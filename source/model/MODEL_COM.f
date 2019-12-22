      MODULE MODEL_COM
!@sum  MODEL_COM Main model variables, independent of resolution
!@auth Original Development Team
!@ver  1.0
      USE RESOLUTION, only : im,jm,lm,ls1,kep,istrat,
     *     psf,pmtop,ptop,psfmpt,pstrat,plbot
      IMPLICIT NONE
      SAVE

!@param IMH half the number of latitudinal boxes
      INTEGER, PARAMETER :: IMH=IM/2
!@param FIM,BYIM real values related to number of long. grid boxes
      REAL*8, PARAMETER :: FIM=IM, BYIM=1./FIM
!@param JEQ grid box zone around or immediately north of the equator
      INTEGER, PARAMETER :: JEQ=1+JM/2

      CHARACTER*132 XLABEL !@var XLABEL=runID+brief description of run
      INTEGER :: LRUNID    !@var Run name stored in XLABEL(1:LRUNID)
      INTEGER :: LHEAD=15  !@var length of crucial beg of module_headers

!**** Vertical resolution dependent variables (set in INPUT)
!@var SIGE sigma levels at layer interfaces (1)
!!!!  Note:   sige(1)=1,  sige(ls1)=0,  sige(lm+1)=-pstrat/psfmpt
      REAL*8, DIMENSION(LM+1) :: SIGE
!@var SIG,DSIG,byDSIG mid point, depth, 1/depth of sigma levels (1)
      REAL*8, DIMENSION(LM) ::
     &     SIG,    ! = (sige(1:lm)+sige(2:lm+1))*0.5d0,
     &     DSIG,   ! =  sige(1:lm)-sige(2:lm+1),
     &     byDSIG  ! =  1./DSIG

!**** Model control parameters:
!@dbparam KOCEAN: if 0 => specified, if 1 => predicted ocean
!@dbparam MFILTR: if 1 => SLP, if 2 => T, if 3 => SLP&T is filtered
      integer :: KOCEAN = 1, MFILTR = 1
!@dbparam dt_UVfilter dU,dV are multiplied by dt/dt_UVfilter
      REAL*8 :: dt_UVfilter = 0. ! UV-filter is NOT used

!@dbparam X_SDRAG.  SDRAG ~X_SDRAG(1)+X_SDRAG(2)*wind_magnitude
      REAL*8, DIMENSION(2) :: X_SDRAG = (/2.5D-4,2.5D-5/)
!@dbparam C_SDRAG.  SDRAG=C_SDRAG (const.) above PTOP
      REAL*8 :: C_SDRAG = 2.5D-5
      REAL*8, DIMENSION(LS1:LM) :: CSDRAGL
!@dbparam P_CSDRAG pressure level above which const.drag is increased
      REAL*8 :: P_CSDRAG=0.
!@dbparam P(P)_SDRAG pressure level above which SDRAG is applied (mb)
      REAL*8 :: P_SDRAG=0., PP_SDRAG = 1.d0 ! (PP_... near poles)
!@var L(P)SDRAG level above which SDRAG is applied (near pole)
      INTEGER :: LSDRAG=LM, LPSDRAG=LM  ! non-polar, polar limit
!@var ANG_SDRAG if =1: ang.momentum lost by SDRAG is added in below PTOP
      INTEGER :: ANG_SDRAG=0  ! default: SDRAG does NOT conserve ang.mom
!@dbparam Wc_JDRAG critical velocity for J.Hansen/Judith Perlwitz drag
      REAL*8 :: Wc_JDRAG=30.d0  !  if 0.: no JDRAG-feature in Sdrag

!**** Diagnostic control parameters
!@dbparam KCOPY: if 1 => acc, if 2 => +rsf, if 3 => +od are saved
!@dbparam NMONAV number of months in a diagnostic accuml. period
!@dbparam Kvflxo if 1 => vert.fluxes into ocean are saved daily
!@dbparam Kradia if -1 save data for, if 1|2 do   inst|adj forcing run
!@dbparam NIPRNT number of instantaneous initial printouts
      integer :: KCOPY=2, NMONAV=1, Kvflxo=0, Kradia=0,iu_rad, NIPRNT=1

C**** (Simplified) Calendar Related Terms
!@param JDperY,JMperY    number of days,months per year
!@var   JDendOfM(0:12)   last Julian day in month
!@var   JDmidOfM(0:13)   middle Julian day in month
      integer, PARAMETER :: JDPERY = 365, JMPERY = 12
      integer :: JDendOfM(0:JMPERY) = (
     *     /0,31,59,90,120,151,181,212,243,273,304,334,365/)
      integer :: JDmidOfM(0:JMPERY+1) = (
     *     /-15,16,47,75,106,136,167,197,228,259,289,320,350,381/)

!@var AMON,AMONTH(0:12)  (3-4 letter) names for current,all months
!@var AMON0  (3-4 letter) name of first month of the current acc-period
      CHARACTER*4 :: AMON='none',AMON0='none', AMONTH(0:12) = (/'IC  ',
     *  'JAN ','FEB ','MAR ','APR ','MAY ','JUNE',
     *  'JULY','AUG ','SEP ','OCT ','NOV ','DEC '/)

!@var NDAY and IYEAR1 relate CALENDAR TIME and INTERNAL TIME Itime :
!@var NDAY number of Internal Time Units per day (1 ITU = DTsrc sec)
!@nlparam IYEAR1  year 1 of internal clock (Itime=0 to 365*NDAY)
      INTEGER :: NDAY,IYEAR1=-1   !@var relate internal to calendar time

!@var ITIME current time in ITUs (1 ITU = DTsrc sec, currently 1 hour)
!@var JDAY,JMON,JDATE,JYEAR,JHOUR current Julian day,month,day,year,hour
      INTEGER :: Itime,JDAY,JMON,JDATE,JYEAR,JHOUR
!@var ItimeI,ItimeE   time at start,end of run
!@var Itime0          time at start of current accumulation period
!@var JMON0,JDATE0,JYEAR0,JHOUR0 date-info about Itime0 (beg.of acc.per)
      INTEGER :: ItimeI,ItimeE,   Itime0,JMON0,JDATE0,JYEAR0,JHOUR0

!@dbparam DTSRC source time step (s)   = 1 ITU
      REAL*8 :: DTsrc = 3600.
!@dbparam DT (atmospheric) dynamics time step (s)
      REAL*8 :: DT    =  450.         ! DT = DTdyn_atm

C**** Time step related multipliers:  N... NI...
C**** general rule:   DTxxx = Nxxx*DTsrc  and  DTxxx = DTsrc/NIxxx
C**** except that the time steps related to NDAa, NDA5k, NDAsf are
C**** slightly larger, to sample all points within the cycle

!@var     NIdyn:  DT atm_dyn  =  DTsrc/NIdyn     (NIdyn=DTsrc/DT)
!@dbparam NIsurf: DT_Surface  =  DTsrc/NIsurf
!@dbparam NRad:   DT_Rad      =  NRad*DTsrc
!@dbparam NFILTR: DT_filter   =  NFILTR*DTsrc
      INTEGER :: NIdyn, NIsurf = 2, NRad = 5 , NFILTR = 1

!@dbparam Ndisk:  DT_saversf    =  Ndisk *DTsrc fort.1/fort.2 saves
!@dbparam Nssw:   DT_checkSsw   =  Nssw  *DTsrc
      INTEGER :: NDisk = 24, Nssw = 1

!@dbparam NDAA:   DT_DiagA    =  NDAA*DTsrc + 2*DT(dyn)
!@dbparam NDA5k:  DT_Diag5k   =  NDA5k*DTsrc + 2*DT(dyn) SpAnal KE
!@dbparam NDA5d:  DT_Diag5d   =  NDA5d*DTsrc     Consrv  SpAnal dyn
!@dbparam NDA5s:  DT_Diag5s   =  NDA5s*DTsrc     Consrv  SpAnal src
!@dbparam NDASf:  DT_DiagSrfc =  NDASf*DTsrc + DTsrc/NIsurf
!@dbparam NDA4:   DT_Diag4    =  NDA4 *DTsrc   Energy history
      INTEGER :: NDAa=7, NDA5d=1, NDA5k=7, NDA5s=1, NDASf=1, NDA4=24

!**** Accounting variables
!@dbparam IRAND last seed used by rand.number generator
!@var KDISK next rsf (fort.)1 or 2 to be written to
!@var NSTEP number of dynamics steps since start of run
!@var MRCH  flags position in dynamics cycle (>0 fw, <0 bw step)
      INTEGER :: IRAND=123456789, KDISK=1, NSTEP,MRCH
!@var MODRD,MODD5K,MODD5S: if MODxxx=0 do xxx, else skip xxx
      INTEGER :: MODRD, MODD5K, MODD5S
!@var MDYN,MCNDS,MRAD,MSURF,MDIAG,MELSE timing-indices
      INTEGER  MDYN,MCNDS,MRAD,MSURF,MDIAG,MELSE
!@param NSAMPL number of diagnostic sampling schemes
      INTEGER, PARAMETER :: NSAMPL = 12
!@var IDACC(NSAMPL) counters for diagn. accumulations
      INTEGER, DIMENSION(NSAMPL) :: IDACC

!**** Boundary condition arrays:
!@var ZATMO,HLAKE Topography arrays: elevation (m), lake depth (m) ???
      REAL*8, DIMENSION(IM,JM) :: ZATMO,HLAKE,
!@var Fxx fraction of gridbox of type xx (land,ocean,...)
     *     FLAND,FOCEAN,FLICE,FLAKE0,FEARTH
!@var VDATA(:,:,k)  fraction of gridbox of veg.type k=1-11
      REAL*8, DIMENSION(IM,JM,11) :: VDATA
!@var WFCS water field capacity of first ground layer (kg/m2)  ???
      REAL*8, DIMENSION(IM,JM) :: WFCS

!**** IO read/write flags used by the io_xyz routines
!@param IOWRITE Flag used for writing normal restart files
!@param IOWRITE_SINGLE Flag used for saving diags in single precision
!@param IOWRITE_MON Flag used for saving restart part only (no diags)
!@param IOREAD Flag used for reading in (composite) restart files
!@param IOREADNT Flag used for reading in restart files (w/o tracers)
!@param IRSFIC Flag used for reading in restart part to start NEW run
!@param IRSFICNT Flag used for reading restart (w/o tracers) for NEW run
!@param IRSFICNO Flag used for reading restart (w/o ocean) for NEW run
!@param IRERUN Flag used for reading in restart part to extend OLD run
      INTEGER, PARAMETER :: ioread=1,ioread_single=2,
     *     irerun=3,irsfic=4,irsficnt=5,ioreadnt=6,irsficno=7,
     *     iowrite=-1,iowrite_single=-2,iowrite_mon=-3

!**** Main model prognostic variables
!@var U,V east-west, and north-south velocities (m/s)
!@var T potential temperature (referenced to 1 mb) (K)
!@var Q specific humidity (kg water vapor/kg air)
!@var WM cloud liquid water amount (kg water/kg air)
      REAL*8, DIMENSION(IM,JM,LM) :: U,V,T,Q,WM
!@var P surface pressure (hecto-Pascals - PTOP)
      REAL*8, DIMENSION(IM,JM) :: P

C**** Define surface types (mostly used for weighting diagnostics)
!@param NTYPE number of different surface types
      INTEGER, PARAMETER :: NTYPE=6   ! orig = 3
!@var FTYPE fractions of each surface type
      REAL*8, DIMENSION(NTYPE,IM,JM) :: FTYPE
!@param ITxx indices of various types (used only when it matters)
      INTEGER, PARAMETER :: ITOCEAN=1, ITOICE=2, ITEARTH=3,
     *                      ITLANDI=4, ITLAKE=5, ITLKICE=6

C**** Variables specific for stratosphere and/or strat diagnostics
!@var AIRX, AIRMX*DXYP(J) Used in stratosphere model. (mb m^2)
      REAL*8 AIRX(IM,JM)
!@var LMC max layer of mc convective mass flux. (Strat model)
      INTEGER, DIMENSION(2,IM,JM) :: LMC
!@var DO_GWDRAG when true, prints Gravity Wave diagnostics
      LOGICAL :: DO_GWDRAG = .false.
!@var iDO_GWDRAG number if AIJ Gravity wave diagnostics
      INTEGER :: iDO_GWDRAG = 0

!@dbparam UOdrag parameter that decides whether ocean.ice velocities
!@+   feed into drag calculation in surface (default = 0)
      INTEGER :: UOdrag = 0

!@nlparam QCHECK TRUE for running diagnostic checks
      LOGICAL :: QCHECK = .FALSE.

!@var stop_on TRUE stops the model (set with "kill -15 PID)
      LOGICAL :: stop_on = .FALSE.

      END MODULE MODEL_COM

      MODULE TIMINGS
!@sum  TIMINGS contains variables for keeping track of computing time
!@auth Gavin Schmidt
!@ver  1.0
      IMPLICIT NONE
      SAVE
!@param NTIMEMAX maximum number of possible time accumulators
      INTEGER, PARAMETER :: NTIMEMAX=10
!@var NTIMEACC actual number of time accumulators
      INTEGER :: NTIMEACC = 0
!@var TIMING array that holds timing info
      INTEGER, DIMENSION(0:NTIMEMAX) :: TIMING
!@var TIMESTR array that holds timing info description
      CHARACTER*12, DIMENSION(NTIMEMAX) :: TIMESTR

      END MODULE TIMINGS

      SUBROUTINE SET_TIMER(STR,MINDEX)
!@sum  SET_TIMER sets an index of TIMING for a particular description
!@auth Gavin Schmidt
!@ver  1.0
      USE TIMINGS
      IMPLICIT NONE
!@var STR string that describes timing accumulator
      CHARACTER*12, INTENT(IN) :: STR
!@var MINDEX index for that accumulator
      INTEGER, INTENT(OUT) :: MINDEX
      INTEGER N

C**** Check whether index has been set
      DO N=1,NTIMEACC
        IF (STR.EQ.TIMESTR(N)) THEN
          MINDEX=N
          RETURN
        END IF
      END DO
C**** Otherwise increase number of indexes
      NTIMEACC = NTIMEACC + 1
      IF (NTIMEACC.gt.NTIMEMAX) call stop_model(
     &     "Too many timing indices: increase NTIMEMAX",255)
      MINDEX = NTIMEACC
      TIMESTR(MINDEX) = STR
C****
      RETURN
      END SUBROUTINE SET_TIMER

      SUBROUTINE TIMER (MNOW,MSUM)
!@sum  TIMER keeps track of elapsed CPU time in hundredths of seconds
!@auth Gary Russell
!@ver  1.0
      USE TIMINGS
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: MNOW   !@var MNOW current CPU time (.01 s)
      INTEGER, INTENT(INOUT) :: MSUM !@var MSUM index for running total
      INTEGER :: MINC                !@var MINC time since last call
      INTEGER, SAVE :: MLAST = 0     !@var MLAST  last CPU time

      CALL GETTIME(MNOW)
      MINC  = MNOW - MLAST
      TIMING(MSUM)  = TIMING(MSUM) + MINC
      MLAST = MNOW
      RETURN
      END SUBROUTINE TIMER

      SUBROUTINE TIMEOUT (MBEGIN,MIN,MOUT)
!@sum  TIMEOUT redistributes timing info between counters
!@auth Gary Russell
!@ver  1.0
      USE TIMINGS
      IMPLICIT NONE
!@var MBEGIN CPU time start of section (.01 s)
      INTEGER, INTENT(IN) :: MBEGIN
      INTEGER, INTENT(INOUT) :: MIN  !@var MIN index to be added to
      INTEGER, INTENT(INOUT) :: MOUT !@var MOUT index to be taken from
      INTEGER :: MINC                !@var MINC time since MBEGIN
      INTEGER :: MNOW                !@var MNOW current CPU time (.01 s)

      CALL GETTIME(MNOW)
      MINC  = MNOW - MBEGIN
      TIMING(MIN)  = TIMING(MIN)  + MINC
      TIMING(MOUT) = TIMING(MOUT) - MINC
      RETURN
      END SUBROUTINE TIMEOUT

      SUBROUTINE io_label(kunit,it,itm,iaction,ioerr)
!@sum  io_label reads and writes label/parameters to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM
      USE TIMINGS, only : ntimemax,ntimeacc,timestr,timing
      USE PARAM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var it input/ouput value of hour
!@var itm maximum hour returned (different from it if post-processing)
      INTEGER, INTENT(INOUT) :: it,itm
!@var LABEL2 content of record 2
      CHARACTER*80 :: LABEL2
!@var NTIM1,TSTR1,TIM1 timing related dummy arrays
      INTEGER NTIM1,TIM1(NTIMEMAX)
      CHARACTER*12 TSTR1(NTIMEMAX)
!@var ITmin,ITmax minimal/maximal time in acc periods to be combined
      INTEGER, SAVE :: ITmax=-1, ITmin=-1 ! to protect against long runs
      INTEGER it1,it0,idac1(12),nd1,iy1,iti1,ite1,it01,im0,jm0,lm0,ls10

      SELECT CASE (IACTION)
      CASE (:IOWRITE)           ! output to end-of-month restart file
        WRITE (kunit,err=10) it,XLABEL,nday,iyear1,itimei,itimee,itime0,
     *       NTIMEACC,TIMING(1:NTIMEACC),TIMESTR(1:NTIMEACC)
C**** doc line: basic model parameters
        write(label2,'(a13,4i4,a)') 'IM,JM,LM,LS1=',im,jm,lm,ls1,' '
        WRITE (kunit,err=10) LABEL2
C**** write parameters database here
        call write_param(kunit)
      CASE (IOREAD:)          ! label always from input file
        READ (kunit,err=10) it,XLABEL,nd1,iy1,iti1,ite1,it01,
     *        NTIM1,TIM1(1:NTIM1),TSTR1(1:NTIM1)
C**** use doc-record to check the basic model parameters
       READ (kunit,err=10) LABEL2
        READ (label2,'(13x,4i4)',err=10) im0,jm0,lm0,ls10
        if (im.ne.im0.or.jm.ne.jm0.or.lm.ne.lm0.or.ls10.ne.ls1) then
          ioerr = 0   ! warning
        end if
        SELECT CASE (IACTION)   ! set model common according to iaction
        CASE (ioread)           ! parameters from rundeck & restart file
          call read_param(kunit,.false.)
          nday=nd1 ; itimei=iti1      ! changeable only at new starts
          itimee=ite1 ; itime0=it01   ! are changed later if appropriate
          if (iyear1.lt.0) iyear1=iy1 ! rarely changes on restarts
          NTIMEACC=NTIM1
          TIMESTR(1:NTIM1)=TSTR1(1:NTIM1)
          TIMING(1:NTIM1)=TIM1(1:NTIM1)
        CASE (IRSFIC,irsficnt,IRSFICNO) ! rundeck/defaults except label
          read(kunit,err=10)          ! skip parameters, dates
          it=it*24/nd1                ! switch itime to ihour
        CASE (IRERUN)           ! parameters from rundeck & restart file
          call read_param(kunit,.false.)
          nday=nd1 ; itimei=iti1      ! changeable only at new starts
          itimee=ite1 ; itime0=it01   ! is changed later if appropriate
          if (iyear1.lt.0) iyear1=iy1 ! rarely changes on restart/reruns
        CASE (IOREAD_SINGLE)    ! parameters/label from 1-many acc files
          call read_param(kunit,.true.)  ! ignore rundeck
          call sync_param( "kradia",kradia)
          nday=nd1 ; iyear1=iy1 ; itime0=it01
          NTIMEACC=NTIM1                 ! use timing from current file
          TIMESTR(1:NTIM1)=TSTR1(1:NTIM1)
          TIMING(1:NTIM1)=TIMING(1:NTIM1)+TIM1(1:NTIM1)

C**** keep track of min/max time over the combined diagnostic period
          if (it.gt.ITmax)                   ITmax = it
          if (ITmin.lt.0 .or. it01.lt.ITmin) ITmin = it01
          itime0 = ITmin

        END SELECT
      END SELECT
      itm = max(it,ITmax)
      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_label

      SUBROUTINE io_model(kunit,iaction,ioerr)
!@sum  io_model reads and writes model variables to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "MODEL01"

      MODULE_HEADER(lhead+1:80) = 'R8 dim(im,jm,lm):u,v,t, p(im,jm),'//
     *  ' dim(im,jm,lm):q,MliqW'

      SELECT CASE (IACTION)
      CASE (:IOWRITE) ! output to end-of-month restart file
        WRITE (kunit,err=10) MODULE_HEADER,U,V,T,P,Q,WM
      CASE (IOREAD:)          ! input from restart file
        READ (kunit,err=10) HEADER,U,V,T,P,Q,WM
        IF (HEADER(1:LHEAD).ne.MODULE_HEADER(1:LHEAD)) THEN
          PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      END SELECT
      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_model

      subroutine getdte(It,Nday,Iyr0,Jyr,Jmn,Jd,Jdate,Jhour,amn)
!@sum  getdte gets julian calendar info from internal timing info
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : hrday
      USE MODEL_COM, only : JDperY,JDendOfM,amonth
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: It,Nday,Iyr0
      INTEGER, INTENT(OUT) :: Jyr,Jmn,Jd,Jdate,Jhour
      CHARACTER*4, INTENT(OUT) :: amn

      Jyr=Iyr0+It/(Nday*JDperY)
      Jd=1+It/Nday-(Jyr-Iyr0)*JDperY
      Jmn=1
      do while (Jd.GT.JDendOfM(Jmn))
        Jmn=Jmn+1
      end do
      Jdate=Jd-JDendOfM(Jmn-1)
      Jhour=nint(mod(It*hrday/Nday,hrday))
      amn=amonth(Jmn)

      return
      end subroutine getdte

