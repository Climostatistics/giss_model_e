#include "rundeck_opts.h"

      SUBROUTINE OCEANS
!@sum  OCEANS integrates ocean source terms and dynamics
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : rhows,grav
      USE MODEL_COM, only : idacc,modd5s,msurf
#ifdef TRACERS_OCEAN
      USE TRACER_COM, only : t_qlimit,ntm
      USE OCEAN, only : trmo,txmo,tymo,tzmo
#endif
      USE OCEAN, only : im,jm,lmo,ndyno,mo,g0m,gxmo,gymo,gzmo,s0m,sxmo,
     *     symo,szmo,dts,dtofs,dto,dtolf,bydxypo,mdyno,msgso
     *     ,ogeoz,ogeoz_sv,opbot,ze,lmm,imaxj
      USE ODIAG, only : oijl,oij,ijl_mo,ijl_g0m,ijl_s0m,ijl_gflx
     *     ,ijl_sflx,ijl_mfu,ijl_mfv,ijl_mfw,ijl_ggmfl,ijl_sgmfl,ij_ssh
     *     ,ij_pb
#ifdef TRACERS_OCEAN
     *     ,toijl,toijl_conc,toijl_tflx,toijl_gmfl
#endif
      USE OCEAN_DYN, only : mmi,smu,smv,smw
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,JM,LMO) :: MM0=0,MM1=0,MMX=0,UM0=0,VM0=0,
     *     UM1=0,VM1=0
      INTEGER NS,I,J,L,mnow,n
C****
C**** Integrate Ocean Dynamics terms
C****
      OGEOZ_SV=OGEOZ

C**** Apply surface fluxes to ocean
      CALL GROUND_OC
         CALL CHECKO('GRNDOC')

C**** Apply ice/ocean and air/ocean stress to ocean
      CALL OSTRES
         CALL CHECKO('OSTRES')
         CALL TIMER (MNOW,MSURF)
         IF (MODD5S.EQ.0) CALL DIAGCA (11)

C**** Calculate vertical diffusion
      CALL OCONV
         CALL CHECKO('OCONV ')
         CALL TIMER (MNOW,MSGSO)

      IDACC(11) = IDACC(11) + 1

!$OMP PARALLEL DO  PRIVATE(L)
      DO L=1,LMO
        SMU(:,:,L) = 0. ; SMV(:,:,L) = 0. ; SMW(:,:,L) = 0.
      END DO
!$OMP END PARALLEL DO

      CALL OVTOM (MMI,UM0,VM0)
      CALL OPGF0
      NS=NDYNO
C**** Initial Forward step,  QMX = QM0 + DT*F(Q0)
      CALL OFLUX  (NS,MMI,.FALSE.)
      CALL OADVM (MM1,MMI,DTOFS)
      CALL OADVV (UM1,VM1,UM0,VM0,DTOFS)
      CALL OPGF  (UM1,VM1,DTOFS)
      CALL OMTOV (MM1,UM1,VM1)
C**** Initial Backward step,  QM1 = QM0 + DT*F(Q0)
      CALL OFLUX  (NS,MMI,.FALSE.)
      CALL OADVM (MM1,MMI,DTO)
      CALL OADVV (UM1,VM1,UM0,VM0,DTO)
      CALL OPGF  (UM1,VM1,DTO)
      CALL OMTOV (MM1,UM1,VM1)
C**** First even leap frog step,  Q2 = Q0 + 2*DT*F(Q1)
      CALL OPGF0
      CALL OFLUX  (NS,MMI,.TRUE.)
      CALL OADVM (MM0,MMI,DTOLF)
      CALL OADVV (UM0,VM0,UM0,VM0,DTOLF)
      CALL OPGF  (UM0,VM0,DTOLF)
      CALL OMTOV (MM0,UM0,VM0)
      NS=NS-1
C**** Odd leap frog step,  Q3 = Q1 + 2*DT*F(Q2)
  420 CALL OPGF0
      CALL OFLUX  (NS,MM1,.FALSE.)
      CALL OADVM (MM1,MM1,DTOLF)
      CALL OADVV (UM1,VM1,UM1,VM1,DTOLF)
      CALL OPGF  (UM1,VM1,DTOLF)
      CALL OMTOV (MM1,UM1,VM1)
      NS=NS-1
C**** Even leap frog step,  Q4 = Q2 + 2*DT*F(Q3)
      CALL OPGF0
      CALL OFLUX  (NS,MM0,.TRUE.)
      CALL OADVM (MM0,MM0,DTOLF)
      CALL OADVV (UM0,VM0,UM0,VM0,DTOLF)
      CALL OPGF  (UM0,VM0,DTOLF)
      CALL OMTOV (MM0,UM0,VM0)
      NS=NS-1
      IF(NS.GT.1)  GO TO 420
!$OMP PARALLEL DO  PRIVATE(L)
      DO L=1,LMO
        OIJL(:,:,L,IJL_MFU) = OIJL(:,:,L,IJL_MFU) + SMU(:,:,L)
        OIJL(:,:,L,IJL_MFV) = OIJL(:,:,L,IJL_MFV) + SMV(:,:,L)
        OIJL(:,:,L,IJL_MFW) = OIJL(:,:,L,IJL_MFW) + SMW(:,:,L)
      END DO
!$OMP END PARALLEL DO
C**** Advection of Potential Enthalpy and Salt
      CALL OADVT (G0M,GXMO,GYMO,GZMO,DTOLF,.FALSE.,OIJL(1,1,1,IJL_GFLX))
      CALL OADVT (S0M,SXMO,SYMO,SZMO,DTOLF,.TRUE.,OIJL(1,1,1,IJL_SFLX))
#ifdef TRACERS_OCEAN
      DO N=1,NTM
        CALL OADVT(TRMO(1,1,1,N),TXMO(1,1,1,N),TYMO(1,1,1,N),
     *       TZMO(1,1,1,N),DTOLF,t_qlimit(n),TOIJL(1,1,1,TOIJL_TFLX,N))
      END DO
#endif
        CALL CHECKO ('OADVT ')
!$OMP PARALLEL DO PRIVATE(L)
        DO L=1,LMO
          OIJL(:,:,L,IJL_MO)  = OIJL(:,:,L,IJL_MO) +  MO(:,:,L)
          OIJL(:,:,L,IJL_G0M) = OIJL(:,:,L,IJL_G0M) + G0M(:,:,L)
          OIJL(:,:,L,IJL_S0M) = OIJL(:,:,L,IJL_S0M) + S0M(:,:,L)
        END DO
!$OMP END PARALLEL DO
        DO J=1,JM
          DO I=1,IMAXJ(J)
            OIJ(I,J,IJ_SSH) = OIJ(I,J,IJ_SSH) + OGEOZ(I,J)
            OIJ(I,J,IJ_PB)  = OIJ(I,J,IJ_PB)  + (OPBOT(I,J)-ZE(LMM(I,J))
     *           *RHOWS*GRAV) 
          END DO
        END DO

#ifdef TRACERS_OCEAN
        DO N=1,NTM
!$OMP PARALLEL DO PRIVATE(L)
          DO L=1,LMO
            TOIJL(:,:,L,TOIJL_CONC,N)=TOIJL(:,:,L,TOIJL_CONC,N)
     *           +TRMO(:,:,L,N)
          END DO
!$OMP END PARALLEL DO
        END DO
#endif
        CALL TIMER (MNOW,MDYNO)

        IF (MODD5S.EQ.0) CALL DIAGCA (12)
C**** Apply Wajowicz horizontal diffusion to UO and VO ocean currents
      CALL ODIFF(DTS)
        CALL CHECKO ('ODIFF ')
C**** Apply GM + Redi tracer fluxes
      CALL GMKDIF
      CALL GMFEXP(G0M,GXMO,GYMO,GZMO,.FALSE.,OIJL(1,1,1,IJL_GGMFL))
      CALL GMFEXP(S0M,SXMO,SYMO,SZMO,.TRUE.,OIJL(1,1,1,IJL_SGMFL))
#ifdef TRACERS_OCEAN
      DO N = 1,NTM
        CALL GMFEXP(TRMO(1,1,1,N),TXMO(1,1,1,N),TYMO(1,1,1,N),
     *       TZMO(1,1,1,N),t_qlimit(n),TOIJL(1,1,1,TOIJL_GMFL,N))
      END DO
#endif
        CALL CHECKO ('GMDIFF')
        CALL TIMER (MNOW,MSGSO)
C****
C**** Acceleration and advection of tracers through ocean straits
C****
      CALL STPGF
      CALL STCONV
      CALL STBDRA
        CALL CHECKO ('STBDRA')
      CALL STADV
        CALL CHECKO ('STADV ')
C**** remove STADVI since it is not really consistent with ICEDYN
c      CALL STADVI
c        CALL CHECKO ('STADVI')
#ifdef TRACERS_OCEAN
        CALL OC_TDECAY
#endif
        CALL TIMER (MNOW,MSGSO)
      CALL TOC2SST

      RETURN
      END SUBROUTINE OCEANS

      SUBROUTINE init_OCEAN(iniOCEAN,istart)
!@sum init_OCEAN initiallises ocean variables
!@auth Original Development Team
!@ver  1.0
      USE FILEMANAGER, only : openunit,closeunit
      USE TIMINGS, only : timing,ntimeacc
      USE PARAM
      USE CONSTANT, only : twopi,radius,by3,grav
      USE MODEL_COM, only : dtsrc,kocean
      USE OCEAN, only : im,jm,lmo,focean,ze1,zerat,sigeo,dsigo,sigo,lmm
     *     ,lmu,lmv,hatmo,hocean,ze,mo,g0m,gxmo,gymo,gzmo,s0m,sxmo
     *     ,symo,szmo,uo,vo,dxypo,ogeoz,dts,dtolf,dto,dtofs,mdyno,msgso
     *     ,ndyno,imaxj,bydxypo,ogeoz_sv,bydts
      USE OCFUNC, only : vgsp,tgsp,hgsp,agsp,bgsp,cgs
      USE SW2OCEAN, only : init_solar
      IMPLICIT NONE
      INTEGER I,J,L,N,iu_OIC,iu_OFTAB,IP1,IM1,LMIJ,I1,J1,I2,J2
     *     ,iu_TOPO
      REAL*4, DIMENSION(IM,JM,LMO):: MO4,G0M4,S0M4,GZM4,SZM4
      CHARACTER*80 TITLE
      REAL*8 FJEQ,SM,SG0,SGZ,SS0,SSZ
      LOGICAL, INTENT(IN) :: iniOCEAN
      INTEGER, INTENT(IN) :: istart
C****
C**** Check that KOCEAN is set correctly
C****
      IF (KOCEAN.eq.0) THEN
        call stop_model(
     &       "Must have KOCEAN > 0 for interactive ocean runs",255)
      END IF
C****
C**** set up time steps from atmospheric model
C****
      call sync_param("DTO",DTO)

      DTS=DTSRC
      BYDTS=1d0/DTS
      NDYNO=2*NINT(.5*DTS/DTO)
      DTO=DTS/NDYNO
      DTOLF=2.*DTO
      DTOFS=2.*DTO*BY3
C**** Set up timing indexes
      CALL SET_TIMER(" OCEAN DYNAM",MDYNO)
      CALL SET_TIMER(" OCEAN PHYS.",MSGSO)
C****
C**** Arrays needed each ocean model run
C****
      CALL GEOMO
C**** Calculate ZE, SIGEO, DSIGO and SIGO
      DO 110 L=0,LMO
  110 ZE(L) = ZE1*(ZERAT**L-1d0)/(ZERAT-1d0)
      SIGEO(0) = 0.
      DO 120 L=1,LMO
      SIGEO(L) = ZE(L)/ZE(LMO)
      DSIGO(L) = SIGEO(L) - SIGEO(L-1)
  120 SIGO(L) = (SIGEO(L) + SIGEO(L-1))*5d-1
C**** Read in table function for specific volume
      CALL openunit("OFTAB",iu_OFTAB,.TRUE.,.TRUE.)
      READ  (iu_OFTAB) TITLE,VGSP
      WRITE (6,*) 'Read from unit ',iu_OFTAB,': ',TITLE
      READ  (iu_OFTAB) TITLE,TGSP
      WRITE (6,*) 'Read from unit ',iu_OFTAB,': ',TITLE
      READ  (iu_OFTAB) TITLE,CGS
      WRITE (6,*) 'Read from unit ',iu_OFTAB,': ',TITLE
      READ  (iu_OFTAB) TITLE,HGSP
      WRITE (6,*) 'Read from unit ',iu_OFTAB,': ',TITLE
      READ  (iu_OFTAB) TITLE,AGSP
      WRITE (6,*) 'Read from unit ',iu_OFTAB,': ',TITLE
      READ  (iu_OFTAB) TITLE,BGSP
      WRITE (6,*) 'Read from unit ',iu_OFTAB,': ',TITLE
      call closeunit(iu_OFTAB)

C**** READ IN LANDMASKS AND TOPOGRAPHIC DATA
      call openunit("TOPO_OC",iu_TOPO,.true.,.true.)
      CALL READT (iu_TOPO,0,FOCEAN,IM*JM,FOCEAN,1) ! Ocean fraction
      CALL READT (iu_TOPO,0,HATMO ,IM*JM,HATMO ,4) ! Atmo. Topography
      CALL READT (iu_TOPO,0,HOCEAN,IM*JM,HOCEAN,1) ! Ocean depths
      call closeunit(iu_TOPO)

C**** Calculate LMM and modify HOCEAN
      DO 170 J=1,JM
      DO 170 I=1,IM
      LMM(I,J) =0
      IF(FOCEAN(I,J).LE.0.)  GO TO 170
      DO 150 L=1,LMO-1
  150 IF(HATMO(I,J)+HOCEAN(I,J) .le. 5d-1*(ZE(L)+ZE(L+1)))  GO TO 160
C     L=LMO
  160 LMM(I,J)=L
      HOCEAN(I,J) = -HATMO(I,J) + ZE(L)
  170 CONTINUE
C**** Calculate LMU
      I=IM
      DO 180 J=1,JM
      DO 180 IP1=1,IM
      LMU(I,J) = MIN(LMM(I,J),LMM(IP1,J))
  180 I=IP1
C**** Calculate LMV
      DO 190 J=1,JM-1
      DO 190 I=1,IM
  190 LMV(I,J) = MIN(LMM(I,J),LMM(I,J+1))
C****
      IF(iniOCEAN) THEN
C**** Initialize a run from ocean initial conditions
      CALL openunit("OIC",iu_OIC,.TRUE.,.TRUE.)
      READ  (iu_OIC,ERR=820) TITLE,MO4,G0M4,GZM4,S0M4,SZM4
      call closeunit(iu_OIC)
      WRITE (6,*) 'OIC read from unit ',iu_OIC,': ',TITLE
C**** Calculate layer mass from column mass and check for consistency
      DO 313 J=1,JM
      DO 313 I=1,IM
      LMIJ=LMM(I,J)
      DO 311 L=1,LMIJ
  311 MO(I,J,L) = MO4(I,J,L)
      DO 312 L=LMIJ+1,LMO
  312 MO(I,J,L) = 0.
      IF((LMM(I,J).GT.0).AND.(ABS(MO(I,J,1)/ZE(1)-1d3).GT.5d1))
     *  WRITE (6,931) I,J,LMIJ,MO(I,J,1),ZE(1)
  313 CONTINUE
C**** Initialize velocity field to zero
      UO=0
      VO=0
C**** Define mean value of mass, potential heat, and salinity at poles
      DO 370 L=1,LMO
      J=JM
      SM  = 0.
      SG0 = 0.
      SGZ = 0.
      SS0 = 0.
      SSZ = 0.
      DO 331 I=1,IM
      SM  = SM  +   MO(I,J,L)
      SG0 = SG0 + G0M4(I,J,L)
      SGZ = SGZ + GZM4(I,J,L)
      SS0 = SS0 + S0M4(I,J,L)
  331 SSZ = SSZ + SZM4(I,J,L)
      DO 332 I=1,IM
      MO(I,J,L)   = SM /IM
      G0M4(I,J,L) = SG0/IM
      GZM4(I,J,L) = SGZ/IM
      S0M4(I,J,L) = SS0/IM
  332 SZM4(I,J,L) = SSZ/IM
C**** Define East-West horizontal gradients
      GXMO=0
      GYMO=0
      SXMO=0
      SYMO=0
      IM1=IM-1
      I=IM
      DO 345 J=2,JM-1
      DO 345 IP1=1,IM
      IF(LMM(I  ,J).LT.L)  GO TO 344
      IF(LMM(IM1,J).GE.L)  GO TO 342
      IF(LMM(IP1,J).LT.L)  GO TO 344
      GXMO(I,J,L) = .5*(G0M4(IP1,J,L)-G0M4(I,J,L))
      SXMO(I,J,L) = .5*(S0M4(IP1,J,L)-S0M4(I,J,L))
      GO TO 344
  342 IF(LMM(IP1,J).GE.L)  GO TO 343
      GXMO(I,J,L) = .5*(G0M4(I,J,L)-G0M4(IM1,J,L))
      SXMO(I,J,L) = .5*(S0M4(I,J,L)-S0M4(IM1,J,L))
      GO TO 344
  343 GXMO(I,J,L) = .25*(G0M4(IP1,J,L)-G0M4(IM1,J,L))
      SXMO(I,J,L) = .25*(S0M4(IP1,J,L)-S0M4(IM1,J,L))
  344 IM1=I
  345 I=IP1
C**** Define North-South horizontal gradients
      DO 354 J=2,JM-1
      DO 354 I=1,IM
      IF(LMM(I,J  ).LT.L)  GO TO 354
      IF(LMM(I,J-1).GE.L)  GO TO 352
      IF(LMM(I,J+1).LT.L)  GO TO 354
      GYMO(I,J,L) = .5*(G0M4(I,J+1,L)-G0M4(I,J,L))
      SYMO(I,J,L) = .5*(S0M4(I,J+1,L)-S0M4(I,J,L))
      GO TO 354
  352 IF(LMM(I,J+1).GE.L)  GO TO 353
      GYMO(I,J,L) = .5*(G0M4(I,J,L)-G0M4(I,J-1,L))
      SYMO(I,J,L) = .5*(S0M4(I,J,L)-S0M4(I,J-1,L))
      GO TO 354
  353 GYMO(I,J,L) = .25*(G0M4(I,J+1,L)-G0M4(I,J-1,L))
      SYMO(I,J,L) = .25*(S0M4(I,J+1,L)-S0M4(I,J-1,L))
  354 CONTINUE
C**** Multiply specific quantities by mass
      DO 360 J=1,JM
      DO 360 I=1,IM
      G0M(I,J,L)  = G0M4(I,J,L)*(MO(I,J,L)*DXYPO(J))
      GXMO(I,J,L) = GXMO(I,J,L)*(MO(I,J,L)*DXYPO(J))
      GYMO(I,J,L) = GYMO(I,J,L)*(MO(I,J,L)*DXYPO(J))
      GZMO(I,J,L) = GZM4(I,J,L)*(MO(I,J,L)*DXYPO(J))
      S0M(I,J,L)  = S0M4(I,J,L)*(MO(I,J,L)*DXYPO(J))
      SXMO(I,J,L) = SXMO(I,J,L)*(MO(I,J,L)*DXYPO(J))
      SYMO(I,J,L) = SYMO(I,J,L)*(MO(I,J,L)*DXYPO(J))
  360 SZMO(I,J,L) = SZM4(I,J,L)*(MO(I,J,L)*DXYPO(J))
  370 CONTINUE
C**** Initiallise geopotential field (needed by KPP)
      OGEOZ = 0.
      OGEOZ_SV = 0.

      END IF

C**** Initialize straits arrays
      call init_STRAITS(iniOCEAN)

C**** Initialize solar radiation penetration arrays
      call init_solar

C**** Initialize ocean diagnostics
      call init_ODIAG

C**** Initialize KPP mixing scheme
      call kmixinit(ZE)

#ifdef TRACERS_OCEAN
C**** Set diagnostics for ocean tracers
      call init_tracer_ocean
#endif

C**** Set atmospheric surface variables
      IF (ISTART.gt.0) CALL TOC2SST

      RETURN
C**** Terminate because of improper start up
  820 call stop_model('init_OCEAN: Error reading ocean IC',255)
C****
  931 FORMAT ('0Inconsistency between LMM and M:',3I4,2F10.1)

      RETURN
C****
      END SUBROUTINE init_OCEAN

      SUBROUTINE daily_OCEAN(end_of_day)
!@sum  daily_OCEAN performs the daily tasks for the ocean module
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : sday
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: end_of_day

C**** Only do this at end of the day
      IF (end_of_day) THEN

C**** Add glacial melt from Antarctica and Greenland
        CALL GLMELT(SDAY)

C**** set gtemp arrays for ocean
        CALL TOC2SST
      END IF
C****
      RETURN
      END SUBROUTINE daily_OCEAN

      SUBROUTINE io_ocean(kunit,iaction,ioerr)
!@sum  io_ocean is a driver to ocean related io routines
!@auth Gavin Schmidt
!@ver  1.0
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR

      call io_ocdyn  (kunit,iaction,ioerr)
      call io_straits(kunit,iaction,ioerr)

      RETURN
C****
      END SUBROUTINE io_ocean

      SUBROUTINE io_ocdyn(kunit,iaction,ioerr)
!@sum  io_ocdyn reads and writes ocean arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,irsficno,irsfic
     *     ,irsficnt,irerun,lhead
      USE OCEAN
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "OCDYN01"
#ifdef TRACERS_OCEAN
!@var TRHEADER Character string label for individual records
      CHARACTER*80 :: TRHEADER, TRMODULE_HEADER = "TROCDYN02"

      write (TRMODULE_HEADER(lhead+1:80),'(a13,i3,a1,i3,a)')
     *     'R8 dim(im,jm,',LMO,',',NTM,'):TRMO,TX,TY,TZ'
#endif

      write (MODULE_HEADER(lhead+1:80),'(a13,i2,a)') 'R8 dim(im,jm,',
     *   LMO,'):M,U,V,G0,GX,GY,GZ,S0,SX,SY,SZ, OGZ,OGZSV'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,MO,UO,VO,G0M,GXMO,GYMO,GZMO
     *     ,S0M,SXMO,SYMO,SZMO,OGEOZ,OGEOZ_SV
#ifdef TRACERS_OCEAN
       WRITE (kunit,err=10) TRMODULE_HEADER,TRMO,TXMO,TYMO,TZMO
#endif
      CASE (IOREAD:)            ! input from restart file
        SELECT CASE (IACTION)
          CASE (IRSFICNO)   ! initial conditions (no ocean data)
            READ (kunit)
          CASE (ioread,irerun,irsfic) ! restarts
            READ (kunit,err=10) HEADER,MO,UO,VO,G0M,GXMO,GYMO,GZMO,S0M
     *           ,SXMO,SYMO,SZMO,OGEOZ,OGEOZ_SV
            IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",HEADER
     *             ,MODULE_HEADER
              GO TO 10
            END IF
#ifdef TRACERS_OCEAN
            READ (kunit,err=10) TRHEADER,TRMO,TXMO,TYMO,TZMO
            IF (TRHEADER(1:LHEAD).NE.TRMODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",TRHEADER
     *             ,TRMODULE_HEADER
              GO TO 10
            END IF
#endif
          CASE (irsficnt) ! restarts (never any tracer data)
            READ (kunit,err=10) HEADER,MO,UO,VO,G0M,GXMO,GYMO,GZMO,S0M
     *           ,SXMO,SYMO,SZMO,OGEOZ,OGEOZ_SV
            IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",HEADER
     *             ,MODULE_HEADER
              GO TO 10
            END IF
          END SELECT
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
C****
      END SUBROUTINE io_ocdyn

      SUBROUTINE CHECKO(SUBR)
!@sum  CHECKO Checks whether Ocean variables are reasonable
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : byrt3,teeny
      USE MODEL_COM, only : qcheck
#ifdef TRACERS_OCEAN
      USE TRACER_COM, only : ntm, trname, t_qlimit
#endif
      USE SEAICE_COM, only : rsi
      USE OCEAN
      IMPLICIT NONE
      REAL*8 SALIM,GO1,SO1,relerr,errmax
      LOGICAL QCHECKO
      INTEGER I,J,L,n,imax,jmax,lmax
!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR

C**** Check for NaN/INF in ocean data
      IF (QCHECK) THEN
      CALL CHECK3(MO  ,IM,JM,LMO,SUBR,'mo ')
      CALL CHECK3(G0M ,IM,JM,LMO,SUBR,'g0m')
      CALL CHECK3(GXMO,IM,JM,LMO,SUBR,'gxm')
      CALL CHECK3(GYMO,IM,JM,LMO,SUBR,'gym')
      CALL CHECK3(GZMO,IM,JM,LMO,SUBR,'gzm')
      CALL CHECK3(S0M ,IM,JM,LMO,SUBR,'s0m')
      CALL CHECK3(SXMO,IM,JM,LMO,SUBR,'sxm')
      CALL CHECK3(SYMO,IM,JM,LMO,SUBR,'sym')
      CALL CHECK3(SZMO,IM,JM,LMO,SUBR,'szm')
      CALL CHECK3(UO  ,IM,JM,LMO,SUBR,'uo ')
      CALL CHECK3(VO  ,IM,JM,LMO,SUBR,'vo ')
#ifdef TRACERS_OCEAN
      CALL CHECK3(TRMO,IM,JM,LMO*NTM,SUBR,'trm')
      CALL CHECK3(TXMO,IM,JM,LMO*NTM,SUBR,'txm')
      CALL CHECK3(TYMO,IM,JM,LMO*NTM,SUBR,'tym')
      CALL CHECK3(TZMO,IM,JM,LMO*NTM,SUBR,'tzm')
#endif

C**** Check for varaibles out of bounds
      QCHECKO=.FALSE.
      DO J=2,JM
      DO I=1,IMAXJ(J)
        IF(FOCEAN(I,J).gt.0.) THEN
C**** Check potential specific enthalpy/salinity
          DO L=1,LMM(I,J)
          GO1 = G0M(I,J,L)/(MO(I,J,L)*DXYPO(J))
          SO1 = S0M(I,J,L)/(MO(I,J,L)*DXYPO(J))
          IF(GO1.lt.-10000. .or. GO1.gt.200000.) THEN
            WRITE (6,*) 'After ',SUBR,': I,J,L,GO=',I,J,L,GO1
            IF (GO1.lt.-20000. .or. GO1.gt.200000.) QCHECKO=.TRUE.
          END IF
          IF(SO1.gt.0.045 .or. SO1.lt.0.) THEN
            WRITE (6,*) 'After ',SUBR,': I,J,L,SO=',I,J,L,1d3*SO1
            IF (SO1.gt.0.05 .or. SO1.lt.0.) QCHECKO=.TRUE.
            IF (I.eq.47.and.J.eq.30) QCHECKO=.FALSE.
          END IF
          END DO
C**** Check all ocean currents
          DO L = 1,LMO
            IF(ABS(UO(I,J,L)).gt.7. .or. ABS(VO(I,J,L)).gt.5) THEN
              WRITE (6,*) 'After ',SUBR,': I,J,L,UO,VO=',I,J,L,UO(I,J,L)
     *             ,VO(I,J,L)
              QCHECKO=.TRUE.
            END IF
          END DO
C**** Check first layer ocean mass
          IF(MO(I,J,1).lt.2000. .or. MO(I,J,1).gt.20000.) THEN
            IF (I.eq.47.and.(J.eq.33.or.J.eq.34)) GOTO 230 ! not Caspian
            WRITE (6,*) 'After ',SUBR,': I,J,MO=',I,J,MO(I,J,1)
            QCHECKO=.TRUE.
          END IF
C**** Check ocean salinity in each eighth box for the first layer
 230      SALIM = .043
          IF(I.eq.47 .and. J.eq.30) GOTO 240 ! not Persian Gulf   !.048
          IF(.5*ABS(SXMO(I,J,1))+.5*ABS(SYMO(I,J,1))+BYRT3*ABS(SZMO(I,J
     *         ,1)).lt.S0M(I,J,1) .and.(.5*ABS(SXMO(I,J,1))+.5
     *         *ABS(SYMO(I,J,1))+BYRT3*ABS(SZMO(I,J,1))+S0M(I,J,1))
     *         /(MO(I,J,1)*DXYPO(J)).lt.SALIM)  GO TO 240
          WRITE (6,*) 'After ',SUBR,': I,J,S0,SX,SY,SZ=',I,J,
     *         S0M(I,J,1)/(MO(I,J,1)*DXYPO(J)),SXMO(I,J,1)/(MO(I,J,1)
     *         *DXYPO(J)),SYMO(I,J,1)/(MO(I,J,1)*DXYPO(J)),SZMO(I,J,1)
     *         /(MO(I,J,1)*DXYPO(J))
          QCHECKO=.TRUE.
 240      CONTINUE
        END IF
      END DO
      END DO

#ifdef TRACERS_OCEAN
      do n=1,ntm
C**** Check for negative tracers
        if (t_qlimit(n)) then
        do l=1,lmo
        do j=1,jm
        do i=1,imaxj(j)
          if (l.le.lmm(i,j)) then
            if (trmo(i,j,l,n).lt.0) then
              write(6,*) "Neg Tracer in ocean after ",subr,i,j,l,
     *             trname(n),trmo(i,j,l,n)
              QCHECKO=.true.
            end if
          end if
        end do
        end do
        end do
        end if
C**** Check conservation of water tracers in ocean
        if (trname(n).eq.'Water') then
          errmax = 0. ; imax=1 ; jmax=1 ; lmax=1
          do l=1,lmo
          do j=1,jm
          do i=1,imaxj(j)
            if (l.le.lmm(i,j)) then
              relerr=max(
     *             abs(trmo(i,j,l,n)-mo(i,j,l)*dxypo(j)+s0m(i,j,l)),
     *             abs(txmo(i,j,l,n)+sxmo(i,j,l)),
     *             abs(tymo(i,j,l,n)+symo(i,j,l)),
     *             abs(tzmo(i,j,l,n)+szmo(i,j,l)))/
     *             (mo(i,j,l)*dxypo(j)-s0m(i,j,l))
            if (relerr.gt.errmax) then
              imax=i ; jmax=j ; lmax=l ; errmax=relerr
            end if
            end if
          end do
          end do
          end do
          print*,"Relative error in ocean fresh water mass after ",subr
     *         ,":",imax,jmax,lmax,errmax,trmo(imax,jmax,lmax,n),mo(imax
     *         ,jmax,lmax)*dxypo(jmax)-s0m(imax,jmax,lmax),txmo(imax
     *         ,jmax,lmax,n),-sxmo(imax,jmax,lmax),tymo(imax,jmax,lmax,n
     *         ),-symo(imax,jmax ,lmax),tzmo(imax,jmax,lmax,n),
     *         -szmo(imax,jmax,lmax)
        end if
      end do
#endif

      IF (QCHECKO)
     &     call stop_model("QCHECKO: Ocean Variables out of bounds",255)

      END IF
C****
      CALL CHECKOST(SUBR)
C****
      END SUBROUTINE CHECKO

      SUBROUTINE conserv_OKE(OKE)
!@sum  conserv_OKE calculates zonal ocean kinetic energy
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : shw
      USE OCEAN, only : im,jm,lmo,fim,imaxj,focean,mo,uo,vo,lmm
      IMPLICIT NONE
!@var OKE zonal ocean kinetic energy per unit area (J/m**2)
      REAL*8, DIMENSION(JM) :: OKE
      INTEGER I,J,L,IP1
      REAL*8 OKEIN

      OKE=0.
      DO J=2,JM-1
        I=IM
        DO IP1=1,IM
          DO L=1,LMM(I,J)
            OKE(J) = OKE(J)+((MO(I,J,L)+MO(IP1,J,L))*UO(I,J,L)*UO(I,J,L)
     *      + MO(I,J,L)*(VO(I,J-1,L)*VO(I,J-1,L) + VO(I,J,L)*VO(I,J,L)))
          END DO
          I=IP1
        END DO
        OKE(J)=OKE(J)*0.25
      END DO
      DO L=1,LMO
        OKEIN = UO(1,JM,L)*UO(1,JM,L)*IM
        DO I=1,IM
          OKEIN  = OKEIN  + VO(I,JM-1,L)*VO(I,JM-1,L)
        END DO
        OKE(JM)= OKE(JM)+ OKEIN*MO(1,JM,L)
      END DO
      OKE(JM)= OKE(JM)*0.5
C****
      RETURN
      END SUBROUTINE conserv_OKE

      SUBROUTINE conserv_OCE(OCEANE)
!@sum  conserv_OCE calculates zonal ocean potential enthalpy(atmos grid)
!@auth Gavin Schmidt
!@ver  1.0
      USE GEOM, only : bydxyp
      USE OCEAN, only : im,jm,fim,imaxj,focean,g0m,lmm
      USE STRAITS, only : nmst,jst,g0mst
      IMPLICIT NONE
!@var OCEANE zonal ocean potential enthalpy (J/m^2)
      REAL*8, DIMENSION(JM) :: OCEANE
      INTEGER I,J,L,N

      DO J=1,JM
        OCEANE(J)=0
        DO I=1,IMAXJ(J)
          DO L=1,LMM(I,J)
            OCEANE(J) = OCEANE(J) + G0M(I,J,L)*FOCEAN(I,J)*BYDXYP(J)
          END DO
        END DO
      END DO
      OCEANE(1) =FIM*OCEANE(1)
      OCEANE(JM)=FIM*OCEANE(JM)
C**** include straits variables
      DO N=1,NMST
        J=JST(N,1)
        OCEANE(J)=OCEANE(J)+SUM(G0MST(:,N))*BYDXYP(J)
      END DO
C****
      RETURN
      END SUBROUTINE conserv_OCE

      SUBROUTINE conserv_OMS(OMASS)
!@sum  conserv_OMS calculates zonal ocean mass (on atmos grid)
!@auth Gavin Schmidt
!@ver  1.0
      USE GEOM, only : bydxyp
      USE OCEAN, only : im,jm,fim,imaxj,focean,mo,g0m,lmm,dxypo
      USE STRAITS, only : nmst,jst,mmst
      IMPLICIT NONE
!@var OMASS zonal ocean mass (kg/m^2)
      REAL*8, DIMENSION(JM) :: OMASS,OMSSV
      COMMON /OCCONS/OMSSV
      INTEGER I,J,L,N

      DO J=1,JM
        OMASS(J)=0
        DO I=1,IMAXJ(J)
          DO L=1,LMM(I,J)
            OMASS(J) = OMASS(J) + MO(I,J,L)*FOCEAN(I,J)
          END DO
        END DO
        OMASS(J)=OMASS(J)*DXYPO(J)*BYDXYP(J)
      END DO
      OMASS(1) =FIM*OMASS(1)
      OMASS(JM)=FIM*OMASS(JM)
C**** include straits variables
      DO N=1,NMST
        J=JST(N,1)
        OMASS(J)=OMASS(J)+SUM(MMST(:,N))*BYDXYP(J)
      END DO
C**** save mass for AM calculation
      OMSSV=OMASS
C****
      RETURN
      END SUBROUTINE conserv_OMS

      SUBROUTINE conserv_OAM(OAM)
!@sum  conserv_OAM calculates zonal ocean angular momentum
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : radius,omega
      USE OCEAN, only : im,jm,fim,imaxj,focean,mo,uo,cosvo,cospo,lmu
      IMPLICIT NONE
!@var OAM ocean angular momentum divided by area (kg/s)
      REAL*8, DIMENSION(JM) :: OAM,OMSSV
      COMMON /OCCONS/OMSSV
      INTEGER I,J,L,IP1
      REAL*8 UMIL

      OAM=0
      DO J=2,JM-1
        UMIL = 0.
        I=IM
        DO IP1=1,IM
          DO L=1,LMU(I,J)
            UMIL = UMIL + UO(I,J,L)*(MO(I,J,L)+MO(IP1,J,L))
          END DO
          I=IP1
        END DO
        OAM(J) = UMIL*COSPO(J) + OMSSV(J)*RADIUS*OMEGA*(COSVO(J-1)
     *       *COSVO(J-1)+COSVO(J)*COSVO(J))
        OAM(J)=0.5*RADIUS*OAM(J)
      END DO
      UMIL = 0.
      DO L=1,LMU(1,JM)
        UMIL = UMIL + UO(1,JM,L)*MO(1,JM,L)
      END DO
      OAM(JM) = UMIL*COSPO(JM)*IM*2. +OMSSV(JM)*RADIUS*OMEGA
     *     *COSVO(JM-1)*COSVO(JM-1)
      OAM(JM)=0.5*RADIUS*OAM(JM)
C****
      RETURN
      END SUBROUTINE conserv_OAM

      SUBROUTINE conserv_OSL(OSALT)
!@sum  conserv_OSL calculates zonal ocean salt on atmos grid
!@auth Gavin Schmidt
!@ver  1.0
      USE GEOM, only : bydxyp
      USE OCEAN, only : im,jm,fim,imaxj,focean,mo,s0m,lmm
      USE STRAITS, only : nmst,jst,s0mst
      IMPLICIT NONE
!@var OSALT zonal ocean salt (kg/m^2)
      REAL*8, DIMENSION(JM) :: OSALT
      INTEGER I,J,L,N

      DO J=1,JM
        OSALT(J)=0
        DO I=1,IMAXJ(J)
          DO L=1,LMM(I,J)
            OSALT(J) = OSALT(J) + S0M(I,J,L)*FOCEAN(I,J)*BYDXYP(J)
          END DO
        END DO
      END DO
      OSALT(1) =FIM*OSALT(1)
      OSALT(JM)=FIM*OSALT(JM)
C**** include straits variables
      DO N=1,NMST
        J=JST(N,1)
        OSALT(J)=OSALT(J)+SUM(S0MST(:,N))*BYDXYP(J)
      END DO
C****
      RETURN
      END SUBROUTINE conserv_OSL

      SUBROUTINE OVTOM (MM,UM,VM)
!@sum  OVTOM converts rho/u/v to mass/momentum in mass units
!@auth Gary Russell
!@ver  1.0
C**** Input:  M (kg/m**2), U (m/s), V (m/s)
      USE OCEAN, only : im,jm,lmo,dxyp=>dxypo,mo,uo,vo
      IMPLICIT NONE
C**** Output: MM (kg), UM (kg*m/s), VM (kg*m/s)
      REAL*8, INTENT(OUT) :: MM(IM,JM,LMO),UM(IM,JM,LMO),VM(IM,JM,LMO)
      INTEGER I,J,L,IP1

      DO 70 L=1,LMO
C**** Convert density to mass
      DO 10 J=1,JM
      DO 10 I=1,IM
      MM(I,J,L) = MO(I,J,L)*DXYP(J)
   10 CONTINUE
C**** Convert U velocity to U momentum
      I=IM
      DO 40 J=2,JM-1
      DO 40 IP1=1,IM
      UM(I, J,L) = UO(I, J,L)*(MM(I,J,L)+MM(IP1,J,L))*5d-1
   40 I=IP1
C     UM(IM,1,L) = UO(IM,1,L)*MM(IM,1,L)*IM
      UM(1,JM,L) = UO(1,JM,L)*MM(1,JM,L)*IM
C**** Convert V velocity to V momentum
      DO 60 I=1,IM
C     VM(I, 1  ,L) = VO(I, 1  ,L)*(MM(1, 1,L)+5d-1*MM(I, 2  ,L))
      VM(I,JM-1,L) = VO(I,JM-1,L)*(MM(1,JM,L)+5d-1*MM(I,JM-1,L))
      DO 60 J=2,JM-2
   60 VM(I,J,L) = VO(I,J,L)*(MM(I,J,L)+MM(I,J+1,L))*5d-1
   70 CONTINUE
      RETURN
      END SUBROUTINE OVTOM

      SUBROUTINE OMTOV (MM,UM,VM)
!@sum  OMTOV converts mass/momentum in mass units to rho/u/v
!@auth Gary Russell
!@ver  1.0
C**** Output: M (kg/m**2), U (m/s), V (m/s)
      USE OCEAN, only : im,jm,lmo,mo,uo,vo,lmm,lmu,lmv,
     *     dxyp=>dxypo,bydxyp=>bydxypo
      IMPLICIT NONE
C**** Input:  MM (kg), UM (kg*m/s), VM (kg*m/s)
      REAL*8, INTENT(IN) :: MM(IM,JM,LMO),UM(IM,JM,LMO),VM(IM,JM,LMO)
      INTEGER I,J,L
      REAL*8 BYD

!$OMP PARALLEL DO  PRIVATE(I,J,L,BYD)
      DO L=1,LMO
C**** Convert mass to density
        DO J=1,JM
          BYD = BYDXYP(J)
          DO I=1,IM
            IF (L.LE.LMM(I,J)) MO(I,J,L) = MM(I,J,L)*BYD
          END DO
        END DO
C**** Convert V momentum to V velocity
        DO J=2,JM-2
          DO I=1,IM
            IF (L.LE.LMV(I,J)) VO(I,J,L) = VM(I,J,L)*2d0/
     *         (MM(I,J,L) + MM(I,J+1,L))
          END DO
        END DO
        DO I=1,IM
          IF (L.LE.LMV(I,JM-1)) VO(I,JM-1,L) = VM(I,JM-1,L)/
     *       (MM(I,JM,L)+5d-1*MM(I,JM-1,L))
        END DO
C**** Convert U momentum to U velocity
        DO J=2,JM-1
          DO I=1,IM-1
            IF (L.LE.LMU(I,J)) UO(I,J,L) = UM(I,J,L)*2d0/
     *         (MM(I,J,L) + MM(I+1,J,L))
          END DO
          IF (L.LE.LMU(IM,J)) UO(IM,J,L) = UM(IM,J,L)*2d0/
     *         (MM(IM,J,L) + MM(1,J,L))
        END DO
      END DO
!$OMP END PARALLEL DO
C     UO(IM,1,L) = UM(IM,1,L)/(MM(IM,1,L)*IM)
      DO L=1,LMU(1,JM)
        UO(1,JM,L) = UM(1,JM,L)/(MM(1,JM,L)*IM)
        DO I=2,IM
           UO(I,JM,L) = UO(1,JM,L)
        END DO
      END DO
      RETURN
      END SUBROUTINE OMTOV

      SUBROUTINE OFLUX (NS,MM,QSAVE)
!@sum OFLUX calculates the fluid fluxes
!@auth Gary Russell
!@ver  1.0
C**** Input:  M (kg/m^2), U (m/s), V (m/s)
C**** Output: MU (kg/s) = DY * U * M
C****         MV (kg/s) = DX * V * M
C****         MW (kg/s) = MW(L-1) + CONV-SCONV*DSIGO + (MM-SMM*DSIGO)
C****       CONV (kg/s) = MU-MU+MV-MV
C****
      USE OCEAN, only : im,jm,lmo,byim,uo,vo,dts,dto,lmm,lmu,lmv
     *     ,dxyp=>dxypo,dxp=>dxpo,dyp=>dypo,dxv=>dxvo,dyv=>dyvo,dsigo
     *     ,sigeo,mo
      USE OCEAN_DYN, only : smu,smv,smw,mu,mv,mw,conv
      IMPLICIT NONE
      REAL*8, PARAMETER :: CHI=.125d0, BY12CHI= 2d0/(1d0-2d0*CHI)
      LOGICAL, INTENT(IN) :: QSAVE
      INTEGER, INTENT(IN) :: NS
      REAL*8, INTENT(IN) :: MM(IM,JM,LMO)
      INTEGER I,J,L,LMIJ
      REAL*8, DIMENSION(IM) :: DMVS,DMVN
      REAL*8 BYNSDTO,MUS,MUN,MVS,MVN,ADMVS,ADMVN,SCONV,SMM,BYSIGEO
C****
      BYNSDTO=1d0/(NS*DTO)
C**** Compute fluid fluxes for the C grid
C****
C**** Smooth the West-East velocity near the poles
!$OMP PARALLEL DO  PRIVATE(I,L)
      DO 110 L=1,LMO
      DO 110 I=IM+1,IM*(JM-1)
  110   MU(I,1,L) = UO(I,1,L)
!$OMP END PARALLEL DO
      CALL OPFIL (MU,LMO)
C**** Compute MU, the West-East mass flux, at non-polar points
!$OMP PARALLEL DO  PRIVATE(I,J,L, MUS,MUN,MVS,MVN, DMVS,DMVN,
!$OMP&                     ADMVS,ADMVN)
      DO 220 L=1,LMO
      DO J=2,JM-1
        DO I=1,IM-1
          MU(I,J,L) = 5d-1*DYP(J)*MU(I,J,L)*(MO(I,J,L)+MO(I+1,J,L))
        END DO
        MU(IM,J,L) = 5d-1*DYP(J)*MU(IM,J,L)*(MO(IM,J,L)+MO(1,J,L))
      END DO
C**** Compute MV, the South-North mass flux
      DO 130 J=1,JM-1
      DO 130 I=1,IM
  130   MV(I,J,L) = 5d-1*DXV(J)*VO(I,J,L)*(MO(I,J,L)+MO(I,J+1,L))
C**** Compute MU*2/(1-2*CHI) at the poles
      MUS = DYP( 1)*UO(1, 1,L)*MO(1, 1,L)
      MUN = DYP(JM)*UO(1,JM,L)*MO(1,JM,L)
      MVS = 0.
      MVN = 0.
      DO 140 I=1,IM
        MVS = MVS + MV(I, 1  ,L)
  140   MVN = MVN + MV(I,JM-1,L)
      MVS = MVS*BYIM
      MVN = MVN*BYIM
      DMVS(1) = 0.
      DMVN(1) = 0.
      DO 150 I=2,IM
        DMVS(I) = DMVS(I-1) + (MV(I, 1  ,L)-MVS)
  150   DMVN(I) = DMVN(I-1) + (MV(I,JM-1,L)-MVN)
      ADMVS = 0.
      ADMVN = 0.
      DO 160 I=1,IM
        ADMVS = ADMVS + DMVS(I)
  160   ADMVN = ADMVN + DMVN(I)
      ADMVS = ADMVS*BYIM
      ADMVN = ADMVN*BYIM
      DO 170 I=1,IM
        MU(I, 1,L) = (ADMVS-DMVS(I) + MUS)*BY12CHI
  170   MU(I,JM,L) = (DMVN(I)-ADMVN + MUN)*BY12CHI
C****
C**** Compute horizontal fluid convergence
C****
      DO J=2,JM-1
        CONV(1,J,L) = MU(IM,J,L)-MU(1,J,L) + (MV(1,J-1,L)-MV(1,J,L))
        DO I=2,IM
          CONV(I,J,L) = MU(I-1,J,L)-MU(I,J,L) + (MV(I,J-1,L)-MV(I,J,L))
        END DO
      END DO
      CONV(IM,1,L) = -MVS
      CONV(1,JM,L) =  MVN
  220 CONTINUE
!$OMP END PARALLEL DO
C****
C**** Compute vertically integrated column convergence and mass
C****
!$OMP PARALLEL DO  PRIVATE(I,L,LMIJ,SCONV,SMM,BYSIGEO)
      DO 440 I=IM,IM*(JM-1)+1
      LMIJ=1
      IF(LMM(I,1).LE.1)  GO TO 420
      LMIJ=LMM(I,1)
      SCONV = 0.
      SMM   = 0.
      DO 310 L=1,LMIJ
        SCONV = SCONV + CONV(I,1,L)
  310   SMM   = SMM   +   MM(I,1,L)
      BYSIGEO=1d0/SIGEO(LMIJ)
      SCONV = SCONV*BYSIGEO
      SMM   = SMM  *BYSIGEO
C****
C**** Compute MW, the downward fluid flux
C****
      MW(I,1,1) = CONV(I,1,1)-SCONV*DSIGO(1) +
     +           (  MM(I,1,1)-  SMM*DSIGO(1))*BYNSDTO
      DO 410 L=2,LMIJ-1
  410   MW(I,1,L) = CONV(I,1,L)-SCONV*DSIGO(L) + MW(I,1,L-1) +
     +           (  MM(I,1,L)-  SMM*DSIGO(L))*BYNSDTO
  420 DO 430 L=LMIJ,LMO-1
  430   MW(I,1,L) = 0.
  440 CONTINUE
!$OMP END PARALLEL DO
C****
C**** Sum mass fluxes to be used for advection of tracers
C****
      IF (QSAVE) THEN
!$OMP PARALLEL DO  PRIVATE(L)
        DO L=1,LMO
          SMU(:,:,L) = SMU(:,:,L) + MU(:,:,L)
          SMV(:,:,L) = SMV(:,:,L) + MV(:,:,L)
          SMW(:,:,L) = SMW(:,:,L) + MW(:,:,L)
        END DO
!$OMP END PARALLEL DO
      END IF

      RETURN
      END SUBROUTINE OFLUX

      SUBROUTINE OPFIL (X,LMAX)
!@sum  OPFIL smoothes X in zonal direction
!@auth Gary Russell
!@ver  1.0
C****
C**** OPFIL smoothes X in zonal direction by reducing coefficients
C**** of its Fourier series for high wave numbers near the poles.
C**** INDM for polar ocean filter set to work with AVR4X5LD.Z12.
C****            Now AVR4X5LD.Z12.gas1
      USE CONSTANT, only : twopi
      USE OCEAN, only : im,jm,lmo,dxyp=>dxypo,dxp=>dxpo,dyp=>dypo,dxv
     *     =>dxvo,dyv=>dyvo,j40s
      USE FILEMANAGER, only : openunit, closeunit
      IMPLICIT NONE
      REAL*8, PARAMETER :: DLON=TWOPI/IM
      INTEGER, INTENT(IN) :: LMAX
      REAL*8, INTENT(INOUT) :: X(IM,JM,LMO)
      INTEGER, PARAMETER :: NMAX=IM/2
      INTEGER, SAVE :: NSEGM,INDM, JXMAX, JNOF
      INTEGER*2, SAVE, DIMENSION(:,:), ALLOCATABLE :: NSEG
      INTEGER*4, SAVE, DIMENSION(:,:), ALLOCATABLE :: INDX
      INTEGER*2, SAVE, DIMENSION(:,:,:), ALLOCATABLE :: IMIN,ILEN
      REAL*4, SAVE, DIMENSION(:), ALLOCATABLE :: REDUCO
      REAL*8, SAVE :: AVCOS(IM,NMAX),AVSIN(IM,NMAX)
      INTEGER, SAVE :: IFIRST=1
      REAL*8  AN(0:NMAX),BN(0:NMAX), Y(IM)
      INTEGER I,J,L,N,JA,JX,NS,I0,IL,IND,K,iu_AVR,IN
      REAL*8 :: REDUC,DRATM,SM
      CHARACTER*80 TITLE
C****
      IF (IFIRST.EQ.1) THEN
      IFIRST=0
      JXMAX=2*J40S-1
      JNOF=JM-2*J40S
C**** Calculate  cos(TWOPI*N*I/IM)  and  sin(TWOPI*N*I/IM)
      DO 10 I=1,NMAX
      AVCOS(I,1) = COS(DLON*I)
      AVSIN(I,1) = SIN(DLON*I)
      AVCOS(I+NMAX,1) = -AVCOS(I,1)
   10 AVSIN(I+NMAX,1) = -AVSIN(I,1)
      DO 20 N=2,NMAX
      DO 20 I=1,IM
      IN = I*N-((I*N-1)/IM)*IM
      AVCOS(I,N) = AVCOS(IN,1)
   20 AVSIN(I,N) = AVSIN(IN,1)
C**** Read in reduction contribution matrices from disk
      call openunit("AVR",iu_AVR,.TRUE.,.TRUE.)
C**** read in size of arrays and allocate space
      READ  (iu_AVR) TITLE,NSEGM,INDM
      ALLOCATE(IMIN(NSEGM,LMO,4:JXMAX))
      ALLOCATE(ILEN(NSEGM,LMO,4:JXMAX))
      ALLOCATE(NSEG(LMO,4:JXMAX))
      ALLOCATE(INDX(IM,2:J40S))
      ALLOCATE(REDUCO(INDM))
C**** read in arrays
      READ  (iu_AVR) TITLE,NSEG,IMIN,ILEN,INDX,REDUCO
      call closeunit (iu_AVR)
      WRITE (6,*) 'Read on unit ',iu_AVR,': ',TITLE
      END IF
C****
C**** Loop over J and L.  J = latitude, JA = absolute latitude.
C****
!$OMP PARALLEL DO  PRIVATE(I,J,JA,JX,K,L,N,NS,I0,IL,IND,Y,REDUC,
!$OMP&                     AN,BN,DRATM,SM)
  100 DO 330 JX=4,JXMAX
      J =JX
      JA=JX
      IF(JX.GT.J40S)  J =JNOF+JX
      IF(JX.GT.J40S)  JA=2*J40S+1-JX
      DO 330 L=1,LMAX
      IF(ILEN(1,L,JX).GE.IM)  GO TO 300
C****
C**** Land boxes exist at this latitude and layer, loop over ocean
C**** segments.
C****
      DO 270 NS=1,NSEG(L,JX)
      I0 = IMIN(NS,L,JX)
      IL = ILEN(NS,L,JX)
      IND= INDX(IL,JA)
      IF(IL+I0.GT.IM)  GO TO 200
C**** Ocean segment does not wrap around the IDL.
C**** Copy X to temporary array and filter X in place.
      DO 110 I=1+I0,IL+I0
  110 Y(I-I0) = X(I,J,L)
      DO 140 I=1+I0,IL+I0
      REDUC = 0.
      DO 130 K=1,IL
      IND=IND+1
  130 REDUC = REDUC + REDUCO(IND)*Y(K)
  140 X(I,J,L) = X(I,J,L) - REDUC
      GO TO 270
C**** Ocean segment wraps around the IDL.
C**** Copy X to temporary array and filter X in place.
  200 DO 210 I=1+I0,IM
  210 Y(I-I0) = X(I,J,L)
      DO 220 I=1,IL+I0-IM
  220 Y(I-I0+IM) = X(I,J,L)
      DO 240 I=1+I0,IM
      REDUC = 0.
      DO 230 K=1,IL
      IND=IND+1
  230 REDUC = REDUC + REDUCO(IND)*Y(K)
  240 X(I,J,L) = X(I,J,L) - REDUC
      DO 260 I=1,IL+I0-IM
      REDUC = 0.
      DO 250 K=1,IL
      IND=IND+1
  250 REDUC = REDUC + REDUCO(IND)*Y(K)
  260 X(I,J,L) = X(I,J,L) - REDUC
  270 CONTINUE
      GO TO 330
C****
C**** No land boxes at this latitude and layer,
C**** perform standard polar filter
C****
  300 CALL FFT (X(1,J,L),AN,BN)
      DRATM = NMAX*DXP(J)/DYP(3)
      DO N=NMAX,INT(DRATM)+1,-1   ! guarantees that SM > 0
        SM = 1. - DRATM/N
        DO I=1,IM
          X(I,J,L) = X(I,J,L) - SM*(AN(N)*AVCOS(I,N)+BN(N)*AVSIN(I,N))
        END DO
      END DO
  330 CONTINUE
!$OMP END PARALLEL DO
      RETURN
      END SUBROUTINE OPFIL

      SUBROUTINE OADVM (MM2,MM0,DT)
!@sum  OADVM calculates the updated column mass
!@auth Gary Russell
!@ver  1.0
C****
C**** Input:  MM0 (kg), CONV (kg/s), MW (kg/s), DT (s)
C**** Output: MM2 (kg) = MM0 + DT*DM*DSIGO
C****
      USE OCEAN, only : im,jm,lmo,lmm
      USE OCEAN_DYN, only : mu,mv,mw,conv
      IMPLICIT NONE

      REAL*8, INTENT(OUT) :: MM2(IM,JM,LMO)
      REAL*8, INTENT(IN)  :: MM0(IM,JM,LMO)
      REAL*8, INTENT(IN)  :: DT
      INTEGER I,J,L,LMIJ
C****
C**** Compute the new mass MM2
C****
!$OMP PARALLEL DO  PRIVATE(I)
      DO 40 I=IM,IM*(JM-1)+1
      LMIJ=LMM(I,1)
      IF(LMIJ-1)  40,10,20
   10 MM2(I,1,1) = MM0(I,1,1) + DT*CONV(I,1,1)
      GO TO 40
   20 MM2(I,1,1) = MM0(I,1,1) + DT*(CONV(I,1,1)-MW(I,1,1))
      DO 30 L=2,LMIJ-1
   30 MM2(I,1,L) = MM0(I,1,L) + DT*(CONV(I,1,L)+(MW(I,1,L-1)-MW(I,1,L)))
      MM2(I,1,LMIJ) = MM0(I,1,LMIJ) + DT*(CONV(I,1,LMIJ)+MW(I,1,LMIJ-1))
   40 CONTINUE
!$OMP END PARALLEL DO
C**** Fill in values at the poles
      DO 80 L=1,LMO
      DO 80 I=1,IM
C     MM2(I, 1,L) = MM2(IM,1,L)
   80 MM2(I,JM,L) = MM2(1,JM,L)
      RETURN
      END SUBROUTINE OADVM

      SUBROUTINE OADVV (UM2,VM2,UM0,VM0,DT1)
!@sum  OADVV advects oceanic momentum (with coriolis force)
!@auth Gary Russell
!@ver  1.0
C****
C**** Input:  MO (kg/m**2), UO (m/s), VO (m/s) = from odd solution
C****          MU (kg/s), MV (kg/s), MW (kg/s) = fluid fluxes
C**** Output: UM2 (kg*m/s) = UM0 + DT*(MU*U-MU*U+MV*U-MV*U+M*CM*V)
C****         VM2 (kg*m/s) = VM0 + DT*(MU*V-MU*V+MV*V-MV*V-M*CM*U)
C****
      USE CONSTANT, only : twopi,omega,sday,radius
      USE OCEAN, only : im,jm,lmo,mo,uo,vo,
     *     dxyp=>dxypo,dxp=>dxpo,dyp=>dypo,dxv=>dxvo,dyv=>dyvo,
     *     cosv=>cosvo,cosp=>cospo
      USE OCEAN_DYN, only : mu,mv,mw
      IMPLICIT NONE
      REAL*8, PARAMETER :: CHI=.125d0
      REAL*8, INTENT(OUT), DIMENSION(IM,JM,LMO) :: UM2,VM2
      REAL*8, INTENT(IN),  DIMENSION(IM,JM,LMO) :: UM0,VM0
      REAL*8, INTENT(IN) :: DT1
      REAL*8, SAVE, DIMENSION(JM) :: DCOSP,TANP
      REAL*8, DIMENSION(IM,JM,LMO) :: DUM,DVM
      REAL*8, DIMENSION(0:JM) :: GX,GXUY,DGDUDN,DGDUUP
      REAL*8, DIMENSION(IM) :: UDCOSY,UYUTAY
      REAL*8, DIMENSION(IM,JM,LMO-1) :: FLUXA
      LOGICAL*4 :: QFIRST = .TRUE.
      INTEGER I,J,L
      REAL*8 DT2,DT4,DTC2,DTD4,UMUX,VMVY,VMUNE,VMUSE,UMVX,VMUY,UMUNE
     *     ,UMUNW,DUMS,DUMN,UMVNE,UMVNW
C****
      IF (QFIRST) THEN
        QFIRST = .FALSE.
        DO J=1,JM
          DCOSP(J) = COSV(J-1)-COSV(J)
          TANP(J)  = DCOSP(J)/COSP(J)
        END DO
      END IF
C****
      DT2  = DT1*5d-1
      DT4  = DT1*2.5d-1
      DTC2 = DT1*CHI*5d-1
      DTD4 = DT1*(1d0-2d0*CHI)*2.5d-1
C**** Zero out momentum changes
      DUM = 0.
      DVM = 0.
C****
C**** Horizontal advection of momentum
C****
!$OMP PARALLEL DO  PRIVATE(I,J,L, UMUX,VMVY,VMUNE,VMUSE,
!$OMP&  UMVX,UMVNE,UMVNW,VMUY, GX,GXUY,DGDUDN,DGDUUP,UDCOSY,UYUTAY)
      DO 480 L=1,LMO
        DO J=2,JM-1
C**** Contribution of West-East flux to U wind
C**** Contribution of South-North and corner fluxes to V wind
C**** (split into first then rest)
          UMUX = DT4*(MU(IM,J,L)+MU(1,J,L))*(UO(IM,J,L)+UO(1,J,L))
          DUM(IM,J,L) = DUM(IM,J,L) - UMUX
          DUM(1 ,J,L) = DUM(1 ,J,L) + UMUX
          VMVY  = DT4 * (MV(IM,J-1,L)+MV(IM,J,L))*
     *                  (VO(IM,J-1,L)+VO(IM,J,L))
          VMUNE = DTC2* MU(IM,J,L) * (VO(IM,J-1,L) + VO(1,J  ,L))
          VMUSE = DTC2* MU(IM,J,L) * (VO(IM,J  ,L) + VO(1,J-1,L))
          DVM(IM,J-1,L) = DVM(IM,J-1,L) - (VMUNE+VMVY)
          DVM(1 ,J-1,L) = DVM(1 ,J-1,L) +  VMUSE
          DVM(IM,J  ,L) = DVM(IM,J  ,L) - (VMUSE-VMVY)
          DVM(1 ,J  ,L) = DVM(1 ,J  ,L) +  VMUNE
          DO I=2,IM     ! not first
            UMUX = DT4*(MU(I-1,J,L)+MU(I,J,L))*(UO(I-1,J,L)+UO(I,J,L))
            DUM(I  ,J,L) = DUM(I  ,J,L) + UMUX
            DUM(I-1,J,L) = DUM(I-1,J,L) - UMUX
            VMVY  = DT4* (MV(I-1,J-1,L)+MV(I-1,J,L))*
     *                   (VO(I-1,J-1,L)+VO(I-1,J,L))
            VMUNE = DTC2* MU(I-1,J,L) * (VO(I-1,J-1,L) + VO(I,J  ,L))
            VMUSE = DTC2* MU(I-1,J,L) * (VO(I-1,J  ,L) + VO(I,J-1,L))
            DVM(I  ,J-1,L) = DVM(I  ,J-1,L) +  VMUSE
            DVM(I-1,J-1,L) = DVM(I-1,J-1,L) - (VMUNE+VMVY)
            DVM(I  ,J  ,L) = DVM(I  ,J  ,L) +  VMUNE
            DVM(I-1,J  ,L) = DVM(I-1,J  ,L) - (VMUSE-VMVY)
          END DO
        END DO
C**** Contribution of South-North and corner fluxes to U wind
C**** Contribution of West-East flux to V wind
C**** (split into first then rest)
      DO J=1,JM-1
        UMVX  = DTD4*(MV(IM,J,L)+MV(1,J,L))*(UO(IM,J,L)+UO(IM,J+1,L))
        UMVNE = DTC2* MV(1,J,L)*(UO(IM,J,L)+UO(1 ,J+1,L))
        UMVNW = DTC2* MV(1,J,L)*(UO(1 ,J,L)+UO(IM,J+1,L))
        DUM(IM,J  ,L) = DUM(IM,J  ,L) - (UMVNE+UMVX)
        DUM(1 ,J  ,L) = DUM(1 ,J  ,L) -  UMVNW
        DUM(IM,J+1,L) = DUM(IM,J+1,L) + (UMVNW+UMVX)
        DUM(1 ,J+1,L) = DUM(1 ,J+1,L) +  UMVNE
        VMUY = DTD4*(MU(IM,J,L)+MU(IM,J+1,L))*(VO(IM,J,L)+VO(1,J,L))
        DVM(IM,J  ,L) = DVM(IM,J  ,L) - VMUY
        DVM(1 ,J  ,L) = DVM(1 ,J  ,L) + VMUY
        DO I=2,IM
          UMVX=DTD4*(MV(I-1,J,L)+MV(I,J,L))*(UO(I-1,J,L)+UO(I-1,J+1,L))
          UMVNE = DTC2* MV(I,J,L)*(UO(I-1,J,L)+UO(I  ,J+1,L))
          UMVNW = DTC2* MV(I,J,L)*(UO(I  ,J,L)+UO(I-1,J+1,L))
          DUM(I  ,J  ,L) = DUM(I  ,J  ,L) -  UMVNW
          DUM(I-1,J  ,L) = DUM(I-1,J  ,L) - (UMVNE+UMVX)
          DUM(I  ,J+1,L) = DUM(I  ,J+1,L) +  UMVNE
          DUM(I-1,J+1,L) = DUM(I-1,J+1,L) + (UMVNW+UMVX)
          VMUY=DTD4*(MU(I-1,J,L)+MU(I-1,J+1,L))*(VO(I-1,J,L)+VO(I,J,L))
          DVM(I  ,J  ,L) = DVM(I  ,J  ,L) + VMUY
          DVM(I-1,J  ,L) = DVM(I-1,J  ,L) - VMUY
        END DO
      END DO
C**** Coriolis force and metric term
C****
C**** U component
          GX(0) = 0.d0
        GXUY(0) = 0.d0
      DGDUDN(0) = 0.d0
      DGDUUP(0) = 0.d0
          GX(JM) = 0.d0
        GXUY(JM) = 0.d0
      DGDUDN(JM) = 0.d0
      DGDUUP(JM) = 0.d0
C**** First I=1,then I=2,IM-1 then I=IM
      DO J=1,JM-1
        GX(J)   =  MV(1,J,L)+MV(2,J,L)
        GXUY(J) = (MV(1,J,L)+MV(2,J,L))*(UO(1,J,L)+UO(1,J+1,L))
        DGDUDN(J) = MV(1,J,L)*(UO(IM,J  ,L)-UO(1,J  ,L)) -
     *              MV(2,J,L)*(UO(1 ,J  ,L)-UO(2,J  ,L))
        DGDUUP(J) = MV(1,J,L)*(UO(IM,J+1,L)-UO(1,J+1,L)) -
     *              MV(2,J,L)*(UO(1 ,J+1,L)-UO(2,J+1,L))
      END DO
      DO J=1,JM
        DUM(1,J,L) = DUM(1,J,L) +
     *    DT2*(OMEGA*RADIUS*(GX(J-1)+GX(J))*DCOSP(J) + TANP(J)*
     *    (2.5d-1*(GXUY(J-1)+GXUY(J))+5d-1*CHI*(DGDUDN(J-1)+DGDUUP(J))))
      END DO
      DO I=2,IM-1
        DO J=1,JM-1
          GX(J)   =  MV(I,J,L)+MV(I+1,J,L)
          GXUY(J) = (MV(I,J,L)+MV(I+1,J,L))*(UO(I,J,L)+UO(I,J+1,L))
          DGDUDN(J) = MV(I  ,J,L)*(UO(I-1,J  ,L)-UO(I  ,J  ,L)) -
     *                MV(I+1,J,L)*(UO(I  ,J  ,L)-UO(I+1,J  ,L))
          DGDUUP(J) = MV(I  ,J,L)*(UO(I-1,J+1,L)-UO(I  ,J+1,L)) -
     *                MV(I+1,J,L)*(UO(I  ,J+1,L)-UO(I+1,J+1,L))
        END DO
        DO J=1,JM
          DUM(I,J,L) = DUM(I,J,L) +
     *    DT2*(OMEGA*RADIUS*(GX(J-1)+GX(J))*DCOSP(J) + TANP(J)*
     *    (2.5d-1*(GXUY(J-1)+GXUY(J))+5d-1*CHI*(DGDUDN(J-1)+DGDUUP(J))))
        END DO
      END DO
      DO J=1,JM-1
        GX(J)   =  MV(IM,J,L)+MV(1,J,L)
        GXUY(J) = (MV(IM,J,L)+MV(1,J,L))*(UO(IM,J,L)+UO(IM,J+1,L))
        DGDUDN(J) = MV(IM,J,L)*(UO(IM-1,J  ,L)-UO(IM,J  ,L)) -
     *              MV(1 ,J,L)*(UO(IM  ,J  ,L)-UO(1 ,J  ,L))
        DGDUUP(J) = MV(IM,J,L)*(UO(IM-1,J+1,L)-UO(IM,J+1,L)) -
     *              MV(1 ,J,L)*(UO(IM  ,J+1,L)-UO(1 ,J+1,L))
      END DO
      DO J=1,JM
        DUM(IM,J,L) = DUM(IM,J,L) +
     *    DT2*(OMEGA*RADIUS*(GX(J-1)+GX(J))*DCOSP(J) + TANP(J)*
     *    (2.5d-1*(GXUY(J-1)+GXUY(J))+5d-1*CHI*(DGDUDN(J-1)+DGDUUP(J))))
      END DO
C**** V component
      DO J=1,JM-1
        DO I=1,IM
          UDCOSY(I) =  UO(I,J,L)*DCOSP(J) + UO(I,J+1,L)*DCOSP(J+1)
          UYUTAY(I) = (UO(I,J,L)* TANP(J) + UO(I,J+1,L)* TANP(J+1))*
     *            (UO(I,J,L)+UO(I,J+1,L))
        END DO
        DVM(1,J,L) = DVM(1,J,L) - DT4*DXV(J)*(MO(1,J,L)+MO(1,J+1,L))*
     *    (OMEGA*RADIUS*(UDCOSY(IM)+UDCOSY(1)) +
     *    2.5d-1*(UYUTAY(IM)+UYUTAY(1)) - 5d-1*CHI*(TANP(J)+TANP(J+1))*
     *    (UO(IM,J,L)-UO(1,J,L))*(UO(IM,J+1,L)-UO(1,J+1,L)))
        DO I=2,IM
          DVM(I,J,L) = DVM(I,J,L) - DT4*DXV(J)*(MO(I,J,L)+MO(I,J+1,L))*
     *      (OMEGA*RADIUS*(UDCOSY(I-1)+UDCOSY(I)) +
     *      2.5d-1*(UYUTAY(I-1)+UYUTAY(I))-5d-1*CHI*(TANP(J)+TANP(J+1))*
     *      (UO(I-1,J,L)-UO(I,J,L))*(UO(I-1,J+1,L)-UO(I,J+1,L)))
        END DO
      END DO
  480 CONTINUE
!$OMP END PARALLEL DO
C****
C**** Vertical advection of momentum
C****
!$OMP PARALLEL DO  PRIVATE(I,J,L)
      DO L=1,LMO-1
        DO J=2,JM-1
          DO I=1,IM-1
            FLUXA(I,J,L) =
     *           DT4*(MW(I,J,L)+MW(I+1,J,L))*(UO(I,J,L)+UO(I,J,L+1))
          END DO
          FLUXA(IM,J,L) =
     *         DT4*(MW(IM,J,L)+MW(1,J,L))*(UO(IM,J,L)+UO(IM,J,L+1))
        END DO
        FLUXA(1,JM,L) = DT2*MW(1,JM,L)*IM*(UO(1,JM,L)+UO(1,JM,L+1))
      END DO
!$OMP END PARALLEL DO
C
!$OMP PARALLEL DO  PRIVATE(L)
      DO L=1,LMO
        IF(L.EQ.1)  THEN
          DUM(:,2:JM-1,L) = DUM(:,2:JM-1,L) - FLUXA(:,2:JM-1,L)
          DUM(1,JM,L)     = DUM(1,JM,L) - FLUXA(1,JM,L)
        ELSE IF(L.LT.LMO)  THEN
          DUM(:,2:JM-1,L) = DUM(:,2:JM-1,L) + FLUXA(:,2:JM-1,L-1)
          DUM(1,JM,L)     = DUM(1,JM,L) + FLUXA(1,JM,L-1)
          DUM(:,2:JM-1,L) = DUM(:,2:JM-1,L) - FLUXA(:,2:JM-1,L)
          DUM(1,JM,L)     = DUM(1,JM,L) - FLUXA(1,JM,L)
        ELSE IF(L.EQ.LMO)  THEN
          DUM(:,2:JM-1,L) = DUM(:,2:JM-1,L) + FLUXA(:,2:JM-1,L-1)
          DUM(1,JM,L)     = DUM(1,JM,L) + FLUXA(1,JM,L-1)
        END IF
      END DO
!$OMP END PARALLEL DO
C
C**** V component
C
!$OMP PARALLEL DO  PRIVATE(I,J,L)
      DO L=1,LMO-1
        DO J=2,JM-2
          DO I=1,IM
            FLUXA(I,J,L) =
     *           DT4*(MW(I,J,L)+MW(I,J+1,L))*(VO(I,J,L)+VO(I,J,L+1))
          END DO
        END DO
        DO I=1,IM
          FLUXA(I,JM-1,L) = DT4*(MW(I,JM-1,L)+2.*MW(1,JM,L))*
     *         (VO(I,JM-1,L)+VO(I,JM-1,L+1))
        END DO
      END DO
!$OMP END PARALLEL DO
C
!$OMP PARALLEL DO  PRIVATE(L)
      DO L=1,LMO
        IF(L.EQ.1)  THEN
          DVM(:,2:JM-1,L) = DVM(:,2:JM-1,L) - FLUXA(:,2:JM-1,L)
        ELSE IF(L.LT.LMO)  THEN
          DVM(:,2:JM-1,L) = DVM(:,2:JM-1,L) + FLUXA(:,2:JM-1,L-1)
          DVM(:,2:JM-1,L) = DVM(:,2:JM-1,L) - FLUXA(:,2:JM-1,L)
        ELSE IF(L.EQ.LMO)  THEN
          DVM(:,2:JM-1,L) = DVM(:,2:JM-1,L) + FLUXA(:,2:JM-1,L-1)
        END IF
      END DO
!$OMP END PARALLEL DO
C****
C**** Add changes to momentum
C****
!$OMP PARALLEL DO  PRIVATE(I,L,DUMS,DUMN)
      DO L=1,LMO
C**** U component
        DO 610 I=IM+1,IM*(JM-1)
  610     UM2(I,1,L) = UM0(I,1,L) + DUM(I,1,L)
        DUMS = 0.
        DUMN = 0.
        DO 620 I=1,IM
          DUMS = DUMS + DUM(I, 1,L)
  620     DUMN = DUMN + DUM(I,JM,L)
        UM2(IM,1,L) = UM0(IM,1,L) + DUMS
        UM2(1,JM,L) = UM0(1,JM,L) + DUMN
C**** V component
        DO 630 I=1,IM*(JM-1)
  630     VM2(I,1,L) = VM0(I,1,L) + DVM(I,1,L)
      END DO
!$OMP END PARALLEL DO
      RETURN
      END SUBROUTINE OADVV

      SUBROUTINE OPGF (UM,VM,DT1)
!@sum  OPGF adds the pressure gradient force to the momentum
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : grav
      USE OCEAN, only : im,jm,lmo,g0m,gzmo,szmo,mo,hocean,lmm,lmu,lmv
     *     ,focean,dypo,dxvo,opress,ogeoz,imaxj
      USE OCEAN_DYN, only : po,phi,dh,dzgdp,vbar
      IMPLICIT NONE
C****
C**** Input: G0M (J), GZMO, S0M (kg), SZMO, DT (s), MO (kg/m**2)
C**** Output: UM (kg*m/s) = UM - DT*(DH*D(P)+MO*D(PHI))*DYP
C****         VM (kg*m/s) = VM - DT*(DH*D(P)+MO*D(PHI))*DXV
C****
      REAL*8, INTENT(OUT), DIMENSION(IM,JM,LMO) :: UM,VM
      REAL*8, INTENT(IN) :: DT1
      REAL*8, DIMENSION(IM,JM,LMO) :: DUM,DVM
      INTEGER I,J,L,IP1
      REAL*8, SAVE :: PHIE
      REAL*8 DT2
C****
C**** SURFCB  OPRESS  Atmospheric and sea ice pressure (Pa-101325)
C****         OGEOZ  Ocean surface geopotential (m^2/s^2)
C****
      DT2   = DT1*5d-1
C****
C**** Calculate the mass weighted pressure P (Pa),
C**** geopotential PHI (m**2/s**2), and layer thickness DH (m)
C****
!$OMP PARALLEL DO  PRIVATE(I,J,L,PHIE)
      DO J=1,JM
      DO I=1,IMAXJ(J)
        IF(FOCEAN(I,J).GT.0.) THEN
C**** Calculate pressure by integrating from the top down
          PO(I,J,1) = OPRESS(I,J) + MO(I,J,1)*5d-1*GRAV
          DO L=2,LMM(I,J)
            PO(I,J,L) = PO(I,J,L-1) + (MO(I,J,L-1)+MO(I,J,L))*5d-1*GRAV
          END DO
C**** Calculate geopotential by integrating from the bottom up
          PHIE = -HOCEAN(I,J)*GRAV
          DO L=LMM(I,J),1,-1
            PHI(I,J,L) = PHIE + DZGDP(I,J,L)*MO(I,J,L)*GRAV
            DH(I,J,L) =  VBAR(I,J,L)*MO(I,J,L)
            PHIE = PHIE + VBAR(I,J,L)*MO(I,J,L)*GRAV
          END DO
          OGEOZ(I,J) = PHIE
        END IF
      END DO
      END DO
!$OMP END PARALLEL DO
C**** Define polar values at all latitudes
      DO L=1,LMM(1,JM)
C        PO(2:IM, 1,L) =  PO(1,IM,L)
C       PHI(2:IM, 1,L) = PHI(1,IM,L)
C        DH(2:IM, 1,L) =  DH(1,IM,L)
         PO(2:IM,JM,L) =  PO(1,JM,L)
        PHI(2:IM,JM,L) = PHI(1,JM,L)
         DH(2:IM,JM,L) =  DH(1,JM,L)
      END DO
C      OGEOZ(2:IM, 1) = OGEOZ(1,1)
      OGEOZ(2:IM,JM) = OGEOZ(1,JM)
C****
C**** Calculate smoothed East-West Pressure Gradient Force
C****
!$OMP PARALLEL DO  PRIVATE(I,IP1,J,L)
      DO 230 J=2,JM-1
      I=IM
      DO 230 IP1=1,IM
      DO 210 L=1,LMU(I,J)
  210 DUM(I,J,L) = -DT2*DYPO(J)*
     *  (( DH(IP1,J,L)+ DH(I,J,L))*( PO(IP1,J,L)- PO(I,J,L)) +
     +   (PHI(IP1,J,L)-PHI(I,J,L))*(MO(IP1,J,L)+MO(I,J,L)))
      DO 220 L=LMU(I,J)+1,LMO
  220 DUM(I,J,L) = 0.
  230 I=IP1
!$OMP END PARALLEL DO
      CALL OPFIL (DUM,LMO)
C****
C**** Calculate North-South Pressure Gradient Force
C****
!$OMP PARALLEL DO  PRIVATE(I,J,L)
      DO 310 J=1,JM-1
      DO 310 I=1,IM
      DO 310 L=1,LMV(I,J)
  310 DVM(I,J,L) = -DT2*DXVO(J)*
     *  (( DH(I,J+1,L)+ DH(I,J,L))*(PO(I,J+1,L)-PO(I,J,L)) +
     +   (PHI(I,J+1,L)-PHI(I,J,L))*(MO(I,J+1,L)+MO(I,J,L)))
!$OMP END PARALLEL DO
C****
C**** Add pressure gradient force to momentum
C****
!$OMP PARALLEL DO  PRIVATE(I,L)
      DO 610 I=IM+1,IM*(JM-1)
      DO 610 L=1,LMU(I,1)
  610 UM(I,1,L) = UM(I,1,L) + DUM(I,1,L)
!$OMP END PARALLEL DO
!$OMP PARALLEL DO  PRIVATE(I,L)
      DO 630 I=1,IM*(JM-1)
      DO 630 L=1,LMV(I,1)
  630 VM(I,1,L) = VM(I,1,L) + DVM(I,1,L)
!$OMP END PARALLEL DO
      RETURN
C****
      END SUBROUTINE OPGF

      SUBROUTINE OPGF0
!@sum  OPGF0 calculates specific volume of each grid box
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
C****
C**** OPGF0 calculates VOLGSP at two locations inside each grid box,
C**** which will be kept constant during the dynamics of momentum,
C**** and calculates two weighted averages from those values
C****
      USE CONSTANT, only : grav,rrt12=>byrt12
      USE OCEAN, only : im,jm,lmo,g0m,gzmo,s0m,szmo,mo,opress,lmm
     *     ,bydxypo,imaxj,opbot
      USE OCEAN_DYN, only : dzgdp,vbar,mmi
      IMPLICIT NONE
      INTEGER I,J,L
      REAL*8, SAVE :: PE
      REAL*8 PUP,PDN,GUP,GDN,SUP,SDN,VUP,VDN,PBAR,BYMMI
      REAL*8 VOLGSP

!$OMP PARALLEL DO  PRIVATE(I,J,L,PE,BYMMI,PBAR,PUP,PDN,GUP,GDN,
!$OMP&                     SUP,SDN,VUP,VDN)
      DO J=1,JM
      DO I=1,IMAXJ(J)
        PE = OPRESS(I,J)
        DO L=1,LMM(I,J)
          BYMMI=1d0/MMI(I,J,L)
          PBAR= PE   + MO(I,J,L)*(GRAV*5d-1)
          PE  = PE   + MO(I,J,L)* GRAV
          PUP = PBAR - MO(I,J,L)*(GRAV*RRT12)
          PDN = PBAR + MO(I,J,L)*(GRAV*RRT12)
          GUP = (G0M(I,J,L)-2d0*RRT12*GZMO(I,J,L))*BYMMI
          GDN = (G0M(I,J,L)+2d0*RRT12*GZMO(I,J,L))*BYMMI
          SUP = (S0M(I,J,L)-2d0*RRT12*SZMO(I,J,L))*BYMMI
          SDN = (S0M(I,J,L)+2d0*RRT12*SZMO(I,J,L))*BYMMI
          VUP = VOLGSP(GUP,SUP,PUP)
          VDN = VOLGSP(GDN,SDN,PDN)
          DZGDP(I,J,L) = (VUP*(5d-1-RRT12) + VDN*(5d-1+RRT12))*5d-1
          VBAR(I,J,L) = (VUP + VDN)*5d-1
        END DO
        OPBOT(I,J) = PE
      END DO
      END DO
!$OMP END PARALLEL DO
      RETURN
      END SUBROUTINE OPGF0

      SUBROUTINE OADVT (RM,RX,RY,RZ,DT,QLIMIT,OIJL)
!@sum  OADVT advects tracers using the linear upstream scheme.
!@auth Gary Russell
!@ver  1.0
C****
C**** Input:  MB (kg) = mass before advection
C****          DT (s) = time step
C****       MU (kg/s) = west to east mass flux
C****       MV (kg/s) = south to north mass flux
C****       MW (kg/s) = downward vertical mass flux
C****          QLIMIT = whether slope limitations should be used
C**** Output: RM (kg) = tracer mass
C****   RX,RY,RZ (kg) = first moments of tracer mass
C****       OIJL (kg) = diagnostic accumulation of tracer mass flux
C****
      USE OCEAN, only : im,jm,lmo
      USE OCEAN_DYN, only : mb=>mmi,smu,smv,smw
      IMPLICIT NONE
      REAL*8, INTENT(INOUT), DIMENSION(IM,JM,LMO) :: RM,RX,RY,RZ
      REAL*8, INTENT(INOUT), DIMENSION(IM,JM,LMO,3) :: OIJL
      REAL*8, DIMENSION(IM,JM,LMO) :: MA
      INTEGER I,J,L
      LOGICAL, INTENT(IN) :: QLIMIT
      REAL*8, INTENT(IN) :: DT
C****
C**** Load mass after advection from mass before advection
C****
      MA = MB
C****
C**** Advect the tracer using the Slopes Scheme
C****
      CALL OADVTX (RM,RX,RY,RZ,MA,SMU,5d-1*DT,QLIMIT,OIJL)
      CALL OADVTY (RM,RX,RY,RZ,MA,SMV,     DT,QLIMIT,OIJL(1,1,1,2))
      CALL OADVTZ (RM,RX,RY,RZ,MA,SMW,     DT,QLIMIT,OIJL(1,1,1,3))
      CALL OADVTX (RM,RX,RY,RZ,MA,SMU,5d-1*DT,QLIMIT,OIJL)
C**** Fill in values at the poles
      DO 20 L=1,LMO
      DO 20 I=1,IM
C     RM(I, 1,L) = RM(IM,1,L)
C     RZ(I, 1,L) = RZ(IM,1,L)
      RM(I,JM,L) = RM(1,JM,L)
   20 RZ(I,JM,L) = RZ(1,JM,L)
      RETURN
      END SUBROUTINE OADVT

      SUBROUTINE OADVTX (RM,RX,RY,RZ,MO,MU,DT,QLIMIT,OIJL)
!@sum  OADVTX advects tracer in x direction using linear upstream scheme
!@auth Gary Russell
!@ver  1.0
C**** If QLIMIT is true, the gradients are
C**** limited to prevent the mean tracer from becoming negative.
C****
C**** Input: DT (s) = time step
C****     MU (kg/s) = west to east mass flux
C****        QLIMIT = whether slope limitations should be used
C**** Input and Output: RM (kg) = tracer mass
C****             RX,RY,RZ (kg) = first moments of tracer mass
C****                    M (kg) = ocean mass
C****
      USE OCEAN, only : im,jm,lmo,lmm
      IMPLICIT NONE
      REAL*8, INTENT(INOUT), DIMENSION(IM,JM,LMO) :: RM,RX,RY,RZ
      REAL*8, INTENT(INOUT), DIMENSION(IM,JM,LMO) :: OIJL
      REAL*8, INTENT(INOUT), DIMENSION(IM,JM,LMO) :: MO
      REAL*8, INTENT(IN), DIMENSION(IM,JM,LMO) :: MU
      LOGICAL*4, INTENT(IN) :: QLIMIT
      REAL*8, INTENT(IN) :: DT
      REAL*8, DIMENSION(IM) :: AM,A,FM,FX,FY,FZ
      REAL*8 RXY
      INTEGER I,J,L,IM1,IP1,ICKERR

C**** Loop over layers and latitudes
      ICKERR=0
!$OMP PARALLEL DO  PRIVATE(I,J,L,IP1,IM1,A,AM,FM,FX,FY,FZ)
!$OMP&             REDUCTION(+:ICKERR)
      DO 320 L=1,LMO
      DO 320 J=2,JM-1
C****
C**** Calculate FM (kg), FX (kg**2), FY (kg) and FZ (kg)
C****
      I=IM
      DO 140 IP1=1,IM
      AM(I) = DT*MU(I,J,L)
      IF(AM(I)) 110,120,130
C**** Ocean mass flux is negative
  110 A(I)  = AM(I)/MO(IP1,J,L)
      IF(A(I).LT.-1d0)  WRITE (6,*) 'A<-1:',I,J,L,A(I),MO(IP1,J,L)
      FM(I) = A(I)*(RM(IP1,J,L)-(1d0+A(I))*RX(IP1,J,L))
      FX(I) = AM(I)*(A(I)*A(I)*RX(IP1,J,L)-3d0*FM(I))
      FY(I) = A(I)*RY(IP1,J,L)
      FZ(I) = A(I)*RZ(IP1,J,L)
      GO TO 140
C**** Ocean mass flux is zero
  120 A(I)  = 0.
      FM(I) = 0.
      FX(I) = 0.
      FY(I) = 0.
      FZ(I) = 0.
      GO TO 140
C**** Ocean mass flux is positive
  130 A(I)  = AM(I)/MO(I,J,L)
      IF(A(I).GT.1d0)  WRITE (6,*) 'A>1:',I,J,L,A(I),MO(I,J,L)
      FM(I) = A(I)*(RM(I,J,L)+(1d0-A(I))*RX(I,J,L))
      FX(I) = AM(I)*(A(I)*A(I)*RX(I,J,L)-3d0*FM(I))
      FY(I) = A(I)*RY(I,J,L)
      FZ(I) = A(I)*RZ(I,J,L)
  140 I=IP1
C****
C**** Modify the tracer moments so that the tracer mass in each
C**** division is non-negative
C****
      IF(.NOT.QLIMIT)  GO TO 300
      IM1=IM
      DO 290 I=1,IM
      IF(A(IM1).GE.0.)  GO TO 240
C**** Water is leaving through the left edge: 2 or 3 divisions
      IF(FM(IM1).LE.0.)  GO TO 210
C**** Left most division is negative, RML = -FM(I-1) < 0: Case 2 or 4
      RX(I,J,L) = RM(I,J,L)/(1d0+A(IM1))
      FM(IM1) = 0.
      FX(IM1) = AM(IM1)*A(IM1)*A(IM1)*RX(I,J,L)
      IF(A(I).LE.0.)  GO TO 290
      FM(I) = A(I)*(RM(I,J,L)+(1d0-A(I))*RX(I,J,L))
      FX(I) = AM(I)*(A(I)*A(I)*RX(I,J,L)-3d0*FM(I))
      GO TO 290
C**** Left most division is non-negative, RML = -FM(I-() > 0:
C**** Case 1, 3 or 5
  210 IF(A(I).LE.0.)  GO TO 230
C**** Water is leaving through the right edge: 3 divisions
      IF(FM(I).GE.0.)  GO TO 290
C**** Right most division is negative, RMR = FM(I) < 0: Case 3 or 5
  220 RX(I,J,L) = -RM(I,J,L)/(1d0-A(I))
      FM(I) = 0.
      FX(I) = AM(I)*A(I)*A(I)*RX(I,J,L)
      FM(IM1) = A(IM1)*(RM(I,J,L)-(1d0+A(IM1))*RX(I,J,L))
      FX(IM1) = AM(IM1)*(A(IM1)*A(IM1)*RX(I,J,L)-3d0*FM(IM1))
      GO TO 290
C**** No water is leaving through the right edge: 2 divisions
  230 IF(RM(I,J,L)+FM(IM1).GE.0.)  GO TO 290
C**** Right most division is negative, RMR = RM(I,J)+FM(I-1) < 0: Case 3
      RX(I,J,L) = RM(I,J,L)/A(IM1)
      FM(IM1) = -RM(I,J,L)
      FX(IM1) = AM(IM1)*(A(IM1)+3d0)*RM(I,J,L)
      GO TO 290
C**** No water is leaving through the left edge: 1 or 2 divisions
  240 IF(A(I).LE.0.)  GO TO 290
C**** Water is leaving through the right edge: 2 divisions
      IF(FM(I).GE.0.)  GO TO 250
C**** Right most division is negative, RMR = FM(I) < 0: Case 3
      RX(I,J,L) = -RM(I,J,L)/(1d0-A(I))
      FM(I) = 0.
      FX(I) = AM(I)*A(I)*A(I)*RX(I,J,L)
      GO TO 290
C**** Right most division is non-negative, RMR = FM(I) > 0: Case 1 or 2
  250 IF(RM(I,J,L)-FM(I).GE.0.)  GO TO 290
C**** Left most division is negative, RML = RM(I,J)-FM(I) < 0: Case 2
      RX(I,J,L) = RM(I,J,L)/A(I)
      FM(I) = RM(I,J,L)
      FX(I) = AM(I)*(A(I)-3d0)*RM(I,J,L)
C****
  290 IM1=I
C****
C**** Calculate new tracer mass and first moments of tracer mass
C****
  300 IM1=IM
      DO 310 I=1,IM
      IF(L.GT.LMM(I,J))  GO TO 310
      RM(I,J,L) = RM(I,J,L) +  (FM(IM1)-FM(I))
      RX(I,J,L) = (RX(I,J,L)*MO(I,J,L) + (FX(IM1)-FX(I))
     *  + 3d0*((AM(IM1)+AM(I))*RM(I,J,L)-MO(I,J,L)*(FM(IM1)+FM(I))))
     *  / (MO(I,J,L)+AM(IM1)-AM(I))
      RY(I,J,L) = RY(I,J,L) + (FY(IM1)-FY(I))
      RZ(I,J,L) = RZ(I,J,L) + (FZ(IM1)-FZ(I))
C****
      if ( QLIMIT ) then ! limit tracer gradients
        RXY = abs(RX(I,J,L)) + abs(RY(I,J,L))
        if ( RXY > RM(I,J,L) ) then
          RX(I,J,L) = RX(I,J,L)*( RM(I,J,L)/(RXY + tiny(RXY)) )
          RY(I,J,L) = RY(I,J,L)*( RM(I,J,L)/(RXY + tiny(RXY)) )
        end if
        if ( abs(RZ(I,J,L)) > RM(I,J,L) )
     *       RZ(I,J,L) = sign(RM(I,J,L), RZ(I,J,L)+0d0)
      end if
C****
      MO(I,J,L) = MO(I,J,L) +  AM(IM1)-AM(I)
         IF(MO(I,J,L).LE.0.)                ICKERR=ICKERR+1
         IF(QLIMIT .AND. RM(I,J,L).LT.0.)   ICKERR=ICKERR+1
         OIJL(I,J,L) = OIJL(I,J,L) + FM(I)
  310 IM1=I
  320 CONTINUE
!$OMP END PARALLEL DO

C**** IF NO ERROR HAS OCCURRED - RETURN, ELSE STOP
      IF(ICKERR.EQ.0)  RETURN
      DO 440 L=1,LMO
      DO 440 J=2,JM-1
      DO 440 I=1,IM
         IF(MO(I,J,L).LE.0.)  GO TO 800
         IF(QLIMIT .AND. RM(I,J,L).LT.0.) GO TO 810
  440 CONTINUE
      WRITE(6,*) 'ERROR CHECK INCONSISTENCY: OADVTX ',ICKERR
      call stop_model("OADVTX",255)

  800 WRITE (6,*) 'MO<0 in OADVTX:',I,J,L,MO(I,J,L)
  810 WRITE (6,*) 'RM in OADVTX:',I,J,L,RM(I,J,L)
      WRITE (6,*) 'A=',(I,A(I),I=1,IM)
      call stop_model("OADVTX",255)
      END SUBROUTINE OADVTX

      SUBROUTINE OADVTY (RM,RX,RY,RZ,MO,MV,DT,QLIMIT,OIJL)
C****
C**** OADVTY advects tracers in the south to north direction using the
C**** linear upstream scheme.  If QLIMIT is true, the gradients are
C**** limited to prevent the mean tracer from becoming negative.
C****
C**** Input: DT (s) = time step
C****     MV (kg/s) = south to north mass flux
C****        QLIMIT = whether slope limitations should be used
C**** Input and Output: RM (kg) = tracer mass
C****             RX,RY,RZ (kg) = first moments of tracer mass
C****                    M (kg) = ocean mass
C****
      USE OCEAN, only : im,jm,lmo,lmm
      IMPLICIT NONE
      REAL*8, INTENT(INOUT), DIMENSION(IM,JM,LMO) :: RM,RX,RY,RZ
      REAL*8, INTENT(INOUT), DIMENSION(IM,JM,LMO) :: OIJL
      REAL*8, INTENT(INOUT), DIMENSION(IM,JM,LMO) :: MO
      REAL*8, INTENT(IN), DIMENSION(IM,JM,LMO) :: MV
      LOGICAL*4, INTENT(IN) :: QLIMIT
      REAL*8, INTENT(IN) :: DT
      REAL*8, DIMENSION(JM) :: BM,B,FM,FX,FY,FZ
      INTEGER I,J,L,ICKERR
      REAL*8 SBMN,SFMN,SFZN,RXY

C**** Loop over layers and longitudes
      ICKERR=0
!$OMP PARALLEL DO  PRIVATE(I,J,L, SBMN,SFMN,SFZN, B,BM,FM,FX,FY,FZ),
!$OMP&             REDUCTION(+:ICKERR)
      DO 350 L=1,LMO
C     SBMS = 0.
C     SFMS = 0.
C     SFZS = 0.
      SBMN = 0.
      SFMN = 0.
      SFZN = 0.
      DO 330 I=1,IM
C****
C**** Calculate FM (kg), FX (kg), FY (kg**2) and FZ (kg) near South Pole
C****
      BM(1) = DT*MV(I,1,L)
C     IF(BM(1)) 101,102,103
C**** Ocean mass flux is negative
C 101 B(1)  = BM(1)/MO(I,2,L)
C     IF(B(1).LT.-1d0)  WRITE (6,*) 'B<-1:',I,1,L,B(1),MO(1,2,L)
C     FM(1) = B(1)*(RM(I,2,L)-(1d0+B(1))*RY(I,2,L))
C     FX(1) = B(1)*RX(I,2,L)
C     FY(1) = BM(1)*(B(1)*B(1)*RY(I,2,L)-3d0*FM(1))
C     FZ(1) = B(1)*RZ(I,2,L)
C     GO TO 110
C**** Ocean mass flux is zero
  102 B(1)  = 0.
      FM(1) = 0.
      FX(1) = 0.
      FY(1) = 0.
      FZ(1) = 0.
C     GO TO 110
C**** Ocean mass flux is positive
C 103 B(1)  = BM(1)/MO(IM,1,L)
C     IF(B(1).GT.1d0)  WRITE (6,*) 'B>1:',I,1,L,B(1),MO(IM,1,L)
C     FM(1) = B(1)*RM(IM,1,L)
C     FX(1) = 0.
C     FY(1) = -3d0*BM(1)*FM(1)
C     FZ(1) = B(1)*RZ(IM,1,L)
C****
C**** Calculate FM (kg), FX (kg), FY (kg**2) and FZ (kg) in the interior
C****
  110 DO 114 J=2,JM-2
      BM(J) = DT*MV(I,J,L)
      IF(BM(J)) 111,112,113
C**** Ocean mass flux is negative
  111 B(J)  = BM(J)/MO(I,J+1,L)
      IF(B(J).LT.-1d0)  WRITE (6,*) 'B<-1:',I,J,L,B(J),MO(I,J+1,L)
      FM(J) = B(J)*(RM(I,J+1,L)-(1d0+B(J))*RY(I,J+1,L))
      FX(J) = B(J)*RX(I,J+1,L)
      FY(J) = BM(J)*(B(J)*B(J)*RY(I,J+1,L)-3d0*FM(J))
      FZ(J) = B(J)*RZ(I,J+1,L)
      GO TO 114
C**** Ocean mass flux is zero
  112 B(J)  = 0.
      FM(J) = 0.
      FX(J) = 0.
      FY(J) = 0.
      FZ(J) = 0.
      GO TO 114
C**** Ocean mass flux is positive
  113 B(J)  = BM(J)/MO(I,J,L)
      IF(B(J).GT.1d0)  WRITE (6,*) 'B>1:',I,J,L,B(J),MO(I,J,L)
      FM(J) = B(J)*(RM(I,J,L)+(1d0-B(J))*RY(I,J,L))
      FX(J) = B(J)*RX(I,J,L)
      FY(J) = BM(J)*(B(J)*B(J)*RY(I,J,L)-3d0*FM(J))
      FZ(J) = B(J)*RZ(I,J,L)
  114 CONTINUE
C****
C**** Calculate FM (kg), FX (kg), FY (kg**2) and FZ (kg) near North Pole
C****
      BM(JM-1) = DT*MV(I,JM-1,L)
      IF(BM(JM-1)) 121,122,123
C**** Ocean mass flux is negative
  121 B(JM-1)  = BM(JM-1)/MO(1,JM,L)
      IF(B(JM-1).LT.-1d0) WRITE(6,*) 'B<-1:',I,JM-1,L,B(JM-1),MO(1,JM,L)
      FM(JM-1) = B(JM-1)*RM(1,JM,L)
      FX(JM-1) = 0.
      FY(JM-1) = -3d0*BM(JM-1)*FM(JM-1)
      FZ(JM-1) = B(JM-1)*RZ(1,JM,L)
      GO TO 200
C**** Ocean mass flux is zero
  122 B(JM-1)  = 0.
      FM(JM-1) = 0.
      FX(JM-1) = 0.
      FY(JM-1) = 0.
      FZ(JM-1) = 0.
      GO TO 200
C**** Ocean mass flux is positive
  123 B(JM-1)  = BM(JM-1)/MO(I,JM-1,L)
      IF(B(JM-1).GT.1d0) WRITE(6,*) 'B>1:',I,JM-1,L,B(JM-1),MO(1,JM-1,L)
      FM(JM-1) = B(JM-1)*(RM(I,JM-1,L)+(1d0-B(JM-1))*RY(I,JM-1,L))
      FX(JM-1) = B(JM-1)*RX(I,JM-1,L)
      FY(JM-1) = BM(JM-1)*(B(JM-1)*B(JM-1)*RY(I,JM-1,L)-3d0*FM(JM-1))
      FZ(JM-1) = B(JM-1)*RZ(I,JM-1,L)
C****
C**** Modify the tracer moments so that the tracer mass in each
C**** division is non-negative
C****
  200 IF(.NOT.QLIMIT)  GO TO 300
      DO 290 J=2,JM-1
      IF(B(J-1).GE.0.)  GO TO 240
C**** Water is leaving through the south edge: 2 or 3 divisions
      IF(FM(J-1).LE.0.)  GO TO 210
C**** South most division is negative, RMS = -FM(J-1) < 0: Case 2 or 4
      RY(I,J,L) = RM(I,J,L)/(1d0+B(J-1))
      FM(J-1) = 0.
      FY(J-1) = BM(J-1)*B(J-1)*B(J-1)*RY(I,J,L)
      IF(B(J).LE.0.)  GO TO 290
      FM(J) = B(J)*(RM(I,J,L)+(1d0-B(J))*RY(I,J,L))
      FY(J) = BM(J)*(B(J)*B(J)*RY(I,J,L)-3d0*FM(J))
      GO TO 290
C**** South most division is non-negative, RMS = -FM(J-1) > 0:
C**** Case 1, 3 or 5
  210 IF(B(J).LE.0.)  GO TO 230
C**** Water is leaving through the north edge: 3 divisions
      IF(FM(J).GE.0.)  GO TO 290
C**** North most division is negative, RMN = FM(J) < 0: Case 3 or 5
  220 RY(I,J,L) = -RM(I,J,L)/(1d0-B(J))
      FM(J) = 0.
      FY(J) = BM(J)*B(J)*B(J)*RY(I,J,L)
      FM(J-1) = B(J-1)*(RM(I,J,L)-(1d0+B(J-1))*RY(I,J,L))
      FY(J-1) = BM(J-1)*(B(J-1)*B(J-1)*RY(I,J,L)-3d0*FM(J-1))
      GO TO 290
C**** No water is leaving through the north edge: 2 divisions
  230 IF(RM(I,J,L)+FM(J-1).GE.0.)  GO TO 290
C**** North most division is negative, RMN = RM(I,J)+FM(J-1) < 0: Case 3
      RY(I,J,L) = RM(I,J,L)/B(J-1)
      FM(J-1) = -RM(I,J,L)
      FY(J-1) = BM(J-1)*(B(J-1)+3d0)*RM(I,J,L)
      GO TO 290
C**** No water is leaving through the south edge: 1 or 2 divisions
  240 IF(B(J).LE.0.)  GO TO 290
C**** Water is leaving through the north edge: 2 divisions
      IF(FM(J).GE.0.)  GO TO 250
C**** North most division is negative, RMN = FM(J) < 0: Case 3
      RY(I,J,L) = -RM(I,J,L)/(1d0-B(J))
      FM(J) = 0.
      FY(J) = BM(J)*B(J)*B(J)*RY(I,J,L)
      GO TO 290
C**** North most division is non-negative, RMN = FM(J) > 0: Case 1 or 2
  250 IF(RM(I,J,L)-FM(J).GE.0.)  GO TO 290
C**** South most division is negative, RMS = RM(I,J)-FM(J) < 0: Case 2
      RY(I,J,L) = RM(I,J,L)/B(J)
      FM(J) = RM(I,J,L)
      FY(J) = BM(J)*(B(J)-3d0)*RM(I,J,L)
C****
  290 CONTINUE
C 300 SBMS = SBMS + BM(1)
C     SFMS = SFMS + FM(1)
C     SFZS = SFZS + FZ(1)
  300 SBMN = SBMN + BM(JM-1)
      SFMN = SFMN + FM(JM-1)
      SFZN = SFZN + FZ(JM-1)
C****
C**** Calculate new tracer mass and first moments of tracer mass
C****
      DO 310 J=2,JM-1
      IF(L.GT.LMM(I,J))  GO TO 310
      RM(I,J,L) = RM(I,J,L) + (FM(J-1)-FM(J))
      RX(I,J,L) = RX(I,J,L) + (FX(J-1)-FX(J))
      RY(I,J,L) = (RY(I,J,L)*MO(I,J,L) + (FY(J-1)-FY(J))
     *  + 3d0*((BM(J-1)+BM(J))*RM(I,J,L)-MO(I,J,L)*(FM(J-1)+FM(J))))
     *  / (MO(I,J,L)+BM(J-1)-BM(J))
      RZ(I,J,L) = RZ(I,J,L) + (FZ(J-1)-FZ(J))
C****
      if ( QLIMIT ) then ! limit tracer gradients
        RXY = abs(RX(I,J,L)) + abs(RY(I,J,L))
        if ( RXY > RM(I,J,L) ) then
          RX(I,J,L) = RX(I,J,L)*( RM(I,J,L)/(RXY + tiny(RXY)) )
          RY(I,J,L) = RY(I,J,L)*( RM(I,J,L)/(RXY + tiny(RXY)) )
        end if
        if ( abs(RZ(I,J,L)) > RM(I,J,L) )
     *       RZ(I,J,L) = sign(RM(I,J,L), RZ(I,J,L)+0d0)
      end if
C****
      MO(I,J,L) = MO(I,J,L) +  BM(J-1)-BM(J)
         IF(MO(I,J,L).LE.0.)             ICKERR=ICKERR+1
         IF(QLIMIT.AND.RM(I,J,L).LT.0.)  ICKERR=ICKERR+1
  310 CONTINUE
         DO 320 J=1,JM-1
  320    OIJL(I,J,L) = OIJL(I,J,L) + FM(J)
  330    CONTINUE
C     IF(L.GT.LMM(IM,1))  GO TO 340
C     RM(IM,1,L) = RM(IM,1,L) - SFMS/IM
C     RZ(IM,1,L) = RZ(IM,1,L) - SFZS/IM
C     MO(IM,1,L) =  MO(IM,1,L) - SBMS/IM
C        IF(MO(IM,1,L).LE.0.)  GO TO 800
C        IF(QLIMIT.AND.RM(IM,1,L).LT.0.)  GO TO 810
  340 IF(L.GT.LMM(1,JM))  GO TO 350
      RM(1,JM,L) = RM(1,JM,L) + SFMN/IM
      RZ(1,JM,L) = RZ(1,JM,L) + SFZN/IM
      MO(1,JM,L) = MO(1,JM,L) + SBMN/IM
         IF(MO(1,JM,L).LE.0.)             ICKERR=ICKERR+1
         IF(QLIMIT.AND.RM(1,JM,L).LT.0.)  ICKERR=ICKERR+1
  350 CONTINUE
!$OMP END PARALLEL DO

C**** IF NO ERROR HAS OCCURRED - RETURN, ELSE STOP
      IF(ICKERR.EQ.0)  RETURN
      DO L=1,LMO
        DO J=2,JM-1
          DO I=1,IM
            IF(MO(I,J,L).LE.0.)  GO TO 800
            IF(QLIMIT .AND. RM(I,J,L).LT.0.) GO TO 810
          END DO
        END DO
        IF(MO(1,JM,L).LE.0.)  GO TO 800
        IF(QLIMIT.AND.RM(1,JM,L).LT.0.)  GO TO 810
      END DO
      WRITE(6,*) 'ERROR CHECK INCONSISTENCY: OADVTY ',ICKERR
      call stop_model("OADVTY",255)

  800 WRITE (6,*) 'MO<0 in OADVTY:',I,J,L,MO(I,J,L)
  810 WRITE (6,*) 'RM in OADVTY:',I,J,L,RM(I,J,L)
      WRITE (6,*) 'B=',(J,B(J),J=1,JM-1)
      call stop_model("OADVTY",255)
      END SUBROUTINE OADVTY

      SUBROUTINE OADVTZ (RM,RX,RY,RZ,MO,MW,DT,QLIMIT,OIJL)
C****
C**** OADVTZ advects tracers in the vertical direction using the
C**** linear upstream scheme.  If QLIMIT is true, the gradients are
C**** limited to prevent the mean tracer from becoming negative.
C****
C**** Input: DT (s) = time step
C****     MW (kg/s) = downward vertical mass flux
C****        QLIMIT = whether slope limitations should be used
C**** Input and Output: RM (kg) = tracer mass
C****             RX,RY,RZ (kg) = first moments of tracer mass
C****                    M (kg) = ocean mass
C****
      USE OCEAN, only : im,jm,lmo,lmm
      IMPLICIT NONE
      REAL*8, INTENT(INOUT), DIMENSION(IM,JM,LMO) :: RM,RX,RY,RZ
      REAL*8, INTENT(INOUT), DIMENSION(IM,JM,LMO) :: OIJL
      REAL*8, INTENT(INOUT), DIMENSION(IM,JM,LMO) :: MO
      REAL*8, INTENT(IN), DIMENSION(IM,JM,LMO-1) :: MW
      LOGICAL*4, INTENT(IN) :: QLIMIT
      REAL*8, INTENT(IN) :: DT
      REAL*8, DIMENSION(0:LMO) :: CM,C,FM,FX,FY,FZ
      INTEGER I,J,L,LMIJ,ICKERR
      REAL*8 SBMN,SFMN,SFZN,RXY
C****
C**** Loop over latitudes and longitudes
      ICKERR=0
!$OMP PARALLEL DO  PRIVATE(I,L,LMIJ, C,CM, FM,FX,FY,FZ)
!$OMP&             REDUCTION(+:ICKERR)
      DO 330 I=IM,IM*(JM-1)+1
      CM(0) = 0.
       C(0) = 0.
      FM(0) = 0.
      FX(0) = 0.
      FY(0) = 0.
      FZ(0) = 0.
      LMIJ=LMM(I,1)
      IF(LMIJ.LE.1)  GO TO 330
      CM(LMIJ) = 0.
       C(LMIJ) = 0.
      FM(LMIJ) = 0.
      FX(LMIJ) = 0.
      FY(LMIJ) = 0.
      FZ(LMIJ) = 0.
C****
C**** Calculate FM (kg), FX (kg), FY (kg) and FZ (kg**2)
C****
      DO 120 L=1,LMIJ-1
      CM(L) = DT*MW(I,1,L)
      IF(CM(L).LT.0.)  GO TO 110
C**** Ocean mass flux is positive
      C(L)  = CM(L)/MO(I,1,L)
      IF(C(L).GT.1d0)  WRITE (6,*) 'C>1:',I,L,C(L),MO(I,1,L)
      FM(L) = C(L)*(RM(I,1,L)+(1d0-C(L))*RZ(I,1,L))
      FX(L) = C(L)*RX(I,1,L)
      FY(L) = C(L)*RY(I,1,L)
      FZ(L) = CM(L)*(C(L)*C(L)*RZ(I,1,L)-3d0*FM(L))
      GO TO 120
C**** Ocean mass flux is negative
  110 C(L)  = CM(L)/MO(I,1,L+1)
      IF(C(L).LT.-1d0)  WRITE (6,*) 'C<-1:',I,L,C(L),MO(I,1,L+1)
      FM(L) = C(L)*(RM(I,1,L+1)-(1d0+C(L))*RZ(I,1,L+1))
      FX(L) = C(L)*RX(I,1,L+1)
      FY(L) = C(L)*RY(I,1,L+1)
      FZ(L) = CM(L)*(C(L)*C(L)*RZ(I,1,L+1)-3d0*FM(L))
  120 CONTINUE
C****
C**** Modify the tracer moments so that the tracer mass in each
C**** division is non-negative
C****
      IF(.NOT.QLIMIT)  GO TO 300
      DO 290 L=1,LMIJ
      IF(C(L-1).GE.0.)  GO TO 240
C**** Water is leaving through the bottom edge: 2 or 3 divisions
      IF(FM(L-1).LE.0.)  GO TO 210
C**** Bottom most division is negative, RMB = -FM(L-1) < 0: Case 2 or 4
      RZ(I,1,L) = RM(I,1,L)/(1d0+C(L-1))
      FM(L-1) = 0.
      FZ(L-1) = CM(L-1)*C(L-1)*C(L-1)*RZ(I,1,L)
      IF(C(L).LE.0.)  GO TO 290
      FM(L) = C(L)*(RM(I,1,L)+(1d0-C(L))*RZ(I,1,L))
      FZ(L) = CM(L)*(C(L)*C(L)*RZ(I,1,L)-3d0*FM(L))
      GO TO 290
C**** Bottom most division is non-negative, RMB = -FM(L-1) > 0:
C**** Case 1, 3 or 5
  210 IF(C(L).LE.0.)  GO TO 230
C**** Water is leaving through the top edge: 3 divisions
      IF(FM(L).GE.0.)  GO TO 290
C**** Top most division is negative, RMT = FM(L) < 0: Case 3 or 5
  220 RZ(I,1,L) = -RM(I,1,L)/(1d0-C(L))
      FM(L) = 0.
      FZ(L) = CM(L)*C(L)*C(L)*RZ(I,1,L)
      FM(L-1) = C(L-1)*(RM(I,1,L)-(1d0+C(L-1))*RZ(I,1,L))
      FZ(L-1) = CM(L-1)*(C(L-1)*C(L-1)*RZ(I,1,L)-3d0*FM(L-1))
      GO TO 290
C**** No water is leaving through the top edge: 2 divisions
  230 IF(RM(I,1,L)+FM(L-1).GE.0.)  GO TO 290
C**** Top most division is negative, RMT = RM(I,1,L)+FM(L-1) < 0: Case 3
      RZ(I,1,L) = RM(I,1,L)/C(L-1)
      FM(L-1) = -RM(I,1,L)
      FZ(L-1) = CM(L-1)*(C(L-1)+3d0)*RM(I,1,L)
      GO TO 290
C**** No water is leaving through the bottom edge: 1 or 2 divisions
  240 IF(C(L).LE.0.)  GO TO 290
C**** Water is leaving through the top edge: 2 divisions
      IF(FM(L).GE.0.)  GO TO 250
C**** Top most division is negative, RMT = FM(L) < 0: Case 3
      RZ(I,1,L) = -RM(I,1,L)/(1d0-C(L))
      FM(L) = 0.
      FZ(L) = CM(L)*C(L)*C(L)*RZ(I,1,L)
      GO TO 290
C**** Top most division is non-negative, RMT = FM(L) > 0: Case 1 or 2
  250 IF(RM(I,1,L)-FM(L).GE.0.)  GO TO 290
C**** Bottom most division is negative, RMB = RM(I,1,L)-FM(L) < 0: Cas 2
      RZ(I,1,L) = RM(I,1,L)/C(L)
      FM(L) = RM(I,1,L)
      FZ(L) = CM(L)*(C(L)-3d0)*RM(I,1,L)
C****
  290 CONTINUE
C****
C**** Calculate new tracer mass and first moments of tracer mass
C****
  300 DO 310 L=1,LMIJ
      RM(I,1,L) = RM(I,1,L) + (FM(L-1)-FM(L))
      RX(I,1,L) = RX(I,1,L) + (FX(L-1)-FX(L))
      RY(I,1,L) = RY(I,1,L) + (FY(L-1)-FY(L))
      RZ(I,1,L) = (RZ(I,1,L)*MO(I,1,L) + (FZ(L-1)-FZ(L))
     *  + 3d0*((CM(L-1)+CM(L))*RM(I,1,L)-MO(I,1,L)*(FM(L-1)+FM(L))))
     *  / (MO(I,1,L)+CM(L-1)-CM(L))
C****
      if ( QLIMIT ) then ! limit tracer gradients
        RXY = abs(RX(I,1,L)) + abs(RY(I,1,L))
        if ( RXY > RM(I,1,L) ) then
          RX(I,1,L) = RX(I,1,L)*( RM(I,1,L)/(RXY + tiny(RXY)) )
          RY(I,1,L) = RY(I,1,L)*( RM(I,1,L)/(RXY + tiny(RXY)) )
        end if
        if ( abs(RZ(I,1,L)) > RM(I,1,L) )
     *       RZ(I,1,L) = sign(RM(I,1,L), RZ(I,1,L)+0d0)
      end if
C****
      MO(I,1,L) = MO(I,1,L) +  CM(L-1)-CM(L)
         IF(MO(I,1,L).LE.0.)              ICKERR=ICKERR+1
         IF(QLIMIT.AND.RM(I,1,L).LT.0.)   ICKERR=ICKERR+1
  310 CONTINUE
         DO 320 L=1,LMIJ-1
  320    OIJL(I,1,L) = OIJL(I,1,L) + FM(L)
  330 CONTINUE
!$OMP END PARALLEL DO

C**** IF NO ERROR HAS OCCURRED - RETURN, ELSE STOP
      IF(ICKERR.EQ.0)  RETURN
      DO I=IM,IM*(JM-1)+1
        LMIJ=LMM(I,1)
        DO L=1,LMIJ
          IF(MO(I,1,L).LE.0.)  GO TO 800
          IF(QLIMIT .AND. RM(I,1,L).LT.0.) GO TO 810
        END DO
      END DO
      WRITE(6,*) 'ERROR CHECK INCONSISTENCY: OADVTZ ',ICKERR
      call stop_model("OAVDTZ",255)

  800 WRITE (6,*) 'MO<0 in OADVTZ:',I,L,MO(I,1,L)
  810 WRITE (6,*) 'RM in OADVTZ:',I,L,RM(I,1,L)
      WRITE (6,*) 'C=',(L,C(L),L=0,LMIJ)
      call stop_model("OADVTZ",255)
      END SUBROUTINE OADVTZ

      SUBROUTINE OBDRAG
!@sum  OBDRAG exerts a drag on the Ocean Model's bottom layer
!@auth Gary Russell
!@ver  1.0
      USE OCEAN, only : im,jm,lmo,mo,uo,vo,lmu,lmv,dts
      IMPLICIT NONE
      INTEGER, PARAMETER :: IQ1=1+IM/4, IQ2=IQ1+IM/4, IQ3=IQ2+IM/4
      REAL*8, PARAMETER :: BDRAGX=1d0, SDRAGX=1d-1
      REAL*8, DIMENSION(IM,JM,LMO) :: UT,VT
      INTEGER I,J,L,IP1,JV,IM1
      REAL*8 WSQ
C****
C**** UO = UO*(1-x/y)  is approximated by  UO*y/(y+x)  for stability
C****
C**** Save UO,VO into UT,VT which will be unchanged
      DO 10 I=1,IM*JM*LMO
      UT(I,1,1) = UO(I,1,1)
   10 VT(I,1,1) = VO(I,1,1)
C****
C**** Reduce West-East ocean current
C****
C**** Bottom drag in the interior
      I=IM
      DO 120 J=2,JM-1
      DO 110 IP1=1,IM
      IF(LMU(I,J).LE.0)  GO TO 110
      L=LMU(I,J)
      WSQ = UT(I,J,L)*UT(I,J,L) + 1d-20 +
     *  2.5d-1*(VT(I,J  ,L)*VT(I,J  ,L) + VT(IP1,J  ,L)*VT(IP1,J  ,L)
     *     + VT(I,J-1,L)*VT(I,J-1,L) + VT(IP1,J-1,L)*VT(IP1,J-1,L))
      UO(I,J,L) = UO(I,J,L) * (MO(I,J,L)+MO(IP1,J,L)) /
     *           (MO(I,J,L)+MO(IP1,J,L) + DTS*BDRAGX*SQRT(WSQ)*2d0)
  110 I=IP1
  120 CONTINUE
C**** Bottom drag at the poles
      J=JM
      JV=JM-1
C     IF(LMU(1,J).LE.0)  GO TO
      L=LMU(1,J)
      WSQ = UT(1,J,L)*UT(1,J,L) + 1d-20 +
     *  2.5d-1*(VT(  1,JV,L)*VT(  1,JV,L) + VT(IQ1,JV,L)*VT(IQ1,JV,L)
     *     + VT(IQ2,JV,L)*VT(IQ2,JV,L) + VT(IQ3,JV,L)*VT(IQ3,JV,L))
      UO(1,J,L) = UO(1,J,L) * MO(1,J,L) /
     *           (MO(1,J,L) + DTS*BDRAGX*SQRT(WSQ))
      DO 130 I=2,IM
  130 UO(I,J,L) = UO(1,J,L)
C****
C**** Reduce South-North ocean current
C****
      IM1=IM
      DO 220 J=1,JM-1
      DO 210 I=1,IM
      IF(LMV(I,J).LE.0)  GO TO 210
      L=LMV(I,J)
      WSQ = VT(I,J,L)*VT(I,J,L) + 1d-20 +
     *  2.5d-1*(UT(IM1,J+1,L)*UT(IM1,J+1,L) + UT(I,J+1,L)*UT(I,J+1,L)
     *     + UT(IM1,J  ,L)*UT(IM1,J  ,L) + UT(I,J  ,L)*UT(I,J  ,L))
      VO(I,J,L) = VO(I,J,L) * (MO(I,J,L)+MO(I,J+1,L)) /
     *           (MO(I,J,L)+MO(I,J+1,L) + DTS*BDRAGX*SQRT(WSQ)*2d0)
  210 IM1=I
  220 CONTINUE
      RETURN
C****
      END SUBROUTINE OBDRAG

      SUBROUTINE OCOAST
!@sum OCOAST reduces the horizontal perpendicular gradients of tracers
!@sum in coastline ocean grid boxes
!@auth Gary Russell
!@ver  1.0
      USE CONSTANT, only : sday
      USE OCEAN, only : im,jm,dts,lmm,gxmo,gymo,sxmo,symo
#ifdef TRACERS_OCEAN
     *     ,ntm,txmo,tymo
#endif
      IMPLICIT NONE
      INTEGER I,IM1,IP1,J,LMIN,L,N
      REAL*8 REDUCE

      REDUCE = 1d0 - DTS/(SDAY*2d1)
C**** Reduce West-East gradient of tracers
      IM1=IM-1
      I=IM
      DO 120 J=2,JM-1
      DO 120 IP1=1,IM
      LMIN = MIN(LMM(IM1,J),LMM(IP1,J)) + 1
      DO 110 L=LMIN,LMM(I,J)
        GXMO(I,J,L) = GXMO(I,J,L)*REDUCE
        SXMO(I,J,L) = SXMO(I,J,L)*REDUCE
#ifdef TRACERS_OCEAN
        DO N = 1,NTM
          TXMO(I,J,L,N) = TXMO(I,J,L,N) *REDUCE
        END DO
#endif
 110  CONTINUE
      IM1=I
  120 I=IP1
C**** Reduce South-North gradient of tracers
      DO 220 J=2,JM-1
      DO 220 I=1,IM
      LMIN = MIN(LMM(I,J-1),LMM(I,J+1)) + 1
      DO 210 L=LMIN,LMM(I,J)
        GYMO(I,J,L) = GYMO(I,J,L)*REDUCE
        SYMO(I,J,L) = SYMO(I,J,L)*REDUCE
#ifdef TRACERS_OCEAN
        DO N = 1,NTM
          TYMO(I,J,L,N) = TYMO(I,J,L,N) *REDUCE
        END DO
#endif
 210  CONTINUE
  220 CONTINUE
      RETURN
      END SUBROUTINE OCOAST

      SUBROUTINE OSTRES
!@sum OSTRES applies the atmospheric surface stress over open ocean
!@sum and the sea ice stress to the layer 1 ocean velocities
!@auth Gary Russell
!@ver  1.0
      USE OCEAN, only : im,jm,uo,vo,mo,dxyso,dxyno,dxyvo,lmu,lmv,cosic
     *     ,sinic,ratoc
      USE FLUXES, only : dmua,dmva,dmui,dmvi
      IMPLICIT NONE
      INTEGER I,J,IP1
C****
C**** All stress now defined for whole box, not just ocn or ice fraction
C**** FLUXCB  DMUA(1)  U momentum downward into open ocean (kg/m*s)
C****         DMVA(1)  V momentum downward into open ocean (kg/m*s)
C****         DMUA(2,JM,1)  polar atmo. mass slowed to zero (kg/m**2)
C****         DMUI     U momentum downward from sea ice (kg/m*s)
C****         DMVI     V momentum downward from sea ice (kg/m*s)

C**** Scale stresses for ocean area
      DO J=1,JM
        DO I=1,IM
          DMUA(I,J,1)=RATOC(J)*DMUA(I,J,1)
          DMVA(I,J,1)=RATOC(J)*DMVA(I,J,1)
          DMUI(I,J)=RATOC(J)*DMUI(I,J)
          DMVI(I,J)=RATOC(J)*DMVI(I,J)
        END DO
      END DO
C****
C**** Surface stress is applied to U component
C****
      I=IM
      DO J=2,JM-1
      DO IP1=1,IM
        IF(LMU(I,J).gt.0.)  UO(I,J,1) = UO(I,J,1) +
     *       (DMUA(I,J,1) + DMUA(IP1,J,1) + 2d0*DMUI(I,J)) /
     *       (  MO(I,J,1) +   MO(IP1,J,1))
        I=IP1
      END DO
      END DO
      UO(1,JM,1) = UO(1,JM,1)*(1d0 - DMUA(2,JM,1)/MO(1,JM,1))
C****
C**** Surface stress is applied to V component
C****
      DO J=2,JM-2
      DO I=1,IM
        IF(LMV(I,J).GT.0.)  VO(I,J,1) = VO(I,J,1) +
     *       (DMVA(I,J  ,1)*DXYNO(J) + DMVA(I,J+1,1)*DXYSO(J+1)
     *       +2d0*DMVI(I,J)*DXYVO(J))
     * / (MO(I,J,1)*DXYNO(J) + MO(I,J+1,1)*DXYSO(J+1))
      END DO
      END DO
C**** Surface stress is applied to V component at the North Pole
      DO I=1,IM
        VO(I,JM-1,1) = VO(I,JM-1,1) +
     *       (DMVA(I,JM-1,1)*DXYNO(JM-1)+
     *       (DMVA(1,JM,1)*COSIC(I) - DMUA(1,JM,1)*SINIC(I))*DXYSO(JM)
     *       + 2d0*DMVI(I,JM-1)*DXYVO(JM-1)) /
     *  (MO(I,JM-1,1)*DXYNO(JM-1) + MO(I,JM,1)*DXYSO(JM))
      END DO
      RETURN
      END SUBROUTINE OSTRES

      SUBROUTINE GROUND_OC
!@sum  GROUND_OC adds vertical fluxes into the ocean
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : grav
      USE GEOM, only : dxyp,bydxyp
      USE OCEAN, only : im,jm,mo,g0m,s0m,focean,gzmo,imaxj,dxypo,bydxypo
     *     ,lmo,lmm,ratoc,rocat,opress
#ifdef TRACERS_OCEAN
     *     ,trmo,ntm
#endif
      USE FLUXES, only : solar,e0,evapor,dmsi,dhsi,dssi,runosi,erunosi
     *     ,flowo,eflowo,srunosi,apress,melti,emelti,smelti
#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
     *     ,trflowo,trevapor,trunosi,trmelti
#endif
     *     ,dtrsi
#endif
      USE SEAICE_COM, only : rsi
      IMPLICIT NONE
      INTEGER I,J
      REAL*8 DXYPJ,BYDXYPJ,RUNO,RUNI,ERUNO,ERUNI,SROX(2),G0ML(LMO)
     *     ,MO1,SO1,ROICE,DMOO,DMOI,DEOO,DEOI,GZML(LMO),SRUNO,SRUNI,DSOO
     *     ,DSOI,POCEAN,POICE
#ifdef TRACERS_OCEAN
      REAL*8, DIMENSION(NTM) :: TRUNO,TRUNI,DTROO,DTROI,TRO1
#endif
C****
C**** Add surface source of fresh water and heat
C****
      DO J=2,JM
        DXYPJ=DXYPO(J)
        BYDXYPJ=BYDXYPO(J)
      DO I=1,IMAXJ(J)
      IF(FOCEAN(I,J).gt.0.) THEN
        ROICE = RSI(I,J)
        POCEAN=FOCEAN(I,J)*(1.-ROICE)
        POICE =FOCEAN(I,J)*ROICE
C**** set mass & energy fluxes (incl. river/sea ice runoff + basal flux)
        RUNO = (FLOWO(I,J)+ MELTI(I,J))/(FOCEAN(I,J)*DXYPJ)-
     *                                   RATOC(J)*EVAPOR(I,J,1)
        RUNI = (FLOWO(I,J)+ MELTI(I,J))/(FOCEAN(I,J)*DXYPJ)+
     *                                   RATOC(J)*RUNOSI(I,J)
        ERUNO=(EFLOWO(I,J)+EMELTI(I,J))/(FOCEAN(I,J)*DXYPJ)+
     *                                   RATOC(J)*E0(I,J,1)
        ERUNI=(EFLOWO(I,J)+EMELTI(I,J))/(FOCEAN(I,J)*DXYPJ)+
     *                                   RATOC(J)*ERUNOSI(I,J)
        SRUNO=SMELTI(I,J)/(FOCEAN(I,J)*DXYPJ)
        SRUNI=SMELTI(I,J)/(FOCEAN(I,J)*DXYPJ)+RATOC(J)*SRUNOSI(I,J)
        G0ML(:) =  G0M(I,J,:)
        GZML(:) = GZMO(I,J,:)
        SROX(1)=SOLAR(1,I,J)*RATOC(J) ! open water
        SROX(2)=SOLAR(3,I,J)*RATOC(J) ! through ice
        MO1 = MO(I,J,1)
        SO1 = S0M(I,J,1)
#ifdef TRACERS_OCEAN
        TRO1(:) = TRMO(I,J,1,:)
#ifdef TRACERS_WATER
        TRUNO(:)=(TRFLOWO(:,I,J)+TRMELTI(:,I,J))/(FOCEAN(I,J)*DXYPJ)-
     *       RATOC(J)*TREVAPOR(:,1,I,J)
        TRUNI(:)=(TRFLOWO(:,I,J)+TRMELTI(:,I,J))/(FOCEAN(I,J)*DXYPJ)+
     *       RATOC(J)*TRUNOSI(:,I,J)
#else
        TRUNO(:)=0. ; TRUNI(:)=0.
#endif
#endif

        CALL OSOURC(ROICE,MO1,G0ML,GZML,SO1,DXYPJ,BYDXYPJ,LMM(I,J),RUNO
     *         ,RUNI,ERUNO,ERUNI,SRUNO,SRUNI,SROX,
#ifdef TRACERS_OCEAN
     *         TRO1,TRUNO,TRUNI,DTROO,DTROI,
#endif
     *         DMOO,DEOO,DMOI,DEOI,DSOO,DSOI)

C**** update ocean variables
          MO(I,J,1) = MO1
         S0M(I,J,1) = SO1
         G0M(I,J,:) = G0ML(:)
        GZMO(I,J,:) = GZML(:)

C**** Store mass and energy fluxes for formation of sea ice
        DMSI(1,I,J)=DMOO*ROCAT(J)
        DMSI(2,I,J)=DMOI*ROCAT(J)
        DHSI(1,I,J)=DEOO*ROCAT(J)
        DHSI(2,I,J)=DEOI*ROCAT(J)
        DSSI(1,I,J)=DSOO*ROCAT(J)
        DSSI(2,I,J)=DSOI*ROCAT(J)
#ifdef TRACERS_OCEAN
        TRMO(I,J,1,:) = TRO1(:)
        DTRSI(:,1,I,J)=DTROO(:)*ROCAT(J)
        DTRSI(:,2,I,J)=DTROI(:)*ROCAT(J)
#endif

C**** Calculate pressure anomaly at ocean surface (and scale for areas)
C**** Updated using latest sea ice (this ensures that total column mass
C**** is consistent for OGEOZ calculation). 
        OPRESS(I,J) = RATOC(J)*(APRESS(I,J))+GRAV*(
     *       (1.-RSI(I,J))*DMSI(1,I,J) + RSI(I,J)*DMSI(2,I,J))

        END IF
      END DO
      END DO

      OPRESS(2:IM,1)  = OPRESS(1,1)
      OPRESS(2:IM,JM) = OPRESS(1,JM)
C****
      RETURN
      END SUBROUTINE GROUND_OC

      SUBROUTINE OSOURC (ROICE,MO,G0ML,GZML,S0M,DXYPJ,BYDXYPJ,LMIJ,RUNO
     *     ,RUNI,ERUNO,ERUNI,SRUNO,SRUNI,SROX,
#ifdef TRACERS_OCEAN
     *     TROM,TRUNO,TRUNI,DTROO,DTROI,
#endif
     *     DMOO,DEOO,DMOI,DEOI,DSOO,DSOI)
!@sum  OSOURC applies fluxes to ocean in ice-covered and ice-free areas
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : shi,lhm
#ifdef TRACERS_OCEAN
      USE TRACER_COM, only : ntm,trname
#endif
      USE SW2OCEAN, only : lsrpd,fsr,fsrz
      USE SEAICE, only : fsss
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: ROICE,DXYPJ,BYDXYPJ,RUNO,RUNI,ERUNO,ERUNI
     *     ,SROX(2),SRUNO,SRUNI
      INTEGER, INTENT(IN) :: LMIJ
      REAL*8, INTENT(INOUT) :: MO,G0ML(LSRPD),GZML(LSRPD),S0M
      REAL*8, INTENT(OUT) :: DMOO,DMOI,DEOO,DEOI,DSOO,DSOI
#ifdef TRACERS_OCEAN
      REAL*8, DIMENSION(NTM), INTENT(INOUT) :: TROM
      REAL*8, DIMENSION(NTM), INTENT(IN) :: TRUNO,TRUNI
      REAL*8, DIMENSION(NTM), INTENT(OUT) :: DTROO,DTROI
      REAL*8, DIMENSION(NTM) :: TMOO,TMOI,FRAC
#ifdef TRACERS_SPECIAL_O18
      REAL*8 fracls
#endif
#endif
      REAL*8 MOO,GOO,GMOO,GMOI,MOI,GOI,SMOO,SMOI,SOO,SOI,GFOO,GFOI,TFOO
     *     ,TFOI,SIOO,SIOI
      REAL*8 GFREZS,TFREZS,TSOL
      INTEGER L,LSR,N

      DMOO=0. ; DEOO=0. ; DMOI=0. ; DEOI=0. ; DSOO=0. ; DSOI=0.
#ifdef TRACERS_OCEAN
      DTROI(:) = 0. ; DTROO(:) = 0.
      do n=1,ntm
#ifdef TRACERS_SPECIAL_O18
        FRAC(n)=fracls(trname(n))
#else
        FRAC(n)=1.
#endif
      end do
#endif

      LSR = MIN(LSRPD,LMIJ)
C****
C**** Open Ocean
C****
      MOO  = MO + RUNO
      GMOO = G0ML(1)*BYDXYPJ + ERUNO
      SMOO = S0M*BYDXYPJ + SRUNO
#ifdef TRACERS_OCEAN
      TMOO(:) = TROM(:)*BYDXYPJ+TRUNO(:)
#endif
      IF (ROICE.lt.1d0) THEN
C**** Remove insolation from layer 1 that goes to lower layers
      IF (LSR.gt.1) GMOO = GMOO - SROX(1)*FSR(2)

      GOO  = GMOO/MOO
      SOO  = SMOO/MOO
      GFOO = GFREZS(SOO)
      IF(GOO.lt.GFOO) THEN
C**** Open ocean is below freezing, calculate
C**** DMOO = mass of ocean that freezes over open fraction from
C**** GOO*MOO = GFOO*(MOO-DMOO) + (TFOO*SHI-LHM*(1-SIOO))*DMOO
        TFOO = TFREZS(SOO)
        SIOO = FSSS*SOO
        DMOO = MOO*(GOO-GFOO)/(TFOO*SHI-LHM*(1.-SIOO)-GFOO)
        DEOO = (TFOO*SHI-LHM*(1.-SIOO))*DMOO
        DSOO = SIOO*DMOO
#ifdef TRACERS_OCEAN
        DTROO(:) = TMOO(:)*FRAC(:)*(DMOO-DSOO)/(MOO-SMOO)
#endif
      END IF
      END IF
C****
C**** Ocean underneath the ice
C****
      MOI  = MO + RUNI
      GMOI = G0ML(1)*BYDXYPJ + ERUNI
      SMOI = S0M*BYDXYPJ + SRUNI
#ifdef TRACERS_OCEAN
      TMOI(:) = TROM(:)*BYDXYPJ+TRUNI(:)
#endif
      IF(ROICE.gt.0.) THEN
C**** Remove insolation from layer 1 that goes to lower layers
        IF (LSR.gt.1) GMOI = GMOI - SROX(2)*FSR(2)

        GOI  = GMOI/MOI
        SOI  = SMOI/MOI
        GFOI = GFREZS(SOI)
        IF(GOI.LT.GFOI) THEN
C**** Ocean underneath the ice is below freezing, calculate
C**** DMOI = mass of ocean that freezes under sea ice fraction from
C**** GOI*MOI = GFOI*(MOI-DMOI) + (TFOI*SHI-LHM*(1-SIOI))*DMOI
          TFOI = TFREZS(SOI)
          SIOI = FSSS*SOI
          DMOI = MOI*(GOI-GFOI)/(TFOI*SHI-LHM*(1.-SIOI)-GFOI)
          DEOI = (TFOI*SHI-LHM*(1.-SIOI))*DMOI
          DSOI = SIOI*DMOI
#ifdef TRACERS_OCEAN
          DTROI(:) = TMOI(:)*FRAC(:)*(DMOI-DSOI)/(MOI-SMOI)
#endif
        END IF
      END IF
C**** Update first layer variables
      MO     =  (MOI-DMOI)*ROICE + (1.-ROICE)*( MOO-DMOO)
      G0ML(1)=((GMOI-DEOI)*ROICE + (1.-ROICE)*(GMOO-DEOO))*DXYPJ
      S0M    =((SMOI-DSOI)*ROICE + (1.-ROICE)*(SMOO-DSOO))*DXYPJ
#ifdef TRACERS_OCEAN
      TROM(:)=((TMOI(:)-DTROI(:))*ROICE + (1.-ROICE)*(TMOO(:)-DTROO(:)))
     *     *DXYPJ
#endif
C**** add insolation to lower layers
      TSOL=(SROX(1)*(1.-ROICE)+SROX(2)*ROICE)*DXYPJ
      DO L=2,LSR-1
        G0ML(L)=G0ML(L)+TSOL*(FSR(L)-FSR(L+1))
        GZML(L)=GZML(L)+TSOL*FSRZ(L)
      END DO
      G0ML(LSR) = G0ML(LSR) + TSOL*FSR (LSR)
      GZML(LSR) = GZML(LSR) + TSOL*FSRZ(LSR)
C****
      RETURN
      END SUBROUTINE OSOURC

      SUBROUTINE PRECIP_OC
!@sum  PRECIP_OC driver for applying precipitation to ocean fraction
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE GEOM, only : dxyp
#ifdef TRACERS_OCEAN
      USE TRACER_COM, only : ntm
#ifdef TRACERS_WATER
      USE FLUXES, only : trpreca=>trprec,trunpsia=>trunpsi
#endif
#endif
      USE FLUXES, only : runpsia=>runpsi,srunpsia=>srunpsi,preca=>prec
     *     ,epreca=>eprec
      USE OCEAN, only : im,jm,mo,g0m,s0m,bydxypo,focean,imaxj
#ifdef TRACERS_OCEAN
     *     ,trmo,dxypo
#endif
      USE SEAICE_COM, only : rsia=>rsi
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,JM) :: PREC,EPREC,RUNPSI,RSI,SRUNPSI
#ifdef TRACERS_OCEAN
      REAL*8, DIMENSION(NTM,IM,JM) :: trprec,trunpsi
#endif
      INTEGER I,J

C**** save surface variables before any fluxes are added
      CALL KVINIT

C**** Convert fluxes on atmospheric grid to oceanic grid
C**** build in enough code to allow a different ocean grid.
C**** Since the geometry differs on B and C grids, some processing
C**** of fluxes is necessary anyway
      DO J=1,JM
        DO I=1,IMAXJ(J)
          PREC   (I,J)=PRECA   (I,J)*DXYP(J)*BYDXYPO(J)  ! kg/m^2
          EPREC  (I,J)=EPRECA  (I,J)*DXYP(J)             ! J
          RUNPSI (I,J)=RUNPSIA (I,J)*DXYP(J)*BYDXYPO(J)  ! kg/m^2
          SRUNPSI(I,J)=SRUNPSIA(I,J)*DXYP(J)             ! kg
          RSI    (I,J)=RSIA    (I,J)
#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
          TRPREC(:,I,J)=TRPRECA(:,I,J)                   ! kg
          TRUNPSI(:,I,J)=TRUNPSIA(:,I,J)*DXYP(J)         ! kg
#else
          TRPREC(:,I,J)=0.  ; TRUNPSI(:,I,J)=0.
#endif
#endif
        END DO
      END DO
C****
      DO J=1,JM
        DO I=1,IMAXJ(J)
          IF(FOCEAN(I,J).gt.0. .and. PREC(I,J).gt.0.)  THEN
            MO (I,J,1)= MO(I,J,1) + (1d0-RSI(I,J))*PREC(I,J) +
     *           RSI(I,J)*RUNPSI(I,J)
            G0M(I,J,1)=G0M(I,J,1) + (1d0-RSI(I,J))*EPREC(I,J)
            S0M(I,J,1)=S0M(I,J,1) + RSI(I,J)*SRUNPSI(I,J)
#ifdef TRACERS_OCEAN
            TRMO(I,J,1,:)=TRMO(I,J,1,:)+(1d0-RSI(I,J))*TRPREC(:,I,J)
     *             +RSI(I,J)*TRUNPSI(:,I,J)
#endif
          END IF
        END DO
      END DO

C**** Convert ocean surface temp to atmospheric SST array
      CALL TOC2SST

      RETURN
      END SUBROUTINE PRECIP_OC

      SUBROUTINE ODIFF (DTDIFF)
!@sum  ODIFF applies Wasjowicz horizontal viscosity to velocities
!@auth Gavin Schmidt
!@ver  1.0
C****
C**** ODIFF calculates horizontal Wasjowicz viscosity terms in momentum
C**** equations implicitly using ADI method and assumes no slip/free
C**** slip conditions at the side. K_h (m^2/s) may vary spatially
C**** based on Munk length though must remain isotropic.
C**** (If longitudinal variation is wanted just make K arrays K(I,J))
C**** FSLIP = 0 implies no slip conditions, = 1 implies free slip
C**** Mass variation is included
C****
      USE CONSTANT, only :twopi,rhows,omega,radius
      USE OCEAN, only : im,jm,lmo,mo,uo,vo,cospo,cosvo,rlat,lmu,lmv
     *     ,dxpo,dypo,dxvo,dyvo,dxyvo,dxypo,bydxypo
      USE OCEAN_DYN, only : dh
      IMPLICIT NONE
      REAL*8, PARAMETER :: AKHMIN=1.5d8, FSLIP=0.
      INTEGER, PARAMETER :: IIP=IM*(JM-2)+1
      REAL*8, DIMENSION(JM) :: KYPXP,KXPYV,KYVXV,KXVYP
      REAL*8, SAVE, DIMENSION(JM) :: BYDXYV,
     *     KHP,KHV,TANP,TANV,BYDXV,BYDXP,BYDYV,BYDYP
      REAL*8, DIMENSION(IM,JM,2) :: DUDX,DUDY,DVDX,DVDY
      REAL*8, DIMENSION(IM,JM) :: FUX,FUY,FVX,FVY,BYMU,BYMV
      INTEGER, SAVE :: IFIRST = 1
C**** Local variables
      REAL*8, DIMENSION(IIP) :: AU,BU,CU,RU,UU,AV,BV,CV,RV,UV
      REAL*8, SAVE, DIMENSION(IM,JM,LMO) :: UXA,UXB,UXC,UYA,UYB,UYC,VXA
     *     ,VXB,VXC,VYA,VYB,VYC
      REAL*8, SAVE, DIMENSION(LMO) :: UYPB
      REAL*8, SAVE, DIMENSION(IM,LMO) :: UYPA
      REAL*8, INTENT(IN) :: DTDIFF
      REAL*8, SAVE :: BYDXYPJM
      REAL*8 DSV,DSP,VLAT,DLAT,DT2,DTU,DTV,VX,VY,VT,UT,UX,UY
      INTEGER I,J,L,IP1,IM1,II
C****
      IF(IFIRST.ne.0)  THEN
      IFIRST = 0
      DO J=2,JM-1
C**** Calculate KH = rho_0 BETA* L_Munk^3 where DX=L_Munk
c      KHP(J)=2d0*RHOWS*OMEGA*COSP(J)*(DXP(J)**3)/RADIUS ! tracer lat
c      KHV(J)=2d0*RHOWS*OMEGA*COSV(J)*(DXV(J)**3)/RADIUS ! (v vel pts)
C**** Calculate KH=rho_0 BETA (sqrt(3) L_Munk/pi)^3, L_Munk=min(DX,DY)
        DSP=MIN(DXPO(J),DYPO(J))*2.*SQRT(3.)/TWOPI  ! tracer lat
        DSV=MIN(DXVO(J),DYVO(J))*2.*SQRT(3.)/TWOPI  ! v vel pts
        KHP(J)=2d0*RHOWS*OMEGA*COSPO(J)*(DSP**3)/RADIUS ! tracer lat
        KHV(J)=2d0*RHOWS*OMEGA*COSVO(J)*(DSV**3)/RADIUS ! (v vel pts)
        KHP(J)=MAX(KHP(J),AKHMIN)
        KHV(J)=MAX(KHV(J),AKHMIN)
        BYDXYV(J)=1D0/DXYVO(J)
        BYDXV(J)=1D0/DXVO(J)
        BYDXP(J)=1D0/DXPO(J)
        BYDYV(J)=1D0/DYVO(J)
        BYDYP(J)=1D0/DYPO(J)
        KYPXP(J)=KHP(J)*DYPO(J)*BYDXP(J)
        KXPYV(J)=KHV(J)*DXPO(J)*BYDYV(J)
        KYVXV(J)=KHV(J)*DYVO(J)*BYDXV(J)
        KXVYP(J)=KHP(J)*DXVO(J)*BYDYP(J)
C**** Discretisation errors need TANP/V to be defined like this
        DLAT = TWOPI*NINT(360d0/(JM-1))/720d0
        TANP(J)=TAN(RLAT(J))*TAN(0.5*DLAT)/(RADIUS*0.5*DLAT)
        VLAT = DLAT*(J+0.5-0.5*(1+JM))
        TANV(J)=TAN(VLAT)*SIN(DLAT)/(DLAT*RADIUS)
      END DO
      KHV(1)=KHV(2)
      BYDXV(1)=1D0/DXVO(1)
      BYDYV(1)=1D0/DYVO(1)
      BYDYP(1)=1D0/DYPO(1)
      BYDXP(JM)=1D0/DXPO(JM)
      BYDXYPJM=1D0/(DXYPO(JM)*IM)
      PRINT*,"Kh: ",(J,KHP(J),J=1,JM)
C****
C**** Calculate operators fixed in time for U and V equations
C****
      UXA=0. ; UXB=0. ; UXC=0. ; UYA=0. ; UYB=0. ; UYC=0.
      VXA=0. ; VXB=0. ; VXC=0. ; VYA=0. ; VYB=0. ; VYC=0.
      UYPB=0.; UYPA=0.

      DO L=1,LMO
C**** Calculate flux operators
C**** i.e DUDX(1) and DUDX(2) are the coefficients of u1 and u2 for
C**** calculating centered difference K_h du/dx
C**** including metric terms in y derivatives
        DUDX=0.
        DUDY=0.
        DVDX=0.
        DVDY=0.
        DO J=2,JM-1
          I=IM
          DO IP1=1,IM
            IF (L.LE.LMU(IP1,J)) DUDX(IP1,J,1) = KYPXP(J)
            IF (L.LE.LMU(I,J)) THEN
              DUDX(IP1,J,2) = -KYPXP(J)
              IF (L.LE.LMU(I,J+1)) THEN
                DUDY(I,J,1) =  KXPYV(J)*(1. +0.5*TANV(J)*DYVO(J))
                DUDY(I,J,2) = -KXPYV(J)*(1. -0.5*TANV(J)*DYVO(J))
              ELSE
                DUDY(I,J,2) = -(1.-FSLIP)*2d0*KXPYV(J)
              END IF
            ELSE
              IF (L.LE.LMU(I,J+1)) DUDY(I,J,1) = (1.-FSLIP)*2d0*KXPYV(J)
            END IF
            IF (L.LE.LMV(I,J+1)) DVDY(I,J+1,1) = KXVYP(J)*(1. +
     *           0.5*TANP(J)*DYPO(J))
            IF (L.LE.LMV(I,J)) THEN
              DVDY(I,J+1,2) = -KXVYP(J)*(1.-0.5*TANP(J)*DYPO(J))
              IF (L.LE.LMV(IP1,J)) THEN
                DVDX(I,J,1) =  KYVXV(J)
                DVDX(I,J,2) = -KYVXV(J)
              ELSE
                DVDX(I,J,2) = -(1.-FSLIP)*2d0*KYVXV(J)
              END IF
            ELSE
              IF (L.LE.LMV(IP1,J)) DVDX(I,J,1) = (1.-FSLIP)*2d0*KYVXV(J)
            END IF
            I=IP1
          END DO
        END DO
C****
C**** Combine to form tri-diagonal operators including first metric term
C****
        DO J=2,JM-1
          IM1=IM
          DO I=1,IM
            IF(L.LE.LMU(IM1,J)) THEN
              UXA(IM1,J,L) = -DUDX(IM1,J,2)                  *BYDXYPO(J)
              UXB(IM1,J,L) = (DUDX(I  ,J,2) -DUDX(IM1,J  ,1))*BYDXYPO(J)
              UXC(IM1,J,L) =  DUDX(I  ,J,1)                  *BYDXYPO(J)
              UYA(IM1,J,L) = -DUDY(IM1,J-1,2)                *BYDXYPO(J)
     *                     + 0.5*TANP(J)*KHP(J)*BYDYP(J)
              UYB(IM1,J,L) =(DUDY(IM1,J  ,2)-DUDY(IM1,J-1,1))*BYDXYPO(J)
     *                     + TANP(J)*TANP(J)*KHP(J)
              UYC(IM1,J,L) = DUDY(IM1,J  ,1)                 *BYDXYPO(J)
     *                     - 0.5*TANP(J)*KHP(J)*BYDYP(J)
            END IF
            IF (L.LE.LMV(I,J)) THEN
              VXA(I,J,L) = -DVDX(IM1,J,2)                 *BYDXYV(J)
              VXB(I,J,L) = (DVDX(I  ,J,2) - DVDX(IM1,J,1))*BYDXYV(J)
              VXC(I,J,L) =  DVDX(I  ,J,1)                 *BYDXYV(J)
              VYA(I,J,L) = -DVDY(I,J  ,2)                 *BYDXYV(J)
     *                   + 0.5*TANV(J)*KHV(J)*BYDYV(J)
              VYB(I,J,L) = (DVDY(I,J+1,2) - DVDY(I  ,J,1))*BYDXYV(J)
     *                   + TANV(J)*TANV(J)*KHV(J)
              VYC(I,J,L) =  DVDY(I,J+1,1)                 *BYDXYV(J)
     *                   - 0.5*TANV(J)*KHV(J)*BYDYV(J)
            END IF
            IM1=I
          END DO
        END DO
C**** At North Pole
        IF(L.LE.LMU(1,JM)) THEN
          DO I=1,IM
            UYPB(L) = UYPB(L) - DUDY(I,JM-1,1)
            UYPA(I,L) = -DUDY(I,JM-1,2)*BYDXYPJM
          END DO
          UYPB(L) = UYPB(L)*BYDXYPJM
        END IF
      END DO
      END IF
C**** Solve diffusion equations semi implicitly
      DT2=DTDIFF*5d-1     ! half time-step
!$OMP PARALLEL DO  PRIVATE(AU,AV, BYMU,BYMV,BU,BV, CU,CV, DTU,DTV,
!$OMP&  FUX,FUY,FVX,FVY, I,IP1,IM1,II, J, L, RU,RV,
!$OMP&  UU,UV,UT,UY,UX, VT,VY,VX)
      DO L=1,LMO
C**** Save (0.5*) mass reciprical for velocity points
      DO J=2,JM-1
        I=IM
        DO IP1=1,IM
          IF (L.LE.LMU(I,J)) BYMU(I,J) = 1./(MO(I,J,L)+MO(IP1,J,L))
          IF (L.LE.LMV(I,J)) BYMV(I,J) = 1./(MO(I,J,L)+MO(I,J+1,L))
          I=IP1
        END DO
      END DO
      IF (L.LE.LMU(1,JM)) BYMU(1,JM) = 1./MO(1,JM,L)
C**** Calculate Wasjowicz boundary terms
C**** Need dv/dy,tv,dv/dx for u equation, du/dy,tu,du/dx for v equation
      FUX=0             ! flux in U equation at the x_+ boundary
      FUY=0             ! flux in U equation at the y_+ boundary
      FVX=0             ! flux in V equation at the x_+ boundary
      FVY=0             ! flux in V equation at the y_+ boundary
      DO J=1,JM-1
        IM1=IM-1
        DO I=1,IM
          UT=0          ! mean u*tan on x_+ boundary for V equation
          UY=0          ! mean du/dx on y_+ boundary for V equation
          UX=0          ! mean du/dy on x_+ boundary for V equation
          IF (L.LE.LMU(I  ,J+1)) THEN
            UT=     UO(I  ,J+1,L)*TANP(J+1)
            UX=     UO(I  ,J+1,L)
            UY=     UO(I  ,J+1,L)
          END IF
          IF (L.LE.LMU(I  ,J  )) THEN
            UT=UT + UO(I  ,J  ,L)*TANP(J  )
            UY=UY - UO(I  ,J  ,L)
          END IF
          IF (L.LE.LMU(IM1,J+1)) UX=UX-    UO(IM1,J+1,L)
          UT=0.5*UT
          UX=UX*BYDXP(J+1)
          UY=UY*BYDYV(J)
C****
          VT=0          ! mean v*tan on x_+ boundary for U equation
          VX=0          ! mean dv/dx on y_+ boundary for U equation
          VY=0          ! mean dv/dy on x_+ boundary for U equation
          IF (L.LE.LMV(I  ,J  )) THEN
            VT=     VO(I  ,J  ,L)*TANV(J  )
            VX=     VO(I  ,J  ,L)
            VY=     VO(I  ,J  ,L)
          END IF
          IF (J.GT.1 .and. L.LE.LMV(I  ,J-1)) THEN
            VT=VT + VO(I  ,J-1,L)*TANV(J-1)
            VY=VY - VO(I  ,J-1,L)
          END IF
          IF (L.LE.LMV(IM1,J  )) VX=VX - VO(IM1,J  ,L)
          VT=0.5*VT
          VY=VY*BYDYP(J)
          VX=VX*BYDXV(J)
C**** Calculate fluxes (including FSLIP condition)
          IF (FSLIP.EQ.1.) THEN
            IF (L.LE.LMV(I,J) .AND. L.LE.LMV(IM1,J))
     *           FUY(IM1,J)=KHV(J)*VX
            IF (L.LE.LMU(I,J) .AND. L.LE.LMU(I,J+1))
     *           FVX(I  ,J)=KHV(J)*(UY + UT)
          ELSE
            FUY(IM1,J)=KHV(J)*VX
            FVX(I  ,J)=KHV(J)*(UY + UT)
          END IF
          FUX(IM1,J)=KHP(J)*(VY + VT)
          IF (J.LT.JM-1) FVY(I,J)=KHP(J+1)*UX
          IM1=I
        END DO
      END DO
C**** Calculate tridiagonal matrix for first semi-implicit step (in x)
C**** Minor complication due to cyclic nature of boundary condition
      AU=0. ; BU=0. ; CU=0. ; RU=0. ; UU=0.
      AV=0. ; BV=0. ; CV=0. ; RV=0. ; UV=0.
      DO J=2,JM-1
        IM1=IM-1
        I=IM
        DO IP1=1,IM
          II=IM*(J-2)+I
          BU(II) = 1d0
          BV(II) = 1d0
          IF (L.LE.LMU(I,J)) THEN
            DTU = DT2*(DH(I,J,L)+DH(IP1,J,L))*BYMU(I,J)
            IF (I.gt.1 ) AU(II) =        - DTU*UXA(I,J,L)
                         BU(II) = BU(II) - DTU*UXB(I,J,L)
            IF (I.lt.IM) CU(II) =        - DTU*UXC(I,J,L)
            RU(II) = UO(I,J,L) + DTU*(UYA(I,J,L)*UO(I,J-1,L)
     *           +UYB(I,J,L)*UO(I,J,L) + UYC(I,J,L)*UO(I,J+1,L))
C**** Make properly tridiagonal by making explicit cyclic terms
            IF (I.eq.1 ) RU(II)=RU(II) + DTU*UXA(I,J,L)*UO(IM,J,L)
            IF (I.eq.IM) RU(II)=RU(II) + DTU*UXC(I,J,L)*UO(1,J,L)
C**** Add Wasjowicz cross-terms to RU + second metric term
            RU(II) = RU(II) + DTU*((DYPO(J)*(FUX(IM1,J) - FUX(I,J))
     *           + DXVO(J)*FUY(I,J) - DXVO(J-1)*FUY(I,J-1))*BYDXYPO(J)
     *           - 0.5*(TANV(J-1)*FUY(I,J-1) + TANV(J)*FUY(I,J)))
          END IF
          IF (L.LE.LMV(I,J)) THEN
            DTV = DT2*(DH(I,J,L)+DH(I,J+1,L))*BYMV(I,J)
            IF (I.gt.1 ) AV(II) =        - DTV*VXA(I,J,L)
                         BV(II) = BV(II) - DTV*VXB(I,J,L)
            IF (I.lt.IM) CV(II) =        - DTV*VXC(I,J,L)
            RV(II) = VO(I,J,L) + DTV*(VYA(I,J,L)*VO(I,J-1,L)
     *           +VYB(I,J,L)*VO(I,J,L) + VYC(I,J,L)*VO(I,J+1,L))
C**** Make properly tridiagonal by making explicit cyclic terms
            IF (I.eq.1 ) RV(II)=RV(II) + DTV*VXA(I,J,L)*VO(IM,J,L)
            IF (I.eq.IM) RV(II)=RV(II) + DTV*VXC(I,J,L)*VO(1,J,L)
C**** Add Wasjowicz cross-terms to RV + second metric term
            RV(II) = RV(II) + DTV*((DYVO(J)*(FVX(I,J) - FVX(IM1,J))
     *           + DXPO(J)*FVY(I,J-1) - DXPO(J+1)*FVY(I,J))*BYDXYV(J)
     *           + 0.5*(TANP(J-1)*FVY(I,J-1) + TANP(J)*FVY(I,J)))
          END IF
          IM1=I
          I=IP1
        END DO
      END DO
C**** At North Pole (no metric terms)
      BU(IIP) = 1d0
      BV(IIP) = 1d0
      IF (L.LE.LMU(1,JM)) THEN
      DTU = DT2*DH(1,JM,L)*BYMU(1,JM)
        RU(IIP) = UO(1,JM,L) +DTU*UYPB(L)*UO(1,JM,L)
        DO I=1,IM       ! include Wasjowicz cross-terms at North Pole
          RU(IIP) = RU(IIP) + DTU*(UYPA(I,L)*UO(I,JM-1,L)
     *                      - DXVO(JM-1)*FUY(I,JM-1)*BYDXYPJM)
        END DO
      END IF
C**** Call tridiagonal solver
      CALL TRIDIAG(AU,BU,CU,RU,UU,IIP)
      CALL TRIDIAG(AV,BV,CV,RV,UV,IIP)
      DO II=1,IIP
c       J= 2 + (II-1)/IM
c       I= II-(J-2)*IM
c       UO(I,J,L) = UU(II)
c       VO(I,J,L) = UV(II)
        UO(II,2,L) = UU(II)   ! this cycles through correctly
        VO(II,2,L) = UV(II)
      END DO
      DO I=2,IM
        UO(I,JM,L) = UO(1,JM,L)
      END DO
C**** Now do semi implicit solution in y
C**** Calc. cross-term fluxes + second metric term (at half time step)
C**** Need dv/dy,tv,dv/dx for u equation, du/dy,tu,du/dx for v equation
      FUX=0             ! flux in U equation at the x_+ boundary
      FUY=0             ! flux in U equation at the y_+ boundary
      FVX=0             ! flux in V equation at the x_+ boundary
      FVY=0             ! flux in V equation at the y_+ boundary
      DO J=1,JM-1
        IM1=IM-1
        DO I=1,IM
          UT=0         ! mean u*tan on x_+ boundary for V equation
          UY=0         ! mean du/dx on y_+ boundary for V equation
          UX=0         ! mean du/dy on x_+ boundary for V equation
          IF (L.LE.LMU(I  ,J+1)) THEN
            UT=     UO(I  ,J+1,L)*TANP(J+1)
            UX=     UO(I  ,J+1,L)
            UY=     UO(I  ,J+1,L)
          END IF
          IF (L.LE.LMU(I  ,J  )) THEN
            UT=UT + UO(I  ,J  ,L)*TANP(J  )
            UY=UY - UO(I  ,J  ,L)
          END IF
          IF (L.LE.LMU(IM1,J+1)) UX=UX-    UO(IM1,J+1,L)
          UT=0.5*UT
          UX=UX*BYDXP(J+1)
          UY=UY*BYDYV(J)
C****
          VT=0         ! mean v*tan on x_+ boundary for U equation
          VX=0         ! mean dv/dx on y_+ boundary for U equation
          VY=0         ! mean dv/dy on x_+ boundary for U equation
          IF (L.LE.LMV(I  ,J  )) THEN
            VT=     VO(I  ,J  ,L)*TANV(J  )
            VX=     VO(I  ,J  ,L)
            VY=     VO(I  ,J  ,L)
          END IF
          IF (J.GT.1 .and. L.LE.LMV(I  ,J-1)) THEN
            VT=VT + VO(I  ,J-1,L)*TANV(J-1)
            VY=VY - VO(I  ,J-1,L)
          END IF
          IF (L.LE.LMV(IM1,J  )) VX=VX - VO(IM1,J  ,L)
          VT=0.5*VT
          VY=VY*BYDYP(J)
          VX=VX*BYDXV(J)
C**** Calculate fluxes (including FSLIP condition)
          IF (FSLIP.EQ.1.) THEN
            IF (L.LE.LMV(I,J) .AND. L.LE.LMV(IM1,J))
     *           FUY(IM1,J)=KHV(J)* VX
            IF (L.LE.LMU(I,J) .AND. L.LE.LMU(I,J+1))
     *           FVX(I,J)=KHV(J)*(UY + UT)
          ELSE
            FUY(IM1,J)=KHV(J)* VX
            FVX(I  ,J)=KHV(J)*(UY + UT)
          END IF
          FUX(IM1,J)=KHP(J)*(VY + VT)
          IF (J.LT.JM-1) FVY(I,J)=KHP(J+1)*UX
          IM1=I
        END DO
      END DO
C**** Calculate tridiagonal matrix for second semi-implicit step (in y)
C**** Minor complication due to singular nature of polar box
      AU=0. ; BU=0. ; CU=0. ; RU=0. ; UU=0.
      AV=0. ; BV=0. ; CV=0. ; RV=0. ; UV=0.
      IM1=IM-1
      I=IM
      DO IP1=1,IM
        DO J=2,JM-1
          II=(JM-2)*(I-1)+(J-1)
          BU(II) = 1d0
          BV(II) = 1d0
          IF (L.LE.LMU(I,J)) THEN
            DTU = DT2*(DH(I,J,L)+DH(IP1,J,L))*BYMU(I,J)
            AU(II) =        - DTU*UYA(I,J,L)
            BU(II) = BU(II) - DTU*UYB(I,J,L)
            IF (J.lt.JM-1) CU(II) =     - DTU*UYC(I,J,L)
            RU(II) = UO(I,J,L) + DTU*(UXA(I,J,L)*UO(IM1,J,L)
     *           +UXB(I,J,L)*UO(I,J,L) + UXC(I,J,L)*UO(IP1,J,L))
C**** Make properly tridiagonal by making explicit polar terms
            IF (J.lt.JM-1) RU(II)=RU(II) + DTU*UYC(I,J,L)*UO(1,JM,L)
C**** Add Wasjowicz cross-terms to RU + second metric term
            RU(II) = RU(II) + DTU*((DYPO(J)*(FUX(IM1,J) - FUX(I,J))
     *           + DXVO(J)*FUY(I,J) - DXVO(J-1)*FUY(I,J-1))*BYDXYPO(J)
     *           - 0.5*(TANV(J-1)*FUY(I,J-1) + TANV(J)*FUY(I,J)))
          END IF
          IF (L.LE.LMV(I,J)) THEN
            DTV = DT2*(DH(I,J,L)+DH(I,J+1,L))*BYMV(I,J)
            AV(II) =        - DTV*VYA(I,J,L)
            BV(II) = BV(II) - DTV*VYB(I,J,L)
            IF (J.lt.JM-1) CV(II) =     - DTV*VYC(I,J,L)
            RV(II) = VO(I,J,L) + DTV*(VXA(I,J,L)*VO(IM1,J,L)
     *           +VXB(I,J,L)*VO(I,J,L) + VXC(I,J,L)*VO(IP1,J,L))
C**** Add Wasjowicz cross-terms to RV + second metric term
            RV(II) = RV(II) + DTV*((DYVO(J)*(FVX(I,J) - FVX(IM1,J))
     *           + DXPO(J)*FVY(I,J-1) - DXPO(J+1)*FVY(I,J))*BYDXYV(J)
     *           + 0.5*(TANP(J-1)*FVY(I,J-1) + TANP(J)*FVY(I,J)))
          END IF
          IM1=I
          I=IP1
        END DO
      END DO
C**** At North Pole (do partly explicitly) no metric terms
      BU(IIP) = 1d0
      BV(IIP) = 1d0
      IF (L.LE.LMU(1,JM)) THEN
        DTU = DT2*DH(1,JM,L)*BYMU(1,JM)
        BU(IIP) = BU(IIP) - DTU*UYPB(L)
        RU(IIP) = UO(1,JM,L)
        DO I=1,IM       ! include Wasjowicz cross-terms at North Pole
          RU(IIP)= RU(IIP) + DTU*(UYPA(I,L)*UO(I,JM-1,L)
     *         - DXVO(JM-1)*FUY(I,JM-1)*BYDXYPJM)
        END DO
      END IF
C**** Call tridiagonal solver
      CALL TRIDIAG(AU,BU,CU,RU,UU,IIP)
      CALL TRIDIAG(AV,BV,CV,RV,UV,IIP)
      DO II=1,IIP-1
        I= 1 + (II-1)/(JM-2)
        J= 1 + II - (I-1)*(JM-2)
        UO(I,J,L) = UU(II)
        VO(I,J,L) = UV(II)
      END DO
      DO I=1,IM
        UO(I,JM,L) = UU(IIP)
      END DO
C****
      END DO
!$OMP END PARALLEL DO
C**** Done!
      RETURN
      END SUBROUTINE ODIFF

      SUBROUTINE TOC2SST
!@sum  TOC2SST convert ocean surface variables into atmospheric sst
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : byshi,lhm
#ifdef TRACERS_WATER
      USE TRACER_COM, only : trw0
#endif
      USE OCEAN, only : im,jm,imaxj,g0m,s0m,mo,dxypo,focean,lmm,ogeoz,uo
     *     ,vo,ogeoz_sv
#ifdef TRACERS_OCEAN
     *     ,trmo
#endif
      USE FLUXES, only : gtemp, sss, mlhc, ogeoza, uosurf, vosurf
#ifdef TRACERS_WATER
     *     ,gtracer
#endif
      IMPLICIT NONE
      INTEGER I,J
      REAL*8 TEMGS,shcgs,GO,SO,GO2,SO2,TO
C****
C**** Note that currently everything is on same grid
C****
      DO J=1,JM
        DO I=1,IMAXJ(J)
          IF (FOCEAN(I,J).gt.0.) THEN
            GO= G0M(I,J,1)/(MO(I,J,1)*DXYPO(J))
            SO= S0M(I,J,1)/(MO(I,J,1)*DXYPO(J))
            TO= TEMGS(GO,SO)
            GTEMP(1,1,I,J)= TO
            SSS(I,J) = 1d3*SO
            MLHC(I,J)= MO(I,J,1)*SHCGS(GO,SO)
            IF (LMM(I,J).gt.1) THEN
              GO2= G0M(I,J,2)/(MO(I,J,2)*DXYPO(J))
              SO2= S0M(I,J,2)/(MO(I,J,2)*DXYPO(J))
              TO= TEMGS(GO2,SO2)
            END IF
            GTEMP(2,1,I,J)= TO
   ! atmospheric grid Ocean height
            OGEOZA(I,J)=0.5*(OGEOZ(I,J)+OGEOZ_SV(I,J))
            UOSURF(I,J)=UO(I,J,1)
            VOSURF(I,J)=VO(I,J,1)
#ifdef TRACERS_WATER
#ifdef TRACERS_OCEAN
            GTRACER(:,1,I,J)=TRMO(I,J,1,:)/(MO(I,J,1)*DXYPO(J)-
     *           S0M(I,J,1))
#else
            GTRACER(:,1,I,J)=trw0(:)
#endif
#endif
          END IF
        END DO
      END DO
C**** do poles
      IF (FOCEAN(1,JM).gt.0) THEN
        VOSURF(1,JM)=0.
        DO I=2,IM
          GTEMP(:,1,I,JM)=GTEMP(:,1,1,JM)
          SSS(I,JM)=SSS(1,JM)
          MLHC(I,JM)=MLHC(1,JM)
          UOSURF(I,JM)=UOSURF(1,JM)
          VOSURF(I,JM)=0.
          OGEOZA(I,JM)=OGEOZA(1,JM)
#ifdef TRACERS_WATER
          GTRACER(:,1,I,JM)=GTRACER(:,1,1,JM)
#endif
        END DO
      END IF
c     IF (FOCEAN(1,1).gt.0) THEN
c       DO I=2,IM
c         GTEMP(:,1,I,1)=GTEMP(:,1,1,1)
c         SSS(I,1)=SSS(1,1)
c         MLHC(I,1)=MLHC(1,1)
#ifdef TRACERS_WATER
c         GTRACER(:,1,I,1)=GTRACER(:,1,1,1)
#endif
c       END DO
c     END IF
      RETURN
C****
      END SUBROUTINE TOC2SST

      SUBROUTINE io_oda(kunit,it,iaction,ioerr)
!@sum  io_oda dummy routine for consistency with uncoupled model
!@auth Gavin Schmidt
!@ver  1.0
      RETURN
      END SUBROUTINE io_oda

      SUBROUTINE ADVSI_DIAG
!@sum ADVSI_DIAG dummy routine for consistency with qflux model
      RETURN
      END

      SUBROUTINE AT2OT(FIELDA,FIELDO,NF,QCONSERV)
!@sum  AT2OT interpolates Atm Tracer grid to Ocean Tracer grid
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ima=>im,jma=>jm
      USE OCEAN, only : imo=>im,jmo=>jm,imaxj,ratoc
      IMPLICIT NONE
!@var QCONSERV true if integrated field must be conserved
      LOGICAL, INTENT(IN) :: QCONSERV
!@var N number of fields
      INTEGER, INTENT(IN) :: NF
!@var FIELDA array on atmospheric tracer grid
      REAL*8, INTENT(IN), DIMENSION(NF,IMA,JMA) :: FIELDA
!@var FIELDO array on oceanic tracer grid
      REAL*8, INTENT(OUT), DIMENSION(NF,IMO,JMO) :: FIELDO
      INTEGER I,J

C**** currently no need for interpolation,
C**** just scaling due to area differences for fluxes
      IF (QCONSERV) THEN
        DO J=1,JMO
          DO I=1,IMAXJ(J)
            FIELDO(:,I,J) = FIELDA(:,I,J)*RATOC(J)
          END DO
        END DO
        DO I=2,IMO
          FIELDO(:,I,JMO)=FIELDO(:,1,JMO)
          FIELDO(:,I,1  )=FIELDO(:,1,  1)
        END DO
      ELSE
        FIELDO = FIELDA
      END IF
C****
      RETURN
      END SUBROUTINE AT2OT

      SUBROUTINE OT2AT(FIELDO,FIELDA,NF,QCONSERV)
!@sum  OT2AT interpolates Ocean Tracer grid to Atm Tracer grid
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ima=>im,jma=>jm
      USE GEOM, only : imaxj
      USE OCEAN, only : imo=>im,jmo=>jm,rocat
      IMPLICIT NONE
!@var QCONSERV true if integrated field must be conserved
      LOGICAL, INTENT(IN) :: QCONSERV
!@var N number of fields
      INTEGER, INTENT(IN) :: NF
!@var FIELDA array on atmospheric tracer grid
      REAL*8, INTENT(OUT), DIMENSION(NF,IMA,JMA) :: FIELDA
!@var FIELDO array on oceanic tracer grid
      REAL*8, INTENT(IN), DIMENSION(NF,IMO,JMO) :: FIELDO
      REAL*8 RAT
      INTEGER I,J

C**** currently no need for interpolation,
C**** just scaling due to area differences for fluxes
      IF (QCONSERV) THEN
        DO J=1,JMA
          DO I=1,IMAXJ(J)
            FIELDA(:,I,J) = FIELDO(:,I,J)*ROCAT(J)
          END DO
        END DO
        DO I=2,IMA
          FIELDA(:,I,JMA)=FIELDA(:,1,JMA)
          FIELDA(:,I,1  )=FIELDA(:,1,  1)
        END DO
      ELSE
        FIELDA = FIELDO
      END IF
C****
      RETURN
      END SUBROUTINE OT2AT

      SUBROUTINE AT2OV(FIELDA,FIELDO,NF,QCONSERV,QU)
!@sum  OT2AT interpolates Atm Tracer grid to Ocean Velocity grid
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ima=>im,jma=>jm
      USE GEOM, only : imaxj
      USE OCEAN, only : imo=>im,jmo=>jm,ramvn,ramvs,ratoc
      IMPLICIT NONE
!@var QCONSERV true if integrated field must be conserved
!@var QU true if u-velocity pts. are wanted (false for v velocity pts.)
      LOGICAL, INTENT(IN) :: QCONSERV, QU
!@var N number of fields
      INTEGER, INTENT(IN) :: NF
!@var FIELDA array on atmospheric tracer grid
      REAL*8, INTENT(IN), DIMENSION(NF,IMA,JMA) :: FIELDA
!@var FIELDO array on oceanic tracer grid
      REAL*8, INTENT(OUT), DIMENSION(NF,IMO,JMO) :: FIELDO
      INTEGER I,J,IP1

C**** Ocean velocities are on C grid
      IF (QU) THEN  ! interpolate onto U points
        IF (QCONSERV) THEN
          DO J=1,JMO-1
            I=IMO
            DO IP1=1,IMO
              FIELDO(:,I,J)=0.5*RATOC(J)*(FIELDA(:,I,J)+FIELDA(:,IP1,J))
              I=IP1
            END DO
          END DO
C**** do poles
          DO I=1,IMO
            FIELDO(:,I,JMO) = FIELDA(:,1,JMO)*RATOC(JMO)
            FIELDO(:,I,  1) = FIELDA(:,1,  1)*RATOC(JMO)
          END DO
        ELSE   ! no area weighting
          DO J=1,JMO-1
            I=IMO
            DO IP1=1,IMO
              FIELDO(:,I,J) = 0.5*(FIELDA(:,I,J)+FIELDA(:,IP1,J))
              I=IP1
            END DO
          END DO
C**** do poles
          DO I=1,IMO
            FIELDO(:,I,JMO) = FIELDO(:,1,JMO)
            FIELDO(:,I,  1) = FIELDO(:,1,  1)
          END DO
        END IF
      ELSE                      ! interpolate onto V points
        IF (QCONSERV) THEN
          DO J=1,JMO-1
            DO I=1,IMO
              FIELDO(:,I,J) = RATOC(J)*RAMVN(J)*FIELDA(:,I,J)+
     *             RATOC(J+1)*RAMVS(J+1)*FIELDA(:,I,J+1)
            END DO
          END DO
        ELSE
          DO J=1,JMO-1
            DO I=1,IMO
              FIELDO(:,I,J) = 0.5*(FIELDA(:,I,J)+FIELDA(:,I,J+1))
            END DO
          END DO
        END IF
        FIELDO(:,:,JMO) = 0.
      END IF
C****
      RETURN
      END SUBROUTINE AT2OV

      SUBROUTINE GLMELT(DT)
!@sum  GLMELT adds glacial melt around Greenland and Antarctica to ocean
!@auth Sukeshi Sheth/Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : dtsrc
      USE OCEAN, only : imaxj,jm,g0m,mo,ze,focean,dxypo
#ifdef TRACERS_OCEAN
     *     ,trmo
#endif
      USE FLUXES, only : gmelt,egmelt
#ifdef TRACERS_WATER
     *     ,trgmelt
#endif
      IMPLICIT NONE
!@var MAXL number of ocean levels over which GMELT is applied
      INTEGER, PARAMETER :: MAXL=5
      REAL*8, INTENT(IN) :: DT  !@var DT timestep for GLMELT call
      REAL*8 DZ
      INTEGER I,J,L

      DO L=1,MAXL
C**** divide over depth and scale for time step
        DZ=DT*(ZE(L)-ZE(L-1))/(DTsrc*ZE(MAXL))
        DO J=1,JM
          DO I=1,IMAXJ(J)
            IF (FOCEAN(I,J).GT.0. .and. GMELT(I,J).gt.0) THEN
              MO(I,J,L) = MO(I,J,L)+GMELT(I,J)*DZ/(DXYPO(J)*FOCEAN(I,J))
              G0M(I,J,L)=G0M(I,J,L)+EGMELT(I,J)*DZ
#ifdef TRACERS_OCEAN
              TRMO(I,J,L,:)=TRMO(I,J,L,:)+TRGMELT(:,I,J)*DZ
#endif
            END IF
          END DO
        END DO
      END DO
C****
      RETURN
      END SUBROUTINE GLMELT
