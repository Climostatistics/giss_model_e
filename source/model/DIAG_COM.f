      MODULE DAGCOM
!@sum  DAGCOM Diagnostic model variables
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : twopi
      USE MODEL_COM, only : im,jm,lm,imh,fim,ntype,kep,istrat
      USE GEOM, only : dlon
      USE RADNCB, only : LM_REQ

      IMPLICIT NONE
      SAVE
C**** Accumulating_period information
      INTEGER, DIMENSION(12) :: MONACC  !@var MONACC(1)=#Januaries, etc
      CHARACTER*12 :: ACC_PERIOD='PARTIAL'    !@var string MONyyr1-yyr2
!@var AMON0,JMON0,JDATE0,JYEAR0,JHOUR0,Itime0  beg.of acc-period

C**** ACCUMULATING DIAGNOSTIC ARRAYS
!@param KAJ number of accumulated zonal budget diagnostics
      INTEGER, PARAMETER :: KAJ=80
!@var AJ zonal budget diagnostics for each surface type
      REAL*8, DIMENSION(JM,KAJ,NTYPE) :: AJ

!@param NREG number of regions for budget diagnostics
      INTEGER, PARAMETER :: NREG=24
!@var AREG regional budget diagnostics
      REAL*8, DIMENSION(NREG,KAJ) :: AREG
!@var TITREG,NAMREG title and names of regions for AREG diagnostics
      CHARACTER*4 TITREG*80,NAMREG(2,23)
!@var JREH lat/lon array defining regions for AREG diagnostics
      INTEGER, DIMENSION(IM,JM) :: JREG

!@param KAPJ number of zonal pressure diagnostics
      INTEGER, PARAMETER :: KAPJ=2
!@var APJ zonal pressure diagnostics
      REAL*8, DIMENSION(JM,KAPJ) :: APJ

!@param KAJL,KAJLX number of AJL diagnostics,KAJLX includes composites
      INTEGER, PARAMETER :: KAJL=70+KEP, KAJLX=KAJL+50
!@var AJL latitude/height diagnostics
      REAL*8, DIMENSION(JM,LM,KAJL) :: AJL

!@param KASJL number of ASJL diagnostics
      INTEGER, PARAMETER :: KASJL=4
!@var ASJL latitude/height supplementary diagnostics (merge with AJL?)
      REAL*8, DIMENSION(JM,LM_REQ,KASJL) :: ASJL

!@param KAIJ,KAIJX number of AIJ diagnostics, KAIJX includes composites
      INTEGER, PARAMETER :: KAIJ=180 , KAIJX=KAIJ+100
!@var AIJ latitude/longitude diagnostics
      REAL*8, DIMENSION(IM,JM,KAIJ) :: AIJ

!@param KAIL number of AIL diagnostics
      INTEGER, PARAMETER :: KAIL=15
!@var AIL longitude/height diagnostics
      REAL*8, DIMENSION(IM,LM,KAIL) :: AIL
!@var J50N,J70N,J5NUV,J5SUV,J5S,J5N special latitudes for AIL diags
      INTEGER, PARAMETER :: J50N  = (50.+90.)*(JM-1)/180.+1.5
      INTEGER, PARAMETER :: J70N  = (70.+90.)*(JM-1)/180.+1.5
      INTEGER, PARAMETER :: J5NUV = (90.+5.)*(JM-1.)/180.+2.
      INTEGER, PARAMETER :: J5SUV = (90.-5.)*(JM-1.)/180.+2.
      INTEGER, PARAMETER :: J5N   = (90.+5.)*(JM-1.)/180.+1.5
      INTEGER, PARAMETER :: J5S   = (90.-5.)*(JM-1.)/180.+1.5

C NEHIST=(TROPO/L STRAT/M STRAT/U STRAT)X(ZKE/EKE/SEKE/ZPE/EPE)X(SH/NH)
!@param NED number of different energy history diagnostics
!@param NEHIST,HIST_DAYS number of energy history columns,rows (max)
      INTEGER, PARAMETER :: NED=10
      INTEGER, PARAMETER :: NEHIST=NED*(2+ISTRAT)
      INTEGER, PARAMETER :: HIST_DAYS=100
!@var ENERGY energy diagnostics
      REAL*8, DIMENSION(NEHIST,HIST_DAYS) :: ENERGY

!@var NPTS number of points at which standard conserv. diags are called
      INTEGER, PARAMETER :: NPTS = 11
!@param NQUANT Number of conserved quantities in conservation diags
      INTEGER, PARAMETER :: NQUANT=22
!@param KCON number of conservation diagnostics
      INTEGER, PARAMETER :: KCON=170
!@var CONSRV conservation diagnostics
      REAL*8, DIMENSION(JM,KCON) :: CONSRV
!@var SCALE_CON scales for conservation diagnostics
      REAL*8, DIMENSION(KCON) :: SCALE_CON
!@var TITLE_CON titles for conservation diagnostics
      CHARACTER*32, DIMENSION(KCON) :: TITLE_CON
!@var NSUM_CON indices for summation of conservation diagnostics
!@var IA_CON IDACC numbers for conservation diagnostics
      INTEGER, DIMENSION(KCON) :: NSUM_CON, IA_CON
!@var NOFM indices for CONSRV array
      INTEGER, DIMENSION(NPTS+1,NQUANT) :: NOFM
!@var icon_xx indexes for conservation quantities
      INTEGER icon_AM,icon_KE,icon_MS,icon_TPE,icon_WM,icon_LKM
     *     ,icon_LKE,icon_EWM,icon_WTG,icon_HTG,icon_OCE,icon_OMSI
     *     ,icon_OHSI,icon_OSSI,icon_LMSI,icon_LHSI,icon_MLI,icon_HLI
!@var KCMX actual number of conservation diagnostics
      INTEGER :: KCMX = 25 ! take up first 25 indexes for special cases
!@var CONPT0 default titles for each point where conserv diags. are done
      CHARACTER*10, DIMENSION(NPTS) :: CONPT0 = (/
     *     "DYNAMICS  ","CONDENSATN","RADIATION ","PRECIPITAT",
     *     "LAND SURFC","SURFACE   ","FILTER    ","OCEAN     ",
     *     "DAILY     ","SRF OCN FL","OCN DYNAM "/)

!@param KSPECA,NSPHER number of spectral diagnostics, and harmonics used
      INTEGER, PARAMETER :: KSPECA=20
      INTEGER, PARAMETER :: NSPHER=4*(2+ISTRAT)
!@var SPECA spectral diagnostics
      REAL*8, DIMENSION((IMH+1),KSPECA,NSPHER) :: SPECA
!@var KLAYER index for dividing up atmosphere into layers for spec.anal.
      INTEGER, DIMENSION(LM) :: KLAYER
!@param PSPEC pressure levels at which layers are seperated and defined
C**** 1000 - 150: troposphere           150 - 10 : low strat.
C****   10 - 1: mid strat               1 and up : upp strat.
      REAL*8, DIMENSION(4), PARAMETER :: PSPEC = (/ 150., 10., 1., 0. /)
!@var LSTR level of interface between low and mid strat. (approx 10 mb)
      INTEGER :: LSTR = LM   ! defaults to model top.

!@param KTPE number of spectral diagnostics for pot. enthalpy
      INTEGER, PARAMETER :: KTPE=8
      integer, parameter :: NHEMI=2
!@var ATPE pot. enthalpy spectral diagnostics
      REAL*8, DIMENSION(KTPE,NHEMI) :: ATPE

!@param HR_IN_DAY hours in day
      INTEGER, PARAMETER :: HR_IN_DAY=24
!@param NDIUVAR number of diurnal diagnostics
      INTEGER, PARAMETER :: NDIUVAR=56
!@param NDIUPT number of points where diurnal diagnostics are kept
      INTEGER, PARAMETER :: NDIUPT=4
!@dbparam IJDD,NAMDD (i,j)-coord.,names of boxes w/diurnal cycle diag
      INTEGER, DIMENSION(2,NDIUPT) :: IJDD
      CHARACTER*4, DIMENSION(NDIUPT) :: NAMDD
      DATA        IJDD    /  63,17,  17,34,  37,27,  13,23 /
      DATA        NAMDD   / 'AUSD', 'MWST', 'SAHL', 'EPAC' /
!@var ADIURN diurnal diagnostics (24 hour cycles at selected points)
      REAL*8, DIMENSION(HR_IN_DAY,NDIUVAR,NDIUPT) :: ADIURN

!@param KAJK number of zonal constant pressure diagnostics
!@param KAJKX number of zonal constant pressure composit diagnostics
      INTEGER, PARAMETER :: KAJK=51, KAJKX=KAJK+100
!@var AJK zonal constant pressure diagnostics
      REAL*8, DIMENSION(JM,LM,KAJK) :: AJK

!@param KAIJK,KAIJX number of lat/lon constant pressure diagnostics
      INTEGER, PARAMETER :: KAIJK=6 , kaijkx=kaijk+100
!@var KAIJK lat/lon constant pressure diagnostics
      REAL*8, DIMENSION(IM,JM,LM,KAIJK) :: AIJK

!@param NWAV_DAG number of components in spectral diagnostics
      INTEGER, PARAMETER :: NWAV_DAG=min(9,imh)
!@param Max12HR_sequ,Min12HR_sequ lengths of time series for wave powers
      INTEGER, PARAMETER :: Max12HR_sequ=2*31, Min12HR_sequ=2*28
!@param RE_AND_IM complex components of wave power diagnostics
      INTEGER, PARAMETER :: RE_AND_IM=2
!@param KWP number of wave power diagnostics
      INTEGER, PARAMETER :: KWP=12
!@var WAVE frequency diagnostics (wave power)
      REAL*8,
     &     DIMENSION(RE_AND_IM,Max12HR_sequ,NWAV_DAG,KWP) :: WAVE

C**** parameters and variables for ISCCP diags
!@param ntau,npress number of ISCCP optical depth,pressure categories
      integer, parameter :: ntau=7,npres=7
!@param nisccp number of ISCCP histogram regions
      integer, parameter :: nisccp = 6
!@var isccp_reg latitudinal index for ISCCP histogram regions
      integer :: isccp_reg(JM)
!@var AISCCP accumlated array of ISCCP histogram
      real*8 :: AISCCP(ntau,npres,nisccp)

!@param KGZ number of pressure levels for some diags
      INTEGER, PARAMETER :: KGZ = 13
!@param kgz_max is the actual number of geopotential heights saved
      INTEGER kgz_max
!@param PMB pressure levels for geopotential heights (extends to strat)
!@param GHT ~mean geopotential heights at PMB level (extends to strat)
!@param PMNAME strings describing PMB pressure levels
      REAL*8, DIMENSION(KGZ), PARAMETER ::
     &     PMB=(/1000d0,850d0,700d0,500d0,300d0,100d0,30d0,10d0,
     *     3.4d0,.7d0,.16d0,.07d0,.03d0/),
     *     GHT=(/0.,1500.,3000.,5600.,9500.,16400.,24000.,30000.,
     *     40000.,50000.,61000.,67000.,72000./)
      CHARACTER*4, DIMENSION(KGZ), PARAMETER :: PMNAME=(/
     *     "1000","850 ","700 ","500 ","300 ","100 ","30  ","10  ",
     *     "3.4 ","0.7 ",".16 ",".07 ",".03 " /)

C**** Instantaneous constant pressure level fields
!@var Z_inst saved instantaneous height field (at PMB levels)
!@var RH_inst saved instantaneous relative hum (at PMB levels)
!@var T_inst saved instantaneous temperature(at PMB levels)
      REAL*8, DIMENSION(KGZ,IM,JM) :: Z_inst,RH_inst,T_inst

!@param KACC total number of diagnostic elements
      INTEGER, PARAMETER :: KACC= JM*KAJ*NTYPE + NREG*KAJ
     *     + JM*KAPJ + JM*LM*KAJL + JM*LM_REQ*KASJL + IM*JM*KAIJ +
     *     IM*LM*KAIL + NEHIST*HIST_DAYS + JM*KCON +
     *     (IMH+1)*KSPECA*NSPHER + KTPE*NHEMI + HR_IN_DAY*NDIUVAR*NDIUPT
     *     + RE_AND_IM*Max12HR_sequ*NWAV_DAG*KWP + JM*LM*KAJK +
     *     IM*JM*LM*KAIJK+ntau*npres*nisccp

      COMMON /ACCUM/ AJ,AREG,APJ,AJL,ASJL,AIJ,AIL,
     &  ENERGY,CONSRV,SPECA,ATPE,ADIURN,WAVE,
     &  AJK,AIJK,AISCCP
      REAL*8, DIMENSION(KACC) :: ACC
      REAL*8, DIMENSION(LM+LM_REQ+1,IM,JM,5) :: AFLX_ST
      EQUIVALENCE (ACC,AJ,AFLX_ST)

!@param KTSF number of freezing temperature diagnostics
      integer, parameter :: ktsf=4
!@var TSFREZ freezing temperature diagnostics
C****   1  FIRST DAY OF GROWING SEASON (JULIAN DAY)
C****   2  LAST DAY OF GROWING SEASON (JULIAN DAY)
C****   3  LAST DAY OF ICE-FREE LAKE (JULIAN DAY)
C****   4  LAST DAY OF ICED-UP LAKE  (JULIAN DAY)
      REAL*8, DIMENSION(IM,JM,KTSF) :: TSFREZ

!@param KTD number of diurnal temperature diagnostics
      INTEGER, PARAMETER :: KTD=9
!@var TDIURN diurnal range temperature diagnostics
C****   1  MIN TG1 OVER EARTH FOR CURRENT DAY (C)
C****   2  MAX TG1 OVER EARTH FOR CURRENT DAY (C)
C****   3  MIN TS OVER EARTH FOR CURRENT DAY (K)
C****   4  MAX TS OVER EARTH FOR CURRENT DAY (K)
C****   5  SUM OF COMPOSITE TS OVER TIME FOR CURRENT DAY (C)
C****   6  MAX COMPOSITE TS FOR CURRENT DAY (K)
C****   7  MAX TG1 OVER OCEAN ICE FOR CURRENT DAY (C)
C****   8  MAX TG1 OVER LAND ICE FOR CURRENT DAY (C)
C****   9  MIN COMPOSITE TS FOR CURRENT DAY (K)
      REAL*8, DIMENSION(IM,JM,KTD) :: TDIURN

!@nlparam KDIAG array of flags to control diagnostics printout
      INTEGER, DIMENSION(12) :: KDIAG

!@param NKEYNR number of key number diagnostics
      INTEGER, PARAMETER :: NKEYNR=42
!@param NKEYMO number of months key diagnostics are saved
      INTEGER, PARAMETER :: NKEYMO=50
!@var KEYNR time-series of key numbers
      INTEGER, DIMENSION(NKEYNR,NKEYMO) :: KEYNR = 0
!@var KEYCT next index in KEYNR to be used (1->nkeymo)
      INTEGER :: KEYCT = 1

!@nlparam IWRITE,JWRITE,ITWRITE control rad.debug output (i,j,amount)
      INTEGER :: IWRITE = 0, JWRITE = 0, ITWRITE = 0
!@nlparam QDIAG TRUE for outputting binary diagnostics
      LOGICAL :: QDIAG = .FALSE.
!@nlparam QDIAG_ratios TRUE for forming ratios if title="q1 x q2"
      LOGICAL :: QDIAG_ratios = .TRUE.

!@var OA generic diagnostic array for ocean heat transport calculations
C****
C****       DATA SAVED IN ORDER TO CALCULATE OCEAN TRANSPORTS
C****
C****       1  ACE1I+SNOWOI  (INSTANTANEOUS AT NOON GMT)
C****       2  MSI2   (INSTANTANEOUS AT NOON GMT)
C****       3  HSIT   (INSTANTANEOUS AT NOON GMT)
C****       4  ENRGP  (INTEGRATED OVER THE DAY)
C****       5  SRHDT  (FOR OCEAN, INTEGRATED OVER THE DAY)
C****       6  TRHDT  (FOR OCEAN, INTEGRATED OVER THE DAY)
C****       7  SHDT   (FOR OCEAN, INTEGRATED OVER THE DAY)
C****       8  EVHDT  (FOR OCEAN, INTEGRATED OVER THE DAY)
C****       9  TRHDT  (FOR OCEAN ICE, INTEGRATED OVER THE DAY)
C****      10  SHDT   (FOR OCEAN ICE, INTEGRATED OVER THE DAY)
C****      11  EVHDT  (FOR OCEAN ICE, INTEGRATED OVER THE DAY)
C****      12  SRHDT  (FOR OCEAN ICE, INTEGRATED OVER THE DAY)
C****
C**** Extra array needed for dealing with advected ice
C****      13  HCHSI  (HORIZ CONV SEA ICE ENRG, INTEGRATED OVER THE DAY)
C****
!@param KOA number of diagnostics needed for ocean heat transp. calcs
      INTEGER, PARAMETER :: KOA = 13  ! 12
      REAL*8, DIMENSION(IM,JM,KOA) :: OA

C****
C**** Information about acc-arrays:
C****      names, indices, units, idacc-numbers, etc.

!@var iparm/dparm int/double global parameters written to acc-file
      integer, parameter :: niparm_max=100
      character(len=20), dimension(niparm_max) :: iparm_name
      integer, dimension(niparm_max) :: iparm
      integer :: niparm=0
      integer, parameter :: ndparm_max=100
      character(len=20), dimension(ndparm_max) :: dparm_name
      REAL*8, dimension(ndparm_max) :: dparm
      integer :: ndparm=0

!@var J_xxx zonal J diagnostic names
      INTEGER :: J_SRINCP0, J_SRNFP0, J_SRNFP1, J_SRABS, J_SRINCG,
     *     J_SRNFG, J_TRNFP0, J_TRNFP1, J_TRHDT, J_RNFP0, J_RNFP1,
     *     J_RHDT, J_SHDT, J_EVHDT, J_HZ1, J_TG2, J_TG1, J_EVAP,
     *     J_PRCP, J_TX, J_TX1, J_TSRF, J_DTSGST, J_DTDGTR, J_RICST,
     *     J_RICTR, J_ROSST, J_ROSTR, J_RSI, J_TYPE, J_RSNOW,
     *     J_OHT, J_DTDJS, J_DTDJT, J_LSTR, J_LTRO, J_EPRCP,
     *     J_RUN, J_ERUN, J_HZ0, J_H2OCH4,
     *     J_RVRD,J_ERVR,J_IMELT, J_HMELT, J_SMELT,J_IMPLM, J_IMPLH,
     *     J_WTR1,J_ACE1, J_WTR2,J_ACE2, J_SNOW, J_BRTEMP, J_HZ2,
     *     J_PCLDSS,J_PCLDMC, J_PCLD,J_CTOPP, J_PRCPSS, J_PRCPMC, J_QP,
     *     J_GAM,J_GAMM, J_GAMC,J_TRINCG, J_FTHERM, J_HSURF, J_HATM,
     *     J_PLAVIS,J_PLANIR,J_ALBVIS, J_ALBNIR, J_SRRVIS, J_SRRNIR,
     *     J_SRAVIS,J_SRANIR,J_CLDDEP, J_CLRTOA, J_CLRTRP, J_TOTTRP
!@var NAME_J,UNITS_J Names/Units of zonal J diagnostics
      character(len=20), dimension(kaj) :: name_j,units_j
!@var LNAME_J Long names of zonal J diagnostics
      character(len=80), dimension(kaj) :: lname_j
!@var STITLE_J short titles for print out for zonal J diagnostics
      character(len=16), dimension(kaj) :: stitle_j
!@var SCALE_J scale for zonal J diagnostics
      real*8, dimension(kaj) :: scale_j
!@var IA_J IDACC indexes for zonal J diagnostics
      integer, dimension(kaj) :: ia_j
!@var k_j_out number of directly printed out budget diags
      integer k_j_out

      character(len=20), dimension(kaj) :: name_reg
      character(len=20), dimension(kapj) :: name_pj,units_pj
      character(len=80), dimension(kapj) :: lname_pj

!@var IJ_xxx AIJ diagnostic names
      INTEGER :: IJ_RSOI, IJ_RSNW, IJ_SNOW, IJ_SHDT, IJ_PREC, IJ_EVAP,
     *     IJ_SSAT, IJ_BETA,  IJ_SLP1,  IJ_P4UV, IJ_PRES, IJ_PHI1K,
     *     IJ_PHI850, IJ_PHI700, IJ_PHI500, IJ_PHI300, IJ_PHI100,
     *     IJ_PHI30, IJ_PHI10, IJ_PHI3p4, IJ_PHI0p7, IJ_PHI0p16,
     *     IJ_PHI0p07, IJ_PHI0p03, IJ_T850, IJ_T500, IJ_T300, IJ_Q850,
     *     IJ_Q500, IJ_Q300, IJ_PMCCLD, IJ_CLDTPPR, IJ_CLDCV, IJ_DSEV,
     *     IJ_CLDTPT, IJ_CLDCV1, IJ_CLDT1T,IJ_CLDT1P,
     *     ij_wtrcld,ij_icecld,ij_optdw,ij_optdi,
     *     IJ_RH1, IJ_RH850, IJ_RH500, IJ_RH300,
     *     IJ_TRNFP0, IJ_SRTR, IJ_NETH, IJ_SRNFP0, IJ_SRINCP0, IJ_SRNFG,
     *     IJ_SRINCG, IJ_TG1, IJ_RSIT, IJ_TDSL, IJ_TDCOMP, IJ_DTDP,
     *     IJ_RUNE, IJ_TS1, IJ_RUNLI, IJ_WS, IJ_TS, IJ_US, IJ_VS,
     *     IJ_SLP, IJ_UJET, IJ_VJET, IJ_PCLDL, IJ_PCLDM, IJ_PCLDH,
     *     IJ_BTMPW, IJ_SRREF, IJ_SRVIS, IJ_TOC2, IJ_TAUS, IJ_TAUUS,
     *     IJ_TAUVS, IJ_GWTR, IJ_QS, IJ_STRNGTS, IJ_ARUNU, IJ_DTGDTS,
     *     IJ_PUQ, IJ_PVQ, IJ_TGO, IJ_MSI2, IJ_TGO2, IJ_EVAPO,
     *     IJ_EVAPI, IJ_EVAPLI,IJ_EVAPE, IJ_F0OC,IJ_F0OI,IJ_F0LI,IJ_F0E,
     *     IJ_F1LI, IJ_SNWF, IJ_TSLI, IJ_ERUN2, IJ_SHDTLI, IJ_EVHDT,
     *     IJ_TRHDT, IJ_TMAX, IJ_TMIN, IJ_TMNMX, IJ_PEVAP, IJ_TMAXE,
     *     IJ_WMSUM, IJ_PSCLD, IJ_PDCLD, IJ_DCNVFRQ, IJ_SCNVFRQ,
     *     IJ_EMTMOM, IJ_SMTMOM, IJ_FMU, IJ_FMV, IJ_SSTABX,
     *     IJ_FGZU, IJ_FGZV, IJ_ERVR, IJ_MRVR,
     *     IJ_LKON, IJ_LKOFF, IJ_LKICE, IJ_PTROP, IJ_TTROP, IJ_TSI,
     *     IJ_SSI1,IJ_SSI2,IJ_SMFX, IJ_MSU2, IJ_MSU2R, IJ_MSU3, IJ_MSU4,
     *     IJ_MLTP,IJ_FRMP, IJ_P850, IJ_CLR_SRINCG
!@var IJ_Gxx names for old AIJG arrays (should be more specific!)
      INTEGER :: IJ_G01,IJ_G02,IJ_G03,IJ_G04,IJ_G05,IJ_G06,IJ_G07,
     *     IJ_G08,IJ_G09,IJ_G10,IJ_G11,IJ_G12,IJ_G13,IJ_G14,IJ_G15,
     *     IJ_G16,IJ_G17,IJ_G18,IJ_G19,IJ_G20,IJ_G21,IJ_G22,IJ_G23,
     *     IJ_G24,IJ_G25,IJ_G26,IJ_G27,IJ_G28,IJ_G29
!@var IJ_GWx names for gravity wave diagnostics
      INTEGER :: IJ_GW1,IJ_GW2,IJ_GW3,IJ_GW4,IJ_GW5,IJ_GW6,IJ_GW7,IJ_GW8
     *     ,IJ_GW9
!@var IJ_xxxI names for ISCCP diagnostics
      INTEGER :: IJ_CTPI,IJ_TAUI,IJ_LCLDI,IJ_MCLDI,IJ_HCLDI,IJ_TCLDI

!@param LEGEND "contour levels" for ij-maps
      CHARACTER(LEN=40), DIMENSION(25), PARAMETER :: LEGEND=(/ !
     1  '0=0,1=5...9=45,A=50...K=100             ', ! ir_pct    fac=.2
     2  '0=0...9=90,A=100...I=180...R=270        ', ! ir_angl       .1
     3  '1=.5...9=4.5,A=5...Z=17.5,+=MORE        ', ! ir_0_18        2
     4  '1=.1...9=.9,A=1...Z=3.5,+=MORE          ', ! ir_0_4        10
     5  '1=2...9=18,A=20...Z=70,+=MORE           ', ! ir_0_71       .5
     6  '1=50...9=450,A=500...Z=1750,+=MORE      ', ! ir_0_1775     .02
     7  '1=100...9=900,A=1000...Z=3500,+=MORE    ', ! ir_0_3550     .01
     8  '1=20...9=180,A=200...Z=700,+=MORE       ', ! ir_0_710      .05
     9  'A=1...Z=26,3=30...9=90,+=100-150,*=MORE ', ! ir_0_26_150    1
     O  '0=0,A=.1...Z=2.6,3=3...9=9,+=10-15      ', ! ir_0_3_15     10
     1  '-=LESS,Z=-78...0=0...9=27,+=MORE        ', ! ir_m80_28     .33
     2  '-=LESS,Z=-260...0=0...9=90,+=MORE       ', ! ir_m265_95    .1
     3  '-=LESS,Z=-520...0=0...9=180,+=MORE      ', ! ir_m530_190   .05
     4  '-=LESS,Z=-1300...0=0...9=450,+=MORE     ', ! ir_m1325_475  .02
     5  '-=LESS,Z=-2600...0=0...9=900,+=MORE     ', ! ir_m2650_950  .01
     6  '-=LESS,Z=-3900...0=0...9=1350,+=MORE    ', ! ir_m3975_1425 .007
     7  '-=LESS,Z=-5200...0=0...9=1800,+=MORE    ', ! ir_m5300_1900 .005
     8  '-=LESS,9=-.9...0=0,A=.1...Z=2.6,+=MORE  ', ! ir_m1_3       10
     9  '-=LESS,9=-45...0=0,A=5...I=45...+=MORE  ', ! ir_m45_130    .2
     O  '-=LESS,9=-90...0=0,A=10...Z=260,+=MORE  ', ! ir_m95_265    .1
     1  '-=LESS,9=-180...A=20...Z=520,+=MORE     ', ! ir_m190_530   .05
     2  '-=LESS,9=-9...0=0,A=1...Z=26,+=MORE     ', ! ir_m9_26       1
     3  '-=LESS,9=-36...0=0,A=4...Z=104,+=MORE   ', ! ir_m38_106    .25
     4  '1=5...9=45,A=50...Z=175,+=MORE          ', ! ir_0_180      .2
     5  '9=-512...1=-2,0=0,A=2,B=4,C=8...+=MORE  '/)! ir_log2       1.
!@var ir_xxxx names for indices to LEGEND indicating the (rounded) range
      integer, parameter :: ir_pct=1, ir_angl=2, ir_0_18=3, ir_0_4=4,
     * ir_0_71=5, ir_0_1775=6, ir_0_3550=7, ir_0_710=8, ir_0_26_150=9,
     * ir_0_3_15=10, ir_m80_28=11, ir_m265_95=12, ir_m530_190=13,
     * ir_m1325_475=14, ir_m2650_950=15, ir_m3975_1425=16,
     * ir_m5300_1900=17, ir_m1_3=18, ir_m45_130=19, ir_m95_265=20,
     * ir_m190_530=21, ir_m9_26=22, ir_m38_106=23, ir_0_180=24,
     * ir_log2=25
!@var fac_legnd = 1/(range_of_1_colorbox)
      real*8, dimension(25) :: fac_legnd=(/
     1      1d0/5,  1d0/10,    2.d0,   10.d0,   1d0/2,
     6     1d0/50, 1d0/100,  1d0/20,    1.d0,   10.d0,
     1      1d0/3,  1d0/10,  1d0/20,  1d0/50, 1d0/100,
     6    1d0/150, 1d0/200,   10.d0,   1d0/5,  1d0/10,
     1     1d0/20,    1.d0,   1d0/4,   1d0/5,     1d0  /)

!@param CBAR "color bars" for ij-maps
      CHARACTER(LEN=38), PARAMETER, DIMENSION(5) :: CBAR=(/
     &     ' 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ+',  ! ib_pos
     &     ' 0123456789ABCDEFGHIJKX               ',  ! ib_pct
     &     '-9876543210ABCDEFGHIJKLMNOPQRSTUVWXYZ+',  ! ib_npp,ib_ntr
     &     ' 0ABCDEFGHIJKLMNOPQRSTUVWXYZ3456789+* ',  ! ib_hyb
     &     '-ZYXWVUTSRQPONMLKJIHGFEDCBA0123456789+'/) ! ib_nnp
!@var ib_xxx indices for color bars
      integer, parameter :: ib_pos=1,ib_pct=2,ib_npp=3,ib_hyb=4,ib_nnp=5
     &     ,ib_ntr=6

!@dbparam isccp_diags: if 1 accumulate ISCCP cloud data (default 0)
      INTEGER :: isccp_diags = 0

!@var SCALE_IJ scaling for weighted AIJ diagnostics
      REAL*8, DIMENSION(KAIJ) :: SCALE_IJ
!@var NAME_IJ,UNITS_IJ Names/Units of lat/lon IJ diagnostics
      character(len=30), dimension(kaijx) :: name_ij,units_ij
!@var LNAME_IJ Long names of lat/lon IJ diagnostics
      character(len=80), dimension(kaijx) :: lname_ij
!@var IW_IJ weighting indices for IJ diagnostics
      integer, dimension(kaij) :: iw_ij
!@var nwts_ij = number of weight-ij-arrays used in IJ-diagnostics
      integer, parameter :: nwts_ij = 7
!@var wt_ij various weight-arrays use in ij-diagnostics
      real*8, dimension(im,jm,nwts_ij) :: wt_ij
!@var IW_xxx index for weight-array
      integer, parameter :: iw_all=1 , iw_ocn=2 , iw_lake=3,
     *   iw_lice=4 , iw_soil=5 , iw_bare=6 , iw_veg=7
!@var IR_IJ range indices for IJ diagnostics
      integer, dimension(kaij) :: ir_ij
!@var IA_IJ IDACC indexes for lat/lon IJ diagnostics
      integer, dimension(kaij) :: ia_ij
!@var jgrid_ij 1=primary grid  2=secondary grid
      integer, dimension(kaij) :: jgrid_ij

!@var JL_xxx names for JL diagnostic indices
      INTEGER ::
     &     jl_mcmflx,jl_srhr,jl_trcr,jl_sshr,jl_trbhr,jl_mchr,jl_totntlh
     *     ,jl_zmfntlh,jl_totvtlh,jl_zmfvtlh,jl_ape,jl_dtdyn,jl_dudfmdrg
     *     ,jl_totcld,jl_dumtndrg,jl_dushrdrg,jl_dumcdrgm10
     *     ,jl_dumcdrgp10,jl_dumcdrgm40,jl_dumcdrgp40,jl_dumcdrgm20
     *     ,jl_dumcdrgp20,jl_sscld,jl_mccld,jl_sdifcoef,jl_dudtsdif
     *     ,jl_gwfirst,jl_dtdtsdrg,jl_epflxv,jl_rhe,jl_epflxn,jl_damdc
     *     ,jl_dammc,jl_40,jl_uepac,jl_vepac,jl_wepac,jl_uwpac,jl_vwpac
     *     ,jl_wwpac,jl_47,jl_zmfntmom,jl_totntmom,jl_mchphas,jl_mcdtotw
     *     ,jl_dudtsdrg,jl_mcldht,jl_trbke,jl_trbdlht,jl_mcheat,jl_mcdry
     *     ,jl_cldmc,jl_cldss,jl_csizmc,jl_csizss
     *     ,jl_wcld,jl_icld,jl_wcod,jl_icod,jl_wcsiz,jl_icsiz

!@var SNAME_JL Names of lat-sigma JL diagnostics
      character(len=30), dimension(kajlx) :: sname_jl
!@var LNAME_JL,UNITS_JL Descriptions/Units of JL diagnostics
      character(len=50), dimension(kajlx) :: lname_jl,units_jl
!@var SCALE_JL printout scaling factors for JL diagnostics
      REAL*8, dimension(kajlx) :: scale_jl
!@var IA_JL,JGRID_JL idacc-numbers,gridtypes for JL diagnostics
      integer, dimension(kajlx) :: ia_jl,jgrid_jl
!@var POW_JL printed output scaled by 10**(-pow_jl)
      integer, dimension(kajlx) :: pow_jl

!@var NAME_SJL Names of radiative-layer-only SJL diagnostics
      character(len=30), dimension(kasjl) :: name_sjl
!@var LNAME_SJL,UNITS_SJL Descriptions/Units of SJL diagnostics
      character(len=50), dimension(kasjl) :: lname_sjl,units_sjl
!@var SCALE_SJL printout scaling factors for SJL diagnostics
      REAL*8, dimension(kasjl) :: scale_sjl
!@var IA_SJL idacc-numbers for SJL diagnostics
      integer, dimension(kasjl) :: ia_sjl

!@var JK_xxx names for JK diagnostic indices
      INTEGER ::
     &     JK_dpa ,JK_dpb ,JK_temp ,JK_hght
     &    ,JK_q ,JK_theta ,JK_rh ,JK_u
     &    ,JK_v ,JK_zmfke ,JK_totke ,JK_zmfntsh
     &    ,JK_totntsh ,JK_zmfntgeo ,JK_totntgeo ,JK_zmfntlh
     &    ,JK_totntlh ,JK_zmfntke ,JK_totntke ,JK_zmfntmom
     &    ,JK_totntmom ,JK_p2kedpgf ,JK_dpsqr ,JK_nptsavg
     &    ,JK_vvel ,JK_zmfvtdse ,JK_totvtdse ,JK_zmfvtlh
     &    ,JK_totvtlh ,JK_vtgeoeddy ,JK_barekegen ,JK_potvort
     &    ,JK_vtpv ,JK_vtpveddy ,JK_nptsavg1 ,JK_totvtke
     &    ,JK_vtameddy ,JK_totvtam ,JK_sheth ,JK_dudtmadv
     &    ,JK_dtdtmadv ,JK_dudttem ,JK_dtdttem ,JK_epflxncp
     &    ,JK_epflxvcp ,JK_uinst ,JK_totdudt ,JK_tinst
     &    ,JK_totdtdt ,JK_eddvtpt ,JK_cldh2o

!@var SNAME_JK Names of lat-pressure JK diagnostics
      character(len=30), dimension(kajkx) :: sname_jk
!@var LNAME_JK,UNITS_JK Descriptions/Units of JK diagnostics
      character(len=50), dimension(kajkx) :: lname_jk,units_jk
!@var SCALE_JK printout scaling factors for JK diagnostics
      REAL*8, dimension(kajkx) :: scale_jk
!@var IA_JK,JGRID_JK idacc-numbers,gridtypes for JK diagnostics
      integer, dimension(kajkx) :: ia_jk,jgrid_jk
!@var POW_JK printed output scaled by 10**(-pow_jk)
      integer, dimension(kajkx) :: pow_jk

!@var IJK_xxx AIJK diagnostic names
      INTEGER :: IJK_U, IJK_V, IJK_DSE, IJK_DP, IJK_T, IJK_Q
!@var SCALE_IJK scaling for weighted AIJK diagnostics
      REAL*8, DIMENSION(KAIJKx) :: SCALE_IJK
!@var OFF_IJK offset for weighted AIJK diagnostics
      REAL*8, DIMENSION(KAIJKx) :: OFF_IJK

!@var NAME_IJK Names of lon-lat-pressure IJK diagnostics
      character(len=30), dimension(kaijkx) :: name_ijk
!@var LNAME_IJK,UNITS_IJK Descriptions/Units of IJK diagnostics
      character(len=50), dimension(kaijkx) :: lname_ijk,units_ijk
!@var jgrid_ijk 1=primary grid  2=secondary grid
      integer, dimension(KAIJKx) :: jgrid_ijk

      character(len=20), dimension(kwp) :: name_wave,units_wave
      character(len=80), dimension(kwp) :: lname_wave

      character(len=20), dimension(kcon) :: name_consrv,units_consrv
      character(len=80), dimension(kcon) :: lname_consrv

      character(len=20), dimension(kail) :: name_il,units_il
      character(len=80), dimension(kail) :: lname_il
      real*8, dimension(kail) :: scale_il
      integer, dimension(kail) :: ia_il
!@var IL_xxx names for longitude height diagnostics
      INTEGER :: IL_UEQ,IL_VEQ,IL_WEQ,IL_TEQ,IL_QEQ,IL_MCEQ,IL_REQ
     *     ,IL_W50N,IL_T50N,IL_R50N,IL_U50N,IL_W70N,IL_T70N,IL_R70N
     *     ,IL_U70N

      character(len=20), dimension(ndiuvar) :: name_dd,units_dd
      character(len=80), dimension(ndiuvar) :: lname_dd
      real*8, dimension(ndiuvar) :: scale_dd

!@var IDD_xxx names for diurnal diagnostics
      INTEGER :: IDD_ISW, IDD_PALB, IDD_GALB, IDD_ABSA, IDD_ECND,
     *     IDD_SPR, IDD_PT5, IDD_PT4, IDD_PT3, IDD_PT2, IDD_PT1, IDD_TS,
     *     IDD_TG1, IDD_Q5, IDD_Q4, IDD_Q3, IDD_Q2, IDD_Q1, IDD_QS,
     *     IDD_QG, IDD_SWG, IDD_LWG, IDD_SH, IDD_LH, IDD_HZ0, IDD_UG,
     *     IDD_VG, IDD_WG, IDD_US, IDD_VS, IDD_WS, IDD_CIA, IDD_RIS,
     *     IDD_RIG, IDD_CM, IDD_CH, IDD_CQ, IDD_EDS, IDD_DBL, IDD_DCF,
     *     IDD_LDC, IDD_PR, IDD_EV, IDD_DMC, IDD_SMC, IDD_CL7, IDD_CL6,
     *     IDD_CL5, IDD_CL4, IDD_CL3, IDD_CL2, IDD_CL1, IDD_W, IDD_CCV,
     *     IDD_SSP, IDD_MCP

!@var tf_xxx tsfrez diagnostic names
      INTEGER :: tf_day1,tf_last,tf_lkon,tf_lkoff
      character(len=20), dimension(ktsf) :: name_tsf,units_tsf
      character(len=80), dimension(ktsf) :: lname_tsf

      character(len=8), dimension(ntype) :: stype_names=
     &     (/ 'OCEAN   ','OCEANICE','EARTH   ',
     &        'LANDICE ','LAKE    ','LAKEICE ' /)

c idacc-indices of various processes
      integer, parameter ::
     &     ia_src=1, ia_rad=2, ia_srf=3, ia_dga=4, ia_d4a=5, ia_d5f=6,
     *     ia_d5d=7, ia_d5s=8, ia_12hr=9, ia_filt=10, ia_ocn=11,
     *     ia_inst=12

!@var PLE,PLM, PLE_DN ref pressures at upper, middle and lower edge
      REAL*8, DIMENSION(LM) :: PLE
      REAL*8, DIMENSION(LM) :: PLE_DN
      REAL*8, DIMENSION(LM+LM_REQ) :: PLM
!@var P1000K scaling to change reference pressure from 1mb to 1000mb
      REAL*8 :: P1000K
!@var inci,incj print increments for i and j, so maps/tables fit on page
      integer, parameter :: inci=(im+35)/36,incj=(JM+23)/24, jmby2=jm/2
!@var linect = current line on page of print out
      integer linect

!@var XWON scale factor for diag. printout needed for Wonderland model
      REAL*8 :: XWON = TWOPI/(DLON*FIM)

!@var LMOMAX max no. of layers in any ocean
      INTEGER, PARAMETER :: LMOMAX=50
!@var ZOC, ZOC1 ocean depths for diagnostics (m) (ONLY FOR DEEP OCEAN)
      REAL*8 :: ZOC(LMOMAX) = 0. , ZOC1(LMOMAX+1) = 0.

      END MODULE DAGCOM

      SUBROUTINE io_diags(kunit,it,iaction,ioerr)
!@sum  io_diag reads and writes diagnostics to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,ioread_single,irerun
     *    ,iowrite,iowrite_mon,iowrite_single,lhead, idacc,nsampl
     *    ,Kradia
      USE DAGCOM
      IMPLICIT NONE
      REAL*4, save :: ACCS(KACC)
      REAL*4 TSFREZS(IM,JM,KTSF),AFLXS(LM+LM_REQ+1,IM,JM,5)
      integer monac1(12),i_ida,i_xtra
!@var Kcomb counts acc-files as they are added up
      INTEGER, SAVE :: Kcomb=0

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER idac1(12)
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "DIAG01"
!@var it input/ouput value of hour
      INTEGER, INTENT(INOUT) :: it

      if(kradia.gt.0) then
        write (MODULE_HEADER(LHEAD+1:80),'(a6,i8,a20,i3,a7)')
     *   '#acc(=',idacc(2),') R8:SU.SD.TU.TD.dT(',lm+lm_req+1,',ijM,5)'

        SELECT CASE (IACTION)
        CASE (IOWRITE)            ! output to standard restart file
          WRITE (kunit,err=10) MODULE_HEADER,idacc(2),AFLX_ST,it
        CASE (IOWRITE_SINGLE)     ! output in single precision
          MODULE_HEADER(LHEAD+18:LHEAD+18) = '4'
          MODULE_HEADER(LHEAD+44:80) = ',monacc(12)'
          WRITE (kunit,err=10) MODULE_HEADER,idacc(2),
     *          REAL(AFLX_ST,KIND=4), monacc,it
        CASE (IOWRITE_MON)        ! output to end-of-month restart file
          MODULE_HEADER(LHEAD+1:80) = 'itime '
          WRITE (kunit,err=10) MODULE_HEADER,it
        CASE (ioread)           ! input from restart file
          READ (kunit,err=10) HEADER,idacc(2),AFLX_ST,it
          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
            GO TO 10
          END IF
        CASE (IOREAD_SINGLE)      !
          READ (kunit,err=10) HEADER,idac1(2),AFLXS,monac1
          AFLX_ST=AFLX_ST+AFLXS
          IDACC(2) = IDACC(2) + IDAC1(2)
          monacc = monacc + monac1
        END SELECT
        return
      end if

C**** The regular model (Kradia le 0)
      write (MODULE_HEADER(LHEAD+1:LHEAD+15),'(a10,i4,a1)')
     *   'I/R8 keys(',1+NKEYNR*NKEYMO,')'             ! keyct,keynr(:,:)
      i_ida = Lhead + 10+4+1 + 10+2+1 + 1
      write (MODULE_HEADER(LHEAD+10+4+1+1:i_ida-1),'(a10,i2,a1)')
     *   ',TSFR(IJM,',KTSF,')'
      write (MODULE_HEADER(i_ida:i_ida+9),'(a7,i2,a1)')
     *   ',idacc(',nsampl,')'
      write (MODULE_HEADER(i_ida+9+1:i_ida+9 + 5+8+1),'(a5,i8,a1)')
     *   ',acc(',kacc,')'
      i_xtra = i_ida+9 + 5+8+1 + 1

      SELECT CASE (IACTION)
      CASE (IOWRITE)            ! output to standard restart file
        write (MODULE_HEADER(i_xtra:80),             '(a7,i2,a)')
     *   ',x(IJM,',KTD+KOA,')'  ! make sure that i_xtra+7+2 < 80
        WRITE (kunit,err=10) MODULE_HEADER,keyct,KEYNR,TSFREZ,
     *     idacc,ACC,
     *     TDIURN,OA,it
      CASE (IOWRITE_SINGLE)     ! output in single precision
        MODULE_HEADER(LHEAD+1:LHEAD+4) = 'I/R4'
        MODULE_HEADER(i_xtra:80) = ',monacc(12)'
        WRITE (kunit,err=10) MODULE_HEADER,keyct,KEYNR,
     *     REAL(TSFREZ,KIND=4),idacc,REAL(ACC,KIND=4),
     *     monacc,it
      CASE (IOWRITE_MON)        ! output to end-of-month restart file
        MODULE_HEADER(i_ida:80) = ',it '
        WRITE (kunit,err=10) MODULE_HEADER,keyct,KEYNR,TSFREZ,it
      CASE (ioread)           ! input from restart file
        READ (kunit,err=10) HEADER,keyct,KEYNR,TSFREZ,
     *      idacc, ACC,
     *      TDIURN,OA,it
        IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
          PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      CASE (IOREAD_SINGLE)      !
        READ (kunit,err=10) HEADER,keyct,KEYNR,TSFREZS,
     *      idac1,ACCS,
     *      monac1
!**** Here we could check the dimensions written into HEADER  ??????
        TSFREZ=TSFREZS
        ACC=ACC+ACCS
        IDACC = IDACC + IDAC1
!@var idacc(5) is the length of a time series (daily energy history).
!****   If combining acc-files, rather than concatenating these series,
!****   we average their beginnings (up to the length of the shortest)
        Kcomb = Kcomb + 1          ! reverse addition, take min instead
        if (Kcomb.gt.1) IDACC(5) = MIN(IDACC(5)-IDAC1(5),IDAC1(5))
        monacc = monacc + monac1
      CASE (irerun)      ! only keynr,tsfrez needed at beg of acc-period
        READ (kunit,err=10) HEADER,keyct,KEYNR,TSFREZ  ! 'it' not read
        IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
          PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_diags

      SUBROUTINE aPERIOD (JMON1,JYR1,months,years,moff,  aDATE,LDATE)
!@sum  aPERIOD finds a 7 or 12-character name for an accumulation period
!@+   if the earliest month is NOT the beginning of the 2-6 month period
!@+   the name will reflect that fact ONLY for 2 or 3-month periods
!@auth Reto A. Ruedy
!@ver  1.0
      USE MODEL_COM, only : AMONTH
      implicit none
!@var JMON1,JYR1 month,year of beginning of period 1
      INTEGER JMON1,JYR1
!@var JMONM,JMONL middle,last month of period
      INTEGER JMONM,JMONL
!@var months,years length of 1 period,number of periods
      INTEGER months,years
!@var moff = # of months from beginning of period to JMON1 if months<12
      integer moff
!@var yr1,yr2 (end)year of 1st and last period
      INTEGER yr1,yr2
!@var aDATE date string: MONyyr1(-yyr2)
      character*12 aDATE
!@var LDATE length of date string (7 or 12)
      INTEGER LDATE

      LDATE = 7                  ! if years=1
      if(years.gt.1) LDATE = 12

      aDATE(1:12)=' '
      aDATE(1:3)=AMONTH(JMON1)        ! letters 1-3 of month IF months=1
      yr1=JYR1
      JMONL=JMON1+months-1
      if(JMONL.GT.12) then
         yr1=yr1+1
         JMONL=JMONL-12
      end if
      if (moff.gt.0.and.months.le.3) then  ! earliest month is NOT month
        JMONL = 1 + mod(10+jmon1,12)       ! 1 of the 2-3 month period
        yr1=JYR1
        if (jmon1.gt.1) yr1=yr1+1
      end if
      yr2=yr1+years-1
      write(aDATE(4:7),'(i4.4)') yr1
      if(years.gt.1) write(aDATE(8:12),'(a1,i4.4)') '-',yr2

      if(months.gt.12) aDATE(1:1)='x'                ! should not happen
      if(months.le.1 .or. months.gt.12) return

!**** 1<months<13: adjust characters 1-3 of aDATE (=beg) if necessary:
!**** beg=F?L where F/L=letter 1 of First/Last month for 2-11 mo.periods
!****    =F+L                                        for 2 month periods
!****    =FML where M=letter 1 of Middle month       for 3 month periods
!****    =FnL where n=length of period if n>3         4-11 month periods
      aDATE(3:3)=AMONTH(JMONL)(1:1)            ! we know: months>1
      IF (months.eq.2) then
        aDATE(2:2)='+'
        return
      end if
      if (months.eq.3) then
        JMONM = JMONL-1
        if (moff.eq.1) jmonm = jmon1+1
        if (jmonm.gt.12) jmonm = jmonm-12
        if (jmonm.le.0 ) jmonm = jmonm+12
        aDATE(2:2)=AMONTH(JMONM)(1:1)
        return
      end if
      if (moff.gt.0) then  ! can't tell non-consec. from consec. periods
        jmon1 = jmon1-moff
        if (jmon1.le.0) jmon1 = jmon1+12
        JMONL=JMON1+months-1
        if (jmonl.gt.12) jmonl = jmonl-12
        aDATE(1:1)=AMONTH(JMON1)(1:1)
        aDATE(3:3)=AMONTH(JMONL)(1:1)
      end if
      IF (months.ge.4.and.months.le.9) write (aDATE(2:2),'(I1)') months
      IF (months.eq.10) aDATE(2:2)='X'         ! roman 10
      IF (months.eq.11) aDATE(2:2)='B'         ! hex   11
      IF (months.eq.6) THEN                    !    exceptions:
         IF (JMON1.eq. 5) aDATE(1:3)='NHW'     ! NH warm season May-Oct
         IF (JMON1.eq.11) aDATE(1:3)='NHC'     ! NH cold season Nov-Apr
      END IF
      IF (months.eq.7) THEN                    !    to avoid ambiguity:
         IF (JMON1.eq. 1) aDATE(1:3)='J7L'     ! Jan-Jul J7J->J7L
         IF (JMON1.eq. 7) aDATE(1:3)='L7J'     ! Jul-Jan J7J->L7J
      END IF
      IF (months.eq.12) THEN
C****    beg=ANn where the period ends with month n if n<10 (except 4)
         aDATE(1:3)='ANN'                      ! regular annual mean
         IF (JMONL.le. 9) WRITE(aDATE(3:3),'(I1)') JMONL
         IF (JMONL.eq. 4) aDATE(1:3)='W+C'     ! NH warm+cold seasons
         IF (JMONL.eq.10) aDATE(1:3)='C+W'     ! NH cold+warm seasons
         IF (JMONL.eq.11) aDATE(1:3)='ANM'     ! meteor. annual mean
      END IF
      return
      end SUBROUTINE aPERIOD
