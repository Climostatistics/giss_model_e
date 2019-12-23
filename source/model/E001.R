E001.R GISS Model E  2002 modelE                 rar  6/20/02

modelE with 12 lyrs, top at 10 mb - 1979 atmosphere/ocean
no gravity wave drag;     uses dry convection (rather than turbulence)
Sdrag: weak linear strat. drag in top layer
       lost ang.mom is added in below 150 mb
sealevel pressure filter applied every hour, U-filter used
6-band oice albedo; Hogstrom(1984) pbl drag
Note: Many of these choices may be changed using the PARAMETERs below.

Preprocessor Options
!#define TRACERS_ON                  ! include tracers code
End Preprocessor Options

Object modules: (in order of decreasing priority)
RES_M12                             ! horiz/vert resolution, 4x5deg, 12 layers -> 10mb
MODEL_COM GEOM_B IORSF              ! model variables and geometry
MODELE                              ! Main and model overhead
PARAM PARSER                        ! parameter database
ATMDYN_COM ATMDYN MOMEN2ND          ! atmospheric dynamics
QUS_COM QUSDEF QUS_DRV              ! advection of tracers
TQUS_DRV                            ! advection of Q
CLOUDS CLOUDS_DRV CLOUDS_COM        ! clouds modules
SURFACE FLUXES                      ! surface calculation and fluxes
GHY_COM GHY_DRV GHY                 ! land surface and soils
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
! pick exactly one of the next 2 choices: ATURB or DRYCNV
! ATURB                             ! turbulence in whole atmosphere
DRYCNV                              ! drycnv
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_DRV                 ! land ice modules
ICEDYN_DRV ICEDYN  ! or: ICEDYN_DUM ! dynamic sea ice modules
OCEAN OCNML                         ! ocean modules
SNOW_DRV SNOW                       ! snow model
RAD_COM RAD_DRV RADIATION           ! radiation modules
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
CONST FFT72 UTILDBL SYSTEM          ! utilities
POUT                                ! post-processing output

Data input files:
    ! the first group of files is specific to prescribed ocean runs
! AIC=1DEC1951.rsfE000   ! or:    ! initial conditions (atm./ground), no GIC, ISTART=8
AIC=AIC.RES_M12.D771201           ! initial conditions (atm.),     needs GIC, ISTART=2
GIC=GIC.rsfB357M12.1DEC1956.1.ext ! initial conditions (ground)
OSST=OST4X5.B.1975-84avg.Hadl1.1  ! prescr. climatological ocean (1 yr of data)
SICE=SICE4X5.B.1975-84avg.Hadl1.1 ! prescr. climatological sea ice
    ! if the prescr. ocean varies from year to year use instead:
! OSST=OST4X5.B.1950.M02.Hadl1.1  ! ocean data   Feb 1950 - 1999
! SICE=SICE4X5.B.1950.M02.Hadl1.1 ! ocean data   Feb 1950 - 1999
    ! the next 3 files are specific to q-flux ocean runs
! AIC=E001/1JAN1956.rsfE001.O250D      ! AIC/OHT made by aux/mkOTSPEC.E001.M250D
! OHT=E001/OTSPEC.E001.M250D.1951-1955 ! horizontal ocean heat transport
OCNML=Z1O.B4X5.cor                ! mixed layer depth (use for post processing)
    ! files needed for all models
CDN=CD4X500S    ! surf.drag coefficient
! VEG=V72X46.1.cor2   ! or:       ! vegetation fractions  (sum=1), need crops_yr=-1
VEG=V72X46.1.cor2_no_crops CROPS=CROPS_72X46N.cor4  ! veg. fractions, crops history
SOIL=S4X50093 TOPO=Z72X46N.cor4   ! soil/topography bdy.conds
REG=REG4X5                        ! special regions-diag
RVR=RD4X525.RVR                   ! river direction file
RADN1=sgpgxg.table8               ! rad.tables and history files
RADN2=kdist33.tautab8
RADN3=miescatpar.abcdv
RADN4=o3Prather1979-80.London1957-70
RADN5=TROAER.1875-1990.Jun2002
RADN6=dust8.tau9x8x13
RADN7=STRATAER.VOL.1850-1999.Apr02
RADN8=cloud.epsilon4.72x46
! RADN9=solar.lean99.uvflux        ! need KSOLAR<2
RADN9=solar.lean02.ann.uvflux      ! need KSOLAR=2
RADNA=O3.1850-2050.rec             ! with recovery of O3 after 2000
!  RADNA=O3.1850-2050.con          ! O3 'constant' after 2000
RADNB=o3WangJacob.1890.1979
RADNE=topcld.trscat8
GHG=GHG.1850-2050.Mar2002
dH2O=dH2O_by_CH4_monthly
TOP_INDEX=top_index_72x46.ij

Label and Namelist:
E001 (1979 atm/ocn)

DTFIX=300

&&PARAMETERS
! parameters set for prescribed ocean runs:
KOCEAN=0        ! ocn is prescribed
Kvflxo=1        ! save VFLXO (daily) if ocn prescribed
ocn_cycl=1      ! =0 if ocean varies from year to year

! parameters usually not changed when switching to q-flux ocean:

X_SDRAG=.00025,.000025  ! used above P(P)_sdrag mb (and in top layer)
C_SDRAG=0.      ! constant SDRAG above PTOP=150mb
P_sdrag=0.      ! linear SDRAG only in top layer (except near poles)
! PP_sdrag=20.  ! linear SDRAG above PP_sdrag mb near poles
ANG_sdrag=1     ! if 1: SDRAG conserves ang.momentum by adding loss below PTOP

xCDpbl=1.
U00ice=.60      ! U00ice up  => nethtz0 down (alb down) goals: nethtz0=0 (ann.
U00wtrX=.80     ! U00wtrX up => nethtz0 up   (alb down)           global mean)
HRMAX=550.      ! HRMAX up   => nethtz0 down (alb up  )        plan.alb 30%

RWCLDOX=1.5  !  wtr cld particle size *3/2 over ocean
RICLDX=.3333 !  ice cld particle size * 1(at 0mb)->1/3(at 1000mb)

CO2X=1.
H2OstratX=1.

H2ObyCH4=1.     ! activates strat.H2O generated by CH4
KVEGA6=3        ! 6-band albedo (Schramm)
KSOLAR=2

! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
crops_yr=1979 ! if -1, crops in VEG-file is used
s0_yr=1979
s0_day=182
ghg_yr=1979
ghg_day=182
volc_yr=1979
volc_day=182
aero_yr=1979
o3_yr=1979

DT_UVfilter=450.  ! usually same as DT (below)

! parameters that may have to be changed in emergencies:
DT=450.         ! from default: DTsrc=3600.,
NIsurf=2        ! increase as layer 1 gets thinner

! parameters that affect at most diagn. output:
Ndisk=24        ! use =240 on halem
SUBDD='SLP'     ! save SLP at sub-daily frequency
NSUBDD=12       ! saving sub-daily diags 12hrly
KCOPY=2         ! saving acc + rsf
isccp_diags=1   ! use =0 to save cpu time
nda5d=1         ! use =7 to save cpu time
nda5s=1         ! use =7 to save cpu time
&&END_PARAMETERS

 &INPUTZ
   YEARI=1949,MONTHI=12,DATEI=1,HOURI=0, ! IYEAR1=YEARI (default)
   YEARE=1956,MONTHE=1,DATEE=1,HOURE=0,     KDIAG=0,2,2,9*0,
   ISTART=2,IRANDI=0, YEARE=1949,MONTHE=12,HOURE=1,
 &END
