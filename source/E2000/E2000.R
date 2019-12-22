E2000.R GISS Model E SI2000 pbl/U00-cld      rar 12/02/02

E2000: modelE version of SI2000 physics
E034LM12: E034E but pbl-scheme and cloud-U00 as in SI2000
E034DM12: E034 but particle size of ice clouds strongly reduced, strong puddl. in SH
E034M12: E033M12 but ocn density increased by 3% and improved
         ice-water heat transports in lakes
E033M12: E032IM12 but saving heat conv. for ice dynamics
E032IM12: E032AM12 but water vapor and CO2 reduced by 20% for radiation
E032AM12: E032M12 but increasing albedos of shallow lakes

Preprocessor Options
!#define TRACERS_ON                  ! include tracers code
End Preprocessor Options

Object modules: (in order of decreasing priority)
RES_M12                             ! horiz/vert resolution
MODEL_COM GEOM_B IORSF              ! model variables and geometry
MODELE                              ! Main and model overhead
PARAM PARSER                        ! parameter database
ATMDYN_COM ATMDYN MOMEN2ND          ! atmospheric dynamics
QUS_COM QUSDEF QUS_DRV              ! advection of tracers
TQUS_DRV                            ! advection of Q
CLOUDS_SI2000 CLOUDS_DRV_SI2000 CLOUDS_COM        ! clouds modules
SURFACE FLUXES                      ! surface calculation and fluxes
GHY_COM GHY_DRV GHY                 ! land surface and soils
PBL_COM PBL_DRV_SI2000 PBL_SI2000          ! atmospheric pbl
! pick exactly one of the next 2 choices: ATURB or DRYCNV
! ATURB                             ! turbulence in whole atmosphere
DRYCNV                              ! drycnv
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_DRV                 ! land ice modules
OCEAN OCNML                         ! ocean modules
ICEDYN_DRV ICEDYN
SNOW                                ! snow model
RAD_COM RAD_DRV_rnd RADIATION          ! radiation modules
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
CONST_SI2000 FFT72 UTILDBL SYSTEM          ! utilities
POUT                                ! post-processing output

Data input files:
! AIC=DEC1958.rsfB394M12.modelE.16 ! model init cond (atm. and ground) ISTART=7
AIC=AIC.RES_M12.D771201   ! observed init cond   (atm. only)       ISTART=2
GIC=GIC.NOV1956.rsfB357M12 ! initial ground conditions              ISTART=2
! OHT=OTSPEC.RunIDM12.M250D ! hor.heat transp.  for q-flux ocean only
OCNML=Z1O.B4X5.cor         ! mixed layer depth,needed for post-processing only
OSST=OST4X5.B.1946-55avg.Hadl1.1 ! prescr. climatological ocean (1 yr of data)
SICE=SICE4X5.B.1946-55avg.Hadl1.1 ! prescr. climatological sea ice
CDN=CD4X500S VEG=V72X46.1.cor
SOIL=S4X50093 TOPO=Z72X46N.cor4 ! bdy.cond
REG=REG4X5           ! special regions-diag
RVR=RD4X525.RVR      ! river direction file
RADN1=sgpgxg.table8    ! rad.tables
RADN2=kdist33.tautabs4
RADN3=miescatpar.abcdv
RADN4=o3Prather1979-80.London1957-70
RADN5=TROAER.1875-1990.Jun2002
RADN6=dust8.tau9x8x13
RADN7=STRATAER.VOL.1850-1999.Apr02
RADN8=cloud.epsilon4.72x46
! RADN9=solar.lean99.uvflux         ! need KSOLAR<2
RADN9=solar.lean02.ann.uvflux  ! need KSOLAR=2
RADNA=O3.1850-2050.depl.con       ! with recovery of O3 after 2000
!  RADNA=O3.1850-2050.depl.con    ! O3 'constant' after 2000
RADNB=o3WangJacob.1890.1979
RADNE=topcld.trscat8
GHG=GHG.1850-2050.Mar2002
dH2O=dH2O_by_CH4
TOP_INDEX=top_index_72x46.ij

Label and Namelist:
E2000 (modelE SI2000)

DTFIX=300
&&PARAMETERS
X_SDRAG=.00025,.000025
C_SDRAG=0.
ANG_sdrag=1     ! if 1: SDRAG conserves ang.momentum by adding loss below PTOP

CO2X=.8
H2OstratX=.8

KOCEAN=0
! Nrad=1
xCDpbl=1.       ! tune surface mom.drag to get reas. SLP (990mb at 65S)
U00ice=.60
U00wtrX=1.
HRMAX=300.  ! tune this to get rad.equilibrium (range: 100.-1500. meters)
snoage_fac_max=0. ! disable Loth-Graf snow aging
KVEGA6=-4       ! 6-band albedo (Schramm)

KSOLAR=2
DT=450.,        ! from default: DTsrc=3600.,
DT_UVfilter=450.
SUBDD='SLP'     ! save SLP at sub-daily frequency
NSUBDD=0        ! saving sub-daily diags 12hrly
Kvflxo=1        ! saving VFLXO (daily)
KCOPY=2         ! saving acc + rsf
isccp_diags=0

! if the params below change, you may have to adjust U00wtrX
s0_yr=1951
ghg_yr=1951
ghg_day=182
s0_day=182
volc_yr=1951
volc_day=182
aero_yr=1951
o3_yr=1951
Ndisk=24
NDAA=1
NIPRNT=12
nrad=1
&&END_PARAMETERS

 &INPUTZ
   YEARI=1949,MONTHI=12,DATEI=1,HOURI=0,
   YEARE=1956,MONTHE=1,DATEE=1,HOURE=0,
   ISTART=2,IRANDI=0, YEARE=1949,MONTHE=12,HOURE=12,
 &END
