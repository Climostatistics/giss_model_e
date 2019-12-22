E001o.R GISS Model E                                 gas 06/00

E001o: new modelE (based on B402A - coupled version)

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
CLOUDS CLOUDS_DRV CLOUDS_COM        ! clouds modules
SURFACE FLUXES                      ! surface calculation and fluxes
GHY_COM GHY_DRV GHY                 ! land surface and soils
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
! pick exactly one of the next 2 choices:
! ATURB                             ! turbulence in whole atmosphere
DRYCNV                              ! drycnv
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_DRV                 ! land ice modules
ODIAG_COM OCEAN_COM OSTRAITS_COM OGEOM ! dynamic ocean modules
OCNDYN OSTRAITS OCNGM OCNKPP           ! dynamic ocean routines
ODIAG_PRT                              ! ocean diagnostic print out
OCNFUNTAB                           ! ocean function look up table
ICEDYN ICEDYN_DRV                   ! ice dynamics
SNOW                                ! snow model
RAD_COM RAD_DRV RADIATION           ! radiation modules
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
CONST FFT72 UTILDBL SYSTEM          ! utilities
POUT                                ! post-processing output

Data input files:
AIC=DEC1958.rsfB394M12.modelE.17.ext
OIC=OIC4X5LD.Z12.gas1.CLEV94.DEC01   ! ocean initial conditions
OFTAB=OFTABLE_NEW                    ! ocean function table
AVR=AVR72X46.L13.gas1.modelE         ! ocean filter
KBASIN=KB4X513.OCN.gas1              ! ocean basin designations
TOPO_OC=Z72X46N_gas.1 ! ocean bdy.cond
CDN=CD4X500S.ext VEG=V72X46.1.cor2.ext
SOIL=S4X50093.ext TOPO=Z72X46N_gas.1 ! bdy.cond
REG=REG4X5           ! special regions-diag
RVR=RD4X525.gas1.RVR      ! river direction file
RADN1=sgpgxg.table8    ! rad.tables
RADN2=kdist33.tautab8
RADN3=miescatpar.abcdv
RADN4=o3Prather1979-80.London1957-70
RADN5=TROAER.1875-1990.Jun2002
RADN6=dust8.tau9x8x13
RADN7=STRATAER.VOL.1850-1999.Apr02
RADN8=cloud.epsilon4.72x46
!  RADNA=O3.1850-2050.depl.con    ! O3 'constant' after 2000
RADN9=solar.lean99.uvflux          ! need KSOLAR<2
! RADN9=solar.lean02.ann.uvflux    ! need KSOLAR=2
RADNA=O3.1850-2050.depl.rec       ! with recovery of O3 after 2000
RADNB=o3WangJacob.1890.1979
RADNE=topcld.trscat8
GHG=GHG.1850-2050.Mar2002
dH2O=dH2O_by_CH4_monthly
TOP_INDEX=top_index_72x46.ij.ext

Label and Namelist:
E001o (new modelE based on B402A - coupled version)
R=00BG/B

&&PARAMETERS
X_SDRAG=.00025,.000025
C_SDRAG=0.
KOCEAN=1
U00ice=.60   ! use same value as corresponding run with prescribed ocean
U00wtr=.80   ! use same value as corresponding run with prescribed ocean
HRMAX=550.   ! use same value as corresponding run with prescribed ocean
KSOLAR=1
DT=450.,        ! from default: DTsrc=3600.,
dt_UVfilter=450.
NSUBDD=0        ! saving SLP 0hrly
Kvflxo=0        ! not saving VFLXO (daily)
KCOPY=2         ! saving acc + rsf
isccp_diags=1
&&END_PARAMETERS

 &INPUTZ
   YEARI=1950,MONTHI=1,DATEI=1,HOURI=0, !  from default: IYEAR1=YEARI
   YEARE=1956,MONTHE=1,DATEE=1,HOURE=0, KDIAG=0,2,2,9*0,
   YEARE=1950,MONTHE=2,
   ISTART=6,IRANDI=0, YEARE=1950,MONTHE=1,HOURE=1,IWRITE=1,JWRITE=1,
 &END
