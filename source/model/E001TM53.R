E001TM53.R GISS Model E                                 gas 06/00

WARNING: The boundary conditions used here may not be what you want
         and no tuning has yet been done.
  Please check and see before running
E001TM53: 53 layer 4x5 model sample tracers
 Air mass, SF6, RN222, CO2, 14CO2, CFC-11, CH4, N2O, linearizedO3

Preprocessor Options
#define TRACERS_ON                  ! include tracers code
#define TRACERS_SPECIAL_Lerner
End Preprocessor Options

Object modules: (in order of decreasing priority)
RES_M53                             ! horiz/vert resolution
MODEL_COM GEOM_B IORSF              ! model variables and geometry
MODELE                              ! Main and model overhead
PARAM PARSER                        ! parameter database
ATMDYN_COM ATMDYN MOMEN2ND          ! atmospheric dynamics
STRATDYN STRAT_DIAG                 ! strospheric dynamics (incl. gw drag)
QUS_COM QUSDEF QUS_DRV              ! advection of tracers
TQUS_DRV                            ! advection of Q
TRACER_COM TRACERS_DRV              ! configurable tracer code
TRACERS                             ! generic tracer code
TRDIAG_COM TRACER_PRT               ! tracer diagnostic printout
! use next line if #define TRACERS_SPECIAL_Lerner
TRACER_SPECIAL_Lerner               ! routines called when TRACERS_SPECIAL_Lerner is activated
CLOUDS CLOUDS_DRV CLOUDS_COM        ! clouds modules
SURFACE FLUXES                              ! surface calculation and fluxes
GHY_COM GHY_DRV GHY                 ! land surface and soils
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
! pick exactly one of the next 2 choices: ATURB or DRYCNV
! ATURB                             ! turbulence in whole atmosphere
DRYCNV                              ! drycnv
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_DRV                 ! land ice modules
OCEAN OCNML                         ! ocean modules
ICEDYN_DRV ICEDYN                   ! ice dynamics modules
SNOW_DRV SNOW                       ! snow model
RAD_COM RAD_DRV RADIATION           ! radiation modules
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
CONST FFT72 UTILDBL SYSTEM          ! utilities
POUT                                ! post-processing output

Data input files:
AIC=AIC.RES_M53.D771201
GIC=GIC.rsfB357M12.1DEC1956.1.ext ! initial conditions (ground)
! OHT=OTSPEC.RB399AM12.M250D ! not needed if KOCEAN=0
OCNML=Z1O.B4X5.cor   ! needed only for postprocessing
OSST=OST4X5.B.1946-55avg.Hadl1.1 SICE=SICE4X5.B.1946-55avg.Hadl1.1 ! ocn
CDN=CD4X500S VEG=V72X46.1.cor2
SOIL=S4X50093 TOPO=Z72X46N.cor4 ! bdy.cond
REG=REG4X5           ! special regions-diag
RVR=RD4X525.RVR      ! river direction file
ZVAR=ZVAR4X5         ! topographic variation for gwdrag
RADN1=sgpgxg.table8    ! rad.tables
RADN2=kdist33.tautab8
RADN3=miescatpar.abcdv
RADN4=o3Prather1979-80.London1957-70
RADN5=TROAER.1875-1990.Jun2002
RADN6=dust8.tau9x8x13
RADN7=STRATAER.VOL.1850-1999.Apr02
RADN8=cloud.epsilon4.72x46
! RADN9=solar.lean99.uvflux       ! need KSOLAR<2
RADN9=solar.lean02.ann.uvflux     ! need KSOLAR=2
RADNA=O3.1850-2050.depl.rec       ! with recovery of O3 after 2000
!  RADNA=O3.1850-2050.depl.con    ! O3 'constant' after 2000
RADNB=o3WangJacob.1890.1979
RADNE=topcld.trscat8
GHG=GHG.1850-2050.Mar2002
dH2O=dH2O_by_CH4_monthly
TOP_INDEX=top_index_72x46.ij
CO2_IC=CO2ijl_IC_Jan1_scale334_M23  !wofsy+B140TQaM9
CO2_FOS_FUEL=CO2_sources/gcm_data/CO2FOS_MRL_4X5
CO2_FERT=CO2_sources/gcm_data/CO2fert01_4X5
CO2_REGROWTH=CO2_sources/gcm_data/CO2_Nforest_4X5
CO2_LAND_USE=CO2_sources/gcm_data/CO2DEF_HOU_4X5
CO2_VEG=CO2_sources/gcm_data/CO2VEG_MON_4X5          ! Monthly source
CO2_OCEAN=CO2_sources/gcm_data/CO2_4X5_Ocean_flux02  ! Monthly source
14CO2_IC_DATA=workshop.14co2                         ! for 14CO2 Oct. 1963
LINOZ_TABLE=chemtab_solfb_31aug01.txt       ! linoz O3 coefficients
N2O_TABLE=N2Oloss.table                     ! Stratosphere tracer forcing
CFC11_TABLE=F11loss.table                   ! Stratosphere tracer forcing
CH4_TABLE=CH4chem.table                     ! Stratosphere tracer forcing
CH4_TROP_FRQ=CLIM.RUN.OHCH4.FRQ      !tropo loss frequency table (9 layers, n-grid)
N2O_IC=N2O_Shindell_Jan9293_M23             !initial conditions
CH4_IC=Wofsy_data_CH4       !wofsy jl initial conditions
CH4_ANIMALS=methane/gcm_data/CH4ANIMLS_4X5      ! Annual
CH4_COALMINE=methane/gcm_data/CH4COAL_4X5       ! Annual
CH4_GASLEAK=methane/gcm_data/CH4GASLEAK_4X5     ! Annual
CH4_GASVENT=methane/gcm_data/CH4GASVENT_4X5     ! Annual
CH4_CITYDUMP=methane/gcm_data/CH4MSW_4X5        ! Annual
CH4_SOIL_ABS=methane/gcm_data/CH4SOILABS_4X5    ! Annual
CH4_TERMITES=methane/gcm_data/CH4TRMITE_4X5     ! Annual
CH4_COALBURN=methane/gcm_data/COAL_BURN_BY_POP84_4X5   ! Annual
CH4_BURN=methane/gcm_data/CH4BURN_4X5           ! Monthly
CH4_RICE=methane/gcm_data/CH4RICEC_4X5          ! Monthly
CH4_WETL=methane/gcm_data/CH4WETL+TUNDRA_4X5    ! Monthly

Label and Namelist:
E001TM53 (4x5, 53 layer model, sample tracers)
R=00BG/B

&&PARAMETERS
X_SDRAG=.0005,.00005  ! used for lin. sdrag above P_SDRAG mb
C_SDRAG=0.     ! no constant sdrag
P_SDRAG=.1     ! lin. sdrag above .1mb (top 2 layers) except near poles
PP_SDRAG=1.    ! lin. sdrag above 1.mb near poles (top 4 layers)
ANG_SDRAG=1    ! if =1: sdrag conserves ang mom.

KOCEAN=0
U00ice=.85   ! tune this first to get reas.alb/cldcvr (range: .4-.6), then
HRMAX=400.   ! tune this to get rad.equilibrium (range: 100.-1500. meters)

H2ObyCH4=1.  ! activates strat.H2O generated by CH4
KVEGA6=3     ! 6-band albedo (Schramm)
KSOLAR=2

LMCM=26              ! max level of moist convection
XCDNST=300.,10000.   ! strat. gw drag parameters
DT=180.,             ! from default: DTsrc=3600.,
dt_UVfilter=180.,
NIsurf=4,            ! number of surface time steps

NSUBDD=0        ! saving sub-daily diags
Kvflxo=0        ! saving VFLXO (daily)
KCOPY=2         ! saving acc + rsf
isccp_diags=0

DEFTHRESH=0.00003   ! deformation threshold (default = 15d-6)
PBREAK=200.         ! p level for breaking gravity waves
CDEF=3.             ! parameter for GW DEF drag
CMTN=.5             ! parameter for GW MTN drag
CMC=0.0000002       ! parameter for GW Moist Convective drag

to_volume_MixRat=1,1,1,1,1,1,1,1,1   ! for tracer printout
&&END_PARAMETERS

 &INPUTZ
   YEARI=1950,MONTHI=1,DATEI=1,HOURI=0,  !  from default: IYEAR1=YEARI
   YEARE=1956,MONTHE=1,DATEE=1,HOURE=0,    KDIAG=0,2,2,9*0,
   YEARE=1950,MONTHE=2,
   ISTART=2,IRANDI=0, YEARE=1950,MONTHE=1,HOURE=1,IWRITE=1,JWRITE=1,
 &END
