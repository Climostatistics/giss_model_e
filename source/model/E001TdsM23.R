E001TdsM23.R GISS Model E                                 gas 06/00

WARNING: The boundary conditions used here may not be what you want
         and no tuning has yet been done.
  Please check and see before running
E001TdsM23: new modelE (based on B402A - strat. version)
! example run deck with Drew Shindell's tropospheric
! chemistry and 15 tracers...
 T1=Ox
 T2=NOx
 T3=N2O5
 T4=HNO3
 T5=H2O2
 T6=CH3OOH
 T7=HCHO
 T8=HO2NO2
 T9=CO
 T10=CH4
 T11=PAN
 T12=Isoprene
 T13=AlkylNit
 T14=Alkenes
 T15=Paraffin

Preprocessor Options
#define TRACERS_ON                  ! include tracers code
#define TRACERS_WATER               ! include water tracers code
#define TRACERS_SPECIAL_Shindell
End Preprocessor Options

Object modules: (in order of decreasing priority)
RES_M23                             ! horiz/vert resolution
MODEL_COM GEOM_B IORSF              ! model variables and geometry
MODELE                              ! Main and model overhead
PARAM PARSER                        ! parameter database
ATMDYN_COM ATMDYN MOMEN2ND          ! atmospheric dynamics
STRATDYN STRAT_DIAG                 ! strospheric dynamics (incl. gw drag)
QUS_COM QUSDEF QUS_DRV              ! advection of tracers
TQUS_DRV                            ! advection of Q
TRACER_COM TRACERS_DRV              ! common and driver for tracers
! ---TRACER SPECIFIC CODES----------
TRACERS_SPECIAL_Shindell            ! routines specific to drew's 15-tracers
TRCHEM_Shindell_COM                 ! Drew Shindell's tracers common
TRCHEM_calc                         ! chemical reaction calculations
TRCHEM_init                         ! chemistry initialization, I/O
TRCHEM_family                       ! tracer family chemistry
TRCHEM_fastj                        ! tracer chem photlysis code/rad transf
TRCHEM_master                       ! trop chem "driver"/strat prescrioption
! ----------------------------------
TRACERS                             ! generic tracer code
TRDIAG_COM TRACER_PRT               ! tracer diagnostic printout
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
OCEAN OCNML                         ! ocean modules
ICEDYN_DRV ICEDYN  ! or: ICEDYN_DUM ! dynamic ice modules
SNOW_DRV SNOW                       ! snow model
RAD_COM RAD_DRV RADIATION           ! radiation modules
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
CONST FFT72 UTILDBL SYSTEM          ! utilities
POUT                                ! post-processing output

Data input files:
AIC=DEC1986.rsfB436TM23.modelE.17
! OHT=OTSPEC.RB399AM12.M250D   !  not needed if KOCEAN=0
OCNML=Z1O.B4X5.cor ! needed for post-processing only
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
! RADN9=solar.lean99.uvflux          ! need KSOLAR<2
RADN9=solar.lean02.ann.uvflux    ! need KSOLAR=2
RADNA=O3.1850-2050.depl.rec       ! with recovery of O3 after 2000
!  RADNA=O3.1850-2050.depl.con    ! O3 'constant' after 2000
RADNB=o3WangJacob.1890.1979
RADNE=topcld.trscat8
GHG=GHG.1850-2050.Mar2002
dH2O=dH2O_by_CH4_monthly
TOP_INDEX=top_index_72x46.ij
!-----------------------------------------------
MOLEC=chem_files/ds3ch4_moleculesXX
JPLRX=chem_files/gs_jpl00_trop_15_fix
JPLPH=chem_files/ds_photlist_trop_15
RATJ=chem_files/ratj.giss_15
SPECFJ=chem_files/jv_spec00_15.dat
ATMFJ=chem_files/jv_atms.dat
Ox_IC=gsin/Ox_init_cond_M23_4x5 !see README in /usr/people/cmrun/gsin
Ox_corr=gsin/corrOx_M23         ! see README in /usr/people/cmrun/gsin
CO_INDUSTRIAL=CO_sources/co_industrial
CO_BIOMASS=CO_sources/co_biomass_new
Alkenes_INDUSTRIAL=gsin/alkenes_industrial_4x5   !monthly  KG/hr/grid
Alkenes_BIOMASS=gsin/alkenes_biomass_4x5         !monthly  KG/hr/grid
Alkenes_VEGETATION=gsin/alkenes_vegetation_4x5   !monthly  KG/hr/grid
Paraffin_INDUSTRIAL=gsin/paraffin_industrial_4x5 !monthly  KG/hr/grid
Paraffin_BIOMASS=gsin/paraffin_biomass_4x5       !monthly  KG/hr/grid
Paraffin_VEGETATION=gsin/paraffin_vegetation_4x5 !monthly  KG/hr/grid
Isoprene_VEGETATION=gsin/isoprene_vegetation_4x5 !monthly  KG/hr/grid
NOx_FOSSIL_FUEL=NOy_sources/fossil_fuel_4x5
NOx_BIOMASS=NOy_sources/bio_burning_4x5
NOx_SOIL=NOy_sources/soil_nox_4x5
NOx_AIRCRAFT=NOy_sources/aircraft_4x5
SULFATE_SA=NOy_sinks/sulfate_fakeM23_M_SA
!-----------------------------------------------
CH4_ANIMALS=methane/gcm_data/CH4ANIMLS_4X5    ! Annual 1.3847
CH4_COALMINE=methane/gcm_data/CH4COAL_4X5      ! Annual 1.0285
CH4_GASLEAK=methane/gcm_data/CH4GASLEAK_4X5   ! Annual 3.904
CH4_GASVENT=methane/gcm_data/CH4GASVENT_4X5   ! Annual 1.659
CH4_CITYDUMP=methane/gcm_data/CH4MSW_4X5       ! Annual 1.233
CH4_SOIL_ABS=methane/gcm_data/CH4SOILABS_4X5   ! Annual 1.194
CH4_TERMITES=methane/gcm_data/CH4TRMITE_4X5    ! Annual 0.999
CH4_COALBURN=methane/gcm_data/COAL_BURN_BY_POP84_4X5  ! Annual 7.2154
CH4_BURN=methane/gcm_data/CH4BURN_4X5      ! Monthly 0.4369
CH4_RICE=methane/gcm_data/CH4RICEC_4X5     ! Monthly 0.7533
CH4_WETL=methane/gcm_data/CH4WETL+TUNDRA_4X5  ! Monthly 0.9818; zonal also

Label and Namelist:
E001TdsM23 (new modelE based on B402A - strat. version)
R=00BG/B

&&PARAMETERS
X_SDRAG=.00025,.000025  ! used for lin. sdrag above P_SDRAG mb
C_SDRAG=0.     ! no constant sdrag
P_SDRAG=.1     ! lin. sdrag above .1mb (top 2 layers) except near poles
PP_SDRAG=1.    ! lin. sdrag above 1.mb near poles (top 4 layers)
ANG_SDRAG=1    ! if =1: sdrag conserves ang mom.
PBREAK = 200.  ! The level for GW breaking above.
DEFTHRESH=0.000035 !the default is 15d-6
PCONPEN=500.   ! penetrating convection defn for GWDRAG
CMC = 0.0000003 

KOCEAN=0
U00ice=.75   ! tune this first to get reas.alb/cldcvr (range: .4-.6), then
HRMAX=500.   ! tune this to get rad.equilibrium (range: 100.-1500. meters)

H2ObyCH4=1.  ! activates strat.H2O generated by CH4
KVEGA6=3     ! 6-band albedo (Schramm)
KSOLAR=2

LMCM=16              ! max level of moist convection
XCDNST=300.,10000.   ! strat. gw drag parameters
DT=180.,             ! from default: DTsrc=3600.,
dt_UVfilter=180.
NIsurf=4,            ! number of surface time steps

NSUBDD=0        ! saving sub-daily diags 12hrly
Kvflxo=0        ! saving VFLXO (daily)
KCOPY=2         ! saving acc + rsf
isccp_diags=0
&&END_PARAMETERS

 &INPUTZ
   YEARI=1950,MONTHI=1,DATEI=1,HOURI=0,  !  from default: IYEAR1=YEARI
   YEARE=1956,MONTHE=1,DATEE=1,HOURE=0,  KDIAG=0,2,2,9*0,
   YEARE=1950,MONTHE=2,
   ISTART=7,IRANDI=0, YEARE=1950,MONTHE=1,HOURE=1,IWRITE=1,JWRITE=1,
 &END
