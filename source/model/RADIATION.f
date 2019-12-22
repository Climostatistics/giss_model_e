        
      MODULE SURF_ALBEDO
!@sum SURF_ALBEDO contains parameters/variables needed for albedo calc
!@auth A. Lacis/V. Oinas (modifications by I. Alienov/G. Schmidt)
      implicit none
      save
!@var NKBAND number of K-bands
      integer, parameter :: NKBAND=33

!@var NVEG number of real vegetation types (not including bare soil)
!@var NV total number of vegetation types
      integer, parameter :: NVEG = 9, NV=11

!@var SRFOAM look up table for ocean foam as a function of wind speed
      real*8, parameter, dimension(25) :: SRFOAM = (/
     *     0.000,0.000,0.000,0.000,0.001,0.002,0.003,0.005,0.007,0.010,
     *     0.014,0.019,0.025,0.032,0.041,0.051,0.063,0.077,0.094,0.112,
     *     0.138,0.164,0.191,0.218,0.246/)

!@var SEASON julian day for start of season (used for veg albedo calc)
C                      1       2       3       4
C                    WINTER  SPRING  SUMMER  AUTUMN
      real*8, parameter, dimension(4)::
     *     SEASON=(/ 15.00,  105.0,  196.0,  288.0/)
C**** parameters used for vegetation albedo
!@var albvnd veg alb by veg type, season and band
      real*8, parameter :: ALBVND(NV,4,6) = RESHAPE( (/
C     (1)  >SRBALB(6) = VIS  (300 - 770 nm)
C        1     2     3     4     5     6     7     8     9    10    11
C      BSAND TNDRA GRASS SHRUB TREES DECID EVERG RAINF CROPS BDIRT ALGAE
     1 .500, .067, .089, .089, .078, .100, .067, .061, .089, .000, .200,
     2 .500, .062, .100, .100, .073, .055, .067, .061, .100, .000, .200,
     3 .500, .085, .091, .139, .085, .058, .083, .061, .091, .000, .200,
     4 .500, .080, .090, .111, .064, .055, .061, .061, .090, .000, .200,
C
C     (2)  >SRBALB(5) = NIR  (770 - 860 nm)    (ANIR=Ref)
C        1     2     3     4     5     6     7     8     9    10    11
C      BSAND TNDRA GRASS SHRUB TREES DECID EVERG RAINF CROPS BDIRT ALGAE
     1 .500, .200, .267, .267, .233, .300, .200, .183, .267, .000, .200,
     2 .500, .206, .350, .300, .241, .218, .200, .183, .350, .000, .200,
     3 .500, .297, .364, .417, .297, .288, .250, .183, .364, .000, .200,
     4 .500, .255, .315, .333, .204, .218, .183, .183, .315, .000, .200,
C
C     (3)  >SRBALB(4) = NIR  (860 -1250 nm)    (ANIR*1.0)
C        1     2     3     4     5     6     7     8     9    10    11
C      BSAND TNDRA GRASS SHRUB TREES DECID EVERG RAINF CROPS BDIRT ALGAE
     1 .500, .200, .267, .267, .233, .300, .200, .183, .267, .000, .200,
     2 .500, .206, .350, .300, .241, .218, .200, .183, .350, .000, .200,
     3 .500, .297, .364, .417, .297, .288, .250, .183, .364, .000, .200,
     4 .500, .255, .315, .333, .204, .218, .183, .183, .315, .000, .200,
C
C     (4)  >SRBALB(3) = NIR  (1250-1500 nm)    (ANIR*0.4)
C        1     2     3     4     5     6     7     8     9    10    11
C      BSAND TNDRA GRASS SHRUB TREES DECID EVERG RAINF CROPS BDIRT ALGAE
     1 .500, .080, .107, .107, .093, .120, .080, .073, .107, .000, .200,
     2 .500, .082, .140, .120, .096, .083, .080, .073, .140, .000, .200,
     3 .500, .119, .145, .167, .119, .115, .100, .073, .145, .000, .200,
     4 .500, .102, .126, .132, .081, .087, .073, .073, .126, .000, .200,
C
C     (5)  >SRBALB(2) = NIR  (1500-2200 nm)    (ANIR*0.5)
C        1     2     3     4     5     6     7     8     9    10    11
C      BSAND TNDRA GRASS SHRUB TREES DECID EVERG RAINF CROPS BDIRT ALGAE
     1 .500, .100, .133, .133, .116, .150, .100, .091, .133, .000, .200,
     2 .500, .103, .175, .150, .120, .109, .100, .091, .175, .000, .200,
     3 .500, .148, .182, .208, .148, .144, .125, .091, .182, .000, .200,
     4 .500, .127, .157, .166, .102, .109, .091, .091, .157, .000, .200,
C
C     (6)  >SRBALB(1) = NIR  (2200-4000 nm)    (ANIR*0.1)
C        1     2     3     4     5     6     7     8     9    10    11
C      BSAND TNDRA GRASS SHRUB TREES DECID EVERG RAINF CROPS BDIRT ALGAE
     1 .500, .020, .027, .027, .023, .030, .020, .018, .027, .000, .200,
     2 .500, .021, .035, .030, .024, .022, .020, .018, .035, .000, .200,
     3 .500, .030, .036, .042, .030, .029, .025, .018, .036, .000, .200,
     4 .500, .026, .032, .033, .020, .022, .018, .018, .032, .000, .200
     *     /), (/11,4,6/) )
C
!@var VTMASK vegetation depth mask by type (kg/m^2)
      real*8, parameter :: VTMASK(NV) = (/
C        1     2     3     4     5     6     7     8     9    10    11
C     BSAND TNDRA GRASS SHRUB TREES DECID EVERG RAINF CROPS BDIRT ALGAE
     * 10.0, 20.0, 20.0, 50.0, 200., 500.,1000.,2500., 20.0, 10.0, .001
     *     /)

!@var ASNALB snow albedo (original version)
!@var AOIALB seaice albedo (original version)
!@var ALIALB land ice albedo (original version)
C****  albedo = asnalb+snowage*snowage_fac so these numbers are the min
      real*8, parameter ::
C                        VIS  NIR1  NIR2  NIR3  NIR4  NIR5    NIR
     *     ASNALB(7)=(/.60d0,.55d0,.55d0,.30d0,.10d0,.05d0, .35d0/),
     *     AOIALB(7)=(/.55d0,.50d0,.45d0,.25d0,.10d0,.05d0, .30d0/),
     *     ALIALB(7)=(/.60d0,.55d0,.50d0,.30d0,.10d0,.05d0, .35d0/),
C**** shorthand for the 2 band version
     *     ASNVIS=ASNALB(1), ASNNIR=ASNALB(7),
     *     AOIVIS=AOIALB(1), AOINIR=AOIALB(7),
     *     ALIVIS=ALIALB(1), ALINIR=ALIALB(7)

C**** variables that control original snow aging calculation
!@var AGEXPF exponent in snowage calculation depends on hemi/surf type
!@var ALBDIF difference in albedo as function of snowage
      REAL*8, PARAMETER, DIMENSION(3,2) ::
     *     AGEXPF = RESHAPE( (/
C          SH EA   SH OC   SH LI   NH EA   NH OC   NH LI
     *     0.2d0,  0.2d0,  0.2d0,  0.2d0,  0.2d0,  0.2d0 /), (/3,2/) ),
     *     ALBDIF = RESHAPE( (/
C          SH EA   SH OC   SH LI   NH EA   NH OC   NH LI
     *     0.35d0, 0.35d0, 0.35d0, 0.35d0, 0.35d0, 0.35d0/), (/3,2/) )

!@var DMOICE, DMLICE masking depth for snow on ice and land ice
      real*8, parameter :: DMOICE = 10., DMLICE = 10.

!@var AOCEAN K-band dependent Thermal radiation characteristics for ocn
      real*8, parameter, dimension(NKBAND) :: AOCEAN = (/
     +        0.04000,0.09566,0.10273,0.10389,0.10464,0.10555,0.10637,
     +        0.10666,0.10697,0.10665,0.10719,0.10728,0.11007,0.04009,
     +        0.04553,0.05554,0.08178,0.09012,0.09464,0.09548,0.09532,
     +        0.09558,0.09558,0.09568,0.09565,0.05771,0.04985,0.04670,
     +        0.04630,0.04575,0.04474,0.04468,0.04500/)

!@var AGSIDV K-band dependent Thermal radiation for other types
C                     AGSNOW  AGLICE  AGROCK  AGVEG
      real*8, parameter, dimension(NKBAND,4) :: AGSIDV = RESHAPE( (/
     +        0.01400,0.09262,0.09170,0.07767,0.07130,0.06603,0.06540,
     +        0.06397,0.06358,0.06361,0.06365,0.06386,0.06564,0.01354,
     +        0.01537,0.02320,0.04156,0.03702,0.03633,0.03417,0.03346,
     +        0.03342,0.03322,0.03350,0.03170,0.01967,0.01845,0.01977,
     +        0.01986,0.01994,0.02013,0.02041,0.02100,
     +        0.01400,0.09262,0.09170,0.07767,0.07130,0.06603,0.06540,
     +        0.06397,0.06358,0.06361,0.06365,0.06386,0.06564,0.01354,
     +        0.01537,0.02320,0.04156,0.03702,0.03633,0.03417,0.03346,
     +        0.03342,0.03322,0.03350,0.03170,0.01967,0.01845,0.01977,
     +        0.01986,0.01994,0.02013,0.02041,0.02100,
     +        0.04500,0.10209,0.08806,0.05856,0.04835,0.04052,0.04001,
     +        0.03775,0.03687,0.03740,0.03637,0.03692,0.03570,0.07001,
     +        0.05665,0.05326,0.05349,0.04356,0.03845,0.03589,0.03615,
     +        0.03610,0.03602,0.03613,0.03471,0.13687,0.14927,0.16484,
     +        0.16649,0.16820,0.17199,0.17484,0.18000,
     +        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     *        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. /),
     *        (/ NKBAND,4 /) )

!@var AVSCAT,ANSCAT,AVFOAM,ANFOAM for ocean albedo calc
      real*8, parameter ::
     *     AVSCAT=0.0156d0, ANSCAT=0d0, AVFOAM=0.2197d0, ANFOAM=0.1514d0

C**** miscellaneous constants
      real*8, parameter ::
!@var WETTRA,WETSRA adjustment factors for wet earth albedo calc
     *     WETTRA=1.0, WETSRA=1.0,
!@var ZOCSRA,ZSNSRA,ZICSRA,ZDSSRA,ZVGSRA adjustment factors for
!@+   solar zenith angle effects (not fully enabled)
     *     ZOCSRA=1.0, ZSNSRA=1.0, ZICSRA=1.0, ZDSSRA=1.0, ZVGSRA=1.0,
!@var EOCTRA,ESNTRA,EICTRA,EDSTRA,EVGTRA adjustment factors for
!@+   thermal radiation effects  (not fully enabled)
     *     EOCTRA=1.0, EDSTRA=1.0, ESNTRA=1.0, EICTRA=1.0, EVGTRA=1.0

C**** ALBVNH is set only once a day then saved
!@var ALBVNH hemispherically varying vegetation albedo
      real*8, dimension(NV,6,2) :: ALBVNH


!@var GZSNOW asymmetry parameter for snow over three types
!@+   from Wiscombe and Warren (1980) JAS
!@+   Note this is used for ice + melt-ponds as well.
      REAL*8, PARAMETER, DIMENSION(7,3,2) :: GZSNOW = RESHAPE( (/
C       VIS     NIR1    NIR2     NIR3     NIR4     NIR5    NIRT
     * 0.95d0, 0.94d0, 0.905d0, 0.896d0, 0.894d0, 0.89d0, 0.91d0,! SH EA
     * 0.95d0, 0.94d0, 0.905d0, 0.896d0, 0.894d0, 0.89d0, 0.91d0,!    OI
     * 0.95d0, 0.94d0, 0.905d0, 0.896d0, 0.894d0, 0.89d0, 0.91d0,!    LI
     * 0.95d0, 0.94d0, 0.905d0, 0.896d0, 0.894d0, 0.89d0, 0.91d0,! NH EA
     * 0.95d0, 0.94d0, 0.905d0, 0.896d0, 0.894d0, 0.89d0, 0.91d0,!    OI
     * 0.95d0, 0.94d0, 0.905d0, 0.896d0, 0.894d0, 0.89d0, 0.91d0 !    LI
     *     /), (/7,3,2/) )

      END MODULE SURF_ALBEDO

      MODULE RE001
!@sum radiation module based originally on rad00b.radcode1.F
!@auth A. Lacis/V. Oinas/R. Ruedy

      INCLUDE 'BR00B.COM'
      REAL*4 yr1S0,yr2S0
!@dbparm snoage_fac_max  max snow age reducing-factor for sea ice albedo
      REAL*8 :: snoage_fac_max=.5d0
!@dbparm kzsnow =1 for snow/ice albedo zenith angle depend. (default=0)
Ccb   INTEGER :: kzsnow=0
      INTEGER :: JYEAR=1980,JDAY=1
      SAVE

      CONTAINS

!!    SUBROUTINE SETNEW(NIPSET,NVALUE,IJYEAR,IJDAY,ILAYER,XVALUE)
!!
!!    !INCLUDE 'BR00B.COM'
!!
!!    IF(NIPSET.EQ.0) THEN
!!
C     ------------------------------------------------------------------
C     MADVEL  Model Add-on Data of Extended Climatology Enable Parameter
C             Each MADVEL digit is ON/OFF switch for corresponding input
C                         123456
C             e.g. MADVEL=111111   (zero digit skips input process)
C
C     MADO3M   =  1   Reads  Makiko's 1951-1997 Ozone climatology RFILEA
C                            Wang-Jacob 1890-1979 Trop O3 change  RFILEB
C     MADAER   =  1   Reads  Aerosol 50y tropospheric climatology RFILE5
C     MADDST   =  1   Reads  Dust-windblown mineral climatology   RFILE6
C     MADVOL   =  1   Reads  Volcanic 1950-00 aerosol climatology RFILE7
C              =  2   Reads  Volcanic 1950-00 aerosol climatology RFILE7
C     MADEPS   =  1   Reads  Epsilon cloud heterogeneity data     RFILE8
C     MADLUV   =  1   Reads  Lean's SolarUV 1882-1998 variability RFILE9
C
C                 Related Model Add-on Data Parameters set in BLOCK DATA
C
C     MADGHG   =  1  Default Enables UPDGHG update. (MADGHG=0),no update
C     MADSUR   =  1   Reads  V72X46N.1.cor Vegetation type data   RFILEC
C                            Z72X46N Ocean fraction, topography   RFILED
C     ------------------------------------------------------------------
!!
!!    MADVEL=NVALUE
!!    MADLUV=MADVEL                      !   Reset MAD Enable Parameters
!!    NORMAL=100000
!!    MADO3M=MADLUV/NORMAL                     !       Ozone Climatology
!!    MADLUV=MADLUV-MADO3M*NORMAL
!!    NORMAL=NORMAL/10
!!    MADAER=MADLUV/NORMAL                     !    Trop Aer Climatology
!!    MADLUV=MADLUV-MADAER*NORMAL
!!    NORMAL=NORMAL/10
!!    MADDST=MADLUV/NORMAL                     !        Dust Climatology
!!    MADLUV=MADLUV-MADDST*NORMAL
!!    NORMAL=NORMAL/10
!!    MADVOL=MADLUV/NORMAL                     !    Volc Aer Climatology
!!    MADLUV=MADLUV-MADVOL*NORMAL
!!    NORMAL=NORMAL/10
!!    MADEPS=MADLUV/NORMAL                     !     Cld Eps Climatology
!!    MADLUV=MADLUV-MADEPS*NORMAL              !    Solar UV Climatology
!!    RETURN
!!    ENDIF

C          Internal Control Parameter (Initialization / Dynamic Control)
C          -------------------------------------------------------------
C          Permits external SETNEW CALL (with access to RADCOM) to reset
C                 FSPARE and ISPARE Time-Trend Parameters KYEARX, KJDAYX
C                 according to:
C                              NIPSET       Process       KYEARX  KJDAYX
C                                1       SolarCon, UV       0       0
C                                2       GH Gas Trend       0       0
C                                3       Ozone Distr        0       0
C                                4       AerClimtolgy       0       0
C                                5       Desert Dust        0       0
C                                6       Volcanic Aer       0       0
C                                7       Epsilon Clds       0       0
C                                8       Refl Surface       0       0
C                          ---------------------------------------------
!!    IF(NIPSET.LT.9) THEN
!!    NYSPAR=30+NIPSET
!!    NDSPAR=40+NIPSET
!!    ISPARE(NYSPAR)=IJYEAR
!!    ISPARE(NDSPAR)=IJDAY
!!    ENDIF

C     IF(NIPSET.EQ.9): Reset all 8 Time-Trend Parameters to IJYEAR,IJDAY
C     ------------------------------------------------------------------

!!    IF(NIPSET.EQ.9) THEN
!!    DO 110 NIP=1,8
!!    NYSPAR=30+NIP
!!    NDSPAR=40+NIP
!!    ISPARE(NYSPAR)=IJYEAR
!!    ISPARE(NDSPAR)=IJDAY
!!110 CONTINUE
!!    ENDIF
C                   ----------------------------------------------------
C                   Permits External SETNEW CALL (with Access to RADCOM)
C
C                   (11)    to Reset KSOLAR to value NVALUE if NIPSET=11
C                   (12)    to Reset KTREND to value NVALUE if NIPSET=12
C                   (13)    to Reset MADGHG to value NVALUE if NIPSET=13
C                   (14)    to Reset MADSUR to value NVALUE if NIPSET=14
C                   ----------------------------------------------------
!!
!!    IF(NIPSET.EQ.11) THEN
!!    KSOLAR=NVALUE
!!    IF(IJYEAR.GT.0) ISPARE(31)=IJYEAR
!!    IF(IJDAY.GT.0)  ISPARE(41)=IJDAY
!!    ENDIF
!!    IF(NIPSET.EQ.12) THEN
!!    KTREND=NVALUE
!!    IF(IJYEAR.GT.0) ISPARE(32)=IJYEAR
!!    IF(IJDAY.GT.0)  ISPARE(42)=IJDAY
!!    ENDIF
!!    IF(NIPSET.EQ.13) THEN
!!    MADGHG=NVALUE
!!    ENDIF
!!    IF(NIPSET.EQ.14) THEN
!!    MADSUR=NVALUE
!!    ENDIF
!!
!!
C     ---------------------------------------------------------------
C     GHG Concentration Direct Specification  (for use with KTREND=0)
C                        Permits direct setting of GHG Default Values
C                        --------------------------------------------
C
C     GAS  NUMBER  1      2     3      4      5       6       7
C                 H2O    CO2    O3    O2     NO2     N2O     CH4
C     DATA PPMVK0/0.0, 337.90, 0.0, 210000., 0.0,  0.3012, 1.5470
C
C     GAS  NUMBER  8           9        10        11         12
C                CCL3F1      CCL2F2     N2       CFC-Y      CFC-Z
C    +        ,0.1666E-03, 0.3003E-03, 0.0,   0.0978D-03, 0.0010D-06/
C
C     ---------------------------------------------------------------
C     Makiko's GHG Trend Compilation  GHG.1850-2050.Dec1999 in GTREND
C     GTREND:  1980.,  337.9,  .3012,  1.547,  .1666,  .3003,  .0978,
C     ---------------------------------------------------------------
!!
!!    IF(NIPSET.EQ.15) THEN
!!    PPMVK0(NVALUE)=XVALUE
!!    ENDIF
!!
!!
C          ----------------------------------------------------------
C          Option to permit multiplying aerosol vertical distribution
C          (VD reference layer) by scaling factor for BCI,OCI,SUI,DST
C          The scaling option is implemented in GETAER and GETDST, so
C          VDF factors can be set and re-set with each CALL to RCOMPX
C
C            Aerosol Vertical Profile (Optical Depth) Scaling Factors
C            --------------------------------------------------------
C     Layr  1    2    3    4    5    6    7    8    9    10   11   12
C     Ptop 934  854  720  550  390  285  255  150  100   60   30   10
C     Pbot 984  934  854  720  550  390  285  255  150  100   60   30
C
C     DATA VDFBCI/             (Reference values, ILAYER=1,2, ... 12)
C    A     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
C    ----------------------------------------------------------------
!!
!!    IF(NIPSET.EQ.16) THEN
!!    VDFBCI(ILAYER)=XVALUE         !   Black Carbon (Industrial) VDF
!!    ENDIF
!!
!!    IF(NIPSET.EQ.17) THEN
!!    VDFOCI(ILAYER)=XVALUE         ! Organic Carbon (Industrial) VDF
!!    ENDIF
!!
!!    IF(NIPSET.EQ.18) THEN
!!    VDFSUI(ILAYER)=XVALUE         !        Sulfate (Industrial) VDF
!!    ENDIF
!!
!!    IF(NIPSET.EQ.19) THEN
!!    VDFDST(ILAYER)=XVALUE         !                 Desert Dust VDF
!!    ENDIF
!!
!!    IF(NIPSET.EQ.20) THEN
!!    VDFBCI(ILAYER)=XVALUE         ! Set BCI,OCI,SUI,DST to same VDF
!!    VDFOCI(ILAYER)=XVALUE
!!    VDFSUI(ILAYER)=XVALUE
!!    VDFDST(ILAYER)=XVALUE
!!    ENDIF
!!
!!
C     ---------------------------------------------------------------
C     Option to directly reset ISPARE(ILAYER) Parameter to new NVALUE
C     ---------------------------------------------------------------
!!
!!    IF(NIPSET.EQ.21) THEN
!!    ISPARE(ILAYER)=NVALUE
!!    ENDIF
!!
C          ISPARE     PARAMTER     DEFAULT     PARAMETER DESCRIPTION
C             1        KVRAER         1        Repartition Aer VDist
C             2        MEANAC         0        Ann-Mean Aer Climtogy
C             3        MEANDD         0        Ann-Mean  Desert Dust
C             4        MEANVA         0        Ann-Mean Volcanic Aer
C             5        NCARO3         0        NCAR London 1976 Ozon
C             6        KUVFAC         0        ON/OFF UV Mult Factor
C             7        KSNORM         0        Norm S0 when KUVFAC=1
C             8        KWTRAB         0        WRITER: Qab,Qex,Qsc,g
C             9        KGGVDF         0        GG Vertical Prof Grad
C            10        KPGRAD         1        Pole-to-Pole GHG Grad
C            11        KLATZ0         1        GG VDist Latitude Dep
C            12        KCLDEM         1        Top Cloud Scatter Cor
C            13        KANORM         0        SGP Atmos Column Norm
C            14        KPFCO2         0        1=MOD CO2PROF: FPXCO2
C            15        KPFOZO         0        1=MOD O3 PROF: FPXOZO
C            16        KVEGA6         1        6 SpecBand Veg Albedo
C            17        KORDER         0        WRITER k-d spec order
C            18        KUFH2O         1        H2O Col AbsorbScaling
C            19        KUFCO2         1        CO2 Col AbsorbScaling
C            20        unused         1
C            21        KCSELF         1        H2O Cont Self-Broaden
C            22        KCFORN         1        CO2 Con Foreign-Broad
C          23-30       unused
!!
C     (Time-Trend Over-ride Factors to Select Year and Day of Trend)
C     --------------------------------------------------------------
C          ISPARE     PARAMTER     DEFAULT     PARAMETER DESCRIPTION
C            31        KYEARS         0        SolarTrend Reset Year
C            32        KYEARG         0        GHG  Trend Reset Year
C            33        KYEARO         0        Ozon Trend Reset Year
C            34        KYEARA         0        Aero Trend Reset Year
C            35        KYEARD         0        Dust Trend Reset Year
C            36        KYEARV         0        Volc Trend Reset Year
C            37        KYEARE         0        CldEpTrend Reset Year
C            38        KYEARR         0        RVegeTrend Reset Year
C            39        unused
C
C            41        KJDAYS         0        SolarTrend Reset Date
C            42        KJDAYG         0        GHG  Trend Reset Date
C            43        KJDAYO         0        Ozon Trend Reset Date
C            44        KJDAYA         0        Aero Trend Reset Date
C            45        KJDAYD         0        Dust Trend Reset Date
C            46        KJDAYV         0        Volc Trend Reset Date
C            47        KJDAYE         0        CldEpTrend Reset Date
C            48        KJDAYR         0        RVegeTrend Reset Date
C          48-98       unused
!!
C     ---------------------------------------------------------------
C     Option to directly reset FSPARE(ILAYER) Parameter to new XVALUE
C     ---------------------------------------------------------------
!!
!!    IF(NIPSET.EQ.22) THEN
!!    FSPARE(ILAYER)=XVALUE
!!    ENDIF
!!
C          FSPARE     PARAMTER     DEFAULT     PARAMETER DESCRIPTION
C           1-3        AGEXPF      3*0.200     SnowAge Exp Factor SH
C           4-6         "  "       3*0.200         (E, OI, LI)    NH
C           7-9        ALBDIF      3*0.200     Snow/Ice Alb Dif   SH
C          10-12        "  "       3*0.200         (E, OI, LI)    NH
!!
C            13         PSOIL                  Bare-Soil PEARTH Frac
C            14        REFF0        0.30           Default VAer Reff
C            15        VEFF0        0.35           Default VAer Veff
C            16        O3WJT0       0.0010     O3 LonDep Extrap Time
C            17        SNOLIM       0.3670     Min SnowMask Exposure
C            18        X0YBCI       0.0001     BCI Aer TimeDep Reset
C            19        X0YOCI       0.0001     OCI Aer TimeDep Reset
C            20        X0YSUI       0.0001     SUI Aer TimeDep Reset
C            21        AVSCAT       0.0156     Ocean VIS  Scat Param
C            22        ANSCAT       0.0002     Ocean NIR  Scat Param
C            23        AVFOAM       0.2197     Foam  VIS  Scat Param
C            24        ANFOAM       0.1514     Foam  NIR  Scat Param
!!
C          25-27       UVWAVL                  UV Sol Flx Partitions
C          28-30       UVFACT                  UV Solar Flux Factors
!!
C            31        PTOPTR       150.00     TropTop SIGMA-lev Pmb
C            32        unused
C            33        unused
C            34        unused
C          35-40       unused
!!
C          (Climatology Aerosol Maximum PI0 Reset Parameters)
C          --------------------------------------------------
C          FSPARE     PARAMETER   PI0VIS(max)     PI0TRA(max)
C          41  61       ACID1       1.00000          0.0
C          42  62       SSALT       1.00000          0.0
C          43  63       SLFT1       1.00000          0.0
C          44  64       SLFT2       1.00000          0.0
C          45  65       BSLT1        .98929          0.0
C          46  66       BSLT2        .95609          0.0
C          47  67       DUST1        .91995          0.0
C          48  68       DUST2        .78495          0.0
C          49  69       DUST3        .63576          0.0
C          50  70       CARB1        .31482          0.0
C          51  71       CARB2        .47513          0.0
!!
C          Vegetated Surface, Bare Soil Albedo Parameters
C          ----------------------------------------------
C          52  72   BVSOIL   TSOIL      (not implemented)
C          53  73   XVSOIL   TVEGE      (not implemented)
C          54  74   BNSOIL  FTRUFS      (not implemented)
C          55  75   XNSOIL  FTRUFV      (not implemented)
C          56  76   BVVEGE  DTRUFS      (not implemented)
C          57  77   XVVEGE  DTRUFV      (not implemented)
C          58  78   BNVEGE  unused      (not implemented)
C          59  79   XNVEGE  unused      (not implemented)
C          60  80   unused  unused
!!
C          (Aerosol Vertical Profile (Optical Depth) Scaling Factors)
C          ----------------------------------------------------------
C          FSPARE VDFBCI  FSPARE VDFOCI  FSPARE VDFSUI  FSPARE VDFDST
C            81    1.0     101    1.0      121   1.0     141    1.0
C            82    1.0     102    1.0      122   1.0     142    1.0
C            83    1.0     103    1.0      123   1.0     143    1.0
C            84    1.0     104    1.0      124   1.0     144    1.0
C            85    1.0     105    1.0      125   1.0     145    1.0
C            86    1.0     106    1.0      126   1.0     146    1.0
C            87    1.0     107    1.0      127   1.0     147    1.0
C            88    1.0     108    1.0      128   1.0     148    1.0
C            89    1.0     109    1.0      129   1.0     149    1.0
C            90    1.0     110    1.0      130   1.0     150    1.0
C            91    1.0     111    1.0      131   1.0     151    1.0
C            92    1.0     112    1.0      132   1.0     152    1.0
!!
C          FSPAER     PARAMETER
C          93-100      unused
C         113-120      unused
C         133-140      unused
C         153-160      unused
!!
C         (Greenhouse Gas 1980 Reference Amts, Vertical,Lat Gradients)
C         ------------------------------------------------------------
C          FSPARE         PPMVK0   FSPARE    PPMVDF  FSPARE   PPGRAD
C            161        .0000000     181    1.00000    201    .00000
C            162     337.9000000     182    1.00000    202    .00000
C            163        .0000000     183    1.00000    203    .00000
C            164  210000.0000000     184    1.00000    204    .00000
C            165        .0000000     185    1.00000    205    .00000
C            166        .3012000     186     .88888    206    .01000
C            167       1.5470000     187     .88888    207    .09000
C            168        .0001666     188     .88888    208    .06000
C            169        .0003003     189     .88888    209    .06000
C            170        .9000000     190    1.00000    210    .00000
C            171        .0000978     191     .88888    211    .06000
C            172        .0000000     192     .88888    212    .06000
!!
C          FSPAER     PARAMETER
C         173-180      unused
C         193-200      unused
C         213-998      unused
!!
!!
C     ----------------------------------------------------------------
C     Option to re-scale (9) Climatology aerosol column optical depths
C     separately for solar and thermal components of each aerosol type
C     by means of FSAERO, FTAERO scaling factors (activated in GETAER)
C
C     Defaults:
C                       TROPOSPHERIC AEROSOL PARAMETERS
C                  BCI  OCI  SUI  SEA  SUN    ANT  OCN  OCB  BCB  SSB
C     DATA FSAERO/ 0.5, 0.5, 1.0, 4.0, 1.0,   0.0, 1.0, 1.0, 1.0, 1.0/
C
C     DATA FTAERO/ 0.5, 0.5, 1.0, 4.0, 1.0,   0.0, 1.0, 1.0, 1.0, 1.0/
C
C     ----------------------------------------------------------------
C                                                (NVALUE=1,2,3, ... 9)
C
C     Climatology Aerosol SW (Solar) Component:
C     ----------------------------------------
!!
!!    IF(NIPSET.EQ.23) THEN
!!    FSAERO(NVALUE)=XVALUE
!!    ENDIF
!!
!!
C     Climatology Aerosol LW (Thermal) Component:
C     ------------------------------------------
!!
!!    IF(NIPSET.EQ.24) THEN
!!    FTAERO(NVALUE)=XVALUE
!!    ENDIF
!!
!!
C     ----------------------------------------------------------------
C     Option to re-scale (8) Desert Dust aerosol column optical depths
C     separately for solar and thermal components of each aerosol size
C     by means of FSDUST, FTDUST scaling factors (activated in SETDST)
C
C     Defaults:
C                          MINERAL DUST PARAMETERS
C                         CLAY                  SILT
C     Dust  size:  0.1, 0.2, 0.4, 0.8,   1.0, 2.0, 4.0, 8.0
C     DATA FSDUST/ 1.0, 1.0, 1.0, 1.0,   1.0, 1.0, 1.0, 1.0/
C
C     DATA FTDUST/ 1.0, 1.0, 1.0, 1.0,   1.0, 1.0, 1.0, 1.0/
C
C     ----------------------------------------------------------------
C                                                (NVALUE=1,2,3, ... 8)
C
C     Desert Dust SW (Solar) Size Component:
C     -------------------------------------
!!
!!    IF(NIPSET.EQ.25) THEN
!!    FSDUST(NVALUE)=XVALUE
!!    ENDIF
!!
!!
C     Desert Dust LW (Thermal) Size Component:
C     ---------------------------------------
!!
!!    IF(NIPSET.EQ.26) THEN
!!    FTDUST(NVALUE)=XVALUE
!!    ENDIF
!!
!!
C     ----------------------------------------------------------------
C     Option to re-scale Desert Dust total optical depth for solar and
C     thermal components by means of the FSDAER,FTDAER scaling factors
C                                              (operational in GETDST)
C     Defaults:
C                FSDAER = FSXAER(4) = 1.0     FTDAER = FTXAER(4) = 1.0
C
C     ----------------------------------------------------------------
!!
!!
C     Desert Dust SW (Solar) Component:
C     --------------------------------
!!
!!    IF(NIPSET.EQ.27) THEN
!!    FSXAER(4)=XVALUE
!!    ENDIF
!!
!!
C     Desert Dust LW (Thermal) Component:
C     ----------------------------------
!!
!!    IF(NIPSET.EQ.28) THEN
!!    FTXAER(4)=XVALUE
!!    ENDIF
!!
!!
C     ----------------------------------------------------------------
C     Option to re-scale  Volcanic aerosol optical depth for solar and
C     thermal components by means of the FSVAER,FTVAER scaling factors
C                                              (operational in GETVOL)
C     Defaults:
C                FSVAER = FSXAER(5) = 1.0     FTVAER = FTXAER(5) = 1.0
C
C     ----------------------------------------------------------------
C
C
C     Volcanic Aerosol SW (Solar) Component:
C     -------------------------------------
!!
!!    IF(NIPSET.EQ.29) THEN
!!    FSXAER(5)=XVALUE
!!    ENDIF
!!
!!
C     Volcanic Aerosol LW (Thermal) Component:
C     ---------------------------------------
!!
!!    IF(NIPSET.EQ.30) THEN
!!    FTXAER(5)=XVALUE
!!    ENDIF
!!
!!    RETURN
!!    END SUBROUTINE SETNEW

      SUBROUTINE RCOMP1(NRFUN)

      !INCLUDE 'BR00B.COM'

CF    CHARACTER*2  NLCHAR
      SAVE  IFIRST,PLB49,PLB14,PLB08,PTOPO3,DLAT14
      CHARACTER*80 EPSTAG,TITLE
CF    CHARACTER*15 LLPATH/'/rad4/gcm00b/'/
CF    CHARACTER*40 RFILE1/'rad00b.radfile1.sgpgxg.table8           '/
CF    CHARACTER*40 RFILE2/'rad00b.radfile2.kdist33.tautables4      '/
CF    CHARACTER*40 RFILE3/'rad00b.radfile3.miescatpar.abcdv        '/
CF    CHARACTER*40 RFILE4/'rad00b.radfile4.Ozone.Data.Pratherxx    '/
CF    CHARACTER*40 RFILEA/'rad00b.radfileA.o3trend.1951-2050       '/
CF    CHARACTER*40 RFILEB/'rad00b.radfileB.o3WangJacob.1890.1979   '/
CF    CHARACTER*40 RFILE5/'rad00b.radfile5.trop8aer.tau5090        '/
CF    CHARACTER*40 RFILE6/'rad00b.radfile6.dust8.tau9x8x13         '/
CF    CHARACTER*40 RFILE7/'rad00b.radfile7.stratvol.1950-2000      '/
CF    CHARACTER*40 RFILE8/'rad00b.radfile8.cloud.epsilon4.72x46    '/
CF    CHARACTER*40 RFILE9/'rad00b.radfile9.solar.lean99.uvflux     '/
CF    CHARACTER*40 RFILEC/'rad00b.radfileC.V72X46.1.cor            '/
CF    CHARACTER*40 RFILED/'rad00b.radfileD.Z72X46N                 '/
CF    CHARACTER*40 RFILEE/'rad00b.radfileE.topcld.trscat8          '/
CF    CHARACTER*40 RFILEF/'rad00b.radfileF.ozone4AB.O3CLIM.forNL   '/
CF    CHARACTER*40 BLNK40/'                                        '/
CF    CHARACTER*80 RFILEN

      REAL*4 OZONLJ(44,46),R72X46(72,46),VTAUR4(1800,24)

C     ------------------------------------------------------------------
C     Solar,GHG Trend, VolcAer Size Selection Parameters:    Defaults
C                                           Process       KYEARX  KJDAYX
c                                         SolarCon, UV       0       0
c                                         GH Gas Trend       0       0
c                                                         REFF0= 0.3
c                                                         VEFF0= 0.35
C     ------------------------------------------------------------------

      REAL*8 OZONNL(LX)
      DIMENSION EJMLAT(47),E20LAT(20)
      DIMENSION SLAT14(14),SLAT46(46),WJLO3(7),FTOPO3(19),ITOPO3(19)
      DIMENSION PLB50(50),OZON50(50),OZONXX(50),TTOPO3(46)
      DIMENSION PLB49(49),PLB14(14),PLB08(8),PTOPO3(25),DLAT14(15)
      DATA PLB49/984.,934.,854.,720.,550.,390.,285.,210.,150.,125.,100.,
     +            80.,60.,55.,50.,45.,40.,35.,30.,25.,20.,15.,10.0,7.00,
     +            5.000,4.000,3.000,2.000,1.500,1.000,0.700,0.500,0.400,
     +            0.300,0.200,0.150,0.100,0.070,0.050,0.040,0.030,0.020,
     +            0.015,0.010,0.007,0.005,0.004,0.003,0.0000001/
      DATA PLB14/984.,934.,854.,720.,550.,390.,285.,210.,150.,125.,100.,
     +            80.,60.,55./
      DATA PLB08/984.,934.,854.,720.,550.,390.,255.,150./
      DATA PTOPO3/5.000,4.000,3.000,2.000,1.500,1.000,0.700,0.500,0.400,
     +            0.300,0.200,0.150,0.100,0.070,0.050,0.040,0.030,0.020,
     +            0.015,0.010,0.007,0.005,0.004,0.003,0.0000001/
      DATA DLAT14/-90.0,-60.0,-50.0,-40.0,-30.0,-20.0,-10.0,0.0,
     +             10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 90.0/


c     NRFUN is now set as an argument from calling routine so that unit
c     numbers can be set automatically
      DIMENSION NRFUN(14)
C          radfile1   2   3   4   5   6   7   8   9   A   B   C   D   E
c     DATA NRFUN/71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84/

C       LPATH1 = Character Length of Path Leader of Input Data File Path
C       ----------------------------------------------------------------

CF    DATA LPATH1/13/
      DATA IFIRST/1/

      IF(IFIRST.LT.1) GO TO 9999

C     ------------------------------------------------------------------
C     Input data are read as specified in the first CALL RCOMP1(MADVEL).
C     Subsequent calls to RCOMP1 can be used to re-initialize parameters
C     in SETXXX subroutines to different values, but no new data is read
C     ------------------------------------------------------------------

CF    RFILEN( 1:40)=BLNK40
CF    RFILEN(41:80)=BLNK40
CF    LP3=LPATH1+1
CF    LP4=LPATH1+40

C     ------------------------------------------------------------------
C     MADVEL  Model Add-on Data of Extended Climatology Enable Parameter
C             Each MADVEL digit is ON/OFF switch for corresponding input
C             e.g. MADVEL=123456   (zero digit skips input process)
C
C     MADO3M   =  1   Reads  Makiko's 1951-1997 Ozone climatology RFILEA
C                            Wang-Jacob 1890-1979 Trop O3 change  RFILEB
C     MADAER   =  2   Reads  Aerosol 50y tropospheric climatology RFILE5
C     MADDST   =  3   Reads  Dust-windblown mineral climatology   RFILE6
C     MADVOL   =  4   Reads  Volcanic 1950-00 aerosol climatology RFILE7
C     MADEPS   =  5   Reads  Epsilon cloud heterogeneity data     RFILE8
C     MADLUV   =  6   Reads  Lean's SolarUV 1882-1998 variability RFILE9
C
C                 Related Model Add-on Data Parameters set in BLOCK DATA
C
C     MADGHG   =  1  Default Enables UPDGHG update. (MADGHG=0),no update
C     MADSUR   =  1   Reads  V72X46N.1.cor Vegetation type data   RFILEC
C                            Z72X46N Ocean fraction, topography   RFILED
C     ------------------------------------------------------------------

!!    MADLUV=MADVEL                    !     Reset MAD Enable Parameters
!!    NORMAL=100000
!!    MADO3M=MADLUV/NORMAL                     !       Ozone Climatology
!!    MADLUV=MADLUV-MADO3M*NORMAL
!!    NORMAL=NORMAL/10
!!    MADAER=MADLUV/NORMAL                     !    Trop Aer Climatology
!!    MADLUV=MADLUV-MADAER*NORMAL
!!    NORMAL=NORMAL/10
!!    MADDST=MADLUV/NORMAL                     !        Dust Climatology
!!    MADLUV=MADLUV-MADDST*NORMAL
!!    NORMAL=NORMAL/10
!!    MADVOL=MADLUV/NORMAL                     !    Volc Aer Climatology
!!    MADLUV=MADLUV-MADVOL*NORMAL
!!    NORMAL=NORMAL/10
!!    MADEPS=MADLUV/NORMAL                     !     Cld Eps Climatology
!!    MADLUV=MADLUV-MADEPS*NORMAL              !    Solar UV Climatology

      IF(MADO3M.EQ.0) MOZONE=0     ! no O3 options, if no O3 climatology

C              Initialize variables that might not otherwise get defined
C              ---------------------------------------------------------

      IF(NLP.GT.LX)   call stop_model('increase LX',255)
CF    IF(LASTVC.GT.0) CALL SETATM
      DO 110 L=1,NLP
      TAUWC(L)=0.D0
      TAUIC(L)=0.D0
      SIZEWC(L)=0.D0
      SIZEIC(L)=0.D0
      CLDEPS(L)=0.D0
      FPXCO2(L)=1.D0
      FPXOZO(L)=1.D0
      TLB(L)=250.D0
      TLT(L)=250.D0
      TLM(L)=250.D0
      SHL(L)=0.D0
      RHL(L)=0.D0
      DO 101 I=1,6
      SRAEXT(L,I)=0.D0
      SRASCT(L,I)=0.D0
      SRAGCB(L,I)=0.D0
!obso SRAPI0(L,I)=0.D0
      SRBEXT(L,I)=0.D0
      SRBSCT(L,I)=0.D0
      SRBGCB(L,I)=0.D0
!obso SRBPI0(L,I)=0.D0
      SRDEXT(L,I)=0.D0
      SRDSCT(L,I)=0.D0
      SRDGCB(L,I)=0.D0
!obso SRDPI0(L,I)=0.D0
      SRVEXT(L,I)=0.D0
      SRVSCT(L,I)=0.D0
      SRVGCB(L,I)=0.D0
!obso SRVPI0(L,I)=0.D0
      DBLEXT(L,I)=0.D0
      DBLSCT(L,I)=0.D0
      DBLGCB(L,I)=0.D0
      DBLPI0(L,I)=0.D0
      SRCEXT(L,I)=0.D0
      SRCSCT(L,I)=0.D0
      SRCGCB(L,I)=0.D0
      SRCPI0(L,I)=0.D0
  101 CONTINUE
      DO 102 I=1,33
      TRAALK(L,I)=0.D0
      TRBALK(L,I)=0.D0
      TRDALK(L,I)=0.D0
      TRVALK(L,I)=0.D0
      TRCALK(L,I)=0.D0
      TRGXLK(L,I)=0.D0
  102 CONTINUE
      DO 103 I=1,12
      U0GAS(L,I)=0.D0
      ULGAS(L,I)=0.D0
  103 CONTINUE
      DO 104 I=1,8
      TRACER(L,I)=0.D0
  104 CONTINUE
  110 CONTINUE

      DO 120 L=1,NLP
      PLB0(L)=PLB(L)
      CALL PHATMO(PLB0(L),HLB0(L),D,TLB(L),O,Q,S,OCM,WCM,1,2)
  120 CONTINUE
      DO 121 L=1,NL
      TLT(L)=TLB(L+1)
      TLM(L)=0.5D0*(TLB(L)+TLT(L))
  121 CONTINUE
c
c     Set Undefined Variables to Zero
c
!obso TRSLCR=0.0
!!    DO 136 I=1,33
!!       TAUSL(I)=0.0
!!       FTAUSL(I)=0.0
!!136 CONTINUE
!     DO 138 L=1,LX
!setatm  PSIG(L)=0.0
! 138 CONTINUE

!obso DO 144 L=1,1800
!obso DO 144 I=1,27
!obso    TAUARC(I,L)=0.0
! 144 CONTINUE

C-----------------------------------------------------------------------
CR(1) Reads GTAU Asymmetry Parameter Conversion Table used within SGPGXG
C
C       (SGPGXG does Multiple Scattering Parameterization used in SOLAR)
C       ----------------------------------------------------------------

      NRFU=NRFUN(1)
CF    RFILEN(  1:LPATH1)=LLPATH
CF    RFILEN(LP3:LP4)=RFILE1
CF    OPEN (NRFU,FILE=RFILEN,FORM='UNFORMATTED',STATUS='OLD')
cg    OPEN (NRFU,FORM='UNFORMATTED',STATUS='OLD')             ! CF
      READ (NRFU) GTAU,TGDATA
cg    CLOSE(NRFU)
      IF(KANORM.gt.1.or.KCNORM.gt.1.or.KZSNOW.gt.0) CALL SETGTS


C-----------------------------------------------------------------------
CR(2)    Reads in Merged k-Distribution Tau Tables for Thermal Radiation
C        CFCs, H2O Continuum Tau Table, Merged k-Distr Planck Flux Table
C
C        (Reads: TAUTBL,TAUWV0,PLANCK,XKCFC,H2OCN8,H2OCF8
c                DUCH4,SDUCH4,DUN2O,SDUN2O,ULOX,DUX      used in TAUGAS)
C       ----------------------------------------------------------------

      NRFU=NRFUN(2)
CF    RFILEN(  1:LPATH1)=LLPATH
CF    RFILEN(LP3:LP4)=RFILE2
CF    OPEN (NRFU,FILE=RFILEN,FORM='UNFORMATTED',STATUS='OLD')
cg    OPEN (NRFU,FORM='UNFORMATTED',STATUS='OLD')             ! CF
      READ(NRFU) TAUTBL,TAUWV0,PLANCK,XKCFC,H2OCN8,H2OCF8
      READ(NRFU) DUCH4,SDUCH4,DUN2O,SDUN2O,ULOX,DUX
cg    CLOSE(NRFU)


C        Define Window Flux to Brightness Temperature Conversion Factors
C        ---------------------------------------------------------------

      WAVNA=850.D0
      WAVNB=900.D0
      DO I=1,630
      PFWI=0.001D0*I
      IF(I.GT.100) PFWI=(0.1D0+0.01D0*(I-100))
      IF(I.GT.190) PFWI=(1.0D0+0.10D0*(I-190))
      TKPFW(I)=TKOFPF(WAVNA,WAVNB,PFWI)
      END DO
      DO I=1,900
      PFWI=I
      TKPFT(I)=TKOFPF(0.D0,10000.D0,PFWI)
      END DO

C                            PLANCK Table interpolation limit parameters
C                            -------------------------------------------
      ITPFT0=123
      ITNEXT=250
C-----------------------------------------------------------------------
CR(3)        Read Mie Scattering Parameters [Qext, Qscat, AsymParameter]
C            (1) Tropospheric Aerosols [11 Background, 8 Trop8 Aerosols]
C            (2) Clouds [5 Water, 5 non-spherical Ice, 5 Mie Ice Clouds]
C            (3) Desert Dust Aerosols  [25 particle sizes - to select 8]
C            (4) Volcanic Aerosols [20 particle sizes, 5 size variances]
C            (5) Sulfate  Aerosols [22 particle sizes, 0.1 - 10. micron]
C            (6) Soot   Aerosols [25 particle sizes, 0.001 - 5.0 micron]
C            -----------------------------------------------------------

      NRFU=NRFUN(3)
CF    RFILEN(  1:LPATH1)=LLPATH
CF    RFILEN(LP3:LP4)=RFILE3
CF    OPEN (NRFU,FILE=RFILEN,FORM='FORMATTED',STATUS='OLD')
cg    OPEN (NRFU,FORM='FORMATTED',STATUS='OLD')               ! CF

C                               GCM 11 background aerosol Mie parameters
C                               ----------------------------------------
      DO 301 N=1,11
      READ (NRFU,3000) TITLE
 3000 FORMAT(A80)
      READ (NRFU,3001) (SRAQEX(K,N),K=1,6)
 3001 FORMAT( 18X,6(F7.5,1X))
      READ (NRFU,3001) (SRAQSC(K,N),K=1,6)
      READ (NRFU,3001) (SRAQCB(K,N),K=1,6)
  301 CONTINUE
      READ (NRFU,3002) (Q55A11(N),N=1,11)
 3002 FORMAT( 18X,6(F7.5,1X)/18X,6(F7.5,1X))
      READ (NRFU,3003) (REFA11(N),N=1,11)
 3003 FORMAT( 18X,6(F7.3,1X)/18X,5(F7.3,1X))
      READ (NRFU,3003) (VEFA11(N),N=1,11)
      DO 302 N=1,11
      READ (NRFU,3000) TITLE
      READ (NRFU,3004) (TRAQEX(K,N),K=1,33)
 3004 FORMAT( 14X,7(F7.5,1X),4(/14X,7(F7.5,1X)))
 3005 FORMAT(/14X,7(F7.5,1X),4(/14X,7(F7.5,1X)))
      READ (NRFU,3005) (TRAQSC(K,N),K=1,33)
      READ (NRFU,3005) (TRAQCB(K,N),K=1,33)
  302 CONTINUE

C                       GCM 9 (of 10) climatology aerosol Mie parameters
C                       ------------------------------------------------
      DO 303 N=1,10
      IF(N.EQ.6) GO TO 303
      READ (NRFU,3000) TITLE
      READ (NRFU,3001) (SRBQEX(K,N),K=1,6)
      READ (NRFU,3001) (SRBQSC(K,N),K=1,6)
      READ (NRFU,3001) (SRBQCB(K,N),K=1,6)
  303 CONTINUE
      READ (NRFU,3002) (Q55B10(N),N=1,5),(Q55B10(N),N=7,10)
      READ (NRFU,3003) (REFB10(N),N=1,5),(REFB10(N),N=7,10)
      READ (NRFU,3003) (VEFB10(N),N=1,5),(VEFB10(N),N=7,10)
      DO 304 N=1,10
      IF(N.EQ.6) GO TO 304
      READ (NRFU,3000) TITLE
      READ (NRFU,3004) (TRBQEX(K,N),K=1,33)
      READ (NRFU,3005) (TRBQSC(K,N),K=1,33)
      READ (NRFU,3005) (TRBQCB(K,N),K=1,33)
  304 CONTINUE


C                               Cloud Water, Ice-non, Ice-Mie parameters
C                               ----------------------------------------
      DO 305 N=1,15
      READ (NRFU,3000) TITLE
      READ (NRFU,3001) (SRCQEX(K,N),K=1,6)
      READ (NRFU,3001) (SRCQSC(K,N),K=1,6)
      READ (NRFU,3001) (SRCQCB(K,N),K=1,6)
  305 CONTINUE
      READ (NRFU,3006) (Q55C15(N),N=1,15)
 3006 FORMAT( 18X,6(F7.5,1X)/18X,6(F7.5,1X)/18X,6(F7.5,1X))
      READ (NRFU,3007) (REFC15(N),N=1,15)
 3007 FORMAT( 18X,6(F7.3,1X)/18X,6(F7.3,1X)/18X,6(F7.3,1X))
      READ (NRFU,3007) (VEFC15(N),N=1,15)
      DO 306 N=1,15
      READ (NRFU,3000) TITLE
      READ (NRFU,3004) (TRCQEX(K,N),K=1,33)
      READ (NRFU,3005) (TRCQSC(K,N),K=1,33)
      READ (NRFU,3005) (TRCQCB(K,N),K=1,33)
      READ (NRFU,3005) (TRCQAL(K,N),K=1,33)
  306 CONTINUE

C                               Desert Dust 25 sizes, Mie parameter data
C                               ----------------------------------------
      DO 307 N=1,25
      READ (NRFU,3001) (SRDQEX(K,N),K=1,6)
      READ (NRFU,3001) (SRDQSC(K,N),K=1,6)
      READ (NRFU,3001) (SRDQCB(K,N),K=1,6)
  307 CONTINUE
      READ (NRFU,3008) (Q55D25(N),N=1,25)
 3008 FORMAT( 18X,5(F7.5,1X),4(/18X,5(F7.5,1X)))
      READ (NRFU,3009) (REFD25(N),N=1,25)
 3009 FORMAT( 18X,12(F3.1,1X)/18X,12(F3.1,1X)/18X,F3.0)
      READ (NRFU,3010) (VEFD25(N),N=1,25)
 3010 FORMAT( 18X,12(F3.1,1X)/18X,12(F3.1,1X)/18X,F3.1)
      DO 308 N=1,25
      READ (NRFU,3000) TITLE
      READ (NRFU,3004) (TRDQEX(K,N),K=1,33)
      READ (NRFU,3005) (TRDQSC(K,N),K=1,33)
      READ (NRFU,3005) (TRDQCB(K,N),K=1,33)
      READ (NRFU,3005) (TRDQAL(K,N),K=1,33)
  308 CONTINUE

C                               Volcanic aerosol Mie size, variance data
C                               ----------------------------------------
      DO 313 M=1,5
      IF(M.EQ.4) GO TO 313
      DO 311 N=1,20
      READ (NRFU,3001) (SRVQEX(K,N,M),K=1,6)
      READ (NRFU,3001) (SRVQSC(K,N,M),K=1,6)
      READ (NRFU,3001) (SRVQCB(K,N,M),K=1,6)
  311 CONTINUE
      READ (NRFU,3011) (Q55V20(N,M),N=1,20)
 3011 FORMAT( 18X,5(F7.5,1X),3(/18X,5(F7.5,1X)))
      READ (NRFU,3012) (REFV20(N,M),N=1,20)
 3012 FORMAT( 18X,12(F3.1,1X)/18X,8(F3.1,1X))
      READ (NRFU,3012) (VEFV20(N,M),N=1,20)
      DO 312 N=1,20
      READ (NRFU,3000) TITLE
      READ (NRFU,3004) (TRVQEX(K,N,M),K=1,33)
      READ (NRFU,3005) (TRVQSC(K,N,M),K=1,33)
      READ (NRFU,3005) (TRVQCB(K,N,M),K=1,33)
      READ (NRFU,3005) (TRVQAL(K,N,M),K=1,33)
  312 CONTINUE
  313 CONTINUE
      DO 316 N=1,20
      DO 314 K=1,6
      SRVQEX(K,N,4)=(SRVQEX(K,N,3)+SRVQEX(K,N,5))/2.D0
      SRVQSC(K,N,4)=(SRVQSC(K,N,3)+SRVQSC(K,N,5))/2.D0
      SRVQCB(K,N,4)=(SRVQCB(K,N,3)+SRVQCB(K,N,5))/2.D0
  314 CONTINUE
      Q55V20(N,4)=(Q55V20(N,3)+Q55V20(N,5))/2.D0
      REFV20(N,4)=(REFV20(N,3)+REFV20(N,5))/2.D0
      VEFV20(N,4)=(VEFV20(N,3)+VEFV20(N,5))/2.D0
      DO 315 K=1,33
      TRVQEX(K,N,4)=(TRVQEX(K,N,3)+TRVQEX(K,N,5))/2.D0
      TRVQSC(K,N,4)=(TRVQSC(K,N,3)+TRVQSC(K,N,5))/2.D0
      TRVQCB(K,N,4)=(TRVQCB(K,N,3)+TRVQCB(K,N,5))/2.D0
      TRVQAL(K,N,4)=(TRVQAL(K,N,3)+TRVQAL(K,N,5))/2.D0
  315 CONTINUE
  316 CONTINUE

C                            Sulfate aerosol, Mie parameter 22-size data
C                            -------------------------------------------
      DO 321 N=1,22
      READ (NRFU,3000) TITLE
      READ (NRFU,3001) (SRUQEX(K,N),K=1,6)
      READ (NRFU,3001) (SRUQSC(K,N),K=1,6)
      READ (NRFU,3001) (SRUQCB(K,N),K=1,6)
  321 CONTINUE
      READ (NRFU,3008) (Q55U22(N),N=1,22)
      READ (NRFU,3013) (REFU22(N),N=1,22)
 3013 FORMAT( 18X,5(F7.3,1X),4(/18X,5(F7.3,1X)))
      READ (NRFU,3013) (VEFU22(N),N=1,22)
      DO 322 N=1,22
      READ (NRFU,3000) TITLE
      READ (NRFU,3004) (TRUQEX(K,N),K=1,33)
      READ (NRFU,3005) (TRUQSC(K,N),K=1,33)
      READ (NRFU,3005) (TRUQCB(K,N),K=1,33)
      READ (NRFU,3005) (TRUQAL(K,N),K=1,33)
  322 CONTINUE

C                               Soot aerosol, Mie parameter 25-size data
C                               ----------------------------------------
      DO 323 N=1,25
      READ (NRFU,3000) TITLE
      READ (NRFU,3001) (SRSQEX(K,N),K=1,6)
      READ (NRFU,3001) (SRSQSC(K,N),K=1,6)
      READ (NRFU,3001) (SRSQCB(K,N),K=1,6)
  323 CONTINUE
      READ (NRFU,3008) (Q55S25(N),N=1,25)
      READ (NRFU,3013) (REFS25(N),N=1,25)
      READ (NRFU,3013) (VEFS25(N),N=1,25)
      DO 324 N=1,25
      READ (NRFU,3000) TITLE
      READ (NRFU,3004) (TRSQEX(K,N),K=1,33)
      READ (NRFU,3005) (TRSQSC(K,N),K=1,33)
      READ (NRFU,3005) (TRSQCB(K,N),K=1,33)
      READ (NRFU,3005) (TRSQAL(K,N),K=1,33)
  324 CONTINUE

C                            Seasalt aerosol, Mie parameter 22-size data
C                            Nitrate aerosol, Mie parameter 22-size data
C                            (Water) aerosol, Mie parameter 22-size data
C                            -------------------------------------------
      N1=23
      DO 326 KK=1,4
      N2=N1+21
      DO 325 N=N1,N2
      READ (NRFU,3000) TITLE
      READ (NRFU,3001) (SRUQEX(K,N),K=1,6)
      READ (NRFU,3001) (SRUQSC(K,N),K=1,6)
      READ (NRFU,3001) (SRUQCB(K,N),K=1,6)
  325 CONTINUE
      READ (NRFU,3008) (Q55U22(N),N=N1,N2)
      READ (NRFU,3013) (REFU22(N),N=N1,N2)
      READ (NRFU,3013) (VEFU22(N),N=N1,N2)
      N1=N2+1
  326 CONTINUE
      N1=23
      DO 328 KK=1,4
      N2=N1+21
      DO 327 N=N1,N2
      READ (NRFU,3000) TITLE
      READ (NRFU,3004) (TRUQEX(K,N),K=1,33)
      READ (NRFU,3005) (TRUQSC(K,N),K=1,33)
      READ (NRFU,3005) (TRUQCB(K,N),K=1,33)
      READ (NRFU,3005) (TRUQAL(K,N),K=1,33)
  327 CONTINUE
      N1=N2+1
  328 CONTINUE

C         Mie scattering Qext,Cosb tables (Nr=1.30-1.60),(Reff=0.0-15.0)
C         used for hygroscopic aerosol relative humidity size dependence
C         --------------------------------------------------------------

      READ (NRFU,3014) (XNR(J),J=1,7)
 3014 FORMAT(13X,F5.3,12F9.3)
      DO 331 I=1,183
      READ (NRFU,3015) J,SIZENR(I),(QEXTNR(I,J),J=1,7)
 3015 FORMAT(I3,F6.2,13F9.5)
  331 CONTINUE
      READ (NRFU,3014) (XNR(J),J=8,13)
      DO 332 I=1,183
      READ (NRFU,3015) J,SIZENR(I),(QEXTNR(I,J),J=8,13)
  332 CONTINUE
      READ (NRFU,3014) (XNR(J),J=1,7)
      DO 333 I=1,183
      READ (NRFU,3015) J,SIZENR(I),(COSBNR(I,J),J=1,7)
  333 CONTINUE
      READ (NRFU,3014) (XNR(J),J=8,13)
      DO 334 I=1,183
      READ (NRFU,3015) J,SIZENR(I),(COSBNR(I,J),J=8,13)
  334 CONTINUE

  399 CONTINUE
cg    CLOSE(NRFU)

C-----------------------------------------------------------------------
CR(4)   Reads Prather (1988) Ozone Latitude and Height Distribution Data
C            + London (1976) Ozone Climatology with Longitude Dependence
C            +Keating (1985) Upper Stratosphere Climatology O3 above 5mb
C       ----------------------------------------------------------------

      IF(MOZONE.EQ.0) GO TO 470
      NRFU=NRFUN(4)
CF    RFILEN(  1:LPATH1)=LLPATH
CF    RFILEN(LP3:LP4)=RFILE4
CF    OPEN (NRFU,FILE=RFILEN,FORM='FORMATTED',STATUS='OLD')
cg    OPEN (NRFU,FORM='FORMATTED',STATUS='OLD')               ! CF

      DO 410 I=1,27
      READ (NRFU,5401) TITLE
 5401 FORMAT(A80)
  410 CONTINUE
      DO 414 M=1,12
      DO 411 I=1,4
      READ (NRFU,5401) TITLE
  411 CONTINUE
      DO 412 L=1,31
      READ (NRFU,5402) PO3(L),(O3DPPM(L,J,M),J=1,8)
 5402 FORMAT(3X,F9.4,4X,8F8.3)
  412 CONTINUE
      READ (NRFU,5401) TITLE
      READ (NRFU,5401) TITLE
      DO 413 L=1,31
      READ (NRFU,5402) PO3(L),(O3DPPM(L,J,M),J=9,16)
  413 CONTINUE
  414 CONTINUE
      DO 415 I=1,13
      READ (NRFU,5401) TITLE
  415 CONTINUE
C                                        London (1976) Ozone Climatology
C                                        -------------------------------
      DO 417 K=1,18
      READ (NRFU,5401) TITLE
      DO 416 J=1,18
      READ (NRFU,5403) (O3AVE(I,J,K),I=1,12)
 5403 FORMAT(11X,12(F4.3,1X))
  416 CONTINUE
  417 CONTINUE
      READ (NRFU,5401) TITLE
      READ (NRFU,5404) (AO3AV1(I),I=1,216)
 5404 FORMAT(30X,7(F5.4,1X),19(/6X,11(F5.4,1X)))


      DO 419 K=1,2
      READ (NRFU,5401) TITLE
      DO 418 J=1,19
      READ (NRFU,5405) (SO3JFS(I,J,K),I=1,11)
 5405 FORMAT(11X,11(F4.1,1X))
  418 CONTINUE
  419 CONTINUE

      DO 420 J=2,19
      E20LAT(J)=-90.D0+(J-1.5D0)*180.D0/18.D0
  420 CONTINUE
      E20LAT( 1)=-90.D0
      E20LAT(20)= 90.D0
      DO 421 J=2,46
      EJMLAT(J)=-90.D0+(J-1.5D0)*180.D0/45.D0
  421 CONTINUE
      EJMLAT( 1)=-90.D0
      EJMLAT(47)= 90.D0
C                              -----------------------------------------
C                              Keating (1985) Upper Strat O3 Climatology
C                              (TOPO3 in ppmv at 25 Pressures above 5mb)
C                              (Tabulated at 19 Lat, Interpolated to 46)
C                              -----------------------------------------

      READ (NRFU,5401) TITLE
      READ (NRFU,5401) TITLE
      DO 427 M=1,12
      READ (NRFU,5401) TITLE
      READ (NRFU,5401) TITLE
      DO 424 LL=1,25
      L=26-LL
      READ (NRFU,5406) (ITOPO3(J),J=1,19)
 5406 FORMAT(4X,19I4)
      DO 422 J=1,19
      FTOPO3(J)=ITOPO3(J)/100.D0
  422 CONTINUE
      CALL RETERP(FTOPO3,E20LAT,20,TTOPO3,EJMLAT,47)
      DO 423 J=1,46
      TOPO3(L,J,M)=TTOPO3(J)
  423 CONTINUE
  424 CONTINUE

C                  -----------------------------------------------------
C                  Scale Height of Homogeneous Atm at 0 deg C = 7.991 km
C                  Conversion of PPM to ACM: ACM per PPM = 7.991D5*1.D-6
C            Thus:
C                  TOPO3 is in ACM: (1 ACM = 1000 DU = 2.687E+19mol/cm2)
C                  -----------------------------------------------------

      DO 426 J=1,46
      DO 425 L=1,24
      ACML=7.991D+05*1.D-06*(PTOPO3(L)-PTOPO3(L+1))/1013.25D0
      TOPO3(L,J,M)=ACML*0.5D0*(TOPO3(L,J,M)+TOPO3(L+1,J,M))
  425 CONTINUE
      TOPO3(25,J,M)=0.D0
  426 CONTINUE
  427 CONTINUE
      DO 431 M=1,12
      READ (NRFU,5401) TITLE
      READ (NRFU,5401) TITLE
      READ (NRFU,5401) TITLE
      DO 430 LL=1,50
      L=51-LL
      READ (NRFU,5407) (SAGEZP(L,J,M),J=1,12)
 5407 FORMAT(8X,12F10.5)
      SLAT14( 1)=3.D0*SAGEZP(L, 1,M)-2.D0*SAGEZP(L, 2,M)
      SLAT14(14)=3.D0*SAGEZP(L,12,M)-2.D0*SAGEZP(L,11,M)
      DO 428 J=1,12
      SLAT14(J+1)=SAGEZP(L,J,M)
  428 CONTINUE
      CALL RETERP(SLAT14,DLAT14,15,SLAT46,EJMLAT,47)
      DO 429 J=1,46
      SZP50(L,J,M)=SLAT46(J)
  429 CONTINUE
  430 CONTINUE
  431 CONTINUE
cg    CLOSE(NRFU)

C-----------------------------------------------------------------------
CR(A)   Reads Makiko's Jan 1951 - Dec 1997 Zonal Ozone Trend Climatology
C             O3CLIM(M,L,J)  M=Month(1-1200), L=Layer(1-44), J=Lat(1-46)
C             Ozone Vertical Distribution is Re-Partitioned to NL Layers
C             and Ozone in DU is converted Ozone in CMA (1 CMA = 1000DU)
C             ----------------------------------------------------------

      IF(MADO3M.LT.1) GO TO 499
      NRFU=NRFUN(10)
CF    RFILEN(  1:LPATH1)=LLPATH
CF    RFILEN(LP3:LP4)=RFILEA
CF    OPEN (NRFU,FILE=RFILEN,FORM='UNFORMATTED',STATUS='OLD')
cg    OPEN (NRFU,FORM='UNFORMATTED',STATUS='OLD')             ! CF

      NLO3=45
      DO 432 L=1,NLO3
      IF(L.LT.15) PLB50(L)=PLB14(L)
      IF(L.GT.14) PLB50(L)=1382.D0*EXP(-(L+5.5D0)/6.2D0)
  432 CONTINUE
      NL1=NL+1
      DO M=1,MO3X
      MM=M-((M-1)/12)*12
      READ (NRFU,END=440) TITLE,OZONLJ
      DO 439 J=1,46
      DO 434 L=1,44
      OZON50(L)=OZONLJ(L,J)*1.D-03
  434 CONTINUE

      IF(MOZONE.GT.4) THEN
C        ---------------------------------------------------------------
C        1.  Replace exponential pressure profile with SAGE2 SZP profile
C            and repartition the SAGE ozone to standard pressure profile
C        2.  Linear transition from tropospheric ozone (below 5.0 mb) to
C            the stratospheric (Keating 1985) ozone profile above 1.0 mb
C            unless MOZONE=6
C        3.  Add (Keating 1985) stratospheric ozone profile above 1.0 mb
C        ---------------------------------------------------------------

      NLO3=49
      DO 435 L=1,31
      PLB50(L+14)=0.5D0*(SZP50(L,J,MM)+SZP50(L+1,J,MM))
  435 CONTINUE
      CALL REPART(OZON50,PLB50,45,OZONXX,PLB49,NLO3)
      DO 436 L=1,NLO3
      PLB50(L)=PLB49(L)
      IF(L.LE.29) OZON50(L)=OZONXX(L)
      IF(L.GT.29) OZON50(L)=TOPO3(L-24,J,MM)
  436 CONTINUE
      END IF
      IF(MOZONE.eq.5) THEN
      DO 437 LL=1,4      ! transition L=26-29
      WBOT=(5-LL)*0.2D0
      WTOP=LL*0.2D0
      OZON50(25+LL)=WBOT*OZONXX(25+LL)+WTOP*TOPO3(LL+1,J,MM)
  437 CONTINUE
      ENDIF
      CALL REPART(OZON50,PLB50,NLO3,OZONNL,PLB,NL1)
      DO 438 L=1,NL
      O3CLIM(M,L,J)=OZONNL(L)
  438 CONTINUE
  439 CONTINUE
      END DO

  440 CONTINUE
cg    CLOSE(NRFU)

C-----------------------------------------------------------------------
CR(B)   Reads Wang-Jacob Tropospheric Ozone: 1890 and 1979 distributions
C             WJ1890(I,J,L,M) Pre-Industrial (L=1,7) 72x46 Monthly Means
C             WJ1979(I,J,L,M) Mid-Industrial (L=1,7) 72x46 Monthly Means
C             Exponentially interpolated longitudinal trend is extracted
C             and superimposed on Makiko's Zonal Ozone Trend Climatology
C             ----------------------------------------------------------

      NRFU=NRFUN(11)
CF    RFILEN(  1:LPATH1)=LLPATH
CF    RFILEN(LP3:LP4)=RFILEB
CF    OPEN (NRFU,FILE=RFILEN,FORM='UNFORMATTED',STATUS='OLD')
cg    OPEN (NRFU,FORM='UNFORMATTED',STATUS='OLD')             ! CF

      DO 444 M=1,12
      DO 443 L=1,7
      READ (NRFU) TITLE,R72X46
      DO 442 J=1,46
      DO 441 I=1,72
      WJ1890(I,J,L,M)=R72X46(I,J)
  441 CONTINUE
  442 CONTINUE
  443 CONTINUE
  444 CONTINUE
      DO 448 M=1,12
      DO 447 L=1,7
      READ (NRFU) TITLE,R72X46
      DO 446 J=1,46
      DO 445 I=1,72
      WJ1979(I,J,L,M)=R72X46(I,J)
  445 CONTINUE
  446 CONTINUE
  447 CONTINUE
  448 CONTINUE

C             Define ozone data polar points as J=2 and J=45 mean values
C             ----------------------------------------------------------

      DO 460 M=1,12
      DO 459 L=1,7
      SUM=0.D0
      DO 451 I=1,72
      SUM=SUM+WJ1890(I,2,L,M)
  451 CONTINUE
      SUM=SUM/72.D0
      DO 452 I=1,72
      WJ1890(I,1,L,M)=SUM
  452 CONTINUE
      SUM=0.D0
      DO 453 I=1,72
      SUM=SUM+WJ1979(I,2,L,M)
  453 CONTINUE
      SUM=SUM/72.D0
      DO 454 I=1,72
      WJ1979(I,1,L,M)=SUM
  454 CONTINUE
      SUM=0.D0
      DO 455 I=1,72
      SUM=SUM+WJ1890(I,45,L,M)
  455 CONTINUE
      SUM=SUM/72.D0
      DO 456 I=1,72
      WJ1890(I,46,L,M)=SUM
  456 CONTINUE
      SUM=0.D0
      DO 457 I=1,72
      SUM=SUM+WJ1979(I,45,L,M)
  457 CONTINUE
      SUM=SUM/72.D0
      DO 458 I=1,72
      WJ1979(I,46,L,M)=SUM
  458 CONTINUE
  459 CONTINUE
  460 CONTINUE
cg    CLOSE(NRFU)

  470 CONTINUE

C                -------------------------------------------------------
C                Repartition  Wang-Jacob Tropospheric Longitudinal Ozone
C                             from 7 reference layers with top at 150 mb
C                             to currently specified PLB pressure levels
C                             ------------------------------------------
      NL1=NL+1
      DO 477 M=1,12
      DO 476 J=1,46
      DO 475 I=1,72
      DO 471 L=1,7
      ACML=7.991D+05*1.D-09*(PLB08(L)-PLB08(L+1))/1013.25D0
      WJLO3(L)=WJ1890(I,J,L,M)*ACML
  471 CONTINUE
      CALL REPART(WJLO3,PLB08,8,OZONXX,PLB,NL1)
      DO 472 L=1,NL
      WJ1890(I,J,L,M)=OZONXX(L)
  472 CONTINUE
      DO 473 L=1,7
      ACML=7.991D+05*1.D-09*(PLB08(L)-PLB08(L+1))/1013.25D0
      WJLO3(L)=WJ1979(I,J,L,M)*ACML
  473 CONTINUE
      CALL REPART(WJLO3,PLB08,8,OZONXX,PLB,NL1)
      DO 474 L=1,NL
      WJ1979(I,J,L,M)=OZONXX(L)
  474 CONTINUE
  475 CONTINUE
  476 CONTINUE
  477 CONTINUE

  499 CONTINUE

C-----------------------------------------------------------------------
CR(5) TROAER: Monthly-Mean Tropospheric Aerosols  (Column Optical Depth)
C                  Decade Maps (N=1,5) for BCI, OCI, SUI (K=1,3) Aerosol
C                  for K=4 (N=1,5) gives SEA, SUN, OCN, OCB, BCB Aerosol
C                  Map: IJ=72x46,  Month: M=1-12, AnnAv: M=13
C
C     VDBSCU: Zonal Mean Vertical Distribution for BCI,SUI-TYPE Aerosols
C     ------------------------------------------------------------------

      IF(MADAER.LT.1) GO TO 599
      NRFU=NRFUN(5)
CF    RFILEN(  1:LPATH1)=LLPATH
CF    RFILEN(LP3:LP4)=RFILE5
CF    OPEN (NRFU,FILE=RFILEN,FORM='UNFORMATTED',STATUS='OLD')
cg    OPEN (NRFU,FORM='UNFORMATTED',STATUS='OLD')             ! CF
      READ (NRFU) TROAER
      READ (NRFU) VDBCSU
cg    CLOSE(NRFU)

C     ------------------------------------------------------------------
C     VDBCSU Input Data Adjustment:
C
C             Renormalize VDBCSU(46,12,13,3) BC and SU Vertical Profiles
C                  and Recompute Annual-Mean BC and SU Vertical Profiles
C     ------------------------------------------------------------------

      DO 515 N=1,3
      DO 514 M=1,12
      DO 513 J=1,46
      SUM=0.D0
      DO 511 L=1,12
      SUM=SUM+VDBCSU(J,L,M,N)
  511 CONTINUE
      DO 512 L=1,12
      VDBCSU(J,L,M,N)=VDBCSU(J,L,M,N)/SUM
  512 CONTINUE
  513 CONTINUE
  514 CONTINUE
  515 CONTINUE
      DO 519 N=1,3
      DO 518 L=1,12
      DO 517 J=1,46
      SUM=0.D0
      DO 516 M=1,12
      SUM=SUM+VDBCSU(J,L,M,N)
  516 CONTINUE
      VDBCSU(J,L,13,N)=SUM/12.D0
  517 CONTINUE
  518 CONTINUE
  519 CONTINUE
  599 CONTINUE

C-----------------------------------------------------------------------
CR(6) DUST:   Monthly-Mean Desert Dust (Clay,Silt) 8-Size Optical Depths
C                  Map: IJ=72x46,  Lay: L=1-9, Siz: S=1-8, Month: M=1-12
C                  -----------------------------------------------------

      IF(MADDST.LT.1) GO TO 699
      NRFU=NRFUN(6)
CF    RFILEN(  1:LPATH1)=LLPATH
CF    RFILEN(LP3:LP4)=RFILE6
CF    OPEN (NRFU,FILE=RFILEN,FORM='UNFORMATTED',STATUS='OLD')
cg    OPEN (NRFU,FORM='UNFORMATTED',STATUS='OLD')             ! CF
      READ (NRFU) TDUST
cg    CLOSE(NRFU)

  699 CONTINUE

C     ------------------------------------------------------------------
C     IFIRST=1     Read in 2 Files of Tropospheric Aerosol and Dust Data
C                  Apply Multiplicative Scaling Factors: for Aerosol and
C                  and Desert Dust Optical Depths, and Single Scattering
C                  Albedos using (FSAERO,FTAERO,PI0MAX, & FSDUST,FTDUST)
C                  defined in BLOCK DATA RADPAR COMMON/AERCOM (defaults)
C
C        NOTE:     Scaling Factors Multiply the Appropriate Aerosol  MIE
C                  Radiative Parameters (QXAERO,QSAERO, & QXDUST,QSDUST)
C                  and NOT  the Column Optical Depth (72x46) Global Maps
C
C     IFIRST=0     Subsequent CALLs to SETAER  can Reset Scaling Factors
C     ------------------------------------------------------------------


C-----------------------------------------------------------------------
CR(7)        Read Makiko's Stratospheric binary data made in April, 2002
C                               (1800 months (1850-1999) x 24 latitudes)
C                               ----------------------------------------

      IF(MADVOL.LT.1) GO TO 799
      NRFU=NRFUN(7)
CF    RFILEN(  1:LPATH1)=LLPATH
CF    RFILEN(LP3:LP4)=RFILE7
CF    OPEN (NRFU,FILE=RFILEN,FORM='UNFORMATTED',STATUS='OLD')
cg    OPEN (NRFU,FORM='UNFORMATTED',STATUS='OLD')             ! CF
      READ (NRFU) TITLE
      IF(TITLE(1:13).eq.'Optical Depth')
     &     call stop_model('use new RADN7',255)
      REWIND (NRFU)
      DO K=1,5
        READ (NRFU) TITLE,VTAUR4
        DO J=1,24
        DO I=1,1800
          V4TAUR(I,J,K)=VTAUR4(I,J)
        END DO
        END DO
      END DO
cg    CLOSE(NRFU)

  799 CONTINUE

C-----------------------------------------------------------------------
CR(8)  ISCCP Derived Cloud Variance (EPSILON) Cloud Optical Depth Factor
C      Low, Mid, High  Cloud Optical Depths are Reduced by (1 - EPSILON)
C
C               INPUT DATA FILE:  UNIT = INFILE
C                                 TAG  = EPSTAG  (CHARACTER*80)
C                                 DATA = EPLMHC  (72,46,12,4) REAL*4
C
C      Data are 72X46 Monthly Mean Low, Mid, High, Column EPSILON Values
C      Cloud Heterogeneity selections used in UPDEPS, GETEPS (in SETCLD)
C
C             EPSCON  Column Cloud Inhomogeneity EPSILON (when KCLDEP=1)
C             KCLDEP  Selects Cloud Inhomogeneity Option (0-4):
C                     KCLDEP =  0  Sets Column CLDEPS to Zero
C                     KCLDEP =  1  Sets Column CLDEPS to EPSCON
C                     KCLDEP =  2  Keeps whatever is specified in CLDEPS
C                     KCLDEP =  3  Uses: Column EPCOL(72,46) Climatology
C                     KCLDEP =  4  Uses: Ht Dep EPLOW, EPMID, EPHIG Data
C               --------------------------------------------------------

      IF(MADEPS.LT.1) GO TO 899
      NRFU=NRFUN(8)
CF    RFILEN(  1:LPATH1)=LLPATH
CF    RFILEN(LP3:LP4)=RFILE8
CF    OPEN (NRFU,FILE=RFILEN,FORM='UNFORMATTED',STATUS='OLD')
cg    OPEN (NRFU,FORM='UNFORMATTED',STATUS='OLD')             ! CF
      READ (NRFU) EPSTAG,EPLMHC
cg    CLOSE(NRFU)

      DO 811 N=1,4
      DO 810 M=1,12
      DO 809 I=1,72
      L=0
  801 CONTINUE
      L=L+1
      IF(EPLMHC(I,L,M,N).LT.0.0) GO TO 801
      IF(L.GT.1) THEN
      DO 802 J=1,L
      EPLMHC(I,J,M,N)=EPLMHC(I,L,M,N)
  802 CONTINUE
      ENDIF
      L=47
  803 CONTINUE
      L=L-1
      IF(EPLMHC(I,L,M,N).LT.0.0) GO TO 803
      IF(L.LT.46) THEN
      DO 804 J=L,46
      EPLMHC(I,J,M,N)=EPLMHC(I,L,M,N)
  804 CONTINUE
      ENDIF
      J=0
  805 CONTINUE
      J=J+1
      IF(J.GT.46) GO TO 808
      IF(EPLMHC(I,J,M,N).GE.0.0) GO TO 805
      K=J-1
      EPK=EPLMHC(I,K,M,N)
  806 CONTINUE
      J=J+1
      IF(J.GT.46) GO TO 808
      IF(EPLMHC(I,J,M,N).LT.0.0) GO TO 806
      L=J
      EPL=EPLMHC(I,L,M,N)
      NN=L-K-1
      DEP=(EPL-EPK)/(L-K)
      DO 807 L=1,NN
      EPLMHC(I,K+L,M,N)=EPK+L*DEP
  807 CONTINUE
      GO TO 805
  808 CONTINUE
  809 CONTINUE
  810 CONTINUE
  811 CONTINUE

  899 CONTINUE


C-----------------------------------------------------------------------
CR(E)
C             KCLDEM  Selects: Top-Cloud (Thermal) Scattering Correction
C                     KCLDEM =  0  Utilizes Non-scattering approximation
C                     KCLDEM =  1  Modifies emission and transmission by
C                                  top cloud (over-rides old correction)
C             ----------------------------------------------------------

      NRFU=NRFUN(14)
CF    RFILEN(  1:LPATH1)=LLPATH
CF    RFILEN(LP3:LP4)=RFILEE
CF    OPEN (NRFU,FILE=RFILEN,FORM='UNFORMATTED',STATUS='OLD')
cg    OPEN (NRFU,FORM='UNFORMATTED',STATUS='OLD')             ! CF
      READ (NRFU) RIJTPG,FDXTPG,FEMTPG
cg    CLOSE(NRFU)


C-----------------------------------------------------------------------
CR(9)         Read Judith Lean's Solar UV and Solar Constant Variability
C                                      Monthly-Mean Solar UV (1882-1998)
C                                      ---------------------------------


      IF(MADLUV.LT.1) GO TO 909
      NRFU=NRFUN(9)
CF    RFILEN(  1:LPATH1)=LLPATH
CF    RFILEN(LP3:LP4)=RFILE9
CF    OPEN (NRFU,FILE=RFILEN,FORM='FORMATTED',STATUS='OLD')
cg    OPEN (NRFU,FORM='FORMATTED',STATUS='OLD')               ! CF

      READ(NRFU,'(a80)') TITLE
      if(ksolar.ge.2 .and. TITLE(1:3).ne.'ANN')
     &     call stop_model('change RADN9',255)
      if(ksolar.lt.2 .and. TITLE(1:3).eq.'ANN')
     &     call stop_model('change RADN9',255)
      READ(NRFU,'(5F14.2)') (WSLEAN(I),I=1,190)
      READ(NRFU,'(a80)') TITLE
      READ(NRFU,'(5E14.3)') (DSLEAN(I),I=1,190)
      DO I=1,190
        WSLEAN(I)=WSLEAN(I)/1000.D0
        DSLEAN(I)=DSLEAN(I)/1000.D0
        W1LEAN(I)=WSLEAN(I)-0.5D0*DSLEAN(I)
      END DO
      READ(NRFU,'(a80)') TITLE
      READ(NRFU,'(a80)') TITLE
      IF(KSOLAR.LT.2) THEN
C****   Read in monthly-mean data
        DO I=1,Ms0X
          READ(NRFU,'(2I6,3F17.6)') IYEAR,IMONTH,TSI1(I),TSI2(I)
          READ(NRFU,'(5E14.6)')     (FSLEAN(K),K=1,190)
          SUM=0.D0
          DO K=1,190
            SUM=SUM+FSLEAN(K)*DSLEAN(K)
          END DO
          SFNORM=TSI1(I)/SUM
          DO K=1,190
            UVLEAN(I,K)=FSLEAN(K)*SFNORM
          END DO
        END DO
      ELSE
C****   Read in annual-mean data
        DO I=1,Ms0X
          READ(NRFU,'(F12.1,2F15.4)',end=908) yr2S0,TSI1(I),TSI2(I)
          if(I.eq.1) yr1S0 = yr2S0
          READ(NRFU,'(5E14.6)')               (UVLEAN(I,K),K=1,190)
        END DO
  908   write(6,*) 'read S0-history: ',yr1S0,' - ',yr2S0
      END IF
cg    CLOSE(NRFU)
      GO TO 949

  909 CONTINUE
      DO I=1,190
        WSLEAN(I)=WSLEAN(I)/1000.D0
        DSLEAN(I)=DSLEAN(I)/1000.D0
        W1LEAN(I)=WSLEAN(I)-0.5D0*DSLEAN(I)
      END DO

  949 CONTINUE


C-----------------------------------------------------------------------
CR(C)     Read:    Elaine Mathews 10 Fractional Vegetation Distributions
C         10 global maps (72x46) depict fractional vegetation/soil types
C         Map-1 (bright sand) + Map-10 (black dirt) define desert albedo
C         (sum of Maps 1-10 over land-area (ILON,JLAT) grid boxes = 1.0)
C
C         Map-11 refers to plankton concentrations over ocean areas that
C         are yet to be implemented.
C         --------------------------------------------------------------


CF    IF(MADSUR.LT.1) GO TO 999
CF    NRFU=NRFUN(12)
CF    RFILEN(  1:LPATH1)=LLPATH
CF    RFILEN(LP3:LP4)=RFILEC
CF    OPEN (NRFU,FILE=RFILEN,FORM='UNFORMATTED',STATUS='OLD')
CF    OPEN (NRFU,FORM='UNFORMATTED',STATUS='OLD')             ! CF

CF    DO 953 K=1,11
CF    READ(NRFU) TITLE,R72X46
CF    DO 952 J=1,46
CF    DO 951 I=1,72
CF    FVEG11(I,J,K)=R72X46(I,J)
CF951 CONTINUE
CF952 CONTINUE
CF953 CONTINUE


C-----------------------------------------------------------------------
CR(D)      Read:   1   FOCEAN   72x46 ocean fraction   (FOCEAN = 0 or 1)
C                  2   FLAKE    72x46 lake  fraction
C                  3   FGRND    72x46 lake  fraction
C                  4   FGICE    72x46 glacial ice fraction
C                               (FLAKE + FGRND + FGICE + FOCEAN = 1.000)
C


C                  5   ZATMO    72x46 topography (ocean = 0.0)
C                  6   HOCEAN   72x46 ocean depth
C                  7   HLAKE    72x46 lake  depth
C                  8   HGICE    72x46 glice depth
C                  9   ZSOLID   72x46 topography of solid ground surface
C                  -----------------------------------------------------
C
C     FOLGIZ is for off-line use only, and is not used in GCM radiation.
C     GCM supplies dynamically changing POCEAN,POICE,PEARTH,PLICE values
C     ------------------------------------------------------------------

CF    NRFU=NRFUN(13)
CF    RFILEN(  1:LPATH1)=LLPATH
CF    RFILEN(LP3:LP4)=RFILED
CF    OPEN (NRFU,FILE=RFILEN,FORM='UNFORMATTED',STATUS='OLD')
CF    OPEN (NRFU,FORM='UNFORMATTED',STATUS='OLD')             ! CF

CF    DO 963 K=1,9
CF    READ(NRFU) TITLE,R72X46
CF    DO 962 J=1,46
CF    DO 961 I=1,72
CF    FOLGIZ(I,J,K)=R72X46(I,J)
CF961 CONTINUE
CF962 CONTINUE
CF963 CONTINUE
CF    CLOSE(NRFU)

  999 CONTINUE

      IFIRST=0
 9999 CONTINUE


C     ---------------------------------------------------------------
C     LASTVC     Initialize:  Default Atmospheric Layering, Structure
C                (for Off-Line use)  as Specified by LASTVC Parameter
C                If LASTVC < 0, GCM defines all Radiation Model Input
C     otherwise:
C                Each LASTVC digit(6) specifies a model configuration
C          e.g.:    LASTVC= 123456
C                L=0,1,..9  Layers NL=  Any,GCM12,GCM23,Pset,Hset,etc
C                A=0,1,..6  Atmosphere  Any,Trop,MLS,MLW,SAS,SAW,Std
C                S=0,1,..9  Surf Types  POCEAN=1,PEARTH=1,POICE=1,etc
C                T=0,1,..9  Tracer Aer  Tau=0,  Tau=0.1 Aer Comp(1-9)
C                V=0,1,..9  Vegetation  Sand,Tundra,Grass,Shrubs, etc
C                C=0,1,..9  Cloud,R=10  Clim Cloud Tau in Layer(1,-9)
C                ----------------------------------------------------

CF    IF(LASTVC.GE.0) CALL SETATM


C             -------------------------------------------------------
C             Set Solar Constant for Default Reference Time: Jan 1950
C             Default used for KSOLAR(=1) is that specified in RADCOM
C             -------------------------------------------------------

      JJDAYS=1
      JYEARS=1950
      IF(KJDAYS.GT.0)             JJDAYS=KJDAYS
      IF(KYEARS.GT.0)             JYEARS=KYEARS
C----------------------------------------------
                      CALL SETSOL(JYEARS,JJDAYS)
C----------------------------------------------


C             -------------------------------------------------------
C             Set Default Greenhouse Gas Reference Year to:  Mid 1980
C             Default used for KTREND(=1) is that specified in RADCOM
C             -------------------------------------------------------

      JJDAYG=184
      JYEARG=1980
C----------------------------------------------
                      CALL SETGHG(JYEARG,JJDAYG)
C----------------------------------------------
      IF(KJDAYG.GT.0)             JJDAYG=KJDAYG
      IF(KYEARG.GT.0)             JYEARG=KYEARG
C----------------------------------------------
                      CALL UPDGHG(JYEARG,JJDAYG)
C----------------------------------------------

C--------------------------------
                      CALL SETO3D
                      CALL SETO2A
                      CALL SETGAS
C
                      CALL SETBAK
      IF(MADAER.GT.0) CALL SETAER
      IF(MADDST.GT.0) CALL SETDST
C--------------------------------


C               -----------------------------------------------------
C               Set Volcanic Aerosol Effective Variance Default Value
C                   Particle Size(REFF0=0.3) when not known from data
C                   (VEFF0=0.35 is value based on thermal ISAMS data)
C                   -------------------------------------------------

C----------------------------------------------
      IF(MADVOL.GT.0) CALL SETVOL
C----------------------------------------------

C--------------------------------
                      CALL SETCLD

C--------------------------------

      CALL SOLAR0

      RETURN
      END SUBROUTINE RCOMP1

      SUBROUTINE RCOMPT

      !INCLUDE 'BR00B.COM'

C-----------------------------------------------------------------------
C
C     Time Trend Selection Parameters and Options:
C     -------------------------------------------
C
C             The Nominal Default Values are KYEARX = 0, and KJDAYX = 0,
C             in which case RADCOM supplied Time JYEAR and JDAY are used
C
C             When Non-Zero Values are specified for  KYEARX and KJDAYX,
C             the JYEAR,JDAY Time Dependence of the Specified Process is
C             over-ridden by the Non-Zero KYEARX and KJDAYX Value.
C             ----------------------------------------------------------
C                 Process       KYEARX  KJDAYX
c   KYEARS   ! SolarCon, UV       0       0
c   KYEARG   ! GH Gas Trend       0       0
c   KYEARO   ! Ozone Distr        0       0
c   KYEARA   ! AerClimtolgy       0       0
c   KYEARD   ! Desert Dust        0       0
c   KYEARV   ! Volcanic Aer       0       0
c   KYEARE   ! Epsilon Clds       0       0
c   KYEARR   ! Refl Surface       0       0

C     ------------------------------------------------------------------
C     MADVEL  Model Add-on Data of Extended Climatology Enable Parameter
C             Each MADVEL digit is ON/OFF switch for corresponding input
C             e.g. MADVEL=123456   (zero digit skips input process)
C
C     MADO3M  =  1  Updates  Makiko's 1951-1997 Ozone climatology RFILEA
C                   Updates  Wang-Jacob 1890-1979 Trop O3 changes RFILEB
C     MADAER  =  2  Updates  Aerosol 50y tropospheric climatology RFILE5
C     MADDST  =  3  Updates  Dust-windblown mineral climatology   RFILE6
C     MADVOL  =  4  Updates  Volcanic 1950-00 aerosol climatology RFILE7
C     MADEPS  =  5  Updates  Epsilon cloud heterogeneity data     RFILE8
C     MADLUV  =  6  Updates  Lean's SolarUV 1882-1998 variability RFILE9
C
C                 Related Model Add-on Data Parameters set in BLOCK DATA
C
C     MADGHG   =  1  Default Enables UPDGHG update. (MADGHG=0),no update
C     MADSUR   =  1          V72X46N.1.cor Vegetation type data   RFILEC
C                            Z72X46N Ocean fraction, topography   RFILED
C     ------------------------------------------------------------------

C                      -------------------------------------------------
C                      Set Seasonal and Time (JDAY) Dependent Quantities
C                      -------------------------------------------------

      JJDAYS=JDAY
      JYEARS=JYEAR
      IF(KJDAYS.GT.0)             JJDAYS=KJDAYS
      IF(KYEARS.GT.0)             JYEARS=KYEARS
C----------------------------------------------
      IF(MADLUV.GT.0) CALL UPDSOL(JYEARS,JJDAYS)
C----------------------------------------------

      JJDAYG=JDAY
      JYEARG=JYEAR
      IF(KJDAYG.GT.0)             JJDAYG=KJDAYG
      IF(KYEARG.GT.0)             JYEARG=KYEARG
C----------------------------------------------
      IF(MADGHG.GT.0) CALL UPDGHG(JYEARG,JJDAYG)
C----------------------------------------------


      JJDAYO=JDAY
      JYEARO=JYEAR
      IF(KJDAYO.GT.0)             JJDAYO=KJDAYO
      IF(KYEARO.GT.0)             JYEARO=KYEARO
C----------------------------------------------
                      CALL UPDO3D(JYEARO,JJDAYO)
C----------------------------------------------

      JJDAYA=JDAY
      JYEARA=JYEAR
      IF(KJDAYA.GT.0)             JJDAYA=KJDAYA
      IF(KYEARA.GT.0)             JYEARA=KYEARA
C----------------------------------------------
      IF(MADAER.GT.0) CALL UPDAER(JYEARA,JJDAYA)
C----------------------------------------------

      JJDAYD=JDAY
      JYEARD=JYEAR
      IF(KJDAYD.GT.0)             JJDAYD=KJDAYD
      IF(KYEARD.GT.0)             JYEARD=KYEARD
C----------------------------------------------
      IF(MADDST.GT.0) CALL UPDDST(JYEARD,JJDAYD)
C----------------------------------------------

      JJDAYV=JDAY
      JYEARV=JYEAR
      IF(KJDAYV.GT.0)             JJDAYV=KJDAYV
      IF(KYEARV.GT.0)             JYEARV=KYEARV
C----------------------------------------------
      IF(MADVOL.GT.0) CALL UPDVOL(JYEARV,JJDAYV)
C----------------------------------------------

      JJDAYE=JDAY
      JYEARE=JYEAR
      IF(KJDAYE.GT.0)             JJDAYE=KJDAYE
      IF(KYEARE.GT.0)             JYEARE=KYEARE
C----------------------------------------------
      IF(MADEPS.GT.0) CALL UPDEPS(JYEARE,JJDAYE)
C----------------------------------------------

      JJDAYR=JDAY
      JYEARR=JYEAR
      IF(KJDAYR.GT.0)             JJDAYR=KJDAYR
      IF(KYEARR.GT.0)             JYEARR=KYEARR
C----------------------------------------------
                      CALL UPDSUR(JYEARR,JJDAYR)
C----------------------------------------------

      RETURN
      END SUBROUTINE RCOMPT

      SUBROUTINE RCOMPX

      !INCLUDE 'BR00B.COM'


C     ------------------------------------------------------------------
C     MADVEL  Model Add-on Data of Extended Climatology Enable Parameter
C             Each MADVEL digit is ON/OFF switch for corresponding input
C             e.g. MADVEL=123456   (zero digit skips process)
C
C     MADO3M  =  1           Makiko's 1951-1997 Ozone climatology RFILEA
C     MADAER  =  2  Updates  Aerosol 50y tropospheric climatology RFILE5
C     MADDST  =  3  Updates  Dust-windblown mineral climatology   RFILE6
C     MADVOL  =  4  Updates  Volcanic 1950-00 aerosol climatology RFILE7
C     MADEPS  =  5           Epsilon cloud heterogeneity data     RFILE8
C     MADLUV  =  6           Lean's SolarUV 1882-1998 variability RFILE9
C
C                 Related Model Add-on Data Parameters set in BLOCK DATA
C
C     MADGHG   =  1  Default Enables UPDGHG update. (MADGHG=0),no update
C     MADSUR   =  1          V72X46N.1.cor Vegetation type data   RFILEC
C                            Z72X46N Ocean fraction, topography   RFILED
C     ------------------------------------------------------------------
C
C      -----------------------------------------------------------------
C      Get Surface, Atmosphere, Sun Angle, Radiative Forcing, etc. Input
C      to compute Solar/Thermal Radiation for given (JLAT,ILON) Grid-box
C
C      The Radiation Model utilizes Data with 72x46 (lon,lat) resolution
C                 for GCM resolution other than 72x46, set JLAT and ILON
C                 to appropriately Sample  (rather than interpolate) the
C                 72x46 aerosol, ozone, cloud heterogeneity data sets
C
C      The Radiation Model can accommodate arbitrary vertical resolution
C      -----------------------------------------------------------------


C--------------------------------
                      CALL GETO3D
                      CALL GETGAS
C--------------------------------


C--------------------------------
      IF(MADBAK.GT.0) THEN ; CALL GETBAK
       ELSE ; SRBEXT=1.d-20 ; SRBSCT=0. ; SRBGCB=0. ; TRBALK=0. ; END IF
      IF(MADAER.GT.0) THEN ; CALL GETAER
       ELSE ; SRAEXT=0.     ; SRASCT=0. ; SRAGCB=0. ; TRAALK=0. ; END IF
      IF(MADDST.GT.0) THEN ; CALL GETDST
       ELSE ; SRDEXT=0.     ; SRDSCT=0. ; SRDGCB=0. ; TRDALK=0. ; END IF
      IF(MADVOL.GT.0) THEN ; CALL GETVOL
       ELSE ; SRVEXT=0.     ; SRVSCT=0. ; SRVGCB=0. ; TRVALK=0. ; END IF
C--------------------------------


C--------------------------------  (GETSUR sets albedo needed by GETCLD)
                      CALL GETSUR
                      CALL GETEPS
                      CALL GETCLD
C--------------------------------

C--------------------------------
                      CALL THERML

                      CALL SOLARM
C--------------------------------

      RETURN
      END SUBROUTINE RCOMPX

      SUBROUTINE SETSOL(JYEARS,JJDAYS)

      !INCLUDE 'BR00B.COM'

      DATA LMOREF/0/ ; SAVE LMOREF

C-----------------------------------------------------------------------
C
C     SETSOL Parameters:
C----------------------
C            KSOLAR    Selects Solar Spectrum, (Lean vs Thekaekara Flux)
C            JYEARS    JYEAR Proxy:  Sets: Solar Constant Reference Year
C            JJDAYS    JDAY  Proxy:  Sets Reference Year Month JDAY/30.5
C                      (Nominal Reference: JYEARS= 1950 JJDAYS= January)
C
C-----------------------------------------------------------------------
C            KSOLAR   SOLSPEC     UVWAVLs       UVFACTs         KUVFAC
C-----------------------------------------------------------------------
C              -1     THEK      Can be set    Can be set   (if KUVFAC=1)
C-----------------------------------------------------------------------
C               0     LEAN      Can be set    Can be set   (if KUVFAC=1)
C-----------------------------------------------------------------------
C               1     LEAN      Can be set    Can be set   (if KUVFAC=1)
C-----------------------------------------------------------------------
C
C                               (Option to Modify Solar UV Fluxes)
C            UVWAVL    Specified Edges of UV Flux Variation SubIntervals
C            UVFACT    Factors to Change the Amplitude of UV Variability
C
C            KUVFAC    ON/OFF switch for activating UV Flux Modification
C            KSNORM    Re-Normalize S0 (VIS) (after UV Amplitude Change)
C                         (Nominal UVWAVLs are: 0.295,0.310,0.366)
C
C-----------------------------------------------------------------------
C     SETSOL Output:
C------------------
C
C        AO3 =   Ozone Absorption Table AO3(460) found in COMMON/BLOCKD/
C                (Solar UV Flux Weighted Absorption Table is used by the
C                FUNCTION AO3ABS(OCM) in SOLAR to compute Ozone Heating)
C                AO3 is the fraction of total Solar Flux absorbed by O3.
C
C     S00WM2 =   Solar Constant Reference Value for Time = JYEARS,JJDAYS
C                (Thekaekara, if KSOLAR=-1, Reference = 1367 WATTS/M**2)
C
C
C     SETSOL  is Generally Called once at Model Initialization to Select
C                Solar Flux (LEAN,THEK), and to Define S00WM2 (RATLS0=1)
C
C-----------------------------------------------------------------------
C NOTE:
C-----
C     S00WM2 = Nominal Reference Solar Constant 1366.448785D0 WATTS/M**2
C                (Spectral Integral: Lean99 Solar Flux for January 1950)
C
C     KSOLAR=-1  Reproduces Thekaekhara Ozone Absorption, e.g., XRAD83XX
C     KSOLAR= 0  Uses Lean99 Solar Flux as set for Time= (JYEARS,JJDAYS)
C     KSOLAR= 1  Sets Lean99 Solar Flux to Current Time= (JYEARS,JJDAYS)
C     KSOLAR= 2  same as 1 but based on annual (not monthly) data
C                (JJDAYS used to select the specified Monthly-Mean Flux)
C
C-----------------------------------------------------------------------
C
C     UPDSOL Parameters:
C----------------------
C            JYEARS    JYEAR Proxy:  Selects Solar Constant Current Year
C            JJDAYS    JDAY  Proxy:  Selects Lean Data Month JJDAYS/30.5
C
C     UPDSOL Output:
C------------------
C
C        AO3 =   Ozone Absorption Table AO3(460) found in COMMON/BLOCKD/
C                (Solar UV Flux Weighted Absorption Table is used by the
C                FUNCTION AO3ABS(OCM) in SOLAR to compute Ozone Heating)
C                AO3 is the fraction of total Solar Flux absorbed by O3.
C
C     RATLS0 =   Ratio:  Current-Time Solar Constant to Reference S00WM2
C
C-----------------------------------------------------------------------
C     Remark:
C
C     UPDSOL  is Called in RCOMPT to Update Solar Constant and Ozone AO3
C                Solar UV Absorption Dependence.  (Monthly-Mean Data are
C                NOT Interpolated in Time, but get Updated with Changing
C                Month, i.e., whenever JDAY/30.5 Reaches Integer Value.)
C
C-----------------------------------------------------------------------


      IF(KSOLAR.GE.0) THEN
C                                           Lean99 Solar Flux, UV Option
C                                           ----------------------------
      if(Ksolar.lt.2) then    ! monthly data
        JMO=1+JJDAYS/30.5D0
        IF(JMO.GT.12) JMO=12
        LMO=(JYEARS-iy1S0)*12+JMO
        IF(LMO.GT.Ms0X) LMO=LMO-mcycs0*((LMO-Ms0X+mcycs0-1)/mcycs0)
        IF(LMO.LT.1) LMO=LMO+mcycs0*((mcycs0-lmo)/mcycs0)
      else                    ! annual data
        Is0x = nint( yr2s0-yr1s0+1 )
        lmo = nint( jyears - yr1s0 + 1.5 )
        IF(LMO.GT.Is0X) LMO=LMO-icycs0*((LMO-Is0X+icycs0-1)/icycs0)
        IF(LMO.LT.1) LMO=LMO+icycs0*((icycs0-lmo)/icycs0)
      end if
      LMOREF=LMO

C                        IF(MADLUV.EQ.0) Default Option is then in force
C                        Default (FRLEAN) = Lean 1950 Jan Solar, UV flux
C                        CORFAC accounts for DSLEAN units in BLOCK DATA,
C                        and TSI1/TSI2 normalization of Lean input data.
C                        -----------------------------------------------

      CORFAC=1366.2911D0/1366.4487855D0
      FLXSUM=0.D0
      DO 110 K=1,190
      IF(MADLUV.EQ.0) FLXSUM=FLXSUM+FRLEAN(K)*DSLEAN(K)*CORFAC
      IF(MADLUV.GT.0) FLXSUM=FLXSUM+UVLEAN(LMO,K)*DSLEAN(K)
  110 CONTINUE
      S00WM2=FLXSUM

      I=0
      DO 120 K=1,50
      I=I+1
      WSOLAR(I)=W1LEAN(K)
      IF(MADLUV.EQ.0) FSOLAR(I)=FRLEAN(K)
      IF(MADLUV.GT.0) FSOLAR(I)=UVLEAN(LMO,K)
      I=I+1
      WSOLAR(I)=W1LEAN(K+1)
      FSOLAR(I)=FSOLAR(I-1)
  120 CONTINUE
      NWSUV=100
C                                          Thekaekhara Solar Flux Option
C                                          -----------------------------
      ELSE
      DO 130 I=1,190
      WSOLAR(I)=WTHEK(I)
      FSOLAR(I)=FTHEK(I)
  130 CONTINUE
      S00WM2=1367.D0
      LMOREF=-111
      NWSUV=190
      ENDIF

C                                         Option to Modify Solar UV Flux
C                                         ------------------------------
      IF(KUVFAC.EQ.1) THEN
      F1FLUX=0.D0
      F2FLUX=0.D0
      F3FLUX=0.D0
      DO 160 I=1,NWSUV
      IF(WSOLAR(I).GT.UVWAVL(1)) GO TO 140
      F1FLUX=F1FLUX+FSOLAR(I)
      FSOLAR(I)=FSOLAR(I)*UVFACT(1)
      GO TO 160
  140 CONTINUE
      IF(WSOLAR(I).GT.UVWAVL(2)) GO TO 150
      F2FLUX=F2FLUX+FSOLAR(I)
      FSOLAR(I)=FSOLAR(I)*UVFACT(2)
      GO TO 160
  150 CONTINUE
      IF(WSOLAR(I).GT.UVWAVL(3)) GO TO 170
      F3FLUX=F3FLUX+FSOLAR(I)
      FSOLAR(I)=FSOLAR(I)*UVFACT(3)
  160 CONTINUE
  170 CONTINUE
      UVNORM=F1FLUX*(1.D0-UVFACT(1))+F2FLUX*(1.D0-UVFACT(2))
     +      +F3FLUX*(1.D0-UVFACT(3))
      IF(KSNORM.EQ.0) S00WM2=S00WM2-UVNORM
      ENDIF
C                  -----------------------------------------------------
C                  When KUVFAC=1 option multiplicative factors UVFACT(I)
C                  are used to change the UV spectral flux distribution,
C                  KSNORM=1 provides the option to keep S00WM2 constant.
C                  -----------------------------------------------------

      RATLS0=1.D0

      DO 190 I=1,460
      II=(I-10)/90-4
      XX=I-((I-10)/90)*90
      OCM=XX*10.D0**II
      DO 180 J=1,226
      TAUK=FUVKO3(J)*OCM
      IF(TAUK.GT.35.D0) TAUK=35.D0
      UVA(J)=1.D0-DEXP(-TAUK)
  180 CONTINUE
      UVWAVA=0.100D0
      UVWAVB=0.400D0
      CALL FXGINT(UVA,XWAVO3,226,FSOLAR,WSOLAR,NWSUV,UVWAVA,UVWAVB,AO33)
      AO3(I)=AO33/S00WM2
  190 CONTINUE

C                       ------------------------------------------------
C                       NOTE:  AO3 is the Ozone-path Absorption Function
C                              AO3 convolves O3 asborption with solar UV
C                              spectral variations by FXGINT integration
C                              AO3 is expressed as the absorbed fraction
C                              of the total solar flux (S00WM2=1366W/m2)
C                              -----------------------------------------
      RETURN


C--------------------------------
      ENTRY UPDSOL(JYEARS,JJDAYS)
C--------------------------------


      IF(KSOLAR.LT.1) GO TO 300
      IF(JYEARS.LT.1) GO TO 300

      if(Ksolar.lt.2) then    ! monthly data       ksolar=1
        JMO=1+JJDAYS/30.5D0
        IF(JMO.GT.12) JMO=12
        LMO=(JYEARS-iy1S0)*12+JMO
        IF(LMO.GT.Ms0X) LMO=LMO-mcycs0*((LMO-Ms0X+mcycs0-1)/mcycs0)
        IF(LMO.LT.1) LMO=LMO+mcycs0*((mcycs0-lmo)/mcycs0)
      else                    ! annual data        ksolar=2
        Is0x = nint( yr2s0-yr1s0+1 )
        lmo = nint( jyears - yr1s0 + 1.5 )
        IF(LMO.GT.Is0X) LMO=LMO-icycs0*((LMO-Is0X+icycs0-1)/icycs0)
        IF(LMO.LT.1) LMO=LMO+icycs0*((icycs0-lmo)/icycs0)
      end if

      IF(LMO.EQ.LMOREF) GO TO 300
      LMOREF=LMO

      IF(KSOLAR.GE.0) THEN
C                                               Select Lean99 Solar Flux
C                                               ------------------------
      FLXSUM=0.D0
      DO 210 K=1,190
      FLXSUM=FLXSUM+UVLEAN(LMO,K)*DSLEAN(K)
  210 CONTINUE

      I=0
      DO 220 K=1,50
      I=I+1
      WSOLAR(I)=W1LEAN(K)
      FSOLAR(I)=UVLEAN(LMO,K)
      I=I+1
      WSOLAR(I)=W1LEAN(K+1)
      FSOLAR(I)=FSOLAR(I-1)
  220 CONTINUE
      NWSUV=100
      ELSE
C                                          Select Thekaekhara Solar Flux
C                                          -----------------------------
      DO 230 I=1,190
      WSOLAR(I)=WTHEK(I)
      FSOLAR(I)=FTHEK(I)
  230 CONTINUE
      FLXSUM=1367.D0
      LMOREF=-111
      NWSUV=190
      ENDIF
C                                         Option to Modify Solar UV Flux
C                                         ------------------------------
      IF(KUVFAC.EQ.1) THEN
      F1FLUX=0.D0
      F2FLUX=0.D0
      F3FLUX=0.D0
      DO 260 I=1,NWSUV
      IF(WSOLAR(I).GT.UVWAVL(1)) GO TO 240
      F1FLUX=F1FLUX+FSOLAR(I)
      FSOLAR(I)=FSOLAR(I)*UVFACT(1)
      GO TO 260
  240 CONTINUE
      IF(WSOLAR(I).GT.UVWAVL(2)) GO TO 250
      F2FLUX=F2FLUX+FSOLAR(I)
      FSOLAR(I)=FSOLAR(I)*UVFACT(2)
      GO TO 260
  250 CONTINUE
      IF(WSOLAR(I).GT.UVWAVL(3)) GO TO 270
      F3FLUX=F3FLUX+FSOLAR(I)
      FSOLAR(I)=FSOLAR(I)*UVFACT(3)
  260 CONTINUE
  270 CONTINUE
      UVNORM=F1FLUX*(1.D0-UVFACT(1))+F2FLUX*(1.D0-UVFACT(2))
     +      +F3FLUX*(1.D0-UVFACT(3))
      IF(KSNORM.EQ.0) FLXSUM=FLXSUM-UVNORM
      ENDIF

      RATLS0=FLXSUM/S00WM2

      DO 290 I=1,460
      II=(I-10)/90-4
      XX=I-((I-10)/90)*90
      OCM=XX*10.D0**II
      DO 280 J=1,226
      TAUK=FUVKO3(J)*OCM
      IF(TAUK.GT.35.D0) TAUK=35.D0
      UVA(J)=1.D0-DEXP(-TAUK)
  280 CONTINUE
      UVWAVA=0.100D0
      UVWAVB=0.400D0
      CALL FXGINT(UVA,XWAVO3,226,FSOLAR,WSOLAR,NWSUV,UVWAVA,UVWAVB,AO33)
      AO3(I)=AO33/FLXSUM
  290 CONTINUE

  300 CONTINUE

      RETURN
      END SUBROUTINE SETSOL

      SUBROUTINE SETGHG(JYEARG,JJDAYG)
C
      !INCLUDE 'BR00B.COM'
C
C     ---------------------------------------------------------------
C     SETGHG  Sets Default Greenhouse Gas Reference Year (for FULGAS)
C
C     Control Parameter:
C                     KTREND (specified in RADCOM) activates GH Trend
C               Default
C     KTREND  =   1
C     Selects   GTREND
C     ---------------------------------------------------------------
C
      TREF=JYEARG+(JJDAYG-0.999D0)/366.D0
C
      IF(KTREND.EQ.0) THEN
        XREF(1)=PPMV80(2)
        XREF(2)=PPMV80(6)
        XREF(3)=PPMV80(7)
        XREF(4)=PPMV80(8)*1000.D0
        XREF(5)=PPMV80(9)*1000.D0
        XREF(6)=PPMV80(11)*1000.D0  ! YREF11=PPMV80(11)*1000.D0
        XREF(7)=PPMV80(12)*1000.D0  ! ZREF12=PPMV80(12)*1000.D0
        RETURN
      END IF

CF    IF(KTREND.EQ.6) THEN
      CALL GTREND(XREF,TREF)     ! finds xref 1-6 (yref11=xx6=xref(6))
      XREF(7)=1.D-13             ! ZREF12=1.D-13
CF    ENDIF
      DO 120 I=1,NGHG
      IF(XREF(I).LT.1.D-06) XREF(I)=1.D-06
  120 CONTINUE
      PPMV80(2)=XREF(1)
      PPMV80(6)=XREF(2)
      PPMV80(7)=XREF(3)
      PPMV80(8)=XREF(4)/1000.D0
      PPMV80(9)=XREF(5)/1000.D0
      PPMV80(11)=XREF(6)/1000.d0   ! YREF11/1000.D0
      PPMV80(12)=XREF(7)/1000.D0   ! ZREF12/1000.D0
      RETURN
C
C--------------------------------
      ENTRY UPDGHG(JYEARG,JJDAYG)
C--------------------------------
C
      TNOW=JYEARG+(JJDAYG-0.999D0)/366.D0
C
      IF(KTREND.EQ.0) THEN
        FULGAS(2)=PPMVK0(2)/XREF(1)
        FULGAS(6)=PPMVK0(6)/XREF(2)
        FULGAS(7)=PPMVK0(7)/XREF(3)
        FULGAS(8)=PPMVK0(8)/XREF(4)
        FULGAS(9)=PPMVK0(9)/XREF(5)
        FULGAS(11)=PPMVK0(11)/XREF(6) ! YREF11
        FULGAS(12)=PPMVK0(12)/XREF(7) ! .../ZREF12
        RETURN
      END IF

CF    IF(KTREND.EQ.1) THEN
      CALL GTREND(XNOW,TNOW) ! finds xnow 1-6 (ynow11=xx6=xnow(6))
      XNOW(7)=1.D-20         ! ZNOW12=1.D-20
CF    ENDIF
      FULGAS(2)=XNOW(1)/XREF(1)
      FULGAS(6)=XNOW(2)/XREF(2)
      FULGAS(7)=XNOW(3)/XREF(3)
      FULGAS(8)=XNOW(4)/XREF(4)
      FULGAS(9)=XNOW(5)/XREF(5)
      FULGAS(11)=XNOW(6)/XREF(6) ! YNOW11/YREF11
      FULGAS(12)=XNOW(7)/XREF(7) ! ZNOW12/ZREF12
C
      RETURN
      END SUBROUTINE SETGHG

      SUBROUTINE SETO3D

      !INCLUDE 'BR00B.COM'

      real*8 D72(72),ATOT(32),O3XCOL(LX),VO3LON(360)
      SAVE  IFIRST,MPI,MPJ,WTMPI,WTMPJ,MMOREF,NL0,NL1,N150MB,NTOPMB
      DATA IFIRST/1/,MMOREF/0/,NL0/0/,N150MB/8/

C     ---------------------------------------------------------------
C     Ozone Data : Time Trend, Longitude Dependence Selection Options
C
C        O3Data   Source    Time Trend   VD Profile      L-Dependence
C     A  O3DPPM   McPeters   Seasonal    0 - 60 km         Zonal
C     B  O3AVE+   London     Seasonal    Column+Profile  Longitudinal
C     C  O3CLIM   Makiko    1850-2050    44-Layer Prof     Zonal
C     D  WJ1890+  Wang-J    1890,1979    Below 150mb     Longitudinal
C     E  TOPO3    Keating    Seasonal    Above 1(or 5)mb   Zonal
C     ---------------------------------------------------------------
C     NOTE:  C,D Data (O3 Clim,W-J Lon dep) are only read if MADO3M>0
C
C     Control Parameters:  MOZONE, KO3LON, O3WJT0, FULGAS(13), NCARO3
C               Defaults:     5       0     0.001     1.0         0
C     Options:
C     -------
C              Zonal       Time          Longitudinal     Additional
C     MOZONE  O3 Prof      Trend          Dependence       Options
C     ------  -------    -----------    ---------------   -----------
C       0       A         Seasonal      B if(KO3LON=1)       (1)
C       1       A        A + Delta C    B if(KO3LON=1)        *
C       2       C            C          B if(KO3LON=1)
C       3       C          C + D        D (below 150mb)      (2)
C       4       C        C + D + B      D (below 150mb)      (2)
C                                     + B (above 150mb)      (3)
C       5      C/E     C/E + D + B      D (below 150mb)      (2)
C                                     + B (above 150mb)      (3) **
C       6    =5 but no lin.transition from E->C between 5mb and 1mb
C
C       8    =6 but uses Wang-Jacob 1890 data below 150 mb
C       9    =6 but uses Wang-Jacob 1979 data below 150 mb
C     ---------------------------------------------------------------
C     NOTE:  (1)  If NCARO3=1, then London's NCAR O3 data is selected
C            (2)  O3WJT0 extends W-J time rate lon dep change (years)
C            (3)  FULGAS(13) scales B lon dep amplitude (above 150mb)
C                 (When FULGAS(13)=0.0, MOZONE=4 is same as MOZONE=3)
C
C     REMARKS:
C              *  If A + Delta C is negative, O3 amount is set to 0.0
C
C             **  MOZONE=5 uses SAGE2 (Latitudinal,Seasonal) pressure
C                 height dependence rather, than global mean profile.
C                 Keating(85) data (E TOPO3) are used above 1mb level
C                 with linear transition to (C O3CLIM) data below 5mb
C     ---------------------------------------------------------------

      IF(IFIRST.EQ.1) THEN
      NL0=NL
      NL1=NL+1
      DO 100 L=1,NL1
      IF(PLB(L).LE.150.D0) GO TO 100
      N150MB=L
  100 CONTINUE
      NTOPMB=N150MB+1

C                     Extrapolate missing O3DPPM data at winter poles
C                     -----------------------------------------------
      DO 112 M=5,7
      DO 111 L=1,31
      O3DPPM(L,1,M)=2.D0*O3DPPM(L,2,M)-O3DPPM(L,3,M)
  111 CONTINUE
  112 CONTINUE
      DO 122 MM=10,12
      M=MM
      IF(MM.EQ.10) M=1
      DO 121 L=1,31
      O3DPPM(L,16,M)=2.D0*O3DPPM(L,15,M)-O3DPPM(L,14,M)
  121 CONTINUE
  122 CONTINUE

C               -----------------------------------------------------
C               Scale Height of Homogeneous Atm at 0 deg D = 7.991 km
C               Conversion of PPM to ACM: ACM per PPM = 7.991D5*1.D-6
C         Thus:
C               O3DPPM is in ACM (1 ACM = 1000 DU = 2.687E+19mol/cm2)
C               (ACM per 2-km layer, at J= -75, -65, ... 65, 75 dLat)
C               -----------------------------------------------------

      DO 131 L=1,31
      ZO3(L)=2.D0*(L-1)
      PO3(L)=1013.25D0*10.D0**(-ZO3(L)/16.D0)
  131 CONTINUE
      PO3(32)=0.D0
      DO 143 M=1,12
      DO 142 J=1,16
      DO 141 L=1,30
      ACML=7.991D5*1.D-6*(PO3(L)-PO3(L+1))/1013.25D0
      O3DPPM(L,J,M)=ACML*0.5D0*(O3DPPM(L,J,M)+O3DPPM(L+1,J,M))
  141 CONTINUE
      L=31
      ACML=7.991D5*1.D-6*(PO3(L)-PO3(L+1))/1013.25D0
      O3DPPM(L,J,M)=ACML*0.5D0*(O3DPPM(L,J,M)+0.D0)
  142 CONTINUE
  143 CONTINUE

C          ----------------------------------------------------------
C          Define ALPHA Wang-Jacob tropospheric ozone trend parameter
C             for exponential growth based on 1979-1890 distributions
C             -------------------------------------------------------

      IF(MADO3M.GT.0) THEN
      DO 154 M=1,12
      DO 153 L=1,N150MB
      DO 152 J=1,46
      DO 151 I=1,72
      ALPHA=LOG(WJ1979(I,J,L,M)/WJ1890(I,J,L,M))/(1979.D0-1890.D0)
      O3WJA(I,J,L,M)=ALPHA
  151 CONTINUE
  152 CONTINUE
  153 CONTINUE
  154 CONTINUE
      ENDIF
      IFIRST=0
      ENDIF

C              ------------------------------------------------------
C              Define Annual Average Ozone Distribution (NL x MLAT46)
C              ------------------------------------------------------

      DO 240 J=1,MLAT46
      XDLAT=DLAT46(J)
      XXLAT=0.1D0*(XDLAT+85.D0)
      JI=XXLAT
      WJJ=XXLAT-JI
      WJI=1.D0-WJJ
      JJ=JI+1
      IF(JI.LT.1) THEN
      JI=1
      JJ=1
      ENDIF
      IF(JI.GT.15) THEN
      JI=16
      JJ=16
      ENDIF

      DO 220 L=1,31
      OSUM=0.D0
      DO 210 M=1,12
      OSUM=OSUM+WJI*O3DPPM(L,JI,M)+WJJ*O3DPPM(L,JJ,M)
  210 CONTINUE
      O3XCOL(L)=OSUM/12.D0
  220 CONTINUE
      CALL REPART(O3XCOL,PO3,32,ATOT,PLB0,NL1)
      SUM=0.D0
      DO 230 L=1,NL0
      O3DLJ(L,J)=ATOT(L)
      SUM=SUM+ATOT(L)
  230 CONTINUE
!obso COLSUM(J)=SUM
      L=NL1
  240 CONTINUE
      DO 250 L=1,NL0
      U0GAS(L,3)=O3DLJ(L,JLAT)
  250 CONTINUE

      CALL SETO3L

      RETURN

C--------------------------------
      ENTRY UPDO3D(JYEARO,JJDAYO)
C--------------------------------

      XJDAY=JJDAYO-0.999D0           ! needed for MOZONE<2 and MOZONE>7
      XPMO=XJDAY/30.5D0+.5D0
      MPI=XPMO
      WTMPJ=XPMO-MPI
      WTMPI=1.D0-WTMPJ
      IF(MPI.LT.1) MPI=12
      MPJ=MPI+1
      IF(MPJ.GT.12) MPJ=1
      IF(MOZONE.LT.2) THEN

C     -----------------------------------------------------------------
C        Ozone Data Latitude Coverage is given at 10 degree Intervals
C     (-75,-65,-55,-45,-35,-25,-15,-05,  5, 15, 25, 35, 45, 55, 65, 75)
C        Latitude Interpolation is Linear between -75 and 75 degrees
C                 NO Extrapolation Pole-ward of -75 or 75
C     -----------------------------------------------------------------

      DO 330 J=1,MLAT46
      XDLAT=DLAT46(J)
      XXLAT=0.1D0*(XDLAT+85.D0)
      JI=XXLAT
      WTJJ=XXLAT-JI
      WTJI=1.D0-WTJJ
      JJ=JI+1
      IF(JI.LT.1) THEN
      JI=1
      JJ=1
      ENDIF
      IF(JI.GT.15) THEN
      JI=16
      JJ=16
      ENDIF

      DO 310 L=1,31
      O3XCOL(L)=(WTMPI*O3DPPM(L,JI,MPI)+WTMPJ*O3DPPM(L,JI,MPJ))*WTJI
     +         +(WTMPI*O3DPPM(L,JJ,MPI)+WTMPJ*O3DPPM(L,JJ,MPJ))*WTJJ
  310 CONTINUE
      CALL REPART(O3XCOL,PO3,32,ATOT,PLB0,NL1)
      SUM=0.D0
      DO 320 L=1,NL0
      O3DLJ(L,J)=ATOT(L)
      SUM=SUM+ATOT(L)
  320 CONTINUE
!obso COLSUM(J)=SUM
      L=NL1
  330 CONTINUE
      IF(MOZONE.LT.1) GO TO 380
      ENDIF

C              ---------------------------------------------------------
C              Include Wang-Jacob tropospheric ozone longitudinal change
C               (imposed on zonal-mean ozone profile below 150 mb level)
C              O3WJT0= time beyond 1979 for (1890-1979) ALPHA to go to 0
C
C              IF(MOZONE.EQ.4) London's longitudinal ozone dependence is
C              also imposed on the zonal-mean ozone profile above 150 mb
C
C              IF(MOZONE.gt.4)  Keating (1985) upper stratospheric ozone
C              is added above 1 mb pressure level (transition 1 - 0.4mb)
C              ---------------------------------------------------------

      IF(MOZONE.GT.2) THEN
      MMO=1+JJDAYO/30.5D0
      IF(MMO.GT.12) MMO=12
      IF(MMO.EQ.MMOREF) GO TO 337
      MMOREF=MMO
      DTYEAR=JYEARO-1890.D0
      KXTRAP=0
      IF(JYEARO.GT.1979) THEN
      KXTRAP=1
      DTYEAR=(O3WJT0-(JYEARO-1979.D0))/(O3WJT0+0.001D0)
      ENDIF
      IF(DTYEAR.LT.0.D0) DTYEAR=0.D0
      DO 334 L=1,N150MB
      DO 333 J=1,46
      SUM72=0.D0
      DO 331 I=1,72
      ALPHA=O3WJA(I,J,L,MMO)
      IF(KXTRAP.EQ.0) D72(I)=WJ1890(I,J,L,MMO)*EXP(ALPHA*DTYEAR)
      IF(KXTRAP.EQ.1) D72(I)=WJ1979(I,J,L,MMO)*EXP(ALPHA*DTYEAR)
      SUM72=SUM72+D72(I)
  331 CONTINUE
      SUM72=SUM72/72.D0
      DO 332 I=1,72
      O3LF(I,J,L)=D72(I)/SUM72
  332 CONTINUE
  333 CONTINUE
  334 CONTINUE
      DO 336 J=1,46
      DO 335 I=1,72
      IF(MOZONE.LT.4) O3LF(I,J,NTOPMB)=1.D0
      IF(MOZONE.GE.4) O3LF(I,J,NTOPMB)=1.D0
     +               +(O3NCAR(I,J,MMO)-1.D0)*FULGAS(13)
  335 CONTINUE
  336 CONTINUE

  337 CONTINUE
      ENDIF

      MONTHZ=12*(JYEARO-iy1O3)
      IF(MONTHZ.LT.0)  MONTHZ=0
      IF(MONTHZ.GT.MO3X-12) MONTHZ=MO3X-12
      XJDAY=JJDAYO-0.999D0
      XMMO=XJDAY/30.5D0+.5D0
      MMI=XMMO
      WTMMJ=XMMO-MMI
      WTMMI=1.D0-WTMMJ
      MMI=MMI+MONTHZ
      MMJ=MMI+1
      IF(MMI.LT.1) MMI=12
      IF(MMI.EQ.MO3X) MMJ=MO3X-11

      IF(MOZONE.GT.1) THEN
      DO 350 J=1,46
      SUM=0.0
      DO 340 L=1,NL0
      O3DLJ(L,J)=WTMMI*O3CLIM(MMI,L,J)+WTMMJ*O3CLIM(MMJ,L,J)
      SUM=SUM+O3DLJ(L,J)
  340 CONTINUE
  350 CONTINUE
      GO TO 380
      ENDIF

      DO 370 J=1,46
      DO 360 L=1,NL0
      DELTAO=WTMMI*O3CLIM(MMI,L,J)+WTMMJ*O3CLIM(MMJ,L,J)
     +      -WTMMI*O3CLIM(MPI,L,J)-WTMMJ*O3CLIM(MPJ,L,J)
      O3DLJ(L,J)=O3DLJ(L,J)+DELTAO
      IF(O3DLJ(L,J).LT.1.D-06) O3DLJ(L,J)=1.D-06
  360 CONTINUE
  370 CONTINUE
  380 CONTINUE
      DO 390 L=1,NL0
      U0GAS(L,3)=O3DLJ(L,JLAT)
  390 CONTINUE

      IF(KO3LON.EQ.1) CALL UPDO3L
      RETURN


C-----------------
      ENTRY GETO3D
C-----------------

      IF(MOZONE.GT.2) THEN
      DO 610 L=1,N150MB
      U0GAS(L,3)=O3DLJ(L,JLAT)*O3LF(ILON,JLAT,L)
  610 CONTINUE
      DO 620 L=NTOPMB,NL0
      U0GAS(L,3)=O3DLJ(L,JLAT)*O3LF(ILON,JLAT,NTOPMB)
  620 CONTINUE
      ENDIF

C          -------------------------------------------------------------
C          Option to include London's NCAR ozone longitudinal dependence
C          - but not the Wang-Jacob tropospheric longitudinal dependence
C          -------------------------------------------------------------

      IF(MOZONE.LT.3) THEN
      FACTOR=1.D0
      IF(KO3LON.EQ.1) THEN
      CALL GETLON(VO3LON,O3Z0AV)
      FACTOR=1.D0+(VO3LON(ILON)-1.D0)*FULGAS(13)
      ENDIF
      DO 710 L=1,NL0
      U0GAS(L,3)=O3DLJ(L,JLAT)*FACTOR
  710 CONTINUE
      ENDIF

C          -------------------------------------------------------------
C          Mozone=8: replace Ozone below 150mb by 1890 Wang-Jacobs data
C          Mozone=9: replace Ozone below 150mb by 1979 Wang-Jacobs data
C          -------------------------------------------------------------

      IF (MOZONE.EQ.8) THEN
        DO L=1,N150MB
        U0GAS(L,3)=WJ1890(ILON,JLAT,L,MPI)*WTMPI +
     +             WJ1890(ILON,JLAT,L,MPJ)*WTMPJ
        END DO
      ELSE IF (MOZONE.EQ.9) THEN
        DO L=1,N150MB
        U0GAS(L,3)=WJ1979(ILON,JLAT,L,MPI)*WTMPI +
     +             WJ1979(ILON,JLAT,L,MPJ)*WTMPJ
        END DO
      END IF

      RETURN
      END SUBROUTINE SETO3D

      SUBROUTINE SETO3L

      !INCLUDE 'BR00B.COM'


C-----------------------------------------------------------------------
C
C        London et al (1976) Jul,1957-Dec,1970 NCAR Atlas of Total Ozone
C
C        Average Global Column Amount -- O3AVE(Month,Latitude,Longitude)
C
C                                  Month=1-12  Jac,Feb,...,Dec
C                                  Lat  =1-18  -85,-75,..., 85
C
C        (these data are used to define Ozone longitudinal dependence)
C-----------------------------------------------------------------------


      SAVE  IFIRST,IMO,JMO,WTIM,WTJM,WTSPR,WTAUT
      SAVE    XJDMO,HKMSPR,HKMAUT,CNCAUT,CNCSPR,DEGLAT,
     *        PLBSO3,SOJDAY,PMLAT,AO3JIM,O3LB,CONCS,CONCA,
     *        BHKMS,BHKMA,WTJLAT,WTJLON,ILATIJ,ILONIJ,
     *        WTLSEP,WTLJAN,LSEPJ,LJANJ, NL0 ,NL1

      real*8 XJDMO(14),HKMSPR(14),HKMAUT(14),D72(72)
      DIMENSION CNCAUT(14),CNCSPR(14),DEGLAT(14)

      DIMENSION AO3JIM(144),O3LB(LX),         WO3LON(360)
      DIMENSION CONCS(144),CONCA(144),BHKMS(144),BHKMA(144)
      DIMENSION WTJLAT(144),WTJLON(144),ILATIJ(144),ILONIJ(144)
      DIMENSION WTLSEP(144),WTLJAN(144),LSEPJ(144),LJANJ(144)
      PARAMETER(ACMMGG=2.37251E-4,ACMPKM=7.1509E-4,H10MB=31.05467)
      DATA IFIRST/1/
      DATA DEGLAT/-85.0,-71.0,-59.0,-47.0,-35.0,-22.0,-9.0,
     +            9.0,22.0,35.0,47.0,59.0,71.0,85.0/
      DATA XJDMO/-15.0,16.0,45.0,75.0,105.0,136.0,166.0,197.0,228.0
     +           ,258.0,289.0,319.0,350.0,381.0/
      DATA HKMSPR/18.5,18.5,19.0,23.5,24.0,24.5,26.5,
     +            26.5,25.0,22.5,21.0,20.0,18.5,16.5/
      DATA HKMAUT/16.5,18.5,20.0,21.0,22.5,25.0,26.5,
     +            26.5,24.5,24.0,23.5,19.0,18.5,18.5/
      DATA CNCSPR/0.0181,0.0212,0.0187,0.0167,0.0162,0.0183,0.0175,
     +            0.0187,0.0200,0.0196,0.0225,0.0291,0.0287,0.0300/
      DATA CNCAUT/0.0300,0.0287,0.0291,0.0225,0.0196,0.0200,0.0187,
     +            0.0175,0.0183,0.0162,0.0167,0.0187,0.0212,0.0181/
C
      DIMENSION PLBSO3(11),SOJDAY(6),PMLAT(6)
      DATA PLBSO3/10.0,7.0,5.0,3.0,2.0,1.5,1.0,0.7,0.5,0.3,0.1/
      DATA SOJDAY/-91.,31.,92.,213.,274.,396./
      DATA PMLAT/1.,1.,-1.,-1.,1.,1./
C-------------------------------------------------------------------
C    Set O3 Vertical Profile Parameters for Latitude GCM Grid Points
C-------------------------------------------------------------------

      IF(IFIRST.EQ.1) THEN
      NL0=NL
      NL1=NL+1
      IFIRST=0
      ENDIF

      DO 103 J=1,MLAT46
      DLATJ=DLAT46(J)
      ILATI=(DLATJ+95.001D0)/10.D0
      IF(ILATI.LT. 1) ILATI= 1
      IF(ILATI.GT.17) ILATI=17
      ILATIJ(J)=ILATI
      LATD=ILATI*10-95
      WTJL=(DLATJ-LATD)*0.1D0
      WTJLAT(J)=WTJL
      DO 101 JJ=2,14
      II=JJ-1
      IF(DLATJ.LE.DEGLAT(JJ)) GO TO 102
  101 CONTINUE
      JJ=14
  102 CONTINUE
      WTJJ=(DLATJ-DEGLAT(II))/(DEGLAT(JJ)-DEGLAT(II))
      WTII=1.-WTJJ
      CONCS(J)=WTII*CNCSPR(II)+WTJJ*CNCSPR(JJ)
      CONCA(J)=WTII*CNCAUT(II)+WTJJ*CNCAUT(JJ)
      BHKMS(J)=WTII*HKMSPR(II)+WTJJ*HKMSPR(JJ)
      BHKMA(J)=WTII*HKMAUT(II)+WTJJ*HKMAUT(JJ)
  103 CONTINUE

      DO 104 I=1,MLON72
      DLONI=DLON72(I)
      ILONG=DLONI/20.D0
      WTJLG=(DLONI-ILONG*20)/20.D0
      WTJLON(I)=WTJLG
      WTILG=1.-WTJLG
      ILONG=ILONG+1
      JLONG=ILONG+1
      IF(ILONG.GT.18) ILONG=18
      IF(ILONG.GT.17) JLONG=1
      ILONIJ(I)=ILONG
  104 CONTINUE
      DO 114 M=1,12
      DO 113 J=1,MLAT46
      ILATI=ILATIJ(J)
      WTJL=WTJLAT(J)
      WTIL=1.D0-WTJL
      JLATI=ILATI+1
      SUM72=0.D0
      DO 111 I=1,MLON72
      ILONG=ILONIJ(I)
      JLONG=ILONG+1
      IF(JLONG.GT.18) JLONG=1
      WTJLG=WTJLON(I)
      WTILG=1.D0-WTJLG
      AO3J=WTIL*(WTILG*O3AVE(M,ILATI,ILONG)
     +          +WTJLG*O3AVE(M,ILATI,JLONG))
     +    +WTJL*(WTILG*O3AVE(M,JLATI,ILONG)
     +          +WTJLG*O3AVE(M,JLATI,JLONG))
      SUM72=SUM72+AO3J
      D72(I)=AO3J
  111 CONTINUE
      SUM72=SUM72/72.D0
      DO 112 I=1,72
      O3NCAR(I,J,M)=D72(I)/SUM72
  112 CONTINUE
  113 CONTINUE
  114 CONTINUE
      RETURN

C-----------------
      ENTRY UPDO3L
C-----------------

      XJDAY=JDAY
      WTAUT=(XJDAY-91.D0)/213.D0
      IF(XJDAY.LT. 91.D0) WTAUT=( 91.D0-XJDAY)/152.D0
      IF(XJDAY.GT.304.D0) WTAUT=(456.D0-XJDAY)/152.D0
      WTSPR=1.D0-WTAUT
      DO 200 JMO=1,14
      XJDMJ=XJDMO(JMO)
      IF(XJDAY.LT.XJDMJ) GO TO 201
      XJDMI=XJDMJ
  200 CONTINUE
      XJDMI=XJDMO(13)
  201 CONTINUE
      DAYMO=XJDMJ-XJDMI
      WTJM=(XJDAY-XJDMI)/DAYMO
      WTIM=1.D0-WTJM
      JMO=JMO-1
      IMO=JMO-1
      IF(IMO.LT.1) IMO=12
      IF(JMO.GT.12) JMO=1
      JJDAY=1
      SJDAY=SOJDAY(JJDAY)
  202 CONTINUE
      JJDAY=JJDAY+1
      SIDAY=SJDAY
      SJDAY=SOJDAY(JJDAY)
      IF(XJDAY.GT.SJDAY) GO TO 202
      WTJAN=(XJDAY-SIDAY)/(SJDAY-SIDAY)
      IF(JJDAY.EQ.3.OR.JJDAY.EQ.5) WTJAN=1.-WTJAN
      WTSEP=1.D0-WTJAN
      DO 203 J=1,MLAT46
      DLATJ=DLAT46(J)
      DLSEP=10.D0+0.099999D0*DLATJ*PMLAT(JJDAY)
      DLJAN=10.D0+0.099999D0*DLATJ*PMLAT(JJDAY-1)
      LSEP=DLSEP
      LJAN=DLJAN
      LJANJ(J)=LJAN
      LSEPJ(J)=LSEP
      WTLSEP(J)=DLSEP-LSEP
      WTLJAN(J)=DLJAN-LJAN
  203 CONTINUE
!obso IF(KO3LON.EQ.1) RETURN
!obso IF(AO3J.GT.1.D-10) GO TO 400
      RETURN

C--------------------------------
      ENTRY GETLON(WO3LON,O3Z0AV)
C--------------------------------

      ILATI=ILATIJ(JLAT)
      WTJL=WTJLAT(JLAT)
      WTIL=1.D0-WTJL
      JLATI=ILATI+1
      LSEP=LSEPJ(JLAT)
      LJAN=LJANJ(JLAT)
      WTLS=WTLSEP(JLAT)
      WTLJ=WTLJAN(JLAT)
      AO3J=WTIM*(WTIL*AO3AVE(ILATI,IMO)+WTJL*AO3AVE(JLATI,IMO))
     +    +WTJM*(WTIL*AO3AVE(ILATI,JMO)+WTJL*AO3AVE(JLATI,JMO))
      BHKMJ=WTSPR*BHKMS(JLAT)+WTAUT*BHKMA(JLAT)
      CONCJ=WTSPR*CONCS(JLAT)+WTAUT*CONCA(JLAT)
      AO3JJ=AO3J

      SUMO=0.D0
      SUMN=0.D0
      DO 300 I=1,MLON72
      ILONG=ILONIJ(I)
      JLONG=ILONG+1
      IF(JLONG.GT.18) JLONG=1
      WTJLG=WTJLON(I)
      WTILG=1.D0-WTJLG
      AO3J=WTIM*(WTIL*(WTILG*O3AVE(IMO,ILATI,ILONG)
     +                +WTJLG*O3AVE(IMO,ILATI,JLONG))
     +          +WTJL*(WTILG*O3AVE(IMO,JLATI,ILONG)
     +                +WTJLG*O3AVE(IMO,JLATI,JLONG)))
     +    +WTJM*(WTIL*(WTILG*O3AVE(JMO,ILATI,ILONG)
     +                +WTJLG*O3AVE(JMO,ILATI,JLONG))
     +          +WTJL*(WTILG*O3AVE(JMO,JLATI,ILONG)
     +                +WTJLG*O3AVE(JMO,JLATI,JLONG)))
      SUMO=SUMO+AO3J
      SUMN=SUMN+1.D0
      AO3JIM(I)=AO3J
  300 CONTINUE
      O3MEAN=SUMO/SUMN
      AO3J=AO3JJ

      O3Z0AV=O3MEAN
      DO 310 I=1,MLON72
      WO3LON(I)=AO3JIM(I)/O3MEAN
  310 CONTINUE
!obso AO3J=AO3JJ+(AO3JIM(ILON)-AO3JJ)*FULGAS(13)
!obso IF(KO3LON.EQ.1) RETURN
      RETURN

!obso 400 CONTINUE
!obso     CKMJ=0.25D0*AO3J/CONCJ
!obso     GTOP=0.D0
!obso     POI=0.D0
!obso     FI=0.D0
!obso     L=NL0
!obso     PLL=PLB0(L)
!obso     J=12
!obso 401 CONTINUE
!obso     J=J-1
!obso     IF(J.LT.1) GO TO 404
!obso     POJ=PLBSO3(J)
!obso     FJ=WTSEP*(WTLS*SO3SO(J,LSEP+1)+(1.D0-WTLS)*SO3SO(J,LSEP))
!obso    +  +WTJAN*(WTLJ*SO3JF(J,LJAN+1)+(1.D0-WTLJ)*SO3JF(J,LJAN))
!obso 402 CONTINUE
!obso     DP=POJ-POI
!obso     IF(POJ.GT.PLL) GO TO 403
!obso     GTOP=GTOP+(FI+FJ)*DP*ACMMGG
!obso     POI=POJ
!obso     FI=FJ
!obso     GO TO 401
!obso 403 CONTINUE
!obso     FF=(FJ-FI)/DP
!obso     DP=PLL-POI
!obso     FF=FI+FF*DP
!obso     GTOP=GTOP+(FI+FF)*DP*ACMMGG
!obso     POI=PLL
!obso     FI=FF
!obso     O3LB(L)=GTOP
!obso     L=L-1
!obso     PLL=PLB0(L)
!obso     GO TO 402
!obso 404 CONTINUE
!obso     FI=FJ*ACMPKM
!obso     HI=H10MB
!obso     HJ=BHKMJ+CKMJ
!obso     XPBC=EXP(-BHKMJ/CKMJ)
!obso     XPHC=EXP(HJ/CKMJ)
!obso     DTERM=1.D0+XPHC*XPBC
!obso     ATERM=(1.D0+XPBC)/DTERM
!obso     FTERM=ATERM/DTERM*XPHC*XPBC/CKMJ
!obso     TTERM=AO3J-GTOP-FI*(HI-HJ)*0.5D0
!obso     AA=TTERM/(FTERM*(HI-HJ)*0.5D0+1.D0-ATERM)
!obso     FJ=AA*FTERM
!obso     GTOPBC=GTOP+(FI+FJ)*(HI-HJ)*0.5D0-AA*ATERM
!obso     TOP=AA*(1.D0+XPBC)
!obso     GO TO 406
!obso 405 CONTINUE
!obso     DH=HI-HJ
!obso     FF=(FJ-FI)/DH
!obso     DH=HI-H
!obso     FF=FI+FF*DH
!obso     GTOP=GTOP+(FI+FF)*DH*0.5D0
!obso     HI=H
!obso     FI=FF
!obso     O3LB(L)=GTOP
!obso     L=L-1
!obso 406 CONTINUE
!obso     H=HLB(L)
!obso     IF(H.GT.HJ) GO TO 405
!obso     O3LB(L)=TOP/(1.+XPBC*EXP(H/CKMJ))+GTOPBC
!obso     L=L-1
!obso     IF(L.GT.0) GO TO 406
!obso     O3LB(NL1)=0.D0
!obso     DO 407 L=1,NL0
!obso     U0GAS(L,3)=(O3LB(L)-O3LB(L+1))
!obso 407 CONTINUE
!obso     RETURN
      END SUBROUTINE SETO3L

      SUBROUTINE SETGAS

      !INCLUDE 'BR00B.COM'

      DIMENSION SINLAT(46)

C-----------------------------------------------------------------------
C     Global   U.S. (1976) Standard Atmosphere  P, T, Geo Ht  Parameters
C-----------------------------------------------------------------------

      DIMENSION P36(36),UFAC36(36)
      DATA P36/
     $1.2000E+03, .9720E+03, .9445E+03, .9065E+03, .8515E+03, .7645E+03,
     $ .6400E+03, .4975E+03, .3695E+03, .2795E+03, .2185E+03, .1710E+03,
     $ .1250E+03, .8500E+02, .6000E+02, .4000E+02, .2500E+02, .1500E+02,
     $ .7500E+01, .4000E+01, .2500E+01, .1500E+01, .7500E+00, .4000E+00,
     $ .2500E+00, .1500E+00, .7810E-01, .4390E-01, .2470E-01, .1390E-01,
     $ .7594E-02, .3623E-02, .1529E-02, .7030E-03, .2059E-03, .0E0/
      DATA UFAC36/
     $ 0.800,0.800,0.800,0.750,0.750,0.750,0.700,0.750,0.846,0.779,
     $ 0.892,0.886,0.881,0.875,0.870,0.846,0.840,0.902,0.880,0.775,
     $ 0.796,0.842,0.866,0.861,0.821,0.903,1.264,1.732,2.000,1.701,
     $ 1.609,1.478,1.253,1.372,1.571,1.571/

      DIMENSION SPLB(8),STLB(8),SHLB(8),SDLB(8)
      DATA SPLB/1013.25,226.32,54.748,8.6801,1.109,.66938,.039564
     +         ,3.7338E-03/
      DATA STLB/288.15,216.65,216.65,228.65,270.65,270.65,214.65,186.87/
      DATA SHLB/0.0,11.0,20.0,32.0,47.0,51.0,71.0,84.852/
      DATA SDLB/-6.5,0.0,1.0,2.8,0.0,-2.8,-2.0,0.0/

      SAVE P36,UFAC36,SPLB,STLB,SHLB,SDLB,SINLAT,NL0,NL1

      PARAMETER(HPCON=34.16319,P0=1013.25,PI=3.141592653589793D0)
      DATA IFIRST/1/

      SAVE  IFIRST

      IF(IFIRST.EQ.1) THEN
      DO 90 I=1,MLAT46
      RADLAT=DLAT46(I)*PI/180.D0
      SINLAT(I)=SIN(RADLAT)
   90 CONTINUE
      NL0=NL
      IFIRST=0
      ENDIF
C                  -----------------------------------------------------
C                  Use PLB to fix Standard Heights for Gas Distributions
C                  -----------------------------------------------------

      NL1=NL0+1
CCC   PS0=PLB0(1)

      DO 100 L=1,NL0
      DPL(L)=PLB0(L)-PLB0(L+1)
      PL(L)=(PLB0(L)+PLB0(L+1))*0.5D0
!obso HLB(L)=HLB0(L)
  100 CONTINUE
!obso HLB(NL1)=HLB0(NL1)
      CALL RETERP(UFAC36,P36,36,FPXCO2,PL,NL0)
cc    IUFAC=1
cc    IF(IUFAC.EQ.0) THEN
cc      DO I=1,NL
cc        FPXCO2(I)=1.
cc      END DO
cc    ENDIF

      NLAY=LASTVC/100000
      NATM=(LASTVC-NLAY*100000)/10000
      IF(NATM.GT.0) GO TO 112

C     ----------------------------------------------------------------
C     Define Default Global Mean Gas Amounts for Off-Line Use Purposes
C
C                                         Global Mean H2O Distribution
C                                         ----------------------------
      RHP=0.77D0
      EST=10.D0**(9.4051D0-2353.D0/TLB(1))
      FWB=0.662D0*RHP*EST/(PLB0(1)-RHP*EST)
      DO 111 L=1,NL0
      PLT=PLB0(L+1)
      DP=PLB0(L)-PLT
      RHP=0.77D0*(PLT/P0-0.02D0)/.98D0
      EST=10.D0**(9.4051D0-2353.D0/TLT(L))
      FWT=0.662D0*RHP*EST/(PLT-RHP*EST)
      IF(FWT.GT.3.D-06) GO TO 110
      FWT=3.D-06
      RHP=FWT*PLT/(EST*(FWT+0.662D0))
  110 CONTINUE
      ULGASL=0.5D0*(FWB+FWT)*DP*1268.75D0
      U0GAS(L,1)=ULGASL
      SHL(L)=ULGASL/(ULGASL+1268.75D0*DP)
      EQ=0.5D0*(PLB0(L)+PLT)*SHL(L)/(0.662D0+0.378D0*SHL(L))
      ES=10.D0**(9.4051D0-2353.D0/TLM(L))
      RHL(L)=EQ/ES
      FWB=FWT
  111 CONTINUE
  112 CONTINUE

C                                         ----------------------------
C                                         Global Mean NO2 Distribution
C                                         ----------------------------
      ACM=0.D0
      HI=0.D0
      FI=CMANO2(1)
      HL=HLB0(2)
      L=1
      J=1
  130 CONTINUE
      J=J+1
      IF(J.GT.42) GO TO 133
      HJ=HI+2.D0
      FJ=CMANO2(J)
  131 CONTINUE
      DH=HJ-HI
      IF(HJ.GT.HL) GO TO 132
      ACM=ACM+(FI+FJ)*DH*0.5D0
      HI=HJ
      FI=FJ
      GO TO 130
  132 CONTINUE
      FF=FI+(FJ-FI)*(HL-HI)/DH
      DH=HL-HI
      ACM=ACM+(FI+FJ)*DH*0.5D0
      U0GAS(L,5)=ACM
      ACM=0.D0
      HI=HL
      FI=FF
      IF(L.EQ.NL0) GO TO 133
      L=L+1
      HL=HLB0(L+1)
      GO TO 131
  133 CONTINUE
      U0GAS(L,5)=ACM
      ACM=0.D0
      L=L+1
      IF(L.LT.NL1) GO TO 133
C                            -----------------------------------------
C                            (CO2,O2) Uniformly Mixed Gas Distribution
C                            -----------------------------------------
      DO 141 K=2,4,2
      DO 140 N=1,NL0
      U0GAS(N,K)=PPMV80(K)*0.8D0*DPL(N)/P0
  140 CONTINUE
  141 CONTINUE
C                -----------------------------------------------------
C                (N20,CH4,F11,F12) Specified Vertical Gas Distribution
C                -----------------------------------------------------
      DO 151 K=6,12
      IF(K.EQ.10) GO TO 151
      DO 150 N=1,NL0
      GGVDF=1.D0-(1.D0-PPMVDF(K))*(1.D0-PLB0(N)/PLB0(1))
      IF(KGGVDF.LT.1) GGVDF=1.D0
      U0GAS(N,K)=PPMV80(K)*0.8D0*DPL(N)/P0*GGVDF
      ZT=(HLB0(N+1)-Z0(K))/ZH(K)
      IF(ZT.LE.0.D0) GO TO 150
      ZB=(HLB0(N)-Z0(K))/ZH(K)
      EXPZT=EXP(-ZT)
      EXPZB=EXP(-ZB)
      IF(ZB.LT.0.D0) EXPZB=1.D0-ZB
      U0GAS(N,K)=U0GAS(N,K)*(EXPZB-EXPZT)/(ZT-ZB)
  150 CONTINUE
  151 CONTINUE
C                         --------------------------------------------
C                         Specification of  FULGAS  Scaled Gas Amounts
C                         --------------------------------------------

      DO 230 L=1,NL0
Cc*** Adjust all water levels
cc    ULGAS(L,1)=U0GAS(L,1)*FULGAS(1)
c**** Only adjust stratospheric water levels (above LS1_loc)
      IF (L.lt.LS1_loc) THEN
        ULGAS(L,1)=U0GAS(L,1)
      ELSE
        ULGAS(L,1)=U0GAS(L,1)*FULGAS(1)
      END IF
C****
      ULGAS(L,3)=U0GAS(L,3)*FULGAS(3)
!obso ULGAS(L,5)=U0GAS(L,5)*FULGAS(5)
  230 CONTINUE

      IF(KPFOZO.EQ.1) THEN
      DO 235 L=1,NL0
      ULGAS(L,3)=ULGAS(L,3)*FPXOZO(L)
  235 CONTINUE
      ENDIF

      PARTTR=(PLB(1)-PTOPTR)/(PLB0(1)-PTOPTR)
      DO 240 L=1,NL0
      IF(PLB(L).LE.PTOPTR) PARTTR=1.D0
      DO 239 K=2,12
      IF(K.EQ.3) GO TO 239
      PARTTG=PARTTR
      IF(KPGRAD.GT.0) PARTTG=PARTTG*(1.D0+0.5D0*PPGRAD(K)*SINLAT(JLAT))
      ULGAS(L,K)=U0GAS(L,K)*FULGAS(K)*PARTTG
  239 CONTINUE
  240 CONTINUE

      IF(KPFCO2.EQ.1) THEN
      DO 245 L=1,NL0
      ULGAS(L,2)=ULGAS(L,2)*FPXCO2(L)
  245 CONTINUE
      ENDIF

      RETURN


C-----------------
      ENTRY GETGAS
C-----------------
C                        ---------------------------------------------
C                        Specify ULGAS: Get Gas Absorption from TAUGAS
C                        ---------------------------------------------

C                -----------------------------------------------------
C                N20,CH4,F11,F12 Specified Latitudinal Z0 Distribution
C                -----------------------------------------------------

      IF(KLATZ0.GT.0) THEN
      PTRO=100.D0
      DL=DLAT46(JLAT)
      DLS=-40.D0
      DLN= 40.D0
      IF(DL.LT.DLS) PTRO=189.D0-(DL+40.D0)*2.22D0
      IF(DL.GT.DLN) PTRO=189.D0+(DL-40.D0)*2.22D0
      DO 249 N=1,NL0
      IF(PLB0(N).GE.PTRO) Z0LAT=HLB0(N)
  249 CONTINUE
      DO 251 K=6,12
      IF(K.EQ.10) GO TO 251
      DO 250 N=1,NL0
      U0GAS(N,K)=PPMV80(K)*0.8D0*(PLB0(N)-PLB0(N+1))/P0
      ZT=(HLB0(N+1)-Z0LAT)/ZH(K)
      IF(ZT.LE.0.D0) GO TO 250
      ZB=(HLB0(N)-Z0LAT)/ZH(K)
      EXPZT=EXP(-ZT)
      EXPZB=EXP(-ZB)
      IF(ZB.LT.0.D0) EXPZB=1.D0-ZB
      U0GAS(N,K)=U0GAS(N,K)*(EXPZB-EXPZT)/(ZT-ZB)
  250 CONTINUE
  251 CONTINUE
      ENDIF

      DO 300 L=1,NL
      DPL(L)=PLB(L)-PLB(L+1)
      PL(L)=(PLB(L)+PLB(L+1))*0.5D0
  300 CONTINUE

      IF(KEEPRH.EQ.1) GO TO 311
      DO 310 L=1,NL
      U0GAS(L,1)=1268.75D0*DPL(L)*SHL(L)/(1.D0-SHL(L))
      EQ=PL(L)*SHL(L)/(0.662D0+0.378D0*SHL(L))
      ES=10.D0**(9.4051D0-2353.D0/TLM(L))
      RHL(L)=EQ/ES
  310 CONTINUE
      GO TO 313
  311 CONTINUE
      DO 312 L=1,NL
      ES=10.D0**(9.4051D0-2353.D0/TLM(L))
      SHL(L)=0.622D0*(RHL(L)*ES)/(PL(L)-0.378D0*(RHL(L)*ES))
      U0GAS(L,1)=1268.75D0*DPL(L)*SHL(L)/(1.D0-SHL(L))
  312 CONTINUE
  313 CONTINUE
      DO 314 L=1,NL
Cc*** Adjust ALL water levels
cc    ULGAS(L,1)=U0GAS(L,1)*FULGAS(1)
c**** Only adjust stratospheric water levels (above LS1_loc)
      IF (L.lt.LS1_loc) THEN
        ULGAS(L,1)=U0GAS(L,1)
      ELSE
        ULGAS(L,1)=U0GAS(L,1)*FULGAS(1)
      END IF
C****
  314 CONTINUE

      DO 330 L=1,NL0
      ULGAS(L,3)=U0GAS(L,3)*FULGAS(3)
!obso ULGAS(L,5)=U0GAS(L,5)*FULGAS(5)
  330 CONTINUE

      IF(KPFOZO.EQ.1) THEN
      DO 335 L=1,NL0
      ULGAS(L,3)=ULGAS(L,3)*FPXOZO(L)
  335 CONTINUE
      ENDIF

      PARTTR=(PLB(1)-PTOPTR)/(PLB0(1)-PTOPTR)
      DO 340 L=1,NL0
      IF(PLB(L).LE.PTOPTR) PARTTR=1.D0
      DO 339 K=2,12
      IF(K.EQ.3) GO TO 339
      PARTTG=PARTTR
      IF(KPGRAD.GT.0) PARTTG=PARTTG*(1.D0+0.5D0*PPGRAD(K)*SINLAT(JLAT))
      ULGAS(L,K)=U0GAS(L,K)*FULGAS(K)*PARTTG
  339 CONTINUE
  340 CONTINUE

      IF(KPFCO2.EQ.1) THEN
      DO 345 L=1,NL0
      ULGAS(L,2)=ULGAS(L,2)*FPXCO2(L)
  345 CONTINUE
      ENDIF

      IF(MRELAY.GT.0) THEN
      IF(NO3COL.GT.0) THEN
      SUMCOL=0.D0
      DO 410 L=1,NL0
      SUMCOL=SUMCOL+U0GAS(L,3)
  410 CONTINUE
      DO 420 L=1,NL0
      ULGAS(L,3)=U0GAS(L,3)*RO3COL/SUMCOL
  420 CONTINUE
      ENDIF
      DO 450 K=2,12
      IF(K.EQ.10.AND.KEEP10.GT.0) GO TO 450
      DO 430 L=1,NL0
      UGAS0(L)=ULGAS(L,K)
  430 CONTINUE
      CALL REPART(UGAS0,PLB0,NL1,UGASR,PLB,NLP)
      DO 440 L=1,NL
      ULGAS(L,K)=UGASR(L)
  440 CONTINUE
  450 CONTINUE
      IF(KEEP10.GT.0) THEN
      IF(KEEP10.LT.10) THEN
      DO 460 L=1,NL
      ULGAS(L,KEEP10)=ULGAS(L,10)
  460 CONTINUE
      ENDIF
      IF(KEEP10.GT.10) THEN
      DO 470 L=1,NL
      ULGAS(L,KEEP10-10)=ULGAS(L,KEEP10-10)+ULGAS(L,10)
  470 CONTINUE
      ENDIF
      ENDIF
      ENDIF

C-----------------
      CALL  TAUGAS
C-----------------

      RETURN
      END SUBROUTINE SETGAS

      SUBROUTINE SETO2A

      !INCLUDE 'BR00B.COM'

      DIMENSION SFWM2(18),SIGMA(18,6)
      DATA SFWM2/
     A 2.196E-03, 0.817E-03, 1.163E-03, 1.331E-03, 1.735E-03, 1.310E-03,
     B 1.311E-03, 2.584E-03, 2.864E-03, 4.162E-03, 5.044E-03, 6.922E-03,
     C 6.906E-03,10.454E-03, 5.710E-03, 6.910E-03,14.130E-03,18.080E-03/
      DATA SIGMA/
     A     2.74E-19, 2.74E-19, 2.74E-19, 2.74E-19, 2.74E-19, 2.74E-19,
     B     4.33E-21, 4.89E-21, 6.63E-21, 1.60E-20, 7.20E-20, 1.59E-18,
     C     2.10E-21, 2.32E-21, 3.02E-21, 6.30E-21, 3.46E-20, 7.52E-19,
     D     5.95E-22, 9.72E-22, 2.53E-21, 7.57E-21, 7.38E-20, 7.44E-19,
     E     3.33E-22, 1.02E-22, 4.09E-21, 1.63E-20, 8.79E-20, 3.81E-19,
     F     1.09E-21, 1.16E-21, 1.45E-21, 3.32E-21, 2.00E-20, 4.04E-19,
     G     1.15E-21, 1.30E-21, 1.90E-21, 4.89E-21, 2.62E-20, 4.08E-19,
     H     3.90E-22, 4.90E-22, 9.49E-22, 3.33E-21, 2.14E-20, 2.39E-19,
     I     1.29E-22, 2.18E-22, 8.28E-22, 3.46E-21, 1.94E-20, 1.06E-19,
     J     6.26E-23, 7.80E-23, 2.62E-22, 1.83E-21, 1.25E-20, 3.95E-20,
     K     2.74E-23, 3.58E-23, 8.64E-23, 4.03E-22, 2.13E-21, 1.95E-20,
     L     1.95E-23, 2.44E-23, 4.89E-23, 2.87E-22, 1.95E-21, 1.36E-20,
     M     1.84E-23, 1.96E-23, 2.71E-23, 8.52E-23, 6.48E-22, 3.89E-21,
     N     1.80E-23, 1.81E-23, 1.87E-23, 2.69E-23, 1.34E-22, 1.52E-21,
     O     1.80E-23, 1.80E-23, 1.82E-23, 2.40E-23, 5.71E-23, 5.70E-22,
     P     1.76E-23, 1.76E-23, 1.76E-23, 1.76E-23, 1.76E-23, 3.50E-23,
     Q     1.71E-23, 1.71E-23, 1.71E-23, 1.71E-23, 1.71E-23, 2.68E-23,
     R     1.00E-23, 1.00E-23, 1.00E-23, 1.00E-23, 1.00E-23, 1.00E-23/

      DIMENSION WTKO2(6)

      DATA WTKO2/0.05,0.20,0.25,0.25,0.20,0.05/

      SAVE  SFWM2,SIGMA,WTKO2,NL0,NL1,IFIRST

CCC   DATA STPMOL/2.68714D+19/,S00/1367.0/
CCC   DATA NW/18/,NZ/11/,NKO2/6/
      PARAMETER(STPMOL=2.68714D+19,S00=1367.0)
      PARAMETER(NW=18,NZ=11,NKO2=6)
      DATA IFIRST/1/

      IF(IFIRST.EQ.1) THEN
      NL0=NL
      DO 10 N=1,NL0
      ULGAS(N,4)=PPMV80(4)*0.8D0*(PLB0(N)-PLB0(N+1))/PLB0(1)
   10 CONTINUE
      IFIRST=0
      ENDIF

      NL1=NL0+1
      FSUM=0.0
      DO 100 I=1,NW
      FSUM=FSUM+SFWM2(I)
  100 CONTINUE
      DO 110 J=1,NZ
      ZTABLE(NL1,J)=FSUM
  110 CONTINUE
      SUMMOL=0.D0
      DO 150 N=1,NL0
      L=NL1-N
      SUMMOL=SUMMOL+ULGAS(L,4)*STPMOL
      DO 140 J=1,NZ
      ZCOS=0.01D0*(1/J)+0.1D0*(J-1)
      FSUM=0.D0
      DO 130 I=1,NW
      WSUM=0.D0
      DO 120 K=1,NKO2
      TAU=SIGMA(I,K)*SUMMOL/ZCOS
      IF(TAU.GT.30.D0) TAU=30.D0
      WSUM=WSUM+WTKO2(K)*EXP(-TAU)
  120 CONTINUE
      FSUM=FSUM+WSUM*SFWM2(I)
  130 CONTINUE
      ZTABLE(L,J)=FSUM
  140 CONTINUE
  150 CONTINUE
      DO 170 J=1,NZ
      DO 160 L=1,NL0
      DLFLUX=ZTABLE(L+1,J)-ZTABLE(L,J)
      ZTABLE(L,J)=DLFLUX/1366.D0
  160 CONTINUE
  170 CONTINUE
      RETURN

C-----------------
      ENTRY GETO2A
C-----------------

C              ---------------------------------------------------------
C              UV absorption by Oxygen is expressed as a fraction of the
C              total solar flux S0. Hence, O2FHRL(L)=ZTABLE(L,J) must be
C              normalized within SOLARM, dividing the GETO2A absorptions
C              O2FHRL(L) and O2FHRB(L) by the fraction of the solar flux
C              within the spectral interval DKS0(15), nominally by 0.05.
C              ---------------------------------------------------------

      ZCOS=1.D0+10.D0*COSZ
      JI=ZCOS
      IF(JI.GT.10) JI=10
      JJ=JI+1
      WTJ=ZCOS-JI
      WTI=1.0-WTJ
      DO 200 LL=1,NL
      L=NL1-LL
      IF(L.LT.1) GO TO 210
      O2FHRL(L)=(WTI*ZTABLE(L,JI)+WTJ*ZTABLE(L,JJ))
      O2FHRB(L)=ZTABLE(L,6)
  200 CONTINUE
  210 CONTINUE
      RETURN
      END SUBROUTINE SETO2A

      SUBROUTINE SETBAK

      !INCLUDE 'BR00B.COM'

C     ------------------------------------------------------------------
C     SETBAK,GETBAK  Initializes Background Aerosol Specification, i.e.,
C                    Aerosol Composition and Distribution that is set in
C                    BLOCK DATA RADPAR by AGOLDH, BGOLDH, CGOLDH Factors
C                    and controlled by FGOLDH ON/OFF Scaling Parameters.
C
C     ------------------------------------------------------------------
C     Tau Scaling Factors:    Solar    Thermal    apply to:
c                             FSTAER   FTTAER  ! Total Aerosol
c                             FSBAER   FTBAER  ! Bgrnd Aerosol
C
C     Control Parameters/Aerosol Scaling (kill) Factors
C                        FSTAER    SW   (All-type) Aerosol Optical Depth
C                        FTTAER    LW   (All-type) Aerosol Optical Depth
C                        FSBAER    SW   SETBAKonly Aerosol Optical Depth
C                        FTBAER    LW   SETBAKonly Aerosol Optical Depth
C                        -----------------------------------------------

      DIMENSION SGOLDH(5),TGOLDH(5)

      DATA IFIRST/1/,NL0/0/

      SAVE  IFIRST,NL0

C     ------------------------------------------------------------------
C     Thermal: Set (5) Aerosol Type Compositions & Vertical Distribution
C     ------------------------------------------------------------------

      IF(IFIRST.EQ.1) THEN
      NL0=NL
      IFIRST=0
      ENDIF

      DO 100 J=1,5
      DO 100 K=1,33
      DO 100 L=1,NL0
      TRAX(L,K,J)=0.D0
  100 CONTINUE

      DO 105 I=1,11
      DO 103 J=1,5
      IF(AGOLDH(I,J).LT.1.D-06) GO TO 103
      C=CGOLDH(I,J)
      BC=EXP(-BGOLDH(I,J)/C)
      ABC=AGOLDH(I,J)*(1.D0+BC)

      HXPB=1.D0
      DO 102 L=1,NL0
      HXPT=HLB0(L+1)/C
      IF(HXPT.GT.80.D0) GO TO 102
      HXPT=EXP(HXPT)
      ABCD=ABC/(1.D0+BC*HXPB)
     +    -ABC/(1.D0+BC*HXPT)
      HXPB=HXPT
      DO 101 K=1,33
      TRAX(L,K,J)=TRAX(L,K,J)+ABCD*(TRAQEX(K,I)-TRAQSC(K,I))
  101 CONTINUE
  102 CONTINUE
  103 CONTINUE
      DO 104 K=1,33
      TRAQAB(K,I)=TRAQEX(K,I)-TRAQSC(K,I)
  104 CONTINUE
  105 CONTINUE

      DO 107 K=1,33
      DO 106 L=1,NL0
      TRBALK(L,K)=0.D0
  106 CONTINUE
  107 CONTINUE

C-----------------------------------------------------------------------
C     SOLAR:   Set (5) Aerosol Type Compositions & Vertical Distribution
C-----------------------------------------------------------------------

      DO 110 J=1,5
      DO 110 K=1,6
      DO 110 L=1,NL0
      SRAX(L,K,J)=1.D-20
      SRAS(L,K,J)=1.D-30
      SRAC(L,K,J)=0.D0
  110 CONTINUE

      DO 114 I=1,11
      DO 113 J=1,5
      IF(AGOLDH(I,J).LT.1.D-06) GO TO 113
      C=CGOLDH(I,J)
      BC=EXP(-BGOLDH(I,J)/C)
      ABC=AGOLDH(I,J)*(1.D0+BC)

      HXPB=1.D0
      DO 112 L=1,NL0
      HXPT=HLB0(L+1)/C
      IF(HXPT.GT.80.D0) GO TO 112
      HXPT=EXP(HXPT)
      ABCD=ABC/(1.D0+BC*HXPB)
     +    -ABC/(1.D0+BC*HXPT)
      HXPB=HXPT
      DO 111 K=1,6
      SRAQX=SRAQEX(K,I)
      SRAQS=SRAQSC(K,I)
      SRAX(L,K,J)=SRAX(L,K,J)+ABCD*SRAQX
      SRAS(L,K,J)=SRAS(L,K,J)+ABCD*SRAQS
      SRAC(L,K,J)=SRAC(L,K,J)+ABCD*SRAQCB(K,I)*SRAQS
  111 CONTINUE
  112 CONTINUE
  113 CONTINUE
  114 CONTINUE

      DO 117 J=1,5
      DO 116 K=1,6
      DO 115 L=1,NL0
      SRAC(L,K,J)=SRAC(L,K,J)/SRAS(L,K,J)
  115 CONTINUE
  116 CONTINUE
  117 CONTINUE

      DO 119 K=1,6
      DO 118 L=1,NL0
      SRBEXT(L,K)=1.D-20
      SRBSCT(L,K)=0.D0
      SRBGCB(L,K)=0.D0
!obso SRBPI0(L,K)=0.D0
  118 CONTINUE
  119 CONTINUE

      RETURN

C-----------------
      ENTRY GETBAK
C-----------------

C     ------------------------------------------------------------------
C     GETBAK   Specifies Background Aerosol Contribution and Initializes
C                    (1) Thermal Radiation Aerosol Coefficient Table:
C                        TRAALK(L,K), for (L=1,NL), (K=1,33)
C
C                    (2) Solar Radiation Coefficient Tables:
C                        SRAEXT(L,K),SRASCT(L,K),SRAGCB(L,K) for (K=1,6)
C                    ---------------------------------------------------

C                                                              (Thermal)
C                                                              ---------
      DO 200 J=1,5
      TGOLDH(J)=FTTAER*FTBAER*FGOLDH(J)
  200 CONTINUE
      DO 203 L=1,NL0
      DO 202 K=1,33
      SUMABS=1.D-20
      DO 201 J=1,5
      SUMABS=SUMABS+TGOLDH(J)*TRAX(L,K,J)
  201 CONTINUE
      TRBALK(L,K)=SUMABS
  202 CONTINUE
  203 CONTINUE

C                                                                (Solar)
C                                                                -------

      DO 210 J=1,5
      SGOLDH(J)=FSTAER*FSBAER*FGOLDH(J)
  210 CONTINUE
      DO 212 K=1,6
      DO 212 L=1,NL0
      EXTSUM=1.D-20
      SCTSUM=1.D-30
      COSSUM=0.D0
      DO 211 J=1,5
      EXTSUM=EXTSUM+SGOLDH(J)*SRAX(L,K,J)
      SCTSUM=SCTSUM+SGOLDH(J)*SRAS(L,K,J)
      COSSUM=COSSUM+SGOLDH(J)*SRAS(L,K,J)*SRAC(L,K,J)
  211 CONTINUE
      SRBEXT(L,K)=EXTSUM
      SRBSCT(L,K)=SCTSUM
      SRBGCB(L,K)=COSSUM/SCTSUM
  212 CONTINUE


C     ------------------------------------------------------------------
C                            Option to repartition aerosol optical depth
C                                  for nonstandard GCM vertical layering
C                                  -------------------------------------
      IF(MRELAY.GT.0) THEN
      DO 240 K=1,6
      DO 220 L=1,NL0
      AERX1(L)=QAERO(L,K)
      AERS1(L)=SAERO(L,K)
      AERG1(L)=CAERO(L,K)*SAERO(L,K)
  220 CONTINUE
      CALL REPART(AERX1,PLB0,NL0,AERX2,PLB,NLP)
      CALL REPART(AERS1,PLB0,NL0,AERS2,PLB,NLP)
      CALL REPART(AERG1,PLB0,NL0,AERG2,PLB,NLP)


      DO 230 L=1,NL
      SRAEXT(L,K)=AERX2(L)
      SRASCT(L,K)=AERS2(L)
      SRAGCB(L,K)=AERG2(L)/(AERS2(L)+1.D-10)
  230 CONTINUE
  240 CONTINUE

      DO 270 K=1,33
      DO 250 L=1,NL0
      AERX1(L)=AAERO(L,K)
  250 CONTINUE
      CALL REPART(AERX1,PLB0,NL0,AERX2,PLB,NLP)
      DO 260 L=1,NL
      TRAALK(L,K)=AERX2(L)
  260 CONTINUE
  270 CONTINUE
      ENDIF

C     ------------------------------------------------------------------
C     Option to add on Tracer Type aerosol thermal & solar contributions
C     ------------------------------------------------------------------

      IF(NTRACE.GT.0) THEN

      DO 303 JJ=1,NTRACE
      J=5+JJ
      I=ITR(JJ)
C                                                              (Thermal)
C                                                              ---------
      DO 302 K=1,33
      SUMEXT=FGOLDH(J)*(TRAQEX(K,I)-TRAQSC(K,I))
      DO 301 L=1,NL
      TRBALK(L,K)=TRBALK(L,K)+SUMEXT*TRACER(L,JJ)
  301 CONTINUE
  302 CONTINUE
  303 CONTINUE

C                                                                (Solar)
C                                                                -------
      DO 305 K=1,6
      DO 305 L=1,NL
      EXTSUM=SRBEXT(L,K)
      SCTSUM=SRBSCT(L,K)
      COSSUM=SRBGCB(L,K)*SRBSCT(L,K)
      DO 304 JJ=1,NTRACE
      J=5+JJ
      I=ITR(JJ)
      SRAQX=SRAQEX(K,I)
      SRAQS=SRAQSC(K,I)
      IF(SRAQS.GT.SRAQX) SRAQS=SRAQX
      EXTSUM=EXTSUM+FGOLDH(J)*TRACER(L,JJ)*SRAQX
      SCTSUM=SCTSUM+FGOLDH(J)*TRACER(L,JJ)*SRAQS
      COSSUM=COSSUM+FGOLDH(J)*TRACER(L,JJ)*SRAQCB(K,I)*SRAQS
  304 CONTINUE
      SRBEXT(L,K)=EXTSUM
      SRBSCT(L,K)=SCTSUM
      SRBGCB(L,K)=COSSUM/SCTSUM
  305 CONTINUE
      ENDIF

      RETURN
      END SUBROUTINE SETBAK

      SUBROUTINE SETAER

      !INCLUDE 'BR00B.COM'

C     ---------------------------------------------------------------
C     GISS MONTHLY-MEAN AEROSOL (1950-1990) & DESERT DUST CLIMATOLOGY
C     ---------------------------------------------------------------
C
C     Modification of Aerosol Composition Amounts and/or of Radiative
C     Parameters is via BLOCK DATA RADPAR Parameters in COMMON/AERCOM
C
C     UPDAER INPUT: via ENTRY UPDAER Arguments:   IYEAR, IMONTH, IDAY
C                   specifies the time of Global Aerosol Distribution
C
C     GETAER INPUT: via ILON and JLAT set in radCOMMON for 72x46 grid
C
C           OUTPUT: via SRAEXT(L,K)  Aerosol Extinction Optical Depth
C                       SRASCT(L,K)  Aerosol Scattering Optical Depth
C                       SRAGCB(L,K)  Aerosol Asymmetry Parameter  g
C                       TRAALK(L,K)  Thermal Absorption Optical Depth
C
C     Defining Parameters are required for GETAER INPUT via radCOMMON
C
C                          NL  = Number of Model Layers
C                         NLP  = Number of Model Layers + 1
C                         PLB  = Layer Edge Pressure
C
C                      MEANAC  = 0  (Default BLOCK DATA RADPAR value)
C                      KVRAER  = 1  (Default BLOCK DATA RADPAR value)
C     REMARKS:
C
C     MEANAC=1  Over-rides IMONTH, IDAY to yield Annual-Mean Aerosols
C     KVRAER=0  Puts 12-Layer Aerosol Amounts in  first 12 PLB layers
C     KVRAER=1  Repartions Aerosol Vertical Distribution according to
C               PLB Pressure levels, i.e., from PLB(1) up to PLB(NLP)
C               (This takes Topography & Model Layering into account)
C     ---------------------------------------------------------------

C     ------------------------------------------------------------------
C     Tau Scaling Factors:    Solar    Thermal    apply to:
c                             FSTAER   FTTAER  ! Total Aerosol
c                             FSAAER   FTAAER  ! AClim Aerosol

C     Control Parameters/Aerosol Scaling (kill) Factors
C                        FSTAER    SW   (All-type) Aerosol Optical Depth
C                        FTTAER    LW   (All-type) Aerosol Optical Depth
C                        FSAAER    SW   AClim Aer  Aerosol Optical Depth
C                        FTAAER    LW   AClim Aer  Aerosol Optical Depth
C                        -----------------------------------------------
C
C            Aerosol Vertical Profile (Optical Depth) Scaling Factors
C            --------------------------------------------------------
C     Layr  1    2    3    4    5    6    7    8    9    10   11   12
C     Ptop 934  854  720  550  390  285  255  150  100   60   30   10
C     Pbot 984  934  854  720  550  390  285  255  150  100   60   30

      DIMENSION QXAERN(25),QSAERN(25),QGAERN(25)
      DIMENSION PLB12L(13),SRAX12(12),SRAS12(12),SRAG12(12),TRAX12(12)

      DATA PLB12L/
     A            984.0, 934.0, 854.0, 720.0, 550.0, 390.0, 285.0
     B           ,255.0, 150.0, 100.0,  60.0,  30.0,  10.0/
      SAVE PLB12L

C     ------------------------------------------------------------------
C     TROAER: Monthly-Mean Tropospheric Aerosols  (Column Optical Depth)
C                  Decade Maps (N=1,5) for BCI, OCI, SUI (K=1,3) aerosol
C                  for K=4 (N=1,5) gives SEA, SUN, OCN, OCB, BCB aerosol
C                  Map: IJ=72x46,  Month: M=1-12, AnnAv: M=13
C
C     VDBSCU: Zonal Mean Vertical Distribution for BCI,SUI-TYPE Aerosols
C
C                  Apply multiplicative scaling factors: for climatology
C                  aerosol optical depths, and single scattering albedos
C                  using (FSAERO,FTAERO,PI0MAX) factors defined in BLOCK
C                  DATA RADPAR in COMMON/AERCOM (contains defaults)
C
C        NOTE:     Scaling factors multiply the appropriate aerosol  MIE
C                  radiative parameters (QXAERO,QSAERO) and NOT the data
C                  of the column optical depth (72x46) global maps
C
C     IFIRST=0     Subsequent CALLs to SETAER  can reset scaling factors
C     ------------------------------------------------------------------
C                  Define aerosol size according to REAERO specification
C                  FRSULF= Sulfate fraction of basic aerosol composition
C
C     Set Size OCI (N=2) = Organic (Anthrop) Aerosol  (Nominal Reff=0.3)
C     Set Size SUI (N=3) = Sulfate (Anthrop) Aerosol  (Nominal Reff=1.0)
C     Set Size SEA (N=4) = SeaSalt (Natural) Aerosol  (Nominal Reff=2.0)
C     Set Size SUN (N=5) = Sulfate (Natural) Aerosol  (Nominal Reff=0.3)
C     Set Size ANT (N=6) = Nitrate (Anthrop) Aerosol  (Nominal Reff=1.0)
C     Set Size OCN (N=7) = Organic (Natural) Aerosol  (Nominal Reff=0.3)
C     Set Size OCB (N=8) = Organic (Anthrop) Aerosol  (Nominal Reff=0.3)
C     ------------------------------------------------------------------

      DO 115 NA=2,8
      N0=0
      IF(NA.EQ.2) N0=88
      IF(NA.EQ.4) N0=22
      IF(NA.EQ.6) N0=44
      IF(NA.EQ.7) N0=88
      IF(NA.EQ.8) N0=88
      AREFF=REAERO(NA)
      DO 112 K=1,6
      DO 111 N=1,22
      NN=N0+N
      WTS=FRSULF(NA)
      WTA=1.D0-WTS
      QXAERN(N)=SRUQEX(K,NN)*WTA+SRUQEX(K,N)*WTS
      QSAERN(N)=SRUQSC(K,NN)*WTA+SRUQSC(K,N)*WTS
      QGAERX=SRUQCB(K,NN)*SRUQSC(K,NN)*WTA+SRUQCB(K,N)*SRUQSC(K,N)*WTS
      QGAERN(N)=QGAERX/QSAERN(N)
  111 CONTINUE
      CALL SPLINE(REFU22,QXAERN,22,AREFF,SRBQEX(K,NA),1.D0,1.D0,1)
      CALL SPLINE(REFU22,QSAERN,22,AREFF,SRBQSC(K,NA),1.D0,1.D0,1)
      CALL SPLINE(REFU22,QGAERN,22,AREFF,SRBQCB(K,NA),1.D0,1.D0,1)
  112 CONTINUE
      DO 114 K=1,33
      DO 113 N=1,22
      NN=N0+N
      WTS=FRSULF(NA)
      WTA=1.D0-WTS
      QXAERN(N)=TRUQEX(K,NN)*WTA+TRUQEX(K,N)*WTS
      QSAERN(N)=TRUQSC(K,NN)*WTA+TRUQSC(K,N)*WTS
      QGAERX=TRUQCB(K,NN)*TRUQSC(K,NN)*WTA+TRUQCB(K,N)*TRUQSC(K,N)*WTS
      QGAERN(N)=QGAERX/(QSAERN(N)+1d-30)
  113 CONTINUE
      CALL SPLINE(REFU22,QXAERN,22,AREFF,TRBQEX(K,NA),1.D0,1.D0,1)
      CALL SPLINE(REFU22,QSAERN,22,AREFF,TRBQSC(K,NA),1.D0,1.D0,1)
      CALL SPLINE(REFU22,QGAERN,22,AREFF,TRBQCB(K,NA),1.D0,1.D0,1)
  114 CONTINUE
      CALL SPLINE(REFU22,Q55U22,22,AREFF,Q55B10(NA),1.D0,1.D0,1)
  115 CONTINUE
      REFB10(6)=AREFF
      VEFB10(6)=0.2D0


C     Set Size BCI (N=1) = Black Carbon (Industrial)  (Nominal Reff=0.1)
C     Set Size BCB (N=9) = Black Carbon (BioBurning)  (Nominal Reff=0.5)
C     ------------------------------------------------------------------

      DO 125 NNA=1,2
      NA=1
      IF(NNA.EQ.2) NA=9
      AREFF=REAERO(NA)
      DO 122 K=1,6
      DO 121 N=1,25
      QXAERN(N)=SRSQEX(K,N)
      QSAERN(N)=SRSQSC(K,N)
      QGAERN(N)=SRSQCB(K,N)
  121 CONTINUE
      CALL SPLINE(REFS25,QXAERN,25,AREFF,SRBQEX(K,NA),1.D0,1.D0,1)
      CALL SPLINE(REFS25,QSAERN,25,AREFF,SRBQSC(K,NA),1.D0,1.D0,1)
      CALL SPLINE(REFS25,QGAERN,25,AREFF,SRBQCB(K,NA),1.D0,1.D0,1)
  122 CONTINUE
      DO 124 K=1,33
      DO 123 N=1,25
      QXAERN(N)=TRSQEX(K,N)
      QSAERN(N)=TRSQSC(K,N)
      QGAERN(N)=TRSQCB(K,N)
  123 CONTINUE
      CALL SPLINE(REFS25,QXAERN,25,AREFF,TRBQEX(K,NA),1.D0,1.D0,1)
      CALL SPLINE(REFS25,QSAERN,25,AREFF,TRBQSC(K,NA),1.D0,1.D0,1)
      CALL SPLINE(REFS25,QGAERN,25,AREFF,TRBQCB(K,NA),1.D0,1.D0,1)
  124 CONTINUE
      CALL SPLINE(REFS25,Q55U22,25,AREFF,Q55B10(NA),1.D0,1.D0,1)
  125 CONTINUE

      DO 150 N=1,10
      DO 130 I=1,6
      PI0F=1.D0
      PI0=SRBQSC(I,N)/SRBQEX(I,N)
      IF(PI0.GT.PI0MAX(N)) PI0F=PI0MAX(N)/PI0
      QXAERO(I,N)=SRBQEX(I,N)
      QSAERO(I,N)=SRBQSC(I,N)*PI0F
      QCAERO(I,N)=SRBQCB(I,N)
  130 CONTINUE
      PIAERO(N)=QSAERO(6,N)/(QXAERO(6,N)+1.D-10)
      DO 140 K=1,33
      ATAERO(K,N)=TRBQEX(K,N)-TRBQSC(K,N)
      TRBQAB(K,N)=TRBQEX(K,N)-TRBQSC(K,N)
  140 CONTINUE
c     IF(N.EQ.3) CALL SETREL(N)
  150 CONTINUE

      RETURN

C--------------------------------
      ENTRY UPDAER(JYEARA,JJDAYA)
C--------------------------------

C                    Aerosol Climatology Time Dependence Selection
C                    ---------------------------------------------
C
C          Time Selection is by UPDAER(IYEAR,IMONTH,IDAY) Argument
C          Integer Year, Month, Day values are converted to XJYEAR
C
C     IF(MEANAC.EQ.0)     Time Dependence is Linearly interpolated
C                         between Monthly-mean Aerosol Data Tables
C
C     IF(MEANAC.EQ.1)     Annual-mean  (72X46) Climatology is used
C
C     MEANAC=0 Default value set in BLOCK DATA RADPAR
C     ------------------------------------------------------------

      XJYEAR=JYEARA+(JJDAYA-0.999D0)/366.D0

      XJYBCI=XJYEAR
      IF(X0YBCI.GT.1.D0) XJYBCI=X0YBCI
      XJYOCI=XJYEAR
      IF(X0YOCI.GT.1.D0) XJYOCI=X0YOCI
      XJYSUI=XJYEAR
      IF(X0YSUI.GT.1.D0) XJYSUI=X0YSUI

C                         ----------------------------------------
C                         Decadal Time Interpolation for BC and SU
C                         Tables from 1950 (IBC=1) to 1990 (JBC=5)
C                         BCWTID, BCWTJD Weights need not sum to 1
C                         Exponential Extrapolation: <1950, 1990>
C                         ----------------------------------------

      CALL BCTAUW(XJYBCI,IBCI,JBCI,BCIWID,BCIWJD)
      CALL BCTAUW(XJYOCI,IOCI,JOCI,OCIWID,OCIWJD)
      CALL SUTAUW(XJYSUI,ISUI,JSUI,SUIWID,SUIWJD)
      XJYEAR=XJYEAR+1.D0/24.D0
      JJYEAR=XJYEAR
      DJYEAR=XJYEAR-JJYEAR
      XMO=DJYEAR*12.D0
      MA=XMO
      MB=MA+1
      WMB=XMO-MA
      WMA=1.D0-WMB
      IF(MA.LT.1)  MA=12
      IF(MB.GT.12) MB=1
      IF(MEANAC.EQ.1) THEN
      WMA=1.D0
      WMB=0.D0
      MA=13
      MB=13
      ENDIF

      DO 220 J=1,46
      DO 210 I=1,72
      TAUBCI=WMA*TROAER(I,J,MA,IBCI,1)+WMB*TROAER(I,J,MB,IBCI,1)
      TAUBCJ=WMA*TROAER(I,J,MA,JBCI,1)+WMB*TROAER(I,J,MB,JBCI,1)
      TAUCOL(I,J,1)=BCIWID*TAUBCI+BCIWJD*TAUBCJ
      TAUOCI=WMA*TROAER(I,J,MA,IOCI,2)+WMB*TROAER(I,J,MB,IOCI,2)
      TAUOCJ=WMA*TROAER(I,J,MA,JOCI,2)+WMB*TROAER(I,J,MB,JOCI,2)
      TAUCOL(I,J,2)=OCIWID*TAUOCI+OCIWJD*TAUOCJ
      TAUSUI=WMA*TROAER(I,J,MA,ISUI,3)+WMB*TROAER(I,J,MB,ISUI,3)
      TAUSUJ=WMA*TROAER(I,J,MA,JSUI,3)+WMB*TROAER(I,J,MB,JSUI,3)
      TAUCOL(I,J,3)=SUIWID*TAUSUI+SUIWJD*TAUSUJ
      TAUCOL(I,J,4)=WMA*TROAER(I,J,MA,1,4)+WMB*TROAER(I,J,MB,1,4)
      TAUCOL(I,J,5)=WMA*TROAER(I,J,MA,2,4)+WMB*TROAER(I,J,MB,2,4)
      TAUCOL(I,J,6)=WMA*TROAER(I,J,MA,2,4)+WMB*TROAER(I,J,MB,2,4)
      TAUCOL(I,J,7)=WMA*TROAER(I,J,MA,3,4)+WMB*TROAER(I,J,MB,3,4)
      TAUCOL(I,J,8)=WMA*TROAER(I,J,MA,4,4)+WMB*TROAER(I,J,MB,4,4)
      TAUCOL(I,J,9)=WMA*TROAER(I,J,MA,5,4)+WMB*TROAER(I,J,MB,5,4)
      TAUCOL(I,J,10)=SSBTAU
  210 CONTINUE
  220 CONTINUE

      DO 240 L=1,12
      DO 230 J=1,46
      VDAERO(J,L,1)=WMA*VDBCSU(J,L,MA,1)+WMB*VDBCSU(J,L,MB,1)
      VDAERO(J,L,2)=WMA*VDBCSU(J,L,MA,2)+WMB*VDBCSU(J,L,MB,2)
      VDAERO(J,L,3)=WMA*VDBCSU(J,L,MA,3)+WMB*VDBCSU(J,L,MB,3)
  230 CONTINUE
  240 CONTINUE
      RETURN


C-----------------
      ENTRY GETAER
C-----------------

C     ------------------------------------------------------------
C     Tropospheric Aerosol and Mineral Dust Data are on 72X46 Grid
C     with a Nominal 5 deg Longitude and 4 deg Latitude Resolution
C     with 12 Layer Vertical Distribution of Aerosol Optical Depth
C
C     Geographic Location is specified by: ILON and JLAT in RADCOM
C     ------------------------------------------------------------
C
C                  -----------------------------------------------
C                  For Model Resolution other than 72X36, Redefine
C                  ILON and JLAT appropriately
C
C     ------------------------------------------------------------
C     Aerosol Vertical Distribution  vs. Model Vertical Resolution
C
C     IF(KVRAER.EQ.0) 12 Layer Vertical Distribution is used as is
C
C     IF(KVRAER.EQ.1) Aerosol Height Distribution is Repartitioned
C                             Using PLB in RADCOM
C                     (This takes Surface Topography into account)
C     ------------------------------------------------------------


C                                                (Solar Component)
C                                                -----------------
      DO 330 K=1,6
      DO 320 L=1,12
      VDTAU=SSBTAU*VDGAER(L,10)*FSAERO(10)
      QSUM=VDTAU*QXAERO(K,10)+2.D-10
      SSUM=VDTAU*QSAERO(K,10)+1.D-10
      CSUM=SSUM*QCAERO(K,10)
      DO 310 N=1,3
      VDTAU=TAUCOL(ILON,JLAT,N)*VDAERO(JLAT,L,N)*VDFAER(L,N)*FSAERO(N)
      QSUM=QSUM+VDTAU*QXAERO(K,N)
      VDQS=VDTAU*QSAERO(K,N)
      SSUM=SSUM+VDQS
      CSUM=CSUM+VDQS*QCAERO(K,N)
  310 CONTINUE
      QAERO(L,K)=QSUM
      SAERO(L,K)=SSUM
      CAERO(L,K)=CSUM/SSUM
  320 CONTINUE
  330 CONTINUE

      DO 360 K=1,6
      DO 350 L=1,6
      QSUM=QAERO(L,K)
      SSUM=SAERO(L,K)
      CSUM=CAERO(L,K)*SSUM
      DO 340 N=4,9
      VDTAU=TAUCOL(ILON,JLAT,N)*VDGAER(L,N)*FSAERO(N)
      QSUM=QSUM+VDTAU*QXAERO(K,N)
      VDQS=VDTAU*QSAERO(K,N)
      SSUM=SSUM+VDQS
      CSUM=CSUM+VDQS*QCAERO(K,N)
  340 CONTINUE
      QAERO(L,K)=QSUM
      SAERO(L,K)=SSUM
      CAERO(L,K)=CSUM/SSUM
  350 CONTINUE
  360 CONTINUE

C                                              (Thermal Component)
C                                              -------------------
      DO 430 K=1,33
      DO 420 L=1,12
      VDTAU=SSBTAU*VDGAER(L,10)*FTAERO(10)
      ASUM=VDTAU*ATAERO(K,10)
      DO 410 N=1,3
      VDTAU=TAUCOL(ILON,JLAT,N)*VDAERO(JLAT,L,N)*VDFAER(L,N)*FTAERO(N)
      ASUM=ASUM+VDTAU*ATAERO(K,N)
  410 CONTINUE
      AAERO(L,K)=ASUM
  420 CONTINUE
  430 CONTINUE

      DO 460 K=1,33
      DO 450 L=1,6
      ASUM=AAERO(L,K)
      DO 440 N=4,9
      VDTAU=TAUCOL(ILON,JLAT,N)*VDGAER(L,N)*FTAERO(N)
      ASUM=ASUM+VDTAU*ATAERO(K,N)
  440 CONTINUE
      AAERO(L,K)=ASUM
  450 CONTINUE
  460 CONTINUE

C                     -------------------------------------------------
C                     Define Aerosol Thermal and Solar Contributions in
C                                       TRAALK and SRAEXT,SRASCT,SRAGCB
C
C                     Apply Solar/Thermal Optical Depth Scaling Factors
C                              AClim Aerosol  Solar  FSXA=FSTAER*FSAAER
C                              AClim Aerosol Thermal FTXA=FSTAER*FTAAER
C                              ----------------------------------------

      FSXTAU=FSTAER*FSAAER+1.D-10
      FTXTAU=FTTAER*FTAAER
      IF(KVRAER.EQ.1) THEN
C                           -------------------------------------------
C                           Option to repartition aerosol optical depth
C                                 for nonstandard GCM vertical layering
C                                 -------------------------------------
      DO 530 K=1,6
      DO 510 L=1,12
      SRAX12(L)=QAERO(L,K)*FSXTAU
      SRAS12(L)=SAERO(L,K)*FSXTAU
      SRAG12(L)=CAERO(L,K)*SAERO(L,K)*FSXTAU
  510 CONTINUE
      CALL REPART(SRAX12,PLB12L,13,SRAXNL,PLB,NLP)
      CALL REPART(SRAS12,PLB12L,13,SRASNL,PLB,NLP)
      CALL REPART(SRAG12,PLB12L,13,SRAGNL,PLB,NLP)
      DO 520 L=1,NL
      SRAEXT(L,K)=SRAXNL(L)
      SRASCT(L,K)=SRASNL(L)
      SRAGCB(L,K)=SRAGNL(L)/(SRASNL(L)+1.D-10)
  520 CONTINUE
  530 CONTINUE

      DO 560 K=1,33
      DO 540 L=1,12
      TRAX12(L)=AAERO(L,K)*FTXTAU
  540 CONTINUE
      CALL REPART(TRAX12,PLB12L,13,TRAXNL,PLB,NLP)
      DO 550 L=1,NL
      TRAALK(L,K)=TRAXNL(L)
  550 CONTINUE
  560 CONTINUE

      ELSE
      DO 620 K=1,6
      DO 610 L=1,NL
      SRAEXT(L,K)=QAERO(L,K)*FSXTAU
      SRASCT(L,K)=SAERO(L,K)*FSXTAU
      SRAGCB(L,K)=CAERO(L,K)
      IF(L.LT.13) GO TO 610
      SRAEXT(L,K)=2.D-10
      SRASCT(L,K)=1.D-10
      SRAGCB(L,K)=0.D0
  610 CONTINUE
  620 CONTINUE

      DO 640 K=1,33
      DO 630 L=1,NL
      TRAALK(L,K)=AAERO(L,K)*FTXTAU
      IF(L.LT.13) GO TO 630
      TRAALK(L,K)=0.D0
  630 CONTINUE
  640 CONTINUE
      ENDIF

      RETURN
      END SUBROUTINE SETAER

      SUBROUTINE SETDST

      !INCLUDE 'BR00B.COM'

C     ---------------------------------------------------------------
C     GISS MONTHLY-MEAN AEROSOL (1950-1990) & DESERT DUST CLIMATOLOGY
C     ---------------------------------------------------------------
C
C     Modification of Aerosol Composition Amounts and/or of Radiative
C     Parameters is via BLOCK DATA RADPAR Parameters in COMMON/AERCOM
C
C     UPDDST INPUT: via ENTRY UPDDST Arguments:   IYEAR, IMONTH, IDAY
C                     specifies the time of Dust Aerosol Distribution
C
C     GETDST INPUT: via ILON and JLAT set in radCOMMON for 72x46 grid
C
C           OUTPUT: via SRDEXT(L,K)   D Dust Extinction Optical Depth
C                       SRDSCT(L,K)   D Dust Scattering Optical Depth
C                       SRDGCB(L,K)   D Dust Asymmetry Parameter  g
C                       TRDALK(L,K)  Thermal Absorption Optical Depth
C
C     Defining Parameters are required for GETdst INPUT via radCOMMON
C
C                          NL  = Number of Model Layers
C                         NLP  = Number of Model Layers + 1
C                         PLB  = Layer Edge Pressure
C
C                      MEANDD  = 0  (Default BLOCK DATA RADPAR value)
C                      KVRAER  = 1  (Default BLOCK DATA RADPAR value)
C     REMARKS:
C
C     MEANDD=1  Over-rides IMONTH, IDAY to yield Annual-Mean Aerosols
C     KVRAER=0  Puts 12-Layer Aerosol Amounts in  first 12 PLB layers
C     KVRAER=1  Repartions Aerosol Vertical Distribution according to
C               PLB Pressure levels, i.e., from PLB(1) up to PLB(NLP)
C               (This takes Topography & Model Layering into account)
C     ---------------------------------------------------------------

C-----------------------------------------------------------------------
C     Tau Scaling Factors:    Solar    Thermal    apply to:
c                             FSTAER   FTTAER  ! Total Aerosol
c                             FSDAER   FTDAER  ! Dust  Aerosol
C
C     Control Parameters/Aerosol Scaling (kill) Factors
C                        FSTAER    SW   (All-type) Aerosol Optical Depth
C                        FTTAER    LW   (All-type) Aerosol Optical Depth
C                        FSDAER    SW   Dust Aer   Aerosol Optical Depth
C                        FTDAER    LW   Dust Aer   Aerosol Optical Depth
C                        -----------------------------------------------


C       Dust Aerosol Vertical Profile (Optical Depth) Scaling Factors
C       -------------------------------------------------------------
C     Layr  1    2    3    4    5    6    7    8    9    10   11   12
C     Ptop 934  854  720  550  390  285  255  150  100   60   30   10
C     Pbot 984  934  854  720  550  390  285  255  150  100   60   30

      DIMENSION PLB12L(13),SRAX12(12),SRAS12(12),SRAG12(12),TRAX12(12)
      DATA PLB12L/
     A            984.0, 934.0, 854.0, 720.0, 550.0, 390.0, 285.0
     B           ,255.0, 150.0, 100.0,  60.0,  30.0,  10.0/

      DIMENSION FRACD1(9),FRACD2(9),FRACDL(9),LFRACD(9)
CCC   DIMENSION FRACD1(9),FRACD2(9),FRACDL(9),LFRACD(9),TAUCON(8)

      DATA FRACD1/105.0, 30.0, 60.0, 50.0, 30.0, 30.0, 20.0, 3.0, 1.0/
      DATA FRACD2/135.0,135.0,105.0, 80.0, 80.0, 60.0, 60.0, 4.0, 6.0/
      DATA LFRACD/  6,    6,    7,    8,    8,    9,    9,    8,  12/

      DATA IFIRST/1/
      SAVE  IFIRST,PLB12L,FRACD1,FRACD2,FRACDL,LFRACD

C     ------------------------------------------------------------------
C     DUST:   Monthly-Mean Desert Dust (Clay,Silt) 8-Size Optical Depths
C                   Map: IJ=72x46, Lay: L=1-9, Siz: S=1-8, Month: M=1-12
C
C     IFIRST=1      Read in Dust Data: (Dust Data are read in by RCOMP1)
C                   Apply multiplicative scaling factors for Desert Dust
C                   optical depths, and single scattering albedos, using
C                   factors FSDUST, FTDUST  defined in BLOCK DATA RADPAR
C                   COMMON/AERCOM (contains default values)
C
C        NOTE:      Scaling factors multiply the appropriate Desert Dust
C                   MIE radiative parameters (QXDUST,QSDUST) and NOT the
C                   column optical depth data of the (72x46) global maps
C
C     IFIRST=0      Subsequent CALLs to SETAER can reset scaling factors
C     ------------------------------------------------------------------


      IF(IFIRST.EQ.1) THEN

      DO 110 L=1,9
      FRACDL(L)=FRACD1(L)/FRACD2(L)
  110 CONTINUE

      DO 120 N=1,8
      FNDUST=FSDUST(N)
      CALL SDDUST
     +(QXDUST(1,N),QSDUST(1,N),QCDUST(1,N),QDST55(N),REDUST(N),FNDUST)
      CALL TDDUST(ATDUST(1,N),REDUST(N),FTDUST(N))
CCC   TAUCON(N)=0.75E+03*QDST55(N)/(RODUST(N)*REDUST(N))
      PIDUST(N)=QSDUST(6,N)/(QXDUST(6,N)+1.D-10)
  120 CONTINUE

      M=13
      DO 160 N=1,8
      DO 150 J=1,46
      DO 140 I=1,72
      DSUM=0.D0
      DO 130 L=1,9
      DSUM=DSUM+TDUST(I,J,L,N,M)
  130 CONTINUE
      DSTCOL(I,J,N)=DSUM
  140 CONTINUE
  150 CONTINUE
  160 CONTINUE

      DO 180 I=1,25
      DO 170 K=1,33
      TRDQAB(K,I)=TRDQEX(K,I)-TRDQSC(K,I)
  170 CONTINUE
  180 CONTINUE
      IFIRST=0
      RETURN
      ENDIF

      DO 190 N=1,8
      FNDUST=FSDUST(N)
      CALL SDDUST
     +(QXDUST(1,N),QSDUST(1,N),QCDUST(1,N),QDST55(N),REDUST(N),FNDUST)
      CALL TDDUST(ATDUST(1,N),REDUST(N),FTDUST(N))
      PIDUST(N)=QSDUST(6,N)/(QXDUST(6,N)+1.D-10)
  190 CONTINUE

      RETURN


C--------------------------------
      ENTRY UPDDST(JYEARD,JJDAYD)
C--------------------------------

      XJYEAR=JYEARD+(JJDAYD-0.999D0)/366.D0

      XJYEAR=XJYEAR+1.D0/24.D0
      JJYEAR=XJYEAR
      DJYEAR=XJYEAR-JJYEAR
      XMO=DJYEAR*12.D0
      MA=XMO
      MB=MA+1
      WMB=XMO-MA
      WMA=1.D0-WMB
      IF(MA.LT.1)  MA=12
      IF(MB.GT.12) MB=1
      IF(MEANDD.EQ.1) THEN
      WMA=1.D0
      WMB=0.D0
      MA=13
      MB=13
      ENDIF

      DO 290 N=1,8
      DO 280 J=1,46
      DO 270 I=1,72
      DSUM=0.D0
      DO 250 L=1,5
      DUSTLN(I,J,L,N)=WMA*TDUST(I,J,L,N,MA)+WMB*TDUST(I,J,L,N,MB)
      DSUM=DSUM+DUSTLN(I,J,L,N)
  250 CONTINUE
      DO 260 L=6,12
      K=L-5
      M=LFRACD(K)
      DUSTTL=WMA*TDUST(I,J,M,N,MA)+WMB*TDUST(I,J,M,N,MB)
      DUSTLN(I,J,L,N)=DUSTTL*FRACDL(K)
      DSUM=DSUM+DUSTLN(I,J,L,N)
  260 CONTINUE
      DUSTLN(I,J,10,N)=DUSTLN(I,J,10,N)+DUSTTL*FRACDL(9)
      DUSTLN(I,J, 7,N)=DUSTLN(I,J, 7,N)+DUSTLN(I,J, 8,N)*FRACDL(8)
      DSUM=DSUM+DUSTLN(I,J,8,N)*FRACDL(8)+DUSTTL*FRACDL(9)
      DSTCOL(I,J,N)=DSUM
  270 CONTINUE
  280 CONTINUE
  290 CONTINUE
      RETURN


C-----------------
      ENTRY GETDST
C-----------------

C     ------------------------------------------------------------
C     Wind-blown Mineral Dust data are represented on a 72X46 Grid
C     with a nominal 5 deg longitude and 4 deg latitude resolution
C     with 12 layer vertical distribution of aerosol optical depth
C
C     Geographic Location is specified by: ILON and JLAT in RADCOM
C
C                  For Model Resolution other than 72X36, Redefine
C                  ILON and JLAT appropriately
C
C     Aerosol Vertical Distribution  vs. Model Vertical Resolution
C
C     IF(KVRAER.EQ.0) 12 Layer Vertical Distribution is used as is
C
C     IF(KVRAER.EQ.1) Aerosol Height Distribution is Repartitioned
C                             Using PLB in RADCOM
C                     (This takes Surface Topography into account)
C     ------------------------------------------------------------


C                                                (Solar Component)
C                                                -----------------

      DO 390 L=1,12
      VDF=VDFDST(L)
      DO 380 K=1,6
      QSUM=0.D0
      SSUM=1.D-10
      CSUM=0.D0
      DO 370 N=1,8
      QSUM=QSUM+QXDUST(K,N)*DUSTLN(ILON,JLAT,L,N)
      SRDN=     QSDUST(K,N)*DUSTLN(ILON,JLAT,L,N)
      SSUM=SSUM+SRDN
      CSUM=CSUM+QCDUST(K,N)*SRDN
  370 CONTINUE
      CDUST(L,K)=CSUM/SSUM
      QDUST(L,K)=QSUM*VDF
      SDUST(L,K)=SSUM*VDF
  380 CONTINUE
  390 CONTINUE

C                                              (Thermal Component)
C                                              -------------------

      DO 490 L=1,12
      VDF=VDFDST(L)
      DO 480 K=1,33
      ASUM=0.D0
      DO 470 N=1,8
      ASUM=ASUM+ATDUST(K,N)*DUSTLN(ILON,JLAT,L,N)
  470 CONTINUE
      ADUST(L,K)=ASUM*VDF
  480 CONTINUE
  490 CONTINUE

C                ------------------------------------------------------
C                Define Mineral Dust Thermal and Solar Contributions to
C                         TRDALK and SRDEXT,SRDSCT,SRDGCB, respectively
C
C                     Apply Solar/Thermal Optical Depth Scaling Factors
C                              Dust Aerosol  Solar   FSXD=FSTAER*FSDAER
C                              Dust Aerosol Thermal  FTXD=FSTAER*FTDAER
C                              ----------------------------------------

      FSXTAU=FSTAER*FSDAER
      FTXTAU=FTTAER*FTDAER
      IF(KVRAER.EQ.1) THEN
C                           -------------------------------------------
C                           Option to repartition aerosol optical depth
C                                 for nonstandard GCM vertical layering
C                                 -------------------------------------
      DO 530 K=1,6
      DO 510 L=1,12
      SRAX12(L)=QDUST(L,K)*FSXTAU
      SRAS12(L)=SDUST(L,K)*FSXTAU
      SRAG12(L)=CDUST(L,K)*SDUST(L,K)*FSXTAU
  510 CONTINUE
      CALL REPART(SRAX12,PLB12L,13,SRAXNL,PLB,NLP)
      CALL REPART(SRAS12,PLB12L,13,SRASNL,PLB,NLP)
      CALL REPART(SRAG12,PLB12L,13,SRAGNL,PLB,NLP)
      DO 520 L=1,NL
      SRDEXT(L,K)=SRAXNL(L)
      SRDSCT(L,K)=SRASNL(L)
      SRDGCB(L,K)=SRAGNL(L)/(SRASNL(L)+1.D-10)
  520 CONTINUE
  530 CONTINUE

      DO 560 K=1,33
      DO 540 L=1,12
      TRAX12(L)=ADUST(L,K)*FTXTAU
  540 CONTINUE
      CALL REPART(TRAX12,PLB12L,13,TRAXNL,PLB,NLP)
      DO 550 L=1,NL
      TRDALK(L,K)=TRAXNL(L)
  550 CONTINUE
  560 CONTINUE

      ELSE
      DO 620 K=1,6
      DO 610 L=1,NL
      SRDEXT(L,K)=QDUST(L,K)*FSXTAU
      SRDSCT(L,K)=SDUST(L,K)*FSXTAU
      SRDGCB(L,K)=CDUST(L,K)
      IF(L.LT.13) GO TO 610
      SRDEXT(L,K)=2.D-10
      SRDSCT(L,K)=1.D-10
      SRDGCB(L,K)=0.D0
  610 CONTINUE
  620 CONTINUE

      DO 640 K=1,33
      DO 630 L=1,NL
      TRDALK(L,K)=ADUST(L,K)*FTXTAU
      IF(L.LT.13) GO TO 630
      TRDALK(L,K)=0.D0
  630 CONTINUE
  640 CONTINUE
      ENDIF

      RETURN
      END SUBROUTINE SETDST

      SUBROUTINE SETVOL

      !INCLUDE 'BR00B.COM'

      DIMENSION E24LAT(25),EJMLAT(47)
      DIMENSION HLATTF(4),HLATKM(5)

      DATA HLATKM/ 15.0, 20.0, 25.0, 30.0, 35.0/,LATVOL/0/

      parameter (htplim=1.d-3)
      SAVE E24LAT,EJMLAT,HLATKM,FSXTAU,FTXTAU,NJ25,NJJM,LATVOL

C     ------------------------------------------------------------------
C     Tau Scaling Factors:    Solar    Thermal    apply to:
c                             FSTAER   FTTAER  ! Total Aerosol
c                             FSVAER   FTVAER  ! SETVOL Aer

C     Control Parameters/Aerosol Scaling (kill) Factors
C                        FSTAER    SW  (All-type) Aerosol Optical Depth
C                        FTTAER    LW  (All-type) Aerosol Optical Depth
C                        FSVAER    SW  SETVOLonly Aerosol Optical Depth
C                        FTVAER    LW  SETVOLonly Aerosol Optical Depth
C                        -----------------------------------------------

C     -----------------------------------------------------------------
C     VEFF0   Selects Size Distribution Variance (this affects Thermal)
C     REFF0   Selects Effective Particle Size for Archive Volcanic Data
C     -----------------------------------------------------------------

      FSXTAU=FSTAER*FSVAER
      FTXTAU=FTTAER*FTVAER

C                   Set Grid-Box Edge Latitudes for Data Repartitioning
C                   ---------------------------------------------------
      NJ25=25
      DO 110 J=2,24
      E24LAT(J)=-90.D0+(J-1.5D0)*180.D0/23.D0
  110 CONTINUE
      E24LAT( 1)=-90.D0
      E24LAT(25)= 90.D0
      NJJM=46+1
      DO 120 J=2,46
      EJMLAT(J)=-90.D0+(J-1.5D0)*180.D0/(MLAT46-1)
  120 CONTINUE
      EJMLAT(   1)=-90.D0
      EJMLAT(NJJM)= 90.D0

!obso DO 140 M=1,13
!obso DO 130 J=1,27
!obso SIZA13(J,M)=REFF0
! 130 CONTINUE
! 140 CONTINUE

      DO 150 L=1,NL
      HTPROF(L)=0.D0
  150 CONTINUE

C                       -----------------------------------------------
C                       Initialize H2SO4 Q,S,C,A Tables for Input VEFF0
C                       -----------------------------------------------
C     ------------------
      CALL SETQVA(VEFF0)
C     ------------------

      RETURN


C--------------------------------
      ENTRY UPDVOL(JYEARV,JDAYVA)
C--------------------------------

C                                          (Makiko's 1850-1999 data)
C                                          -------------------------
      XYYEAR=JYEARV+JDAYVA/366.D0
      IF(XYYEAR.LT.1850.D0) XYYEAR=1850.D0
      XYI=(XYYEAR-1850.D0)*12.D0+1.D0
      IF(XYI.GT.1799.999D0) XYI=1799.999D0
      MI=XYI
      WMJ=XYI-MI
      WMI=1.D0-WMJ
      MJ=MI+1
      DO 250 J=1,24
      GDATA(J)=WMI*V4TAUR(MI,J,5)+WMJ*V4TAUR(MJ,J,5)
  250 CONTINUE
      CALL RETERP(GDATA,E24LAT,NJ25,SIZLAT,EJMLAT,NJJM)
      DO 270 K=1,4
      DO 260 J=1,24
      FDATA(J)=WMI*V4TAUR(MI,J,K)+WMJ*V4TAUR(MJ,J,K)
  260 CONTINUE
      CALL RETERP(FDATA,E24LAT,NJ25,HTFLAT(1,K),EJMLAT,NJJM)
  270 CONTINUE
      DO 290 J=1,24
      SUM=0
      DO 280 K=1,4
      SUM=SUM+HTFLAT(J,K)
  280 CONTINUE
      TAULAT(J)=SUM
  290 CONTINUE

      RETURN


C-----------------
      ENTRY GETVOL
C-----------------

      IF(MRELAY.GT.0)    GO TO 300
      IF(JLAT.EQ.LATVOL) GO TO 350

C                      Set JLAT Dependent Aerosol Distribution and Size
C                      ------------------------------------------------
  300 CONTINUE

      DO 310 K=1,4
      HLATTF(K)=HTFLAT(JLAT,K)
  310 CONTINUE
      CALL REPART(HLATTF,HLATKM,5,HTPROF,HLB0,NLP)
ccc   LHPMAX=0       ! not used
ccc   LHPMIN=NL      ! not used
      DO 320 L=1,NL
      N=NLP-L
ccc   IF(HTPROF(L).GE.HTPLIM) LHPMAX=L
ccc   IF(HTPROF(N).GE.HTPLIM) LHPMIN=N
  320 CONTINUE
      SUMHTF=1.D-10
      DO 330 L=1,NL
      IF(HTPROF(L).LT.HTPLIM) HTPROF(L)=0.D0
      SUMHTF=SUMHTF+HTPROF(L)
  330 CONTINUE

!obso TAUVOL=TAULAT(JLAT)
      SIZVOL=SIZLAT(JLAT)

C                        Select H2SO4 Q,S,C,A Tables for  Size = SIZVOL
C                        ----------------------------------------------

C------------------------
      CALL GETQVA(SIZVOL)
C------------------------

      LATVOL=JLAT
  350 CONTINUE
C                                  ------------------------------------
C                                  H2SO4 Thermal Contribution in TRVALK
C                                  ------------------------------------
      DO 420 K=1,33
      DO 410 L=1,NL
      TRVALK(L,K)=HTPROF(L)*AVH2S(K)*FTXTAU
  410 CONTINUE
  420 CONTINUE

C                      H2SO4 Solar Contribution in SRVEXT,SRVSCT,SRVGCB
C                      ------------------------------------------------

      DO 440 K=1,6
      DO 430 L=1,NL
      SRVEXT(L,K)=QVH2S(K)*HTPROF(L)*FSXTAU
      SRVSCT(L,K)=SVH2S(K)*HTPROF(L)*FSXTAU*PIVMAX
      SRVGCB(L,K)=GVH2S(K)
  430 CONTINUE
  440 CONTINUE

      RETURN
      END SUBROUTINE SETVOL

      SUBROUTINE SETQVA(VEFF)

      !INCLUDE 'BR00B.COM'

C     ------------------------------------------------------------------
C     SETQVA   Selects (interpolates) H2SO4 Mie Parameters for specified
C              Variance VEFF for subsequent Size interpolation by GETQVA
C     ------------------------------------------------------------------

ceq   DIMENSION SRQV( 6,20),SRSV( 6,20),SRGV( 6,20),Q55V(   20),REFV(20)
ceq   DIMENSION TRQV(33,20),TRSV(33,20),TRGV(33,20),TRAV(33,20),VEFV(20)
      DIMENSION TRAB(33,20),V5(5),Q5(5),RV20(20),QV20(20)

      DATA V5/0.1D0,0.2D0,0.3D0,0.4D0,0.5D0/
      SAVE TRAB,V5

C     ------------------------------------------------------------------
C     SRVQEX Volcanic Aerosol sizes (Reff) range from 0.1 to 5.0 microns
C     To utilize equal interval interpolation, Reff N=9,20 are redefined
C     so Volcanic Aerosol sizes have effective range of 0.1-2.0 microns.
C     ------------------------------------------------------------------

      DO 130 N=1,20
      RV20(N)=REFV20(N,1)
      VEFV(N)=VEFF
      REFV(N)=N/10.D0
      DO 110 I=1,5
      Q5(I)=Q55V20(N,I)
  110 CONTINUE
      CALL SPLINE(V5,Q5,5,VEFF,Q55V(N),1.D0,1.D0,1)
      DO 115 K=1,6
      DO 111 I=1,5
      Q5(I)=SRVQEX(K,N,I)
  111 CONTINUE
      CALL SPLINE(V5,Q5,5,VEFF,SRQV(K,N),1.D0,1.D0,1)
      DO 112 I=1,5
      Q5(I)=SRVQSC(K,N,I)
  112 CONTINUE
      CALL SPLINE(V5,Q5,5,VEFF,SRSV(K,N),1.D0,1.D0,1)
      DO 113 I=1,5
      Q5(I)=SRVQCB(K,N,I)
  113 CONTINUE
      CALL SPLINE(V5,Q5,5,VEFF,SRGV(K,N),1.D0,1.D0,1)
  115 CONTINUE
      DO 120 K=1,33
      DO 116 I=1,5
      Q5(I)=TRVQEX(K,N,I)
  116 CONTINUE
      CALL SPLINE(V5,Q5,5,VEFF,TRQV(K,N),1.D0,1.D0,1)
      DO 117 I=1,5
      Q5(I)=TRVQSC(K,N,I)
  117 CONTINUE
      CALL SPLINE(V5,Q5,5,VEFF,TRSV(K,N),1.D0,1.D0,1)
      DO 118 I=1,5
      Q5(I)=TRVQCB(K,N,I)
  118 CONTINUE
      CALL SPLINE(V5,Q5,5,VEFF,TRGV(K,N),1.D0,1.D0,1)
      DO 119 I=1,5
      Q5(I)=TRVQAL(K,N,I)
  119 CONTINUE
      CALL SPLINE(V5,Q5,5,VEFF,TRAV(K,N),1.D0,1.D0,1)
  120 CONTINUE
      DO 125 K=1,33
      TRAB(K,N)=TRQV(K,N)-TRSV(K,N)
  125 CONTINUE
  130 CONTINUE
      DO 131 N=1,20
      QV20(N)=Q55V(N)
  131 CONTINUE
      DO 132 N=9,20
      REFN=REFV(N)
      CALL SPLINE(RV20,QV20,20,REFN,Q55V(N),1.D0,1.D0,1)
  132 CONTINUE
      DO 140 K=1,6
      DO 133 N=1,20
      QV20(N)=SRQV(K,N)
  133 CONTINUE
      DO 134 N=9,20
      REFN=REFV(N)
      CALL SPLINE(RV20,QV20,20,REFN,SRQV(K,N),1.D0,1.D0,1)
  134 CONTINUE
      DO 135 N=1,20
      QV20(N)=SRSV(K,N)
  135 CONTINUE
      DO 136 N=9,20
      REFN=REFV(N)
      CALL SPLINE(RV20,QV20,20,REFN,SRSV(K,N),1.D0,1.D0,1)
  136 CONTINUE
      DO 137 N=1,20
      QV20(N)=SRGV(K,N)
  137 CONTINUE
      DO 138 N=9,20
      REFN=REFV(N)
      CALL SPLINE(RV20,QV20,20,REFN,SRGV(K,N),1.D0,1.D0,1)
  138 CONTINUE
  140 CONTINUE
      DO 150 K=1,33
      DO 141 N=1,20
      QV20(N)=TRQV(K,N)
  141 CONTINUE
      DO 142 N=9,20
      REFN=REFV(N)
      CALL SPLINE(RV20,QV20,20,REFN,TRQV(K,N),1.D0,1.D0,1)
  142 CONTINUE
      DO 143 N=1,20
      QV20(N)=TRSV(K,N)
  143 CONTINUE
      DO 144 N=9,20
      REFN=REFV(N)
      CALL SPLINE(RV20,QV20,20,REFN,TRSV(K,N),1.D0,1.D0,1)
  144 CONTINUE
      DO 145 N=1,20
      QV20(N)=TRGV(K,N)
  145 CONTINUE
      DO 146 N=9,20
      REFN=REFV(N)
      CALL SPLINE(RV20,QV20,20,REFN,TRGV(K,N),1.D0,1.D0,1)
  146 CONTINUE
      DO 147 N=1,20
      QV20(N)=TRAV(K,N)
  147 CONTINUE
      DO 148 N=9,20
      REFN=REFV(N)
      CALL SPLINE(RV20,QV20,20,REFN,TRAV(K,N),1.D0,1.D0,1)
  148 CONTINUE
      DO 149 N=9,20
      TRAB(K,N)=TRQV(K,N)-TRSV(K,N)
  149 CONTINUE
  150 CONTINUE

      RETURN

C-------------------------
      ENTRY GETQVA(SIZVOL)
C-------------------------

C     ------------------------------------------------------------------
C     Volcanic Aerosol sizes have effective range of  0.1 - 2.0 microns.
C     ------------------------------------------------------------------

      RADX=SIZVOL*10.D0
      IF(RADX.LT.1.000001D0) RADX=1.000001D0
      IF(RADX.GT.19.99999D0) RADX=19.99999D0
      JRXLO=RADX
      WTJHI=RADX-JRXLO
      WTJLO=1.D0-WTJHI
      JRXHI=JRXLO+1

      DO 210 I=1,6
      QVH2S(I)=WTJLO*SRQV(I,JRXLO)+WTJHI*SRQV(I,JRXHI)
      SVH2S(I)=WTJLO*SRSV(I,JRXLO)+WTJHI*SRSV(I,JRXHI)
      GVH2S(I)=WTJLO*SRGV(I,JRXLO)+WTJHI*SRGV(I,JRXHI)
  210 CONTINUE

      Q55H2S=WTJLO*Q55V(JRXLO)+WTJHI*Q55V(JRXHI)

      DO 220 I=1,33
      AVH2S(I)=WTJLO*TRAB(I,JRXLO)+WTJHI*TRAB(I,JRXHI)
  220 CONTINUE

      RETURN
      END SUBROUTINE SETQVA

      SUBROUTINE SETCLD

      !INCLUDE 'BR00B.COM'

C-----------------------------------------------------------------------
C     Control Parameters used in SETCLD,GETCLD,GETEPS: defined in RADCOM
C
C             ICE012  Selects Water, Non-Mie, Mie  Ice Cloud Qex,Qsc,Pi0
C             TAUWC0  Minimum Optical Depth for Water Clouds
C             TAUIC0  Minimum Optical Depth for   Ice Clouds
C             FCLDTR  Scaling Factor for Thermal Cloud Optical Depth
C             FCLDSR  Scaling Factor for  Solar  Cloud Optical Depth
C             EPSCON  Column Cloud Inhomogeneity EPSILON (when KCLDEP=1)
C             KCLDEP  Selects Cloud Inhomogeneity Option (0-4):
C                     KCLDEP =  0  Sets Column CLDEPS to Zero
C                     KCLDEP =  1  Sets Column CLDEPS to EPSCON
C                     KCLDEP =  2  Keeps whatever is specified in CLDEPS
C                     KCLDEP =  3  Uses: Column EPCOL(72,46) Climatology
C                     KCLDEP =  4  Uses: Ht Dep EPLOW, EPMID, EPHIG Data
C
C-----------------------------------------------------------------------
C                    Define Cloud Absorption Cross-Sections
C
C     Selected by:   ICE012 = 0    Liquid Water Droplets   (N =  1 -  5)
C                    ICE012 = 1    Ice - Non-Spherical     (N =  6 - 10)
C                    ICE012 = 2    Ice - Mie (Spherical)   (N = 11 - 15)
C
C     Define Solar,Thermal Cloud Single Scattering Albedo: SRCQPI( 6,15)
C                                                          TRCQPI(33,15)
C-----------------------------------------------------------------------

      DO 120 N=1,15
      DO 110 K=1,33
      TRCQAB(K,N)=TRCQEX(K,N)-TRCQSC(K,N)
      TRCQPI(K,N)=TRCQSC(K,N)/TRCQEX(K,N)
  110 CONTINUE
  120 CONTINUE

      DO 140 N=1,15
      DO 130 K=1,6
      SRCQPI(K,N)=SRCQSC(K,N)/SRCQEX(K,N)
  130 CONTINUE
  140 CONTINUE

C                          Initialize  GETCLD Output Parameters to Zero
C                          --------------------------------------------
CCC   LBOTCW=0
CCC   LTOPCW=0
CCC   LBOTCI=0
CCC   LTOPCI=0
      DO 150 K=1,33
      TRCTCA(K)=0.D0
  150 CONTINUE
      DO 180 L=1,NL
      DO 160 K=1,33
      TRCALK(L,K)=0.D0
  160 CONTINUE
      DO 170 K=1,6
      SRCEXT(L,K)=1.D-20
      SRCSCT(L,K)=0.D0
      SRCGCB(L,K)=0.D0
  170 CONTINUE
  180 CONTINUE

      RETURN
C-----------------
      ENTRY GETCLD
C-----------------

C-----------------------------------------------------------------------
C           Define:    TRCALK(LX,33)  Thermal Radiation Cloud Absorption
C                      TRCTCA(33)     Thermal Radiation Top Cloud Albedo
C
C                      SRCEXT(LX,6)   Solar Radiation Cloud Ext Op Depth
C                      SRCSCT(LX,6)   Solar Radiation Cloud Sct Op Depth
C                      SRCGCB(LX,6)   Solar Radiation Cloud Asym Param g
C
C                         LTOPCL      Top Cloud Layer Location
C                         LBOTCL      Bot Cloud Layer Location
C
C                         LTOPCW      Top Water Cloud Layer Location
C                         LBOTCW      Bot Water Cloud Layer Location
C
C                         LTOPCI      Top Ice Cloud Layer Location
C                         LBOTCI      Bot Ice Cloud Layer Location
C
C-----------------------------------------------------------------------

      LBOTCW=0
      LTOPCW=0
      LBOTCI=0
      LTOPCI=0
      DO 210 K=1,33
      TRCTCA(K)=0.D0
  210 CONTINUE
      DO 280 L=1,NL
      DO 220 K=1,33
      TRCALK(L,K)=0.D0
  220 CONTINUE
      DO 230 K=1,6
      SRCEXT(L,K)=1.D-20
      SRCSCT(L,K)=1.D-30
      SRCGCB(L,K)=0.D0
      SRCPI0(L,K)=0.D0
  230 CONTINUE
C                                         Water Cloud Size Interpolation
C                                         ------------------------------

      IF(TAUWC(L).GT.TAUWC0) THEN
      SIZWCL=SIZEWC(L)
      LTOPCW=L
      IF(LBOTCW.EQ.0) LBOTCW=L
      IF(SIZWCL.LT.15.D0) THEN
      IF(SIZWCL.LT.3.0D0) SIZWCL=3.0D0
      IRWAT=2
      XRW=SIZWCL/10.0D0-1.00D0
      ELSE
      IF(SIZWCL.GT.25.D0) SIZWCL=25.D0
      IRWAT=4
      XRW=SIZWCL/10.0D0-2.00D0
      ENDIF
      XMW=1.D0-XRW-XRW
      XPW=1.D0+XRW+XRW
      EPS=CLDEPS(L)
      VEP=EPS/(1.D0-EPS)
      VEP1=1.D0+VEP
      TAUWCL=TAUWC(L)
      DO 240 K=1,33
      QAWATK=XMW*XPW*TRCQAB(K,IRWAT)
     +      -XMW*XRW*TRCQAB(K,IRWAT-1)+XPW*XRW*TRCQAB(K,IRWAT+1)
      QPWATK=XMW*XPW*TRCQPI(K,IRWAT)
     +      -XMW*XRW*TRCQPI(K,IRWAT-1)+XPW*XRW*TRCQPI(K,IRWAT+1)
      VEPP=VEP*QPWATK
      TRCALK(L,K)=TRCALK(L,K)+TAUWCL*QAWATK/(VEP1-VEPP)
  240 CONTINUE
      SRCGFW=SRCGSF(1)
      DO 250 K=1,6
      QXWATK=XMW*XPW*SRCQEX(K,IRWAT)
     +      -XMW*XRW*SRCQEX(K,IRWAT-1)+XPW*XRW*SRCQEX(K,IRWAT+1)
      QSWATK=XMW*XPW*SRCQSC(K,IRWAT)
     +      -XMW*XRW*SRCQSC(K,IRWAT-1)+XPW*XRW*SRCQSC(K,IRWAT+1)
      QGWATK=XMW*XPW*SRCQCB(K,IRWAT)
     +      -XMW*XRW*SRCQCB(K,IRWAT-1)+XPW*XRW*SRCQCB(K,IRWAT+1)
      QPWATK=XMW*XPW*SRCQPI(K,IRWAT)
     +      -XMW*XRW*SRCQPI(K,IRWAT-1)+XPW*XRW*SRCQPI(K,IRWAT+1)
      QGWATK=QGWATK*SRCGFW
      VEPP=VEP*QPWATK
      VEP2=VEP1-VEPP
      SRCEXT(L,K)=SRCEXT(L,K)+TAUWCL*QXWATK/VEP1
      SRCSCT(L,K)=TAUWCL*QSWATK/(VEP1*VEP2)
      SRCGCB(L,K)=QGWATK*VEP2/(VEP1-VEPP*QGWATK)
  250 CONTINUE
      ENDIF
C                                           Ice Cloud Size Interpolation
C                                           ----------------------------
      IF(TAUIC(L).GT.TAUIC0) THEN
      SIZICL=SIZEIC(L)
      LTOPCI=L
      IF(LBOTCI.EQ.0) LBOTCI=L
      IF(SIZICL.LT.25.D0) THEN
      IF(SIZICL.LT.3.0D0) SIZICL=3.0D0
      IRICE=2+ICE012*5
      XRI=SIZICL/20.D0-0.75D0
      ELSE
      IF(SIZICL.GT.75.D0) SIZICL=75.D0
      IRICE=4+ICE012*5
      XRI=SIZICL/50.D0-1.00D0
      ENDIF
      XMI=1.D0-XRI-XRI
      XPI=1.D0+XRI+XRI
      EPS=CLDEPS(L)
      VEP=EPS/(1.D0-EPS)
      VEP1=1.D0+VEP
      TAUICL=TAUIC(L)
      DO 260 K=1,33
      QAICEK=XMI*XPI*TRCQAB(K,IRICE)
     +      -XMI*XRI*TRCQAB(K,IRICE-1)+XPI*XRI*TRCQAB(K,IRICE+1)
      QPICEK=XMI*XPI*TRCQPI(K,IRICE)
     +      -XMI*XRI*TRCQPI(K,IRICE-1)+XPI*XRI*TRCQPI(K,IRICE+1)
      VEPP=VEP*QPICEK
      TRCALK(L,K)=TRCALK(L,K)+TAUICL*QAICEK/(VEP1-VEPP)
  260 CONTINUE

      SRCGFC=SRCGSF(2)
      IF(ICE012.EQ.2) SRCGFC=SRCGSF(3)
      DO 270 K=1,6
      QXICEK=XMI*XPI*SRCQEX(K,IRICE)
     +      -XMI*XRI*SRCQEX(K,IRICE-1)+XPI*XRI*SRCQEX(K,IRICE+1)
      QSICEK=XMI*XPI*SRCQSC(K,IRICE)
     +      -XMI*XRI*SRCQSC(K,IRICE-1)+XPI*XRI*SRCQSC(K,IRICE+1)
      QGICEK=XMI*XPI*SRCQCB(K,IRICE)
     +      -XMI*XRI*SRCQCB(K,IRICE-1)+XPI*XRI*SRCQCB(K,IRICE+1)
      QPICEK=XMI*XPI*SRCQPI(K,IRICE)
     +      -XMI*XRI*SRCQPI(K,IRICE-1)+XPI*XRI*SRCQPI(K,IRICE+1)
      QGICEK=QGICEK*SRCGFC
      VEPP=VEP*QPICEK
      VEP2=VEP1-VEPP
      SRCEXT(L,K)=SRCEXT(L,K)+TAUICL*QXICEK/VEP1
      SCTTAU=TAUICL*QSICEK/(VEP1*VEP2)
      GCBICE=QGICEK*VEP2/(VEP1-VEPP*QGICEK)
      SCTGCB=SRCSCT(L,K)*SRCGCB(L,K)+SCTTAU*GCBICE
      SRCSCT(L,K)=SRCSCT(L,K)+SCTTAU
      SRCGCB(L,K)=SCTGCB/SRCSCT(L,K)
  270 CONTINUE
      ENDIF
  280 CONTINUE

C     ------------------------------------------------------------------
C     Identify Top Cloud (LTOPCL) and define top cloud albedo correction
C
C     Full Scattering Correction:      KCLDEM=1   ECLTRA=1.0   (default)
C     Partial(rad99a) Correction:      KCLDEM=0   ECLTRA=1.0
C       No Scattering Correction:      KCLDEM=0   ECLTRA=0.0
C
C     KCLDEM=1 Top-cloud scattering correction uses TXCTPG,TSCTPG,TGCTPG
C              to generate correction (over-rides old ECLTRA correction)
C              (KCLDEM correction is computed in THRMAL at LTOPCL level)
C     ------------------------------------------------------------------

      LTOPCL=LTOPCI
      IF(LTOPCI.GT.LTOPCW) GO TO 330
      IF(LTOPCW.LT.1) GO TO 350
      LTOPCL=LTOPCW
      TCTAUW=TAUWC(LTOPCL)
      DO 310 K=1,33
      ALWATK=XMW*XPW*TRCQAL(K,IRWAT)
     +      -XMW*XRW*TRCQAL(K,IRWAT-1)+XPW*XRW*TRCQAL(K,IRWAT+1)
      QXWATK=XMW*XPW*TRCQEX(K,IRWAT)
     +      -XMW*XRW*TRCQEX(K,IRWAT-1)+XPW*XRW*TRCQEX(K,IRWAT+1)
      TRCTCA(K)=(1.D0-EXP(-TAUWC(LTOPCL)*QXWATK))*ALWATK*ECLTRA
      QSWATK=XMW*XPW*TRCQSC(K,IRWAT)
     +      -XMW*XRW*TRCQSC(K,IRWAT-1)+XPW*XRW*TRCQSC(K,IRWAT+1)
      QGWATK=XMW*XPW*TRCQCB(K,IRWAT)
     +      -XMW*XRW*TRCQCB(K,IRWAT-1)+XPW*XRW*TRCQCB(K,IRWAT+1)
      TXCTPG(K)=QXWATK*TCTAUW
      TSCTPG(K)=QSWATK*TCTAUW
      TGCTPG(K)=QGWATK
  310 CONTINUE
      LBOTCL=LBOTCW
      IF(LBOTCI.LT.1) GO TO 360
      IF(LBOTCI.LE.LBOTCW) LBOTCL=LBOTCI
      IF(LTOPCI.EQ.LTOPCW) THEN
      TCTAUW=TAUWC(LTOPCL)
      TCTAUC=TAUIC(LTOPCL)
      WTI=TAUIC(LTOPCL)/(TAUIC(LTOPCL)+TAUWC(LTOPCL))
      WTW=TAUWC(LTOPCL)/(TAUIC(LTOPCL)+TAUWC(LTOPCL))
      DO 320 K=1,33
      ALICEK=XMI*XPI*TRCQAL(K,IRICE)
     +      -XMI*XRI*TRCQAL(K,IRICE-1)+XPI*XRI*TRCQAL(K,IRICE+1)
      QXICEK=XMI*XPI*TRCQEX(K,IRICE)
     +      -XMI*XRI*TRCQEX(K,IRICE-1)+XPI*XRI*TRCQEX(K,IRICE+1)
      TRCTCI=(1.D0-EXP(-TAUIC(LTOPCL)*QXICEK))*ALICEK*ECLTRA
      TRCTCA(K)=WTW*TRCTCA(K)+WTI*TRCTCI
      QSICEK=XMI*XPI*TRCQSC(K,IRICE)
     +      -XMI*XRI*TRCQSC(K,IRICE-1)+XPI*XRI*TRCQSC(K,IRICE+1)
      QGICEK=XMI*XPI*TRCQCB(K,IRICE)
     +      -XMI*XRI*TRCQCB(K,IRICE-1)+XPI*XRI*TRCQCB(K,IRICE+1)
      TXCTPG(K)=TXCTPG(K)+QXICEK*TCTAUC
      SCTGCB=TSCTPG(K)*TGCTPG(K)+QSICEK*TCTAUC*QGICEK
      TSCTPG(K)=TSCTPG(K)+QSICEK*TCTAUC
      TGCTPG(K)=SCTGCB/(1.D-10+TSCTPG(K))
  320 CONTINUE
      ENDIF
      GO TO 360
  330 CONTINUE
      LTOPCL=LTOPCI
      TCTAUC=TAUIC(LTOPCL)
      DO 340 K=1,33
      ALICEK=XMI*XPI*TRCQAL(K,IRICE)
     +      -XMI*XRI*TRCQAL(K,IRICE-1)+XPI*XRI*TRCQAL(K,IRICE+1)
      QXICEK=XMI*XPI*TRCQEX(K,IRICE)
     +      -XMI*XRI*TRCQEX(K,IRICE-1)+XPI*XRI*TRCQEX(K,IRICE+1)
      TRCTCA(K)=(1.D0-EXP(-TAUIC(LTOPCL)*QXICEK))*ALICEK*ECLTRA
      QSICEK=XMI*XPI*TRCQSC(K,IRICE)
     +      -XMI*XRI*TRCQSC(K,IRICE-1)+XPI*XRI*TRCQSC(K,IRICE+1)
      QGICEK=XMI*XPI*TRCQCB(K,IRICE)
     +      -XMI*XRI*TRCQCB(K,IRICE-1)+XPI*XRI*TRCQCB(K,IRICE+1)
      TXCTPG(K)=QXICEK*TCTAUC
      TSCTPG(K)=QSICEK*TCTAUC
      TGCTPG(K)=QGICEK
  340 CONTINUE
      LBOTCL=LBOTCI
      IF(LBOTCW.EQ.0) GO TO 360
      LBOTCL=LBOTCW
      GO TO 360
  350 CONTINUE
      LBOTCL=0
      LTOPCL=0
  360 CONTINUE

      RETURN


C--------------------------------
      ENTRY UPDEPS(JYEARE,JJDAYE)
C--------------------------------

C                 Select ISCCP-Based Cloud Heterogeneity Time Dependence
C                 ------------------------------------------------------

      XJDAY=JJDAYE-0.999D0
      XMO=XJDAY/30.5D0+.5D0
      MI=XMO
      WTMJ=XMO-MI
      WTMI=1.D0-WTMJ
      IF(MI.LT.1) MI=12
      MJ=MI+1
      IF(MJ.GT.12) MJ=1

      DO 420 J=1,46
      DO 410 I=1,72
      EPLOW(I,J)=WTMI*EPLMHC(I,J,MI,1)+WTMJ*EPLMHC(I,J,MJ,1)
      EPMID(I,J)=WTMI*EPLMHC(I,J,MI,2)+WTMJ*EPLMHC(I,J,MJ,2)
      EPHIG(I,J)=WTMI*EPLMHC(I,J,MI,3)+WTMJ*EPLMHC(I,J,MJ,3)
      EPCOL(I,J)=WTMI*EPLMHC(I,J,MI,4)+WTMJ*EPLMHC(I,J,MJ,4)
  410 CONTINUE
  420 CONTINUE

      RETURN


C-----------------
      ENTRY GETEPS
C-----------------

C             ----------------------------------------------------------
C                     Select Cloud Heterogeneity CLDEPS Options
C             EPSCON  Column Cloud Inhomogeneity EPSILON (when KCLDEP=1)
C             KCLDEP  Selects Cloud Inhomogeneity Option (0-4):
C                     KCLDEP =  0  Sets Column CLDEPS to Zero
C                     KCLDEP =  1  Sets Column CLDEPS to EPSCON
C                     KCLDEP =  2  Keeps whatever is specified in CLDEPS
C                     KCLDEP =  3  Uses: Column EPCOL(72,46) Climatology
C                     KCLDEP =  4  Uses: Ht Dep EPLOW, EPMID, EPHIG Data
C                     --------------------------------------------------

      IF(KCLDEP.EQ.0) THEN
      DO 510 L=1,NL
      CLDEPS(L)=0.D0
  510 CONTINUE
      ENDIF
      IF(KCLDEP.EQ.1) THEN
      DO 520 L=1,NL
      CLDEPS(L)=EPSCON
  520 CONTINUE
      ENDIF
      I=ILON
      J=JLAT
      IF(KCLDEP.EQ.3) THEN
      DO 530 L=1,NL
      CLDEPS(L)=EPCOL(I,J)
  530 CONTINUE
      ENDIF
      IF(KCLDEP.EQ.4) THEN
CCC   LOWCLD=0
CCC   MIDCLD=0
      DO 540 L=1,NL
      CLDEPS(L)=EPMID(I,J)
      IF(PLB(L).GT.750.D0) CLDEPS(L)=EPLOW(I,J)
      IF(PLB(L).LT.430.D0) CLDEPS(L)=EPHIG(I,J)
  540 CONTINUE
      ENDIF

      RETURN
      END SUBROUTINE SETCLD


      SUBROUTINE TAUGAS

      !INCLUDE 'BR00B.COM'

C     ----------------------------------------------------------
C     TAUGAS INPUT REQUIRES:  NL,PL,DPL,TLM,ULGAS
C                             TAUTBL,TAUWV0,XKCFC,H2OCN8,H2OCF8
C                             DUCH4,SDUCH4,DUN2O,SDUN2O,ULOX,DUX
C     TAUGAS OUTPUT DATA IS:  TRGXLK
C     ----------------------------------------------------------

      PARAMETER(NTX=8, TLOX=181.d0, DTX=23.d0, NPX=19, NGUX=968,
     *       NPUX=19, NPU0=14, NPU2=14, DELCH4=0.10D0, DELN2O=0.05D0,
     *       NPU=5, P0=1013.25d0)

      DIMENSION IGASX(18),KGX(18),NUX(13),IGUX(13),NGX(4),IG1X(4)
     +     ,PX(19),XKCFCW(8,2),MLGAS(18),PDPU2(14),PU(5)
      DATA       PX/1000.,750.,500.,300.,200.,100.,50.,20.,10.,5.,
     *                 2.,1.,.5,.2,.1,.03,.01,.003,.001/
C
      DATA NGX/12,12,08,33/, IG1X/2,14,26,1/
      DATA PDPU2/1.E4,1.E5,2.E5,5.E5,1.E6,2.E6,5.E6,1.E7,2.E7,
     * 5.E7,1.E8,2.E8,5.E8,1.E9/
      DATA    PU/  50.,200.,800.,3200.,12800./

      DATA IGASX/  1,  2,  3,  1,  1,  2,  2,  3,  3,  6,  6,  6,  7,
     *             8,  8,  9,  9,  1/
      DATA   KGX/  1,  2,  3,  2,  3,  1,  3,  1,  2,  1,  2,  3,  1,
     *             2,  3,  2,  3,  4/
      DATA   NUX/ 25,  9,  9,  9,  9,  5,  5,  5,  5,  2,  2,  2,  2/
      DATA  IGUX/  0,300,408,480,588,660,720,760,820,880,904,928,944/
C
      DATA XKCFCW/
     $ 12.4414,11.7842,11.3630,10.8109,10.3200, 9.8900, 9.3916, 8.8933,
     $  5.3994, 5.6429, 5.8793, 6.1687, 6.2300, 6.5200, 6.8650, 7.2100/

      DIMENSION RCO2(5),FCO2(5),FH2O(5),RH2O(8),XFUL(8)
      DATA RCO2/0.25,0.50,1.00,2.00,4.00/
      DATA FCO2/1.07,1.03,1.00,0.96,0.92/
      DATA FH2O/1.02,1.01,1.00,1.00,1.00/
      DATA RH2O/0.0625,0.125,0.25,0.50,1.00,2.00,4.00,8.00/
      DATA XFUL/0.4600,0.540,0.60,0.64,0.66,0.69,0.72,0.72/


      SAVE            IGASX,KGX,NUX,IGUX,NGX,IG1X,PX,XKCFCW,MLGAS,
     *                PDPU2,PU,RCO2,FCO2,FH2O,RH2O,XFUL

C----------------------------------------------------------------------
C     N2O and CH4 TOA flux change normalization to off-line LBL results
C          January 20, 2000 parameterization: DUCH4,SDUCH4,DUN2O,SDUN2O
C----------------------------------------------------------------------

      DO 100 K=1,33
      DO 100 L=1,NL
      TRGXLK(L,K)=0.D0
  100 CONTINUE

      DO 110 I=1,18
      MLGAS(I)=1
  110 CONTINUE

C              KWVCON = ON/OFF flag for water vapor continuum absorption
C              ---------------------------------------------------------
      IF(KWVCON.LT.1) MLGAS(18)=0

      UH2O=0.D0
      UCO2=0.D0
      UO3 =0.D0
      UN2O=0.D0
      UCH4=0.D0
      DO 120 IP=1,NL
      UH2O=UH2O+ULGAS(IP,1)
      UCO2=UCO2+ULGAS(IP,2)
      UO3 =UO3 +ULGAS(IP,3)
      UN2O=UN2O+ULGAS(IP,6)
      UCH4=UCH4+ULGAS(IP,7)
  120 CONTINUE
      UIN2O=UN2O/DELN2O+1.D0
      UICH4=UCH4/DELCH4+1.D0
      IN2O=UIN2O
      ICH4=UICH4
      IF(IN2O.GT.149) IN2O=149
      IF(ICH4.GT.149) ICH4=149
      DIN2O=UIN2O-IN2O
      DICH4=UICH4-ICH4
      TUN2OX=UN2O+DUN2O(IN2O)*(1.D0-DIN2O)+DUN2O(IN2O+1)*DIN2O
      TUCH4X=UCH4+DUCH4(ICH4)*(1.D0-DICH4)+DUCH4(ICH4+1)*DICH4
      TXUN2O=TUN2OX/UN2O
      TXUCH4=TUCH4X/UCH4

C--------------------------------------------------------------------
C        Absorption (TAU) interpolation for gas amounts in ULGAS(N,K)
C--------------------------------------------------------------------


      FULMLT=1.D0
      XU1=1.D0
      XU2=1.D0
      IF(KUFH2O.EQ.0) GO TO 170
      SUM1=0.D0
      SUM2=0.D0
      DO 130 L=1,NL
      SUM1=SUM1+ULGAS(L,1)
      SUM2=SUM2+ULGAS(L,2)
  130 CONTINUE
      URAT1=SUM1/3252.D0
      FULMLT=0.46D0
      IF(URAT1.LT.RH2O(1)) GO TO 150
      I=1
  140 CONTINUE
      IF(URAT1.LE.RH2O(I+1)) THEN
      A=RH2O(I+1)-URAT1
      B=URAT1-RH2O(I)
      DENO=RH2O(I+1)-RH2O(I)
      FULMLT=(XFUL(I)*A+XFUL(I+1)*B)/DENO
      ELSE
      I=I+1
      IF(I.LT.8) GO TO 140
      FULMLT=0.72D0
      END IF
  150 CONTINUE
      IF(KUFCO2.EQ.0) GO TO 170
      URAT2=SUM2/262.52D0
      XU1=1.02D0
      XU2=1.07D0
      IF(URAT2.LT.RCO2(1)) GO TO 170
      I=1
  160 CONTINUE
      IF(URAT2.LE.RCO2(I+1)) THEN
      A=RCO2(I+1)-URAT2
      B=URAT2-RCO2(I)
      DENO=RCO2(I+1)-RCO2(I)
      XU1=(FH2O(I)*A+FH2O(I+1)*B)/DENO
      XU2=(FCO2(I)*A+FCO2(I+1)*B)/DENO
      ELSE
      I=I+1
      IF(I.LT.5) GO TO 160
      XU1=1.00D0
      XU2=0.92D0
      ENDIF
  170 CONTINUE
      H2ORAT=MIN( UH2O*FULMLT*XU1/3671.D0 , 1.D0)
      CO2RAT=MIN( UCO2*XU2/252.D0         , 1.D0)
      O3RAT =MIN( UO3/.343D0              , 1.D0)

      IPX=2
      DO 500 IP=1,NL
  200 CONTINUE
      WPB = (PL(IP)-PX(IPX))/(PX(IPX-1)-PX(IPX))
      IF(WPB.GE.0.D0.OR.IPX.GE.NPX) GO TO 210
      IPX = IPX+1
      GO TO 200
  210 CONTINUE
      WTB = (TLM(IP)-TLOX)/DTX
      ITX = MIN0(MAX0(INT(WTB),0),NTX-2)
      WTB = WTB-FLOAT(ITX)

      WBB = WPB*WTB
      WBA = WPB-WBB
      WAB = WTB-WBB
      WAA = 1.D0-(WBB+WBA+WAB)

      IAA = NGUX*(ITX+NTX*(IPX-1))
      IBA = IAA-NGUX*NTX

      ITX1=ITX+1
      ITX2=ITX+2
      IPXM1=IPX-1

      IULOW=0
      DO 400 IGAS=1,18
      IF(MLGAS(IGAS).LT.1) GO TO 400
      NG = NGX(KGX(IGAS))
      UGAS = ULGAS(IP,IGASX(IGAS))
      KK=IG1X(KGX(IGAS))
C                                   Apply absorber scaling for H2O & CO2
C                                   ------------------------------------
      IF(IGAS.LT.18) THEN
      IF(IGASX(IGAS).EQ.1) UGAS=UGAS*FULMLT*XU1
      IF(IGASX(IGAS).EQ.2) UGAS=UGAS*FPXCO2(IP)*XU2

C                                         Modified scaling for CH4 & N2O
C                                         ------------------------------
      IF(PL(IP).GT.200.D0) THEN
      FPL=DMIN1((PL(IP)-200.D0)/300.D0,1.D0)
      IF(IGAS.EQ.10.OR.IGAS.EQ.13) FPL=FPL*H2ORAT
      IF(IGAS.EQ.11) FPL=FPL*CO2RAT*FPXCO2(IP)
      IF(IGAS.EQ.12) FPL=FPL*O3RAT
      IF(IGAS.EQ.13) UGAS=UGAS*(1.D0+FPL*(1.8D0*TXUCH4-1.D0))
      IF(IGAS.GT.9.AND.IGAS.LT.13)
     +               UGAS=UGAS*(1.D0+FPL*(1.3D0*TXUN2O-1.D0))
      ENDIF
      ENDIF
C                                 Apply water vapor continuum absorption
C                                 --------------------------------------
      IF(IGAS.EQ.18) THEN
      UGASW=UGAS*1.15D0
      UGS0=UGAS*FULMLT*XU1
      UGAS=UGS0*1.036D0

C                  KWSELF = ON/FF flag for H2O self broadening continuum
C                  -----------------------------------------------------
      IF(KCSELF.GT.0) THEN
      U=UGASW
      IK0=1
      IKF=1
      DO 245 I=1,3
      IF(I.EQ.2) THEN
      U=UGS0
      IK0=2
      IKF=13
      ENDIF
      IF(I.EQ.3) THEN
      U=UGAS
      IK0=14
      IKF=33
      ENDIF
      PU2=PL(IP)/DPL(IP)*U**2
      IF(PU2.GT.PDPU2(1)) THEN
      DO 220 IPU=2,NPU2
      IPUI=IPU
      IF(PU2.LE.PDPU2(IPU)) GO TO 230
  220 CONTINUE
  230 CONTINUE
      IPU=IPUI
      IPU1=IPU-1
      WTPU=(PU2-PDPU2(IPU1))/(PDPU2(IPU)-PDPU2(IPU1))
      ELSE
      WTPU=PU2/PDPU2(1)
      ENDIF

      DO 240 IK=IK0,IKF
      IF(PU2.GT.PDPU2(1)) THEN
      TAUT1=WTPU*(H2OCN8(IK,ITX1,IPU)-H2OCN8(IK,ITX1,IPU1))+
     $      H2OCN8(IK,ITX1,IPU1)
      TAUT2=WTPU*(H2OCN8(IK,ITX2,IPU)-H2OCN8(IK,ITX2,IPU1))+
     $      H2OCN8(IK,ITX2,IPU1)
      ELSE
      TAUT1=WTPU*H2OCN8(IK,ITX1,1)
      TAUT2=WTPU*H2OCN8(IK,ITX2,1)
      ENDIF

      TAUHCN=WTB*(TAUT2-TAUT1)+TAUT1
      TRGXLK(IP,KK)=TRGXLK(IP,KK)+TAUHCN
      KK=KK+1
  240 CONTINUE
  245 CONTINUE
      ENDIF
C               KWFORN = ON/FF flag for H2O foreign broadening continuum
C               --------------------------------------------------------

      IF(KCFORN.LT.1) GO TO 400
      KK=IG1X(KGX(IGAS))
      U=UGASW
      IK0=1
      IKF=1
      DO 275 I=1,3
      IF(I.EQ.2) THEN
      U=UGS0
      IK0=2
      IKF=13
      ENDIF
      IF(I.EQ.3) THEN
      U=UGAS
      IK0=14
      IKF=33
      ENDIF
      UP=PL(IP)/P0*U
      IF(UP.GT.PU(1)) THEN
      DO 250 IPU=2,NPU
      IPUI=IPU
      IF(UP.LE.PU(IPU)) GO TO 260
  250 CONTINUE
  260 CONTINUE
      IPU=IPUI
      IPU1=IPU-1
      WTPU=(UP-PU(IPU1))/(PU(IPU)-PU(IPU1))
      ELSE
      WTPU=UP/PU(1)
      ENDIF
      DO 270 IK=IK0,IKF
      IF(UP.GT.PU(1)) THEN
      TAUT1=WTPU*(H2OCF8(IK,ITX1,IPU)-H2OCF8(IK,ITX1,IPU1))+
     +      H2OCF8(IK,ITX1,IPU1)
      TAUT2=WTPU*(H2OCF8(IK,ITX2,IPU)-H2OCF8(IK,ITX2,IPU1))+
     +      H2OCF8(IK,ITX2,IPU1)
      ELSE
      TAUT1=WTPU*H2OCF8(IK,ITX1,1)
      TAUT2=WTPU*H2OCF8(IK,ITX2,1)
      ENDIF
      TAUHFB=WTB*(TAUT2-TAUT1)+TAUT1
      TRGXLK(IP,KK)=TRGXLK(IP,KK)+TAUHFB
      KK=KK+1
  270 CONTINUE
  275 CONTINUE
      GO TO 400
      ENDIF

      IF(IGAS.GT.13) THEN
      IGCFC=IGAS-13
      DO 280 IK=1,NG
      XA=WTB*(XKCFC(IK,ITX2,IPXM1,IGCFC)-XKCFC(IK,ITX1,IPXM1,IGCFC))+
     +        XKCFC(IK,ITX1,IPXM1,IGCFC)
      XB=WTB*(XKCFC(IK,ITX2,IPX,IGCFC)-XKCFC(IK,ITX1,IPX,IGCFC))+
     +        XKCFC(IK,ITX1,IPX,IGCFC)
      XK=WPB*(XA-XB)+XB
      TAUCF=XK*UGAS
      TRGXLK(IP,KK)=TRGXLK(IP,KK)+TAUCF
      KK=KK+1
  280 CONTINUE
      GO TO 400
      ENDIF
      IU = IPX + NPUX*(IGAS-1)
      NU = NUX(IGAS)
      IF(NU.GT.1) GO TO 290
      XUA = 0.D0
      XUB = 0.D0
      GO TO 300
  290 CONTINUE
      XUA = (UGAS-ULOX(IU))/DUX(IU)
      XUB = (UGAS-ULOX(IU-1))/DUX(IU-1)
  300 CONTINUE
      IUA = INT(XUA)
      IUB = INT(XUB)

      IF(IGAS.EQ.1) THEN
      IF(XUA.LT.0..AND.XUB.LT.0.) IULOW=1
      ENDIF

      QAA = 1.D0
      QAB = 1.D0
      IF(XUA.GT.0.D0.AND.IUA.LT.NU-1) GO TO 310
      FNU1=NU-1
      XUA = DMIN1(DMAX1(XUA,0.D0),FNU1)
      IUA = MIN0(INT(XUA),NU-2)
      QAA = UGAS/(ULOX(IU)+DUX(IU)*FLOAT(IUA))
      QAB = UGAS/(ULOX(IU)+DUX(IU)*FLOAT(IUA+1))
  310 CONTINUE
      QBA = 1.D0
      QBB = 1.D0
      IF(XUB.GT.0.D0.AND.IUB.LT.NU-1) GO TO 320
      FNU1=NU-1
      XUB = DMIN1(DMAX1(XUB,0.D0),FNU1)
      IUB = MIN0(INT(XUB),NU-2)
      QBA = UGAS/(ULOX(IU-1)+DUX(IU-1)*FLOAT(IUB))
      QBB = UGAS/(ULOX(IU-1)+DUX(IU-1)*FLOAT(IUB+1))
  320 CONTINUE
      UAB = XUA-FLOAT(IUA)
      UBB = XUB-FLOAT(IUB)
      UAA = 1.D0-UAB
      UBA = 1.D0-UBB

      WAAA = WAA*UAA*QAA
      WAAB = WAA*UAB*QAB
      WABA = WAB*UAA*QAA
      WABB = WAB*UAB*QAB
      WBAA = WBA*UBA*QBA
      WBAB = WBA*UBB*QBB
      WBBA = WBB*UBA*QBA
      WBBB = WBB*UBB*QBB

      IAAA = IAA+IGUX(IGAS) + NG*IUA
      IAAB = IAAA+NG
      IABA = IAAA+NGUX
      IABB = IABA+NG
      IBAA = IBA+IGUX(IGAS) + NG*IUB
      IBAB = IBAA+NG
      IBBA = IBAA+NGUX
      IBBB = IBBA+NG

      IH2O0=0
      IF(IGAS.EQ.6.OR.IGAS.EQ.8.OR.IGAS.EQ.10.OR.IGAS.EQ.13) THEN
      IF(IULOW.EQ.1) IH2O0=1
      ENDIF

      DO 330 IG=1,NG
      IF(IH2O0.EQ.0) THEN
      TAUIPG=
     +       WAAA*TAUTBL(IAAA+IG)
     +     + WAAB*TAUTBL(IAAB+IG)
     +     + WABA*TAUTBL(IABA+IG)
     +     + WABB*TAUTBL(IABB+IG)
     +     + WBAA*TAUTBL(IBAA+IG)
     +     + WBAB*TAUTBL(IBAB+IG)
     +     + WBBA*TAUTBL(IBBA+IG)
     +     + WBBB*TAUTBL(IBBB+IG)
      ELSE
      TAUIPG=
     +       WAAA*TAUWV0(IAAA+IG)
     +     + WAAB*TAUWV0(IAAB+IG)
     +     + WABA*TAUWV0(IABA+IG)
     +     + WABB*TAUWV0(IABB+IG)
     +     + WBAA*TAUWV0(IBAA+IG)
     +     + WBAB*TAUWV0(IBAB+IG)
     +     + WBBA*TAUWV0(IBBA+IG)
     +     + WBBB*TAUWV0(IBBB+IG)
      ENDIF

      TAUSUM=TRGXLK(IP,KK)+TAUIPG
      IF(TAUSUM.GT.0.D0) TRGXLK(IP,KK)=TAUSUM
      KK=KK+1
  330 CONTINUE
  400 CONTINUE
C                               CFC11 and CFC12 Window Absorption (1997)
C                               ----------------------------------------

      IF(MLGAS(14).EQ.1.OR.MLGAS(15).EQ.1) THEN
      XK=WTB*(XKCFCW(ITX2,1)-XKCFCW(ITX1,1))+XKCFCW(ITX1,1)
      TAU11=XK*(ULGAS(IP,8)+ULGAS(IP,11))
      TRGXLK(IP,1)=TRGXLK(IP,1)+TAU11
      ENDIF
      IF(MLGAS(16).EQ.1.OR.MLGAS(17).EQ.1) THEN
      XK=WTB*(XKCFCW(ITX2,2)-XKCFCW(ITX1,2))+XKCFCW(ITX1,2)
      TAU12=XK*(ULGAS(IP,9)+ULGAS(IP,12))
      TRGXLK(IP,1)=TRGXLK(IP,1)+TAU12
      ENDIF
  500 CONTINUE

      RETURN
      END SUBROUTINE TAUGAS


      SUBROUTINE THERML

      !INCLUDE 'BR00B.COM'

      PARAMETER (R6=.16666667D0,R24=4.1666667D-02)
      PARAMETER (A=0.3825D0,B=0.5742D0,C=0.0433D0)

C     ------------------------------------------------------------------
C             Top-cloud Thermal Scattering Correction Control Parameters
C             ----------------------------------------------------------
C
C             ECLTRA = 1.0  Scattering correction is enabled
C             with KCLDEM = 1, Rigorous scattering correction is applied
C             with KCLDEM = 0, Approximate scattering correction is used
C
C             ECLTRA = 0.0  No scattering correction is used
C                                          (Independent of KCLDEM value)
C
C     ------------------------------------------------------------------
C                                   Lower Edge Temperature Interpolation
C                                   ------------------------------------
C     TLGRAD=1.0  (Default)
C                 Layer-mean temperatures (TLM) supplied by GCM are used
C                 to define the layer edge temperature TLT (top) and TLB
C                 (bottom) using overall atmospheric temperature profile
C                 to establish temperature gradient within each layer so
C                 as to minimize the temperature discontinuities between
C                 layer edges and to conserve layer thermal energy.
C
C     TLGRAD=0.0  This results in isothermal layers with TLT = TLB = TLM
C
C     TLGRAD<0.0  TLT and TLB are used as specified, without any further
C                 adjustments.  This is mainly for off-line use when the
C                 temperature profile (TLM,TLT,TLB) can be fully defined
C                 from a continuous temperature profile.
C
C     NOTE:       TLGRAD can also accommodate values between 0.0 and 1.0
C
C     PTLISO      (Default PTLISO=2.5mb)
C                 Pressure level above which model layers are defined to
C                 be isothermal.  This is appropriate for optically thin
C                 layers where emitted flux depends on mean temperature.
C     ------------------------------------------------------------------


      IF(TLGRAD.LT.0.D0) GO TO 130
      TA=TLM(1)
      TB=TLM(2)
      P1=PLB(1)
      P2=PLB(2)
      P3=PLB(3)
      DT1CPT=0.5D0*TA*(PLB(1)**0.286D0-PLB(2)**0.286D0)/PL(1)**0.286D0
      DTHALF=(TA-TB)*(P1-P2)/(P1-P3)
      IF(DTHALF.GT.DT1CPT) DTHALF=DT1CPT
      TLB(1)=TA+DTHALF*TLGRAD
      TLT(1)=TA-DTHALF*TLGRAD
      DO 110 L=3,NL
      TC=TLM(L)
      P4=PLB(L+1)
      DTHALF=0.5D0*((TA-TB)/(P1-P3)+(TB-TC)/(P2-P4))*(P2-P3)*TLGRAD
      TLB(L-1)=TB+DTHALF
      TLT(L-1)=TB-DTHALF
      TA=TB
      TB=TC
      P1=P2
      P2=P3
      P3=P4
  110 CONTINUE
      DTHALF=(TA-TB)*(P2-P3)/(P1-P3)*TLGRAD
      TLB(NL)=TC+DTHALF
      TLT(NL)=TC-DTHALF
      DO 120 L=NL,1,-1
      IF(PLB(L).GT.PTLISO) GO TO 130
      TLT(L)=TLM(L)
      TLB(L)=TLM(L)
  120 CONTINUE
  130 CONTINUE

C     ------------------------------------------------------------------
C                   WEIGHT ASSIGNMENTS FOR PLANCK FUNCTION INTERPOLATION
C                    (Effective range is from TK = 123 K to TK = 373 K)
C     ------------------------------------------------------------------

      DO 140 L=1,NL
      ITL=TLB(L)
      WTLB(L)=TLB(L)-ITL
      ITLB(L)=ITL-ITPFT0
      IF(ITLB(L).LT.1) ITLB(L)=1
      IF(ITLB(L).GT.249) ITLB(L)=249
      ITL=TLT(L)
      WTLT(L)=TLT(L)-ITL
      ITLT(L)=ITL-ITPFT0
      IF(ITLT(L).LT.1) ITLT(L)=1
      IF(ITLT(L).GT.249) ITLT(L)=249
  140 CONTINUE

      DO 180 L=1,NL
      IF(L.EQ.LTOPCL) THEN
      DO 170 K=1,33
      CLTAUX=TXCTPG(K)+TRGXLK(L,K)+1.D-10
      CLTAUS=TSCTPG(K)
      CLCOSB=TGCTPG(K)
      CLPI0K=CLTAUS*ECLTRA/CLTAUX
      CLPI0(K)=CLPI0K
      CTX=CLTAUX*10.D0
      IF(CLTAUX.GE.3.D0) THEN
      CTX=CLTAUX*2.D0+24.D0
      IF(CTX.GT.47.999999D0) CTX=47.999999D0
      ENDIF
      ICT=CTX
      DT2=CTX-ICT
      DT1=1.D0-DT2
      IT1=ICT+1
      IT2=ICT+2
      CTP=CLPI0K*20.D0
      ICP=CTP
      DP2=CTP-ICP
      DP1=1.D0-DP2
      IP1=ICP+1
      IP2=ICP+2
      CTG=CLCOSB*20.D0
      ICG=CTG
      DG2=CTG-ICG
      DG1=1.D0-DG2
      IG1=ICG+1
      IG2=ICG+2
      WT1=DT1*DP1*DG1
      WT2=DT2*DP1*DG1
      WT3=DT2*DP2*DG1
      WT4=DT1*DP2*DG1
      WT5=DT1*DP1*DG2
      WT6=DT2*DP1*DG2
      WT7=DT2*DP2*DG2
      WT8=DT1*DP2*DG2
      DO 150 II=1,6
      RIJTCK(II,K)=WT1*RIJTPG(II,IT1,IP1,IG1)+WT2*RIJTPG(II,IT2,IP1,IG1)
     +            +WT3*RIJTPG(II,IT2,IP2,IG1)+WT4*RIJTPG(II,IT1,IP2,IG1)
     +            +WT5*RIJTPG(II,IT1,IP1,IG2)+WT6*RIJTPG(II,IT2,IP1,IG2)
     +            +WT7*RIJTPG(II,IT2,IP2,IG2)+WT8*RIJTPG(II,IT1,IP2,IG2)
  150 CONTINUE
      DO 160 II=1,3
      FEMTCK(II,K)=WT1*FEMTPG(II,IT1,IP1,IG1)+WT2*FEMTPG(II,IT2,IP1,IG1)
     +            +WT3*FEMTPG(II,IT2,IP2,IG1)+WT4*FEMTPG(II,IT1,IP2,IG1)
     +            +WT5*FEMTPG(II,IT1,IP1,IG2)+WT6*FEMTPG(II,IT2,IP1,IG2)
     +            +WT7*FEMTPG(II,IT2,IP2,IG2)+WT8*FEMTPG(II,IT1,IP2,IG2)
      FDXTCK(II,K)=WT1*FDXTPG(II,IT1,IP1,IG1)+WT2*FDXTPG(II,IT2,IP1,IG1)
     +            +WT3*FDXTPG(II,IT2,IP2,IG1)+WT4*FDXTPG(II,IT1,IP2,IG1)
     +            +WT5*FDXTPG(II,IT1,IP1,IG2)+WT6*FDXTPG(II,IT2,IP1,IG2)
     +            +WT7*FDXTPG(II,IT2,IP2,IG2)+WT8*FDXTPG(II,IT1,IP2,IG2)
  160 CONTINUE
  170 CONTINUE
      ENDIF
  180 CONTINUE
      DO 190 L=1,NLP
      TRDFLB(L)=0.D0
      TRUFLB(L)=0.D0
  190 CONTINUE

      BG=BGFEMT(1)
      TOTLZF(1)=0.D0
      TOTLZF(2)=0.D0
      TOTLZF(3)=0.D0
!!    TRSLTS=0.D0
!!    TRSLTG=0.D0
!!    TRSLWV=0.D0
!!    TRSLBS=0.D0

      K=0
      IMOL=0
  200 CONTINUE
      ITK0=K*ITNEXT
      K=K+1
      IF(K.GT.33) GO TO 300
      BG=BGFEMT(K)
      DFLB(NLP,K)=0.D0
      DNACUM=0.D0
      DNBCUM=0.D0
      DNCCUM=0.D0
      L=NL
  205 CONTINUE
      TAUAG=TRGXLK(L,K)
      TAUAP=TRCALK(L,K)+TRAALK(L,K)+TRBALK(L,K)+TRDALK(L,K)+TRVALK(L,K)
      TAUAX=TAUAG+TAUAP
      IF(TAUAX.GT.1.D-06) GO TO 215
      DFLB(L,K)=0.D0
      ENA(L)=0.D0
      DNA(L)=0.D0
      TRA(L)=1.D0
      ENB(L)=0.D0
      DNB(L)=0.D0
      TRB(L)=1.D0
      ENC(L)=0.D0
      DNC(L)=0.D0
      TRC(L)=1.D0
      L=L-1
      IF(L.GT.0) GO TO 205
  210 CONTINUE
      L=L+1
      UFLB(L,K)=BG
      TRUFLB(L)=TRUFLB(L)+BG
      IF(L.LT.NLP) GO TO 210
      TOTLZF(1)=TOTLZF(1)+BG
      TOTLZF(2)=TOTLZF(2)+BG
      TOTLZF(3)=TOTLZF(3)+BG
      GO TO 200
  215 CONTINUE
  220 CONTINUE
      ITL=ITLT(L)+ITK0
      BTOP=PLANCK(ITL)-(PLANCK(ITL)-PLANCK(ITL+1))*WTLT(L)
      ITL=ITLB(L)+ITK0
      BBOT=PLANCK(ITL)-(PLANCK(ITL)-PLANCK(ITL+1))*WTLB(L)
      TAUAG=TRGXLK(L,K)
      TAUAP=TRCALK(L,K)+TRAALK(L,K)+TRBALK(L,K)+TRDALK(L,K)+TRVALK(L,K)
      TAUAX=TAUAG+TAUAP

C               Optically thin limit emission/transmission approximation
C               --------------------------------------------------------

      IF(TAUAX.LT.1.D-04) THEN
      TAUBX=TAUAX+TAUAX
      TAUCX=10.D0*TAUAX
      BBAR=0.5D0*(BTOP+BBOT)
      TRA(L)=1.D0-TAUAX
      ENA(L)=BBAR*TAUAX
      DNA(L)=ENA(L)
      DNACUM=DNACUM*TRA(L)+DNA(L)
      TRB(L)=1.D0-TAUBX
      ENB(L)=BBAR*TAUBX
      DNB(L)=ENB(L)
      DNBCUM=DNBCUM*TRB(L)+DNB(L)
      TRC(L)=1.D0-TAUCX
      ENC(L)=BBAR*TAUCX
      DNC(L)=ENC(L)
      DNCCUM=DNCCUM*TRC(L)+DNC(L)
      GO TO 230
      ENDIF

C                     TAUB absorber-dependent extinction path adjustment
C                     --------------------------------------------------

      PLBN=PLB(L)
      TAUBG=TAUAG+TAUAG
      TAUCG=10.D0*TAUAG
      IF(TAUAG.GT.0.1) THEN
         IF(IMOL.EQ.1) THEN
         IF(PLBN.GT.250.0) THEN
      F=0.75
      IF(TAUAG.LT.3.0) F=0.92-0.053*TAUAG
      TAUBG=TAUBG*F
         ELSE
      F=0.70
      IF(TAUAG.LT.2.5) F=0.90-0.073*TAUAG
      TAUBG=TAUBG*F
      ENDIF
      ENDIF
         IF(IMOL.EQ.2) THEN
         IF(PLBN.GT.250.0) THEN
      F=0.58
      IF(TAUAG.LT.3.5) F=0.93-0.097*TAUAG
      TAUBG=TAUBG*F
         ELSE
      F=0.70
      IF(TAUAG.LT.3.5) F=0.92-0.062*TAUAG
      TAUBG=TAUBG*F
      ENDIF
      ENDIF
         IF(IMOL.EQ.3) THEN
         IF(PLBN.GT.250.0) THEN
      F=0.95
      IF(TAUAG.LT.0.5) F=0.99-0.016*TAUAG
      TAUBG=TAUBG*F
         ELSE
      F=0.75
      IF(TAUAG.LT.3.7) F=0.97-0.060*TAUAG
      TAUBG=TAUBG*F
      ENDIF
      ENDIF
      ENDIF

C                     TAUC absorber-dependent extinction path adjustment
C                     --------------------------------------------------

      IF(TAUAG.GT.0.01) THEN
         IF(IMOL.EQ.1) THEN
         IF(PLBN.GT.250.0) THEN
      F=0.65
      IF(TAUAG.LT.0.37) F=0.96-0.67*TAUAG
      TAUCG=TAUCG*F
         ELSE
      F=0.50
      IF(TAUAG.LT.0.47) F=0.87-0.71*TAUAG
      TAUCG=TAUCG*F
      ENDIF
      ENDIF
         IF(IMOL.EQ.2) THEN
         IF(PLBN.GT.250.0) THEN
      F=0.50
      IF(TAUAG.LT.0.75) F=0.95-0.32*TAUAG
      TAUCG=TAUCG*F
         ELSE
      F=0.50
      IF(TAUAG.LT.0.70) F=0.90-0.59*TAUAG
      TAUCG=TAUCG*F
      ENDIF
      ENDIF
         IF(IMOL.EQ.3) THEN
         IF(PLBN.GT.250.0) THEN
      F=0.95
      IF(TAUAG.LT.0.50) F=0.98-.039*TAUAG
      TAUCG=TAUCG*F
         ELSE
      F=0.75
      IF(TAUAG.LT.0.70) F=0.98-0.29*TAUAG
      TAUCG=TAUCG*F
      ENDIF
      ENDIF
      ENDIF

      TAUBP=TAUAP+TAUAP
      TAUCP=10.D0*TAUAP
      TAUA=TAUAG+TAUAP
      TAUB=TAUBG+TAUBP
      TAUC=TAUCG+TAUCP

      IF(L.EQ.LTOPCL.AND.KCLDEM.EQ.1) GO TO 225

      BDIF=BBOT-BTOP
      BBTA=BDIF/TAUA
      BBTB=BDIF/TAUB
      BBTC=BDIF/TAUC

C            Optically thick limit non-scattering emission approximation
C            -----------------------------------------------------------

      IF(TAUA.GT.9.D0) THEN
      TRA(L)=0.D0
      TRB(L)=0.D0
      TRC(L)=0.D0
      ENA(L)=BTOP+BBTA
      ENB(L)=BTOP+BBTB
      ENC(L)=BTOP+BBTC
      DNA(L)=BBOT-BBTA
      DNB(L)=BBOT-BBTB
      DNC(L)=BBOT-BBTC
      DNACUM=BBOT-BBTA
      DNBCUM=BBOT-BBTB
      DNCCUM=BBOT-BBTC
      GO TO 230
      ENDIF

      IF(TAUA.LT.0.5D0) THEN
      TAUA2=TAUA*TAUA
      TRANA=1.D0-TAUA+(0.5D0-R6*TAUA+R24*TAUA2)*TAUA2
      ELSE
      TRANA=EXP(-TAUA)
      ENDIF
      IF(TAUB.LT.0.5D0) THEN
      TAUB2=TAUB*TAUB
      TRANB=1.D0-TAUB+(0.5D0-R6*TAUB+R24*TAUB2)*TAUB2
      ELSE
      TRANB=EXP(-TAUB)
      ENDIF
      IF(TAUC.LT.0.5D0) THEN
      TAUC2=TAUC*TAUC
      TRANC=1.D0-TAUC+(0.5D0-R6*TAUC+R24*TAUC2)*TAUC2
      ELSE
      TRANC=EXP(-TAUC)
      ENDIF

      TRA(L)=TRANA
      ENA(L)=BTOP+BBTA-(BBOT+BBTA)*TRANA
      DNA(L)=BBOT-BBTA-(BTOP-BBTA)*TRANA
      DNACUM=DNACUM*TRANA+DNA(L)
      TRB(L)=TRANB
      ENB(L)=BTOP+BBTB-(BBOT+BBTB)*TRANB
      DNB(L)=BBOT-BBTB-(BTOP-BBTB)*TRANB
      DNBCUM=DNBCUM*TRANB+DNB(L)
      TRC(L)=TRANC
      ENC(L)=BTOP+BBTC-(BBOT+BBTC)*TRANC
      DNC(L)=BBOT-BBTC-(BTOP-BBTC)*TRANC
      DNCCUM=DNCCUM*TRANC+DNC(L)
      GO TO 230

C                          ---------------------------------------------
C                          Top-cloud multiple scattering corrections for
C                          emitted, transmitted, and reflected radiances
C                          and fluxes at the top-cloud (L=LTOPCL) level.
C                          ---------------------------------------------

  225 CONTINUE
      TRA(L)=EXP(-TAUAG-TAUAP*FDXTCK(3,K))
      TRB(L)=EXP(-TAUBG-TAUBP*FDXTCK(2,K))
      TRC(L)=EXP(-TAUCG-TAUCP*FDXTCK(1,K))
      DEC=C*DNCCUM*RIJTCK(1,K)+B*DNBCUM*RIJTCK(2,K)+A*DNACUM*RIJTCK(3,K)
      DEB=C*DNCCUM*RIJTCK(2,K)+B*DNBCUM*RIJTCK(4,K)+A*DNACUM*RIJTCK(5,K)
      DEA=C*DNCCUM*RIJTCK(3,K)+B*DNBCUM*RIJTCK(5,K)+A*DNACUM*RIJTCK(6,K)
      ALBTCK(1,K)=C*RIJTCK(1,K)+B*RIJTCK(2,K)+A*RIJTCK(3,K)
      ALBTCK(2,K)=C*RIJTCK(2,K)+B*RIJTCK(4,K)+A*RIJTCK(5,K)
      ALBTCK(3,K)=C*RIJTCK(3,K)+B*RIJTCK(5,K)+A*RIJTCK(6,K)
      COALB1=1.D0-ALBTCK(1,K)
      COALB2=1.D0-ALBTCK(2,K)
      COALB3=1.D0-ALBTCK(3,K)
      TAUA=TAUAG+TAUAP*FEMTCK(3,K)
      TAUB=TAUBG+TAUBP*FEMTCK(2,K)
      TAUC=TAUCG+TAUCP*FEMTCK(1,K)
      TRANA=EXP(-TAUA)
      TRANB=EXP(-TAUB)
      TRANC=EXP(-TAUC)
      BDIF=BBOT-BTOP
      BBTA=BDIF/TAUA
      BBTB=BDIF/TAUB
      BBTC=BDIF/TAUC
      ENA(L)=(BTOP+BBTA-(BBOT+BBTA)*TRANA)*COALB3
      DNA(L)=(BBOT-BBTA-(BTOP-BBTA)*TRANA)*COALB3
      DNACUM=DNACUM*TRA(L)+DNA(L)
      ENB(L)=(BTOP+BBTB-(BBOT+BBTB)*TRANB)*COALB2
      DNB(L)=(BBOT-BBTB-(BTOP-BBTB)*TRANB)*COALB2
      DNBCUM=DNBCUM*TRB(L)+DNB(L)
      ENC(L)=(BTOP+BBTC-(BBOT+BBTC)*TRANC)*COALB1
      DNC(L)=(BBOT-BBTC-(BTOP-BBTC)*TRANC)*COALB1
      DNCCUM=DNCCUM*TRC(L)+DNC(L)
      ENC(L)=ENC(L)+DEC
      ENB(L)=ENB(L)+DEB
      ENA(L)=ENA(L)+DEA
  230 CONTINUE
      FDNABC=A*DNACUM+B*DNBCUM+C*DNCCUM
      TRDFLB(L)=TRDFLB(L)+FDNABC
      DFLB(L,K)=FDNABC
      L=L-1
      IF(L.GT.0) GO TO 220

C             Old form of scattering correction is skipped when KCLDEM=1
C             ----------------------------------------------------------

      IF(KCLDEM.EQ.1) GO TO 235
      IF(LTOPCL.LT.1) GO TO 235
      ENA(LTOPCL)=ENA(LTOPCL)*(1.0-TRCTCA(K))+TRCTCA(K)*DFLB(LTOPCL+1,K)
      ENB(LTOPCL)=ENB(LTOPCL)*(1.0-TRCTCA(K))+TRCTCA(K)*DFLB(LTOPCL+1,K)
      ENC(LTOPCL)=ENC(LTOPCL)*(1.0-TRCTCA(K))+TRCTCA(K)*DFLB(LTOPCL+1,K)
  235 CONTINUE

C     ------------------------------------------------------------------
C                                         SURFACE LAYER FLUX COMPUTATION
C     with TAUSL,FTAUSL=0 defaults, surface layer calculation is skipped
C     ------------------------------------------------------------------

      L=1
      DFSL(K)=FDNABC
!!    TAUA=TAUSL(K)+FTAUSL(K)
!!    IF(TAUA.GT.1.D-06) GO TO 240
      BG=BG+FDNABC*TRGALB(K)
      UNA=BG
      UNB=BG
      UNC=BG
      FUNABC=BG
!!    GO TO 245
!!240 CONTINUE
!!    ITS=TSL
!!    WTS=TSL-ITS
!!    WTS1=1.D0-WTS
!!    ITS=ITS-ITPFT0
!!    ITS=ITS+(K-1)*ITNEXT
!!    BS=PLANCK(ITS)*WTS1+PLANCK(ITS+1)*WTS
!!    TA=EXP(-TAUA)
!!    TB=TA*TA
!!    TC=(TB*TB*TA)**2
!!    DNA(1)=(DNA(1)-BS)*TA+BS
!!    DNB(1)=(DNB(1)-BS)*TB+BS
!!    DNC(1)=(DNC(1)-BS)*TC+BS
!!    FDNABC=A*DNA(1)+B*DNB(1)+C*DNC(1)
!!    BG=BGFEMT(K)+FDNABC*TRGALB(K)
!!    UNA=(BG-BS)*TA+BS
!!    UNB=(BG-BS)*TB+BS
!!    UNC=(BG-BS)*TC+BS
!!    FUNABC=A*UNA+B*UNB+C*UNC
!!    BSP=PLANCK(ITS+1)*WTS1+PLANCK(ITS+2)*WTS
!!    BSM=PLANCK(ITS-1)*WTS1+PLANCK(ITS  )*WTS
!!    SLABS=1.D0-A*TA-B*TB-C*TC
!!    TRSLTS=TRSLTS+(BSP-BSM)*SLABS
!!    TRSLTG=TRSLTG+BGFEMD(K)*SLABS
!!    TRSLBS=TRSLBS+BS*SLABS

C     ------------------------------------------------------------------
C                                                UPWARD FLUX COMPUTATION
C     ------------------------------------------------------------------

  245 CONTINUE
      TRUFLB(L)=TRUFLB(L)+FUNABC
      UFLB(L,K)=FUNABC

C       ----------------------------------------------------------------
C       At top-cloud level, compute component of upwelling flux relected
C       downward by cloud bottom and add to downwelling flux below cloud
C       ----------------------------------------------------------------

      IF(L.EQ.LTOPCL.AND.KCLDEM.EQ.1) THEN
      DEC=C*UNC*RIJTCK(1,K)+B*UNB*RIJTCK(2,K)+A*UNA*RIJTCK(3,K)
      DEB=C*UNC*RIJTCK(2,K)+B*UNB*RIJTCK(4,K)+A*UNA*RIJTCK(5,K)
      DEA=C*UNC*RIJTCK(3,K)+B*UNB*RIJTCK(5,K)+A*UNA*RIJTCK(6,K)
      N=L
  250 CONTINUE
      DNA(N)=DNA(N)+DEA
      DNB(N)=DNB(N)+DEB
      DNC(N)=DNC(N)+DEC
      DDFLUX=A*DEA+B*DEB+C*DEC
      TRDFLB(N)=TRDFLB(N)+DDFLUX
      DFLB(N,K)=DFLB(N,K)+DDFLUX
      IF(N.LT.2) GO TO 255
      N=N-1
      DEA=DEA*TRA(N)
      DEB=DEB*TRB(N)
      DEC=DEC*TRC(N)
      GO TO 250
  255 CONTINUE
      ENDIF
      UNA=UNA*TRA(L)+ENA(L)
      UNB=UNB*TRB(L)+ENB(L)
      UNC=UNC*TRC(L)+ENC(L)
      FUNABC=A*UNA+B*UNB+C*UNC

      L=L+1
      IF(L.LT.NLP) GO TO 245
      IF(K.EQ.1) THEN
      TRUFTW=FUNABC
      TRDFGW=TRDFLB(1)
      TRUFGW=BG
      WINDZF(1)=UNA
      WINDZF(2)=UNB
      WINDZF(3)=UNC
      ENDIF

      TRUFLB(NLP)=TRUFLB(NLP)+FUNABC
      UFLB(NLP,K)=FUNABC
      UFSL(K)=UFLB(1,K)
      TOTLZF(1)=TOTLZF(1)+UNA
      TOTLZF(2)=TOTLZF(2)+UNB
      TOTLZF(3)=TOTLZF(3)+UNC

      IMOL=2
      IF(K.LT.13) IMOL=1
      IF(K.GT.24) IMOL=3
      GO TO 200
  300 CONTINUE

      DO 310 L=1,NLP
      TRNFLB(L)=TRUFLB(L)-TRDFLB(L)
  310 CONTINUE
      DO 320 L=1,NL
      TRFCRL(L)=TRNFLB(L+1)-TRNFLB(L)
  320 CONTINUE

C**** Window region and spectr. integrated total flux diagnostics
      PFW=max( 1.0001d-2, 10.D0*TRUFTW )
      IF(PFW.GT.629.999d0) PFW=629.999d0
      IPF=PFW
      IF(IPF.LT.10) GO TO 330
      DPF=PFW-IPF
      IPF=IPF+180
      GO TO 350
  330 CONTINUE
      PFW=10.D0*PFW
      IPF=PFW
      IF(IPF.LT.10) GO TO 340
      DPF=PFW-IPF
      IPF=IPF+90
      GO TO 350
  340 CONTINUE
      PFW=10.D0*PFW
      IPF=PFW
C     IF(IPF.LT.1) IPF=1
  350 CONTINUE
      BTEMPW=TKPFW(IPF)+DPF*(TKPFW(IPF+1)-TKPFW(IPF))

      DO 390 II=1,3
      PFW=max( 1.d0, TOTLZF(II) )
      IF(PFW.GT.899.999d0) PFW=899.999d0
      IPF=PFW
      DPF=PFW-IPF
C     IF(IPF.LT.1) IPF=1
C     IF(IPF.GT.899) IPF=899
      TOTLZT(II)=TKPFT(IPF)+DPF*(TKPFT(IPF+1)-TKPFT(IPF))

      PFW=max( 1.0001d-2, 10.D0*WINDZF(II) )
      IF(PFW.GT.629.999d0) PFW=629.999d0
      IPF=PFW
      IF(IPF.LT.10) GO TO 360
      DPF=PFW-IPF
      IPF=IPF+180
      GO TO 380
  360 CONTINUE
      PFW=10.D0*PFW
      IPF=PFW
      IF(IPF.LT.10) GO TO 370
      DPF=PFW-IPF
      IPF=IPF+90
      GO TO 380
  370 CONTINUE
      PFW=10.D0*PFW
      IPF=PFW
  380 CONTINUE
      WINDZT(II)=TKPFW(IPF)+DPF*(TKPFW(IPF+1)-TKPFW(IPF))
  390 CONTINUE
      RETURN
      END SUBROUTINE THERML

      SUBROUTINE SOLAR0

      !INCLUDE 'BR00B.COM'

      DIMENSION NMKWAV(17),LORDER(16)
      DATA NMKWAV/ 200, 360, 770, 795, 805, 810, 860,1250,1500,1740,
     +            2200,3000,3400,3600,3800,4000,9999/
      DATA LORDER/15,14, 8, 7, 6, 5, 13, 12, 4, 3, 2, 1, 11, 10, 9,16/

      N2=1
      DO 1 I=1,30
      N2=N2*2
      DBLN(I)=N2
    1 CONTINUE
      DO 2 I=1,16
      SKDKS0(I)=DKS0(I)
      NORDER(I)=LORDER(I)
      NMWAVA(I)=NMKWAV(I)
      NMWAVB(I)=NMKWAV(I+1)
    2 CONTINUE

      TCLMIN=MIN(TAUIC0,TAUWC0)

      RETURN
      END SUBROUTINE SOLAR0

      SUBROUTINE SOLARM

      !INCLUDE 'BR00B.COM'

C     ------------------------------------------------------------------
C     SOLARM Returns:
C                      SRDFLB   Solar downward flux at layer bottom edge
C                      SRUFLB   Solar  upward  flux at layer bottom edge
C                      SRNFLB   Solar net downward flux  (in Watts/m**2)
C                      SRFHRL   Solar heating rate/layer (in Watts/m**2)
C                      FSRNFG   Solar flux abs at ground by surface-type
C                               (see explanatory note at end of SOLARM)
C               Also:
C                      TOA:   SRIVIS SROVIS PLAVIS  SRINIR SRONIR PLANIR
C                      BOA:   SRDVIS SRUVIS ALBVIS  SRDNIR SRUNIR ALBNIR
C                      ATM:   SRTVIS SRRVIS SRAVIS  SRTNIR SRRNIR SRANIR
C                             SRXVIS SRXNIR (Direct beam only at ground)
C
C           Spectral:  (by k-distribution/pseudo-spectral) breakdown:
C                      SKDFLB   Solar downward flux at layer bottom edge
C                      SKUFLB   Solar  upward  flux at layer bottom edge
C                      SKNFLB   Solar net downward flux  (in Watts/m**2)
C                      SKFHRL   Solar heating rate/layer (in Watts/m**2)
C
C                      SRKALB   Planetary albedo (by spectral breakdown)
C                      SRKINC   Incident fluxedo (by spectral breakdown)
C                      SRKGAX   Direct  k-d flux absorbed by ground-type
C                      SRKGAD   Diffuse k-d flux absorbed by ground-type
C     ------------------------------------------------------------------
C     Remarks:
C              NORMS0=1 Incident (TOA) Solar flux normalized to equal S0
C                       (COSZ dependence included in calculated results)
C                        The returned solar fluxes have to be multiplied
C                       by COSZ to yield actual atmospheric heating rate
C
C              NMKWAV   Spectral/k-distribution subdivisions are nominal
C                       (due to spectral trading of absorption features)
C
C                 VIS   Designates solar visible wavelengths (   <770nm)
C                 NIR   Designates solar near-IR wavelengths (770>   nm)
C                       VIS comprises .53 of S0, NIR comprises .47 of S0
C     ------------------------------------------------------------------
C
CCC   SAVE
C     ------------------------------------------------------------------
C         Fractional solar flux k-distribution/pseudo-spectral intervals
C
C         KSLAM=    1     1     2     2     5     5     5     5
C             K=    1     2     3     4     5     6     7     8
CCC   DATA DKS0/ .010, .030, .040, .040, .040, .002, .004, .013,
C         KSLAM=    1     1     1     3     4     6     6     1
C             K=    9    10    11    12    13    14    15    16
CCC  +           .002, .003, .003, .072, .200, .480, .050, .011/
C
C     ------------------------------------------------------------------
C     The nominal spectral order for k-dist/pseudo-spectral intervals is
C     (WavA and WavB designate approximate spectral interval boundaries)
C
C             L=   12    11    10     9     6     5     4     3
C     WavA (nm)= 3000  2200  1740  1500   810   805   795   770
C     WavB (nm)= 3400  3000  2200  1740   860   810   805   795
C             K=    1     2     3     4     5     6     7     8
C     DATA DKS0/ .010, .030, .040, .040, .040, .002, .004, .013,
C
C             L=   15    14    13     8     7     2     1    16
C     WavA (nm)= 3800  3500  3400  1250   860   360   200  4000
C     WavB (nm)= 4000  3800  3600  1500  1250   770   360  9999
C             K=    9    10    11    12    13    14    15    16
C    +           .002, .003, .003, .072, .200, .480, .050, .011/
C
C     ------------------------------------------------------------------
C     6 spectral intervals overlap the 16 solar k-distribution intervals
C
C     Cloud and aerosol Mie scattering parameters (also surface albedos)
C     are averaged over these spectral intervals. These intervals are in
C     reverse spectral order. Thus spectral interval 6 refers to visible
C     (VIS) wavelengths, intervals 1-5 refer to nearIR (NIR) wavelengths
C     KSLAM designates the spectral interval of first 14 k-distributions
C     (K=15 for UV ozone absorption refers to (VIS) spectral interval 6)
C     (K=16 represents strong absorbing spectral regions via interval 1)
C
C     The nominal Mie scattering spectral band subdivisions are:
C
C                         -------------NIR------------      VIS
C                     L=    1     2     3     4     5        6
C             WavA (nm)=  2200  1500  1250   860   770      300
C             WavB (nm)=  4000  2200  1500  1250   860      770
C
C     ------------------------------------------------------------------

      DIMENSION COLEXT(6),COLSCT(6),COLGCB(6)   ! ,ALLGCB(6)

C                            -------------------------------------------
C                            NO2, O3 Chappuis Band, Rayleigh, parameters
C                            -------------------------------------------
      PARAMETER(XCMNO2=5.465,XCMO3=.0399623,TOTRAY=0.000155)

      S0COSZ=S0
      IF(NORMS0.EQ.0) S0COSZ=S0*COSZ

      DO 10 N=1,NLP
      SRNFLB(N)=0.D0
      SRDFLB(N)=0.D0
      SRUFLB(N)=0.D0
      SRFHRL(N)=0.D0
      SKDFLB(N,16)=0.D0
      SKUFLB(N,16)=0.D0
      SRKALB(N)=0.D0  ! for WRITER only
   10 CONTINUE
      dblext=0. ; dblsct=0. ; dblgcb=0. ; dblpi0=0. ! for writer only
      skdflb=0. ; sknflb=0. ; skuflb=0.             ! for writer only
      skfhrl=0. ; srkgax=0. ; srkgad=0.             ! for writer only
C                     TOA solar flux VIS/NIR subdivision
C                     (incident, outgoing, plane albedo)
C                     ----------------------------------
      SRIVIS=0.D0
      SROVIS=0.D0
      PLAVIS=1.D0
      SRINIR=0.D0
      SRONIR=0.D0
      PLANIR=1.D0
C                     BOA solar flux VIS/NIR subdivision
C                     (incident, upward, surface albedo)
C                     ----------------------------------
      SRDVIS=0.D0
      SRUVIS=0.D0
      ALBVIS=1.D0
      SRDNIR=0.D0
      SRUNIR=0.D0
      ALBNIR=1.D0
C                     Fractional atmos only flux VIS/NIR subdivision
C                     (fractions reflected, transmitted, & absorbed)
C                     ----------------------------------------------
      SRRVIS=1.D0
      SRTVIS=0.D0
      SRAVIS=0.D0
      SRRNIR=0.D0
      SRTNIR=0.D0
      SRANIR=0.D0
C                     Direct beam, fractional S0 VIS/NIR subdivision
C                     ----------------------------------------------
      SRXVIS=0.D0
      SRXNIR=0.D0
!obso SRSLHR=0.D0
C                     Ground surface absorbed solar flux subdivision
C                     according to 4 fractional surface-type albedos
C                     ----------------------------------------------
      FSRNFG(1)=0.D0
      FSRNFG(2)=0.D0
      FSRNFG(3)=0.D0
      FSRNFG(4)=0.D0

      IF(COSZ.LT.0.001D0) RETURN
      COSMAG=35.D0/SQRT(1224.D0*COSZ*COSZ+1.D0)
      SECZ=1.D0/COSZ

      TAURAY=TOTRAY*FRAYLE

      DO 90 K=1,6
      RTAU=1.D-10
      IF(K.EQ.6) RTAU=TAURAY
      COLEXT(K)=0.D0
      COLSCT(K)=0.D0
      COLGCB(K)=0.D0
      DO 30 L=1,NL
      RTAUL=RTAU*(PLB(L)-PLB(L+1))
      SUMEXT=RTAUL+SRCEXT(L,K)+SRAEXT(L,K)+SRBEXT(L,K)
     +                        +SRDEXT(L,K)+SRVEXT(L,K)
      SUMSCT=RTAUL+SRCSCT(L,K)+SRASCT(L,K)+SRBSCT(L,K)
     +                        +SRDSCT(L,K)+SRVSCT(L,K)
      SUMGCB=      SRCSCT(L,K)*SRCGCB(L,K)+SRASCT(L,K)*SRAGCB(L,K)
     +            +SRBSCT(L,K)*SRBGCB(L,K)+SRDSCT(L,K)*SRDGCB(L,K)
     +            +SRVSCT(L,K)*SRVGCB(L,K)
      DBLEXT(L,K)=SUMEXT
      DBLSCT(L,K)=SUMSCT
      DBLGCB(L,K)=SUMGCB/(SUMSCT+1.D-10)
      DBLPI0(L,K)=SUMSCT/(SUMEXT+1.D-10)
      COLEXT(K)=COLEXT(K)+DBLEXT(L,K)
      COLSCT(K)=COLSCT(K)+DBLSCT(L,K)
      COLGCB(K)=COLGCB(K)+DBLSCT(L,K)*DBLGCB(L,K)
   30 CONTINUE
      COLGCB(K)=COLGCB(K)/(COLSCT(K)+1.D-10)

      IF(KANORM.GT.0) THEN
C     -----------------------------------------------------------------
C     KANORM (default = 0)  Option to renormalize aerosol column albedo
C                           to make column albedo less dependent on the
C                           number of model layers due to SGP treatment
C
C                           KANORM=1  aerosol column only is normalized
C
C                           KANORM=2  aerosol plus ground is normalized
C                                     with Tau equivalent ground albedo
C                                     ---------------------------------
      COLPFG=COLGCB(K)
      SURFBB=SRBALB(K)
      TAUSBB=0.D0
      IF(KANORM.GT.1) CALL GTSALB(XXG,XXT,SURX,SURFBB,COLPFG,TAUSBB,2)
      DBLEXT(NLP,K)=TAUSBB
      ALLTAU=TAUSBB+COLEXT(K)
      CALL SGPGXG(COSZ,ALLTAU,COLPFG,SGPG)
cc    ALLGCB(K)=SGPG
      DO 40 L=1,NL
      DBLGCB(L,K)=SGPG
   40 CONTINUE
      ELSE

      DO 50 L=1,NL
      TAULAY=DBLEXT(L,K)
      GCBLAY=DBLGCB(L,K)
      CALL SGPGXG(COSZ,TAULAY,GCBLAY,SGPG)
      DBLGCB(L,K)=SGPG
   50 CONTINUE
      ENDIF

      IF(LTOPCL.LT.1) GO TO 90
      RTAU=1.D-10
      IF(K.EQ.6) RTAU=TAURAY
      COLEXT(K)=0.D0
      COLSCT(K)=0.D0
      COLGCB(K)=0.D0
      DO 60 L=1,NL
      IF(SRCEXT(L,K).LT.TCLMIN) GO TO 60
      RTAUL=RTAU*(PLB(L)-PLB(L+1))
      SUMEXT=RTAUL+SRCEXT(L,K)+SRAEXT(L,K)+SRBEXT(L,K)
     +                        +SRDEXT(L,K)+SRVEXT(L,K)
      SUMSCT=RTAUL+SRCSCT(L,K)+SRASCT(L,K)+SRBSCT(L,K)
     +                        +SRDSCT(L,K)+SRVSCT(L,K)
      SUMGCB=      SRCSCT(L,K)*SRCGCB(L,K)+SRASCT(L,K)*SRAGCB(L,K)
     +            +SRBSCT(L,K)*SRBGCB(L,K)+SRDSCT(L,K)*SRDGCB(L,K)
     +            +SRVSCT(L,K)*SRVGCB(L,K)
      DBLEXT(L,K)=SUMEXT
      DBLSCT(L,K)=SUMSCT
      DBLGCB(L,K)=SUMGCB/(SUMSCT+1.D-10)
      DBLPI0(L,K)=SUMSCT/(SUMEXT+1.D-10)
      COLEXT(K)=COLEXT(K)+DBLEXT(L,K)
      COLSCT(K)=COLSCT(K)+DBLSCT(L,K)
      COLGCB(K)=COLGCB(K)+DBLSCT(L,K)*DBLGCB(L,K)
   60 CONTINUE
      COLGCB(K)=COLGCB(K)/(COLSCT(K)+1.D-10)

      IF(KCNORM.GT.0) THEN
C     -----------------------------------------------------------------
C     KCNORM (default = 0)  Option to renormalize  cloud  column albedo
C                           to make column albedo less dependent on the
C                           number of model layers due to SGP treatment
C
C                           KCNORM=1    cloud column only is normalized
C
C                           KCNORM=2    cloud plus ground is normalized
C                                     with Tau equivalent ground albedo
C                                     ---------------------------------
      COLPFG=COLGCB(K)
      SURFBB=SRBALB(K)
      TAUSBB=0.D0
      IF(KCNORM.GT.1) CALL GTSALB(XXG,XXT,SURX,SURFBB,COLPFG,TAUSBB,2)
      DBLEXT(NLP,K)=TAUSBB
      ALLTAU=TAUSBB+COLEXT(K)
      CALL SGPGXG(COSZ,ALLTAU,COLPFG,SGPG)
cc    ALLGCB(K)=SGPG
      DO 70 L=1,NL
      IF(SRCEXT(L,K).LT.TCLMIN) GO TO 70
      DBLGCB(L,K)=SGPG
   70 CONTINUE
      ELSE
      DO 80 L=1,NL
      IF(SRCEXT(L,K).LT.TCLMIN) GO TO 80
      TAULAY=DBLEXT(L,K)
      GCBLAY=DBLGCB(L,K)
      CALL SGPGXG(COSZ,TAULAY,GCBLAY,SGPG)
      DBLGCB(L,K)=SGPG
   80 CONTINUE
      ENDIF
   90 CONTINUE

      K = 0
  300 CONTINUE
      K = K+1

      KLAM=KSLAM(K)
      DKS0X=DKS0(K)*S0COSZ
      RBNB=SRBALB(KLAM)
      RBNX=SRXALB(KLAM)
      RCNB=0.D0
      RCNX=0.D0
      SRKINC(K)=DKS0X

      N = 0
  200 CONTINUE
      N = N+1

      SRB(N)=RBNB
      SRX(N)=RBNX
      TLN=TLM(N)
      PLN=PL(N)
      ULN=ULGAS(N,1)

C     Select parameterized k-distribution gas absorption by H2O, O2, CO2
C     ------------------------------------------------------------------

      GO TO (101,102,103,104,105,106,107,108,109,110,111,112,113,114),K
 101  CONTINUE
C--------K=6-------H2O       DS0=.01
      TERMA=(35.66+TLN*(.0416-.0004622*TLN+.001057*PLN))*(1.+.04286*PLN)
      TERMB=(1.+.00171*ULN)*(1.+PLN*(189.088+.1316*PLN))
      TAU1 =TERMA/TERMB
      IF(TAU1.GT.0.02343) TAU1=0.02343
      TAU=TAU1*ULN
      GO TO 120
 102  CONTINUE
C--------K=5-------H2O       DS0=.03
      TERMA=(2.792+TLN*(.0914-.0002848*TLN+.0003395*PLN))
     +     *(1.+.02964*PLN)
      TERMB=(1.0+.000657*ULN)*(1.+PLN*(240.70+.13847*PLN))
      TAU1 =TERMA/TERMB
      IF(TAU1.GT.0.00520) TAU1=0.00520
      TAU=TAU1*ULN
      GO TO 120
 103  CONTINUE
C--------K=4-------H2O       DS0=.04
      TERMA=(.4768+.467E-04*PLN*TLN)*(1.+TLN*(.00191-.719E-05*TLN))
      TERMB=(1.+.717E-04*ULN)*(1.+PLN*(130.56+.0876*PLN))/(1.+.0266*PLN)
      TAU1 =TERMA/TERMB
      IF(TAU1.GT.0.00150) TAU1=0.0015
      TAU=TAU1*ULN
      GO TO 120
 104  CONTINUE
C--------K=3-------H2O       DS0=.04
      TERMA=(.000247*TLN-.091+PLN*(.00035+.78E-06*TLN))*(1.+.2847*PLN)
      TERMB=(1.+.2066E-04*ULN)*(1.+PLN*(137.17+.16132*PLN))
      TAU  =(TERMA/TERMB)*ULN
      GO TO 120
 105  CONTINUE
C--------K=2-------H2O       DS0=.04
      TERMA=(PLN*(1.974/TLN+.0001117*TLN)-10.713)*(1.+.005788*TLN)
     +     *(1.+.001517*PLN)
      TERMB=(1.+.3218E-04*ULN)*(1.+PLN*(863.44+.2048*PLN))
      TAU  =(TERMA/TERMB)*ULN
      GO TO 120
 106  CONTINUE
C--------K=4-------O2        DS0=.002
      ULN=ULGAS(N,4)
      TERMA=(.2236E-05-.1181E-09*TLN)*(1.+PLN*(.6364E-05*PLN+.001168))
      TERMB=1.+.1521E-05*ULN
      TAU  =(TERMA/TERMB)*ULN
      GO TO 120
 107  CONTINUE
C--------K=3-------O2        DS0=.004
      ULN=ULGAS(N,4)
      TERMA=(.3179E-06-.9263E-11*TLN)*(1.+PLN*(.8832E-05*PLN+.0005292))
      TERMB=1.+.1968E-06*ULN
      TAU  =(TERMA/TERMB)*ULN
      GO TO 120
 108  CONTINUE
C--------K=2-------O2        DS0=.013
      ULN=ULGAS(N,4)
      TERMA=(.2801E-07-.1638E-12*TLN)*(1.+PLN*(.1683E-04*PLN-.001721))
      TERMB=1.+.8097E-07*ULN
      TAU  =(TERMA/TERMB)*ULN
      GO TO 120
 109  CONTINUE
C--------K=4-------CO2       DS0=.002
      ULN=ULGAS(N,2)
      TERMA=(50.73-.03155*TLN-PLN*(.5543+.00091*TLN))*(1.-.1004*PLN)
      TERMB=(1.+.006468*ULN)*(1.+PLN*(49.51+.8285*PLN))
      TAU  =(TERMA/TERMB)*ULN
      IF(PLN.LT.175.0) TAU=(.00018*PLN+0.00001)*ULN
      GO TO 120
 110  CONTINUE
C--------K=3-------CO2       DS0=.003
      ULN=ULGAS(N,2)
      TERMA=(1.+.01319*TLN)*(PLN*(.008001*ULN+.4589E-03)-.8396*ULN)
      TERMB=ULN*(PLN+295.7+1.967*ULN)+.15126*PLN
      TAU  =(TERMA/TERMB)*ULN
      GO TO 120
 111  CONTINUE
C--------K=2-------CO2       DS0=.003
      ULN=ULGAS(N,2)
      TERMA=(1.+.02257*TLN)*(PLN*(.002295*ULN-.5489E-04)-.7571*ULN)
      TERMB=ULN*(PLN+803.9+2.477*ULN)-.09899*PLN
      TAU  =(TERMA/TERMB)*ULN
      GO TO 120
 112  CONTINUE
      TAU=0.D0
      GO TO 120
 113  CONTINUE
      TAU=0.D0
      GO TO 120
 114  CONTINUE
      TAU=XCMNO2*ULGAS(N,5)+XCMO3*ULGAS(N,3)
 120  CONTINUE

C     With 10 doublings to get to Tau=1.0, maximum seed tau is < 1/1024.
C     ------------------------------------------------------------------

      IF(TAU.LT.0.D0) TAU=0.D0

      TAU=TAU+DBLEXT(N,KLAM)
      IF(TAU.LT.1.D-06) GO TO 180
      PIZERO=DBLSCT(N,KLAM)/TAU
      IF(PIZERO.LT.0.001D0) GO TO 180

      PFF=DBLGCB(N,KLAM)

      NDBLS=0
      PR=1.D0-PFF
      PT=1.D0+PFF
      IF(TAU.GT.0.0019531D0) THEN
      DBLS=10.D0+1.44269D0*LOG(TAU)
      NDBLS=DBLS
      TAU=TAU/DBLN(NDBLS)
      ENDIF

C     Set optically thin limit values of R,T,X using PI0 renormalization
C     ------------------------------------------------------------------

      XANB=EXP(-TAU-TAU)
      XANX=EXP(-TAU*SECZ)
      TANB=PT*XANB
      XXT=(SECZ-2.D0)*TAU
      TANX=PT*SECZ
     +    *(.5D0+XXT*(.25D0+XXT*(.0833333D0+XXT*(.0208333D0+XXT))))*XANX
      RASB=PR*(1.D0-TAU*(2.D0-2.66667D0*TAU*(1.D0-TAU)))
      XXT=(SECZ+2.D0)*TAU
      RASX=PR*SECZ
     +    *(.5D0-XXT*(.25D0-XXT*(.0833333D0-XXT*(.0208333D0-XXT))))
      BNORM=(1.D0-XANB)/(RASB+TANB)*PIZERO
      XNORM=(1.D0-XANX)/(RASX+TANX)*PIZERO
      RASB=RASB*BNORM
      RASX=RASX*XNORM
      TANB=TANB*BNORM
      TANX=TANX*XNORM

C     Compute and record R,T,X atmospheric layer doubling/adding results
C     ------------------------------------------------------------------

      IF(NDBLS.LT.1) GO TO 170
      DO 160 NN=1,NDBLS
      RARB=RASB*RASB
      RARX=XANX*RASX
      XATB=XANB+TANB
      DENOM=1.D0-RARB
      DB=(TANB+XANB*RARB)/DENOM
      DX=(TANX+RARX*RASB)/DENOM
      UB=RASB*(XANB+DB)
      UX=RARX+RASB*DX
      RASB=RASB+XATB*UB
      RASX=RASX+XATB*UX
      TANB=XANB*TANB+XATB*DB
      TANX=XANX*TANX+XATB*DX
      XANB=XANB*XANB
      XANX=XANX*XANX
  160 CONTINUE
  170 CONTINUE
      RARB=RASB*RBNB
      RARX=RASB*RBNX
      XATB=XANB+TANB
      DENOM=1.D0-RARB
      DB=(TANB+XANB*RARB)/DENOM
      DX=(TANX+XANX*RARX)/DENOM
      UB=RBNB*(XANB+DB)
      UX=RBNX*XANX+RBNB*DX
      RBNB=RASB+XATB*UB
      RBNX=RASX+XATB*UX
      XATC=XATB/(1.D0-RASB*RCNB)
      RCNX=RASX+(XANX*RCNX+TANX*RCNB)*XATC
      RCNB=RASB+RCNB*XATB*XATC
      GO TO 190
  180 CONTINUE
      RASB=0.D0
      RASX=0.D0
      TANB=0.D0
      TANX=0.D0
      XANB=EXP(-TAU-TAU)
      XANX=EXP(-TAU*SECZ)
      DX=0.D0
      UX=RBNX*XANX
      RBNB=RBNB*XANB*XANB
      RBNX=UX*XANB
      RCNB=RCNB*XANB*XANB
      RCNX=RCNX*XANX*XANB
  190 CONTINUE
      RNB(N)=RASB
      RNX(N)=RASX
      TNB(N)=TANB
      TNX(N)=TANX
      XNB(N)=XANB
      XNX(N)=XANX
      IF(N.LT.NL) GO TO 200

C     Record fluxes, spectral components at TOA, & top-layer bottom edge
C     ------------------------------------------------------------------

      SRDFLB(NLP)=SRDFLB(NLP)+DKS0X
      SRUFLB(NLP)=SRUFLB(NLP)+DKS0X*RBNX
      SRDFLB(NL)=SRDFLB(NL)+DKS0X*(XANX+DX)
      SRUFLB(NL)=SRUFLB(NL)+DKS0X*UX
      SKDFLB(NLP,K)=DKS0X
      SKUFLB(NLP,K)=DKS0X*RBNX
      SKDFLB(NL,K)=DKS0X*(XANX+DX)
      SKUFLB(NL,K)=DKS0X*UX
      RBXTOA=RBNX
      SRKALB(K)=RBNX

C     Add successively layer N (at bottom) to form upper composite layer
C     ------------------------------------------------------------------

      DO 230 M=2,NL
      N=NLP-M
      XBNB=XNB(N)
      XBNX=XNX(N)
      RBNX=RNX(N)
      IF(RBNX.GT.1.D-05) GO TO 210
      RASB=RASB*XBNB*XBNB
      TANX=TANX*XBNB
      GO TO 220
  210 CONTINUE
      RBNB=RNB(N)
      TBNB=TNB(N)
      TBNX=TNX(N)
      RARB=RASB*RBNB
      XBTB=XBNB+TBNB
      DENOM=1.D0-RARB
      TANX=TBNX*XANX+XBTB*(TANX+XANX*RBNX*RASB)/DENOM
      RASB=RBNB+XBTB*XBTB*RASB/DENOM
  220 CONTINUE
      XANX=XANX*XBNX
      RBNB=SRB(N)
      RBNX=SRX(N)
      DX=(TANX+XANX*RBNX*RASB)/(1.D0-RASB*RBNB)
      UX=RBNX*XANX+RBNB*DX
      SRUFLB(N)=SRUFLB(N)+DKS0X*UX
      SRDFLB(N)=SRDFLB(N)+DKS0X*(XANX+DX)
      SKUFLB(N,K)=DKS0X*UX
      SKDFLB(N,K)=DKS0X*(XANX+DX)
  230 CONTINUE

C     Record absorbed spectral flux at ground for surface type fractions
C     ------------------------------------------------------------------

      DO 231 I=1,4
      SRKGAX(K,I)=DKS0X*XANX*(1.D0-PRNX(KLAM,I))
      SRKGAD(K,I)=DKS0X*  DX*(1.D0-PRNB(KLAM,I))
  231 CONTINUE

      IF(K.EQ.NKSLAM) GO TO 301
      SRINIR=SRINIR+DKS0X
      SRONIR=SRONIR+DKS0X*RBXTOA
      SRDNIR=SRDNIR+SKDFLB(1,K)
      SRUNIR=SRUNIR+SKUFLB(1,K)
      SRRNIR=SRRNIR+DKS0X*RCNX
      SRTNIR=SRTNIR+DKS0X*(TANX+XANX)
      SRXNIR=SRXNIR+DKS0X*XANX
      GO TO 300

  301 CONTINUE
      SRIVIS=DKS0X
      SROVIS=DKS0X*RBXTOA
      SRDVIS=SKDFLB(1,K)
      SRUVIS=SKUFLB(1,K)
      SRRVIS=DKS0X*RCNX
      SRTVIS=DKS0X*(TANX+XANX)
      SRXVIS=SRXVIS+DKS0X*XANX

C     ------------------------------------------------------------------
C     UV absorption by O3 and O2 within solar spectral band DKS0(15)=.05
C     ------------------------------------------------------------------

      K=15
      DKS0X=DKS0(K)*S0COSZ
      SRKINC(K)=DKS0X

      N=NLP
      ATOPX=0.D0
      ATOPD=0.D0
      O3CMX=0.D0
      O3CMD=0.D0
  302 CONTINUE
      N=N-1
      O3CMX=O3CMX+COSMAG*ULGAS(N,3)
      O3CMD=O3CMD+1.90D0*ULGAS(N,3)
      CALL AO3ABS(O3CMX,ABOTX)
      CALL AO3ABS(O3CMD,ABOTD)
      AO3X(N)=(ABOTX-ATOPX)/DKS0(15)
      AO3D(N)=(ABOTD-ATOPD)/DKS0(15)
      ATOPX=ABOTX
      ATOPD=ABOTD
      IF(N.GT.1) GO TO 302
  303 CONTINUE
      O3CMX=O3CMX+1.90D0*ULGAS(N,3)
      O3CMD=O3CMD+1.90D0*ULGAS(N,3)
      CALL AO3ABS(O3CMX,ATOPX)
      CALL AO3ABS(O3CMD,ATOPD)
      AO3UXN=(ATOPX-ABOTX)/DKS0(15)
      AO3UDN=(ATOPD-ABOTD)/DKS0(15)
      AO3U(N)=XNX(N)*AO3UXN+(1.D0-XNX(N))*AO3UDN
      ABOTX=ATOPX
      ABOTD=ATOPD
      N=N+1
      IF(N.LT.NLP) GO TO 303
      RBNB=SRBALB(KLAM)
      RBNX=SRXALB(KLAM)
      RCNB=0.D0
      RCNX=0.D0
C                                  -------------------------------------
C                                  Get Oxygen UV absorption contribution
C                                  -------------------------------------
C----------------
      CALL GETO2A
C----------------
C                 ------------------------------------------------------
C                 Add Layers from Ground up. Retain Composite RBNB, RBNX
C             R,T,X of "A" (above) layer are corrected for O3 absorption
C             ----------------------------------------------------------

      DO 304 N=1,NL
      O2FHRL(N)=O2FHRL(N)/DKS0(15)*FULGAS(4)
      O2FHRB(N)=O2FHRB(N)/DKS0(15)*FULGAS(4)
      SRB(N)=RBNB
      SRX(N)=RBNX
      XANX=XNX(N)*(1.D0-AO3X(N)-O2FHRL(N))
      XANB=XNB(N)*(1.D0-AO3D(N)-O2FHRB(N))
      RASX=RNX(N)*(1.D0-AO3U(N))
      RASB=RNB(N)*(1.D0-AO3U(N))
      TANX=TNX(N)*(1.D0-AO3D(N))
      TANB=TNB(N)*(1.D0-AO3D(N))
      ABSRTX=1.D0-XANX-TANX-RASX
      ABSRTB=1.D0-XANB-TANB-RASB
      RARB=RASB*RBNB
      RARX=RASB*RBNX
      XATB=XANB+TANB
      DENOM=1.D0-RARB
      DB=(TANB+XANB*RARB)/DENOM
      DX=(TANX+XANX*RARX)/DENOM
      UB=RBNB*(XANB+DB)
      UX=RBNX*XANX+RBNB*DX
      RBNB=RASB+XATB*UB
      RBNX=RASX+XATB*UX
      XATC=XATB/(1.D0-RASB*RCNB)
      RCNX=RASX+(XANX*RCNX+TANX*RCNB)*XATC
      RCNB=RASB+RCNB*XATB*XATC
  304 CONTINUE
      VRD(NLP)=1.D0
      VRU(NLP)=RBNX
      SRKALB(15)=RBNX
      N=NL
      VRD(N)=XANX+DX
      VRU(N)=UX
  310 CONTINUE
      N=N-1
      XBNX=XNX(N)*(1.D0-AO3X(N)-O2FHRL(N))
      XBNB=XNB(N)*(1.D0-AO3D(N)-O2FHRB(N))
      RBNX=RNX(N)*(1.D0-AO3U(N))
      RBNB=RNB(N)*(1.D0-AO3U(N))
      TBNX=TNX(N)*(1.D0-AO3D(N))
      TBNB=TNB(N)*(1.D0-AO3D(N))

C     Add successively layer N (at bottom) to form upper composite layer
C     ------------------------------------------------------------------

      RARB=RASB*RBNB
      XBTB=XBNB+TBNB
      DENOM=1.D0/(1.D0-RARB)
      TANX=TBNX*XANX+XBTB*(TANX+XANX*RBNX*RASB)*DENOM
      RASB=RBNB+XBTB*XBTB*RASB*DENOM
      XANX=XANX*XBNX

C     Add upper & bottom composite layers to get flux at layer interface
C     ------------------------------------------------------------------

      RBNB=SRB(N)
      RBNX=SRX(N)
      DX=(TANX+XANX*RBNX*RASB)/(1.D0-RASB*RBNB)
      UX=RBNX*XANX+RBNB*DX
      VRD(N)=XANX+DX
      VRU(N)=UX
      IF(N.GT.1) GO TO 310
      DO 321 I=1,4
      SRKGAX(15,I)=DKS0X*XANX*(1.D0-PRNX(6,I))
      SRKGAD(15,I)=DKS0X*  DX*(1.D0-PRNB(6,I))
  321 CONTINUE

      DO 325 N=1,NLP
      VRD(N)=VRD(N)*DKS0X
      VRU(N)=VRU(N)*DKS0X
      SKDFLB(N,K)=VRD(N)
      SKUFLB(N,K)=VRU(N)
  325 CONTINUE
      SRIVIS=SRIVIS+VRD(NLP)
      SROVIS=SROVIS+VRU(NLP)
      PLAVIS=SROVIS/SRIVIS
      SRDVIS=SRDVIS+VRD(1)
      SRUVIS=SRUVIS+VRU(1)
      ALBVIS=SRUVIS/(SRDVIS+1.D-10)
      SRRVIS=SRRVIS+DKS0X*RCNX
      SRTVIS=SRTVIS+DKS0X*(TANX+XANX)
      SRXVIS=SRXVIS+DKS0X*XANX
      SRAVIS=1.D0-SRRVIS-SRTVIS

C     K16 strong absorbing contributions are computed without scattering
C     ------------------------------------------------------------------

      K=16
      DKS0X=DKS0(16)*S0COSZ
      SRKINC(16)=DKS0X
      SRKA16=0.D0
      DO 344 I=1,4
      SRKGAX(16,I)=0.D0
      SRKGAD(16,I)=0.D0
  344 CONTINUE
      DO 345 KK=1,3
      IF(KK.EQ.1) DKS0XX=DKS0X*0.002D0/0.011D0
      IF(KK.EQ.2) DKS0XX=DKS0X*0.008D0/0.011D0
      IF(KK.EQ.3) DKS0XX=DKS0X*0.001D0/0.011D0
      N=NLP
      TRNC=1.D0
      DO 330 M=1,NL
      N=N-1
      PLN=PL(N)
      CLX=DBLEXT(N,1)-DBLSCT(N,1)

C--------K=5-------CO2       DS0=.002
      IF(KK.EQ.1) THEN
      TRN1=0.D0
      ULN=ULGAS(N,2)*SECZ
      IF(ULN.GT.7.D0) ULN=7.D0
      TERMA=.003488*PLN*(1.+39.59*EXP(-8.769*ULN/(1.+4.419*ULN)))
     +     *(1.+ULN*(.001938*PLN-.00503*ULN))
      TERMB=(1.+.04712*PLN*(1.+.4877*ULN))
      TAUG=TERMA/TERMB*ULN
      TAU1=TAUG+CLX*SECZ
      IF(TAU1.LT.10.0) TRN1=EXP(-TAU1)
      FAC(N)=TRN1
      ENDIF

C--------K=7-------H2O       DS0=.008
      IF(KK.EQ.2) THEN
      TRN2=0.D0
      ULN=ULGAS(N,1)*SECZ
      TERMA=.001582*PLN*(1.+6.769*EXP(-9.59*ULN/(1.+5.026*ULN)))
     +     *(1.+ULN*(.2757E-03*PLN+.001429*ULN))
      TERMB=(1.+.003683*PLN*(1.+1.187*ULN))
      TAUG=TERMA/TERMB*ULN
      TAU2=TAUG+CLX*SECZ
      IF(TAU2.LT.10.0) TRN2=EXP(-TAU2)
      FAC(N)=TRN2
      ENDIF

C--------K=5-------O2        DS0=.001
      IF(KK.EQ.3) THEN
      TRN3=0.D0
      ULN=ULGAS(N,4)*SECZ
      TERMA=(.1366E-03-.2203E-07*TLN)*(1.+PLN*(.1497E-06*ULN+.001261))
      TERMB=(1.+.3867E-03*ULN)/(1.+.2075E-04*ULN)
      TAUG=TERMA/TERMB*ULN
      TAU3=TAUG+CLX*SECZ
      IF(TAU3.LT.10.0) TRN3=EXP(-TAU3)
      FAC(N)=TRN3
      ENDIF

      TRNC=TRNC*FAC(N)
      SRDFLB(N)=SRDFLB(N)+DKS0XX*TRNC
      SKDFLB(N,K)=SKDFLB(N,K)+DKS0XX*TRNC
  330 CONTINUE
      SRDFLB(NLP)=SRDFLB(NLP)+DKS0XX
      SRUFLB(1)=SRUFLB(1)+DKS0XX*TRNC*SRXALB(1)
      SKDFLB(NLP,K)=SKDFLB(NLP,K)+DKS0XX
      SKUFLB(1,K)=SKUFLB(1,K)+DKS0XX*TRNC*SRXALB(1)

C     For completeness, any incident flux at ground is relflected upward
C     ------------------------------------------------------------------

      TRNU=TRNC
      DO 340 N=2,NLP
      TRNU=TRNU*FAC(N-1)
      SRUFLB(N)=SRUFLB(N)+DKS0XX*TRNC*SRXALB(1)*TRNU
      SKUFLB(N,K)=SKUFLB(N,K)+DKS0XX*TRNC*SRXALB(1)*TRNU
  340 CONTINUE
      DO 341 I=1,4
      SRKGAX(16,I)=SRKGAX(16,I)+DKS0XX*TRNC*(1.D0-PRNX(1,I))
  341 CONTINUE
      SRKA16=SRKA16+TRNU*SRXALB(1)

      SRINIR=SRINIR+DKS0XX
      SRONIR=SRONIR+DKS0XX*TRNU*SRXALB(1)
      SRDNIR=SRDNIR+SKDFLB(1,K)
      SRUNIR=SRUNIR+SKUFLB(1,K)
  345 CONTINUE
      PLANIR=SRONIR/SRINIR
      ALBNIR=SRUNIR/(SRDNIR+1.D-10)
      SRKALB(16)=SRKA16/DKS0X

      DO 350 N=1,NLP
      SRDFLB(N)=SRDFLB(N)+VRD(N)
      SRUFLB(N)=SRUFLB(N)+VRU(N)
      SRNFLB(N)=SRDFLB(N)-SRUFLB(N)
  350 CONTINUE
      DO 360 N=1,NL
      SRFHRL(N)=SRNFLB(N+1)-SRNFLB(N)
  360 CONTINUE
      SRRNIR=SRRNIR+DKS0X*RCNX
      SRTNIR=SRTNIR+DKS0X*(TANX+XANX)
      SRXNIR=SRXNIR+DKS0X*XANX

      S0VIS=0.53D0*S0
      SRTVIS=SRTVIS/S0VIS
      SRRVIS=SRRVIS/S0VIS
      SRXVIS=SRXVIS/S0VIS
      SRAVIS=1.D0-SRTVIS-SRRVIS

      S0NIR=0.47D0*S0
      SRTNIR=SRTNIR/S0NIR
      SRRNIR=SRRNIR/S0NIR
      SRXNIR=SRXNIR/S0NIR
      SRANIR=1.D0-SRTNIR-SRRNIR


C     ------------------------------------------------------------------
C     FSRNFG defines the total solar flux absorbed at the ground surface
C              taking into account the albedo of different surface types
C     Thus:
C              SRNFLB(1)=POCEAN*FSRNFG(1)+PEARTH*FSRNFG(2)
C                       + POICE*FSRNFG(3)+ PLICE*FSRNFG(4)
C
C     NOTE:    If any surface type POCEAN, PEARTH, POICE, PLICE are Zero
C              the corresponding FSRNFG(I) absorbed solar flux at ground
C              is computed with that surface-type albedo set equal to 0.
C              ---------------------------------------------------------

      DO 420 I=1,4
      SUMX=0.D0
      SUMD=0.D0
      DO 410 K=1,16
      SUMX=SUMX+SRKGAX(K,I)
      SUMD=SUMD+SRKGAD(K,I)
  410 CONTINUE
      FSRNFG(I)=SUMX+SUMD
  420 CONTINUE


      DO 520 L=1,NLP
      DO 510 K=1,16
      SKNFLB(L,K)=SKDFLB(L,K)-SKUFLB(L,K)
  510 CONTINUE
  520 CONTINUE

      DO 540 L=1,NL
      DO 530 K=1,16
      SKFHRL(L,K)=SKNFLB(L+1,K)-SKNFLB(L,K)
  530 CONTINUE
  540 CONTINUE

      DO 560 L=1,NLP
      SUMD=0.D0
      SUMU=0.D0
      SUMN=0.D0
      SUMH=0.D0
      DO 550 K=1,16
      SUMD=SUMD+SKDFLB(L,K)
      SUMU=SUMU+SKUFLB(L,K)
      SUMN=SUMN+SKNFLB(L,K)
      SUMH=SUMH+SKFHRL(L,K)
  550 CONTINUE
      SKDFLB(L,17)=SUMD
      SKUFLB(L,17)=SUMU
      SKNFLB(L,17)=SUMN
      SKFHRL(L,17)=SUMH
  560 CONTINUE

      RETURN
      END SUBROUTINE SOLARM

      SUBROUTINE SETGTS
CCC   SUBROUTINE GTSALB(GIN,TAUIN,RBBOUT,RBBIN,EGIN,TAUOUT,KGTAUR)

      !INCLUDE 'BR00B.COM'

      DIMENSION FFKG(4,3),RBBK(3)
      DIMENSION GVALUE(14)
      DATA GVALUE/.0,.25,.45,.50,.55,.60,.65,.70,.75,.80,.85,.90,.95,1./
      SAVE GVALUE
CCC   DATA IFIRST/1/

CCC   IF(IFIRST.EQ.1) THEN
      DO 100 I=1,122
      TAUTGD(I)=(I-1)*0.1D0
      IF(I.GT.24) TAUTGD(I)=(I-24)*0.2D0+2.2D0
      IF(I.GT.48) TAUTGD(I)=(I-48)*0.5D0+7.0D0
      IF(I.GT.72) TAUTGD(I)=(I-72)+19.0D0
      IF(I.GT.96) TAUTGD(I)=(I-96)*5.0D0+40.0D0
      IF(I.GT.112) TAUTGD(I)=(I-112)*100.0D0+100.0D0
      IF(I.EQ.121) TAUTGD(I)=9999.99D0
      IF(I.EQ.122) TAUTGD(I)=12000.0D0
  100 CONTINUE

      DO 110 I=1,768
      IF(I.LT.602) TAUTGS(I)=(I-1)*0.05D0
      IF(I.GT.601) TAUTGS(I)=(I-601)*0.50D0+30.0D0
      IF(I.GT.741) TAUTGS(I)=(I-741)*50.0D0+100.D0
      IF(I.GT.758) TAUTGS(I)=(I-758)*1000.D0
  110 CONTINUE

      DO 130 J=1,13
      DO 120 I=1,768
      CWM=0.5
      CWE=0.5
      IF(I.GT.759) CWM=0.0
      IF(I.GT.759) CWE=0.0
      TIJ=TAUTGS(I)
      CALL SPLINE(TAUTGD,TGDATA(1,J),122,TIJ,RBBI,CWM,CWE,0)
      SALBTG(I,J)=RBBI
  120 CONTINUE
  130 CONTINUE
      DO 150 J=1,13
      DO 140 I=2,1000
      RBB=(I-1)*0.001D0
      CWM=0.5
      CWE=0.5
      CALL SPLINE(SALBTG(1,J),TAUTGS,768,RBB,BTAU,CWM,CWE,0)
      TAUGSA(I,J)=BTAU
  140 CONTINUE
  150 CONTINUE
      DO 160 J=1,14
      SALBTG(1,J)=0.D0
      TAUGSA(1,J)=0.D0
      TAUGSA(1001,J)=10000.D0
  160 CONTINUE

      DO 170 I=1,768
      SALBTG(I,14)=SALBTG(I,13)*2.D0-SALBTG(I,12)
  170 CONTINUE
      DO 180 I=1,1001
CCC   SALBI=(I-1)*0.001D0
      TAUGSA(I,14)=TAUGSA(I,13)*2.D0-TAUGSA(I,12)
  180 CONTINUE
CCC   IFIRST=0
CCC   ENDIF
      RETURN

      ENTRY GTSALB(GIN,TAUIN,RBBOUT,RBBIN,EGIN,TAUOUT,KGTAUR)

      KTERPL=0
      CUSPWM=0.5
      CUSPWE=0.5

      G=GIN
      TAU=TAUIN
      RBB=RBBIN
      EG=EGIN

      RBBOUT=0.0
      TAUOUT=0.0
C                                           ---------------------------
C                                           OPTICAL DEPTH INTERPOLATION
C                                           0.05 ON (0.00 < TAU < 30.0)
C                                           0.50 ON (30.0 < TAU < 100.)
C                                           50.0 ON (100. < TAU < 1000)
C                                           ---------------------------

      IF(KGTAUR.EQ.2) GO TO 300

  200 CONTINUE

      DELTAU=0.05D0
      TI=TAU/DELTAU
      IT=TI
      IF(IT.GT.599) GO TO 210
      WTJ=TI-IT
      WTI=1.D0-WTJ
      IT=IT+1
      GO TO 240
  210 CONTINUE
      DELTAU=0.50D0
      TI=TAU/DELTAU
      IT=TI
      IF(IT.GT.199) GO TO 220
      WTJ=TI-IT
      WTI=1.0-WTJ
      IT=IT+541
      GO TO 240
  220 CONTINUE
      DELTAU=50.0D0
      TI=TAU/DELTAU
      IT=TI
      IF(IT.GT.19) GO TO 230
      WTJ=TI-IT
      WTI=1.0-WTJ
      IT=IT+649
      GO TO 240
  230 CONTINUE
      DELTAU=1000.0D0
      TI=TAU/DELTAU
      IT=TI
      WTJ=TI-IT
      WTI=1.0-WTJ
      IT=IT+758
  240 CONTINUE
      JT=IT+1

C                                    ---------------------------------
C                                    ASYMMETRY PARAMETER INTERPOLATION
C                                    0.05 CUBIC SPLINE (0.5 < G < 0.9)
C                                    0.25 QUADRATIC ON (0.0 < G < 0.5)
C                                    LINEAR EXTRAP FOR (.95 < G < 1.0)
C                                    ---------------------------------

      GI=G*20.D0
      IF(GI.GT.10.0) GO TO 250
      IG=2
      JG=3
      ITERPL=1
      GO TO 260

  250 CONTINUE
      ITERPL=4
      IG=GI
      WGJ=GI-IG
      WGI=1.D0-WGJ
      IG=IG-6
      IF(IG.GT.12) THEN
      ITERPL=2
      IG=12
      ENDIF
      JG=IG+1

  260 CONTINUE

      IGM=IG-1
      JGP=JG+1

      K=0
      DO 270 KG=IGM,JGP
      K=K+1
      F1=SALBTG(IT-1,KG)
      F2=SALBTG(IT  ,KG)
      F3=SALBTG(JT  ,KG)
      F4=SALBTG(JT+1,KG)
      IF(IT.EQ.1) F1=-F3
      F21=(F2-F1)
      F32=(F3-F2)
      F43=(F4-F3)
      F3221=(F32+F21)*0.5D0
      F4332=(F43+F32)*0.5D0
      A=F2
      B=F3221
      C=3.D0*F32-F3221-F3221-F4332
      D=(F3221+F4332-F32-F32)
      XF=WTJ
      FFCUSP=A+XF*(B+XF*(C+XF*D))
      XE=1.D0-XF-XF
      IF(XE.LT.0.0) XE=-XE
      XEXM=XE**2
      CUSPWT=(1.0-XEXM)*CUSPWM+XEXM*CUSPWE
      FFLINR=A+XF*F32
      FFKG(K,1)=FFCUSP*CUSPWT+FFLINR*(1.D0-CUSPWT)
      FFKG(K,2)=F2
      FFKG(K,3)=F3
  270 CONTINUE

      IF(ITERPL.LT.4) GO TO 290

      DO 280 K=1,3
      F1=FFKG(1,K)
      F2=FFKG(2,K)
      F3=FFKG(3,K)
      F4=FFKG(4,K)
      F21=(F2-F1)
      F32=(F3-F2)
      F43=(F4-F3)
      F3221=(F32+F21)*0.5D0
      F4332=(F43+F32)*0.5D0
      A=F2
      B=F3221
      C=3.D0*F32-F3221-F3221-F4332
      D=(F3221+F4332-F32-F32)
      XF=WGJ
      FFCUSP=A+XF*(B+XF*(C+XF*D))
      XE=1.D0-XF-XF
      IF(XE.LT.0.0) XE=-XE
      XEXM=XE**2
      CUSPWT=(1.0-XEXM)*CUSPWM+XEXM*CUSPWE
      FFLINR=A+XF*F32
      RBBK(K)=FFCUSP*CUSPWT+FFLINR*(1.D0-CUSPWT)
  280 CONTINUE
      RBB=RBBK(1)
      RB2=RBBK(2)
      RB3=RBBK(3)
      TBB=TAU
      TB2=TAUTGS(IT)
      TB3=TAUTGS(JT)
      IF(KGTAUR.EQ.1) RETURN
      IF(KTERPL.EQ.1) GO TO 400
      GO TO 300

  290 CONTINUE
      XG=G*2.D0-0.5D0
      IF(ITERPL.EQ.2) XG=G*10.D0-9.D0
      XM=1.D0-XG-XG
      XP=1.D0+XG+XG
      RBB=XM*XP*FFKG(ITERPL+1,1)-XG*XM*FFKG(ITERPL,1)+XG*XP*FFKG(4,1)
      RB2=XM*XP*FFKG(ITERPL+1,2)-XG*XM*FFKG(ITERPL,2)+XG*XP*FFKG(4,2)
      RB3=XM*XP*FFKG(ITERPL+1,3)-XG*XM*FFKG(ITERPL,3)+XG*XP*FFKG(4,3)

      IF(KGTAUR.EQ.1) RETURN
      IF(KTERPL.EQ.1) GO TO 400

  300 CONTINUE
      RBBB=RBB

      RI=RBB*1000.D0
      IR=RI
      WRJ=RI-IR
      WRI=1.D0-WRJ
      IR=IR+1
      JR=IR+1

      EI=EG*20.D0
      IF(EI.GT.10.0) GO TO 310
      IE=2
      JE=3
      ITERPL=1
      GO TO 320
  310 CONTINUE

      ITERPL=4
      IE=EI
      WEJ=EI-IE
      WEI=1.D0-WEJ
      IE=IE-6
      IF(IE.GT.12) THEN
      ITERPL=2
      IE=12
      ENDIF
      JE=IE+1
  320 CONTINUE

      DELALB=0.001D0
      IEM=IE-1
      JEP=JE+1
      K=0
      DO 330 KG=IEM,JEP
      K=K+1
      F1=TAUGSA(IR-1,KG)
      F2=TAUGSA(IR  ,KG)
      F3=TAUGSA(JR  ,KG)
      F4=TAUGSA(JR+1,KG)
      IF(IR.EQ.1) F1=-F3
      F21=(F2-F1)
      F32=(F3-F2)
      F43=(F4-F3)
      F3221=(F32+F21)*0.5D0
      F4332=(F43+F32)*0.5D0
      A=F2
      B=F3221
      C=3.D0*F32-F3221-F3221-F4332
      D=(F3221+F4332-F32-F32)
      XF=WRJ
      FFCUSP=A+XF*(B+XF*(C+XF*D))
      XE=1.D0-XF-XF
      IF(XE.LT.0.0) XE=-XE
      XEXM=XE**2
      CUSPWT=(1.0-XEXM)*CUSPWM+XEXM*CUSPWE
      FFLINR=A+XF*F32
      FFKG(K,1)=FFCUSP*CUSPWT+FFLINR*(1.D0-CUSPWT)
      FFKG(K,2)=F2
      FFKG(K,3)=F3
  330 CONTINUE
      X1=GVALUE(IE-1)
      X2=GVALUE(IE  )
      X3=GVALUE(JE  )
      X4=GVALUE(JE+1)
      XX=WEJ
      IF(ITERPL.LT.4) GO TO 350

      DO 340 K=1,3
      F1=FFKG(1,K)
      F2=FFKG(2,K)
      F3=FFKG(3,K)
      F4=FFKG(4,K)
      F21=(F2-F1)
      F32=(F3-F2)
      F43=(F4-F3)
      F3221=(F32+F21)*0.5D0
      F4332=(F43+F32)*0.5D0
      A=F2
      B=F3221
      C=3.D0*F32-F3221-F3221-F4332
      D=(F3221+F4332-F32-F32)
      XF=WEJ
      FFCUSP=A+XF*(B+XF*(C+XF*D))
      XE=1.D0-XF-XF
      IF(XE.LT.0.0) XE=-XE
      XEXM=XE**2
      CUSPWT=(1.0-XEXM)*CUSPWM+XEXM*CUSPWE
      FFLINR=A+XF*F32
      RBBK(K)=FFCUSP*CUSPWT+FFLINR*(1.D0-CUSPWT)
  340 CONTINUE
      TBB=RBBK(1)
      TB2=RBBK(2)
      TB3=RBBK(3)

      IF(KTERPL.EQ.1) GO TO 400
      GO TO 390

  350 CONTINUE
      XG=EG*2.D0-0.5D0
      IF(ITERPL.EQ.2) XG=G*10.D0-9.D0
      XM=1.D0-XG-XG
      XP=1.D0+XG+XG
      TBB=XM*XP*FFKG(ITERPL+1,1)-XG*XM*FFKG(ITERPL,1)+XG*XP*FFKG(4,1)
      TB2=XM*XP*FFKG(ITERPL+1,2)-XG*XM*FFKG(ITERPL,2)+XG*XP*FFKG(4,2)
      TB3=XM*XP*FFKG(ITERPL+1,3)-XG*XM*FFKG(ITERPL,3)+XG*XP*FFKG(4,3)
      IF(KTERPL.EQ.1) GO TO 400
  390 CONTINUE
      KTERPL=1
      TAU=TBB
      G=EGIN
      GO TO 200
  400 CONTINUE
      IF(ABS(WTI*WTJ).LT.0.1D0) DTAU=(RBBB-RB2)/(RB3-RB2)
      IF(ABS(WTI*WTJ).GE.0.1D0) THEN
      C=(RB3-RBB)/WTI-(RBB-RB2)/WTJ
      B=(RBB-RB2)/WTJ-WTJ*C
      A=RB2
      BB=B*B+4.D0*C*(RBBB-A)
      IF(BB.GT.0.D0) DTAU=(SQRT(BB)-B)/(C+C)
      ENDIF
      TAUOUT=(IT-1+DTAU)*DELTAU
      RBBOUT=RBBB

      RETURN
CCC   END SUBROUTINE GTSALB
      END SUBROUTINE SETGTS

      SUBROUTINE SGPGXG(XMU,TAU,G,GG)

      !INCLUDE 'BR00B.COM'

C     ----------------------------------------------------------------
C     COSBAR ADJUSTMENT TO REPRODUCE THE SOLAR ZENITH ANGLE DEPENDENCE
C     FOR AEROSOL ALBEDO FOR OPTICAL THICKNESSES [0.0 < TAU < 10000.0]
C     ----------------------------------------------------------------


C                          -------------------------------------------
C                          XMU (COSZ) SOLAR ZENITH ANGLE INTERPOLATION
C                          DATA INTERVAL:  0.02  ON  [0.0 < XMU < 1.0]
C                          -------------------------------------------

      XI=XMU*50.D0+0.999999D0  ! safer: +1.  ?????
      IX=XI
      JX=IX+1
      WXJ=XI-IX
      WXI=1.D0-WXJ

C                                      -------------------------------
C                                      COSBAR DEPENDENCE INTERPOLATION
C                                         0.10 ON [0.0 < COSBAR < 1.0]
C                                      -------------------------------

      GI=G*10.D0
      IG=GI
      WGJ=GI-IG
      WGI=1.D0-WGJ
      IG=IG+1
      JG=IG+1

C                            -----------------------------------------
C                               AEROSOL TAU INTERPOLATION INTERVALS
C                            -----------------------------------------
C                            dTau      1      1  (Lin Int)  61     62
C                            0.10 ON [0.00 , 0.00 < TAU < 6.00 , 6.10]
C                                     63     64            92     93
C                            0.50 ON [5.50 , 6.00 < TAU < 20.0 , 20.5]
C                                     94     95           111    112
C                            5.00 ON [15.0 , 20.0 < TAU < 100. , 105.]
C                                     113    114          132    133
C                            50.0 ON [50.0 , 100. < TAU < 1000 , 1050]
C                                            134          143
C                            1000 ON [     , 1000 < TAU < 10000,     ]
C                            -----------------------------------------

      IF(TAU.LT.6.D0) THEN
      TI=TAU*10.D0+1.
      IT=TI
      WTJ=TI-IT
      IT0=0
CCC   GO TO 100
      ELSEIF(TAU.LT.20.D0) THEN
      TI=(TAU-6.D0)*2.00D0+2.D0
      IT=TI
      WTJ=TI-IT
      IT0=62
CCC   GO TO 100
      ELSEIF(TAU.LT.100.D0) THEN
      TI=(TAU-20.D0)*0.20D0+2.D0
      IT=TI
      WTJ=TI-IT
      IT0=93
CCC   GO TO 100
      ELSEIF(TAU.LT.1000.D0) THEN
      TI=(TAU-100.D0)*0.02D0+2.D0
      IT=TI
      WTJ=TI-IT
      IT0=112
CCC   GO TO 100
      ELSE
      TI=TAU*0.001D0+1.D-6
      IT=TI
      WTJ=TI-IT
      IF(IT.GT.9) IT=9
      IT0=133
      ENDIF
C 100 CONTINUE
      WTI=1.D0-WTJ
      IT=IT+IT0
      JT=IT+1

      GG=WGI*(WTI*(WXI*GTAU(IX,IG,IT)+WXJ*GTAU(JX,IG,IT))
     +      + WTJ*(WXI*GTAU(IX,IG,JT)+WXJ*GTAU(JX,IG,JT)))
     +  +WGJ*(WTI*(WXI*GTAU(IX,JG,IT)+WXJ*GTAU(JX,JG,IT))
     +      + WTJ*(WXI*GTAU(IX,JG,JT)+WXJ*GTAU(JX,JG,JT)))

      RETURN
      END SUBROUTINE SGPGXG

      SUBROUTINE BCTAUW(XJYEAR,IDEC,JDEC,BCWTID,BCWTJD)
      IMPLICIT REAL*8 (A-H,O-Z)

C-------------------------------------------------------------------
C     Black Carbon interdecadal TAU interpolation is based on linear
C     TAU trend (between decadal global TAUmaps) with a superimposed
C     intra-decadal time dependence scaled to the Black Carbon Total
C     emission rate.
C
C        INPUT:  XJYEAR  (Fractional Julian year)
C                                                 by 25     by 10
C       OUTPUT:  IDEC    (Map Index: I= -3,4-7  (0) -1925,1950-1980)
C                JDEC    (Map Index: J=1-4,5-8  1875-1950,1960-1990)
C
C              BCWTID    (Multiplicative Weight for BC TAU-Map IDEC)
C              BCWTJD    (Multiplicative Weight for BC TAU-Map JDEC)
C
C-------------------------------------------------------------------

      DIMENSION EYEAR(50),BCEHC(50),BCEBC(50),BCEDI(50),BCTOT(50)

      DIMENSION BCE(5,45)

C    Global Annual Emissions of BC U   Emission (Mt/yr)

      DATA BCE/
C      Year    Hard_Coal    Brown_Coal      Diesel        Total
     A 50.0, 2.280581713, 0.4449132979, 0.1599090248, 2.885536671,
     B 51.0, 2.443193913, 0.4855868816, 0.1884280443, 3.117194653,
     C 52.0, 2.473641872, 0.5115299225, 0.2027695477, 3.187930107,
     D 53.0, 2.481340885, 0.5448409319, 0.2149295360, 3.241089582,
     E 54.0, 2.505670071, 0.5780177116, 0.2343477309, 3.317960978,
     F 55.0, 2.698692560, 0.6238067150, 0.2733324766, 3.595800638,
     G 56.0, 2.855226278, 0.6531309485, 0.3043369055, 3.812692404,
     H 57.0, 2.975781679, 0.6821750998, 0.3207367063, 3.978575468,
     I 58.0, 3.341105223, 0.7035279870, 0.3370627165, 4.381746292,
     J 59.0, 3.638528824, 0.7075053453, 0.3695519567, 4.715488434,
     K 60.0, 3.770926714, 0.7416650057, 0.3832504749, 4.896034241,
     L 61.0, 3.392980337, 0.7805693150, 0.4217525721, 4.595387459,
     M 62.0, 3.288835049, 0.8179932237, 0.4603823125, 4.567360401,
     N 63.0, 3.359177589, 0.8604368567, 0.5090782642, 4.728550911,
     O 64.0, 3.432664871, 0.8952696323, 0.5388473868, 4.866865158,
     P 65.0, 3.529418945, 0.8819132447, 0.5785927773, 4.989773750,
     Q 66.0, 3.577459812, 0.8817394972, 0.6323299408, 5.091631413,
     R 67.0, 3.418204546, 0.8635972142, 0.6592246890, 4.941041946,
     S 68.0, 3.452457905, 0.8943673372, 0.7338049412, 5.080585003,

ceq   DATA BCE2/
     A 69.0, 3.626069546, 0.9298774004, 0.7889106274, 5.344810009,
     B 70.0, 3.264039755, 0.9229136109, 0.8880128860, 5.074741840,
     C 71.0, 3.437611580, 0.9374827743, 0.9531223178, 5.328329086,
     D 72.0, 3.473345757, 0.7836616039, 1.0180075170, 5.274850368,
     E 54.0, 2.505670071, 0.5780177116, 0.2343477309, 3.317960978,
     F 74.0, 3.506143808, 0.8251076341, 1.0828053950, 5.413989067,
     G 75.0, 3.906814098, 0.8527192473, 1.0454736950, 5.804963112,
     H 76.0, 4.005736828, 0.8900613785, 1.1400985720, 6.035901546,
     I 77.0, 4.236912251, 0.9103702307, 1.2190728190, 6.366260529,
     J 78.0, 4.459666252, 0.9303293228, 1.2408012150, 6.630728722,
     K 79.0, 4.697422504, 0.9856286645, 1.3019220830, 6.984815121,
     L 80.0, 4.796229839, 0.9959300756, 1.2336660620, 7.026207924,
     M 81.0, 4.789204121, 1.0459070210, 1.1664049630, 7.001126766,
     N 82.0, 4.872739315, 1.0975246430, 1.1601715090, 7.130136490,
     O 83.0, 4.983223438, 1.1424025300, 1.1732926370, 7.298912525,
     P 84.0, 5.265352249, 1.2178678510, 1.2251536850, 7.708741188,
     Q 85.0, 5.763637543, 1.2965050940, 1.2428865430, 8.303324699,
     R 86.0, 5.924767494, 1.3386499880, 1.2930148840, 8.556744576,
     S 87.0, 6.155550480, 1.3738890890, 1.3162037130, 8.845513344,

ceq   DATA BCE3/
     A 88.0, 6.379704475, 1.3670797350, 1.3813229800, 9.127896309,
     B 89.0, 6.594299316, 1.4169263840, 1.4029121400, 9.414231300,
     C 90.0, 6.566919804, 1.4685817960, 1.4224120380, 9.458042145,
     D 91.0, 6.661097050, 1.2067918780, 1.4163945910, 9.284657478,
     E 92.0, 7.737902641, 1.3509917260, 1.4471185210, 10.53625107,
     F 93.0, 7.393332005, 1.2448183300, 1.4543261530, 10.09271908,
     G 94.0, 7.515841007, 1.2333894970, 1.4780857560, 10.22745800/

C     DATA POST90/50.D0/
C     DATA PRE50/50.D0/
      PARAMETER(POST90=-250.D0)
      DATA IFIRST/1/

      SAVE IFIRST,EYEAR,BCEHC,BCEBC,BCEDI,BCTOT,BCE

      IF(IFIRST.EQ.1) THEN
      DO 110 I=1,45
      EYEAR(I)=BCE(1,I)+1900.D0
      BCEHC(I)=BCE(2,I)
      BCEBC(I)=BCE(3,I)
      BCEDI(I)=BCE(4,I)
      BCTOT(I)=BCE(5,I)
  110 CONTINUE
      IFIRST=0
      ENDIF

      IF(XJYEAR.LT.1950.D0) THEN
C**** 0 until 1850, then lin. interpolate obs. data (every 25 years)
        XDEC=(XJYEAR-1850.d0)/25.d0
        IF (XDEC.lt.0.) Xdec=0.
        IDEC=XDEC
        DDEC=XDEC-IDEC
        BCWTID=1.-DDEC
        JDEC=IDEC+1
        BCWTJD=DDEC
        IF(IDEC.LT.1) THEN
          IDEC=1
          JDEC=1
          BCWTID=XDEC
          BCWTJD=0.
        END IF
      ELSE IF(XJYEAR.GE.1990.D0) THEN
C**** Slow reduction after 1990 (POST90=-250 years e-folding time)
C**** Actually we will use no reduction after 1990
        DYEAR=XJYEAR-1990.D0
        BCWTID=0.D0
        BCWTJD=1.D0 !  BCWTJD=1.D0*EXP(DYEAR/POST90)
        IDEC=8
        JDEC=8
      ELSE
C**** lin. interpolate obs. data (every 10 years) 1950-1990
        DYEAR=XJYEAR-1950.D0
        IYEAR=DYEAR
        DELTAY=DYEAR-IYEAR
        DELDEC=DYEAR/10.D0-IYEAR/10
        IDEC=(IYEAR+10)/10
        JDEC=IDEC+1
        IYDI=1+(IDEC-1)*10
        IYDJ=IYDI+10
        IYYI=IYEAR+1
        IYYJ=IYYI+1
        BCED=BCTOT(IYDI)+DELDEC*(BCTOT(IYDJ)-BCTOT(IYDI))
        BCEY=BCTOT(IYYI)+DELTAY*(BCTOT(IYYJ)-BCTOT(IYYI))
        RATYD=BCEY/BCED
        BCWTID=RATYD*(1.D0-DELDEC)
        BCWTJD=RATYD*DELDEC
        IDEC=IDEC+3
        JDEC=JDEC+3
      END IF

      RETURN
      END SUBROUTINE BCTAUW

      SUBROUTINE SUTAUW(XJYEAR,IDEC,JDEC,SUWTID,SUWTJD)
      IMPLICIT REAL*8 (A-H,O-Z)

C-------------------------------------------------------------------
C     Anthropogenic Sulfate inter-decadal TAU interpolation is based
C     on a linear TAU trend (between decadal global TAU-maps) with a
C     superimposed intradecadal time dependence scaled in proportion
C     to the Anthropogenic Sulfate global emission rate.
C
C        INPUT:  XJYEAR  (Fractional Julian year)
C                                                 by 25     by 10
C       OUTPUT:  IDEC    (Map Index: I= -3,4-7  (0) -1925,1950-1980)
C                JDEC    (Map Index: J=1-4,5-8  1875-1950,1960-1990)
C
C              SUWTID    (Multiplicative Weight for SU TAU-Map IDEC)
C              SUWTJD    (Multiplicative Weight for SU TAU-Map JDEC)
C
C-------------------------------------------------------------------

      DIMENSION EYEAR(50),SUANT(50),SUNAT(50)

      DIMENSION SUE(3,41)

C     Global Emission of Sulfate

C     Emission (Mt/yr)
C               year      Anthropogenic_Sulfate Natural_Sulfate
      DATA SUE/
     A          1950.0,     30.46669769,           14.4,
     B          1951.0,     32.38347244,           14.4,
     C          1952.0,     32.18632889,           14.4,
     D          1953.0,     32.83379745,           14.4,
     E          1954.0,     32.79270935,           14.4,
     F          1955.0,     35.79611969,           14.4,
     G          1956.0,     39.93603897,           14.4,
     H          1957.0,     38.68806839,           14.4,
     I          1958.0,     39.35904312,           14.4,
     J          1959.0,     41.06065369,           14.4,
     K          1960.0,     42.67050934,           14.4,
     L          1961.0,     41.32410431,           14.4,
     M          1962.0,     41.80470276,           14.4,
     N          1963.0,     43.26312637,           14.4,
     O          1964.0,     44.68368530,           14.4,
     P          1965.0,     45.81701660,           14.4,
     Q          1966.0,     46.61584091,           14.4,
     R          1967.0,     46.42276001,           14.4,
     S          1968.0,     47.77438354,           14.4,
ceq   DATA SUE2/
     A          1969.0,     49.30817032,           14.4,
     B          1970.0,     52.81050873,           14.4,
     C          1971.0,     52.95043945,           14.4,
     D          1972.0,     54.10167694,           14.4,
     E          1973.0,     55.93037415,           14.4,
     F          1974.0,     57.31056213,           14.4,
     G          1975.0,     58.52788162,           14.4,
     H          1976.0,     59.71361542,           14.4,
     I          1977.0,     62.59599304,           14.4,
     J          1978.0,     61.98198318,           14.4,
     K          1979.0,     64.71042633,           14.4,
     L          1980.0,     65.28986359,           14.4,
     M          1981.0,     63.23768234,           14.4,
     N          1982.0,     62.88000488,           14.4,
     O          1983.0,     61.45023346,           14.4,
     P          1984.0,     63.85008621,           14.4,
     Q          1985.0,     66.47412872,           14.4,
     R          1986.0,     68.00902557,           14.4,
     S          1987.0,     69.87956238,           14.4,
ceq   DATA SUE3/
     A          1988.0,     70.52937317,           14.4,
     B          1989.0,     72.06355286,           14.4,
     C          1990.0,     71.29174805,           14.4/

C     DATA POST90/100.D0/
C     DATA PRE50/50.D0/
      PARAMETER(POST90=-250.D0)
      DATA IFIRST/1/

      SAVE  IFIRST,EYEAR,SUANT,SUNAT,SUE

      IF(IFIRST.EQ.1) THEN
      DO 110 I=1,41
      EYEAR(I)=SUE(1,I)
      SUANT(I)=SUE(2,I)
      SUNAT(I)=SUE(3,I)
  110 CONTINUE
      IFIRST=0
      ENDIF

      IF(XJYEAR.LT.1950.D0) THEN
C**** 0 until 1850, then lin. interpolate obs. data (every 25 years)
        XDEC=(XJYEAR-1850.d0)/25.d0
        IF (XDEC.lt.0.) Xdec=0.
        IDEC=XDEC
        DDEC=XDEC-IDEC
        SUWTID=1.-DDEC
        JDEC=IDEC+1
        SUWTJD=DDEC
        IF(IDEC.LT.1) THEN
          IDEC=1
          JDEC=1
          SUWTID=XDEC
          SUWTJD=0.
        END IF
      ELSE IF(XJYEAR.GE.1990.D0) THEN
C**** Slow reduction after 1990 (POST90=-250 years e-folding time)
C**** Actually we will use no reduction after 1990
        DYEAR=XJYEAR-1990.D0
        SUWTID=0.D0
        SUWTJD=1.D0  !  SUWTJD=1.D0*EXP(DYEAR/POST90)
        IDEC=8
        JDEC=8
      ELSE
C**** lin. interpolate obs. data (every 10 years) 1950-1990
        DYEAR=XJYEAR-1950.D0
        IYEAR=DYEAR
        DELTAY=DYEAR-IYEAR
        DELDEC=DYEAR/10.D0-IYEAR/10
        IDEC=(IYEAR+10)/10
        JDEC=IDEC+1
        IYDI=1+(IDEC-1)*10
        IYDJ=IYDI+10
        IYYI=IYEAR+1
        IYYJ=IYYI+1
        SUED=SUANT(IYDI)+DELDEC*(SUANT(IYDJ)-SUANT(IYDI))
        SUEY=SUANT(IYYI)+DELTAY*(SUANT(IYYJ)-SUANT(IYYI))
        RATYD=SUEY/SUED
        SUWTID=RATYD*(1.D0-DELDEC)
        SUWTJD=RATYD*DELDEC
        IDEC=IDEC+3
        JDEC=JDEC+3
      END IF

      RETURN
      END SUBROUTINE SUTAUW

!!    SUBROUTINE STDAER(QXAER,QSAER,QCAER,ATAER,PIAMAX,NA)
!!
!!    !INCLUDE 'BR00B.COM'
!!
C     ------------------------------------------------------------------
C     STDAER   Selects Mie Scattering Parameters for Aerosol Climatology
C                     and applies PIAMAX, FSAERA, FTAERA Scaling Factors
C                     --------------------------------------------------
!!
!!    DIMENSION QXAER(6),QSAER(6),QCAER(6),ATAER(33)
!!    DIMENSION FSAERA(10),FTAERA(10),PIAMAX(10)
!!    FSAERN=1.D0                ! FSAERA,FTAERA Scaling moved to GETAER
!!    FTAERN=1.D0
!!
!!    DO 110 I=1,6
!!    PI0F=1.D0
!!    PI0=SRBQSC(I,NA)/SRBQEX(I,NA)
!!    IF(PI0.GT.PIAMAX(NA)) PI0F=PIAMAX(NA)/PI0
!!    QXAER(I)=SRBQEX(I,NA)*FSAERN
!!    QSAER(I)=SRBQSC(I,NA)*FSAERN*PI0F
!!    QCAER(I)=SRBQCB(I,NA)
!!110 CONTINUE
!!
!!    DO 120 I=1,33
!!    ATAER(I)=(TRBQEX(I,NA)-TRBQSC(I,NA))*FTAERN
!!120 CONTINUE
!!
!!    RETURN
!!    END SUBROUTINE STDAER
!!
      SUBROUTINE SDDUST(QQDUST,SSDUST,GGDUST,Q55DST,SIZDST,FNDUST)

      !INCLUDE 'BR00B.COM'

C     ------------------------------------------------------------------
C     SDDUST     Selects Solar Mie Scattering Parameters for Desert Dust
C                                     and applies FSDUST Scaling Factors
C                                     ----------------------------------

      DIMENSION QQDUST(6),SSDUST(6),GGDUST(6)

      RDI=0.0
      DO 110 J=1,25
      RDJ=REFD25(J)
      IF(SIZDST.LT.RDJ) GO TO 120
      RDI=RDJ
  110 CONTINUE
      WTJ=1.0
      J=25
      I=25
      GO TO 130
  120 WTJ=(SIZDST-RDI)/(RDJ-RDI)
      I=J-1
      IF(I.LT.1) I=1
  130 WTI=1.0-WTJ

      DO 140 K=1,6
      QQDUST(K)=(WTI*SRDQEX(K,I)+WTJ*SRDQEX(K,J))*FNDUST
      SSDUST(K)=(WTI*SRDQSC(K,I)+WTJ*SRDQSC(K,J))*FNDUST
      GGDUST(K)=(WTI*SRDQCB(K,I)+WTJ*SRDQCB(K,J))
  140 CONTINUE
      Q55DST=WTI*Q55D25(I)+WTJ*Q55D25(J)

      RETURN
      END SUBROUTINE SDDUST

      SUBROUTINE TDDUST(ADDUST,SIZDST,FNDUST)

      !INCLUDE 'BR00B.COM'

C     ------------------------------------------------------------------
C     TDDUST   Selects Thermal Mie Scattering Parameters for Desert Dust
C                                     and applies FTDUST Scaling Factors
C                                     ----------------------------------

      DIMENSION QDDUST(33),SDDUST(33),           ADDUST(33)
CCC   DIMENSION QDDUST(33),SDDUST(33),GDDUST(33),ADDUST(33)

      RDI=0.0
      DO 110 J=1,25
      RDJ=REFD25(J)
      IF(SIZDST.LT.RDJ) GO TO 120
      RDI=RDJ
  110 CONTINUE
      WTJ=1.0
      J=25
      I=25
      GO TO 130
  120 WTJ=(SIZDST-RDI)/(RDJ-RDI)
      I=J-1
      IF(I.LT.1) I=1
  130 WTI=1.0-WTJ

      DO 140 K=1,33
      QDDUST(K)=(WTI*TRDQEX(K,I)+WTJ*TRDQEX(K,J))*FNDUST
      SDDUST(K)=(WTI*TRDQSC(K,I)+WTJ*TRDQSC(K,J))*FNDUST
CCC   GDDUST(K)=(WTI*TRDQCB(K,I)+WTJ*TRDQCB(K,J))
      ADDUST(K)=QDDUST(K)-SDDUST(K)
  140 CONTINUE

      RETURN
      END SUBROUTINE TDDUST

      SUBROUTINE AO3ABS(OCM,O3ABS)

      !INCLUDE 'BR00B.COM'

C              ---------------------------------------------------------
C              UV absorption by Ozone  is expressed as a fraction of the
C              total solar flux S0. Hence O3ABS (fraction of total solar
C              flux absored by OCM cm ofozone) must be normalized within
C              SOLARM by dividing O3ABS by the corresponding fraction of
C              the solar flux within the spectral interval DKS0(15)=0.05
C              ---------------------------------------------------------

      O3ABS=AO3(460)
      IP=0
      XX=OCM*1.D+04
      IX=XX
      IF(IX.GT.99) GO TO 110
      IF(IX.LT.1 ) GO TO 130
      GO TO 120
  110 CONTINUE
      IP=IP+90
      XX=XX*0.1D0
      IX=XX
      IF(IX.GT.99) GO TO 110
  120 CONTINUE
      DX=XX-IX
      IX=IX+IP
      IF(IX.GT.459) GO TO 140
      O3ABS=AO3(IX)+DX*(AO3(IX+1)-AO3(IX))
      GO TO 140
  130 CONTINUE
      O3ABS=XX*AO3(1)
  140 CONTINUE

      RETURN
      END SUBROUTINE AO3ABS

      SUBROUTINE WRITER(KWRU,INDEX)
C
      USE SURF_ALBEDO, only : AVSCAT, ANSCAT, AVFOAM, ANFOAM,
     *     WETTRA, WETSRA, ZOCSRA, ZSNSRA, ZICSRA, ZDSSRA, ZVGSRA,
     *     EOCTRA, ESNTRA, EICTRA, EDSTRA, EVGTRA, AGEXPF, ALBDIF
      !INCLUDE 'BR00B.COM'   ! use RE001
C
C     ------------------------------------------------------------------
C     WRITER Radiative Input/Output Cloud/Aerosol Data/Conrol Parameters
C
C         INDEX
C           0     control parameter defaults in RADPAR
C           1     RADCOM Radiative conrol/scaling parameter GHG defaults
C           2     RADCOM Atmospheric composition P,H,T,Cld,Aer  profiles
C           3     RADCOM Computed LW SW fluxes cooling and heating rates
C           4     Aerosol and Cloud: Mie scattering radiative parameters
C                 A  SW aerosol Mie scattering Qx,Qs,g in use parameters
C                 B  SW  cloud  Mie scattering Qx,Qs,g in use parameters
C                 C  SW cld+aer Mie scattering Qx,Qs,g in use parameters
C                 D  SW LW aerosol 11-compositon  Mie Qx,Qs,g parameters
C                 E  SW LW aerosol 10-compositon  Mie Qx,Qs,g parameters
C                 F  SW LW aerosol 8-size D dust  Mie Qx,Qs,g parameters
C                 G  SW LW  cloud  15-size/phase  Mie Qx,Qs,g parameters
C           5     LW cld,aer,gas total optical k-distribution extinction
C           6     LW gas absorb: total optical k-distribution extinction
C           7     A  LW cloud   TRCALK optical k-distribution extinction
C                 B  LW aerosol TRAALK optical k-distribution extinction
C           8     SW Spectral/k-dist flux, albedo, absorption components
C                 A  Spectral components of downward & upward solar flux
C                 B  Spectral components of net solar flux, heating rate
C           9     LW flux contribution from each k-distribution interval
C                 1  Downward LW flux  from each k-distribution interval
C                 2  Upward   LW flux  from each k-distribution interval
C                 3  Net (Up) LW flux  from each k-distribution interval
C                 4  Flux cooling rate from each k-distribution interval
C                 5  Fraction coolrate from each k-distribution interval
C        NOTE:
C                 KWTRAB sets LW Mie parameters in 4-D,E,F,G
C                        KWTRAB=0 (default) sets LW output to be Mie Qab
C                        KWTRAB=1           sets LW output to be Mie Qex
C                        KWTRAB=2           sets LW output to be Mie Qsc
C                        KWTRAB=3           sets LW output to be Mie Qcb
C                        KWTRAB=4           sets LW output to be Mie Pi0
C
C                 INDEX  0-9 : show item 'INDEX' only
C                 INDEX 11-19: show items 1->last digit of 'INDEX'
C                 INDEX 21-29: show items 0->last digit of 'INDEX'
C                 KWRU directs the output to selected (KWRU) file number
C     ------------------------------------------------------------------
C
      INTEGER, INTENT(IN) :: INDEX
      CHARACTER*8 FTYPE
      CHARACTER*6 GHG
      CHARACTER*3 TRABCD,TRAXSG
      DIMENSION TRABCD(5),TRAXSG(5),TRPI0K(25)
      DATA TRABCD/'TRA','TRB','TRC','TRD','TRE'/
      DATA TRAXSG/'QAB','QEX','QSC','QCB','PI0'/
C
      DIMENSION TKEFF(3)
C
      DIMENSION BGFLUX(33),BGFRAC(33),TAUSUM(33)
      DIMENSION SUM0(20),SUM1(LX),SUM2(LX),SUM3(LX),FTYPE(5),GHG(12)
      DIMENSION WSREXT(LX,6),WSRSCT(LX,6),WSRGCB(LX,6),WSRPI0(LX,6)
      DIMENSION FSR1(17),FSR2(17),ISR1(16),KSLAMW(16),IORDER(16)
      DATA KSLAMW/1,1,2,2,5,5,5,5,1,1,1,3,4,6,6,1/
      DATA IORDER/12,11,10, 9, 6, 5, 4, 3,15,14,13, 8, 7, 2, 1,16/
C
      DATA FTYPE/'DOWNWARD','  UPWARD','UPWD NET','COOLRATE','FRACTION'/
      DATA GHG/'   H2O','   CO2','    O3','    O2','   NO2','   N2O'
     +        ,'   CH4','CCL3P1','CCL2P2','    N2',' CFC-Y',' CFC-Z'/
      character*1,dimension(4) :: AUXGAS = (/'0','L','X','X'/)
CCC   DATA P0/1013.25/
      PARAMETER(P0=1013.25)
      SAVE  KSLAMW,IORDER,FTYPE,AUXGAS
C
      DO 20 K=1,6
      DO 10 L=1,NL
      WSREXT(L,K)=SRAEXT(L,K)+SRBEXT(L,K)
     +           +SRDEXT(L,K)+SRVEXT(L,K)
      WSRSCT(L,K)=SRASCT(L,K)+SRBSCT(L,K)
     +           +SRDSCT(L,K)+SRVSCT(L,K)
      WSRGCB(L,K)=SRASCT(L,K)*SRAGCB(L,K)+SRBSCT(L,K)*SRBGCB(L,K)
     +           +SRDSCT(L,K)*SRDGCB(L,K)+SRVSCT(L,K)*SRVGCB(L,K)
      WSRPI0(L,K)=WSRSCT(L,K)/(WSREXT(L,K)+1.E-10)
      WSRGCB(L,K)=WSRGCB(L,K)/(WSRSCT(L,K)+1.D-10)
   10 CONTINUE
   20 CONTINUE
C
      ACOLX=0.0
      BCOLX=0.0
      DCOLX=0.0
      VCOLX=0.0
      DO 30 L=1,NL
      ACOLX=ACOLX+SRAEXT(L,6)
      BCOLX=BCOLX+SRBEXT(L,6)
      DCOLX=DCOLX+SRDEXT(L,6)
      VCOLX=VCOLX+SRVEXT(L,6)
   30 CONTINUE
      TCOLX=ACOLX+BCOLX+DCOLX+VCOLX
C
      KW=KWRU
      INDJ=MOD(INDEX,10)
      IF(INDJ.LT.1.and.INDEX.gt.0) INDJ=10
      INDI=1
      IF(INDEX.GT.20.or.INDEX.eq.0) INDI=0
      IF(INDEX.LT.11) INDI=INDJ
      DO 9999 INDX=INDI,INDJ
C
      KPAGE=1
      IF(INDX.EQ.0) GO TO 90
      GO TO (100,200,300,400,500,600,700,800,900,1000),INDX
C
C-------------
   90 CONTINUE
C-------------
C
      WRITE(KW,6000)
 6000 FORMAT(' CALL WRITER(KW,0) :',2X,'PAGE 1/2  '
     +          ,'CONTROL PARAMS   DEFINITIONS'/
     +      /' CONTROL PARAMTER      DEFAULT  PARAMETER DESCRIPTION')
      WRITE(KW,6001)  KVRAER,MEANAC,MEANDD,MEANVA,NCARO3,KUVFAC,KSNORM
     + ,KWTRAB,KGGVDF,KPGRAD,KLATZ0,KCLDEM,KANORM,KPFCO2,KPFOZO,KVEGA6
     + ,KORDER,KUFH2O,KUFCO2,KCSELF,KCFORN
 6001 FORMAT(7X,'   KVRAER = ',I1,'     1      Repartition Aer VDist'
     2      /7X,'   MEANAC = ',I1,'     0      Use Ann-Mean Aer Clim'
     3      /7X,'   MEANDD = ',I1,'     0      Use Ann-Mean Des Dust'
     4      /7X,'   MEANVA = ',I1,'     0      Use Ann-Mean Volc Aer'
     5      /7X,'   NCARO3 = ',I1,'     0      NCAR London 1976 Ozon'
     6      /7X,'   KUVFAC = ',I1,'     0      ON/OFF UV Mult Factor'
     7      /7X,'   KSNORM = ',I1,'     0      Norm S0 when KUVFAC=1'
     8      /7X,'   KWTRAB = ',I1,'     0      WRITER: Qab,Qex,Qsc,g'
     9      /7X,'   KGGVDF = ',I1,'     0      Use GHG VertProf Grad'
     A      /7X,'   KPGRAD = ',I1,'     1      Pole-to-Pole GHG Grad'
     1      /7X,'   KLATZ0 = ',I1,'     1      Use GHG VDist Lat Dep'
     2      /7X,'   KCLDEM = ',I1,'     1      Use TopCloud Scat Cor'
     3      /7X,'   KANORM = ',I1,'     0      Use SGP Atmo Col Norm'
     4      /7X,'   KPFCO2 = ',I1,'     0      1=MOD CO2PROF: FPXCO2'
     5      /7X,'   KPFOZO = ',I1,'     0      1=MOD O3 PROF: FPXOZO'
     6      /7X,'   KVEGA6 = ',I1,'     0      Schramm"s ocn ice alb'
     7      /7X,'   KORDER = ',I1,'     0      WRITER k-d spec order'
     8      /7X,'   KUFH2O = ',I1,'     1      Col Absorber Scal H2O'
     +      /7X,'   KUFCO2 = ',I1,'     1      Col Absorber Scal CO2'
     1      /7X,'   KCSELF = ',I1,'     1      H2O Cont Self-Broaden'
     2      /7X,'   KCFORN = ',I1,'     1      H2O Con Foreign-Broad'
     R      )
C
      WRITE(KW,6004)
 6004 FORMAT(/' CONTROL PARAMTER    DEFAULT       SNOW/ICE FACTORS')

      WRITE(KW,6005) agexpf,albdif
 6005 FORMAT(7X   ,'   AGEXPF = ',F7.3,'      SNOWAGE XPFACTOR SH EARTH'
     A      /7X   ,'    SH O  = ',F7.3,'         "       "     "  OCICE'
     B      /7X   ,'    SH L  = ',F7.3,'         "       "     "  LDICE'
     C      /7X   ,'    NH E  = ',F7.3,'         "       "     NH EARTH'
     D      /7X   ,'    NH O  = ',F7.3,'         "       "     "  OCICE'
     E      /7X   ,'    NH L  = ',F7.3,'         "       "     "  LDICE'
     F      /7X   ,'   ALBDIF = ',F7.3,'      SNOW/ICE ALBDIF  SH EARTH'
     G      /7X   ,'    SH O  = ',F7.3,'         "       "     "  OCICE'
     H      /7X   ,'    SH L  = ',F7.3,'         "       "     "  LDICE'
     I      /7X   ,'    NH E  = ',F7.3,'         "       "     NH EARTH'
     J      /7X   ,'    NH O  = ',F7.3,'         "       "     "  OCICE'
     K      /7X   ,'    NH L  = ',F7.3,'         "       "     "  LDICE'
     L      )
C
      WRITE(KW,6006)
 6006 FORMAT('0CONTROL PARAMTER    VALUE',16X,' DEFAULT')
      WRITE(KW,6007) PTOPTR,REFF0 ,VEFF0 ,O3WJT0,X0YBCI,X0YOCI,X0YSUI
     +              ,AVSCAT,ANSCAT,AVFOAM,ANFOAM
 6007 FORMAT(7X,'   PTOPTR = ',F7.1,6X,'TropTop (SIGMA lev) Pressure'
     A      /7X,'   REFF0  = ',F7.3,'                 0.300         '
     B      /7X,'   VEFF0  = ',F7.3,'                 0.350         '
     C      /7X,'   O3WJT0 = ',F7.5,'                 0.001         '
     E      /7X,'   X0YBCI = ',F7.3,'                 0.001         '
     F      /7X,'   X0YOCI = ',F7.3,'                 0.001         '
     G      /7X,'   X0YSUI = ',F7.3,'                 0.0001        '
     H      /7X,'   AVSCAT = ',F7.5,'                 0.01560       '
     I      /7X,'   ANSCAT = ',F7.5,'                 0.00020       '
     J      /7X,'   AVFOAM = ',F7.5,'                 0.21970       '
     K      /7X,'   ANFOAM = ',F7.5,'                 0.15140       '
     X      )
      WRITE(KW,6008)
 6008 FORMAT(/10X,'UV Solar Flux Spectral Partitions and Factors')
      WRITE(KW,6009) UVWAVL,UVFACT
 6009 FORMAT(10X,'UVWAVL = ',F7.5,2F8.5/10X,'UVFACT = ',F7.5,2F8.5)
C
      WRITE(KW,6013)
 6013 FORMAT(/' CONTROL PARAMETER  PI0VIS     PI0TRA      DEFAULT')
      WRITE(KW,6014) PI0VIS,PI0TRA
 6014 FORMAT(7X,'   ACID1 = ',F8.6,F11.6,'     1.0                 '
     A      /7X,'   SSALT = ',F8.6,F11.6,'     1.0                 '
     B      /7X,'   SLFT1 = ',F8.6,F11.6,'     1.0                 '
     C      /7X,'   SLFT2 = ',F8.6,F11.6,'     1.0                 '
     D      /7X,'   BSLT1 = ',F8.6,F11.6,'      .98929             '
     E      /7X,'   BSLT2 = ',F8.6,F11.6,'      .95609             '
     F      /7X,'   DUST1 = ',F8.6,F11.6,'      .91995             '
     G      /7X,'   DUST2 = ',F8.6,F11.6,'      .78495             '
     H      /7X,'   DUST3 = ',F8.6,F11.6,'      .63576             '
     I      /7X,'   CARB1 = ',F8.6,F11.6,'      .31482             '
     J      /7X,'   CARB2 = ',F8.6,F11.6,'      .47513             '
     K      )

      WRITE(KW,6017)
 6017 FORMAT(/'  layer  VDFBCI VDFOCI VDFSUI VDFDST')
      WRITE(KW,6018) (I,VDFBCI(I),VDFOCI(I),VDFSUI(I),VDFDST(I),I=1,12)
 6018 FORMAT(3X,I4,' :',4F6.2)

      WRITE(KW,6019)
 6019 FORMAT(/'  GHGAS',9X,'PPMVK0    PPMVDF    PPGRAD')
      WRITE(KW,6020) (ghg(I),PPMVK0(I),PPMVDF(I),PPGRAD(I),I=1,12)
 6020 FORMAT(1X,a6,' ',F15.7,F10.5,F10.5)
      GO TO 9999
C
C-------------
  100 CONTINUE
C-------------
C
      NPAGE=1
      IF(INDEX.LT.11) NPAGE=KPAGE
      WRITE(KW,6101)
      WRITE(KW,6102)
      FACTOR=P0/(PLB(1)-PLB(2))*1.25
      PPMCO2=ULGAS(1,2)*FACTOR
      PPMO2 =ULGAS(1,4)*FACTOR
      PPMN2O=ULGAS(1,6)*FACTOR
      PPMCH4=ULGAS(1,7)*FACTOR
      PPMF11=ULGAS(1,8)*FACTOR
      PPMF12=ULGAS(1,9)*FACTOR
      PPMY11=ULGAS(1,11)*FACTOR
      PPMZ12=ULGAS(1,12)*FACTOR
      WRITE(KW,6103) (FULGAS(I),I=1,3),(FULGAS(I),I=6,9)
     +              ,FULGAS(11),FULGAS(12),    (FGOLDH(I),I=1,5)
C     IF(KGASSR.GT.0)
C    +WRITE(KW,6104) (FULGAS(I+9),I=1,2),(FULGAS(I+9),I=4,9)
C    +              ,FULGAS(11),FULGAS(12),   (FGOLDH(I+9),I=1,5)
      WRITE(KW,6105) PPMCO2,PPMN2O,PPMCH4,PPMF11,PPMF12,PPMY11,PPMZ12
     +             ,(FGOLDH(I),I=6,9),PPMV80(2),(PPMV80(I),I=6,9)
     +             ,(PPMV80(I),I=11,12),KTREND,JYEAR,JDAY,LASTVC
      WRITE(KW,6106) TAUWC0,FCLDTR,EOCTRA,ZOCSRA,KZSNOW,KCLDEM,NTRACE
     +             ,FSAAER,FTTAER,MOZONE,KCLDEP,NL
      WRITE(KW,6107) TAUIC0,FCLDSR,ESNTRA,ZSNSRA,WETTRA,KVEGA6,ITR(1)
     +             ,ITR(5),FSBAER,FTBAER,KO3LON,KEEPAL,NLP
      WRITE(KW,6108)FULGAS(13),FRAYLE,EICTRA,ZICSRA,WETSRA,KCNORM,ITR(2)
     +             ,ITR(6),FSAAER,FTAAER,KVRAER,KEEP10,MLAT46
      WRITE(KW,6109) TLGRAD,ECLTRA,EDSTRA,ZDSSRA,KANORM,KPGRAD,ITR(3)
     +             ,ITR(7),FSDAER,FTDAER,KWVCON,ICE012,MLON72
      WRITE(KW,6110) PTLISO,ZCLSRA,EVGTRA,ZVGSRA,KEEPRH,KLATZ0,ITR(4)
     +             ,ITR(8),FSVAER,FTVAER,KSOLAR,NORMS0
C
 6101 FORMAT(' (1)FUL:  1',7X,'2',8X,'3',7X,'6',7X,'7',8X,'8',8X,'9'
     +      ,8X,'11',7X,'12',4X,'RADCOM 1/F: (Control/Default'
     +      ,'/Scaling Parameters)')
 6102 FORMAT(4X,'GAS: ','H2O',5X,'CO2',7X,'O3'
     +      ,5X,'N2O',5X,'CH4',5X,'CFC-11',3X,'CFC-12'
     +      ,3X,'CFY-11',3X,'CFZ-12'
     +      ,2X,'Aerosol Global    Ocean     Land  Desert    Haze')
 6103 FORMAT(1X,'FULGAS=',F5.3,F10.5,F7.3,F9.5,F8.5,4F9.5
     +      ,2X,'FGOLDH=',F7.5,2F9.6,2F8.5)
 6104 FORMAT(1H+,T84,'T'
     +      /1X,'FULGAS=',1P,1E7.1,1P,2E8.1,1P,2E8.1,1P,4E9.1
     +      ,' S','FGOLDH=',1P,1E7.1,1P,2E9.2,1P,2E8.1)
 6105 FORMAT(1X,'PPM(1)=(now)',2X,F8.3,8X,F8.5,F8.5,4(1X,F8.7)
     +      ,2X,'TRACER=',F7.5,2F9.6,F8.5
     +      /' PPMV80=(ref)=',0P,F9.3,8X,2F8.5,4(1X,F8.7),2X
     +      ,'KTREND=',I1,2X,'JYEAR=',I4,' JDAY=',I3,5X,'LASTVC=',I7)
 6106 FORMAT(1X,'TAUWC0=',1P,E6.0,' FCLDTR=',0P,F4.2,' EOCTRA=',F3.1
     +      ,1X,'ZOCSRA=',   F3.1,' KZSNOW=',     I4,' KCLDEM=',  I3
     +      ,1X,'NTRACE=',    I3,2X,'FSTAER=',  F3.1,' FTTAER=',F3.1
     +      ,1X,'MOZONE=',    I1,1X,'KCLDEP=',    I1,'     NL=',  I2)
 6107 FORMAT(1X,'TAUIC0=',1P,E6.0,' FCLDSR=',0P,F4.2,' ESNTRA=',F3.1
     +      ,1X,'ZSNSRA=',  F3.1,1X,'WETTRA=',  F4.2,' KVEGA6=',  I3
     +      ,1X,'ITR(1)=',    2I2,1X,'FSBAER=', F3.1,' FTBAER=',F3.1
     +      ,1X,'K03LON=',    I1,1X,'KEEPAL=',    I1,'    NLP=',  I2)
 6108 FORMAT(1X,'FO3LON=',   F6.2,' FRAYLE=',0P,F4.1,' EICTRA=',F3.1
     +      ,1X,'ZICSRA=',  F3.1,1X,'WETSRA=',  F4.2,' KCNORM=',  I3
     +      ,1X,'ITR(2)=',    2I2,1X,'FSAAER=', F3.1,' FTAAER=',F3.1
     +      ,1X,'KVRAER=',    I1,1X,'KEEP10=',    I1,' MLAT46=',  I2)
 6109 FORMAT(1X,'TLGRAD=',   F6.2,' ECLTRA=',0P,F4.2,' EDSTRA=',F3.1
     +      ,1X,'ZDSSRA=',  F3.1,1X,'KANORM=',  I4  ,' KPGRAD=',  I3
     +      ,1X,'ITR(3)=',    2I2,1X,'FSDAER=', F3.1,' FTDAER=',F3.1
     +      ,1X,'KWVCON=',    I1,1X,'ICE012=',    I1,' MLON72=',  I2)
 6110 FORMAT(1X,'PTLISO=',  F6.1,1X,'ZCLSRA=',  F4.2,' EVGTRA=',F3.1
     +      ,1X,'ZVGSRA=',  F3.1,1X,'KEEPRH=',  I4  ,' KLATZ0=',  I3
     +      ,1X,'ITR(4)=',    2I2,1X,'FSVAER=', F3.1,' FTVAER=',F3.1
     +      ,1X,'KSOLAR=',    I1,1X,'NORMS0=',    I1,'        ')
      GO TO 9999
C
C-------------
  200 CONTINUE
C-------------
C
      NPAGE=0
      LUXGAS=0
      IF(INDEX.LT.11) NPAGE=KPAGE
      WRITE(KW,6201) AUXGAS(LUXGAS+1),S00WM2,S0,COSZ
      DO 202 K=1,9
      DO 201 L=1,NL
      UXGAS(L,K)=ULGAS(L,K)
  201 CONTINUE
  202 CONTINUE
      IF(LUXGAS.LT.2) GO TO 205
      LGS=(LUXGAS-2)*9
      DO 203 L=1,NL
      UXGAS(L,1)=U0GAS(L,1)*FULGAS(1+LGS)
      UXGAS(L,3)=U0GAS(L,3)*FULGAS(3+LGS)
  203 UXGAS(L,5)=U0GAS(L,5)*FULGAS(5+LGS)
C
      DO 204 L=1,NL
      UXGAS(L,2)=U0GAS(L,2)*FULGAS(2+LGS)
      UXGAS(L,4)=U0GAS(L,4)*FULGAS(4+LGS)
      UXGAS(L,6)=U0GAS(L,6)*FULGAS(6+LGS)
      UXGAS(L,7)=U0GAS(L,7)*FULGAS(7+LGS)
      UXGAS(L,8)=U0GAS(L,8)*FULGAS(8+LGS)
  204 UXGAS(L,9)=U0GAS(L,9)*FULGAS(9+LGS)
  205 CONTINUE
      DO 206 L=NL,1,-1
      EPS=CLDEPS(L)
      TAER=WSREXT(L,6)
      IPI0=WSRPI0(L,6)*1000.D0+1.D-05
      HLM=0.5D0*(HLB0(L+1)+HLB0(L))
      TLAPS=(TLT(L)-TLB(L))/(HLB0(L+1)-HLB0(L))
      IRHL=RHL(L)*100.0
      IF(PL(L).LT.1.D0) THEN
      WRITE(KW,6212) L,PL(L),HLM,TLM(L),TLAPS,SHL(L),IRHL
     +     ,(UXGAS(L,K),K=1,3),(UXGAS(L,K),K=6,9),UXGAS(L,5)
     +     ,SIZEWC(L),SIZEIC(L),TAUWC(L),TAUIC(L),EPS,TAER,IPI0
      ELSE
      IF(UXGAS(L,1).GE.1.D0) THEN
      WRITE(KW,6202) L,PL(L),HLM,TLM(L),TLAPS,SHL(L),IRHL
     +     ,(UXGAS(L,K),K=1,3),(UXGAS(L,K),K=6,9),UXGAS(L,5)
     +     ,SIZEWC(L),SIZEIC(L),TAUWC(L),TAUIC(L),EPS,TAER,IPI0
      ELSE
      WRITE(KW,6211) L,PL(L),HLM,TLM(L),TLAPS,SHL(L),IRHL
     +     ,(UXGAS(L,K),K=1,3),(UXGAS(L,K),K=6,9),UXGAS(L,5)
     +     ,SIZEWC(L),SIZEIC(L),TAUWC(L),TAUIC(L),EPS,TAER,IPI0
      ENDIF
      ENDIF
  206 CONTINUE
      DO 207 I=1,16
  207 SUM0(I)=0.
      DO 210 L=1,NL
      DO 208 I=1,9
      SUM0(I)=SUM0(I)+ULGAS(L,I)
  208 CONTINUE
      DO 209 I=1,4
      SUM0(12+I)=SUM0(12+I)+TRACER(L,I)
  209 CONTINUE
      SUM0(10)=SUM0(10)+TAUWC(L)
      SUM0(11)=SUM0(11)+TAUIC(L)
  210 CONTINUE
      TAU55=0.0
      DO 211 L=1,NL
      TAU55=TAU55+WSREXT(L,6)
  211 CONTINUE
      SUM0(12)=TAU55
      TGMEAN=POCEAN*TGO**4+PEARTH*TGE**4+PLICE*TGLI**4+POICE*TGOI**4
      TGMEAN=SQRT(TGMEAN)
      TGMEAN=SQRT(TGMEAN)
      WRITE(KW,6203) (SUM0(I),I=1,3),(SUM0(I),I=6,9),SUM0(5)
     +               ,SUM0(10),SUM0(11),SUM0(12)
      WRITE(KW,6204) POCEAN,TGO,AGESN(1),WMAG,SUM0(13),JYEAR
     +             ,BXA(4:5),LASTVC
      WRITE(KW,6205) PEARTH,TGE,SNOWE,WEARTH,SUM0(14),JDAY,BXA(6:7)
      WRITE(KW,6206) POICE,TGOI,SNOWOI,      SUM0(15),JLAT
     +             ,(SRBALB(I),I=1,6)
      WRITE(KW,6207) PLICE,TGLI,SNOWLI,FULGAS(10),SUM0(16),ILON
     +             ,(SRXALB(I),I=1,6)
      PSUM=POCEAN+PEARTH+POICE+PLICE
      WRITE(KW,6208) TGMEAN,FULGAS(4),FULGAS(13)
     +               ,PSUM,TSL,FULGAS(5),(PVT(I),I=1,11)
      write(kw,*) 'hsn,hin,hmp,fmp, flags, LS1_loc',
     * hsn,hin,hmp,fmp, flags, LS1_loc
      WRITE(KW,6209) (PRNB(1:2,I),PRNX(1:2,I),I=1,4),BXA(1:3)
      WRITE(KW,6210)
 6201 FORMAT(' (2) RADCOM G/L: (Input Data)'
     +      ,2X,'Absorber Amount per Layer:'
     +      ,'  U',1A1,'GAS(L,K) in cm**3(STP)/cm**2',2X,'S00WM2=',F9.4
     +      ,1X,'S0=',F9.4,2X,'COSZ=',F6.4
     +      /' LN     PL  HLM   TLM  TLAP   SHL  .RH     '
     +      ,'H2O   CO2    O3    N2O    CH4  CFC-11'
     +      ,'  CFC-12   NO2   WC.SIZ.IC  WC.TAU.IC CLEP  A TAU PI0')
 6202 FORMAT(1X,I2,F7.2,F5.1,F7.2,F5.1,1X,F7.6,I3
     +      ,F8.2,F6.2,1X,F6.5,1X,F5.4
     +      ,F7.4,1P,3E8.1,0P,2F5.1,F6.2,F5.2,1X,F4.3,F6.3,I5)
 6203 FORMAT(24X,' Column Amount',F7.1,F7.2,1X,F6.5
     +       ,1X,F5.4,F7.4,1P,3E8.1,0P,10X,F6.2,F5.2,5X,F6.3)
 6204 FORMAT( 1X,'POCEAN=',F6.4,'    TGO=' ,F6.2,1X,' AGESN=',F6.3
     +      , 1X,'  WMAG=',F6.3,' TRACER 1=',F5.3,' JYEAR=',I4
     +      , 3X,'BSNVIS=',F6.4,' BSNNIR=' ,F6.4,7X,'LASTVC=',I7)
 6205 FORMAT(    ' PEARTH=',F6.4,'    TGE=',F6.2,'  SNOWE=',F6.3
     +      ,    ' WEARTH=',F6.3,'  Sums: 2=',F5.3
     +      ,     '  JDAY=',I4  ,2X,' XSNVIS=',F6.4,' XSNNIR=',F6.4
     +      , 8X,'NIRALB VISALB')
 6206 FORMAT(    '  POICE=',F6.4,'   TGOI=',F6.2,' SNOWOI=',F6.3
     +      ,                14X,'        3=',F5.3
     +      ,     '  JLAT=',I4,  2X,' SRBALB=',F6.4
     +      ,4F7.4,F7.4)
 6207 FORMAT(    '  PLICE=',F6.4,'   TGLI=',F6.2,' SNOWLI=',F6.3
     +      ,    ' FULGAS(10)=',F3.1,'       4=',F5.3
     +      ,     '  ILON=',I4,  2X,' SRXALB=',F6.4
     +      ,4F7.4,F7.4)
 6208 FORMAT(8X,6('-'),' TGMEAN=',F6.2,' FULGAS(4)=',F3.1,' FULGAS(13)='
     +      ,F3.1,' BSAND TUNDRA GRASSL SHRUBS  TREES DECIDF'
     +       ,' EVERGF','  RAINF','  CROPS','  BDIRT','  ALGAE'
     +      /    '   PSUM=',F6.4,'    TSL=',F6.2,' FULGAS(5)=',F3.1
     +       ,9X,T54,'PVT=',F6.4,10F7.4)
 6209 FORMAT(' BOCVIS BOCNIR XOCVIS XOCNIR]BEAVIS BEANIR XEAVIS XEANIR'
     +      ,'|BOIVIS BOINIR XOIVIS XOINIR]BLIVIS BLINIR XLIVIS XLINIR'
     +      ,'|EXPSNE]EXPSNO]EXPSNL'/1X,F6.4,18F7.4)
 6210 FORMAT(' ')
 6211 FORMAT(1X,I2,F7.2,F5.1,F7.2,F5.1,1X,F7.6,I3
     +      ,F8.5,F6.2,1X,F6.5,1X,F5.4
     +      ,F7.4,1P,3E8.1,0P,2F5.1,F6.2,F5.2,1X,F4.3,F6.3,I5)
 6212 FORMAT(1X,I2,F7.4,F5.1,F7.2,F5.1,1X,F7.6,I3
     +      ,F8.5,F6.2,1X,F6.5,1X,F5.4
     +      ,F7.4,1P,3E8.1,0P,2F5.1,F6.2,F5.2,1X,F4.3,F6.3,I5)
C
      GO TO 9999
C
C-------------
  300 CONTINUE
C-------------
C
      NPAGE=0
      IF(INDEX.LT.11) NPAGE=KPAGE
      IF(NL.GT.13) NPAGE=1
      L=NLP
      SRALB =SRUFLB(L)/(SRDFLB(L)+1.E-10)
      STNFLB=SRNFLB(L)-TRNFLB(L)
      WRITE(KW,6301) NORMS0
      WRITE(KW,6302) L,PLB(L),HLB0(L),TLT(L-1) ! TLB(LN+1) unused/set
     +             ,TRDFLB(L),TRUFLB(L),TRNFLB(L)
     +             ,SRDFLB(L),SRUFLB(L),SRNFLB(L),STNFLB,SRALB
      DO 301 N=1,NL
      L=NLP-N
      CRHRF=8.4167/(PLB(L)-PLB(L+1))
      STNFLB=SRNFLB(L)-TRNFLB(L)
      STFHR =SRFHRL(L)-TRFCRL(L)
      TRDCR =TRFCRL(L)*CRHRF
      SRDHR =SRFHRL(L)*CRHRF
      STDHR=STFHR*CRHRF
      SRALB =SRUFLB(L)/(SRDFLB(L)+1.E-10)
      SRXVIS=SRXATM(1)
      SRXNIR=SRXATM(2)
      IF(PLB(L).LT.1.D0) THEN
      WRITE(KW,6313) L,PLB(L),HLB0(L),TLB(L),TLT(L)
     +             ,TRDFLB(L),TRUFLB(L),TRNFLB(L),TRFCRL(L)
     +             ,SRDFLB(L),SRUFLB(L),SRNFLB(L),SRFHRL(L)
     +             ,STNFLB,STFHR,STDHR,TRDCR,SRDHR,SRALB
      ELSE
      WRITE(KW,6303) L,PLB(L),HLB0(L),TLB(L),TLT(L)
     +             ,TRDFLB(L),TRUFLB(L),TRNFLB(L),TRFCRL(L)
     +             ,SRDFLB(L),SRUFLB(L),SRNFLB(L),SRFHRL(L)
     +             ,STNFLB,STFHR,STDHR,TRDCR,SRDHR,SRALB
      ENDIF
  301 CONTINUE
C
      DO 302 II=1,3
      PFW=TRDFLB(1)
      IF(II.EQ.2) PFW=TRUFLB(1)
      IF(II.EQ.3) PFW=TRUFLB(NLP)
      IPF=PFW
      DPF=PFW-IPF
      IF(IPF.LT.1) IPF=1
      IF(IPF.GT.899) IPF=899
  302 TKEFF(II)=TKPFT(IPF)+DPF*(TKPFT(IPF+1)-TKPFT(IPF))
C
      WRITE(KW,6304) WINDZF(1),WINDZT(1),TOTLZF(1),TOTLZT(1),
     +              (FSRNFG(I),I=1,4),LTOPCL,JLAT,JYEAR
      WRITE(KW,6305) WINDZF(2),WINDZT(2),TOTLZF(2),TOTLZT(2),
     +              (FTRUFG(I),I=1,4),LBOTCL,ILON,JDAY
      IF(KORDER.EQ.0) WRITE(KW,6306) WINDZF(3),WINDZT(3),TOTLZF(3)
     +              ,TOTLZT(3),(I,I=1,16)
      IF(KORDER.EQ.1) WRITE(KW,6307) WINDZF(3),WINDZT(3),TOTLZF(3)
     +              ,TOTLZT(3),(I,I=1,16)
      FRACSL=0.D0
      IF(KORDER.EQ.0) WRITE(KW,6308) TKEFF(1),TKEFF(2),TKEFF(3)
     +               ,(SRKALB(NORDER(I)),I=1,16),BTEMPW,TRUFTW,SRIVIS
     +               ,SROVIS,PLAVIS,SRINIR,SRONIR,PLANIR
      IF(KORDER.EQ.1) WRITE(KW,6308) TKEFF(1),TKEFF(2),TKEFF(3)
     +               ,(SRKALB(I),I=1,16),BTEMPW,TRUFTW,SRIVIS
     +               ,SROVIS,PLAVIS,SRINIR,SRONIR,PLANIR
      WRITE(KW,6309) TRDFGW,TRUFGW,SRDVIS,SRUVIS,ALBVIS,SRDNIR,SRUNIR
     +             ,ALBNIR
      WRITE(KW,6310) SRXVIS,SRXNIR,SRTVIS,SRRVIS,SRAVIS,SRTNIR,SRRNIR
     +             ,SRANIR
C
C
 6301 FORMAT(/' (3) RADCOM M/S: (Output Data)'
     +      ,T37,'Thermal Fluxes (W/M**2)',4X,'Solar Fluxes (W/M**2)'
     +      ,1X,'NORMS0=',I1,'  Energy  Input  Heat/Cool Deg/Day Alb'
     +      ,'do'/' LN     PLB   HLB    TLB    TLT '
     +      ,'  TRDFLB TRUFLB TRNFLB TRFCRL   SRDFLB  SRUFLB  SRNFLB'
     +      ,' SRFHRL  STNFLB  STFHR  SR-TR TR=CR SR=HR SRALB')
 6302 FORMAT(1X,I2,F9.3,F6.2,1X,F6.2,8X,3F7.2,8X,3F8.2,7X,F8.2,26X,F6.4)
 6303 FORMAT(1X,I2,F9.3,F6.2,2F7.2,1X,3F7.2,F7.2,1X,3F8.2,F7.2,1X,F7.2
     +      ,1X,F6.2,1X,3F6.2,1X,F5.4)
 6304 FORMAT( 1X,'XMU  WINDZF WINDZT   TOTLZF TOTLZT'/
     +        1X,'1.0',1X,F7.3,F7.2,2X,F7.3,F7.2,2X,'FR.SRNLB1'
     +      ,' OCEAN=',F7.2,' EARTH=',F7.2,'  OICE=',F7.2,'   LICE='
     +      ,F7.2,1X,' LTOPCL=',I2,' JLAT=',I2,' JYEAR=',I4)
 6305 FORMAT( 1X,'0.5',1X,F7.3,F7.2,2X,F7.3,F7.2,2X,'FR.TRULB1'
     +      ,' OCEAN=',F7.4,' EARTH=',F7.4,'  OICE=',F7.4,'   LICE='
     +      ,F7.4,1X,' LBOTCL=',I2,' ILON=',I2,'  JDAY=',I4)
 6306 FORMAT( 1X,'0.1',1X,F7.3,F7.2,2X,F7.3,F7.2,2X,'L=',I3,15I6)
 6307 FORMAT( 1X,'0.1',1X,F7.3,F7.2,2X,F7.3,F7.2,2X,'K=',I3,15I6)
 6308 FORMAT(' TKeff= ',F6.2,2F7.2,'  SRKALB=',16F6.4/
     +        1X,'At Top of Atm: ',' BTEMPW=',F6.2,1X,' TRUFTW=',F6.3
     +      , 2X,' SRIVIS=',F6.2,' SROVIS=',F6.2,   ' PLAVIS=',F6.4
     +      , 2X,' SRINIR=',F6.2,' SRONIR=',F6.2,   ' PLANIR=',F6.4)
 6309 FORMAT( 1X,'At Bot of Atm: ',' TRDFGW=',F6.3,1X,' TRUFGW=',F6.3
     +      , 2X,' SRDVIS=',F6.2,' SRUVIS=',F6.2,   ' ALBVIS=',F6.4
     +      , 2X,' SRDNIR=',F6.2,' SRUNIR=',F6.2,   ' ALBNIR=',F6.4)
 6310 FORMAT( 1X,'In Atmosphere: ',' SRXVIS=',F6.4,1X,' SRXNIR=',F6.4
     +      , 2X,' SRTVIS=',F6.4,' SRRVIS=',F6.4,   ' SRAVIS=',F6.4
     +      , 2X,' SRTNIR=',F6.4,' SRRNIR=',F6.4,   ' SRANIR=',F6.4)
 6311 FORMAT(' ')
 6313 FORMAT(1X,I2,F9.5,F6.2,2F7.2,1X,F7.4,2F7.2,F7.4,1X,3F8.2,F7.4
     +      ,1X,F7.2,F7.4,1X,3F6.2,1X,F5.4)
      GO TO 9999
C
C-------------
  400 CONTINUE
C-------------
C
C                                (4A)  Total Aerosol Qx, Qs, g, Pi0
C                                ----------------------------------
      NPAGE=1
      IF(INDEX.LT.11) NPAGE=KPAGE
      WRITE(KW,6401)
      DO 402 K=1,6
      SUM1(K)=0.
      SUM2(K)=0.
      SUM3(K)=0.
      DO 401 L=1,NL
      SUM1(K)=SUM1(K)+WSREXT(L,K)
      SUM2(K)=SUM2(K)+WSRSCT(L,K)
      SUM3(K)=SUM3(K)+WSRSCT(L,K)*WSRGCB(L,K)
  401 CONTINUE
      SUM3(K)=SUM3(K)/(SUM2(K)+1.D-10)
      SUM0(K)=SUM2(K)/(SUM1(K)+1.D-10)
  402 CONTINUE
      WRITE(KW,6402) (K,K=1,6),(K,K=1,6)
      DO 403 N=1,NL
      L=NLP-N
      WRITE(KW,6403) L,PLB(L),HLB0(L)
     +              ,(WSREXT(L,J),J=1,6),(WSRSCT(L,J),J=1,6)
  403 CONTINUE
      WRITE(KW,6404) (SUM1(K),K=1,6),(SUM2(K),K=1,6)
      NPAGE=0
      IF(NL.GT.13) NPAGE=1
      WRITE(KW,6405) KANORM
      WRITE(KW,6406) (K,K=1,6),(K,K=1,6)
      DO 404 N=1,NL
      L=NLP-N
      WRITE(KW,6407) L,PL(L),DPL(L)
     +              ,(WSRGCB(L,J),J=1,6),(WSRPI0(L,J),J=1,6)
  404 CONTINUE
      WRITE(KW,6408) (SUM3(K),K=1,6),(SUM0(K),K=1,6)
C     WRITE(KW,6420) (SRBALB(K),K=1,6)
C     WRITE(KW,6421) (SRXALB(K),K=1,6)
C     WRITE(KW,6422)
      SUM=0.
      DO 406 J=1,5
      TAU55=0.
      DO 405 I=1,11  !  NAERO
  405 TAU55=TAU55+AGOLDH(I,J)*FGOLDH(J)
      WRITE(KW,6423) J,FGOLDH(J),TAU55
  406 SUM=SUM+TAU55
      WRITE(KW,6424) SUM
      WRITE(KW,6424) BCOLX,ACOLX,DCOLX,VCOLX,TCOLX
      DO 407 I=1,8
      WRITE(KW,6425)
  407 CONTINUE
C
C                                (4B)  Water/Ice Cloud Qx, Qs, g, Pi0
C                                ------------------------------------
      NPAGE=1
      IF(INDEX.LT.11) NPAGE=KPAGE
      WRITE(KW,6411)
      DO 412 K=1,6
      SUM1(K)=0.
      SUM2(K)=0.
      SUM3(K)=0.
      DO 411 L=1,NL
      SUM1(K)=SUM1(K)+SRCEXT(L,K)
      SUM2(K)=SUM2(K)+SRCSCT(L,K)
      SUM3(K)=SUM3(K)+SRCSCT(L,K)*SRCGCB(L,K)
  411 SRCPI0(L,K)=SRCSCT(L,K)/(SRCEXT(L,K)+1.D-10)
      SUM3(K)=SUM3(K)/(SUM2(K)+1.D-10)
      SUM0(K)=SUM2(K)/(SUM1(K)+1.D-10)
  412 CONTINUE
      WRITE(KW,6412) (K,K=1,6),(K,K=1,6)
      DO 413 N=1,NL
      L=NLP-N
      WRITE(KW,6413) L,PLB(L),HLB0(L)
     +              ,(SRCEXT(L,J),J=1,6),(SRCSCT(L,J),J=1,6)
  413 CONTINUE
      WRITE(KW,6414) (SUM1(K),K=1,6),(SUM2(K),K=1,6)
      NPAGE=0
      IF(NL.GT.13) NPAGE=1
      WRITE(KW,6415) KANORM
      WRITE(KW,6416) (K,K=1,6),(K,K=1,6)
      DO 414 N=1,NL
      L=NLP-N
      WRITE(KW,6417) L,PL(L),DPL(L)
     +              ,(SRCGCB(L,J),J=1,6),(SRCPI0(L,J),J=1,6)
  414 CONTINUE
      WRITE(KW,6418) (SUM3(K),K=1,6),(SUM0(K),K=1,6)
      WRITE(KW,6420) (SRBALB(K),K=1,6)
      WRITE(KW,6421) (SRXALB(K),K=1,6)
      WRITE(KW,6422)
      SUM=0.
      DO 416 J=1,5
      TAU55=0.
      DO 415 I=1,11  !  NAERO
  415 TAU55=TAU55+AGOLDH(I,J)*FGOLDH(J)
      WRITE(KW,6423) J,FGOLDH(J),TAU55
  416 SUM=SUM+TAU55
      WRITE(KW,6424) SUM
      DO 417 I=1,2
      WRITE(KW,6425)
  417 CONTINUE
C
C                                (4C)  Aerosol + Cloud Qx, Qs, g, Pi0
C                                ------------------------------------
      NPAGE=1
      IF(INDEX.LT.11) NPAGE=KPAGE
      WRITE(KW,6426)
      DO 419 K=1,6
      SUM1(K)=0.
      SUM2(K)=0.
      SUM3(K)=0.
      DO 418 L=1,NL
      SUM1(K)=SUM1(K)+DBLEXT(L,K)
      SUM2(K)=SUM2(K)+DBLSCT(L,K)
      SUM3(K)=SUM3(K)+DBLSCT(L,K)*DBLGCB(L,K)
  418 DBLPI0(L,K)=DBLSCT(L,K)/(DBLEXT(L,K)+1.E-10)
      SUM3(K)=SUM3(K)/(SUM2(K)+1.E-10)
      SUM0(K)=SUM2(K)/(SUM1(K)+1.E-10)
  419 CONTINUE
      WRITE(KW,6427) (K,K=1,6),(K,K=1,6)
      DO 420 N=1,NL
      L=NLP-N
      WRITE(KW,6428) L,PLB(L),HLB0(L)
     +              ,(DBLEXT(L,J),J=1,6),(DBLSCT(L,J),J=1,6)
  420 CONTINUE
      WRITE(KW,6429) (SUM1(K),K=1,6),(SUM2(K),K=1,6)
      NPAGE=0
      IF(NL.GT.13) NPAGE=1
      WRITE(KW,6430) KANORM
      WRITE(KW,6431) (K,K=1,6),(K,K=1,6)
      DO 421 N=1,NL
      L=NLP-N
      WRITE(KW,6432) L,PL(L),DPL(L)
     +              ,(DBLGCB(L,J),J=1,6),(DBLPI0(L,J),J=1,6)
  421 CONTINUE
      WRITE(KW,6433) (SUM3(K),K=1,6),(SUM0(K),K=1,6)
      WRITE(KW,6434) (SRBALB(K),K=1,6)
      WRITE(KW,6435) (SRXALB(K),K=1,6)
      WRITE(KW,6436)
      SUM=0.
      DO 423 J=1,5
      TAU55=0.
      DO 422 I=1,11  !  NAERO
  422 TAU55=TAU55+AGOLDH(I,J)*FGOLDH(J)
      WRITE(KW,6437) J,FGOLDH(J),TAU55
  423 SUM=SUM+TAU55
      WRITE(KW,6438) SUM
      DO 424 I=1,2
      WRITE(KW,6439)
  424 CONTINUE
C
C                                (4D)  11-Comp Aerosol Qx, Qs, g, Pi0
C                                ------------------------------------
      NPAGE=1
      IF(INDEX.LT.11) NPAGE=KPAGE
      WRITE(KW,6440)  KWTRAB,(N,N=1,11)
      WRITE(KW,6441)
      DO 425 K=1,6
      WRITE(KW,6442) K,(SRAQEX(K,N),N=1,11)
  425 CONTINUE
      WRITE(KW,6443)
      DO 426 K=1,6
      WRITE(KW,6442) K,(SRAQSC(K,N),N=1,11)
  426 CONTINUE
      WRITE(KW,6444)
      DO 427 K=1,6
      WRITE(KW,6442) K,(SRAQCB(K,N),N=1,11)
  427 CONTINUE
      WRITE(KW,6445) TRABCD(1),TRAXSG(KWTRAB+1)
      DO 428 K=1,33
      IF(KWTRAB.EQ.0) WRITE(KW,6442) K,(TRAQAB(K,N),N=1,11)
      IF(KWTRAB.EQ.1) WRITE(KW,6442) K,(TRAQEX(K,N),N=1,11)
      IF(KWTRAB.EQ.2) WRITE(KW,6442) K,(TRAQSC(K,N),N=1,11)
      IF(KWTRAB.EQ.3) WRITE(KW,6442) K,(TRAQCB(K,N),N=1,11)
      IF(KWTRAB.EQ.4) THEN
      DO N=1,11
      TRPI0K(N)=TRAQSC(K,N)/(1.D-10+TRAQEX(K,N))
      END DO
      WRITE(KW,6442) K,(TRPI0K(N),N=1,11)
      ENDIF
  428 CONTINUE
      DO 429 I=1,1
      WRITE(KW,6446)
  429 CONTINUE
C
C
C                                (4E)  10-Comp Aerosol Qx, Qs, g, Pi0
C                                ------------------------------------
      WRITE(KW,6450) KWTRAB,(N,N=1,10),(REAERO(N),N=1,10)
      WRITE(KW,6451)
      DO 435 K=1,6
      WRITE(KW,6452) K,(SRBQEX(K,N),N=1,10)
  435 CONTINUE
      WRITE(KW,6453)
      DO 436 K=1,6
      WRITE(KW,6452) K,(SRBQSC(K,N),N=1,10)
  436 CONTINUE
      WRITE(KW,6454)
      DO 437 K=1,6
      WRITE(KW,6452) K,(SRBQCB(K,N),N=1,10)
  437 CONTINUE
      WRITE(KW,6455) TRABCD(2),TRAXSG(KWTRAB+1)
      DO 438 K=1,33
      IF(KWTRAB.EQ.0) WRITE(KW,6442) K,(TRBQAB(K,N),N=1,10)
      IF(KWTRAB.EQ.1) WRITE(KW,6442) K,(TRBQEX(K,N),N=1,10)
      IF(KWTRAB.EQ.2) WRITE(KW,6442) K,(TRBQSC(K,N),N=1,10)
      IF(KWTRAB.EQ.3) WRITE(KW,6442) K,(TRBQCB(K,N),N=1,10)
      IF(KWTRAB.EQ.4) THEN
      DO N=1,10
      TRPI0K(N)=TRBQSC(K,N)/(1.D-10+TRBQEX(K,N))
      END DO
      WRITE(KW,6442) K,(TRPI0K(N),N=1,10)
      ENDIF
  438 CONTINUE
      DO 439 I=1,1
      WRITE(KW,6456)
  439 CONTINUE
C
C                             (4F  8-size Dust Aerosol Qx, Qs, g, Pi0
C                             ---------------------------------------
      NPAGE=1
      IF(INDEX.LT.11) NPAGE=KPAGE
      WRITE(KW,6460) KWTRAB,(N,N=1,8),(REDUST(N),N=1,8)
      WRITE(KW,6461)
      DO 445 K=1,6
      WRITE(KW,6462) K,(SRAQEX(K,N),N=1,8)
  445 CONTINUE
      WRITE(KW,6463)
      DO 446 K=1,6
      WRITE(KW,6462) K,(SRAQSC(K,N),N=1,8)
  446 CONTINUE
      WRITE(KW,6464)
      DO 447 K=1,6
      WRITE(KW,6462) K,(SRAQCB(K,N),N=1,8)
  447 CONTINUE
      WRITE(KW,6465) TRABCD(4),TRAXSG(KWTRAB+1)
      DO 448 K=1,33
      IF(KWTRAB.EQ.0) WRITE(KW,6442) K,(TRDQAB(K,N),N=1,8)
      IF(KWTRAB.EQ.1) WRITE(KW,6442) K,(TRDQEX(K,N),N=1,8)
      IF(KWTRAB.EQ.2) WRITE(KW,6442) K,(TRDQSC(K,N),N=1,8)
      IF(KWTRAB.EQ.3) WRITE(KW,6442) K,(TRDQCB(K,N),N=1,8)
      IF(KWTRAB.EQ.4) THEN
      DO N=1,8
      TRPI0K(N)=TRDQSC(K,N)/(1.D-10+TRDQEX(K,N))
      END DO
      WRITE(KW,6442) K,(TRPI0K(N),N=1,8)
      ENDIF
  448 CONTINUE
      DO 449 I=1,1
      WRITE(KW,6466)
  449 CONTINUE
C
C                           (4G  15-Size/phase Cloud Qx, Qs, g, Pi0
C                           ---------------------------------------
      NPAGE=1
      IF(INDEX.LT.11) NPAGE=KPAGE
      WRITE(KW,6470) KWTRAB,(N,N=1,15)
      WRITE(KW,6471)
      DO 430 K=1,6
      WRITE(KW,6472) K,(SRCQEX(K,N),N=1,15)
  430 CONTINUE
      WRITE(KW,6473)
      DO 431 K=1,6
      WRITE(KW,6472) K,(SRCQSC(K,N),N=1,15)
  431 CONTINUE
      WRITE(KW,6474)
      DO 432 K=1,6
      WRITE(KW,6472) K,(SRCQCB(K,N),N=1,15)
  432 CONTINUE
      WRITE(KW,6475) TRABCD(3),TRAXSG(KWTRAB+1)
      DO 433 K=1,33
      IF(KWTRAB.EQ.0) WRITE(KW,6472) K,(TRCQAB(K,N),N=1,15)
      IF(KWTRAB.EQ.1) WRITE(KW,6472) K,(TRCQEX(K,N),N=1,15)
      IF(KWTRAB.EQ.2) WRITE(KW,6472) K,(TRCQSC(K,N),N=1,15)
      IF(KWTRAB.EQ.3) WRITE(KW,6472) K,(TRCQCB(K,N),N=1,15)
      IF(KWTRAB.EQ.4) THEN
      DO N=1,15
      TRPI0K(N)=TRCQSC(K,N)/(1.D-10+TRCQEX(K,N))
      END DO
      WRITE(KW,6442) K,(TRPI0K(N),N=1,15)
      ENDIF
  433 CONTINUE
      DO 434 I=1,2
      WRITE(KW,6476)
  434 CONTINUE
C
 6401 FORMAT(' (4A) Aerosol Input for Solar Radiation:'
     +      ,' Aerosol Radiative Parameters'
     +      ,T81,'LIST: SRAEXT(L,K),SRASCT(L,K),SRAGCB(L,K),SRAPI0(L,K)'
     +      //T42,'TAU -- EXTINCTION',T99,'TAU -- SCATTERING'
     +      ,/T24,53('-'),4X,53('-'))
 6402 FORMAT(' LN    PLB     HLB     K=',I3,5I9,7X,'K=',I3,5I9)
 6403 FORMAT(1X,I2,2F8.3,3X,6F9.6,3X,6F9.6)
 6404 FORMAT(/1X,T7,'COLUMN AMOUNT=',2X,6F9.6,3X,6F9.6)
 6405 FORMAT(6X,'KANORM=',1I1/T48,'COSBAR',T105,'PIZERO'
     +      ,/T24,53('-'),4X,53('-'))
 6406 FORMAT(' LN     PL     DPL     K=',I3,5I9,7X,'K=',I3,5I9)
 6407 FORMAT(1X,I2,2F8.3,3X,6F9.6,3X,6F9.6)
 6408 FORMAT(/1X,T7,'COLUMN   MEAN=',2X,6F9.6,3X,6F9.6)
C
 6411 FORMAT(' (4B) Cloud Input for Solar Radiation:'
     +      ,'   Cloud Radiative Parameters'
     +      ,T81,'LIST: SRCEXT(L,K),SRCSCT(L,K),SRCGCB(L,K),SRCPI0(L,K)'
     +      //T42,'TAU -- EXTINCTION',T99,'TAU -- SCATTERING'
     +      ,/T24,53('-'),4X,53('-'))
 6412 FORMAT(' LN    PLB     HLB     K=',I3,5I9,7X,'K=',I3,5I9)
 6413 FORMAT(1X,I2,2F8.3,3X,6F9.6,3X,6F9.6)
 6414 FORMAT(/1X,T7,'COLUMN AMOUNT=',2X,6F9.6,3X,6F9.6)
 6415 FORMAT(6X,'KANORM=',1I1/T48,'COSBAR',T105,'PIZERO'
     +      ,/T24,53('-'),4X,53('-'))
 6416 FORMAT(' LN     PL     DPL     K=',I3,5I9,7X,'K=',I3,5I9)
 6417 FORMAT(1X,I2,2F8.3,3X,6F9.6,3X,6F9.6)
 6418 FORMAT(/1X,T7,'COLUMN   MEAN=',2X,6F9.6,3X,6F9.6)
C
 6420 FORMAT(/1X,T7,'ALBEDO RSURFB=',2X,6F9.6,3X,6F9.6)
 6421 FORMAT( 1X,T7,'ALBEDO RSURFX=',2X,6F9.6,3X,6F9.6)
 6422 FORMAT(///T44,'AEROSOL COMPOSITION AND TYPE MIX:'
     +      ,T81,'FACTOR',6X,'VALUE',T107,'TAU(0.55)'/)
 6423 FORMAT(T81,'FGOLDH(',I1,') =',1P,E9.2,5X,0P,F7.4)
 6424 FORMAT(/T11,'SUM COLUMN TAU(0.55) =    BkGrnd    ClimAer   D Dust'
     +       ,'    VolAer    TotAer'/T33,5F10.5)
 6425 FORMAT(' ')
C
 6426 FORMAT(' (4C) Cloud+Aerosol  Output from SOLARM/SGPGXG:'
     +      ,'   Cloud+Aerosol Rad Parameters'
     +      ,T81,'LIST: DBLEXT(L,K),DBLSCT(L,K),DBLGCB(L,K),DBLPI0(L,K)'
     +      //T42,'TAU -- EXTINCTION',T99,'TAU -- SCATTERING'
     +      ,/T24,53('-'),4X,53('-'))
 6427 FORMAT(' LN    PLB     HLB     K=',I3,5I9,7X,'K=',I3,5I9)
 6428 FORMAT(1X,I2,2F8.3,3X,6F9.6,3X,6F9.6)
 6429 FORMAT(/1X,T7,'COLUMN AMOUNT=',2X,6F9.6,3X,6F9.6)
 6430 FORMAT(6X,'KANORM=',1I1/T48,'COSBAR',T105,'PIZERO'
     +      ,/T24,53('-'),4X,53('-'))
 6431 FORMAT(' LN     PL     DPL     K=',I3,5I9,7X,'K=',I3,5I9)
 6432 FORMAT(1X,I2,2F8.3,3X,6F9.6,3X,6F9.6)
 6433 FORMAT(/1X,T7,'COLUMN   MEAN=',2X,6F9.6,3X,6F9.6)
 6434 FORMAT(/1X,T7,'ALBEDO RSURFB=',2X,6F9.6,3X,6F9.6)
 6435 FORMAT( 1X,T7,'ALBEDO RSURFX=',2X,6F9.6,3X,6F9.6)
 6436 FORMAT(///T44,'AEROSOL COMPOSITION AND TYPE MIX:'
     +      ,T81,'FACTOR',6X,'VALUE',T107,'TAU(0.55)'/)
 6437 FORMAT(T81,'FGOLDH(',I1,') =',1P,E9.2,5X,0P,F7.4)
 6438 FORMAT(/T81,'SUM COLUMN TAU(0.55) =',F10.4)
 6439 FORMAT(' ')
C
 6440 FORMAT(' (4D) Background Aerosol Solar and Thermal Mie '
     +      ,'Scattering Parameters:'
     +      ,T81,'List: SRAQEX(L,K),SRAQST(L,K),SRAQCB(L,K), TRAB Q S G'
     +      /'      KWTRAB=',I1/7X,11I8/
     +        '   AEROSOL  ACID1   SSALT   SLFT1   SLFT2   BSLT1'
     +        ,'   BSLT2   DUST1   DUST2   DUST3   CARB1   CARB2'/
     +        '   SIZE      0.5     2.0     0.3     1.0     0.5 '
     +        ,'    2.0     0.5     2.0     8.0     0.1     0.5 ')
 6441 FORMAT('  K  SRAQEX')
 6442 FORMAT(I3,6X,15F8.5)
 6443 FORMAT('  K  SRAQSC')
 6444 FORMAT('  K  SRAQCB')
 6445 FORMAT('  K  ',2A3)
 6446 FORMAT(' ')
C
 6450 FORMAT(' (4E) Climatology Aerosol Solar and Thermal Mie '
     +      ,'Scattering Parameters:'
     +      ,T81,'List: SRBQEX(L,K),SRBQST(L,K),SRBQCB(L,K), TRAB Q S G'
     +      /'      KWTRAB=',I1/7X,10I8/
     +        '   AEROSOL   BCI     OCI     SUI     SEA     SUN '
     +        ,'    ANT     OCN     OCB     BCB     SSB         '/
     +        '   SIZE ',10F8.1)
 6451 FORMAT('  K  SRAQEX')
 6452 FORMAT(I3,6X,15F8.5)
 6453 FORMAT('  K  SRAQSC')
 6454 FORMAT('  K  SRAQCB')
 6455 FORMAT('  K  ',2A3)
 6456 FORMAT(' ')
C
 6460 FORMAT(' (4F) Desert Dust Aerosol Solar and Thermal Mie '
     +      ,'Scattering Parameters:'
     +      ,T81,'List: SRDQEX(L,K),SRDQST(L,K),SRDQCB(L,K), TRAB Q S G'
     +      /'      KWTRAB=',I1/7X,8I8/
     +        '   AEROSOL  CLAY1   CLAY2   CLAY3   CLAY4   SILT1'
     +        ,'   SILT2   SILT3   SILT4                        '/
     +        '   SIZE ',8F8.1)
 6461 FORMAT('  K  SRAQEX')
 6462 FORMAT(I3,6X,15F8.5)
 6463 FORMAT('  K  SRAQSC')
 6464 FORMAT('  K  SRAQCB')
 6465 FORMAT('  K  ',2A3)
 6466 FORMAT(' ')
C
 6470 FORMAT(' (4G) Cloud Input for Solar, Thermal Radiation:'
     +      ,'  Mie Cloud Radiative Properties'
     +      ,T81,'List: SRCQEX(L,K),SRCQST(L,K),SRCQCB(L,K), TRAB Q S G'
     +      /'      KWTRAB=',I1/7X,15I8/
     +        ' WIM CLOUD  WAT05   WAT10   WAT15   WAT20   WAT25'
     +                ,'   ICE05   ICE15   ICE25   ICE50   ICE75'
     +                ,'   MIC05   MIC15   MIC25   MIC50   MIC75')
 6471 FORMAT('  K  SRCQEX')
 6472 FORMAT(I3,6X,15F8.5)
 6473 FORMAT('  K SRCQSC')
 6474 FORMAT('  K SRCQCB')
 6475 FORMAT('  K ',2A3)
 6476 FORMAT(' ')
      GO TO 9999
C
C-------------
  500 CONTINUE
C-------------
C
      NPAGE=1
      IF(INDEX.LT.11) NPAGE=KPAGE
      SIGMA=5.6697D-08
      TGMEAN=POCEAN*TGO**4+PEARTH*TGE**4+PLICE*TGLI**4+POICE*TGOI**4
      TGMEAN=SQRT(TGMEAN)
      TGMEAN=SQRT(TGMEAN)
      SIGT4=SIGMA*TGMEAN**4
      ITG=TGMEAN
      WTG=TGMEAN-ITG
      ITG=ITG-ITPFT0
      SUMK=0.0
      DO 501 K=1,33
      BGFLUX(K)=PLANCK(ITG)-(PLANCK(ITG)-PLANCK(ITG+1))*WTG
      BGFRAC(K)=BGFLUX(K)/SIGT4
      SUMK=SUMK+BGFLUX(K)
      ITG=ITG+ITNEXT
  501 CONTINUE
      LK=0
      DO 503 K=1,33
      TAUSUM(K)=0. !! TAUSL(K)
      DO 502 L=1,NL
      TRTAUK(L,K)=TRGXLK(L,K)+TRCALK(L,K)+TRAALK(L,K)
      TAUSUM(K)=TAUSUM(K)+TRGXLK(L,K)+TRCALK(L,K)+TRAALK(L,K)
  502 CONTINUE
  503 CONTINUE
      WRITE(KW,6501)
      WRITE(KW,6502) (K,K=1,13)
      DO 504 N=1,NL
      L=NLP-N
      WRITE(KW,6503) L,PL(L),TLM(L),(TRTAUK(L,K),K=1,13)
  504 CONTINUE
!!    WRITE(KW,6504) (TAUSL(K),K=1,13)
      WRITE(KW,6505) (TAUSUM(K),K=1,13)
      WRITE(KW,6506) SUMK,(BGFLUX(K),K=1,13)
      WRITE(KW,6507) TGMEAN,SIGT4,(BGFRAC(K),K=1,13)
      NPAGE=0
      IF(NL.GT.13)  NPAGE=1
      WRITE(KW,6508) NPAGE
      WRITE(KW,6509) (K,K=14,33)
      DO 505 N=1,NL
      L=NLP-N
      WRITE(KW,6510) L,(TRTAUK(L,K),K=14,33)
  505 CONTINUE
!!    WRITE(KW,6511) ( TAUSL(K),K=14,33)
      WRITE(KW,6512) (TAUSUM(K),K=14,33)
      WRITE(KW,6513) (BGFLUX(K),K=14,33)
      WRITE(KW,6514) (BGFRAC(K),K=14,33)
      DO 506 I=1,10
      WRITE(KW,6515)
  506 CONTINUE
C
 6501 FORMAT(' (5) TAU TABLE FOR THERMAL RADIATION: CONTAINS'
     +      ,' TOTAL SPECIFIED GAS, CLOUD & AEROSOL ABSORPTION'
     +      ,T99,'TRGXLK(L,K),TRCALK(L,K),TRCAAK(L,K)'/
     +      ,/1X,'K-DIST BREAKDOWN:',T23,'WINDOW'
     +      ,3X,'WATER VAPOR:',T71,'PRINCIPAL ABSORBER REGION'
     +      ,/T23,6('-'),3X,101('-'))
 6502 FORMAT(' LN     PL     TLM      K=',I1,4X,'K=',I2,9I9,3I8)
 6503 FORMAT(1X,I2,F8.3,F7.2,1X,10F9.4,3F8.3)
!!6504 FORMAT(/4X,'SURFACE LAYER= ',10F9.4,3F8.3)
 6505 FORMAT( 4X,'COLUMN AMOUNT= ',10F9.4,3F8.3)
 6506 FORMAT(/1X,'PF W/M**2= '  ,F6.2,1X,10F9.3,3F8.3)
 6507 FORMAT( 1X,'TG=',F6.2,'= ',F6.2,1X,10F9.4,3F8.3)
 6508 FORMAT(1I1/4X,'CARBON DIOXIDE:',T36,'PRINCIPAL ABSORBER REGION'
     +      ,T83,'OZONE:',T100,'PRINCIPAL ABSORBER REGION'
     +      /4X,76('-'),2X,50('-'))
 6509 FORMAT(1X,'LN  K=',I2,5I7,6I6,3X,'K=',I2,3I7,6I6)
 6510 FORMAT( 1X,  I2,6F7.4,2F6.3,3F6.2,1F6.1,4F7.4,3F6.3,F6.2)
!! 6511 FORMAT(/1X,'SL',6F7.4,2F6.3,3F6.2,1F6.1,4F7.4,3F6.3,F6.2)
 6512 FORMAT( 1X,'CA',5F7.4,1F7.3,3F6.2,2F6.1,1F6.0,4F7.4,2F6.3,2F6.2)
 6513 FORMAT(/1X,'PF',1F7.4,5F7.3,1F6.2,3F6.3,2F6.3,2F7.3,2F7.4,4F6.3)
 6514 FORMAT( 1X,'FR',6F7.4,2F6.3,3F6.3,1F6.3,4F7.4,3F6.3,F6.3)
 6515 FORMAT(' ')
      GO TO 9999
C
C-------------
  600 CONTINUE
C-------------
C
      NPAGE=1
      IF(INDEX.LT.11) NPAGE=KPAGE
      SIGMA=5.6697D-08
      TGMEAN=POCEAN*TGO**4+PEARTH*TGE**4+PLICE*TGLI**4+POICE*TGOI**4
      TGMEAN=SQRT(TGMEAN)
      TGMEAN=SQRT(TGMEAN)
      SIGT4=SIGMA*TGMEAN**4
      ITG=TGMEAN
      WTG=TGMEAN-ITG
      ITG=ITG-ITPFT0
      SUMK=0.0
      DO 601 K=1,33
      BGFLUX(K)=PLANCK(ITG)-(PLANCK(ITG)-PLANCK(ITG+1))*WTG
      BGFRAC(K)=BGFLUX(K)/SIGT4
      SUMK=SUMK+BGFLUX(K)
      ITG=ITG+ITNEXT
  601 CONTINUE
      WRITE(KW,6601)
      WRITE(KW,6602) (K,K=1,13)
      DO 602 N=1,NL
      L=NLP-N
      WRITE(KW,6603) L,PL(L),TLM(L),(TRGXLK(L,K),K=1,13)
  602 CONTINUE
      LK=0
      DO 604 K=1,33
      TAUSUM(K)=0. !!  TAUSL(K)
      DO 603 L=1,NL
  603 TAUSUM(K)=TAUSUM(K)+TRGXLK(L,K)
  604 CONTINUE
!!    WRITE(KW,6604) (TAUSL(K),K=1,13)
      WRITE(KW,6605) (TAUSUM(K),K=1,13)
      WRITE(KW,6606) SUMK,(BGFLUX(K),K=1,13)
      WRITE(KW,6607) TGMEAN,SIGT4,(BGFRAC(K),K=1,13)
      NPAGE=0
      IF(NL.GT.13)  NPAGE=1
      WRITE(KW,6608)
      WRITE(KW,6609) (K,K=14,33)
      DO 605 N=1,NL
      L=NLP-N
      WRITE(KW,6610) L,(TRGXLK(L,K),K=14,33)
  605 CONTINUE
!!    WRITE(KW,6611) ( TAUSL(K),K=14,33)
      WRITE(KW,6612) (TAUSUM(K),K=14,33)
      WRITE(KW,6613) (BGFLUX(K),K=14,33)
      WRITE(KW,6614) (BGFRAC(K),K=14,33)
      DO 606 I=1,10
      WRITE(KW,6615)
  606 CONTINUE
C
 6601 FORMAT(' (6) TAU TABLE FOR THERMAL RADIATION: INCLUDES ANY'
     +      ,' SPECIFIED OVERLAP, CLOUD & AEROSOL ABSORPTION'
     +      ,T114,'TRGXLK(L,K),TAUSL(L)'/
     +      ,/1X,'K-DIST BREAKDOWN:',T23,'WINDOW'
     +      ,3X,'WATER VAPOR:',T71,'PRINCIPAL ABSORBER REGION'
     +      ,/T23,6('-'),3X,101('-'))
 6602 FORMAT(' LN     PL     TLM      K=',I1,4X,'K=',I2,9I9,3I8)
 6603 FORMAT(1X,I2,F8.3,F7.2,1X,10F9.4,3F8.3)
 6604 FORMAT(/4X,'SURFACE LAYER= ',10F9.4,3F8.3)
 6605 FORMAT( 4X,'COLUMN AMOUNT= ',10F9.4,3F8.3)
 6606 FORMAT(/1X,'PF W/M**2= '  ,F6.2,1X,10F9.3,3F8.3)
 6607 FORMAT( 1X,'TG=',F6.2,'= ',F6.2,1X,10F9.4,3F8.3)
 6608 FORMAT(/4X,'CARBON DIOXIDE:',T36,'PRINCIPAL ABSORBER REGION'
     +      ,T83,'OZONE:',T100,'PRINCIPAL ABSORBER REGION'
     +      /4X,76('-'),2X,50('-'))
 6609 FORMAT(1X,'LN  K=',I2,5I7,6I6,3X,'K=',I2,3I7,6I6)
 6610 FORMAT( 1X,  I2,6F7.4,2F6.3,3F6.2,1F6.1,4F7.4,3F6.3,F6.2)
!! 6611 FORMAT(/1X,'SL',6F7.4,2F6.3,3F6.2,1F6.1,4F7.4,3F6.3,F6.2)
 6612 FORMAT( 1X,'CA',5F7.4,1F7.3,3F6.2,2F6.1,1F6.0,4F7.4,2F6.3,2F6.2)
 6613 FORMAT(/1X,'PF',1F7.4,5F7.3,1F6.2,3F6.3,2F6.3,2F7.3,2F7.4,4F6.3)
 6614 FORMAT( 1X,'FR',6F7.4,2F6.3,3F6.3,1F6.3,4F7.4,3F6.3,F6.3)
 6615 FORMAT(' ')
      GO TO 9999
C
C-------------
  700 CONTINUE
C-------------
C
      SIGMA=5.6697D-08
      TGMEAN=POCEAN*TGO**4+PEARTH*TGE**4+PLICE*TGLI**4+POICE*TGOI**4
      TGMEAN=SQRT(TGMEAN)
      TGMEAN=SQRT(TGMEAN)
      SIGT4=SIGMA*TGMEAN**4
      ITG=TGMEAN
      WTG=TGMEAN-ITG
      ITG=ITG-ITPFT0
      SUMK=0.0
      DO 701 K=1,33
      BGFLUX(K)=PLANCK(ITG)-(PLANCK(ITG)-PLANCK(ITG+1))*WTG
      BGFRAC(K)=BGFLUX(K)/SIGT4
      SUMK=SUMK+BGFLUX(K)
      ITG=ITG+ITNEXT
  701 CONTINUE
      WRITE(KW,6701)
      WRITE(KW,6702) (K,K=1,13)
      DO 702 N=1,NL
      L=NLP-N
      WRITE(KW,6703) L,PL(L),TLM(L),(TRCALK(L,K),K=1,13)
  702 CONTINUE
      LK=0
      DO 704 K=1,33
      TAUSUM(K)=0.0
      DO 703 L=1,NL
      LK=LK+1
  703 TAUSUM(K)=TAUSUM(K)+TRCALK(L,K)
  704 CONTINUE
      WRITE(KW,6704) (TAUSUM(K),K=1,13),(TRCTCA(K),K=1,13)
      WRITE(KW,6705)
      WRITE(KW,6706)  SUMK,(BGFLUX(K),K=1,13)
      WRITE(KW,6707) TGMEAN,SIGT4,(BGFRAC(K),K=1,13)
C
      WRITE(KW,6708)
      WRITE(KW,6709) (K,K=14,33)
      DO 705 N=1,NL
      L=NLP-N
      WRITE(KW,6710) L,(TRCALK(L,K),K=14,33)
  705 CONTINUE
      WRITE(KW,6711) (TAUSUM(K),K=14,33),(TRCTCA(K),K=14,33)
      WRITE(KW,6712) (BGFLUX(K),K=14,33)
      WRITE(KW,6713) (BGFRAC(K),K=14,33)
      DO 706 I=1,8
      WRITE(KW,6714)
  706 CONTINUE
C
      SIGMA=5.6697D-08
      TGMEAN=POCEAN*TGO**4+PEARTH*TGE**4+PLICE*TGLI**4+POICE*TGOI**4
      TGMEAN=SQRT(TGMEAN)
      TGMEAN=SQRT(TGMEAN)
      SIGT4=SIGMA*TGMEAN**4
      ITG=TGMEAN
      WTG=TGMEAN-ITG
      ITG=ITG-ITPFT0
      SUMK=0.0
      DO 711 K=1,33
      BGFLUX(K)=PLANCK(ITG)-(PLANCK(ITG)-PLANCK(ITG+1))*WTG
      BGFRAC(K)=BGFLUX(K)/SIGT4
      SUMK=SUMK+BGFLUX(K)
      ITG=ITG+ITNEXT
  711 CONTINUE
      WRITE(KW,6721)
      WRITE(KW,6722) (K,K=1,13)
      DO 712 N=1,NL
      L=NLP-N
      WRITE(KW,6723) L,PL(L),TLM(L),(TRAALK(L,K),K=1,13)
  712 CONTINUE
      DO 714 K=1,33
      TAUSUM(K)=0.0
      DO 713 L=1,NL
  713 TAUSUM(K)=TAUSUM(K)+TRAALK(L,K)
  714 CONTINUE
      WRITE(KW,6724) (TAUSUM(K),K=1,13)
      WRITE(KW,6725)
      WRITE(KW,6726)         SUMK,(BGFLUX(K),K=1,13)
      WRITE(KW,6727) TGMEAN,SIGT4,(BGFRAC(K),K=1,13)
      NPAGE=0
      IF(NL.GT.13)  NPAGE=1
      WRITE(KW,6728) NPAGE
      WRITE(KW,6729) (K,K=14,33)
      DO 715 N=1,NL
      L=NLP-N
      WRITE(KW,6730) L,(TRAALK(L,K),K=14,33)
  715 CONTINUE
      WRITE(KW,6731) (TAUSUM(K),K=14,33)
      WRITE(KW,6732) (BGFLUX(K),K=14,33)
      WRITE(KW,6733) (BGFRAC(K),K=14,33)
      DO 716 I=1,12
      WRITE(KW,6734)
  716 CONTINUE
C
 6701 FORMAT(' (7A) TRCALK TABLE FOR THERMAL RADIATION: CONTAINS'
     +      ,' 33 KD CLOUD ABSORPTION OPTICAL DEPTHS AT'
     +      ,' THERMAL WAVELENGTHS ',T117,'LIST: TRCALK(L,K)'/
     +      ,/1X,'K-DIST BREAKDOWN:',T23,'WINDOW'
     +      ,3X,'WATER VAPOR:',T71,'PRINCIPAL ABSORBER REGION'
     +      ,/T23,6('-'),3X,101('-'))
 6702 FORMAT(' LN     PL     TLM      K=',I1,6X,I2,9I9,3I8)
 6703 FORMAT(1X,I2,F8.3,F7.2,1X,9F9.5,4F8.5)
 6704 FORMAT(/4X,'COLUMN AMOUNT= ',9F9.4,4F8.5
     +       /4X,'TOPCLD ALBEDO= ',9F9.4,4F8.5)
 6705 FORMAT(/' K-INTERVAL CONTRIBUTIONS:'/' COMPARE WITH GROUND FLUX:')
 6706 FORMAT( 1X,'PF W/M**2= '  ,F6.2,1X,10F9.3,3F8.3)
 6707 FORMAT( 1X,'TG=',F6.2,'= ',F6.2,1X,10F9.4,3F8.3)
 6708 FORMAT(/T25,'CARBON DIOXIDE:   PRINCIPAL ABSORBER REGION'
     +      ,T93,'OZONE:   PRINCIPAL ABSORBER REGION'
     +      /4X,77('-'),1X,51('-'))
 6709 FORMAT(1X,'LN  K=',I2,5I7,6I6,3X,'K=',I2,3I7,6I6)
 6710 FORMAT( 1X,  I2,6F7.4,5F6.4,F6.4,4F7.4,3F6.4,F6.4)
 6711 FORMAT(/1X,'CA',5F7.4,1F7.3,3F6.2,2F6.1,1F6.0,4F7.4,2F6.3,2F6.2
     +       /1X,'TA',5F7.4,1F7.4,3F6.3,2F6.3,1F6.3,4F7.4,2F6.3,2F6.3)
 6712 FORMAT(/1X,'PF',1F7.4,5F7.3,1F6.2,3F6.3,2F6.3,2F7.3,2F7.4,4F6.3)
 6713 FORMAT( 1X,'FR',6F7.4,2F6.3,3F6.3,1F6.3,4F7.4,3F6.3,F6.3)
 6714 FORMAT(' ')
C
 6721 FORMAT(' (7B) AEROSOL TAU TABLE FOR THERMAL RADIATION:'
     +      ,'  AEROSOL ABSORPTION OPTICAL DEPTH AT THERMAL WAVELENGTHS'
     +      ,T116,'LIST:  TRAALK(L,K)'/
     +      ,/1X,'K-DIST BREAKDOWN:',T23,'WINDOW'
     +      ,3X,'WATER VAPOR:',T71,'PRINCIPAL ABSORBER REGION'
     +      ,/T23,6('-'),3X,101('-'))
 6722 FORMAT(' LN     PL     TLM      K=',I1,6X,I2,9I9,3I8)
 6723 FORMAT(1X,I2,F8.3,F7.2,1X,10F9.5,3F8.5)
 6724 FORMAT(/4X,'COLUMN AMOUNT= ',10F9.5,3F8.5)
 6725 FORMAT(' K-INTERVAL CONTRIBUTIONS:'/' COMPARE WITH GROUND FLUX:')
 6726 FORMAT( 1X,'PF W/M**2= '  ,F6.2,1X,10F9.3,3F8.3)
 6727 FORMAT( 1X,'TG=',F6.2,'= ',F6.2,1X,10F9.4,3F8.3)
 6728 FORMAT(1I1/4X,'CARBON DIOXIDE:',T36,'PRINCIPAL ABSORBER REGION'
     +      ,T83,'OZONE:',T100,'PRINCIPAL ABSORBER REGION'
     +      /4X,76('-'),2X,50('-'))
 6729 FORMAT(1X,'LN  K=',I2,5I7,6I6,3X,'K=',I2,3I7,6I6)
 6730 FORMAT( 1X,  I2,6F7.5,2F6.4,3F6.4,F6.4,4F7.4,3F6.4,F6.4)
 6731 FORMAT( 1X,'CA',5F7.5,1F7.5,3F6.4,2F6.4,1F6.4,4F7.4,2F6.4,2F6.4)
 6732 FORMAT(/1X,'PF',1F7.4,5F7.3,1F6.2,3F6.3,2F6.3,2F7.3,2F7.4,4F6.3)
 6733 FORMAT( 1X,'FR',6F7.4,2F6.3,3F6.3,1F6.3,4F7.4,3F6.3,F6.3)
 6734 FORMAT(' ')
      GO TO 9999
C
C-------------
  800 CONTINUE
C-------------
C
      WRITE(KW,6800)
      DO 801 K=1,16
      ISR1(K)=NORDER(K)
      IF(KORDER.EQ.1) ISR1(K)=K
  801 CONTINUE
      WRITE(KW,6801) (ISR1(K),K=1,16)
      SUMK=0.0
      DO 802 K=1,16
      FSR1(K)=SKDKS0(NORDER(K))
      IF(KORDER.EQ.1) FSR1(K)=SKDKS0(K)
      SUMK=SUMK+FSR1(K)
  802 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6802) (FSR1(K),K=1,17)
      DO 803 K=1,16
      ISR1(K)=NMWAVA(K)
      IF(KORDER.EQ.1) ISR1(K)=NMWAVA(IORDER(K))
  803 CONTINUE
      WRITE(KW,6803) (ISR1(K),K=1,16)
      DO 804 K=1,16
      ISR1(K)=NMWAVB(K)
      IF(KORDER.EQ.1) ISR1(K)=NMWAVB(IORDER(K))
  804 CONTINUE
      WRITE(KW,6804) (ISR1(K),K=1,16)
      IF(KORDER.EQ.0) WRITE(KW,6805)
      IF(KORDER.EQ.1) WRITE(KW,6806)
      DO 805 K=1,16
      ISR1(K)=K
      IF(KORDER.EQ.1)ISR1(K)=IORDER(K)
  805 CONTINUE
      WRITE(KW,6807) (ISR1(K),K=1,16)
      DO 807 N=1,NLP
      L=NLP+1-N
      SUMK=0.0
      DO 806 K=1,16
      FSR1(K)=SKDFLB(L,NORDER(K))
      IF(KORDER.EQ.1) FSR1(K)=SKDFLB(L,K)
      SUMK=SUMK+FSR1(K)
  806 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6808) L,(FSR1(K),K=1,17)
  807 CONTINUE
      DO 808 K=1,16
      ISR1(K)=K
      IF(KORDER.EQ.1) ISR1(K)=IORDER(K)
  808 CONTINUE
      WRITE(KW,6809) (ISR1(K),K=1,16)
      DO 810 N=1,NLP
      L=NLP+1-N
      SUMK=0.0
      DO 809 K=1,16
      FSR1(K)=SKUFLB(L,NORDER(K))
      IF(KORDER.EQ.1) FSR1(K)=SKUFLB(L,K)
      SUMK=SUMK+FSR1(K)
  809 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6810) L,(FSR1(K),K=1,17)
  810 CONTINUE
      DO 811 K=1,16
      ISR1(K)=K
      IF(KORDER.EQ.1) ISR1(K)=IORDER(K)
  811 CONTINUE
      WRITE(KW,6811) (ISR1(K),K=1,16)
      SUMT=0.D0
      SUMK=0.D0
      DO 812 K=1,16
      FSR1(K)=SRKALB(NORDER(K))
      FSR2(K)=SKDKS0(NORDER(K))
      IF(KORDER.EQ.1) FSR1(K)=SRKALB(K)
      IF(KORDER.EQ.1) FSR2(K)=SKDKS0(K)
      SUMK=SUMK+FSR1(K)*FSR2(K)
  812 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6812) (FSR1(K),K=1,17)
      SUMT=SUMT+FSR1(17)
      SUMK1=0.D0
      SUMK2=0.D0
      DO 813 K=1,16
      FSR1(K)=SKNFLB(NLP,NORDER(K))-SKNFLB(1,NORDER(K))
      FSR2(K)=SKDFLB(NLP,NORDER(K))
      IF(KORDER.EQ.1) FSR1(K)=SKNFLB(NLP,K)-SKNFLB(1,K)
      IF(KORDER.EQ.1) FSR2(K)=SKDFLB(NLP,K)
      SUMK1=SUMK1+FSR1(K)
      SUMK2=SUMK2+FSR2(K)
      FSR1(K)=FSR1(K)/(FSR2(K)+1.d-20)
  813 CONTINUE
      FSR1(17)=SUMK1/(SUMK2+1.d-20)
      WRITE(KW,6813) (FSR1(K),K=1,17)
      SUMT=SUMT+FSR1(17)
      SUMK1=0.D0
      SUMK2=0.D0
      DO 814 K=1,16
      FSR1(K)=SKNFLB(1,NORDER(K))
      FSR2(K)=SKDFLB(NLP,NORDER(K))
      IF(KORDER.EQ.1) FSR1(K)=SKNFLB(1,K)
      IF(KORDER.EQ.1) FSR2(K)=SKDFLB(NLP,K)
      SUMK1=SUMK1+FSR1(K)
      SUMK2=SUMK2+FSR2(K)
      FSR1(K)=FSR1(K)/(FSR2(K)+1.d-20)
  814 CONTINUE
      FSR1(17)=SUMK1/(SUMK2+1.d-20)
      WRITE(KW,6814) (FSR1(K),K=1,17)
      SUMT=SUMT+FSR1(17)
      DO 815 K=1,16
      ISR1(K)=KSLAMW(NORDER(K))
      IF(KORDER.EQ.1) ISR1(K)=KSLAMW(K)
  815 CONTINUE
      WRITE(KW,6815) SUMT,(ISR1(K),K=1,16)
      SUMK=0.D0
      DO 816 K=1,16
      KK=KSLAMW(NORDER(K))
      IF(KORDER.EQ.1) KK=KSLAMW(K)
      FSR1(K)=SRBXAL(KK,1)
      FSR2(K)=SRBXAL(KK,2)
  816 CONTINUE
      WRITE(KW,6816) (FSR1(K),K=1,16)
      WRITE(KW,6817) (FSR2(K),K=1,16)
      WRITE(KW,6818)   COSZ,SRIVIS,SROVIS,PLAVIS,SRINIR,SRONIR,PLANIR
      WRITE(KW,6819) SRXVIS,SRXNIR,SRDVIS,SRUVIS,ALBVIS,SRDNIR
     +                            ,SRUNIR,ALBNIR
      WRITE(KW,6820) SRTVIS,SRRVIS,SRAVIS,SRTNIR,SRRNIR,SRANIR
      DO 817 I=1,1
      IF(KORDER.EQ.1) WRITE(KW,6821)
  817 CONTINUE
C
      WRITE(KW,6840)
      DO 821 K=1,16
      ISR1(K)=NORDER(K)
      IF(KORDER.EQ.1) ISR1(K)=K
  821 CONTINUE
      WRITE(KW,6841) (ISR1(K),K=1,16)
      SUMK=0.0
      DO 822 K=1,16
      FSR1(K)=SKDKS0(NORDER(K))
      IF(KORDER.EQ.1) FSR1(K)=SKDKS0(K)
      SUMK=SUMK+FSR1(K)
  822 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6842) (FSR1(K),K=1,17)
      IF(KORDER.EQ.0) WRITE(KW,6843)
      IF(KORDER.EQ.1) WRITE(KW,6844)
      DO 825 K=1,16
      ISR1(K)=K
      IF(KORDER.EQ.1)ISR1(K)=IORDER(K)
  825 CONTINUE
      WRITE(KW,6845) (ISR1(K),K=1,16)
      DO 827 N=1,NLP
      L=NLP+1-N
      SUMK=0.0
      DO 826 K=1,16
      FSR1(K)=SKNFLB(L,NORDER(K))
      IF(KORDER.EQ.1) FSR1(K)=SKNFLB(L,K)
      SUMK=SUMK+FSR1(K)
  826 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6846) L,(FSR1(K),K=1,17)
  827 CONTINUE
      DO 828 K=1,16
      ISR1(K)=K
      IF(KORDER.EQ.1) ISR1(K)=IORDER(K)
  828 CONTINUE
      WRITE(KW,6847) (ISR1(K),K=1,16)
      DO 830 N=1,NL
      L=NL+1-N
      SUMK=0.0
      DO 829 K=1,16
      FSR1(K)=SKFHRL(L,NORDER(K))
      IF(KORDER.EQ.1) FSR1(K)=SKFHRL(L,K)
      SUMK=SUMK+FSR1(K)
  829 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6848) L,(FSR1(K),K=1,17)
  830 CONTINUE
      DO 831 K=1,16
      ISR1(K)=K
      IF(KORDER.EQ.1)ISR1(K)=IORDER(K)
  831 CONTINUE
      WRITE(KW,6849) (ISR1(K),K=1,16)
      DO 833 N=1,4
      SUMK=0.0
      DO 832 K=1,16
      FSR1(K)=SRKGAX(NORDER(K),N)
      IF(KORDER.EQ.1) FSR1(K)=SRKGAX(K,N)
      SUMK=SUMK+FSR1(K)
  832 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6850) L,(FSR1(K),K=1,17)
  833 CONTINUE
      DO 834 K=1,16
      ISR1(K)=K
      IF(KORDER.EQ.1)ISR1(K)=IORDER(K)
  834 CONTINUE
      WRITE(KW,6851)
      DO 836 N=1,4
      SUMK=0.0
      DO 835 K=1,16
      FSR1(K)=SRKGAD(NORDER(K),N)
      IF(KORDER.EQ.1) FSR1(K)=SRKGAD(K,N)
      SUMK=SUMK+FSR1(K)
  835 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6852) N,(FSR1(K),K=1,17)
  836 CONTINUE
      WRITE(KW,6853)
      DO 838 N=1,4
      SUMK=0.0
      DO 837 K=1,16
      FSR1(K)=SRKGAX(NORDER(K),N)
      FSR2(K)=SRKGAD(NORDER(K),N)
      IF(KORDER.EQ.1) FSR1(K)=SRKGAX(K,N)
      IF(KORDER.EQ.1) FSR2(K)=SRKGAD(K,N)
      FSR1(K)=FSR1(K)+FSR2(K)
      SUMK=SUMK+FSR1(K)
  837 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6854) N,(FSR1(K),K=1,17)
  838 CONTINUE
      WRITE(KW,6855) SRNFLB(1),POCEAN,FSRNFG(1),PEARTH,FSRNFG(2)
     +                        ,POICE ,FSRNFG(3),PLICE ,FSRNFG(4)
C
 6800 FORMAT(' (8A)     SPECTRAL/k-DISTRIBUTION COMPONENT BREAKDOWN'
     +      ,' FOR DOWNWARD AND UPWARD SOLAR RADIATIVE FLUXES'
     +      ,T108,'SKDFLB(L,K)  SKUFLB(L,K)  SRKALB(K)'/)
 6801 FORMAT('     K=',I5,2I8,2I7,I8,I9,7I8,I7,I8,'       Total')
 6802 FORMAT('  DKS0=',F6.3,2F8.3,2F7.3,F8.3,F9.3,7F8.3,F7.3,F8.3,F11.3)
 6803 FORMAT('NMWAVA=',I6,2I8,2I7,I8,I9,7I8,I7,I8)
 6804 FORMAT('NMWAVB=',I6,2I8,2I7,I8,I9,7I8,I7,I8)
 6805 FORMAT('ABSORB'/'   GAS= O3,O2  O3,NO2      O2     O2     O2'
     +        ,'     H2O',22X,'H2O     H2O     H2O     H2O     CO2'
     +        ,'     CO2    CO2  CO2,H2O,O2'
     +        /'SKDFLB (Downward Spectral Flux)'
     +        /6X,6('-'),'VIS',6('-'),2X,46('-'),'NIR',59('-'))
 6806 FORMAT('ABSORB'/'   GAS=   H2O     H2O     H2O    H2O    H2O'
     +        ,'      O2       O2      O2     CO2     CO2     CO2',18X
     +        ,'O3,NO2  O3,O2 CO2,H2O,O2'/'SKDFLB  (Downard Spectral'
     +        ,' Flux)',T110,6('-'),'VIS',5('-'))
 6807 FORMAT(' N   L=',I4,I9,I8,2I7,I8,I9,7I8,I7,I8,'       Total')
 6808 FORMAT(I2,1X,2F9.3,F8.3,2F7.3,F8.3,F9.3,7F8.3,F7.3,F8.3,F11.3)
 6809 FORMAT(/'SKUFLB (Upward Spectral Flux)'/' N   L='
     +         ,I4,I9,I8,2I7,I8,I9,7I8,I7,I8,'       Total')
 6810 FORMAT(I2,1X,2F9.3,F8.3,2F7.3,F8.3,F9.3,7F8.3,F7.3,F8.3,F11.3)
 6811 FORMAT(/'SRKALB '
     +         ,I4,I9,I8,2I7,I8,I9,7I8,I7,I8,'       Total')
 6812 FORMAT('   TOA='
     +       ,F5.4,F9.4,F8.4,2F7.4,F8.4,F9.4,7F8.4,F7.4,F8.4,F11.4)
 6813 FORMAT('ABSORB'/'  ATMO='
     +       ,F5.4,F9.4,F8.4,2F7.4,F8.4,F9.4,7F8.4,F7.4,F8.4,F11.4)
 6814 FORMAT('ABSORB'/'  SURF='
     +       ,F5.4,F9.4,F8.4,2F7.4,F8.4,F9.4,7F8.4,F7.4,F8.4,F11.4)
 6815 FORMAT('ALSURF',T133,'Sum=',F6.4
     +       /' KSLAM= ',I3,I9,I8,2I7,I8,I9,7I8,I7,I8)
 6816 FORMAT('   SRX='
     +       ,F5.4,F9.4,F8.4,2F7.4,F8.4,F9.4,7F8.4,F7.4,F8.4,F11.4)
 6817 FORMAT('   SRB='
     +       ,F5.4,F9.4,F8.4,2F7.4,F8.4,F9.4,7F8.4,F7.4,F8.4,F11.4)
 6818 FORMAT(/'At Top of Atm:  ',' COSZ  =',F6.4,14X
     +      , 2X,' SRIVIS=',F7.3,'  SROVIS=',F7.3,   '   PLAVIS=',F6.4
     +      , 2X,' SRINIR=',F7.3,'  SRONIR=',F7.3,   '   PLANIR=',F6.4)
 6819 FORMAT( 'At Bot of Atm:  ',' SRXVIS=',F6.4,1X,' SRXNIR=',F6.4
     +      , 1X,' SRDVIS=',F7.3,'  SRUVIS=',F7.3,   '   ALBVIS=',F6.4
     +      , 2X,' SRDNIR=',F7.3,'  SRUNIR=',F7.3,   '   ALBNIR=',F6.4)
 6820 FORMAT( 'In Atmosphere:  ',' (VIS=0.53*S0)',2X,'(NIR=0.47*S0)'
     +      , 1X,' SRTVIS=',F7.5,'  SRRVIS=',F7.5,   '   SRAVIS=',F6.4
     +      , 2X,' SRTNIR=',F7.5,'  SRRNIR=',F7.5,   '   SRANIR=',F6.4)
 6821 FORMAT(' ')

 6840 FORMAT(' (8B)  SPECTRAL/k-DISTRIBUTION COMPONENT BREAKDOWN'
     +      ,' FOR NET DOWNWARD SOLAR FLUX & HEATING RATE'
     +      ,T106,'SKNFLB(L,K)  SKFHRL(L,K)  SRKGAX(L,I)'/)
 6841 FORMAT('     K=',I5,2I8,2I7,I8,I9,7I8,I7,I8,'       Total')
 6842 FORMAT('  DKS0=',F6.3,2F8.3,2F7.3,F8.3,F9.3,7F8.3,F7.3,F8.3,F11.3)
 6843 FORMAT('   GAS= O3,O2  O3,NO2      O2     O2     O2'
     +        ,'     H2O',22X,'H2O     H2O     H2O     H2O     CO2',5X
     +        ,'CO2    CO2  CO3,H2O,O2'/'SKNFLB (Spectral Net Flux)')
 6844 FORMAT('   GAS=   H2O     H2O     H2O    H2O    H2O'
     +        ,'      O2       O2      O2     CO2     CO2     CO2',18X
     +        ,'O3,NO2  O3,O2 CO2,H2O,O2'
     +        /'SKDFLB (Spectral Net Flux)',T110,6('-'),'VIS',5('-'))
 6845 FORMAT(' N   L=',I4,I9,I8,2I7,I8,I9,7I8,I7,I8,'       Total')
 6846 FORMAT(I2,1X,2F9.3,F8.3,2F7.3,F8.3,F9.3,7F8.3,F7.3,F8.3,F11.3)
 6847 FORMAT(/'SKFHRL (Spectral Heating Rate)'/' N   L='
     +       ,I4,I9,I8,2I7,I8,I9,7I8,I7,I8,'       Total')
 6848 FORMAT(I2,1X,2F9.3,F8.3,2F7.3,F8.3,F9.3,7F8.3,F7.3,F8.3,F11.3)
 6849 FORMAT(/'SRKGAX (Direct Beam Spectral Absorption at Ground)'
     +       /' N   L=',I4,8I8,I9,4I8,I7,I8,'      Total')
 6850 FORMAT(I2,1X,2F9.3,F8.3,2F7.3,F8.3,F9.3,7F8.3,F7.3,F8.3,F11.3)
 6851 FORMAT( 'SRKGAD (Diffuse Spectral Absorption at Ground)')
 6852 FORMAT(I2,1X,2F9.3,F8.3,2F7.3,F8.3,F9.3,7F8.3,F7.3,F8.3,F11.3)
 6853 FORMAT( 'SRKGAD (Total Spectral Absorption at Ground)')
 6854 FORMAT(I2,1X,2F9.3,F8.3,2F7.3,F8.3,F9.3,7F8.3,F7.3,F8.3,F11.3)
 6855 FORMAT(/'Absorption at Ground by Surface-type'
     +       ,T39,'SRNFLB(1) = POCEAN * FSRNFG(1) + PEARTH * FSRNFG(2) '
     +                     ,'+  POICE * FSRNFG(3) +  PLICE * FSRNFG(4) '
     +       /T39,F7.3,'   = ',F6.4,' *',F8.3,'   + ',F6.4,' *',F8.3
     +                ,'   + ',F6.4,' *',F8.3,'   + ',F6.4,' *',F8.3)
      GO TO 9999
C-------------
  900 CONTINUE
C-------------
C
      SIGMA=5.6697D-08
      TGMEAN=POCEAN*TGO**4+PEARTH*TGE**4+PLICE*TGLI**4+POICE*TGOI**4
      TGMEAN=SQRT(TGMEAN)
      TGMEAN=SQRT(TGMEAN)
      SIGT4=SIGMA*TGMEAN**4
      ITG=TGMEAN
      WTG=TGMEAN-ITG
      ITG=ITG-ITPFT0
      DO 901 K=1,33
      BGFLUX(K)=PLANCK(ITG)-(PLANCK(ITG)-PLANCK(ITG+1))*WTG
      BGFRAC(K)=BGFLUX(K)/SIGT4
      ITG=ITG+ITNEXT
  901 CONTINUE
      DO 910 NW=1,5
      DO 903 K=1,33
      DO 902 L=1,NLP
      IF(NW.EQ.1) WFLB(L,K)=DFLB(L,K)
      IF(NW.EQ.2) WFLB(L,K)=UFLB(L,K)
      IF(NW.EQ.3) WFLB(L,K)=UFLB(L,K)-DFLB(L,K)
      IF(NW.GT.3.AND.L.GT.NL) GO TO 902
      IF(NW.EQ.4) WFLB(L,K)=WFLB(L+1,K)-WFLB(L,K)
      IF(NW.EQ.5.AND.ABS(TRFCRL(L)).LT.1.E-10) WFLB(L,K)=1.E-30
      IF(NW.EQ.5) WFLB(L,K)=WFLB(L,K)/(ABS(TRFCRL(L))+1.E-10)
  902 CONTINUE
      IF(NW.EQ.1) WFSL(K)=DFSL(K)
      IF(NW.EQ.2) WFSL(K)=UFSL(K)
      IF(NW.EQ.3) WFSL(K)=UFSL(K)-DFSL(K)
      IF(NW.EQ.4) WFSL(K)=WFSL(K)-UFLB(1,K)+DFLB(1,K)
!obso IF(NW.EQ.5.AND.ABS(TRSLCR).LT.1.E-10) WFSL(K)=1.E-30
      IF(NW.EQ.5) WFSL(K)=0.   !obsolete: =WFSL(K)/(ABS(TRSLCR)+1.E-10)
  903 CONTINUE
      DO 907 L=1,NLP
      IF(L.GT.NL.AND.NW.GT.3) GO TO 907
      ASUM1=0.
      BSUM1=0.
      CSUM1=0.
      DSUM1=0.
      ESUM1=0.
      FSUM1=0.
      SUM=0.
      DO 904 K=2,13
      ASUM1=ASUM1+  WFSL(K)
      BSUM1=BSUM1+BGFEMT(K)
      CSUM1=CSUM1+BGFLUX(K)
      DSUM1=DSUM1+BGFRAC(K)
      ESUM1=ESUM1+TRCTCA(K)
      FSUM1=FSUM1+TRGALB(K)
  904 SUM=SUM+WFLB(L,K)
      SUM1(L)=SUM
      ASUM2=0.
      BSUM2=0.
      CSUM2=0.
      DSUM2=0.
      ESUM2=0.
      FSUM2=0.
      SUM=0.
      DO 905 K=14,25
      ASUM2=ASUM2+  WFSL(K)
      BSUM2=BSUM2+BGFEMT(K)
      CSUM2=CSUM2+BGFLUX(K)
      DSUM2=DSUM2+BGFRAC(K)
      ESUM2=ESUM2+TRCTCA(K)
      FSUM2=FSUM2+TRGALB(K)
  905 SUM=SUM+WFLB(L,K)
      SUM2(L)=SUM
      ASUM3=0.
      BSUM3=0.
      CSUM3=0.
      DSUM3=0.
      ESUM3=0.
      FSUM3=0.
      SUM=0.
      DO 906 K=26,33
      ASUM3=ASUM3+  WFSL(K)
      BSUM3=BSUM3+BGFEMT(K)
      CSUM3=CSUM3+BGFLUX(K)
      DSUM3=DSUM3+BGFRAC(K)
      ESUM3=ESUM3+TRCTCA(K)
      FSUM3=FSUM3+TRGALB(K)
  906 SUM=SUM+WFLB(L,K)
      SUM3(L)=SUM
  907 CONTINUE
C
      NPAGE=1
      WRITE(KW,6901) NW,FTYPE(NW)
      WRITE(KW,6902) (K,K=1,13)
      DO 908 N=1,NLP
      L=NLP+1-N
      IF(L.GT.NL.AND.NW.GT.3) GO TO 908
      SUML=SUM1(L)+SUM2(L)+SUM3(L)+WFLB(L,1)
      WRITE(KW,6903) L,PL(L),SUML,SUM1(L),SUM2(L),SUM3(L)
     +               ,(WFLB(L,K),K=1,13)
  908 CONTINUE
      SUMA=ASUM1+ASUM2+ASUM3+  WFSL(1)
      SUMB=BSUM1+BSUM2+BSUM3+BGFEMT(1)
      SUMC=CSUM1+CSUM2+CSUM3+BGFLUX(1)
      SUMD=DSUM1+DSUM2+DSUM3+BGFRAC(1)
      SUME=ESUM1+ESUM2+ESUM3+TRCTCA(1)
      SUMF=FSUM1+FSUM2+FSUM3+TRGALB(1)
      WRITE(KW,6904) SUMA,ASUM1,ASUM2,ASUM3,(  WFSL(K),K=1,13)
      WRITE(KW,6905) SUMB,BSUM1,BSUM2,BSUM3,(BGFEMT(K),K=1,13)
      WRITE(KW,6906) SUMC,CSUM1,CSUM2,CSUM3,(BGFLUX(K),K=1,13)
      WRITE(KW,6907) SUMD,DSUM1,DSUM2,DSUM3,(BGFRAC(K),K=1,13)
      WRITE(KW,6908) SUME,ESUM1,ESUM2,ESUM3,(TRCTCA(K),K=1,13)
      WRITE(KW,6909) SUMF,FSUM1,FSUM2,FSUM3,(TRGALB(K),K=1,13)
      NPAGE=0
      WRITE(KW,6910) NPAGE
      WRITE(KW,6911) (K,K=14,33)
      DO 909 N=1,NLP
      L=NLP+1-N
      IF(L.GT.NL.AND.NW.GT.3) GO TO 909
      WRITE(KW,6912) L,(WFLB(L,K),K=14,33)
  909 CONTINUE
      WRITE(KW,6913) (  WFSL(K),K=14,33)
      WRITE(KW,6914) (BGFEMT(K),K=14,33)
      WRITE(KW,6915) (BGFLUX(K),K=14,33)
      WRITE(KW,6916) (BGFRAC(K),K=14,33)
      WRITE(KW,6917) (TRCTCA(K),K=14,33)
      WRITE(KW,6918) (TRGALB(K),K=14,33)
      LINFIL=2
      IF(NW.GT.3) LINFIL=4
      DO 911 I=1,LINFIL
      WRITE(KW,6919)
  911 CONTINUE
  910 CONTINUE
C
 6901 FORMAT(' (9.',I1,') THERMAL RADIATION: K-DISTRIBUTION'
     +       ,' BREAKDOWN FOR  ',1A8,' FLUX'/
     +       /T21,'PRINCIPAL REGION SUM',2X,'WINDOW'
     +       ,T52,'WATER VAPOR:',T76,'PRINCIPAL ABSORBER REGION'
     +       /20X,20('-'),2X,6('-'),3X,81('-'))
 6902 FORMAT(1X,'LN     PL    TOTAL    H2O   CO2    O3     K='
     +       ,I2,4X,'K=',I2,12I7)
 6903 FORMAT(1X,I2,2F8.2,3F7.2,F8.3,1X,12F7.3)
 6904 FORMAT(/' SL',  9X,4F7.2,F8.3,1X,12F7.3)
 6905 FORMAT(/' BG',  9X,4F7.2,F8.3,1X,12F7.3)
 6906 FORMAT( ' PF',  9X,4F7.2,F8.3,1X,12F7.3)
 6907 FORMAT( ' FR',  9X,4F7.2,F8.3,1X,12F7.3)
 6908 FORMAT(/' AC',  9X,4F7.2,F8.3,1X,12F7.3)
 6909 FORMAT( ' AG',  9X,4F7.2,F8.3,1X,12F7.3)
 6910 FORMAT(1I1/5X,'CARBON DIOXIDE:',T36,'PRINCIPAL ABSORBER REGION'
     +       ,T85,'OZONE:',T101,'PRINCIPAL ABSORBER REGION'
     +       /5X,76('-'),3X,48('-'))
 6911 FORMAT(1X,'LN  K=',I2,6I7,5I6,4X,'K=',I2,1I7,6I6)
 6912 FORMAT( 1X,I2,7F7.3,5F6.3,1X,2F7.3,6F6.3)
 6913 FORMAT(/' SL',7F7.3,5F6.3,1X,2F7.3,6F6.3)
 6914 FORMAT(/' BG',7F7.3,5F6.3,1X,2F7.3,6F6.3)
 6915 FORMAT( ' PF',7F7.3,5F6.3,1X,2F7.3,6F6.3)
 6916 FORMAT( ' FR',7F7.3,5F6.3,1X,2F7.3,6F6.3)
 6917 FORMAT(/' AC',7F7.3,5F6.3,1X,2F7.3,6F6.3)
 6918 FORMAT( ' AG',7F7.3,5F6.3,1X,2F7.3,6F6.3)
 6919 FORMAT(' ')
      RETURN
C-------------
 1000 CONTINUE
C-------------
C
 9999 CONTINUE
      RETURN
      END SUBROUTINE WRITER

      SUBROUTINE WRITET(KWRU,INDEX,JYRREF,JYRNOW,JMONTH,KLIMIT)
C
      !INCLUDE 'BR00B.COM'   ! use RE001
C
C     ------------------------------------------------------------------
C     WRITET  GHG, Solar UV, Ozone, Aerosol Trend Diagnostic Information
C
C         INDEX
C           1     GHG DT0 Trends / FULGAS Ratios for CO2,NO2,CH4,F11,F12
C           2     GHG DF  Change / Ann Increase Rate CO2,NO2,CH4,F11,F12
C           3     Lean Solar Constant, UV Spectral Variation Time Trends
C           4     Ozone Zonal-mean (Latitude and Vertical) Distributions
C           5     Ozone Surface-150mb, 150mb-TOA, Column Longitude Distr
C                 A  O3 (Wang-Jacobs) Relative Longitudinal Distribution
C                 B  O3 (London-NCAR) Relative Longitudinal Distribution
C                 C  O3 (W-J, London) Relative Longitudinal Distribution
C           6     Tropospheric Climatology Aerosol Latitude/Height Distr
C                 A  Zonal-mean Extinction Optical Depth
C                 B  Zonal-mean Single Scattering Albedo
C                 C  Zonal-mean Asymmetry Parameter
C           7     Tropospheric Desert Dust Aerosol Latitude/Height Distr
C                 A  Zonal-mean Extinction Optical Depth
C                 B  Zonal-mean Single Scattering Albedo
C                 C  Zonal-mean Asymmetry Parameter
C           8     Stratospheric (Volcanic) Aerosol Latitude/Height Distr
C                 A  Zonal-mean Extinction Optical Depth
C                 B  Zonal-mean Single Scattering Albedo
C                 C  Zonal-mean Asymmetry Parameter
C           9     Total Column Atmospheric Aerosol Latitude/Height Distr
C                 A  Zonal-mean Extinction Optical Depth
C                 B  Zonal-mean Single Scattering Albedo
C                 C  Zonal-mean Asymmetry Parameter
C        NOTE:
C                 Time Trend (year) Specification is by JYRREF to JYRNOW
C                 Time Specification (O3,Aerosol) is by JYRREF to JMONTH
C                      (If JMONTH = 0, JDAY is used)
C
C                 INDEX < 10 is selective, INDEX > 10 is digit inclusive
C                 KLIMIT = 0 full output,  KLIMIT > 0 abbreviated output
C                 KWRU directs the output to selected (KWRU) file number
C     ------------------------------------------------------------------
C
      CHARACTER*32 CHAER(4)

      DIMENSION WREF(7),WDAT(7),WPPM(7),DWM2(7),XDT0(5),XRAT(5)
      DIMENSION QX(49,LX),QS(49,LX),QG(49,LX),QP(49,LX),O3(49,LX)
      DIMENSION QXCOL(49),QSCOL(49),QGCOL(49),QPCOL(49),O3COL(49)
      DIMENSION SFL0(5),SFLX(5),DFLX(5),RFLX(5),O3L(46,72),LO3(36)
C
      DATA NSW1/24/, NSW2/32/, NSW3/40/, NSW4/48/
C
      DATA CHAER/'Tropospheric Climatology Aerosol',
     +           'Tropospheric Desert Dust Aerosol',
     +           'Stratospheric (Volcanic) Aerosol',
     +           'Total Column Atmospheric Aerosol'/
      SAVE  NSW1,NSW2,NSW3,NSW4,CHAER
C
      KW=KWRU
      INDJ=MOD(INDEX,10)
      IF(INDJ.LT.1) INDJ=10
      INDI=1
      IF(INDEX.EQ.0)  INDJ=1
      IF(INDEX.LT.11) INDI=INDJ
      DO 9999 INDX=INDI,INDJ
C
      GO TO (100,100,300,400,500,600,600,600,600,1000),INDX
C
C-------------
  100 CONTINUE
C-------------
C
      KINDEX=INDX
      DO 110 I=1,5
      WREF(I)=XREF(I)
      WPPM(I)=PPMV80(I+4)
  110 CONTINUE
      WREF(6)=PPMV80(11)*1000.D0
      WREF(7)=PPMV80(12)*1000.D0
      WPPM(1)=PPMV80(2)
      WPPM(6)=PPMV80(11)
      WPPM(7)=PPMV80(12)
C
      YREF11=PPMV80(11)
      ZREF12=PPMV80(12)
C
      JJDAYG=184
C
      IF(KINDEX.EQ.1) THEN
      WRITE(KW,6101) JJDAYG
 6101 FORMAT(/1X,'(1)=INDEX'
     +      ,T12,'JDAY=',I3,'   RCM RAD EQUIL NO-FEEDBACK DT0'
     +      ,T55,'PRESENT TREND UPDGHG INPUT DATA TO GCM'
     +      ,T96,'FULGAS FACTOR RELATIVE TO 1980 AMOUNTS')
       WRITE(KW,6102) KTREND
 6102 FORMAT(1X,'KTREND=',I2,1X,40('-'),3X,38('-'),3X,38('-')
     +      /1X,'YEAR DTSUM  *DTCO2   DTN2O   DTCH4   DTF11   DTF12'
     +         ,         '   PPMCO2  PPMN20  PPMCH4  PPBF11  PPBF12'
     +         ,         '   FULCO2  FULN2O  FULCH4  FULF11  FULF12')
      ENDIF
C
      IF(KINDEX.EQ.2) THEN
      WRITE(KW,6201) JJDAYG
 6201 FORMAT(/1X,'(2)=INDEX'
     +      ,T12,'JDAY=',I3,'   RCM EQ NO-FEEDBACK DFLUX W/M2'
     +      ,T55,'PRESENT TREND UPDGHG INPUT DATA TO GCM'
     +      ,T96,'ANNUAL CHANGE RATE OF TRACE GAS AMOUNT')
      WRITE(KW,6202) KTREND
 6202 FORMAT(1X,'KTREND=',I2,1X,40('-'),3X,38('-'),3X,38('-')
     +      /1X,'YEAR DTSUM  *DTCO2   DTN2O   DTCH4   DTF11   DTF12'
     +         ,         '   PPMCO2  PPMN20  PPMCH4  PPBF11  PPBF12'
     +         ,         '   RATCO2  RATN2O  RATCH4  RATF11  RATF12')
       ENDIF
C
      JYEARG=JYRREF-1
      CALL UPDGHG(JYEARG,JJDAYG)
      YNOW11=FULGAS(11)*YREF11
      ZNOW12=FULGAS(12)*ZREF12
      XDT0(4)=XDT0(4)+0.066*(YNOW11-YREF11)
      XDT0(5)=XDT0(5)+0.084*(ZNOW12-ZREF12)
C
      DO 120 I=1,5
      WDAT(I)=XNOW(I)
  120 CONTINUE
C
      DO 230 J=JYRREF,JYRNOW
      KWSKIP=0
      IF(J.GT.JYRREF) KWSKIP=KLIMIT
      IF(J.EQ.1980)   KWSKIP=0
      IF(J.EQ.JYRNOW) KWSKIP=0
      JYEARG=J
      CALL UPDGHG(JYEARG,JJDAYG)
      YNOW11=FULGAS(11)*YREF11
      ZNOW12=FULGAS(12)*ZREF12
      XDT0(4)=XDT0(4)+0.066*(YNOW11-YREF11)
      XDT0(5)=XDT0(5)+0.084*(ZNOW12-ZREF12)
      SDT0=0.D0
      DFL0=0.D0
      DO 220 I=1,5
      SDT0=SDT0+XDT0(I)
      XRAT(I)=(XNOW(I)-WDAT(I))/(1.D-10+WDAT(I))
      IF(XRAT(I).GT.9.9999) XRAT(I)=9.9999
      WDAT(I)=XNOW(I)
      DWM2(I)=XDT0(I)*3.300D0
      DFL0=DFL0+DWM2(I)
  220 CONTINUE
      IYEAR=JYEARG
      IF(KINDEX.EQ.1) THEN
      IF(KWSKIP.EQ.0)
     +WRITE(KW,6103) IYEAR,SDT0,(XDT0(I),I=1,5),(XNOW(I),I=1,5)
     +              ,FULGAS(2),(FULGAS(I),I=6,9)
 6103 FORMAT(1X,I4,F6.3,5F8.4,1X,F8.2,4F8.4,1X,5F8.4)
      ENDIF
      IF(KINDEX.EQ.2) THEN
      IF(KWSKIP.EQ.0)
     +WRITE(KW,6203) IYEAR,DFL0,(DWM2(I),I=1,5),(XNOW(I),I=1,5)
     +              ,(XRAT(I),I=1,5)
 6203 FORMAT(1X,I4,F6.3,5F8.4,1X,F8.2,4F8.4,1X,5F8.4)
      ENDIF
      NSPACE=IYEAR-(IYEAR/10)*10
      IF(KLIMIT.GT.0) GO TO 230
      IF(NSPACE.EQ.0) WRITE(KW,6104)
 6104 FORMAT(' ')
  230 CONTINUE
      GO TO 9999
C
C-------------
  300 CONTINUE
C-------------
C
      if(ksolar.lt.0) go to 9999
      LMO=(1950-iy1S0)*12+1
      if(ksolar.gt.1) LMO=nint(1950 - yr1s0 + 1.5)
      DO 310 I=1,5
      SFL0(I)=0.D0
  310 CONTINUE
      DO 320 K=1,190
      IF(K.LE.NSW1)              SFL0(1)=SFL0(1)+UVLEAN(LMO,K)*DSLEAN(K)
      IF(K.GT.NSW1.AND.K.LE.NSW2)SFL0(2)=SFL0(2)+UVLEAN(LMO,K)*DSLEAN(K)
      IF(K.GT.NSW2.AND.K.LE.NSW3)SFL0(3)=SFL0(3)+UVLEAN(LMO,K)*DSLEAN(K)
      IF(K.GT.NSW3.AND.K.LE.NSW4)SFL0(4)=SFL0(4)+UVLEAN(LMO,K)*DSLEAN(K)
                                 SFL0(5)=SFL0(5)+UVLEAN(LMO,K)*DSLEAN(K)
  320 CONTINUE
C
      if(ksolar.eq.2)
     *   WRITE(KW,6299) int(yr1s0),int(yr2s0),JYRREF,JYRNOW,SFL0(5)
      if(ksolar.lt.2) WRITE(KW,6300) JYRREF,JYRNOW,SFL0(5)
 6299 FORMAT(/' (3)=INDEX  Annual-mean Solar flux (from J.Lean annual'
     +      ,I6,'-',I4,' data) for JYRREF=',I4,' to JYRNOW=',I4,'  mid'
     +      ,' 1950 Ref S00WM2=',F9.4/12X,'Solar UV Spectral Flux W/m2'
     +      ,T57,'Delta Solar UV Spectral Flux W/m2'
     +      ,T97,'Solar UV Spectral Flux Ratios'
     +      /'  YEAR    0-280 280-320 320-360 360-400   Total '
     +      ,6X,'0-280 280-320 320-360 360-400   Total '
     +      ,4X,'0-280 280-320 320-360 360-400   Total ')
 6300 FORMAT(/' (3)=INDEX  Annual-mean Solar flux (from J.Lean monthly'
     +      ,' 1882-1998 data) for JYRREF=',I4,' to JYRNOW=',I4,'  Jan'
     +      ,' 1950 Ref S00WM2=',F9.4/12X,'Solar UV Spectral Flux W/m2'
     +      ,T57,'Delta Solar UV Spectral Flux W/m2'
     +      ,T97,'Solar UV Spectral Flux Ratios'
     +      /'  YEAR    0-280 280-320 320-360 360-400   Total '
     +      ,6X,'0-280 280-320 320-360 360-400   Total '
     +      ,4X,'0-280 280-320 320-360 360-400   Total ')
C
      if(ksolar.lt.2) then
        mavg = 12
        iyr1 = iy1s0
        lmax = ms0x
        icyc = mcycs0
      else
        mavg = 1
        iyr1 = yr1s0
        lmax = nint( yr2s0-yr1s0+1 )
        icyc = icycs0
      end if
      DO 370 J=JYRREF,JYRNOW
      KWSKIP=0
      IF(J.GT.JYRREF) KWSKIP=KLIMIT
      IF(J.EQ.JYRNOW) KWSKIP=0
      JYEARS=J
      DO 330 I=1,5
      SFLX(I)=0.D0
  330 CONTINUE
      LMO=(JYEARS-iyr1)*mavg
      DO 350 M=1,mavg
      LMO=LMO+1
      IF(LMO.GT.lmax) LMO=LMO-icyc*((LMO-lmax+icyc-1)/icyc)
      IF(LMO.LT.1) LMO=LMO+icyc*((icyc-LMO)/icyc)
      DO 340 K=1,190
      IF(K.LE.NSW1)              SFLX(1)=SFLX(1)+UVLEAN(LMO,K)*DSLEAN(K)
      IF(K.GT.NSW1.AND.K.LE.NSW2)SFLX(2)=SFLX(2)+UVLEAN(LMO,K)*DSLEAN(K)
      IF(K.GT.NSW2.AND.K.LE.NSW3)SFLX(3)=SFLX(3)+UVLEAN(LMO,K)*DSLEAN(K)
      IF(K.GT.NSW3.AND.K.LE.NSW4)SFLX(4)=SFLX(4)+UVLEAN(LMO,K)*DSLEAN(K)
                                 SFLX(5)=SFLX(5)+UVLEAN(LMO,K)*DSLEAN(K)
  340 CONTINUE
  350 CONTINUE
      DO 360 I=1,5
      SFLX(I)=SFLX(I)/mavg
      DFLX(I)=SFLX(I)-SFL0(I)
      RFLX(I)=SFLX(I)/SFL0(I)
  360 CONTINUE
      IF(KWSKIP.EQ.0)
     +WRITE(KW,6301) JYEARS,(SFLX(I),I=1,5),(DFLX(I),I=1,5)
     +                     ,(RFLX(I),I=1,5)
 6301 FORMAT(2X,I4,1X,4F8.4,F10.4,2X,5F8.4,2X,5F8.5)
      NSPACE=JYEARS-(JYEARS/10)*10
      IF(KLIMIT.GT.0) GO TO 370
      IF(NSPACE.EQ.0) WRITE(KW,6302)
 6302 FORMAT(' ')
  370 CONTINUE
      GO TO 9999
C
C-------------
  400 CONTINUE
C-------------
C
      JJDAYO=JMONTH*30-15
      IF(JMONTH.LT.1) JJDAYO=JDAY
      CALL UPDO3D(JYRREF,JJDAYO)
      DO 450 J=1,46
      DO 410 L=1,NL
      O3(J,L)=0.D0
  410 CONTINUE
      JLAT=J
      DO 430 I=1,72
      ILON=I
      CALL GETO3D
      DO 420 L=1,NL
      O3(J,L)=O3(J,L)+U0GAS(L,3)/72.D0
  420 CONTINUE
  430 CONTINUE
      SUMO3=0.D0
      DO 440 L=1,NL
      SUMO3=SUMO3+O3(J,L)
  440 CONTINUE
      O3COL(J)=SUMO3
  450 CONTINUE
      CALL BOXAV1(DLAT46,O3COL,46, 1,23,QOSH)
      CALL BOXAV1(DLAT46,O3COL,46,24,46,QONH)
      CALL BOXAV1(DLAT46,O3COL,46, 1,46,QOGL)
      O3COL(47)=QOSH
      O3COL(48)=QONH
      O3COL(49)=QOGL
      DO 460 L=1,NL
      CALL BOXAV1(DLAT46,O3(1,L),46, 1,23,QOSH)
      CALL BOXAV1(DLAT46,O3(1,L),46,24,46,QONH)
      CALL BOXAV1(DLAT46,O3(1,L),46, 1,46,QOGL)
      O3(47,L)=QOSH
      O3(48,L)=QONH
      O3(49,L)=QOGL
  460 CONTINUE
C
      IF(KLIMIT.GT.0)
     +WRITE(KW,6400) JYRREF,JJDAYO,JMONTH,MOZONE,(L,L=2,15)
 6400 FORMAT(/' (4)=INDEX  JYRREF=',I5,'  JDAY=',I3,'   JMONTH=',I2
     +      ,T50,' Ozone: Zonal-mean Vertical Distribution (cmSTP)'
     +      ,T126,'MOZONE=',I1/'  JLAT DLAT46   COLUMN  L =   1',14I7)
      IF(KLIMIT.LT.1)
     +WRITE(KW,7400) JYRREF,JJDAYO,JMONTH,MOZONE
     +              ,(PLB0(I),I=1,15),(L,L=2,15)
 7400 FORMAT(/' (4)=INDEX  JYRREF=',I5,'  JDAY=',I3,'   JMONTH=',I2
     +      ,T50,' Ozone: Zonal-mean Vertical Distribution (cmSTP)'
     +      ,T126,'MOZONE=',I1//21X,'PLB0 =',F6.1,9F7.1,5F7.2
     +      /'  JLAT DLAT46   COLUMN  L =   1',14I7)
C
      DO 470 JJ=1,46
      J=47-JJ
      IF(KLIMIT.GT.0) GO TO 470
      WRITE(KW,6401) J,DLAT46(J),O3COL(J),(O3(J,L),L=1,15)
 6401 FORMAT(I5,F8.2,F9.5,4X,15F7.5)
  470 CONTINUE
      IF(KLIMIT.LT.1) WRITE(KW,6402)
 6402 FORMAT(' ')
      WRITE(KW,6403) O3COL(48),(O3(48,L),L=1,15)
 6403 FORMAT(11X,'NH',F9.5,4X,15F7.5)
      IF(KLIMIT.LT.1) WRITE(KW,6402)
      WRITE(KW,6404) O3COL(47),(O3(47,L),L=1,15)
 6404 FORMAT(11X,'SH',F9.5,4X,15F7.5)
      IF(KLIMIT.LT.1) WRITE(KW,6402)
      WRITE(KW,6405) O3COL(49),(O3(49,L),L=1,15)
 6405 FORMAT( 7X,'GLOBAL',F9.5,4X,15F7.5)
      GO TO 9999
C
C
C-------------
  500 CONTINUE
C-------------
C
      JJDAYO=JMONTH*30-15
      IF(JMONTH.LT.1) JJDAYO=JDAY
      CALL UPDO3D(JYRREF,JJDAYO)
      DO 590 N=1,3
      N1=1
      N2=8
      IF(N.EQ.2) N1=9
      IF(N.GT.1) N2=NL
      DO 530 J=1,46
      JLAT=J
      DO 520 I=1,72
      ILON=I
      CALL GETO3D
      SUMO3=0.D0
      DO 510 L=N1,N2
      SUMO3=SUMO3+U0GAS(L,3)
  510 CONTINUE
      O3L(J,I)=SUMO3
  520 CONTINUE
  530 CONTINUE
      DO 560 J=1,46
      SUMO3=0.D0
      DO 540 I=1,72
      SUMO3=SUMO3+O3L(J,I)/72.D0
  540 CONTINUE
      DO 550 I=1,72
      O3L(J,I)=O3L(J,I)/SUMO3
  550 CONTINUE
  560 CONTINUE
C
      IF(N.EQ.1)
     +WRITE(KW,6510) JYRREF,JJDAY,JMONTH,MOZONE,(I,I=10,310,10)
 6510 FORMAT(/'  5A=INDEX  JYEAR=',I5,'  JDAY=',I3,'   JMONTH=',I2
     +      ,T50,' Ozone Longitudinal Variation:  Troposphere'
     +      ,' (Wang-Jacobs) Surf to 150 mb',T126,'MOZONE=',I1
     +      /'  J LON=0',31I4)
      IF(N.EQ.2)
     +WRITE(KW,6520) JYRREF,JJDAY,JMONTH,MOZONE,(I,I=10,310,10)
 6520 FORMAT(/'  5B=INDEX  JYEAR=',I5,'  JDAY=',I3,'   JMONTH=',I2
     +      ,T50,' Ozone Longitudinal Variation:  Stratosphere'
     +      ,' (London-NCAR) 150 mb to TOA',T126,'MOZONE=',I1
     +      /'  J LON=0',31I4)
      IF(N.EQ.3.AND.KLIMIT.LT.1)
     +WRITE(KW,6530) JYRREF,JJDAY,JMONTH,MOZONE,(I,I=10,310,10)
 6530 FORMAT(/'  5C=INDEX  JYEAR=',I5,'  JDAY=',I3,'   JMONTH=',I2
     +      ,T50,' Ozone Longitudinal Variation:  Total Column'
     +      ,' (W-J/London) Surface to TOA',T126,'MOZONE=',I1
     +      /'  J LON=0',31I4)
      IF(KLIMIT.LT.1) WRITE(KW,6540)
 6540 FORMAT(' ')
C
      IF(N.EQ.3.AND.KLIMIT.GT.0) GO TO 590
      DO 580 JJ=1,46
      J=47-JJ
      KWSKIP=KLIMIT
      IF(J.EQ.36) KWSKIP=0
      IF(J.EQ.24) KWSKIP=0
      IF(J.EQ.12) KWSKIP=0
      DO 570 I=1,36
      II=I*2-1
      LO3(I)=O3L(J,II)*100.D0+0.5D0
  570 CONTINUE
      IF(KWSKIP.EQ.0)
     +WRITE(KW,6501) J,(LO3(I),I=1,32)
 6501 FORMAT(I4,1X,36I4)
  580 CONTINUE
  590 CONTINUE
C
      GO TO 9999
C
C-------------
  600 CONTINUE
C-------------
C
      KAEROS=4
      IF(INDX.EQ.6) KAEROS=1
      IF(INDX.EQ.7) KAEROS=2
      IF(INDX.EQ.8) KAEROS=3
      JJDAY=JMONTH*30-15
      IF(JMONTH.LT.1) JJDAY=JDAY
      K=6
      IF(KAEROS.EQ.1.OR.KAEROS.GT.3) CALL UPDAER(JYRREF,JJDAY)
      IF(KAEROS.EQ.2.OR.KAEROS.GT.3) CALL UPDDST(JYRREF,JJDAY)
      IF(KAEROS.EQ.3.OR.KAEROS.GT.3) CALL UPDVOL(JYRREF,JJDAY)
C
      DO 650 J=1,46
      DO 610 L=1,NL
      QX(J,L)=0.D0
      QS(J,L)=0.D0
      QG(J,L)=0.D0
  610 CONTINUE
      JLAT=J
      DO 630 I=1,72
      ILON=I
      IF(KAEROS.EQ.1.OR.KAEROS.GT.3) CALL GETAER
      IF(KAEROS.EQ.2.OR.KAEROS.GT.3) CALL GETDST
      IF(KAEROS.EQ.3.OR.KAEROS.GT.3) CALL GETVOL
      DO 620 L=1,NL
      IF(KAEROS.EQ.1.OR.KAEROS.GT.3) QX(J,L)=QX(J,L)+SRAEXT(L,K)/72.D0
      IF(KAEROS.EQ.2.OR.KAEROS.GT.3) QX(J,L)=QX(J,L)+SRDEXT(L,K)/72.D0
      IF(KAEROS.EQ.3.OR.KAEROS.GT.3) QX(J,L)=QX(J,L)+SRVEXT(L,K)/72.D0
      IF(KAEROS.EQ.1.OR.KAEROS.GT.3) QS(J,L)=QS(J,L)+SRASCT(L,K)/72.D0
      IF(KAEROS.EQ.2.OR.KAEROS.GT.3) QS(J,L)=QS(J,L)+SRDSCT(L,K)/72.D0
      IF(KAEROS.EQ.3.OR.KAEROS.GT.3) QS(J,L)=QS(J,L)+SRVSCT(L,K)/72.D0
      IF(KAEROS.EQ.1.OR.KAEROS.GT.3) QG(J,L)=QG(J,L)+SRAGCB(L,K)
     +                                              *SRASCT(L,K)/72.D0
      IF(KAEROS.EQ.2.OR.KAEROS.GT.3) QG(J,L)=QG(J,L)+SRDGCB(L,K)
     +                                              *SRDSCT(L,K)/72.D0
      IF(KAEROS.EQ.3.OR.KAEROS.GT.3) QG(J,L)=QG(J,L)+SRVGCB(L,K)
     +                                              *SRVSCT(L,K)/72.D0
  620 CONTINUE
  630 CONTINUE
      SUMXL=1.D-10
      SUMSL=1.D-20
      SUMGL=1.D-20
      DO 640 L=1,NL
      SUMXL=SUMXL+QX(J,L)
      SUMSL=SUMSL+QS(J,L)
      SUMGL=SUMGL+QG(J,L)
      QG(J,L)=(1.D-20+QG(J,L))/(1.D-10+QS(J,L))
      QP(J,L)=(1.D-20+QS(J,L))/(1.D-10+QX(J,L))
      IF(QP(J,L).GT.0.99999D0) QP(J,L)=0.99999D0
  640 CONTINUE
      QXCOL(J)=SUMXL
      QSCOL(J)=SUMSL
      QGCOL(J)=(1.D-15+SUMGL)/(1.D-05+SUMSL)
      QPCOL(J)=(1.D-20+SUMSL)/(1.D-10+SUMXL)
  650 CONTINUE
      CALL BOXAV1(DLAT46,QXCOL,46, 1,23,QXSH)
      CALL BOXAV1(DLAT46,QXCOL,46,24,46,QXNH)
      CALL BOXAV1(DLAT46,QXCOL,46, 1,46,QXGL)
      QXCOL(47)=QXSH
      QXCOL(48)=QXNH
      QXCOL(49)=QXGL
      CALL BOXAV1(DLAT46,QSCOL,46, 1,23,QSSH)
      CALL BOXAV1(DLAT46,QSCOL,46,24,46,QSNH)
      CALL BOXAV1(DLAT46,QSCOL,46, 1,46,QSGL)
      QSCOL(47)=QSSH
      QSCOL(48)=QSNH
      QSCOL(49)=QSGL
      CALL BOXAV2(DLAT46,QXCOL,QPCOL,46, 1,23,QPSH)
      CALL BOXAV2(DLAT46,QXCOL,QPCOL,46,24,46,QPNH)
      CALL BOXAV2(DLAT46,QXCOL,QPCOL,46, 1,46,QPGL)
      QPCOL(47)=QPSH
      QPCOL(48)=QPNH
      QPCOL(49)=QPGL
      CALL BOXAV2(DLAT46,QSCOL,QGCOL,46, 1,23,QGSH)
      CALL BOXAV2(DLAT46,QSCOL,QGCOL,46,24,46,QGNH)
      CALL BOXAV2(DLAT46,QSCOL,QGCOL,46, 1,46,QGGL)
      QGCOL(47)=QGSH
      QGCOL(48)=QGNH
      QGCOL(49)=QGGL
      DO 660 L=1,NL
      CALL BOXAV1(DLAT46,QX(1,L),46, 1,23,QXSH)
      CALL BOXAV1(DLAT46,QX(1,L),46,24,46,QXNH)
      CALL BOXAV1(DLAT46,QX(1,L),46, 1,46,QXGL)
      QX(47,L)=QXSH
      QX(48,L)=QXNH
      QX(49,L)=QXGL
      CALL BOXAV1(DLAT46,QS(1,L),46, 1,23,QSSH)
      CALL BOXAV1(DLAT46,QS(1,L),46,24,46,QSNH)
      CALL BOXAV1(DLAT46,QS(1,L),46, 1,46,QSGL)
      QS(47,L)=QSSH
      QS(48,L)=QSNH
      QS(49,L)=QSGL
      CALL BOXAV2(DLAT46,QX(1,L),QP(1,L),46, 1,23,QPSH)
      CALL BOXAV2(DLAT46,QX(1,L),QP(1,L),46,24,46,QPNH)
      CALL BOXAV2(DLAT46,QX(1,L),QP(1,L),46, 1,46,QPGL)
      QP(47,L)=QPSH
      QP(48,L)=QPNH
      QP(49,L)=QPGL
      CALL BOXAV2(DLAT46,QS(1,L),QG(1,L),46, 1,23,QGSH)
      CALL BOXAV2(DLAT46,QS(1,L),QG(1,L),46,24,46,QGNH)
      CALL BOXAV2(DLAT46,QS(1,L),QG(1,L),46, 1,46,QGGL)
      QG(47,L)=QGSH
      QG(48,L)=QGNH
      QG(49,L)=QGGL
  660 CONTINUE
C
      KA=KAEROS
      IF(KLIMIT.GT.0)
     +WRITE(KW,6600) INDX,JYRREF,JJDAY,JMONTH,CHAER(KA),(L,L=2,15)
 6600 FORMAT(/I3,'A=INDEX JYEAR=',I5,'  JDAY=',I3,'   JMONTH=',I2
     +      ,T50,' ZONAL MEAN AEROSOL OPTICAL DEPTH',T100,A32
     +      /'  JLAT DLAT46   COLUMN  L =   1',14I7)
      IF(KLIMIT.LT.1)
     +WRITE(KW,7600) INDX,JYRREF,JJDAY,JMONTH,CHAER(KA)
     +              ,(PLB0(I),I=1,15),(L,L=2,15)
 7600 FORMAT(/I3,'A=INDEX JYEAR=',I5,'  JDAY=',I3,'   JMONTH=',I2
     +      ,T50,' ZONAL MEAN AEROSOL OPTICAL DEPTH',T100,A32
     +      //21X,'PLB0 =',F6.1,9F7.1,5F7.2
     +      /'  JLAT DLAT46   COLUMN  L =   1',14I7)
C
      IF(KLIMIT.LT.1) THEN
      DO 670 JJ=1,46
      J=47-JJ
      WRITE(KW,6601) J,DLAT46(J),QXCOL(J),(QX(J,L),L=1,15)
 6601 FORMAT(I5,F8.2,F9.5,4X,15F7.5)
  670 CONTINUE
      WRITE(KW,6602) QXCOL(48),(QX(48,L),L=1,15)
 6602 FORMAT(/11X,'NH',F9.5,4X,15F7.5)
      WRITE(KW,6603) QXCOL(47),(QX(47,L),L=1,15)
 6603 FORMAT(/11X,'SH',F9.5,4X,15F7.5)
      WRITE(KW,6604) QXCOL(49),(QX(49,L),L=1,15)
 6604 FORMAT(/7X,'GLOBAL',F9.5,4X,15F7.5)
      ENDIF
      IF(KLIMIT.GT.0) THEN
      WRITE(KW,6605) QXCOL(48),(QX(48,L),L=1,15)
 6605 FORMAT( 11X,'NH',F9.5,4X,15F7.5)
      WRITE(KW,6606) QXCOL(47),(QX(47,L),L=1,15)
 6606 FORMAT( 11X,'SH',F9.5,4X,15F7.5)
      WRITE(KW,6607) QXCOL(49),(QX(49,L),L=1,15)
 6607 FORMAT( 7X,'GLOBAL',F9.5,4X,15F7.5)
      ENDIF
C
      IF(KLIMIT.GT.0) GO TO 699
      WRITE(KW,6610) INDX,JYRREF,JJDAY,JMONTH,CHAER(KA)
     +              ,(PLB0(I),I=1,15),(L,L=2,15)
 6610 FORMAT(/I3,'B=INDEX JYEAR=',I5,'  JDAY=',I3,'   JMONTH=',I2
     +      ,T50,' ZONAL MEAN AEROSOL SINGLE SCATTERING ALBEDO'
     +      ,T100,A32//21X,'PLB0 =',F6.1,9F7.1,5F7.2
     +      /'  JLAT DLAT46   COLUMN  L =   1',14I7)
C
      DO 680 JJ=1,46
      J=47-JJ
      WRITE(KW,6611) J,DLAT46(J),QPCOL(J),(QP(J,L),L=1,15)
 6611 FORMAT(I5,F8.2,F9.5,4X,15F7.5)
  680 CONTINUE
      WRITE(KW,6612) QPCOL(48),(QP(48,L),L=1,15)
 6612 FORMAT(/11X,'NH',F9.5,4X,15F7.5)
      WRITE(KW,6613) QPCOL(47),(QP(47,L),L=1,15)
 6613 FORMAT(/11X,'SH',F9.5,4X,15F7.5)
      WRITE(KW,6614) QPCOL(49),(QP(49,L),L=1,15)
 6614 FORMAT(/7X,'GLOBAL',F9.5,4X,15F7.5)
C
      WRITE(KW,6620) INDX,JYRREF,JJDAY,JMONTH,CHAER(KA)
     +              ,(PLB0(I),I=1,15),(L,L=2,15)
 6620 FORMAT(/I3,'C=INDEX JYEAR=',I5,'  JDAY=',I3,'   JMONTH=',I2
     +      ,T50,' ZONAL MEAN AEROSOL ASYMMETRY PARAMETER'
     +      ,T100,A32//21X,'PLB0 =',F6.1,9F7.1,5F7.2
     +      /'  JLAT DLAT46   COLUMN  L =   1',14I7)
C
      DO 690 JJ=1,46
      J=47-JJ
      WRITE(KW,6621) J,DLAT46(J),QGCOL(J),(QG(J,L),L=1,15)
 6621 FORMAT(I5,F8.2,F9.5,4X,15F7.5)
  690 CONTINUE
      WRITE(KW,6622) QGCOL(48),(QG(48,L),L=1,15)
 6622 FORMAT(/11X,'NH',F9.5,4X,15F7.5)
      WRITE(KW,6623) QGCOL(47),(QG(47,L),L=1,15)
 6623 FORMAT(/11X,'SH',F9.5,4X,15F7.5)
      WRITE(KW,6624) QGCOL(49),(QG(49,L),L=1,15)
 6624 FORMAT(/7X,'GLOBAL',F9.5,4X,15F7.5)
  699 CONTINUE
      GO TO 9999
C
 1000 CONTINUE
C
 9999 CONTINUE
C
      RETURN
      END SUBROUTINE WRITET

      SUBROUTINE BOXAV1(DEGLAT,TAULAT,NLAT,JALIM,JBLIM,TAU)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION DEGLAT(NLAT),TAULAT(NLAT)
C
C--------------------------------------------------------------------
C     BOXAV1  Performs:
C                       Latitudinal average (area-weighted) of TAULAT
C
C              DEGLAT   Center latitude of grid-box variable (TAULAT)
C                       of the form:  DEGLAT = -90+(J-1)*180/(NLAT-1)
C
C              TAULAT   Zonal average value is constant over grid-bos
C
C        JALIM, JBLIM   Latitude boxes for which (TAULAT) is averaged
C
C                 TAU   Area-weighted (TAULAT) latitude average value
C--------------------------------------------------------------------
C
      ASUM=0.D0
      TSUM=0.D0
      PI=ACOS(-1.D0)
      RADIAN=180.D0/PI
      J1=JALIM-1
      IF(J1.LT.1) J1=1
      RLAT1=(0.5D0*(DEGLAT(J1)+DEGLAT(JALIM))+90.D0)/RADIAN
      ALAT1=SIN(RLAT1)
      DO 110 J=JALIM,JBLIM
      J2=J+1
      IF(J2.GT.NLAT) J2=NLAT
      RLAT2=(0.5D0*(DEGLAT(J)+DEGLAT(J2))+90.D0)/RADIAN
      ALAT2=SIN(RLAT2)
      ALATJ=0.5D0*(ALAT1+ALAT2)/(RLAT2-RLAT1)
      ASUM=ASUM+ALATJ
      TSUM=TSUM+ALATJ*TAULAT(J)
      RLAT1=RLAT2
      ALAT1=ALAT2
  110 CONTINUE
      TAU=TSUM/ASUM
      RETURN
      END SUBROUTINE BOXAV1

      SUBROUTINE BOXAV2(DEGLAT,TAULAT,SIZLAT,NLAT,JALIM,JBLIM,SIZ)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION DEGLAT(NLAT),TAULAT(NLAT),SIZLAT(NLAT)
C
C--------------------------------------------------------------------
C     BOXAV2  Performs:
C                       TAULAT-weighted latitudinal average of SIZLAT
C
C              DEGLAT   Center latitude of grid-box variable (TAULAT)
C                       of the form:  DEGLAT = -90+(J-1)*180/(NLAT-1)
C
C              TAULAT   Zonal average value is constant over grid-box
C              SIZLAT   Zonal average value is constant over grid-box
C
C        JALIM, JBLIM   Latitude boxes for which variable is averaged
C
C                 SIZ   TAULAT-weighted latitudinal average of SIZLAT
C--------------------------------------------------------------------
C
      ASUM=0.D0
      TSUM=0.D0
      PI=ACOS(-1.D0)
      RADIAN=180.D0/PI
      J1=JALIM-1
      IF(J1.LT.1) J1=1
      RLAT1=(0.5D0*(DEGLAT(J1)+DEGLAT(JALIM))+90.D0)/RADIAN
      ALAT1=SIN(RLAT1)
      DO 110 J=JALIM,JBLIM
      J2=J+1
      IF(J2.GT.NLAT) J2=NLAT
      RLAT2=(0.5D0*(DEGLAT(J)+DEGLAT(J2))+90.D0)/RADIAN
      ALAT2=SIN(RLAT2)
      ALATJ=0.5D0*(ALAT1+ALAT2)/(RLAT2-RLAT1)
      ASUM=ASUM+ALATJ*TAULAT(J)
      TSUM=TSUM+ALATJ*TAULAT(J)*SIZLAT(J)
      RLAT1=RLAT2
      ALAT1=ALAT2
  110 CONTINUE
      SIZ=(1.D-20+TSUM)/(1.D-10+ASUM)
      RETURN
      END SUBROUTINE BOXAV2
      END MODULE RE001


      SUBROUTINE GTREND(XNOW,TNOW)
C
      USE RE001, only: nghg,ghgyr1,ghgyr2,ghgam
!alt  integer nghg,nyrsghg,ghgyr1,ghgyr2        ! w/o module RE001
!alt  parameter (nyrsghg=2050-1850+1,nghg=6)
!alt  common/ghgcom/ghgyr1,ghgyr2,ghgam(nghg,nyrsghg)
      real*8 xnow(nghg),tnow,year,dy,frac
      integer iy,n
C
C-------------------------------------------------------------
C        Makiko's GHG Trend Compilation  GHG.1850-2050.Dec1999
C
C        Annual-Mean      Greenhouse Gas Mixing Ratios
C-------------------------------------------------------------
C                 CO2     N2O     CH4   CFC-11  CFC-12  others
C        Year     ppm     ppm     ppm     ppb     ppb     ppb
C-------------------------------------------------------------
C     Read from external file - outside table: use value from
C                                      years ghgyr1 or ghgyr2
      YEAR=TNOW
      IF(TNOW.LE.ghgyr1+.5D0) YEAR=ghgyr1+.5D0
      IF(TNOW.gE.ghgyr2+.49999D0) YEAR=ghgyr2+.49999D0
      DY=YEAR-(ghgyr1+.5D0)
      IY=DY
      frac=DY-IY
      IY=IY+1
C
C     CO2 N2O CH4 CFC-11 CFC-12 other_GHG  SCENARIO
C--------------------------------------------------
C
      do n=1,nghg
        XNOW(N)=GHGAM(N,IY)+frac*(GHGAM(N,IY+1)-GHGAM(N,IY))
      end do
C
      RETURN
      END SUBROUTINE GTREND

      SUBROUTINE PHATMO(P,H,D,T,O,Q,S,OCM,WCM,NPHD,NATM)
C     ------------------------------------------------------------------
C     -------------     MCCLATCHY (1972) ATMOSPHERE DATA     -----------
C     ------------------------------------------------------------------
C
C        INPUT DATA
C------------------
C                  NATM=0  GIVES ABREVIATED DATA FOR  STANDARD ATMOSPHER
C                                (INPUT: P OR H) (RETURNS: H OR P & D,T)
C
C                  NATM=1  GIVES ATMOSPHERE DATA FOR  TROPICAL LATITUDES
C                  NATM=2  GIVES ATMOSPHERE DATA FOR  MIDLATITUDE SUMMER
C                  NATM=3  GIVES ATMOSPHERE DATA FOR  MIDLATITUDE WINTER
C                  NATM=4  GIVES ATMOSPHERE DATA FOR  SUBARCTIC SUMMER
C                  NATM=5  GIVES ATMOSPHERE DATA FOR  SUBARCTIC WINTER
C                  NATM=6  GIVES ATMOSPHERE DATA FOR  STANDARD ATMOSPHER
C
C                  NPHD=1  RETURNS H,D,T,O,Q,S DATA FOR GIVEN PRESSURE P
C                  NPHD=2  RETURNS P,D,T,O,Q,S DATA FOR GIVEN   HEIGHT H
C                  NPHD=3  RETURNS P,H,T,O,Q,S DATA FOR GIVEN  DENSITY D
C
C       OUTPUT DATA
C------------------
C                  P = PRESSURE IN MILLIBARS
C                  H = HEIGHT IN KILOMETERS
C                  D = DENSITY IN GRAMS/METER**3
C                  T = TEMPERATURE (ABSOLUTE)
C                  O = OZONE MIXING RATIO (GRAMS OZONE)/(GRAMS AIR)
C                  Q = SPECIFIC HUMIDITY (GRAMS WATER VAPOR)/(GRAMS AIR)
C                  S = SATURATION RATIO (GRAMS WATER VAPOR)/(GRAMS AIR)
C                  OCM = OZONE (CM-STP) ABOVE GIVEN HEIGHT
C                  WCM = WATER VAPOR (CM-STP) ABOVE GIVEN HEIGHT
C
C           REMARKS
C------------------
C                  INPUT P,H,D PARAMETERS ARE NOT ALTERED
C                  P,D INTERPOLATION IS EXPONENTIAL WITH HEIGHT
C                  NO EXTRAPOLATION IS MADE OUTSIDE 0-100 KM INTERVAL
C                  S  IS NOT COMPUTED ABOVE 40 KM (FORMULA NOT ACCURATE)
C
C                  R = Q/S          GIVES RELATIVE HUMIDITY
C                  W = Q/(1-Q)      GIVES WATER VAPOR MIXING RATIO
C                  N = D*2.079E 16  GIVES NUMBER DENSITY PER CM**3
C

      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION    PRS1(33),PRS2(33),PRS3(33),PRS4(33),PRS5(33),PRS6(33)
     1            ,DNS1(33),DNS2(33),DNS3(33),DNS4(33),DNS5(33),DNS6(33)
     2            ,TMP1(33),TMP2(33),TMP3(33),TMP4(33),TMP5(33),TMP6(33)
     3            ,WVP1(33),WVP2(33),WVP3(33),WVP4(33),WVP5(33),WVP6(33)
     4            ,OZO1(33),OZO2(33),OZO3(33),OZO4(33),OZO5(33),OZO6(33)
      DIMENSION   PRES(33,6),DENS(33,6),TEMP(33,6),WVAP(33,6),OZON(33,6)

      EQUIVALENCE
     +     (PRES(1,1),PRS1(1)),(DENS(1,1),DNS1(1)),(TEMP(1,1),TMP1(1))
     +    ,(PRES(1,2),PRS2(1)),(DENS(1,2),DNS2(1)),(TEMP(1,2),TMP2(1))
     +    ,(PRES(1,3),PRS3(1)),(DENS(1,3),DNS3(1)),(TEMP(1,3),TMP3(1))
     +    ,(PRES(1,4),PRS4(1)),(DENS(1,4),DNS4(1)),(TEMP(1,4),TMP4(1))
     +    ,(PRES(1,5),PRS5(1)),(DENS(1,5),DNS5(1)),(TEMP(1,5),TMP5(1))
     +    ,(PRES(1,6),PRS6(1)),(DENS(1,6),DNS6(1)),(TEMP(1,6),TMP6(1))
      EQUIVALENCE (WVAP(1,1),WVP1(1)),(OZON(1,1),OZO1(1))
      EQUIVALENCE (WVAP(1,2),WVP2(1)),(OZON(1,2),OZO2(1))
      EQUIVALENCE (WVAP(1,3),WVP3(1)),(OZON(1,3),OZO3(1))
      EQUIVALENCE (WVAP(1,4),WVP4(1)),(OZON(1,4),OZO4(1))
      EQUIVALENCE (WVAP(1,5),WVP5(1)),(OZON(1,5),OZO5(1))
      EQUIVALENCE (WVAP(1,6),WVP6(1)),(OZON(1,6),OZO6(1))

      DIMENSION HTKM(33)
      DATA HTKM/1.0E-09, 1., 2., 3., 4., 5., 6., 7., 8., 9.,10.,11.
     1         ,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.
     2         ,25.,30.,35.,40.,45.,50.,70.,99.9/


C----------------------------------------------------------------------
C0000 GLOBAL   U.S. (1976) STANDARD ATMOSPHERE   P, T, GEO H  PARAMETERS
C----------------------------------------------------------------------

      DIMENSION SPLB(8),STLB(8),SHLB(8),SDLB(8)
      DATA SPLB/1013.25,226.32,54.748,8.6801,1.109,.66938,.039564
     +         ,3.7338E-03/
      DATA STLB/288.15,216.65,216.65,228.65,270.65,270.65,214.65,186.87/
      DATA SHLB/0.0,11.0,20.0,32.0,47.0,51.0,71.0,84.852/
      DATA SDLB/-6.5,0.0,1.0,2.8,0.0,-2.8,-2.0,0.0/
      DATA HPCON/34.16319/


C-----------------------------------------------------------------------
C1111 TROPICAL LATITUDES      MCCLATCHY (1972) ATMOSPHERE DATA VS HEIGHT
C-----------------------------------------------------------------------

      DATA PRS1/      1.013E 03,9.040E 02,8.050E 02,7.150E 02,6.330E 02,
     1      5.590E 02,4.920E 02,4.320E 02,3.780E 02,3.290E 02,2.860E 02,
     2      2.470E 02,2.130E 02,1.820E 02,1.560E 02,1.320E 02,1.110E 02,
     3      9.370E 01,7.890E 01,6.660E 01,5.650E 01,4.800E 01,4.090E 01,
     4      3.500E 01,3.000E 01,2.570E 01,1.220E 01,6.000E 00,3.050E 00,
     5      1.590E 00,8.540E-01,5.790E-02,3.000E-04/
      DATA DNS1/      1.167E 03,1.064E 03,9.689E 02,8.756E 02,7.951E 02,
     1      7.199E 02,6.501E 02,5.855E 02,5.258E 02,4.708E 02,4.202E 02,
     2      3.740E 02,3.316E 02,2.929E 02,2.578E 02,2.260E 02,1.972E 02,
     3      1.676E 02,1.382E 02,1.145E 02,9.515E 01,7.938E 01,6.645E 01,
     4      5.618E 01,4.763E 01,4.045E 01,1.831E 01,8.600E 00,4.181E 00,
     5      2.097E 00,1.101E 00,9.210E-02,5.000E-04/
      DATA TMP1/  300.0,294.0,288.0,284.0,277.0,270.0,264.0,257.0,250.0,
     1244.0,237.0,230.0,224.0,217.0,210.0,204.0,197.0,195.0,199.0,203.0,
     2207.0,211.0,215.0,217.0,219.0,221.0,232.0,243.0,254.0,265.0,270.0,
     3  219.0,210.0/
      DATA WVP1/1.9E 01,1.3E 01,9.3E 00,4.7E 00,2.2E 00,1.5E 00,8.5E-01,
     1  4.7E-01,2.5E-01,1.2E-01,5.0E-02,1.7E-02,6.0E-03,1.8E-03,1.0E-03,
     2  7.6E-04,6.4E-04,5.6E-04,5.0E-04,4.9E-04,4.5E-04,5.1E-04,5.1E-04,
     3  5.4E-04,6.0E-04,6.7E-04,3.6E-04,1.1E-04,4.3E-05,1.9E-05,6.3E-06,
     4  1.4E-07,1.0E-09/
      DATA OZO1/5.6E-05,5.6E-05,5.4E-05,5.1E-05,4.7E-05,4.5E-05,4.3E-05,
     1  4.1E-05,3.9E-05,3.9E-05,3.9E-05,4.1E-05,4.3E-05,4.5E-05,4.5E-05,
     2  4.7E-05,4.7E-05,6.9E-05,9.0E-05,1.4E-04,1.9E-04,2.4E-04,2.8E-04,
     3  3.2E-04,3.4E-04,3.4E-04,2.4E-04,9.2E-05,4.1E-05,1.3E-05,4.3E-06,
     4  8.6E-08,4.3E-11/

C-----------------------------------------------------------------------
C2222 MIDLATITUDE SUMMER      MCCLATCHY (1972) ATMOSPHERE DATA VS HEIGHT
C-----------------------------------------------------------------------

      DATA PRS2/      1.013E 03,9.020E 02,8.020E 02,7.100E 02,6.280E 02,
     1      5.540E 02,4.870E 02,4.260E 02,3.720E 02,3.240E 02,2.810E 02,
     2      2.430E 02,2.090E 02,1.790E 02,1.530E 02,1.300E 02,1.110E 02,
     3      9.500E 01,8.120E 01,6.950E 01,5.950E 01,5.100E 01,4.370E 01,
     4      3.760E 01,3.220E 01,2.770E 01,1.320E 01,6.520E 00,3.330E 00,
     5      1.760E 00,9.510E-01,6.710E-02,3.000E-04/
      DATA DNS2/      1.191E 03,1.080E 03,9.757E 02,8.846E 02,7.998E 02,
     1      7.211E 02,6.487E 02,5.830E 02,5.225E 02,4.669E 02,4.159E 02,
     2      3.693E 02,3.269E 02,2.882E 02,2.464E 02,2.104E 02,1.797E 02,
     3      1.535E 02,1.305E 02,1.110E 02,9.453E 01,8.056E 01,6.872E 01,
     4      5.867E 01,5.014E 01,4.288E 01,1.322E 01,6.519E 00,3.330E 00,
     5      1.757E 00,9.512E-01,6.706E-02,5.000E-04/
      DATA TMP2/  294.0,290.0,285.0,279.0,273.0,267.0,261.0,255.0,248.0,
     1242.0,235.0,229.0,222.0,216.0,216.0,216.0,216.0,216.0,216.0,217.0,
     2218.0,219.0,220.0,222.0,223.0,224.0,234.0,245.0,258.0,270.0,276.0,
     3  218.0,210.0/
      DATA WVP2/1.4E 01,9.3E 00,5.9E 00,3.3E 00,1.9E 00,1.0E 00,6.1E-01,
     1  3.7E-01,2.1E-01,1.2E-01,6.4E-02,2.2E-02,6.0E-03,1.8E-03,1.0E-03,
     2  7.6E-04,6.4E-04,5.6E-04,5.0E-04,4.9E-04,4.5E-04,5.1E-04,5.1E-04,
     3  5.4E-04,6.0E-04,6.7E-04,3.6E-04,1.1E-04,4.3E-05,1.9E-05,6.3E-06,
     4  1.4E-07,1.0E-09/
      DATA OZO2/6.0E-05,6.0E-05,6.0E-05,6.2E-05,6.4E-05,6.6E-05,6.9E-05,
     1  7.5E-05,7.9E-05,8.6E-05,9.0E-05,1.1E-04,1.2E-04,1.5E-04,1.8E-04,
     2  1.9E-04,2.1E-04,2.4E-04,2.8E-04,3.2E-04,3.4E-04,3.6E-04,3.6E-04,
     3  3.4E-04,3.2E-04,3.0E-04,2.0E-04,9.2E-05,4.1E-05,1.3E-05,4.3E-06,
     4  8.6E-08,4.3E-11/

C-----------------------------------------------------------------------
C3333 MIDLATITUDE WINTER      MCCLATCHY (1972) ATMOSPHERE DATA VS HEIGHT
C-----------------------------------------------------------------------

      DATA PRS3/      1.018E 03,8.973E 02,7.897E 02,6.938E 02,6.081E 02,
     1      5.313E 02,4.627E 02,4.016E 02,3.473E 02,2.992E 02,2.568E 02,
     2      2.199E 02,1.882E 02,1.610E 02,1.378E 02,1.178E 02,1.007E 02,
     3      8.610E 01,7.350E 01,6.280E 01,5.370E 01,4.580E 01,3.910E 01,
     4      3.340E 01,2.860E 01,2.430E 01,1.110E 01,5.180E 00,2.530E 00,
     5      1.290E 00,6.820E-01,4.670E-02,3.000E-04/
      DATA DNS3/      1.301E 03,1.162E 03,1.037E 03,9.230E 02,8.282E 02,
     1      7.411E 02,6.614E 02,5.886E 02,5.222E 02,4.619E 02,4.072E 02,
     2      3.496E 02,2.999E 02,2.572E 02,2.206E 02,1.890E 02,1.620E 02,
     3      1.388E 02,1.188E 02,1.017E 02,8.690E 01,7.421E 01,6.338E 01,
     4      5.415E 01,4.624E 01,3.950E 01,1.783E 01,7.924E 00,3.625E 00,
     5      1.741E 00,8.954E-01,7.051E-02,5.000E-04/
      DATA TMP3/  272.2,268.7,265.2,261.7,255.7,249.7,243.7,237.7,231.7,
     1225.7,219.7,219.2,218.7,218.2,217.7,217.2,216.7,216.2,215.7,215.2,
     2215.2,215.2,215.2,215.2,215.2,215.2,217.4,227.8,243.2,258.5,265.7,
     3  230.7,210.2/
      DATA WVP3/3.5E 00,2.5E 00,1.8E 00,1.2E 00,6.6E-01,3.8E-01,2.1E-01,
     1  8.5E-02,3.5E-02,1.6E-02,7.5E-03,6.9E-03,6.0E-03,1.8E-03,1.0E-03,
     2  7.6E-04,6.4E-04,5.6E-04,5.0E-04,4.9E-04,4.5E-04,5.1E-04,5.1E-04,
     3  5.4E-04,6.0E-04,6.7E-04,3.6E-04,1.1E-04,4.3E-05,1.9E-05,6.3E-06,
     4  1.4E-07,1.0E-09/
      DATA OZO3/6.0E-05,5.4E-05,4.9E-05,4.9E-05,4.9E-05,5.8E-05,6.4E-05,
     1  7.7E-05,9.0E-05,1.2E-04,1.6E-04,2.1E-04,2.6E-04,3.0E-04,3.2E-04,
     2  3.4E-04,3.6E-04,3.9E-04,4.1E-04,4.3E-04,4.5E-04,4.3E-04,4.3E-04,
     3  3.9E-04,3.6E-04,3.4E-04,1.9E-04,9.2E-05,4.1E-05,1.3E-05,4.3E-06,
     4  8.6E-08,4.3E-11/

C-----------------------------------------------------------------------
C4444 SUBARCTIC SUMMER        MCCLATCHY (1972) ATMOSPHERE DATA VS HEIGHT
C-----------------------------------------------------------------------

      DATA PRS4/      1.010E 03,8.960E 02,7.929E 02,7.000E 02,6.160E 02,
     1      5.410E 02,4.730E 02,4.130E 02,3.590E 02,3.107E 02,2.677E 02,
     2      2.300E 02,1.977E 02,1.700E 02,1.460E 02,1.250E 02,1.080E 02,
     3      9.280E 01,7.980E 01,6.860E 01,5.890E 01,5.070E 01,4.360E 01,
     4      3.750E 01,3.227E 01,2.780E 01,1.340E 01,6.610E 00,3.400E 00,
     5      1.810E 00,9.870E-01,7.070E-02,3.000E-04/
      DATA DNS4/      1.220E 03,1.110E 03,9.971E 02,8.985E 02,8.077E 02,
     1      7.244E 02,6.519E 02,5.849E 02,5.231E 02,4.663E 02,4.142E 02,
     2      3.559E 02,3.059E 02,2.630E 02,2.260E 02,1.943E 02,1.671E 02,
     3      1.436E 02,1.235E 02,1.062E 02,9.128E 01,7.849E 01,6.750E 01,
     4      5.805E 01,4.963E 01,4.247E 01,1.338E 01,6.614E 00,3.404E 00,
     5      1.817E 00,9.868E-01,7.071E-02,5.000E-04/
      DATA TMP4/  287.0,282.0,276.0,271.0,266.0,260.0,253.0,246.0,239.0,
     1232.0,225.0,225.0,225.0,225.0,225.0,225.0,225.0,225.0,225.0,225.0,
     2225.0,225.0,225.0,225.0,226.0,228.0,235.0,247.0,262.0,274.0,277.0,
     3  216.0,210.0/
      DATA WVP4/9.1E 00,6.0E 00,4.2E 00,2.7E 00,1.7E 00,1.0E 00,5.4E-01,
     1  2.9E-01,1.3E-02,4.2E-02,1.5E-02,9.4E-03,6.0E-03,1.8E-03,1.0E-03,
     2  7.6E-04,6.4E-04,5.6E-04,5.0E-04,4.9E-04,4.5E-04,5.1E-04,5.1E-04,
     3  5.4E-04,6.0E-04,6.7E-04,3.6E-04,1.1E-04,4.3E-05,1.9E-05,6.3E-06,
     4  1.4E-07,1.0E-09/
      DATA OZO4/4.9E-05,5.4E-05,5.6E-05,5.8E-05,6.0E-05,6.4E-05,7.1E-05,
     1  7.5E-05,7.9E-05,1.1E-04,1.3E-04,1.8E-04,2.1E-04,2.6E-04,2.8E-04,
     2  3.2E-04,3.4E-04,3.9E-04,4.1E-04,4.1E-04,3.9E-04,3.6E-04,3.2E-04,
     3  3.0E-04,2.8E-04,2.6E-04,1.4E-04,9.2E-05,4.1E-05,1.3E-05,4.3E-06,
     4  8.6E-08,4.3E-11/

C-----------------------------------------------------------------------
C5555 SUBARCTIC WINTER        MCCLATCHY (1972) ATMOSPHERE DATA VS HEIGHT
C-----------------------------------------------------------------------

      DATA PRS5/      1.013E 03,8.878E 02,7.775E 02,6.798E 02,5.932E 02,
     1      5.158E 02,4.467E 02,3.853E 02,3.308E 02,2.829E 02,2.418E 02,
     2      2.067E 02,1.766E 02,1.510E 02,1.291E 02,1.103E 02,9.431E 01,
     3      8.058E 01,6.882E 01,5.875E 01,5.014E 01,4.277E 01,3.647E 01,
     4      3.109E 01,2.649E 01,2.256E 01,1.020E 01,4.701E 00,2.243E 00,
     5      1.113E 00,5.719E-01,4.016E-02,3.000E-04/
      DATA DNS5/      1.372E 03,1.193E 03,1.058E 03,9.366E 02,8.339E 02,
     1      7.457E 02,6.646E 02,5.904E 02,5.226E 02,4.538E 02,3.879E 02,
     2      3.315E 02,2.834E 02,2.422E 02,2.071E 02,1.770E 02,1.517E 02,
     3      1.300E 02,1.113E 02,9.529E 01,8.155E 01,6.976E 01,5.966E 01,
     4      5.100E 01,4.358E 01,3.722E 01,1.645E 01,7.368E 00,3.330E 00,
     5      1.569E 00,7.682E-01,5.695E-02,5.000E-04/
      DATA TMP5/  257.1,259.1,255.9,252.7,247.7,240.9,234.1,227.3,220.6,
     1217.2,217.2,217.2,217.2,217.2,217.2,217.2,216.6,216.0,215.4,214.8,
     2214.1,213.6,213.0,212.4,211.8,211.2,216.0,222.2,234.7,247.0,259.3,
     3  245.7,210.0/
      DATA WVP5/1.2E 00,1.2E 00,9.4E-01,6.8E-01,4.1E-01,2.0E-01,9.8E-02,
     1  5.4E-02,1.1E-02,8.4E-03,5.5E-03,3.8E-03,2.6E-03,1.8E-03,1.0E-03,
     2  7.6E-04,6.4E-04,5.6E-04,5.0E-04,4.9E-04,4.5E-04,5.1E-04,5.1E-04,
     3  5.4E-04,6.0E-04,6.7E-04,3.6E-04,1.1E-04,4.3E-05,1.9E-05,6.3E-06,
     4  1.4E-07,1.0E-09/
      DATA OZO5/4.1E-05,4.1E-05,4.1E-05,4.3E-05,4.5E-05,4.7E-05,4.9E-05,
     1  7.1E-05,9.0E-05,1.6E-04,2.4E-04,3.2E-04,4.3E-04,4.7E-04,4.9E-04,
     2  5.6E-04,6.2E-04,6.2E-04,6.2E-04,6.0E-04,5.6E-04,5.1E-04,4.7E-04,
     3  4.3E-04,3.6E-04,3.2E-04,1.5E-04,9.2E-05,4.1E-05,1.3E-05,4.3E-06,
     4  8.6E-08,4.3E-11/

C----------------------------------------------------------------------
C6666 GLOBAL   U.S. (1976) STANDARD ATMOSPHERE   P, T, GEO H  PARAMETERS
C----------------------------------------------------------------------

      DATA PRS6/    1.01325E+03,8.987E+02,7.950E+02,7.011E+02,6.164E+02,
     1      5.402E+02,4.718E+02,4.106E+02,3.560E+02,3.074E+02,2.644E+02,
     2      2.263E+02,1.933E+02,1.651E+02,1.410E+02,1.204E+02,1.029E+02,
     3      8.787E+01,7.505E+01,6.410E+01,5.475E+01,4.678E+01,4.000E+01,
     4      3.422E+01,2.931E+01,2.511E+01,1.172E+01,5.589E+00,2.775E+00,
     5      1.431E+00,7.594E-01,4.634E-02,2.384E-04/
      DATA DNS6/      1.225E+03,1.112E+03,1.006E+03,9.091E+02,8.191E+02,
     1      7.361E+02,6.597E+02,5.895E+02,5.252E+02,4.663E+02,4.127E+02,
     2      3.639E+02,3.108E+02,2.655E+02,2.268E+02,1.937E+02,1.654E+02,
     3      1.413E+02,1.207E+02,1.031E+02,8.803E+01,7.487E+01,6.373E+01,
     4      5.428E+01,4.627E+01,3.947E+01,1.801E+01,8.214E+00,3.851E+00,
     5      1.881E+00,9.775E-01,7.424E-02,4.445E-04/
      DATA TMP6/
     1         288.150,281.650,275.150,268.650,262.150,255.650,249.150,
     2         242.650,236.150,229.650,223.150,216.650,216.650,216.650,
     3         216.650,216.650,216.650,216.650,216.650,216.650,216.650,
     4         217.650,218.650,219.650,220.650,221.650,226.650,237.050,
     5         251.050,265.050,270.650,217.450,186.870/
      DATA WVP6/      1.083E+01,6.323E+00,3.612E+00,2.015E+00,1.095E+00,
     1      5.786E-01,2.965E-01,1.469E-01,7.021E-02,3.226E-02,1.419E-02,
     2      5.956E-03,5.002E-03,4.186E-03,3.490E-03,2.896E-03,2.388E-03,
     3      1.954E-03,1.583E-03,1.267E-03,9.967E-04,8.557E-04,7.104E-04,
     4      5.600E-04,4.037E-04,2.406E-04,5.404E-05,2.464E-05,1.155E-05,
     5      5.644E-06,2.932E-06,2.227E-07,1.334E-09/
      DATA OZO6/      7.526E-05,3.781E-05,6.203E-05,3.417E-05,5.694E-05,
     1      3.759E-05,5.970E-05,4.841E-05,7.102E-05,6.784E-05,9.237E-05,
     2      9.768E-05,1.251E-04,1.399E-04,1.715E-04,1.946E-04,2.300E-04,
     3      2.585E-04,2.943E-04,3.224E-04,3.519E-04,3.714E-04,3.868E-04,
     4      3.904E-04,3.872E-04,3.728E-04,2.344E-04,9.932E-05,3.677E-05,
     5      1.227E-05,4.324E-06,5.294E-08,1.262E-10/


      IF(NATM.GT.0) GO TO 200
      O=1.E-10
      Q=1.E-10
      S=1.E-10
      OCM=1.E-10
      WCM=1.E-10
      IF(NPHD.LT.2) GO TO 150
      DO 110 N=2,8
      IF(H.LT.SHLB(N)) GO TO 120
  110 CONTINUE
      N=9
  120 CONTINUE
      N=N-1
      IF(ABS(SDLB(N)).LT.1.E-04) GO TO 130
      P=SPLB(N)*(1.+SDLB(N)/STLB(N)*(H-SHLB(N)))**(-HPCON/SDLB(N))
      GO TO 140
  130 CONTINUE
      P=SPLB(N)*EXP(-HPCON/STLB(N)*(H-SHLB(N)))
  140 CONTINUE
      T=STLB(N)+SDLB(N)*(H-SHLB(N))
      D=P/T*28.9644E 05/8.31432E 03
      RETURN

  150 CONTINUE
      DO 160 N=2,8
      IF(P.GT.SPLB(N)) GO TO 170
  160 CONTINUE
      N=9
  170 CONTINUE
      N=N-1
      IF(ABS(SDLB(N)).LT.1.E-04) GO TO 180
      H=SHLB(N)+STLB(N)/SDLB(N)*((SPLB(N)/P)**(SDLB(N)/HPCON)-1.)
      GO TO 190
  180 CONTINUE
      H=SHLB(N)+STLB(N)/HPCON*LOG(SPLB(N)/P)
  190 CONTINUE
      T=STLB(N)+SDLB(N)*(H-SHLB(N))
      D=P/T*28.9644E 05/8.31432E 03
      RETURN

 200  CONTINUE
      IF(NPHD.EQ.1) GO TO 240
      IF(NPHD.EQ.2) GO TO 220
      XX=D
      XI=DENS(1,NATM)
      IF(D.GT.XI) XX=XI
      IF(D.LT.5.0E-04) GO TO 280
      DO 210 J=2,33
      XJ=DENS(J,NATM)
      IF(XX.GT.XJ) GO TO 260
      XI=XJ
  210 CONTINUE
  220 CONTINUE
      XX=H
      XI=HTKM(1)
      IF(H.LT.XI) XX=XI
      IF(H.GT.99.9) GO TO 280
      DO 230 J=2,33
      XJ=HTKM(J)
      IF(XX.LT.XJ) GO TO 260
      XI=XJ
  230 CONTINUE
  240 CONTINUE
      XX=P
      XI=PRES(1,NATM)
      IF(P.GT.XI) XX=XI
      IF(P.LT.3.0E-04) GO TO 280
      DO 250 J=2,33
      XJ=PRES(J,NATM)
      IF(XX.GT.XJ) GO TO 260
      XI=XJ
  250 CONTINUE
  260 CONTINUE
      DELTA=(XX-XI)/(XJ-XI)
      I=J-1
      IF(NPHD.NE.2) THEN
      RAT=LOG(XX/XI)/LOG(XJ/XI)
      H=HTKM(I)+(HTKM(J)-HTKM(I))*RAT
      T=TEMP(I,NATM)+(TEMP(J,NATM)-TEMP(I,NATM))*RAT
      ENDIF
      PI=PRES(I,NATM)
      PJ=PRES(J,NATM)
      DI=DENS(I,NATM)
      DJ=DENS(J,NATM)
      IF(NPHD.EQ.1) D=DI+DELTA*(DJ-DI)
      IF(NPHD.EQ.3) P=PI+DELTA*(PJ-PI)
      IF(NPHD.EQ.2) THEN
      P=PI*(PJ/PI)**DELTA
      D=DI*(DJ/DI)**DELTA
      T=TEMP(I,NATM)+DELTA*(TEMP(J,NATM)-TEMP(I,NATM))
      ENDIF
      O=OZON(I,NATM)/DI+DELTA*(OZON(J,NATM)/DJ-OZON(I,NATM)/DI)
      Q=WVAP(I,NATM)/DI+DELTA*(WVAP(J,NATM)/DJ-WVAP(I,NATM)/DI)
      ES=10.D0**(9.4051D0-2353.D0/T)
      IF(P.LT.PI) PI=P
      S=1.D+06
      RS=(PI-ES+0.622*ES)/(0.622*ES)
      IF(RS.GT.1.E-06) S=1./RS
      OI=O
      QI=Q
      OCM=0.D0
      WCM=0.D0
      DO 270 K=J,33
      PJ=PRES(K,NATM)
      DJ=DENS(K,NATM)
      OJ=OZON(K,NATM)/DJ
      QJ=WVAP(K,NATM)/DJ
      DP=PI-PJ
      OCM=OCM+0.5D0*(OI+OJ)*DP
      WCM=WCM+0.5D0*(QI+QJ)*DP
      OI=OJ
      QI=QJ
      PI=PJ
  270 CONTINUE
      WCM=WCM/0.980D0*22420.7D0/18.D0
      OCM=OCM/0.980D0*22420.7D0/48.D0
      RETURN
  280 CONTINUE
      T=210.D0
      IF(NATM.EQ.6) T=186.87
      O=1.D-10
      Q=1.D-10
      S=1.D-10
      OCM=1.D-10
      WCM=1.D-10
      IF(NPHD.NE.1) P=1.D-05
      IF(NPHD.NE.2) H=99.99
      IF(NPHD.NE.3) D=2.D-05
      RETURN
      END SUBROUTINE PHATMO

      FUNCTION PFOFTK(WAVNA,WAVNB,TK)
C     ------------------------------------------------------------------
C
C        INPUT DATA
C                  WAVNA,WAVNB  SPECLTRAL INTERVAL IN WAVENUMBERS
C                               (ORDER OF WAVNA,WAVNB NOT IMPORTANT)
C
C                  TK           ABSOLUTE TEMPERATURE IN DEGREES KELVIN
C
C       OUTPUT DATA
C                  PFOFTK       PLANCK FLUX (W/M**2)
C
C
C           REMARKS
C                   PLANCK INTENSITY (W/M**2/STER) IS GIVEN BY PFOFTK/PI
C
C     ------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION BN(21),BD(21)
      DATA BN/1.D0,-1.D0,1.D0,-1.D0,1.D0,-1.D0,5.D0,-691.D0,7.D0
     1,-3617.D0,43867.D0,-174611.D0,854513.D0,-236364091.D0
     2,8553103.D0,-23749461029.D0,8615841276005.D0,-7709321041217.D0
     3,2577687858367.D0,-2631527155305348D 04,2929993913841559D0/
      DATA BD/1.D0,2.D0,6.D0,30.D0,42.D0,30.D0,66.D0,2730.D0,6.D0
     1,510.D0,798.D0,330.D0,138.D0,2730.D0,6.D0,870.D0,14322.D0
     2,510.D0,6.D0,1919190.D0,6.D0/
      DATA PI4/97.40909103400244D0/
C     DATA  PI/3.141592653589793D0/
      DATA HCK/1.43879D0/
      DATA DGXLIM/1.D-06/
      PFOFTK=0.D0
      IF(TK.LT.1.D-06) RETURN
      DO 160 II=1,2
      IF(II.EQ.1) X=HCK*WAVNA/TK
      IF(II.EQ.2) X=HCK*WAVNB/TK
      IF(X.GT.2.3D0) GO TO 120
      XX=X*X
      GSUM=1.D0/3.D0-X/8.D0+XX/60.D0
      NB=3
      XNF=XX/2.D0
      DO 100 N=4,38,2
      NB=NB+1
      NNB=NB
      B=BN(NB)/BD(NB)
      XN3=N+3
      XNM=N*(N-1)
      XNF=XNF*(XX/XNM)
      DG=B/XN3*XNF
      GSUM=GSUM+DG
      DGB=DG
      IF(ABS(DG).LT.DGXLIM) GO TO 110
  100 CONTINUE
  110 CONTINUE
      GX=GSUM*XX*X
      GO TO 150
  120 CONTINUE
      GSUM=PI4/15.D0
      DO 130 N=1,20
      NNB=N
      XN=N
      XNN=XN*XN
      XNX=XN*X
      IF(XNX.GT.100.D0) GO TO 140
      GTERM=(X*X*(3.D0+XNX)+6.D0*(1.D0+XNX)/XNN)/XNN
      DG=GTERM*EXP(-XNX)
      GSUM=GSUM-DG
      DGB=DG
      IF(DG.LT.DGXLIM) GO TO 140
  130 CONTINUE
  140 CONTINUE
      GX=GSUM
  150 CONTINUE
      IF(II.EQ.1) GXA=GX
      IF(II.EQ.2) GXB=GX
  160 CONTINUE
      PNORM=15.D0/PI4
      PFOFTK=ABS(GXB-GXA)*PNORM
      PFOFTK=PFOFTK*5.6692D-08*TK**4
      RETURN
      END FUNCTION PFOFTK

      FUNCTION TKOFPF(WAVNA,WAVNB,FLUXAB)
C     ------------------------------------------------------------------
C
C        INPUT DATA
C------------------
C                  WAVNA,WAVNB  SPECLTRAL INTERVAL IN WAVENUMBERS
C                               (ORDER OF WAVNA,WAVNB NOT IMPORTANT)
C                  FLUXAB       PLANCK FLUX (W/M**2) IN INTERVAL
C                                                       (WAVNA,WAVNB)
C
C       OUTPUT DATA
C------------------
C                  TK           BRIGHTNESS TEMPERATURE IN DEGREES KELVIN
C
C
C           REMARKS
C------------------
C                   TKOFPF IS INVERSE FUNCTION OF PFOFTK(WAVNA,WAVNB,TK)
C                   THE OUTPUT OF TKOFPF SATISFIES THE IDENTITY
C                                 FLUXAB=PFOFTK(WAVNA,WAVNB,TK)
C                   (UNITS FOR FLUXAB AND PFOFTK MUST BE IDENTICAL)
C
C     ------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LOGFIT
      DATA DELFIT/1.D-06/
      DATA NMAX/20/
      IF(FLUXAB.LE.0.D0) RETURN
      LOGFIT=.FALSE.
      NFIT=0
      PF=FLUXAB
      XA=0.D0
      YA=0.D0
      XB=250.D0
      YB=PFOFTK(WAVNA,WAVNB,XB)
      XX=PF*XB/YB
      YY=PFOFTK(WAVNA,WAVNB,XX)
      IF(ABS(YY-PF).LT.DELFIT) GO TO 200
      IF((YY/PF).LT.0.5D0) GO TO 150
      IF((YY/PF).GT.2.0D0) GO TO 170
      IF(XX.GT.XB) GO TO 110
      XC=XB
      YC=YB
      XB=XX
      YB=YY
      GO TO 120
  110 CONTINUE
      XC=XX
      YC=YY
  120 CONTINUE
      XBA=XB-XA
      XCA=XC-XA
      XBC=XB-XC
      YBA=YB-YA
      YCA=YC-YA
      YBC=YB-YC
      NFIT=NFIT+1
      IF(NFIT.GT.NMAX) GO TO 200
      YXBA=YBA/XBA
      YXCA=YCA/XCA
      C=(YXBA-YXCA)/XBC
      B=YXBA-(XB+XA)*C
      A=YA-XA*(B+XA*C)
      ROOT=SQRT(B*B+4.D0*C*(PF-A))
      XX=0.5D0*(ROOT-B)/C
      IF(XX.LT.XA.OR.XX.GT.XC) XX=-0.5D0*(ROOT+B)/C
      YY=PFOFTK(WAVNA,WAVNB,XX)
      IF(LOGFIT) YY=LOG(YY)
      IF(ABS(YY-PF).LT.DELFIT) GO TO 200
      IF(XX.GT.XB) GO TO 130
      XC=XB
      YC=YB
      GO TO 140
  130 CONTINUE
      XA=XB
      YA=YB
  140 CONTINUE
      XB=XX
      YB=YY
      GO TO 120
  150 CONTINUE
      XA=XX
      YA=YY
  160 CONTINUE
      XC=XB
      YC=YB
      XB=XB/2.D0
      YB=PFOFTK(WAVNA,WAVNB,XB)
      IF(YB.LT.YA) GO TO 190
      IF(YB.GT.PF) GO TO 160
      XA=XB
      YA=YB
      GO TO 190
  170 CONTINUE
      XC=XX
      YC=YY
  180 CONTINUE
      XA=XB
      YA=YB
      XB=XB*2.D0
      YB=PFOFTK(WAVNA,WAVNB,XB)
      IF(YB.GT.YC) GO TO 190
      IF(YB.LT.PF) GO TO 180
      XC=XB
      YC=YB
  190 CONTINUE
      XB=XA+(PF-YA)*(XC-XA)/(YC-YA)
      YB=PFOFTK(WAVNA,WAVNB,XB)
      XX=XB
      IF(ABS(YB-PF).LT.DELFIT) GO TO 200
      PF=LOG(PF)
      YA=LOG(YA)
      YB=LOG(YB)
      YC=LOG(YC)
      LOGFIT=.TRUE.
      GO TO 120
  200 CONTINUE
      TKOFPF=XX
      RETURN
      END FUNCTION TKOFPF

      SUBROUTINE SPLINE(X,F,NXF,XX,FF,CUSPWM,CUSPWE,KXTRAP)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(NXF),F(NXF)

C---------------------------------------------------------------------
C
C    SPLINE locates XX between points (F2,X2)(F3,X3) on 4-point spread
C       and returns 4-point Cubic Spline interpolated value FF = F(XX)
C
C    Quadratic Derivatives of Spline are continuous at (F2,X2),(F3,X3)
C    (X-Coordinate may be specified in increasing or decreasing order)
C
C---------------------------------------------------------------------
C
C    Custom Control Parameters:  CUSPWM,CUSPWE,KXTRAP
C------------------------------
C
C    In cases where data points are unevenly spaced and/or data points
C    exhibit abrupt changes in value, Spline Interpolation may produce
C    undesirable bulging of interpolated values. In more extreme cases
C    Linear Interpolation may be less problematic to use.
C
C    Interpolation can be weighted between: Cubic Spline and Linear by
C    adjusting weights CUSPWM and CUSPWE to values between 1.0 and 0.0
C
C    CUSPWM = Cubic Spline Weight at the (X2-X3) Interval Mid-point
C    CUSPWE = Cubic Spline Weight at the (X2-X3) Interval End-points
C
C    For example, with:
C
C    CUSPWM=1.0,CUSPWE=1.0  FF returns Cubic Spline interpolated value
C    CUSPWM=0.0,CUSPWE=0.0  FF returns   Linearly   interpolated value
C
C---------------------------------------------------------------------
C
C     Extrapolation for XX outside of defined interval:  X(1)<->X(NXF)
C
C               KXTRAP = 0    No Extrapolation  (i.e., sets F(XX)=0.0)
C                        1    Fixed Extrapolation (F(XX) = edge value)
C                        2    Linear Extrapolation using 2 edge points
C
C---------------------------------------------------------------------

      K=2
      X2=X(K)
      X3=X(NXF-1)
      BETW=(XX-X2)*(X3-XX)
      IF(BETW.LE.0.D0) GO TO 120

  100 CONTINUE
      K=K+1
      X3=X(K)
      BETW=(XX-X2)*(X3-XX)
      IF(BETW.GE.0.D0) GO TO 110
      X2=X3
      GO TO 100

  110 CONTINUE
      F3=F(K)
      F4=F(K+1)
      X4=X(K+1)
      F2=F(K-1)
      X2=X(K-1)
      F1=F(K-2)
      X1=X(K-2)
      X21=X2-X1
      X31=X3-X1
      X32=X3-X2
      X43=X4-X3
      X42=X4-X2
      F21=(F2-F1)/(X21*X21)
      F32=(F3-F2)/(X32*X32)
      F43=(F4-F3)/(X43*X43)
      F3221=(F32+F21)/X31*X21
      F4332=(F43+F32)/X42*X43
      A=F2
      B=X32*F3221
      C=3.D0*F32-F3221-F3221-F4332
      D=(F3221+F4332-F32-F32)/X32
      XF=XX-X2

C                             FFCUSP= Cubic Spline Interpolation Result
C                             -----------------------------------------

      FFCUSP=A+XF*(B+XF*(C+XF*D))
      XE=(X3+X2-XX-XX)/X32
      IF(XE.LT.0.D0) XE=-XE
      XEXM=XE**2
      CUSPWT=(1.D0-XEXM)*CUSPWM+XEXM*CUSPWE

C                                   FFLINR= Linear Interpolation Result
C                                   -----------------------------------
      FFLINR=A+XF*F32*X32
      FF=FFCUSP*CUSPWT+FFLINR*(1.D0-CUSPWT)
      GO TO 160

C                Edge Point Interval Interpolation and/or Extrapolation
C                ------------------------------------------------------
  120 CONTINUE
      BETW=(X2-XX)*(X3-X2)
      IF(BETW.LT.0.D0) GO TO 140

C                          X(1),X(2)  Edge Point Interval Interpolation
C                          --------------------------------------------
      X1=X(1)
      F1=F(1)
      F2=F(2)
      X21=X2-X1
      F21=(F2-F1)/X21
      XF=XX-X1
      BETW=(X2-XX)*XF
      IF(BETW.LT.0.D0) GO TO 130
      F3=F(3)
      X3=X(3)
      X32=X3-X2
      X31=X3-X1
      C=((F3-F2)/X32-F21)/X31
      B=F21-X21*C
      A=F1
      FFCUSP=A+XF*(B+XF*C)
      FFLINR=A+XF*F21
      XE=1.D0-2.D0*XF/X21
      IF(XE.LT.0.D0) XE=-XE
      XEXM=XE**2
      CUSPWT=(1.D0-XEXM)*CUSPWM+XEXM*CUSPWE
      FF=FFCUSP*CUSPWT+FFLINR*(1.D0-CUSPWT)
      GO TO 160

  130 CONTINUE
C                  Extrapolation for XX Outside of Interval X(1) - X(2)
C                  ----------------------------------------------------
C                  IF(KXTRAP.EQ.0)  (No Extrapolation:  sets F(XX)=0.0)
C                  IF(KXTRAP.EQ.1)  (Extrapolation at Fixed Edge Value)
C                  IF(KXTRAP.EQ.2)  (2 Edge Point Linear Extrapolation)

      IF(KXTRAP.EQ.0) FF=0.D0
      IF(KXTRAP.EQ.1) FF=F1
      IF(KXTRAP.EQ.2) FF=F1+XF*F21
      GO TO 160

  140 CONTINUE
C                    X(NXF-1),X(NXF)  Edge Point Interval Interpolation
C                    --------------------------------------------------
      F3=F(NXF)
      X3=X(NXF)
      F2=F(NXF-1)
      X2=X(NXF-1)
      X32=X3-X2
      F32=(F3-F2)/X32
      XF=XX-X3
      BETW=(X2-XX)*(XX-X3)
      IF(BETW.LT.0.D0) GO TO 150
      F1=F(NXF-2)
      X1=X(NXF-2)
      X21=X2-X1
      X31=X3-X1
      F21=(F2-F1)/X21
      XF=XX-X2

C                    3-Point Quadratic Interpolation for Edge Intervals
C                    --------------------------------------------------
C
C      (Edge Option)     ----------------------------------------------
C                        For Linear Interpolation within Edge Intervals
C                        between X(1),X(2), and between X(NXF-1),X(NXF)
C                        set the value of coefficient C below, to C=0.0
C                        ----------------------------------------------

      C=(F32-F21)/X31
      B=F21+X21*C
      A=F2
      FFCUSP=A+XF*(B+XF*C)
      FFLINR=A+XF*F32
      XE=1.D0-2.D0*XF/X32
      IF(XE.LT.0.D0) XE=-XE
      XEXM=XE**2
      CUSPWT=(1.D0-XEXM)*CUSPWM+XEXM*CUSPWE
      FF=FFCUSP*CUSPWT+FFLINR*(1.D0-CUSPWT)
      GO TO 160

  150 CONTINUE
C              Extrapolation for X Outside of Interval  X(NXF-1)-X(NXF)
C              --------------------------------------------------------
C                  IF(KXTRAP.EQ.0)  (No Extrapolation:  sets F(XX)=0.0)
C                  IF(KXTRAP.EQ.1)  (Extrapolation at Fixed Edge Value)
C                  IF(KXTRAP.EQ.2)  (2 Edge Point Linear Extrapolation)

      IF(KXTRAP.EQ.0) FF=0.D0
      IF(KXTRAP.EQ.1) FF=F3
      IF(KXTRAP.EQ.2) FF=F3+XF*(F3-F2)/(X3-X2)

  160 CONTINUE
      RETURN
      END SUBROUTINE SPLINE

      SUBROUTINE REPART(FXL,XLB,NXB,GYL,YLB,NYB)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION FXL(NXB),XLB(NXB),GYL(NYB),YLB(NYB)

C     ------------------------------------------------------------------
C
C     REPART   Repartitions FXL (a histogram-tyoe distribution function)
C              where XLB depicts the NXB partitions that define FXL data
C              FXL is assumed to be constant between XLB(N) AND XLB(N+1)
C
C              GYL(N) is the new histogram distribution function defined
C              by the (input) NYB partitions YLB(N).  The YLB partitions
C              can differ from XLB in number and/or spacing, or in upper
C              and lower limits.  XLB and YLB coordinates are assumed to
C              be linear both increasing or decreasing in the same sense
C
C     RETURNS  GYL as a histogram-type distribution function, with NYB-1
C              values assumed to be constant between YLB(N) and YLB(N+1)
C
C       NOTE:  The column amount of a vertically distributed quantity is
C              conserved, within the Repartition Interval that is common
C              to to XLB and YLB (layer bottom edge and top edge) limits
C
C     ------------------------------------------------------------------

      NXF=NXB-1
      NYG=NYB-1
      SUMG=0.D0
      DO 50 I=1,NYG
      GYL(I)=0.D0
   50 CONTINUE
      SUMY=0.D0
      I=1
      XA=XLB(I)
      J=1
      YA=YLB(J)
      XB=XLB(I+1)
      IF(XB.LT.XA) GO TO 200
  100 CONTINUE
      YB=YLB(J+1)
      IF(YB.GT.XA) GO TO 110
      GYL(J)=0.D0
      J=J+1
      IF(J.GT.NYG) GO TO 160
      YA=YB
      GO TO 100
  110 CONTINUE
      XB=XLB(I+1)
      IF(XB.GT.YA) GO TO 120
      I=I+1
      IF(I.GT.NXF) GO TO 160
      XA=XB
      GO TO 110
  120 CONTINUE
      XAYA=XA
      IF(YA.GT.XA) XAYA=YA
      IF(YB.GT.XB) GO TO 130
      PART=(YB-XAYA)/(XB-XA)
      SUMG=SUMG+PART*FXL(I)
      SUMY=SUMY+PART
      GYL(J)=SUMG
      J=J+1
      IF(J.GT.NYG) GO TO 160
      SUMG=0.D0
      SUMY=0.D0
      YA=YB
      YB=YLB(J+1)
      GO TO 120
  130 CONTINUE
      PART=(XB-XAYA)/(XB-XA)
      SUMG=SUMG+PART*FXL(I)
      SUMY=SUMY+PART
      I=I+1
      IF(I.GT.NXF) GO TO 140
      XA=XB
      XB=XLB(I+1)
      GO TO 120
  140 CONTINUE
      GYL(J)=SUMG
  150 CONTINUE
      J=J+1
      IF(J.GT.NYG) GO TO 160
      GYL(J)=0.D0
      GO TO 150
  160 CONTINUE
      GO TO 300

  200 CONTINUE
      YB=YLB(J+1)
      IF(YB.LT.XA) GO TO 210
      GYL(J)=0.D0
      J=J+1
      IF(J.GT.NYG) GO TO 260
      YA=YB
      GO TO 200
  210 CONTINUE
      XB=XLB(I+1)
      IF(XB.LT.YA) GO TO 220
      I=I+1
      IF(I.GT.NXF) GO TO 260
      XA=XB
      GO TO 210
  220 CONTINUE
      XAYA=XA
      IF(YA.LT.XA) XAYA=YA
      IF(YB.LT.XB) GO TO 230
      PART=(YB-XAYA)/(XB-XA)
      SUMG=SUMG+PART*FXL(I)
      SUMY=SUMY+PART
      GYL(J)=SUMG
      J=J+1
      IF(J.GT.NYG) GO TO 260
      SUMG=0.D0
      SUMY=0.D0
      YA=YB
      YB=YLB(J+1)
      GO TO 220
  230 CONTINUE
      PART=(XB-XAYA)/(XB-XA)
      SUMG=SUMG+PART*FXL(I)
      SUMY=SUMY+PART
      I=I+1
      IF(I.GT.NXF) GO TO 240
      XA=XB
      XB=XLB(I+1)
      GO TO 220
  240 CONTINUE
      GYL(J)=SUMG
  250 CONTINUE
      J=J+1
      IF(J.GT.NYG) GO TO 260
      GYL(J)=0.D0
      GO TO 250
  260 CONTINUE

  300 CONTINUE
      RETURN
      END SUBROUTINE REPART

      SUBROUTINE RETERP(FXL,XLB,NXB,GYL,YLB,NYB)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION FXL(NXB),XLB(NXB),GYL(NYB),YLB(NYB)

C     ------------------------------------------------------------------
C     RETERP
C              Interpolates FXL (a histogram-type distribution function)
C              where XLB depicts the NXB partitions that define FXL data
C              FXL is assumed to be constant between XLB(N) AND XLB(N+1)
C
C              GYL(N) is the new histogram distribution function defined
C              by the (input) NYB partitions YLB(N).  The YLB partitions
C              can differ from XLB in number and/or spacing, or in upper
C              and lower limits.  XLB and YLB coordinates are assumed to
C              be linear both increasing or decreasing in the same sense
C
C     RETURNS  GYL as a histogram-type distribution function, with NYB-1
C              values assumed to be constant between YLB(N) and YLB(N+1)
C
C     ------------------------------------------------------------------

      NXF=NXB-1
      NYG=NYB-1
      SUMG=0.D0
      DO 50 I=1,NYG
      GYL(I)=0.D0
   50 CONTINUE
      SUMY=0.D0
      I=1
      XA=XLB(I)
      J=1
      YA=YLB(J)
      XB=XLB(I+1)
      IF(XB.LT.XA) GO TO 200
  100 CONTINUE
      YB=YLB(J+1)
      IF(YB.GT.XA) GO TO 110
      GYL(J)=0.D0
      J=J+1
      IF(J.GT.NYG) GO TO 160
      YA=YB
      GO TO 100
  110 CONTINUE
      XB=XLB(I+1)
      IF(XB.GT.YA) GO TO 120
      I=I+1
      IF(I.GT.NXF) GO TO 160
      XA=XB
      GO TO 110
  120 CONTINUE
      XAYA=XA
      IF(YA.GT.XA) XAYA=YA
      IF(YB.GT.XB) GO TO 130
      PART=(YB-XAYA)/(XB-XA)
      SUMG=SUMG+PART*FXL(I)
      SUMY=SUMY+PART
      GYL(J)=SUMG/SUMY
      J=J+1
      IF(J.GT.NYG) GO TO 160
      SUMG=0.D0
      SUMY=0.D0
      YA=YB
      YB=YLB(J+1)
      GO TO 120
  130 CONTINUE
      PART=(XB-XAYA)/(XB-XA)
      SUMG=SUMG+PART*FXL(I)
      SUMY=SUMY+PART
      I=I+1
      IF(I.GT.NXF) GO TO 140
      XA=XB
      XB=XLB(I+1)
      GO TO 120
  140 CONTINUE
      GYL(J)=SUMG/SUMY
  150 CONTINUE
      J=J+1
      IF(J.GT.NYG) GO TO 160
      GYL(J)=0.D0
      GO TO 150
  160 CONTINUE
      GO TO 300

  200 CONTINUE
      YB=YLB(J+1)
      IF(YB.LT.XA) GO TO 210
      GYL(J)=0.D0
      J=J+1
      IF(J.GT.NYG) GO TO 260
      YA=YB
      GO TO 200
  210 CONTINUE
      XB=XLB(I+1)
      IF(XB.LT.YA) GO TO 220
      I=I+1
      IF(I.GT.NXF) GO TO 260
      XA=XB
      GO TO 210
  220 CONTINUE
      XAYA=XA
      IF(YA.LT.XA) XAYA=YA
      IF(YB.LT.XB) GO TO 230
      PART=(YB-XAYA)/(XB-XA)
      SUMG=SUMG+PART*FXL(I)
      SUMY=SUMY+PART
      GYL(J)=SUMG/SUMY
      J=J+1
      IF(J.GT.NYG) GO TO 260
      SUMG=0.D0
      SUMY=0.D0
      YA=YB
      YB=YLB(J+1)
      GO TO 220
  230 CONTINUE
      PART=(XB-XAYA)/(XB-XA)
      SUMG=SUMG+PART*FXL(I)
      SUMY=SUMY+PART
      I=I+1
      IF(I.GT.NXF) GO TO 240
      XA=XB
      XB=XLB(I+1)
      GO TO 220
  240 CONTINUE
      GYL(J)=SUMG/SUMY
  250 CONTINUE
      J=J+1
      IF(J.GT.NYG) GO TO 260
      GYL(J)=0.D0
      GO TO 250
  260 CONTINUE

  300 CONTINUE
      RETURN
      END SUBROUTINE RETERP

      SUBROUTINE FXGINT(F,X,NX,G,Y,NY,ALIM,BLIM,ABINT)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(NX),X(NX),G(NY),Y(NY)
CCC   DATA DELTA/1.D-07/
      PARAMETER (DELTA=1.D-07)

C     ------------------------------------------------------------------
C     FXGINT  PERFORMS NUMERICAL INTEGRATION (AREA UNDER CURVE) OF  F*G
C             BETWEEN THE LIMITS X=ALIM AND X=BLIM  (WITH BLIM GT ALIM)
C
C       F(X)  IS DEFINED BY CONNECTING SUCCESSIVE F(X) DATA POINTS USING
C             STRAIGHT-LINE SEGMENTS, I.E. F(X) IS PIECE-WISE CONTINUOUS
C             THE  X  COORDINATE CAN BE IN ASCENDING OR DESCENDING ORDER
C
C       G(Y)  IS DEFINED BY CONNECTING SUCCESSIVE G(Y) DATA POINTS USING
C             STRAIGHT-LINE SEGMENTS, I.E. G(Y) IS PIECE-WISE CONTINUOUS
C             THE  Y  COORDINATE CAN BE IN ASCENDING OR DESCENDING ORDER
C
C             (X,Y ARE THE SAME LINEAR COORDINATE INDEPENDENTLY DEFINED)
C
C             (F(X) IS ZERO OUTSIDE THE INTERVAL BETWEEN X(1) AND X(NX))
C             (G(Y) IS ZERO OUTSIDE THE INTERVAL BETWEEN Y(1) AND Y(NY))
C     ------------------------------------------------------------------

      ABINT=0.D0
      JX=1
      JY=1
      KX=1
      KY=1
      XA=X(JX)
      YA=Y(JY)
      XB=X(NX)
      YB=Y(NY)
      XX=XA
      IF(XB.GT.XA) GO TO 100
      XA=XB
      XB=XX
      JX=NX
      KX=-1
  100 CONTINUE
      XX=YA
      IF(YB.GT.YA) GO TO 120
      YA=YB
      YB=XX
      JY=NY
      KY=-1
  120 CONTINUE
      XMIN=DMAX1(XA,YA)
      XMAX=DMIN1(XB,YB)
      IF(XMIN.GE.BLIM) GO TO 190
      IF(XMAX.LE.ALIM) GO TO 190
      IF(XMIN.LT.ALIM) XMIN=ALIM
      IF(XMAX.GT.BLIM) XMAX=BLIM
  130 CONTINUE
      JX=JX+KX
      XJ=X(JX)
      IF(XJ.LE.XMIN) GO TO 130
      IX=JX-KX
      XI=X(IX)
      IF((XJ-XI).LT.DELTA) GO TO 130
      FI=F(IX)
      FJ=F(JX)
      BF=(FJ-FI)/(XJ-XI)
      AF=FJ-BF*XJ
  140 CONTINUE
      JY=JY+KY
      YJ=Y(JY)
      IF(YJ.LE.XMIN) GO TO 140
      IY=JY-KY
      YI=Y(IY)
      IF((YJ-YI).LT.DELTA) GO TO 140
      GI=G(IY)
      GJ=G(JY)
      BG=(GJ-GI)/(YJ-YI)
      AG=GJ-BG*YJ
      X2=XMIN
  160 CONTINUE
      X1=X2
      X2=DMIN1(XJ,YJ)
      IF(X2.GT.XMAX) X2=XMAX
      DINT=(AF*AG)*(X2-X1)
     *    +(AF*BG+BF*AG)*(X2**2-X1**2)/2.D0
     *    +(BF*BG)*(X2**3-X1**3)/3.D0
      ABINT=ABINT+DINT
      IF(DABS(X2-XMAX).LT.DELTA) GO TO 190
      IF((XJ-X2).GT.DELTA) GO TO 180
  170 CONTINUE
      XI=XJ
      FI=FJ
      IX=JX
      JX=JX+KX
      XJ=X(JX)
      FJ=F(JX)
      IF(DABS(XJ-XI).LT.DELTA) GO TO 170
      BF=(FJ-FI)/(XJ-XI)
      AF=FJ-BF*XJ
  180 CONTINUE
      IF(YJ.GT.X2) GO TO 160
      YI=YJ
      GI=GJ
      IY=JY
      JY=JY+KY
      YJ=Y(JY)
      GJ=G(JY)
      IF(DABS(YJ-YI).LT.DELTA) GO TO 180
      BG=(GJ-GI)/(YJ-YI)
      AG=GJ-BG*YJ
      GO TO 160
  190 CONTINUE
      RETURN
      END SUBROUTINE FXGINT

      BLOCK DATA RADPAR
C
      INCLUDE 'BR00B.COM'
C
C
C     Lean Spectral irradiance on following wavelength (nm) grid centers
C         with the following wavelength bands (nm) centered on above wls
cg    DATA WSLEAN/
cg   A  117.50, 122.50, 127.50, 135.00, 145.00, 155.00, 165.00, 175.00,
cg   B  185.00, 195.00, 205.00, 215.00, 222.50, 227.50, 232.50, 237.50,
cg   C  242.50, 247.50, 252.50, 257.50, 262.50, 267.50, 272.50, 277.50,
cg   D  282.50, 287.50, 292.50, 297.50, 302.50, 307.50, 312.50, 317.50,
cg   E  322.50, 327.50, 332.50, 337.50, 342.50, 347.50, 352.50, 357.50,
cg   F  362.50, 367.50, 372.50, 377.50, 382.50, 387.50, 392.50, 397.50,
cg   G  402.50, 407.50, 412.50, 417.50, 422.50, 427.50, 432.50, 437.50,
cg   H  442.50, 447.50, 452.50, 457.50, 462.50, 467.50, 472.50, 477.50,
cg   I  482.50, 487.50, 492.50, 497.50, 502.50, 507.50, 512.50, 517.50,
cg   J  522.50, 527.50, 532.50, 537.50, 542.50, 547.50, 552.50, 557.50,
cg   K  562.50, 567.50, 572.50, 577.50, 582.50, 587.50, 592.50, 597.50,
cg   L  602.50, 607.50, 615.00, 625.00, 635.00, 645.00, 655.00, 665.00,
cg   M  675.00, 685.00, 695.00, 705.00, 715.00, 725.00, 735.00, 745.00,
cg   N  755.00, 765.00, 775.00, 785.00, 795.00, 805.00, 815.00, 825.00,
cg   O  835.00, 845.00, 855.00, 865.00, 875.00, 885.00, 895.00, 905.00,
cg   P  915.00, 925.00, 935.00, 945.00, 955.00, 965.00, 975.00, 985.00,
cg   Q  995.00,1025.00,1075.00,1125.00,1175.00,1225.00,1275.00,1325.00,
cg   R 1375.00,1425.00,1475.00,1525.00,1575.00,1625.00,1675.00,1725.00,
cg   S 1775.00,1825.00,1875.00,1925.00,1975.00,2050.00,2150.00,2250.00,
C
ceq   DATA WLEAN2/
cg   A 2350.00,2450.00,2550.00,2650.00,2750.00,2850.00,2950.00,3050.00,
cg   B 3150.00,3250.00,3350.00,3450.00,3550.00,3650.00,3750.00,3850.00,
cg   C 3950.00,4050.00,4150.00,4250.00,4350.00,4450.00,4550.00,4650.00,
cg   D 4750.00,4850.00,4950.00,5500.00,6500.00,7500.00,8500.00,9500.00,
cg   E 10500.0,11500.0,12500.0,13500.0,14500.0,57500.0/
C
cg    DATA DSLEAN/
cg   A   5.0,  5.0,  5.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,
cg   B  10.0, 10.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,
cg   C   5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,
cg   D   5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,
cg   E   5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,
cg   F   5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,
cg   G   5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,
cg   H   5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,
cg   I   5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,
cg   J  10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,
cg   K  10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,
cg   L  10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,
cg   M  10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 50.0,
cg   N  50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
cg   O  50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 100.,
cg   P  100., 100., 100., 100., 100., 100., 100., 100., 100., 100.,
cg   Q  100., 100., 100., 100., 100., 100., 100., 100., 100., 100.,
cg   R  100., 100., 100., 100., 100., 100., 100., 100., 100.,1000.,
cg   S 1000.,1000.,1000.,1000.,1000.,1000.,1000.,1000.,1000.,8.5E+04/
C
C      Lean Solar UV Flux  1950     1      1366.291100      1366.448785
cg    DATA FRLEAN/
cg   A 1.463600E-03,1.526467E+00,2.927199E-02,7.005412E-02,6.927072E-02,
cg   B 1.668700E-01,3.785660E-01,1.120390E+00,2.706021E+00,5.377562E+00,
cg   C 1.171729E+01,3.687586E+01,5.222137E+01,4.619777E+01,4.697743E+01,
cg   D 4.752818E+01,5.774523E+01,5.098351E+01,5.023283E+01,1.083320E+02,
cg   E 1.425212E+02,2.479949E+02,2.079845E+02,1.800741E+02,2.359987E+02,
cg   F 3.226166E+02,5.530246E+02,5.127549E+02,5.201678E+02,5.984619E+02,
cg   G 7.012117E+02,7.082256E+02,7.617171E+02,9.915579E+02,9.935029E+02,
cg   H 9.177540E+02,9.600278E+02,9.351981E+02,1.050501E+03,9.218711E+02,
cg   I 9.948324E+02,1.204492E+03,1.063533E+03,1.195645E+03,9.637095E+02,
cg   J 1.043356E+03,1.067941E+03,1.255184E+03,1.806284E+03,1.705940E+03,
cg   K 1.699230E+03,1.758670E+03,1.758020E+03,1.650341E+03,1.553274E+03,
cg   L 1.767777E+03,1.929360E+03,1.987496E+03,2.075871E+03,2.094106E+03,
cg   M 2.095205E+03,2.030999E+03,2.060222E+03,2.104337E+03,2.075792E+03,
cg   N 1.880081E+03,1.987674E+03,1.980032E+03,1.911138E+03,1.980752E+03,
cg   O 1.892300E+03,1.744249E+03,1.859837E+03,1.829749E+03,1.893766E+03,
cg   P 1.880606E+03,1.852021E+03,1.880300E+03,1.880174E+03,1.814568E+03,
cg   Q 1.803343E+03,1.795452E+03,1.813090E+03,1.799145E+03,1.821609E+03,
cg   R 1.749601E+03,1.765703E+03,1.764881E+03,1.739113E+03,1.725504E+03,
cg   S 1.668046E+03,1.662454E+03,1.637083E+03,1.593666E+03,1.520814E+03,
ceq   DATA F2LEAN/
cg   A 1.538301E+03,1.500148E+03,1.475486E+03,1.460724E+03,1.429136E+03,
cg   B 1.376901E+03,1.345483E+03,1.314613E+03,1.283486E+03,1.271272E+03,
cg   C 1.243894E+03,1.198048E+03,1.161738E+03,1.133216E+03,1.116448E+03,
cg   D 1.101408E+03,1.070334E+03,1.040272E+03,1.009878E+03,9.452869E+02,
cg   E 9.636011E+02,9.627610E+02,9.396095E+02,9.293125E+02,9.004694E+02,
cg   F 8.922053E+02,8.606506E+02,8.539306E+02,8.280872E+02,8.124155E+02,
cg   G 7.956104E+02,7.779253E+02,7.611612E+02,7.472715E+02,6.951995E+02,
cg   H 6.224917E+02,5.687959E+02,5.203566E+02,4.779112E+02,4.352621E+02,
cg   I 3.917702E+02,3.566886E+02,3.290026E+02,3.042735E+02,2.808124E+02,
cg   J 2.586596E+02,2.364040E+02,2.154338E+02,1.943942E+02,1.761917E+02,
cg   K 1.594696E+02,1.447463E+02,1.318530E+02,1.206954E+02,1.048143E+02,
cg   L 8.809857E+01,7.489415E+01,6.361460E+01,5.467346E+01,4.723338E+01,
cg   M 4.113513E+01,3.584235E+01,3.141150E+01,2.763573E+01,2.431808E+01,
cg   N 2.154878E+01,1.912381E+01,1.700619E+01,1.523811E+01,1.367492E+01,
cg   O 1.226657E+01,1.098629E+01,1.000669E+01,9.063802E+00,8.205545E+00,
cg   P 7.499935E+00,6.828093E+00,6.078227E+00,5.463104E+00,4.993065E+00,
cg   Q 4.607049E+00,4.285606E+00,3.973604E+00,3.678151E+00,2.540683E+00,
cg   R 1.341226E+00,7.647217E-01,4.665436E-01,3.005465E-01,2.023893E-01,
cg   S 1.411388E-01,1.014034E-01,7.471365E-02,5.625803E-02,2.896050E-03/
C
      DATA WTHEK/.115,.120,.125,.130,.140,.150,.160,.170,.180,.190,.200,
     1 .210,.220,.225,.230,.235,.240,.245,.250,.255,.260,.265,.270,.275,
     2      .280,.285,.290,.295,.300,.305,.310,.315,.320,.325,.330,.335,
     3           .340,.345,.350,.355,.360,.365,.370,.375,.380,.385,.390,
     4           .395,.400,.405,.410,.415,.420,.425,.430,.435,.440,.445,
     5           .450,.455,.460,.465,.470,.475,.480,.485,.490,.495,.500,
     6           .505,.510,.515,.520,.525,.530,.535,.540,.545,.550,.555,
     7           .560,.565,.570,.575,.580,.585,.590,.595,.600,.605,.610,
     8           .620,.630,.640,.650,.660,.670,.680,.690,.700,.710,.720,
     9           .730,.740,.750,.760,.770,.780,.790,.800,.810,.820,.830,
     A .840,.850,.860,.870,.880,.890,.900,.910,.920,.930,.940,.950,.960,
     B 0.97,0.98,0.99,1.00,1.05,1.10,1.15,1.20,1.25,1.30,1.35,1.40,1.45,
     C 1.50,1.55,1.60,1.65,1.70,1.75,1.80,1.85,1.90,1.95,2.00,2.10,2.20,
     D 2.30,2.40,2.50,2.60,2.70,2.80,2.90,3.00,3.10,3.20,3.30,3.40,3.50,
     E 3.60,3.70,3.80,3.90,4.00,4.10,4.20,4.30,4.40,4.50,4.60,4.70,4.80,
     F      4.9, 5.0, 6.0, 7.0, 8.0, 9.0,10.0,11.0,12.0,13.0,14.0,15.00/
C
      DATA FTHEK/.007,.900,.007,.007,.030,.070,.230,.630,1.25,2.71,10.7,
     1 22.9,57.5,64.9,66.7,59.3,63.0,72.3,70.4,104.,130.,185.,232.,204.,
     2    222.,315.,482.,584.,514.,603.,689.,764.,830.,975.,1059.,1081.,
     31074.,1069.,1093.,1083.,1068.,1132.,1181.,1157.,1120.,1098.,1098.,
     41189.,1429.,1644.,1751.,1774.,1747.,1693.,1639.,1663.,1810.,1922.,
     52006.,2057.,2066.,2048.,2033.,2044.,2074.,1976.,1950.,1960.,1942.,
     61920.,1882.,1833.,1833.,1852.,1842.,1818.,1783.,1754.,1725.,1720.,
     71695.,1705.,1712.,1719.,1715.,1712.,1700.,1682.,1666.,1647.,1635.,
     81602.,1570.,1544.,1511.,1486.,1456.,1427.,1402.,1389.,1344.,1314.,
     91290.,1260.,1235.,1211.,1185.,1159.,1134.,1109.,1085.,1060.,1036.,
     A1013.,990.,968.,947.,926.,908.,891.,880.,869.,858.,847.,837.,820.,
     B 803.,785.,767.,748.,668.,593.,535.,485.,438.,397.,358.,337.,312.,
     C 288.,267.,245.,223.,202.,180.,159.,142.,126.,114.,103., 90., 79.,
     D 69.0,62.0,55.0,48.0,43.0,39.0,35.0,31.0,26.0,22.6,19.2,16.6,14.6,
     E 13.5,12.3,11.1,10.3, 9.5,8.70,7.80,7.10,6.50,5.92,5.35,4.86,4.47,
     F   4.11,3.79,1.82,0.99,.585,.367,.241,.165,.117,.0851,.0634,.0481/
C
      DATA XWAVO3/.2002,.2012,.2022,.2032,.2042,.2052,.2062,.2072,.2082,
     A.2092,.2102,.2112,.2122,.2132,.2142,.2152,.2162,.2172,.2182,.2192,
     B.2202,.2212,.2222,.2232,.2242,.2252,.2262,.2272,.2282,.2292,.2302,
     C.2312,.2322,.2332,.2342,.2352,.2362,.2372,.2382,.2392,.2400,.2402,
     D.2412,.2422,.2432,.2438,.2444,.2452,.2458,.2463,.2472,.2478,.2482,
     E.2490,.2492,.2500,.2508,.2519,.2527,.2539,.2543,.2553,.2562,.2566,
     F.2571,.2575,.2579,.2587,.2597,.2604,.2617,.2624,.2635,.2643,.2650,
     G.2654,.2662,.2669,.2675,.2682,.2692,.2695,.2702,.2712,.2718,.2722,
     H.2732,.2742,.2746,.2752,.2762,.2772,.2782,.2792,.2802,.2812,.2822,
     I.2830,.2842,.2852,.2862,.2872,.2882,.2892,.2902,.2912,.2922,.2932,
     J.2942,.2952,.2962,.2972,.2982,.2992,.2998,
     &            .3004,.3016,.3021,.3029,.3036,.3037,.3051,.3053,.3059,
     A.3061,.3066,.3075,.3077,.3083,.3085,.3092,.3098,.3100,.3104,.3106,
     B.3109,.3112,.3130,.3135,.3146,.3148,.3151,.3154,.3167,.3170,.3173,
     C.3176,.3190,.3194,.3199,.3200,.3209,.3210,.3216,.3220,.3223,.3226,
     D.3239,.3242,.3245,.3248,.3253,.3255,.3269,.3272,.3275,.3279,.3292,
     E.3295,.3299,.3303,.3309,.3312,.3328,.3332,.3334,.3338,.3357,.3365,
     F.3369,.3372,.3391,.3395,.3398,.3401,.3417,.3421,.3426,.3430,.3437,
     G.3439,.3451,.3455,.3460,.3463,.3466,.3472,.3481,.3485,.3489,.3493,
     H.3499,.3501,.3506,.3514,.3521,.3523,.3546,.3550,.3554,.3556,.3561,
     I.3567,.3572,.3573,.3588,.3594,.3599,.3600,.3604,.3606,.3639,.3647,
     J.3650,.3654,.3660/
      DATA FUVKO3/  8.3,  8.3,  8.1,  8.3,  8.6,  9.0,  9.7, 10.8, 11.7,
     A 13.0, 14.3, 16.0, 18.0, 20.6, 23.0, 26.1, 29.3, 32.6, 36.9, 40.8,
     B 46.9, 51.4, 56.7, 63.4, 69.1, 76.6, 84.0, 91.4, 99.9,110.0,118.0,
     C126.0,136.0,145.0,154.0,164.0,175.0,186.0,192.0,201.0,210.0,212.0,
     D221.0,230.0,239.0,248.0,250.0,259.0,264.0,264.0,273.0,277.0,275.0,
     E283.0,283.0,290.0,283.0,297.0,290.0,300.0,290.0,302.0,295.0,283.0,
     F293.0,290.0,286.0,297.0,281.0,280.0,271.0,275.0,254.0,264.0,250.0,
     G248.0,242.0,228.0,230.0,216.0,213.0,211.0,199.0,188.0,188.0,178.0,
     H169.0,153.0,155.0,148.0,136.0,127.0,117.0,108.0, 97.0, 88.7, 81.3,
     I 78.7, 67.9, 61.4, 54.3, 49.6, 43.1, 38.9, 34.6, 30.2, 27.5, 23.9,
     J 21.0, 18.6, 16.2, 14.2, 12.3, 10.7,  9.5,
     &            8.880,7.520,6.960,6.160,5.810,5.910,4.310,4.430,4.130,
     A4.310,4.020,3.330,3.390,3.060,3.100,2.830,2.400,2.490,2.330,2.320,
     B2.120,2.200,1.436,1.595,1.074,1.138,1.068,1.262,0.818,0.948,0.860,
     C1.001,0.543,0.763,0.665,0.781,0.382,0.406,0.373,0.608,0.484,0.601,
     D0.209,0.276,0.259,0.470,0.319,0.354,0.131,0.223,0.185,0.339,0.080,
     E0.093,0.079,0.184,0.139,0.214,0.053,0.074,0.068,0.152,0.038,0.070,
     F.0540000,.1030000,.0240000,.0382500,.0292500,.0550000,.0135000,
     G.0155250,.0127500,.0188250,.0167250,.0262500,.0115500,.0140250,
     H.0099750,.0115500,.0081000,.0104250,.0050100,.0057000,.0046650,
     I.0073425,.0051825,.0055275,.0040575,.0077700,.0048900,.0054600,
     J.0015375,.0017775,.0013275,.0014100,.0011550,.0023325,.0018825,
     K.0019650,.0009600,.0013650,.0011925,.0013200,.0008925,.0009825,
     L.0001350,.0006300,.0004500,.0006225,0.0/
C
C
C-----------------------------------------------------------------------
C          NO2 Trace Gas Vertical Distribution and Concentration Profile
C-----------------------------------------------------------------------
C
      DATA CMANO2/
     1  8.66E-06,5.15E-06,2.85E-06,1.50E-06,9.89E-07,6.91E-07,7.17E-07,
     2  8.96E-07,3.67E-06,4.85E-06,5.82E-06,6.72E-06,7.77E-06,8.63E-06,
     3  8.77E-06,8.14E-06,6.91E-06,5.45E-06,4.00E-06,2.67E-06,1.60E-06,
     4  8.36E-07,3.81E-07,1.58E-07,6.35E-08,2.57E-08,1.03E-08,4.18E-09,
     5  1.66E-09,6.57E-10,2.58E-10,1.02E-10,4.11E-11,1.71E-11,7.73E-12,
     6  9.07E-12,4.63E-12,2.66E-12,1.73E-12,1.28E-12,1.02E-12,1.00E-30/
C
      DATA PVT/8*0.125,3*0.0/
C
      DATA SRBXAL/30*0.0/
C
      DATA DLAT46/
     A    -90.,-86.,-82.,-78.,-74.,-70.,-66.,-62.,-58.,-54.,-50.,-46.,
     B    -42.,-38.,-34.,-30.,-26.,-22.,-18.,-14.,-10., -6., -2.,  2.,
     C      6., 10., 14., 18., 22., 26., 30., 34., 38., 42., 46., 50.,
     D     54., 58., 62., 66., 70., 74., 78., 82., 86., 90./
C
      DATA DLON72/
     A      0.,  5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55.,
     B     60., 65., 70., 75., 80., 85., 90., 95.,100.,105.,110.,115.,
     C    120.,125.,130.,135.,140.,145.,150.,155.,160.,165.,170.,175.,
     D    180.,185.,190.,195.,200.,205.,210.,215.,220.,225.,230.,235.,
     E    240.,245.,250.,255.,260.,265.,270.,275.,280.,285.,290.,295.,
     F    300.,305.,310.,315.,320.,325.,330.,335.,340.,345.,350.,355./
C
!obso DATA PLB/984.,934.,854.,720.,550.,390.,285.,210.,150.,100.,60.,
!obso+          30.,10.,5.,2.,0.001,41*0.0/
C
!obso DATA HLB/
!obso+          0.25, 0.70, 1.47, 2.89, 5.06, 7.65, 9.90,11.97,14.12,
!obso+         16.67,19.95,24.47,31.97,36.98,44.00,99.99,41*99.999/
C
C     DATA PLB/
C    +     1013.2500, 961.7485, 879.3460, 741.3219, 566.2166, 401.4117,
C    +      262.3575, 154.2043,  71.8018,  10.0000,   5.0000,   2.0000,
C    +        1.E-05,  27*0.00/
C
C     DATA HLB/
C    +         1.E-10,1.0,2.0,3.0,6.0,11.0,20.0,32.0,47.0,51.0,71.0
C    +        ,84.852,99.99,27*99.999/
C
!!    DATA TAUSL/33*0.0/, FTAUSL/33*0.0/
C
C-----------------------------------------------------------------------
C     TRACE GAS REFERENCE AMOUNTS & DISTRIBUTIONS ARE DEFINED IN  SETGAS
C-----------------------------------------------------------------------
C
C
C                  H2O  CO2   O3   O2  NO2  N2O  CH4  F11  F12
C                   1    2    3    4    5    6    7    8    9
      DATA FULGAS/ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
C
C                  N2C   CFC11+  CFC12+  FO3LON
     +           , 1.0,   1.0,    1.0,    1.0/
C
C
C               GLOBAL OCEAN LAND DESERT  HAZE   TR1  TR2  TR3  TR4
C                   1    2    3     4      5      6    7    8    9
C
C                                                TR5  TR6  TR7  TR8
C                                                 10   11   12   13
C     DATA FGOLDH/13*0.0/
C
      DATA FGOLDH/ 1.0, .68, .32, 1.d-20,1.d-20, 0.0, 0.0, 0.0, 0.0
     +           , 0.0, 0.0, 0.0, 0.0/
C
      DATA    LASTVC/-123456/
C
      DATA   FSXAER/5*1.D0/,  FTXAER/5*1.D0/
      DATA   PIVMAX/1.0000/
C
      DATA   TLGRAD/ 1.0/,    RO3COL/ 1.0/, ECLTRA/1.00/
      DATA   PTLISO/2.5E+00/, ZCLSRA/1.00/
C
      DATA   FCLDTR/1.0/,     FCLDSR/1.0/,  NTRACE/0/
      DATA                    ITR/8*0/,     FRAYLE/1.0/,  EPSCON/0.0/
C
      DATA   KSOLAR/1/,       KTREND/1/,    MADBAK/0/
      DATA   KEEPRH/0/,       KEEP10/0/,    NO3COL/0/,    KZSNOW/ 1/
      DATA   KCNORM/0/,       KEEPAL/0/,    MRELAY/0/
      DATA   KCLDEP/4/,       KO3LON/0/,    MOZONE/6/
      DATA   MADGHG/1/,       MADSUR/1/,    ICE012/1/
      DATA   KWVCON/1/,       KSNORM/0/,    NORMS0/1/
C
      DATA                    JLAT/ 7/,   S00WM2/1366.2911/, S0/1366.0/
      DATA                    ILON/18/,   COSZ/0.5000/
C
cg    DATA   POCEAN/0.700/,   TGO/288.15/,   AGESN/3*1./,    WMAG/2.00/
cg    DATA   PEARTH/0.000/,   TGE/288.15/,   SNOWE/0.30/,  WEARTH/0.00/
cg    DATA    POICE/0.200/,  TGOI/288.15/,  SNOWOI/0.10/,
cg    DATA    PLICE/0.100/,  TGLI/288.15/,  SNOWLI/0.20/
cg    DATA                    TSL/288.15/
C
      DATA KVRAER/1/,  MEANAC/0/,  MEANDD/0/,  MEANVA/0/,  NCARO3/0/
      DATA KUVFAC/0/,  KWTRAB/0/   KCLDEM/1/,  KANORM/0/
C
      DATA KPFCO2/0/,  KPFOZO/0/,  KVEGA6/0/,  KORDER/0/
      DATA KUFH2O/1/,  KUFCO2/1/,              KCSELF/1/,   KCFORN/1/
C
      DATA KYEARS/    0/, KJDAYS/  0/
      DATA KYEARG/    0/, KJDAYG/  0/
      DATA KYEARO/    0/, KJDAYO/  0/
      DATA KYEARA/    0/, KJDAYA/  0/
      DATA KYEARD/    0/, KJDAYD/  0/
      DATA KYEARV/    0/, KJDAYV/  0/
      DATA KYEARE/    0/, KJDAYE/  0/
      DATA KYEARR/    0/, KJDAYR/  0/
C
      DATA                              REFF0/0.30d0/,  VEFF0/0.35d0/
C
      DATA                              O3WJT0/0.001/
C
      DATA              X0YBCI/0.001/,  X0YOCI/0.001/,  X0YSUI/0.001/
C
      DATA UVWAVL/0.0010,0.0020,0.0030/
      DATA UVFACT/1.0000,1.0000,1.0000/
C
      DATA PTOPTR/150.0/
C
C
C-----------------------------------------------------------------------
C      COMPOSITION & VERTICAL DISTRIBUTION FOR 5 SPECIFIED AEROSOL TYPES
C-----------------------------------------------------------------------
C TYPE
C    1   STRATOSPHERIC GLOBAL AEROSOL  A,B,C ARE GLOBAL AVERAGE VALUES
C    2    TROPOSPHERIC  OCEAN AEROSOL  A,B,C ARE GLOBAL AVERAGE VALUES
C    3    TROPOSPHERIC   LAND AEROSOL  A,B,C ARE GLOBAL AVERAGE VALUES
C    4    TROPOSPHERIC DESERT AEROSOL  A,B,C ARE  LOCAL AVERAGE VALUES
C    5    TROPOSPHERIC   HAZE AEROSOL  A,B,C ARE  LOCAL AVERAGE VALUES
C
C        1     2     3     4     5     6     7     8     9    10    11
C      ACID1 SSALT SLFT1 SLFT2 BSLT1 BSLT2 DUST1 DUST2 DUST3 CARB1 CARB2
      DATA AGOLDH/
     1 .005,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,
     2   .0, .020, .010, .010, .005,   .0, .010,   .0,   .0, .005,   .0,
     3   .0,   .0,   .0, .020, .005,   .0, .010, .010,   .0,   .0, .015,
     4   .0,   .0,   .0,   .0,   .0,   .0,   .0, .020, .010,   .0,   .0,
     5   .0,   .0,   .0, .010,   .0,   .0,   .0,   .0,   .0,   .0, .005/
      DATA BGOLDH/
     1 20.0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,
     2   .0, 1.00, 4.00, 1.00, 4.00, 1.00, 4.00,   .0,   .0, 1.00,   .0,
     3   .0,   .0,   .0, 0.00, 2.00,   .0, 4.00, 2.00,   .0,   .0, 0.00,
     4   .0,   .0,   .0,   .0,   .0,   .0,   .0, 2.00, 0.00,   .0,   .0,
     5   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0, 0.00/
      DATA CGOLDH/
     1 3.00,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,
     2   .0, 1.00, 3.00, 2.00, 3.00, 1.00, 2.00,   .0,   .0, 1.00,   .0,
     3   .0,   .0,   .0, 1.00, 3.00,   .0, 1.00, 1.00,   .0,   .0, 1.00,
     4   .0,   .0,   .0,   .0,   .0,   .0,   .0, 1.00, 1.00,   .0,   .0,
     5   .0,   .0,   .0, 1.00,   .0,   .0,   .0,   .0,   .0,   .0, 1.00/
C
C
C           1          2          3          4          5          6
C         ACID1      SSALT      SLFT1      SLFT2      BSLT1      BSLT2
      DATA PI0VIS/
     1   1.00000,   1.00000,   1.00000,   1.00000,   0.98929,   0.95609,
C
C           7          8          9         10         11
C         DUST1      DUST2      DUST3      CARB1      CARB2
     2   0.91995,   0.78495,   0.63594,   0.31482,   0.47513/
C
C
C          ------------------------------------------------------------
C          Multiplicative Scaling Factors for Cloud Asymmetry Parameter
C          ------------------------------------------------------------
C
C                  Water    Ice    MieIce
      DATA SRCGSF/ 1.000,  1.000,  1.000/
C
      DATA TAUWC0/0.0010/, TAUIC0/0.0010/
C
C                       TROPOSPHERIC AEROSOL PARAMETERS
C                  BCI  OCI  SUI  SEA  SUN    ANT  OCN  OCB  BCB  SSB
      DATA REAERO/ 0.1, 0.3, 0.3, 2.0, 0.3,   1.0, 0.3, 0.3, 0.2, 0.5/

      DATA VEAERO/ 0.2, 0.2, 0.2, 0.2, 0.2,   0.2, 0.2, 0.2, 0.2, 0.2/

      DATA ROAERO/ 2.0, 2.0, 2.0, 2.0, 2.0,   2.0, 2.0, 2.0, 2.0, 2.0/

C     DATA PI0MAX/ 1.0, 1.0, 1.0, 1.0, 1.0,   1.0, 1.0, 1.0, 1.0, 1.0/
C       (OCI,OCN,OCB PI0MAXs are based on Novakov(1998) Ni=0.005)
      DATA PI0MAX/ 1.0, .96, 1.0, 1.0, 1.0,   1.0, .98, .93, 1.0, 1.0/

      DATA FSAERO/ 0.5, 0.5, 1.0, 4.0, 1.0,   0.0, 1.0, 1.0, 1.0, 1.0/

      DATA FTAERO/ 0.5, 0.5, 1.0, 4.0, 1.0,   0.0, 1.0, 1.0, 1.0, 1.0/

      DATA FRSULF/ 0.0, .33, 0.0, 0.0, 0.0,   0.0, .33, .33, 0.0, 0.0/

      DATA SSBTAU/0.005/

C                          MINERAL DUST PARAMETERS
C                         CLAY                  SILT
      DATA REDUST/ 0.1, 0.2, 0.4, 0.8,   1.0, 2.0, 4.0, 8.0/
C
      DATA VEDUST/ 0.2, 0.2, 0.2, 0.2,   0.2, 0.2, 0.2, 0.2/
C
      DATA RODUST/ 2.5, 2.5, 2.5, 2.5,   2.6, 2.6, 2.6, 2.6/
C
      DATA FSDUST/ 1.0, 1.0, 1.0, 1.0,   1.0, 1.0, 1.0, 1.0/
C
      DATA FTDUST/ 1.0, 1.0, 1.0, 1.0,   1.0, 1.0, 1.0, 1.0/
C
C
C              TROPOSPHERIC AEROSOL GLOBAL MEAN VERTICAL DISTRIBUTION
C              ------------------------------------------------------
C     Layr  1    2    3    4    5    6    7    8    9    10   11   12
C     Ptop 934  854  720  550  390  285  255  150  100   60   30   10
C     Pbot 984  934  854  720  550  390  285  255  150  100   60   30
      DATA VDGAER/
C       L=   1    2    3    4    5    6    7    8    9   10   11   12
     A     .10, .10, .10, .10, .10, .10, .10, .10, .10, .10, .00, .00,
     B     .10, .10, .10, .10, .10, .10, .10, .10, .10, .10, .00, .00,
     C     .10, .10, .10, .10, .10, .10, .10, .10, .10, .10, .00, .00,
     D     .30, .20, .20, .15, .10, .05, .00, .00, .00, .00, .00, .00,
     E     .20, .20, .20, .20, .10, .10, .00, .00, .00, .00, .00, .00,
     F     .20, .20, .20, .20, .10, .10, .00, .00, .00, .00, .00, .00,
     G     .20, .20, .20, .20, .10, .10, .00, .00, .00, .00, .00, .00,
     H     .20, .20, .20, .20, .10, .10, .00, .00, .00, .00, .00, .00,
     I     .20, .20, .20, .20, .10, .10, .00, .00, .00, .00, .00, .00,
     J     .00, .00, .00, .00, .00, .00, .25, .20, .20, .20, .10, .05/
C
C            Aerosol Vertical Profile (Optical Depth) Scaling Factors
C            --------------------------------------------------------
C     Layr  1    2    3    4    5    6    7    8    9    10   11   12
C     Ptop 934  854  720  550  390  285  255  150  100   60   30   10
C     Pbot 984  934  854  720  550  390  285  255  150  100   60   30
      DATA VDFBCI/
     A     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
      DATA VDFOCI/
     A     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
      DATA VDFSUI/
     A     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
      DATA VDFDST/
     A     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
C
C
C-----------------------------------------------------------------------
C     GHG 1980 Reference Concentrations and Vertical Profile Definitions
C-----------------------------------------------------------------------
C
C     NUMBER   1    2    3    4  5    6    7    8     9   10   11  12
C             H2O  CO2  O3   O2 NO2  N2O  CH4 CFC11 CFC12 N2 CF-Y  CF-Z
C     FULGAS/ 1.0, 1.0,1.0, 1.0,1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0 /
      DATA Z0/0.0, 0.0,0.0, 0.0,0.0, 16., 16., 16., 16., 0.0, 16., 16.0/
      DATA ZH/8.0, 8.0,8.0, 8.0,8.0, 30., 50., 30., 30., 0.0, 30., 30.0/
C
C     GAS  NUMBER  1      2     3      4      5       6       7
C                 H2O    CO2    O3    O2     NO2     N2O     CH4
      DATA PPMV80/0.0, 337.90, 0.0, 210000., 0.0,  0.3012, 1.5470
C
C     GAS  NUMBER   8           9        10        11         12
C                 CCL3F1      CCL2F2     N2       CFC-Y      CFC-Z
     +         ,0.1666E-03, 0.3003E-03, 0.0,   0.0978D-03, 0.0010D-10/
C
C     KTREND=0  Direct Set Default Values
C     -----------------------------------
C
C     GAS  NUMBER  1      2     3      4      5       6       7
C                 H2O    CO2    O3    O2     NO2     N2O     CH4
      DATA PPMVK0/0.0, 337.90, 0.0, 210000., 0.0,  0.3012, 1.5470
C
C     GAS  NUMBER   8           9        10        11         12
C                 CCL3F1      CCL2F2     N2       CFC-Y      CFC-Z
     +         ,0.1666E-03, 0.3003E-03, 0.0,   0.0978D-03, 0.0010D-10/
C
C     Makiko's GHG Trend Compilation  GHG.1850-2050.Dec1999 in GTREND
C     GTREND:  1980.,  337.9,  .3012,  1.547,  .1666,  .3003,  .0978,
C     ---------------------------------------------------------------
C
C     ---------------------------------------------------------------
C     Minschwaner et al JGR (1998) CH4, N2O, CFC-12 Vertical profiles
C       IF(KGGVDF.GT.0) Then:
C     Gas decreases are linear with pressure, from unity at ground to
C     the fractional value PPMVDF(NGAS) at the top of the atmosphere.
C     Exponential decrease by EXP(-(Z-Z0)/H) is superimposed on this.
C       IF(KLATZ0.GT.0) Then: Z0 depends on latitude, KGGVDF not used
C     Pole-to-Pole latitudinal gradient (PPGRAD) is also superimposed
C     ---------------------------------------------------------------
C
      DATA KGGVDF/0/, KPGRAD/1/, KLATZ0/1/
C
C     GAS  NUMBER  1     2    3    4    5       6         7
C                 H2O   CO2   O3   O2  NO2     N2O       CH4
      DATA PPMVDF/1.0,  1.0, 1.0, 1.0, 1.0,  0.88888,  0.88888,
C
C     GAS  NUMBER   8          9      10       11        12
C                 CCL3F1    CCL2F2    N2     CFC-Y      CFC-Z
     +           0.88888,  0.88888,  1.0,   0.88888,   0.88888/
C
C     GAS  NUMBER  1     2    3    4    5        6         7
C                 H2O   CO2   O3   O2  NO2      N2O       CH4
      DATA PPGRAD/0.0,  0.0, 0.0, 0.0, 0.0,   0.0100,   0.0900,
C
C     GAS  NUMBER   8          9       10      11        12
C                 CCL3F1    CCL2F2     N2     CFC-Y     CFC-Z
     +            0.0600,   0.0600,   0.0,   0.0600,   0.0600/
C
      DATA DKS0/ .010, .030, .040, .040, .040, .002, .004, .013,
     +           .002, .003, .003, .072, .200, .480, .050, .011/

      DATA NKSLAM/14/, KSLAM/1,1,2,2,5,5,5,5,1,1,1,3,4,6,6,1/

      END BLOCK DATA RADPAR

c*********************************************************************
c  The following part computes the albedo. It has been taken out
c  from the radiation module with the intention to convert it later
c  into a separate module.
c*********************************************************************

      SUBROUTINE UPDSUR(JYEARR,JJDAYR)
!@sum UPDSUR updates variables for surface albedo once a day
!@auth A. Lacis/V. Oinas (modifications by I. Alienov/G. Schmidt)
      use SURF_ALBEDO, only : albvnh,albvnd,season,nv
      implicit none
!@var jyearr radiation year (not used for anything)
!@var jjdayr julian day (used for seasonality of veg albedo)
      integer, intent(in) :: jyearr, jjdayr
      integer K,KS1,KS2,KN1,KN2,L
      real*8 XJDAY,SEASN1,SEASN2,WT2,WT1
C
C                      Define Seasonal Albedo Dependence
C                      ---------------------------------
C
      XJDAY=JJDAYR
      SEASN1=-77.0D0
      DO K=1,4
        SEASN2=SEASON(K)
        IF(XJDAY.LE.SEASN2) GO TO 120
        SEASN1=SEASN2
      END DO
      K=1
      SEASN2=380.0D0
  120 CONTINUE
      WT2=(XJDAY-SEASN1)/(SEASN2-SEASN1)
      WT1=1.D0-WT2
      KS1=1+MOD(K,4)
      KS2=1+MOD(K+1,4)
      KN1=1+MOD(K+2,4)
      KN2=K

      DO K=1,NV
      DO L=1,6
C     -------------------
C     Southern Hemisphere
C     -------------------
        ALBVNH(K,L,1)=WT1*ALBVND(K,KS1,L)+WT2*ALBVND(K,KS2,L)
C     -------------------
C     Northern Hemisphere
C     -------------------
        ALBVNH(K,L,2)=WT1*ALBVND(K,KN1,L)+WT2*ALBVND(K,KN2,L)
      END DO
      END DO
      RETURN
      END SUBROUTINE UPDSUR

      SUBROUTINE GETSUR
!@sum GETSUR computes surface albedo for each grid box
!@auth A. Lacis/V. Oinas (modifications by I. Alienov/G. Schmidt)
      use SURF_ALBEDO
      use RE001, only:
C**** config data
     *     MLAT46,jnorth,KEEPAL,KVEGA6,snoage_fac_max,KZSNOW,
C**** input from radiation
     *     COSZ,PLANCK,ITNEXT,ITPFT0,
C**** input from driver
     *     AGESN,ILON,JLAT,POCEAN,POICE,PEARTH,PLICE,PLAKE,LKDEPTH,
     *     TGO,TGOI,TGE,TGLI,TSL,HIN,FLAGS,FMP,HSN,HMP,
     *     SNOWOI,SNOWE,SNOWLI,SNOW_FRAC,WEARTH,WMAG,PVT,
C**** output
     *     BXA,PRNB,PRNX,SRBALB,SRXALB,TRGALB,BGFEMD,BGFEMT,
     *     DTRUFG,FTRUFG
      implicit none
      integer K,J,I,L,JH,IWM,JWM,ITOC,ITEA,ITOI,ITLI
      real*8 ASNAGE,FSNAGE,BOCSUM,BEASUM,BOISUM,BLISUM,AVSCUM,ANSCUM,WMJ
     *     ,WMI,FRFOAM,AV,BV,WTOC,BOCM,BOCP,TRAPOC,BOCM1,BOCP1,BOC
     *     ,DSFRAC,VGFRAC,SEAVIS,SEANIR,VTFRAC,WTEA,BEAM,BEAP,TRAPEA
     *     ,BEAM1,BEAP1,BEA,WTOI,BOIM,BOIP,TRAPOI,BOIM1,BOIP1,BOI,AMEAN
     *     ,WTLI,BLIM,BLIP,BGF,TRAPLI,BLIM1,BLIP1,BLI,KKZSNO

C**** local variables for albedos for each surface type
      real*8 BOCVIS,BEAVIS,BOIVIS,BLIVIS,BOCNIR,BEANIR,BOINIR,BLINIR,
     +       XOCVIS,XEAVIS,XOIVIS,XLIVIS,XOCNIR,XEANIR,XOINIR,XLINIR,
     +       EXPSNE,EXPSNO,EXPSNL,BSNVIS,BSNNIR,XSNVIS,XSNNIR,
     +       BVSOIL,BNSOIL,BVVEGE,BNVEGE,XVSOIL,XNSOIL,XVVEGE,XNVEGE,
     +       BVSURF,BNSURF,XVSURF,XNSURF

C**** arrays needed if 6 band albedo is used
      real*8, dimension(6) :: BOCVN,BEAVN,BOIVN,BLIVN,BSNVN,BVNSUR,
     *                        XOCVN,XEAVN,XOIVN,XLIVN,XSNVN,XVNSUR

C**** Equivalence 2 band variables to 6 band array for easier passing
      EQUIVALENCE
     *     (BOCVN(1),BOCVIS),(BOCVN(2),BOCNIR),
     *     (BEAVN(1),BEAVIS),(BEAVN(2),BEANIR),
     *     (BOIVN(1),BOIVIS),(BOIVN(2),BOINIR),
     *     (BLIVN(1),BLIVIS),(BLIVN(2),BLINIR),
     *     (BSNVN(1),BSNVIS),(BSNVN(2),BSNNIR),
     *     (XOCVN(1),XOCVIS),(XOCVN(2),XOCNIR),
     *     (XEAVN(1),XEAVIS),(XEAVN(2),XEANIR),
     *     (XOIVN(1),XOIVIS),(XOIVN(2),XOINIR),
     *     (XLIVN(1),XLIVIS),(XLIVN(2),XLINIR),
     *     (XSNVN(1),XSNVIS),(XSNVN(2),XSNNIR),
     *     (BVNSUR(1),BVSURF),(BVNSUR(2),BNSURF),
     *     (XVNSUR(1),XVSURF),(XVNSUR(2),XNSURF)

C**** variables used for sea ice albedo calculation (for 4 bands)
      real*8, dimension(4) :: almp,alsd,alsf,ali,albtf,albtr
      real*8 :: patchy,snagfac

C     -----------------------------------------------------------------
C     Ocean Albedo Dependence on Zenith Angle and Wind Speed
C
      real*8 BVH2O, XVH2O, WMAG1, X
      BVH2O(WMAG1)=.0488D0+.0974D0/(5.679D0+WMAG1)+
     &     .0004D0/(.3333D0+WMAG1)
      XVH2O(WMAG1,X)=.021D0+X*X*(.0421D0+X*(.1283D0+X*(-.04D0+X*(3.117D0
     +              /(5.679D0+WMAG1)+X*.025D0/(.3333D0+WMAG1)))))
C     -----------------------------------------------------------------
C
C-----------------------------------------------------------------------
C     Select albedo computations and fixups using KVEGA6
C     KVEGA6=-5  2-band albedo, Antarc/Greenl alb=.8, we puddling:SI95
C     KVEGA6=-4  2-band albedo, Antarc/Greenl alb=.8, ws puddling:R00BF
C     KVEGA6=-3  2-band albedo, Antarc/Greenl alb=.8, wk puddling:R00BG
C     KVEGA6=-2  2-band albedo, Antarc/Greenl alb=.8, no puddling
C     KVEGA6=-1  2-band albedo - no 'fixups'
C     KVEGA6= 0  Schramm oi.alb, Antarc/Greenl alb=.8
C     KVEGA6= 1  6-band albedo - no 'fixups'
C     KVEGA6= 2  6-band albedo, Antarc/Greenl alb=.8, no puddling
C     KVEGA6= 3  6-band Schramm oi.alb, Antarc/Greenl alb=.8
C
C           Get Albedo, Thermal Flux, Flux Derivative for each Surf Type
C           ------------------------------------------------------------
C
      JH=1
      IF(JLAT.GT.JNORTH) JH=2
      KKZSNO=KZSNOW
      IF(COSZ.LT.0.001) KKZSNO=0
C
      EXPSNE=1.D0 ; EXPSNO=1.D0 ;  EXPSNL=1.D0
C
      DO K=1,NKBAND
        TRGALB(K)=0.D0
        BGFEMD(K)=0.D0
        BGFEMT(K)=0.D0
      END DO
C
      BOCSUM=0.D0
      BEASUM=0.D0
      BOISUM=0.D0
      BLISUM=0.D0
      DO K=1,4
        DTRUFG(K)=0.D0
      END DO
C
      AVSCUM=0.D0
      ANSCUM=0.D0
C
      BOCVN=0. ; BEAVN=0. ; BOIVN=0. ; BLIVN=0.
      XOCVN=0. ; XEAVN=0. ; XOIVN=0. ; XLIVN=0.
C
C                                             --------------------------
C                                             Ocean Albedo Specification
C                                             --------------------------
C
      IF(POCEAN.LT.1.D-04) GO TO 400
      X=0.5D0+(0.5D0-COSZ)*ZOCSRA
      BOCVIS=BVH2O(WMAG)  +AVSCAT+AVSCUM
      XOCVIS=XVH2O(WMAG,X)+AVSCAT+AVSCUM
      BOCNIR=BVH2O(WMAG)  +ANSCAT+ANSCUM
      XOCNIR=XVH2O(WMAG,X)+ANSCAT+ANSCUM
C
      IWM=WMAG
      IF(IWM.LT.1) IWM=1
      IF(IWM.GT.24) IWM=24
      JWM=IWM+1
      WMJ=WMAG-IWM
      WMI=1.D0-WMJ
      FRFOAM=WMI*SRFOAM(IWM)+WMJ*SRFOAM(JWM)
C
      BOCVIS=BOCVIS*(1.D0-FRFOAM)+FRFOAM*AVFOAM
      XOCVIS=XOCVIS*(1.D0-FRFOAM)+FRFOAM*AVFOAM
      BOCNIR=BOCNIR*(1.D0-FRFOAM)+FRFOAM*ANFOAM
      XOCNIR=XOCNIR*(1.D0-FRFOAM)+FRFOAM*ANFOAM
      IF(KVEGA6.GT.0) THEN      ! fill in higher bands
        DO L=3,6                ! 1/2 already equivalenced
          BOCVN(L)=BOCNIR
          XOCVN(L)=XOCNIR
        END DO
      ENDIF

C**** For lakes increase albedo if lakes are very shallow
C**** This is a fix to prevent lakes from overheating when they
C**** are not allowed to evaporate away. We assume that departing
C**** lakes would leave bare soil behind.
      if (PLAKE.gt.0 .and. lkdepth.lt.1.) then  ! < 1m
        DO L=1,6
          BOCVN(L)=(BOCVN(L)*(lkdepth-0.4d0)+
     *         (1.-lkdepth)*ALBVNH(1,L,1))/0.6d0
          XOCVN(L)=(XOCVN(L)*(lkdepth-0.4d0)+
     *         (1.-lkdepth)*ALBVNH(1,L,1))/0.6d0
        END DO
      end if
C
      X=1.D0/(1.D0+WMAG)
      AV=(-.0147087D0*X*X+.0292266D0*X-.0081079D0)*EOCTRA
      BV=(1.01673D0-0.0083652D0*WMAG)*EOCTRA
C
      ITOC=TGO
      WTOC=TGO-ITOC
      ITOC=ITOC-ITPFT0
      IF(ITOC.LT.124) ITOC=124
      BOCSUM=0.D0
      BOCM=0.D0
      BOCP=0.D0
C
      DO K=1,NKBAND
        TRAPOC=AV+BV*AOCEAN(K)
        BOCM1 =(PLANCK(ITOC-1)-(PLANCK(ITOC-1)-PLANCK(ITOC  ))*WTOC)
     +       *(1.D0-TRAPOC)
        BOCM  =BOCM+BOCM1
        BOCP1 =(PLANCK(ITOC+1)-(PLANCK(ITOC+1)-PLANCK(ITOC+2))*WTOC)
     +       *(1.D0-TRAPOC)
        BOCP  =BOCP+BOCP1
        BOC   =(PLANCK(ITOC  )-(PLANCK(ITOC  )-PLANCK(ITOC+1))*WTOC)
     +       *(1.D0-TRAPOC)
        BOCSUM=BOCSUM+BOC
        ITOC=ITOC+ITNEXT
C
        TRGALB(K)=TRGALB(K)+POCEAN*TRAPOC
        BGFEMD(K)=BGFEMD(K)+POCEAN*(BOCP1-BOCM1)
        BGFEMT(K)=BGFEMT(K)+POCEAN*BOC
      END DO
      DTRUFG(1)=0.5D0*(BOCP-BOCM)
C
  400 CONTINUE
      IF(PEARTH.LT.1.D-04) GO TO 500
C
C                                         ------------------------------
C                                         Land Snow Albedo Specification
C                                         ------------------------------
      ASNAGE=ALBDIF(1,JH)*EXP(-AGEXPF(1,JH)*AGESN(1))
      IF(KVEGA6.LE.0) THEN      ! 2 band
        BSNVIS=ASNVIS+ASNAGE
        BSNNIR=ASNNIR+ASNAGE
C**** Set zenith angle dependence if required
        IF (KKZSNO.GT.0) THEN
          CALL RXSNOW(BSNVIS,COSZ,GZSNOW(1,1,JH),XSNVIS)
          CALL RXSNOW(BSNNIR,COSZ,GZSNOW(7,1,JH),XSNNIR)
        ELSE
          XSNVIS=BSNVIS
          XSNNIR=BSNNIR
        END IF
      ELSE                      ! 6 band
        DO L=1,6
          FSNAGE=1.D0
          IF(L.GT.2) FSNAGE=2.0D0/L
          BSNVN(L)=ASNALB(L)+ASNAGE*FSNAGE
C**** Set zenith angle dependence if required
          IF (KKZSNO.GT.0) THEN
            CALL RXSNOW(BSNVN(L),COSZ,GZSNOW(L,1,JH),XSNVN(L))
          ELSE
            XSNVN(L)=BSNVN(L)
          END IF
        END DO
      ENDIF

C                                          -----------------------------
C                                          Soil/Veg Albedo Specification
C                                          -----------------------------
c**** In the following code when computing albedo for snow covered soil
c**** we use snow_frac(1:2) which is the snow fraction cover for
c**** bare/vegetated soil. It is computed in GHY_DRV.f in accordance
c**** with the surface topography.
c**** The final snow cover is minimum of snow_frac and the snow fraction
c**** obtained using the vegetation masking.
      DSFRAC=PVT(1)+PVT(10)
      VGFRAC=1.D0-DSFRAC
      IF(KVEGA6.LE.0) THEN                                      ! 2-band
      IF(SNOWE .LE.1.D-04) THEN
        BEAVIS=PVT(1)*ALBVNH(1,1,JH)*(1.D0-0.5D0*WEARTH*WETSRA)
        BEANIR=PVT(1)*ALBVNH(1,2,JH)*(1.D0-0.5D0*WEARTH*WETSRA)
        BVSOIL=BEAVIS
        BNSOIL=BEANIR
        DO K=2,NVEG
          BEAVIS=BEAVIS+PVT(K)*ALBVNH(K,1,JH)
          BEANIR=BEANIR+PVT(K)*ALBVNH(K,2,JH)
        END DO
        SEAVIS=BEAVIS
        SEANIR=BEANIR
        BVVEGE=BVSOIL
        BNVEGE=BNSOIL
        IF(VGFRAC.GT.0.001D0) THEN
          BVVEGE=(BEAVIS-BVSOIL*DSFRAC)/VGFRAC
          BNVEGE=(BEANIR-BNSOIL*DSFRAC)/VGFRAC
        ENDIF
      ELSE
        VTFRAC=PVT(1)*MAX((1.d0-snow_frac(1)),EXP(-SNOWE/VTMASK(1)))
        EXPSNE=VTFRAC +
     &       PVT(10)*MAX((1.d0-snow_frac(1)),EXP(-SNOWE/VTMASK(10)))
        DSFRAC=EXPSNE
        BEAVIS=VTFRAC*ALBVNH(1,1,JH)*(1.D0-0.5D0*WEARTH*WETSRA)
        BEANIR=VTFRAC*ALBVNH(1,2,JH)*(1.D0-0.5D0*WEARTH*WETSRA)
        DO K=2,NVEG
          VTFRAC=PVT(K)*MAX((1.d0-snow_frac(2)),EXP(-SNOWE/VTMASK(K)))
          BEAVIS=BEAVIS+VTFRAC*ALBVNH(K,1,JH)
          BEANIR=BEANIR+VTFRAC*ALBVNH(K,2,JH)
          EXPSNE=EXPSNE+VTFRAC
        END DO
      END IF
      XEAVIS=BEAVIS
      XEANIR=BEANIR
      BEAVIS=BEAVIS+BSNVIS*(1.D0-EXPSNE)
      BEANIR=BEANIR+BSNNIR*(1.D0-EXPSNE)
      XEAVIS=XEAVIS+XSNVIS*(1.D0-EXPSNE)
      XEANIR=XEANIR+XSNNIR*(1.D0-EXPSNE)
      VGFRAC=EXPSNE-DSFRAC
      XVSOIL=BVSOIL
      XNSOIL=BNSOIL
      XVVEGE=BVVEGE
      XNVEGE=BNVEGE
C
      ELSE                                                      ! 6-band
      IF(SNOWE .LE.1.D-04) THEN
        DO L=1,6
          BEAVN(L)=PVT(1)*ALBVNH(1,L,JH)*(1.D0-0.5D0*WEARTH*WETSRA)
        END DO
        BVSOIL=BEAVN(1)
        BNSOIL=BEAVN(2)
        DO K=2,NVEG
          DO L=1,6
            BEAVN(L)=BEAVN(L)+PVT(K)*ALBVNH(K,L,JH)
          END DO
        END DO
        SEAVIS=BEAVN(1)
        SEANIR=BEAVN(2)
        BVVEGE=BVSOIL
        BNVEGE=BNSOIL
        IF(VGFRAC.GT.0.001D0) THEN
          BVVEGE=(BEAVN(1)-BVSOIL*DSFRAC)/VGFRAC
          BNVEGE=(BEAVN(2)-BNSOIL*DSFRAC)/VGFRAC
        ENDIF
      ELSE
        VTFRAC=PVT(1)*MAX((1.d0-snow_frac(1)),EXP(-SNOWE/VTMASK(1)))
        EXPSNE=VTFRAC +
     &       PVT(10)*MAX((1.d0-snow_frac(1)),EXP(-SNOWE/VTMASK(10)))
        DSFRAC=EXPSNE
        DO L=1,6
          BEAVN(L)=VTFRAC*ALBVNH(1,L,JH)*(1.D0-0.5D0*WEARTH*WETSRA)
        END DO
        DO K=2,NVEG
          VTFRAC=PVT(K)*MAX((1.d0-snow_frac(2)),EXP(-SNOWE/VTMASK(K)))
          DO L=1,6
            BEAVN(L)=BEAVN(L)+VTFRAC*ALBVNH(K,L,JH)
          END DO
          EXPSNE=EXPSNE+VTFRAC
        END DO
      END IF
      DO L=1,6
        XEAVN(L)=BEAVN(L)
      END DO
      DO L=1,6
        BEAVN(L)=BEAVN(L)+BSNVN(L)*(1.D0-EXPSNE)
        XEAVN(L)=XEAVN(L)+XSNVN(L)*(1.D0-EXPSNE)
      END DO
      VGFRAC=EXPSNE-DSFRAC
      XVSOIL=BVSOIL
      XNSOIL=BNSOIL
      XVVEGE=BVVEGE
      XNVEGE=BNVEGE
      ENDIF                                                 ! end 6-band
C
      ITEA=TGE
      WTEA=TGE-ITEA
      ITEA=ITEA-ITPFT0
      BEASUM=0.D0
      BEAM=0.D0
      BEAP=0.D0
C
      DO K=1,NKBAND
        TRAPEA=AGSIDV(K,1)*(1.D0-EXPSNE)
     +        +AGSIDV(K,3)*DSFRAC*EDSTRA*(1.D0-WETTRA*WEARTH)
     +        +AGSIDV(K,4)*VGFRAC
        BEAM1 =(PLANCK(ITEA-1)-(PLANCK(ITEA-1)-PLANCK(ITEA  ))*WTEA)
     +       *(1.D0-TRAPEA)
        BEAM  =BEAM+BEAM1
        BEAP1 =(PLANCK(ITEA+1)-(PLANCK(ITEA+1)-PLANCK(ITEA+2))*WTEA)
     +       *(1.D0-TRAPEA)
        BEAP  =BEAP+BEAP1
        BEA   =(PLANCK(ITEA  )-(PLANCK(ITEA  )-PLANCK(ITEA+1))*WTEA)
     +       *(1.D0-TRAPEA)
        BEASUM=BEASUM+BEA
        ITEA=ITEA+ITNEXT
C
        TRGALB(K)=TRGALB(K)+PEARTH*TRAPEA
        BGFEMD(K)=BGFEMD(K)+PEARTH*(BEAP1-BEAM1)
        BGFEMT(K)=BGFEMT(K)+PEARTH*BEA
      END DO
      DTRUFG(2)=0.5D0*(BEAP-BEAM)
C
  500 CONTINUE
      IF(POICE.LT.1.D-04) GO TO 600
C
C                                         ------------------------------
C                                         Ocean Ice Albedo Specification
C                                         ------------------------------
      IF(KVEGA6.EQ.0 .or. KVEGA6.eq.3) then
C**** This albedo specification comes from Schramm et al 96 (4 spectral
C**** bands). Depending on KVEGA6 we either average to 2 or 6 bands
C**** Bare ice:
        if(hin.gt.0. .and. hin.lt.1.)then
          ali(1)=.76d0+.14d0*log(hin)
          ali(2)=.247d0+.029d0*log(hin)
          ali(3)=.055d0
          ali(4)=.036d0
        elseif (hin.ge.1. .and. hin.lt.2.) then
          ali(1)=.77d0+.018d0*(hin-1.)
          ali(2)=.247d0+.196d0*(hin-1.)
          ali(3)=.055d0
          ali(4)=.036d0
        elseif (hin.ge.2.) then
          ali(1)=.778d0
          ali(2)=.443d0
          ali(3)=.055d0
          ali(4)=.036d0
        endif
        albtf(1:4)=ali(1:4)
        albtr(1:4)=ali(1:4)
C**** Snow:
        patchy=0.
        if(hsn.gt.0.)then
          if(hsn.ge.0.1d0)then
            patchy=1d0
          else
            patchy=hsn/0.1d0
          endif
          if(flags)then         ! wet snow
            alsf(1)=.871d0
            alsf(2)=.702d0
            alsf(3)=.079d0
            alsf(4)=.001d0
            alsd(1:4)=alsf(1:4)
          else                  ! dry snow
            alsf(1)=.975d0
            alsf(2)=.832d0
            alsf(3)=.25d0
            alsf(4)=.025d0
            alsd(1)=.98d0-.008d0*cosz
            alsd(2)=.902d0-.116d0*cosz
            alsd(3)=.384d0-.222d0*cosz
            alsd(4)=.053d0-.0047d0*cosz
          endif
C**** consider snow aging based on Loth and Graf (1998)
C****  Dry, Wet(thick), Wet(thin) snow decreases by
C**** 0.006,  0.015 and 0.071 per day, respectively (for mean)
C**** assume decrease for each band is proportional
          if (flags) then
            if (hsn.gt.0.25) then
              snagfac = 0.015d0/0.7d0
            else
              snagfac = 0.071d0/0.7d0
            end if
          else
            snagfac = 0.006d0/0.82d0
          end if
C**** make sure it doesn't get too low!
          snagfac=min(snoage_fac_max,snagfac*AGESN(2))
          alsf(1:4)=alsf(1:4)*(1.-snagfac)
          alsd(1:4)=alsd(1:4)*(1.-snagfac)
C****
          albtf(1:4)=albtf(1:4)*(1.-patchy)+alsf(1:4)*patchy
          albtr(1:4)=albtr(1:4)*(1.-patchy)+alsd(1:4)*patchy
        endif
C**** Melt ponds:
        almp(1)=.15d0+exp(-8.1d0*hmp-.47d0)
        almp(2)=.054d0+exp(-31.8d0*hmp-.94d0)
        almp(3)=.033d0+exp(-2.6d0*hmp-3.82d0)
        almp(4)=.03d0
c**** combined sea ice albedo
        albtf(1:4)=albtf(1:4)*(1.-fmp)+almp(1:4)*fmp
C**** The zenith angle dependence from Schramm is only used is KKSNOW
C**** is zero. This only includes dry snow. Otherwise, Andy's
C**** calculation is used for all surface types (dry/wet snow, ice,
C**** meltponds)

C**** Zenith angle dependence for dry snow only
        albtr(1:4)=albtr(1:4)*(1.-fmp)+almp(1:4)*fmp
C**** Uncomment code below for zenith angle dependence for all types
C**** based on Dickinson (1981) (only possible if KKSNOW=0)
C****          a = a1                                   cosz>0.5
C****              a1 + (1-a1)*0.5 * (3/(1+4*cosz) -1) 0<cosz<.5
C**** ==> a_diff = 0.84 a1 + 0.16 (integrating over cosz)
C**** ==> a1= (a_diff-0.16)/0.84
c        if (cosz.gt.0.5) then
c          albtr(1:4) = (albtf(1:4)-0.16d0)/0.84d0
c        else
c          albtr(1:4) = (1.5*(albtf(1:4)-0.16d0)/0.84d0+1.-2.*cosz)/
c     *         (4.*cosz+1.)
c        end if

        IF(KVEGA6.GT.0) THEN
C**** 6 band albedo: map 4 Schramm wavelength intervals to 6 GISS ones
C**** Band#  range  %solar(grnd)  range   %solar(grnd)  band#
C****  (1)  250-690  (49.3%)   -->  300-770    (58.5%) (1)
C****  (2) 690-1190  (34.9%)   -->  770-860     (8.6%) (2)
C****                          -->  860-1250   (19.5%) (3)
C****  (3) 1190-2380 (14.8%)   --> 1250-1500    (3.8%) (4)
C****                          --> 1500-2200    (7.6%) (5)
C****  (4) 2380-4000 (1.0%)    --> 2200-4000    (2.0%) (6)
C**** Adjust weighting to force same broadband albedo
          BOIVN(1)=albtf(1)*.493d0/.585d0
          BOIVN(2:3)=albtf(2)*.349d0/.281d0
          BOIVN(4:5)=albtf(3)*.148d0/.114d0
          BOIVN(6)=albtf(4)*.01d0/.02d0

C**** set zenith angle dependence if required
          IF (KKZSNO.GT.0) THEN
            DO L=1,6
              CALL RXSNOW(BOIVN(L),COSZ,GZSNOW(L,2,JH),XOIVN(L))
            END DO
          ELSE ! use Schramm values (only for dry snow)
            XOIVN(1)=albtr(1)*.493d0/.585d0
            XOIVN(2:3)=albtr(2)*.349d0/.281d0
            XOIVN(4:5)=albtr(3)*.148d0/.114d0
            XOIVN(6)=albtr(4)*.01d0/.02d0
          END IF
        ELSE
C**** 2 band albedo: weight the 3 NIR bands by the solar irradiance to
C**** create a composite NIR value.
C**** Adjust weighting to force same broadband albedo
          BOIVIS=albtf(1)*.493d0/.585d0
          BOINIR=(.349d0*albtf(2)+.148d0*albtf(3)+.01d0*albtf(4))/.415d0
C**** set zenith angle dependence if required
          IF (KKZSNO.GT.0) THEN
            CALL RXSNOW(BOIVIS,COSZ,GZSNOW(1,2,JH),XOIVIS)
            CALL RXSNOW(BOINIR,COSZ,GZSNOW(7,2,JH),XOINIR)
          ELSE ! use Schramm values (only calculated for dry snow)
            XOIVIS=albtr(1)*.493d0/.585d0
            XOINIR=(.349d0*albtr(2)+.148d0*albtr(3)+.01d0*albtr(4))/
     *           .415d0
          END IF
        END IF
        EXPSNO=1.-patchy
C**** end of Schramm's version
      else
C**** original version

      EXPSNO=EXP(-SNOWOI/DMOICE)
C**** Set snow albedo over sea ice
      ASNAGE=ALBDIF(2,JH)*EXP(-AGEXPF(2,JH)*AGESN(2))
      IF(KVEGA6.LE.0) THEN      ! 2 band
        BSNVIS=ASNVIS+ASNAGE
        BSNNIR=ASNNIR+ASNAGE
      ELSE                      ! 6 band
        DO L=1,6
          FSNAGE=1.D0
          IF(L.GT.2) FSNAGE=2.0D0/L
          BSNVN(L)=ASNALB(L)+ASNAGE*FSNAGE
        END DO
      ENDIF

C**** set ice albedo
      IF (KVEGA6.le.0) THEN     ! 2 band
      BOIVIS=AOIVIS*EXPSNO+BSNVIS*(1.D0-EXPSNO)
      BOINIR=AOINIR*EXPSNO+BSNNIR*(1.D0-EXPSNO)
c**** Puddlings: weak in both Hemispheres, i.e. if Ts > 0C, then
c**** set albedos indep. of snow to .3/.15 up to .55/.3 as Ts grows
      if (kvega6.eq.-3) then
        IF(TSL.GT.273.16d0) THEN
          BOIVIS=.3d0           !   Ts > 10C
          BOINIR=.15d0
          IF(TSL.LT.283.16d0) THEN
            BOIVIS=AOIVIS-(TSL-273.16d0)*.1d0*(AOIVIS-.30d0) !  0<Ts<10C
            BOINIR=AOINIR-(TSL-273.16d0)*.1d0*(AOINIR-.15d0)
          END IF
        END IF
c**** Puddlings: weak in NH, strong (or extreme) in SH
      else if (kvega6.le.-4) then
        if (jlat.gt.jnorth) then    ! NH weak puddling
          IF(TSL.GT.273.16d0) THEN
            BOIVIS=.3d0                                         ! Ts>10C
            BOINIR=.15d0
            IF(TSL.LT.283.16d0) THEN
              BOIVIS=AOIVIS-(TSL-273.16d0)*.1d0*(AOIVIS-.30d0) !0<Ts<10C
              BOINIR=AOINIR-(TSL-273.16d0)*.1d0*(AOINIR-.15d0)
            END IF
          END IF
        else if (tgoi.gt.273.06.or.kvega6.eq.-5) then ! SH strong puddl
          BOIVIS=.25d0
          BOINIR=.1d0
        end if
      end if
c**** End of puddling section

C**** Set zenith angle dependence if required
        IF (KKZSNO.GT.0) THEN
          CALL RXSNOW(BOIVIS,COSZ,GZSNOW(1,2,JH),XOIVIS)
          CALL RXSNOW(BOINIR,COSZ,GZSNOW(7,2,JH),XOINIR)
        ELSE
          XOIVIS=BOIVIS
          XOINIR=BOINIR
        END IF
      ELSE                      ! end of 2-band original version
        DO L=1,6                !        6-band original version
          BOIVN(L)=AOIALB(L)*EXPSNO+BSNVN(L)*(1.D0-EXPSNO)
C**** Set zenith angle dependence if required
          IF (KKZSNO.GT.0) THEN
            CALL RXSNOW(BOIVN(L),COSZ,GZSNOW(L,2,JH),XOIVN(L))
          ELSE
            XOIVN(L)=BOIVN(L)
          END IF
        END DO
      ENDIF
      endif                     ! end of original version
C
      ITOI=TGOI
      WTOI=TGOI-ITOI
      ITOI=ITOI-ITPFT0
      BOISUM=0.D0
      BOIM=0.D0
      BOIP=0.D0
C
      DO K=1,NKBAND
        TRAPOI=AGSIDV(K,1)*ESNTRA*(1.-EXPSNO)
     +        +AGSIDV(K,2)*EICTRA*EXPSNO
        BOIM1 =(PLANCK(ITOI-1)-(PLANCK(ITOI-1)-PLANCK(ITOI  ))*WTOI)
     +       *(1.D0-TRAPOI)
        BOIM  =BOIM+BOIM1
        BOIP1 =(PLANCK(ITOI+1)-(PLANCK(ITOI+1)-PLANCK(ITOI+2))*WTOI)
     +       *(1.D0-TRAPOI)
        BOIP  =BOIP+BOIP1
        BOI   =(PLANCK(ITOI  )-(PLANCK(ITOI  )-PLANCK(ITOI+1))*WTOI)
     +       *(1.D0-TRAPOI)
        BOISUM=BOISUM+BOI
        ITOI=ITOI+ITNEXT
C
        TRGALB(K)=TRGALB(K)+POICE*TRAPOI
        BGFEMD(K)=BGFEMD(K)+POICE*(BOIP1-BOIM1)
        BGFEMT(K)=BGFEMT(K)+POICE*BOI
      END DO
      DTRUFG(3)=0.5D0*(BOIP-BOIM)
C
  600 CONTINUE
      IF(PLICE.LT.1.E-04) GO TO 700
C                                          -----------------------------
C                                          Land Ice Albedo Specification
C                                          -----------------------------
C**** Set snow albedo over land ice
      ASNAGE=ALBDIF(3,JH)*EXP(-AGEXPF(3,JH)*AGESN(3))
      IF(KVEGA6.LE.0) THEN      ! 2 band
        BSNVIS=ASNVIS+ASNAGE
        BSNNIR=ASNNIR+ASNAGE
        XSNVIS=BSNVIS
        XSNNIR=BSNNIR
      ELSE                      ! 6 band
        DO L=1,6
          FSNAGE=1.D0
          IF(L.GT.2) FSNAGE=2.0D0/L
          BSNVN(L)=ASNALB(L)+ASNAGE*FSNAGE
        END DO
        XSNVN(1:6)=BSNVN(1:6)
      ENDIF

      EXPSNL=EXP(-SNOWLI/DMLICE)
      IF (KVEGA6.le.0) THEN     ! 2 band
        BLIVIS=ALIVIS*EXPSNL+BSNVIS*(1.D0-EXPSNL)
        BLINIR=ALINIR*EXPSNL+BSNNIR*(1.D0-EXPSNL)
      ELSE                      ! 6-band
        DO L=1,6
          BLIVN(L)=ALIALB(L)*EXPSNL+BSNVN(L)*(1.D0-EXPSNL)
        END DO
      END IF
C**** For KVEGA6 != 1 or -1:
C**** Specify the Albedo for Antarctica and Greenland: vis.alb = 95%
C**** and mean albedo=80%, i.e. AMEAN = .585*BLIVIS+.415*BLINIR = .80
C****
      if (abs(kvega6).ne.1) then
      IF( JLAT.LT.NINT(MLAT46/6.) .OR.
     *   (JLAT.LT.45.AND.JLAT.GT.38.AND.ILON.LT.33.AND.ILON.GT.23)) THEN
        AMEAN=.8d0
        BLIVIS=.95d0
        BLINIR=(AMEAN-.585d0*BLIVIS)/.415d0
        IF (KVEGA6.gt.0) THEN   ! 6 band
          DO L=3,6  ! fill in higher bands
            BLIVN(L)=BLINIR
          END DO
        END IF
      END IF
      end if

C**** zenith angle dependence if required
      IF (KVEGA6.le.0) THEN     ! 2 band
        IF (KKZSNO.GT.0) THEN
          CALL RXSNOW(BLIVIS,COSZ,GZSNOW(1,3,JH),XLIVIS)
          CALL RXSNOW(BLINIR,COSZ,GZSNOW(7,3,JH),XLINIR)
        ELSE
          XLIVIS=BLIVIS
          XLINIR=BLINIR
        END IF
      ELSE  ! 6 band
        DO L=1,6
          IF (KKZSNO.GT.0) THEN
            CALL RXSNOW(BLIVN(L),COSZ,GZSNOW(L,3,JH),XLIVN(L))
          ELSE
            XLIVN(L)=BLIVN(L)
          END IF
        END DO
      END IF
C
      ITLI=TGLI
      WTLI=TGLI-ITLI
      ITLI=ITLI-ITPFT0
      BLISUM=0.D0
      BLIM=0.D0
      BLIP=0.D0
      BGF=0.D0
      DO K=1,NKBAND
        TRAPLI=AGSIDV(K,1)*ESNTRA*(1.-EXPSNL)
     +        +AGSIDV(K,2)*EICTRA*EXPSNL
        BLIM1 =(PLANCK(ITLI-1)-(PLANCK(ITLI-1)-PLANCK(ITLI  ))*WTLI)
     +       *(1.D0-TRAPLI)
        BLIM  =BLIM+BLIM1
        BLIP1 =(PLANCK(ITLI+1)-(PLANCK(ITLI+1)-PLANCK(ITLI+2))*WTLI)
     +       *(1.D0-TRAPLI)
        BLIP  =BLIP+BLIP1
        BLI   =(PLANCK(ITLI  )-(PLANCK(ITLI  )-PLANCK(ITLI+1))*WTLI)
     +       *(1.D0-TRAPLI)
        BLISUM=BLISUM+BLI
        ITLI=ITLI+ITNEXT
        TRGALB(K)=TRGALB(K)+PLICE*TRAPLI
        BGFEMD(K)=BGFEMD(K)+PLICE*(BLIP1-BLIM1)
        BGFEMT(K)=BGFEMT(K)+PLICE*BLI
      END DO
      DTRUFG(4)=0.5D0*(BLIP-BLIM)
  700 CONTINUE

C**** if 2-band, fill in rest of 6-band arrays
      IF(KVEGA6.LT.1) THEN                                      ! 2-band
        DO L=3,6
          BOCVN(L)=BOCVN(2)
          BEAVN(L)=BEAVN(2)
          BOIVN(L)=BOIVN(2)
          BLIVN(L)=BLIVN(2)
          XOCVN(L)=XOCVN(2)
          XEAVN(L)=XEAVN(2)
          XOIVN(L)=XOIVN(2)
          XLIVN(L)=XLIVN(2)
        END DO
      END IF
C**** write some BXA for diagnostic output (in WRITER) (replaces equiv.)
      BXA(1)=EXPSNE ; BXA(2)=EXPSNO ; BXA(3)=EXPSNL
      BXA(4)=BSNVIS ; BXA(5)=BSNNIR ; BXA(6)=XSNVIS ; BXA(7)=XSNNIR

C**** calculate final variables always over 6-bands
      DO L=1,6
        BVNSUR(L)=POCEAN*BOCVN(L)+PEARTH*BEAVN(L)
     +           + POICE*BOIVN(L)+ PLICE*BLIVN(L)
        XVNSUR(L)=POCEAN*XOCVN(L)+PEARTH*XEAVN(L)
     +           + POICE*XOIVN(L)+ PLICE*XLIVN(L)
      END DO
      DO L=1,6
        J=7-L
        PRNB(J,1)=BOCVN(L)
        PRNB(J,2)=BEAVN(L)
        PRNB(J,3)=BOIVN(L)
        PRNB(J,4)=BLIVN(L)
        PRNX(J,1)=XOCVN(L)
        PRNX(J,2)=XEAVN(L)
        PRNX(J,3)=XOIVN(L)
        PRNX(J,4)=XLIVN(L)
      END DO
      IF(KEEPAL.NE.1) THEN
        DO J=1,6
          L=7-J
          SRBALB(J)=BVNSUR(L)
          SRXALB(J)=XVNSUR(L)
        END DO
      ENDIF
C
C                     --------------------------------------------------
C                     Define each Surface Flux Factors, Flux Derivatives
C                     --------------------------------------------------
      BGF=0.D0
      DO K=1,NKBAND
        BGFEMD(K)=BGFEMD(K)*0.5D0
        BGF=BGF+BGFEMT(K)
      END DO
C
!obso BGM=BOCM*POCEAN+BEAM*PEARTH+BOIM*POICE+BLIM*PLICE
!obso BGP=BOCP*POCEAN+BEAP*PEARTH+BOIP*POICE+BLIP*PLICE
!obso TTRUFG=0.5D0*(BGP-BGM)
      FTRUFG(1)=BOCSUM/BGF
      FTRUFG(2)=BEASUM/BGF
      FTRUFG(3)=BOISUM/BGF
      FTRUFG(4)=BLISUM/BGF
C
      RETURN
      END SUBROUTINE GETSUR

      SUBROUTINE RXSNOW(RBSNO,XCOSZ,GGSNO,RXSNO)
!@sum RXSNOW calculate zenith angle dependence for snow/ice albedo
!@auth A. Lacis (modified by G. Schmidt)
      USE RE001, only : gtsalb,sgpgxg
      IMPLICIT NONE
!@var RBSNO diffuse albedo
      REAL*8, INTENT(IN) :: RBSNO
!@var XCOSZ zenith angle
      REAL*8, INTENT(IN) :: XCOSZ
!@var GGSNO Asymmetry parameter for snow
      REAL*8, INTENT(IN) :: GGSNO
!@var RXSNO direct albedo
      REAL*8, INTENT(OUT) :: RXSNO
      INTEGER NDBLS,NN
      REAL*8 XXG,XXT,GGSN,RBSN,FRTOP,TAU,TAUSN,GPFF,PR,PT,DBLS,SECZ,XANB
     *     ,XANX,TANB,TANX,RASB,RASX,BNORM,XNORM,RARB,RARX,XATB,DENOM,DB
     *     ,DX,UB,UX,DRBRAT,RBBOUT

      IF(RBSNO.LT.0.05D0) THEN
        RXSNO=RBSNO
        RETURN
      ENDIF
      call sys_flush(6)
      XXG=0.D0
      XXT=0.D0
      GGSN=GGSNO
      IF(GGSNO.GT.0.9D0) GGSN=0.9D0
      RBSN=RBSNO
      FRTOP=1.D0
      IF(RBSNO.GT.0.5D0) THEN
        RBSN=0.5D0
        FRTOP=((1.D0-RBSNO)/0.5D0)**2
      ENDIF

      CALL GTSALB(XXG,XXT,RBBOUT,RBSN,GGSN,TAUSN,2)
      CALL SGPGXG(XCOSZ,TAUSN,GGSN,GPFF)
      PR=1.D0-GPFF
      PT=1.D0+GPFF
      DBLS=10.D0+1.44269D0*LOG(TAUSN)
      NDBLS=DBLS
      TAU=TAUSN/2**NDBLS
C     Set optically thin limit values of R,T,X using PI0 renormalization
C     ------------------------------------------------------------------
C
      SECZ=1.D0/XCOSZ
      XANB=EXP(-TAU-TAU)
      XANX=EXP(-TAU*SECZ)
      TANB=PT*XANB
      XXT=(SECZ-2.D0)*TAU
      TANX=PT*SECZ
     +    *(.5D0+XXT*(.25D0+XXT*(.0833333D0+XXT*(.0208333D0+XXT))))*XANX
      RASB=PR*(1.D0-TAU*(2.D0-2.66667D0*TAU*(1.D0-TAU)))
      XXT=(SECZ+2.D0)*TAU
      RASX=PR*SECZ
     +    *(.5D0-XXT*(.25D0-XXT*(.0833333D0-XXT*(.0208333D0-XXT))))
      BNORM=(1.D0-XANB)/(RASB+TANB)
      XNORM=(1.D0-XANX)/(RASX+TANX)
      RASB=RASB*BNORM
      RASX=RASX*XNORM
      TANB=TANB*BNORM
      TANX=TANX*XNORM
      DO NN=1,NDBLS
        RARB=RASB*RASB
        RARX=XANX*RASX
        XATB=XANB+TANB
        DENOM=1.D0-RARB
        DB=(TANB+XANB*RARB)/DENOM
        DX=(TANX+RARX*RASB)/DENOM
        UB=RASB*(XANB+DB)
        UX=RARX+RASB*DX
        RASB=RASB+XATB*UB
        RASX=RASX+XATB*UX
        TANB=XANB*TANB+XATB*DB
        TANX=XANX*TANX+XATB*DX
        XANB=XANB*XANB
        XANX=XANX*XANX
      END DO
      DRBRAT=RASX/RBSN-1.D0
      RXSNO=RBSNO*(1.D0+DRBRAT*FRTOP)
      RETURN
      END SUBROUTINE RXSNOW

