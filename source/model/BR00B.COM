      IMPLICIT REAL*8(A-H,O-Z)

!@var LX max.number of rad. vertical layer edges = model_LM+4 or more !!
!@var MLAT46,MLON72 horizontal grid dimensions the radiation refers to
!@var JNORTH latitude index defining northern hemisphere : jlat>jnorth
      PARAMETER (LX=53+4, MLAT46=46,MLON72=72, JNORTH=MLAT46/2)

!@var iy1O3,MO3X first year, max.number of months of prescr. O3 history
      PARAMETER (iy1O3=1850, MO3X=12*(2050-iy1O3+1),
!@var iy1S0,MS0X first year, max.number of months of prescr. S0 history
!@var icycs0,mcycs0 solar cycle in yrs,months used to extend S0 history
     *           iy1S0=1882, MS0X=12*(1998-iy1S0+1),
     *           icycs0=11,  mcycs0=icycs0*12)

!@var nghg nr. of well-mixed GHgases: CO2 N2O CH4 CFC-11 CFC-12 others
      integer nghg,nyrsghg,ghgyr1,ghgyr2
      parameter (nyrsghg=2050-1850+1,nghg=6)
!@var nyrsghg max.number of years of prescr. greenhouse gas history
!@var ghgyr1,ghgyr2 first and last year of GHG history
!@var ghgam GHG-mixing ratios in ppm,ppm,ppm,ppb,ppb,ppb
      real*8 GHGAM(nghg,nyrsghg),XREF(nghg+1),XNOW(nghg+1)
      common/ghgcom/ghgyr1,ghgyr2,ghgam,xref,xnow

      REAL*4 TROAER,VDBCSU,TDUST,EPLMHC,UVLEAN
!offl REAL*4 FVEG11,FOLGIZ      !  needed only for offline use ->setIC
      REAL*4 O3CLIM

      COMMON/RADCOM_BLOCKD_CONTROL/  ! BLOCK DATA initialized parameters
!alb A              DMOICE,DMLICE,WETTRA,WETSRA,NV ! now in surf_albedo
!alb A             ,ASHZOI,ANHZOI ! not implemented; no ZOICE in ij_data
!alb A             ,FZSRA(6),FEMTRA(6)         ! gone except for ...(1):
     A              ZCLSRA,  ECLTRA
     B             ,EPSCON,RO3COL,FCLDTR,FCLDSR,PTLISO,TLGRAD,FGOLDH(13)
     C             ,FSXAER(5),FTXAER(5), TAUWC0,TAUIC0,FULGAS(13)
     D             ,S00WM2,RATLS0,S0,SSBTAU,FRAYLE,PIVMAX,FSPARE(998)
     E             ,KWVCON,KEEPAL,KEEPRH,KEEP10,KCNORM,KCLDEP,ICE012
     F             ,MRELAY,MOZONE,KO3LON,NO3COL,NORMS0,KSOLAR,KZSNOW
     G             ,KTREND,NTRACE,ITR(8),LASTVC,ISPARE(98)
     H             ,MADO3M,MADAER,MADDST,MADVOL,MADEPS,MADLUV,MADGHG
     I             ,MADBAK,ITPFT0,ITNEXT,MADSUR         ! ,JYEAR,JDAY

      COMMON/RADCOM_INPUT_IJDATA/    !              Input data to RCOMPX
     A              PLB(LX),HLB(LX),TLB(LX),TLT(LX),TLM(LX),ULGAS(LX,12)
     B             ,TAUWC(LX),TAUIC(LX),SIZEWC(LX),SIZEIC(LX),CLDEPS(LX)
     C             ,SHL(LX),RHL(LX),TRACER(LX,8),SRBXAL(15,2)
     D             ,PVT(11),AGESN(3),SNOWE,SNOWOI,SNOWLI,WEARTH
     E             ,POCEAN,PEARTH,POICE,PLICE,TGO,TGE,TGOI,TGLI,TSL,COSZ
!alb  writer-arrays       ,BXA(7),PRNB(6,4),PRNX(6,4) replace BXA(153)
     F             ,WMAG  ,BXA(7),PRNB(6,4),PRNX(6,4)
!alb          new variables mostly for alt. sea ice albedo
     X             ,hsn,hin,hmp,fmp,snow_frac(2),lkdepth, PLAKE
C             integer variables start here, followed by logicals
     Y             ,JLAT,ILON,NL,NLP, LS1_loc,flags ! the last 2 are new
!@var LS1_loc local tropopause level
C**** Additional info passed for new prescription of sea ice albedo
!@var hsn depth of snow (m)
!@var hin depth of ice (m)
!@var hmp depth of melt pond (m)
!@var fmp fraction of melt pond area
!@var lkdepth lake depth (m)
      REAL*8 :: hsn,hin,hmp,fmp,lkdepth
!@var flags true is snow is wet
      LOGICAL*4 :: flags
!$OMP  THREADPRIVATE(/RADCOM_INPUT_IJDATA/)

C     INPUT data ! some local, some not: repeat this section in driver
      REAL*8 U0GAS
      COMMON/RADCOM_hybrid/U0GAS(LX,12)
!$OMP  THREADPRIVATE(/RADCOM_hybrid/)
C     end of section to be repeated in driver (needed for 'copyin')

C     OUTPUT DATA
      COMMON/RADCOM_OUTPUT_IJDATA/
     A              TRDFLB(LX),TRUFLB(LX),TRNFLB(LX),TRFCRL(LX)
     B             ,SRDFLB(LX),SRUFLB(LX),SRNFLB(LX),SRFHRL(LX)
     C             ,SRIVIS,SROVIS,PLAVIS,SRINIR,SRONIR,PLANIR,SRXATM(4)
     D             ,SRDVIS,SRUVIS,ALBVIS,SRDNIR,SRUNIR,ALBNIR,FSRNFG(4)
     E             ,SRTVIS,SRRVIS,SRAVIS,SRTNIR,SRRNIR,SRANIR,FTRUFG(4)
     F             ,TRDFGW,TRUFGW,TRUFTW,BTEMPW,DTRUFG(4) ! ,TTRUFG gone
     G             ,WINDZF(3),WINDZT(3),TOTLZF(3),TOTLZT(3),SRKINC(16)
     I             ,SRKALB(16),SRKGAX(16,4),SRKGAD(16,4),SKDFLB(LX,17)
     J             ,SKUFLB(LX,17),SKNFLB(LX,17),SKFHRL(LX,17)
!obsoK             ,FTAUSL(33),TAUSL(33),TRSLTS,TRSLTG,TRSLWV,TRSLBS
!obsoK             ,TRDFSL,TRUFSL,TRSLCR,SRSLHR
     L             ,LBOTCL,LTOPCL
!$OMP  THREADPRIVATE(/RADCOM_OUTPUT_IJDATA/)

      COMMON/RADCOM_BLOCKD_DATA/     !       BLOCK DATA initialized data
     A              AGOLDH(11,5),BGOLDH(11,5),CGOLDH(11,5)
!alb  ->surf_albedo AGSIDV,AOCEAN,ASN/AOI/ALIALB,ALBVND,V6ALB,GZSNOW
     C             ,CMANO2(42),PPMV80(12),Z0(12),ZH(12)
     D             ,SRCGSF(3),DLAT46(46),DLON72(72)
     E             ,WTHEK(190),FTHEK(190)    ,WSOLAR(190),FSOLAR(190) !?
     F   ,WSLEAN(190),FSLEAN(190),DSLEAN(190),W1LEAN(190),FRLEAN(190) !?
     G                         ,XWAVO3(226),FUVKO3(226),UVA(226)
     H           ,REAERO(10),VEAERO(10),ROAERO(10),FSAERO(10),FTAERO(10)
     I             ,FRSULF(10),REDUST(8),VEDUST(8),RODUST(8),FSDUST(8)
     J             ,FTDUST(8),VDGAER(12,10),PI0MAX(10),PIDUST(8)

      COMMON/WORKDATA/          !          Temp data generated by RCOMPX
     A              SRAEXT(LX,6),SRASCT(LX,6),SRAGCB(LX,6),TRCALK(LX,33)
     B             ,SRBEXT(LX,6),SRBSCT(LX,6),SRBGCB(LX,6),TRAALK(LX,33)
     C             ,SRDEXT(LX,6),SRDSCT(LX,6),SRDGCB(LX,6),TRBALK(LX,33)
     D             ,SRVEXT(LX,6),SRVSCT(LX,6),SRVGCB(LX,6),TRTAUK(LX,33)
     E             ,DBLEXT(LX,6),DBLSCT(LX,6),DBLGCB(LX,6),DBLPI0(LX,6)
     F             ,SRCEXT(LX,6),SRCSCT(LX,6),SRCGCB(LX,6),SRCPI0(LX,6)
     G             ,TRDALK(LX,33),TRVALK(LX,33),TRGXLK(LX,33),TRCTCA(33)
     H             ,DFLB(LX,33),UFLB(LX,33),WFLB(LX,33),PL(LX),DPL(LX)
     I             ,TRGALB(LX),BGFEMT(LX),BGFEMD(LX),WTLB(LX),WTLT(LX)
     J             ,ENA(LX),ENB(LX),ENC(LX),TRA(LX),TRB(LX),TRC(LX)
     K             ,DFSL(33),UFSL(33),WFSL(33),ITLB(LX),ITLT(LX)
     L             ,AERX1(LX),AERS1(LX),AERG1(LX),TRAXNL(LX)
     M             ,AERX2(LX),AERS2(LX),AERG2(LX),UGAS0(LX),UGASR(LX)
     N             ,RNB(LX),RNX(LX),TNB(LX),TNX(LX),XNB(LX),XNX(LX)
     O             ,SRB(LX),SRX(LX),VRU(LX),VRD(LX),FAC(LX)
     P             ,UXGAS(LX,9),TAUN(33*LX)
!alb P             ,BVSURF,BNSURF,XVSURF,XNSURF -> surf_albedo
     Q             ,DNA(LX),DNB(LX),DNC(LX),Q55H2S
!
!     COMMON/WORKDAT2/          !          Temp data generated by RCOMPX
     A             ,RIJTCK(6,33),FDXTCK(3,33),ALBTCK(3,33),CLPI0(33)
     B             ,FEMTCK(3,33),TXCTPG(33),TSCTPG(33),TGCTPG(33)
     C             ,QAERO(LX,6),SAERO(LX,6),CAERO(LX,6),AAERO(LX,33)
     D             ,QDUST(LX,6),SDUST(LX,6),CDUST(LX,6),ADUST(LX,33)
     E             ,O2FHRL(LX),SRAXNL(LX),SRASNL(LX),SRAGNL(LX),AO3X(LX)
     F             ,O2FHRB(LX),AO3D(LX),AO3U(LX)  ! ,ALBVNH(11,6,2)->alb
     G             ,HTPROF(LX),QVH2S(6),SVH2S(6),GVH2S(6),AVH2S(33)
!$OMP  THREADPRIVATE(/WORKDATA/)

      COMMON/RADINFO_WRITER_DATA/   !   Temp data used by WRITER, WRITET
     A              SRCQPI( 6,15),TRCQPI(33,15)
     C         ,TRAQAB(33,11),TRBQAB(33, 9),TRCQAB(33,15),TRDQAB(33,25)
     E             ,SKDKS0(16),NORDER(16),NMWAVA(16),NMWAVB(16)

      COMMON/RADDAT_SETUP1_GENERATED   ! Model files generated by RCOMP1
     H             /HLB0(LX),PLB0(LX),TLM0(LX),U0GAS3(LX)
     A             ,TKPFW(630),TKPFT(900),AO3(460),O3DLJ(LX,46)
     D             ,O3DLJM(LX,46,12),FPXCO2(LX),FPXOZO(LX),PIAERO(10)
     C             ,SRAX(LX,6,5),SRAS(LX,6,5),SRAC(LX,6,5),ZTABLE(LX,11)
     D             ,QXAERO(6,10),QSAERO(6,10),QCAERO(6,10),ATAERO(33,10)
     E             ,QXDUST(6,8),QSDUST(6,8),QCDUST(6,8),ATDUST(33,8)
     D             ,QDST55(8),ZO3(32),TRAX(LX,33,5)
     E             ,DKS0(16),DBLN(30),TCLMIN,KSLAM(16),NKSLAM
!    I             ! ,psig(lx) -> setatm

      COMMON/RADDAT_TR_SGP_TABLES                  !  radfile1, radfile2
     A             /TAUTBL(148000),PLANCK(8250),XKCFC(12,8,19,4)
     B             ,TAUWV0(148000),H2OCN8(33,8,14),H2OCF8(33,8,5)
     C             ,DUCH4(150),SDUCH4(150),DUN2O(150),SDUN2O(150)
     D             ,ULOX(247),DUX(247),GTAU(51,11,143),TGDATA(122,13)
     E           ,SALBTG(768,14),TAUGSA(1001,14),TAUTGD(122),TAUTGS(768)

      COMMON/RADDAT_OZONE_MCPETERS             !                radfile4
     A             /O3DPPM(31,16,12),PO3(32),O3AVE(12,18,18),AO3AV1(216)
     B             ,SO3JFS(11,19,2),TOPO3(30,46,12),SAGEZP(50,12,12)
     C             ,SZP50(50,46,12),O3NCAR(72,46,12)

      COMMON/RADDAT_AERCLD_MIEPAR              !                radfile3
     A             /SRAQEX( 6,11),SRAQSC( 6,11),SRAQCB( 6,11),Q55A11(11)
     B             ,TRAQEX(33,11),TRAQSC(33,11),TRAQCB(33,11),REFA11(11)
     C             ,SRBQEX( 6,10),SRBQSC( 6,10),SRBQCB( 6,10),Q55B10(10)
     D             ,TRBQEX(33,10),TRBQSC(33,10),TRBQCB(33,10),REFB10(10)
     E             ,SRCQEX( 6,15),SRCQSC( 6,15),SRCQCB( 6,15),Q55C15(15)
     F             ,TRCQEX(33,15),TRCQSC(33,15),TRCQCB(33,15),REFC15(15)
     G             ,TRCQAL(33,15),VEFC15(15)   ,VEFA11(   11),VEFB10(10)
     H             ,SRDQEX( 6,25),SRDQSC( 6,25),SRDQCB( 6,25),Q55D25(25)
     I             ,TRDQEX(33,25),TRDQSC(33,25),TRDQCB(33,25),REFD25(25)
     J             ,TRDQAL(33,25),VEFD25(25),SIZENR(183),XNR(13)
     K     ,SRVQEX( 6,20,6),SRVQSC( 6,20,6),SRVQCB( 6,20,6),Q55V20(20,6)
     L     ,TRVQEX(33,20,6),TRVQSC(33,20,6),TRVQCB(33,20,6),REFV20(20,6)
     M     ,TRVQAL(33,20,6),VEFV20(20,6),QEXTNR(183,13),COSBNR(183,13)
     N         ,SRUQEX( 6,120),SRUQSC( 6,120),SRUQCB( 6,120),Q55U22(120)
     O         ,TRUQEX(33,120),TRUQSC(33,120),TRUQCB(33,120),REFU22(120)
     P         ,TRUQAL(33,120),VEFU22(120),TRSQAL(33,25),VEFS25(25)
     Q             ,SRSQEX( 6,25),SRSQSC( 6,25),SRSQCB( 6,25),Q55S25(25)
     R             ,TRSQEX(33,25),TRSQSC(33,25),TRSQCB(33,25),REFS25(25)

      COMMON/RADDAT_CLDCOR_TRSCAT              !                radfileE
     A         /RIJTPG(6,49,17,21),FDXTPG(3,49,17,21),FEMTPG(3,49,17,21)

      COMMON/RADMAD1_OZONE_DECADAL/  !  (user SETO3D)  radfileA,radfileB
     A        O3CLIM(MO3X,LX,46),WJ1890(72,46,LX,12),WJ1979(72,46,LX,12)
     B       ,O3WJA(72,46,LX,12),O3LF(72,46,LX)

      COMMON/RADMAD2_TROPAER_DECADAL/      !  (user SETAER)     radfile5
     A               TROAER(72,46,13,8,4),VDBCSU(46,12,13,3)
     B              ,TAUCOL(72,46,10),VDAERO(46,12,10)

      COMMON/RADMAD3_DUST_SEASONAL/        !  (user SETDST)     radfile6
     A            TDUST(72,46,9,8,13),DUSTLN(72,46,12,8),DSTCOL(72,46,8)

      COMMON/RADMAD4_VOLCAER_DECADAL/      !  (user SETVOL)     radfile7
     A               V4TAUR(1800,24,5),FDATA(80),GDATA(80)
     C              ,HTFLAT(49,4),TAULAT(49),SIZLAT(49)

      COMMON/RADMAD5_CLDEPS_3D_SEASONAL/   !  (user SETCLD)     radfile8
     A               EPLMHC(72,46,12,4),EPLOW(72,46),EPMID(72,46)
     B              ,EPHIG(72,46),EPCOL(72,46)

      COMMON/RADMAD6_SOLARUV_DECADAL/      !  (user SETSOL)     radfile9
     A               UVLEAN(Ms0X,190),TSI1(Ms0X),TSI2(Ms0X)

!offl COMMON/RADMAD7_VEGTOP/       !    (user SETSUR)  radfileC,radfileD
!offlC               FVEG11(72,46,11),FOLGIZ(72,46,9)

      COMMON/RADMAD8_RELHUM_AERDATA/       !  (user SETREL)     radfileH
     A TRHQAB(33,168,4),SRHQEX(6,168,4),SRHQSC(6,168,4),SRHQCB( 6,168,4)


C     EQUIVALENCE (FSPARE(13),PSOIL)
      EQUIVALENCE (FSPARE(14),REFF0)
      EQUIVALENCE (FSPARE(15),VEFF0)
      EQUIVALENCE (FSPARE(16),O3WJT0)
      EQUIVALENCE (FSPARE(18),X0YBCI)
      EQUIVALENCE (FSPARE(19),X0YOCI)
      EQUIVALENCE (FSPARE(20),X0YSUI)
      DIMENSION UVWAVL(3)
      EQUIVALENCE (FSPARE(25),UVWAVL(1)) ! -FSPARE(27)
      DIMENSION UVFACT(3)
      EQUIVALENCE (FSPARE(28),UVFACT(1)) ! -FSPARE(30)
      EQUIVALENCE (FSPARE(31),PTOPTR)
      DIMENSION PI0VIS(11)
      EQUIVALENCE (FSPARE(41),PI0VIS(1)) ! -FSPARE(51)
      DIMENSION PI0TRA(11)
      EQUIVALENCE (FSPARE(61),PI0TRA(1)) ! -FSPARE(71)
      EQUIVALENCE (FSPARE(72),TSOIL)
      EQUIVALENCE (FSPARE(73),TVEGE)
      EQUIVALENCE (FSPARE(74),FTRUFS)
      EQUIVALENCE (FSPARE(75),FTRUFV)
      EQUIVALENCE (FSPARE(76),DTRUFS)
      EQUIVALENCE (FSPARE(77),DTRUFV)
      DIMENSION VDFAER(20,3)
      EQUIVALENCE (FSPARE(81),VDFAER(1,1)) ! -FSPARE(140)
      DIMENSION VDFBCI(12)
      EQUIVALENCE (FSPARE(81),VDFBCI(1))   ! -FSPARE(92)
      DIMENSION VDFOCI(12)
      EQUIVALENCE (FSPARE(101),VDFOCI(1))  ! -FSPARE(112)
      DIMENSION VDFSUI(12)
      EQUIVALENCE (FSPARE(121),VDFSUI(1))  ! -FSPARE(132)
      DIMENSION VDFDST(12)
      EQUIVALENCE (FSPARE(141),VDFDST(1))  ! -FSPARE(152)
      DIMENSION PPMVK0(12)
      EQUIVALENCE (FSPARE(161),PPMVK0(1))  ! -FSPARE(162)
      DIMENSION PPMVDF(12)
      EQUIVALENCE (FSPARE(181),PPMVDF(1))  ! -FSPARE(192)
      DIMENSION PPGRAD(12)
      EQUIVALENCE (FSPARE(201),PPGRAD(1))  ! -FSPARE(212)

      EQUIVALENCE (ISPARE( 1),KVRAER)
      EQUIVALENCE (ISPARE( 2),MEANAC)
      EQUIVALENCE (ISPARE( 3),MEANDD)
      EQUIVALENCE (ISPARE( 4),MEANVA)
      EQUIVALENCE (ISPARE( 5),NCARO3)
      EQUIVALENCE (ISPARE( 6),KUVFAC)
      EQUIVALENCE (ISPARE( 7),KSNORM)
      EQUIVALENCE (ISPARE( 8),KWTRAB)
      EQUIVALENCE (ISPARE( 9),KGGVDF)
      EQUIVALENCE (ISPARE(10),KPGRAD)
      EQUIVALENCE (ISPARE(11),KLATZ0)
      EQUIVALENCE (ISPARE(12),KCLDEM)
      EQUIVALENCE (ISPARE(13),KANORM)
      EQUIVALENCE (ISPARE(14),KPFCO2)
      EQUIVALENCE (ISPARE(15),KPFOZO)
      EQUIVALENCE (ISPARE(16),KVEGA6)
      EQUIVALENCE (ISPARE(17),KORDER)
      EQUIVALENCE (ISPARE(18),KUFH2O)
      EQUIVALENCE (ISPARE(19),KUFCO2)
C     EQUIVALENCE (ISPARE(20),??????)    ! unused
      EQUIVALENCE (ISPARE(21),KCSELF)
      EQUIVALENCE (ISPARE(22),KCFORN)
C     EQUIVALENCE (ISPARE(23-30),??????) ! unused
      EQUIVALENCE (ISPARE(31),KYEARS)
      EQUIVALENCE (ISPARE(32),KYEARG)
      EQUIVALENCE (ISPARE(33),KYEARO)
      EQUIVALENCE (ISPARE(34),KYEARA)
      EQUIVALENCE (ISPARE(35),KYEARD)
      EQUIVALENCE (ISPARE(36),KYEARV)
      EQUIVALENCE (ISPARE(37),KYEARE)
      EQUIVALENCE (ISPARE(38),KYEARR)
      EQUIVALENCE (ISPARE(41),KJDAYS)
      EQUIVALENCE (ISPARE(42),KJDAYG)
      EQUIVALENCE (ISPARE(43),KJDAYO)
      EQUIVALENCE (ISPARE(44),KJDAYA)
      EQUIVALENCE (ISPARE(45),KJDAYD)
      EQUIVALENCE (ISPARE(46),KJDAYV)
      EQUIVALENCE (ISPARE(47),KJDAYE)
      EQUIVALENCE (ISPARE(48),KJDAYR)


      DIMENSION AO3AVE(18,12)
      EQUIVALENCE (AO3AV1(1),AO3AVE(1,1))
      DIMENSION SO3JF(11,19)
      EQUIVALENCE (SO3JFS(1,1,1),SO3JF(1,1))
      DIMENSION SO3SO(11,19)
      EQUIVALENCE (SO3JFS(1,1,2),SO3SO(1,1))

      EQUIVALENCE (FSXAER(1),FSTAER),  (FTXAER(1),FTTAER)
      EQUIVALENCE (FSXAER(2),FSBAER),  (FTXAER(2),FTBAER)
      EQUIVALENCE (FSXAER(3),FSAAER),  (FTXAER(3),FTAAER)
      EQUIVALENCE (FSXAER(4),FSDAER),  (FTXAER(4),FTDAER)
      EQUIVALENCE (FSXAER(5),FSVAER),  (FTXAER(5),FTVAER)

      DIMENSION SRQV( 6,20),SRSV( 6,20),SRGV( 6,20),Q55V(   20),REFV(20)
      DIMENSION TRQV(33,20),TRSV(33,20),TRGV(33,20),TRAV(33,20),VEFV(20)
      EQUIVALENCE (SRVQEX(1,1,6),SRQV(1,1)),(SRVQSC(1,1,6),SRSV(1,1))
      EQUIVALENCE (SRVQCB(1,1,6),SRGV(1,1)),(Q55V20(1,6),Q55V(1))
      EQUIVALENCE (TRVQEX(1,1,6),TRQV(1,1)),(TRVQSC(1,1,6),TRSV(1,1))
      EQUIVALENCE (TRVQCB(1,1,6),TRGV(1,1)),(TRVQAL(1,1,6),TRAV(1,1))
      EQUIVALENCE (REFV20(1,6),REFV(1)),(VEFV20(1,6),VEFV(1))

      DIMENSION    SRBALB(6),SRXALB(6)
      EQUIVALENCE (SRBXAL(1,1),SRBALB(1)),(SRBXAL(1,2),SRXALB(1))

      EQUIVALENCE (SRXATM(1),SRXVIS),(SRXATM(2),SRXNIR)
      EQUIVALENCE (SRXATM(3),XXAVIS),(SRXATM(4),XXANIR)

