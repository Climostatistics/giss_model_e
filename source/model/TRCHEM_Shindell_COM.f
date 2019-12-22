      MODULE TRCHEM_Shindell_COM
!@sum  TRCHEM_Shindell_COM declares variables for tracer chemistry
!@+    and sources.
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on various chemistry modules of B436Tds3YM23 model)
c
      USE RESOLUTION, only : im,jm,lm,ls1,istrat,psf,
     &                       ptop,sig,sige,dsig,bydsig
      USE MODEL_COM, only  : dtsrc,Itime,ItimeI,T,JEQ
      USE CONSTANT, only   : pi, mair, mwat, radian
      USE DYNAMICS, only   : am, byam, PMID, PK
      USE GEOM, only       : BYDXYP,dxyp
      USE RE001, only : TAUWC, TAUIC
      USE TRACER_COM, only : ntm, trm, TR_MM
c
      IMPLICIT NONE
      SAVE
c
C**************  P  A  R  A  M  E  T  E  R  S  *******************
C
!@param p_1 number of reactants per reaction
!@param p_2 unknown for now
!@param p_3 unknown for now
!@param p_4 unknown for now
!@param p_5 unclear what this is, for now'
!@param p_6 unclear what this is, for now'
!@param p_7 unclear what this is, for now'
!@param n_rx maximum number of chemical reactions
!@param n_bi maximum number of bimolecular reactions
!@param n_tri maximum number of trimolecular reactions
!@param n_bnd1 maximum number of spectral bands 1
!@param n_bnd2 maximum number of spectral bands 2
!@param n_bnd3 maximum number of spectral bands 3
!@param n_nst maximum number of reverse reactions
!@param n_fam maximum number of chemical families
!@param n_spc maximum number of chemical species
!@param n_oig max number of optically-important gases
!@param n_srb max number of schumann-runge bands
!@param N__ Number of levels in Mie grid: 2*(2*lpar+2+jaddto(1))+3
!@param M__ Number of Gauss points used
!@param NLFASTJ maximum number levels after inserting extra Mie levels
!@param NS maximum number of species which require J-values calculating
!@param NWFASTJ maximum number of wavelength bins that can be used
!@param JPPJ number of photolysis reactions
!@param JPNL number of photolysis levels
!@param NJVAL Number of species for which to calculate J-values
!@param szamax max Zenith Angle(98 deg at 63 km;99 degrees at 80 km)
!@param dtaumax max optical depth above which must instert new level
!@param ZZHT Scale height above top of atmosphere (cm)
!@param odmax Maximum allowed optical depth, above which they're scaled
!@param luselb Use reflective photolysis boundary treatment
!@param zlbatm Optical depth above which to set lower boundary
!@param NFASTJ number of quadrature points in OPMIE
!@param MFIT expansion of phase function in OPMIE
!@param MFASTJ ?
!@param CMEQ1 ?
!@param nc "number of molecules considered" ?
!@param ny "number of calculated gases" ?
!@param numfam number of families
!@param n_igas uknown: number of prod/dest gases?
!@param n_phot how often to do photolysis (in increments of DTsrc)
!@param O3MULT =2.14E-2    !is 1.E10*2.69E13*48./6.02E26
!@param BYO3MULT = 1/O3MULT
!@param pfix_O2 fixed ratio of O2/M
!@param pfix_H2 fixed ratio of H2/M
!@param pfix_Aldehyde fixed ratio of Aldehyde/M for initial conditions
!@param pfix_CH4_S fixed ratio of CH4/M in South. Hemis. (if used)
!@param pfix_CH4_S fixed ratio of CH4/M in South. Hemis. (if used)
!@param MWabyMWw ratio of molecular weights of air/water
!@param O3_1_fact factor to alter surface O3 that is passed to FASTJ
!@+     this is fastj level 1, not model level 1.  Currently, it is 
!@+     decreased by a factor of (972/1000)mb
!@param RGAMMASULF unknown
!@param RKBYPIM=8.*RBOLTZ/pi/MASSN2O55=8.*1.38062E-23/3.14159/1.793E-25
!@param cboltz Boltzman's Constant = 1.3806d-19
!@param dlogp 10.d0**(-2./16.)
!@param byradian 1/radian = conversion from radians to degrees
      INTEGER, PARAMETER ::
     & p_1   =     2,
     & p_2   =   111,
     & p_3   =   200,
     & p_4   =    70,
     & p_5   =    14,
     & p_6   =   116,
     & p_7   =    32,     
     & n_rx  =   101,
     & n_bi  =    96,
     & n_tri =    14,
     & n_bnd1=    31,
     & n_nst =     5,
     & n_fam =     4,
     & n_spc =    35,
     & n_bnd2=    87,
     & n_bnd3=   107,
     & n_oig =     3,
     & n_srb =    18,
C ----------------------------------------------     
c     & n_Ox=        1,    ! note, these
c     & n_NOx=       2,    ! first 15 species are
c     & n_N2O5=      3,    ! tracers, and therefore
c     & n_HNO3=      4,    ! these parameters are
c     & n_H2O2=      5,    ! be defined in 
c     & n_CH3OOH=    6,    ! TRACERS_Air.f.
c     & n_HCHO=      7,    ! Note the UNDERSCORE !!!
c     & n_HO2NO2=    8,    !  T
c     & n_CO=        9,    !  R
c     & n_CH4=      10,    !  A
c     & n_PAN=      11,    !  C
c     & n_Isoprene= 12,    !  E
c     & n_AlkylNit= 13,    !  R
c     & n_Alkenes=  14,    !  S
c     & n_Paraffin= 15,    ! ---------------
     & nC2O3=     16,
     & nXO2=      17,    
     & nXO2N=     18,
     & nRXPAR=    19,
     & nROR=      20,
     & nAldehyde= 21,
     & nH2O=      22,
     & nCH3O2=    23,
     & nH2=       24,
     & nOH=       25,
     & nHO2=      26,
     & nO3=       27,
     & nO=        28,
     & nO1D=      29,
     & nNO=       30,
     & nNO2=      31,
     & nNO3=      32,
     & nHONO=     33,
     & nO2=       34,
     & nM=        35,     !you must always put nM last (highest number)
C ----------------------------------------------   
     & N__=     1800,     !jan00, was 450, then 900 in Nov99
     & M__=        4,
     & NLFASTJ=  300,     !300 is arbitrary for now GSF
     & NS     =   51,
     & NWFASTJ=   15, 
     & JPPJ   =   16,
     & JPNL   =   12,
     & NJVAL  =   16,     !formerly read in from jv_spec00_15.dat
     & MFIT   =    8,
     & NFASTJ =    4,
     & MFASTJ =    1,
     & nc     =   35,     !formerly in param sub
     & ny     =   33,     !formerly in param sub  
     & numfam =    2,     !formerly in param sub  
     & n_igas =   43,
     & n_phot=     2      ! currently means 2 hours
C     
      REAL*8, PARAMETER ::  O3MULT       = 2.14D-2,
     &                      BYO3MULT     = 1./O3MULT,
     &                      pfix_O2      = 0.209476,
     &                      pfix_H2      = 560.E-9,
     &                      pfix_Aldehyde= 2.E-9,
     &                      pfix_CH4_S   = 1.75E-6,
     &                      pfix_CH4_N   = 1.855E-6,
     &                      MWabyMWw     = mair/mwat,
     &                      O3_1_fact    = 0.972,
     &                      RGAMMASULF   = 0.1,
     &                      RKBYPIM      = 1.961E2,
     &                      cboltz       = 1.3806d-19,
     &                      dlogp        = 7.49894209E-1, !=10^(-.125)
     &                      szamax       = 98.0d0,
     &                      dtaumax      = 1.0d0,
     &                      ZZHT         = 5.d5,
     &                      odmax        = 200.d0,
     &                      zlbatm       = 4.d0,
     &                      CMEQ1        = 0.25D0,
     &                      byradian     = 1.d0/radian
C  
      LOGICAL, PARAMETER :: luselb       = .false.
c
C**************  V  A  R  I  A  B  L  E  S *******************  
!@var nn reactant's number in mol list, first index reactant 1 or 2,
!@+      second - reaction number
!@var nnr reaction product's number in mol list, indicies like nn
!@var kss mollst number for product gases from photolysis
!@var ks mollst number for source gas in photolysis reaction
!@var nps reaction numbers by molecule, photolytic production
!@var nds reaction numbers by molecule, photolytic destruction
!@var npnr reaction numbers by molecule, photolytic production
!@var ndnr reaction numbers by molecule, photolytic destruction
!@var kps reaction numbers by molecule, chemical production
!@var kds reaction numbers by molecule, chemical destruction
!@var kpnr reaction numbers by molecule, chemical production
!@var kdnr reaction numbers by molecule, chemical destruction
!@var fam ___?
!@var nst reverse reaction number for dissociation reactions
!@var lbeg beginning of spectral range for photodissociation for ind.
!@+   gas (18 Sch-runge bands, then 200-730 nm in 5 nm steps)
!@var nir length of photodissociation spectra of gas (# of non-zero
!@+   absorption cross sections)
!@var lprn,jprn,iprn l, j, and i point for chemistry debugging
!@var ay name of gas being considered, in 2 four character strings
!@var y concentration of gas, 1st index=gas number, 2nd=verticle level
!@var rr rate constant of chemical reaction, first index - reaction
!@+   number, 2nd is i, 3rd is j, 4th is verticle level
!@var ss photodissociation coefficient, indexed as rr
!@var pe rate constant for bimolecular chemical reaction
!@var ea activation energy constant for bimolecular chemical reactions
!@var ro,r1,sn,sb rate parameters for trimolecular reactions
!@var sigg effective absorption cross section for radiation, first
!@+   index spectral interval number, second photolysis reaction number
!@var conc concentration of optically important gases (O2 & O3), first
!@+   vertivle level, second=gas number (1=O2,2=O3)
!@var qfu flux of solar radiation in the upper atmosphere in 18
!@+   schumann-runge bands 175 - 200 nm
!@var qf flux of solar radiation, first index for spectral interval
!@+   number, second - verticle level (photons/cm^2*c)
!@var wlt wavelength from 200 to 730 nm in 5 nm steps (107 intervals)
!@var sO3 absorption cross section of ozone (cm^2)
!@var sO2 absorption cross section of oxygen (cm^2)
!@var sech cross section of optically important gases, first index
!@+   gas number (1=O2,2=O3), second - spectral interval number 
!@var f8 flux of solar radiation outside the atmosphere, index
!@+   corresponds to the number of the spectral interval (cm-2*c-1)
!@var TXL temperature profile
!@var prnrts logical: print rate of each chemical reaction?
!@var prnchg logical: print chemical changes?
!@var prnls logical: print reaction lists by species?
!@var epot,egr,oz unclear for now
!@var yNO3,pHOx,pNOx,pOx,yCH3O2,yC2O3,yROR,yXO2,yAldehyde,yXO2N,yRXPAR?
!@var NCFASTJ number of levels in the fastj atmosphere
!@var title_aer_pf titles read from aerosol phase function file
!@var TITLE0 blank title read in I think
!@var TITLEJ titles read from O2, O3, and other species X-sections
!@var jlabel Reference label identifying appropriate J-value to use
!@var jind mapping index for jvalues? 
!@var jndlev Levels at which we want J-values (centre of CTM levels)
!@var jaddlv Additional levels associated with each level
!@var jaddto Cumulative total of new levels to be added
!@var NW1, NW2 beginning, ending wavelength for wavelength "bins"
!@var MIEDX Type of aerosol scattering, currently 6 set up:
!@+   1=Rayly 2=iso 30iso-equiv 4=bkgrd-sulf,5=volc-sulf,6=liq water
!@var NAA Number of categories for scattering phase functions
!@var npdep Number of pressure dependencies
!@var jpdep Index of cross sections requiring P dependence
!@var PFASTJ pressure sent to FASTJ
!@var nss this is a copy of JPPJ that is read in from a file
!@var jfacta Quantum yield (or multiplication factor) for photolysis
!@var WBIN Boundaries of wavelength bins
!@var WL Centres of wavelength bins - 'effective wavelength'
!@var NWWW Number of wavelength bins, from NW1:NW2
!@var FL Solar flux incident on top of atmosphere (cm-2.s-1)
!@var   Rayleigh parameters (effective cross-section) (cm2)
!@var QBC Black Carbon abs. extinct. (specific cross-sect.m2/g)
!@var QO2      O2 cross-sections
!@var QO3      O3 cross-sections
!@var Q1D      O3 => O(1D) quantum yield
!@var TQQ      Temperature for supplied cross sections
!@var QQQ      Supplied cross sections in each wavelength bin (cm2)
!@var QAAFASTJ Aerosol scattering phase functions
!@var NK       Number of wavelengths at which functions are supplied
!@var WAAFASTJ Wavelengths for the NK supplied phase functions
!@var PAA Scaling for extinctions
!@var zpdep    Pressure dependencies by wavelength bin
!@var lpdep    Label for pressure dependence
!@var OREF     O3 reference profile
!@var TREF     temperature reference profile
!@var BREF     black carbon reference profile
!@var U0 cosine of the solar zenith angle
!@var ZFASTJ Altitude of each fastj pressure level (approx.) (cm)
!@var RFLECT Surface albedo (Lamertian) in fastj
!@var odtmp Optical depth (temporary array)
!@var odsum Column optical depth
!@var nlbatm Level of lower photolysis boundary - usually surface ('1')
!@var aer obs fastj aerosol profile?
!@var O3J Ozone profile on photolysis grid
!@var TJ Temperature profile on photolysis grid
!@var DBC Mass of Black Carbon at each pressure level (g.cm-3)
!@var FFF Actinic flux at each level for each wavelength bin and level
!@var DMFASTJ Total number density at each pressure level (cm-3)
!@var DO3 Ozone number density at each pressure level (cm-3
!@var XQO2   Absorption cross-section of O2
!@var XQO3   Absorption cross-section of O3
!@var DTAUDZ   Local extinction at each point
!@var PIRAY    Contribution of Rayleigh scattering to extinction
!@var PIAER    Contribution of Aerosol scattering to extinction
!@var TTAU     Opt depth of air vert'y above each point(to top of atm)
!@var XLTAU    TTAU along the slant path
!@var FTAU     Attenuation of solar beam
!@var RZ      Distance from centre of Earth to each point (cm)
!@var RQ      Square of distance ratios
!@var TANHT   Tangent height for the current SZA
!@var XL      Slant path between points
!@var WTAU    Weighted slant path - each side of point
!@var dpomega   change in pomega per increment  (linear)
!@var POMEGAJ  Scattering phase function. the 2nd dimension on POMEGAJ
!@+   was 30 in the 9-layer code, which is (2*LM+2)+10. But it seems to
!@+   only need +1 (I.e. NCFASTJ+1): GSF
!@var POMEGA,ZTAU,FZ,ZREFL,jndlv,FJFASTJ,EMU,ZFLUX,ZREFL,ZU0,WFASTJ ?
!@var PM0,PM,BFASTJ,AFASTJ,AAFASTJ,WTFASTJ,CC,HFASTJ,C1,SFASTJ,U1,V1 ?
!@var RR2 former RR from fastj ?
!@var nfam number of beginning molecule of each chemical family
!@var ZJ photodissociation coefficient? (level,reaction)
!@var RCLOUDFJ cloudiness (optical depth) parameter, radiation to fastj
!@var SALBFJ surface albedo parameter from radiation to fastj
!@var O3DLJI model ozone to radiation (calculated)
!@var O3DLJI_clim model ozone to radiation (climatological)
!@var COlat carbon monoxide latitude distribution for init cond (ppbv)
!@var COalt adjustments of COlat by altitude (unitless)
!@var CH4FACT, FACTJ,r179m2v for setting CH4 ICs and strat distribution
!@var BYFJM = 1/JM
!@var avg67 average of lev 6 and 7 CH4 (used for IC and strat)
!@var mass2vol local array to convert between mass and volume units.
!@var bymass2vol local array to convert between mass and volume units.
!@var MODPHOT if MODPHOT=0 do photolysis, else skip it
!@var TX temperature variable for master chem
!@var ta, pres local arrays to hold temperature,pressure
!@var TFASTJ temperature sent to FASTJ
!@var O3_FASTJ ozone sent to fastj
!@var FASTJLAT,FASTJLON latitude & LONGITUDE (degrees) for use in fastj
!@var SZA the solar zenith angle (degrees)
!@var ILIMIT maximum of main I-loop (1 or IM)
!@var I,J,L,N,igas,inss,LL,Lqq,JJ,J3 dummy loop variables
!@var JFASTJ photolysis rates
!@var RVELN2O5, prod_sulf unknown
!@var sulfate N2O5 sulfate sink (formerly SRC(I,J,L,20) variable)   
!@var prod_sulfate  N2O5 change by sulfate reactions in mass units
!@var wprod_sulf N2O5 change by sulfate reactions in molecules/cm3/s
!@var DT2 variable chemical time step, set in masterchem
!@var nr total number of        reactions read in from gs_jpl00_trop_15
!@var nr3 #of trimolecular      reactions read in from gs_jpl00_trop_15
!@var nr2 #of mono+bi-molecular reactions read in from gs_jpl00_trop_15
!@var nmm #of monomolecular     reactions read in from gs_jpl00_trop_15
!@var nhet #of heterogenous     reactions read in from gs_jpl00_trop_15
!@var pfactor to convert units on species chemical changes
!@var bypfactor to convert units on species chemical changes
!@var dNO3,gwprodHNO3,gprodHNO3,gwprodN2O5,changeAldehyde,
!@+   changeAlkenes,changeIsoprene,changeHCHO,changeAlkylNit,
!@+   changeHNO3,changeNOx,changeN2O5,wprodHCHO working variables to 
!@+   calculate nighttime chemistry changes
!@var rlossN,rprodN,ratioN variables for nitrogen conservation
!@var change change due to chemistry in mass/time
!@var chemrate,photrate ?   
!@var MDOFM cumulative days at end of each month
!@var OxIC initial conditions for Ox in KG read from file (Jan 1st?)
!@var corrOx correction factor to tweak inital Ox in stratosphere 
!@var F0 reactivity factor for oxidation of biological substances
!@var HSTAR Henry's Law const for dry dep. why different from wet dep?
      INTEGER nr,nr2,nr3,nmm,nhet,MODPHOT,ILIMIT,igas,LL,I,J,L,N,
     & inss,Lqq,JJ,J3,lprn,jprn,iprn,NW1,NW2,MIEDX,NAA,npdep,nss,
     & NWWW,NK,nlbatm,NCFASTJ                                 
      INTEGER, DIMENSION(n_fam)        :: nfam
      INTEGER, DIMENSION(p_1,p_2)      :: nn, nnr, kss
      INTEGER, DIMENSION(p_2)          :: ks
      INTEGER, DIMENSION(p_3)          :: nps, nds, npnr, ndnr
      INTEGER, DIMENSION(p_4)          :: kps, kds, kpnr, kdnr
      INTEGER, DIMENSION(n_nst)        :: nst
      INTEGER, DIMENSION(n_bnd1)       :: lbeg, nir
      INTEGER, DIMENSION(LM)           :: jndlv,jndlev
      INTEGER, DIMENSION(JPPJ)         :: jind
      INTEGER, DIMENSION(NLFASTJ)      :: jaddlv, jaddto
      INTEGER, DIMENSION(NJVAL)        :: jpdep  
      INTEGER, DIMENSION(12)           :: MDOFM    
C      
      REAL*4, DIMENSION(2,n_spc)       :: ay
C      
      REAL*8 ZFLUX,ZREFL,ZU0,U0,RFLECT,odsum,XLTAU,TANHT,CH4FACT,FACTj,
     & r179m2v,BYFJM,FASTJLAT,FASTJLON,SZA,RVELN2O5,
     & prod_sulf,DT2,wprod_sulf,dNO3,gwprodHNO3,gprodHNO3,gwprodN2O5,
     & changeAldehyde,changeAlkenes,changeIsoprene,changeHCHO,
     & wprodHCHO,changeAlkylNit,changeHNO3,changeNOx,changeN2O5,
     & wprodCO,rlossN,rprodN,ratioN,pfactor,bypfactor
      REAL*8, DIMENSION(JM,4,12)       :: corrOx ! JM,(L=12,15),month
      REAL*8, DIMENSION(n_spc,LM)      :: y
      REAL*8, DIMENSION(n_rx,IM,JM,LM) :: rr, ss
      REAL*8, DIMENSION(n_bi)          :: pe, ea
      REAL*8, DIMENSION(n_tri)         :: ro, r1, sn, sb
      REAL*8, DIMENSION(n_bnd2,n_rx)   :: sigg
      REAL*8, DIMENSION(LM,n_oig)      :: conc
      REAL*8, DIMENSION(n_srb,p_5)     :: qfu
      REAL*8, DIMENSION(n_bnd3,LM)     :: qf
      REAL*8, DIMENSION(LM,n_bnd3)     :: epot, egr
      REAL*8, DIMENSION(n_bnd3)        :: wlt, sO3, sO2
      REAL*8, DIMENSION(n_oig,n_bnd3)  :: sech
      REAL*8, DIMENSION(p_6)           :: f8
      REAL*8, DIMENSION(p_7)           :: oz
      REAL*8, DIMENSION(IM,JM,LM)   :: yNO3,pHOx,pNOx,pOx,yCH3O2,yC2O3,
     &                                yROR,yXO2,yAldehyde,yXO2N,yRXPAR,
     &                                    RCLOUDFJ,TX,sulfate,OxIC
      REAL*8, DIMENSION(M__)         :: AFASTJ,C1,HFASTJ,V1,WTFASTJ,EMU
      REAL*8, DIMENSION(M__,M__)  :: BFASTJ,AAFASTJ,CC,SFASTJ,WFASTJ,U1
      REAL*8, DIMENSION(M__,2*M__)     :: PM
      REAL*8, DIMENSION(M__,M__,N__)   :: DD
      REAL*8, DIMENSION(M__,N__)       :: RR2    
      REAL*8, DIMENSION(2*M__)         :: PM0,dpomega
      REAL*8, DIMENSION(2*M__,N__)     :: POMEGA
      REAL*8, DIMENSION(N__)           :: ZTAU,FZ
      REAL*8, DIMENSION(2*M__,2*LM+2+1):: POMEGAJ
      REAL*8, DIMENSION(2*LM+2)        :: PFASTJ
      REAL*8, DIMENSION(NWFASTJ+1)     :: WBIN
      REAL*8, DIMENSION(NWFASTJ)       :: WL, FL, QRAYL, QBC
      REAL*8, DIMENSION(NWFASTJ,3)     :: QO3, QO2, Q1D, zpdep
      REAL*8, DIMENSION(3,NS)          :: TQQ
      REAL*8, DIMENSION(4,9)           :: QAAFASTJ, WAAFASTJ
      REAL*8, DIMENSION(8,4,9)         :: PAA
      REAL*8, DIMENSION(31,18,12)      :: OREF
      REAL*8, DIMENSION(41,18,12)      :: TREF
      REAL*8, DIMENSION(41)            :: BREF
      REAL*8, DIMENSION(NLFASTJ)       ::aer,ZFASTJ,O3J,TJ,DBC,DMFASTJ,
     &                                    XQO3,XQO2,DTAUDZ,TTAU,FTAU,
     &                                    PIAER,RZ,RQ,DO3,PIRAY
      REAL*8, DIMENSION(LM)            ::odtmp,TXL,COalt,ta,pres,TFASTJ
      REAL*8, DIMENSION(NS)            :: VALJ
      REAL*8, DIMENSION(N__)           :: FJFASTJ
      REAL*8, DIMENSION(NWFASTJ,2,NS-3):: QQQ
      REAL*8, DIMENSION(NWFASTJ,jpnl)  :: FFF
      REAL*8, DIMENSION(JPPJ)          :: jfacta
      REAL*8, DIMENSION(JPNL,JPPJ)      :: zj, JFASTJ
      REAL*8, DIMENSION(p_2,LM)         :: chemrate, photrate 
      REAL*8, DIMENSION(IM,JM)          :: SALBFJ
      REAL*8, DIMENSION(40,JM,IM)       :: O3DLJI, O3DLJI_clim
      REAL*8, DIMENSION(2*(LS1-1))      :: O3_FASTJ
      REAL*8, DIMENSION(19)             :: COlat
      REAL*8, DIMENSION(n_igas,LM)      :: dest, prod
      REAL*8, DIMENSION(NTM)            :: mass2vol,bymass2vol,F0,HSTAR
      REAL*8, DIMENSION(IM,JM,LM,ntm)   :: change
      REAL*8, DIMENSION(NLFASTJ,NLFASTJ):: WTAU
      REAL*8, DIMENSION(IM,(JM-3)/3)    :: avg67 
CGSF :: somehow need to do RCLOUDFJ(I,J,L)=TAUWC(L)+TAUIC(L) in RADIA
C     
      LOGICAL                         fam,prnrts,prnchg,prnls      
C      
      CHARACTER*20, DIMENSION(9)   :: title_aer_pf !formerly TITLEA( )
      CHARACTER*78                    TITLE0
      CHARACTER*7, DIMENSION(3,NS) :: TITLEJ
      CHARACTER*7, DIMENSION(JPPJ) :: jlabel
      CHARACTER*7, DIMENSION(3)    :: lpdep
C
c
C**************  D   A   T   A *******************  
C
      DATA EMU/.06943184420297D0,.33000947820757D0,.66999052179243D0,
     $          .93056815579703D0/
      DATA WTFASTJ/.17392742256873D0,.32607257743127D0,
     $          .32607257743127D0,.17392742256873D0/
      DATA nfam /27,30,0,0/ ! why 4?
CGSF [CO] ppbv based on 10deg lat-variation Badr & Probert 1994 fig 9:
      DATA COlat/40.,40.,40.,40.,45.,50.,60.,70.,80.,90.,110.
     &           ,125.,140.,165.,175.,180.,170.,165.,150./
CGSF multiplier of free trop [CO] by layer Badr & Probert 94 fig10&11,
CGSF Lopez-Valverde et al 93 fig 3, and Warneck 88, ch1 fig14.
      DATA COalt/2.,1.5625,1.375,1.25,1.125,1.0625,1.,1.,1.,
     &           1.,1.,.5,.375,.20,.20,.20,.20,.20,.25,.4,2.5,12.,60./
      DATA MDOFM/31,59,90,120,151,181,212,243,273,304,334,365/
      DATA F0/1.,0.1,0.,  0.,1.,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0./ 
C Note: it is not clear that these are the right HSTAR numbers to use.
C See the notes in DB396Tds3M23.f and the values in the file: 
C /u/cmrun/chem_files/henrys_law_tables.pdf  Example: why shouldn't
C HSTAR exactly match tr_RKD, the wet dep version of the Henry's
C Law constants? In any case, the ones below are in mole/(L atm). The
C conversion to mole/Joule, which I used in wet dep is:
C  mole/Joule = mole/(Pa*m3) X (101325. Pa/atm) X (1.E-6 m3/cm3)
C  X (1.E3 cm3 / L) = mole/(L*atm)           GSF 5/22/02
      DATA HSTAR/1.D-2,1.D-2,0.D0,  1.D5,1.D5,3.D2,  6.D3,0.D0,0.D0,
     &         0.D0,3.68D0,1.3D-2,  0.D0,0.D0,0.D0/    
C
      END MODULE TRCHEM_Shindell_COM
