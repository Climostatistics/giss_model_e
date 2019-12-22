#include "rundeck_opts.h"

      MODULE CLOUDS
!@sum  CLOUDS column physics of moist conv. and large-scale condensation
!@auth M.S.Yao/A. Del Genio (modifications by Gavin Schmidt)
!@ver  1.0 (taken from CB265)
!@cont MSTCNV,LSCOND
      USE CONSTANT, only : rgas,grav,lhe,lhs,lhm,sha,bysha
     *     ,by3,tf,bytf,rvap,bygrav,deltx,bymrat,teeny,gamd
      USE MODEL_COM, only : im,lm,dtsrc
      USE QUSDEF, only : nmom,xymoms,zmoms,zdir
#ifdef TRACERS_ON
      USE TRACER_COM, only: ntm
#ifdef TRACERS_WATER
     & ,nGAS, nPART, nWATER, tr_wd_TYPE, tr_RKD, tr_DHD, tr_evap_fact
#endif
#endif
CCC   USE RANDOM
      IMPLICIT NONE
      SAVE
C**** parameters and constants
      REAL*8, PARAMETER :: TI=233.16d0   !@param TI pure ice limit
!@param WMU critical cloud water content for rapid conversion (g m**-3)
      REAL*8, PARAMETER :: WMU=.25
      REAL*8, PARAMETER :: WMUL=.5       !@param WMUL WMU over land
      REAL*8, PARAMETER :: WMUI=.1d0     !@param WMUI WMU for ice clouds
      REAL*8, PARAMETER :: BRCLD=.2d0    !@param BRCLD for cal. BYBR
      REAL*8, PARAMETER :: SLHE=LHE*BYSHA
      REAL*8, PARAMETER :: SLHS=LHS*BYSHA
!@param FCLW fraction of condensate in plume that remains as CLW
!@param CCMUL multiplier for convective cloud cover
!@param CCMUL1 multiplier for deep anvil cloud cover
!@param CCMUL2 multiplier for shallow anvil cloud cover
!@param COETAU multiplier for convective cloud optical thickness
      REAL*8, PARAMETER :: FCLW=0.5
      REAL*8, PARAMETER :: CCMUL=1.,CCMUL1=5.,CCMUL2=3.,COETAU=.08d0
!@var RcldXL layer-dep. adj.factors for ice cloud particle size (1)
      real*8, dimension(lm) :: RcldXL

      REAL*8 :: BYBR,BYDTsrc,XMASS
!@var BYBR factor for converting cloud particle radius to effect. radius
!@var XMASS dummy variable

C**** Set-able variables
!@dbparam LMCM max level for originating MC plumes
      INTEGER :: LMCM = -1 ! defaults to LS1-1 if not set in rundeck
!@dbparam ISC integer to turn on computation of stratocumulus clouds
      INTEGER :: ISC = 0  ! set ISC=1 to compute stratocumulus clouds
!@dbparam U00wtrX multiplies U00ice for critical humidity for water clds
      REAL*8 :: U00wtrX = 1.0d0     ! default
!@dbparam U00ice critical humidity for ice cloud condensation
      REAL*8 :: U00ice = .7d0       ! default
!@dbparam HRMAX maximum distance an air parcel rises from surface
      REAL*8 :: HRMAX = 1000.d0     ! default (m)
!@dbparam RWCldOX multiplies part.size of water clouds over ocean
      REAL*8 :: RWCldOX=1.d0
!@dbparam RICldX multiplies part.size of ice clouds at 1000mb
!@+       RICldX changes linearly to 1 as p->0mb
      REAL*8 :: RICldX=1.d0 , xRICld

#ifdef TRACERS_ON
!@var ntx,NTIX: Number and Indices of active tracers used in convection
      integer, dimension(ntm) :: ntix
      integer ntx
#endif

C**** input variables
      LOGICAL DEBUG
!@var RA ratio of primary grid box to secondary gridbox
      REAL*8, DIMENSION(IM) :: RA
!@var UM,VM,U_0,V_0 velocity related variables (UM,VM)=(U,V)*AIRM
      REAL*8, DIMENSION(IM,LM) :: UM,VM
      REAL*8, DIMENSION(IM,LM) :: U_0,V_0

!@var Miscellaneous vertical arrays set in driver
!@var PLE pressure at layer edge
!@var LHP array of precip phase ! may differ from LHX
      REAL*8, DIMENSION(LM+1) :: PLE,LHP
      REAL*8, DIMENSION(LM) :: PL,PLK,AIRM,BYAM,ETAL,TL,QL,TH,RH,WMX
     *     ,VSUBL,MCFLX,DGDSM,DPHASE,DTOTW,DQCOND,DGDQM,AQ,DPDT,RH1
!@var PL layer pressure (mb)
!@var PLK PL**KAPA
!@var AIRM the layer's pressure depth (mb)
!@var BYAM 1./AIRM
!@var ETAL fractional entrainment rate
!@var TL, QL temperature, specific humidity of the layer
!@var TH potential temperature (K)
!@var RH relative humidity
!@var RH1 relative humidity to compare with the threshold humidity
!@var WMX cloud water mixing ratio (kg/Kg)
!@var VSUBL downward vertical velocity due to cumulus subsidence (cm/s)
!@var MCFLX, DGDSM, DPHASE, DQCOND, DGDQM dummy variables
!@var AQ time change rate of specific humidity (s**-1)
!@var DPDT time change rate of pressure (mb/s)
      REAL*8, DIMENSION(LM+1) :: PRECNVL
!@var PRECNVL convective precip entering the layer top
C**** new arrays must be set to model arrays in driver (before MSTCNV)
      REAL*8, DIMENSION(LM) :: SDL,WML
!@var SDL vertical velocity in sigma coordinate
!@var WML cloud water mixing ratio (kg/Kg)
C**** new arrays must be set to model arrays in driver (after MSTCNV)
      REAL*8, DIMENSION(LM) :: TAUMCL,SVLATL,CLDMCL,SVLHXL,SVWMXL
!@var TAUMCL convective cloud optical thickness
!@var SVLATL saved LHX for convective cloud
!@var CLDMCL convective cloud cover
!@var SVLHXL saved LHX for large-scale cloud
!@var SVWMXL saved detrained convective cloud water
      REAL*8, DIMENSION(LM) :: CSIZEL
!@var CSIZEL cloud particle radius (micron)
C**** new arrays must be set to model arrays in driver (before LSCOND)
      REAL*8, DIMENSION(LM) :: TTOLDL,CLDSAVL
!@var TTOLDL previous potential temperature
!@var CLDSAVL saved large-scale cloud cover
C**** new arrays must be set to model arrays in driver (after LSCOND)
      REAL*8, DIMENSION(LM) :: SSHR,DCTEI,TAUSSL,CLDSSL
!@var SSHR,DCTEI height diagnostics of dry and latent heating by MC
!@var TAUSSL large-scale cloud optical thickness
!@var CLDSSL large-scale cloud cover

!@var SM,QM Vertical profiles of T/Q
      REAL*8, DIMENSION(LM) :: SM,QM
      REAL*8, DIMENSION(NMOM,LM) :: SMOM,QMOM

#ifdef TRACERS_ON
!@var TM Vertical profiles of tracers
      REAL*8, DIMENSION(LM,NTM) :: TM
      REAL*8, DIMENSION(nmom,lm,ntm) :: TMOM
      COMMON/CLD_TRCCOM/TM,TMOM
!$OMP  THREADPRIVATE (/CLD_TRCCOM/)
#ifdef TRACERS_WATER
!@var TRWML Vertical profile of liquid water tracers (kg)
!@var TRSVWML New liquid water tracers from m.c. (kg)
      REAL*8, DIMENSION(NTM,LM) :: TRWML, TRSVWML
!@var TRPRSS super-saturated tracer precip (kg)
!@var TRPRMC moist convective tracer precip (kg)
      REAL*8, DIMENSION(NTM)    :: TRPRSS,TRPRMC
!@var FQCONDT fraction of tracer that condenses
!@var FQEVPT  fraction of tracer that evaporates (in downdrafts)
!@var FPRCPT fraction of tracer that evaporates (in net re-evaporation)
!@var FWASHT  fraction of tracer scavenged by below-cloud precipitation
      REAL*8 :: FQCONDT, FWASHT, FPRCPT, FQEVPT
!@var WMXTR available water mixing ratio for tracer condensation ( )?
!@var b_beta_DT precipitating gridbox fraction from lowest precipitating
!@+   layer. The name was chosen to correspond to Koch et al. p. 23,802.
!@var precip_mm precipitation (mm) from the grid box above for washout
      REAL*8 WMXTR, b_beta_DT, precip_mm
      COMMON/CLD_WTRTRCCOM/TRWML, TRSVWML,TRPRSS,TRPRMC
     *  ,FQCONDT, FWASHT, FPRCPT, FQEVPT,WMXTR, b_beta_DT, precip_mm
!$OMP  THREADPRIVATE (/CLD_WTRTRCCOM/)
#endif
#endif

!@var KMAX index for surrounding velocity
!@var LP50 50mb level
      INTEGER ::  KMAX,LP50
!@var PEARTH fraction of land in grid box
!@var TS average surface temperture (C)
!@var RIS, RI1, RI2 Richardson numbers
      REAL*8 :: PEARTH,TS,QS,US,VS,RIS,RI1,RI2,DXYPJ
!@var DCL max level of planetary boundary layer
      INTEGER :: DCL

C**** output variables
      REAL*8 :: PRCPMC,PRCPSS,HCNDSS,WMSUM
!@var PRCPMC precip due to moist convection
!@var PRCPSS precip due to large-scale condensation
!@var HCNDSS heating due to large-scale condensation
!@var WMSUM cloud liquid water path
      REAL*8 :: CLDSLWIJ,CLDDEPIJ
!@var CLDSLWIJ shallow convective cloud cover
!@var CLDDEPIJ deep convective cloud cover
      INTEGER :: LMCMAX,LMCMIN
!@var LMCMAX upper-most convective layer
!@var LMCMIN lowerest convective layer
!@var AIRXL is convective mass flux (kg/m*m)
      REAL*8 AIRXL
!@var RNDSSL stored random number sequences
      REAL*8  RNDSSL(3,LM)
CCOMP  does not work yet:
CCOMP  THREADPRIVATE (RA,UM,VM,U_0,V_0,PLE,PL,PLK,AIRM,BYAM,ETAL
CCOMP*  ,TL,QL,TH,RH,WMX,VSUBL,MCFLX,SSHR,DGDSM,DPHASE
CCOMP*  ,DTOTW,DQCOND,DCTEI,DGDQM,dxypj
CCOMP*  ,AQ,DPDT,PRECNVL,SDL,WML,SVLATL,SVLHXL,SVWMXL,CSIZEL,RH1
CCOMP*  ,TTOLDL,CLDSAVL,TAUMCL,CLDMCL,TAUSSL,CLDSSL,RNDSSL
CCOMP*  ,SM,QM,SMOM,QMOM,PEARTH,TS,QS,US,VS,DCL,RIS,RI1,RI2, AIRXL
CCOMP*  ,PRCPMC,PRCPSS,HCNDSS,WMSUM,CLDSLWIJ,CLDDEPIJ,LMCMAX
CCOMP*  ,LMCMIN,KMAX,DEBUG)
      COMMON/CLDPRV/RA,UM,VM,U_0,V_0,PLE,PL,PLK,AIRM,BYAM,ETAL
     *  ,TL,QL,TH,RH,WMX,VSUBL,MCFLX,SSHR,DGDSM,DPHASE,LHP
     *  ,DTOTW,DQCOND,DCTEI,DGDQM,DXYPJ
     *  ,AQ,DPDT,PRECNVL,SDL,WML,SVLATL,SVLHXL,SVWMXL,CSIZEL,RH1
     *  ,TTOLDL,CLDSAVL,TAUMCL,CLDMCL,TAUSSL,CLDSSL,RNDSSL
     *  ,SM,QM,SMOM,QMOM,PEARTH,TS,QS,US,VS,RIS,RI1,RI2, AIRXL
     *  ,PRCPMC,PRCPSS,HCNDSS,WMSUM,CLDSLWIJ,CLDDEPIJ,DEBUG
     *  ,LMCMAX,LMCMIN,KMAX,DCL     ! integers last (alignment)
!$OMP  THREADPRIVATE (/CLDPRV/)

      CONTAINS

      SUBROUTINE MSTCNV(IERR,LERR)
!@sum  MSTCNV moist convective processes (precip, convective clouds,...)
!@auth M.S.Yao/A. Del Genio (modularisation by Gavin Schmidt)
!@ver  1.0 (taken from CB265)
!@calls adv1d,QSAT,DQSATDT,THBAR
      IMPLICIT NONE
      REAL*8 LHX,MPLUME,MCLOUD,MPMAX,SENV,QENV
!@var LHX latent heat of evaporation (J/Kg)
!@var MPLUME mass of convective plume (mb)
!@var MCLOUD air mass available for re-evaporation of precip
!@var MPMAX convective plume at the detrainment level
!@var SENV,QENV dummy variables

C**** functions
      REAL*8 :: QSAT, DQSATDT
!@var QSAT saturation humidity
!@var DQSATDT dQSAT/dT

      REAL*8, DIMENSION(0:LM) :: CM     !@var CM air mass of subsidence
      REAL*8, DIMENSION(IM) :: UMP,VMP,UMDN,VMDN
!@var UMP, VMP momentum carried by convective plumes
!@var UMDN,VMDN dummy variables
!@var DQM,DSM,DQMR,DSMR Vertical profiles of T/Q and changes
      REAL*8, DIMENSION(LM) ::
     * SMOLD,QMOLD, DQM,DSM,DQMR,DSMR
!@var SMOLD,QMOLD profiles prior to any moist convection
      REAL*8, DIMENSION(LM) :: F,CMNEG
      REAL*8, DIMENSION(NMOM,LM) :: FMOM
      REAL*8, DIMENSION(NMOM,LM) :: SMOMOLD,QMOMOLD
      REAL*8, DIMENSION(NMOM,LM) :: DSMOM,DQMOM,DSMOMR,DQMOMR
      REAL*8, DIMENSION(NMOM) ::
     &     SMOMP,QMOMP, SMOMPMAX,QMOMPMAX, SMOMDN,QMOMDN

#ifdef TRACERS_ON
!@var TMOLD: old TM (tracer mass)
!@var DTM,DTMR: Vertical profiles of Tracers mass and changes
      REAL*8, DIMENSION(LM,NTM)      :: TMOLD,  DTM,  DTMR,TM1
      REAL*8, DIMENSION(NMOM,LM,NTM) :: TMOMOLD,DTMOM,DTMOMR
      REAL*8, DIMENSION(NTM)      :: TMP,  TMPMAX,  TMDN, TENV
      REAL*8, DIMENSION(NMOM,NTM) :: TMOMP,TMOMPMAX,TMOMDN
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(NTM)      :: TRPRCP
      REAL*8, DIMENSION(NTM,LM)   :: TRCOND
      REAL*8 DTSUM,HEFF
#endif
#endif

      REAL*8, DIMENSION(LM) ::
     *     DM,COND,CDHEAT,CCM,SM1,QM1,DMR,ML,SMT,QMT,TPSAV,SVTP,DDM
!@var DM change in air mass
!@var COND condensate
!@var CDHEAT heating due to condensation
!@var CCM convective plume mass
!@var SM1, QM1 dummy variables
!@var DMR change of air mass
!@var ML layer air mass
!@var SMT, QMT dummy variables
!@var TPSAV, SVTP  arrays to save plume temperature
!@var DDM downdraft mass

!@var IERRT,LERRT error reports from advection
      INTEGER :: IERRT,LERRT
      INTEGER LDRAFT,LMAX,LMIN,MCCONT,MAXLVL
     *     ,MINLVL,ITER,IC,LFRZ,NSUB,LDMIN
!@var LDRAFT the layer the downdraft orginates
!@var LEVAP the layer evaporation of precip starts
!@var LMAX, LMIN the lowest, the highest layer of a convective event
!@var MCCONT interger to count convective events
!@var MAXLVL, MINLVL the lowest, the highest layer of convective events
!@var ITER number for iteration
!@var IC interger for cloud types
!@var LFRZ freezing level
!@var nsub = LMAX - LMIN + 1
!@var LDMIN the lowest layer the downdraft drops
      REAL*8 TERM1,FMP0,SMO1
     *     ,QMO1,SMO2,QMO2,SDN,QDN,SUP,QUP,SEDGE,QEDGE,WMDN,WMUP,SVDN
     *     ,SVUP,WMEDG,SVEDG,DMSE,FPLUME,DFP,FMP2,FRAT1,FRAT2,SMN1
     *     ,QMN1,SMN2,QMN2,SMP,QMP,TP,GAMA,DQSUM,TNX,TNX1
     *     ,DQ,DMSE1,FCTYPE,BETAU,ALPHAU
     *     ,CDHDRT,DDRAFT,DELTA
     *     ,ALPHA,BETA,CDHM,CDHSUM,CLDM,CLDREF,CONSUM,DQEVP
     *     ,DQRAT,EPLUME,ETADN,ETAL1,EVPSUM,FCDH
     *     ,FCDH1,FCLD,FCLOUD,FDDL,FDDP,FENTR,FENTRA,FEVAP,FLEFT
     *     ,FQCOND,FQEVP,FPRCP,FSEVP,FSSUM,HEAT1
     *     ,PRHEAT,PRCP
     *     ,QMDN,QMIX,QMPMAX,QMPT,QNX,QSATC,QSATMP
     *     ,RCLD,RCLDE,SLH,SMDN,SMIX,SMPMAX,SMPT,SUMAJ
     *     ,SUMDP,DDRUP,EDRAFT
     *     ,TOLD,TOLD1,TEMWM,TEM,WTEM,WCONST,WORK,DELTX1
!@var TERM1 contribution to non-entraining convective cloud
!@var FMP0 non-entraining convective mass
!@var SMO1,QMO1,SMO2,QMO2,SDN,QDN,SUP,QUP,SEDGE,QEDGE dummy variables
!@var WMDN,WMUP,SVDN,SVUP,WMEDG,SVEDG,DDRUP dummy variables
!@var DMSE difference in moist static energy
!@var FPLUME fraction of convective plume
!@var DFP an iterative increment
!@var FMP2,FRAT1,FRAT2,SMN1,QMN1,SMN2,QMN2 dummy variables
!@var SMP, QMP plume's SM, QM
!@var TP plume's temperature (K)
!@var GAMA,DQSUM,TNX,DQ,CONSUM,BETAU,ALPHAU dummy variables
!@var DMSE1 difference in moist static energy
!@var FCTYPE fraction for convective cloud types
!@var CDHDRT,ALPHA,BETA,CDHM,CDHSUM,CLDREF dummy variables
!@var DDRAFT downdraft mass
!@var DELTA fraction of plume stays in the layer
!@var CLDM subsidence due to convection
!@var DQEVP amount of condensate that evaporates in downdrafts
!@var DQRAT fraction for the condensate to evaporate
!@var EPLUME air mass of entrainment
!@var ETADN fraction of the downdraft
!@var ETAL1 fractional entrainment rate
!@var EVPSUM,FCDH,FCDH1,FCLD,FCLOUD,FDDL,FDDP dummy variables
!@var FENTR fraction of entrainment
!@var FENTRA fraction of entrainment
!@var FEVAP fraction of air available for precip evaporation
!@var FLEFT fraction of plume after removing downdraft mass
!@var FQCOND fraction of water vapour that condenses in plume
!@var FQEVP fraction of water vapour that evaporates in downdraft
!@var FPRCP fraction of evaporated precipitation
!@var FSEVP fraction of energy lost to evaporate
!@var FSSUM fraction of energy lost to evaporate
!@var HEAT1 heating due to phase change
!@var PRHEAT energy of condensate
!@var PRCP precipipation
!@var QMDN,QMIX,SMDN,SMIX dummy variables
!@var QMPMAX,SMPMAX detrainment of QMP, SMP
!@var QMPT,SMPT dummy variables
!@var QNX,SUMAJ,SUMDP dummy variables
!@var QSATC saturation vapor mixing ratio
!@var QSATMP plume's saturation vapor mixing ratio
!@var RCLD,RCLDE cloud particle's radius, effective radius
!@var SLH LHX/SHA
!@var EDRAFT entrainment into downdrafts
!@var TOLD,TOLD1 old temperatures
!@var TEMWM,TEM,WTEM,WCONST dummy variables
!@var WORK work done on convective plume

      LOGICAL MC1  !@var MC1 true for the first convective event

      REAL*8,  PARAMETER :: CK1 = 1.       !@param CK1 a tunning const.

      INTEGER K,L,N  !@var K,L,N loop variables
      INTEGER ITYPE  !@var convective cloud types
!@var IERR,LERR error reports from advection
      INTEGER, INTENT(OUT) :: IERR,LERR
!@var DUM, DVM changes of UM,VM
      REAL*8, DIMENSION(IM,LM) :: DUM,DVM

      REAL*8 THBAR  !@var THBAR virtual temperature at layer edge
!@var BELOW_CLOUD logical- is the current level below cloud?
      LOGICAL BELOW_CLOUD
C****
C**** MOIST CONVECTION
C****
C**** CONVECTION OCCURS AT THE LOWEST MOIST CONVECTIVELY UNSTABLE
C**** LEVEL AND CONTINUES UNTIL A STABLE LAYER PAIR IS REACHED.  RE-
C**** EVAPORATION AND PRECIPITATION ARE COMPUTED AT THE END OF THIS
C**** CYCLE.  THE WHOLE PROCESS MAY BE REPEATED FROM A NEW LOWEST
C**** UNSTABLE LEVEL.
C****
      ierr=0 ; lerr=0
      LMCMIN=0
      LMCMAX=0
      MCCONT=0
C**** initiallise arrays of computed output
      TAUMCL=0
      SVWMXL=0
      SVLATL=0
      VSUBL=0
      PRECNVL=0
      CLDMCL=0
      CLDSLWIJ=0
      CLDDEPIJ=0
      PRCPMC=0.
      SVTP=0
      TPSAV=0
      CSIZEL=10.    !  effective droplet radius in stem of convection
#ifdef TRACERS_WATER
      trsvwml = 0.
      TRPRCP = 0.
      TRPRMC = 0.
#endif
C**** zero out diagnostics
         MCFLX =0.
         DGDSM=0.
         DPHASE=0.
         DTOTW=0.
         DQCOND=0.
         DGDQM=0.
C**** save initial values (which will be updated after subsid)
      SM1=SM
      QM1=QM
#ifdef TRACERS_ON
      TM1(:,1:NTX) = TM(:,1:NTX)
#endif
C**** OUTER MC LOOP OVER BASE LAYER
      DO 600 LMIN=1,LMCM-1
      MAXLVL=0
      MINLVL=LM
C****
C**** COMPUTE THE CONVECTIVE MASS OF THE NON-ENTRAINING PART
C****
      TERM1=-10.*CK1*SDL(LMIN+1)*BYGRAV
      FMP0=TERM1*XMASS
      IF(FMP0.LE.0.) FMP0=0.
C**** CREATE A PLUME IN THE BOTTOM LAYER
C****
C**** ITERATION TO FIND FPLUME WHICH RESTORES THE ATM TO NEUTRAL STATE
C****
      DELTX1=.608d0
      SMO1=SM(LMIN)
      QMO1=QM(LMIN)
      SMO2=SM(LMIN+1)
      QMO2=QM(LMIN+1)
      SDN=SMO1*BYAM(LMIN)
      SUP=SMO2*BYAM(LMIN+1)
      SEDGE=THBAR(SUP,SDN)
      QDN=QMO1*BYAM(LMIN)
      QUP=QMO2*BYAM(LMIN+1)
      WMDN=WML(LMIN)
      WMUP=WML(LMIN+1)
      SVDN=SDN*(1.+DELTX1*QDN-WMDN)
      SVUP=SUP*(1.+DELTX1*QUP-WMUP)
      QEDGE=.5*(QUP+QDN)
      WMEDG=.5*(WMUP+WMDN)
      SVEDG=SEDGE*(1.+DELTX1*QEDGE-WMEDG)
      LHX=LHE
      SLH=LHX*BYSHA
      DMSE=(SVUP-SVEDG)*PLK(LMIN+1)+(SVEDG-SVDN)*PLK(LMIN)+
     *  SLHE*(QSAT(SUP*PLK(LMIN+1),LHX,PL(LMIN+1))-QDN)
      IF(DMSE.GT.-1d-10) GO TO 600
C****
      FPLUME=.25
      DFP = .25
      DO ITER=1,8
      DFP=DFP*0.5
      FMP2=FPLUME*AIRM(LMIN)
      FRAT1=FMP2*BYAM(LMIN+1)
      FRAT2=FMP2*BYAM(LMIN+2)
      SMN1=SMO1*(1.-FPLUME)+FRAT1*SMO2
      QMN1=QMO1*(1.-FPLUME)+FRAT1*QMO2
      SMN2=SMO2*(1.-FRAT1)+FRAT2*SM(LMIN+2)
      QMN2=QMO2*(1.-FRAT1)+FRAT2*QM(LMIN+2)
      SMP=SMO1*FPLUME
      QMP=QMO1*FPLUME
      TP=SMO1*PLK(LMIN+1)*BYAM(LMIN)
      QSATMP=FMP2*QSAT(TP,LHX,PL(LMIN+1))
      GAMA=SLH*QSATMP*DQSATDT(TP,LHX)/FMP2
      DQSUM=(QMP-QSATMP)/(1.+GAMA)
      IF(DQSUM.LE.0.) GO TO 205
      FEVAP=.5*FPLUME
      MCLOUD=FEVAP*AIRM(LMIN+1)
      TNX=SMO2*PLK(LMIN+1)*BYAM(LMIN+1)
      QNX=QMO2*BYAM(LMIN+1)
      QSATC=QSAT(TNX,LHX,PL(LMIN+1))
      DQ=MCLOUD*(QSATC-QNX)/(1.+SLH*QSATC*DQSATDT(TNX,LHX))
      IF(DQ.GT.DQSUM) DQ=DQSUM
      SMN2=SMN2-SLH*DQ/PLK(LMIN+1)
      QMN2=QMN2+DQ
      DQSUM=DQSUM-DQ
      IF(DQSUM.LE.0.) GO TO 205
      MCLOUD=FEVAP*AIRM(LMIN)
      TNX=SMO1*PLK(LMIN)*BYAM(LMIN)
      QNX=QMO1*BYAM(LMIN)
      QSATC=QSAT(TNX,LHX,PL(LMIN))
      DQ=MCLOUD*(QSATC-QNX)/(1.+SLH*QSATC*DQSATDT(TNX,LHX)) ! new func'n
      IF(DQ.GT.DQSUM) DQ=DQSUM
      SMN1=SMN1-SLH*DQ/PLK(LMIN)
      QMN1=QMN1+DQ
  205 SDN=SMN1*BYAM(LMIN)
      SUP=SMN2*BYAM(LMIN+1)
      SEDGE=THBAR(SUP,SDN)
      QDN=QMN1*BYAM(LMIN)
      QUP=QMN2*BYAM(LMIN+1)
      SVDN=SDN*(1.+DELTX1*QDN-WMDN)
      SVUP=SUP*(1.+DELTX1*QUP-WMUP)
      QEDGE=.5*(QUP+QDN)
      SVEDG=SEDGE*(1.+DELTX1*QEDGE-WMEDG)
      DMSE1=(SVUP-SVEDG)*PLK(LMIN+1)+(SVEDG-SVDN)*PLK(LMIN)+
     *  SLHE*(QSAT(SUP*PLK(LMIN+1),LHX,PL(LMIN+1))-QDN)
      IF (ABS(DMSE1).LE.1.d-3) GO TO 411
      IF(DMSE1.GT.1.d-3) FPLUME=FPLUME-DFP
      IF(DMSE1.LT.-1.d-3) FPLUME=FPLUME+DFP
      END DO
  411 IF(FPLUME.LE..001) GO TO 600
C****
C**** ITERATION THROUGH CLOUD TYPES
C****
      ITYPE=2                        ! always 2 types of clouds:
C     IF(LMIN.LE.2) ITYPE=2          ! entraining & non-entraining
      FCTYPE=1.
C**** SAVE ORIG PROFILES
      SMOLD(:) = SM(:)
      SMOMOLD(:,:) = SMOM(:,:)
      QMOLD(:) = QM(:)
      QMOMOLD(:,:) = QMOM(:,:)
#ifdef TRACERS_ON
      TMOLD(:,1:NTX) = TM(:,1:NTX)
      TMOMOLD(:,:,1:NTX) = TMOM(:,:,1:NTX)
#endif

      DO 570 IC=1,ITYPE
C**** INITIALLISE VARIABLES USED FOR EACH TYPE
      DO L=1,LM
        COND(L)=0.
        CDHEAT(L)=0.
        DM(L)=0.
        DMR(L)=0.
        CCM(L)=0.
        DDM(L)=0.
      END DO
      DUM(1:KMAX,:)=0.
      DVM(1:KMAX,:)=0.
      DSM(:) = 0.
      DSMOM(:,:) = 0.
      DSMR(:) = 0.
      DSMOMR(:,:) = 0.
      DQM(:) = 0.
      DQMOM(:,:) = 0.
      DQMR(:) = 0.
      DQMOMR(:,:) = 0.
#ifdef TRACERS_ON
      DTM(:,1:NTX) = 0.
      DTMOM(:,:,1:NTX) = 0.
      DTMR(:,1:NTX) = 0.
      DTMOMR(:,:,1:NTX) = 0.
#endif
#ifdef TRACERS_WATER
      TRCOND = 0.
#endif
      MC1=.FALSE.
      LHX=LHE
      MPLUME=MIN(1.*AIRM(LMIN),1.*AIRM(LMIN+1))
      IF(MPLUME.GT.FMP2) MPLUME=FMP2
      IF(ITYPE.EQ.2) THEN
      FCTYPE=1.
      IF(MPLUME.GT.FMP0) FCTYPE=FMP0/MPLUME
      IF(IC.EQ.2) FCTYPE=1.-FCTYPE
      IF(FCTYPE.LT.0.001) GO TO 570
      END IF
      MPLUME=MPLUME*FCTYPE
C     FPLUM0=FMP1*BYAM(LMIN)
      FPLUME=MPLUME*BYAM(LMIN)
      SMP  =  SMOLD(LMIN)*FPLUME
      SMOMP(xymoms)=SMOMOLD(xymoms,LMIN)*FPLUME
      QMP  =  QMOLD(LMIN)*FPLUME
      QMOMP(xymoms)=QMOMOLD(xymoms,LMIN)*FPLUME
      TPSAV(LMIN)=SMP*PLK(LMIN)/MPLUME
      DMR(LMIN)=-MPLUME
        DSMR(LMIN)=-SMP
      DSMOMR(xymoms,LMIN)=-SMOMP(xymoms)
      DSMOMR(zmoms,LMIN)=-SMOMOLD(zmoms,LMIN)*FPLUME
        DQMR(LMIN)=-QMP
      DQMOMR(xymoms,LMIN)=-QMOMP(xymoms)
      DQMOMR(zmoms,LMIN)=-QMOMOLD(zmoms,LMIN)*FPLUME
#ifdef TRACERS_ON
      TMP(1:NTX) = TMOLD(LMIN,1:NTX)*FPLUME
      TMOMP(xymoms,1:NTX)=TMOMOLD(xymoms,LMIN,1:NTX)*FPLUME
        DTMR(LMIN,1:NTX)=-TMP(1:NTX)
      DTMOMR(xymoms,LMIN,1:NTX)=-TMOMP(xymoms,1:NTX)
      DTMOMR( zmoms,LMIN,1:NTX)=-TMOMOLD(zmoms,LMIN,1:NTX)*FPLUME
#endif
      DO K=1,KMAX
         UMP(K)=UM(K,LMIN)*FPLUME
         DUM(K,LMIN)=-UMP(K)
         VMP(K)=VM(K,LMIN)*FPLUME
         DVM(K,LMIN)=-VMP(K)
      ENDDO
C****
C**** RAISE THE PLUME TO THE TOP OF CONVECTION AND CALCULATE
C**** ENTRAINMENT, CONDENSATION, AND SECONDARY MIXING
C****
      CDHSUM=0.
      CDHDRT=0.
      ETADN=0.
      LDRAFT=LM
      EVPSUM=0.
      DDRAFT=0.
      LFRZ=0
      LMAX=LMIN
 220  L=LMAX+1
      SVTP(L)=SMP*PLK(L)/(MPLUME+teeny)
C**** TEST FOR SUFFICIENT AIR, MOIST STATIC STABILITY AND ENERGY
C     IF(L.GT.LMIN+1.AND.SDL(L).GT.0.) GO TO 340
      IF(MPLUME.LE..001*AIRM(L)) GO TO 340
      SDN=SMP/MPLUME
      SUP=SM1(L)*BYAM(L)
      QDN=QMP/MPLUME
      QUP=QM1(L)*BYAM(L)
      WMDN=0.
      WMUP=WML(L)
      SVDN=SDN*(1.+DELTX1*QDN-WMDN)
      SVUP=SUP*(1.+DELTX1*QUP-WMUP)
      IF(LMAX.GT.LMIN) THEN
      SEDGE=THBAR(SUP,SDN)
      QEDGE=.5*(QUP+QDN)
      WMEDG=.5*(WMUP+WMDN)
      SVEDG=SEDGE*(1.+DELTX1*QEDGE-WMEDG)
      LHX=LHE
      DMSE=(SVUP-SVEDG)*PLK(L)+(SVEDG-SVDN)*PLK(L-1)+
     *  SLHE*(QSAT(SUP*PLK(L),LHX,PL(L))-QDN)
      IF(DMSE.GT.-1d-10) GO TO 340
      END IF
      IF(PLK(L-1)*(SVUP-SVDN)+SLHE*(QUP-QDN).GE.0.) GO TO 340
C**** TEST FOR CONDENSATION ALSO DETERMINES IF PLUME REACHES UPPER LAYER
      TP=SMP*PLK(L)/MPLUME
      TPSAV(L)=TP
      IF(TPSAV(L-1).GE.TF.AND.TPSAV(L).LT.TF) LFRZ=L-1
      IF(TP.LT.TI) LHX=LHS
      QSATMP=MPLUME*QSAT(TP,LHX,PL(L))
      IF(QMP.LT.QSATMP) GO TO 340
      IF(TP.GE.TF.OR.LHX.EQ.LHS) GO TO 290
      LHX=LHS
      QSATMP=MPLUME*QSAT(TP,LHX,PL(L))
  290 SLH=LHX*BYSHA
C****
C**** DEPOSIT PART OF THE PLUME IN LOWER LAYER
C****
      IF(MPLUME.GT..95*AIRM(L)) THEN
        DELTA=(MPLUME-.95*AIRM(L))/MPLUME
        DM(L-1)=DM(L-1)+DELTA*MPLUME
        MPLUME=.95*AIRM(L)

        DSM(L-1)=  DSM(L-1)+DELTA*SMP
        SMP = SMP  *(1.-DELTA)
        DSMOM(xymoms,L-1)=DSMOM(xymoms,L-1)+DELTA*SMOMP(xymoms)
        SMOMP(xymoms) = SMOMP(xymoms)*(1.-DELTA)

        DQM(L-1)=  DQM(L-1)+DELTA*QMP
        QMP = QMP  *(1.-DELTA)
        DQMOM(xymoms,L-1)=DQMOM(xymoms,L-1)+DELTA*QMOMP(xymoms)
        QMOMP(xymoms) = QMOMP(xymoms)*(1.-DELTA)

        DO K=1,KMAX
          DUM(K,L-1)=DUM(K,L-1)+UMP(K)*DELTA
          DVM(K,L-1)=DVM(K,L-1)+VMP(K)*DELTA
          UMP(K)=UMP(K)-UMP(K)*DELTA
          VMP(K)=VMP(K)-VMP(K)*DELTA
        ENDDO

#ifdef TRACERS_ON
        DTM(L-1,1:NTX) = DTM(L-1,1:NTX)+DELTA*TMP(1:NTX)
        DTMOM(xymoms,L-1,1:NTX)=DTMOM(xymoms,L-1,1:NTX)+DELTA
     *       *TMOMP(xymoms,1:NTX)
        TMP(1:NTX) = TMP(1:NTX)*(1.-DELTA)
        TMOMP(xymoms,1:NTX) = TMOMP(xymoms,1:NTX)*(1.-DELTA)
#endif
      END IF
C****
C**** CONVECTION IN UPPER LAYER   (WORK DONE COOLS THE PLUME)
C****
      WORK=MPLUME*(SUP-SDN)*(PLK(L-1)-PLK(L))/PLK(L-1)
C     SMP=SMP-WORK
      DSM(L-1)=DSM(L-1)-WORK
      CCM(L-1)=MPLUME
C****
C**** ENTRAINMENT
C****
      IF(IC.EQ.2.OR.(IC.EQ.1.AND.PL(L).GE.800.)) THEN
      FENTR=ETAL(L)*FPLUME
      IF(FENTR+FPLUME.GT.1.) FENTR=1.-FPLUME
      ETAL1=FENTR/(FPLUME+teeny)
      FPLUME=FPLUME+FENTR
      EPLUME=MPLUME*ETAL1
C**** Reduce EPLUME so that mass flux is less than mass in box
      IF (EPLUME.GT.AIRM(L)*0.975d0-MPLUME) THEN
        EPLUME=AIRM(L)*0.975d0-MPLUME
      END IF
      MPLUME=MPLUME+EPLUME
      FENTRA = EPLUME*BYAM(L)
      DSMR(L)=DSMR(L)-EPLUME*SUP        ! = DSM(L)-SM(L)*FENTRA
      DSMOMR(:,L)=DSMOMR(:,L)-SMOM(:,L)*FENTRA
      DQMR(L)=DQMR(L)-EPLUME*QUP        ! = DQM(L)-QM(L)*FENTRA
      DQMOMR(:,L)=DQMOMR(:,L)-QMOM(:,L)*FENTRA
      DMR(L)=DMR(L)-EPLUME
      SMP=SMP+EPLUME*SUP
      SMOMP(xymoms)= SMOMP(xymoms)+ SMOM(xymoms,L)*FENTRA
      QMP=QMP+EPLUME*QUP
      QMOMP(xymoms)= QMOMP(xymoms)+ QMOM(xymoms,L)*FENTRA
#ifdef TRACERS_ON
      DTMR(L,1:NTX) = DTMR(L,1:NTX)-TM(L,1:NTX)*FENTRA
      DTMOMR(:,L,1:NTX) = DTMOMR(:,L,1:NTX)-TMOM(:,L,1:NTX)*FENTRA
      TMP(1:NTX) = TMP(1:NTX)+TM(L,1:NTX)*FENTRA
      TMOMP(xymoms,1:NTX) = TMOMP(xymoms,1:NTX)+TMOM(xymoms,L,1:NTX)
     *     *FENTRA
#endif
      DO K=1,KMAX
         UMP(K)=UMP(K)+U_0(K,L)*EPLUME
         DUM(K,L)=DUM(K,L)-U_0(K,L)*EPLUME
         VMP(K)=VMP(K)+V_0(K,L)*EPLUME
         DVM(K,L)=DVM(K,L)-V_0(K,L)*EPLUME
      ENDDO
      END IF
C****
C**** CHECK THE DOWNDRAFT POSSIBILITY
C****
      IF(L-LMIN.LE.1) GO TO 291
      IF(ETADN.GT.1d-10) GO TO 291
      SMIX=.5*(SUP+SMP/MPLUME)
      QMIX=.5*(QUP+QMP/MPLUME)
C     WMIX=.5*(WMUP+WMDN)
C     SVMIX=SMIX*(1.+DELTX1*QMIX-WMIX)
C     DMMIX=(SVUP-SVMIX)*PLK(L)+
C    *  SLHE*(QSAT(SUP*PLK(L),LHX,PL(L))-QMIX)
C     IF(DMMIX.LT.1d-10) GO TO 291
      IF(SMIX.GE.SUP) GO TO 291
      IF(PL(L).GT.700.) GO TO 291
      LDRAFT=L
      ETADN=BY3
      FLEFT=1.-.5*ETADN
      DDRAFT=ETADN*MPLUME
      FDDP = .5*DDRAFT/MPLUME
      FDDL = .5*DDRAFT*BYAM(L)
      MPLUME=FLEFT*MPLUME
      SMDN=DDRAFT*SMIX         ! = SM(L)*FDDL +  SMP*FDDP
      SMOMDN(xymoms)=SMOM(xymoms,L)*FDDL +  SMOMP(xymoms)*FDDP
      SMP=FLEFT*SMP
      SMOMP(xymoms)= SMOMP(xymoms)*FLEFT
      QMDN=DDRAFT*QMIX         ! = QM(L)*FDDL +  QMP*FDDP
      QMOMDN(xymoms)=QMOM(xymoms,L)*FDDL +  QMOMP(xymoms)*FDDP
      QMP=FLEFT*QMP
      QMOMP(xymoms)= QMOMP(xymoms)*FLEFT
      DMR(L) = DMR(L)-.5*DDRAFT
      DSMR(L)=DSMR(L)-.5*DDRAFT*SUP        ! = DSM(L)-SM(L)*FDDL
      DSMOMR(:,L)=DSMOMR(:,L) - SMOM(:,L)*FDDL
      DQMR(L)=DQMR(L)-.5*DDRAFT*QUP        ! = DQM(L)-QM(L)*FDDL
      DQMOMR(:,L)=DQMOMR(:,L) - QMOM(:,L)*FDDL
#ifdef TRACERS_ON
      Tmdn(1:NTX) = tm(l,1:NTX)*fddl+Tmp(1:NTX)*fddp
      tmomdn(xymoms,1:NTX) = tmom(xymoms,l,1:NTX)*fddl+
     *     tmomp(xymoms,  1:NTX)*fddp
      dtmr    (l,1:NTX) = dtmr    (l,1:NTX)-fddl *tm    (l,1:NTX)
      dtmomr(:,l,1:NTX) = dtmomr(:,l,1:NTX)-fddl *tmom(:,l,1:NTX)
      Tmp         (1:NTX) = Tmp         (1:NTX)*fleft
      tmomp(xymoms,1:NTX) = tmomp(xymoms,1:NTX)*fleft
#endif
      DO K=1,KMAX
         UMDN(K)=.5*(ETADN*UMP(K)+DDRAFT*U_0(K,L))
         UMP(K)=UMP(K)*FLEFT
         DUM(K,L)=DUM(K,L)-.5*DDRAFT*U_0(K,L)
         VMDN(K)=.5*(ETADN*VMP(K)+DDRAFT*V_0(K,L))
         VMP(K)=VMP(K)*FLEFT
         DVM(K,L)=DVM(K,L)-.5*DDRAFT*V_0(K,L)
      ENDDO
C****
C**** CONDENSE VAPOR IN THE PLUME AND ADD LATENT HEAT
C****
  291 DQSUM=0.
      SMPT=SMP
      QMPT=QMP
      DO 292 N=1,3
      TP=SMP*PLK(L)/MPLUME
      QSATMP=MPLUME*QSAT(TP,LHX,PL(L))
      GAMA=SLH*QSATMP*DQSATDT(TP,LHX)/MPLUME
      DQ=(QMP-QSATMP)/(1.+GAMA)
      SMP=SMP+SLH*DQ/PLK(L)
      QMP=QMP-DQ
  292 DQSUM=DQSUM+DQ
      FQCOND = 0
      IF(DQSUM.GE.0.) THEN
        IF (QMPT.gt.teeny) FQCOND = DQSUM/QMPT
        QMOMP(xymoms) =  QMOMP(xymoms)*(1.-FQCOND)
      ELSE  ! no change
        DQSUM=0.
        SMP=SMPT
        QMP=QMPT
      END IF
      COND(L)=DQSUM
#ifdef TRACERS_WATER
C**** CONDENSING TRACERS
      WMXTR=DQSUM*BYAM(L)
      DO N=1,NTX
        CALL GET_COND_FACTOR(L,N,WMXTR,TP,LHX,FPLUME,FQCOND
     *       ,FQCONDT)
        TRCOND(N,L) = FQCONDT * TMP(N)
        TMP(N)         = TMP(N)         *(1.-FQCONDT)
        TMOMP(xymoms,N)= TMOMP(xymoms,N)*(1.-FQCONDT)
      END DO
#endif
      TAUMCL(L)=TAUMCL(L)+DQSUM
      CDHEAT(L)=SLH*COND(L)
      CDHSUM=CDHSUM+CDHEAT(L)
      IF(ETADN.GT.1d-10) CDHDRT=CDHDRT+SLH*COND(L)
C****
C**** UPDATE ALL QUANTITIES CARRIED BY THE PLUME
C****
      MCCONT=MCCONT+1
      IF(MCCONT.EQ.1) MC1=.TRUE.
      IF(MC1.AND.PLE(LMIN)-PLE(L+2).GE.450.) SVLATL(L)=LHX
      SMPMAX=SMP
      SMOMPMAX(xymoms) =  SMOMP(xymoms)
      QMPMAX=QMP
      QMOMPMAX(xymoms) =  QMOMP(xymoms)
#ifdef TRACERS_ON
C**** Tracers at top of plume
      TMPMAX(1:NTX) = TMP(1:NTX)
      TMOMPMAX(xymoms,1:NTX) = TMOMP(xymoms,1:NTX)
#endif
      MPMAX=MPLUME
      LMAX = LMAX + 1
      IF (LMAX.LT.LM) GO TO 220   ! CHECK FOR NEXT POSSIBLE LMAX
C**** UPDATE CHANGES CARRIED BY THE PLUME IN THE TOP CLOUD LAYER
  340 IF(LMIN.EQ.LMAX) GO TO 600
      IF(TPSAV(LMAX).GE.TF) LFRZ=LMAX
      DM(LMAX)=DM(LMAX)+MPMAX
      DSM(LMAX)=DSM(LMAX)+SMPMAX
      DSMOM(xymoms,LMAX)=DSMOM(xymoms,LMAX) + SMOMPMAX(xymoms)
      DQM(LMAX)=DQM(LMAX)+QMPMAX
      DQMOM(xymoms,LMAX)=DQMOM(xymoms,LMAX) + QMOMPMAX(xymoms)
#ifdef TRACERS_ON
C**** Add plume tracers at LMAX
      DTM(LMAX,1:NTX) = DTM(LMAX,1:NTX) + TMPMAX(1:NTX)
      DTMOM(xymoms,LMAX,1:NTX) =
     *   DTMOM(xymoms,LMAX,1:NTX) + TMOMPMAX(xymoms,1:NTX)
#endif
      CCM(LMAX)=0.
#ifdef TRACERS_SPECIAL_Shindell
      CALL calc_lightning(LMAX,LFRZ,IC)
#endif
      DO K=1,KMAX
         DUM(K,LMAX)=DUM(K,LMAX)+UMP(K)
         DVM(K,LMAX)=DVM(K,LMAX)+VMP(K)
      ENDDO
      CDHM=0.
      IF(MINLVL.GT.LMIN) MINLVL=LMIN
      IF(MAXLVL.LT.LMAX) MAXLVL=LMAX
      IF(LMCMIN.EQ.0) LMCMIN=LMIN
      IF(LMCMAX.LT.MAXLVL) LMCMAX=MAXLVL
C****
C**** PROCESS OF DOWNDRAFTING
C****
      LDMIN=LMIN
      EDRAFT=0.
      IF(ETADN.GT.1d-10) THEN
      CONSUM=0.
      DO 347 L=LMIN,LDRAFT
  347 CONSUM=CONSUM+COND(L)
      TNX=SMDN*PLK(LMIN)/DDRAFT
      QNX=QMDN/DDRAFT
      LHX=LHE
      IF(TPSAV(LMIN).LT.TF) LHX=LHS
      SLH=LHX*BYSHA
      QSATC=QSAT(TNX,LHX,PL(LMIN))
      DQ=(QSATC-QNX)/(1.+SLH*QSATC*DQSATDT(TNX,LHX))
      DQRAT=DQ*DDRAFT/(CONSUM+teeny)
      DO 346 L=LDRAFT,1,-1
      LHX=LHE
      IF (L.GE.LMIN) THEN  ! avoids reference to uninitiallised value
        IF(TPSAV(L).LT.TF) LHX=LHS
      END IF
      SLH=LHX*BYSHA
      DQEVP=DQRAT*COND(L)
      IF(DQEVP.GT.COND(L)) DQEVP=COND(L)
      IF (L.LT.LMIN) DQEVP=0.
      FSEVP = 0
      IF (ABS(PLK(L)*SMDN).gt.teeny) FSEVP = SLH*DQEVP/(PLK(L)*SMDN)
      SMDN=SMDN-SLH*DQEVP/PLK(L)
      SMOMDN(xymoms)=SMOMDN(xymoms)*(1.-FSEVP)
      FQEVP = 0
      IF (COND(L).gt.0.) FQEVP = DQEVP/COND(L)
      QMDN=QMDN+DQEVP
      COND(L)=COND(L)-DQEVP
      TAUMCL(L)=TAUMCL(L)-DQEVP
      CDHEAT(L)=CDHEAT(L)-DQEVP*SLH
      EVPSUM=EVPSUM+DQEVP*SLH
#ifdef TRACERS_WATER
C**** RE-EVAPORATION OF TRACERS IN DOWNDRAFTS
C**** (If 100% evaporation, allow all tracers to evaporate completely.)
      DO N=1,NTX
        IF(FQEVP.eq.1.) THEN                 ! total evaporation
          TMDN(N)     = TMDN(N) + TRCOND(N,L)
          TRCOND(N,L) = 0.D0
        ELSE ! otherwise, tracers evaporate dependent on type of tracer
          CALL GET_EVAP_FACTOR(N,TNX,LHX,.FALSE.,1d0,FQEVP,FQEVPT)
          TMDN(N)     = TMDN(N)     + FQEVPT * TRCOND(N,L)
          TRCOND(N,L) = TRCOND(N,L) - FQEVPT * TRCOND(N,L)
        END IF
      END DO
#endif
C**** ENTRAINMENT OF DOWNDRAFTS
cc    IF(L.LT.LDRAFT.AND.L.GT.1) THEN
cc      DDRUP=DDRAFT
cc      DDRAFT=DDRAFT*(1.+ETAL(L))
cc      IF(DDRUP.GT.DDRAFT) DDRAFT=DDRUP
cc      IF(DDRAFT.GT..95d0*(AIRM(L-1)+DMR(L-1)))
cc   *    DDRAFT=.95d0*(AIRM(L-1)+DMR(L-1))
cc      EDRAFT=DDRAFT-DDRUP
cc      IF (EDRAFT.gt.0) THEN  ! usual case, entrainment into downdraft
cc        FENTRA=EDRAFT*BYAM(L)
cc        SENV=SM1(L)/AIRM(L)
cc        QENV=QM1(L)/AIRM(L)
cc        SMDN=SMDN+EDRAFT*SENV
cc        QMDN=QMDN+EDRAFT*QENV
cc        SMOMDN(xymoms)= SMOMDN(xymoms)+ SMOM(xymoms,L)*FENTRA
cc        QMOMDN(xymoms)= QMOMDN(xymoms)+ QMOM(xymoms,L)*FENTRA
cc        DSMR(L)=DSMR(L)-EDRAFT*SENV
cc        DSMOMR(:,L)=DSMOMR(:,L)-SMOM(:,L)*FENTRA
cc        DQMR(L)=DQMR(L)-EDRAFT*QENV
cc        DQMOMR(:,L)=DQMOMR(:,L)-QMOM(:,L)*FENTRA
cc        DMR(L)=DMR(L)-EDRAFT
#ifdef TRACERS_ON
cc        Tenv(1:NTX)=tm1(l,1:NTX)/airm(l)
cc        TMDN(1:NTX)=TMDN(1:NTX)+EDRAFT*Tenv(1:NTX)
cc        TMOMDN(xymoms,1:NTX)= TMOMDN(xymoms,1:NTX)+ TMOM(xymoms,L
cc   *         ,1:NTX)*FENTRA
cc        DTMR(L,1:NTX)=DTMR(L,1:NTX)-EDRAFT*TENV(1:NTX)
cc        DTMOMR(:,L,1:NTX)=DTMOMR(:,L,1:NTX)-TMOM(:,L,1:NTX)*FENTRA
#endif
cc      ELSE  ! occasionally detrain into environment if ddraft too big
cc        FENTRA=EDRAFT/DDRUP  ! < 0
cc        DSM(L)=DSM(L)-FENTRA*SMDN
cc        DSMOM(xymoms,L)=DSMOM(xymoms,L)-SMOMDN(xymoms)*FENTRA
cc        DQM(L)=DQM(L)-FENTRA*QMDN
cc        DQMOM(xymoms,L)=DQMOM(xymoms,L)-QMOMDN(xymoms)*FENTRA
cc        SMDN=SMDN*(1+FENTRA)
cc        QMDN=QMDN*(1+FENTRA)
cc        SMOMDN(xymoms)= SMOMDN(xymoms)*(1+FENTRA)
cc        QMOMDN(xymoms)= QMOMDN(xymoms)*(1+FENTRA)
cc        DM(L)=DM(L)-EDRAFT
#ifdef TRACERS_ON
cc        DTM(L,1:NTX)=DTM(L,1:NTX)-FENTRA*TMDN(1:NTX)
cc        DTMOM(xymoms,L,1:NTX)=DTMOM(xymoms,L,1:NTX)-
cc   *         TMOMDN(xymoms,1:NTX)*FENTRA
cc        TMDN(1:NTX)=TMDN(1:NTX)*(1.+FENTRA)
cc        TMOMDN(xymoms,1:NTX)= TMOMDN(xymoms,1:NTX)*(1.+FENTRA)
#endif
cc      END IF
cc    ENDIF
      IF(L.GT.1) DDM(L-1)=DDRAFT
      LDMIN=L
C**** ALLOW FOR DOWNDRAFT TO DROP BELOW LMIN, IF IT'S NEGATIVE BUOYANT
      IF (L.LE.LMIN.AND.L.GT.1) THEN
        GO TO 345   ! force exit at LMIN
cc      SMIX=SMDN/DDRAFT
cc      IF (SMIX.GE.(SM1(L-1)/AIRM(L-1))) GO TO 345
      ENDIF
  346 CONTINUE
  345 DSM(LDMIN)=DSM(LDMIN)+SMDN
      DSMOM(xymoms,LDMIN)=DSMOM(xymoms,LDMIN) + SMOMDN(xymoms)
      DQM(LDMIN)=DQM(LDMIN)+QMDN
      DQMOM(xymoms,LDMIN)=DQMOM(xymoms,LDMIN) + QMOMDN(xymoms)
#ifdef TRACERS_ON
      DTM(LDMIN,1:NTX) = DTM(LDMIN,1:NTX) + TMDN(1:NTX)
      DTMOM(xymoms,LDMIN,1:NTX) = DTMOM(xymoms,LDMIN,1:NTX) +
     *     TMOMDN(xymoms,1:NTX)
#endif
      DO K=1,KMAX
      DUM(K,LDMIN)=DUM(K,LDMIN)+UMDN(K)
      DVM(K,LDMIN)=DVM(K,LDMIN)+VMDN(K)
      ENDDO
      DM(LDMIN)=DM(LDMIN)+DDRAFT
      END IF
C****
C**** SUBSIDENCE AND MIXING
C****
C**** Calculate vertical mass fluxes (Note CM for subsidence is defined
C**** in opposite sense than normal (positive is down))
      DO L=0,LDMIN-1
        CM(L) = 0.
      END DO
      DO L=LDMIN,LMAX
        CM(L) = CM(L-1) - DM(L) - DMR(L)
        SMT(L)=SM(L)    ! Save profiles for diagnostics
        QMT(L)=QM(L)
      END DO
      DO L=LMAX+1,LM
        CM(L) = 0.
      END DO
C**** simple upwind scheme for momentum
      ALPHA=0.
      DO 380 L=LDMIN,LMAX
      CLDM=CCM(L)
      IF(L.LT.LDRAFT.AND.ETADN.GT.1d-10) CLDM=CCM(L)-DDM(L)
      IF(MC1) VSUBL(L)=100.*CLDM*RGAS*TL(L)/(PL(L)*GRAV*DTsrc)
      BETA=CLDM*BYAM(L+1)
      IF(CLDM.LT.0.) BETA=CLDM*BYAM(L)
      BETAU=BETA
      ALPHAU=ALPHA
      IF(BETA.LT.0.) BETAU=0.
      IF(ALPHA.LT.0.) ALPHAU=0.
      DO K=1,KMAX
       UM(K,L)=
     *    UM(K,L)+RA(K)*(-ALPHAU*UM(K,L)+BETAU*UM(K,L+1)+DUM(K,L))
       VM(K,L)=
     *    VM(K,L)+RA(K)*(-ALPHAU*VM(K,L)+BETAU*VM(K,L+1)+DVM(K,L))
      ENDDO
  380 ALPHA=BETA
C**** Subsidence uses Quadratic Upstream Scheme for QM and SM
      ML(LDMIN:LMAX) = AIRM(LDMIN:LMAX) +   DMR(LDMIN:LMAX)
      SM(LDMIN:LMAX) =   SM(LDMIN:LMAX) +  DSMR(LDMIN:LMAX)
      SMOM(:,LDMIN:LMAX) =  SMOM(:,LDMIN:LMAX) + DSMOMR(:,LDMIN:LMAX)
C****
      nsub = lmax-ldmin+1
      cmneg(ldmin:lmax)=-cm(ldmin:lmax)
      cmneg(lmax)=0. ! avoid roundoff error (esp. for qlimit)
      call adv1d(sm(ldmin),smom(1,ldmin), f(ldmin),fmom(1,ldmin),
     &     ml(ldmin),cmneg(ldmin), nsub,.false.,1, zdir,ierrt,lerrt)
      SM(LDMIN:LMAX) =   SM(LDMIN:LMAX) +   DSM(LDMIN:LMAX)
      SMOM(:,LDMIN:LMAX) =  SMOM(:,LDMIN:LMAX) +  DSMOM(:,LDMIN:LMAX)
      ierr=max(ierrt,ierr) ; lerr=max(lerrt,lerr)
C****
      ML(LDMIN:LMAX) = AIRM(LDMIN:LMAX) +   DMR(LDMIN:LMAX)
      QM(LDMIN:LMAX) =   QM(LDMIN:LMAX) +  DQMR(LDMIN:LMAX)
      QMOM(:,LDMIN:LMAX) =  QMOM(:,LDMIN:LMAX) + DQMOMR(:,LDMIN:LMAX)
      call adv1d(qm(ldmin),qmom(1,ldmin), f(ldmin),fmom(1,ldmin),
     &     ml(ldmin),cmneg(ldmin), nsub,.true.,1, zdir,ierrt,lerrt)
      QM(LDMIN:LMAX) =   QM(LDMIN:LMAX) +   DQM(LDMIN:LMAX)
      QMOM(:,LDMIN:LMAX) =  QMOM(:,LDMIN:LMAX) +  DQMOM(:,LDMIN:LMAX)
      ierr=max(ierrt,ierr) ; lerr=max(lerrt,lerr)
C**** diagnostics
      DO L=LDMIN,LMAX
        FCDH=0.
        IF(L.EQ.LMAX) FCDH=CDHSUM-(CDHSUM-CDHDRT)*.5*ETADN+CDHM
        FCDH1=0.
        IF(L.EQ.LDMIN) FCDH1=(CDHSUM-CDHDRT)*.5*ETADN-EVPSUM
        MCFLX(L)=MCFLX(L)+CCM(L)
        DGDSM(L)=DGDSM(L)+PLK(L)*(SM(L)-SMT(L))-FCDH-FCDH1
        DTOTW(L)=DTOTW(L)+SLHE*(QM(L)-QMT(L)+COND(L))
        DGDQM(L)=DGDQM(L)+SLHE*(QM(L)-QMT(L))
      END DO
#ifdef TRACERS_ON
C**** Subsidence of tracers by Quadratic Upstream Scheme
      DO N=1,NTX
      ML(LDMIN:LMAX) =  AIRM(LDMIN:LMAX) +    DMR(LDMIN:LMAX)
      TM(LDMIN:LMAX,N) =  TM(LDMIN:LMAX,N) + DTMR(LDMIN:LMAX,N)
      TMOM(:,LDMIN:LMAX,N) = TMOM(:,LDMIN:LMAX,N)+DTMOMR(:,LDMIN:LMAX,N)
      call adv1d(tm(ldmin,n),tmom(1,ldmin,n), f(ldmin),fmom(1,ldmin),
     &     ml(ldmin),cmneg(ldmin), nsub,.true.,1, zdir,ierrt,lerrt)
      TM(LDMIN:LMAX,N) = TM(LDMIN:LMAX,N) +   DTM(LDMIN:LMAX,N)
      TMOM(:,LDMIN:LMAX,N) = TMOM(:,LDMIN:LMAX,N) +DTMOM(:,LDMIN:LMAX,N)
      ierr=max(ierrt,ierr) ; lerr=max(lerrt,lerr)
      END DO
#endif
C**** save new 'environment' profile for static stability calc.
      DO L=1,LM
        SM1(L)=SM(L)
        QM1(L)=QM(L)
      END DO
#ifdef TRACERS_ON
      TM1(:,1:NTX) = TM(:,1:NTX)
#endif
C****
C**** REEVAPORATION AND PRECIPITATION
C****
      IF(MC1.AND.PLE(LMIN)-PLE(LMAX+1).GE.450.) THEN
        DO L=LMAX,LMIN+4,-1
c          IF (PLE(L+1).GT.550.) EXIT
          SVWMXL(L)=FCLW*COND(L)*BYAM(L)
          COND(L)=COND(L)*(1.-FCLW)
#ifdef TRACERS_WATER
C**** Apportion cloud tracers and condensation
C**** Note that TRSVWML is in mass units unlike SVWMX
          TRSVWML(1:NTX,L) = FCLW*TRCOND(1:NTX,L)
          TRCOND(1:NTX,L) = (1.-FCLW)*TRCOND(1:NTX,L)
#endif
        END DO
      END IF
      PRCP=COND(LMAX)
      PRHEAT=CDHEAT(LMAX)
#ifdef TRACERS_WATER
C**** Tracer precipiation
C Note that all of the tracers that condensed do not precipitate here,
C since a fraction (FCLW) of TRCOND was removed above.
      TRPRCP(1:NTX) = TRCOND(1:NTX,LMAX)
#endif
         DPHASE(LMAX)=DPHASE(LMAX)+CDHSUM-(CDHSUM-CDHDRT)*.5*ETADN+CDHM
      DO 540 L=LMAX-1,1,-1
      IF(PRCP.LE.0.) GO TO 530
      FCLOUD=CCMUL*CCM(L)*BYAM(L+1)
      IF(PLE(LMIN)-PLE(L+2).GE.450.) FCLOUD=CCMUL1*CCM(L)*BYAM(L+1)
      IF(L.LT.LMIN) FCLOUD=CCMUL*CCM(LMIN)*BYAM(LMIN+1)
      IF(PLE(LMIN)-PLE(LMAX+1).LT.450.) THEN
        IF(L.EQ.LMAX-1) FCLOUD=CCMUL2*CCM(L)*BYAM(L+1)
        IF(L.LT.LMIN) FCLOUD=0.
      ENDIF
      IF(FCLOUD.GT.1.) FCLOUD=1.
      FEVAP=.5*CCM(L)*BYAM(L+1)
      IF(L.LT.LMIN) FEVAP=.5*CCM(LMIN)*BYAM(LMIN+1)
      IF(FEVAP.GT..5) FEVAP=.5
      CLDMCL(L+1)=CLDMCL(L+1)+FCLOUD
      CLDREF=CLDMCL(L+1)
      IF(PLE(LMAX+1).GT.700..AND.CLDREF.GT.CLDSLWIJ)
     *  CLDSLWIJ=CLDREF
      IF(PLE(LMIN)-PLE(LMAX+1).GE.450..AND.CLDREF.GT.CLDDEPIJ)
     *  CLDDEPIJ=CLDREF
C**** FORWARD STEP COMPUTES HUMIDITY CHANGE BY RECONDENSATION
C**** Q = Q + F(TOLD,PRHEAT,QOLD+EVAP)
      PRECNVL(L+1)=PRECNVL(L+1)+PRCP*BYGRAV
      MCLOUD=0.
      IF(L.LE.LMIN) MCLOUD=2.*FEVAP*AIRM(L)
      TOLD=SMOLD(L)*PLK(L)*BYAM(L)
      TOLD1=SMOLD(L+1)*PLK(L+1)*BYAM(L+1)
C**** decide whether to melt snow/ice
      HEAT1=0.
      IF(L.EQ.LFRZ.OR.(L.LE.LMIN.AND.TOLD.GE.TF.AND.TOLD1.LT.TF.AND.L.GT
     *     .LFRZ)) HEAT1=LHM*PRCP*BYSHA
C**** and deal with possible inversions and re-freezing of rain
      IF (TOLD.LT.TF.AND.TOLD1.GT.TF) HEAT1=-LHM*PRCP*BYSHA
      TNX=TOLD
      QNX=QMOLD(L)*BYAM(L)
      LHX=LHE
      IF(TNX.LT.TF) LHX=LHS
      IF(L.GT.LMIN) THEN    ! needed to avoid unintiallised value
        IF (TPSAV(L).GE.TF) LHX=LHE
      END IF
      SLH=LHX*BYSHA
      DQSUM=0.
      DO 510 N=1,3
      QSATC=QSAT(TNX,LHX,PL(L))
      DQ=(QSATC-QNX)/(1.+SLH*QSATC*DQSATDT(TNX,LHX))
      TNX=TNX-SLH*DQ
      QNX=QNX+DQ
  510 DQSUM=DQSUM+DQ*MCLOUD
      IF(DQSUM.LT.0.) DQSUM=0.
      IF(DQSUM.GT.PRCP) DQSUM=PRCP
      FPRCP=DQSUM/PRCP
      PRCP=PRCP-DQSUM
C**** UPDATE TEMPERATURE AND HUMIDITY DUE TO NET REVAPORATION IN CLOUDS
      FSSUM = 0
      IF (ABS(PLK(L)*SM(L)).gt.teeny) FSSUM = (SLH*DQSUM+HEAT1)/
     *     (PLK(L)*SM(L))
      SM(L)=SM(L)-(SLH*DQSUM+HEAT1)/PLK(L)
      SMOM(:,L) =  SMOM(:,L)*(1.-FSSUM)
      QM(L)=QM(L)+DQSUM
#ifdef TRACERS_WATER
C**** Tracer net re-evaporation
C**** (If 100% evaporation, allow all tracers to evaporate completely.)
      BELOW_CLOUD = L.lt.LMIN
      IF(FPRCP.eq.1.) THEN      !total evaporation
        DO N=1,NTX
          TM(L,N)   = TM(L,N)  + TRPRCP(N)
          TRPRCP(N) = 0.D0
        END DO
      ELSE ! otherwise, tracers evaporate dependent on type of tracer
C**** estimate effective humidity
        if (below_cloud) then
          TNX1=(SM(L)*PLK(L)-SLH*DQSUM*(1./(2.*MCLOUD)-1.))*BYAM(L)
          HEFF=MIN(1d0,(QM(L)+DQSUM*(1./(2.*MCLOUD)-1.))*BYAM(L)
     *         /QSAT(TNX1,LHX,PL(L)))
        else
          heff=1.
        end if
        DO N=1,NTX
          CALL GET_EVAP_FACTOR(N,TNX,LHX,BELOW_CLOUD,HEFF,FPRCP,FPRCPT)
          TM(L,N) = TM(L,N)     + FPRCPT*TRPRCP(N)
          TRPRCP(N) = TRPRCP(N) - FPRCPT*TRPRCP(N)
        END DO
      END IF
C**** CONDENSING and WASHOUT of TRACERS BELOW CLOUD
      IF(BELOW_CLOUD) THEN ! BELOW CLOUD
        WMXTR = PRCPMC*BYAM(L)
        precip_mm = PRCPMC*100.*bygrav
        b_beta_DT = FPLUME
        DO N=1,NTX
          CALL GET_COND_FACTOR(L,N,WMXTR,TNX,LHX,FPLUME,0d0,FQCONDT)
          CALL GET_WASH_FACTOR(N,b_beta_DT,precip_mm,FWASHT)
          TRCOND(N,L) = FPLUME * FQCONDT * TM(L,N)
          TRPRCP(N)=TRPRCP(N) + TM(L,N)*FWASHT
          TM(L,N)=TM(L,N) * (1.-FPLUME*FQCONDT-FWASHT)
          TMOM(xymoms,L,N)=TMOM(xymoms,L,N) *
     &                       (1.-FPLUME*FQCONDT-FWASHT)
        END DO
      END IF
#endif
         FCDH1=0.
         IF(L.EQ.LDMIN) FCDH1=(CDHSUM-CDHDRT)*.5*ETADN-EVPSUM
         DPHASE(L)=DPHASE(L)-SLH*DQSUM+FCDH1-HEAT1
         DQCOND(L)=DQCOND(L)-SLH*DQSUM
C**** ADD PRECIPITATION AND LATENT HEAT BELOW
  530 PRHEAT=CDHEAT(L)+SLH*PRCP
      PRCP=PRCP+COND(L)
#ifdef TRACERS_WATER
      TRPRCP(1:NTX) = TRPRCP(1:NTX) + TRCOND(1:NTX,L)
#ifdef TRACERS_SPECIAL_O18
C**** Isotopic equilibration of liquid precip with water vapour
      IF (LHX.eq.LHE .and. PRCP.gt.0) THEN
        DO N=1,NTX
          CALL ISOEQUIL(NTIX(N),TNX,.TRUE.,QM(L),PRCP,TM(L,N),TRPRCP(N)
     *         ,0.5d0)
        END DO
      END IF
#endif
#endif
  540 CONTINUE
C****
      IF(PRCP.GT.0.) CLDMCL(1)=CLDMCL(1)+CCM(LMIN)*BYAM(LMIN+1)
      PRCPMC=PRCPMC+PRCP
#ifdef TRACERS_WATER
      TRPRMC(1:NTX) = TRPRMC(1:NTX) + TRPRCP(1:NTX)
#endif
      IF(LMCMIN.GT.LDMIN) LMCMIN=LDMIN
C****
C**** END OF LOOP OVER CLOUD TYPES
C****
  570 CONTINUE
C****
C**** END OF OUTER LOOP OVER CLOUD BASE
C****
  600 CONTINUE
      IF(LMCMIN.GT.0) THEN
C**** ADJUSTMENT TO CONSERVE CP*T
        SUMAJ=0.
        SUMDP=0.
        DO L=LMCMIN,LMCMAX
          SUMDP=SUMDP+AIRM(L)
          SUMAJ=SUMAJ+DGDSM(L)
        END DO
        DO L=LMCMIN,LMCMAX
          DGDSM(L)=DGDSM(L)-SUMAJ*AIRM(L)/SUMDP
          SM(L)=SM(L)-SUMAJ*AIRM(L)/(SUMDP*PLK(L))
        END DO
C**** LOAD MASS EXCHANGE ARRAY FOR GWDRAG
        AIRXL = 0.
        DO L=LMCMIN,LMCMAX
          AIRXL = AIRXL+MCFLX(L)
        END DO
      END IF
C**** CALCULATE OPTICAL THICKNESS
      WCONST=WMU*(1.-PEARTH)+WMUL*PEARTH
      WMSUM=0.
      DO L=1,LMCMAX
         TL(L)=(SM(L)*BYAM(L))*PLK(L)
         TEMWM=(TAUMCL(L)-SVWMXL(L)*AIRM(L))*1.d2*BYGRAV
         IF(TL(L).GE.TF) WMSUM=WMSUM+TEMWM
         IF(CLDMCL(L).GT.0.) THEN
               TAUMCL(L)=AIRM(L)*COETAU
            IF(L.EQ.LMCMAX .AND. PLE(LMCMIN)-PLE(LMCMAX+1).LT.450.)
     *         TAUMCL(L)=AIRM(L)*.02d0
            IF(L.LE.LMCMIN .AND. PLE(LMCMIN)-PLE(LMCMAX+1).GE.450.)
     *         TAUMCL(L)=AIRM(L)*.02d0
         END IF
         IF(SVLATL(L).EQ.0.) THEN
            SVLATL(L)=LHE
            IF ( (SVTP(L).gt.0. .and. SVTP(L).LT.TF) .or.
     *           (SVTP(L).eq.0. .and. TL(L).lt.TF) ) SVLATL(L)=LHS
         ENDIF
         IF(SVWMXL(L).GT.0.) THEN
            FCLD=CLDMCL(L)+1.E-20
            TEM=1.d5*SVWMXL(L)*AIRM(L)*BYGRAV
            WTEM=1.d5*SVWMXL(L)*PL(L)/(FCLD*TL(L)*RGAS)
            IF(SVLATL(L).EQ.LHE.AND.SVWMXL(L)/FCLD.GE.WCONST*1.d-3)
     *           WTEM=1d2*WCONST*PL(L)/(TL(L)*RGAS)
            IF(SVLATL(L).EQ.LHS.AND.SVWMXL(L)/FCLD.GE.WMUI*1.d-3)
     *           WTEM=1d2*WMUI*PL(L)/(TL(L)*RGAS)
            IF(WTEM.LT.1.d-10) WTEM=1.d-10
            IF(SVLATL(L).EQ.LHE)  THEN
               RCLD=(10.*(1.-PEARTH)+7.0*PEARTH)*(WTEM*4.)**BY3
             ELSE
               RCLD=25.0*(WTEM/4.2d-3)**BY3 * 500./(pl(l)+500.)
            END IF
            RCLDE=RCLD/BYBR   !  effective droplet radius in anvil
            CSIZEL(L)=RCLDE   !  effective droplet radius in anvil
            TAUMCL(L)=1.5*TEM/(FCLD*RCLDE+1.E-20)
            IF(TAUMCL(L).GT.100.) TAUMCL(L)=100.
         END IF
      END DO

      RETURN
      END SUBROUTINE MSTCNV

      SUBROUTINE LSCOND(IERR,WMERR,LERR)
!@sum  LSCOND column physics of large scale condensation
!@auth M.S.Yao/A. Del Genio (modularisation by Gavin Schmidt)
!@ver  1.0 (taken from CB265)
!@calls CTMIX,QSAT,DQSATDT,THBAR
      IMPLICIT NONE

!@var IERR,WMERR,LERR error reporting
      INTEGER, INTENT(OUT) :: IERR,LERR
      REAL*8, INTENT(OUT) :: WMERR
      REAL*8 LHX

C**** functions
      REAL*8 :: QSAT, DQSATDT,ERMAX
!@var QSAT saturation humidity
!@var DQSATDT dQSAT/dT

!@param CM00 upper limit for autoconversion rate
!@param AIRM0 scaling factor for computing rain evaporation
!@param HEFOLD e-folding length for computing HRISE
!@param COEFM coefficient for ratio of cloud water amount and WCONST
!@param COEFT coefficient used in computing PRATM
!@param COESIG coefficient for equ. 23 of Del Genio et al. (1996)
!@param COEEC coefficient for computing cloud evaporation
C**** Adjust COEFT and COEFM to change proportion of super-cooled rain
C**** to snow. Increasing COEFT reduces temperature range of super
C**** -cooled rain, increasing COEFM enhances probability of snow.
      REAL*8, PARAMETER :: CM00=1.d-4, AIRM0=100.d0, GbyAIRM0=GRAV/AIRM0
      REAL*8, PARAMETER :: HEFOLD=500.,COEFM=10.,COEFT=2.5
      REAL*8, PARAMETER :: COESIG=1d-3,COEEC=1000.
      REAL*8, DIMENSION(IM) :: UMO1,UMO2,UMN1,UMN2 !@var dummy variables
      REAL*8, DIMENSION(IM) :: VMO1,VMO2,VMN1,VMN2 !@var dummy variables
!@var Miscellaneous vertical arrays
      REAL*8, DIMENSION(LM) ::
     *     QSATL,RHF,ATH,SQ,ER,QHEAT,
     *     CAREA,PREP,RH00,EC,WMXM
!@var QSATL saturation water vapor mixing ratio
!@var RHF environmental relative humidity
!@var ATH change in potential temperature
!@var SQ ERMAX dummy variables
!@var ER evaporation of precip
!@var QHEAT change of latent heat
!@var CAREA fraction of clear region
!@var PREP precip conversion rate
!@var RH00 threshold relative humidity
!@var EC cloud evaporation rate
!@var WMXM cloud water mass (mb)
      REAL*8, DIMENSION(LM+1) :: PREBAR,PREICE
!@var PREBAR,PREICE precip entering layer top for total, snow

#ifdef TRACERS_WATER
!@var TRPRBAR tracer precip entering layer top for total (kg)
      REAL*8, DIMENSION(NTM,LM+1) :: TRPRBAR
!@var DTPR change of tracer by precip (kg)
!@var DTER change of tracer by evaporation (kg)
!@var DTQW change of tracer by condensation (kg)
!@var FWTOQ fraction of CLW that goes to water vapour
!@var FPR fraction of CLW that precipitates
!@var FER fraction of precipitate that evaporates
      REAL*8 DTPR,DTER,DTQW,TWMTMP,DTSUM,FWTOQ,FPR,FER
!@var DTPRT tracer-specific change of tracer by precip (kg)
!@var DTERT tracer-specific change of tracer by evaporation (kg)
!@var DTWRT tracer-specific change of tracer by washout (kg)
!@var DTQWT tracer-specific change of tracer by condensation (kg)
!@var FWTOQT tracer-specific fraction of tracer in CLW that evaporates
!@var FQTOWT tracer-specific fraction of gas tracer condensing in CLW
!@var FPRT tracer-specific fraction of tracer in CLW that precipitates
!@var FERT tracer-specific fraction of tracer in precipitate evaporating
      REAL*8 ::DTWRT,DTPRT,DTERT,DTQWT,FWTOQT,FQTOWT,FPRT,FERT,PRLIQ
!@var BELOW_CLOUD logical- is the current level below cloud?
!@var CLOUD_YET logical- in L loop, has there been any cloud so far?
      LOGICAL BELOW_CLOUD,CLOUD_YET
#endif

      REAL*8 AIRMR,BETA,BMAX
     *     ,CBF,CBFC0,CK,CKIJ,CK1,CK2,CKM,CKR,CM,CM0,CM1,DFX,DQ,DQSDT
     *     ,DQSUM,DQUP,DRHDT,DSE,DSEC,DSEDIF,DWDT,DWDT1,DWM,ECRATE,EXPST
     *     ,FCLD,FMASS,FMIX,FPLUME,FPMAX,FQTOW,FRAT,FUNI,HRISE
     *     ,FUNIL,FUNIO,HCHANG,HDEP,HPHASE,OLDLAT,OLDLHX,PFR,PMI,PML
     *     ,HDEP1,PRATIO,QCONV,QHEATC,QLT1,QLT2,QMN1,QMN2,QMO1,QMO2,QNEW
     *     ,QNEWU,QOLD,QOLDU,QSATC,QSATE,RANDNO,RCLDE,RHI,RHN,RHO,RHT1
     *     ,RHW,SEDGE,SIGK,SLH,SMN1,SMN2,SMO1,SMO2,TEM,TEMP,TEVAP,THT1
     *     ,THT2,TLT1,TNEW,TNEWU,TOLD,TOLDU,TOLDUP,VDEF,WCONST,WMN1,WMN2
     *     ,WMNEW,WMO1,WMO2,WMT1,WMT2,WMX1,WTEM,VVEL,XY,RCLD,FCOND,HDEPx
     *     ,PRATW,PRATM
!@var BETA,BMAX,CBFC0,CKIJ,CK1,CK2,HRISE,PRATW,PRATM dummy variables
!@var AIRMR
!@var CBF enhancing factor for precip conversion
!@var CK ratio of cloud top jumps in moist static energy and total water
!@var CKM CTEI threshold in MacVean and Mason theory
!@var CKR CTEI threshold in Randall theory
!@var CM conversion rate for large cloud water content
!@var CM0,CM1 limiting autoconversion rate for large cloud water content
!@var DFX iteration increment
!@var DQ condensed water vapor
!@var DQSDT derivative of saturation vapor pressure w.r.t. temperature
!@var DQSUM,DQUP dummy variables
!@var DRHDT time change of relative humidity
!@var DSE moist static energy jump at cloud top
!@var DSEC critical DSE for CTEI to operate
!@var DSFDIF DSE-DSEC
!@var DWDT,DWDT1 time change rates of cloud water
!@var ECRATE cloud droplet evaporation rate
!@var EXPST exponential term in determining the fraction in CTEI
!@var FCLD cloud fraction
!@var FMIX,FRAT fraction of mixing in CTEI
!@var FMASS mass of mixing in CTEI
!@var FPLUME fraction of mixing in CTEI
!@var FPMAX max fraction of mixing in CTEI
!@var FQTOW fraction of water vapour that goes to CLW
!@var FUNI the probablity for ice cloud to form
!@var FUNIL,FUNIO FUNI over land, ocean
!@var HCHANG,HPHASE latent heats for changing phase
!@var HDEP,HDEP1 layer depth (m)  (note change of unit km-->m)
!@var OLDLAT,OLDLHX previous LHX
!@var PFR PROBABLITY OF GLACIATION OF SUPER-COOLED WATER
!@var PMI icy precip entering the layer top
!@var PML layer's cloud water devided by GRAV
!@var PRATIO PMI/PML
!@var QCONV convergence of latent heat
!@var QHEATC,QLT1,QLT2,QMN1,QMN2,QMO1,QMO2 dummy variables
!@var QNEW,QNEWU updated specific humidity
!@var QOLD,QOLDU previous specific humidity
!@var QSATC saturation vapor mixing ratio
!@var QSATE saturation vapor mixing ratio w.r.t. water
!@var RANDNO random number
!@var RCLD,RCLDE cloud particle's radius, effective radius
!@var RHI relative humidity w.r.t. ice
!@var RHN dummy variable
!@var RHO air density
!@var RHT1,RHW,SIGK dummy variables
!@var SEDGE potential temperature at layer edge
!@var SMN1,SMN2,SMO1,SMO2,TEM,TEMP,TEVAP,THT1,THT2,TLT1 dummy variables
!@var TNEW,TNEWU updated tempertures
!@var TOLD,TOLDU,TOLDUP previous temperatures
!@var VDEF = VVEL - VSUB
!@var WCONST,WMN1,WMN2,WMO1,WMO2,WMT1,WMT2,WMX1 dummy variables
!@var WMNEW updated cloud water mixing ratio
!@var WTEM cloud water density (g m**-3)
!@var VVEL vertical velocity (cm/s)
!@var XY,FCOND dummy variables
      INTEGER LN,ITER !@var LN,ITER loop variables
      LOGICAL BANDF  !@var BANDF true if Bergero-Findeisen proc. occurs
      LOGICAL FORM_CLOUDS !@var FORM_CLOUDS true if clouds are formed

      INTEGER K,L,N  !@var K,L,N loop variables

      REAL*8 THBAR !@var THBAR potential temperature at layer edge

C****
C**** LARGE-SCALE CLOUDS AND PRECIPITATION
C**** THE LIQUID WATER CONTENT IS PREDICTED
C****
      IERR=0
C****
      PRCPSS=0.
      HCNDSS=0.
      CKIJ=1.
C**** initialise vertical arrays
      ER=0.
      EC=0.
      PREP=0.
      PREBAR=0.
      LHP=0.
      QHEAT=0.
      CLDSSL=0
      TAUSSL=0
#ifdef TRACERS_WATER
      TRPRSS = 0.
      TRPRBAR = 0.
      BELOW_CLOUD=.false.
      CLOUD_YET=.false.
#endif
      DO L=1,LP50
        CAREA(L)=1.-CLDSAVL(L)
        IF(WMX(L).LE.0.) CAREA(L)=1.
      END DO
      DQUP=0.
      TOLDUP=TL(LP50)
      PREICE(LP50+1)=0.
      WCONST=WMU*(1.-PEARTH)+WMUL*PEARTH
         SSHR=0.
         DCTEI=0.
C****
C**** MAIN L LOOP FOR LARGE-SCALE CONDENSATION, PRECIPITATION AND CLOUDS
C****
      DO L=LP50,1,-1
      TOLD=TL(L)
      QOLD=QL(L)
      OLDLHX=SVLHXL(L)
      OLDLAT=SVLATL(L)
C**** COMPUTE VERTICAL VELOCITY IN CM/S
      TEMP=100.*RGAS*TL(L)/(PL(L)*GRAV)
      IF(L.EQ.1)  THEN
         VVEL=-SDL(L+1)*TEMP
      ELSE IF(L.EQ.LM)  THEN
         VVEL=-SDL(L)*TEMP
      ELSE
         VVEL=-.5*(SDL(L)+SDL(L+1))*TEMP
      END IF
C**** COMPUTE THE LIMITING AUTOCONVERSION RATE FOR CLOUD WATER CONTENT
      CM0=CM00
      VDEF=VVEL-VSUBL(L)
      IF(VDEF.GT.0.) CM0=CM00*10.**(-VDEF)
      FCLD=1.-CAREA(L)+1.E-20
C**** COMPUTE THE PROBABILITY OF ICE FORMATION, FUNI, AND
C**** THE PROBABLITY OF GLACIATION OF SUPER-COOLED WATER, PFR
C**** DETERMINE THE PHASE MOISTURE CONDENSES TO
C**** DETERMINE THE POSSIBILITY OF B-F PROCESS
      BANDF=.FALSE.
      LHX=LHE
      CBF=1. + EXP(-((TL(L)-258.16d0)/10.)**2)

      IF (TL(L).LE.TI) THEN     ! below -40: force ice
        LHX=LHS
      ELSEIF (TL(L).GE.TF) THEN ! above freezing: force water
        LHX=LHE
      ELSE                      ! in between: compute probability
        IF(TL(L).GT.269.16) THEN ! OC/SI/LI clouds: water above -4
          FUNIO=0.
        ELSE
          FUNIO=1.-EXP(-((TL(L)-269.16d0)/15.)**2)
        END IF
        IF(TL(L).GT.263.16) THEN ! land clouds water: above -10
          FUNIL=0.
        ELSE
          FUNIL=1.-EXP(-((TL(L)-263.16d0)/15.)**2)
        END IF
        FUNI=FUNIO*(1.-PEARTH)+FUNIL*PEARTH
        RANDNO=RNDSSL(1,L)       !  RANDNO=RANDU(XY)
        IF(RANDNO.LT.FUNI) LHX=LHS

C**** special case 1) if ice previously then stay as ice (if T<Tf)
        IF((OLDLHX.EQ.LHS.OR.OLDLAT.EQ.LHS).AND.TL(L).LT.TF) THEN
          IF(LHX.EQ.LHE) BANDF=.TRUE.
          LHX=LHS
        ENDIF

        IF (L.LT.LP50) THEN
C**** Decide whether precip initiates B-F process
          PML=WMX(L)*AIRM(L)*BYGRAV
          PMI=PREICE(L+1)*DTsrc
          RANDNO=RNDSSL(2,L)     !  RANDNO=RANDU(XY)
C**** Calculate probability of ice precip seeding a water cloud
          IF (LHX.EQ.LHE.AND.PMI.gt.0) THEN
            PRATIO=MIN(PMI/(PML+1.E-20),10d0)
            CBFC0=.5*CM0*CBF*DTsrc
            PFR=(1.-EXP(-(PRATIO*PRATIO)))*(1.-EXP(-(CBFC0*CBFC0)))
            IF(PFR.GT.RANDNO) THEN
              BANDF=.TRUE.
              LHX=LHS
            END IF
          END IF
C**** If liquid rain falls into an ice cloud, B-F must occur
          IF (LHP(L+1).EQ.LHE .AND. LHX.EQ.LHS .AND. PML.GT.0.)
     *         BANDF=.TRUE.
        END IF
      END IF
      IF(LHX.EQ.LHS .AND. (OLDLHX.EQ.LHE.OR.OLDLAT.EQ.LHE)) BANDF=.TRUE.

C**** COMPUTE RELATIVE HUMIDITY
      QSATL(L)=QSAT(TL(L),LHX,PL(L))
      RH1(L)=QL(L)/QSATL(L)
      IF(LHX.EQ.LHS) THEN
        QSATE=QSAT(TL(L),LHE,PL(L))
        RHW=.00536d0*TL(L)-.276d0
        RH1(L)=QL(L)/QSATE
        IF(TL(L).LT.238.16) RH1(L)=QL(L)/(QSATE*RHW)
      END IF
C**** PHASE CHANGE OF CLOUD WATER CONTENT
      HCHANG=0.
      IF(OLDLHX.EQ.LHE.AND.LHX.EQ.LHS) HCHANG= WML(L)*LHM
      IF(OLDLHX.EQ.LHS.AND.LHX.EQ.LHE) HCHANG=-WML(L)*LHM
      IF(OLDLAT.EQ.LHE.AND.LHX.EQ.LHS) HCHANG=HCHANG+SVWMXL(L)*LHM
      IF(OLDLAT.EQ.LHS.AND.LHX.EQ.LHE) HCHANG=HCHANG-SVWMXL(L)*LHM
      SVLHXL(L)=LHX
      TL(L)=TL(L)+HCHANG/SHA
      TH(L)=TL(L)/PLK(L)
      ATH(L)=(TH(L)-TTOLDL(L))*BYDTsrc
C**** COMPUTE RH IN THE CLOUD-FREE AREA, RHF
      RHI=QL(L)/QSAT(TL(L),LHS,PL(L))
    ! this formulation is used for consistency with current practice
      RH00(L)=U00wtrX*U00ice
      IF(LHX.EQ.LHS) RH00(L)=U00ice
      IF(L.EQ.1) THEN
        HDEP=AIRM(L)*TL(L)*RGAS/(GRAV*PL(L))
        RH00(L)=1.-GAMD*LHE*HDEP/(RVAP*TS*TS)
        IF(DCL.LE.1) THEN
         IF(RIS.GT.1.) HDEP1=10d0
         IF(RIS.LE.1..AND.RI1.GT.1.) HDEP1=50d0
         IF(RIS.LE.1..AND.RI1.LE.1..AND.RI2.GT.1.) HDEP1=100d0
         IF(RIS.LE.1..AND.RI1.LE.1..AND.RI2.LE.1.) HDEP1=HDEP
         RH00(L)=1.-GAMD*LHE*HDEP1/(RVAP*TS*TS)
        ENDIF
        IF(RH00(L).LT.0.) RH00(L)=0.
      ENDIF
C**** Special formulation for PBL layers
cc    IF(L.GT.1.AND.L.LE.DCL) THEN
C**** integrated HDEP over pbl depth
cc      HDEP=0.
cc      DO LN=1,DCL
cc        HDEP=HDEP+AIRM(LN)*TL(LN)*RGAS/(GRAV*PL(LN))
cc      END DO
c       IF(L.EQ.DCL) HDEP=HDEP+0.5*AIRM(L)*TL(L)*RGAS/(GRAV*PL(L))
cc      HRISE=HRMAX
cc      IF(HDEP.LT.10.*HEFOLD) HRISE=HRMAX*(1.-EXP(-HDEP/HEFOLD))
cc      IF(HRISE.GT.HDEP) HRISE=HDEP
C**** hdep is simply layer dependent (not used: resolution sensitive)
c        HDEP=AIRM(L)*TL(L)*RGAS/(GRAV*PL(L))
c        IF(L.EQ.DCL) HDEP=0.5*HDEP
cc      RH00(L)=max( RH00(L) , 1.-GAMD*LHE*HRISE/(RVAP*TS*TS) )
c       IF(RH00(L).LT.0.) RH00(L)=0.
cc    ENDIF
C****
      IF(RH00(L).GT.1.) RH00(L)=1.
      RHF(L)=RH00(L)+(1.-CAREA(L))*(1.-RH00(L))
C**** Set precip phase to be the same as the cloud, unless precip above
C**** is ice and temperatures are still below freezing
      LHP(L)=LHX
      IF (LHP(L+1).eq.LHS .and. TL(L).lt.TF) LHP(L)=LHP(L+1)
C**** COMPUTE THE AUTOCONVERSION RATE OF CLOUD WATER TO PRECIPITATION
      IF(WMX(L).GT.0.) THEN
        RHO=1d5*PL(L)/(RGAS*TL(L))
        TEM=RHO*WMX(L)/(WCONST*FCLD+teeny)
        IF(LHX.EQ.LHS ) TEM=RHO*WMX(L)/(WMUI*FCLD+teeny)
        TEM=TEM*TEM
        IF(TEM.GT.10.) TEM=10.
        CM1=CM0
        IF(BANDF) CM1=CM0*CBF
        IF(LHX.EQ.LHS) CM1=CM0
        CM=CM1*(1.-1./EXP(TEM*TEM))+1.*100.*(PREBAR(L+1)+
     *       PRECNVL(L+1)*BYDTsrc)
        IF(CM.GT.BYDTsrc) CM=BYDTsrc
        PREP(L)=WMX(L)*CM
        IF(TL(L).LT.TF.AND.LHX.EQ.LHE) THEN ! check snowing pdf
          PRATM=1.d5*COEFM*WMX(L)*PL(L)/(WCONST*FCLD*TL(L)*RGAS+teeny)
          PRATM=MIN(PRATM,1.D0)*(1.-EXP(MAX(-1d2,(TL(L)-TF)/COEFT)))
          IF(PRATM.GT.RNDSSL(3,L)) LHP(L)=LHS
        END IF
      ELSE
        CM=0.
      END IF
C**** DECIDE WHETHER TO FORM CLOUDS
C**** FORM CLOUDS ONLY IF RH GT RH00
      IF (RH1(L).LT.RH00(L)) THEN
        FORM_CLOUDS=.FALSE.
      ELSE   ! COMPUTE THE CONVERGENCE OF AVAILABLE LATENT HEAT
        SQ(L)=LHX*QSATL(L)*DQSATDT(TL(L),LHX)*BYSHA
        TEM=-LHX*DPDT(L)/PL(L)
        QCONV=LHX*AQ(L)-RH(L)*SQ(L)*SHA*PLK(L)*ATH(L)
     *       -TEM*QSATL(L)*RH(L)
        FORM_CLOUDS= (QCONV.GT.0. .OR. WMX(L).GT.0.)
      END IF
C****
      ERMAX=LHX*PREBAR(L+1)*GRAV*BYAM(L)
      IF (FORM_CLOUDS) THEN
C**** COMPUTE EVAPORATION OF RAIN WATER, ER
        RHN=MIN(RH(L),RHF(L))
        IF(WMX(L).GT.0.)  THEN
          ER(L)=(1.-RHN)*LHX*PREBAR(L+1)*GRAV*BYAM(L) !*GbyAIRM0 ! GRAV/AIRM0
        ELSE                    !  WMX(l).le.0.
          IF(PREICE(L+1).GT.0..AND.TL(L).LT.TF)  THEN
            ER(L)=(1.-RHI)*LHX*PREBAR(L+1)*GRAV*BYAM(L) ! *GbyAIRM0 ! GRAV/AIRM0
          ELSE
            ER(L)=(1.-RH(L))*LHX*PREBAR(L+1)*GRAV*BYAM(L) !*GbyAIRM0 ! GRAV/AIRM0
          END IF
        END IF
        ER(L)=MAX(0d0,MIN(ER(L),ERMAX))
C**** COMPUTATION OF CLOUD WATER EVAPORATION
        IF (CAREA(L).GT.0.) THEN
          WTEM=1d5*WMX(L)*PL(L)/(FCLD*TL(L)*RGAS+teeny)
          IF(LHX.EQ.LHE.AND.WMX(L)/FCLD.GE.WCONST*1d-3)
     *         WTEM=1d2*WCONST*PL(L)/(TL(L)*RGAS)
          IF(LHX.EQ.LHS.AND.WMX(L)/FCLD.GE.WMUI*1d-3)
     *         WTEM=1d2*WMUI*PL(L)/(TL(L)*RGAS)
          IF(WTEM.LT.1d-10) WTEM=1d-10
          IF(LHX.EQ.LHE)  THEN
            RCLD=1d-6*(10.*(1.-PEARTH)+7.0*PEARTH)*(WTEM*4.)**BY3
          ELSE
            RCLD=25.d-6*(WTEM/4.2d-3)**BY3 * 500./(pl(l)+500.)
          END IF
          CK1=1000.*LHX*LHX/(2.4d-2*RVAP*TL(L)*TL(L))
          CK2=1000.*RGAS*TL(L)/(2.4d-3*QSATL(L)*PL(L))
          TEVAP=COEEC*(CK1+CK2)*RCLD*RCLD
          WMX1=WMX(L)-PREP(L)*DTsrc
          ECRATE=(1.-RHF(L))/(TEVAP*FCLD+teeny)
          IF(ECRATE.GT.BYDTsrc) ECRATE=BYDTsrc
          EC(L)=WMX1*ECRATE*LHX
        END IF
C**** COMPUTE NET LATENT HEATING DUE TO STRATIFORM CLOUD PHASE CHANGE,
C**** QHEAT, AND NEW CLOUD WATER CONTENT, WMNEW
        DRHDT=2.*CAREA(L)*CAREA(L)*(1.-RH00(L))*(QCONV+ER(L)+EC(L))/LHX/
     *       (WMX(L)/(FCLD+teeny)+2.*CAREA(L)*QSATL(L)*(1.-RH00(L))
     *       +teeny)
        IF(ER(L).EQ.0..AND.EC(L).EQ.0..AND.WMX(L).LE.0.) DRHDT=0.
        QHEAT(L)=(QCONV-LHX*DRHDT*QSATL(L))/(1.+RH(L)*SQ(L))
        DWDT=QHEAT(L)/LHX-PREP(L)+CAREA(L)*ER(L)/LHX
        WMNEW =WMX(L)+DWDT*DTsrc
        IF(WMNEW.LT.0.) THEN
          WMNEW=0.
          QHEAT(L)=(-WMX(L)*BYDTsrc+PREP(L))*LHX-CAREA(L)*ER(L)
        END IF
      ELSE
C**** UNFAVORABLE CONDITIONS FOR CLOUDS TO EXIT, PRECIP OUT CLOUD WATER
        IF(WMX(L).GT.0.) PREP(L)=WMX(L)*BYDTsrc
        ER(L)=(1.-RH(L))*LHX*PREBAR(L+1)*GRAV*BYAM(L) !*GbyAIRM0 ! GRAV/AIRM0
        IF(PREICE(L+1).GT.0..AND.TL(L).LT.TF)
     *       ER(L)=(1.-RHI)*LHX*PREBAR(L+1)*GRAV*BYAM(L) !*GbyAIRM0 ! GRAV/AIRM0
        ER(L)=MAX(0d0,MIN(ER(L),ERMAX))
        QHEAT(L)=-CAREA(L)*ER(L)
        WMNEW=0.
      END IF

C**** PHASE CHANGE OF PRECIPITATION, FROM ICE TO WATER
C**** This occurs if 0 C isotherm is crossed, or ice is falling into
C**** a super-cooled water cloud (that has not had B-F occur).
C**** Note: on rare occasions we have ice clouds even if T>0
C**** In such a case, no energy of phase change is needed.
      HPHASE=0.
      IF (LHP(L+1).EQ.LHS.AND.LHP(L).EQ.LHE.AND.PREICE(L+1).GT.0) THEN
c      IF(L.LT.LM.AND.TL(L).GT.TF) THEN    ! back to SI2000
        HPHASE=HPHASE+LHM*PREICE(L+1)*GRAV*BYAM(L)
        PREICE(L+1)=0.
      ENDIF
C**** PHASE CHANGE OF PRECIP, FROM WATER TO ICE
      IF (LHP(L+1).EQ.LHE.AND.LHP(L).EQ.LHS.AND.PREBAR(L+1).GT.0)
     *     HPHASE=HPHASE-LHM*PREBAR(L+1)*GRAV*BYAM(L)
C**** Make sure energy is conserved for transfers between P and CLW
      IF (LHP(L).NE.LHX)
     *     HPHASE=HPHASE+(ER(L)*CAREA(L)/LHX-PREP(L))*LHM
C**** COMPUTE THE PRECIP AMOUNT ENTERING THE LAYER TOP
      IF (ER(L).eq.ERMAX) THEN ! to avoid round off problem
        PREBAR(L)=PREBAR(L+1)*(1.-CAREA(L))+AIRM(L)*PREP(L)*BYGRAV
      ELSE
        PREBAR(L)=MAX(0d0,PREBAR(L+1)+
     *       AIRM(L)*(PREP(L)-ER(L)*CAREA(L)/LHX)*BYGRAV)
      END IF
C**** UPDATE NEW TEMPERATURE AND SPECIFIC HUMIDITY
      QNEW =QL(L)-DTsrc*QHEAT(L)/LHX
      IF(QNEW.LT.0.) THEN
        QNEW=0.
        QHEAT(L)=QL(L)*LHX*BYDTsrc
        DWDT1=QHEAT(L)/LHX-PREP(L)+CAREA(L)*ER(L)/LHX
        WMNEW=WMX(L)+DWDT1*DTsrc
C**** IF WMNEW .LT. 0., THE COMPUTATION IS UNSTABLE
        IF(WMNEW.LT.0.) THEN
          IERR=1
          LERR=L
          WMERR=WMNEW
          WMNEW=0.
        END IF
      END IF
C**** Only Calculate fractional changes of Q to W
#ifdef TRACERS_WATER
      FPR=0.
      IF (WMX(L).gt.0.) FPR=PREP(L)*DTsrc/WMX(L)              ! CLW->P
      FPR=MIN(1d0,FPR)
      FER=0.
      IF (PREBAR(L+1).gt.0.) FER=CAREA(L)*ER(L)*AIRM(L)/(
     *     GRAV*LHX*PREBAR(L+1))                              ! P->Q
      FER=MIN(1d0,FER)
      FWTOQ=0.                                                ! CLW->Q
#endif
      FQTOW=0.                                                ! Q->CLW
      IF (QHEAT(L)+CAREA(L)*ER(L).gt.0) THEN
        IF (LHX*QL(L)+DTsrc*CAREA(L)*ER(L).gt.0.) FQTOW=(QHEAT(L)
     *       +CAREA(L)*ER(L))*DTsrc/(LHX*QL(L)+DTsrc*CAREA(L)*ER(L))
#ifdef TRACERS_WATER
      ELSE
        IF (WMX(L)-PREP(L)*DTsrc.gt.0.) FWTOQ=-(QHEAT(L)+CAREA(L)
     *       *ER(L))*DTsrc/(LHX*(WMX(L)-PREP(L)*DTsrc))
        FWTOQ=MIN(1d0,FWTOQ)
#endif
      END IF
      QL(L)=QNEW
C**** adjust gradients down if Q decreases
      QMOM(:,L)= QMOM(:,L)*(1.-FQTOW)
      WMX(L)=WMNEW
      TL(L)=TL(L)+DTsrc*(QHEAT(L)-HPHASE)/SHA
      TH(L)=TL(L)/PLK(L)
      TNEW=TL(L)
      QSATC=QSAT(TL(L),LHX,PL(L))
      RH(L)=QL(L)/QSATC
#ifdef TRACERS_WATER
C**** update tracers from cloud formation (in- and below-cloud
C****    precipitation, evaporation, condensation, and washout)
      DO N=1,NTX
c ---------------------- initialize fractions ------------------------
        FPRT  =0.
        FERT  =0.
        FWASHT=0.
        FQTOWT=0.
        FWTOQT=0.
c ----------------------- calculate fractions --------------------------
c precip. tracer evap
        CALL GET_EVAP_FACTOR(N,TL(L),LHP(L),.FALSE.,1d0,FER,FERT)
c clw tracer evap
        CALL GET_EVAP_FACTOR(N,TL(L),LHX,.FALSE.,1d0,FWTOQ,FWTOQT)
        IF(BELOW_CLOUD) THEN
          precip_mm = PREBAR(L+1)*100.*dtsrc
          CALL GET_WASH_FACTOR(N,b_beta_DT,precip_mm,FWASHT) !washout
          WMXTR = PREBAR(L+1)*grav*BYAM(L)*dtsrc
        ELSE
          WMXTR = WMX(L)
c         b_beta_DT is needed at the lowest precipitating level,
c         so saving it here for below cloud case:
          b_beta_DT = FCLD*CM*dtsrc
        END IF
        CALL GET_PREC_FACTOR(N,BELOW_CLOUD,CM,FCLD,FPR,FPRT) !precip CLW
        CALL GET_COND_FACTOR(L,N,WMXTR,TL(L),LHX,FCLD,FQTOW
     *       ,FQTOWT)  !condens
c ---------------------- calculate fluxes ------------------------
        DTWRT = FWASHT*TM(L,N)
        DTERT = FERT  *TRPRBAR(N,L+1)
        DTPRT = FPRT  *TRWML(N,L)
        DTQWT =
     &   FQTOWT*(TM(L,N)+DTERT) - FWTOQT*TRWML(N,L)*(1.-FPRT)
c ---------------------- apply fluxes ------------------------
        TRWML(N,L) = TRWML(N,L)*(1.-FPRT)                   + DTQWT
        TM(L,N)    = TM(L,N)                + DTERT - DTWRT - DTQWT
        TRPRBAR(N,L)=TRPRBAR(N,L+1)*(1.-FERT)+DTPRT + DTWRT
        IF (PREBAR(L).eq.0) TRPRBAR(N,L)=0.  ! remove round off error
        IF (WMX(L).eq.0) TRWML(N,L)=0.       ! remove round off error
c
        TMOM(:,L,N)  = TMOM(:,L,N)*(1. - FQTOWT - FWASHT)
#ifdef TRACERS_SPECIAL_O18
C**** Isotopic equilibration of the CLW and water vapour
        IF (LHX.eq.LHE .and. WMX(L).gt.0) THEN  ! only if liquid
          CALL ISOEQUIL(NTIX(N),TL(L),.TRUE.,QL(L),WMX(L),TM(L,N)
     *         ,TRWML(N,L),1d0)
        END IF
C**** Isotopic equilibration of Precip (if liquid) and water vapour
C**** Note that precip is either all water or all ice
        IF (LHP(L).eq.LHE .AND. PREBAR(L).gt.0 .AND. QL(L).gt.0) THEN
          PRLIQ=PREBAR(L)*DTSrc*BYAM(L)*GRAV
          CALL ISOEQUIL(NTIX(N),TL(L),.TRUE.,QL(L),PRLIQ,TM(L,N)
     *         ,TRPRBAR(N,L),1d0)
        END IF
#endif
      END DO
#endif
C**** CONDENSE MORE MOISTURE IF RELATIVE HUMIDITY .GT. 1
      IF(RH(L).GT.1.) THEN
      SLH=LHX*BYSHA
      DQSUM=0.
      DO N=1,3
        IF(N.NE.1) QSATC=QSAT(TL(L),LHX,PL(L))
        DQ=(QL(L)-QSATC)/(1.+SLH*QSATC*DQSATDT(TL(L),LHX))
        TL(L)=TL(L)+SLH*DQ
        QL(L)=QL(L)-DQ
        DQSUM=DQSUM+DQ
      END DO
      IF(DQSUM.GT.0.) THEN
      WMX(L)=WMX(L)+DQSUM
      FCOND=DQSUM/QNEW
C**** adjust gradients down if Q decreases
      QMOM(:,L)= QMOM(:,L)*(1.-FCOND)
#ifdef TRACERS_WATER
C**** CONDENSING MORE TRACERS
      IF(BELOW_CLOUD) THEN
        WMXTR = PREBAR(L+1)*grav*BYAM(L)*dtsrc
      ELSE
        WMXTR = WMX(L)
      END IF
      DO N=1,NTX
        CALL GET_COND_FACTOR(L,N,WMXTR,TL(L),LHX,FCLD,FCOND
     *       ,FQCONDT)
        TRWML(N,L)  =TRWML(N,L)+ FQCONDT*TM(L,N)
        TM(L,N)     =TM(L,N)    *(1.-FQCONDT)
        TMOM(:,L,N) =TMOM(:,L,N)*(1.-FQCONDT)
      END DO
#endif
      ELSE
      TL(L)=TNEW
      QL(L)=QNEW
      END IF
      RH(L)=QL(L)/QSAT(TL(L),LHX,PL(L))
      TH(L)=TL(L)/PLK(L)
      TNEW=TL(L)
      END IF
      IF(RH(L).LE.1.) CAREA(L)=DSQRT((1.-RH(L))/(1.-RH00(L)+teeny))
      IF(CAREA(L).GT.1.) CAREA(L)=1.
      IF(RH(L).GT.1.) CAREA(L)=0.
      IF(WMX(L).LE.0.) CAREA(L)=1.
      IF(CAREA(L).LT.0.) CAREA(L)=0.
cSI2K      RHF(L)=RH00(L)+(1.-CAREA(L))*(1.-RH00(L))
cSI2K      IF(RH(L).LE.RHF(L).AND.RH(L).LT..999999.AND.WMX(L).gt.0.) THEN
      IF(RH(L).LE.RHF(L)) THEN
C**** PRECIP OUT CLOUD WATER IF RH LESS THAN THE RH OF THE ENVIRONMENT
#ifdef TRACERS_WATER
        TRPRBAR(1:NTX,L) = TRPRBAR(1:NTX,L) + TRWML(1:NTX,L)
        TRWML(1:NTX,L) = 0.
#endif
        PREBAR(L)=PREBAR(L)+WMX(L)*AIRM(L)*BYGRAV*BYDTsrc
        IF(LHP(L).EQ.LHS .AND. LHX.EQ.LHE) THEN
          HCHANG=WMX(L)*LHM
          TL(L)=TL(L)+HCHANG/SHA
          TH(L)=TL(L)/PLK(L)
        END IF
        WMX(L)=0.
      END IF
C**** set phase of condensation for next box down
      PREICE(L)=0.
      IF (PREBAR(L).gt.0 .AND. LHP(L).EQ.LHS) PREICE(L)=PREBAR(L)
      IF (PREBAR(L).le.0) LHP(L)=0.
C**** COMPUTE THE LARGE-SCALE CLOUD COVER
      IF(RH(L).LE.1.) CAREA(L)=DSQRT((1.-RH(L))/(1.-RH00(L)+teeny))
      IF(CAREA(L).GT.1.) CAREA(L)=1.
      IF(RH(L).GT.1.) CAREA(L)=0.
      IF(WMX(L).LE.0.) CAREA(L)=1.
      IF(CAREA(L).LT.0.) CAREA(L)=0.
      CLDSSL(L)=1.-CAREA(L)
#ifdef TRACERS_WATER
      IF(CLDSSL(L).gt.0.) CLOUD_YET=.true.
      IF(CLOUD_YET.and.CLDSSL(L).eq.0.) BELOW_CLOUD=.true.
#endif
      TOLDUP=TOLD
C**** ACCUMULATE SOME DIAGNOSTICS
         HCNDSS=HCNDSS+(TNEW-TOLD)*AIRM(L)
         SSHR(L)=SSHR(L)+(TNEW-TOLD)*AIRM(L)
      END DO  ! end of loop over L

      PRCPSS=MAX(0d0,PREBAR(1)*GRAV*DTsrc) ! fix small round off err
#ifdef TRACERS_WATER
      TRPRSS(1:NTX)=MAX(0d0,TRPRBAR(1:NTX,1))
#endif
C****
C**** CLOUD-TOP ENTRAINMENT INSTABILITY
C****
      DO L=LP50-1,1,-1
        SM(L)=TH(L)*AIRM(L)*PLK(L)
        QM(L)=QL(L)*AIRM(L)
        WMXM(L)=WMX(L)*AIRM(L)
        SM(L+1)=TH(L+1)*AIRM(L+1)*PLK(L+1)
        QM(L+1)=QL(L+1)*AIRM(L+1)
        WMXM(L+1)=WMX(L+1)*AIRM(L+1)
        TOLD=TL(L)
        TOLDU=TL(L+1)
        QOLD=QL(L)
        QOLDU=QL(L+1)
        FCLD=1.-CAREA(L)+1d-30
        IF(CAREA(L).EQ.1. .OR. (CAREA(L).LT.1..AND.CAREA(L+1).LT.1.))
     *       CYCLE
        SEDGE=THBAR(TH(L+1),TH(L))
        DSE=(TH(L+1)-SEDGE)*PLK(L+1)+(SEDGE-TH(L))*PLK(L)+
     *       SLHE*(QL(L+1)-QL(L))
        DWM=QL(L+1)-QL(L)+(WMX(L+1)-WMX(L))/FCLD
        DQSDT=DQSATDT(TL(L),LHE)*QL(L)/(RH(L)+1d-30)
        BETA=(1.+BYMRAT*TL(L)*DQSDT)/(1.+SLHE*DQSDT)
        CKM=(1.+SLHE*DQSDT)*(1.+.392d0*TL(L)/SLHE)/
     *       (2.+(1.+BYMRAT*TL(L)/SLHE)*SLHE*DQSDT)
        CKR=TL(L)/(BETA*SLHE)
        CK=DSE/(SLHE*DWM)
        SIGK=0.
        IF(CKR.GT.CKM) CYCLE
        IF(CK.GT.CKR) SIGK=COESIG*((CK-CKR)/(CKM-CKR+teeny))**5
        EXPST=EXP(-SIGK*DTsrc)
        IF(L.LE.1) CKIJ=EXPST
        DSEC=DWM*TL(L)/BETA
        IF(CK.LT.CKR) CYCLE
        FPMAX=MIN(1d0,1.-EXPST)
        IF(FPMAX.LE.0.) CYCLE
        IF(DSE.GE.DSEC) CYCLE
C**** MIXING TO REMOVE CLOUD-TOP ENTRAINMENT INSTABILITY
        AIRMR=(AIRM(L+1)+AIRM(L))*BYAM(L+1)*BYAM(L)
        SMO1=SM(L)
        QMO1=QM(L)
        WMO1=WMXM(L)
        SMO2=SM(L+1)
        QMO2=QM(L+1)
        WMO2=WMXM(L+1)
        DO K=1,KMAX
          UMO1(K)=UM(K,L)
          VMO1(K)=VM(K,L)
          UMO2(K)=UM(K,L+1)
          VMO2(K)=VM(K,L+1)
        ENDDO
        FPLUME=FPMAX
        DFX=FPMAX
        DO ITER=1,9
          DFX=DFX*0.5
          FMIX=FPLUME*FCLD
          FMASS=FMIX*AIRM(L)
          FMASS=MIN(FMASS,(AIRM(L+1)*AIRM(L))/(AIRM(L+1)+AIRM(L)))
          FMIX=FMASS*BYAM(L)
          FRAT=FMASS*BYAM(L+1)
          SMN1=SMO1*(1.-FMIX)+FRAT*SMO2
          QMN1=QMO1*(1.-FMIX)+FRAT*QMO2
          WMN1=WMO1*(1.-FMIX)+FRAT*WMO2
          SMN2=SMO2*(1.-FRAT)+FMIX*SMO1
          QMN2=QMO2*(1.-FRAT)+FMIX*QMO1
          WMN2=WMO2*(1.-FRAT)+FMIX*WMO1
          THT1=SMN1*BYAM(L)/PLK(L)
          QLT1=QMN1*BYAM(L)
          TLT1=THT1*PLK(L)
          LHX=SVLHXL(L)
          RHT1=QLT1/(QSAT(TLT1,LHX,PL(L)))
          WMT1=WMN1*BYAM(L)
          THT2=SMN2*BYAM(L+1)/PLK(L+1)
          QLT2=QMN2*BYAM(L+1)
          WMT2=WMN2*BYAM(L+1)
          SEDGE=THBAR(THT2,THT1)
          DSE=(THT2-SEDGE)*PLK(L+1)+(SEDGE-THT1)*PLK(L)+SLHE*(QLT2-QLT1)
          DWM=QLT2-QLT1+(WMT2-WMT1)/FCLD
          DQSDT=DQSATDT(TLT1,LHE)*QLT1/(RHT1+1d-30)
          BETA=(1.+BYMRAT*TLT1*DQSDT)/(1.+SLHE*DQSDT)
          CKM=(1.+SLHE*DQSDT)*(1.+.392d0*TLT1/SLHE)/
     *         (2.+(1.+BYMRAT*TLT1/SLHE)*SLHE*DQSDT)
          DSEC=DWM*TLT1/BETA
          DSEDIF=DSE-DSEC
          IF(DSEDIF.GT.1d-3) FPLUME=FPLUME-DFX
          IF(DSEDIF.LT.-1d-3) FPLUME=FPLUME+DFX
          IF(ABS(DSEDIF).LE.1d-3.OR.FPLUME.GT.FPMAX) EXIT
        END DO
C**** UPDATE TEMPERATURE, SPECIFIC HUMIDITY AND MOMENTUM DUE TO CTEI
        TH(L)=SMN1*BYAM(L)/PLK(L)
        TL(L)=TH(L)*PLK(L)
        QL(L)=QMN1*BYAM(L)
        LHX=SVLHXL(L)
        RH(L)=QL(L)/QSAT(TL(L),LHX,PL(L))
        WMX(L)=WMN1*BYAM(L)
        TH(L+1)=SMN2*BYAM(L+1)/PLK(L+1)
        QL(L+1)=QMN2*BYAM(L+1)
        WMX(L+1)=WMN2*BYAM(L+1)
C**** need to scale SM moments for conservation purposes
        SMOM(:,L)=SMOM(:,L)*PLK(L)
        SMOM(:,L+1)=SMOM(:,L+1)*PLK(L+1)
        CALL CTMIX (SM(L),SMOM(1,L),FMASS*AIRMR,FMIX,FRAT)
        SMOM(:,L)=SMOM(:,L)/PLK(L)
        SMOM(:,L+1)=SMOM(:,L+1)/PLK(L+1)
C****
        CALL CTMIX (QM(L),QMOM(1,L),FMASS*AIRMR,FMIX,FRAT)
#ifdef TRACERS_ON
        DO N=1,NTX
          CALL CTMIX (TM(L,N),TMOM(1,L,N),FMASS*AIRMR,FMIX,FRAT)
#ifdef TRACERS_WATER
C**** mix cloud liquid water tracers as well
          TWMTMP      = TRWML(N,L  )*(1.-FMIX)+FRAT*TRWML(N,L+1)
          TRWML(N,L+1)= TRWML(N,L+1)*(1.-FRAT)+FMIX*TRWML(N,L  )
          TRWML(N,L)  = TWMTMP
#endif
        END DO
#endif
        DO K=1,KMAX
          UMN1(K)=(UMO1(K)*(1.-FMIX)+FRAT*UMO2(K))
          VMN1(K)=(VMO1(K)*(1.-FMIX)+FRAT*VMO2(K))
          UMN2(K)=(UMO2(K)*(1.-FRAT)+FMIX*UMO1(K))
          VMN2(K)=(VMO2(K)*(1.-FRAT)+FMIX*VMO1(K))
          UM(K,L)=UM(K,L)+(UMN1(K)-UMO1(K))*RA(K)
          VM(K,L)=VM(K,L)+(VMN1(K)-VMO1(K))*RA(K)
          UM(K,L+1)=UM(K,L+1)+(UMN2(K)-UMO2(K))*RA(K)
          VM(K,L+1)=VM(K,L+1)+(VMN2(K)-VMO2(K))*RA(K)
        END DO
C**** RE-EVAPORATION OF CLW IN THE UPPER LAYER
        QL(L+1)=QL(L+1)+WMX(L+1)
        TH(L+1)=TH(L+1)-(LHX*BYSHA)*WMX(L+1)/PLK(L+1)
        TL(L+1)=TH(L+1)*PLK(L+1)
        RH(L+1)=QL(L+1)/QSAT(TL(L+1),LHX,PL(L+1))
        WMX(L+1)=0.
#ifdef TRACERS_WATER
        TM(L+1,1:NTX)=TM(L+1,1:NTX)+TRWML(1:NTX,L+1)
        TRWML(1:NTX,L+1)=0.
#endif
        IF(RH(L).LE.1.) CAREA(L)=DSQRT((1.-RH(L))/(1.-RH00(L)+teeny))
        IF(CAREA(L).GT.1.) CAREA(L)=1.
        IF(RH(L).GT.1.) CAREA(L)=0.
        CLDSSL(L)=1.-CAREA(L)
        TNEW=TL(L)
        TNEWU=TL(L+1)
        QNEW=QL(L)
        QNEWU=QL(L+1)
        HCNDSS=HCNDSS+(TNEW-TOLD)*AIRM(L)+(TNEWU-TOLDU)*AIRM(L+1)
        SSHR(L)=SSHR(L)+(TNEW-TOLD)*AIRM(L)
        SSHR(L+1)=SSHR(L+1)+(TNEWU-TOLDU)*AIRM(L+1)
        DCTEI(L)=DCTEI(L)+(QNEW-QOLD)*AIRM(L)*LHX*BYSHA
        DCTEI(L+1)=DCTEI(L+1)+(QNEWU-QOLDU)*AIRM(L+1)*LHX*BYSHA
      END DO

C**** COMPUTE CLOUD PARTICLE SIZE AND OPTICAL THICKNESS
      WMSUM=0.
      DO L=1,LP50
        FCLD=CLDSSL(L)+teeny
        WTEM=1.d5*WMX(L)*PL(L)/(FCLD*TL(L)*RGAS+teeny)
        LHX=SVLHXL(L)
        IF(LHX.EQ.LHS.AND.WMX(L)/FCLD.GE.WMUI*1d-3)
     *       WTEM=1d5*WMUI*1.d-3*PL(L)/(TL(L)*RGAS)
        IF(WTEM.LT.1d-10) WTEM=1.d-10
        IF(LHX.EQ.LHE) THEN
          RCLD=(10.*(1.-PEARTH)+7.0*PEARTH)*(WTEM*4.)**BY3
          QHEATC=(QHEAT(L)+CAREA(L)*(EC(L)+ER(L)))/LHX
          IF(RCLD.GT.20..AND.PREP(L).GT.QHEATC) RCLD=20.
        ELSE
          RCLD=25.0*(WTEM/4.2d-3)**BY3 ! * 500./(pl(l)+500.)
        ENDIF
        RCLDE=RCLD/BYBR
        CSIZEL(L)=RCLDE
        TEM=AIRM(L)*WMX(L)*1.d2*BYGRAV
        TAUSSL(L)=1.5d3*TEM/(FCLD*RCLDE+teeny)
        IF(TAUSSL(L).GT.100.) TAUSSL(L)=100.
        IF(LHX.EQ.LHE) WMSUM=WMSUM+TEM
      END DO

C**** CALCULATE OPTICAL THICKNESS
      DO L=1,LP50
        CLDSAVL(L)=CLDSSL(L)
        IF(WMX(L).LE.0.) SVLHXL(L)=0.
        IF(TAUMCL(L).LE.0..OR.CKIJ.NE.1.) THEN
          BMAX=1.-EXP(-(CLDSAVL(L)/.3d0))
          IF(CLDSAVL(L).GE..95d0) BMAX=CLDSAVL(L)
          IF(L.EQ.1) THEN
            CLDSSL(L)=CLDSSL(L)+(BMAX-CLDSSL(L))*CKIJ
            TAUSSL(L)=TAUSSL(L)*CLDSAVL(L)/(CLDSSL(L)+teeny)
          ENDIF
          IF(TAUSSL(L).LE.0.) CLDSSL(L)=0.
          IF(L.gt.1.and.TAUMCL(L).LE.0.) THEN
            CLDSSL(L)=CLDSSL(L)**(2.*BY3)
            TAUSSL(L)=TAUSSL(L)*CLDSAVL(L)**BY3
          END IF
        END IF
      END DO

      RETURN
      END SUBROUTINE LSCOND

      END MODULE CLOUDS


C-----------------------------------------------------------------------

      SUBROUTINE ISCCP_CLOUD_TYPES(pfull,phalf,qv,cc,conv,dtau_s,dtau_c
     *  ,skt,at,dem_s,dem_c,itrop,fq_isccp,ctp,tauopt,nbox,jerr)
!@sum  ISCCP_CLOUD_TYPES calculate isccp cloud diagnostics in a column
!@auth G. Tselioudis (modifications by Gavin Schmidt)
!@ver  1.0
      USE CONSTANT, only : bygrav, wtmair=>mair, bymrat
      USE RANDOM, only : rinit,rfinal,randu
      USE MODEL_COM, only : nlev=>lm,qcheck
      implicit none
!@var  overlap type: 1=max, 2=rand,  3=max/rand
!@var  top_height 1 = adjust top height, that is compute infrared
!@+     brightness temperature and adjust cloud top pressure accordingly
!@+               2 = do not adjust top height, that is cloud top
!@+     pressure is the actual cloud top pressure in the model
      INTEGER, PARAMETER :: top_height=1, overlap=3
!@var  emsfc_lw    longwave emissivity of surface at 10.5 microns
      REAL*8, PARAMETER :: emsfc_lw=0.99d0
      INTEGER,PARAMETER :: ncol =20    !@var ncol number of subcolumns
      REAL*8, PARAMETER :: byncol = 1d0/ncol
      REAL*8, PARAMETER :: Navo = 6.023d23 !@var Navo Avogardos Numbers
!@var pc1bylam Planck constant c1 by wavelength (10.5 microns)
      REAL*8, PARAMETER :: pc1bylam = 1.439d0/10.5d-4
!@var t0 ave temp  (K)
      REAL*8, PARAMETER :: t0 = 296.
!@var t0bypstd ave temp by sea level press
      REAL*8, PARAMETER :: t0bypstd = t0/1.013250d6
      real*8, parameter :: bywc = 1./2.56d0 , byic= 1./2.13d0
!     -----
!     Input
!     -----
!@var pfull pressure of full model levels (Pascals)
      REAL*8 pfull(nlev)  ! pfull(1) is top level, pfull(nlev) is bottom
!@var phalf pressure of half model levels (Pascals)
      REAL*8 phalf(nlev+1) ! phalf(1) is top, phalf(nlev+1) is surf pres
!@var qv  water vapor specific humidity (kg vapor/ kg air)
      REAL*8 qv(nlev)
!@var cc  input cloud cover in each model level (fraction)
      REAL*8 cc(nlev)           ! NOTE: This is the HORIZONTAL area of
                                ! each grid box covered by clouds
!@var conv input convective cloud cover in each model level(frac)
      REAL*8 conv(nlev)         ! NOTE: This is the HORIZONTAL area of
                                ! each box covered by convective clouds
!@var dtau_s mean 0.67 micron optical depth of stratiform cloud level
!@var dtau_c mean 0.67 micron optical depth of convective cloud level
      REAL*8 dtau_s(nlev), dtau_c(nlev)
           !  NOTE:  these cloud optical depths are only for the
           !         cloudy part of the grid box, they are not weighted
           !         with the 0 cloud optical depth of the clear
           !         part of the grid box

      INTEGER :: itrop    !@var itrop tropopause level (WMO definition)

!     The following input variables are used only if top_height = 1
      REAL*8 skt            !@var skt skin Temperature (K)
      REAL*8 at(nlev)       !@var at temperature in each model level (K)
!@var dem_s 10.5 micron longwave emissivity of stratiform clouds
!@var dem_c 10.5 micron longwave emissivity of convective clouds
      REAL*8 dem_s(nlev),dem_c(nlev)
!     ------
!     Output
!     ------
!@var fq_isccp the fraction of the model grid box covered by each
!@+   of the 49 ISCCP D level cloud types
      REAL*8 fq_isccp(7,7)
!@var ctp cloud top pressure averaged over grid box (mb)
      REAL*8 ctp
!@var tauopt cloud optical thickness averaged over grid box
      REAL*8 tauopt
!@var nbox number of boxes with clouds (used as a flag)
      INTEGER nbox
!     ------
! Working variables
!     ------
!@var frac_out boxes gridbox divided up into
      REAL*8 frac_out(ncol,nlev) ! Equivalent of BOX in orig. version,
                                 ! but indexed by column then row,
                                 ! rather than by row then column
!@var tca total cloud cover (fraction) with extra layer of zeros on top
      REAL*8 tca(ncol,0:nlev)  ! in this version this just contains the
                      ! values input from cc but replicated across ncol
!@var cca convect. cloud cover each model level(frac)
      REAL*8 cca(ncol,nlev)  ! from conv but replicated across ncol
!@var threshold pointer to position in gridbox
!@var maxocc Flag for max overlapped conv cld
!@var maxpsc Flag for max overlapped strat cld
!@var boxpos ordered pointer to position in gridbox
!@var threshold_min min. allowed value of threshold
      REAL*8, dimension(ncol) :: threshold,maxocc,maxosc,boxpos,
     *     threshold_min
!@var dem,bb,bbs working variables for 10.5 micron longwave emissivity
!@+ in part of gridbox under consideration
      real*8 dem(ncol),bb(nlev),bbs
      integer seed  !@var seed saved value for random number generator
!@var dtautmp,demtmp temporary variables for dtau,dem
      real*8, dimension(ncol) :: dtautmp, demtmp
      real*8 ptrop,attrop,atmax,atmin,btcmin,transmax
      integer ilev,ibox,ipres,itau,ilev2
      integer acc(nlev,ncol),match(nlev-1),nmatch,levmatch(ncol)

      !variables needed for water vapor continuum absorption
      real*8 fluxtop_clrsky,trans_layers_above_clrsky,taumin
      real*8 dem_wv(nlev)
      real*8 press, dpress, atmden, rvh20, wk, rhoave, rh20s, rfrgn
      real*8 tmpexp,tauwv, XX

      real*8, dimension(ncol) :: tau,tb,ptop,emcld,fluxtop,
     *     trans_layers_above
      real*8, parameter :: isccp_taumin = 0.1d0

      INTEGER  JERR

      character*1 cchar(6),cchar_realtops(6)
      data cchar / ' ','-','1','+','I','+'/
      data cchar_realtops / ' ',' ','1','1','I','I'/

      JERR=0

!     assign 2d tca array using 1d input array cc
      do ilev=0,nlev
        do ibox=1,ncol
          if (ilev.eq.0) then
            tca(ibox,ilev)=0
          else
            tca(ibox,ilev)=cc(ilev)
          endif
        enddo
      enddo
!     assign 2d cca array using 1d input array conv
      do ilev=1,nlev
        do ibox=1,ncol
          cca(ibox,ilev)=conv(ilev)
        enddo
      enddo

C**** In order to ensure that the model does not go down a different
C**** path depending on whether this routine is used, we save current
C**** seed and reset it afterwards
cc    CALL RFINAL(SEED)

      if (top_height .eq. 1) then
        ptrop = pfull(itrop)
        attrop = at(itrop)
        atmin = 400.
        atmax = 0.
        do ilev=1,nlev-1
          if (at(ilev) .gt. atmax) atmax=at(ilev)
          if (at(ilev) .lt. atmin) atmin=at(ilev)
        end do
      end if

!     ---------------------------------------------------!
!     find unpermittable data.....
!
      do ilev=1,nlev
        if (cc(ilev) .lt. 0.) then
          print*, ' error = cloud fraction less than zero'
          JERR=1 ; return; ! stop
        end if
        if (cc(ilev) .gt. 1.) then
          print*,' error = cloud fraction greater than 1'
          JERR=1 ; return; ! stop
        end if
        if (conv(ilev) .lt. 0.) then
          print*,' error = convective cloud fraction less than zero'
          JERR=1 ; return; ! stop
        end if
        if (conv(ilev) .gt. 1.) then
          print*,' error = convective cloud fraction greater than 1'
          JERR=1 ; return; ! stop
        end if

        if (dtau_s(ilev) .lt. 0.) then
          print*,' error = stratiform cloud opt. depth less than zero'
          JERR=1 ; return; ! stop
        end if
        if (dem_s(ilev) .lt. 0.) then
          print*,' error = stratiform cloud emissivity less than zero'
          JERR=1 ; return; ! stop
        end if
        if (dem_s(ilev) .gt. 1.) then
          print*,' error = stratiform cloud emissivity greater than 1'
          JERR=1 ; return; ! stop
        end if

        if (dtau_c(ilev) .lt. 0.) then
          print*, ' error = convective cloud opt. depth less than zero'
          JERR=1 ; return; ! stop
        end if
        if (dem_c(ilev) .lt. 0.) then
          print*,' error = convective cloud emissivity less than zero'
          JERR=1 ; return; ! stop
        end if
        if (dem_c(ilev) .gt. 1.) then
          print*,' error = convective cloud emissivity greater than 1'
          JERR=1 ; return; ! stop
        end if
      end do

!     ---------------------------------------------------!
!     Initialise working variables
!     ---------------------------------------------------!

!     Initialised frac_out to zero
      frac_out(:,:)=0.0
      do ibox=1,ncol
        boxpos(ibox)=(ibox-.5)/ncol
      enddo

!     ---------------------------------------------------!
!     ALLOCATE CLOUD INTO BOXES, FOR NCOLUMNS, NLEVELS
!     frac_out is the array that contains the information
!     where 0 is no cloud, 1 is a stratiform cloud and 2 is a
!     convective cloud

      !loop over vertical levels
      do ilev = 1,nlev

!     Initialise threshold

        if (ilev.eq.1) then
          do ibox=1,ncol
            ! if max overlap
            if (overlap.eq.1) then
              ! select pixels spread evenly
              ! across the gridbox
              threshold(ibox)=boxpos(ibox)
            else
              ! select random pixels from the non-convective
              ! part the gridbox ( some will be converted into
              ! convective pixels below )
              threshold(ibox)=
     *        cca(ibox,ilev)+(1-cca(ibox,ilev))*randu(xx)
            endif
          enddo
        endif

        do ibox=1,ncol
          ! All versions
          if (boxpos(ibox).le.cca(ibox,ilev)) then
            maxocc(ibox) = 1
          else
            maxocc(ibox) = 0
          end if

          select case (overlap)
          case (1)              ! Max overlap
            threshold_min(ibox)=cca(ibox,ilev)
            maxosc(ibox)=1
          case (2)              ! Random overlap
            threshold_min(ibox)=cca(ibox,ilev)
            maxosc(ibox)=0
          case (3)              ! Max/Random overlap
            threshold_min(ibox)=max(cca(ibox,ilev),
     &           min(tca(ibox,ilev-1),tca(ibox,ilev)))
            if (threshold(ibox).lt.min(tca(ibox,ilev-1),tca(ibox,ilev))
     &           .and.(threshold(ibox).gt.cca(ibox,ilev))) then
              maxosc(ibox)= 1
            else
              maxosc(ibox)= 0
            end if
          end select
          ! Reset threshold
          threshold(ibox)=
                                !if max overlapped conv cloud
     &         maxocc(ibox) * boxpos(ibox) +
                                !else
     &         (1-maxocc(ibox)) * (
                                !if max overlapped strat cloud
     &         (maxosc(ibox)) * (
                                !threshold=boxpos
     &         threshold(ibox) ) +
                                !else
     &         (1-maxosc(ibox)) * (
                                !threshold_min=random[thrmin,1]
     &         threshold_min(ibox)+(1-threshold_min(ibox))*RANDU(XX)))
        end do

!       Fill frac_out with 1's where tca is greater than the threshold

        do ibox=1,ncol
          if (tca(ibox,ilev).gt.threshold(ibox)) then
            frac_out(ibox,ilev)=1
          else
            frac_out(ibox,ilev)=0
          end if
        end do

!    Code to partition boxes into stratiform and convective parts
!    goes here

        do ibox=1,ncol
          if (threshold(ibox).le.cca(ibox,ilev)) then
                              ! = 2 IF threshold le cca(ibox)
            frac_out(ibox,ilev) = 2
          else
                              ! = the same IF NOT threshold le cca(ibox)
            frac_out(ibox,ilev) = frac_out(ibox,ilev)
          end if
        end do
      end do                  !loop over nlev
!
!     ---------------------------------------------------!
!     COMPUTE CLOUD OPTICAL DEPTH FOR EACH COLUMN and
!     put into vector tau

      !initialize tau to zero
      tau(:)=0.

      !compute total cloud optical depth for each column
      do ilev=1,nlev
            !increment tau for each of the boxes
        do ibox=1,ncol
          if (frac_out(ibox,ilev).eq.1) then
            dtautmp(ibox)= dtau_s(ilev)
          else if (frac_out(ibox,ilev).eq.2) then
            dtautmp(ibox)= dtau_c(ilev)
          else
            dtautmp(ibox)= 0.
          end if

          tau(ibox)=tau(ibox)+dtautmp(ibox)
        end do
      end do
!     ---------------------------------------------------!
!     COMPUTE INFRARED BRIGHTNESS TEMPERATURES
!     AND CLOUD TOP TEMPERATURE SATELLITE SHOULD SEE
!
!     again this is only done if top_height = 1
!
!     fluxtop is the 10.5 micron radiance at the top of the
!              atmosphere
!     trans_layers_above is the total transmissivity in the layers
!             above the current layer
!     fluxtop_clrsky and trans_layers_above_clrsky are the clear
!             sky versions of these quantities.
      if (top_height .eq. 1) then

        !---------------------------------------------------------------
        !
        !             DO CLEAR SKY RADIANCE CALCULATION FIRST
        !
        !compute water vapor continuum emissivity
        !this treatment follows Schwarkzopf and Ramasamy
        !JGR 1999,vol 104, pages 9467-9499.
        !the emissivity is calculated at a wavenumber of 955 cm-1,
        !or 10.47 microns
        do ilev=1,nlev
                            !press and dpress are dyne/cm2 = Pascals *10
          press = pfull(ilev)*10.
          dpress = (phalf(ilev+1)-phalf(ilev))*10.
                                !atmden = g/cm2 = kg/m2 / 10
          atmden = dpress*bygrav
          rvh20 = qv(ilev)*bymrat    !wtmair/wtmh20
          wk = rvh20*Navo*atmden/wtmair
c          rhoave = (press/pstd)*(t0/at(ilev))
          rhoave = (press/at(ilev))*t0bypstd
          rh20s = rvh20*rhoave
          rfrgn = rhoave-rh20s
          tmpexp = exp(-0.02d0*(at(ilev)-t0))
          tauwv = wk*1d-20*( (0.0224697d0*rh20s*tmpexp) +
     &         (3.41817d-7*rfrgn)         )*0.98d0
          dem_wv(ilev) = 1. - exp(-tauwv)
        end do

        !initialize variables
        fluxtop_clrsky = 0.
        trans_layers_above_clrsky=1.

        do ilev=1,nlev
            ! Black body emission at temperature of the layer
          bb(ilev)=1 / ( exp(pc1bylam/at(ilev)) - 1. )

                          ! increase TOA flux by flux emitted from layer
                          !    times total transmittance in layers above
          fluxtop_clrsky = fluxtop_clrsky
     &         + dem_wv(ilev) * bb(ilev) * trans_layers_above_clrsky

                         ! update trans_layers_above with transmissivity
                         !     from this layer for next time around loop

          trans_layers_above_clrsky=
     &         trans_layers_above_clrsky*(1.-dem_wv(ilev))
        end do                  !loop over level

        !add in surface emission
        bbs=1/( exp(pc1bylam/skt) - 1. )

        fluxtop_clrsky = fluxtop_clrsky + emsfc_lw * bbs
     &       * trans_layers_above_clrsky
!           END OF CLEAR SKY CALCULATION
!
!----------------------------------------------------------------
        fluxtop(:)=0.
        trans_layers_above(:)=1.

        do ilev=1,nlev
          do ibox=1,ncol
                                ! emissivity for point in this layer
            if (frac_out(ibox,ilev).eq.1) then
              dem(ibox)= 1. -
     &             ( (1. - dem_wv(ilev)) * (1. -  dem_s(ilev)) )
            else if (frac_out(ibox,ilev).eq.2) then
              dem(ibox)= 1. -
     &             ( (1. - dem_wv(ilev)) * (1. -  dem_c(ilev)) )
            else
              dem(ibox)=  dem_wv(ilev)
            end if
                ! increase TOA flux by flux emitted from layer
         ! times total transmittance in layers above
            fluxtop(ibox) = fluxtop(ibox)
     &           + dem(ibox) * bb(ilev)
     &           * trans_layers_above(ibox)
                ! update trans_layers_above with transmissivity
         ! from this layer for next time around loop
            trans_layers_above(ibox)=
     &           trans_layers_above(ibox)*(1.-dem(ibox))

          end do                ! ibox
        end do                  ! ilev

                                !add in surface emission
        do ibox=1,ncol
          fluxtop(ibox) = fluxtop(ibox)
     &         + emsfc_lw * bbs * trans_layers_above(ibox)
        end do

        do ibox=1,ncol
            !now that you have the top of atmosphere radiance account
            !for ISCCP procedures to determine cloud top temperature

            !account for partially transmitting cloud recompute flux
            !ISCCP would see assuming a single layer cloud
            !note choice here of 2.13, as it is primarily ice
            !clouds which have partial emissivity and need the
            !adjustment performed in this section
            !
            !Note that this is discussed on pages 85-87 of
            !the ISCCP D level documentation (Rossow et al. 1996)

            !compute minimum brightness temperature and optical depth
          btcmin = 1. /  ( exp(pc1bylam/(attrop-5.)) - 1. )
          transmax = (fluxtop(ibox)-btcmin)/(fluxtop_clrsky-btcmin)
          taumin = -log(max(min(transmax,0.9999999d0),0.001d0))

          if (transmax .gt. 0.001 .and. transmax .le. 0.9999999) then
            emcld(ibox) = 1. - exp(-tau(ibox)*byic)
            fluxtop(ibox) = fluxtop(ibox) -
     &           ((1.-emcld(ibox))*fluxtop_clrsky)
            fluxtop(ibox)=max(1d-6,(fluxtop(ibox)/emcld(ibox)))
          end if

          if (tau(ibox) .gt.  1d-7) then
                !cloudy box
            tb(ibox)= pc1bylam/(log(1. + (1./fluxtop(ibox))))

            if (tau(ibox) .lt. taumin) then
              tb(ibox) = attrop - 5.
            end if
          else
                !clear sky brightness temperature
            tb(ibox) = pc1bylam/(log(1.+(1./fluxtop_clrsky)))
          end if
        enddo                   ! ibox
      end if
!
!     ---------------------------------------------------!
!     DETERMINE CLOUD TOP PRESSURE
!
!     again the 2 methods differ according to whether
!     or not you use the physical cloud top pressure (top_height = 2)
!     or the radiatively determined cloud top pressure (top_height = 1)
!
      !compute cloud top pressure
      do ibox=1,ncol

               !segregate according to optical thickness
        if (tau(ibox) .le. 1d-7) then
          ptop(ibox)=0.
          levmatch(ibox)=0
        else
          if (top_height .eq. 1) then
                               ! find level whose temperature most
                               ! closely matches brightness temperature
            nmatch=0
            do ilev=1,nlev-1
              if ((at(ilev)   .ge. tb(ibox) .and.
     &             at(ilev+1) .lt. tb(ibox)) .or.
     &             (at(ilev) .le. tb(ibox) .and.
     &             at(ilev+1) .gt. tb(ibox))) then

                nmatch=nmatch+1
                if(abs(at(ilev)-tb(ibox)) .lt.
     &               abs(at(ilev+1)-tb(ibox))) then
                  match(nmatch)=ilev
                else
                  match(nmatch)=ilev+1
                end if
              end if
            end do

            if (nmatch .ge. 1) then
              ptop(ibox)=pfull(match(nmatch))
              levmatch(ibox)=match(nmatch)
            else
              if (tb(ibox) .lt. atmin) then
                ptop(ibox)=ptrop
                levmatch(ibox)=itrop
              end if
              if (tb(ibox) .gt. atmax) then
                ptop(ibox)=pfull(nlev)
                levmatch(ibox)=nlev
              end if
            end if
          else
            ptop(ibox)=0.
            ilev=1
            do while(ptop(ibox) .eq. 0.
     &           .and. ilev .lt. nlev+1)
              if (frac_out(ibox,ilev) .ne. 0) then
                ptop(ibox)=pfull(ilev)
                levmatch(ibox)=ilev
              end if
              ilev=ilev+1
            end do
          end if
        end if
      end do
!
!     ---------------------------------------------------!
!     DETERMINE ISCCP CLOUD TYPE FREQUENCIES
!
!     Now that ptop and tau have been determined,
!     determine amount of each of the 49 ISCCP cloud
!     types
!
      !compute isccp frequencies
      fq_isccp(:,:)=0.
      ctp = 0.
      tauopt = 0.
      nbox = 0

      do ibox=1,ncol
            !convert ptop to millibars
        ptop(ibox)=ptop(ibox)*1d-2

        if (tau(ibox) .gt. 1d-7 .and. ptop(ibox) .gt. 0.) then
            !reset itau, ipres
          itau = 0
          ipres = 0
            !determine optical depth category
          if (tau(ibox) .lt. isccp_taumin) then
            itau=1
          else if (tau(ibox) .ge. isccp_taumin
     &           .and. tau(ibox) .lt. 1.3d0) then
            itau=2
          else if (tau(ibox) .ge. 1.3d0 .and. tau(ibox) .lt. 3.6d0) then
            itau=3
          else if (tau(ibox) .ge. 3.6d0 .and. tau(ibox) .lt. 9.4d0) then
            itau=4
          else if (tau(ibox) .ge. 9.4d0 .and. tau(ibox) .lt. 23.) then
            itau=5
          else if (tau(ibox) .ge. 23. .and. tau(ibox) .lt. 60.) then
            itau=6
          else if (tau(ibox) .ge. 60.) then
            itau=7
          end if

            !determine cloud top pressure category
          if (    ptop(ibox) .gt. 0.  .and.ptop(ibox) .lt. 180.) then
            ipres=1
          else if(ptop(ibox) .ge. 180..and.ptop(ibox) .lt. 310.) then
            ipres=2
          else if(ptop(ibox) .ge. 310..and.ptop(ibox) .lt. 440.) then
            ipres=3
          else if(ptop(ibox) .ge. 440..and.ptop(ibox) .lt. 560.) then
            ipres=4
          else if(ptop(ibox) .ge. 560..and.ptop(ibox) .lt. 680.) then
            ipres=5
          else if(ptop(ibox) .ge. 680..and.ptop(ibox) .lt. 800.) then
            ipres=6
          else if(ptop(ibox) .ge. 800.) then
            ipres=7
          end if

            !update frequencies
          if(ipres .gt. 0.and.itau .gt. 0) then
            fq_isccp(itau,ipres)=fq_isccp(itau,ipres)+byncol ! 1.d0/ncol
          end if

C**** accumulate ptop/tauopt over columns for output
          if (itau.gt.1) then
            ctp   = ctp  +ptop(ibox)
            tauopt=tauopt+ tau(ibox)
            nbox = nbox + 1
          end if
        end if
      end do

      if (nbox.gt.0) then
        ctp = ctp/REAL(nbox,KIND=8)
        tauopt=tauopt/REAL(nbox,KIND=8)
      end if
cc    CALL RINIT(SEED)   ! reset seed to original value

!
!     ---------------------------------------------------!
!     OPTIONAL PRINTOUT OF DATA TO CHECK PROGRAM
!
      RETURN   ! this is too much information!
      if (QCHECK) then
        do ilev=1,nlev
          do ibox=1,ncol
            acc(ilev,ibox)=frac_out(ibox,ilev)*2
            if (levmatch(ibox) .eq. ilev)
     &           acc(ilev,ibox)=acc(ilev,ibox)+1
          end do
        end do

        do ilev=1,nlev
          write(6,'(i2,1X,8(f7.2,1X),50(a1))')
     &         ilev,pfull(ilev)/100.,at(ilev),
     &         cc(ilev)*100.0,dem_wv(ilev),dem_s(ilev),dtau_s(ilev)
     *         ,dem_c(ilev),dtau_c(ilev),(cchar(acc(ilev,ibox)+1),ibox=1
     *         ,ncol)
        end do
        print*, 'skt = ',skt
        write (6,'(8I7)') (ibox,ibox=1,ncol)
        write (6,'(a)') 'tau:'
        write (6,'(8f7.2)') (tau(ibox),ibox=1,ncol)
        write (6,'(a)') 'tb:'
        write (6,'(8f7.2)') (tb(ibox),ibox=1,ncol)
        write (6,'(a)') 'ptop:'
        write (6,'(8f7.2)') (ptop(ibox),ibox=1,ncol)
        print*, 'ctp,tauopt,nbox =',ctp,tauopt,nbox
      end if
      return
      end

