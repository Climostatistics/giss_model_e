#include "rundeck_opts.h"

      SUBROUTINE SURFCE
!@sum SURFCE calculates the surface fluxes which include
!@+   sensible heat, evaporation, thermal radiation, and momentum
!@+   drag.  It also calculates instantaneous surface temperature,
!@+   surface specific humidity, and surface wind components.
!@auth Nobody will claim responsibilty
      USE CONSTANT, only : grav,rgas,kapa,sday,lhm,lhe,lhs,twopi
     *     ,sha,tf,rhow,rhoi,shv,shw,shi,rvap,stbo,bygrav,by6,byshi
     *     ,byrhoi,deltx,byrt3,teeny
      USE MODEL_COM, only : im,jm,lm,fim,dtsrc,nisurf,u,v,t,p,q
     *     ,idacc,dsig,jday,ndasf,jeq,fland,flice,focean
     *     ,fearth,nday,modrd,itime,jhour,sige,byim,itocean
     *     ,itoice,itlake,itlkice,itlandi,qcheck,UOdrag
      USE GEOM, only : dxyp,imaxj,bydxyp,idjj,idij,rapj,kmaxj,sinip
     *     ,cosip
      USE SOMTQ_COM, only : tmom,qmom,mz,nmom
      USE DYNAMICS, only : pmid,pk,pedn,pek,pdsig,plij,am,byam
      USE RADNCB, only : trhr,fsf,cosz1
#ifdef TRACERS_ON
      USE TRACER_COM, only : ntm,itime_tr0,needtrs,trm,trmom,ntsurfsrc
#ifdef TRACERS_WATER
     *     ,nWATER,nGAS,nPART,tr_wd_TYPE,trname,trw0
#endif
#endif
C**** Interface to PBL
      USE SOCPBL, only : zgs,ZS1,TGV,TKV,QG_SAT,HEMI,DTSURF,POLE
     &     ,US,VS,WS,WSM,WSH,TSV,QSRF,PSI,DBL,KHS,KQS !,PPBL ! ,KMS
     &     ,UG,VG,WG,ZMIX
      USE PBLCOM, only : ipbl,cmgs,chgs,cqgs,tsavg,dclev
      USE PBL_DRV, only : pbl,evap_max,fr_sat,uocean,vocean
#ifdef TRACERS_ON
     *     ,trtop,trs,trsfac,trconstflx,ntx,ntix
#ifdef TRACERS_WATER
     *     ,tr_evap_max
#endif
#endif
      USE DAGCOM, only : oa,aij,tdiurn,aj,areg,adiurn,ndiupt,jreg
     *     ,ij_tsli,ij_shdtli,ij_evhdt,ij_trhdt,ij_shdt,ij_trnfp0
     *     ,ij_srtr,ij_neth,ij_ws,ij_ts,ij_us,ij_vs,ij_taus,ij_tauus
     *     ,ij_tauvs,ij_qs,ij_tg1,ij_evap,ij_evapo,ij_tgo,ij_f0oc
     *     ,ij_f0oi,ij_evapi,ij_f0li,ij_evapli,j_evap,j_evhdt
     *     ,j_tsrf,j_shdt,j_trhdt,j_type,j_tg1,j_tg2,ijdd,idd_spr
     *     ,idd_pt5,idd_pt4,idd_pt3,idd_pt2,idd_pt1,idd_ts,idd_tg1
     *     ,idd_q5,idd_q4,idd_q3,idd_q2,idd_q1,idd_qs,idd_qg,idd_swg
     *     ,idd_lwg,idd_sh,idd_lh,idd_hz0,idd_ug,idd_vg,idd_wg,idd_us
     *     ,idd_vs,idd_ws,idd_cia,idd_cm,idd_ch,idd_cq,idd_eds,idd_dbl
     *     ,idd_ev,idd_ldc,idd_dcf
      USE LANDICE, only : hc2li,z1e,z2li,hc1li
      USE LANDICE_COM, only : snowli
      USE SEAICE, only : xsi,z1i,ace1i,hc1i,alami,byrli,byrls,
     *     solar_ice_frac,tfrez
      USE SEAICE_COM, only : rsi,msi,snowi,flag_dsws
      USE LAKES_COM, only : mwl,mldlk,gml,flake
      USE LAKES, only : minmld
      USE FLUXES, only : dth1,dq1,e0,e1,evapor,runoe,erunoe,sss
     *     ,solar,dmua,dmva,gtemp,nstype,uflux1,vflux1,tflux1,qflux1
     *     ,uosurf,vosurf,uisurf,visurf
#ifdef TRACERS_ON
     *     ,trsrfflx,trsource
#ifdef TRACERS_WATER
     *     ,trevapor,trunoe,gtracer
#endif
      USE TRACER_DIAG_COM, only : taijn,tij_surf
#ifdef TRACERS_WATER
     *     ,tij_evap,tij_grnd,tajls,jls_source
#endif
#endif
      USE SOIL_DRV, only: earth
      IMPLICIT NONE

      INTEGER I,J,K,KR,JR,NS,NSTEPS,MODDSF,MODDD,ITYPE,IH,IDTYPE,IM1
      REAL*8 PLAND,PLICE,POICE,POCEAN,PIJ,PS,P1K,PGK,PKDN
     *     ,BETA,ELHX,ACE2,CDTERM,CDENOM,dF1dTG,HCG1,HCG2,EVHDT,F1DT
     *     ,CM,CH,CQ,DHGS,DQGS,DGS,BETAUP,EVHEAT,F0,F1,DSHDTG,DQGDTG
     *     ,DEVDTG,DTRDTG,DF0DTG,DFDTG,DTG,dSNdTG,dEVdQS,HSDEN,HSCON
     *     ,HSMUL,dHS,dQS,dT2,dTS,DQ1X,EVHDT0,EVAP,F0DT,FTEVAP,PWATER
     *     ,PSK,Q1,THV1,PTYPE,TG1,SRHEAT,SNOW,TG2
     *     ,SHDT,TRHDT,TG,TS,RHOSRF,RCDMWS,RCDHWS,RCDQWS,SHEAT,TRHEAT
     *     ,QSDEN,QSCON,QSMUL,T2DEN,T2CON,T2MUL,TGDEN,FQEVAP
     *     ,Z1BY6L,QZ1,EVAPLIM,F2,FSRI(2),HTLIM

      REAL*8 MSUM, MA1, MSI1
      REAL*8, DIMENSION(NSTYPE,IM,JM) :: TGRND,TGRN2
      REAL*8, PARAMETER :: qmin=1.d-12
      REAL*8, PARAMETER :: S1BYG1 = BYRT3,
     *     Z1IBYL=Z1I/ALAMI, Z2LI3L=Z2LI/(3.*ALAMI), Z1LIBYL=Z1E/ALAMI
      REAL*8 QSAT,DQSATDT
      REAL*8 AREGIJ(7,3,IM,JM)
c
#ifdef TRACERS_ON
      real*8 rhosrf0, totflux
      integer n,nx,nsrc
#ifdef TRACERS_WATER
      real*8, dimension(ntm) :: tevaplim,trgrnd
      real*8  TEV,dTEVdTQS,tevap,dTQS,TDP,TDT1
#ifdef TRACERS_SPECIAL_O18
     *     ,FRACVL,FRACVS,FRACLK,frac
#endif
#endif
#endif

      NSTEPS=NIsurf*ITime
      DTSURF=DTsrc/NIsurf
      IH=JHOUR+1

C**** ZERO OUT ENERGY AND EVAPORATION FOR GROUND AND INITIALIZE TGRND
      DO J=1,JM
      DO I=1,IM
        TGRND(2,I,J)=GTEMP(1,2,I,J)
        TGRND(3,I,J)=GTEMP(1,3,I,J)
C       TGRND(4,I,J)=GTEMP(1,4,I,J)
        TGRN2(2,I,J)=GTEMP(2,2,I,J)
        TGRN2(3,I,J)=GTEMP(2,3,I,J)
      END DO
      END DO

C**** Zero out fluxes summed over type and surface time step
      E0=0. ; E1=0. ; EVAPOR=0. ; RUNOE=0. ; ERUNOE=0.
      DMUA=0. ; DMVA=0. ; SOLAR=0.
      AREGIJ = 0.
#ifdef TRACERS_WATER
      TREVAPOR = 0. ; TRUNOE = 0.
#endif

C****
C**** OUTSIDE LOOP OVER TIME STEPS, EXECUTED NIsurf TIMES EVERY HOUR
C****
      DO NS=1,NIsurf
         MODDSF=MOD(NSTEPS+NS-1,NDASF*NIsurf+1)
         IF(MODDSF.EQ.0) IDACC(3)=IDACC(3)+1
         MODDD=MOD(1+ITime/NDAY+NS,NIsurf)   ! 1+ not really needed ??
C**** ZERO OUT FLUXES ACCUMULATED OVER SURFACE TYPES
         DTH1=0. ;  DQ1 =0. ;  uflux1=0. ; vflux1=0.
#ifdef TRACERS_WATER
         trsrfflx = 0.
#endif

      call loadbl
C****
C**** OUTSIDE LOOP OVER J AND I, EXECUTED ONCE FOR EACH GRID POINT
C****
!$OMP   PARALLEL DO PRIVATE (ACE2, BETA,BETAUP,   CM,CH,CQ,
!$OMP*  CDTERM,CDENOM, DGS,DSHDTG,DQGDTG,DEVDTG,DTRDTG,
!$OMP*  DF0DTG,DFDTG,DTG,DQ1X,DF1DTG,DHGS,DQGS,DSNDTG,DEVDQS,
!$OMP*  DHS,DQS,DT2,DTS, EVAP,EVAPLIM,ELHX,EVHDT,EVHEAT,EVHDT0,
!$OMP*  F0DT,F1DT,F0,F1,F2,FSRI, HCG1,HCG2,HSDEN,HSCON,
!$OMP*  HSMUL,HTLIM,I,ITYPE,IDTYPE,IM1, J,K,
!$OMP*  KR, MSUM,MA1,MSI1, PS,P1K,PLAND,PWATER,
!$OMP*  PLICE,PIJ,POICE,POCEAN,PGK,PKDN,PTYPE,PSK, Q1,QSDEN,
!$OMP*  QSCON,QSMUL, RHOSRF,RCDMWS,RCDHWS,RCDQWS, SHEAT,SRHEAT,
!$OMP*  SNOW,SHDT, T2DEN,T2CON,T2MUL,TGDEN,TS,
!$OMP*  THV1,TG,TG1,TG2,TRHDT,TRHEAT,Z1BY6L
#ifdef TRACERS_ON
!$OMP*  ,n,nx,nsrc,rhosrf0,totflux
#ifdef TRACERS_WATER
!$OMP*  ,tevaplim,tevap,trgrnd,TEV,dTEVdTQS,dTQS,TDP,TDT1
#endif
#endif
!$OMP*  )
!$OMP*  SCHEDULE(DYNAMIC,2)
C
      DO J=1,JM
      HEMI=1.
      IF(J.LE.JM/2) HEMI=-1.
      POLE= (J.EQ.1 .or. J.EQ.JM)

      IM1=IM
      DO I=1,IMAXJ(J)

      EVAPLIM = 0. ; HTLIM=0.  ! need initialisation
#ifdef TRACERS_WATER
      tevaplim = 0.
#endif
C****
C**** DETERMINE SURFACE CONDITIONS
C****
      PLAND=FLAND(I,J)
      PWATER=1.-PLAND
      PLICE=FLICE(I,J)
      POICE=RSI(I,J)*PWATER
      POCEAN=PWATER-POICE
      PIJ=P(I,J)
      PS=PEDN(1,I,J)
      PSK=PEK(1,I,J)
      P1K=PK(1,I,J)
      Q1=Q(I,J,1)
      THV1=T(I,J,1)*(1.+Q1*deltx)
CCC   JR=JREG(I,J)
      MA1=AM(1,I,J) !@var MA1 mass of lowest atmospheric layer (kg/m^2)
      MSUM = (PS*100.)/GRAV !@var MSUM total mass of atmosphere (kg/m^2)
      PGK = (PS*100.)**KAPA
      PKDN = (GRAV*(MSUM-MA1*0.25))**KAPA
#ifdef TRACERS_ON
C**** Set up tracers for PBL calculation if required
      nx=0
      do n=1,ntm
        if (itime_tr0(n).le.itime .and. needtrs(n)) then
          nx=nx+1
          ntix(nx) = n
C**** Calculate first layer tracer concentration
          trtop(nx)=trm(i,j,1,n)*byam(1,i,j)*bydxyp(j)
        end if
      end do
      ntx = nx
#endif
C**** QUANTITIES ACCUMULATED HOURLY FOR DIAGDD
         IF(MODDD.EQ.0) THEN
         DO KR=1,NDIUPT
           IF(I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
             ADIURN(IH,IDD_SPR,KR)=ADIURN(IH,IDD_SPR,KR)+PS
             ADIURN(IH,IDD_PT5,KR)=ADIURN(IH,IDD_PT5,KR)+PSK*T(I,J,5)
             ADIURN(IH,IDD_PT4,KR)=ADIURN(IH,IDD_PT4,KR)+PSK*T(I,J,4)
             ADIURN(IH,IDD_PT3,KR)=ADIURN(IH,IDD_PT3,KR)+PSK*T(I,J,3)
             ADIURN(IH,IDD_PT2,KR)=ADIURN(IH,IDD_PT2,KR)+PSK*T(I,J,2)
             ADIURN(IH,IDD_PT1,KR)=ADIURN(IH,IDD_PT1,KR)+PSK*T(I,J,1)
             ADIURN(IH,IDD_Q5,KR)=ADIURN(IH,IDD_Q5,KR)+Q(I,J,5)
             ADIURN(IH,IDD_Q4,KR)=ADIURN(IH,IDD_Q4,KR)+Q(I,J,4)
             ADIURN(IH,IDD_Q3,KR)=ADIURN(IH,IDD_Q3,KR)+Q(I,J,3)
             ADIURN(IH,IDD_Q2,KR)=ADIURN(IH,IDD_Q2,KR)+Q(I,J,2)
             ADIURN(IH,IDD_Q1,KR)=ADIURN(IH,IDD_Q1,KR)+Q1
           END IF
         END DO
         END IF
C****
      DO ITYPE=1,3       ! no earth type
      ipbl(i,j,itype)=0
!!!      SELECT CASE (ITYPE)
C****
C**** OPEN OCEAN/LAKE
C****
!!!      CASE (1)
      if ( ITYPE == 1 ) then

      PTYPE=POCEAN
      IF (PTYPE.gt.0) THEN
      TG1=GTEMP(1,1,I,J)
      TG2=GTEMP(2,1,I,J)   ! diagnostic only
      IF (FLAKE(I,J).gt.0) THEN
C**** limit evap/cooling if between MINMLD and 40cm, no evap below 40cm
        IF (MWL(I,J).lt.MINMLD*RHOW*FLAKE(I,J)*DXYP(J)) THEN
          EVAPLIM=MAX(0.5*(MWL(I,J)/(FLAKE(I,J)*DXYP(J))-0.4d0*RHOW),
     *         0d0)
        ELSE
          EVAPLIM=MWL(I,J)/(FLAKE(I,J)*DXYP(J))-(0.5*MINMLD+0.2d0)*RHOW
        END IF
#ifdef TRACERS_WATER
C**** limit on tracer evporation from lake
        TEVAPLIM(1:NTX)=EVAPLIM*GTRACER(NTIX(1:NTX),1,I,J)
#endif
        HTLIM = GML(I,J)/(FLAKE(I,J)*DXYP(J)) + 0.5*LHM*EVAPLIM
        IDTYPE=ITLAKE
        uocean = 0. ; vocean = 0. ! no velocities for lakes
      ELSE
        IDTYPE=ITOCEAN
        if (UOdrag.eq.1) then ! use uocean for drag calculation
C**** Convert UOSURF,VOSURF from C grid to A grid
C**** Note that uosurf,vosurf start with j=1, (not j=2 as in atm winds)
          if (pole) then
            uocean = 0. ; vocean = 0.
            do k=1,kmaxj(j)
              uocean = uocean + rapj(k,j)*(uosurf(idij(k,i,j),idjj(k,j)
     *             -1)*cosip(k)-hemi*vosurf(idij(k,i,j),idjj(k,j)-1)
     *             *sinip(k))
              vocean = vocean + rapj(k,j)*(vosurf(idij(k,i,j),idjj(k,j)
     *             -1)*cosip(k)+hemi*uosurf(idij(k,i,j),idjj(k,j)-1)
     *             *sinip(k))
            end do
          else
            uocean = 0.5*(uosurf(i,j)+uosurf(im1,j))
            vocean = 0.5*(vosurf(i,j)+vosurf(i,j-1))
          end if
        else
          uocean=0. ; vocean=0.
        end if
      END IF
      SRHEAT=FSF(ITYPE,I,J)*COSZ1(I,J)
      SOLAR(1,I,J)=SOLAR(1,I,J)+DTSURF*SRHEAT
            OA(I,J,5)=OA(I,J,5)+SRHEAT*DTSURF
      BETA=1.
      ELHX=LHE
      END IF
C****
C**** OCEAN/LAKE ICE
C****
!!!      CASE (2)
      else if ( ITYPE == 2 ) then

      PTYPE=POICE
      IF (PTYPE.gt.0) THEN
      IF (FLAKE(I,J).gt.0) THEN
        IDTYPE=ITLKICE
        uocean = 0. ; vocean = 0. ! no dynamic ice for lakes
      ELSE
        IDTYPE=ITOICE
        if (UOdrag.eq.1) then ! use ice velcoities in drag calculation
C**** Convert UISURF,VISURF from C grid to A grid
C**** Note that uisurf,visurf start with j=1, (not j=2 as in atm winds)
          if (pole) then
            uocean = 0. ; vocean = 0.
            do k=1,kmaxj(j)
              uocean = uocean + rapj(k,j)*(uisurf(idij(k,i,j),idjj(k,j)
     *             -1)*cosip(k)-hemi*visurf(idij(k,i,j),idjj(k,j)-1)
     *             *sinip(k))
              vocean = vocean + rapj(k,j)*(visurf(idij(k,i,j),idjj(k,j)
     *             -1)*cosip(k)+hemi*uisurf(idij(k,i,j),idjj(k,j)-1)
     *             *sinip(k))
            end do
          else
            uocean = 0.5*(uisurf(i,j)+uisurf(im1,j))
            vocean = 0.5*(visurf(i,j)+visurf(i,j-1))
          end if
        else
          uocean = 0. ; vocean =0.
        end if
      END IF
      SNOW=SNOWI(I,J)
      TG1=TGRND(2,I,J)
      TG2=TGRN2(2,I,J)
      MSI1 = SNOW+ACE1I ! snow and first layer ice mass (kg/m^2)
      ACE2=MSI(I,J) ! second (physical) layer ice mass (kg/m^2)
      dF1dTG = 2./(ACE1I*BYRLI+SNOW*BYRLS)
      HCG1 = SHI*XSI(1)*MSI1 ! heat capacity of top ice layer (J/C*m^2)
      HCG2 = SHI*XSI(2)*MSI1 ! heat capacity of second layer ice
      SRHEAT=FSF(ITYPE,I,J)*COSZ1(I,J)
      SOLAR(2,I,J)=SOLAR(2,I,J)+DTSURF*SRHEAT
C**** fraction of solar radiation leaving layer 1 and 2
      IF (SRHEAT.gt.0) THEN ! only bother if there is sun
        call solar_ice_frac(SNOW,ACE2,FLAG_DSWS(I,J),FSRI,2)
      ELSE
        FSRI(1:2) = 0
      END IF
            OA(I,J,12)=OA(I,J,12)+SRHEAT*DTSURF
      BETA=1.
      ELHX=LHS

      END IF
C****
C**** LAND ICE
C****
!!!      CASE (3)
      else if ( ITYPE == 3 ) then

      PTYPE=PLICE
      IF (PTYPE.gt.0) THEN
      IDTYPE=ITLANDI
      SNOW=SNOWLI(I,J)
      TG1=TGRND(3,I,J)
      TG2=GTEMP(2,3,I,J)
      SRHEAT=FSF(ITYPE,I,J)*COSZ1(I,J)
      Z1BY6L=(Z1LIBYL+SNOW*BYRLS)*BY6
      CDTERM=TG2
      CDENOM=1./(2.*Z1BY6L+Z2LI3L)
      HCG1=HC1LI+SNOW*SHI
      BETA=1.
      ELHX=LHS
      uocean = 0. ; vocean = 0. ! no land ice velocity
      END IF
!!!      END SELECT
      endif
C****
      IF (PTYPE.gt.0) THEN
C****
C**** BOUNDARY LAYER INTERACTION
C****
      TKV=THV1*PSK
      ZS1=.5*DSIG(1)*RGAS*BYGRAV*TKV*PIJ/PS
      SHDT=0.
      EVHDT=0.
      TRHDT=0.
      F1DT=0.

      TG=TG1+TF
      QG_SAT=QSAT(TG,ELHX,PS)
      IF (ITYPE.eq.1 .and. focean(i,j).gt.0) QG_SAT=0.98d0*QG_SAT
      TGV=TG*(1.+QG_SAT*deltx)
#ifdef TRACERS_ON
C**** Set up b.c. for tracer PBL calculation if required
      do nx=1,ntx
        n=ntix(nx)
C**** Set surface boundary conditions for tracers depending on whether
C**** they are water or another type of tracer
#ifdef TRACERS_WATER
        tr_evap_max(nx)=1.
C**** The select is used to distinguish water from gases or particle
!!!        select case (tr_wd_TYPE(n))
!!!        case (nWATER)
        if ( tr_wd_TYPE(n) == nWATER ) then
#ifdef TRACERS_SPECIAL_O18
          if (ITYPE.eq.1) then   ! liquid water => fractionation
        trgrnd(nx)=gtracer(n,itype,i,j)*QG_SAT*FRACVL(TG1,trname(n))
          else   ! ice, no fractionation
            trgrnd(nx)=gtracer(n,itype,i,j)*QG_SAT
          end if
#else
            trgrnd(nx)=gtracer(n,itype,i,j)*QG_SAT
#endif
C**** trsfac and trconstflx are multiplied by cq*wsh in PBL
          trsfac(nx)=1.
          trconstflx(nx)=trgrnd(nx)
!!!        case (nGAS, nPART)
        else if ( tr_wd_TYPE(n) == nGAS .or.
     &         tr_wd_TYPE(n) == nPART ) then
#endif
C**** For non-water tracers (i.e. if TRACERS_WATER is not set, or there
C**** is a non-soluble tracer mixed in.)
C**** Calculate trsfac (set to zero for const flux)
          trsfac(nx)=0.
C**** Calculate trconstflx (m/s * conc) (could be dependent on itype)
          rhosrf0=100.*ps/(rgas*tgv) ! estimated surface density
          totflux=0.
          do nsrc=1,ntsurfsrc(n)
            totflux = totflux+trsource(i,j,nsrc,n)
          end do
          trconstflx(nx)=totflux/(dxyp(j)*rhosrf0)
#ifdef TRACERS_WATER
!!!        end select
        endif
#endif
      end do
#endif
C =====================================================================
      fr_sat = 1. ! entire surface is saturated
      evap_max = 1.
      CALL PBL(I,J,ITYPE,PTYPE)
      CM = cmgs(i,j,itype)
      CH = chgs(i,j,itype)
      CQ = cqgs(i,j,itype)
      DHGS=(ZMIX-ZGS)*CH*WSH
      DQGS=(ZMIX-ZGS)*CQ*WSH
      DGS =DQGS
C =====================================================================
      TS=TSV/(1.+QSRF*deltx)
C**** CALCULATE RHOSRF*CM*WS AND RHOSRF*CH*WS
      RHOSRF=100.*PS/(RGAS*TSV)
      RCDMWS=CM*WSM*RHOSRF
      RCDHWS=CH*WSH*RHOSRF
      RCDQWS=CQ*WSH*RHOSRF
C**** CALCULATE FLUXES OF SENSIBLE HEAT, LATENT HEAT, THERMAL
C****   RADIATION, AND CONDUCTION HEAT (WATTS/M**2)
      SHEAT=SHA*RCDHWS*(TS-TG)
      BETAUP = BETA
      IF (QSRF .GT. QG_SAT) BETAUP = 1.
      EVHEAT=(LHE+TG1*SHV)*BETAUP*RCDQWS*(QSRF-QG_SAT)
      TRHEAT=TRHR(0,I,J)-STBO*(TG*TG)*(TG*TG)
C****
!!!      SELECT CASE (ITYPE)

!!!      CASE (1) ! FLUXES USING IMPLICIT TIME STEP FOR OCEAN POINTS
      if ( ITYPE == 1) then
        DSHDTG=-RCDHWS*SHA
        dEVdQS = LHE*RCDQWS
        dHS = -(1.+2.*S1BYG1)*DTSURF*PGK*SHEAT/
     A       ((1.+2.*S1BYG1)*DTSURF*PGK*RCDHWS+MA1*PKDN)
        dTS = -(1.+2.*S1BYG1)*DTSURF*PGK*SHEAT/
     A       (MA1*PKDN*SHA-(1.+2.*S1BYG1)*DTSURF*PGK*DSHDTG)
        dQS = -(1.+2.*S1BYG1)*DTSURF*EVHEAT/
     A       ((1.+2.*S1BYG1)*DTSURF*dEVdQS+MA1*LHE)
        SHDT = DTSURF*(SHEAT-dTS*DSHDTG)
        EVHDT=DTSURF*(EVHEAT+dQS*dEVdQS) ! latent heat flux
        TRHDT=DTSURF*TRHEAT
C****
!!!      CASE (2) ! FLUXES USING IMPLICIT TIME STEP FOR ICE POINTS
      else if ( ITYPE == 2 ) then

! heat flux on first/second/third layers (W/m^2)
        F1 = (TG1-TG2)*dF1dTG + SRHEAT*FSRI(1)
        F2 = SRHEAT*FSRI(2)
        EVHEAT = LHE*RCDQWS*(QSRF-QG_SAT) ! latent heat flux (W/m^2)
        F0=SRHEAT+TRHEAT+SHEAT+EVHEAT
        dSNdTG=-RCDHWS*KHS*SHA/(DHGS+KHS)
        dQGdTG=QG_SAT*DQSATDT(TG,ELHX) ! d(QG)/dTG
        dEVdTG = -dQGdTG*LHE*RCDQWS*KHS/(DGS+KHS) ! d(EVHEAT)/dTG
        dTRdTG = -4*STBO*TG*TG*TG ! d(TRHEAT)/dTG
        dF0dTG = dSNdTG+dEVdTG+dTRdTG ! d(F0)/dTG
C       dSNdHS = RCDHWS ! d(SHEAT)/dHS - kg/(sec*m^2)
        dEVdQS = LHE*RCDQWS     ! d(EVHEAT)/dQS
        HSDEN = -(1.+2.*S1BYG1)*DTSURF*PGK*BETA*dSNdTG+MA1*PKDN*SHA
        HSCON = -(1.+2.*S1BYG1)*DTSURF*PGK*SHEAT/HSDEN ! (J*sec)/kg
        HSMUL=-(1.+2.*S1BYG1)*DTSURF*PGK*BETA*dSNdTG/HSDEN ! J/(kg*degC)
        QSDEN = (1.+2.*S1BYG1)*BETA*DTSURF*dEVdQS+MA1*LHE
        QSCON = -(1.+2.*S1BYG1)*DTSURF*EVHEAT/QSDEN
        QSMUL = -(1.+2.*S1BYG1)*DTSURF*BETA*dEVdTG/QSDEN
        T2DEN = HCG2+BETA*DTSURF*dF1dTG
        T2CON = DTSURF*(F1-F2)/T2DEN
        T2MUL = BETA*DTSURF*dF1dTG/T2DEN
        TGDEN = HCG1-BETA*DTSURF*(dF0dTG-dF1dTG-
     A       HSMUL*dSNdTG+QSMUL*dEVdQS+T2MUL*dF1dTG) ! W/(m^2*degC)
        dTG = DTSURF*(F0-F1+BETA*
     A       (QSCON*dEVdQS-HSCON*dSNdTG+T2CON*dF1dTG))/TGDEN ! degC
        IF (TG1+dTG .GT. 0.) dTG = -TG1
        dHS = HSCON+HSMUL*dTG   ! (J*sec)/kg
        dQS = QSCON+QSMUL*dTG
        dT2 = T2CON+T2MUL*dTG
        SHDT = DTSURF*(SHEAT+BETA*((dTG-dHS)*dSNdTG)) ! sensible
        EVHDT = DTSURF*(EVHEAT+BETA*(dTG*dEVdTG+dQS*dEVdQS)) ! latent
        TRHDT = DTSURF*(TRHEAT+BETA*dTG*dTRdTG) ! thermal flux (J/m^2)
        F1DT = DTSURF*(F1+BETA*(dTG*dF1dTG-dT2*dF1dTG))
        TG1 = TG1+dTG          ! first layer sea ice temperature (degC)
        TG2 = TG2+dT2          ! second layer sea ice temperature (degC)
        TGRN2(ITYPE,I,J) = TG2
C****
!!!      CASE (3) ! IMPLICIT TIME STEP OVER LANDICE
      else if ( ITYPE == 3 ) then

        F0=SRHEAT+TRHEAT+SHEAT+EVHEAT
        F1=(TG1-CDTERM-F0*Z1BY6L)*CDENOM
        DSHDTG=-RCDHWS*KHS*SHA/(DHGS+KHS)
        DQGDTG=QG_SAT*DQSATDT(TG,ELHX)
        DEVDTG=-RCDQWS*KHS*LHE*BETAUP*DQGDTG/(DQGS+KHS)
        DTRDTG=-4.*STBO*TG*TG*TG
        DF0DTG=DSHDTG+DEVDTG+DTRDTG
        DFDTG=DF0DTG-(1.-DF0DTG*Z1BY6L)*CDENOM
        DTG=(F0-F1)*DTSURF/(HCG1-DTSURF*DFDTG)
        SHDT=DTSURF*(SHEAT+DTG*DSHDTG)
        EVHDT=DTSURF*(EVHEAT+DTG*DEVDTG)
        TRHDT=DTSURF*(TRHEAT+DTG*DTRDTG)
        F1DT=DTSURF*(TG1-CDTERM-(F0+DTG*DFDTG)*Z1BY6L)*CDENOM
        TG1=TG1+DTG

!!!      END SELECT
      endif

C**** CALCULATE EVAPORATION
      DQ1X =EVHDT/((LHE+TG1*SHV)*MA1)
      EVHDT0=EVHDT
C**** Limit evaporation if lake mass is at minimum
      IF (ITYPE.EQ.1 .and. FLAKE(I,J).GT.0 .and.
     *     (EVAPOR(I,J,1)-DQ1X*MA1).gt.EVAPLIM) THEN
        if (QCHECK) WRITE(99,*) "Lake EVAP limited: I,J,EVAP,MWL",I,J
     *     ,EVAPOR(I,J,1)-DQ1X*MA1,MWL(I,J)/(RHOW*FLAKE(I,J)*DXYP(J))
        DQ1X=(EVAPOR(I,J,1)-EVAPLIM)*BYAM(1,I,J)
      ELSEIF (DQ1X.GT.Q1+DQ1(I,J)) THEN
        DQ1X=(Q1+DQ1(I,J))
      ELSE
        GO TO 3720
      END IF
      EVHDT=DQ1X*(LHE+TG1*SHV)*MA1
      IF (ITYPE.NE.1) TG1=TG1+(EVHDT-EVHDT0)/HCG1
 3720 EVAP=-DQ1X*MA1
#ifdef TRACERS_WATER
C****
C**** Calculate Water Tracer Evaporation
C****
      DO NX=1,NTX
        N=NTIX(NX)
        if (tr_wd_TYPE(n).eq.nWATER) THEN
          IF (ITYPE.EQ.1) THEN  ! OCEAN
C**** do calculation implicitly for TQS
#ifdef TRACERS_SPECIAL_O18
            TEV=-RCDQWS*(trs(nx)-trgrnd(nx))*FRACLK(WSM,trname(n))
            dTEVdTQS =-RCDQWS*FRACLK(WSM,trname(n))
#else
            TEV=-RCDQWS*(trs(nx)-trgrnd(nx))
            dTEVdTQS =-RCDQWS
#endif
            dTQS = -(1.+2.*S1BYG1)*DTSURF*TEV/
     *           ((1.+2.*S1BYG1)*DTSURF*dTEVdTQS-MA1)
            TEVAP=DTSURF*(TEV+dTQS*dTEVdTQS)
          ELSE                  ! ICE AND LAND ICE
C**** tracer flux is set by source tracer concentration
            IF (EVAP.GE.0) THEN ! EVAPORATION
              TEVAP=EVAP*trgrnd(nx)/QG_SAT
            ELSE                ! DEW (fractionates)
#ifdef TRACERS_SPECIAL_O18
              IF (TG1.gt.0) THEN
                frac=FRACVL(TG1,trname(n))
              ELSE
                frac=FRACVS(TG1,trname(n))
              END IF
              TEVAP=EVAP*trs(nx)/(QSRF*frac+teeny)
#else
              TEVAP=EVAP*trs(nx)/(QSRF+teeny)
#endif
            END IF
          END IF
C**** Limit evaporation if lake mass is at minimum
          IF (ITYPE.EQ.1 .and. FLAKE(I,J).GT.0 .and.
     *         (TREVAPOR(n,1,I,J)+TEVAP.gt.TEVAPLIM(NX))) THEN
            IF (QCHECK) WRITE(99,*) "Lake TEVAP limited: I,J,TEVAP,TMWL"
     *           ,N,TREVAPOR(n,1,I,J)+TEVAP,TEVAPLIM(NX)
            TEVAP= TEVAPLIM(NX)-TREVAPOR(n,1,I,J)
          END IF
          TDP = TEVAP*DXYP(J)*ptype
          TDT1 = trsrfflx(I,J,n)*DTSURF
          IF (TRM(I,J,1,n)+TDT1+TDP.lt.0..and.tdp.lt.0) THEN
            IF (QCHECK) WRITE(99,*) "LIMITING TRDEW",I,J,N,TDP,TRM(I,J,1
     *           ,n),TDT1
            TEVAP = -(TRM(I,J,1,n)+TDT1)/(DXYP(J)*ptype)
            trsrfflx(I,J,n)= - TRM(I,J,1,n)/DTSURF
          ELSE
            trsrfflx(I,J,n)=trsrfflx(I,J,n)+TDP/DTSURF
          END IF
          TREVAPOR(n,ITYPE,I,J)=TREVAPOR(n,ITYPE,I,J)+TEVAP
        END IF
      END DO
#endif
C**** ACCUMULATE SURFACE FLUXES AND PROGNOSTIC AND DIAGNOSTIC QUANTITIES
      F0DT=DTSURF*SRHEAT+TRHDT+SHDT+EVHDT
C**** Limit heat fluxes out of lakes if near minimum depth
      IF (ITYPE.eq.1 .and. FLAKE(I,J).gt.0 .and.
     *     E0(I,J,1)+F0DT+HTLIM.lt.0) THEN
        if (QCHECK) write(6,*) "Limiting heat flux from lake",i,j,SHDT
     *       ,F0DT,E0(I,J,1),DTSURF*SRHEAT,TRHDT,EVHDT,HTLIM
        SHDT = -(HTLIM+E0(I,J,1)+DTSURF*SRHEAT+TRHDT+EVHDT)
        F0DT = -E0(I,J,1)-HTLIM
        if (QCHECK) write(6,*) "New SHDT,F0DT",i,j,SHDT,F0DT
      END IF
      E0(I,J,ITYPE)=E0(I,J,ITYPE)+F0DT
      E1(I,J,ITYPE)=E1(I,J,ITYPE)+F1DT
      EVAPOR(I,J,ITYPE)=EVAPOR(I,J,ITYPE)+EVAP
      TGRND(ITYPE,I,J)=TG1
      DTH1(I,J)=DTH1(I,J)-SHDT*PTYPE/(SHA*MA1*P1K)
      DQ1(I,J) =DQ1(I,J) -DQ1X*PTYPE
      DMUA(I,J,ITYPE)=DMUA(I,J,ITYPE)+PTYPE*DTSURF*RCDMWS*(US-UOCEAN)
      DMVA(I,J,ITYPE)=DMVA(I,J,ITYPE)+PTYPE*DTSURF*RCDMWS*(VS-VOCEAN)
      uflux1(i,j)=uflux1(i,j)+PTYPE*RCDMWS*(US-UOCEAN)
      vflux1(i,j)=vflux1(i,j)+PTYPE*RCDMWS*(VS-UOCEAN)
C****
C**** ACCUMULATE DIAGNOSTICS FOR EACH SURFACE TIME STEP AND ITYPE
C****
        AJ(J,J_EVAP ,IDTYPE)=AJ(J,J_EVAP ,IDTYPE)+ EVAP*PTYPE
        AJ(J,J_EVHDT,IDTYPE)=AJ(J,J_EVHDT,IDTYPE)+EVHDT*PTYPE
        AJ(J,J_SHDT ,IDTYPE)=AJ(J,J_SHDT ,IDTYPE)+ SHDT*PTYPE
        AJ(J,J_TRHDT,IDTYPE)=AJ(J,J_TRHDT,IDTYPE)+TRHDT*PTYPE
        IF(MODDSF.EQ.0) THEN
          AJ(J,J_TSRF,IDTYPE)=AJ(J,J_TSRF,IDTYPE)+(TS-TF)*PTYPE
          AJ(J,J_TYPE,IDTYPE)=AJ(J,J_TYPE,IDTYPE)+        PTYPE
          AJ(J,J_TG1 ,IDTYPE)=AJ(J,J_TG1 ,IDTYPE)+    TG1*PTYPE
          AJ(J,J_TG2 ,IDTYPE)=AJ(J,J_TG2 ,IDTYPE)+    TG2*PTYPE
        END IF
C**** QUANTITIES ACCUMULATED FOR REGIONS IN DIAGJ
CCC     AREG(JR,J_TRHDT)=AREG(JR,J_TRHDT)+TRHDT*PTYPE*DXYP(J)
CCC     AREG(JR,J_SHDT )=AREG(JR,J_SHDT )+SHDT *PTYPE*DXYP(J)
CCC     AREG(JR,J_EVHDT)=AREG(JR,J_EVHDT)+EVHDT*PTYPE*DXYP(J)
CCC     AREG(JR,J_EVAP )=AREG(JR,J_EVAP )+EVAP *PTYPE*DXYP(J)
CCC     IF(MODDSF.EQ.0)
CCC  *       AREG(JR,J_TSRF)=AREG(JR,J_TSRF)+(TS-TF)*PTYPE*DXYP(J)
        AREGIJ(1,ITYPE,I,J)=TRHDT*PTYPE*DXYP(J)
        AREGIJ(2,ITYPE,I,J)=SHDT *PTYPE*DXYP(J)
        AREGIJ(3,ITYPE,I,J)=EVHDT*PTYPE*DXYP(J)
        AREGIJ(4,ITYPE,I,J)=EVAP *PTYPE*DXYP(J)
        IF(MODDSF.EQ.0) THEN
          AREGIJ(5,ITYPE,I,J)=(TS-TF)*PTYPE*DXYP(J)
          AREGIJ(6,ITYPE,I,J)=    TG1*PTYPE*DXYP(J)
          AREGIJ(7,ITYPE,I,J)=    TG2*PTYPE*DXYP(J)
        END IF
C
C**** QUANTITIES ACCUMULATED FOR LATITUDE-LONGITUDE MAPS IN DIAGIJ
        AIJ(I,J,IJ_SHDT)=AIJ(I,J,IJ_SHDT)+SHDT*PTYPE
        IF(MODRD.EQ.0) AIJ(I,J,IJ_TRNFP0)=AIJ(I,J,IJ_TRNFP0)+TRHDT
     *       *PTYPE/DTSRC
        AIJ(I,J,IJ_SRTR)=AIJ(I,J,IJ_SRTR)+(SRHEAT*DTSURF+TRHDT)*PTYPE
        AIJ(I,J,IJ_NETH)=AIJ(I,J,IJ_NETH)+(SRHEAT*DTSURF+TRHDT+SHDT
     *       +EVHDT)*PTYPE
        AIJ(I,J,IJ_EVAP)=AIJ(I,J,IJ_EVAP)+EVAP*PTYPE
        IF(MODDSF.EQ.0) THEN
          AIJ(I,J,IJ_WS)=AIJ(I,J,IJ_WS)+WS*PTYPE
          AIJ(I,J,IJ_TS)=AIJ(I,J,IJ_TS)+(TS-TF)*PTYPE
          AIJ(I,J,IJ_US)=AIJ(I,J,IJ_US)+US*PTYPE
          AIJ(I,J,IJ_VS)=AIJ(I,J,IJ_VS)+VS*PTYPE
          AIJ(I,J,IJ_TAUS)=AIJ(I,J,IJ_TAUS)+RCDMWS*WSM*PTYPE
          AIJ(I,J,IJ_TAUUS)=AIJ(I,J,IJ_TAUUS)+RCDMWS*(US-UOCEAN)*PTYPE
          AIJ(I,J,IJ_TAUVS)=AIJ(I,J,IJ_TAUVS)+RCDMWS*(VS-VOCEAN)*PTYPE
          AIJ(I,J,IJ_QS)=AIJ(I,J,IJ_QS)+QSRF*PTYPE
          AIJ(I,J,IJ_TG1)=AIJ(I,J,IJ_TG1)+TG1*PTYPE
        END IF
C**** QUANTITIES ACCUMULATED HOURLY FOR DIAGDD
        IF(MODDD.EQ.0) THEN
          DO KR=1,NDIUPT
            IF(I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
              ADIURN(IH,IDD_TS,KR)=ADIURN(IH,IDD_TS,KR)+TS*PTYPE
              ADIURN(IH,IDD_TG1,KR)=ADIURN(IH,IDD_TG1,KR)+(TG1+TF)*PTYPE
              ADIURN(IH,IDD_QS,KR)=ADIURN(IH,IDD_QS,KR)+QSRF*PTYPE
              ADIURN(IH,IDD_QG,KR)=ADIURN(IH,IDD_QG,KR)+QG_SAT*PTYPE
              ADIURN(IH,IDD_SWG,KR)=ADIURN(IH,IDD_SWG,KR)+SRHEAT*DTSURF
     *             *PTYPE
              ADIURN(IH,IDD_LWG,KR)=ADIURN(IH,IDD_LWG,KR)+TRHDT*PTYPE
              ADIURN(IH,IDD_SH,KR)=ADIURN(IH,IDD_SH,KR)+SHDT*PTYPE
              ADIURN(IH,IDD_LH,KR)=ADIURN(IH,IDD_LH,KR)+EVHDT*PTYPE
              ADIURN(IH,IDD_HZ0,KR)=ADIURN(IH,IDD_HZ0,KR)
     *             +(SRHEAT*DTSURF+TRHDT+SHDT+EVHDT)*PTYPE
              ADIURN(IH,IDD_UG,KR)=ADIURN(IH,IDD_UG,KR)+UG*PTYPE
              ADIURN(IH,IDD_VG,KR)=ADIURN(IH,IDD_VG,KR)+VG*PTYPE
              ADIURN(IH,IDD_WG,KR)=ADIURN(IH,IDD_WG,KR)+WG*PTYPE
              ADIURN(IH,IDD_US,KR)=ADIURN(IH,IDD_US,KR)+US*PTYPE
              ADIURN(IH,IDD_VS,KR)=ADIURN(IH,IDD_VS,KR)+VS*PTYPE
              ADIURN(IH,IDD_WS,KR)=ADIURN(IH,IDD_WS,KR)+WS*PTYPE
              ADIURN(IH,IDD_CIA,KR)=ADIURN(IH,IDD_CIA,KR)+PSI*PTYPE
              ADIURN(IH,IDD_CM,KR)=ADIURN(IH,IDD_CM,KR)+CM*PTYPE
              ADIURN(IH,IDD_CH,KR)=ADIURN(IH,IDD_CH,KR)+CH*PTYPE
              ADIURN(IH,IDD_CQ,KR)=ADIURN(IH,IDD_CQ,KR)+CQ*PTYPE
              ADIURN(IH,IDD_EDS,KR)=ADIURN(IH,IDD_EDS,KR)+KHS*PTYPE
              ADIURN(IH,IDD_DBL,KR)=ADIURN(IH,IDD_DBL,KR)+DBL*PTYPE
              ADIURN(IH,IDD_EV,KR)=ADIURN(IH,IDD_EV,KR)+EVAP*PTYPE
            END IF
          END DO
        END IF
C****
#ifdef TRACERS_ON
C**** Save surface tracer concentration whether calculated or not
      nx=0
      do n=1,ntm
        if (itime_tr0(n).le.itime) then
          if (needtrs(n)) then
            nx=nx+1
            taijn(i,j,tij_surf,n) = taijn(i,j,tij_surf,n)+trs(nx)*ptype
          else
            taijn(i,j,tij_surf,n) = taijn(i,j,tij_surf,n)
     *           +max((trm(i,j,1,n)-trmom(mz,i,j,1,n))*byam(1,i,j)
     *           *bydxyp(j),0d0)*ptype
          end if
#ifdef TRACERS_WATER
          taijn(i,j,tij_evap,n)=taijn(i,j,tij_evap,n)+
     *         trevapor(n,itype,i,j)*ptype
          taijn(i,j,tij_grnd,n)=taijn(i,j,tij_grnd,n)+
     *         gtracer(n,itype,i,j)*ptype/nisurf
          tajls(j,1,jls_source(1,n))=tajls(j,1,jls_source(1,n))
     *         +trevapor(n,itype,i,j)*ptype
          if (focean(i,j).gt.0) tajls(j,1,jls_source(2,n))=tajls(j,1
     *         ,jls_source(2,n))+trevapor(n,itype,i,j)*ptype
#endif
        end if
      end do
#endif
C****
C**** SAVE SOME TYPE DEPENDENT FLUXES/DIAGNOSTICS
C****
!!!      SELECT CASE (ITYPE)
!!!      CASE (1)  ! ocean
      if ( ITYPE == 1 ) then
        OA(I,J,6)=OA(I,J,6)+TRHDT
        OA(I,J,7)=OA(I,J,7)+SHDT
        OA(I,J,8)=OA(I,J,8)+EVHDT
        IF (MODDSF.eq.0) AIJ(I,J,IJ_TGO)=AIJ(I,J,IJ_TGO)+TG1
        AIJ(I,J,IJ_EVAPO)=AIJ(I,J,IJ_EVAPO)+EVAP*PTYPE
        IF (FOCEAN(I,J).gt.0) AIJ(I,J,IJ_F0OC)=AIJ(I,J,IJ_F0OC) +F0DT
     *       *PTYPE
C****
!!!      CASE (2)  ! seaice
      else if ( ITYPE == 2 ) then
        IF (TG1.GT.TDIURN(I,J,7)) TDIURN(I,J,7) = TG1
        OA(I,J,9)=OA(I,J,9)+TRHDT
        OA(I,J,10)=OA(I,J,10)+SHDT
        OA(I,J,11)=OA(I,J,11)+EVHDT
        AIJ(I,J,IJ_F0OI) =AIJ(I,J,IJ_F0OI) +F0DT*PTYPE
        AIJ(I,J,IJ_EVAPI)=AIJ(I,J,IJ_EVAPI)+EVAP*PTYPE
C****
!!!      CASE (3) ! land ice
      else if ( ITYPE == 3 ) then
        IF (TG1.GT.TDIURN(I,J,8)) TDIURN(I,J,8) = TG1
        IF (MODDSF.eq.0) AIJ(I,J,IJ_TSLI)=AIJ(I,J,IJ_TSLI)+(TS-TF)
        AIJ(I,J,IJ_SHDTLI)=AIJ(I,J,IJ_SHDTLI)+SHDT
        AIJ(I,J,IJ_EVHDT)=AIJ(I,J,IJ_EVHDT)+EVHDT
        AIJ(I,J,IJ_TRHDT)=AIJ(I,J,IJ_TRHDT)+TRHDT
        AIJ(I,J,IJ_F0LI)=AIJ(I,J,IJ_F0LI)+F0DT
        AIJ(I,J,IJ_EVAPLI)=AIJ(I,J,IJ_EVAPLI)+EVAP
C****
!!!      END SELECT
      endif
C****
      END IF
      END DO   ! end of itype loop
      IM1=I
      END DO   ! end of I loop

      END DO   ! end of J loop
!$OMP  END PARALLEL DO

      DO J=1,JM
      DO I=1,IMAXJ(J)
        JR=JREG(I,J)
        DO K=1,3
          AREG(JR,J_TRHDT)=AREG(JR,J_TRHDT)+AREGIJ(1,K,I,J)
          AREG(JR,J_SHDT )=AREG(JR,J_SHDT )+AREGIJ(2,K,I,J)
          AREG(JR,J_EVHDT)=AREG(JR,J_EVHDT)+AREGIJ(3,K,I,J)
          AREG(JR,J_EVAP )=AREG(JR,J_EVAP )+AREGIJ(4,K,I,J)
          IF(MODDSF.EQ.0) THEN
            AREG(JR,J_TSRF )=AREG(JR,J_TSRF )+AREGIJ(5,K,I,J)
            AREG(JR,J_TG1)=AREG(JR,J_TG1)+AREGIJ(6,K,I,J)
            AREG(JR,J_TG2)=AREG(JR,J_TG2)+AREGIJ(7,K,I,J)
          END IF
        END DO
      END DO
      END DO
C****
C**** EARTH
C****
      CALL EARTH(NS,MODDSF,MODDD)
C****
C**** UPDATE FIRST LAYER QUANTITIES
C****
!$OMP  PARALLEL DO PRIVATE (I,J,FTEVAP,FQEVAP,P1K)
!$OMP*          SCHEDULE(DYNAMIC,2)
      DO J=1,JM
      DO I=1,IMAXJ(J)
        FTEVAP=0
        IF (DTH1(I,J)*T(I,J,1).lt.0) FTEVAP=-DTH1(I,J)/T(I,J,1)
        FQEVAP=0
        IF (DQ1(I,J).lt.0.and.Q(I,J,1).gt.0) FQEVAP=-DQ1(I,J)/Q(I,J,1)
! Z-moments should be set from PBL
        TMOM(:,I,J,1) = TMOM(:,I,J,1)*(1.-FTEVAP)
        QMOM(:,I,J,1) = QMOM(:,I,J,1)*(1.-FQEVAP)
        IF ( Q(I,J,1)+DQ1(I,J) .LT. qmin ) THEN
          WRITE(99,*)
     &         ITime,'I,J:',I,J,' Q1:',Q(I,J,1)+DQ1(I,J),'->',qmin
          dq1(i,j)=qmin-q(i,j,1)
          QMOM(:,I,J,1)=0.
        ENDIF
c****   retrieve fluxes
        P1K=PK(1,I,J)
        tflux1(i,j)=dth1(i,j)*(-AM(1,I,J)*P1K)/(dtsurf)
        qflux1(i,j)=dq1(i,j)*(-AM(1,I,J))/(dtsurf)
C**** Diurnal cycle of temperature diagnostics
        tdiurn(i,j,5)=tdiurn(i,j,5)+(tsavg(i,j)-tf)
        if(tsavg(i,j).gt.tdiurn(i,j,6)) tdiurn(i,j,6)=tsavg(i,j)
        if(tsavg(i,j).lt.tdiurn(i,j,9)) tdiurn(i,j,9)=tsavg(i,j)
      END DO
      END DO
!$OMP  END PARALLEL DO
#ifdef TRACERS_ON
C****
C**** Apply tracer surface sources and sinks
C****
      call apply_tracer_2Dsource(dtsurf)
#endif
c****
c**** apply surface fluxes to the first layer of the atmosphere
c****  (replaced with dummy subroutine when ATURB is used)
c****
      call apply_fluxes_to_atm(dtsurf)

C**** Call dry convection or aturb depending on rundeck
      CALL ATM_DIFFUS(1,1,dtsurf)
C****
C**** ACCUMULATE SOME ADDITIONAL BOUNDARY LAYER DIAGNOSTICS
C****
      IF(MODDD.EQ.0) THEN
        DO KR=1,NDIUPT
C**** CHECK IF DRY CONV HAS HAPPENED FOR THIS DIAGNOSTIC
          IF(DCLEV(IJDD(1,KR),IJDD(2,KR)).GT.1.) THEN
            ADIURN(IH,IDD_DCF,KR)=ADIURN(IH,IDD_DCF,KR)+1.
            ADIURN(IH,IDD_LDC,KR)=ADIURN(IH,IDD_LDC,KR)
     *           +DCLEV(IJDD(1,KR),IJDD(2,KR))
          END IF
        END DO
      END IF
C****
      END DO   ! end of surface time step
      RETURN
C****
      END SUBROUTINE SURFCE

