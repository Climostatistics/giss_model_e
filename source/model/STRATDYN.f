!@sum  STRATDYN stratospheric only routines
!@auth Bob Suozzo/Jean Lerner/Gavin Schmidt
!@ver  1.0

C**** TO DO:
C****   i) A-grid <-> B-grid  should be done with indexes etc.
C****  ii) PK type variables should be done in dynamics and used here

      MODULE STRAT
!@sum  STRAT local stratospheric variables for GW drag etc.
!@auth Bob Suozzo/Jean Lerner
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm
      IMPLICIT NONE
      SAVE
!@dbparam XCDNST parameters for GW drag (in param. database)
      REAL*8, DIMENSION(2) :: XCDNST(2)
!@dbparam CMTN parameter for GW MTN drag (in param. database)
!@dbparam CDEF parameter for GW DEF drag (in param. database)
!@dbparma CMC parameter for GW M. Convective drag (in param. database)
C**** (used to be FMC)
      REAL*8 :: CMTN = .5, CDEF = 3., CMC = 2d-7
!@dbparam PBREAK p. level above which GW drag acts (in param. database)
      REAL*8 :: PBREAK = 500.   ! default is 500mb
!@dbparam PCONPEN level of penetrating moist conv (in param. database)
      REAL*8 :: PCONPEN = 400.   ! default is 400mb
!@dbparam DEFTHRESH threshold for deformation wave (1/s)
      REAL*8 :: DEFTHRESH = 15d-6  ! default is 15x10^-6 s^-1
!@dbparam PBREAKTOP p. level to force GW breaking in top layer
C**** This should be set to 100. (or something similar) to force
C**** breaking of remaining gravity waves in top layer. Otherwise,
C**** momentum passes through model top.
      REAL*8 :: PBREAKTOP = 0.05d0   ! default is 0.05mb

!@var ZVART,ZVARX,ZVARY,ZWT topogrpahic variance
C**** (must be in common due to read statement)
      REAL*8, DIMENSION(IM,JM) :: ZVART,ZVARX,ZVARY,ZWT
      COMMON/ZVARCB/ZVART,ZVARX,ZVARY,ZWT
!@var DEFRM deformation field
      REAL*8, DIMENSION(IM,JM) :: DEFRM
!@var LDEF,LDEFM deformation levels
      INTEGER LDEF,LDEFM
!@var LBREAK,LSHR,LD2 levels for various GW drag terms
      INTEGER :: LBREAK,LSHR,LD2 = LM   ! need default for LD2

!@dbparam QGWMTN =1 turns on GW Mountain Wave drag terms
!@dbparam QGWSHR =1 turns on GW Shear drag terms
!@dbparam QGWDEF =1 turns on GW Deformation drag terms
!@dbparam QGWCNV =1 turns on GW Convective drag terms
      INTEGER :: QGWMTN = 1, QGWSHR = 1, QGWDEF = 1, QGWCNV = 1

!@var PK local P**Kapa array - should be done by DYNAMICS?
      REAL*8, DIMENSION(IM,JM,LM) :: PK
!@param NM number of gravity wave drag sources
      INTEGER, PARAMETER :: NM=9
!@var Arrays needed for GWDRAG
      REAL*8 :: EK(NM,JM),PKS(LM)

      END MODULE


      SUBROUTINE init_GWDRAG
!@sum init_GWDRAG
!@auth Jean Lerner
C**** DO_GWDRAG=true activates the printing of the diagnostics
C**** accumulated in the routines contained herein
      USE FILEMANAGER
      USE PARAM
      USE CONSTANT, only : twopi,kapa
      USE MODEL_COM, only : im,jm,lm,ls1,do_gwdrag,ptop,sig,psfmpt,sige
      USE GEOM, only : areag,dxyv,dlat_dg
      USE STRAT, only : xcdnst, qgwmtn, qgwshr, qgwdef, qgwcnv,lbreak
     *     ,ld2,lshr,ldef,ldefm,zvarx,zvary,zvart,zwt,pks,nm,ek, cmtn
     *     ,cdef,cmc,pbreak,pbreaktop,defthresh,pconpen
      IMPLICIT NONE
      REAL*8 PLEV,PLEVE,EKS,EK1,EK2,EKX
      INTEGER I,J,L,iu_zvar

C**** define flag for optional diagnostics
      DO_GWDRAG = .true.

C**** sync gwdrag parameters from input
      call sync_param( "XCDNST", XCDNST, 2 )
      call sync_param( "CMTN", CMTN)
      call sync_param( "CDEF", CDEF)
      call sync_param( "CMC", CMC)
      call sync_param( "PBREAK", PBREAK)
      call sync_param( "PCONPEN", PCONPEN)
      call sync_param( "PBREAKTOP", PBREAKTOP)
      call sync_param( "DEFTHRESH", DEFTHRESH)

C**** sync more gwdrag parameters from input
      call sync_param( "QGWMTN", QGWMTN)
      call sync_param( "QGWSHR", QGWSHR)
      call sync_param( "QGWDEF", QGWDEF)
      call sync_param( "QGWCNV", QGWCNV)

C**** Calculate levels for deformation etc.
C**** Note: these levels work for the 23 layer model, but may
C**** need testing for other resolutions
      DO L=1,LM
        PLEV=PSFMPT*SIG(L)+PTOP
        PLEVE=PSFMPT*SIGE(L)+PTOP
        IF (PLEV.GE.700) LDEF=L
        IF (PLEVE.GE.PBREAK) LBREAK=L+1
        IF (PLEVE.GE.300.) LSHR=L
        IF (PLEV.GE.200.) LDEFM=L
        IF (PLEV.GE.0.2d0) LD2=L
      END DO
      WRITE (*,*) ' LEVEL FOR DEFORMATION IS: LDEF,PDEF= ',LDEF,PSFMPT
     *     *SIG(LDEF)+PTOP,' LDEFM=',LDEFM
      WRITE (*,*) ' LEVELS FOR WIND SHEAR GENERATION: LSHR,LD2= ',LSHR
     *     ,LD2
      WRITE (*,*) ' LBREAK=',LBREAK
C****
C**** TOPOGRAPHY VARIANCE FOR MOUNTAIN WAVES
C****
      call openunit("ZVAR",iu_ZVAR,.true.,.true.)
      CALL READT (iu_ZVAR,0,ZVART,IM*JM*4,ZVART,1)
      call closeunit(iu_ZVAR)

      DO J=JM-1,1,-1
      DO I=1,IM
        ZVART(I,J+1)=ZVART(I,J)
        ZVARX(I,J+1)=ZVARX(I,J)
        ZVARY(I,J+1)=ZVARY(I,J)
        ZWT(I,J+1)=ZWT(I,J)
      END DO
      END DO
      DO L=LS1,LM
        PKS(L)=(PSFMPT*SIG(L)+PTOP)**KAPA
      END DO

C**** define wave number array EK for GWDRAG
C**** EKX is the mean wave number for wave lengths between a 1x1 degree
C**** box and a model grid box weighted by 1/EK; wave_length=root(area)
      EKS=0.
      DO J=2,JM
        EK1=TWOPI/SQRT(DXYV(J))                ! 2pi/grid_box_size
        EK2=EK1*SQRT((360./IM)*DLAT_DG)        ! 2pi/1x1deg_box_size
        EKX=0.
        if(EK2.gt.EK1) EKX=(EK2-EK1)/LOG(EK2/EK1) ! weighted mean
        EKS=EKS+EKX*DXYV(J)
        EK(1,J)=EKX
        EK(2,J)=EKX
      END DO
      EKS=EKS*IM/AREAG
      EK(3:NM,1:JM)=EKS
      WRITE (6,970) (J,EK(1,J),J=2,JM)
      WRITE (6,971) EKS
  970 FORMAT ('0  J,EK:',9X,1P,7(I4,E12.2)/,9(1X,8(I4,E12.2)/))
  971 FORMAT ('   AVG EK: ',4X,E12.2)

      END SUBROUTINE init_GWDRAG


      SUBROUTINE VDIFF (P,U,V,UT,VT,T,DT1)
!@sum VDIFF Vertical Diffusion in stratosphere
!@auth Bob Suozzo/Jean Lerner
!@ver  1.0
C****
C**** Vertical diffusion coefficient depends on wind shear.
C**** Uses TRIDIAG for implicit scheme (MU=1) as in diffuse53.
C**** This version only does diffusion for lowest LDIFM layers.
C****
      USE CONSTANT, only : rgas,grav,twopi,kapa,sha
      USE MODEL_COM, only : im,jm,lm,psfmpt,sig,ptop,ls1
     *     ,sige,mrch
      USE GEOM, only : sini=>siniv,cosi=>cosiv,imaxj,rapvn,rapvs,dxyv
     *     ,kmaxj,idij,idjj,rapj
      USE PBLCOM, only : tsurf=>tsavg,qsurf=>qsavg,usurf=>usavg,
     *     vsurf=>vsavg
      USE DAGCOM, only : ajl,jl_dudtsdif,JL_dTdtsdrg
      USE STRAT, only : defrm,pk
      IMPLICIT NONE
      INTEGER, PARAMETER :: LDIFM=LM
      REAL*8, PARAMETER :: BYRGAS = 1./RGAS
      REAL*8, DIMENSION(IM,JM,LM+1) :: VKEDDY
      REAL*8, DIMENSION(IM,JM,LM) :: RHO, DUT, DVT, DKE
      REAL*8, DIMENSION(0:LDIFM+1) :: UL,VL,TL,PL,RHOL
      REAL*8, DIMENSION(LDIFM) :: AIRM,AM,AL,AU,B,DU,DV,DTEMP,DQ
      REAL*8, DIMENSION(LDIFM+1) :: TE,PLE,RHOE,DPE,DFLX,KMEDGE,KHEDGE
     *     ,LMEDGE,LHEDGE
      REAL*8, PARAMETER :: MU=1.
      REAL*8, INTENT(INOUT), DIMENSION(IM,JM,LM) :: U,V,T
      REAL*8, INTENT(IN), DIMENSION(IM,JM,LM) :: UT,VT
      REAL*8, INTENT(INOUT), DIMENSION(IM,JM) :: P
      REAL*8, INTENT(IN) :: DT1
      REAL*8 G2DT,PIJ,TPHYS,ediff
      INTEGER I,J,K,L,IP1,NDT,N

      G2DT=GRAV*GRAV*DT1
C**** Fill in USURF,VSURF at poles (Shouldn't this be done already?)
      DO 60 I=2,IM
      USURF(I,1)=USURF(1,1)*COSI(I)-VSURF(1,1)*SINI(I)
      VSURF(I,1)=VSURF(1,1)*COSI(I)+USURF(1,1)*SINI(I)
      USURF(I,JM)=USURF(1,JM)*COSI(I)+VSURF(1,JM)*SINI(I)
   60 VSURF(I,JM)=VSURF(1,JM)*COSI(I)-USURF(1,JM)*SINI(I)
C**** Calculate RHO(I,J,L)
      DO L=1,LS1-1
      DO J=1,JM
      DO I=1,IMAXJ(J)
        RHO(I,J,L)=   BYRGAS*(P(I,J)*SIG(L)+PTOP)/(T(I,J,L)*PK(I,J,L))
      END DO
      END DO
      END DO

      DO L=LS1,LM
      DO J=1,JM
      DO I=1,IMAXJ(J)
        RHO(I,J,L)=   BYRGAS*(PSFMPT*SIG(L)+PTOP)/(T(I,J,L)*PK(I,J,L))
      END DO
      END DO
      END DO
C**** Fill in T,RHO at poles (again shouldn't this be done already?)
      DO L=1,LM
      DO I=1,IM
        T(I,1,L)=T(1,1,L)
        T(I,JM,L)=T(1,JM,L)
        RHO(I,1,L)=RHO(1,1,L)
        RHO(I,JM,L)=RHO(1,JM,L)
      END DO
      END DO
C**** Get Vertical Diffusion Coefficient for this timestep
      CALL GETVK (U,V,VKEDDY,LDIFM)
C****
C**** U,V Diffusion
C****
      DUT=0 ; DVT=0.
      DO 300 J=2,JM
      I=IM
      DO 300 IP1=1,IM
C**** Surface values are used for F(0)
C**** Note area weighting for four point means
      PIJ=(P(I,J-1)+P(IP1,J-1))*RAPVN(J-1)+(P(I,J)+P(IP1,J))*RAPVS(J)
      PL(0)=(PIJ+PTOP)
      UL(0)=(USURF(I  ,J-1) + USURF(IP1,J-1))*RAPVN(J-1) +
     *      (USURF(I  ,J)   + USURF(IP1,J  ))*RAPVS(J)
      VL(0)=(VSURF(I  ,J-1) + VSURF(IP1,J-1))*RAPVN(J-1) +
     *      (VSURF(I  ,J  ) + VSURF(IP1,J  ))*RAPVS(J)
      TPHYS=(TSURF(I  ,J-1) + TSURF(IP1,J-1))*RAPVN(J-1) +
     *      (TSURF(I  ,J)   + TSURF(IP1,J  ))*RAPVS(J)
      TL(0)=TPHYS*PL(0)**KAPA
      RHOL(0)=PL(0)/(RGAS*TPHYS)
      DO L=1,MIN(LDIFM+1,LM)
        IF (L.GE.LS1) PIJ=PSFMPT
        PL(L)=PIJ*SIG(L)+PTOP
        UL(L)=U(I,J,L)
        VL(L)=V(I,J,L)
        RHOL(L)=(RHO(I,J-1,L)+RHO(IP1,J-1,L))*RAPVN(J-1)+
     *          (RHO(I,J  ,L)+RHO(IP1,J  ,L))*RAPVS(J)
        PLE(L)=PIJ*SIGE(L)+PTOP
      END DO
C**** Edge values at LM+1 don't matter since diffusiv flx=F(L)-F(L-1)=0
      IF (L.EQ.LM+1) THEN
         PL(L)   =PIJ*SIGE(L)+PTOP
         UL(L)   =UL(L-1)
         VL(L)   =VL(L-1)
         RHOL(L) =RHOL(L-1)
         PLE(L)  =PIJ*SIGE(L)+PTOP
       ENDIF
      DO L=1,LDIFM
        AIRM(L)  =PLE(L)-PLE(L+1)
      END DO
      DO L=1,LDIFM+1
        KMEDGE(L)=VKEDDY(I,J,L)
        KHEDGE(L)=VKEDDY(I,J,L)
        DPE(L)   =PL(L-1)-PL(L)
      END DO
C**** RHOE is obtained from average of T (=P/RHO)
      DO L=1,LDIFM+1
        RHOE(L)=2.*PLE(L)*RHOL(L-1)*RHOL(L)/
     *       (PL(L-1)*RHOL(L)+PL(L)*RHOL(L-1))
      END DO
C**** Calculate diffusive flux and number of timesteps NDT
      DO L=1,LDIFM+1
        DFLX(L)=G2DT*RHOE(L)**2*KMEDGE(L)/DPE(L)
      END DO
      NDT=1

      DO N=1,NDT
C**** dq/dt by diffusion as tridiagonal matrix
        CALL DFUSEQ(AIRM,DFLX,UL(0),MU,AM,AL,AU,B,LDIFM)
        CALL TRIDIAG(Al,Am,AU,B,DU,LDIFM)
        CALL DFUSEQ(AIRM,DFLX,VL(0),MU,AM,AL,AU,B,LDIFM)
        CALL TRIDIAG(Al,Am,AU,B,DV,LDIFM)
C**** Update model winds
        IF (MRCH.GT.0) THEN
          DO L=1,LM
            AJL(J,L,JL_DUDTSDIF) = AJL(J,L,JL_DUDTSDIF) + DU(L)
          END DO
        ENDIF
        DO L=1,LM
          DUT(I,J,L) = DUT(I,J,L) + DU(L)*AIRM(L)*DXYV(J)
          DVT(I,J,L) = DVT(I,J,L) + DV(L)*AIRM(L)*DXYV(J)
          DKE(I,J,L) = DU(L)*(U(I,J,L)+0.5*DU(L))+
     *                 DV(L)*(V(I,J,L)+0.5*DV(L))
          U(I,J,L) = U(I,J,L) + DU(L)
          V(I,J,L) = V(I,J,L) + DV(L)
        END DO
      END DO
  300 I=IP1
C**** conservation diagnostic
      IF (MRCH.gt.0) CALL DIAGCD (6,UT,VT,DUT,DVT,DT1)

C**** PUT THE KINETIC ENERGY BACK IN AS HEAT
!$OMP  PARALLEL DO PRIVATE(I,J,L,K,ediff)
        DO L=1,LM
          DO J=1,JM
            DO I=1,IMAXJ(J)
              ediff=0.
              DO K=1,KMAXJ(J)   ! loop over surrounding vel points
                ediff=ediff+DKE(IDIJ(K,I,J),IDJJ(K,J),L)*RAPJ(K,J)
              END DO
              ediff=ediff/(SHA*PK(I,J,L))
              T(I,J,L)=T(I,J,L)-ediff
              AJL(J,L,JL_dTdtsdrg)=AJL(J,L,JL_dTdtsdrg)-ediff
            END DO
          END DO
        END DO
!$OMP  END PARALLEL DO
      RETURN
C****
      END SUBROUTINE VDIFF

      SUBROUTINE DFUSEQ(AIRM,DFLX,F,MU,AM,AL,AU,B,LM)
!@sum  DFUSEQ calculate tridiagonal terms
!@auth Bob Suozzo/Jean Lerner
!@ver  1.0
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LM
      REAL*8, INTENT(IN), DIMENSION(LM) :: AIRM
      REAL*8, INTENT(IN), DIMENSION(0:LM+1) :: F
      REAL*8, INTENT(IN), DIMENSION(LM+1) :: DFLX
      REAL*8, INTENT(OUT), DIMENSION(LM) :: AM,AL,AU,B
      REAL*8, INTENT(IN) :: MU
      INTEGER L
      REAL*8 BYAIRM

      DO L=1,LM
        BYAIRM = 1./AIRM(L)
        B(L)=     BYAIRM*(DFLX(L+1)*(F(L+1)-F(L))
     *       - DFLX(L)*(F(L)-F(L-1)))
        AM(L)=1.+ (MU*BYAIRM)*(DFLX(L+1)+DFLX(L))
        if (l.lt.lm) AL(L+1)=-(MU*BYAIRM)*DFLX(L)
        AU(L)=-(MU*BYAIRM)*DFLX(L+1)
      END DO
      RETURN
C****
      END SUBROUTINE DFUSEQ

      SUBROUTINE GETVK (U,V,VKEDDY,LDIFM)
!@sum GETVK calculate vertical diff. coefficient (use high wind shear)
!@auth Bob Suozzo/Jean Lerner
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm
      USE PBLCOM, only : usurf=>usavg,vsurf=>vsavg
      IMPLICIT NONE
!@var Vert. Diffusion coefficent
      REAL*8, INTENT(OUT), DIMENSION(IM,JM,LM+1) :: VKEDDY
      REAL*8, INTENT(IN), DIMENSION(IM,JM,LM) :: U,V
      INTEGER, INTENT(IN) :: LDIFM
      REAL*8, PARAMETER :: XEDDY = 10., DV2MAX = 25.**2
      INTEGER I,J,L,IP1
      REAL*8 US1,VS1,DELV2

C**** Calculate surface winds on velocity grid (rewrite!)
C**** Is this calculated in PBL?
      DO J=2,JM
        I=IM
        DO IP1=1,IM
          US1=.25*(USURF(I,J-1)+USURF(IP1,J-1)+USURF(I,J)+USURF(IP1,J))
          VS1=.25*(VSURF(I,J-1)+VSURF(IP1,J-1)+VSURF(I,J)+VSURF(IP1,J))
          DELV2   =(U(I,J,1)-US1)**2 + (V(I,J,1)-VS1)**2
          VKEDDY(I,J,1)=0.
          IF (DELV2.GT.DV2MAX) VKEDDY(I,J,1)=XEDDY
          I=IP1
        END DO
      END DO

      DO L=2,MIN(LDIFM+1,LM)
      DO J=2,JM
      DO I=1,IM
        DELV2   =(U(I,J,L)-U(I,J,L-1))**2 + (V(I,J,L)-V(I,J,L-1))**2
        VKEDDY(I,J,L)=0.
        IF (DELV2.GT.DV2MAX) VKEDDY(I,J,L)=XEDDY
      END DO
      END DO
      END DO
      IF (LDIFM.EQ.LM) THEN
        DO J=2,JM
          DO I=1,IM
            VKEDDY(I,J,LM+1)=0.
          END DO
        END DO
      ENDIF
      DO J=2,JM
        DO I=1,IM
          VKEDDY(I,J,1)=0.
          VKEDDY(I,J,2)=0.
        END DO
      END DO
      RETURN
C****
      END SUBROUTINE GETVK


      SUBROUTINE GWDRAG (P,U,V,UT,VT,T,SZ,DT1)
!@sum  GWDRAG puts a momentum drag in the stratosphere
!@auth Bob Suozzo/Jean Lerner
!@ver  1.0
C****
C**** GWDRAG is called from DYNAM with arguments:
C****      P = Pressure (mb) at end of timestep
C****      U,V,T = wind and Potential Temp. at beginning of timestep
C****      DT1 = timestep (s)
C****
C**** LD: incident layer - lowest layer at which wave is allowed to
C****   break - set to LM+1 if wave doesn't exist
C**** MU: momentum flux of wave - conserved if the wave doesn't break
C**** MUB: breaking momentum flux - proportional to -(W-CN)**3
C****    MU is the Greek letter mu, for momentum.
C****    MU and MUB must be negative.
C**** CN: wave phase speed (m/s)
C**** W-CN: (local wind in direction of wave) - (phase speed of wave)
C****    calculated as WMC = (U*UR + V*VR - CN)  must be positive!!!
C**** For propagation, we need MUB < MU < 0. When MUB > 0, the wave
C****    cannot propagate, instead, it breaks and deposits all its
C****    momentum.  The level where MUB changes sign is a critical
C****    level.  W-CN also changes sign.
C**** (UR,VR): unit vector in direction of the wave
C**** There is some trickiness in the way (W-CN), MU, and UR,VR
C**** are gotten for MC waves. (UR,VR) is in the direction of the mean
C**** wind in the source region. But adding +/- 10, 20, and 40 m/s means
C**** sometimes W-CN is negative, unless the direction of (UR,VR) is
C**** changed.  7159-7159.9 fixes this so that (UR,VR) agrees with CN.
C**** NM:
C**** 1   Mountain waves
C**** 2   Shear waves
C**** 3-8 Convective waves
C**** 9   Deformation wave
C****
      USE CONSTANT, only : grav,sha,twopi,kapa,rgas
      USE MODEL_COM, only : im,jm,lm,byim,airx,lmc,nidyn,sig,sige
     *     ,dsig,psfmpt,ptop,ls1,mrch,zatmo
      USE STRAT, only : nm,xcdnst,defrm,zvart,zvarx,zvary,zwt,ldef,ldefm
     *     ,lbreak,ld2,lshr,pk,ek,pks, qgwmtn, qgwshr, qgwdef, qgwcnv
     *     ,cmtn,cdef,cmc,pbreaktop,defthresh,pconpen
C**** Do I need to put the common decalaration here also?
      USE GEOM, only : dxyv,bydxyv,fcor,imaxj,ravpn,ravps,rapvn,rapvs
     *     ,kmaxj,rapj,idij,idjj
      USE DAGCOM, only : aij,ajl,ij_gw1,ij_gw2,ij_gw3,ij_gw4,ij_gw5
     *     ,ij_gw6,ij_gw7,ij_gw8,ij_gw9
     &     ,jl_sdifcoef,jl_dtdtsdrg,JL_gwFirst
      IMPLICIT NONE
!@var BVF(LMC1) is Brunt-Vaissala frequency at top of convection
!@var CLDHT is height of cloud = 8000*LOG(P(cloud bottom)/P(cloud top)
      REAL*8, INTENT(INOUT), DIMENSION(IM,JM,LM) :: T,U,V,SZ
      REAL*8, INTENT(IN), DIMENSION(IM,JM,LM) :: UT,VT
      REAL*8, INTENT(IN), DIMENSION(IM,JM) :: P
      REAL*8, INTENT(IN) :: DT1
      REAL*8, PARAMETER :: ERR=1d-20, H0=8000., XFROUD=1.
!@var XLIMIT per timestep limit on mixing and drag
      REAL*8, PARAMETER :: XLIMIT=.1d0
!@var ROTK should this be set from CONSTANT?
      REAL*8, PARAMETER :: ROTK = 1.5, RKBY3= ROTK*ROTK*ROTK

      REAL*8, DIMENSION(IM,JM,LM) :: DUT3,DVT3,DKE,TLS,THLS,BVS
      INTEGER, DIMENSION(IM,2:JM) :: LDRAGA
      REAL*8, DIMENSION(IM,LM) :: UIL,VIL,TLIL,THIL,BVIL
      REAL*8, DIMENSION(LM,JM) :: DUJL
      REAL*8, DIMENSION(LM) :: PL,DP,TL,THL,RHO,BVF,UL,VL,DL,DUT,DVT,
     *     DQT,DTT,RDI,DFTL,DFM,DFR,WMC,UEDGE,VEDGE,BYFACS,CN
      REAL*8 MUB(LM+1,NM),PLE(LM+1),MU(NM),UR(NM),VR(NM),WT(NM)
     *     ,VARXS(IM),VARXN(IM),VARYS(IM),VARYN(IM)
      DATA CN(1)/0./
      REAL*8 :: DTHR,BYDT1
      INTEGER LD(NM)
      INTEGER I,J,L,N,K,LN,LMC0,LMC1,NMX,LDRAG,LD1,LTOP,IP1
      REAL*8 FCORU,PIJ,SP,U0,V0,W0,BV0,ZVAR,P0,DU,DV,DW,CU,USRC,VSRC
     *     ,AIRX4,AIRXS,CLDDEP,FPLUME,CLDHT,WTX,TEDGE,WMCE,BVEDGE,DFMAX
     *     ,EXCESS,ALFA,XDIFF,DFT,DWT,FDEFRM,WSRC,WCHECK,DUTN,PDN
     *     ,YDN,FLUXUD,FLUXVD,PUP,YUP,DX,DLIMIT,FLUXU,FLUXV,MDN
     *     ,MUP,MUR,BVFSQ,ediff
C****
      DTHR=2./NIdyn
      BYDT1 = 1./DT1

C**** Start main loop
!$OMP  PARALLEL DO PRIVATE(I,J,L)
      DO L=1,LM
        IF (L.GE.LS1) THEN
C**** P**KAPA IN THE STRATOSPHERE
          DO J=1,JM
          DO I=1,IMAXJ(J)
            PK(I,J,L)=PKS(L)
          END DO
          END DO
        ELSE
C**** P**KAPA IN THE TROPOSPHERE
          DO J=1,JM
          DO I=1,IMAXJ(J)
            PK(I,J,L)=(P(I,J)*SIG(L)+PTOP)**KAPA
          END DO
          END DO
        END IF
      END DO
!$OMP  END PARALLEL DO
C****
C**** FILL IN QUANTITIES AT POLES
C****
      DO I=2,IM
        AIRX(I,1)=AIRX(1,1)
        LMC(1,I,1)=LMC(1,1,1)
        LMC(2,I,1)=LMC(2,1,1)
        AIRX(I,JM)=AIRX(1,JM)
        LMC(1,I,JM)=LMC(1,1,JM)
        LMC(2,I,JM)=LMC(2,1,JM)
      END DO
      DO L=1,LM
      DO I=2,IM
        T(I,1,L)=T(1,1,L)
        SZ(I,1,L)=SZ(1,1,L)
        PK(I,1,L)=PK(1,1,L)
        T(I,JM,L)=T(1,JM,L)
        SZ(I,JM,L)=SZ(1,JM,L)
        PK(I,JM,L)=PK(1,JM,L)
      END DO
      END DO
C****
C**** DEFORMATION
C****
      IF(MRCH.EQ.0)  CALL DEFORM (P,U,V)
!$OMP  PARALLEL DO PRIVATE(L)
      DO L=1,LM
        DUT3(:,:,L)=0. ; DVT3(:,:,L)=0.
      END DO
!$OMP  END PARALLEL DO
!$OMP  PARALLEL DO PRIVATE(I,IP1,J,L)
      DO L=1,LM
      DO J=2,JM
      I=IM
      DO IP1=1,IM
         TLS(I,J,L) =
     *    .25*(PK(I,J-1,L)*T(I,J-1,L)+PK(IP1,J-1,L)*T(IP1,J-1,L)+
     *     PK(I,J,L)*T(I,J,L)+PK(IP1,J,L)*T(IP1,J,L))
         THLS(I,J,L) =
     *    .25*(T(I,J-1,L)+T(IP1,J-1,L)+T(I,J,L)+T(IP1,J,L))
         BVS(I,J,L) =
     *    .5*(SZ(I,J-1,L)+SZ(IP1,J-1,L)+SZ(I,J,L)+SZ(IP1,J,L))
         I=IP1
       END DO
       END DO
       END DO
!$OMP  END PARALLEL DO
C****
C**** BEGINNING OF OUTER LOOP OVER I,J
C****
!$OMP  PARALLEL DO PRIVATE(ALFA,AIRX4,BVF,BVFSQ,BV0,BVEDGE,BYFACS,
!$OMP*  CN,CU,DP,DL,DUT,DVT,DU,DV,DW,DFM,DFR,DFT,DFMAX,DFTL,DWT,DUTN,DX,
!$OMP*  DLIMIT,EXCESS,FCORU,FDEFRM,FLUXU,FLUXV,FLUXUD,FLUXVD,
!$OMP*  I,IP1,J,L,LD,LD1,LMC0,LMC1,LN,LDRAG,LTOP, MU,MUB,MDN,MUP,MUR,
!$OMP*  N,NMX, PIJ,PLE,PL,P0,PDN,PUP, RHO,SP,TL,THL,TEDGE,
!$OMP*  UL,U0,UR,USRC,UEDGE, VL,V0,VR,VSRC,VEDGE, W0,WSRC,WCHECK,
!$OMP*  WMCE,WMC, XDIFF, YDN,YUP, ZVAR,WTX,WT,AIRXS,CLDDEP,CLDHT,FPLUME,
!$OMP*  UIL,VIL, TLIL,THIL,BVIL)
      DO J=2,JM
C**** parallel reductions
         UIL(:,:)=U(:,J,:)
         VIL(:,:)=V(:,J,:)
         TLIL(:,:)=TLS(:,J,:)
         THIL(:,:)=THLS(:,J,:)
         BVIL(:,:)=BVS(:,J,:)

      CN(1)=0.
      FCORU=(ABS(FCOR(J-1))+ABS(FCOR(J)))*BYDXYV(J)
      I=IM
      DO IP1=1,IM
C****
C**** CALCULATE VERTICAL ARRAYS
C****
      PIJ=(P(I,J-1)+P(IP1,J-1))*RAPVN(J-1)+(P(I,J)+P(IP1,J))*RAPVS(J)
      SP=PIJ
      DO L=1,LM
      IF (L.GE.LS1) PIJ=PSFMPT
c      TL(L)=.25*(PK(I,J-1,L)*T(I,J-1,L)+PK(IP1,J-1,L)*T(IP1,J-1,L)+
c     *  PK(I,J,L)*T(I,J,L)+PK(IP1,J,L)*T(IP1,J,L))
c      THL(L)=.25*(T(I,J-1,L)+T(IP1,J-1,L)+T(I,J,L)+T(IP1,J,L))
      TL(L)=TLIL(I,L)
      THL(L)=THIL(I,L)
      PLE(L)=PIJ*SIGE(L)+PTOP
      PL(L)=PIJ*SIG(L)+PTOP
      DP(L)=PIJ*DSIG(L)
      RHO(L)=PL(L)/(RGAS*TL(L))
c      BVFSQ=.5*(SZ(I,J-1,L)+SZ(IP1,J-1,L)+SZ(I,J,L)+SZ(IP1,J,L))/
c     *  (DP(L)*THL(L))*GRAV*GRAV*RHO(L)
      BVFSQ=BVIL(I,L)/(DP(L)*THL(L))*GRAV*GRAV*RHO(L)
      IF (PL(L).GE..4d0) THEN
         BVF(L)=SQRT(MAX(BVFSQ,1.d-10))
        ELSE
         BVF(L)=SQRT(MAX(BVFSQ,1.d-4))
      END IF
CRAD  RDI(L)=960.*960./(TL(L)*TL(L))*EXP(-960./TL(L))
      DL(L)=0.
      DUT(L)=0.
      DVT(L)=0.
      UL(L)=UIL(I,L)   ! U(I,J,L)
      VL(L)=VIL(I,L)   ! V(I,J,L)
      END DO
      PLE(LM+1)=PIJ*SIGE(LM+1)+PTOP
      PIJ=SP
C****
C**** INITIALIZE THE MOMENTUM FLUX FOR VARIOUS WAVES
C****
      MU(:)=0.
      LD(:)=LM+1                ! set drag level = LM+1 (turns off drag)
      WT(:)=1.                  ! initialize area weight to 1
      MUB(:,:)=0.
C**** MOUNTAIN WAVES generate at 1 s.d. above topography ...
      IF (QGWMTN.eq.1) THEN
        LD(1)=LBREAK
        U0=(UL(1)*DSIG(1)+UL(2)*DSIG(2))/(SIGE(1)-SIGE(3))
        V0=(VL(1)*DSIG(1)+VL(2)*DSIG(2))/(SIGE(1)-SIGE(3))
        W0=SQRT(U0*U0+V0*V0)
        UR(1)=U0/(W0+ ERR)
        VR(1)=V0/(W0+ ERR)
        BV0=(BVF(1)*DSIG(1)+BVF(2)*DSIG(2))/(SIGE(1)-SIGE(3))
        ZVAR=ABS(UR(1))*ZVARX(I,J)+ABS(VR(1))*ZVARY(I,J)
        IF(ZVAR.LT.XCDNST(1)*XCDNST(1)) ZVAR=0.
C.... if Froude number (U0/BV0*ZSD) > 1
C.... limit ZSD to be consistent with Froude no. (U0/BV0*ZSD) > 1
        IF (ZVAR.GT.(XFROUD*W0/BV0)**2) ZVAR=(XFROUD*W0/BV0)**2
        P0=(SIG(1)*DSIG(1)+SIG(2)*DSIG(2))/(SIGE(1)-SIGE(3))*PIJ+PTOP
        WT(1)=ZWT(I,J)
        MU(1)=-CMTN*EK(1,J)/(H0*ROTK)*P0*BV0*W0*ZVAR
        IF(MU(1)*(UL(LBREAK)*UR(1)+VL(LBREAK)*VR(1)).GE.0.) MU(1)=0.
      END IF
C**** DEFORMATION WAVE (X    d cm-2)
      IF (QGWDEF.eq.1) THEN
        IF (ZATMO(I,J)/GRAV.GT.1000.) GO TO 155
        IF (DEFRM(I,J).LT. DEFTHRESH) GO TO 155 !  deform. Threshold
        FDEFRM=- CDEF*GRAV/(1000.*DEFTHRESH)*DEFRM(I,J)
        DU=UL(LDEF +1)-UL(LDEF )
        DV=VL(LDEF +1)-VL(LDEF )
        DW=SQRT(DU**2        + DV**2       )
        UR(9)=DU       /(DW+ ERR)
        VR(9)=DV       /(DW+ ERR)
        CN(9)=UL(LDEF )*UR(9)+VL(LDEF )*VR(9)
        MU(9)=FDEFRM
        LD(9)=LBREAK
      END IF
C**** WIND SHEAR: USE SHEAR BETWEEN 7 AND 8 UNLESS CRIT. LEVEL ABOVE..
 155  IF (QGWSHR.eq.1) THEN
        LN=LSHR
 160    L=LN
        CU=.5*(UL(L)+UL(L+1))
 170    LN=LN+1
        IF (LN.GE.LD2) GO TO 172
        IF ((UL(LN)-CU)*(UL(LN+1)-CU).LT.0.) GO TO 160
        GO TO 170
 172    LD(2)=L+2
        DU=UL(L+1)-UL(L)
        DV=VL(L+1)-VL(L)
        DW=SQRT(DU*DU+DV*DV)
        UR(2)=DU/(DW+ ERR)
        VR(2)=DV/(DW+ ERR)
        CN(2)=.5*((UL(L+1)+UL(L))*UR(2)+(VL(L+1)+VL(L))*VR(2))
        MU(2)=-FCORU*PLE(L+1)*DW*DW/(240.*H0*BVF(L+1))
      END IF
C**** MOIST CONVECTIVE MASS FLUX BEGINS TWO LEVELS ABOVE CLOUD...
C**** AMPLITUDE DEPENDS ON |U(SOURCE)-C|.
C**** NOTE:  NM LE 2, NM EQ 4, NM GE 8  ARE ALLOWED FOR MC DRAG
      IF (QGWCNV.eq.1) THEN
        USRC=0.
        VSRC=0.
        AIRX4=AIRX(I,J-1)+AIRX(IP1,J-1)+AIRX(I,J)+AIRX(IP1,J)
        AIRXS= ((AIRX(I,J-1)+AIRX(IP1,J-1))*RAVPN(J-1)
     *       +  (AIRX(I,J)+AIRX(IP1,J))*RAVPS(J))
        IF (AIRX4.LE.0.) GO TO 200
        LMC0=.5+(LMC(1,I,J-1)*AIRX(I,J-1)+LMC(1,IP1,J-1)*AIRX(IP1,J-1)+
     *       LMC(1,I,J)*AIRX(I,J)+LMC(1,IP1,J)*AIRX(IP1,J))/(AIRX4+ ERR)
C**** Note: LMC1 was defined in CB245M31 as LMAX+1
        LMC1=.5+(LMC(2,I,J-1)*AIRX(I,J-1)+LMC(2,IP1,J-1)*AIRX(IP1,J-1)+
     *       LMC(2,I,J)*AIRX(I,J)+LMC(2,IP1,J)*AIRX(IP1,J))/(AIRX4+ ERR)
C**** Do nothing for shallow convection (below 800 mb)
        IF (PLE(LMC1-1).gt.800.) GO TO 200
        NMX=4
CCC     CLDDEP=PIJ*(SIGE(LMC0)-SIGE(LMC1))
CCC     FPLUME=AIRXS/(DXYV(J)*CLDDEP)
CCC     CLDHT=H0*LOG((PIJ*SIGE(LMC0)+PTOP)/(PIJ*SIGE(LMC1)+PTOP))
        IF(LMC1.LT.LS1)  THEN
           CLDDEP=PIJ*(SIGE(LMC0)-SIGE(LMC1))
           FPLUME=AIRXS/(DXYV(J)*CLDDEP)
           CLDHT=H0*LOG((PIJ*SIGE(LMC0)+PTOP)/(PIJ*SIGE(LMC1)+PTOP))
         ELSE
           CLDDEP=PIJ*SIGE(LMC0) - PSFMPT*SIGE(LMC1)
           FPLUME=AIRXS/(DXYV(J)*CLDDEP)
           CLDHT=H0*LOG((PIJ*SIGE(LMC0)+PTOP)/(PSFMPT*SIGE(LMC1)+PTOP))
        END IF
        WTX=FPLUME
        IF (WTX.GT.1. OR. WTX.LT.0.) THEN
          PRINT *, 'WARNING IN GWDRAG, WTX INCORRECT',WTX,FPLUME,CLDDEP
          IF (WTX.GT.1.) WTX=1.
          IF (WTX.LT.0.) call stop_model(' WTX <0 IN GWDRAG',255)
        END IF
        DO L=LMC0,LMC1-1
          USRC=USRC+UL(L)
          VSRC=VSRC+VL(L)
        END DO
        USRC=USRC/(LMC1-LMC0)
        VSRC=VSRC/(LMC1-LMC0)
        WSRC=SQRT(USRC*USRC+VSRC*VSRC)
        UR(3)=USRC/(WSRC+ ERR)
        VR(3)=VSRC/(WSRC+ ERR)
        MU(3)=-EK(3,J)*CMC*BVF(LMC1-1)*PL(LMC1-1)*CLDHT**2
        MU(4)=MU(3)
        CN(3)=WSRC-10.
        CN(4)=WSRC+10.
        UR(4)=UR(3)
        VR(4)=VR(3)
        LD(3)=10
        LD(4)=10
        WT(3)=WTX
        WT(4)=WTX
C**** If convection is penetrating (i.e. above PCONPEN) do second set
        IF (PLE(LMC1).LT.PCONPEN .AND. NM.GE.8) THEN
          NMX=8
          DO N=3,NMX
            WT(N)=WTX
            LD(N)=LMC1+1
          END DO
          CN(5)=WSRC-40.
          CN(6)=WSRC+40.
          CN(7)=WSRC-20.
          CN(8)=WSRC+20.
          DO N=5,NMX
            MU(N)=MU(3)
            UR(N)=UR(3)
            VR(N)=VR(3)
          END DO
        ENDIF
        WCHECK=UL(LD(3))*UR(3)+VL(LD(3))*VR(3)
        DO N=3,NMX
          IF (WCHECK.GT.CN(N)) CYCLE
          UR(N)=-UR(N)
          VR(N)=-VR(N)
          CN(N)=-CN(N)
        END DO
      END IF
C****
C**** BREAKING MOMENTUM FLUX AT LAYER EDGES
C****
  200 CONTINUE
      DO L=2,LM
        UEDGE(L)=.5*(UL(L-1)+UL(L))
        VEDGE(L)=.5*(VL(L-1)+VL(L))
        TEDGE=.5*(TL(L-1)+TL(L))
        BVEDGE=.5*(BVF(L-1)+BVF(L))
        BYFACS(L)=-.5*GRAV*PLE(L)/(RGAS*RKBY3*BVEDGE*TEDGE)
      END DO
      DO N=1,NM
        DO L=LD(N),LM
          WMCE=UEDGE(L)*UR(N)+VEDGE(L)*VR(N) - CN(N)
          IF (WMCE.GE.0.) THEN
            MUB(L,N)=BYFACS(L)*EK(N,J)*WMCE**3
            IF (MUB(L,N).GE.-ERR) MUB(L,N)=0.
          END IF
        END DO
      END DO
C**** IF TOP LEVEL IS ABOVE PBREAKTOP FORCE BREAKING GW
      IF (PLE(LM).lt.PBREAKTOP) THEN
C**** DEPOSIT REMAINING MOM. FLUX IN TOP LAYER
        MUB(LM+1,:)=0.
C**** DISTRIBUTE CRIT LEVEL NEAR TOP
        DO N=1,NM
          IF (MUB(LM,N).EQ.0.) MUB(LM,N)=.3d0*MUB(LM-1,N)
        END DO
      END IF
Cw*** APPLY AREA WEIGHTING OF MTN WAVE TO BREAKING MOM. FLUX
Cw    DO 225 L=LD(1),LM
Cw225 MUB(L,1)=MUB(L,1)*ZWT(I,J)
Cw*** APPLY AREA WEIGHTING OF MC WAVE TO BREAKING MOM. FLUX
Cw    DO 230 N=3,NMX
Cw    DO 230 L=LD(N),LM
Cw    MUB(L,N)=MUB(L,N)*WT(N)
Cw230 CONTINUE
C****
C**** DETERMINE INCIDENT MOMENTUM FLUX
C****
C**** INCIDENT FLUX FOR MTN WAVES
      LN=LD(1)
      IF (MU(1).LT.MUB(LN,1)) MU(1)=MUB(LN,1)
      IF (LN.LE.LM.AND.BVF(LN).LT.1.E-5) MU(1)=0.
C**** INCIDENT FLUX FOR SHR AND MC WAVES
      DO N=2,NM
        LN=LD(N)
        IF (LN.GT.LM) CYCLE
        IF (MU(N).LT.MUB(LN,N)) MU(N)=MUB(LN,N)
      END DO
C**** LOCATE LDRAG - LOWEST LAYER TO APPLY DRAG
      LDRAG=LM+2    ! better default
      DO N=1,NM
        IF (MU(N).LE.-ERR) LDRAG=MIN(LD(N),LDRAG)
      END DO
C**** DIAGNOSTICS
         IF (MRCH.EQ.2) THEN
         AIJ(I,J,IJ_GW1)=AIJ(I,J,IJ_GW1)+MU(9)*UR(9)*DTHR
         AIJ(I,J,IJ_GW2)=AIJ(I,J,IJ_GW2)+MU(1)*UR(1)*DTHR  *WT(1)
         AIJ(I,J,IJ_GW3)=AIJ(I,J,IJ_GW3)+MU(2)*UR(2)*DTHR
         AIJ(I,J,IJ_GW4)=AIJ(I,J,IJ_GW4)+MU(3)*UR(3)*DTHR  *WT(3)
         AIJ(I,J,IJ_GW5)=AIJ(I,J,IJ_GW5)+MU(7)*UR(7)*DTHR  *WT(7)
         AIJ(I,J,IJ_GW6)=AIJ(I,J,IJ_GW6)+MU(5)*UR(5)*DTHR  *WT(5)
         AIJ(I,J,IJ_GW7)=AIJ(I,J,IJ_GW7)+CN(2)*UR(2)*DTHR
         AIJ(I,J,IJ_GW8)=AIJ(I,J,IJ_GW8)+USRC*DTHR
         END IF
C****
C**** CALCULATE THE DRAG
C****
      LDRAGA(I,J)=LDRAG
      IF (LDRAG.GT.LM) GO TO 500
      DO 400 N=1,NM
      IF (LD(N).GT.LM) GO TO 400
      LD1=LD(N)
      WMC(LD1-1)=UL(LD1-1)*UR(N)+VL(LD1-1)*VR(N)-CN(N)
      DO 260 L=LD1,LM-1
      DFM(L)=0.
      DFR(L)=0.
      WMC(L)=UL(L)*UR(N)+VL(L)*VR(N)-CN(N)
      IF (WMC(L).GT..01) GO TO 255
      WMC(L)=.5*(WMC(L-1)+WMC(L))
      IF (WMC(L).LT..01) WMC(L)=.01d0
  255 LTOP=L
      IF (MUB(L+1,N).EQ.0.) GO TO 270
  260 CONTINUE
      LTOP=L
      DFM(L)=0.
      DFR(L)=0.
      WMC(L)=UL(L)*UR(N)+VL(L)*VR(N)-CN(N)
      IF (WMC(L).GT..01) GO TO 270
      WMC(L)=.5*(WMC(L-1)+WMC(L))
      IF (WMC(L).LT..01) WMC(L)=.01d0
  270 CONTINUE
      DO 300 L=LD1,LTOP
      IF (L.EQ.LM.AND.MRCH.EQ.2)
     *   AIJ(I,J,IJ_GW9)=AIJ(I,J,IJ_GW9)+MU(N)*UR(N)*DTHR  *WT(N)
C**** RADIATIVE DAMPING
      MUR=MU(N)
CRAD  DTW=DP(L)*BVF(L)/(EK(N,J)*WMC(L)*WMC(L)*RHO(L)*GRAV+ ERR)
CRAD  WAVEM=1.E3*BVF(L)/(ABS(WMC(L))+ ERR)
CRAD  IF (WAVEM.GT.1.6) WAVEM=1.6d0
CRAD  WM3BY2=WAVEM*SQRT(WAVEM)
CRAD  DTR=86400./(RDI(L)*(AZ(L)+BZ(L)*WAVEM*WAVEM/(CZ(L)+WM3BY2))+ERR)
CRAD  MUR=MU(N)*EXP(-2.*DTW/(DTR+ ERR))
      DFR(L)=MU(N)-MUR
C**** MECHANICAL (TURBULENT) DRAG
      MU(N)=MUR
      IF (N.NE.1.AND.N.NE.9) THEN    ! not mtn or deformation
      IF (MUR.LT.MUB(L+1,N)) MU(N)=MUB(L+1,N)  ! saturation drag
      ELSE
         IF (MUR.LT.MUB(L+1,N)) MU(N)=0.     ! sensitivity test
      ENDIF
      DFM(L)=MUR-MU(N)
  300 CONTINUE
C**** LIMIT THE CONVERGENCE TO XLIMIT*(U-C)
  320 DO L=LTOP,LD1,-1
      DFT=DFR(L)+DFM(L)
      DFMAX=-XLIMIT*WMC(L)*DP(L)*BYDT1
      IF (DFT.GE.DFMAX) CYCLE
      EXCESS=DFT-DFMAX
      ALFA=DFM(L)/DFT
      DFR(L)=DFR(L)-EXCESS*(1.-ALFA)
      DFM(L)=DFM(L)-EXCESS*ALFA
      IF (L.GT.LD1) THEN
        DFR(L-1)=DFR(L-1)+EXCESS*(1.-ALFA)
        DFM(L-1)=DFM(L-1)+EXCESS*ALFA
      ENDIF
      END DO
C**** COMBINE MECHANICAL AND RADIATIVE DRAG
      DO L=LD1,LTOP
        DFTL(L)=DFM(L)+DFR(L)
      END DO
C**** ACCUMULATE SATURATION FLUX FOR MTN WAVES DURING BRKING
C     IF (N.NE.1.OR.MRCH.NE.2) GO TO 375
C     DO 370 L=LD1,LTOP
C     IF (DFM(L).EQ.0.) GO TO 370
C        AJL(J,L,JL_30)=AJL(J,L,JL_30)+ZWT(I,J)
C        AJL(J,L,JL_DUDFMDRG)=
C    *     AJL(J,L,JL_DUDFMDRG)+MUB(L+1,1)*UR(1)*ZWT(I,J)
C 370 CONTINUE
C**** CALCULATE DIFFUSION COEFFICIENT AND ADD DRAG TO THE WINDS
C**** (DEFORMATION DIFFUSION IS * XDIFF)
      XDIFF=1.
      DO L=LD1,LTOP
        DL(L)=DL(L)+ABS(DFM(L)*WMC(L))*XDIFF*WT(N)
        DFT=DFTL(L)
        DWT=DFT*DT1/DP(L)
        IF (DFT.GT.-ERR) CYCLE
C**** SHEAR AND MOUNTAIN ORIENTATION
        DUTN=DWT*UR(N)*WT(N)
        DUT(L)=DUT(L)+DUTN
        DVT(L)=DVT(L)+DWT*VR(N)
        IF (MRCH.EQ.2) AJL(J,L,N+JL_gwFirst-1)=AJL(J,L,N+JL_gwFirst-1)
     *       +DUTN
      END DO
  400 CONTINUE
C****
C**** MOMENTUM DIFFUSION   (DOUBLED)
C****    (limited to XLIMIT per timestep)
      PDN=PL(LDRAG-1)
      MDN=DP(LDRAG-1)
      YDN=0.
      FLUXUD=0.
      FLUXVD=0.
      DO L=LDRAG,LM
      PUP=PL(L)
      MUP=DP(L)
      YUP=GRAV*GRAV*DT1*RHO(L)*RHO(L)*DL(L)/(DP(L)*BVF(L)*BVF(L))
      DX=   (YDN+YUP)/(PDN-PUP)     ! double diffusion coefficient
      DLIMIT=XLIMIT*MIN(MUP,MDN)
      IF (DX.GT.DLIMIT) THEN
C       IF (DX.GT.DLIMIT)
C    *   WRITE (99,9430) MRCH,I,J,L,DX,MUP,MDN,UL(L),UL(L-1)
C9430 FORMAT (' 7267 DIFFX > LIMIT:',4I3,1P,3E10.2,0P,2F7.1)
         DX=DLIMIT
      ENDIF
      FLUXU=DX/(1.+DX*(MUP+MDN)/(MUP*MDN))*(UL(L)-UL(L-1))
      FLUXV=DX/(1.+DX*(MUP+MDN)/(MUP*MDN))*(VL(L)-VL(L-1))
      DUT(L-1)=DUT(L-1)-(FLUXUD-FLUXU)/MDN
      DVT(L-1)=DVT(L-1)-(FLUXVD-FLUXV)/MDN
C        IF (MRCH.EQ.2) AJL(J,L-1,JL_DUDTSDIF)=
C    *      AJL(J,L-1,JL_DUDTSDIF)-(FLUXUD-FLUXU)/MDN
      FLUXUD=FLUXU
      FLUXVD=FLUXV
      YDN=YUP
      MDN=MUP
      PDN=PUP
      END DO
C**** DIFFUSION IN THE TOP LAYER COMES ONLY FROM BELOW.
      DUT(LM)=DUT(LM)-FLUXUD/MDN
      DVT(LM)=DVT(LM)-FLUXVD/MDN

C**** ACCUMULATE DIAGNOSTICS  (DU, DIFFUSION COEFFICIENT)
      IF (MRCH.EQ.2) THEN
C        AJL(J,LM,JL_DUDTSDIF)=AJL(J,LM,JL_DUDTSDIF)-FLUXUD/MDN
        DO L=LDRAG,LM
cc          AJL(J,L,JL_SDIFCOEF)=AJL(J,L,JL_SDIFCOEF)+
cc     &         DL(L)/(BVF(L)*BVF(L))*DTHR
          DUJL(L,J)=DL(L)/(BVF(L)*BVF(L))*DTHR
        END DO
C****
C**** Save KE change and diffusion coefficient on A-grid
C****
        DO L=LDRAG-1,LM
          DKE(I,J,L)=.5*(((DUT(L)+UIL(I,L))**2+(DVT(L)+VIL(I,L))**2)-
     *         (UIL(I,L)**2+VIL(I,L)**2))
          DUT3(I,J,L) = DUT(L)*DP(L)*DXYV(J)
          DVT3(I,J,L) = DVT(L)*DP(L)*DXYV(J)
        END DO
      END IF
C****
C**** UPDATE THE U AND V WINDS
C****
      DO L=LDRAG-1,LM
        UIL(I,L)=UIL(I,L)+DUT(L)
        VIL(I,L)=VIL(I,L)+DVT(L)
      END DO

C**** END OF LOOP OVER I
  500 I=IP1
      END DO

         U(:,J,:)=UIL(:,:)
         V(:,J,:)=VIL(:,:)
C**** END OF LOOP OVER J
      END DO
!$OMP  END PARALLEL DO
C****

      IF (MRCH.EQ.2) THEN
C**** conservation diagnostic
        CALL DIAGCD (6,UT,VT,DUT3,DVT3,DT1)

C**** PUT THE KINETIC ENERGY BACK IN AS HEAT
!$OMP  PARALLEL DO PRIVATE(I,J,L,K,ediff)
        DO L=1,LM
          DO J=1,JM
            DO I=1,IMAXJ(J)
              ediff=0.
              DO K=1,KMAXJ(J)   ! loop over surrounding vel points
                IF (L.ge.LDRAGA(IDIJ(K,I,J),IDJJ(K,J))-1) ediff=ediff
     *               + DKE(IDIJ(K,I,J),IDJJ(K,J),L)*RAPJ(K,J)
              END DO
              ediff=ediff/(SHA*PK(I,J,L))
              T(I,J,L)=T(I,J,L)-ediff
              AJL(J,L,JL_dTdtsdrg)=AJL(J,L,JL_dTdtsdrg)-ediff
            END DO
            AJL(J,L,JL_SDIFCOEF)=AJL(J,L,JL_SDIFCOEF)+ DUJL(L,J)
          END DO
        END DO
!$OMP  END PARALLEL DO

      END IF
C****
      RETURN
  911 FORMAT ('0TOPOGRAPHY VARIANCE    J=',I6,' :')
  912 FORMAT (1X,1P,12E11.1)
  993 FORMAT (1H0,'N=993','I,J,N,LN,MU,MUB=',4I4,2P,2F10.4)
  994 FORMAT (1H0,'N=994  I,J=',2I4,'   L    U  V   MUB(1)  MUB(2)...',
     *  14X,'PLE',6X,'BVF'/(21X,I4,2F8.1,2P,6F11.4,0P,F11.3,F11.4))
  995 FORMAT (' DUT TOO BIG: I,J,L,N,UL,CN,DWTT,MUR,MU,DP=',4I4,3F8.1,
     *  2P,2F8.4,0P,F10.4)
  996 FORMAT (1H0,'MTN WAVE GEN.: I,J,L1,U0,M,U(L1),P0=',
     *  3I4,F7.1,2P,F8.4,0P,F7.1,F11.1)
  997 FORMAT (1X,'SHR WAVE GEN.: I,J,LUP,LDN,DU2K,M,BVF=',
     *  4I4,F7.1,2P,F8.4,1P,2E11.1)
  998 FORMAT (1X,'MC GEN.: I,J,L0,1,HT,W3,M3,C3,US,VS=',
     *  4I4,2X,-3P,F7.2,0P,F7.3,1P,E12.2,1P,F7.2,2X,2F7.2)
  920 FORMAT (1X)
C****
      END SUBROUTINE GWDRAG

      SUBROUTINE DEFORM (P,U,V)
!@sum  DEFORM calculate defomation terms
!@auth Bob Suozzo/Jean Lerner
!@ver  1.0
C****
C**** Deformation terms  DEFRM1=du/dx-dv/dy   DEFRM2=du/dy+dv/dx
C**** For spherical coordinates, we are treating DEFRM1 like DIV
C**** and DEFRM2 like CURL (i.e., like FLUX and CIRCULATION),
C**** except the "V" signs are switched.  DEFRM is RMS on u,v grid
C****
      USE MODEL_COM, only : im,jm,lm,psfmpt,ptop,sig,dsig
      USE DYNAMICS, only : pu,pv
      USE GEOM, only : bydxyv,dxyv,dxv,dyp
      USE STRAT, only : ldef,ldefm,defrm
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,JM,LM), INTENT(INOUT) :: U,V
      REAL*8, DIMENSION(IM,JM), INTENT(INOUT) :: P
      REAL*8, DIMENSION(IM) :: UDXS,DUMS1,DUMS2,DUMN1,DUMN2
      REAL*8, DIMENSION(IM,JM) ::  DEFRM1,DEFRM2,DEF1A,DEF2A
      CHARACTER*80 TITLE
      INTEGER I,J,L,IP1,IM1
      REAL*8 UDXN

C****
C**** Deformation terms  DEFRM1=du/dx-dv/dy   DEFRM2=du/dy+dv/dx
C**** For spherical coordinates, we are treating DEFRM1 like DIV
C**** and DEFRM2 like CURL (i.e., like FLUX and CIRCULATION),
C**** except the "V" signs are switched
C****
      L=LDEF
C**** U-terms
      DO 91 J=1,JM
      IM1=IM
      DO 91 I=1,IM
      DEFRM1(I,J) = PU(I,J,L)-PU(IM1,J,L)
   91 IM1=I
      DO 92 I=1,IM
   92 UDXS(I)=0.
      IM1=IM
      DO 93 J=1,JM-1
      DO 93 I=1,IM
      UDXN=.5*(U(IM1,J+1,L)+U(I,J+1,L))*DXV(J+1)
      DEFRM2(I,J ) = UDXN-UDXS(I)
      UDXS(I)=UDXN
   93 IM1=I
      DO 94 I=1,IM
   94 DEFRM2(I,JM) =     -UDXS(I)
C**** V-terms
      IM1=IM
      DO 98 I=1,IM
      DO 96 J=1,JM-1
   96 DEFRM1(I,J ) = DEFRM1(I,J ) + (PV(I,J+1,L)-PV(I,J ,L))
      DEFRM1(I,1 ) = DEFRM1(I,1 ) + (PV(I,2  ,L)           )
      DEFRM1(I,JM) = DEFRM1(I,JM) + (           -PV(I,JM,L))
      DO 97 J=2,JM-1
   97 DEFRM2(I,J ) = DEFRM2(I,J ) +
     *  .5*((V(I,J,L)+V(I,J+1,L))-(V(IM1,J,L)+V(IM1,J+1,L)))*DYP(J)
      DEFRM2(I,1 ) = DEFRM2(I,1 ) + (V(I,2 ,L)-V(IM1,2 ,L))*DYP(1)
      DEFRM2(I,JM) = DEFRM2(I,JM) + (V(I,JM,L)-V(IM1,JM,L))*DYP(JM)
   98 IM1=I
C**** Convert to UV-grid
      I=IM
      DO 99 IP1=1,IM
      DUMS1(I)=DEFRM1(I,1)+DEFRM1(IP1,1)
      DUMS2(I)=DEFRM2(I,1)+DEFRM2(IP1,1)
   99 I=IP1
      DO 110 J=2,JM
      I=IM
      DO 100 IP1=1,IM
      DUMN1(I)=DEFRM1(I,J)+DEFRM1(IP1,J)
      DUMN2(I)=DEFRM2(I,J)+DEFRM2(IP1,J)
  100 I=IP1
      DO 105 I=1,IM
      DEFRM1(I,J)=DUMS1(I)+DUMN1(I)
  105 DEFRM2(I,J)=DUMS2(I)+DUMN2(I)
      DO 110 I=1,IM
      DUMS1(I)=DUMN1(I)
      DUMS2(I)=DUMN2(I)
  110 CONTINUE
      DO 120 J=2,JM
      I=IM
      DO 120 IP1=1,IM
      DEFRM1(I,J)=DEFRM1(I,J)/
     *   (DXYV(J)*DSIG(L)*(P(I,J)+P(IP1,J)+P(I,J-1)+P(IP1,J-1)))
      DEFRM2(I,J)=.25*DEFRM2(I,J)*BYDXYV(J)
      DEFRM(I,J)=SQRT(DEFRM1(I,J)**2+DEFRM2(I,J)**2)
  120 I=IP1
C**** Set deformation to zero near the poles
      DO 130 I=1,IM
      DEFRM(I,2)=0.
      DEFRM(I,3)=0.
      DEFRM(I,JM-1)=0.
  130 DEFRM(I,JM  )=0.
C****
C**** Print Deformation terms
C****
CWF   IF (MOD(TAU,8.D0).EQ.0.)  RETURN
CW    D1MAX=0.
CW    D2MAX=0.
CWF   DFMAX=0.
CWF   DO 220 J=4,JM-2
CWF   DO 220 I=1,IM
CW    DEF1A(I,J)=ABS(DEFRM1(I,J))
CW    DEF2A(I,J)=ABS(DEFRM2(I,J))
CW    IF(D1MAX.LT.DEF1A(I,J)) D1MAX=DEF1A(I,J)
CW    IF(D2MAX.LT.DEF2A(I,J)) D2MAX=DEF2A(I,J)
CWF   IF(DFMAX.LT.DEFRM(I,J)) DFMAX=DEFRM(I,J)
  220 CONTINUE
CWF   IHR=NINT(TAU)
CW    SCALE=100./10.**(INT(LOG10(D1MAX)))
CW    TITLE='ABS DEFORMATION TERM 1 = DU/DX - DV/DY (XXXXXXXX S-1)'
CW    WRITE (TITLE(41:48),'(1PE8.0)') SCALE
CCW   CALL MAP0(IM,JM,IHR,TITLE,DEF1A,SCALE,0., 0)
CW    SCALE=100./10.**(INT(LOG10(D2MAX)))
CW    TITLE='ABS DEFORMATION TERM 2 = DU/DY - DV/DX (XXXXXXXX S-1)'
CW    WRITE (TITLE(41:48),'(1PE8.0)') SCALE
CCW   CALL MAP0(IM,JM,IHR,TITLE,DEF2A,SCALE,0., 0)
CWF   SCALE=100./10.**(INT(LOG10(DFMAX)))
CWF   TITLE='RMS OF DEFORMATION TERMS 1 AND 2     (XXXXXXXX S-1)'
CWF   WRITE (TITLE(38:45),'(1PE8.0)') SCALE
CWF   CALL MAP0(IM,JM,IHR,TITLE,DEFRM,SCALE,0., 0)
      RETURN
C****
      END SUBROUTINE DEFORM

      SUBROUTINE io_strat(kunit,iaction,ioerr)
!@sum  io_strat reads and writes strat. model variables to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "STRAT01"

      MODULE_HEADER(lhead+1:80) = 'R8: airx(im,jm), I: lmc(2,im,jm)'

      SELECT CASE (IACTION)
      CASE (:IOWRITE) ! output to end-of-month restart file
        WRITE (kunit,err=10) MODULE_HEADER,AIRX,LMC
      CASE (IOREAD:)          ! input from restart file
        READ (kunit,err=10) HEADER,AIRX,LMC
        IF (HEADER(1:lhead).ne.MODULE_HEADER(1:lhead)) THEN
          PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      END SELECT
      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_strat

