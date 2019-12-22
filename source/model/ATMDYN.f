#include "rundeck_opts.h"

      SUBROUTINE DYNAM
!@sum  DYNAM Integrate dynamic terms
!@auth Original development team
!@ver  1.0
      USE CONSTANT, only : by3,sha,mb2kg
      USE MODEL_COM, only : im,jm,lm,u,v,t,p,q,wm,dsig,NIdyn,dt,MODD5K
     *     ,NSTEP,NDA5K,ndaa,mrch,psfmpt,ls1,byim,dt_UVfilter,psf
      USE GEOM, only : dyv,dxv,dxyp,areag
      USE SOMTQ_COM, only : tmom,qmom,mz
      USE DYNAMICS, only : ptold,pu,pv,pit,sd,phi,dut,dvt
     &    ,pua,pva,sda,ps,mb,pk
      USE DAGCOM, only : aij,ij_fmv,ij_fmu,ij_fgzu,ij_fgzv
      IMPLICIT NONE

      REAL*8, DIMENSION(IM,JM) :: PRAT
      REAL*8, DIMENSION(IM,JM,LM) :: UT,VT,TT,TZ,TZT,MA
      REAL*8, DIMENSION(IM,JM,LM) :: UX,VX,PIJL
      REAL*8 PA(IM,JM),PB(IM,JM),PC(IM,JM),FPEU(IM,JM),FPEV(IM,JM)

      REAL*8 DTFS,DTLF,PP,UU,VV
      INTEGER I,J,L,IP1,IM1   !@var I,J,L,IP1,IM1  loop variables
      INTEGER NS, NSOLD,MODDA    !? ,NIdynO

      REAL*8, DIMENSION(JM) :: KEJ,PEJ
      REAL*8 ediff,TE0,TE

!?    NIdynO=MOD(NIdyn,2)   ! NIdyn odd is currently not an option
      DTFS=DT*2./3.
      DTLF=2.*DT
      NS=0
      NSOLD=0                            ! strat
      PTOLD(:,:) = P(:,:)
C**** Initialise total energy (J/m^2)
      call conserv_PE(PEJ)
      call conserv_KE(KEJ)
      TE0=(sum(PEJ(:)*DXYP(:))+sum(KEJ(2:JM)))/AREAG
C**** Initialize mass fluxes used by tracers and Q
      PS (:,:)   = P(:,:)
      PUA(:,:,:) = 0.
      PVA(:,:,:) = 0.
      SDA(:,:,:) = 0.
C**** Leap-frog re-initialization: IF (NS.LT.NIdyn)
  300 CONTINUE
c     UX(:,:,:)  = U(:,:,:)
c     UT(:,:,:)  = U(:,:,:)
c     VX(:,:,:)  = V(:,:,:)
c     VT(:,:,:)  = V(:,:,:)
!     copy z-moment of temperature into contiguous memory
c     tz(:,:,:) = tmom(mz,:,:,:)
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         UX(:,:,L)  = U(:,:,L)
         UT(:,:,L)  = U(:,:,L)
         VX(:,:,L)  = V(:,:,L)
         VT(:,:,L)  = V(:,:,L)
         TZ(:,:,L)  = TMOM(MZ,:,:,L)
      ENDDO
!$OMP  END PARALLEL DO
C
      PA(:,:) = P(:,:)
      PB(:,:) = P(:,:)
      PC(:,:) = P(:,:)
C**** INITIAL FORWARD STEP, QX = Q + .667*DT*F(Q)
      MRCH=0
C     CALL DYNAM (UX,VX,TX,PX,Q,U,V,T,P,Q,DTFS)
      CALL CALC_PIJL(LM,P,PIJL)
      CALL AFLUX (U,V,PIJL)
      CALL ADVECM (P,PB,DTFS)
      CALL GWDRAG (PB,UX,VX,U,V,T,TZ,DTFS)   ! strat
      CALL VDIFF (PB,UX,VX,U,V,T,DTFS)       ! strat
      CALL ADVECV (P,UX,VX,PB,U,V,Pijl,DTFS)  !P->pijl
      CALL PGF (UX,VX,PB,U,V,T,TZ,Pijl,DTFS)
      if (dt_UVfilter.gt.0.) CALL FLTRUV(UX,VX,U,V)
C**** INITIAL BACKWARD STEP IS ODD, QT = Q + DT*F(QX)
      MRCH=-1
C     CALL DYNAM (UT,VT,TT,PT,QT,UX,VX,TX,PX,Q,DT)
      CALL CALC_PIJL(LS1-1,PB,PIJL)
      CALL AFLUX (UX,VX,PIJL)
      CALL ADVECM (P,PA,DT)
      CALL GWDRAG (PA,UT,VT,UX,VX,T,TZ,DT)   ! strat
      CALL VDIFF (PA,UT,VT,UX,VX,T,DT)       ! strat
      CALL ADVECV (P,UT,VT,PA,UX,VX,Pijl,DT)   !PB->pijl
      CALL PGF (UT,VT,PA,UX,VX,T,TZ,Pijl,DT)
      if (dt_UVfilter.gt.0.) CALL FLTRUV(UT,VT,UX,VX)
      GO TO 360
C**** ODD LEAP FROG STEP, QT = QT + 2*DT*F(Q)
  340 MRCH=-2
C     CALL DYNAM (UT,VT,TT,PT,QT,U,V,T,P,Q,DTLF)
      CALL CALC_PIJL(LS1-1,P,PIJL)
      CALL AFLUX (U,V,PIJL)
      CALL ADVECM (PA,PB,DTLF)
      CALL GWDRAG (PB,UT,VT,U,V,T,TZ,DTLF)   ! strat
      CALL VDIFF (PB,UT,VT,U,V,T,DTLF)       ! strat
      CALL ADVECV (PA,UT,VT,PB,U,V,Pijl,DTLF)   !P->pijl
      CALL PGF (UT,VT,PB,U,V,T,TZ,Pijl,DTLF)
      if (dt_UVfilter.gt.0.) CALL FLTRUV(UT,VT,U,V)
      PA(:,:) = PB(:,:)     ! LOAD PB TO PA
C**** EVEN LEAP FROG STEP, Q = Q + 2*DT*F(QT)
  360 NS=NS+2
         MODD5K=MOD(NSTEP+NS-NIdyn+NDA5K*NIdyn+2,NDA5K*NIdyn+2)
      MRCH=2
C     CALL DYNAM (U,V,T,P,Q,UT,VT,TT,PT,QT,DTLF)
      CALL CALC_PIJL(LS1-1,PA,PIJL)
      CALL AFLUX (UT,VT,PIJL)
      CALL ADVECM (PC,P,DTLF)
      CALL GWDRAG (P,U,V,UT,VT,T,TZ,DTLF)   ! strat
      CALL ADVECV (PC,U,V,P,UT,VT,Pijl,DTLF)     !PA->pijl
         MODDA=MOD(NSTEP+NS-NIdyn+NDAA*NIdyn+2,NDAA*NIdyn+2)  ! strat
         IF(MODDA.LT.MRCH) CALL DIAGA0   ! strat
C**** ACCUMULATE MASS FLUXES FOR TRACERS and Q
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
        PUA(:,:,L)=PUA(:,:,L)+PU(:,:,L)
        PVA(:,:,L)=PVA(:,:,L)+PV(:,:,L)
        IF (L.LE.LM-1) SDA(:,:,L)=SDA(:,:,L)+SD(:,:,L)
      END DO
!$OMP  END PARALLEL DO
C**** ADVECT Q AND T
CCC   TT(:,:,:) = T(:,:,:)
CCC   TZT(:,:,:)= TZ(:,:,:)
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         TT(:,:,L)  = T(:,:,L)
         TZT(:,:,L) = TZ(:,:,L)
      ENDDO
!$OMP  END PARALLEL DO
      call calc_amp(pc,ma)
      CALL AADVT (MA,T,TMOM, SD,PU,PV, DTLF,.FALSE.,FPEU,FPEV)
!     save z-moment of temperature in contiguous memory for later
CCC   tz(:,:,:) = tmom(mz,:,:,:)
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         TZ(:,:,L) = TMOM(MZ,:,:,L)
      ENDDO
!$OMP  END PARALLEL DO
      CALL VDIFF (P,U,V,UT,VT,T,DTLF)          ! strat
      PC(:,:)    = .5*(P(:,:)+PC(:,:))
CCC   TT(:,:,:)  = .5*(T(:,:,:)+TT(:,:,:))
CCC   TZT(:,:,:) = .5*(TZ(:,:,:)+TZT(:,:,:))
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         TT(:,:,L)  = .5*(T(:,:,L)+TT(:,:,L))
         TZT(:,:,L) = .5*(TZ(:,:,L)+TZT(:,:,L))
      ENDDO
!$OMP  END PARALLEL DO
         DO L=1,LM
           AIJ(:,2:JM,IJ_FMV)  = AIJ(:,2:JM,IJ_FMV )+PV(:,2:JM,L)*DTLF
           AIJ(:,1,IJ_FMU)  = AIJ(:, 1,IJ_FMU )+PU(:, 1,L)*DTLF*BY3
           AIJ(:,JM,IJ_FMU) = AIJ(:,JM,IJ_FMU )+PU(:,JM,L)*DTLF*BY3
           AIJ(:,2:JM-1,IJ_FMU)=AIJ(:,2:JM-1,IJ_FMU)+PU(:,2:JM-1,L)*DTLF
         END DO
      CALL CALC_PIJL(LS1-1,PC,PIJL)
      CALL PGF (U,V,P,UT,VT,TT,TZT,Pijl,DTLF)    !PC->pijl
         DO L=1,LM
         DO J=2,JM
         IM1=IM
         DO I=1,IM
           PP=.5*(PHI(I,J-1,L)+PHI(I,J,L))
           UU=.5*(U(I,J,L)+U(IM1,J,L))
           AIJ(I,J,IJ_FGZU)=AIJ(I,J,IJ_FGZU)+PP*UU*DYV(J)*DTLF
           VV=.5*(V(I,J,L)+V(IM1,J,L))
           AIJ(I,J,IJ_FGZV)=AIJ(I,J,IJ_FGZV)+PP*VV*DXV(J)*DTLF
         IM1=I
         END DO
         END DO
         END DO
      CALL CALC_AMPK(LS1-1)
      if (dt_UVfilter.gt.0.) CALL FLTRUV(U,V,UT,VT)
      PC(:,:) = P(:,:)      ! LOAD P TO PC
      CALL SDRAG (DTLF)
         IF (MOD(NSTEP+NS-NIdyn+NDAA*NIdyn+2,NDAA*NIdyn+2).LT.MRCH) THEN
           CALL DIAGA
           CALL DIAGB
           CALL EPFLUX (U,V,T,P)
         ENDIF
C**** Restart after 8 steps due to divergence of solutions
      IF (NS-NSOLD.LT.8 .AND. NS.LT.NIdyn) GO TO 340
      NSOLD=NS
      IF (NS.LT.NIdyn) GO TO 300

C**** This fix adjusts thermal energy to conserve total energy TE=KE+PE
C**** Currently energy is put in uniformly weighted by mass
      call conserv_PE(PEJ)
      call conserv_KE(KEJ)
      TE=(sum(PEJ(:)*DXYP(:))+sum(KEJ(2:JM)))/AREAG
      ediff=(TE-TE0)/(PSF*SHA*mb2kg)        ! C
!$OMP  PARALLEL DO PRIVATE (L)
      do l=1,lm
        T(:,:,L)=T(:,:,L)-ediff/PK(L,:,:)
      end do
!$OMP  END PARALLEL DO

C**** Scale WM mixing ratios to conserve liquid water
      PRAT(:,:)=PTOLD(:,:)/P(:,:)
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LS1-1
        WM(:,:,L)=WM(:,:,L)*PRAT(:,:)
      END DO
!$OMP  END PARALLEL DO

      RETURN
      END SUBROUTINE DYNAM


      SUBROUTINE QDYNAM
!@sum  QDYNAM is the driver to integrate dynamic terms by the method
!@+          of pre-computing Courant limits using mean fluxes
!@+    It replaces CALL AADVT (MA,Q,QMOM, SD,PU,PV, DTLF,.TRUE.,
!@auth J. Lerner
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm,q,dt,byim
      USE SOMTQ_COM, only : tmom,qmom
      USE DAGCOM, only: ajl,jl_totntlh,jl_zmfntlh,jl_totvtlh,jl_zmfvtlh
      USE DYNAMICS, only: ps,mb,ma
      USE TRACER_ADV, only:
     *    AADVQ,AADVQ0,sbf,sbm,sfbm,scf,scm,sfcm,ncyc
      IMPLICIT NONE
      REAL*8 DTLF,byncyc,byma
      INTEGER I,J,L   !@var I,J,L loop variables

      DTLF=2.*DT
      CALL CALC_AMP(PS,MB)
      CALL AADVQ0 (DTLF)  ! uses the fluxes pua,pva,sda from DYNAM
C****
C**** convert from concentration to mass units
C****
!$OMP PARALLEL DO PRIVATE (L,J,I)
      DO L=1,LM
      DO J=1,JM
      DO I=1,IM
        Q(I,J,L)=Q(I,J,L)*MB(I,J,L)
        QMOM(:,I,J,L)=QMOM(:,I,J,L)*MB(I,J,L)
      enddo; enddo; enddo
!$OMP END PARALLEL DO
C**** ADVECT
        sfbm = 0.; sbm = 0.; sbf = 0.
        sfcm = 0.; scm = 0.; scf = 0.
      CALL AADVQ (Q,QMOM, .TRUE. ,'q       ')
        byncyc = 1./ncyc
        AJL(:,:,jl_totntlh) = AJL(:,:,jl_totntlh) + sbf(:,:)
        AJL(:,:,jl_zmfntlh) = AJL(:,:,jl_zmfntlh)
     &    + sbm(:,:)*sfbm(:,:)*byim*byncyc
        AJL(:,:,jl_totvtlh) = AJL(:,:,jl_totvtlh) + scf(:,:)
        AJL(:,:,jl_zmfvtlh)  = AJL(:,:,jl_zmfvtlh)
     &    + scm(:,:)*sfcm(:,:)*byim*byncyc
C****
C**** convert from mass to concentration units (using updated MA)
C****
!$OMP PARALLEL DO PRIVATE (L,I,J,BYMA)
      DO L=1,LM
      DO J=1,JM
      DO I=1,IM
        BYMA = 1.D0/MA(I,J,L)
        Q(I,J,L)=Q(I,J,L)*BYMA
        QMOM(:,I,J,L)=QMOM(:,I,J,L)*BYMA
      enddo; enddo; enddo
!$OMP END PARALLEL DO

      RETURN
      END SUBROUTINE QDYNAM


#ifdef TRACERS_ON
      SUBROUTINE TrDYNAM
!@sum  TrDYNAM is the driver to integrate tracer dynamic terms
!@auth J. Lerner
!@ver  1.0
      USE MODEL_COM, only: im,jm,lm,itime,dt,byim
      USE TRACER_COM, only: itime_tr0,trm,trmom,trname,t_qlimit,ntm
      USE TRACER_ADV
      USE TRACER_DIAG_COM, only:
     &  jlnt_nt_tot,jlnt_nt_mm,jlnt_vt_tot,jlnt_vt_mm,TAJLN
      IMPLICIT NONE
      REAL*8 DTLF,byncyc
      INTEGER N

C**** uses the fluxes pua,pva,sda from DYNAM and QDYNAM
      DO N=1,NTM
        IF (itime.LT.itime_tr0(N)) cycle
        sfbm = 0.; sbm = 0.; sbf = 0.
        sfcm = 0.; scm = 0.; scf = 0.
        CALL AADVQ (TRM(1,1,1,n),TrMOM(1,1,1,1,n),t_qlimit(n),trname(n))
        byncyc = 1./ncyc
        TAJLN(:,:,jlnt_nt_tot,n) = TAJLN(:,:,jlnt_nt_tot,n) + sbf(:,:)
        TAJLN(:,:,jlnt_nt_mm, n) = TAJLN(:,:,jlnt_nt_mm, n)
     &    + sbm(:,:)*sfbm(:,:)*byim*byncyc
        TAJLN(:,:,jlnt_vt_tot,n) = TAJLN(:,:,jlnt_vt_tot,n) + scf(:,:)
        TAJLN(:,:,jlnt_vt_mm, n) = TAJLN(:,:,jlnt_vt_mm, n)
     &    + scm(:,:)*sfcm(:,:)*byim*byncyc
      ENDDO
      RETURN
      END SUBROUTINE TrDYNAM
#endif

      SUBROUTINE AFLUX (U,V,PIJL)
!@sum  AFLUX Calculates horizontal/vertical air mass fluxes
!@+    Input: U,V velocities, PIJL pressure
!@+    Output: PIT  pressure tendency (mb m^2/s)
!@+            SD   sigma dot (mb m^2/s)
!@+            PU,PV horizontal mass fluxes (mb m^2/s)
!@+            CONV  horizontal mass convergence (mb m^2/s)
!@+            SPA
!@auth Original development team
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm,ls1,psfmpt,dsig,bydsig,byim
     &     ,zatmo,sige
      USE GEOM
      USE DYNAMICS, only : pit,sd,conv,pu,pv,sd_clouds,spa
      IMPLICIT NONE
C**** CONSTANT PRESSURE AT L=LS1 AND ABOVE, PU,PV CONTAIN DSIG
!@var U,V input velocities (m/s)
!@var PIJL input 3-D pressure field (mb) (no DSIG)
      REAL*8, INTENT(IN), DIMENSION(IM,JM,LM) :: U,V,PIJL
      REAL*8, DIMENSION(IM) :: DUMMYS,DUMMYN
      INTEGER I,J,L,IP1,IM1
      REAL*8 PUS,PUN,PVS,PVN,PBS,PBN,SDNP,SDSP
      REAL*8 WT,DXDSIG,DYDSIG,PVSA(LM),PVNA(LM),xx
C****
C**** BEGINNING OF LAYER LOOP
C****
!$OMP  PARALLEL DO PRIVATE (I,J,L,IM1,IP1,DXDSIG,DYDSIG,DUMMYS,DUMMYN,
!$OMP*                      PUS,PUN,PVS,PVN,PBS,PBN)
      DO 2000 L=1,LM
C****
C**** COMPUTATION OF MASS FLUXES     P,T  PU     PRIMARY GRID ROW
C**** ARAKAWA'S SCHEME B             PV   U,V    SECONDARY GRID ROW
C****
C**** COMPUTE PU, THE WEST-EAST MASS FLUX, AT NON-POLAR POINTS
      DO 2154 J=2,JM-1
      DO 2154 I=1,IM
 2154 SPA(I,J,L)=U(I,J,L)+U(I,J+1,L)
      CALL AVRX (SPA(1,1,L))
      I=IM
      DO 2166 J=2,JM-1
      DYDSIG = 0.25D0*DYP(J)*DSIG(L)
      DO 2165 IP1=1,IM
      PU(I,J,L)=DYDSIG*SPA(I,J,L)*(PIJL(I,J,L)+PIJL(IP1,J,L))
 2165 I=IP1
 2166 CONTINUE
C**** COMPUTE PV, THE SOUTH-NORTH MASS FLUX
      IM1=IM
      DO 2172 J=2,JM
      DXDSIG = 0.25D0*DXV(J)*DSIG(L)
      DO 2170 I=1,IM
      PV(I,J,L)=DXDSIG*(V(I,J,L)+V(IM1,J,L))*(PIJL(I,J,L)+PIJL(I,J-1,L))

 2170 IM1=I
 2172 CONTINUE
C**** COMPUTE PU*3 AT THE POLES
      PUS=0.
      PUN=0.
      PVS=0.
      PVN=0.
      DO 1110 I=1,IM
      PUS=PUS+U(I,2,L)
      PUN=PUN+U(I,JM,L)
      PVS=PVS+PV(I,2,L)
 1110 PVN=PVN+PV(I,JM,L)
      PUS=.25*DYP(2)*PUS*PIJL(1,1,L)*BYIM
      PUN=.25*DYP(JM-1)*PUN*PIJL(1,JM,L)*BYIM
      PVS=PVS*BYIM
      PVN=PVN*BYIM
      PVSA(L)=PVS
      PVNA(L)=PVN
      DUMMYS(1)=0.
      DUMMYN(1)=0.
      DO 1120 I=2,IM
      DUMMYS(I)=DUMMYS(I-1)+(PV(I,2,L) -PVS)*BYDSIG(L)
 1120 DUMMYN(I)=DUMMYN(I-1)+(PV(I,JM,L)-PVN)*BYDSIG(L)
      PBS=0.
      PBN=0.
      DO 1130 I=1,IM
      PBS=PBS+DUMMYS(I)
 1130 PBN=PBN+DUMMYN(I)
      PBS=PBS*BYIM
      PBN=PBN*BYIM
      DO 1140 I=1,IM
      SPA(I,1,L)=4.*(PBS-DUMMYS(I)+PUS)/(DYP(2)*PIJL(1,1,L))
      SPA(I,JM,L)=4.*(DUMMYN(I)-PBN+PUN)/(DYP(JM-1)*PIJL(1,JM,L))
      PU(I,1,L)=3.*(PBS-DUMMYS(I)+PUS)*DSIG(L)
 1140 PU(I,JM,L)=3.*(DUMMYN(I)-PBN+PUN)*DSIG(L)
C****
C**** CONTINUITY EQUATION
C****
C**** COMPUTE CONV, THE HORIZONTAL MASS CONVERGENCE
c     DO 1510 J=2,JM-1
c     IM1=IM
c     DO 1510 I=1,IM
c     CONV(I,J,L)=(PU(IM1,J,L)-PU(I,J,L)+PV(I,J,L)-PV(I,J+1,L))
c1510 IM1=I
c     CONV(1,1,L)=-PVS
c     CONV(1,JM,L)=PVN
 2000 CONTINUE
!$OMP  END PARALLEL DO
C****
C**** END OF HORIZONTAL ADVECTION LAYER LOOP
C****
c
c modify uphill air mass fluxes around steep topography
      do 2015 j=2,jm-1
      i = im
      do 2010 ip1=1,im
         xx = zatmo(ip1,j)-zatmo(i,j)
         if(xx.eq.0.0)  go to 2007
         DO 2005 L=1,LS1-1
ccc         if((zatmo(ip1,j)-zatmo(i,j))*pu(i,j,l).gt.0.) then
            if(xx*pu(i,j,l).gt.0.) then
               if(pu(i,j,l).gt.0.) then
                  wt = (pijl(ip1,j,l)/pijl(i,j,l)-sige(l+1))/dsig(l)
               else
                  wt = (pijl(i,j,l)/pijl(ip1,j,l)-sige(l+1))/dsig(l)
               endif
               if(wt.lt.1.) then
                  if(wt.lt.0.) wt=0.
                  pu(i,j,l+1) = pu(i,j,l+1) + pu(i,j,l)*(1.-wt)
                  pu(i,j,l) = pu(i,j,l)*wt
               else
                  go to 2007
               endif
            endif
 2005    CONTINUE
 2007    CONTINUE
         i = ip1
 2010 CONTINUE
 2015 CONTINUE
ccc   do j=2,jm
      do 2035 j=3,jm-1
      do 2035 i=1,im
         xx = zatmo(i,j)-zatmo(i,j-1)
         if(xx.eq.0.0)  go to 2035
         DO 2020 L=1,LS1-1
ccc         if((zatmo(i,j)-zatmo(i,j-1))*pv(i,j,l).gt.0.) then
            if(xx*pv(i,j,l).gt.0.) then
               if(pv(i,j,l).gt.0.) then
                  wt = (pijl(i,j,l)/pijl(i,j-1,l)-sige(l+1))/dsig(l)
               else
                  wt = (pijl(i,j-1,l)/pijl(i,j,l)-sige(l+1))/dsig(l)
               endif
               if(wt.lt.1.) then
                  if(wt.lt.0.) wt=0.
                  pv(i,j,l+1) = pv(i,j,l+1) + pv(i,j,l)*(1.-wt)
                  pv(i,j,l) = pv(i,j,l)*wt
               else
                  go to 2035
               endif
            endif
 2020    CONTINUE
 2035 CONTINUE
C
C     Now Really Do  CONTINUITY EQUATION
C
C     COMPUTE CONV, THE HORIZONTAL MASS CONVERGENCE
C
!$OMP  PARALLEL DO PRIVATE (I,J,L,IM1)
      DO 2400 L=1,LM
      DO 1510 J=2,JM-1
      IM1=IM
      DO 1510 I=1,IM
      CONV(I,J,L)=(PU(IM1,J,L)-PU(I,J,L)+PV(I,J,L)-PV(I,J+1,L))
 1510 IM1=I
      CONV(1,1,L)=-PVSA(L)
      CONV(1,JM,L)=PVNA(L)
 2400 CONTINUE
!$OMP  END PARALLEL DO
C
C**** COMPUTE PIT, THE PRESSURE TENDENCY
CC    PIT(I,J)=CONV(I,J,1)
C     DO 2420 L=LM,2,-1
C     PIT(1,1)=PIT(1,1)+CONV(1,1,L)
C     PIT(1,JM)=PIT(1,JM)+CONV(1,JM,L)
C     DO 2420 J=2,JM-1
C     DO 2420 I=1,IM
C2420 PIT(I,J)=PIT(I,J)+CONV(I,J,L)
!$OMP  PARALLEL DO PRIVATE(I,J,L)
      DO 2420 J=1,JM
         DO 2410 I=1,IMAXJ(J)
         DO 2410 L=LM-1,1,-1
            PIT(I,J)=PIT(I,J)+SD(I,J,L)
 2410    CONTINUE
 2420 CONTINUE
!$OMP  END PARALLEL DO
C**** COMPUTE SD, SIGMA DOT                        -------
C     SD(1, 1,LM-1)=CONV(1, 1,LM)                     |
C     SD(1,JM,LM-1)=CONV(1,JM,LM)             completely wasteful
C     DO 2430 J=2,JM-1                                |
C     DO 2430 I=1,IM                                  |
C2430 SD(I,J,LM-1)=CONV(I,J,LM)                    -------
C     DO 2435 L=LM-2,LS1-1,-1
C     SD(1, 1,L)=SD(1, 1,L+1)+CONV(1, 1,L+1)
C     SD(1,JM,L)=SD(1,JM,L+1)+CONV(1,JM,L+1)
C     DO 2435 J=2,JM-1
C     DO 2435 I=1,IM
C     SD(I, J,L)=SD(I, J,L+1)+CONV(I, J,L+1)
C2435 CONTINUE
!$OMP  PARALLEL DO PRIVATE(I,J,L)
      DO 2435 J=1,JM
         DO 2430 I=1,IMAXJ(J)
         DO 2430 L=LM-2,LS1-1,-1
            SD(I,J,L)=SD(I,J,L+1)+SD(I,J,L)
 2430    CONTINUE
 2435 CONTINUE
!$OMP  END PARALLEL DO
C     DO 2440 L=LS1-2,1,-1
C     SD(1, 1,L)=SD(1, 1,L+1)+CONV(1, 1,L+1)-DSIG(L+1)*PIT(1, 1)
C     SD(1,JM,L)=SD(1,JM,L+1)+CONV(1,JM,L+1)-DSIG(L+1)*PIT(1,JM)
C     DO 2440 J=2,JM-1
C     DO 2440 I=1,IM
C     SD(I, J,L)=SD(I, J,L+1)+CONV(I, J,L+1)-DSIG(L+1)*PIT(I, J)
C2440 CONTINUE
!$OMP  PARALLEL DO PRIVATE(I,J,L)
      DO 2440 J=1,JM
         DO 2438 I=1,IMAXJ(J)
         DO 2438 L=LS1-2,1,-1
            SD(I,J,L)=SD(I,J,L+1)+SD(I,J,L)-DSIG(L+1)*PIT(I,J)
 2438    CONTINUE
 2440 CONTINUE
!$OMP  END PARALLEL DO
      DO 2450 L=1,LM-1
      DO 2450 I=2,IM
      SD(I,1,L)=SD(1,1,L)
 2450 SD(I,JM,L)=SD(1,JM,L)
C**** temporary fix for CLOUDS module
      SD_CLOUDS(:,:,1)    = PIT
!$OMP PARALLEL DO PRIVATE (L)
      DO L=2,LM
        SD_CLOUDS(:,:,L) = SD(:,:,L-1)
      END DO
!$OMP END PARALLEL DO
C****
      RETURN
      END SUBROUTINE AFLUX


      SUBROUTINE ADVECM (P,PA,DT1)
!@sum  ADVECM Calculates updated column pressures using mass fluxes
!@auth Original development team
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm,mrch,zatmo,u,v,t,q,ptop
      USE GEOM, only : bydxyp,imaxj
      USE DYNAMICS, only : pit
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: P(IM,JM)
      REAL*8, INTENT(OUT) :: PA(IM,JM)
      REAL*8, INTENT(IN) :: DT1
      INTEGER I,J,L,K  !@var I,J,L,K  loop variables

C**** COMPUTE PA, THE NEW SURFACE PRESSURE
      DO J=1,JM
        DO I=1,IMAXJ(J)
          PA(I,J)=P(I,J)+(DT1*PIT(I,J)*BYDXYP(J))
          IF (PA(I,J)+PTOP.GT.1160. .or. PA(I,J)+PTOP.LT.350.) THEN
            WRITE (6,990) I,J,MRCH,P(I,J),PA(I,J),ZATMO(I,J),DT1,
     *           (U(I-1,J,L),U(I,J,L),U(I-1,J+1,L),U(I,J+1,L),
     *            V(I-1,J,L),V(I,J,L),V(I-1,J+1,L),V(I,J+1,L),
     *            T(I,J,L),Q(I,J,L),L=1,LM)
            write(6,*) "Pressure diagnostic error"
            IF (PA(I,J)+PTOP.lt.250. .or. PA(I,J)+PTOP.GT.1200.)
     &           call stop_model('ADVECM: Pressure diagnostic error',11)
          END IF
        END DO
      END DO
      PA(2:IM, 1)=PA(1,1)
      PA(2:IM,JM)=PA(1,JM)
C****
      RETURN
  990 FORMAT (/'0PRESSURE DIAGNOSTIC     I,J,MRCH,P,PA=',3I4,2F10.2/
     *  '     ZATMO=',F10.3,' DT=',F6.1/
     *  '0    U(I-1,J)     U(I,J)   U(I-1,J+1)    U(I,J+1)    V(I-1,J)',
     *   '     V(I,J)   V(I-1,J+1)    V(I,J+1)     T(I,J)     Q(I,J)'/
     *  (1X,9F12.3,F12.6))
      END SUBROUTINE ADVECM


      SUBROUTINE PGF (UT,VT,PB,U,V,T,SZ,P,DT1)
!@sum  PGF Adds pressure gradient forces to momentum
!@auth Original development team
!@ver  1.0
      USE CONSTANT, only : grav,rgas,kapa,bykapa,bykapap1,bykapap2
      USE MODEL_COM, only : im,jm,lm,ls1,mrch,dsig,psfmpt,sige,ptop
     *     ,zatmo,sig,modd5k,bydsig
      USE GEOM, only : imaxj,dxyv,dxv,dyv,dxyp,dyp,dxp
      USE DYNAMICS, only : gz,pu,pit,phi,spa,dut,dvt
      IMPLICIT NONE

      REAL*8, DIMENSION(IM,JM,LM) :: U,V,T
      REAL*8, DIMENSION(IM,JM) :: FD,RFDUX
      REAL*8 UT(IM,JM,LM),VT(IM,JM,LM),TT(IM,JM,LM),
     *  PA(IM,JM),PB(IM,JM),QT(IM,JM,LM),P(IM,JM,LM)
      REAL*8, DIMENSION(IM,JM,1) :: PU0

      REAL*8 PKE(LS1:LM+1)
      REAL*8 SZ(IM,JM,LM),DT4,DT1
      REAL*8 PIJ,PDN,PKDN,PKPDN,PKPPDN,PUP,PKUP,PKPUP,PKPPUP,DP,P0,X
     *     ,BYDP
      REAL*8 TZBYDP,FLUX,FDNP,FDSP,RFDU,PHIDN,FACTOR
      INTEGER I,J,L,IM1,IP1  !@var I,J,IP1,IM1,L loop variab.
C****
      DT4=DT1/4.
      DO L=LS1,LM+1
        PKE(L)=(PSFMPT*SIGE(L)+PTOP)**KAPA
      END DO
C****
C**** VERTICAL DIFFERENCING
C****
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=LS1,LM
      SPA(:,:,L)=0.
      END DO
!$OMP  END PARALLEL DO
!$OMP  PARALLEL DO PRIVATE(I,J,L,DP,P0,PIJ,PHIDN,TZBYDP,X,
!$OMP*             BYDP,PDN,PKDN,PKPDN,PKPPDN,PUP,PKUP,PKPUP,PKPPUP)
      DO J=1,JM
      DO I=1,IMAXJ(J)
        PIJ=P(I,J,1)
        PDN=PIJ+PTOP
        PKDN=PDN**KAPA
        PHIDN=ZATMO(I,J)
C**** LOOP OVER THE LAYERS
        DO L=1,LM
          PKPDN=PKDN*PDN
          PKPPDN=PKPDN*PDN
          IF(L.GE.LS1) THEN
            DP=DSIG(L)*PSFMPT
            BYDP=1./DP
            P0=SIG(L)*PSFMPT+PTOP
            TZBYDP=2.*SZ(I,J,L)*BYDP
            X=T(I,J,L)+TZBYDP*P0
            PUP=SIGE(L+1)*PSFMPT+PTOP
            PKUP=PKE(L+1)
            PKPUP=PKUP*PUP
            PKPPUP=PKPUP*PUP
          ELSE
            DP=DSIG(L)*PIJ
            BYDP=1./DP
            P0=SIG(L)*PIJ+PTOP
            TZBYDP=2.*SZ(I,J,L)*BYDP
            X=T(I,J,L)+TZBYDP*P0
            PUP=SIGE(L+1)*PIJ+PTOP
            PKUP=PUP**KAPA
            PKPUP=PKUP*PUP
            PKPPUP=PKPUP*PUP
C****   CALCULATE SPA, MASS WEIGHTED THROUGHOUT THE LAYER
            SPA(I,J,L)=RGAS*((X+TZBYDP*PTOP)*(PKPDN-PKPUP)*BYKAPAP1
     *      -X*PTOP*(PKDN-PKUP)*BYKAPA-TZBYDP*(PKPPDN-PKPPUP)*BYKAPAP2)
     *      *BYDP
          END IF
C**** CALCULATE PHI, MASS WEIGHTED THROUGHOUT THE LAYER
          PHI(I,J,L)=PHIDN+RGAS*(X*PKDN*BYKAPA-TZBYDP*PKPDN*BYKAPAP1
     *      -(X*(PKPDN-PKPUP)*BYKAPA-TZBYDP*(PKPPDN-PKPPUP)*BYKAPAP2)
     *      *BYDP*BYKAPAP1)
C**** CALULATE PHI AT LAYER TOP (EQUAL TO BOTTOM OF NEXT LAYER)
          PHIDN=PHIDN+RGAS*(X*(PKDN-PKUP)*BYKAPA-TZBYDP*(PKPDN-PKPUP)
     *     *BYKAPAP1)
          PDN=PUP
          PKDN=PKUP
        END DO
      END DO
      END DO
!$OMP END PARALLEL DO
C**** SET POLAR VALUES FROM THOSE AT I=1
      DO L=1,LM
        SPA(2:IM,1,L)=SPA(1,1,L)
        SPA(2:IM,JM,L)=SPA(1,JM,L)
        PHI(2:IM,1,L)=PHI(1,1,L)
        PHI(2:IM,JM,L)=PHI(1,JM,L)
      END DO
!$OMP  PARALLEL DO PRIVATE(L)
      DO L=1,LM
        GZ(:,:,L)=PHI(:,:,L)
      END DO
!$OMP END PARALLEL DO
C****
C**** PRESSURE GRADIENT FORCE
C****
C**** NORTH-SOUTH DERIVATIVE AFFECTS THE V-COMPONENT OF MOMENTUM
C
!$OMP  PARALLEL DO PRIVATE(I,IM1,J,L,FACTOR,FLUX)
      DO 3236 L=1,LM
      DO 3236 J=2,JM
      FACTOR = DT4*DXV(J)*DSIG(L)
      IM1=IM
      DO 3234 I=1,IM
      FLUX=    ((P(I,J,L)+P(I,J-1,L))*(PHI(I,J,L)-PHI(I,J-1,L))+
     *  (SPA(I,J,L)+SPA(I,J-1,L))*(P(I,J,L)-P(I,J-1,L)))*FACTOR
      DVT(I,J,L)  =DVT(I,J,L)  -FLUX
      DVT(IM1,J,L)=DVT(IM1,J,L)-FLUX
 3234 IM1=I
 3236 CONTINUE
!$OMP  END PARALLEL DO
C
C**** SMOOTHED EAST-WEST DERIVATIVE AFFECTS THE U-COMPONENT
C
!$OMP  PARALLEL DO PRIVATE(I,IP1,J,L,FACTOR)
      DO 3300 L=1,LM
      PU(:,1,L)=0.
      PU(:,JM,L)=0.
      I=IM
      DO 3290 J=2,JM-1
      DO 3280 IP1=1,IM
      PU(I,J,L)=(P(IP1,J,L)+P(I,J,L))*(PHI(IP1,J,L)-PHI(I,J,L))+
     *  (SPA(IP1,J,L)+SPA(I,J,L))*(P(IP1,J,L)-P(I,J,L))
 3280 I=IP1
 3290 CONTINUE
      CALL AVRX (PU(1,1,L))
      DO 3294 J=2,JM
      FACTOR = -DT4*DYV(J)*DSIG(L)
      DO 3294 I=1,IM
 3294 DUT(I,J,L)=DUT(I,J,L)+FACTOR*(PU(I,J,L)+PU(I,J-1,L))
 3300 CONTINUE
!$OMP  END PARALLEL DO
C
C**** CALL DIAGNOSTICS
      IF(MRCH.GT.0) THEN
         IF(MODD5K.LT.MRCH) CALL DIAG5D (6,MRCH,DUT,DVT)
         CALL DIAGCD (3,U,V,DUT,DVT,DT1)
      ENDIF
C****
C****
C**** UNDO SCALING PERFORMED AT BEGINNING OF DYNAM
C****
      DO 3410 J=2,JM-1
      DO 3410 I=1,IM
 3410 FD(I,J)=PB(I,J)*DXYP(J)
      FDSP=PB(1, 1)*DXYP( 1)
      FDNP=PB(1,JM)*DXYP(JM)
      FDSP=FDSP+FDSP
      FDNP=FDNP+FDNP
      DO 3520 I=1,IM
      FD(I, 1)=FDSP
 3520 FD(I,JM)=FDNP
C
!$OMP  PARALLEL DO PRIVATE(I,IP1,J)
      DO 3530 J=2,JM
      I=IM
      DO 3525 IP1=1,IM
      RFDUX(I,J)=4./(FD(I,J)+FD(IP1,J)+FD(I,J-1)+FD(IP1,J-1))
 3525 I = IP1
 3530 CONTINUE
!$OMP  END PARALLEL DO
C
!$OMP  PARALLEL DO PRIVATE(I,J,L,RFDU)
      DO 3550 L=1,LM
      DO 3550 J=2,JM
      RFDU=1./(PSFMPT*DXYV(J)*DSIG(L))
      DO 3540 I=1,IM
      IF(L.LT.LS1) RFDU=RFDUX(I,J)*BYDSIG(L)
      VT(I,J,L)=VT(I,J,L)+DVT(I,J,L)*RFDU
      UT(I,J,L)=UT(I,J,L)+DUT(I,J,L)*RFDU
 3540 CONTINUE
 3550 CONTINUE
!$OMP  END PARALLEL DO
C
      RETURN
      END SUBROUTINE PGF

      SUBROUTINE AVRX0
!@sum  AVRX Smoothes zonal mass flux and geopotential near the poles
!@auth Original development team
!@ver  1.0
      USE MODEL_COM, only : im,jm,imh
      USE GEOM, only : dlon,dxp,dyp,bydyp
      USE DYNAMICS, only : xAVRX
C**** THIS VERSION OF AVRX DOES SO BY TRUNCATING THE FOURIER SERIES.
      IMPLICIT NONE
      REAL*8, INTENT(INOUT) :: X(IM,JM)
CCC   REAL*8, SAVE :: SM(IMH,JM),DRAT(JM)
CCC   REAL*8, SAVE, DIMENSION(IMH) :: BYSN
      REAL*8  SM(IMH,JM),DRAT(JM),BYSN(IMH)
      REAL*8, DIMENSION(0:IMH) :: AN,BN
CCC   INTEGER, SAVE :: NMIN(JM)
CCC   INTEGER, SAVE :: IFIRST = 1
      INTEGER NMIN(JM)
      INTEGER J,N
C
      COMMON /AVRXX/ SM,DRAT,NMIN
C
CCC   IF (IFIRST.EQ.1) THEN
CCC   IFIRST=0
C     CALL FFT0(IM)
      DO N=1,IMH
        BYSN(N)=xAVRX/SIN(.5*DLON*N)
      END DO
      DO 50 J=2,JM-1
        DRAT(J) = DXP(J)*BYDYP(3)
        DO 40 N=IMH,1,-1
          SM(N,J) = BYSN(N)*DRAT(J)
          IF(SM(N,J).GT.1.) THEN
            NMIN(J) = N+1
            GO TO 50
          ENDIF
 40     CONTINUE
 50   CONTINUE
CCC   END IF
      RETURN
C****
      ENTRY AVRX (X)
C****
      DO 140 J=2,JM-1
      IF (DRAT(J).GT.1) GO TO 140
      CALL FFT (X(1,J),AN,BN)
      DO N=NMIN(J),IMH-1
        AN(N)=SM(N,J)*AN(N)
        BN(N)=SM(N,J)*BN(N)
      END DO
      AN(IMH) = SM(IMH,J)*AN(IMH)
      CALL FFTI(AN,BN,X(1,J))
  140 CONTINUE
      RETURN
      END SUBROUTINE AVRX0

      SUBROUTINE FILTER
!@sum  FILTER Performs 8-th order shapiro filter in zonal direction
!@auth Original development team
!@ver  1.0
!@calls SHAP1D
C****
C**** MFILTR=1  SMOOTH P USING SEA LEVEL PRESSURE FILTER
C****        2  SMOOTH T USING TROPOSPHERIC STRATIFICATION OF TEMPER
C****        3  SMOOTH P AND T
C****
      USE CONSTANT, only : bbyg,gbyrb,kapa,sha,mb2kg
      USE MODEL_COM, only : im,jm,lm,ls1,t,p,q,wm,mfiltr,zatmo,ptop
     *     ,byim,sig,itime,psf
      USE GEOM, only : areag,dxyp
      USE SOMTQ_COM, only : tmom,qmom
      USE DYNAMICS, only : pk
#ifdef TRACERS_ON
      USE TRACER_COM, only: ntm,trname,ITIME_TR0,trm,trmom
#endif
      USE PBLCOM, only : tsavg
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,JM) :: X,Y
      REAL*8 PSUMO(JM)

      REAL*8 POLD(IM,JM),PRAT(IM,JM)
      REAL*8 PSUMN,PDIF,AKAP,TE0,TE,ediff
      INTEGER I,J,L,N  !@var I,J,L  loop variables
      REAL*8, DIMENSION(JM) :: KEJ,PEJ

      IF (MOD(MFILTR,2).NE.1) GO TO 200
C**** Initialise total energy (J/m^2)
      call conserv_PE(PEJ)
      call conserv_KE(KEJ)
      TE0=(sum(PEJ(:)*DXYP(:))+sum(KEJ(2:JM)))/AREAG
C****
C**** SEA LEVEL PRESSURE FILTER ON P
C****
!$OMP  PARALLEL DO PRIVATE(I,J)
      DO J=2,JM-1
        PSUMO(J)=0.
        DO I=1,IM
          PSUMO(J)=PSUMO(J)+P(I,J)
          POLD(I,J)=P(I,J)      ! Save old pressure
          Y(I,J)=(1.+BBYG*ZATMO(I,J)/TSAVG(I,J))**GBYRB
          X(I,J)=(P(I,J)+PTOP)*Y(I,J)
        END DO
      END DO
!$OMP  END PARALLEL DO
      CALL SHAP1D (8,X)
!$OMP  PARALLEL DO PRIVATE(I,J,PSUMN,PDIF)
      DO J=2,JM-1
        PSUMN=0.
        DO I=1,IM
          P(I,J)=X(I,J)/Y(I,J)-PTOP
C**** reduce large variations (mainly due to topography)
          P(I,J)=MIN(MAX(P(I,J),0.99d0*POLD(I,J)),1.01d0*POLD(I,J))
          PSUMN=PSUMN+P(I,J)
        END DO
        PDIF=(PSUMN-PSUMO(J))*BYIM
        DO I=1,IM
          P(I,J)=P(I,J)-PDIF
        END DO
      END DO
!$OMP  END PARALLEL DO
C**** Scale mixing ratios (incl moments) to conserve mass/heat
!$OMP  PARALLEL DO PRIVATE(I,J)
      DO J=2,JM-1
        DO I=1,IM
          PRAT(I,J)=POLD(I,J)/P(I,J)
        END DO
      END DO
!$OMP  END PARALLEL DO
!$OMP  PARALLEL DO PRIVATE (I,J,L)
      DO L=1,LS1-1
      DO J=2,JM-1
      DO I=1,IM
         Q(I,J,L)= Q(I,J,L)*PRAT(I,J)
         T(I,J,L)= T(I,J,L)*PRAT(I,J)
        WM(I,J,L)=WM(I,J,L)*PRAT(I,J)
        QMOM(:,I,J,L)=QMOM(:,I,J,L)*PRAT(I,J)
        IF (PRAT(I,J).lt.1.) THEN
          TMOM(:,I,J,L)=TMOM(:,I,J,L)*PRAT(I,J)
        END IF
      END DO
      END DO
      END DO
!$OMP  END PARALLEL DO
#ifdef TRACERS_ON
C**** In general, only an air tracer is affected by the filter
C**** This fix conserves tracer concentration, BUT NOT MASS!
C****   It is not wanted for most tracers.
C**** Thus, this code slows model and should be removed if not wanted
C**** Instead of if(trname...) could use n=n_air and avoid the IF-test
C**** But if n_air=0 this will cause problems...
      do n=1,ntm
      if (trname(n).ne.'Air') cycle
!     if (itime.lt.itime_tr0(n)) cycle   !probably not needed
!$OMP  PARALLEL DO PRIVATE (I,J,L)
      DO L=1,LS1-1
        DO J=2,JM-1
          DO I=1,IM
             trm(I,J,L,n)=  trm(I,J,L,n)/PRAT(I,J)
             trmom(:,I,J,L,n)=trmom(:,I,J,L,n)/PRAT(I,J)
      end do; end do; end do
!$OMP  END PARALLEL DO
      end do
#endif
      CALL CALC_AMPK(LS1-1)

C**** This fix adjusts thermal energy to conserve total energy TE=KE+PE
      call conserv_PE(PEJ)
      call conserv_KE(KEJ)
      TE=(sum(PEJ(:)*DXYP(:))+sum(KEJ(2:JM)))/AREAG
      ediff=(TE-TE0)/(PSF*SHA*mb2kg)        ! C
!$OMP  PARALLEL DO PRIVATE (L)
      do l=1,lm
        T(:,:,L)=T(:,:,L)-ediff/PK(L,:,:)
      end do
!$OMP  END PARALLEL DO

  200 IF (MFILTR.LT.2) RETURN
C****
C**** TEMPERATURE STRATIFICATION FILTER ON T
C****
      AKAP=KAPA-.205d0    ! what is this number?
!$OMP  PARALLEL DO PRIVATE (J,L,X,Y)
      DO L=1,LM
        IF(L.LT.LS1) THEN
          DO J=2,JM-1
            Y(:,J)=(SIG(L)*P(:,J)+PTOP)**AKAP
            X(:,J)=T(:,J,L)*Y(:,J)
          END DO
          CALL SHAP1D (8,X)
          DO J=2,JM-1
            T(:,J,L)=X(:,J)/Y(:,J)
          END DO
        ELSE
          DO J=2,JM-1
            X(:,J)=T(:,J,L)
          END DO
          CALL SHAP1D (8,X)
          DO J=2,JM-1
            T(:,J,L)=X(:,J)
          END DO
        END IF
      END DO
!$OMP  END PARALLEL DO
C
      RETURN
      END SUBROUTINE FILTER


      SUBROUTINE FLTRUV(U,V,UT,VT)
!@sum  FLTRUV Filters 2 gridpoint noise from the velocity fields
!@auth Original development team
!@ver  1.0
      USE CONSTANT, only : sha
      USE MODEL_COM, only : im,jm,lm,byim,mrch,dt,dt_UVfilter,t
      USE GEOM, only : dxyn,dxys,idij,idjj,rapj,imaxj,kmaxj
      USE DYNAMICS, only : pdsig,pk
C**********************************************************************
C**** FILTERING IS DONE IN X-DIRECTION WITH A 8TH ORDER SHAPIRO
C**** FILTER. THE EFFECT OF THE FILTER IS THAT OF DISSIPATION AT
C**** THE SMALLEST SCALES. (Y-dir and 2D filter commented out 'cq')
C**********************************************************************
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,JM,LM), INTENT(INOUT) :: U,V
      REAL*8, DIMENSION(IM,JM,LM), INTENT(IN) :: UT,VT
      REAL*8, DIMENSION(IM,JM,LM) :: DUT,DVT,USAVE,VSAVE,DKE
      REAL*8 X(IM),Y(0:JM+1),F2D(IM,JM),DP,Xby4toN
      REAL*8 :: DT1=0.
cq    LOGICAL*4 QFILY,QFIL2D
      INTEGER I,J,K,L,N,IP1  !@var I,J,L,N  loop variables
      REAL*8 Y4TO8,YJ,YJM1,X1,XI,XIM1
      INTEGER, PARAMETER :: NSHAP=8, ISIGN=(-1.0)**(NSHAP-1)
      REAL*8, PARAMETER :: BY16=1./16., by4toN=1./(4.**NSHAP)
      REAL*8 ediff
C****
cq       QFILY  = .FALSE.
cq       QFIL2D = .FALSE.
C****
      IF (MRCH.eq.2) THEN
        USAVE=U ; VSAVE=V
      END IF
      Xby4toN=(DT/DT_UVfilter)*by4toN
C****
C**** Filtering in east-west direction
C****
!$OMP  PARALLEL DO PRIVATE (I,J,L,N,X,X1,XI,XIM1,F2D)
      DO 350 L=1,LM
C**** Filter U component of momentum
cq    IF(QFIL2D) GOTO 250
      DO 240 J=2,JM
      DO 210 I=1,IM
  210 X(I) = U(I,J,L)
      DO 230 N=1,NSHAP
      X1   = X(1)
      XIM1 = X(IM)
      DO 220 I=1,IM-1
      XI   = X(I)
      X(I) = XIM1-XI-XI+X(I+1)
  220 XIM1 = XI
  230 X(IM)= XIM1-X(IM)-X(IM)+X1
      DO 240 I=1,IM
  240 U(I,J,L) = U(I,J,L) + ISIGN*X(I)*Xby4toN
cq    GOTO 270
cq250 DO 255,J=2,JM-1
cq    F2D(1,J)=.125*(U(IM,J,L)+U(1,J+1,L)+U(2,J,L)+
cq   *                  U(1,J-1,L)-4.*U(1,J,L))+
cq   *        BY16*(U(IM,J-1,L)+U(IM,J+1,L)+U(2,J+1,L)+
cq   *              U(2,J-1,L)-4.*U(1,J,L))
cq    DO 260,I=2,IM-1
cq    F2D(I,J)=.125*(U(I-1,J,L)+U(I,J+1,L)+U(I+1,J,L)+
cq   *                  U(I,J-1,L)-4.*U(I,J,L))+
cq   *        BY16*(U(I-1,J-1,L)+U(I-1,J+1,L)+U(I+1,J+1,L)+
cq   *              U(I+1,J-1,L)-4.*U(I,J,L))
cq260 CONTINUE
cq    F2D(IM,J)=.125*(U(IM-1,J,L)+U(IM,J+1,L)+U(1,J,L)+
cq   *                  U(IM,J-1,L)-4.*U(IM,J,L))+
cq   *        BY16*(U(IM-1,J-1,L)+U(IM-1,J+1,L)+U(1,J+1,L)+
cq   *              U(1,J-1,L)-4.*U(IM,J,L))
cq255 CONTINUE
cq    DO 265,J=2,JM-1
cq    DO 265,I=1,IM
cq      U(I,J,L)=U(I,J,L)+F2D(I,J)
cq265 CONTINUE
cq270 CONTINUE
C**** Filter V component of momentum
      DO 340 J=2,JM
      DO 310 I=1,IM
  310 X(I) = V(I,J,L)
      DO 330 N=1,NSHAP
      X1   = X(1)
      XIM1 = X(IM)
      DO 320 I=1,IM-1
      XI   = X(I)
      X(I) = XIM1-XI-XI+X(I+1)
  320 XIM1 = XI
  330 X(IM)= XIM1-X(IM)-X(IM)+X1
      DO 340 I=1,IM
  340 V(I,J,L) = V(I,J,L) + ISIGN*X(I)*Xby4toN
  350 CONTINUE
!$OMP  END PARALLEL DO
C****
C**** Filtering in north-south direction
C****
cq    IF(QFILY) THEN
cq    Y4TO8 = 1./(4.**NSHAP)
C**** Filter U component of momentum
cq!$OMP  PARALLEL DO PRIVATE (I,J,L,N,Y,YJ,YJM1)
cq    DO 650 L=1,LM
cq    DO 540 I=1,IM
cq    DO 510 J=1,JM
cq510 Y(J) = U(I,J,L)
cq    DO 530 N=1,NSHAP
cq       Y(0)    = Y(2)
cq       Y(JM+1) = Y(JM-1)
cq       YJM1 = Y(0)
cq       DO 520 J=1,JM
cq          YJ   = Y(J)
cq          Y(J) = YJM1-YJ-YJ+Y(J+1)
cq520    YJM1 = YJ
cq530 CONTINUE
cq    DO 540 J=1,JM
cq540 U(I,J,L) = U(I,J,L) + ISIGN*Y(J)*Y4TO8
C**** Make U winds at poles to be uniform
cq    DO 560 J=1,JM,JM-1
cq    YJ = 0.
cq    DO 550 I=1,IM
cq550 YJ = YJ + U(I,J,L)
cq    DO 560 I=1,IM
cq560 U(I,J,L) = YJ*BYIM
C**** Filter V component of momentum
cq    DO 640 I=1,IM
cq    DO 610 J=1,JM-1
cq610 Y(J) = V(I,J,L)
cq    DO 630 N=1,NSHAP
cq    YJM1 = Y(1)
cq    Y(1) = Y(2)-Y(1)
cq    DO 620 J=2,JM-2
cq    YJ   = Y(J)
cq    Y(J) = YJM1-YJ-YJ+Y(J+1)
cq620 YJM1 = YJ
cq630 Y(JM-1)= YJM1-Y(JM-1)
cq    DO 640 J=1,JM-1
cq640 V(I,J,L) = V(I,J,L) + ISIGN*Y(J)*Y4TO8
cq650 CONTINUE
cq!$OMP  END PARALLEL DO
cq    END IF
C****
      IF (MRCH.eq.2) THEN
!$OMP  PARALLEL DO PRIVATE (I,IP1,J,L,DP)
        DO L=1,LM
        DO J=2,JM
        I=IM
        DO IP1=1,IM
          DP=0.5*((PDSIG(L,IP1,J-1)+PDSIG(L,I,J-1))*DXYN(J-1)
     *           +(PDSIG(L,IP1,J  )+PDSIG(L,I,J  ))*DXYS(J  ))
          DKE(I,J,L)=0.5*(U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L)-USAVE(I,J
     *         ,L)*USAVE(I,J,L)-VSAVE(I,J,L)*VSAVE(I,J,L))
          DUT(I,J,L)=(U(I,J,L)-USAVE(I,J,L))*DP
          DVT(I,J,L)=(V(I,J,L)-VSAVE(I,J,L))*DP
          I=IP1
        END DO
        END DO
        END DO
!$OMP  END PARALLEL DO
        CALL DIAGCD(5,UT,VT,DUT,DVT,DT1)

C***** Add in dissipiated KE as heat locally
!$OMP  PARALLEL DO PRIVATE(I,J,L,ediff,K)
        DO L=1,LM
          DO J=1,JM
            DO I=1,IMAXJ(J)
              ediff=0.
              DO K=1,KMAXJ(J)   ! loop over surrounding vel points
                ediff=ediff+DKE(IDIJ(K,I,J),IDJJ(K,J),L)*RAPJ(K,J)
              END DO
              T(I,J,L)=T(I,J,L)-ediff/(SHA*PK(L,I,J))
            END DO
          END DO
        END DO
!$OMP  END PARALLEL DO

      END IF
      RETURN
      END SUBROUTINE FLTRUV


      SUBROUTINE SHAP1D (NORDER,X)
!@sum  SHAP1D Smoothes in zonal direction use n-th order shapiro filter
!@auth Original development team
!@ver  1.0
      USE MODEL_COM, only :im,jm
      IMPLICIT NONE
!@var NORDER order of shapiro filter (must be even)
      INTEGER, INTENT(IN) :: NORDER
      REAL*8, INTENT(INOUT), DIMENSION(IM,JM) :: X
      REAL*8, DIMENSION(IM)::XS
      REAL*8 by4toN,XS1,XSIM1,XSI
      INTEGER I,J,N   !@var I,J,N  loop variables

      by4toN=1./4.**NORDER
      DO J=2,JM-1
        XS(:) = X(:,J)
        DO N=1,NORDER
          XS1=XS(1)
          XSIM1=XS(IM)
          DO I=1,IM-1
            XSI=XS(I)
            XS(I)=XSIM1-XSI-XSI+XS(I+1)
            XSIM1=XSI
          END DO
          XS(IM)=XSIM1-XS(IM)-XS(IM)+XS1
        END DO
        X(:,J)=X(:,J)-XS(:)*by4toN
      END DO
      RETURN
      END SUBROUTINE SHAP1D


      SUBROUTINE CALC_PIJL(lmax,p,pijl)
!@sum  CALC_PIJL Fills in P as 3-D
!@auth Jean Lerner
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm,ls1,psfmpt
C****
      implicit none
      REAL*8, dimension(im,jm) :: p
      REAL*8, dimension(im,jm,lm) :: pijl
      integer :: l,lmax

      do l=1,ls1-1
        pijl(:,:,l) = p(:,:)
      enddo
      do l=ls1,lmax
        pijl(:,:,l) = PSFMPT
      enddo
      return
      end subroutine calc_pijl


      SUBROUTINE CALC_AMPK(LMAX)
!@sum  CALC_AMPK calculate air mass and pressure functions
!@auth Jean Lerner/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : bygrav,kapa
      USE MODEL_COM, only : im,jm,lm,ls1,p,dsig,sig,sige,ptop,psfmpt
      USE DYNAMICS, only : plij,pdsig,pmid,pk,pedn,pek,sqrtp,am,byam
      IMPLICIT NONE

      INTEGER :: I,J,L  !@var I,J,L  loop variables
      INTEGER, INTENT(IN) :: LMAX !@var LMAX max. level for update

C**** Calculate air mass, layer pressures, P**K, and sqrt(P)
C**** Note that only layers LS1 and below vary as a function of surface
C**** pressure. Routine should be called with LMAX=LM at start, and
C**** subsequentaly with LMAX=LS1-1
C**** Note Air mass is calculated in (kg/m^2)

C**** Fill in polar boxes
      P(2:IM,1) = P(1,1)
      P(2:IM,JM)= P(1,JM)

!$OMP  PARALLEL DO PRIVATE (I,J,L)
      DO J=1,JM
        DO I=1,IM
          DO L=1,LS1-1
            PLIJ(L,I,J) = P(I,J)
            PDSIG(L,I,J) = P(I,J)*DSIG(L)
            PMID(L,I,J) = SIG(L)*P(I,J)+PTOP
            PK  (L,I,J) = PMID(L,I,J)**KAPA
            PEDN(L,I,J) = SIGE(L)*P(I,J)+PTOP
            PEK (L,I,J) = PEDN(L,I,J)**KAPA
            AM  (L,I,J) = P(I,J)*DSIG(L)*1d2*BYGRAV
            BYAM(L,I,J) = 1./AM(L,I,J)
          END DO
          DO L=LS1,LMAX
            PLIJ(L,I,J) = PSFMPT
            PDSIG(L,I,J) = PSFMPT*DSIG(L)
            PMID(L,I,J) = SIG(L)*PSFMPT+PTOP
            PK  (L,I,J) = PMID(L,I,J)**KAPA
            PEDN(L,I,J) = SIGE(L)*PSFMPT+PTOP
            PEK (L,I,J) = PEDN(L,I,J)**KAPA
            AM  (L,I,J) = PSFMPT*DSIG(L)*1d2*BYGRAV
            BYAM(L,I,J) = 1./AM(L,I,J)
          END DO
          IF (LMAX.eq.LM) THEN
            PEDN(LM+1,I,J) = SIGE(LM+1)*PSFMPT+PTOP
            PEK (LM+1,I,J) = PEDN(LM+1,I,J)**KAPA
          END IF
          SQRTP(I,J) = SQRT(P(I,J))
        END DO
      END DO
!$OMP  END PARALLEL DO

      RETURN
      END SUBROUTINE CALC_AMPK

      SUBROUTINE CALC_AMP(p,amp)
!@sum  CALC_AMP Calc. AMP: kg air*grav/100, incl. const. pressure strat
!@auth Jean Lerner/Max Kelley
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm,ls1,dsig,psf,ptop
      USE GEOM, only : dxyp
C****
      implicit none
      REAL*8, dimension(im,jm) :: p
      REAL*8, dimension(im,jm,lm) :: amp
      integer :: j,l
C
!$OMP  PARALLEL DO PRIVATE(J,L)
      DO L=1,LM
        IF(L.LT.LS1) THEN
ccc   do l=1,ls1-1
          do j=1,jm
            amp(:,j,l) = p(:,j)*dxyp(j)*dsig(l)
          enddo
ccc   enddo
        ELSE
ccc   do l=ls1,lm
          do j=1,jm
            amp(:,j,l) = (psf-ptop)*dxyp(j)*dsig(l)
          enddo
        END IF
ccc   enddo
      enddo
!$OMP  END PARALLEL DO
C
      return
C****
      end subroutine calc_amp

      SUBROUTINE PGRAD_PBL
!@sum  PGRAD_PBL calculates surface/layer 1 pressure gradients for pbl
!@auth Ye Cheng
!@ver  1.0
C**** As this is written, it must be called after the call to CALC_AMPK
C**** after DYNAM (since it uses pk/pmid). It would be better if it used
C**** SPA and PU directly from the dynamics. (Future work).
      USE CONSTANT, only : rgas
      USE MODEL_COM, only : im,jm,t,p,zatmo,sig,byim
      USE GEOM, only : bydyp,bydxp,cosip,sinip
      USE DYNAMICS, only : phi,dpdy_by_rho,dpdy_by_rho_0,dpdx_by_rho
     *     ,dpdx_by_rho_0,pmid,pk
      IMPLICIT NONE
      REAL*8 by_rho1,dpx1,dpy1,dpx0,dpy0,hemi
      INTEGER I,J,K,IP1,IM1,J1

C**** (Pressure gradient)/density at first layer and surface
C**** to be used in the PBL, at the promary grids

      ! for dPdy/rho at non-pole grids

      DO I=1,IM
        DO J=2,JM-1
          by_rho1=(rgas*t(I,J,1)*pk(1,I,J))/(100.*pmid(1,I,J))
          DPDY_BY_RHO(I,J)=(100.*(P(I,J+1)-P(I,J-1))*SIG(1)*by_rho1
     2         +PHI(I,J+1,1)-PHI(I,J-1,1))*BYDYP(J)*.5d0
          DPDY_BY_RHO_0(I,J)=(100.*(P(I,J+1)-P(I,J-1))*by_rho1
     2         +ZATMO(I,J+1)-ZATMO(I,J-1))*BYDYP(J)*.5d0
        END DO
      END DO

      ! for dPdx/rho at non-pole grids

      DO J=2,JM-1
        IM1=IM-1
        I=IM
        DO IP1=1,IM
          by_rho1=(rgas*t(I,J,1)*pk(1,I,J))/(100.*pmid(1,I,J))
          DPDX_BY_RHO(I,J)=(100.*(P(IP1,J)-P(IM1,J))*SIG(1)*by_rho1
     2         +PHI(IP1,J,1)-PHI(IM1,J,1))*BYDXP(J)*.5d0
          DPDX_BY_RHO_0(I,J)=(100.*(P(IP1,J)-P(IM1,J))*by_rho1
     2         +ZATMO(IP1,J)-ZATMO(IM1,J))*BYDXP(J)*.5d0
          IM1=I
          I=IP1
        END DO
      END DO

      ! at poles

      DO J=1,JM,JM-1
        if (J.eq.1) then
          hemi=-1.; J1=2
        else
          hemi= 1.; J1=JM-1
        endif
        dpx1=0. ; dpy1=0.
        dpx0=0. ; dpy0=0.
        DO K=1,IM
          dpx1=dpx1+(DPDX_BY_RHO(K,J1)*COSIP(K)
     2         -hemi*DPDY_BY_RHO(K,J1)*SINIP(K))
          dpy1=dpy1+(DPDY_BY_RHO(K,J1)*COSIP(K)
     2         +hemi*DPDX_BY_RHO(K,J1)*SINIP(K))
          dpx0=dpx0+(DPDX_BY_RHO_0(K,J1)*COSIP(K)
     2         -hemi*DPDY_BY_RHO_0(K,J1)*SINIP(K))
          dpy0=dpy0+(DPDY_BY_RHO_0(K,J1)*COSIP(K)
     2         +hemi*DPDX_BY_RHO_0(K,J1)*SINIP(K))
        END DO
        DPDX_BY_RHO(1,j)  =dpx1*BYIM
        DPDY_BY_RHO(1,j)  =dpy1*BYIM
        DPDX_BY_RHO_0(1,j)=dpx0*BYIM
        DPDY_BY_RHO_0(1,j)=dpy0*BYIM
      END DO

      END SUBROUTINE PGRAD_PBL

      SUBROUTINE SDRAG(DT1)
!@sum  SDRAG puts a drag on the winds in the top layers of atmosphere
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : grav,rgas,sha
      USE MODEL_COM, only : im,jm,lm,ls1,psfmpt,u,v,sige,bydsig,ptop,t
     *  ,q,x_sdrag,csdragl,lsdrag,lpsdrag,ang_sdrag,itime
     *  ,Wc_Jdrag
      USE GEOM, only : cosv,dxyn,dxys,imaxj,kmaxj,idij,idjj,rapj
      USE DAGCOM, only : ajl,jl_dudtsdrg
      USE DYNAMICS, only : pk,pdsig
      IMPLICIT NONE

!@var DT1 time step (s)
      REAL*8, INTENT(IN) :: DT1
!@var L(P)SDRAG lowest level at which SDRAG_lin is applied (near poles)
C**** SDRAG_const is applied above PTOP (150 mb) and below the SDRAG_lin
C**** regime (but not above P_CSDRAG)
      REAL*8 WL,TL,RHO,CDN,X,BYPIJU,DP,DPL(LM),du
!@var DUT,DVT change in momentum (mb m^3/s)
!@var DKE change in kinetic energy (m^2/s^2)
      REAL*8, DIMENSION(IM,JM,LM) :: DUT,DVT,DKE
      INTEGER I,J,L,IP1,K,Lmax
      logical cd_lin
!@var ang_mom is the sum of angular momentun at layers LS1 to LM
      REAL*8, DIMENSION(IM,JM)    :: ang_mom, sum_airm
!@param wmax imposed limit for stratospheric winds (m/s)
      real*8, parameter :: wmax = 200.d0 , wmaxp = wmax*3.d0/4.d0
      real*8 ediff,wmaxj,xjud

      ang_mom=0. ;  sum_airm=0. ; dke=0. ; dut=0.
C*
      BYPIJU=1./PSFMPT
      DUT=0. ; DVT=0.
      DO L=LS1,LM
      DO J=2,JM
      cd_lin=.false.
      IF( L.ge.LSDRAG .or.
     *   (L.ge.LPSDRAG.and.COSV(J).LE..15) ) cd_lin=.true.
      wmaxj=wmax
      if(COSV(J).LE..15) wmaxj=wmaxp
      I=IM
      DO IP1=1,IM
        TL=T(I,J,L)*PK(L,I,J)
C**** check T to make sure it stayed within physical bounds
        if (TL.lt.100..or.TL.gt.373.) then
          write(99,*) 'SDRAG:',itime,i,j,l,'T,U,V=',TL,U(I,J,L),V(I,J,L)
          call stop_model('Stopped in ATMDYN::SDRAG',11)
        end if
        RHO=(PSFMPT*SIGE(L+1)+PTOP)/(RGAS*TL)
        WL=SQRT(U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L))
        xjud=1.
        if(Wc_JDRAG.gt.0.) xjud=(Wc_JDRAG/(Wc_JDRAG+min(WL,wmaxj)))**2
C**** WL is restricted to Wmax by adjusting X, if necessary;
C**** the following is equivalent to first reducing (U,V), if necessary,
C**** then finding the drag and applying it to the reduced winds
                    CDN=CSDRAGl(l)*xjud
        IF (cd_lin) CDN=(X_SDRAG(1)+X_SDRAG(2)*min(WL,wmaxj))*xjud
        X=DT1*RHO*CDN*min(WL,wmaxj)*GRAV*BYDSIG(L)*BYPIJU
        if (wl.gt.wmaxj) X = 1. - (1.-X)*wmaxj/wl
C**** adjust diags for possible difference between DT1 and DTSRC
        AJL(J,L,JL_DUDTSDRG) = AJL(J,L,JL_DUDTSDRG)-U(I,J,L)*X
        DP=0.5*((PDSIG(L,IP1,J-1)+PDSIG(L,I,J-1))*DXYN(J-1)
     *         +(PDSIG(L,IP1,J  )+PDSIG(L,I,J  ))*DXYS(J  ))
        ang_mom(i,j) = ang_mom(i,j)+U(I,J,L)*X*DP
        DUT(I,J,L)=-X*U(I,J,L)*DP
        DVT(I,J,L)=-X*V(I,J,L)*DP
        DKE(I,J,L)=0.5*(U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L))*(X*X-2.*X)
        U(I,J,L)=U(I,J,L)*(1.-X)
        V(I,J,L)=V(I,J,L)*(1.-X)
        I=IP1
      END DO
      END DO
      END DO
C*
C***  Add the lost angular momentum uniformly back in if ang_sdrag>0
C***  only below 150mb if ang_sdrag=1, into whole column if ang_sdrag>1
C*
      if (ang_sdrag.gt.0) then
        lmax=ls1-1
        if (ang_sdrag.gt.1) lmax=lm
        do j = 2,jm
        I=IM
        do ip1 = 1,im
          do l = 1,lmax
            DPL(L)=0.5*((PDSIG(L,IP1,J-1)+PDSIG(L,I,J-1))*DXYN(J-1)
     *        +(PDSIG(L,IP1,J  )+PDSIG(L,I,J  ))*DXYS(J  ))
            sum_airm(i,j) = sum_airm(i,j)+DPL(L)
          end do
C*
          do l = 1,lmax
            du = ang_mom(i,j)/sum_airm(i,j)
            DUT(I,J,L) = DUT(I,J,L) + du*dpl(l)
            U(I,J,L)=U(I,J,L) + du
            AJL(J,L,JL_DUDTSDRG) = AJL(J,L,JL_DUDTSDRG) + du
            dke(i,j,l) = dke(i,j,l) + .5*du*(2*u(i,j,l)+du)
          end do
          I=IP1
        end do
        end do
      end if

C***** Add in dissipiated KE as heat locally
!$OMP  PARALLEL DO PRIVATE(I,J,L,ediff,K)
      DO L=1,LM
        DO J=1,JM
          DO I=1,IMAXJ(J)
            ediff=0.
            DO K=1,KMAXJ(J)     ! loop over surrounding vel points
              ediff=ediff+DKE(IDIJ(K,I,J),IDJJ(K,J),L)*RAPJ(K,J)
            END DO
            T(I,J,L)=T(I,J,L)-ediff/(SHA*PK(L,I,J))
          END DO
        END DO
      END DO
!$OMP  END PARALLEL DO

C**** conservation diagnostic
C**** (technically we should use U,V from before but this is ok)
      CALL DIAGCD (4,U,V,DUT,DVT,DT1)
      RETURN
      END SUBROUTINE SDRAG

      SUBROUTINE CALC_TROP
!@sum  CALC_TROP (to calculate tropopause height and layer)
!@auth J. Lerner
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm,t
      USE GEOM, only : imaxj
      USE DAGCOM, only : aij, ij_ptrop, ij_ttrop
      USE DYNAMICS, only : pk, pmid, PTROPO, LTROPO
      IMPLICIT NONE
      INTEGER I,J,L,IERR
      REAL*8, DIMENSION(LM) :: TL

C**** Find WMO Definition of Tropopause to Nearest L
!$OMP  PARALLEL DO PRIVATE (I,J,L,TL,IERR)
      do j=1,jm
      do i=1,imaxj(j)
        do l=1,lm
          TL(L)=T(I,J,L)*PK(L,I,J)
        end do
        CALL TROPWMO(TL,PMID(1,I,J),PK(1,I,J),PTROPO(I,J),LTROPO(I,J)
     *       ,IERR)
        IF (IERR.gt.0) print*,"TROPWMO error: ",i,j
        AIJ(I,J,IJ_PTROP)=AIJ(I,J,IJ_PTROP)+PTROPO(I,J)
        AIJ(I,J,IJ_TTROP)=AIJ(I,J,IJ_TTROP)+TL(LTROPO(I,J))
      end do
      end do
!$OMP  END PARALLEL DO
      PTROPO(2:IM,1) = PTROPO(1,1)
      LTROPO(2:IM,1) = LTROPO(1,1)
      PTROPO(2:IM,JM)= PTROPO(1,JM)
      LTROPO(2:IM,JM)= LTROPO(1,JM)

      END SUBROUTINE CALC_TROP

      SUBROUTINE DISSIP
!@sum DISSIP adds in dissipated KE as heat locally
!@auth Gavin Schmidt
      USE CONSTANT, only : sha
      USE MODEL_COM, only : jm,lm,t
      USE GEOM, only : imaxj,kmaxj,idij,idjj,rapj
      USE DYNAMICS, only : dke,pk
      IMPLICIT NONE
      INTEGER I,J,L,K
      REAL*8 ediff

C**** DKE (m^2/s^2) is saved from surf,dry conv,aturb and m.c

!$OMP  PARALLEL DO PRIVATE(I,J,L,ediff,K)
      DO L=1,LM
        DO J=1,JM
          DO I=1,IMAXJ(J)
            ediff=0.
            DO K=1,KMAXJ(J)     ! loop over surrounding vel points
              ediff=ediff+DKE(IDIJ(K,I,J),IDJJ(K,J),L)*RAPJ(K,J)
            END DO
            T(I,J,L)=T(I,J,L)-ediff/(SHA*PK(L,I,J))
          END DO
        END DO
      END DO
!$OMP  END PARALLEL DO

      END SUBROUTINE DISSIP

      subroutine tropwmo(ptm1, papm1, pk, ptropo, ltropp,ierr)
!@sum  tropwmo calculates tropopasue height according to WMO formula
!@auth D. Nodorp/T. Reichler/C. Land
!@+    GISS Modifications by Jean Lerner/Gavin Schmidt
!@ver  1.0
!@alg  WMO Tropopause Definition
!@+
!@+ From A Temperature Lapse Rate Definition of the Tropopause Based on
!@+ Ozone, J. M. Roe and W. H. Jasperson, 1981
!@+
!@+ In the following discussion the lapse rate is defined as -dT/dz.
!@+
!@+ The main features of the WMO tropopause definition are as follows:
!@+ * The first tropopause (i.e., the conventional tropopause) is
!@+   defined as the lowest level at which the lapse rate decreases to 2
!@+   K/km or less, and the average lapse rate from this level to any
!@+   level within the next higher 2 km does not exceed 2 K/km.
!@+ * If above the first tropopause the average lapse rate between any
!@+   level and all higher levels within 1 km exceed 3 K/km, then a
!@+   second tropopause is defined by the same criterion as under the
!@+   statement above. This tropopause may be either within or above the
!@+   1 km layer.
!@+ * A level otherwise satisfying the definition of tropopause, but
!@+   occuring at an altitude below that of the 500 mb level will not be
!@+   designated a tropopause unless it is the only level satisfying the
!@+   definition and the average lapse rate fails to exceed 3 K/km over
!@+   at least 1 km in any higher layer.
!@+ * (GISS failsafe) Some cases occur when the lapse rate never falls
!@+   below 2 K/km. In such cases the failsafe level is that where the
!@+   lapse rate first falls below 3 K/km. If this still doesn't work
!@+   (ever?), the level is set to the pressure level below 30mb.
!@+
      USE MODEL_COM, only : klev=>lm
      USE CONSTANT, only : zkappa=>kapa,zzkap=>bykapa,grav,rgas
      implicit none

      real*8, intent(in), dimension(klev) :: ptm1, papm1, pk
      real*8, intent(out) :: ptropo
      integer, intent(out) :: ltropp,ierr
      real*8, dimension(klev) :: zpmk, zpm, za, zb, ztm, zdtdz
!@param zgwmo min lapse rate (* -1) needed for trop. defn. (-K/km)
!@param zgwmo2 GISS failsafe minimum lapse rate (* -1) (-K/km)
!@param zdeltaz distance to check for lapse rate changes (km)
!@param zfaktor factor for caluclating height from pressure (-rgas/grav)
!@param zplimb min pressure at which to define tropopause (mb)
      real*8, parameter :: zgwmo  = -2d-3, zgwmo2=-3d-3,
     *     zdeltaz = 2000.0, zfaktor = -GRAV/RGAS, zplimb=500.
      real*8 zptph, zp2km, zag, zbg, zasum, zaquer, zptf
      integer iplimb,iplimt, jk, jj, kcount, ltset,l
      logical ldtdz
c****
c****  2. Calculate the height of the tropopause
c****  -----------------------------------------
      ltset = -999
      ierr=0
      iplimb=1
c**** set limits based on pressure
      do jk=2,klev-1
        if (papm1(jk-1).gt.600d0) then
          iplimb=jk
        else
          if (papm1(jk).lt.30d0) exit
        end if
      end do
      iplimt=jk
c****
c****  2.1 compute dt/dz
c****  -----------------
c****       ztm  lineare Interpolation in p**kappa
c****     gamma  dt/dp = a * kappa + papm1(jx,jk)**(kappa-1.)

      do jk=iplimb+1,iplimt       ! -1 ?????
        zpmk(jk)=0.5*(pk(jk-1)+pk(jk))

        zpm(jk)=zpmk(jk)**zzkap ! p mitte

        za(jk)=(ptm1(jk-1)-ptm1(jk))/(pk(jk-1)-pk(jk))
        zb(jk) = ptm1(jk)-(za(jk)*pk(jk))

        ztm(jk)=za(jk)*zpmk(jk)+zb(jk) ! T mitte
        zdtdz(jk)=zfaktor*zkappa*za(jk)*zpmk(jk)/ztm(jk)
      end do
c****
c****  2.2 First test: valid dt/dz ?
c****  -----------------------------
c****
      do 1000 jk=iplimb+1,iplimt-1

c**** GISS failsafe test
        if (zdtdz(jk).gt.zgwmo2.and.ltset.ne.1) then
          ltropp=jk
          ltset =1
        end if
c****
        if (zdtdz(jk).gt.zgwmo .and. ! dt/dz > -2K/km
     &       zpm(jk).le.zplimb) then ! zpm not too low
          ltropp = jk
          ltset = 1
c****
c****  2.3 dtdz is valid > something in German
c****  ----------------------------------------
c****    1.lineare in p^kappa (= Dieters neue Methode)

          zag = (zdtdz(jk)-zdtdz(jk+1))/
     &         (zpmk(jk)-zpmk(jk+1)) ! a-gamma
          zbg = zdtdz(jk+1) - zag*zpmk(jk+1) ! b-gamma
          if(((zgwmo-zbg)/zag).lt.0.) then
            zptf=0.
          else
            zptf=1.
          end if
          zptph = zptf*abs((zgwmo-zbg)/zag)**zzkap
          ldtdz=zdtdz(jk+1).lt.zgwmo
          if(.not.ldtdz) zptph=zpm(jk)
c****
c****  2.4 2nd test: dt/dz above 2km must not be lower than -2K/km
c****  -----------------------------------------------------------
c****
          zp2km = zptph + zdeltaz*zpm(jk)
     &         / ztm(jk)*zfaktor ! p at ptph + 2km
          zasum = 0.0           ! zdtdz above
          kcount = 0            ! number of levels above
c****
c****  2.5 Test until pm < p2km
c****  --------------------------
c****
          do jj=jk,iplimt-1
            if(zpm(jj).gt.zptph) cycle ! doesn't happen
            if(zpm(jj).lt.zp2km) goto 2000 ! ptropo valid
            zasum = zasum+zdtdz(jj)
            kcount = kcount+1
            zaquer = zasum/float(kcount) ! dt/dz mean
            if(zaquer.le.zgwmo) goto 1000 ! dt/dz above < 2K/1000
                                          ! discard it
          end do                ! test next level
          goto 2000
        endif
 1000 continue                  ! next level
 2000 continue

      if (ltset.eq.-999) then
        ltropp=iplimt-1  ! default = last level below 30mb
        print*,"In tropwmo ltropp not set, using default: ltropp ="
     *       ,ltropp
        write(6,'(12(I4,5F10.5,/))') (l,ptm1(l),papm1(l),pk(l),zdtdz(l)
     *       ,zpm(l),l=iplimb+1,iplimt-1)
        ierr=1
      end if
      ptropo = papm1(ltropp)
c****
      return
      end subroutine tropwmo
