C**** QUSEM12 E001M12 SOMTQ QUSB261AM12
C**** QUSBM9=QUSB140M9 with correction to second order moment calc.
C**** Changes for constant pressure above LS1 + REAL*8
C**** QUS   is Russell quadratic upstream scheme for temperature
C**** and water vapor advection, with limits applied to water vapor.
C**** Changes for constant pressure above LS1
C**** FQU,FQV for additional diagnostics
C**** Routines included: AADVT, AADVTX, AADVTY, AADVTZ

      MODULE QUSCOM
!@sum  QUSCOM contains gcm-specific advection parameters/workspace
!@auth Maxwell Kelley
      USE QUSDEF
      IMPLICIT NONE
      SAVE
      INTEGER :: IM,JM,LM
      INTEGER :: XSTRIDE,YSTRIDE,ZSTRIDE
      REAL*8 :: BYIM
C**** AIR MASS FLUXES
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: MFLX
C**** WORKSPACE FOR AADVTX
cc    REAL*8, DIMENSION(:), ALLOCATABLE :: AM,F_I
cc    REAL*8, DIMENSION(:,:), ALLOCATABLE :: FMOM_I
C**** WORKSPACE FOR AADVTY
cc    REAL*8, DIMENSION(:), ALLOCATABLE :: BM,F_J
cc    REAL*8, DIMENSION(:,:), ALLOCATABLE :: FMOM_J
C**** WORKSPACE FOR AADVTZ
cc    REAL*8, DIMENSION(:), ALLOCATABLE :: CM,F_L
cc    REAL*8, DIMENSION(:,:), ALLOCATABLE :: FMOM_L
      END MODULE QUSCOM

      SUBROUTINE init_QUS(IM_GCM,JM_GCM,LM_GCM)
!@sum  init_QUS sets gcm-specific advection parameters/workspace
!@auth Maxwell Kelley
      use QUSCOM
      USE PARAM
      INTEGER, INTENT(IN) :: IM_GCM,JM_GCM,LM_GCM
C**** SET RESOLUTION
      IM = IM_GCM
      JM = JM_GCM
      LM = LM_GCM
      BYIM = 1.D0/REAL(IM,KIND=8)
      XSTRIDE = 1
      YSTRIDE = IM
      ZSTRIDE = IM*JM
C**** ALLOCATE SPACE FOR AIR MASS FLUXES
      ALLOCATE(MFLX(IM,JM,LM))
C**** ALLOCATE WORKSPACE FOR AADVTX
cc    ALLOCATE(AM(IM),F_I(IM),FMOM_I(NMOM,IM))
C**** ALLOCATE WORKSPACE FOR AADVTY
cc    ALLOCATE(BM(JM),F_J(JM),FMOM_J(NMOM,JM))
C**** ALLOCATE WORKSPACE FOR AADVTZ
cc    ALLOCATE(CM(LM),F_L(LM),FMOM_L(NMOM,LM))

      call sync_param("prather_limits",prather_limits)

      RETURN
      END SUBROUTINE init_QUS

      SUBROUTINE AADVT (MA,RM,RMOM,SD,PU,PV,DT,QLIMIT,FQU,FQV)
!@sum  AADVT advection driver
!@auth G. Russell, modified by Maxwell Kelley
c****
c**** AADVT advects tracers using the Quadradic Upstream Scheme.
c****
c**** input:
c****  pu,pv,sd (kg/s) = east-west,north-south,vertical mass fluxes
c****      qlimit = whether moment limitations should be used
C****         DT (s) = time step
c****
c**** input/output:
c****     rm = tracer concentration
c****   rmom = moments of tracer concentration
c****     ma (kg) = fluid mass
c****
      USE QUSDEF
      USE QUSCOM, ONLY : IM,JM,LM, MFLX
      IMPLICIT NONE

      REAL*8, dimension(im,jm,lm) :: rm,ma
      REAL*8, dimension(NMOM,IM,JM,LM) :: rmom

      REAL*8, INTENT(IN) :: DT
      REAL*8, dimension(im,jm,lm), intent(in) :: pu,pv
      REAL*8, dimension(im,jm,lm-1), intent(in) :: sd
      LOGICAL, INTENT(IN) :: QLIMIT

      REAL*8, dimension(im,jm), intent(inout) :: fqu,fqv

      INTEGER :: I,J,L,N
      REAL*8 :: BYMA

C**** Initialise diagnostics
      FQU=0.  ; FQV=0.

C**** Fill in values at the poles
!$OMP  PARALLEL DO PRIVATE(I,L,N)
      DO L=1,LM
         DO I=2,IM
           RM(I,1 ,L) =   RM(1,1 ,L)
           RM(I,JM,L) =   RM(1,JM,L)
           DO N=1,NMOM
             RMOM(N,I,1 ,L) =  RMOM(N,1,1 ,L)
             RMOM(N,I,JM,L) =  RMOM(N,1,JM,L)
         enddo
         enddo
      enddo
!$OMP  END PARALLEL DO
C****
C**** convert from concentration to mass units
C****
!$OMP  PARALLEL DO PRIVATE(I,J,L)
      DO L=1,LM
      DO J=1,JM
      DO I=1,IM
         RM(I,J,L)=RM(I,J,L)*MA(I,J,L)
         RMOM(:,I,J,L)=RMOM(:,I,J,L)*MA(I,J,L)
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO
C****
C**** Advect the tracer using the quadratic upstream scheme
C****
CC    mflx(:,:,:)=pu(:,:,:)*(.5*dt)
!$OMP  PARALLEL DO PRIVATE(L)
       DO L=1,LM
          mflx(:,:,l)=pu(:,:,l)*(.5*dt)
       ENDDO
!$OMP  END PARALLEL DO
      CALL AADVTX (RM,RMOM,MA,MFLX,QLIMIT,FQU)
CC    mflx(:,1:jm-1,:)=pv(:,2:jm,:)*dt
CC    mflx(:,jm,:)=0.
!$OMP  PARALLEL DO PRIVATE(L)
       DO L=1,LM
          mflx(:,1:jm-1,l)=pv(:,2:jm,l)*dt
          mflx(:,jm,l)=0.
       ENDDO
!$OMP  END PARALLEL DO
      CALL AADVTY (RM,RMOM,MA,MFLX,QLIMIT,FQV)
CC    mflx(:,:,1:lm-1)=sd(:,:,1:lm-1)*(-dt)
CC    mflx(:,:,lm)=0.
!$OMP  PARALLEL DO PRIVATE(L)
      DO L=1,LM
         IF(L.NE.LM)  THEN
            MFLX(:,:,L)=SD(:,:,L)*(-DT)
         ELSE
            MFLX(:,:,L)=0.
         END IF
      ENDDO
!$OMP  END PARALLEL DO
      CALL AADVTZ (RM,RMOM,MA,MFLX,QLIMIT)
CC    mflx(:,:,:)=pu(:,:,:)*(.5*dt)
!$OMP  PARALLEL DO PRIVATE(L)
       DO L=1,LM
          mflx(:,:,l)=pu(:,:,l)*(.5*dt)
       ENDDO
!$OMP  END PARALLEL DO
      CALL AADVTX (RM,RMOM,MA,MFLX,QLIMIT,FQU)
C****
C**** convert from mass to concentration units
C****
!$OMP  PARALLEL DO PRIVATE(I,J,L,BYMA)
      DO L=1,LM
      DO J=1,JM
      DO I=1,IM
         BYMA = 1.D0/MA(I,J,L)
         RM(I,J,L)=RM(I,J,L)*BYMA
         RMOM(:,I,J,L)=RMOM(:,I,J,L)*BYMA
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO
      RETURN
      END

      subroutine aadvtx(rm,rmom,mass,mu,qlimit,fqu)
!@sum  AADVTX advection driver for x-direction
!@auth Maxwell Kelley
c****
c**** aadvtx advects tracers in the west to east direction using the
c**** quadratic upstream scheme.  if qlimit is true, the moments are
c**** limited to prevent the mean tracer from becoming negative.
c****
c**** input:
c****     mu (kg) = west-east mass flux, positive eastward
c****      qlimit = whether moment limitations should be used
c****
c**** input/output:
c****     rm (kg) = tracer mass
c****   rmom (kg) = moments of tracer mass
c****   mass (kg) = fluid mass
c****
      use QUSDEF
ccc   use QUSCOM, only : im,jm,lm, xstride,am,f_i,fmom_i
      use QUSCOM, only : im,jm,lm, xstride
      implicit none
      REAL*8, dimension(im,jm,lm) :: rm,mass,mu,hfqu
      REAL*8, dimension(NMOM,IM,JM,LM) :: rmom
      logical ::  qlimit
      REAL*8, INTENT(OUT), DIMENSION(IM,JM) :: FQU
      REAL*8  AM(IM), F_I(IM), FMOM_I(NMOM,IM)
      integer :: i,j,l,ierr,nerr,ICKERR
c**** loop over layers and latitudes
      ICKERR=0
!$OMP  PARALLEL DO PRIVATE(J,L,AM,F_I,FMOM_I,IERR,NERR)
!$OMP* SHARED(IM,QLIMIT,XSTRIDE)
!$OMP* REDUCTION(+:ICKERR)
      do l=1,lm
      do j=2,jm-1
      am(:) = mu(:,j,l)
c****
c**** call 1-d advection routine
c****
      call adv1d(rm(1,j,l),rmom(1,1,j,l), f_i,fmom_i, mass(1,j,l),
     &        am, im, qlimit,xstride,xdir,ierr,nerr)
      if (ierr.gt.0) then
        write(6,*) "Error in aadvtx: i,j,l=",nerr,j,l
        if (ierr.eq.2) then
          write(0,*) "Error in qlimit: abs(a) > 1"
CCC       call stop_model('Error in qlimit: abs(a) > 1',11)
          ICKERR=ICKERR+1
        end if
      end if
c****
c**** store tracer flux in fqu array
c****
CCC   fqu(:,j)  = fqu(:,j) + f_i(:)
      hfqu(:,j,l)  = f_i(:)
      enddo ! j
      enddo ! l
!$OMP  END PARALLEL DO
c
c     now sum into fqu
c
      do l=1,lm
      do j=2,jm-1
         fqu(:,j)  = fqu(:,j) + hfqu(:,j,l)
      enddo ! j
      enddo ! l
C
      IF(ICKERR.GT.0)  CALL stop_model('Stopped in aadvtx',11)
C
      return
c****
      end subroutine aadvtx

      subroutine aadvty(rm,rmom,mass,mv,qlimit,fqv)
!@sum  AADVTY advection driver for y-direction
!@auth Maxwell Kelley
c****
c**** aadvty advects tracers in the south to north direction using the
c**** quadratic upstream scheme.  if qlimit is true, the moments are
c**** limited to prevent the mean tracer from becoming negative.
c****
c**** input:
c****     mv (kg) = north-south mass flux, positive northward
c****      qlimit = whether moment limitations should be used
c****
c**** input/output:
c****     rm (kg) = tracer mass
c****   rmom (kg) = moments of tracer mass
c****   mass (kg) = fluid mass
c****
      use QUSDEF
ccc   use QUSCOM, only : im,jm,lm, ystride,bm,f_j,fmom_j, byim
      use QUSCOM, only : im,jm,lm, ystride,               byim
      implicit none
      REAL*8, dimension(im,jm,lm) :: rm,mass,mv
      REAL*8, dimension(NMOM,IM,JM,LM) :: rmom
      logical ::  qlimit
      REAL*8, intent(out), dimension(im,jm) :: fqv
      REAL*8  HFQV(IM,JM,LM),BM(JM),F_J(JM),FMOM_J(NMOM,JM)
      integer :: i,j,l,ierr,nerr,ICKERR
      REAL*8 ::
     &     m_sp,m_np,rm_sp,rm_np,rzm_sp,rzm_np,rzzm_sp,rzzm_np
c**** loop over layers
      ICKERR=0
!$OMP  PARALLEL DO PRIVATE(I,L,M_SP,M_NP,RM_SP,RM_NP,RZM_SP,RZZM_SP,
!$OMP*             F_J,FMOM_J,RZM_NP,RZZM_NP,BM,IERR,NERR)
!$OMP* SHARED(JM,QLIMIT,YSTRIDE)
!$OMP* REDUCTION(+:ICKERR)
      do l=1,lm
c**** scale polar boxes to their full extent
      mass(:,1:jm:jm-1,l)=mass(:,1:jm:jm-1,l)*im
      m_sp = mass(1,1 ,l)
      m_np = mass(1,jm,l)
      rm(:,1:jm:jm-1,l)=rm(:,1:jm:jm-1,l)*im
      rm_sp = rm(1,1 ,l)
      rm_np = rm(1,jm,l)
      do i=1,im
         rmom(:,i,1 ,l)=rmom(:,i,1 ,l)*im
         rmom(:,i,jm,l)=rmom(:,i,jm,l)*im
      enddo
      rzm_sp  = rmom(mz ,1,1 ,l)
      rzzm_sp = rmom(mzz,1,1 ,l)
      rzm_np  = rmom(mz ,1,jm,l)
      rzzm_np = rmom(mzz,1,jm,l)
c**** loop over longitudes
      do i=1,im
c****
c**** load 1-dimensional arrays
c****
      bm   (:) = mv(i,:,l) !/nstep
      bm(jm)= 0.
      rmom(ihmoms,i,1 ,l) = 0.! horizontal moments are zero at pole
      rmom(ihmoms,i,jm,l) = 0.
c****
c**** call 1-d advection routine
c****
      call adv1d(rm(i,1,l),rmom(1,i,1,l), f_j,fmom_j, mass(i,1,l),
     &     bm, jm,qlimit,ystride,ydir,ierr,nerr)
      if (ierr.gt.0) then
        write(6,*) "Error in aadvty: i,j,l=",i,nerr,l
        if (ierr.eq.2) then
          write(0,*) "Error in qlimit: abs(b) > 1"
ccc       call stop_model('Error in qlimit: abs(b) > 1',11)
          ICKERR=ICKERR+1
        endif
      end if
c**** store tracer flux in fqv array
ccc   fqv(i,:) = fqv(i,:) + f_j(:)
ccc   fqv(i,jm) = 0.   ! play it safe
      hfqv(i,:,l) = f_j(:)
      rmom(ihmoms,i,1 ,l) = 0.! horizontal moments are zero at pole
      rmom(ihmoms,i,jm,l) = 0.
      enddo ! end loop over longitudes
c**** average and unscale polar boxes
      mass(:,1 ,l) = (m_sp + sum(mass(:,1 ,l)-m_sp))*byim
      mass(:,jm,l) = (m_np + sum(mass(:,jm,l)-m_np))*byim
      rm(:,1 ,l) = (rm_sp + sum(rm(:,1 ,l)-rm_sp))*byim
      rm(:,jm,l) = (rm_np + sum(rm(:,jm,l)-rm_np))*byim
      rmom(mz ,:,1 ,l) = (rzm_sp  + sum(rmom(mz ,:,1 ,l)-rzm_sp ))*byim
      rmom(mzz,:,1 ,l) = (rzzm_sp + sum(rmom(mzz,:,1 ,l)-rzzm_sp))*byim
      rmom(mz ,:,jm,l) = (rzm_np  + sum(rmom(mz ,:,jm,l)-rzm_np ))*byim
      rmom(mzz,:,jm,l) = (rzzm_np + sum(rmom(mzz,:,jm,l)-rzzm_np))*byim
      enddo ! end loop over levels
!$OMP  END PARALLEL DO
c
c     sum into fqv
c
      do l=1,lm
      do i=1,im
         fqv(i,:) = fqv(i,:) + hfqv(i,:,l)
         fqv(i,jm) = 0.
      enddo
      enddo
C
      IF(ICKERR.GT.0)  CALL stop_model('Stopped in aadvty',11)
C
      return
c****
      end subroutine aadvty


      subroutine aadvtz(rm,rmom,mass,mw,qlimit)
!@sum  AADVTZ advection driver for z-direction
!@auth Maxwell Kelley
c****
c**** aadvtz advects tracers in the upward vertical direction using the
c**** quadratic upstream scheme.  if qlimit is true, the moments are
c**** limited to prevent the mean tracer from becoming negative.
c****
c**** input:
c****     mw (kg) = vertical mass flux, positive upward
c****      qlimit = whether moment limitations should be used
c****
c**** input/output:
c****     rm (kg) = tracer mass
c****   rmom (kg) = moments of tracer mass
c****   mass (kg) = fluid mass
c****
      use QUSDEF
ccc   use QUSCOM, only : im,jm,lm, zstride,cm,f_l,fmom_l
      use QUSCOM, only : im,jm,lm, zstride
      implicit none
      REAL*8, dimension(im,jm,lm) :: rm,mass,mw
      REAL*8, dimension(NMOM,IM,JM,LM) :: rmom
      logical ::  qlimit
      REAL*8  CM(LM),F_L(LM),FMOM_L(NMOM,LM)
      integer :: i,j,l,ierr,nerr,ICKERR
c**** loop over latitudes and longitudes
      ICKERR=0
!$OMP  PARALLEL DO PRIVATE(I,J,CM,F_L,FMOM_L,IERR,NERR)
!$OMP* SHARED(LM,QLIMIT,ZSTRIDE)
!$OMP* REDUCTION(+:ICKERR)
      do j=1,jm
      do i=1,im
      cm(:) = mw(i,j,:)
      cm(lm)= 0.
c****
c**** call 1-d advection routine
c****
      call adv1d(rm(i,j,1),rmom(1,i,j,1),f_l,fmom_l,mass(i,j,1),
     &        cm,lm,qlimit,zstride,zdir,ierr,nerr)
      if (ierr.gt.0) then
        write(6,*) "Error in aadvtz: i,j,l=",i,j,nerr
        if (ierr.eq.2) then
          write(0,*) "Error in qlimit: abs(c) > 1"
ccc       call stop_model('Error in qlimit: abs(c) > 1',11)
          ICKERR=ICKERR+1
        endif
      end if
      enddo ! i
      enddo ! j
!$OMP  END PARALLEL DO
C
      IF(ICKERR.GT.0) call stop_model('Stopped in aadvtz',11)
      return
c****
      end subroutine aadvtz
