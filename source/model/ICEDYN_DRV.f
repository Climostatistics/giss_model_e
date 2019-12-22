#include "rundeck_opts.h"

      MODULE ICEDYN_COM
!@sum  ICEDYN_COM holds global variables for dynamic sea ice
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,jm
      USE ICEDYN, only : imic,jmic
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm
#endif
      IMPLICIT NONE
      SAVE

C**** variables used in advection (on ICE grid)
!@var RSIX,RSIY first order moments for seaice concentration
!@var USI,VSI east-west, and north-south sea ice velocities (m/s)
      REAL*8, DIMENSION(IMIC,JMIC) :: RSIX,RSIY,USI,VSI

!@var USIDT,VSIDT sea ice fluxes, saved for advection (m)
      REAL*8, DIMENSION(IMIC,JMIC) :: USIDT,VSIDT

C**** Needed for ADVSI (on ATM grid)
!@var RSISAVE saved value of sea ice concentration before DYNSI
      REAL*8, DIMENSION(IM,JM) :: RSISAVE

C**** Ice advection diagnostics
      INTEGER, PARAMETER :: KICIJ=12
!@var IJ_xxx Names for ICIJ diagnostics
      INTEGER IJ_USI,IJ_VSI,IJ_DMUI,IJ_DMVI,IJ_PICE,IJ_MUSI,IJ_MVSI
     *     ,IJ_HUSI,IJ_HVSI,IJ_SUSI,IJ_SVSI,IJ_RSI
!@var ICIJ lat-lon ice dynamic diagnostics (on atm grid)
      REAL*8, DIMENSION(IMIC,JMIC,KICIJ)  :: ICIJ
!@var lname_icij Long names for ICIJ diagnostics
      CHARACTER*50, DIMENSION(KICIJ) :: LNAME_ICIJ
!@var sname_icij Short names for ICIJ diagnostics
      CHARACTER*30, DIMENSION(KICIJ) :: SNAME_ICIJ
!@var units_icij Units for ICIJ diagnostics
      CHARACTER*50, DIMENSION(KICIJ) :: UNITS_ICIJ
!@var ia_icij IDACC numbers for ICIJ diagnostics
      INTEGER, DIMENSION(KICIJ) :: IA_ICIJ
!@var scale_icij scales for ICIJ diagnostics
      REAL*8, DIMENSION(KICIJ) :: SCALE_ICIJ
!@var ijgrid_icij Grid descriptor for ICIJ diagnostics
       INTEGER, DIMENSION(KICIJ) :: IJGRID_ICIJ
#ifdef TRACERS_WATER
!@var KTICIJ number of lat/lon ice dynamic tracer diagnostics
      INTEGER, PARAMETER :: KTICIJ=2
!@var TICIJ  lat/lon ice dynamic tracer diagnostics
      REAL*8, DIMENSION(IMIC,JMIC,KTICIJ,NTM)  :: TICIJ
!@var ticij_xxx indices for TICIJ diags
      INTEGER :: tICIJ_tusi,tICIJ_tvsi
!@var lname_ticij Long names for TICIJ diagnostics
      CHARACTER*50, DIMENSION(KTICIJ) :: LNAME_TICIJ
!@var sname_ticij Short names for TICIJ diagnostics
      CHARACTER*30, DIMENSION(KTICIJ) :: SNAME_TICIJ
!@var units_ticij Units for TICIJ diagnostics
      CHARACTER*50, DIMENSION(KTICIJ) :: UNITS_TICIJ
!@var ia_ticij IDACC numbers for TICIJ diagnostics
      INTEGER, DIMENSION(KTICIJ) :: IA_TICIJ
!@var scale_ticij scales for TICIJ diagnostics
      REAL*8, DIMENSION(KTICIJ) :: SCALE_TICIJ
!@var ijgrid_ticij Grid descriptor for TICIJ diagnostics
      INTEGER, DIMENSION(KTICIJ) :: IJGRID_TICIJ
#endif

      END MODULE ICEDYN_COM

      SUBROUTINE io_icedyn(kunit,iaction,ioerr)
!@sum  io_icedyn reads and writes dynamic ice arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,irsfic,irsficno,irsficnt
     *     ,irerun,lhead 
      USE ICEDYN_COM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "ICEDYN01"

      write(MODULE_HEADER(lhead+1:80),'(a7,i3,a1,i3,a)')
     *     'R8 dim(',imic,',',jmic,'):RSIX,RSIY,USI,VSI'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,RSIX,RSIY,USI,VSI
      CASE (IOREAD:)            ! input from restart file
        SELECT CASE (IACTION)
        CASE (IRSFICNO)           ! initial conditions (no ocean)
        CASE (ioread,irerun,irsfic,irsficnt)    ! restarts
          READ (kunit,err=10) HEADER,RSIX,RSIY,USI,VSI
          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
            GO TO 10
          END IF
        END SELECT
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
C****
      END SUBROUTINE io_icedyn

      SUBROUTINE io_icdiag(kunit,it,iaction,ioerr)
!@sum  io_icdiag reads and writes ice dynamic diagnostic arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,iowrite_mon,iowrite_single
     *     ,irsfic,irsficnt,irerun,ioread_single,lhead
      USE ICEDYN_COM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "ICDIAG01"
!@var it input/ouput value of hour
      INTEGER, INTENT(INOUT) :: it
!@var ICIJ4 dummy arrays for reading diag. files
      REAL*4, DIMENSION(IMIC,JMIC,KICIJ)  :: ICIJ4
#ifdef TRACERS_WATER
      REAL*4, DIMENSION(IMIC,JMIC,KTICIJ,NTM)  :: TICIJ4
!@var TR_HEADER Character string label for individual tracer records
      CHARACTER*80 :: TR_HEADER, TR_MODULE_HEADER = "TRICDIAG01"
      write(TR_MODULE_HEADER(lhead+1:80),'(a9,i3,a1,i3,a1,i2,a1,i2,a4)')
     *     'R8 Ticij(',imic,',',jmic,',',kticij,',',ntm,'),it'
#endif
      write(MODULE_HEADER(lhead+1:80),'(a8,i3,a1,i3,a1,i2,a4)')
     *     'R8 ICij(',imic,',',jmic,',',kicij,'),it'

      SELECT CASE (IACTION)
      CASE (IOWRITE,IOWRITE_MON)  ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,ICIJ,it
#ifdef TRACERS_WATER
        WRITE (kunit,err=10) TR_MODULE_HEADER,TICIJ,it
#endif
      CASE (IOWRITE_SINGLE)    ! output to acc file
        MODULE_HEADER(LHEAD+1:LHEAD+2) = 'R4'
        WRITE (kunit,err=10) MODULE_HEADER,REAL(ICIJ,KIND=4),it
#ifdef TRACERS_WATER
        TR_MODULE_HEADER(LHEAD+1:LHEAD+2) = 'R4'
        WRITE (kunit,err=10) TR_MODULE_HEADER,REAL(TICIJ,KIND=4),it
#endif
      CASE (IOREAD:)            ! input from restart file
        SELECT CASE (IACTION)
        CASE (ioread_single)    ! accumulate diagnostic files
          READ (kunit,err=10) HEADER,ICIJ4,it
C**** accumulate diagnostics
          ICIJ=ICIJ+ICIJ4
          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",HEADER
     *           ,MODULE_HEADER
            GO TO 10
          END IF
#ifdef TRACERS_WATER
          READ (kunit,err=10) TR_HEADER,TICIJ4,it
C**** accumulate diagnostics
          TICIJ=TICIJ+TICIJ4
          IF (TR_HEADER(1:LHEAD).NE.TR_MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",TR_HEADER
     *           ,TR_MODULE_HEADER
            GO TO 10
          END IF
#endif
        CASE (ioread,irerun)    ! restarts
          READ (kunit,err=10) HEADER,ICIJ,it
          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",HEADER
     *           ,MODULE_HEADER
            GO TO 10
          END IF
#ifdef TRACERS_WATER
          READ (kunit,err=10) TR_HEADER,TICIJ,it
          IF (TR_HEADER(1:LHEAD).NE.TR_MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",TR_HEADER
     *           ,TR_MODULE_HEADER
            GO TO 10
          END IF
#endif
        CASE (IRSFIC)  ! initial conditions 
          READ (kunit)
#ifdef TRACERS_WATER
          READ (kunit)
#endif
        CASE (IRSFICNT)  ! initial conditions (with no tracers)
          READ (kunit)
        END SELECT
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
C****
      END SUBROUTINE io_icdiag

      SUBROUTINE reset_icdiag
!@sum reset_icdiag resets ice dynamic diagnostic arrays
!@auth Gavin Schmidt
      USE ICEDYN_COM
      IMPLICIT NONE

      ICIJ=0.
#ifdef TRACERS_WATER
      TICIJ=0.
#endif
      RETURN
      END SUBROUTINE reset_icdiag

      SUBROUTINE DYNSI
!@sum  DYNSI calculate ice velocites
!@+    Note that the ice velocities are calculated on the ice grid
!@auth Jiping Liu/Gavin Schmidt (based on code from J. Zhang)
!@ver  1.0
      USE CONSTANT, only : rhoi,grav,omega,rhows
      USE MODEL_COM, only : im,jm,p,ptop,dts=>dtsrc,focean
      USE GEOM, only : dxyn,dxys,dxyv,dxyp,bydxyp,dxp,dyv,imaxj
      USE ICEDYN, only : imic,jmic,nx1,ny1,press,heffm,uvm,dwatn,cor
     *     ,sinwat,coswat,bydts,sinen,uice,vice,heff,area,gairx,gairy
     *     ,gwatx,gwaty,pgfub,pgfvb,amass,uicec,vicec,uib,vib,dmu,dmv
      USE ICEDYN_COM, only : usi,vsi,usidt,vsidt,rsisave,icij,ij_usi
     *     ,ij_vsi,ij_dmui,ij_dmvi,ij_pice,ij_rsi
      USE FLUXES, only : dmua,dmva,dmui,dmvi,UI2rho,ogeoza,uosurf,vosurf
     *     ,apress,uisurf,visurf
      USE SEAICE, only : ace1i
      USE SEAICE_COM, only : rsi,msi,snowi
      IMPLICIT NONE
      SAVE
C**** intermediate calculation for pressure gradient terms
      REAL*8, DIMENSION(IM,JM) :: PGFU,PGFV
C****
      REAL*8, PARAMETER :: BYRHOI=1D0/RHOI
      REAL*8 :: hemi
      INTEGER I,J,n,k,ip1,im1,l
      REAL*8 USINP,DMUINP,duA,dvA

C**** Start main loop
C**** Replicate polar boxes
      RSI(2:IM,JM)=RSI(1,JM)
      MSI(2:IM,JM)=MSI(1,JM)
      DMUA(2:IM,JM,2) = DMUA(1,JM,2)
      DMVA(2:IM,JM,2) = DMVA(1,JM,2)

C**** save current value of sea ice concentration for ADVSI
C**** RSISAVE is on atmospheric grid
      RSISAVE(:,:)=RSI(:,:)

C**** Pressure anomaly at surface APRESS: calculated by sea ice routines
C**** APRESS is on atmospheric grid. We are now no longer using this as
C**** a forcing in the sea ice dynamics because it is partially
C**** included already in the internal ice pressure gradient. The
C**** atmospheric pressure gradient term does not in general produce a
C**** horizontal force in a solid (such as ice).

C**** calculate sea surface tilt on atmospheric C grid
C**** (using OGEOZA on atmospheric grid plus displacement of free
C**** surface due to presence of ice). This is ignored in favour of
C**** geostrophy if osurf_tilt=0.
C**** PGF is an accelaration
      PGFU(1:IM,JM)=0
      DO J=2,JM-1
        I=IM
        DO IP1=1,IM
          IF(FOCEAN(I,J).gt.0 .and. FOCEAN(IP1,J).gt.0. .and.
     *         RSI(I,J)+RSI(IP1,J).gt.0.) THEN
c            PGFU(I,J)=-((APRESS(IP1,J)-APRESS(I,J))*BYRHOI
c     *           +OGEOZA(IP1,J)-OGEOZA(I,J))/DXP(J)
            PGFU(I,J)=-(OGEOZA(IP1,J)-OGEOZA(I,J)+
     *           (RSI(IP1,J)*(MSI(IP1,J)+SNOWI(IP1,J)+ACE1I)
     *           -RSI(I,J)*(MSI(I,J)+SNOWI(I,J)+ACE1I))*GRAV/RHOWS )
     *           /DXP(J)
          ELSE
            PGFU(I,J)=0.
          END IF
          I=IP1
        END DO
      END DO
      DO J=1,JM-1
        DO I=1,IM
          IF(FOCEAN(I,J+1).gt.0 .and. FOCEAN(I,J).gt.0. .and.
     *         RSI(I,J)+RSI(I,J+1).gt.0.) THEN
c            PGFV(I,J)=-((APRESS(I,J+1)-APRESS(I,J))*BYRHOI
c     *           +OGEOZA(I,J+1)-OGEOZA(I,J))/DYV(J+1)
            PGFV(I,J)=-(OGEOZA(I,J+1)-OGEOZA(I,J)+
     *           (RSI(I,J+1)*(MSI(I,J+1)+SNOWI(I,J+1)+ACE1I) -
     *           RSI(I,J)*(MSI(I,J)+SNOWI(I,J)+ACE1I))*GRAV/RHOWS  )
     *           /DYV(J+1)
          ELSE
            PGFV(I,J)=0.
          END IF
        END DO
      END DO

C**** DMUA is defined over the whole box (not just over ptype)
C**** Convert to stress over ice fraction only (on atmospheric grid)
      DO J=1,JM
        DO I=1,IM
          IF (FOCEAN(I,J)*RSI(I,J).gt.0) THEN
            DMUA(I,J,2) = DMUA(I,J,2)/(FOCEAN(I,J)*RSI(I,J))
            DMVA(I,J,2) = DMVA(I,J,2)/(FOCEAN(I,J)*RSI(I,J))
          ELSE
            DMUA(I,J,2) = 0.
            DMVA(I,J,2) = 0.
          END IF
        END DO
      END DO

C**** Set up ice grid variables
C**** HEFF,AREA on primary (tracer) grid for ice
      do j=1,NY1
        do i=2,NX1-1
          HEFF(I,J)=RSI(I-1,J)*(ACE1I+MSI(I-1,J))*BYRHOI
          AREA(I,J)=RSI(I-1,J)
c if heffm(i,j).eq.1 .and. focean(i-1,j).eq.0 interpolate.....
c          if (heffm(i,j).eq.1 .and. focean(i-1,j).eq.0) then
c     rsifix=(rsi(i-2,j)*focean(i-2,j)+rsi(i,j)*focean(i,j)+rsi(i-1,j+1)
c    *  *focean(i-1,j+1)+rsi(i-1,j-1)*focean(i-1,j-1))/(focean(i-2,j)
c    *  +focean(i,j)+focean(i-1,j+1)+focean(i-1,j-1))
c     msifix=(msi(i-2,j)*focean(i-2,j)+msi(i,j)*focean(i,j)+msi(i-1,j+1)
c    *  *focean(i-1,j+1)+msi(i-1,j-1)*focean(i-1,j-1))/(focean(i-2,j)
c    *  +focean(i,j)+focean(i-1,j+1)+focean(i-1,j-1))
c      HEFF(I,J)=rsifix*(ACE1I+msifix)*BYRHOI
c      AREA(I,J)=rsifix
c          end if
          HEFF(I,J)=HEFF(I,J)*HEFFM(I,J)
        enddo
      enddo
C**** fill in overlap regions
      DO J=1,NY1
        HEFF(1,J)=HEFF(NX1-1,J)
        AREA(1,J)=AREA(NX1-1,J)
        HEFF(NX1,J)=HEFF(2,J)
        AREA(NX1,J)=AREA(2,J)
      END DO

C****
C**** Set up mass per unit area and coriolis term (on ice grid A)
C****
      DO J=1,NY1-1
      DO I=1,NX1-1
        AMASS(I,J)=RHOI*0.25*(HEFF(I,J)
     *       +HEFF(I+1,J)+HEFF(I,J+1)+HEFF(I+1,J+1))
        COR(I,J)=AMASS(I,J)*2.0*OMEGA*SINEN(I,J)
      END DO
      END DO

c**** interpolate air, current and ice velocity from C grid to B grid
C**** This should be more generally from ocean grid to ice grid
C**** NOTE: UOSURF, VOSURF are expected to be on the C-grid
      do j=1,jm-1
        im1=im
        do i=1,im
          UIB  (i,j)=0.5*(USI (im1,j)  +USI (im1,j+1))   ! iceC--> iceB
          GWATX(i,j)=0.5*(UOSURF(im1,j)+UOSURF(im1,j+1)) ! ocnC--> iceB
          PGFUB(i,j)=0.5*(PGFU(im1,j)  +PGFU(im1,j+1))   ! atmC--> iceB
          VIB  (i,j)=0.5*(VSI (im1,j)  +VSI (i,j))
          GWATY(i,j)=0.5*(VOSURF(im1,j)+VOSURF(i,j))
          PGFVB(i,j)=0.5*(PGFV(im1,j)  +PGFV(i,j))
          im1=i
        enddo
      enddo
c**** set north pole
      do i=1,im
        UIB  (i,jm)=USI(1,jm)
        GWATX(i,jm)=UOSURF(1,jm)
        PGFUB(i,jm)=PGFU(1,jm)
        VIB  (i,jm)=0.
        GWATY(i,jm)=0.
        PGFVB(i,jm)=0.
      enddo
      DO J=1,NY1
        UIB(nx1-1,J)=UIB(1,J)
        VIB(nx1-1,J)=VIB(1,J)
        GWATX(nx1-1,J)=GWATX(1,J)
        GWATY(nx1-1,J)=GWATY(1,J)
        PGFUB(nx1-1,J)=PGFUB(1,J)
        PGFVB(nx1-1,J)=PGFVB(1,J)
        UIB(nx1,J)=UIB(2,J)
        VIB(nx1,J)=VIB(2,J)
        GWATX(nx1,J)=GWATX(2,J)
        GWATY(nx1,J)=GWATY(2,J)
        PGFUB(nx1,J)=PGFUB(2,J)
        PGFVB(nx1,J)=PGFVB(2,J)
      ENDDO

c**** interpolate air stress from A grid in atmos, to B grid in ice
C**** change of unit from change of momentum, to flux
      do j=1,jm-1
        im1=im
        do i=1,im
          GAIRX(i,j)=0.25*(dmua(i,j,2)+dmua(im1,j,2)+dmua(im1,j+1,2)
     #         +dmua(i,j+1,2))*bydts
          GAIRY(i,j)=0.25*(dmva(i,j,2)+dmva(im1,j,2)+dmva(im1,j+1,2)
     #         +dmva(i,j+1,2))*bydts
          im1=i
        enddo
      enddo
      GAIRX(1:im,jm)=dmua(1,jm,2)*bydts
      GAIRY(1:im,jm)=dmva(1,jm,2)*bydts

      do j=1,ny1
       GAIRX(nx1-1,j)=GAIRX(1,j)
       GAIRY(nx1-1,j)=GAIRY(1,j)
       GAIRX(nx1,j)=GAIRX(2,j)
       GAIRY(nx1,j)=GAIRY(2,j)
      enddo

c**** read in sea ice velocity
      DO J=1,NY1
      DO I=1,NX1
       UICE(I,J,1)=UIB(I,J)
       VICE(I,J,1)=VIB(I,J)
       UICE(I,J,2)=0.
       VICE(I,J,2)=0.
       UICE(I,J,3)=0.
       VICE(I,J,3)=0.
      END DO
      END DO

C**** do the looping over pseudo-timesteps
      CALL VPICEDYN

C**** Calculate stress on ice velocity grid (B grid)
      DO J=1,NY1
        hemi=-1.
        if (J.gt.NY1/2) hemi=1.
        DO I=1,NX1
          DMU(i,j)=DTS*dwatn(i,j)*(COSWAT*(UICE(i,j,1)-GWATX(i,j))-
     *         HEMI*SINWAT*(VICE(i,j,1)-GWATY(i,j)))
          DMV(i,j)=DTS*dwatn(i,j)*(HEMI*SINWAT*(UICE(i,j,1)-GWATX(i,j))
     *         +COSWAT*(VICE(i,j,1)-GWATY(i,j)))
        END DO
      END DO

C**** interpolate ice velocity and stress from B grid to C grid in atm
      do j=2,jm
        i=im
        do ip1=1,im
          usi(i,j)=0.5*(uice(i+1,j-1,1)+uice(i+1,j,1))
          IF (abs(USI(I,J)).lt.1d-10) USI(I,J)=0
          DMUI(I,J) = 0.5*(dmu(i+1,j-1)+dmu(i+1,j))
C**** Rescale DMUI to be net momentum into ocean
          DMUI(I,J) = 0.5*DMUI(I,J)*(FOCEAN(I,J)*RSI(I,J)+FOCEAN(IP1,J)
     *         *RSI(IP1,J))
          i=ip1
        enddo
      enddo
      do j=1,jm-1
        do i=1,im
          vsi(i,j)=0.5*(vice(i,j,1)+vice(i+1,j,1))
          IF (abs(VSI(I,J)).lt.1d-10) VSI(I,J)=0
          DMVI(I,J) = 0.5*(dmv(i,j)+dmv(i+1,j))
C**** Rescale DMVI to be net momentum into ocean
          IF (J.lt.JM-1) DMVI(I,J) = 0.5*DMVI(I,J)*(FOCEAN(I,J)*RSI(I,J)
     *         *DXYN(J)+FOCEAN(I,J+1)*RSI(I,J+1)*DXYS(J+1))/DXYV(J+1)
          IF (J.eq.JM-1) DMVI(I,JM-1) = 0.5*DMVI(I,JM-1)*(FOCEAN(I,JM-1)
     *         *RSI(I,JM-1)*DXYN(JM-1)+FOCEAN(1,JM)*RSI(1,JM)*DXYS(JM))
     *         /DXYV(JM)
        enddo
      enddo
      usi(1:im,1)=0.
      vsi(1:im,jm)=0.
      dmui(1:im,1)=0.
      dmvi(1:im,jm)=0.

C**** Calculate ustar*2*rho for ice-ocean fluxes on atmosphere grid
C**** UI2rho = | tau |
      do j=1,jm
        do i=1,imaxj(j)
          UI2rho(i,j)=0
          if (FOCEAN(I,J)*RSI(i,j).gt.0) THEN
C**** calculate 4 point average of B grid values of stresses
            duA = 0.5*(DXYN(J)*(dmu(i+1,j)+dmu(i,j))+DXYS(j)*(dmu(i+1
     *           ,j-1)+dmu(i,j-1)))*BYDXYP(J)
            dvA = 0.5*(DXYN(J)*(dmv(i+1,j)+dmv(i,j))+DXYS(j)*(dmv(i+1
     *           ,j-1)+dmv(i,j-1)))*BYDXYP(J)
            UI2rho(i,j)= sqrt (duA**2 + dvA**2) * bydts
          end if
        end do
      end do

C**** set north pole in C grid (atm)
      USINP=0.
      DMUINP=0.
      do i=1,im
        USINP = USINP + USI(i,jm)
        DMUINP = DMUINP + DMUI(i,jm)
      enddo
      USINP=USINP/IM
      DMUINP=DMUINP/IM
      USI(1:im,jm)=USINP
      DMUI(1:im,jm)=DMUINP
      VSI(1:im,jm)=0.
      DMVI(1:im,jm)=0.

C**** calculate mass fluxes for the ice advection
      DO J=1,JM-1
        I=IM
        DO IP1=1,IM
          USIDT(I,J)=0.
          IF (FOCEAN(I,J).gt.0 .and. FOCEAN(IP1,J).gt.0. .and.
     *         RSI(I,J)+RSI(IP1,J).gt.1d-4) THEN
            USIDT(I,J)=USI(I,J)*DTS
            ICIJ(I,J,IJ_USI) =ICIJ(I,J,IJ_USI) +(RSI(I,J)+RSI(IP1,J))
     *           *USI(i,j)
            ICIJ(I,J,IJ_DMUI)=ICIJ(I,J,IJ_DMUI)+DMUI(i,j)
          END IF
          VSIDT(I,J)=0.
          IF (FOCEAN(I,J+1).gt.0 .and. FOCEAN(I,J).gt.0. .and.
     *         RSI(I,J)+RSI(I,J+1).gt.1d-4) THEN
            VSIDT(I,J)=VSI(I,J)*DTS
            ICIJ(I,J,IJ_VSI) =ICIJ(I,J,IJ_VSI) +(RSI(I,J)+RSI(I,J+1))
     *           *VSI(i,j)
            ICIJ(I,J,IJ_DMVI)=ICIJ(I,J,IJ_DMVI)+DMVI(i,j)
          END IF
          ICIJ(I,J,IJ_PICE)=ICIJ(I,J,IJ_PICE)+ RSI(I,J)*press(i+1,j)
          ICIJ(I,J,IJ_RSI) =ICIJ(I,J,IJ_RSI) + RSI(I,J)
          I=IP1
        END DO
      END DO
      VSIDT(1:IM,JM)=0.
      USIDT(1:IM,JM)=USI(1,JM)*DTS
      ICIJ(1,JM,IJ_USI) =ICIJ(1,JM,IJ_USI) +RSI(1,JM)*USI(1,JM)
      ICIJ(1,JM,IJ_DMUI)=ICIJ(1,JM,IJ_DMUI)+DMUI(1,JM)
      ICIJ(1,JM,IJ_RSI) =ICIJ(1,JM,IJ_RSI) +RSI(1,JM)
      ICIJ(1,JM,IJ_PICE)=ICIJ(1,JM,IJ_PICE)+RSI(1,JM)*press(1,JM)
C**** Set uisurf,visurf for use in atmospheric drag calculations
      uisurf=usi
      visurf=vsi
C****
      END SUBROUTINE DYNSI

      SUBROUTINE ADVSI
!@sum  ADVSI advects sea ice
!@+    Currently set up to advect ice on AGCM grid (i.e. usidt/vsidt are
!@+    on the AGCM grid, and RSI/MSI/HSI etc. are unchanged)
!@+    At some point this will change (USIDT/VSIDT on ice grid, and RSI
!@+    etc. will need to be interpolated back and forth).
!@auth Gary Russell/Gavin Schmidt
      USE CONSTANT, only : byshi,lhm,grav
      USE MODEL_COM, only : im,jm,focean,p,ptop,kocean
      USE GEOM, only : dxyp,dyp,dxp,dxv,bydxyp,imaxj
c      USE ICEGEOM, only : dxyp,dyp,dxp,dxv,bydxyp ?????
      USE ICEDYN_COM, only : usidt,vsidt,rsix,rsiy,rsisave,icij,ij_musi
     *     ,ij_mvsi,ij_husi,ij_hvsi,ij_susi,ij_svsi
#ifdef TRACERS_WATER
     *     ,ticij,ticij_tusi,ticij_tvsi
#endif
      USE SEAICE, only : ace1i,xsi
      USE SEAICE_COM, only : rsi,msi,snowi,hsi,ssi,lmi
#ifdef TRACERS_WATER
     *     ,trsi,ntm
#endif
      USE FLUXES, only : gtemp,apress,msicnv,fwsim
#ifdef TRACERS_WATER
     *     ,gtracer
#endif
      USE DAGCOM, only : oa
      IMPLICIT NONE
      REAL*8, DIMENSION(IM) :: FAW,FASI,FXSI,FYSI
!@var NTRICE max. number of tracers to be advected (mass/heat/salt+)
#ifndef TRACERS_WATER
      INTEGER, PARAMETER :: NTRICE=2+2*LMI
#else
      INTEGER, PARAMETER :: NTRICE=2+(2+NTM)*LMI
      INTEGER ITR
#endif
      REAL*8 FMSI(NTRICE,IM),SFMSI(NTRICE),AMSI(NTRICE),BYFOA(IM,JM)
      INTEGER I,J,L,IM1,IP1,K
      REAL*8 SFASI,DMHSI,ASI,YRSI,XRSI,FRSI
!@var MHS mass/heat/salt content of sea ice
      REAL*8, DIMENSION(NTRICE,IM,JM) :: MHS
C****
C**** FLUXCB  USIDT  U compon of time integrated sea ice velocity (m)
C****         VSIDT  V compon of time integrated sea ice velocity (m)
C****
C****         FAW    flux of surface water area (m^2) = USIDT*DYP
C****         FASI   flux of sea ice area (m^2) = USIDT*DYP*RSIedge
C****         FMSI   flux of sea ice mass (kg) or heat (J) or salt (kg)

C**** Regularise ice concentration gradients to prevent advection errors
      DO J=2,JM-1
      DO I=1,IM
        IF (RSI(I,J).gt.1d-4) THEN
          IF (RSISAVE(I,J).gt.RSI(I,J)) THEN ! reduce gradients
            FRSI=(RSISAVE(I,J)-RSI(I,J))/RSISAVE(I,J)
            RSIX(I,J)=RSIX(I,J)*(1.-FRSI)
            RSIY(I,J)=RSIY(I,J)*(1.-FRSI)
          END IF
          IF(RSI(I,J)-RSIX(I,J).lt.0.)  RSIX(I,J) =    RSI(I,J)
          IF(RSI(I,J)+RSIX(I,J).lt.0.)  RSIX(I,J) =   -RSI(I,J)
          IF(RSI(I,J)-RSIX(I,J).gt.1d0) RSIX(I,J) =    RSI(I,J)-1d0
          IF(RSI(I,J)+RSIX(I,J).gt.1d0) RSIX(I,J) =1d0-RSI(I,J)
          IF(RSI(I,J)-RSIY(I,J).lt.0.)  RSIY(I,J) =    RSI(I,J)
          IF(RSI(I,J)+RSIY(I,J).lt.0.)  RSIY(I,J) =   -RSI(I,J)
          IF(RSI(I,J)-RSIY(I,J).gt.1d0) RSIY(I,J) =    RSI(I,J)-1d0
          IF(RSI(I,J)+RSIY(I,J).gt.1d0) RSIY(I,J) =1d0-RSI(I,J)
        ELSE
          RSIX(I,J) = 0.  ; RSIY(I,J) = 0.
        END IF
      END DO
      END DO

C**** set up local MHS array to contain all advected quantities
C**** MHS(1:2) = MASS, MHS(3:6) = HEAT, MHS(7:10)=SALT
C**** Currently this is on atmospheric grid
      MHS(1,:,:) = ACE1I + SNOWI
      MHS(2,:,:) = MSI
      DO L=1,LMI
        MHS(L+2,:,:) = HSI(L,:,:)
        MHS(L+2+LMI,:,:) = SSI(L,:,:)
#ifdef TRACERS_WATER
C**** add tracers to advected arrays
        DO ITR=1,NTM
          MHS(L+2+(1+ITR)*LMI,:,:)=TRSI(ITR,L,:,:)
        END DO
#endif
      END DO
C**** define inverse area array
      DO J=1,JM
      DO I=1,IM
        IF (FOCEAN(I,J).gt.0) THEN
          BYFOA(I,J)=BYDXYP(J)/FOCEAN(I,J)
        ELSE
          BYFOA(I,J)=0.
        END IF
      END DO
      END DO
C****
C**** North-South Advection of Sea Ice
C****
      SFASI  = 0.
      SFMSI(1:NTRICE) = 0.
      DO 340 I=1,IM
C****
C**** Calculate south-north sea ice fluxes at grid box edges
C****
      DO 120 J=2,JM-2
      IF(VSIDT(I,J).eq.0.)  GO TO 120
      FAW(J) = VSIDT(I,J)*DXV(J+1) ! be careful with atm. grid index
      IF(VSIDT(I,J).le.0.) THEN
C**** Sea ice velocity is southward at grid box edge
        FASI(J)=FAW(J)*(RSI(I,J+1)-(1d0+FAW(J)*BYDXYP(J+1))*RSIY(I,J+1))
     *       *FOCEAN(I,J+1)
        FXSI(J)=FAW(J)*RSIX(I,J+1)*FOCEAN(I,J+1)
        FYSI(J)=FAW(J)*(FAW(J)*BYDXYP(J+1)*
     *       FAW(J)*RSIY(I,J+1)*FOCEAN(I,J+1) - 3d0*FASI(J))
        FMSI(1:NTRICE,J) = FASI(J)*MHS(1:NTRICE,I,J+1)
      ELSE
C**** Sea ice velocity is northward at grid box edge
        FASI(J)=FAW(J)*(RSI(I,J)+(1d0-FAW(J)*BYDXYP(J))*RSIY(I,J))
     *       *FOCEAN(I,J)
        FXSI(J)=FAW(J)*RSIX(I,J)*FOCEAN(I,J)
        FYSI(J)=FAW(J)*(FAW(J)*BYDXYP(J)*FAW(J)*RSIY(I,J)*FOCEAN(I,J)
     *       -3d0*FASI(J))
        FMSI(1:NTRICE,J) = FASI(J)*MHS(1:NTRICE,I,J)
      END IF
         ICIJ(I,J,IJ_MVSI)=ICIJ(I,J,IJ_MVSI)+SUM(FMSI(1:2,J))
         ICIJ(I,J,IJ_HVSI)=ICIJ(I,J,IJ_HVSI)+SUM(FMSI(3:2+LMI,J))
         ICIJ(I,J,IJ_SVSI)=ICIJ(I,J,IJ_SVSI)+SUM(FMSI(3+LMI:2+2*LMI,J))
#ifdef TRACERS_WATER
         DO ITR=1,NTM
           TICIJ(I,J,TICIJ_TVSI,ITR)=TICIJ(I,J,TICIJ_TVSI,ITR)+
     *          SUM(FMSI(3+(1+ITR)*LMI:2+(2+ITR)*LMI,J))
         END DO
#endif
  120 CONTINUE
C****
C**** Calculate south-north sea ice fluxes near North Pole
C****
      IF(VSIDT(I,JM-1).eq.0.)  GO TO 200
      FAW(JM-1) = VSIDT(I,JM-1)*DXV(JM) ! be careful with atm.grid index
      IF(VSIDT(I,JM-1).le.0.) THEN
C**** Sea ice velocity is southward from North Pole box
        FASI(JM-1) = FAW(JM-1)*RSI(1,JM)*FOCEAN(1,JM)
        FXSI(JM-1) = 0.
        FYSI(JM-1) = -FAW(JM-1)*FASI(JM-1)
        FMSI(1:NTRICE,JM-1) = FASI(JM-1)*MHS(1:NTRICE,1,JM)
      ELSE
C**** Sea ice velocity is northward into North Pole box
        FASI(JM-1) = FAW(JM-1)*FOCEAN(I,JM-1)*
     *       (RSI(I,JM-1)+(1d0-FAW(JM-1)*BYDXYP(JM-1))*RSIY(I,JM-1))
        FXSI(JM-1) = FAW(JM-1)*RSIX(I,JM-1)*FOCEAN(I,JM-1)
        FYSI(JM-1) = FAW(JM-1)*(FAW(JM-1)*BYDXYP(JM-1)*
     *       FAW(JM-1)*RSIY(I,JM-1)*FOCEAN(I,JM-1)-3d0*FASI(JM-1))
        FMSI(1:NTRICE,JM-1) = FASI(JM-1)*MHS(1:NTRICE,I,JM-1)
      END IF
C**** Accumulate sea ice leaving and entering North Pole box
      SFASI = SFASI + FASI(JM-1)
      SFMSI(1:NTRICE) = SFMSI(1:NTRICE) + FMSI(1:NTRICE,JM-1)
       ICIJ(I,JM-1,IJ_MVSI)=ICIJ(I,JM-1,IJ_MVSI)+SUM(FMSI(1:2,JM-1))
       ICIJ(I,JM-1,IJ_HVSI)=ICIJ(I,JM-1,IJ_HVSI)+SUM(FMSI(3:2+LMI,JM-1))
       ICIJ(I,JM-1,IJ_SVSI)=ICIJ(I,JM-1,IJ_SVSI)+
     *      SUM(FMSI(3+LMI:2+2*LMI,JM-1))
#ifdef TRACERS_WATER
         DO ITR=1,NTM
           TICIJ(I,JM-1,TICIJ_TVSI,ITR)=TICIJ(I,JM-1,TICIJ_TVSI,ITR)+
     *          SUM(FMSI(3+(1+ITR)*LMI:2+(2+ITR)*LMI,JM-1))
         END DO
#endif
C****
C**** Update sea ice variables due to south-north fluxes
C****
  200 DO 330 J=2,JM-1
      IF(VSIDT(I,J-1)) 240,210,280
C**** VSIDT(J-1)=0.
  210 IF(VSIDT(I,J))  220,330,230
C**** VSIDT(J-1)=0, VSIDT(J)<0.
  220 ASI = RSI(I,J)*DXYP(J)*FOCEAN(I,J) -  FASI(J)
      DO 225 K=1,NTRICE
  225 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J)*FOCEAN(I,J) - FMSI(K,J)
      IF(ASI.gt.DXYP(J)*FOCEAN(I,J))  GO TO 320
      YRSI = (RSIY(I,J)*DXYP(J)*DXYP(J)*FOCEAN(I,J) - FYSI(J)
     *    + 3d0*(FAW(J)*ASI-DXYP(J)*FASI(J))) / (DXYP(J)-FAW(J))
      RSI(I,J)  = ASI*BYFOA(I,J)
      RSIY(I,J) = YRSI*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J) - FXSI(J)*BYFOA(I,J)
      IF (ASI.gt.0) MHS(1:NTRICE,I,J) = AMSI(1:NTRICE)/ASI
      GO TO 310
C**** VSIDT(J-1)=0, VSIDT(J)>0.
  230 RSI(I,J)  =  RSI(I,J) -  FASI(J)*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J)*(1d0-FAW(J)*BYDXYP(J))
      RSIY(I,J) = RSIY(I,J)*(1d0-FAW(J)*BYDXYP(J))**2
      GO TO 310
C**** VSIDT(J-1)<0.
  240 IF(VSIDT(I,J))  260,250,270
C**** VSIDT(J-1)<0, VSIDT(J)=0.
  250 RSI(I,J)  =  RSI(I,J) +  FASI(J-1)*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J)*(1d0+FAW(J-1)*FOCEAN(I,J-1)*BYFOA(I,J))
      RSIY(I,J) = RSIY(I,J)*(1d0+FAW(J-1)*FOCEAN(I,J-1)*BYFOA(I,J))**2
      GO TO 310
C**** VSIDT(J-1)<0, VSIDT(J)<0  or  VSIDT(J-1)>0, VSIDT(J) not 0.
  260 ASI = RSI(I,J)*DXYP(J)*FOCEAN(I,J) + ( FASI(J-1)- FASI(J))
      DO 265 K=1,NTRICE
  265 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J)*FOCEAN(I,J) +
     *       (FMSI(K,J-1)-FMSI(K,J))
      IF(ASI.gt.DXYP(J)*FOCEAN(I,J))  GO TO 320
      YRSI = (RSIY(I,J)*DXYP(J)*DXYP(J)*FOCEAN(I,J)+(FYSI(J-1)-FYSI(J))
     *    + 3d0*((FAW(J-1)+FAW(J))*ASI-DXYP(J)*(FASI(J-1)+FASI(J))))
     *    / (DXYP(J) + (FAW(J-1)-FAW(J)))
      RSI(I,J)  = ASI*BYFOA(I,J)
      RSIY(I,J) = YRSI*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J) + (FXSI(J-1)-FXSI(J))*BYFOA(I,J)
      IF (ASI.gt.0) MHS(1:NTRICE,I,J) = AMSI(1:NTRICE)/ASI
      GO TO 310
C**** VSIDT(J-1)<0, VSIDT(J)>0.
  270 RSI(I,J)  =  RSI(I,J) + (FASI(J-1)-FASI(J))*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J)*(1d0+(FAW(J-1)*FOCEAN(I,J-1)
     *     -FAW(J)*FOCEAN(I,J))*BYFOA(I,J))
      RSIY(I,J) = RSIY(I,J)*(1d0+(FAW(J-1)*FOCEAN(I,J-1)
     *     -FAW(J)*FOCEAN(I,J))*BYFOA(I,J))**2
      GO TO 310
C**** VSIDT(J-1)>0.
  280 IF(VSIDT(I,J).ne.0.)  GO TO 260
C**** VSIDT(J-1)>0, VSIDT(J)=0.
      ASI = RSI(I,J)*DXYP(J)*FOCEAN(I,J) + FASI(J-1)
      DO 285 K=1,NTRICE
  285 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J)*FOCEAN(I,J) + FMSI(K,J-1)
      IF(ASI.gt.DXYP(J)*FOCEAN(I,J))  GO TO 320
      YRSI = (RSIY(I,J)*DXYP(J)*DXYP(J)*FOCEAN(I,J) + FYSI(J-1)
     *    + 3d0*(FAW(J-1)*ASI-DXYP(J)*FASI(J-1))) / (DXYP(J)+FAW(J-1))
      RSI(I,J)  = ASI*BYFOA(I,J)
      RSIY(I,J) = YRSI*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J) + FXSI(J-1)*BYFOA(I,J)
      IF (ASI.gt.0) MHS(1:NTRICE,I,J) = AMSI(1:NTRICE)/ASI
C**** Limit RSIX and RSIY so that sea ice is positive at the edges
  310 RSI(I,J) = MAX(0d0,RSI(I,J))
      IF(RSI(I,J)-RSIX(I,J).lt.0.)  RSIX(I,J) =    RSI(I,J)
      IF(RSI(I,J)+RSIX(I,J).lt.0.)  RSIX(I,J) =   -RSI(I,J)
      IF(RSI(I,J)-RSIX(I,J).gt.1d0) RSIX(I,J) =    RSI(I,J)-1d0
      IF(RSI(I,J)+RSIX(I,J).gt.1d0) RSIX(I,J) =1d0-RSI(I,J)
      IF(RSI(I,J)-RSIY(I,J).lt.0.)  RSIY(I,J) =    RSI(I,J)
      IF(RSI(I,J)+RSIY(I,J).lt.0.)  RSIY(I,J) =   -RSI(I,J)
      IF(RSI(I,J)-RSIY(I,J).gt.1d0) RSIY(I,J) =    RSI(I,J)-1d0
      IF(RSI(I,J)+RSIY(I,J).gt.1d0) RSIY(I,J) =1d0-RSI(I,J)
      GO TO 330
C**** Sea ice crunches into itself and completely covers grid box
  320 RSI(I,J)   = 1d0
      RSIX(I,J)  = 0.
      RSIY(I,J)  = 0.
      MHS(1,I,J) = AMSI(1)/ASI
      MHS(2,I,J) =(AMSI(1)+AMSI(2))*BYFOA(I,J) - MHS(1,I,J)
      DO K=1,(NTRICE-2)/LMI
        MHS(3+LMI*(K-1),I,J) = AMSI(3+LMI*(K-1)) / ASI
        MHS(4+LMI*(K-1),I,J) = AMSI(4+LMI*(K-1)) / ASI
        DMHSI = (AMSI(3+LMI*(K-1))+AMSI(4+LMI*(K-1))+AMSI(5+LMI*(K-1))
     *       +AMSI(6+LMI*(K-1)))*(BYFOA(I,J) -1d0 / ASI )
        MHS(5+LMI*(K-1),I,J) = AMSI(5+LMI*(K-1)) / ASI +
     *       XSI(3)*DMHSI
        MHS(6+LMI*(K-1),I,J) = AMSI(6+LMI*(K-1)) / ASI +
     *       XSI(4)*DMHSI
      END DO
C**** End of loop over J
  330 CONTINUE
C**** End of loop over I
  340 CONTINUE
C****
C**** Advection of Sea Ice leaving and entering North Pole box
C****
      ASI = RSI(1,JM)*DXYP(JM)*FOCEAN(1,JM) + SFASI/IM
      DO 345 K=1,NTRICE
  345 AMSI(K) = RSI(1,JM)*DXYP(JM)*MHS(K,1,JM)*FOCEAN(1,JM)+SFMSI(K)/IM
      IF(ASI.gt.DXYP(JM)*FOCEAN(1,JM))  GO TO 350
      RSI(1,JM)   = ASI*BYFOA(1,JM)
      IF (ASI.gt.0) MHS(1:NTRICE,1,JM) = AMSI(1:NTRICE)/ASI
      GO TO 400
C**** Sea ice crunches into itself at North Pole box
  350 RSI(1,JM)   = 1d0
      MHS(1,1,JM) = AMSI(1)/ASI
      MHS(2,1,JM) =(AMSI(1)+AMSI(2))*BYFOA(1,JM)-MHS(1,1,JM)
      DO K=1,(NTRICE-2)/LMI
        MHS(3+LMI*(K-1),1,JM) = AMSI(3+LMI*(K-1)) / ASI
        MHS(4+LMI*(K-1),1,JM) = AMSI(4+LMI*(K-1)) / ASI
        DMHSI = (AMSI(3+LMI*(K-1))+AMSI(4+LMI*(K-1))+AMSI(5+LMI*(K-1))
     *       +AMSI(6+LMI*(K-1)))*(BYFOA(1,JM) -1d0/ ASI)
        MHS(5+LMI*(K-1),1,JM) = AMSI(5+LMI*(K-1)) / ASI +
     *       XSI(3)*DMHSI
        MHS(6+LMI*(K-1),1,JM) = AMSI(6+LMI*(K-1)) / ASI +
     *       XSI(4)*DMHSI
      END DO
C****
C**** East-West Advection of Sea Ice
C****
  400 DO 640 J=2,JM-1
C****
C**** Calculate west-east sea ice fluxes at grid box edges
C****
      I=IM
      DO IP1=1,IM
      IF(USIDT(I,J).eq.0.)  GO TO 420
      FAW(I) = USIDT(I,J)*DYP(J)
      IF(USIDT(I,J).le.0.) THEN
C**** Sea ice velocity is westward at grid box edge
        FASI(I)=FAW(I)*(RSI(IP1,J)-(1d0+FAW(I)*BYDXYP(J))*RSIX(IP1,J))
     *       *FOCEAN(IP1,J)
        FXSI(I)=FAW(I)*(FAW(I)*BYDXYP(J)*
     *       FAW(I)*RSIX(IP1,J)*FOCEAN(IP1,J)-3d0*FASI(I))
        FYSI(I)=FAW(I)*RSIY(IP1,J)*FOCEAN(IP1,J)
        FMSI(1:NTRICE,I) = FASI(I)*MHS(1:NTRICE,IP1,J)
      ELSE
C**** Sea ice velocity is eastward at grid box edge
        FASI(I)=FAW(I)*(RSI(I,J)+(1d0-FAW(I)*BYDXYP(J))*RSIX(I,J))
     *       *FOCEAN(I,J)
        FXSI(I)=FAW(I)*(FAW(I)*BYDXYP(J)*FAW(I)*RSIX(I,J)*FOCEAN(I,J)
     *       -3d0*FASI(I))
        FYSI(I)=FAW(I)*RSIY(I,J)*FOCEAN(I,J)
        FMSI(1:NTRICE,I) = FASI(I)*MHS(1:NTRICE,I,J)
      END IF
         ICIJ(I,J,IJ_MUSI)=ICIJ(I,J,IJ_MUSI)+SUM(FMSI(1:2,I))
         ICIJ(I,J,IJ_HUSI)=ICIJ(I,J,IJ_HUSI)+SUM(FMSI(3:2+LMI,I))
         ICIJ(I,J,IJ_SUSI)=ICIJ(I,J,IJ_SUSI)+SUM(FMSI(3+LMI:2+2*LMI,I))
#ifdef TRACERS_WATER
         DO ITR=1,NTM
           TICIJ(I,J,TICIJ_TUSI,ITR)=TICIJ(I,J,TICIJ_TUSI,ITR)+
     *          SUM(FMSI(3+(1+ITR)*LMI:2+(2+ITR)*LMI,I))
         END DO
#endif
  420 I=IP1
      END DO
C****
C**** Update sea ice variables due to west-east fluxes
C****
      IM1=IM
      DO 630 I=1,IM
      IF(USIDT(IM1,J)) 540,510,580
C**** USIDT(IM1)=0.
  510 IF(USIDT(I,J))  520,630,530
C**** USIDT(IM1)=0, USIDT(I)<0.
  520 ASI = RSI(I,J)*DXYP(J)*FOCEAN(I,J) -  FASI(I)
      DO 525 K=1,NTRICE
  525 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J)*FOCEAN(I,J) - FMSI(K,I)
      IF(ASI.gt.DXYP(J)*FOCEAN(I,J))  GO TO 620
      XRSI = (RSIX(I,J)*DXYP(J)*DXYP(J)*FOCEAN(I,J) - FXSI(I)
     *    + 3d0*(FAW(I)*ASI-DXYP(J)*FASI(I))) / (DXYP(J)-FAW(I))
      RSI(I,J)  = ASI*BYFOA(I,J)
      RSIX(I,J) = XRSI*BYFOA(I,J)
      RSIY(I,J) = RSIY(I,J) - FYSI(I)*BYFOA(I,J)
      IF (ASI.gt.0) MHS(1:NTRICE,I,J) = AMSI(1:NTRICE)/ASI
      GO TO 610
C**** USIDT(IM1)=0, USIDT(I)>0.
  530 RSI(I,J)  =  RSI(I,J) -  FASI(I)*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J)*(1d0-FAW(I)*BYDXYP(J))**2
      RSIY(I,J) = RSIY(I,J)*(1d0-FAW(I)*BYDXYP(J))
      GO TO 610
C**** USIDT(IM1)<0.
  540 IF(USIDT(I,J))  560,550,570
C**** USIDT(IM1)<0, USIDT(I)=0.
  550 RSI(I,J)  =  RSI(I,J) +  FASI(IM1)*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J)*(1d0+FAW(IM1)*FOCEAN(IM1,J)*BYFOA(I,J))**2
      RSIY(I,J) = RSIY(I,J)*(1d0+FAW(IM1)*FOCEAN(IM1,J)*BYFOA(I,J))
      GO TO 610
C**** USIDT(IM1)<0, USIDT(I)<0  or  USIDT(IM1)>0, USIDT(I) not 0.
  560 ASI = RSI(I,J)*DXYP(J)*FOCEAN(I,J) + (FASI(IM1)- FASI(I))
      DO 565 K=1,NTRICE
  565 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J)*FOCEAN(I,J) +
     *       (FMSI(K,IM1)-FMSI(K,I))
      IF(ASI.gt.DXYP(J)*FOCEAN(I,J))  GO TO 620
      XRSI = (RSIX(I,J)*DXYP(J)*DXYP(J)*FOCEAN(I,J)+(FXSI(IM1)-FXSI(I))
     *    + 3d0*((FAW(IM1)+FAW(I))*ASI-DXYP(J)*(FASI(IM1)+FASI(I))))
     *    / (DXYP(J) + (FAW(IM1)-FAW(I)))
      RSI(I,J)  = ASI*BYFOA(I,J)
      RSIX(I,J) = XRSI*BYFOA(I,J)
      RSIY(I,J) = RSIY(I,J) + (FYSI(IM1)-FYSI(I))*BYFOA(I,J)
      IF (ASI.gt.0) MHS(1:NTRICE,I,J) = AMSI(1:NTRICE)/ASI
      GO TO 610
C**** USIDT(IM1)<0, USIDT(I)>0.
  570 RSI(I,J)  =  RSI(I,J) + (FASI(IM1)-FASI(I))*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J)*(1d0+(FAW(IM1)*FOCEAN(IM1,J)
     *     -FAW(I)*FOCEAN(I,J))*BYFOA(I,J))**2
      RSIY(I,J) = RSIY(I,J)*(1d0+(FAW(IM1)*FOCEAN(IM1,J)
     *     -FAW(I)*FOCEAN(I,J))*BYFOA(I,J))
      GO TO 610
C**** USIDT(IM1)>0.
  580 IF(USIDT(I,J).ne.0.)  GO TO 560
C**** USIDT(IM1)>0, USIDT(I)=0.
      ASI = RSI(I,J)*DXYP(J)*FOCEAN(I,J) + FASI(IM1)
      DO 585 K=1,NTRICE
  585 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J)*FOCEAN(I,J) + FMSI(K,IM1)
      IF(ASI.gt.DXYP(J)*FOCEAN(I,J))  GO TO 620
      XRSI = (RSIX(I,J)*DXYP(J)*DXYP(J)*FOCEAN(I,J) + FXSI(IM1)
     *    + 3d0*(FAW(IM1)*ASI-DXYP(J)*FASI(IM1))) / (DXYP(J)+FAW(IM1))
      RSI(I,J)  = ASI*BYFOA(I,J)
      RSIX(I,J) = XRSI*BYFOA(I,J)
      RSIY(I,J) = RSIY(I,J) + FYSI(IM1)*BYFOA(I,J)
      IF (ASI.gt.0) MHS(1:NTRICE,I,J) = AMSI(1:NTRICE)/ASI
C**** Limit RSIX and RSIY so that sea ice is positive at the edges
  610 RSI(I,J) = MAX(0d0,RSI(I,J))
      IF(RSI(I,J)-RSIX(I,J).lt.0.)  RSIX(I,J) =    RSI(I,J)
      IF(RSI(I,J)+RSIX(I,J).lt.0.)  RSIX(I,J) =   -RSI(I,J)
      IF(RSI(I,J)-RSIX(I,J).gt.1d0) RSIX(I,J) =    RSI(I,J)-1d0
      IF(RSI(I,J)+RSIX(I,J).gt.1d0) RSIX(I,J) =1d0-RSI(I,J)
      IF(RSI(I,J)-RSIY(I,J).lt.0.)  RSIY(I,J) =    RSI(I,J)
      IF(RSI(I,J)+RSIY(I,J).lt.0.)  RSIY(I,J) =   -RSI(I,J)
      IF(RSI(I,J)-RSIY(I,J).gt.1d0) RSIY(I,J) =    RSI(I,J)-1d0
      IF(RSI(I,J)+RSIY(I,J).gt.1d0) RSIY(I,J) =1d0-RSI(I,J)
      GO TO 630
C**** Sea ice crunches into itself and completely covers grid box
  620 RSI(I,J)   = 1d0
      RSIX(I,J)  = 0.
      RSIY(I,J)  = 0.
      MHS(1,I,J) = AMSI(1)/ASI
      MHS(2,I,J) =(AMSI(1)+AMSI(2))*BYFOA(I,J) - MHS(1,I,J)
      DO K=1,(NTRICE-2)/LMI
        MHS(3+LMI*(K-1),I,J) = AMSI(3+LMI*(K-1)) / ASI
        MHS(4+LMI*(K-1),I,J) = AMSI(4+LMI*(K-1)) / ASI
        DMHSI = (AMSI(3+LMI*(K-1))+AMSI(4+LMI*(K-1))+AMSI(5+LMI*(K-1))
     *       +AMSI(6+LMI*(K-1)))*(BYFOA(I,J) -1d0/ ASI)
        MHS(5+LMI*(K-1),I,J) = AMSI(5+LMI*(K-1)) / ASI +
     *       XSI(3)*DMHSI
        MHS(6+LMI*(K-1),I,J) = AMSI(6+LMI*(K-1)) / ASI +
     *       XSI(4)*DMHSI
      END DO
C**** End of loop over I
  630 IM1=I
C**** End of loop over J
  640 CONTINUE

      IF (KOCEAN.eq.1) THEN ! full ocean calculation, adjust sea ice
C**** set global variables from local array
C**** Currently on atmospheric grid, so no interpolation necessary
        DO J=1,JM
          DO I=1,IMAXJ(J)
C**** Fresh water sea ice mass convergence (needed for qflux model)
            MSICNV(I,J) = RSI(I,J)*(MHS(1,I,J)+MHS(2,I,J)-SUM(MHS(3
     *           +LMI:2*LMI+2,I,J))) - RSISAVE(I,J)*(ACE1I+SNOWI(I,J)
     *           +MSI(I,J)-SUM(SSI(:,I,J)))
C**** sea ice prognostic variables
            SNOWI(I,J)= MAX(0d0,MHS(1,I,J) - ACE1I)
            MSI(I,J)  = MHS(2,I,J)
            DO L=1,LMI
              HSI(L,I,J) = MHS(L+2,I,J)
              SSI(L,I,J) = MHS(L+2+LMI,I,J)
#ifdef TRACERS_WATER
              DO ITR=1,NTM
                TRSI(ITR,L,I,J)=MHS(L+2+(1+ITR)*LMI,I,J)
              END DO
#endif
            END DO
            FWSIM(I,J)=RSI(I,J)*(ACE1I+SNOWI(I,J)+MSI(I,J)-
     *           SUM(SSI(1:LMI,I,J)))
          END DO
        END DO
C**** Set atmospheric arrays
        DO J=1,JM
          DO I=1,IMAXJ(J)
            IF (FOCEAN(I,J).gt.0) THEN
C**** set total atmopsheric pressure anomaly in case needed by ocean
              APRESS(I,J) = 100.*(P(I,J)+PTOP-1013.25d0)+RSI(I,J)
     *             *(SNOWI(I,J)+ACE1I+MSI(I,J))*GRAV
              GTEMP(1:2,2,I,J)=((HSI(1:2,I,J)-SSI(1:2,I,J)*LHM)/
     *             (XSI(1:2)*MHS(1,I,J))+LHM)*BYSHI
#ifdef TRACERS_WATER
              GTRACER(:,2,I,J)=TRSI(:,1,I,J)/(XSI(1)*MHS(1,I,J)
     *             -SSI(1,I,J))
#endif
            END IF
          END DO
        END DO
        DO I=2,IM               ! North pole
          APRESS(I,JM)=APRESS(1,JM)
          GTEMP(1:2,2,I,JM)= GTEMP(1:2,2,1,JM)
#ifdef TRACERS_WATER
          GTRACER(:,2,I,JM)=GTRACER(:,2,1,JM)
#endif
        END DO
      ELSE          ! fixed SST case, save implied heat convergence
        DO J=1,JM
          DO I=1,IMAXJ(J)
            IF (FOCEAN(I,J).gt.0) THEN
              OA(I,J,13)=OA(I,J,13)+(RSI(I,J)*SUM(MHS(3:2+LMI,I,J))
     *             -RSISAVE(I,J)*SUM(HSI(1:LMI,I,J)))
C**** reset sea ice concentration
              RSI(I,J)=RSISAVE(I,J)
            END IF
          END DO
        END DO
      END IF
C****
      RETURN
      END SUBROUTINE ADVSI

      SUBROUTINE init_icedyn(iniOCEAN)
!@sum  init_icedyn initializes ice dynamics variables
!@auth Gavin Schmidt
      USE MODEL_COM, only : im,jm,dtsrc,foceanA=>focean
      USE DAGCOM, only : ia_src
      USE ICEDYN_COM
      USE ICEDYN, only : setup_icedyn_grid,focean,osurf_tilt
      USE FLUXES, only : uisurf,visurf
      USE PARAM
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: iniOCEAN
      INTEGER k,i,j

C**** setup ice dynamics grid
C**** Currently using ATM grid:
C**** -------------------------
C**** If a different grid (but still lat/lon) is required, edit
C**** definition of focean here and the definition of imic,jmic in
C**** ICEDYN.f, and obviously the code in DYNSI.
C**** -------------------------
C**** Set land masks for ice dynamics
      focean(:,:)=foceanA(:,:)         ! EDIT FOR GRID CHANGE

      call setup_icedyn_grid

C**** Initiallise ice dynamics if ocean model needs initialising
      if (iniOCEAN) THEN
        RSIX=0.
        RSIY=0.
        USI=0.
        VSI=0.
      end if

C**** set uisurf,visurf for atmopsherice drag calculations
      uisurf=usi
      visurf=vsi

C**** set properties for ICIJ diagnostics
      k=0

      k=k+1
      IJ_USI=k
      lname_icij(k)="Sea ice EW velocity x POICEU"
      sname_icij(k)="icij_usi"
      units_icij(k)="m/s"
      ia_icij(k)=ia_src
      scale_icij(k)=1.
      ijgrid_icij(k)=2

      k=k+1
      IJ_VSI=k
      lname_icij(k)="Sea ice NS velocity x POICEV"
      sname_icij(k)="icij_vsi"
      units_icij(k)="m/s"
      ia_icij(k)=ia_src
      scale_icij(k)=1.
      ijgrid_icij(k)=2

      k=k+1
      IJ_DMUI=k
      lname_icij(k)="Ice-ocean EW stress"
      sname_icij(k)="icij_dmui"
      units_icij(k)="kg/m s^2"
      ia_icij(k)=ia_src
      scale_icij(k)=1./dtsrc
      ijgrid_icij(k)=1

      k=k+1
      IJ_DMVI=k
      lname_icij(k)="Ice-ocean NS stress"
      sname_icij(k)="icij_dmvi"
      units_icij(k)="kg/m s^2"
      ia_icij(k)=ia_src
      scale_icij(k)=1./dtsrc
      ijgrid_icij(k)=1

      k=k+1
      IJ_PICE=k
      lname_icij(k)="Sea ice internal pressure x POICE"
      sname_icij(k)="icij_psi"
      units_icij(k)="10^3 kg/m s^2"
      ia_icij(k)=ia_src
      scale_icij(k)=1d-3
      ijgrid_icij(k)=1

      k=k+1
      IJ_MUSI=k
      lname_icij(k)="Sea ice NS mass flux"
      sname_icij(k)="icij_musi"
      units_icij(k)="10^7 kg/s"
      ia_icij(k)=ia_src
      scale_icij(k)=1d-7/dtsrc
      ijgrid_icij(k)=2

      k=k+1
      IJ_MVSI=k
      lname_icij(k)="Sea ice EW mass flux"
      sname_icij(k)="icij_mvsi"
      units_icij(k)="10^7 kg/s"
      ia_icij(k)=ia_src
      scale_icij(k)=1d-7/dtsrc
      ijgrid_icij(k)=2

      k=k+1
      IJ_HUSI=k
      lname_icij(k)="Sea ice NS heat flux"
      sname_icij(k)="icij_husi"
      units_icij(k)="10^12 W"
      ia_icij(k)=ia_src
      scale_icij(k)=1d-12/dtsrc
      ijgrid_icij(k)=2

      k=k+1
      IJ_HVSI=k
      lname_icij(k)="Sea ice EW heat flux"
      sname_icij(k)="icij_hvsi"
      units_icij(k)="10^12 W"
      ia_icij(k)=ia_src
      scale_icij(k)=1d-12/dtsrc
      ijgrid_icij(k)=2

      k=k+1
      IJ_SUSI=k
      lname_icij(k)="Sea ice NS salt flux"
      sname_icij(k)="icij_susi"
      units_icij(k)="10^3 kg/s"
      ia_icij(k)=ia_src
      scale_icij(k)=1d-3/dtsrc
      ijgrid_icij(k)=2

      k=k+1
      IJ_SVSI=k
      lname_icij(k)="Sea ice EW salt flux"
      sname_icij(k)="icij_svsi"
      units_icij(k)="10^3 kg/s"
      ia_icij(k)=ia_src
      scale_icij(k)=1d-3/dtsrc
      ijgrid_icij(k)=2

      k=k+1
      IJ_RSI=k
      lname_icij(k)="Ocean ice fraction (ice dynamic grid)"
      sname_icij(k)="icij_rsi"
      units_icij(k)="%"
      ia_icij(k)=ia_src
      scale_icij(k)=100.
      ijgrid_icij(k)=1

      if (k.gt.KICIJ) then
        write(6,*) "Too many ICIJ diags: increase KICIJ to at least"
     *       ,k
        call stop_model("ICIJ diagnostic error",255)
      end if

#ifdef TRACERS_WATER
C**** simple tracer diags same description for all tracers
C**** set properties for TICIJ diagnostics
      k=0

      k=k+1
      TICIJ_TUSI=k
      lname_ticij(k)="Sea ice NS tracer flux"
      sname_ticij(k)="ticij_tusi"
      units_ticij(k)="kg/s"
      ia_ticij(k)=ia_src
      scale_ticij(k)=1./dtsrc
      ijgrid_ticij(k)=2

      k=k+1
      TICIJ_TVSI=k
      lname_ticij(k)="Sea ice EW tracer flux"
      sname_ticij(k)="ticij_tvsi"
      units_ticij(k)="kg/s"
      ia_ticij(k)=ia_src
      scale_ticij(k)=1./dtsrc
      ijgrid_ticij(k)=2

      if (k.gt.KTICIJ) then
        write(6,*) "Too many TICIJ diags: increase KTICIJ to at least"
     *       ,k
        call stop_model("TICIJ diagnostic error",255)
      end if
#endif

      RETURN
      END SUBROUTINE init_icedyn

      SUBROUTINE diag_ICEDYN
!@sum  diag_ICEDYN prints out diagnostics for ice dynamics
!@auth Gavin Schmidt
      USE CONSTANT, only : undef,teeny
      USE MODEL_COM, only : xlabel,lrunid,jmon0,jyear0,idacc,jdate0
     *     ,amon0,jdate,amon,jyear
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm,trname,ntrocn
#endif
      USE DAGCOM, only : qdiag,acc_period
      USE ICEDYN_COM
      USE ICEDYN, only : focean
      USE FILEMANAGER, only : openunit
      IMPLICIT NONE
      REAL*8, DIMENSION(IMIC,JMIC) :: Q,ADENOM
      INTEGER I,J,L,N,KXLB,ijgrid,IP1,k1,k
      CHARACTER XLB*30
      CHARACTER TITLE*80,lname*50,sname*30,units*50
      character*50 :: unit_string
      REAL*8 QJ(JM),QSUM
      REAL*8 byiacc

      IF (.not. QDIAG) RETURN
C**** determine label to be added to all titles
      KXLB = INDEX(XLABEL(1:11),'(')-1
      IF(KXLB.le.0) KXLB = 10
      XLB = ' '
      XLB(1:13)=acc_period(1:3)//' '//acc_period(4:12)
      XLB = TRIM(XLB)//" "//XLABEL(1:KXLB)

C**** Open output files
      call open_ij(trim(acc_period)//'.icij'//
     *     XLABEL(1:LRUNID),imic,jmic)
C****
C**** Simple scaled ICIJ diagnostics
C****
      DO K=1,KICIJ
        byiacc=1./(IDACC(IA_ICIJ(K))+teeny)
        adenom=1.
        lname=lname_icij(k)
        k1 = index(lname,' x ')
        if (k1 .gt. 0) then
          if (index(lname,' x POICE ') .gt. 0) then
            do j=1,jmic
            do i=1,imic
              adenom(i,j)=icij(i,j,ij_rsi) * byiacc
            end do
            end do
          else if (index(lname,' x POICEU') .gt. 0) then
            do j=1,jmic
            i=imic
            do ip1=1,imic
              adenom(i,j)=0.5*(icij(i,j,ij_rsi)+icij(ip1,j,ij_rsi))
     *             * byiacc
              i=ip1
            end do
            end do
          else if (index(lname,' x POICEV') .gt. 0) then
            do j=1,jmic-1
            do i=1,imic
              adenom(i,j)=0.5*(icij(i,j,ij_rsi)+icij(i,j+1,ij_rsi))
     *             * byiacc
            end do
            end do
          end if
          lname(k1:50) = ' '
        end if

        Q=UNDEF ; QJ=UNDEF ; QSUM=UNDEF
        DO J=1,JMIC
          DO I=1,IMIC
            IF (ADENOM(I,J).gt.0 .and. FOCEAN(I,J).gt.0.5)
     *           Q(I,J)=SCALE_ICIJ(K)*ICIJ(I,J,K)*byiacc/adenom(i,j)
          END DO
        END DO
        Q(2:IMIC,JMIC)=Q(1,JMIC)
        Q(2:IMIC,1)=Q(1,1)
        TITLE=trim(LNAME)//" ("//trim(UNITS_ICIJ(K))//") "
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME_ICIJ(K),LNAME_ICIJ(K),UNITS_ICIJ(K),Q
     *       ,QJ,QSUM,IJGRID_ICIJ(K))

      END DO

#ifdef TRACERS_WATER
C**** simple tracer diags (no need for weighting)
C**** Name and scale are tracer dependent
      DO K=1,KTICIJ
        byiacc=1./(IDACC(IA_TICIJ(K))+teeny)

        DO N=1,NTM
        lname=trim(trname(n))//" "//lname_ticij(k)
        sname=trim(trname(n))//"_"//sname_ticij(k)
        Q=UNDEF
        DO J=1,JMIC
          DO I=1,IMIC
            IF (FOCEAN(I,J).gt.0.5) Q(I,J)=10**(-ntrocn(n))*
     *           SCALE_TICIJ(K)*TICIJ(I,J,K,N)*byiacc
          END DO
        END DO
        Q(2:IMIC,JMIC)=Q(1,JMIC)
        Q(2:IMIC,1)=Q(1,1)
        UNITS=unit_string(ntrocn(n),UNITS_TICIJ(K))
        TITLE=trim(LNAME)//" ("//trim(UNITS)//")"
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID_TICIJ(K))
      END DO
      END DO
#endif
      call close_ij
C****
      RETURN
      END SUBROUTINE diag_ICEDYN

