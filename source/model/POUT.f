!@sum  POUT default output routines for standard GISS formats
!@auth Gavin Schmidt
!@ver  1.0

C****
C**** If other formats are desired, please replace these routines
C**** with ones appropriate for your chosen output format, but with the
C**** same interface
C****
C**** Note: it would be nice to amalgamate IL and JL, but that will
C**** have to wait.

      module gissout
!@sum gissout contains variables for outputting GISS format binaries
!@auth G. Schmidt
      implicit none
!@var iu_ij,iu_jl,iu_il,iu_j !  units for selected diag. output
      integer iu_ij,iu_ijk,iu_il,iu_j,iu_jl
!@var im,jm,lm,lm_req local dimensions set in open_* routines
      integer :: im,jm,lm,lm_req
!@var JMMAX maximum conceivable JM
      INTEGER, PARAMETER :: JMMAX=200
!@var LAT_DG latitude of mid points of primary and sec. grid boxs (deg)
      REAL*8, DIMENSION(JMMAX,2) :: LAT_DG

      end module gissout

      subroutine open_ij(filename,im_gcm,jm_gcm)
!@sum  OPEN_IJ opens the lat-lon binary output file
!@auth M. Kelley
!@ver  1.0
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
!@var IM_GCM,JM_GCM dimensions for ij output
      INTEGER, INTENT(IN) :: im_gcm,jm_gcm

C**** open unit for ij
      call openunit(filename,iu_ij,.true.,.false.)

C**** set dimensions
      im=im_gcm
      jm=jm_gcm

      return
      end subroutine open_ij

      subroutine close_ij
!@sum  CLOSE_IJ closes the lat-lon binary output file
!@auth M. Kelley
!@ver  1.0
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
      call closeunit(iu_ij)
      return
      end subroutine close_ij

      subroutine POUT_IJ(TITLE,SNAME,LNAME,UNITS,XIJ,XJ,XSUM,IJGRID)
!@sum  POUT_IJ output lat-lon binary records
!@auth Gavin Schmidt
!@ver  1.0
      USE GISSOUT
      IMPLICIT NONE
!@var TITLE 80 byte title including description and averaging period
      CHARACTER, INTENT(IN) :: TITLE*80
!@var SNAME short name of field
      CHARACTER, INTENT(IN) :: SNAME*30
!@var LNAME long name of field
      CHARACTER, INTENT(IN) :: LNAME*50
!@var UNITS units of field
      CHARACTER, INTENT(IN) :: UNITS*50
!@var XIJ lat/lon output field
      REAL*8, DIMENSION(IM,JM), INTENT(IN) :: XIJ
!@var XJ lat sum/mean of output field
      REAL*8, DIMENSION(JM), INTENT(IN) :: XJ
!@var XSUM global sum/mean of output field
      REAL*8, INTENT(IN) :: XSUM
!@var IJGRID = 1 for primary lat-lon grid, 2 for secondary lat-lon grid
      INTEGER, INTENT(IN) :: IJGRID

      WRITE(iu_ij) TITLE,REAL(XIJ,KIND=4),REAL(XJ,KIND=4),
     *             REAL(XSUM,KIND=4)
      return
      end subroutine POUT_IJ

      subroutine open_jl(filename,jm_gcm,lm_gcm,lm_req_gcm,lat_dg_gcm)
!@sum  OPEN_JL opens the lat-height binary output file
!@auth M. Kelley
!@ver  1.0
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
!@var LM_GCM,JM_GCM,lm_req_gcm dimensions for jl output
      INTEGER, INTENT(IN) :: lm_gcm,jm_gcm,lm_req_gcm
!@var lat_dg_gcm latitude of mid points of grid boxs (deg)
      REAL*8, INTENT(IN), DIMENSION(JM_GCM,2) :: lat_dg_gcm

      call openunit(filename,iu_jl,.true.,.false.)

C**** set dimensions
      jm=jm_gcm
      lm=lm_gcm
      lm_req=lm_req_gcm
      lat_dg(1:JM,:)=lat_dg_gcm(1:JM,:)

      return
      end subroutine open_jl

      subroutine close_jl
!@sum  CLOSE_JL closes the lat-height binary output file
!@auth M. Kelley
!@ver  1.0
      USE GISSOUT
      IMPLICIT NONE
      close(iu_jl)
      return
      end subroutine close_jl

      subroutine POUT_JL(TITLE,LNAME,SNAME,UNITS_IN,
     &     J1,KLMAX,XJL,PM,CX,CY)
!@sum  POUT_JL output lat-height binary records
!@auth Gavin Schmidt
!@ver  1.0
      USE GISSOUT
      IMPLICIT NONE
!@var TITLE 80 byte title including description and averaging period
      CHARACTER, INTENT(IN) :: TITLE*80
!@var LNAME long name of field
      CHARACTER, INTENT(IN) :: LNAME*50
!@var SNAME short name of field
      CHARACTER, INTENT(IN) :: SNAME*30
!@var UNITS units of field
      CHARACTER, INTENT(IN) :: UNITS_IN*50
!@var KLMAX max level to output
!@var J1 minimum j value to output (needed for secondary grid fields)
      INTEGER, INTENT(IN) :: KLMAX,J1
!@var XJL output field
!@+       (J1:JM,1:KLMAX) is field
!@+       (JM+1:JM+3,1:KLMAX) are global/NH/SH average over L
!@+       (J1:JM+3,LM+LM_REQ+1) are averages over J
      REAL*8, DIMENSION(JM+3,LM+LM_REQ+1), INTENT(IN) :: XJL
!@var PM pressure levels (MB)
      REAL*8, DIMENSION(LM+LM_REQ), INTENT(IN) :: PM

      CHARACTER*16, INTENT(IN) :: CX,CY
      CHARACTER*16, PARAMETER :: CBLANK = '                '
      REAL*8 XCOOR(JM)
      INTEGER J,L,JXMAX

      JXMAX = JM-J1+1
      XCOOR(1:JXMAX) = LAT_DG(J1:JM,J1)

      WRITE (iu_jl) TITLE,JXMAX,KLMAX,1,1,
     *     ((REAL(XJL(J1+J-1,L),KIND=4),J=1,JXMAX),L=1,KLMAX)
     *     ,(REAL(XCOOR(J),KIND=4),J=1,JXMAX)
     *     ,(REAL(PM(L),KIND=4),L=1,KLMAX)
     *     ,REAL(1.,KIND=4),REAL(1.,KIND=4)
     *     ,CX,CY,CBLANK,CBLANK,'NASAGISS'
     *     ,(REAL(XJL(J,LM+LM_REQ+1),KIND=4),J=J1,JM+3)
     *     ,((REAL(XJL(J,L),KIND=4),J=JM+1,JM+3),L=1,KLMAX)

      return
      end

      subroutine open_il(filename,im_gcm,lm_gcm,lm_req_gcm)
!@sum  OPEN_IL opens the lon-height binary output file
!@auth M. Kelley
!@ver  1.0
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
!@var IM_GCM,LM_GCM,lm_req_gcm dimensions for il output
      INTEGER, INTENT(IN) :: im_gcm,lm_gcm,lm_req_gcm

      call openunit(filename,iu_il,.true.,.false.)

C**** set units
      im=im_gcm
      lm=lm_gcm
      lm_req=lm_req_gcm

      return
      end subroutine open_il

      subroutine close_il
!@sum  CLOSE_IL closes the lon-height binary output file
!@auth M. Kelley
!@ver  1.0
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
      call closeunit(iu_il)
      return
      end subroutine close_il

      subroutine POUT_IL(TITLE,sname,lname,unit,I1,ISHIFT,KLMAX,XIL
     *     ,PM,CX,CY,ASUM,GSUM,ZONAL)
!@sum  POUT_IL output lon-height binary records
!@auth Gavin Schmidt
!@ver  1.0
      USE GISSOUT
      USE GEOM, only : lon_dg
      IMPLICIT NONE
!@var TITLE 80 byte title including description and averaging period
      CHARACTER, INTENT(IN) :: TITLE*80
!@var KLMAX max level to output
!@var ISHIFT flag for secondary grid
!@var I1 coordinate index associated with first long. (for wrap-around)
      INTEGER, INTENT(IN) :: KLMAX,ISHIFT,I1
!@var XIL output field
      REAL*8, DIMENSION(IM,LM+LM_REQ+1), INTENT(IN) :: XIL
!@var PM pressure levels (MB)
      REAL*8, DIMENSION(LM+LM_REQ), INTENT(IN) :: PM
!@var ASUM vertical mean/sum
      REAL*8, DIMENSION(IM), INTENT(IN) :: ASUM
!@var GSUM total sum/mean
      REAL*8, INTENT(IN) :: GSUM
!@var ZONAL zonal sum/mean
      REAL*8, DIMENSION(LM+LM_REQ), INTENT(IN) :: ZONAL

      character(len=20), intent(in) :: sname,unit
      character(len=80), intent(in) :: lname
      CHARACTER*16, INTENT(IN) :: CX,CY
      CHARACTER*16, PARAMETER :: CBLANK = ' '
      REAL*8 XCOOR(IM)
      INTEGER I,L

C**** Allow for the possibility of wrap-around arrays
      XCOOR(1:IM-I1+1) = LON_DG(I1:IM,ISHIFT)
      IF (I1.gt.1) XCOOR(IM-I1+2:IM) = LON_DG(1:I1-1,ISHIFT)

      WRITE (iu_il) TITLE,IM,KLMAX,1,1,
     *     ((REAL(XIL(I,L),KIND=4),I=1,IM),L=1,KLMAX)
     *     ,(REAL(XCOOR(I),KIND=4),I=1,IM)
     *     ,(REAL(PM(L),KIND=4),L=1,KLMAX)
     *     ,REAL(0.,KIND=4),REAL(0.,KIND=4)
     *     ,CX,CY,CBLANK,CBLANK,'NASAGISS'
     *     ,(REAL(ASUM(I),KIND=4),I=1,IM),REAL(GSUM,KIND=4)
     *     ,(REAL(ZONAL(L),KIND=4),L=1,KLMAX)

      return
      end subroutine POUT_IL

      subroutine close_j
!@sum  CLOSE_J closes the latitudinal budget-page ascii output file
!@auth M. Kelley
!@ver  1.0
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
      call closeunit(iu_j)
      return
      end subroutine close_j

      subroutine open_j(filename,ntypes,jm_gcm,lat_dg_gcm)
!@sum  OPEN_J opens the latitudinal budget-page ascii output file
!@auth M. Kelley
!@ver  1.0
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
!@var ntypes number of surface types to be output
      integer, intent(in) :: ntypes
!@var JM_GCM dimensions for j output
      INTEGER, INTENT(IN) :: jm_gcm
!@var lat_dg_gcm latitude of mid points of grid boxs (deg)
      REAL*8, INTENT(IN), DIMENSION(JM_GCM,2) :: lat_dg_gcm

      call openunit(filename,iu_j,.false.,.false.)

C**** set dimensions
      jm=jm_gcm
      lat_dg(1:JM,:)=lat_dg_gcm(1:JM,:)

      return
      end subroutine open_j

      subroutine POUT_J(TITLE,SNAME,LNAME,UNITS,BUDG,KMAX,TERRAIN,
     *     iotype)
!@sum  POUT_J output zonal budget ascii file (aplot format)
!@auth Gavin Schmidt
!@ver  1.0
      USE GISSOUT
      USE DAGCOM, only : KAJ
      IMPLICIT NONE
      CHARACTER*16, DIMENSION(KAJ),INTENT(INOUT) :: TITLE
      CHARACTER*16 :: NEWTIT
!@var LNAME,SNAME,UNITS dummy strings
      CHARACTER*50, DIMENSION(KAJ),INTENT(IN) :: LNAME
      CHARACTER*30, DIMENSION(KAJ),INTENT(IN) :: SNAME
      CHARACTER*50, DIMENSION(KAJ),INTENT(IN) :: UNITS
      CHARACTER*16, INTENT(IN) :: TERRAIN
      REAL*8, DIMENSION(JM+3,KAJ), INTENT(IN) :: BUDG
      INTEGER, INTENT(IN) :: KMAX,iotype
      INTEGER K,N,J,n1

C**** Convert spaces in TITLE to underscore
C**** Try simply removing spaces for compactness
      DO K=1,KMAX
        newtit=' '
        n1=1
        do n=2,len_trim(title(K))    ! skip control character
          if (title(K)(n:n).ne.' ') then
            n1=n1+1
            newtit(n1:n1)=title(K)(n:n)
          end if
        end do
        title(K)=newtit
      END DO

      WRITE(iu_j,*) "Zonal Budgets for surface type ",TERRAIN
      WRITE(iu_j,*) "Latitude"
      WRITE(iu_j,*) "Zonal Average"
      WRITE(iu_j,'(A4,100A)') "Lat",(TRIM(TITLE(K)(1:14)),K=1,KMAX)

      DO J=1,JM
        WRITE(iu_j,'(I4,100(1X,F8.3))') NINT(LAT_DG(J,1)),
     *       (BUDG(J,K),K=1,KMAX)
      END DO
      WRITE(iu_j,*)
C**** output hemispheric and global means
      WRITE(iu_j,'(A4,100F8.3)') "NH",(BUDG(JM+1,K),K=1,KMAX)
      WRITE(iu_j,'(A4,100F8.3)') "SH",(BUDG(JM+2,K),K=1,KMAX)
      WRITE(iu_j,'(A4,100F8.3)') "GLOB",(BUDG(JM+3,K),K=1,KMAX)
      WRITE(iu_j,*)

      return
      end

      subroutine open_ijk(filename,im_gcm,jm_gcm,lm_gcm)
!@sum  OPEN_IJK opens the lat-lon-height binary output file
!@auth M. Kelley
!@ver  1.0
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
!@var IM_GCM,JM_GCM,LM_GCM dimensions for ij output
      INTEGER, INTENT(IN) :: im_gcm,jm_gcm,lm_gcm

      call openunit(filename,iu_ijk,.true.,.false.)

C**** set dimensions
      im=im_gcm
      jm=jm_gcm
      lm=lm_gcm

      return
      end subroutine open_ijk

      subroutine close_ijk
!@sum  CLOSE_IJK closes the lat-lon-height binary output file
!@auth M. Kelley
!@ver  1.0
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
      call closeunit(iu_ijk)
      return
      end subroutine close_ijk

      subroutine POUT_IJK(TITLE,SNAME,LNAME,UNITS,XIJK,XJK,XK,IJGRID)
!@sum  POUT_IJK outputs lat-lon-height binary records
!@auth M. Kelley
!@ver  1.0
      USE GISSOUT
      IMPLICIT NONE
!@var TITLE 80 byte title including description and averaging period
      CHARACTER, DIMENSION(LM), INTENT(IN) :: TITLE*80
!@var SNAME short name of field
      CHARACTER, INTENT(IN) :: SNAME*30
!@var LNAME long name of field
      CHARACTER, INTENT(IN) :: LNAME*50
!@var UNITS units of field
      CHARACTER, INTENT(IN) :: UNITS*50
!@var XIJK lat/lon/height output field
      REAL*8, DIMENSION(IM,JM,LM), INTENT(IN) :: XIJK
!@var XJK lat sum/mean of output field
      REAL*8, DIMENSION(JM,LM), INTENT(IN) :: XJK
!@var XK global sum/mean of output field
      REAL*8, DIMENSION(LM), INTENT(IN) :: XK
!@var IJGRID = 1 for primary lat-lon grid, 2 for secondary lat-lon grid
      INTEGER, INTENT(IN) :: IJGRID
      INTEGER :: I,K

      DO K=1,LM
         WRITE(iu_ijk) TITLE(K), REAL(XIJK(:,:,K),KIND=4),
     &     REAL(XJK(:,K),KIND=4), REAL(XK(K),KIND=4)
      ENDDO
      return
      end subroutine POUT_IJK
