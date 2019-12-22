      MODULE SOMTQ_COM
!@sum  SOMTQ_COM contains the arrays containing second order moments
!@auth Gary Russell
!@ver  1.0
      USE QUSDEF
      USE MODEL_COM, only : im,jm,lm
      IMPLICIT NONE
      SAVE
      REAL*8, DIMENSION(NMOM,IM,JM,LM) :: TMOM,QMOM

      END MODULE SOMTQ_COM

      SUBROUTINE io_somtq(kunit,iaction,ioerr)
!@sum  io_somtq reads and writes second order moments to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,lhead
      USE SOMTQ_COM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "QUS01"

      write (MODULE_HEADER(lhead+1:80),'(a7,i2,a)')
     * 'R8 dim(',nmom,',im,jm,lm):Tmom,Qmom'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)           ! output to standard restart file
        WRITE (KUNIT,ERR=10) MODULE_HEADER,TMOM,QMOM
      CASE (IOREAD:)            ! input from restart file
        READ (KUNIT,ERR=10) HEADER,TMOM,QMOM
        IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
          PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_somtq

      subroutine tq_zmom_init(t,q)
      USE MODEL_COM, only : im,jm,lm,sige,sig
      USE SOMTQ_COM
      implicit none
      REAL*8, dimension(im,jm,lm) :: t,q
      integer :: i,j,l
      REAL*8 :: rdsig
C**** INITIALIZES VERTICAL SLOPES OF T,Q
      DO J=1,JM
      DO I=1,IM
         RDSIG=(SIG(1)-SIGE(2))/(SIG(1)-SIG(2))
         TMOM(MZ,I,J,1)=(T(I,J,2)-T(I,J,1))*RDSIG
         QMOM(MZ,I,J,1)=(Q(I,J,2)-Q(I,J,1))*RDSIG
         IF(Q(I,J,1)+QMOM(MZ,I,J,1).LT.0.) QMOM(MZ,I,J,1)=-Q(I,J,1)
         DO L=2,LM-1
            RDSIG=(SIG(L)-SIGE(L+1))/(SIG(L-1)-SIG(L+1))
            TMOM(MZ,I,J,L)=(T(I,J,L+1)-T(I,J,L-1))*RDSIG
            QMOM(MZ,I,J,L)=(Q(I,J,L+1)-Q(I,J,L-1))*RDSIG
            IF(Q(I,J,L)+QMOM(MZ,I,J,L).LT.0.) QMOM(MZ,I,J,L)=-Q(I,J,L)
         END DO
         RDSIG=(SIG(LM)-SIGE(LM+1))/(SIG(LM-1)-SIG(LM))
         TMOM(MZ,I,J,LM)=(T(I,J,LM)-T(I,J,LM-1))*RDSIG
         QMOM(MZ,I,J,LM)=(Q(I,J,LM)-Q(I,J,LM-1))*RDSIG
         IF(Q(I,J,LM)+QMOM(MZ,I,J,LM).LT.0.) QMOM(MZ,I,J,LM)=-Q(I,J,LM)
      END DO
      END DO
      return
      end subroutine tq_zmom_init
