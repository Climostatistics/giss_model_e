#include "rundeck_opts.h"

      MODULE LAKES_COM
!@sum  LAKES_COM model variables for Lake/Rivers module
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : IM,JM,ioread,iowrite,lhead,irerun,irsfic
     *     ,irsficno
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm
#endif
      IMPLICIT NONE
      SAVE
!@var MWL mass of lake water (kg)
      REAL*8, DIMENSION(IM,JM) :: MWL
!@var GML total enthalpy of lake (J)
      REAL*8, DIMENSION(IM,JM) :: GML
!@var TLAKE temperature of lake (C)
      REAL*8, DIMENSION(IM,JM) :: TLAKE
!@var MLDLK mixed layer depth in lake (m)
      REAL*8, DIMENSION(IM,JM) :: MLDLK
!@var FLAKE variable lake fraction (1)
      REAL*8, DIMENSION(IM,JM) :: FLAKE
!@var TANLK tan(alpha) = slope for conical lake (1)
      REAL*8, DIMENSION(IM,JM) :: TANLK

#ifdef TRACERS_WATER
!@var TRLAKE tracer amount in each lake level (kg)      
      REAL*8, DIMENSION(NTM,2,IM,JM) :: TRLAKE
#endif

      END MODULE LAKES_COM

      SUBROUTINE io_lakes(kunit,iaction,ioerr)
!@sum  io_lakes reads and writes lake arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE LAKES_COM
      IMPLICIT NONE
      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "LAKE01"
#ifdef TRACERS_WATER
!@var TRHEADER Character string label for individual records
      CHARACTER*80 :: TRHEADER, TRMODULE_HEADER = "TRLAK01"

      write (TRMODULE_HEADER(lhead+1:80)
     *     ,'(a7,i3,a)')'R8 dim(',NTM,',2,im,jm):TRLAKE'
#endif

      MODULE_HEADER(lhead+1:80) = 'R8 dim(im,jm):MixLD,MWtr,Tlk,Enth'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,MLDLK,MWL,TLAKE,GML !,FLAKE
#ifdef TRACERS_WATER
        WRITE (kunit,err=10) TRMODULE_HEADER,TRLAKE
#endif
      CASE (IOREAD:)            ! input from restart file
        READ (kunit,err=10) HEADER,MLDLK,MWL,TLAKE,GML   !,FLAKE
        IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
          PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
          GO TO 10
        END IF
#ifdef TRACERS_WATER
        SELECT CASE (IACTION)
        CASE (IRERUN,IOREAD,IRSFIC,IRSFICNO)    ! reruns/restarts
          READ (kunit,err=10) TRHEADER,TRLAKE
          IF (TRHEADER(1:LHEAD).NE.TRMODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",TRHEADER
     *           ,TRMODULE_HEADER
            GO TO 10
          END IF
        END SELECT
#endif
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_lakes


