      MODULE CLOUDS_COM
!@sum  CLOUDS_COM model variables for moist convction and
!@+          large-scale condensation
!@auth M.S.Yao/T. Del Genio (modularisation by Gavin Schmidt)
!@ver  1.0 (taken from CB265)
      USE MODEL_COM, only : IM,JM,LM

      IMPLICIT NONE
      SAVE
!@var TTOLD,QTOLD previous potential temperature, humidity
      REAL*8, DIMENSION(LM,IM,JM) :: TTOLD,QTOLD
!@var SVLHX,SVLAT previous latent heat of evaporation
      REAL*8, DIMENSION(LM,IM,JM) :: SVLHX,SVLAT
!@var RHSAV previous relative humidity
      REAL*8, DIMENSION(LM,IM,JM) :: RHSAV
!@var CLDSAV previous cloud cover area (percent)
      REAL*8, DIMENSION(LM,IM,JM) :: CLDSAV
!@var PBLTOP height of PBL (m) (NOT USED)
      REAL*8, DIMENSION(IM,JM) :: PBLTOP

C**** variables saved for radiation calculations
!@var TAUSS optical depth from super-saturated clouds
      REAL*8, DIMENSION(LM,IM,JM) :: TAUSS
!@var TAUMC optical depth from moist-convective clouds
      REAL*8, DIMENSION(LM,IM,JM) :: TAUMC
!@var CLDSS super-saturated cloud cover area (percent)
      REAL*8, DIMENSION(LM,IM,JM) :: CLDSS
!@var CLDMC moist convective cloud cover area (percent)
      REAL*8, DIMENSION(LM,IM,JM) :: CLDMC
!@var CSIZMC,CSIZSS mc,ss effective cloud droplet radius (microns)
      REAL*8, DIMENSION(LM,IM,JM) :: CSIZMC,CSIZSS

!@var LLOW,LMID,LHI max levels for low, mid and high clouds
      INTEGER LLOW,LMID,LHI

      END MODULE CLOUDS_COM

      SUBROUTINE io_clouds(kunit,iaction,ioerr)
!@sum  io_clouds reads and writes cloud arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,lhead
      USE CLOUDS_COM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "CLD01"
!@var UDUM1,UDUM2,LDUM dummy variables
      REAL*8 UDUM1,UDUM2
      INTEGER LDUM1

      write(MODULE_HEADER(lhead+1:80),'(a)')
     *  'R8 dim(im,jm,lm):potT,Hum,LatHeat,RHum,CldCv (all old)'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)           ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,
     *     TTOLD,QTOLD,SVLHX,RHSAV,CLDSAV
      CASE (IOREAD:)            ! input from restart file
        READ (kunit,err=10) HEADER,
     *     TTOLD,QTOLD,SVLHX,RHSAV,CLDSAV
        IF (HEADER(1:15).NE.MODULE_HEADER(1:15)) THEN
          PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_clouds
