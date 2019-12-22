#include "rundeck_opts.h"

      MODULE PBLCOM
!@sum  PBLCOM contains the arrays used by the Boundary Layer code
!@auth Greg Hartke/Ye Cheng
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm
#ifdef TRACERS_ON
      USE TRACER_COM, only : ntm
#endif
      USE SOCPBL, only : npbl=>n
      IMPLICIT NONE
      SAVE
!@var uabl boundary layer profile for zonal wind
!@var vabl boundary layer profile for meridional wind
!@var tabl boundary layer profile for temperature
!@var qabl boundary layer profile for humidity
      real*8, dimension(npbl,im,jm,4) :: uabl,vabl,tabl,qabl
!@var eabl boundary layer profile for turbulent KE (calc. on sec. grid)
      real*8, dimension(npbl,im,jm,4) :: eabl
#ifdef TRACERS_ON
!@var trabl boundary layer profile for tracers
      real*8, dimension(npbl,ntm,im,jm,4) :: trabl
#endif

!@var cmgs drag coefficient (dimensionless surface momentum flux)
!@var chgs Stanton number   (dimensionless surface heat flux)
!@var cqgs Dalton number    (dimensionless surface moisture flux)
      real*8, dimension(im,jm,4) :: cmgs,chgs,cqgs

!@var ipbl flag for whether pbl properties were found at last timestep
      integer, dimension(im,jm,4) :: ipbl

!@var ROUGHL log10(zgs/roughness length), prescribed with zgs=30 m.
      REAL*8, dimension(im,jm) :: roughl
!@var WSAVG     COMPOSITE SURFACE WIND MAGNITUDE (M/S)
!@var TSAVG     COMPOSITE SURFACE AIR TEMPERATURE (K)
!@var QSAVG     COMPOSITE SURFACE AIR SPECIFIC HUMIDITY (1)
!@var DCLEV     LAYER TO WHICH DRY CONVECTION MIXES (1)
!@var USAVG     COMPOSITE SURFACE U WIND
!@var VSAVG     COMPOSITE SURFACE V WIND
!@var TAUAVG    COMPOSITE SURFACE MOMENTUM TRANSFER (TAU)
!@var w2_l1     COMPOSITE vertical component of t.k.e. at gcm layer 1
!@var USTAR_pbl friction velocity (sqrt of srfc mom flux) (m/s)
      REAL*8, dimension(im,jm) ::
     &     wsavg,tsavg,qsavg,dclev,usavg,vsavg,tauavg,tgvavg,qgavg
     &    ,w2_l1
      REAL*8, dimension(im,jm,4) :: ustar_pbl

!@var egcm  3-d turbulent kinetic energy in the whole atmosphere
!@var w2gcm vertical component of egcm
!@var t2gcm 3-d turbulent temperature variance in the whole atmosphere
      real*8, dimension(lm,im,jm) :: egcm,w2gcm,t2gcm

!@var uflux surface turbulent u-flux (=-<uw>)
!@var vflux surface turbulent v-flux (=-<vw>)
!@var tflux surface turbulent t-flux (=-<tw>)
!@var qflux surface turbulent q-flux (=-<qw>)
      real*8, dimension(im,jm) :: uflux,vflux,tflux,qflux

      END MODULE PBLCOM

      SUBROUTINE io_pbl(kunit,iaction,ioerr)
!@sum  io_pbl reads and writes model variables to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,irsfic,irerun,iowrite,irsficno,lhead
      USE PBLCOM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "PBL01"
#ifdef TRACERS_ON
!@var TR_HEADER Character string label for tracer record
      CHARACTER*80 :: TR_HEADER, TR_MODULE_HEADER = "TRPBL01"
      write (TR_MODULE_HEADER(lhead+1:80),'(a7,i2,i3,a)') 'R8 dim(',npbl
     *     ,ntm,',ijm,4):TRt'
#endif
      write (MODULE_HEADER(lhead+1:80),'(a7,i2,a)') 'R8 dim(',npbl,
     *  ',ijm,4):Ut,Vt,Tt,Qt,Et dim(ijm,4,3):Cmhq, I:Ipb(ijm,4)'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        WRITE (KUNIT,ERR=10) MODULE_HEADER,UABL,VABL,TABL,QABL,EABL,
     *     CMGS,CHGS,CQGS,IPBL
#ifdef TRACERS_ON
        WRITE (KUNIT,ERR=10) TR_MODULE_HEADER,TRABL
#endif
      CASE (IOREAD:)            ! input from restart file or restart
        READ (KUNIT,ERR=10) HEADER,UABL,VABL,TABL,QABL,EABL,CMGS,CHGS,
     *     CQGS,IPBL
        IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
          PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
          GO TO 10
        END IF
#ifdef TRACERS_ON
        SELECT CASE (IACTION)
        CASE (IOREAD,IRERUN,IRSFIC,IRSFICNO)    ! restarts
          READ (KUNIT,ERR=10) TR_HEADER,TRABL
          IF (TR_HEADER(1:LHEAD).NE.TR_MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in tracer module version ",TR_HEADER
     *           ,TR_MODULE_HEADER
            GO TO 10
          END IF
        END SELECT
#endif
      END SELECT
      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_pbl

      SUBROUTINE io_bldat(kunit,iaction,ioerr)
!@sum  io_bldat reads and writes boundary layer data to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,lhead
      USE PBLCOM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "BLD02"

      MODULE_HEADER(lhead+1:80) = 'R8 dim(ijm):ws,ts,qs,'//
     *  'LvlDC,us,vs,tau,u*(.,4),ke;w2(lijm),tgv,qg'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,wsavg,tsavg,qsavg,dclev
     *       ,usavg,vsavg,tauavg,ustar_pbl,egcm,w2gcm,tgvavg,qgavg
      CASE (IOREAD:)            ! input from restart file
        READ (kunit,err=10) HEADER,wsavg,tsavg,qsavg,dclev,usavg
     *       ,vsavg,tauavg,ustar_pbl,egcm,w2gcm,tgvavg,qgavg
        IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
          PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_bldat

