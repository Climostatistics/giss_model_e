#include "rundeck_opts.h"

      SUBROUTINE io_rsf(kunit,it,iaction,ioerr)
!@sum   io_rsf controls the reading and writing of the restart files
!@auth  Gavin Schmidt
!@ver   1.0
!@calls io_model,io_ocean,io_lakes,io_seaice,io_earth,io_soils,io_snow
!@+     io_landice,io_bldat,io_pbl,io_clouds,io_somtq,io_rad,io_diags
!@+     io_ocdiag,io_icedyn,io_icdiag
      USE MODEL_COM, only : ioread_single,iowrite_single,Kradia

      IMPLICIT NONE
!@var iaction flag for reading or writing rsf file
!@var kunit Fortran unit number of file i/o
      INTEGER, INTENT(IN) :: iaction,kunit
!@var it hour of model run
      INTEGER, INTENT(INOUT) :: it
!@var IOERR (1,0,-1) if there (is, is maybe, is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var IT1 hour for correct reading check
!@var ITM maximum hour for post-processing 
      INTEGER IT1,itm

      ioerr=-1
      rewind kunit

C**** For all iaction < 0  ==> WRITE, For all iaction > 0  ==> READ
C**** Particular values may produce variations in indiv. i/o routines

C**** Calls to individual i/o routines
      call io_label  (kunit,it,itm,iaction,ioerr)
      it1=it
      if (Kradia.gt.0) then
        if (Kradia.gt.1 .and. iaction.ne.ioread_single .and.
     *   iaction.ne.iowrite_single) call io_rad (kunit,iaction,ioerr)
        call io_diags  (kunit,it,iaction,ioerr)
        if(it1.ne.it .or. ioerr.eq.1)
     &       call stop_model('restart problem',255)
        return
      end if
      if(iaction.ne.ioread_single.and.iaction.ne.iowrite_single) then
        call io_model  (kunit,iaction,ioerr)
        call io_strat  (kunit,iaction,ioerr)
        call io_ocean  (kunit,iaction,ioerr)
        call io_lakes  (kunit,iaction,ioerr)
        call io_seaice (kunit,iaction,ioerr)
        call io_earth  (kunit,iaction,ioerr)
        call io_soils  (kunit,iaction,ioerr)
        call io_snow   (kunit,iaction,ioerr)
        call io_landice(kunit,iaction,ioerr)
        call io_bldat  (kunit,iaction,ioerr)
        call io_pbl    (kunit,iaction,ioerr)
        call io_clouds (kunit,iaction,ioerr)
        call io_somtq  (kunit,iaction,ioerr)
        call io_rad    (kunit,iaction,ioerr)
        call io_icedyn (kunit,iaction,ioerr)
#ifdef TRACERS_ON
        call io_tracer (kunit,iaction,ioerr)
#endif
      end if
      call io_diags  (kunit,it,iaction,ioerr)
      call io_ocdiag (kunit,it,iaction,ioerr)
      call io_icdiag (kunit,it,iaction,ioerr)
#ifdef TRACERS_ON
      call io_trdiag (kunit,it,iaction,ioerr)
#endif

      if (it1.ne.it) THEN
        WRITE(6,*) "TIMES DO NOT MATCH READING IN RSF FILE",it,it1
        ioerr=1
      END IF
      if (ioerr.eq.1) WRITE(6,*) "I/O ERROR IN RESTART FILE: KUNIT="
     *     ,kunit
      close (kunit)

C**** return maximum time
      it=itm

      RETURN
      END SUBROUTINE io_rsf
