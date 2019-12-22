      module ncinp
      implicit none

      include '/usr/local/netcdf-3.4/include/netcdf.inc'

      private
      public ::
     &     accfile
     &    ,open_acc,get_dim_size,close_acc
     &    ,status,acc_fid
     &    ,xlabel
     &    ,ia,idacc,im,jm,lm
     &    ,psf,ptop,dtsrc,sday,grav,omega,rgas,kapa,sha,tf,ls1
     &    ,radius,dlat,dlon
     &    ,twopi

      character(len=80) :: accfile
      integer :: status,acc_fid

! GCM variables and parameters
      character(len=132) :: xlabel
      integer, dimension(12) :: idacc ! 12 hard-coded for now
      integer :: ia ! idacc-number of current variable

      real :: psf,ptop,dtsrc,sday,grav,omega,rgas,kapa,sha,tf
      real :: radius,dlat,dlon

      integer :: im,jm,lm
      integer :: ls1

! some other parameters
      real, parameter :: twopi=6.2831853

      contains

      subroutine open_acc
      logical :: ex
      character(len=80) :: err_str
      integer :: tmp_id
!
      inquire(file=accfile,exist=ex)
      if(.not.ex) then
         err_str = 'nonexistent file: '//trim(accfile)
         stop trim(err_str)
      endif
      status = nf_open(trim(accfile),nf_nowrite,acc_fid)
      if(status.ne.nf_noerr) then ! is it a netcdf file?
         write(6,*) 'not a netcdf file: '//trim(accfile)
         stop
      endif
! check whether it really is an acc-file
      status = nf_get_att_text(acc_fid,nf_global,'XLABEL',xlabel)
      if(status.ne.nf_noerr) then
         write(6,*) 'missing XLABEL, not an acc file: '//trim(accfile)
         stop
      endif
! preliminaries:

! get IDACC array
      status = nf_inq_varid(acc_fid,'IDACC',tmp_id)
      status = nf_get_var_int(acc_fid,tmp_id,idacc)

! get real parameters
      status = nf_get_att_real(acc_fid,nf_global,'PSF',psf)
      status = nf_get_att_real(acc_fid,nf_global,'PTOP',ptop)
      status = nf_get_att_real(acc_fid,nf_global,'DTsrc',dtsrc)
      status = nf_get_att_real(acc_fid,nf_global,'SDAY',sday)
      status = nf_get_att_real(acc_fid,nf_global,'GRAV',grav)
      status = nf_get_att_real(acc_fid,nf_global,'OMEGA',omega)
      status = nf_get_att_real(acc_fid,nf_global,'RGAS',rgas)
      status = nf_get_att_real(acc_fid,nf_global,'KAPA',kapa)
      status = nf_get_att_real(acc_fid,nf_global,'TF',tf)

! get im,jm,lm dimensions
      status = nf_get_att_int(acc_fid,nf_global,'IM0',im)
      status = nf_get_att_int(acc_fid,nf_global,'JM0',jm)
      status = nf_get_att_int(acc_fid,nf_global,'LM0',lm)
! get ls1
      status = nf_get_att_int(acc_fid,nf_global,'LS1',ls1)

! initialize certain constants
      sha=rgas/kapa
      dlat=twopi/float(jm-1) ! temporary
      dlon=twopi/float(im) ! temporary
      radius=6375000. ! temporary

      return
      end subroutine open_acc

      subroutine get_dim_size(dim_name,dim_size)
      character(len=20) :: dim_name
      integer :: dim_size

      integer :: tmp_id
      status = nf_inq_dimid(acc_fid,trim(dim_name),tmp_id)
      status = nf_inq_dimlen(acc_fid,tmp_id,dim_size)
      return
      end subroutine get_dim_size

      subroutine close_acc
      status = nf_close(acc_fid)
      return
      end subroutine close_acc

      end module ncinp

      subroutine getacc(acc_name,acc)
      use ncinp
      implicit none
      include '/usr/local/netcdf-3.4/include/netcdf.inc'
      integer :: varid
      character(len=20) :: acc_name
      real :: acc
      status = nf_inq_varid(acc_fid,trim(acc_name),varid)
      status = nf_get_var_real(acc_fid,varid,acc)
      status = nf_get_att_int(acc_fid,varid,'idacc',ia)
      return
      end subroutine getacc

      subroutine gettxt(acc_name,text)
      use ncinp
      implicit none
      include '/usr/local/netcdf-3.4/include/netcdf.inc'
      integer :: varid
      character(len=20) :: acc_name
      character :: text
      status = nf_inq_varid(acc_fid,trim(acc_name),varid)
      status = nf_get_var_text(acc_fid,varid,text)
      return
      end subroutine gettxt
