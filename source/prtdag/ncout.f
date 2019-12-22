      module ncout
      implicit none

      include '/usr/local/netcdf-3.4/include/netcdf.inc'

      private
      public ::
     &     outfile
     &    ,open_out,def_dim_out,set_dim_out,close_out
     &    ,status_out,varid_out,out_fid
     &    ,ndims_out,dimids_out,file_dimlens
     &    ,units,long_name,missing

      character(len=80) :: outfile
      integer :: status_out,out_fid
      integer :: ndims_out
      integer :: varid_out
      integer, dimension(7) :: dimids_out
      character(len=20) :: units=''
      character(len=80) :: long_name=''

      integer, parameter :: ndfmax=100 ! max dims in file
      integer :: ndims_file=0
      integer, dimension(ndfmax) :: file_dimids
      integer, dimension(ndfmax) :: file_dimlens
      character(len=20), dimension(ndfmax) :: file_dimnames

      real :: missing=-1.e30

      contains
      
      subroutine open_out
      status_out = nf_create (trim(outfile), nf_clobber, out_fid)
      if(status_out.ne.nf_noerr) then
         write(6,*) 'cannot create: '//trim(outfile)
         stop
      endif
      return
      end subroutine open_out

      subroutine def_dim_out(dim_name,dimlen)
      character(len=20) :: dim_name
      integer :: dimlen
      integer :: tmp_id
      if(ndims_file.eq.ndfmax)
     &     stop 'def_dim_out: too many dimensions in file'
      ndims_file = ndims_file + 1
      file_dimlens(ndims_file) = dimlen
      file_dimnames(ndims_file)=dim_name
      status_out = nf_def_dim(out_fid,trim(dim_name),dimlen,tmp_id)
      file_dimids(ndims_file) = tmp_id
      return
      end subroutine def_dim_out

      subroutine set_dim_out(dim_name,dim_num)
      character(len=20) :: dim_name
      integer :: dim_num
      integer :: idim
      if(ndims_file.eq.0) stop 'set_dim_out: no dims defined'
      if(dim_num.gt.ndims_out) stop 'set_dim_out: invalid dim #'
      idim=1
      do while(idim.le.ndims_file)
         if(dim_name.eq.file_dimnames(idim)) then
            dimids_out(dim_num)=file_dimids(idim)
            exit
         endif
         idim = idim + 1
      enddo
      if(idim.gt.ndims_file) stop 'set_dim_out: invalid dim_name'
      return
      end subroutine set_dim_out

      subroutine close_out
      status_out = nf_close(out_fid)
      return
      end subroutine close_out

      end module ncout

      subroutine wrtarr(var_name,var)
      use ncout
      implicit none
      include '/usr/local/netcdf-3.4/include/netcdf.inc'
      character(len=20) :: var_name
      real :: var(1)
      integer :: var_nelems,n
      status_out = nf_redef(out_fid)
      status_out=nf_def_var(out_fid,trim(var_name),nf_real,
     &     ndims_out,dimids_out,varid_out)
      if(len_trim(units).gt.0) status_out =
     &     nf_put_att_text(out_fid,varid_out,
     &     'units',len_trim(units),units)
      if(len_trim(long_name).gt.0) status_out =
     &     nf_put_att_text(out_fid,varid_out,
     &     'long_name',len_trim(long_name),long_name)
c define missing value attribute only if there are missing values
      var_nelems = product(file_dimlens(dimids_out(1:ndims_out)))
      do n=1,var_nelems
         if(var(n).eq.missing) then
            status_out = nf_put_att_real(out_fid,varid_out,
     &           'missing_value',nf_float,1,missing)
            exit
         endif
      enddo
      status_out = nf_enddef(out_fid)
      status_out = nf_put_var_real(out_fid,varid_out,var)
      units=''
      long_name=''
      return
      end subroutine wrtarr
