c f90 -c -64 NCACC.f
      subroutine write_nc_acc(fileout)
      use MODEL_COM
      use GEOM
      use DAGCOM
      implicit none
      include 'netcdf.inc'

      integer :: varid,status,ncid
      character(len=80) :: fileout
      integer, dimension(4) :: dimids

      integer ::
     &     lon_did,lat_did,lonb_did,latb_did,sig_did,sige_did,
     &     str_levid,plm_did,ple_did,ntype_did

      character(len=20) :: var_name
      REAL*8, dimension(jm,ntype) :: aj_tmp

      integer :: i,j,k

c 1st elem of ple =  glob mean edge pres between layers 1 and 2
c 1st elem of ple_shift = glob mean surf pres
      REAL*8, dimension(lm+1) :: ple_shift

      ple_shift(1) = psf
      ple_shift(2:lm+1) = ple(:)

c refresh some global parameters that change over the course of a run
!!! should the following lines will be really removed?
!!! there are no dparm_defs, iparm_defs in DEFACC.f any more.
!!! Those who familiar with this code please check. I.A.
!!!      call iparm_defs
!!!      call dparm_defs

      status = nf_create (trim(fileout), nf_clobber, ncid)
c-----------------------------------------------------------------------
c write global header
c-----------------------------------------------------------------------
      status = nf_put_att_text(ncid,nf_global,'XLABEL',132,xlabel)
      status = nf_put_att_text(ncid,nf_global,'AMON',  4,amon)
      status = nf_put_att_text(ncid,nf_global,'AMON0',  4,amon0)
      call ncwrt_iparm(iparm_name,ncid,iparm,100)
      call ncwrt_dparm(dparm_name,ncid,dparm,100)
c-----------------------------------------------------------------------
c define some popular dimensions
c-----------------------------------------------------------------------
      status = nf_def_dim(ncid, 'longitude', im,          lon_did)
      status = nf_def_dim(ncid, 'latitude' ,jm,          lat_did)
      status = nf_def_dim(ncid, 'lonb', im,          lonb_did)
      status = nf_def_dim(ncid, 'latb' ,jm,        latb_did)
      status = nf_def_dim(ncid, 'sig',   lm,      sig_did)
      status = nf_def_dim(ncid, 'sige', lm+1,      sige_did)
      status = nf_def_dim(ncid, 'p',   lm,      plm_did)
      status = nf_def_dim(ncid, 'ple',   lm+1,      ple_did)
      status = nf_def_dim(ncid, 'NTYPE',ntype, ntype_did)
c-----------------------------------------------------------------------
c define some supplemental variables not currently present in acc-file
c-----------------------------------------------------------------------
c longitudes at a-grid grid box centers
      status = nf_def_var (ncid, 'longitude', nf_real, 1, lon_did,varid)
c latitudes at a-grid grid box centers
      status = nf_def_var (ncid, 'latitude', nf_real, 1, lat_did, varid)
c a-grid grid box areas
      status = nf_def_var (ncid, 'area', nf_real, 1, lat_did, varid)
c longitudes at b-grid grid box centers
      status = nf_def_var (ncid, 'lonb', nf_real, 1, lonb_did, varid)
c latitudes at b-grid grid box centers
      status = nf_def_var (ncid, 'latb', nf_real, 1, latb_did, varid)
c b-grid grid box areas
      status = nf_def_var (ncid, 'areab', nf_real, 1, latb_did, varid)
c b-grid grid box east-west length
      status = nf_def_var (ncid, 'dxv', nf_real, 1, latb_did, varid)
c sigmas
      status = nf_def_var(ncid, 'sig', nf_real, 1, sig_did, varid)
c edge sigmas
      status = nf_def_var(ncid, 'sige', nf_real, 1, sige_did, varid)
c global mean layer midpoint pressures
      status = nf_def_var(ncid, 'p', nf_real, 1, plm_did, varid)
c global mean layer edge pressures
      status = nf_def_var(ncid, 'ple', nf_real, 1, ple_did, varid)

c-----------------------------------------------------------------------
c idacc
c-----------------------------------------------------------------------
      status = nf_def_dim(ncid, 'NDACC', 12,  dimids(1))
      status = nf_def_var(ncid, 'IDACC', nf_int, 1, dimids, varid)
c-----------------------------------------------------------------------
c surface type names
c-----------------------------------------------------------------------
      status=nf_def_dim(ncid,'STYPE_CLEN',len(stype_names(1)),dimids(1))
      dimids(2) = ntype_did
      status = nf_def_var(ncid, 'STYPE_NAMES', nf_char, 2,dimids,varid)
c-----------------------------------------------------------------------
c keynr
c-----------------------------------------------------------------------
      status = nf_def_dim(ncid, 'NKEYNR',42, dimids(1))
      status = nf_def_dim(ncid, 'KEYNR2',50, dimids(2))
      status = nf_def_var(ncid, 'KEYNR', nf_int, 2,dimids,varid)
c-----------------------------------------------------------------------
c tsfrez
c-----------------------------------------------------------------------
      dimids(1) = lon_did
      dimids(2) = lat_did
      call ncdefarr(name_tsf,lname_tsf,units_tsf,
     &     ncid,ktsf,nf_real,2,dimids)
c-----------------------------------------------------------------------
c aj
c-----------------------------------------------------------------------
      dimids(1)=lat_did
      dimids(2)=ntype_did
      call ncdefarr(name_j,lname_j,units_j,
     &     ncid,kaj,nf_real,2,dimids)
      call ncput_idacc(name_j,ia_j,ncid,kaj)
      status = nf_def_dim(ncid, 'reg',nreg, dimids(1))
      call ncdefarr(name_reg,lname_j,units_j,
     &     ncid,kaj,nf_real,1,dimids)
      call ncput_idacc(name_reg,ia_j,ncid,kaj)

c-----------------------------------------------------------------------
c apj
c-----------------------------------------------------------------------
      dimids(1)=lat_did
      call ncdefarr(name_pj,lname_pj,units_pj,
     &     ncid,kapj,nf_real,1,dimids)
c-----------------------------------------------------------------------
c ajl
c-----------------------------------------------------------------------
      dimids(1)=lat_did
      dimids(2)=sig_did
      call ncdefarr(sname_jl,lname_jl,units_jl,
     &     ncid,kajlx,nf_real,2,dimids)
c-----------------------------------------------------------------------
c asjl
c-----------------------------------------------------------------------
      status = nf_def_dim(ncid, 'str_level', lm_req,  str_levid)
      dimids(1)=lat_did
      dimids(2)=str_levid
      call ncdefarr(name_sjl,lname_sjl,units_sjl,
     &     ncid,kasjl,nf_real,2,dimids)

c-----------------------------------------------------------------------
c aij
c-----------------------------------------------------------------------
      dimids(1)=lon_did
      dimids(2)=lat_did
      call ncdefarr(name_ij,lname_ij,units_ij,
     &     ncid,kaij,nf_real,2,dimids)
      call ncput_idacc(name_ij,ia_ij,ncid,kaij)
c-----------------------------------------------------------------------
c ail
c-----------------------------------------------------------------------
      dimids(1)=lon_did
      dimids(2)=sig_did
      call ncdefarr(name_il,lname_il,units_il,
     &     ncid,kail,nf_real,2,dimids)
c-----------------------------------------------------------------------
c energy
c-----------------------------------------------------------------------
      status = nf_def_dim(ncid, 'ENERGY1', 20, dimids(1))
      status = nf_def_dim(ncid, 'ENERGY2',100, dimids(2))
      status = nf_def_var(ncid, 'ENERGY' ,nf_real,2,dimids,varid)
c-----------------------------------------------------------------------
c consrv
c-----------------------------------------------------------------------
      dimids(1)=lat_did
      call ncdefarr(name_consrv,lname_consrv,units_consrv,
     &     ncid,kcon,nf_real,1,dimids)
c-----------------------------------------------------------------------
c speca
c-----------------------------------------------------------------------
      status = nf_def_dim(ncid, 'IM_BY2_PLUS1',im/2+1, dimids(1))
      status = nf_def_dim(ncid, 'SP_QUANT'    ,    20, dimids(2))
      status = nf_def_dim(ncid, 'SP_NSPHER'   ,     8, dimids(3))
      status = nf_def_var(ncid, 'SPECA'       ,nf_real,3,dimids,varid)
c-----------------------------------------------------------------------
c atpe
c-----------------------------------------------------------------------
      status = nf_def_dim(ncid, 'ATPE1', 8, dimids(1)) ! m=1->8 in d5
      status = nf_def_dim(ncid, 'NHEMI', 2, dimids(2))
      status = nf_def_var(ncid, 'ATPE', nf_real, 2, dimids, varid)
c-----------------------------------------------------------------------
c adiurn
c-----------------------------------------------------------------------
      status = nf_def_dim(ncid, 'HR_IN_DAY', hr_in_day, dimids(1))
      status = nf_def_dim(ncid, 'NDIUVAR',   ndiuvar, dimids(2))
      status = nf_def_dim(ncid, 'NDIUPT', ndiupt, dimids(3))
      status = nf_def_var(ncid, 'ADIURN',nf_real,3,dimids,varid)
      status = nf_put_att_text(ncid,varid,'NAMDD',4*ndiupt,namdd)
      status = nf_put_att_int(ncid,varid,'IJDD',nf_int,2*ndiupt,IJDD)
c-----------------------------------------------------------------------
c wave power
c-----------------------------------------------------------------------
      status = nf_def_dim(ncid, 'RE_AND_IM', RE_AND_IM, dimids(1))
      status = nf_def_dim(ncid, 'Max12hr_sequ', Max12hr_sequ, dimids(2))
      status = nf_def_dim(ncid, 'NWAV_DAG', NWAV_DAG, dimids(3))
      call ncdefarr(name_wave,lname_wave,units_wave,
     &     ncid,kwp,nf_real,3,dimids)
c-----------------------------------------------------------------------
c ajk
c-----------------------------------------------------------------------
      dimids(1)=latb_did
      dimids(2)=sig_did
      call ncdefarr(sname_jk,lname_jk,units_jk,
     &     ncid,kajkx,nf_real,2,dimids)
c-----------------------------------------------------------------------
c aijk
c-----------------------------------------------------------------------
      dimids(1)=lonb_did
      dimids(2)=latb_did
      dimids(3)=sig_did
      call ncdefarr(name_ijk,lname_ijk,units_ijk,
     &     ncid,kaijk,nf_real,3,dimids)
c-----------------------------------------------------------------------
c end netcdf definitions
c-----------------------------------------------------------------------

      status = nf_enddef (ncid)

c-----------------------------------------------------------------------
c write data to netcdf file
c-----------------------------------------------------------------------
      var_name='longitude'; call ncwrtdbl1(var_name,ncid,lon_dg(1,1))
      var_name='latitude';    call ncwrtdbl1(var_name,ncid,lat_dg(1,1))
      var_name='area';   call ncwrtdbl1(var_name,ncid,dxyp)
      var_name='lonb';   call ncwrtdbl1(var_name,ncid,lon_dg(1,2))
      var_name='latb';   call ncwrtdbl1(var_name,ncid,lat_dg(1,2))
      var_name='areab';  call ncwrtdbl1(var_name,ncid,dxyv)
      var_name='dxv';  call ncwrtdbl1(var_name,ncid,dxv)
      var_name='sig';    call ncwrtdbl1(var_name,ncid,sig)
      var_name='sige';   call ncwrtdbl1(var_name,ncid,sige)
      var_name='p';   call ncwrtdbl1(var_name,ncid,plm)
      var_name='ple';   call ncwrtdbl1(var_name,ncid,ple_shift)
      var_name='IDACC';  call ncwrtint1(var_name,ncid,idacc)
      var_name='STYPE_NAMES';  call ncwrtchar(var_name,ncid,stype_names)
      var_name='KEYNR';  call ncwrtint1(var_name,ncid,keynr)
      call ncwrtdbl(name_tsf,ncid,im*jm,ktsf,tsfrez)
      do k=1,kaj
         aj_tmp=aj(:,k,:)
         call ncwrtdbl1(name_j(k),ncid,aj_tmp)
      enddo
      call ncwrtdbl(name_reg,ncid,nreg,kaj,areg)
      call ncwrtdbl(name_pj,ncid,jm,kapj,apj)
      call ncwrtdbl(sname_jl,ncid,jm*lm,kajlx,ajl)
      call ncwrtdbl(name_sjl,ncid,jm*lm_req,kasjl,asjl)
      call ncwrtdbl(name_ij,ncid,im*jm,kaij,aij)
      call ncwrtdbl(name_il,ncid,im*lm,kail,ail)
      var_name='ENERGY'; call ncwrtdbl1(var_name,ncid,energy)
      call ncwrtdbl(name_consrv,ncid,jm,kcon,consrv)
      var_name='SPECA';  call ncwrtdbl1(var_name,ncid,speca)
      var_name='ATPE';   call ncwrtdbl1(var_name,ncid,atpe)
      var_name='ADIURN'; call ncwrtdbl1(var_name,ncid,adiurn)
      call ncwrtdbl(name_wave,ncid,re_and_im*Max12hr_sequ*NWAV_DAG,
     &     kwp,wave)
      call ncwrtdbl(sname_jk,ncid,jm*lm,kajkx,ajk)
      call ncwrtdbl(name_ijk,ncid,im*jm*lm,kaijk,aijk)
c-----------------------------------------------------------------------
c close netcdf file
c-----------------------------------------------------------------------
      status = nf_close (ncid)
      return
      end subroutine write_nc_acc

      subroutine ncdefarr(arr_name,arr_lname,arr_units,
     &     ncid,arr_size,nf_type,ndims,dimids)
      implicit none
      include 'netcdf.inc'
      integer :: ncid,arr_size,nf_type,ndims
      integer, dimension(ndims) :: dimids
      integer :: status
      character(len=20), dimension(arr_size) :: arr_name,arr_units
      character(len=80), dimension(arr_size) :: arr_lname
      integer :: n,varid
      do n=1,arr_size
         status = nf_def_var(
     &        ncid,trim(arr_name(n)),nf_type,ndims,dimids,varid)
         status = nf_put_att_text(ncid,varid,'long_name',
     &        len_trim(arr_lname(n)),trim(arr_lname(n)))
         status = nf_put_att_text(ncid,varid,'units',
     &        len_trim(arr_units(n)),trim(arr_units(n)))
      enddo
      return
      end subroutine ncdefarr

      subroutine ncput_idacc(arr_name,arr_idacc,ncid,arr_size)
      implicit none
      include 'netcdf.inc'
      integer :: ncid,arr_size
      character(len=20), dimension(arr_size) :: arr_name
      integer, dimension(arr_size) :: arr_idacc
      integer :: status,n,varid
      do n=1,arr_size
         status = nf_inq_varid(ncid,trim(arr_name(n)),varid)
         status = nf_put_att_int(ncid,varid,'idacc',nf_int,1,
     &        arr_idacc(n))
      enddo
      return
      end subroutine ncput_idacc

      subroutine ncwrtdbl(arr_name,ncid,inner_size,arr_size,arr)
      implicit none
      include 'netcdf.inc'
      integer :: ncid,inner_size,arr_size
      integer :: status
      character(len=20), dimension(arr_size) :: arr_name
      REAL*8, dimension(inner_size,arr_size) :: arr
      integer :: n,varid
      do n=1,arr_size
         status = nf_inq_varid(ncid,trim(arr_name(n)),varid)
         status = nf_put_var_double(ncid,varid,arr(1,n))
      enddo
      return
      end subroutine ncwrtdbl

      subroutine ncgetdbl(arr_name,ncid,inner_size,arr_size,arr)
      implicit none
      include 'netcdf.inc'
      integer :: ncid,inner_size,arr_size
      integer :: status
      character(len=20), dimension(arr_size) :: arr_name
      REAL*8, dimension(inner_size,arr_size) :: arr
      integer :: n,varid
      do n=1,arr_size
         status = nf_inq_varid(ncid,trim(arr_name(n)),varid)
         status = nf_get_var_double(ncid,varid,arr(1,n))
      enddo
      return
      end subroutine ncgetdbl

      subroutine ncwrtdbl1(var_name,ncid,var)
      implicit none
      include 'netcdf.inc'
      integer :: ncid
      integer :: status,varid
      character(len=20) :: var_name
      REAL*8 :: var
      status = nf_inq_varid(ncid,trim(var_name),varid)
      status = nf_put_var_double(ncid,varid,var)
      return
      end subroutine ncwrtdbl1

      subroutine ncgetdbl1(var_name,ncid,var)
      implicit none
      include 'netcdf.inc'
      integer :: ncid
      integer :: status,varid
      character(len=20) :: var_name
      REAL*8 :: var
      status = nf_inq_varid(ncid,trim(var_name),varid)
      status = nf_get_var_double(ncid,varid,var)
      return
      end subroutine ncgetdbl1

      subroutine ncwrtint1(var_name,ncid,var)
      implicit none
      include 'netcdf.inc'
      integer :: ncid
      integer :: status,varid
      character(len=20) :: var_name
      integer :: var
      status = nf_inq_varid(ncid,trim(var_name),varid)
      status = nf_put_var_int(ncid,varid,var)
      return
      end subroutine ncwrtint1

      subroutine ncgetint1(var_name,ncid,var)
      implicit none
      include 'netcdf.inc'
      integer :: ncid
      integer :: status,varid
      character(len=20) :: var_name
      integer :: var
      status = nf_inq_varid(ncid,trim(var_name),varid)
      status = nf_get_var_int(ncid,varid,var)
      return
      end subroutine ncgetint1

      subroutine ncwrtchar(var_name,ncid,var)
      implicit none
      include 'netcdf.inc'
      integer :: ncid
      integer :: status,varid
      character(len=20) :: var_name
      character :: var
      status = nf_inq_varid(ncid,trim(var_name),varid)
      status = nf_put_var_text(ncid,varid,var)
      return
      end subroutine ncwrtchar

      subroutine ncgetchar(var_name,ncid,var)
      implicit none
      include 'netcdf.inc'
      integer :: ncid
      integer :: status,varid
      character(len=20) :: var_name
      character :: var
      status = nf_inq_varid(ncid,trim(var_name),varid)
      status = nf_get_var_text(ncid,varid,var)
      return
      end subroutine ncgetchar

      subroutine ncwrt_iparm(iparm_name,ncid,iparm,kjc)
      implicit none
      include 'netcdf.inc'
      character(len=20), dimension(kjc) :: iparm_name
      integer, dimension(kjc) :: iparm
      integer :: ncid,kjc
      integer :: k,status
      do k=1,kjc
         if(iparm_name(k).eq.'nowrite') cycle
         status = nf_put_att_int(ncid,nf_global,trim(iparm_name(k)),
     &                           nf_int,1,iparm(k))
      enddo
      return
      end subroutine ncwrt_iparm

      subroutine ncget_iparm(iparm_name,ncid,iparm,kjc)
      implicit none
      include 'netcdf.inc'
      character(len=20), dimension(kjc) :: iparm_name
      integer, dimension(kjc) :: iparm
      integer :: ncid,kjc
      integer :: k,status
      do k=1,kjc
         if(iparm_name(k).eq.'nowrite') cycle
         status = nf_get_att_int(
     &        ncid,nf_global,trim(iparm_name(k)),iparm(k))
      enddo
      return
      end subroutine ncget_iparm

      subroutine ncwrt_dparm(dparm_name,ncid,dparm,krc)
      implicit none
      include 'netcdf.inc'
      character(len=20), dimension(krc) :: dparm_name
      REAL*8, dimension(krc) :: dparm
      integer :: ncid,krc
      integer :: k,status
      do k=1,krc
         if(dparm_name(k).eq.'nowrite') cycle
         status = nf_put_att_double(ncid,nf_global,trim(dparm_name(k)),
     &                           nf_real,1,dparm(k))
      enddo
      return
      end subroutine ncwrt_dparm

      subroutine ncget_dparm(dparm_name,ncid,dparm,krc)
      implicit none
      include 'netcdf.inc'
      character(len=20), dimension(krc) :: dparm_name
      REAL*8, dimension(krc) :: dparm
      integer :: ncid,krc
      integer :: k,status
      do k=1,krc
         if(dparm_name(k).eq.'nowrite') cycle
         status=nf_get_att_double(
     &        ncid,nf_global,trim(dparm_name(k)),dparm(k))
      enddo
      return
      end subroutine ncget_dparm
