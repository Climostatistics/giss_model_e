      program ijkdag
C****
C****
      use ncinp
      use ncout
      implicit none
      integer :: nargs,iargc
      character(len=80) :: stop_str,cmd_str
      real, dimension(:), allocatable :: lon_gcm,lat_gcm,plm_gcm
      real, dimension(:), allocatable :: lona,lata,ple_gcm
      real, dimension(:), allocatable :: dxyp,dxv,dyp
      real, dimension(:,:,:), allocatable :: arrn,arrb,arrb2,dpb
      character(len=20) :: var_name,acc_name,dim_name
      integer :: i,j,l,im1,ip1
      real :: dpaccg
      call getarg(0,cmd_str)
      nargs = iargc()
      if(nargs.ne.2) then
         stop_str = 'usage: '//trim(cmd_str)//' acc_file out_file'
         stop trim(stop_str)
      endif
      call getarg(1,accfile)
      call getarg(2,outfile)
      if(accfile.eq.outfile) stop 'cannot overwrite input file'

! open the acc file
      call open_acc

! allocate space using these dimensions
      allocate(lon_gcm(im),lona(im))
      allocate(lat_gcm(jm),lata(jm))
      allocate(plm_gcm(lm),ple_gcm(lm+1))
      allocate(dxyp(jm),dxv(jm),dyp(jm))
      allocate(arrn(im,jm,lm))
      allocate(dpb(im,jm,lm))
      allocate(arrb(im,jm-1,lm),arrb2(im,jm-1,lm))
! get lon,lat,p coordinates
      acc_name='lonb'; call getacc(acc_name,lon_gcm)
      acc_name='latb'; call getacc(acc_name,lat_gcm)
      lat_gcm(1:jm-1) = lat_gcm(2:jm) ! shift off nonexistent 1st latitude
      acc_name='p'; call getacc(acc_name,plm_gcm)
      acc_name='longitude'; call getacc(acc_name,lona)
      acc_name='latitude'; call getacc(acc_name,lata)
      acc_name='ple'; call getacc(acc_name,ple_gcm)
! get geometry
      acc_name='area'; call getacc(acc_name,dxyp)
      acc_name='dxv'; call getacc(acc_name,dxv)
      acc_name='dyp'; call getacc(acc_name,dyp)
      dyp(:)=radius*dlat; dyp(1)=.5*dyp(1); dyp(jm)=.5*dyp(jm)
! define output file
      call open_out
      dim_name='longitude'; call def_dim_out(dim_name,im)
      dim_name='latitude'; call def_dim_out(dim_name,jm-1)
      dim_name='p'; call def_dim_out(dim_name,lm)
      dim_name='lona'; call def_dim_out(dim_name,im)
      dim_name='lata'; call def_dim_out(dim_name,jm)
      dim_name='ple'; call def_dim_out(dim_name,lm)

      ndims_out = 1

c b-grid
      dim_name='longitude'; call set_dim_out(dim_name,1)
      units='degrees_east'
      var_name='longitude';call wrtarr(var_name,lon_gcm)
      dim_name='latitude'; call set_dim_out(dim_name,1)
      units='degrees_north'
      var_name='latitude';call wrtarr(var_name,lat_gcm)
      dim_name='p'; call set_dim_out(dim_name,1)
      units='mb'; long_name='Pressure'
      var_name='p';call wrtarr(var_name,plm_gcm)

c a-grid (for omega)
      dim_name='lona'; call set_dim_out(dim_name,1)
      units='degrees_east'
      var_name='lona';call wrtarr(var_name,lona)
      dim_name='lata'; call set_dim_out(dim_name,1)
      units='degrees_north'
      var_name='lata';call wrtarr(var_name,lata)
      dim_name='ple'; call set_dim_out(dim_name,1)
      units='mb'; long_name='Pressure'
      var_name='ple';call wrtarr(var_name,ple_gcm(2))

      ndims_out = 3
      dim_name='longitude'; call set_dim_out(dim_name,1)
      dim_name='latitude'; call set_dim_out(dim_name,2)
      dim_name='p'; call set_dim_out(dim_name,3)

! get DPB array from acc file
      acc_name='DPB'; call getacc(acc_name,dpb)

! u-velocity
      acc_name='UDPB'; call getacc(acc_name,arrn)
      call scale(arrb,arrn,dpb,im,jm,lm,missing,1.)
      var_name='u'; call wrtarr(var_name,arrb)

! v-velocity
      acc_name='VDPB'; call getacc(acc_name,arrn)
      call scale(arrb,arrn,dpb,im,jm,lm,missing,1.)
      var_name='v'; call wrtarr(var_name,arrb)

! temperature
      acc_name='TDPBx4'; call getacc(acc_name,arrn)
      call scale(arrb,arrn,dpb,im,jm,lm,missing,.25)
      var_name='t'; call wrtarr(var_name,arrb)
! and height
      acc_name='DSEDPBx4'; call getacc(acc_name,arrn)
      where(dpb(:,2:jm,:).gt.0.)
         arrb(:,:,:) = (.25*arrn(:,2:jm,:)/dpb(:,2:jm,:) -
     &        sha*arrb(:,:,:))/grav
      elsewhere
         arrb(:,:,:) = missing
      end where
      var_name='z'; call wrtarr(var_name,arrb)

! humidity
      acc_name='QDPBx4'; call getacc(acc_name,arrn)
      call scale(arrb,arrn,dpb,im,jm,lm,missing,.25)
      var_name='q'; call wrtarr(var_name,arrb)


      ndims_out = 3
      dim_name='lona'; call set_dim_out(dim_name,1)
      dim_name='lata'; call set_dim_out(dim_name,2)
      dim_name='ple'; call set_dim_out(dim_name,3)

! omega-velocity
      acc_name='UDPB'; call getacc(acc_name,arrn)
      arrb(:,:,:) = arrn(:,2:jm,:)
      acc_name='VDPB'; call getacc(acc_name,arrn)
      arrb2(:,:,:)= arrn(:,2:jm,:)
c zero flux at the top
      arrn(:,:,lm) = 0.
c loop downward through layers
      do l=lm-1,1,-1
c**** calculate non-polar vertical winds
         do j=2,jm-1
            im1=im
            do i=1,im
               arrn(i,j,l) = arrn(i,j,l+1) +
     &              .5/(dxyp(j)) *
     &              ( dyp(j) * ( (arrb(im1,j-1,l+1) + arrb(im1,j,l+1))
     &              - (arrb(i,j-1,l+1) + arrb(i,j,l+1)))
     &              + dxv(j) * (arrb2(im1,j-1,l+1) + arrb2(i,j-1,l+1))
     &              - dxv(j+1) * (arrb2(im1,j,l+1) + arrb2(i,j,l+1)) )
               im1=i
            enddo
         enddo
c**** calculate polar vertical winds
         arrn(:, 1,l) = arrn(:, 1,l+1) +
     &        dxv( 2)*sum(arrb2(:,1   ,l+1))/(float(im)*dxyp(1))
         arrn(:,jm,l) = arrn(:,jm,l+1) -
     &        dxv(jm)*sum(arrb2(:,jm-1,l+1))/(float(im)*dxyp(jm))
      enddo
      arrn = arrn*(-1.e5/idacc(4))
      do l=1,lm
         dpaccg=idacc(4)*(ple_gcm(l)-ple_gcm(l+1))
         do j=2,jm
            i=im
            do ip1=1,im
               if(dpb(i,j,l).lt.dpaccg) then
                  arrn(i  ,j-1,l) = missing
                  arrn(i  ,j  ,l) = missing
                  arrn(ip1,j-1,l) = missing
                  arrn(ip1,j  ,l) = missing
               endif
               i=ip1
            enddo
         enddo
      enddo
      var_name='omega'; call wrtarr(var_name,arrn)

      call close_acc
      call close_out

      end program ijkdag

      subroutine scale(arrb,arrn,dpb,im,jm,lm,missing,fac)
      implicit none
      real, dimension(im,jm-1,lm) :: arrb
      real, dimension(im,jm,lm) :: arrn,dpb
      integer :: im,jm,lm
      real :: missing,fac
      where(dpb(:,2:jm,:).gt.0.)
         arrb(:,:,:) = fac*arrn(:,2:jm,:)/dpb(:,2:jm,:)
      elsewhere
         arrb(:,:,:) = missing
      end where
      return
      end subroutine scale
