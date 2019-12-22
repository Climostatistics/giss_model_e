      program jldag
C****
C****
      use ncinp
      use ncout
      implicit none
      integer :: lm_req
      integer :: nargs,iargc
      character(len=80) :: stop_str,cmd_str
      character(len=80) :: title
      real :: dtsec,dtday,bygrav
c from diagjl/k
      real, dimension(:,:,:), allocatable :: dpab,dsjk
      real, dimension(:,:), allocatable :: ax,axb,bx,arqx,brqx
      real, dimension(:), allocatable ::
     &     DXYP,BYDXYP,DXYV,DXV,BYP,BYPV,DXYPPO,ONES,ONESPO
      real, dimension(:), allocatable ::
     &     BYDPS,BYPKS,SIG,DSIG,BYDSIG,SIGE,PLE,PLM,PME

      INTEGER :: J,L,N,J1,K,KM,LR,M,N

      REAL ::
     &     BYIACN,BYIADA,BYIARD,BYIMDA,
     &     SCALE,XWON,BY100G,BYN,P1000K,SDDP,BYIM

      REAL, PARAMETER :: ONE=1.,P1000=1000.

      real, dimension(:), allocatable :: lats,latsb
      character(len=20) :: var_name,acc_name,dim_name

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

      dim_name='str_level'; call get_dim_size(dim_name,lm_req)

! allocate space using these dimensions
      allocate(lats(jm),latsb(jm))
      allocate(sig(lm),dsig(lm),sige(lm+1),ple(lm+1),plm(lm+lm_req))
      allocate(ax(jm,lm),bx(jm,lm),dpab(jm,lm,2),dsjk(jm,lm,2))
      allocate(axb(jm-1,lm))
      allocate(arqx(jm,lm+lm_req),brqx(jm,lm+lm_req))
      allocate(DXYP(jm),BYDXYP(jm),DXYV(jm),DXV(jm))
      allocate(BYP(jm),BYPV(jm),DXYPPO(jm),ONES(jm),ONESPO(jm))
      allocate(BYDPS(lm_req),BYPKS(lm_req),PME(LM),BYDSIG(LM))

! get lat,ht coordinates
      acc_name='latitude'; call getacc(acc_name,lats)
      acc_name='latb'; call getacc(acc_name,latsb)
      acc_name='sig'; call getacc(acc_name,sig)
      acc_name='sige'; call getacc(acc_name,sige)
      latsb(1:jm-1)=latsb(2:jm) ! shift gcm's bgrid latitudes
      dsig(1:lm)=sige(1:lm)-sige(2:lm+1)
      bydsig=1./dsig
! get gridbox areas
      acc_name='area'; call getacc(acc_name,dxyp)
      acc_name='areab'; call getacc(acc_name,dxyv)
      acc_name='dxv'; call getacc(acc_name,dxv)
      bydxyp=1./dxyp

! define output file
      call open_out
      dim_name='latitude'; call def_dim_out(dim_name,jm)
      dim_name='latb'; call def_dim_out(dim_name,jm-1)
      dim_name='p'; call def_dim_out(dim_name,lm)
      dim_name='sig'; call def_dim_out(dim_name,lm)
      dim_name='sige'; call def_dim_out(dim_name,lm)
      dim_name='prqt'; call def_dim_out(dim_name,lm+lm_req)
      dim_name='pe'; call def_dim_out(dim_name,lm)

c**** initialize certain quantities
      byim=1./dble(im)
      xwon=twopi/(dlon*dble(im)) ! for wonderland model?
      km=lm
      bygrav=1./grav
      by100g=.01*bygrav
      p1000k=p1000**kapa
      ple(:)=sige(:)*(psf-ptop)+ptop
      plm(1:lm)=sig(1:lm)*(psf-ptop)+ptop
      pme(:)=(psf-ptop)*sige(1:lm)+ptop
      plm(lm+1)=.75*ple(lm+1)
      plm(lm+2)=.35*ple(lm+1)
      plm(lm+3)=.1*ple(lm+1)
      bydps(1)=1./(.5*ple(lm+1))
      bydps(2)=1./(.3*ple(lm+1))
      bydps(3)=1./(.2*ple(lm+1))
      bypks(1)=1./(.75*ple(lm+1))**kapa
      bypks(2)=1./(.35*ple(lm+1))**kapa
      bypks(3)=1./(.1*ple(lm+1))**kapa
      do j=1,jm
         dxyppo(j)=dxyp(j)
         ones(j)=1.
         onespo(j)=1.
      enddo
      dxyppo(jm)=dxyp(jm)*dble(im)
      dxyppo(1)=dxyp(1)*dble(im)
      onespo(1)=dble(im)
      onespo(jm)=dble(im)

! put lat,ht into output file
      ndims_out = 1
      dim_name='latitude'; call set_dim_out(dim_name,1)
      units='degrees_north'
      var_name='latitude'; call wrtarr(var_name,lats)
      dim_name='latb'; call set_dim_out(dim_name,1)
      units='degrees_north'
      var_name='latb'; call wrtarr(var_name,latsb)
      dim_name='p'; call set_dim_out(dim_name,1)
      units='mb'
      var_name='p'; call wrtarr(var_name,plm)
      dim_name='sig'; call set_dim_out(dim_name,1)
      units='1'
      var_name='sig'; call wrtarr(var_name,sig)
      dim_name='sige'; call set_dim_out(dim_name,1)
      units='1'
      var_name='sige'; call wrtarr(var_name,sige(2))
      dim_name='prqt'; call set_dim_out(dim_name,1)
      units='mb'
      var_name='prqt'; call wrtarr(var_name,plm)
      dim_name='pe'; call set_dim_out(dim_name,1)
      units='mb'
      var_name='pe'; call wrtarr(var_name,ple)


! find total integration time in seconds, days, and other units
      dtsec = idacc(1)*3600.
      dtday = idacc(1)/24.

      BYIACN=1./(IDACC(1)+1.D-20)
      BYIARD=1./(IDACC(2)+1.D-20)
      BYIADA=1./(IDACC(4)+1.D-20)
      BYIMDA=BYIADA*BYIM
      acc_name='APJ1'; call getacc(acc_name,byp)
      DO J=1,JM
         BYP(J)=1./(BYP(J)+1.D-20)
      ENDDO

      ndims_out = 2

! (re)set shape of output arrays
      dim_name='latitude'; call set_dim_out(dim_name,1)
      dim_name='prqt'; call set_dim_out(dim_name,2)
c****
      title='SOLAR RADIATION HEATING RATE (DEGREES KELVIN/DAY)'
      acc_name='AJL09'; call getacc(acc_name,arqx(1,1))
      SCALE=1.D-2*GRAV*SDAY*IDACC(4)*BYIARD/SHA
      call scale2d(arqx(1,1),scale,byp,bydsig,jm,lm)
      acc_name='ASJL03'; call getacc(acc_name,arqx(1,lm+1))
      SCALE=1.D-2*GRAV*SDAY*BYIM*BYIARD/SHA
      call scale2d(arqx(1,lm+1),scale,onespo,bydps,jm,lm_req)
      var_name='srheat'; call wrtarr(var_name,arqx)
c****
      title='THERMAL RADIATION COOLING RATE (DEGREES KELVIN/DAY)'
      acc_name='AJL10'; call getacc(acc_name,brqx(1,1))
      SCALE=-1.D-2*GRAV*SDAY*IDACC(4)*BYIARD/SHA
      call scale2d(brqx(1,1),scale,byp,bydsig,jm,lm)
      acc_name='ASJL04'; call getacc(acc_name,brqx(1,lm+1))
      SCALE=-1.D-2*GRAV*SDAY*BYIM*BYIARD/SHA
      call scale2d(brqx(1,lm+1),scale,onespo,bydps,jm,lm_req)
      var_name='trcool'; call wrtarr(var_name,brqx)
c****
      title='TOTAL RADIATION COOLING RATE (DEGREES KELVIN/DAY)'
      ARQX(:,:)=-ARQX(:,:)+BRQX(:,:)
      var_name='rdcool'; call wrtarr(var_name,arqx)
c****
! (re)set shape of output arrays
      dim_name='latitude'; call set_dim_out(dim_name,1)
      dim_name='sig'; call set_dim_out(dim_name,2)

      title='TOTAL CLOUD COVER'
      acc_name='AJL19'; call getacc(acc_name,ax)
      SCALE=BYIARD*BYIM
      call scale2d(ax,scale,onespo,ones,jm,lm)
      var_name='pcloud'; call wrtarr(var_name,ax)
c****
      title='SUPER SATURATION CLOUD COVER'
      acc_name='AJL28'; call getacc(acc_name,ax)
      SCALE=BYIARD*BYIM
      call scale2d(ax,scale,onespo,ones,jm,lm)
      var_name='pcldss'; call wrtarr(var_name,ax)
c****
      title='MOIST CONVECTIVE CLOUD COVER'
      acc_name='AJL29'; call getacc(acc_name,ax)
      SCALE=BYIARD*BYIM
      call scale2d(ax,scale,onespo,ones,jm,lm)
      var_name='pcldmc'; call wrtarr(var_name,ax)
c****
      title='LARGE SCALE CONDENSATION HEATING (DEGREES KELVIN/DAY)'
      acc_name='AJL11'; call getacc(acc_name,ax)
      scale=IDACC(4)*BYIACN*SDAY/DTsrc
      call scale2d(ax,scale,byp,ones,jm,ls1-1)
      scale=byim*BYIACN*SDAY/(DTsrc*(psf-ptop))
      call scale2d(ax(1,ls1),scale,ones,ones,jm,lm-ls1+1)
      var_name='cdheat'; call wrtarr(var_name,ax)
c****
      title='HEATING BY DRY CONVECTION (DEGREES KELVIN/DAY)'
      acc_name='AJL12'; call getacc(acc_name,ax)
      scale=IDACC(4)*BYIACN*SDAY/DTsrc
      call scale2d(ax,scale,byp,ones,jm,ls1-1)
      scale=byim*BYIACN*SDAY/(DTsrc*(psf-ptop))
      call scale2d(ax(1,ls1),scale,ones,ones,jm,lm-ls1+1)
      var_name='dcheat'; call wrtarr(var_name,ax)
c****
      title='CHANGE OF LATENT HEAT BY DRY CONV. (WATTS/UNIT SIGMA/M2)'
      acc_name='AJL55'; call getacc(acc_name,ax)
      SCALE=(100./grav)*SHA*BYIACN*BYIM/DTsrc
      call scale2d(ax,scale,onespo,bydsig,jm,lm)
      var_name='dcdry'; call wrtarr(var_name,ax)
c****
      title='TOTAL HEATING BY MOIST CONVECTION (DEGREES KELVIN/DAY)'
      acc_name='AJL56'; call getacc(acc_name,ax)
      scale=IDACC(4)*BYIACN*SDAY/DTsrc
      call scale2d(ax,scale,byp,ones,jm,ls1-1)
      scale=byim*BYIACN*SDAY/(DTsrc*(psf-ptop))
      call scale2d(ax(1,ls1),scale,ones,ones,jm,lm-ls1+1)
      var_name='mcheat'; call wrtarr(var_name,ax)
c****
      title='TOTAL DRYING BY MOIST CONVECTION (DEG K/DAY EQUIVALENT)'
      acc_name='AJL57'; call getacc(acc_name,ax)
      scale=IDACC(4)*BYIACN*SDAY/DTsrc
      call scale2d(ax,scale,byp,ones,jm,ls1-1)
      scale=byim*BYIACN*SDAY/(DTsrc*(psf-ptop))
      call scale2d(ax(1,ls1),scale,ones,ones,jm,lm-ls1+1)
      var_name='mcdry'; call wrtarr(var_name,ax)
c****
! (re)set shape of output arrays
      dim_name='latitude'; call set_dim_out(dim_name,1)
      dim_name='sige'; call set_dim_out(dim_name,2)

      title='MOIST CONVECTION MASS FLUX (KG/SECOND/M2)'
      acc_name='AJL08'; call getacc(acc_name,ax)
      SCALE=(100./grav)*BYIACN*BYIM/DTsrc
      call scale2d(ax,scale,onespo,ones,jm,lm)
      var_name='mcmflx'; call wrtarr(var_name,ax)
c****
C**** BEGINNING OF DIAGJK
! note than ajk01 is on a-grid
      acc_name='AJK01'; call getacc(acc_name,dpab(1,1,1))
      acc_name='AJK02'; call getacc(acc_name,dpab(1,1,2))

      DO J=1,JM
      BYP(J)=1./(SUM(DPAB(J,1:KM,1))+1.D-20)
      BYPV(J)=1./(SUM(DPAB(J,1:KM,2))+1.D-20)
      ENDDO

cC**** WTJ: area weighting for JKMAP, JLMAP hemispheres
c      DO J=1,JM
c        WTJ(J,1,1)=1.
c        WTJ(J,2,1)=2.*FIM*DXYP(J)/AREAG
c      END DO
c      DO J=2,JM
c        WTJ(J,1,2)=1.
c        WTJ(J,2,2)=2.*FIM*DXYV(J)/AREAG
c      END DO
cc**** equatorial velocity point belongs in neither hemisphere
c      WTJ(JM/2+1,1,2)=.5
c      WTJ(JM/2+1,2,2)=WTJ(JM/2+1,2,2)/2.
C****
C**** INITIALIZE DELTA SIGMA IN PRESSURE COORDINATES
C****
      do j1=1,2
         do j=j1,jm
         DSJK(J,1:KM,J1)=DPAB(J,1:KM,J1)/(SUM(DPAB(J,1:KM,J1))+1.D-20)
         enddo
      enddo
      where(dsjk.le.0.) dsjk=missing
! (re)set shape of output arrays
      dim_name='latb'; call set_dim_out(dim_name,1)
      dim_name='sig'; call set_dim_out(dim_name,2)

C**** # OF GRIDPOINTS, DELTA P, S.D. OF DELTA P
      title='NUMBER OF GRIDPOINTS INCLUDED IN AVERAGE CP'
      acc_name='AJK24'; call getacc(acc_name,ax)
      SCALE=XWON*BYIADA
      call scale2d(ax,scale,ones,ones,jm,lm)
      var_name='ngdpts'; call wrtarr(var_name,ax)
c****
      title='PRESSURE DIFFERENCES  (MB) CP'
      SCALE=BYIMDA
      call scale2d(ax,scale,ones,ones,jm,lm)
      var_name='dpbgrid'; call wrtarr(var_name,ax)
c****
      title='STANDARD DEVIATION OF PRESSURE DIFFERENCES  (MB)'
      acc_name='AJK23'; call getacc(acc_name,ax)
      acc_name='AJK24'; call getacc(acc_name,bx)
      DO J=2,JM
         DO K=1,KM
            BYN=1./(BX(J,K)+1.D-10)
            AX(J,K)=0.
            SDDP=(AX(J,K)-DPAB(J,K,2)*DPAB(J,K,2)*BYN)*BYN
            IF (SDDP.GT.0.) AX(J,K)=SQRT(SDDP)
         ENDDO
      ENDDO
      SCALE=1.
      call scale2d(ax,scale,ones,ones,jm,lm)
      var_name='sddp'; call wrtarr(var_name,ax)
c****
! (re)set shape of output arrays
      dim_name='latitude'; call set_dim_out(dim_name,1)
      dim_name='prqt'; call set_dim_out(dim_name,2)

      title='TEMPERATURE (DEGREES CENTIGRADE)'
      acc_name='AJK03'; call getacc(acc_name,arqx(1,1))
      call scalejk(arqx(1,1),1.,byp,ones,dsjk(1,1,1),jm,lm,missing)
      acc_name='ASJL01'; call getacc(acc_name,arqx(1,lm+1))
      SCALE=BYIMDA
      call scale2d(arqx(1,lm+1),scale,onespo,ones,jm,lm_req)
      var_name='temp'; call wrtarr(var_name,arqx)
c****
      title='POTENTIAL TEMPERATURE AT 1000 MB (DEGREES KELVIN)'
      acc_name='AJK06'; call getacc(acc_name,arqx(1,1))
      SCALE=P1000K
      call scalejk(arqx(1,1),scale,byp,ones,dsjk(1,1,1),jm,lm,missing)
      acc_name='ASJL01'; call getacc(acc_name,arqx(1,lm+1))
      ARQX(:,lm+1:lm+lm_req)=ARQX(:,lm+1:lm+lm_req)*
     &     BYIMDA*SPREAD(ONESPO,DIM=2,NCOPIES=LM_REQ)+TF
      call scale2d(arqx(1,lm+1),scale,ones,bypks,jm,lm_req)
      var_name='pt1000'; call wrtarr(var_name,arqx)
c****
      title='HEIGHT (HUNDREDS OF METERS)'
      acc_name='AJK04'; call getacc(acc_name,arqx(1,1))
      SCALE=BY100G
      call scalejk(arqx(1,1),scale,byp,ones,dsjk(1,1,1),jm,lm,missing)
      acc_name='ASJL02'; call getacc(acc_name,arqx(1,lm+1))
      SCALE=BYIMDA*BY100G
      call scale2d(arqx(1,lm+1),scale,onespo,ones,jm,lm_req)
      var_name='height'; call wrtarr(var_name,arqx)
c****
! (re)set shape of output arrays
      dim_name='latitude'; call set_dim_out(dim_name,1)
      dim_name='p'; call set_dim_out(dim_name,2)

      title='SPECIFIC HUMIDITY (KG H2O/KG AIR)'
      acc_name='AJK05'; call getacc(acc_name,ax)
      call scalejk(ax,1.,byp,ones,dsjk(1,1,1),jm,lm,missing)
      var_name='q'; call wrtarr(var_name,ax)
c****
      title='RELATIVE HUMIDITY'
      acc_name='AJK07'; call getacc(acc_name,ax)
      call scalejk(ax,1.,byp,ones,dsjk(1,1,1),jm,lm,missing)
      var_name='rh'; call wrtarr(var_name,ax)
c****
      title='TOTAL CLOUD WATER CONTENT (KG/KG)'
      acc_name='AJK51'; call getacc(acc_name,ax)
      call scalejk(ax,1.,byp,ones,dsjk(1,1,1),jm,lm,missing)
      var_name='cldh2o'; call wrtarr(var_name,ax)
c****
! (re)set shape of output arrays
      dim_name='latb'; call set_dim_out(dim_name,1)
      dim_name='p'; call set_dim_out(dim_name,2)

      title='ZONAL WIND (U COMPONENT)'
      acc_name='AJK08'; call getacc(acc_name,ax)
      call scalejk(ax,1.,bypv,ones,dsjk(1,1,2),jm,lm,missing)
      axb(:,:)=ax(2:jm,:)
      var_name='u'; call wrtarr(var_name,axb)
c****
      title='MERIDIONAL WIND (V COMPONENT)'
      acc_name='AJK09'; call getacc(acc_name,ax)
      call scalejk(ax,1.,bypv,ones,dsjk(1,1,2),jm,lm,missing)
      axb(:,:)=ax(2:jm,:)
      var_name='v'; call wrtarr(var_name,axb)
c****
! (re)set shape of output arrays
      dim_name='latb'; call set_dim_out(dim_name,1)
      dim_name='sige'; call set_dim_out(dim_name,2)

      title='STREAM FUNCTION (10**9 KILOGRAMS/SECOND) CP'
      acc_name='AJK09'; call getacc(acc_name,ax)
      DO J=2,JM
         AX(J,1)=AX(J,1)
         DO K=2,KM
            AX(J,K)=AX(J,K)+AX(J,K-1)
         ENDDO
      ENDDO
      SCALE=100.D-9*XWON*BYIADA*BYGRAV
      call scale2d(ax,scale,dxv,ones,jm,lm)
      axb(:,:)=ax(2:jm,:)
      var_name='stfunc'; call wrtarr(var_name,axb)
c****
! (re)set shape of output arrays
      dim_name='latitude'; call set_dim_out(dim_name,1)
      dim_name='pe'; call set_dim_out(dim_name,2)

      title='VERTICAL VELOCITY (MB/SECOND)'
      acc_name='AJK25'; call getacc(acc_name,ax)
      SCALE=-BYIMDA
      call scale2d(ax,scale,bydxyp,ones,jm,lm)
      var_name='vvel'; call wrtarr(var_name,ax)
c****
c close input/output files
      call close_acc
      call close_out
c****
      end program jldag

      subroutine scale2d(arr,scale,scalej,scalel,jm,lm)
      implicit none
      integer :: jm,lm
      real, dimension(jm,lm) :: arr
      real :: scale
      real, dimension(jm) :: scalej
      real, dimension(lm) :: scalel
      integer :: j,l
      do l=1,lm
         do j=1,jm
            arr(j,l)=arr(j,l)*scale*scalej(j)*scalel(l)
         enddo
      enddo
      return
      end subroutine scale2d

      subroutine scalejk(arr,scale,scalej,scalel,massjl,jm,lm,missing)
      implicit none
      integer :: jm,lm
      real, dimension(jm,lm) :: arr,massjl
      real :: scale,missing
      real, dimension(jm) :: scalej
      real, dimension(lm) :: scalel
      integer :: j,l
      do l=1,lm
         do j=1,jm
            arr(j,l)=arr(j,l)*scale*scalej(j)*scalel(l)
         enddo
      enddo
      where(massjl.ne.missing)
         arr(:,:)=arr(:,:)/massjl(:,:)
      elsewhere
         arr(:,:)=missing
      end where
      return
      end subroutine scalejk
