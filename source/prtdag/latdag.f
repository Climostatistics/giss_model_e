c
      program latdag

      use ncinp
      use ncout
      implicit none

      integer :: nargs,iargc,iarg
      character(len=80) :: stop_str,cmd_str

      real :: scale,A1BYA2,A2BYA1,BYA1
      real :: dtday,dtsec

      integer :: ntype,stype_clen
      integer :: itype,ic
      real, dimension(:), allocatable :: lats,sige,dxyp
      real, dimension(:,:), allocatable :: acc
      real, dimension(:), allocatable :: acc1,acc2,acc3,acc4,acc5,accsum
      real, dimension(:), allocatable :: fj,sptype
      real, dimension(:), allocatable :: wtype
      character, dimension(:,:), allocatable :: stype_names
      real, parameter :: P1000=1000.
      character(len=20) :: var_name,acc_name,dim_name,type_name

      call getarg(0,cmd_str)
      nargs = iargc()
      if(nargs.lt.2) then
         write(6,*) 'usage: '//trim(cmd_str)//' acc_file out_file'//
     &        ' [type-name1 [type-name2 ... for composite] ]'
         write(6,*) 'omit type names for composite over all types'
         stop
      endif
      call getarg(1,accfile)
      call getarg(2,outfile)
      if(accfile.eq.outfile) stop 'cannot overwrite input file'

! open the acc file
      call open_acc

! get ntype dimensions
      dim_name='STYPE_CLEN'; call get_dim_size(dim_name,stype_clen)
      dim_name='NTYPE'; call get_dim_size(dim_name,ntype)
      if(nargs-2.gt.ntype) stop 'too many surface types in arguments'
! allocate space using these dimensions
      allocate(lats(jm),dxyp(jm))
      allocate(sige(lm+1))
      allocate(acc(jm,ntype))
      allocate (acc1(jm),acc2(jm),acc3(jm),acc4(jm),acc5(jm),accsum(jm))
      allocate(fj(jm),sptype(jm))
      allocate(wtype(ntype))
      allocate(stype_names(stype_clen,ntype))
! get lon,ht coordinates
      acc_name='latitude'; call getacc(acc_name,lats)
      acc_name='area'; call getacc(acc_name,dxyp)
      acc_name='sige'; call getacc(acc_name,sige)

! get names of surface types
      acc_name='STYPE_NAMES'; call gettxt(acc_name,stype_names)

C**** determine surface type weightings from arguments
      wtype=0.
      if(nargs.eq.2) then ! composite of all types
         wtype = 1.
      else
         do iarg=3,nargs
            type_name=''
            call getarg(iarg,type_name)
            itype = 1
            do while(itype.le.ntype)
               ic=1
               do while(ic.le.stype_clen .and.
     &                  type_name(ic:ic).eq.stype_names(ic,itype))
                  ic = ic +1
               enddo
               if(ic.le.stype_clen) then
                  itype = itype + 1
               else
                  exit
               endif
            enddo
            if(itype.gt.ntype) then
               write(6,*) 'invalid type: ',trim(type_name)
               stop
            endif
            wtype(itype)=1.
         enddo
      endif

! define output file
      call open_out
      dim_name='latitude'; call def_dim_out(dim_name,jm)

      ndims_out = 1
      dim_name='latitude'; call set_dim_out(dim_name,1)
      units='degrees_north'
      var_name='latitude';call wrtarr(var_name,lats)

! find total integration time in seconds, days, and other units
      dtsec = idacc(1)*3600.
      dtday = idacc(1)/24.

C****
      BYA1=1./(IDACC(1)+1.E-20)
      A2BYA1=DBLE(IDACC(2))/DBLE(IDACC(1))
      A1BYA2=IDACC(1)/(IDACC(2)+1.E-20)

C****
c  AJ69
      var_name = 'ptype'
      long_name = 'SURF TYPE FRACT'
      units = '1'
      acc_name='PTYPE'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1. / idacc(ia)
      sptype = accsum*scale
      fj = sptype
      fj(2:jm-1) = fj(2:jm-1)/im
      call wrtarr(var_name,fj)
c grid areas
      var_name='area'
      long_name = 'area of latitude band'
      units = 'm^2'
      call wrtarr(var_name,dxyp)
c  AJ01
      var_name = 'inc_sw'
      long_name = 'INC SW'
      units = 'W/m^2'
      acc_name='SRINCP0'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ02
      var_name = 'srnfp0'
      long_name = 'SW ABS BELOW P0'
      units = '1'
      acc_name='SRNFP0'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ03
      var_name = 'srnfp1'
      long_name = 'SW ABS BELOW P1'
      units = 'W/m^2'
      acc_name='SRNFP1'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ04
      var_name = 'swabsatm'
      long_name = 'SW ABS BY ATMOS'
      units = 'W/m^2'
      acc_name='SRNFP0'; call getaj(acc_name,acc,acc1,jm,wtype,ntype)
      acc_name='SRNFG'; call getaj(acc_name,acc,acc2,jm,wtype,ntype)
      accsum = acc1 - acc2
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ05
      var_name = 'sw_inc_z0'
      long_name = 'SW INC ON ZO'
      units = 'W/m^2'
      acc_name='SRINCG'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ06
      var_name = 'sw_abs_z0'
      long_name = 'SW ABS AT ZO'
      units = 'W/m^2'
      acc_name='SRNFG'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ07
      var_name = 'lw_net_p0'
      long_name = 'NET LW AT P0'
      units = 'W/m^2'
      acc_name='HSURF'; call getaj(acc_name,acc,acc1,jm,wtype,ntype)
      acc_name='TRHDT'; call getaj(acc_name,acc,acc2,jm,wtype,ntype)
      accsum  = acc1+A2BYA1*acc2/DTSRC
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ08
      var_name = 'lw_net_p1'
      long_name = 'NET LW AT P1'
      units = 'W/m^2'
      acc_name='HATM'; call getaj(acc_name,acc,acc1,jm,wtype,ntype)
      acc_name='TRHDT'; call getaj(acc_name,acc,acc2,jm,wtype,ntype)
      accsum = acc1 +A2BYA1*acc2/DTSRC
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ09
      var_name = 'lw_net_z0'
      long_name = 'NET LW AT ZO (positive downward)'
      units = 'W/m^2'
      acc_name='TRHDT'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1./DTSRC / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ10
      var_name = 'rnet_p0'
      long_name = 'NET RAD AT P0'
      units = 'W/m^2'
      acc_name='SRNFP0'; call getaj(acc_name,acc,acc1,jm,wtype,ntype)
      acc_name='TRNFP0'; call getaj(acc_name,acc,acc2,jm,wtype,ntype)
      accsum = acc1+acc2
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ11
      var_name = 'rnet_p1'
      long_name = 'NET RAD AT P1'
      units = 'W/m^2'
      acc_name='SRNFP1'; call getaj(acc_name,acc,acc1,jm,wtype,ntype)
      acc_name='TRNFP1'; call getaj(acc_name,acc,acc2,jm,wtype,ntype)
      accsum = acc1+acc2
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ12
      var_name = 'rnet_z0'
      long_name = 'NET RAD AT Z0'
      units = 'W/m^2'
      acc_name='SRNFG'; call getaj(acc_name,acc,acc1,jm,wtype,ntype)
      acc_name='TRHDT'; call getaj(acc_name,acc,acc2,jm,wtype,ntype)
      accsum = A1BYA2*acc1*DTSRC+acc2
      scale = 1./DTSRC / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ13
      var_name = 'shflx'
      long_name = 'SENS HEAT FLUX'
      units = 'W/m^2'
      acc_name='SHDT'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1./DTSRC / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ14
      var_name = 'lhflx'
      long_name = 'LATENT HEAT FLUX'
      units = 'W/m^2'
      acc_name='EVHDT'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1./DTSRC / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ15
      var_name = 'f2dt'
      long_name = 'CONDC AT -Z1-Z2'
      units = 'W/m^2'
      acc_name='F2DT'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1./DTSRC / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ16
      var_name = 'htz1'
      long_name = 'NET HEAT AT -Z1'
      units = 'W/m^2'
      acc_name='EDIFS'; call getaj(acc_name,acc,acc1,jm,wtype,ntype)
      acc_name='F1DT'; call getaj(acc_name,acc,acc2,jm,wtype,ntype)
      accsum = acc1+acc2
      scale = 1./DTSRC / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ17
      var_name = 'tg2'
      long_name = 'TG2'
      units = 'degC'
      acc_name='TG2'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ18
      var_name = 'tg1'
      long_name = 'TG1'
      units = 'degC'
      acc_name='TG1'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ19
      var_name = 'evap'
      long_name = 'EVAP'
      units = 'mm/day'
      acc_name='EVAP'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = SDAY/DTSRC / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ20
      var_name = 'prec'
      long_name = 'PREC'
      units = 'mm/day'
      acc_name='PRCPSS'; call getaj(acc_name,acc,acc1,jm,wtype,ntype)
      acc_name='PRCPMC'; call getaj(acc_name,acc,acc2,jm,wtype,ntype)
      accsum = acc1+acc2
      scale = 100.*SDAY/(DTSRC*GRAV) / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ21
      var_name = 'tair'
      long_name = 'T AIR'
      units = 'degC'
      acc_name='TX'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ22
      var_name = 't1'
      long_name = 'T1'
      units = 'degC'
      acc_name='TX1'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ23
      var_name = 'tsurf'
      long_name = 'T SURF'
      units = 'degC'
      acc_name='TSRF'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ24
      var_name = 'stsb_st'
      long_name = 'STAT STB(STRAT)'
      units = 'unknown'
      acc_name='DTSGST'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1.D3*GRAV*(P1000)**KAPA / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ25
      var_name = 'stsb_tr'
      long_name = 'STAT STB(TROPO)'
      units = 'unknown'
      acc_name='DTDGTR'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1.D3*GRAV*(P1000)**KAPA / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ26
      var_name = 'rich_st'
      long_name = 'RICH NUM (STRATOSPHERE)'
      units = '1'
      acc_name='RICST'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 16.*RGAS / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ27
      var_name = 'rich_tr'
      long_name = 'RICH NUM (TROPOSPHERE)'
      units = '1'
      acc_name='RICTR'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 16.*RGAS / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ28
      var_name = 'ross_st'
      long_name = 'ROSS NUM (STRATOSPHERE)'
      units = '1'
      acc_name='ROSST'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = .5/(2.*OMEGA*DBLE(IM)) / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ29
      var_name = 'ross_tr'
      long_name = 'ROSS NUM (TROPOSPHERE)'
      units = '1'
      acc_name='ROSTR'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = .5/(2.*OMEGA*DBLE(IM)) / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ30
      var_name = 'poice'
      long_name = 'OC/LK ICE COVER'
      units = '1'
      acc_name='RSI'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ31
      var_name = 'psnow'
      long_name = 'SNOW COVER'
      units = '1'
      acc_name='RSNOW'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ33
      var_name = 'oht'
      long_name = 'OCEAN TRANSPORT'
      units = 'W/m^2'
      acc_name='OHT'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1./DTSRC / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ34
      var_name = 'tg3'
      long_name = 'OCEAN TEMP AT MAX. MIXED LAYER DEPTH'
      units = 'degC'
      acc_name='OMLT'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ35
      var_name = 'dtdlt_st'
      long_name = 'DT/DLAT(STRAT)'
      units = 'K'
      acc_name='DTDJS'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = .5D2*(JM-1)/((SIGE(LS1)-SIGE(LM+1)+1.E-20)*180.)/idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ36
      var_name = 'dtdlt_tr'
      long_name = 'DT/DLAT(TROPO)'
      units = 'K'
      acc_name='DTDJT'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = .5E2*(JM-1.)/((SIGE(1)-SIGE(LS1))*180.) / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ37
      var_name = 'lross_st'
      long_name = 'L(STRAT)(10**5)'
      units = 'm'
      acc_name='LSTR'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1.D-5*SQRT(RGAS)/(2.*OMEGA) / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ38
      var_name = 'lross_tr'
      long_name = 'L(TROP) (10**5)'
      units = 'm'
      acc_name='LTRO'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1.D-5*SQRT(RGAS)/(2.*OMEGA) / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ39
      var_name = 'ENERGP'
      long_name = 'PRECIP HEAT FLX'
      units = 'W/m^2'
      acc_name='EPRCP'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1./DTSRC / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ40
      var_name = 'ERUN1'
      long_name = 'HEAT RUNOFF Z0'
      units = 'W/m^2'
      acc_name='ERUN1'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1./DTSRC / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ41
      var_name = 'EDIFS'
      long_name = 'HT WTR DIFS -Z1'
      units = 'W/m^2'
      acc_name='EDIFS'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1./DTSRC / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ42
      var_name = 'F1DT'
      long_name = 'CONDUCTN AT -Z1'
      units = 'W/m^2'
      acc_name='F1DT'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1./DTSRC / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ43
      var_name = 'ERUN2'
      long_name = 'ICE ENRG -Z1-Z2'
      units = 'W/m^2'
      acc_name='ERUN2'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1./DTSRC / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ44
      var_name='hnet_z0'
      long_name = 'NET HEAT AT Z0'
      units = 'W/m^2'
      acc_name='RHDT'; call getaj(acc_name,acc,acc1,jm,wtype,ntype)
      acc_name='SHDT'; call getaj(acc_name,acc,acc2,jm,wtype,ntype)
      acc_name='EVHDT'; call getaj(acc_name,acc,acc3,jm,wtype,ntype)
      acc_name='EPRCP'; call getaj(acc_name,acc,acc4,jm,wtype,ntype)
      acc_name='ERUN1'; call getaj(acc_name,acc,acc5,jm,wtype,ntype)
      accsum = acc1+acc2+acc3+acc4-acc5
      scale = 1./DTSRC / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ45
      var_name = 'DIFS'
      long_name = 'H2O DIFS AT -Z1'
      units = 'mm/day'
      acc_name='DIFS'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = SDAY/DTSRC / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ46
      var_name = 'ice_z1z2'
      long_name = 'ICE THRU -Z1-Z2'
      units = 'mm/day'
      acc_name='IMELT'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = SDAY/DTSRC / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ47
      var_name = 'RUN2'
      long_name = 'WATR RUNOFF MLD'
      units = 'mm/day'
      acc_name='RUN2'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = SDAY/DTSRC / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ48
      var_name = 'DWTR2'
      long_name = 'HEAT RUNOFF THROUGH THE MIXED LAYER DEPTH'
      units = 'W/m^2'
      acc_name='DWTR2'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1./DTSRC / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ49
      var_name = 'WTR1'
      long_name = 'WATER IN G1'
      units = 'kg/m^2'
      acc_name='WTR1'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ50
      var_name = 'ICE1'
      long_name = 'ICE IN G1'
      units = 'kg/m^2'
      acc_name='ACE1'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ51
      var_name = 'WTR2'
      long_name = 'WATER IN G2'
      units = 'kg/m^2'
      acc_name='WTR2'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ52
      var_name = 'ICE2'
      long_name = 'ICE IN G2'
      units = 'kg/m^2'
      acc_name='ACE2'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ53
      var_name = 'SNOWDP'
      long_name = 'SNOW DEPTH'
      units = 'kg/m^2'
      acc_name='SNOW'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ54
      var_name = 'RUN1'
      long_name = 'WATER RUNOFF Z0'
      units = 'mm/day'
      acc_name='RUN1'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = SDAY/DTSRC / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ55
      var_name = 'BTEMPW'
      long_name = 'LW WINDOW BTEMP'
      units = 'degC'
      acc_name='BRTEMP'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ56
      var_name = 'htz2'
      long_name = 'NET HEAT -Z1-Z2'
      units = 'W/m^2'
      acc_name='F2DT'; call getaj(acc_name,acc,acc1,jm,wtype,ntype)
      acc_name='ERUN2'; call getaj(acc_name,acc,acc2,jm,wtype,ntype)
      accsum = acc1+acc2
      scale = 1./DTSRC / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ57
      var_name = 'sscld'
      long_name = 'TOT SUP SAT CLD'
      units = '1'
      acc_name='PCLDSS'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ58
      var_name = 'mccld'
      long_name = 'TOT MST CNV CLD'
      units = '1'
      acc_name='PCLDMC'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ59
      var_name = 'pcld'
      long_name = 'TOTAL CLD COVER'
      units = '1'
      acc_name='PCLD'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ60
      var_name = 'cldtopmc'
      long_name = 'MC CLD DPTH'
      units = '100 PA'
      acc_name='CDLDEP'; call getaj(acc_name,acc,acc1,jm,wtype,ntype)
      acc_name='PCLDMC'; call getaj(acc_name,acc,acc2,jm,wtype,ntype)
      accsum = acc1/(acc2+1.E-20)
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ61
      var_name = 'PRCPSS'
      long_name = 'SS PRECIP'
      units = 'mm/day'
      acc_name='PRCPSS'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 100.*SDAY/(DTSRC*GRAV) / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ62
      var_name = 'PRCPMC'
      long_name = 'MC PRECIP'
      units = 'mm/day'
      acc_name='PRCPMC'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 100.*SDAY/(DTSRC*GRAV) / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ63
      var_name = 'atmh2o'
      long_name = 'H2O OF ATM (MM)'
      units = 'kg/m^2'
      acc_name='QP'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 100./GRAV / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ64
      var_name = 'GAM'
      long_name = 'GAM'
      units = 'K/KM'
      acc_name='GAM'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1.D3*GRAV / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ65
      var_name = 'GAMM'
      long_name = 'GAMM'
      units = 'K/KM'
      acc_name='GAMM'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1.D3*.0098/(SIGE(1)-SIGE(LS1)) / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ66
      var_name = 'GAMC'
      long_name = 'GAMC'
      units = 'K/KM'
      acc_name='GAMC'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1.D3 / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ67
      var_name = 'lw_inc_z0'
      long_name = 'LW INC ON ZO'
      units = 'W/m^2'
      acc_name='TRINCG'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 1. / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
c  AJ68
      var_name = 'HT_THRMCL'
      long_name = 'ENERGY DIFFUSION INTO THERMOCLINE'
      units = 'W/m^2'
      acc_name='FTHERM'; call getaj(acc_name,acc,accsum,jm,wtype,ntype)
      scale = 2.E3*4185./SDAY / idacc(ia)
      fj = accsum*scale/(sptype+1.d-20)
      call wrtarr(var_name,fj)
C****
C**** CALCULATE AND PRINT ALBEDOS
C****
      acc_name='SRINCP0';call getaj(acc_name,acc,acc2,jm,wtype,ntype)
      acc_name='SRINCG'; call getaj(acc_name,acc,acc3,jm,wtype,ntype)
      scale = 1.
c AJ02/AJ01
      var_name = 'palb'
      long_name = 'PLANETARY ALBDO'
      units = '1'
      acc_name='SRNFP0'; call getaj(acc_name,acc,acc1,jm,wtype,ntype)
      fj = 1.-acc1/(acc2+1.E-20)
      call wrtarr(var_name,fj)
c AJ72/AJ01
      var_name = 'palb_vis'
      long_name = 'PLAN ALB VISUAL'
      units = '1'
      acc_name='PLAVIS'; call getaj(acc_name,acc,acc1,jm,wtype,ntype)
      fj = acc1/(acc2+1.E-20)
      call wrtarr(var_name,fj)
c AJ73/AJ01
      var_name = 'palb_nir'
      long_name = 'PLAN ALB NEARIR'
      units = '1'
      acc_name='PLANIR'; call getaj(acc_name,acc,acc1,jm,wtype,ntype)
      fj = acc1/(acc2+1.E-20)
      call wrtarr(var_name,fj)
c AJ06/AJ05
      var_name = 'galb'
      long_name = 'SURFACE G ALBDO'
      units = '1'
      acc_name='SRNFG'; call getaj(acc_name,acc,acc1,jm,wtype,ntype)
      fj = 1.-acc1/(acc3+1.E-20)
      call wrtarr(var_name,fj)
c AJ74/AJ01
      var_name = 'galb_vis'
      long_name = 'SURF ALB VISUAL'
      units = '1'
      acc_name='ALBVIS'; call getaj(acc_name,acc,acc1,jm,wtype,ntype)
      fj = acc1/(acc2+1.E-20)
      call wrtarr(var_name,fj)
c AJ75/AJ01
      var_name = 'galb_nir'
      long_name = 'SURF ALB NEARIR'
      units = '1'
      acc_name='ALBNIR'; call getaj(acc_name,acc,acc1,jm,wtype,ntype)
      fj = acc1/(acc2+1.E-20)
      call wrtarr(var_name,fj)
c AJ76/AJ01
      var_name = 'aalb_vis'
      long_name = 'ATMO ALB VISUAL'
      units = '1'
      acc_name='SRRVIS'; call getaj(acc_name,acc,acc1,jm,wtype,ntype)
      fj = acc1/(acc2+1.E-20)
      call wrtarr(var_name,fj)
c AJ77/AJ01
      var_name = 'aalb_nir'
      long_name = 'ATMO ALB NEARIR'
      units = '1'
      acc_name='SRRNIR'; call getaj(acc_name,acc,acc1,jm,wtype,ntype)
      fj = acc1/(acc2+1.E-20)
      call wrtarr(var_name,fj)
c AJ78/AJ01
      var_name = 'aabs_vis'
      long_name = 'ATMO ABS VISUAL'
      units = '1'
      acc_name='SRAVIS'; call getaj(acc_name,acc,acc1,jm,wtype,ntype)
      fj = acc1/(acc2+1.E-20)
      call wrtarr(var_name,fj)
c AJ79/AJ01
      var_name = 'aabs_nir'
      long_name = 'ATMO ABS NEARIR '
      units = '1'
      acc_name='SRANIR'; call getaj(acc_name,acc,acc1,jm,wtype,ntype)
      fj = acc1/(acc2+1.E-20)
      call wrtarr(var_name,fj)
c
      call close_acc
      call close_out

      end program latdag

      subroutine getaj(acc_name,acc,acc_sum,dim1,wtype,ntype)
      use ncinp
      implicit none
      integer :: dim1,ntype
      character(len=20) :: acc_name
      real, dimension(dim1,ntype) :: acc
      real, dimension(dim1) :: acc_sum
      real, dimension(ntype) :: wtype
      integer :: itype
      character(len=20) :: aj_name

      aj_name='J_'//trim(acc_name)
      call getacc(aj_name,acc)

      acc_sum = 0.
      do itype=1,ntype
         acc_sum = acc_sum + wtype(itype)*acc(:,itype)
      enddo

      return
      end subroutine getaj

C****
C**** THIS SUBROUTINE PRODUCES AREA WEIGHTED STATISTICS OF
C****
C   K   N
C****
C***1   1  SOLAR RADIATION INCIDENT ON PLANET (W/M**2)
C****
C**1A   2/1  PLANETARY ALBEDO (10**-2)
C**1B  72/1  PLANETARY ALBEDO VISUAL (10**-2)
C**1C  73/1  PLANETARY ALBEDO NEAR IR (10**-2)
C**1D   6/5  GROUND ALBEDO (10**-2)
C**1E  74/1  GROUND ALBEDO VISUAL (10**-2)
C**1F  75/1  GROUND ALBEDO NEAR IR (10**-2)
C**1G  76/1  ATMOSPHERIC ALBEDO VISUAL (10**-2)
C**1H  77/1  ATMOSPHERIC ALBEDO NEAR IR (10**-2)
C**1I  78/1  ATMOSPHERIC ABSORPTION VISUAL (10**-2)
C**1J  79/1  ATMOSPHERIC ABSORPTION NEAR IR (10**-2)
C****
C***2   2  SOLAR RADIATION ABSORBED BY PLANET (W/M**2)
C***3   3  SOLAR RADIATION ABSORBED BELOW PTOP (W/M**2)
C***4   4  SOLAR RADIATION ABSORBED BY ATMOSPHERE (W/M**2)
C***5   5  SOLAR RADIATION INCIDENT ON GROUND (W/M**2)
C***6   6  SOLAR RADIATION ABSORBED BY GROUND (W/M**2)
C***7  32  SOLAR RADIATION WATER CORRECTION
C***8   7  THERMAL RADIATION EMITTED BY PLANET (W/M**2)
C***9   8  THERMAL RADIATION AT PTOP (W/M**2)
C**10   9  THERMAL RADIATION EMITTED BY GROUND (W/M**2)
C****
C**11  67  THERMAL RADIATION INCIDENT ON GROUND (W/M**2)
C****  55  BRIGHTNESS TEMPERATURE THROUGH WINDOW REGION (K-TF)
C****  10  NET RADIATION ABSORBED BY PLANET (W/M**2)
C****  11  NET RADIATION ABSORBED BELOW PTOP (W/M**2)
C****  12  NET RADIATION ABSORBED BY GROUND (W/M**2)
C****  13  SENSIBLE HEAT FLUX INTO THE GROUND (W/M**2)
C****  14  EVAPORATION HEAT FLUX INTO THE GROUND (W/M**2)
C****  39  PRECIPITATION HEAT FLUX INTO THE GROUND (W/M**2)
C****  40  HEAT RUNOFF FROM FIRST GROUND LAYER (W/M**2)
C****  44  NET HEATING AT Z0 (W/M**2)
C****
C**21  42  CONDUCTION AT -Z1 (W/M**2)
C****  41  HEAT OF WATER OR ICE DUFFUSION AT -Z1 (W/M**2)
C****  16  NET HEATING AT -Z1 (W/M**2)
C****  15  CONDUCTION AT -Z1-Z2 (W/M**2)
C****  43  ENERGY OF ICE MELTING (OR TRANSPORTING) AT -Z1-Z2 (W/M**2)
C****  56  NET HEATING AT -Z1-Z2 (W/M**2)
C****  33  OCEAN TRANSPORT (W/M**2)
C****  48  HEAT RUNOFF THROUGH THE MIXED LAYER DEPTH (W/M**2)
C****  68  ENERGY DIFFUSION INTO THE THERMOCLINE (W/M**2)
C****  18  MEAN TEMPERATURE OF FIRST GROUND LAYER (.1 K-TF)
C****
C**31  17  MEAN TEMPERATURE OF SECOND GROUND LAYER (.1 K-TF)
C****  34  OCEAN TEMPERATURE AT THE MAXIMUM MIXED LAYER DEPTH
C****  23  SURFACE AIR TEMPERATURE (.1 K-TF)
C****  22  FIRST LAYER AIR TEMPERATURE (.1 K-TF)
C****  21  COMPOSITE AIR TEMPERATURE (.1 K-TF)
C****  35  STRATO TEMPERATURE CHANGE PER DEGREE LATITUDE (10**-2 K)
C****  36  TROPO TEMPERATURE CHANGE PER DEGREE LATITUDE (10**-2 K)
C****  24  STRATOSPHERIC STATIC STABILITY (10**-3 K/M)
C****  25  TROPOSPHERIC STATIC STABILITY (10**-3 K/M)
C****  26  STRATOSPHERIC RICHARDSON NUMBER (1)
C****
C**41  27  TROPOSPHERIC RICHARDSON NUMBER (1)
C****  28  STRATOSPHERIC ROSSBY NUMBER (1)
C****  29  TROPOSPHERIC ROSSBY NUMBER (1)
C****  37  L IN THE STRATOSPHERE (10**5 M)
C****  38  L IN THE TROPOSPHERE (10**5 M)
C****  64  GAM  (10**-3 K/M)
C****  65  GAMM  (10**-3 K/M)
C****  66  GAMC  (10**-3 K/M)
C****  57  INTEGRATED SUPER-SATURATION CLOUD COVER (10**-2)
C****  58  INTEGRATED MOIST CONVECTIVE CLOUD COVER (10**-2)
C****
C**51  59  INTEGRATED TOTAL CLOUD COVER (10**-2)
C****  60  MOIST CONVECTIVE CLOUD DEPTH (100 N)
C****  61  SUPER SATURATION PRECIPITATION (KG/M**2/86400 S)
C****  62  MOIST CONVECTIVE PRECIPITATION (KG/M**2/86400 S)
C****  20  PRECIPITATION (KG/M**2/86400 S)
C****  19  EVAPORATION (KG/M**2/86400 S)
C****  63  WATER CONTENT OF ATMOSPHERE (KG/M**2)
C****  54  WATER RUNOFF AT Z0 (KG/M**2/86400 S)
C****  45  WATER OR ICE DIFFUSION AT -Z1 (KG/M**2/86400 S)
C****  46  ICE MELTING (OR TRANSPORTING) AT -Z1-Z2 (KG/M**2/86400 S)
C****
C**61  47  WATER RUNOFF THROUGH MIXED LAYER DEPTH (KG/M**2/86400 S)
C****  49  WATER CONTAINED IN FIRST GROUND LAYER (KG/M**2)
C****  50  ICE CONTAINED IN FIRST GROUND LAYER (KG/M**2)
C****  51  WATER CONTAINED IN SECOND GROUND LAYER (KG/M**2)
C****  52  ICE CONTAINED IN SECOND GROUND LAYER (KG/M**2)
C****  53  SNOW DEPTH (KG/M**2)
C****  31  SNOW COVER (10**-2)
C**68  30  OCEAN ICE COVER (10**-2)
C****

cc  AJ32
c      scale = 1.
c      var_name = 'SW CORRECTION'
c      long_name = 'SW CORRECTION'
c      units = 'unknown'
c      acc_name='SWCOR'; call getaj(acc_name,acc,accsum,ia)
