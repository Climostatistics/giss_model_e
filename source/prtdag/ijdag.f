      program ijdag
C****
C****
      use ncinp
      use ncout
      implicit none
      integer :: nargs,iargc
      character(len=80) :: stop_str,cmd_str
      real :: dtsec,dtday
      real, dimension(:), allocatable :: lons,lats
      real, dimension(:,:), allocatable :: acc,acc2
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

! allocate space
      allocate(lons(im))
      allocate(lats(jm))
      allocate(acc(im,jm),acc2(im,jm))

! get lon,lat coordinates
      acc_name='longitude'; call getacc(acc_name,lons)
      acc_name='latitude'; call getacc(acc_name,lats)

! define output file
      call open_out
      dim_name='longitude'; call def_dim_out(dim_name,im)
      dim_name='latitude'; call def_dim_out(dim_name,jm)

      ndims_out = 1
      dim_name='longitude'; call set_dim_out(dim_name,1)
      units='degrees_east'
      var_name='longitude';call wrtarr(var_name,lons)
      dim_name='latitude'; call set_dim_out(dim_name,1)
      units='degrees_north'
      var_name='latitude';call wrtarr(var_name,lats)

! find total integration time in seconds, days, and other units
      dtsec = idacc(1)*3600.
      dtday = idacc(1)/24.

! set shape of output arrays
      ndims_out = 2
      dim_name='longitude'; call set_dim_out(dim_name,1)
      dim_name='latitude'; call set_dim_out(dim_name,2)

! ACCUMULATED OCEAN ICE FRACTION units = '1'
      acc_name='RSOI'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='oicefr'; call wrtarr(var_name,acc)

! ACCUMULATED SNOW FRACTION units = '1'
      acc_name='RSNW'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='snowfr'; call wrtarr(var_name,acc)

! ACCUMULATED SNOW MASS units = 'kg/m^2'
      acc_name='SNOW'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='snowdp'; call wrtarr(var_name,acc)

! INTEGRATED SENSIBLE HEAT FLUX units = 'J/m^2'
      acc_name='SHDT'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/dtsec
      var_name='sensht'; call wrtarr(var_name,acc)

! ACCUMULATED PRECIPITATION units = 'kg/m^2'
      acc_name='PREC'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/dtday
      var_name='precip'; call wrtarr(var_name,acc)

! ACCUMULATED EVAPORATION units = 'kg/m^2'
      acc_name='EVAP'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/dtday
      var_name='evapor'; call wrtarr(var_name,acc)

! ACCUMULATED EVAPORATION EFFICIENCY units = '1'
      acc_name='BETA'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='evbeta'; call wrtarr(var_name,acc)

! ACCUMULATED SURFACE PRESSURE units = '100 PA'
      acc_name='PRES'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='prsurf'; call wrtarr(var_name,acc)

!  PHI1000 units = 'm^2/s^2'
      acc_name='PHI1K'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='phi1000'; call wrtarr(var_name,acc)

!  PHI850 units = 'm^2/s^2'
      acc_name='PHI850'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='phi850'; call wrtarr(var_name,acc)

! PHI700 units = 'm^2/s^2'
      acc_name='PHI700'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='phi700'; call wrtarr(var_name,acc)

! PHI500 units = 'm^2/s^2'
      acc_name='PHI500'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='phi500'; call wrtarr(var_name,acc)

! PHI300 units = 'm^2/s^2'
      acc_name='PHI300'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='phi300'; call wrtarr(var_name,acc)

! PHI100 units = 'm^2/s^2'
      acc_name='PHI100'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='phi100'; call wrtarr(var_name,acc)

! PHI30 units = 'm^2/s^2'
      acc_name='PHI30'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='phi30'; call wrtarr(var_name,acc)

! T850 units = 'degC'
      acc_name='T850'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='t850'; call wrtarr(var_name,acc)

! ACCUMULATED CONVECTIVE CLOUD COVER units = '1'
      acc_name='PMCCLD'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='pcldmc'; call wrtarr(var_name,acc)

! ACCUMULATED CLOUD TOP PRESSURE units = '100 PA'
      acc_name='CLDTPPR'; call getacc(acc_name,acc)
      acc_name='CLDCV'; call getacc(acc_name,acc2)
      where(acc2.gt.0.)
         acc(:,:) = acc(:,:)/acc2(:,:)
      elsewhere
         acc(:,:) = missing
      end where
      var_name='cldtppr'; call wrtarr(var_name,acc)

! ACCUMULATED CLOUD COVER units = '1'
      acc_name='CLDCV'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='pcld'; call wrtarr(var_name,acc)

! 16*P4*(SHA*T4+Z4)*V1*DSIG*DXV units = '100 W*m/s^2'
      acc_name='PEV'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='pev'; call wrtarr(var_name,acc)

! ACCUMULATED NET LONGWAVE FLUX, TOA units = 'W/m^2'
      acc_name='TRNFP0'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='trnfp0'; call wrtarr(var_name,acc)

! ACCUMULATED RADIATION ABS BY SURF units = 'J/m^2'
      acc_name='SRTR'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/dtsec
      var_name='srtr'; call wrtarr(var_name,acc)

! ACCUMLATED NET SURFACE HEATING units = 'J/m^2'
      acc_name='NETH'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/dtsec
      var_name='neth'; call wrtarr(var_name,acc)

! ACCUMULATED NET SHORTWAVE FLUX, TOA units = 'W/m^2'
      acc_name='SRNFP0'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='srnfp0'; call wrtarr(var_name,acc)

! ACCUMULATED INCIDENT SHORTWAVE FLUX, TOA units = 'W/m^2'
      acc_name='SRINCP0'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='srincp0'; call wrtarr(var_name,acc)

! ACCUMULATED NET SHORTWAVE FLUX, SURF units = 'W/m^2'
      acc_name='SRNFG'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='srnfg'; call wrtarr(var_name,acc)

! ACCUMULATED INCIDENT SHORTWAVE FLUX, SURF units = 'W/m^2'
      acc_name='SRINCG'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='srincg'; call wrtarr(var_name,acc)

! ACCUMULATED GROUND TEMPERATURE units = 'degC'
      acc_name='TG1'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='tgrnd'; call wrtarr(var_name,acc)

! POICE+PLICE+(IF SNOW)PEARTH units = '1'
      acc_name='RSIT'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='rsit'; call wrtarr(var_name,acc)

! DIURNAL DELTA TS  OVER SOIL units = 'K'
      acc_name='TDSL'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='tdsl'; call wrtarr(var_name,acc)

! DTHETA/DPHI IN TROPOSPHERE units = 'K s^2/m^2'
      acc_name='DTDP'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='dtdp'; call wrtarr(var_name,acc)

! RUN1 OVER EARTH units = 'mm/day'
      acc_name='RUNE'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/dtday
      var_name='rune'; call wrtarr(var_name,acc)

! RUN1 OVER LAND ICE units = 'mm/day'
      acc_name='RUNLI'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/dtday
      var_name='runli'; call wrtarr(var_name,acc)

! ACCUMULATED SURFACE WIND SPEED units = 'm/s'
      acc_name='WS'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='wspeed'; call wrtarr(var_name,acc)

! ACCUMULATED SURFACE AIR TEMPERATURE units = 'degC'
      acc_name='TS'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='tsurf'; call wrtarr(var_name,acc)

! ACCUMULATED EASTWARD COMP. OF SURFACE WIND units = 'm/s'
      acc_name='US'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='usurf'; call wrtarr(var_name,acc)

! ACCUMULATED NORTHWARD COMP. OF SURFACE WIND units = 'm/s'
      acc_name='VS'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='vsurf'; call wrtarr(var_name,acc)

! PSL USING TS units = 'mb-1000'
      acc_name='SLP'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='slpres'; call wrtarr(var_name,acc)

! UJET units = 'm/s'
      acc_name='UJET'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='ujet'; call wrtarr(var_name,acc)

! VJET units = 'm/s'
      acc_name='VJET'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='vjet'; call wrtarr(var_name,acc)

! ACCUMULATED LOW CLOUD COVER units = '1'
      acc_name='PCLDL'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='pcldl'; call wrtarr(var_name,acc)

! ACCUMULATED MIDDLE CLOUD COVER units = '1'
      acc_name='PCLDM'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='pcldm'; call wrtarr(var_name,acc)

! ACCUMULATED HIGH CLOUD COVER units = '1'
      acc_name='PCLDH'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='pcldh'; call wrtarr(var_name,acc)

! ACCUMULATED BRIGHTNESS TEMPERATURE IN WINDOW units = 'degC'
      acc_name='BTMPW'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='btmpw'; call wrtarr(var_name,acc)

! ACCUMULATED INC. SHORTWAVE x VISIBLE ALBEDO units = 'W/m^2'
      acc_name='SRREF'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='srref'; call wrtarr(var_name,acc)

! TGO2=ODATA(4) units = 'degC'
      acc_name='TOC2'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='toc2'; call wrtarr(var_name,acc)

! ACCUMULATED MAG. OF SURFACE STRESS units = 'kg/m*s^2'
      acc_name='TAUS'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='taumag'; call wrtarr(var_name,acc)

! ACCUMULATED EASTWARD COMP. OF SURFACE STRESS units = 'kg/m*s^2'
      acc_name='TAUUS'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='tauus'; call wrtarr(var_name,acc)

! ACCUMULATED NORTHWARD COMP. OF SURFACE STRESS units = 'kg/m*s^2'
      acc_name='TAUVS'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='tauvs'; call wrtarr(var_name,acc)

! WATER1+WATER2+ICE1+ICE2 units = 'kg/m^2'
      acc_name='GWTR'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='gwtr'; call wrtarr(var_name,acc)

! ACCUMULATED SURFACE AIR SPECIFIC HUMIDITY units = '1'
      acc_name='QS'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='qsurf'; call wrtarr(var_name,acc)

! MAX(0,33-1.8*DAILY MEAN ON TS IN C) units = 'Fahrenheit*day'
      acc_name='STRNGTS'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='strngts'; call wrtarr(var_name,acc)

! UNDERGROUND RUNOFF units = 'mm/day'
      acc_name='ARUNU'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/dtday
      var_name='arunu'; call wrtarr(var_name,acc)

! 18*(DEL(TG)/DEL(TS)-1) units = 'unknown'
      acc_name='DTGDTS'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='dtgdts'; call wrtarr(var_name,acc)

! 8*P*U*Q (VERTICALLY INTEGRATED) units = '12.5 PA*m/s'
      acc_name='PUQ'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='puq'; call wrtarr(var_name,acc)

! 8*P*V*Q (VERTICALLY INTEGRATED) units = '12.5 PA*m/s'
      acc_name='PVQ'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='pvq'; call wrtarr(var_name,acc)

! TGO=ODATA(1) units = 'degC'
      acc_name='TGO'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='tgo'; call wrtarr(var_name,acc)

! ACE2OI=ODATA(3)*POICE units = 'kg/m^2'
      acc_name='MSI2'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='msi2'; call wrtarr(var_name,acc)

! WIND SPEED IN TOP LAYER units = 'm/s'
      acc_name='WLM'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='wlm'; call wrtarr(var_name,acc)

! TGO12=ODATA(5) units = 'degC'
      acc_name='TGO2'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='tgo2'; call wrtarr(var_name,acc)

! EVAP*POCEAN units = 'mm/day'
      acc_name='EVAPO'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/dtday
      var_name='evapo'; call wrtarr(var_name,acc)

! EVAP*POICE units = 'mm/day'
      acc_name='EVAPI'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/dtday
      var_name='evapi'; call wrtarr(var_name,acc)

! EVAP OVER LAND ICE units = 'mm/day'
      acc_name='EVAPLI'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/dtday
      var_name='evapli'; call wrtarr(var_name,acc)

! EVAP OVER EARTH units = 'mm/day'
      acc_name='EVAPE'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/dtday
      var_name='evape'; call wrtarr(var_name,acc)

! F0DT*POCEAN, NET HEAT AT Z0 units = 'J/m^2'
      acc_name='F0OC'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='f0oc'; call wrtarr(var_name,acc)

! F0DT*POICE, NET HEAT AT Z0 units = 'J/m^2'
      acc_name='F0OI'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='f0oi'; call wrtarr(var_name,acc)

! F0DT, NET HEAT AT Z0 OVER LAND ICE units = 'J/m^2'
      acc_name='F0LI'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='f0li'; call wrtarr(var_name,acc)

!  F0DT, NET HEAT AT Z0 OVER EARTH units = 'J/m^2'
      acc_name='F0E'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='f0e'; call wrtarr(var_name,acc)

! F1DT OVER LAND ICE units = 'J/m^2'
      acc_name='F1LI'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='f1li'; call wrtarr(var_name,acc)

! SNOW FALL (H2O EQUIV) units = 'mm/day'
      acc_name='SNWF'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/dtday
      var_name='snwf'; call wrtarr(var_name,acc)

! SURF AIR TEMP OVER LAND ICE units = 'degC'
      acc_name='TSLI'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='tsli'; call wrtarr(var_name,acc)

! F2DT OVER LAND ICE units = 'J/m^2'
      acc_name='ERUN2'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='erun2'; call wrtarr(var_name,acc)

! SENS HEAT FLUX OVER LAND ICE units = 'W/m^2'
      acc_name='SHDTLI'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='shdtli'; call wrtarr(var_name,acc)

! LATENT HEAT FLUX OVER LAND ICE units = 'W/m^2'
      acc_name='EVHDT'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='evhdtli'; call wrtarr(var_name,acc)

! TRHDT OVER LAND ICE units = 'J/m^2'
      acc_name='TRHDT'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='trhdtli'; call wrtarr(var_name,acc)

! MAX(COMPOSITE TS) units = 'degC'
      acc_name='TMAX'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='tmax'; call wrtarr(var_name,acc)

! MIN(COMPOSITE TS) units = 'degC'
      acc_name='TMIN'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='tmin'; call wrtarr(var_name,acc)

! MIN(DIURNAL MAX OF COMPOSITE TS) units = 'degC'
      acc_name='TMNMX'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='tmnmx'; call wrtarr(var_name,acc)

! POTENTIAL EVAPORATION units = 'mm/day'
      acc_name='PEVAP'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/dtday
      var_name='pevap'; call wrtarr(var_name,acc)

! MAX TS OVER EARTH FOR CURRENT DAY units = 'K'
      acc_name='TMAXE'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='tmaxe'; call wrtarr(var_name,acc)

! LIQUID WATER PATH units = 'kg/m^2'
      acc_name='WMSUM'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='wmsum'; call wrtarr(var_name,acc)

! ACCUMULATED SHALLOW CONVECTIVE CLOUD COVER units = '1'
      acc_name='PSCLD'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='pscld'; call wrtarr(var_name,acc)

! ACCUMULATED DEEP CONVECTIVE CLOUD COVER units = '1'
      acc_name='PDCLD'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='pdcld'; call wrtarr(var_name,acc)

! ACCUMULATED DEEP CONVECTIVE CLOUD OCCURENCE units = '1'
      acc_name='DCNVFRQ'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='dcnvfrq'; call wrtarr(var_name,acc)

! ACCUMULATED SHALLOW CONVECTIVE CLOUD OCCURENCE units = '1'
      acc_name='SCNVFRQ'; call getacc(acc_name,acc)
      acc(:,:) = acc(:,:)/idacc(ia)
      var_name='scnvfrq'; call wrtarr(var_name,acc)

! albedos
      acc_name='SRNFP0'; call getacc(acc_name,acc)
      acc_name='SRINCP0'; call getacc(acc_name,acc2)
      where(acc2.gt.0.)
         acc(:,:) = 1.-acc(:,:)/acc2(:,:)
      elsewhere
         acc(:,:) = missing
      end where
      var_name='palb'; call wrtarr(var_name,acc)

      acc_name='SRNFG'; call getacc(acc_name,acc)
      acc_name='SRINCG'; call getacc(acc_name,acc2)
      where(acc2.gt.0.)
         acc(:,:) = 1.-acc(:,:)/acc2(:,:)
      elsewhere
         acc(:,:) = missing
      end where
      var_name='galb'; call wrtarr(var_name,acc)

      call close_acc
      call close_out

      end program ijdag
