      subroutine def_acc
c-----------------------------------------------------------------------
c define acc names, units, etc
c-----------------------------------------------------------------------
      use DAGCOM
      implicit none
      integer :: k
      call tsf_defs
      call j_defs
      name_reg=name_j
      do k=1,kaj ! to avoid naming conflicts, put a prefix
         name_j(k) = 'J_'//trim(name_j(k))
         name_reg(k) = 'reg_'//trim(name_reg(k))
      enddo
      call pj_defs
      call jl_defs
      call sjl_defs
      call ij_defs
      call il_defs
      call wave_defs
      call jk_defs
      call ijk_defs
      call diurn_defs
      return
      end subroutine def_acc

      subroutine j_defs
!@sum j_defs definitions for j_xx zonal budget diagnostics
!@+   diags are printed out in the order they are defined
!@auth G. Schmidt/M. Kelley
      use CONSTANT, only : grav,sday,shw,rgas,omega,bygrav,gamd
      use MODEL_COM, only : jm,lm,ls1,dtsrc,fim,sige,kocean
      use DAGCOM
      implicit none
      integer :: k,kk
c
      do k=1,kaj
         write(name_j(k),'(a2,i3.3)') 'AJ',k
         lname_j(k) = 'unused'
         units_j(k) = 'unused'
         stitle_j(k)= 'no output'
         scale_j(k) = 1.
         ia_j(k)    = 1.
      enddo
c
      k=0
c
      k=k+1
      J_SRINCP0= k ! SRINCP0 (W/M**2)                              2 RD
      name_j(k) = 'inc_sw'
      lname_j(k) = 'SOLAR RADIATION INCIDENT ON PLANET'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' INC SW (W/m^2)'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_SRNFP0=  k ! SRNFP0 (W/m**2)                               2 RD
      name_j(k) = 'sw_abs_p0'
      lname_j(k) = 'SOLAR RADIATION ABSORBED BY PLANET'
      units_j(k) = 'W/m^2'
      stitle_j(k)= '0SW ABS BELOW P0'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_SRNFP1=  k ! SRNFP1 (W/m**2)                               2 RD
      name_j(k) = 'sw_abs_p1'
      lname_j(k) = 'SOLAR RADIATION ABSORBED BELOW PTOP'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' SW ABS BELOW P1'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_SRABS =  k ! SRABSATM=AJ(SRNFP0)-AJ(SRNFG) (W/m**2)        2 D1
      name_j(k) = 'sw_abs_atm'
      lname_j(k) = 'SOLAR RADIATION ABSORBED BY ATMOSPHERE'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' SW ABS BY ATMOS'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_SRINCG=  k ! SRINCG (W/m**2)                               2 RD
      name_j(k) = 'sw_inc_z0'
      lname_j(k) = 'SOLAR RADIATION INCIDENT ON GROUND'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' SW INC ON Z0   '
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_SRNFG =  k ! SRNFG (W/m**2)                                2 RD
      name_j(k) = 'sw_abs_z0'
      lname_j(k) = 'SOLAR RADIATION ABSORBED BY GROUND'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' SW ABS AT Z0   '
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_TRNFP0=  k ! TRNFP0=AJ(ALBVIS)+A2BYA1*AJ(TRHDT)/DTS(W/m**2)2 D1
      name_j(k) = 'net_lw_p0'
      lname_j(k) = 'THERMAL RADIATION EMITTED BY PLANET'
      units_j(k) = 'W/m^2'
      stitle_j(k)= '0NET LW AT P0   '
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_TRNFP1=  k ! TRNFP1=AJ(ALBNIR)+A2BYA1*AJ(TRHDT)/DTS(W/m**2)2 D1
      name_j(k) = 'net_lw_p1'
      lname_j(k) = 'NET THERMAL RADIATION AT PTOP'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' NET LW AT P1   '
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_TRHDT =  k ! TRHDT (J/m**2)                                1 SF
      name_j(k) = 'net_lw_z0'
      lname_j(k) = 'NET THERMAL RADIATION AT GROUND'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' NET LW AT Z0   '
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c
      k=k+1
      J_TRINCG= k ! TRINCG (W/m**2)                               2 RD
      name_j(k) = 'lw_inc_z0'
      lname_j(k) = 'THERMAL RADIATION INCIDENT ON GROUND'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' LW INC ON Z0   '
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_BRTEMP= k ! BTEMPW-TF                                     2 RD
      name_j(k) = 'btemp_window'
      lname_j(k) = 'BRIGHTNESS TEMP THROUGH WINDOW REGION'
      units_j(k) = 'C'
      stitle_j(k)= ' LW WINDOW BTEMP'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_RNFP0 = k ! RNFP0=AJ(SRNFP0)+AJ(TRNFP0) (W/m**2)          2 D1
      name_j(k) = 'net_rad_p0'
      lname_j(k) = 'NET RADIATION OF PLANET'
      units_j(k) = 'W/m^2'
      stitle_j(k)= '0NET RAD AT P0  '
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_CLRTOA= k ! clear sky radiative forcing top of atmosphere
      name_j(k) = 'net_clr_toa'
      lname_j(k) = 'NET CLEAR SKY RADIATION AT P0'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' NET CLR RAD  P0'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_RNFP1 = k ! RNFP1=AJ(SRNFP1)+AJ(TRNFP1) (W/m**2)          2 D1
      name_j(k) = 'net_rad_p1'
      lname_j(k) = 'NET RADIATION BELOW PTOP'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' NET RAD AT P1  '
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_CLRTRP= k ! clear sky radiative forcing tropopause
      name_j(k) = 'net_clr_trp'
      lname_j(k) = 'NET CLEAR SKY RADIATION AT TROPOPAUSE (WMO)'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' NET CLR RAD TRP'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_TOTTRP= k ! total radiative forcing tropopause
      name_j(k) = 'net_tot_trp'
      lname_j(k) = 'NET RADIATION AT TROPOPAUSE (WMO)'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' NET RAD (TROPP)'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_RHDT  = k ! RHDT=A1BYA2*AJ(SRNFG)*DTS+AJ(TRHDT)(J/m^2)    1 D1
      name_j(k) = 'net_rad_z0'
      lname_j(k) = 'NET RADIATION ABSORBED BY GROUND'
      units_j(k) = 'J/m^2'
      stitle_j(k)= ' NET RAD AT Z0  '
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c
      k=k+1
      J_SHDT  = k ! SHEATDT (J/m**2)                              1 SF
      name_j(k) = 'snsht_flx'
      lname_j(k) = 'SENSIBLE HEAT FLUX INTO THE GROUND'
      units_j(k) = 'W/m^2'
      stitle_j(k)= '0SENSBL HEAT FLX'
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c
      k=k+1
      J_EVHDT = k ! EVHDT (J/m**2)                                1 SF
      name_j(k) = 'evht_flx'
      lname_j(k) = 'LATENT HEAT FLUX INTO THE GROUND'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' EVAPOR HEAT FLX'
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c
      k=k+1
      J_EPRCP = k ! ENERGP (J/m**2)                               1 CN
      name_j(k) = 'prec_ht_flx'
      lname_j(k) = 'PRECIPITATION HEAT FLUX INTO THE GROUND'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' PRECIP HEAT FLX'
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c
      k=k+1
      J_HZ0   = k ! HEATZ0=RHDT+SHDT+EVHDT+EPRCP (J/m**2)   1 D1
      name_j(k) = 'nt_ht_z0'
      lname_j(k) = 'NET HEATING AT GROUND SURFACE'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' NET HEAT AT Z0 '
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c
      k=k+1
      J_ERVR = k ! Energy in river discharge                     1 GP
      name_j(k) = 'ht_rvr_disch'
      lname_j(k) = 'HEAT IN RIVER DISCHARGE'
      units_j(k) = 'W/m^2'
      stitle_j(k)= '0HT RVR DISCH   '
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c
      k=k+1
      J_HZ1   = k ! Net heating at ocean surface                  1 D1
      name_j(k) = 'net_ht_z1'
      lname_j(k) = 'NET HEAT AT SURFACE (INCL RIVERS)'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' NET HEAT Z0+RVR'
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c
      k=k+1
      J_ERUN = k ! ERUNOFF (J/m**2)                                1 GP
      name_j(k) = 'ht_runoff'
      lname_j(k) = 'HEAT RUNOFF'
      units_j(k) = 'W/m^2'
      stitle_j(k)= '0HEAT RUNOFF '
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c
      k=k+1
      J_HMELT = k ! net amount of energy associated with ice melt/form
      name_j(k) = 'ht_ice_melt'
      lname_j(k) = 'NET HEAT OF ICE MELT/FORMATION'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' HT ICE MLT/FORM'
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c
C**** Note this is used for ice in fixed SST runs, but for ocean in
C**** qflux runs. Over land, it is always used for landice changes.
      k=k+1
      J_IMPLH = k               !                                 1 GP
      name_j(k) = 'impl_ht'
      lname_j(k) = 'DOWNWARD IMPLICIT HEAT FLUX AT ICE BASE/OCN ML'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' DWN IMPL HT FLX'
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c
      k=k+1
      J_HZ2   = k !
      name_j(k) = 'net_ht_hz2'
      lname_j(k) = 'NET HEAT CONVERGENCE INTO GROUND TYPE'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' NET HT CNV GRND'
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c
      k=k+1
      J_PRCPSS= k ! PRCPSS (100 PA)                               1 CN
      name_j(k) = 'ssprec'
      lname_j(k) = 'SUPER SATURATION PRECIPITATION'
      units_j(k) = 'mm/day'
      stitle_j(k)= '1SS PRECIP(MM/D)'
      scale_j(k) = 100.*SDAY/(DTsrc*GRAV)
      ia_j(k) = ia_src
c
      k=k+1
      J_PRCPMC= k ! PRCPMC (100 PA)                               1 CN
      name_j(k) = 'mcprec'
      lname_j(k) = 'MOIST CONVECTIVE PRECIPITATION'
      units_j(k) = 'mm/day'
      stitle_j(k)= ' MC PRECIP(MM/D)'
      scale_j(k) = 100.*SDAY/(DTsrc*GRAV)
      ia_j(k) = ia_src
c
      k=k+1
      J_PRCP  = k ! PRCP=AJ(PRCPSS)+AJ(PRCPMC) (100 PA)           1 D1
      name_j(k) = 'prec'
      lname_j(k) = 'PRECIPITATION'
      units_j(k) = 'mm/day'
      stitle_j(k)= ' PRECIP (MM/DAY)'
      scale_j(k) = 100.*SDAY/(DTsrc*GRAV)
      ia_j(k) = ia_src
c
      k=k+1
      J_EVAP  = k ! EVAP (KG/m**2)                                1 GD
      name_j(k) = 'evap'
      lname_j(k) = 'EVAPORATION'
      units_j(k) = 'mm/day'
      stitle_j(k)= ' EVAPOR (MM/DAY)'
      scale_j(k) = SDAY/DTSRC
      ia_j(k) = ia_src
c
      k=k+1
      J_RUN  = k ! RUNOFF (KG/m**2)                                1 GP
      name_j(k) = 'wat_runoff'
      lname_j(k) = 'WATER RUNOFF AT GROUND SURFACE'
      units_j(k) = 'mm/day'
      stitle_j(k)= '0WATER RUNOFF'
      scale_j(k) = SDAY/DTSRC
      ia_j(k) = ia_src
c
      k=k+1
      J_RVRD  = k ! RIVER DISCHARGE                           1 GP
      name_j(k) = 'river_discharge'
      lname_j(k) = 'RIVER DISCHARGE'
      units_j(k) = 'mm/day'
      stitle_j(k)= ' RVR DISCH(MM/D)'
      scale_j(k) = SDAY/DTSRC
      ia_j(k) = ia_src
c
      k=k+1
      J_IMELT = k ! net amount ice melt/formation in oc/lk  1 GP
      name_j(k) = 'ice_melt'
      lname_j(k) = 'NET ICE MELTING/FORMATION'
      units_j(k) = 'mm/day'
      stitle_j(k)= ' ICE MELT/FORM  '
      scale_j(k) = SDAY/DTSRC
      ia_j(k) = ia_src
c
C**** Note this is used for ice in fixed SST runs, but for ocean in
C**** qflux runs. Over land, it is always used for landice changes.
      k=k+1
      J_IMPLM = k               !                                 1 GP
      name_j(k) = 'impl_m_flux'
      lname_j(k) = 'DOWNWARD IMPLICIT MASS FLUX AT ICE BASE/OCN ML'
      units_j(k) = 'mm/day'
      stitle_j(k)= ' DWN IMPL WT FLX'
      scale_j(k) = SDAY/DTSRC
      ia_j(k) = ia_src
c
      k=k+1
      J_H2OCH4 = k               !                                 1 GP
      name_j(k) = 'h2o_from_ch4'
      lname_j(k) = 'WATER DERIVED FROM CH4 OXIDATION IN STRATOSPHERE'
      units_j(k) = '10^6 mm/day'
      stitle_j(k)= ' H2O BY CH4(x1M)'
      scale_j(k) = 2d6
      ia_j(k) = ia_12hr
c
      k=k+1
      J_SMELT= k ! salt flux associated with ice melt/formation   1 GD
      name_j(k) = 's_ice_melt'
      lname_j(k) = 'SALT IN ICE MELT/FORMATION'
      units_j(k) = '10^-3 kg/m^2/day'
      stitle_j(k)= '0SALT MELT (x1K)'
      scale_j(k) = 1000.*SDAY/DTSRC
      ia_j(k) = ia_src
c
      k=k+1
      J_TG1   = k ! TG1 (K-TF)                                    1 GD
      name_j(k) = 'tg1'
      lname_j(k) = 'TEMPERATURE OF GROUND LAYER 1'
      units_j(k) = '.1 C'
      stitle_j(k)= '0TG1 (.1 C)     '
      scale_j(k) = 10.
      ia_j(k) = ia_srf
c
      k=k+1
      J_TG2   = k ! TG2 (K-TF)                                    1 GD
      name_j(k) = 'tg2'
      lname_j(k) = 'TEMPERATURE OF GROUND LAYER 2'
      units_j(k) = '.1 C'
      stitle_j(k)= ' TG2 (.1 C)     '
      scale_j(k) = 10.
      ia_j(k) = ia_srf
c
      k=k+1
      J_TSRF  = k ! TS (K-TF)                                     3 SF
      name_j(k) = 'tsurf'
      lname_j(k) = 'SURFACE AIR TEMPERATURE'
      units_j(k) = '.1 C'
      stitle_j(k)= ' T SURF (.1 C)  '
      scale_j(k) = 10.
      ia_j(k) = ia_srf
c
      k=k+1
      J_TX    = k ! TX (K-TF)  (INTEGRAL OVER ATMOSPHERE OF)      4 DA
      name_j(k) = 'tair'
      lname_j(k) = 'AIR TEMPERATURE'
      units_j(k) = '.1 C'
      stitle_j(k)= ' T AIR (.1 C)   '
      scale_j(k) = 10.
      ia_j(k) = ia_dga
c
      k=k+1
      J_TX1   = k ! TX1 (K-TF)                                    4 DA
      name_j(k) = 't1'
      lname_j(k) = 'TEMPERATURE OF AIR LAYER 1'
      units_j(k) = '.1 C'
      stitle_j(k)= ' T1 (.1 C)      '
      scale_j(k) = 10.
      ia_j(k) = ia_dga
c
      k=k+1
      J_DTDJS = k ! T(J+1)-T(J-1)  (SUM OVER STRATOSPHERE OF)     4 DA
      name_j(k) = 'dtdlat_strat'
      lname_j(k) = 'STRATO TEMP CHANGE PER DEGREE LATITUDE'
      units_j(k) = 'deg C/deg lat'
      stitle_j(k)= '0DT/DLAT(STRAT) '
      scale_j(k) = .5D2*(JM-1.)/((SIGE(LS1)-SIGE(LSTR+1)+1d-12)*180.)
      ia_j(k) = ia_dga
c
      k=k+1
      J_DTDJT = k ! T(J+1)-T(J-1)  (SUM OVER TROPOSPHERE OF)      4 DA
      name_j(k) = 'dtdlat_trop'
      lname_j(k) = 'TROPO TEMP CHANGE PER DEGREE LATITUDE'
      units_j(k) = 'deg C/deg lat'
      stitle_j(k)= ' DT/DLAT(TROPO) '
      scale_j(k) = .5d2*(JM-1.)/((SIGE(1)-SIGE(LS1))*180.)
      ia_j(k) = ia_dga
c
      k=k+1
      J_DTSGST= k ! DTH/DPHI  (STRATOSPHERE)                      4 DA
      name_j(k) = 'sstab_strat'
      lname_j(k) = 'STRATOSPHERIC STATIC STABILITY'
      units_j(k) = 'C/km'
      stitle_j(k)= '0STAT STB(STRAT)'
      scale_j(k) = 1.D3*GRAV*P1000K
      ia_j(k) = ia_dga
c
      k=k+1
      J_DTDGTR= k ! DTH/DPHI  (TROPOSPHERE)                       4 DA
      name_j(k) = 'sstab_trop'
      lname_j(k) = 'TROPOSPHERIC STATIC STABILITY'
      units_j(k) = 'C/km'
      stitle_j(k)= ' STAT STB(TROPO)'
      scale_j(k) = 1.D3*GRAV*P1000K
      ia_j(k) = ia_dga
c
      k=k+1
      J_RICST = k ! .0625*DTH*DLNP/(DU*DU+DV*DV)  (STRATOSPHERE)  4 DA
      name_j(k) = 'rich_num_strat'
      lname_j(k) = 'STRATOSPHERIC RICHARDSON NUMBER'
      units_j(k) = '1'
      stitle_j(k)= '0RICH NUM(STRAT)'
      scale_j(k) = 16.*RGAS
      ia_j(k) = ia_dga
c
      k=k+1
      J_RICTR = k ! .0625*DTH*DLNP/(DU*DU+DV*DV)  (TROPOSPHERE)   4 DA
      name_j(k) = 'rich_num_trop'
      lname_j(k) = 'TROPOSPHERIC RICHARDSON NUMBER'
      units_j(k) = '1'
      stitle_j(k)= ' RICH NUM(TROPO)'
      scale_j(k) = 16.*RGAS
      ia_j(k) = ia_dga
c
      k=k+1
      J_ROSST = k ! 4*UMAX/(DX*SINJ)  (STRATOSPHERE)              4 DA
      name_j(k) = 'ross_num_strat'
      lname_j(k) = 'STRATOSPHERIC ROSSBY NUMBER'
      units_j(k) = '1'
      stitle_j(k)= ' ROSS NUM(STRAT)'
      scale_j(k) = .5/(2.*OMEGA*FIM)
      ia_j(k) = ia_dga
c
      k=k+1
      J_ROSTR = k ! 4*UMAX/(DX*SINJ)  (TROPOSPHERE)               4 DA
      name_j(k) = 'ross_num_trop'
      lname_j(k) = 'TROPOSPHERIC ROSSBY NUMBER'
      units_j(k) = '1'
      stitle_j(k)= ' ROSS NUM(TROPO)'
      scale_j(k) = .5/(2.*OMEGA*FIM)
      ia_j(k) = ia_dga
c
      k=k+1
      J_LSTR  = k ! SQRT(DTH/DLNP)/SINJ  (STRATOSPHERE)           4 DA
      name_j(k) = 'ross_radius_strat'
      lname_j(k) = 'ROSSBY RADIUS IN THE STRATOSPHERE'
      units_j(k) = '10**5 m'
      stitle_j(k)= ' L(STRAT)(10**5)'
      scale_j(k) = 1d-5*SQRT(RGAS)/(2.*OMEGA)
      ia_j(k) = ia_dga
c
      k=k+1
      J_LTRO  = k ! SQRT(DTH/DLNP)/SINJ  (TROPOSPHERE)            4 DA
      name_j(k) = 'ross_radius_trop'
      lname_j(k) = 'ROSSBY RADIUS IN THE TROPOSPHERE'
      units_j(k) = '10**5 m'
      stitle_j(k)=  ' L(TROP) (10**5)'
      scale_j(k) = 1d-5*SQRT(RGAS)/(2.*OMEGA)
      ia_j(k) = ia_dga
c
      k=k+1
      J_GAM   = k ! GAM  (K/m)  (*SIG(TROPOSPHERE)/GRAV)          4 DA
      name_j(k) = 'lapse_rate'
      lname_j(k) = 'MEAN LAPSE RATE'
      units_j(k) = 'K/km'
      stitle_j(k)= '0GAM(K/KM)      '
      scale_j(k) = 1d3*GRAV
      ia_j(k) = ia_dga
c
      k=k+1
      J_GAMM  = k ! GAMM  (K-S**2/m**2)  (SIG(TROPOSPHERE)/GAMD)  4 DA
      name_j(k) = 'lapse_rate_m'
      lname_j(k) = 'MOIST ADIABATIC LAPSE RATE'
      units_j(k) = 'K/km'
      stitle_j(k)= ' GAMM(K/KM)     '
      scale_j(k) = 1.D3*GAMD/(SIGE(1)-SIGE(LS1))
      ia_j(k) = ia_dga
c
      k=k+1
      J_GAMC  = k ! GAMC  (K/m)                                   4 DA
      name_j(k) = 'lapse_rate_c'
      lname_j(k) = 'GAMC'
      units_j(k) = 'K/Km'
      stitle_j(k)= ' GAMC(K/KM)     '
      scale_j(k) = 1d3
      ia_j(k) = ia_dga
c
      k=k+1
      J_PCLDSS= k ! PCLDSS (1)  (COMPOSITE OVER ATMOSPHERE)       2 RD
      name_j(k) = 'sscld'
      lname_j(k) = 'SUPER SATURATION CLOUD COVER'
      units_j(k) = '%'
      stitle_j(k)= '0TOT SUP SAT CLD'
      scale_j(k) = 100.
      ia_j(k) = ia_rad
c
      k=k+1
      J_PCLDMC= k ! PCLDMC (1)  (COMPOSITE OVER ATMOSPHERE)       2 RD
      name_j(k) = 'mccld'
      lname_j(k) = 'MOIST CONVECTIVE CLOUD COVER'
      units_j(k) = '%'
      stitle_j(k)= ' TOT MST CNV CLD'
      scale_j(k) = 100.
      ia_j(k) = ia_rad
c
      k=k+1
      J_PCLD  = k ! PCLD (1)  (COMPOSITE OVER ATMOSPHERE)         2 RD
      name_j(k) = 'totcld'
      lname_j(k) = 'TOTAL CLOUD COVER'
      units_j(k) = '%'
      stitle_j(k)= ' TOTAL CLD COVER'
      scale_j(k) = 100.
      ia_j(k) = ia_rad
c
      k=k+1
      J_QP    = k ! Q*P (100 PA)  (INTEGRAL OVER ATMOSPHERE OF)   4 DA
      name_j(k) = 'atmh2o'
      lname_j(k) = 'WATER CONTENT OF ATMOSPHERE'
      units_j(k) = 'mm'
      stitle_j(k)= ' H2O OF ATM (MM)'
      scale_j(k) = 100.*BYGRAV
      ia_j(k) = ia_dga
c
      k=k+1
      J_WTR1  = k ! WTR1 (KG/m**2)                                1 GD
      name_j(k) = 'wat_g1'
      lname_j(k) = 'WATER IN GROUND LAYER 1'
      units_j(k) = 'kg/m^2'
      stitle_j(k)= '0WATER IN G1    '
      scale_j(k) = 1.
      ia_j(k) = ia_src
c
      k=k+1
      J_ACE1  = k ! ACE1 (KG/m**2)                                1 GD
      name_j(k) = 'ice_g1'
      lname_j(k) = 'ICE IN GROUND LAYER 1'
      units_j(k) = 'kg/m^2'
      stitle_j(k)= ' ICE IN G1      '
      scale_j(k) = 1.
      ia_j(k) = ia_src
c
      k=k+1
      J_WTR2  = k ! WTR2 (KG/m**2)                                1 GD
      name_j(k) = 'wat_g2'
      lname_j(k) = 'WATER IN GROUND LAYER 2'
      units_j(k) = '10^2 kg/m^2'
      stitle_j(k)= ' WATER G2 x0.01 '
      scale_j(k) = 1d-2
      ia_j(k) = ia_src
c
      k=k+1
      J_ACE2  = k ! ACE2 (KG/m**2)                                1 GD
      name_j(k) = 'ice_g2'
      lname_j(k) = 'ICE IN GROUND LAYER 2'
      units_j(k) = '10^2 kg/m^2'
      stitle_j(k)= ' ICE G2   x0.01 '
      scale_j(k) = 1d-2
      ia_j(k) = ia_src
c
      k=k+1
      J_SNOW  = k ! SNOW (KG/m**2)                                1 GD
      name_j(k) = 'snowdp'
      lname_j(k) = 'SNOW DEPTH'
      units_j(k) = 'kg/m^2'
      stitle_j(k)= ' SNOW DEPTH     '
      scale_j(k) = 1.
      ia_j(k) = ia_src
c
      k=k+1
      J_RSNOW = k ! PSNOW (1)                                     4 DA
      name_j(k) = 'snow_cover'
      lname_j(k) = 'SNOW COVER'
      units_j(k) = '%'
      stitle_j(k)= ' SNOW COVER     '
      scale_j(k) = 100.
      ia_j(k) = ia_src
c
      k=k+1
      J_RSI   = k ! RSI (1)                                       1 GD
      name_j(k) = 'ocn_lak_ice_frac'
      lname_j(k) = 'OCEAN/LAKE ICE COVER'
      units_j(k) = '%'
      stitle_j(k)= ' OC/LK ICE COVER'
      scale_j(k) = 100.
      ia_j(k) = ia_src
c
      IF (KOCEAN.eq.1) THEN ! only for non-fixed SST runs
      k=k+1
      J_OHT   = k ! OCEAN TRANSPORT                               1 GD
      name_j(k) = 'ocn_ht_trans'
      lname_j(k) = 'CONVERGED OCEAN HEAT TRANSPORT'
      units_j(k) = 'W/m^2'
      stitle_j(k)= '0OCEAN TRNS CONV'
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c
      k=k+1
      J_FTHERM= k ! ENERGY DIFFUSION INTO THERMOCLINE (W/m**2) .5*9 MN
      name_j(k) = 'ht_thermocline'
      lname_j(k) = 'ENERGY DIFFUSION INTO THE THERMOCLINE'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' HT INTO THRMOCL'
      scale_j(k) = 2d3*SHW/SDAY
      ia_j(k) = ia_12hr
      END IF
c
      k=k+1
      J_TYPE  = k ! PTYPE                                         1 GD
      name_j(k) = 'surf_type_frac'
      lname_j(k) = 'SURF TYPE FRACT'
      units_j(k) = '%'
      stitle_j(k)= ' SURF TYPE FRACT'
      scale_j(k) = 100.
      ia_j(k) = ia_srf
c set the number of directly output budget diagnostics
      k_j_out=k
c None of the following will be printed out:
      k=k+1
      J_HSURF = k ! TRNFP0-TRNFG (W/m**2)                         2 RD
      name_j(k) = 'HSURF'
      lname_j(k) = 'SURFACE THERMAL HEATING'
      units_j(k) = 'W/m^2'
      stitle_j(k)= 'no output'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_HATM  = k ! TRNFP1-TRNFG (W/m**2)                         2 RD
      name_j(k) = 'HATM'
      lname_j(k) = 'ATMOSPHERIC THERMAL HEATING'
      units_j(k) = 'W/m^2'
      stitle_j(k)= 'no output'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
C**** Note: next eight diagnostics must remain in order
      J_PLAVIS= k ! PLAVIS*S0*COSZ (W/m**2)                       2 RD
      name_j(k) = 'plan_refl_vis'
      lname_j(k) = 'PLANETARY REFLECTED RAD. IN VISUAL'
      units_j(k) = 'W/m^2'
      stitle_j(k)= 'no output'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_PLANIR= k ! PLANIR*S0*COSZ (W/m**2)                       2 RD
      name_j(k) = 'plan_refl_nir'
      lname_j(k) = 'PLANETARY REFLECTED RAD. IN NEAR IR'
      units_j(k) = 'W/m^2'
      stitle_j(k)= 'no output'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_ALBVIS= k ! ALBVIS*S0*COSZ (W/m**2)                       2 RD
      name_j(k) = 'surf_refl_vis'
      lname_j(k) = 'GROUND REFLECTED RAD. IN VISUAL'
      units_j(k) = 'W/m^2'
      stitle_j(k)= 'no output'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_ALBNIR= k ! ALBNIR*S0*COSZ (W/m**2)                       2 RD
      name_j(k) = 'surf_refl_nir'
      lname_j(k) = 'GROUND REFLECTED RAD. IN NEAR IR'
      units_j(k) = 'W/m^2'
      stitle_j(k)= 'no output'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_SRRVIS= k ! SRRVIS*S0*COSZ (W/m**2)                       2 RD
      name_j(k) = 'atm_refl_vis'
      lname_j(k) = 'ATMOSPHERIC REFLECTED RAD. IN VISUAL'
      units_j(k) = 'W/m^2'
      stitle_j(k)= 'no output'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_SRRNIR= k ! SRRNIR*S0*COSZ (W/m**2)                       2 RD
      name_j(k) = 'atm_refl_nir'
      lname_j(k) = 'ATMOSPHERIC REFLECTED RAD. IN NEAR IR'
      units_j(k) = 'W/m^2'
      stitle_j(k)= 'no output'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_SRAVIS= k ! SRAVIS*S0*COSZ (W/m**2)                       2 RD
      name_j(k) = 'atm_abs_vis'
      lname_j(k) = 'ATMOSPHERIC ABSORPTION IN VISUAL'
      units_j(k) = 'W/m^2'
      stitle_j(k)= 'no output'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_SRANIR= k ! SRANIR*S0*COSZ (W/m**2)                       2 RD
      name_j(k) = 'atm_abs_nir'
      lname_j(k) = 'ATMOSPHERIC ABSORPTION IN NEAR IR'
      units_j(k) = 'W/m^2'
      stitle_j(k)= 'no output'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_CLDDEP= k ! PBOTMC-PTOPMC (100 PA)                        2 RD
      name_j(k) = 'CLDDEP'
      lname_j(k) = 'MOIST CONVECTIVE CLOUD DEPTH'
      units_j(k) = 'mb'
      stitle_j(k)= 'no output'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      if (k .gt. kaj) then
        write (6,*) 'j_defs: Increase kaj=',kaj,' to at least ',k
        call stop_model( 'kaj too small', 255 )
      end if
      write (6,*) 'Number of AJ diagnostics defined: kajmax=',k
      do kk=1,k
        write (6,'(i4,'':'',a)') kk,trim(lname_j(kk))
      end do

      return
      end subroutine j_defs

      subroutine ij_defs
      use constant
      use MODEL_COM
      use DAGCOM
      implicit none
      integer :: k,kk
c
      do k=1,kaij
         write(name_ij(k),'(a3,i3.3)') 'AIJ',k
         lname_ij(k) = 'unused'
         units_ij(k) = 'unused'
         ia_ij(k) = ia_src
         scale_ij(k) = 1.
         jgrid_ij(k) = 1
         iw_ij(k) = iw_all
         ir_ij(k) = ir_pct
      enddo
c
      k=0
C**** AIJ diagnostic names:
C**** NAME     NO.    DESCRIPTION   (SCALE)*IDACC  LOCATION
C**** ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c
      k=k+1 !
      IJ_RSOI = k ! POICE (1)            1 GD
      lname_ij(k) = 'OCEAN/LAKE ICE COVERAGE'
      units_ij(k) = '%'
      name_ij(k) = 'RSOI'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_RSNW = k ! PSNOW (1)            1 GD
      lname_ij(k) = 'SNOW COVERAGE'
      units_ij(k) = '%'
      name_ij(k) = 'RSNW'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_SNOW = k ! SNOW (KG/m**2)       1 GD
      lname_ij(k) = 'SNOW DEPTH'    ! 'SNOW MASS'
      units_ij(k) = 'mm H2O'
      name_ij(k) = 'SNOW'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
cc    iw_ij(k) = iw_firm
      ir_ij(k) = ir_0_26_150
c
      k=k+1 !
      IJ_SHDT = k ! SHDT (J/m**2)        1 SF
      lname_ij(k) = 'SENSIBLE HEAT FLUX'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'SHDT'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      ir_ij(k) = ir_m265_95
c
      k=k+1 !
      IJ_PREC = k ! PREC (KG/m**2)       1 CN
      lname_ij(k) = 'PRECIPITATION'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'PREC'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
      ir_ij(k) = ir_0_3_15
c
      k=k+1 !
      IJ_EVAP = k ! EVAP (KG/m**2)       1 SF
      lname_ij(k) = 'EVAPORATION'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'EVAP'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
      ir_ij(k) = ir_0_3_15
c
      k=k+1 !
      IJ_BETA = k ! BETA (1)             1 GD
      lname_ij(k) = 'GROUND WETNESS (VEG ROOTS)'
      units_ij(k) = '%'
      name_ij(k) = 'BETA'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
      iw_ij(k) = iw_veg
c
      k=k+1 !
      IJ_PRES = k ! PIJ (100 PA)  (NO PRINTOUT)  4 DA
      lname_ij(k) = 'SURFACE PRESSURE  - PTOP'
      units_ij(k) = 'mb'
      name_ij(k) = 'PRES'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_1775
c
      k=k+1 !
      IJ_PHI1K = k ! PHI1000 (M**2/S**2) 4 DA
      lname_ij(k) = '1000mb HEIGHT'
      units_ij(k) = 'm'
      name_ij(k) = 'PHI1K'
      ia_ij(k) = ia_dga
      scale_ij(k) = BYGRAV
      ir_ij(k) = ir_m190_530
c
      k=k+1 !
      IJ_PHI850 = k ! PHI850 (M**2/S**2-1500*GRAV) 4 DA
      lname_ij(k) = '850 mb HEIGHT'
      units_ij(k) = 'm-1500'
      name_ij(k) = 'PHI850'
      ia_ij(k) = ia_dga
      scale_ij(k) = BYGRAV
      ir_ij(k) = ir_m265_95
c
      k=k+1 !
      IJ_PHI700 = k ! PHI700-3000*GRAV  4 DA
      lname_ij(k) = '700 mb HEIGHT'
      units_ij(k) = 'm-3000'
      name_ij(k) = 'PHI700'
      ia_ij(k) = ia_dga
      scale_ij(k) = BYGRAV
      ir_ij(k) = ir_m530_190
c
      k=k+1 !
      IJ_PHI500 = k ! PHI500-5600*GRAV  4 DA
      lname_ij(k) = '500 mb HEIGHT'
      units_ij(k) = 'm-5600'
      name_ij(k) = 'PHI500'
      ia_ij(k) = ia_dga
      scale_ij(k) = BYGRAV
      ir_ij(k) = ir_m1325_475
c
      k=k+1 !
      IJ_PHI300 = k ! PHI300-9500*GRAV  4 DA
      lname_ij(k) = '300 mb HEIGHT'
      units_ij(k) = 'm-9500'
      name_ij(k) = 'PHI300'
      ia_ij(k) = ia_dga
      scale_ij(k) = BYGRAV
      ir_ij(k) = ir_m2650_950
c
      k=k+1 !
      IJ_PHI100 = k ! PHI100-16400*GRAV 4 DA
      lname_ij(k) = '100 mb HEIGHT'
      units_ij(k) = 'm-16400'
      name_ij(k) = 'PHI100'
      ia_ij(k) = ia_dga
      scale_ij(k) = BYGRAV
      ir_ij(k) = ir_m2650_950
c
      k=k+1 !
      IJ_PHI30 = k ! PHI30-24000*GRAV   4 DA
      lname_ij(k) = '30 mb HEIGHT'
      units_ij(k) = 'm-24000'
      name_ij(k) = 'PHI30'
      ia_ij(k) = ia_dga
      scale_ij(k) = BYGRAV
      ir_ij(k) = ir_m3975_1425
c
      if (kgz_max.gt.k-IJ_PHI1K+1) then
      k=k+1
      IJ_PHI10 = k ! PHI10-30000*GRAV   4 DA
      lname_ij(k) = '10 mb HEIGHT'
      units_ij(k) = 'm-30000'
      name_ij(k) = 'PHI10'
      ia_ij(k) = ia_dga
      scale_ij(k) = BYGRAV
      ir_ij(k) = ir_m5300_1900
      end if
c
      if (kgz_max.gt.k-IJ_PHI1K+1) then
      k=k+1
      IJ_PHI3p4 = k ! PHI3.4-40000*GRAV   4 DA
      lname_ij(k) = '3.4 mb HEIGHT'
      units_ij(k) = 'm-40000'
      name_ij(k) = 'PHI3.4'
      ia_ij(k) = ia_dga
      scale_ij(k) = BYGRAV
      ir_ij(k) = ir_m5300_1900
      end if
c
      if (kgz_max.gt.k-IJ_PHI1K+1) then
      k=k+1
      IJ_PHI0p7 = k ! PHI0.7-50000*GRAV   4 DA
      lname_ij(k) = '0.7 mb HEIGHT'
      units_ij(k) = 'm-50000'
      name_ij(k) = 'PHI0.7'
      ia_ij(k) = ia_dga
      scale_ij(k) = BYGRAV
      ir_ij(k) = ir_m5300_1900
      end if
c
      if (kgz_max.gt.k-IJ_PHI1K+1) then
      k=k+1
      IJ_PHI0p16 = k ! PHI0.16-61000*GRAV   4 DA
      lname_ij(k) = '0.16 mb HEIGHT'
      units_ij(k) = 'm-61000'
      name_ij(k) = 'PHI0.16'
      ia_ij(k) = ia_dga
      scale_ij(k) = BYGRAV
      ir_ij(k) = ir_m5300_1900
      end if
c
      if (kgz_max.gt.k-IJ_PHI1K+1) then
      k=k+1
      IJ_PHI0p07 = k ! PHI0.07-67000*GRAV   4 DA
      lname_ij(k) = '0.07 mb HEIGHT'
      units_ij(k) = 'm-67000'
      name_ij(k) = 'PHI0.07'
      ia_ij(k) = ia_dga
      scale_ij(k) = BYGRAV
      ir_ij(k) = ir_m5300_1900
      end if
c
      if (kgz_max.gt.k-IJ_PHI1K+1) then
      k=k+1
      IJ_PHI0p03 = k ! PHI0.03-72000*GRAV   4 DA
      lname_ij(k) = '0.03 mb HEIGHT'
      units_ij(k) = 'm-72000'
      name_ij(k) = 'PHI0.03'
      ia_ij(k) = ia_dga
      scale_ij(k) = BYGRAV
      ir_ij(k) = ir_m5300_1900
      end if
c
      k=k+1 !
      IJ_P850 = k !
      lname_ij(k) = 'FREQUENCY OF 850mb PRESSURE'  ! weighting function
      units_ij(k) = '%'
      name_ij(k) = 'P850'
      ia_ij(k) = ia_dga
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_T300 = k !
      lname_ij(k) = 'TEMPERATURE AT 300mb'
      units_ij(k) = 'C'
      name_ij(k) = 'T300'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
c
      k=k+1 !
      IJ_T500 = k !
      lname_ij(k) = 'TEMPERATURE AT 500mb'
      units_ij(k) = 'C'
      name_ij(k) = 'T500'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
c
      k=k+1 !
      IJ_T850 = k !
      lname_ij(k) = 'TEMPERATURE AT 850mb x P850'
      units_ij(k) = 'C'
      name_ij(k) = 'T850'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
c
      k=k+1
      IJ_Q300 = k
      lname_ij(k) = 'SPECIFIC HUMIDITY AT 300mb'
      units_ij(k) = 'g/kg'
      name_ij(k) = 'Q300'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1d3
      ir_ij(k) = ir_0_18
c
      k=k+1
      IJ_Q500 = k
      lname_ij(k) = 'SPECIFIC HUMIDITY AT 500mb'
      units_ij(k) = 'g/kg'
      name_ij(k) = 'Q500'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1d3
      ir_ij(k) = ir_0_18
c
      k=k+1
      IJ_Q850 = k
      lname_ij(k) = 'SPECIFIC HUMIDITY AT 850mb x P850'
      units_ij(k) = 'g/kg'
      name_ij(k) = 'Q850'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1d3
      ir_ij(k) = ir_0_18
c
      k=k+1 !
      IJ_RH1 = k !
      lname_ij(k) = 'LAYER 1 RELATIVE HUMIDITY'
      units_ij(k) = '%'
      name_ij(k) = 'RH1'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1d2
      ir_ij(k) = ir_pct
c
      k=k+1
      IJ_RH300 = k
      lname_ij(k) = 'RELATIVE HUMIDITY AT 300mb'
      units_ij(k) = '%'
      name_ij(k) = 'RH300'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1d2
      ir_ij(k) = ir_pct
c
      k=k+1
      IJ_RH500 = k
      lname_ij(k) = 'RELATIVE HUMIDITY AT 500mb'
      units_ij(k) = '%'
      name_ij(k) = 'RH500'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1d2
      ir_ij(k) = ir_pct
c
      k=k+1
      IJ_RH850 = k
      lname_ij(k) = 'RELATIVE HUMIDITY AT 850mb x P850'
      units_ij(k) = '%'
      name_ij(k) = 'RH850'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1d2
      ir_ij(k) = ir_pct
c
      k=k+1 !
      IJ_PMCCLD = k ! PCLDMC (1)  (COMPOSITE OVER ATMOSPHERE) 2 RD
      lname_ij(k) = 'CONVECTIVE CLOUD COVER'
      units_ij(k) = '%'
      name_ij(k) = 'PMCCLD'
      ia_ij(k) = ia_rad
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_CLDTPPR = k ! P-CLOUD TOP   (100 PA)                  2 RD
      lname_ij(k) = 'CLOUD TOP PRESSURE x TOTAL CLOUD COVER'
      units_ij(k) = 'mb'
      name_ij(k) = 'CLDTPPR'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
cc    iw_ij(k) = iw_cldcv  ! built in
      ir_ij(k) = ir_0_1775
c
      k=k+1 !
      IJ_CLDCV = k ! PCLD (1)  (COMPOSITE OVER ATMOSPHERE)   2 RD
      lname_ij(k) = 'TOTAL CLOUD COVER'
      units_ij(k) = '%'
      name_ij(k) = 'CLDCV'
      ia_ij(k) = ia_rad
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_DSEV  = k
       ! 16*P4*(SHA*T4+Z4)*V1*DSIG*DXV (100 W*M/S**2) (UV GRID) 4 DA
      lname_ij(k) = 'TOTAL NT DRY STAT ENRGY' !  NT: NORTHWARD TRANSPORT
      units_ij(k) = '10^14 W'
      name_ij(k) = 'DSEV'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1.d-14*100.*BYGRAV/16.
      jgrid_ij(k) = 2
      ir_ij(k) = ir_m95_265
c
      k=k+1 !
      IJ_TRNFP0 = k ! TRNFP0 (W/m**2)                         2 RS
      lname_ij(k) = 'NET THERMAL RADIATION, TOA'   ! >0 if down !
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'TRNFP0'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_m530_190
c
      k=k+1 !
      IJ_SRTR = k ! SRHDT+TRHDT (J/m**2)                    1 RD/SF
      lname_ij(k) = 'NET RADIATION AT GROUND'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'SRTR'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      ir_ij(k) = ir_m95_265
c
      k=k+1 !
      IJ_NETH = k    ! SRHDT+TRHDT+SHDT+EVHDT+ENRGP (J/m**2)   1 SC
      lname_ij(k) = 'NET HEATING AT GROUND'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'NETH'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      ir_ij(k) = ir_m530_190
c
      k=k+1 !
      IJ_SRNFP0 = k ! SRNFP0 (W/m**2)                         2 RD
      lname_ij(k) = 'NET SOLAR RADIATION, TOA'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'SRNFP0'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
c
      k=k+1 !
      IJ_SRINCP0 = k ! SRINCP0 (W/m**2)                        2 RD
      lname_ij(k) = 'INCIDENT SOLAR RADIATION, TOA'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'SRINCP0'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
c
      k=k+1 !
      IJ_SRNFG = k   ! SRNFG (W/m**2)                          2 RD
      lname_ij(k) = 'NET SOLAR RADIATION, SURF'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'SRNFG'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
c
      k=k+1 !
      IJ_SRINCG = k ! SRINCG (W/m**2)                         2 RD
      lname_ij(k) = 'INCIDENT SOLAR RADIATION, SURF'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'SRINCG'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
c
      k=k+1 !
      IJ_TG1  = k ! TG1 (K-TF)                                1 GD
      lname_ij(k) = 'GROUND TEMPERATURE'
      units_ij(k) = 'C'
      name_ij(k) = 'TG1'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
c
      k=k+1 !
      IJ_RSIT = k ! POICE+PLICE+(IF SNOW)PEARTH               4 DA
      lname_ij(k) = 'SNOW AND ICE COVERAGE'
      units_ij(k) = '%'
      name_ij(k) = 'RSIT'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_TDSL = k ! DIURNAL DELTA TS (K) OVER SOIL (NO PRT) .5*9 MN
      lname_ij(k) = 'DIURNAL SURF AIR TEMP RANGE OVER SOIL'
      units_ij(k) = 'K'
      name_ij(k) = 'TDSL'
      ia_ij(k) = ia_12hr
      scale_ij(k) = 2.
      iw_ij(k) = iw_soil
      ir_ij(k) = ir_0_18
c
      k=k+1
      IJ_TDCOMP = k
      lname_ij(k) = 'DIURNAL SURF AIR TEMP RANGE' ! composite
      units_ij(k) = 'C'
      name_ij(k) = 'TDCOMP'
      ia_ij(k) = ia_12hr
      scale_ij(k) = 2.
      ir_ij(k) = ir_0_18
c
      k=k+1 !
      IJ_DTDP = k ! DTHETA/DPHI (K S**2/m**2) IN TROPOSPHERE  4 DA
      lname_ij(k) = 'TROP STATIC STABILITY'
      units_ij(k) = 'C/km'
      name_ij(k) = 'DTDP'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1000.*GRAV*P1000K
      ir_ij(k) = ir_0_18
c
      k=k+1 !
      IJ_SSTABX = k ! PEAK DTHETA/DPHI (K S**2/m**2) IN PBL   1 CL
      lname_ij(k) = 'PEAK STATIC STABILITY IN PBL'
      units_ij(k) = 'C/km'
      name_ij(k) = 'SSTABX'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.*GRAV*P1000K
      ir_ij(k) = ir_0_18
c
      k=k+1 !
      IJ_RUNE = k ! RUN1 OVER EARTH  (KG/m**2)                1 PG
      lname_ij(k) = 'GROUND RUNOFF OVER SOIL'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'RUNE'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
      iw_ij(k) = iw_soil
      ir_ij(k) = ir_m1_3
c
      k=k+1 !
      IJ_RUNLI = k ! RUN1 OVER LAND ICE  (KG/m**2) (NO PRT)    1 PG
      lname_ij(k) = 'SURFACE RUNOFF OVER LAND ICE'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'RUNLI'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
      iw_ij(k) = iw_lice
      ir_ij(k) = ir_m1_3
c
      k=k+1 !
      IJ_WS   = k ! SURFACE WIND SPEED (M/S)                  3 SF
      lname_ij(k) = 'SURFACE WIND SPEED'
      units_ij(k) = 'm/s'
      name_ij(k) = 'WS'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_18
c
      k=k+1 !
      IJ_TS   = k ! TS (K-TF)                                 3 SF
      lname_ij(k) = 'SURFACE AIR TEMPERATURE'
      units_ij(k) = 'C'
      name_ij(k) = 'TS'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
c
      k=k+1 !
      IJ_US   = k ! US (M/S)                                  3 SF
      lname_ij(k) = 'U COMPONENT OF SURFACE AIR WIND'
      units_ij(k) = 'm/s'
      name_ij(k) = 'US'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m9_26
c
      k=k+1 !
      IJ_VS   = k ! VS (M/S)                                  3 SF
      lname_ij(k) = 'V COMPONENT OF SURFACE AIR WIND'
      units_ij(k) = 'm/s'
      name_ij(k) = 'VS'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m9_26
c
      k=k+1 !
      IJ_SLP  = k ! PSL (100 PA-1000)  (USING TS)             4 DA
      lname_ij(k) = 'SEA LEVEL PRESSURE'
      units_ij(k) = 'mb-1000'
      name_ij(k) = 'SLP'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1.
      ir_ij(k) = ir_m9_26
c
      k=k+1 !
      IJ_UJET = k ! UJET (M/S)                                4 DA
      lname_ij(k) = 'U COMPONENT OF JET WINDS'
      units_ij(k) = 'm/s'
      name_ij(k) = 'UJET'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1.
      jgrid_ij(k) = 2
      ir_ij(k) = ir_m38_106
c
      k=k+1 !
      IJ_VJET = k ! VJET (M/S)                                4 DA
      lname_ij(k) = 'V COMPONENT OF JET WINDS'
      units_ij(k) = 'm/s'
      name_ij(k) = 'VJET'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1.
      jgrid_ij(k) = 2
      ir_ij(k) = ir_m38_106
c
      k=k+1 !
      IJ_PCLDL = k ! PCLD(LOW) (1)                            2 RD
      lname_ij(k) = 'LOW LEVEL CLOUDINESS'
      units_ij(k) = '%'
      name_ij(k) = 'PCLDL'
      ia_ij(k) = ia_rad
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_PCLDM = k ! PCLD(MID) (1)                            2 RD
      lname_ij(k) = 'MIDDLE LEVEL CLOUDINESS'
      units_ij(k) = '%'
      name_ij(k) = 'PCLDM'
      ia_ij(k) = ia_rad
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_PCLDH = k ! PCLD(HIGH) (1)                           2 RD
      lname_ij(k) = 'HIGH LEVEL CLOUDINESS'
      units_ij(k) = '%'
      name_ij(k) = 'PCLDH'
      ia_ij(k) = ia_rad
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_BTMPW = k ! BTEMPW-TF (K-TF)                         2 RD
      lname_ij(k) = 'BRIGHTNESS TEMP THRU WNDW' ! window region
      units_ij(k) = 'C'
      name_ij(k) = 'BTMPW'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
c
      k=k+1 !
      IJ_SRREF = k ! PLAVIS*S0*COSZ (W/m**2)                  2 RD
      lname_ij(k) = 'REFLECTED SOLAR RADIATION IN VISUAL'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'SRREF'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
c
      k=k+1
      IJ_SRVIS = k ! ALBVIS*S0*COSZ (W/m**2)                  2 RD
      lname_ij(k) = 'REFLECTED SOLAR RADIATION IN VISUAL AT SURF'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'SRVIS'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
c
      k=k+1 !
      IJ_TOC2 = k ! TGO2= TOCEAN(2)  (C)                   .5*9 MN
      lname_ij(k) = 'OCEAN TEMPERATURE BELOW MIXED LAYER' ! lyr 2
      units_ij(k) = 'C'
      name_ij(k) = 'TOC2'
      ia_ij(k) = ia_12hr
      scale_ij(k) = 2.
      iw_ij(k) = iw_ocn
      ir_ij(k) = ir_m9_26
c
      k=k+1 !
      IJ_TAUS = k ! TAUS  (MOM. SURF. DRAG) (kg/m**2) (NO PRT)  3 SF
      lname_ij(k) = 'MAG OF MOMENTUM SURFACE DRAG'
      units_ij(k) = 'g/m*s^2'
      name_ij(k) = 'TAUS'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1000.
      ir_ij(k) = ir_0_710
c
      k=k+1 !
      IJ_TAUUS = k ! TAUUS (MOM. SURF. DRAG) (kg/m**2) (NO PRT)  3 SF
      lname_ij(k) = 'U COMPON OF MOMENTUM SRF DRAG'
      units_ij(k) = 'g/m*s^2'
      name_ij(k) = 'TAUUS'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1000.
      ir_ij(k) = ir_m2650_950
c
      k=k+1 !
      IJ_TAUVS = k ! TAUVS (MOM. SURF. DRAG) (kg/m**2) (NO PRT)  3 SF
      lname_ij(k) = 'V COMPON OF MOMENTUM SRF DRAG'
      units_ij(k) = 'g/m*s^2'
      name_ij(k) = 'TAUVS'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1000.
      ir_ij(k) = ir_m2650_950
c
      k=k+1 !
      IJ_GWTR = k ! WATER1+WATER2+ICE1+ICE2 (EARTH POINTS ONLY) 1 GD
      lname_ij(k) = 'TOTAL EARTH WATER' ! includes ice
      units_ij(k) = 'kg/m^2'
      name_ij(k) = 'GWTR'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      iw_ij(k) = iw_soil
      ir_ij(k) = ir_0_710
c
      k=k+1 !
      IJ_QS   = k ! QS                                (NO PRT)  3 SF
      lname_ij(k) = 'SURFACE AIR SPECIFIC HUMIDITY'
      units_ij(k) = '10^-4 g/g'
      name_ij(k) = 'QS'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.d4
      ir_ij(k) = ir_0_180
c
      k=k+1 !
      IJ_STRNGTS = k ! MAX(0,65F-TS_daily_avg in F)          .5*9 MN
      lname_ij(k) = 'MONTHLY HEATING' ! monthly heating need ?
      units_ij(k) = 'degF days'
      name_ij(k) = 'STRNGTS'
      ia_ij(k) = ia_12hr
      scale_ij(k) = 2.*30.
      ir_ij(k) = ir_0_3550
c
      k=k+1 !
      IJ_ARUNU = k ! ARUNU                                      1 EA
      lname_ij(k) = 'UNDERGROUND RUNOFF OVER SOIL'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'ARUNU'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
      iw_ij(k) = iw_soil
      ir_ij(k) = ir_m1_3
c
      k=k+1 !
      IJ_DTGDTS = k ! 18*(DEL(TG)/DEL(TS)-1),DEL= DIURN MX-MN .5*9 MN
      lname_ij(k) = 'PLANT WATER STRESS'
      units_ij(k) = '1'
      name_ij(k) = 'DTGDTS'
      ia_ij(k) = ia_12hr
      scale_ij(k) = 2.*30.
      ir_ij(k) = ir_m190_530
c
      k=k+1 !
      IJ_PUQ  = k ! 8*P*U*Q (VERT. INTEGRATED) (12.5 PA*M/S) 4 DA
      lname_ij(k) = 'EAST-WEST HUMIDITY FLUX (VERT SUM)'
      units_ij(k) = 'mb*m/s'
      name_ij(k) = 'PUQ'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1./8.
      ir_ij(k) = ir_m45_130
c
      k=k+1 !
      IJ_PVQ  = k ! 8*P*V*Q (VERT. INTEGRATED) (12.5 PA*M/S) 4 DA
      lname_ij(k) = 'NORTH-SOUTH HUMIDITY FLUX (VERT SUM)'
      units_ij(k) = 'mb*m/s'
      name_ij(k) = 'PVQ'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1./8.
      ir_ij(k) = ir_m45_130
c
      k=k+1 !
      IJ_TGO  = k               !      3 SF
      lname_ij(k) = 'SEA SURFACE TEMPERATURE'    ! layer 1
      units_ij(k) = 'C'
      name_ij(k) = 'TGO'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      iw_ij(k) = iw_ocn
      ir_ij(k) = ir_m9_26
c
      k=k+1 !
      IJ_MSI2 = k ! ACE2OI= MSI2*POICE  (KG/m**2)        1 GD
      lname_ij(k) = 'LAYER 2 OCEAN ICE MASS x POICE'
      units_ij(k) = '10^3 kg/m^2'
      name_ij(k) = 'MSI2'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-3
cc    iw_ij(k) = iw_oic !! built in
      ir_ij(k) = ir_0_4
c
c     k=k+1 !
c     IJ_WLM  = k ! WIND SPEED IN TOP LAYER (M/S) before SDRAG 1 SD
c     lname_ij(k) = 'WIND SPEED IN TOP LAYER'
c     units_ij(k) = 'm/s'
c     name_ij(k) = 'WLM'
c     ia_ij(k) = ia_src
c     scale_ij(k) = 1.
c     jgrid_ij(k) = 2
c     ir_ij(k) = ir_0_26_150
c
      k=k+1 !
      IJ_TGO2 = k ! TGO12= TOCEAN(3) (C)                  .5*9 MN
      lname_ij(k) = 'OCEAN TEMPERATURE AT ANN-MAX MIXED-LAYER' ! layer 3
      units_ij(k) = 'C'
      name_ij(k) = 'TGO2'
      ia_ij(k) = ia_12hr
      scale_ij(k) = 2.
      iw_ij(k) = iw_ocn
      ir_ij(k) = ir_m9_26
c
      k=k+1 !
      IJ_EVAPO = k ! EVAP*POCEAN  (KG/m**2)                  1 GD
      lname_ij(k) = 'OCEAN EVAPORATION x POCEAN'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'EVAPO'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
c     iw built-in
      ir_ij(k) = ir_0_3_15
c
      k=k+1 !
      IJ_EVAPI = k ! EVAP*POICE  (KG/m**2)                   1 GD
      lname_ij(k) = 'OCEAN ICE EVAPORATION x POICE'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'EVAPI'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
c     iw built-in
      ir_ij(k) = ir_0_3_15
c
      k=k+1 !
      IJ_EVAPLI = k ! EVAP OVER LAND ICE  (KG/m**2)          1 GD
      lname_ij(k) = 'LAND ICE EVAPORATION x PLICE'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'EVAPLI'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
c     iw built-in
      ir_ij(k) = ir_0_3_15
c
      k=k+1 !
      IJ_EVAPE = k ! EVAP OVER EARTH  (KG/m**2)              1 GD
      lname_ij(k) = 'SOIL EVAPORATION x PSOIL'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'EVAPE'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
c     iw built-in
      ir_ij(k) = ir_0_3_15
c
      k=k+1 !
      IJ_F0OC = k ! F0DT*POCEAN, NET HEAT AT Z0  (J/m**2)    1 GD
      lname_ij(k) = 'NET HEAT INTO OCEAN x POCEAN'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'F0OC'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
c     iw built-in
      ir_ij(k) = ir_m530_190
c
      k=k+1 !
      IJ_F0OI = k ! F0DT*POICE, NET HEAT AT Z0  (J/m**2)     1 GD
      lname_ij(k) = 'NET HEAT INTO OCEAN ICE x POICE'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'F0OI'
      ia_ij(k) = ia_src
      scale_ij(k) = (1./DTsrc)
c     iw built-in
      ir_ij(k) = ir_m530_190
c
      k=k+1 !
      IJ_F0LI = k ! F0DT, NET HEAT AT Z0 OVER LAND ICE  (J/m**2) 1 GD
      lname_ij(k) = 'NET HEAT INTO LAND ICE'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'F0LI'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      iw_ij(k) = iw_lice
      ir_ij(k) = ir_m530_190
c
      k=k+1 !
      IJ_F0E  = k ! F0DT, NET HEAT AT Z0 OVER EARTH  (J/m**2) 1 GD
      lname_ij(k) = 'NET HEAT INTO SOIL'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'F0E'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      iw_ij(k) = iw_soil
      ir_ij(k) = ir_m530_190
c
      k=k+1 !
      IJ_F1LI = k ! F1DT OVER LAND ICE  (J/m**2)             1 PG
      lname_ij(k) = 'CONDUCTION AT LYR 1 BOTTOM OVER LAND ICE'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'F1LI'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      iw_ij(k) = iw_lice
      ir_ij(k) = ir_m530_190
c
      k=k+1 !
      IJ_SNWF = k ! SNOW FALL  (KG/m**2)                     1 PR
      lname_ij(k) = 'SNOW FALL (H2O EQUIV)'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'SNWF'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
      ir_ij(k) = ir_0_3_15
c
      k=k+1 !
      IJ_TSLI = k ! SURF AIR TEMP OVER LAND ICE  (C)  NISURF*1 SF
      lname_ij(k) = 'SURF AIR TEMP OVER LAND ICE'
      units_ij(k) = 'C'
      name_ij(k) = 'TSLI'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d0/NIsurf
      iw_ij(k) = iw_lice
      ir_ij(k) = ir_m80_28
c
      k=k+1 !
      IJ_ERUN2 = k ! F2DT OVER LAND ICE  (J/m**2)            1 PG
      lname_ij(k) = 'TRANSPORT OF ENERGY AT LYR 2 BOTTOM OVER LAND ICE'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'ERUN2'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      iw_ij(k) = iw_lice
      ir_ij(k) = ir_m80_28
c
      k=k+1 !
      IJ_SHDTLI = k ! SHDT OVER LAND ICE  (J/m**2)           1 SF
      lname_ij(k) = 'SENS HEAT FLUX OVER LAND ICE'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'SHDTLI'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      iw_ij(k) = iw_lice
      ir_ij(k) = ir_m38_106
c
      k=k+1 !
      IJ_EVHDT = k ! EVHDT OVER LAND ICE  (J/m**2)           1 SF
      lname_ij(k) = 'LATENT HEAT FLUX OVER LAND ICE'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'EVHDT'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      iw_ij(k) = iw_lice
      ir_ij(k) = ir_m38_106
c
      k=k+1 !
      IJ_TRHDT = k ! TRHDT OVER LAND ICE  (J/m**2)           1 SF
      lname_ij(k) = 'NET THERMAL RADIATION INTO LAND ICE'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'TRHDT'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      iw_ij(k) = iw_lice
      ir_ij(k) = ir_m38_106
c
      k=k+1 !
      IJ_TMNMX  = k ! MIN(DIURNAL MAX OF COMPOSITE TS)      12 MN
      lname_ij(k) = 'SURFC AIR TEMPERATURE: LOWEST DIURNAL HIGH  + TF'
      units_ij(k) = 'C' ! after offset (TF=273.16)
      name_ij(k) = 'TMNMX'
      ia_ij(k) = ia_inst
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
c
      k=k+1 !
      IJ_PEVAP = k ! POTENTIAL EVAPORATION (KG/m**2)         1 EA
      lname_ij(k) = 'POTENTIAL EVAPORATION'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'PEVAP'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
      ir_ij(k) = ir_0_26_150
c
      k=k+1 !
      IJ_TMAXE = k ! MAX TS OVER EARTH FOR CURRENT DAY (C).5*9 MN
      lname_ij(k) = 'SURFACE AIR TEMPERATURE: DIURNAL HIGH/SOIL'
      units_ij(k) = 'C'
      name_ij(k) = 'TMAXE'
      ia_ij(k) = ia_12hr
      scale_ij(k) = 2.
      iw_ij(k) = iw_soil
      ir_ij(k) = ir_m80_28
c
      k=k+1 !
      IJ_WMSUM = k ! LIQUID WATER PATH (kg/m**2)             1 CL
      lname_ij(k) = 'LIQUID WATER PATH'
      units_ij(k) = '.1 kg/m^2'
      name_ij(k) = 'WMSUM'
      ia_ij(k) = ia_src
      scale_ij(k) = 10.
      ir_ij(k) = ir_0_18
c
      k=k+1 !
      IJ_PSCLD = k ! SHALLOW CONVECTIVE CLOUD COVER  (1)     1 CL
      lname_ij(k) = 'SHALLOW CONVECTIVE CLOUD COVER'
      units_ij(k) = '%'
      name_ij(k) = 'PSCLD'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_PDCLD = k ! DEEP CONVECTIVE CLOUD COVER     (1)     1 CL
      lname_ij(k) = 'DEEP CONVECTIVE CLOUD COVER'
      units_ij(k) = '%'
      name_ij(k) = 'PDCLD'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_DCNVFRQ = k ! DEEP CONVECTIVE CLOUD OCCURRENCE (1)   1 CL
      lname_ij(k) = 'DEEP CONVECTIVE CLOUD FREQUENCY'
      units_ij(k) = '%'
      name_ij(k) = 'DCNVFRQ'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_SCNVFRQ = k ! SHALLOW CONVECTIVE CLOUD OCCURRENCE (1) 1 CL
      lname_ij(k) = 'SHALLOW CONV CLOUD FREQUENCY'
      units_ij(k) = '%'
      name_ij(k) = 'SCNVFRQ'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
c     k=k+1 !
cfree IJ_EMTMOM = k ! INCIDENT MTN EAST MOM. FLUX (MB-M/S**2)  1 SD
cfree lname_ij(k) = 'INCIDENT MTN EAST MOMENTUM FLUX'
cfree units_ij(k) = 'mb m/s^2'
cfree name_ij(k) = 'EMTMOM'
cfree ia_ij(k) = ia_src
cfree scale_ij(k) = 1./DTsrc
c
c     k=k+1 !
cfree IJ_SMTMOM = k ! INCIDENT MTN SOUTH MOM. FLUX (MB-M/S**2) 1 SD
cfree lname_ij(k) = 'INCIDENT MTN SOUTH MOMENTUM FLUX'
cfree units_ij(k) = 'mb m/s^2'
cfree name_ij(k) = 'SMTMOM'
cfree ia_ij(k) = ia_src
cfree scale_ij(k) = 1./DTsrc
c
      k=k+1 !
      IJ_FMU  = k ! EAST-WEST MASS FLUX (KG/S) 100./GRAV/3600.*1 DY
      lname_ij(k) = 'EAST-WEST MASS FLUX'
      units_ij(k) = '10^10 kg/s'
      name_ij(k) = 'FMU'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-10*100.*BYGRAV/DTsrc
      ir_ij(k) = ir_m38_106
c
      k=k+1 !
      IJ_FMV  = k ! NORTH-SOUTH MASS FLUX (KG/S) 100./GRAV/3600.*1 DY
      lname_ij(k) = 'NORTH-SOUTH MASS FLUX'
      units_ij(k) = '10^10 kg/s'
      name_ij(k) = 'FMV'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-10*100.*BYGRAV/DTsrc
c
      k=k+1 !
      IJ_FGZU = k ! EAST-WEST GEOPOTENTIAL FLUX (W)    /3600.*1 DY
      lname_ij(k) = 'EAST-WEST GEOPOTENTIAL FLUX'
      units_ij(k) = '10^10 W'
      name_ij(k) = 'FGZU'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-10/DTsrc
      ir_ij(k) = ir_m1325_475
c
      k=k+1 !
      IJ_FGZV = k ! NORTH-SOUTH GEOPOTENTIAL FLUX (W)  /3600.*1 DY
      lname_ij(k) = 'NORTH-SOUTH GEOPOTENTIAL FLUX'
      units_ij(k) = '10^10 W'
      name_ij(k) = 'FGZV'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-10/DTsrc
      ir_ij(k) = ir_m1325_475
c
      k=k+1 !
      IJ_ERVR = k ! Energy Outflow by Rivers (10**10 W) E-10/DTS*1 RV
      lname_ij(k) = 'Energy Outflow by Rivers'
      units_ij(k) = '10^10 W'
      name_ij(k) = 'ERVR'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-10/DTsrc
      ir_ij(k) = ir_m1325_475
c
      k=k+1 !
      IJ_MRVR = k ! Mass Outflow by Rivers (10**5 kg/s)  E-5/DTS*1 RV
      lname_ij(k) = 'Mass Outflow by Rivers'
      units_ij(k) = '10^5 kg/s'
      name_ij(k) = 'MRVR'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-5/DTsrc
      ir_ij(k) = ir_m1325_475
c
c     k=k+1 !
c     IJ_SDRAG = k ! DU/DT BY SDRAG (M S-2)                       1 SD
c     lname_ij(k) = 'STRATOSPHERIC DRAG'
c     units_ij(k) = '10^-6 m/s^2'
c     name_ij(k) = 'SDRAG'
c     ia_ij(k) = ia_src
c     scale_ij(k) = 1.d6/DTsrc
c     ir_ij(k) = ir_m190_530
c
      k=k+1 !
      IJ_LKON = k
      lname_ij(k) = 'LAST ICE-FREE DAY (SH-58,NH-242) x LKICE'
      units_ij(k) = 'JULIAN DAY'
      name_ij(k) = 'LKONDAY'
      ia_ij(k) = ia_inst
c
      k=k+1 !
      IJ_LKOFF = k
      lname_ij(k) = 'LAST ICED-UP DAY (SH-58,NH-242) x LKICE'
      units_ij(k) = 'JULIAN DAY'
      name_ij(k) = 'LKOFFDAY'
      ia_ij(k) = ia_inst
c
      k=k+1 !
      IJ_LKICE = k
      lname_ij(k) = 'LAKE ICE WEIGHTING' ! for lake freeze/thaw diags
      units_ij(k) = ' '
      name_ij(k) = 'LKICEWT'
      ia_ij(k) = ia_inst
c
C**** Here I am adding all the previous AIJG to AIJ
C**** actual number for Gxx = 100 + xx
      k=k+1
      IJ_G01 = k
      name_ij(k) = 'bs_wlay1' !
      lname_ij(k) = 'LAYER 1 BARE SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      iw_ij(k) = iw_bare
      ir_ij(k) = ir_0_710
c
      k=k+1
      IJ_G02 = k
      name_ij(k) = 'bs_wlay2' !
      lname_ij(k) = 'LAYER 2 BARE SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      iw_ij(k) = iw_bare
      ir_ij(k) = ir_0_710
c
      k=k+1
      IJ_G03 = k
      name_ij(k) = 'bs_wlay3' !
      lname_ij(k) = 'LAYER 3 BARE SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      iw_ij(k) = iw_bare
      ir_ij(k) = ir_0_710
c
      k=k+1
      IJ_G04 = k
      name_ij(k) = 'bs_wlay4' !
      lname_ij(k) = 'LAYER 4 BARE SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      iw_ij(k) = iw_bare
      ir_ij(k) = ir_0_710
c
      k=k+1
      IJ_G05 = k
      name_ij(k) = 'bs_beta' !
      lname_ij(k) = 'BARE SOIL WETNESS, BETA'
      units_ij(k) = '%'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
      iw_ij(k) = iw_bare
c
      k=k+1
      IJ_G06 = k
      name_ij(k) = 'beta_pen' !
      lname_ij(k) = 'PENMAN SOIL WETNESS, BETA'
      units_ij(k) = '%'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
      iw_ij(k) = iw_soil
c
      k=k+1
      IJ_G07 = k
      name_ij(k) = 'vs_wcan' !
      lname_ij(k) = 'VEGETATION CANOPY SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      iw_ij(k) = iw_veg
      ir_ij(k) = ir_0_710
c
      k=k+1
      IJ_G08 = k
      name_ij(k) = 'vs_wlay1' !
      lname_ij(k) = 'LAYER 1 VEGETATED SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      iw_ij(k) = iw_veg
      ir_ij(k) = ir_0_710
c
      k=k+1
      IJ_G09 = k
      name_ij(k) = 'vs_wlay2' !
      lname_ij(k) = 'LAYER 2 VEGETATED SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      iw_ij(k) = iw_veg
      ir_ij(k) = ir_0_710
c
      k=k+1
      IJ_G10 = k
      name_ij(k) = 'vs_wlay3' !
      lname_ij(k) = 'LAYER 3 VEGETATED SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      iw_ij(k) = iw_veg
      ir_ij(k) = ir_0_710
c
      k=k+1
      IJ_G11 = k
      name_ij(k) = 'bvs_wet' !
      lname_ij(k) = 'BARE & VEGETATED SOIL WETNESS'
      units_ij(k) = '%'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
      iw_ij(k) = iw_soil
c
      k=k+1
      IJ_G12 = k
      name_ij(k) = 'cond_atm' !
      lname_ij(k) = 'CONDUCTANCE OF ATMOSPHERE'
      units_ij(k) = 'm/s'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      iw_ij(k) = iw_soil
      ir_ij(k) = ir_0_71
c
      k=k+1
      IJ_G13 = k
      name_ij(k) = 'cond_can' !
      lname_ij(k) = 'CONDUCTANCE OF CANOPY'
      units_ij(k) = '.01 m/s'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
      iw_ij(k) = iw_veg
      ir_ij(k) = ir_0_71
c
      k=k+1
      IJ_G14 = k
      name_ij(k) = 'pev_pen' !
      lname_ij(k) = 'PENMAN POTENTIAL EVAPORATION'
      units_ij(k) = 'mm/day'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
      iw_ij(k) = iw_soil
      ir_ij(k) = ir_0_26_150
c
      k=k+1
      IJ_G15 = k
      name_ij(k) = 'bs_tlay1' !
      lname_ij(k) = 'BARE SOIL LAYER 1 TEMPERATURE'
      units_ij(k) = 'C'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./NIsurf
      iw_ij(k) = iw_bare
      ir_ij(k) = ir_m80_28
c
      k=k+1
      IJ_G16 = k
      name_ij(k) = 'bs_tlay2' !
      lname_ij(k) = 'BARE SOIL LAYER 2 TEMPERATURE'
      units_ij(k) = 'C'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./NIsurf
      iw_ij(k) = iw_bare
      ir_ij(k) = ir_m80_28
c
      k=k+1
      IJ_G17 = k
      name_ij(k) = 'bs_tlay3' !
      lname_ij(k) = 'BARE SOIL LAYER 3 TEMPERATURE'
      units_ij(k) = 'C'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./NIsurf
      iw_ij(k) = iw_bare
      ir_ij(k) = ir_m80_28
c
      k=k+1
      IJ_G18 = k
      name_ij(k) = 'bs_evap' !
      lname_ij(k) = 'BARE SOIL EVAPORATION'
      units_ij(k) = 'mm/day'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
      iw_ij(k) = iw_bare
      ir_ij(k) = ir_0_3_15
c
      k=k+1
      IJ_G19 = k
      name_ij(k) = 'drycan_evap' !
      lname_ij(k) = 'DRY CANOPY EVAPORATION'
      units_ij(k) = 'mm/day'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
      iw_ij(k) = iw_veg
      ir_ij(k) = ir_0_3_15
c
      k=k+1
      IJ_G20 = k
      name_ij(k) = 'wetcan_evap' !
      lname_ij(k) = 'WET CANOPY EVAPORATION'
      units_ij(k) = 'mm/day'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
      iw_ij(k) = iw_veg
      ir_ij(k) = ir_0_3_15
c
      k=k+1
      IJ_G21 = k
      name_ij(k) = 'can_temp' !
      lname_ij(k) = 'CANOPY TEMPERATURE'
      units_ij(k) = 'C'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./NIsurf
      iw_ij(k) = iw_veg
      ir_ij(k) = ir_m80_28
c
      k=k+1
      IJ_G22 = k
      name_ij(k) = 'vs_tlay1' !
      lname_ij(k) = 'VEGETATED SOIL LAYER 1 TEMPERATURE'
      units_ij(k) = 'C'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./NIsurf
      iw_ij(k) = iw_veg
      ir_ij(k) = ir_m80_28
c
      k=k+1
      IJ_G23 = k
      name_ij(k) = 'vs_tlay2' !
      lname_ij(k) = 'VEGETATED SOIL LAYER 2 TEMPERATURE'
      units_ij(k) = 'C'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./NIsurf
      iw_ij(k) = iw_veg
      ir_ij(k) = ir_m80_28
c
      k=k+1
      IJ_G24 = k
      name_ij(k) = 'vs_tlay3' !
      lname_ij(k) = 'VEGETATED SOIL LAYER 3 TEMPERATURE'
      units_ij(k) = 'C'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./NIsurf
      iw_ij(k) = iw_veg
      ir_ij(k) = ir_m80_28
c
      k=k+1
      IJ_G25 = k
      name_ij(k) = 'wtbl_depth' !
      lname_ij(k) = 'AVERAGE WATER TABLE DEPTH'
      units_ij(k) = 'm'
      ia_ij(k) = ia_src
      scale_ij(k) = -1./NIsurf
      iw_ij(k) = iw_soil
      ir_ij(k) = ir_0_3_15
c
      k=k+1
      IJ_G26 = k
      name_ij(k) = 'vs_wetness' !
      lname_ij(k) = 'VEGETATED SOIL WETNESS'
      units_ij(k) = '%'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
      iw_ij(k) = iw_veg
c
      k=k+1
      IJ_G27 = k
      name_ij(k) = 'beta_trans' !
      lname_ij(k) = 'TRANSPIRATION EFFICIENCY, BETAT'
      units_ij(k) = '%'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
      iw_ij(k) = iw_veg
c
      k=k+1
      IJ_G28 = k
      name_ij(k) = 'bs_snowdp' !
      lname_ij(k) = 'SNOW DEPTH OVER BARE SOIL'
      units_ij(k) = 'mm H2O'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      iw_ij(k) = iw_veg
      ir_ij(k) = ir_0_3550
c
      k=k+1
      IJ_G29 = k
      name_ij(k) = 'vs_snowdp' !
      lname_ij(k) = 'SNOW DEPTH OVER VEG SOIL'
      units_ij(k) = 'mm H2O'
      scale_ij(k) = 1000.
      ia_ij(k) = ia_src
      iw_ij(k) = iw_veg
      ir_ij(k) = ir_0_3550
c
c Gravity Wave diagnostics
      iDO_GWDRAG = 0
      if (DO_GWDRAG) then
c
      k=k+1                ! ij diags from gwdrag calculations
      IJ_GW1 = k
      name_ij(k) = 'ij_def_drag_mom_flux'
      lname_ij(k) = 'DEFORM. DRAG MOM FLUX'
      units_ij(k) = '.1 N/m**2'
      scale_ij(k) = 10.*100.*BYGRAV
      ia_ij(k) = ia_src
      ir_ij(k) = ir_m1_3
c
      k=k+1                ! ij diags from gwdrag calculations
      IJ_GW2 = k
      name_ij(k) = 'ij_mtn_wave_mom_flux'
      lname_ij(k) = 'MTN WAVE MOM. FLUX'
      units_ij(k) = '.1 N/m**2'  ! dynes/cm^2
      scale_ij(k) = 10.*100.*BYGRAV
      ia_ij(k) = ia_src
      ir_ij(k) = ir_m1_3
c
      k=k+1                ! ij diags from gwdrag calculations
      IJ_GW3 = k
      name_ij(k) = 'ij_shr_wave_mom_flux'
      lname_ij(k) = 'SHEAR WAVE MOM. FLUX'
      units_ij(k) = '.001 N/m**2'
      scale_ij(k) = 1000.*100.*BYGRAV
      ia_ij(k) = ia_src
      ir_ij(k) = ir_m1_3
c
      k=k+1                ! ij diags from gwdrag calculations
      IJ_GW4 = k
      name_ij(k) = 'ij_mc_c_m10r_mom_flux'
      lname_ij(k) = 'MC C=-10R MOM. FLUX'
      units_ij(k) = '.001 N/m**2'
      scale_ij(k) = 1000.*100.*BYGRAV
      ia_ij(k) = ia_src
      ir_ij(k) = ir_m1_3
c
      k=k+1                ! ij diags from gwdrag calculations
      IJ_GW5 = k
      name_ij(k) = 'ij_mc_c_m20r_mom_flux'
      lname_ij(k) = 'MC C=-20R MOM. FLUX'
      units_ij(k) = '.001 N/m**2'
      scale_ij(k) = 1000.*100.*BYGRAV
      ia_ij(k) = ia_src
      ir_ij(k) = ir_m1_3
c
      k=k+1                ! ij diags from gwdrag calculations
      IJ_GW6 = k
      name_ij(k) = 'ij_mc_c_m40r_mom_flux'
      lname_ij(k) = 'MC C=-40R MOM. FLUX'
      units_ij(k) = '.001 N/m**2'
      scale_ij(k) = 1000.*100.*BYGRAV
      ia_ij(k) = ia_src
      ir_ij(k) = ir_m1_3
c
      k=k+1                ! ij diags from gwdrag calculations
      IJ_GW7 = k
      name_ij(k) = 'ij_phase_speed_of_wind_shear'
      lname_ij(k) = 'PHASE SPEED OF SHEAR WAVE'
      units_ij(k) = 'm/s'
      ia_ij(k) = ia_src
      ir_ij(k) = ir_m45_130
c
      k=k+1                ! ij diags from gwdrag calculations
      IJ_GW8 = k
      name_ij(k) = 'ij_source_speed_of_mc'
      lname_ij(k) = 'MC SOURCE WIND SPEED'
      units_ij(k) = 'm/s'
      ia_ij(k) = ia_src
      ir_ij(k) = ir_m45_130
c
      k=k+1                ! ij diags from gwdrag calculations
      IJ_GW9 = k
      name_ij(k) = 'ij_exit_tot_mom_flux'
      lname_ij(k) = 'EXIT TOT. MOM. FLUX'
      units_ij(k) = '.0001 N/m**2'
      scale_ij(k) = 10000.*100.*BYGRAV
      ia_ij(k) = ia_src
      ir_ij(k) = ir_m1_3
      iDO_GWDRAG = k-IJ_GW1+1
      END IF

      if (isccp_diags.eq.1) then
      k=k+1 !
      IJ_LCLDI = k
      lname_ij(k) = 'LOW LEVEL CLOUDINESS (ISCCP)'
      units_ij(k) = '%'
      name_ij(k) = 'LCLDI'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_MCLDI = k
      lname_ij(k) = 'MIDDLE LEVEL CLOUDINESS (ISCCP)'
      units_ij(k) = '%'
      name_ij(k) = 'MCLDI'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_HCLDI = k
      lname_ij(k) = 'HIGH LEVEL CLOUDINESS (ISCCP)'
      units_ij(k) = '%'
      name_ij(k) = 'HCLDI'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
C**** Note this diagnostic is NOT the total cloud cover (that is got by
C**** summing the low+ mid+high diagnostics). Instead, this is the
C**** fraction of time that a cloud appears in the grid box (which may
C**** well cover less than 100% of the box). This is needed for
C**** weighting the cloud top pressure and optical depth
      k=k+1 !
      IJ_TCLDI = k
      lname_ij(k) = 'FRACTION OF TIME FOR ISCCP CLOUD'
      units_ij(k) = '%'
      name_ij(k) = 'TCLDI'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_CTPI = k
      lname_ij(k) = 'CLOUD TOP PRESSURE (ISCCP) x TOTAL ISCCP CLOUD'
      units_ij(k) = 'mb'
      name_ij(k) = 'CTPI'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_1775
c
      k=k+1 !
      IJ_TAUI = k
      lname_ij(k) = 'CLOUD OPTICAL DEPTH (ISCCP) x TOTAL ISCCP CLOUD'
      units_ij(k) = ''
      name_ij(k) = 'TAUI'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.

C**** Also include MSU radiation diagnotsics here

c      k=k+1 !
c      IJ_MSU2 = k
c      lname_ij(k) = 'MSU CHANNEL 2'
c      units_ij(k) = 'C'
c      name_ij(k) = 'MSU2'
c      ia_ij(k) = ia_inst
c      ir_ij(k) = ir_m80_28
c      scale_ij(k) = 1.

c      k=k+1 !
c      IJ_MSU2R = k
c      lname_ij(k) = 'MSU CHANNEL 2R'
c      units_ij(k) = 'C'
c      name_ij(k) = 'MSU2R'
c      ia_ij(k) = ia_inst
c      ir_ij(k) = ir_m80_28
c      scale_ij(k) = 1.

c      k=k+1 !
c      IJ_MSU3 = k
c      lname_ij(k) = 'MSU CHANNEL 3'
c      units_ij(k) = 'C'
c      name_ij(k) = 'MSU3'
c      ia_ij(k) = ia_inst
c      ir_ij(k) = ir_m80_28
c      scale_ij(k) = 1.

c      k=k+1 !
c      IJ_MSU4 = k
c      lname_ij(k) = 'MSU CHANNEL 4'
c      units_ij(k) = 'C'
c      name_ij(k) = 'MSU4'
c      ia_ij(k) = ia_inst
c      ir_ij(k) = ir_m80_28
c      scale_ij(k) = 1.

      end if

      k=k+1
      IJ_PTROP = k
      lname_ij(k) = 'TROPOPAUSE PRESSURE (WMO)'
      units_ij(k) = 'mb'
      name_ij(k) = 'PTROP'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.

      k=k+1
      IJ_TTROP = k
      lname_ij(k) = 'TROPOPAUSE TEMPERATURE (WMO)'
      units_ij(k) = 'K'
      name_ij(k) = 'TTROP'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.

      k=k+1
      IJ_TSI = k
      lname_ij(k) = 'SEA ICE TEMPERATURE (MASS LAYER 2) x POICE'
      units_ij(k) = 'C'
      name_ij(k) = 'TEMPSI'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.

      k=k+1
      IJ_SSI1 = k
      lname_ij(k) = 'SEA ICE SALINITY (MASS LAYER 1) x POICE'
      units_ij(k) = 'psu'
      name_ij(k) = 'SSI1'
      ia_ij(k) = ia_src
      scale_ij(k) = 1d3

      k=k+1
      IJ_SSI2 = k
      lname_ij(k) = 'SEA ICE SALINITY (MASS LAYER 2) x POICE'
      units_ij(k) = 'psu'
      name_ij(k) = 'SSI2'
      ia_ij(k) = ia_src
      scale_ij(k) = 1d3

      k=k+1
      IJ_MLTP = k
      lname_ij(k) = 'SEA ICE MELT POND MASS x POICE'
      units_ij(k) = 'kg/m^2'
      name_ij(k) = 'MLTP'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.

      k=k+1
      IJ_FRMP = k
      lname_ij(k) = 'SEA ICE MELT POND FRACTION x POICE'
      units_ij(k) = '%'
      name_ij(k) = 'FRMP'
      ia_ij(k) = ia_rad
      scale_ij(k) = 100.

      IF (KOCEAN.eq.0) THEN
        k=k+1
        IJ_SMFX = k
        lname_ij(k) = 'SEA ICE IMPLICIT MASS FLUX'
        units_ij(k) = 'kg/m^2'
        name_ij(k) = 'SIMSFX'
        ia_ij(k) = ia_12hr
        scale_ij(k) = 2.
      END IF

      k=k+1 !
      IJ_CLR_SRINCG = k ! SRINCG*CLRSKY (W/m**2)            2 RD
      lname_ij(k) = 'CLR SKY INCIDENT SOLAR RADIATION, SRF x CLRSKY'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'CLR_SRINCG'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
c
      k=k+1 !
      IJ_CLDCV1 = k ! PCLD (1)  (COMPOSITE OVER ATMOSPHERE)   2 RD
      lname_ij(k) = 'TAU>1 CLOUD COVER'
      units_ij(k) = '%'
      name_ij(k) = 'CLDCV1'
      ia_ij(k) = ia_rad
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_CLDT1P  = k ! P-CLOUD TOP   (100 PA)                  2 RD
      lname_ij(k) = 'CLOUD TAU=1 PRESSURE x TAU>1 CLOUD COVER'
      units_ij(k) = 'mb'
      name_ij(k) = 'CLDT1P'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
cc    iw_ij(k) = iw_cldcv  ! built in
      ir_ij(k) = ir_0_1775
c
      k=k+1 !
      IJ_CLDTPT = k !
      lname_ij(k) = 'CLOUD TOP TEMPERATURE x TOTAL CLOUD COVER'
      units_ij(k) = 'C'
      name_ij(k) = 'CLDTPT'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
cc    iw_ij(k) = iw_cldcv  ! built in
      ir_ij(k) = ir_m80_28
c
      k=k+1 !
      IJ_CLDT1T = k !
      lname_ij(k) = 'CLOUD TAU=1 TEMPERATURE x TAU>1 CLOUD COVER'
      units_ij(k) = 'C'
      name_ij(k) = 'CLDT1T'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
cc    iw_ij(k) = iw_cldcv1 ! built in
      ir_ij(k) = ir_m80_28
c
      k=k+1 !
      IJ_WTRCLD = k ! PCLD (1)  (COMPOSITE OVER ATMOSPHERE)   2 RD
      lname_ij(k) = 'WATER CLOUD COVER'
      units_ij(k) = '%'
      name_ij(k) = 'WTRCLD'
      ia_ij(k) = ia_rad
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_ICECLD = k ! PCLD (1)  (COMPOSITE OVER ATMOSPHERE)   2 RD
      lname_ij(k) = 'ICE CLOUD COVER'
      units_ij(k) = '%'
      name_ij(k) = 'ICECLD'
      ia_ij(k) = ia_rad
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_OPTDW = k
      lname_ij(k) = 'WATER CLOUD OPTICAL DEPTH x WATER CLOUD COVER'
      units_ij(k) = ''
      name_ij(k) = 'optdw'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
c
      k=k+1 !
      IJ_OPTDI = k
      lname_ij(k) = 'ICE CLOUD OPTICAL DEPTH x ICE CLOUD COVER'
      units_ij(k) = ''
      name_ij(k) = 'optdi'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.

      if (k .gt. kaij) then
        write (6,*) 'ij_defs: Increase kaij=',kaij,' to at least ',k
        call stop_model( 'kaij too small', 255 )
      end if

      write (6,*) 'Number of AIJ diagnostics defined: kaijmax=',k
      do kk=1,k
        write (6,'(i4,'':'',a)') kk,trim(lname_ij(kk))
      end do
      return
      end subroutine ij_defs

      subroutine il_defs
      USE CONSTANT, only : grav,rgas,by3,sha,bygrav
      USE MODEL_COM, only : dtsrc,jeq
      USE GEOM, only : dxyp
      use DAGCOM
      implicit none
      real*8 :: bydj,bydjuv,daeq
      integer :: k,j,kk
c
      do k=1,kail
         write(name_il(k),'(a3,i3.3)') 'AIL',k
         lname_il(k) = 'unused'
         units_il(k) = 'unused'
         scale_il(k) = 1.
         ia_il(k)    = 0.
      enddo

C**** some scaling numbers for the equatorial diags.
      bydj   = 1./REAL(j5n-j5s+1,KIND=8)
      bydjuv = 1./REAL(j5nuv-j5suv+1,KIND=8)
      daeq=0.
      do j=j5s,j5n
        daeq=daeq+DXYP(J)
      end do
C****
      k=0
c
      k = k + 1
      IL_UEQ=k
      name_il(k) = 'u_equator'
      lname_il(k) = 'ZONAL WIND (U COMPONENT) AROUND +/- 5 DEG'
      units_il(k) = 'm/s'
      scale_il(k) = bydjuv
      ia_il(k)    = ia_dga
      k = k + 1
      IL_VEQ=k
      name_il(k) = 'v_equator'
      lname_il(k) = 'MERIDIONAL WIND (V COMPONENT) AROUND +/- 5 DEG'
      units_il(k) = 'm/s'
      scale_il(k) = bydjuv
      ia_il(k)    = ia_dga
      k = k + 1
      IL_WEQ=k
      name_il(k) = 'vvel_equator'
      lname_il(k) = 'VERTICAL VELOCITY AROUND +/- 5 DEG'
      units_il(k) = '10**-4 m/s'
      scale_il(k) = -1d4*RGAS*BYGRAV/daeq
      ia_il(k)    = ia_dga
      k = k + 1
      IL_TEQ=k
      name_il(k) = 'temp_equator'
      lname_il(k) = 'TEMPERATURE AROUND +/- 5 DEG'
      units_il(k) = 'C'
      scale_il(k) = bydj
      ia_il(k)    = ia_dga
      k = k + 1
      IL_QEQ=k
      name_il(k) = 'rh_equator'
      lname_il(k) = 'RELATIVE HUMIDITY AROUND +/- 5 DEG'
      units_il(k) = '%'
      scale_il(k) = 1d2*bydj
      ia_il(k)    = ia_dga
      k = k + 1
      IL_MCEQ=k
      name_il(k) = 'mcheat_equator'
      lname_il(k) = 'MOIST CONVECTIVE HEATING AROUND +/- 5 DEG'
      units_il(k) = '10**13 WATTS/DSIG'
      scale_il(k) = 100d-13*SHA/(GRAV*DTsrc)
      ia_il(k)    = ia_src
      k = k + 1
      IL_REQ=k
      name_il(k) = 'rad_cool_equator'
      lname_il(k) = 'TOTAL RADIATIVE COOLING AROUND +/- 5 DEG'
      units_il(k) = '10**13 WATTS/DSIG'
      scale_il(k) = -1d-13
      ia_il(k)    = ia_rad
      k = k + 1
      IL_W50N=k
      name_il(k) = 'vvel_50N'
      lname_il(k) = 'VERTICAL VELOCITY AT 50 N'
      units_il(k) = '10**-4 m/s'
      scale_il(k) = 1d4*RGAS/(GRAV*DXYP(J50N))
      ia_il(k)    = ia_dga
      k = k + 1
      IL_T50N=k
      name_il(k) = 'temp_50N'
      lname_il(k) = 'TEMPERATURE AT 50 N'
      units_il(k) = 'C'
      scale_il(k) = 1.
      ia_il(k)    = ia_dga
      k = k + 1
      IL_R50N=k
      name_il(k) = 'rad_cool_50N'
      lname_il(k) = 'TOTAL RADIATIVE COOLING AT 50 N'
      units_il(k) = '10**13 WATTS/UNIT SIGMA'
      scale_il(k) = 1d-13
      ia_il(k)    = ia_rad
      k = k + 1
      IL_U50N=k
      name_il(k) = 'u_50N'
      lname_il(k) = 'ZONAL WIND AT 50 N'
      units_il(k) = 'm/s'
      scale_il(k) = 0.5
      ia_il(k)    = ia_dga
      k = k + 1
      IL_W70N=k
      name_il(k) = 'vvel_70N'
      lname_il(k) = 'VERTICAL VELOCITY AT 70 N'
      units_il(k) = '10**-4 m/s'
      scale_il(k) = -1d4*RGAS/(GRAV*DXYP(J70N))
      ia_il(k)    = ia_dga
      k = k + 1
      IL_T70N=k
      name_il(k) = 'temp_70N'
      lname_il(k) = 'TEMPERATURE AT 70 N'
      units_il(k) = 'C'
      scale_il(k) = 1.
      ia_il(k)    = ia_dga
      k = k + 1
      IL_R70N=k
      name_il(k) = 'rad_cool_70N'
      lname_il(k) = 'TOTAL RADIATIVE COOLING AT 70 N'
      units_il(k) = '10**13 WATTS/UNIT SIGMA'
      scale_il(k) = -1d-13
      ia_il(k)    = ia_rad
      k = k + 1
      IL_U70N=k
      name_il(k) = 'u_70N'
      lname_il(k) = 'ZONAL WIND AT 70 N'
      units_il(k) = 'm/s'
      scale_il(k) = 0.5
      ia_il(k)    = ia_dga
c
      write (6,*) 'Number of IL diagnostics defined: kailmax=',k
      do kk=1,k
        write (6,'(i4,'':'',a)') kk,trim(lname_il(kk))
      end do
      return
      end subroutine il_defs


      subroutine jl_defs
      use CONSTANT, only : sday,grav,twopi,sha,rgas,bygrav,radius,lhe
      use MODEL_COM, only : fim,dtsrc,nidyn,byim
      use GEOM, only : dlon
      use DAGCOM
      implicit none
      integer :: k,kk
c
      do k=1,kajlx
         write(sname_jl(k),'(a3,i3.3)') 'AJL',k
         lname_jl(k) = 'unused'
         units_jl(k) = 'unused'
         pow_jl(k) = 0
      enddo
c
      k=0
c
      k=k+1
      jl_mcmflx = k
      sname_jl(k) = 'mc_mflx' !                  'FMX(MC)*P'
      lname_jl(k) = 'VERTICAL MASS EXCHANGE FROM MOIST CONVECTION'
      units_jl(k) = 'mb/s'
      scale_jl(k) = 1./(FIM*DTsrc)
      pow_jl(k) = -5
      ia_jl(k) = ia_src
      jgrid_jl(k) = 1
c
      k=k+1
      jl_srhr = k
      sname_jl(k) = 'srad_heat' !
      lname_jl(k) = 'SOLAR RADIATION HEATING RATE' !'SRHR'
      units_jl(k) = 'K/DAY' !'W/m^2'
      pow_jl(k) = -2
      scale_jl(k) = 1.D-2*GRAV*SDAY/SHA
      ia_jl(k) = ia_rad
      jgrid_jl(k) = 1
c
      k=k+1
      jl_trcr = k
      sname_jl(k) = 'trad_cool' !
      lname_jl(k) = 'THERMAL RADIATION COOLING RATE' !'TRHR'
      units_jl(k) = 'K/DAY' !'W/m^2'
      pow_jl(k) = -2
      scale_jl(k) = -1.D-2*GRAV*SDAY/SHA
      ia_jl(k) = ia_rad
      jgrid_jl(k) = 1
c
      k=k+1
      jl_sshr = k
      sname_jl(k) = 'lscond_heat' !
      lname_jl(k) = 'HEATING BY LARGE SCALE CONDENSATION' !'DTX(SS)*P'
      units_jl(k) = 'W/(m^2*mb)'
      pow_jl(k) = -2
      scale_jl(k) = 100.*BYGRAV*SHA/DTsrc
      ia_jl(k) = ia_src
      jgrid_jl(k) = 1
c
      k=k+1
      jl_trbhr = k
      sname_jl(k) = 'turb_heat' !
      lname_jl(k) = 'HEATING BY TURBULENCE' !'DT(DC)*P'
      units_jl(k) = 'W/(m^2*mb)'
      pow_jl(k) = -2
      scale_jl(k) = 100.*BYGRAV*SHA/DTsrc
      ia_jl(k) = ia_src
      jgrid_jl(k) = 1
c
      k=k+1
      jl_mchr = k
      sname_jl(k) = 'AJL13'
      lname_jl(k) = 'DT(MC)*P  DRY HEATING'
      units_jl(k) = '100 PA*K'
c
      k=k+1
      jl_ape = k
      sname_jl(k) = 'avail_pe' !
      lname_jl(k) = 'AVAILABLE POTENTIAL ENERGY'
      units_jl(k) = 'm^2/s^2'
      scale_jl(k) = .5*RGAS*BYIM
      ia_jl(k) = ia_dga
      jgrid_jl(k) = 1
c
      k=k+1
      jl_dtdyn = k
      sname_jl(k) = 'DT_DYNAMICS' !
      lname_jl(k) = 'DTEMP/DT BY DYNAMICS'
      units_jl(k) = 'K/DAY'
      pow_jl(k) = -1
      scale_jl(k) = SDAY*NIDYN/(FIM*7200.)
      ia_jl(k) = ia_dga
      jgrid_jl(k) = 1
c
      k=k+1
      jl_totcld = k
      sname_jl(k) = 'totcld' !
      lname_jl(k) = 'TOTAL CLOUD COVER' !'PCLD*P (TOTAL)'
      units_jl(k) = '%'
      scale_jl(k) = 100.*BYIM
      ia_jl(k) = ia_rad
      jgrid_jl(k) = 1
c
      k=k+1
      jl_gwFirst = k   ! The next consecutive 9 are Gravity Wave Diags
      jl_dumtndrg = k
      sname_jl(k) = 'dudt_mtndrg' !
      lname_jl(k) = 'DU/DT BY STRAT MTN DRAG'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./(FIM*DTsrc)
      ia_jl(k) = ia_src
      jgrid_jl(k) = 2
c
      k=k+1
      jl_dushrdrg = k
      sname_jl(k) = 'dudt_shrdrg'
      lname_jl(k) = 'DU/DT BY STRAT SHR DRAG'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./(FIM*DTsrc)
      ia_jl(k) = ia_src
      jgrid_jl(k) = 2
c
      k=k+1
      jl_dumcdrgm10 = k
      sname_jl(k) = 'dudt_mcdrgm10' !
      lname_jl(k) = 'DU/DT BY STRAT MC DRAG C=-10'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./(FIM*DTsrc)
      ia_jl(k) = ia_src
      jgrid_jl(k) = 2
c
      k=k+1
      jl_dumcdrgp10 = k
      sname_jl(k) = 'dudt_mcdrgp10' !
      lname_jl(k) = 'DU/DT BY STRAT MC DRAG C=+10'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./(FIM*DTsrc)
      ia_jl(k) = ia_src
      jgrid_jl(k) = 2
c
      k=k+1
      jl_dumcdrgm40 = k
      sname_jl(k) = 'dudt_mcdrgm40' !
      lname_jl(k) = 'DU/DT BY STRAT MC DRAG C=-40'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./(FIM*DTsrc)
      ia_jl(k) = ia_src
      jgrid_jl(k) = 2
c
      k=k+1
      jl_dumcdrgp40 = k
      sname_jl(k) = 'dudt_mcdrgp40' !
      lname_jl(k) = 'DU/DT BY STRAT MC DRAG C=+40'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./(FIM*DTsrc)
      ia_jl(k) = ia_src
      jgrid_jl(k) = 2
c
      k=k+1
      jl_dumcdrgm20 = k
      sname_jl(k) = 'dudt_mcdrgm20' !
      lname_jl(k) = 'DU/DT BY STRAT MC DRAG C=-20'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./(FIM*DTsrc)
      ia_jl(k) = ia_src
      jgrid_jl(k) = 2
c
      k=k+1
      jl_dumcdrgp20 = k
      sname_jl(k) = 'dudt_mcdrgp20' !
      lname_jl(k) = 'DU/DT BY STRAT MC DRAG C=+20'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./(FIM*DTsrc)
      ia_jl(k) = ia_src
      jgrid_jl(k) = 2
c Last of the Gravity Wave JL's
      k=k+1
      jl_dudfmdrg = k
      sname_jl(k) = 'dudt_dfmdrg' !
      lname_jl(k) = 'DU/DT BY STRAT DEFORM DRAG'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./(FIM*DTsrc)
      ia_jl(k) = ia_src
      jgrid_jl(k) = 2
c
      k=k+1
      jl_sscld = k
      sname_jl(k) = 'sscld' !
      lname_jl(k) = 'SUPER SATURATION CLOUD COVER' !'PCLD*P (SS)'
      units_jl(k) = '%'
      scale_jl(k) = 100.*BYIM
      ia_jl(k) = ia_rad
      jgrid_jl(k) = 1
c
      k=k+1
      jl_mccld = k
      sname_jl(k) = 'mccld' !
      lname_jl(k) = 'MOIST CONVECTIVE CLOUD COVER' !'PCLD*P (MC)'
      units_jl(k) = '%'
      scale_jl(k) = 100.*BYIM
      ia_jl(k) = ia_rad
      jgrid_jl(k) = 1
c
      k=k+1
      jl_sdifcoef = k
      sname_jl(k) = 'strat_diff_coeff' !
      lname_jl(k) = 'STRAT. DIFFUSION COEFF'
      units_jl(k) = 'm^2/s'
c
      k=k+1
      jl_dudtsdif = k
      sname_jl(k) = 'dudt_sdiff' !     ! gwdrag
      lname_jl(k) = 'DU/DT  STRATOSPHERIC DIFFUSION'
      units_jl(k) = 'm/s^2'
c
      k=k+1
      jl_dtdtsdrg = k
      sname_jl(k) = 'DT_SDRAG' !
      lname_jl(k) = 'DTEMP/DT BY STRATOSPHERIC DRAG'
      units_jl(k) = 'K/DAY'
      pow_jl(k) = -1
      scale_jl(k) = SDAY/(FIM*DTsrc)
      ia_jl(k) = ia_src
      jgrid_jl(k) = 1
c
      k=k+1
      jl_epflxv = k
      sname_jl(k) = 'epflx_vert' !
      lname_jl(k) = 'VERTICAL ELIASSEN-PALM FLUX'
      units_jl(k) = 'm^2/s^2'
      pow_jl(k) = -3
      scale_jl(k) = .125*100.*BYGRAV*BYIM
      ia_jl(k) = ia_dga
      jgrid_jl(k) = 1
c
      k=k+1
      jl_epflxn = k
      sname_jl(k) = 'epflx_north' !
      lname_jl(k) = 'NORTHWARD ELIASSEN-PALM FLUX'
      units_jl(k) = 'm^2/s^2'
      scale_jl(k) = 1.
      ia_jl(k) = ia_dga
      jgrid_jl(k) = 2
c
      k=k+1
      jl_damdc = k
      sname_jl(k) = 'dudt_dc' !
      lname_jl(k) = 'CHANGE OF U-WIND BY TURBULENCE'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./DTsrc
      ia_jl(k) = ia_src
      jgrid_jl(k) = 2
c
      k=k+1
      jl_dammc = k
      sname_jl(k) = 'dudt_mc' !          'DU(MC)*P (UV GRID)'
      lname_jl(k) = 'CHANGE OF U-WIND BY MOIST CONV'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./DTsrc
      ia_jl(k) = ia_src
      jgrid_jl(k) = 2
c
      k=k+1
      jl_40 = k
      sname_jl(k) = 'AJL40' !DU(ED)*P*(DTSURF*DSIG*ED/DZ**2)  (UV GRID)
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      jl_uepac = k
      sname_jl(k) = 'u_epac' !
      lname_jl(k) = 'U WIND AVERAGED OVER I,EAST PACIFIC'
      units_jl(k) = '10^-1 m/s'
      scale_jl(k) = .2E+1
      ia_jl(k) = ia_dga
      jgrid_jl(k) = 2
c
      k=k+1
      jl_vepac = k
      sname_jl(k) = 'v_epac' !
      lname_jl(k) = 'V WIND AVERAGED OVER EAST PACIFIC'
      units_jl(k) = '10^-1 m/s'
      scale_jl(k) = .2E+1
      ia_jl(k) = ia_dga
      jgrid_jl(k) = 2
c
      k=k+1
      jl_wepac = k
      sname_jl(k) = 'vvel_epac' !
      lname_jl(k) = 'VERTICAL VELOCITY FOR EAST PACIFIC'
      units_jl(k) = '10**-5 m/s'
      scale_jl(k) = -1.D5*RGAS/(5.*GRAV)
      ia_jl(k) = ia_dga
      jgrid_jl(k) = 1
c
      k=k+1
      jl_uwpac = k
      sname_jl(k) = 'u_wpac' !
      lname_jl(k) = 'U WIND AVERAGED OVER WEST PACIFIC'
      units_jl(k) = '10^-1 m/s'
      scale_jl(k) = .2E+1
      ia_jl(k) = ia_dga
      jgrid_jl(k) = 2
c
      k=k+1
      jl_vwpac = k
      sname_jl(k) = 'v_wpac' !
      lname_jl(k) = 'V WIND AVERAGED OVER WEST PACIFIC'
      units_jl(k) = 'm/s'
      scale_jl(k) = .2d0
      pow_jl(k) = -1
      ia_jl(k) = ia_dga
      jgrid_jl(k) = 2
c
      k=k+1
      jl_wwpac = k
      sname_jl(k) = 'vvel_wpac' !
      lname_jl(k) = 'VERTICAL VELOCITY FOR WEST PACIFIC'
      units_jl(k) = '10**-5 m/s'
      scale_jl(k) = -1.D5*RGAS/(5.*GRAV)
      ia_jl(k) = ia_dga
      jgrid_jl(k) = 1
c
      k=k+1
      jl_47 = k
      sname_jl(k) = 'AJL47' !V-V*  =D((V-VI)*(T-TI)/DTHDP)/DP
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      jl_zmfntmom = k
      sname_jl(k) = 'zmf_nt_mom' !
      lname_jl(k) = 'NORTH TRANS ZON. MOM. BY ZON. MEAN FLOW'
      units_jl(k) = 'unknown'
c
      k=k+1
      jl_totntmom = k
      sname_jl(k) = 'tot_nt_mom' !
      lname_jl(k) = 'TOTAL NORTH TRANS ZON. MOM.'
      units_jl(k) = 'unknown'
c
      k=k+1
      jl_mchphas = k
      sname_jl(k) = 'AJL50'
      lname_jl(k) = 'DT(MC)*P  CHANGE OF PHASE'
      units_jl(k) = '100 PA*K'
c
      k=k+1
      jl_mcdtotw = k
      sname_jl(k) = 'mc_del_tot_wat' !'CLHE*DQ(MC BEFORE COND)*P'
      lname_jl(k) = 'CHANGE IN TOTAL WATER BY MOIST CONV'
      units_jl(k) = '100 PA*K'
c
      k=k+1
      jl_dudtsdrg = k
      sname_jl(k) = 'dudt_sdrag' !
      lname_jl(k) = 'DU/DT BY SDRAG'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./(FIM*DTsrc)
      ia_jl(k) = ia_src
      jgrid_jl(k) = 2
c
      k=k+1
      jl_mcldht = k
      sname_jl(k) = 'moist_lat' !
      lname_jl(k) = 'CHANGE OF LATENT HEAT BY CTEI'
      units_jl(k) = 'W/(m^2*mb)'
      pow_jl(k) = -2
      scale_jl(k) = 100.*BYGRAV*SHA/DTsrc
      ia_jl(k) = ia_src
      jgrid_jl(k) = 1
c
      k=k+1
      jl_trbke = k
      sname_jl(k) = 'tke' !
      lname_jl(k) = 'TURBULENT KINETIC ENERGY'
      units_jl(k) = 'm^2/s^2'
      scale_jl(k) = 1.
      ia_jl(k) = ia_src
      jgrid_jl(k) = 1
c
      k=k+1
      jl_trbdlht = k
      sname_jl(k) = 'turb_lat' !
      lname_jl(k) = 'CHANGE OF LATENT HEAT BY TURBULENCE'
      units_jl(k) = 'W/(m^2*mb)'
      pow_jl(k) = -2
      scale_jl(k) = 100.*BYGRAV*SHA/DTsrc
      ia_jl(k) = ia_src
      jgrid_jl(k) = 1
c
      k=k+1
      jl_mcheat = k
      sname_jl(k) = 'tot_ht_mc' !
      lname_jl(k) = 'TOTAL HEATING BY MOIST CONVECTION (Q1)'
      units_jl(k) = 'W/(m^2*mb)'
      pow_jl(k) = -2
      scale_jl(k) = 100.*BYGRAV*SHA/DTsrc
      ia_jl(k) = ia_src
      jgrid_jl(k) = 1
c
      k=k+1
      jl_mcdry = k
      sname_jl(k) = 'tot_dry_mc' !
      lname_jl(k) = 'TOTAL DRYING BY MOIST CONVECTION (Q2)'
      units_jl(k) = 'W/(m^2*mb)'
      pow_jl(k) = -2
      scale_jl(k) = 100.*BYGRAV*SHA/DTsrc
      ia_jl(k) = ia_src
      jgrid_jl(k) = 1
c
      k=k+1
      jl_zmfntlh = k
      sname_jl(k) = 'jl_zmf_nt_lh'    ! used in DIAGJK but not printed
      lname_jl(k) = 'MEAN MERIDIONAL NORTHWARD TRANS. OF LATENT HEAT'
      units_jl(k) = 'W/mb'
      pow_jl(k) = 9
      scale_jl(k) = 100.*bygrav*LHE*XWON*fim/DTsrc
      ia_jl(k) = ia_src
      jgrid_jl(k) = 2
c
      k=k+1
      jl_totntlh = k
      sname_jl(k) = 'jl_tot_nt_lh'
      lname_jl(k) = 'TOTAL NORTHWARD TRANSPORT OF LATENT HEAT'
      units_jl(k) = 'W/mb'
      pow_jl(k) = 10
      scale_jl(k) = 100.*bygrav*LHE*XWON*fim/DTsrc
      ia_jl(k) = ia_src
      jgrid_jl(k) = 2
c
      k=k+1
      jl_zmfvtlh = k     ! used in DIAGJK but not printed
      sname_jl(k) = 'jl_zmf_vt_lh'
      lname_jl(k) = 'MEAN MERIDIONAL VERTICAL TRANS. OF LATENT HEAT'
      units_jl(k) = 'W/m^2'
      scale_jl(k) = 100.*BYGRAV*LHE*XWON*byim/DTsrc
      ia_jl(k) = ia_src
      jgrid_jl(k) = 1
c
      k=k+1
      jl_totvtlh = k
      sname_jl(k) = 'jl_tot_vt_lh'
      lname_jl(k) = 'TOTAL VERTICAL TRANSPORT OF LATENT HEAT'
      units_jl(k) = 'W/m^2'
      scale_jl(k) = 100.*BYGRAV*LHE*XWON*byim/DTsrc
      ia_jl(k) = ia_src
      jgrid_jl(k) = 1
c
      k=k+1
      jl_rhe= k
      sname_jl(k) = 'rhe'
      lname_jl(k) = 'EFFECTIVE RELATIVE HUMIDITY' ! from cloud scheme
      units_jl(k) = '%'
      scale_jl(k) = 100.*byim
      ia_jl(k) = ia_src
      jgrid_jl(k) = 1
c
      k=k+1
      jl_cldmc= k
      sname_jl(k) = 'cldmc'   ! no output
      lname_jl(k) = 'MOIST CONVECTIVE CLOUD FRACTION' !from cloud scheme
      units_jl(k) = '%'
      scale_jl(k) = 100.*byim
      ia_jl(k) = ia_src
      jgrid_jl(k) = 1
c
      k=k+1
      jl_cldss= k
      sname_jl(k) = 'cldss'  ! no output
      lname_jl(k) = 'LARGE-SCALE CLOUD FRACTION' ! from cloud scheme
      units_jl(k) = '%'
      scale_jl(k) = 100.*byim
      ia_jl(k) = ia_src
      jgrid_jl(k) = 1
c
      k=k+1
      jl_csizmc= k
      sname_jl(k) = 'csizmc'
      lname_jl(k) = 'MOIST CONVECTIVE EFFECTIVE CLOUD PARTICLE SIZE'
      units_jl(k) = 'micron'
      scale_jl(k) = byim
      ia_jl(k) = ia_src
      jgrid_jl(k) = 1
c
      k=k+1
      jl_csizss= k
      sname_jl(k) = 'csizss'
      lname_jl(k) = 'LARGE-SCALE EFFECTIVE CLOUD PARTICLE SIZE'
      units_jl(k) = 'micron'
      scale_jl(k) = byim
      ia_jl(k) = ia_src
      jgrid_jl(k) = 1
c
      k=k+1
      jl_wcld = k
      sname_jl(k) = 'wcld' !
      lname_jl(k) = 'WATER CLOUD COVER' !'PCLD*P'
      units_jl(k) = '%'
      scale_jl(k) = 100.*BYIM
      ia_jl(k) = ia_rad
      jgrid_jl(k) = 1
c
      k=k+1
      jl_icld = k
      sname_jl(k) = 'icld' !
      lname_jl(k) = 'ICE CLOUD COVER' !'PCLD*P'
      units_jl(k) = '%'
      scale_jl(k) = 100.*BYIM
      ia_jl(k) = ia_rad
      jgrid_jl(k) = 1
c
      k=k+1
      jl_wcod = k
      sname_jl(k) = 'wcod' !
      lname_jl(k) = 'WATER CLOUD OPTICAL DEPTH' ! od*wcldcv
      units_jl(k) = '/1000mb'
c     scale_jl(k) = 1000.*BYIM ! /dp  built in
c     ia_jl(k) = ia_rad
      jgrid_jl(k) = 1
c
      k=k+1
      jl_icod = k
      sname_jl(k) = 'icod' !
      lname_jl(k) = 'ICE CLOUD OPTICAL DEPTH'   ! od*icldcv
      units_jl(k) = '/1000mb'
c     scale_jl(k) = 1000.*BYIM ! /dp  built in
c     ia_jl(k) = ia_rad
      jgrid_jl(k) = 1
c
      k=k+1
      jl_wcsiz= k
      sname_jl(k) = 'wcsiz'
      lname_jl(k) = 'EFFECTIVE WATER CLOUD PARTICLE SIZE' ! SIZ*OPT.DPTH
      units_jl(k) = 'micron'
      scale_jl(k) = byim
      ia_jl(k) = ia_rad
      jgrid_jl(k) = 1
c
      k=k+1
      jl_icsiz= k
      sname_jl(k) = 'icsiz'
      lname_jl(k) = 'EFFECTIVE ICE CLOUD PARTICLE SIZE' ! SIZ*OPT.DPTH
      units_jl(k) = 'micron'
      scale_jl(k) = byim
      ia_jl(k) = ia_rad
      jgrid_jl(k) = 1

      if (k .gt. kajl) then
        write (6,*) 'jl_defs: Increase kajl=',kajl,' to at least ',k
        call stop_model( 'kajl too small', 255 )
      end if

      write (6,*) 'Number of AJL diagnostics defined: kajlmax=',k
      do kk=1,k
        write (6,'(i4,'':'',a)') kk,trim(lname_jl(kk))
      end do
      return
      end subroutine jl_defs

      subroutine sjl_defs
      use CONSTANT, only : grav,sday,sha,bygrav
      use MODEL_COM, only : byim
      use DAGCOM
      implicit none
      integer :: k,kk
c
      do k=1,kasjl
         write(name_sjl(k),'(a4,i3.3)') 'ASJL',k
         lname_sjl(k) = 'unused'
         units_sjl(k) = 'unused'
      enddo
c
      k=0
c
      k=k+1
      name_sjl(k) = 'ASJL01'
      lname_sjl(k) = 'TX'
      units_sjl(k) = 'C'
      scale_sjl(k) = BYIM
      ia_sjl(k) = ia_dga
c
      k=k+1
      name_sjl(k) = 'ASJL02'
      lname_sjl(k) = 'HEIGHT'
      units_sjl(k) = '100 M'
      scale_sjl(k) = BYIM*.01*BYGRAV
      ia_sjl(k) = ia_dga
c
      k=k+1
      name_sjl(k) = 'srad_heat' !'ASJL03'
      lname_sjl(k) = 'SOLAR RADIATION HEATING RATE' !'SRHR'
      units_sjl(k) = '10**-2 K/DAY' !'W/m^2'
      scale_sjl(k) = 100.D-2*GRAV*SDAY*BYIM/SHA
      ia_sjl(k) = ia_rad
c
      k=k+1
      name_sjl(k) = 'trad_cool' !'ASJL04'
      lname_sjl(k) = 'THERMAL RADIATION COOLING RATE' !'TRHR'
      units_sjl(k) = '10**-2 K/DAY' !'W/m^2'
      scale_sjl(k) = -100.D-2*GRAV*SDAY*BYIM/SHA
      ia_sjl(k) = ia_rad
c
      write (6,*) 'Number of ASJL diagnostics defined: kasjlmax=',k
      do kk=1,k
        write (6,'(i4,'':'',a)') kk,trim(lname_sjl(kk))
      end do
      return
      end subroutine sjl_defs

      subroutine jk_defs
      use CONSTANT, only : sday,twopi,rgas,lhe,bygrav,bymrat
      use MODEL_COM, only : fim,byim,dt
      use GEOM, only : dlon
      use DAGCOM
      implicit none
      integer :: k,kk
c
      do k=1,kajkx
         write(sname_jk(k),'(a3,i3.3)') 'AJK',k
         lname_jk(k) = 'unused'
         units_jk(k) = 'unused'
         pow_jk(k) = 0
      enddo
c
      k=0
c
      k=k+1
      jk_dpa = k             !'AJK01'
      sname_jk(k) = 'dp_cp1' !   DP=MIN(PM(K),PS)-MIN(PM(K+1),PS)
      lname_jk(k) =  'PRESSURE DIFFERENCES (CP,PT)' ! DP (PT GRID)
      units_jk(k) = 'mb'
      scale_jk(k) = byim
      ia_jk(k) = ia_dga
      jgrid_jk(k) = 1
c
      k=k+1
      jk_dpb = k
      sname_jk(k) = 'dp_cp2' !'AJK02'
      lname_jk(k) = 'PRESSURE DIFFERENCES (CP,UV)' ! DP (UV GRID)
      units_jk(k) = 'mb'
      scale_jk(k) = byim
      ia_jk(k) = ia_dga
      jgrid_jk(k) = 2
c
      k=k+1
      jk_temp = k
      sname_jk(k) = 'temp' !'AJK03'
      lname_jk(k) = 'TEMPERATURE' !'(TX-273.16)*DP'
      units_jk(k) = 'C'
      scale_jk(k) = 1.
      jgrid_jk(k) = 1
c
      k=k+1
      jk_hght = k
      sname_jk(k) = 'height' !'AJK04'
      lname_jk(k) = 'HEIGHT' !'PHI*DP'
      units_jk(k) = 'm'
      pow_jk(k) = 2
      scale_jk(k) = BYGRAV
      jgrid_jk(k) = 1
c
      k=k+1
      jk_q = k
      sname_jk(k) = 'q' !'AJK05'
      lname_jk(k) = 'SPECIFIC HUMIDITY' !'Q*DP'
      units_jk(k) = 'ppmv'
      pow_jk(k) = -1
      scale_jk(k) = 1.e6*bymrat
      jgrid_jk(k) = 1
c
      k=k+1
      jk_theta = k
      sname_jk(k) = 'pot_temp' !'AJK06'
      lname_jk(k) = 'POTENTIAL TEMPERATURE' !'TH*DP'
      units_jk(k) = 'K'
      scale_jk(k) = p1000k
      jgrid_jk(k) = 1
c
      k=k+1
      jk_rh = k
      sname_jk(k) = 'rh' !'AJK07'
      lname_jk(k) = 'RELATIVE HUMIDITY' !'RH*DP'
      units_jk(k) = '%'
      scale_jk(k) = 100.
      jgrid_jk(k) = 1
c
      k=k+1
      jk_u = k
      sname_jk(k) = 'u' !'AJK08'
      lname_jk(k) = 'ZONAL WIND (U COMPONENT)' !'U*DP4  (UV GRID)'
      units_jk(k) = 'm/s' !'100 PA*m/s'
      pow_jk(k) = -1
      scale_jk(k) = 1.
      jgrid_jk(k) = 2
c
      k=k+1
      jk_v = k
      sname_jk(k) = 'v' !'AJK09'
      lname_jk(k) = 'MERIDIONAL WIND (V COMPONENT)' !'V*DP4  (UV GRID)'
      units_jk(k) = 'm/s' !'100 PA*m/s'
      pow_jk(k) = -2
      scale_jk(k) = 1.
      jgrid_jk(k) = 2
c
      k=k+1
      jk_zmfke = k
      sname_jk(k) = 'zmf_ke' !'AJK10'
      lname_jk(k) = 'KINETIC ENERGY OF ZONAL MEAN FLOW'
      units_jk(k) = 'm^2/s^2'
c
      k=k+1
      jk_totke = k
      sname_jk(k) = 'tot_ke' !'AJK11'
      lname_jk(k) = 'TOTAL KINETIC ENERGY'
      units_jk(k) = 'm^2/s^2'
      scale_jk(k) = .5
      jgrid_jk(k) = 2
c
      k=k+1
      jk_zmfntsh = k
      sname_jk(k) = 'zmf_nt_sh' !'AJK12'
      lname_jk(k) = 'NORTH. TRANS. SENS. HT. BY ZON. MEAN FLOW'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_totntsh = k
      sname_jk(k) = 'tot_nt_sh' !'AJK13'
      lname_jk(k) = 'TOT NORTH. TRANS. SENS. HT'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_zmfntgeo = k
      sname_jk(k) = 'zmf_nt_geo' !'AJK14'
      lname_jk(k) = 'NORTH. TRANS. GEOPOT. BY ZON. MEAN FLOW'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_totntgeo = k
      sname_jk(k) = 'tot_nt_geo' !'AJK15'
      lname_jk(k) = 'TOT NORTH. TRANS. GEOPOT.'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_zmfntlh = k
      sname_jk(k) = 'zmf_nt_lh' !'AJK16'
      lname_jk(k) = 'NORTH TRANS LAT HT BY ZON. MEAN FLOW'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_totntlh = k
      sname_jk(k) = 'tot_nt_lh' !'AJK17'
      lname_jk(k) = 'TOTAL NORTHWARD TRANSPORT OF LATENT HEAT'
      units_jk(k) = 'W/mb'
      pow_jk(k) = 10
      scale_jk(k) = .25*LHE*XWON*FIM*100.*BYGRAV
      jgrid_jk(k) = 2
c
      k=k+1
      jk_zmfntke = k
      sname_jk(k) = 'zmf_nt_ke' !'AJK18'
      lname_jk(k) = 'NORTH TRANS KE BY ZON. MEAN FLOW'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_totntke = k
      sname_jk(k) = 'tot_nt_ke' !'AJK19'
      lname_jk(k) = 'TOTAL NORTHWARD TRANSPORT OF KINETIC ENERGY'
      units_jk(k) = 'W/mb'
      pow_jk(k) = 8
      scale_jk(k) = .5*XWON*FIM
      jgrid_jk(k) = 2
c
      k=k+1
      jk_zmfntmom = k
      sname_jk(k) = 'zmf_nt_mom' !'AJK20'
      lname_jk(k) = 'NORTH TRANS ZON. MOM. BY ZON. MEAN FLOW'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_totntmom = k
      sname_jk(k) = 'tot_nt_mom' !'AJK21'
      lname_jk(k) = 'TOT NORTH. TRANS. ZON. MOM.'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_p2kedpgf = k
      sname_jk(k) = 'p2k_eddy_pgf' !'AJK22'
      lname_jk(k) = 'P-K BY EDDY PRESSURE GRADIENT FORCE'
      units_jk(k) = 'W/(m^2*mb)'
      scale_jk(k) = .5*1d2*BYGRAV/DT
      jgrid_jk(k) = 2
      pow_jk(k) = -4
c
      k=k+1
      jk_dpsqr = k
      sname_jk(k) = 'dp_sqr' !'AJK23'
      lname_jk(k) = 'SQUARE OF PRESSURE DIFFERENCES'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_nptsavg = k
      sname_jk(k) = 'npts_avg' !'AJK24'
      lname_jk(k) = 'NUMBER OF GRIDPOINTS INCLUDED IN AVERAGE (CP)'
      units_jk(k) = '1'
      scale_jk(k) = TWOPI/(DLON*FIM)
      ia_jk(k) = ia_dga
      jgrid_jk(k) = 2
c
      k=k+1
      jk_vvel = k
      sname_jk(k) = 'vvel' !'AJK25'
      lname_jk(k) = 'VERTICAL VELOCITY (POSITIVE UPWARD)'
      units_jk(k) = 'mb/s'
      pow_jk(k) = -5
      scale_jk(k) = -byim
      ia_jk(k) = ia_dga
      jgrid_jk(k) = 1
c
      k=k+1
      jk_zmfvtdse = k
      sname_jk(k) = 'zmf_vt_dse' !'AJK26'
      lname_jk(k) = 'VERT TRANS DRY STAT. ENER. BY ZON. MEAN FLOW'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_totvtdse = k
      sname_jk(k) = 'tot_vt_dse' !'AJK27'
      lname_jk(k) = 'TOTAL LGE SCALE VERT. TRANS. OF DRY STAT ENRG (CP)'
      units_jk(k) = 'W/m^2'
      scale_jk(k) = -100.*BYGRAV*BYIM
      ia_jk(k) = ia_dga
      pow_jk(k) = 1
      jgrid_jk(k) = 1
c
      k=k+1
      jk_zmfvtlh = k
      sname_jk(k) = 'zmf_vt_lh' !'AJK28'
      lname_jk(k) = 'VERT TRANS LATENT HEAT BY ZON. MEAN FLOW'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_totvtlh = k
      sname_jk(k) = 'tot_vt_lh' !'AJK29'
      lname_jk(k) = 'TOTAL LGE SCALE VERT. TRANS. OF LATENT HEAT (CP)'
      units_jk(k) = 'W/m^2'
      scale_jk(k) = -100.*BYGRAV*BYIM*LHE
      ia_jk(k) = ia_dga
      jgrid_jk(k) = 1
c
      k=k+1
      jk_vtgeoeddy = k
      sname_jk(k) = 'vt_geopot_eddy' !'AJK30'
      lname_jk(k) = 'VERT. TRANS. OF GEOPOT. ENERGY BY EDDIES (CP)'
      units_jk(k) = 'W/m^2'
      scale_jk(k) = -100.*BYGRAV*BYIM
      ia_jk(k) = ia_dga
      pow_jk(k) = -2
      jgrid_jk(k) = 1
c
      k=k+1
      jk_barekegen = k
      sname_jk(k) = 'baroc_eddy_ke_gen' !'AJK31'
      lname_jk(k) = 'BAROCLINIC EDDY KINETIC ENERGY GEN.'
      units_jk(k) = 'W/(m^2*mb)'
      scale_jk(k) = .5*RGAS*1d2*BYGRAV
      pow_jk(k) = -4
      jgrid_jk(k) = 1
c
      k=k+1
      jk_potvort = k
      sname_jk(k) = 'pot_vort' !'AJK32'
      lname_jk(k) = 'POTENTIAL VORTICITY (CP)'
      units_jk(k) = 'K/(mb*s)'
      pow_jk(k) = -6
      scale_jk(k) = byim*p1000k
      ia_jk(k) = ia_dga
      jgrid_jk(k) = 1
c
      k=k+1
      jk_vtpv = k
      sname_jk(k) = 'vt_pv' !'AJK33'
      lname_jk(k) = 'VERT. TRANSPORT OF POTENTIAL VORTICITY (CP)'
      units_jk(k) = 'K/s^2'
      pow_jk(k) = -9
      scale_jk(k) = -.25*P1000K
      ia_jk(k) = ia_dga
      jgrid_jk(k) = 1
c
      k=k+1
      jk_vtpveddy = k
      sname_jk(k) = 'vt_pv_eddy' !'AJK34'
      lname_jk(k) = 'VERT. TRANS. OF POT. VORT. BY EDDIES (CP)'
      units_jk(k) = 'K/S^2'
      pow_jk(k) = -9
      scale_jk(k) = -.25*P1000K
      ia_jk(k) = ia_dga
      jgrid_jk(k) = 1
c
      k=k+1
      jk_nptsavg1 = k
      sname_jk(k) = 'npts_avg1' !'AJK35'
      lname_jk(k) = 'NUMBER OF GRIDPOINTS IN AVERAGE (PT GRID)'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_totvtke = k
      sname_jk(k) = 'tot_vt_ke' !'AJK36'
      lname_jk(k) ='TOTAL LGE SCALE VERT. TRANS. OF KINETIC ENRG (CP)'
      units_jk(k) = 'W/m^2'
      pow_jk(k) = -1
      scale_jk(k) = -.125*100.*BYGRAV*BYIM
      ia_jk(k) = ia_dga
      jgrid_jk(k) = 2
c
      k=k+1
      jk_vtameddy = k
      sname_jk(k) = 'vt_u_eddy' !'AJK37'
      lname_jk(k) = 'EDDY VERTICAL ZONAL MOM. FLUX (CP)'
      units_jk(k) = 'm^2/s^2'
      pow_jk(k) = -3
      scale_jk(k) = -.25*100.*BYGRAV*BYIM
      ia_jk(k) = ia_dga
      jgrid_jk(k) = 2
c
      k=k+1
      jk_totvtam = k
      sname_jk(k) = 'tot_vt_u' !'AJK38'
      lname_jk(k) = 'TOTAL VERTICAL ZONAL MOM. FLUX (CP)'
      units_jk(k) = 'm^2/s^2'
      pow_jk(k) = -2
      scale_jk(k) = -.25*100.*BYGRAV*BYIM
      ia_jk(k) = ia_dga
      jgrid_jk(k) = 2
c
      k=k+1
      jk_sheth = k
      sname_jk(k) = 'sheth' !'AJK39'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_dudtmadv = k
      sname_jk(k) = 'dudt_mean_advec' !'AJK40'
      lname_jk(k) = 'DU/DT BY MEAN ADVECTION (CP)'
      units_jk(k) = 'm/s^2'
      pow_jk(k) = -6
      scale_jk(k) = 1.
      ia_jk(k) = ia_dga
      jgrid_jk(k) = 2
c
      k=k+1
      jk_dtdtmadv = k
      sname_jk(k) = 'dtempdt_mean_advec' !'AJK41'
      lname_jk(k) = 'DTEMP/DT BY MEAN ADVECTION (CP)'
      units_jk(k) = 'K/DAY'
      pow_jk(k) = -1
      scale_jk(k) = SDAY
      ia_jk(k) = ia_dga
      jgrid_jk(k) = 1
c
      k=k+1
      jk_dudttem = k
      sname_jk(k) = 'dudt_advec_tem' !'AJK42'
      lname_jk(k) = 'DU/DT BY TRANSFORMED ADVECTION (CP)'
      units_jk(k) = 'm/s^2'
      pow_jk(k) = -6
      scale_jk(k) = 1.
      ia_jk(k) = ia_dga
      jgrid_jk(k) = 2
c
      k=k+1
      jk_dtdttem = k
      sname_jk(k) = 'dtempdt_advec_tem' !'AJK43'
      lname_jk(k) = 'DTEMP/DT BY TRANSFORMED ADVECTION (CP)'
      units_jk(k) = 'K/DAY'
      pow_jk(k) = -1
      scale_jk(k) = SDAY
      ia_jk(k) = ia_dga
      jgrid_jk(k) = 1
c
      k=k+1
      jk_epflxncp = k
      sname_jk(k) = 'epflx_north_cp' !'AJK44'
      lname_jk(k) = 'NORTHWARD COMP. OF ELIASSEN-PALM FLUX (CP)'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_epflxvcp = k
      sname_jk(k) = 'epflx_vert_cp' !'AJK45'
      lname_jk(k) = 'NORTHWARD COMP. OF ELIASSEN-PALM FLUX (CP)'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_uinst = k
      sname_jk(k) = 'u_inst' !'AJK46'
      lname_jk(k) = 'INSTANTANEOUS ZONAL AVERAGE OF ZONAL WIND'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_totdudt = k
      sname_jk(k) = 'tot_dudt' !'AJK47'
      lname_jk(k) = 'DU/DT   TOTAL CHANGE (CP)'
      units_jk(k) = 'm/s^2'
      pow_jk(k) = -6
      jgrid_jk(k) = 2
c
      k=k+1
      jk_tinst = k
      sname_jk(k) = 't_inst' !'AJK48'
      lname_jk(k) = 'INSTANTANEOUS ZONAL AVERAGE OF TEMPERATURE'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_totdtdt = k
      sname_jk(k) = 'dtempdt' !'AJK49'
      lname_jk(k) = 'DTEMP/DT   TOTAL CHANGE (CP)'
      units_jk(k) = 'K/DAY'
      pow_jk(k) = -1
      jgrid_jk(k) = 1
c
      k=k+1
      jk_eddvtpt = k
      sname_jk(k) = 'edd_vt_pt'
      lname_jk(k) = 'EDDY VERTICAL TRANSPORT OF POT. TEMP.'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_cldh2o = k
      sname_jk(k) = 'cldh2o' !'AJK51'
      lname_jk(k) = 'TOTAL CLOUD WATER CONTENT'
      units_jk(k) = 'kg/kg'
      pow_jk(k) = -6
      scale_jk(k) = 1.
      jgrid_jk(k) = 1
c
      write (6,*) 'Number of AJK diagnostics defined: kajkmax=',k
      do kk=1,k
        write (6,'(i4,'':'',a)') kk,trim(lname_jk(kk))
      end do
      return
      end subroutine jk_defs

      subroutine ijk_defs
      use CONSTANT, only : bygrav,tf
      use DAGCOM
      implicit none
      integer :: k,kk
c
      do k=1,kaijkx
         write(name_ijk(k),'(a4,i3.3)') 'AIJK',k
         lname_ijk(k) = 'unused'
         units_ijk(k) = 'unused'
         scale_ijk(k) = 1.
         jgrid_ijk(k) = 2
         off_ijk(k)   = 0.
      enddo
c
      k=0
c
      k=k+1
      IJK_U=k
      name_ijk(k) = 'u' !'UDPB'
      lname_ijk(k) = 'U-WIND'!            x delta p, b-grid'
      units_ijk(k) = 'm/s'
      scale_ijk(k) = 1.
      off_ijk(k)   = 0.
c
      k=k+1
      IJK_V=k
      name_ijk(k) = 'v' !'VDPB'
      lname_ijk(k) = 'V-WIND'!            x delta p, b-grid'
      units_ijk(k) = 'm/s'
      scale_ijk(k) = 1.
      off_ijk(k)   = 0.
c
      k=k+1
      IJK_DSE=k
      name_ijk(k) = 'dse' !'DSEDPB'
      lname_ijk(k) = 'DRY STAT. ENERGY'!  x delta p x 4, b-grid'
      units_ijk(k) = 'm^2/s^2'
      scale_ijk(k) = .25
      off_ijk(k)   = 0.
c
      k=k+1
      IJK_DP=k
      name_ijk(k) = 'dp' !'DPB'
      lname_ijk(k) = 'DELTA-P'!           b-grid'
      units_ijk(k) = '100 PA'
      scale_ijk(k) = 1.
      off_ijk(k)   = 0.
c
      k=k+1
      IJK_T=k
      name_ijk(k) = 't' !'TDPB'
      lname_ijk(k) = 'TEMPERATURE'!       x delta p x 4, b-grid'
      units_ijk(k) = 'C'
      scale_ijk(k) = 0.25
      off_ijk(k)   = -TF
c
      k=k+1
      IJK_Q=k
      name_ijk(k) = 'q' !'QDPB'
      lname_ijk(k) = 'SPECIFIC HUMIDITY'! x delta p x 4, b-grid'
      units_ijk(k) = 'kg/kg' !'10**-5'
      scale_ijk(k) = 0.25    !*1d5
      off_ijk(k)   = 0.
c
      write (6,*) 'Number of AIJK diagnostics defined: kaijkmax=',k
      do kk=1,k
        write (6,'(i4,'':'',a)') kk,trim(lname_ijk(kk))
      end do
      return
      end subroutine ijk_defs


      subroutine wave_defs
      use DAGCOM
      implicit none
      integer :: k,kk
c
      do k=1,kwp
         write(name_wave(k),'(a5,i3.3)') 'WAVEP',k
         lname_wave(k) = 'unused'
         units_wave(k) = 'unused'
      enddo
c
      k=0
c
      k=k+1
      name_wave(k) = 'U850_EQ'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'V850_EQ'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'U300_EQ'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'V300_EQ'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'U050_EQ'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'V050_EQ'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'Z922_50N'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'Z700_50N'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'Z500_50N'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'Z300_50N'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'Z100_50N'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'Z010_50N'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      write (6,*) 'Number of Wave diagnostics defined: kwavemax=',k
      do kk=1,k
        write (6,'(i4,'':'',a)') kk,trim(name_wave(kk))
      end do
      return
      end subroutine wave_defs


      subroutine pj_defs
      use DAGCOM
      implicit none
      integer :: k
c
      do k=1,kapj
         write(name_pj(k),'(a3,i3.3)') 'APJ',k
         lname_pj(k) = 'unused'
         units_pj(k) = 'unused'
      enddo
c
      k=0
c
      k=k+1
      name_pj(k) = 'APJ1'
      lname_pj(k) = 'P'
      units_pj(k) = '100 PA'
c
      k=k+1
      name_pj(k) = 'APJ2'
      lname_pj(k) = '4*P4I (UV GRID)'
      units_pj(k) = '100 PA'
c
      return
      end subroutine pj_defs

      subroutine tsf_defs
      use DAGCOM
      implicit none
      integer :: k
c
      do k=1,ktsf
         write(name_tsf(k),'(a6,i3.3)') 'TSFREZ',k
         lname_tsf(k) = 'unused'
         units_tsf(k) = 'unused'
      enddo
c
      k=0
c
      k=k+1
      name_tsf(k) = 'TSFREZ1'
      lname_tsf(k) = 'FIRST DAY OF GROWING SEASON'
      units_tsf(k) = 'JULIAN DAY'
      tf_day1 = k
c
      k=k+1
      name_tsf(k) = 'TSFREZ2'
      lname_tsf(k) = 'LAST DAY OF GROWING SEASON'
      units_tsf(k) = 'JULIAN DAY'
      tf_last = k
c
      k=k+1
      name_tsf(k) = 'LKICEON'
      lname_tsf(k) = 'LAST DAY OF ICE-FREE LAKE'
      units_tsf(k) = 'JULIAN DAY'
      tf_lkon = k
c
      k=k+1
      name_tsf(k) = 'LKICEOFF'
      lname_tsf(k) = 'LAST DAY OF ICED-UP LAKE'
      units_tsf(k) = 'JULIAN DAY'
      tf_lkoff = k
c
      return
      end subroutine tsf_defs

      subroutine diurn_defs
!@sum  diurn_defs definitions for diurnal diagnostic accumulated arrays
!@auth G. Schmidt
!@ver  1.0
      use CONSTANT, only : sha,rgas,twopi,sday,grav
      use MODEL_COM, only : dtsrc,nisurf
      use DAGCOM
      implicit none
      integer :: k,kk
c
      do k=1,NDIUVAR
         write(name_dd(k),'(a5,i3.3)') 'DIURN',k
         lname_dd(k) = 'unused'
         units_dd(k) = 'unused'
         scale_dd(k) = 1.
      enddo
c
      k=0
c
      k=k+1
      IDD_ISW=k
      name_dd(k)='INCSW'
      units_dd(k)='W/m^2'
      scale_dd(k)=1.
      lname_dd(k)=' INC SW RADIATION'
c
      k=k+1
      IDD_PALB=k
      name_dd(k)='PLALB'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' P ALBD '
c
      k=k+1
      IDD_GALB=k
      name_dd(k)='GRDALB'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' G ALBD '
c
      k=k+1
      IDD_ABSA=k
      name_dd(k)='ABSATM'
      units_dd(k)='W/m^2'
      scale_dd(k)=1.
      lname_dd(k)=' ABS ATM'
c
      k=k+1
      IDD_ECND=k
      name_dd(k)='ENRGCND'
      units_dd(k)='W/m^2'
      scale_dd(k)=100.*SHA/(GRAV*DTsrc)
      lname_dd(k)=' E CNDS '
c
      k=k+1
      IDD_SPR=k
      name_dd(k)='SRFPRS'
      units_dd(k)='mb'
      scale_dd(k)=1.
      lname_dd(k)=' SRF PRS'
c
      k=k+1
      IDD_PT5=k
      name_dd(k)='POTT5'
      units_dd(k)='K'
      scale_dd(k)=1.
      lname_dd(k)=' PT 5   '
c
      k=k+1
      IDD_PT4=k
      name_dd(k)='POTT4'
      units_dd(k)='K'
      scale_dd(k)=1.
      lname_dd(k)=' PT 4   '
c
      k=k+1
      IDD_PT3=k
      name_dd(k)='POTT3'
      units_dd(k)='K'
      scale_dd(k)=1.
      lname_dd(k)=' PT 3   '
c
      k=k+1
      IDD_PT2=k
      name_dd(k)='POTT2'
      units_dd(k)='K'
      scale_dd(k)=1.
      lname_dd(k)=' PT 2   '
c
      k=k+1
      IDD_PT1=k
      name_dd(k)='POTT1'
      units_dd(k)='K'
      scale_dd(k)=1.
      lname_dd(k)=' PT 1   '
c
      k=k+1
      IDD_TS=k
      name_dd(k)='TSURF'
      units_dd(k)='K'
      scale_dd(k)=1.
      lname_dd(k)=' TS     '
c
      k=k+1
      IDD_TG1=k
      name_dd(k)='TGRND'
      units_dd(k)='K'
      scale_dd(k)=1.
      lname_dd(k)=' TG1    '
c
      k=k+1
      IDD_Q5=k
      name_dd(k)='Q5'
      units_dd(k)='1d-5 kg/kg'
      scale_dd(k)=1d5
      lname_dd(k)=' Q 5    '
c
      k=k+1
      IDD_Q4=k
      name_dd(k)='Q4'
      units_dd(k)='1d-5 kg/kg'
      scale_dd(k)=1d5
      lname_dd(k)=' Q 4    '
c
      k=k+1
      IDD_Q3=k
      name_dd(k)='Q3'
      units_dd(k)='1d-5 kg/kg'
      scale_dd(k)=1d5
      lname_dd(k)=' Q 3    '
c
      k=k+1
      IDD_Q2=k
      name_dd(k)='Q2'
      units_dd(k)='1d-5 kg/kg'
      scale_dd(k)=1d5
      lname_dd(k)=' Q 2    '
c
      k=k+1
      IDD_Q1=k
      name_dd(k)='Q1'
      units_dd(k)='1d-5 kg/kg'
      scale_dd(k)=1d5
      lname_dd(k)=' Q 1    '
c
      k=k+1
      IDD_QS=k
      name_dd(k)='QSURF'
      units_dd(k)='1d-5 kg/kg'
      scale_dd(k)=1d5
      lname_dd(k)=' QS     '
c
      k=k+1
      IDD_QG=k
      name_dd(k)='QGRND'
      units_dd(k)='1d-5 kg/kg'
      scale_dd(k)=1d5
      lname_dd(k)=' QG     '
c
      k=k+1
      IDD_SWG=k
      name_dd(k)='SWGRND'
      units_dd(k)='W/m^2'
      scale_dd(k)=NISURF/DTsrc
      lname_dd(k)=' SW ON G'
c
      k=k+1
      IDD_LWG=k
      name_dd(k)='LWGRND'
      units_dd(k)='W/m^2'
      scale_dd(k)=NISURF/DTsrc
      lname_dd(k)=' LW AT G'
c
      k=k+1
      IDD_SH=k
      name_dd(k)='SENSHT'
      units_dd(k)='W/m^2'
      scale_dd(k)=NISURF/DTsrc
      lname_dd(k)=' SNSB HT'
c
      k=k+1
      IDD_LH=k
      name_dd(k)='LATHT'
      units_dd(k)='W/m^2'
      scale_dd(k)=NISURF/DTsrc
      lname_dd(k)=' LAT HT '
c
      k=k+1
      IDD_HZ0=k
      name_dd(k)='NETHT'
      units_dd(k)='W/m^2'
      scale_dd(k)=NISURF/DTsrc
      lname_dd(k)=' HEAT Z0'
c
      k=k+1
      IDD_UG=k
      name_dd(k)='UGEOS'
      units_dd(k)='0.1 m/s'
      scale_dd(k)=10.
      lname_dd(k)=' UG*10  '
c
      k=k+1
      IDD_VG=k
      name_dd(k)='VGEOS'
      units_dd(k)='0.1 m/s'
      scale_dd(k)=10.
      lname_dd(k)=' VG*10  '
c
      k=k+1
      IDD_WG=k
      name_dd(k)='WGEOS'
      units_dd(k)='0.1 m/s'
      scale_dd(k)=10.
      lname_dd(k)=' WG*10  '
c
      k=k+1
      IDD_US=k
      name_dd(k)='USURF'
      units_dd(k)='0.1 m/s'
      scale_dd(k)=10.
      lname_dd(k)=' US*10  '
c
      k=k+1
      IDD_VS=k
      name_dd(k)='VSURF'
      units_dd(k)='m/s'
      scale_dd(k)=10.
      lname_dd(k)=' VS*10  '
c
      k=k+1
      IDD_WS=k
      name_dd(k)='WSURF'
      units_dd(k)='m/s'
      scale_dd(k)=10.
      lname_dd(k)=' WS*10  '
c
      k=k+1
      IDD_CIA=k
      name_dd(k)='CRISOANG'
      units_dd(k)=''
      scale_dd(k)=360./TWOPI
      lname_dd(k)=' ALPHA0 '
c
      k=k+1
      IDD_RIS=k
      name_dd(k)='RIS1'
      units_dd(k)='0.01'
      scale_dd(k)=100.
      lname_dd(k)=' RIS1*E2'
c
      k=k+1
      IDD_RIG=k
      name_dd(k)='RIGS'
      units_dd(k)='0.01'
      scale_dd(k)=100.
      lname_dd(k)=' RIGS*E2'
c
      k=k+1
      IDD_CM=k
      name_dd(k)='CM'
      units_dd(k)='1d-4'
      scale_dd(k)=1d4
      lname_dd(k)=' CM*E4  '
c
      k=k+1
      IDD_CH=k
      name_dd(k)='CH'
      units_dd(k)='1d-4'
      scale_dd(k)=1d4
      lname_dd(k)=' CH*E4  '
c
      k=k+1
      IDD_CQ=k
      name_dd(k)='CQ'
      units_dd(k)='1d-4'
      scale_dd(k)=1d4
      lname_dd(k)=' CQ*E4  '
c
      k=k+1
      IDD_EDS=k
      name_dd(k)='EDS'
      units_dd(k)='0.1'
      scale_dd(k)=10.
      lname_dd(k)=' EDS1*10'
c
      k=k+1
      IDD_DBL=k
      name_dd(k)='DBL'
      units_dd(k)='m'
      scale_dd(k)=1.
      lname_dd(k)=' DBL    '
c
      k=k+1
      IDD_DCF=k
      name_dd(k)='DCFREQ'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' DC FREQ'
c
      k=k+1
      IDD_LDC=k
      name_dd(k)='LDC'
      units_dd(k)='0.1 L'
      scale_dd(k)=10.
      lname_dd(k)=' LDC*10 '
c
      k=k+1
      IDD_PR=k
      name_dd(k)='PREC'
      units_dd(k)='0.01 mm/day'
      scale_dd(k)=100.*100.*SDAY/(DTsrc*GRAV)
      lname_dd(k)=' PRC*100'  ! check scale
c
      k=k+1
      IDD_EV=k
      name_dd(k)='EVAP'
      units_dd(k)='0.01 mm/day'
      scale_dd(k)=100.*SDAY*NISURF/DTsrc
      lname_dd(k)=' EVP*100'
c
      k=k+1
      IDD_DMC=k
      name_dd(k)='DEEPMC'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' DEEP MC'
c
      k=k+1
      IDD_SMC=k
      name_dd(k)='SHLWMC'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' SHLW MC'
c
      k=k+1
      IDD_CL7=k
      name_dd(k)='CLD7'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' CLD 7  '
c
      k=k+1
      IDD_CL6=k
      name_dd(k)='CLD6'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' CLD 6  '
c
      k=k+1
      IDD_CL5=k
      name_dd(k)='CLD5'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' CLD 5  '
c
      k=k+1
      IDD_CL4=k
      name_dd(k)='CLD4'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' CLD 4  '
c
      k=k+1
      IDD_CL3=k
      name_dd(k)='CLD3'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' CLD 3  '
c
      k=k+1
      IDD_CL2=k
      name_dd(k)='CLD2'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' CLD 2  '
c
      k=k+1
      IDD_CL1=k
      name_dd(k)='CLD1'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' CLD 1  '
c
      k=k+1
      IDD_W=k
      name_dd(k)='VERTVEL'
      units_dd(k)='1d-5 m/s'
      scale_dd(k)=1.
      lname_dd(k)=' W TO-5 '
c
      k=k+1
      IDD_CCV=k
      name_dd(k)='CLDCOV'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' C COVER'
c
      k=k+1
      IDD_SSP=k
      name_dd(k)='SSPREC'
      units_dd(k)='0.01 mm/day'
      scale_dd(k)=100.*100.*SDAY/(DTsrc*GRAV)
      lname_dd(k)=' SSP*100'
c
      k=k+1
      IDD_MCP=k
      name_dd(k)='MCPREC'
      units_dd(k)='0.01 mm/day'
      scale_dd(k)=100.*100.*SDAY/(DTsrc*GRAV)
      lname_dd(k)=' MCP*100'

      write (6,*) 'Number of Diurn diagnostics defined: kaddmax=',k
      do kk=1,k
        write (6,'(i4,'':'',a)') kk,trim(lname_dd(kk))
      end do
      return
      end subroutine diurn_defs
