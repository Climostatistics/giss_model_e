#include "rundeck_opts.h"

#ifdef TRACERS_ON
!@sum  TRACERS: generic tracer routines used for all tracers
!@+    Routines included: 
!@+      Generic diags: set_generic_tracer_diags
!@+      Apply previously set sources: apply_tracer_sources 
!@+      Radioactive Decay: tdecay 
!@+      Gravitaional Settling: trgrav
!@+      Check routine: checktr
!@auth Jean Lerner/Gavin Schmidt
!@ver  1.0

      SUBROUTINE set_generic_tracer_diags
!@sum set_generic_tracer_diags init trace gas attributes and diagnostics
!@auth J. Lerner
!@calls sync_param
      USE CONSTANT, only: mair,sday
      USE MODEL_COM, only: dtsrc,nisurf
      USE DAGCOM, only: ia_src,ia_12hr,ir_log2
      USE TRACER_COM
      USE TRACER_DIAG_COM
      USE PARAM
      implicit none
      integer :: l,k,n
      character*20 sum_unit(ntm),inst_unit(ntm)   ! for conservation
      character*10, DIMENSION(NTM) :: CMR,CMRWT
      logical :: qcon(KTCON-1), T=.TRUE. , F=.FALSE.
      character*50 :: unit_string

C**** Get itime_tr0 from rundeck if it exists
      call sync_param("itime_tr0",itime_tr0,ntm)
C**** Get to_volume_MixRat from rundecks if it exists
      call sync_param("to_volume_MixRat",to_volume_MixRat,ntm)
#ifdef TRACERS_WATER
C**** Decide on water tracer conc. units from rundeck if it exists
      call sync_param("to_per_mil",to_per_mil,ntm)
#endif
#ifdef TRACERS_SPECIAL_O18
C**** set super saturation parameter for isotopes if needed
      call sync_param("supsatfac",supsatfac)
#endif

C**** Get factor to convert from mass mixing ratio to volume mr
      do n=1,ntm
        if (to_volume_MixRat(n) .eq.1) then
          MMR_to_VMR(n) = mair/tr_mm(n)
          cmr(n) = 'V/V air'
        else
          MMR_to_VMR(n) = 1.d0
          cmr(n) = 'kg/kg air'
        endif
#ifdef TRACERS_WATER
        if (to_per_mil(n) .eq.1) then
          cmrwt(n) = 'per mil'
          cmr(n) = 'per mil'     ! this overrides to_volume_ratio
          MMR_to_VMR(n) = 1.
        else
          cmrwt(n) = 'kg/m^2'
        end if
#endif
      end do
C****
C**** TAJLN(J,L,KQ,N)  (SUM OVER LONGITUDE AND TIME OF)
C****
C**** jlq_power Exponent associated with a physical process 
C****      (for printing tracers).   (ntm_power+jlq_power)
C**** jls_power Exponent associated with a source/sink (for printing)
C**** Defaults for JLN
      scale_jlq(:) = 1./DTsrc
      jgrid_jlq(:) = 1
      ia_jlq(:) = ia_src

C**** Tracer concentration
      do n=1,ntm
        k = 1        ! <<<<< Be sure to do this
        jlnt_conc = k
        sname_jln(k,n) = trim(trname(n))//'_CONCENTRATION' 
        lname_jln(k,n) = trim(trname(n))//' CONCENTRATION' 
        jlq_power(k) = 0.
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),cmr(n))
        scale_jlq(k) = 1.d0
        scale_jln(n) = MMR_to_VMR(n)
#ifdef TRACERS_WATER
        if (to_per_mil(n).gt.0) units_jln(k,n) = unit_string(0,cmr(n))
#endif

C**** Physical processes affecting tracers
C****   F (TOTAL NORTHWARD TRANSPORT OF TRACER MASS)  (kg)
        k = k + 1 
        jlnt_nt_tot = k
        sname_jln(k,n) = 'tr_nt_tot_'//trname(n)
        lname_jln(k,n) = 'TOTAL NORTHWARD TRANSPORT OF '//
     &     trim(trname(n))//' MASS'
        jlq_power(k) = 11.
        jgrid_jlq(k) = 2
C****   STM/SM (MEAN MERIDIONAL N.T. OF TRACER MASS)  (kg)
        k = k + 1 
        jlnt_nt_mm = k
        sname_jln(k,n) = 'tr_nt_mm_'//trname(n)
        lname_jln(k,n) = 'NORTHWARD TRANS. OF '//
     &     trim(trname(n))//' MASS BY MERIDIONAL CIRC.'
        jlq_power(k) = 10.
        jgrid_jlq(k) = 2
C****   F (TOTAL VERTICAL TRANSPORT OF TRACER MASS)  (kg)
        k = k + 1
        jlnt_vt_tot = k
        sname_jln(k,n) = 'tr_vt_tot_'//trname(n)
        lname_jln(k,n) = 'TOTAL VERTICAL TRANSPORT OF '//
     &     trim(trname(n))//' MASS'
        jlq_power(k) = 11.
C****   STM/SM (MEAN MERIDIONAL V.T. OF TRACER MASS)  (kg)
        k = k + 1
        jlnt_vt_mm = k
        sname_jln(k,n) = 'tr_vt_mm_'//trname(n)
        lname_jln(k,n) = 'VERTICAL TRANS. OF '//
     &     trim(trname(n))//' MASS BY MERIDIONAL CIRC.'
        jlq_power(k) = 10.
C****   TMBAR-TM (CHANGE OF TRACER MASS BY MOIST CONVEC)(kg)
        k = k + 1
        jlnt_mc = k
        sname_jln(k,n) = 'tr_mc_'//trname(n)
        lname_jln(k,n) = 'CHANGE OF '//
     &     trim(trname(n))//' MASS BY MOIST CONVECTION'
        jlq_power(k) = 10.
C****   TMBAR-TM (CHANGE OF TRACER MASS BY Large-scale CONDENSE)  (kg)
        k = k + 1
        jlnt_lscond = k
        sname_jln(k,n) = 'tr_lscond'//trname(n)
        lname_jln(k,n) ='CHANGE OF '//
     &     trim(trname(n))//' MASS BY LARGE-SCALE CONDENSE'
        jlq_power(k) = 10.
C****   TMBAR-TM (CHANGE OF TRACER MASS BY DRY CONVEC)  (kg)
        k = k + 1
        jlnt_turb = k
        sname_jln(k,n) = 'tr_turb_'//trname(n)
        lname_jln(k,n) = 'CHANGE OF '//
     &     trim(trname(n))//' MASS BY TURBULENCE/DRY CONVECTION'
        jlq_power(k) = 10.
      end do

      if (k.gt. ktajl) then
        write (6,*) 
     &   'tjl_defs: Increase ktajl=',ktajl,' to at least ',k
        call stop_model('ktajl too small',255)
      end if

C**** Construct UNITS string for output
      do n=1,ntm
      do k=2,ktajl
      units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),'kg/s')
      end do; end do

C**** CONTENTS OF TAIJLN(I,J,LM,N)  (SUM OVER TIME OF)
C****        TML (M*M * KG TRACER/KG AIR)
C**** Set defaults that are true for all tracers and layers
      ia_ijt    = ia_src
      ir_ijt(:) = ir_log2   !n
      ijtm_power(:) = ntm_power(:)+4   !n for integrated mass
      do n=1,ntm
        ijtc_power(n) = ntm_power(n)+1 !n for concentration
#ifdef TRACERS_WATER
        if (to_per_mil(n) .eq.1) ijtc_power(n) = 0
#endif
      end do
C**** Tracer concentrations (TAIJLN)
      do n=1,ntm
      do l=1,lm
        write(sname_ijt(l,n),'(a,i2.2)') trim(TRNAME(n))//'_L_',l 
        write(lname_ijt(l,n),'(a,i2)')   trim(TRNAME(n))//' L ',l
        units_ijt(l,n) = unit_string(ijtc_power(n),cmr(n))
        scale_ijt(l,n) = MMR_to_VMR(n)*10.**(-ijtc_power(n))
      end do
      end do

C**** AIJN
C****     1  TM (SUM OVER ALL LAYERS) (M*M * KG TRACER/KG AIR)
C****     2  TRS (SURFACE TRACER CONC.) (M*M * KG TRACER/KG AIR)
C****     3  TM (SUM OVER ALL LAYERS) (M*M * KG TRACER)
      do n=1,ntm
C**** Summation of mass over all layers
      k = 1        ! <<<<< Be sure to do this
      tij_mass = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_Total_Mass'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//' Total Mass'
        units_tij(k,n) = unit_string(ijtm_power(n),'kg/m^2')
        scale_tij(k,n) = 10.**(-ijtm_power(n))
C**** Average concentration over layers
      k = k+1
      tij_conc = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_Average'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//' Average'
        units_tij(k,n) = unit_string(ijtc_power(n),cmr(n))
        scale_tij(k,n) = MMR_to_VMR(n)*10.**(-ijtc_power(n))
C**** Surface concentration
      k = k+1
      tij_surf = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_At_Surface'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//' At Surface'
        units_tij(k,n) = unit_string(ijtc_power(n),cmr(n))
        scale_tij(k,n)=MMR_to_VMR(n)*10.**(-ijtc_power(n))/
     *                 REAL(NIsurf,KIND=8)
#ifdef TRACERS_WATER
C**** the following diagnostics are set assuming that the particular
C**** tracer exists in water. 
C**** Tracers in precipitation
      k = k+1
      tij_prec = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_in_prcp'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' in Precip'
        units_tij(k,n)=unit_string(ijtc_power(n),cmrwt(n))
        scale_tij(k,n)=10.**(-ijtc_power(n))
C**** Tracers in evaporation
      k = k+1
      tij_evap = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_in_evap'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' in Evaporation'
        units_tij(k,n)=unit_string(ijtc_power(n),cmrwt(n))
        scale_tij(k,n)=10.**(-ijtc_power(n))/REAL(NIsurf,KIND=8)
C**** Tracers in river runoff
      k = k+1
      tij_rvr = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_in_rvr'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' in River Outflow'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
          scale_tij(k,n)=1.
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,'kg/kg')
          scale_tij(k,n)=10.**(-ijtc_power(n)-3)
        end if
C**** Tracers in sea ice
      k = k+1
      tij_seaice = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_in_ice'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' in Sea Ice  x POICE'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
          scale_tij(k,n)=1.
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,cmrwt(n))
          scale_tij(k,n)=10.**(-ijtc_power(n)-3)
        end if
C**** Tracers conc. in ground component (ie. water or ice surfaces)
      k = k+1
      tij_grnd = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_at_Grnd'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' at Ground'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,'kg/kg wat')
        end if
        scale_tij(k,n)=10.**(-ijtc_power(n)-3)/REAL(NIsurf,KIND=8)
C**** Tracers conc. in lakes (layer 1)
      k = k+1
      tij_lk1 = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_Lake1'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' Lakes layer 1'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,'kg/kg wat')
        end if
        scale_tij(k,n)=10.**(-ijtc_power(n)-3)/REAL(NIsurf,KIND=8)
C**** Tracers conc. in lakes (layer 2)
      k = k+1
      tij_lk2 = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_Lake2'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' Lakes layer 2'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,'kg/kg wat')
        end if
        scale_tij(k,n)=10.**(-ijtc_power(n)-3)/REAL(NIsurf,KIND=8)
C**** Tracers conc. in soil water 
      k = k+1
      tij_soil = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_in_Soil'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' Soil Water'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,'kg/kg wat')
        end if
        scale_tij(k,n)=10.**(-ijtc_power(n)-3)/REAL(NIsurf,KIND=8)
#endif
      end do

      if (k .gt. ktaij) then
        write (6,*) 
     &   'tij_defs: Increase ktaij=',ktaij,' to at least ',k
        call stop_model('ktaij too small',255)
      end if

      END SUBROUTINE set_generic_tracer_diags


      SUBROUTINE apply_tracer_2Dsource(dtstep)
!@sum apply_tracer_2Dsource adds non-interactive surface sources 
!@+       to tracers
!@auth Jean Lerner/Gavin Schmidt
      USE MODEL_COM, only : jm
      USE GEOM, only : imaxj
      USE QUSDEF, only : mz,mzz
      USE TRACER_COM, only : ntm,trm,trmom,ntsurfsrc
      USE FLUXES, only : trsource,trflux1,trsrfflx
      USE TRACER_DIAG_COM, only : taijs,tajls,ijts_source,jls_source
     *     ,itcon_surf
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: dtstep
      INTEGER n,ns,naij,najl,j,i
      REAL*8, DIMENSION(JM) :: dtracer
      REAL*8 ftr1

C**** This is tracer independent coding designed to work for all
C**** surface sources.
C**** Note that tracer flux is added to first layer either implicitly
C**** in ATURB or explcitly in 'apply_fluxes_to_atm' call in SURFACE.

      do n=1,ntm

        trflux1(:,:,n) = 0.
C**** Non-interactive sources
        do ns=1,ntsurfsrc(n)
C**** diagnostics
          naij = ijts_source(ns,n)
          taijs(:,:,naij) = taijs(:,:,naij) + trsource(:,:,ns,n)*dtstep
          najl = jls_source(ns,n)
          do j=1,jm
            tajls(j,1,najl) = tajls(j,1,najl)+sum(trsource(:,j,ns,n))
     *           *dtstep
            dtracer(j)=0.
            do i=1,imaxj(j)
              dtracer(j)=dtracer(j)+trsource(i,j,ns,n)*dtstep
            end do
          end do
          call DIAGTCB(dtracer,itcon_surf(ns,n),n)
C**** trflux1 is total flux into first layer
          trflux1(:,:,n) = trflux1(:,:,n)+trsource(:,:,ns,n)
        end do
C**** Interactive sources
C**** diagnostics
c        naij = ijts_source(ns,n)  ????
c        taijs(:,:,naij) = taijs(:,:,naij) + trsrfflx(:,:,n)*dtstep
c        najl = jls_source(ns,n)   ????
c        do j=1,jm
c          tajls(j,1,najl) = tajls(j,1,najl)+sum(trsrfflx(:,j,n))
c     *         *dtstep
c        end do
c        call DIAGTCA(itcon_surf(ns,n),n)  ????

C**** modify vertical moments (only from non-interactive sources)
        trmom( mz,:,:,1,n) = trmom( mz,:,:,1,n)-1.5*trflux1(:,:,n)
     *       *dtstep
        trmom(mzz,:,:,1,n) = trmom(mzz,:,:,1,n)+0.5*trflux1(:,:,n)
     *       *dtstep

C**** Accumulate interactive sources as well
        trflux1(:,:,n) = trflux1(:,:,n)+trsrfflx(:,:,n)

C**** Technically speaking the vertical moments should be modified here 
C**** as well. But for consistency with water vapour we only modify
C**** moments for dew.
        do j=1,jm
          do i=1,imaxj(j)
            if (trsrfflx(i,j,n).lt.0 .and. trm(i,j,1,n).gt.0) then
              ftr1=-trsrfflx(i,j,n)*dtstep/trm(i,j,1,n)
              trmom(:,i,j,1,n)=trmom(:,i,j,1,n)*(1.-ftr1)
            end if
          end do
        end do
      end do
C****
      RETURN
      END SUBROUTINE apply_tracer_2Dsource


      SUBROUTINE apply_tracer_3Dsource(ns,n)
!@sum apply_tracer_3Dsource adds 3D sources to tracers
!@auth Jean Lerner/Gavin Schmidt
      USE CONSTANT, only : teeny
      USE MODEL_COM, only : jm,im,lm,dtsrc
      USE GEOM, only : imaxj
      USE QUSDEF, only: nmom
      USE TRACER_COM, only : ntm,trm,trmom
      USE FLUXES, only : tr3Dsource
      USE TRACER_DIAG_COM, only : tajls,jls_3Dsource,itcon_3Dsrc
      IMPLICIT NONE
      REAL*8 fr3d
      INTEGER n,ns,najl,i,j,l

C**** This is tracer independent coding designed to work for all
C**** 3D sources.
C**** Modify tracer amount, moments, and diagnostics
C**** tracer moments are modified elsewhere
      najl = jls_3Dsource(ns,n)
      do l=1,lm
      do j=1,jm
      do i=1,imaxj(j)
C**** calculate fractional loss
        if (tr3Dsource(i,j,l,ns,n).lt.0.) then
          fr3d = -tr3Dsource(i,j,l,ns,n)*dtsrc/(trm(i,j,l,n)+teeny)
          trmom(1:nmom,i,j,l,n) = trmom(1:nmom,i,j,l,n)*(1.-fr3d)
        end if
C**** update tracer mass and diagnostics
        trm(i,j,l,n) = trm(i,j,l,n)+tr3Dsource(i,j,l,ns,n)*dtsrc
        tajls(j,l,najl) = tajls(j,l,najl)+
     *                              tr3Dsource(i,j,l,ns,n)*dtsrc
      end do; end do; end do
      call DIAGTCA(itcon_3Dsrc(ns,n),n)
C****
      RETURN
      END SUBROUTINE apply_tracer_3Dsource


      SUBROUTINE TDECAY
!@sum TDECAY decays radioactive tracers every source time step
!@auth Gavin Schmidt/Jean Lerner
      USE MODEL_COM, only : im,jm,lm,itime,dtsrc
      USE TRACER_COM, only : ntm,trm,trmom,trdecay,itime_tr0
#ifdef TRACERS_WATER
     *     ,trwm
      USE SEAICE_COM, only : trsi
      USE LAKES_COM, only : trlake
      USE LANDICE_COM, only : trlndi,trsnowli
      USE GHYCOM, only : tr_wbare,tr_wvege,tr_wsn_ij
#endif
      USE TRACER_DIAG_COM, only : tajls,jls_decay,itcon_decay
      IMPLICIT NONE
      real*8, save, dimension(ntm) :: expdec = 1.
      real*8, dimension(im,jm,lm) :: told
      logical, save :: ifirst=.true.
      integer n,najl,j,l

      if (ifirst) then               
        do n=1,ntm
          if (trdecay(n).gt.0.0) expdec(n)=exp(-trdecay(n)*dtsrc)
        end do
        ifirst = .false.
      end if

      do n=1,ntm
        if (trdecay(n).gt.0. .and. itime.ge.itime_tr0(n)) then
C**** Atmospheric decay
          told(:,:,:)=trm(:,:,:,n)
#ifdef TRACERS_WATER
     *               +trwm(:,:,:,n)
          trwm(:,:,:,n)   = expdec(n)*trwm(:,:,:,n)
#endif
          trm(:,:,:,n)    = expdec(n)*trm(:,:,:,n)
          trmom(:,:,:,:,n)= expdec(n)*trmom(:,:,:,:,n)
#ifdef TRACERS_WATER
C**** Note that ocean tracers are dealt with by separate ocean code.
C**** Decay sea ice tracers
          trsi(n,:,:,:)   = expdec(n)*trsi(n,:,:,:)
C**** ...lake tracers
          trlake(n,:,:,:) = expdec(n)*trlake(n,:,:,:)
C**** ...land surface tracers
          tr_wbare(n,:,:,:) = expdec(n)*tr_wbare(n,:,:,:)
          tr_wvege(n,:,:,:) = expdec(n)*tr_wvege(n,:,:,:)
          tr_wsn_ij(n,:,:,:,:)= expdec(n)*tr_wsn_ij(n,:,:,:,:)
          trsnowli(n,:,:) = expdec(n)*trsnowli(n,:,:)
          trlndi(n,:,:)   = expdec(n)*trlndi(n,:,:)
#endif
C**** atmospheric diagnostics
          najl = jls_decay(n)
          do l=1,lm
          do j=1,jm
            tajls(j,l,najl)=tajls(j,l,najl)+sum(trm(:,j,l,n)
#ifdef TRACERS_WATER
     *           +trwm(:,j,l,n)
#endif
     *           -told(:,j,l))
          enddo 
          enddo
          call DIAGTCA(itcon_decay(n),n)
        end if
      end do
C****
      return
      end subroutine tdecay


      SUBROUTINE TRGRAV
!@sum TRGRAV gravitationally settles particular tracers 
!@auth Gavin Schmidt/Reha Cakmur
      USE CONSTANT, only : visc_air,grav
      USE MODEL_COM, only : im,jm,lm,itime,dtsrc,zatmo
      USE GEOM, only : imaxj
      USE SOMTQ_COM, only : mz,mzz,mzx,myz,zmoms
      USE DYNAMICS, only : gz
      USE TRACER_COM, only : ntm,trm,trmom,itime_tr0,trradius,trpdens
      USE TRACER_DIAG_COM, only : tajls,jls_grav,itcon_grav
      USE FLUXES, only : trgrdep
      IMPLICIT NONE
      real*8, save, dimension(ntm) :: stokevdt = 0.
      logical, save :: ifirst=.true.
      real*8, dimension(im,jm,lm) :: told
      real*8 fgrfluxd,fgrfluxu
      integer n,najl,i,j,l

      if (ifirst) then               
C**** Calculate settling velocity based on Stokes' Law using particle
C**** density and effective radius
        do n=1,ntm
          if (trradius(n).gt.0.0) stokevdt(n)=dtsrc*2.*grav*trpdens(n)
     *         *trradius(n)**2/(9.*visc_air)
        end do
        ifirst = .false.
      end if

      do n=1,ntm
        if (trradius(n).gt.0. .and. itime.ge.itime_tr0(n)) then
C**** Gravitional settling 
          do l=1,lm
          do j=1,jm
          do i=1,imaxj(j)
            told(i,j,l)=trm(i,j,l,n)
C**** Calculate height differences using geopotential
            if (l.eq.1) then   ! layer 1 calc
              fgrfluxd=stokevdt(n)*grav/(gz(i,j,l)-zatmo(i,j))
              trgrdep(i,j,n)=fgrfluxd*trm(i,j,l,n)
            else               ! above layer 1
              fgrfluxd=stokevdt(n)*grav/(gz(i,j,l)-gz(i,j,l-1))
            end if
            if (l.lt.lm) then  ! below top layer
              fgrfluxu=stokevdt(n)*grav/(gz(i,j,l+1)-gz(i,j,l))
            else               ! top layer
              fgrfluxu=0.
            end if
            trm(i,j,l,n)=trm(i,j,l  ,n)*(1.-fgrfluxd)
     *                 + trm(i,j,l+1,n)*    fgrfluxu
            trmom(zmoms,i,j,l,n)=trmom(zmoms,i,j,l,n)*(1.-fgrfluxd)
          end do
          end do
          end do

          najl = jls_grav(n)
          do l=1,lm
          do j=1,jm
          tajls(j,l,najl)=tajls(j,l,najl)+sum(trm(:,j,l,n)-told(:,j,l))
          enddo 
          enddo
          call DIAGTCA(itcon_grav(n),n)
        end if
      end do
C****
      return
      end subroutine trgrav


      SUBROUTINE read_monthly_sources(iu,jdlast,tlca,tlcb,data,
     *  frac,imon)
!@sum Read in monthly sources and interpolate to current day
!@+   Calling routine must have the lines:
!@+      real*8 tlca(im,jm,nm),tlcb(im,jm,nm)
!@+      integer imon(nm)   ! nm=number of files that will be read
!@+      data jdlast /0/
!@+      save jdlast,tlca,tlcb,imon
!@+   Input: iu, the fileUnit#; jdlast
!@+   Output: interpolated data array + two monthly data arrays
!@auth Jean Lerner and others
      USE MODEL_COM, only: jday,im,jm,idofm=>JDmidOfM
      implicit none
      real*8 tlca(im,jm),tlcb(im,jm),data(im,jm),frac
      integer imon,iu,jdlast

      if (jdlast.EQ.0) then ! NEED TO READ IN FIRST MONTH OF DATA
        imon=1          ! imon=January
        if (jday.le.16)  then ! JDAY in Jan 1-15, first month is Dec
          call readt(iu,0,tlca,im*jm,tlca,12)
          rewind iu
        else            ! JDAY is in Jan 16 to Dec 16, get first month
  120     imon=imon+1
          if (jday.gt.idofm(imon) .AND. imon.le.12) go to 120
          call readt(iu,0,tlca,im*jm,tlca,imon-1)
          if (imon.eq.13)  rewind iu
        end if
      else              ! Do we need to read in second month?
        if (jday.ne.jdlast+1) then ! Check that data is read in daily
          if (jday.ne.1 .OR. jdlast.ne.365) then
            write(6,*)
     *      'Incorrect values in Tracer Source:JDAY,JDLAST=',JDAY,JDLAST
            call stop_model('stopped in TRACERS.f',255)
          end if
          imon=imon-12  ! New year
          go to 130
        end if
        if (jday.le.idofm(imon)) go to 130
        imon=imon+1     ! read in new month of data
        tlca(:,:) = tlcb(:,:)
        if (imon.eq.13) rewind iu
      end if
      call readt(iu,0,tlcb,im*jm,tlcb,1)
  130 continue
c**** Interpolate two months of data to current day
      frac = float(idofm(imon)-jday)/(idofm(imon)-idofm(imon-1))
      data(:,:) = tlca(:,:)*frac + tlcb(:,:)*(1.-frac)
      return
      end subroutine read_monthly_sources


      subroutine checktr(subr)
!@sum  CHECKTR Checks whether atmos tracer variables are reasonable
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : teeny
      USE MODEL_COM, only : ls1,im,jm,lm,q,wm
      USE GEOM, only : dxyp,imaxj
      USE SOMTQ_COM, only : qmom
      USE DYNAMICS, only : am
      USE FLUXES, only : gtracer
      USE TRACER_COM
      IMPLICIT NONE
      LOGICAL QCHECKT
      INTEGER I,J,L,N,m, imax,jmax,lmax
      REAL*8 relerr, errmax,errsc
!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR

      CALL CHECK3(gtracer(1,1,1,1),NTM,4,IM*JM,SUBR,'GTRACE')
      do n=1,ntm
        CALL CHECK3(trmom(1,1,1,1,n),NMOM,IM,JM*LM,SUBR,
     *       'X'//trname(n))
        CALL CHECK3(trm(1,1,1,n),IM,JM,LM,SUBR,trname(n))
#ifdef TRACERS_WATER
        CALL CHECK3(trwm(1,1,1,n),IM,JM,LM,SUBR,'WM'//trname(n))
#endif

C**** check for negative tracer amounts (if t_qlimit is set)
        if (t_qlimit(n)) then
          QCHECKT=.false.
          do l=1,lm
          do j=1,jm
          do i=1,imaxj(j)
            if (trm(i,j,l,n).lt.0) then
              write(6,*) "Negative mass for ",trname(n),i,j,l,trm(i,j,l
     *             ,n)," after ",SUBR,"." 
              QCHECKT=.true.
            end if
          end do
          end do
          end do
          if (QCHECKT)
     &         call stop_model("CHECKTR: Negative tracer amount",255)
        end if

C**** check whether air mass is conserved
        
        if (trname(n).eq.'Air') then
          errmax = 0. ; lmax=1 ; imax=1 ; jmax=1 
          do l=1,lm
          do j=1,jm
          do i=1,imaxj(j)
            relerr=abs(trm(i,j,l,n)-am(l,i,j)*dxyp(j))/
     *           (am(l,i,j)*dxyp(j))
            if (relerr.gt.errmax) then
              lmax=l ; imax=i ; jmax=j ; errmax=relerr
            end if
          end do
          end do
          end do
          print*,"Relative error in air mass after ",trim(subr),":",imax
     *         ,jmax,lmax,errmax,trm(imax,jmax,lmax,n),am(lmax,imax
     *         ,jmax)*dxyp(jmax)
        end if
      
#ifdef TRACERS_WATER
        if (trname(n).eq.'Water') then
          errmax = 0. ; lmax=1 ; imax=1 ; jmax=1
          do l=1,lm
          do j=1,jm
          do i=1,imaxj(j)
            errsc=(q(i,j,l)+sum(abs(qmom(:,i,j,l))))*am(l,i,j)*dxyp(j)
            if (errsc.eq.0.) errsc=1.
            relerr=abs(trm(i,j,l,n)-q(i,j,l)*am(l,i,j)*dxyp(j))/errsc
            if (wm(i,j,l).gt.0 .and. trwm(i,j,l,n).gt.1.) relerr
     *           =max(relerr,(trwm(i,j,l,n)-wm(i,j,l)*am(l,i,j)*dxyp(j))
     *           /(wm(i,j,l)*am(l,i,j)*dxyp(j)))
            if ((wm(i,j,l).eq.0 .and.trwm(i,j,l,n).gt.1) .or. (wm(i,j,l)
     *           .gt.teeny .and.trwm(i,j,l,n).eq.0))
     *           print*,"Liquid water mismatch: ",subr,i,j,l,trwm(i,j,l
     *           ,n),wm(i,j,l)*am(l,i,j)*dxyp(j) 
            do m=1,nmom
              relerr=max(relerr,(trmom(m,i,j,l,n)-qmom(m,i,j,l)*am(l,i,j
     *             )*dxyp(j))/errsc)
            end do
            if (relerr.gt.errmax) then
              lmax=l ; imax=i ; jmax=j ; errmax=relerr
            end if
          end do
          end do
          end do
          print*,"Relative error in water mass after ",trim(subr),":"
     *         ,imax,jmax,lmax,errmax,trm(imax,jmax,lmax,n),q(imax,jmax
     *         ,lmax)*am(lmax,imax,jmax)*dxyp(jmax),trwm(imax,jmax,lmax
     *         ,n),wm(imax,jmax,lmax)*am(lmax,imax,jmax)*dxyp(jmax)
     *         ,(trmom(m,imax,jmax,lmax,n),qmom(m,imax,jmax,lmax)
     *         *am(lmax,imax,jmax)*dxyp(jmax),m=1,nmom)
        end if
#endif
      end do

      return
      end subroutine checktr

#endif


      SUBROUTINE io_tracer(kunit,iaction,ioerr)
!@sum  io_tracer reads and writes tracer variables to file
!@auth Jean Lerner
!@ver  1.0
#ifdef TRACERS_ON      
      USE MODEL_COM, only: ioread,iowrite,irsfic,irsficno,irerun,lhead
      USE TRACER_COM
#ifdef TRACERS_SPECIAL_Shindell
      USE TRCHEM_Shindell_COM, only: yNO3,pHOx,pNOx,pOx,yCH3O2,
     & yC2O3,yROR,yXO2,yAldehyde,yXO2N,yRXPAR,OxIC,corrOx
#endif
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
#ifdef TRACERS_WATER
      CHARACTER*80 :: HEADER, MODULE_HEADER = "TRACERW01"

      write (MODULE_HEADER(lhead+1:80),'(a,i2,a,a,i1,a,i2,a,i2,a)')
     *     'R8 TRM(im,jm,lm,',NTM,')',
     *     ',TRmom(',NMOM,',im,jm,lm,',NTM,'),trwm(im,jm,lm,',NTM,')' 
#else
      CHARACTER*80 :: HEADER, MODULE_HEADER = "TRACER01"

      write (MODULE_HEADER(lhead+1:80),'(a,i2,a,a,i1,a,i2,a)')
     *           'R8 TRM(im,jm,lm,',NTM,')',
     *  ',TRmom(',NMOM,',im,jm,lm,',NTM,')'

#endif

      SELECT CASE (IACTION)

      CASE (:IOWRITE) ! output to end-of-month restart file
        WRITE (kunit,err=10) MODULE_HEADER,TRM,TRmom
#ifdef TRACERS_WATER
     *     ,TRWM
#endif
#ifdef TRACERS_SPECIAL_Shindell
     *     ,yNO3,pHOx,pNOx,pOx,yCH3O2,yC2O3,yROR,yXO2
     *     ,yAldehyde,yXO2N,yRXPAR,OxIC,corrOx
#endif
      CASE (IOREAD:)          ! input from restart file
        SELECT CASE (IACTION)
        CASE (ioread,irerun,irsfic,irsficno) ! restarts
          READ (kunit,err=10) HEADER,TRM,TRmom
#ifdef TRACERS_WATER
     *       ,TRWM
#endif
#ifdef TRACERS_SPECIAL_Shindell
     *     ,yNO3,pHOx,pNOx,pOx,yCH3O2,yC2O3,yROR,yXO2
     *     ,yAldehyde,yXO2N,yRXPAR,OxIC,corrOx
#endif
          IF (HEADER(1:lhead).ne.MODULE_HEADER(1:lhead)) THEN
            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
            GO TO 10
          END IF
        END SELECT
      END SELECT

      RETURN
 10   IOERR=1
#endif
      RETURN
      END SUBROUTINE io_tracer
