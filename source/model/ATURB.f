#include "rundeck_opts.h"

      subroutine atm_diffus(lbase_min,lbase_max,dtime)
!@sum  atm_diffus updates u,v,t,q due to
!@+  turbulent transport throughout all GCM layers
!@+  using a second order closure (SOC)
!@+  turbulence model developed at GISS, 2000.
!@auth Ye Cheng/G. Hartke (modifications by G. Schmidt)
!@ver  1.0 (from diffB347D6M20)
!@cont atm_diffus,getdz,dout,de_solver_main,de_solver_edge,lgcm,kgcm,
!@+    find_pbl_top,ave_uv_to_agrid,ave_uv_to_bgrid,ave_s_to_bgrid,
!@+    tom_BUV,tom_B
!@var lbase_min/max levels through which to apply turbulence (dummy)
!@var dtime time step
!@var qmin minimum value of specific humidity
!@var itest longitude at which to call dout
!@var jtest latitude at which to call dout
!@var call_diag logical variable whether dout is called
!@var level the turb. model level, 25 is to solve e differential eqn
!@var non_local logical variable, true is to turn on non-local calc.

      USE CONSTANT, only : grav,deltx,lhe,sha,by3,teeny
      USE MODEL_COM, only :
     *      im,jm,lm,sig,sige,u_3d=>u,v_3d=>v,t_3d=>t,q_3d=>q,p,itime
cc      USE QUSDEF, only : nmom,zmoms,xymoms
cc      USE SOMTQ_COM, only : tmom,qmom
      USE GEOM, only : imaxj,kmaxj,ravj,idij,idjj,bydxyp,dxyp
      USE DYNAMICS, only : pk,pdsig,plij,pek,byam,am,dke
      USE DAGCOM, only : ajl,jl_trbhr,jl_damdc,jl_trbke,jl_trbdlht
#ifdef TRACERS_ON
      USE TRACER_COM, only : ntm,itime_tr0,trm  !,trmom
      USE TRACER_DIAG_COM, only: tajln,jlnt_turb
      USE FLUXES, only : trflux1
#endif
      USE SOCPBL, only : g0,g5,g6,g7,b1,b123,b2,prt,kappa,zgs
      USE PBLCOM, only : tsavg,qsavg,dclev,uflux,vflux,tflux,qflux
     *     ,e_3d=>egcm,w2_3d=>w2gcm,t2_3d=>t2gcm
      USE FLUXES, only : uflux1,vflux1,tflux1,qflux1


      IMPLICIT NONE

      integer, intent(in) :: lbase_min,lbase_max
      real*8, intent(in) :: dtime

      real*8, parameter :: qmin=1.d-20
      integer, parameter :: level=25,itest= 1,jtest=11
      logical, parameter :: non_local=.false.,call_diag=.false.
      integer, SAVE :: ifirst=0

      real*8, dimension(lm) :: u,v,t,q,e,t2,u0,v0,t0,q0,e0,t20
     &    ,tau,dudz,dvdz,dtdz,dqdz,g_alpha,as2,an2
     &    ,rhoebydz,bydzerho,rho,rhoe,dz,dze,gm,gh
     &    ,km,kh,kq,ke,kt2,kwt,gc_wt,gc_wq,gc_ew,gc_w2t,gc_wt2
     &    ,lscale,qturb,p2,p3,p4,rhobydze,bydzrhoe,uw,vw,w2,wt,wt0
     &    ,gc_wt_by_t2,wq,wq0

      real*8, dimension(lm,im,jm) :: u_3d_old,rho_3d,rhoe_3d,dz_3d
     &    ,dze_3d,u_3d_agrid,v_3d_agrid,t_3d_virtual,km_3d,km_3d_bgrid
     &    ,dz_3d_bgrid,dze_3d_bgrid,rho_3d_bgrid,rhoe_3d_bgrid
     &    ,wt_3d,v_3d_old
      real*8, dimension(im,jm) :: tvsurf,uflux_bgrid,vflux_bgrid
cc      real*8, dimension(nmom,lm) :: tmomij,qmomij

      real*8 :: uflx,vflx,tvflx,qflx,tvs
     &   ,ustar2,dbll,t0ijl,rak,alpha1,ustar
     &   ,flux_bot,flux_top,x_surf,dgcdz
      integer :: idik,idjk,
     &    i,j,l,k,n,iter !@i,j,l,k,n,iter loop variable
      real*8, save :: cc1,cc2
#ifdef TRACERS_ON
!@var tr0ij initial vertical tracer concentration profile (kg/kg)
!@var trij vertical tracer concentration profile (kg/kg)
!@var trmomij vertical tracer concentration moment profile (kg/kg)
      real*8, dimension(lm,ntm) :: tr0ij,trij
cc      real*8, dimension(nmom,lm,ntm) :: trmomij
!@var trflx surface tracer flux (-w tr) (kg/kg m/s)
      real*8, dimension(ntm) :: trflx
#endif

      ! Note that lbase_min/max are here for backwards compatibility
      ! with original drycnv. They are only used to determine where the
      ! routine has been called from.

      if (lbase_min.eq.2) return       ! quit if called from main

      if(ifirst.eq.0) then
          ifirst=1
          cc1=(1.-g7)/(5.*g5)
          cc2=cc1*0.5d0*(g6+g7)
      endif

      !  convert input T to virtual T
      do j=1,jm
        do i=1,imaxj(j)
          !@var tvsurf(i,j) surface virtual temperature
          !@var tsavg(i,j) COMPOSITE SURFACE AIR TEMPERATURE (K)
          tvsurf(i,j)=tsavg(i,j)*(1.d0+deltx*qsavg(i,j))
          do l=1,lm
            ! t_3d_virtual is virtual potential temp. referenced at 1 mb
            t_3d_virtual(l,i,j)=t_3d(i,j,l)*(1.d0+deltx*q_3d(i,j,l))
          end do
        end do
      end do

      ! integrate equations other than u,v at agrids

      ! get u_3d_agrid and v_3d_agrid
      call ave_uv_to_agrid(u_3d,v_3d,u_3d_agrid,v_3d_agrid,im,jm,lm)

      call getdz(t_3d_virtual,p,dz_3d,dze_3d,rho_3d,rhoe_3d,tvsurf
     &    ,im,jm,lm)

      loop_j_tq: do j=1,jm
        loop_i_tq: do i=1,imaxj(j)

          do l=1,lm
            u(l)=u_3d_agrid(l,i,j)
            v(l)=v_3d_agrid(l,i,j)
            ! virtual potential temp. referenced at 1 mb
            t(l)=t_3d_virtual(l,i,j)
            q(l)=q_3d(i,j,l)
cc            qmomij(:,l)=qmom(:,i,j,l)
cc            tmomij(:,l)=tmom(:,i,j,l) ! vert. grad. should virtual ?
            if(q(l).lt.qmin) q(l)=qmin
            e(l)=e_3d(l,i,j) !e_3d was called egcM
            t2(l)=t2_3d(l,i,j)
            rho(l)=rho_3d(l,i,j)
            rhoe(l)=rhoe_3d(l,i,j)
            t0(l)=t(l)
            q0(l)=q(l)
            e0(l)=e(l)
            t20(l)=t2(l)
            qturb(l)=sqrt(2.d0*e(l))
            dze(l)=dze_3d(l,i,j)
            bydzerho(l)=1.d0/(dze(l)*rho(l))
            rhobydze(l)=rho(l)/dze(l)
#ifdef TRACERS_ON
            do n=1,ntm
              if (itime_tr0(n).le.itime) then
                trij(l,n)=trm(i,j,l,n)*byam(l,i,j)*bydxyp(j)
cc                trmomij(:,l,n)=trmom(:,i,j,l,n)
                tr0ij(l,n)=trij(l,n)
              end if
            end do
#endif
          end do
          do l=1,lm-1
            dz(l)=dz_3d(l,i,j)
            bydzrhoe(l+1)=1.d0/(dz(l)*rhoe(l+1))
            rhoebydz(l+1)=rhoe(l+1)/dz(l)
          end do

          ! tvs is surface virtual potential temp. referenced at 1 mb
          tvs=tvsurf(i,j)/pek(1,i,j)
          uflx=uflux1(i,j)/rhoe(1)
          vflx=vflux1(i,j)/rhoe(1)
          qflx=qflux1(i,j)/rhoe(1)
          ! tvflx is virtual, potential temp. flux referenced at 1 mb
          tvflx=tflux1(i,j)*(1.d0+deltx*qsavg(i,j))/(rhoe(1)*pek(1,i,j))
     &         +deltx*tsavg(i,j)/pek(1,i,j)*qflx
          ! redefine uflux1,vflux1 for later use
          uflux1(i,j)=uflx
          vflux1(i,j)=vflx
#ifdef TRACERS_ON
          do n=1,ntm
            if (itime_tr0(n).le.itime) then
C**** minus sign needed for ATURB conventions
              trflx(n)=-trflux1(i,j,n)*bydxyp(j)/rhoe(1)
            end if
          end do
#endif

          if(abs(uflx).lt.teeny) uflx=sign(teeny,uflx)
          if(abs(vflx).lt.teeny) vflx=sign(teeny,vflx)
          ustar=(uflx*uflx+vflx*vflx)**(0.25d0)
          ustar2=ustar*ustar
          alpha1=atan2(vflx,uflx)

          ! calculate z-derivatives at the surface

          ! @var zgs height of surface layer (m), imported from SOCPBL
          dudz(1)=ustar/(kappa*zgs)*cos(alpha1)
          dvdz(1)=ustar/(kappa*zgs)*sin(alpha1)
          dtdz(1)=(tvflx*prt/ustar)/(kappa*zgs)
          dqdz(1)=(qflx*prt/ustar)/(kappa*zgs)

          g_alpha(1)=grav/tvs

          ! calculate z-derivatives on the edges of the layers

          do l=2,lm
            dudz(l)=(u(l)-u(l-1))/dz(l-1)
            dvdz(l)=(v(l)-v(l-1))/dz(l-1)
            dtdz(l)=(t(l)-t(l-1))/dz(l-1)
            dqdz(l)=(q(l)-q(l-1))/dz(l-1)
            g_alpha(l)=grav*2.d0/(t(l)+t(l-1))
          end do

          !@var an2 brunt-vassala frequency
          !@var as2 shear number squared

          do l=1,lm
            an2(l)=g_alpha(l)*dtdz(l)
            as2(l)=dudz(l)*dudz(l)+dvdz(l)*dvdz(l)
          end do

          ! calculate turbulence length scale lscale
          call lgcm(lscale,qturb,as2,an2,dze,rhoe,lm)

          ! calculate turbulent diffusivities km,kh,kq,ke and kt2
          call kgcm(km,kh,kq,ke,kt2,kwt,gc_wt,gc_wq,gc_ew,gc_w2t,gc_wt2
     1        ,gc_wt_by_t2,gm,gh
     2        ,uw,vw,w2,wt,wq,tau,u,v,t,e,qturb,t2,dudz,dvdz,as2,dtdz
     3        ,dqdz,g_alpha,an2,lscale,dz,dze,non_local,level,lm)

          ! integrate differential eqn for e
          p2(1)=0. ; p2(lm)=0.
          p3(1)=0. ; p3(lm)=0.
          p4(1)=0. ; p4(lm)=0.
          do l=2,lm-1
              dgcdz=(gc_ew(l)-gc_ew(l-1))/dz(l-1)
              p2(l)=0.d0
              p3(l)=2.d0*(qturb(l)/(b1*lscale(l)))
              p4(l)=-uw(l)*dudz(l)-vw(l)*dvdz(l)+g_alpha(l)*wt(l)-dgcdz
          end do
          x_surf=0.5d0*b123*ustar2
          call de_solver_edge(e,e0,ke,p2,p3,p4,
     &        rhobydze,bydzrhoe,x_surf,dtime,lm)

          do l=1,lm
              qturb(l)=sqrt(2.d0*e(l))
          end do

          if(level.eq.3) then

              ! integrate differential eqn for t2
              do l=2,lm-1
                  dgcdz=(gc_wt2(l)-gc_wt2(l-1))/dz(l-1)
                  p3(l)=2.d0*(qturb(l)/(b2*lscale(l))
     &                 +gc_wt_by_t2(l)*dtdz(l))
c                 p4(l)=-2.d0*wt(l)*dtdz(l)-dgcdz
                  p4(l)=2.d0*kh(l)*dtdz(l)*dtdz(l)-dgcdz
              end do
              x_surf=tvflx*tvflx*b2*prt/(ustar2*b1**by3)
              call de_solver_edge(t2,t20,kt2,p2,p3,p4,
     &        rhobydze,bydzrhoe,x_surf,dtime,lm)

c         elseif(level.eq.31) then

c             ! integrate differential eqn for wt
c             do l=2,lm-1
c                 dgcdz=(gc_w2t(l)-gc_w2t(l-1))/dz(l-1)
c                 p3(l)=g5/tau(l)-cc2*tau(l)*(dudz(l)+dvdz(l))
c                 p4(l)=-w2(l)*dtdz(l)+g0*g_alpha(l)*t2(l)
c    &            +cc1*tau(l)*(uw(l)+vw(l))*dtdz(l)-dgcdz
c             end do
c             x_surf=-tvflx
c             call de_solver_edge(wt,wt0,kwt,p2,p3,p4,
c    &            rhobydze,bydzrhoe,x_surf,dtime,lm)

          endif

          ! integrate differential eqn for T
          do l=2,lm-1
              p3(l)=0.d0
              p4(l)=-(gc_wt(l+1)-gc_wt(l))/dze(l)
          end do
          flux_bot=rhoe(1)*(tvflx+gc_wt(2)-gc_wt(1))
          flux_top=rhoe(lm)*gc_wt(lm)
          call de_solver_main(t,t0,kh,p2,p3,p4,
     &        rhoebydz,bydzerho,flux_bot,flux_top,dtime,lm)

C**** also diffuse moments
cc        call diff_mom(tmomij)

          ! integrate differential eqn for Q
          do l=2,lm-1
              p3(l)=0.d0
              p4(l)=-(gc_wq(l+1)-gc_wq(l))/dze(l)
          end do
          flux_bot=rhoe(1)*(qflx+gc_wq(2)-gc_wq(1))
          flux_top=rhoe(lm)*gc_wq(lm)
          call de_solver_main(q,q0,kq,p2,p3,p4,
     &        rhoebydz,bydzerho,flux_bot,flux_top,dtime,lm)
          do l=1,lm
              if(q(l).lt.qmin) q(l)=qmin
          end do

C**** also diffuse moments
cc        call diff_mom(qmomij)

#ifdef TRACERS_ON
C**** Use q diffusion coefficient for tracers
C**** Note that non-local terms for tracers are not properly defined
C**** and so will not work properly if used.
          do n=1,ntm
            if (itime_tr0(n).le.itime) then
              p2=0. ; p3=0. ; p4=0.
              flux_bot=rhoe(1)*trflx(n) !tr0ij(1,n)
              flux_top=0.     !rhoe(lm)*gc_wq(lm)
              call de_solver_main(trij(1,n),tr0ij(1,n),kq,p2,p3,p4,
     &             rhoebydz,bydzerho,flux_bot,flux_top,dtime,lm)
c     c        call diff_mom(trmomij)
            end if
          end do
#endif

          !@var dbll PBL top layer number counted from below, real*8
          call find_pbl_top(e,dbll,lm)
          dclev(i,j)=dbll

          do l=1,lm
            ! update 3-d t,q,e,t2 and km
            t0ijl=t_3d(i,j,l)
            t_3d(i,j,l)=t(l)/(1.d0+deltx*q(l))
C**** moment variation to be added
cc            qmom(:,i,j,l)=qmomij(:,l)
cc            tmom(:,i,j,l)=tmomij(:,l)

            q_3d(i,j,l)=q(l)
            e_3d(l,i,j)=e(l)
            w2_3d(l,i,j)=w2(l)
            t2_3d(l,i,j)=t2(l)
            km_3d(l,i,j)=km(l)
            ! ACCUMULATE DIAGNOSTICS for t and q
            AJL(J,L,JL_TRBHR)=AJL(J,L,JL_TRBHR)
     2                 +(T_3d(I,J,L)-t0ijl)*PK(L,I,J)*PLIJ(L,I,J)
            AJL(J,L,JL_TRBDLHT)=AJL(J,L,JL_TRBDLHT)
     2                 +(Q_3d(I,J,L)-q0(l))*PDSIG(L,I,J)*LHE/SHA
            AJL(J,L,JL_TRBKE)=AJL(J,L,JL_TRBKE)+e(l)
#ifdef TRACERS_ON
            do n=1,ntm
              if (itime_tr0(n).le.itime) then
                tajln(j,l,jlnt_turb,n)=tajln(j,l,jlnt_turb,n) +
     &               (trij(l,n)-tr0ij(l,n))*am(l,i,j)*dxyp(j)
                trm(i,j,l,n)=trij(l,n)*am(l,i,j)*dxyp(j)
cc                trmom(:,i,j,l,n)=trmomij(:,l,n)
              end if
            end do
#endif

          end do

          ! Write out diagnostics if at selected grid point:

          if (call_diag.and.(i.eq.itest).and.(j.eq.jtest)) then
            call dout(u,v,t,q,e,t2,dz,dze,dudz,dvdz,as2,dtdz,g_alpha,
     1      an2,dqdz,rho,rhoe,km,kh,kq,ke,kt2,gc_wt,gc_wq,gc_ew,
     2      gc_w2t,gc_wt2,gm,gh,lscale,tvs,uflx,vflx,tvflx,qflx,
     3      i,j,iter,non_local,level,lm)
          endif

        end do loop_i_tq
      end do loop_j_tq

      ! integrate U,V equations at bgrids

      call ave_uv_to_bgrid(uflux1,vflux1,uflux_bgrid,vflux_bgrid,
     &                     im,jm,1)
      call ave_s_to_bgrid(km_3d,km_3d_bgrid,im,jm,lm)
      call ave_s_to_bgrid(dz_3d,dz_3d_bgrid,im,jm,lm)
      call ave_s_to_bgrid(dze_3d,dze_3d_bgrid,im,jm,lm)
      call ave_s_to_bgrid(rho_3d,rho_3d_bgrid,im,jm,lm)
      call ave_s_to_bgrid(rhoe_3d,rhoe_3d_bgrid,im,jm,lm)

      loop_j_uv: do j=2,jm
        loop_i_uv: do i=1,im

          do l=1,lm
            u(l)=u_3d(i,j,l)
            v(l)=v_3d(i,j,l)
            rho(l)=rho_3d_bgrid(l,i,j)
            rhoe(l)=rhoe_3d_bgrid(l,i,j)
            u0(l)=u(l)
            v0(l)=v(l)
            u_3d_old(l,i,j)=u(l)
            v_3d_old(l,i,j)=v(l)
            km(l)=km_3d_bgrid(l,i,j)
            dze(l)=dze_3d_bgrid(l,i,j)
            bydzerho(l)=1.d0/(dze(l)*rho(l))
          end do
          do l=1,lm-1
            dz(l)=dz_3d_bgrid(l,i,j)
            rhoebydz(l+1)=rhoe(l+1)/dz(l)
          end do

          uflx=uflux_bgrid(i,j)
          vflx=vflux_bgrid(i,j)

          ! integrate differential eqns for U and V
          do l=2,lm-1
              p3(l)=0.d0
              p4(l)=0.d0
          end do
          flux_bot=rhoe(1)*uflx
          flux_top=0.d0
          call de_solver_main(u,u0,km,p2,p3,p4,
     &        rhoebydz,bydzerho,flux_bot,flux_top,dtime,lm)
          flux_bot=rhoe(1)*vflx
          call de_solver_main(v,v0,km,p2,p3,p4,
     &        rhoebydz,bydzerho,flux_bot,flux_top,dtime,lm)

          ! update u and v
          do l=1,lm
            u_3d(i,j,l)= u(l)
            v_3d(i,j,l)= v(l)
          end do

        end do loop_i_uv
      end do loop_j_uv

      ! ACCUMULATE DIAGNOSTICS for u and v
      DO J=1,JM
        DO I=1,IMAXJ(J)
          DO L=1,LM
            DO K=1,KMAXJ(J)
              RAK=RAVJ(K,J)
              IDIK=IDIJ(K,I,J)
              IDJK=IDJJ(K,J)
              AJL(IDJK,L,JL_DAMDC)=AJL(IDJK,L,JL_DAMDC)
     &        +(U_3d(IDIK,IDJK,L)-u_3d_old(L,IDIK,IDJK))*PLIJ(L,I,J)*RAK
            ENDDO
          ENDDO
        ENDDO
      ENDDO

C**** Save additional changes in KE for addition as heat later
!$OMP  PARALLEL DO PRIVATE (L,I,J)
      DO L=1,LM
      DO J=2,JM
      DO I=1,IM
        DKE(I,J,L)=DKE(I,J,L)+0.5*(U_3d(I,J,L)*U_3d(I,J,L)
     *       +V_3d(I,J,L)*V_3d(I,J,L)-U_3d_old(L,I,J)*U_3d_old(L,I,J)
     *       -V_3d_old(L,I,J)*V_3d_old(L,I,J))
      END DO
      END DO
      END DO
!$OMP  END PARALLEL DO

      return
      end subroutine atm_diffus

      subroutine getdz(tv,p,dz,dze,rho,rhoe,tvsurf,im,jm,lm)
!@sum  getdz computes the 3-d finite difference dz and dze
!@+    as well as the 3-d density rho and rhoe
!@+    called at the primary grid (A-grid)
!@auth Ye Cheng/G. Hartke
!@ver  1.0
!@var  dz(l,i,j) z(l+1,i,j) - z(l,i,j)
!@var  dze(l,i,j) ze(l+1,i,j) - ze(l,i,j)
!@var  rho,rhoe air density at the main/edge grids
!@var  tvsurf surface virtual temperature
!@var  im,jm,lm 3-d grids
!@var  z vertical coordinate associated with SIG(l)
!@var  ze vertical coordinate associated with SIGE(l)
!@var  temp0 virtual temperature (K) at (i,j) and SIG(l)
!@var  temp1 virtual temperature (K) at (i,j) and SIG(l+1)
!@var  temp1e average of temp0 and temp1

c
c     Grids:
c
c                   -------------------------  lm+1
c                lm  - - - - - - - - - - - - -
c                   -------------------------  lm
c                l+1 - - - - - - - - - - - - -
c                    -------------------------  l+1
C     (main)     l   - - - - - - - - - - - - -
c                    -------------------------  l     (edge)
c                l-1 - - - - - - - - - - - - -
c                    -------------------------  l-1
c                2   - - - - - - - - - - - - -
c                    -------------------------    2
c                1   - - - - - - - - - - - - -
c                    -------------------------    1
c           rhoe(l+1,i,j)=100.d0*(pl-pl1)/(grav*dz(l,i,j))
c           rho(l,i,j)=100.d0*(ple-pl1e)/(grav*dze(l,i,j))
c
c     at main: u,v,tv,q,ke
c     at edge: e,lscale,km,kh,gm,gh
c
      USE CONSTANT, only : grav,rgas
      USE GEOM, only : imaxj
      USE DYNAMICS, only : pmid,pk,pedn

      implicit none

      integer, intent(in) :: im,jm,lm
      real*8, dimension(im,jm), intent(in) :: p,tvsurf
      real*8, dimension(lm,im,jm), intent(in) :: tv
      real*8, dimension(lm,im,jm), intent(out) :: rho,rhoe,dz,dze

      real*8 :: temp0,temp1,temp1e,pl1,pl,pl1e,ple,plm1e
      integer :: i,j,l  !@var i,j,l loop variable

      do j=1,jm
        do i=1,imaxj(j)
          do l=1,lm-1

            pl1 =pmid(l+1,i,j)
            pl  =pmid(l,i,j)
            pl1e=pedn(l+1,i,j)
            ple =pedn(l,i,j)
            temp0 =tv(l,i,j)*pk(l,i,j)
            temp1 =tv(l+1,i,j)*pk(l+1,i,j)
            temp1e=0.5d0*(temp0+temp1)
            dz(l,i,j)    =-(rgas/grav)*temp1e*log(pl1/pl)
            dze(l,i,j)=-(rgas/grav)*temp0 *log(pl1e/ple)
            rhoe(l+1,i,j)=100.d0*(pl-pl1)/(grav*dz(l,i,j))
            rho(l,i,j)=100.d0*(ple-pl1e)/(grav*dze(l,i,j))
            if(l.eq.1) then
              rhoe(1,i,j)=100.d0*ple/(tvsurf(i,j)*rgas)
c             rhoe(1,i,j)=2.d0*rho(1,i,j)-rhoe(2,i,j)
            endif
            if(l.eq.lm-1) then
c             rho(lm,i,j)=100.d0*pl1/(temp1*rgas)
              plm1e=pedn(lm+1,i,j)
              dze(lm,i,j)=-(rgas/grav)*temp1 *log(plm1e/pl1e)
              rho(lm,i,j)=100.d0*(pl1e-plm1e)/(grav*dze(lm,i,j))
            endif

          end do
        end do
      end do

      return
      end subroutine getdz

      subroutine dout(u,v,t,q,e,t2,dz,dze,dudz,dvdz,as2,dtdz,g_alpha,
     1      an2,dqdz,rho,rhoe,km,kh,kq,ke,kt2,gc_wt,gc_wq,gc_ew,
     2      gc_w2t,gc_wt2,gm,gh,lscale,tvs,uflx,vflx,tvflx,qflx,
     3      i,j,iter,non_local,level,n)
!@sum dout writes out diagnostics at i,j
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
!@var u  west-east   velocity component
!@var v  south-north velocity component
!@var t  virt. pot. temperature (referenced to 1mb)
!@var p  pressure at main grid z
!@var pe  pressure at secondary grid ze
!@var q  specific humidity
!@var e  turbulent kinetic energy
!@var t2  turbulent virt. pot. temperature variance
!@var dz(j) z(j+1)-z(j)
!@var dze(j) ze(j+1)-ze(j)
!@var dudz z-derivative of u at secondary grid ze
!@var dvdz z-derivative of v at secondary grid ze
!@var as2 dudz^2+dvdz^2
!@var dtdz z-derivative of t at secondary grid ze
!@var g_alpha grav*(thermal expansion coefficient)
!@var an2 g_alpha*dtdz, brunt-vasala frequency
!@var dqdz z-derivative of q at secondary grid ze
!@var rho  density at z
!@var rhoe  density at ze
!@var km turbulent viscosity for u and v equations
!@var kh turbulent diffusivity for t
!@var kq turbulent diffusivity for q
!@var ke turbulent diffusivity for e
!@var kt2 turbulent diffusivity for t2
!@var gc countergradient term is vertical heat flux (at model level 3)
!@var gm normalized velocity gradient, tau**2*as2
!@var gh normalized temperature gradient, tau**2*an2
!@var lscale  turbulent length scale
!@var tvs surface virtual temperature
!@var uflx momentun flux -uw at surface, ze(1)
!@var vflx momentun flux -vw at surface, ze(1)
!@var tvflx heat flux -wt at surface, ze(1)
!@var qflx moisture flux -wq at surface, ze(1)
!@var i/j location at which the output is wriiten
!@var n number of vertical main layers

      USE CONSTANT, only : grav,rgas
      USE MODEL_COM, only : sig,sige,ptop
      USE DYNAMICS, only : pmid,pedn,plij,pk,pek
      USE SOCPBL, only : b1,b2

      implicit none

      integer, intent(in) :: level,n,iter,i,j
      real*8, dimension(n), intent(in) ::
     &   u,v,t,q,e,t2,dz,dze,dudz,dvdz,as2,dtdz,g_alpha,
     &   an2,dqdz,rho,rhoe,km,kh,kq,ke,kt2,gc_wt,gc_wq,
     &   gc_ew,gc_w2t,gc_wt2,gm,gh,lscale
      real*8, intent(in) :: tvs,uflx,vflx,tvflx,qflx
      logical, intent(in) :: non_local

      real*8, dimension(n) :: p,pe,t2_out
      real*8 :: z,ze,dmdz,uf,hf,qf,ri,rif,sigmat,wt1
      integer :: l  !@var l loop variable

      do l=1,n
          p(l)=100.d0*pmid(l,i,j)
          pe(l)=100.d0*pedn(l,i,j)
          if(level.eq.25) then
            t2_out(l)=b2*lscale(l)/sqrt(2.d0*e(l))*kh(l)*dtdz(l)**2
          else
            t2_out(l)=t2(l)
          endif
      end do
      Write (67,1000) "iter=",iter, "i=",i,"j=",j
      Write (67,1100) "pe(1)=",pe(1)
      Write (67,1100) "uflx=",uflx,"vflx=",vflx
      Write (67,1100) "tvflx=",tvflx*pek(1,i,j),"qflx=",qflx

      ! Fields on main vertical grid:
      z=(rgas/grav)*0.5d0*(tvs*pek(1,i,j)+t(1)*pk(1,i,j))
     2  *log(pe(1)/p(1))+10.d0
      write (67,1500)
      do l=1,n-1
        write (67,2000) l,z,p(l),dz(l),u(l),v(l),t(l)*pk(l,i,j),
     2                  q(l),ke(l),gc_ew(l),gc_w2t(l),gc_wt2(l)
        z=z+dz(l)
      end do
      write (67,2500) l,z,p(n),u(n),v(n),t(n)*pk(n,i,j),q(n),
     2                0.,0.,0.,0.
      write (67,*)

      ! Fields on secondary vertical grid:

      write (67,3000)
      l=1
      ze=10.d0
      write (67,2100) l,ze,pe(l),lscale(l),e(l),t2_out(l)*pek(l,i,j)**2
      ze=10.d0+dze(1)
      do l=2,n
        wt1=-kh(l)*dtdz(l)*pek(l,i,j)
        write (67,2000) l,ze,pe(l),gc_wt(l),wt1,kh(l),kt2(l),
     2       gm(l),gh(l),lscale(l),e(l),t2_out(l)*pek(l,i,j)**2
        ze=ze+dze(l)
      end do
      write (67,*)

      ! Fluxes on the secondary grid:

      ze=10.d0
      write (67,4000)
      do l=2,n
        ze=ze+dze(l-1)
        dmdz=sqrt(as2(l))
        uf=km(l)*dmdz
        hf=kh(l)*dtdz(l)*pek(l,i,j)
        qf=kq(l)*dqdz(l)
        ri=an2(l)/as2(l)
        sigmat=km(l)/kh(l)
        rif=ri/sigmat
        write (67,2000) l,ze,uf,hf,qf,dmdz,dtdz(l)*pek(l,i,j),
     2                  dqdz(l),ri,rif,sigmat,gc_wq(l)
      end do

      write (67,*) "------------"
      write (67,*)

      return

1000  format(3(4x,a,i4))
1001  format(9(1pe14.4))
1100  format(2(4x,a,1pe14.4))
1500  format (1h ,' l',1x,
     1            '     z     ',1x,
     2            '     p     ',1x,'     dz    ',1x,'     u     ',1x,
     3            '     v     ',1x,'     t     ',1x,'     q     ',1x,
     4            '    ke     ',1x,'    gc_ew ',1x,'   gc_w2t ',1x,
     5            '   gc_wt2  ',/)
2000  format (1h ,i2,1x,1pe11.4,9(1x,1pe11.4),1x,1pe10.3)
2100  format (1h ,i2,2(1x,1pe11.4),24x,1x,11x,36x,2(1x,1pe11.4),
     2           ,1x,1pe10.3)
2500  format (1h ,i2,2(1x,1pe11.4),12x,7(1x,1pe11.4),1x,1pe10.3)
3000  format (1h ,' l',1x,
     2            '  z (edge) ',1x,
     2            '  p (edge) ',1x,'   gc_wt   ',1x,' -kh*dtdz  ',1x,
     3            '     kh    ',1x,'    kt2    ',1x,'     gm    ',1x,
     4            '     gh    ',1x,'   lscale  ',1x,'     e     ',1x,
     5            '     t2    ',/)
4000  format (1h ,' l',1x,
     2            '  z (edge) ',1x,
     3            '    uflux  ',1x,'    hflux  ',1x,'    qflux  ',1x,
     4            '     dmdz  ',1x,'    dtdz   ',1x,'    dqdz   ',1x,
     5            '     Ri    ',1x,'     Rif   ',1x,'   Sigmat  ',1x,
     6            '    gc_wq  ',/)

      end subroutine dout

      subroutine de_solver_main(x,x0,p1,p2,p3,p4,
     &    rhoebydz,bydzerho,flux_bot,flux_top,dtime,n)
!@sum differential eqn solver for x using tridiagonal method
!@+   d/dt x = d/dz (P1 d/dz x) - P2 d/dz x - P3 x + P4
!@+   where p2=0 and x is at the main grid
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
!@var x the unknown to be solved (at main drid)
!@var x0 x at previous time step
!@var p1,p2,p3,p4 coeff. of the d.e.
!@var rhoebydz(j) rhoe(j+1)/dz(j)
!@var bydzerho(j) 1/(dze(j)*rho(j))
!@var flux_bot==rhoe(1)*(flux+gc(2)-gc(1))
!@var flux_top=rhoe(n)*gc(n)
!@var dtime time step
!@var n number of vertical edge grid

      implicit none

      integer, intent(in) :: n
      real*8, dimension(n), intent(in) ::
     &        x0,p1,p2,p3,p4,rhoebydz,bydzerho
      real*8, intent(in) :: flux_bot,flux_top,dtime
      real*8, dimension(n), intent(out) :: x

      real*8, dimension(n) :: sub,dia,sup,rhs
      real*8 :: alpha
      integer :: j  !@var j loop variable

c     sub(j)*x_jm1_kp1+dia(j)*x_j_kp1+sup(j)*x_jp1_kp1 = rhs(j)
c     k refers to time step, j refers to main grid
c     except p1(j) which is defined on the edge grid

      do j=2,n-1
          sub(j)=-dtime*p1(j)*rhoebydz(j)*bydzerho(j)
          sup(j)=-dtime*p1(j+1)*rhoebydz(j+1)*bydzerho(j)
          dia(j)=1.d0-(sub(j)+sup(j))+dtime*p3(j)
          rhs(j)=x0(j)+dtime*p4(j)
      end do

c     Lower boundary conditions(x=T) :
c     d/dt T = -d/dz wt where
c     d/dt T = (T(1)-T0(1))/dtime
c     -d/dz wt = -(wt(2)-wt(1))/dze(1), dze(1)=ze(2)-ze(1)
c     wt(2)=-p1(2)*(T(2)-T(1))/dz(1)+gc(2), dz(1)=z(2)-z(1)
c     wt(1)=-tvflx+gc(1)
c     flux_bot==rhoe(1)*(tvflx+gc(2)-gc(1))
c     rhoebydz and bydzerho are in place to balance mass flux

      alpha=dtime*p1(2)*rhoebydz(2)*bydzerho(1)
      dia(1)=1.d0+alpha
      sup(1)=-alpha
      rhs(1)=x0(1)-dtime*bydzerho(1)*flux_bot

c     Upper boundary conditions:

c     d/dt T = -d/dz wt where
c     d/dt T = (T(n)-T0(n))/dtime
c     -d/dz wt = -(wt(n+1)-wt(n))/dze(n), dze(n)=ze(n+1)-ze(n)
c     wt(n)=-p1(n)*(T(n)-T(n-1))/dz(n-1)+gc(n), dz(n-1)=z(n)-z(n-1)
c     wt(n+1)=0,gc(n+1)=0.
c     flux_top=rhoe(n)*flux_top

      alpha=dtime*p1(n)*rhoebydz(n)*bydzerho(n)
      sub(n)=-alpha
      dia(n)=1.d0+alpha
      rhs(n)=x0(n)+dtime*bydzerho(n)*flux_top

      call tridiag(sub,dia,sup,rhs,x,n)

      return
      end subroutine de_solver_main

      subroutine de_solver_edge(x,x0,p1,p2,p3,p4,
     &    rhobydze,bydzrhoe,x_surf,dtime,n)
!@sum differential eqn solver for x using tridiagonal method
!@+   d/dt x = d/dz (P1 d/dz x) - P2 d/dz x - P3 x + P4
!@+   where p2=0 and x is at the edge grid
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
!@var x the unknown to be solved (at edge drid)
!@var x0 x at previous time step
!@var p1,p2,p3,p4 coeff. of the d.e.
!@var rhobydze(j) rho(j)/dze(j)
!@var bydzrhoe(j) 1/(dz(j-1)*rhoe(j))
!@var x_surf surface value of x
!@var dtime time step
!@var n number of vertical edge grid

      use CONSTANT, only : teeny
      implicit none

      integer, intent(in) :: n
      real*8, dimension(n), intent(in) ::
     &        x0,p1,p2,p3,p4,rhobydze,bydzrhoe
      real*8, intent(in) :: x_surf,dtime
      real*8, dimension(n), intent(out) :: x

      real*8, dimension(n) :: sub,dia,sup,rhs
      real*8 :: alpha
      integer :: j  !@var j loop variable

c     sub(j)*x_jm1_kp1+dia(j)*x_j_kp1+sup(j)*x_jp1_kp1 = rhs(j)
c     k refers to time step, j refers to edge grid
c     except p1(j) which is defined on the main grid

      do j=2,n-1
          sub(j)=-dtime*p1(j-1)*rhobydze(j-1)*bydzrhoe(j)
          sup(j)=-dtime*p1(j)*rhobydze(j)*bydzrhoe(j)
          dia(j)=1.d0-(sub(j)+sup(j))+dtime*p3(j)
          rhs(j)=x0(j)+dtime*p4(j)
      end do

c     Boundary conditions:

      dia(1)=1.d0
      sup(1)=0.d0
      rhs(1)=x_surf

      sub(n)=-1.d0
      dia(n)=1.d0
      rhs(n)=0.d0

      call tridiag(sub,dia,sup,rhs,x,n)
c
      do j=1,n
         if(x(j).lt.teeny) x(j)=teeny
      end do

      return
      end subroutine de_solver_edge

      subroutine lgcm(lscale,qturb,as2,an2,dze,rhoe,n)
!@sum lgcm calculates the z-profle of the turbulent length scale
!@auth Ye Cheng/G. Hartke
!@ver  1.0
!@var lscale z-profle of turbulent length scale
!@var qturb z-profle of sqrt(2*e)
!@var as2 z-profle of dudz^2+dvdz^2
!@var an2 z-profle of g_alpha*dtdz
!@var dze(j) ze(j+1)-ze(j)
!@var rhoe the z-profile of the density at the layer edge
!@var n number of GCM layers
!@var zgs height of surface layer (m), imported from SOCPBL

      USE CONSTANT, only : grav,teeny
      USE SOCPBL, only : kappa,zgs

      implicit none

      integer, intent(in) :: n
      real*8, dimension(n), intent(in) :: qturb,as2,an2,dze,rhoe
      real*8, dimension(n), intent(out) :: lscale

      real*8, dimension(n) :: ze
      real*8, parameter :: alpha0=0.1d0
      real*8 :: sum1,sum2,qj,qjm1,l0,l1,kz,lmax,lmax2
      integer :: j  !@var j loop variable

      ze(1)=zgs
      do j=2,n
          ze(j)=ze(j-1)+dze(j-1)
      end do

c     integration of monotonically tabulated function by
c     trapezoidal rule

      sum1=0.d0
      sum2=0.d0
      do j=2,n
        qj=qturb(j)
        qjm1=qturb(j-1)
        sum1=sum1+dze(j-1)*(qj*rhoe(j)+qjm1*rhoe(j-1))
        sum2=sum2+dze(j-1)*(qj*rhoe(j)*ze(j)+
     &                         qjm1*rhoe(j-1)*ze(j-1))
      end do
      l0=alpha0*sum2/sum1

      kz=kappa*ze(1)
      lscale(1)=l0*kz/(l0+kz)
      if (lscale(1).gt.dze(1)) lscale(1)=dze(1)
      do j=2,n
        kz=kappa*ze(j)
        l1=l0*kz/(l0+kz)
        if (an2(j).gt.0.d0) then  !must NOT use this for level3
          lmax  =0.53d0*qturb(j)/(sqrt(an2(j))+teeny)
          lmax2 =1.95d0*qturb(j)/(sqrt(as2(j))+teeny)
          lmax=min(lmax,lmax2)
          if (l1.gt.lmax) l1=lmax
        endif
        if (l1.gt.dze(j)) l1=dze(j)
        lscale(j)=l1
      end do

      return
      end subroutine lgcm

      subroutine kgcm(km,kh,kq,ke,kt2,kwt,gc_wt,gc_wq,gc_ew,gc_w2t
     2    ,gc_wt2,gc_wt_by_t2
     3    ,gm,gh,uw,vw,w2,wt,wq,tau,u,v,t,e,qturb,t2,dudz,dvdz,as2
     4    ,dtdz,dqdz,g_alpha,an2,lscale,dz,dze,non_local,level,n)
c     dz(j)==z(j+1)-z(j), dze(j)==ze(j+1)-ze(j)
c     at main: u,v,t,q,ke
c     at edge: e,lscale,km,kh,gm,gh
!@sum kgcm computes the turbulent stability functions Km, Kh and Ke
!@+   using the GISS second order closure model (2000)
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
!@var u,v,t,e,lscale z-profiles
!@var dz(j) z(j+1)-z(j)
!@var km turbulent viscosity for u and v equations
!@var kh turbulent conductivity for t and q equations
!@var ke turbulent diffusivity for e equation
!@var gc countergradient term in heat flux (level 3)
!@var gm normalized velocity gradient, tau**2*as2
!@var gh normalized temperature gradient, tau**2*an2
!@var n number of GCM layers
!@var tau B1*lscale/sqrt(2*e)
!@var as2 shear squared, (dudz)**2+(dvdz)**2
!@var an2 Brunt-Vaisala frequency, grav/T*dTdz
!@var se stability constant for e, adjustable
!@var st2 stability constant for t2, adjustable
      USE CONSTANT, only : grav,teeny,by3
      USE SOCPBL, only :prt,ghmin,ghmax,gmmax0,d1,d2,d3,d4,d5
     *     ,s0,s1,s2,s4,s5,s6,b1,b2
     *     ,d1_3,d2_3,d3_3,d4_3,d5_3
     *     ,s0_3,s1_3,s2_3,s3_3,s4_3,s5_3,s6_3
     *     ,g0,g1,g2,g3,g4,g5,g6,g7,g8

      implicit none

      integer, intent(in) :: level,n    !@var n  array dimension
      real*8, dimension(n), intent(in) :: u,v,t,e,qturb,dudz,dvdz,as2
     &    ,dtdz,dqdz,g_alpha,an2,lscale,dz,dze
      real*8, dimension(n), intent(inout) :: t2,kwt
      real*8, dimension(n), intent(out) ::
     &     km,kh,kq,ke,kt2,gc_wt,gc_wq,gc_ew,gc_w2t,gc_wt2
     &    ,gc_wt_by_t2,gm,gh,uw,vw,w2,wt,wq,tau
      logical, intent(in) :: non_local

      real*8, parameter :: se=0.1d0,st2=0.06d0,kmax=600.d0
     &  ,kmmin=1.5d-5,khmin=2.5d-5,kqmin=2.5d-5
     &  ,kemin=1.5d-5,kt2min=1.5d-5

      real*8, dimension(n) :: u2,v2,uv,ut,vt
      real*8 :: ell,byden,ghj,gmj,gmmax
     &    ,sm,sh,sq,taue,e_lpbl
     &    ,kh_canuto,sig,sw,tpj,tpjm1,tppj,w3pj,taupj,m
     &    ,g_alphaj,tauj,dudzj,dvdzj,as2j,an2j,dtdzj
     &    ,uwj,vwj,w2j,wtj,du2dz,dv2dz,dw2dz
     &    ,duvdz,duwdz,dvwdz,dutdz,dvtdz,dwtdz,dt2dz,ke0
     &    ,x,dedz,aa,bb,cc,ghmin3,bydz,tmp,sq_by_sh
      integer :: j  !@var j loop variable
      integer, save :: ifirst=0
      real*8, save :: c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13
     &    ,c14,c15,c16,c17,g9,twoby3,byg5,b2byb1

      if(ifirst.eq.0) then
          ifirst=1
          twoby3=2.*by3
          byg5=1./g5
          b2byb1=b2/b1
          c0=0.5d0*g1
          c1=g4*(g3+by3*g2+0.5d0*(g6+g7)/g5)
          c2=g4/g5
          c3=(g3**2-by3*g2**2)
          c4=(g2+3*g3)*by3
          c5=2*g2*by3
          c6=2*g4*by3
          c7=(3*g3-g2)*by3
          c8=4*g4*by3
          c9=(g2+g3)*0.5d0
          c10=(g6+g7)*0.5d0
          ! c11-c14 are for level 3
          c11=g0*byg5
          c12=2*d3_3+s5_3
          c13=2*d1_3+s4_3
          c14=c2*c11
c         g9=0.5d0*g8   ! then sq=sh
c         c15=(g8-g9)/g5
c         c16=g9/g5
c         c17=(g6*g6-g7*g7)/(4*g5*g5)

      endif

      do j=1,n
        ell=lscale(j)
        tauj=B1*ell/(qturb(j)+teeny)
        tau(j)=tauj
        ghj=tauj*tauj*an2(j)
        gmj=tauj*tauj*as2(j)
        if(level.eq.25) then
          if(ghj.lt.ghmin) ghj=ghmin
          if(ghj.gt.ghmax) ghj=ghmax
          gmmax=(1+d1*ghj+d3*ghj*ghj)/(d2+d4*ghj)
          gmmax=min(gmmax,gmmax0)
          if(gmj.gt.gmmax) gmj=gmmax
          byden=1./(1+d1*ghj+d2*gmj+d3*ghj*ghj+d4*ghj*gmj+d5*gmj*gmj)
          sm=(s0+s1*ghj+s2*gmj)*byden
          sh=(s4+s5*ghj+s6*gmj)*byden
c         sq_by_sh =(1+c15*ghj-c17*gmj)/(1+c16*ghj-c17*gmj)
c         sq=sq_by_sh*sh
          sq=sh
          gc_wt(j)=0.d0
          gc_wq(j)=0.d0
          gc_ew(j)=0.d0
          gc_w2t(j)=0.d0
          gc_wt2(j)=0.d0
        elseif(level.eq.3) then
          tmp=g_alpha(j)*tauj
          tmp=tmp*tmp*t2(j)/(e(j)+teeny)
          aa=c12
          bb=c13-c14*tmp
          cc=2-c11*tmp
          ghmin3=(-bb+sqrt(bb*bb-4*aa*cc))/(2*aa)
c         ghmin3=int(ghmin3*10000.)/10000.
          if(ghj.lt.ghmin3) ghj=ghmin3
          if(ghj.gt.ghmax) ghj=ghmax
          gmmax=(1+d1_3*ghj+d3_3*ghj*ghj)/(d2_3+d4_3*ghj)
          if(gmmax.lt.teeny) gmmax=teeny
          gmmax=min(gmmax,gmmax0)
          if(gmj.gt.gmmax) gmj=gmmax
          byden=1./(1+d1_3*ghj+d2_3*gmj+d3_3*ghj*ghj+d4_3*ghj*gmj
     &         +d5_3*gmj*gmj)
          sm=(s0_3+s1_3*ghj+s2_3*gmj+s3_3*tmp)*byden
          sh=(s4_3+s5_3*ghj+s6_3*gmj)*byden
          gc_wt_by_t2(j)=c11*(1+c2*ghj+c3*gmj)*byden*g_alpha(j)*tauj
          gc_wt(j)=gc_wt_by_t2(j)*t2(j) ! used in d.e. for T
c         gc_wq(j)=gc_wt(j)*dqdz(j)/(dtdz(j)+teeny) ! used in d.e. for Q
c         gc_wq(j)=-gc_wt_by_t2(j)*prt*tauj*wt(j)*dqdz(j)
          gc_wq(j)=0.d0
          gc_ew(j)=0.d0               ! used in d.e. for e
          gc_wt2(j)=0.d0              ! used in d.e. for t2
          gc_w2t(j)=0.d0              ! used in d.e. for wt

c         when e and wt are solved prognostically
c         ref: /u/acyxc/papers/2ndOrder/2001/uiujhi_2001_wt
c         Sm=(c0+c1*tau(j)*g_alpha(j)*wt(j)/e(j))/(1+c2*ghj+c3*gmj)
        endif
        taue=tauj*e(j)
        km(j)=min(max(taue*sm,kmmin),kmax)
        kh(j)=min(max(taue*sh,khmin),kmax)
c       kq(j)=min(max(taue*sq,kqmin),kmax)
        if(level.eq.25) then
          kq(j)=kh(j)
        elseif(level.eq.3) then
          byden=1./(1+d1*ghj+d2*gmj+d3*ghj*ghj+d4*ghj*gmj+d5*gmj*gmj)
          sh=(s4+s5*ghj+s6*gmj)*byden
          kq(j)=min(max(taue*sh,kqmin),kmax)
        endif
c       ke(j)=min(max(taue*se,kemin),kmax)
        ke(j)=5.*km(j)
        kt2(j)=min(max(taue*st2,kt2min),kmax)
        gm(j)=gmj
        gh(j)=ghj
        uw(j)=-km(j)*dudz(j)
        vw(j)=-km(j)*dvdz(j)
        wt(j)=-kh(j)*dtdz(j)+gc_wt(j)
        wq(j)=-kq(j)*dqdz(j)+gc_wq(j)
        w2j=twoby3*e(j)+tauj*(c7*(
     2    dudz(j)*uw(j)+dvdz(j)*vw(j))+c8*g_alpha(j)*wt(j))
        w2(j)=max(0.24d0*e(j),w2j) ! Mellor-Yamada 1982
      enddo
      if(level.eq.25) then
        do j=1,n
          t2(j)=b2byb1*tau(j)*kh(j)*dtdz(j)*dtdz(j)
        end do
      endif

        do j=1,n-1
          ke(j)=0.5d0*(ke(j)+ke(j+1))    !defined on main levels
          kt2(j)=0.5d0*(kt2(j)+kt2(j+1)) !defined on main levels
        end do

      if(non_local) then

        do j=1,n
          tauj=tau(j)
          u2(j)=twoby3*e(j)-tauj*(c4*dudz(j)*uw(j)
     2        -c5*dvdz(j)*vw(j)+c6*g_alpha(j)*wt(j))
          v2(j)=twoby3*e(j)-tauj*(c4*dvdz(j)*vw(j)
     2        -c5*dudz(j)*uw(j)+c6*g_alpha(j)*wt(j))
          w2(j)=twoby3*e(j)+tauj*(c7*(
     2    dudz(j)*uw(j)+dvdz(j)*vw(j))+c8*g_alpha(j)*wt(j))
          uv(j)=-c9*tauj*(dvdz(j)*uw(j)+dudz(j)*vw(j))
          ut(j)=-byg5*tauj*(dtdz(j)*uw(j)+c10*dudz(j)*wt(j))
          vt(j)=-byg5*tauj*(dtdz(j)*vw(j)+c10*dvdz(j)*wt(j))
        end do

        do j=1,n-1  ! on main grids
          g_alphaj=0.5d0*(g_alpha(j)+g_alpha(j+1))
          tauj=0.5d0*(tau(j)+tau(j+1))
          dudzj=0.5d0*(dudz(j)+dudz(j+1))
          dvdzj=0.5d0*(dvdz(j)+dvdz(j+1))
          dtdzj=0.5d0*(dtdz(j)+dtdz(j+1))
          as2j=dudzj*dudzj+dvdzj*dvdzj
          an2j=g_alphaj*dtdzj
c         if(an2j.lt.0.d0) tauj=tauj/(1.d0-0.04d0*an2j*tauj*tauj)
          uwj=0.5d0*(uw(j)+uw(j+1))
          vwj=0.5d0*(vw(j)+vw(j+1))
          w2j=0.5d0*(w2(j)+w2(j+1))
          wtj=0.5d0*(wt(j)+wt(j+1))
          bydz=1./dze(j)
          du2dz=(u2(j+1)-u2(j))*bydz
          dv2dz=(v2(j+1)-v2(j))*bydz
          dw2dz=(w2(j+1)-w2(j))*bydz
          duvdz=(uv(j+1)-uv(j))*bydz
          duwdz=(uw(j+1)-uw(j))*bydz
          dvwdz=(vw(j+1)-vw(j))*bydz
          dutdz=(ut(j+1)-ut(j))*bydz
          dvtdz=(vt(j+1)-vt(j))*bydz
          dwtdz=(wt(j+1)-wt(j))*bydz
          dt2dz=(t2(j+1)-t2(j))*bydz
          dedz=(e(j+1)-e(j))*bydz
c         dedz=0.5d0*(du2dz+dv2dz+dw2dz)

          call tom_BUV(g_alphaj,tauj,dudzj,dvdzj,as2j,an2j,
     &      uwj,vwj,w2j,wtj,du2dz,dv2dz,dw2dz,
     &      duvdz,duwdz,dvwdz,dutdz,dvtdz,dwtdz,dt2dz,
     &      Ke(j),gc_ew(j))
          kt2(j)=0.d0
          gc_wt2(j)=0.d0
c         tmp=0.d0
c         call tom_BUV(g_alphaj,tauj,tmp,tmp,tmp,an2j,
c    &      tmp,tmp,w2j,wtj,du2dz,dv2dz,dw2dz,
c    &      tmp,tmp,tmp,tmp,tmp,dwtdz,dt2dz,
c    &      Ke(j),gc_ew(j),ew(j))
c          x=tauj*tauj*an2j
c         call tom_B(g_alphaj,tauj,x,w2j,wtj,dw2dz,dwtdz,dt2dz,dedz,
c    &        ke(j),kwt(j),kt2(j),gc_ew(j),gc_w2t(j),gc_wt2(j))
c         call tom_B(g_alphaj,tauj,x,w2j,wtj,dw2dz,dwtdz,dt2dz,dedz,
c    &        ke(j),tmp,kt2(j),gc_ew(j),tmp,gc_wt2(j))
c         call tom_B(g_alphaj,tauj,x,w2j,wtj,dw2dz,dwtdz,dt2dz,dedz,
c    &        ke(j),tmp,tmp,gc_ew(j),tmp,tmp)
          ke(j)=min(max(ke(j),kemin),kmax)
          kt2(j)=min(max(kt2(j),kt2min),kmax)
        end do

      endif

      return
      end subroutine kgcm

      subroutine find_pbl_top(e,dbll,n)
!@sum find_pbl_top Find the pbl top (at main level lpbl)
!@+   if e(j) le certain fraction of e(1), real(j) is pbl top
!@auth  Ye Cheng
!@ver   1.0
!@var dbll the (real*8) layer number corresponds to the top of the pbl

      implicit none

      integer, intent(in) :: n
      real*8, dimension(n), intent(in) :: e
      real*8, intent(out) :: dbll

      real*8, parameter :: fraction = 0.1d0
      real*8 :: e1p    ! certain percent of e_1
      integer j

      e1p=fraction*e(1)
      do j=2,n
        if (e(j).lt.e1p) exit
      end do
      dbll=j-1   ! dbll is real*8
      return
      end subroutine find_pbl_top

      subroutine ave_uv_to_agrid(u,v,u_a,v_a,im,jm,lm)
!@sum Computes u_a,v_a from u and v
!@var u x-component at secondary grids (B_grid)
!@var v y-component at secondary grids (B_grid)
!@var u_a x-component at primary grids (A_grid)
!@var v_a y-component at primary grids (A_grid)
!@auth Ye Cheng
!@ver  1.0
      USE GEOM, only : imaxj,idij,idjj,kmaxj,rapj,cosiv,siniv
      implicit none

      integer, intent(in) :: im,jm,lm
      real*8, dimension(im,jm,lm), intent(in) ::  u,v
      real*8, dimension(lm,im,jm), intent(out) :: u_a,v_a

      real*8 :: HEMI,u_t,v_t
      integer :: i,j,l,k

c     polar boxes
      DO J=1,JM,JM-1
        HEMI=1.
        IF(J.LE.JM/2) HEMI=-1.
        DO I=1,IMAXJ(J)
          DO L=1,LM
            u_t=0.d0; v_t=0.d0
            DO K=1,KMAXJ(J)
              u_t=u_t+rapj(k,j)*(u(idij(k,i,j),idjj(k,j),L)*cosiv(k)-
     2                      hemi*v(idij(k,i,j),idjj(k,j),L)*siniv(k))
              v_t=v_t+rapj(k,j)*(v(idij(k,i,j),idjj(k,j),L)*cosiv(k)+
     2                      hemi*u(idij(k,i,j),idjj(k,j),L)*siniv(k))
            END DO
            u_a(l,i,j)=u_t
            v_a(l,i,j)=v_t
          END DO
        END DO
      END DO
c     non polar boxes
      DO J=2,JM-1
        DO I=1,IMAXJ(J)
          DO L=1,LM
            u_t=0.d0; v_t=0.d0
            DO K=1,KMAXJ(J)
              u_t=u_t+u(idij(k,i,j),idjj(k,j),L)*rapj(k,j)
              v_t=v_t+v(idij(k,i,j),idjj(k,j),L)*rapj(k,j)
            END DO
            u_a(l,i,j)=u_t
            v_a(l,i,j)=v_t
          END DO
        END DO
      END DO
C****
      return
      end subroutine ave_uv_to_agrid

      subroutine ave_uv_to_bgrid(u_a,v_a,u,v,im,jm,lm)
!@sum Computes u and v from u_a and v_a
!@var u_a x-component of wind at primary grids (A_grid)
!@var v_a y-component of wind at primary grids (A_grid)
!@var u x-component of wind at secondary grids (B_grid)
!@var v y-component of wind at secondary grids (B_grid)
!@auth Ye Cheng
!@ver  1.0

      USE GEOM, only : imaxj,idij,idjj,kmaxj,ravj,cosiv,siniv

      implicit none

      integer, intent(in) :: im,jm,lm
      real*8, dimension(lm,im,jm), intent(in) ::  u_a,v_a
      real*8, dimension(im,jm,lm), intent(out) :: u,v

      real*8 :: HEMI
      integer :: i,j,l,k

      u=0.d0; v=0.d0
c     polar boxes
      DO J=1,JM,JM-1
        HEMI=1.
        IF(J.LE.JM/2) HEMI=-1.
        DO I=1,IMAXJ(J)
        DO L=1,LM
        DO K=1,KMAXJ(J)
          U(IDIJ(K,I,J),IDJJ(K,J),L)=U(IDIJ(K,I,J),IDJJ(K,J),L)
     *      +RAVJ(K,J)*(U_A(L,I,J)*COSIV(K)+V_A(L,I,J)*SINIV(K)*HEMI)
          V(IDIJ(K,I,J),IDJJ(K,J),L)=V(IDIJ(K,I,J),IDJJ(K,J),L)
     *      +RAVJ(K,J)*(V_A(L,I,J)*COSIV(K)-U_A(L,I,J)*SINIV(K)*HEMI)
        END DO
        END DO
        END DO
      END DO
c     non polar boxes
      DO J=2,JM-1
        DO I=1,IMAXJ(J)
        DO L=1,LM
        DO K=1,KMAXJ(J)
          U(IDIJ(K,I,J),IDJJ(K,J),L)=U(IDIJ(K,I,J),IDJJ(K,J),L)
     *          +RAVJ(K,J)*U_A(L,I,J)
          V(IDIJ(K,I,J),IDJJ(K,J),L)=V(IDIJ(K,I,J),IDJJ(K,J),L)
     *          +RAVJ(K,J)*V_A(L,I,J)
        END DO
        END DO
        END DO
      END DO

      return
      end subroutine ave_uv_to_bgrid

      subroutine ave_s_to_bgrid(s,s_b,im,jm,lm)
!@sum Computes s_b from s
!@var s scalar at primary grid (A_grid)
!@var s_b scalar at secondary grid (B_grid)
!@auth Ye Cheng
!@ver  1.0

      USE GEOM, only : imaxj,idij,idjj,kmaxj,ravj

      implicit none

      integer, intent(in) :: im,jm,lm
      real*8, dimension(lm,im,jm), intent(in) ::  s
      real*8, dimension(lm,im,jm), intent(out) :: s_b

      integer :: i,j,l,k

      S_B=0.d0
      DO J=1,JM
        DO I=1,IMAXJ(J)
        DO L=1,LM
        DO K=1,KMAXJ(J)
          S_B(L,IDIJ(K,I,J),IDJJ(K,J))=S_B(L,IDIJ(K,I,J),IDJJ(K,J))
     *          +RAVJ(K,J)*S(L,I,J)
        END DO
        END DO
        END DO
      END DO
1003  format(4(i4,1x),3(1pe14.4))
      return
      end subroutine ave_s_to_bgrid

      subroutine tom_BUV(ga,tau,dudz,dvdz,as2,an2,uw,vw,w2,wt0,
     &   du2dz,dv2dz,dw2dz,duvdz,duwdz,dvwdz,dutdz0,dvtdz0,dwtdz0,
     &   dt2dz0,Ke,gc_ew)
c
c      q2w  = - Ke*dq2dz + gc_q2w
c      ew  = - Ke*dedz + gc_ew, gc_ew = gc_q2w/2
c output of 3m_eqns,3m_solve_sb0,3m_solve_sb0_more,more2,more3,more31,
c more32, more33 on kirk:/u/acyxc/papers/third/3m_publication
c tau=q2/epsilon=B1*ell/q, ell->0.4*z for small z (height)
c each t obsorbs a lamda=g*alpha*tau, ga==g*alpha       --- 04/6/00
c
      USE CONSTANT, only : by3

      implicit none

      real*8, intent(in) :: ga,tau,dudz,dvdz,as2,an2,uw,vw,w2
     &     ,du2dz,dv2dz,dw2dz,duvdz,duwdz,dvwdz
     &     ,wt0,dutdz0,dvtdz0,dwtdz0,dt2dz0
      real*8, intent(out) :: Ke,gc_ew

      real*8, parameter :: c=1.d0/8.d0
      real*8, SAVE :: c0,c1,c2,c3,c4,c5,c6,c7,c8,c9
     &    ,c10,c11,c12,c13,c14,c15
     &    ,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27
     &    ,d0,d1,d2,d3,d4,d5,d6,d7,d8

      real*8 :: U,V,n0,n1,n2,s1,s2,b1,b2,b3,b4,b5,d
     &    ,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,num_rest,gc_q2w
     &    ,tmp,wt,dutdz,dvtdz,dwtdz,dt2dz,dq2dz,dedz
      integer, SAVE :: ifirst=0

      U=dudz*tau
      V=dvdz*tau
      S2=as2*tau*tau
      N2=an2*tau*tau
      tmp=ga*tau
      wt    = tmp * wt0
      dutdz = tmp * dutdz0
      dvtdz = tmp * dvtdz0
      dwtdz = tmp * dwtdz0
      dt2dz = tmp*tmp * dt2dz0

      dq2dz=du2dz+dv2dz+dw2dz
      dedz=0.5d0*dq2dz

      if (ifirst.eq.0) then

          ifirst=1
          c0 = -2*c**2*(13*c+9)
          c1 = -900
          c2 = -30*c*(44*c + 3)
          c3 = -54*c**3
          c4 = 25*c
          c5 = 5*c**3
          c6 = 25*c + 15
          c7 =  5*c**3
          c8 =  6*c**2*(3*c+2)
          c9   =  1500*(5*c+3)
          c10   = 750*c**2*(7*c+3)
          c11   = 50*c*(268*c**2+183*c+9)
          c12   = 25*c**3*(146*c**2+63*c+9)
          c13   = 30*c**3*(242*c**2+181*c+15)
          c14   = 30*c**5*(57*c**2+28*c-3)
          c15   = 36*c**5*(3*c+2)*(10*c+1)
          c16   = 4500*(5*c+3)
          c17   = 1500*c**2*(17*c+9)
          c18   = 50*c*(394*c**2+303*c+27)
          c19   = 50*c**3*(-17*c**2+24*c+27)
          c20   = 30*c**3*( 84*c**2+125*c+45 )
          c21   = 45*c**5*(7*c+3)
          c22   = 108*c**5*(3*c+2)
          c23   = 2050*c+1230
          c24   = 90*c**2*(29*c+16)
          c25   = 6*c**2*(321*c+209)
          c26   = 9*c**4*(57*c+31)
          c27   = 108*c**4*(3*c+2)

          d0 = 7500*(5*c+3)**2
          d1 = 3750*c**3*(23*c+18)
          d2 = 250*c*(5*c+3)*(292*c**2+237*c+9)
          d3 = 125*c**4*(526*c**2+645*c+54)
          d4 = 150*c**3*(1254*c**3+2081*c**2+933*c+72)
          d5 = 1350*c**6*(19*c**2+17*c+1)
          D6 = 90*c**5*(180*c**3+438*c**2+299*c+57)
          d7 = 405*c**8
          d8 = 324*c**7*(3*c+2)
      endif

      n0 = c1 + c2*N2 + c3*N2**2
      n1 = c4 + c5*N2
      s1 = -2*c5*U*V
      b1 = c6 +  c7*V**2 + c8*N2
      b2 = c6 +  c7*U**2 + c8*N2
      b3 = c9+c10*S2+(c11+c12*S2)*N2+(c13+c14*S2)*N2**2+c15*N2**3
      b4 = c16+c17*S2+(c18+c19*S2)*N2+(c20+c21*S2)*N2**2+c22*N2**3
      b5 = c23 + c24*S2 + (c25 + c26*S2)*N2 + c27*N2**2
      D = (d0 + d1*S2) + (d2 + d3*S2)*N2 + (d4 + d5*S2)*N2**2 +
     &    (d6 + d7*S2)*N2**3 + d8*N2**4

      m1 = n0*(n1*c*(3*U*uw+V*vw)+b1*(5*w2+2*c*wt))
      m2 = n0*(n1*c*(3*V*vw+U*uw)+b2*(5*w2+2*c*wt))

      m3 = -n0*n1*(4*c+3)*( U*uw+V*vw )-45*b3*w2-6*b4*c*wt
      m4 = n0*(2*n1*c*(V*uw+U*vw)+s1*(5*w2+2*c*wt))

      m5 = n0*(10*b1*uw+5*s1*vw-2*n1*(4*c+3)*U*w2+5*c0*U*wt)
      m6 = n0*(10*b2*vw+5*s1*uw-2*n1*(4*c+3)*V*w2+5*c0*V*wt)

      m7 = n0*(2*c*(2*b1*uw+s1*vw)+c0*U*(5*w2+4*c*wt))
      m8 = n0*(2*c*(2*b2*vw+s1*uw)+c0*V*(5*w2+4*c*wt))

      m9 = 5*n0*c0*(U*uw+V*vw)+12*c*(-b4*w2-10*b5*c*wt)
      M10 = 2*c*(n0*c0*(U*uw+V*vw)-6*c*b5*(5*w2+wt))

c     Num=m1*du2dz+m2*dv2dz+m3*dw2dz+m4*duvdz+m5*duwdz+
c         m6*dvwdz+m7*dutdz+m8*dvtdz+m9*dwtdz+m10*dt2dz
c        =m1*(du2dz-dq2dz/3)+m2*(dv2dz-dq2dz/3)+m3*(dw2dz-dq2dz/3)
c        +m4*duvdz+m5*duwdz+m6*dvwdz+m7*dutdz+m8*dvtdz
c        +m9*dwtdz+m10*dt2dz+(m1+m2+m3)/3*dq2dz
c     q2w = tau*c/(2*D)*Num
c         = - Ke*dq2dz + gc_q2w
c     ew  = - Ke*dedz + gc_ew , where gc_ew=1/2*gc_q2w
      Ke=-tau*c/(6.*D)*(m1+m2+m3)
      gc_q2w=tau*c/(2.*D)*
     & (m1*(du2dz-dq2dz*by3)+m2*(dv2dz-dq2dz*by3)+m3*(dw2dz-dq2dz*by3)
     & +m4*duvdz+m5*duwdz+m6*dvwdz+m7*dutdz+m8*dvtdz
     & +m9*dwtdz+m10*dt2dz)
      gc_ew=0.5d0*gc_q2w
c     ew=-Ke*dedz+gc_ew

      return
      end subroutine tom_BUV

      subroutine tom_B(ga,tau,x,w2,wt0,dw2dz,dwtdz0,dt2dz0,dedz,
     2  Ke,Kwt,Kt2,gc_ew,gc_w2t,gc_wt2)
c
c     algebraic third moment expressions in terms of
c     first and second moments with pure bouynacy
c     c is adjustable from 1/5 to 1/9
c     input:
c         tau = 2*TKE/epsilon
c         x    = g*alpha*dT/dz*tau**2 = (N*tau)**2
c         w2,dw2dz,dedz are in the usual sense
c         wt0,dwtdz0,dt2dz0 are wt,dwtdz,dt2dz
c     output:
c         Ke,Kwt,Kt2 diffusivities
c         gc_ew,gc_w2t,gc_wt2 the rest of the tom expressions
c         ew = -Ke*d/dz e + gc_ew
c         w2t = -Kwt*d/dz <wt> + gc_w2t
c         wt2 = -Kt2*d/dz <t2> + gc_wt2
c     date:
c         12-06-01
c     reference:
c       /u/acyxc/papers/third/with_V/
c       3rdm_eqns,3rdm_solve_B,3rdm_solve_B_more,test_tom_B.f
c
      USE CONSTANT, only : by3
      implicit none

      real*8, parameter :: c=1.d0/8.d0
      real*8, intent(in) :: ga,tau,x,w2,wt0,dw2dz,dwtdz0,dt2dz0,dedz
      real*8, intent(out) :: Ke,Kwt,Kt2,gc_ew,gc_w2t,gc_wt2

      real*8, save :: d0,d1,d2,d3,
     & a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,
     & a13,a14,a15,a16,a17,a18,a19,
     & b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,
     & c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,twoby3
      real*8 :: A_w2,B_w2,C_w2,D_w2,A_wt,B_wt,C_wt,D_wt,
     & x2,tau_by_Den,tmp,wt,dwtdz,dt2dz
      integer, save :: ifirst=0

c     wt    == g*alpha*tau*<wt>, dimension of <w2>
c     dwtdz == g*alpha*tau*d<wt>/dz, dimension of d<w2>/dz
c     dt2dz == (g*alpha*tau)**2*d<t2>/dz, dimension of d<w2>/dz
      tmp=ga*tau
      wt=wt0*tmp
      dwtdz=dwtdz0*tmp
      dt2dz=dt2dz0*tmp*tmp

      if(ifirst.eq.0) then ! calc. for the constants only once

        ifirst=1

        twoby3=2.d0*by3
        d0 = 1500*(5*c+3)*(2*c+1)
        d1 = 25*c*(345*c+216*c**3+656*c**2+27)
        d2 = 15*c**3*(132*c**2+36+115*c)
        d3 = 81*c**5

        a0 = -2250*c*(2*c+1)
        a1 = -75*c**2*(52*c+9)*0.5d0
        a2 = -675*c**5
        a3 = -900*c**2*(2*c+1)
        a4 = 45*c**3*(10*c-3)
        a5 = -2700*c**2*(2*c+1)
        a6 = -15*c**3*(28*c+72*c**2+27)
        a7 = -81*c**5
        a8 = -2460*c**3
        a9 = -540*c**5
        a10 = -1230*c**3
        a11 = -270*c**5
        a12 = -369*c**3
        a13 = -81*c**5

        a14 = -2250*c*(2*c+1)
        a15 = -75.*0.5d0*c**2*(9+72*c**2+88*c)
        a16 = -405.*0.5d0*c**4
        a17 = -900*c**2*(2*c+1)
        a18 = -15*c**3*(9+72*c**2+88*c)
        a19 = -81*c**5

        b0 = 3375*c**2*(2*c+1)*(c+1)
        b1 = 675.*0.25d0*c**3*(8*c**2+7*c+3)
        b2 = 405.*0.25d0*c**5
        b3 = -750*c*(2*c+1)*(5*c+3)
        b4 = -450.*0.25d0*c**2*(8*c**2+9*c+3)
        b5 = -270.*0.25d0*c**4
        b6 = -1500*c*(2*c+1)*(5*c+3)
        b7 = -225*c**2*(3+12*c**2+8*c**3+9*c)
        b8 = -135*c**4*(c+1)
        b9 = -1300*c**2*(5*c+3)
        b10 = -60*c**4*(13+15*c)
        b11 = 0.5d0
        b12 = 0.15d0
        b13 = -2250*c**3*(2*c+1)
        b14 = -675.*0.5d0*c**4
        b15 = 0.4d0*c

        c0 = -3375*c**3*(c+1)
        c1 = -675*c**5
        c2 = 750*c**2*(5*c+3)
        c3 = 450*c**4
        c4 = 1500*c**2*(5*c+3)
        c5 = 900*c**4*(c+1)
        c6 = -1500*c*(5*c+3)
        c7 = -900*c**3*(4+3*c)
        c8 = -540*c**5
        c9 = 0.5d0
        c10 = 0.15d0
        c11 = 2250*c**4
        c12 = 900*c**5

      endif

      X2=x*x
      tau_by_Den = tau/(d0+d1*x+d2*x2+d3*x2*x)

c ew:

      A_W2 = a0+a1*x+a2*x2
      A_wt = a3+a4*x
      B_w2 = a5+a6*x+a7*x2
      B_wt = a8+a9*x
      C_w2 = a10+a11*x
      C_wt = a12+a13*x
      D_w2 = a14+a15*x+a16*x2
      D_wt = a17+a18*x+a19*x2

c     Ke = -tau_by_Den*(D_w2*w2+D_wt*wt)
c     ew = tau_by_Den*
c    &  ((A_w2*w2+A_wt*wt)*dw2dz+(B_w2*w2+B_wt*wt)*dwtdz
c    &  +(C_w2*w2+C_wt*wt)*dt2dz+(D_w2*w2+D_wt*wt)*dedz)
c        = tau_by_Den*
c    &  ((A_w2*w2+A_wt*wt)*(dw2dz-twoby3*dedz)+(B_w2*w2+B_wt*wt)*dwtdz
c    &  +(C_w2*w2+C_wt*wt)*dt2dz+(D_w2*w2+D_wt*wt)*dedz
c    &    +2./3*(A_w2*w2+A_wt*wt)*dedz)

c     Ke=-tau_by_Den*(D_w2*w2+D_wt*wt+twoby3*(A_w2*w2+A_wt*wt))
      Ke=-tau_by_Den*((D_w2+twoby3*A_w2)*w2+(D_wt+twoby3*A_wt)*wt)
      gc_ew=tau_by_Den*
     &  ((A_w2*w2+A_wt*wt)*(dw2dz-twoby3*dedz)+(B_w2*w2+B_wt*wt)*dwtdz
     &  +(C_w2*w2+C_wt*wt)*dt2dz)
c     ew=-Ke*dedz+gc_ew

c w2t:

c     A_w2 = (b0+b1*x+b2*x2)*x
c     A_wt = b3+b4*x+b5*x2
c     B_w2 = b6+b7*x+b8*x2
c     B_wt = b9+b10*x
c     C_w2 = b11*B_wt
c     C_wt = b12*B_wt
c     D_w2 = (b13+b14*x)*x
c     D_wt = b15*D_w2

c     Kwt = -tau_by_Den*(B_w2*w2+B_wt*wt)
c     gc_w2t = tau_by_Den*
c    &  ((A_w2*w2+A_wt*wt)*dw2dz
c    &  +(C_w2*w2+C_wt*wt)*dt2dz+(D_w2*w2+D_wt*wt)*dedz)
c     gc_w2t = gc_w2t/tmp
cc    w2t = -Kwt*dwtdz0 + gc_w2t
      Kwt=0.
      gc_w2t=0.

c wt2:

      A_w2 = (c0+c1*x)*x2
      A_wt = (c2+c3*x)*x
      B_w2 = (c4+c5*x)*x
      B_wt = c6+c7*x+c8*x2
      C_w2 = c9*B_wt
      C_wt = c10*B_wt
      D_w2 = c11*x2
      D_wt = c12*x2

      Kt2 = -tau_by_Den*(C_w2*w2+C_wt*wt)
      gc_wt2 = tau_by_Den*
     &  ((A_w2*w2+A_wt*wt)*dw2dz+(B_w2*w2+B_wt*wt)*dwtdz
     &  +(D_w2*w2+D_wt*wt)*dedz)
      gc_wt2=gc_wt2/(tmp*tmp)
c     wt2 = -Kt2*dt2dz0 + gc_wt2

      return
      end subroutine tom_B

      subroutine apply_fluxes_to_atm
!@sum dummy subroutine - replaces the real one needed by DRYCNV
!@auth I. Aleinov
!@ver  1.0
      return
      end subroutine apply_fluxes_to_atm
