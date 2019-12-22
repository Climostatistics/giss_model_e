#include "rundeck_opts.h"

#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
!@sum  OCN_TRACER: tracer-dependent routines for GISS Ocean tracers
!@+    Routines included:
!@+      Those that MUST EXIST for all tracers:
!@+        Diagnostic specs: init_tracer_ocean (if TRACERS_OCEAN true)
!@+        Tracer initialisation + sources: tracer_ic_ocean
!@+
!@auth Jean Lerner/Gavin Schmidt
!@ver  1.0

      subroutine tracer_ic_ocean
!@sum tracer_ic_ocean initialise ocean tracers
!@auth Gavin Schmidt
!@ver 1.0
      USE MODEL_COM, only: itime
      USE TRACER_COM, only : itime_tr0, ntm, trname, trw0
      USE OCEAN, only : im,jm,lmo,focean,dxypo,mo,lmm,imaxj
#ifdef TRACERS_OCEAN
     *     ,trmo,txmo,tymo,tzmo,s0m,sxmo,symo,szmo
#endif
      USE SEAICE, only : xsi,lmi
      USE STRAITS, only : nmst,msist,ssist
#ifdef TRACERS_OCEAN
     *     ,lmst,ist,jst,xst,yst,mmst,s0mst,sxmst,szmst,trmst,txmst
     *     ,tzmst
#endif
#ifdef TRACERS_WATER
     *     ,trsist
      USE FLUXES, only : gtracer
#endif
      USE FILEMANAGER, only : openunit,closeunit
      IMPLICIT NONE
      integer n,i,j,l,k,nst,i1,j1,i2,j2
      real*8 t01,t02
#ifdef TRACERS_SPECIAL_O18
      integer iu_O18ic,ip1,im1
      character*80 title
      real*4, dimension(im,jm,lmo) :: t0m4,tzm4
      real*8 fwm,trsum,wsum,tdiff,afac
#endif


C**** Note that only sea ice related arrays are initialised if
C**** only TRACERS_WATER is true.

      do n=1,ntm

        if (itime.eq.itime_tr0(n)) then
        select case (trname(n))

        case default
#ifdef TRACERS_OCEAN
          trmo(:,:,:,n)=0.
          txmo(:,:,:,n)=0.
          tymo(:,:,:,n)=0.
          tzmo(:,:,:,n)=0.
#endif
#ifdef TRACERS_WATER
          do j=1,jm
          do i=1,im
            if (focean(i,j).gt.0) gtracer(n,1,i,j)=0.
          end do
          end do
#endif
        case ('Water', 'H2O18', 'HDO')
#if (defined TRACERS_OCEAN) && (defined TRACERS_SPECIAL_O18)
C**** Open ic file for isotope tracers
          call openunit("H2O18ic",iu_O18ic,.true.,.true.)
#endif

#ifdef TRACERS_OCEAN
#ifndef TRACERS_SPECIAL_O18
C**** main ocean variabiles and gradients
          do j=1,jm
          do i=1,im
            do l=1,lmm(i,j)
              trmo(i,j,l,n)=trw0(n)*(mo(i,j,l)*dxypo(j)-s0m(i,j,l))
              txmo(i,j,l,n)=-trw0(n)*sxmo(i,j,l)
              tymo(i,j,l,n)=-trw0(n)*symo(i,j,l)
              tzmo(i,j,l,n)=-trw0(n)*szmo(i,j,l)
            end do
          end do
          end do
#else
C**** read in initial conditions for isotopes
C**** search through for correct tracer (since there is no guarantee
C**** that they will be in same order as tracer indices).
C**** data are now in 'per mil'
          rewind (iu_O18ic)
 10       read  (iu_O18ic,err=800,end=810) title,t0m4,tzm4
          if (index(title,trim(trname(n))).eq.0) goto 10
          write (6,*) 'Read from H2O18ic: ',title
          call closeunit(iu_O18ic)

C**** Turn per mil data into mass ratios (using current standard)
          t0m4(:,:,:)=(t0m4(:,:,:)*1d-3+1.)*trw0(n)
          tzm4(:,:,:)=0.    ! tzm4(:,:,:)*1d-3*trw0(n) corrupted?
C****
          do l=1,lmo
            txmo(:,:,l,n) = 0.
            tymo(:,:,l,n) = 0.
C**** Define East-West horizontal gradients
            im1=im-1
            i=im
            do j=2,jm-1
              do ip1=1,im
                if (lmm(i,j).ge.l) then
                  if (lmm(im1,j).ge.l) then
                    if (lmm(ip1,j).ge.l) then
                      txmo(i,j,l,n)=2.5d-1*(t0m4(ip1,j,l)-t0m4(im1,j,l))
                    else
                      txmo(i,j,l,n)=  5d-1*(t0m4(  i,j,l)-t0m4(im1,j,l))
                    end if
                  else
                    if (lmm(ip1,j).ge.l)
     *                   txmo(i,j,l,n)=5d-1*(t0m4(ip1,j,l)-t0m4(i,j,l))
                  end if
                end if
                im1=i
                i=ip1
              end do
            end do
C**** Define North-South horizontal gradients
            do j=2,jm-1
              do i=1,im
                if (lmm(i,j).ge.l)  then
                  if (lmm(i,j-1).ge.l)  then
                    if (lmm(i,j+1).ge.l)  then
                      tymo(i,j,l,n)=2.5d-1*(t0m4(i,j+1,l)-t0m4(i,j-1,l))
                    else
                      tymo(i,j,l,n)=  5d-1*(t0m4(i,  j,l)-t0m4(i,j-1,l))
                    end if
                  else
                    if (lmm(i,j+1).ge.l)
     *                   tymo(i,j,l,n)=5d-1*(t0m4(i,j+1,l)-t0m4(i,j,l))
                  end if
                end if
              end do
            end do
C**** Multiply ratios by freshwater mass
            do j=1,jm
            do i=1,im
            if (focean(i,j).gt.0) then
              if (l.le.lmm(i,j)) then
                fwm = mo(i,j,l)*dxypo(j)-s0m(i,j,l)
                trmo(i,j,l,n)=t0m4(i,j,l)*fwm
                txmo(i,j,l,n)=txmo(i,j,l,n)*fwm-sxmo(i,j,l)*t0m4(i,j,l)
                tymo(i,j,l,n)=tymo(i,j,l,n)*fwm-symo(i,j,l)*t0m4(i,j,l)
                tzmo(i,j,l,n)=tzm4(i,j,l)  *fwm-szmo(i,j,l)*t0m4(i,j,l)
              else
                trmo(i,j,l,n)=0.
                txmo(i,j,l,n)=0.
                tymo(i,j,l,n)=0.
                tzmo(i,j,l,n)=0.
              end if
            end if
            end do
            end do
          end do
C**** Balance tracers so that average concentration is TRW0
          trsum = 0
          do j=2,jm
            afac = 1
            if (j.eq.jm) afac = im
            do i=1,imaxj(j)
              do l=1,lmm(i,j)
                trsum = trsum + afac*focean(i,j)*trmo(i,j,l,n)
              end do
            end do
          end do
          if (trname(n).eq.'Water') then
            wsum = trsum
          else
            tdiff = trsum - wsum*trw0(n)
            write(6,*) "Average oceanic tracer concentration ",
     *           trname(n),(trsum/(wsum*trw0(n))-1d0)*1d3,tdiff,trsum
     *           ,wsum
            do l=1,lmo
              do j=2,jm
                do i=1,imaxj(j)
                  trmo(i,j,l,n) = trmo(i,j,l,n) * (1d0 - tdiff/trsum)
                end do
              end do
            end do
          end if
#endif
          do j=1,jm
          do i=1,im
#ifdef TRACERS_WATER
            if (focean(i,j).gt.0) gtracer(n,1,i,j)=trmo(i,j,1,n)/(mo(i,j
     *           ,1)*dxypo(j)-s0m(i,j,1))
#endif
          do l=lmm(i,j)+1,lmo
            trmo(i,j,l,n)=0. ; txmo(i,j,l,n)=0.
            tymo(i,j,l,n)=0. ; tzmo(i,j,l,n)=0.
          end do
          end do
          end do
#endif

C**** Initiallise strait values based on adjacent ocean boxes
          do nst=1,nmst
#ifdef TRACERS_OCEAN
            i1=ist(nst,1)
            j1=jst(nst,1)
            i2=ist(nst,2)
            j2=jst(nst,2)
            do l=1,lmst(nst)
              t01=trmo(i1,j1,l,n)/(mo(i1,j1,l)*dxypo(j1)-s0m(i1,j1,l))
              t02=trmo(i2,j2,l,n)/(mo(i2,j2,l)*dxypo(j2)-s0m(i2,j2,l))
              trmst(l,nst,n) = 5d-1*(mmst(l,nst)-s0mst(l,nst))*(t01+t02)
              txmst(l,nst,n) = -5d-1*sxmst(l,nst)*(t01+t02)
              tzmst(l,nst,n) = -5d-1*szmst(l,nst)*(t01+t02)
            end do
            do l=lmst(nst)+1,lmo
              trmst(l,nst,n) = 0.
              txmst(l,nst,n) = 0.
              tzmst(l,nst,n) = 0.
            end do
#endif
#ifdef TRACERS_WATER
            trsist(n,1:2,nst) = trw0(n)*(msist(1,nst)*xsi(1:2)
     *           -ssist(1:2,nst))
            trsist(n,3:lmi,nst)=trw0(n)*(msist(2,nst)*xsi(3:lmi)
     *           -ssist(3:lmi,nst))
#endif
          end do
C****
        end select
        write(6,*) trname(n)," tracer initialised in ocean"
        end if
      end do

      return
 800  write(6,*) "Error reading input file H2O18ic"
      call stop_model('stopped in OCN_TRACER.f',255)
 810  write(6,*) "Tracer ",trname(n)," not found in file H2O18ic"
      call stop_model('stopped in OCN_TRACER.f',255)
      end subroutine tracer_ic_ocean
#endif

#ifdef TRACERS_OCEAN
      subroutine init_tracer_ocean
!@sum Initialise ocean tracer diagnostics (none so far)
      return
      end subroutine init_tracer_ocean
#endif

#ifdef TRACERS_OCEAN
      SUBROUTINE OC_TDECAY
!@sum OC_TDECAY decays radioactive tracers in ocean
!@auth Gavin Schmidt/Jean Lerner
      USE MODEL_COM, only : itime,dtsrc
      USE TRACER_COM, only : ntm,trdecay,itime_tr0
      USE OCEAN, only : trmo,txmo,tymo,tzmo
      IMPLICIT NONE
      real*8, save, dimension(ntm) :: expdec = 1.
      logical, save :: ifirst=.true.
      integer n

      if (ifirst) then               
        do n=1,ntm
          if (trdecay(n).gt.0.0) expdec(n)=exp(-trdecay(n)*dtsrc)
        end do
        ifirst = .false.
      end if

      do n=1,ntm
        if (trdecay(n).gt.0. .and. itime.ge.itime_tr0(n)) then
C**** Oceanic decay
          trmo(:,:,:,n)   = expdec(n)*trmo(:,:,:,n)
          txmo(:,:,:,n)   = expdec(n)*txmo(:,:,:,n)
          tymo(:,:,:,n)   = expdec(n)*tymo(:,:,:,n)
          tzmo(:,:,:,n)   = expdec(n)*tzmo(:,:,:,n)
        end if
      end do
C****
      return
      end subroutine oc_tdecay
#endif
