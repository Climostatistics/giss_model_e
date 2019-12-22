c     Family chemistry calculations:
c     equil values are production/loss from rxns within family only
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
      SUBROUTINE Oxinit(lmax,I,J)
!@sum Oxinit Find O,O1D and Ox initial conc assuming equilibrium with O3
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on ds3ch4_famchem_apr1902_M23)
c
C**** GLOBAL parameters and variables:
c
      USE TRACER_COM, only : n_CH4, n_Ox
      USE TRCHEM_Shindell_COM, only:ss,rr,y,nO2,nM,nH2O,nO,nO1D,nO3,pOx
C
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
C
!@var az,bz,P1 dummy variables
!@var L dummy loop variable
!@var I,J passed horizontal position indicies
!@var lmax maximum altitude for chemistry (usually LS1-1 ~ tropopause)
      real*8 az, bz, P1
      integer, intent(IN) :: lmax,I,J
      integer L
C
      do L=1,lmax
c       for concentration of O:
        az=(ss(2,I,J,L)+ss(3,I,J,L))/(rr(47,I,J,L)*y(nO2,L))
c       for concentration of O(1D):
        bz=ss(2,I,J,L)/
     &  (rr(8,I,J,L)*y(nO2,L)+rr(9,I,J,L)*y(nM,L)+
     &  rr(10,I,J,L)*y(nH2O,L)+rr(11,I,J,L)*y(n_CH4,L))
        P1=1./(1.+az+bz)
        y(nO,L)=P1*az*y(n_Ox,L)
        y(nO1D,L)=P1*bz*y(n_Ox,L)
        y(nO3,L)=y(n_Ox,L)-y(nO,L)-y(nO1D,L)
        if(y(nO,L).lt.0.)y(nO,L)=0.0
        if(y(nO1D,L).lt.0.)y(nO1D,L)=0.0
        if(y(nO3,L).lt.1.)y(nO3,L)=1.0
        if(y(n_Ox,L).lt.1.)y(n_Ox,L)=1.0
        pOx(I,J,L)=y(nO3,L)/y(n_Ox,L)
      enddo
c
      return
      END SUBROUTINE OXinit
C
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
      SUBROUTINE NOxfam(lmax,I,J)
!@sum NOxfam Find NOx family (NO,NO2,NO3,HONO) partitioning assuming
!@+   equilibrium at a given concentration of NOx. Only called during
!@+   daylight, then assume NO3=HONO=1.
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on ds3ch4_famchem_apr1902_M23)
c
C**** GLOBAL parameters and variables:
c
      USE TRACER_COM, only : n_NOx
      USE TRCHEM_Shindell_COM, only:rr,y,yNO3,nO3,nHO2,yCH3O2,nO,nC2O3,
     &                           ta,nXO2,ss,nNO,nNO2,pNOx,nNO3,nHONO
C
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var b,c,p1,p2 dummy variables
!@var L dummy loop variable
!@var I,J passed horizontal position indicies
!@var lmax maximum altitude for chemistry (usually LS1-1 ~ tropopause)
      real*8 b,c,p1,p2
      integer L
      integer, intent(IN) :: lmax,I,J

      do L=1,lmax
c       If dawn then set NO3 back to zero:
        IF(yNO3(I,J,L).GT.0.)yNO3(I,J,L)=0.
c       B is for NO->NO2 reactions :
        B=rr(5,I,J,L)*y(nO3,L)+rr(6,I,J,L)*y(nHO2,L)
     &   +rr(20,I,J,L)*yCH3O2(I,J,L) + rr(48,I,J,L)*y(nO,L)
     &   +rr(39,I,J,L)*y(nC2O3,L)+4.2E-12*exp(180./ta(L))*y(nXO2,L)
C       C is for NO2->NO reactions :
        C=ss(1,I,J,L)
        p2=B/(B+C)
        p1=1-p2
        y(nNO,L)=p1*y(n_NOx,L)
        y(nNO2,L)=p2*y(n_NOx,L)
C       Set limits on NO, NO2, NOx:
        if(y(nNO,L).lt.1.) y(nNO,L) = 1.0
        if(y(nNO2,L).lt.1.)y(nNO2,L)= 1.0
        if(y(n_NOx,L).lt.1.)y(n_NOx,L)= 1.0
        pNOx(I,J,L)=y(nNO2,L)/y(n_NOx,L)
        y(nNO3,L) =1.0
        y(nHONO,L)=1.0
      enddo

      return
      END SUBROUTINE NOxfam
C
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
      SUBROUTINE HOxfam(lmax,I,J)
!@sum HOxfam Find HOx family (OH,HO2) partitioning assuming equilibrium
!@+   concentration of HOx.
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on ds3ch4_famchem_apr1902_M23)
c
C**** GLOBAL parameters and variables:
c
      USE TRACER_COM, only : n_CH4,n_HNO3,n_CH3OOH,n_H2O2,n_HCHO,n_CO,
     &                       n_Paraffin,n_Alkenes,n_Isoprene,n_AlkylNit
      USE TRCHEM_Shindell_COM, only:pHOx,rr,y,nNO2,nNO,yCH3O2,nH2O,nO3,
     &                           nO2,nM,nHO2,nOH,nH2,nAldehyde,nXO2,
     &                           nXO2N,ta,ss,nC2O3,nROR
C
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var aqqz,bqqz,cqqz,cz,dz,sqroot,temp_yHOx dummy variables
!@var L dummy loop variable
!@var I,J passed horizontal position indicies
!@var lmax maximum altitude for chemistry (usually LS1-1 ~ tropopause)
      real*8 aqqz, bqqz, cqqz, cz, dz, sqroot, temp_yHOx
      integer L
      integer, intent(IN) :: lmax,I,J
C
      do L=1,lmax   ! >> beginning of altitude loop <<
c
c      First calculate equilibrium amount of HOx
c      A: loss rxns with HOx**2, B: loss rxns linear in HOx, C: prod
c       equations are in terms of HO2 (so *pHOx when OH is reactant)
C
C >>> Drew: temporarily, I omitted all your debugging writes in this
C subroutine. If you need them back, let me know.  Thanks, Greg <<<
C
       aqqz=2.*(pHOx(I,J,L)*rr(1,I,J,L) + pHOx(I,J,L)*pHOx(I,J,L)*
     & (rr(3,I,J,L)+rr(49,I,J,L)) + rr(15,I,J,L))
C
       bqqz=pHOx(I,J,L)*(rr(12,I,J,L)*y(n_CH4,L)+rr(16,I,J,L)*
     & y(n_HNO3,L)+rr(23,I,J,L)*y(n_CH3OOH,L)+rr(50,I,J,L)
     & *y(nNO2,L)+rr(53,I,J,L)*y(nNO,L))+rr(22,I,J,L)*yCH3O2(I,J,L)
     & +pHOx(I,J,L)*(rr(38,I,J,L)*y(nAldehyde,L)+rr(37,I,J,L)
     & *y(n_Paraffin,L)*0.89+rr(34,I,J,L)*y(n_Alkenes,L)
     & +rr(30,I,J,L)*y(n_Isoprene,L)*0.15+rr(33,I,J,L)*y(n_AlkylNit,L))
     & +rr(43,I,J,L)*y(nXO2,L)+y(nXO2N,L)*
     & (rr(44,I,J,L)*rr(43,I,J,L)/(4.2E-12*exp(180./ta(L))))
C
       cqqz=(2.*ss(4,I,J,L)*y(n_H2O2,L)+ss(9,I,J,L)*y(n_HNO3,L)+
     & ss(13,I,J,L)*y(n_HCHO,L)+ss(14,I,J,L)*y(n_CH3OOH,L)+
     & (rr(20,I,J,L)*y(nNO,L)+0.66*rr(27,I,J,L)*yCH3O2(I,J,L))
     & *yCH3O2(I,J,L))
C
       cqqz=cqqz+
     & ((2.*rr(10,I,J,L)*y(nH2O,L)+rr(11,I,J,L)*y(n_CH4,L))*
     & ss(2,I,J,L)*y(nO3,L))/
     & (rr(8,I,J,L)*y(nO2,L)+rr(9,I,J,L)*y(nM,L)+
     & rr(10,I,J,L)*y(nH2O,L)+rr(11,I,J,L)*y(n_CH4,L))
     & +ss(16,I,J,L)*y(nAldehyde,L)*2.+(rr(39,I,J,L)*y(nNO,L)
     & +rr(40,I,J,L)*y(nC2O3,L)*2.)*y(nC2O3,L)
     & +(rr(42,I,J,L)*0.94+1.6E3)*y(nROR,L)+rr(35,I,J,L)*y(n_Alkenes,L)
     & *y(nO3,L)*0.65+rr(31,I,J,L)*y(n_Isoprene,L)*y(nO3,L)*0.58
c
       sqroot=sqrt(bqqz*bqqz+4.*aqqz*cqqz)
       y(nHO2,L)=(sqroot-bqqz)/(2.*aqqz)
       y(nOH,L)=pHOx(I,J,L)*y(nHO2,L)
       temp_yHOx=y(nOH,L)+y(nHO2,L)
c
c      Now partition HOx into OH and HO2:
c      CZ: OH->HO2 reactions :
       cz=rr(2,I,J,L)*y(nO3,L)+rr(13,I,J,L)*y(n_CO,L)
     & +rr(14,I,J,L)*y(n_H2O2,L)+rr(19,I,J,L)*y(nH2,L)
     & +rr(21,I,J,L)*y(n_HCHO,L)+rr(37,I,J,L)*y(n_Paraffin,L)*
     & *0.11+rr(30,I,J,L)*y(n_Isoprene,L)*0.85
     & +rr(34,I,J,L)*y(n_Alkenes,L)
C
C      DZ: HO2->OH reactions :
       dz=rr(4,I,J,L)*y(nO3,L)+rr(6,I,J,L)*y(nNO,L)
     & +rr(41,I,J,L)*0.79*y(nC2O3,L)+rr(15,I,J,L)*y(nHO2,L)
     & *(2.*ss(4,I,J,L)/(2.*ss(4,I,J,L)+y(nOH,L)*rr(14,I,J,L)))
c      Previous two lines additional OH production via rxn 41,
c      which also produces HO2, and R15 then S4/(S4+S14) fraction
C
       if(cz+dz.gt.0.)then
        y(nOH,L)=(dz/(cz+dz))*temp_yHOx
        if(y(nOH,L).gt.temp_yHOx)y(nOH,L)=temp_yHOx-1.0
c---->  warning: OH caps follow
        if(j.le.3.and.y(nOH,L).ge.3.E5) y(nOH,L)=3.E5
        if(j.ge.44.and.y(nOH,L).ge.3.E5)y(nOH,L)=3.E5
       else
        y(nOH,L)=1.0
       endif
       y(nHO2,L)=(temp_yHOx-y(nOH,L))
C      Some limits on OH, HO2:
       if(y(nOH,L).lt.1.)y(nOH,L)=1.0
       if(y(nHO2,L).lt.1.)y(nHO2,L)=1.0
       if(y(nHO2,L).gt.1.E9)y(nHO2,L)=1.E9
       pHOx(I,J,L)=y(nOH,L)/y(nHO2,L)
C
      enddo  ! >> end of altitude loop <<

      return
      END SUBROUTINE HOxfam

