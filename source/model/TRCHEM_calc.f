      SUBROUTINE chemstep(I,J,change)
!@sum chemstep Calculate new concentrations after photolysis & chemistry
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on ds3ch4_chem_calc_jun1202_M23)
!@calls rates,chem1,chem1prn
c
C**** GLOBAL parameters and variables:
C
      USE RESOLUTION, only    : im,jm,lm,ls1
      USE DYNAMICS, only      : am, byam
      USE GEOM, only          : BYDXYP,dxyp
      USE TRACER_COM, only: n_CH4,n_CH3OOH,n_Paraffin,n_PAN,n_Isoprene,
     &                   n_AlkylNit,n_Alkenes,n_N2O5,n_NOx,n_HO2NO2,
     &                   n_Ox,n_HNO3,n_H2O2,n_CO,n_HCHO,trm,ntm
      USE TRCHEM_Shindell_COM, only: chemrate,photrate,mass2vol,
     &                   yCH3O2,yC2O3,yXO2,yXO2N,yRXPAR,yAldehyde,
     &                   yROR,nCH3O2,nC2O3,nXO2,nXO2N,nRXPAR,
     &                   nAldehyde,nROR,nr,nn,dt2,nss,ks,dest,prod,
     &                   ny,nhet,rr,nO1D,nOH,nNO,nHO2,ta,nM,ss,
     &                   nO3,nNO2,nNO3,prnrts,jprn,iprn,lprn,ay,
     &                   prnchg,y,bymass2vol,kss,nps,kps,nds,kds,
     &                   npnr,nnr,ndnr,kpnr,kdnr,nH2O,changeAldehyde
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var change change due to chemistry in mass/time
!@var I,J passed horizontal spatial indicies
!@var L,iter dummy loop variable
!@var maxl highest level with chemistry
!@var qqqCH3O2,CH3O2loss,XO2_NO,XO2N_HO2,RXPAR_PAR,ROR_CH2,C2O3prod,
!@+   C2O3dest,XO2prod,XO2dest,XO2_XO2,XO2Nprod,XO2Ndest,RXPARprod,
!@+   RXPARdest,Aldehydeprod,Aldehydedest,RORprod,RORdest,total,
!@+   rnewval,sumN,dNOx,ratio,sumD,newD,ratioD,newP,ratioP,changeA
!@+   sumP dummy temp variables
!@var tempiter,tempiter2 temp vars for equilibrium calcs iterations
      REAL*8, DIMENSION(IM,JM,LM,ntm)     :: change
      INTEGER, INTENT(IN) :: I,J
      INTEGER L,iter,maxl,igas
      REAL*8 qqqCH3O2,CH3O2loss,XO2_NO,XO2N_HO2,RXPAR_PAR,ROR_CH2,
     & C2O3prod,C2O3dest,XO2prod,XO2dest,XO2_XO2,XO2Nprod,XO2Ndest,
     & RXPARprod,RXPARdest,Aldehydeprod,Aldehydedest,RORprod,RORdest,
     & total,rnewval,sumN,dNOx,ratio,sumD,newD,ratioD,newP,ratioP,
     & changeA,sumP,tempiter,tempiter2
C
C     TROPOSPHERIC CHEMISTRY ONLY:
      maxl=LS1-1

      do L=1,maxl
       y(nCH3O2,L)   =yCH3O2(I,J,L)
       y(nC2O3,L)    =yC2O3(I,J,L)
       y(nXO2,L)     =yXO2(I,J,L)
       y(nXO2N,L)    =yXO2N(I,J,L)
       y(nRXPAR,L)   =yRXPAR(I,J,L)
       y(nAldehyde,L)=yAldehyde(I,J,L)
       y(nROR,L)     =yROR(I,J,L)
      enddo
C
cc    calculate reaction rates with present concentrations
      call rates(maxl,I,J)
c
c     chem1 call sample:
c     (klist,l,numel,nlist,ndlist,rate,change,multip)
c              numel=number of elements in reaction list nlist (1 or 2)
c              change=dest or prod array
c              multip=1(prod) or -1(dest)
c     chemical destruction:
      call chem1(kdnr,maxl,2,nn,ndnr,chemrate,dest,-1)
c     chemical production:
      call chem1(kpnr,maxl,2,nnr,npnr,chemrate,prod,1)
c     photolytic destruction:
      call chem1(kds,maxl,1,ks,nds,photrate,dest,-1)
c     photolytic production:
      call chem1(kps,maxl,2,kss,nps,photrate,prod,1)
C
c     Oxidation of Isoprene and Alkenes produces less than one
c      HCHO, Alkenes, and CO per rxn, correct here following Houweling
      do L=1,maxl
       prod(n_CO,L)=prod(n_CO,L)-0.63*rr(35,I,J,L)*y(n_Alkenes,L)*
     &  y(nO3,L)*dt2
       prod(n_HCHO,L)=prod(n_HCHO,L)-0.36*rr(35,I,J,L)*y(n_Alkenes,L)*
     &  y(nO3,L)*dt2
       prod(n_HCHO,L)=prod(n_HCHO,L)-0.39*rr(30,I,J,L)*y(n_Isoprene,L)*
     &  y(nOH,L)*dt2
       prod(n_Alkenes,L)=prod(n_Alkenes,L)-0.42*rr(30,I,J,L)*
     &  y(n_Isoprene,L)*y(nOH,L)*dt2
       prod(n_HCHO,L)=prod(n_HCHO,L)-0.10*rr(31,I,J,L)*y(n_Isoprene,L)*
     &  y(nO3,L)*dt2
       prod(n_Alkenes,L)=prod(n_Alkenes,L)-0.45*rr(31,I,J,L)*
     &  y(n_Isoprene,L)*y(nO3,L)*dt2
      enddo
c
c     set CH3O2 values (concentration = production/specific loss)
      do L=1,maxl
        iter=1
        qqqCH3O2=(rr(11,I,J,L)*y(nO1D,L)+rr(12,I,J,L)*y(nOH,L))
     &  *y(n_CH4,L)+rr(23,I,J,L)*y(n_CH3OOH,L)*y(nOH,L)

       tempiter=rr(20,I,J,L)*y(nNO,L)+rr(22,I,J,L)*y(nHO2,L)
 10    CH3O2loss=tempiter+rr(27,I,J,L)*yCH3O2(I,J,L)
       if(CH3O2loss.gt.1E-7)then
          y(nCH3O2,L)=qqqCH3O2/CH3O2loss
        else
          y(nCH3O2,L)=1.0
        endif
        yCH3O2(I,J,L)=y(nCH3O2,L)
        iter=iter+1
        if(iter.le.7)goto10 ! replace with while loop?
      enddo
c
      do L=1,maxl
c
c       Set C2O3, XO2, XO2N, RXPAR, Aldehyde & ROR values
c
c        First set various specific loss rates
         XO2_NO=y(nNO,L)*4.2E-12*exp(180./ta(L))
         XO2N_HO2=y(nHO2,L)*y(nNO,L)*rr(44,I,J,L)*
     &   rr(43,I,J,L)/XO2_NO
         RXPAR_PAR=y(n_Paraffin,L)*8.E-11
         ROR_CH2=1.6E3
c
c       Set value for C2O3
        iter=1
        C2O3prod=rr(38,I,J,L)*yAldehyde(I,J,L)*y(nOH,L)+
     &  (rr(29,I,J,L)*y(nM,L)+ss(15,I,J,L))*y(n_PAN,L)+0.15*
     &  rr(31,I,J,L)*y(nO3,L)*y(n_Isoprene,L)
        tempiter=rr(39,I,J,L)*y(nNO,L)+rr(54,I,J,L)*y(nNO2,L)+
     &   rr(41,I,J,L)*y(nHO2,L)
 20     C2O3dest=tempiter+rr(40,I,J,L)*yC2O3(I,J,L)
        if(C2O3dest.gt.1E-7)then
          y(nC2O3,L)=(C2O3prod/C2O3dest)
        else
          y(nC2O3,L)=1.0
        endif
        yC2O3(I,J,L)=y(nC2O3,L)
        iter=iter+1
        if(iter.le.7)goto20     ! replace with while loop?
c
c       Set value for XO2
        iter=1
        XO2prod=ss(16,I,J,L)*yAldehyde(I,J,L)+
     &  y(nC2O3,L)*(rr(39,I,J,L)*y(nNO2,L)+rr(40,I,J,L)*
     &  y(nC2O3,L)*2.+rr(41,I,J,L)*y(nHO2,L))
     &  +rr(42,I,J,L)*yROR(I,J,L)*0.96
     &  +y(nOH,L)*(rr(37,I,J,L)*y(n_Paraffin,L)*0.87+rr(34,I,J,L)*
     &  y(n_Alkenes,L)+rr(30,I,J,L)*y(n_Isoprene,L)*0.85+
     &  rr(33,I,J,L)*y(n_AlkylNit,L))+
     &  y(nO3,L)*(rr(35,I,J,L)*y(n_Alkenes,L)*0.29+
     &  rr(31,I,J,L)*y(n_Isoprene,L)*0.18)
        tempiter=XO2_NO+rr(43,I,J,L)*y(nHO2,L)
        tempiter2=1.7E-14*exp(1300./ta(L))
 30     XO2_XO2=yXO2(I,J,L)*tempiter2
        XO2dest=tempiter+XO2_XO2
        if(XO2dest.gt.1E-7.and.ss(16,I,J,L).gt.5E-5)then
          y(nXO2,L)=(XO2prod/XO2dest)
        else
          y(nXO2,L)=1.0
        endif
        yXO2(I,J,L)=y(nXO2,L)
        iter=iter+1
        if(iter.le.7)goto30   ! replace with while loop?
c
c       Set value for XO2N
        XO2Nprod=rr(37,I,J,L)*y(n_Paraffin,L)*y(nOH,L)*0.13+
     &  rr(42,I,J,L)*yROR(I,J,L)*0.04+rr(30,I,J,L)*y(n_Isoprene,L)*
     &  y(nOH,L)*0.15
        XO2Ndest=XO2N_HO2+rr(44,I,J,L)*y(nNO,L)
        if(XO2Ndest.gt.1E-7)then
          y(nXO2N,L)=(XO2Nprod/XO2Ndest)
        else
          y(nXO2N,L)=1.0
        endif
        yXO2N(I,J,L)=y(nXO2N,L)
c
c       Set value for RXPAR
        RXPARprod=rr(37,I,J,L)*y(n_Paraffin,L)*y(nOH,L)*0.11+
     &  rr(34,I,J,L)*yROR(I,J,L)*2.1+rr(35,I,J,L)*y(n_Alkenes,L)*
     &  y(nO3,L)*0.9
        RXPARdest=RXPAR_PAR
        if(RXPARdest.gt.0.)then
          y(nRXPAR,L)=(RXPARprod/RXPARdest)
        else
          y(nRXPAR,L)=1.0
        endif
        yRXPAR(I,J,L)=y(nRXPAR,L)
c
c       Set value for Aldehyde
        Aldehydeprod=rr(37,I,J,L)*y(n_Paraffin,L)*y(nOH,L)*0.11+
     &  rr(34,I,J,L)*y(n_Alkenes,L)*y(nOH,L)+
     &  rr(42,I,J,L)*yROR(I,J,L)*1.1+rr(35,I,J,L)*y(n_Alkenes,L)*
     &  y(nO3,L)*0.44
        Aldehydedest=rr(38,I,J,L)*y(nOH,L)+ss(16,I,J,L)

c       Check for equilibrium
        if(Aldehydedest*y(nAldehyde,L)*dt2.lt.y(nAldehyde,L))then
         changeAldehyde=
     &   (Aldehydeprod-y(nAldehyde,L)*Aldehydedest)*dt2
         if(changeAldehyde.gt.y(nAldehyde,L))
     &   changeAldehyde=y(nAldehyde,L)
         y(nAldehyde,L)=y(nAldehyde,L)+changeAldehyde
         if(y(nAldehyde,L).lt.0.)y(nAldehyde,L)=1.0
        else
          y(nAldehyde,L)=(Aldehydeprod/(Aldehydedest+0.5E-5))
        endif
        yAldehyde(I,J,L)=y(nAldehyde,L)
c
c       Set value for ROR
        RORprod=rr(37,I,J,L)*y(n_Paraffin,L)*y(nOH,L)*0.76+
     &  rr(42,I,J,L)*yROR(I,J,L)*0.02
        RORdest=rr(42,I,J,L)+ROR_CH2
        if(RORdest.gt.0.)then
          y(nROR,L)=(RORprod/RORdest)
        else
          y(nROR,L)=1.0
        endif
        yROR(I,J,L)=y(nROR,L)
c
c       Add in parrafin loss term via rxpar reaction
        dest(n_Paraffin,L)=dest(n_Paraffin,L)-y(nRXPAR,L)*RXPAR_PAR*dt2

c       Add in CH3OOH production via XO2N + HO2
        prod(n_CH3OOH,L)=prod(n_CH3OOH,L)+XO2N_HO2*y(nXO2N,L)*dt2
c
      end do  ! L
c
c     If NOx in equil with N2O5, HO2NO2 or PAN, remove from changes
      do L=1,maxl
       if(-dest(n_N2O5,L).ge.y(n_N2O5,L))then
        dest(n_NOx,L)=dest(n_NOx,L)+2.*chemrate(52,L)
        prod(n_NOx,L)=prod(n_NOx,L)-2.*(chemrate(46,L)+photrate(7,L))
       endif
       if(-dest(n_HO2NO2,L).ge.y(n_HO2NO2,L))then
        dest(n_NOx,L)=dest(n_NOx,L)+chemrate(51,L)
        prod(n_NOx,L)=prod(n_NOx,L)-(chemrate(18,L)+chemrate(45,L)+
     &   photrate(10,L)+photrate(11,L))
       endif
       if(-dest(n_PAN,L).ge.y(n_PAN,L))then
        dest(n_NOx,L)=dest(n_NOx,L)+chemrate(54,L)
        prod(n_NOx,L)=prod(n_NOx,L)-(chemrate(29,L)+photrate(15,L))
       endif
      enddo
c
c     print rxtn rates if desired (chem1prn : argument before multip is
c     index = number of call):
c
       if(prnrts.and.J.eq.jprn.and.I.eq.iprn)then
        do igas=1,ntm
         total=0.
         write(6,108) ' Species: ',ay(1,igas),ay(2,igas)
c
         call chem1prn
     *   (kdnr,2,nn,ndnr,chemrate,1,-1,igas,total,maxl,I,J)
c
         if(igas.eq.n_Ox)write(6,109)
     *   ' losses from 10 and 11 reset by fraction of O1D via rxn,',
     *   'dy = ?'
         if(igas.eq.n_NOx.and.-dest(n_HO2NO2,lprn).ge.y(n_HO2NO2,lprn))
     *    write(6,110)'loss by reaction 51 (HO2NO2 formation) removed',
     *    rr(51,I,J,lprn)*y(nHO2,lprn)*y(nNO2,lprn)*dt2
         if(igas.eq.n_NOx.and.-dest(n_N2O5,lprn).ge.y(n_N2O5,lprn))
     *    write(6,110)'losses by reaction 52 (N2O5 formation) removed',
     *    2*rr(52,I,J,lprn)*y(nNO3,lprn)*y(nNO2,lprn)*dt2
         if(igas.eq.n_NOx.and.-dest(n_PAN,lprn).ge.y(n_PAN,lprn))
     *    write(6,110)'losses by reaction 54 (PAN formation) removed',
     *    chemrate(54,lprn)
c
         call chem1prn
     &   (kpnr,2,nnr,npnr,chemrate,2,1,igas,total,maxl,I,J)
c
         if(igas.eq.n_NOx.and.-dest(n_HO2NO2,lprn).ge.y(n_HO2NO2,lprn))
     *    write(6,110)'gain by reactions 18 & 45 (from HO2NO2) rmoved',
     *    (rr(18,I,J,lprn)*y(nOH,L)+rr(45,I,J,lprn)*y(nM,lprn)+
     *    ss(10,I,J,lprn)+ss(11,I,J,lprn))*y(n_HO2NO2,lprn)*dt2
         if(igas.eq.n_NOx.and.-dest(n_N2O5,lprn).ge.y(n_N2O5,lprn))
     *        write(6,110)
     *        'gains by reaction 46 (N2O5 decomposition) removed',
     *        2.*(rr(46,I,J,lprn)*y(nM,lprn))*y(n_N2O5,lprn)*dt2
         if(igas.eq.n_NOx.and.-dest(n_PAN,lprn).ge.y(n_PAN,lprn))
     *    write(6,110)'gain by reaction 29 (from PAN) removed',
     *    rr(29,I,J,lprn)*y(nM,lprn)*y(n_PAN,lprn)*dt2
c
         call chem1prn(kds,1,ks,nds,photrate,3,-1,igas,total,maxl,I,J)
c
         call chem1prn(kps,2,kss,nps,photrate,4,1,igas,total,maxl,I,J)
c
         if(igas.eq.n_NOx.and.-dest(n_N2O5,lprn).ge.y(n_N2O5,lprn))
     *    write(6,110)'gains by reaction 7 (N2O5 photolysis) removed',
     *    ss(7,I,J,lprn)*y(n_N2O5,lprn)*dt2
         if(igas.eq.n_NOx.and.-dest(n_N2O5,lprn).ge.y(n_N2O5,lprn))
     *        write(6,110)
     *        'net change due to N2O5 is ',
     *        2.*(y(n_N2O5,lprn)-(rr(52,I,J,lprn)*y(nNO3,lprn)*
     *        y(nNO2,lprn))/
     *        (rr(46,I,J,lprn)*y(nM,lprn)+ss(7,I,J,lprn)))
         if(igas.eq.n_NOx.and.-dest(n_HO2NO2,lprn).ge.y(n_HO2NO2,lprn))
     *    write(6,110)'gain by rxns 10 & 11 (HO2NO2 photolysis) rmoved'
     *    ,(ss(10,I,J,lprn)+ss(11,I,J,lprn))*y(n_HO2NO2,lprn)*dt2
         if(igas.eq.n_NOx.and.-dest(n_HO2NO2,lprn).ge.y(n_HO2NO2,lprn))
     *        write(6,110)'net change due to HO2NO2 is ',
     *        y(n_HO2NO2,lprn)-((rr(51,I,J,lprn)*y(nHO2,lprn)*
     *       y(nNO2,lprn))/(rr(18,I,J,lprn)*
     *     y(nOH,lprn)+rr(45,I,J,lprn)*y(nM,lprn)+ss(10,I,J,lprn)
     *     +ss(11,I,J,lprn)))
         if(igas.eq.n_NOx.and.-dest(n_PAN,lprn).ge.y(n_PAN,lprn))
     *        write(6,110)'net change due to PAN is ',
     *        y(n_PAN,lprn)-((rr(54,I,J,lprn)*y(nC2O3,lprn)*
     *        y(nNO2,lprn))/
     *        (rr(29,I,J,lprn)*y(nM,lprn)+ss(15,I,J,lprn)))
         if(igas.eq.n_Ox.or.igas.eq.n_NOx)total=
     *    100.*(dest(igas,lprn)+prod(igas,lprn))/y(igas,lprn)
         if(igas.eq.n_CH3OOH) write(6,'(a48,a6,e10.3)')
     *    'production from XO2N + HO2 ','dy = ',
     *    y(nHO2,lprn)*y(nNO,lprn)*rr(44,I,J,lprn)*rr(43,I,J,lprn)/
     *    (y(nNO,lprn)*4.2E-12*exp(180./ta(lprn)))
     *    *y(nXO2N,lprn)*dt2
         if(igas.eq.n_Paraffin) write(6,'(a48,a6,e10.3)')
     *    'destruction from RXPAR ',
     *    'dy = ',-y(nRXPAR,lprn)*y(n_Paraffin,lprn)*8.E-11*dt2
         write(6,118) ' Total change in ',ay(1,igas),ay(2,igas),
     *  ' is ',total,' percent; dy= ',dest(igas,lprn)+prod(igas,lprn)
         write(6,*)
        enddo ! igas
       endif  ! chem diags
 108  format(a10,2x,2a4)
 109  format(a55,7x,a6,e10.3)
 110  format(a68,e10.3)
 118  format(a17,2a4,a4,f10.0,a14,e12.3)
 119  format(a45,7x,f10.0,a5,e12.3)
c
      if(prnchg.and.J.eq.jprn.and.I.eq.iprn)write(6,'(a35,3(2x,i2))')
     * ' Total change by species at I, J, L',i,j,lprn
c
cc    Loop to calculate tracer changes
      do igas=1,ntm
       do L=1,maxl
c
         change(I,J,L,igas)=
     &   (dest(igas,L)+prod(igas,L))*dxyp(J)*AM(L,I,J)*
     &   bymass2vol(igas)/y(nM,L)
c
c        Set N2O5 to equilibrium when necessary (near ground,
c        N2O5 is thermally unstable, has a very short lifetime)
         if(igas.eq.n_N2O5.and.-dest(igas,L).ge.y(n_N2O5,L))then
           rnewval=(rr(52,I,J,L)*y(nNO3,L)*y(nNO2,L))/
     &     (rr(46,I,J,L)*y(nM,L)+ss(7,I,J,L))
           if(rnewval.lt.1.)rnewval=1.
           change(I,J,L,igas)=(rnewval-y(n_N2O5,L))*dxyp(J)*AM(L,I,J)*
     &     bymass2vol(igas)/y(nM,L)
c          endif
         endif
c
c        Conserve NOx with respect to N2O5
         if(igas.eq.n_NOx.and.-dest(n_N2O5,L).ge.y(n_N2O5,L))then
          rnewval=(rr(52,I,J,L)*y(nNO3,L)*y(nNO2,L))/
     &    (rr(46,I,J,L)*y(nM,L)+ss(7,I,J,L))
          change(I,J,L,igas)=
     &    change(I,J,L,igas)+2.*(y(n_N2O5,L)-rnewval)
     &    *dxyp(J)*AM(L,I,J)*bymass2vol(igas)/y(nM,L)
         endif
c
c        Set HO2NO2 to equil when necessary
         if(igas.eq.n_HO2NO2.and.-dest(igas,L).ge.y(n_HO2NO2,L))then
          rnewval=(rr(51,I,J,L)*y(nHO2,L)*y(nNO2,L))/
     *    (rr(18,I,J,L)*y(nOH,L)+rr(45,I,J,L)*y(nM,L)+ss(10,I,J,L)
     *    +ss(11,I,J,L))
          if(rnewval.lt.1.)rnewval=1.
          change(I,J,L,igas)=(rnewval-y(n_HO2NO2,L))*dxyp(J)*AM(L,I,J)*
     &    bymass2vol(igas)/y(nM,L)
         endif
c
c        Conserve NOx with respect to HO2NO2
         if(igas.eq.n_NOx.and.-dest(n_HO2NO2,L).ge.y(n_HO2NO2,L))then
          rnewval=(rr(51,I,J,L)*y(nHO2,L)*y(nNO2,L))/
     *    (rr(18,I,J,L)*y(nOH,L)+rr(45,I,J,L)*y(nM,L)+ss(10,I,J,L)
     *    +ss(11,I,J,L))
          change(I,J,L,igas)=
     &    change(I,J,L,igas)+(y(n_HO2NO2,L)-rnewval)*
     &    dxyp(J)*AM(L,I,J)*bymass2vol(igas)/y(nM,L)
         endif
c
c        Set PAN to equilibrium when necessary (near ground,
c         PAN is thermally unstable, has a very short lifetime)
         if(igas.eq.n_PAN.and.-dest(igas,L).ge.y(n_PAN,L))then
          rnewval=(rr(54,I,J,L)*y(nC2O3,L)*y(nNO2,L))/
     *    (rr(29,I,J,L)*y(nM,L)+ss(15,I,J,L))
          if(rnewval.lt.1.)rnewval=1.
          change(I,J,L,igas)=(rnewval-y(n_PAN,L))*dxyp(J)*AM(L,I,J)
     &    *bymass2vol(igas)/y(nM,L)
         endif
c
c        Conserve NOx with respect to PAN
         if(igas.eq.n_NOx.and.-dest(n_PAN,L).ge.y(n_PAN,L))then
          rnewval=(rr(54,I,J,L)*y(nC2O3,L)*y(nNO2,L))/
     *    (rr(29,I,J,L)*y(nM,L)+ss(15,I,J,L))
          if(rnewval.lt.1.)rnewval=1.
          change(I,J,L,igas)=change(I,J,L,igas)+(y(n_PAN,L)-rnewval)*
     &    dxyp(J)*AM(L,I,J)*bymass2vol(igas)/y(nM,L)
         endif
c
       end do ! L
      end do  ! igas
c
cc    Insure nitrogen conservation
      if(prnchg.and.J.eq.jprn.and.I.eq.iprn)then
       write(*,*) 'changes before nitrogen conservation routine'
       write(*,*) 'NOx, N2O5, HO2NO2, HNO3, PAN, AlkylNit'
       write(*,*) change(I,J,lprn,n_NOx),change(I,J,lprn,n_N2O5),
     & change(I,J,lprn,n_HO2NO2),change(I,J,lprn,n_HNO3),
     & change(I,J,lprn,n_PAN),change(I,J,lprn,n_AlkylNit)
      endif
c
      do L=1,maxl
c
c       First check for nitrogen loss > 100%
        if(-change(I,J,L,n_NOx).gt.trm(I,J,L,n_NOx))
     &  change(I,J,L,n_NOx)=1.-trm(I,J,L,n_NOx)
        if(-change(I,J,L,n_N2O5).gt.trm(I,J,L,n_N2O5))
     &  change(I,J,L,n_N2O5)=1.-trm(I,J,L,n_N2O5)
        if(-change(I,J,L,n_HO2NO2).gt.trm(I,J,L,n_HO2NO2))
     &  change(I,J,L,n_HO2NO2)=1.-trm(I,J,L,n_HO2NO2)
        if(-change(I,J,L,n_HNO3).gt.trm(I,J,L,n_HNO3))
     &  change(I,J,L,n_HNO3)=1.-trm(I,J,L,n_HNO3)
        if(-change(I,J,L,n_PAN).gt.trm(I,J,L,n_PAN))
     &  change(I,J,L,n_PAN)=1.-trm(I,J,L,n_PAN)
        if(-change(I,J,L,n_AlkylNit).gt.trm(I,J,L,n_AlkylNit))
     &  change(I,J,L,n_AlkylNit)=1.-trm(I,J,L,n_AlkylNit)
c
c       Next insure balance between dNOx and sum of dOthers
        sumN=(2.*change(I,J,L,n_N2O5))*mass2vol(n_N2O5)+
     *  (change(I,J,L,n_HNO3))*mass2vol(n_HNO3)+
     *  (change(I,J,L,n_HO2NO2))*mass2vol(n_HO2NO2)+
     *  (change(I,J,L,n_PAN))*mass2vol(n_PAN)+
     *  (change(I,J,L,n_AlkylNit))*mass2vol(n_AlkylNit)
        dNOx=change(I,J,L,n_NOx)*mass2vol(n_NOx)
        if(prnchg.and.J.eq.jprn.and.I.eq.iprn.and.L.eq.lprn)
     &  write(*,*) 'sumN, dNOx = ',sumN,dNOx
        ratio=-sumN/dNOx
        IF(ratio.le.0.999 .OR. ratio.ge.1.001) THEN
         if(dNOx.gt.0.)then
c         NOx being produced (net positive change)
          if (ratio.gt.1.)then
           sumD=0.
c          reduce N destruction to match NOx prodcution
           if(change(I,J,L,n_N2O5).lt.0.)   sumD=sumD+
     *     2.*change(I,J,L,n_N2O5)*mass2vol(n_N2O5)
           if(change(I,J,L,n_HO2NO2).lt.0.) sumD=sumD+
     *     change(I,J,L,n_HO2NO2)*mass2vol(n_HO2NO2)
           if(change(I,J,L,n_HNO3).lt.0.)   sumD=sumD+
     *     change(I,J,L,n_HNO3)*mass2vol(n_HNO3)
           if(change(I,J,L,n_PAN).lt.0.)    sumD=sumD+
     *     change(I,J,L,n_PAN)*mass2vol(n_PAN)
           if(change(I,J,L,n_AlkylNit).lt.0.)sumD=sumD+
     *     change(I,J,L,n_AlkylNit)*mass2vol(n_AlkylNit)
           newD=(sumN/ratio)+sumD-sumN
           ratioD=newD/sumD
           if(change(I,J,L,n_N2O5).lt.0.)    change(I,J,L,n_N2O5)=
     *     change(I,J,L,n_N2O5)    *ratioD
           if(change(I,J,L,n_HO2NO2).lt.0.)  change(I,J,L,n_HO2NO2)=
     *     change(I,J,L,n_HO2NO2)  *ratioD
           if(change(I,J,L,n_HNO3).lt.0.)    change(I,J,L,n_HNO3)=
     *     change(I,J,L,n_HNO3)    *ratioD
           if(change(I,J,L,n_PAN).lt.0.)     change(I,J,L,n_PAN)=
     *     change(I,J,L,n_PAN)     *ratioD
           if(change(I,J,L,n_AlkylNit).lt.0.)change(I,J,L,n_AlkylNit)=
     *     change(I,J,L,n_AlkylNit)*ratioD
          endif
c
          if (ratio.le.1..and.ratio.gt.0.)then
c          reduce NOx production to match N loss
           change(I,J,L,n_NOx)=change(I,J,L,n_NOx)*ratio
          endif
         else
c         NOx being destroyed (net change is negative)
          if (ratio.gt.1.)then
           sumP=0
c          reduce N production to match NOx loss
           if(change(I,J,L,n_N2O5).gt.0.)    sumP=sumP+
     *     2.*change(I,J,L,n_N2O5)*mass2vol(n_N2O5)
           if(change(I,J,L,n_HO2NO2).gt.0.)  sumP=sumP+
     *     change(I,J,L,n_HO2NO2)*mass2vol(n_HO2NO2)
           if(change(I,J,L,n_HNO3).gt.0.)    sumP=sumP+
     *     change(I,J,L,n_HNO3)*mass2vol(n_HNO3)
           if(change(I,J,L,n_PAN).gt.0.)     sumP=sumP+
     *     change(I,J,L,n_PAN)*mass2vol(n_PAN)
           if(change(I,J,L,n_AlkylNit).gt.0.)sumP=sumP+
     *     change(I,J,L,n_AlkylNit)*mass2vol(n_AlkylNit)
           newP=(sumN/ratio)+sumP-sumN
           ratioP=newP/sumP
           if(change(I,J,L,n_N2O5).gt.0)     change(I,J,L,n_N2O5)=
     *        change(I,J,L,n_N2O5)*ratioP
           if(change(I,J,L,n_HO2NO2).gt.0.)  change(I,J,L,n_HO2NO2)=
     *        change(I,J,L,n_HO2NO2)*ratioP
           if(change(I,J,L,n_HNO3).gt.0.)    change(I,J,L,n_HNO3)=
     *        change(I,J,L,n_HNO3)*ratioP
           if(change(I,J,L,n_PAN).gt.0.)     change(I,J,L,n_PAN)=
     *        change(I,J,L,n_PAN)*ratioP
           if(change(I,J,L,n_AlkylNit).gt.0.)change(I,J,L,n_AlkylNit)=
     *        change(I,J,L,n_AlkylNit)*ratioP
          endif
          if (ratio.le.1.and.ratio.gt.0)then
c          increase N production to match NOx loss (12/4/2001)
           sumP=0.
           if(change(I,J,L,n_N2O5).gt.0.)    sumP=sumP+
     *     2.*change(I,J,L,n_N2O5)*mass2vol(n_N2O5)
           if(change(I,J,L,n_HO2NO2).gt.0.)  sumP=sumP+
     *     change(I,J,L,n_HO2NO2)*mass2vol(n_HO2NO2)
           if(change(I,J,L,n_HNO3).gt.0.)    sumP=sumP+
     *     change(I,J,L,n_HNO3)*mass2vol(n_HNO3)
           if(change(I,J,L,n_PAN).gt.0.)     sumP=sumP+
     *     change(I,J,L,n_PAN)*mass2vol(n_PAN)
           if(change(I,J,L,n_AlkylNit).gt.0.)sumP=sumP+
     *     change(I,J,L,n_AlkylNit)*mass2vol(n_AlkylNit)
           newP=(sumN/ratio)+sumP-sumN
           ratioP=newP/sumP
           if(change(I,J,L,n_N2O5).gt.0.)    change(I,J,L,n_N2O5)=
     *     change(I,J,L,n_N2O5)*ratioP
           if(change(I,J,L,n_HO2NO2).gt.0.)  change(I,J,L,n_HO2NO2)=
     *     change(I,J,L,n_HO2NO2)*ratioP
           if(change(I,J,L,n_HNO3).gt.0.)    change(I,J,L,n_HNO3)=
     *     change(I,J,L,n_HNO3)*ratioP
           if(change(I,J,L,n_PAN).gt.0.)     change(I,J,L,n_PAN)=
     *     change(I,J,L,n_PAN)*ratioP
           if(change(I,J,L,n_AlkylNit).gt.0.)change(I,J,L,n_AlkylNit)=
     *     change(I,J,L,n_AlkylNit)*ratioP
          endif
         endif

        END IF ! skipped section above if ratio very close to one
C
        if(prnchg.and.J.eq.jprn.and.I.eq.iprn.and.L.eq.lprn)
     &  write(*,*) 'ratio for conservation =',ratio
C
      end do ! L
C
cc    Print chemical changes in a particular grid box if desired
      if(prnchg.and.J.eq.jprn.and.I.eq.iprn)then
        do igas=1,ntm
         changeA=change(I,J,Lprn,igas)*y(nM,lprn)*mass2vol(igas)*
     *   bydxyp(J)*byam(lprn,I,J)
         if(y(igas,lprn).eq.0.)then
            write(6,156) ay(1,igas),ay(2,igas),': ',changeA
     *            ,' molecules;  y=0'
         else
            write(6,155) ay(1,igas),ay(2,igas),': ',changeA
     *           ,' molecules produced; ',
     *   (100.*changeA)/y(igas,lprn),' percent of'
     *   ,y(igas,lprn),'(',1.E9*y(igas,lprn)/y(nM,lprn),' ppbv)'
         endif
c
         if(igas.eq.ntm)then
          write(6,'(a10,58x,e13.3,6x,f10.3,a5)')
     *    ' H2O     :',y(nH2O,LPRN),(y(nH2O,LPRN)/
     *    y(nM,LPRN))*1.E6,' ppmv'
          write(6,'(a10,58x,e13.3,6x,f10.3,a5)')
     *    ' CH3O2   :',yCH3O2(I,J,LPRN),(yCH3O2(I,J,LPRN)/
     *    y(nM,LPRN))*1.E9,' ppbv'
          write(6,'(a10,58x,e13.3,6x,f10.3,a5)')
     *    ' C2O3    :',y(nC2O3,LPRN),(y(nC2O3,LPRN)/y(nM,LPRN))*1.E9,
     *    ' ppbv'
          write(6,'(a10,58x,e13.3,6x,f10.3,a5)')
     *    ' XO2     :',y(nXO2,LPRN),(y(nXO2,LPRN)/y(nM,LPRN))*1.E9,
     *    ' ppbv'
          write(6,'(a10,58x,e13.3,6x,f10.3,a5)')
     *    ' XO2N    :',y(nXO2N,LPRN),(y(nXO2N,LPRN)/
     *    y(nM,LPRN))*1.E9,' ppbv'
          write(6,'(a10,58x,e13.3,6x,f10.3,a5)')
     *    ' RXPAR   :',y(nRXPAR,LPRN),(y(nRXPAR,LPRN)/
     *    y(nM,LPRN))*1.E9,' ppbv'
          write(6,'(a10,58x,e13.3,6x,f10.3,a5)')
     *    ' Aldehyde:',y(nAldehyde,LPRN),(y(nAldehyde,LPRN)/
     *    y(nM,LPRN))*1.E9,' ppbv'
          write(6,'(a10,58x,e13.3,6x,f10.3,a5)')
     *    ' ROR     :',y(nROR,LPRN),(y(nROR,LPRN)/
     *    y(nM,LPRN))*1.E9,' ppbv'
         endif
        enddo
      endif  !end of prnchg loop
c
c*** tracer masses & slopes are now updated in apply_tracer_3Dsource ***
      do igas=1,ntm
       do L=1,maxl
        if(change(I,J,L,igas).gt.1.E20 .OR.
     &    -change(I,J,L,igas).gt.trm(I,J,L,igas)) THEN
           WRITE(6,*)'>>change set to 0 in chemstep: I,J,L,igas,change'
     &    ,I,J,L,igas,change(I,J,L,igas)
           change(I,J,L,igas) = 0.
        endif
c
c Leaving this here to reinstate it later in some form. GSF 2/14/02 :
c         if(change(I,J,L,igas).le.1.E20) THEN
c           if(L.eq.1.and.igas.eq.1)then
c             changeA=(change(I,J,1,igas)*y(nM,1)*mass2vol(igas))*
c     *       bydxyp(J)*byam(1,I,J)
c             TAIJS(I,J,1)=TAIJS(I,J,1)+
c     *       1.E9*(changeA)/y(nM,L) !O3 change in ppbv
c           endif
c           if(igas.eq.ntm)then
c            if(y(nOH,L).gt.0..and.y(nOH,L).lt.1.E20)tempJLS(J,L,nOH)=
c     *       tempJLS(J,L,nOH)+y(nOH,L)
c            tempJLS(J,L,nH2O)=tempJLS(J,L,nH2O)+(y(nH2O,L)/y(nM,L))
c           endif
c         endif
c
       end do    ! L
      end do     ! igas

 155  format(1x,2a4,a2,e13.3,a21,f10.0,a11,2x,e13.3,3x,a1,f12.5,a6)
 156  format(1x,2a4,a2,e13.3,a16)
 157  format(1x,2a4,a32,e13.3,a11)
c
      return
      end SUBROUTINE chemstep
cc    __________________________________________________________________
      SUBROUTINE rates(maxl,I,J)
!@sum rates calculate reaction rates with present concentrations
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on ds3ch4_chem_calc_jun1202_M23)
c
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only: nr,chemrate,photrate,rr,y,nn,dt2,
     &                          ss,ks,ny,dest,prod,JPPJ
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var kalt local dummy L-loop variable
!@var maxl passed highest chemistry level
!@var ireac,igas dummy loop variables
!@var I,J passed horixontal spatial indicies
      INTEGER kalt, ireac, igas, maxl
      INTEGER, INTENT(IN) :: I,J
c
      do kalt=1,maxl
        do ireac=1,nr
          chemrate(ireac,kalt)=rr(ireac,I,J,kalt)*y(nn(1,ireac),kalt)*
     &    y(nn(2,ireac),kalt)*dt2
        enddo
        do ireac=1,JPPJ
          photrate(ireac,kalt)=ss(ireac,I,J,kalt)*y(ks(ireac),kalt)*dt2
        enddo
c       Initialize change arrays
        do igas=1,ny
          dest(igas,kalt)=0.0
          prod(igas,kalt)=0.0
        enddo
      enddo
      return
      end SUBROUTINE rates

cc    __________________________________________________________________
      SUBROUTINE chem1(kdnr,maxl,numel,nn,ndnr,chemrate,dest,multip)
!@sum chem1 calculate chemical destruction/production
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on ds3ch4_chem_calc_jun1202_M23)
c
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only: p_2, p_3, p_4, n_igas, numfam,nfam
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var maxl passed highest chemistry level
!@var numel first index of nn array (shouldn't code like this anymore?)
!@var kdnr kdnr,kpnr,kds, or kps    passed from chemstep
!@var nn nn,nnr,ks, or kss          passed from chemstep
!@var ndnr ndnr,npnr,nds, or nps    passed from chemstep
!@var chemrate chemrate or photrate passed from chemstep
!@var dest dest or prod             passed from chemstep
!@var multip -1 for destruction, +1 for production
!@var i,ireac,igas,ial,nbeg,nend dummy loop variables
!@var dk dummy variable
      INTEGER ireac,igas,ial,i,dk,nbeg,nend
      INTEGER, INTENT(IN)            :: maxl, numel,multip
      INTEGER, DIMENSION(p_4)        :: kdnr
      INTEGER, DIMENSION(numel,p_2)  :: nn       ! can't do anymore?
      INTEGER, DIMENSION(p_3)        :: ndnr
      REAL*8,  DIMENSION(p_2,maxl)   :: chemrate ! can't do anymore?
      REAL*8,  DIMENSION(n_igas,maxl):: dest     ! can't do anymore?
C
      ireac=0
c     Reactive families
      do igas=1,numfam
        dk=kdnr(igas+1)-kdnr(igas)
        if(dk.ge.1) then
          do i=1,dk
           ireac=ireac+1
           do ial=1,maxl
            if(nn(1,ndnr(ireac)).ge.nfam(igas).and.nn(1,ndnr(ireac)).lt.
     *       nfam(igas+1))then
             dest(igas,ial)=dest(igas,ial)+multip
     &       *chemrate(ndnr(ireac),ial)
c            Save change array for individual family elements
             dest(nn(1,ndnr(ireac)),ial)=dest(nn(1,ndnr(ireac)),ial)+
     *          multip*chemrate(ndnr(ireac),ial)
            endif
            if(numel.eq.2)then
             if(nn(2,ndnr(ireac)).ge.nfam(igas).and.nn(2,ndnr(ireac))
     *        .lt.nfam(igas+1))then
              dest(igas,ial)=dest(igas,ial)+
     *        multip*chemrate(ndnr(ireac),ial)
              dest(nn(2,ndnr(ireac)),ial)=dest(nn(2,ndnr(ireac)),ial)+
     *        multip*chemrate(ndnr(ireac),ial)
             endif
            endif
           end do ! ial
          end do  ! i
        end if
      end do      ! igas
C
c     INDIVIDUAL SPECIES
      nbeg=numfam+1
      nend=nfam(1)-1
      do igas=nbeg,nend
        dk=kdnr(igas+1)-kdnr(igas)
        if(dk.ge.1) then
          do i=1,dk
           ireac=ireac+1
           do ial=1,maxl
            dest(igas,ial)=dest(igas,ial)+multip*
     &      chemrate(ndnr(ireac),ial)
           end do
          end do
        end if
      end do
      return
      end SUBROUTINE chem1

cc    __________________________________________________________________
      SUBROUTINE chem1prn(kdnr,numel,nn,ndnr,chemrate,
     &                    index,multip,igas,total,maxl,I,J)
!@sum chem1prn for printing out the chemical reactions
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on ds3ch4_chem_calc_jun1202_M23)
c
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only: ay, lprn, nfam, p_4, numfam, y,
     &                              p_2, p_3
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var kdnr kdnr,kpnr,kds, or kps from chemstep
!@var numel first index of nn array (shouldn't code like this anymore?)
!@var nn nn,nnr,ks, or kss from chemstep
!@var ndnr ndnr,npnr,nds, or nps from chemstep
!@var chemrate chemrate or photrate from chemstep
!@var index passed index to know which call this is... {1,2,3,4}
!@var multip 1 for production, -1 for destruction
!@var igas passed index for gas number
!@var total dummy summation
!@var maxl highest chemistry level
!@var I,J passed horizontal spatial indicies
!@var label character string for printing
!@var irec dummy loop variables
!@var per dummy temp variable
      INTEGER, INTENT(IN) :: igas,I,J,maxl,multip,index,numel
      INTEGER ireac
      REAL*8 total,per
      REAL*8,  DIMENSION(p_2,maxl)   :: chemrate !can't do this anymore?
      INTEGER, DIMENSION(p_3)        :: ndnr
      INTEGER, DIMENSION(numel,p_2)  :: nn       !can't do this anymore?
      INTEGER, DIMENSION(p_4)        :: kdnr
      character*17 label

c     skip this section during (within) family chemisty calls
      if(igas.le.numfam) then
c
c     FAMILIES
      if(kdnr(igas+1)-kdnr(igas).lt.1)goto40
      do ireac=kdnr(igas),kdnr(igas+1)-1
         if(index.le.2)label=' chem reaction # '
         if(index.gt.2)label=' phot reaction # '
         if(nn(1,ndnr(ireac)).ge.nfam(igas).and.nn(1,ndnr(ireac)).lt.
     *     nfam(igas+1))then
           per=0.
           if(y(igas,lprn).ne.0.)per=multip*100.*
     *     chemrate(ndnr(ireac),lprn)/y(igas,lprn)
           write(6,177) label,ndnr(ireac),' percent change from ',
     *     ay(1,nn(1,ndnr(ireac))),ay(2,nn(1,ndnr(ireac))),' = ',per,
     *     ' dy=',multip*chemrate(ndnr(ireac),lprn)
           total=total+per
         endif
         if(numel.eq.2)then
          if(nn(2,ndnr(ireac)).ge.nfam(igas).and.nn(2,ndnr(ireac)).lt.
     *     nfam(igas+1))then
           per=0.
           if(y(igas,lprn).ne.0.)per=multip*100.*
     *     chemrate(ndnr(ireac),lprn)/y(igas,lprn)
           write(6,177) label,ndnr(ireac),' percent change from ',
     *     ay(1,nn(2,ndnr(ireac))),ay(2,nn(2,ndnr(ireac))),' = ',per,
     *     ' dy=',multip*chemrate(ndnr(ireac),lprn)
           total=total+per
          endif
         endif  ! end of numel=2 loop
      end do  !end of ireac loop
      goto40

      end if ! non family chem only loop

c     INDIVIDUAL SPECIES
      if(kdnr(igas+1)-kdnr(igas).lt.1)goto40
      do ireac=kdnr(igas),kdnr(igas+1)-1
        if(index.le.2)label=' chem reaction # '
        if(index.gt.2)label=' phot reaction # '
c       skip same reaction if written twice
        if(ireac.eq.1 . or. ndnr(ireac).ne.ndnr(ireac-1) .or.
     &  total.eq.0.) then
         if(nn(1,ndnr(ireac)).eq.igas)then
           per=0.
           if(y(igas,lprn).ne.0.)per=100.*multip*
     *     chemrate(ndnr(ireac),lprn)/y(igas,lprn)
           write(6,106) label,ndnr(ireac),' percent change = ',per,
     *     ' dy=',multip*chemrate(ndnr(ireac),lprn)
           total=total+per
         endif
         if(numel.eq.2)then
          if(nn(2,ndnr(ireac)).eq.igas)then
           per=0.
           if(y(igas,lprn).ne.0.)per=100.*multip*
     *     chemrate(ndnr(ireac),lprn)/y(igas,lprn)
           write(6,106) label,ndnr(ireac),' percent change = ',per,
     *     ' dy=',multip*chemrate(ndnr(ireac),lprn)
           total=total+per
          endif
         endif  !end numel=2 loop
 106     format(a17,i3,a18,f10.0,a4,e12.3)
 177     format(a17,i3,a21,2a4,a3,f10.0,a4,e12.3)
  20    end if
  30  end do ! ireac
  40  continue
      return
      end SUBROUTINE chem1prn
