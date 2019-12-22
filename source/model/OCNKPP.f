#include "rundeck_opts.h"

      MODULE KPP_COM
!@sum  KPP_COM holds variables related to the KPP mixing scheme
!@auth Gavin Schmidt
!@ver  1.0
#ifdef TRACERS_OCEAN
      USE TRACER_COM, only : ntm
#endif
      USE OCEAN, only : im,jm,lmo
      USE SW2OCEAN, only : lsrpd
      IMPLICIT NONE
      SAVE
!@var KPL level to which mixed layer descends (1)
      INTEGER, DIMENSION(IM,JM) :: KPL

      REAL*8, DIMENSION(IM,JM,LSRPD) :: G0M1
      REAL*8, DIMENSION(IM,JM) :: MO1,GXM1,GYM1,S0M1,SXM1,SYM1,UO1,VO1
#ifdef TRACERS_OCEAN
      REAL*8, DIMENSION(NTM,IM,JM) :: TRMO1,TXMO1,TYMO1
#endif

      END MODULE KPP_COM
C****
      MODULE KPPE
!@sum  KPPE contains variables and routines for KPP mixing scheme
!@auth NCAR (modifications by Gavin Schmidt)
!@ver  1.0
c====================== include file "KPP_1D.COM" =====================
c
      USE OCEAN, only : lmo
      USE SW2OCEAN, only : lsrpd,fsr,fsrz,dfsrdz,dfsrdzb
      IMPLICIT NONE
      SAVE
c set max. number of ocean levels from ocean code
!@var km is max number of KPP layers
      INTEGER, PARAMETER :: km=LMO

c     rules for parameter constants
c
c use prefix of "c" for whole real numbers (ie: c57 for 57.0)
c use "m" after prefix to designate negative values (minus sign)
c      (ie: cm7 for -7.0)
c use prefix of "p" for non repeating fractions (ie: p5 for 0.5)
c use prefix of "r" for reciprocals (ie: r3 for 1/3.0)
c combine use of prefix above and "e" for scientific notation, with
c       (ie: c5e4 for 5.0e4, c1em10 for 1.0e-10)
c
      real*8,parameter :: c0=0.0, c1=1.0, c2=2.0, c4=4.0, c5=5.0, c8=8.0
      real*8,parameter :: c16=16.0, c360=360.0
      real*8,parameter :: p25=0.25, p5=0.5, p75=0.75
      real*8,parameter :: epsln=1.0e-20
c
      real*8,parameter :: c24=24.0, c60=60.0, c1440=1440.0
      real*8,parameter :: r24=c1/c24, r60=c1/c60, r1440=c1/c1440
      real*8,parameter :: secday=c1/(c60*c1440)
c
      real*8,parameter :: c1p5 = 1.5 , c3   = 3.     , c35   = 35.
      real*8,parameter :: c10 = 10. , c100 = 100.   , c1000 = 1000.
      real*8,parameter :: c10000 = 10000.
      real*8,parameter :: r3 = c1/c3 , r10  = c1/c10 , r100 = r10/c10
      real*8,parameter :: r1000  = c1/c1000, r10000 = r10*r1000

c====================== include file CVMIX_CLEAN.COM ==================

c     variables used for vertical diffusion
c
c inputs: (set through namelist)
c
!@var fkph   = vertical diffusion coefficient (cm**2/sec)
!@var fkpm   = vertical viscosity coefficient (cm**2/sec)
!@var bvdc   = background vertical diffusion constant
!@var bvvc   = background vertical viscosity constant
!@var vvclim = vertical viscosity coefficient limit
!@var vdclim = vertical diffusion coefficient limit
c
!@var vvcric = maximum viscosity   due to shear instability    (cm**2/s)
c              (based on local richardson number)
!@var vdcric = maximum diffusivity due to shear instability    (cm**2/s)
c
c arrays used for vertical diffusion in "k-profile" mixing scheme,
c computed in "kmix.F" and "kmixs.F".
c note, that in this scheme temperature diffusvity might be
c different from the diffusivity of all other tracers due to double
c diffusion.
c
      real*8, parameter :: vvcric=50d0, vdcric=50d0, fkpm=10d0, fkph=0
     *     .3d0,vvclim = 1000d0, vdclim = 1000d0

c====================== include file "KMIX_CLEAN.COM" =================
c     Define various parameters and common blocks for kmix vertical-
c     mixing scheme; used in "kmixs.F" subroutines
c
c parameters for several subroutines
c
!@var  epsl    = epsln in "pconst.h" (set in "kmixinit")    = 1.0e-20
!@var  epsilon = nondimensional extent of the surface layer = 0.1
!@var  vonk    = von Karman's constant                      = 0.4
!@var  conc1,conam,concm,conc2,zetam,conas,concs,conc3,zetas
c          = scalar coefficients

      real*8, parameter :: epsl=epsln, epsilon=0.1d0, vonk=0.4d0, conc1
     *     =5d0,conam=1.257d0, concm=8.380d0, conc2=16d0,zetam= -0.2d0
     *     ,conas=-28.86d0,concs=98.96d0,conc3=16d0,zetas= -1d0

c     parameters for subroutine "bldepth"

c to compute depth of boundary layer:
c
!@var  Ricr    = critical bulk Richardson Number            = 0.3
!@var  cekman  = coefficient for ekman depth                = 0.7
!@var  cmonob  = coefficient for Monin-Obukhov depth        = 1.0
!@var  concv   = ratio of interior buoyancy frequency to
c               buoyancy frequency at entrainment depth    = 1.8
!@var  hbf     = fraction of bounadry layer depth to
c               which absorbed solar radiation
c               contributes to surface buoyancy forcing    = 1.0
!@var  Vtc     = non-dimensional coefficient for velocity
c               scale of turbulant velocity shear
c               (=function of concv,concs,epsilon,vonk,Ricr)
c
      real*8, parameter :: Ricr= 0.3d0, cekman= 0.7d0, cmonob=1d0, concv
     *     =1.8d0,hbf=1d0
      real*8 Vtc
      common /kmixcbd/ Vtc

c     parameters and common arrays for subroutines "kmixinit" and
c     "wscale" to compute turbulent velocity scales:
c
!@var  nni     = number of values for zehat in the look up table
!@var  nnj     = number of values for ustar in the look up table
c
!@var  wmt     = lookup table for wm, the turbulent velocity scale
c               for momentum
!@var  wst     = lookup table for ws, the turbulent velocity scale
c               for scalars
!@var  deltaz  = delta zehat in table
!@var  deltau  = delta ustar in table
!@var  rdeltaz  = recipricol of delta zehat in table
!@var  rdeltau  = recipricol of delta ustar in table
!@var  zmin    = minimum limit for zehat in table (m3/s3)
!@var  zmax    = maximum limit for zehat in table
!@var  umin    = minimum limit for ustar in table (m/s)
!@var  umax    = maximum limit for ustar in table

      integer, parameter :: nni = 890, nnj = 480
c
      real*8 wmt,wst
      common /kmixcws/ wmt(0:nni+1,0:nnj+1),wst(0:nni+1,0:nnj+1)
      real*8, parameter :: zmin=-4d-7,zmax=0.,umin=0.,umax=4d-2,
     *     deltaz = (zmax-zmin)/(nni+1), deltau = (umax-umin)/(nnj+1),
     *     rdeltaz = 1./deltaz, rdeltau = 1./deltau
c
c     parameters for subroutine "ri_iwmix"
c     to compute vertical mixing coefficients below boundary layer:
c
!@var  Riinfty = local Richardson Number limit
c              for shear instability                      = 0.7
!@var  rRiinfty = 1/Riinfty
c    (note: the vertical mixing coefficients are defined
c    in (m2/s) units in "kmixinit". they are initialized
c    in "kmixbdta" and read via namelist "eddy" in (cm2/s)
c    units using their c-g-s names)
c   (m2/s)                                                 (cm2/s)
!@var difm0   = viscosity max due to shear instability     = vvcric
!@var difs0   = diffusivity ..                             = vdcric
!@var difmiw  = viscosity background due to internal waves = fkpm
!@var difsiw  = diffusivity ..                             = fkph
!@var difmcon = viscosity due to convective instability    = vvclim
!@var difscon = diffusivity ..                             = vdclim
!@var BVSQcon = value of N^2 where the convection coefs first become max
!@var rBVSQcon = 1/BVSQcon
!@var num_v_smooth_Ri = number of vertical smoothings of Ri

      real*8, parameter :: Riinfty= 0.7d0, rRiinfty=1./Riinfty, BVSQcon=
     *     -1d-7,rBVSQcon=1./BVSQcon, num_v_smooth_Ri=1
      real*8, parameter :: difm0   = vvcric * r10000, difs0   = vdcric
     *     * r10000,difmiw  = fkpm   * r10000, difsiw  = fkph   * r10000
     *     ,difmcon = vvclim * r10000, difscon = vdclim * r10000

c     parameters for subroutine "ddmix"
c     to compute additional diffusivity due to double diffusion:
!@var  Rrho0   = limit for double diffusive density ratio
!@var  dsfmax  = maximum diffusivity in case of salt fingering (m2/s)

      real*8, parameter :: Rrho0 = 1.9d0, dsfmax = 0.001d0

c     parameters for subroutine "blmix"
c     to compute mixing within boundary layer:
!@var  cstar   = proportionality coefficient for nonlocal transport
!@var  cg      = non-dimensional coefficient for counter-gradient term

      real*8, parameter :: cstar = 10d0
      real*8 cg
      common /kmixcbm/ cg

c add variables for depth dependent mixing due to rough topography
!@var diftop diffusion scale for topographic mixing
!@var fz500 diffusion scale as function of height above topography
      real*8, parameter :: diftop=0d0    ! 10d0 * r10000)
      real*8 fz500
      common /topomix/FZ500(km,km)

      END MODULE KPPE

      SUBROUTINE KPPMIX(LDD   , ZE    ,
     $                  zgrid , hwide , kmtj  , Shsq   , dVsq  ,
     $                  ustar , Bo    , Bosol , alphaDT, betaDS,
     $                  dbloc , Ritop , Coriol, byhwide,
     $                  visc  , difs  , dift  , ghats  , hbl
     $                                                 , kbl )

!@sum KPPMIX Main driver subroutine for kpp vertical mixing scheme and
!@+   interface to greater ocean model
c
c     written  by: bill large,    june  6, 1994
c     modified by: jan morzel,    june 30, 1994
c                  bill large,  august 11, 1994
c                  bill large, january 25, 1995 : "dVsq" and 1d code
c     modified for GISS by Gavin Schmidt, march 1998
c              for ModelE                 march 2001
c
      USE KPPE
!@var mdiff number of diffusivities for local arrays
      integer, parameter :: mdiff = 3
c input
      real*8 ZE(0:LMO)      !@var ZE GISS vertical layering (m)
      real*8 zgrid(0:km+1)  !@var zgrid vertical grid (<= 0) (m)
      real*8 hwide(0:km+1)  !@var hwide layer thicknesses    (m)
      real*8 byhwide(0:km+1)!@var byhwide 1/layer thicknesses (1/m)
      integer kmtj    !@var kmtj number of vertical layers on this row
      real*8 Shsq(km) !@var Shsq (local velocity shear)^2  (m/s)^2
      real*8 dVsq(km) !@var dVsq (velocity shear re sfc)^2 (m/s)^2
      real*8 ustar    !@var ustar surface friction velocity (m/s)
      real*8 Bo       !@var Bo surface turbulent buoy. forcing (m^2/s^3)
      real*8 Bosol    !@var Bosol radiative buoyancy forcing (m^2/s^3)
!@var alphaDT alpha * DT  across interfaces (kg/m^3)
      real*8 alphaDT(km)
!@var betaDS beta  * DS  across interfaces (kg/m^3)
      real*8 betaDS(km)
!@var dbloc local delta buoyancy across interfaces (m/s^2)
      real*8 dbloc(km)
!@var Ritop numerator of bulk Richardson Number (m/s)^2
c          Ritop = (-z - -zref)* delta buoyancy w/ respect to sfc
      real*8 Ritop(km)
      real*8 Coriol   !@var Coriol Coriolis parameter            (1/s)
      logical LDD     !@var LDD = TRUE for double diffusion
c output
!@var visc vertical viscosity coefficient (m^2/s)
      real*8 visc(0:km+1)
!@var difs vertical scalar diffusivity (m^2/s)
      real*8 difs(0:km+1)
!@var dift vertical temperature diffusivity (m^2/s)
      real*8 dift(0:km+1)
      real*8 ghats(km)  !@var ghats nonlocal transport (s/m^2)
      real*8 hbl        !@var hbl boundary layer depth (m)
      real*8 byhbl      !@var byhbl 1/boundary layer depth (1/m)
c local
      real*8 bfsfc      !@var bfsfc surface buoyancy forcing (m^2/s^3)
      real*8 ws         !@var ws momentum velocity scale
      real*8 wm         !@var wm scalar   velocity scale
      real*8 caseA      !@var caseA = 1 in case A; =0 in case B
      real*8 stable     !@var stable 1 in stable forcing; 0 in unstable
      real*8 dkm1(mdiff)!@var dkm1 boundary layer difs at kbl-1 level
      real*8 gat1(mdiff)!@var gat1 shape function at sigma=1
!@var dat1 derivative of shape function at sigma=1
      real*8 dat1(mdiff)
      real*8 blmc(km,mdiff)!@var blmc boundary layer mixing coefficients
      real*8 sigma      !@var sigma normalized depth (d / hbl)
      real*8 Rib(2)     !@var Rib bulk Richardson number
      integer kbl       !@var kbl index of first grid level below hbl
      integer kmax  !@var kmax minimum of LSRPD and kmtj, used in swfrac
c local ri_iwmix
      real*8 Rigg       !@var Rigg local richardson number
      real*8 fri,fcon   !@var fri,fcon function of Rig
      real*8 ftop       !@var ftop function of topography
c local enhance
      real*8 delta  !@var delta fraction hbl lies beteen zgrid neighbors
c local wscale
      real*8 zehat      !@var zehat = zeta *  ustar**3
      INTENT (IN) LDD,ZE,zgrid,hwide,kmtj,Shsq,dVsq,ustar,Bo,Bosol
     *            ,alphaDT,betaDS,dbloc,Ritop,Coriol,byhwide
      INTENT (OUT) visc, difs, dift, ghats, hbl, kbl
      INTEGER ki,mr,ka,ku,kl,iz,izp1,ju,jup1,ksave,kn,kt
      REAL*8 ratio,zdiff,zfrac,fzfrac,wam,wbm,was,wbs,u3,bvsq
     *     ,delhat,dvdzup,dvdzdn,viscp,diftp,visch,difsh,f1,bywm,byws
     *     ,sig,a1,a2,a3,gm,gs,gt,dstar,udiff,ufrac,vtsq,r,difsp,difth
     *     ,dkmp5

      kmax = MIN(LSRPD,kmtj)
c compute interior mixing coefficients everywhere, due to constant
c internal wave activity, static instability, and local shear
c instability.
c Zero surface values

c      call ri_iwmix ( km, km+1, kmtj, Shsq, dbloc, zgrid,
c    *                visc, difs , dift  )
c     compute interior viscosity diffusivity coefficients due
c     to shear instability (dependent on a local richardson number),
c     to background internal wave activity, and
c     to static instability (local richardson number < 0).

c     compute interior gradient Ri at all interfaces ki=1,km, (not surfa
ce)
c       use visc(ki=1,km) as temporary storage to be smoothed
c       use dift(ki=1,km) as temporary storage of N^2 = BVSQ
c       set values at bottom and below to nearest value above bottom

      do ki = 1, kmtj
         visc(ki)  = dbloc(ki) * (zgrid(ki)-zgrid(ki+1)) /
     $        ( Shsq(ki) + epsl)
         dift(ki) = dbloc(ki) * byhwide(ki+1)
      end do

c vertically smooth Ri num_v_smooth_Ri times

      do mr = 1,num_v_smooth_Ri
        call z121(visc,kmtj,km)
      end do

      do ki = 1, kmtj

c evaluate f of Vaisala squared    for convection        store in fcon
c evaluate f of   smooth Ri (fri) for shear instability store in fri

         Rigg  = DMAX1( dift(ki) , BVSQcon )
         ratio = DMIN1( (BVSQcon-Rigg) * rBVSQcon , c1 )
         fcon  = (c1 - ratio*ratio)
         fcon  = fcon * fcon * fcon

         Rigg  = DMAX1( visc(ki) , c0 )
         ratio = DMIN1( Rigg * rRiinfty , c1 )
         fri   = (c1 - ratio*ratio)
         fri   = fri  * fri   * fri
c add increased tracer mixing over rough topography. As a first cut,
c assume that topography is 'rough' if kmtj != LMO.
c functional form kv_top = diftop * exp (-z/500) where z is
c distance from the bottom.
c The array FZ500 is zero if kmtj=LMO

          ftop =  FZ500(ki,kmtj)

c evaluate diffusivities and viscosity
c mixing due to internal waves, and shear and static instability

         visc(ki) = (difmiw + fcon * difmcon + fri * difm0)
         difs(ki) = (difsiw + fcon * difscon + fri * difs0
     *                      + ftop * diftop)
         dift(ki) = difs(ki)
      end do

c set surface values to 0.0

      visc(0)    = c0
      dift(0)    = c0
      difs(0)    = c0

c add double diffusion if desired

      if (LDD) call ddmix (alphaDT,betaDS,visc,difs,dift,kmtj)

c Zero values at seafloor and below for blmix

      visc(kmtj:km+1) = c0
      difs(kmtj:km+1) = c0
      dift(kmtj:km+1) = c0

c compute boundary layer mixing coefficients:
c diagnose the new boundary layer depth

c       call bldepth (km, km+1, zgrid, hwide, kmtj, dVsq,
c      $             dbloc, Ritop, ustar, Bo, Bosol, Coriol,
c      $             hbl, bfsfc, stable, caseA, kbl,
c      $             Rib, sigma, wm, ws)

c     the oceanic planetray boundary layer depth, hbl, is determined as
c     the shallowest depth where the bulk richardson number is
c     equal to the critical value, Ricr.
c     Bulk richardson numbers are evaluated by computing velocity and
c     buoyancy differences between values at zgrid(kl) < 0 and surface
c     reference values.
c     in this configuration, the reference values are equal to the
c     values in the surface layer.
c     when using a very fine vertical grid, these values should be
c     computed as the vertical average of velocity and buoyancy from
c     the surface down to epsilon*zgrid(kl).
c     When the bulk richardson number at k exceeds Ricr, hbl is
c     linearly interpolated between grid levels zgrid(k) and zgrid(k-1).
c     The water column and the surface forcing are diagnosed for
c     stable/ustable forcing conditions, and where hbl is relative
c     to grid points (caseA), so that conditional branches can be
c     avoided in later subroutines.

c find bulk Richardson number at every grid level until > Ricr
c
c note: the reference depth is -epsilon/2.*zgrid(k), but the reference
c       u,v,t,s values are simply the surface layer values,
c       and not the averaged values from 0 to 2*ref.depth,
c       which is necessary for very fine grids(top layer < 2m thickness)
c note: max values when Ricr never satisfied are
c       kbl=kmtj and hbl=-zgrid(kmtj)

c indices for array Rib(k), the bulk Richardson number.
      ka = 1
      ku = 2

c initialize hbl and kbl to bottomed out values
      Rib(ka) = 0.0
      kbl    = kmtj
      hbl    = -zgrid(kbl)

      kl = 1
 10   kl = kl + 1

c compute bfsfc = sw fraction at hbf * zgrid
c      call swfrac(zgrid(kl),ZE,kl,kmtj,kmax,bfsfc)  ! inlined

C**** Get actual k (since k is on staggered grid)
      kt = kl + NINT(0.5d0 + SIGN(0.5d0,-(ZE(kl)+zgrid(kl))))

C**** calculate fraction
      if (kt.gt.kmax) then
         bfsfc = 0.
      else
         if (kt.eq.kmax) then
            bfsfc = FSR(kt) + (zgrid(kl)+ZE(kt-1))*dFSRdZB(kt)
         else
            bfsfc = FSR(kt) + (zgrid(kl)+ZE(kt-1))*dFSRdZ(kt)
         end if
      end if

c use caseA as temporary array
         caseA  = -zgrid(kl)

c compute bfsfc= Bo + radiative contribution down to hbf * hbl
         bfsfc  = Bo + Bosol * (1. - bfsfc)

         stable = 0.5 + SIGN( 5d-1, bfsfc )
         sigma  = stable * 1. + (1.-stable) * epsilon

c compute velocity scales at sigma, for hbl= caseA = -zgrid(kl)
c         call wscale(sigma, caseA, ustar, bfsfc,   wm, ws)
c     compute turbulent velocity scales.
c     use a 2D-lookup table for wm and ws as functions of ustar and
c     zetahat (=vonk*sigma*hbl*bfsfc).
c use lookup table for zehat < zmax only; otherwise use
c stable formulae

      zehat = vonk * sigma * caseA * bfsfc

      if (zehat.le.zmax) then
         zdiff  = zehat-zmin
         iz = int( zdiff * rdeltaz )
         iz = min( iz , nni )
         iz = max( iz , 0  )
         izp1=iz+1

         udiff  = ustar-umin
         ju = int( udiff * rdeltau)
         ju = min( ju , nnj )
         ju = max( ju , 0  )
         jup1=ju+1

         zfrac = zdiff*rdeltaz - float(iz)
         ufrac = udiff*rdeltau - float(ju)

         fzfrac= 1.-zfrac
         wam   = (fzfrac)  * wmt(iz,jup1) + zfrac*wmt(izp1,jup1)
         wbm   = (fzfrac)  * wmt(iz,ju  ) + zfrac*wmt(izp1,ju  )
         wm = (1.-ufrac)* wbm          + ufrac*wam

         was   = (fzfrac)  * wst(iz,jup1) + zfrac*wst(izp1,jup1)
         wbs   = (fzfrac)  * wst(iz,ju  ) + zfrac*wst(izp1,ju  )
         ws = (1.-ufrac)* wbs          + ufrac*was
      else
         u3    = ustar*ustar*ustar
         wm = vonk * ustar * u3 / ( u3 + conc1*zehat )
         ws = wm
      endif

c compute the turbulent shear contribution to Rib
         bvsq =0.5*
     $        ( dbloc(kl-1) * byhwide(kl) +
     $          dbloc(kl  ) * byhwide(kl+1) )
         Vtsq = - zgrid(kl) * ws * sqrt(abs(bvsq)) * Vtc

c           compute bulk Richardson number at new level, dunder
c           note: Ritop needs to be zero on land and ocean bottom
c           points so that the following if statement gets triggered
c           correctly. otherwise, hbl might get set to (big) negative
c           values, that might exceed the limit for the "exp" function
c           in "swfrac"

         Rib(ku) = Ritop(kl) / (dVsq(kl)+Vtsq+epsl)

c           linearly interpolate to find hbl where Rib = Ricr
         if((kbl.eq.kmtj).and.(Rib(ku).gt.Ricr)) then
            hbl = -zgrid(kl-1) + (zgrid(kl-1)-zgrid(kl)) *
     $           (Ricr - Rib(ka)) / (Rib(ku)-Rib(ka))
            kbl = kl
         else
            ksave = ka
            ka    = ku
            ku    = ksave
            if (kl.lt.kmtj) go to 10
         end if

c find stability and buoyancy forcing for boundary layer
c      call swfrac(-hbl,ZE,kbl-1,kmtj,bfsfc,kmax)  ! inlined

C**** Get actual k (since k is on staggered grid)
      kt = kbl-1+ NINT(0.5d0 + SIGN(0.5d0,-(ZE(kbl-1)-hbl)))

C**** calculate fraction
      if (kt.gt.kmax) then
         bfsfc = 0.
      else
         if (kt.eq.kmax) then
            bfsfc = FSR(kt) + (ZE(kt-1)-hbl)*dFSRdZB(kt)
         else
            bfsfc = FSR(kt) + (ZE(kt-1)-hbl)*dFSRdZ(kt)
         end if
      end if

      bfsfc  = Bo + Bosol * (1. - bfsfc)
      stable = 0.5 + SIGN( 5d-1, bfsfc )
      bfsfc  = bfsfc + stable * epsl !ensures bfsfc never=0

c check hbl limits for hekman or hmonob (NOT USED)

c      if(bfsfc.gt.0.0) then
c         hekman = cekman * ustar / (abs(Coriol)+epsl)
c         hmonob = cmonob * ustar*ustar*ustar
c     $        /(vonk * (bfsfc+epsl) )
c         hlimit = stable     * DMIN1(hekman,hmonob) +
c     $        (stable-1.) * zgrid(km)
c         hbl = DMIN1(hbl,hlimit)
c         hbl = DMAX1(hbl,-zgrid(1))
c      endif
c      kbl = kmtj

c find new kbl
c      do kl=2,kmtj
c         if((kbl.eq.kmtj).and.(-zgrid(kl).gt.hbl)) kbl = kl
c      end do

c find stability and buoyancy forcing for final hbl values
c      call swfrac(-hbl,ZE,kbl-1,kmtj,kmax,bfsfc)

c      bfsfc  = Bo + Bosol * (1. - bfsfc)
c      stable = 0.5 + SIGN( 5d-1, bfsfc )
c      bfsfc  = bfsfc + stable * epsl

c determine caseA and caseB
      caseA  = 0.5 + SIGN(5d-1,-zgrid(kbl) - 0.5*hwide(kbl) - hbl)

      byhbl = 1d0/hbl
c compute boundary layer diffusivities

c       call blmix    (km   , km+1 , mdiff , zgrid, hwide ,
c      $               ustar, bfsfc, hbl  , stable, caseA,
c      $               visc , difs , dift , kbl   ,
c      $               gat1 , dat1 , dkm1 , blmc  , ghats,
c      $               sigma, wm   , ws   )
c     mixing coefficients within boundary layer depend on surface
c     forcing and the magnitude and gradient of interior mixing below
c     the boundary layer ("matching").
c     Caution: if mixing bottoms out at hbl = -zgrid(km) then
c     fictitious layer at km+1 is needed with small but finite width
c     hwide(km+1) (eg. epsl = 1.e-20).

c compute velocity scales at hbl

      sigma = stable * 1.0 + (1.-stable) * epsilon

c      call wscale(sigma, hbl, ustar, bfsfc,   wm, ws)
c     compute turbulent velocity scales.
c     use a 2D-lookup table for wm and ws as functions of ustar and
c     zetahat (=vonk*sigma*hbl*bfsfc).
c use lookup table for zehat < zmax only; otherwise use
c stable formulae

      zehat = vonk * sigma * hbl * bfsfc

      if (zehat.le.zmax) then
         zdiff  = zehat-zmin
         iz = int( zdiff * rdeltaz )
         iz = min( iz , nni )
         iz = max( iz , 0  )
         izp1=iz+1

         udiff  = ustar-umin
         ju = int( udiff * rdeltau)
         ju = min( ju , nnj )
         ju = max( ju , 0  )
         jup1=ju+1

         zfrac = zdiff*rdeltaz - float(iz)
         ufrac = udiff*rdeltau - float(ju)

         fzfrac= 1.-zfrac
         wam   = (fzfrac)  * wmt(iz,jup1) + zfrac*wmt(izp1,jup1)
         wbm   = (fzfrac)  * wmt(iz,ju  ) + zfrac*wmt(izp1,ju  )
         wm = (1.-ufrac)* wbm          + ufrac*wam

         was   = (fzfrac)  * wst(iz,jup1) + zfrac*wst(izp1,jup1)
         wbs   = (fzfrac)  * wst(iz,ju  ) + zfrac*wst(izp1,ju  )
         ws = (1.-ufrac)* wbs          + ufrac*was
      else
         u3    = ustar*ustar*ustar
         wm = vonk * ustar * u3 / ( u3 + conc1*zehat )
         ws = wm
      endif

      kn    = int(caseA+epsl) *(kbl -1) +
     $     (1-int(caseA+epsl)) * kbl

c find the interior viscosities and derivatives at hbl
      delhat = 0.5*hwide(kn) - zgrid(kn) - hbl
      R      = 1.0 - delhat * byhwide(kn)
      dvdzup = (visc(kn-1) - visc(kn)) * byhwide(kn)
      dvdzdn = (visc(kn)   - visc(kn+1)) * byhwide(kn+1)
      viscp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+
     $     R  * (dvdzdn + abs(dvdzdn)) )

      dvdzup = (difs(kn-1) - difs(kn)) * byhwide(kn)
      dvdzdn = (difs(kn)   - difs(kn+1)) * byhwide(kn+1)
      difsp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+
     $     R  * (dvdzdn + abs(dvdzdn)) )

      dvdzup = (dift(kn-1) - dift(kn)) * byhwide(kn)
      dvdzdn = (dift(kn)   - dift(kn+1)) * byhwide(kn+1)
      diftp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+
     $     R  * (dvdzdn + abs(dvdzdn)) )

      visch  = visc(kn) + viscp * delhat
      difsh  = difs(kn) + difsp * delhat
      difth  = dift(kn) + diftp * delhat

      f1 = stable * conc1 * bfsfc / (ustar**4+epsl)

      bywm = 1./(wm+epsl)
      byws = 1./(ws+epsl)

      gat1(1) = visch * byhbl * bywm
      dat1(1) = -viscp * bywm + f1 * visch
      dat1(1) = min(dat1(1),0d0)

      gat1(2) = difsh  * byhbl * byws
      dat1(2) = -difsp * byws  + f1 * difsh
      dat1(2) = min(dat1(2),0d0)

      gat1(3) = difth * byhbl * byws
      dat1(3) = -diftp * byws + f1 * difth
      dat1(3) = min(dat1(3),0d0)

      do ki = 1,kbl-1

c compute turbulent velocity scales on the interfaces
         sig     = (-zgrid(ki) + 0.5 * hwide(ki)) * byhbl
         sigma= stable*sig + (1.-stable)*DMIN1(sig,epsilon)

c         call wscale(sigma, hbl, ustar, bfsfc,   wm,  ws)
c     compute turbulent velocity scales.
c     use a 2D-lookup table for wm and ws as functions of ustar and
c     zetahat (=vonk*sigma*hbl*bfsfc).
c use lookup table for zehat < zmax only; otherwise use
c stable formulae

      zehat = vonk * sigma * hbl * bfsfc

      if (zehat.le.zmax) then
         zdiff  = zehat-zmin
         iz = int( zdiff * rdeltaz )
         iz = min( iz , nni )
         iz = max( iz , 0  )
         izp1=iz+1

         udiff  = ustar-umin
         ju = int( udiff * rdeltau)
         ju = min( ju , nnj )
         ju = max( ju , 0  )
         jup1=ju+1

         zfrac = zdiff*rdeltaz - float(iz)
         ufrac = udiff*rdeltau - float(ju)

         fzfrac= 1.-zfrac
         wam   = (fzfrac)  * wmt(iz,jup1) + zfrac*wmt(izp1,jup1)
         wbm   = (fzfrac)  * wmt(iz,ju  ) + zfrac*wmt(izp1,ju  )
         wm = (1.-ufrac)* wbm          + ufrac*wam

         was   = (fzfrac)  * wst(iz,jup1) + zfrac*wst(izp1,jup1)
         wbs   = (fzfrac)  * wst(iz,ju  ) + zfrac*wst(izp1,ju  )
         ws = (1.-ufrac)* wbs          + ufrac*was
      else
         u3    = ustar*ustar*ustar
         wm = vonk * ustar * u3 / ( u3 + conc1*zehat )
         ws = wm
      endif

c compute the dimensionless shape functions at the interfaces
         sig = (-zgrid(ki) + 0.5 * hwide(ki)) * byhbl
         a1 = sig - 2.
         a2 = 3.-2.*sig
         a3 = sig - 1.

         Gm = a1 + a2 * gat1(1) + a3 * dat1(1)
         Gs = a1 + a2 * gat1(2) + a3 * dat1(2)
         Gt = a1 + a2 * gat1(3) + a3 * dat1(3)

c compute boundary layer diffusivities at the interfaces
         blmc(ki,1) = hbl * wm * sig * (1. + sig * Gm)
         blmc(ki,2) = hbl * ws * sig * (1. + sig * Gs)
         blmc(ki,3) = hbl * ws * sig * (1. + sig * Gt)

c nonlocal transport term = ghats * <ws>o
         ghats(ki) = (1.-stable) * cg * byws * byhbl
      end do

c find diffusivities at kbl-1 grid level
      sig      =  -zgrid(kbl-1)  * byhbl
      sigma =stable * sig + (1.-stable) * MIN(sig,epsilon)

c      call wscale(sigma, hbl, ustar, bfsfc,   wm, ws)
c     compute turbulent velocity scales.
c     use a 2D-lookup table for wm and ws as functions of ustar and
c     zetahat (=vonk*sigma*hbl*bfsfc).
c use lookup table for zehat < zmax only; otherwise use
c stable formulae

      zehat = vonk * sigma * hbl * bfsfc

      if (zehat.le.zmax) then
         zdiff  = zehat-zmin
         iz = int( zdiff * rdeltaz )
         iz = min( iz , nni )
         iz = max( iz , 0  )
         izp1=iz+1

         udiff  = ustar-umin
         ju = int( udiff * rdeltau)
         ju = min( ju , nnj )
         ju = max( ju , 0  )
         jup1=ju+1

         zfrac = zdiff*rdeltaz - float(iz)
         ufrac = udiff*rdeltau - float(ju)

         fzfrac= 1.-zfrac
         wam   = (fzfrac)  * wmt(iz,jup1) + zfrac*wmt(izp1,jup1)
         wbm   = (fzfrac)  * wmt(iz,ju  ) + zfrac*wmt(izp1,ju  )
         wm = (1.-ufrac)* wbm          + ufrac*wam

         was   = (fzfrac)  * wst(iz,jup1) + zfrac*wst(izp1,jup1)
         wbs   = (fzfrac)  * wst(iz,ju  ) + zfrac*wst(izp1,ju  )
         ws = (1.-ufrac)* wbs          + ufrac*was
      else
         u3    = ustar*ustar*ustar
         wm = vonk * ustar * u3 / ( u3 + conc1*zehat )
         ws = wm
      endif

      sig = -zgrid(kbl-1) * byhbl
      a1= sig - 2.
      a2 = 3.-2.*sig
      a3 = sig - 1.
      Gm = a1 + a2 * gat1(1) + a3 * dat1(1)
      Gs = a1 + a2 * gat1(2) + a3 * dat1(2)
      Gt = a1 + a2 * gat1(3) + a3 * dat1(3)
      dkm1(1) = hbl * wm * sig * (1. + sig * Gm)
      dkm1(2) = hbl * ws * sig * (1. + sig * Gs)
      dkm1(3) = hbl * ws * sig * (1. + sig * Gt)

c       call enhance  (km   , km+1 , mdiff , dkm1 , visc  ,
c      $               difs , dift , hbl  , kbl   , zgrid, caseA ,
c      $               blmc ,ghats )
c     enhance the diffusivity at the kbl-.5 interface

      if (kbl.le.kmtj) then
         ki = kbl-1
         delta = (hbl+zgrid(ki)) * byhwide(ki+1)

         dkmp5 = caseA*visc(ki) + (1.-caseA)*blmc(ki,1)
         dstar = (1.-delta)**2 * dkm1(1) + delta**2 * dkmp5
         blmc(ki,1) = (1.-delta) * visc(ki) + delta * dstar

         dkmp5 = caseA*difs(ki) + (1.-caseA)*blmc(ki,2)
         dstar = (1.-delta)**2 * dkm1(2) + delta**2 * dkmp5
         blmc(ki,2) = (1.-delta) * difs(ki) + delta * dstar

         dkmp5 = caseA*dift(ki) + (1.-caseA)*blmc(ki,3)
         dstar = (1.-delta)**2 * dkm1(3) + delta**2 * dkmp5
         blmc(ki,3) = (1.-delta) * dift(ki) + delta * dstar

         ghats(ki) = (1.-caseA) * ghats(ki)
      end if

c combine interior and boundary layer coefficients and nonlocal term

      visc(1:kbl-1)=blmc(1:kbl-1,1)
      difs(1:kbl-1)=blmc(1:kbl-1,2)
      dift(1:kbl-1)=blmc(1:kbl-1,3)
      ghats(kbl:km)=0.

      return
      end subroutine kppmix

      subroutine wscale(sigma, hbl, ustar, bfsfc, wm , ws)

c     compute turbulent velocity scales.
c     use a 2D-lookup table for wm and ws as functions of ustar and
c     zetahat (=vonk*sigma*hbl*bfsfc).
c
c     note: the lookup table is only used for unstable conditions
c     (zehat.le.0), in the stable domain wm (=ws) gets computed
c     directly.
      USE KPPE
      INTENT (IN) sigma,hbl,ustar,bfsfc
      INTENT (OUT) wm,ws

c  input
      real*8 sigma      ! normalized depth (d/hbl)
      real*8 hbl        ! boundary layer depth (m)
      real*8 ustar      ! surface friction velocity         (m/s)
      real*8 bfsfc    ! total surface buoyancy flux       (m^2/s^3)
c  output
      real*8 wm,ws ! turbulent velocity scales at sigma
c local
      real*8 zehat           ! = zeta *  ustar**3
      real*8 zdiff,udiff,ufrac,fzfrac,wam,was,wbs,wbm,u3,zfrac
      integer iz,izp1,ju,jup1
c use lookup table for zehat < zmax only; otherwise use
c stable formulae

      zehat = vonk * sigma * hbl * bfsfc

      if (zehat.le.zmax) then
         zdiff  = zehat-zmin
         iz = int( zdiff * rdeltaz )
         iz = min( iz , nni )
         iz = max( iz , 0  )
         izp1=iz+1

         udiff  = ustar-umin
         ju = int( udiff * rdeltau)
         ju = min( ju , nnj )
         ju = max( ju , 0  )
         jup1=ju+1

         zfrac = zdiff*rdeltaz - float(iz)
         ufrac = udiff*rdeltau - float(ju)

         fzfrac= 1.-zfrac
         wam   = (fzfrac)  * wmt(iz,jup1) + zfrac*wmt(izp1,jup1)
         wbm   = (fzfrac)  * wmt(iz,ju  ) + zfrac*wmt(izp1,ju  )
         wm = (1.-ufrac)* wbm          + ufrac*wam

         was   = (fzfrac)  * wst(iz,jup1) + zfrac*wst(izp1,jup1)
         wbs   = (fzfrac)  * wst(iz,ju  ) + zfrac*wst(izp1,ju  )
         ws = (1.-ufrac)* wbs          + ufrac*was
      else
         u3    = ustar*ustar*ustar
         wm = vonk * ustar * u3 / ( u3 + conc1*zehat )
         ws = wm
      endif

      return
      end subroutine wscale

      subroutine ddmix (alphaDT, betaDS, visc, difs, dift, kmtj)

c     Rrho dependent interior flux parameterization.
c     Add double-diffusion diffusivities to Ri-mix values at blending
c     interface and below.

      USE KPPE
      INTENT (IN) alphaDT,betaDS,kmtj
      INTENT (INOUT) visc, difs, dift
c input
      real*8 alphaDT(km)   ! alpha * DT  across interfaces
      real*8 betaDS(km)    ! beta  * DS  across interfaces
c output
      real*8 visc(0:km+1)  ! interior viscosity           (m^2/s)
      real*8 dift(0:km+1)  ! interior thermal diffusivity (m^2/s)
      real*8 difs(0:km+1)  ! interior scalar  diffusivity (m^2/s)
c local
      real*8 diffdd            ! double diffusion diffusivity scale
      real*8 prandtl           ! prandtl number
      integer kmtj,ki
      real*8 Rrho

      do 100 ki= 1, kmtj

c salt fingering case

         if((alphaDT(ki).gt.betaDS(ki)).and.
     $        (betaDS (ki).gt.0.         )) then

            Rrho       = MIN(alphaDT(ki) / betaDS(ki) , Rrho0)
c           diffdd     = dsfmax*(1.0-((Rrho-1)/(Rrho0-1))**2)**pexp2
            diffdd     =         1.0-((Rrho-1)/(Rrho0-1))**2
            diffdd     = dsfmax*diffdd*diffdd*diffdd
            dift(ki) = dift(ki) + 0.7*diffdd
            difs(ki) = difs(ki) + diffdd

c diffusive convection

         else if ((alphaDT(ki).lt.0.0).and.(betaDS(ki).lt.0.0)
     $           .and.(alphaDT(ki).lt.betaDS(ki)) ) then

            Rrho    = alphaDT(ki) / betaDS(ki)
            diffdd=1.5d-6*9.0*0.101d0*exp(4.6d0*exp(-0.54d0*(1/Rrho-1)))
            prandtl = 0.15*Rrho
            if (Rrho.gt.0.5d0) prandtl = (1.85d0-0.85d0/Rrho)*Rrho
            dift(ki) = dift(ki) + diffdd
            difs(ki) = difs(ki) + prandtl*diffdd
         endif
 100  continue
      return
      end subroutine ddmix

      subroutine kmixinit(ZE)
c     initialize some constants for kmix subroutines, and initialize
c     for kmix subroutine "wscale" the 2D-lookup table for wm and ws
c     as functions of ustar and zetahat (=vonk*sigma*hbl*bfsfc).

      USE KPPE
c local
      real*8 zehat                        ! = zeta *  ustar**3
      real*8 zeta                       ! = stability parameter d/L
      real*8 ZE(0:LMO)                  ! GCM vertical grid
      integer lbot,j,i,l
      real*8 usta

c define some non-dimensional constants and
c the vertical mixing coefficients in m-k-s units

c      epsl    = epsln
      Vtc     = concv * sqrt(0.2/concs/epsilon) / vonk**2 / Ricr
      cg      = cstar * vonk * (concs * vonk * epsilon)**(1./3.)

c      difm0   = vvcric * r10000
c      difs0   = vdcric * r10000
c      difmiw  = fkpm   * r10000
c      difsiw  = fkph   * r10000
c      difmcon = vvclim * r10000
c      difscon = vdclim * r10000

c construct the wm and ws lookup tables

c      deltaz = (zmax-zmin)/(nni+1)
c      deltau = (umax-umin)/(nnj+1)
c      rdeltaz = 1./deltaz
c      rdeltau = 1./deltau

      do 100 i=0,nni+1
         zehat = deltaz*(i) + zmin
         do 90 j=0,nnj+1
            usta = deltau*(j) + umin
            zeta = zehat/(usta**3+epsl)

            if(zehat.ge.0.) then
               wmt(i,j) = vonk*usta/(1.+conc1*zeta)
               wst(i,j) = wmt(i,j)
            else
               if(zeta.gt.zetam) then
                  wmt(i,j) = vonk* usta * sqrt(sqrt(1.-conc2*zeta))
               else
                  wmt(i,j) = vonk* (conam*usta**3-concm*zehat)**(1./3.)
               endif
               if(zeta.gt.zetas) then
                  wst(i,j) = vonk* usta * sqrt(1.-conc3*zeta)
               else
                  wst(i,j) = vonk* (conas*usta**3-concs*zehat)**(1./3.)
               endif
            endif
 90      continue
 100  continue

c  set up depth dependent factor for mixing due to rough topography
c  assume FZ500 = exp (-z/500)
c  for two levels above bottom, starting at level km/2

      FZ500=0.
      do lbot=LMO/2,LMO-1
         do l=lbot-2,lbot-1
            FZ500(l,lbot)=exp(-(ZE(lbot)-ZE(l))/500.)
         end do
      end do

      return
      end subroutine kmixinit

      subroutine swfrac(z,ZE,k,kmtj,kmax,bfsfc)
!@sum swfrac Calculate fraction of solar energy penetrating to depth z
!@+   using linear interpolation for depths between grid levels
!@+   There is a slight error since 'z' is scaled by free surface height
!@+   and ZE is not.
C****
!@var  k is grid box containing point
!@var  kmtj is number of grid points in column
!@var  kmax is the number of grid points that recieve solar radiation
C****
      USE KPPE
      INTENT (IN) z,k,kmtj,kmax,ZE
      INTENT (OUT) bfsfc
      real*8 z,bfsfc
      real*8 ZE(0:LMO)
      integer k,l,kmtj,kmax,kt

C**** Get actual k (since k is on staggered grid)
      kt = k + NINT(0.5d0 + SIGN(0.5d0,-(ZE(k)+z)))

C**** calculate fraction
      if (kt.gt.kmax) then
         bfsfc = 0.
      else
         if (kt.eq.kmax) then
            bfsfc = FSR(kt) + (z+ZE(kt-1))*dFSRdZB(kt)
         else
            bfsfc = FSR(kt) + (z+ZE(kt-1))*dFSRdZ(kt)
         end if
      end if

      return
      end subroutine swfrac

      subroutine z121 (v,kmtj,km)
!@sum z121 Apply 121 smoothing in k to 2-d array V(k=1,km)
!@+   top (0) value is used as a dummy
!@+   bottom (km+1) value is set to input value from above.
      IMPLICIT NONE
      REAL*8, PARAMETER :: p5=5d-1, p25=2.5d-1
      INTEGER, INTENT (IN) :: kmtj,km
      REAL*8, INTENT (INOUT) :: V(0:km+1)  ! 2-D array to be smoothed
      INTEGER K
      REAL*8 tmp

      V(0)      =  p25 * V(1)
      V(kmtj+1) =        V(kmtj)

      do k=1,kmtj
         tmp      =  V(k)
         V(k)   =  V(0)  + p5 * V(k) + p25 * V(k+1)
         V(0)   =  p25 * tmp
      end do
      return
      end subroutine z121

C**** Is this still necessary now that fluxes are saved?
      SUBROUTINE KVINIT
!@sum KVINIT Initialise KMIX and save pre-source term surface values of
!@+   enthalpy, mass and horizontal gradients for kppmix calculation
      USE OCEAN, only : im,jm,lmo,gxmo,sxmo,gymo,symo,g0m,s0m,mo,uo,vo
#ifdef TRACERS_OCEAN
     *     ,trmo,txmo,tymo
#endif
      USE KPP_COM, only : G0M1,MO1,GXM1,GYM1,S0M1,SXM1,SYM1,UO1,VO1
     *     ,lsrpd
#ifdef TRACERS_OCEAN
     *     ,trmo1,txmo1,tymo1
#endif
      IMPLICIT NONE
      INTEGER I,J
C**** Save surface values
      DO J=1,JM
         DO I=1,IM
            G0M1(I,J,1:LSRPD) = G0M(I,J,1:LSRPD)
            GXM1(I,J) = GXMO(I,J,1)
            GYM1(I,J) = GYMO(I,J,1)
            S0M1(I,J) = S0M(I,J,1)
            SXM1(I,J) = SXMO(I,J,1)
            SYM1(I,J) = SYMO(I,J,1)
            MO1(I,J)  = MO(I,J,1)
            UO1(I,J)  = UO(I,J,1)
            VO1(I,J)  = VO(I,J,1)
#ifdef TRACERS_OCEAN
            TRMO1(:,I,J) = TRMO(I,J,1,:)
            TXMO1(:,I,J) = TXMO(I,J,1,:)
            TYMO1(:,I,J) = TYMO(I,J,1,:)
#endif
         END DO
      END DO
      RETURN
      END

      SUBROUTINE OCONV
!@sum  OCONV does vertical mixing using coefficients from KPP scheme
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0
      USE CONSTANT, only : grav,omega
#ifdef TRACERS_OCEAN
      USE TRACER_COM, only : t_qlimit
      USE OCEAN, only : trmo,txmo,tymo,tzmo,ntm
#endif
      USE OCEAN, only : im,jm,lmo,g0m,s0m,gxmo,sxmo,symo,gymo,szmo,gzmo
     *     ,ogeoz,hocean,ze,bydxypo,mo,sinpo,dts,lmm,lmv,lmu,ramvs
     *     ,dxypo,cosic,sinic,uo,vo,ramvn,bydts
      USE SEAICE_COM, only : rsi
      USE ODIAG, only : oijl,oij,ij_hbl,ij_bo,ij_bosol,ij_ustar,ijl_kvm
     *     ,ijl_kvg,ijl_wgfl,ijl_wsfl,ol,l_rho,l_temp,l_salt  !ij_ogeoz
      USE KPP_COM, only : g0m1,s0m1,mo1,gxm1,gym1,sxm1,sym1,uo1,vo1,kpl
#ifdef TRACERS_OCEAN
     *     ,trmo1,txmo1,tymo1
#endif
      USE FLUXES, only : solar,dmua,dmva,dmui,dmvi
      USE SW2OCEAN, only : fsr,lsrpd
      IMPLICIT NONE

      LOGICAL*4 QPOLE
      REAL*8 UT(IM,JM,LMO),VT(IM,JM,LMO),UL(LMO,IM+2),MML(LMO),
     *     G0ML(LMO,2,2),S0ML(LMO,2,2),GZML(LMO,2,2),SZML(LMO,2,2),
     *     BYMML(LMO),DTBYDZ(LMO),BYDZ2(LMO),RAVM(IM+1),RAMV(IM+1),
     *     UL0(LMO,IM+2),G0ML0(LMO,2,2),S0ML0(LMO,2,2),MML0(LMO),
     *     BYMML0(LMO),MMLT(LMO),BYMMLT(LMO),
     *     AKVM(0:LMO+1),AKVG(0:LMO+1),AKVS(0:LMO+1),GHATM(LMO),
     *     GHATG(LMO),GHATS(LMO),FLG(LMO),FLS(LMO),TXY
#ifdef TRACERS_OCEAN
     *     ,TRML(LMO,NTM,2,2),TZML(LMO,NTM,2,2),TRML1(NTM,2,2)
     *     ,DELTATR(NTM),GHATT(LMO,NTM),FLT(LMO,NTM)
      INTEGER NSIGT,N
#endif
      INTEGER LMUV(IM+1)

C**** CONV parameters: BETA controls degree of convection (default 0.5).
      REAL*8, PARAMETER :: BETA=5d-1, BYBETA=1d0/BETA
C**** KPP variables
      REAL*8, PARAMETER :: epsln=1d-20
      REAL*8 zgrid(0:LMO+1),hwide(0:LMO+1),Shsq(LMO),dVsq(LMO)
     *     ,talpha(LMO),sbeta(LMO),dbloc(LMO),dbsfc(LMO),Ritop(LMO),
     *     alphaDT(LMO),betaDS(LMO),ghat(LMO),byhwide(0:LMO+1)
      REAL*8 G(LMO),S(LMO),TO(LMO),BYRHO(LMO),RHO(LMO),PO(LMO)
      REAL*8 UKJM(LMO,IM+1),UKM(LMO,4,2,2,IM,2:JM-1),OLJ(3,LMO,JM)
      LOGICAL, PARAMETER :: LDD = .FALSE.
      INTEGER I,J,K,L,IQ,JQ,LMIJ,KMUV,IM1,ITER,NSIGG,NSIGS,KBL,II
      REAL*8 CORIOL,UISTR,VISTR,U2rho,DELTAM,DELTAE,DELTASR,ANSTR
     *     ,RJ,RI,ZSCALE,HBL,HBLP,Ustar,BYSHC,B0,Bosol,R,R2,DTBYDZ2,DM
     *     ,RHOM,RHO1,Bo,DELTAS
      REAL*8 VOLGSP,ALPHAGSP,BETAGSP,TEMGSP,SHCGS,TEMGS

C**** Load UO,VO into UT,VT.  UO,VO will be updated, while UT,VT
C**** will be fixed during convection.
!$OMP PARALLEL DO PRIVATE(L)
      DO L=1,LMO
        UT(:,:,L) = UO(:,:,L)
        VT(:,:,L) = VO(:,:,L)
      END DO
!$OMP END PARALLEL DO
C****
C**** Outside loop over J
C**** Processes are checked and applied on every horizontal quarter box.
C****
!$OMP PARALLEL DO  PRIVATE(ANSTR,AKVM,AKVS,AKVG,ALPHADT, BYMML,BYMMLT,
!$OMP&  BYMML0,BYHWIDE,BYRHO,BYSHC,BO,BOSOL,BETADS,BYDZ2, CORIOL,DBLOC,
!$OMP&  DBSFC,DELTAE,DELTAM,DELTAS,DELTASR,DM,DTBYDZ,DTBYDZ2,DVSQ,
!$OMP&  FLG,FLS,G,G0ML0,G0ML,GHAT,GHATM,GHATG,GHATS,GZML, HBL,HBLP,
!$OMP&  HWIDE,I,II,IM1,ITER,IQ,J,JQ,K,KBL,KMUV,L,LMIJ,LMUV,MML,MML0,
!$OMP&  MMLT, NSIGG,NSIGS, PO, QPOLE, R,R2,RAVM,RAMV,RJ,RI,RITOP,
!$OMP&  RHO,RHO1,RHOM, S,S0ML0,S0ML,SBETA,SHSQ,SZML, TO,TALPHA,TXY,
#ifdef TRACERS_OCEAN
!$OMP&  N,TRML,TZML,TRML1,DELTATR,GHATT,FLT,NSIGT,
#endif
!$OMP&  UL,UL0, U2RHO,UISTR,USTAR, VISTR, ZGRID,ZSCALE)
      DO 790 J=2,JM
C**** coriolis parameter, defined at tracer point
      Coriol = 2d0*OMEGA*SINPO(J)
C**** initialise diagnostics saved over quarter boxes and longitude
      OLJ(1:3,1:LMO,J) = 0.

      QPOLE = J.EQ.JM
      IF (QPOLE)  THEN
C****
C**** Define parameters and currents at the North Pole
C****
      LMIJ=LMM(1,JM)
      IQ=1
      JQ=1
      KMUV=IM+1
      DO I=1,IM
        LMUV(I) = LMV(I,JM-1)
        RAVM(I) = 1d0/IM
        RAMV(I) = RAMVS(JM)
        UL0(1,I) = VT(I,JM-1,1)
        UL(1,I)  = VO1(I,JM-1)
        DO L=2,LMIJ
          UL0(L,I) = VT(I,JM-1,L)
          UL(L,I)  = UL0(L,I)
        END DO
      END DO
      LMUV(IM+1) = LMU(1,JM)
      RAVM(IM+1) = 1d0
      RAMV(IM+1) = 1d0
      UL0(1,IM+1) = UT(1,JM,1)
      MML0(1)     = MO(1,JM,1)*DXYPO(JM)
      BYMML0(1)   = 1d0/MML0(1)
      DTBYDZ(1)   = DTS/MO(1,JM,1)
      G0ML0(1,1,1)= G0M(1,JM,1)
      S0ML0(1,1,1)= S0M(1,JM,1)
      GZML(1,1,1) = GZMO(1,JM,1)
      SZML(1,1,1) = SZMO(1,JM,1)
      UL(1,IM+1)  = UO1(1,JM)
      MMLT(1)     = MO1(1,JM)*DXYPO(JM)
      BYMMLT(1)   = 1d0/MMLT(1)
      S0ML(1,1,1) = S0M1(1,JM)
      DO L=2,LMIJ
        UL0(L,IM+1) = UT(1,JM,L)
        MML0(L)     = MO(1,JM,L)*DXYPO(JM)
        BYMML0(L)   = 1d0/MML0(L)
        DTBYDZ(L)   = DTS/MO(1,JM,L)
        BYDZ2(L-1)  = 2d0/(MO(1,JM,L-1)+MO(1,JM,L))
        G0ML0(L,1,1)= G0M(1,JM,L)
        S0ML0(L,1,1)= S0M(1,JM,L)
        GZML(L,1,1) = GZMO(1,JM,L)
        SZML(L,1,1) = SZMO(1,JM,L)
        UL(L,IM+1)  = UL0(L,IM+1)
        MMLT(L)     = MML0(L)
        BYMMLT(L)   = BYMML0(L)
        S0ML(L,1,1) = S0ML0(L,1,1)
      END DO
      DO L=1,MIN(LSRPD,LMIJ)
        G0ML(L,1,1) = G0M1(1,JM,L)
      END DO
      DO L=LSRPD,LMIJ
        G0ML(L,1,1) = G0ML0(L,1,1)
      END DO
#ifdef TRACERS_OCEAN
      TRML1(:,1,1)=TRMO1(:,1,JM)
      DO L = 1,LMIJ
        TRML(L,:,1,1) = TRMO(1,JM,L,:)
        TZML(L,:,1,1) = TZMO(1,JM,L,:)
      END DO
#endif
C**** Calculate whole box turbulent stresses sqrt(|tau_0|/rho_0)
C**** except for RHO factor and sqrt: U2rho = Ustar**2 * rho
C**** DM[UV]A are defined on t-grid. DM[UV]I defined on u,v grid
C**** Note that rotational ice stresses are not calculated at the pole
C**** DMUA/I now defined over whole box, not just surface type
      UISTR = 0
      VISTR = 0
      DO I=1,IM
        UISTR = UISTR + DMVI(I,JM-1)*COSIC(I)
        VISTR = VISTR - DMVI(I,JM-1)*SINIC(I)
      END DO
      UISTR = UISTR/IM
      VISTR = VISTR/IM
      U2rho = SQRT(
     *     (DMUA(1,JM,1) + UISTR)**2+
     *     (DMVA(1,JM,1) + VISTR)**2)*BYDTS
C**** Calculate surface mass, salt and heat fluxes
      DELTAM = (MO(1,JM,1) -  MO1(1,JM))*BYDTS
      DELTAE = (G0ML0(1,1,1) - G0ML(1,1,1))*BYDXYPO(JM)*BYDTS
      DELTAS = (S0ML0(1,1,1) - S0ML(1,1,1))*BYDXYPO(JM)*BYDTS
      DELTASR= (SOLAR(1,1,JM)*(1d0-RSI(1,JM))+SOLAR(3,1,JM)*RSI(1,JM))
     *     *BYDTS               ! W/m^2
#ifdef TRACERS_OCEAN
      DELTATR(:)=(TRML(1,:,1,1)-TRML1(:,1,1))*BYDXYPO(JM)*BYDTS !kg/m2/s
#endif
      KPL(1,JM) = 1  ! Initialize mixed layer depth
      I=1
      GOTO 500
      END IF
C****
C**** Define parameters and currents away from the poles
C****
      IM1=IM
      I=1
  210 IF(LMM(I,J).LE.1)  GO TO 730
      LMIJ=LMM(I,J)
      KMUV=4
      LMUV(1) = LMU(IM1,J)
      LMUV(2) = LMU(I  ,J)
      LMUV(3) = LMV(I,J-1)
      LMUV(4) = LMV(I,J  )
      MML0(1)   = 2.5d-1*MO(I,J,1)*DXYPO(J)
      BYMML0(1) = 1d0/MML0(1)
      DTBYDZ(1) = DTS/MO(I,J,1)
      MMLT(1)   = 2.5d-1*MO1(I,J)*DXYPO(J)
      BYMMLT(1) = 1d0/MMLT(1)
      DO L=2,LMIJ
        MML0(L)   = 2.5d-1*MO(I,J,L)*DXYPO(J)
        BYMML0(L) = 1d0/MML0(L)
        DTBYDZ(L) = DTS/MO(I,J,L)
        BYDZ2(L-1)= 2d0/(MO(I,J,L)+MO(I,J,L-1))
        MMLT(L)   = MML0(L)
        BYMMLT(L) = BYMML0(L)
      END DO
C**** Calculate whole box turbulent stresses sqrt(|tau_0|/rho_0)
C**** except for RHO factor and sqrt: U2rho = Ustar**2 * rho
C**** DM[UV]A are defined on t-grid. DM[UV]I defined on u,v grid
C**** DMUA/I now defined over whole box, not just surface type
      UISTR = 0
      VISTR = 0
      IF (RSI(I,J).gt.0) THEN
        IF (LMU(I,J).gt.0)   UISTR = UISTR + DMUI(I,J)
        IF (LMU(IM1,J).gt.0) UISTR = UISTR + DMUI(IM1,J)
        ANSTR=1.-SIGN(0.25,LMU(I,J)-0.5)-SIGN(0.25,LMU(IM1,J)-0.5)
        UISTR = UISTR*ANSTR
        IF (LMV(I,J).gt.0)   VISTR = VISTR + DMVI(I,J)
        IF (LMV(I,J-1).gt.0) VISTR = VISTR + DMVI(I,J-1)
        ANSTR=1.-SIGN(0.25,LMV(I,J)-0.5)-SIGN(0.25,LMV(I,J-1)-0.5)
        VISTR = VISTR*ANSTR
      END IF
      U2rho = SQRT((DMUA(I,J,1) + UISTR)**2 +
     *             (DMVA(I,J,1) + VISTR)**2)*BYDTS
C**** Calculate surface mass flux and Solar forcing
      DELTAM = (MO(I,J,1) -  MO1(I,J))*BYDTS ! kg/m^2 s
      DELTASR = (SOLAR(1,I,J)*(1d0-RSI(I,J))+SOLAR(3,I,J)*RSI(I,J))
     *     *BYDTS               ! W/m^2
      KPL(I,J) = 1  ! Initialize mixed layer depth

C**** Loop over quarter boxes
      JQ=1
  230 RJ = BETA*(2d0*JQ-3d0)
      RAVM(3) = (1.25d0-5d-1*JQ)
      RAMV(3) = (1.25d0-5d-1*JQ)*RAMVS(J)*5d-1
      RAVM(4) = (5d-1*JQ-2.5d-1)
      RAMV(4) = (5d-1*JQ-2.5d-1)*RAMVN(J)*5d-1
      IQ=1
  240 RI = BETA*(2d0*IQ-3d0)
      RAVM(1) = (1.25d0 - 5d-1*IQ)
      RAMV(1) = (1.25d0 - 5d-1*IQ)*2.5d-1
      RAVM(2) = (5d-1*IQ-2.5d-1)
      RAMV(2) = (5d-1*IQ-2.5d-1)*2.5d-1
      UL0(1,1) = UT(IM1,J,1)
      UL0(1,2) = UT(I  ,J,1)
      UL0(1,3) = VT(I,J-1,1)
      UL0(1,4) = VT(I,J  ,1)
      G0ML0(1,IQ,JQ)=2.5d-1*(G0M(I,J,1)+RI*GXMO(I,J,1)+RJ*GYMO(I,J,1))
      S0ML0(1,IQ,JQ)=2.5d-1*(S0M(I,J,1)+RI*SXMO(I,J,1)+RJ*SYMO(I,J,1))
      GZML(1,IQ,JQ)=2.5d-1* GZMO(I,J,1)
      SZML(1,IQ,JQ)=2.5d-1* SZMO(I,J,1)
      UL(1,1) = UO1(IM1,J)
      UL(1,2) = UO1(I  ,J)
      UL(1,3) = VO1(I,J-1)
      UL(1,4) = VO1(I,J  )
      S0ML(1,IQ,JQ)=2.5d-1*(S0M1(I,J) + RI*SXM1(I,J) + RJ*SYM1(I,J))
      if (abs(SZML(1,IQ,JQ))>S0ML(1,IQ,JQ))
     \     SZML(1,IQ,JQ) = sign(S0ML(1,IQ,JQ),SZML(1,IQ,JQ))
      DO 250 L=2,LMIJ
      UL0(L,1) = UT(IM1,J,L)
      UL0(L,2) = UT(I  ,J,L)
      UL0(L,3) = VT(I,J-1,L)
      UL0(L,4) = VT(I,J  ,L)
      G0ML0(L,IQ,JQ)=2.5d-1*(G0M(I,J,L)+RI*GXMO(I,J,L)+RJ*GYMO(I,J,L))
      S0ML0(L,IQ,JQ)=2.5d-1*(S0M(I,J,L)+RI*SXMO(I,J,L)+RJ*SYMO(I,J,L))
      GZML(L,IQ,JQ)=2.5d-1* GZMO(I,J,L)
      SZML(L,IQ,JQ)=2.5d-1* SZMO(I,J,L)
      UL(L,1) = UL0(L,1)
      UL(L,2) = UL0(L,2)
      UL(L,3) = UL0(L,3)
      UL(L,4) = UL0(L,4)
      S0ML(L,IQ,JQ) = S0ML0(L,IQ,JQ)
      if (abs(SZML(L,IQ,JQ))>S0ML(L,IQ,JQ))
     \     SZML(L,IQ,JQ) = sign(S0ML(L,IQ,JQ),SZML(L,IQ,JQ))
  250 CONTINUE
      G0ML(1,IQ,JQ)=2.5d-1*(G0M1(I,J,1)+ RI*GXM1(I,J) + RJ*GYM1(I,J))
      DO L=2,MIN(LSRPD,LMIJ)
        G0ML(L,IQ,JQ)=2.5d-1*(G0M1(I,J,L)+RI*GXMO(I,J,L)+RJ*GYMO(I,J,L))
      END DO
      DO L=LSRPD,LMIJ
        G0ML(L,IQ,JQ) = G0ML0(L,IQ,JQ)
      END DO
#ifdef TRACERS_OCEAN
      DO N=1,NTM
        TRML1(N,IQ,JQ) = 2.5d-1*(TRMO1(N,I,J)+RI*TXMO1(N,I,J)
     *       + RJ*TYMO1(N,I,J))
        DO L=1,LMIJ
          TRML(L,N,IQ,JQ) = 2.5d-1*(TRMO(I,J,L,N)+RI*TXMO(I,J,L,N)
     *         + RJ*TYMO(I,J,L,N))
          TZML(L,N,IQ,JQ) = 2.5d-1* TZMO(I,J,L,N)
C****
          if (t_qlimit(N)) then
            TRML(L,N,IQ,JQ) = MAX(0d0,TRML(L,N,IQ,JQ))
            if (abs(TZML(L,N,IQ,JQ))>TRML(L,N,IQ,JQ))
     \           TZML(L,N,IQ,JQ) = sign(TRML(L,N,IQ,JQ),TZML(L,N,IQ,JQ))
          end if
        END DO
      END DO
#endif

C**** Calculate surface heat and salt flux: dE(W/m^2);dS,dTr(kg/m^2/s)
      DELTAE=4d0*(G0ML0(1,IQ,JQ)-G0ML(1,IQ,JQ))*BYDXYPO(J)*BYDTS 
      DELTAS=4d0*(S0ML0(1,IQ,JQ)-S0ML(1,IQ,JQ))*BYDXYPO(J)*BYDTS 
#ifdef TRACERS_OCEAN
      DELTATR(:)=4d0*(TRML(1,:,IQ,JQ)-TRML1(:,IQ,JQ))*BYDXYPO(J)*BYDTS 
#endif
C****
C**** Vertical mixing dependent on KPP boundary layer scheme
C****
  500 IF(LMIJ.LE.1)  GO TO 725
C**** Do not recalculate quantities good for all quarter boxes
      IF(IQ+JQ.EQ.2)  THEN

C**** Z0 = OGEOZ/GRAV+HOCEAN Ocean depth (m)
C****      scale depths using Z0/ZE(LMIJ)

      ZSCALE = (OGEOZ(I,J) + GRAV*HOCEAN(I,J))/(ZE(LMIJ)*GRAV)
c       OIJ(I,J,IJ_OGEOZ)=OIJ(I,J,IJ_OGEOZ)+OGEOZ(I,J)
      zgrid(0) = epsln
      hwide(0) = epsln
      byhwide(0) = 0.
      PO(1) = 5d-1*MO(I,J,1)*GRAV
      DO L=1,LMO-1
         zgrid(L) = -0.5*ZSCALE*(ZE(L-1) + ZE(L)) ! tracer level
         hwide(L) = zgrid(L-1) - zgrid(L)       ! tracer level
         byhwide(L) = 1d0/hwide(L)
         PO(L+1) = PO(L) + 5d-1*GRAV*(MO(I,J,L)+MO(I,J,L+1))
      END DO
      zgrid(LMO)   = -0.5*ZSCALE*(ZE(LMO-1) + ZE(LMO)) ! tracer level
      hwide(LMO)   = zgrid(LMO-1) - zgrid(LMO) ! tracer level
      byhwide(LMO) = 1d0/hwide(LMO)
      zgrid(LMO+1) = -ZE(LMO)*ZSCALE
      hwide(LMO+1) = epsln
      byhwide(LMO+1) = 0.

      END IF

C**** If swfrac is not used add solar from next two levels to Bo
C      DELTAE = DELTAE + FSR(2)*DELTASR
C**** else remove solar from Bo
      DELTAE = DELTAE - (1d0-FSR(2))*DELTASR

C**** Start iteration for diffusivities
      MML = MMLT   ! start with pre-source mass
      BYMML = BYMMLT
      ITER = 0
      HBL = 0
  510 ITER =ITER + 1
      HBLP = HBL
      IF (ITER.eq.2) THEN
        MML = MML0      ! continue with post-source mass
        BYMML = BYMML0
      END IF
C**** Initiallise fields
      Shsq=0.
      dVsq=0.
      dbloc=0.
      dbsfc=0.
      Ritop=0.

C**** velocity shears
      DO L=1,LMIJ-1
        DO K=1,KMUV
          Shsq(L) = Shsq(L) + RAVM(K)*(UL(L,K)-UL(L+1,K))**2 !interface
          dVsq(L) = dVsq(L) + RAVM(K)*(UL(1,K)-UL(L  ,K))**2 !tracer pts
         END DO
      END DO
      DO K=1,KMUV
        dVsq(LMIJ)=dVsq(LMIJ)+RAVM(K)*(UL(1,K)-UL(LMIJ,K))**2
      END DO

C****    density related quantities
C****
C****    density of surface layer                                (kg/m3)
C****          rho   = rho{G(1),S(1),PO(1)}
C****    local buoyancy gradient at km interfaces:                (m/s2)
C****          dbloc = g/rho{k+1,k+1} * [ rho{k,k+1}-rho{k+1,k+1} ]
C****    buoyancy difference with respect to "zref", the surface  (m/s2)
C****          dbsfc = g/rho{k,k} * [ rho{k,k}-rho{1,k} ] (CORRECT)
C****    thermal expansion coefficient without 1/rho factor    (kg/m3/C)
C****          talpha= d(rho{k,k})/d(T(k))
C****    salt expansion coefficient without 1/rho factor     (kg/m3/PSU)
C****          sbeta = d(rho{k,k})/d(S(k))

      DO L=1,LMIJ
         G(L)     = G0ML(L,IQ,JQ)*BYMML(L)
         S(L)     = S0ML(L,IQ,JQ)*BYMML(L)
         BYRHO(L) = VOLGSP(G(L),S(L),PO(L))
         RHO(L)   = 1d0/BYRHO(L)
         IF (L.ge.2) THEN
           RHOM = 1d0/VOLGSP(G(L-1),S(L-1),PO(L)) ! L-1 density wrt P_L
           RHO1 = 1d0/VOLGSP(G(1),S(1),PO(L)) ! surf density wrt P_L
           dbsfc(L)   = GRAV*(1d0-RHO1*BYRHO(L))
           dbloc(L-1) = GRAV*(1d0-RHOM*BYRHO(L))
C**** numerator of bulk richardson number on grid levels
           Ritop(L) = (zgrid(1)-zgrid(L)) * dbsfc(L)
         END IF
      END DO
      talpha(1) =  ALPHAGSP(G(1),S(1),PO(1)) ! >0
      sbeta(1)  =   BETAGSP(G(1),S(1),PO(1))

C**** surface wind turbulent friction speed (m/s) = sqrt( |tau_0|/rho )
      Ustar = SQRT(U2rho*BYRHO(1))

C**** surface buoyancy forcing turbulent and solar (m^2/s^3)
C****    Bo    = -g*(talpha * wtsfc - sbeta * wssfc)/rho
C****    Bosol = -g*(talpha * wtsol                )/rho
C**** Bo includes all buoyancy and heat forcing

      BYSHC = 1d0/SHCGS(G(1),S(1))

      Bo    = - GRAV*BYRHO(1)**2 *(
     *       sbeta(1)*DELTAS*1d3 + talpha(1)*BYSHC*DELTAE -
     *     ( sbeta(1)*S(1)*1d3   + talpha(1)*BYSHC*G(1) )*DELTAM)
      Bosol = - GRAV*BYRHO(1)**2 * talpha(1)*BYSHC*DELTASR

C**** Double diffusive option (ddmix) needs alphaDT,betaDS
C**** alphaDT  = mean -talpha * delta(temp.) at interfaces (kg/m3)
C**** betaDS   = mean sbeta  * delta(salt)  at interfaces (kg/m3)

      if (LDD) then
        TO(1)     =    TEMGSP(G(1),S(1),PO(1))
        do L=2,LMIJ
          TO(L)     =   TEMGSP(G(L),S(L),PO(L))
          talpha(L)= ALPHAGSP(G(L),S(L),PO(L)) ! >0
          sbeta(L) =  BETAGSP(G(L),S(L),PO(L))
          alphaDT(L-1) = - 5d-1 * (talpha(L-1) + talpha(L))
     $         * (TO(L-1) - TO(L))
          betaDS (L-1) = 5d-1 * (sbeta (L-1) + sbeta(L))
     $         * (S(L-1) - S(L))*1d3
        end do
      end if

      CALL KPPMIX(LDD,ZE,zgrid,hwide,LMIJ,Shsq,dVsq,Ustar,Bo
     *     ,Bosol ,alphaDT,betaDS,dbloc,Ritop,Coriol,byhwide,
     *     AKVM,AKVS,AKVG,GHAT,HBL,KBL)

C**** Calculate non-local transport terms for each scalar
C****        ghat[sg] = kv * ghat * <w[sg]0>   (J,kg)
C**** Correct units for diffusivities (m^2/s) => (kg^2/m^4 s)
C****                            ghat (s/m^2) => (s m^4/kg^2)
      DO L=1,LMIJ-1
         R = 5d-1*(RHO(L)+RHO(L+1))
         R2 = R**2
         GHATM(L) = 0.          ! no non-local momentum transport
C**** GHAT terms must be zero for consistency with OSOURC
         GHATG(L) = 0. !AKVG(L)*GHAT(L)*DELTAE*DXYPO(J)
         GHATS(L) = 0. !AKVS(L)*GHAT(L)*DELTAM*S(1)*DXYPO(J) CHECK!
         AKVM(L) = AKVM(L)*R2
         AKVG(L) = AKVG(L)*R2
         AKVS(L) = AKVS(L)*R2
      END DO

C**** For each field (U,G,S + TRACERS) call OVDIFF
C**** Momentum
      DO K=1,KMUV
        CALL OVDIFF(UL(1,K),AKVM(1),GHATM,DTBYDZ,BYDZ2
     *       ,LMIJ,UL0(1,K))
      END DO
C**** Enthalpy
      CALL OVDIFFS(G0ML(1,IQ,JQ),AKVG(1),GHATG,DTBYDZ,BYDZ2,DTS
     *     ,LMIJ,G0ML0(1,IQ,JQ),FLG)
C**** Salinity
      CALL OVDIFFS(S0ML(1,IQ,JQ),AKVS(1),GHATS,DTBYDZ,BYDZ2,DTS
     *     ,LMIJ,S0ML0(1,IQ,JQ),FLS)
      IF ((ITER.eq.1  .or. ABS(HBLP-HBL).gt.(ZE(KBL)-ZE(KBL-1))*0.25)
     *     .and. ITER.lt.4) GO TO 510
#ifdef TRACERS_OCEAN
C**** Tracers are diffused after iteration and follow salinity
      DO L=1,LMIJ
        GHATT(L,:)=0.         !AKVS(L)*GHAT(L)*DELTATR(:)*DXYP(J)*0.25
      END DO
      DO N=1,NTM
        CALL OVDIFFS(TRML(1,N,IQ,JQ),AKVS(1),GHATT(1,N),DTBYDZ,BYDZ2
     *       ,DTS,LMIJ,TRML(1,N,IQ,JQ),FLT(1,N))
      END DO
#endif
C**** Gradients of scalars
C**** Implicitly apply interpolated KV to linear profile
C**** Surface tracer + mass fluxes included (no solar flux)
C**** Note that FL[GS] are upward fluxes.
      DTBYDZ2 = 12d0*DTBYDZ(1)**2*BYDTS
      DM = DELTAM*DTBYDZ(1)
      GZML(1,IQ,JQ)=(GZML(1,IQ,JQ)+3d0*(FLG(1)
     *     +0.25*DTS*DELTAE*DXYPO(J)))/(1d0+DTBYDZ2*AKVG(1))
      SZML(1,IQ,JQ)=(SZML(1,IQ,JQ)+3d0*(FLS(1)-DELTAS*(1-DM)*0.25
     *     *DXYPO(J)*DTS+DM*2.5d-1*(S0M1(I,J) + RI*SXM1(I,J) + RJ*SYM1(I
     *     ,J))))/(1d0+DTBYDZ2*AKVS(1))
#ifdef TRACERS_OCEAN
      TZML(1,:,IQ,JQ)=(TZML(1,:,IQ,JQ)+3d0*(FLT(1,:)-
     *     DELTATR(:)*(1-DM)*0.25*DXYPO(J)*DTS+DM*2.5d-1*(TRMO1(:,I,J)
     *     +RI*TXMO1(:,I,J)+ RJ*TYMO1(:,I,J))))/(1d0+DTBYDZ2*AKVS(1))
#endif
      DO L=2,LMIJ-1
        DTBYDZ2 = 6d0*DTBYDZ(L)**2*BYDTS
        GZML(L,IQ,JQ)=(GZML(L,IQ,JQ)+3d0*(FLG(L-1)+FLG(L)))
     *       /(1d0+DTBYDZ2*(AKVG(L-1)+AKVG(L)))
        SZML(L,IQ,JQ)=(SZML(L,IQ,JQ)+3d0*(FLS(L-1)+FLS(L)))
     *       /(1d0+DTBYDZ2*(AKVS(L-1)+AKVS(L)))
#ifdef TRACERS_OCEAN
        TZML(L,:,IQ,JQ)=(TZML(L,:,IQ,JQ)+3d0*(FLT(L-1,:)+
     *       FLT(L,:))) /(1d0+DTBYDZ2*(AKVS(L-1)+AKVS(L)))
#endif
      END DO
      DTBYDZ2 = 12d0*DTBYDZ(LMIJ)**2*BYDTS
      GZML(LMIJ,IQ,JQ)=(GZML(LMIJ,IQ,JQ)+3d0*FLG(LMIJ-1))
     *     /(1d0+DTBYDZ2*AKVG(LMIJ-1))
      SZML(LMIJ,IQ,JQ)=(SZML(LMIJ,IQ,JQ)+3d0*FLS(LMIJ-1))
     *     /(1d0+DTBYDZ2*AKVS(LMIJ-1))
#ifdef TRACERS_OCEAN
      TZML(LMIJ,:,IQ,JQ)=(TZML(LMIJ,:,IQ,JQ)+3d0*FLT(LMIJ-1,:))
     *     /(1d0+DTBYDZ2*AKVS(LMIJ-1))
#endif

C**** Diagnostics for non-local transport and vertical diffusion
       DO L=1,LMIJ-1
         OIJL(I,J,L,IJL_KVM) = OIJL(I,J,L,IJL_KVM) + AKVM(L)
         OIJL(I,J,L,IJL_KVG) = OIJL(I,J,L,IJL_KVG) + AKVG(L)
         OIJL(I,J,L,IJL_WGFL)= OIJL(I,J,L,IJL_WGFL) + FLG(L) ! heat flux
         OIJL(I,J,L,IJL_WSFL)= OIJL(I,J,L,IJL_WSFL) + FLS(L) ! salt flux
c         OIJL(I,J,L,IJL_KVS) = OIJL(I,J,L,IJL_KVS) + AKVS(L)
c         OIJL(I,J,L,IJL_KVGG) = OIJL(I,J,L,IJL_KVGG) + AKVG(L)*GHATG(L)
c         OIJL(I,J,L,IJL_KVSG) = OIJL(I,J,L,IJL_KVSG) + AKVS(L)*GHATS(L)
       END DO
C**** Also set vertical diagnostics
       DO L=1,LMIJ
CCC         OL(L,L_RHO) = OL(L,L_RHO) + RHO(L)
CCC         OL(L,L_TEMP)= OL(L,L_TEMP)+ TEMGS(G(L),S(L))
CCC         OL(L,L_SALT)= OL(L,L_SALT)+ S(L)
         OLJ(1,L,J)= OLJ(1,L,J) + RHO(L) ! L_RHO
         OLJ(2,L,J)= OLJ(2,L,J) + TEMGS(G(L),S(L)) ! L_TEMP
         OLJ(3,L,J)= OLJ(3,L,J) + S(L) ! L_SALT
       END DO
C**** Set diagnostics
       OIJ(I,J,IJ_HBL) = OIJ(I,J,IJ_HBL) + HBL ! boundary layer depth
       OIJ(I,J,IJ_BO) = OIJ(I,J,IJ_BO) + Bo ! surface buoyancy forcing
       OIJ(I,J,IJ_BOSOL) = OIJ(I,J,IJ_BOSOL) + Bosol ! solar buoy frcg
       OIJ(I,J,IJ_USTAR) = OIJ(I,J,IJ_USTAR) + Ustar ! turb fric speed
       IF(KBL.gt.KPL(I,J)) KPL(I,J)=KBL ! save max. mixed layer depth
C****
C**** Update current prognostic variables
C****
CCC      IF (QPOLE) THEN
CCC        DO II=1,IM
CCC          VO(II,JM-1,1:LMUV(II))=VO(II,JM-1,1:LMUV(II))+RAMV(II)*
CCC     *         (UL(1:LMUV(II),II)-VT(II,JM-1,1:LMUV(II)))
CCC        END DO
CCC        UO(1,JM,1:LMIJ)=UO(1,JM,1:LMIJ) + RAMV(IM+1)*
CCC     *       (UL(1:LMIJ,IM+1)-UT(1,JM,1:LMIJ))
CCC      ELSE
CCC        UO(IM1,J,1:LMUV(1))=UO(IM1,J,1:LMUV(1)) + RAMV(1)*
CCC     *       (UL(1:LMUV(1),1)-UT(IM1,J,1:LMUV(1)))
CCC        UO(I  ,J,1:LMUV(2))=UO(I  ,J,1:LMUV(2)) + RAMV(2)*
CCC     *       (UL(1:LMUV(2),2)-UT(I  ,J,1:LMUV(2)))
CCC        VO(I,J-1,1:LMUV(3))=VO(I,J-1,1:LMUV(3)) + RAMV(3)*
CCC     *       (UL(1:LMUV(3),3)-VT(I,J-1,1:LMUV(3)))
CCC        VO(I,J  ,1:LMUV(4))=VO(I,J  ,1:LMUV(4)) + RAMV(4)*
CCC     *       (UL(1:LMUV(4),4)-VT(I,J  ,1:LMUV(4)))
CCC      END IF

      IF(QPOLE)  THEN
         DO II=1,IM
           UKJM(1:LMUV(II),II) = RAMV(II)*(UL(1:LMUV(II),II)-
     *          VT(II,JM-1,1:LMUV(II)))
         END DO
         UKJM(1:LMIJ,IM+1)= RAMV(IM+1)*(UL(1:LMIJ,IM+1)-
     *        UT(1,JM,1:LMIJ))
      ELSE
        UKM(1:LMUV(1),1,IQ,JQ,I,J) = RAMV(1)*(UL(1:LMUV(1),1)
     *          -UT(IM1,J,1:LMUV(1)))
        UKM(1:LMUV(2),2,IQ,JQ,I,J) = RAMV(2)*(UL(1:LMUV(2),2)
     *          -UT(I  ,J,1:LMUV(2)))
        UKM(1:LMUV(3),3,IQ,JQ,I,J) = RAMV(3)*(UL(1:LMUV(3),3)
     *          -VT(I,J-1,1:LMUV(3)))
        UKM(1:LMUV(4),4,IQ,JQ,I,J) = RAMV(4)*(UL(1:LMUV(4),4)
     *          -VT(I,J  ,1:LMUV(4)))
      END IF

  725 IF(QPOLE)  GO TO 750
C**** End of quarter box loops
      IQ=IQ+1
      IF(IQ.LE.2)  GO TO 240
      JQ=JQ+1
      IF(JQ.LE.2)  GO TO 230
C****
C**** Recreate main prognostic variables of potential enthalpy and
C**** salt from those on the quarter boxes
C****
      DO L=1,LMIJ
      G0M(I,J,L)= G0ML(L,2,2)+G0ML(L,2,1)+G0ML(L,1,2)+G0ML(L,1,1)
      NSIGG = EXPONENT(G0M(I,J,L)) -2 - 42
      GXMO(I,J,L)=(G0ML(L,2,2)+G0ML(L,2,1)-G0ML(L,1,2)-G0ML(L,1,1))
     *     *BYBETA
      GYMO(I,J,L)=(G0ML(L,2,2)-G0ML(L,2,1)+G0ML(L,1,2)-G0ML(L,1,1))
     *     *BYBETA
      CALL REDUCE_FIG(NSIGG,GXMO(I,J,L))
      CALL REDUCE_FIG(NSIGG,GYMO(I,J,L))
      GZMO(I,J,L)= GZML(L,2,2)+GZML(L,2,1)+GZML(L,1,2)+GZML(L,1,1)
      S0M(I,J,L)= S0ML(L,2,2)+S0ML(L,2,1)+S0ML(L,1,2)+S0ML(L,1,1)
      NSIGS = EXPONENT(S0M(I,J,L)) - 2 - 42 + 4
      SXMO(I,J,L)=(S0ML(L,2,2)+S0ML(L,2,1)-S0ML(L,1,2)-S0ML(L,1,1))
     *     *BYBETA
      SYMO(I,J,L)=(S0ML(L,2,2)-S0ML(L,2,1)+S0ML(L,1,2)-S0ML(L,1,1))
     *     *BYBETA
      CALL REDUCE_FIG(NSIGS,SXMO(I,J,L))
      CALL REDUCE_FIG(NSIGS,SYMO(I,J,L))
      SZMO(I,J,L)= SZML(L,2,2)+SZML(L,2,1)+SZML(L,1,2)+SZML(L,1,1)
C**** limit salinity gradients
      TXY = abs(SXMO(I,J,L)) + abs(SYMO(I,J,L))
      if ( TXY > S0M(I,J,L) ) then
        SXMO(I,J,L)=SXMO(I,J,L)*(S0M(I,J,L)/(TXY + tiny(TXY)))
        SYMO(I,J,L)=SYMO(I,J,L)*(S0M(I,J,L)/(TXY + tiny(TXY)))
      end if
      if ( abs(SZMO(I,J,L)) > S0M(I,J,L) )
     *     SZMO(I,J,L) = sign(S0M(I,J,L),SZMO(I,J,L)+0d0)
C****
      END DO
#ifdef TRACERS_OCEAN
      DO N = 1,NTM
        DO L = 1,LMIJ
        TRMO(I,J,L,N) = TRML(L,N,2,2) + TRML(L,N,2,1) +
     *                  TRML(L,N,1,2) + TRML(L,N,1,1)
        NSIGT = EXPONENT(TRMO(I,J,L,N)) - 2 - 42
        TXMO(I,J,L,N) =(TRML(L,N,2,2) + TRML(L,N,2,1) -
     *                  TRML(L,N,1,2) - TRML(L,N,1,1))*BYBETA
        TYMO(I,J,L,N) =(TRML(L,N,2,2) - TRML(L,N,2,1) +
     *                  TRML(L,N,1,2) - TRML(L,N,1,1))*BYBETA
        CALL REDUCE_FIG(NSIGT,TXMO(I,J,L,N))
        CALL REDUCE_FIG(NSIGT,TYMO(I,J,L,N))
        TZMO(I,J,L,N) = TZML(L,N,2,2) + TZML(L,N,2,1) +
     *                  TZML(L,N,1,2) + TZML(L,N,1,1)
C****
        if (t_qlimit(n)) then   ! limit gradients
          TXY = abs(TXMO(I,J,L,N)) + abs(TYMO(I,J,L,N))
          if ( TXY > TRMO(I,J,L,N) ) then
            TXMO(I,J,L,N) = TXMO(I,J,L,N)
     *           *( TRMO(I,J,L,N)/(TXY + tiny(TXY)) )
            TYMO(I,J,L,N) = TYMO(I,J,L,N)
     *           *( TRMO(I,J,L,N)/(TXY + tiny(TXY)) )
          end if
          if ( abs(TZMO(I,J,L,N)) > TRMO(I,J,L,N) )
     *         TZMO(I,J,L,N) = sign(TRMO(I,J,L,N),TZMO(I,J,L,N)+0d0)
        end if
C****
      END DO
      END DO
#endif
C**** End of I loop
  730 IM1=I
      I=I+1
      IF(I.LE.IM)  GO TO 210
      GO TO 790
C****
C**** Load the prognostic variables of potential enthalpy and
C**** salt from the column arrays at the poles
C****
  750 DO L=1,LMIJ
      G0M(1,JM,L) = G0ML(L,1,1)
      GZMO(1,JM,L) = GZML(L,1,1)
      S0M(1,JM,L) = S0ML(L,1,1)
      SZMO(1,JM,L) = SZML(L,1,1)
#ifdef TRACERS_OCEAN
      TRMO(1,JM,L,:) = TRML(L,:,1,1)
      TZMO(1,JM,L,:) = TZML(L,:,1,1)
#endif
      END DO
C**** End of outside J loop
  790 CONTINUE
!$OMP END PARALLEL DO

C**** Update velocities outside parallel region

C**** North pole
      DO I=1,IM
        VO(I,JM-1,1:LMV(I,JM-1))=VO(I,JM-1,1:LMV(I,JM-1))+
     *       UKJM(1:LMV(I,JM-1),I)
      END DO
      UO(1,JM,1:LMU(1,JM))=UO(1,JM,1:LMU(1,JM)) + UKJM(1:LMU(1,JM),IM+1)
C**** Everywhere else
      DO J=2,JM-1
        IM1=IM
        DO I=1,IM
          DO IQ=1,2
          DO JQ=1,2
            UO(IM1,J,1:LMU(IM1,J))=UO(IM1,J,1:LMU(IM1,J)) +
     *           UKM(1:LMU(IM1,J),1,IQ,JQ,I,J)
            UO(I  ,J,1:LMU(I  ,J))=UO(I  ,J,1:LMU(I  ,J)) +
     *           UKM(1:LMU(I  ,J),2,IQ,JQ,I,J)
            VO(I,J-1,1:LMV(I,J-1))=VO(I,J-1,1:LMV(I,J-1)) +
     *           UKM(1:LMV(I,J-1),3,IQ,JQ,I,J)
            VO(I,J  ,1:LMV(I,J  ))=VO(I,J  ,1:LMV(I,J  )) +
     *           UKM(1:LMV(I,J  ),4,IQ,JQ,I,J)
          END DO
          END DO
          IM1=I
        END DO
      END DO
C**** sum global mean diagnostics
      DO J=2,JM
        DO L=1,LMO
          OL(L,L_RHO) = OL(L,L_RHO) + OLJ(1,L,J)
          OL(L,L_TEMP)= OL(L,L_TEMP)+ OLJ(2,L,J)
          OL(L,L_SALT)= OL(L,L_SALT)+ OLJ(3,L,J)
        END DO
      END DO
C****
      RETURN
      END SUBROUTINE OCONV

      SUBROUTINE STCONV
!@sum  STCONV uses vertical diffusion coefficients from KPP schmeme
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0
      USE CONSTANT, only : grav,omega
#ifdef TRACERS_OCEAN
      USE TRACER_COM, only : t_qlimit
#endif
      USE OCEAN,only : lmo,dts,ze,sinpo
      USE STRAITS, only : must,mmst,g0mst,gzmst,gxmst,s0mst,szmst,sxmst
     *     ,lmst,nmst,dist,wist,jst
#ifdef TRACERS_OCEAN
     *     ,trmst,txmst,tzmst,ntm
#endif
      USE ODIAG, only : olnst,ln_kvm,ln_kvg,ln_wgfl,ln_wsfl
      IMPLICIT NONE
      REAL*8 MMLT,MML0
      REAL*8, DIMENSION(LMO,2) :: UL,G0ML,S0ML,GZML,SZML
      REAL*8, DIMENSION(LMO) :: MML,BYMML,DTBYDZ,BYDZ2,UL0,G0ML0,S0ML0
      REAL*8, DIMENSION(0:LMO+1) :: AKVM,AKVG,AKVS
      REAL*8, DIMENSION(LMO) :: G,S,TO,BYRHO,RHO,PO,GHAT,FLG,FLS
#ifdef TRACERS_OCEAN
      REAL*8 TRML(LMO,NTM,2),TZML(LMO,NTM,2),FLT(LMO,NTM)
      INTEGER ITR,NSIGT
#endif
C**** CONV parameters: BETA controls degree of convection (default 0.5).
      REAL*8, PARAMETER :: BETA=5d-1,BYBETA=1d0/BETA
      LOGICAL*4, PARAMETER :: LDD=.FALSE.
      REAL*8, SAVE :: zgrid(0:LMO+1),hwide(0:LMO+1),byhwide(0:LMO+1)
      REAL*8 Shsq(LMO),dVsq(LMO)
     *     ,talpha(LMO),sbeta(LMO),dbloc(LMO),dbsfc(LMO),Ritop(LMO),
     *     alphaDT(LMO),betaDS(LMO)
      REAL*8, PARAMETER :: epsln=1d-20
      REAL*8, SAVE :: Bo,Bosol,bydts,ustar
      REAL*8 U2rho,RI,Coriol,HBL,HBLP,RHOM,RHO1,R,R2,DTBYDZ2
      REAL*8 VOLGSP,ALPHAGSP,BETAGSP,TEMGSP,SHCGS
      INTEGER, SAVE :: IFIRST = 1
      INTEGER I,L,N,LMIJ,IQ,ITER,NSIGG,NSIGS,KBL

      IF (IFIRST.eq.1) THEN
        IFIRST=0
        zgrid(0) = epsln
        hwide(0) = epsln
        byhwide(0) = 0.
        DO L=1,LMO
          zgrid(L) = -0.5*(ZE(L-1) + ZE(L)) ! tracer level
          hwide(L) = zgrid(L-1) - zgrid(L) ! tracer level
          byhwide(L) = 1d0/hwide(L)
        END DO
        zgrid(LMO+1) = -ZE(LMO)
        hwide(LMO+1) = epsln
        byhwide(LMO+1) = 0.

        Ustar = 0.   ! No wind/ice stress in straits
        Bo = 0       ! No buoyancy forcing
        Bosol = 0
        BYDTS = 1d0/DTS
      END IF
C****
C**** Outside loop over straits
C****
      DO 790 N=1,NMST
C****
C**** Define parameters and currents
C****
      LMIJ=LMST(N)
      PO(1) = 5d-1*GRAV*MMST(1,N)/(DIST(N)*WIST(N))
      DO 220 L=1,LMIJ-1
      PO(L+1) = PO(L  ) + 5d-1*GRAV*(MMST(L,N)+MMST(L+1,N))
     *       /(DIST(N)*WIST(N))
      DTBYDZ(L) = DTS*DIST(N)*WIST(N)/MMST(L,N)
      BYDZ2(L)  = 2d0*DIST(N)*WIST(N)/(MMST(L+1,N)+MMST(L,N))
      MML(L)    = 5d-1*MMST(L,N)
  220 BYMML(L)  =  1d0/MML(L)
      MML(LMIJ) = 5d-1*MMST(LMIJ,N)
      BYMML(LMIJ) = 1d0/MML(LMIJ)
      DTBYDZ(LMIJ) = DTS*DIST(N)*WIST(N)/MMST(LMIJ,N)
C**** Loop over half boxes
      IQ=1
  240 RI = BETA*(2d0*IQ-3d0)
      DO L=1,LMIJ
        UL(L,IQ) = 5d-1*  MUST(L,N)*DIST(N)*BYMML(L)
      G0ML(L,IQ) = 5d-1*(G0MST(L,N) + RI*GXMST(L,N))
      S0ML(L,IQ) = 5d-1*(S0MST(L,N) + RI*SXMST(L,N))
      GZML(L,IQ) = 5d-1* GZMST(L,N)
      SZML(L,IQ) = 5d-1* SZMST(L,N)
#ifdef TRACERS_OCEAN
      TRML(L,:,IQ)=5d-1*(TRMST(L,N,:) + RI*TXMST(L,N,:))
      TZML(L,:,IQ)=5d-1* TZMST(L,N,:)
#endif
      END DO
C****
C**** Vertical mixing derived from KPP scheme
C****
  500 IF(LMIJ.LE.1)  GO TO 700

C**** Coriolis parameter, defined at strait end
      Coriol = 2d0*OMEGA*SINPO(JST(N,IQ))
C**** Save original values
        UL0 = UL(1:LMO,IQ)
      G0ML0 = G0ML(1:LMO,IQ)
      S0ML0 = S0ML(1:LMO,IQ)
C**** Start iteration for diffusivities
      ITER = 0
      HBL = 0
  510 ITER =ITER + 1
      HBLP = HBL
C**** Initiallise fields
      Shsq=0.
      dVsq=0.
      dbloc=0.
      dbsfc=0.
      Ritop=0.
      alphaDT=0. ; betaDS=0.

C**** velocity shears

      DO L=1,LMIJ-1
         Shsq(L)  = (UL(L,IQ)-UL(L+1,IQ))**2 !interface
         dvsq(L)  = (UL(1,IQ)-UL(L  ,IQ))**2 !tracer pts
      END DO
      dVsq(LMIJ) = (UL(1,IQ)-UL(LMIJ,IQ))**2 ! tracer pts

C****    density related quantities
C****
C****    density of surface layer                                (kg/m3)
C****          rho   = rho{G1),S(1),zt(1)}
C****    local buoyancy gradient at km interfaces:                (m/s2)
C****          dbloc = g/rho{k+1,k+1} * [ rho{k,k+1}-rho{k+1,k+1} ]
C****    buoyancy difference with respect to "zref", the surface  (m/s2)
C****          dbsfc = g/rho{k,k} * [ rho{k,k}-rho{1,k} ] CORRECT
C****    thermal expansion coefficient without 1/rho factor    (kg/m3/C)
C****          talpha= d(rho{k,k})/d(T(k))
C****    salt expansion coefficient without 1/rho factor     (kg/m3/PSU)
C****          sbeta =  d(rho{k,k})/d(S(k))

      G(1)     = G0ML(1,IQ)*BYMML(1)
      S(1)     = S0ML(1,IQ)*BYMML(1)
      BYRHO(1) = VOLGSP(G(1),S(1),PO(1))
      RHO(1)   = 1d0/BYRHO(1)
      DO L=2,LMIJ
         G(L)     = G0ML(L,IQ)*BYMML(L)
         S(L)     = S0ML(L,IQ)*BYMML(L)
         BYRHO(L) = VOLGSP(G(L),S(L),PO(L))
         RHO(L)   = 1d0/BYRHO(L)
         RHOM = 1d0/VOLGSP(G(L-1),S(L-1),PO(L)) ! L-1 density wrt PO(L)
         RHO1 = 1d0/VOLGSP(G(1),S(1),PO(L)) ! surf density wrt PO(L)
         dbsfc(L)   = GRAV*(1d0-RHO1*BYRHO(L))
         dbloc(L-1) = GRAV*(1d0-RHOM*BYRHO(L))
C**** numerator of bulk richardson number on grid levels
         Ritop(L) = (zgrid(1)-zgrid(L)) * dbsfc(L)
      END DO

C**** Double diffusive option (ddmix) needs alphaDT,betaDS
C**** alphaDT = mean -talpha * delta(temp.)    at interfaces  (kg/m3)
C**** betaDS  = mean sbeta  * delta(salt)     at interfaces  (kg/m3)

      if (LDD) then
        TO(1)      =    TEMGSP(G(1),S(1),PO(1))
        talpha(1) =  ALPHAGSP(G(1),S(1),PO(1)) ! >0
        sbeta(1)  =   BETAGSP(G(1),S(1),PO(1))
        do L=2,LMIJ
          TO(L)      =    TEMGSP(G(L),S(L),PO(L))
          talpha(L) =  ALPHAGSP(G(L),S(L),PO(L)) ! >0
          sbeta(L)  =   BETAGSP(G(L),S(L),PO(L))
          alphaDT(L-1) = - 5d-1 * (talpha(L-1) + talpha(L))
     $         * (TO(L-1) - TO(L))
          betaDS (L-1) = 5d-1 * (sbeta (L-1) + sbeta(L))
     $         * (S(L-1) - S(L))*1d3
        end do
      end if

C**** Get diffusivities for the whole column
      CALL KPPMIX(LDD,ZE,zgrid,hwide,LMIJ,Shsq,dVsq,Ustar,Bo
     *     ,Bosol ,alphaDT,betaDS,dbloc,Ritop,Coriol,byhwide,
     *     AKVM,AKVS,AKVG,GHAT,HBL,KBL)

C**** Correct units for diffusivities (m^2/s) => (kg^2/m^4 s)
      DO L=1,LMIJ-1
         R = 5d-1*(RHO(L)+RHO(L+1))
         R2 = R**2
         AKVM(L) = AKVM(L)*R2
         AKVG(L) = AKVG(L)*R2
         AKVS(L) = AKVS(L)*R2
         GHAT(L) = 0.  ! no non-local transports since no surface fluxes
      END DO

C**** For each field (U,G,S + TRACERS) call OVDIFF
C**** Momentum
      CALL OVDIFF(UL(1,IQ),AKVM(1),GHAT,DTBYDZ,BYDZ2,LMIJ,UL0)
C**** Enthalpy
      CALL OVDIFFS(G0ML(1,IQ),AKVG(1),GHAT,DTBYDZ,BYDZ2,DTS
     *     ,LMIJ,G0ML0,FLG)
C**** Salinity
      CALL OVDIFFS(S0ML(1,IQ),AKVS(1),GHAT,DTBYDZ,BYDZ2,DTS
     *     ,LMIJ,S0ML0,FLS)
      IF ((ITER.eq.1  .or. ABS(HBLP-HBL).gt.(ZE(KBL)-ZE(KBL-1))*0.25)
     *     .and. ITER.lt.4) GO TO 510
#ifdef TRACERS_OCEAN
C**** Tracers are diffused after iteration (GHAT always zero)
      DO ITR = 1,NTM
        CALL OVDIFFS(TRML(1,ITR,IQ),AKVS(1),GHAT,DTBYDZ,BYDZ2
     *       ,DTS,LMIJ,TRML(1,ITR,IQ),FLT(1,ITR))
      END DO
#endif
C**** Implicitly apply interpolated KV to linear profile
C**** No surface fluxes
      DTBYDZ2 = 12d0*DTBYDZ(1)**2*BYDTS
      GZML(1,IQ)=(GZML(1,IQ)+3d0*FLG(1))/(1d0+DTBYDZ2*AKVG(1))
      SZML(1,IQ)=(SZML(1,IQ)+3d0*FLS(1))/(1d0+DTBYDZ2*AKVS(1))
#ifdef TRACERS_OCEAN
      TZML(1,:,IQ)=(TZML(1,:,IQ)+3d0*FLT(1,:))/(1d0+DTBYDZ2*AKVS(1))
#endif
      DO L=2,LMIJ-1
        DTBYDZ2 = 6d0*DTBYDZ(L)**2*BYDTS
        GZML(L,IQ)=(GZML(L,IQ)+3d0*(FLG(L-1)+FLG(L)))
     *       /(1d0+DTBYDZ2*(AKVG(L-1)+AKVG(L)))
        SZML(L,IQ)=(SZML(L,IQ)+3d0*(FLS(L-1)+FLS(L)))
     *       /(1d0+DTBYDZ2*(AKVS(L-1)+AKVS(L)))
#ifdef TRACERS_OCEAN
        TZML(L,:,IQ)=(TZML(L,:,IQ)+3d0*(FLT(L-1,:)+FLT(L,:)))
     *       /(1d0+DTBYDZ2*(AKVS(L-1)+AKVS(L)))
#endif
      END DO
      DTBYDZ2 = 12d0*DTBYDZ(LMIJ)**2*BYDTS
      GZML(LMIJ,IQ)=(GZML(LMIJ,IQ)+3d0*FLG(LMIJ-1))
     *     /(1d0+DTBYDZ2*AKVG(LMIJ-1))
      SZML(LMIJ,IQ)=(SZML(LMIJ,IQ)+3d0*FLS(LMIJ-1))
     *     /(1d0+DTBYDZ2*AKVS(LMIJ-1))
#ifdef TRACERS_OCEAN
      TZML(LMIJ,:,IQ)=(TZML(LMIJ,:,IQ)+3d0*FLT(LMIJ-1,:))
     *     /(1d0+DTBYDZ2*AKVS(LMIJ-1))
#endif
C****
      DO L=1,LMIJ-1
        OLNST(L,N,LN_KVG) = OLNST(L,N,LN_KVG) + AKVG(L)
        OLNST(L,N,LN_KVM) = OLNST(L,N,LN_KVM) + AKVM(L)
C       OLNST(L,N,LN_KVS) = OLNST(L,N,LN_KVS) + AKVS(L)
        OLNST(L,N,LN_WGFL)= OLNST(L,N,LN_WGFL) + FLG(L)
        OLNST(L,N,LN_WSFL)= OLNST(L,N,LN_WSFL) + FLS(L)
      END DO

C**** End of half box loops
  700 IQ=IQ+1
      IF(IQ.LE.2)  GO TO 240
C****
C**** Recreate main prognostic variables of potential enthalpy and
C**** salt from those on the half boxes
C****
      DO L=1,LMIJ
       MUST(L,N) = (  UL(L,2) +   UL(L,1))*MML(L)/DIST(N)
      G0MST(L,N) =  G0ML(L,2) + G0ML(L,1)
      NSIGG = EXPONENT(G0MST(L,N)) -1 - 42
      GXMST(L,N) = (G0ML(L,2) - G0ML(L,1))*BYBETA
      CALL REDUCE_FIG(NSIGG,GXMST(L,N))
      GZMST(L,N) =  GZML(L,2) + GZML(L,1)
      S0MST(L,N) =  S0ML(L,2) + S0ML(L,1)
      NSIGS = EXPONENT(S0MST(L,N)) -1 - 42 + 4
      SXMST(L,N) = (S0ML(L,2) - S0ML(L,1))*BYBETA
      CALL REDUCE_FIG(NSIGS,SXMST(L,N))
      SZMST(L,N) =  SZML(L,2) + SZML(L,1)
C****  limit salinity gradients
      if ( abs(SXMST(L,N)) > S0MST(L,N) )
     *     SXMST(L,N) = sign(S0MST(L,N),SXMST(L,N))
      if ( abs(SZMST(L,N)) > S0MST(L,N) )
     *     SZMST(L,N) = sign(S0MST(L,N),SZMST(L,N))
C****
#ifdef TRACERS_OCEAN
      DO ITR=1,NTM
        TRMST(L,N,ITR) = TRML(L,ITR,2) + TRML(L,ITR,1)
        NSIGT = EXPONENT(TRMST(L,N,ITR)) -1 - 42
        TXMST(L,N,ITR) =(TRML(L,ITR,2) - TRML(L,ITR,1))*BYBETA
        CALL REDUCE_FIG(NSIGT,TXMST(L,N,ITR))
        TZMST(L,N,ITR) = TZML(L,ITR,2) + TZML(L,ITR,1)
C****
        if (t_qlimit(itr)) then  ! limit gradients
          if ( abs(TXMST(L,N,ITR)) > TRMST(L,N,ITR) )
     *         TXMST(L,N,ITR) = sign(TRMST(L,N,ITR),TXMST(L,N,ITR))
          if ( abs(TZMST(L,N,ITR)) > TRMST(L,N,ITR) )
     *         TZMST(L,N,ITR) = sign(TRMST(L,N,ITR),TZMST(L,N,ITR))
        end if
C****
      END DO
#endif
      END DO
C**** End of outside loop over straits
  790 CONTINUE
      RETURN
      END SUBROUTINE STCONV

      SUBROUTINE OVDIFF(U,K,GHAT,DTBYDZ,BYDZ2,LMIJ,U0)
!@sum  OVDIFF Implicit vertical diff + non local transport for velocity
!@auth Gavin Schmidt
!@ver  1.0
      USE OCEAN, only : LMO
      IMPLICIT NONE
      REAL*8, DIMENSION(LMO), INTENT(IN) :: U0,K,GHAT,DTBYDZ,BYDZ2
      REAL*8, DIMENSION(LMO), INTENT(OUT) :: U
      INTEGER, INTENT(IN) :: LMIJ
      REAL*8, DIMENSION(LMO) :: A,B,C,R
      INTEGER L
C****  U0,U  input and output field (velocity or concentration)
C****     K  vertical diffusivity ((z)^2/ s)
C****  GHAT  non local transport of scalar
C****           kv * ghats * surface flux
C**** DTBYDZ  DT/DZ_L
C**** BYDZ2  1d0/DZ_L+1/2
C****    DT  timestep (s)
C**** Boundary conditions assumed to be no-flux at Z=0, Z=Z(LMIJ)
C**** Calculate operators for tridiagonal solver
      A(1) = 0
      B(1) = 1d0   + DTBYDZ(1)*BYDZ2(1)*K(1)
      C(1) =       - DTBYDZ(1)*BYDZ2(1)*K(1)
      R(1) = U0(1) - DTBYDZ(1)*GHAT(1)
      DO L=2,LMIJ-1
        A(L) =       - DTBYDZ(L)* BYDZ2(L-1)*K(L-1)
        B(L) = 1d0   + DTBYDZ(L)*(BYDZ2(L-1)*K(L-1)+BYDZ2(L)*K(L))
        C(L) =       - DTBYDZ(L)*                   BYDZ2(L)*K(L)
        R(L) = U0(L) + DTBYDZ(L)*(GHAT(L-1) - GHAT(L))
      END DO
      A(LMIJ) =          - DTBYDZ(LMIJ)*BYDZ2(LMIJ-1)*K(LMIJ-1)
      B(LMIJ) = 1d0      + DTBYDZ(LMIJ)*BYDZ2(LMIJ-1)*K(LMIJ-1)
      C(LMIJ) = 0
      R(LMIJ) = U0(LMIJ) + DTBYDZ(LMIJ)*GHAT(LMIJ-1)

      CALL TRIDIAG(A,B,C,R,U,LMIJ)

      RETURN
      END

      SUBROUTINE OVDIFFS(U,K,GHAT,DTBYDZ,BYDZ2,DT,LMIJ,U0,FL)
!@sum  OVDIFFS Implicit vertical diff + non local transport for tracers
!@auth Gavin Schmidt
!@ver  1.0
      USE OCEAN, only : LMO
      IMPLICIT NONE
      REAL*8, DIMENSION(LMO), INTENT(IN) :: U0,K,GHAT,DTBYDZ,BYDZ2
      REAL*8, DIMENSION(LMO), INTENT(OUT) :: U,FL
      REAL*8, INTENT(IN) :: DT
      INTEGER, INTENT(IN) :: LMIJ
      REAL*8, DIMENSION(LMO) :: A,B,C,R
      INTEGER L
C****  U0,U  input and output fields (total tracer)
C****     K  vertical diffusivity ((z)^2/ s)
C****  GHAT  non local transport of scalar
C****            kv * ghats * surface flux
C**** DTBYDZ  DT/DZ_L
C**** BYDZ2  1d0/DZ_L+1/2
C****    DT  timestep (s)
C****    FL  diffusive flux at boundary in units of total tracer
C**** Boundary conditions assumed to be no-flux at Z=0, Z=Z(LMIJ)
C**** Calculate operators for tridiagonal solver
      A(1) = 0
      B(1) = 1d0   + DTBYDZ(1)*BYDZ2(1)*K(1)
      C(1) =       - DTBYDZ(2)*BYDZ2(1)*K(1)
      R(1) = U0(1) - DT * GHAT(1)
      DO L=2,LMIJ-1
        A(L) =       - DTBYDZ(L-1)* BYDZ2(L-1)*K(L-1)
        B(L) = 1d0   + DTBYDZ(L  )*(BYDZ2(L-1)*K(L-1)+BYDZ2(L)*K(L))
        C(L) =       - DTBYDZ(L+1)*                   BYDZ2(L)*K(L)
        R(L) = U0(L) + DT * (GHAT(L-1) - GHAT(L))
      END DO
      A(LMIJ) =          - DTBYDZ(LMIJ-1)*BYDZ2(LMIJ-1)*K(LMIJ-1)
      B(LMIJ) = 1d0      + DTBYDZ(LMIJ  )*BYDZ2(LMIJ-1)*K(LMIJ-1)
      C(LMIJ) = 0
      R(LMIJ) = U0(LMIJ) + DT * GHAT(LMIJ-1)

      CALL TRIDIAG(A,B,C,R,U,LMIJ)

      DO L=1,LMIJ-1
        FL(L)=K(L)*(DTBYDZ(L+1)*U(L+1)-DTBYDZ(L)*U(L))*BYDZ2(L)
      END DO
C****
      RETURN
      END

      SUBROUTINE REDUCE_FIG(NSIG,RX)
!@sub reduce_fig reduce significant figures if calculation is garbage
!@+   made separate to avoid OMP compiler bug (due to NINT)
!@auth Gavin Schmidt
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NSIG
      REAL*8, INTENT(INOUT) :: RX

      IF (NSIG+30.gt.EXPONENT(RX)) RX =
     *     SCALE(REAL(NINT(SCALE(RX,-NSIG)),KIND=8),NSIG)

      RETURN
      END SUBROUTINE
