      MODULE DYNAMICS
!@sum  DYNAMICS contains all the pressure and momentum related variables
!@auth Original development team
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm
      IMPLICIT NONE
      SAVE
C**** Some helpful arrays (arrays should be L first)
!@var  PLIJ  Surface pressure: P(I,J) or PSF-PTOP (mb)
      REAL*8, DIMENSION(LM,IM,JM) :: PLIJ
!@var  PDSIG  Surface pressure * DSIG(L) (mb)
      REAL*8, DIMENSION(LM,IM,JM) :: PDSIG
!@var  AM  Air mass of each box (kg/m^2)
      REAL*8, DIMENSION(LM,IM,JM) :: AM     ! PLIJ*DSIG(L)*100/grav
!@var  BYAM  1/Air mass (m^2/kg)
      REAL*8, DIMENSION(LM,IM,JM) :: BYAM
!@var  PMID  Pressure at mid point of box (mb)
      REAL*8, DIMENSION(LM,IM,JM) :: PMID    ! SIG(L)*PLIJ+PTOP
!@var  PK   PMID**KAPA
      REAL*8, DIMENSION(LM,IM,JM) :: PK
!@var  PEUP  Pressure at lower edge of box (incl. surface) (mb)
      REAL*8, DIMENSION(LM+1,IM,JM) :: PEDN  ! SIGE(L)*PLIJ+PTOP
!@var  PEK  PEUP**KAPA
      REAL*8, DIMENSION(LM+1,IM,JM) :: PEK
!@var  SQRTP  square root of P (used in diagnostics)
      REAL*8, DIMENSION(IM,JM) :: SQRTP
!@var  PTROPO  Pressure at mid point of tropopause level (mb)
      REAL*8, DIMENSION(IM,JM) :: PTROPO
!@var  LTROPO  Tropopause layer
      INTEGER, DIMENSION(IM,JM) :: LTROPO

C**** module should own dynam variables used by other routines
!@var PTOLD pressure at beginning of dynamic time step (for clouds)
      REAL*8, DIMENSION(IM,JM)    :: PTOLD
!@var SD_CLOUDS vert. integrated horizontal convergence (for clouds)
      REAL*8, DIMENSION(IM,JM,LM) :: SD_CLOUDS
!@var GZ geopotential height (for Clouds and Diagnostics)
      REAL*8, DIMENSION(IM,JM,LM) :: GZ
!@var DPDX_BY_RHO,DPDY_BY_RHO (pressure gradients)/density at L=1
      REAL*8, DIMENSION(IM,JM)  :: DPDX_BY_RHO,DPDY_BY_RHO
!@var DPDX_BY_RHO_0,DPDY_BY_RHO_0 surface (pressure gradients)/density
      REAL*8, DIMENSION(IM,JM)  :: DPDX_BY_RHO_0,DPDY_BY_RHO_0

      REAL*8, DIMENSION(IM,JM,LM) :: PU,PV,CONV
      REAL*8, DIMENSION(IM,JM,LM-1) :: SD
!@var PIT  pressure tendency (mb m^2/s)
      REAL*8, DIMENSION(IM,JM) :: PIT
      EQUIVALENCE (SD(1,1,1),CONV(1,1,2))
      EQUIVALENCE (PIT(1,1),CONV(1,1,1))

      REAL*8, DIMENSION(IM,JM,LM) :: PHI,SPA
      REAL*8, DIMENSION(IM,JM,LM) :: DUT,DVT
!@var xAVRX scheme-depend. coefficient for AVRX: 1,byrt2 (2nd,4th order)
      REAL*8 xAVRX

!@var PUA,PVA,SDA,PS save PU,PV,SD,P for hourly tracer advection
!@var MB Air mass array for tracers (before advection)
!@var MA Air mass array for tracers (updated during advection)
      REAL*8, DIMENSION(IM,JM,LM) :: PUA,PVA,SDA,MB,MA
      REAL*8, DIMENSION(IM,JM) :: PS

!@var DKE change in KE due to dissipation (SURF/DC/MC) (m^2/s^2)
      REAL*8, DIMENSION(IM,JM,LM) :: DKE

      END MODULE DYNAMICS

