#include "rundeck_opts.h"

      MODULE OCEAN
!@sum  OCEAN dynamic ocean related variables
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
C**** Note that we currently use the same horizontal grid as for the
C**** atmosphere. However, we can redefine im,jm if necessary.
      USE MODEL_COM, only : im,jm,fim,byim
      USE GEOM, only : imaxj
#ifdef TRACERS_OCEAN
      USE TRACER_COM, only : ntm
#endif
      IMPLICIT NONE
      SAVE
      INTEGER, PARAMETER :: LMO=13  !@param LMO no. of ocean levels
!@param LSRPD level of solar radiation penetration depth
      INTEGER, PARAMETER :: LSRPD=3 
!@param ZE1,ZERAT control the grid spacing in the vertical
      REAL*8, PARAMETER :: ZE1=12d0,ZERAT=1.5d0

!@var MO mass of ocean (kg/m^2)
!@var UO E-W velocity on C-grid (m/s)
!@var VO N-S velocity on C-grid (m/s)
!@var G0M,GXMO,GYMO,GZMO pot. enthalpy of ocean (+moments) (J)
!@var S0M,SXMO,SYMO,SZMO salinity of ocean (+moments) (kg)
      REAL*8, DIMENSION(IM,JM,LMO) :: MO, UO, VO,
     *     G0M,GXMO,GYMO,GZMO, S0M,SXMO,SYMO,SZMO

C**** ocean geometry (should this be in a separate module?)
      REAL*8, DIMENSION(JM) :: DXYPO,DXPO,DYPO,DXVO,DYVO
     *     ,COSPO,SINPO,DXYVO,DXYSO,DXYNO,RAMVS,RAMVN,RLAT,BYDXYPO
      REAL*8, DIMENSION(0:JM) :: COSVO
      REAL*8, DIMENSION(IM) :: SINIC,COSIC
      REAL*8, DIMENSION(LMO) :: DSIGO,SIGO
      REAL*8, DIMENSION(0:LMO) :: SIGEO,ZE
      INTEGER, DIMENSION(IM,JM) :: LMM,LMU,LMV
!@var RATOC,ROCAT Ratio of areas for converting atm. fluxes to ocean
      REAL*8, DIMENSION(JM) :: RATOC,ROCAT
!@var J40S max. grid box below 40S (used in OPFIL)
      INTEGER :: J40S 

      REAL*8, DIMENSION(IM,JM) :: HATMO,HOCEAN,FOCEAN
C**** ocean related parameters
      INTEGER NDYNO,MDYNO,MSGSO
!@dbparam DTO timestep for ocean dynamics (s)
      REAL*8 :: DTO=450.        ! default. setable parameter
      REAL*8 DTOFS,DTOLF,DTS,BYDTS

!@var OPRESS Anomalous pressure at surface of ocean (under ice) (Pa)
      REAL*8, DIMENSION(IM,JM) :: OPRESS
!@var OGEOZ ocean geopotential at surface (m^2/s^2)
      REAL*8, DIMENSION(IM,JM) :: OGEOZ
      REAL*8, DIMENSION(IM,JM) :: OGEOZ_SV
!@var OPBOT ocean bottom pressure (diagnostic only) (Pa)
      REAL*8, DIMENSION(IM,JM) :: OPBOT

#ifdef TRACERS_OCEAN
!@var TRMO,TXMO,TYMO,TZMO tracer amount (+moments) in ocean (kg)      
      REAL*8, DIMENSION(IM,JM,LMO,NTM) :: TRMO,TXMO,TYMO,TZMO
#endif

      END MODULE OCEAN

      MODULE OCEAN_DYN
!@sum  OCEAN_DYN contains variables used in ocean dynamics
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0
      USE OCEAN, only : im,jm,lmo
      IMPLICIT NONE
      SAVE
C**** variables for the pressure gradient terms
!@var PO ocean pressure field
!@var PHI ocean geopotential
!@var DH height of each ocean layer
!@var DZGDP 
!@var VBAR mean specific volume of each layer
      REAL*8,DIMENSION(IM,JM,LMO) :: PO,PHI,DH,DZGDP,VBAR

C**** momentum and mass fluxes
!@var MMI initial mass field (kg)
      REAL*8, DIMENSION(IM,JM,LMO) :: MMI
!@var SMU,SMV,SMW integrated mass fluxes
      REAL*8, DIMENSION(IM,JM,LMO) :: SMU,SMV,SMW
!@var CONV mass flux convergence 
      REAL*8, DIMENSION(IM,JM,LMO) :: CONV
!@var MU,MV,MW instantaneous mass fluxes
      REAL*8, DIMENSION(IM,JM,LMO) :: MU,MV,MW
C****
      END MODULE OCEAN_DYN

      MODULE SW2OCEAN
!@sum  SW2OCEAN variables for putting solar radiation into ocean
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0
      USE OCEAN, only : ze
      IMPLICIT NONE
      SAVE
      INTEGER, PARAMETER :: LSRPD = 3
      REAL*8, DIMENSION(LSRPD) :: FSR,FSRZ,dFSRdZ,dFSRdZB
      REAL*8, PARAMETER :: RFRAC=.62d0, ZETA1=1.5d0, ZETA2=2d1

      CONTAINS

      SUBROUTINE init_solar
!@sum  init_solar calculates penetration of solar radiation
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0
      REAL*8 EF,EFZ,Z
      INTEGER L
      EF(Z) = RFRAC*EXP(-Z/ZETA1) + (1d0-RFRAC)*EXP(-Z/ZETA2)
      EFZ(Z)=ZETA1*RFRAC*EXP(-Z/ZETA1)+ZETA2*(1d0-RFRAC)*EXP(-Z/ZETA2)
C****
C**** Calculate the fraction of solar energy absorbed in each layer
C****
      do l=1,LSRPD
         FSR(l) = EF(ZE(l-1))
         FSRZ(l) = -3d0*(EF(ZE(l-1))+EF(ZE(l))) +
     *        6d0*(EFZ(ZE(l-1))-EFZ(ZE(l)))/(ZE(l)-ZE(l-1))
         dFSRdZB(l) = FSR(l)/(ZE(l)-ZE(l-1))
      end do
      do l=1,LSRPD-1
         dFSRdZ(l) = (FSR(l)-FSR(l+1))/(ZE(l)-ZE(l-1))
      end do
      dFSRdZ (LSRPD) = 0.
C****
      END SUBROUTINE init_solar
C****
      END MODULE SW2OCEAN
