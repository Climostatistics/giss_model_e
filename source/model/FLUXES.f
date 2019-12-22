#include "rundeck_opts.h"

      MODULE FLUXES
!@sum  FLUXES contains the fluxes between various components
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      USE TRACER_COM, only: ntm
#endif
#ifdef TRACERS_ON
     *     ,ntsurfsrcmax,nt3Dsrcmax
#endif
      IMPLICIT NONE

!@var RUNOSI run off from sea/lake ice after surface (kg/m^2)
!@var ERUNOSI energy of run off from sea/lake ice after surface (J/m^2)
!@var SRUNOSI salt in run off from sea/lake ice after surface (kg/m^2)
      REAL*8, DIMENSION(IM,JM) :: RUNOSI, ERUNOSI, SRUNOSI
!@var RUNPSI run off from sea/lake ice after precip (kg/m^2)
!@var ERUNPSI energy of run off from sea/lake ice after precip (J/m^2)
!@var SRUNPSI salt in run off from sea/lake ice after precip (kg/m^2)
      REAL*8, DIMENSION(IM,JM) :: RUNPSI, SRUNPSI   ! ,ERUNPSI (not yet)
!@var RUNOE run off from earth (kg/m^2)
!@var ERUNOE energy of run off from earth (J/m^2)
      REAL*8, DIMENSION(IM,JM) :: RUNOE, ERUNOE
C**** DMSI,DHSI,DSSI are fluxes for ice formation within water column
!@var DMSI mass flux of sea ice 1) open water and 2) under ice (kg/m^2)
!@var DHSI energy flux of sea ice 1) open water and 2) under ice (J/m^2)
!@var DSSI salt flux in sea ice 1) open water and 2) under ice (kg/m^2)
      REAL*8, DIMENSION(2,IM,JM) :: DMSI, DHSI, DSSI
!@var fmsi_io,fhsi_io,fssi_io basal ice-ocean fluxes (kg or J/m^2)
      REAL*8, DIMENSION(IM,JM) :: fmsi_io,fhsi_io,fssi_io
!@var RUNOLI run off from land ice (kg/m^2) (Energy always=0)
      REAL*8, DIMENSION(IM,JM) :: RUNOLI

C**** surface energy fluxes defined over type
!@param NSTYPE number of surface types for radiation purposes
      INTEGER, PARAMETER :: NSTYPE=4
!@var E0 net energy flux at surface for each type (J/m^2)
!@var E1 net energy flux at layer 1 for each type (J/m^2)
      REAL*8, DIMENSION(IM,JM,NSTYPE) :: E0,E1
!@var EVAPOR evaporation over each type (kg/m^2) 
      REAL*8, DIMENSION(IM,JM,NSTYPE) :: EVAPOR
!@var SOLAR absorbed solar radiation (J/m^2)
!@+   SOLAR(1)  absorbed by open water
!@+   SOLAR(2)  absorbed by ice
!@+   SOLAR(3)  absorbed by water under the ice
      REAL*8, DIMENSION(3,IM,JM) :: SOLAR

C**** Momemtum stresses are calculated as if they were over whole box
!@var DMUA,DMVA momentum flux from atmosphere for each type (kg/m s) 
!@+   On atmospheric A grid (tracer point)
      REAL*8, DIMENSION(IM,JM,NSTYPE) :: DMUA,DMVA
!@var DMUI,DMVI momentum flux from sea ice to ocean (kg/m s)
!@+   On atmospheric C grid 
      REAL*8, DIMENSION(IM,JM) :: DMUI,DMVI
!@var UI2rho Ustar*2*rho ice-ocean friction velocity on atmospheric grid
      REAL*8, DIMENSION(IM,JM) :: UI2rho
!@var OGEOZA ocean surface height geopotential (m^2/s^2) (on ATM grid)
      REAL*8, DIMENSION(IM,JM) :: OGEOZA

C**** currently saved - should be replaced by fluxed quantities
!@var DTH1,DQ1 heat/water flux from atmos. summed over type (C, kg/kg)
      REAL*8, DIMENSION(IM,JM) :: DTH1,DQ1

!@var uflux1 surface turbulent u-flux (=-<uw>) 
!@var vflux1 surface turbulent v-flux (=-<vw>)
!@var tflux1 surface turbulent t-flux (=-<tw>)
!@var qflux1 surface turbulent q-flux (=-<qw>)
      real*8, dimension(im,jm) :: uflux1,vflux1,tflux1,qflux1

C**** The E/FLOWO and E/S/MELTI arrays are used to flux quantities to 
C**** the ocean that are not tied to the open water/ice covered 
C**** fractions. This is done separately for river flow and complete
C**** sea ice melt.
!@var FLOWO,EFLOWO mass, energy from rivers into ocean (kg, J)
      REAL*8, DIMENSION(IM,JM) :: FLOWO,EFLOWO
!@var MELTI,EMELTI,SMELTI mass,energy,salt from simelt into ocean (kg,J)
      REAL*8, DIMENSION(IM,JM) :: MELTI,EMELTI,SMELTI

!@var PREC precipitation (kg/m^2)
      REAL*8, DIMENSION(IM,JM) :: PREC
!@var EPREC energy of preciptiation (J/m^2)
      REAL*8, DIMENSION(IM,JM) :: EPREC
!@var PRECSS precipitation from super-saturation (kg/m^2)
      REAL*8, DIMENSION(IM,JM) :: PRECSS

!@var GTEMP ground temperature (upper two levels) over surface type (C)
      REAL*8, DIMENSION(2,NSTYPE,IM,JM) :: GTEMP
!@var SSS sea surface salinity on atmospheric grid (ppt)
      REAL*8, DIMENSION(IM,JM) :: SSS
!@var MLHC ocean mixed layer heat capacity (J/m^2 C) 
      REAL*8, DIMENSION(IM,JM) :: MLHC
!@var UOSURF, VOSURF ocean surface velocity (Atm C grid) (m/s)
      REAL*8, DIMENSION(IM,JM) :: UOSURF,VOSURF
!@var UISURF, VISURF dynamic ice surface velocity (Atm C grid) (m/s)
      REAL*8, DIMENSION(IM,JM) :: UISURF,VISURF
!@var APRESS total atmos + sea ice pressure (at base of sea ice) (Pa)
      REAL*8, DIMENSION(IM,JM) :: APRESS
!@var FWSIM fresh water sea ice mass (kg/m^2) (used for qflux model)
      REAL*8, DIMENSION(IM,JM) :: FWSIM
!@var MSICNV fresh water sea ice mass convergence after advsi (kg/m^2)
!@+   (used for qflux model)
      REAL*8, DIMENSION(IM,JM) :: MSICNV

#ifdef TRACERS_ON
!@var TRSOURCE non-interactive surface sources/sinks for tracers (kg/s)
      REAL*8, DIMENSION(IM,JM,ntsurfsrcmax,NTM):: trsource
!@var TRSRFFLX interactive surface sources/sinks for tracers (kg/s)
      REAL*8, DIMENSION(IM,JM,NTM):: trsrfflx
!@var TRFLUX1 total surface flux for each tracer (kg/s)
      REAL*8, DIMENSION(IM,JM,NTM):: trflux1
!@var TRGRDEP gravitationally settled tracers at surface (kg)
      REAL*8, DIMENSION(NTM,IM,JM):: TRGRDEP
!@var GTRACER ground concentration of tracer on atmospheric grid (kg/kg)
      REAL*8, DIMENSION(NTM,NSTYPE,IM,JM):: GTRACER
!@var TR3DSOURCE 3D sources/sinks for tracers (kg/s)
      REAL*8, DIMENSION(IM,JM,LM,nt3Dsrcmax,NTM):: tr3Dsource

#ifdef TRACERS_WATER
!@var TRPREC tracers in precip (kg)
      REAL*8, DIMENSION(NTM,IM,JM):: TRPREC
!@var TREVAPOR tracer evaporation over each type (kg/m^2) 
      REAL*8, DIMENSION(NTM,NSTYPE,IM,JM) :: TREVAPOR
!@var TRUNPSI tracer in run off from sea/lake ice after precip (kg/m^2)
!@var TRUNOSI tracer in run off from sea/lake ice after surface (kg/m^2)
!@var TRUNOE tracer runoff from earth (kg/m^2)
!@var TRUNOLI tracer runoff from land ice (kg/m^2)
      REAL*8, DIMENSION(NTM,IM,JM):: TRUNPSI, TRUNOSI, TRUNOE, TRUNOLI
!@var TRFLOWO tracer in river runoff into ocean (kg)
      REAL*8, DIMENSION(NTM,IM,JM) :: TRFLOWO
!@var TRMELTI tracer from simelt into ocean (kg)
      REAL*8, DIMENSION(NTM,IM,JM) :: TRMELTI
!@var ftrsi_io ice-ocean tracer fluxes under ice (kg/m^2)
      REAL*8, DIMENSION(NTM,IM,JM) :: ftrsi_io
#endif
#endif

#if (defined TRACERS_OCEAN) || (defined TRACERS_WATER)
!@var DTRSI tracer flux in sea ice under ice and on open water (kg/m^2)
      REAL*8, DIMENSION(NTM,2,IM,JM) :: DTRSI
#endif

      END MODULE FLUXES

