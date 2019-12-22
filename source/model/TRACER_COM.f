#include "rundeck_opts.h"

!@sum  TRACER_COM: Exists alone to minimize the number of dependencies
!@+    This version for simple trace gases/chemistry/isotopes
!@auth Jean Lerner/Gavin Schmidt
!@ver  1.0

      MODULE TRACER_COM
!@sum  TRACER_COM tracer variables
!@auth Jean Lerner
!@ver  1.0
#ifdef TRACERS_ON
      USE QUSDEF, only: nmom
      USE MODEL_COM, only: im,jm,lm
#endif
      IMPLICIT NONE
      SAVE

C**** Each tracer has a variable name and a unique index
!@param NTM number of tracers
!@var TRNAME: Name for each tracer >>> MUST BE LEFT-JUSTIFIED <<<
#ifdef TRACERS_SPECIAL_O18
      integer, parameter :: ntm=3
      character*8, parameter :: trname(ntm)=(/
     *     'Water   ','H2O18   ','HDO     '/)  !,'HTO     '/)
#else     
#ifdef TRACERS_SPECIAL_Lerner     
      integer, parameter :: ntm=9
      character*8, parameter :: trname(ntm)= (/
     *     'Air     ','SF6     ','Rn222   ','CO2     ','N2O     ',
     *     'CFC11   ','14CO2   ','CH4     ','O3      '/)
#else
#ifdef TRACERS_SPECIAL_Shindell
      integer, parameter :: ntm=15
      character*8, parameter :: trname(ntm)=(/
     *    'Ox      ','NOx     ','N2O5    ','HNO3    ','H2O2    ',
     *    'CH3OOH  ','HCHO    ','HO2NO2  ','CO      ','CH4     ',
     *    'PAN     ','Isoprene','AlkylNit','Alkenes ','Paraffin'/)
#else
#ifdef TRACERS_WATER
      integer, parameter :: ntm=2
      character*8, parameter :: trname(ntm)=(/'Air     ','Water   '/)
#else
#ifdef TRACERS_OCEAN
      integer, parameter :: ntm=1
      character*8, parameter :: trname(ntm)=(/'Water   '/)
#else ! default for TRACERS_ON
      integer, parameter :: ntm=1
      character*8, parameter :: trname(ntm)=(/'Air     '/)
#endif
#endif
#endif
#endif
#endif
!@var N_XXX: variable names of indices for tracers
      integer :: 
     *     n_Air,    n_SF6,   n_Rn222, n_CO2,      n_N2O,
     *     n_CFC11,  n_14CO2, n_CH4,   n_O3,       n_water,
     *     n_H2O18,  n_HDO,   n_HTO,   n_Ox,       n_NOx, 
     *     n_N2O5,   n_HNO3,  n_H2O2,  n_CH3OOH,   n_HCHO,
     *     n_HO2NO2, n_CO,    n_PAN,   n_Isoprene, n_AlkylNit,
     *     n_Alkenes,n_Paraffin 

C****    The following are set in tracer_IC
!@var T_QLIMIT: if t_qlimit=.true. tracer is maintained as positive
      logical, dimension(ntm) :: t_qlimit
!@var trdecay radioactive decay constant (1/s) (=0 for stable tracers)
      real*8, dimension(ntm) :: trdecay
!@dbparam ITIME_TR0: start time for each tracer (hours)
      integer, dimension(ntm) :: itime_tr0
!@var MTRACE: timing index for tracers
      integer mtrace

#ifdef TRACERS_ON
!@var NTM_POWER: Power of 10 associated with each tracer (for printing)
      integer, dimension(ntm) :: ntm_power
!@var TR_MM: molecular mass of each tracer (g/mole)
      real*8, dimension(ntm) :: tr_mm
!@var needtrs: true if surface tracer value from PBL is required
      logical, dimension(ntm) :: needtrs

C**** Note units for these parameters!
C**** Example: clay dust; trpdens=2.5d3, trradius=0.73d-6 
C****          silt dust; trpdens=2.65d3, trradius=6.1d-6 
C****
!@var trpdens tracer particle density (kg/m^3) 
!@+               (=0 for non-particle tracers)
      real*8, dimension(ntm) :: trpdens
!@var trradius tracer effective radius (m) (=0 for non particle tracers)
      real*8, dimension(ntm) :: trradius

!@var TRM: Tracer array (kg)
      real*8, dimension(im,jm,lm,ntm) :: trm
!@var TRMOM: Second order moments for tracers (kg)
      real*8, dimension(nmom,im,jm,lm,ntm) :: trmom

!@var ntsurfsrcmax maximum number of surface 2D sources/sinks
      integer, parameter :: ntsurfsrcmax=14
!@var ntsurfsrc no. of non-interactive surface sources for each tracer
      integer, dimension(ntm) :: ntsurfsrc 
!@var nt3Dsrcmax maximum number of 3D tracer sources/sinks
#ifdef TRACERS_SPECIAL_Shindell
      integer, parameter :: nt3Dsrcmax=4
#else
      integer, parameter :: nt3Dsrcmax=3
#endif    
#endif

#ifdef TRACERS_WATER
!@param nWD_TYPES number of tracer types for wetdep purposes
!@param nGAS   index for wetdep tracer type = gas
!@param nPART  index for wetdep tracer type = particle/aerosol
!@param nWATER index for wetdep tracer type = water 
!@+       (original condense method)
      integer, parameter :: nWD_TYPES=3, nGAS=1, nPART=2, nWATER=3
!@param tr_evap_fact fraction of re-evaporation by tracer type
C note, tr_evap_fact is not dimensioned as NTM:
      REAL*8, parameter, dimension(nWD_TYPES) :: tr_evap_fact=
     *     (/1.d0, 0.5d0,  1.d0/)

!@var tr_wd_TYPE: tracer wet dep type (gas, particle, water)
      integer, dimension(ntm) :: tr_wd_TYPE
!@var tr_RKD: Henry's Law coefficient (in mole/Joule please !)
      real*8, dimension(ntm) :: tr_RKD
!@var tr_DHD: coefficient of temperature-dependence term of Henry's
!@+   Law coefficient (in Joule/mole please !)
      real*8, dimension(ntm) :: tr_DHD
!@var tr_H2ObyCH4 conc. of tracer in water from methane oxidation 
      real*8, dimension(ntm) :: tr_H2ObyCH4
!@var TRWM tracer in cloud liquid water amount (kg)
      real*8, dimension(im,jm,lm,ntm) :: trwm
#ifdef TRACERS_SPECIAL_O18
!@dbparam supsatfac factor controlling super saturation for isotopes
      real*8 :: supsatfac = 2d-3
#endif
#endif

#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
!@var TRW0 default tracer concentration in water (kg/kg)
      real*8, dimension(ntm) :: trw0
!@var NTROCN scaling power factor for ocean/ice tracer concentrations
      integer, dimension(ntm) :: ntrocn
#endif

#ifdef TRACERS_OCEAN
!@var TRGLAC tracer ratio in glacial runoff to ocean (kg/kg)
      real*8, dimension(ntm) :: trglac
#endif

      END MODULE TRACER_COM

