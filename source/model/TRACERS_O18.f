#include "rundeck_opts.h"

!@sum  TRACERS_O18 water isotope specific routines and functions
!@auth Gavin Schmidt

      FUNCTION FRACVL(TEMP,trname)
!@sum FRACVL Calculate vapor-->liquid equilibrium fractionation factor
!@auth Gavin Schmidt
      USE CONSTANT, only : tf
      IMPLICIT NONE
!@var TEMP  temperature (deg C)
      REAL*8, INTENT(IN) :: TEMP
      CHARACTER, INTENT(IN) :: trname*8
      INTEGER, PARAMETER :: NTSPM=4
c**** quadratic fit to Majoube (1971)
c      REAL*8, PARAMETER ::
c     *     A(NTSPM) = (/  0d0, -3.75d-7, -6.375d-6, -6.875d-6  /),
c     *     B(NTSPM) = (/  0d0,  1.025d-4, 1.2475d-3, 1.7d-3    /),
c     *     C(NTSPM) = (/  1d0,  0.9884d0, 0.9001d0 , 0.86975d0 /)
c**** exponentials 
      REAL*8, PARAMETER ::
     *     A(NTSPM) = (/  0d0,  1137d0   ,  24844d0  , 46480d0  /),
     *     B(NTSPM) = (/  0d0, -0.4156d0 , -76.248d0 ,-103.87d0 /),
     *     C(NTSPM) = (/  0d0, -2.0667d-3,  52.612d-3, 0d0      /)

!@var ITR   species number of tracer
C****       1: Fresh Water 2: o18 3: Deu  4: Tritium
      INTEGER ITR
      REAL*8 FRACVL,TK
C****
      select case (trname)
      case ('Water')
        ITR=1
      case ('H2O18')
        ITR=2
      case ('HDO')
        ITR=3
      case ('HTO')   ! tritium data not yet confirmed
        ITR=4
      case default
        write(6,*) "Tracer name ",trname," not defined in FRACVL"
        call stop_model('Tracer name not defined in FRACVL',255)
      end select
C**** Quadratic fit
c      FRACVL=C(ITR) + TEMP*(B(ITR) + TEMP*A(ITR))
C**** Exponential
      TK=TEMP+TF
      FRACVL=EXP(-A(ITR)/TK**2 - B(ITR)/TK - C(ITR))
C****
      RETURN
      END FUNCTION FRACVL

      FUNCTION FRACVS(TEMP,trname)
!@sum FRACVS Calculate vapour --> solid (ice) equil. fractionation fact.
!@auth Gavin Schmidt
      USE CONSTANT, only : tf
      IMPLICIT NONE
!@var TEMP  temperature (deg C)
      REAL*8, INTENT(IN) :: TEMP
      CHARACTER, INTENT(IN) :: trname*8
      INTEGER, PARAMETER :: NTSPM=4
C**** linear fit
c      REAL*8, PARAMETER ::
c     *     A(NTSPM) = (/  0d0,  1.36d-4,   1.46d-3, 0d0/),
c     *     B(NTSPM) = (/  1d0,  0.9850d0, 0.8834d0, 0.7845d0/)
C**** exponential
      REAL*8, PARAMETER ::
     *     A(NTSPM) = (/  0d0,  0d0       , 16288d0 , 46480d0 /),
     *     B(NTSPM) = (/  0d0,  11.839d0  , 0d0     ,-103.87d0/),
     *     C(NTSPM) = (/  0d0, -0.028244d0,-0.0934d0, 0d0     /)

!@var ITR   species number of tracer
C****       1: Fresh Water 2: o18 3: Deu  4: Tritium
      INTEGER ITR
      REAL*8 FRACVS,TK
C****
      select case (trname)
      case ('Water')
        ITR=1
      case ('H2O18')
        ITR=2
      case ('HDO')
        ITR=3
      case ('HTO')   ! tritium data not yet confirmed
        ITR=4
      case default
        write(6,*) "Tracer name ",trname," not defined in FRACVS"
        call stop_model('Tracer name not defined in FRACVS',255)
      end select
C**** linear fit
c      FRACVS=B(ITR) + A(ITR)*TEMP
C**** Exponential
      TK=TEMP+TF
      FRACVS=EXP(-A(ITR)/TK**2 - B(ITR)/TK -C(ITR))
C****
      RETURN
      END FUNCTION FRACVS

      FUNCTION FRACLS(trname)
!@sum FRACLS Calculate liquid --> solid equilibrium fractionation factor
!@auth Gavin Schmidt
      IMPLICIT NONE
      CHARACTER, INTENT(IN) :: trname*8
      INTEGER, PARAMETER :: NTSPM=4
      REAL*8, PARAMETER ::
     *     A(NTSPM) = (/ 1d0, 1.0035d0, 1.0208d0,  1.04d0/)
!@var ITR   species number of tracer
C****       1: Fresh Water 2: o18 3: Deu  4: Tritium
      INTEGER ITR
      REAL*8 FRACLS
C****
      select case (trname)
      case ('Water')
        ITR=1
      case ('H2O18')
        ITR=2
      case ('HDO')
        ITR=3
      case ('HTO')
        ITR=4
      case default
        FRACLS=1.    ! no fractionation
        RETURN
      end select
      FRACLS=A(ITR)
C****
      RETURN
      END FUNCTION FRACLS

      FUNCTION FRACLK(WS,trname)
!@sum FRACLK calculates the liquid/vapor kinetic fractionation factor
!@+          from either (Merlivat and Jouzel,1972) or as a constant
!@auth Gavin Schmidt
      IMPLICIT NONE
!@var WS surface wind speed (m/s)
      REAL*8, INTENT(IN) :: WS
      CHARACTER, INTENT(IN) :: trname*8
      INTEGER, PARAMETER :: NTSPM=4
      REAL*8, PARAMETER ::
     *     A(NTSPM) = (/ 1d0, 0.994d0,  0.99472d0,  0.98944d0  /),
     *     B(NTSPM) = (/ 0d0, 0.285d-3, 0.2508d-3,  0.5016d-3  /),
     *     C(NTSPM) = (/ 0d0, 0.82d-3,  0.7216d-3,  0.14432d-2 /)
!@var ITR   species number of tracer
C****       1: Fresh Water 2: o18 3: Deu  4: Tritium
      INTEGER ITR
      REAL*8 FRACLK
C****
C**** Calculate kinetic fractionation factor:
C****
      select case (trname)
      case ('Water')
        ITR=1
      case ('H2O18')
        ITR=2
      case ('HDO')
        ITR=3
      case ('HTO')   ! tritium data not yet confirmed
        ITR=4
      case default
        write(6,*) "Tracer name ",trname," not defined in FRACLK"
        call stop_model('Tracer name not defined in FRACLK',255)
      end select
      FRACLK=A(ITR)
      IF(WS.GE.7.) FRACLK=1d0-(B(ITR)*WS+C(ITR))
C****
      RETURN
      END FUNCTION FRACLK

#ifdef TRACERS_SPECIAL_O18
      SUBROUTINE ISOEQUIL(N,TEMP,LIQU,QMV,QML,TRMV,TRML,FEQ)
!@sum ISOEQUIL equilibrates isotopes in vapor & liquid/solid reservoirs
!@auth Gavin Schmidt/Georg Hoffmann
      USE CONSTANT, only : tf
      USE TRACER_COM, only: trname,tr_wd_TYPE,nWATER
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N
      LOGICAL, INTENT(IN) :: LIQU !@var LIQU true if QML is liquid
      REAL*8, INTENT(IN) :: TEMP  !@var TEMP temperature (K)
!@var QMV,QML vapour and liquid water masses (kg or kg/m^2)
      REAL*8, INTENT(IN) :: QMV,QML
!@var TRMV,TRML vapour and liquid tracer mass (kg or kg/m^2)
      REAL*8, INTENT(INOUT) :: TRMV,TRML
!@var FEQ fraction of condensate equilibrated (1. = all, 0. = none)
      REAL*8, INTENT(IN) :: FEQ
      REAL*8 TDEGC,ZALPH,ZDELEQU,ZXFAC,FRACVL,FRACVS

      SELECT CASE(tr_wd_TYPE(N))
      CASE(nWATER)
        TDEGC = TEMP - TF
        IF (QMV.GT.0.) THEN
          IF (LIQU) THEN
            ZALPH = 1./FRACVL(TDEGC,trname(N))
          ELSE
            ZALPH = 1./FRACVS(TDEGC,trname(N))
          END IF
          ZXFAC = FEQ*ZALPH*QML/QMV
          ZDELEQU = (FEQ*TRML - ZXFAC*TRMV)/(1.+ZXFAC)
          TRML = TRML - ZDELEQU
          TRMV = TRMV + ZDELEQU
        END IF
      END SELECT

      RETURN
C****
      END SUBROUTINE ISOEQUIL

      FUNCTION KIN_COND_ICE(ALPH,SUPSAT,trname)
!@sum calculate kinetic fractionation when condensing to ice in 
!@+   super-saturated conditions
!@auth Gavin Schmidt/Georg Hoffmann
      USE CONSTANT, only : tf
      IMPLICIT NONE
      CHARACTER, INTENT(IN) :: trname*8
!@var SUPSAT super_saturation factor from cloud scheme
      REAL*8, INTENT(IN) :: SUPSAT
!@var ALPH equilibrium fractionation
      REAL*8, INTENT(IN) :: ALPH
!@var ZDIFREL = inverse ratio of diffusion coeffs w.r.t normal water
      real*8 :: ZDIFREL(4) = (/ 1d0 ,1.0285d0, 1.0251d0, 1.0331d0/)
      integer itr
      real*8 kin_cond_ice
C****
C**** Calculate kinetic condensation when condensing to ice
C****
      select case (trname)
      case ('Water')
        ITR=1
      case ('H2O18')
        ITR=2
      case ('HDO')
        ITR=3
      case ('HTO')
        ITR=4
      case default
        write(6,*) "Tracer name ",trname," not defined in KIN_COND"
        call stop_model('Tracer name not defined in KIN_COND',255)
      end select
      KIN_COND_ICE=alph*SUPSAT/(1.+(SUPSAT-1.)*alph*ZDIFREL(ITR))

      return
      end function kin_cond_ice

      FUNCTION KIN_EVAP_PREC(ALPH,HEFF,trname)
!@sum calculate kinetic fractionation when evaporating into 
!@+   undersaturated environment
!@auth Gavin Schmidt/Georg Hoffmann
      USE CONSTANT, only : tf
      IMPLICIT NONE
      CHARACTER, INTENT(IN) :: trname*8
!@var HEFF effective relative humidity from cloud scheme
      REAL*8, INTENT(IN) :: HEFF
!@var alph equilibrium fractionation 
      REAL*8, INTENT(IN) :: ALPH
!@var ZDIFRELGAM = ZDIFREL^0.58 
      real*8 :: ZDIFRELGAM(4) = (/ 1d0, 1.0164d0, 1.0145d0, 1.0191d0/)
      integer itr
      real*8 kin_evap_prec
C****
C**** Calculate kinetic condensation when evaporating below clouds
C****
      select case (trname)
      case ('Water')
        ITR=1
      case ('H2O18')
        ITR=2
      case ('HDO')
        ITR=3
      case ('HTO')
        ITR=4
      case default
        write(6,*) "Tracer name ",trname," not defined in KIN_EVAP"
        call stop_model('Tracer name not defined in KIN_EVAP',255)
      end select
      KIN_EVAP_PREC=alph*HEFF/(1.+(HEFF-1.)*alph*ZDIFRELGAM(ITR))

      return
      end function kin_evap_prec

      FUNCTION delta(W,T,n)
!@sum calculate per mil values
      use TRACER_COM, only : trw0
      implicit none
      real*8, intent(in) :: W,T
      integer, intent(in) :: n
      real*8 delta

      delta=0.
      if (trw0(n).gt.0 .and. W.gt.0) delta=(T/(W*trw0(n))-1.)*1d3

      return
      end function delta
#endif
