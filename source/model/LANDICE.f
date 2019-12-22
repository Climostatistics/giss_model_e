#include "rundeck_opts.h"

      MODULE LANDICE
!@sum  LANDICE contains variables and routines dealing with land ice
!@auth Original Development Team
!@ver  1.0
!@cont PRECLI,LNDICE
      USE CONSTANT, only : lhm,rhoi,shi
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm
#endif
      IMPLICIT NONE
      SAVE
!@var Z1E,Z2LI first and second layer thicknesses for land ice
      REAL*8, PARAMETER :: Z1E=.1d0, Z2LI=2.9d0
!@var ACE1LI first layer land ice mass (kg/m^2)
      REAL*8, PARAMETER :: ACE1LI = Z1E*RHOI
!@var HC1LI heat capacity of first layer land ice (J/m^2)
      REAL*8, PARAMETER :: HC1LI = ACE1LI*SHI
!@var ACE2LI second layer land ice mass (kg/m^2)
      REAL*8, PARAMETER :: ACE2LI = Z2LI*RHOI
!@var HC2LI heat capacity of second layer land ice (J/m^2)
      REAL*8, PARAMETER :: HC2LI = ACE2LI*SHI
!@dbparam glmelt_on determines whether glacial melt is used for oceans
      INTEGER :: glmelt_on = 0   ! default is 'off'

      CONTAINS

      SUBROUTINE PRECLI(SNOW,TG1,TG2,PRCP,ENRGP,
#ifdef TRACERS_WATER
     *     TRSNOW,TRLI,TRPRCP,TRDIFS,TRUN0,
#endif
     *     EDIFS,DIFS,ERUN2,RUN0)
!@sum  PRECLI apply precipitation to land ice fraction
!@auth Original Development Team
!@ver  1.0
      REAL*8, INTENT(INOUT) :: SNOW,TG1,TG2
      REAL*8, INTENT(IN) :: PRCP,ENRGP
      REAL*8, INTENT(OUT) :: EDIFS,DIFS,ERUN2,RUN0
#ifdef TRACERS_WATER
!@var TRSNOW tracer amount in snow (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(INOUT) :: TRSNOW
!@var TRLI tracer amount in land ice (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(INOUT) :: TRLI
!@var TRPRCP tracer amount in precip (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(IN) :: TRPRCP
!@var TRUN0 tracer runoff from ice (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(OUT) :: TRUN0
!@var TRDIFS implicit tracer flux at base of ice (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(OUT) :: TRDIFS
#endif

      REAL*8 HC1,DWATER
C**** initiallize output
      EDIFS=0. ; DIFS=0. ; ERUN2=0. ; RUN0=0.
#ifdef TRACERS_WATER
      TRDIFS(:)=0. ; TRUN0(:)=0.
#endif

      HC1=HC1LI+SNOW*SHI
      IF (ENRGP.GE.0.) THEN
        IF (ENRGP.GE.-TG1*HC1) THEN
C**** RAIN HEATS UP TG1 TO FREEZING POINT AND MELTS SOME SNOW OR ICE
          DWATER=(TG1*HC1+ENRGP)/LHM
          TG1=0.
          RUN0=DWATER+PRCP
          IF (DWATER.LE.SNOW) THEN
C**** RAIN MELTS SOME SNOW
#ifdef TRACERS_WATER
            TRUN0(:)=TRPRCP(:)+DWATER*TRSNOW(:)/SNOW
            TRSNOW(:)=TRSNOW(:)*(1.-DWATER/SNOW)
#endif
            SNOW=SNOW-DWATER
          ELSE
C**** RAIN MELTS ALL SNOW AND SOME ICE, ICE MOVES UP THROUGH THE LAYERS
            DIFS=SNOW-DWATER  ! < 0
            SNOW=0.
            TG1=-TG2*DIFS/ACE1LI
            EDIFS=DIFS*(TG2*SHI-LHM)
            ERUN2=EDIFS
#ifdef TRACERS_WATER
            TRUN0(:)=TRPRCP(:)+TRSNOW(:)-DIFS*TRLI(:)/(ACE2LI+ACE1LI)
            TRDIFS(:)=DIFS*TRLI(:)/(ACE2LI+ACE1LI)
            TRSNOW(:)=0.
#endif
          END IF
        ELSE
C**** RAIN COOLS TO FREEZING POINT AND HEATS UP TG1
          TG1=TG1+ENRGP/HC1
          RUN0=PRCP
#ifdef TRACERS_WATER
          TRUN0(:)=TRPRCP(:)
#endif
        END IF
      ELSE
C**** SNOW INCREASES SNOW AMOUNT AND SNOW TEMPERATURE RECOMPUTES TG1
        TG1=(TG1*HC1+ENRGP+LHM*PRCP)/(HC1+PRCP*SHI)
        SNOW=SNOW+PRCP
#ifdef TRACERS_WATER
        TRSNOW(:)=TRSNOW(:)+TRPRCP(:)
#endif
        IF (SNOW.gt.ACE1LI) THEN
C**** SNOW IS COMPACTED INTO ICE, ICE MOVES DOWN THROUGH THE LAYERS
          DIFS=SNOW-.9d0*ACE1LI
#ifdef TRACERS_WATER
          TRDIFS(:)=TRLI(:)*DIFS/(ACE2LI+ACE1LI)
          TRLI(:)=TRLI(:)*(1.-DIFS/(ACE2LI+ACE1LI))+TRSNOW(:)*DIFS/SNOW
          TRSNOW(:)=TRSNOW(:)*(1.-DIFS/SNOW)
#endif
          SNOW=.9d0*ACE1LI
          EDIFS=DIFS*(TG1*SHI-LHM)
          ERUN2=DIFS*(TG2*SHI-LHM)
          TG2=TG2+(TG1-TG2)*DIFS/ACE2LI
        END IF
      END IF
      RETURN
      END SUBROUTINE PRECLI

      SUBROUTINE LNDICE(SNOW,TG1,TG2,F0DT,F1DT,EVAP,
#ifdef TRACERS_WATER
     *     TRSNOW,TRLI,TREVAP,TRDIFS,TRUN0,
#endif
     *     EDIFS,DIFS,RUN0)
!@sum  LNDICE apply surface fluxes to land ice fraction
!@auth Original Development team
!@ver  1.0
      REAL*8, INTENT(INOUT) :: SNOW,TG1,TG2
      REAL*8, INTENT(IN) :: F0DT,F1DT,EVAP
      REAL*8, INTENT(OUT) :: EDIFS,DIFS,RUN0
#ifdef TRACERS_WATER
C**** Tracer concentration in ice assumed constant
!@var TRSNOW tracer amount in snow (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(INOUT) :: TRSNOW
!@var TRLI tracer amount in land ice (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(INOUT) :: TRLI
!@var TREVAP tracer amount in evaporation (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(IN) :: TREVAP
!@var TRUN0 tracer runoff from ice (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(OUT) :: TRUN0
!@var TRDIFS implicit tracer flux at base of ice (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(OUT) :: TRDIFS
#endif
      REAL*8 :: SNANDI,HC1,ENRG1
      INTEGER N

C**** initiallize output
      EDIFS=0. ; DIFS=0. ; RUN0=0.
#ifdef TRACERS_WATER
      TRDIFS(:)=0. ; TRUN0(:)=0.
#endif

C**** CALCULATE TG1
      SNANDI=SNOW+ACE1LI-EVAP
      HC1=SNANDI*SHI
      ENRG1=F0DT+EVAP*(TG1*SHI-LHM)-F1DT
      IF (ENRG1.GT.-TG1*HC1) THEN
C**** FLUXES HEAT UP TG1 TO FREEZING POINT AND MELT SOME SNOW AND ICE
        RUN0=(ENRG1+TG1*HC1)/LHM
        TG1=0.
        SNANDI=SNANDI-RUN0
#ifdef TRACERS_WATER
        IF (SNOW-EVAP.gt.RUN0) THEN  ! only snow melts
          TRUN0(:)=RUN0*(TRSNOW(:)-TREVAP(:))/(SNOW-EVAP)
          TRSNOW(:)=TRSNOW(:)-TREVAP(:)-TRUN0(:)
        ELSE ! all snow + some ice melts 
          TRUN0(:)=TRSNOW(:)-TREVAP(:)+(ACE1LI-SNANDI)*TRLI(:)/(ACE2LI
     *         +ACE1LI)
          TRSNOW(:)=0.
c         TRLI(:)=TRLI(:)*(1.- ((ACE1LI-SNANDI)/(ACE2LI+ACE1LI)))
          do n=1,ntm
            TRUN0(N)=MAX(TRUN0(N),0d0)
          end do
        END IF
#endif
      ELSE
C**** FLUXES RECOMPUTE TG1 WHICH IS BELOW FREEZING POINT
        TG1=TG1+ENRG1/HC1
#ifdef TRACERS_WATER
        if (SNOW.gt.EVAP) THEN
          TRSNOW(:)=TRSNOW(:)-TREVAP(:)
        else
c         TRLI(:)=TRLI(:)-(TREVAP(:)-TRSNOW(:))
          TRSNOW(:)=0.
        end if
#endif
      END IF
      IF (SNANDI.LT.ACE1LI) THEN
C**** SOME ICE HAS MELTED OR EVAPORATED, TAKE IT FROM G2
        SNOW=0.
        DIFS=SNANDI-ACE1LI
        TG1=(TG1*SNANDI-TG2*DIFS)/ACE1LI
        EDIFS=DIFS*(TG2*SHI-LHM)
#ifdef TRACERS_WATER
c       TRLI(:)=TRLI(:)*(ACE2LI+ACE1LI)/(ACE2LI+ACE1LI-DIFS) 
        TRDIFS(:)=TRLI(:)*(DIFS/(ACE2LI+ACE1LI))
#endif
      ELSE
        SNOW=SNANDI-ACE1LI
      END IF
C**** CALCULATE TG2
      TG2=TG2+F1DT/HC2LI

      RETURN
      END SUBROUTINE LNDICE

      END MODULE LANDICE

      MODULE LANDICE_COM
!@sum  LANDICE_COM contains the model arrays for land ice
!@auth Gavin Schmidt
!@ver  1.0
!@cont io_landice
      USE MODEL_COM, only : im,jm
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm
#endif
      IMPLICIT NONE
      SAVE
!@var SNOWLI snow amount on land ice (kg/m^2)
      REAL*8, DIMENSION(IM,JM) :: SNOWLI
!@var TLANDI temperature of each land ice layer (C)
      REAL*8, DIMENSION(2,IM,JM) :: TLANDI

#ifdef TRACERS_WATER
!@var TRSNOWLI tracer amount in land ice snow (kg/m^2)
      REAL*8, DIMENSION(NTM,IM,JM) :: TRSNOWLI
!@var TRLNDI tracer amount in land ice (kg/m^2)
      REAL*8, DIMENSION(NTM,IM,JM) :: TRLNDI
!@var TRLI0 default tracer conc. for land ice (kg/kg)
      REAL*8, DIMENSION(NTM) :: TRLI0 
#endif

      END MODULE LANDICE_COM

      SUBROUTINE io_landice(kunit,iaction,ioerr)
!@sum  io_landice reads and writes landice variables to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,lhead,irsfic,irsficno,irerun
      USE LANDICE_COM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "GLAIC01"
#ifdef TRACERS_WATER
!@var TRHEADER Character string label for individual records
      CHARACTER*80 :: TRHEADER, TRMODULE_HEADER = "TRGLAC01"

      write (TRMODULE_HEADER(lhead+1:80)
     *     ,'(a7,i3,a)')'R8 dim(',NTM,',im,jm):TRSNOWLI,TRLNDI'
#endif

      MODULE_HEADER(lhead+1:80) = 'R8 SNOW(im,jm),T(2,im,jm)'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,SNOWLI,TLANDI
#ifdef TRACERS_WATER
        WRITE (kunit,err=10) TRMODULE_HEADER,TRSNOWLI,TRLNDI
#endif
      CASE (IOREAD:)            ! input from restart file
        READ (kunit,err=10) HEADER,SNOWLI,TLANDI
        IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
          PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
          GO TO 10
        END IF
#ifdef TRACERS_WATER
        SELECT CASE (IACTION)
        CASE (IRERUN,IOREAD,IRSFIC,IRSFICNO)    ! from reruns/restarts
          READ (kunit,err=10) TRHEADER,TRSNOWLI,TRLNDI
          IF (TRHEADER(1:LHEAD).NE.TRMODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",TRHEADER
     *           ,TRMODULE_HEADER
            GO TO 10
          END IF
        END SELECT
#endif
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_landice


