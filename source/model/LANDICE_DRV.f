#include "rundeck_opts.h"

!@sum  LANDICE_DRV contains drivers for LANDICE related routines
!@auth Gavin Schmidt
!@ver  1.0
!@cont init_LI,PRECIP_LI,GROUND_LI

      SUBROUTINE init_LI
!@sum  init_ice initialises landice arrays
!@auth Original Development Team
!@ver  1.0
      USE MODEL_COM, only : im,jm,flice
      USE LANDICE, only: ace1li,ace2li
      USE LANDICE_COM, only : tlandi,snowli
#ifdef TRACERS_WATER
     *     ,trsnowli,trlndi,trli0
#endif
      USE FLUXES, only : gtemp
#ifdef TRACERS_WATER
     *     ,gtracer
#endif
      USE DAGCOM, only : npts,icon_MLI,icon_HLI,title_con,conpt0
      IMPLICIT NONE
      LOGICAL :: QCON(NPTS), T=.TRUE. , F=.FALSE.
      INTEGER I,J

C**** set GTEMP array for landice
      DO J=1,JM
        DO I=1,IM
          IF (FLICE(I,J).gt.0) THEN
            GTEMP(1:2,3,I,J)=TLANDI(1:2,I,J)
#ifdef TRACERS_WATER
            IF (SNOWLI(I,J).gt.1d-5) THEN
              GTRACER(:,3,I,J)=TRSNOWLI(:,I,J)/SNOWLI(I,J)
            ELSE
              GTRACER(:,3,I,J)=trli0(:) !TRLNDI(:,I,J)/(ACE1LI+ACE2LI)
            END IF
#endif
          END IF
        END DO
      END DO

C**** Set conservation diagnostics for land ice mass, energy
      QCON=(/ F, F, F, T, T, F, F, F, T, F, F/)
      CALL SET_CON(QCON,CONPT0,"LAND ICE","(KG/M^2)        ",
     *     "(10**-9 KG/SM^2)",1d0,1d9,icon_MLI)
      QCON=(/ F, F, F, T, T, F, F, F, T, F, F/)
      CALL SET_CON(QCON,CONPT0,"LNDI ENR","(10**6 J/M^2)   ",
     *     "(10**-3 W/M^2)  ",1d-6,1d3,icon_HLI)
C****
      END SUBROUTINE init_LI

      SUBROUTINE PRECIP_LI
!@sum  PRECIP_LI driver for applying precipitation to land ice fraction
!@auth Original Development team
!@ver  1.0
!@calls LANDICE:PRECLI
      USE MODEL_COM, only : im,jm,flice,itlandi
      USE GEOM, only : imaxj,dxyp,bydxyp
      USE FLUXES, only : runoli,prec,eprec,gtemp
#ifdef TRACERS_WATER
     *     ,trunoli,trprec,gtracer
#endif
      USE LANDICE, only: ace1li,ace2li,precli
      USE LANDICE_COM, only : snowli,tlandi
#ifdef TRACERS_WATER
     *     ,trsnowli,trlndi,ntm,trli0
#endif
      USE DAGCOM, only : aj,areg,aij,jreg,ij_f0li,ij_f1li,ij_erun2
     *     ,ij_runli,j_run,j_implh,j_implm
      IMPLICIT NONE

      REAL*8 SNOW,TG1,TG2,PRCP,ENRGP,EDIFS,DIFS,ERUN2,RUN0,PLICE,DXYPJ
#ifdef TRACERS_WATER
!@var TRSNOW tracer amount in snow (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRSNOW
!@var TRLI tracer amount in land ice (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRLI
!@var TRPRCP tracer amount in precip (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRPRCP
!@var TRUN0 tracer runoff from ice (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRUN0
!@var TRDIFS implicit tracer flux at base of ice (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRDIFS
#endif
      INTEGER I,J,JR

      DO J=1,JM
      DXYPJ=DXYP(J)
      DO I=1,IMAXJ(J)
      PLICE=FLICE(I,J)
      PRCP=PREC(I,J)
      JR=JREG(I,J)
      RUNOLI(I,J)=0
#ifdef TRACERS_WATER
      TRUNOLI(:,I,J)=0.
#endif
      IF (PLICE.gt.0 .and. PRCP.gt.0) THEN

        ENRGP=EPREC(I,J)      ! energy of precipitation
        SNOW=SNOWLI(I,J)
        TG1=TLANDI(1,I,J)
        TG2=TLANDI(2,I,J)
#ifdef TRACERS_WATER
        TRLI(:)=TRLNDI(:,I,J)
        TRSNOW(:)=TRSNOWLI(:,I,J)
        TRPRCP(:)=TRPREC(:,I,J)*BYDXYP(J)
#endif
        AIJ(I,J,IJ_F0LI)=AIJ(I,J,IJ_F0LI)+ENRGP

        CALL PRECLI(SNOW,TG1,TG2,PRCP,ENRGP,
#ifdef TRACERS_WATER
     *       TRSNOW,TRLI,TRPRCP,TRDIFS,TRUN0,
#endif
     *       EDIFS,DIFS,ERUN2,RUN0)

C**** RESAVE PROGNOSTIC QUANTITIES AND FLUXES
        SNOWLI(I,J)=SNOW
        TLANDI(1,I,J)=TG1
        TLANDI(2,I,J)=TG2
        RUNOLI(I,J)  =RUN0
        GTEMP(1:2,3,I,J)=TLANDI(1:2,I,J)
#ifdef TRACERS_WATER
        TRLNDI(:,I,J)=TRLI(:)
        TRSNOWLI(:,I,J)=TRSNOW(:)
        TRUNOLI(:,I,J)=TRUN0(:)
c       TRDIFS(:)     !  diagnostic?
        IF (SNOW.gt.1d-5) THEN
          GTRACER(:,3,I,J)=TRSNOW(:)/SNOW
        ELSE
          GTRACER(:,3,I,J)=trli0(:)  !TRLI(:)/(ACE1LI+ACE2LI)
        END IF
#endif
C**** ACCUMULATE DIAGNOSTICS
c        AJ(J,J_TYPE, ITLANDI)=AJ(J,J_TYPE, ITLANDI)+      PLICE
        AJ(J,J_RUN,  ITLANDI)=AJ(J,J_RUN,  ITLANDI)+RUN0 *PLICE
C       AJ(J,J_ERUN ,ITLANDI)=AJ(J,J_ERUN ,ITLANDI)+ERUN0*PLICE ! (Tg=0)
        AJ(J,J_IMPLM,ITLANDI)=AJ(J,J_IMPLM,ITLANDI)+DIFS *PLICE
        AJ(J,J_IMPLH,ITLANDI)=AJ(J,J_IMPLH,ITLANDI)+ERUN2*PLICE
        AREG(JR,J_RUN  )=AREG(JR,J_RUN  )+RUN0 *PLICE*DXYPJ
c       AREG(JR,J_ERUN )=AREG(JR,J_ERUN )+ERUN0*PLICE*DXYPJ ! (Tg=0)
        AREG(JR,J_IMPLM)=AREG(JR,J_IMPLM)+DIFS *PLICE*DXYPJ
        AREG(JR,J_IMPLH)=AREG(JR,J_IMPLH)+ERUN2*PLICE*DXYPJ
        AIJ(I,J,IJ_F1LI) =AIJ(I,J,IJ_F1LI) +EDIFS
        AIJ(I,J,IJ_ERUN2)=AIJ(I,J,IJ_ERUN2)+ERUN2
        AIJ(I,J,IJ_RUNLI)=AIJ(I,J,IJ_RUNLI)+RUN0
      END IF
      END DO
      END DO
      END SUBROUTINE PRECIP_LI

      SUBROUTINE GROUND_LI
!@sum  GROUND_LI driver for applying surface fluxes to land ice fraction
!@auth Original Development team
!@ver  1.0
!@calls LANDICE:LNDICE
      USE MODEL_COM, only : im,jm,flice,itlandi
      USE GEOM, only : imaxj,dxyp
      USE DAGCOM, only : aj,areg,aij,jreg,ij_runli,ij_f1li,ij_erun2
     *     ,j_wtr1,j_ace1,j_wtr2,j_ace2,j_snow,j_run
     *     ,j_implh,j_implm,j_rsnow,ij_rsnw,ij_rsit,ij_snow
      USE LANDICE_COM, only : snowli,tlandi
#ifdef TRACERS_WATER
     *     ,ntm,trsnowli,trlndi,trli0
#endif
      USE LANDICE, only : lndice,ace1li,ace2li
      USE FLUXES, only : e0,e1,evapor,gtemp,runoli
#ifdef TRACERS_WATER
     *     ,trunoli,trevapor,gtracer
#endif
      IMPLICIT NONE

      REAL*8 SNOW,TG1,TG2,F0DT,F1DT,EVAP,EDIFS,DIFS,RUN0,PLICE,DXYPJ
     *     ,SCOVLI
#ifdef TRACERS_WATER
!@var TRSNOW tracer amount in snow (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRSNOW
!@var TRLI tracer amount in land ice (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRLI
!@var TREVAP tracer amount in evaporation (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TREVAP
!@var TRUN0 tracer runoff from ice (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRUN0
!@var TRDIFS implicit tracer flux at base of ice (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRDIFS
#endif
      INTEGER I,J,JR

      DO J=1,JM
      DXYPJ=DXYP(J)
      DO I=1,IMAXJ(J)
      PLICE=FLICE(I,J)
      JR=JREG(I,J)
      RUNOLI(I,J)=0.
      IF (PLICE.gt.0) THEN

        SNOW=SNOWLI(I,J)
        TG1=TLANDI(1,I,J)
        TG2=TLANDI(2,I,J)
        F0DT=E0(I,J,3)
        F1DT=E1(I,J,3)
        EVAP=EVAPOR(I,J,3)
#ifdef TRACERS_WATER
        TRLI(:)=TRLNDI(:,I,J)
        TRSNOW(:)=TRSNOWLI(:,I,J)
        TREVAP(:)=TREVAPOR(:,3,I,J)
#endif

        CALL LNDICE(SNOW,TG1,TG2,F0DT,F1DT,EVAP,
#ifdef TRACERS_WATER
     *     TRSNOW,TRLI,TREVAP,TRDIFS,TRUN0,
#endif
     *     EDIFS,DIFS,RUN0)

C**** RESAVE PROGNOSTIC QUANTITIES AND FLUXES
        SNOWLI(I,J)=SNOW
        TLANDI(1,I,J)=TG1
        TLANDI(2,I,J)=TG2
        RUNOLI(I,J) = RUN0
        GTEMP(1:2,3,I,J)=TLANDI(1:2,I,J)
#ifdef TRACERS_WATER
        TRLNDI(:,I,J)=TRLI(:)
        TRSNOWLI(:,I,J)=TRSNOW(:)
        TRUNOLI(:,I,J)=TRUN0(:)
c       TRDIFS(:)     !  diagnostic?
        IF (SNOW.gt.1d-5) THEN
          GTRACER(:,3,I,J)=TRSNOW(:)/SNOW
        ELSE
          GTRACER(:,3,I,J)=trli0(:) !TRLI(:)/(ACE1LI+ACE2LI)
        END IF
#endif
C**** ACCUMULATE DIAGNOSTICS
        SCOVLI=0
        IF (SNOWLI(I,J).GT.0.) SCOVLI=PLICE
        AJ(J,J_RSNOW,ITLANDI)=AJ(J,J_RSNOW,ITLANDI)+SCOVLI
        AREG(JR,J_RSNOW)=AREG(JR,J_RSNOW)+SCOVLI*DXYPJ
        AIJ(I,J,IJ_RSNW)=AIJ(I,J,IJ_RSNW)+SCOVLI
        AIJ(I,J,IJ_SNOW)=AIJ(I,J,IJ_SNOW)+SNOW*PLICE
        AIJ(I,J,IJ_RSIT)=AIJ(I,J,IJ_RSIT)+PLICE

        AJ(J,J_RUN,ITLANDI)  =AJ(J,J_RUN,ITLANDI)  +RUN0 *PLICE
        AJ(J,J_SNOW,ITLANDI) =AJ(J,J_SNOW,ITLANDI) +SNOW *PLICE
        AJ(J,J_ACE1,ITLANDI) =AJ(J,J_ACE1,ITLANDI)+ACE1LI*PLICE
        AJ(J,J_ACE2,ITLANDI) =AJ(J,J_ACE2,ITLANDI)+ACE2LI*PLICE
        AJ(J,J_IMPLH,ITLANDI)=AJ(J,J_IMPLH,ITLANDI)+EDIFS*PLICE
        AJ(J,J_IMPLM,ITLANDI)=AJ(J,J_IMPLM,ITLANDI)+DIFS *PLICE

        AREG(JR,J_RUN)  =AREG(JR,J_RUN)  +RUN0  *PLICE*DXYPJ
        AREG(JR,J_SNOW) =AREG(JR,J_SNOW) +SNOW  *PLICE*DXYPJ
        AREG(JR,J_ACE1) =AREG(JR,J_ACE1) +ACE1LI*PLICE*DXYPJ
        AREG(JR,J_ACE2) =AREG(JR,J_ACE2) +ACE2LI*PLICE*DXYPJ
        AREG(JR,J_IMPLH)=AREG(JR,J_IMPLH)+EDIFS *PLICE*DXYPJ
        AREG(JR,J_IMPLM)=AREG(JR,J_IMPLM)+DIFS  *PLICE*DXYPJ

        AIJ(I,J,IJ_F1LI) =AIJ(I,J,IJ_F1LI) +EDIFS+F1DT
        AIJ(I,J,IJ_RUNLI)=AIJ(I,J,IJ_RUNLI)+RUN0
        AIJ(I,J,IJ_ERUN2)=AIJ(I,J,IJ_ERUN2)+EDIFS
C****
      END IF
      END DO
      END DO
      END SUBROUTINE GROUND_LI

      SUBROUTINE conserv_MLI(ICE)
!@sum  conserv_MLI calculates total amount of snow and ice in land ice
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,jm,fim,flice
      USE GEOM, only : imaxj
      USE LANDICE_COM, only : snowli
      USE LANDICE, only : lndice,ace1li,ace2li
      IMPLICIT NONE
!@var ICE total land ice snow and ice mass (kg/m^2)
      REAL*8, DIMENSION(JM) :: ICE
      INTEGER I,J

      DO J=1,JM
        ICE(J)=0
        DO I=1,IMAXJ(J)
          ICE(J)=ICE(J)+FLICE(I,J)*(ACE1LI+ACE2LI+SNOWLI(I,J))
        END DO
      END DO
      ICE(1) =FIM*ICE(1)
      ICE(JM)=FIM*ICE(JM)
      RETURN
C****
      END SUBROUTINE conserv_MLI

      SUBROUTINE conserv_HLI(EICE)
!@sum  conserv_HLI calculates total land ice energy 
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : shi,lhm
      USE MODEL_COM, only : im,jm,fim,flice
      USE GEOM, only : imaxj
      USE LANDICE_COM, only : snowli,tlandi
      USE LANDICE, only : ace1li,ace2li
      IMPLICIT NONE
!@var EICE total land ice energy (J/m^2)
      REAL*8, DIMENSION(JM) :: EICE
      INTEGER I,J

      DO J=1,JM
        EICE(J)=0
        DO I=1,IMAXJ(J)
          EICE(J)=EICE(J)+FLICE(I,J)*((TLANDI(1,I,J)*SHI-LHM)*(ACE1LI
     *         +SNOWLI(I,J))+(TLANDI(2,I,J)*SHI-LHM)*ACE2LI)
        END DO
      END DO
      EICE(1) =FIM*EICE(1)
      EICE(JM)=FIM*EICE(JM)
      RETURN
C****
      END SUBROUTINE conserv_HLI
