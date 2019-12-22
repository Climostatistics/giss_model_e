#include "rundeck_opts.h"

      MODULE LAKES
!@sum  LAKES subroutines for Lakes and Rivers
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0 (based on LB265)
      USE CONSTANT, only : grav,bygrav,shw,rhow,lhm,shi,teeny
      USE MODEL_COM, only : im,jm
#ifdef TRACERS_WATER
      USE TRACER_COM, only : trname,ntm
#endif
      IMPLICIT NONE
      SAVE
C****
C**** Changes from Model III: MO -> MWL (kg), G0M -> GML (J),
C****                         GZM -> TLAKE (deg C)
!@var KDIREC directions for river flow
C**** (0 no flow, 1-8 anti-clockwise from top RH corner
      INTEGER, DIMENSION(IM,JM) :: KDIREC
!@var IDPOLE,JDPOLE special directions for south pole flow
      INTEGER :: IDPOLE,JDPOLE
!@var RATE rate of river flow downslope (fraction)
!@var DHORZ horizontal distance to downstream box (m)
      REAL*8, DIMENSION(IM,JM) :: RATE,DHORZ
!@var IFLOW,JFLOW grid box indexes for downstream direction
      INTEGER IFLOW(IM,JM),JFLOW(IM,JM)
!@param NRVRMX Max No. of specified rivers
      INTEGER, PARAMETER :: NRVRMX = 42 
!@var NRVR actual No. of specified rivers
      INTEGER :: NRVR
!@var IRVRMTH,JRVRMTH indexes for specified river mouths
      INTEGER, DIMENSION(NRVRMX) :: IRVRMTH,JRVRMTH
!@var NAMERVR Names of specified rivers
      CHARACTER*8, DIMENSION(NRVRMX) :: NAMERVR

!@param MINMLD minimum mixed layer depth in lake (m)
      REAL*8, PARAMETER :: MINMLD = 1.
!@param TMAXRHO temperature of maximum density (pure water) (C)
      REAL*8, PARAMETER :: TMAXRHO = 4.
!@param KVLAKE lake diffusion constant at mixed layer depth (m^2/s)
      REAL*8, PARAMETER :: KVLAKE = 1d-5
!@param TFL freezing temperature for lakes (=0 C)
      REAL*8, PARAMETER :: TFL = 0.
!@param AC1LMIN, AC2LMIN minimum ice thickness for lake ice (kg/m^2)
      REAL*8, PARAMETER :: AC1LMIN = 0.1, AC2LMIN=0.1  ! (not used)
!@param FLEADLK lead fraction for lakes
      REAL*8, PARAMETER :: FLEADLK = 0.
!@param BYZETA reciprocal of solar rad. extinction depth for lake (1/m)
      REAL*8, PARAMETER :: BYZETA = 1./0.35d0

      CONTAINS

      SUBROUTINE LKSOURC (I0,J0,ROICE,MLAKE,ELAKE,RUN0,FODT,FIDT,SROX
     *     ,FSR2,
#ifdef TRACERS_WATER
     *     TRLAKEL,TRUN0,TREVAP,TRO,TRI,
#endif
     *     EVAPO,ENRGFO,ACEFO,ACEFI,ENRGFI)
!@sum  LKSOURC applies fluxes to lake in ice-covered and ice-free areas
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : qcheck
      IMPLICIT NONE
!@var MLAKE,ELAKE mass and energy in lake layers (kg,J /m^2)
      REAL*8, INTENT(INOUT), DIMENSION(2) :: MLAKE,ELAKE
      INTEGER, INTENT(IN) :: I0,J0
      REAL*8, INTENT(IN) :: ROICE, EVAPO, RUN0
      REAL*8, INTENT(IN) :: FODT, FIDT, SROX(2)
      REAL*8, INTENT(OUT) :: ENRGFO, ACEFO, ENRGFI, ACEFI
#ifdef TRACERS_WATER
      REAL*8, INTENT(INOUT), DIMENSION(NTM,2) :: TRLAKEL
      REAL*8, INTENT(IN), DIMENSION(NTM) :: TRUN0,TREVAP
      REAL*8, INTENT(OUT), DIMENSION(NTM) :: TRO,TRI
      REAL*8, DIMENSION(NTM) :: DTR2,TRUNO,TRUNI,TRF1,TRF2,FRAC
#ifdef TRACERS_SPECIAL_O18
      REAL*8 fracls
      INTEGER N
#endif
#endif
!@var emin min energy deficit required before ice forms (J/m^2)
      REAL*8, PARAMETER :: emin=-1d-10 
      REAL*8 ENRGF1, ACEF1, ENRGF2, ACEF2, FHO, FHI, FH0, FH1, FH2, FSR2
      REAL*8 ENRGI, ENRGI2, ENRGO, ENRGO2, RUNO, RUNI, TLK2, DM2, DH2
      REAL*8 FRATO,FRATI,E2O,E2I
C**** initiallize output
      ENRGFO=0. ; ACEFO=0. ; ACEFI=0. ; ENRGFI=0.

C**** Calculate heat and mass fluxes to lake
      ENRGO = FODT-SROX(1)*FSR2 ! in open water
      ENRGO2=     +SROX(1)*FSR2 ! in open water, second layer
      ENRGI = FIDT-SROX(2)*FSR2 ! under ice
      ENRGI2=     +SROX(2)*FSR2 ! under ice, second layer
      RUNO  =-EVAPO
      RUNI  = RUN0
#ifdef TRACERS_WATER
      TRUNO(:)=-TREVAP(:)
      TRUNI(:)= TRUN0(:)
      FRAC(:)=1.
#ifdef TRACERS_SPECIAL_O18
      do n=1,ntm
        FRAC(n)=fracls(trname(n)) ! fractionation when freezing
      end do
#endif
#endif

C**** Bring up mass from second layer if required/allowed
      IF (MLAKE(1)+RUNO.lt.MINMLD*RHOW.and.MLAKE(2).gt.0) THEN
        DM2 = MIN(MLAKE(2),MINMLD*RHOW-(MLAKE(1)+RUNO))
        DH2 = DM2*(ELAKE(2)+(1.-ROICE)*ENRGO2+ROICE*ENRGI2)/MLAKE(2)
#ifdef TRACERS_WATER
        DTR2(:) = DM2*TRLAKEL(:,2)/MLAKE(2)
#endif
      ELSE
        DM2 = 0.
        DH2 = 0.
#ifdef TRACERS_WATER
        DTR2(:) = 0.
#endif
      END IF

C**** Apply fluxes to 2nd layer
      IF (DM2.lt.MLAKE(2)) THEN
        MLAKE(2)=MLAKE(2) - DM2
        ELAKE(2)=ELAKE(2) - DH2 + (1.-ROICE)*ENRGO2 + ROICE*ENRGI2
#ifdef TRACERS_WATER
        TRLAKEL(:,2)=TRLAKEL(:,2) - DTR2(:)
#endif
      ELSE
        MLAKE(2)=0.
        ELAKE(2)=0.
#ifdef TRACERS_WATER
        TRLAKEL(:,2)=0.
#endif
      END IF

      E2O = 0. ; E2I = 0.

C**** Calculate energy in mixed layer (open ocean)
      IF (ROICE.LT.1d0) THEN
        FHO=ELAKE(1)+ENRGO+DH2-(MLAKE(1)+DM2+RUNO)*TFL*SHW
        IF (FHO.LT.emin) THEN ! FLUXES COOL WATER TO FREEZING, FORM ICE
          ACEFO =FHO/(TFL*(SHI-SHW)-LHM)
          ACEFO =MIN(ACEFO,MAX(MLAKE(1)+DM2+RUNO-MINMLD*RHOW,0d0))
          ENRGFO=ACEFO*(TFL*SHI-LHM)
          E2O=FHO-ENRGFO
        END IF
      END IF

      IF (ROICE.GT.0) THEN
C**** Calculate energy in mixed layer (under ice)
        FHI=ELAKE(1)+DH2+ENRGI-(MLAKE(1)+DM2+RUNI)*TFL*SHW
        IF (FHI.LT.emin) THEN ! FLUXES COOL WATER TO FREEZING, FORM ICE
          ACEFI =FHI/(TFL*(SHI-SHW)-LHM)
          ACEFI =MIN(ACEFI,MAX(MLAKE(1)+DM2+RUNI-MINMLD*RHOW,0d0))
          ENRGFI=ACEFI*(TFL*SHI-LHM)
          E2I=FHI-ENRGFI
        END IF
      END IF
#ifdef TRACERS_WATER
      TRO(:)=ACEFO*FRAC(:)*TRLAKEL(:,1)/MLAKE(1)
      TRI(:)=ACEFI*FRAC(:)*TRLAKEL(:,1)/MLAKE(1)
#endif

C**** Update first layer variables
      MLAKE(1)=MLAKE(1)+DM2+(1.-ROICE)*(RUNO -ACEFO)+ROICE*(RUNI-ACEFI)
      ELAKE(1)=ELAKE(1)+DH2+(1.-ROICE)*(ENRGO-ENRGFO)+
     *                                             ROICE*(ENRGI-ENRGFI)
#ifdef TRACERS_WATER
      TRLAKEL(:,1)=TRLAKEL(:,1)+DTR2(:)+(1.-ROICE)*(TRUNO(:)-TRO(:))+
     *     ROICE*(TRUNI(:)-TRI(:))
#endif

      ACEF1=0. ; ACEF2=0. ; ENRGF1=0. ; ENRGF2=0.
C**** Take remaining energy and start to freeze second layer
      FH2= ELAKE(1)-MLAKE(1)*TFL*SHW
      IF (FH2.LT.emin) THEN
        IF (MLAKE(2).gt.0) THEN
C**** FH2=-ACEF2*(TLK2-TFL)*SHW+ACEF2*LHM
          TLK2    =ELAKE(2)/(MLAKE(2)*SHW)
          ACEF2   =-FH2/(TLK2*SHW-TFL*SHI+LHM)
          ACEF2   =MIN(ACEF2,MLAKE(2))
          ENRGF2  =ACEF2*(TFL*SHI-LHM)
          ELAKE(1)=ELAKE(1)+ACEF2*TLK2*SHW-ENRGF2
          ELAKE(2)=ELAKE(2)-ACEF2*TLK2*SHW
          MLAKE(2)=MLAKE(2)-ACEF2
        END IF
        FH1= ELAKE(1)-MLAKE(1)*TFL*SHW
        IF (FH1.lt.emin) THEN      ! all layer 2 froze, freeze layer 1
          ACEF1   =FH1/(TFL*(SHI-SHW)-LHM)
C**** limit freezing if lake is between 50 and 20cm depth
          IF (MLAKE(1).lt.0.5d0*RHOW) THEN
            ACEF1=MIN(ACEF1,MAX(0.5*(MLAKE(1)-0.2d0*RHOW),0d0))
            if (qcheck) print*,"Lake freezing limited",ACEF1/RHOW
     *           ,MLAKE(1)/RHOW
          END IF
          ENRGF1  =ACEF1*(TFL*SHI-LHM)
          ELAKE(1)=ELAKE(1)-ENRGF1
          MLAKE(1)=MLAKE(1)-ACEF1
          FH0     =ELAKE(1)-MLAKE(1)*TFL*SHW
          IF (FH0.lt.-1d-8) THEN ! max. amount of lake frozen, cool ice
            if (qcheck) WRITE(6,*)
     *           "Minimum lake level reached: rsi,mlake,elake",i0,j0
     *           ,roice,mlake(1)/rhow,elake(1)
            ENRGF1  =ENRGF1+FH0
            ELAKE(1)=MLAKE(1)*TFL*SHW
          END IF
        END IF
      END IF
#ifdef TRACERS_WATER
      TRF1(:) = ACEF1*FRAC(:)*TRLAKEL(:,1)/(MLAKE(1)+ACEF1)
      TRLAKEL(:,1)=TRLAKEL(:,1)-TRF1(:)
      IF (MLAKE(2).gt.0) THEN
        TRF2(:) = ACEF2*FRAC(:)*TRLAKEL(:,2)/(MLAKE(2)+ACEF2)
        TRLAKEL(:,2)=TRLAKEL(:,2)-TRF2(:)
      ELSE ! possibility of complete freezing (and so no frac)
        TRF2(:) = TRLAKEL(:,2)
        TRLAKEL(:,2) = 0.
      END IF
#endif

C**** combine mass and energy fluxes for output
C**** Note that output fluxes are over open water/ice covered fractions
C**** distribute ice fluxes according to flux amounts
      FRATO = 1d0
      FRATI = 1d0
      IF (E2I+E2O.lt.0) THEN
        FRATO = E2O/(E2I*ROICE+E2O*(1.-ROICE))
        FRATI = E2I/(E2I*ROICE+E2O*(1.-ROICE))
      END IF
      ACEFO =ACEFO + (ACEF1 +ACEF2 )*FRATO
      ACEFI =ACEFI + (ACEF1 +ACEF2 )*FRATI
      ENRGFO=ENRGFO+ (ENRGF1+ENRGF2)*FRATO
      ENRGFI=ENRGFI+ (ENRGF1+ENRGF2)*FRATI
#ifdef TRACERS_WATER
      TRO(:)=TRO(:) + (TRF1(:) + TRF2(:))* FRATO
      TRI(:)=TRI(:) + (TRF1(:) + TRF2(:))* FRATI
#endif

      RETURN
      END SUBROUTINE LKSOURC

      SUBROUTINE LKMIX(MLAKE,ELAKE,
#ifdef TRACERS_WATER
     *     TRLAKEL,
#endif
     *     HLAKE,TKE,ROICE,DTSRC)
!@sum  LKMIX calculates mixing and entrainment in lakes
!@auth Gavin Schmidt
!@ver  1.0
      IMPLICIT NONE
!@var MLAKE,ELAKE mass and energy in lake layers (kg,J /m^2)
      REAL*8, INTENT(INOUT), DIMENSION(2) :: MLAKE,ELAKE
!@var TKE turbulent kinetic energy input at surface of lake (J/m^2)
!@var ROICE ice fraction
      REAL*8, INTENT(IN) :: TKE,ROICE
!@var HLAKE sill depth for lake (m)
      REAL*8, INTENT(IN) :: HLAKE
!@var DTSRC source time step (s)
      REAL*8, INTENT(IN) :: DTSRC
#ifdef TRACERS_WATER
!@var TRLAKEL tracer mass in lake layers (kg/m^2)
      REAL*8, INTENT(INOUT), DIMENSION(NTM,2) :: TRLAKEL
      REAL*8, DIMENSION(NTM) :: DTML,TR1N,TR2N,TRLT
#endif
!@param MAXRHO,RHO0,BFAC freshwater density function approximation
      REAL*8, PARAMETER :: MAXRHO=1d3, RHO0=999.842594d0,
     *     BFAC=(MAXRHO-RHO0)/16d0

      REAL*8 TLK1, TLK2, HLT, MLT, DTK, E1N, E2N, ATKE, H1, H2,
     *      DRHO, DML, DHML

C**** Only mix if there is a second layer!
      IF (MLAKE(2).gt.0) THEN
        TLK1=ELAKE(1)/(MLAKE(1)*SHW)
        TLK2=ELAKE(2)/(MLAKE(2)*SHW)
        HLT=ELAKE(1)+ELAKE(2)
        MLT=MLAKE(1)+MLAKE(2)
#ifdef TRACERS_WATER
        TRLT(:)=TRLAKEL(:,1)+TRLAKEL(:,2)
#endif
C**** Test for static stability
        IF ((TMAXRHO-TLK1)*(TLK2-TLK1).lt.0) THEN
C**** mix uniformly and set MLD to minimum
          MLAKE(1)=MAX(MINMLD*RHOW,(MLAKE(1)+MLAKE(2))-HLAKE*RHOW)
          MLAKE(2)=MLT-MLAKE(1)
          ELAKE(1)=HLT*MLAKE(1)/MLT
          ELAKE(2)=HLT*MLAKE(2)/MLT
#ifdef TRACERS_WATER
          TRLAKEL(:,1)=TRLT(:)*MLAKE(1)/MLT
          TRLAKEL(:,2)=TRLT(:)*MLAKE(2)/MLT
#endif
        ELSE ! not unstable, implicitly diffuse heat + entrain
C**** reduce mixing if there is ice cover, no mixing if temperature
C**** gradient is negative and deep water is at max density.
          IF (TLK2.lt.TLK1 .and. TLK2.gt.TMAXRHO) THEN
            DTK=2.*KVLAKE*(1.-ROICE)*DTSRC*RHOW**2
            E1N=(ELAKE(1)+DTK*HLT/(MLT*MLAKE(2)))/
     *           (1.+DTK/(MLAKE(1)*MLAKE(2)))
            E2N=(ELAKE(2)+DTK*HLT/(MLT*MLAKE(1)))/
     *           (1.+DTK/(MLAKE(1)*MLAKE(2)))
            ELAKE(1)=E1N
            ELAKE(2)=E2N
#ifdef TRACERS_WATER
C**** diffuse tracers using same KV as for heat?
            TR1N(:)=(TRLAKEL(:,1)+DTK*TRLT(:)/(MLT*MLAKE(2)))/
     *           (1.+DTK/(MLAKE(1)*MLAKE(2)))
            TR2N(:)=(TRLAKEL(:,2)+DTK*TRLT(:)/(MLT*MLAKE(1)))/
     *           (1.+DTK/(MLAKE(1)*MLAKE(2)))
            TRLAKEL(:,1)=TR1N(:)
            TRLAKEL(:,2)=TR2N(:)
#endif
          END IF
C**** entrain deep water if there is available TKE
C**** take a factor of TKE and calculate change in PE
          IF (TKE.gt.0) THEN
            ATKE=0.2d0*TKE      ! 20% of TKE is available for mixing
            H1=MLAKE(1)/RHOW
            H2=MLAKE(2)/RHOW
C**** DRHO=RHO(TLK2)-RHO(TLK1)~=(TLK2-TLK1)*dRHOdT(TLK1)
C**** Assumes a parabolic density function going through MAXRHO at
C**** TMAXRHO, and RHO0 at T=0. (reasonable up to about 12 C)
            DRHO=(TLK2-TLK1)*2d0*BFAC*(TMAXRHO-TLK1)
            DML=ATKE*BYGRAV/(DRHO*0.5*H1)
            IF (DML*RHOW.lt.MLAKE(2)) THEN
              DHML=DML*ELAKE(2)/H2
              ELAKE(1)=ELAKE(1)+DHML
              ELAKE(2)=ELAKE(2)-DHML
              MLAKE(1)=MLAKE(1)+DML*RHOW
              MLAKE(2)=MLAKE(2)-DML*RHOW
#ifdef TRACERS_WATER
              DTML(:)=DML*TRLAKEL(:,2)/H2
              TRLAKEL(:,1)=TRLAKEL(:,1)+DTML(:)
              TRLAKEL(:,2)=TRLAKEL(:,2)-DTML(:)
#endif
            ELSE                ! entire second layer is entrained
              MLAKE(1)=MLT
              MLAKE(2)=0
              ELAKE(1)=HLT
              ELAKE(2)=0
#ifdef TRACERS_WATER
              TRLAKEL(:,1)=TRLT(:)
              TRLAKEL(:,2)=0.
#endif
            END IF
          END IF
        END IF
      END IF
      RETURN
C****
      END SUBROUTINE LKMIX

      END MODULE LAKES

      SUBROUTINE init_LAKES(inilake)
!@sum  init_LAKES initiallises lake variables
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE FILEMANAGER
      USE CONSTANT, only : rhow,shw,tf,pi,grav
      USE MODEL_COM, only : im,jm,flake0,zatmo,dtsrc,flice,hlake
     *     ,focean,fearth,fland,jday
      USE GEOM, only : dxyp,dxv,dyv,dxp,dyp,imaxj
#ifdef TRACERS_WATER
      USE TRACER_COM, only : trw0
      USE FLUXES, only : gtracer
#endif
      USE FLUXES, only : gtemp,mlhc
      USE SEAICE_COM, only : rsi
      USE PBLCOM, only : tsavg
      USE LAKES
      USE LAKES_COM
      USE DAGCOM, only : npts,icon_LKM,icon_LKE,title_con,conpt0
      IMPLICIT NONE

      LOGICAL inilake
!@var I,J,I72,IU,JU,ID,JD loop variables
      INTEGER I,J,I72,IU,JU,ID,JD,INM,KD
      INTEGER iu_RVR  !@var iu_RVR unit number for river direction file
      CHARACTER TITLEI*80, CDIREC(IM,JM)*1, CONPT(NPTS)*10
      REAL*8 SPMIN,SPMAX,SPEED0,SPEED,DZDH,DZDH1,MLK1
      LOGICAL :: QCON(NPTS), T=.TRUE. , F=.FALSE.
C****
C**** LAKECB  MWL      Mass of water in lake (kg)
C****         GML      Liquid lake enthalpy (J)
C****         TLAKE    Temperature of lake surface (C)
C****         HLAKE    Lake sill depth (m)
C****         ALAKE    Lake surface area (m^2)
C****         TANLK    Lake slope (tan(alpha)) (1)
C****
C**** FIXDCB  FLAKE0    Original Lake fraction (1)
C****
      FLAKE = FLAKE0   ! this is only here until FLAKE is variable

C**** Ensure that HLAKE is a minimum of 1m for FLAKE>0
      DO J=1,JM
        DO I=1,IM
          IF (FLAKE(I,J).gt.0 .and. HLAKE(I,J).lt.1.) THEN
            print*,"Warning: Fixing HLAKE",i,j,FLAKE(I,J),HLAKE(I,J)
     *           ,"--> 1m"
            HLAKE(I,J)=1.
          END IF
        END DO
      END DO

      IF (INILAKE) THEN
C**** Set lake variables from surface temperature
C**** This is just an estimate for the initiallisation
        DO J=1,JM
          DO I=1,IM
c            FLAKE(I,J) = FLAKE0(I,J)
            IF (FLAKE(I,J).gt.0) THEN
              IF (HLAKE(I,J).lt.1.) print*,
     *             "HLAKE too small for lake fraction",i,j,HLAKE(I,J)
     *             ,FLAKE(I,J)
              TLAKE(I,J) = MAX(0d0,TSAVG(I,J)-TF)
              MWL(I,J)   = RHOW*HLAKE(I,J)*FLAKE(I,J)*DXYP(J)
              MLK1       = MINMLD*RHOW*FLAKE(I,J)*DXYP(J)
              GML(I,J)   = SHW*(MLK1*TLAKE(I,J)
     *             +(MWL(I,J)-MLK1)*MAX(TLAKE(I,J),4d0))
              MLDLK(I,J) = MINMLD
#ifdef TRACERS_WATER
              TRLAKE(:,1,I,J)=MLK1*TRW0(:)
              TRLAKE(:,2,I,J)=(MWL(I,J)-MLK1)*TRW0(:)
#endif
            ELSE
              TLAKE(I,J) = 0.
              MWL(I,J)   = 0.
              GML(I,J)   = 0.
              MLDLK(I,J) = MINMLD
#ifdef TRACERS_WATER
              TRLAKE(:,:,I,J)=0.
#endif
            END IF
          END DO
        END DO
      END IF

C**** Set geometric variables
C**** TANLK=TAN(ALPHA) = R/H for a conical lake of equivalent volume
      DO J=1,JM
        DO I=1,IM
          IF (FLAKE0(I,J).gt.0) THEN
            TANLK(I,J) = SQRT(FLAKE0(I,J)*DXYP(J)/PI)/(3d0*HLAKE(I,J))
          ELSE
            TANLK(I,J) = 2d3   ! reasonable average value
          END IF
        END DO
      END DO

      CALL PRINTLK("IN")
C**** Set GTEMP arrays for lakes
      DO J=1,JM
        DO I=1,IM
          IF (FLAKE(I,J).gt.0) THEN
            GTEMP(1,1,I,J)=TLAKE(I,J)
            IF (MWL(I,J).gt.(1d-10+MLDLK(I,J))*RHOW*FLAKE(I,J)*DXYP(J))
     *           THEN
              GTEMP(2,1,I,J)=(GML(I,J)-TLAKE(I,J)*SHW*MLDLK(I,J)*RHOW
     *             *FLAKE(I,J)*DXYP(J))/(SHW*(MWL(I,J)-MLDLK(I,J)
     *             *RHOW*FLAKE(I,J)*DXYP(J)))
            ELSE
              GTEMP(2,1,I,J)=TLAKE(I,J)
            END IF
#ifdef TRACERS_WATER
            GTRACER(:,1,I,J)=TRLAKE(:,1,I,J)/(MLDLK(I,J)*RHOW*FLAKE(I,J)
     *           *DXYP(J))
#endif
            MLHC(I,J)= SHW*MLDLK(I,J)*RHOW
          END IF
        END DO
      END DO

C****
C**** Always initiallise River direction and Rate
C**** Read in CDIREC: Number = octant direction, Letter = river mouth
      call openunit("RVR",iu_RVR,.false.,.true.)
      READ  (iu_RVR,910) TITLEI
      WRITE (6,*) 'River Direction file read: ',TITLEI
      READ  (iu_RVR,910)
      DO I72=1,1+(IM-1)/72
        DO J=JM,1,-1
          READ  (iu_RVR,911) (CDIREC(I,J),I=72*(I72-1)+1,MIN(IM,I72*72))
        END DO
      END DO
C**** read in named rivers (if any)
      READ (iu_RVR,*,END=10)
      READ (iu_RVR,'(A80)',END=10) TITLEI
      READ (iu_RVR,*,END=10)
      IF (TITLEI.eq."Named River Mouths:") THEN
        DO I=1,NRVRMX,5
          READ(iu_RVR,'(5(A8,X))') NAMERVR(I:MIN(NRVRMX,I+4))
        END DO
      END IF
 10   call closeunit (iu_RVR)
C**** Create integral direction array KDIREC from CDIREC
      INM=0
      DO J=1,JM
      DO I=1,IM
C**** KD: -16 = blank, 0-8 directions >8 named rivers
        KD= ICHAR(CDIREC(I,J)) - 48
C**** If land but no ocean, and no direction, print warning
        IF (FLAND(I,J).gt.0 .and. FOCEAN(I,J).le.0 .and.
     *       (KD.gt.8 .or. KD.lt.0)) THEN
          WRITE(6,*) "Land box has no river direction I,J: ",I,J
     *     ,FOCEAN(I,J),FLAND(I,J),FLICE(I,J),FLAKE0(I,J),FEARTH(I,J)
        END IF
C**** Default direction is down (if ocean box), or no outlet (if not)
C**** Also ensure that all ocean boxes are done properly
        IF ((KD.lt.0 .or. KD.gt.8) .or. FOCEAN(I,J).eq.1.) THEN
          KDIREC(I,J)=0.
        ELSE
          KDIREC(I,J) = KD
        END IF
C**** Check for specified river mouths
        IF (KD.GE.17 .AND. KD.LE.42) THEN
          INM=INM+1
          IRVRMTH(INM)=I
          JRVRMTH(INM)=J
          IF (CDIREC(I,J).ne.NAMERVR(INM)(1:1)) THEN
            WRITE(6,*) "Warning: Named river in RVR does not correspond"
     *           //" with letter in direction file. Please check"
            WRITE(6,*) "INM, CDIREC, NAMERVR = ",INM,CDIREC(I,J)
     *           ," ",NAMERVR(INM)
            NAMERVR(INM)=CDIREC(I,J)  ! set default
          END IF
          IF (FOCEAN(I,J).le.0) THEN
            WRITE(6,*) "Warning: Named river outlet must be in ocean",i
     *           ,j,NAMERVR(INM),FOCEAN(I,J),FLAND(I,J),FLICE(I,J)
     *           ,FLAKE0(I,J),FEARTH(I,J)
          END IF
        END IF
      END DO
      END DO
      NRVR=INM
C****
C**** From each box calculate the downstream river box
C****
      DO J=2,JM-1
        DO I=1,IM
          SELECT CASE (KDIREC(I,J))
          CASE (0)
            IFLOW(I,J) = I
            JFLOW(I,J) = J
            DHORZ(I,J) = 0.5*SQRT(DXP(J)*DXP(J)+DYP(J)*DYP(J))
          CASE (1)
            IFLOW(I,J) = I+1
            JFLOW(I,J) = J+1
            DHORZ(I,J) = SQRT(DXV(J+1)*DXV(J+1)+DYV(J+1)*DYV(J+1))
                                ! SQRT(DXV(J)*DXV(J)+DYV(J)*DYV(J))
            IF(I.eq.IM)  IFLOW(I,J) = 1
          CASE (2)
            IFLOW(I,J) = I
            JFLOW(I,J) = J+1
            DHORZ(I,J) = DYV(J+1)   ! DYV(J)
          CASE (3)
            IFLOW(I,J) = I-1
            JFLOW(I,J) = J+1
            DHORZ(I,J) = SQRT(DXV(J+1)*DXV(J+1)+DYV(J+1)*DYV(J+1))
                                ! SQRT(DXV(J)*DXV(J)+DYV(J)*DYV(J))
            IF(I.eq.1)  IFLOW(I,J) = IM
          CASE (4)
            IFLOW(I,J) = I-1
            JFLOW(I,J) = J
            DHORZ(I,J) = DXP(J)
            IF(I.eq.1)  IFLOW(I,J) = IM
          CASE (5)
            IFLOW(I,J) = I-1
            JFLOW(I,J) = J-1
            DHORZ(I,J) = SQRT(DXV(J)*DXV(J)+DYV(J)*DYV(J))
                            ! SQRT(DXV(J-1)*DXV(J-1)+DYV(J-1)*DYV(J-1))
            IF(I.eq.1)  IFLOW(I,J) = IM
          CASE (6)
            IFLOW(I,J) = I
            JFLOW(I,J) = J-1
            DHORZ(I,J) = DYV(J)    ! DYV(J-1)
          CASE (7)
            IFLOW(I,J) = I+1
            JFLOW(I,J) = J-1
            DHORZ(I,J) = SQRT(DXV(J)*DXV(J)+DYV(J)*DYV(J))
                            ! SQRT(DXV(J-1)*DXV(J-1)+DYV(J-1)*DYV(J-1))
            IF(I.eq.IM)  IFLOW(I,J) = 1
          CASE (8)
            IFLOW(I,J) = I+1
            JFLOW(I,J) = J
            DHORZ(I,J) = DXP(J)
            IF(I.eq.IM)  IFLOW(I,J) = 1
          END SELECT
        END DO
      END DO
C**** South Pole is a special case
      DO I=1,IM
        IF(KDIREC(I,1).eq.2)  THEN
          IDPOLE     = I
          JDPOLE     = 2
          IFLOW(1,1) = I
          JFLOW(1,1) = 2
          DHORZ(1,1) = DYV(2)    ! DYV(1)
        END IF
        IF(KDIREC(I,2).eq.6)  THEN
          IFLOW(I,2) = 1
          JFLOW(I,2) = 1
        END IF
      END DO
C****
C**** Calculate river flow RATE (per source time step)
C****
      SPEED0= .35d0
      SPMIN = .15d0
      SPMAX = 5.
      DZDH1 = .00005
      DO JU=1,JM-1
        DO IU=1,IMAXJ(JU)
          IF(KDIREC(IU,JU).gt.0) THEN
            JD=JFLOW(IU,JU)
            ID=IFLOW(IU,JU)
            DZDH  = (ZATMO(IU,JU)-ZATMO(ID,JD)) / (GRAV*DHORZ(IU,JU))
          ELSE
            DZDH  = ZATMO(IU,JU) / (GRAV*DHORZ(IU,JU))
          END IF
          SPEED = SPEED0*DZDH/DZDH1
          IF(SPEED.lt.SPMIN)  SPEED = SPMIN
          IF(SPEED.gt.SPMAX)  SPEED = SPMAX
          RATE(IU,JU) = DTsrc*SPEED/DHORZ(IU,JU)
        END DO
      END DO

C**** Set conservation diagnostics for Lake mass and energy
      CONPT=CONPT0
      CONPT(4)="PREC+LAT M"
      CONPT(5)="SURFACE"   ; CONPT(8)="RIVERS"
      QCON=(/ F, F, F, T, T, F, F, T, T, F, F/)
      CALL SET_CON(QCON,CONPT,"LAK MASS","(KG/M^2)      ",
     *     "(10**-9 KG/SM^2)",1d0,1d9,icon_LKM)
      QCON=(/ F, F, F, T, T, F, F, T, T, F, F/)
      CALL SET_CON(QCON,CONPT,"LAK ENRG","(10**3 J/M^2) ",
     *     "(10**-3 W/M^2)  ",1d-3,1d3,icon_LKE)

      RETURN
C****
 910  FORMAT (A72)
 911  FORMAT (72A1)
      END SUBROUTINE init_LAKES

      SUBROUTINE RIVERF
!@sum  RIVERF transports lake water from each grid box downstream
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0 (based on LB265)
      USE CONSTANT, only : grav,shw,rhow,teeny
      USE MODEL_COM, only : im,jm,focean,zatmo,hlake,itlake,itlkice
     *     ,itocean,itoice,fland
      USE GEOM, only : dxyp,bydxyp
      USE DAGCOM, only : aij,ij_ervr,ij_mrvr,ij_f0oc,aj,areg,jreg,j_rvrd
     *     ,j_ervr
#ifdef TRACERS_WATER
      USE TRACER_DIAG_COM, only : taijn,tij_rvr
      USE FLUXES, only : trflowo,gtracer
#endif
      USE FLUXES, only : flowo,eflowo,gtemp,mlhc
      USE LAKES, only : kdirec,idpole,jdpole,rate,iflow,jflow
      USE LAKES_COM, only : tlake,gml,mwl,mldlk,flake
#ifdef TRACERS_WATER
     *     ,trlake,ntm
#endif
      USE SEAICE_COM, only : rsi
      IMPLICIT NONE

!@var I,J,IU,JU,ID,JD loop variables
      INTEGER I,J,IU,JU,ID,JD,JR,ITYPE
      REAL*8 MWLSILL,DMM,DGM,HLK1,DPE
      REAL*8, DIMENSION(IM,JM) :: FLOW,EFLOW
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(NTM) :: DTM
      REAL*8, DIMENSION(NTM,IM,JM) :: TRFLOW
#endif
C****
C**** LAKECB  MWL  Liquid lake mass  (kg)
C****         GML  Liquid lake enthalpy  (J)
C****         TLAKE  Lake surface temperature (C)
C****
C**** Calculate net mass and energy changes due to river flow
C****
      FLOW = 0. ; EFLOW = 0.
      FLOWO = 0. ; EFLOWO = 0.
#ifdef TRACERS_WATER
      TRFLOW = 0.
      TRFLOWO = 0.
#endif
      DO JU=2,JM-1
        DO IU=1,IM
C**** Also allow flow into ocean fraction of same box if KDIREC=0
          IF (KDIREC(IU,JU).gt.0 .or.
     *       (FLAND(IU,JU).gt.0 .and. FOCEAN(IU,JU).gt.0))  THEN
            JD=JFLOW(IU,JU)
            ID=IFLOW(IU,JU)
C**** Only overflow if lake mass is above sill height (HLAKE (m))
            MWLSILL = RHOW*HLAKE(IU,JU)*FLAKE(IU,JU)*DXYP(JU)
            IF(MWL(IU,JU).gt.MWLSILL) THEN
              DMM = (MWL(IU,JU)-MWLSILL)*RATE(IU,JU)
              IF (MWL(IU,JU)-DMM.lt.1d-20) DMM=MWL(IU,JU)
c              IF (FLAKE(IU,JU).gt.0) THEN
c                MLM=RHOW*MLDLK(IU,JU)*FLAKE(IU,JU)*DXYP(JU)
c                DMM=MIN(DMM,MLM)   ! not necessary since MLM>TOTD-HLAKE
c              END IF
              DGM=TLAKE(IU,JU)*DMM*SHW  ! TLAKE always defined
              FLOW(IU,JU) =  FLOW(IU,JU) - DMM
              EFLOW(IU,JU) = EFLOW(IU,JU) - DGM
#ifdef TRACERS_WATER
              if (flake(iu,ju).gt.0) then
                DTM(:) = DMM*GTRACER(:,1,IU,JU)
              else
                DTM(:) = DMM*TRLAKE(:,1,IU,JU)/MWL(IU,JU)
              end if
              TRFLOW(:,IU,JU) = TRFLOW(:,IU,JU) - DTM(:)
#endif
              IF(FOCEAN(ID,JD).le.0.) THEN
                DPE=0.  ! DMM*(ZATMO(IU,JU)-ZATMO(ID,JD))
                FLOW(ID,JD) =  FLOW(ID,JD) + DMM
                EFLOW(ID,JD) = EFLOW(ID,JD) + DGM+DPE
#ifdef TRACERS_WATER
                TRFLOW(:,ID,JD) = TRFLOW(:,ID,JD) + DTM(:)
#endif
              ELSE ! Save river mouth flow to for output to oceans
C**** DPE: also add potential energy change to ocean.
C**** Normally ocean is at sea level (Duh!), but in some boxes ZATMO
C**** may not be zero if there is land as well, while in the Caspian,
C**** the ocean level is below zero.
C**** Note: this is diasabled until PE of precip is properly calculated
C**** in atmosphere as well. Otherwise, there is an energy imbalance.
                DPE=0.  ! DMM*(ZATMO(IU,JU)-MIN(0d0,ZATMO(ID,JD)))
                FLOWO(ID,JD)=FLOWO(ID,JD)+DMM
                EFLOWO(ID,JD)=EFLOWO(ID,JD)+DGM+DPE
#ifdef TRACERS_WATER
                TRFLOWO(:,ID,JD) = TRFLOWO(:,ID,JD) + DTM(:)
#endif
C**** accumulate river runoff diags (moved from ground)
                AJ(JD,J_RVRD,ITOCEAN)=AJ(JD,J_RVRD,ITOCEAN)+
     *               (1.-RSI(ID,JD))*DMM*BYDXYP(JD)
                AJ(JD,J_ERVR,ITOCEAN)=AJ(JD,J_ERVR,ITOCEAN)+
     *               (1.-RSI(ID,JD))*(DGM+DPE)*BYDXYP(JD)
                AJ(JD,J_RVRD,ITOICE)=AJ(JD,J_RVRD,ITOICE) +
     *               RSI(ID,JD)*DMM*BYDXYP(JD)
                AJ(JD,J_ERVR,ITOICE)=AJ(JD,J_ERVR,ITOICE) +
     *               RSI(ID,JD)*(DGM+DPE)*BYDXYP(JD)
                AIJ(ID,JD,IJ_F0OC)=AIJ(ID,JD,IJ_F0OC)+
     *               (1.-RSI(ID,JD))*(DGM+DPE)*BYDXYP(JD)
              END IF
              JR=JREG(ID,JD)
              AREG(JR,J_RVRD)=AREG(JR,J_RVRD)+DMM
              AREG(JR,J_ERVR)=AREG(JR,J_ERVR)+DGM+DPE
              AIJ(ID,JD,IJ_MRVR)=AIJ(ID,JD,IJ_MRVR) + DMM
              AIJ(ID,JD,IJ_ERVR)=AIJ(ID,JD,IJ_ERVR) + DGM+DPE
#ifdef TRACERS_WATER
              TAIJN(ID,JD,TIJ_RVR,:)=TAIJN(ID,JD,TIJ_RVR,:)+DTM(:)
     *             *BYDXYP(JD)
#endif
            END IF
          END IF
        END DO
      END DO
C****
C**** Calculate river flow at the South Pole
C****
       FLOW(1,1) =  FLOW(1,1)/IM
      EFLOW(1,1) = EFLOW(1,1)/IM
#ifdef TRACERS_WATER
      TRFLOW(:,1,1) = TRFLOW(:,1,1)/IM
#endif
C**** Only overflow if lake mass is above sill height (HLAKE (m))
      MWLSILL = RHOW*HLAKE(1,1)*FLAKE(1,1)*DXYP(1)
      IF(MWL(1,1).gt.MWLSILL) THEN
        DMM = (MWL(1,1)-MWLSILL)*RATE(1,1)
c        IF (FLAKE(1,1).gt.0) THEN
c         MLM=RHOW*MLDLK(1,1)*FLAKE(1,1)*DXYP(1)
c         DMM=MIN(DMM,MLM)   ! not necessary since MLM>TOTD-HLAKE
c        END IF
        DGM=TLAKE(1,1)*DMM*SHW ! TLAKE always defined
        DPE=0.  ! DMM*(ZATMO(1,1)-ZATMO(IDPOLE,JDPOLE))
        FLOW(1,1) =  FLOW(1,1) - DMM
        EFLOW(1,1) = EFLOW(1,1) - DGM
#ifdef TRACERS_WATER
        if (flake(iu,ju).gt.0) then
          DTM(:) = DMM*GTRACER(:,1,1,1)
        else
          DTM(:) = DMM*TRLAKE(:,1,1,1)/(MWL(1,1)+teeny)
        end if
        TRFLOW(:,1,1) = TRFLOW(:,1,1) - DTM(:)
        TAIJN(IDPOLE,JDPOLE,TIJ_RVR,:)=TAIJN(IDPOLE,JDPOLE,TIJ_RVR,:) +
     *       DTM(:)*BYDXYP(JDPOLE)
#endif
        AIJ(IDPOLE,JDPOLE,IJ_MRVR)=AIJ(IDPOLE,JDPOLE,IJ_MRVR) + DMM
        AIJ(IDPOLE,JDPOLE,IJ_ERVR)=AIJ(IDPOLE,JDPOLE,IJ_ERVR) + DGM+DPE
         FLOW(IDPOLE,JDPOLE) =  FLOW(IDPOLE,JDPOLE) + IM*DMM
        EFLOW(IDPOLE,JDPOLE) = EFLOW(IDPOLE,JDPOLE) + IM*(DGM+DPE)
#ifdef TRACERS_WATER
        TRFLOW(:,IDPOLE,JDPOLE) = TRFLOW(:,IDPOLE,JDPOLE) + IM*DTM(:)
#endif
      END IF
C****
C**** Apply net river flow to continental reservoirs
C****
      DO J=2,JM-1
        DO I=1,IM
          IF(FLAND(I,J).gt.0.) THEN
            MWL(I,J) = MWL(I,J) +  FLOW(I,J)
            GML(I,J) = GML(I,J) + EFLOW(I,J)
#ifdef TRACERS_WATER
            TRLAKE(:,1,I,J) = TRLAKE(:,1,I,J) + TRFLOW(:,I,J)
#endif
C**** remove pathologically small values
            IF (MWL(I,J).lt.1d-20) THEN
              MWL(I,J)=0.
              GML(I,J)=0.
#ifdef TRACERS_WATER
              TRLAKE(:,1:2,I,J) = 0.
#endif
            END IF
            IF (FLAKE(I,J).gt.0) THEN
              HLK1=(MLDLK(I,J)*RHOW)*TLAKE(I,J)*SHW
              MLDLK(I,J)=MLDLK(I,J)+FLOW(I,J)/(RHOW*FLAKE(I,J)*DXYP(J))
              TLAKE(I,J)=(HLK1*FLAKE(I,J)*DXYP(J)+EFLOW(I,J))
     *             /(MLDLK(I,J)*RHOW*FLAKE(I,J)*DXYP(J)*SHW)
C**** accumulate some diagnostics
              AJ(J,J_RVRD,ITLAKE) =AJ(J,J_RVRD,ITLAKE) + FLOW(I,J)*
     *             BYDXYP(J)*(1.-RSI(I,J))
              AJ(J,J_ERVR,ITLAKE) =AJ(J,J_ERVR,ITLAKE) +EFLOW(I,J)*
     *             BYDXYP(J)*(1.-RSI(I,J))
              AJ(J,J_RVRD,ITLKICE)=AJ(J,J_RVRD,ITLKICE)+ FLOW(I,J)*
     *             BYDXYP(J)*RSI(I,J)
              AJ(J,J_ERVR,ITLKICE)=AJ(J,J_ERVR,ITLKICE)+EFLOW(I,J)*
     *             BYDXYP(J)*RSI(I,J)
            ELSE
              TLAKE(I,J)=GML(I,J)/(SHW*MWL(I,J)+teeny)
C**** accounting fix to ensure river flow with no lakes is counted
              AJ(J,J_RVRD,ITLAKE)=AJ(J,J_RVRD,ITLAKE)+ FLOW(I,J)
     *             *BYDXYP(J)
              AJ(J,J_ERVR,ITLAKE)=AJ(J,J_ERVR,ITLAKE)+EFLOW(I,J)
     *             *BYDXYP(J)
            END IF
          END IF
        END DO
      END DO
      MWL(1,1) = MWL(1,1) +  FLOW(1,1)
      GML(1,1) = GML(1,1) + EFLOW(1,1)
#ifdef TRACERS_WATER
      TRLAKE(:,1,1,1) = TRLAKE(:,1,1,1) + TRFLOW(:,1,1)
#endif
      IF (MWL(1,1).lt.1d-20) THEN
        MWL(1,1)=0.
        GML(1,1)=0.
#ifdef TRACERS_WATER
        TRLAKE(:,1:2,1,1) = 0.
#endif
      END IF
      IF (FLAKE(1,1).gt.0) THEN
        HLK1=(MLDLK(1,1)*RHOW)*TLAKE(1,1)*SHW
        MLDLK(1,1)=MLDLK(1,1)+FLOW(1,1)/(RHOW*FLAKE(1,1)*DXYP(1))
        TLAKE(1,1)=(HLK1*FLAKE(1,1)*DXYP(1)+EFLOW(1,1))
     *       /(MLDLK(1,1)*RHOW*FLAKE(1,1)*DXYP(1)*SHW)
C**** accumulate some diagnostics
        AJ(1,J_RVRD,ITLAKE) =AJ(1,J_RVRD,ITLAKE) + FLOW(1,1)*
     *             BYDXYP(1)*(1.-RSI(1,1))
        AJ(1,J_ERVR,ITLAKE) =AJ(1,J_ERVR,ITLAKE) +EFLOW(1,1)*
     *       BYDXYP(1)*(1.-RSI(1,1))
        AJ(1,J_RVRD,ITLKICE)=AJ(1,J_RVRD,ITLKICE)+ FLOW(1,1)*
     *       BYDXYP(1)*RSI(1,1)
        AJ(1,J_ERVR,ITLKICE)=AJ(1,J_ERVR,ITLKICE)+EFLOW(1,1)*
     *       BYDXYP(1)*RSI(1,1)
      ELSE
        TLAKE(1,1)=GML(1,1)/(SHW*MWL(1,1)+teeny)
        AJ(1,J_RVRD,ITLAKE)=AJ(1,J_RVRD,ITLAKE)+ FLOW(1,1)*BYDXYP(1)
        AJ(1,J_ERVR,ITLAKE)=AJ(1,J_ERVR,ITLAKE)+EFLOW(1,1)*BYDXYP(1)
      END IF

      CALL PRINTLK("RV")
C**** Set GTEMP array for lakes
      DO J=1,JM
        DO I=1,IM
          IF (FLAKE(I,J).gt.0) THEN
            GTEMP(1,1,I,J)=TLAKE(I,J)
#ifdef TRACERS_WATER
            GTRACER(:,1,I,J)=TRLAKE(:,1,I,J)/(MLDLK(I,J)*RHOW*FLAKE(I,J)
     *           *DXYP(J))
#endif
            MLHC(I,J) = SHW*MLDLK(I,J)*RHOW
          END IF
        END DO
      END DO

      RETURN
C****
      END SUBROUTINE RIVERF

      SUBROUTINE diag_RIVER
!@sum  diag_RIVER prints out the river outflow for various rivers
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : rhow,sday,teeny,undef
      USE MODEL_COM, only : jyear0,amon0,jdate0,jhour0,jyear,amon
     *     ,jdate,jhour,itime,dtsrc,idacc,itime0,nday,jdpery,jmpery
      USE GEOM, only : bydxyp
      USE DAGCOM, only : aij,ij_mrvr
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm,trname,trw0,n_water,itime_tr0
     *     ,tr_wd_type,nwater
      USE TRACER_DIAG_COM, only : taijn,tij_rvr,to_per_mil,units_tij
     *     ,scale_tij
#endif
      USE LAKES, only : irvrmth,jrvrmth,namervr,nrvr
      IMPLICIT NONE
      REAL*8 RVROUT(6), SCALERVR, DAYS
      INTEGER INM,I,N
#ifdef TRACERS_WATER
      REAL*8 TRRVOUT(6,NTM)
#endif

      DAYS=(Itime-Itime0)/DFLOAT(nday)
      WRITE(6,900) JYEAR0,AMON0,JDATE0,JHOUR0,JYEAR,AMON,JDATE,JHOUR
     *     ,ITIME,DAYS
C**** convert kg/(source time step) to km^3/mon
      SCALERVR = 1d-9*SDAY*JDPERY/(JMPERY*RHOW*DTSRC)
      DO INM=1,NRVR,6
        DO I=1,6
          RVROUT(I) = SCALERVR*AIJ(IRVRMTH(I-1+INM),JRVRMTH(I-1+INM)
     *         ,IJ_MRVR)/IDACC(1)
        END DO
        WRITE(6,901) (NAMERVR(I-1+INM),RVROUT(I),I=1,6)
      END DO

#ifdef TRACERS_WATER
      DO N=1,NTM
        if (itime.ge.itime_tr0(n) .and. tr_wd_TYPE(n).eq.nWater) then
          WRITE(6,*) "River outflow tracer concentration "
     *         ,trim(units_tij(tij_rvr,n)),":",TRNAME(N)
          DO INM=1,NRVR,6
            DO I=1,6
              IF (AIJ(IRVRMTH(I-1+INM),JRVRMTH(I-1+INM),IJ_MRVR).gt.0)
     *             THEN
              if (to_per_mil(n).gt.0) then
                TRRVOUT(I,N)=1d3*(TAIJN(IRVRMTH(I-1+INM),JRVRMTH(I-1
     *               +INM),TIJ_RVR,N)/(trw0(n)*AIJ(IRVRMTH(I-1+INM)
     *               ,JRVRMTH(I-1+INM),IJ_MRVR)*BYDXYP(JRVRMTH(I-1+INM
     *               ))) -1.)
              else
                TRRVOUT(I,N)=scale_tij(TIJ_RVR,n)*TAIJN(IRVRMTH(I-1+INM)
     *               ,JRVRMTH(I-1+INM),TIJ_RVR,N)/(AIJ(IRVRMTH(I-1+INM)
     *               ,JRVRMTH(I-1+INM),IJ_MRVR)*BYDXYP(JRVRMTH(I-1+INM))
     *               +teeny)
              end if
              ELSE
                TRRVOUT(I,N)=undef
              END IF
            END DO
            WRITE(6,901) (NAMERVR(I-1+INM),TRRVOUT(I,N),I=1,6)
          END DO
        end if
      END DO
#endif

      RETURN
C****
 900  FORMAT ('1* River Outflow (km^3/mon) **  From:',I6,A6,I2,',  Hr'
     *     ,I3,6X,'To:',I6,A6,I2,', Hr',I3,'  Model-Time:',I9,5X
     *     ,'Dif:',F7.2,' Days')
 901  FORMAT (' ',A8,':',F8.3,5X,A8,':',F8.3,5X,A8,':',F8.3,5X,
     *            A8,':',F8.3,5X,A8,':',F8.3,5X,A8,':',F8.3)
      END SUBROUTINE diag_RIVER

      SUBROUTINE CHECKL (SUBR)
!@sum  CHECKL checks whether the lake variables are reasonable.
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0 (based on LB265)
      USE CONSTANT, only : rhow
      USE MODEL_COM, only : im,jm,hlake,fearth,qcheck
      USE GEOM, only : dxyp,imaxj
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm, trname, t_qlimit
#endif
      USE LAKES
      USE LAKES_COM
      IMPLICIT NONE

      INTEGER I,J,N !@var I,J loop variables
      CHARACTER*6, INTENT(IN) :: SUBR
      LOGICAL QCHECKL
#ifdef TRACERS_WATER
      integer :: imax,jmax
      real*8 relerr,errmax
#endif

C**** Check for NaN/INF in lake data
      CALL CHECK3(MWL,IM,JM,1,SUBR,'mwl')
      CALL CHECK3(GML,IM,JM,1,SUBR,'gml')
      CALL CHECK3(MLDLK,IM,JM,1,SUBR,'mld')
      CALL CHECK3(TLAKE,IM,JM,1,SUBR,'tlk')

      QCHECKL = .FALSE.
      DO J=2,JM-1
      DO I=1,IM
        IF(FEARTH(I,J).gt.0.) THEN
C**** check for negative mass
          IF (MWL(I,J).lt.0 .or. MLDLK(I,J).lt.0) THEN
            WRITE(6,*) 'After ',SUBR,': I,J,TSL,MWL,GML,MLD=',I,J
     *           ,TLAKE(I,J),MWL(I,J),GML(I,J),MLDLK(I,J)
            QCHECKL = .TRUE.
          END IF
C**** check for reasonable lake surface temps
          IF (TLAKE(I,J).ge.50 .or. TLAKE(I,J).lt.-0.5) THEN
            WRITE(6,*) 'After ',SUBR,': I,J,TSL=',I,J,TLAKE(I,J)
c            QCHECKL = .TRUE.
          END IF
        END IF
C**** Check total lake mass ( <0.4 m, >20x orig depth)
        IF(FLAKE(I,J).gt.0.) THEN
          IF(MWL(I,J).lt.0.4d0*RHOW*DXYP(J)*FLAKE(I,J)) THEN
            WRITE (6,*) 'After ',SUBR,
     *           ': I,J,FLAKE,HLAKE,lake level low=',I,J,FLAKE(I,J),
     *           HLAKE(I,J),MWL(I,J)/(RHOW*DXYP(J)*FLAKE(I,J))
          END IF
          IF(MWL(I,J).gt.RHOW*MAX(20.*HLAKE(I,J),3d1)*DXYP(J)*FLAKE(I,J)
     *         )THEN
            WRITE (6,*) 'After ',SUBR,
     *           ': I,J,FLAKE,HLAKE,lake level high=',I,J,FLAKE(I,J),
     *           HLAKE(I,J),MWL(I,J)/(RHOW*DXYP(J)*FLAKE(I,J))
          END IF
        END IF
      END DO
      END DO

#ifdef TRACERS_WATER
      do n=1,ntm
C**** Check for neg tracers in lake
        if (t_qlimit(n)) then
          do j=1,jm
          do i=1,imaxj(j)
            if (fearth(i,j).gt.0) then
              if (trlake(n,1,i,j).lt.0 .or. trlake(n,2,i,j).lt.0) then
                print*,"Neg tracer in lake after ",SUBR,i,j,trname(n)
     *               ,trlake(n,:,i,j) 
                QCHECKL=.TRUE.
              end if
            end if
          end do
          end do
        end if
C**** Check conservation of water tracers in lake
        if (trname(n).eq.'Water') then
          errmax = 0. ; imax=1 ; jmax=1
          do j=1,jm
          do i=1,imaxj(j)
            if (fearth(i,j).gt.0) then
              if (flake(i,j).gt.0) then
                relerr=max(
     *               abs(trlake(n,1,i,j)-mldlk(i,j)*rhow*flake(i,j)
     *               *dxyp(j))/trlake(n,1,i,j),abs(trlake(n,1,i,j)
     *               +trlake(n,2,i,j)-mwl(i,j))/(trlake(n,1,i,j)
     *               +trlake(n,2,i,j)))
              else
                if ((mwl(i,j).eq.0 .and. trlake(n,1,i,j)+trlake(n,2,i,j)
     *               .gt.0) .or. (mwl(i,j).gt.0 .and. trlake(n,1,i,j)
     *               +trlake(n,2,i,j).eq.0))  then
                  print*,"CHECKL ",SUBR,i,j,mwl(i,j),trlake(n,1:2,i,j)
                  relerr=0.
                else
                  if (mwl(i,j).gt.1d-20) then
                    relerr=abs(trlake(n,1,i,j)
     *                 +trlake(n,2,i,j)-mwl(i,j))/(trlake(n,1,i,j)
     *                 +trlake(n,2,i,j))
                  else
                    if (mwl(i,j).gt.0) print*,"CHECKL2 ",SUBR,i,j,mwl(i
     *                   ,j),trlake(n,1:2,i,j)
                    relerr=0.
                  end if
                end if
              end if
              if (relerr.gt.errmax) then
                imax=i ; jmax=j ; errmax=relerr
              end if
            end if
          end do
          end do
          print*,"Relative error in lake mass after ",trim(subr),":"
     *         ,imax,jmax,errmax,trlake(n,:,imax,jmax),mldlk(imax,jmax)
     *         *rhow*flake(imax,jmax)*dxyp(jmax),mwl(imax,jmax)
     *         -mldlk(imax,jmax)*rhow*flake(imax,jmax)*dxyp(jmax)
        end if
      end do
#endif

      IF (QCHECKL)
     &     call stop_model('CHECKL: Lake variables out of bounds',255)
      RETURN
C****
      END SUBROUTINE CHECKL

      SUBROUTINE daily_LAKE
!@sum  daily_LAKE does lake things at the beginning of every day
!@auth G. Schmidt
!@ver  1.0
      IMPLICIT NONE

C**** Experimental code: not yet functional
C**** Update lake fraction as a function of lake mass at end of day
C**** Assume lake is conical
C****   => A = pi*(h*tanlk)^2, M=(1/3)*pi*rho*h*(h*tanlk)^2
C****
c        PRINT*,"TEST FLAKE CHANGE"
c        DO J=1,JM
c        DO I=1,IMAXJ(J)
c          IF (FLAKE(I,J).gt.0) THEN
c            PRINT*,"tan(a)",I,J,TANLK(I,J)
c            PRINT*,FLAKE(I,J),(9d0*PI*(TANLK(I,J)*MWL(I,J)/RHOW)**2)
c     *           **BY3/DXYP(J)
C**** remove lakes that are too small ??
c            IF (FLAKE(I,J).lt.0.005) THEN
c          PRINT*,"Lake->0 (too small)",I,J,MWL(I,J),GML(I,J),MLDLK(I,J)
c          MLDLK(I,J)=MIN(MLDLK(I,J),MWL(I,J)/(RHOW*FLAKE(I,J)*DXYP(J)))
c              FLAKE(I,J) = 0.
c          END IF
c        END DO
c        END DO
C****
      CALL PRINTLK("DY")
C****
      RETURN
      END SUBROUTINE daily_LAKE

      SUBROUTINE PRECIP_LK
!@sum  PRECIP_LK driver for applying precipitation/melt to lake fraction
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : rhow,shw,teeny
      USE MODEL_COM, only : im,jm,fland,flice,itlake,itlkice
      USE GEOM, only : imaxj,dxyp,bydxyp
      USE SEAICE_COM, only : rsi
      USE LAKES_COM, only : mwl,gml,tlake,mldlk,flake
#ifdef TRACERS_WATER
     *     ,trlake,ntm
#endif
      USE FLUXES, only : runpsi,runoli,prec,eprec,gtemp,melti,emelti
#ifdef TRACERS_WATER
     *     ,trunpsi,trunoli,trprec,gtracer,trmelti
#endif
      USE DAGCOM, only : aj,j_run
      IMPLICIT NONE

      REAL*8 PRCP,ENRGP,PLICE,PLKICE,RUN0,ERUN0,POLAKE,HLK1
      INTEGER I,J,ITYPE
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(NTM) :: TRUN0
#endif

      CALL PRINTLK("PR")

      DO J=1,JM
      DO I=1,IMAXJ(J)
      IF (FLAKE(I,J)+FLICE(I,J).gt.0) THEN
        POLAKE=(1.-RSI(I,J))*FLAKE(I,J)
        PLKICE=RSI(I,J)*FLAKE(I,J)
        PLICE=FLICE(I,J)
        PRCP=PREC(I,J)
        ENRGP=EPREC(I,J)        ! energy of precipitation

C**** calculate fluxes over whole box
        RUN0 =POLAKE*PRCP  + PLKICE* RUNPSI(I,J) + PLICE* RUNOLI(I,J)
        ERUN0=POLAKE*ENRGP ! PLKICE*ERUNPSI(I,J) + PLICE*ERUNOLI(I,J) =0

C**** simelt is given as kg, so divide by area
        IF (FLAKE(I,J).gt.0) THEN
          RUN0  =RUN0+ MELTI(I,J)*BYDXYP(J)
          ERUN0=ERUN0+EMELTI(I,J)*BYDXYP(J)
        END IF

        MWL(I,J) = MWL(I,J) +  RUN0*DXYP(J)
        GML(I,J) = GML(I,J) + ERUN0*DXYP(J)
#ifdef TRACERS_WATER
        TRUN0(:) = POLAKE*TRPREC(:,I,J)*BYDXYP(J)
     *       + PLKICE*TRUNPSI(:,I,J) + PLICE *TRUNOLI(:,I,J)
        IF (FLAKE(I,J).gt.0) TRUN0(:)=TRUN0(:)+TRMELTI(:,I,J)*BYDXYP(J)
        TRLAKE(:,1,I,J)=TRLAKE(:,1,I,J) + TRUN0(:)*DXYP(J)
#endif

        IF (FLAKE(I,J).gt.0) THEN
          HLK1=TLAKE(I,J)*MLDLK(I,J)*RHOW*SHW
          MLDLK(I,J)=MLDLK(I,J) + RUN0/(FLAKE(I,J)*RHOW)
          TLAKE(I,J)=(HLK1*FLAKE(I,J)+ERUN0)/(MLDLK(I,J)*FLAKE(I,J)
     *         *RHOW*SHW)
          GTEMP(1,1,I,J)=TLAKE(I,J)
#ifdef TRACERS_WATER
          GTRACER(:,1,I,J)=TRLAKE(:,1,I,J)/(MLDLK(I,J)*RHOW*FLAKE(I,J)
     *         *DXYP(J))
#endif
          AJ(J,J_RUN,ITLAKE) =AJ(J,J_RUN,ITLAKE) -PLICE*RUNOLI(I,J)
     *         *(1.-RSI(I,J))
          AJ(J,J_RUN,ITLKICE)=AJ(J,J_RUN,ITLKICE)-PLICE*RUNOLI(I,J)
     *         *RSI(I,J)
        ELSE
          TLAKE(I,J)=GML(I,J)/(MWL(I,J)*SHW+teeny)
C**** accounting fix to ensure runoff with no lakes is counted
C**** no regional diagnostics required
          AJ(J,J_RUN,ITLAKE) =AJ(J,J_RUN,ITLAKE)-PLICE*RUNOLI(I,J)
        END IF

      END IF
      END DO
      END DO
      RETURN
C****
      END SUBROUTINE PRECIP_LK

      SUBROUTINE GROUND_LK
!@sum  GROUND_LK driver for applying surface fluxes to lake fraction
!@auth Gavin Schmidt
!@ver  1.0
!@calls
      USE CONSTANT, only : rhow,shw,teeny
      USE MODEL_COM, only : im,jm,flice,fland,hlake
     *     ,fearth,dtsrc,itlake,itlkice
      USE GEOM, only : imaxj,dxyp
      USE FLUXES, only : runosi, erunosi, e0, evapor, dmsi, dhsi, dssi,
     *     runoli, runoe, erunoe, solar, dmua, dmva, gtemp
#ifdef TRACERS_WATER
     *     ,trunoli,trunoe,trevapor,dtrsi,trunosi,gtracer
#endif
      USE SEAICE_COM, only : rsi
      USE DAGCOM, only : aj,areg,jreg,j_wtr1,j_wtr2,j_run,j_erun
      USE LAKES_COM, only : mwl,gml,tlake,mldlk,flake
#ifdef TRACERS_WATER
     *     ,trlake,ntm
      USE TRACER_DIAG_COM,only: taijn,tij_lk1,tij_lk2
#endif
      USE LAKES, only : lkmix,lksourc,byzeta,minmld
      IMPLICIT NONE
C**** grid box variables
      REAL*8 ROICE, POLAKE, PLKICE, PEARTH, PLICE
!@var MLAKE,ELAKE mass and energy /m^2 for lake model layers
      REAL*8, DIMENSION(2) :: MLAKE,ELAKE
C**** fluxes
      REAL*8 EVAPO, FIDT, FODT, RUN0, ERUN0, RUNLI, RUNE, ERUNE,
     *     HLK1,TLK1,TLK2,TKE,SROX(2),FSR2, U2RHO
C**** output from LKSOURC
      REAL*8 ENRGFO, ACEFO, ACEFI, ENRGFI
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(NTM) :: TRUN0,TRO,TRI,TREVAP,TOTTRL
      REAL*8, DIMENSION(NTM,2) :: TRLAKEL
#endif
      INTEGER I,J,JR

      CALL PRINTLK("GR")

      DO J=1,JM
      DO I=1,IMAXJ(J)
      JR=JREG(I,J)
      ROICE=RSI(I,J)
      PLKICE=FLAKE(I,J)*ROICE
      POLAKE=FLAKE(I,J)*(1.-ROICE)
C**** Add land ice and surface runoff to lake variables
      IF (FLAND(I,J).gt.0) THEN
        PLICE =FLICE(I,J)
        PEARTH=FEARTH(I,J)
        RUNLI=RUNOLI(I,J)
        RUNE =RUNOE(I,J)
        ERUNE=ERUNOE(I,J)
C**** calculate flux over whole box
        RUN0 =RUNLI*PLICE + RUNE*PEARTH
        ERUN0=             ERUNE*PEARTH
        MWL(I,J) = MWL(I,J) + RUN0*DXYP(J)
        GML(I,J) = GML(I,J) +ERUN0*DXYP(J)
#ifdef TRACERS_WATER
        TRLAKE(:,1,I,J)=TRLAKE(:,1,I,J)+
     *       (TRUNOLI(:,I,J)*PLICE+TRUNOE(:,I,J)*PEARTH)*DXYP(J)
#endif
        IF (FLAKE(I,J).gt.0) THEN
          HLK1=TLAKE(I,J)*MLDLK(I,J)*RHOW*SHW
          MLDLK(I,J)=MLDLK(I,J) + RUN0/(FLAKE(I,J)*RHOW)
          TLAKE(I,J)=(HLK1*FLAKE(I,J)+ERUN0)/(MLDLK(I,J)*FLAKE(I,J)
     *         *RHOW*SHW)
#ifdef TRACERS_WATER
          GTRACER(:,1,I,J)=TRLAKE(:,1,I,J)/(MLDLK(I,J)*RHOW*FLAKE(I,J)
     *         *DXYP(J))
#endif
          AJ(J,J_RUN,ITLAKE)=AJ(J,J_RUN,ITLAKE)-
     *         (RUNE*PEARTH+RUNLI*PLICE)*(1.-RSI(I,J))
          AJ(J,J_RUN,ITLKICE)=AJ(J,J_RUN,ITLKICE)-
     *         (RUNE*PEARTH+RUNLI*PLICE)*RSI(I,J)
          AJ(J,J_ERUN,ITLAKE )=AJ(J,J_ERUN,ITLAKE )-ERUNE*PEARTH*
     *         (1.-RSI(I,J))
          AJ(J,J_ERUN,ITLKICE)=AJ(J,J_ERUN,ITLKICE)-ERUNE*PEARTH*
     *         RSI(I,J)
        ELSE
          TLAKE(I,J)=GML(I,J)/(MWL(I,J)*SHW+teeny)
C**** accounting fix to ensure runoff with no lakes is counted
C**** no regional diagnostics required
          AJ(J,J_RUN,ITLAKE)=AJ(J,J_RUN,ITLAKE)-
     *         (RUNE*PEARTH+RUNLI*PLICE)
          AJ(J,J_ERUN,ITLAKE )=AJ(J,J_ERUN,ITLAKE )-ERUNE*PEARTH
        END IF
      END IF

      IF (FLAKE(I,J).gt.0) THEN
        TLK1 =TLAKE(I,J)
        EVAPO=EVAPOR(I,J,1)     ! evap/dew over open lake (kg/m^2)
        FODT =E0(I,J,1)         ! net heat over open lake (J/m^2)
        SROX(1)=SOLAR(1,I,J)      ! solar radiation open lake (J/m^2)
        SROX(2)=SOLAR(3,I,J)      ! solar radiation through ice (J/m^2)
        FSR2 =EXP(-MLDLK(I,J)*BYZETA)
C**** get ice-ocean fluxes from sea ice routine (over ice fraction)
        RUN0 =RUNOSI(I,J) ! includes ACE2M + basal term
        FIDT =ERUNOSI(I,J)
C**** calculate kg/m^2, J/m^2 from saved variables
        MLAKE(1)=MLDLK(I,J)*RHOW
        MLAKE(2)=MAX(MWL(I,J)/(FLAKE(I,J)*DXYP(J))-MLAKE(1),0d0)
        ELAKE(1)=TLK1*SHW*MLAKE(1)
        ELAKE(2)=GML(I,J)/(FLAKE(I,J)*DXYP(J))-ELAKE(1)
#ifdef TRACERS_WATER
        TRLAKEL(:,:)=TRLAKE(:,:,I,J)/(FLAKE(I,J)*DXYP(J))
        TRUN0(:)=TRUNOSI(:,I,J)
        TREVAP(:)=TREVAPOR(:,1,I,J)
#endif
        IF (MLAKE(2).lt.1d-10) THEN
          MLAKE(2)=0.
          ELAKE(2)=0.
#ifdef TRACERS_WATER
          TRLAKEL(:,2)=0.
#endif
        END IF

C**** Limit FSR2 in the case of thin second layer
        FSR2=MIN(FSR2,MLAKE(2)/(MLAKE(1)+MLAKE(2)))

C**** Apply fluxes and calculate the amount of frazil ice formation
        CALL LKSOURC (I,J,ROICE,MLAKE,ELAKE,RUN0,FODT,FIDT,SROX,FSR2,
#ifdef TRACERS_WATER
     *       TRLAKEL,TRUN0,TREVAP,TRO,TRI,
#endif
     *       EVAPO,ENRGFO,ACEFO,ACEFI,ENRGFI)

C**** Mixing and entrainment
C**** Calculate turbulent kinetic energy for lake
c       U2rho=(1.-ROICE)*SQRT(DMUA(I,J,1)**2+DMVA(I,J,1)**2)/DTSRC
c       TKE=0.5 * (19.3)^(2/3) * U2rho /rhoair ! (m/s)^2
        TKE=0.  ! 3.6d0*U2rho/rhoair*MLAKE(1)  ! (J/m^2)

        CALL LKMIX (MLAKE,ELAKE,
#ifdef TRACERS_WATER
     *       TRLAKEL,
#endif
     *       HLAKE(I,J),TKE,ROICE,DTSRC)

C**** Resave prognostic variables
        MWL(I,J)  =(MLAKE(1)+MLAKE(2))*(FLAKE(I,J)*DXYP(J))
        GML(I,J)  =(ELAKE(1)+ELAKE(2))*(FLAKE(I,J)*DXYP(J))
        MLDLK(I,J)= MLAKE(1)/RHOW
        IF (MLAKE(2).eq.0.) MLDLK(I,J)=MIN(MINMLD,MLDLK(I,J))
        TLAKE(I,J)= ELAKE(1)/(SHW*MLAKE(1))
        IF (MLAKE(2).gt.0) THEN
          TLK2    = ELAKE(2)/(SHW*MLAKE(2))
        ELSE
          TLK2    = TLAKE(I,J)
        END IF
#ifdef TRACERS_WATER
        IF (MLAKE(2).eq.0. .and. MLAKE(1)-MLDLK(I,J)*RHOW.gt.1d-10) THEN
          TOTTRL(:)=TRLAKEL(:,1)
          TRLAKEL(:,2)=(MLAKE(1)-MLDLK(I,J)*RHOW)*TRLAKEL(:,1)/MLAKE(1)
          TRLAKEL(:,1)=TOTTRL(:)-TRLAKEL(:,2)
        END IF
        TRLAKE(:,:,I,J)=TRLAKEL(:,:)*(FLAKE(I,J)*DXYP(J))
        GTRACER(:,1,I,J)=TRLAKEL(:,1)/(MLDLK(I,J)*RHOW)
#endif
        GTEMP(1,1,I,J)=TLAKE(I,J)
        GTEMP(2,1,I,J)=TLK2       ! diagnostic only

C**** Open lake diagnostics
        AJ(J,J_WTR1, ITLAKE)=AJ(J,J_WTR1, ITLAKE)+MLAKE(1)*POLAKE
        AJ(J,J_WTR2, ITLAKE)=AJ(J,J_WTR2, ITLAKE)+MLAKE(2)*POLAKE
C**** Ice-covered ocean diagnostics
        AJ(J,J_WTR1, ITLKICE)=AJ(J,J_WTR1, ITLKICE)+MLAKE(1)*PLKICE
        AJ(J,J_WTR2, ITLKICE)=AJ(J,J_WTR2, ITLKICE)+MLAKE(2)*PLKICE
C**** regional diags
        AREG(JR,J_WTR1)=AREG(JR,J_WTR1)+MLAKE(1)*FLAKE(I,J)*DXYP(J)
        AREG(JR,J_WTR2)=AREG(JR,J_WTR2)+MLAKE(2)*FLAKE(I,J)*DXYP(J)
#ifdef TRACERS_WATER
C**** tracer diagnostics
        TAIJN(I,J,tij_lk1,:)=TAIJN(I,J,tij_lk1,:)+TRLAKEL(:,1) !*PLKICE?
        TAIJN(I,J,tij_lk2,:)=TAIJN(I,J,tij_lk2,:)+TRLAKEL(:,2) !*PLKICE?
#endif

C**** Store mass and energy fluxes for formation of sea ice
        DMSI(1,I,J)=ACEFO
        DMSI(2,I,J)=ACEFI
        DHSI(1,I,J)=ENRGFO
        DHSI(2,I,J)=ENRGFI
        DSSI(:,I,J)=0.     ! always zero salinity
#ifdef TRACERS_WATER
        DTRSI(:,1,I,J)=TRO(:)
        DTRSI(:,2,I,J)=TRI(:)
#endif
      END IF
      END DO
      END DO

      CALL PRINTLK("G2")

      RETURN
C****
      END SUBROUTINE GROUND_LK


      SUBROUTINE PRINTLK(STR)
!@sum  PRINTLK print out selected diagnostics from specified lakes
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : lhm,byshi,rhow,shw
      USE MODEL_COM, only : qcheck
      USE GEOM, only : dxyp
      USE LAKES_COM, only : tlake,mwl,mldlk,gml,flake
      USE SEAICE, only : xsi,ace1i,rhoi
      USE SEAICE_COM, only : rsi,hsi,msi,snowi
      IMPLICIT NONE
      CHARACTER*2, INTENT(IN) :: STR
      INTEGER, PARAMETER :: NDIAG=1
      INTEGER I,J,N
      INTEGER, DIMENSION(NDIAG) :: IDIAG = (/24/),
     *                             JDIAG = (/15/)
      REAL*8 HLK2,TLK2, TSIL(4)

      IF (.NOT.QCHECK) RETURN

      DO N=1,NDIAG
        I=IDIAG(N)
        J=JDIAG(N)
        HLK2 = MWL(I,J)/(RHOW*FLAKE(I,J)*DXYP(J)) - MLDLK(I,J)
        IF (HLK2.gt.0) THEN
        TLK2 = (GML(I,J)/(SHW*RHOW*FLAKE(I,J)*DXYP(J)) -
     *       TLAKE(I,J)*MLDLK(I,J))/HLK2
        ELSE
          TLK2=0.
        END IF
        TSIL(:)=0.
        IF (RSI(I,J).gt.0) THEN
          TSIL(1:2) = (HSI(1:2,I,J)/(XSI(1:2)*(ACE1I+SNOWI(I,J)))+LHM)
     *         *BYSHI
          TSIL(3:4) = (HSI(3:4,I,J)/(XSI(3:4)*MSI(I,J))+LHM)*BYSHI
        END IF
        WRITE(99,*) STR,I,J,FLAKE(I,J),TLAKE(I,J),TLK2,MLDLK(I,J),HLK2
     *       ,RSI(I,J),MSI(I,J)/RHOI,SNOWI(I,J)/RHOW,TSIL(1:4)
      END DO

      RETURN
      END  SUBROUTINE PRINTLK

      SUBROUTINE conserv_LKM(LKM)
!@sum  conserv_LKM calculates zonal lake mass
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,jm,fland,fearth,fim
      USE GEOM, only : imaxj,bydxyp
      USE LAKES_COM, only : mwl
      IMPLICIT NONE
      REAL*8, DIMENSION(JM) :: LKM
      INTEGER :: I,J
C****
C**** LAKE MASS (kg/m^2)
C****
      DO J=1,JM
        LKM(J)=0.
        DO I=1,IMAXJ(J)
          IF (FLAND(I,J).gt.0) LKM(J)=LKM(J)+MWL(I,J)
        END DO
        LKM(J)=LKM(J)*BYDXYP(J)
      END DO
      LKM(1) =FIM*LKM(1)
      LKM(JM)=FIM*LKM(JM)
      RETURN
      END SUBROUTINE conserv_LKM

      SUBROUTINE conserv_LKE(LKE)
!@sum  conserv_LKE calculates zonal lake energy
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,jm,zatmo,fim,fland,fearth
      USE GEOM, only : imaxj,bydxyp
      USE LAKES_COM, only : gml,mwl
      IMPLICIT NONE
      REAL*8, DIMENSION(JM) :: LKE
      INTEGER :: I,J
      REAL*8 AREA
C****
C**** LAKE ENERGY (J/m^2) (includes potential energy (DISABLED))
C****
      DO J=1,JM
        LKE(J)=0.
        DO I=1,IMAXJ(J)
          IF (FLAND(I,J).gt.0) LKE(J)=LKE(J)+GML(I,J)
c     *         +ZATMO(I,J)*MWL(I,J)
        END DO
        LKE(J)=LKE(J)*BYDXYP(J)
      END DO
      LKE(1) =FIM*LKE(1)
      LKE(JM)=FIM*LKE(JM)
      RETURN
      END SUBROUTINE conserv_LKE
