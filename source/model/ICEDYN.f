C*************************************************************
C PLEASE KEEP THIS NOTE OF MODEL-DEVELOPMENT HISTORY
C Matrix solve uses Thomas algorithm, 10/1991, Jinlun Zhang
C Spherical coordinate system, 10/27/93, Jinlun Zhang
C Latest finite differencing scheme for treatment of NP,
C 9/9/1996,Jinlun Zhang
C Alternating direction implicit (ADI) method is used, 10/1998,
C Jinlun Zhang
C For details about ADI dynamics model, see Zhang and Rothrock,
C "Modeling
C Arctic Sea ice with an efficient plastic solution",
C submitted to JGR, 1999
C Adapted for GISS coupled model Jiping Liu/Gavin Schmidt 2000
C - Further modularised May 2001
C*************************************************************

      MODULE ICEDYN
!@sum  ICEDYN holds local variables for dynamic sea ice
!@auth Gavin Schmidt (based on code from Jinlun Zhang)
      USE CONSTANT, only : radian,radius
      USE MODEL_COM, only : im,jm
      USE SEAICE, only : osurf_tilt
      IMPLICIT NONE
      SAVE
C**** Definition for ice advection grid (EDIT FOR ADVSI GRID CHANGE)
      INTEGER, PARAMETER :: IMIC=IM, JMIC=JM

C**** local grid variables for ice rheology scheme
C**** Edit the definition of nx1,ny1 to change the grid for the
C**** rheology calculations without changing ADVSI grid.

!@var nx1 number of grid points in the longitudinal direction
!@+  (calculated points from 2 through nx1-1. End points are boundaries)
!@var ny1 number of grid points in the latitudinal direction
!@+  (calculated points from 2 through ny1-1. End points are boundaries)
      integer, parameter :: nx1=imic+2, ny1=jmic
      INTEGER, parameter :: NYPOLE=NY1-1,NXLCYC=NX1-1
      integer :: NPOL=1,LCYC=1

!@var FOCEAN land/ocean mask on ice dynamic grid
      REAL*8, DIMENSION(NX1-2,NY1) :: FOCEAN

C**** input
!@var HEFFM ice mass mask (1/0)
!@var UVM ice velocity mask (1/0)
!@var COR coriolis term for ice dynamic equation
!@var GAIRX,GAIRY atmosphere-ice stress (B grid)
!@var GWATX,GWATY ocean velocities (B grid) (m/s)
!@var PGFUB,PGFVB pressure accelaration force
!@var AMASS ice mass (kg/m^2)
!@var UICEC,VICEC velocity arrays  (m/s)
!@var UICE,VICE velocity arrays  (m/s)
!@var HEFF ice thickness (mean over box) (m)
!@var AREA ice area (frac)
!@var UIB,VIB velocity arrays (m/s)  (????)
C**** internal variables
!@var PRESS ice internal pressure (Pa)
!@var FORCEX,FORCEY external force 
!@var DRAGS,DRAGA symmetric/anti-symmetric drag terms
!@var ZMAX,ZMIN max,min values of ZETA
!@var ETA,ZETA viscosities
C**** output
!@var DWATN non-linear water drag term
!@var DMU,DMV ice-ocean stress
      REAL*8, DIMENSION(NX1,NY1) :: PRESS,HEFFM,UVM,DWATN,COR
     *     ,ZMAX,ZMIN,ETA,ZETA,DRAGS,DRAGA,GAIRX,GAIRY
     *     ,GWATX,GWATY,PGFUB,PGFVB,FORCEX,FORCEY,AMASS,UICEC,VICEC,UIB
     *     ,VIB,DMU,DMV,HEFF,AREA
      REAL*8, DIMENSION(NX1,NY1,3) :: UICE,VICE
C**** Geometry 
!@var SINEN sin(phi)
!@var BYDXDY 
!@var DXT,DXU x-direction distances on tracer and velocity grid
!@var DYT,DYU y-direction distances on tracer and velocity grid
      REAL*8, DIMENSION(NX1,NY1) :: SINEN,BYDXDY
      REAL*8, DIMENSION(NX1) :: DXT,DXU,BYDX2,BYDXR
      REAL*8, DIMENSION(NY1) :: DYT,DYU,BYDY2,BYDYR,CST,CSU,TNGT,TNG
     *     ,BYCSU

!@var OIPHI ice-ocean turning angle (25 degrees)
!@var ECCEN value of eccentricity for yield curve ellipse
      REAL*8, PARAMETER :: ECCEN=2.0, OIPHI=25d0*radian
!@var SINWAT,COSWAT sin and cos of ice-ocean turning angle
      REAL*8 SINWAT,COSWAT

!@var PSTAR maximum sea ice pressure (Pa)
      REAL*8, PARAMETER :: PSTAR=2.75d4

!@var BYDTS reciprocal of timestep in ice dynamics code
      REAL*8 :: BYDTS

      CONTAINS

      SUBROUTINE FORM
!@sum  FORM calculates ice dynamics input parameters for relaxation
!@auth Jiping Liu/Gavin Schmidt (based on code from J. Zhang)
!@ver  1.0
      IMPLICIT NONE
      INTEGER I,J
      REAL*8 AAA
C****
C**** Set up non linear water drag
C****
      DO J=1,NY1-1
      DO I=1,NX1-1
        DWATN(I,J)=5.5*SQRT((UICE(I,J,1)-GWATX(I,J))**2
     1       +(VICE(I,J,1)-GWATY(I,J))**2)
      END DO
      END DO
C NOW SET UP SYMMETRIC DRAG
      DO J=1,NY1-1
      DO I=1,NX1-1
        DRAGS(I,J)=DWATN(I,J)*COSWAT
      END DO
      END DO
C NOW SET UP ANTI SYMMETRIC DRAG PLUS CORIOLIS
      DO J=1,NY1
      DO I=1,NX1
        IF(J.GT.NY1/2) THEN
          DRAGA(I,J)=DWATN(I,J)*SINWAT+COR(I,J)
        ELSE
          DRAGA(I,J)=DWATN(I,J)*(-SINWAT)+COR(I,J)
        END IF
      END DO
      END DO
C NOW SET UP FORCING FIELD
      DO J=1,NY1
      DO I=1,NX1

C FIRST DO WIND
        FORCEX(I,J)=GAIRX(i,j)
        FORCEY(I,J)=GAIRY(i,j)

C NOW ADD IN CURRENT FORCE
        IF(J.GT.NY1/2) THEN
          FORCEX(I,J)=FORCEX(I,J)+DWATN(I,J)*(COSWAT*GWATX(I,J)
     1         -SINWAT*GWATY(I,J))
          FORCEY(I,J)=FORCEY(I,J)+DWATN(I,J)*(SINWAT*GWATX(I,J)
     1         +COSWAT*GWATY(I,J))
        ELSE
          FORCEX(I,J)=FORCEX(I,J)+DWATN(I,J)*(COSWAT*GWATX(I,J)
     1         +SINWAT*GWATY(I,J))
          FORCEY(I,J)=FORCEY(I,J)+DWATN(I,J)*(-SINWAT*GWATX(I,J)
     1         +COSWAT*GWATY(I,J))
        END IF
        
C     NOW ADD IN TILT 
        if (osurf_tilt.eq.1) then
C**** This assumes explicit knowledge of sea surface tilt
          FORCEX(I,J)=FORCEX(I,J)+AMASS(I,J)*PGFUB(I,J)
          FORCEY(I,J)=FORCEY(I,J)+AMASS(I,J)*PGFVB(I,J)
        else
C**** Otherwise estimate tilt using geostrophy
          FORCEX(I,J)=FORCEX(I,J)-COR(I,J)*GWATY(I,J)
          FORCEY(I,J)=FORCEY(I,J)+COR(I,J)*GWATX(I,J)
        end if
      END DO
      END DO

C NOW SET UP ICE PRESSURE AND VISCOSITIES
      DO J=1,NY1
      DO I=1,NX1
        PRESS(I,J)=PSTAR*HEFF(I,J)*EXP(-20.0*(1.0-AREA(I,J)))
        ZMAX(I,J)=(5d12/2d4)*PRESS(I,J)
c       ZMIN(I,J)=0.0D+00
        ZMIN(I,J)=4d8
      END DO
      END DO

      CALL PLAST

      AAA=0.0
      DO I=2,NX1-1
        AAA=AAA+PRESS(I,NY1-1)
      END DO
      AAA=AAA/FLOAT(NX1-2)
      DO I=1,NX1
        PRESS(I,NY1)=AAA
      END DO

 8481 CONTINUE

      DO J=1,NY1
        PRESS(1,J)=PRESS(NX1-1,J)
        PRESS(NX1,J)=PRESS(2,J)
      END DO

C NOW SET VISCOSITIES AND PRESSURE EQUAL TO ZERO AT OUTFLOW PTS

      DO J=1,NY1
      DO I=1,NX1
        PRESS(I,J)=PRESS(I,J)*HEFFM(I,J)
        ETA(I,J)=ETA(I,J)*HEFFM(I,J)
        ZETA(I,J)=ZETA(I,J)*HEFFM(I,J)
      END DO
      END DO

C NOW CALCULATE PRESSURE FORCE AND ADD TO EXTERNAL FORCE
      DO J=1,NY1-1
      DO I=1,NX1-1
      FORCEX(I,J)=FORCEX(I,J)-(0.25/(DXU(I)*CSU(J)))
     1  *(PRESS(I+1,J)+PRESS(I+1,J+1)-PRESS(I,J)-PRESS(I,J+1))
      FORCEY(I,J)=FORCEY(I,J)-0.25/DYU(J)
     1  *(PRESS(I,J+1)+PRESS(I+1,J+1)-PRESS(I,J)-PRESS(I+1,J))
C NOW PUT IN MINIMAL MASS FOR TIME STEPPING CALCULATIONS
      END DO
      END DO

      DO J=1,NY1
        FORCEX(1,J)=FORCEX(NX1-1,J)
        FORCEY(1,J)=FORCEY(NX1-1,J)
        FORCEX(NX1,J)=FORCEX(2,J)
        FORCEY(NX1,J)=FORCEY(2,J)
      END DO

      RETURN
      END SUBROUTINE FORM

      SUBROUTINE PLAST
!@sum  PLAST Calculates strain rates and viscosity for dynamic ice
!@auth Jiping Liu/Gavin Schmidt (based on code from J. Zhang)
!@ver  1.0
      IMPLICIT NONE
      REAL*8, DIMENSION(NX1,NY1) :: E11,E22,E12
c      REAL*8 :: SS11
      REAL*8, PARAMETER :: ECM2 = 1.0/(ECCEN**2),GMIN=1d-20
      REAL*8 DELT,DELT1,AAA
      INTEGER I,J

C EVALUATE STRAIN RATES
      DO J=2,NY1-1
        DO I=2,NX1-1
          E11(I,J)=0.5/(DXT(I)*CST(J))*(UICE(I,J,1)+UICE(I,J-1,1)
     *         -UICE(I-1,J,1)-UICE(I-1,J-1,1))-0.25*(VICE(I,J,1)+
     *         VICE(I-1,J,1)+VICE(I-1,J-1,1)+VICE(I,J-1,1))*TNGT(J)
     *         /RADIUS
          E22(I,J)=0.5/DYT(J)*(VICE(I,J,1)+VICE(I-1,J,1)
     *         -VICE(I,J-1,1)-VICE(I-1,J-1,1))
          E12(I,J)=0.5*(0.5/DYT(J)*(UICE(I,J,1)+UICE(I-1,J,1)-
     *         UICE(I,J-1,1)-UICE(I-1,J-1,1))+0.5/(DXT(I)*CST(J))*
     *         (VICE(I,J,1)+VICE(I,J-1,1)-VICE(I-1,J,1)-VICE(I-1,J-1,1))
     *         +0.25*(UICE(I,J,1)+UICE(I-1,J,1)+UICE(I-1,J-1,1)+UICE(I,J
     *         -1,1))*TNGT(J)/RADIUS)
C NOW EVALUATE VISCOSITIES
          DELT=(E11(I,J)**2+E22(I,J)**2)*(1.0+ECM2)+4.0*ECM2*E12(I,J)**2
     *         +2.0*E11(I,J)*E22(I,J)*(1.0-ECM2)
          DELT1=SQRT(DELT)
          DELT1=MAX(GMIN,DELT1)
          ZETA(I,J)=0.5*PRESS(I,J)/DELT1
        END DO
      END DO

C NOW PUT MIN AND MAX VISCOSITIES IN
      DO J=1,NY1
        DO I=1,NX1
          ZETA(I,J)=MIN(ZMAX(I,J),ZETA(I,J))
          ZETA(I,J)=MAX(ZMIN(I,J),ZETA(I,J))
        END DO
      END DO

      AAA=0.0
      DO I=2,NX1-1
        AAA=AAA+ZETA(I,NY1-1)
      END DO
      AAA=AAA/FLOAT(NX1-2)
      DO I=1,NX1
        ZETA(I,NY1)=AAA
      END DO

      DO J=1,NY1
        ZETA(1,J)=ZETA(NX1-1,J)
        ZETA(NX1,J)=ZETA(2,J)
      END DO

      DO J=1,NY1
        DO I=1,NX1
          ETA(I,J)=ECM2*ZETA(I,J)
c         E11(I,J)=E11(I,J)*HEFFM(I,J)
c         E22(I,J)=E22(I,J)*HEFFM(I,J)
c         E12(I,J)=E12(I,J)*HEFFM(I,J)
c         SS11=(ZETA(I,J)-ETA(I,J))*(E11(I,J)+E22(I,J))-PRESS(I,J)*0.5
        END DO
      END DO
      RETURN
      END SUBROUTINE PLAST

      SUBROUTINE RELAX
!@sum  RELAX calculates ice dynamics relaxation method
!@auth Jiping Liu/Gavin Schmidt (based on code from J. Zhang)
!@ver  1.0
      IMPLICIT NONE
      REAL*8, DIMENSION(NX1,NY1) :: AU,BU,CU,FXY,FXY1
      REAL*8, DIMENSION(NY1,NX1) :: AV,BV,CV,FXYa,FXY1a
      REAL*8, DIMENSION(NX1) :: CUU,URT         !CUU,
      REAL*8, DIMENSION(NY1) :: CVV,VRT,U_tmp   !CVV,
      REAL*8, PARAMETER :: BYRAD2 = 1./(RADIUS*RADIUS)
      INTEGER I,J,J1,J2,IMD,JMD
      REAL*8 DELXY,DELXR,DELX2,DELY2,DELYR,ETAMEAN,ZETAMEAN,AA1,AA2
     *     ,AA3,AA4,AA5,AA6,AA9

      DO J=1,NY1
        DO I=1,NX1
          FORCEX(I,J)=FORCEX(I,J)*UVM(I,J)
          FORCEY(I,J)=FORCEY(I,J)*UVM(I,J)
        END DO
      END DO
C MUST UPDATE HEFF BEFORE CALLING RELAX
C FIRST SET U(2)=U(1)
      DO J=1,NY1
        DO I=1,NX1
C NOW MAKE SURE BDRY PTS ARE EQUAL TO ZERO
          UICE(I,J,2)=UICE(I,J,1)
          VICE(I,J,2)=VICE(I,J,1)
          UICE(I,J,1)=UICE(I,J,3)*UVM(I,J)
          VICE(I,J,1)=VICE(I,J,3)*UVM(I,J)
        END DO
      END DO

      DO I=1,NX1/2
      UICE(I,NY1,1)=-UICEC(I+(NX1-2)/2,NY1-1)
      VICE(I,NY1,1)=-VICEC(I+(NX1-2)/2,NY1-1)
      UICE(I,NY1,3)=-UICEC(I+(NX1-2)/2,NY1-1)
      VICE(I,NY1,3)=-VICEC(I+(NX1-2)/2,NY1-1)

      UICEC(I,NY1)=-UICEC(I+(NX1-2)/2,NY1-1)
      VICEC(I,NY1)=-VICEC(I+(NX1-2)/2,NY1-1)
      END DO

      DO I=NX1/2+1,NX1-1
      UICE(I,NY1,1)=-UICEC(I-(NX1-2)/2,NY1-1)
      VICE(I,NY1,1)=-VICEC(I-(NX1-2)/2,NY1-1)
      UICE(I,NY1,3)=-UICEC(I-(NX1-2)/2,NY1-1)
      VICE(I,NY1,3)=-VICEC(I-(NX1-2)/2,NY1-1)

      UICEC(I,NY1)=-UICEC(I-(NX1-2)/2,NY1-1)
      VICEC(I,NY1)=-VICEC(I-(NX1-2)/2,NY1-1)
      END DO

      DO J=1,NY1
      UICE(1,J,1)=UICEC(NX1-1,J)
      VICE(1,J,1)=VICEC(NX1-1,J)
      UICE(NX1,J,1)=UICEC(2,J)
      VICE(NX1,J,1)=VICEC(2,J)
      UICE(1,J,3)=UICEC(NX1-1,J)
      VICE(1,J,3)=VICEC(NX1-1,J)
      UICE(NX1,J,3)=UICEC(2,J)
      VICE(NX1,J,3)=VICEC(2,J)

      UICEC(1,J)=UICEC(NX1-1,J)
      VICEC(1,J)=VICEC(NX1-1,J)
      UICEC(NX1,J)=UICEC(2,J)
      VICEC(NX1,J)=VICEC(2,J)
      END DO

C FIRST DO UICE
C THE FIRST HALF

      DO J=2,NYPOLE
      DO I=2,NXLCYC
      DELXY=BYDXDY(I,J)    ! 0.5/(DXU(I)*DYU(J))
      DELXR=BYDXR(I)       ! 0.5/(DXU(I)*RADIUS)
      DELX2=BYDX2(I)       ! 0.5/(DXU(I)*DXU(I))
      ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))
      ZETAMEAN=0.25*(ZETA(I,J+1)+ZETA(I+1,J+1)+ZETA(I,J)+ZETA(I+1,J))

      FXY(I,J)=DRAGA(I,J)*VICEC(I,J)+FORCEX(I,J)
     3+0.5*(ZETA(I+1,J+1)*(VICEC(I+1,J+1)+VICEC(I,J+1)
     3-VICEC(I+1,J)-VICEC(I,J))+ZETA(I+1,J)*(VICEC(I+1,J)
     3+VICEC(I,J)-VICEC(I+1,J-1)-VICEC(I,J-1))+ZETA(I,J+1)
     3*(VICEC(I,J)+VICEC(I-1,J)-VICEC(I,J+1)-VICEC(I-1,J+1))
     3+ZETA(I,J)*(VICEC(I,J-1)+VICEC(I-1,J-1)-VICEC(I,J)
     3-VICEC(I-1,J)))*DELXY*BYCSU(J)
     3
     4-0.5*(ETA(I+1,J+1)*(VICEC(I+1,J+1)+VICEC(I,J+1)
     4-VICEC(I+1,J)-VICEC(I,J))+ETA(I+1,J)*(VICEC(I+1,J)
     4+VICEC(I,J)-VICEC(I+1,J-1)-VICEC(I,J-1))+ETA(I,J+1)
     4*(VICEC(I,J)+VICEC(I-1,J)-VICEC(I,J+1)-VICEC(I-1,J+1))
     4+ETA(I,J)*(VICEC(I,J-1)+VICEC(I-1,J-1)-VICEC(I,J)
     4-VICEC(I-1,J)))*DELXY*BYCSU(J)
     4
     5+0.5*(VICEC(I+1,J)-VICEC(I-1,J))*(ETA(I,J+1)+ETA(I+1,J+1)
     5-ETA(I,J)-ETA(I+1,J))*DELXY*BYCSU(J)+0.5*ETAMEAN*((VICEC(I+1,J+1)
     5-VICEC(I-1,J+1))*BYCSU(J+1)-(VICEC(I+1,J-1)-VICEC(I-1,J-1))
     5*BYCSU(J-1))*DELXY
     5
     6-((ZETA(I+1,J+1)+ZETA(I+1,J)-ZETA(I,J)-ZETA(I,J+1))
     6+(ETA(I+1,J+1)+ETA(I+1,J)-ETA(I,J)-ETA(I,J+1)))
     6*TNG(J)*VICEC(I,J)*DELXR*BYCSU(J)
     6-(ETAMEAN+ZETAMEAN)*TNG(J)*(VICEC(I+1,J)-VICEC(I-1,J))
     6*DELXR*BYCSU(J)
     6
     7-ETAMEAN*2.0*TNG(J)*(VICEC(I+1,J)-VICEC(I-1,J))*DELXR*BYCSU(J)

      END DO
      END DO

      DO J=2,NYPOLE
      DO I=2,NXLCYC
      DELX2=BYDX2(I)       ! 0.5/(DXU(I)*DXU(I))
      DELY2=BYDY2(J)       ! 0.5/(DYU(J)*DYU(J))
      DELYR=BYDYR(J)       ! 0.5/(DYU(J)*RADIUS)
      ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))
      AA1=((ETA(I+1,J)  +ZETA(I+1,J)  )*BYCSU(J)+
     *     (ETA(I+1,J+1)+ZETA(I+1,J+1))*BYCSU(J))*BYCSU(J)
      AA2=((ETA(I,J)+ZETA(I,J))*BYCSU(J)+(ETA(I,J+1)+ZETA(I,J+1))
     &     *BYCSU(J))*BYCSU(J)
      AA3=ETA(I,J+1)+ETA(I+1,J+1)
      AA4=ETA(I,J)+ETA(I+1,J)
      AA5=-(ETA(I,J+1)+ETA(I+1,J+1)-ETA(I,J)-ETA(I+1,J))*TNG(J)
      AA6=2.0*ETAMEAN*TNG(J)*TNG(J)
      AU(I,J)=-AA2*DELX2*UVM(I,J)
      BU(I,J)=((AA1+AA2)*DELX2+AA6*BYRAD2
     &+AMASS(I,J)*BYDTS*2.0+DRAGS(I,J))*UVM(I,J)+(1.0-UVM(I,J))
      CU(I,J)=-AA1*DELX2*UVM(I,J)
      END DO
      END DO
      DO J=2,NYPOLE
      AU(2,J)=0.0
      CU(NXLCYC,J)=0.0
c      CU(2,J)=CU(2,J)/BU(2,J)  ! absorbed into TRIDIAG
      END DO

      DO 1200 J=2,NYPOLE
      DO I=2,NXLCYC

      DELX2=BYDX2(I)       ! 0.5/(DXU(I)*DXU(I))
      DELY2=BYDY2(J)       ! 0.5/(DYU(J)*DYU(J))
      DELYR=BYDYR(J)       ! 0.5/(DYU(J)*RADIUS)
      ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))

      AA1=((ETA(I+1,J)  +ZETA(I+1,J)  )*BYCSU(J)+
     *     (ETA(I+1,J+1)+ZETA(I+1,J+1))*BYCSU(J))*BYCSU(J)
      AA2=((ETA(I,J)+ZETA(I,J))*BYCSU(J)+(ETA(I,J+1)+ZETA(I,J+1))
     &     *BYCSU(J))*BYCSU(J)
      AA3=ETA(I,J+1)+ETA(I+1,J+1)
      AA4=ETA(I,J)+ETA(I+1,J)
      AA5=-(ETA(I,J+1)+ETA(I+1,J+1)-ETA(I,J)-ETA(I+1,J))*TNG(J)
      AA6=2.0*ETAMEAN*TNG(J)*TNG(J)

      IF(I.EQ.2) THEN
      AA9=AA2*DELX2*UICEC(I-1,J)*UVM(I,J)*FLOAT(LCYC-0)
      ELSE IF(I.EQ.NXLCYC) THEN
      AA9=AA1*DELX2*UICEC(I+1,J)*UVM(I,J)*FLOAT(LCYC-0)
      ELSE
      AA9=0.0
      END IF

      URT(I)=AA9+FXY(I,J)-AA5*DELYR*UICE(I,J,2)
     1-(AA3+AA4)*DELY2*UICE(I,J,2)
     1+(ETA(I,J+1)+ETA(I+1,J+1))*UICE(I,J+1,2)*DELY2
     2+(ETA(I,J)+ETA(I+1,J))*UICE(I,J-1,2)*DELY2
     3+ETAMEAN*DELYR*(UICE(I,J+1,2)*TNG(J+1)-UICE(I,J-1,2)*TNG(J-1))
     4-ETAMEAN*DELYR*2.0*TNG(J)*(UICE(I,J+1,2)-UICE(I,J-1,2))
      URT(I)=(URT(I)+AMASS(I,J)*BYDTS*UICE(I,J,2)*2.0)*UVM(I,J)
      END DO

      CALL TRIDIAG(AU(2,J),BU(2,J),CU(2,J),URT(2),UICE(2,J,1),NXLCYC-1)

c      DO I=2,NXLCYC
c       CUU(I)=CU(I,J)
c      END DO
c      URT(2)=URT(2)/BU(2,J)
c      DO I=3,NXLCYC
c       IMD=I-1
c       CUU(I)=CUU(I)/(BU(I,J)-AU(I,J)*CUU(IMD))
c       URT(I)=(URT(I)-AU(I,J)*URT(IMD))/(BU(I,J)-AU(I,J)*CUU(IMD))
c      END DO
c      DO I=1,NXLCYC-2
c       J1=NXLCYC-I
c       J2=J1+1
c       URT(J1)=URT(J1)-CUU(J1)*URT(J2)
c      END DO
c      DO I=2,NXLCYC
c       UICE(I,J,1)=URT(I)
c      END DO

 1200 CONTINUE

      DO I=2,NXLCYC
      DO J=2,NYPOLE
      UICE(I,J,3)=UICE(I,J,1)
      END DO
      END DO

C NOW THE SECOND HALF
      DO I=2,NXLCYC
      DO J=2,NYPOLE
      DELX2=BYDX2(I)       ! 0.5/(DXU(I)*DXU(I))
      DELY2=BYDY2(J)       ! 0.5/(DYU(J)*DYU(J))
      DELYR=BYDYR(J)       ! 0.5/(DYU(J)*RADIUS)
      ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))
      ZETAMEAN=0.25*(ZETA(I,J+1)+ZETA(I+1,J+1)+ZETA(I,J)+ZETA(I+1,J))

      AA1=ETA(I,J+1)+ETA(I+1,J+1)
      AA2=ETA(I,J)+ETA(I+1,J)
      AA5=-(ETA(I,J+1)+ETA(I+1,J+1)-ETA(I,J)-ETA(I+1,J))*TNG(J)
      AA6=2.0*ETAMEAN*TNG(J)*TNG(J)

      AV(J,I)=(-AA2*DELY2+ETAMEAN*DELYR*(TNG(J-1)-2.0*TNG(J)))*UVM(I,J)
      BV(J,I)=((AA1+AA2)*DELY2+AA5*DELYR+AA6*BYRAD2
     &+AMASS(I,J)*BYDTS*2.0+DRAGS(I,J))*UVM(I,J)+(1.0-UVM(I,J))
      CV(J,I)=(-AA1*DELY2-ETAMEAN*DELYR*(TNG(J+1)-2.0*TNG(J)))*UVM(I,J)
      END DO
      END DO
      DO I=2,NXLCYC
      AV(2,I)=0.0
      CV(NYPOLE,I)=0.0
c      CV(2,I)=CV(2,I)/BV(2,I)  ! absorbed into TRIDIAG
      END DO

      DO I=2,NXLCYC
      DO J=2,NYPOLE
      DELX2=BYDX2(I)       ! 0.5/(DXU(I)*DXU(I))
      DELY2=BYDY2(J)       ! 0.5/(DYU(J)*DYU(J))
      DELYR=BYDYR(J)       ! 0.5/(DYU(J)*RADIUS)
      ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))

      AA1=((ETA(I+1,J)  +ZETA(I+1,J)  )*BYCSU(J)+
     *     (ETA(I+1,J+1)+ZETA(I+1,J+1))*BYCSU(J))*BYCSU(J)
      AA2=((ETA(I,J)+ZETA(I,J))*BYCSU(J)+(ETA(I,J+1)+ZETA(I,J+1))
     &     *BYCSU(J))*BYCSU(J)

      IF(J.EQ.NYPOLE) THEN
      AA9=( (ETA(I,J+1)+ETA(I+1,J+1))*DELY2*UICEC(I,J+1)
     &     +ETAMEAN*DELYR*(TNG(J+1)-2.0*TNG(J))*UICEC(I,J+1) )*UVM(I,J)
     &*FLOAT(NPOL-0)
      ELSE
      AA9=0.0
      END IF

      FXY1a(J,I)=AA9+AMASS(I,J)*BYDTS*UICE(I,J,1)*2.0
     5-(AA1+AA2)*DELX2*UICE(I,J,1)
     6+((ETA(I+1,J)+ZETA(I+1,J)+ETA(I+1,J+1)+ZETA(I+1,J+1))
     6*UICE(I+1,J,1)
     6+(ETA(I,J)+ZETA(I,J)+ETA(I,J+1)+ZETA(I,J+1))*UICE(I-1,J,1))
     6*DELX2*BYCSU(J)*BYCSU(J)

      END DO
      END DO

      DO 1300 I=2,NXLCYC
      DO J=2,NYPOLE
      VRT(J)=FXY(I,J)+FXY1a(J,I)
      VRT(J)=VRT(J)*UVM(I,J)
      END DO

      CALL TRIDIAG(AV(2,I),BV(2,I),CV(2,I),VRT(2),U_tmp(2),NYPOLE-1)

c      DO J=2,NYPOLE
c      CVV(J)=CV(J,I)
c      END DO
c      VRT(2)=VRT(2)/BV(2,I)
c      DO J=3,NYPOLE
c      JMD=J-1
c      CVV(J)=CVV(J)/(BV(J,I)-AV(J,I)*CVV(JMD))
c      VRT(J)=(VRT(J)-AV(J,I)*VRT(JMD))/(BV(J,I)-AV(J,I)*CVV(JMD))
c      END DO
c      DO J=1,NYPOLE-2
c      J1=NYPOLE-J
c      J2=J1+1
c      VRT(J1)=VRT(J1)-CVV(J1)*VRT(J2)
c      END DO
      DO J=2,NYPOLE
        UICE(I,J,1)=U_tmp(J)     ! VRT(J)   !
      END DO
 1300 CONTINUE

C NOW DO VICE
C THE FIRST HALF

      DO I=2,NXLCYC
      DO J=2,NYPOLE
      DELXY=BYDXDY(I,J)    ! 0.5/(DXU(I)*DYU(J))
      DELXR=BYDXR(I)       ! 0.5/(DXU(I)*RADIUS)
      DELY2=BYDY2(J)       ! 0.5/(DYU(J)*DYU(J))
      DELYR=BYDYR(J)       ! 0.5/(DYU(J)*RADIUS)
      ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))
      ZETAMEAN=0.25*(ZETA(I,J+1)+ZETA(I+1,J+1)+ZETA(I,J)+ZETA(I+1,J))

      FXYa(J,I)=-DRAGA(I,J)*UICEC(I,J)+FORCEY(I,J)
     3+(0.5*(UICEC(I+1,J)-UICEC(I-1,J))*(ZETA(I,J+1)+ZETA(I+1,J+1)
     3-ZETA(I,J)-ZETA(I+1,J))*DELXY*BYCSU(J)+0.5*ZETAMEAN*
     3((UICEC(I+1,J+1)
     3-UICEC(I-1,J+1))*BYCSU(J+1)-(UICEC(I+1,J-1)-UICEC(I-1,J-1))
     3*BYCSU(J-1))*DELXY)
     3
     4-(0.5*(UICEC(I+1,J)-UICEC(I-1,J))*(ETA(I,J+1)+ETA(I+1,J+1)
     4-ETA(I,J)-ETA(I+1,J))*DELXY*BYCSU(J)+0.5*ETAMEAN*((UICEC(I+1,J+1)
     4-UICEC(I-1,J+1))*BYCSU(J+1)-(UICEC(I+1,J-1)-UICEC(I-1,J-1))
     4*BYCSU(J-1))*DELXY)
     4
     5+0.5*(ETA(I+1,J+1)*(UICEC(I+1,J+1)+UICEC(I,J+1)
     5-UICEC(I+1,J)-UICEC(I,J))+ETA(I+1,J)*(UICEC(I+1,J)
     5+UICEC(I,J)-UICEC(I+1,J-1)-UICEC(I,J-1))+ETA(I,J+1)
     5*(UICEC(I,J)+UICEC(I-1,J)-UICEC(I,J+1)-UICEC(I-1,J+1))
     5+ETA(I,J)*(UICEC(I,J-1)+UICEC(I-1,J-1)-UICEC(I,J)
     5-UICEC(I-1,J)))*DELXY*BYCSU(J)
     5
     6+(ETA(I+1,J+1)+ETA(I+1,J)-ETA(I,J)-ETA(I,J+1))
     6*TNG(J)*UICEC(I,J)*DELXR*BYCSU(J)
     6+ETAMEAN*TNG(J)*(UICEC(I+1,J)-UICEC(I-1,J))*DELXR*BYCSU(J)
     6
     7+ETAMEAN*2.0*TNG(J)*(UICEC(I+1,J)-UICEC(I-1,J))*DELXR*BYCSU(J)

      END DO
      END DO

      DO I=2,NXLCYC
      DO J=2,NYPOLE
      DELX2=BYDX2(I)       ! 0.5/(DXU(I)*DXU(I))
      DELY2=BYDY2(J)       ! 0.5/(DYU(J)*DYU(J))
      DELYR=BYDYR(J)       ! 0.5/(DYU(J)*RADIUS)
      ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))
      ZETAMEAN=0.25*(ZETA(I,J+1)+ZETA(I+1,J+1)+ZETA(I,J)+ZETA(I+1,J))
      AA1=ETA(I,J+1)+ZETA(I,J+1)+ETA(I+1,J+1)+ZETA(I+1,J+1)
      AA2=ETA(I,J)+ZETA(I,J)+ETA(I+1,J)+ZETA(I+1,J)
      AA3=(ETA(I+1,J)*BYCSU(J)+ETA(I+1,J+1)*BYCSU(J))*BYCSU(J)
      AA4=(ETA(I,J)*BYCSU(J)+ETA(I,J+1)*BYCSU(J))*BYCSU(J)
      AA5=((ZETA(I,J+1)-ETA(I,J+1))+(ZETA(I+1,J+1)-ETA(I+1,J+1))
     &-(ZETA(I,J)-ETA(I,J))-(ZETA(I+1,J)-ETA(I+1,J)))*TNG(J)
      AA6=2.0*ETAMEAN*TNG(J)*TNG(J)

      AV(J,I)=(-AA2*DELY2-(ZETAMEAN-ETAMEAN)*TNG(J-1)*DELYR
     &-ETAMEAN*2.0*TNG(J)*DELYR)*UVM(I,J)
      BV(J,I)=((AA1+AA2)*DELY2+AA5*DELYR+AA6*BYRAD2
     &+AMASS(I,J)*BYDTS*2.0+DRAGS(I,J))*UVM(I,J)+(1.0-UVM(I,J))
      CV(J,I)=(-AA1*DELY2+(ZETAMEAN-ETAMEAN)*TNG(J+1)*DELYR
     &+ETAMEAN*2.0*TNG(J)*DELYR)*UVM(I,J)
      END DO
      END DO
      DO I=2,NXLCYC
      AV(2,I)=0.0
      CV(NYPOLE,I)=0.0
c      CV(2,I)=CV(2,I)/BV(2,I)  ! absorbed into TRIDIAG
      END DO

      DO 1301 I=2,NXLCYC
      DO J=2,NYPOLE
      DELX2=BYDX2(I)       ! 0.5/(DXU(I)*DXU(I))
      DELY2=BYDY2(J)       ! 0.5/(DYU(J)*DYU(J))
      DELYR=BYDYR(J)       ! 0.5/(DYU(J)*RADIUS)

      ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))
      ZETAMEAN=0.25*(ZETA(I,J+1)+ZETA(I+1,J+1)+ZETA(I,J)+ZETA(I+1,J))

      AA1=ETA(I,J+1)+ZETA(I,J+1)+ETA(I+1,J+1)+ZETA(I+1,J+1)
      AA2=ETA(I,J)+ZETA(I,J)+ETA(I+1,J)+ZETA(I+1,J)
      AA3=(ETA(I+1,J)*BYCSU(J)+ETA(I+1,J+1)*BYCSU(J))*BYCSU(J)
      AA4=(ETA(I,J)*BYCSU(J)+ETA(I,J+1)*BYCSU(J))*BYCSU(J)
      AA5=((ZETA(I,J+1)-ETA(I,J+1))+(ZETA(I+1,J+1)-ETA(I+1,J+1))
     &-(ZETA(I,J)-ETA(I,J))-(ZETA(I+1,J)-ETA(I+1,J)))*TNG(J)
      AA6=2.0*ETAMEAN*TNG(J)*TNG(J)

      IF(J.EQ.NYPOLE) THEN
      AA9=(AA1*DELY2-(ZETAMEAN-ETAMEAN)*TNG(J+1)*DELYR
     &-ETAMEAN*2.0*TNG(J)*DELYR)*VICEC(I,J+1)*UVM(I,J)*FLOAT(NPOL-0)
      ELSE
      AA9=0.0
      END IF

      VRT(J)=AA9+FXYa(J,I)-(AA3+AA4)*DELX2*VICE(I,J,2)
     6+((ETA(I+1,J)*BYCSU(J)+ETA(I+1,J+1)*BYCSU(J))*VICE(I+1,J,2)*DELX2
     7     +(ETA(I,J)*BYCSU(J)+ETA(I,J+1)*BYCSU(J))*VICE(I-1,J,2)*DELX2)
     *     *BYCSU(J)
      VRT(J)=(VRT(J)+AMASS(I,J)*BYDTS*VICE(I,J,2)*2.0)*UVM(I,J)
      END DO

      CALL TRIDIAG(AV(2,I),BV(2,I),CV(2,I),VRT(2),U_tmp(2),NYPOLE-1)

c      DO J=2,NYPOLE
c      CVV(J)=CV(J,I)
c      END DO
c      VRT(2)=VRT(2)/BV(2,I)
c      DO J=3,NYPOLE
c      JMD=J-1
c      CVV(J)=CVV(J)/(BV(J,I)-AV(J,I)*CVV(JMD))
c      VRT(J)=(VRT(J)-AV(J,I)*VRT(JMD))/(BV(J,I)-AV(J,I)*CVV(JMD))
c      END DO
c      DO J=1,NYPOLE-2
c      J1=NYPOLE-J
c      J2=J1+1
c      VRT(J1)=VRT(J1)-CVV(J1)*VRT(J2)
c      END DO
      DO J=2,NYPOLE
        VICE(I,J,1)=U_tmp(J)   ! VRT(J)   !
      END DO
 1301 CONTINUE

      DO I=2,NXLCYC
      DO J=2,NYPOLE
      VICE(I,J,3)=VICE(I,J,1)
      END DO
      END DO

C NOW THE SECOND HALF

      DO J=2,NYPOLE
      DO I=2,NXLCYC
      DELX2=BYDX2(I)       ! 0.5/(DXU(I)*DXU(I))
      DELY2=BYDY2(J)       ! 0.5/(DYU(J)*DYU(J))
      DELYR=BYDYR(J)       ! 0.5/(DYU(J)*RADIUS)
      ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))
      ZETAMEAN=0.25*(ZETA(I,J+1)+ZETA(I+1,J+1)+ZETA(I,J)+ZETA(I+1,J))

      AA1=ETA(I,J+1)+ZETA(I,J+1)+ETA(I+1,J+1)+ZETA(I+1,J+1)
      AA2=ETA(I,J)+ZETA(I,J)+ETA(I+1,J)+ZETA(I+1,J)
      AA3=(ETA(I+1,J)*BYCSU(J)+ETA(I+1,J+1)*BYCSU(J))*BYCSU(J)
      AA4=(ETA(I,J)*BYCSU(J)+ETA(I,J+1)*BYCSU(J))*BYCSU(J)
      AA6=2.0*ETAMEAN*TNG(J)*TNG(J)

      AU(I,J)=-AA4*DELX2*UVM(I,J)
      BU(I,J)=((AA3+AA4)*DELX2+AA6*BYRAD2
     &+AMASS(I,J)*BYDTS*2.0+DRAGS(I,J))*UVM(I,J)+(1.0-UVM(I,J))
      CU(I,J)=-AA3*DELX2*UVM(I,J)
      END DO
      END DO
      DO J=2,NYPOLE
      AU(2,J)=0.0
      CU(NXLCYC,J)=0.0
c      CU(2,J)=CU(2,J)/BU(2,J)   ! absorbed into TRIDIAG
      END DO

      DO J=2,NYPOLE
      DO I=2,NXLCYC

      DELX2=BYDX2(I)       ! 0.5/(DXU(I)*DXU(I))
      DELY2=BYDY2(J)       ! 0.5/(DYU(J)*DYU(J))
      DELYR=BYDYR(J)       ! 0.5/(DYU(J)*RADIUS)
      ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))
      ZETAMEAN=0.25*(ZETA(I,J+1)+ZETA(I+1,J+1)+ZETA(I,J)+ZETA(I+1,J))

      AA1=ETA(I,J+1)+ZETA(I,J+1)+ETA(I+1,J+1)+ZETA(I+1,J+1)
      AA2=ETA(I,J)+ZETA(I,J)+ETA(I+1,J)+ZETA(I+1,J)
      AA3=(ETA(I+1,J)*BYCSU(J)+ETA(I+1,J+1)*BYCSU(J))*BYCSU(J)
      AA4=(ETA(I,J)*BYCSU(J)+ETA(I,J+1)*BYCSU(J))*BYCSU(J)
      AA5=((ZETA(I,J+1)-ETA(I,J+1))+(ZETA(I+1,J+1)-ETA(I+1,J+1))
     &-(ZETA(I,J)-ETA(I,J))-(ZETA(I+1,J)-ETA(I+1,J)))*TNG(J)
      AA6=2.0*ETAMEAN*TNG(J)*TNG(J)

      IF(I.EQ.2) THEN
      AA9=AA4*DELX2*VICEC(I-1,J)*UVM(I,J)*FLOAT(LCYC-0)
      ELSE IF(I.EQ.NXLCYC) THEN
      AA9=AA3*DELX2*VICEC(I+1,J)*UVM(I,J)*FLOAT(LCYC-0)
      ELSE
      AA9=0.0
      END IF
      FXY1(I,J)=AA9+AMASS(I,J)*BYDTS*VICE(I,J,1)*2.0
     1-AA5*DELYR*VICE(I,J,1)
     1-(AA1+AA2)*DELY2*VICE(I,J,1)
     1+AA1*DELY2*VICE(I,J+1,1)-((ZETAMEAN-ETAMEAN)*TNG(J+1)*DELYR
     1+ETAMEAN*2.0*TNG(J)*DELYR)*VICE(I,J+1,1)
     2+AA2*DELY2*VICE(I,J-1,1)+((ZETAMEAN-ETAMEAN)*TNG(J-1)*DELYR
     2+ETAMEAN*2.0*TNG(J)*DELYR)*VICE(I,J-1,1)
      END DO
      END DO

      DO 1201 J=2,NYPOLE
      DO I=2,NXLCYC
      URT(I)=FXYa(J,I)+FXY1(I,J)
      URT(I)=URT(I)*UVM(I,J)
      END DO

      CALL TRIDIAG(AU(2,J),BU(2,J),CU(2,J),URT(2),VICE(2,J,1),NXLCYC-1)

c      DO I=2,NXLCYC
c      CUU(I)=CU(I,J)
c      END DO
c      URT(2)=URT(2)/BU(2,J)
c      DO I=3,NXLCYC
c      IMD=I-1
c      CUU(I)=CUU(I)/(BU(I,J)-AU(I,J)*CUU(IMD))
c      URT(I)=(URT(I)-AU(I,J)*URT(IMD))/(BU(I,J)-AU(I,J)*CUU(IMD))
c      END DO
c      DO I=1,NXLCYC-2
c      J1=NXLCYC-I
c      J2=J1+1
c      URT(J1)=URT(J1)-CUU(J1)*URT(J2)
c      END DO
c      DO I=2,NXLCYC
c      VICE(I,J,1)=URT(I)
c      END DO
 1201 CONTINUE

      DO J=2,NYPOLE
      DO I=2,NXLCYC
      UICE(I,J,1)=UICE(I,J,1)*UVM(I,J)
      VICE(I,J,1)=VICE(I,J,1)*UVM(I,J)
      END DO
      END DO

      RETURN
      END SUBROUTINE RELAX

      SUBROUTINE setup_icedyn_grid
      USE MODEL_COM, only : dts=>dtsrc
      IMPLICIT NONE
      REAL*8 :: dlat,dlon,phit,phiu,hemi,rms
      INTEGER I,J,n,k,kki,sumk,l
C****
C**** calculate grid and initialise arrays
c****
      dlat=nint(180./(ny1-1))*radian
      dlon=nint(360./(nx1-2))*radian
      bydts = 1./dts
c****
      do j = 1,ny1
        dyt(j) = dlat*radius
        dyu(j) = dlat*radius
      enddo
      do i=1,nx1
        dxt(i) = dlon*radius
        dxu(i) = dlon*radius
      enddo
      dxt(1) = dxt(nx1-1)
      dxt(nx1) = dxt(2)
      dxu(1) = dxu(nx1-1)
      dxu(nx1) = dxu(2)

      do i=1,nx1
        bydx2(i)=0.5/(dxu(i)*dxu(i))
        bydxr(i)=0.5/(dxu(i)*radius)
      end do
      do j=1,ny1
        bydy2(j)=0.5/(dyu(j)*dyu(j))
        bydyr(j)=0.5/(dyu(j)*radius)
        do i=1,nx1
          bydxdy(i,j) = 0.5/(dxu(i)*dyu(j))
        end do
      end do

      do j = 1,ny1
       phit = (-90.+(j-1)*4.)*radian
       phiu = (-88.+(j-1)*4.)*radian
       cst(j) = cos(phit)
       csu(j) = cos(phiu)
       bycsu(j) = 1./csu(j)
       tng(j) = sin(phiu)/csu(j)
       TNGT(J)=SIN(PHIT)/CST(J)
       DO I=1,NX1
         SINEN(I,J)=SIN(PHIU)
       enddo
      enddo
      TNGT(NY1)=TNGT(NY1-1)
      TNG(NY1)=TNG(NY1-1)
      CSU(NY1)=CSU(NY1-1)
      bycsu(NY1) = 1./csu(NY1)

C**** sin/cos ice-ocean turning angle
      SINWAT=SIN(OIPHI)
      COSWAT=COS(OIPHI)

C**** Set land masks for tracer and velocity points
      do j=1,ny1
        do i=2,nx1-1
          heffm(i,j)=nint(focean(i-1,j))
        enddo
        heffm(1,j)=heffm(nx1-1,j)
        heffm(nx1,j)=heffm(2,j)  
      enddo
C**** define velocity points (including exterior corners)
      do j=1,ny1-1
        do i=1,nx1-1
c          sumk=heffm(i,j)+heffm(i+1,j)+heffm(i,j+1)+heffm(i+1,j+1)
c          if (sumk.ge.3) uvm(i,j)=1  ! includes exterior corners
          uvm(i,j) = nint(min(heffm(i,j), heffm(i+1,j), heffm(i,j+1),
     *         heffm(i+1,j+1)))
        end do
      end do
C**** reset tracer points to surround velocity points (except for single
      do j=2,ny1-1
        do i=2,nx1-1
          k = nint(max (uvm(i,j), uvm(i-1,j), uvm(i,j-1), uvm(i-1,j-1)))
c         sumk = nint(uvm(i,j)+uvm(i+1,j)+uvm(i,j+1)+uvm(i+1,j+1))
c set to k except if an island 
c         if (.not. (sumk.eq.4.and.focean(i-1,j).eq.0) ) then
            heffm(i,j) = k
c         end if
        enddo
      enddo
C**** final sweep to reinstate islands
c      do j=2,ny1-1
c        do i=2,nx1-1
c          sumk = nint(uvm(i,j)+uvm(i+1,j)+uvm(i,j+1)+uvm(i+1,j+1))
c          if (sumk.eq.4.and.heffm(i,j).eq.0.) then
c            uvm(i,j)=0 ; uvm(i+1,j)=0 ; uvm(i,j+1)=0 ; uvm(i+1,j+1)=0
c          end if
c        enddo
c      enddo
c set lateral boundary conditions
      do j=1,ny1
        heffm(1,j)   = heffm(nx1-1,j)
        heffm(nx1,j) = heffm(2,j)
      enddo
      do j=1,ny1-1
        do i=1,nx1-1
          uvm(i,j) = nint(min(heffm(i,j), heffm(i+1,j), heffm(i,j+1),
     *         heffm(i+1,j+1)))
        end do
      end do
c set cyclic conditions on eastern and western boundary
      do j=1,ny1
        uvm(1,j) = uvm(nx1-1,j)
        uvm(nx1,j) = uvm(2,j)
      enddo

      RETURN
      END SUBROUTINE setup_icedyn_grid

      END MODULE ICEDYN

      SUBROUTINE VPICEDYN
!@sum  vpicedyn is the entry point into the viscous-plastic ice 
!@+    dynamics code
!@auth Gavin Schmidt (based on code from J. Zhang)
      USE ICEDYN, only : nx1,ny1,form,relax,uice,vice,uicec,vicec
      IMPLICIT NONE
      REAL*8, DIMENSION(NX1,NY1) :: USAVE,VSAVE
      REAL*8 rms
      INTEGER kki,i,j

      rms=0.
C KKI LOOP IS FOR PSEUDO-TIMESTEPPING
      KKI=0.
 10   KKI=KKI+1

C FIRST DO PREDICTOR
      DO J=1,NY1
      DO I=1,NX1
       UICE(I,J,3)=UICE(I,J,1)
       VICE(I,J,3)=VICE(I,J,1)
       UICEC(I,J)=UICE(I,J,1)
       VICEC(I,J)=VICE(I,J,1)
      END DO
      END DO

      CALL FORM
      CALL RELAX

      DO J=1,NY1
       UICE(1,J,1)=UICE(NX1-1,J,1)
       VICE(1,J,1)=VICE(NX1-1,J,1)
       UICE(NX1,J,1)=UICE(2,J,1)
       VICE(NX1,J,1)=VICE(2,J,1)
      END DO

C NOW DO REGULAR TIME STEP
C NOW DO MODIFIED EULER STEP
c
      DO J=1,NY1
      DO I=1,NX1
       UICE(I,J,1)=0.5*(UICE(I,J,1)+UICE(I,J,2))
       VICE(I,J,1)=0.5*(VICE(I,J,1)+VICE(I,J,2))
      END DO
      END DO

      CALL FORM

C NOW SET U(1)=U(2) AND SAME FOR V
      DO J=1,NY1
      DO I=1,NX1
       UICE(I,J,3)=UICE(I,J,1)
       VICE(I,J,3)=VICE(I,J,1)
       UICEC(I,J)=UICE(I,J,1)
       VICEC(I,J)=VICE(I,J,1)
       UICE(I,J,1)=UICE(I,J,2)
       VICE(I,J,1)=VICE(I,J,2)
      END DO
      END DO

      CALL RELAX

      DO J=1,NY1
       UICE(1,J,1)=UICE(NX1-1,J,1)
       VICE(1,J,1)=VICE(NX1-1,J,1)
       UICE(NX1,J,1)=UICE(2,J,1)
       VICE(NX1,J,1)=VICE(2,J,1)
      END DO

      if (kki.gt.1) then ! test convergence
        rms=0.
        do i=1,nx1
          do j=1,ny1
            rms=rms+(USAVE(i,j)-UICE(i,j,1))**2+(VSAVE(i,j)-VICE(i,j,1))
     *           **2
          end do
        end do
      end if

      if (kki.eq.20) then
        write(6,*) "Too many iterations in VPICEDYN. kki:",kki,rms
      elseif (kki.eq.1 .or. rms.gt.0.01d0) then
        USAVE=UICE(:,:,1)
        VSAVE=VICE(:,:,1)
        goto 10 
      end if

      RETURN
      END SUBROUTINE VPICEDYN

