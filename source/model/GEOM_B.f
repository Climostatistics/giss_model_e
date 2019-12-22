      MODULE GEOM
!@sum  GEOM contains spherical geometric variables and arrays
!@auth Original development team
!@ver  1.0 (B grid version)
!@cont GEOM_B
      USE CONSTANT, only : OMEGA,RADIUS,TWOPI,SDAY,radian
      USE MODEL_COM, only : IM,JM,LM,FIM,BYIM
      IMPLICIT NONE
C**** The primary grid is the A grid (including both poles)
C**** The secondary grid is for the B grid velocities, located on the
C**** vertices of the A grid (Note: first velocity point is J=2)
C**** Polar boxes can have different latitudinal size and are treated
C**** as though they were 1/IM of their actual area
      SAVE
!@param  DLON grid spacing in longitude (deg)
      REAL*8, PARAMETER :: DLON = TWOPI*BYIM
C**** For the wonderland model set DLON=DLON/3
c      REAL*8, PARAMETER :: DLON=TWOPI/(IM*3)
!@param  DLAT,DLAT_DG grid spacing in latitude (rad,deg)
      REAL*8  :: DLAT,DLAT_DG
!@param  FJEQ equatorial value of J
      REAL*8, PARAMETER :: FJEQ=.5*(1+JM)
!@var  J1U index of southernmost latitude (currently 2, later 1)
      INTEGER, parameter :: J1U = 2
!@var  JRANGE_HEMI lowest,highest lat index for SH,NH for A,B grid
      INTEGER, parameter, dimension(2,2,2) :: JRANGE_HEMI = reshape(
     *  (/1,JM/2,  1+JM/2,JM,  J1U,J1U-1+JM/2, J1U-1+JM/2,JM+J1U-2/),
     *  (/2,2,2/))
!@var  LAT latitude of mid point of primary grid box (radians)
      REAL*8, DIMENSION(JM) :: LAT
!@var  LAT_DG latitude of mid points of primary and sec. grid boxs (deg)
      REAL*8, DIMENSION(JM,2) :: LAT_DG
!@var  LON longitude of mid points of primary grid box (radians)
      REAL*8, DIMENSION(IM) :: LON
!@var  LON_DG longitude of mid points of prim. and sec. grid boxes (deg)
      REAL*8, DIMENSION(IM,2) :: LON_DG
!@var  DXYP,BYDXYP area of grid box (+inverse) (m^2)
C**** Note that this is not the exact area, but is what is required for
C**** some B-grid conservation quantities
      REAL*8, DIMENSION(JM) :: DXYP,BYDXYP
!@var AREAG global integral of area (m^2)
      REAL*8 :: AREAG
!@var WTJ area weighting used in JLMAP, JKMAP (for hemispheric means)
      REAL*8, DIMENSION(JM,2,2) :: WTJ

!@var DXYV,BYDXYV area of grid box around velocity point (recip.)(m^2)
      REAL*8, DIMENSION(JM) :: DXYV,BYDXYV

!@var  DXP,DYP,BYDXP,BYDYP distance between points on primary grid
!@+     (+inverse)
      REAL*8, DIMENSION(JM) :: DXP,DYP,BYDXP,BYDYP
!@var  DXV,DYV distance between velocity points (secondary grid)
      REAL*8, DIMENSION(JM) :: DXV,DYV
!@var  DXYN,DXYS half box areas to the North,South of primary grid point
      REAL*8, DIMENSION(JM) :: DXYS,DXYN
!@var  SINP sin of latitude at primary grid points
      REAL*8, DIMENSION(JM) :: SINP
!@var  COSP, COSV cos of latitude at primary, secondary latitudes
      REAL*8, DIMENSION(JM) :: COSP,COSV
!@var  RAPVS,RAPVN,RAVPS,RAVPN area scalings for primary and sec. grid
      REAL*8, DIMENSION(JM) :: RAPVS,RAPVN,RAVPS,RAVPN

!@var SINIV,COSIV,SINIP,COSIP longitud. sin,cos for wind,pressure grid
      REAL*8, DIMENSION(IM) :: SINIV,COSIV,SINIP,COSIP
!@var  RAVJ scaling for A grid U/V to B grid points (func. of lat. j)
!@var  RAPJ scaling for B grid -> A grid conversion (1/4,1/im at poles)
      REAL*8, DIMENSION(IM,JM) :: RAPJ,RAVJ
!@var  IDJJ J index of adjacent U/V points for A grid (func. of lat. j)
      INTEGER, DIMENSION(IM,JM) :: IDJJ
!@var  IDIJ I index of adjacent U/V points for A grid (func. of lat/lon)
      INTEGER, DIMENSION(IM,IM,JM) :: IDIJ
!@var  KMAXJ varying number of adjacent velocity points
      INTEGER, DIMENSION(JM) :: KMAXJ
!@var  IMAXJ varying number of used longitudes
      INTEGER, DIMENSION(JM) :: IMAXJ
!@var  FCOR latitudinally varying coriolis parameter
      REAL*8, DIMENSION(JM) :: FCOR

      CONTAINS

      SUBROUTINE GEOM_B
!@sum  GEOM_B Calculate spherical geometry for B grid
!@auth Original development team (modifications by G. Schmidt)
!@ver  1.0 (B grid version)
      IMPLICIT NONE
      REAL*8, PARAMETER :: EDPERD=1.,EDPERY = 365.

      INTEGER :: I,J,K,IM1  !@var I,J,K,IM1  loop variables
      INTEGER :: JVPO,JMHALF
      REAL*8  :: RAVPO

C**** latitudinal spacing depends on whether you have even spacing or
C**** a partial box at the pole
      DLAT_DG=180./JM                     ! even spacing (default)
      IF (JM.eq.46) DLAT_DG=180./(JM-1)   ! 1/2 box at pole for 4x5
cc    IF (JM.eq.24) DLAT_DG=180./(JM-1)   ! 1/2 box at pole, orig 8x10
      IF (JM.eq.24) DLAT_DG=180./(JM-1.5) ! 1/4 box at pole, 'real' 8x10
      DLAT=DLAT_DG*radian
      LAT(1)  = -.25*TWOPI
      LAT(JM) = -LAT(1)
      SINP(1)  = -1.
      SINP(JM) = 1.
      COSP(1)  = 0.
      COSP(JM) = 0.
      DXP(1)  = 0.
      DXP(JM) = 0.
      DO J=2,JM-1
        LAT(J)  = DLAT*(J-FJEQ)
        SINP(J) = SIN(LAT(J))
        COSP(J) = COS(LAT(J))
        DXP(J)  = RADIUS*DLON*COSP(J)
      END DO
      BYDXP(2:JM-1) = 1.D0/DXP(2:JM-1)
      DO J=2,JM
        COSV(J) = .5*(COSP(J-1)+COSP(J))
        DXV(J)  = .5*(DXP(J-1)+DXP(J))
        DYV(J)  = RADIUS*(LAT(J)-LAT(J-1))
      END DO
      DYP(1)  = RADIUS*(LAT(2)-LAT(1)-0.5*DLAT)
      DYP(JM) = RADIUS*(LAT(JM)-LAT(JM-1)-0.5*DLAT)
      DXYP(1) = .5*DXV(2)*DYP(1)
      BYDXYP(1) = 1./DXYP(1)
      DXYP(JM)= .5*DXV(JM)*DYP(JM)
      BYDXYP(JM) = 1./DXYP(JM)
      DXYS(1)  = 0.
      DXYS(JM) = DXYP(JM)
      DXYN(1)  = DXYP(1)
      DXYN(JM) = 0.
      AREAG = DXYP(1)+DXYP(JM)
      DO J=2,JM-1
        DYP(J)  = .5*(DYV(J)+DYV(J+1))
        DXYP(J) = .5*(DXV(J)+DXV(J+1))*DYP(J)
        BYDXYP(J) = 1./DXYP(J)
        DXYS(J) = .5*DXYP(J)
        DXYN(J) = .5*DXYP(J)
        AREAG = AREAG+DXYP(J)
      END DO
      BYDYP(:) = 1.D0/DYP(:)
      AREAG = AREAG*FIM
      RAVPS(1)  = 0.
      RAPVS(1)  = 0.
      RAVPN(JM) = 0.
      RAPVN(JM) = 0.
      DO J=2,JM
        DXYV(J) = DXYN(J-1)+DXYS(J)
        BYDXYV(J) = 1./DXYV(J)
        RAPVS(J)   = .5*DXYS(J)/DXYV(J)
        RAPVN(J-1) = .5*DXYN(J-1)/DXYV(J)
        RAVPS(J)   = .5*DXYS(J)/DXYP(J)
        RAVPN(J-1) = .5*DXYN(J-1)/DXYP(J-1)
      END DO
C**** LONGITUDES (degrees); used in ILMAP
      LON_DG(1,1) = -180.+360./(2.*FLOAT(IM))
      LON_DG(1,2) = -180.+360./    FLOAT(IM)
      DO I=2,IM
        LON_DG(I,1) = LON_DG(I-1,1)+360./FLOAT(IM)
        LON_DG(I,2) = LON_DG(I-1,2)+360./FLOAT(IM)
      END DO
C**** LATITUDES (degrees); used extensively in the diagn. print routines
      LAT_DG(1,1:2)=-90.
      LAT_DG(JM,1)=90.
      DO J=2,JM-1
        LAT_DG(J,1)=DLAT_DG*(J-FJEQ)    ! primary (tracer) latitudes
      END DO
      DO J=2,JM
        LAT_DG(J,2)=DLAT_DG*(J-JM/2-1)  ! secondary (velocity) latitudes
      END DO
C**** WTJ: area weighting for JKMAP, JLMAP hemispheres
      JMHALF= JM/2
      DO J=1,JM
        WTJ(J,1,1)=1.
        WTJ(J,2,1)=2.*FIM*DXYP(J)/AREAG
      END DO
      DO J=2,JM
        WTJ(J,1,2)=1.
        WTJ(J,2,2)=2.*FIM*DXYV(J)/AREAG
      END DO
      WTJ(JMHALF+1,1,2)=.5
      WTJ(JMHALF+1,2,2)=WTJ(JMHALF+1,2,2)/2.
      WTJ(1,1,2)=0.
      WTJ(1,2,2)=0.
C**** CALCULATE CORIOLIS PARAMETER
c      OMEGA = TWOPI*(EDPERD+EDPERY)/(EDPERD*EDPERY*SDAY)
      FCOR(1)  = -RADIUS*OMEGA*.5*COSP(2)*DXV(2)
      FCOR(JM) = -FCOR(1)
      DO J=2,JM-1
        FCOR(J) = OMEGA*(DXV(J)*DXV(J)-DXV(J+1)*DXV(J+1))/DLON
      END DO

C**** Set indexes and scalings for the influence of A grid points on
C**** adjacent velocity points

C**** Calculate relative directions of polar box to nearby U,V points
      DO I=1,IM
        SINIV(I)=SIN((I-1)*DLON)
        COSIV(I)=COS((I-1)*TWOPI*BYIM) ! DLON)
        LON(I)=DLON*(I-.5)
        SINIP(I)=SIN(LON(I))
        COSIP(I)=COS(LON(I))
      END DO

C**** Conditions at the poles
      DO J=1,JM,JM-1
        IF(J.EQ.1) THEN
          JVPO=2
          RAVPO=2.*RAPVN(1)
        ELSE
          JVPO=JM
          RAVPO=2.*RAPVS(JM)
        END IF
        KMAXJ(J)=IM
        IMAXJ(J)=1
        RAVJ(1:KMAXJ(J),J)=RAVPO
        RAPJ(1:KMAXJ(J),J)=BYIM
        IDJJ(1:KMAXJ(J),J)=JVPO
        DO K=1,KMAXJ(J)
          IDIJ(K,1:IM,J)=K
        END DO
      END DO
C**** Conditions at non-polar points
      DO J=2,JM-1
        KMAXJ(J)=4
        IMAXJ(J)=IM
        DO K=1,2
          RAVJ(K,J)=RAPVS(J)
          RAPJ(K,J)=RAVPS(J)    ! = .25
          IDJJ(K,J)=J
          RAVJ(K+2,J)=RAPVN(J)
          RAPJ(K+2,J)=RAVPN(J)  ! = .25
          IDJJ(K+2,J)=J+1
        END DO
        IM1=IM
        DO I=1,IM
          IDIJ(1,I,J)=IM1
          IDIJ(2,I,J)=I
          IDIJ(3,I,J)=IM1
          IDIJ(4,I,J)=I
          IM1=I
        END DO
      END DO
      
      RETURN
      END SUBROUTINE GEOM_B

      END MODULE GEOM

