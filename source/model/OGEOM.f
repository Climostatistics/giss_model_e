C****
      SUBROUTINE GEOMO
!@sum  GEOMO Calculates the spherical geometry for the C grid
!@auth Gary Russell
!@ver  1.0  
      USE CONSTANT, only : twopi,radius,omega
      USE GEOM, only : dxyp,bydxyp
      USE OCEAN
      IMPLICIT NONE
      INTEGER I,J
      REAL*8 FJEQ,DLON,DLAT,VLATS,SINVS,VLATN,SINVN

      DLON   = TWOPI/IM
      DLAT   = NINT(180./(JM-1))*TWOPI/360.
C**** Geometric parameters defined at secondary latitudes
      FJEQ    = .5*(1+JM)
      RLAT(1) = -TWOPI/4.
      SINPO(1) = SIN(RLAT(1))
      DO 110 J=1,JM-1
      RLAT(J+1) = DLAT*(J+1.-FJEQ)
      IF(J.GE.JM-1)  RLAT(J+1) = TWOPI/4.
      SINPO(J+1) = SIN(RLAT(J+1))
      COSVO(J) = (SINPO(J+1)-SINPO(J))/(RLAT(J+1)-RLAT(J))
      DXVO(J)  = DLON*RADIUS*COSVO(J)
      DYVO(J)  = (RLAT(J+1)-RLAT(J))*RADIUS
  110 DXYVO(J) = DLON*RADIUS*RADIUS*(SINPO(J+1)-SINPO(J))
      DXVO(JM) = 0.
      DYVO(JM) = 0.
      DXYVO(JM)= 0.
      COSVO(0) = 0.
      COSVO(JM)= 0.
C**** Geometric parameters defined at primary latitudes
      VLATS = -TWOPI/4.
      SINVS = SIN(VLATS)
      DO 120 J=1,JM
      VLATN = DLAT*(J+.5-FJEQ)
      IF(J.EQ.JM)  VLATN = TWOPI/4.
      SINVN   = SIN(VLATN)
      DXYPO(J) = DLON*RADIUS*RADIUS*(SINVN-SINVS)
      BYDXYPO(J) = 1d0/DXYPO(J)
      IF (J.gt.1) DXPO(J)  = .5*(DXVO(J-1)+DXVO(J))
      DYPO(J)  = (VLATN-VLATS)*RADIUS
      COSPO(J) = .5*(COSVO(J-1)+COSVO(J))
      DXYSO(J) = .5*DXYPO(J)
      DXYNO(J) = .5*DXYPO(J)
      VLATS   = VLATN
  120 SINVS   = SINVN
      DXPO(1)  = .5*DXVO(1)
      DXPO(JM) = .5*DXVO(JM-1)
      DXYSO(1) = 0.
      DXYSO(JM)= DXYPO(JM)
      DXYNO(1) = DXYPO(1)
      DXYNO(JM)= 0.
C**** Calculate area ratios for applying source terms to V winds
      DO 130 J=1,JM-1
      RAMVN(J  ) = DXYNO(J  )/(DXYNO(J)+DXYSO(J+1))
  130 RAMVS(J+1) = DXYSO(J+1)/(DXYNO(J)+DXYSO(J+1))
      RAMVS(1  ) = 0.
      RAMVN(JM ) = 0.
C**** Calculate area ratios for converting atmospheric grid to ocean
      DO J=1,JM
        RATOC(J)=DXYP(J)*BYDXYPO(J)
        ROCAT(J)=DXYPO(J)*BYDXYP(J)
      END DO
C**** calculate SINI/COSI terms for C grid
      DO I=1,IM
        SINIC(I) = SIN((I-.5)*TWOPI/IM)
        COSIC(I) = COS((I-.5)*TWOPI/IM)
      END DO
C**** calculate the grid box at 40S (for use with OPFIL)
      DO J=1,JM/2
        IF (RLAT(J).gt.-40.*TWOPI/360.) THEN
          J40S=J-1
          EXIT
        END IF
      END DO
C****
      RETURN
      END
