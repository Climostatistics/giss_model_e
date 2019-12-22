      SUBROUTINE ICEDYN_DUM
!@sum ICEDYN_DUM dummy routines to replace ice dynamics
      ENTRY io_icedyn
      ENTRY io_icdiag
      ENTRY reset_icdiag
      ENTRY ADVSI
      ENTRY init_icedyn
      ENTRY diag_ICEDYN
      RETURN
      END SUBROUTINE ICEDYN_DUM

      SUBROUTINE DYNSI
!@sum DYNSI simple coding to estimate ice-ocean friction velocity
!@auth Gavin Schmidt
      USE CONSTANT, only : rhows
      USE MODEL_COM, only : im,jm,kocean,focean,dtsrc
      USE GEOM, only : imaxj
      USE SEAICE, only : oi_ustar0
      USE SEAICE_COM, only : rsi
      USE FLUXES, only : UI2rho,dmua,dmva

      IMPLICIT NONE
      INTEGER I,J
      REAL*8 ustar1

      IF (KOCEAN.eq.1) THEN
        DO J=1,JM
        DO I=1,IMAXJ(J)
c          UI2rho(I,J) = rhows*(oi_ustar0)**2  ! default
C**** with wind stress dependence
          if (rsi(i,j)*focean(i,j).gt.0) then
            ustar1= SQRT(SQRT(DMUA(I,J,2)**2+DMVA(I,J,2)**2)/(RSI(i,j)
     *           *focean(i,j)*DTSRC*RHOWS))
            UI2rho(I,J)=rhows*(oi_ustar0*max(1d0,1d3*ustar1))**2
          end if
        END DO
        END DO
      END IF

      RETURN
      END SUBROUTINE DYNSI

