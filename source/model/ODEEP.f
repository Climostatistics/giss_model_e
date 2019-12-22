!@sum ODEEP contains routines used for Qflux mixed layer with deep diff.
!@auth G. Schmidt/G. Russell
!@ver  1.0
      MODULE ODEEP_COM
!@sum  ODEEP_COM defines the variables for deep diffusing Qflux model
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0
      USE MODEL_COM, only : im,jm
      IMPLICIT NONE
      SAVE

!@param LMOM number of layers for deep ocean diffusion
c      INTEGER, PARAMETER :: LMOM = 9    ! good for 1000m
      INTEGER, PARAMETER :: LMOM = 12    ! good for 5015m
!@var TG3M Monthly accumulation of temperatures at base of mixed layer
      REAL*8, DIMENSION(IM,JM,12) :: TG3M
!@var RTGO Temperature anomaly in thermocline
      REAL*8, DIMENSION(LMOM,IM,JM) :: RTGO
!@var STG3 accumulated temperature at base of mixed layer
      REAL*8, DIMENSION(IM,JM) :: STG3
!@var DTG3 accumulated temperature diff. from initial monthly values
      REAL*8, DIMENSION(IM,JM) :: DTG3
!@var EDO ocean vertical diffusion (m^2/s)
      REAL*8, DIMENSION(IM,JM) :: EDO
!@var DZ thermocline layer thickness (m)
      REAL*8, DIMENSION(LMOM) :: DZ
!@var DZO,BYDZO distance between centres in thermocline layer (m)
      REAL*8, DIMENSION(LMOM-1) :: DZO,BYDZO

      END MODULE ODEEP_COM

      SUBROUTINE init_ODEEP(iniOCEAN)
!@sum  init_ODEEP initialise deep ocean arrays
!@auth G. Schmidt
!@ver  1.0
      USE FILEMANAGER, only : openunit,closeunit
      USE MODEL_COM, only : im,jm
      USE ODEEP_COM, only : tg3m,stg3,dtg3,rtgo,dz,dzo,bydzo,edo,lmom
      USE STATIC_OCEAN, only : tocean
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: iniOCEAN
      INTEGER :: iu_tg3m,iu_EDDY,L
      CHARACTER*80 TITLE

!@param FAC ratio of adjacent deep ocean layers
C**** NOTE: For LMOM is 9 this value gives a total depth of 1000m
C**** For any different number of layers, the effective depth is given
C**** by the equation Z=10*(1-x^(LMOM-1))/(1-x).
C**** In particular, for LMOM=12, the total depth is 5015m
      REAL*8, PARAMETER :: FAC=1.705357255658901d0

C**** READ IN EDDY DIFFUSIVITY AT BASE OF MIXED LAYER
      CALL openunit("EDDY",iu_EDDY,.TRUE.,.TRUE.)
      CALL READT (iu_EDDY,0,EDO,IM*JM,EDO,1)
      call closeunit(iu_EDDY)

C**** DEFINE THE VERTICAL LAYERING EVERYWHERE EXCEPT LAYER 1 THICKNESS
      DZ(2)=10.
      DZO(1)=0.5*DZ(2)          ! 10./SQRT(FAC)
      BYDZO(1)=1./DZO(1)
      DO L=2,LMOM-1
        DZ(L+1)=DZ(L)*FAC
        DZO(L)=0.5*(DZ(L+1)+DZ(L)) !DZO(L-1)*FAC
        BYDZO(L)=1./DZO(L)
      END DO

C**** read in initial conditions and climatology for the temperatures at
C**** base of mixed layer, initialise deep arrays for start of run
      if (iniOCEAN) then
        stg3=0. ; dtg3=0. ; rtgo=0.

        call openunit("TG3M",iu_tg3m,.true.,.true.)
        READ(iu_tg3m) TITLE,TG3M
        WRITE(6,*) "Read from TG3M",TITLE
        call closeunit (iu_tg3m)

      end if
      return
C****
      end subroutine init_odeep

      SUBROUTINE io_ocean(kunit,iaction,ioerr)
!@sum  io_ocean reads and writes ocean arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,irsficno,lhead
      USE STATIC_OCEAN
      USE ODEEP_COM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "OCNDEEP01"

      WRITE(MODULE_HEADER(lhead+1:80),'(a45,i2,a)') 'R8 To(3,ijm)'//
     *     ', dim(ijm): MixLD,Stg3,dtg3,rtgo(',lmom,',ijm),'//
     *     'tg3m(12,ijm)'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,TOCEAN,Z1O,STG3,DTG3,RTGO,
     *     TG3M
      CASE (IOREAD:)            ! input from restart file
        SELECT CASE (IACTION)
        CASE (IRSFICNO)         ! no ocean initial conditions
C**** Note this is for reading in a full rsf file, but from a qflux run
C**** We do not check HEADER here because it will be wrong. The other
C**** data MUST be initialised by setting iniOCEAN=.TRUE. in init_ODEEP.
          READ (kunit) HEADER,TOCEAN,Z1O
        CASE DEFAULT            ! restart file
          READ (kunit,err=10) HEADER,TOCEAN,Z1O,STG3,DTG3,RTGO,TG3M
          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
            GO TO 10
          END IF
        END SELECT
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
C****
      END SUBROUTINE io_ocean

      SUBROUTINE io_ocdiag(kunit,it,iaction,ioerr)
!@sum  io_ocdiag reads and writes ocean diagnostic arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,irsfic,irerun,iowrite_single
     *     ,ioread_single,lhead,im,jm
      USE ODEEP_COM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "OCDIAGDEEP01"
!@var it input/ouput value of hour
      INTEGER, INTENT(INOUT) :: it
!@var RTGO4 dummy variable for reading in single precision
      REAL*4 RTGO4(LMOM,IM,JM)
      INTEGER I,J

C**** no output required for rsf files. Only acc files
      write(MODULE_HEADER(lhead+1:80),'(a,i2,a)')
     *     'R4 RTGO(',lmom,'im,jm)'

      SELECT CASE (IACTION)
      CASE (IOWRITE_SINGLE)     ! output to acc file
        WRITE (kunit,err=10) MODULE_HEADER,REAL(RTGO,KIND=4)
      CASE (IOREAD:)            ! input from restart file
        SELECT CASE (IACTION)
        CASE (ioread_single)    ! read in from acc file
          READ (kunit,err=10) HEADER,RTGO4
C**** sum RTGO over input files
          RTGO=RTGO+RTGO4
          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",HEADER
     *           ,MODULE_HEADER
            GO TO 10
          END IF
        END SELECT
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
C****
      END SUBROUTINE io_ocdiag

      SUBROUTINE reset_odiag(isum)
!@sum reset_odiag zeros out ocean diagnostics if needed
!@auth G. Schmidt
!@ver  1.0
      USE ODEEP_COM, only : rtgo
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: isum

C**** there is a confusion of definition for rtgo
C****   i) for the model run, it is the actual temperature anomaly
C****  ii) for post-processing, it is the accumulatated anomaly
C**** Thus it is only initiallised here for case ii).
      if (isum.eq.1) rtgo=0

      return
      end subroutine reset_odiag

      SUBROUTINE conserv_OCE(OCEANE)
!@sum  conserv_OCE calculates zonal ocean energy for Qflux ocean
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : shw,rhows
      USE MODEL_COM, only : im,jm,fim,focean
      USE GEOM, only : imaxj
      USE STATIC_OCEAN, only : tocean,z1o,z12o
      USE ODEEP_COM, only : dz,rtgo,lmom
      IMPLICIT NONE
!@var OCEANE zonal ocean energy (J/M^2)
      REAL*8, DIMENSION(JM) :: OCEANE
      INTEGER I,J,L

      OCEANE=0
      DO J=1,JM
        DO I=1,IMAXJ(J)
          IF (FOCEAN(I,J).gt.0) THEN
            OCEANE(J)=OCEANE(J)+(TOCEAN(1,I,J)*Z1O(I,J)
     *           +TOCEAN(2,I,J)*(Z12O(I,J)-Z1O(I,J)))*SHW*RHOWS
            DO L=2,LMOM
              OCEANE(J)=OCEANE(J)+(RTGO(L,I,J)*DZ(L)*SHW*RHOWS)
            END DO
          END IF
        END DO
      END DO
      OCEANE(1) =FIM*OCEANE(1)
      OCEANE(JM)=FIM*OCEANE(JM)
C****
      END SUBROUTINE conserv_OCE

      SUBROUTINE ODIFS
!@sum  ODFIS calculates heat diffusion at the base of the mixed layer
!@+    compares that to the control run's temperature, calls odffus,
!@+    and reduces the upper ocean temperatures by the amount of heat
!@+    that is diffused into the thermocline
!@auth Gary Russell/G. Schmidt
!@ver  1.0
!@calls ODFFUS
      USE FILEMANAGER
      USE CONSTANT, only : sday
      USE MODEL_COM, only : im,jm,focean,jmon,jday,jdate,itocean
     *     ,itoice
      USE GEOM, only : imaxj
      USE ODEEP_COM, only : tg3m,rtgo,stg3,dtg3,edo,dz,dzo,bydzo,lmom
      USE SEAICE_COM, only : rsi
      USE DAGCOM, only : aj,j_ftherm
      USE FLUXES, only : gtemp
      USE STATIC_OCEAN, only : z12o,tocean
      IMPLICIT NONE
      REAL*8, PARAMETER :: PERDAY=1./365d0
!@param ALPHA degree of implicitness (1 fully implicit,0 fully explicit)
      REAL*8, PARAMETER :: ALPHA=.5d0
      REAL*8 :: ADTG3
      INTEGER I,J,L
C****
C**** ACCUMULATE OCEAN TEMPERATURE AT MAXIMUM MIXED LAYER
C****
      DO J=1,JM
        DO I=1,IMAXJ(J)
          STG3(I,J)=STG3(I,J)+TOCEAN(3,I,J)
        END DO
      END DO
C****
C**** AT THE END OF EACH MONTH, UPDATE THE OCEAN TEMPERATURE
C**** DIFFERENCE AND REPLACE THE MONTHLY SUMMED TEMPERATURE
C****
      IF(JDATE.EQ.1) THEN
      DO J=1,JM
        DO I=1,IMAXJ(J)
          DTG3(I,J)=DTG3(I,J)+(STG3(I,J)-TG3M(I,J,JMON))
          TG3M(I,J,JMON)=STG3(I,J)
          STG3(I,J)=0.
        END DO
      END DO
      END IF
C****
C**** DIFFUSE THE OCEAN TEMPERATURE DIFFERENCE OF THE UPPER LAYERS
C**** INTO THE THERMOCLINE AND REDUCE THE UPPER TEMPERATURES BY THE
C**** HEAT THAT IS DIFFUSED DOWNWARD
C****
      DO J=1,JM
        DO I=1,IMAXJ(J)
          IF(FOCEAN(I,J).GT.0.) THEN

            ADTG3=DTG3(I,J)*PERDAY
            RTGO(1,I,J)=ADTG3
C**** Set first layer thickness
            DZ(1)=Z12O(I,J)

            CALL ODFFUS (SDAY,ALPHA,EDO(I,J),DZ,BYDZO,RTGO(1,I,J),LMOM)

            DO L=1,3
              TOCEAN(L,I,J)=TOCEAN(L,I,J)+(RTGO(1,I,J)-ADTG3)
            END DO
            AJ(J,J_FTHERM,ITOCEAN)=AJ(J,J_FTHERM,ITOCEAN)-(RTGO(1,I,J)
     *           -ADTG3)*Z12O(I,J)*FOCEAN(I,J)*(1.-RSI(I,J))
            AJ(J,J_FTHERM,ITOICE )=AJ(J,J_FTHERM,ITOICE )-(RTGO(1,I,J)
     *           -ADTG3)*Z12O(I,J)*FOCEAN(I,J)*RSI(I,J)
            GTEMP(1:2,1,I,J) = TOCEAN(1:2,I,J)
          END IF
        END DO
      END DO

      RETURN
      END SUBROUTINE ODIFS

      SUBROUTINE ODFFUS (DT,ALPHA,ED,DZ,BYDZO,R,LMIJ)
!@sum  ODFFUS calculates the vertical mixing of a tracer
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0
!@calls TRIDIAG
      IMPLICIT NONE
!@var LMIJ IS THE NUMBER OF VERTICAL LAYERS
      INTEGER, INTENT(IN) :: LMIJ
!@var ED diffusion coefficient between adjacent layers (m**2/s)
!@var ALPHA determines the time scheme (0 explicit,1 fully implicit)
!@var DT time step (s)
      REAL*8, INTENT(IN) :: ED,ALPHA,DT
!@var DZ the depth of the layers (m)
!@var BYDZO is the inverse of depth between layer centers (1/m)
      REAL*8, INTENT(IN) :: DZ(LMIJ),BYDZO(LMIJ-1)
!@var R tracer concentration
      REAL*8, INTENT(INOUT) :: R(LMIJ)

      REAL*8 AM(LMIJ),BM(LMIJ),CM(LMIJ),DM(LMIJ)
      INTEGER L
C**** SET UP TRIDIAGONAL MATRIX ENTRIES AND RIGHT HAND SIDE
      AM(1)=0
      BM(1)=DZ(1)+ALPHA*DT*ED*BYDZO(1)
      CM(1)=     -ALPHA*DT*ED*BYDZO(1)
      DM(1)=DZ(1)*R(1)-(1.-ALPHA)*DT*ED*(R(1)-R(2))*BYDZO(1)

      DO L=2,LMIJ-1
        AM(L)=     -ALPHA*DT* ED*BYDZO(L-1)
        BM(L)=DZ(L)+ALPHA*DT*(ED*BYDZO(L-1)+ED*BYDZO(L))
        CM(L)=     -ALPHA*DT*               ED*BYDZO(L)
        DM(L)=DZ(L)*R(L)+(1.-ALPHA)*DT*(ED*(R(L-1)-R(L))*BYDZO(L-1)
     *                                 -ED*(R(L)-R(L+1))*BYDZO(L))
      END DO

      AM(LMIJ)=        -ALPHA*DT*ED*BYDZO(LMIJ-1)
      BM(LMIJ)=DZ(LMIJ)+ALPHA*DT*ED*BYDZO(LMIJ-1)
      CM(LMIJ)=0.
      DM(LMIJ)=DZ(LMIJ)*R(LMIJ)+(1.-ALPHA)*DT*ED*
     *         (R(LMIJ-1)-R(LMIJ))*BYDZO(LMIJ-1)

      CALL TRIDIAG(AM,BM,CM,DM,R,LMIJ)

      RETURN
      END SUBROUTINE ODFFUS

      SUBROUTINE CHECKO(SUBR)
!@sum  CHECKO Checks whether deep ocean values are reasonable
!@auth Original Development Team
!@ver  1.0
      USE MODEL_COM, only : im,jm
      USE ODEEP_COM, only : lmom,stg3,dtg3,tg3m,rtgo
      USE STATIC_OCEAN, only : tocean
      IMPLICIT NONE

!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR
      LOGICAL QCHECKO
      INTEGER I,J

C**** Check for NaN/INF in ocean data
      CALL CHECK3(TOCEAN,3 ,IM,JM,SUBR,'toc')
      CALL CHECK3(DTG3  ,IM,JM, 1,SUBR,'dtg3')
      CALL CHECK3(TG3M  ,12,IM,JM,SUBR,'tg3m')
      CALL CHECK3(STG3  ,IM,JM,1 ,SUBR,'stg3')
      CALL CHECK3(RTGO,LMOM,IM,JM,SUBR,'rtgo')

      QCHECKO = .FALSE.
C**** Check for reasonable values for ocean variables
      DO J=1,JM
        DO I=1,IM
          IF (TOCEAN(1,I,J).lt.-2. .or. TOCEAN(1,I,J).gt.50.) THEN
            WRITE(6,*) 'After ',SUBR,': I,J,TOCEAN=',I,J,TOCEAN(1:3,I,J)
            QCHECKO = .TRUE.
          END IF
       END DO
      END DO
      IF (QCHECKO)
     *     call stop_model("CHECKO: Ocean variables out of bounds",255)

      END SUBROUTINE CHECKO

      SUBROUTINE diag_OCEAN
!@sum  diag_OCEAN prints out diagnostics for ocean
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : jm,lrunid,xlabel,idacc
      USE GEOM, only : imaxj,lat_dg
      USE ODEEP_COM, only : lmom,rtgo,dz
      USE DAGCOM, only : acc_period,qdiag,zoc
      IMPLICIT NONE
      CHARACTER LNAME*50,SNAME*30,UNITS*50
      INTEGER I,J,L
      REAL*8 ATGO(JM,LMOM),SCALED,ONES(JM)

C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
      IF(QDIAG) call open_jl(trim(acc_period)//'.o'//XLABEL(1:LRUNID),
     *  jm,lmom,0,lat_dg)
      LNAME="Zonally averaged deep ocean temperature anomaly"
      SNAME="tgo_deep_anom"
      UNITS="DEGREES C"
C**** calculate zonal average
      DO L=1,LMOM
        DO J=1,JM
          ATGO(J,L)=0.
          DO I=1,IMAXJ(J)
            ATGO(J,L)=ATGO(J,L)+RTGO(L,I,J)
          END DO
        END DO
      END DO
C**** depths are calculated from base of the mixed layer
      ZOC(1)=0.
      DO L=2,LMOM
        ZOC(L)=ZOC(L-1)+DZ(L)
      END DO
      SCALED=1./IDACC(12)
      ONES(1:JM)=1.
C**** Print out a depth/latitude plot of the deep ocean temp anomaly
      CALL JLMAP(LNAME,SNAME,UNITS,1,ZOC,ATGO,SCALED,ONES,ONES,LMOM,2,1)
C****
      if(qdiag) call close_jl

      RETURN
      END SUBROUTINE diag_OCEAN

