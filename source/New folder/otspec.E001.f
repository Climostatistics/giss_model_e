      program otspec
!@sum  otspec reads the ocean vertical fluxes saved from a GCM
!@+   run, calculates the ocean energy from the climatological
!@+   ocean file, and accumulates the zeroeth and first harmonics
!@+   of each of those quantities.  The spectral coefficients are
!@+   used to calculate the zeroeth and first harmonics of the
!@+   ocean transports.
C****
C**** Input: RSFIC - restart file of a run with prescribed sst
C****        VFLX - flux data from fixed sst run
C****        OSST - climatological ocean data
C****        MLMAX - ann max mixed layer depths
C****        TOPO - topography
C****        XCORR - XCORR produced by program (from vertflux.f)
C****        OCNOUT - ocean IC for start of new run (from dec31.f)
C**** Output:
C****        OHT - OTSPEC.RB150.M250D
C****      OHTLP - line plottable zonal northward heat transports
C****     RSFNEW - augmented restart file
C****
      USE PARAM
      USE CONSTANT, only : twopi,sday,rhows
      USE MODEL_COM, only: im,jm,lm,iowrite_mon,irerun
      USE TIMINGS, only : ntimeacc,timing,timestr
      USE STATIC_OCEAN
      USE DAGCOM, only : oa,koa
      USE SEAICE_COM, only : rsi,snowi,msi,ssi
      USE SEAICE, only : ace1i,ac2oim
      USE FLUXES, only : sss
      USE GEOM
      USE FILEMANAGER
      implicit none
      integer first_month, first_year, last_month,
     *        last_year, years, year, i, j, k,m,iok
      integer months, monthe, itime1, month, kday, last_day
      REAL*8 COT(IM,JM),AOT(IM,JM,4),BOT(IM,JM,4)
      REAL*4 AMPOT(IM,JM),PHAOT(IM,JM),COTS(IM,JM),TAU4
      REAL*8 ARG,COSDAY(4),SINDAY(4),VFX, XCORR, SYEARS, OMEG, OE,
     *  CV(IM,JM),AV(IM,JM,4),BV(IM,JM,4),AE(IM,JM,4),BE(IM,JM,4)
      CHARACTER*80 TITLE(5),TITLE0, RunID, file_name
      REAL*4 month_day(12)
      INTEGER iu_AIC,ioerr,iu_TOPO,iu_MLMAX,iu_VFLX,iu_OHTLP,iu_OCNOUT
     *     ,iu_XCORR,iu_OHT
      INTEGER ItimeX
      REAL*8 onht(jm),toceansv(3,IM,JM)
      real*8 z1ox(im,jm),z12o_max
C****
      character*4 month_name(12), tmonth, tyear
      data month_name /'JAN ','FEB ','MAR ','APR ','MAY ','JUN ',
     *                 'JUL ','AUG ','SEP ','OCT ','NOV ','DEC '/
      data month_day /31,28,31,30,31,30,31,31,30,31,30,31/

      DATA TITLE/
     1  'Amplitude of Ocean Transport (W/m**2) - 1st Harmonic',
     2  'Phase of Ocean Transport (days) - 1st Harmonic',
     3  'Mean Ocean Transport (W/m**2)',
     4  'Cosine coefficient of Ocean Transport (W/m**2)',
     5  'Sine coefficient of Ocean Transport (W/m**2)'/
C****
C**** Ocean vertical flux data saved from run
C****    K= 1  SNOWOI (instantaneous at noon GMT)
C****       2  FWSIM  (instantaneous at noon GMT)
C****       3  HSIT   (instantaneous at noon GMT)
C****       4  ENRGP  (integrated over the day)
C****       5  SRHDT  (integrated over the day)
C****       6  TRHDT  (for ocean, integrated over the day)
C****       7  SHDT   (for ocean, integrated over the day)
C****       8  EVHDT  (for ocean, integrated over the day)
C****       9  TRHDT  (for ocean ice, integrated over the day)
C****      10  SHDT   (for ocean ice, integrated over the day)
C****      11  EVHDT  (for ocean ice, integrated over the day)
C****      12  SRHDT  (for ocean ice, integrated over the day)
C****
C**** Extra array needed for dealing with advected ice
C****      13  HCHSI  (HORIZ CONV SEA ICE ENRG, INTEGRATED OVER THE DAY)
C****
      call getarg(1,RunID )
      call getarg(2,title0)
      read(title0,*) first_month
      call getarg(3,title0)
      read(title0,*) first_year
      call getarg(4,title0)
      read(title0,*) last_month
      call getarg(5,title0)
      read(title0,*) last_year
      years = last_year-first_year+1
      call getarg(6,title0)
      read(title0,*) z12o_max
C*
      call openunit("XCORR",iu_XCORR,.true.,.true.)
      READ(iu_XCORR) XCORR
      WRITE(6,*) ' XCORR=',XCORR
      call closeunit(iu_XCORR)
C****
C**** Read in input files
C****
      call openunit("SICE",iu_SICE,.true.,.true.)
      CALL READT (iu_SICE,0,DM,IM*JM,DM,1)
C*
C**** Read in FOCEAN - ocean fraction
C*
      call openunit("TOPO",iu_TOPO,.true.,.true.)
      CALL READT (iu_TOPO,0,FOCEAN,IM*JM,FOCEAN,1) ! Ocean fraction
      call closeunit(iu_TOPO)
C*
C**** Read in ocean data below mixed layer on December 31
C*
      call openunit("OCNOUT",iu_OCNOUT,.true.,.true.)
      READ (iu_OCNOUT) TOCEAN,Z1O
      WRITE(6,*)'TOCEAN,Z1O read from OCNOUT',TOCEAN(:,1,45),Z1O(1,45)
     *     ,TOCEAN(:,1,23),Z1O(1,23)
      call closeunit(iu_OCNOUT)
C*
C**** Zero out spectral coefficients
C*
      CV = 0. ; AV = 0. ; BV = 0. ; AE = 0. ; BE = 0.
      FLAND = 1. - FOCEAN
C****
C**** Calculate spherical geometry
C****
      call GEOM_B
C**** set up unit numbers for ocean climatologies
      call openunit("OSST",iu_OSST,.true.,.true.)
C**** Set up unit number of mixed layer depth climatogies
      call openunit("OCNML",iu_OCNML,.true.,.true.)
C**** find and limit ocean ann max mix layer depths
      z12o = 0.
      do m=1,jmpery
        CALL READT (iu_OCNML,0,z1ox,IM*JM,z1ox,1)
        do j=1,jm
        do i=1,im
ccc       z12o(i,j)=min( z12o_max , max(z12o(i,j),z1ox(i,j)) )
ccc   the above line could substitute for next 3 lines: same results ?
          if (focean(i,j).gt.0. .and. z1ox(i,j).gt.z12o_max)
     *        z1ox(i,j)=z12o_max
          if (z1ox(i,j).gt.z12o(i,j)) z12o(i,j)=z1ox(i,j)
        end do
        end do
      end do
      rewind iu_OCNML
      write(6,*) 'Mixed Layer Depths limited to',z12o_max
C**** define sea surface salinity (needed for OCLIM)
      sss(:,:)=sss0
C**** Initialise sea ice mass etc.
      do j=1,jm
        do i=1,im
          msi(i,j)=ac2oim
          ssi(1:2,i,j)=ssi0*ace1i*xsi(1:2)
          ssi(3:4,i,j)=ssi0*msi(i,j)*xsi(3:4)
        end do
      end do
C****
C**** Loop over years of data and days in the year
C****
      DO year=first_year,last_year
        months = 1
        monthe = 12
        if (year .eq. first_year) months = first_month
        if (year .eq. last_year)  monthe = last_month
        jday = JDendOfM(first_month-1)
        do month = months, monthe
          tmonth = month_name(month)
          write (tyear, '(i4.4)') year
          file_name = 'VFLXO'//trim(tmonth)//tyear
          call openunit(file_name,iu_VFLX,.true.,.true.)
          last_day = month_day(month)
          do kday = 1,last_day
            jday = jday + 1
            ARG  = JDAY*TWOPI/365.
            DO K=1,4
              COSDAY(K) = DCOS(K*ARG)
              SINDAY(K) = DSIN(K*ARG)
            END DO
C*
            call READi (iu_VFLX, itime,OA,itime1,IM*JM*KOA*2,iok)
            if(iok.gt.0) go to 555
C**** Interpolate daily ocean data from the monthly climatology
C*
            kocean = 0
            jmon = month
            jdate = kday
C*
            DO J = 1,JM
              DO I = 1,IM
                SNOWI(I,J) = OA(I,J,1)
              END DO
            END DO

            CALL OCLIM (.true.)

C***  Read in the ocean mixed layer depth data
C***  and interpolate them for the current day
            kocean = 1
            CALL OCLIM (.true.)

            CALL OSTRUC(.false.)
C****
C**** Calculate the vertical flux (J/m**2) and ocean energy (J/m**2)
C****
            DO J=1,JM
              DO I=1,IMAXJ(J)
                IF(FOCEAN(I,J).LE.0.)  CYCLE
                VFX = OA(I,J,4)
     *               + (1.-RSI(I,J))*(OA(I,J,6)+OA(I,J,7)+OA(I,J,8)
     *               +XCORR*OA(I,J,5))+ RSI(I,J)*(OA(I,J,9)+OA(I,J,10)
     *               +OA(I,J,11)+XCORR*OA(I,J,12))
     *               +OA(I,J,13)
C*
                OE = RSI(I,J)*OA(I,J,3)
     *               + ((Z1O(I,J)*RHOWS-OA(I,J,2))*TOCEAN(1,I,J)+
     *               (Z12O(I,J)-Z1O(I,J))*RHOWS*TOCEAN(2,I,J))*SHW
C*
C**** Accumulate the spectral coefficients
C*
                CV(I,J) = CV(I,J) + VFX
                DO K=1,4
                  AV(I,J,K) = AV(I,J,K) + VFX*COSDAY(K)
                  BV(I,J,K) = BV(I,J,K) + VFX*SINDAY(K)
                  AE(I,J,K) = AE(I,J,K) + OE *COSDAY(K)
                  BE(I,J,K) = BE(I,J,K) + OE *SINDAY(K)
                END DO
              END DO
            END DO
          end do
          call closeunit(iu_VFLX)
        end do
      END DO
      SYEARS = SDAY*365.*years
C****
C**** SCALE AV TO W/M**2 , AE TO J/M**2 TO CALCULATE SPECTRAL COEFF
C****
      DO J=1,JM
      DO I=1,IM
        CV(I,J) = CV(I,J)/SYEARS
        DO K=1,4
          AV(I,J,K) = AV(I,J,K)*2./SYEARS
          BV(I,J,K) = BV(I,J,K)*2./SYEARS
          AE(I,J,K) = AE(I,J,K)*2./(365.*years)
          BE(I,J,K) = BE(I,J,K)*2./(365.*years)
        END DO
      END DO
      END DO
C****
C**** Calculate the ocean transports spectral coefficients
C****
      OMEG = TWOPI/(SDAY*365.)
      DO J=1,JM
      DO I=1,IMAXJ(J)
        IF(FOCEAN(I,J).LE.0.)  CYCLE
        COT(I,J) =  -CV(I,J)
        DO K=1,4
          AOT(I,J,K) =  BE(I,J,K)*K*OMEG - AV(I,J,K)
          BOT(I,J,K) = -AE(I,J,K)*K*OMEG - BV(I,J,K)
        END DO
      END DO
      END DO
C**** Compute phase and amplitude of ocean transports
      DO J=1,JM
      DO I=1,IMAXJ(J)
        IF(FOCEAN(I,J).LE.0.) CYCLE
        AMPOT(I,J) = SQRT(AOT(I,J,1)*AOT(I,J,1)+BOT(I,J,1)*BOT(I,J,1))
        PHAOT(I,J) = ATAN2(BOT(I,J,1),AOT(I,J,1))*365./TWOPI
        COTS(I,J)  = COT(I,J)
      END DO
      END DO
      TAU4 = itime
      CALL MAP1 (IM,JM,ITIME,TITLE(1),AMPOT,SNGL(FOCEAN),1.,0.,0)
      CALL MAP1 (IM,JM,ITIME,TITLE(2),PHAOT,SNGL(FOCEAN),1.,0.,0)
      CALL MAP1 (IM,JM,ITIME,TITLE(3),COTS ,SNGL(FOCEAN),1.,0.,0)
C****
C**** Write ocean transports spectral coefficients to disk
C****
      call openunit("OHT",iu_OHT,.true.,.false.)
      WRITE (iu_OHT) BOT,AOT,COT,z12o_max
      print*,"OHT comp:",BOT(71,23,1:4),AOT(71,23,1:4),COT(71,23)
      call closeunit(iu_OHT)

C**** save tocean data
      TOCEANSV=TOCEAN

C**** Combine final restart file of PRESCRIBED OCEAN DATA model run
C**** with mean & bottom temperature of 2nd ocean layer to create
C**** Initial Conditions for a PREDICTED OCEAN DATA model run.
C****
      call openunit("RSFIC",iu_AIC,.true.,.true.)
      call io_rsf(iu_AIC,ItimeX,irerun,ioerr)
      call closeunit(iu_AIC)
C*
C****  Set the ocean temperature below the mixed layer
C*
      tocean(2:3,:,:) = toceansv(2:3,:,:)
C*
      ntimeacc = 1
      timing = 0
      timestr = " "
C*
      call openunit("RSFNEW",iu_AIC,.true.,.false.)
      call io_rsf(iu_AIC,ItimeX,iowrite_mon,ioerr)
      call closeunit(iu_AIC)

C**** Output aplot format file of ocean heat transports
      print*,"Calculating global northward heat transport..."

      onht(jm)=cot(1,jm)*im*dxyp(jm)*focean(1,jm)
      write(*,*) cot(1,jm),1.d-15*onht(jm),dxyp(jm)
      do j=jm-1,2,-1
        onht(j)=onht(j+1)
        do i=1,im
          onht(j)=onht(j)+COT(i,j)*dxyp(j)*focean(i,j)
        end do
      end do

      call openunit("OHTLP",iu_OHTLP,.false.,.false.)

      write(iu_OHTLP,*) 'Global Northward Ocean Heat Transport '
      write(iu_OHTLP,*) 'Latitude'
      write(iu_OHTLP,*) '10**15 W'
      write(iu_OHTLP,*) ' lat  ',RunID
      do j=2,jm
        write(iu_OHTLP,*) lat_dg(j,2),1d-15*onht(j)
      end do
      write(iu_OHTLP,*) ' '

      call closeunit(iu_OHTLP)
      WRITE(0,*) ' NORMAL END'

      STOP
 555  write (*,*) ' Reached end of file ',file_name
      call exit_rc (11)
      END


