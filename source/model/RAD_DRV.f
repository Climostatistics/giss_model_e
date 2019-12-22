!@sum RAD_DRV contains drivers for the radiation related routines
!@ver  1.0
!@cont COSZ0, init_RAD, RADIA
C**** semi-random cloud overlap (computed opt.d+diagn)
C**** to be used with R99E or later radiation  routines.  carbon/2
C****
      SUBROUTINE COSZ0
!@sum  COSZ0 calculates Earth's zenith angle, weighted by time/sunlight
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : twopi
      USE MODEL_COM
      USE GEOM, only : dlat,dlon,lon,sinip,cosip
      USE RADNCB, only : cosd,sind,sinj,cosj
      IMPLICIT NONE
      SAVE
      REAL*8, DIMENSION(IM) :: LT1,LT2,SLT1,SLT2,S2LT1,S2LT2
      REAL*8, DIMENSION(IM,JM) :: COSZ,COSZA
      COMMON/WORK5/LT1,LT2,SLT1,SLT2,S2LT1,S2LT2
C**** ZERO1 HAS TO EQUAL THE CUT-OFF VALUE FOR COSZ USED IN SOLAR
C**** COSZS WORKS CORRECTLY ONLY IF ZERO1 >> 1.D-3
      REAL*8, PARAMETER :: ZERO1=1.D-2
      INTEGER I,J,L
      REAL*8 S2DAWN,S2DUSK,ECOSZ,ECOSQZ,CLT1,CLT2,ZERO2,CDUSK,DUSK,DAWN
     *     ,SDUSK,SDAWN,CJCD,SJSD,SR1,CR1,SR2,CR2,ROT1,ROT2,DROT
C****
      ENTRY COSZT (ROT1,ROT2,COSZ)
C****
C**** THIS ENTRY COMPUTES THE ZENITH ANGLE WEIGHTED BY DAYTIME
C**** HOURS FROM ROT1 TO ROT2, GREENWICH MEAN TIME IN RADIANS.  ROT1
C**** MUST BE BETWEEN 0 AND 2*PI.  ROT2 MUST BE BETWEEN ROT1 AND
C**** ROT1+2*PI.  I=1 MUST LIE ON THE INTERNATIONAL DATE LINE.
C****
      DROT=ROT2-ROT1
C**** COMPUTE THE SINES AND COSINES OF THE INITIAL AND FINAL GMT'S
      SR1=SIN(ROT1)
      CR1=COS(ROT1)
      SR2=SIN(ROT2)
      CR2=COS(ROT2)
C**** COMPUTE THE INITIAL AND FINAL LOCAL TIMES (MEASURED FROM NOON TO
C****   NOON) AND THEIR SINES AND COSINES
      DO I=1,IM
        LT1(I)=ROT1+LON(I)
        SLT1(I)=SR1*COSIP(I)+CR1*SINIP(I)
        LT2(I)=ROT2+LON(I)
        SLT2(I)=SR2*COSIP(I)+CR2*SINIP(I)
      END DO
C****
C**** CALCULATION FOR POLAR GRID BOXES
C****
      DO J=1,JM,JM-1
        SJSD=SINJ(J)*SIND
        CJCD=COSJ(J)*COSD
        IF (SJSD+CJCD.GT.ZERO1) THEN
          IF (SJSD-CJCD.LT.0.) THEN
C**** AVERAGE COSZ FROM DAWN TO DUSK NEAR THE POLES
            DUSK=ACOS(-SJSD/CJCD)
            SDUSK=SQRT(CJCD*CJCD-SJSD*SJSD)/CJCD
            DAWN=-DUSK
            SDAWN=-SDUSK
            COSZ(1,J)=(SJSD*(DUSK-DAWN)+CJCD*(SDUSK-SDAWN))/TWOPI
          ELSE
C**** CONSTANT DAYLIGHT NEAR THE POLES
            COSZ(1,J)=SJSD
          END IF
        ELSE
C**** CONSTANT NIGHTIME NEAR THE POLES
          COSZ(1,J)=0.
        END IF
      END DO
C****
C**** LOOP OVER NON-POLAR LATITUDES
C****
      DO 500 J=2,JM-1
      SJSD=SINJ(J)*SIND
      CJCD=COSJ(J)*COSD
      IF (SJSD+CJCD.GT.ZERO1) THEN
      IF (SJSD-CJCD.LT.0.) THEN
C**** COMPUTE DAWN AND DUSK (AT LOCAL TIME) AND THEIR SINES
      DUSK=ACOS(-SJSD/CJCD)
      SDUSK=SQRT(CJCD*CJCD-SJSD*SJSD)/CJCD
      DAWN=-DUSK
      SDAWN=-SDUSK
C**** NEITHER CONSTANT DAYTIME NOR CONSTANT NIGHTIME AT THIS LATITUDE,
C**** LOOP OVER LONGITUDES
      ZERO2=ZERO1/CJCD
      DO 400 I=1,IM
C**** FORCE DUSK TO LIE BETWEEN LT1 AND LT1+2*PI
      IF (DUSK.GT.LT1(I)+ZERO2) GO TO 220
      DUSK=DUSK+TWOPI
      DAWN=DAWN+TWOPI
  220 IF (DAWN.LT.LT2(I)-ZERO2) GO TO 240
C**** CONTINUOUS NIGHTIME FROM INITIAL TO FINAL TIME
      COSZ(I,J)=0.
      GO TO 400
  240 IF (DAWN.GE.LT1(I)) GO TO 300
      IF (DUSK.LT.LT2(I)) GO TO 260
C**** CONTINUOUS DAYLIGHT FROM INITIAL TIME TO FINAL TIME
      COSZ(I,J)=SJSD+CJCD*(SLT2(I)-SLT1(I))/DROT
      GO TO 400
  260 IF (DAWN+TWOPI.LT.LT2(I)-ZERO2) GO TO 280
C**** DAYLIGHT AT INITIAL TIME AND NIGHT AT FINAL TIME
      COSZ(I,J)=(SJSD*(DUSK-LT1(I))+CJCD*(SDUSK-SLT1(I)))/DROT
      GO TO 400
C**** DAYLIGHT AT INITIAL AND FINAL TIMES WITH NIGHTIME IN BETWEEN
  280 COSZ(I,J)=(SJSD*(LT2(I)-DAWN-TWOPI+DUSK-LT1(I))+CJCD*
     *  (SLT2(I)-SDAWN+SDUSK-SLT1(I)))/DROT
      GO TO 400
  300 IF (DUSK.LT.LT2(I)) GO TO 320
C**** NIGHT AT INITIAL TIME AND DAYLIGHT AT FINAL TIME
      COSZ(I,J)=(SJSD*(LT2(I)-DAWN)+CJCD*(SLT2(I)-SDAWN))/DROT
      GO TO 400
C**** NIGHTIME AT INITIAL AND FINAL TIMES WITH DAYLIGHT IN BETWEEN
  320 COSZ(I,J)=(SJSD*(DUSK-DAWN)+CJCD*(SDUSK-SDAWN))/DROT
  400 CONTINUE
      ELSE
C**** CONSTANT DAYLIGHT AT THIS LATITUDE
        DO I=1,IM
          COSZ(I,J)=SJSD+CJCD*(SLT2(I)-SLT1(I))/DROT
        END DO
      END IF
      ELSE
C**** CONSTANT NIGHTIME AT THIS LATITUDE
        COSZ(1:IM,J)=0.
      END IF
  500 CONTINUE
      RETURN
C****
C****
      ENTRY COSZS (ROT1,ROT2,COSZ,COSZA)
C****
C**** THIS ENTRY COMPUTES THE ZENITH ANGLE TWICE, FIRST WEIGHTED BY THE
C**** DAYTIME HOURS FROM ROT1 TO ROT2 AND SECONDLY WEIGHTED BY THE
C**** INCIDENT SUN LIGHT FROM ROT1 TO ROT2.  COSZT MUST HAVE BEEN
C**** CALLED JUST PREVIOUSLY.
C****
      DROT=ROT2-ROT1
C**** COMPUTE THE SINES AND COSINES OF THE INITIAL AND FINAL GMT'S
      SR1=SIN(ROT1)
      CR1=COS(ROT1)
      SR2=SIN(ROT2)
      CR2=COS(ROT2)
C**** COMPUTE THE INITIAL AND FINAL LOCAL TIMES (MEASURED FROM NOON TO
C****   NOON) AND THEIR SINES AND COSINES
      DO I=1,IM
        LT1(I)=ROT1+LON(I)
        SLT1(I)=SR1*COSIP(I)+CR1*SINIP(I)
        CLT1=CR1*COSIP(I)-SR1*SINIP(I)
        S2LT1(I)=2.*SLT1(I)*CLT1
        LT2(I)=ROT2+LON(I)
        SLT2(I)=SR2*COSIP(I)+CR2*SINIP(I)
        CLT2=CR2*COSIP(I)-SR2*SINIP(I)
        S2LT2(I)=2.*SLT2(I)*CLT2
      END DO
C****
C**** CALCULATION FOR POLAR GRID BOXES
C****
      DO J=1,JM,JM-1
        SJSD=SINJ(J)*SIND
        CJCD=COSJ(J)*COSD
        IF (SJSD+CJCD.GT.ZERO1) THEN
          IF (SJSD-CJCD.LT.0.) THEN
C**** AVERAGE COSZ FROM DAWN TO DUSK NEAR THE POLES
            CDUSK=-SJSD/CJCD
            DUSK=ACOS(CDUSK)
            SDUSK=SQRT(CJCD*CJCD-SJSD*SJSD)/CJCD
            S2DUSK=2.*SDUSK*CDUSK
            DAWN=-DUSK
            SDAWN=-SDUSK
            S2DAWN=-S2DUSK
            ECOSZ=SJSD*(DUSK-DAWN)+CJCD*(SDUSK-SDAWN)
            ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SDUSK-SDAWN)+
     *           .5*CJCD*(DUSK-DAWN+.5*(S2DUSK-S2DAWN)))
            COSZ(1,J)=ECOSZ/TWOPI
            COSZA(1,J)=ECOSQZ/ECOSZ
          ELSE
C**** CONSTANT DAYLIGHT NEAR THE POLES
            ECOSZ=SJSD*TWOPI
            ECOSQZ=SJSD*ECOSZ+.5*CJCD*CJCD*TWOPI
            COSZ(1,J)=ECOSZ/TWOPI
            COSZA(1,J)=ECOSQZ/ECOSZ
          END IF
        ELSE
C**** CONSTANT NIGHTIME NEAR THE POLES
          COSZ(1,J)=0.
          COSZA(1,J)=0.
        END IF
      END DO
C****
C**** LOOP OVER NON-POLAR LATITUDES
C****
      DO 900 J=2,JM-1
      SJSD=SINJ(J)*SIND
      CJCD=COSJ(J)*COSD
      IF (SJSD+CJCD.GT.ZERO1) THEN
      IF (SJSD-CJCD.LT.0.) THEN
C**** COMPUTE DAWN AND DUSK (AT LOCAL TIME) AND THEIR SINES
      CDUSK=-SJSD/CJCD
      DUSK=ACOS(CDUSK)
      SDUSK=SQRT(CJCD*CJCD-SJSD*SJSD)/CJCD
      S2DUSK=2.*SDUSK*CDUSK
      DAWN=-DUSK
      SDAWN=-SDUSK
      S2DAWN=-S2DUSK
C**** NEITHER CONSTANT DAYTIME NOR CONSTANT NIGHTIME AT THIS LATITUDE,
C**** LOOP OVER LONGITUDES
      ZERO2=ZERO1/CJCD
      DO 800 I=1,IM
C**** FORCE DUSK TO LIE BETWEEN LT1 AND LT1+2*PI
      IF (DUSK.GT.LT1(I)+ZERO2) GO TO 620
      DUSK=DUSK+TWOPI
      DAWN=DAWN+TWOPI
  620 IF (DAWN.LT.LT2(I)-ZERO2) GO TO 640
C**** CONTINUOUS NIGHTIME FROM INITIAL TO FINAL TIME
      COSZ(I,J)=0.
      COSZA(I,J)=0.
      GO TO 800
  640 IF (DAWN.GE.LT1(I)) GO TO 700
      IF (DUSK.LT.LT2(I)) GO TO 660
C**** CONTINUOUS DAYLIGHT FROM INITIAL TIME TO FINAL TIME
      ECOSZ=SJSD*DROT+CJCD*(SLT2(I)-SLT1(I))
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SLT2(I)-SLT1(I))+
     *  .5*CJCD*(DROT+.5*(S2LT2(I)-S2LT1(I))))
      COSZ(I,J)=ECOSZ/DROT
      COSZA(I,J)=ECOSQZ/ECOSZ
      GO TO 800
  660 IF (DAWN+TWOPI.LT.LT2(I)-ZERO2) GO TO 680
C**** DAYLIGHT AT INITIAL TIME AND NIGHT AT FINAL TIME
      ECOSZ=SJSD*(DUSK-LT1(I))+CJCD*(SDUSK-SLT1(I))
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SDUSK-SLT1(I))+
     *  .5*CJCD*(DUSK-LT1(I)+.5*(S2DUSK-S2LT1(I))))
      COSZ(I,J)=ECOSZ/DROT
      COSZA(I,J)=ECOSQZ/ECOSZ
      GO TO 800
C**** DAYLIGHT AT INITIAL AND FINAL TIMES WITH NIGHTIME IN BETWEEN
  680 ECOSZ=SJSD*(DROT-DAWN-TWOPI+DUSK)+
     *  CJCD*(SLT2(I)-SDAWN+SDUSK-SLT1(I))
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SDUSK-SLT1(I)+SLT2(I)-SDAWN)+
     *  .5*CJCD*(DUSK+DROT-DAWN-TWOPI+
     *  .5*(S2DUSK-S2LT1(I)+S2LT2(I)-S2DAWN)))
      COSZ(I,J)=ECOSZ/DROT
      COSZA(I,J)=ECOSQZ/ECOSZ
      GO TO 800
  700 IF (DUSK.GE.LT2(I)) THEN
C**** NIGHT AT INITIAL TIME AND DAYLIGHT AT FINAL TIME
        ECOSZ=SJSD*(LT2(I)-DAWN)+CJCD*(SLT2(I)-SDAWN)
        ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SLT2(I)-SDAWN)+
     *       .5*CJCD*(LT2(I)-DAWN+.5*(S2LT2(I)-S2DAWN)))
        COSZ(I,J)=ECOSZ/DROT
        COSZA(I,J)=ECOSQZ/ECOSZ
      ELSE
C**** NIGHTIME AT INITIAL AND FINAL TIMES WITH DAYLIGHT IN BETWEEN
        ECOSZ=SJSD*(DUSK-DAWN)+CJCD*(SDUSK-SDAWN)
        ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SDUSK-SDAWN)+
     *       .5*CJCD*(DUSK-DAWN+.5*(S2DUSK-S2DAWN)))
        COSZ(I,J)=ECOSZ/DROT
        COSZA(I,J)=ECOSQZ/ECOSZ
      END IF
  800 CONTINUE
      ELSE
C**** CONSTANT DAYLIGHT AT THIS LATITUDE
        DO I=1,IM
          ECOSZ=SJSD*DROT+CJCD*(SLT2(I)-SLT1(I))
          ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SLT2(I)-SLT1(I))+
     *         .5*CJCD*(DROT+.5*(S2LT2(I)-S2LT1(I))))
          COSZ(I,J)=ECOSZ/DROT
          COSZA(I,J)=ECOSQZ/ECOSZ
        END DO
      END IF
C**** CONSTANT NIGHTIME AT THIS LATITUDE
      ELSE
        COSZ(1:IM,J)=0.
        COSZA(1:IM,J)=0.
      END IF
  900 CONTINUE
      RETURN
      END

      SUBROUTINE init_RAD
!@sum  init_RAD initialises radiation code
!@auth Original Development Team
!@ver  1.0
!@calls RE001:RCOMP1, ORBPAR
      USE FILEMANAGER
      USE PARAM
      USE CONSTANT, only : grav,bysha,twopi
      USE MODEL_COM, only : jm,lm,ls1,dsig,sige,psfmpt,ptop,dtsrc,nrad
     *     ,kradia
      USE GEOM, only : dlat,lat_dg
      USE RE001, only : rcomp1,writer,writet       ! routines
     &     ,FULGAS ,PTLISO ,KTREND ,NL ,NLP, PLB, PTOPTR
     *     ,KCLDEM,KVEGA6,MOZONE,KSOLAR, SHL, snoage_fac_max, KZSNOW
     *     ,KYEARS,KJDAYS,MADLUV, KYEARG,KJDAYG,MADGHG
     *     ,KYEARO,KJDAYO,MADO3M, KYEARA,KJDAYA,MADAER
     *     ,KYEARD,KJDAYD,MADDST, KYEARV,KJDAYV,MADVOL
     *     ,KYEARE,KJDAYE,MADEPS, KYEARR,KJDAYR
      USE RADNCB, only : s0x,co2x,ch4x,h2ostratx,s0_yr,s0_day
     *     ,ghg_yr,ghg_day,volc_yr,volc_day,aero_yr,O3_yr,crops_yr
     *     ,lm_req,coe,sinj,cosj,H2ObyCH4,dH2O
     *     ,obliq,eccn,omegt,obliq_def,eccn_def,omegt_def
     *     ,calc_orb_par,paleo_orb_yr
     *     ,PLB0,shl0  ! saved to avoid OMP-copyin of input arrays
      USE DAGCOM, only : iwrite,jwrite,itwrite
      IMPLICIT NONE

      INTEGER J,L,LR,MADVEL
      REAL*8 COEX,SPHIS,CPHIS,PHIN,SPHIN,CPHIN,PHIM,PHIS,PLBx(LM+1)
     *     ,pyear
!@var NRFUN indices of unit numbers for radiation routines
      INTEGER NRFUN(14),IU
!@var RUNSTR names of files for radiation routines
      CHARACTER*5 :: RUNSTR(14) = (/"RADN1","RADN2","RADN3",
     *     "RADN4","RADN5","RADN6","RADN7","RADN8",
     *     "RADN9","RADNA","RADNB","RADNC","RADND",
     *     "RADNE"/)
!@var QBIN true if files for radiation input files are binary
      LOGICAL :: QBIN(14)=(/.TRUE.,.TRUE.,.FALSE.,.FALSE.,.TRUE.,.TRUE.
     *     ,.TRUE.,.TRUE.,.FALSE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE./)

C**** sync radiation parameters from input
      call sync_param( "S0X", S0X )
      call sync_param( "CO2X", CO2X )
      call sync_param( "CH4X", CH4X )
      call sync_param( "H2OstratX", H2OstratX )
      call sync_param( "H2ObyCH4", H2ObyCH4 )
      call sync_param( "S0_yr", S0_yr )
      call sync_param( "ghg_yr", ghg_yr )
      call sync_param( "ghg_day", ghg_day )
      call sync_param( "S0_day", S0_day )
      call sync_param( "volc_yr", volc_yr )
      call sync_param( "volc_day", volc_day )
      call sync_param( "aero_yr", aero_yr )
      call sync_param( "O3_yr", O3_yr )
      call sync_param( "crops_yr", crops_yr )
      call sync_param( "MOZONE", MOZONE )
      call sync_param( "KSOLAR", KSOLAR )
      call sync_param( "KVEGA6", KVEGA6 )
      call sync_param( "KZSNOW", KZSNOW )
      call sync_param( "calc_orb_par", calc_orb_par )
      call sync_param( "paleo_orb_yr", paleo_orb_yr )
      call sync_param( "snoage_fac_max", snoage_fac_max )
      if(snoage_fac_max.lt.0. .or. snoage_fac_max.gt.1.) then
        write(*,*) 'set 0<snoage_fac_max<1, not',snoage_fac_max
        call stop_model('init_RAD: snoage_fac_max out of range',255)
      end if

C**** Set orbital parameters appropriately
      if (calc_orb_par.eq.1) then ! calculate from paleo-year
        pyear=1950.-paleo_orb_yr ! since 0 BP is defined as 1950CE
        call orbpar(pyear, eccn, obliq, omegt)
        write(6,*)
        write(6,*) " Orbital Parameters Calculated:"
        write(6,'(a,f8.0,a,f8.0,a)') "   Paleo-year: ",pyear," (CE);",
     *       paleo_orb_yr," (BP)"
        write(6,'(a,f8.7,a,f8.7,a)') "   Eccentricity: ",eccn,
     *       " (default = ",eccn_def,")"
        write(6,'(a,f9.6,a,f9.6,a)') "   Obliquity (degs): ",obliq,
     *       " (default = ",obliq_def,")"
        write(6,'(a,f7.3,a,f7.3,a)') "   Precession (degs from ve): ",
     *       omegt," (default = ",omegt_def,")"
        write(6,*)
      else  ! set from defaults (defined in CONSTANT module)
        omegt=omegt_def
        obliq=obliq_def
        eccn=eccn_def
      end if

C**** COMPUTE THE AREA WEIGHTED LATITUDES AND THEIR SINES AND COSINES
      PHIS=-.25*TWOPI
      SPHIS=-1.
      CPHIS=0.
      DO J=1,JM-1
        PHIN=DLAT*(J-.5*JM)
        SPHIN=SIN(PHIN)
        CPHIN=COS(PHIN)
        PHIM=(PHIN*SPHIN+CPHIN-PHIS*SPHIS-CPHIS)/(SPHIN-SPHIS)
        SINJ(J)=SIN(PHIM)
        COSJ(J)=COS(PHIM)
        PHIS=PHIN
        SPHIS=SPHIN
        CPHIS=CPHIN
      END DO
      PHIN=.25*TWOPI
      SPHIN=1.
      CPHIN=0.
      PHIM=(PHIN*SPHIN+CPHIN-PHIS*SPHIS-CPHIS)/(SPHIN-SPHIS)
      SINJ(JM)=SIN(PHIM)
      COSJ(JM)=COS(PHIM)
C****
C**** SET THE CONTROL PARAMETERS FOR THE RADIATION (need mean pressures)
C****
      NL=LM+LM_REQ
      NLP=NL+1
      COEX=1d-2*GRAV*BYSHA
      DO L=1,LM
        COE(L)=DTsrc*COEX/DSIG(L)
        PLB(L)=SIGE(L)*PSFMPT+PTOP
        PLBx(L)=PLB(L)           ! needed for CH4 prod. H2O
      END DO
      PLB(LM+1)=SIGE(LM+1)*PSFMPT+PTOP
      PLB(LM+2)=.5*PLB(LM+1)
      PLB(NL)=.2d0*PLB(LM+1)
      PLB(NL+1)=1d-5
      PTOPTR=PTOP ! top of sigma-coord.system
      DO LR=LM+1,NL
        COE(LR)=DTsrc*NRAD*COEX/(PLB(LR)-PLB(LR+1))
        PLB0(LR-LM) = PLB(LR+1)
      END DO
      if (kradia.gt.1) then
        do l=1,ls1-1
          COE(L)=DTsrc*nrad*COEX/DSIG(L)
        end do
        do l=ls1,lm
          COE(L)=DTsrc*NRAD*COEX/(PLB(L)-PLB(L+1))
        end do
      end if
      PTLISO=15.
      KTREND=1   !  GHgas trends are determined by input file
!note KTREND=0 is a possible but virtually obsolete option
C****
C             Model Add-on Data of Extended Climatology Enable Parameter
C     MADO3M  =  1   Reads            1951-1997 Ozone climatology RFILEA
C     MADAER  =  1   Reads   Aerosol 50y tropospheric climatology RFILE5
C     MADDST  =  1   Reads   Dust-windblown mineral climatology   RFILE6
C     MADVOL  =  1   Reads   Volcanic 1950-00 aerosol climatology RFILE7
C     MADEPS  =  1   Reads   Epsilon cloud heterogeniety data     RFILE8
C     MADLUV  =  1   Reads   Lean's SolarUV 1882-1998 variability RFILE9
C**** Radiative forcings are either constant = obs.value at given yr/day
C****    or time dependent (year=0); if day=0 an annual cycle is used
C****                                         even if the year is fixed
      KYEARS=s0_yr   ; KJDAYS=s0_day ;  MADLUV=1   ! solar 'constant'
      KYEARG=ghg_yr  ; KJDAYG=ghg_day              ! well-mixed GHGases
      if(ghg_yr.gt.0)  MADGHG=0                    ! skip GHG-updating
      KYEARO=O3_yr   ; KJDAYO=0 ;       MADO3M=1   ! ozone (ann.cycle)
      KYEARA=Aero_yr ; KJDAYA=0 ;       MADAER=1 !trop.aeros (ann.cycle)
      KYEARV=Volc_yr ; KJDAYV=Volc_day; MADVOL=1   ! Volc. Aerosols
C**** NO time history (yet), except for ann.cycle, for forcings below;
C****  if arg3=day0 (1->365), data from that day are used all year
      KYEARD=0       ; KJDAYD=0 ;       MADDST=1   ! Desert dust
      KYEARE=0       ; KJDAYE=0 ;       MADEPS=1   !cloud Epsln - KCLDEP
      KYEARR=0       ; KJDAYR=0           ! surf.reflectance (ann.cycle)
C**** New options (currently not used)
      KCLDEM=0  ! 0:old 1:new LW cloud scattering scheme
      if (ktrend.ne.0) then
C****   Read in time history of well-mixed greenhouse gases
        call openunit('GHG',iu,.false.,.true.)
        call ghghst(iu)
        call closeunit(iu)
        if (H2ObyCH4.gt.0..and.Kradia.le.0) then
C****     Read in dH2O: H2O prod.rate in kg/m^2 per day and ppm_CH4
          call openunit('dH2O',iu,.false.,.true.)
          call getqma(iu,lat_dg,plbx,dh2o,lm,jm)
          call closeunit(iu)
        end if
      end if
C**** set up unit numbers for 14 more radiation input files
      DO IU=1,14
        IF (IU.EQ.12.OR.IU.EQ.13) CYCLE ! not used in GCM
        call openunit(RUNSTR(IU),NRFUN(IU),QBIN(IU),.true.)
      END DO
C***********************************************************************
C     Main Radiative Initializations
C     ------------------------------------------------------------------
      CALL RCOMP1 (NRFUN)
      CALL WRITER(6,0)  ! print out ispare/fspare control parameters
C***********************************************************************
      DO IU=1,14
        IF (IU.EQ.12.OR.IU.EQ.13) CYCLE ! not used in GCM
        call closeunit(NRFUN(IU))
      END DO
C**** Save initial (currently permanent and global) Q in rad.layers
      do LR=1,LM_REQ
        shl0(LR) = shl(LM+LR)
      end do
      write(6,*) 'spec.hum in rad.equ.layers:',shl0
C**** Optionally scale selected greenhouse gases
      IF(ghg_yr.gt.0) FULGAS(2)=FULGAS(2)*CO2X
      if(ghg_yr.gt.0) FULGAS(7)=FULGAS(7)*CH4X
      IF(H2OstratX.GE.0.) FULGAS(1)=FULGAS(1)*H2OstratX
C**** write trend table for forcing 'itwrite' for years iwrite->jwrite
C**** itwrite: 1-2=GHG 3=So 4-5=O3 6-9=aerosols: Trop,DesDust,Volc,Total
      if(jwrite.gt.1500) call writet (6,itwrite,iwrite,jwrite,1,0)
C****
      END SUBROUTINE init_RAD

      SUBROUTINE RADIA
!@sum  RADIA adds the radiation heating to the temperatures
!@auth Original Development Team
!@ver  1.0
!@calls tropwmo, RE001:rcompt, RE001:rcompx, RE001:writer, coszs, coszt
      USE CONSTANT, only : sday,lhe,lhs,twopi,tf,stbo,rhow,mair,grav
     *     ,kapa
      USE MODEL_COM
      USE GEOM
      USE RE001
     &  , only : writer,rcompx,rcompt ! routines
     &          ,lx  ! for threadprivate copyin common block
     &          ,tauwc0,tauic0 ! set in radpar block data
C     INPUT DATA         ! not (i,j) dependent
     X          ,S00WM2,RATLS0,S0,JYEARR=>JYEAR,JDAYR=>JDAY
C     INPUT DATA  (i,j) dependent
     &             ,JLAT,ILON,nl,nlp, PLB ,TLB,TLM ,SHL, ltopcl
     &             ,TAUWC ,TAUIC ,SIZEWC ,SIZEIC
     &             ,POCEAN,PEARTH,POICE,PLICE,PLAKE,COSZ,PVT
     &             ,TGO,TGE,TGOI,TGLI,TSL,WMAG,WEARTH
     &             ,AGESN,SNOWE,SNOWOI,SNOWLI
     &             ,hsn,hin,hmp,fmp,flags,LS1_loc,snow_frac,lkdepth
C     OUTPUT DATA
     &          ,TRDFLB ,TRNFLB ,TRUFLB, TRFCRL
     &          ,SRDFLB ,SRNFLB ,SRUFLB, SRFHRL
     &          ,PLAVIS ,PLANIR ,ALBVIS ,ALBNIR ,FSRNFG
     &          ,SRRVIS ,SRRNIR ,SRAVIS ,SRANIR ,SRXVIS
     &          ,BTEMPW
      USE RADNCB, only : rqt,srhr,trhr,fsf,cosz1,s0x,rsdist,lm_req
     *     ,coe,plb0,shl0,tchg,alb,fsrdir
      USE RANDOM
      USE CLOUDS_COM, only : tauss,taumc,svlhx,rhsav,svlat,cldsav,
     *     cldmc,cldss,csizmc,csizss,llow,lmid,lhi
      USE PBLCOM, only : wsavg,tsavg
      USE DAGCOM, only : aj,areg,jreg,aij,ail,ajl,asjl,adiurn,
     *     iwrite,jwrite,itwrite,ndiupt,j_pcldss,j_pcldmc,ij_pmccld,
     *     j_clddep,j_pcld,ij_cldcv,ij_pcldl,ij_pcldm,ij_pcldh,
     *     ij_cldtppr,j_srincp0,j_srnfp0,j_srnfp1,j_srincg,
     *     j_srnfg,j_brtemp,j_trincg,j_hsurf,j_hatm,j_plavis,ij_trnfp0,
     *     ij_srnfp0,ij_srincp0,ij_srnfg,ij_srincg,ij_btmpw,ij_srref
     *     ,ij_srvis,j50n,j70n,j_clrtoa,j_clrtrp,j_tottrp,il_req,il_r50n
     *     ,il_r70n,ijdd,idd_cl7,idd_cl6,idd_cl5,idd_cl4,idd_cl3,idd_cl2
     *     ,idd_cl1,idd_ccv,idd_isw,idd_palb,idd_galb,idd_absa,j5s,j5n
     *     ,jl_srhr,jl_trcr,jl_totcld,jl_sscld,jl_mccld,ij_frmp
     *     ,jl_wcld,jl_icld,jl_wcod,jl_icod,jl_wcsiz,jl_icsiz
     *     ,ij_clr_srincg,ij_CLDTPT,ij_cldt1t,ij_cldt1p,ij_cldcv1
     *     ,ij_wtrcld,ij_icecld,ij_optdw,ij_optdi
     *     ,AFLX_ST
      USE DYNAMICS, only : pk,pedn,plij,pmid,pdsig,ltropo,am,byam
      USE SEAICE, only : rhos,ace1i,rhoi
      USE SEAICE_COM, only : rsi,snowi,pond_melt,msi,flag_dsws
      USE GHYCOM, only : snowe_com=>snowe,snoage,wearth_com=>wearth
     *     ,aiearth,fr_snow_rad_ij
      USE LANDICE_COM, only : snowli_com=>snowli
      USE LAKES_COM, only : flake,mwl
      USE FLUXES, only : gtemp
      IMPLICIT NONE
C
C     INPUT DATA   partly (i,j) dependent, partly global
      REAL*8 U0GAS,taulim
      COMMON/RADCOM_hybrid/U0GAS(LX,12)
!$OMP  THREADPRIVATE(/RADCOM_hybrid/)

      REAL*8, DIMENSION(IM,JM) :: COSZ2,COSZA,TRINCG,BTMPW,WSOIL,fmp_com
      REAL*8, DIMENSION(4,IM,JM) :: SNFS,TNFS
      REAL*8, DIMENSION(LM_REQ,IM,JM) :: TRHRS,SRHRS
      REAL*8, DIMENSION(0:LM+LM_REQ,IM,JM) :: TRHRA,SRHRA ! for adj.frc
      REAL*8, DIMENSION(LM) :: TOTCLD

      INTEGER, SAVE :: JDLAST = -9
      INTEGER I,J,L,K,KR,LR,JR,IH,INCH,JK,IT,iy,iend,icc1
      REAL*8 ROT1,ROT2,PLAND,PIJ,CSS,CMC,DEPTH,QSS,TAUSSL,RANDSS
     *     ,TAUMCL,ELHX,CLDCV,DXYPJ,SRNFLG,X,OPNSKY,CSZ2,tauup,taudn
     *     ,taucl,wtlin,MSTRAT,STRATQ,STRJ,MSTJ,optdw,optdi
     *     ,QR(LM,IM,JM),CLDinfo(LM,3,IM,JM)
      REAL*8 QSAT
      LOGICAL NO_CLOUD_ABOVE
C
      REAL*8  RDSS(LM,IM,JM),RDMC(IM,JM), AREGIJ(7,IM,JM)
      INTEGER ICKERR,JCKERR
C
C****
C**** FLAND     LAND COVERAGE (1)
C**** FLICE     LAND ICE COVERAGE (1)
C****
C**** GTEMP(1)  GROUND TEMPERATURE ARRAY OVER ALL SURFACE TYPES (C)
C****   RSI  RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)
C****
C**** VDATA  1-11 RATIOS FOR THE 11 VEGETATION TYPES (1)
C****

C**** limit optical cloud depth from below: taulim
      taulim=min(tauwc0,tauic0) ! currently both .001
      tauwc0 = taulim ; tauic0 = taulim
C**** Calculate mean cosine of zenith angle for the current physics step
      ROT1=(TWOPI*MOD(ITIME,NDAY))/NDAY  ! MOD(ITIME,NDAY)*TWOPI/NDAY ??
      ROT2=ROT1+TWOPI*DTsrc/SDAY
      CALL COSZT (ROT1,ROT2,COSZ1)

      if (kradia.gt.0) then    ! read in all rad. input data (frc.runs)
        iend = 1
        it = itime-1           ! make sure, at least 1 record is read
        do while (mod(itime-it,8760).ne.0)
          read(iu_rad,end=10,err=10) it,T,RQT,TsAvg,QR,P,CLDinfo
     *     ,rsi,msi,(((GTEMP(1,k,i,j),k=1,4),i=1,im),j=1,jm),wsoil,wsavg
     *     ,snowi,snowli_com,snowe_com,snoage,fmp_com,flag_dsws,ltropo
     *     ,srhra,trhra,iy  ! original output data (for adj.frc. only)
          if (qcheck) write(6,*) 'reading RADfile at Itime',Itime,it,iy
        end do
        iend = 0
   10   if (it.ne.iy.or.iend.eq.1) then
          write(6,*) 'RAD input file bad or too short:',itime,it,iy,iend
          call stop_model('RADIA: input file bad or too short',255)
        end if
C****   Find arrays derived from P : PEdn and PK (forcing experiments)
        do j=1,jm
        do i=1,imaxj(j)
          pedn(LM+1,i,j) = SIGE(LM+1)*PSFMPT+PTOP
        do l=1,lm
          pij=p(i,j)
          if(l.ge.ls1) pij=psfmpt
          pedn(l,i,j) = SIGE(L)*PIJ+PTOP
          pk(l,i,j)   = (SIG(L)*PIJ+PTOP)**KAPA
        end do ; end do ; end do
      end if

      IF (MODRD.NE.0) GO TO 900
      IDACC(2)=IDACC(2)+1
C****
C**** Interface with radiation routines, done only every NRAD time steps
C****
C**** Calculate mean cosine of zenith angle for the full radiation step
      ROT2=ROT1+TWOPI*NRAD*DTsrc/SDAY
      CALL COSZS (ROT1,ROT2,COSZ2,COSZA)
      JDAYR=JDAY
      JYEARR=JYEAR
C*********************************************************
C     Update time dependent radiative parameters each day
      IF(JDAY.NE.JDLAST) CALL RCOMPT
C*********************************************************
      JDLAST=JDAY
      S0=S0X*S00WM2*RATLS0/RSDIST

         if(kradia.le.0) then
      IF (QCHECK) THEN
C****   Calculate mean strat water conc
        STRATQ=0.
        MSTRAT=0.
        DO J=1,JM
          STRJ=0.
          MSTJ=0.
          DO I=1,IMAXJ(J)
            DO L=LTROPO(I,J)+1,LM
              STRJ=STRJ+Q(I,J,L)*AM(L,I,J)*DXYP(J)
              MSTJ=MSTJ+AM(L,I,J)*DXYP(J)
            END DO
          END DO
          IF (J.eq.1 .or. J.eq.JM) THEN
            STRJ=STRJ*FIM
            MSTJ=MSTJ*FIM
          END IF
          STRATQ=STRATQ+STRJ
          MSTRAT=MSTRAT+MSTJ
        END DO
        PRINT*,"Strat water vapour (ppmv), mass (mb)",1d6*STRATQ*mair
     *       /(18.*MSTRAT),PMTOP+1d-2*GRAV*MSTRAT/AREAG
      END IF
C
C**** GET THE RANDOM NUMBERS OUTSIDE PARALLEL REGIONS
C**** but keep MC calculation seperate from SS clouds
C**** MC clouds are considered as a block for each I,J grid point
      DO J=1,JM                    ! complete overlap
      DO I=1,IMAXJ(J)
        RDMC(I,J) = RANDU(X)
      END DO
      END DO
C**** SS clouds are considered as a block for each continuous cloud
      DO J=1,JM                    ! semi-random overlap
      DO I=1,IMAXJ(J)
        NO_CLOUD_ABOVE = .TRUE.
        DO L=LM,1,-1
          IF(TAUSS(L,I,J).le.taulim) CLDSS(L,I,J)=0.
          IF(TAUMC(L,I,J).le.taulim) CLDMC(L,I,J)=0.
          IF(CLDSS(L,I,J).GT.0.) THEN
            IF (NO_CLOUD_ABOVE) THEN
              RANDSS = RANDU(X)
              NO_CLOUD_ABOVE = .FALSE.
            END IF
          ELSE
            RANDSS = 1.
            NO_CLOUD_ABOVE = .TRUE.
          END IF
          RDSS(L,I,J) = RANDSS
        END DO
      END DO
      END DO
      end if                    ! kradia le 0
C****
C**** MAIN J LOOP
C****
      ICKERR=0
      JCKERR=0
!$OMP  PARALLEL PRIVATE(CSS,CMC,CLDCV, DEPTH,OPTDW,OPTDI, ELHX,
!$OMP*   I,INCH,IH,IT, J, K,KR, L,LR,icc1, OPNSKY, CSZ2,
!$OMP*   PLAND,PIJ, QSS, TOTCLD,TAUSSL,TAUMCL)
!$OMP*   COPYIN(/RADCOM_hybrid/)
!$OMP*   SHARED(ITWRITE)
!$OMP    DO SCHEDULE(DYNAMIC,2)
!$OMP*   REDUCTION(+:ICKERR,JCKERR)
      DO 600 J=1,JM
      NL=LM+LM_REQ ; NLP=NL+1   ! radiation allows var. # of layers
C**** Radiation input files use a 72x46 grid independent of IM and JM
C**** (ilon,jlat) is the 4x5 box containing the center of box (i,j)
      JLAT=INT(1.+(J-1.)*45./(JM-1.)+.5)  !  lat_index w.r.to 72x46 grid
C****
C**** MAIN I LOOP
C****
      DO I=1,IMAXJ(J)
      ILON=INT(.5+(I-.5)*72./IM+.5)       !  lon_index w.r.to 72x46 grid
CCC      JR=JREG(I,J)
C**** DETERMINE FRACTIONS FOR SURFACE TYPES AND COLUMN PRESSURE
      PLAND=FLAND(I,J)
      POICE=RSI(I,J)*(1.-PLAND)
      POCEAN=(1.-PLAND)-POICE
      PLAKE=FLAKE(I,J)
      PLICE=FLICE(I,J)
      PEARTH=FEARTH(I,J)
C****
      LS1_loc=LTROPO(I,J)+1  ! define stratosphere for radiation
      if (kradia.gt.0) then     ! rad forcing model
        do l=1,lm
          tlm(l) = T(i,j,l)*pk(l,i,j)
          shl(l) = QR(l,i,j)
          tauwc(l) = CLDinfo(l,1,i,j)
          tauic(l) = CLDinfo(l,2,i,j)
          SIZEWC(L)= CLDinfo(l,3,i,j)
          SIZEIC(L)= SIZEWC(L)
        end do
      else                      ! full model
C****
C**** DETERMINE CLOUDS (AND THEIR OPTICAL DEPTHS) SEEN BY RADIATION
C****
      CSS=0. ; CMC=0. ; CLDCV=0. ; DEPTH=0. ; OPTDW=0. ; OPTDI=0.
      DO L=1,LM
        PIJ=PLIJ(L,I,J)
        QSS=Q(I,J,L)/(RHSAV(L,I,J)+1.D-20)
        shl(L)=QSS
        IF(CLDSAV(L,I,J).LT.1.)
     *       shl(L)=(Q(I,J,L)-QSS*CLDSAV(L,I,J))/(1.-CLDSAV(L,I,J))
        TLm(L)=T(I,J,L)*PK(L,I,J)
        TAUSSL=0.
        TAUMCL=0.
        TAUWC(L)=0.
        TAUIC(L)=0.
        SIZEWC(L)=0.
        SIZEIC(L)=0.
        TOTCLD(L)=0.
C**** Determine large scale and moist convective cloud cover for radia
        IF (CLDSS(L,I,J).GT.RDSS(L,I,J)) THEN
          TAUSSL=TAUSS(L,I,J)
          shl(L)=QSS
          CSS=1.
          AJL(J,L,JL_SSCLD)=AJL(J,L,JL_SSCLD)+CSS
        END IF
        IF (CLDMC(L,I,J).GT.RDMC(I,J)) THEN
          CMC=1.
          AJL(J,L,JL_MCCLD)=AJL(J,L,JL_MCCLD)+CMC
          DEPTH=DEPTH+PDSIG(L,I,J)
          IF(TAUMC(L,I,J).GT.TAUSSL) THEN
            TAUMCL=TAUMC(L,I,J)
            ELHX=LHE
            IF(TLm(L).LE.TF) ELHX=LHS
            shl(L)=QSAT(TLm(L),ELHX,PMID(L,I,J))
          END IF
        END IF
        IF(TAUSSL+TAUMCL.GT.0.) THEN
             CLDCV=1.
          TOTCLD(L)=1.
          AJL(J,L,JL_TOTCLD)=AJL(J,L,JL_TOTCLD)+1.
          IF(TAUMCL.GT.TAUSSL) THEN
            SIZEWC(L)=CSIZMC(L,I,J)
            SIZEIC(L)=CSIZMC(L,I,J)
            IF(SVLAT(L,I,J).EQ.LHE) THEN
              TAUWC(L)=TAUMCL
              OPTDW=OPTDW+TAUWC(L)
              AJL(j,l,jl_wcld)=AJL(j,l,jl_wcld)+1.
            ELSE
              TAUIC(L)=TAUMCL
              OPTDI=OPTDI+TAUIC(L)
              AJL(j,l,jl_icld)=AJL(j,l,jl_icld)+1.
            END IF
          ELSE
            SIZEWC(L)=CSIZSS(L,I,J)
            SIZEIC(L)=CSIZSS(L,I,J)
            IF(SVLHX(L,I,J).EQ.LHE) THEN
              TAUWC(L)=TAUSSL
              OPTDW=OPTDW+TAUWC(L)
              AJL(j,l,jl_wcld)=AJL(j,l,jl_wcld)+1.
            ELSE
              TAUIC(L)=TAUSSL
              OPTDI=OPTDI+TAUIC(L)
              AJL(j,l,jl_icld)=AJL(j,l,jl_icld)+1.
            END IF
          END IF
          AJL(j,l,jl_wcod) =AJL(j,l,jl_wcod)+tauwc(l)
          AJL(j,l,jl_icod) =AJL(j,l,jl_icod)+tauic(l)
          AJL(j,l,jl_wcsiz)=AJL(j,l,jl_wcsiz)+sizewc(l)*tauwc(l)
          AJL(j,l,jl_icsiz)=AJL(j,l,jl_icsiz)+sizeic(l)*tauic(l)
        END IF
      END DO
C**** effective cloud cover diagnostics
         OPNSKY=1.-CLDCV
         DO IT=1,NTYPE
           AJ(J,J_PCLDSS,IT)=AJ(J,J_PCLDSS,IT)+CSS  *FTYPE(IT,I,J)
           AJ(J,J_PCLDMC,IT)=AJ(J,J_PCLDMC,IT)+CMC  *FTYPE(IT,I,J)
           AJ(J,J_CLDDEP,IT)=AJ(J,J_CLDDEP,IT)+DEPTH*FTYPE(IT,I,J)
           AJ(J,J_PCLD  ,IT)=AJ(J,J_PCLD  ,IT)+CLDCV*FTYPE(IT,I,J)
         END DO
CCC      AREG(JR,J_PCLDSS)=AREG(JR,J_PCLDSS)+CSS  *DXYP(J)
CCC      AREG(JR,J_PCLDMC)=AREG(JR,J_PCLDMC)+CMC  *DXYP(J)
CCC      AREG(JR,J_CLDDEP)=AREG(JR,J_CLDDEP)+DEPTH*DXYP(J)
CCC      AREG(JR,J_PCLD)  =AREG(JR,J_PCLD)  +CLDCV*DXYP(J)
         AREGIJ(1,I,J)=CSS  *DXYP(J)
         AREGIJ(2,I,J)=CMC  *DXYP(J)
         AREGIJ(3,I,J)=DEPTH*DXYP(J)
         AREGIJ(4,I,J)=CLDCV*DXYP(J)
         AIJ(I,J,IJ_PMCCLD)=AIJ(I,J,IJ_PMCCLD)+CMC
         AIJ(I,J,IJ_CLDCV) =AIJ(I,J,IJ_CLDCV) +CLDCV
         DO L=1,LLOW
           IF (TOTCLD(L).NE.1.) cycle
           AIJ(I,J,IJ_PCLDL)=AIJ(I,J,IJ_PCLDL)+1.
           exit
         end do
         DO L=LLOW+1,LMID
           IF (TOTCLD(L).NE.1.) cycle
           AIJ(I,J,IJ_PCLDM)=AIJ(I,J,IJ_PCLDM)+1.
           exit
         end do
         DO L=LMID+1,LHI
           IF (TOTCLD(L).NE.1.) cycle
           AIJ(I,J,IJ_PCLDH)=AIJ(I,J,IJ_PCLDH)+1.
           exit
         end do

         if(optdw.gt.0.) then
            AIJ(I,J,IJ_optdw)=AIJ(I,J,IJ_optdw)+optdw
            AIJ(I,J,IJ_wtrcld)=AIJ(I,J,IJ_wtrcld)+1.
         end if
         if(optdi.gt.0.) then
            AIJ(I,J,IJ_optdi)=AIJ(I,J,IJ_optdi)+optdi
            AIJ(I,J,IJ_icecld)=AIJ(I,J,IJ_icecld)+1.
         end if

         DO KR=1,NDIUPT
           IF (I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
             DO INCH=1,NRAD
               IH=1+MOD(JHOUR+INCH-1,24)
               ADIURN(IH,IDD_CL7,KR)=ADIURN(IH,IDD_CL7,KR)+TOTCLD(7)
               ADIURN(IH,IDD_CL6,KR)=ADIURN(IH,IDD_CL6,KR)+TOTCLD(6)
               ADIURN(IH,IDD_CL5,KR)=ADIURN(IH,IDD_CL5,KR)+TOTCLD(5)
               ADIURN(IH,IDD_CL4,KR)=ADIURN(IH,IDD_CL4,KR)+TOTCLD(4)
               ADIURN(IH,IDD_CL3,KR)=ADIURN(IH,IDD_CL3,KR)+TOTCLD(3)
               ADIURN(IH,IDD_CL2,KR)=ADIURN(IH,IDD_CL2,KR)+TOTCLD(2)
               ADIURN(IH,IDD_CL1,KR)=ADIURN(IH,IDD_CL1,KR)+TOTCLD(1)
               ADIURN(IH,IDD_CCV,KR)=ADIURN(IH,IDD_CCV,KR)+CLDCV
             END DO
           END IF
         END DO
      end if ! kradia le 0 (full model)
C****
C**** SET UP VERTICAL ARRAYS OMITTING THE I AND J INDICES
C****
C**** EVEN PRESSURES
      PLB(LM+1)=PEDN(LM+1,I,J)
      DO L=1,LM
        PLB(L)=PEDN(L,I,J)
C**** TEMPERATURES
C---- TLm(L)=T(I,J,L)*PK(L,I,J)     ! already defined
        IF(TLm(L).LT.130..OR.TLm(L).GT.370.) THEN
          WRITE(99,*) 'In Radia: Time,I,J,L,TL',ITime,I,J,L,TLm(L)
          WRITE(99,*) 'GTEMP:',GTEMP(:,:,I,J)
CCC       STOP 'In Radia: Temperature out of range'
          ICKERR=ICKERR+1
        END IF
C**** MOISTURE VARIABLES
C---- shl(L)=Q(I,J,L)        ! already defined
      END DO
C**** Radiative Equilibrium Layer data
      DO K=1,LM_REQ
        IF(RQT(K,I,J).LT.130..OR.RQT(K,I,J).GT.370.) THEN
        WRITE(99,*) 'In RADIA: Time,I,J,L,TL',ITime,I,J,LM+K,RQT(K,I,J)
CCC     STOP 'In Radia: RQT out of range'
        JCKERR=JCKERR+1
        END IF
        TLm(LM+K)=RQT(K,I,J)
        PLB(LM+k+1) = PLB0(k)
        shl(LM+k)    = shl0(k)
        tauwc(LM+k) = 0.
        tauic(LM+k) = 0.
        sizewc(LM+k)= 0.
        sizeic(LM+k)= 0.
      END DO
      if (kradia.gt.1) then
        do l=1,lm+lm_req
          tlm(l) = tlm(l) + Tchg(l,i,j)
          AFLX_ST(L,I,J,5)=AFLX_ST(L,I,J,5)+Tchg(L,I,J)
        end do
      end if
C**** Zenith angle and GROUND/SURFACE parameters
      COSZ=COSZA(I,J)
      TGO =GTEMP(1,1,I,J)+TF
      TGOI=GTEMP(1,2,I,J)+TF
      TGLI=GTEMP(1,3,I,J)+TF
      TGE =GTEMP(1,4,I,J)+TF
      TSL=TSAVG(I,J)
      SNOWOI=SNOWI(I,J)
      SNOWLI=SNOWLI_COM(I,J)
      SNOWE=SNOWE_COM(I,J)                    ! snow depth (kg/m**2)
      snow_frac(:) = fr_snow_rad_ij(:,i,j)    ! snow cover (1)
      AGESN(1)=SNOAGE(3,I,J)    ! land         ! ? why are these numbers
      AGESN(2)=SNOAGE(1,I,J)    ! ocean ice        so confusing ?
      AGESN(3)=SNOAGE(2,I,J)    ! land ice
c      print*,"snowage",i,j,SNOAGE(1,I,J)
C**** set up parameters for new sea ice and snow albedo
      hsn=snowoi/rhos
      if (poice.gt.0.) then
        hin=(ace1i+msi(i,j))/rhoi
        flags=flag_dsws(i,j)
        if (kradia .le. 0) then
          fmp=min(1.118d0*sqrt(pond_melt(i,j)/rhow),1d0)
             AIJ(I,J,IJ_FRMP) = AIJ(I,J,IJ_FRMP) + fmp*POICE
        else
          fmp = fmp_com(i,j)
        end if
        hmp=min(0.8d0*fmp,0.9d0*hin)
      else
        hin=0. ; flags=.FALSE. ; fmp=0. ; hmp=0.
      endif
C**** set up new lake depth parameter to incr. albedo for shallow lakes
      if (plake.gt.0) then
        lkdepth = MWL(I,J)/(RHOW*PLAKE*DXYP(J))
      end if
C****
      if (kradia .le. 0) then
        WEARTH=(WEARTH_COM(I,J)+AIEARTH(I,J))/(WFCS(I,J)+1.D-20)
        if (wearth.gt.1.) wearth=1.
      else                            ! rad.frc. model
        wearth = wsoil(i,j)
      end if
      DO K=1,11
        PVT(K)=VDATA(I,J,K)
      END DO
      WMAG=WSAVG(I,J)
C****
C*****************************************************
C     Main RADIATIVE computations, SOLAR and THERMAL
      CALL RCOMPX
C*****************************************************
      IF(I.EQ.IWRITE.AND.J.EQ.JWRITE) CALL WRITER(6,ITWRITE)
      CSZ2=COSZ2(I,J)
      if (kradia.gt.0) then  ! rad. forc. model; acc diagn
        do L=1,LM+LM_REQ+1
          AFLX_ST(L,I,J,1)=AFLX_ST(L,I,J,1)+SRUFLB(L)*CSZ2
          AFLX_ST(L,I,J,2)=AFLX_ST(L,I,J,2)+SRDFLB(L)*CSZ2
          AFLX_ST(L,I,J,3)=AFLX_ST(L,I,J,3)+TRUFLB(L)
          AFLX_ST(L,I,J,4)=AFLX_ST(L,I,J,4)+TRDFLB(L)
        end do
        if(kradia.eq.1) cycle
        do l=LS1_loc,ls1-1
          tchg(l,i,j) = tchg(l,i,j) + ( srfhrl(l)*csz2-srhra(l,i,j) +
     *      (-trfcrl(l)-trhra(l,i,j)) ) * coe(l)/p(i,j)
        end do
        do l=max(ls1,LS1_loc),lm+lm_req
          tchg(l,i,j) = tchg(l,i,j) + ( srfhrl(l)*csz2-srhra(l,i,j) +
     *      (-trfcrl(l)-trhra(l,i,j)) ) * coe(l)
        end do
        cycle
      else if (kradia.lt.0) then ! save i/o data for frc.runs
        fmp_com(i,j) = fmp                  ! input data
        wsoil(i,j) = wearth
        do L=1,LM
          QR(L,I,J) = shl(L)
          CLDinfo(L,1,I,J) = tauwc(L)
          CLDinfo(L,2,I,J) = tauic(L)
          CLDinfo(L,3,I,J) = sizeic(L)  ! sizeic=sizewc currently
        end do
        SRHRA(0,I,J)=SRNFLB(1)*CSZ2      ! output data (for adj frc)
        TRHRA(0,I,J)=-TRNFLB(1)
        do L=1,LM+LM_REQ
          SRHRA(L,I,J)=SRFHRL(L)*CSZ2
          TRHRA(L,I,J)=-TRFCRL(L)
        end do
      end if
C****
C**** Save relevant output in model arrays
C****
      FSF(1,I,J)=FSRNFG(1)   !  ocean
      FSF(2,I,J)=FSRNFG(3)   !  ocean ice
      FSF(3,I,J)=FSRNFG(4)   !  land ice
      FSF(4,I,J)=FSRNFG(2)   !  soil
      SRHR(0,I,J)=SRNFLB(1)
      TRHR(0,I,J)=STBO*(POCEAN*TGO**4+POICE*TGOI**4+PLICE*TGLI**4
     *  +PEARTH*TGE**4)-TRNFLB(1)
      DO L=1,LM
        SRHR(L,I,J)=SRFHRL(L)
        TRHR(L,I,J)=-TRFCRL(L)
      END DO
      DO LR=1,LM_REQ
        SRHRS(LR,I,J)= SRFHRL(LM+LR)
        TRHRS(LR,I,J)=-TRFCRL(LM+LR)
      END DO
      DO K=1,4
        SNFS(K,I,J)=SRNFLB(K+LM)
        TNFS(K,I,J)=TRNFLB(K+LM)-TRNFLB(1)
      END DO
      TRINCG(I,J)=TRDFLB(1)
      BTMPW(I,J)=BTEMPW-TF
      ALB(I,J,1)=SRNFLB(1)/(SRDFLB(1)+1.D-20)
      ALB(I,J,2)=PLAVIS
      ALB(I,J,3)=PLANIR
      ALB(I,J,4)=ALBVIS
      ALB(I,J,5)=ALBNIR
      ALB(I,J,6)=SRRVIS
      ALB(I,J,7)=SRRNIR
      ALB(I,J,8)=SRAVIS
      ALB(I,J,9)=SRANIR
      FSRDIR(I,J)=SRXVIS                  ! added by adf
C**** Save clear sky/tropopause diagnostics here
        AIJ(I,J,IJ_CLR_SRINCG)=AIJ(I,J,IJ_CLR_SRINCG)+
     +                                    OPNSKY*SRDFLB(1)*CSZ2
      DO IT=1,NTYPE
        AJ(J,J_CLRTOA,IT)=AJ(J,J_CLRTOA,IT)+OPNSKY*(SRNFLB(LM+LM_REQ+1)
     *     *CSZ2-TRNFLB(LM+LM_REQ+1))*FTYPE(IT,I,J)
        AJ(J,J_CLRTRP,IT)=AJ(J,J_CLRTRP,IT)+OPNSKY*(SRNFLB(LTROPO(I,J))
     *     *CSZ2-TRNFLB(LTROPO(I,J)))*FTYPE(IT,I,J)
        AJ(J,J_TOTTRP,IT)=AJ(J,J_TOTTRP,IT)+(SRNFLB(LTROPO(I,J))
     *     *CSZ2-TRNFLB(LTROPO(I,J)))*FTYPE(IT,I,J)
      END DO
CCC   AREG(JR,J_CLRTOA)=AREG(JR,J_CLRTOA)+OPNSKY*(SRNFLB(LM+LM_REQ+1)
CCC  *     *CSZ2-TRNFLB(LM+LM_REQ+1))*DXYP(J)
CCC   AREG(JR,J_CLRTRP)=AREG(JR,J_CLRTRP)+OPNSKY*
CCC  *     (SRNFLB(LTROPO(I,J))*CSZ2-TRNFLB(LTROPO(I,J)))*DXYP(J)
CCC   AREG(JR,J_TOTTRP)=AREG(JR,J_TOTTRP)+
CCC  *     (SRNFLB(LTROPO(I,J))*CSZ2-TRNFLB(LTROPO(I,J)))*DXYP(J)
      AREGIJ(5,I,J)=OPNSKY*(SRNFLB(LM+LM_REQ+1)
     *     *CSZ2-TRNFLB(LM+LM_REQ+1))*DXYP(J)
      AREGIJ(6,I,J)=OPNSKY*
     *     (SRNFLB(LTROPO(I,J))*CSZ2-TRNFLB(LTROPO(I,J)))*DXYP(J)
      AREGIJ(7,I,J)=
     *     (SRNFLB(LTROPO(I,J))*CSZ2-TRNFLB(LTROPO(I,J)))*DXYP(J)
C**** Save cloud top diagnostics here
      if (CLDCV.le.0.) go to 590
      AIJ(I,J,IJ_CLDTPPR)=AIJ(I,J,IJ_CLDTPPR)+plb(ltopcl+1)
      AIJ(I,J,IJ_CLDTPT)=AIJ(I,J,IJ_CLDTPT)+(tlb(ltopcl+1) - tf)
C**** Save cloud tau=1 related diagnostics here (opt.depth=1 level)
      tauup=0.
      DO L=LM,1,-1
        taucl=tauwc(l)+tauic(l)
        taudn=tauup+taucl
        if (taudn.gt.1.) then
          aij(i,j,ij_cldcv1)=aij(i,j,ij_cldcv1)+1.
          wtlin=(1.-tauup)/taucl
          aij(i,j,ij_cldt1t)=aij(i,j,ij_cldt1t)+( tlb(l+1)-tf +
     +          (tlb(l)-tlb(l+1))*wtlin )
          aij(i,j,ij_cldt1p)=aij(i,j,ij_cldt1p)+( plb(l+1)+
     +          (plb(l)-plb(l+1))*wtlin )
          go to 590
        end if
      end do
  590 continue

      END DO
C****
C**** END OF MAIN LOOP FOR I INDEX
C****
  600 CONTINUE
C****
C**** END OF MAIN LOOP FOR J INDEX
C****
!$OMP  END DO
!$OMP  END PARALLEL
CcOMP  END PARALLEL DO
      if(kradia.gt.0) return
C**** Stop if temperatures were out of range
      IF(ICKERR.GT.0)
     &     call stop_model('In Radia: Temperature out of range',11)
      IF(JCKERR.GT.0)  call stop_model('In Radia: RQT out of range',11)

C**** save all input data to disk if kradia<0
      if (kradia.lt.0) write(iu_rad) itime,T,RQT,TsAvg,QR,P,CLDinfo
     *  ,rsi,msi,(((GTEMP(1,k,i,j),k=1,4),i=1,im),j=1,jm),wsoil,wsavg
     *  ,snowi,snowli_com,snowe_com,snoage,fmp_com,flag_dsws,ltropo
     *  ,SRHRA,TRHRA,itime
C****
C**** ACCUMULATE THE RADIATION DIAGNOSTICS
C****
C       delayed accumulation to preserve order of summation
         DO J=1,JM
         DO I=1,IMAXJ(J)
           JR=JREG(I,J)
           AREG(JR,J_PCLDSS)=AREG(JR,J_PCLDSS)+AREGIJ(1,I,J)
           AREG(JR,J_PCLDMC)=AREG(JR,J_PCLDMC)+AREGIJ(2,I,J)
           AREG(JR,J_CLDDEP)=AREG(JR,J_CLDDEP)+AREGIJ(3,I,J)
           AREG(JR,J_PCLD)  =AREG(JR,J_PCLD)  +AREGIJ(4,I,J)
           AREG(JR,J_CLRTOA)=AREG(JR,J_CLRTOA)+AREGIJ(5,I,J)
           AREG(JR,J_CLRTRP)=AREG(JR,J_CLRTRP)+AREGIJ(6,I,J)
           AREG(JR,J_TOTTRP)=AREG(JR,J_TOTTRP)+AREGIJ(7,I,J)
         END DO
         END DO
C
         DO 780 J=1,JM
         DXYPJ=DXYP(J)
         DO L=1,LM
           DO I=1,IMAXJ(J)
             AJL(J,L,JL_SRHR)=AJL(J,L,JL_SRHR)+SRHR(L,I,J)*COSZ2(I,J)
             AJL(J,L,JL_TRCR)=AJL(J,L,JL_TRCR)+TRHR(L,I,J)
           END DO
         END DO
         DO 770 I=1,IMAXJ(J)
         CSZ2=COSZ2(I,J)
         JR=JREG(I,J)
         DO LR=1,LM_REQ
           ASJL(J,LR,3)=ASJL(J,LR,3)+SRHRS(LR,I,J)*CSZ2
           ASJL(J,LR,4)=ASJL(J,LR,4)+TRHRS(LR,I,J)
         END DO
         DO KR=1,NDIUPT
           IF (I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
             DO INCH=1,NRAD
               IH=1+MOD(JHOUR+INCH-1,24)
               ADIURN(IH,IDD_PALB,KR)=ADIURN(IH,IDD_PALB,KR)+
     *              (1.-SNFS(4,I,J)/S0)
               ADIURN(IH,IDD_GALB,KR)=ADIURN(IH,IDD_GALB,KR)+
     *              (1.-ALB(I,J,1))
               ADIURN(IH,IDD_ABSA,KR)=ADIURN(IH,IDD_ABSA,KR)+
     *              ((SNFS(4,I,J)-SNFS(1,I,J))*CSZ2-TNFS(4,I,J)
     *              +TNFS(1,I,J))
             END DO
           END IF
         END DO

         DO IT=1,NTYPE
         AJ(J,J_SRINCP0,IT)=AJ(J,J_SRINCP0,IT)+(S0*CSZ2)*FTYPE(IT,I,J)
         AJ(J,J_SRNFP0 ,IT)=AJ(J,J_SRNFP0 ,IT)+(SNFS(4,I,J)*CSZ2)*
     *          FTYPE(IT,I,J)
         AJ(J,J_SRINCG ,IT)=AJ(J,J_SRINCG ,IT)+(SRHR(0,I,J)*CSZ2/
     *          (ALB(I,J,1)+1.D-20))*FTYPE(IT,I,J)
         AJ(J,J_BRTEMP ,IT)=AJ(J,J_BRTEMP ,IT)+BTMPW(I,J) *FTYPE(IT,I,J)
         AJ(J,J_TRINCG ,IT)=AJ(J,J_TRINCG ,IT)+TRINCG(I,J)*FTYPE(IT,I,J)
         AJ(J,J_HSURF  ,IT)=AJ(J,J_HSURF  ,IT)-TNFS(4,I,J)*FTYPE(IT,I,J)
         AJ(J,J_SRNFP1 ,IT)=AJ(J,J_SRNFP1 ,IT)+SNFS(1,I,J)*CSZ2
     *          *FTYPE(IT,I,J)
         AJ(J,J_HATM   ,IT)=AJ(J,J_HATM   ,IT)-TNFS(1,I,J)*FTYPE(IT,I,J)
         END DO
         AREG(JR,J_SRINCP0)=AREG(JR,J_SRINCP0)+(S0*CSZ2)*DXYPJ
         AREG(JR,J_SRNFP0)=AREG(JR,J_SRNFP0)+(SNFS(4,I,J)*CSZ2)*DXYPJ
         AREG(JR,J_SRNFP1)=AREG(JR,J_SRNFP1)+(SNFS(1,I,J)*CSZ2)*DXYPJ
         AREG(JR,J_SRINCG)=AREG(JR,J_SRINCG)+
     *     (SRHR(0,I,J)*CSZ2/(ALB(I,J,1)+1.D-20))*DXYPJ
C**** Note: confusing because the types for radiation are a subset
         AJ(J,J_SRNFG,ITOCEAN)=AJ(J,J_SRNFG,ITOCEAN)+(FSF(1,I,J)*CSZ2)
     *        *FOCEAN(I,J)*(1.-RSI(I,J))
         AJ(J,J_SRNFG,ITLAKE) =AJ(J,J_SRNFG,ITLAKE) +(FSF(1,I,J)*CSZ2)
     *        * FLAKE(I,J)*(1.-RSI(I,J))
         AJ(J,J_SRNFG,ITEARTH)=AJ(J,J_SRNFG,ITEARTH)+(FSF(4,I,J)*CSZ2)
     *        *FEARTH(I,J)
         AJ(J,J_SRNFG,ITLANDI)=AJ(J,J_SRNFG,ITLANDI)+(FSF(3,I,J)*CSZ2)
     *        * FLICE(I,J)
         AJ(J,J_SRNFG,ITOICE )=AJ(J,J_SRNFG,ITOICE )+(FSF(2,I,J)*CSZ2)
     *        *FOCEAN(I,J)*RSI(I,J)
         AJ(J,J_SRNFG,ITLKICE)=AJ(J,J_SRNFG,ITLKICE)+(FSF(2,I,J)*CSZ2)
     *        * FLAKE(I,J)*RSI(I,J)
C****
         AREG(JR,J_HATM)  =AREG(JR,J_HATM)  - TNFS(1,I,J)      *DXYPJ
         AREG(JR,J_SRNFG) =AREG(JR,J_SRNFG) +(SRHR(0,I,J)*CSZ2)*DXYPJ
         AREG(JR,J_HSURF) =AREG(JR,J_HSURF) - TNFS(4,I,J)      *DXYPJ
         AREG(JR,J_BRTEMP)=AREG(JR,J_BRTEMP)+  BTMPW(I,J)      *DXYPJ
         AREG(JR,J_TRINCG)=AREG(JR,J_TRINCG)+ TRINCG(I,J)      *DXYPJ
         DO K=2,9
           JK=K+J_PLAVIS-2     ! accumulate 8 radiation diags.
           DO IT=1,NTYPE
             AJ(J,JK,IT)=AJ(J,JK,IT)+(S0*CSZ2)*ALB(I,J,K)*FTYPE(IT,I,J)
           END DO
           AREG(JR,JK)=AREG(JR,JK)+(S0*CSZ2)*ALB(I,J,K)*DXYPJ
         END DO
         AIJ(I,J,IJ_SRNFG)  =AIJ(I,J,IJ_SRNFG)  +(SRHR(0,I,J)*CSZ2)
         AIJ(I,J,IJ_BTMPW)  =AIJ(I,J,IJ_BTMPW)  +BTMPW(I,J)
         AIJ(I,J,IJ_SRREF)  =AIJ(I,J,IJ_SRREF)  +S0*CSZ2*ALB(I,J,2)
         AIJ(I,J,IJ_SRVIS)  =AIJ(I,J,IJ_SRVIS)  +S0*CSZ2*ALB(I,J,4)
         AIJ(I,J,IJ_TRNFP0) =AIJ(I,J,IJ_TRNFP0) - TNFS(4,I,J)
         AIJ(I,J,IJ_SRNFP0) =AIJ(I,J,IJ_SRNFP0) +(SNFS(4,I,J)*CSZ2)
         AIJ(I,J,IJ_SRINCG) =AIJ(I,J,IJ_SRINCG) +(SRHR(0,I,J)*CSZ2/
     *        (ALB(I,J,1)+1.D-20))
         AIJ(I,J,IJ_SRINCP0)=AIJ(I,J,IJ_SRINCP0)+(S0*CSZ2)
  770    CONTINUE
  780    CONTINUE
         DO L=1,LM
           DO I=1,IM
             DO J=J5S,J5N
               AIL(I,L,IL_REQ)=AIL(I,L,IL_REQ)+
     *              (SRHR(L,I,J)*COSZ2(I,J)+TRHR(L,I,J))*DXYP(J)
             END DO
             AIL(I,L,IL_R50N)=AIL(I,L,IL_R50N)+(SRHR(L,I,J50N)*COSZ2(I
     *            ,J50N)+TRHR(L,I,J50N))*DXYP(J50N)
             AIL(I,L,IL_R70N)=AIL(I,L,IL_R70N)+(SRHR(L,I,J70N)*COSZ2(I
     *            ,J70N)+TRHR(L,I,J70N))*DXYP(J70N)
           END DO
         END DO
C****
C**** Update radiative equilibrium temperatures
C****
      DO J=1,JM
        DO I=1,IMAXJ(J)
          DO LR=1,LM_REQ
            RQT(LR,I,J)=RQT(LR,I,J)+(SRHRS(LR,I,J)*COSZ2(I,J)
     *           +TRHRS(LR,I,J))*COE(LR+LM)
          END DO
        END DO
      END DO
C****
C**** Update other temperatures every physics time step
C****
  900 DO J=1,JM
        DO I=1,IMAXJ(J)
          DO L=1,LM
            T(I,J,L)=T(I,J,L)+(SRHR(L,I,J)*COSZ1(I,J)+TRHR(L,I,J))*
     *           COE(L)/(PLIJ(L,I,J)*PK(L,I,J))
          END DO
        END DO
      END DO
C**** daily diagnostics
      IH=1+JHOUR
      DO KR=1,NDIUPT
        ADIURN(IH,IDD_ISW,KR)=ADIURN(IH,IDD_ISW,KR)+
     *       S0*COSZ1(IJDD(1,KR),IJDD(2,KR))
      END DO

      RETURN
      END

      SUBROUTINE GHGHST(iu)
!@sum  reads history for nghg well-mixed greenhouse gases
!@auth R. Ruedy
!@ver  1.0

      USE RE001, only : nghg,nyrsghg,ghgyr1,ghgyr2,ghgam
      USE RADNCB, only : ghg_yr
      IMPLICIT NONE
      INTEGER iu,n,k
      CHARACTER*80 title

      write(6,*)
      do n=1,5
      read(iu,'(a)') title
      write(6,'(1x,a80)') title
      end do
      if(title(1:2).eq.'--') then                 ! older format
        read(iu,'(a)') title
        write(6,'(1x,a80)') title
      end if

      read(title,*) ghgyr1,(ghgam(k,1),k=1,nghg)
      ghgyr2=ghgyr1
      do n=2,nyrsghg
        read(iu,'(a)',end=20) title
        read(title,*) ghgyr2,(ghgam(k,n),k=1,nghg)
        if(ghg_yr.eq.0.or.abs(ghg_yr-ghgyr2).le.1)
     *      write(6,'(1x,a80)') title
        do k=1,nghg
          if(ghgam(k,n).lt.0.) ghgam(k,n)=ghgam(k,n-1)
        end do
      end do
   20 continue
      if(ghg_yr.ne.0.and.ghg_yr.ne.ghgyr2) write(6,'(1x,a80)') title
      write(*,*) 'read GHG table for years',ghgyr1,' - ',ghgyr2
      return
      end SUBROUTINE GHGHST

      subroutine getqma (iu,dglat,plb,dh2o,lm,jm)
!@sum  reads H2O production rates induced by CH4 (Tim Hall)
!@auth R. Ruedy
!@ver  1.0
      implicit none
      integer, parameter:: jma=18,lma=24
      integer m,iu,jm,lm,j,j1,j2,l,ll,ldn(lm),lup(lm)
      real*8 PLB(lm+1),dH2O(jm,lm,12),dglat(jm)
      real*4 pb(0:lma+1),h2o(jma,0:lma),xlat(jma),z(lma),dz(0:lma)
      character*100 title
      real*4 pdn,pup,w1,w2,dh,fracl

C**** read headers/latitudes
      read(iu,'(a)') title
      write(6,'(''0'',a100)') title
      read(iu,'(a)') title
      write(6,'(1x,a100)') title
      read(iu,'(a)') title
      write(6,'(1x,a100)') title
      read(title(10:100),*) (xlat(j),j=1,jma)

C**** read heights z(km) and data (kg/km^3/year)
      do m=1,12
        read(iu,'(a)') title
        write(6,'(1x,a100)') title
        do l=lma,1,-1
          read(iu,'(a)') title
          write(6,'(1x,a100)') title
          read(title,*) z(l),(H2O(j,l),j=1,jma)
        end do
        do j=1,jma
          h2o(j,0) = 0.
        end do

C**** Find edge heights and pressures
        dz(0) = 0.
        dz(1) = z(2)-z(1)
        do l=2,lma-1
           dz(l)=.5*(z(l+1)-z(l-1))
        end do
        dz(lma) = z(lma)-z(lma-1)

        pb(0) = plb(1)
        do l=1,lma
           Pb(l)=1000.*10.**(-(z(l)-.5*dz(l))/16.)
        end do
C**** extend both systems vertically to p=0
        pb(lma+1)=0.
        plb(lm+1)=0.

C**** Interpolate vertical resolution to model layers
        ldn(:) = 0
        do l=1,lm
          do while (pb(ldn(l)+1).ge.plb(l) .and. ldn(l).lt.lma)
            ldn(l)=ldn(l)+1
          end do
          lup(l)=ldn(l)
          do while (pb(lup(l)+1).gt.plb(l+1) .and. lup(l).lt.lma)
            lup(l)=lup(l)+1
          end do
        end do

C**** Interpolate (extrapolate) horizontally and vertically
        j2=2
        do j=1,jm
C**** coeff. for latitudinal linear inter/extrapolation
          do while (j2.lt.jma .and. dglat(j).gt.xlat(j2))
            j2 = j2+1
          end do
          j1 = j2-1
          w1 = (xlat(j2)-dglat(j))/(xlat(j2)-xlat(j1))
C**** for extrapolations, only use half the slope
          if(w1.gt.1.) w1=.5+.5*w1
          if(w1.lt.0.) w1=.5*w1
          w2 = 1.-w1
          do l=1,lm
            dh = 0.
            pdn = plb(l)
            if (lup(l).gt.0) then
              do ll=ldn(l),lup(l)
                pup = max(REAL(pb(ll+1),KIND=8),plb(l+1))
                fracl= (pdn-pup)/(pb(ll)-pb(ll+1))
                dh = dh+(w1*h2o(j1,ll)+w2*h2o(j2,ll))*fracl*dz(ll)
                pdn = pup
              end do
            end if
            dh2o(j,l,m) = 1.d-6*dh/1.74d0/365. !->(kg/m^2/ppm_CH4/day)
          end do
        end do
      end do
      return
      end subroutine getqma

      SUBROUTINE ORBPAR (YEAR, ECCEN,OBLIQ,OMEGVP)
!@sum ORBPAR calculates the three orbital parameters as a function of
!@+   YEAR.  The source of these calculations is: Andre L. Berger,
!@+   1978, "Long-Term Variations of Daily Insolation and Quaternary
!@+   Climatic Changes", JAS, v.35, p.2362.  Also useful is: Andre L.
!@+   Berger, May 1978, "A Simple Algorithm to Compute Long Term
!@+   Variations of Daily Insolation", published by Institut
!@+   D'Astronomie de Geophysique, Universite Catholique de Louvain,
!@+   Louvain-la Neuve, No. 18.
!@auth Gary Russell (with extra terms from D. Thresher)
C****
C**** Tables and equations refer to the first reference (JAS).  The
C**** corresponding table or equation in the second reference is
C**** enclosed in parentheses.
C****
      USE CONSTANT, only : twopi,PI180=>radian
      IMPLICIT NONE
C**** Input:
!@var YEAR   = years C.E. are positive, B.C.E are -ve (i.e 4BCE = -3)
      REAL*8, INTENT(IN) :: YEAR
C**** Output:
!@var ECCEN  = eccentricity of orbital ellipse
!@var OBLIQ  = latitude of Tropic of Cancer in degrees
!@var OMEGVP = longitude of perihelion =
!@+          = spatial angle from vernal equinox to perihelion
!@+            in degrees with sun as angle vertex
      REAL*8, INTENT(OUT) :: ECCEN,OBLIQ,OMEGVP
C**** Table 1 (2).  Obliquity relative to mean ecliptic of date: OBLIQD
      REAL*8, PARAMETER, DIMENSION(3,47) :: TABL1 = RESHAPE( (/
     1  -2462.22D0,  31.609970D0,  251.9025D0,
     2   -857.32D0,  32.620499D0,  280.8325D0,
     3   -629.32D0,  24.172195D0,  128.3057D0,
     4   -414.28D0,  31.983780D0,  292.7251D0,
     5   -311.76D0,  44.828339D0,   15.3747D0,
     6    308.94D0,  30.973251D0,  263.7952D0,
     7   -162.55D0,  43.668243D0,  308.4258D0,
     8   -116.11D0,  32.246689D0,  240.0099D0,
     9    101.12D0,  30.599442D0,  222.9725D0,
     O    -67.69D0,  42.681320D0,  268.7810D0,
     1     24.91D0,  43.836456D0,  316.7998D0,
     2     22.58D0,  47.439438D0,  319.6023D0,
     3    -21.16D0,  63.219955D0,  143.8050D0,
     4    -15.65D0,  64.230484D0,  172.7351D0,
     5     15.39D0,   1.010530D0,   28.9300D0,
     6     14.67D0,   7.437771D0,  123.5968D0,
     7    -11.73D0,  55.782181D0,   20.2082D0,
     8     10.27D0,   0.373813D0,   40.8226D0,
     9      6.49D0,  13.218362D0,  123.4722D0,
     O      5.85D0,  62.583237D0,  155.6977D0,
     1     -5.49D0,  63.593765D0,  184.6277D0,
     2     -5.43D0,  76.438309D0,  267.2771D0,
     3      5.16D0,  45.815262D0,   55.0196D0,
     4      5.08D0,   8.448301D0,  152.5268D0,
     5     -4.07D0,  56.792709D0,   49.1382D0,
     6      3.72D0,  49.747849D0,  204.6609D0,
     7      3.40D0,  12.058272D0,   56.5233D0,
     8     -2.83D0,  75.278214D0,  200.3284D0,
     9     -2.66D0,  65.241013D0,  201.6651D0,
     O     -2.57D0,  64.604294D0,  213.5577D0,
     1     -2.47D0,   1.647247D0,   17.0374D0,
     2      2.46D0,   7.811584D0,  164.4194D0,
     3      2.25D0,  12.207832D0,   94.5422D0,
     4     -2.08D0,  63.856659D0,  131.9124D0,
     5     -1.97D0,  56.155991D0,   61.0309D0,
     6     -1.88D0,  77.448837D0,  296.2073D0,
     7     -1.85D0,   6.801054D0,  135.4894D0,
     8      1.82D0,  62.209412D0,  114.8750D0,
     9      1.76D0,  20.656128D0,  247.0691D0,
     O     -1.54D0,  48.344406D0,  256.6113D0,
     1      1.47D0,  55.145462D0,   32.1008D0,
     2     -1.46D0,  69.000534D0,  143.6804D0,
     3      1.42D0,  11.071350D0,   16.8784D0,
     4     -1.18D0,  74.291306D0,  160.6835D0,
     5      1.18D0,  11.047742D0,   27.5932D0,
     6     -1.13D0,   0.636717D0,  348.1074D0,
     7      1.09D0,  12.844549D0,   82.6496D0/), (/3,47/) )
C**** Table 2 (4).  Precessional parameter: ECCEN sin(omega) (unused)
      REAL*8, PARAMETER, DIMENSION(3,46) :: TABL2 = RESHAPE(  (/
     1     .0186080D0,  54.646484D0,   32.012589D0,
     2     .0162752D0,  57.785370D0,  197.181274D0,
     3    -.0130066D0,  68.296539D0,  311.699463D0,
     4     .0098883D0,  67.659821D0,  323.592041D0,
     5    -.0033670D0,  67.286011D0,  282.769531D0,
     6     .0033308D0,  55.638351D0,   90.587509D0,
     7    -.0023540D0,  68.670349D0,  352.522217D0,
     8     .0014002D0,  76.656036D0,  131.835892D0,
     9     .0010070D0,  56.798447D0,  157.536392D0,
     O     .0008570D0,  66.649292D0,  294.662109D0,
     1     .0006499D0,  53.504456D0,  118.253082D0,
     2     .0005990D0,  67.023102D0,  335.484863D0,
     3     .0003780D0,  68.933258D0,  299.806885D0,
     4    -.0003370D0,  56.630219D0,  149.162415D0,
     5     .0003334D0,  86.256454D0,  283.915039D0,
     6     .0003334D0,  23.036499D0,  320.110107D0,
     7     .0002916D0,  89.395340D0,   89.083817D0,
     8     .0002916D0,  26.175385D0,  125.278732D0,
     9     .0002760D0,  69.307068D0,  340.629639D0,
     O    -.0002330D0,  99.906509D0,  203.602081D0,
     1    -.0002330D0,  36.686569D0,  239.796982D0,
     2     .0001820D0,  67.864838D0,  155.484787D0,
     3     .0001772D0,  99.269791D0,  215.494690D0,
     4     .0001772D0,  36.049850D0,  251.689606D0,
     5    -.0001740D0,  56.625275D0,  130.232391D0,
     6    -.0001240D0,  68.856720D0,  214.059708D0,
     7     .0001153D0,  87.266983D0,  312.845215D0,
     8     .0001153D0,  22.025970D0,  291.179932D0,
     9     .0001008D0,  90.405869D0,  118.013870D0,
     O     .0001008D0,  25.164856D0,   96.348694D0,
     1     .0000912D0,  78.818680D0,  160.318298D0,
     2     .0000912D0,  30.474274D0,   83.706894D0,
     3    -.0000806D0, 100.917038D0,  232.532120D0,
     4    -.0000806D0,  35.676025D0,  210.866943D0,
     5     .0000798D0,  81.957565D0,  325.487061D0,
     6     .0000798D0,  33.613159D0,  248.875565D0,
     7    -.0000638D0,  92.468735D0,   80.005234D0,
     8    -.0000638D0,  44.124329D0,    3.393823D0,
     9     .0000612D0, 100.280319D0,  244.424728D0,
     O     .0000612D0,  35.039322D0,  222.759552D0,
     1    -.0000603D0,  98.895981D0,  174.672028D0,
     2    -.0000603D0,  35.676025D0,  210.866943D0,
     3     .0000597D0,  87.248322D0,  342.489990D0,
     4     .0000597D0,  24.028381D0,   18.684967D0,
     5     .0000559D0,  86.630264D0,  324.737793D0,
     6     .0000559D0,  22.662689D0,  279.287354D0/), (/3,46/) )
C**** Table 3 (5).  Eccentricity: ECCEN (unused)
      REAL*8, PARAMETER, DIMENSION(3,42) :: TABL3 = RESHAPE( (/
     1     .01102940D0,   3.138886D0,  165.168686D0,
     2    -.00873296D0,  13.650058D0,  279.687012D0,
     3    -.00749255D0,  10.511172D0,  114.518250D0,
     4     .00672394D0,  13.013341D0,  291.579590D0,
     5     .00581229D0,   9.874455D0,  126.410858D0,
     6    -.00470066D0,   0.636717D0,  348.107422D0,
     7    -.00254464D0,  12.639528D0,  250.756897D0,
     8     .00231485D0,   0.991874D0,   58.574905D0,
     9    -.00221955D0,   9.500642D0,   85.588211D0,
     O     .00201868D0,   2.147012D0,  106.593765D0,
     1    -.00172371D0,   0.373813D0,   40.822647D0,
     2    -.00166112D0,  12.658154D0,  221.112030D0,
     3     .00145096D0,   1.010530D0,   28.930038D0,
     4     .00131342D0,  12.021467D0,  233.004639D0,
     5     .00101442D0,   0.373813D0,   40.822647D0,
     6    -.00088343D0,  14.023871D0,  320.509521D0,
     7    -.00083395D0,   6.277772D0,  330.337402D0,
     8     .00079475D0,   6.277772D0,  330.337402D0,
     9     .00067546D0,  27.300110D0,  199.373871D0,
     O    -.00066447D0,  10.884985D0,  155.340912D0,
     1     .00062591D0,  21.022339D0,  229.036499D0,
     2     .00059751D0,  22.009552D0,   99.823303D0,
     3    -.00053262D0,  27.300110D0,  199.373871D0,
     4    -.00052983D0,   5.641055D0,  342.229980D0,
     5    -.00052983D0,   6.914489D0,  318.444824D0,
     6     .00052836D0,  12.002811D0,  262.649414D0,
     7     .00051457D0,  16.788940D0,   84.855621D0,
     8    -.00050748D0,  11.647654D0,  192.181992D0,
     9    -.00049048D0,  24.535049D0,   75.027847D0,
     O     .00048888D0,  18.870667D0,  294.654541D0,
     1     .00046278D0,  26.026688D0,  223.159103D0,
     2     .00046212D0,   8.863925D0,   97.480820D0,
     3     .00046046D0,  17.162750D0,  125.678268D0,
     4     .00042941D0,   2.151964D0,  125.523788D0,
     5     .00042342D0,  37.174576D0,  325.784668D0,
     6     .00041713D0,  19.748917D0,  252.821732D0,
     7    -.00040745D0,  21.022339D0,  229.036499D0,
     8    -.00040569D0,   3.512699D0,  205.991333D0,
     9    -.00040569D0,   1.765073D0,  124.346024D0,
     O    -.00040385D0,  29.802292D0,   16.435165D0,
     1     .00040274D0,   7.746099D0,  350.172119D0,
     2     .00040068D0,   1.142024D0,  273.759521D0/), (/3,42/) )
C**** Table 4 (1).  Fundamental elements of the ecliptic: ECCEN sin(pi)
      REAL*8, PARAMETER, DIMENSION(3,19) :: TABL4 = RESHAPE( (/
     1     .01860798D0,   4.207205D0,   28.620089D0,
     2     .01627522D0,   7.346091D0,  193.788772D0,
     3    -.01300660D0,  17.857263D0,  308.307024D0,
     4     .00988829D0,  17.220546D0,  320.199637D0,
     5    -.00336700D0,  16.846733D0,  279.376984D0,
     6     .00333077D0,   5.199079D0,   87.195000D0,
     7    -.00235400D0,  18.231076D0,  349.129677D0,
     8     .00140015D0,  26.216758D0,  128.443387D0,
     9     .00100700D0,   6.359169D0,  154.143880D0,
     O     .00085700D0,  16.210016D0,  291.269597D0,
     1     .00064990D0,   3.065181D0,  114.860583D0,
     2     .00059900D0,  16.583829D0,  332.092251D0,
     3     .00037800D0,  18.493980D0,  296.414411D0,
     4    -.00033700D0,   6.190953D0,  145.769910D0,
     5     .00027600D0,  18.867793D0,  337.237063D0,
     6     .00018200D0,  17.425567D0,  152.092288D0,
     7    -.00017400D0,   6.186001D0,  126.839891D0,
     8    -.00012400D0,  18.417441D0,  210.667199D0,
     9     .00001250D0,   0.667863D0,   72.108838D0/), (/3,19/) )
C**** Table 5 (3).  General precession in longitude: psi
      REAL*8, PARAMETER, DIMENSION(3,78) :: TABL5 = RESHAPE( (/
     1    7391.0225890d0,  31.609974d0,   251.9025d0,
     2    2555.1526947d0,  32.620504d0,   280.8325d0,
     3    2022.7629188d0,  24.172203d0,   128.3057d0,
     4   -1973.6517951d0,   0.636717d0,   348.1074d0,
     5    1240.2321818d0,  31.983787d0,   292.7252d0,
     6     953.8679112d0,   3.138886d0,   165.1686d0,
     7    -931.7537108d0,  30.973257d0,   263.7951d0,
     8     872.3795383d0,  44.828336d0,    15.3747d0,
     9     606.3544732d0,   0.991874d0,    58.5749d0,
     O    -496.0274038d0,   0.373813d0,    40.8226d0,
     1     456.9608039d0,  43.668246d0,   308.4258d0,
     2     346.9462320d0,  32.246691d0,   240.0099d0,
     3    -305.8412902d0,  30.599444d0,   222.9725d0,
     4     249.6173246d0,   2.147012d0,   106.5937d0,
     5    -199.1027200d0,  10.511172d0,   114.5182d0,
     6     191.0560889d0,  42.681324d0,   268.7809d0,
     7    -175.2936572d0,  13.650058d0,   279.6869d0,
     8     165.9068833d0,   0.986922d0,    39.6448d0,
     9     161.1285917d0,   9.874455d0,   126.4108d0,
     O     139.7878093d0,  13.013341d0,   291.5795d0,
     1    -133.5228399d0,   0.262904d0,   307.2848d0,
     2     117.0673811d0,   0.004952d0,    18.9300d0,
     3     104.6907281d0,   1.142024d0,   273.7596d0,
     4      95.3227476d0,  63.219948d0,   143.8050d0,
     5      86.7824524d0,   0.205021d0,   191.8927d0,
     6      86.0857729d0,   2.151964d0,   125.5237d0,
     7      70.5893698d0,  64.230478d0,   172.7351d0,
     8     -69.9719343d0,  43.836462d0,   316.7998d0,
     9     -62.5817473d0,  47.439436d0,   319.6024d0,
     O      61.5450059d0,   1.384343d0,    69.7526d0,
     1     -57.9364011d0,   7.437771d0,   123.5968d0,
     2      57.1899832d0,  18.829299d0,   217.6432d0,
     3     -57.0236109d0,   9.500642d0,    85.5882d0,
     4     -54.2119253d0,   0.431696d0,   156.2147d0,
     5      53.2834147d0,   1.160090d0,    66.9489d0,
     6      52.1223575d0,  55.782177d0,    20.2082d0,
     7     -49.0059908d0,  12.639528d0,   250.7568d0,
     8     -48.3118757d0,   1.155138d0,    48.0188d0,
     9     -45.4191685d0,   0.168216d0,     8.3739d0,
     O     -42.2357920d0,   1.647247d0,    17.0374d0,
     1     -34.7971099d0,  10.884985d0,   155.3409d0,
     2      34.4623613d0,   5.610937d0,    94.1709d0,
     3     -33.8356643d0,  12.658184d0,   221.1120d0,
     4      33.6689362d0,   1.010530d0,    28.9300d0,
     5     -31.2521586d0,   1.983748d0,   117.1498d0,
     6     -30.8798701d0,  14.023871d0,   320.5095d0,
     7      28.4640769d0,   0.560178d0,   262.3602d0,
     8     -27.1960802d0,   1.273434d0,   336.2148d0,
     9      27.0860736d0,  12.021467d0,   233.0046d0,
     O     -26.3437456d0,  62.583231d0,   155.6977d0,
     1      24.7253740d0,  63.593761d0,   184.6277d0,
     2      24.6732126d0,  76.438310d0,   267.2772d0,
     3      24.4272733d0,   4.280910d0,    78.9281d0,
     4      24.0127327d0,  13.218362d0,   123.4722d0,
     5      21.7150294d0,  17.818769d0,   188.7132d0,
     6     -21.5375347d0,   8.359495d0,   180.1364d0,
     7      18.1148363d0,  56.792707d0,    49.1382d0,
     8     -16.9603104d0,   8.448301d0,   152.5268d0,
     9     -16.1765215d0,   1.978796d0,    98.2198d0,
     O      15.5567653d0,   8.863925d0,    97.4808d0,
     1      15.4846529d0,   0.186365d0,   221.5376d0,
     2      15.2150632d0,   8.996212d0,   168.2438d0,
     3      14.5047426d0,   6.771027d0,   161.1199d0,
     4     -14.3873316d0,  45.815258d0,    55.0196d0,
     5      13.1351419d0,  12.002811d0,   262.6495d0,
     6      12.8776311d0,  75.278220d0,   200.3284d0,
     7      11.9867234d0,  65.241008d0,   201.6651d0,
     8      11.9385578d0,  18.870667d0,   294.6547d0,
     9      11.7030822d0,  22.009553d0,    99.8233d0,
     O      11.6018181d0,  64.604291d0,   213.5577d0,
     1     -11.2617293d0,  11.498094d0,   154.1631d0,
     2     -10.4664199d0,   0.578834d0,   232.7153d0,
     3      10.4333970d0,   9.237738d0,   138.3034d0,
     4     -10.2377466d0,  49.747842d0,   204.6609d0,
     5      10.1934446d0,   2.147012d0,   106.5938d0,
     6     -10.1280191d0,   1.196895d0,   250.4676d0,
     7      10.0289441d0,   2.133898d0,   332.3345d0,
     8     -10.0034259d0,   0.173168d0,    27.3039d0/), (/3,78/) )
C****
      REAL*8 :: YM1950,SUMC,ARG,ESINPI,ECOSPI,PIE,PSI,FSINFD
      INTEGER :: I
C****
      YM1950 = YEAR-1950.
C****
C**** Obliquity from Table 1 (2):
C****   OBLIQ# = 23.320556 (degrees)             Equation 5.5 (15)
C****   OBLIQ  = OBLIQ# + sum[A cos(ft+delta)]   Equation 1 (5)
C****
      SUMC = 0.
      DO I=1,47
        ARG  = PI180*(YM1950*TABL1(2,I)/3600.+TABL1(3,I))
        SUMC = SUMC + TABL1(1,I)*COS(ARG)
      END DO
      OBLIQ = 23.320556D0 + SUMC/3600.
!      OBLIQ  = OBLIQ*PI180 ! not needed for output in degrees
C****
C**** Eccentricity from Table 4 (1):
C****   ECCEN sin(pi) = sum[M sin(gt+beta)]           Equation 4 (1)
C****   ECCEN cos(pi) = sum[M cos(gt+beta)]           Equation 4 (1)
C****   ECCEN = ECCEN sqrt[sin(pi)^2 + cos(pi)^2]
C****
      ESINPI = 0.
      ECOSPI = 0.
      DO I=1,19
        ARG    = PI180*(YM1950*TABL4(2,I)/3600.+TABL4(3,I))
        ESINPI = ESINPI + TABL4(1,I)*SIN(ARG)
        ECOSPI = ECOSPI + TABL4(1,I)*COS(ARG)
      END DO
      ECCEN  = SQRT(ESINPI*ESINPI+ECOSPI*ECOSPI)
C****
C**** Perihelion from Equation 4,6,7 (9) and Table 4,5 (1,3):
C****   PSI# = 50.439273 (seconds of degree)         Equation 7.5 (16)
C****   ZETA =  3.392506 (degrees)                   Equation 7.5 (17)
C****   PSI = PSI# t + ZETA + sum[F sin(ft+delta)]   Equation 7 (9)
C****   PIE = atan[ECCEN sin(pi) / ECCEN cos(pi)]
C****   OMEGVP = PIE + PSI + 3.14159                 Equation 6 (4.5)
C****
      PIE = ATAN2(ESINPI,ECOSPI)
      FSINFD = 0.
      DO I=1,78
        ARG    = PI180*(YM1950*TABL5(2,I)/3600.+TABL5(3,I))
        FSINFD = FSINFD + TABL5(1,I)*SIN(ARG)
      END DO
      PSI    = PI180*(3.392506D0+(YM1950*50.439273D0+FSINFD)/3600.)
      OMEGVP = MOD(PIE+PSI+.5*TWOPI,TWOPI)
      IF(OMEGVP.lt.0.)  OMEGVP = OMEGVP + TWOPI
      OMEGVP = OMEGVP/PI180  ! for output in degrees
C****
      RETURN
      END SUBROUTINE ORBPAR

      SUBROUTINE ORBIT (OBLIQ,ECCN,OMEGT,DAY,SDIST,SIND,COSD,LAMBDA)
!@sum ORBIT receives the orbital parameters and time of year, and
!@+   returns the distance from the sun and its declination angle.
!@+   The reference for the following calculations is: V.M.Blanco
!@+   and S.W.McCuskey, 1961, "Basic Physics of the Solar System",
!@+   pages 135 - 151.
!@auth Gary L. Russell and Robert J. Suozzo, 12/13/85

C**** Input
!@var OBLIQ = latitude of tropics in degrees
!@var ECCEN = eccentricity of the orbital ellipse
!@var OMEGT = angle from vernal equinox to perihelion in degrees
!@var DAY   = day of the year in days; 0 = Jan 1, hour 0

C**** Constants:
!@param VERQNX = occurence of vernal equinox = day 79 = Mar 21 hour 0

C**** Intermediate quantities:
!@var PERIHE = perihelion during the year in temporal radians
!@var MA     = mean anomaly in temporal radians = 2 JDAY/365 - PERIHE
!@var EA     = eccentric anomaly in radians
!@var TA     = true anomaly in radians
!@var BSEMI  = semi minor axis in units of the semi major axis
!@var GREENW = longitude of Greenwich in the Earth's reference frame

C**** Output:
!@var SDIST = square of distance to the sun in units of semi major axis
!@var SIND = sine of the declination angle
!@var COSD = cosine of the declination angle
!@var LAMBDA = sun longitude in Earth's rotating reference frame (OBS)
      USE CONSTANT, only : pi,radian,edpery
      IMPLICIT NONE
      REAL*8, PARAMETER :: VERQNX = 79.
      REAL*8, INTENT(IN) :: OBLIQ,ECCN,OMEGT,DAY
      REAL*8, INTENT(OUT) :: SIND,COSD,SDIST,LAMBDA

      REAL*8 MA,OMEGA,DOBLIQ,ECCEN,PERIHE,EA,DEA,BSEMI,COSEA
     *     ,SINEA,TA,SUNX,SUNY,GREENW,SINDD
C****
      OMEGA=OMEGT*radian
      DOBLIQ=OBLIQ*radian
      ECCEN=ECCN
C****
C**** Determine time of perihelion using Kepler's equation:
C**** PERIHE-VERQNX = OMEGA - ECCEN sin(OMEGA)
C****
      PERIHE = OMEGA-ECCEN*SIN(OMEGA)+VERQNX*2.*PI/EDPERY
C     PERIHE = DMOD(PERIHE,2.*PI)
      MA = 2.*PI*DAY/EDPERY - PERIHE
      MA = DMOD(MA,2.*PI)
C****
C**** Numerically solve Kepler's equation: MA = EA - ECCEN sin(EA)
C****
      EA = MA+ECCEN*(SIN(MA)+ECCEN*SIN(2.*MA)/2.)
  110 DEA = (MA-EA+ECCEN*SIN(MA))/(1.-ECCEN*COS(EA))
      EA = EA+DEA
      IF (DABS(DEA).GT.1.D-8)  GO TO 110
C****
C**** Calculate the distance to the sun and the true anomaly
C****
      BSEMI = DSQRT(1.-ECCEN*ECCEN)
      COSEA = COS(EA)
      SINEA = SIN(EA)
      SDIST  = (1.-ECCEN*COSEA)*(1.-ECCEN*COSEA)
      TA = DATAN2(SINEA*BSEMI,COSEA-ECCEN)
C****
C**** Change the reference frame to be the Earth's equatorial plane
C**** with the Earth at the center and the positive x axis parallel to
C**** the ray from the sun to the Earth were it at vernal equinox.
C**** The distance from the current Earth to that ray (or x axis) is:
C**** DIST sin(TA+OMEGA).  The sun is located at:
C****
C**** SUN    = (-DIST cos(TA+OMEGA),
C****           -DIST sin(TA+OMEGA) cos(OBLIQ),
C****            DIST sin(TA+OMEGA) sin(OBLIQ))
C**** SIND   = sin(TA+OMEGA) sin(OBLIQ)
C**** COSD   = sqrt(1-SIND**2)
C**** LAMBDA = atan[tan(TA+OMEGA) cos(OBLIQ)] - GREENW
C**** GREENW = 2*3.14159 DAY (EDPERY-1)/EDPERY
C****
      SINDD = SIN(TA+OMEGA)*SIN(DOBLIQ)
      COSD = DSQRT(1.-SINDD*SINDD)
      SIND = SINDD
C     GREENW = 2.*PI*(DAY-VERQNX)*(EDPERY+1.)/EDPERY
C     SUNX = -COS(TA+OMEGA)
C     SUNY = -SIN(TA+OMEGA)*COS(DOBLIQ)
      LAMBDA = 0. ! just to keep the compiler happy
C     LAMBDA = DATAN2(SUNY,SUNX)-GREENW
C     LAMBDA = DMOD(LAMBDA,2.*PI)
C****
      RETURN
      END SUBROUTINE ORBIT


      subroutine updveg (year)
!@sum  reads appropriate crops data and updates the vegetation file
!@auth R. Ruedy
!@ver  1.0
      USE FILEMANAGER
      USE MODEL_COM, only : im,jm,vdata      
      USE GEOM, only : imaxj
      implicit none
      integer, intent(in) :: year

      real*8 wt,crop(im,jm),crops         ! temporary vars

      integer :: year1,year2,year_old=-1, iu, i,j,k
      real*8 vdata0(im,jm,11),crop1(im,jm),crop2(im,jm)  ! to limit i/o
      save    year1,year2,year_old,vdata0,crop1,crop2    ! to limit i/o

      character*80 title
      real*4 crop4(im,jm)

C**** check whether update is needed
      if (year.eq.year_old) return

C**** first iteration actions:
      if (year_old.lt.0) then
C****     check whether a no-crops vege-file was used
        do j=2,jm-1
        do i=1,im
           if(vdata(i,j,9).gt.0.)
     *     call stop_model('updveg: use no_crops_VEG_file',255)
        end do
        end do
C****     open and read input file
        call openunit('CROPS',iu,.true.,.true.)
        read(iu) title,crop4
        read(title,*) year1
        crop1=crop4 ; crop2=crop4 ; year2=year1
        if (year1.ge.year)          year2=year+1
C****     save orig. (no-crop) vdata to preserve restart-independence
        vdata0 = vdata
      end if

      wt=0.
      do while (year2.lt.year)
         year1 = year2 ; crop1 = crop2
         read (iu,end=10) title,crop4
         read(title,*) year2
         crop2 = crop4
      end do
      wt = (year-year1)/(real(year2-year1,kind=8))
   10 continue
      write(6,*) 'Using crops data from year',year1+wt*(year2-year1)

C**** Modify the vegetation fractions
      do j=1,jm
      do i=1,imaxj(j)
         if (crop1(i,j).ge.0.) then
            crops = crop1(i,j) + wt*(crop2(i,j)-crop1(i,j))
            do k=1,11
              vdata(i,j,k) = vdata0(i,j,k)*(1.-crops)
            end do
              vdata(i,j,9) = crops
         end if
      end do
      end do

      year_old = year

      return
      end subroutine updveg

