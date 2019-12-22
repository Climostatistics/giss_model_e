#include "rundeck_opts.h"

      SUBROUTINE diag_OCEAN
!@sum  diag_OCEAN prints out diagnostics for ocean
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0
C**** Note this is an incorporation and modification of the stand alone
C**** ocean diagnostic programs from Gary. All diagnostics are on the
C**** ocean grid.
      USE CONSTANT, only : undef,teeny
      USE MODEL_COM, only : xlabel,lrunid,jmon0,jyear0,idacc,jdate0
     *     ,amon0,jdate,amon,jyear
      USE OCEAN, only : im,jm,lmo,ndyno,dts,dto,imaxj,lmm,ze
      USE DAGCOM, only : qdiag,acc_period
      USE ODIAG
      USE FILEMANAGER, only : openunit
      IMPLICIT NONE
      INTEGER I,J,L,N,NOL(LMO),JEQ,JDLAT,KXLB
      REAL*8 DLON

C**** Calculate latitudes
      JEQ = JM/2
      JDLAT = NINT(180./(JM-1))
      DO J=1,JM
        FLAT(J,1)  = JDLAT*(J-JEQ-0.5)  ! primary grid
      END DO
      FLAT(1,2)=-90.
      DO J=1,JM-1
        FLAT(J+1,2)  = JDLAT*(J-JEQ)  ! secondary grid  (shifted by 1)
      END DO
C**** Calculate longitudes
      DLON=360./REAL(IM,KIND=8)
      DO I=1,IM
        FLON(I,1)  = -180.+DLON*(I-0.5)  ! primary grid
      END DO
      DO I=1,IM
        FLON(I,2)  = -180.+DLON*I        ! secondary grid
      END DO

C**** determine label to be added to all titles
      KXLB = INDEX(XLABEL(1:11),'(')-1
      IF(KXLB.le.0) KXLB = 10
      XLB = ' '
      XLB(1:13)=acc_period(1:3)//' '//acc_period(4:12)
      XLB = TRIM(XLB)//" "//XLABEL(1:KXLB)

C**** Open output files
      IF(QDIAG) then
        call open_ij(trim(acc_period)//'.oij'//XLABEL(1:LRUNID),im,jm)
        call open_jl(trim(acc_period)//'.ojl'//XLABEL(1:LRUNID),jm,lmo,0
     *       ,flat)
        call open_il(trim(acc_period)//'.oil'//XLABEL(1:LRUNID),im,lmo,0
     *       )
        call openunit(trim(acc_period)//'.otj'//XLABEL(1:LRUNID),iu_otj,
     *       .FALSE.,.FALSE.)
      END IF
C**** Call diagnostic calculating routines
      CALL OIJOUT  ! lat-lon diags
      CALL OSFOUT  ! overturning streamfunction diags
      CALL OJLOUT  ! lat-height and lon-height sections and means
      WRITE (6,*)
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
      CALL OTJOUT  ! transports
      CALL STROUT  ! strait diags
C**** Miscellaneous diags for print out

C****
C**** vertical mean diagnostics
C****
      WRITE (6,*)
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
C**** Calculate number of points in average
      NOL = 0.
      DO J=1,JM
      DO I=1,IMAXJ(J)
      DO L=1,LMM(I,J)
        NOL(L)=NOL(L)+1
      END DO
      END DO
      END DO

      WRITE(6,*) " Ocean Mean quantities:"
      WRITE(6,*) " Level   Rho0       Temp    Salinity"
      DO L=1,LMO
        WRITE(6,'(2X,I4,10F10.3)') L,OL(L,L_RHO)/(4.*NOL(L)*IDACC(1)),
     *       OL(L,L_TEMP)/(4.*NOL(L)*IDACC(1)),
     *       1d3*OL(L,L_SALT)/(4.*NOL(L)*IDACC(1))
      END DO
      call sys_flush(6)
C****
      IF (QDIAG) THEN
        call close_ij
        call close_jl
        call close_il
      END IF
C****
      RETURN
  907 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)
      END SUBROUTINE diag_OCEAN

      SUBROUTINE OIJOUT
!@sum  OIJOUT prints out lat-lon diagnostics for ocean
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0
      USE CONSTANT, only : undef,teeny
      USE MODEL_COM, only : xlabel,lrunid,jmon0,jyear0,idacc,jdate0
     *     ,amon0,jdate,amon,jyear
#ifdef TRACERS_OCEAN
      USE TRACER_COM, only : ntm,trw0,trname,ntrocn
#endif
      USE OCEAN, only : im,jm,lmo,focean,dxypo,ndyno,dts,dto
     *     ,imaxj,lmm,ze,dxvo,dypo
      USE DAGCOM, only : qdiag,acc_period
      USE STRAITS, only : nmst,wist,dist,lmst,name_st
      USE ODIAG
#ifdef TRACERS_OCEAN
      USE TRACER_DIAG_COM, only : to_per_mil
#endif
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,JM) :: Q,SFIJM,SFIJS,ADENOM
      REAL*8, DIMENSION(IM,JM,LMO) :: Q3
      REAL*8, DIMENSION(LMO,NMST) :: AS
      INTEGER I,J,K,L,NS,N,KB,IJGRID,IP1,k1,KK
     *     ,LMSTMIN(LMO),SUMORMN(LMO),JEQ,JDLAT,KXLB
      INTEGER :: LMINEF=1, LMAXEF=LMO, KVMF(3) = (/ 3, 6, 9/),
     *           LMINMF=1, LMAXMF=1,   KCMF(3) = (/ 3, 6, 9/),
     *           LMINSF=1, LMAXSF=LMO, KVDC(3) = (/ 3, 6, 9/)
      REAL*8 GOS,SOS,FAC,FACST,GSMAX,GSMIN,CKMIN,CKMAX,ACMIN,ACMAX
     *     ,TSUM,TEMGS,QJ(JM),QSUM,MQ,DLON,byiacc
      CHARACTER NAME(KOLNST)*40,TITLE*80,lname*50,sname*30,units*50
      character*50 :: unit_string

      QJ=0.
      QSUM=0.
      IJGRID=1
      Q=0.
      Q3=0.

C**** lat/lon diagnostics
      WRITE (6,*)
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
  907 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)

      IF(QDIAG) then
C****
C**** Ocean Potential Temperature (C)
C****
      LNAME="OCEAN POTENTIAL TEMPERATURE"
      UNITS="C"
      SNAME="oc_temp"
      TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
      TITLE(51:80)=XLB
C**** Loop over layers
      DO L=1,LMO
      DO J=1,JM
      DO I=1,IMAXJ(J)
        Q(I,J) = UNDEF
        IF(FOCEAN(I,J).gt..5 .and. OIJL(I,J,L,IJL_MO).gt.0.)  THEN
          GOS = OIJL(I,J,L,IJL_G0M) / (OIJL(I,J,L,IJL_MO)*DXYPO(J))
          SOS = OIJL(I,J,L,IJL_S0M) / (OIJL(I,J,L,IJL_MO)*DXYPO(J))
          Q(I,J) = TEMGS(GOS,SOS)
        END IF
      END DO
      END DO
      Q(2:IM,JM)=Q(1,JM)
      Q(2:IM,1)=Q(1,1)
      WRITE (LNAME(40:47),'(A5,I3)') 'Level',L
      WRITE (TITLE(40:47),'(A5,I3)') 'Level',L
      CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID)
      END DO
C****
C**** Ocean Salinity (per mil)
C****
      LNAME="OCEAN SALINITY"
      UNITS="psu"
      SNAME="oc_salt"
      TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
      TITLE(51:80)=XLB
C**** Loop over layers
      DO L=1,LMO
      DO J=1,JM
      DO I=1,IMAXJ(J)
        Q(I,J) = UNDEF
        IF(FOCEAN(I,J).gt..5 .and. OIJL(I,J,L,IJL_MO).gt.0.)
     *       Q(I,J) = 1d3*OIJL(I,J,L,IJL_S0M) / (OIJL(I,J,L,IJL_MO)
     *       *DXYPO(J))
      END DO
      END DO
      Q(2:IM,JM)=Q(1,JM)
      Q(2:IM,1)=Q(1,1)
      WRITE (LNAME(40:47),'(A5,I3)') 'Level',L
      WRITE (TITLE(40:47),'(A5,I3)') 'Level',L
      CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID)
      END DO
#ifdef TRACERS_OCEAN
C****
C**** Ocean Tracers 
C****
      DO N=1,NTM
      LNAME="OCEAN "//trname(n)
      if (to_per_mil(n).gt.0) THEN
        UNITS="per mil"
      ELSE
        UNITS=unit_string(ntrocn(n),'kg/kg')
      END IF
      SNAME="oc_"//trim(trname(n))
      TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
      TITLE(51:80)=XLB
C**** Loop over layers
      DO L=1,LMO
      DO J=1,JM
      DO I=1,IMAXJ(J)
        Q(I,J) = UNDEF
        IF(FOCEAN(I,J).gt..5 .and. OIJL(I,J,L,IJL_MO).gt.0.) THEN
        if (to_per_mil(n).gt.0) THEN
          Q(I,J)=1d3*(TOIJL(I,J,L,TOIJL_CONC,N)/((OIJL(I,J,L,IJL_MO)
     *         *DXYPO(J)-OIJL(I,J,L,IJL_S0M))*trw0(n))-1.)
        else
          Q(I,J)=10.**(-ntrocn(n))*TOIJL(I,J,L,TOIJL_CONC,N)/
     *         (OIJL(I,J,L,IJL_MO)*DXYPO(J))
        end if
        END IF
      END DO
      END DO
      Q(2:IM,JM)=Q(1,JM)
      Q(2:IM,1)=Q(1,1)
      WRITE (LNAME(40:47),'(A5,I3)') 'Level',L
      WRITE (TITLE(40:47),'(A5,I3)') 'Level',L
      CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID)
      END DO
      END DO
#endif      
C****
C**** East-West or North-South Velocities (cm/s)
C****
  300 K=IJL_MFU
      LNAME="EAST-WEST VELOCITY"
      UNITS="cm/s"
      SNAME="uvel"//char(l+48)
      Q = 0.
      DO J=1,JM
        I=IM
        DO IP1=1,IMAXJ(J)
          MQ = 0.
          DO L=LMINMF,LMAXMF
            MQ = MQ + (OIJL(I,J,L,IJL_MO)+OIJL(IP1,J,L,IJL_MO))
            Q(I,J) = Q(I,J) + OIJL(I,J,L,K)
          END DO
          Q(I,J) = 1d2 * 4.* Q(I,J) / (MQ*NDYNO*DYPO(J)+teeny)
          I=IP1
        END DO
      END DO
      Q(2:IM,JM)=Q(1,JM)
      Q(2:IM,1)=Q(1,1)
      TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
      IF(LMINMF.eq.LMAXMF) WRITE (TITLE(39:41),'(I3)') LMINMF
      IF(LMINMF.lt.LMAXMF) WRITE (TITLE(39:46),'(I3,A2,I3)') LMINMF," -"
     *     ,LMAXMF
      IF(LMINMF.eq.LMAXMF) WRITE (LNAME(39:41),'(I3)') LMINMF
      IF(LMINMF.lt.LMAXMF) WRITE (LNAME(39:46),'(I3,A2,I3)') LMINMF," -"
     *     ,LMAXMF
      TITLE(51:80)=XLB
      CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID)

      K=IJL_MFV
      LNAME="NORTH-SOUTH VELOCITY"
      UNITS="cm/s"
      SNAME="vvel"//char(l+48)
      Q = 0.
      DO J=1,JM-1
        DO I=1,IMAXJ(J)
          MQ = 0.
          DO L=LMINMF,LMAXMF
            MQ = MQ + (OIJL(I,J,L,IJL_MO)+OIJL(I,J+1,L,IJL_MO))
            Q(I,J) = Q(I,J) + OIJL(I,J,L,K)
          END DO
          Q(I,J) = 1d2 * 4.* Q(I,J) / (MQ*NDYNO*DXVO(J)+teeny)
        END DO
      END DO
      Q(2:IM,1)=Q(1,1)
      TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
      IF(LMINMF.eq.LMAXMF) WRITE (TITLE(39:41),'(I3)') LMINMF
      IF(LMINMF.lt.LMAXMF) WRITE (TITLE(39:46),'(I3,A2,I3)') LMINMF," -"
     *     ,LMAXMF
      IF(LMINMF.eq.LMAXMF) WRITE (LNAME(39:41),'(I3)') LMINMF
      IF(LMINMF.lt.LMAXMF) WRITE (LNAME(39:46),'(I3,A2,I3)') LMINMF," -"
     *     ,LMAXMF
      TITLE(51:80)=XLB
      CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,2.)
C****
C**** Vertical Mass Flux (10^-2 kg/s*m^2)
C****
      DO K=1,3  !LMO
c     IF(KVMF(K).le.0)  GO TO 370
        L =KVMF(K)
        LNAME="VERTICAL MASS FLUX"
        UNITS="10^-2 kg/s*m^2"
        SNAME="vert_mflx"//char(L+48)
        DO J=1,JM
          DO I=1,IMAXJ(J)
            Q(I,J) = 2d2*OIJL(I,J,L,IJL_MFW) / (IDACC(1)*NDYNO*DXYPO(J))
          END DO
        END DO
        Q(2:IM,JM)=Q(1,JM)
        Q(2:IM,1)=Q(1,1)
        TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
        WRITE (TITLE(40:47),'(A5,I3)') "Level",L
        WRITE (LNAME(40:47),'(A5,I3)') "Level",L
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID)
      END DO
C****
C**** East-West or North-South Heat Flux (10^15 W)
C****
      DO K=IJL_GFLX,IJL_GFLX+1
c      IF(.not.QL(K))  GO TO 440
        IF (K.eq.IJL_GFLX) THEN
          LNAME="EAST-WEST HEAT FLUX"
          UNITS="10^15 W"
          SNAME="ew_hflx"
        ELSE
          LNAME="NORTH-SOUTH HEAT FLUX"
          UNITS="10^15 W"
          SNAME="ns_hflx"
        END IF
        Q = 0.
        DO J=1,JM
        DO I=1,IMAXJ(J)
          DO L=LMINEF,LMAXEF
            Q(I,J) = Q(I,J) + OIJL(I,J,L,K)
          END DO
          Q(I,J) = 1d-15*Q(I,J) / (IDACC(1)*DTS)
        END DO
        END DO
        Q(2:IM,JM)=Q(1,JM)
        Q(2:IM,1)=Q(1,1)
        TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
        IF(LMINEF.eq.LMAXEF) WRITE (TITLE(39:41),'(I3)') LMINEF
        IF(LMINEF.lt.LMAXEF) WRITE (TITLE(39:46),'(I3,A2,I3)') LMINEF,
     *       " -",LMAXEF
        IF(LMINEF.eq.LMAXEF) WRITE (LNAME(39:41),'(I3)') LMINEF
        IF(LMINEF.lt.LMAXEF) WRITE (LNAME(39:46),'(I3,A2,I3)') LMINEF,
     *       " -",LMAXEF
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID)
      END DO
C****
C**** East-West or North-South Salt Flux (10^6 kg/s)
C****
      DO K=IJL_SFLX,IJL_SFLX+1
c      IF(.not.QL(K))  GO TO 540
        IF (K.eq.IJL_SFLX) THEN
          LNAME="EAST-WEST SALT FLUX"
          UNITS="10^6 kg/s"
          SNAME="ew_sflx"
        ELSE
          LNAME="NORTH-SOUTH SALT FLUX"
          UNITS="10^6 kg/s"
          SNAME="ns_sflx"
        END IF
        Q = 0.
        DO J=1,JM
        DO I=1,IMAXJ(J)
          DO L=LMINSF,LMAXSF
            Q(I,J) = Q(I,J) + OIJL(I,J,L,K)
          END DO
          Q(I,J) = 1d-6*Q(I,J) / (IDACC(1)*DTS)
        END DO
        END DO
        Q(2:IM,JM)=Q(1,JM)
        Q(2:IM,1)=Q(1,1)
        TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
        IF(LMINSF.eq.LMAXSF) WRITE (TITLE(39:41),'(I3)') LMINSF
        IF(LMINSF.lt.LMAXSF) WRITE (TITLE(39:46),'(I3,A2,I3)') LMINSF,
     *       " -",LMAXSF
        IF(LMINSF.eq.LMAXSF) WRITE (LNAME(39:41),'(I3)') LMINSF
        IF(LMINSF.lt.LMAXSF) WRITE (LNAME(39:46),'(I3,A2,I3)') LMINSF,
     *       " -",LMAXSF
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID)
      END DO
C****
C**** Gent-McWilliams fluxes (10^-2 kg/s*m )
C****
      DO K=1,3 
        L =KCMF(K)
        DO KK=0,2 
          SELECT CASE (KK)
          CASE (0)      ! E-W fluxes
            LNAME="GM/EDDY E-W HEAT FLUX"
            UNITS="10^-15 W"
            SNAME="gm_ew_hflx"//char(l+48)
            DO J=1,JM
              DO I=1,IMAXJ(J)
                Q(I,J) = 1d-15*OIJL(I,J,L,KK+IJL_GGMFL)/(IDACC(1)*DTS)
              END DO
            END DO    
          CASE (1)  ! N-S fluxes
            LNAME="GM/EDDY N-S HEAT FLUX"
            UNITS="10^-15 W"
            SNAME="gm_ns_hflx"//char(l+48)
            DO J=1,JM
              DO I=1,IMAXJ(J)
                Q(I,J) = 1d-15*OIJL(I,J,L,KK+IJL_GGMFL)/(IDACC(1)*DTS)
              END DO
            END DO    
          CASE (2)    !  Vertical fluxes
            LNAME="GM/EDDY VERTICAL HEAT FLUX"
            UNITS="W/m^2"
            SNAME="gm_vt_hflx"//char(l+48)
            DO J=1,JM
              DO I=1,IMAXJ(J)
                Q(I,J)=OIJL(I,J,L,KK+IJL_GGMFL)/(IDACC(1)*DTS*DXYPO(J))
              END DO
            END DO    
          END SELECT

          Q(2:IM,JM)=Q(1,JM)
          Q(2:IM,1)=Q(1,1)
          TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
          WRITE (TITLE(40:47),'(A5,I3)') "Level",L
          WRITE (LNAME(40:47),'(A5,I3)') "Level",L
          TITLE(51:80)=XLB
          CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID)
        END DO
      END DO
C****
C**** Vertical Diffusion Coefficients (cm /s)
C****
      DO K=1,3 
        L =KCMF(K)
        LNAME="VERT. MOM. DIFF."
        UNITS="cm^2/s"
        SNAME="kvm"//char(l+48)
        Q=UNDEF
        DO J=1,JM
        DO I=1,IMAXJ(J)
          IF (OIJL(I,J,L+1,IJL_MO).gt.0)
     *       Q(I,J)=1d4*.00097d0**2*OIJL(I,J,L,IJL_KVM)/(IDACC(1)*4.)
        END DO
        END DO
        Q(2:IM,JM)=Q(1,JM)
        Q(2:IM,1)=Q(1,1)
        TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
        WRITE (TITLE(40:47),'(A5,I3)') "Level",L
        WRITE (LNAME(40:47),'(A5,I3)') "Level",L
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,2)
      END DO
      
      DO K=1,3 
        L =KCMF(K)
        LNAME="VERT. HEAT DIFF."
        UNITS="cm^2/s"
        SNAME="kvh"//char(l+48)
        Q=UNDEF
        DO J=1,JM
        DO I=1,IMAXJ(J)
          IF (OIJL(I,J,L+1,IJL_MO).gt.0)
     *       Q(I,J)=1d4*.00097**2*OIJL(I,J,L,IJL_KVG)/(IDACC(1)*4.)
        END DO
        END DO
        Q(2:IM,JM)=Q(1,JM)
        Q(2:IM,1)=Q(1,1)
        TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
        WRITE (TITLE(40:47),'(A5,I3)') "Level",L
        WRITE (LNAME(40:47),'(A5,I3)') "Level",L
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,2)
      END DO
C****
C**** Simple scaled OIJ diagnostics
C****
      DO K=1,KOIJ
        byiacc=1./(IDACC(IA_OIJ(K))+teeny)
        adenom=1.
        lname=lname_oij(k)
        k1 = index(lname,' x ')
        if (k1 .gt. 0) then
          if (index(lname,' x PO4') .gt. 0) then
            adenom(1,jm)=0.25
            adenom(1, 1)=0.25
          end if
          lname(k1:50) = ' '
        end if

        Q=UNDEF
        DO J=1,JM
          DO I=1,IMAXJ(J)
            IF (ADENOM(I,J).gt.0 .and. FOCEAN(I,J).gt.0.5)
     *           Q(I,J)=SCALE_OIJ(K)*OIJ(I,J,K)*byiacc/adenom(i,j)
          END DO
        END DO
        Q(2:IM,JM)=Q(1,JM)
        Q(2:IM,1)=Q(1,1)
        TITLE=trim(LNAME)//" ("//trim(UNITS_OIJ(K))//") "
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME_OIJ(K),LNAME_OIJ(K),UNITS_OIJ(K),Q,QJ
     *       ,QSUM,IJGRID_OIJ(K))

      END DO

      END IF
C****
C**** Calculate Horizontal Mass Stream Function and write it
C****
      LNAME= 'HORIZONTAL MASS TRANSPORT STREAMFUNCTION'
      UNITS= 'Sv'
      SNAME= 'osfij'
      FAC   = -2d-9/(IDACC(1)*NDYNO)
      FACST = -1d-9/(IDACC(1)*NDYNO*DTO)
      CALL STRMIJ (OIJL(1,1,1,IJL_MFU),FAC,OLNST(1,1,LN_MFLX),FACST
     *     ,SFIJM)

      TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
      TITLE(51:80)=XLB
      IF (QDIAG) CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,SFIJM,QJ,QSUM
     *     ,IJGRID)
C****
C**** Calculate Salt Stream Function and write it
C****
c      FAC   = -1d-6/(IDACC(1)*NDYNO*DTO)
c      FACST = -.5E-6/(IDACC(1)*NDYNO*DTO)
c      CALL STRMIJ (OIJL(1,1,1,IJL_SFLX),FAC,OLNST(1,1,LN_SFLX),FACST
c     *     ,SFIJS)
C**** Subtract .035 times the Mass Stream Function
c      DO J=1,JM
c        DO I=1,IM
c        IF (SFIJS(I,J).ne.UNDEF) SFIJS(I,J) = SFIJS(I,J)-35.*SFIJM(I,J)
c        END DO
c      END DO
c      TITLE = NAME(2)//'  Run '//XLABEL(1:6)
c      WRITE(TITLE(63:80),'(A6,I4)') JMON0,JYEAR0

C**** Output Key diagnostics: Gulf Stream, ACC, Kuroshio
      WRITE(6,'(A)') " Key horizontal mass stream function diags:"
      GSMAX=0
      GSMIN=100.
      DO J=30,33
        DO I=20,24
          IF (SFIJM(I,J).GT.GSMAX) GSMAX=SFIJM(I,J)
          IF (SFIJM(I,J).LT.GSMIN) GSMIN=SFIJM(I,J)
        END DO
      END DO
      WRITE(6,'(a,F10.3)') " Gulf Stream (MAX in (20-24,30-33)):", GSMAX
     *     -GSMIN
      CKMAX=0
      CKMIN=100.
      DO J=30,33
        DO I=63,67
          IF (SFIJM(I,J).GT.CKMAX) CKMAX=SFIJM(I,J)
          IF (SFIJM(I,J).LT.CKMIN) CKMIN=SFIJM(I,J)
        END DO
      END DO
      WRITE(6,'(a,F10.3)') " Kuroshio    (MAX in (63-67,30-33)):", CKMAX
     *     -CKMIN
      ACMAX=0
      ACMIN=100.
      I=23
      DO J=6,10
        IF (SFIJM(I,J).GT.ACMAX) ACMAX=SFIJM(I,J)
        IF (SFIJM(I,J).LT.ACMIN) ACMIN=SFIJM(I,J)
      END DO
      WRITE(6,'(a,F10.3)') " ACC               (Drakes Passage):", ACMAX
     *     -ACMIN
C****
      RETURN
      END SUBROUTINE OIJOUT

      SUBROUTINE OSFOUT
!@sum  OSFOUT prints out streamfunction diagnostics for ocean
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0
      USE CONSTANT, only : undef
      USE MODEL_COM, only : xlabel,lrunid,jmon0,jyear0,idacc,jdate0
     *     ,amon0,jdate,amon,jyear
      USE OCEAN, only : im,jm,lmo,dxypo,ndyno,dts,dto,imaxj
      USE DAGCOM, only : qdiag,acc_period,zoc1
      USE ODIAG
      IMPLICIT NONE
      REAL*8, DIMENSION(JM-1,0:LMO,0:4) :: SFM,SFS
      REAL*8, DIMENSION(JM+3,LMO+1) :: XJL
      INTEGER I,J,K,L,NS,N,KB,IP1
     *     ,LMSTMIN(LMO),SUMORMN(LMO),JEQ,JDLAT,KXLB
      REAL*8 FAC,FACST
      CHARACTER TITLE*80,lname*50,sname*30,units*50

      TITLE(51:80)=XLB
      LNAME=""
      XJL=undef
C****
C**** Calculate Mass Stream Function and write it
C****
      FAC   = -2d-9/(IDACC(1)*NDYNO)
      FACST = -1d-9/(IDACC(1)*NDYNO*DTO)
      CALL STRMJL (OIJL(1,1,1,IJL_MFU),FAC,OLNST(1,1,LN_MFLX),FACST,SFM)
      TITLE(9:50)=" Mass Stream Function (Sv)"
      DO KB=1,4
        TITLE(1:8)=TRIM(BASIN(KB))
        lname(1:31)=TITLE(1:31)
        sname='sf_'//trim(BASIN(KB))
        units='Sv'
        XJL(2:JM,1:LMO+1)=SFM(1:JM-1,0:LMO,KB)
        IF (QDIAG) CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,2,LMO+1,XJL
     *       ,ZOC1,"Latitude","Depth (m)")
C****
C**** Write data to PRinT file
C****
c        WRITE (6,921) 'Southern Hemisphere',(NINT(FLAT(J,2)),J=1,JEQ)
c        DO L=0,LMO
c          WRITE (6,922) ZOC1(L+1),(NINT(SFM(J,L,KB)),J=1,JEQ)
c        END DO
c        WRITE (6,921) 'Northern Hemisphere',(NINT(FLAT(J,2)),J=JEQ,
c    *     JM-1)
c        DO L=0,LMO
c          WRITE (6,922) ZOC1(L+1),(NINT(SFM(J,L,KB)),J=JEQ,JM-1)
c        END DO
      END DO

C**** Output Key diagnostics: North Atl. overturning
      WRITE(6,'(A46,F6.2)') " North Atlantic overturning: 900m 48N: ",
     *     SFM(35,9,1)
C**** Output Key diagnostics: North Pacific overturning
      WRITE(6,'(A46,F6.2)') " North Pacific overturning:  900m 48N: ",
     *     SFM(35,9,2)
C**** Output Key diagnostics: AABW production
      WRITE(6,'(A46,F6.2)') " Antarctic Bottom Water production:"
     *     //" 3000m 54S: ",SFM(10,12,4)

C****
C**** Calculate Salt Stream Function and write it
C****
c      FAC   = -1d-6/(IDACC(1)*NDYNO*DTO)
c      FACST = -.5E-6/(IDACC(1)*NDYNO*DTO)
c     CALL STRMJL(OIJL(1,1,1,IJL_SFLX),FAC,OLNST(1,1,LN_SFLX),FACST,SFS)
cC**** Subtract .035 times the Mass Stream Function
c      TITLE(9:50)=" Salt Stream Function (Sv)"
c      DO KB=1,4
c        DO L=0,LMO-1
c        DO J=1,JM-1
c          IF (SFS(J,L,KB).ne.UNDEF) SFS(J,L,KB) = SFS(J,L,KB) -
c     *         35.*SFM(J,L,KB)
c        END DO
c        END DO
c        TITLE(1:8)=TRIM(BASIN(KB))
c        lname(1:31)=TITLE(1:31)
c        sname='salt_sf_'//trim(BASIN(KB))
c        units='Sv'
c        XJL(2:JM,1:LMO+1)=SFS(1:JM-1,0:LMO,KB)
c        IF (QDIAG) CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,2,LMO+1,XJL
c     *       ,ZOC1,"Latitude","Depth (m)")
c      END DO
C****
      RETURN
  920 FORMAT ('1',A80/)
  921 FORMAT ('0',A20 / '0',6X,23I5 / 8X,23('-----'))
  922 FORMAT (F6.0,2X,23I5)

      END SUBROUTINE OSFOUT

      SUBROUTINE STROUT
!@sum  STROUT prints out strait diagnostics for ocean
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0
      USE CONSTANT, only : undef,teeny
      USE MODEL_COM, only : xlabel,lrunid,jmon0,jyear0,idacc,jdate0
     *     ,amon0,jdate,amon,jyear
      USE OCEAN, only : im,jm,lmo,focean,dxypo,ndyno,dts,dto
     *     ,imaxj,lmm,ze,dxvo,dypo
      USE DAGCOM, only : acc_period
      USE STRAITS, only : nmst,wist,dist,lmst,name_st
      USE ODIAG
      IMPLICIT NONE
      REAL*8, DIMENSION(LMO,NMST) :: AS
      INTEGER I,J,K,L,NS,N,KB,IP1
     *     ,LMSTMIN(LMO),SUMORMN(LMO),JEQ,JDLAT,KXLB
      REAL*8 GOS,SOS,FAC,FACST,SCALEO(10),TSUM,TEMGS,MQ
      CHARACTER*40, DIMENSION(KOLNST) :: NAME = ""
      CHARACTER :: TITLE*80,lname*50,sname*30,units*50
C****
C**** Strait diagnostics
C****
      TITLE=""
      TITLE(51:80)=XLB
      SUMORMN(:)=1
      NAME(LN_MFLX)='Strait Transport of Mass (10^6 kg/s)'
      SCALEO(LN_MFLX) = 1.D- 6 / DTS
      NAME(LN_GFLX)='Strait Trans of Potential Enthalpy (10^11 W)'
      SCALEO(LN_GFLX) = .5D-11 / DTS
      NAME(LN_SFLX)='Strait Transport of Salt (10^5 kg/s)'
      SCALEO(LN_SFLX) = .5D- 5 / DTS
      NAME(LN_KVM)= 'Vertical Diffusion Coefficient (m^2/s)'
      SCALEO(LN_KVM) = .5
      SUMORMN(LN_KVM)=2
      NAME(LN_ICFL)='Sea Ice Mass Flux (10^6 kg/s)'
      SCALEO(LN_ICFL) = .5D-6 / DTS
      SUMORMN(LN_ICFL)=0
      DO N=1,KOLNST
        if (n.eq.LN_MFLX .or. n.eq.LN_GFLX.or. n.eq.LN_SFLX.or. n.eq
     *       .LN_KVM.or. n.eq.LN_ICFL) THEN
          AS = 0.
          TITLE(1:40) = NAME(N)
          DO NS=1,NMST
            DO L=1,LMST(NS)
              AS(L,NS) = OLNST(L,NS,N)*SCALEO(N)/IDACC(1)
            END DO
          END DO
          IF (N.eq.LN_ICFL) THEN
c            LMSTMIN(:)=1
c            CALL STABLE (LMSTMIN,name_st,AS,TITLE,SUMORMN(N))
          ELSE
            CALL STABLE (LMST,name_st,AS,TITLE,SUMORMN(N))
          END IF
        END IF
      END DO
C****
      RETURN
      END SUBROUTINE STROUT

      SUBROUTINE STABLE (LMST,STRAIT,AX,TITLE,SUMORMN)
!@sum STABLE produces a layer by strait table on the line printer.
C****
C**** Input:
!@var LMST = number of layers for each strait
!@var STRAIT = names of straits
!@var AX = two dimensional input array
!@var TITLE = title of run
!@var SUMORMN flag for whether sum or mean or neither ar printed
C****
      USE OCEAN, only : lmo
      USE STRAITS, only : nmst
      IMPLICIT NONE
      INTEGER, PARAMETER :: LMAXST=6
      INTEGER, INTENT(IN), DIMENSION(NMST) :: LMST
      INTEGER, INTENT(IN) :: SUMORMN
      REAL*8, DIMENSION(LMO,NMST), INTENT(IN) :: AX
      CHARACTER*20, INTENT(IN), DIMENSION(NMST) :: STRAIT
      CHARACTER*80, INTENT(IN) :: TITLE
      INTEGER, SAVE :: LINECT = 0.
      REAL*8 SUML
      INTEGER NS,L,LMAX
C****
C**** Produce leading title lines
C****
      LINECT = LINECT + 19
      IF(LINECT.LE.63)  THEN
        WRITE (6,900) TITLE
      ELSE
        LINECT = 8
        WRITE (6,901) TITLE
      ENDIF
      WRITE (6,902) ('------',L=1,LMAXST),(L,L=1,LMAXST),
     *              ('------',L=1,LMAXST)
C****
C**** Calculate sums and print the table
C****
      DO NS=1,NMST
        LMAX = LMST(NS)
        SUML = 0.
        DO L=1,LMAX
          SUML = SUML + AX(L,NS)
        END DO
        SELECT CASE (SUMORMN)
        CASE (0)                ! no sum or mean
          WRITE (6,'(1X,A20,9X,10I6)') STRAIT(NS),NINT(AX(1:LMAX,NS))
        CASE (1)
          WRITE (6,'(1X,A20,F8.1,1X,10I6)') STRAIT(NS),SUML
     *         ,NINT(AX(1:LMAX,NS))
        CASE (2)
          SUML=SUML/(LMAX-1)
          WRITE (6,'(1X,A20,F8.1,1X,10I6)') STRAIT(NS),SUML
     *         ,NINT(AX(1:LMAX,NS))
        END SELECT
      END DO
C****
      RETURN
  900 FORMAT ('0'//1X,A80)
  901 FORMAT ('1',A80)
  902 FORMAT ('0',29('-'),6A6 / ' Strait',15X,'Sum/Mean',6I6 /
     *        ' ',29('-'),6A6)
  904 FORMAT (1X,A20,F8.1,1X,10I6)
      END SUBROUTINE STABLE

      SUBROUTINE STRMJL (OIJL,FAC,OLNST,FACST,SF)
!@sum  STRMJL calculates the latitude by layer stream function
!@auth G. Russell/G. Schmidt
!@ver  1.0 (from OSFJL064)
C****
C**** Input:
!@var OIJL  = west-east and south-north tracer fluxes (kg/s)
!@var FAC   = global scaling factor
!@var OLNST = strait mass flux (kg/s)
!@var FACST = global scaling factor for straits
C**** Output:
!@var    SF = stream function (kg/s)
C****
      USE CONSTANT, only : undef
      USE OCEAN, only : im,jm,lmo
      USE STRAITS, only : nmst
      USE ODIAG, only : kbasin
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: FAC,FACST
      REAL*8, INTENT(IN), DIMENSION(IM,JM,LMO,2) :: OIJL
      REAL*8, INTENT(IN), DIMENSION(LMO,NMST) :: OLNST
      REAL*8, INTENT(OUT), DIMENSION(JM-1,0:LMO,0:4) :: SF
      REAL*8, DIMENSION(4) :: SUMB
      INTEGER :: KB,I,J,L,K
C****
C**** Zero out the Stream Function at the ocean bottom
C****
      SF(:,:,:) = 0.
C****
C**** Integrate the Stream Function upwards from the ocean bottom
C****
      DO L=LMO-1,0,-1
      DO J=1,JM-1
        DO KB=1,4
          SF(J,L,KB) = SF(J,L+1,KB)
        END DO
        DO I=1,IM
          KB = KBASIN(I,J+1)
          SF(J,L,KB) = SF(J,L,KB) + OIJL(I,J,L+1,2)*FAC
        END DO
      END DO
C****
C**** Add south-north flow between basins
C****
C**** Add Bering Strait flow into Atlantic from Pacific
C      DO 310 J=1,38
C      SF(J,L,1) = SF(J,L,1) + OIJL(3,39,L+1,2)*FAC
C  310 SF(J,L,2) = SF(J,L,2) - OIJL(3,39,L+1,2)*FAC
C**** Add Indonesian Straits flow into Pacific from Indian
C      DO 320 J=1,20
C      SF(J,L,2) = SF(J,L,2) + OIJL(59,21,L+1,2)*FAC
C      SF(J,L,2) = SF(J,L,2) + OIJL(60,21,L+1,2)*FAC
C      SF(J,L,3) = SF(J,L,3) - OIJL(59,21,L+1,2)*FAC
C  320 SF(J,L,3) = SF(J,L,3) - OIJL(60,21,L+1,2)*FAC
C      DO 330 J=1,19
C      SF(J,L,2) = SF(J,L,2) + OIJL(62,20,L+1,2)*FAC
C  330 SF(J,L,3) = SF(J,L,3) - OIJL(62,20,L+1,2)*FAC
C****
C**** Add west-east flow between basins
C****
C**** Add Drake Passage flow into Atlantic from Pacific
C      FLOW = 0.
C      DO 410 J=8,1,-1
C      FLOW = FLOW + OIJL(23,J+1,L+1,1)*FAC
C      SF(J,L,1) = SF(J,L,1) + FLOW
C  410 SF(J,L,2) = SF(J,L,2) - FLOW
C**** Add flow south of Africa into Indian from Atlantic
C      FLOW = 0.
C      DO 420 J=13,1,-1
C      FLOW = FLOW + OIJL(40,J+1,L+1,1)*FAC
C      SF(J,L,3) = SF(J,L,3) + FLOW
C  420 SF(J,L,1) = SF(J,L,1) - FLOW
C**** Add flow south of Australia into Pacific from Indian
C      FLOW = 0.
C      DO 430 J=12,1,-1
C      FLOW = FLOW + OIJL(65,J+1,L+1,1)*FAC
C      SF(J,L,2) = SF(J,L,2) + FLOW
C  430 SF(J,L,3) = SF(J,L,3) - FLOW
C****
C**** Add strait flow from to the Stream Function
C****
C**** Fury & Hecla: (19,42) to (20,40)
      SF(41,L,1) = SF(41,L,1) - OLNST(L+1,1)*FACST
      SF(40,L,1) = SF(40,L,1) - OLNST(L+1,1)*FACST
C**** Nares: (22,43) to (24,44)
      SF(43,L,1) = SF(43,L,1) + OLNST(L+1,2)*FACST
C**** Gibralter: (35,32) to (37,33)
      SF(32,L,1) = SF(32,L,1) + OLNST(L+1,3)*FACST
C**** English: (36,36) to (37,37)
      SF(36,L,1) = SF(36,L,1) + OLNST(L+1,4)*FACST
C**** Bosporous: (42,33) to (43,34)
      SF(33,L,1) = SF(33,L,1) + OLNST(L+1,6)*FACST
C**** Red Sea: (44,29) to (45,28)
      SF(28,L,3) = SF(28,L,3) - OLNST(L+1,7)*FACST
C**** Bab-al-Mandab: (45,28) to (46,27)
      SF(27,L,3) = SF(27,L,3) - OLNST(L+1,8)*FACST
C**** Hormuz: (47,30) to (49,29)
      SF(29,L,3) = SF(29,L,3) - OLNST(L+1,9)*FACST
C**** Korea: (62,32) to (63,33)
      SF(32,L,2) = SF(32,L,2) + OLNST(L+1,11)*FACST
C**** Soya: (64,34) to (65,35)
      SF(34,L,2) = SF(34,L,2) + OLNST(L+1,12)*FACST
C**** Malacca: (56,25) to (58,24), from Indian Ocean to Pacific Ocean
      SF(24,L,2) = SF(24,L,2) - OLNST(L+1,10)*FACST  !*.5
c      DO 510 J=1,23
c      SF(J,L,2) = SF(J,L,2) - OLNST(L+1,10)*FACST*.5
c  510 SF(J,L,3) = SF(J,L,3) + OLNST(L+1,10)*FACST*.5
      END DO
C****
C**** Calculate global Stream Function by summing it over 3 oceans
C****
      SF(:,:,4) = SF(:,:,1) + SF(:,:,2) + SF(:,:,3) + SF(:,:,4)
C****
C**** Mask streamfunction so that topography is put to skip value
C****
      DO J=1,JM
        SUMB = 0
        DO I=1,IM
          IF (KBASIN(I,J).gt.0) SUMB(KBASIN(I,J))=1.
        END DO
        SUMB(4)=SUMB(1)+SUMB(2)+SUMB(3)+SUMB(4)
        DO K=1,4
          IF (SUMB(K).eq.0) THEN
            SF(J,0:LMO,K) = UNDEF   ! UNDEF areas of no integral
          ELSE
            L=0
  620       L=L+1       ! UNDEF below bottom values
            IF (L.LE.LMO) THEN
              IF (SF(J,L,K).eq.0 .and. SF(J,L-1,K).eq.0) THEN
                SF(J,L:LMO,K) = UNDEF
              ELSE
                GOTO 620
              END IF
            END IF
          END IF
        END DO
      END DO
C****
      RETURN
      END

      SUBROUTINE OBASIN
!@sum  OBASIN Read in KBASIN: 0=continent,1=Atlantic,2=Pacific,3=Indian
!@auth G. Russell
!@ver  1.0
      USE OCEAN, only : IM,JM,focean
      USE ODIAG, only : kbasin
      USE FILEMANAGER
      IMPLICIT NONE
      CHARACTER TITLE*72, CBASIN(IM,JM)
      INTEGER J,I,iu_KB
C****
C**** read in basin data
      call openunit('KBASIN',iu_KB,.false.,.true.)

      READ  (iu_KB,900) TITLE
      WRITE (6,*) 'Read on unit ',iu_KB,': ',TITLE
      READ  (iu_KB,900)
      DO J=JM,1,-1
        READ  (iu_KB,901) (CBASIN(I,J),I=1,72)
      END DO
      call closeunit(iu_KB)

      DO J=1,JM
      DO I=1,IM
        SELECT CASE (CBASIN(I,J))
        CASE DEFAULT
          KBASIN(I,J) = 0
          IF (FOCEAN(I,J).gt.0) WRITE(6,*)
     *         "Warning: Ocean box not defined in KBASIN ",i,j
        CASE ('A')
          KBASIN(I,J) = 1
        CASE ('P')
          KBASIN(I,J) = 2
        CASE ('I')
          KBASIN(I,J) = 3
        CASE ('G')
          KBASIN(I,J) = 4
        END SELECT
      END DO
      END DO
C****
      RETURN
  900 FORMAT (A72)
  901 FORMAT (72A1)
      END SUBROUTINE OBASIN

      SUBROUTINE STRMIJ (OIJL,FAC,OLNST,FACST,SF)
!@sum STRMIJ calculates the latitude by longitude stream function for a
!@+   given quantity.
C****
C**** Input:
!@var OIJL  = west-east and south-north tracer fluxes (kg/s)
!@var   FAC = global scaling factor
!@var OLNST = strait mass flux (kg/s)
!@var FACST = global scaling factor for straits
C**** Output:
!@var    SF = stream function (kg/s)
C****
      USE OCEAN, only : im,jm,lmo,lmm
      USE STRAITS, only : nmst,lmst
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: FAC,FACST
      REAL*8, INTENT(IN), DIMENSION(IM,JM,LMO) :: OIJL
      REAL*8, INTENT(IN), DIMENSION(LMO,NMST) :: OLNST
      REAL*8, INTENT(OUT), DIMENSION(IM,JM) :: SF
      REAL*8, DIMENSION(4) :: SUMB
      INTEGER :: I,J,L
      REAL*8 TSUM
C****
C**** Integrate up from South Pole
C****
      SF=0
      DO J=2,JM-1
      DO I=1,IM
        SF(I,J) = SF(I,J-1)
        DO L=1,LMM(I,J)
          SF(I,J) = SF(I,J) + OIJL(I,J,L)*FAC
        END DO
      END DO
C****
C**** Add strait flow from to the Stream Function
C****
C**** Fury & Hecla: (19,42) to (20,40)
      IF (J.eq.41) SF(19,41) = SF(19,41) + SUM(OLNST(1:LMST(1),1))*FACST
C**** Nares: (22,43) to (24,44)
      IF (J.eq.44) SF(22,44) = SF(22,44) + SUM(OLNST(1:LMST(2),2))*FACST
      IF (J.eq.44) SF(23,44) = SF(23,44) + SUM(OLNST(1:LMST(2),2))*FACST
C**** Gibrater: (35,32) to (37,33)
      IF (J.eq.33) SF(35,33) = SF(35,33) + SUM(OLNST(1:LMST(3),3))*FACST
      IF (J.eq.33) SF(36,33) = SF(36,33) + SUM(OLNST(1:LMST(3),3))*FACST
C**** Engish: (36,36) to (37,37)
      IF (J.eq.37) SF(36,37) = SF(36,37) + SUM(OLNST(1:LMST(4),4))*FACST
C**** Bosporous: (42,33) to (43,34)
      IF (J.eq.34) SF(42,34) = SF(42,34) + SUM(OLNST(1:LMST(6),6))*FACST
C**** Red Sea: (44,29) to (45,28)
      IF (J.eq.29) SF(44,29) = SF(44,29) + SUM(OLNST(1:LMST(7),7))*FACST
C**** Bab-al-Mandab: (45,28) to (46,27)
      IF (J.eq.28) SF(45,28) = SF(45,28) + SUM(OLNST(1:LMST(8),8))*FACST
C**** Hormuz: (47,30) to (49,29)
      IF (J.eq.30) SF(47,30) = SF(47,30) + SUM(OLNST(1:LMST(9),9))*FACST
      IF (J.eq.30) SF(48,30) = SF(48,30) + SUM(OLNST(1:LMST(9),9))*FACST
C**** Korea: (62,32) to (63,33)
      IF (J.eq.33) SF(62,33) = SF(62,33)+SUM(OLNST(1:LMST(11),11))*FACST
C**** Soya: (64,34) to (65,35)
      IF (J.eq.35) SF(64,35) = SF(64,35)+SUM(OLNST(1:LMST(12),12))*FACST
C**** Malacca: (56,25) to (58,24),
      IF (J.eq.25) SF(56,25) = SF(56,25)+SUM(OLNST(1:LMST(10),10))*FACST
      IF (J.eq.25) SF(57,25) = SF(57,25)+SUM(OLNST(1:LMST(10),10))*FACST
      END DO
C**** Correct SF for mean E-W drift (SF over topography --> 0)
      TSUM=0
      DO I=14,17   ! (mid N. America as example)
        TSUM=TSUM+SF(I,34)
      END DO
      TSUM=TSUM/4.
      DO J=1,JM-1
        DO I=1,IM
          SF(I,J)=SF(I,J)-TSUM
        END DO
      END DO
      RETURN
      END

      SUBROUTINE OJLOUT
!@sum OJLOUT calculates basin means, lat. and long. sections
!@+   and advective tracer fluxes
!@auth Gavin Schmidt/Gary Russell
      USE CONSTANT, only : undef,teeny
      USE MODEL_COM, only : idacc
#ifdef TRACERS_OCEAN
      USE TRACER_COM, only : ntm,trw0,trname,ntrocn
#endif
      USE OCEAN, only : im,jm,lmo,ze,imaxj,focean,ndyno,dypo,dts,dxvo
     *     ,dxypo
      USE DAGCOM, only : qdiag,zoc
      USE ODIAG
#ifdef TRACERS_OCEAN
      USE TRACER_DIAG_COM, only : to_per_mil
#endif
      IMPLICIT NONE
      INTEGER, PARAMETER :: NSEC=3, NBAS=4
      REAL*8, DIMENSION(JM+3,LMO+1) :: XJL
      REAL*8 XB0(JM,LMO,0:NBAS),X0(IM,LMO),XS(IM,LMO),
     *     XBG(JM,LMO,0:NBAS),XBS(JM,LMO,0:NBAS),XG(IM,LMO)
#ifdef TRACERS_OCEAN
      REAL*8 XBT(JM,LMO,0:NBAS,NTM),XT(IM,LMO,NTM)
#endif
      CHARACTER TITLE*80,EW*1,NS*1,LNAME*50,SNAME*50,UNITS*50
      CHARACTER LABI*16,LABJ*16
      character*50 :: unit_string
C**** these parameters control the output (Resolution dependent!)
C****  latitudinal sections: at ILON from 1 to JM
C****      Defaults:  3 (Pacific), 30 (Atlantic), 49 (Indian)
C****  longitudinal sections: at JLAT starting at I1 (with wrap around)
C****      Defaults:  7 (South. Oc), 23 (Eq), 35 (sub-trop gyre)
C****      overlaps across the dateline to avoid splitting Pacific
      INTEGER*4 :: JLAT(NSEC) = (/ 7,23,33/), ILON(NSEC) = (/ 3,30,49/),
     *     I1=41, KITR(3)
      INTEGER I,J,L,KB,ISEC,II,ILAT,JLON,N
      REAL*8 GOS,SOS,TEMGS,ASUM(IM),GSUM,ZONAL(LMO)

      IF (QDIAG) then
      XJL=0.
C**** Calculate basin and zonal averages of ocean mass
      XB0 = 0.
      DO L=1,LMO
        DO J= 1,JM
          DO I=1,IMAXJ(J)
            KB=KBASIN(I,J)
            XB0(J,L,KB) = XB0(J,L,KB) + OIJL(I,J,L,IJL_MO)
          END DO
        END DO
      END DO
      XB0(:,:,4) = XB0(:,:,1) + XB0(:,:,2) + XB0(:,:,3)
      TITLE(51:80) = XLB
C****
C**** Quantites that are calculated
C****
C**** Basin and Zonal averages of tracers
C****
      XBG = 0.  ;  XBS = 0.
#ifdef TRACERS_OCEAN
      XBT = 0.
#endif
      DO L=1,LMO
        DO J= 1,JM
          DO I=1,IMAXJ(J)
            KB=KBASIN(I,J)
            IF(FOCEAN(I,J).gt..5)  THEN
              XBG(J,L,KB) = XBG(J,L,KB) + OIJL(I,J,L,IJL_G0M)
              XBS(J,L,KB) = XBS(J,L,KB) + OIJL(I,J,L,IJL_S0M)
#ifdef TRACERS_OCEAN
              XBT(J,L,KB,:) = XBT(J,L,KB,:) + TOIJL(I,J,L,TOIJL_CONC,:)
#endif
            END IF
          END DO
        END DO
      END DO
      XBG(:,:,4) = XBG(:,:,1) + XBG(:,:,2) + XBG(:,:,3)
      XBS(:,:,4) = XBS(:,:,1) + XBS(:,:,2) + XBS(:,:,3)
#ifdef TRACERS_OCEAN
      XBT(:,:,4,:)= XBT(:,:,1,:)+XBT(:,:,2,:)+XBT(:,:,3,:)
#endif
      DO KB=1,4
        DO J=1,JM
          DO L=1,LMO
            IF (XB0(J,L,KB).ne.0.) THEN
#ifdef TRACERS_OCEAN
              do n=1,ntm
              if (to_per_mil(n).gt.0) then
                XBT(j,l,kb,n)= 1d3*(XBT(j,l,kb,n)/
     *               ((XB0(J,L,KB)*DXYPO(J)-XBS(J,L,KB))*trw0(n))-1.)
              else
                XBT(j,l,kb,n)= 10.**(-ntrocn(n))*XBT(j,l,kb,n)/
     *               (XB0(J,L,KB)*DXYPO(J))
              end if
              end do
#endif
              GOS = XBG(J,L,KB)/(XB0(J,L,KB)*DXYPO(J))
              SOS = XBS(J,L,KB)/(XB0(J,L,KB)*DXYPO(J))
              XBG(J,L,KB) = TEMGS(GOS,SOS)
              XBS(J,L,KB) = 1d3*SOS
            ELSE
              XBG(J,L,KB)=undef
              XBS(J,L,KB)=undef
#ifdef TRACERS_OCEAN
              XBT(j,l,kb,:)= undef
#endif
            END IF
          END DO
        END DO
        if (KB.le.3) then
          TITLE(21:50)=" "//TRIM(BASIN(KB))//" Basin Average"
          SNAME(5:30)="_"//BASIN(KB)(1:3)
        else
          TITLE(21:50)=" Zonal Average"
          SNAME(5:30)="_Zonal"
        end if
        LNAME(21:50)=TITLE(21:50)
        TITLE(1:20)="Temperature (C)"
        LNAME(1:20)="Temperature"
        SNAME(1:4)="Temp"
        UNITS="C"
        XJL(1:JM,1:LMO)=XBG(1:JM,1:LMO,KB)
        CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,1,LMO,XJL,ZOC
     *       ,"Latitude","Depth (m)")
        TITLE(1:20)="Salinity  (ppt)"
        LNAME(1:20)="Salinity"
        SNAME(1:4)="Salt"
        UNITS="ppt"
        XJL(1:JM,1:LMO)=XBS(1:JM,1:LMO,KB)
        CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,1,LMO,XJL,ZOC
     *       ,"Latitude","Depth (m)")
#ifdef TRACERS_OCEAN
        DO N=1,NTM
          if (to_per_mil(n).gt.0) then 
            UNITS="permil"
          else
            UNITS=unit_string(ntrocn(n),'kg/kg')
          end if
          TITLE(1:20)=trim(trname(n))//" ("//trim(UNITS)//")"
          LNAME(1:20)=trim(trname(n))
          SNAME(1:4)=trname(n)(1:4)
          XJL(1:JM,1:LMO)=XBT(1:JM,1:LMO,KB,N)
          CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,1,LMO,XJL,ZOC
     *         ,"Latitude","Depth (m)")
        END DO
#endif
      END DO
C****
C**** Longitudinal sections
C****
      DO ISEC=1,NSEC
        XB0 = undef  ; XBS = undef ;  XBG = undef
#ifdef TRACERS_OCEAN
        XBT = undef
#endif
        IF (ILON(ISEC).ne.-1) THEN
          I=ILON(ISEC)
          IF (I.le.IM/2) THEN
            EW="W"
          ELSE
            EW="E"
          END IF
          DO L=1,LMO
            DO J= 1,JM
              IF (OIJL(I,J,L,IJL_MO).ne.0) THEN
                GOS = OIJL(I,J,L,IJL_G0M)/(OIJL(I,J,L,IJL_MO)*DXYPO(J))
                SOS = OIJL(I,J,L,IJL_S0M)/(OIJL(I,J,L,IJL_MO)*DXYPO(J))
                XBG(J,L,1)= TEMGS(GOS,SOS)
                XBS(J,L,1)= 1d3*SOS
#ifdef TRACERS_OCEAN
                do n=1,ntm
                  if (to_per_mil(n).gt.0) then
                    XBT(j,l,1,n)=1d3*(TOIJL(I,J,L,TOIJL_CONC,n)/
     *              ((OIJL(I,J,L,IJL_MO)*DXYPO(J)-OIJL(I,J,L,IJL_S0M))
     *                   *trw0(n))-1.) 
                  else
                    XBT(j,l,1,n)= 10.**(-ntrocn(n))*TOIJL(I,J,L
     *                   ,TOIJL_CONC,n)/(OIJL(I,J,L,IJL_MO)*DXYPO(J))
                  end if
                end do
#endif
              ELSE
                XBG(J,L,1)=undef ; XBS(J,L,1)=undef
#ifdef TRACERS_OCEAN
                XBT(J,L,1,:)=undef
#endif
              END IF
            END DO
          END DO
          ILAT=ABS(NINT(-182.5+ILON(ISEC)*5.))
          WRITE(LABI,'(I3,A1)') ILAT,EW
          TITLE(1:50)="Temperature Section         (C)"
          WRITE(TITLE(21:25),'(I3,A1)') ILAT,EW
          LNAME(1:25)=TITLE(1:25)
          SNAME="Temp_"//adjustl(labi)
          UNITS="C"
          XJL(1:JM,1:LMO)=XBG(1:JM,1:LMO,1)
          CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,1,LMO,XJL,ZOC
     *         ,"Latitude","Depth (m)")
          TITLE(1:50)="Salinity Section            (ppt)"
          WRITE(TITLE(21:25),'(I3,A1)') ILAT,EW
          LNAME(1:25)=TITLE(1:25)
          SNAME="Salt_"//adjustl(labi)
          UNITS="ppt"
          XJL(1:JM,1:LMO)=XBS(1:JM,1:LMO,1)
          CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,1,LMO,XJL,ZOC
     *         ,"Latitude","Depth (m)")
#ifdef TRACERS_OCEAN
          DO N=1,NTM
          if (to_per_mil(n).gt.0) then 
            UNITS="permil"
          else
            UNITS=unit_string(ntrocn(n),'kg/kg')
          end if
          TITLE(1:50)=trname(n)//" Section          ("//TRIM(UNITS)//")"
          WRITE(TITLE(21:25),'(I3,A1)') ILAT,EW
          LNAME(1:25)=TITLE(1:25)
          SNAME=trim(trname(n))//"_"//adjustl(labi)
          XJL(1:JM,1:LMO)=XBT(1:JM,1:LMO,1,N)
          CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,1,LMO,XJL,ZOC
     *         ,"Latitude","Depth (m)")
          END DO
#endif

C**** EW fluxes of mass, heat and salt
C**** Note that ocean fluxes are defined 1:JM-1 and need to be moved up
C**** one place for compatibility with AGCM output format routines
          DO L=1,LMO
            DO J= 1,JM-1
C**** GM fluxes are also saved, so add the GM heat and salt fluxes here
              IF (OIJL(I,J,L,IJL_MFU).ne.0) THEN
                XB0(J,L,1)= 2d0*OIJL(I,J,L,IJL_MFU)/(IDACC(1)*NDYNO
     *               *DYPO(J) *(ZE(L)-ZE(L-1)))
                XBG(J,L,1)= 1d-6*(OIJL(I,J,L,IJL_GFLX)+OIJL(I,J,L
     *               ,IJL_GGMFL))/(IDACC(1)*DTS*DYPO(J)*(ZE(L)-ZE(L-1)))
                XBS(J,L,1)= 1d-6*(OIJL(I,J,L,IJL_MFU)+OIJL(I,J,L
     *               ,IJL_SGMFL))/(IDACC(1)*DTS*DYPO(J)*(ZE(L)-ZE(L-1)))
              ELSE
                XB0(J,L,1)=undef
                XBG(J,L,1)=undef
                XBS(J,L,1)=undef
              END IF
            END DO
          END DO
          ILAT=ABS(NINT(-180.+ILON(ISEC)*5.))
          WRITE(LABI,'(I3,A1)') ILAT,EW
          TITLE(1:50)=" EW Mass Flux            (kg/s m^2)"
          WRITE(TITLE(21:25),'(I3,A1)') ILAT,EW
          LNAME(1:25)=TITLE(1:25)
          SNAME="EWmflx_"//adjustl(labi)
          UNITS="kg/s m^2"
          XJL(2:JM,1:LMO)=XB0(1:JM-1,1:LMO,1)
          CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,2,LMO,XJL,ZOC
     *         ,"Latitude","Depth (m)")
          TITLE(1:50)=" EW Heat Flux            (10^6 J/s m^2)"
          WRITE(TITLE(21:25),'(I3,A1)') ILAT,EW
          LNAME(1:25)=TITLE(1:25)
          SNAME="EWhflx_"//adjustl(labi)
          UNITS="10^6 J/s m^2"
          XJL(2:JM,1:LMO)=XBG(1:JM-1,1:LMO,1)
          CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,2,LMO,XJL,ZOC
     *         ,"Latitude","Depth (m)")
          TITLE(1:50)=" EW Salt Flux            (10^6 kg/s m^2)"
          WRITE(TITLE(21:25),'(I3,A1)') ILAT,EW
          LNAME(1:25)=TITLE(1:25)
          SNAME="EWsflx_"//adjustl(labi)
          UNITS="10^6 kg/s m^2"
          XJL(2:JM,1:LMO)=XBS(1:JM-1,1:LMO,1)
          CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,2,LMO,XJL,ZOC
     *         ,"Latitude","Depth (m)")
        END IF
      END DO
C****
C**** Latitudinal sections
C****
      ASUM=0. ; GSUM=0. ; ZONAL=0.
      DO ISEC=1,NSEC
        X0=undef ; XS=undef ; XG=undef
#ifdef TRACERS_OCEAN
        XT=undef
#endif
        IF (JLAT(ISEC).ne.-1) THEN
          IF (JLAT(ISEC).le.23) THEN
            NS="S"
          ELSE
            NS="N"
          END IF
          J=JLAT(ISEC)
          DO L=1,LMO
            DO II= 1,IM
              I = II+I1-1
              IF (I.gt.IM) I=I-IM
              IF (OIJL(I,J,L,IJL_MO).ne.0) THEN
                GOS = OIJL(I,J,L,IJL_G0M)/(OIJL(I,J,L,IJL_MO)*DXYPO(J))
                SOS = OIJL(I,J,L,IJL_S0M)/(OIJL(I,J,L,IJL_MO)*DXYPO(J))
                XG(II,L)= TEMGS(GOS,SOS)
                XS(II,L)= 1d3*SOS
#ifdef TRACERS_OCEAN
                do n=1,ntm
                  if (to_per_mil(n).gt.0) then
                    XT(II,l,n)= 1d3*(TOIJL(I,J,L,TOIJL_CONC,n)/
     *              ((OIJL(I,J,L,IJL_MO)*DXYPO(J)-OIJL(I,J,L,IJL_S0M))
     *                   *trw0(n))-1.)
                  else
                    XT(II,l,n)= 10.**(-ntrocn(n))*TOIJL(I,J,L,TOIJL_CONC
     *                   ,n)/(OIJL(I,J,L,IJL_MO)*DXYPO(J))
                  end if
                end do
#endif
              ELSE
                XG(II,L)=undef ; XS(II,L)=undef
#ifdef TRACERS_OCEAN
                XT(II,L,:)=undef
#endif
              END IF
            END DO
          END DO
          JLON=ABS(NINT(-94.+JLAT(ISEC)*4.))
          WRITE(LABJ,'(I3,A1)') JLON,NS
          TITLE(1:50)="Temperature Section              (C)"
          WRITE(TITLE(21:25),'(I3,A1)') JLON,NS
          LNAME(1:25)=TITLE(1:25)
          SNAME="Temp_"//adjustl(labj)
          UNITS="C"
          CALL POUT_IL(TITLE,sname,lname,units,I1,1,LMO,XG
     *         ,ZOC,"Longitude","Depth (m)",ASUM,GSUM,ZONAL)
          TITLE(1:50)="Salinity Section                 (ppt)"
          WRITE(TITLE(21:25),'(I3,A1)') JLON,NS
          LNAME(1:25)=TITLE(1:25)
          SNAME="Salt_"//adjustl(labj)
          UNITS="ppt"
          CALL POUT_IL(TITLE,sname,lname,units,I1,1,LMO,XS
     *         ,ZOC,"Longitude","Depth (m)",ASUM,GSUM,ZONAL)
#ifdef TRACERS_OCEAN
          DO N=1,NTM
          if (to_per_mil(n).gt.0) then 
            UNITS="permil"
          else
            UNITS=unit_string(ntrocn(n),'kg/kg')
          end if
          TITLE(1:50)=trname(n)//" Section          ("//TRIM(UNITS)//")"
          WRITE(TITLE(21:25),'(I3,A1)') JLON,NS
          LNAME(1:25)=TITLE(1:25)
          SNAME=trim(trname(n))//"_"//adjustl(labj)
          CALL POUT_IL(TITLE,sname,lname,units,I1,1,LMO,XT(1,1,N)
     *         ,ZOC,"Longitude","Depth (m)",ASUM,GSUM,ZONAL)
          END DO
#endif
C**** Fluxes
          DO L=1,LMO
            DO II= 1,IM
              I = II+I1-1
              IF (I.gt.IM) I=I-IM
              IF (OIJL(I,J,L,IJL_MFV).ne.0) THEN
                X0(II,L)= 2d0*OIJL(I,J,L,IJL_MFV)/(IDACC(1)*NDYNO
     *               *DXVO(J) *(ZE(L)-ZE(L-1)))
                XG(II,L)= 1d-6*(OIJL(I,J,L,IJL_GFLX+1)+OIJL(I,J,L
     *               ,IJL_GGMFL+1))/(IDACC(1)*DTS*DXVO(J)*(ZE(L)-ZE(L-1)
     *               ))
                XS(II,L)= 1d-6*(OIJL(I,J,L,IJL_SFLX+1)+OIJL(I,J,L
     *               ,IJL_SGMFL+1))/(IDACC(1)*DTS*DXVO(J)*(ZE(L)-ZE(L-1)
     *               ))
              ELSE
                X0(II,L)=undef ; XG(II,L)=undef ; XS(II,L)=undef
              END IF
            END DO
          END DO
          JLON=ABS(NINT(-92.+JLAT(ISEC)*4.))
          WRITE(LABJ,'(I3,A1)') JLON,NS
          TITLE(1:50)=" NS Mass Flux            (kg/s m^2)"
          WRITE(TITLE(21:25),'(I3,A1)') JLON,NS
          LNAME(1:25)=TITLE(1:25)
          SNAME="NSmflx_"//adjustl(labj)
          UNITS="kg/s m^2"
          CALL POUT_IL(TITLE,sname,lname,units,I1,2,LMO,X0
     *         ,ZOC,"Longitude","Depth (m)",ASUM,GSUM,ZONAL)
          TITLE(1:50)=" NS Heat Flux            (10^6  J/s m^2)"
          WRITE(TITLE(21:25),'(I3,A1)') JLON,NS
          LNAME(1:25)=TITLE(1:25)
          SNAME="NShflx_"//adjustl(labj)
          UNITS="10^6 W/m^2"
          CALL POUT_IL(TITLE,sname,lname,units,I1,2,LMO,XG
     *         ,ZOC,"Longitude","Depth (m)",ASUM,GSUM,ZONAL)
          TITLE(1:50)=" NS Salt Flux            (10^6 kg/s m^2)"
          WRITE(TITLE(21:25),'(I3,A1)') JLON,NS
          LNAME(1:25)=TITLE(1:25)
          SNAME="NSsflx_"//adjustl(labj)
          UNITS="10^6 kg/s m^2"
          CALL POUT_IL(TITLE,sname,lname,units,I1,2,LMO,XS
     *         ,ZOC,"Longitude","Depth (m)",ASUM,GSUM,ZONAL)
        END IF
      END DO
      END IF
C****
      RETURN
      END SUBROUTINE OJLOUT

      SUBROUTINE OTJOUT
!@sum OTJOUT calculate vertically integrated basin and zonal
!@+   northward transports
!@auth Gavin Schmidt/Gary Russell
      USE MODEL_COM, only : idacc
      USE OCEAN, only : im,jm,lmo,ndyno,dts
      USE STRAITS, only : nmst,wist,dist,lmst,name_st
      USE DAGCOM, only : qdiag
      USE ODIAG
      IMPLICIT NONE
      CHARACTER TITLE*80, NAME(7)*50,YAXIS(3)*8
      REAL*8 VLAT(0:JM),X(0:JM,4,6),SCALEM(3),SCALES(3),SOLNST(NMST),
     *       MV(4),MT(4), XSPEC(3),YSPEC(3,3)
      INTEGER, PARAMETER :: INC=1+(JM-1)/24
      INTEGER IS(4)
      DATA NAME/
     1  'Northward Transport of Mass (10^9 kg/sec)',
     2  'North. Trans. of Potential Enthalpy (10^15 W)',
     3  'North. Trans. of Salt - .035*Mass (10^6 kg/sec)',
     4  'North. Trans. of Heat in Atlantic Ocean (10^15 W)',
     5  'North. Trans. of Heat in Pacific Ocean (10^15 W)',
     6  'North. Trans. of Heat in Indian Ocean (10^15 W)',
     7  'North. Trans. of Heat in Global Ocean (10^15 W)'/
      DATA YAXIS /'  Mass  ', 'Enthalpy', '  Salt  '/
      DATA XSPEC /-90.,90.,30./
      DATA YSPEC /-4.,4.,1., -5.,6.,1., -50.,50.,10./
      INTEGER I,J,L,K,KQ,KB,NOIJL,NOIJLGM,NOLNST,NS,IDLAT
      REAL*8 SOIJL,FJEQ,SOIJLGM
C****
C****  N          Contents of OIJL(I,J,L,N)       Scaling factor
C****  -          -------------------------       --------------
C**** IJL_MV      South-north mass flux (kg/s)    2/IDACC(1)*NDYNO
C**** IJL_GFLX+1  South-north heat flux (J/s)     1/IDACC(1)*DTS
C**** IJL_SFLX+1  South-north salt flux (kg/s)    1/IDACC(1)*DTS
C****
C****  N       Contents of OLNST(L,N)             Scaling factor
C****  -       ----------------------             --------------
C**** LN_MFLX  Strait mass flux (kg/s)            1/IDACC(1)*DTS
C**** LN_GFLX  Strait heat flux (J/s)            .5/IDACC(1)*DTS
C**** LN_SFLX  Strait salt flux (kg/s)           .5/IDACC(1)*DTS
C****
      SCALEM(1) = 2.D- 9/(IDACC(1)*NDYNO)
      SCALEM(2) = 1.D-15/(IDACC(1)*DTS)
      SCALEM(3) = 1.D- 6/(IDACC(1)*DTS)
      SCALES(1) = 1.D- 9/(IDACC(1)*DTS)
      SCALES(2) = .5D-15/(IDACC(1)*DTS)
      SCALES(3) = .5D- 6/(IDACC(1)*DTS)
C****
      X = 0.
C****
C**** Loop over quantites
C****
      DO 620 KQ=1,3
C****
C**** Calculate the transport entering each domain from the
C**** southern boundary
C****
      SELECT CASE (KQ)
      CASE (1)
        NOIJL=IJL_MFV
        NOIJLGM=0
      CASE (2)
        NOIJL=IJL_GFLX+1
        NOIJLGM=IJL_GGMFL+1
      CASE (3)
        NOIJL=IJL_SFLX+1
        NOIJLGM=IJL_SGMFL+1
      END SELECT

      DO 120 J=1,JM-1
      DO 120 I=1,IM
      IF(OIJL(I,J,1,IJL_MFV).eq.0)  GO TO 120
      KB = KBASIN(I,J+1)
      SOIJL = 0.
      SOIJLGM = 0.
      DO L=1,LMO
        SOIJL = SOIJL + OIJL(I,J,L,NOIJL)
        IF (NOIJLGM.gt.0) SOIJLGM = SOIJLGM + OIJL(I,J,L,NOIJLGM)
      END DO
      X(J,KB,KQ) = X(J,KB,KQ) + (SOIJL+SOIJLGM)*SCALEM(KQ)
      X(J, 4,KQ) = X(J, 4,KQ) + (SOIJL+SOIJLGM)*SCALEM(KQ)
      IF (KQ.eq.2) THEN
        X(J,KB,5) = X(J,KB,5) + SOIJLGM*SCALEM(KQ)
        X(J, 4,5) = X(J, 4,5) + SOIJLGM*SCALEM(KQ)
      END IF
  120 CONTINUE
      IF(KQ.ne.2)  GO TO 200
C****
C**** Accumulate northward transport of heat by overturning and by
C**** horizontal gyre
C****
      DO 190 J=1,JM-1
      DO 180 L=1,LMO
        IS=0 ; MV=0 ; MT=0
      DO 170 I=1,IM
      IF(OIJL(I,J,L,IJL_MFV).eq.0)  GO TO 170
      KB = KBASIN(I,J+1)
      IS(KB) = IS(KB) + 1
      MV(KB) = MV(KB) + OIJL(I,J,L,IJL_MFV)
      MT(KB) = MT(KB) + OIJL(I,J,L,IJL_GFLX+1)/OIJL(I,J,L,IJL_MFV)
  170 CONTINUE
      IS(4) = IS(1) + IS(2) + IS(3)
      MV(4) = MV(1) + MV(2) + MV(3)
      MT(4) = MT(1) + MT(2) + MT(3)
      DO 180 KB=1,4
      IF(IS(KB).eq.0)  GO TO 180
      X(J,KB,4) = X(J,KB,4) + SCALEM(2)*MV(KB)*MT(KB)/IS(KB)
  180 CONTINUE
      DO 190 KB=1,4
  190 X(J,KB,6) = X(J,KB,2) - X(J,KB,4) - X(J,KB,5)
  200 CONTINUE
C**** I choose to define basins only where they can be integrated
C**** across from edge to edge. Gary has another convention.
C**** If that is required, uncomment these sections.
cC****
cC**** Calculate the transport entering each domain from an east-
cC**** west boundary that is north of the southern boundary
cC****
cC**** Bering Strait affects Atlantic and Pacific
c      SOIJL = 0.
c      DO L=1,LMO
c        SOIJL = SOIJL + OIJL(3,39,L,NOIJL)
c        IF (NOIJLGM.gt.0) SOIJL = SOIJL + OIJL(3,39,L,NOIJLGM)
c      END DO
c      DO 220 J=1,39-1
c      X(J,1,KQ) = X(J,1,KQ) + SOIJL*SCALEM(KQ)
c  220 X(J,2,KQ) = X(J,2,KQ) - SOIJL*SCALEM(KQ)
cC**** Lombock Passage affects Indian and Pacific
c      SOIJL = 0.
c      DO L=1,LMO
c      DO I=59,60
c        SOIJL = SOIJL + OIJL(I,21,L,NOIJL)
c        IF (NOIJLGM.gt.0) SOIJL = SOIJL + OIJL(I,21,L,NOIJLGM)
c      END DO
c      END Do
c      DO 240 J=1,21-1
c      X(J,2,KQ) = X(J,2,KQ) + SOIJL*SCALEM(KQ)
c  240 X(J,3,KQ) = X(J,3,KQ) - SOIJL*SCALEM(KQ)
cC**** Borneo Passage affects Indian and Pacific
c      SOIJL = 0.
c      DO L=1,LMO
c        SOIJL = SOIJL + OIJL(62,20,L,NOIJL)
c        IF (NOIJLGM.gt.0) SOIJL = SOIJL + OIJL(62,20,L,NOIJLGM)
c      END DO
c      DO 260 J=1,20-1
c      X(J,2,KQ) = X(J,2,KQ) + SOIJL*SCALEM(KQ)
c  260 X(J,3,KQ) = X(J,3,KQ) - SOIJL*SCALEM(KQ)
cC****
cC**** Calculate the transport entering each domain from a north-
cC**** south boundary that is north of the southern boundary
cC****
c      SELECT CASE (KQ)
c      CASE (1)
c        NOIJL=IJL_MFU
c        NOIJLGM=0
c      CASE (2)
c        NOIJL=IJL_GFLX
c        NOIJLGM=IJL_GGMFL
c      CASE (3)
c        NOIJL=IJL_SFLX
c        NOIJLGM=IJL_SGMFL
c      END SELECT
cC**** Drake Passage affects Pacific and Atlantic
c      SOIJL = 0.
c      DO 320 J=8,1,-1
c      DO L=1,LMO
c        SOIJL = SOIJL + OIJL(23,J+1,L,NOIJL)
c        IF (NOIJLGM.gt.0) SOIJL = SOIJL + OIJL(23,J+1,L,NOIJLGM)
c      END DO
c      X(J,1,KQ) = X(J,1,KQ) + SOIJL*SCALEM(KQ)
c  320 X(J,2,KQ) = X(J,2,KQ) - SOIJL*SCALEM(KQ)
cC**** Line south of Africa affects Atlantic and Indian
c      SOIJL = 0.
c      DO 340 J=13,1,-1
c      DO L=1,LMO
c        SOIJL = SOIJL + OIJL(40,J+1,L,NOIJL)
c        IF (NOIJLGM.gt.0) SOIJL = SOIJL + OIJL(40,J+1,L,NOIJLGM)
c      END DO
c      X(J,1,KQ) = X(J,1,KQ) - SOIJL*SCALEM(KQ)
c  340 X(J,3,KQ) = X(J,3,KQ) + SOIJL*SCALEM(KQ)
cC**** Line south of Austalia affects Indian and Pacific
c      SOIJL = 0.
c      DO 360 J=12,1,-1
c      DO L=1,LMO
c        SOIJL = SOIJL + OIJL(65,J+1,L,NOIJL)
c        IF (NOIJLGM.gt.0) SOIJL = SOIJL + OIJL(65,J+1,L,NOIJLGM)
c      END DO
c      X(J,2,KQ) = X(J,2,KQ) + SOIJL*SCALEM(KQ)
c  360 X(J,3,KQ) = X(J,3,KQ) - SOIJL*SCALEM(KQ)
C****
C**** Calculate transport through straits from latitude to another
C**** within the same basin
C****
      SELECT CASE (KQ)
      CASE (1)
        NOLNST=LN_MFLX
      CASE (2)
        NOLNST=LN_GFLX
      CASE (3)
        NOLNST=LN_SFLX
      END SELECT

      DO 400 NS=1,NMST
      SOLNST(NS) = 0.
      DO 400 L=1,LMO
  400 SOLNST(NS) = SOLNST(NS) + OLNST(L,NS,NOLNST)
C**** Fury & Hecla: (19,42) to (20,40)
      X(40,1,KQ) = X(40,1,KQ) - SOLNST(1)*SCALES(KQ)
      X(40,4,KQ) = X(40,4,KQ) - SOLNST(1)*SCALES(KQ)
      X(41,1,KQ) = X(41,1,KQ) - SOLNST(1)*SCALES(KQ)
      X(41,4,KQ) = X(41,4,KQ) - SOLNST(1)*SCALES(KQ)
C**** Nares: (22,43) to (24,44)
      X(43,1,KQ) = X(43,1,KQ) + SOLNST(2)*SCALES(KQ)
      X(43,4,KQ) = X(43,4,KQ) + SOLNST(2)*SCALES(KQ)
C**** Gibralter: (35,32) to (37,33)
      X(32,1,KQ) = X(32,1,KQ) + SOLNST(3)*SCALES(KQ)
      X(32,4,KQ) = X(32,4,KQ) + SOLNST(3)*SCALES(KQ)
C**** English: (36,36) to (37,37)
      X(36,1,KQ) = X(36,1,KQ) + SOLNST(4)*SCALES(KQ)
      X(36,4,KQ) = X(36,4,KQ) + SOLNST(4)*SCALES(KQ)
C**** Bosporous: (42,33) to (43,34)
      X(33,1,KQ) = X(33,1,KQ) + SOLNST(6)*SCALES(KQ)
      X(33,4,KQ) = X(33,4,KQ) + SOLNST(6)*SCALES(KQ)
C**** Red Sea: (44,29) to (45,28)
      X(28,3,KQ) = X(28,3,KQ) - SOLNST(7)*SCALES(KQ)
      X(28,4,KQ) = X(28,4,KQ) - SOLNST(7)*SCALES(KQ)
C**** Bab-al-Mandab: (45,28) to (46,27)
      X(27,3,KQ) = X(27,3,KQ) - SOLNST(8)*SCALES(KQ)
      X(27,4,KQ) = X(27,4,KQ) - SOLNST(8)*SCALES(KQ)
C**** Hormuz: (47,30) to (49,29)
      X(29,3,KQ) = X(29,3,KQ) - SOLNST(9)*SCALES(KQ)
      X(29,4,KQ) = X(29,4,KQ) - SOLNST(9)*SCALES(KQ)
C**** Korea: (62,32) to (63,33)
      X(32,2,KQ) = X(32,2,KQ) + SOLNST(11)*SCALES(KQ)
      X(32,4,KQ) = X(32,4,KQ) + SOLNST(11)*SCALES(KQ)
C**** Soya: (64,34) to (65,35)
      X(34,2,KQ) = X(34,2,KQ) + SOLNST(12)*SCALES(KQ)
      X(34,4,KQ) = X(34,4,KQ) + SOLNST(12)*SCALES(KQ)
C****
C**** Calculate transport through straits from one basin to another
C****
C**** Malacca: (56,25) to (58,24)
c      DO 510 J=1,23
c      X( J,2,KQ) = X( J,2,KQ) + SOLNST(10)*SCALES(KQ)
c  510 X( J,3,KQ) = X( J,3,KQ) - SOLNST(10)*SCALES(KQ)
      X(24,3,KQ) = X(24,3,KQ) - SOLNST(10)*SCALES(KQ)
      X(24,4,KQ) = X(24,4,KQ) - SOLNST(10)*SCALES(KQ)
C**** Fill in south polar value
      DO 610 KB=1,4
  610 X(0,KB,KQ)  = X(1,KB,KQ)
  620 CONTINUE
C**** Replace salt by salt - .035 times mass
      DO 630 KB=1,4
      DO 630 J=0,JM
  630 X(J,KB,3) = X(J,KB,3) - .035d3*X(J,KB,1)
C****
      IDLAT= NINT(180./(JM-1))
      FJEQ = JM/2.
      DO 660 J=1,JM-1
  660 VLAT( J) = IDLAT*(J-FJEQ)
      VLAT( 0) = -90.
      VLAT(JM) =  90.
C****
C**** Write titles and data to disk.
C****
      TITLE(51:80)=XLB
      DO KQ=1,3
        TITLE(1:50)=NAME(KQ)
        WRITE (6,907) TITLE(1:72)
C**** print out truncated series to PRT file
        WRITE (6,903) (NINT(VLAT(J)),J=JM,0,-INC)
        DO KB=1,4
          WRITE (6,906) BASIN(KB),(X(J,KB,KQ),J=JM,0,-INC)
        END DO
        WRITE(6,904)
        IF (QDIAG) THEN
          WRITE (iu_otj,*) TITLE
          WRITE (iu_otj,*) 'Latitude'
          WRITE (iu_otj,*) YAXIS(KQ)
          WRITE (iu_otj,*)
     *         ' Lat   Atlantic   Pacific    Indian     Global'
          WRITE (iu_otj,971) (VLAT(J),(X(J,KB,KQ),KB=1,4),J=0,JM)
          WRITE (iu_otj,*)
        END IF
      END DO
C****
C**** Write titles and data to disk for northward transport of heat
C**** by components
C****
      WRITE(6,908)
      DO KB=1,4
      TITLE(1:50)=NAME(3+KB)
      WRITE (6,907) TITLE(1:72)
C**** print out truncated series to PRT file
      WRITE(6,903) (NINT(VLAT(J)),J=JM,0,-INC)
      WRITE(6,906)BASIN(KB)(1:3)//" Overtrn",(X(J,KB,4),J=JM,0,-INC)
      WRITE(6,906)BASIN(KB)(1:3)//" GM flx ",(X(J,KB,5),J=JM,0,-INC)
      WRITE(6,906)BASIN(KB)(1:3)//" Hor gyr",(X(J,KB,6),J=JM,0,-INC)
      WRITE(6,906)BASIN(KB)(1:3)//" Total  ",(X(J,KB,2),J=JM,0,-INC)
      WRITE(6,904)
      IF (QDIAG) THEN
        WRITE (iu_otj,*) TITLE
        WRITE (iu_otj,*) 'Latitude'
        WRITE (iu_otj,*) YAXIS(2)
        WRITE (iu_otj,*)
     *       ' Lat      Total    Overturn    GM_flux    Hor_Gyre'
        WRITE (iu_otj,976)(VLAT(J),X(J,KB,2),X(J,KB,4),X(J,KB,5),X(J,KB
     *       ,6),J=0,JM)
        WRITE (iu_otj,*)
      END IF
      END DO
      RETURN
C****
  970 FORMAT (A50,A6,2X,2A4)
  971 FORMAT (F6.0,4F10.3)
  972 FORMAT (A,1P,3E15.4)
  973 FORMAT (A,2I6)
  976 FORMAT (F6.0,4F10.3)
  903 FORMAT (' ',131('-')/,' Latitude  ',24I5)
  904 FORMAT (' ',131('-'))
  906 FORMAT (' ',A11,24F5.1)
  907 FORMAT ('0',A)
  908 FORMAT ('1')

      END SUBROUTINE OTJOUT
