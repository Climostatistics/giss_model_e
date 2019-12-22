#include "rundeck_opts.h"

!@sum TRACER_PRT Tracer diagnostics.  These routines are generic for
!@+    all tracers
!@+   TRACEA and DIAGTCA are called throughout the day
!@+   The other routines are for printing
!@auth Jean Lerner (with ideas stolen from G.Schmidt, R. Ruedy, etc.)

#ifdef TRACERS_ON
      SUBROUTINE TRACEA
!@sum TRACEA accumulates tracer concentration diagnostics (IJL, JL)
!@auth J.Lerner
!@ver  1.0
      USE MODEL_COM, only: im,jm,lm,itime
      USE GEOM, only: imaxj,bydxyp
      USE SOMTQ_COM, only: mz  
      USE TRACER_COM 
      USE TRACER_DIAG_COM 
      USE DYNAMICS, only: am,byam 
      implicit none

      integer i,j,l,k,n
      real*8 tsum,asum

C****
C**** Accumulate concentration for all tracers
C****
      do 600 n=1,ntm
      IF (itime.lt.itime_tr0(n)) cycle
C**** Latitude-longitude by layer concentration
!$OMP PARALLEL DO PRIVATE (L)
      do l=1,lm
        taijln(:,:,l,n) = taijln(:,:,l,n) + trm(:,:,l,n)*byam(l,:,:)
      end do
!$OMP END PARALLEL DO
C**** Average concentration; surface concentration; total mass
!$OMP PARALLEL DO PRIVATE (J,I,TSUM,ASUM)
      do j=1,jm
      do i=1,im
        tsum = sum(trm(i,j,:,n))*bydxyp(j)  !sum over l
        asum = sum(am(:,i,j))     !sum over l
        taijn(i,j,tij_mass,n) = taijn(i,j,tij_mass,n)+tsum  !MASS 
        taijn(i,j,tij_conc,n) = taijn(i,j,tij_conc,n)+tsum/asum
      enddo; enddo
!$OMP END PARALLEL DO
C**** Zonal mean concentration
!$OMP PARALLEL DO PRIVATE (L,J,TSUM,ASUM)
      do l=1,lm
      do j=1,jm
        tsum = sum(trm(1:imaxj(j),j,l,n)) !sum over i
        asum = sum(am(l,1:imaxj(j),j))    !sum over i
        tajln(j,l,jlnt_conc,n) = tajln(j,l,jlnt_conc,n)+tsum/asum
      enddo; enddo
!$OMP END PARALLEL DO
  600 continue
      return
      end SUBROUTINE TRACEA


      SUBROUTINE DIAGTCA (M,NT)
!@sum  DIAGTCA Keeps track of the conservation properties of tracers
!@auth Gary Russell/Gavin Schmidt/Jean Lerner
!@ver  1.0
      USE MODEL_COM, only: im,jm,lm,fim
      USE GEOM, only: imaxj
      USE TRACER_COM
      USE TRACER_DIAG_COM, only: tconsrv,nofmt,title_tcon
      IMPLICIT NONE

!@var M index denoting which process changed the tracer
      INTEGER, INTENT(IN) :: m
!@var NT index denoting tracer number
      INTEGER, INTENT(IN) :: nt
!@var TOTAL amount of conserved quantity at this time
      REAL*8, DIMENSION(JM) :: total
      INTEGER :: i,j,l,nm,ni
      REAL*8 sstm,stm
C****
C**** THE PARAMETER M INDICATES WHEN DIAGCA IS BEING CALLED
C**** M=1,2...12:  See DIAGCA in DIAG.f
C****   13+ AFTER Sources and Sinks
C****
C**** NOFMT contains the indexes of the TCONSRV array where each
C**** change is to be stored for each quantity. If NOFMT(M,NT)=0,
C**** no calculation is done.
C**** NOFMT(1,NT) is the index for the instantaneous value.

      if (nofmt(m,nt).gt.0) then
C**** Calculate current value TOTAL
!$OMP PARALLEL DO PRIVATE (J,L,I,SSTM,STM)
        do j=1,jm
          sstm = 0.
          do l=1,lm
            stm = 0.
            do i=1,imaxj(j)
              stm = stm+trm(i,j,l,nt)
#ifdef TRACERS_WATER
     *             +trwm(i,j,l,nt)
#endif
            end do
            sstm = sstm+stm
          end do
          total(j) = sstm
        end do
!$OMP END PARALLEL DO
        total(1) = fim*total(1)
        total(jm)= fim*total(jm)
        nm=nofmt(m,nt)
        ni=nofmt(1,nt)
c**** Accumulate difference from last time in TCONSRV(NM)
        if (m.gt.1) then
          tconsrv(:,nm,nt) = 
     &    tconsrv(:,nm,nt)+(total(:)-tconsrv(:,ni,nt))  !do 1,jm
        end if
C**** Save current value in TCONSRV(NI)
        tconsrv(:,ni,nt)=total(:)   !do 1,jm
      end if
      return
      end subroutine diagtca

      SUBROUTINE DIAGTCB (DTRACER,M,NT)
!@sum  DIAGTCB Keeps track of the conservation properties of tracers
!@+    This routine takes an already calculated difference
!@auth Gary Russell/Gavin Schmidt/Jean Lerner
!@ver  1.0
      USE MODEL_COM, only: jm,fim
      USE TRACER_DIAG_COM, only: tconsrv,nofmt,title_tcon
      IMPLICIT NONE

!@var M index denoting which process changed the tracer
      INTEGER, INTENT(IN) :: m
!@var NT index denoting tracer number
      INTEGER, INTENT(IN) :: nt
!@var DTRACER change of conserved quantity at this time
      REAL*8, DIMENSION(JM) :: DTRACER
      INTEGER :: j,nm
      REAL*8 stm
C****
C**** THE PARAMETER M INDICATES WHEN DIAGCA IS BEING CALLED
C**** M=1,2...12:  See DIAGCA in DIAG.f
C****   13+ AFTER Sources and Sinks
C****
C**** NOFMT contains the indexes of the TCONSRV array where each
C**** change is to be stored for each quantity. If NOFMT(M,NT)=0,
C**** no calculation is done.

      if (nofmt(m,nt).gt.0) then
C**** Calculate latitudinal mean of chnage DTRACER
        dtracer(1) = fim*dtracer(1)
        dtracer(jm)= fim*dtracer(jm)
        nm=nofmt(m,nt)
c**** Accumulate difference in TCONSRV(NM)
        if (m.gt.1) then
          tconsrv(:,nm,nt) = tconsrv(:,nm,nt) + dtracer(:)
        end if
C**** No need to save current value
      end if
      return
      end subroutine diagtcb

      SUBROUTINE DIAGTCP
!@sum  DIAGCP produces tables of the conservation diagnostics
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only: teeny, twopi
      USE MODEL_COM, only:
     &     jm,fim,idacc,jhour,jhour0,jdate,jdate0,amon,amon0,
     &     jyear,jyear0,nday,jeq,itime,itime0,xlabel
      USE GEOM, only:
     &     areag,dlon,dxyp,lat_dg
      USE TRACER_COM, only: ntm ,itime_tr0
      USE TRACER_DIAG_COM, only:
     &     TCONSRV,ktcon,scale_tcon,title_tcon,nsum_tcon,ia_tcon,nofmt
      USE DAGCOM, only: inc=>incj,xwon,kdiag
      IMPLICIT NONE

      INTEGER, DIMENSION(JM) :: MAREA
      REAL*8, DIMENSION(KTCON) :: FGLOB
      REAL*8, DIMENSION(2,KTCON) :: FHEM
      REAL*8, DIMENSION(JM+3,KTCON,NTM) :: CNSLAT
      CHARACTER*4, PARAMETER :: HEMIS(2) = (/' SH ',' NH '/),
     *     DASH = ('----')
      INTEGER :: j,jhemi,jnh,jp1,jpm,jsh,jx,n,k,KTCON_max
      REAL*8 :: aglob,ahem,feq,fnh,fsh,days

      if (kdiag(8).ge.2) return
C**** CALCULATE SCALING FACTORS
      IF (IDACC(12).LT.1) IDACC(12)=1
C**** Calculate areas
      DO J=1,JM
        MAREA(J)=1.D-10*XWON*FIM*DXYP(J)+.5
      END DO
      AGLOB=1.D-10*AREAG*XWON
      AHEM=1.D-10*(.5*AREAG)*XWON
C****
C**** Outer loop over tracers
C****
      DO 900 N=1,NTM
      IF (itime.LT.itime_tr0(N)) cycle
C**** CALCULATE SUM OF CHANGES
C**** LOOP BACKWARDS SO THAT INITIALIZATION IS DONE BEFORE SUMMATION!
      DO J=1,JM
        DO K=KTCON,1,-1
          IF (NSUM_TCON(K,N).eq.0) THEN
            TCONSRV(J,K,N)=0.
          ELSEIF (NSUM_TCON(K,N).gt.0) THEN
            TCONSRV(J,NSUM_TCON(K,N),N)=
     *      TCONSRV(J,NSUM_TCON(K,N),N)+TCONSRV(J,K,N)
     *           *SCALE_TCON(K,N)*IDACC(12)/(IDACC(IA_TCON(K,N))+teeny)
          KTCON_max = NSUM_TCON(K,N)
          END IF
        END DO
      END DO
C**** CALCULATE ALL OTHER CONSERVED QUANTITIES ON TRACER GRID
      DO K=1,KTCON_max
        FGLOB(K)=0.
        FHEM(1,K)=0.
        FHEM(2,K)=0.
        DO JSH=1,JEQ-1
          JNH=1+JM-JSH
          FSH=TCONSRV(JSH,K,N)*SCALE_TCON(K,N)/(IDACC(IA_TCON(K,N))
     *         +teeny)
          FNH=TCONSRV(JNH,K,N)*SCALE_TCON(K,N)/(IDACC(IA_TCON(K,N))
     *         +teeny)
          FGLOB (K)=FGLOB (K)+(FSH+FNH)
          FHEM(1,K)=FHEM(1,K)+FSH
          FHEM(2,K)=FHEM(2,K)+FNH
          CNSLAT(JSH,K,N)=FSH/(FIM*DXYP(JSH))
          CNSLAT(JNH,K,N)=FNH/(FIM*DXYP(JNH))
        END DO
        FGLOB (K)=FGLOB (K)/AREAG
        FHEM(1,K)=FHEM(1,K)/(.5*AREAG)
        FHEM(2,K)=FHEM(2,K)/(.5*AREAG)
        CNSLAT(JM+1,K,N)=FHEM(1,K)
        CNSLAT(JM+2,K,N)=FHEM(2,K)
        CNSLAT(JM+3,K,N)=FGLOB (K)
      END DO
C**** LOOP OVER HEMISPHERES
      DAYS=(Itime-Itime0)/DFLOAT(nday)
      WRITE (6,'(''1'')')
      DO JHEMI=2,1,-1
        WRITE (6,'(''0'',A)') XLABEL
        WRITE (6,902) JYEAR0,AMON0,JDATE0,JHOUR0,
     *       JYEAR,AMON,JDATE,JHOUR,ITIME,DAYS
        JP1=1+(JHEMI-1)*(JEQ-1)
        JPM=JHEMI*(JEQ-1)
C**** PRODUCE TABLES 
        WRITE (6,903) (DASH,J=JP1,JPM,INC)
        WRITE (6,904) HEMIS(JHEMI),(NINT(LAT_DG(JX,1)),JX=JPM,JP1,-INC)
        WRITE (6,903) (DASH,J=JP1,JPM,INC)
        DO K=1,KTCON_max
        WRITE (6,905) TITLE_TCON(K,N),FGLOB(K),FHEM(JHEMI,K),
     *         (NINT(CNSLAT(JX,K,N)),JX=JPM,JP1,-INC)
        END DO
        WRITE (6,906) AGLOB,AHEM,(MAREA(JX),JX=JPM,JP1,-INC)
      END DO
  900 CONTINUE
      RETURN
C****
  902 FORMAT ('0Conservation Quantities       From:',
     *  I6,A6,I2,',  Hr',I3,  6X,  'To:',I6,A6,I2,', Hr',I3,
     *  '  Model-Time:',I9,5X,'Dif:',F7.2,' Days')
  903 FORMAT (1X,28('--'),13(A4,'--'))
  904 FORMAT (41X,'GLOBAL',A8,2X,13I6)
  905 FORMAT (A38,2F9.2,1X,13I6)
  906 FORMAT ('0AREA (10^10 M^2)',F30.1,F9.1,1X,13I6)
      END SUBROUTINE DIAGTCP


      MODULE BDjlt
!@sum  stores info for outputting lat-sig/pressure diags for tracers
!@auth J. Lerner
      IMPLICIT NONE
      SAVE
!@param names of derived JLt/JLs tracer output fields
      INTEGER jlnt_nt_eddy,jlnt_vt_eddy  
      END MODULE BDjlt


      SUBROUTINE JLt_TITLEX
!@sum JLt_TITLEX sets up titles, etc. for composit JL output (tracers)
!@auth J. Lerner
!@ver  1.0
      USE TRACER_COM
      USE TRACER_DIAG_COM
      USE BDjlt
      IMPLICIT NONE
      character*50 :: unit_string
      integer k,n,kpmax

!@var jlnt_xx Names for TAJL diagnostics
      do n=1,ntm
      k = ktajl +1
        jlnt_nt_eddy = k
        sname_jln(k,n) = 'tr_nt_eddy_'//trname(n)
        lname_jln(k,n) = 'NORTHWARD TRANS. OF '//
     &     trim(trname(n))//' MASS BY EDDIES'
        jlq_power(k) = 10.
        jgrid_jlq(k) = 2
      k = k+1
        jlnt_vt_eddy = k
        sname_jln(k,n) = 'tr_vt_eddy_'//trname(n)
        lname_jln(k,n) = 'VERTICAL TRANS. OF '//
     &     trim(trname(n))//' MASS BY EDDIES'
        jlq_power(k) = 10.
      end do

      kpmax = k

      if (k .gt. ktajlx) then
        write (6,*) 'JLt_TITLEX: Increase ktajlx=',ktajlx,
     &         ' to at least ',k
        call stop_model('ktajlx too small',255)
      end if
C**** Construct UNITS string for output
      do n=1,ntm
      do k=ktajl+1,kpmax
      units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),' kg/s')
      end do; end do

      END SUBROUTINE JLt_TITLEX


      SUBROUTINE DIAGJLT
!@sum DIAGJLT controls calls to the pressure/lat print routine for 
!@+      tracers
!@auth J. Lerner
!@calls open_jl, JLt_TITLEX, JLMAP_t
C****
C**** THIS ROUTINE PRODUCES LATITUDE BY LAYER TABLES OF TRACERS
C****
      USE CONSTANT, only : undef
      USE MODEL_COM, only: jm,lm,fim,itime,idacc,xlabel,lrunid,psfmpt
     &   ,sige,ptop
      USE GEOM, only: bydxyp,dxyp,lat_dg
      USE TRACER_COM
      USE DAGCOM, only: linect,plm,acc_period,qdiag,lm_req
      USE TRACER_DIAG_COM
      USE BDJLT
      IMPLICIT NONE

      REAL*8, DIMENSION(JM) :: ONESPO
      REAL*8, DIMENSION(JM+LM) :: ONES
      REAL*8, DIMENSION(JM,LM) :: A
      REAL*8, DIMENSION(LM) :: PM
      REAL*8 :: scalet
      INTEGER :: J,L,N,K,jtpow
#ifdef TRACERS_SPECIAL_O18
      INTEGER :: n1,n2
      CHARACTER :: lname*80,sname*30,units*50
      REAL*8 :: dD, d18O
#endif

C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
      IF(QDIAG) call open_jl(trim(acc_period)//'.jlt'//XLABEL(1:LRUNID)
     *  ,jm,lm,0,lat_dg)
!?   *  ,jm,lm,lm_req,lat_dg)
C****
C**** INITIALIZE CERTAIN QUANTITIES
C****
      call JLt_TITLEX

      do l=1,lm
        pm(l)=psfmpt*sige(l+1)+ptop
      end do
      onespo(:)  = 1.d0
      onespo(1)  = fim
      onespo(jm) = fim
      ones(:) = 1.d0
      linect = 65
C****
C**** LOOP OVER TRACERS
C****
      DO 400 N=1,NTM
      IF (itime.LT.itime_tr0(N)) cycle
C****
C**** TRACER CONCENTRATION
C****
      k = jlnt_conc

#ifdef TRACERS_WATER
      if (to_per_mil(n).gt.0) then
C**** Note permil concentrations REQUIRE trw0 and n_water to be defined!
      scalet = 1.
      jtpow = 0.
      do l=1,lm
      do j=1,jm
        if (tajln(j,l,k,n_water).gt.0) then
          a(j,l)=1d3*(tajln(j,l,k,n)/(trw0(n)*tajln(j,l,k,n_water))-1.)
        else
          a(j,l)=undef 
        end if
      end do
      end do
      CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *     plm,a,scalet,ones,ones,lm,2,jgrid_jlq(k))
      else
#endif
      scalet = scale_jln(n)*scale_jlq(k)/idacc(ia_jlq(k))
      jtpow = ntm_power(n)+jlq_power(k)
      scalet = scalet*10.**(-jtpow)
      CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *     plm,tajln(1,1,k,n),scalet,bydxyp,ones,lm,2,jgrid_jlq(k))
#ifdef TRACERS_WATER
      end if
#endif
C****
C**** NORTHWARD TRANSPORTS: Total and eddies
C****
      a(:,:) = 0.
      k = jlnt_nt_tot
      do 205 l=1,lm
      do 205 j=1,jm-1
  205 a(j+1,l) = tajln(j,l,k,n)
      scalet = scale_jlq(k)/idacc(ia_jlq(k))
      jtpow = ntm_power(n)+jlq_power(k)
      scalet = scalet*10.**(-jtpow)
      CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *  plm,a,scalet,ones,ones,lm,1,jgrid_jlq(k))
      do 210 l=1,lm
      do 210 j=1,jm-1
  210 a(j+1,l) = tajln(j,l,k,n) -tajln(j,l,jlnt_nt_mm,n)
      k = jlnt_nt_eddy
      scalet = scale_jlq(k)/idacc(ia_jlq(k))
      jtpow = ntm_power(n)+jlq_power(k)
      scalet = scalet*10.**(-jtpow)
      CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *  plm,a,scalet,ones,ones,lm,1,jgrid_jlq(k))
C****
C**** VERTICAL TRANSPORTS: Total and eddies
C****
      k = jlnt_vt_tot
      scalet = scale_jlq(k)/idacc(ia_jlq(k))
      jtpow = ntm_power(n)+jlq_power(k)
      scalet = scalet*10.**(-jtpow)
      CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *  pm,tajln(1,1,k,n),scalet,ones,ones,lm-1,1,jgrid_jlq(k))
      a(:,:) = 0.
      do 260 l=1,lm-1
      do 260 j=1,jm
  260 a(j,l) = tajln(j,l,k,n)-tajln(j,l,jlnt_vt_mm,n)
      k = jlnt_vt_eddy
      scalet = scale_jlq(k)/idacc(ia_jlq(k))
      jtpow = ntm_power(n)+jlq_power(k)
      scalet = scalet*10.**(-jtpow)
      CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *  pm,a,scalet,ones,ones,lm-1,1,jgrid_jlq(k))
C****
C**** Convective Processes
C****
C**** MOIST CONVECTION 
      k = jlnt_mc
      scalet = scale_jlq(k)/idacc(ia_jlq(k))
      jtpow = ntm_power(n)+jlq_power(k)
      scalet = scalet*10.**(-jtpow)
      CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *  plm,tajln(1,1,k,n),scalet,ones,ones,lm,1,Jgrid_jlq(k))
C**** TURBULENCE (or Dry Convection)
      k = jlnt_turb
      scalet = scale_jlq(k)/idacc(ia_jlq(k))
      jtpow = ntm_power(n)+jlq_power(k)
      scalet = scalet*10.**(-jtpow)
      CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *  plm,tajln(1,1,k,n),scalet,ones,ones,lm,1,jgrid_jlq(k))
C**** LARGE-SCALE CONDENSATION
      k = jlnt_lscond
      scalet = scale_jlq(k)/idacc(ia_jlq(k))
      jtpow = ntm_power(n)+jlq_power(k)
      scalet = scalet*10.**(-jtpow)
      CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *  plm,tajln(1,1,k,n),scalet,ones,ones,lm,1,jgrid_jlq(k))
  400 CONTINUE
C****
C**** Sources and sinks
C****
      do k=1,ktajls
        n = jls_index(k)
        if ( n.eq.0 ) cycle
        if ( itime.LT.itime_tr0(n) )  cycle
        scalet = scale_jls(k)*10.**(-jls_power(k))/idacc(ia_jls(k))
        CALL JLMAP_t (lname_jls(k),sname_jls(k),units_jls(k),
     *    plm,tajls(1,1,k),scalet,ones,ones,jls_ltop(k),1,jgrid_jls(k))
      end do

#ifdef TRACERS_SPECIAL_O18
C****
C**** Calculations of deuterium excess (d=dD-8*d18O)
C**** Note: some of these definitions probably belong in JLt_TITLEX
C****
      if (n_H2O18.gt.0 .and. n_HDO.gt.0) then
        n1=n_H2O18
        n2=n_HDO
        k=jlnt_conc
        scalet = 1.
        jtpow = 0.
        
        do l=1,lm
          do j=1,jm
            if (tajln(j,l,k,n_water).gt.0) then
              d18O=1d3*(tajln(j,l,k,n1)/(trw0(n1)*tajln(j,l,k,n_water))
     *             -1.)
              dD=1d3*(tajln(j,l,k,n2)/(trw0(n2)*tajln(j,l,k,n_water))
     *             -1.)
              a(j,l)=dD-8.*d18O
            else
              a(j,l)=undef 
            end if
          end do
        end do
        lname="Deuterium excess"
        sname="dexcess"
        units="per mil"
        CALL JLMAP_t (lname,sname,units,plm,a,scalet,ones,ones,lm,2
     *       ,jgrid_jlq(k)) 
      end if
#endif

      if (qdiag) call close_jl
      RETURN
      END SUBROUTINE DIAGJLT
#endif

      SUBROUTINE JLMAP_t (LNAME,SNAME,UNITS,
     &     PL,AX,SCALET,SCALEJ,SCALEL,LMAX,JWT,JG)
C****
C**** THIS ROUTINE PRODUCES LAYER BY LATITUDE TABLES ON THE LINE
C**** PRINTER.  THE INTERIOR NUMBERS OF THE TABLE ARE CALCULATED AS
C****               AX * SCALET * SCALEJ * SCALEL.
C**** WHEN JWT=1, THE INSIDE NUMBERS ARE NOT AREA WEIGHTED AND THE
C****    HEMISPHERIC AND GLOBAL NUMBERS ARE SUMMATIONS.
C**** WHEN JWT=2, ALL NUMBERS ARE PER UNIT AREA.
C**** JG INDICATES PRIMARY OR SECONDARY GRID.
C**** THE BOTTOM LINE IS CALCULATED AS THE SUMMATION OF DSIG TIMES THE
C**** NUMBERS ABOVE 
C****
      USE MODEL_COM, only: jm,lm,jdate,jdate0,amon,amon0,jyear,jyear0
     *     ,xlabel,dsig,sige 
      USE GEOM, only: wtj,jrange_hemi,lat_dg
      USE DAGCOM, only: qdiag,acc_period,inc=>incj,linect,jmby2,lm_req
      IMPLICIT NONE

!@var units string containing output field units
      CHARACTER(LEN=50) :: UNITS
!@var lname string describing output field
      CHARACTER(LEN=50) :: LNAME
!@var sname string referencing output field
      CHARACTER(LEN=30) :: SNAME
!@var title string, formed as concatentation of lname//units
      CHARACTER(LEN=64) :: TITLE

      INTEGER :: JG,JWT,LMAX
      REAL*8 :: SCALET
      REAL*8, DIMENSION(JM,LM) :: AX
      REAL*8, DIMENSION(JM) :: SCALEJ
      REAL*8, DIMENSION(LM) :: SCALEL,PL

      CHARACTER*4 DASH,WORD(4)
      DATA DASH/'----'/,WORD/'SUM','MEAN',' ','.1*'/

      INTEGER, DIMENSION(JM) :: MLAT
      REAL*8, DIMENSION(JM) :: ASUM
      REAL*8, DIMENSION(2) :: FHEM,HSUM
      INTEGER :: IWORD,J,JHEMI,K,L
      REAL*8 :: FGLOB,FLATJ,GSUM,SDSIG

!?    REAL*8, DIMENSION(JM+3,LM+LM_REQ+1) :: XJL ! for binary output
      REAL*8, DIMENSION(JM+3,LM+1) :: XJL ! for binary output
      CHARACTER XLB*16,CLAT*16,CPRES*16,CBLANK*16,TITLEO*80
      DATA CLAT/'LATITUDE'/,CPRES/'PRESSURE (MB)'/,CBLANK/' '/

C form title string
      title = trim(lname)//' ('//trim(units)//')'

C****
C**** WRITE XLABEL ON THE TOP OF EACH OUTPUT PAGE
C****
      LINECT = LINECT+LMAX+7
      IF(LINECT.LE.60) GO TO 20
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
      LINECT = LMAX+8
   20 CONTINUE
      WRITE (6,901) TITLE,(DASH,J=JG,JM,INC)
      WRITE (6,904) WORD(JWT),(NINT(LAT_DG(J,JG)),J=JM,JG,-INC)
      WRITE (6,905) (DASH,J=JG,JM,INC)
C****
C**** CALCULATE TABLE NUMBERS AND WRITE THEM TO THE LINE PRINTER
C****
         XJL(:,:) = -1.E30
      SDSIG = 1.-SIGE(LMAX+1)
      ASUM(:) = 0.
      HSUM(:) = 0.
      GSUM = 0.
      DO 240 L=LMAX,1,-1
      FGLOB = 0.
      DO 230 JHEMI=1,2
      FHEM(JHEMI) = 0.
      DO 220 J=JRANGE_HEMI(1,JHEMI,JG),JRANGE_HEMI(2,JHEMI,JG)
      FLATJ = AX(J,L)*SCALET*SCALEJ(J)*SCALEL(L)
         XJL(J,L) = FLATJ
      MLAT(J) = NINT(FLATJ)
      IF (JWT.EQ.1) THEN     
        ASUM(J) = ASUM(J)+FLATJ  !!!!most
      ELSE
        ASUM(J) = ASUM(J)+FLATJ*DSIG(L)/SDSIG  !!!! concentration
      ENDIF
      FHEM(JHEMI) = FHEM(JHEMI)+FLATJ*WTJ(J,JWT,JG)
  220 CONTINUE
  230 FGLOB = FGLOB+FHEM(JHEMI)/JWT
         XJL(JM+3,L)=FHEM(1)   ! SOUTHERN HEM
         XJL(JM+2,L)=FHEM(2)   ! NORTHERN HEM
         XJL(JM+1,L)=FGLOB     ! GLOBAL
      WRITE (6,902) PL(L),FGLOB,FHEM(2),FHEM(1),(MLAT(J),J=JM,JG,-INC)
      IF (JWT.EQ.1) THEN
        HSUM(1) = HSUM(1)+FHEM(1)
        HSUM(2) = HSUM(2)+FHEM(2)
        GSUM    = GSUM   +FGLOB
      ELSE
        HSUM(1) = HSUM(1)+FHEM(1)*DSIG(L)/SDSIG
        HSUM(2) = HSUM(2)+FHEM(2)*DSIG(L)/SDSIG
        GSUM    = GSUM   +FGLOB  *DSIG(L)/SDSIG
      ENDIF
  240 CONTINUE
      WRITE (6,905) (DASH,J=JG,JM,INC)
      ASUM(jmby2+1) = ASUM(jmby2+1)/JG
         DO 180 J=JG,JM
! 180    XJL(J   ,LM+LM_REQ+1)=ASUM(J)
!        XJL(JM+3,LM+LM_REQ+1)=HSUM(1)   ! SOUTHERN HEM
!        XJL(JM+2,LM+LM_REQ+1)=HSUM(2)   ! NORTHERN HEM
!        XJL(JM+1,LM+LM_REQ+1)=GSUM      ! GLOBAL
  180    XJL(J   ,LM+1)=ASUM(J)
         XJL(JM+3,LM+1)=HSUM(1)   ! SOUTHERN HEM
         XJL(JM+2,LM+1)=HSUM(2)   ! NORTHERN HEM
         XJL(JM+1,LM+1)=GSUM      ! GLOBAL
         XLB=' '//acc_period(1:3)//' '//acc_period(4:12)//'  '
         TITLEO=TITLE//XLB
         IF(QDIAG) CALL POUT_JL(TITLEO,LNAME,SNAME,UNITS,
     *        JG,LMAX,XJL,PL,CLAT,CPRES)
      IF (LMAX.EQ.1) THEN
         LINECT = LINECT-1
         RETURN
      ENDIF
      IF (LMAX.EQ.1) RETURN
      WRITE (6,903) WORD(JWT),GSUM,HSUM(2),HSUM(1),
     *    (NINT(ASUM(J)),J=JM,JG,-INC)
      RETURN
C****
  901 FORMAT ('0',30X,A64/2X,32('-'),24A4)
  902 FORMAT (1X,F8.3,3F8.1,1X,24I4)
  903 FORMAT (1X,A6,2X,3F8.1,1X,24I4)
  904 FORMAT ('  P(MB)    ',A4,' G     NH      SH  ',24I4)
  905 FORMAT (2X,32('-'),24A4)
  907 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)
      END

#ifdef TRACERS_ON
      SUBROUTINE DIAGIJt
!@sum  DIAGIJt produces lat-lon fields as maplets (6/page) or full-page
!@+    digital maps, and binary (netcdf etc) files (if qdiag=true)
!@auth Jean Lerner (adapted from work of G. Russell,M. Kelley,R. Ruedy)
!@ver   1.0
      USE MODEL_COM, only: im,jm,lm,jhour,jhour0,jdate,jdate0,amon,amon0
     *     ,jyear,jyear0,nday,itime,itime0,xlabel,lrunid
      USE TRACER_COM
      USE DAGCOM
      USE TRACER_DIAG_COM
      IMPLICIT NONE

      integer, parameter :: ktmax = (lm+ktaij)*ntm+ktaijs
#ifdef TRACERS_SPECIAL_O18
     *     + 2+lm      ! include dexcess diags
#endif
!@var Qk: if Qk(k)=.true. field k still has to be processed
      logical, dimension (ktmax) :: Qk
!@var Iord: Index array, fields are processed in order Iord(k), k=1,2,..
!@+     only important for fields 1->nmaplets which appear in printout
!@+     Iord(k)=0 indicates that a blank space replaces a maplet
      INTEGER nmaplets
      INTEGER, DIMENSION(ktmax) :: nt,ijtype,Iord,irange,iacc
      CHARACTER, DIMENSION(ktmax) :: lname*80,name*30,units*30
      REAL*8, DIMENSION(ktmax) :: scale
      REAL*8, DIMENSION(IM,JM,ktmax) :: aij1,aij2
      REAL*8, DIMENSION(IM,JM) :: SMAP
      REAL*8, DIMENSION(JM) :: SMAPJ
      CHARACTER xlb*32,title*48
!@var LINE virtual half page (with room for overstrikes)
      CHARACTER*133 LINE(53)
      INTEGER ::  I,J,K,kx,L,M,N,kcolmn,nlines,jgrid,n1,n2
      REAL*8 :: DAYS,gm

      if (kdiag(8).ge.1) return
C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
      IF(QDIAG)call open_ij(trim(acc_period)//'.ijt'//XLABEL(1:LRUNID) 
     *     ,im,jm)

c**** always skip unused fields
      Qk = .true.

C**** Fill in the undefined pole box duplicates
      do i=2,im
        taijln(i,1,:,:) = taijln(1,1,:,:)
        taijln(i,jm,:,:) = taijln(1,jm,:,:)
        taijn(i,1,:,:) = taijn(1,1,:,:)
        taijn(i,jm,:,:) = taijn(1,jm,:,:)
        taijs (i,1,:) = taijs(1,1,:)
        taijs (i,jm,:) = taijs(1,jm,:)
      end do

C**** Fill in maplet indices for tracer concentrations
      k = 0
      do n=1,ntm
      if (itime.lt.itime_tr0(n)) cycle
      do l=1,lm
        k = k+1
        iord(k) = l
        nt(k) = n
        ijtype(k) = 1
        name(k) = sname_ijt(l,n) 
        lname(k) = lname_ijt(l,n) 
        units(k) = units_ijt(l,n)
        irange(k) = ir_ijt(n)
        iacc(k) = ia_ijt
        aij1(:,:,k) = taijln(:,:,l,n)
        aij2(:,:,k) = 1.
        scale(k) = scale_ijt(l,n)
#ifdef TRACERS_WATER
        if (to_per_mil(n).gt.0) then
        aij1(:,:,k)=1d3*(taijln(:,:,l,n)-taijln(:,:,l,n_water)*trw0(n))
        aij2(:,:,k)=taijln(:,:,l,n_water)*trw0(n)
        ijtype(k) = 3
        end if
#endif
        if (index(lname_ijt(l,n),'unused').gt.0) Qk(k) = .false.
      end do
C**** Fill in maplet indices for tracer sums/means and ground conc
      do l=1,ktaij
        if (index(lname_tij(l,n),'unused').gt.0) cycle
        k = k+1
        iord(k) = l
        nt(k) = n
        name(k) = sname_tij(l,n)
        lname(k) = lname_tij(l,n)
        units(k) = units_tij(l,n)
        irange(k) = ir_ijt(n)
        scale(k) = scale_tij(l,n)
        iacc(k) = ia_ijt
        ijtype(k) = 2
        aij1(:,:,k) = taijn(:,:,l,n)
        aij2(:,:,k) = 1.
#ifdef TRACERS_WATER
        if (to_per_mil(n).gt.0 .and. l.ne.tij_mass) then
        aij1(:,:,k)=1d3*(taijn(:,:,l,n)-taijn(:,:,l,n_water)*trw0(n))
        aij2(:,:,k)=taijn(:,:,l,n_water)*trw0(n)
        ijtype(k) = 3
        end if
#endif
      end do
      end do

C**** Fill in maplet indices for sources and sinks
      do kx=1,ktaijs
        if (index(lname_ijts(kx),'unused').gt.0) cycle
        n = ijts_index(kx)
        if (itime.lt.itime_tr0(n)) cycle
        k = k+1
        iord(k) = kx
        nt(k) = kx
        ijtype(k) = 1
        name(k) = sname_ijts(kx)
        lname(k) = lname_ijts(kx)
        units(k) = units_ijts(kx)
        irange(k) = ir_ijt(1)
        iacc(k) = ia_ijts(kx)
        aij1(:,:,k) = taijs(:,:,kx)
        aij2(:,:,k) = 1.
        scale(k) = scale_ijts(kx)
      end do

#ifdef TRACERS_SPECIAL_O18
C****
C**** Calculations of deuterium excess (d=dD-d18O)
C****
      if (n_H2O18.gt.0 .and. n_HDO.gt.0) then
        n1=n_H2O18
        n2=n_HDO
C**** precipitation
        k=k+1
        n=tij_prec
        ijtype(k) = 3
        name(k) = "prec_ij_dex"
        lname(k) = "Deuterium excess in precip"
        units(k) = "per mil"
        irange(k) = ir_m45_130
        iacc(k) = ia_src
        iord(k) = 1
        aij1(:,:,k) = 1d3*(taijn(:,:,tij_prec,n2)/trw0(n2)-
     *                8.*taijn(:,:,tij_prec,n1)/trw0(n1)+
     *                7.*taijn(:,:,tij_prec,n_water))
        aij2(:,:,k) = taijn(:,:,tij_prec,n_water)
        scale(k) = 1.
C**** ground concentration
        k=k+1
        n=tij_grnd
        ijtype(k) = 3
        name(k) = "grnd_ij_dex"
        lname(k) = "Deuterium excess at Ground"
        units(k) = "per mil"
        irange(k) = ir_m45_130
        iacc(k) = ia_src
        iord(k) = 2
        aij1(:,:,k) = 1d3*(taijn(:,:,tij_grnd,n2)/trw0(n2)-
     *                8.*taijn(:,:,tij_grnd,n1)/trw0(n1)+
     *                7.*taijn(:,:,tij_grnd,n_water))
        aij2(:,:,k) = taijn(:,:,tij_grnd,n_water)
        scale(k) = 1.
C**** water vapour
        do l=1,lm
          k=k+1
          ijtype(k) = 3
          lname(k) = "Deuterium excess water vapour Level"
          write(lname(k)(36:38),'(I3)') l
          name(k) = "wvap_ij_dex"//adjustl(lname(k)(36:38))
          units(k) = "per mil"
          irange(k) = ir_m45_130
          iacc(k) = ia_src
          iord(k) = l+2
          aij1(:,:,k) = 1d3*(taijln(:,:,l,n2)/trw0(n2)-
     *         8.*taijln(:,:,l,n1)/trw0(n1)+
     *         7.*taijln(:,:,l,n_water))
          aij2(:,:,k) = taijln(:,:,l,n_water)
          scale(k) = 1.
        end do
      end if
#endif

!     nmaplets = ktmax   ! (lm+ktaij)*ntm+ktaijs
      nmaplets = k
      Qk(k+1:ktmax)=.false.

      xlb=acc_period(1:3)//' '//acc_period(4:12)//' '//XLABEL(1:LRUNID)
C****
      DAYS=(Itime-Itime0)/DFLOAT(nday)
C**** Collect the appropriate weight-arrays in WT_IJ
      wt_ij(:,:,:) = 1.

C**** Print out 6-map pages
      do n=1,nmaplets
        if (mod(n-1,6) .eq. 0) then
c**** print header lines
          write (6,'("1",a)') xlabel
          write (6,902) jyear0,amon0,jdate0,jhour0,
     *      jyear,amon,jdate,jhour,itime,days
        end if
        kcolmn = 1 + mod(n-1,3)
        if (kcolmn .eq. 1) line=' '
c**** Find, then display the appropriate array
        if (Iord(n) .gt. 0 .and. Qk(n)) then
          call ijt_mapk (ijtype(n),nt(n),Iord(n),aij1(1,1,n),aij2(1,1,n)
     *         ,smap,smapj,gm,jgrid,scale(n),iacc(n),irange(n),name(n)
     *         ,lname(n),units(n)) 
          title=trim(lname(n))//' ('//trim(units(n))//')'
          call maptxt(smap,smapj,gm,irange(n),title,line,kcolmn,nlines)
          if(qdiag) call pout_ij(title//xlb,name(n),lname(n),units(n),
     *                            smap,smapj,gm,jgrid)
          Qk(n) = .false.
        end if
c**** copy virtual half-page to paper if appropriate
        if (kcolmn.eq.3 .or. n.eq.nmaplets) then
          do k=1,nlines
            write (6,'(a133)') line(k)
          end do
        end if
      end do

      if (.not.qdiag) RETURN
C**** produce binary files of remaining fields if appropriate
      do n=1,ktmax
        if (Qk(n)) then
          call ijt_mapk (ijtype(n),nt(n),Iord(n),aij1(1,1,n),aij2(1,1,n)
     *         ,smap,smapj,gm,jgrid,scale(n),iacc(n),irange(n),name(n)
     *         ,lname(n),units(n)) 
          title=trim(lname(n))//' ('//trim(units(n))//')'
          call pout_ij(title//xlb,name,lname(n),units(n),smap,smapj,gm
     *         ,jgrid)
        end if
      end do
      if(qdiag) call close_ij

      RETURN
C****
  902 FORMAT ('0',15X,'From:',I6,A6,I2,',  Hr',I3,
     *  6X,'To:',I6,A6,I2,', Hr',I3,'  Model-Time:',I9,5X,
     *  'Dif:',F7.2,' Days')
      END SUBROUTINE DIAGIJt
#endif

      subroutine IJt_MAPk(nmap,n,k,aij1,aij2,smap,smapj,gm,jgrid
     *     ,scale,iacc,irange,name,lname,units)
!@sum ijt_MAPk returns the map data and related terms for the k-th field
!@+   for tracers and tracer sources/sinks
      USE CONSTANT, only: teeny
      USE MODEL_COM, only:im,jm, idacc
      USE GEOM, only: dxyp
      USE TRACER_COM
      USE DAGCOM

      IMPLICIT NONE

      REAL*8, DIMENSION(IM,JM) :: anum,adenom,smap,aij1,aij2
      REAL*8, DIMENSION(JM) :: smapj
      integer i,j,k,iwt,jgrid,irange,n,nmap,iacc
      character(len=30) name,units
      character(len=80) lname
      real*8 :: gm,nh,sh, off, byiacc,scale
!@var isumz,isumg = 1 or 2 if zon,glob sums or means are appropriate
      integer isumz,isumg,k1

      isumz = 2 ; isumg = 2  !  default: in most cases MEANS are needed
      iwt = 1 ; jgrid = 1

      adenom = 1.
      byiacc = 1./(idacc(iacc)+teeny)
c**** tracer amounts (divide by area) and Sources and sinks
      if (nmap.eq.1) then
        do j=1,jm
        do i=1,im
          anum(i,j)=aij1(i,j)*byiacc*scale/dxyp(j)
        end do
        end do
c**** tracer sums and means (no division by area)
      else if (nmap.eq.2) then
        if (index(lname,' x POICE') .gt. 0) then  ! weight by sea ice
          k1 = index(lname,' x ')
          adenom(:,:)=aij(:,:,ij_rsoi)*byiacc
          lname(k1:80)=''
        end if
        do j=1,jm
        do i=1,im
          anum(i,j)=aij1(i,j)*byiacc*scale
        end do
        end do
c**** ratios (i.e. per mil diags)
      else if (nmap.eq.3) then
        if (index(lname,' x POICE') .gt. 0) then  ! ignore weighting
          k1 = index(lname,' x ')
          lname(k1:80)=''
        end if
        do j=1,jm
        do i=1,im
          anum(i,j)=aij1(i,j)
          adenom(i,j)=aij2(i,j)
        end do
        end do
C**** PROBLEM
      else  ! should not happen
        write (6,*) 'no field defined for ijt_index',n
        call stop_model(
     &       'ijt_mapk: undefined extra ij_field for tracers',255)
      end if

c**** Find final field and zonal, global, and hemispheric means
  100 continue
      call ij_avg (anum,adenom,wt_ij(1,1,iwt),jgrid,isumz,isumg, ! in
     *             smap,smapj,gm,nh,sh)                    ! out
      return
      end subroutine ijt_mapk

#ifdef TRACERS_ON
      SUBROUTINE io_trdiag(kunit,it,iaction,ioerr)
!@sum  io_trdiag reads and writes tracer diagnostics arrays to file
!@auth Jean Lerner
!@ver  1.0
      USE MODEL_COM, only: ioread,iowrite,iowrite_mon,iowrite_single
     *     ,irerun,ioread_single,lhead
      USE TRACER_DIAG_COM, only: tacc,ktacc
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "TRdiag01"
!@var it input/ouput value of hour
      INTEGER, INTENT(INOUT) :: it
!@var TACC4(KTACC) dummy array for reading diagnostics files
      REAL*4 TACC4(KTACC)

      write (MODULE_HEADER(lhead+1:80),'(a,i8,a)')
     *   'R8 TACC(',ktacc,'),it'

      SELECT CASE (IACTION)
      CASE (IOWRITE,IOWRITE_MON) ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER, TACC,it
      CASE (IOWRITE_SINGLE)    ! output to acc file
        MODULE_HEADER(LHEAD+1:LHEAD+2) = 'R4'
        WRITE (kunit,err=10) MODULE_HEADER, REAL(TACC,KIND=4),it
      CASE (IOREAD:)          ! input from restart file
        SELECT CASE (IACTION)
        CASE (ioread_single)    ! accumulate diagnostic files
          READ (kunit,err=10) HEADER,TACC4
          TACC = TACC+TACC4     ! accumulate diagnostics
          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",HEADER
     *           ,MODULE_HEADER
            GO TO 10
          END IF
        CASE (ioread,irerun)  ! restarts
          READ (kunit,err=10) HEADER, TACC,it
          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",HEADER
     *           ,MODULE_HEADER
            GO TO 10
          end IF
        END SELECT
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_trdiag
#endif


