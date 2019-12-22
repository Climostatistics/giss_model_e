      PROGRAM mkdeep
!@sum  mkdeep creates the mixed layer temperature climatology needed of
!@+    deep ocean diffusion runs
!@auth G. Russell/L. Nazarenko/G. Schmidt
C**** Note this must be compiled using the deep ocean gcmlib
      USE MODEL_COM, only : ioread,Itime,im,jm,amonth
      USE DAGCOM, only : ij_tgo2,aij
      USE FILEMANAGER, only : closeunit,openunit
      IMPLICIT NONE
      CHARACTER TITLE*80,FNAME*60,RUNID*20
      REAL*8 TG3M(IM,JM,12) !! ,TOCEANS(3,IM,JM),RSIS(IM,JM),MSIS(IM,JM)
      INTEGER M,MONACC(12),K,MTG3,ioerr,IYR,IYRE,IYRI,NFILES,IARGC
     *     ,iu_oda,lf

C**** ZERO OUT ACCUMULATING ARRAYS
      TG3M=0.
C**** Initiallise AIJ array (needed due to the way io_oda works)
      IJ_TGO2 = 1
      AIJ(:,:,1) = 0.

C**** LOOP OVER INPUT FILES
      NFILES=IARGC()
      IF (NFILES.le.0) THEN
        PRINT*,"mkdeep: make climatology to initialise deep ocean runs"
        PRINT*,"Usage: mkdeep *.oda.*"
        STOP
      END IF

      IYRE=0
      IYRI=9999
      monacc=0

C**** loop over input files
      DO K=1,NFILES
        CALL GETARG(K,FNAME)
C**** FNAME = monYYYY.odaRUNID
C**** get relevant month
        MTG3=0
        DO M=1,12
          IF (FNAME(1:3).eq.AMONTH(M)(1:3)) THEN
            monacc(M)=monacc(M)+1
            MTG3=M
          END IF
        END DO
        if (MTG3.eq.0) then
          PRINT*,trim(FNAME)," non-standard, month=",FNAME(1:3)
          stop
        end if
C**** years
        READ(FNAME(4:7),'(I4)') IYR
        IYRE=MAX(IYRE,IYR)
        IYRI=MIN(IYRI,IYR)

        call openunit(FNAME,iu_oda,.true.,.true.)
        call io_oda(iu_oda,itime,ioread,ioerr)
        if (ioerr.eq.1) goto 910
        call closeunit(iu_oda)

C**** ACCUMULATE TG3M SEPARATELY FOR EACH MONTH
        TG3M(:,:,MTG3)=TG3M(:,:,MTG3)+AIJ(:,:,IJ_TGO2)

      END DO
      RUNID=FNAME(12:len_trim(FNAME))

C**** DIVIDE EACH MONTH BY NUMBER OF ACCUMULATIONS
      DO M=1,12
        TG3M(:,:,M)=TG3M(:,:,M)/monacc(M)
      END DO

      TITLE="mean daily sums of TG3 for each month, TG3M"
      FNAME="TG3M."//trim(RUNID)
      lf=len_trim(FNAME)
      WRITE(FNAME(lf+1:lf+10),'(a1,I4.0,a1,I4.0)') '.',IYRI,'-',IYRE

      call openunit(FNAME,iu_oda,.true.,.false.)
      WRITE(iu_oda) TITLE,TG3M
      call closeunit(iu_oda)

      STOP
  910 WRITE(6,*) 'READ ERROR ON FILE ',FNAME
      STOP
      END
