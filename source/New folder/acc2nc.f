c compilation:
c f90 -c -64 -mips4 -O2 acc2nc.f
c f90 -64 -mips4 -O2 MODEL_COM.o DAGCOM.o RADNCB.o NCACC.o DEFACC.o
c             GEOM_B.o acc2nc.o /usr/local/netcdf-3.4/lib64/libnetcdf.a
c f90 -64 -mips4 -O2 acc2nc.o ../model/gcmlib.a
c             /usr/local/netcdf-3.4/lib64/libnetcdf.a -o acc2nc
      program acc2nc
      use MODEL_COM
      use DAGCOM
      use GEOM
      implicit none
      character(len=80) :: accfile,outfile,cmdstr,stopstr
      INTEGER :: IARGC
      INTEGER :: IUNIT
      INTEGER :: K,L
      LOGICAL :: EX
      INTEGER :: IT,ioerr
      REAL :: DUMMY4
      CHARACTER(len=8) :: MODULE_HEADER
      integer :: k01,k02,k03,k04,k05,k06,k07,k08,k09,k10,
     &           k11,k12,k13,k14,k15,k16,k17,k18,k19
      real ::
     &     aj_odd,areg_odd,apj_odd,ajl_odd,asjl_odd,aij_odd,
     &     ail_odd,energy_odd,consrv_odd,speca_odd,
     &     atpe_odd,adiurn_odd,wave_odd,ajk_odd,aijk_odd,
     &     tsfrez_odd

      call getarg(0,cmdstr)
      if(iargc().ne.2) then
         stopstr=trim(cmdstr)//' infile outfile'
         stop trim(stopstr)
      endif
      call getarg(1,accfile)
      inquire(file=accfile,exist=ex)
      if(.not.ex) then
         stopstr='file does not exist: '//trim(accfile)
         stop trim(stopstr)
      endif
      call getarg(2,outfile)
C****
C**** READ THE SINGLE-PRECISION ACC FILE
C****

      iunit=30
      open(iunit,file=accfile,form='unformatted')
      call io_rsf(iunit,it,ioread_single,ioerr)
      if (ioerr.eq.-1) go to 840
      close(iunit)

C**** Convert ACC single to double precision

      k01 = size(aj)
      k02 = size(areg)
      k03 = size(apj)
      k04 = size(ajl)
      k05 = size(asjl)
      k06 = size(aij)
      k07 = size(ail)
      k09 = size(energy)
      k10 = size(consrv)
      k11 = size(speca)
      k12 = size(atpe)
      k13 = size(adiurn)
      k14 = size(wave)
      k15 = size(ajk)
      k16 = size(aijk)
      k19 = size(tsfrez)

      call CNV428(aj      ,aj_odd      ,k01)
      call CNV428(areg    ,areg_odd    ,k02)
      call CNV428(apj     ,apj_odd     ,k03)
      call CNV428(ajl     ,ajl_odd     ,k04)
      call CNV428(asjl    ,asjl_odd    ,k05)
      call CNV428(aij     ,aij_odd     ,k06)
      call CNV428(ail     ,ail_odd     ,k07)
      call CNV428(energy  ,energy_odd  ,k09)
      call CNV428(consrv  ,consrv_odd  ,k10)
      call CNV428(speca   ,speca_odd   ,k11)
      call CNV428(atpe    ,atpe_odd    ,k12)
      call CNV428(adiurn  ,adiurn_odd  ,k13)
      call CNV428(wave    ,wave_odd    ,k14)
      call CNV428(ajk     ,ajk_odd     ,k15)
      call CNV428(aijk    ,aijk_odd    ,k16)
      CALL CNV428(tsfrez  ,tsfrez_odd  ,k19)

      call geom_b
      call init_DIAG
      call write_nc_acc(outfile)

      stop
C****
C**** TERMINATE BECAUSE OF IMPROPER PICK-UP OR PARAMETER CONFLICTS
C****
  840 WRITE(6,*) 'READ ERROR ON RESTART FILE TAPE'
c     IF(IM.NE.JC(1).OR.JM.NE.JC(2).OR.LM.NE.JC(3)) then
c        WRITE(6,*) 'IM,JM,LM mismatch: ',IM,JM,LM,' vs ',JC(1:3)
c     ENDIF
c     IF(KACC.NE.JC(5)) WRITE(6,*) 'KACC mismatch: ',KACC,' vs ',JC(5)
      STOP 3
c 850 WRITE(6,'('' END OF TAPE REACHED...LAST TAU ON TAPE '',I6)')IT
c     STOP 3

      end program acc2nc

      subroutine cnv428(arr,arr_odd,n)
      implicit none
      integer :: n
      double precision, dimension(n) :: arr
      real*4 :: arr_odd
      integer :: i
      if(mod(n,2).eq.1) arr(n)=dble(arr_odd)
      do i=n/2,1,-1
         arr(2*i-1:2*i)=dble(transfer(arr(i),arr_odd,size=2))
      enddo
      return
      end subroutine cnv428
