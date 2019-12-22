!@sum This file contains architecture specific code for SGI, IBM, Linux, DEC
      MODULE RANDOM
!@sum   RANDOM generates random numbers: 0<RANDom_nUmber<1
!@auth  Reto Ruedy
!@ver   1.0 (SGI,IBM,Linux,DEC)
!@cont  RANDU, RINIT, RFINAL
      IMPLICIT NONE
      INTEGER, SAVE :: IX            !@var IX     random number seed

      CONTAINS

#if defined(MACHINE_SGI) || defined(MACHINE_Linux) || defined(MACHINE_DEC)
      FUNCTION RANDU (X)
!@sum   RANDU calculates a random number based on the seed IX
!@calls RAN
      REAL*8 X                       !@var X      dummy variable
      REAL*4 RAN                     !@fun RAN    SGI intrinsic func.
      REAL*8 :: RANDU                !@var RANDU  random number
      RANDU=RAN(IX)
      RETURN
      END FUNCTION RANDU
#elif defined( MACHINE_IBM )
      FUNCTION RANDU (X)
!@sum   RANDU calculates a random number based on the seed IX
      REAL*8 X                       !@var X      dummy variable
      REAL*8 :: RANDU                !@var RANDU  random number
      INTEGER :: IY                  !@var IY     dummy integer
   10 IY=IX*65539
      SELECT CASE (IY)
      CASE (:-1)
         IY=(IY+2147483647)+1
      CASE (0)
         IX=1
         GO TO 10
      END SELECT
      IX=IY
      RANDU=DFLOAT(IY)*.465661287308D-9
      RETURN
      END FUNCTION RANDU
#else
      None of supported architectures was specified.
      This will crash the compiling process.
#endif

      SUBROUTINE RINIT (INIT)
!@sum   RINIT sets the initial seed IX
      INTEGER, INTENT(IN)  :: INIT   !@var INIT   first random no. seed
      IX=INIT
      RETURN
      END SUBROUTINE RINIT

      SUBROUTINE RFINAL (IFINAL)
!@sum   RFINAL retrieves final seed value
      INTEGER, INTENT(OUT) :: IFINAL !@var IFINAL last random no. seed
      IFINAL=IX
      RETURN
      END SUBROUTINE RFINAL

      END MODULE RANDOM

#if defined( MACHINE_SGI ) || defined( MACHINE_IBM )
      SUBROUTINE GETTIME (MNOW)
!@sum  GETTIME returns current CPU time
!@auth Gary Russell
!@ver  1.0 (SGI, IBM)
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: MNOW !@var MNOW current CPU time (.01 s)
      INTEGER :: MCLOCK            !@var MCLOCK intrinsic function
      MNOW = MCLOCK()
      RETURN
      END SUBROUTINE GETTIME
#elif defined( MACHINE_Linux )
      SUBROUTINE GETTIME (MNOW)
!@sum  GETTIME returns current CPU time
!@auth Gary Russell
!@ver  1.0 (Absoft version)
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: MNOW !@var MNOW current CPU time (.01 s)
      REAL*4 :: ETIME, TARR(2)       !@var ETIME intrinsic function
      MNOW = NINT(ETIME(TARR)*100.)
      RETURN
      END SUBROUTINE GETTIME
#elif defined( MACHINE_DEC )
      SUBROUTINE GETTIME (MNOW)
!@sum  GETTIME returns current CPU time
!@auth RIck Healy
!@ver  1.0 (DEC version)
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: MNOW !@var MNOW current CPU time (.01 s)
      REAL*4 :: SECNDS             !@var SECNDS intrinsic function
      MNOW = INT(100*SECNDS(0.0))
      RETURN
      END SUBROUTINE GETTIME
#else
      None of supported architectures was specified.
      This will crash the compiling process.
#endif

      SUBROUTINE exit_rc (code)
!@sum  exit_rc stops the run and sets a return code
!@auth Reto A Ruedy
!@ver  1.0 (SGI,IBM,Linux,DEC)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: code !@var code return code set by user
#if defined(MACHINE_SGI) || defined(MACHINE_Linux) || defined(MACHINE_DEC)
      call exit(code) !!! should check if it works for Absoft and DEC
#elif defined( MACHINE_IBM )
      call exit_(code)
#else
      None of supported architectures was specified.
      This will crash the compiling process.
#endif
      RETURN
      END SUBROUTINE exit_rc

      SUBROUTINE sys_flush (unit)
!@sum system call to flush corresponding I/O unit
!@auth I. Aleinov
!@ver  1.0 (SGI,IBM,Linux,DEC)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: unit !@var unit 
#if defined(MACHINE_SGI) || defined(MACHINE_Linux) || defined(MACHINE_DEC)
      call flush(unit) !!! should check if it works for Absoft and DEC
#elif defined( MACHINE_IBM )
      call flush_(unit)
#else
      None of supported architectures was specified.
      This will crash the compiling process.
#endif
      RETURN
      END SUBROUTINE sys_flush

      SUBROUTINE sys_signal (sig, prog)
!@sum system call to "signal"
!@auth I. Aleinov
!@ver  1.0 (SGI,IBM,Linux,Dec) !! should check if works with DEC !!
      IMPLICIT NONE
!@var unit signal number to catch
      INTEGER, INTENT(IN) :: sig
!@var prog handler subroutine for given signal
      EXTERNAL prog
#if defined(MACHINE_SGI) || defined(MACHINE_Linux) || defined(MACHINE_DEC)
      call signal( sig, prog, -1 ) 
#elif defined( MACHINE_IBM )
      call signal( sig, prog )
#else
      None of supported architectures was specified.
      This will crash the compiling process.
#endif
      RETURN
      END SUBROUTINE sys_signal
