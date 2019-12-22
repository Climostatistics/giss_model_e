      module param
!@sum  param does all the bookkeeping for input parameters
!@auth I. Aleinov
!@ver 1.0
!@usage This module has the following public subroutines. Most of them
!@+ are Fortramn 90 interfaces, i.e. they recognize the type of the
!@+ variables implicitly and call corresponding subroutine.
!@+
!@+ Simple copy routines to copy parameters to/from database:
!@+     set_param( name, value, dim, opt ) - put a parameter with the
!@+       name <name> and the value <value> to the database
!@+     get_param( name, value, dim ) - copy the value of the parameter
!@+       <name> from the database to the variable <value>
!@+
!@+ Query logical function:
!@+     is_set_param( name ) - returns .true. if parameter <name> is
!@+       present in the database, .false. otherwise
!@+     sync_param( name, value, dim ) - puts parameter <name> into the
!@+       database if it is not yet there, otherwise gets it
!@+
!@+ Subroutines to work with pointers:
!@+     alloc_param( name, pvalue, initval, dim ) - allocates space in
!@+       the database for the parameter <name>, fills it with data
!@+       provided in <initval> and returns pointer to it in <pvalue>
!@+     get_pparam( name, pvalue, dim ) - returns pointer <pvalue> to
!@+       the parameter data of the parameter <name> in the database
!@+
!@+ Reading/writing subroutines:
!@+     read_param( kunit, ovrwrt ) - reads the parameter database from
!@+       the unit <kunit>
!@+     write_param( kunit ) - writes the parameter database to the unit
!@+       <kunit>
!@+
!@+ Other useful subroutines:
!@+     print_param( kunit ) - does formatted output to the unit <kunit>
!@+       in a way similar to namelist
!@+     query_param( n, name, dim, ptype ) - returns information about
!@+       the parameter by its number in the database <n>. It returns
!@+       'EMPTY' in the <name> if no parameter with such an <n> exists
!@+       (i.e. if <n> is bigger than the number of parameters)
!@+
!@+ The formal arguments in the subroutines are:
!@+     name - character*(*) - the name of the parameter which is a
!@+       character string no longer than 32 bytes
!@+     value - a scalar variable or a linear array of type: integer,
!@+       real*8,character*1 to character*64
!@+     dim - integer - dimension of an array; omit 'dim' for scalars
!@+     opt - character*1 - an optional "option" (opt='o' means
!@+       "overwrite")
!@+     kunut - integer - unit number for reading/writing
!@+     ptype - character*1, intent(out) - returns the type of the
!@+       parameter: 'i' for integer, 'r' for real*8, 'c' for character
!@+     ovrwrt - logical - if .true. then reading overwrites those
!@+       parameters that are already in the database. If .false. then
!@+       those parameters which are already in the database are left
!@+       unchanged and only new parameters are added.
!@+
!@+   Read FAQ's for the full description.
!@+
!@+               CHANGE LOG:
!@+ 04/18/02 added 3 bytes to ParamStr so that its size is
!@+ divisible by 4 (needed for portability SGI,LINUX <-> IBM,COMPAQ).
!@+ Header renamed to "PARAM02 "

      implicit none
      save
      private

      public set_param, get_param, get_pparam, read_param, write_param
      public is_set_param, alloc_param, sync_param, print_param
      public query_param

      integer, parameter :: MAX_PARAMS = 256
      integer, parameter :: MAX_RPARAMS = 128
      integer, parameter :: MAX_IPARAMS = 128
      integer, parameter :: MAX_CPARAMS = 64
      integer, parameter :: MAX_NAME_LEN = 32
      integer, parameter :: MAX_CHAR_LEN = 64

      character*80 :: MODULE_HEADER='PARAM02 '

      type ParamStr
        character(MAX_NAME_LEN) name  ! parameter name
        integer indx                 ! storage for its value
        integer dim                   ! number of elements
        character(1) attrib            ! type: real ('r') or int ('i')
        character(1) reserved(3)
      end type ParamStr

      type (ParamStr), target :: Params(MAX_PARAMS)
      real*8, target :: Rdata(MAX_RPARAMS)
      integer, target :: Idata(MAX_IPARAMS)
      character*(MAX_CHAR_LEN), target :: Cdata(MAX_CPARAMS)
      integer :: num_param = 0
      integer :: num_rparam = 0
      integer :: num_iparam = 0
      integer :: num_cparam = 0

      interface set_param
       module procedure set_iparam, set_rparam, set_cparam
       module procedure  set_aiparam, set_arparam, set_acparam
      end interface

      interface get_param
       module procedure get_iparam, get_rparam, get_cparam
       module procedure get_aiparam, get_arparam, get_acparam
      end interface

      interface get_pparam
       module procedure get_piparam, get_prparam, get_pcparam
       module procedure get_paiparam, get_parparam, get_pacparam
      end interface

      interface alloc_param
      module procedure alloc_iparam, alloc_rparam, alloc_cparam
      module procedure alloc_aiparam, alloc_arparam, alloc_acparam
      end interface

      interface sync_param
      module procedure sync_iparam, sync_rparam, sync_cparam
      module procedure sync_aiparam, sync_arparam, sync_acparam
      end interface

      contains

      function is_set_param( name_in )
      implicit none
      character*(*), intent(in) :: name_in
      logical is_set_param
      integer n
      character*(MAX_NAME_LEN) name

      name = name_in
      call lowcase( name )

      do n=1,num_param
        if ( Params(n)%name == name ) exit
      enddo

      if ( n > num_param  ) then
        is_set_param = .false.
      else
        is_set_param = .true.
      endif
      return
      end function is_set_param

      subroutine set_pstr( name_in, dim, attrib, PStr, flag )
      implicit none
      character*(*), intent(in) :: name_in
      integer, intent(in) :: dim
      character*1, intent(in) ::  attrib
      type (ParamStr), pointer :: PStr
      logical, intent(in) :: flag
      integer n
      character*(MAX_NAME_LEN) name

      if ( len(name_in) > MAX_NAME_LEN ) then
        print *, 'PARAM: parameter name too long: ', name_in
        print *, 'PARAM: maximal length allowed: ', MAX_NAME_LEN
        call stop_model('PARAM: parameter name too long: ',255)
      endif

      name = name_in
      call lowcase( name )

      call get_pstr( name, dim, attrib, PStr )

      if ( associated( PStr ) ) then
        if ( .not. flag ) then
          print *, 'PARAM: attempt to set param which is already set'
          print *, 'name: ', name
          call stop_model(
     &         'PARAM: attempt to set param which is already set',255)
        else
          return  ! return PStr found by get_pstr
        endif
      endif

      if ( num_param >= MAX_PARAMS ) then
        print *, 'PARAM: Maximal number of parameters exceeded'
        print *, 'PARAM: Please recompile param with bigger MAX_PARAMS'
        call stop_model(
     &       'PARAM: Maximal number of parameters exceeded',255)
      endif

      num_param = num_param + 1
      PStr => Params(num_param)
      PStr%name = name
      PStr%attrib = attrib
      PStr%dim = dim
      PStr%reserved(1:3) = (/ ' ', ' ', ' ' /)

      select case (attrib)
        case ('i')
          if ( num_iparam+dim >= MAX_IPARAMS ) then
            print *, 'PARAM: Maximal number of int parameters exceeded'
            print *, 'PARAM: Recompile param with bigger MAX_IPARAMS'
            call stop_model(
     &           'PARAM: Maximal number of int parameters exceeded',255)
          endif
          PStr%indx = num_iparam + 1
          num_iparam = num_iparam + dim
        case ('r')
          if ( num_rparam+dim >= MAX_RPARAMS ) then
            print *, 'PARAM: Maximal number of real parameters exceeded'
            print *, 'PARAM: Recompile param with bigger MAX_RPARAMS'
            call stop_model(
     &          'PARAM: Maximal number of real parameters exceeded',255)
          endif
          PStr%indx = num_rparam + 1
          num_rparam = num_rparam + dim
        case ('c')
          if ( num_cparam+dim >= MAX_CPARAMS ) then
            print *, 'PARAM: Maximal number of char parameters exceeded'
            print *, 'PARAM: Recompile param with bigger MAX_CPARAMS'
            call stop_model(
     &          'PARAM: Maximal number of char parameters exceeded',255)
          endif
          PStr%indx = num_cparam + 1
          num_cparam = num_cparam + dim
        end select

      return
      end subroutine set_pstr


      subroutine get_pstr( name_in, dim, attrib, PStr )
      implicit none
      character*(*), intent(in) :: name_in
      integer, intent(in) :: dim
      character*1, intent(in) ::  attrib
      type (ParamStr), pointer :: PStr
      integer n
      character*(MAX_NAME_LEN) name

      name = name_in
      call lowcase( name )

      nullify( PStr )

      do n=1,num_param
        if ( Params(n)%name == name ) exit
      enddo

      if ( n > num_param  ) return   ! not found - return NULL

      if ( Params(n)%attrib /= attrib .or. Params(n)%dim /= dim ) then
        print *, 'PARAM: wrong type or dim of parameter: ', name
        print *, 'ATT: set: ', Params(n)%attrib, ' called: ', attrib
        print *, 'DIM: set: ', Params(n)%dim, ' called: ', dim
        call stop_model('PARAM: wrong type or dim of parameter',255)
      endif

      PStr => Params(n)
      return
      end subroutine get_pstr



      !***** integers ******!


      subroutine set_iparam( name, value, opt )
      implicit none
      character*(*), intent(in) :: name
      integer, intent(in) :: value
      character*1, optional, intent(in) :: opt
      integer v(1)
      v(1) = value
      call set_aiparam( name, v, 1, opt )
      return
      end subroutine set_iparam


      subroutine set_aiparam( name, value, np, opt )
      implicit none
      character*(*), intent(in) :: name
      integer, intent(in) :: np
      integer, intent(in) :: value(np)
      character*1, optional, intent(in) :: opt
      type (ParamStr), pointer :: PStr
      logical flag

      flag = .false.
      if ( present(opt) ) then
        if ( opt=='o' ) flag = .true.
      endif

      call set_pstr( name, np, 'i', PStr, flag )
      Idata( PStr%indx : PStr%indx+np-1 ) = value(1:np)
      return
      end subroutine set_aiparam


      subroutine get_iparam( name, value )
      implicit none
      character*(*), intent(in) ::  name
      integer, intent(out) ::  value
      integer v(1)

      call get_aiparam( name, v, 1 )
      value = v(1)
      return
      end subroutine get_iparam


      subroutine get_aiparam( name, value, np )
      implicit none
      character*(*), intent(in) ::  name
      integer, intent (in) :: np
      integer, intent(out) ::  value(np)
      type (ParamStr), pointer :: PStr

      call get_pstr( name, np, 'i', PStr )
      if ( .not. associated( PStr) ) then
        print *, 'PARAM: Can''t get - not in database : ', name
        call stop_model(
     &       'PARAM: Can''t get parameter - not in database',255)
      endif
      value(1:np) = Idata( PStr%indx : PStr%indx+np-1 )
      return
      end subroutine get_aiparam


      subroutine get_piparam( name, pvalue )
      implicit none
      character*(*), intent(in) ::  name
      integer, pointer ::  pvalue
      type (ParamStr), pointer :: PStr

      call get_pstr( name, 1, 'i', PStr )
      if ( .not. associated( PStr) ) then
        print *, 'PARAM: Can''t get - not in database : ', name
        call stop_model(
     &       'PARAM: Can''t get parameter - not in database',255)
      endif
      pvalue => Idata( PStr%indx )
      return
      end subroutine get_piparam

      subroutine get_paiparam( name, pvalue, np )
      implicit none
      character*(*), intent(in) ::  name
      integer, pointer ::  pvalue(:)
      integer, intent(in) :: np
      type (ParamStr), pointer :: PStr

      call get_pstr( name, np, 'i', PStr )
      if ( .not. associated( PStr) ) then
        print *, 'PARAM: Can''t get - not in database : ', name
        call stop_model(
     &       'PARAM: Can''t get parameter - not in database',255)
      endif
      pvalue => Idata( PStr%indx:PStr%indx+np-1 )
      return
      end subroutine get_paiparam



      !***** reals ******!

      subroutine set_rparam( name, value, opt )
      implicit none
      character*(*), intent(in) :: name
      real*8, intent(in) :: value
      character*1, optional, intent(in) :: opt
      real*8 v(1)
      v(1) = value
      call set_arparam( name, v, 1, opt )
      return
      end subroutine set_rparam


      subroutine set_arparam( name, value, np, opt )
      implicit none
      character*(*), intent(in) :: name
      integer, intent(in) :: np
      real*8, intent(in) :: value(np)
      character*1, optional, intent(in) :: opt
      type (ParamStr), pointer :: PStr
      logical flag

      flag = .false.
      if ( present(opt) ) then
        if ( opt=='o' ) flag = .true.
      endif

      call set_pstr( name, np, 'r', PStr, flag )
      Rdata( PStr%indx : PStr%indx+np-1 ) = value(1:np)
      return
      end subroutine set_arparam


      subroutine get_rparam( name, value )
      implicit none
      character*(*), intent(in) ::  name
      real*8, intent(out) ::  value
      real*8 v(1)

      call get_arparam( name, v, 1 )
      value = v(1)
      return
      end subroutine get_rparam


      subroutine get_arparam( name, value, np )
      implicit none
      character*(*), intent(in) ::  name
      integer, intent (in) :: np
      real*8, intent(out) ::  value(np)
      type (ParamStr), pointer :: PStr

      call get_pstr( name, np, 'r', PStr )
      if ( .not. associated( PStr) ) then
        print *, 'PARAM: Can''t get - not in database : ', name
        call stop_model(
     &       'PARAM: Can''t get parameter - not in database',255)
      endif
      value(1:np) = Rdata( PStr%indx : PStr%indx+np-1 )
      return
      end subroutine get_arparam


      subroutine get_prparam( name, pvalue )
      implicit none
      character*(*), intent(in) ::  name
      real*8, pointer ::  pvalue
      type (ParamStr), pointer :: PStr

      call get_pstr( name, 1, 'r', PStr )
      if ( .not. associated( PStr) ) then
        print *, 'PARAM: Can''t get - not in database : ', name
        call stop_model(
     &       'PARAM: Can''t get parameter - not in database',255)
      endif
      pvalue => Rdata( PStr%indx )
      return
      end subroutine get_prparam

      subroutine get_parparam( name, pvalue, np )
      implicit none
      character*(*), intent(in) ::  name
      real*8, pointer ::  pvalue(:)
      integer, intent(in) :: np
      type (ParamStr), pointer :: PStr

      call get_pstr( name, np, 'r', PStr )
      if ( .not. associated( PStr) ) then
        print *, 'PARAM: Can''t get - not in database : ', name
        call stop_model(
     &       'PARAM: Can''t get parameter - not in database',255)
      endif
      pvalue => Rdata( PStr%indx:PStr%indx+np-1 )
      return
      end subroutine get_parparam


      !***** Chars ******!

      subroutine set_cparam( name, value, opt )
      implicit none
      character*(*), intent(in) :: name
      character*(*), intent(in) :: value
      character*1, optional, intent(in) :: opt
      character*(MAX_CHAR_LEN) v(1)
      if ( len(value) > MAX_CHAR_LEN ) then
        print *, 'PARAM: Char string too long. MAX = ', MAX_CHAR_LEN
        call stop_model('PARAM: Char string too long',255)
      endif
      v(1) = value
      call set_acparam( name, v, 1, opt )
      return
      end subroutine set_cparam


      subroutine set_acparam( name, value, np, opt )
      implicit none
      character*(*), intent(in) :: name
      integer, intent(in) :: np
      character*(*), intent(in) :: value(np)
      character*1, optional, intent(in) :: opt
      type (ParamStr), pointer :: PStr
      integer n
      logical flag

      flag = .false.
      if ( present(opt) ) then
        if ( opt=='o' ) flag = .true.
      endif

      do n=1,np
        if ( len(value(n)) > MAX_CHAR_LEN ) then
          print *, 'PARAM: Char string too long. MAX = ', MAX_CHAR_LEN
          print *, 'You submitted LEN = ', len(value(n))
          call stop_model('PARAM: Char string too long',255)
        endif
      enddo
      call set_pstr( name, np, 'c', PStr, flag )
      Cdata( PStr%indx : PStr%indx+np-1 ) = value(1:np)
      return
      end subroutine set_acparam


      subroutine get_cparam( name, value )
      implicit none
      character*(*), intent(in) ::  name
      character*(*), intent(out) ::  value
      character*(MAX_CHAR_LEN) v(1)

      call get_acparam( name, v, 1 )
      value = v(1)
      return
      end subroutine get_cparam


      subroutine get_acparam( name, value, np )
      implicit none
      character*(*), intent(in) ::  name
      integer, intent (in) :: np
      character*(*), intent(out) ::  value(np)
      type (ParamStr), pointer :: PStr

      call get_pstr( name, np, 'c', PStr )
      if ( .not. associated( PStr) ) then
        print *, 'PARAM: Can''t get - not in database : ', name
        call stop_model(
     &       'PARAM: Can''t get parameter - not in database',255)
      endif
      value(1:np) = Cdata( PStr%indx : PStr%indx+np-1 )
      return
      end subroutine get_acparam


      subroutine get_pcparam( name, pvalue )
      implicit none
      character*(*), intent(in) ::  name
      character*(*), pointer ::  pvalue
      type (ParamStr), pointer :: PStr

      call get_pstr( name, 1, 'c', PStr )
      if ( .not. associated( PStr) ) then
        print *, 'PARAM: Can''t get - not in database : ', name
        call stop_model(
     &       'PARAM: Can''t get parameter - not in database',255)
      endif
      pvalue => Cdata( PStr%indx )
      return
      end subroutine get_pcparam

      subroutine get_pacparam( name, pvalue, np )
      implicit none
      character*(*), intent(in) ::  name
      character*(*), pointer ::  pvalue(:)
      integer, intent(in) :: np
      type (ParamStr), pointer :: PStr

      call get_pstr( name, np, 'c', PStr )
      if ( .not. associated( PStr) ) then
        print *, 'PARAM: Can''t get - not in database : ', name
        call stop_model(
     &       'PARAM: Can''t get parameter - not in database',255)
      endif
      pvalue => Cdata( PStr%indx:PStr%indx+np-1 )
      return
      end subroutine get_pacparam

      !***** alloc fuctions ******

      subroutine alloc_iparam( name, pvalue, initval )
      implicit none
      character*(*), intent(in) ::  name
      integer, pointer ::  pvalue
      integer, intent(in) :: initval
      call set_param( name, initval )
      call get_pparam( name, pvalue )
      end subroutine alloc_iparam


      subroutine alloc_aiparam( name, pvalue, initval, dim )
      implicit none
      character*(*), intent(in) ::  name
      integer, pointer ::  pvalue(:)
      integer, intent(in) :: initval(:)
      integer, intent(in) :: dim
      call set_param( name, initval, dim )
      call get_pparam( name, pvalue, dim )
      end subroutine alloc_aiparam


      subroutine alloc_rparam( name, pvalue, initval )
      implicit none
      character*(*), intent(in) ::  name
      real*8, pointer ::  pvalue
      real*8, intent(in) :: initval
      call set_param( name, initval )
      call get_pparam( name, pvalue )
      end subroutine alloc_rparam


      subroutine alloc_arparam( name, pvalue, initval, dim )
      implicit none
      character*(*), intent(in) ::  name
      real*8, pointer ::  pvalue(:)
      real*8, intent(in) :: initval(:)
      integer, intent(in) :: dim
      call set_param( name, initval, dim )
      call get_pparam( name, pvalue, dim )
      end subroutine alloc_arparam


      subroutine alloc_cparam( name, pvalue, initval )
      implicit none
      character*(*), intent(in) ::  name
      character*(*), pointer ::  pvalue
      character*(*), intent(in) :: initval
      call set_param( name, initval )
      call get_pparam( name, pvalue )
      end subroutine alloc_cparam


      subroutine alloc_acparam( name, pvalue, initval, dim )
      implicit none
      character*(*), intent(in) ::  name
      character*(*), pointer ::  pvalue(:)
      character*(*), intent(in) :: initval(:)
      integer, intent(in) :: dim
      call set_param( name, initval, dim )
      call get_pparam( name, pvalue, dim )
      end subroutine alloc_acparam


      !***** sync functions ******!

      subroutine sync_iparam( name, value )
      implicit none
      character*(*), intent(in) :: name
      integer, intent(inout) :: value

      if ( is_set_param( name ) ) then
        call get_param( name, value )
      else
        call set_param( name, value )
      endif
      end subroutine sync_iparam


      subroutine sync_aiparam( name, value, np )
      implicit none
      character*(*), intent(in) :: name
      integer, intent(in) :: np
      integer, intent(inout) :: value(np)

      if ( is_set_param( name ) ) then
        call get_param( name, value, np )
      else
        call set_param( name, value, np )
      endif
      end subroutine sync_aiparam


      subroutine sync_rparam( name, value )
      implicit none
      character*(*), intent(in) :: name
      real*8, intent(inout) :: value

      if ( is_set_param( name ) ) then
        call get_param( name, value )
      else
        call set_param( name, value )
      endif
      end subroutine sync_rparam


      subroutine sync_arparam( name, value, np )
      implicit none
      character*(*), intent(in) :: name
      integer, intent(in) :: np
      real*8, intent(inout) :: value(np)

      if ( is_set_param( name ) ) then
        call get_param( name, value, np )
      else
        call set_param( name, value, np )
      endif
      end subroutine sync_arparam


      subroutine sync_cparam( name, value )
      implicit none
      character*(*), intent(in) :: name
      character*(*), intent(inout) :: value

      if ( is_set_param( name ) ) then
        call get_param( name, value )
      else
        call set_param( name, value )
      endif
      end subroutine sync_cparam


      subroutine sync_acparam( name, value, np )
      implicit none
      character*(*), intent(in) :: name
      integer, intent(in) :: np
      character*(*), intent(inout) :: value(np)

      if ( is_set_param( name ) ) then
        call get_param( name, value, np )
      else
        call set_param( name, value, np )
      endif
      end subroutine sync_acparam



      !***** input / output ******!

      subroutine read_param( kunit, ovrwrt )
      implicit none
      integer, intent(in) :: kunit
      logical, intent(in) :: ovrwrt
      integer n, np
      type (ParamStr), save :: LParams(MAX_PARAMS)
      real*8, save :: LRdata(MAX_RPARAMS)
      integer, save :: LIdata(MAX_IPARAMS)
      character*(MAX_CHAR_LEN), save :: LCdata(MAX_CPARAMS)
      integer lnum_param, lnum_rparam, lnum_iparam, lnum_cparam
      character*80 HEADER

      read( kunit, err=10 ) HEADER
      backspace kunit

      if (HEADER(1:8).ne.MODULE_HEADER(1:8)) then
        if (HEADER(1:8).eq.'PARAM01 ') then
          print *, 'WARNING: PARAM: Old format of parameter data.'
          call read_param_comp01( kunit, ovrwrt )
          return
        else
          print * , 'PARAM: No parameter header in input data'
          call stop_model(
     &         'PARAM: No parameter header in input data',255)
        endif
      endif

      read( kunit, err=10 ) HEADER,
     *     lnum_param, lnum_rparam, lnum_iparam, lnum_cparam,
     *     ( LParams(n), n=1,min(lnum_param,MAX_PARAMS) ),
     *     ( LRdata(n), n=1,min(lnum_rparam,MAX_RPARAMS) ),
     *     ( LIdata(n), n=1,min(lnum_iparam,MAX_IPARAMS) ),
     *     ( LCdata(n), n=1,min(lnum_cparam,MAX_CPARAMS) )

      if (     lnum_param  > MAX_PARAMS
     *     .or. lnum_rparam > MAX_RPARAMS
     *     .or. lnum_iparam > MAX_IPARAMS
     *     .or. lnum_cparam > MAX_CPARAMS
     *     ) then
        print *, 'PARAM: parameter list in input file too long'
        print *, 'PARAM: please recompile param with bigger MAX_?PARAMS'
        print *, 'PARAM: ',num_param,num_rparam,num_iparam,num_cparam
        call stop_model(
     &       'PARAM: parameter list in input file too long',255)
      endif

      if ( lnum_param < 1 ) return   ! no parameters in the records

#ifdef MACHINE_DEC
      ! COMPAQ needs big/little endian conversion
      do n=1,lnum_param
        call swap_bytes_4( LParams(n)%indx, 1 )
        call swap_bytes_4( LParams(n)%dim,  1 )
      enddo
#endif

      ! checking big/little endian format, just in case
      if ( LParams(1)%dim > 65536 .or. LParams(1)%dim < 0 ) then
        print *, 'PARAM: wrong big/little endian format in LParams.'
        call stop_model(
     &       'PARAM: wrong big/little endian format in LParams',255)
      endif

      ! now merge the data just read with existing database
      do n=1,lnum_param
        if ( (.not. is_set_param(LParams(n)%name)) .or. ovrwrt ) then
          select case( LParams(n)%attrib )
          case ('i')
            call set_aiparam( LParams(n)%name, LIdata(LParams(n)%indx),
     *           LParams(n)%dim, 'o' )
          case ('r')
            call set_arparam( LParams(n)%name, LRdata(LParams(n)%indx),
     *           LParams(n)%dim, 'o' )
          case ('c')
            call set_acparam( LParams(n)%name, LCdata(LParams(n)%indx),
     *           LParams(n)%dim, 'o' )
          end select
        endif
      enddo
      return
 10   print *, 'PARAM: Error reading, unit = ', kunit
      call stop_model('PARAM: Error reading',255)
      end subroutine read_param


      subroutine write_param( kunit )
      implicit none
      integer, intent(in) :: kunit
      integer n

      write (MODULE_HEADER(9:80),'(i10,a)')
     *  num_param,' is the current number of parameters in database DB'

#ifdef MACHINE_DEC
      ! converting it manually to big-endian for COMPAQ compiler
      do n=1,num_param
        call swap_bytes_4( Params(n)%indx, 1 )
        call swap_bytes_4( Params(n)%dim,  1 )
      enddo
#endif

      write( kunit, err=20 ) MODULE_HEADER,
     *     num_param, num_rparam, num_iparam, num_cparam,
     *     ( Params(n), n=1,min(num_param,MAX_PARAMS) ),
     *     ( Rdata(n), n=1,min(num_rparam,MAX_RPARAMS) ),
     *     ( Idata(n), n=1,min(num_iparam,MAX_IPARAMS) ),
     *     ( Cdata(n), n=1,min(num_cparam,MAX_CPARAMS) )

#ifdef MACHINE_DEC
      ! and back to little-endian ...
      do n=1,num_param
        call swap_bytes_4( Params(n)%indx, 1 )
        call swap_bytes_4( Params(n)%dim,  1 )
      enddo
#endif
      return

 20   print *, 'PARAM: Error writing, unit = ', kunit
      call stop_model('PARAM: Error writing',255)
      end subroutine write_param

      subroutine print_param1( kunit )
      implicit none
      integer, intent(in) :: kunit
      integer, parameter :: nf = 7
      integer n, i

      write( kunit, * ) '&&PARAMETERS'
      do n=1, num_param
        select case( Params(n)%attrib )
        case ('i')
          write( kunit, '(1x,a16,a3,8i16)' )
     $         Params(n)%name, ' = ',
     $        ( Idata(Params(n)%indx+i), i=0,min(Params(n)%dim,nf)-1 )
          if ( Params(n)%dim > nf )
     $         write( kunit, '(20x,8i16)' )
     $        ( Idata(Params(n)%indx+i), i=0,Params(n)%dim-nf-1 )
        case ('r')
          write( kunit, '(1x,a16,a3,8g16.6)' )
     $         Params(n)%name, ' = ',
     $        ( Rdata(Params(n)%indx+i), i=0,min(Params(n)%dim,nf)-1 )
          if ( Params(n)%dim > nf )
     $         write( kunit, '(20x,8g16.6)' )
     $        ( Rdata(Params(n)%indx+i), i=0,Params(n)%dim-nf-1 )
        case ('c')
          write( kunit, '(1x,a16,a3,8a64)' )
     $         Params(n)%name, ' = ',
     $        ( Cdata(Params(n)%indx+i), i=0,min(Params(n)%dim,nf)-1 )
          if ( Params(n)%dim > nf )
     $         write( kunit, '(20x,8a64)' )
     $        ( Cdata(Params(n)%indx+i), i=0,Params(n)%dim-nf-1 )
        end select
      enddo
      write( kunit, * ) '&&END_PARAMETERS'
      end subroutine print_param1


      subroutine print_param( kunit )
      implicit none
      integer, intent(in) :: kunit
      integer, parameter :: nf = 7
      integer n, i

      write( kunit, * ) '&&PARAMETERS'
      do n=1, num_param
        select case( Params(n)%attrib )
        case ('i')
          write( kunit, * )
     $         trim(Params(n)%name), ' = ',
     $        ( Idata(Params(n)%indx+i), i=0,Params(n)%dim-1 )
       case ('r')
          write( kunit, * )
     $         trim(Params(n)%name), ' = ',
     $        ( Rdata(Params(n)%indx+i), i=0,Params(n)%dim-1 )
        case ('c')
          write( kunit, * )
     $         trim(Params(n)%name), ' = ',
     $        ( Cdata(Params(n)%indx+i), i=0,Params(n)%dim-1 )
        end select
      enddo
      write( kunit, * ) '&&END_PARAMETERS'
      end subroutine print_param


      subroutine query_param( n, name, dim, ptype )
      integer, intent(in) :: n
      character*(*), intent(out) :: name
      integer, intent(out) :: dim
      character*1, intent(out) :: ptype

      if ( n>0 .and. n<=num_param ) then
        name = Params(n)%name
        dim = Params(n)%dim
        ptype = Params(n)%attrib
      else
        name = 'EMPTY'
        dim = 0
        ptype = 'U'
      endif

      end subroutine query_param


      subroutine lowcase( str )
      ! converts string str to lower case
      implicit none
      character*(*) str
      integer n, i
      integer A, Z, shift, c

      A = iachar( 'A' )
      Z = iachar( 'Z' )
      shift = iachar( 'a' ) - iachar( 'A' )
      n = len_trim(str)
      do i=1,n
         c = iachar( str(i:i) )
        if ( c>=A .and. c<=Z ) str(i:i) = achar( c + shift )
      enddo
      end subroutine lowcase

!**** the code below is included for compatibility with older     ****
!**** versions; it may be removed later when not needed any more  ****

      subroutine read_param_comp01( kunit, ovrwrt )
      implicit none
      type ParamStr_comp01
        character(MAX_NAME_LEN) name ! parameter name
        integer indx                 ! storage for its value
        integer dim                   ! number of elements
        character*1 attrib            ! type: real ('r') or int ('i')
      end type ParamStr_comp01
      integer, intent(in) :: kunit
      logical, intent(in) :: ovrwrt
      integer n, np
      type (ParamStr_comp01), save :: LParams(MAX_PARAMS)
      real*8, save :: LRdata(MAX_RPARAMS)
      integer, save :: LIdata(MAX_IPARAMS)
      character*(MAX_CHAR_LEN), save :: LCdata(MAX_CPARAMS)
      integer lnum_param, lnum_rparam, lnum_iparam, lnum_cparam
      character*80 HEADER

      read( kunit, err=10 ) HEADER
      backspace kunit

      if (HEADER(1:8).ne.'PARAM01 ') then
        print *, 'PARAM: No parameter header in input data'
        call stop_model('PARAM: No parameter header in input data',255)
      endif

      read( kunit, err=10 ) HEADER,
     *     lnum_param, lnum_rparam, lnum_iparam, lnum_cparam,
     *     ( LParams(n), n=1,min(lnum_param,MAX_PARAMS) ),
     *     ( LRdata(n), n=1,min(lnum_rparam,MAX_RPARAMS) ),
     *     ( LIdata(n), n=1,min(lnum_iparam,MAX_IPARAMS) ),
     *     ( LCdata(n), n=1,min(lnum_cparam,MAX_CPARAMS) )

      if (     lnum_param  > MAX_PARAMS
     *     .or. lnum_rparam > MAX_RPARAMS
     *     .or. lnum_iparam > MAX_IPARAMS
     *     .or. lnum_cparam > MAX_CPARAMS
     *     ) then
        print *, 'PARAM: parameter list in input file too long'
        print *, 'PARAM: please recompile param with bigger MAX_?PARAMS'
        print *, 'PARAM: ',num_param,num_rparam,num_iparam,num_cparam
        call stop_model(
     &       'PARAM: parameter list in input file too long',255)
      endif

      ! now merge the data just read with existing database
      do n=1,lnum_param
        if ( (.not. is_set_param(LParams(n)%name)) .or. ovrwrt ) then
          select case( LParams(n)%attrib )
          case ('i')
            call set_aiparam( LParams(n)%name, LIdata(LParams(n)%indx),
     *           LParams(n)%dim, 'o' )
          case ('r')
            call set_arparam( LParams(n)%name, LRdata(LParams(n)%indx),
     *           LParams(n)%dim, 'o' )
          case ('c')
            call set_acparam( LParams(n)%name, LCdata(LParams(n)%indx),
     *           LParams(n)%dim, 'o' )
          end select
        endif
      enddo
      return
 10   print *, 'PARAM: Error reading, unit = ', kunit
      call stop_model('PARAM: Error reading',255)
      end subroutine read_param_comp01

      end module param

!**** this should be put somewhere else, but since it is used only ****
!**** in this module I put it here for a while ...                 ****

      subroutine swap_bytes_4( value, dim )
      integer dim
      integer(4) value(dim)
      integer(4) a
      character(1) c(4), c1, c2
      equivalence (a,c)
      integer n
!@sum  does conversion big<->little - endian for 4 byte data
!@auth I. Aleinov
!@ver 1.0
      do n=1,dim
        a = value(n)
        c1 = c(1)
        c2 = c(2)
        c(1) = c(4)
        c(2) = c(3)
        c(3) = c2
        c(4) = c1
        ! the next line is needed for stupid optimizers - do not delete
        if ( n > 16384 ) print *, a, c
        value(n) = a
      enddo
      end subroutine swap_bytes_4

