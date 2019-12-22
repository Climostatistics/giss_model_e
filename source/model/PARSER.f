      module PARSER
!@sum programs to read parameters from the rundeck into the database
!@auth I. Aleinov
!@ver 1.0     

      
      contains

      subroutine strip_comment( str )
      ! remove comment at the end of the line. comments symbols: !#
      implicit none
      character*(*) str
      integer n

      n = scan( str, '!#' )
      if ( n == 1 ) str = '' ! empty string
      if ( n > 1 ) str = str(:n-1)

      return
      end subroutine strip_comment
      

      subroutine skip_junk( str )
      implicit none
      character*(*) str
      integer n

      do while ( len_trim( str ) > 0 .and. scan( str, ' =,' ) == 1 )
        str = str(2:)
      enddo
      return
      end subroutine skip_junk


      subroutine sread_int( str, value )
      implicit none
      character*(*) str
      integer value
      integer n

      read ( str, * ) value

      ! remove chars till the next [ =,]
      n = scan( str, ' =,' )
      str = str(n+1:)

      call skip_junk( str )

      return
      end subroutine sread_int


      subroutine sread_real( str, value )
      implicit none
      character*(*) str
      real*8 value
      integer n

      read ( str, * ) value

      ! remove chars till the next [ =,]
      n = scan( str, ' =,' )
      str = str(n+1:)

      call skip_junk( str )

      return
      end subroutine sread_real


      subroutine sread_char( str, value )
      implicit none
      character*(*) str
      character*(*) value
      !character*256 tstr
      integer n, n1

      ! replace '=' with space if not quoted
      n1 = scan( str, '''' )
      n  = scan( str, '=' )
      if ( n>0 .and. ( n1==0 .or. n<n1 ) ) str(n:n) = ' '

      read ( str, * ) value

      if ( scan( str, '''' ) == 1 ) then  ! quated string
        str = str(2:)
        n = scan( str, '''' )
        str = str(n+1:)
      else  ! remove chars till the next [ =,]
        n = scan( str, ' =,' ) 
        str = str(n+1:)
      endif

      call skip_junk( str )

      return
      end subroutine sread_char


      subroutine parse_params( kunit )
      use PARAM
      implicit none
      integer, parameter :: MAXDIM=64
      integer, intent(in) :: kunit
      character buf(256)
      character*256 bufs
      character*256 ttt
      character*16 name
      character*1 type
      integer np, i
      integer ivars(MAXDIM)
      real*8 rvars(MAXDIM)
      character*64 cvars(MAXDIM)

      ! skip unrelated stuff
      do
        read( kunit, '(a256)', err=666, end=667 ) bufs
        if ( len_trim(bufs) < 1 ) cycle
        read( bufs, * ) name
        if ( name == '&&PARAMETERS' ) exit
      enddo

      do
        read( kunit, '(a256)', err=666, end=666 ) bufs

        if ( len_trim(bufs) < 1 ) cycle

        call strip_comment( bufs )
        call skip_junk( bufs )

        if ( len_trim(bufs) < 1 ) cycle

        !read the name of the variable
        call sread_char( bufs, name )

        if ( name == '&&END_PARAMETERS' ) exit  ! end of list 

        if ( len_trim(bufs) < 1 ) then
          print *,'PARSER: no values were given to param: ', name
          call stop_model('PARSER error',255)
        endif

        ! now check the type of variables
        if ( scan( bufs, '''' ) > 0 ) then
          type = 'c'
        else if ( scan( bufs, '.' ) > 0 ) then
          type = 'r'
        else
          type = 'i'
        endif

        select case ( type )
        case ('i')
          np = 0
          do while ( len_trim(bufs) > 0 )
            np = np+1
            call sread_int( bufs, ivars(np) )
          end do
          call set_param( name, ivars, np, 'o' )
        case ('r')
          np = 0
          do while ( len_trim(bufs) > 0 )
            np = np+1
            call sread_real( bufs, rvars(np) )
          end do
          call set_param( name, rvars, np, 'o' )
        case ('c')
          np = 0
          do while ( len_trim(bufs) > 0 )
            np = np+1
            call sread_char( bufs, cvars(np) )
          end do
          call set_param( name, cvars, np, 'o' )
        end select

      enddo

      return
 666  print *, 'PARSER: Error reading params'
      call stop_model( 'PARSER: Error reading params', 255 )
 667  print *, 'PARSER: No &&PARAMETERS or &&END_PARAMETERS found'
      call stop_model(
     &     'PARSER: No &&PARAMETERS or &&END_PARAMETERS found',255)
      end subroutine parse_params

      end module PARSER









