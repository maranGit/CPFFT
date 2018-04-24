c     ****************************************************************
c     *                                                              *
c     *                     subroutine trlist                        *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/26/2011                 *
c     *                                                              *
c     *     replacement for trlist to support user-defined lists     *  
c     *                                                              *
c     ****************************************************************
c
      subroutine trlist( list, mlist, iall, nlist, ierr )
c
c          scan action to input a list of integer terms.
c          the input can be a conventional integerlist or
c          a string containing the name of a previously defined 
c          user list.
c
c          see trscan_list for all details on a conventional
c          integerlist.
c
c          dummy arguments
c              list      (output) - list of parsed input - as described
c              mlist     (input)  - allowable size of list
c              iall      (input)  - value of 'all'
c                                   = 0 - 'all' is not acceptable
c              nlist     (output) - number of terms stored in list
c              ierr      (output) - error code
c                                   = 1 - no error
c                                   = 2 - parse rules failded
c                                   = 3 - list overflow
c                                   = 4 - list not found
c
c
c          parsing rules:
c           - on entry, the calling routine has made the current scan entity
c             either the start of a conventional integer list or a string
c           - on exit we must have scan entity be the next item on line 
c             after the string. this is compatible with the processing of
c             conventional integerlists. 
c           - this routine does not touch the internal scan logical
c             flag tracking tru/false tests. here we do not know
c             what the user code expects so we leave exactly like simple
c             intergerlist
c
c     use main_data, only : user_lists
c
c
c                 used defined, named lists of integers.
c                 usually node numbers or element numbers
c
c
      implicit integer(a-z)
      type :: ulist
        character(len=24) :: name
        integer :: length_list
        integer, allocatable, dimension(:) :: list
      end type
c
      type (ulist), dimension(100) :: user_lists ! 100 is set in param_def  
      include 'param_def'
      dimension list(*)
      character lname*24, name*80
      logical isstring, scanms, debug
      data debug / .false. /
c
c          if we have a regular <integerlist> just process as before and
c          return
c
      call trscan_list( list, mlist, iall, nlist, ierr )
      if( ierr .ne. 4 ) return
c
c          <integerlist> not found. check for user defined list name 
c          in a string. scan entity is
c          already the string if it is there. use a scan function that
c          does not touch the internal scanner flag (next)
c
      if( .not. isstring(idummy) ) return
c
c          user defined list. find list in table, insert into list 
c          space passed in. advance scanner to next entity on return to 
c          match behavior of conventional integerlist.

      lname(1:24) = ' '; call entits( name, nchars )
      if( nchars > 24 ) nchars = 24; lname(1:nchars) = name(1:nchars)
      if( debug )  write(*,*) "... list id: ", lname
c      
      list_col  = 0
      do i = 1, max_user_lists
       if( scanms( user_lists(i)%name, lname, 24 ) ) then
         list_col = i; go to 100
       end if
      end do 
c
      ierr = 4
      return
c
c          list found. check for overflow of space provided for
c          list. extract values from stored lists and return.
c          put the next line entity into the scanner.
c
 100  continue
      stored_length = user_lists(i)%length_list
      if( stored_length == 0 ) then
         ierr = 2
         call ulist_error( 27 )
         call scan
         return
      end if 
      if( mlist < stored_length ) then
        ierr = 3; return
      end if
      nlist = stored_length
      list(1:nlist) = user_lists(list_col)%list(1:nlist)
      ierr = 1
      if( debug ) then
        write(*,*) '.. list_col, stored_length ', list_col, nlist
        write(*,*) 'list: ', list(1:nlist)
      end if 
      call scan
      return
c
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine ulist_error                  *
c     *                                                              *
c     *                       written by : dw and rhd                *
c     *                                                              *
c     *      error message routine for user-define integer lists     *
c     *                                                              *
c     ****************************************************************
c
      subroutine ulist_error( message )
c     use main_data, only : user_lists
      implicit integer(a-z)
      include 'common.main'
      type :: ulist
        character(len=24) :: name
        integer :: length_list
        integer, allocatable, dimension(:) :: list
      end type
c
      type (ulist), dimension(100) :: user_lists ! 100 is set in param_def  
      character(len=80) :: string
c
      select case( message )
      case( 27 )
         call entits( string, strlng )
         write(out,9027) string(1:strlng)
      case default
        write(out,9999)
        stop
      end select
c
      return
c
 9027 format(/1x,'>>>>> warning: the user-list named: ',a,
     & /16x,'has no entries...')
 9999 format(/1x,'>>>>> Fatal Error: routine ulist_error.',
     &   /16x,   'should have not reach this point.')
c
      end
