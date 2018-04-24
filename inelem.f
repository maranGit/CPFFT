c     ****************************************************************
c     *                                                              *
c     *                      subroutine inelem                       *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *     this subroutine supervises and conducts the input of     *
c     *     element type and properties.                             * 
c     *                                                              *
c     ****************************************************************
c
      subroutine inelem()
      use fft, only: matList, mat_props
      implicit none
      include 'common.main'

c                 global

c                 local
      integer, allocatable, dimension(:) :: intlst
      integer :: lenlst, errnum
      integer :: dum
      character :: dums*8
      real :: dumr
      double precision :: dumd

      allocate( intlst(mxlsz) )
c
c                       translate the list of elements input on
c                       this line.
c     
      errnum = 0
      do while (.true.)
        call readsc()
        call scan()
        call trlist(intlst,mxlsz,noelem,lenlst,errnum)
        if( errnum.eq.2 .or. errnum.eq.3 ) then
          call errmsg(9,dum,dums,dumr,dumd)
          cycle
        endif
        if( errnum.eq.4 ) exit

        ! successfully read an element list
        call backsp(1)
        call ElemProps(intlst,lenlst)
      enddo
      call backsp(1)

      deallocate( intlst )
      return
c
      contains
c     ****************************************************************
c     *                      subroutine ElemProps                    *
c     *                   read element properties                    * 
c     ****************************************************************
      subroutine ElemProps(intlst,lenlst)
      implicit none

c                         global
      integer :: intlst(*), lenlst

c                         local
      integer :: nc, currMat, iplist, icn, elem
      character :: unknown*24,lname*80,mname*24
      logical, external :: endcrd, matchs, label, scanms
      logical :: matFound
c
      do while ( .not. endcrd(dum) )

        ! detect and store material
        if( matchs('material',8) ) then

          if( .not. label(dum) ) then
            call errmsg(11,dum,dums,dumr,dumd)
            cycle
          end if

          ! search stored material array for current material
          lname  = ' '
          mname = ' '
          call entits( lname, nc )
          if( nc .gt. 24 ) nc = 24
          mname(1:nc) = lname(1:nc)
          
          currMat = 1
          matFound = .false.
          do while( mat_props(currMat)%assigned )
            if( scanms( mat_props(currMat)%matnam, mname, 24 ) ) then
              matFound = .true.
              exit
            end if
            currMat = currMat + 1
          end do
          if( .not. matFound ) then
            call errmsg(12,dum,mname,dumr,dumd)
            return
          end if

          ! store in matList with current material ID
          icn    = 0
          iplist = 1
          do while( .true. )
            call trxlst(intlst,lenlst,iplist,icn,elem)
            if( iplist.eq.0 ) exit
c           matList(elem) = currMat
          end do

        ! unknown terms in element properties
        else
          call entits(unknown,nc)
          call errmsg(10,dum,unknown(1:nc),dumr,dumd)
        end if
      end do
      return
      end subroutine ! ElemProps()
      end subroutine ! inelem()
