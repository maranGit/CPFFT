c     ****************************************************************
c     *                                                              *
c     *                      subroutine inlod                        *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *              translate and store loading definitions         *
c     *                                                              *
c     ****************************************************************
c
      subroutine inlod()
      implicit none
      include 'param_def'
c                         global
c                         local
      integer, allocatable, dimension(:) :: intlst
      integer :: lenlst, errnum
      integer :: dum
      character :: dums*8
      real :: dumr
      double precision :: dumd
      double precision, parameter :: zero = 0.0D0
      logical, external :: matchs

      allocate( intlst(mxstep) )
      mults = zero

      do while ( .true. )
        call readsc()
        if ( .not. matchs('step',4) ) exit

        ! read step list
        call scan()
        call trlist(intlst,mxlsz,noelem,lenlst,errnum)
        if    ( errnum .eq. 2 ) then
          call errmsg(14,dum,dums,dumr,dumd)
          cycle
        elseif( errnum .eq. 3 ) then
          call errmsg(14,dum,dums,dumr,dumd)
          cycle
        elseif( errnum .eq. 4 ) then
          call errmsg(14,dum,dums,dumr,dumd)
          cycle
        else
          ! successfully read an element list
          call backsp(1)
          call LoadProps(intlst,lenlst)
        endif
      end do

      deallocate( intlst )
      return

      contains
c     ****************************************************************
c     *                      subroutine LoadProps                    *
c     *                   read loading properties                    * 
c     ****************************************************************
      subroutine LoadProps
      use fft, only: mults
      implicit none
c                         global
c                         local
      integer :: iplist, icn, step, param
      double precision :: constraints
      logical, external :: numd

c                      read constraints
      if ( matchs( 'constraints',6 ) ) call splunj
      if ( .not. numd( constraints ) ) then
        call errmsg()
        return
      end if

c                     store constraints in step list
      iplist = 1
      icn = 0
      do while ( iplist .ne. 0 )
        call trxlst(intlst,lenlst,iplist,icn,step)

        if( step .gt. mxstep ) then
          param = mxstep
          call errmsg(61,param,dums,dumr,dumd)
          exit
        end if
c
        if( step .lt. 0 ) then
          param = step
          call errmsg(69,param,dums,dumr,dumd)
          exit
        end if

        mults(step) = constraints
      end do
      end subroutine ! LoadProps
      end subroutine ! inlod
