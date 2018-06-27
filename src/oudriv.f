c     ****************************************************************
c     *                                                              *
c     *                      subroutine oudrive                      *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified : 4/25/2018 RM               *
c     *                                                              *
c     *     drive output of any and all quantities requested by the  *
c     *     user. all phases of output except output of residual     *
c     *     loads are driven from this section of code.              *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine oudrive()
      implicit none
c                          local
      logical, external :: matchs_exact, matchs

      if ( matchs_exact('data') ) call DataOut()
      if ( matchs_exact('model') ) call ModelOut()
      if ( matchs('results',6) ) call ResultOut()

      return
      end subroutine
c
c     ****************************************************************
c     *                      subroutine DataOut                      *
c     *                       written by : RM                        *
c     *                      write input data                        *
c     ****************************************************************
      subroutine DataOut()
      use fft
      implicit none
      include 'common.main'
c                          local
      integer :: ii

c     number of grids
      write(out,1001) N

c     material parameters
      ii = 1
      do while ( mat_assigned(ii) )
        write(out,1002) ii, matnam(ii)
        ii = ii + 1
      end do

c     dimension
      write(out,1003) l_x, l_y, l_z

c     element definition
      do ii = 1, N3
        write(out,1004) ii, matList(ii)
      end do

c     total strains
      do ii = 1, 9
        write(out,1005) ii, FP_max(ii)
      end do

c     constraints of each step
      do ii = 1, nstep
        write(out,1006) ii, mults(ii)
      enddo

c     nonlinear analysis parameters
      return
 1001 format(1x,'Number of grid points: ',i5)
 1002 format(1x,'Name of the ',i5,'th material: ',a24)
 1003 format(1x,'size of x direction: ',i5,
     &  '; y direction: ',i5,'; z direction: ',i5)
 1004 format(1x,'Element ',i5,', material ',i5)
 1005 format(1x,'F',i5,', value ',f7.2)
 1006 format(1x,'Step ',i5,', constraint ',f7.2)
      end subroutine
c
c     ****************************************************************
c     *                      subroutine ResultOut                    *
c     *                       written by : RM                        *
c     *                      write result data                       *
c     ****************************************************************
      subroutine ResultOut()
      use fft, only: out_step
      implicit none
      include 'common.main'
c
c                             local
c
      integer, allocatable, dimension(:) :: intlst
      integer   :: iplist, icn, step, param
      integer   :: lenlst, errnum
      integer   :: dum, dummy
      real      :: dumr
      real(8)   :: dumd
      character :: dums*8
      logical, external :: matchs
      logical :: debug

      debug = .false.
c
c                      read step list
c
      if( .not. matchs('steps',4) ) then
        call errmsg()
        go to 9999
      end if

      allocate( intlst(mxstep) )
      call scan()
      call trlist(intlst,mxlsz,mxstep,lenlst,errnum)

      if    ( errnum .eq. 2 ) then ! syntax error
        call errmsg(14,dum,dums,dumr,dumd)
        go to 9999
      elseif( errnum .eq. 3 ) then ! list overflow
        call errmsg(14,dum,dums,dumr,dumd)
        go to 9999
      elseif( errnum .eq. 4 ) then ! no integer list was found
        call errmsg(14,dum,dums,dumr,dumd)
        go to 9999
      end if
c
c          errnum = 1, successfully read an element list
c                  store steps in output list
c
      call backsp(1)
      out_step(1:mxstep) = .false.
c                     store constraints in step list
      iplist = 1
      icn = 0
      do while ( iplist .ne. 0 )

        call trxlst(intlst,lenlst,iplist,icn,step)

        if( step .ge. mxstep ) then
          param = mxstep
          call errmsg(16,param,dums,dumr,dumd) ! too large step
          go to 9999
        end if
c
        if( step .le. 0 ) then
          param = step
          call errmsg(16,param,dums,dumr,dumd) ! negative step
          go to 9999
        end if

        out_step(step) = .true.

      end do

      if(debug) then
        write(out,9002)
        do dum = 1, 10
          write(out,9001) dum, out_step(dum)
        end do
      end if
c
      deallocate( intlst )
c
 9999 return
 9002 format(1x,'Now checking output results',/)
 9001 format(1x,'IO status of step', i3, ':', l2,/)

      end subroutine 
