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
      logical, external :: matchs_exact

      if ( matchs_exact('data') ) call DataOut()
      if ( matchs_exact('model') ) call ModelOut()
      if ( matchs_exact('result') ) call ResultOut()

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
      do while ( mat_props(ii)%assigned )
        write(out,1002) ii, mat_props(ii)%matnam
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
        write(out,1005) ii, F_total(ii)
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
      implicit none

      return
      end subroutine  
