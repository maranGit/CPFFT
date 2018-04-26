c     ****************************************************************
c     *                                                              *
c     *                      subroutine indypm                       *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *     input parameters controlling how the solution is         *
c     *     performed for analysis                                   *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine indypm()
      use fft, only: maxIter, tstep, tolPCG, tolNR
      implicit none
      include 'param_def'
c                           global
c                           local
      integer :: dum
      character :: dums*24
      real :: dumr
      double precision :: dumd
      logical, external :: matchs, numd, integr, endcrd, matchs_exact

      do while ( .true. )

        call readsc()

        if ( matchs('maximum',7) ) then
          if( matchs('iteration',4) ) call splunj
          if( .not. integr(maxIter) ) 
     &      call errmsg(17,dum,'iteration',dumr,dumd)

        elseif( matchs('convergence',8) ) then
          if( matchs('tolerance',3) ) call splunj
          ! read tolerance for PCG and NR loop
          do while ( .not. endcrd(dum) )
            if( matchs_exact('CG') ) then
              ! read tolerance for PCG
              if( .not. numd(tolPCG) ) 
     &          call errmsg(17,dum,'PCG tolerance',dumr,dumd)
            elseif( matchs_exact('NR') ) then
              ! read tolerance for NR iteration
              if( .not. numd(tolNR) ) 
     &          call errmsg(17,dum,'NR tolerance',dumr,dumd)
            else
              ! Unknow option, read next line
              call errmsg(18,dum,dums,dumr,dumd)
              exit
            end if
          end do ! loop over token

        elseif( matchs('time',4) ) then
          if( matchs('step',4) ) call splunj
          ! read time step
          if( .not. numd(tstep) ) 
     &      call errmsg(17,dum,'time step',dumr,dumd)

        else
          ! find next high-level command
          exit

        end if

      end do ! loop over line
      return
      end subroutine
