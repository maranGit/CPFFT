c     ****************************************************************
c     *                                                              *
c     *                      subroutine inlodcase                    *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                      Read one load case                      *
c     *                                                              *
c     ****************************************************************
c
      subroutine inlodcase
      use fft, only: FP_max, isNBC
      implicit none
c
      integer :: dum
      real :: dumr
      real(8) :: dumd
      character :: dums*8
      real(8), parameter :: zero = 0.0D0
      logical, external :: matchs, numd
c
      FP_max = zero
c
      do while ( .true. )
        call readsc()
c
c                 read essential boundary condition
c
        if     ( matchs('F_xx',4) ) then
          if ( .not. numd( FP_max(1) ) ) then
            call errmsg(13,dum,dums,dumr,dumd)
          else
            isNBC(1) = .false.
          end if
        elseif ( matchs('F_xy',4) ) then
          if ( .not. numd( FP_max(2) ) ) then
            call errmsg(13,dum,dums,dumr,dumd)
          else
            isNBC(2) = .false.
          end if
        elseif ( matchs('F_xz',4) ) then
          if ( .not. numd( FP_max(3) ) ) then
            call errmsg(13,dum,dums,dumr,dumd)
          else
            isNBC(3) = .false.
          end if
        elseif ( matchs('F_yx',4) ) then
          if ( .not. numd( FP_max(4) ) ) then
            call errmsg(13,dum,dums,dumr,dumd)
          else
            isNBC(4) = .false.
          end if
        elseif ( matchs('F_yy',4) ) then
          if ( .not. numd( FP_max(5) ) ) then
            call errmsg(13,dum,dums,dumr,dumd)
          else
            isNBC(5) = .false.
          end if
        elseif ( matchs('F_yz',4) ) then
          if ( .not. numd( FP_max(6) ) ) then
            call errmsg(13,dum,dums,dumr,dumd)
          else
            isNBC(6) = .false.
          end if
        elseif ( matchs('F_zx',4) ) then
          if ( .not. numd( FP_max(7) ) ) then
            call errmsg(13,dum,dums,dumr,dumd)
          else
            isNBC(7) = .false.
          end if
        elseif ( matchs('F_zy',4) ) then
          if ( .not. numd( FP_max(8) ) ) then
            call errmsg(13,dum,dums,dumr,dumd)
          else
            isNBC(8) = .false.
          end if
        elseif ( matchs('F_zz',4) ) then
          if ( .not. numd( FP_max(9) ) ) then
            call errmsg(13,dum,dums,dumr,dumd)
          else
            isNBC(9) = .false.
          end if
c
c                 read natural boundary condition
c
        elseif ( matchs('P_xx',4) ) then
          if ( .not. numd( FP_max(1) ) ) then
            call errmsg(13,dum,dums,dumr,dumd)
          else
            isNBC(1) = .true.
          end if
        elseif ( matchs('P_xy',4) ) then
          if ( .not. numd( FP_max(2) ) ) then
            call errmsg(13,dum,dums,dumr,dumd)
          else
            isNBC(2) = .true.
          end if
        elseif ( matchs('P_xz',4) ) then
          if ( .not. numd( FP_max(3) ) ) then
            call errmsg(13,dum,dums,dumr,dumd)
          else
            isNBC(3) = .true.
          end if
        elseif ( matchs('P_yx',4) ) then
          if ( .not. numd( FP_max(4) ) ) then
            call errmsg(13,dum,dums,dumr,dumd)
          else
            isNBC(4) = .true.
          end if
        elseif ( matchs('P_yy',4) ) then
          if ( .not. numd( FP_max(5) ) ) then
            call errmsg(13,dum,dums,dumr,dumd)
          else
            isNBC(5) = .true.
          end if
        elseif ( matchs('P_yz',4) ) then
          if ( .not. numd( FP_max(6) ) ) then
            call errmsg(13,dum,dums,dumr,dumd)
          else
            isNBC(6) = .true.
          end if
        elseif ( matchs('P_zx',4) ) then
          if ( .not. numd( FP_max(7) ) ) then
            call errmsg(13,dum,dums,dumr,dumd)
          else
            isNBC(7) = .true.
          end if
        elseif ( matchs('P_zy',4) ) then
          if ( .not. numd( FP_max(8) ) ) then
            call errmsg(13,dum,dums,dumr,dumd)
          else
            isNBC(8) = .true.
          end if
        elseif ( matchs('P_zz',4) ) then
          if ( .not. numd( FP_max(9) ) ) then
            call errmsg(13,dum,dums,dumr,dumd)
          else
            isNBC(9) = .true.
          end if
c
c                 higher level command
c
        else
          exit
        end if
      end do
      return
      end subroutine
