      program main
      use file_info
      use fft
      implicit none
      include 'common.main'

      integer :: dum
      real :: dumr
      real(8) :: dumd
      character :: dums*8
      logical :: readnew
      logical, external :: matchs, matchs_exact, endfil

      readnew = .true.
      call readsc()
      do while ( .true. )
        if(readnew) call readsc()
        if(matchs('project',7)) then
          call ProcessInput(1)
          cycle
        elseif(matchs('number',6)) then
          call ProcessInput(2)
          call FFT_init()
          cycle
        elseif(matchs('material',8)) then
          call ProcessInput(3)
          cycle
        elseif(matchs('sizes',5)) then
          call ProcessInput(4)
          cycle
        elseif(matchs('elements',8)) then
          call ProcessInput(5)
          cycle
        elseif(matchs('blocking',8)) then
          call ProcessInput(6)
          cycle
        elseif(matchs('strains',7)) then
          call ProcessInput(7)
          cycle
        elseif(matchs('loading',7)) then
          call ProcessInput(8)
          cycle
        elseif(matchs('nonlinear',9)) then
          call ProcessInput(9)
          cycle
        elseif(matchs('compute',7)) then
          call ProcessInput(10)
          call drive_eps_sig( 1, 0 )
          call FFT_nr3()
          cycle
        elseif(matchs('stop',4)) then
          call ProcessInput(11)
          call fftAllocate( 2 )
          exit
        elseif(endfil(dum)) then
          call ProcessInput(12)
          exit
        else
          write(out,*) ">>>Error: Unrecogonized command"
          exit
        end if
      end do
c
      contains
c
c              ****************************************************
c              *  contains: ProcessInput()                        *
c              *    assign task for each input command            *
c              ****************************************************
c
      subroutine ProcessInput(isw)
      implicit none
      integer :: isw

      select case (isw)
      case (1)
      case (2)
      case (3)
      case (4)
      case (5)
      case (6)
      case (7)
      case (8)
      case (9)
      case (10)
      case (11)
      case (12)
      case default
      end select

      return
      end subroutine
      end program
