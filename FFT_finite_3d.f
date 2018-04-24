      program main
      use file_info
      implicit none
      include 'common.main'

      integer :: dum, dummy
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
        elseif(matchs('number',6)) then
          call ProcessInput(2)
        elseif(matchs('material',8)) then
          call ProcessInput(3)
        elseif(matchs('sizes',5)) then
          call ProcessInput(4)
        elseif(matchs('elements',8)) then
          call ProcessInput(5)
        elseif(matchs('blocking',8)) then
          call ProcessInput(6)
        elseif(matchs('strains',7)) then
          call ProcessInput(7)
        elseif(matchs('loading',7)) then
          call ProcessInput(8)
        elseif(matchs('nonlinear',9)) then
          call ProcessInput(9)
        elseif(matchs('compute',7)) then
          call ProcessInput(10)
        elseif(matchs('stop',4)) then
          call ProcessInput(11)
          exit
        elseif(endfil(dum)) then
          call ProcessInput(12)
          exit
        else
          call errmsg(4,dum,dums,dumr,dumd)
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
      use fft
      implicit none
      integer :: isw
      logical, external :: label, endcrd, integr, numd

      select case (isw)
      case (1)
        if(matchs('name',4)) call splunj
        if(label(dummy)) then
c         name = 1
        else
          call errmsg(1,dum,dums,dumr,dumd)
        endif
        readnew = .true.
      case (2)
        do while(.not.endcrd(dum))
          if(matchs('of',2)) call splunj
          if(matchs('grid',4)) then
            if(.not.integr(N)) call errmsg(2,dum,dums,dumr,dumd)
          elseif(matchs('materials',4)) then
            if(.not.integr(nummat)) call errmsg(3,dum,dums,dumr,dumd)
          else
            call errmsg(4,dum,dums,dumr,dumd)
          endif
        enddo
        call FFT_init()
        readnew = .true.
      case (3)
        call inmat()
        readnew = .true.
      case (4)
        if(matchs('of',2)) call splunj
        do while (.not.endcrd(dum))
          if(matchs('x_direction',4)) then
            if(.not.numd(l_x)) call errmsg(4,dum,dums,dumr,dumd)
          elseif(matchs('y_direction',4)) then
            if(.not.numd(l_y)) call errmsg(4,dum,dums,dumr,dumd)
          elseif(matchs('z_direction',4)) then
            if(.not.numd(l_z)) call errmsg(4,dum,dums,dumr,dumd)
          else
            call errmsg(4,dum,dums,dumr,dumd)
          endif
        end do
        readnew = .true.
      case (5)
        call inelem()
      case (6)
      case (7)
      case (8)
      case (9)
      case (10)
        call drive_eps_sig( 1, 0 )
        call FFT_nr3()
      case (11)
        call fftAllocate( 2 )
      case (12)
      case default
      end select

      return
      end subroutine
      end program
