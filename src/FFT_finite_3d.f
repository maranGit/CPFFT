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
      do while ( .true. )
        if(readnew) call readsc()
        if(matchs('project',7)) then
          call ProcessInput(1)
        elseif(matchs('number',6)) then
          call ProcessInput(2)
        elseif(matchs('material',8)) then
          call ProcessInput(3)
        elseif(matchs('crystal',7)) then
          call ProcessInput(14)
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
        elseif(matchs('output',6)) then
          call ProcessInput(11)
        elseif(matchs('stop',4)) then
          call ProcessInput(12)
          exit
        elseif(endfil(dum)) then
          call ProcessInput(13)
          exit
        else
          call errmsg(4,dum,dums,dumr,dumd)
          readnew = .true.
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

c                        global
      integer :: isw
c                        local
      logical, external :: label, endcrd, integr, numd
      integer :: nc
      character(len=80) :: name

      select case (isw)
      case (1)
        if(matchs('name',4)) call splunj
        if(label(dummy)) then
          stname = ' '
          name= ' '
          call entits(name,nc)
          if(nc.gt.8) nc=8
          stname(1:nc)= name(1:nc)
        else
          call errmsg(1,dum,dums,dumr,dumd)
        endif
        readnew = .true.
      case (2)
        do while(.not.endcrd(dum))
          if(matchs('of',2)) call splunj
          if(matchs('grid',4)) then
            if(.not.integr(N)) call errmsg(2,dum,dums,dumr,dumd)
          else
            call errmsg(4,dum,dums,dumr,dumd)
            call scan()
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
            call scan()
          endif
        end do
        readnew = .true.
      case (5)
        call inelem()
        readnew = .false.
      case (6)
        if(matchs('automatic',4)) call splunj
        call inelbk( matList )
        call fftAllocate( 3 )
        readnew = .true.
      case (7)
        do while ( .true. )
          call readsc()
          if     ( matchs('F_xx',4) ) then
            if ( .not. numd( F_total(1) ) ) 
     &         call errmsg(13,dum,dums,dumr,dumd)
          elseif ( matchs('F_xy',4) ) then
            if ( .not. numd( F_total(2) ) ) 
     &         call errmsg(13,dum,dums,dumr,dumd)
          elseif ( matchs('F_xz',4) ) then
            if ( .not. numd( F_total(3) ) ) 
     &         call errmsg(13,dum,dums,dumr,dumd)
          elseif ( matchs('F_yx',4) ) then
            if ( .not. numd( F_total(4) ) ) 
     &         call errmsg(13,dum,dums,dumr,dumd)
          elseif ( matchs('F_yy',4) ) then
            if ( .not. numd( F_total(5) ) ) 
     &         call errmsg(13,dum,dums,dumr,dumd)
          elseif ( matchs('F_yz',4) ) then
            if ( .not. numd( F_total(6) ) ) 
     &         call errmsg(13,dum,dums,dumr,dumd)
          elseif ( matchs('F_zx',4) ) then
            if ( .not. numd( F_total(7) ) ) 
     &         call errmsg(13,dum,dums,dumr,dumd)
          elseif ( matchs('F_zy',4) ) then
            if ( .not. numd( F_total(8) ) ) 
     &         call errmsg(13,dum,dums,dumr,dumd)
          elseif ( matchs('F_zz',4) ) then
            if ( .not. numd( F_total(9) ) ) 
     &         call errmsg(13,dum,dums,dumr,dumd)
          else
            exit
          end if
        end do
        readnew = .false.
      case (8)
        call inlod()
        readnew = .false.
      case (9)
        if(matchs('analysis',8)) call splunj
        if(matchs('parameters',5)) call splunj
        call indypm()
        readnew = .false.
      case (10)
        call drive_eps_sig( 1, 0 )
        call compute_checks
        call FFT_nr3()
        readnew = .true.
      case (11)
        call oudrive()
        readnew = .true.
      case (12)
        call fftAllocate( 2 )
      case (13)
        call fftAllocate( 2 )
      case (14)
        call incrystal(out)
        readnew = .true.
      case default
        call errmsg(4,dum,dums,dumr,dumd)
      end select

      return
      end subroutine
      end program
