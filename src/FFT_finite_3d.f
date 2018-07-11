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
        elseif(endfil(dum)) then
          call ProcessInput(13)
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
        call setstarttime
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
            if(.not.numd(l_x)) 
     &        call errmsg(23,dum,'x_direction',dumr,dumd)
          elseif(matchs('y_direction',4)) then
            if(.not.numd(l_y)) 
     &        call errmsg(23,dum,'y_direction',dumr,dumd)
          elseif(matchs('z_direction',4)) then
            if(.not.numd(l_z)) 
     &        call errmsg(23,dum,'z_direction',dumr,dumd)
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
c       call fftAllocate( 3 )
c                       Set up the crystal element properties as well
        call read_crystal_data
        call read_simple_angles
        call avg_cry_elast_props
        readnew = .true.
      case (7)
        call inlodcase()
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
        call compute_checks
c
c           before analysis of loca (time) step 1, perform
c           system-wide array initialization.
c
        call fftAllocate(3)
        call drive_eps_sig( 1, 0 )
        call FFT_nr3()
        readnew = .true.
      case (11)
        call oudrive()
        readnew = .true.
      case (12)
        call CPFFT_normal_stop
      case (13)
        call CPFFT_normal_stop
      case (14)
        call incrystal(out)
        readnew = .true.
      case default
        call errmsg(4,dum,dums,dumr,dumd)
      end select

      return
      end subroutine
      end program
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine CPFFT_normal_stop                 *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified : 7/11/2018 rhd              *
c     *                                                              *
c     *     execute a normal shutdown with messages, etc.            *
c     *                                                              *
c     ****************************************************************
c
c

      subroutine CPFFT_normal_stop
c
      use file_info
      implicit none
      include 'common.main'
c
      integer :: i
      real :: t1, dumr, warptime, pardiso_time
      real, external :: wcputime
      character(len=8) :: stcnam, dums, sdate_*24
      character(len=80) :: name, stflnm, rtflnm
      logical :: hilcmd, sbflg1, sbflg2
      logical :: endcrd, label, matchs, debug1, debug2, debug, endfil,
     &           string, matchs_exact
      logical, parameter :: output_1_thread_cpu_times = .false.
c
c
c                       cleanup some allocs first
c
      call cleanup_crystal
      call fftAllocate( 2 )
c
      write(out,*)
      write(out,*)
      t1 = wcputime ( 1 )
      write(out,'(">> total job wall time (secs): ", f12.2)') t1
c
c                      close input and output files
c
      close (out)
      if( filcnt .gt. 1 ) then
         do i = filcnt, 2, -1
            close(i)
         end do
      end if
c
c         MPI:
c            tell workers we are ending now and then stop.
c         threads - Fortran stop
c
      call die_gracefully
c
      end
