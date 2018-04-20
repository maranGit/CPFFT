      program main
      use file_info
      use fft
      implicit none
      include 'common.main'

      logical :: readnew
      logical, external :: matchs, matchs_exact

      readnew = .true.
      call readsc()
      do while ( .true. )
        if(readnew) call readsc()
        if(matchs('project',7)) then
          cycle
        elseif(matchs('number',6)) then
          call FFT_init()
          cycle
        elseif(matchs('material',8)) then
          cycle
        elseif(matchs('sizes',5)) then
          cycle
        elseif(matchs('elements',8)) then
          cycle
        elseif(matchs('blocking',8)) then
          cycle
        elseif(matchs('strains',7)) then
          cycle
        elseif(matchs('loading',7)) then
          cycle
        elseif(matchs('nonlinear',9)) then
          cycle
        elseif(matchs('compute',7)) then
          call drive_eps_sig( 1, 0 )
          call FFT_nr3()
          cycle
        elseif(matchs('stop',4)) then
          call fftAllocate( 2 )
          exit
        else
          write(out,*) ">>>Error: Unrecogonized command"
          exit
        end if
      end do
c
      end program

