      program main
      implicit none
      include 'common.main'
      
      out = 11
      open(out, file = "RanOut.out")

c               initialize all global arrays
      call FFT_init()
      
c               extract initial stiffness
      call drive_eps_sig( 1, 0 )

c              global Newton-Raphson loop
      call FFT_nr3()

c              deallocate all arrays
      call fftAllocate( 2 )
      
      close(out)
c
      end program

