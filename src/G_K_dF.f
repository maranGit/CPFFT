c     ****************************************************************
c     *                                                              *
c     *                     subroutine   G_K_dF                      *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 1/22/18                     *
c     *                                                              *
c     ****************************************************************
c
      subroutine G_K_dF(F, GKF, flgK)
      use fft, only: N3, ndim2, Ghat4, K4, real1, real2, real3
      implicit none
c
c                        global
c
      real(8), intent(in)  :: F(N3, *)
      real(8), intent(out) :: GKF(N3, *)
      logical, intent(in)  :: flgK
c
c                        local
c
      logical :: local_debug
      integer :: ii, now_thread
      integer, external :: omp_get_thread_num
c
      data local_debug  / .false. /
c
c                    multiply by K4 if required
c
      if (flgK) then
        call ddot42n(K4(1,1), F(1,1), real1(1,1), GKF(1,1), N3)
      else
        call dcopy(N3*ndim2, F(1,1), 1, real1(1,1), 1)
      endif
c
c                               FFT
c
c     In the original algorithm, this step should be
c       ifftshift( fftn( fftshift( x ) ) ).
c     Here, multiplication with Fourier coefficient is employed
c       as pre-processing equivalent to fftshift and ifftshift.
c     There should be two multiplication before FFT and after FFT.
c     The second one is omitted because it cancled with the 
c       multiplication before IFFT.
c     Therefore, the intermediate result should be different from
c       my corresponding MATLAB code.
c
      call omp_set_dynamic( .false. )
c$OMP PARALLEL DO ORDERED 
c$OMP&         PRIVATE( ii, now_thread )
c$OMP&         SHARED( real1, real2, real3, ndim2 )
      do ii = 1, ndim2
        now_thread = omp_get_thread_num() + 1
        call fftfem3d(real1(1,ii), real2(1,ii), real3(1,ii))
      enddo
c$OMP END PARALLEL DO
c
c                     multiply by G_hat matrix
c             original complex array is (real2, real3)
c             after convolution,  it is (real1, real2)
c             GKF serves as a temporary array
c
      call ddot42n(Ghat4(1,1), real2(1,1), real1(1,1), GKF(1,1), N3)
      call ddot42n(Ghat4(1,1), real3(1,1), real2(1,1), GKF(1,1), N3)
c
c                              IFFT
c
c     In the original algorithm, this step should be
c       ifftshift( ifftn( fftshift( x ) ) ).
c     Here, multiplication with Fourier coefficient is employed
c       as pre-processing equivalent to fftshift and ifftshift.
c     There should be two multiplication before IFFT and after IFFT.
c     The first one is omitted because it cancled with the 
c       multiplication after FFT.
c
c$OMP PARALLEL DO ORDERED 
c$OMP&         PRIVATE( ii, now_thread )
c$OMP&         SHARED( real1, real2, GKF, ndim2 )
      do ii = 1, ndim2
        now_thread = omp_get_thread_num() + 1
        call ifftfem3d(GKF(1,ii), real1(1,ii), real2(1,ii))
      enddo
c$OMP END PARALLEL DO
c
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine       fftfem3d                    *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 1/22/18                     *
c     *                                                              *
c     *                transform real1 to cplx1                      *
c     *                                                              *
c     ****************************************************************
c
      subroutine fftfem3d(in_re, out_re, out_im)
      use fft, only: N3, N, Nhalf, ndim1, dims, coeffs1, coeffs2
      use MKL_DFTI
      implicit none
c
c                        global
c
      real(8) :: in_re(*), out_re(*), out_im(*)

c
c                        local
c
      type(DFTI_DESCRIPTOR), pointer :: MyHandle
      integer :: istrides(4), ostrides(4), fftStatus, tmp
      real(8), parameter :: zero = 0.0D0, one = 1.0D0
      logical :: debug

c     frequently used array
      istrides = ( / 0, 1, N, N * N / )
      ostrides = ( / 0, 1, N, N * N / )
      debug = .false.
c
c                          perform fftshift
c
      call vdmul( N3, in_re(1), coeffs1(1), out_re(1) )
      call vdmul( N3, in_re(1), coeffs2(1), out_im(1) )
c
c                  create and personalize descriptor
c
      fftStatus = DftiCreateDescriptor
     &            (MyHandle, DFTI_DOUBLE, DFTI_COMPLEX, ndim1, dims)
      fftStatus = DftiSetValue
     &            (MyHandle, DFTI_COMPLEX_STORAGE, DFTI_REAL_REAL)
      fftStatus = DftiSetValue
     &            (MyHandle, DFTI_PLACEMENT, DFTI_INPLACE)
      fftStatus = DftiSetValue
     &            (MyHandle, DFTI_INPUT_STRIDES, istrides)
      fftStatus = DftiSetValue
     &            (MyHandle, DFTI_OUTPUT_STRIDES, ostrides)
      fftStatus = DftiSetValue(MyHandle, DFTI_FORWARD_SCALE, one)
      fftStatus = DftiSetValue(MyHandle, DFTI_BACKWARD_SCALE,
     &                         one/dble(N3))
c
c                         commit descriptor
c
      fftStatus = DftiCommitDescriptor(MyHandle)
c
c                         compute a forward FFT
c
      fftStatus = DftiComputeForward(MyHandle, out_re, out_im)
c
c                          free descriptor
c
      fftStatus = DftiFreeDescriptor(MyHandle)
c
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine       ifftfem3d                   *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 1/22/18                     *
c     *                                                              *
c     *                  transform cplx1 to real1                    *
c     *                                                              *
c     ****************************************************************
c
      subroutine ifftfem3d(in_re, out_re, out_im)
      use fft, only: N3, N, Nhalf, ndim1, dims, coeffs1, coeffs2
      use MKL_DFTI
      implicit none
c
c                        global
c
      real(8) :: in_re(*), out_re(*), out_im(*)
c
c                        local
c
      type(DFTI_DESCRIPTOR), pointer :: MyHandle
      integer :: istrides(4), ostrides(4), fftStatus, tmp
      real(8), parameter :: zero = 0.0D0, one = 1.0D0
      logical :: debug

c     frequently used array
      istrides = ( / 0, 1, N, N * N / )
      ostrides = ( / 0, 1, N, N * N / )
      debug = .false.
c
c                  create and personalize descriptor
c
      fftStatus = DftiCreateDescriptor
     &            (MyHandle, DFTI_DOUBLE, DFTI_COMPLEX, ndim1, dims)
      fftStatus = DftiSetValue
     &            (MyHandle, DFTI_COMPLEX_STORAGE, DFTI_REAL_REAL)
      fftStatus = DftiSetValue
     &            (MyHandle, DFTI_PLACEMENT, DFTI_INPLACE)
      fftStatus = DftiSetValue
     &            (MyHandle, DFTI_INPUT_STRIDES, istrides)
      fftStatus = DftiSetValue
     &            (MyHandle, DFTI_OUTPUT_STRIDES, ostrides)
      fftStatus = DftiSetValue(MyHandle, DFTI_FORWARD_SCALE, one)
      fftStatus = DftiSetValue(MyHandle, DFTI_BACKWARD_SCALE,
     &                         one/dble(N3))
c
c                         commit descriptor
c
      fftStatus = DftiCommitDescriptor(MyHandle)
c
c                         compute a forward FFT
c
      fftStatus = DftiComputeBackward(MyHandle, out_re, out_im)
c
c                          perform ifftshift
c
      call vdmul( N3, out_re(1), coeffs1(1), in_re(1)  )
      call vdmul( N3, out_im(1), coeffs2(1), out_re(1) )
      call daxpy( N3, one      , out_re(1) , 1, in_re(1), 1 )
c
c                          free descriptor
c
      fftStatus = DftiFreeDescriptor(MyHandle)
c
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                    subroutine   ddot42n                      *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                    last modified: 8/7/18                     *
c     *                                                              *
c     *                       should be inlined                      *
c     *                                                              *
c     ****************************************************************
c
      subroutine ddot42n(A4, B2, C2, tmp, n)
c
c                  A4_ijkl * B2_kl = C2_ij, n rows
c
      implicit none
c
c                    global
c
      integer, intent(in)  :: n
      real(8), intent(in)  :: A4(n,*), B2(n,*)
      real(8), intent(out) :: C2(n,*), tmp(n,*)
c
c                    local
c
      integer :: i, j
      real(8), parameter :: one = 1.0D0
c
      do i = 1, 9
        j = i * 9 - 8
        call vdmul( n*9, A4(1,j), B2(1,1), tmp(1,1) )
        call daxpy( n*4, one, tmp(1,6), 1, tmp(1,2), 1 )
        call daxpy( n*2, one, tmp(1,4), 1, tmp(1,2), 1 )
        call daxpy( n  , one, tmp(1,3), 1, tmp(1,2), 1 )
        call vdadd( n, tmp(1,1), tmp(1,2), C2(1,i) )
      end do
c
      return
      end subroutine
