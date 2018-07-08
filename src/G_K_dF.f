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

c     input variables
      real(8), intent(in) :: F(N3, ndim2)
      real(8), intent(out) :: GKF(N3, ndim2)
      logical flgK
c
c               local
      logical :: local_debug
      integer :: ii

      local_debug = .false.

c     multiply by K4
      if (flgK) then
        call ddot42n(K4(1,1), F(1,1), real1(1,1), N3)
      else
        call dcopy(N3*ndim2, F(1,1), 1, real1(1,1), 1)
      endif

c     fft
      do ii = 1, ndim2
        call fftfem3d(real1(1,ii), real2(1,ii), real3(1,ii))
      enddo
c
c                     multiply by G_hat matrix
c             original complex array is (real2, real3)
c             after convolution, it is (real1, real2)
c
      call ddot42n(Ghat4(1,1), real2(1,1), real1(1,1), N3)
      call ddot42n(Ghat4(1,1), real3(1,1), real2(1,1), N3)
c
c     inverse fft
      do ii = 1, ndim2
c       call ifftfem3d(tmpCplx(1,ii), GKF(1,ii))
        call ifftfem3d(GKF(1,ii), real1(1,ii), real2(1,ii))
      enddo
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
      call vdmul( N3, in_re, coeffs1, out_re )
      call vdmul( N3, in_re, coeffs2, out_im )
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
      call vdmul( N3, out_re, coeffs1, in_re )
      call vdmul( N3, out_im, coeffs2, out_re )
      call daxpy( N3, one, out_re, 1, in_re, 1 )
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
c     *             subroutine   trans2, ddot42n, ddot42n_cpmlx      *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 1/22/18                     *
c     *                                                              *
c     *                       should be inlined                      *
c     *                                                              *
c     ****************************************************************
c
      subroutine ddot42n(A4, B2, C2, n)
c
c                  A4_ijkl * B2_kl = C2_ij, n rows
c
      implicit none
c
c                    global
c
      integer, intent(in)  :: n
      real(8), intent(in)  :: A4(n,*), B2(n,*)
      real(8), intent(out) :: C2(n,*)
c
c                    local
c
      integer :: i, j, k, l
      real(8), parameter :: zero = 0.0D0
c
      C2(1:n,1:9) = zero

      do j = 1, 9
        l = ( j - 1 ) * 9

        do k = 1, 9

          do i = 1, n

            C2(i,j) = C2(i,j) + A4(i,l+k) * B2(i,k)

          end do

        end do

      end do

      return
      end subroutine
c
c     subroutine ddot42n_cmplx(A4, B2, n)
c                 A4_ijkl * B2_kl = C2_ij, n rows
c                 A4 is real, B2 is complex, C2 is complex
c     implicit none
c
c                    global
c
c     integer, intent(in)     :: n
c     real(8), intent(in)     :: A4( n, * )
c     complex(8)              :: B2( n, * )
c
c                    local
c
c     complex(8), allocatable :: C2( :, : )
c     integer                 :: i, j, k, l
c     real(8), parameter      :: zero = 0.0D0
c
c     allocate( C2( n, 9 ) )
c     C2 = (zero, zero)
c     D2 = [ 1, 4, 7, 2, 5, 8, 3, 6, 9 ]
c     do ii = 1, 9
c       do jj = 1, 9
c         p = D2(jj) + (ii - 1) * 9
c         C2(:, ii) = C2(:, ii) + dcmplx(A4(:, p),0.0D0) * B2(:, jj)
c       enddo
c     enddo
c
c     do j = 1, 9
c       l = ( j - 1 ) * 9
c
c       do k = 1, 9
c
c         do i = 1, n
c
c           C2(i,j) = C2(i,j) + dcmplx(A4(i,l+k),0.0D0) * B2(i,k)
c
c         end do
c
c       end do
c
c     end do
c
c     do i = 1, 9
c       B2(1:n,i) = C2(1:n,i)
c     end do
c
c     deallocate( C2 )
c
c     return
c     end subroutine
