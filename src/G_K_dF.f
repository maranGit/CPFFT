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
      use fft, only: N3, ndim2, Ghat4, K4, tmpReal, tmpCplx
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
      tmpReal = F
      if (flgK) then
        call trans2(tmpReal, N3)
        call ddot42n(K4, tmpReal, N3)
        call trans2(tmpReal, N3)
      endif

c     fft
      do ii = 1, ndim2
        call fftfem3d(tmpReal(:,ii), tmpCplx(:,ii))
      enddo

      if(local_debug) then
        write(*,*) "Now checking tmpCplx"
        write(*,*) tmpCplx(1,1)
      endif

c     multiply by G_hat matrix
      call ddot42n_cmplx(Ghat4, tmpCplx, N3)

c     inverse fft
      do ii = 1, ndim2
        call ifftfem3d(tmpCplx(:,ii), tmpReal(:,ii))
      enddo

c     return
      GKF = tmpReal

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
      subroutine fftfem3d(real3, cplx3)
      use fft, only: real1, cplx3half, cplx1half, N3, N, Nhalf,
     &               ndim1, dims
      use MKL_DFTI
      implicit none

c     input variables
      real(8) :: real3(N, N, N)
      complex(8) :: cplx3(N, N, N) ! return as a N3x1 vector

c     internal variables
      integer :: istrides(4), ostrides(4)
      type(DFTI_DESCRIPTOR), pointer :: MyHandle
      integer :: fftStatus, tmp
      logical debug
      integer, dimension(:), allocatable :: ffts, iffts
      complex(8), parameter :: zero = (0.0D0, 0.0D0)
      real :: fftrange, fftmax, fftmin

c     allocate local variables
      allocate( ffts(N) )
      allocate( iffts(N) )

c     frequently used array
      ffts = [ (Nhalf + 1) : N , 1 : Nhalf ]
      iffts = [ Nhalf : N, 1 : (Nhalf - 1) ]
      istrides = [0, 1, N    , N*N    ]
      ostrides = [0, 1, Nhalf, Nhalf*N]
      debug = .false.

c     ifftshift
      real3 = real3(iffts, iffts, iffts)
      real1 = reshape(real3, (/N3/))

c     create and personalize descriptor
      fftStatus = DftiCreateDescriptor(MyHandle, DFTI_DOUBLE,
     &                                 DFTI_REAL, ndim1, dims)
      fftStatus = DftiSetValue(MyHandle, DFTI_CONJUGATE_EVEN_STORAGE,
     &                         DFTI_COMPLEX_COMPLEX)
      fftStatus = DftiSetValue(MyHandle, DFTI_PACKED_FORMAT,
     &                         DFTI_CCE_FORMAT)
      fftStatus = DftiSetValue(MyHandle, DFTI_PLACEMENT,
     &                         DFTI_NOT_INPLACE)
      fftStatus = DftiSetValue(MyHandle, DFTI_INPUT_STRIDES, istrides)
      fftStatus = DftiSetValue(MyHandle, DFTI_OUTPUT_STRIDES, ostrides)

c     commit descriptor
      fftStatus = DftiCommitDescriptor(MyHandle)

c     compute a forward FFT
      fftStatus = DftiComputeForward(MyHandle, real1, cplx1half)
      cplx3half = reshape(cplx1half, (/Nhalf, N, N/))

c     free descriptor
      fftStatus = DftiFreeDescriptor(MyHandle)

c     construct the whole matrix and perform fftshift
      cplx3( 1:Nhalf  , :  , :   ) = cplx3half
      cplx3( Nhalf+1:N, 1  , 1   ) = 
     &                 dconjg( cplx3half( Nhalf:2:-1, 1     , 1      ) )
      cplx3( Nhalf+1:N, 1  , 2:N ) = 
     &                 dconjg( cplx3half( Nhalf:2:-1, 1     , N:2:-1 ) )
      cplx3( Nhalf+1:N, 2:N, 1   ) = 
     &                 dconjg( cplx3half( Nhalf:2:-1, N:2:-1, 1      ) )
      cplx3( Nhalf+1:N, 2:N, 2:N ) = 
     &                 dconjg( cplx3half( Nhalf:2:-1, N:2:-1, N:2:-1 ) )
      cplx3                        = cplx3(ffts, ffts, ffts)

c     deallocate local variables
      deallocate( ffts, iffts )

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
      subroutine ifftfem3d(cplx3, real3)
      use fft, only: real1, cplx3half, cplx1half, N3, N, Nhalf,
     &               ndim1, dims
      use MKL_DFTI
      implicit none

c     input variables
      complex(8) :: cplx3(N,N,N)
      real(8) :: real3(N,N,N)

c     internal variables
      integer :: istrides(4), ostrides(4)
      type(DFTI_DESCRIPTOR), pointer :: MyHandle
      integer :: fftStatus, tmp
      logical debug
      real(8), parameter :: one = 1.0D0
      integer, dimension(:), allocatable :: ffts, iffts

c     allocate local variables
      allocate( ffts(N) )
      allocate( iffts(N) )

c     initial parameters for the descriptor
      ffts = [ (Nhalf + 1) : N , 1 : Nhalf ]
      iffts = [ Nhalf : N, 1 : (Nhalf - 1) ]
      istrides = [0, 1, Nhalf, Nhalf*N]
      ostrides = [0, 1, N    , N*N    ]
      debug = .false.

c     ifftshift
      cplx3 = cplx3(iffts, iffts, iffts)

c     prepare for even-storate input
      cplx3half = cplx3(1:Nhalf, :, :)
      cplx1half = reshape(cplx3half, (/Nhalf*N*N/))

c     create and personalize descriptor
      fftStatus = DftiCreateDescriptor(MyHandle, DFTI_DOUBLE,
     &                                 DFTI_REAL, ndim1, dims)
      fftStatus = DftiSetValue(MyHandle, DFTI_CONJUGATE_EVEN_STORAGE,
     &                         DFTI_COMPLEX_COMPLEX)
      fftStatus = DftiSetValue(MyHandle, DFTI_PACKED_FORMAT,
     &                         DFTI_CCE_FORMAT)
      fftStatus = DftiSetValue(MyHandle, DFTI_PLACEMENT,
     &                         DFTI_NOT_INPLACE)
      fftStatus = DftiSetValue(MyHandle, DFTI_INPUT_STRIDES, istrides)
      fftStatus = DftiSetValue(MyHandle, DFTI_OUTPUT_STRIDES, ostrides)
      fftStatus = DftiSetValue(MyHandle, DFTI_FORWARD_SCALE, one)
      fftStatus = DftiSetValue(MyHandle, DFTI_BACKWARD_SCALE,
     &                         one/dble(N3))

c     commit descriptor
      fftStatus = DftiCommitDescriptor(MyHandle)

c     compute a forward FFT
      fftStatus = DftiComputeBackward(MyHandle, cplx1half, real1)

c     free descriptor
      fftStatus = DftiFreeDescriptor(MyHandle)

      real3 = reshape( real1, (/N, N, N/) )
      real3 = real3(ffts, ffts, ffts)

c     deallocate local variables
      deallocate( ffts, iffts )

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
      subroutine trans2(A2, n)
      implicit none
c     B2_ji = A2_ij, n rows
      integer, intent(in) :: n
      real(8) :: A2( n, 9 )
      A2 = A2( :, [ 1, 4, 7, 2, 5, 8, 3, 6, 9 ] )
      return
      end subroutine
c
      subroutine ddot42n(A4, B2, n)
      implicit none
c                  A4_ijkl * B2_lk = C2_ij, n rows
      integer, intent(in) :: n
      real(8), intent(in) :: A4( n, 81 )
      real(8) :: B2( n, 9 ), C2( n, 9 )
      integer :: ii, jj, p
      integer :: D2( 9 )

      C2 = 0.0D0
      D2 = [ 1, 4, 7, 2, 5, 8, 3, 6, 9 ]
      do ii = 1, 9
        do jj = 1, 9
          p = D2(jj) + (ii - 1) * 9
          C2(:, ii) = C2(:, ii) + A4(:, p) * B2(:, jj)
        enddo
      enddo

      B2 = C2

      return
      end subroutine
c
      subroutine ddot42n_cmplx(A4, B2, n)
c                 A4_ijkl * B2_lk = C2_ij, n rows
c                 A4 is real, B2 is complex, C2 is complex
      implicit none
      
      integer, intent(in) :: n
      real(8), intent(in) :: A4( n, 81 )
      complex(8) :: B2( n, 9 ), C2( n, 9 )
      integer :: ii, jj, p
      integer :: D2( 9 )

      C2 = (0.0D0, 0.0D0)
      D2 = [ 1, 4, 7, 2, 5, 8, 3, 6, 9 ]
      do ii = 1, 9
        do jj = 1, 9
          p = D2(jj) + (ii - 1) * 9
          C2(:, ii) = C2(:, ii) + dcmplx(A4(:, p),0.0D0) * B2(:, jj)
        enddo
      enddo

      B2 = C2

      return
      end subroutine
