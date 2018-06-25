c     ****************************************************************
c     *                                                              *
c     *                   subroutine  tangent_homo                   *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *               last modified : 6/21/2018  (RM)                *
c     *                                                              *
c     ****************************************************************
c
      subroutine tangent_homo( C_homo )
      use fft, only: N3, K4, b, tolPCG
      implicit none
      include 'common.main'
c
c                global
c
      real(8) :: C_homo(*)
c
c                local
c
      integer              :: i, j
      real(8), parameter   :: zero = 0.0D0, one = 1.0D0, mone = -1.0D0
      real(8), allocatable :: C_ij(:,:)
c
c     Idea from this paper:
c       Gokuzum, Felix Selim, and Marc‐Andre Keip. "An algorithmically 
c       consistent macroscopic tangent operator for FFT‐based computa-
c       tional homogenization." International Journal for Numerical 
c       Methods in Engineering 113.4 (2018): 581-600.
c
c     Use G operator and CG solver to compute d(delta_F) / d(F_bar), 
c     where delta_F is the F perturbation and F_bar is the mean 
c     deformation gradient. Then d(P) / d(F_bar) and d(P_bar) / d(F_bar)
c     is right at hand.
c
      allocate( C_ij(N3, nstrs) )
      do i = 1, nstrs
c
c                form b = Ghat4 : K4
c
        call dcopy(N3, K4(1, i     ), 1, C_ij(1, 1), 1)
        call dcopy(N3, K4(1, i +  9), 1, C_ij(1, 2), 1)
        call dcopy(N3, K4(1, i + 18), 1, C_ij(1, 3), 1)
        call dcopy(N3, K4(1, i + 27), 1, C_ij(1, 4), 1)
        call dcopy(N3, K4(1, i + 36), 1, C_ij(1, 5), 1)
        call dcopy(N3, K4(1, i + 45), 1, C_ij(1, 6), 1)
        call dcopy(N3, K4(1, i + 54), 1, C_ij(1, 7), 1)
        call dcopy(N3, K4(1, i + 63), 1, C_ij(1, 8), 1)
        call dcopy(N3, K4(1, i + 72), 1, C_ij(1, 9), 1)
        call G_K_dF(C_ij(1,1), b(1,1), .false.)
        call dscal(N3*nstrs, mone, b(1,1), 1)

c          use CG solver for d(delta_F) / d(F_bar)

        call fftPcg(b, C_ij, tolPCG, out)

c     d(F)/d(F_bar) = delta(i,k) * delta(j,l) + d(delta_F) / d(F_bar)

        C_ij(1:N3,i) = C_ij(1:N3,i) + one

c        C_homo = mean{ [ d(P)/d(F) ] : [ d(F)/d(F_bar) ] }

        do j = 1, nstrs
          call tangent_average(C_ij(1,j), C_homo(9*j+i-9), N3)
        end do

      end do
c
      deallocate( C_ij )
c                
      return
      end subroutine

      subroutine tangent_average(a, b, n)
      implicit none
      real(8), intent(in)  :: a(*)
      integer, intent(in)  :: n
      real(8), intent(out) :: b
      integer              :: i
      real(8), parameter   :: zero = 0.0D0

      b = zero
      do i = 1, n
        b = b + a(i)
      enddo
      b = b / dble(n)

      return
      end subroutine
