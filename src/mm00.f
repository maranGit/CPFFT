c     ****************************************************************
c     *                                                              *
c     *                 subroutine   constitutive                    *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 1/22/18                     *
c     *                                                              *
c     *         original material subroutine of FFT program          *
c     *                                                              *
c     ****************************************************************
c
      subroutine constitutive(F, phase, P, K4)
      implicit none
      real(8) :: F(9), P(9), K4(81)
c     internal variables
      integer :: trans(9)
      real(8) :: K, mu, E(9), S(9), C4(81)
      real(8) :: I2(9), I4(81), I4rt(81), I4s(81), II(81)
      real(8) :: tmp1(81), tmp2(81)
      real(8), parameter :: zero = 0.0D0, one = 1.0D0, half = 0.5D0
      real(8), parameter :: two = 2.0D0, three = 3.0D0
      integer :: phase

c     frequently used identity tensor
      trans = [ 1, 4, 7, 2, 5, 8, 3, 6, 9 ]
      I2 = zero
      I2(1) = one
      I2(5) = one
      I2(9) = one
      I4 = zero
      I4( [1, 13, 25, 29, 41, 53, 57, 69, 81] ) = one
      I4rt = zero
      I4rt( [1, 11, 21, 31, 41, 51, 61, 71, 81] ) = one
      I4s = half * ( I4 + I4rt )
      II = zero
      II( [1, 5, 9, 37, 41, 45, 73, 77, 81] ) = one

c     loop over gauss point to update stress and consistent stiffness
      P = zero
      K4 = zero
      if ( phase ) then
        K = 8.33D0;
        mu = 3.86D0;
      else
        K = 0.833D0;
        mu = 0.386D0;
      endif
      call dot22( F(trans), F, E )
      E = half * ( E - I2 )
      C4 = K * II + two * mu * ( I4s - II / three )
      call ddot42( C4, E, S )
      call dot22( F, S, P )
      call dot24( F, C4, tmp1 )
      call dot42( tmp1, F(trans), tmp2 )
      call ddot44( I4rt, tmp2, tmp1 )
      call ddot44( tmp1, I4rt, tmp2 )
      call dot24( S, I4, tmp1 )
      K4 = tmp1 + tmp2
      return
      end subroutine ! consititutive
c
c     ****************************************************************
c     *                                                              *
c     *               subroutine dot22, dot24, dot42, ddot44         *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 1/22/18                     *
c     *                                                              *
c     *                      should be inlined                       *
c     *                                                              *
c     ****************************************************************
c
      subroutine dot22(A2, B2, C2)
      implicit none
c     C2_ik = A2_ij * B2_jk
      real(8), intent(in) :: A2( 9 ), B2( 9 )
      real(8) :: C2( 9 )
      C2(1) = A2(1) * B2(1) + A2(2) * B2(4) + A2(3) * B2(7)
      C2(2) = A2(1) * B2(2) + A2(2) * B2(5) + A2(3) * B2(8)
      C2(3) = A2(1) * B2(3) + A2(2) * B2(6) + A2(3) * B2(9)
      C2(4) = A2(4) * B2(1) + A2(5) * B2(4) + A2(6) * B2(7)
      C2(5) = A2(4) * B2(2) + A2(5) * B2(5) + A2(6) * B2(8)
      C2(6) = A2(4) * B2(3) + A2(5) * B2(6) + A2(6) * B2(9)
      C2(7) = A2(7) * B2(1) + A2(8) * B2(4) + A2(9) * B2(7)
      C2(8) = A2(7) * B2(2) + A2(8) * B2(5) + A2(9) * B2(8)
      C2(9) = A2(7) * B2(3) + A2(8) * B2(6) + A2(9) * B2(9)
      return
      end subroutine

      subroutine dot24( A2, B4, C4 )
      implicit none
c     A2_ij * B4_jkmn = C4_ikmn
      real(8), intent(in) :: A2( 9 ), B4( 81 )
      real(8) :: C4( 81 )
      integer :: ii, jj, p
      C4 = 0.0D0
      do ii = 1, 3
        do jj = 1, 27
          p = ( ii - 1 ) * 27 + jj
          C4(p) = A2(ii*3-2) * B4(jj) + A2(ii*3 - 1) * B4(jj + 27)
     &          + A2(ii * 3) * B4(jj + 54)
        enddo
      enddo
      return
      end subroutine

      subroutine dot42( A4, B2, C4 )
      implicit none
c     A4_ijkl * B2_lm = C4_ijkm
      real(8), intent(in) :: A4( 81 ), B2( 9 )
      real(8) :: C4( 81 )
      integer :: ii, jj, p
      C4 = 0.0D0
      do ii = 1, 27
        do jj = 1, 3
          p = ( ii - 1 ) * 3 + jj
          C4( p ) = A4( ii * 3 - 2 ) * B2( jj )
     &            + A4(ii*3-1) * B2(jj+3) + A4(ii*3) * B2(jj+6)
        enddo
      enddo
      return
      end subroutine

      subroutine ddot42( A4, B2, C2 )
      implicit none
c     A4_ijkl * B2_lk = C2_ij
      real(8), intent(in) :: A4( 81 ), B2( 9 )
      real(8) :: C2( 9 )
      integer :: ii, jj, p
      integer :: D2( 9 )

      C2 = 0.0D0
      D2 = [ 1, 4, 7, 2, 5, 8, 3, 6, 9 ]
      do ii = 1, 9
        do jj = 1, 9
          p = D2( jj ) + ( ii - 1 ) * 9
          C2( ii ) = C2( ii ) + A4( p ) * B2( jj )
        enddo
      enddo
      return
      end subroutine

      subroutine ddot44(A4, B4, C4)
      implicit none
c     A4_ijkl * B4_lkmn = C4_ijmn
      real(8) :: A4(81), B4(81), C4(81)
      integer :: tmp1(9), tmp2(9), p, ii, jj
      tmp1 = [ 1,  2,  3, 4,  5,  6,  7,  8,  9 ]
      tmp2 = [ 0, 27, 54, 9, 36, 63, 18, 45, 72 ]
      do ii = 1, 9
        do jj = 1, 9
          p = ( ii - 1 ) * 9 + jj
          C4(p) = dot_product( A4(tmp1 + ii*9-9), B4(tmp2 + jj) )
        enddo
      enddo
      return
      end subroutine
