c     ****************************************************************
c     *                                                              *
c     *                      subroutine cep2A                        *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                 last modified : 6/12/2018 RM                 *
c     *                                                              *
c     *       pull back from dt/dd to dP/dF, exact analytical        *
c     *       formulation, verified against finite difference        *
c     *                                                              *
c     ****************************************************************
c
      subroutine cep2A( local_work, cep, rnh, detF, detFh,
     &                  fnhinv, fn1inv, dPdF )
      implicit none
      include 'param_def'
      include 'include_sig_up'
c
c                           global
      real(8), intent(in)  :: cep(mxvl,nstr,*)
      real(8), intent(in)  :: rnh(mxvl,3,*)
      real(8), intent(in)  :: detF(*), detFh(*)
      real(8), intent(in)  :: fnhinv(mxvl,3,*), fn1inv(mxvl,3,*)
      real(8), intent(out) :: dPdF(mxvl,*)
c
c                           local
      integer :: span, felem, ii
      real(8) :: Fn33(3,3), t33(3,3), C66(6,6), Rh33(3,3), detFnh
      real(8) :: fnhinv33(3,3), R33(3,3), Fn133(3,3), finv33(3,3)
      real(8) :: detFn1, dPdF81(81)
      real(8), parameter :: zero = 0.0D0, one = 1.0D0
c
c                        preparation
      span  = local_work%span
      felem = local_work%felem
c
      do ii = 1, span
        Fn33 = local_work%fn(ii,1:3,1:3)
        t33(1,1) = local_work%urcs_blk_n1(ii,1,1)
        t33(2,1) = local_work%urcs_blk_n1(ii,4,1)
        t33(3,1) = local_work%urcs_blk_n1(ii,6,1)
        t33(1,2) = local_work%urcs_blk_n1(ii,4,1)
        t33(2,2) = local_work%urcs_blk_n1(ii,2,1)
        t33(3,2) = local_work%urcs_blk_n1(ii,5,1)
        t33(1,3) = local_work%urcs_blk_n1(ii,6,1)
        t33(2,3) = local_work%urcs_blk_n1(ii,5,1)
        t33(3,3) = local_work%urcs_blk_n1(ii,3,1)
        C66 = cep(ii,1:6,1:6)
        Rh33 = rnh(ii,1:3,1:3)
        detFnh = detFh(ii)
        fnhinv33 = fnhinv(ii,1:3,1:3)
        R33(1,1) = local_work%rot_blk_n1(ii,1,1)
        R33(2,1) = local_work%rot_blk_n1(ii,2,1)
        R33(3,1) = local_work%rot_blk_n1(ii,3,1)
        R33(1,2) = local_work%rot_blk_n1(ii,4,1)
        R33(2,2) = local_work%rot_blk_n1(ii,5,1)
        R33(3,2) = local_work%rot_blk_n1(ii,6,1)
        R33(1,3) = local_work%rot_blk_n1(ii,7,1)
        R33(2,3) = local_work%rot_blk_n1(ii,8,1)
        R33(3,3) = local_work%rot_blk_n1(ii,9,1)
        Fn133 = local_work%fn1(ii,1:3,1:3)
        finv33 = fn1inv(ii,1:3,1:3)
        detFn1 = detF(ii)
        dPdF81 = zero
        call cep2A_a(Fn33, t33, C66, Rh33, detFnh, fnhinv33, R33, Fn133,
     &               finv33, detFn1, dPdF81)
        dPdF(ii,1:81) = dPdF81(1:81)
      end do
c
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine cep2A_a                      *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                 last modified : 6/12/2018 RM                 *
c     *                                                              *
c     *         Pull back from dt/dd to dP/dF on ONE ELEMENT         *
c     *                                                              *
c     ****************************************************************
c
      subroutine cep2A_a(Fn, t, C, Rh, detFnh, fnhinv, R, Fn1, 
     &                   finv, detFn1, dPdF)
      implicit none
c
c                           global
c
c                     reference  configuration
      real(8), intent(in)  :: Fn(3,3)
c
c                   intermediate configuration
      real(8), intent(in)  :: t(3,3), C(6,6), Rh(3,3), fnhinv(3,3)
      real(8), intent(in)  :: detFnh
c
c                      current   configuration
      real(8), intent(in)  :: R(3,3), Fn1(3,3), finv(3,3), detFn1
c
c                           output
      real(8), intent(out) :: dPdF(81)
c
c                           local
c
      integer :: indices(4,81), index_voigt(9)
      integer :: m,n,i,j,k,l,a,p,q,tmp
      real(8) :: temp(3,3)
c
c              polar decomposition in both configuration
      real(8) :: transR(3,3) , U(3,3) , traceU , Y(3,3) , RYR(3,3)
      real(8) :: detY_1 , RY(3,3)
      real(8) :: transRh(3,3), Uh(3,3), traceUh, Yh(3,3), RYRh(3,3)
      real(8) :: detYh_1, RYh(3,3)
c
c                        tangent stiffness
      real(8) :: cep99(9,9),dJdF(3,3),dRdF(81),dRhdF(81),dLdF(81)
      real(8) :: dtdF(81),dtdF_1(81),dtdF_2(81),dtdF_3(81)
      real(8) :: dPdF0(81),dPdF1(81),dPdF2(81),dPdF3(81),dPdF4(81)
c
c                        stress
      real(8) :: sigma(3,3),t_Uinv(3,3),R_t(3,3),RtRF(3,3)
c
c                        strain
      real(8) :: dFn(3,3),Fnh(3,3),FpFinv(3,3)
      real(8) :: L33(3,3),D2(3,3), Uinv(3,3)
c
      real(8), parameter :: zero = 0.0D0, one = 1.0D0
      real(8), parameter :: half = 0.5D0, two = 2.0D0

c                         indices array
      indices = 0
      indices(1,1:27) = 1
      indices(1,28:54) = 2
      indices(1,55:81) = 3
c
      indices(2,1:9) = 1
      indices(2,10:18) = 2
      indices(2,19:27) = 3
      indices(2,28:54) = indices(2,1:27)
      indices(2,55:81) = indices(2,1:27)
c
      indices(3,1:81) = 1
      indices(3,4:81:9) = 2
      indices(3,5:81:9) = 2
      indices(3,6:81:9) = 2
      indices(3,7:81:9) = 3
      indices(3,8:81:9) = 3
      indices(3,9:81:9) = 3
c
      indices(4,1:81:3) = 1
      indices(4,2:81:3) = 2
      indices(4,3:81:3) = 3
c
c                         frequently used tensor
      call cep2A_transpose33(R,transR)
      call cep2A_multiply33(transR,Fn1,U)
      traceU = U(1,1) + U(2,2) + U(3,3)
      Y = zero
      Y(1,1) = traceU
      Y(2,2) = traceU
      Y(3,3) = traceU
      Y = Y - U
      call cep2A_multiply33(R,t,temp)
      call cep2A_multiply33(temp,transR,sigma)
      call cep2A_multiply33(finv,R,Uinv)
      call cep2A_multiply33(t,Uinv,t_Uinv)
      call cep2A_multiply33(R,t,R_t)
      call cep2A_multiply33(R,Y,RY)
      call cep2A_multiply33(RY,transR,RYR)
c     detY_1 = one / ( traceU * traceU * traceU - detFn1 )
      call cep2A_det33(Y,detY_1)
      detY_1 = one / detY_1
      call cep2A_multiply33(R_t,Uinv,RtRF)
      call cep2A_transpose33(finv,dJdF)
      dJdF = dJdF * detFn1
      dFn = Fn1 - Fn
      Fnh = half * ( Fn1 + Fn )
      call cep2A_transpose33(Rh,transRh)
      call cep2A_multiply33(transRh,Fnh,Uh)
      traceUh = Uh(1,1) + Uh(2,2) + Uh(3,3)
      Yh = zero
      Yh(1,1) = traceUh
      Yh(2,2) = traceUh
      Yh(3,3) = traceUh
      Yh = Yh - Uh
      call cep2A_multiply33(Rh,Yh,temp)
      call cep2A_multiply33(temp,transRh,RYRh)
c     detYh_1 = one / ( traceUh * traceUh * traceUh - detFnh )
      call cep2A_det33(Yh,detYh_1)
      detYh_1 = one / detYh_1
      call cep2A_multiply33(Rh,Yh,RYh)
      FpFinv = half * fnhinv
      call cep2A_multiply33(dFn,FPFinv,L33)
      call cep2A_transpose33(L33,temp)
      D2 = L33 + temp
      
c                         cep99
      index_voigt = [1,4,6,4,2,5,6,5,3]
      cep99(1:9,1:9) = C(index_voigt,index_voigt)

c                         dR(i,j) / dF(k,l)
      dRdF = zero
      dRhdF = zero
      do m = 1,81
        i = indices(1,m)
        j = indices(2,m)
        k = indices(3,m)
        l = indices(4,m)
        dRdF(m) = detY_1 * (RYR(i,k)*Y(l,j) - RY(i,l)*RY(k,j))
        dRhdF(m) = detYh_1 * (RYRh(i,k)*Yh(l,j) - RYh(i,l)*RYh(k,j))
      end do
      dRhdF = dRhdF * half

c                         dL(i,j) / dF(k,l)
      dLdF = zero
      do m = 1,81
          i = indices(1,m)
          j = indices(2,m)
          k = indices(3,m)
          l = indices(4,m)
          dLdF(m) = -L33(i,k) * FpFinv(l,j)
          if ( i .eq. k ) dLdF(m) = dLdF(m) + FpFinv(l,j)
      end do
      
c                         dt(i,j) / dF(p,q)
      dtdF_1 = zero
      dtdF_2 = zero
      dtdF_3 = zero
      do a = 1,81
        i = indices(1,a)
        j = indices(2,a)
        p = indices(3,a)
        q = indices(4,a)
        do m = 1,3
          do n = 1,3
            tmp = (m-1)*27 + (i-1)*9 + (p-1)*3 + q
            dtdF_1(a) = dtdF_1(a) + dRhdF(tmp) * D2(m,n) * Rh(n,j)
            tmp = (m-1)*27 + (n-1)*9 + (p-1)*3 + q
            dtdF_2(a) = dtdF_2(a) + Rh(m,i) * dLdF(tmp) * Rh(n,j)
            tmp = (n-1)*27 + (m-1)*9 + (p-1)*3 + q
            dtdF_2(a) = dtdF_2(a) + Rh(m,i) * dLdF(tmp) * Rh(n,j)
            tmp = (n-1)*27 + (j-1)*9 + (p-1)*3 + q
            dtdF_3(a) = dtdF_3(a) + Rh(m,i) * D2(m,n) * dRhdF(tmp)
          end do
        end do
      end do
      dtdF = dtdF_1 + dtdF_2 + dtdF_3
      dtdF_1 = dtdF
      call cep2A_ddot44(cep99, dtdF_1, dtdF)
      
c                         loop over tensor component
      dPdF0 = zero
      dPdF1 = zero
      dPdF2 = zero
      dPdF3 = zero
      dPdF4 = zero
      do p = 1,81
        i = indices(1,p)
        j = indices(2,p)
        k = indices(3,p)
        l = indices(4,p)
        dPdF0(p) = RtRF(i,j) * dJdF(k,l)
        do m = 1,3
          dPdF1(p) = dPdF1(p) 
     &             - detFn1 * sigma(i,m) * finv(j,k) * finv(l,m)
          tmp = (i-1)*27 + (m-1)*9 + (k-1)*3 + l
          dPdF2(p) = dPdF2(p) +  detFn1 * dRdF(tmp) * t_Uinv(m,j)
          do n = 1,3
            tmp = (m-1)*27 + (n-1)*9 + (k-1)*3 + l
            dPdF3(p) = dPdF3(p)
     &               + detFn1 * R(i,m) * dtdF(tmp) * Uinv(n,j)
            tmp = (n-1)*27 + (m-1)*9 + (k-1)*3 + l
            dPdF4(p) = dPdF4(p)
     &               + detFn1 * R_t(i,m) * dRdF(tmp) * finv(j,n)
          end do
        end do
      end do
      dPdF = dPdF0 + dPdF1 + dPdF2 + dPdF3 + dPdF4
      
      return

      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *      subroutine multiply33, transpose33, ddot44, det33       *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *      frequently used mathematical subroutines by cep2A       *
c     *                                                              *
c     *                       should be inlined                      *
c     *                                                              *
c     ****************************************************************
c
      subroutine cep2A_multiply33(A, B, C)
      implicit none
      real(8), intent(in)  :: A(3,3), B(3,3)
      real(8), intent(out) :: C(3,3)

      C(1,1) = A(1,1)*B(1,1) + A(1,2)*B(2,1) + A(1,3)*B(3,1)
      C(2,1) = A(2,1)*B(1,1) + A(2,2)*B(2,1) + A(2,3)*B(3,1)
      C(3,1) = A(3,1)*B(1,1) + A(3,2)*B(2,1) + A(3,3)*B(3,1)
      C(1,2) = A(1,1)*B(1,2) + A(1,2)*B(2,2) + A(1,3)*B(3,2)
      C(2,2) = A(2,1)*B(1,2) + A(2,2)*B(2,2) + A(2,3)*B(3,2)
      C(3,2) = A(3,1)*B(1,2) + A(3,2)*B(2,2) + A(3,3)*B(3,2)
      C(1,3) = A(1,1)*B(1,3) + A(1,2)*B(2,3) + A(1,3)*B(3,3)
      C(2,3) = A(2,1)*B(1,3) + A(2,2)*B(2,3) + A(2,3)*B(3,3)
      C(3,3) = A(3,1)*B(1,3) + A(3,2)*B(2,3) + A(3,3)*B(3,3)

      return
      end subroutine
      
      subroutine cep2A_transpose33(A, B)
      implicit none
      real(8), intent(in)  :: A(3,3)
      real(8), intent(out) :: B(3,3)

      B(1,1) = A(1,1)
      B(2,1) = A(1,2)
      B(3,1) = A(1,3)
      B(1,2) = A(2,1)
      B(2,2) = A(2,2)
      B(3,2) = A(2,3)
      B(1,3) = A(3,1)
      B(2,3) = A(3,2)
      B(3,3) = A(3,3)

      return
      end subroutine

      subroutine cep2A_ddot44(A4, B4, C4)
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

      subroutine cep2A_det33(A2,detA2)
      implicit none
      real(8),intent(in) :: A2(3,3)
      real(8),intent(out) :: detA2
      real(8),parameter :: zero = 0.0D0

      detA2 = zero
      detA2 = A2(1,1) * A2(2,2) * A2(3,3)
      detA2 = detA2 + A2(1,2) * A2(2,3) * A2(3,1)
      detA2 = detA2 + A2(1,3) * A2(2,1) * A2(3,2)
      detA2 = detA2 - A2(1,1) * A2(2,3) * A2(3,2)
      detA2 = detA2 - A2(1,2) * A2(2,1) * A2(3,3)
      detA2 = detA2 - A2(1,3) * A2(2,2) * A2(3,1)
      return
      end subroutine

