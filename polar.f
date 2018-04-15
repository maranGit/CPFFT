
c     ****************************************************************
c     *                                                              *
c     *                      subroutine rtcmp1                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 06/30/91                   *
c     *                                                              *
c     *     this subroutine computes the polar decompostion of the   *
c     *     deformation gradient into the rotation tensor [R] and a  *
c     *     deformation tensor [U] for a block of solid elements     *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine rtcmp1(span,f,r)
      implicit none
      include 'param_def'
      integer :: span, i
      double precision
     &     f(mxvl,ndim,*),r(mxvl,ndim,*),ui(mxvl,nstr)
c
c                       compute the inverse of the right
c                       stretch tensor.
c
      call irscp1( span, f, ui )
c
c                       compute the rotation tensor.
c
      do i= 1,span
         r(i,1,1)= f(i,1,1)*ui(i,1)+f(i,1,2)*ui(i,2)+f(i,1,3)*ui(i,4)
         r(i,1,2)= f(i,1,1)*ui(i,2)+f(i,1,2)*ui(i,3)+f(i,1,3)*ui(i,5)
         r(i,1,3)= f(i,1,1)*ui(i,4)+f(i,1,2)*ui(i,5)+f(i,1,3)*ui(i,6)
         r(i,2,1)= f(i,2,1)*ui(i,1)+f(i,2,2)*ui(i,2)+f(i,2,3)*ui(i,4)
         r(i,2,2)= f(i,2,1)*ui(i,2)+f(i,2,2)*ui(i,3)+f(i,2,3)*ui(i,5)
         r(i,2,3)= f(i,2,1)*ui(i,4)+f(i,2,2)*ui(i,5)+f(i,2,3)*ui(i,6)
         r(i,3,1)= f(i,3,1)*ui(i,1)+f(i,3,2)*ui(i,2)+f(i,3,3)*ui(i,4)
         r(i,3,2)= f(i,3,1)*ui(i,2)+f(i,3,2)*ui(i,3)+f(i,3,3)*ui(i,5)
         r(i,3,3)= f(i,3,1)*ui(i,4)+f(i,3,2)*ui(i,5)+f(i,3,3)*ui(i,6)
      end do
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine irscp1                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 06/30/91                   *
c     *                                                              *
c     *     this subroutine computes the inverse of the right        *
c     *     stretch tensor. the computations are for a gauss         *
c     *     point for a block of solid elements                      *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine irscp1( span, f, ui )
c     implicit integer (a-z)
      implicit none
      include 'param_def'

      integer :: span, i
c
c                       parameter declarations
c
      double precision
     &   f(mxvl,ndim,*),ui(mxvl,*)
c
c                       locally allocated arrays
c
      double precision
     &   c(mxvl,6), cc(mxvl,6),
     &   iu(mxvl), iiu(mxvl), iiiu(mxvl), a2(mxvl), b2(mxvl),
     &   c2(mxvl),d2(mxvl), one, two
c
      data one, two / 1.0d00, 2.0d00 /
c
c                       ui is in symmetric upper triangular form.
c
c                       compute the invariants of the right
c                       stretch tensor, the metric tensor, and
c                       its square.
c
      call ivcmp1( span, f, c, cc, iu, iiu, iiiu )
c
c                       compute multipliers.
c
      do i = 1, span
         a2(i)= one/(iiiu(i)*(iu(i)*iiu(i)-iiiu(i)))
         b2(i)= iu(i)*iiu(i)*iiu(i)-iiiu(i)*(iu(i)*iu(i)+iiu(i))
         c2(i)= -iiiu(i)-iu(i)*(iu(i)*iu(i)-two*iiu(i))
         d2(i)= iu(i)
      end do
c
c                       compute the inverse of the right
c                       stretch tensor.
c
      do i = 1, span
         ui(i,1)= a2(i) * ( b2(i) + c2(i)*c(i,1) + d2(i)*cc(i,1) )
         ui(i,2)= a2(i) * (         c2(i)*c(i,2) + d2(i)*cc(i,2) )
         ui(i,3)= a2(i) * ( b2(i) + c2(i)*c(i,3) + d2(i)*cc(i,3) )
         ui(i,4)= a2(i) * (         c2(i)*c(i,4) + d2(i)*cc(i,4) )
         ui(i,5)= a2(i) * (         c2(i)*c(i,5) + d2(i)*cc(i,5) )
         ui(i,6)= a2(i) * ( b2(i) + c2(i)*c(i,6) + d2(i)*cc(i,6) )
      end do
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine ivcmp1                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 9/6/2016 rhd               *
c     *                                                              *
c     *     this subroutine computes the invariants of the right     *
c     *     stretch tensor, the metric tensor, and its square.       *
c     *     the computations are for a gauss point in a bblock of    *
c     *     solid elements                                           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ivcmp1( span, f, c, cc, iu, iiu, iiiu )
      implicit none
      include 'param_def'
c
c               parameter declarations
c
      integer :: span
      double precision ::
     & f(mxvl,ndim,*), c(mxvl,6), cc(mxvl,6), iu(*), iiu(*), iiiu(*)
c
c               locally allocated arrays
c
      integer :: i
      logical, parameter :: new = .true.
      double precision :: ct(mxvl,nstr), ev(mxvl,ndim)
c
c              c and cc are in symmetric upper triangular form.
c              compute the metric tensor.
c
      do i = 1, span
       c(i,1)= f(i,1,1)*f(i,1,1)+f(i,2,1)*f(i,2,1)+f(i,3,1)*f(i,3,1)
       c(i,2)= f(i,1,1)*f(i,1,2)+f(i,2,1)*f(i,2,2)+f(i,3,1)*f(i,3,2)
       c(i,3)= f(i,1,2)*f(i,1,2)+f(i,2,2)*f(i,2,2)+f(i,3,2)*f(i,3,2)
       c(i,4)= f(i,1,1)*f(i,1,3)+f(i,2,1)*f(i,2,3)+f(i,3,1)*f(i,3,3)
       c(i,5)= f(i,1,2)*f(i,1,3)+f(i,2,2)*f(i,2,3)+f(i,3,2)*f(i,3,3)
       c(i,6)= f(i,1,3)*f(i,1,3)+f(i,2,3)*f(i,2,3)+f(i,3,3)*f(i,3,3)
      end do
c
c              compute the square of the metric tensor
c
      do i = 1, span
       cc(i,1)= c(i,1)*c(i,1)+c(i,2)*c(i,2)+c(i,4)*c(i,4)
       cc(i,2)= c(i,1)*c(i,2)+c(i,2)*c(i,3)+c(i,4)*c(i,5)
       cc(i,3)= c(i,2)*c(i,2)+c(i,3)*c(i,3)+c(i,5)*c(i,5)
       cc(i,4)= c(i,1)*c(i,4)+c(i,2)*c(i,5)+c(i,4)*c(i,6)
       cc(i,5)= c(i,2)*c(i,4)+c(i,3)*c(i,5)+c(i,5)*c(i,6)
       cc(i,6)= c(i,4)*c(i,4)+c(i,5)*c(i,5)+c(i,6)*c(i,6)
      end do
c
c              old or new algorithm to get eivenvalues. old
c              uses vectroized Givens rotations to diagonalize
c              the 3x3 symmetric, real matrix. New uses closed
c              form Cardano extraction of eigenvales for this
c              specific type & size of matrix. New is 
c              considerable faster.
c
      if( new ) call evcmp1_new( span, mxvl, c, ev )
      if( .not. new ) then   
c
c              copy the metric tensor to stress vector
c              form then get principal values.
c
          do i = 1, span
             ct(i,1)= c(i,1)
             ct(i,2)= c(i,3)
             ct(i,3)= c(i,6)
             ct(i,4)= c(i,2)
             ct(i,5)= c(i,5)
             ct(i,6)= c(i,4)
          end do
          call evcmp1( span, ct, ev )
      end if
c
c              set the principal values.
c
      do i = 1, span
         ev(i,1)= sqrt(ev(i,1))
         ev(i,2)= sqrt(ev(i,2))
         ev(i,3)= sqrt(ev(i,3))
      end do
c
c              invariants of right stretch tensor.
c
      do i = 1, span
       iu(i)  = ev(i,1)+ev(i,2)+ev(i,3)
       iiu(i) = ev(i,1)*ev(i,2)+ev(i,2)*ev(i,3)+ev(i,1)*ev(i,3)
       iiiu(i)= ev(i,1)*ev(i,2)*ev(i,3)
      end do
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine evcmp1_new                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 09/6/2016                  *
c     *                                                              *
c     *     eigenvalues of metric tensor. symmetric, real 3x3
c     *                                                              *
c     ****************************************************************
c
      subroutine evcmp1_new( span, mxvl, c, lamda )
      implicit none
c
c                 parameter declarations
c
      integer :: span, mxvl
      double precision :: c(mxvl,6), lamda(mxvl,3)
c
c                 locals
c
      integer :: bel
      double precision ::
     &  m11, m12, m13, m22, m23, m33, e1, e2, e3,
     &  swap1, swap2, swap3, zero,
     &  de, dd, ee, ff, m, c1, c0,p, q, sqrtp,phi, cphi, sphi,
     &  one, two, three, third, oneptfive, thirteenptfive,
     &  twentyseven, quarter, sixpt75, oneroot3
      data zero, one, two, three / 0.0d0, 1.0d0, 2.0d0, 3.0d0 /
      data third, oneptfive / 0.3333333333333333333d0, 1.5d0 /
      data thirteenptfive, twentyseven / 13.5d0, 27.0d0 /
      data quarter, sixpt75, oneroot3
     &     / 0.25d0, 6.75d0, 0.5773502691896258d0 /
c      
c
c              calculates the eigenvalues of a symmetric 3x3 matrix 
c              using Cardano's analytical algorithm.
c              Only the diagonal and upper triangular parts of matrix
c              are accessed. The access is read-only.
c              Copyright (C) 2006  Joachim Kopp. Avialble under 
c              GNU Lesser General Public License
c
      do bel = 1, span
       m11 = c(bel,1)
       m12 = c(bel,2)
       m13 = c(bel,4)
       m22 = c(bel,3)
       m23 = c(bel,5)
       m33 = c(bel,6)
       de  = m12 * m23
       dd  = m12**2
       ee  = m23**2
       ff  = m13**2
       m   = m11 + m22 + m33
       c1  = ( m11*m22 + m11*m33 + m22*m33 ) - (dd + ee + ff)
       c0  = m33*dd + m11*ee + m22*ff - m11*m22*m33 - two * m13*de
       p   = m*m - three * c1
       q   = m*(p - oneptfive*c1) - thirteenptfive*c0
       sqrtp = sqrt(abs(p))
       phi = twentyseven * ( quarter * c1*c1 * (p - c1)
     &          + c0 * (q + sixpt75 * c0) )
       phi = third * atan2(sqrt(abs(phi)), q)
       cphi = sqrtp * cos(phi)
       sphi = oneroot3 * sqrtp * sin(phi)
       e2 = third * (m - cphi)
       e3 = e2 + sphi
       e1 = e2 + cphi
       e2 = e2 - sphi
c       
       if( e2 .lt. e1 ) then
            swap1 = e1
            e1    = e2
            e2    = swap1
       end if
c
       if( e3 .lt. e1 ) then
            swap2 = e1
            e1    = e3
            e3    = swap2
       end if
c
       if( e3 .lt. e2 ) then
            swap3 = e2
            e2    = e3
            e3    = swap3
       end if
c         
       lamda(bel,1) = e1
       lamda(bel,2) = e2
       lamda(bel,3) = e3
c
      end do
c       
      return
      end
      

c     ****************************************************************
c     *                                                              *
c     *                      subroutine evcmp1                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 09/5/2016 rhd              *
c     *                                                              *
c     *     eigenvalues of a 3x3 positive definite matrix in         *
c     *     stress vector form for a block of solid elements.        *
c     *                                                              *
c     ****************************************************************
c
      subroutine evcmp1( span, k, lamda )
c     implicit integer (a-z)
      implicit none
      include 'param_def'
c
c                 parameter declarations
c
      double precision
     &  k(mxvl,6), lamda(mxvl,3)
c
c                 local arrays allocated
c
      double precision
     &  m(mxvl,ndim),kbari(mxvl),
     &  kbarj(mxvl), kbar(mxvl), ki(mxvl), kj(mxvl),mi(mxvl),
     &  mj(mxvl), scale(mxvl), alpha(mxvl), gamma(mxvl),x(mxvl),
     &  xsign(mxvl), rad(mxvl), errork(mxvl), swap(mxvl),
     &  ratiok(mxvl), sqtol, thold
      integer iexp(mxvl)
      integer :: bel, maxswp, span, swpnum
      logical cvgtst
      double precision
     &  jactol, one, four, ten, ten_thouth, zero, two
      data maxswp/15/,zero, one, two, jactol, four, ten, ten_thouth
     &   / 0.0d00, 1.0d00, 2.0d00, 1.0d-08,
     &     4.0d00, 10.0d00, 0.0001d00 /
c
c              initialize lamda, m, sweep parameters.
c
      swpnum = 0
c
      do bel = 1, span
c
         m(bel,1)= one
         m(bel,2)= one
         m(bel,3)= one
         lamda(bel,1) = k(bel,1)
         lamda(bel,2) = k(bel,2)
         lamda(bel,3) = k(bel,3)
c
c              scale [k] and ^0m^2 to avoid problems with exponential
c              overflow and underflow.
c
c              find the max and min terms on the diagonal of [k] & ^0m^2
c
         kj(bel) = k(bel,1)
         kj(bel) = min( k(bel,2),kj(bel) )
         kj(bel) = min( k(bel,3),kj(bel) )
         ki(bel) = k(bel,1)
         ki(bel) = max( k(bel,2),ki(bel) )
         ki(bel) = max( k(bel,3),ki(bel) )
         mj(bel) = one
         mi(bel) = one
c
c              compute the scale factor and do the scaling
c
         iexp(bel) = idint( ( log10(kj(bel))+log10(ki(bel))+
     &                      log10(mj(bel))+log10(mi(bel)) ) / four )
         scale(bel) = one / ( ten ** iexp(bel) )
         m(bel,1) = m(bel,1) * scale(bel)
         m(bel,2) = m(bel,2) * scale(bel)
         m(bel,3) = m(bel,3) * scale(bel)
         k(bel,1) = k(bel,1) * scale(bel)
         k(bel,4) = k(bel,4) * scale(bel)
         k(bel,2) = k(bel,2) * scale(bel)
         k(bel,6) = k(bel,6) * scale(bel)
         k(bel,5) = k(bel,5) * scale(bel)
         k(bel,3) = k(bel,3) * scale(bel)
c
      end do
c
c              begin a new sweep
c
      do !  rotate iterations to eliminate off diagonals
c            
      swpnum = swpnum + 1
      thold  = ten_thouth ** swpnum
      sqtol  = jactol * jactol
      if( thold < sqtol ) thold = sqtol
c
c              enter sweep loop -- work on lower triangle only
c                                          ( i > j )
c
c              rows are done from top to bottom
c              columns are done from left to right.
c
c
c           ***************************************
c           *                                     *
c           *           row 2 and column 1.       *
c           *                                     *
c           ***************************************
c
c
      do bel = 1, span
c
c                       check if term is within threshold
c
         ratiok(bel) = (k(bel,4)*k(bel,4))/(k(bel,2)*k(bel,1))
         if( ratiok(bel) < thold )  cycle
c
c                      compute the rotatiom matrix:  an identity
c                      matrix with alpha at position (2,1) and
c                      gamma at position (1,2).
c
            kbari(bel) = -m(bel,2)*k(bel,4)
            kbarj(bel) = -m(bel,1)*k(bel,4)
            kbar(bel)  = k(bel,2)*m(bel,1)-k(bel,1)*m(bel,2)
            rad(bel)   = (kbar(bel)*kbar(bel)/four) + 
     &                   kbari(bel)*kbarj(bel)
c
            xsign(bel) = one
            x(bel) = kbar(bel)/two+sign(xsign(bel),kbar(bel))*
     &               sqrt(rad(bel))
c
            if( (abs(x(bel)).lt.jactol*abs(kbarj(bel))).or.
     &          (abs(x(bel)).lt.jactol*abs(kbari(bel)))    ) then
               alpha(bel) = zero
               gamma(bel) = -k(bel,4)/k(bel,2)
            else
               alpha(bel) = kbarj(bel)/x(bel)
               gamma(bel) = -kbari(bel)/x(bel)
            end if
c
c                       perform the rotation.
c
c                       row 3, column 2
c                       row 3, column 1
c
            ki(bel)  = k(bel,5)
            kj(bel)  = k(bel,6)
            k(bel,5) = ki(bel)+gamma(bel)*kj(bel)
            k(bel,6) = kj(bel)+alpha(bel)*ki(bel)
c
c                       term (2,1) and diagonal terms (2,2) and (1,1).
c
            kj(bel)  = k(bel,1)
            mj(bel)  = m(bel,1)
            ki(bel)  = k(bel,2)
            mi(bel)  = m(bel,2)
            k(bel,1) = kj(bel)+alpha(bel)*alpha(bel)*ki(bel)+
     &                 two*alpha(bel)*k(bel,4)
            m(bel,1) = mj(bel)+alpha(bel)*alpha(bel)*mi(bel)
            k(bel,2) = ki(bel)+gamma(bel)*gamma(bel)*kj(bel)+
     &                 two*gamma(bel)*k(bel,4)
            m(bel,2) = mi(bel)+gamma(bel)*gamma(bel)*mj(bel)
            k(bel,4) = zero
c
      end do
c
c           ***************************************
c           *                                     *
c           *           row 3 and column 1.       *
c           *                                     *
c           ***************************************
c
      do bel = 1, span
c
c                       check if term is within threshold
c
         ratiok(bel) = (k(bel,6)*k(bel,6))/(k(bel,3)*k(bel,1))
         if( ratiok(bel) < thold ) cycle
c
c                       compute the rotatiom matrix:  an identity
c                       matrix with alpha at position (3,1) and
c                       gamma at position (1,3).
c
            kbari(bel) = -m(bel,3)*k(bel,6)
            kbarj(bel) = -m(bel,1)*k(bel,6)
            kbar(bel)  = k(bel,3)*m(bel,1)-k(bel,1)*m(bel,3)
            rad(bel)   = (kbar(bel)*kbar(bel)/four) + 
     &                   kbari(bel)*kbarj(bel)
c
            xsign(bel) = one
            x(bel) = kbar(bel)/two+sign(xsign(bel),kbar(bel))*
     &              sqrt(rad(bel))
            if( (abs(x(bel)).lt.jactol*abs(kbarj(bel))).or.
     &          (abs(x(bel)).lt.jactol*abs(kbari(bel)))    ) then
               alpha(bel) = zero
               gamma(bel) = -k(bel,6)/k(bel,3)
            else
               alpha(bel) =  kbarj(bel)/x(bel)
               gamma(bel) = -kbari(bel)/x(bel)
            end if
c
c                       perform the rotation.
c
c                       row 3, column 2
c                       row 2, column 1
c
            ki(bel)  = k(bel,5)
            kj(bel)  = k(bel,4)
            k(bel,5) = ki(bel)+gamma(bel)*kj(bel)
            k(bel,4) = kj(bel)+alpha(bel)*ki(bel)
c
c                       term (3,1) and diagonal terms (3,3) and (1,1).
c
            kj(bel)  = k(bel,1)
            mj(bel)  = m(bel,1)
            ki(bel)  = k(bel,3)
            mi(bel)  = m(bel,3)
            k(bel,1) = kj(bel)+alpha(bel)*alpha(bel)*ki(bel)+
     &                 two*alpha(bel)*k(bel,6)
            m(bel,1) = mj(bel)+alpha(bel)*alpha(bel)*mi(bel)
            k(bel,3) = ki(bel)+gamma(bel)*gamma(bel)*kj(bel)+
     &                 two*gamma(bel)*k(bel,6)
            m(bel,3) = mi(bel)+gamma(bel)*gamma(bel)*mj(bel)
            k(bel,6) = zero
c
      end do
c
c           ***************************************
c           *                                     *
c           *           row 3 and column 2.       *
c           *                                     *
c           ***************************************
c
      do bel = 1, span
c
c                       check if term is within threshold
c
         ratiok(bel) = (k(bel,5)*k(bel,5))/(k(bel,3)*k(bel,2))
         if( ratiok(bel) < thold ) cycle
c
c                       compute the rotatiom matrix:  an identity
c                       matrix with alpha at position (3,2) and
c                       gamma at position (2,3).
c
            kbari(bel) = -m(bel,3)*k(bel,5)
            kbarj(bel) = -m(bel,2)*k(bel,5)
            kbar(bel)  =  k(bel,3)*m(bel,2)-k(bel,2)*m(bel,3)
            rad(bel)   = (kbar(bel)*kbar(bel)/four)
     &                   + kbari(bel)*kbarj(bel)
c
            xsign(bel) = one
            x(bel) = kbar(bel)/two+sign(xsign(bel),kbar(bel))*
     &               sqrt(rad(bel))
            if( (abs(x(bel)).lt.jactol*abs(kbarj(bel))).or.
     &          (abs(x(bel)).lt.jactol*abs(kbari(bel)))    ) then
               alpha(bel) = zero
               gamma(bel) = -k(bel,5)/k(bel,3)
            else
               alpha(bel) =  kbarj(bel)/x(bel)
               gamma(bel) = -kbari(bel)/x(bel)
            end if
c
c                       perform the rotation.
c
c                       row 3, column 1
c                       row 2, column 1
c
            ki(bel)  = k(bel,6)
            kj(bel)  = k(bel,4)
            k(bel,6) = ki(bel)+gamma(bel)*kj(bel)
            k(bel,4) = kj(bel)+alpha(bel)*ki(bel)
c
c                       term (3,2) and diagonal terms (3,3) and (2,2).
c
            kj(bel)  = k(bel,2)
            mj(bel)  = m(bel,2)
            ki(bel)  = k(bel,3)
            mi(bel)  = m(bel,3)
            k(bel,2) = kj(bel)+alpha(bel)*alpha(bel)*ki(bel)+
     &                 two*alpha(bel)*k(bel,5)
            m(bel,2) = mj(bel)+alpha(bel)*alpha(bel)*mi(bel)
            k(bel,3) = ki(bel)+gamma(bel)*gamma(bel)*kj(bel)+
     &                 two*gamma(bel)*k(bel,5)
            m(bel,3) = mi(bel)+gamma(bel)*gamma(bel)*mj(bel)
            k(bel,5) = zero
c
      end do
c
c              end sweep
c
c              check off-diagonal elements for convergence
c
      cvgtst = .true.
c
      do bel = 1, span
c
         errork(bel) = k(bel,4)*k(bel,4)/(k(bel,2)*k(bel,1))
         if( errork(bel) .gt. sqtol ) then
           cvgtst = .false.
           exit
         end if  
c
         errork(bel) = k(bel,6)*k(bel,6)/(k(bel,3)*k(bel,1))
         if( errork(bel) .gt. sqtol ) then
                 cvgtst = .false.
                 exit
         end if
c
         errork(bel) = k(bel,5)*k(bel,5)/(k(bel,3)*k(bel,2))
         if( errork(bel) .gt. sqtol ) then
            cvgtst = .false.
            exit
         end if
c
      end do
c
      if( cvgtst ) exit
      if( swpnum .lt. maxswp ) cycle
c      
      end do   ! over rotation iterations
c
c              update eigenvalue vector 
c
      do bel = 1, span
         lamda(bel,1) = k(bel,1) / m(bel,1)
         lamda(bel,2) = k(bel,2) / m(bel,2)
         lamda(bel,3) = k(bel,3) / m(bel,3)
      end do
c
c             reorder the eigenvalues. small to big
c
c
      do bel = 1, span
c
         if( lamda(bel,2) .lt. lamda(bel,1) ) then
            swap(bel)    = lamda(bel,1)
            lamda(bel,1) = lamda(bel,2)
            lamda(bel,2) = swap(bel)
         end if
c
         if( lamda(bel,3) .lt. lamda(bel,1) ) then
            swap(bel)    = lamda(bel,1)
            lamda(bel,1) = lamda(bel,3)
            lamda(bel,3) = swap(bel)
         end if
c
         if( lamda(bel,3) .lt. lamda(bel,2) ) then
            swap(bel)    = lamda(bel,2)
            lamda(bel,2) = lamda(bel,3)
            lamda(bel,3) = swap(bel)
         end if
c
      end do
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine getrm1                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 08/23/12 rhd               *
c     *                                                              *
c     *     computes the rotation matrix taking one                  *
c     *     tensor to its corresponding value obtained by adding or  *
c     *     removing the material rotation, for a gauss point of an  *
c     *     element in a block of similar solid elements             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine getrm1( span, q, r, opt )
      implicit none
      include 'param_def'
c
c           parameter declarations
c
      integer :: span, opt
      double precision :: q(mxvl,nstr,*), r(mxvl,ndim,*)
     &  
c
c           locals
c
      integer :: i
      double precision
     & two, rbar(mxvl,3,3)
      data two / 2.0d00 /
c
c           compute q. branch on quantity & direction of rotation.
c
      if ( opt .eq. 1 ) then
c
c      unrotated rate of deformation vector {d} = [q] *
c       (rotated) rate of deformation vector {D}. in tensor form:
c
c                 [d] = trans([R]) [D] [R]
c
c       both [D] and [d] are symmetric, [R] is orthogonal rotation.
c       vector forms for {d} and {D} use engineering shear strains.
c       vector ordering is {x,y,z,xy,yz,xz}
c

         do i = 1, span
            q(i,1,1)= r(i,1,1)**2
            q(i,1,2)= r(i,2,1)**2
            q(i,1,3)= r(i,3,1)**2
            q(i,1,4)= r(i,1,1)*r(i,2,1)
            q(i,1,5)= r(i,3,1)*r(i,2,1)
            q(i,1,6)= r(i,1,1)*r(i,3,1)
            q(i,2,1)= r(i,1,2)**2
            q(i,2,2)= r(i,2,2)**2
            q(i,2,3)= r(i,3,2)**2
            q(i,2,4)= r(i,1,2)*r(i,2,2)
            q(i,2,5)= r(i,3,2)*r(i,2,2)
            q(i,2,6)= r(i,1,2)*r(i,3,2)
            q(i,3,1)= r(i,1,3)**2
            q(i,3,2)= r(i,2,3)**2
            q(i,3,3)= r(i,3,3)**2
            q(i,3,4)= r(i,1,3)*r(i,2,3)
            q(i,3,5)= r(i,3,3)*r(i,2,3)
            q(i,3,6)= r(i,1,3)*r(i,3,3)
            q(i,4,1)= two*r(i,1,1)*r(i,1,2)
            q(i,4,2)= two*r(i,2,1)*r(i,2,2)
            q(i,4,3)= two*r(i,3,1)*r(i,3,2)
            q(i,4,4)= r(i,1,1)*r(i,2,2)+r(i,1,2)*r(i,2,1)
            q(i,4,5)= r(i,2,1)*r(i,3,2)+r(i,3,1)*r(i,2,2)
            q(i,4,6)= r(i,1,1)*r(i,3,2)+r(i,3,1)*r(i,1,2)
            q(i,5,1)= two*r(i,1,2)*r(i,1,3)
            q(i,5,2)= two*r(i,2,3)*r(i,2,2)
            q(i,5,3)= two*r(i,3,2)*r(i,3,3)
            q(i,5,4)= r(i,1,2)*r(i,2,3)+r(i,2,2)*r(i,1,3)
            q(i,5,5)= r(i,2,2)*r(i,3,3)+r(i,2,3)*r(i,3,2)
            q(i,5,6)= r(i,1,2)*r(i,3,3)+r(i,3,2)*r(i,1,3)
            q(i,6,1)= two*r(i,1,1)*r(i,1,3)
            q(i,6,2)= two*r(i,2,1)*r(i,2,3)
            q(i,6,3)= two*r(i,3,1)*r(i,3,3)
            q(i,6,4)= r(i,1,1)*r(i,2,3)+r(i,2,1)*r(i,1,3)
            q(i,6,5)= r(i,2,1)*r(i,3,3)+r(i,3,1)*r(i,2,3)
            q(i,6,6)= r(i,1,1)*r(i,3,3)+r(i,1,3)*r(i,3,1)
        end do
        return
      end if
c
      if ( opt .eq. 2 ) then
c
c       cauchy stress {T} = [q] * (rotated) cauchy stress {t}.
c       in tensor form:
c
c                 [T] = [R] [t] trans([R])
c
c       both [T] and [t] are symmetric, [R] is orthogonal rotation.
c       vector ordering is {x,y,z,xy,yz,xz}. this [q] matrix
c       is the transpose of the one above.
c

         do i = 1, span
            q(i,1,1)= r(i,1,1)**2
            q(i,1,2)= r(i,1,2)**2
            q(i,1,3)= r(i,1,3)**2
            q(i,1,4)= two*r(i,1,1)*r(i,1,2)
            q(i,1,5)= two*r(i,1,3)*r(i,1,2)
            q(i,1,6)= two*r(i,1,1)*r(i,1,3)
            q(i,2,1)= r(i,2,1)**2
            q(i,2,2)= r(i,2,2)**2
            q(i,2,3)= r(i,2,3)**2
            q(i,2,4)= two*r(i,2,1)*r(i,2,2)
            q(i,2,5)= two*r(i,2,3)*r(i,2,2)
            q(i,2,6)= two*r(i,2,1)*r(i,2,3)
            q(i,3,1)= r(i,3,1)**2
            q(i,3,2)= r(i,3,2)**2
            q(i,3,3)= r(i,3,3)**2
            q(i,3,4)= two*r(i,3,1)*r(i,3,2)
            q(i,3,5)= two*r(i,3,3)*r(i,3,2)
            q(i,3,6)= two*r(i,3,1)*r(i,3,3)
            q(i,4,1)= r(i,1,1)*r(i,2,1)
            q(i,4,2)= r(i,1,2)*r(i,2,2)
            q(i,4,3)= r(i,1,3)*r(i,2,3)
            q(i,4,4)= r(i,1,1)*r(i,2,2)+r(i,2,1)*r(i,1,2)
            q(i,4,5)= r(i,1,2)*r(i,2,3)+r(i,1,3)*r(i,2,2)
            q(i,4,6)= r(i,1,1)*r(i,2,3)+r(i,1,3)*r(i,2,1)
            q(i,5,1)= r(i,2,1)*r(i,3,1)
            q(i,5,2)= r(i,3,2)*r(i,2,2)
            q(i,5,3)= r(i,2,3)*r(i,3,3)
            q(i,5,4)= r(i,2,1)*r(i,3,2)+r(i,2,2)*r(i,3,1)
            q(i,5,5)= r(i,2,2)*r(i,3,3)+r(i,3,2)*r(i,2,3)
            q(i,5,6)= r(i,2,1)*r(i,3,3)+r(i,2,3)*r(i,3,1)
            q(i,6,1)= r(i,1,1)*r(i,3,1)
            q(i,6,2)= r(i,1,2)*r(i,3,2)
            q(i,6,3)= r(i,1,3)*r(i,3,3)
            q(i,6,4)= r(i,1,1)*r(i,3,2)+r(i,1,2)*r(i,3,1)
            q(i,6,5)= r(i,1,2)*r(i,3,3)+r(i,1,3)*r(i,3,2)
            q(i,6,6)= r(i,1,1)*r(i,3,3)+r(i,3,1)*r(i,1,3)
         end do
         return
      end if
c
c
      if ( opt .eq. 3 ) then
c
c       unrotated cauchy stress {t} = [q] * cauchy stress {T}.
c       in tensor form:
c
c                 [t] = trans([R]) [T] [R]
c
c       want to use code above for opt = 2. Set rbar = trans([R])
c       and compute [q]. We are computing the
c
c       both [T] and [t] are symmetric, [R] is orthogonal rotation.
c       vector ordering is {x,y,z,xy,yz,xz}. this [q] matrix
c       is the transpose of the one above.
c

         do i = 1, span
          rbar(i,1,1) = r(i,1,1)
          rbar(i,1,2) = r(i,2,1)
          rbar(i,1,3) = r(i,3,1)
          rbar(i,2,1) = r(i,1,2)
          rbar(i,2,2) = r(i,2,2)
          rbar(i,2,3) = r(i,3,2)
          rbar(i,3,1) = r(i,1,3)
          rbar(i,3,2) = r(i,2,3)
          rbar(i,3,3) = r(i,3,3)
         end do
c

         do i = 1, span
           q(i,1,1)= rbar(i,1,1)**2
           q(i,1,2)= rbar(i,1,2)**2
           q(i,1,3)= rbar(i,1,3)**2
           q(i,1,4)= two*rbar(i,1,1)*rbar(i,1,2)
           q(i,1,5)= two*rbar(i,1,3)*rbar(i,1,2)
           q(i,1,6)= two*rbar(i,1,1)*rbar(i,1,3)
           q(i,2,1)= rbar(i,2,1)**2
           q(i,2,2)= rbar(i,2,2)**2
           q(i,2,3)= rbar(i,2,3)**2
           q(i,2,4)= two*rbar(i,2,1)*rbar(i,2,2)
           q(i,2,5)= two*rbar(i,2,3)*rbar(i,2,2)
           q(i,2,6)= two*rbar(i,2,1)*rbar(i,2,3)
           q(i,3,1)= rbar(i,3,1)**2
           q(i,3,2)= rbar(i,3,2)**2
           q(i,3,3)= rbar(i,3,3)**2
           q(i,3,4)= two*rbar(i,3,1)*rbar(i,3,2)
           q(i,3,5)= two*rbar(i,3,3)*rbar(i,3,2)
           q(i,3,6)= two*rbar(i,3,1)*rbar(i,3,3)
           q(i,4,1)= rbar(i,1,1)*rbar(i,2,1)
           q(i,4,2)= rbar(i,1,2)*rbar(i,2,2)
           q(i,4,3)= rbar(i,1,3)*rbar(i,2,3)
           q(i,4,4)= rbar(i,1,1)*rbar(i,2,2)+rbar(i,2,1)*rbar(i,1,2)
           q(i,4,5)= rbar(i,1,2)*rbar(i,2,3)+rbar(i,1,3)*rbar(i,2,2)
           q(i,4,6)= rbar(i,1,1)*rbar(i,2,3)+rbar(i,1,3)*rbar(i,2,1)
           q(i,5,1)= rbar(i,2,1)*rbar(i,3,1)
           q(i,5,2)= rbar(i,3,2)*rbar(i,2,2)
           q(i,5,3)= rbar(i,2,3)*rbar(i,3,3)
           q(i,5,4)= rbar(i,2,1)*rbar(i,3,2)+rbar(i,2,2)*rbar(i,3,1)
           q(i,5,5)= rbar(i,2,2)*rbar(i,3,3)+rbar(i,3,2)*rbar(i,2,3)
           q(i,5,6)= rbar(i,2,1)*rbar(i,3,3)+rbar(i,2,3)*rbar(i,3,1)
           q(i,6,1)= rbar(i,1,1)*rbar(i,3,1)
           q(i,6,2)= rbar(i,1,2)*rbar(i,3,2)
           q(i,6,3)= rbar(i,1,3)*rbar(i,3,3)
           q(i,6,4)= rbar(i,1,1)*rbar(i,3,2)+rbar(i,1,2)*rbar(i,3,1)
           q(i,6,5)= rbar(i,1,2)*rbar(i,3,3)+rbar(i,1,3)*rbar(i,3,2)
           q(i,6,6)= rbar(i,1,1)*rbar(i,3,3)+rbar(i,3,1)*rbar(i,1,3)
         end do
         return
      end if

      write(*,*) ' ** Fatal Error: call to getrm1 with opt /= 1,2,3'
      write(*,*) '    is invalid. program terminated'
      call die_abort
      stop
c
      end
