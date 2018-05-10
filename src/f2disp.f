c     ****************************************************************
c     *                                                              *
c     *                      subroutine f2disp                       *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified : 5/07/2018 RM               *
c     *                                                              *
c     *           compute nodal displacement from deformation        *
c     *           gradient using finite element shape function.      *
c     *           The weak form makes sure that grad u = F.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine f2disp(step)
      use fft, only: Fn1, N, incmap, incid, N3
      implicit none
      include 'common.main'
c
c                          global
c
      integer :: step
c
c                          local
c
      integer :: ndofs, ndofs_1, nen, nen1, currElem, currNode
      real(8), allocatable, dimension(:,:) :: AA, bb, x_current
      logical :: debug

      integer :: ii, kk, mm, nn, ll
      logical :: bf, der
      integer :: ix_local(8), lint
      real(8) :: Coordinate_local(3,8), A_local(8,8), b_local(8,3)
      real(8) :: F_local_T(3,3), B_matrix(3,8), BT_matrix(8,3)
      real(8) :: Wgt, ss(3) ! intpntb
      real(8) :: shl(1,8), shld(3,8), shls(6,8), be(4) ! shlb
      real(8) :: QXY(3,8), shgs(6,8), Jdet, xs(3,3) ! shgb
      real(8), parameter :: zero = 0.0D0
c
c                        problem size
c
      debug   = .false.
      nen     = 8
      nen1    = 9
      ndofs   = ( N + 1 ) * ( N + 1 ) * ( N + 1 )
      ndofs_1 = ndofs - 1
      allocate( AA(ndofs_1,ndofs_1), bb(ndofs_1,3) )
c
c                  global stiffness and force vector
c
      AA = zero
      bb = zero
      bf = .false.
      der = .false.

      do ii = 1, N3

c                 8-point integration to get A matrix
        lint = 8
        currElem = incmap(ii)
        do kk = 1, nen
          ix_local(kk) = incid(currElem+kk-1)
          currNode = ix_local(kk) * 3 - 2
          Coordinate_local( 1, kk ) = c( currNode )
          Coordinate_local( 2, kk ) = c( currNode + 1 )
          Coordinate_local( 3, kk ) = c( currNode + 2 )
        end do

        A_local = zero
        b_local = zero

        if(debug) write(out,9001)

        do ll = 1, lint

          ! get value of B at integration point
          call intpntb(ll,lint,0,Wgt,ss)
          call shlb(ss,nen,nen,der,bf,shl,shld,shls,be)
          call shgb(Coordinate_local,nen,shld,shls,nen,nen,bf,der,
     &              Jdet,QXY,shgs,be,xs)

          ! B matrix
          B_matrix = QXY
          BT_matrix = transpose(QXY)

          ! Gauss integration to get local A matrix
          A_local = A_local + Wgt * Jdet * matmul(BT_matrix, B_matrix)

        end do

        if(debug) write(out,9002)
c
c                 1-point integration to get b vector
c
        lint = 1

        if(debug) write(out,9003)

        do ll = 1, lint

          ! get value of B at integration point
          call intpntb(ll,lint,0,Wgt,ss)
          call shlb(ss,nen,nen,der,bf,shl,shld,shls,be)
          call shgb(Coordinate_local,nen,shld,shls,nen,nen,bf,der,
     &              Jdet,QXY,shgs,be,xs)

          ! B matrix
          B_matrix = QXY
          BT_matrix = transpose(QXY)

          ! Gauss integration to get local b vector
          F_local_T(1,1) = Fn1(ii,1)
          F_local_T(2,1) = Fn1(ii,2)
          F_local_T(3,1) = Fn1(ii,3)
          F_local_T(1,2) = Fn1(ii,4)
          F_local_T(2,2) = Fn1(ii,5)
          F_local_T(3,2) = Fn1(ii,6)
          F_local_T(1,3) = Fn1(ii,7)
          F_local_T(2,3) = Fn1(ii,8)
          F_local_T(3,3) = Fn1(ii,9)
          b_local = b_local + Wgt * Jdet * matmul(BT_matrix, F_local_T)

        end do

        if(debug) write(out,9004)
c
c                           assemble
c
        do kk = 1, nen
          mm = ix_local(kk) - 1
          if( mm .eq. 0 ) cycle
          bb(mm,1:3) = bb(mm,1:3) + b_local(kk,1:3)
          do ll = 1, nen
            nn = ix_local(ll) - 1
            if( nn .eq. 0 ) cycle
            AA(mm,nn) = AA(mm,nn) + A_local(kk,ll)
          end do
        end do
        
      end do ! loop over element
c
c                    solve for displacement
c
      allocate( x_current( ndofs_1, 3 ) )
      x_current = zero
      call drive_pardiso_solver(AA(1,1),bb(1,1),ndofs_1,x_current(1,1))
c
c              store nodal displacement in common.main
c
      u(1:3) = zero
      do ii = 1, ndofs_1
        u(ii*3+1) = x_current(ii,1) - c(ii*3+1)
        u(ii*3+2) = x_current(ii,2) - c(ii*3+2)
        u(ii*3+3) = x_current(ii,3) - c(ii*3+3)
      end do

      deallocate( AA, bb, x_current )
      return
 9001 format(1x,'>>> Entering A construction',/)
 9002 format(1x,'>>> Leaving A construction',/)
 9003 format(1x,'>>> Entering b construction',/)
 9004 format(1x,'>>> Leaving b construction',/)
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine  drive_pardiso_solver             *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 5/2/18                      *
c     *                                                              *
c     *                  solve for nodal displacement                *
c     *                  three linear equations                      *
c     *                                                              *
c     ****************************************************************
c
      subroutine drive_pardiso_solver(A, b, neq, x)
      implicit none
c
c                          global
c
      integer :: neq
      real(8) :: A(neq,*), b(neq,*), x(neq,*)
c
c                          local
c
      integer(8) :: pt(64) ! long int for 64 bit machine
      integer :: maxfct, mnum, mtype, phase, nrhs, iparm(64)
      integer :: msglvl, error, A_size
      integer, allocatable, dimension(:) :: ia, ja
      real(8), allocatable :: A_CSR(:)
      integer :: idum, i, j, k
      real(8) :: ddum
      logical :: local_debug
      logical, external :: isNonZero

      local_debug = .false.

c               use CSR format to store A matrix
      pt(1:64) = 0
      maxfct = 1
      mnum = 1
      mtype = 2 ! symmetric positive definite
      nrhs = 3 ! 3 set of linear equations
      msglvl = 0
      error = 0

      iparm(1:64) = 0 ! use default options

c              number of non-zero terms in A (upper triangular)
      if(local_debug) write(*,*) ">>> Counting non-zero terms in A"

      A_size = neq
      do i = 1, ( neq - 1 )
        do j = ( i + 1 ), neq
          if( isNonZero( A(i,j) ) ) A_size = A_size + 1
        end do
      end do

      allocate( A_CSR(A_size), ia(neq+1), ja(A_size) )

      if(local_debug) write(*,*) ">>> counting finished"
c
c             use CSR format to store A matrix
c
      if(local_debug) write(*,*) ">>> condensing A matrix"

      k = 1

      do i = 1, ( neq - 1 )

        ja(k) = i
        A_CSR(k) = A(i,i)
        ia(i) = k
        k = k + 1

        do j = ( i + 1 ), neq

          if( isNonZero( A(i,j) ) ) then
            ja(k) = j
            A_CSR(k) = A(i,j)
            k = k + 1
          end if

        end do

      end do

      ja(A_size) = neq
      A_CSR(A_size) = A(neq, neq)
      ia(neq) = A_size
      ia(neq+1) = A_size + 1

      if(local_debug) write(*,*) ">>> condensing finished"
c
c             pardiso solver
c
      call pardisoinit(pt, mtype, iparm)

      if(local_debug) write(*,*) ">>> Entering pardiso solver"
      phase = 13
      call pardiso(pt, maxfct, mnum, mtype, phase, neq, A_CSR(1), ia,
     &             ja, idum, nrhs, iparm, msglvl, b(1,1), x(1,1), error)
      if(local_debug) write(*,*) ">>> Leaving pardiso solver"
      phase = -1
      if(local_debug) write(*,*) ">>> Now releasing memory"
      call pardiso(pt, maxfct, mnum, mtype, phase, neq, ddum, idum,
     &             idum, idum, nrhs, iparm, msglvl, ddum, ddum, error)

      deallocate( A_CSR, ia, ja )

      return
      end subroutine ! drive_pardiso_solver
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine intpntb                      *
c     *                      subroutine shlb                         *
c     *                      subroutine shgb                         *
c     *                                                              *
c     *                   last modified: 5/2/18  RM                  *
c     *                                                              *
c     *      generate shape function, integration point, ...         *
c     *                                                              *
c     ****************************************************************
C-----------------------------------------------------------------------
      subroutine intpntb(l,lint,ib,w,ss)
C-----------------------------------------------------------------------
C      Purpose: Gauss quadrature for 3-d Brick element
c         Revised by Timothy Truster, 02/20/2012
C      Inputs:
C         ll     - Number of points/direction
C      Outputs:
C         lint   - Total number of quadrature points
C         s(4,*) - Gauss points (1-3) and weights (4)
C-----------------------------------------------------------------------
      implicit none

      integer   l1,l2,lint, i,j,k,l, ib
      real*8    s(4,125),gauss(10),weight(10)
      real*8    ss(3),w,i1,j1,one,two,three,g
      real*8    ls(9),lt(9)

      integer iin,iout,irsin,irsout
      common /iounit/ iin,iout,irsin,irsout

      data ls /-1.d0,1.d0,1.d0,-1.d0,0.d0,1.d0,0.d0,-1.d0,0.d0/
      data lt /-1.d0,-1.d0,1.d0,1.d0,-1.d0,0.d0,1.d0,0.d0,0.d0/

        one    = 1.0d0
        two    = 2.0d0
        three  = 3.0d0

      if(ib.eq.0) then

c     1 pt. quadrature
      if(lint.eq.1) then

        do i = 1,3
          ss(i) = 0.d0
        end do ! i
        w = 8.0d0

c     2 x 2 x 2 pt. quadrature
      elseif(lint.eq.8)         then
        l1 = 2
        gauss(1)  =-0.577350269189626D0
        gauss(2)  = 0.577350269189626D0
        weight(1) = 1.D0
        weight(2) = 1.D0
        i1 = dble(l)/dble(l1**2)
        i = ceiling(i1)
        l2 = (l-(i-1)*l1**2)
        j1 = dble(l2)/dble(l1)
        j = ceiling(j1)
        k = (l2-(j-1)*l1)
        ss(1)=gauss(k)
        ss(2)=gauss(j)
        ss(3)=gauss(i)
        w=weight(i)*weight(j)*weight(k)

c        for 3*3*3 pts quadrature
      elseif(lint.eq.27)         then
        l1        = 3
        gauss(1)  =-0.774596669241483d0
        gauss(2)  = 0.d0
        gauss(3)  = 0.774596669241483d0
        weight(1) = 0.555555555555556d0
        weight(2) = 0.888888888888889d0
        weight(3) = 0.555555555555556d0
        i1 = dble(l)/dble(l1**2)
        i = ceiling(i1)
        l2 = (l-(i-1)*l1**2)
        j1 = dble(l2)/dble(l1)
        j = ceiling(j1)
        k = (l2-(j-1)*l1)
        ss(1)=gauss(k)
        ss(2)=gauss(j)
        ss(3)=gauss(i)
        w=weight(i)*weight(j)*weight(k)

c        for 4*4*4 pts quadrature
      elseif(lint.eq.64)         then
        l1        = 4
        gauss(1)  =-0.861136311594953d0
        gauss(2)  =-0.339981043584856d0
        gauss(3)  = 0.339981043584856d0
        gauss(4)  = 0.861136311594953d0
        weight(1) = 0.347854845137454d0
        weight(2) = 0.652145154862546d0
        weight(3) = 0.652145154862546d0
        weight(4) = 0.347854845137454d0
        i1 = dble(l)/dble(l1**2)
        i = ceiling(i1)
        l2 = (l-(i-1)*l1**2)
        j1 = dble(l2)/dble(l1)
        j = ceiling(j1)
        k = (l2-(j-1)*l1)
        ss(1)=gauss(k)
        ss(2)=gauss(j)
        ss(3)=gauss(i)
        w=weight(i)*weight(j)*weight(k)

c      For 5*5*5 pts quadrature
       elseif(lint.eq.125)         then
        l1        = 5
        gauss(1)  =-0.906179845938664d0
        gauss(2)  =-0.538469310105683d0
        gauss(3)  = 0.d0
        gauss(4)  =+0.538469310105683d0
        gauss(5)  =+0.906179845938664d0
        weight(1) = 0.236926885056189d0
        weight(2) = 0.478628670449366d0
        weight(3) = 0.568888888888889d0
        weight(4) = 0.478628670449366d0
        weight(5) = 0.236926885056189d0
        i1 = dble(l)/dble(l1**2)
        i = ceiling(i1)
        l2 = (l-(i-1)*l1**2)
        j1 = dble(l2)/dble(l1)
        j = ceiling(j1)
        k = (l2-(j-1)*l1)
        ss(1)=gauss(k)
        ss(2)=gauss(j)
        ss(3)=gauss(i)
        w=weight(i)*weight(j)*weight(k)

c      For 10*10*10 pts quadrature
       elseif(lint.eq.1000)         then
        l1        = 1000

        gauss(1)  = -0.973906528517172
        gauss(2)  = -0.865063366688985
        gauss(3)  = -0.679409568299024
        gauss(4)  = -0.433395394129247
        gauss(5)  = -0.148874338981631
        gauss(6)  =  0.148874338981631
        gauss(7)  =  0.433395394129247
        gauss(8)  =  0.679409568299024
        gauss(9)  =  0.865063366688985
        gauss(10) =  0.973906528517172

        weight(1) = 0.0666713443086881
        weight(2) = 0.149451349150581
        weight(3) = 0.219086362515982
        weight(4) = 0.269266719309996
        weight(5) = 0.295524224714753
        weight(6) = 0.295524224714753
        weight(7) = 0.269266719309996
        weight(8) = 0.219086362515982
        weight(9) = 0.149451349150581
        weight(10)= 0.0666713443086881

        i1 = dble(l)/dble(l1**2)
        i = ceiling(i1)
        l2 = (l-(i-1)*l1**2)
        j1 = dble(l2)/dble(l1)
        j = ceiling(j1)
        k = (l2-(j-1)*l1)
        ss(1)=gauss(k)
        ss(2)=gauss(j)
        ss(3)=gauss(i)
        w=weight(i)*weight(j)*weight(k)

      endif

      else

c
c     2x2 integration.
c     ----------------
      if (lint .eq. 4) then
        l1 = 2
        gauss(1)  =-0.577350269189626D0
        gauss(2)  = 0.577350269189626D0
        weight(1) = 1.D0
        weight(2) = 1.D0
        j1 = dble(l)/dble(l1)
        j = ceiling(j1)
        k = (l-(j-1)*l1)
        ss(1)=gauss(k)
        ss(2)=gauss(j)
        ss(3)=-1.d0
        w=weight(j)*weight(k)
      endif

      endif

      return
      end

c----------------------------------------------------------------------
c
      subroutine shlb(ss,nel,nen,der,bf,shl,shld,shls,bubble)
c
c      3 dimensional 8 noded element + bubble point
c      By Kaiming Xia 09/20/2001
c         Second derivatives of shape functions for cartesian coordinates
c         By Kaiming Xia 11/03/2001
c      * * F E A P * * A Finite Element Analysis Program
c
c      3 dimensional 27 noded element
c      By Kaiming Xia 09/20/2001
c         Second derivatives of shape functions for cartesian coordinates
c         By Kaiming Xia 11/03/2001
c
c      Revised by Ramon Calderer 03/07/2008
c      Revised by Timothy Truster, 02/20/2012
c
c----------------------------------------------------------------------
c
c      Purpose: Compute 3-d isoparametric shape
c               functions and their derivatives w/r x,y,z

c      Inputs:
c         ss(3)     - Natural coordinates of point
c         xl(ndm,*) - Nodal coordinates for element
c         ndm       - Spatial dimension of mesh
c             ideriv          - Flag for activation of second derivatives (default=0
c                     don't compute second derivative
c      Outputs:
c         xsj       - Jacobian determinant at point
c         shp(4,*)  - Shape functions and derivatives at local point
c         shp2(6,*) - Shape functions and derivatives at local point
c                     cartd(1,i) = dN_i/dx
c                     cartd(2,i) = dN_i/dy
c                     cartd(3,i) = dN_i/dz
c                     cartd(4,i) =  N_i
c                     cartd2(1,i) = d^2N_i/dx^2
c                     cartd2(2,i) = d^2N_i/dy^2
c                     cartd2(3,i) = d^2N_i/dz^2
c                     cartd2(4,i) = d^2N_i/dxy
c                     cartd2(5,i) = d^2N_i/dyz
c                     cartd2(6,i) = d^2N_i/dzx
c
c----------------------------------------------------------------------

      implicit none

!     Input Variables
      integer nel,nen
      logical der,bf
      real*8  ss(3)

!     Output Variables
      real*8 shl(1,nen),shld(3,nen),shls(6,nen),bubble(4)

!     Local Variables
      integer   i , j , k
      real*8    onPLz8,onMIz8,onPLy8,onMIy8,x,y,z

      x=ss(1)
      y=ss(2)
      z=ss(3)
      onPLz8=(1.d0+z)/8.d0
      onMIz8=(1.d0-z)/8.d0
      onPLy8=(1.d0+y)/8.d0
      onMIy8=(1.d0-y)/8.d0

      if(nel.eq.8) then
c for 8 noded 3 dimensional brick element
        shld(1,1)=-(1.d0-y)*onMIz8
        shld(2,1)=-(1.d0-x)*onMIz8
        shld(3,1)=-(1.d0-x)*onMIy8
        shl(1,1)=(1.d0-x)*(1.d0-y)*onMIz8
        shld(1,2)=(1.d0-y)*onMIz8
        shld(2,2)=-(1.d0+x)*onMIz8
        shld(3,2)=-(1.d0+x)*onMIy8
        shl(1,2)=(1.d0+x)*(1.d0-y)*onMIz8
        shld(1,3)=(1.d0+y)*onMIz8
        shld(2,3)=(1.d0+x)*onMIz8
        shld(3,3)=-(1.d0+x)*onPLy8
        shl(1,3)=(1.d0+x)*(1.d0+y)*onMIz8
        shld(1,4)=-(1.d0+y)*onMIz8
        shld(2,4)=(1.d0-x)*onMIz8
        shld(3,4)=-(1.d0-x)*onPLy8
        shl(1,4)=(1.d0-x)*(1.d0+y)*onMIz8
        shld(1,5)=-(1.d0-y)*onPLz8
        shld(2,5)=-(1.d0-x)*onPLz8
        shld(3,5)=(1.d0-x)*onMIy8
        shl(1,5)=(1.d0-x)*(1.d0-y)*onPLz8
        shld(1,6)=(1.d0-y)*onPLz8
        shld(2,6)=-(1.d0+x)*onPLz8
        shld(3,6)=(1.d0+x)*onMIy8
        shl(1,6)=(1.d0+x)*(1.d0-y)*onPLz8
        shld(1,7)=(1.d0+y)*onPLz8
        shld(2,7)=(1.d0+x)*onPLz8
        shld(3,7)=(1.d0+x)*onPLy8
        shl(1,7)=(1.d0+x)*(1.d0+y)*onPLz8
        shld(1,8)=-(1.d0+y)*onPLz8
        shld(2,8)=(1.d0-x)*onPLz8
        shld(3,8)=(1.d0-x)*onPLy8
        shl(1,8)=(1.d0-x)*(1.d0+y)*onPLz8

c        Bubble point
        if(bf) then
          bubble(1)=(-2.d0*x)*(1.d0-y*y)*(1.d0-z*z)
          bubble(2)=(1.d0-x*x)*(-2.d0*y)*(1.d0-z*z)
          bubble(3)=(1.d0-x*x)*(1.d0-y*y)*(-2.d0*z)
          bubble(4)=(1.d0-x*x)*(1.d0-y*y)*(1.d0-z*z)
        endif

c    calculate second derivatives for shape funtions in local coordinate

        if(der) then
        shls(1,1)=0.d0
        shls(2,1)=0.d0
        shls(3,1)=0.d0
        shls(4,1)=onMIz8
        shls(5,1)=(1.d0-x)/8.d0
        shls(6,1)=onMIy8
        shls(1,2)=0.d0
        shls(2,2)=0.d0
        shls(3,2)=0.d0
        shls(4,2)=-onMIz8
        shls(5,2)=(1.d0+x)/8.d0
        shls(6,2)=-onMIy8
        shls(1,3)=0.d0
        shls(2,3)=0.d0
        shls(3,3)=0.d0
        shls(4,3)=onMIz8
        shls(5,3)=-(1.d0+x)/8.d0
        shls(6,3)=-onPLy8
        shls(1,4)=0.d0
        shls(2,4)=0.d0
        shls(3,4)=0.d0
        shls(4,4)=-onMIz8
        shls(5,4)=-(1.d0-x)/8.d0
        shls(6,4)=onPLy8
        shls(1,5)=0.d0
        shls(2,5)=0.d0
        shls(3,5)=0.d0
        shls(4,5)=onPLz8
        shls(5,5)=-(1.d0-x)/8.d0
        shls(6,5)=-onMIy8
        shls(1,6)=0.d0
        shls(2,6)=0.d0
        shls(3,6)=0.d0
        shls(4,6)=-onPLz8
        shls(5,6)=-(1.d0+x)/8.d0
        shls(6,6)=onMIy8
        shls(1,7)=0.d0
        shls(2,7)=0.d0
        shls(3,7)=0.d0
        shls(4,7)=onPLz8
        shls(5,7)=(1.d0+x)/8.d0
        shls(6,7)=onPLy8
        shls(1,8)=0.d0
        shls(2,8)=0.d0
        shls(3,8)=0.d0
        shls(4,8)=-onPLz8
        shls(5,8)=(1.d0-x)/8.d0
        shls(6,8)=-onPLy8
        end if

      endif
c
c for 27 noded 3 dimensional brick element
c
      if(nel.eq.27) then

        shld(1,1)=(2.d0*x-1.d0)*y*(y-1.d0)*z*(z-1.d0)/8.d0
        shld(2,1)=x*(x-1.d0)*(2.d0*y-1.d0)*z*(z-1.d0)/8.d0
        shld(3,1)=x*(x-1.d0)*y*(y-1.d0)*(2.d0*z-1.d0)/8.d0
        shl(1,1)=x*(x-1.d0)*y*(y-1.d0)*z*(z-1.d0)/8.d0

        shld(1,2)=(2.d0*x+1.d0)*y*(y-1.d0)*z*(z-1.d0)/8.d0
        shld(2,2)=x*(x+1.d0)*(2.d0*y-1.d0)*z*(z-1.d0)/8.d0
        shld(3,2)=x*(x+1.d0)*y*(y-1.d0)*(2.d0*z-1.d0)/8.d0
        shl(1,2)=x*(x+1.d0)*y*(y-1.d0)*z*(z-1.d0)/8.d0

        shld(1,3)=(2.d0*x+1.d0)*y*(y+1.d0)*z*(z-1.d0)/8.d0
        shld(2,3)=x*(x+1.d0)*(2.d0*y+1.d0)*z*(z-1.d0)/8.d0
        shld(3,3)=x*(x+1.d0)*y*(y+1.d0)*(2.d0*z-1.d0)/8.d0
        shl(1,3)=x*(x+1.d0)*y*(y+1.d0)*z*(z-1.d0)/8.d0

        shld(1,4)=(2.d0*x-1.d0)*y*(y+1.d0)*z*(z-1.d0)/8.d0
        shld(2,4)=x*(x-1.d0)*(2.d0*y+1.d0)*z*(z-1.d0)/8.d0
        shld(3,4)=x*(x-1.d0)*y*(y+1.d0)*(2.d0*z-1.d0)/8.d0
        shl(1,4)=x*(x-1.d0)*y*(y+1.d0)*z*(z-1.d0)/8.d0

        shld(1,5)=(2.d0*x-1.d0)*y*(y-1.d0)*z*(z+1.d0)/8.d0
        shld(2,5)=x*(x-1.d0)*(2.d0*y-1.d0)*z*(z+1.d0)/8.d0
        shld(3,5)=x*(x-1.d0)*y*(y-1.d0)*(2.d0*z+1.d0)/8.d0
        shl(1,5)=x*(x-1.d0)*y*(y-1.d0)*z*(z+1.d0)/8.d0

        shld(1,6)=(2.d0*x+1.d0)*y*(y-1.d0)*z*(z+1.d0)/8.d0
        shld(2,6)=x*(x+1.d0)*(2.d0*y-1.d0)*z*(z+1.d0)/8.d0
        shld(3,6)=x*(x+1.d0)*y*(y-1.d0)*(2.d0*z+1.d0)/8.d0
        shl(1,6)=x*(x+1.d0)*y*(y-1.d0)*z*(z+1.d0)/8.d0

        shld(1,7)=(2.d0*x+1.d0)*y*(y+1.d0)*z*(z+1.d0)/8.d0
        shld(2,7)=x*(x+1.d0)*(2.d0*y+1.d0)*z*(z+1.d0)/8.d0
        shld(3,7)=x*(x+1.d0)*y*(y+1.d0)*(2.d0*z+1.d0)/8.d0
        shl(1,7)=x*(x+1.d0)*y*(y+1.d0)*z*(z+1.d0)/8.d0

        shld(1,8)=(2.d0*x-1.d0)*y*(y+1.d0)*z*(z+1.d0)/8.d0
        shld(2,8)=x*(x-1.d0)*(2.d0*y+1.d0)*z*(z+1.d0)/8.d0
        shld(3,8)=x*(x-1.d0)*y*(y+1.d0)*(2.d0*z+1.d0)/8.d0
        shl(1,8)=x*(x-1.d0)*y*(y+1.d0)*z*(z+1.d0)/8.d0

        shld(1,17)=(2.d0*x-1.d0)*y*(y-1.d0)*(1.d0-z*z)/4.d0
        shld(2,17)=x*(x-1.d0)*(2.d0*y-1.d0)*(1.d0-z*z)/4.d0
        shld(3,17)=x*(x-1.d0)*y*(y-1.d0)*(-2.d0*z)/4.d0
        shl(1,17)=x*(x-1.d0)*y*(y-1.d0)*(1.d0-z*z)/4.d0

        shld(1,18)=(2.d0*x+1.d0)*y*(y-1.d0)*(1.d0-z*z)/4.d0
        shld(2,18)=x*(x+1.d0)*(2.d0*y-1.d0)*(1.d0-z*z)/4.d0
        shld(3,18)=x*(x+1.d0)*y*(y-1.d0)*(-2.d0*z)/4.d0
        shl(1,18)=x*(x+1.d0)*y*(y-1.d0)*(1.d0-z*z)/4.d0

        shld(1,19)=(2.d0*x+1.d0)*y*(y+1.d0)*(1.d0-z*z)/4.d0
        shld(2,19)=x*(x+1.d0)*(2.d0*y+1.d0)*(1.d0-z*z)/4.d0
        shld(3,19)=x*(x+1.d0)*y*(y+1.d0)*(-2.d0*z)/4.d0
        shl(1,19)=x*(x+1.d0)*y*(y+1.d0)*(1.d0-z*z)/4.d0

        shld(1,20)=(2.d0*x-1.d0)*y*(y+1.d0)*(1.d0-z*z)/4.d0
        shld(2,20)=x*(x-1.d0)*(2.d0*y+1.d0)*(1.d0-z*z)/4.d0
        shld(3,20)=x*(x-1.d0)*y*(y+1.d0)*(-2.d0*z)/4.d0
        shl(1,20)=x*(x-1.d0)*y*(y+1.d0)*(1.d0-z*z)/4.d0

        shld(1,9)=(-2.d0*x)*y*(y-1.d0)*z*(z-1.d0)/4.d0
        shld(2,9)=(1.d0-x*x)*(2.d0*y-1.d0)*z*(z-1.d0)/4.d0
        shld(3,9)=(1.d0-x*x)*y*(y-1.d0)*(2.d0*z-1.d0)/4.d0
        shl(1,9)=(1.d0-x*x)*y*(y-1.d0)*z*(z-1.d0)/4.d0

        shld(1,10)=(2.d0*x+1.d0)*(1.d0-y*y)*z*(z-1.d0)/4.d0
        shld(2,10)=x*(x+1.d0)*(-2.d0*y)*z*(z-1.d0)/4.d0
        shld(3,10)=x*(x+1.d0)*(1.d0-y*y)*(2.d0*z-1.d0)/4.d0
        shl(1,10)=x*(x+1.d0)*(1.d0-y*y)*z*(z-1.d0)/4.d0

        shld(1,11)=(-2.d0*x)*y*(y+1.d0)*z*(z-1.d0)/4.d0
        shld(2,11)=(1.d0-x*x)*(2.d0*y+1.d0)*z*(z-1.d0)/4.d0
        shld(3,11)=(1.d0-x*x)*y*(y+1.d0)*(2.d0*z-1.d0)/4.d0
        shl(1,11)=(1.d0-x*x)*y*(y+1.d0)*z*(z-1.d0)/4.d0

        shld(1,12)=(2.d0*x-1.d0)*(1.d0-y*y)*z*(z-1.d0)/4.d0
        shld(2,12)=x*(x-1.d0)*(-2.d0*y)*z*(z-1.d0)/4.d0
        shld(3,12)=x*(x-1.d0)*(1.d0-y*y)*(2.d0*z-1.d0)/4.d0
        shl(1,12)=x*(x-1.d0)*(1.d0-y*y)*z*(z-1.d0)/4.d0

        shld(1,21)=(-2.d0*x)*(1.d0-y*y)*z*(z-1.d0)/2.d0
        shld(2,21)=(1.d0-x*x)*(-2.d0*y)*z*(z-1.d0)/2.d0
        shld(3,21)=(1.d0-x*x)*(1.d0-y*y)*(2.d0*z-1.d0)/2.d0
        shl(1,21)=(1.d0-x*x)*(1.d0-y*y)*z*(z-1.d0)/2.d0

        shld(1,13)=(-2.d0*x)*y*(y-1.d0)*z*(z+1.d0)/4.d0
        shld(2,13)=(1.d0-x*x)*(2.d0*y-1.d0)*z*(z+1.d0)/4.d0
        shld(3,13)=(1.d0-x*x)*y*(y-1.d0)*(2.d0*z+1.d0)/4.d0
        shl(1,13)=(1.d0-x*x)*y*(y-1.d0)*z*(z+1.d0)/4.d0

        shld(1,14)=(2.d0*x+1.d0)*(1.d0-y*y)*z*(z+1.d0)/4.d0
        shld(2,14)=x*(x+1.d0)*(-2.d0*y)*z*(z+1.d0)/4.d0
        shld(3,14)=x*(x+1.d0)*(1.d0-y*y)*(2.d0*z+1.d0)/4.d0
        shl(1,14)=x*(x+1.d0)*(1.d0-y*y)*z*(z+1.d0)/4.d0

        shld(1,15)=(-2.d0*x)*y*(y+1.d0)*z*(z+1.d0)/4.d0
        shld(2,15)=(1.d0-x*x)*(2.d0*y+1.d0)*z*(z+1.d0)/4.d0
        shld(3,15)=(1.d0-x*x)*y*(y+1.d0)*(2.d0*z+1.d0)/4.d0
        shl(1,15)=(1.d0-x*x)*y*(y+1.d0)*z*(z+1.d0)/4.d0

        shld(1,16)=(2.d0*x-1.d0)*(1.d0-y*y)*z*(z+1.d0)/4.d0
        shld(2,16)=x*(x-1.d0)*(-2.0*y)*z*(z+1.d0)/4.d0
        shld(3,16)=x*(x-1.d0)*(1.d0-y*y)*(2.d0*z+1.d0)/4.d0
        shl(1,16)=x*(x-1.d0)*(1.d0-y*y)*z*(z+1.d0)/4.d0

        shld(1,22)=(-2.d0*x)*(1.d0-y*y)*z*(z+1.d0)/2.d0
        shld(2,22)=(1.d0-x*x)*(-2.d0*y)*z*(z+1.d0)/2.d0
        shld(3,22)=(1.d0-x*x)*(1.d0-y*y)*(2.d0*z+1.d0)/2.d0
        shl(1,22)=(1.d0-x*x)*(1.d0-y*y)*z*(z+1.d0)/2.d0

        shld(1,25)=(-2.d0*x)*y*(y-1.d0)*(1.d0-z*z)/2.d0
        shld(2,25)=(1.d0-x*x)*(2.d0*y-1.d0)*(1.d0-z*z)/2.d0
        shld(3,25)=(1.d0-x*x)*y*(y-1.d0)*(-2.d0*z)/2.d0
        shl(1,25)=(1.d0-x*x)*y*(y-1.d0)*(1.d0-z*z)/2.d0

        shld(1,24)=(2.d0*x+1.d0)*(1.d0-y*y)*(1.d0-z*z)/2.d0
        shld(2,24)=x*(x+1.d0)*(-2.d0*y)*(1.d0-z*z)/2.d0
        shld(3,24)=x*(x+1.d0)*(1.d0-y*y)*(-2.d0*z)/2.d0
        shl(1,24)=x*(x+1.d0)*(1.d0-y*y)*(1.d0-z*z)/2.d0

        shld(1,26)=(-2.d0*x)*y*(y+1.d0)*(1.d0-z*z)/2.d0
        shld(2,26)=(1.d0-x*x)*(2.d0*y+1.d0)*(1.d0-z*z)/2.d0
        shld(3,26)=(1.d0-x*x)*y*(y+1.d0)*(-2.d0*z)/2.d0
        shl(1,26)=(1.d0-x*x)*y*(y+1.d0)*(1.d0-z*z)/2.d0

        shld(1,23)=(2.d0*x-1.d0)*(1.d0-y*y)*(1.d0-z*z)/2.d0
        shld(2,23)=x*(x-1.d0)*(-2.d0*y)*(1.d0-z*z)/2.d0
        shld(3,23)=x*(x-1.d0)*(1.d0-y*y)*(-2.d0*z)/2.d0
        shl(1,23)=x*(x-1.d0)*(1.d0-y*y)*(1.d0-z*z)/2.d0

        shld(1,27)=(-2.d0*x)*(1.d0-y*y)*(1.d0-z*z)
        shld(2,27)=(1.d0-x*x)*(-2.d0*y)*(1.d0-z*z)
        shld(3,27)=(1.d0-x*x)*(1.d0-y*y)*(-2.d0*z)
        shl(1,27)=(1.d0-x*x)*(1.d0-y*y)*(1.d0-z*z)

c       Bubble node
        if(bf) then
          bubble(1)=(-4.d0*x*x*x)*(1.d0-y*y*y*y)*(1.d0-z*z*z*z)
          bubble(2)=(1.d0-x*x*x*x)*(-4.d0*y*y*y)*(1.d0-z*z*z*z)
          bubble(3)=(1.d0-x*x*x*x)*(1.d0-y*y*y*y)*(-4.d0*z*z*z)
          bubble(4)=(1.d0-x*x*x*x)*(1.d0-y*y*y*y)*(1.d0-z*z*z*z)
        endif

c     calculate second derivatives for shape funtions in local coordinate

      if(der) then

        shls(1,1)=2.d0*y*(y-1.d0)*z*(z-1.d0)/8.d0
        shls(2,1)=x*(x-1.d0)*2.d0*z*(z-1.d0)/8.d0
        shls(3,1)=x*(x-1.d0)*y*(y-1.d0)*2.d0/8.d0
        shls(4,1)=(2.d0*x-1.d0)*(2.d0*y-1.d0)*z*(z-1.d0)/8.d0
        shls(5,1)=x*(x-1.d0)*(2.d0*y-1.d0)*(2.d0*z-1.d0)/8.d0
        shls(6,1)=(2.d0*x-1.d0)*y*(y-1.d0)*(2.d0*z-1.d0)/8.d0

        shls(1,2)=2.d0*y*(y-1.d0)*z*(z-1.d0)/8.d0
        shls(2,2)=x*(x+1.d0)*2.d0*z*(z-1.d0)/8.d0
        shls(3,2)=x*(x+1.d0)*y*(y-1.d0)*2.d0/8.d0
        shls(4,2)=(2.d0*x+1.d0)*(2.d0*y-1.d0)*z*(z-1.d0)/8.d0
        shls(5,2)=x*(x+1.d0)*(2.d0*y-1.d0)*(2.d0*z-1.d0)/8.d0
        shls(6,2)=(2.d0*x+1.d0)*y*(y-1.d0)*(2.d0*z-1.d0)/8.d0

        shls(1,3)=2.d0*y*(y+1.d0)*z*(z-1.d0)/8.d0
        shls(2,3)=x*(x+1.d0)*2.d0*z*(z-1.d0)/8.d0
        shls(3,3)=x*(x+1.d0)*y*(y+1.d0)*2.d0/8.d0
        shls(4,3)=(2.d0*x+1.d0)*(2.d0*y+1.d0)*z*(z-1.d0)/8.d0
        shls(5,3)=x*(x+1.d0)*(2.d0*y+1.d0)*(2.d0*z-1.d0)/8.d0
        shls(6,3)=(2.d0*x+1.d0)*y*(y+1.d0)*(2.d0*z-1.d0)/8.d0

        shls(1,4)=2.d0*y*(y+1.d0)*z*(z-1.d0)/8.d0
        shls(2,4)=x*(x-1.d0)*2.d0*z*(z-1.d0)/8.d0
        shls(3,4)=x*(x-1.d0)*y*(y+1.d0)*2.d0/8.d0
        shls(4,4)=(2.d0*x-1.d0)*(2.d0*y+1.d0)*z*(z-1.d0)/8.d0
        shls(5,4)=x*(x-1.d0)*(2.d0*y+1.d0)*(2.d0*z-1.d0)/8.d0
        shls(6,4)=(2.d0*x-1.d0)*y*(y+1.d0)*(2.d0*z-1.d0)/8.d0

        shls(1,5)=2.d0*y*(y-1.d0)*z*(z+1.d0)/8.d0
        shls(2,5)=x*(x-1.d0)*2.d0*z*(z+1.d0)/8.d0
        shls(3,5)=x*(x-1.d0)*y*(y-1.d0)*2.d0/8.d0
        shls(4,5)=(2.d0*x-1.d0)*(2.d0*y-1.d0)*z*(z+1.d0)/8.d0
        shls(5,5)=x*(x-1.d0)*(2.d0*y-1.d0)*(2.d0*z+1.d0)/8.d0
        shls(6,5)=(2.d0*x-1.d0)*y*(y-1.d0)*(2.d0*z+1.d0)/8.d0

        shls(1,6)=2.d0*y*(y-1.d0)*z*(z+1.d0)/8.d0
        shls(2,6)=x*(x+1.d0)*2.d0*z*(z+1.d0)/8.d0
        shls(3,6)=x*(x+1.d0)*y*(y-1.d0)*2.d0/8.d0
        shls(4,6)=(2.d0*x+1.d0)*(2.d0*y-1.d0)*z*(z+1.d0)/8.d0
        shls(5,6)=x*(x+1.d0)*(2.d0*y-1.d0)*(2.d0*z+1.d0)/8.d0
        shls(6,6)=(2.d0*x+1.d0)*y*(y-1.d0)*(2.d0*z+1.d0)/8.d0

        shls(1,7)=2.d0*y*(y+1.d0)*z*(z+1.d0)/8.d0
        shls(2,7)=x*(x+1.d0)*2.d0*z*(z+1.d0)/8.d0
        shls(3,7)=x*(x+1.d0)*y*(y+1.d0)*2.d0/8.d0
        shls(4,7)=(2.d0*x+1.d0)*(2.d0*y+1.d0)*z*(z+1.d0)/8.d0
        shls(5,7)=x*(x+1.d0)*(2.d0*y+1.d0)*(2.d0*z+1.d0)/8.d0
        shls(6,7)=(2.d0*x+1.d0)*y*(y+1.d0)*(2.d0*z+1.d0)/8.d0

        shls(1,8)=2.d0*y*(y+1.d0)*z*(z+1.d0)/8.d0
        shls(2,8)=x*(x-1.d0)*2.d0*z*(z+1.d0)/8.d0
        shls(3,8)=x*(x-1.d0)*y*(y+1.d0)*2.d0/8.d0
        shls(4,8)=(2.d0*x-1.d0)*(2.d0*y+1.d0)*z*(z+1.d0)/8.d0
        shls(5,8)=x*(x-1.d0)*(2.d0*y+1.d0)*(2.d0*z+1.d0)/8.d0
        shls(6,8)=(2.d0*x-1.d0)*y*(y+1.d0)*(2.d0*z+1.d0)/8.d0

        shls(1,17)=2.d0*y*(y-1.d0)*(1.d0-z*z)/4.d0
        shls(2,17)=x*(x-1.d0)*2.d0*(1.d0-z*z)/4.d0
        shls(3,17)=x*(x-1.d0)*y*(y-1.d0)*(-2.d0)/4.d0
        shls(4,17)=(2.d0*x-1.d0)*(2.d0*y-1.d0)*(1.d0-z*z)/4.d0
        shls(5,17)=x*(x-1.d0)*(2.d0*y-1.d0)*(-2.d0*z)/4.d0
        shls(6,17)=(2.d0*x-1.d0)*y*(y-1.d0)*(-2.d0*z)/4.d0

        shls(1,18)=2.d0*y*(y-1.d0)*(1.d0-z*z)/4.d0
        shls(2,18)=x*(x+1.d0)*2.d0*(1.d0-z*z)/4.d0
        shls(3,18)=x*(x+1.d0)*y*(y-1.d0)*(-2.d0)/4.d0
        shls(4,18)=(2.d0*x+1.d0)*(2.d0*y-1.d0)*(1.d0-z*z)/4.d0
        shls(5,18)=x*(x+1.d0)*(2.d0*y-1.d0)*(-2.d0*z)/4.d0
        shls(6,18)=(2.d0*x+1.d0)*y*(y-1.d0)*(-2.d0*z)/4.d0

        shls(1,19)=2.d0*y*(y+1.d0)*(1.d0-z*z)/4.d0
        shls(2,19)=x*(x+1.d0)*2.d0*(1.d0-z*z)/4.d0
        shls(3,19)=x*(x+1.d0)*y*(y+1.d0)*(-2.d0)/4.d0
        shls(4,19)=(2.d0*x+1.d0)*(2.d0*y+1.d0)*(1.d0-z*z)/4.d0
        shls(5,19)=x*(x+1.d0)*(2.d0*y+1.d0)*(-2.d0*z)/4.d0
        shls(6,19)=(2.d0*x+1.d0)*y*(y+1.d0)*(-2.d0*z)/4.d0

        shls(1,20)=2.d0*y*(y+1.d0)*(1.d0-z*z)/4.d0
        shls(2,20)=x*(x-1.d0)*2.d0*(1.d0-z*z)/4.d0
        shls(3,20)=x*(x-1.d0)*y*(y+1.d0)*(-2.d0)/4.d0
        shls(4,20)=(2.d0*x-1.d0)*(2.d0*y+1.d0)*(1.d0-z*z)/4.d0
        shls(5,20)=x*(x-1.d0)*(2.d0*y+1.d0)*(-2.d0*z)/4.d0
        shls(6,20)=(2.d0*x-1.d0)*y*(y+1.d0)*(-2.d0*z)/4.d0

        shls(1,9)=(-2.d0)*y*(y-1.d0)*z*(z-1.d0)/4.d0
        shls(2,9)=(1.d0-x*x)*2.d0*z*(z-1.d0)/4.d0
        shls(3,9)=(1.d0-x*x)*y*(y-1.d0)*2.d0/4.d0
        shls(4,9)=(-2.d0*x)*(2.d0*y-1.d0)*z*(z-1.d0)/4.d0
        shls(5,9)=(1.d0-x*x)*(2.d0*y-1.d0)*(2.d0*z-1.d0)/4.d0
        shls(6,9)=(-2.d0*x)*y*(y-1.d0)*(2.d0*z-1.d0)/4.d0

        shls(1,10)=2.d0*(1.d0-y*y)*z*(z-1.d0)/4.d0
        shls(2,10)=x*(x+1.d0)*(-2.d0)*z*(z-1.d0)/4.d0
        shls(3,10)=x*(x+1.d0)*(1.d0-y*y)*2.d0/4.d0
        shls(4,10)=(2.d0*x+1.d0)*(-2.d0*y)*z*(z-1.d0)/4.d0
        shls(5,10)=x*(x+1.d0)*(-2.d0*y)*(2.d0*z-1.d0)/4.d0
        shls(6,10)=(2.D0*x+1.d0)*(1.d0-y*y)*(2.d0*z-1.d0)/4.d0

        shls(1,11)=(-2.d0)*y*(y+1.d0)*z*(z-1.d0)/4.d0
        shls(2,11)=(1.d0-x*x)*2.d0*z*(z-1.d0)/4.d0
        shls(3,11)=(1.d0-x*x)*y*(y+1.d0)*2.d0/4.d0
        shls(4,11)=(-2.d0*x)*(2.D0*y+1.d0)*z*(z-1.d0)/4.d0
        shls(5,11)=(1.d0-x*x)*(2.d0*y+1.d0)*(2.D0*z-1.d0)/4.d0
        shls(6,11)=(-2.D0*x)*y*(y+1.d0)*(2.d0*z-1.d0)/4.d0

        shls(1,12)=2.d0*(1.d0-y*y)*z*(z-1.d0)/4.d0
        shls(2,12)=x*(x-1.d0)*(-2.d0)*z*(z-1.d0)/4.d0
        shls(3,12)=x*(x-1.d0)*(1.d0-y*y)*2.d0/4.d0
        shls(4,12)=(2.d0*x-1.d0)*(-2.D0*y)*z*(z-1.d0)/4.d0
        shls(5,12)=x*(x-1.d0)*(-2.d0*y)*(2.D0*z-1.d0)/4.d0
        shls(6,12)=(2.D0*x-1.d0)*(1.d0-y*y)*(2.d0*z-1.d0)/4.d0

        shls(1,21)=(-2.d0)*(1.d0-y*y)*z*(z-1.d0)/2.d0
        shls(2,21)=(1.d0-x*x)*(-2.d0)*z*(z-1.d0)/2.d0
        shls(3,21)=(1.d0-x*x)*(1.d0-y*y)*2.d0/2.d0
        shls(4,21)=(-2.d0*x)*(-2.D0*y)*z*(z-1.d0)/2.d0
        shls(5,21)=(1.d0-x*x)*(-2.d0*y)*(2.D0*z-1.d0)/2.d0
        shls(6,21)=(-2.D0*x)*(1.d0-y*y)*(2.d0*z-1.d0)/2.d0

        shls(1,13)=(-2.d0)*y*(y-1.d0)*z*(z+1.d0)/4.d0
        shls(2,13)=(1.d0-x*x)*2.d0*z*(z+1.d0)/4.d0
        shls(3,13)=(1.d0-x*x)*y*(y-1.d0)*2.d0/4.d0
        shls(4,13)=(-2.d0*x)*(2.D0*y-1.d0)*z*(z+1.d0)/4.d0
        shls(5,13)=(1.d0-x*x)*(2.d0*y-1.d0)*(2.D0*z+1.d0)/4.d0
        shls(6,13)=(-2.D0*x)*y*(y-1.d0)*(2.d0*z+1.d0)/4.d0

        shls(1,14)=2.d0*(1.d0-y*y)*z*(z+1.d0)/4.d0
        shls(2,14)=x*(x+1.d0)*(-2.d0)*z*(z+1.d0)/4.d0
        shls(3,14)=x*(x+1.d0)*(1.d0-y*y)*2.d0/4.d0
        shls(4,14)=(2.d0*x+1.d0)*(-2.D0*y)*z*(z+1.d0)/4.d0
        shls(5,14)=x*(x+1.d0)*(-2.d0*y)*(2.D0*z+1.d0)/4.d0
        shls(6,14)=(2.D0*x+1.d0)*(1.d0-y*y)*(2.d0*z+1.d0)/4.d0

        shls(1,15)=(-2.d0)*y*(y+1.d0)*z*(z+1.d0)/4.d0
        shls(2,15)=(1.d0-x*x)*2.d0*z*(z+1.d0)/4.d0
        shls(3,15)=(1.d0-x*x)*y*(y+1.d0)*2.d0/4.d0
        shls(4,15)=(-2.d0*x)*(2.D0*y+1.d0)*z*(z+1.d0)/4.d0
        shls(5,15)=(1.d0-x*x)*(2.d0*y+1.d0)*(2.D0*z+1.d0)/4.d0
        shls(6,15)=(-2.D0*x)*y*(y+1.d0)*(2.d0*z+1.d0)/4.d0

        shls(1,16)=2.d0*(1.d0-y*y)*z*(z+1.d0)/4.d0
        shls(2,16)=x*(x-1.d0)*(-2.0)*z*(z+1.d0)/4.d0
        shls(3,16)=x*(x-1.d0)*(1.d0-y*y)*2.d0/4.d0
        shls(4,16)=(2.d0*x-1.d0)*(-2.d0*y)*z*(z+1.d0)/4.d0
        shls(5,16)=x*(x-1.d0)*(-2.0*y)*(2.d0*z+1.d0)/4.d0
        shls(6,16)=(2.d0*x-1.d0)*(1.d0-y*y)*(2.d0*z+1.d0)/4.d0

        shls(1,22)=(-2.d0)*(1.d0-y*y)*z*(z+1.d0)/2.d0
        shls(2,22)=(1.d0-x*x)*(-2.d0)*z*(z+1.d0)/2.d0
        shls(3,22)=(1.d0-x*x)*(1.d0-y*y)*2.d0/2.d0
        shls(4,22)=(-2.d0*x)*(-2.d0*y)*z*(z+1.d0)/2.d0
        shls(5,22)=(1.d0-x*x)*(-2.d0*y)*(2.d0*z+1.d0)/2.d0
        shls(6,22)=(-2.d0*x)*(1.d0-y*y)*(2.d0*z+1.d0)/2.d0

        shls(1,25)=(-2.d0)*y*(y-1.d0)*(1.d0-z*z)/2.d0
        shls(2,25)=(1.d0-x*x)*2.d0*(1.d0-z*z)/2.d0
        shls(3,25)=(1.d0-x*x)*y*(y-1.d0)*(-2.d0)/2.d0
        shls(4,25)=(-2.d0*x)*(2.d0*y-1.d0)*(1.d0-z*z)/2.d0
        shls(5,25)=(1.d0-x*x)*(2.d0*y-1.d0)*(-2.d0*z)/2.d0
        shls(6,25)=(-2.d0*x)*y*(y-1.d0)*(-2.d0*z)/2.d0

        shls(1,24)=2.d0*(1.d0-y*y)*(1.d0-z*z)/2.d0
        shls(2,24)=x*(x+1.d0)*(-2.d0)*(1.d0-z*z)/2.d0
        shls(3,24)=x*(x+1.d0)*(1.d0-y*y)*(-2.d0)/2.d0
        shls(4,24)=(2.d0*x+1.d0)*(-2.d0*y)*(1.d0-z*z)/2.d0
        shls(5,24)=x*(x+1.d0)*(-2.d0*y)*(-2.d0*z)/2.d0
        shls(6,24)=(2.d0*x+1.d0)*(1.d0-y*y)*(-2.d0*z)/2.d0

        shls(1,26)=(-2.d0)*y*(y+1.d0)*(1.d0-z*z)/2.d0
        shls(2,26)=(1.d0-x*x)*2.d0*(1.d0-z*z)/2.d0
        shls(3,26)=(1.d0-x*x)*y*(y+1.d0)*(-2.d0)/2.d0
        shls(4,26)=(-2.d0*x)*(2.d0*y+1.d0)*(1.d0-z*z)/2.d0
        shls(5,26)=(1.d0-x*x)*(2.d0*y+1.d0)*(-2.d0*z)/2.d0
        shls(6,26)=(-2.d0*x)*y*(y+1.d0)*(-2.d0*z)/2.d0

        shls(1,23)=2.d0*(1.d0-y*y)*(1.d0-z*z)/2.d0
        shls(2,23)=x*(x-1.d0)*(-2.d0)*(1.d0-z*z)/2.d0
        shls(3,23)=x*(x-1.d0)*(1.d0-y*y)*(-2.d0)/2.d0
        shls(4,23)=(2.d0*x-1.d0)*(-2.d0*y)*(1.d0-z*z)/2.d0
        shls(5,23)=x*(x-1.d0)*(-2.d0*y)*(-2.d0*z)/2.d0
        shls(6,23)=(2.d0*x-1.d0)*(1.d0-y*y)*(-2.d0*z)/2.d0

        shls(1,27)=(-2.d0)*(1.d0-y*y)*(1.d0-z*z)
        shls(2,27)=(1.d0-x*x)*(-2.d0)*(1.d0-z*z)
        shls(3,27)=(1.d0-x*x)*(1.d0-y*y)*(-2.d0)
        shls(4,27)=(-2.d0*x)*(-2.d0*y)*(1.d0-z*z)
        shls(5,27)=(1.d0-x*x)*(-2.d0*y)*(-2.d0*z)
        shls(6,27)=(-2.d0*x)*(1.d0-y*y)*(-2.d0*z)

      endif

      endif

      return
      end

c----------------------------------------------------------------------
c
      subroutine shgb(xl,nel,shp,shp2,neg,nen,bf,der,xsj,cartd1,
     &                cartd2,bubble,sx)
c
c      8 noded, 27 noded brick element + bubble point
c
c      Written by Kaiming Xia, 09/20/2001
c         Revised by Ramon Calderer, 03/07/2008
c         Revised by Timothy Truster, 02/20/2012
c
c----------------------------------------------------------------------
c
c      Purpose: Compute 3-d isoparametric shape
c               functions and their derivatives w/r x,y,z

c      Inputs:
c         ss(3)     - Natural coordinates of point
c         xl(ndm,*) - Nodal coordinates for element
c         ndm       - Spatial dimension of mesh
c             ideriv          - Flag for activation of second derivatives (default=0
c                     don't compute second derivative
c      Outputs:
c         xsj       - Jacobian determinant at point
c         shp(4,*)  - Shape functions and derivatives at local point
c         shp2(6,*) - Shape functions and derivatives at local point
c                     cartd(1,i) = dN_i/dx
c                     cartd(2,i) = dN_i/dy
c                     cartd(3,i) = dN_i/dz
c                     cartd(4,i) =  N_i
c                     cartd2(1,i) = d^2N_i/dx^2
c                     cartd2(2,i) = d^2N_i/dy^2
c                     cartd2(3,i) = d^2N_i/dz^2
c                     cartd2(4,i) = d^2N_i/dxy
c                     cartd2(5,i) = d^2N_i/dyz
c                     cartd2(6,i) = d^2N_i/dzx
c
c----------------------------------------------------------------------

      implicit none

!     Input Variables
      integer nel,nen,neg
      logical der,bf
      real*8 xl(3,nel),shp(3,nen),shp2(6,nen),bubble(4)

!     Output Variables
      real*8 xsj, cartd1(3,nen),cartd2(6,nen),sx(3,3)

!     Local Variables
      integer   i , j , k, inode, n
      real*8    rxsj,c1(6,3),tc(6,3),tcj(6,3),at1(6),at2(6)
      real*8    t2(6,6),temp(3)
      real*8    xs(3,3),ad(3,3),u

      integer iin,iout,irsin,irsout
      common /iounit/ iin,iout,irsin,irsout

c     Compute jacobian transformation
c        Calculate c1 matrix simultaneously

        do 30 i=1,3
        do 30 j=1,3
        xs(i,j)=0.d0
        c1(i,j)=0.d0
        do 30 inode=1,nel
        xs(i,j)=xs(i,j)+shp(i,inode)*xl(j,inode)
        c1(i,j)=c1(i,j)+shp2(i,inode)*xl(j,inode)
30        continue
        do 40 i=4,6
        do 40 j=1,3
        c1(i,j)=0.d0
        do 40 inode=1,nel
        c1(i,j)=c1(i,j)+shp2(i,inode)*xl(j,inode)
40        continue

! Output transpose which agrees with customary array definition
      do i = 1,3
        do j = 1,3
          sx(i,j) = xs(j,i)
        enddo
      enddo

c below is to calculate determinant and inverswe of jacobian matrix
c
c     Compute adjoint to jacobian

      ad(1,1) = xs(2,2)*xs(3,3) - xs(2,3)*xs(3,2)
      ad(1,2) = xs(3,2)*xs(1,3) - xs(3,3)*xs(1,2)
      ad(1,3) = xs(1,2)*xs(2,3) - xs(1,3)*xs(2,2)
      ad(2,1) = xs(2,3)*xs(3,1) - xs(2,1)*xs(3,3)
      ad(2,2) = xs(3,3)*xs(1,1) - xs(3,1)*xs(1,3)
      ad(2,3) = xs(1,3)*xs(2,1) - xs(1,1)*xs(2,3)
      ad(3,1) = xs(2,1)*xs(3,2) - xs(2,2)*xs(3,1)
      ad(3,2) = xs(3,1)*xs(1,2) - xs(3,2)*xs(1,1)
      ad(3,3) = xs(1,1)*xs(2,2) - xs(1,2)*xs(2,1)

c     Compute determinant of jacobian

      xsj  = xs(1,1)*ad(1,1) + xs(1,2)*ad(2,1) + xs(1,3)*ad(3,1)
        if (xsj.le.(0.d0)) then
c         write(iout,1000) nen,neg
          write(*,1000) nen,neg
          stop
        endif

      rxsj = 1.d0/xsj
c     Compute jacobian inverse
      do j = 1,3
        do i = 1,3
          xs(i,j) = ad(i,j)*rxsj
        end do
      end do

        do 55 n=1,nel
        do 55 i=1,3
        cartd1(i,n)=0.d0
        do 55 k=1,3
55        cartd1(i,n)=cartd1(i,n)+xs(i,k)*shp(k,n)


c     Initialize second derivatives with respect to global coords.
        do i=1,6
          do n=1,nel
            cartd2(i,n)=0.0d0
          end do
        end do

      if (bf) then

          temp(1) = bubble(1)
          temp(2) = bubble(2)
          temp(3) = bubble(3)
          do i = 1,3
            bubble(i) = temp(1)*xs(i,1) + temp(2)*xs(i,2)
     &                + temp(3)*xs(i,3)
          enddo

      endif


        IF(der) then
C---------------------------------------------------------------
C        Below is to calculate second derivative of shape functions

C        Calculate T2 matrix
        t2(1,1)=xs(1,1)**2
        t2(1,2)=xs(1,2)**2
        t2(1,3)=xs(1,3)**2
        t2(1,4)=2.d0*xs(1,1)*xs(1,2)
        t2(1,5)=2.d0*xs(1,2)*xs(1,3)
        t2(1,6)=2.d0*xs(1,3)*xs(1,1)

        t2(2,1)=xs(2,1)**2
        t2(2,2)=xs(2,2)**2
        t2(2,3)=xs(2,3)**2
        t2(2,4)=2.d0*xs(2,1)*xs(2,2)
        t2(2,5)=2.d0*xs(2,2)*xs(2,3)
        t2(2,6)=2.d0*xs(2,3)*xs(2,1)

        t2(3,1)=xs(3,1)**2
        t2(3,2)=xs(3,2)**2
        t2(3,3)=xs(3,3)**2
        t2(3,4)=2.d0*xs(3,1)*xs(3,2)
        t2(3,5)=2.d0*xs(3,2)*xs(3,3)
        t2(3,6)=2.d0*xs(3,3)*xs(3,1)

        t2(4,1)=xs(1,1)*xs(2,1)
        t2(4,2)=xs(1,2)*xs(2,2)
        t2(4,3)=xs(1,3)*xs(2,3)
        t2(4,4)=xs(1,1)*xs(2,2)+xs(1,2)*xs(2,1)
        t2(4,5)=xs(1,2)*xs(2,3)+xs(1,3)*xs(2,2)
        t2(4,6)=xs(1,1)*xs(2,3)+xs(1,3)*xs(2,1)

        t2(5,1)=xs(2,1)*xs(3,1)
        t2(5,2)=xs(2,2)*xs(3,2)
        t2(5,3)=xs(2,3)*xs(3,3)
        t2(5,4)=xs(2,1)*xs(3,2)+xs(2,2)*xs(3,1)
        t2(5,5)=xs(2,2)*xs(3,3)+xs(2,3)*xs(3,2)
        t2(5,6)=xs(2,1)*xs(3,3)+xs(2,3)*xs(3,1)

        t2(6,1)=xs(1,1)*xs(3,1)
        t2(6,2)=xs(1,2)*xs(3,2)
        t2(6,3)=xs(1,3)*xs(3,3)
        t2(6,4)=xs(1,1)*xs(3,2)+xs(1,2)*xs(3,1)
        t2(6,5)=xs(1,2)*xs(3,3)+xs(1,3)*xs(3,2)
        t2(6,6)=xs(1,1)*xs(3,3)+xs(1,3)*xs(3,1)

C        CALCULATE  T1 MATRIX
        do 60 i=1,6
        do 60 j=1,3
        tc(i,j)=0.d0
        do 60 k=1,6
        tc(i,j)=tc(i,j)-t2(i,k)*c1(k,j)
60        continue
        do 70 i=1,6
        do 70 j=1,3
        tcj(i,j)=0.d0
        do 70 k=1,3
        tcj(i,j)=tcj(i,j)+tc(i,k)*xs(k,j)
70        continue

        do 80 n=1,nel
         do 85 i=1,6
         at1(i)=0.d0
         do 85 k=1,3
         at1(i)=at1(i)+tcj(i,k)*shp(k,n)
85     continue
         do 90 i=1,6
         at2(i)=0.d0
         do 90 k=1,6
90         at2(i)=at2(i)+t2(i,k)*shp2(k,n)
         do 100 i=1,6
       cartd2(i,n)=at1(i)+at2(i)
100        continue
80        continue

       endif

 1000 format('1','non-positive determinant in element number  ',i5,
     &          ' in element material group  ',i5)

        return
      end

c
c     ****************************************************************
c     *                                                              *
c     *                    function    iszero                        *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 5/2/18                      *
c     *                                                              *
c     *           check if x is zero, should be inlined              *
c     *                                                              *
c     ****************************************************************
c
      function isNonZero(x)
      implicit none
      real(8) :: x
      real(8), parameter :: tolp = 1.0D-10, tolm = -1.0D-10
      logical :: isNonZero

      isNonZero = .false.
      if ( x .gt. tolp .or. x .lt. tolm ) isNonZero = .true.

      return
      end
