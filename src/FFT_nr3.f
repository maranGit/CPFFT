c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine     FFT_nr3                       *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 6/27/18                     *
c     *                                                              *
c     *                  global Newton-Raphson loop                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine FFT_nr3()
      use fft, only: b, Fn1, Fn, Pn1, Pn, dFm, isNBC, BC_all,
     &               straininc, tolPCG, tolNR, maxIter, nstep,
     &               out_step, N3
      implicit none
      include 'common.main'

c                    local
      real(8), parameter :: zero = 0.0D0, one = 1.0D0, mone = -1.0D0
      real(8) :: P_bar(9), resPbar(3), DbarF(9), barF(9), barF_t(9)
      real(8) :: FBC(9), PBC(9), C_homo(81), Fnorm, resfft, N3_dbl
      integer :: step, ii, jj, iiter_EBC, iiter_NBC
      logical :: debug, existNBC
      real(8), external :: dnrm2
c
      data barF     / one,zero,zero,zero, one,zero,zero,zero, one/
      data barF_t   / one,zero,zero,zero, one,zero,zero,zero, one/
      data DbarF    /zero,zero,zero,zero,zero,zero,zero,zero,zero/
      data P_bar    /zero,zero,zero,zero,zero,zero,zero,zero,zero/
      data debug    /.false./
      data existNBC /.false./
c
      N3_dbl = dble(N3)
c
c                     flat indicating if NBC exists
c            b is borrowed as intermediate array to save memory!
c
      do ii = 1, nstrs
        if ( isNBC(ii) ) existNBC = .true.
      end do
c
c                 initial homogenized tangent stiffness
c
      call tangent_homo( C_homo )
c
c                         global step increment
c 
      do step = 1, nstep

        write(out,1000) step
        ltmstp = step
c
c                         extract boundary condition
c
        PBC = zero
        FBC = zero
        DbarF = zero
        do ii = 1, nstrs
          if ( isNBC(ii) ) then
            PBC( ii ) = BC_all( ii, step )
          else
            FBC( ii ) = BC_all( ii, step )
            DbarF( ii ) = FBC( ii ) - barF_t( ii )
          end if
        end do
c
c                       if NBC exists, update DbarF
c
        if ( existNBC ) then
          call NBC_update(C_homo(1),DbarF(1),P_bar(1),PBC(1),isNBC(1))
        end if
        barF = barF_t + DbarF
        iiter_NBC = 0
c
c                    NR loop for natural boundary condition
c
        do while ( .true. )
c
c         initial residual: distribute "barF" over grid
c         dFm = repmat( DbarF, N3, 1 )
c
          do ii = 1, 9
            dFm(1:N3, ii) = DbarF(ii)
          end do
          Fnorm = dnrm2(N3*nstrs, Fn1(1,1), 1)
          call daxpy(N3*nstrs, one, dFm(1,1), 1, Fn1(1,1), 1)
c
c         get b, solve for first dFm
c
          call G_K_dF(dFm(1,1), b(1,1), .true.)
          call dscal(N3*nstrs, mone, b(1,1), 1)
          call fftPcg(b(1,1), dFm(1,1), tolPCG, out)
          call daxpy(N3*nstrs, one, dFm(1,1), 1, Fn1(1,1), 1)
c
c         initialize for most current boundary condition
c
          resfft = dnrm2(N3*nstrs, dFm(1,1), 1) / Fnorm
          write(out,1003) resfft
          iiter_EBC = 0
          resfft = one ! guarantee at least one iteration
c
c         iterate as long as iterative update does not vanish
c
          do while ( resfft .gt. tolNR )

            call drive_eps_sig( step, iiter_EBC )
            call G_K_dF(Pn1(1,1), b(1,1), .false.)
            call dscal(N3*nstrs, mone, b(1,1), 1)
            call fftPcg(b(1,1), dFm(1,1), tolPCG, out)
            call daxpy(N3*nstrs, one, dFm(1,1), 1, Fn1(1,1), 1)
            resfft = dnrm2(N3*nstrs, dFm(1,1), 1) / Fnorm
            write(out,1001) iiter_EBC, resfft
            if ( iiter_EBC .eq. maxIter ) then
              write(out,9999)
              call die_abort()
            endif
            iiter_EBC = iiter_EBC + 1

          enddo ! end N-R loop
          call drive_eps_sig( step, iiter_EBC )
c
c                   check if P_bar converges to PBC
c
          resPbar = zero
          P_bar = zero
          do ii = 1, nstrs
            do jj = 1, N3
              P_bar(ii) = P_bar(ii) + Pn1(jj,ii)
            end do
            P_bar(ii) = P_bar(ii) / N3_dbl
            resPbar(2) = resPbar(2) + P_bar(ii) * P_bar(ii)
c
c           || Pbar - PBC || only when PBC(ii) is prescribed
c
            if( .not. isNBC(ii) ) cycle
            resPbar(1) = resPbar(1) 
     &                 + (P_bar(ii)-PBC(ii)) * (P_bar(ii)-PBC(ii))
c
          end do
          if ( resPbar(2) .lt. 1.0D-8 ) then
            resPbar(3) = sqrt( resPbar(1) )
          else
            resPbar(3) = sqrt( resPbar(1) / resPbar(2) )
          end if
          write(out,1002) iiter_NBC, resPbar(3)
          if ( resPbar(3) .le. tolNR ) exit
c
c                 check if exceed maximum iteration
c
          if ( iiter_NBC .gt. maxIter ) then
            write(out,9998)
            call die_abort()
          endif
c
c            if P_bar does not converge to PBC, update DbarF
c            compute homogenized tangent stiffness
c            b is borrowed as intermediate array to save memory!
c
          call tangent_homo( C_homo )
          DbarF = zero
          call NBC_update(C_homo(1),DbarF(1),P_bar(1),PBC(1),isNBC(1))
          barF = barF + DbarF
c
          iiter_NBC = iiter_NBC + 1

        end do
c
c            update state variables
c
        barF_t(1:9) = barF(1:9)
        call dcopy(N3*nstrs, Fn1(1,1), 1, Fn(1, 1), 1)
        call dcopy(N3*nstrs, Pn1(1,1), 1, Pn(1, 1), 1)
c
c            update global variables in eleblocks
c
        call update
c
c            output file
c
        if(out_step(step)) call ouresult(step)

      enddo ! end all steps
c
      if ( debug ) then
        write(out,*) "check 1st PK stress"
        write(out,*) Pn1(4000:29000:5000, 1:9)
      end if
      
 9999 format(/,'>> Error: Newton loop does not converge.')
 9998 format(/,'>> Error: Prescribed stress cannot be reached ',
     &         'within given maximum iteration.'/)
 1000 format(/,4x,'--------------------------------------------------',
     &          '-----------------------',/,5x,'Now starting step: ',i7)
 1001 format(7x,'Iteration ', i5, 7x, '  residual ', d10.3)
 1002 format(7x,'Stress iteration ',i5, '  residual ',d10.3)
 1003 format(7x,'Initial residual ', 16x, d10.3)
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine       fftPcg                      *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 1/22/18                     *
c     *                                                              *
c     *                  conjugate gradient algorithm                *
c     *                                                              *
c     ****************************************************************
c
      subroutine fftPcg(b, x, tol, out)
      use fft, only: veclen, tmpPcg
      implicit none

c     input variables
      integer, intent(in)  :: out
      real(8), intent(in)  :: b(*), tol
      real(8), intent(out) :: x(*)

c     internal variables
      integer, parameter :: maxIter = 1000
      integer, parameter :: nipar = 128, ndpar = 128
      real(8), parameter :: zero = 0.0D0, one = 1.0D0
      integer :: ipar( nipar ), itercount, RCI_request
      real(8) :: dpar( ndpar ), n2b, resnorm, relres
      real(8) :: tolb
      logical :: passflg, debug
      real(8), external :: dnrm2
c
      data debug /.false./
c
c     initialize parameters
c
      x(1:veclen) = zero
      ipar = 0
      itercount = 0
      dpar = zero
      passflg = .false.
      n2b = dnrm2( veclen, b(1), 1 )
      tolb = tol * n2b
c
c                improper tolerance
c
      if ( tol .le. epsilon(one) .or. tol .ge. one ) then
        write(out,9000) tol
        call die_abort
      endif
c
c               check for all zero right hand side vector
c                   all zero solution
c
      if ( n2b .le. epsilon(one) ) then
        x(1:veclen) = zero
        write(out,1002)
        return
      endif
c
c     initialize the solver
c     I cannot remember why I use this order ...
c
      dpar( 1 ) = tol ! specifies the relative tolerance.
      call dcg_init( veclen, x(1), b(1), RCI_request,
     &               ipar, dpar, tmpPcg(1,1) )
      ipar( 1 ) = veclen ! length
      ipar( 2 ) = 6 ! the error and warning written to screen
      ipar( 3 ) = 1 ! current stage of the RCI CG computations
      ipar( 5 ) = maxIter ! maximum iteration
      ipar( 8 ) = 1 ! performs stopping test for the maximum iterations
      ipar( 9 ) = 0 ! dcg does not provide stopping test
      ipar( 10 ) = 1 ! I provide stopping test
c
c     Checks the consistency and correctness of the user defined data
c
      call dcg_check( veclen, x(1), b(1), RCI_request,
     &                ipar, dpar, tmpPcg(1,1) )
c
      if (RCI_request .ne. 0) then
        write(out,9001) RCI_request
        call die_abort
      endif

      do while ( .true. )
c
c       Computes the approximate solution vector
c
        call dcg(veclen, x(1), b(1), RCI_request,
     &           ipar, dpar, tmpPcg(1,1) )
c
c       3 tasks according to RCI_request
c
        select case (RCI_request)
  
        case (1) ! task 1: update solution, A*tmpPcg(:,1) = tmpPcg(:,2)
c
          call G_K_dF(tmpPcg(1, 1), tmpPcg(1, 2), .true.)
          cycle
c
        case (2) ! task 2: perform the stopping tests
c
          resnorm = dnrm2( veclen, tmpPcg(1, 3), 1 )
          passflg = ( resnorm .le. tolb .or. resnorm .le. tol )
          if ( passflg ) exit
c
        case (3) ! task 3: apply the preconditioner
c
          write(out,9002)
          call MKL_FREE_BUFFERS
          call die_abort
c
        case (0) ! successfully converged
c
          exit
c
        case default
c
c         dcg gives error message, stop the program
          write(out,9999) RCI_request
          call MKL_FREE_BUFFERS
          call die_abort
c
        end select
c
      enddo

c     Retrieves the number of the current iteration
      call dcg_get(veclen, x(1), b(1), RCI_request,
     &             ipar, dpar, tmpPcg(1,1), itercount)
c
      relres = resnorm / n2b
      if ( debug ) write(out,1001) itercount, relres
      if ( itercount .ge. maxIter ) then
        write(out,9003) maxIter
        call MKL_FREE_BUFFERS
        call die_abort()
      end if
c
c              release memory in dcg, otherwise memory leak
c
      call MKL_FREE_BUFFERS
      return
c
 9000 format(5x,'>>> pcg: improper tolerance',D9.2)
 9001 format(5x,'>>>fftPcg: dcg_check failed'
     &     /,10x,'returned the ERROR code:      ',i2)
 9002 format(5x,'>>>fftPcg: precondition is not available now')
 9003 format(5x,'>>>fftPcg: fail to converge within ',i6,' iterations')
 9999 format(5x,'>>> This example FAILED as the solver has',
     &     /,10x,'returned the ERROR code:      ',i2)
 1001 format(/1x,'>>> pcg converged at iteration ', i4, 
     &            ' to a solution with relative residual ',e14.6)
 1002 format(/1x,'>>> WARNING: Initial RHS is too small. ',
     &           'Returning zero solution ...')
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine   NBC_update                      *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 6/28/18                     *
c     *                                                              *
c     *              update delta_F for prescribed P                 *
c     *              P_bar is old P, PBC is desired P                *
c     *                                                              *
c     ****************************************************************
c
      subroutine NBC_update( C_homo, DbarF, P_bar, PBC, isNBC )
      implicit none
      include 'common.main'
c
c                        global
c
      real(8), intent(in)  :: C_homo(nstrs,*), P_bar(*), PBC(*)
      logical, intent(in)  :: isNBC(*)
      real(8), intent(out) :: DbarF(*)
c
c                        local
c
      real(8) :: AAA(nstrs,nstrs), bbb(nstrs)
      integer :: ipiv(nstrs), info, i
      real(8), parameter :: zero = 0.0D0, one = 1.0D0
c
c                   form A matrix
c
c     do i = 1, nstrs
c       AAA(1:nstrs,i) = C_homo(i,1:nstrs)
c     end do
      do i = 1, nstrs
c       if ( isNBC(i) ) cycle
c       AAA(1:nstrs, i) = zero
c       AAA(i, 1:nstrs) = zero
c       AAA(i, i)       = one
        if ( isNBC(i) ) then
          AAA(i,1:nstrs) = C_homo(1:nstrs,i)
          bbb(i) = PBC(i) - P_bar(i)
        else
          AAA(i,1:nstrs) = zero
          AAA(i, i) = one
          bbb(i) = DbarF(i)
        end if
      end do
c
c                   form b vector
c
c     bbb(1:nstrs) = PBC(1:nstrs) - P_bar(1:nstrs)
c
c                      solve
c
      call dgesv(nstrs,1,AAA(1,1),nstrs,ipiv(1),bbb(1),nstrs,info)
      if ( info .ne. 0 ) then
        write(out,*) ">>> Error: P_bar update failed"
        call die_abort()
      end if
c
c                   return DbarF
c
c     DbarF(1:nstrs) = zero
c     do i = 1, nstrs
c       if ( .not. isNBC(i) ) cycle
c       DbarF(i) = bbb(i)
c     end do
      DbarF(1:nstrs) = bbb(1:nstrs)
c
      return
      end subroutine
