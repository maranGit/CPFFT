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
      logical :: debug
      real(8), external :: dnrm2
c
      data barF   /one,zero,zero,zero,one,zero,zero,zero,one/
      data barF_t /one,zero,zero,zero,one,zero,zero,zero,one/
      data DbarF  /zero,zero,zero,zero,zero,zero,zero,zero,zero/
      data debug  /.false./
c
      N3_dbl = dble(N3)
c
c                         global step increment
c 
      do step = 1, nstep

        write(*,*) "Now starting step ",step
        ltmstp = step
c
c                         extract boundary condition
c
        PBC = zero
        FBC = zero
        do ii = 1, nstrs
          if ( isNBC(ii) ) then
            PBC( ii ) = BC_all( ii, step )
          else
            FBC( ii ) = BC_all( ii, step )
c
c           DbarF(EBC) = FBC, DbarF(NBC) = extrapolation
            DbarF( ii ) = FBC( ii ) - barF_t( ii )
          end if
        end do
        barF = barF_t + DbarF
        iiter_NBC = 0
c
c                    NR loop for natural boundary condition
c
        do while ( .true. )
c
c         initial residual: distribute "barF" over grid using K4
c
          b = zero
c
c         dFm = repmat( DbarF, N3, 1 )
          dFm = zero
          do ii = 1, 9
            dFm(1:N3, ii) = DbarF(ii)
          end do
          call daxpy(N3*nstrs, one, dFm(1,1), 1, Fn1(1,1), 1)
c
c         get b, solve for first dFm
c
          call G_K_dF(dFm, b, .true.)
          call dscal(N3*nstrs, mone, b(1,1), 1)
          call fftPcg(b, dFm, tolPCG, out)
          call daxpy(N3*nstrs, one, dFm(1,1), 1, Fn1(1,1), 1)
c
c         initialize for most current boundary condition
c
          Fnorm = dnrm2(N3*nstrs, Fn1(1,1), 1)
          resfft = one
          iiter_EBC = 0
c
c         iterate as long as iterative update does not vanish
c
          do while ( resfft .gt. tolNR )

            call drive_eps_sig( step, iiter_EBC )
            call G_K_dF(Pn1, b, .false.)
            call dscal(N3*nstrs, mone, b(1,1), 1)
            call fftPcg(b, dFm, tolPCG, out)
            call daxpy(N3*nstrs, one, dFm(1,1), 1, Fn1(1,1), 1)
            resfft = dnrm2(N3*nstrs, dFm(1,1), 1) / Fnorm
            write(*,*) resfft
            if ( iiter_EBC .eq. maxIter ) then
              write(out,9999)
              call die_abort()
            endif
            iiter_EBC = iiter_EBC + 1

          enddo ! end N-R loop
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
c           || Pbar - PBC || only when PBC(ii) is prescribed
            if( .not. isNBC(ii) ) cycle
            resPbar(1) = resPbar(1) 
     &                 + (P_bar(ii)-PBC(ii)) * (P_bar(ii)-PBC(ii))
c
          end do
          resPbar(3) = sqrt( resPbar(1) / resPbar(2) )
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
          call NBC_update(C_homo(1),DbarF(1),P_bar(1),PBC(1),isNBC(1))
c
          iiter_NBC = iiter_NBC + 1

        end do
c
c            update state variables
c
        DbarF(1:9) = barF(1:9) - barF_t(1:9)
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
      integer :: out
      real(8), intent(in) :: b(veclen), tol
      real(8), intent(out) :: x(veclen)

c     internal variables
      integer, parameter :: nipar = 128, ndpar = 128
      real(8), parameter :: zero = 0.0D0, one = 1.0D0
      integer :: ipar( nipar ), itercount, maxit, RCI_request
      real(8) :: dpar( ndpar ), n2b, resnorm, relres, res( veclen )
      real(8) :: tolb
      logical :: passflg, debug
      real(8), external :: dnrm2
    
c     initialize parameters
      debug = .true.
      maxit = veclen ! maximum number of iteration
      ipar = 0
      itercount = 0
      dpar = zero
      tmpPcg = zero
      passflg = .false.
      x = zero
      n2b = dnrm2( veclen, b, 1 )
      tolb = tol * n2b
      
c                improper tolerance
      if ( tol .le. epsilon(one) .or. tol .ge. one ) then
        write(out,9000) tol
        call die_abort
      endif
c
c               check for all zero right hand side vector
c                   all zero solution
      if ( n2b .eq. zero ) then
        x = zero
        relres = zero
        write(out,1001) itercount, relres
        return
      endif
c
c     initialize the solver
      dpar( 1 ) = tol ! specifies the relative tolerance.
      call dcg_init( veclen, x,b, RCI_request, 
     &               ipar, dpar, tmpPcg )
      ipar( 1 ) = veclen ! length
      ipar( 2 ) = 6 ! the error and warning written to screen
      ipar( 3 ) = 1 ! current stage of the RCI CG computations
      ipar( 5 ) = 1000 ! maximum iteration
      ipar( 8 ) = 1 ! performs stopping test for the maximum iterations
      ipar( 9 ) = 0 ! dcg does not provide stopping test
      ipar( 10 ) = 1 ! I provide stopping test

c     Checks the consistency and correctness of the user defined data
      call dcg_check( veclen, x, b, RCI_request,
     &                ipar, dpar, tmpPcg )
      if (RCI_request .ne. 0) then
        write(out,9001) RCI_request
        call die_abort
      endif

      do while ( .true. )
c       Computes the approximate solution vector
        call dcg(veclen, x, b, RCI_request, ipar, dpar, tmpPcg )
c
c       3 tasks according to RCI_request
c
        select case (RCI_request)
  
        case (1) ! task 1: update solution, A*tmpPcg(:,1) = tmpPcg(:,2)
          call G_K_dF(tmpPcg(:, 1), tmpPcg(:, 2), .true.)
          cycle

        case (2) ! task 2: perform the stopping tests
          res = tmpPcg(:, 3)
          resnorm = dnrm2( veclen, res, 1 )
          if (resnorm.gt.tolb .and. resnorm.gt.tol) then
            passflg = .false.
          else
            call G_K_dF(x, res, .true.)
            ! call DAXPY(veclen, -1.D0, b, 1, res)
            res = b - res
            resnorm = dnrm2( veclen, res, 1 )
            passflg = ( resnorm.le.tolb .or. resnorm.le.tol )
          endif
c
          if ( .not. passflg ) then
c           proceed with CG iterations
            cycle
          else
c           stop CG iterations
            exit
          endif
        
        case (3) ! task 3: apply the preconditioner
          write(out,9002)
          call MKL_FREE_BUFFERS
          call die_abort

        case (0) ! successfully converged
          exit

        case default
c         dcg gives error message, stop the program
          write(out,9999) RCI_request
          call MKL_FREE_BUFFERS
          call die_abort
        end select
      enddo

c     Retrieves the number of the current iteration
      call dcg_get(veclen, x, b, RCI_request,
     &             ipar, dpar, tmpPcg, itercount)
c
      relres = resnorm / n2b
      write(out,1001) itercount, relres
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
 9999 format(5x,'>>> This example FAILED as the solver has',
     &     /,10x,'returned the ERROR code:      ',i2)
 1001 format(/1x,'>>> pcg converged at iteration ', i4, 
     &            ' to a solution with relative residual ',e14.6)
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
      do i = 1, nstrs
        AAA(1:nstrs,i) = C_homo(i,1:nstrs)
      end do
      do i = 1, nstrs
        if ( isNBC(i) ) cycle
        AAA(1:nstrs, i) = zero
        AAA(i, 1:nstrs) = zero
        AAA(i, i)       = one
      end do
c
c                   form b vector
c
      bbb(1:nstrs) = PBC(1:nstrs) - P_bar(1:nstrs)
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
      DbarF(1:nstrs) = zero
      do i = 1, nstrs
        if ( .not. isNBC(i) ) cycle
        DbarF(i) = bbb(i)
      end do
c
      return
      end subroutine
