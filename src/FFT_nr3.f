c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine     FFT_nr3                       *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 1/22/18                     *
c     *                                                              *
c     *                  global Newton-Raphson loop                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine FFT_nr3()
      use fft, only: barF, DbarF, barF_t, b, Fn1, Fn, Pn1, Pn, dFm, 
     &                straininc, tolPCG, tolNR, maxIter, nstep,
     &                mults, F_total, out_step
      implicit none
      include 'common.main'

c                    local
      real(8), parameter :: zero = 0.0D0, one = 1.0D0
      real(8) :: Fnorm, resfft
      integer :: step, iiter, ii
      logical :: debug

      debug = .false.
c 
c                         global step increment
c 
      do step = 1, nstep

        write(*,*) "Now starting step ",step
        ltmstp = step
c
c              average strain increment
c
        DbarF = zero
        do ii = 1, nstrs
          DbarF( :, ii ) = F_total( ii ) * mults( step )
        end do
c
c                      average strain
        barF = barF_t + DbarF

c     initial residual: distribute "barF" over grid using K4
        b = zero
        call G_K_dF(DbarF, b, .true.)
        b = -b

        Fn1 = Fn1 + DbarF
        Fnorm = sqrt(sum(Fn1*Fn1))
        resfft = Fnorm
        iiter = 0

c     iterate as long as iterative update does not vanish
        do while ( .true. )
          call fftPcg(b, dFm, tolPCG, out) ! results stored in dFm in mod_fft
          Fn1 = Fn1 + dFm
          
          call drive_eps_sig( step, iiter )
          
          call G_K_dF(Pn1, b, .false.)
          b = -b
          resfft = sqrt(sum(dFm*dFm))
          write(*,*) resfft/Fnorm
          if ((resfft/Fnorm < tolNR) .and. (iiter > 0)) exit
          if ( iiter .eq. maxIter ) then
            write(out,9999)
            call die_abort()
          endif
          iiter = iiter + 1
          
        enddo ! end N-R loop
c
c            update state variables
c
        Fn = Fn1
        barF_t = barF
        Pn = Pn1
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
