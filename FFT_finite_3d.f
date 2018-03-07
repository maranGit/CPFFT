      program main
      use fft
      use elem_block_data
      implicit none
      include 'common.main'
      
      logical debug
      
c     temporary variables
      integer :: blk, iiter, step, nstep
      real(8) :: Fnorm, resfft
      integer :: phase(31,31,31)
      integer :: now_thread
      integer :: alloc_stat
      integer, external :: omp_get_thread_num
      real(8), parameter :: zero = 0.0D0, one = 1.0D0
      real(8), parameter :: tolNR = 1.0D-5, tolPCG = 1.0D-8
      
c     problem size
      ndim1 = 3 ! defined in param_def
      N = 31
      nstep = 1
      phase = 0
      phase(23:31, 1:9, 23:31) = 1
      
c     initialize parameters in common block
      step = 0
      ndim2 = ndim1 * ndim1
      ndim3 = ndim2 * ndim1
      ndim4 = ndim3 * ndim1
      N3 = N * N * N
      noelem = N3
      if ( mod(N, 2) .eq. 1 ) Nhalf = (N + 1) / 2
      if ( mod(N, 2) .eq. 0 ) Nhalf = N / 2 + 1
      veclen = N3 * 9
      dims = [N, N, N]
      debug = .true.
      open(out, file = "RanOut.out")
c
c     autoblock for openmp, store in elblks (elem_block_data.mod)
      matList = reshape( phase, (/N3/))
      call inelbk( matList )
c
c     allocate variables in module
c
      call fftAllocate( 1 )
c
c     initialize deformation gradient
      Fn = zero
      Fn(:,[1,5,9]) = one
      Fn1 = Fn
      Pn = zero
      Pn1 = zero
c
c     form G_hat_4 matrix and store in Ghat4
c
      call formG()
c
c     NEWTON ITERATIONS
c     loop over blocks to update stress and tangent stiffness matrix
c$OMP PARALLEL DO PRIVATE( blk, now_thread )
c$OMP&         SHARED( Fn1, matList, Pn1, K4, N3, nelblk, iiter, step )
      do blk = 1, nelblk
        now_thread = omp_get_thread_num() + 1
        call do_nleps_block( blk, iiter, step )
c        call do_nleps_block( blk, iter, step, step_cut_flags(blk),
c     &                       block_energies(blk),
c     &                       block_plastic_work(blk) )
      enddo
c$OMP END PARALLEL DO
      
      do step = 1, nstep
c     macroscopic loading
        DbarF = zero
        DbarF(:,2) = one
        b = zero

c     initial residual: distribute "barF" over grid using K4
        call G_K_dF(DbarF, b, .true.)
        b = -b
        Fn1 = Fn1 + DbarF
        Fnorm = sqrt(sum(Fn1*Fn1))
        resfft = Fnorm
        iiter = 0

c     iterate as long as iterative update does not vanish
        do while ( resfft/Fnorm .gt. tolNR )
          call fftPcg(b, dFm, tolPCG, out) ! results stored in dFm in mod_fft
          Fn1 = Fn1 + dFm
c$OMP PARALLEL DO PRIVATE( blk, now_thread )
c$OMP&         SHARED( Fn1, matList, Pn1, K4, N3, nelblk, iiter, step )
          do blk = 1, nelblk
            now_thread = omp_get_thread_num() + 1
            call do_nleps_block( blk, iiter, step )
c            call do_nleps_block( blk, iter, step, step_cut_flags(blk),
c     &                           block_energies(blk),
c     &                           block_plastic_work(blk) )
          enddo
c$OMP END PARALLEL DO
          call G_K_dF(Pn1, b, .false.)
          b = -b
          resfft = sqrt(sum(dFm*dFm))
          write(*,*) resfft/Fnorm
          iiter = iiter + 1
          
        enddo ! end N-R loop
c
c            update state variables
        Fn = Fn1
        Pn = Pn1
c            update global variables in eleblocks
      call update

      enddo ! end all steps
c
c     Ran check results
c
      if ( debug ) then
        write(out,*) "check 1st PK stress"
        write(out,*) Pn1(4000:29000:5000, 1:9)
      end if
c
c     deallocate all arrays
c
      call fftAllocate( 2 )
c
c     close output file
c
      close(out)
      
 1001 format( 27f8.4 )
c
      end program
c     ****************************************************************
c     *                                                              *
c     *                      subroutine do_nleps_block               *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified : 02/26/2018 RM              *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine do_nleps_block( blk, iter, step )
      use fft, only: Fn1, Pn1, Fn, Pn, K4, matList
      use elem_block_data
      implicit none
      include 'common.main'
      include 'include_sig_up'
      
      integer, intent(in) :: blk, iter, step
c
c             locals
c
      integer :: ii, span, felem, currElem
      integer :: info_vector(4)
      logical :: local_debug
      
      integer :: ngp, cep_size, block_size, gpn
      integer :: hist_size
      real(8), parameter :: zero = 0.0D0, one = 1.0D0, half = 0.5D0
      double precision :: gp_dtemps(mxvl)
      integer :: alloc_stat
      integer :: mat_type
      logical :: geo_non_flg
      real(8), allocatable, dimension(:,:,:) :: rnh, fnh, dfn, fnhinv
      real(8), allocatable, dimension(:) :: detF
      real(8), allocatable, dimension(:,:) :: ddt, uddt
      real(8), allocatable, dimension(:,:,:) :: qnhalf, qn1
      real(8), allocatable, dimension(:,:) :: P_ref, cep_ref
c
c             cauchu stress and rotation matrix @ n+1
      real(8) :: qtn1(mxvl,nstr,nstr), cs_blk_n1(mxvl,nstr)

c             allocate and initialization
      allocate( detF(mxvl) )
      allocate( rnh(mxvl,ndim,ndim), fnh(mxvl,ndim,ndim), 
     &          dfn(mxvl,ndim,ndim), fnhinv(mxvl,ndim,ndim) )
c             warp3d original allocate and initialization
c             nstr = 6; nstrs = 9;
      allocate( ddt(mxvl,nstr), uddt(mxvl,nstr),
     &          qnhalf(mxvl,nstr,nstr), qn1(mxvl,nstr,nstr) )
      allocate( P_ref(mxvl,nstrs), cep_ref(mxvl,nstrs*nstrs) )
      
      local_debug = .false.
!DIR$ VECTOR ALIGNED
      ddt    = zero
!DIR$ VECTOR ALIGNED
      uddt   = zero
!DIR$ VECTOR ALIGNED
      qnhalf = zero
!DIR$ VECTOR ALIGNED
      qn1    = zero
!DIR$ VECTOR ALIGNED
      detF   = zero
!DIR$ VECTOR ALIGNED
      rnh    = zero
!DIR$ VECTOR ALIGNED
      fnh    = zero
!DIR$ VECTOR ALIGNED
      dfn    = zero
!DIR$ VECTOR ALIGNED
      fnhinv = zero
!DIR$ VECTOR ALIGNED
      P_ref  = zero
!DIR$ VECTOR ALIGNED
      cep_ref= zero
c
c     initialize local_work based on elblks(0,blk) and elblks(1,blk)
c
      ngp = 1
      gpn = 1
      mat_type = 1
      geo_non_flg = .true.
      span = elblks(0, blk)
      felem = elblks(1, blk)
      local_work%blk = blk
      local_work%span = span
      local_work%felem = felem
      local_work%num_int_points = ngp
      local_work%gpn = gpn
      local_work%step = step
      local_work%iter = iter ! global newton iteration (>1)
      local_work%geo_non_flg = geo_non_flg
      local_work%is_cohes_elem = .false.
      local_work%mat_type = mat_type
      local_work%segmental = .false.
      local_work%number_points = 0
      local_work%curve_set_number = 0
      local_work%fgm_enode_props = .false.
      local_work%killed_status_vec = .false.
      local_work%iout = out
      local_work%int_order = 1
      local_work%num_enodes = 8
      local_work%num_enode_dof = 3
      local_work%lnelas_vec = .false.
      local_work%bbar_flg = .true.
      local_work%material_cut_step = .false.
      local_work%adaptive_flag = .false.
      local_work%eps_bbar = 0
c
c     allocate memory for local_work according to material model
c
      call mm01_set_sizes( info_vector )
      hist_size = info_vector(1)
      cep_size = info_vector(2)
      block_size = span * ngp * cep_size
      local_work%hist_size_for_blk = hist_size
      call recstr_allocate( local_work )
      local_work%elem_type = 2 ! lsdisop element
c
c     grab material parameters, state variables and history
c     due(nodal displacement from n to n+1)
c     ue(nodal displacement from 0 to n)
c     ce_0, cd_n, cd_mid, cd_n1: nodal coordinate
c                 ***hard coded for now***
c
      call mm01_hardCoded ! only parameters and temperature
      call dupstr_blocked( blk, span, felem, mat_type, 
     & geo_non_flg, step, iter, local_work )
c
c                grab global Fn and Fn1 to local_work
c
!DIR$ LOOP COUNT MAX=128
!DIR$ VECTOR ALIGNED
      do ii = 1, span
        local_work%fn(ii,1,1) = Fn(ii,1)
        local_work%fn(ii,1,2) = Fn(ii,2)
        local_work%fn(ii,1,3) = Fn(ii,3)
        local_work%fn(ii,2,1) = Fn(ii,4)
        local_work%fn(ii,2,2) = Fn(ii,5)
        local_work%fn(ii,2,3) = Fn(ii,6)
        local_work%fn(ii,3,1) = Fn(ii,7)
        local_work%fn(ii,3,2) = Fn(ii,8)
        local_work%fn(ii,3,3) = Fn(ii,9)
      end do
!DIR$ LOOP COUNT MAX=128
!DIR$ VECTOR ALIGNED
      do ii = 1, span
        local_work%fn1(ii,1,1) = Fn1(ii,1)
        local_work%fn1(ii,1,2) = Fn1(ii,2)
        local_work%fn1(ii,1,3) = Fn1(ii,3)
        local_work%fn1(ii,2,1) = Fn1(ii,4)
        local_work%fn1(ii,2,2) = Fn1(ii,5)
        local_work%fn1(ii,2,3) = Fn1(ii,6)
        local_work%fn1(ii,3,1) = Fn1(ii,7)
        local_work%fn1(ii,3,2) = Fn1(ii,8)
        local_work%fn1(ii,3,3) = Fn1(ii,9)
      end do
      fnh = half * ( local_work%fn + local_work%fn1 )
      dfn = local_work%fn1 - local_work%fn
c
c              compute rotation tensor
c
      call rtcmp1( span, fnh, rnh )
      call rtcmp1( span, local_work%fn1,
     &             local_work%rot_blk_n1(1,1,gpn) )
c
c              compute ddt( detF is a dummy argument )
c
      call inv33( span, felem, fnh, fnhinv, detF )
      call mul33( span, felem, dfn, fnhinv, ddt )
c
c              compute uddt
c
      call getrm1( span, qnhalf, rnh, 1 ) ! form qnhalf
      call qmply1( span, mxvl, nstr, qnhalf, ddt, uddt )
      call rstgp1_update_strains( span, mxvl, nstr, uddt,
     &                            local_work%ddtse(1,1,gpn) )
c
c     recover stress and stiffness
c
      call rstgp1( local_work, uddt )
      if(local_debug) write(out,*) uddt(:,1)
c
c     scatter local variables to global
c
      call rplstr( span, felem, ngp, mat_type, iter,
     &             geo_non_flg, local_work, blk )

c     compute cs_blk_n1(Cauchy stress) from urcs_blk_n1
c
      call getrm1( span, qtn1, local_work%rot_blk_n1(1,1,gpn), 2 )
      call qmply1( span, mxvl, nstr, qtn1,
     &             local_work%urcs_blk_n1(1,1,gpn),
     &             cs_blk_n1 )
c     pull back cs_blk_n1(Cauchy stress) to 1st PK stress
c     P = J*sigma*F^{-T}
      call inv33( span, felem, fnh, fnhinv, detF )
      call mul33( span, felem, dfn, fnhinv, ddt )

c     pull back cep_blocks to dP/dF
c     cep_blocks in elem_block_data.mod
c
      ! do nothing for now
c
c     update global stress and stiffness
c
      ! do nothing for now
c
c     get cut-step flag and update plastic work
c
      do ii = 1, span
        currElem = felem + ii - 1
        call constitutive(Fn1(currElem,:), matList(currElem), 
     &                     Pn1(currElem,:), K4(currElem,:))
      enddo
c
      deallocate( rnh, fnh, dfn, fnhinv )
      deallocate( ddt, uddt, qnhalf, qn1 )
      deallocate( P_ref, cep_ref )
      call recstr_deallocate( local_work )
      contains
c     ========
c
c     ****************************************************************
c     *                                                              *
c     *                  subroutine mm01_hardCoded                   *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified : 2/07/2018 rhd              *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine mm01_hardCoded
      implicit none 
c
c                hard code uddt, uncs and gp_dtemps
c
c     do ii = 1, span
c       local_work%urcs_blk_n(ii,:,gpn) = 
c    &            [63.369052D0,      -1.571854D-003, -1.571854D-003,
c    &             -7.384044D-005,  9.452432D-006, -7.384048D-005,
c    &              0.615002D0,       0.548073D0,       1.112306D-002]
c       uddt(ii,:) = 
c    &            [6.452949D-003, -1.940361D-003, -1.940361D-003,
c    &             5.700588D-006,  2.689669D-007,  5.700588D-006]
c     enddo
      
      gp_dtemps = zero ! temperature change over step
c     local_work%urcs_blk_n1(1,1,gpn) = zero ! cartesian stress for step n+1
      
      ! material parameters
      local_work%e_vec = 30000.0D0 ! young's modulus in the block
      local_work%nu_vec = 0.3D0 ! poisson's ratio
      local_work%beta_vec = one ! isotropic/kinematic fractional factor for elements in the block
      local_work%h_vec = 300.0D0 ! plastic hardening modulus
      local_work%sigyld_vec = 60.0D0 ! uniaxial yield stress
      
      ! all 11 history variables
c     local_work%elem_hist(:,1,:) = 7.9899D-3 ! lamda * dt @ n
c     local_work%elem_hist(:,1,:) = 7.9868D-003
c     local_work%elem_hist(:,2,:) = 36.587D0
c     local_work%elem_hist(:,3,:) = 1.1123D-002
c     local_work%elem_hist(:,4,:) = 8.4488D-004
c     local_work%elem_hist(:,5,:) = 303.03D0
c     local_work%elem_hist(:,6,:) = zero
c     local_work%elem_hist(:,7,:) = zero
c     local_work%elem_hist(:,8,:) = zero
c     local_work%elem_hist(:,9,:) = zero
c     local_work%elem_hist(:,10,:) = zero
c     local_work%elem_hist(:,11,:) = zero
      local_work%rtse(1,1,gpn) = zero ! "relative" trial elastic stress rate
      local_work%e_vec_n = 30000.0D0 ! young's modulus at step n
      local_work%nu_vec_n = 0.3D0 ! poisson's ratio at step n
        
c     other properties required by rstgp1
      local_work%block_energy = zero
      local_work%block_plastic_work = zero
      local_work%beta_fact = zero
      
c     other properties required by drive_01_udate
      local_work%dt = one
      local_work%temperatures = 297.0D0
      local_work%temperatures_ref = 297.0D0

      end subroutine ! mm01_hardCoded
      end subroutine ! do_nleps_block
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine      fftAllocate                  *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 2/3/18                      *
c     *                                                              *
c     *                  allocate array in mod_fft.f                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine fftAllocate( action )
      use fft
      use elem_block_data
      implicit none
      include 'common.main'
      
c     temporary variables
      integer :: blk
      integer :: action, iok, N3ndim2
      
      N3ndim2 = N3 * ndim2
      
      select case ( action )
      case( 1 )
c
c     allocate variables in module
c
        allocate ( Ghat4(N3, ndim4), stat=iok )
        allocate ( K4(N3, ndim4), stat=iok )
        allocate ( Fn(N3, ndim2), stat=iok )
        allocate ( Fn1(N3, ndim2), stat=iok )
        allocate ( Pn(N3, ndim2) )
        allocate ( Pn1(N3, ndim2) )
        allocate ( DbarF(N3, ndim2) )
        allocate ( b(N3, ndim2) )
        allocate ( dFm(N3, ndim2) )
  
c     allocate FFT related variables
        allocate ( real1(N3) )
        allocate ( cplx3half(Nhalf,N,N) )
        allocate ( cplx1half(Nhalf*N*N) )
  
c     allocate pcg related variables
        allocate ( tmpPcg(N3ndim2, 4) )
  
c     allocate internal variables
        allocate( tmpReal(N3, ndim2) )
        allocate( tmpCplx(N3, ndim2) )

c     allocate global arrays in Warp3d
        call history_cep_init()
        call stresses_init()
        call rotation_init()
        call strains_init()
c
      case( 2 )
c
c     deallocate variables in module
c
        deallocate( Ghat4, K4, Fn, Fn1 )
        deallocate( Pn, Pn1, DbarF, b, dFm )
        deallocate( real1, cplx3half, cplx1half )
        deallocate( tmpPcg )
        deallocate( tmpReal, tmpCplx )
        deallocate( history_blk_list )
        do blk = 1, nelblk
          deallocate( cep_blocks(blk)%vector )
        end do
        deallocate( cep_blocks )

      case default
        write(*,*) ">>Error: invalid option. Job terminated"
      end select

      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine      inelbk                       *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 1/22/18                     *
c     *                                                              *
c     *     automatic generation of element blocking: simple version *
c     *     for threaded w/o vectorized blocking. optional           *
c     *     assignment of blocks to domains w/ simple algorithm      *
c     *                                                              *
c     ****************************************************************
c
      subroutine inelbk( matList )
      implicit none
      include 'common.main'
c
      integer :: auto_size
      logical :: display
      integer :: felem, current_size, element, param, i
      logical :: newblk, compatible
      integer :: blk_matmodel, ele_matmodel
      integer, intent(in) :: matList(noelem)
c
c     Ran: hard coded variables for autoblock
c
      auto_size    = mxvl  !  from param_def
      display      = .true.
c
c                     first generation of automatic assignment of
c                     elements to blocks.
c
c                     a) sequential pass thru elements
c                     b) assign to current block or open new
c                        block if full or element/material
c                        combinations are not compatible with
c                        current block.
c                     c) no element renumbering or vectorized
c                        blocking in this first set of
c                        features.
c
c
      nelblk       = 1   ! in common main
      current_size = 1
      felem        = 1
      elblks(1,1)  = 1  ! first element in block
      elblks(0,1)  = 1  ! number elements in block
      if( noelem == 1 ) return
c
      blk_matmodel = matList( 1 )
c
c                     1. element fits in current block?
c                     2. if yes, is it compatible with elements now
c                        in the block?
c                     3. if yes, update current number of elements
c                        in the block, next element
c                     4. otherwise start a new block. set first
c                        element in the block, init block size,
c                        load props for first element in block for
c                        subsequent comparisons
c
      do element = 2, noelem
         newblk       = .false.
         current_size = current_size + 1
         if( current_size .gt. auto_size ) newblk = .true.
         if( .not. newblk ) then
           ele_matmodel = matList( element )
           compatible = ( blk_matmodel .eq. ele_matmodel )
           if( .not. compatible ) newblk = .true.
         end if
         if( .not. newblk ) then
           elblks(0,nelblk)  = current_size
           cycle
         endif
         nelblk = nelblk + 1
         if( nelblk .gt. mxnmbl ) then
            param = nelblk
c      call errmsg(74,param,dums,dumr,dumd)
c      call die_abort
            write(out,*) "too many element blocks required"
            stop
         end if
         felem             = element
         elblks(1,nelblk)  = felem
         elblks(0,nelblk)  = 1
         current_size       = 1
         blk_matmodel = matList( felem )
      end do ! on element
c
c                     display blocking table if requested
c
      if( .not. display ) return
      write(out,9000) nelblk, auto_size
      write(out,9010)
      do i = 1, nelblk
        write(out,9020) i, elblks(1,i),  elblks(0,i), elblks(2,i)
      end do
      write(out,*) ' '
c
      return
c
 9000 format(/,'>> Generated element blocking table:',
     & /,      1x,'  number of blocks, target size: ',i7,i5)
 9010 format(/,5x,
     &'block        1st element in blk        # elements in block',
     &5x,'assigned domain' )
 9020 format(1x,i8, 10x,i9,25x,i4,10x,i6)

      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine       formG                       *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 1/22/18                     *
c     *                                                              *
c     *                  form G_hat_4 matrix                         *
c     *                                                              *
c     ****************************************************************
c
      subroutine formG()
      use fft, only: N, Nhalf, Ghat4, N3, ndim1, ndim2
      implicit none

c     internal variables
      integer :: ii, jj, kk, tmp, Gnon0(27)
      real(8), dimension(:,:), allocatable :: q
      real(8), dimension(:), allocatable :: q_dot_q

c     allocate internal variables
      allocate ( q(N3, ndim1) )
      allocate ( q_dot_q(N3) )

      tmp = 1
      do ii = 1, N
        do jj = 1, N
          do kk = 1, N
            q(tmp, :) = dble( [ ii-Nhalf, jj-Nhalf, kk-Nhalf ] )
            tmp = tmp + 1
          enddo
        enddo
      enddo
      q_dot_q = q(:,1)*q(:,1) + q(:,2)*q(:,2) + q(:,3)*q(:,3)

      do ii = 1, 27
        Gnon0(ii) = ii * 3 - 2
      enddo
      Gnon0( 10:18 ) = Gnon0( 10:18 ) + 1
      Gnon0( 19:27 ) = Gnon0( 19:27 ) + 2
      do kk = 1, 3
        do ii = 1, 3
          do jj = 1, 3
            tmp = (kk - 1) * 9 + (ii - 1) * 3 + jj
            Ghat4(:,Gnon0(tmp)) = q(:,ii) * q(:,jj)
          enddo
        enddo
      enddo
      do ii = 1, N3
        if ( abs( q_dot_q(ii) ) .le. 1D-6 ) then
          Ghat4(ii, :) = 0.0D0
        else
          Ghat4(ii, :) = Ghat4(ii, :) / q_dot_q(ii)
        endif
      enddo

c     deallocate internal variables
      deallocate( q, q_dot_q )
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
      real(8) :: b(veclen), x(veclen), tol

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
      ipar( 5 ) = veclen ! maximum iteration
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
          if (resnorm .gt. tolb) then
            passflg = .false.
          else
            call G_K_dF(x, res, .true.)
            ! call DAXPY(veclen, -1.D0, b, 1, res)
            res = b - res
            resnorm = dnrm2( veclen, res, 1 )
            passflg = ( resnorm .le. tolb )
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
c     *                 subroutine   G_K_dF                          *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 1/22/18                     *
c     *                                                              *
c     *                                                              *
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

c     internal variables
      integer :: ii

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

c     multiply by G_hat matrix
      call ddot42n_cmplx(Ghat4, tmpCplx, N3)

c     inverse fft
      do ii = 1, ndim2
        call ifftfem3d(tmpCplx(:,ii), tmpReal(:,ii))
      enddo

c     return
      GKF = tmpReal

      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine       constitutive                *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 1/22/18                     *
c     *                                                              *
c     *                                                              *
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
      end subroutine ! consititutive
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
     &                 conjg( cplx3half( Nhalf:2:-1, 1     , 1      ) )
      cplx3( Nhalf+1:N, 1  , 2:N ) = 
     &                 conjg( cplx3half( Nhalf:2:-1, 1     , N:2:-1 ) )
      cplx3( Nhalf+1:N, 2:N, 1   ) = 
     &                 conjg( cplx3half( Nhalf:2:-1, N:2:-1, 1      ) )
      cplx3( Nhalf+1:N, 2:N, 2:N ) = 
     &                 conjg( cplx3half( Nhalf:2:-1, N:2:-1, N:2:-1 ) )
      cplx3                        = cplx3(ffts, ffts, ffts)

c     deallocate local variables
      deallocate( ffts, iffts )

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

      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine trans2, dot22, dot24, dot42       *
c     *                            ddot42, ddot42n, ddot44           *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 1/22/18                     *
c     *                                                              *
c     *                  frequently used tensor operation            *
c     *                                                              *
c     ****************************************************************
c
      subroutine trans2(A2, n)
      implicit none
c     B2_ji = A2_ij, n rows
      integer, intent(in) :: n
      real(8) :: A2( n, 9 )
      A2 = A2( :, [ 1, 4, 7, 2, 5, 8, 3, 6, 9 ] )
      end subroutine

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
      end subroutine

      subroutine ddot42n(A4, B2, n)
      implicit none
c     A4_ijkl * B2_lk = C2_ij, n rows
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

      end subroutine

      subroutine ddot42n_cmplx(A4, B2, n)
      implicit none
c     A4_ijkl * B2_lk = C2_ij, n rows
c     A4 is real, B2 is complex, C2 is complex
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
      end subroutine
c     ****************************************************************
c     *                                                              *
c     *                      subroutine inv33                        *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                                                              *
c     *                compute inverse matrice of                    *
c     *                   en element block                           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine inv33( span, felem, jac, gama, dj )
      implicit none
c
      include 'param_def'
c                   parameters
c
      double precision ::
     &  jac(mxvl,3,3), gama(mxvl,3,3)
      integer :: span, felem
c
      double precision :: dj(mxvl) 
c
c                   local work arrays (on stack)
c
      integer :: i, row, col
      logical :: local_debug
      double precision ::
     &   j1(mxvl), j2(mxvl), j3(mxvl)
c
      double precision ::
     &       zero, zero_check, one, half
c
c                 hard coded material parameters
c
      integer :: gpn
c
      data zero, zero_check, one, half, local_debug
     &    / 0.0d0,    1.0d-20, 1.0d0,  0.5d0, .false. /
!DIR$ VECTOR ALIGNED
      gama = zero
c
c           calculate the determinate of the jacobian matrix
c
      gpn = 1
!DIR$ LOOP COUNT MAX=128  
!DIR$ VECTOR ALIGNED
      do i = 1, span
          j1(i)= jac(i,2,2)*jac(i,3,3)-jac(i,2,3)*jac(i,3,2)
          j2(i)= jac(i,2,1)*jac(i,3,3)-jac(i,2,3)*jac(i,3,1)
          j3(i)= jac(i,2,1)*jac(i,3,2)-jac(i,2,2)*jac(i,3,1)
          dj(i)= jac(i,1,1)*j1(i)-jac(i,1,2)*j2(i)+jac(i,1,3)*j3(i)
      end do
c
c           check to insure a positive determinate.
c
!DIR$ LOOP COUNT MAX=128  
!DIR$ IVDEP
      do i = 1, span
       if( dj(i) .le. zero_check ) then
         write(*,9169) gpn,felem+i-1, dj(i)
       end if
      end do
c
c           calculate the inverse of the jacobian matrix
c
!DIR$ LOOP COUNT MAX=128  
!DIR$ VECTOR ALIGNED
      do i = 1, span
         gama(i,1,1)=  j1(i)/dj(i)
         gama(i,2,1)= -j2(i)/dj(i)
         gama(i,3,1)=  j3(i)/dj(i)
         gama(i,1,2)= (jac(i,3,2)*jac(i,1,3)-
     &                 jac(i,1,2)*jac(i,3,3))/dj(i)
         gama(i,2,2)= (jac(i,1,1)*jac(i,3,3)-
     &                 jac(i,3,1)*jac(i,1,3))/dj(i)
         gama(i,3,2)= (jac(i,1,2)*jac(i,3,1)-
     &                 jac(i,1,1)*jac(i,3,2))/dj(i)
         gama(i,1,3)= (jac(i,1,2)*jac(i,2,3)-
     &                 jac(i,1,3)*jac(i,2,2))/dj(i)
         gama(i,2,3)= (jac(i,1,3)*jac(i,2,1)-
     &                 jac(i,1,1)*jac(i,2,3))/dj(i)
         gama(i,3,3)= (jac(i,1,1)*jac(i,2,2)-
     &                 jac(i,1,2)*jac(i,2,1))/dj(i)
      end do
c
      if ( local_debug ) then
        do i=1, span
          write(*,*)'       Jacobian matrix, elem #',i
          do row=1,3
            write(*,9000)(jac(i,row,col),col=1,3)
          end do
          write(*,*)'       determinant of the Jacobian matrix =',dj(i)
          write(*,*)'       inverse Jacobian matrix'
          do row=1,3
            write(*,9000)(gama(i,row,col),col=1,3)
          end do
        end do
      end if
c
      return
c
 9000 format(3x,3f10.5)
 9900 format(I10,3F10.4)
 9169 format(/1x,'>>>>> warning: the determinant of the jacobian',
     &           ' matrix for gauss point ',i6,/7x,'of element ',
     &           i6,' is non-positive. current value: ',e12.5,/)
c
      end

c              subroutine to compute C = A * B
      subroutine mul33( span, felem, AA, BB, CC )
      implicit none
c
      include 'common.main'

c              global variables
      integer, intent(in)  :: span, felem
      real(8), intent(in)  :: AA(mxvl,ndim,ndim), BB(mxvl,ndim,ndim)
      real(8), intent(out) :: CC(mxvl,nstr) ! voigt notation
c
c              local variables
      integer :: i
      logical :: local_debug

      local_debug = .false.
c
c          calculate matrix multiplication 
c          AA and BB: full matrix
c          CC: xx, yy, zz, xy, yz, xz
c
!DIR$ LOOP COUNT MAX=128  
!DIR$ VECTOR ALIGNED
      do i = 1, span
      CC(i,1) = AA(i,1,1)*BB(i,1,1) + AA(i,1,2)*BB(i,2,1)
     &        + AA(i,1,3)*BB(i,3,1) ! xx
      CC(i,4) = AA(i,2,1)*BB(i,1,1) + AA(i,2,2)*BB(i,2,1)
     &        + AA(i,2,3)*BB(i,3,1) + AA(i,1,1)*BB(i,1,2)
     &        + AA(i,1,2)*BB(i,2,2) + AA(i,1,3)*BB(i,3,2) ! yx
      CC(i,6) = AA(i,3,1)*BB(i,1,1) + AA(i,3,2)*BB(i,2,1)
     &        + AA(i,3,3)*BB(i,3,1) + AA(i,1,1)*BB(i,1,3)
     &        + AA(i,1,2)*BB(i,2,3) + AA(i,1,3)*BB(i,3,3) ! zx
c     CC(i,1,2) = AA(i,1,1)*BB(i,1,2) + AA(i,1,2)*BB(i,2,2)
c    &         + AA(i,1,3)*BB(i,3,2) ! xy
      CC(i,2) = AA(i,2,1)*BB(i,1,2) + AA(i,2,2)*BB(i,2,2)
     &        + AA(i,2,3)*BB(i,3,2) ! yy
      CC(i,5) = AA(i,3,1)*BB(i,1,2) + AA(i,3,2)*BB(i,2,2)
     &        + AA(i,3,3)*BB(i,3,2) + AA(i,2,1)*BB(i,1,3)
     &        + AA(i,2,2)*BB(i,2,3) + AA(i,2,3)*BB(i,3,3) ! zy
c     CC(i,1,3) = AA(i,1,1)*BB(i,1,3) + AA(i,1,2)*BB(i,2,3)
c    &         + AA(i,1,3)*BB(i,3,3) ! xz
c     CC(i,2,3) = AA(i,2,1)*BB(i,1,3) + AA(i,2,2)*BB(i,2,3)
c    &         + AA(i,2,3)*BB(i,3,3) ! yz
      CC(i,3) = AA(i,3,1)*BB(i,1,3) + AA(i,3,2)*BB(i,2,3)
     &        + AA(i,3,3)*BB(i,3,3) ! zz
      end do

      if ( local_debug ) then
        write(out,9000)
        write(out,9001), AA(1,:,:)
        write(out,9001), BB(1,:,:)
        write(out,9002), CC(1,:)
      end if
c
 9000 format(3x,'Now checking matrix multiplication')
 9001 format(9e10.3)
 9002 format(6e10.3) 
c
      end subroutine 
