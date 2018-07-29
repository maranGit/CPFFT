c     ****************************************************************
c     *                                                              *
c     *                  subroutine  drive_eps_sig                   *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified : 4/5/2018 RM                *
c     *                                                              *
c     *      recovers all the strains, stresses                      *
c     *      for all the elements in the structure at state (n+1).   *
c     *                                                              *
c     *      Blocks are processed in parallel; using threads (OMP)   *
c     *                                                              *
c     ****************************************************************
c
      subroutine drive_eps_sig( step, iiter )
      implicit none
      include 'common.main'
      
c                         global
      integer :: step, iiter
      
c                         local
      integer :: blk, now_thread
      logical :: debug
      integer, external :: omp_get_thread_num
c
      debug = .false.
      call thyme( 2, 1 )
c
      call omp_set_dynamic( .false. )
c$OMP PARALLEL DO
c$OMP&         PRIVATE( blk, now_thread )
c$OMP&         SHARED( nelblk, iiter, step )
      do blk = 1, nelblk
        now_thread = omp_get_thread_num() + 1
        call do_nleps_block( blk, iiter, step )
      enddo
c$OMP END PARALLEL DO
c
      call thyme( 2, 2 )
c
      return
c
      end subroutine
c
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
      subroutine do_nleps_block( blk, iter, step )
      use fft, only: Fn1, Pn1, Fn, K4, matList, tstep, matprp
      use elem_block_data
      implicit none
      include 'common.main'
      include 'include_sig_up'
      
c             global
      integer, intent(in) :: blk, iter, step
c
c             locals
c
      integer :: ii, span, felem, currElem, currmat
      integer :: info_vector(4)
      logical :: local_debug, include_qbar
      
      integer :: ngp, gpn
      real(8), parameter :: zero = 0.0D0, one = 1.0D0, half = 0.5D0
      double precision :: gp_dtemps(mxvl)
      integer :: alloc_stat
      integer :: mat_type
      logical :: geo_non_flg
      integer :: now_thread
      integer, external :: OMP_GET_THREAD_NUM

c                       try using stack
      real(8) :: cep_blk_n1(mxvl,nstr,nstr), P_blk_n1(mxvl,nstrs)
      real(8) :: fn1inv(mxvl,ndim,ndim),rnh(mxvl,ndim,ndim)
      real(8) :: detF(mxvl), detFh(mxvl)
      real(8) :: fnh(mxvl,ndim,ndim),dfn(mxvl,ndim,ndim)
      real(8) :: fnhinv(mxvl,ndim,ndim),ddt(mxvl,nstr)
      real(8) :: uddt(mxvl,nstr),qnhalf(mxvl,nstr,nstr)
      real(8) :: qn1(mxvl,nstr,nstr),A_blk_n1(mxvl,nstrs*nstrs)
c     real(8), allocatable, dimension(:,:,:) :: rnh, fnh, dfn, fnhinv
c     real(8), allocatable, dimension(:) :: detF 
c     real(8), allocatable, dimension(:,:,:) :: fn1inv, cep_blk_n1
c     real(8), allocatable, dimension(:,:) :: ddt, uddt
c     real(8), allocatable, dimension(:,:,:) :: qnhalf, qn1
c     real(8), allocatable, dimension(:,:) :: P_blk_n1
c     real(8), allocatable, dimension(:,:) :: A_blk_n1
c
c             cauchy stress and rotation matrix @ n+1
      real(8) :: qtn1(mxvl,nstr,nstr), cs_blk_n1(mxvl,nstr)

c             allocate and initialization
c     allocate( detF(mxvl) )
c     allocate( fn1inv(mxvl,ndim,ndim) )
c     allocate( cep_blk_n1(mxvl,nstr,nstr) )
c     allocate( rnh(mxvl,ndim,ndim), fnh(mxvl,ndim,ndim), 
c    &          dfn(mxvl,ndim,ndim), fnhinv(mxvl,ndim,ndim) )
c             warp3d original allocate and initialization
c             nstr = 6; nstrs = 9;
c     allocate( ddt(mxvl,nstr), uddt(mxvl,nstr),
c    &          qnhalf(mxvl,nstr,nstr), qn1(mxvl,nstr,nstr) )
c     allocate( P_blk_n1(mxvl,nstrs) )
c     allocate( A_blk_n1(mxvl,nstrs*nstrs) )
      
      local_debug = .false.
      ddt    = zero
      uddt   = zero
      qnhalf = zero
      qn1    = zero
      detF   = zero
      rnh    = zero
      fnh    = zero
      dfn    = zero
      fnhinv = zero
      P_blk_n1 = zero
      cep_blk_n1 = zero
      fn1inv = zero
c
c     initialize local_work based on elblks(0,blk) and elblks(1,blk)
c
      ngp = fftngp
      gpn = 1
      geo_non_flg = .true.
      span = elblks(0, blk)
      felem = elblks(1, blk)
      currmat = matList(felem)
      mat_type = matprp(9,currmat)
      adaptive_flag = .false.
c
c
c                        initialize local work
c
      local_work%dt = tstep
      local_work%blk = blk
      local_work%span = span
      local_work%felem = felem
      local_work%num_int_points = ngp
      local_work%gpn = gpn
      local_work%step = step
      local_work%iter = iter ! global newton iteration (>1)
      local_work%geo_non_flg = geo_non_flg
      local_work%is_cohes_elem = .false.
      local_work%block_has_nonlocal_solids = .false.
      local_work%mat_type = mat_type
      local_work%matnum = currmat
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
      local_work%adaptive_flag = adaptive_flag
      local_work%eps_bbar = zero
      local_work%is_solid_matl = .true.
      local_work%is_umat = .false.
      local_work%umat_stress_type = 1
      local_work%is_crys_pls = (mat_type .eq. 10)
      local_work%temperatures = 297.0D0
      local_work%temperatures_ref = 297.0D0
c
c     allocate memory for local_work according to material model
c
c     call mm01_set_sizes( info_vector )
c     hist_size = info_vector(1)
c     cep_size = info_vector(2)
c     block_size = span * ngp * cep_size
c     local_work%hist_size_for_blk = hist_size
      call recstr_allocate( local_work )
      local_work%elem_type = 2 ! lsdisop element
c
c     grab material parameters, state variables and history
c     due(nodal displacement from n to n+1)
c     ue(nodal displacement from 0 to n)
c     ce_0, cd_n, cd_mid, cd_n1: nodal coordinate
c                 ***hard coded for now***
c
c                hard code material parameters
c     call mm01_hardCoded
      call rknstr_set_up_materials
c
c                grab global variables to local_work
c
      call dupstr_blocked( blk, span, felem, ngp, mat_type, 
     & geo_non_flg, step, iter, local_work )
c
c                grab global Fn and Fn1 to local_work
c
      do ii = 1, span
        currElem = felem + ii - 1
        local_work%fn(ii,1,1)= Fn(currElem,1)
        local_work%fn(ii,1,2)= Fn(currElem,2)
        local_work%fn(ii,1,3)= Fn(currElem,3)
        local_work%fn(ii,2,1)= Fn(currElem,4)
        local_work%fn(ii,2,2)= Fn(currElem,5)
        local_work%fn(ii,2,3)= Fn(currElem,6)
        local_work%fn(ii,3,1)= Fn(currElem,7)
        local_work%fn(ii,3,2)= Fn(currElem,8)
        local_work%fn(ii,3,3)= Fn(currElem,9)
      end do
      do ii = 1, span
        currElem = felem + ii - 1
        local_work%fn1(ii,1,1)= Fn1(currElem,1)
        local_work%fn1(ii,1,2)= Fn1(currElem,2)
        local_work%fn1(ii,1,3)= Fn1(currElem,3)
        local_work%fn1(ii,2,1)= Fn1(currElem,4)
        local_work%fn1(ii,2,2)= Fn1(currElem,5)
        local_work%fn1(ii,2,3)= Fn1(currElem,6)
        local_work%fn1(ii,3,1)= Fn1(currElem,7)
        local_work%fn1(ii,3,2)= Fn1(currElem,8)
        local_work%fn1(ii,3,3)= Fn1(currElem,9)
      end do
      do ii = 1, span
        fnh(ii,1,1) = half * ( local_work%fn(ii,1,1) 
     &              + local_work%fn1(ii,1,1) )
        fnh(ii,1,2) = half * ( local_work%fn(ii,1,2) 
     &              + local_work%fn1(ii,1,2) )
        fnh(ii,1,3) = half * ( local_work%fn(ii,1,3) 
     &              + local_work%fn1(ii,1,3) )
        fnh(ii,2,1) = half * ( local_work%fn(ii,2,1) 
     &              + local_work%fn1(ii,2,1) )
        fnh(ii,2,2) = half * ( local_work%fn(ii,2,2) 
     &              + local_work%fn1(ii,2,2) )
        fnh(ii,2,3) = half * ( local_work%fn(ii,2,3) 
     &              + local_work%fn1(ii,2,3) )
        fnh(ii,3,1) = half * ( local_work%fn(ii,3,1) 
     &              + local_work%fn1(ii,3,1) )
        fnh(ii,3,2) = half * ( local_work%fn(ii,3,2) 
     &              + local_work%fn1(ii,3,2) )
        fnh(ii,3,3) = half * ( local_work%fn(ii,3,3) 
     &              + local_work%fn1(ii,3,3) )
      enddo
      dfn = local_work%fn1 - local_work%fn
c
c              compute rotation tensor
c
      call rtcmp1( span, fnh, rnh )
      call rtcmp1( span, local_work%fn1,
     &             local_work%rot_blk_n1(1,1,gpn) )
c
c              compute ddt
c
      call inv33( span, felem, fnh, fnhinv, detFh )
      call mul33( span, felem, dfn, fnhinv, ddt, out )
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
c
c           For CP model, calculate the gradient of the elastic
c           rotations at the element level by linear curve fit.
c
c           For linear models this will just be based on the plastic rotations
c           which may or may not be a realistic assumption
c
      if( local_work%mat_type .eq. 10 ) call rknstr_finish_cp
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
c
c     pull back cs_blk_n1(Cauchy stress) to 1st PK stress
c     P = J*sigma*F^{-T}
c
      call inv33( span, felem, local_work%fn1, fn1inv, detF )
      call cs2p( span, felem, cs_blk_n1, fn1inv, detF, P_blk_n1 )
      if(local_debug .and. blk .eq. 1) then
        write(*,*) "current iteration: ", iter
      end if
c
c         rotate from ( d_urcs / d_uddt ) to ( Green-Naghdi / D )
c     (1) extract stiffness from mod_eleblocks
c     (2) push forward to current configuration
c     (3) store in cep_blk_n1(mxvl,nstr,nstr)
c     update: do not need to (2) push forward now. Because I can get
c             dP/dF from dt/dd directly. Now the dP/dF is exactly the
c             same compared with finite difference method.
c
      if( geo_non_flg .and. local_work%is_solid_matl ) then
        include_qbar = .false.
        call gptns1( local_work, cep_blk_n1, qtn1 )
      end if
c
c                 pull back cep_blk_n1 to dP/dF
c
      call cep2A( local_work, cep_blk_n1, rnh, detF, 
     &            detFh, fnhinv, fn1inv, A_blk_n1 )
c
c     update global P and tangent stiffness
c
      do ii = 1, span
        currElem = felem + ii - 1
        Pn1(currElem, 1:9) = P_blk_n1(ii, 1:9)
        K4(currElem, 1:81) = A_blk_n1(ii, 1:81)
      end do
c
c     get cut-step flag and update plastic work
c
c     do ii = 1, span
c       currElem = felem + ii - 1
c       call constitutive(Fn1(currElem,:), matList(currElem), 
c    &                     Pn1(currElem,:), K4(currElem,:))
c     enddo
c
c 100 deallocate( rnh, fnh, dfn, fnhinv, fn1inv )
c     deallocate( ddt, uddt, qnhalf, qn1 )
c     deallocate( P_blk_n1 )
c     deallocate( A_blk_n1 )
      call recstr_deallocate( local_work )
      return
      contains
c     ========
c     ****************************************************************
c     *                                                              *
c     *                   subroutine rknstr_set_up_materials         *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/26/2017 rhd              *
c     *                                                              *
c     *     make calls to the specific material model for block      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rknstr_set_up_materials
      implicit none
c
      logical :: adaptive
      integer :: iout
c
      adaptive = adaptive_flag .and. step .gt. 1
      iout = local_work%iout
c
      select case ( mat_type )
c
        case ( 1 )
          call setup_mm01_rknstr( span, adaptive, local_work )
c         call setup_segmental(  span, props, lprops, iprops,
c    &                           local_work )
c
c       case ( 2 )
c         call setup_mm02_rknstr( span, props, lprops, iprops,
c    &                            local_work )
c
c       case ( 3 )
c         call setup_mm03_rknstr( span, props, lprops, iprops,
c    &                            adaptive, local_work )
c         call setup_segmental(  span, props, lprops, iprops,
c    &                           local_work )
c
c       case ( 4 )
c         call mm04_init( iout, span, felem, props, lprops, iprops,
c    &                    local_work%cohes_type,
c    &                    local_work%intf_prp_block,
c    &                    matprp(1,local_work%matnum),
c    &                    local_work%cohes_rot_block )
c
c       case ( 5 )
c         call setup_mm05_rknstr( span, props, lprops, iprops,
c    &                            adaptive, local_work )
c         call setup_segmental(  span, props, lprops, iprops,
c    &                           local_work )
c
c       case ( 6 )
c         call setup_mm06_rknstr( span, props, lprops, iprops,
c    &                            adaptive, local_work )
c
c       case ( 7 )
c         call setup_mm07_rknstr( span, props, lprops, iprops,
c    &                            adaptive, local_work )
c
c       case ( 8 )
c         call setup_umat_rknstr( span, props, lprops, iprops,
c    &                            adaptive, local_work )
c
        case ( 10 )
          call setup_mm10_rknstr( span, adaptive, local_work )
c
        case default
          write(iout,*) '>>> invalid material model number'
          write(iout,*) '    in rknstr. abort execution'
          call die_abort
          stop
c
      end select
c
c           If this material is going to have interface damage applied to it,
c           call the setup function for mm11
c
c     if( local_work%is_inter_dmg ) then
c       call setup_mm11_rknstr( span, props, lprops, iprops,
c    &                            adaptive, local_work )
c     end if
c
      return
      end subroutine rknstr_set_up_materials
c
c     ****************************************************************
c     *                                                              *
c     *                   subroutine rknstr_finish_cp                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/26/2017 rhd              *
c     *                                                              *
c     *     make calls to the specific material model for block      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rknstr_finish_cp
      implicit none
c
      integer :: sh, eh, i, iblkrow, rs, re
c
c              calculate the gradient of the elastic
c              rotations at the element level by linear curve fit.
c
c     sh  = indexes_common(2,1) ! first index of grad_fe
c     eh  = indexes_common(2,2) ! last index of grad_fe
c
c     do i = 1, local_work%span
c       iblkrow = i
c       if( local_work%ncrystals(i) .gt. 1 )  then
c            local_work%elem_hist1(i,sh:eh,1:ngp) = zero
c       else
c        rs = index_crys_hist(1,3,1) ! first index of Rp, 1st crystal
c        re = index_crys_hist(1,3,2) ! last index of Rp, 1st crystal
c        call mm10_calc_grads( ngp, elem_type, order, geonl,
c    &       local_work%rot_blk_n1(1,1,1),
c    &       local_work%jac(1,1,1),
c    &       local_work%elem_hist(1,rs,1),
c    &       local_work%elem_hist1(1,sh,1), iout, iblkrow, span,
c    &       local_work%hist_size_for_blk )
c       end if
c     end do
c
      return
      end subroutine rknstr_finish_cp
      end subroutine ! do_nleps_block
c
c
c     ****************************************************************
c     *                                                              *
c     *                   subroutine setup_mm01_rknstr               *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified : 5/17/2018 RM               *
c     *                                                              *
c     *     set up material model #1 (bilinear mises) for stress     *
c     *     updating                                                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine setup_mm01_rknstr( span, adaptive, local_work )
      use fft, only: matprp, lmtprp, imatprp, dmatprp, smatprp
      implicit none
      include 'param_def'
      include 'include_sig_up'
      
c                     global
      integer :: span
      logical :: adaptive

c                     local
      integer :: currmat, i
      real(8) :: ym, nu, beta, tan_e, yld, alpha, rho
      
      currmat = local_work%matnum
      
      do i = 1, span
        ym                       = matprp(1,currmat)
        nu                       = matprp(2,currmat)
        yld                      = matprp(5,currmat)
        tan_e                    = matprp(4,currmat)
        rho                      = matprp(7,currmat)
        alpha                    = matprp(6,currmat)
        beta                     = matprp(3,currmat)
        local_work%e_vec(i)      = ym
        local_work%nu_vec(i)     = nu
        local_work%beta_vec(i)   = beta
        local_work%tan_e_vec(i)  =  tan_e
        local_work%sigyld_vec(i) = yld
        local_work%h_vec(i)      = tan_e*ym/(ym - tan_e) 
        local_work%e_vec_n(i)    = ym
        local_work%nu_vec_n(i)   = nu
      end do

      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                   subroutine setup_mm10_rknstr               *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 1/9/2017 rhd               *
c     *                                                              *
c     *     set up material model #10 (crystal plasticity)           *
c     *     for stress updating: values constant across all g. pts.  *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine setup_mm10_rknstr( span, adaptive, local_work )
c     use segmental_curves
      use fft, only : matprp, lmtprp, imatprp, dmatprp, smatprp
      use crystal_data, only : c_array, angle_input, crystal_input,
     &                              data_offset
c
      implicit none
      include 'param_def'
c
c                    parameter declarations
c
c     integer :: iprops(mxelpr,mxvl), span
c     real    :: props(mxelpr,mxvl)
c     logical :: lprops(mxelpr,mxvl), adaptive
      integer :: span
      logical :: adaptive
      include 'include_sig_up'
c
c                    local
c
      integer :: i, out, matnum, ctotal, c, cnum, s, elnum, osn
      double precision :: angles(3), bs(3), ns(3), temp66(6,6),
     &                    rot(3,3), trans_rot(3,3),
     &                    temp_vec3(3),
     &                    A(3,3), trans_A(3,3), A_symm(3,3),
     &                    A_asymm(3,3),
     &                    Rstiff(6,6), trans_Rstiff(6,6)
      double precision, parameter :: half = 0.5d00, zero = 0.0d00
      character :: aconv*5, atype*7
c
      ctotal = 0
      out = local_work%iout
      matnum = local_work%matnum
c
      do i = 1, span
c
       call setup_mm10_rknstr_a( 1 )  ! selected props -> local_work
c
c              will need to (in the near future) extract crystal
c              and orientation information.  Also possibly change
c              to allow for variable number of crystals in each
c              element block.
c
       elnum = local_work%felem+i-1
c
       do c = 1, local_work%ncrystals(i)
c
         call setup_mm10_rknstr_a( 2 ) ! get data for this crystal
         call setup_mm10_rknstr_a( 3 ) ! put props into local_work
         call setup_mm10_rknstr_a( 4 ) ! get crystal->reference rot
c
c                     set up and rotate our orientation tensors
c
         rot       = local_work%c_props(i,c)%rotation_g
         call set_up_mm10_rknstr_c( trans_rot, rot, 3 )
         call setup_mm10_rknstr_a( 5 ) ! get crystal->reference rots
c
c                      rotate forward our stiffness tensor
c
         call mm10_RT2RVE( trans_rot, Rstiff ) ! makes Rstiff
         call set_up_mm10_rknstr_c( trans_Rstiff, Rstiff, 6)
         temp66 = matmul( local_work%c_props(i,c)%init_elast_stiff,
     &                       trans_Rstiff )
         local_work%c_props(i,c)%init_elast_stiff = matmul(
     &               Rstiff, temp66 )   ! 6x6 * 6xx
       end do ! on crystals
c
           ctotal = ctotal + local_work%ncrystals(i)
c
      end do ! on span

c
c                   determine if material model can call for a
c                   reduction in the adaptive step size
c
      local_work%allow_cut = adaptive .and. lmtprp(22,matnum)
c
c
      return

 9501 format(/,1x,'>>>> Not implemented yet in rknstr setup!'/)
 9503 format(/,1x,'>>>> System error: unexpected angle type in rknstr!',
     &            ' Aborting.'/)
 9504 format(/,1x,'>>>> System error: unexpected angle conv in rknstr!',
     &            ' Aborting.'/)

      contains
c     ========
      subroutine set_up_mm10_rknstr_b( a, b, c )
      implicit none
      double precision :: a(3), b(3,3), c(3)
c
c                     a = [b] * c
c
      a(1) = b(1,1)*c(1) + b(1,2)*c(2) + b(1,3)*c(3)
      a(2) = b(2,1)*c(1) + b(2,2)*c(2) + b(2,3)*c(3)
      a(3) = b(3,1)*c(1) + b(3,2)*c(2) + b(3,3)*c(3)
c
      return
      end subroutine set_up_mm10_rknstr_b

      subroutine set_up_mm10_rknstr_c( a, b, n )
      implicit none
      integer :: n, i, j
      double precision :: a(n,n), b(n,n)
c
c                     a = trans[ b ]
c
      if( n .eq. 3 ) then
         do i = 1, 3
!DIR$ VECTOR ALIGNED
            do j = 1, 3
               a(i,j) = b(j,i)
            end do
         end do
         return
      end if
c
       do i = 1, 6
!DIR$ VECTOR ALIGNED
         do j = 1, 6
            a(i,j) = b(j,i)
         end do
      end do
c
      return
      end subroutine set_up_mm10_rknstr_c

      subroutine set_up_mm10_rknstr_d( trans_A, A, bs, ns )
      implicit none
c
      integer :: i, j
      double precision :: trans_A(3,3), A(3,3), bs(3), ns(3)
c
      double precision :: spreadbs(3,3), spreadns(3,3)
c
      spreadbs(1,1) = bs(1)
      spreadbs(2,1) = bs(2)
      spreadbs(3,1) = bs(3)
      spreadbs(1,2) = bs(1)
      spreadbs(2,2) = bs(2)
      spreadbs(3,2) = bs(3)
      spreadbs(1,3) = bs(1)
      spreadbs(2,3) = bs(2)
      spreadbs(3,3) = bs(3)
c
      spreadns(1,1) = ns(1)
      spreadns(2,1) = ns(1)
      spreadns(3,1) = ns(1)
      spreadns(1,2) = ns(2)
      spreadns(2,2) = ns(2)
      spreadns(3,2) = ns(2)
      spreadns(1,3) = ns(3)
      spreadns(2,3) = ns(3)
      spreadns(3,3) = ns(3)
c
      do i = 1, 3
!DIR$ VECTOR ALIGNED
        do j = 1, 3
         A(i,j) = spreadbs(i,j) * spreadns(i,j)
         trans_A(j,i) = A(i,j)
        end do
      end do
c
      return
c
      end  subroutine set_up_mm10_rknstr_d




      subroutine setup_mm10_rknstr_a( dowhat )
      implicit none
c
      integer :: dowhat

      select case( dowhat )

      case( 1 )
c
        local_work%alpha_vec(i,1) = matprp(26,matnum)
        local_work%alpha_vec(i,2) = matprp(27,matnum)
        local_work%alpha_vec(i,3) = matprp(28,matnum)
        local_work%alpha_vec(i,4) = matprp(29,matnum)
        local_work%alpha_vec(i,5) = matprp(30,matnum)
        local_work%alpha_vec(i,6) = matprp(31,matnum)
c
        local_work%alpha_vec_n(i,1) = matprp(26,matnum)
        local_work%alpha_vec_n(i,2) = matprp(27,matnum)
        local_work%alpha_vec_n(i,3) = matprp(28,matnum)
        local_work%alpha_vec_n(i,4) = matprp(29,matnum)
        local_work%alpha_vec_n(i,5) = matprp(30,matnum)
        local_work%alpha_vec_n(i,6) = matprp(31,matnum)
c
        local_work%debug_flag(i) = lmtprp(13,matnum)
        local_work%local_tol(i) = dmatprp(100,matnum)
c
c              may eventually change this to allow for
c              different # of crystals in block
c
        local_work%ncrystals(i) = imatprp(101,matnum)
c
        local_work%angle_convention(i) = imatprp(102,matnum)
        local_work%angle_type(i) = imatprp(103,matnum)
c
      case( 2 )
c
c              get the local crystal number
c
        if( imatprp(104,matnum) .eq. 1 ) then
              cnum = imatprp(105,matnum)
        elseif( imatprp(104,matnum) .eq. 2 ) then
              osn = data_offset(elnum)
              cnum = crystal_input(osn,c)
c
c             couldn't do this earlier, so check here
c
              if( (cnum .gt. max_crystals) .or.
     &              (cnum .lt. 0) ) then
               write (out,'("Crystal ", i3, " not valid")')
     &              cnum
                    call die_gracefully
              end if
        else
              write(out,9502)
              call die_gracefully
        end if
c
c             get the local orientation
c
        if (imatprp(107,matnum ) .eq. 1) then
             angles(1) = dmatprp(108,matnum)
             angles(2) = dmatprp(109,matnum)
             angles(3) = dmatprp(110,matnum)
        elseif( imatprp(107,matnum) .eq. 2 ) then
             osn = data_offset(elnum)
             angles(1:3) = angle_input(osn,c,1:3)
        else
             write (out,9502)
             call die_gracefully
        end if
c
      case( 3 )
c
c             we have the properties, we just need to extract
c             into our local structure
c
        local_work%c_props(i,c)%init_elast_stiff =
     &        c_array(cnum)%elast_stiff
        local_work%c_props(i,c)%init_angles = angles
        local_work%c_props(i,c)%nslip = c_array(cnum)%nslip
        local_work%c_props(i,c)%rateN = c_array(cnum)%harden_n
        local_work%c_props(i,c)%tauHat_y = c_array(cnum)%tau_hat_y
        local_work%c_props(i,c)%Go_y = c_array(cnum)%g_o_y
        local_work%c_props(i,c)%tauHat_v = c_array(cnum)%tau_hat_v
        local_work%c_props(i,c)%Go_v = c_array(cnum)%g_o_v
        local_work%c_props(i,c)%tau_a = c_array(cnum)%tau_a
        local_work%c_props(i,c)%burgers = c_array(cnum)%b
        local_work%c_props(i,c)%p_v = c_array(cnum)%p_v
        local_work%c_props(i,c)%q_v = c_array(cnum)%q_v
        local_work%c_props(i,c)%p_y = c_array(cnum)%p_y
        local_work%c_props(i,c)%q_y = c_array(cnum)%q_y
        local_work%c_props(i,c)%boltzman = c_array(cnum)%boltz
        local_work%c_props(i,c)%theta_o = c_array(cnum)%theta_o
        local_work%c_props(i,c)%eps_dot_o_v = c_array(cnum)%eps_dot_o_v
        local_work%c_props(i,c)%eps_dot_o_y = c_array(cnum)%eps_dot_o_y
        local_work%c_props(i,c)%mu_o = c_array(cnum)%mu_o
        local_work%c_props(i,c)%D_o  = c_array(cnum)%D_o
c
        local_work%c_props(i,c)%t_o = c_array(cnum)%t_o
        local_work%c_props(i,c)%tau_a = c_array(cnum)%tau_a
        local_work%c_props(i,c)%k_o = c_array(cnum)%k_o
        local_work%c_props(i,c)%h_type = c_array(cnum)%h_type
        local_work%c_props(i,c)%s_type = c_array(cnum)%slip_type
        local_work%c_props(i,c)%cnum = cnum
        local_work%c_props(i,c)%num_hard = c_array(cnum)%num_hard
        local_work%c_props(i,c)%tang_calc = c_array(cnum)%tang_calc
        local_work%c_props(i,c)%u1 = c_array(cnum)%u1
        local_work%c_props(i,c)%u2 = c_array(cnum)%u2
        local_work%c_props(i,c)%u3 = c_array(cnum)%u3
        local_work%c_props(i,c)%u4 = c_array(cnum)%u4
        local_work%c_props(i,c)%u5 = c_array(cnum)%u5
        local_work%c_props(i,c)%u6 = c_array(cnum)%u6
        local_work%c_props(i,c)%u7 = c_array(cnum)%u7
        local_work%c_props(i,c)%u8 = c_array(cnum)%u8
        local_work%c_props(i,c)%u9 = c_array(cnum)%u9
        local_work%c_props(i,c)%u10 = c_array(cnum)%u10
        local_work%c_props(i,c)%tau_y = c_array(cnum)%tau_y
        local_work%c_props(i,c)%tau_v = c_array(cnum)%tau_v
        local_work%c_props(i,c)%voche_m = c_array(cnum)%voche_m
        local_work%c_props(i,c)%iD_v = c_array(cnum)%iD_v
c
        local_work%c_props(i,c)%cp_001 = c_array(cnum)%cp_001
        local_work%c_props(i,c)%cp_002 = c_array(cnum)%cp_002
        local_work%c_props(i,c)%cp_003 = c_array(cnum)%cp_003
        local_work%c_props(i,c)%cp_004 = c_array(cnum)%cp_004
        local_work%c_props(i,c)%cp_005 = c_array(cnum)%cp_005
        local_work%c_props(i,c)%cp_006 = c_array(cnum)%cp_006
        local_work%c_props(i,c)%cp_007 = c_array(cnum)%cp_007
        local_work%c_props(i,c)%cp_008 = c_array(cnum)%cp_008
        local_work%c_props(i,c)%cp_009 = c_array(cnum)%cp_009
        local_work%c_props(i,c)%cp_010 = c_array(cnum)%cp_010
        local_work%c_props(i,c)%cp_011 = c_array(cnum)%cp_011
        local_work%c_props(i,c)%cp_012 = c_array(cnum)%cp_012
        local_work%c_props(i,c)%cp_013 = c_array(cnum)%cp_013
        local_work%c_props(i,c)%cp_014 = c_array(cnum)%cp_014
        local_work%c_props(i,c)%cp_015 = c_array(cnum)%cp_015
        local_work%c_props(i,c)%cp_016 = c_array(cnum)%cp_016
        local_work%c_props(i,c)%cp_017 = c_array(cnum)%cp_017
        local_work%c_props(i,c)%cp_018 = c_array(cnum)%cp_018
        local_work%c_props(i,c)%cp_019 = c_array(cnum)%cp_019
        local_work%c_props(i,c)%cp_020 = c_array(cnum)%cp_020
        local_work%c_props(i,c)%cp_021 = c_array(cnum)%cp_021
        local_work%c_props(i,c)%cp_022 = c_array(cnum)%cp_022
        local_work%c_props(i,c)%cp_023 = c_array(cnum)%cp_023
        local_work%c_props(i,c)%cp_024 = c_array(cnum)%cp_024
        local_work%c_props(i,c)%cp_025 = c_array(cnum)%cp_025
        local_work%c_props(i,c)%cp_026 = c_array(cnum)%cp_026
        local_work%c_props(i,c)%cp_027 = c_array(cnum)%cp_027
        local_work%c_props(i,c)%cp_028 = c_array(cnum)%cp_028
        local_work%c_props(i,c)%cp_029 = c_array(cnum)%cp_029
        local_work%c_props(i,c)%cp_030 = c_array(cnum)%cp_030
        local_work%c_props(i,c)%cp_031 = c_array(cnum)%cp_031
        local_work%c_props(i,c)%cp_032 = c_array(cnum)%cp_032
        local_work%c_props(i,c)%cp_033 = c_array(cnum)%cp_033
        local_work%c_props(i,c)%cp_034 = c_array(cnum)%cp_034
        local_work%c_props(i,c)%cp_035 = c_array(cnum)%cp_035
        local_work%c_props(i,c)%cp_036 = c_array(cnum)%cp_036
        local_work%c_props(i,c)%cp_037 = c_array(cnum)%cp_037
        local_work%c_props(i,c)%cp_038 = c_array(cnum)%cp_038
        local_work%c_props(i,c)%cp_039 = c_array(cnum)%cp_039
        local_work%c_props(i,c)%cp_040 = c_array(cnum)%cp_040
        local_work%c_props(i,c)%cp_041 = c_array(cnum)%cp_041
        local_work%c_props(i,c)%cp_042 = c_array(cnum)%cp_042
        local_work%c_props(i,c)%cp_043 = c_array(cnum)%cp_043
        local_work%c_props(i,c)%cp_044 = c_array(cnum)%cp_044
        local_work%c_props(i,c)%cp_045 = c_array(cnum)%cp_045
        local_work%c_props(i,c)%cp_046 = c_array(cnum)%cp_046
        local_work%c_props(i,c)%cp_047 = c_array(cnum)%cp_047
        local_work%c_props(i,c)%cp_048 = c_array(cnum)%cp_048
        local_work%c_props(i,c)%cp_049 = c_array(cnum)%cp_049
        local_work%c_props(i,c)%cp_050 = c_array(cnum)%cp_050
        local_work%c_props(i,c)%cp_051 = c_array(cnum)%cp_051
        local_work%c_props(i,c)%cp_052 = c_array(cnum)%cp_052
        local_work%c_props(i,c)%cp_053 = c_array(cnum)%cp_053
        local_work%c_props(i,c)%cp_054 = c_array(cnum)%cp_054
        local_work%c_props(i,c)%cp_055 = c_array(cnum)%cp_055
        local_work%c_props(i,c)%cp_056 = c_array(cnum)%cp_056
        local_work%c_props(i,c)%cp_057 = c_array(cnum)%cp_057
        local_work%c_props(i,c)%cp_058 = c_array(cnum)%cp_058
        local_work%c_props(i,c)%cp_059 = c_array(cnum)%cp_059
        local_work%c_props(i,c)%cp_060 = c_array(cnum)%cp_060
        local_work%c_props(i,c)%cp_061 = c_array(cnum)%cp_061
        local_work%c_props(i,c)%cp_062 = c_array(cnum)%cp_062
        local_work%c_props(i,c)%cp_063 = c_array(cnum)%cp_063
        local_work%c_props(i,c)%cp_064 = c_array(cnum)%cp_064
        local_work%c_props(i,c)%cp_065 = c_array(cnum)%cp_065
        local_work%c_props(i,c)%cp_066 = c_array(cnum)%cp_066
        local_work%c_props(i,c)%cp_067 = c_array(cnum)%cp_067
        local_work%c_props(i,c)%cp_068 = c_array(cnum)%cp_068
        local_work%c_props(i,c)%cp_069 = c_array(cnum)%cp_069
        local_work%c_props(i,c)%cp_070 = c_array(cnum)%cp_070
        local_work%c_props(i,c)%cp_071 = c_array(cnum)%cp_071
        local_work%c_props(i,c)%cp_072 = c_array(cnum)%cp_072
        local_work%c_props(i,c)%cp_073 = c_array(cnum)%cp_073
        local_work%c_props(i,c)%cp_074 = c_array(cnum)%cp_074
        local_work%c_props(i,c)%cp_075 = c_array(cnum)%cp_075
        local_work%c_props(i,c)%cp_076 = c_array(cnum)%cp_076
        local_work%c_props(i,c)%cp_077 = c_array(cnum)%cp_077
        local_work%c_props(i,c)%cp_078 = c_array(cnum)%cp_078
        local_work%c_props(i,c)%cp_079 = c_array(cnum)%cp_079
        local_work%c_props(i,c)%cp_080 = c_array(cnum)%cp_080
        local_work%c_props(i,c)%cp_081 = c_array(cnum)%cp_081
        local_work%c_props(i,c)%cp_082 = c_array(cnum)%cp_082
        local_work%c_props(i,c)%cp_083 = c_array(cnum)%cp_083
        local_work%c_props(i,c)%cp_084 = c_array(cnum)%cp_084
        local_work%c_props(i,c)%cp_085 = c_array(cnum)%cp_085
        local_work%c_props(i,c)%cp_086 = c_array(cnum)%cp_086
        local_work%c_props(i,c)%cp_087 = c_array(cnum)%cp_087
        local_work%c_props(i,c)%cp_088 = c_array(cnum)%cp_088
        local_work%c_props(i,c)%cp_089 = c_array(cnum)%cp_089
        local_work%c_props(i,c)%cp_090 = c_array(cnum)%cp_090
        local_work%c_props(i,c)%cp_091 = c_array(cnum)%cp_091
        local_work%c_props(i,c)%cp_092 = c_array(cnum)%cp_092
        local_work%c_props(i,c)%cp_093 = c_array(cnum)%cp_093
        local_work%c_props(i,c)%cp_094 = c_array(cnum)%cp_094
        local_work%c_props(i,c)%cp_095 = c_array(cnum)%cp_095
        local_work%c_props(i,c)%cp_096 = c_array(cnum)%cp_096
        local_work%c_props(i,c)%cp_097 = c_array(cnum)%cp_097
        local_work%c_props(i,c)%cp_098 = c_array(cnum)%cp_098
        local_work%c_props(i,c)%cp_099 = c_array(cnum)%cp_099
        local_work%c_props(i,c)%cp_100 = c_array(cnum)%cp_100
c
c          flags to control CP internal solvers
c
        local_work%c_props(i,c)%solver = c_array(cnum)%solver
        local_work%c_props(i,c)%strategy = c_array(cnum)%strategy
        local_work%c_props(i,c)%gpall = c_array(cnum)%gpall
        local_work%c_props(i,c)%gpp = c_array(cnum)%gpp
        local_work%c_props(i,c)%st_it(1:3) = c_array(cnum)%st_it(1:3)
        local_work%c_props(i,c)%method = c_array(cnum)%method
        local_work%c_props(i,c)%miter = c_array(cnum)%miter
        local_work%c_props(i,c)%atol = c_array(cnum)%atol
        local_work%c_props(i,c)%atol1 = c_array(cnum)%atol1
        local_work%c_props(i,c)%rtol = c_array(cnum)%rtol
        local_work%c_props(i,c)%rtol1 = c_array(cnum)%rtol1
        local_work%c_props(i,c)%xtol = c_array(cnum)%xtol
        local_work%c_props(i,c)%xtol1 = c_array(cnum)%xtol1
c
c          alternative model flag
c
        local_work%c_props(i,c)%alter_mode = c_array(cnum)%alter_mode
c
      case( 4 )
c
c          reference rotation
c
        if (local_work%angle_type(i) .eq. 1) then
              atype = "degrees"
        elseif (local_work%angle_type(i) .eq. 2) then
              atype = "radians"
        else
              write(out,9503)
              call die_gracefully
        end if
        if (local_work%angle_convention(i) .eq. 1) then
              aconv="kocks"
        else
              write(out,9504)
              call die_gracefully
        end if
        call mm10_rotation_matrix( local_work%c_props(i,c)%init_angles,
     &        aconv, atype, local_work%c_props(i,c)%rotation_g, out )
c
      case( 5 )
c
c          rotation on each slip system for this crystal
c
        do s = 1, local_work%c_props(i,c)%nslip
            temp_vec3 = c_array(cnum)%bi(s,1:3)
            call set_up_mm10_rknstr_b( bs, trans_rot, temp_vec3 )
            temp_vec3 =  c_array(cnum)%ni(s,1:3)
            call set_up_mm10_rknstr_b( ns, trans_rot, temp_vec3 )
            local_work%c_props(i,c)%ns(1:3,s) = ns
            call set_up_mm10_rknstr_d( trans_A, A, bs, ns )
            A_symm  = half * ( A + trans_A )
            A_asymm = half * ( A - trans_A )
            call mm10_ET2EV( A_symm,  local_work%c_props(i,c)%ms(1,s) )
            call mm10_WT2WV( A_asymm, local_work%c_props(i,c)%qs(1,s) )
        end do ! on slip systems

      end select
c
      return
 9502 format(/,1x,'>>>> System error: unexpected input type in rknstr!',
     &            ' Aborting.'/)
 9503 format(/,1x,'>>>> System error: unexpected angle type in rknstr!',
     &            ' Aborting.'/)
 9504 format(/,1x,'>>>> System error: unexpected angle conv in rknstr!',
     &            ' Aborting.'/)

c
      end subroutine setup_mm10_rknstr_a

c
      end subroutine setup_mm10_rknstr
c
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
      gama = zero
c
c           calculate the determinate of the jacobian matrix
c
      gpn = 1
      do i = 1, span
          j1(i)= jac(i,2,2)*jac(i,3,3)-jac(i,2,3)*jac(i,3,2)
          j2(i)= jac(i,2,1)*jac(i,3,3)-jac(i,2,3)*jac(i,3,1)
          j3(i)= jac(i,2,1)*jac(i,3,2)-jac(i,2,2)*jac(i,3,1)
          dj(i)= jac(i,1,1)*j1(i)-jac(i,1,2)*j2(i)+jac(i,1,3)*j3(i)
      end do
c
c           check to insure a positive determinate.
c
      do i = 1, span
       if( dj(i) .le. zero_check ) then
         write(*,9169) gpn,felem+i-1, dj(i)
       end if
      end do
c
c           calculate the inverse of the jacobian matrix
c
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
      subroutine mul33( span, felem, AA, BB, CC, iout )
      implicit none
c
      include 'param_def'

c              global variables
      integer, intent(in)  :: span, felem, iout
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
        write(iout,9000)
        write(iout,9001), AA(1,:,:)
        write(iout,9001), BB(1,:,:)
        write(iout,9002), CC(1,:)
      end if
      return
c
 9000 format(3x,'Now checking matrix multiplication')
 9001 format(9e10.3)
 9002 format(6e10.3) 
c
      end subroutine 
c     ****************************************************************
c     *                                                              *
c     *                      subroutine cs2p                         *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 3/7/2018                    *
c     *                                                              *
c     *                pull back Cauchy stress to                    *
c     *                      1st P-K stress                          *
c     *                    P = J*sigma*F^{-T}                        *
c     *                 P(11,12,13,21,22,23,31,32,33)                *
c     *                                                              *
c     ****************************************************************
c
      subroutine cs2p( span, felem, cs, finv, detF, P)
      implicit none
      include 'param_def'
c                          global
      integer :: span, felem
      real(8) :: cs(mxvl,*), finv(mxvl,ndim,*), detF(*)
      real(8) :: P(mxvl,*)
c                          local
      integer :: ii
c     
c                   cs_blk_n1 * transpose(fn1inv)
c       voigt notation: 1-11, 2-22, 3-33, 4-12, 5-23, 6-13
      do ii = 1, span
        ! P11 = cs(11)*finv(11) + cs(12)*finv(12) + cs(13)*finv(13)
        P(ii,1) = detF(ii) * ( cs(ii,1)*finv(ii,1,1) 
     &          + cs(ii,4)*finv(ii,1,2) + cs(ii,6)*finv(ii,1,3) )
        ! P12 = cs(11)*finv(21) + cs(12)*finv(22) + cs(13)*finv(23)
        P(ii,2) = detF(ii) * ( cs(ii,1)*finv(ii,2,1) 
     &          + cs(ii,4)*finv(ii,2,2) + cs(ii,6)*finv(ii,2,3) )
        ! P13 = cs(11)*finv(31) + cs(12)*finv(32) + cs(13)*finv(33)
        P(ii,3) = detF(ii) * ( cs(ii,1)*finv(ii,3,1) 
     &          + cs(ii,4)*finv(ii,3,2) + cs(ii,6)*finv(ii,3,3) )
        ! P21 = cs(21)*finv(11) + cs(22)*finv(12) + cs(23)*finv(13)
        P(ii,4) = detF(ii) * ( cs(ii,4)*finv(ii,1,1) 
     &          + cs(ii,2)*finv(ii,1,2) + cs(ii,5)*finv(ii,1,3) )
        ! P22 = cs(21)*finv(21) + cs(22)*finv(22) + cs(23)*finv(23)
        P(ii,5) = detF(ii) * ( cs(ii,4)*finv(ii,2,1) 
     &          + cs(ii,2)*finv(ii,2,2) + cs(ii,5)*finv(ii,2,3) )
        ! P23 = cs(21)*finv(31) + cs(22)*finv(32) + cs(23)*finv(33)
        P(ii,6) = detF(ii) * ( cs(ii,4)*finv(ii,3,1) 
     &          + cs(ii,2)*finv(ii,3,2) + cs(ii,5)*finv(ii,3,3) )
        ! P31 = cs(31)*finv(11) + cs(32)*finv(12) + cs(33)*finv(13)
        P(ii,7) = detF(ii) * ( cs(ii,6)*finv(ii,1,1) 
     &          + cs(ii,5)*finv(ii,1,2) + cs(ii,3)*finv(ii,1,3) )
        ! P32 = cs(31)*finv(21) + cs(32)*finv(22) + cs(33)*finv(23)
        P(ii,8) = detF(ii) * ( cs(ii,6)*finv(ii,2,1) 
     &          + cs(ii,5)*finv(ii,2,2) + cs(ii,3)*finv(ii,2,3) )
        ! P33 = cs(31)*finv(31) + cs(32)*finv(32) + cs(33)*finv(33)
        P(ii,9) = detF(ii) * ( cs(ii,6)*finv(ii,3,1) 
     &          + cs(ii,5)*finv(ii,3,2) + cs(ii,3)*finv(ii,3,3) )
      enddo
      return
      end subroutine
