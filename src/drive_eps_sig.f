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
c
      call omp_set_dynamic( .false. )
c$OMP PARALLEL DO ORDERED 
c$OMP&         PRIVATE( blk, now_thread )
c$OMP&         SHARED( nelblk, iiter, step )
      do blk = 1, nelblk
        now_thread = omp_get_thread_num() + 1
        call do_nleps_block( blk, iiter, step )
      enddo
c$OMP END PARALLEL DO

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
      use fft, only: Fn1, Pn1, Fn, K4, matList, mat_props, tstep
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
      real(8) :: detF(mxvl),fn1inv(mxvl,ndim,ndim),rnh(mxvl,ndim,ndim)
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
      mat_type = mat_props(currmat)%matnum
c
c
c                        initialize local work
c
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
      local_work%eps_bbar = zero
      local_work%is_solid_matl = .true.
      local_work%is_umat = .false.
      local_work%umat_stress_type = 1
      local_work%is_crys_pls = .false.
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
      call mm01_hardCoded
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
c              compute ddt( detF is a dummy argument )
c
      call inv33( span, felem, fnh, fnhinv, detF )
      call mul33( span, felem, dfn, fnhinv, ddt, out )
c
c              compute uddt
c
      call getrm1( span, qnhalf, rnh, 1 ) ! form qnhalf
      call qmply1( span, mxvl, nstr, qnhalf, ddt, uddt )
      if(blk.eq.1) then
        write(*,*) "uddt: "
        write(*,'(6D12.4)') uddt(1,:)
      endif
      call rstgp1_update_strains( span, mxvl, nstr, uddt,
     &                            local_work%ddtse(1,1,gpn) )
c
c     recover stress and stiffness
c
      call rstgp1( local_work, uddt )
      if(blk.eq.1) then
        write(*,*) "Now checking cep: "
        write(*,*) cep_blocks(blk)%vector(10)
      endif
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
c
      if( geo_non_flg .and. local_work%is_solid_matl ) then
        include_qbar = .false.
        if(local_debug) write(*,*) "Entering ctran1()"
        call ctran1( span, felem, blk, cep_blk_n1, qtn1, cs_blk_n1,
     &               include_qbar, detF, local_work%weights,
     &               local_work%is_umat, local_work%umat_stress_type,
     &               local_work%is_crys_pls, local_debug )
        if(local_debug) write(*,*) "Leaving ctran1()"
      end if
c
c     pull back cep_blocks to dP/dF
c     assume that Green-Naghdi rate is close to Lie derivative
c     see Simo & Hughes, chapter 7
c     input:  cep_blk_n1(mxvl,nstr,nstr)
c     output: A_blk_n1(mxvl,nstrs*nstrs)
      call cep2A( span, cs_blk_n1, cep_blk_n1, 
     &            fn1inv, detF, A_blk_n1, out)
      if(blk.eq.1) write(*,*) "A: ", A_blk_n1(1,11)
      if(local_debug .and. blk .eq. 2) then
        write(*,*) " Now checking P of 2nd block:"
        write(*,'(9D12.4)') P_blk_n1(1,:)
      endif
c
c     update global P and tangent stiffness
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
      real(8) :: ym, nu, beta, tan_e, yld, alpha_x, alpha_y, alpha_z
      
      gp_dtemps = zero ! temperature change over step
      
c                    grab material parameters from global to local
      do ii = 1, span
        currmat = matList(felem + ii - 1)
        ym = mat_props(currmat)%dmatprp(1)
        nu = mat_props(currmat)%dmatprp(2)
        yld = mat_props(currmat)%dmatprp(3)
        tan_e = mat_props(currmat)%dmatprp(4)
        alpha_x = mat_props(currmat)%dmatprp(6)
        alpha_y = mat_props(currmat)%dmatprp(7)
        alpha_z = mat_props(currmat)%dmatprp(8)
        beta = mat_props(currmat)%dmatprp(9)
        local_work%e_vec(ii) = ym 
        local_work%nu_vec(ii) = nu
        local_work%beta_vec(ii) = beta
        local_work%tan_e_vec(ii) =  tan_e
        local_work%sigyld_vec(ii) = yld
        local_work%h_vec(ii) = tan_e*ym/(ym - tan_e) 
        local_work%e_vec_n(ii) = ym
        local_work%nu_vec_n(ii) = nu
      end do

      local_work%rtse = zero ! "relative" trial elastic stress rate
        
c     other properties required by rstgp1
      local_work%block_energy = zero
      local_work%block_plastic_work = zero
      local_work%beta_fact = zero
      
c     other properties required by drive_01_udate
      local_work%dt = tstep
      local_work%temperatures = 297.0D0
      local_work%temperatures_ref = 297.0D0

      return
      end subroutine ! mm01_hardCoded
      end subroutine ! do_nleps_block
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
