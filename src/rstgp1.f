c     ****************************************************************
c     *                                                              *
c     *                      subroutine rstgp1                       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/2/2016 rhd               *
c     *                                                              *
c     *     supervise the computation of strains, stresses and       *
c     *     accompaning stress data at an integration point          *
c     *     for a block of similar elements that use the same        *
c     *     material model code                                      *
c     *                                                              *
c     *              ** geometric nonlinear version **               *
c     *                                                              *
c     ****************************************************************
c
c
c      subroutine rstgp1( props, lprops, iprops, local_work )
      subroutine rstgp1( local_work, uddt_out )
c     implicit integer (a-z)
      implicit none
c     use segmental_curves, only : max_seg_points
      integer, parameter :: max_seg_points=20
      integer :: span, felem, type, order, gpn, ngp, nnode
      integer :: ndof, step, iter, mat_type, iout, i, j, k
c
      include 'param_def'
c
c                      parameter declarations
c
c      real    :: props(mxelpr,*)   ! all 3 are same by read-only
c      logical :: lprops(mxelpr,*)
c      integer :: iprops(mxelpr,*)
      include 'include_sig_up'
c
c                       locally defined variables
c
      double precision ::
     &  internal_energy, beta_fact, eps_bbar, plastic_work, zero
      double precision,
     & allocatable :: ddt(:,:), uddt(:,:), qnhalf(:,:,:),
     &                qn1(:,:,:)
      double precision :: uddt_out(mxvl,nstr)
c
      logical :: cut_step,
     &           adaptive, geonl, bbar, material_cut_step,
     &           local_debug, adaptive_flag
c
      data local_debug, zero / .false., 0.0d00 / 
c
      internal_energy   = local_work%block_energy
      plastic_work      = local_work%block_plastic_work
      beta_fact         = local_work%beta_fact
      span              = local_work%span
      felem             = local_work%felem
      type              = local_work%elem_type
      order             = local_work%int_order
      gpn               = local_work%gpn
      ngp               = local_work%num_int_points
      nnode             = local_work%num_enodes
      ndof              = local_work%num_enode_dof
      geonl             = local_work%geo_non_flg
      step              = local_work%step
      iter              = local_work%iter
      bbar              = local_work%bbar_flg
      mat_type          = local_work%mat_type
      material_cut_step = local_work%material_cut_step
      adaptive_flag     = local_work%adaptive_flag
      eps_bbar          = local_work%eps_bbar
      adaptive          = adaptive_flag .and. step .gt. 1
      iout              = local_work%iout
      if( local_debug ) write(iout,*) '... in rstgp1'
c
c        allocate and zero. only span rows are used but
c        array operators (e.g uddt = ..) will operate on
c        full content and access uninitialized values.
c
      allocate( ddt(mxvl,nstr), uddt(mxvl,nstr),
     &          qnhalf(mxvl,nstr,nstr), qn1(mxvl,nstr,nstr) )
      ddt    = zero
      uddt   = uddt_out
      qnhalf = zero
      qn1    = zero
c
c        process cohesive elements separately
c
c      if( local_work%is_cohes_elem ) then
c        if ( local_debug ) write(*,*) '>> calling gtlsn2...'
c        call gtlsn2( span, nnode,
c     &               local_work%due, uddt,
c     &               local_work%ddtse(1,1,gpn),
c     &               local_work%cohes_rot_block,
c     &               local_work%shape(1,gpn),
c     &               local_work%elem_type,gpn, felem, iout )
c        go to 7000
c      end if
c
c        compute deformation gradients (F), perform polar decompositions,
c        F=RU, etc. to get strain increment over the step on the
c        unrotated configuration.
c
c        calculate the element displacements at n+1/2 and n+1
c
c      call rstgp1_a( ndof, nnode, span, local_work%ue,
c     &               local_work%due, local_work%uenh,
c     &               local_work%uen1, mxvl )
c
c        find the deformation gradients F=RU and stress/strain
c        transformation matrices at n+1/2, n+1. [R,n+1] is computed
c        and stored in block structure rot_blk_n1. [qnhalf], [qn1]
c        and {dfn1} are returned. dfn1 is det[F,n+1] for energy
c        integration. if we get a bad deformation jacobian (det <= 0,
c        terminate strain computations and request
c        immediate step size reduction if possible.
c
c        we store F at n and n+1 in local_work for use by WARP3D UMAT
c        and crystal plasticity if they need them
c        (for large displacement analysis).
c
c        qn1 is computed at present only for UMAT and crystal plasticity
c
c      error = 0
c      call gtmat1( qnhalf, qn1, error, local_work )
c      if( error .eq. 1 ) then
c         if( adaptive ) then
c            material_cut_step = .true.
c            local_work%material_cut_step = material_cut_step
c            go to 9999
c         else
c            write(iout,9820)
c            call abort_job
c         end if
c      end if
c
c         compute the deformation tensor increment ddt (often called D).
c         ddt returned in vector (6x1) form. we use the linear [B] matrix
c         evaluated at n+1/2 configuration * the step displacement increment.
c         (displacement increment also equals the velocity * dt)
c
c      call gtlsn1( span, nnode,
c     &             local_work%due, ddt,
c     &             local_work%gama_mid(1,1,1,gpn),
c     &             local_work%nxi(1,gpn), local_work%neta(1,gpn),
c     &             local_work%nzeta(1,gpn),
c     &             local_work%vol_block, bbar, eps_bbar,
c     &             local_work%b )
c
c         compute the unrotated increment (uddt) of the deformation tensor
c         (also called delta-d) (the "d" rate * dt ). uddt = [qnhalf] * ddt.
c         the stress update is driven by uddt.
c         add delta-d to acumulated (integral) over all steps. The unrotated
c         increments and total all refer to the fixed (global) coordinate
c         axes.
c
c      call qmply1( span, mxvl, nstr, qnhalf, ddt, uddt )
c      call rstgp1_update_strains( span, mxvl, nstr, uddt,
c     &                            local_work%ddtse(1,1,gpn) )
c
c------------------------------------------------------------------------
c
 7000 continue
      select case ( mat_type )
      case ( 1 )
c
c                vectorized mises plasticity model
c
c       call drive_01_update( gpn, props, lprops, iprops,
c     &                       local_work, uddt, iout )
       call drive_01_update( gpn, local_work, uddt, iout )
c
      case( 2 )
c
c                linear+power law deformation plasticity model.
c                not supported for finite strain solutions,
c
       write(iout,9000)
       call die_abort
       stop
c
      case( 3 )
c
c                general mises/gurson flow theory model.
c
c       call drive_03_update( gpn, props, lprops, iprops,
c     &                       local_work, uddt, iout )
c       call drive_03_update( gpn, local_work, uddt, iout )
c
      case( 4 )
c
c                linear and non-linear cohesive model
c
c       call drive_04_update( gpn, props, lprops, iprops,
c     &                       local_work, uddt, iout )
c       call drive_04_update( gpn, local_work, uddt, iout )
c
      case ( 5 )
c
c                cyclic plasticity model
c
c       call drive_05_update( gpn, props, lprops, iprops,
c     &                       local_work, uddt, iout )
c       call drive_05_update( gpn, local_work, uddt, iout )
      case ( 6 )
c
c                creep model
c
c       call drive_06_update( gpn, props, lprops, iprops,
c     &                       local_work, uddt, iout )
c       call drive_06_update( gpn, local_work, uddt, iout )
      case ( 7 )
c
c                mises + hydrogen effects
c
c       call drive_07_update( gpn, props, lprops, iprops,
c     &                       local_work, uddt, iout )
c       call drive_07_update( gpn, local_work, uddt, iout )
c
      case ( 8 )
c
c                general UMAT
c
c       call drive_umat_update( gpn, local_work, uddt, qn1, iout )
c
      case ( 9 )
c
c                ALCOA anisotropic plasticty
c
c       call drive_09_update( gpn, props, lprops, iprops,
c     &                       local_work, uddt, iout )
c       call drive_09_update( gpn, local_work, uddt, iout )
c
      case ( 10 )
c
c               CP model
c
c      call drive_10_update( gpn, props, lprops, iprops,
c     &                       local_work, uddt, iout)
c       call drive_10_update( gpn, local_work, uddt, iout )
c
      case default
        write(iout,*) '>>> invalid material model number'
        write(iout,*) '    in rstgp1'
        call die_abort
c
      end select
c
c              If we're actually an interface damage model, 
c              call the interface damage calculations
c
c      if( local_work%is_inter_dmg )
c     &       call drive_11_update(gpn, props, lprops, iprops,
c     &            local_work, uddt, iout)
c
c --------------------------------------------------------------------
c
c          calculate the internal energy and plastic work
c          y integrating the densities over the deformed volume of the.
c          elment. urcs_blk_n1(..,7) is really the current (total) energy
c          density per unit deformed volume - look above...
c          increment of plastic work density stored in plastic_work_incr
c
c      if( iter .ne. 0 ) then
c        call rstgp1_b( span, internal_energy, plastic_work,
c     &                 local_work%urcs_blk_n1(1,7,gpn),
c     &                 local_work%urcs_blk_n1(1,8,gpn),
c     &                 local_work%det_j(1,gpn), local_work%dfn1, 1 )
c        internal_energy = internal_energy * beta_fact *
c     &                    local_work%weights(gpn)
c        plastic_work    = plastic_work * beta_fact *
c     &                    local_work%weights(gpn)
c        local_work%block_energy       = internal_energy
c        local_work%block_plastic_work = plastic_work
c      end if
c
      if( local_debug ) then
        write(iout,*) '>> rstgp1 .. gauss point: ', gpn
        write (iout,9500) internal_energy,
     &                 plastic_work
        write(iout,9110)
        do i = 1, span
         write(iout,9100) i, (local_work%urcs_blk_n(i,k,gpn),k=1,7),
     &                 (local_work%urcs_blk_n1(i,k,gpn),k=1,7),
     &                 (ddt(i,k),k=1,6),
     &                 (uddt(i,k),k=1,6)
        end do
      end if
c
 9999 continue
      deallocate( ddt, uddt, qnhalf, qn1 )
c
      return
c
 9000 format('>>> Fatal Error: the nonlinear elastic material',
     &     /,'                 is not compatible with large',
     &     /,'                 displacement elements.',
     &     /,'                 job terminated....' )
 9100 format(i5,7f15.6,/,5x,7f15.6,/,5x,6f15.6,/,5x,6f15.6)
 9110 format(1x,'Elem    /',20('-'),
     &        ' unrot. Cauchy @ n, unrot. Cauchy @ n+1,',
     &        ' ddt, uddt', 20('-'),'/')
 9500 format('  Internal energy inside of (rstgp1)    = ',e16.6,
     &     /,'  Plasstic work inside of (rstgp1)      = ',e16.6)
 9820 format(///,
     &       '>> FATAL ERROR: strain computation routines requested',
     &     /,'                immediate step size reduction due to',
     &     /,'                invalid determinant of a deformation',
     &     /,'                jacobian. the user has not allowed step',
     &     /,'                size reductions. analysis terminated.',
     &     /// )
      end
c
c     ****************************************************************
c     *                                                              *
c     *                subroutine drive_01_update                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *              last modified : 12/21/2015 rhd                  *
c     *                                                              *
c     *     drives material model #1 (bilinear) to                   *
c     *     update stresses and history for all elements in the      *
c     *     block for gauss point gpn                                *
c     *                                                              *
c     ****************************************************************
c
c
c      subroutine drive_01_update( gpn, props, lprops, iprops,
c     &                            local_work, uddt_displ, iout )
      subroutine drive_01_update( gpn, local_work, uddt_displ, iout )
c      use segmental_curves, only : max_seg_points
      use elem_block_data, only : gbl_cep_blocks => cep_blocks
c      use main_data, only : extrapolated_du, non_zero_imposed_du 
c
      implicit none
      include 'param_def'

      integer, parameter :: max_seg_points=20
c
c                      parameter declarations
c
c      real ::    props(mxelpr,*)   ! all same but read only
c      logical :: lprops(mxelpr,*)
c      integer :: iprops(mxelpr,*)
      integer :: gpn, iout
      double precision ::  uddt_displ(mxvl,nstr)
      include 'include_sig_up'
c
c                       locally defined variables
c
      integer :: span, felem, type, order, ngp, nnode, ndof, step,
     &           iter, now_blk, mat_type, number_points, curve_set,
     &           hist_size_for_blk, curve_type, elem_type, i
c
      double precision ::
     &  dtime, gp_temps(mxvl), gp_rtemps(mxvl), gp_dtemps(mxvl),
     &  zero,  ddummy(1), gp_alpha, ymfgm, et, uddt_temps(mxvl,nstr),
     &  uddt(mxvl,nstr), cep(mxvl,6,6)
c
      logical :: geonl, local_debug, temperatures, segmental,
     &           temperatures_ref, fgm_enode_props
c
      data zero / 0.0d0 /
c
c           vectorized mises plasticity model with constant hardening
c           modulus. the model supports temperature dependence of
c           the elastic modulus, nu, hprime, and thermal
c           expansion alpha can vary. temperature dependent
c           properties enter through segmental curves.
c
      dtime             = local_work%dt
      span              = local_work%span
      felem             = local_work%felem
      type              = local_work%elem_type
      order             = local_work%int_order
      ngp               = local_work%num_int_points
      nnode             = local_work%num_enodes
      ndof              = local_work%num_enode_dof
      geonl             = local_work%geo_non_flg
      step              = local_work%step
      iter              = local_work%iter
      now_blk           = local_work%blk
      mat_type          = local_work%mat_type
      temperatures      = local_work%temperatures
      temperatures_ref  = local_work%temperatures_ref
      segmental         = local_work%segmental
      number_points     = local_work%number_points
      curve_set         = local_work%curve_set_number
      fgm_enode_props   = local_work%fgm_enode_props
      hist_size_for_blk = local_work%hist_size_for_blk
c
      local_debug       = .false. ! felem .eq. 1 .and. gpn .eq. 3
      if( local_debug ) then
        write(iout,9000) felem, gpn, span
        write(iout,9010) dtime, type, order, nnode, ndof, geonl, step, 
     &                   iter, now_blk, mat_type,
     &                   temperatures, temperatures_ref, segmental,
     &                   number_points, curve_set, 
     &                   fgm_enode_props, hist_size_for_blk 
      end if
c
c          determine if the material elastic and properties are
c          described by segmental curve(s) and if they are temperature
c          or strain rate dependent [=0 no dependence, =1 temperature
c          dependent, =2 strain-rate dependent (not used here)]
c
c      curve_type = -1
c      if( segmental ) call set_segmental_type( curve_set, curve_type,
c     &                                          local_work%eps_curve )
c
c          determine if some material properties are specified using
c          values at model nodes to define fgms. interpolate values
c          at the current gauss point (overwrite the constant values).
c          the material cannot be both segmental & fgm (already
c          checked). interpolated values are the same at n and n+1.
c          at present only e, nu, tan_e, sig_yld and
c          isotropic alpha can be specified as fgm properties
c          interpolated from nodal values. note we have to
c          recompute h-prime after the correct e and tan_e are found.
c
c      if( fgm_enode_props ) then
c
c          call  set_fgm_solid_props_for_block(
c     &            span, felem, type, gpn, nnode,
c     &            local_work%e_vec_n, local_work%shape(1,gpn),
c     &            local_work%enode_mat_props, 1,
c     &            local_work%fgm_flags(1,1) )
c          call  set_fgm_solid_props_for_block(
c     &            span, felem, elem_type, gpn, nnode,
c     &            local_work%nu_vec_n, local_work%shape(1,gpn),
c     &            local_work%enode_mat_props, 2,
c     &            local_work%fgm_flags(1,2) )
c          call  set_fgm_solid_props_for_block(
c     &            span, felem, elem_type, gpn, nnode,
c     &            local_work%alpha_vec_n(1,1), local_work%shape(1,gpn),
c     &            local_work%enode_mat_props, 3,
c     &            local_work%fgm_flags(1,3) )
c          call  set_fgm_solid_props_for_block(
c     &            span, felem, elem_type, gpn, nnode,
c     &            local_work%tan_e_vec(1), local_work%shape(1,gpn),
c     &            local_work%enode_mat_props, 6,
c     &            local_work%fgm_flags(1,6) )
c          call  set_fgm_solid_props_for_block(
c     &            span, felem, elem_type, gpn, nnode,
c     &            local_work%sigyld_vec(1), local_work%shape(1,gpn),
c     &            local_work%enode_mat_props, 7,
c     &            local_work%fgm_flags(1,7) )
c
c          do i = 1, span
c            local_work%e_vec(i)  = local_work%e_vec_n(i)
c            local_work%nu_vec(i) = local_work%nu_vec_n(i)
c            gp_alpha  = local_work%alpha_vec_n(i,1)
c            local_work%alpha_vec(i,1)   = gp_alpha
c            local_work%alpha_vec(i,2)   = gp_alpha
c            local_work%alpha_vec(i,3)   = gp_alpha
c            local_work%alpha_vec(i,4)   = zero
c            local_work%alpha_vec(i,5)   = zero
c            local_work%alpha_vec(i,6)   = zero
c            local_work%alpha_vec_n(i,2) = gp_alpha
c            local_work%alpha_vec_n(i,3) = gp_alpha
c            local_work%alpha_vec_n(i,4) = zero
c            local_work%alpha_vec_n(i,5) = zero
c            local_work%alpha_vec_n(i,6) = zero
c            ymfgm = local_work%e_vec_n(i)
c            et    = local_work%tan_e_vec(i)
c            local_work%h_vec(i) = (ymfgm*et)/(ymfgm-et)
c          end do
c
c          if( local_debug ) write(iout,9020)
c      end if
c
c           get increment of temperature at gauss point for elements
c           in the block, the temperature at end of step and the
c           reference temperature.
c
c      call gauss_pt_temps(
c     &        local_work%dtemps_node_blk, gpn, type, span, order,
c     &        nnode, gp_dtemps, local_work%temps_node_blk,
c     &        gp_temps, temperatures, local_work%temps_node_to_process,
c     &        temperatures_ref, local_work%temps_ref_node_blk,
c     &        gp_rtemps )
c      if( local_debug ) write(iout,9030)
c
c          get temperature dependent young's modulus, poisson's
c          ratio, uniaxial yield stress and constant plastic
c          modulus for the temperature at n and n+1. Some
c          properties are not for this model but need to be passed
c          to satisfy syntax of call.
c
c      if( curve_type .eq. 0 .or. curve_type .eq. 1 ) then
c          call set_up_segmental( span, gp_temps, local_work%e_vec,
c     &        local_work%nu_vec, local_work%alpha_vec,
c     &        local_work%e_vec_n, local_work%nu_vec_n,
c     &        local_work%alpha_vec_n,
c     &        local_work%gp_sig_0_vec,
c     &        local_work%gp_h_u_vec,
c     &        local_work%gp_beta_u_vec,
c     &        local_work%gp_delta_u_vec,
c     &        local_work%gp_sig_0_vec_n,
c     &        local_work%gp_h_u_vec_n,
c     &        local_work%gp_beta_u_vec_n,
c     &        local_work%gp_delta_u_vec_n,
c     &        gp_dtemps,
c     &        ddummy, gpn, mxvl )
c
c         call set_up_h_prime( span, local_work%h_vec,
c     &                        local_work%sigyld_vec, felem )
c         if( local_debug ) write(iout,9040)
c      end if
c
c          get the thermal strain increment (actually negative of increment)
c
      uddt_temps = zero
c      if ( temperatures ) then
c        call gp_temp_eps( span, uddt_temps,
c     &                    local_work%alpha_vec, gp_dtemps ,
c     &                    gp_temps, gp_rtemps,
c     &                    local_work%alpha_vec_n )
c        if( local_debug ) write(iout,9050)
c      end if
c    
c            uddt_displ - strain increment due to displacement increment
c            uddt_temps - (negative) of strain increment just due 
c                         to temperature change
c
      uddt = uddt_displ + uddt_temps
      cep  = zero
c      
      do i = 1, span
       if( local_work%killed_status_vec(i) ) uddt(i,1:nstr) = zero
      end do
c    
c            now standard update process. use nonlinear update and [D]
c            computation for iter = 0 and extrapolation or iter > 1
c            for iter = 0 and no extrapolation, use linear-elastic [D]
c            with props at n+1.

c      if( iter >= 1 .or. extrapolated_du ) then !nonlinear update
c
c     ****************************************************************
c     *                                                              *
c     *                 Ran modify this line                         *
c     *                 no need to call drive_01_update              *
c     *                                                              *
c     ****************************************************************
c     if( iter >= 1 ) then !nonlinear update
      if( .true. ) then 
       if( local_debug ) write(iout,9060)
c       call mm01( span, felem, gpn, step, iter, local_work%e_vec,
c     &           local_work%nu_vec, local_work%beta_vec,
c     &           local_work%h_vec, local_work%lnelas_vec,
c     &           local_work%sigyld_vec,
c     &           local_work%urcs_blk_n(1,1,gpn),
c     &           local_work%urcs_blk_n1(1,1,gpn),
c     &           uddt, local_work%elem_hist(1,1,gpn),
c     &           local_work%elem_hist1(1,1,gpn),
c     &           local_work%rtse(1,1,gpn), gp_dtemps,
c     &           local_work%e_vec_n, local_work%nu_vec_n  )
       if(local_debug .and. local_work%blk .eq. 2) then
         write(*,*) "ym_n1: ", local_work%e_vec(1:span)
         write(*,*) "nu_n1: ", local_work%nu_vec(1:span)
         write(*,*) "beta: ", local_work%beta_vec(1:span)
         write(*,*) "hprime_n1: ", local_work%h_vec(1:span)
         write(*,*) "yld_n1: ", local_work%sigyld_vec(1:span)
         write(*,*) "cgn", local_work%urcs_blk_n(1,:,1)
         write(*,*) "uddt_voigt", uddt(1,:)
         write(*,*) "currhist", local_work%elem_hist(1,:,1)
         write(*,*) "ym_n", local_work%e_vec_n(1:span)
         write(*,*) "nu_n", local_work%nu_vec_n(1:span)
       endif
       call mm01( span, felem, gpn, step, iter, local_work%e_vec,
     &           local_work%nu_vec, local_work%beta_vec,
     &           local_work%h_vec, local_work%lnelas_vec,
     &           local_work%sigyld_vec,
     &           local_work%urcs_blk_n(1,1,gpn),
     &           local_work%urcs_blk_n1(1,1,gpn),
     &           uddt, local_work%elem_hist(1,1,gpn),
     &           local_work%elem_hist1(1,1,gpn),
     &           local_work%rtse(1,1,gpn), gp_dtemps,
     &           local_work%e_vec_n, local_work%nu_vec_n,
     &           local_work%iout  )
       call cnst1( span, cep, local_work%rtse(1,1,gpn),
     &            local_work%nu_vec,
     &            local_work%e_vec, local_work%elem_hist1(1,2,gpn),
     &            local_work%elem_hist1(1,5,gpn), local_work%beta_vec,
     &            local_work%elem_hist1(1,1,gpn), 
     &            local_work%elem_hist1(1,4,gpn), felem, iout )
      else  ! linear-elastic update
        if( local_debug ) write(iout,9070)
        call drive_01_update_a
      end if
c
c          save the [D] matrices (lower-triangle)
c
      call rstgp1_store_cep( span, mxvl, gpn, 
     &         gbl_cep_blocks(now_blk)%vector, 21, cep )
      if( local_debug ) write(iout,9080)
c      
      return
c
 9000 format(1x,'.... debug mm01. felem, gpn, span: ',i7,i3,i3)
 9010 format(10x,'...dtime, type, order, nnode, ndof:',e14.6,4i5,
     &     /,10x,'...geonl, step, iter, now_blk, mat_type: ',l2,4i5,
     &     /,10x,'...temperatures, temperatures_ref: ',
     &               2l2,
     &     /,10x,'...segmental, number_points, curve_set: ',l2,i3,i3,
     &     /,10x,'...fgm_enode_props, hist_size_for_blk: ',
     &    l3,i4 ) 
 9610 format(' >> rate iterations to converge: ',i3 )
 9020 format(10x,'...fgm properties determined...')             
 9030 format(10x,'...temperatures computed at integration point...')             
 9040 format(10x,'...temperatures dependent properties computed...')             
 9050 format(10x,'...thermal strains computed...') 
 9060 format(10x,'...update stresses nonlinear procedure...' )            
 9070 format(10x,'...update stresses use linear [D]...' )            
 9080 format(10x,'...[D]s saved to global structure...')
c
      contains
c     ========      
c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_01_update_a                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/21/2015 rhd             *
c     *                                                              *
c     ****************************************************************
c
      subroutine drive_01_update_a
      implicit none
c
      integer :: k, m, i
      double precision :: one, two, e, nu, c1, c2, c3, c4
      data one, two / 1.0d00, 2.0d00 /
c
c              get linear-elastic [D] with potentially temperature
c              dependent properties
c      
      do i = 1, span
c         if( local_work%killed_status_vec(i) ) cycle
         e  = local_work%e_vec(i)
         nu = local_work%nu_vec(i)
         c1 = (e/((one+nu)*(one-two*nu)))
         c2 = (one-nu)*c1
         c3 = ((one-two*nu)/two)*c1
         c4 = nu*c1
         cep(i,1,1)= c2
         cep(i,2,2)= c2
         cep(i,3,3)= c2
         cep(i,4,4)= c3
         cep(i,5,5)= c3
         cep(i,6,6)= c3
         cep(i,1,2)= c4
         cep(i,1,3)= c4
         cep(i,2,1)= c4
         cep(i,3,1)= c4
         cep(i,2,3)= c4
         cep(i,3,2)= c4
      end do
c
c              stresses at n+1 using linear-elastic [D]
c 
       call drive_01_update_b( span, mxvl, uddt, cep, 
     &                        local_work%urcs_blk_n(1,1,gpn),
     &                        local_work%urcs_blk_n1(1,1,gpn),
     &                        local_work%killed_status_vec )
c     
      return
      end subroutine drive_01_update_a
c
      end subroutine drive_01_update
c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_01_update_b                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/21/2015 rhd             *
c     *                                                              *
c     *     support routine for mm01 material driver.                *
c     *     should be inlined                                        *
c     *                                                              *
c     ****************************************************************

      subroutine drive_01_update_b( span, mxvl, uddt,
     &                              local_cep, stress_n, stress_np1,
     &                              killed_status )                
      implicit none
c      
      integer :: span, mxvl
      logical :: killed_status(*)
      double precision ::
     &  local_cep(mxvl,6,6), stress_n(mxvl,6), stress_np1(mxvl,6),
     &  uddt(mxvl,6), zero
      data zero / 0.0d00 /
c
      integer i, k, m
c
c              for each element in block, update stresses by
c              [D-elastic] * uddt. uddt contains thermal increment +
c              increment from imposed nodal displacements
c
      stress_np1 = stress_n
c      
      do k = 1, 6
       do m = 1, 6
         do i = 1, span
           stress_np1(i,k) = stress_np1(i,k) + 
     &                       local_cep(i,m,k) * uddt(i,m)
         end do
       end do
      end do   
c
      do i = 1, span
        if( killed_status(i) ) stress_np1(i,1:6) = zero
      end do
c         
      return
      end

c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine material_model_info               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 04/28/12                   *
c     *                                                              *
c     *     call the material model specific set up routine          *
c     *     to get a vector of various data sizes, parameters        *
c     *     return specific requested value to the calling routine   *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine material_model_info( element_no, block_no, info_type,
     &                                 value )
      use fft, only: matList, mat_props
      implicit none
      include 'common.main'
c                      global data
      integer :: element_no, block_no, info_type, value
c
c                      local data
c
      integer :: local_element_no, mat_type
      integer :: info_vector(10)
      integer :: inter_mat
      logical :: is_inter_dmg
      integer :: currmat  
c
      is_inter_dmg = .false.
c
c                      get the material model type associated with
c                      the element number or the elements in the
c                      block number. calling routine cannot set both
c                      element and block numbers.
c
c     info_vector:
c         1        number of history values per integration
c                  point. Abaqus calles these "statev". Values
c                  double or single precsion based on hardware.
c
c         2        number of values in the symmetric part of the
c                  [D] for each integration point. for solid
c                  elements this is 21, for cohesive elements this 6.
c
c         3        = 0, the material model returns "unrotated"
c                       Cauchy stresses at n+1
c                  = 1, the material model returns the standard
c                       Cauchy stresses at n+1
c
c         4        number of state variables per point to be output
c                  when user requests this type of results
c
      info_vector(1) = 0
      info_vector(2) = 21
      info_vector(3) = 0
      info_vector(4) = 0
c
      if( element_no .gt. 0 .and. block_no .gt. 0 ) then
         write(out,9000) 1
         call die_gracefully
      end if
c
      local_element_no = element_no
      if( element_no .le. 0 ) then
        if( block_no .gt. nelblk .or. block_no .le. 0 ) then
          write(out,9000) 2
          call die_gracefully
        else
          local_element_no = elblks(1,block_no)
        end if
      end if
c
      if( local_element_no .le. 0 .or.
     &    local_element_no .gt. noelem ) then
             write(out,9000) 3
             call die_gracefully
      end if
c
c      mat_type = iprops(25,local_element_no)
      currmat = matList(element_no)
      mat_type = mat_props(currmat)%matnum
c
c              See if we're actually a interface-damaged model
c
c      if( iprops(42, local_element_no) .ne. -1 ) then
c        inter_mat = iprops(42,local_element_no)
c        is_inter_dmg = .true.
c      end if
c
      select case( mat_type )
      case( 1 )
        call mm01_set_sizes( info_vector )
c     case( 2 )
c       call mm02_set_sizes( info_vector )
c     case( 3 )
c       call mm03_set_sizes( info_vector )
c     case( 4 )
c       call mm04_set_sizes( info_vector )
c     case( 5 )
c       call mm05_set_sizes( info_vector )
c     case( 6 )
c       call mm06_set_sizes( info_vector )
c     case( 7 )
c       call mm07_set_sizes( info_vector )
c     case( 8 )
c       call umat_set_features( info_vector )
c     case( 9 )
c       call mm09_set_sizes( info_vector )
c     case(10 )
c       call mm10_set_sizes_special( info_vector, local_element_no )
      case default
        write(out,9000) 4
        call die_gracefully
      end select
c
c              change history length if we are actually an 
c              interface damaged material
c
      if( is_inter_dmg )
     &   call mm11_set_sizes_special(inter_mat,info_vector, 
     &                               local_element_no)
c
      if( info_type .gt. 0 .and. info_type .le. 4 ) then
         value = info_vector(info_type)
      else
         write(out,9000) 5
         call die_gracefully
      end if
c
      return
c
 9000 format(/," SYSTEM Error: material_model_sizes.", /,
     &         "               error type: ", i5, /,
     &         "               Job aborted" )
      end


c     ****************************************************************
c     *                                                              *
c     *             subroutine rstgp1_update_strains                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/28/2015 rhd              *
c     *                                                              *
c     *      computed total strains at n+1 by including the          *
c     *      increment over the current step: n+1 = n + deps         *
c     *      for geometric nonlinear theory, these are the unrotated *
c     *      strains. for linear theory, just the usual small-       *
c     *      strains. strains are in vector (not tensor) form at     *
c     *      at this point. compiler should inline this routine      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rstgp1_update_strains( span, mxvl, nstr, deps,
     &                                  strain_np1 )
      implicit integer (a-z)
c
c                      parameter declarations
c
      double precision ::
     & deps(mxvl,*), strain_np1(mxvl,*)
c
       do i = 1, span
         strain_np1(i,1) = strain_np1(i,1) + deps(i,1)
         strain_np1(i,2) = strain_np1(i,2) + deps(i,2)
         strain_np1(i,3) = strain_np1(i,3) + deps(i,3)
         strain_np1(i,4) = strain_np1(i,4) + deps(i,4)
         strain_np1(i,5) = strain_np1(i,5) + deps(i,5)
         strain_np1(i,6) = strain_np1(i,6) + deps(i,6)
      end do
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                  subroutine gauss_pt_coords                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 09/28/2015 rhd             *
c     *                                                              *
c     *     compute (x,y,z) coordinates for all elements in block    *
c     *     at this integration point                                *
c     *                                                              *
c     ****************************************************************
c
c
c      subroutine gauss_pt_coords(
c     &       gpn, etype, span, int_order,
c     &       nnodel, gp_coords, node_coords, iout  )
c      implicit integer (a-z)
c      include 'param_def'
c
c                      parameter declarations
c
c      double precision ::
c     &  gp_coords(mxvl,3), node_coords(mxvl,*)
c
c                     locally defined arrays-variables
c
c      double precision ::
c     &  sf(mxndel), xi, eta, zeta, weight, zero
c      logical :: local_debug
c      data zero, local_debug / 0.0d00, .false. /
c
c      if( local_debug ) write(iout,*) '... in gauss_pt_coords'
c
c                     get the parametric coordinates for
c                     this integration point. then get the nodal
c                     shape functions evaluated at the point.
c
c      call getgpts( etype, int_order, gpn, xi, eta, zeta, weight )
c      call shapef( etype, xi, eta, zeta, sf(1) )
c
c      do i = 1, span
c         gp_coords(i,1) = zero
c         gp_coords(i,2) = zero
c         gp_coords(i,3) = zero
c      end do
c
c                     interpolate (x,y,z) global coords at this gpn for
c                     each element in block. rows of node_coords have
c                     coords for element nodes -- all x-coord, then
c                     all y-coord, then all z-coord
c
c      ky = nnodel
c      kz = ky + nnodel
c      do enode = 1, nnodel
c        do i = 1, span
c          gp_coords(i,1) = gp_coords(i,1)  +
c     &                      sf(enode) * node_coords(i,enode)
c          gp_coords(i,2) = gp_coords(i,2)  +
c     &                      sf(enode) *  node_coords(i,ky+enode)
c          gp_coords(i,3) = gp_coords(i,3)  +
c     &                      sf(enode) *  node_coords(i,kz+enode)
c        end do
c      end do
c
c      if( .not. local_debug ) return
c         write(iout,*) '>> in  gauss_pt_coords'
c         write(iout,*) 'xi, eta, zeta:'
c         write(iout,9000) xi, eta, zeta
c         write(iout,*) 'coords for element 1 of blk'
c         write(iout,9000) node_coords(1,1:3*nnodel)
c      return
c
c 9000 format(1x,f15.6 )
c      end
c     ****************************************************************
c     *                                                              *
c     *             subroutine rstgp1_make_symmetric_store           *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 10/27/20015 rhd            *
c     *                                                              *
c     *                   support for umat processing                *
c     *                                                              *
c     ****************************************************************
c
      subroutine rstgp1_make_symmetric_store( matrix, symm_vector )
      implicit none
      
      double precision :: matrix(6,6), symm_vector(21)
c
      double precision :: tp(6,6), symm_version(6,6), half
      integer :: i, j, k, map(6)
      data half / 0.5d00 /
      data map / 1,2,3,4,6,5 /
c
c         1. compute transpose of 6 x 6 matrix
c         2. compute symmetrized version
c         3. swap rows, cols 5 & 6 to make shear ordering
c            compatible with WARP3D
c         4. store 21 terms in lower triangle by row
c
      tp = transpose( matrix )
c
      do j = 1, 6
        do i = 1, 6
          symm_version(map(i),map(j)) = half * ( matrix(i,j) +
     &                                  tp(i,j) )
        end do
      end do
c
c      k = 1
c      do i = 1, 6
c        do j = 1, i
c          symm_vector(k) = symm_version(i,j)
c          k = k + 1
c        end do
c      end do

      symm_vector(1)  = symm_version(1,1)
      symm_vector(2)  = symm_version(2,1)
      symm_vector(3)  = symm_version(2,2)
      symm_vector(4)  = symm_version(3,1)
      symm_vector(5)  = symm_version(3,2)
      symm_vector(6)  = symm_version(3,3)
      symm_vector(7)  = symm_version(4,1)
      symm_vector(8)  = symm_version(4,2)
      symm_vector(9)  = symm_version(4,3)
      symm_vector(10) = symm_version(4,4)
      symm_vector(11) = symm_version(5,1)
      symm_vector(12) = symm_version(5,2)
      symm_vector(13) = symm_version(5,3)
      symm_vector(14) = symm_version(5,4)
      symm_vector(15) = symm_version(5,5)
      symm_vector(16) = symm_version(6,1)
      symm_vector(17) = symm_version(6,2)
      symm_vector(18) = symm_version(6,3)
      symm_vector(19) = symm_version(6,4)
      symm_vector(20) = symm_version(6,5)
      symm_vector(21) = symm_version(6,6)
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *    supporting routines for rstgp1 (to be inlined)            *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 10/18/15                       *
c     *                                                              *
c     *  support routines for rstgp1. include here so they can be    *
c     *  inlined.                                                    *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rstgp1_a( ndof, nnode, span, ue, due, uenh, uen1,
     &                     mxvl )
      integer :: span
      double precision ::
     &  ue(mxvl,*), due(mxvl,*), uenh(mxvl,*), uen1(mxvl,*),
     &  half
      data half / 0.5d00 /
c
      do j = 1, ndof*nnode
        do i = 1, span
           uenh(i,j) = ue(i,j) + half*due(i,j)
           uen1(i,j) = ue(i,j) + due(i,j)
        end do
      end do
c
      return
      end

      subroutine rstgp1_b( span, internal_energy, plastic_work,
     &                     gp_energies, gp_plast_work, det_j,
     &                     dfn1, itype )
      integer :: span
      double precision ::
     &  internal_energy, plastic_work, gp_energies(*),
     &  det_j(*), dfn1(*), gp_plast_work(*)
c
      if( itype .ne. 1 ) go to 100
      do i = 1, span
        internal_energy = internal_energy + gp_energies(i) *
     &                     dfn1(i) * det_j(i)
         plastic_work    = plastic_work +
     &                     gp_plast_work(i) * dfn1(i) * det_j(i)
      end do
      return
c
 100  continue
      do i = 1, span
         internal_energy = internal_energy + gp_energies(i) *
     &                       det_j(i)
         plastic_work    = plastic_work + gp_plast_work(i) * det_j(i)
      end do
      return
c
      end

c     ****************************************************************
c     *                                                              *
c     *                 subroutine recstr_cep_uddt_for_block         *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 06/24/12                   *
c     *                                                              *
c     *     compute updated stresses for elements in this block at   *
c     *     this integration point. Uses the stored [Dt] * deps      *
c     *     already adjusted for incremental temperature             *
c     *                                                              *
c     ****************************************************************

      subroutine recstr_cep_uddt_for_block( mxvl, span, ceps_blk,
     &     deps_blk, stress_n, stress_np1, type, nrow_ceps_blk,
     &     gpn )
      implicit none
c
c                      parameter declarations
c
      integer :: mxvl, span, type, nrow_ceps_blk, gpn
      double precision ::
     &  ceps_blk(nrow_ceps_blk,span,*), deps_blk(mxvl,*), 
     &  stress_n(mxvl,6), stress_np1(mxvl,6)
c
c                      locals

      integer :: ielem, i, j, k
      double precision ::
     & full_cep(mxvl,6,6), zero
      data zero  / 0.0d00 /
c
c              handle solid elements (type = 1) and cohesive elements
c              (type = 2 ) to let compiler optimize loops.
c
      if( type .eq. 2 ) go to 1000
c
c              expand compressed (symmetric) [Dts] to full 6x6
c              for simplicity in coding next loop.
c
      k = 1
      do i = 1, 6
       do j = 1, i
         do ielem = 1, span
          full_cep(ielem,i,j) = ceps_blk(k,ielem,gpn)
          full_cep(ielem,j,i) = full_cep(ielem,i,j)
         end do
        k = k + 1
       end do
      end do
c
c              compute stress @ n+1 = stress @ n + [Dt]* deps for
c              each element in block at this integration point
c
      stress_np1(1:mxvl,1:6) = stress_n(1:mxvl,1:6)
c
      do i = 1, 6
       do k = 1, 6
         do ielem = 1, span
           stress_np1(ielem,i) = stress_np1(ielem,i) +
     &         full_cep(ielem,i,k) * deps_blk(ielem,k)
         end do
       end do
      end do
c
      return
c
c              cohesive elements have 3x3 [Dt]
c
 1000 continue
      k = 1
      do i = 1, 3
       do j = 1, i
         do ielem = 1, span
          full_cep(ielem,i,j) = ceps_blk(k,ielem,gpn)
          full_cep(ielem,j,i) = full_cep(ielem,i,j)
         end do
        k = k + 1
       end do
      end do
c
c              compute stress @ n+1 = stress @ n + [Dt]* deps for each
c              element in block at this integration point
c
      stress_np1(1:mxvl,1:3) = stress_n(1:mxvl,1:3)
c
      do i = 1, 3
       do k = 1, 3
         do ielem = 1, span
           stress_np1(ielem,i) =  stress_np1(ielem,i) +
     &         full_cep(ielem,i,k) * deps_blk(ielem,k)
         end do
       end do
      end do
c      
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                 subroutine rstgp1_store_cep                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 10/25/2015 rhd             *
c     *                                                              *
c     *  store full cep from mm.. routine into symmetric global data *
c     *  structure for elements in block for gauss point gpn         *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rstgp1_store_cep( span, mxvl, gpn, gbl_ceps_blk, 
     &                             nrow_ceps_blk, local_cep )                
      implicit none
c
      integer :: span, mxvl, gpn, nrow_ceps_blk 
      double precision ::
     &  gbl_ceps_blk(nrow_ceps_blk,span,*), local_cep(mxvl,6,6)
c
      integer i, k, ii, jj
c
      if( nrow_ceps_blk .eq. 21 ) then ! symmetric [D] 6x6
        do i = 1, span
          gbl_ceps_blk(1,i,gpn)  = local_cep(i,1,1)
          gbl_ceps_blk(2,i,gpn)  = local_cep(i,2,1)
          gbl_ceps_blk(3,i,gpn)  = local_cep(i,2,2)
          gbl_ceps_blk(4,i,gpn)  = local_cep(i,3,1)
          gbl_ceps_blk(5,i,gpn)  = local_cep(i,3,2)
          gbl_ceps_blk(6,i,gpn)  = local_cep(i,3,3)
          gbl_ceps_blk(7,i,gpn)  = local_cep(i,4,1)
          gbl_ceps_blk(8,i,gpn)  = local_cep(i,4,2)
          gbl_ceps_blk(9,i,gpn)  = local_cep(i,4,3)
          gbl_ceps_blk(10,i,gpn) = local_cep(i,4,4)
       end do
       do i = 1, span   
          gbl_ceps_blk(11,i,gpn) = local_cep(i,5,1)
          gbl_ceps_blk(12,i,gpn) = local_cep(i,5,2)
          gbl_ceps_blk(13,i,gpn) = local_cep(i,5,3)
          gbl_ceps_blk(14,i,gpn) = local_cep(i,5,4)
          gbl_ceps_blk(15,i,gpn) = local_cep(i,5,5)
          gbl_ceps_blk(16,i,gpn) = local_cep(i,6,1)
          gbl_ceps_blk(17,i,gpn) = local_cep(i,6,2)
          gbl_ceps_blk(18,i,gpn) = local_cep(i,6,3)
          gbl_ceps_blk(19,i,gpn) = local_cep(i,6,4)
          gbl_ceps_blk(20,i,gpn) = local_cep(i,6,5)
          gbl_ceps_blk(21,i,gpn) = local_cep(i,6,6)
        end do
      elseif( nrow_ceps_blk .eq. 6 ) then ! symmetric [D] 3x3
        do i = 1, span
          gbl_ceps_blk(1,i,gpn)  = local_cep(i,1,1)
          gbl_ceps_blk(2,i,gpn)  = local_cep(i,2,1)
          gbl_ceps_blk(3,i,gpn)  = local_cep(i,2,2)
          gbl_ceps_blk(4,i,gpn)  = local_cep(i,3,1)
          gbl_ceps_blk(5,i,gpn)  = local_cep(i,3,2)
          gbl_ceps_blk(6,i,gpn)  = local_cep(i,3,3)
        end do
      elseif( nrow_ceps_blk .eq. 36 ) then ! non-symmetric [D] 6x6
          k = 1 
          do i = 1, span
            do ii = 1, 6
              do jj = 1, 6
                gbl_ceps_blk(k,i,gpn)  = local_cep(i,ii,jj)
                k = k +1
              end do
            end do   
          end do
      else
       write(*,9000)
       call die_abort
      end if
c
      return
 9000 format(/,3x,">>> FATAL ERROR: wrong. size. rstgp1_store_cep",
     &       /,3x,"                 job aborted" )       
      end
