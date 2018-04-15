c     ****************************************************************
c     *                                                              *
c     *                   subroutine recstr_allocate                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified :  3/15/2017 rhd             *
c     *                                                              *
c     *     allocate data structure in local_work for updating       *
c     *     strains-stresses-internal forces.                        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine recstr_allocate( local_work )
c      use segmental_curves, only : max_seg_points, max_seg_curves
      use elem_block_data, only: history_blk_list
      implicit none

      include 'common.main'
      include 'include_sig_up'

c     use segmental_curves
      integer, parameter :: max_seg_points=20, max_seg_curves=20

c     original declaration
      integer :: local_mt, error, span, blk, ngp, hist_size, nlsize
      double precision :: zero
      data zero / 0.0d00 /
c
      local_mt = local_work%mat_type
c
      allocate(
     &   local_work%ce_0(mxvl,mxecor),
     &   local_work%ce_n(mxvl,mxecor),
     &   local_work%ce_mid(mxvl,mxecor),
     &   local_work%ce_n1(mxvl,mxecor), stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 1
         call die_abort
      end if
c
      allocate( local_work%trnmte(mxvl,mxedof,mxndof) )
c        
      allocate(
     1 local_work%det_j(mxvl,mxgp),
     2 local_work%det_j_mid(mxvl,mxgp),
     3 local_work%nxi(mxndel,mxgp),
     4 local_work%neta(mxndel,mxgp),
     5 local_work%nzeta(mxndel,mxgp),
     6 local_work%gama(mxvl,3,3,mxgp),
     7 local_work%gama_mid(mxvl,3,3,mxgp), stat=error  )
      if( error .ne. 0 ) then
         write(out,9000) 2
         call die_abort
      end if
c
      if( local_work%geo_non_flg ) then
        allocate( local_work%fn(mxvl,3,3),
     &           local_work%fn1(mxvl,3,3),
     &           local_work%dfn1(mxvl) )
        local_work%fn1 = zero
        local_work%dfn1 = zero
      end if
      
      
      allocate(
     &  local_work%vol_block(mxvl,8,3),
     &  local_work%volume_block(mxvl),
     &  local_work%volume_block_0(mxvl),
     &  local_work%volume_block_n(mxvl),
     &  local_work%volume_block_n1(mxvl),
     &  local_work%jac(mxvl,3,3),
     &  local_work%b(mxvl,mxedof,nstr),
     &  local_work%ue(mxvl,mxedof),
     &  local_work%due(mxvl,mxedof),
     &  local_work%uenh(mxvl,mxedof), stat=error  )
      if( error .ne. 0 ) then
         write(out,9000) 3
         call die_abort
      end if
      
      local_work%b = zero      
c
      allocate( local_work%uen1(mxvl,mxedof),
     1  local_work%urcs_blk_n(mxvl,nstrs,mxgp),
     2  local_work%urcs_blk_n1(mxvl,nstrs,mxgp),
     3  local_work%rot_blk_n1(mxvl,9,mxgp),
     4  local_work%rtse(mxvl,nstr,mxgp), stat=error  )
      if( error .ne. 0 ) then
         write(out,9000) 4
         call die_abort
      end if
      local_work%det_j       = zero
      local_work%rot_blk_n1  = zero
      local_work%urcs_blk_n1 = zero
      local_work%urcs_blk_n  = zero
c
      allocate( local_work%ddtse(mxvl,nstr,mxgp),
     1   local_work%strain_n(mxvl,nstr,mxgp),
     2   local_work%dtemps_node_blk(mxvl,mxndel),
     3   local_work%temps_ref_node_blk(mxvl,mxndel),
     4   local_work%temps_node_blk(mxvl,mxndel),
     5   local_work%temps_node_ref_blk(mxvl,mxndel),
     6   local_work%nu_vec(mxvl),
     7   local_work%beta_vec(mxvl),
     8   local_work%h_vec(mxvl),
     9   local_work%tan_e_vec(mxvl),
     a   local_work%e_vec(mxvl), stat=error  )
      if( error .ne. 0 ) then
         write(out,9000) 5
         call die_abort
      end if
c
      allocate( local_work%sigyld_vec(mxvl),
     1   local_work%alpha_vec(mxvl,6),
     2   local_work%e_vec_n(mxvl),
     3   local_work%nu_vec_n(mxvl),
     4   local_work%gp_sig_0_vec(mxvl),
     5   local_work%gp_sig_0_vec_n(mxvl),
     6   local_work%gp_h_u_vec(mxvl),
     7   local_work%gp_h_u_vec_n(mxvl),
     8   local_work%gp_beta_u_vec(mxvl),
     9   local_work%gp_beta_u_vec_n(mxvl), stat=error  )
      if( error .ne. 0 ) then
         write(out,9000) 6
         call die_abort
      end if
c
      allocate( local_work%gp_delta_u_vec(mxvl),
     1   local_work%gp_delta_u_vec_n(mxvl),
     2   local_work%alpha_vec_n(mxvl,6),
     3   local_work%h_vec_n(mxvl),
     4   local_work%n_power_vec(mxvl) )
c     
      local_mt = local_work%mat_type
c      
      if( local_mt .eq. 3 ) then
        allocate( 
     1   local_work%f0_vec(mxvl),
     2   local_work%eps_ref_vec(mxvl),
     3   local_work%m_power_vec(mxvl),
     4   local_work%q1_vec(mxvl),
     5   local_work%q2_vec(mxvl),
     6   local_work%q3_vec(mxvl),
     7   local_work%nuc_s_n_vec(mxvl),
     8   local_work%nuc_e_n_vec(mxvl),
     9   local_work%nuc_f_n_vec(mxvl), stat=error  )
        if( error .ne. 0 ) then
         write(out,9000) 7
         call die_abort
        end if
      end if
c
      allocate( local_work%eps_curve(max_seg_points),
     1    local_work%shape(mxndel,mxgp),
     5    local_work%enode_mat_props(mxndel,mxvl,mxndpr),
     6    local_work%fgm_flags(mxvl,mxndpr), stat=error  )
      if( error .ne. 0 ) then
         write(out,9000) 9
         call die_abort
      end if
c
      if( local_mt .eq. 5 ) allocate( local_work%mm05_props(mxvl,10) )
      if( local_mt .eq. 6 ) allocate( local_work%mm06_props(mxvl,10) )
      if( local_mt .eq. 7 ) allocate( local_work%mm07_props(mxvl,10) )
      if( local_work%is_umat ) then
        allocate(local_work%umat_props(mxvl,50),
     1           local_work%characteristic_length(mxvl) )
      end if

c
      allocate( local_work%trne(mxvl,mxndel), stat=error  )
      if( error .ne. 0 ) then
         write(out,9000) 11
         call die_abort
      end if
      local_work%trne = .false.
c
      if( local_mt .eq. 10 .or. local_mt .eq. 11 ) then
        allocate( local_work%debug_flag(mxvl),
     1    local_work%local_tol(mxvl),
     2    local_work%ncrystals(mxvl),
     3    local_work%angle_type(mxvl),
     4    local_work%angle_convention(mxvl),
     5    local_work%c_props(mxvl,max_crystals),
     6    local_work%nstacks(mxvl),
     7    local_work%nper(mxvl))
      end if
c
      span                         = local_work%span
      blk                          = local_work%blk
      ngp                          = local_work%num_int_points
      hist_size                    = history_blk_list(blk)
      local_work%hist_size_for_blk = hist_size
c
      allocate( local_work%elem_hist1(span,hist_size,ngp),
     &          local_work%elem_hist(span,hist_size,ngp), stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 12
         call die_abort
      end if
c
      if( local_work%is_cohes_nonlocal ) then
         nlsize = nonlocal_shared_state_size
         allocate( local_work%top_surf_solid_stresses_n(mxvl,nstrs),
     &      local_work%bott_surf_solid_stresses_n(mxvl,nstrs),
     &      local_work%top_surf_solid_eps_n(mxvl,nstr),
     &      local_work%bott_surf_solid_eps_n(mxvl,nstr),
     &      local_work%top_surf_solid_elements(mxvl),
     &      local_work%bott_surf_solid_elements(mxvl),
     &      local_work%nonlocal_stvals_bott_n(mxvl,nlsize),
     &      local_work%nonlocal_stvals_top_n(mxvl,nlsize),
     &      stat=error )
         if( error .ne. 0 ) then
           write(out,9000) 13
           call die_abort
         end if
         allocate( local_work%top_solid_matl(mxvl),
     &      local_work%bott_solid_matl(mxvl), stat=error )
         if( error .ne. 0 ) then
           write(out,9000) 14
           call die_abort
         end if
      end if
c
c             always allocate cohes_rot_block even if not
c             a block of cohesive elements. the
c             array is passed in lots of places
c             that process both solid and interface elements.
c             So it needs to exist !

      allocate( local_work%cohes_rot_block(mxvl,3,3), stat=error ) 
      if( error .ne. 0 ) then
         write(out,9000) 15
         call die_abort
      end if
c         
      if( local_work%is_cohes_elem ) then
         allocate( local_work%cohes_temp_ref(mxvl),
     1      local_work%cohes_dtemp(mxvl),
     2      local_work%cohes_temp_n(mxvl),
     3      local_work%intf_prp_block(mxvl,max_interface_props),
     4      stat=error )
         if( error .ne. 0 ) then
           write(out,9000) 16
           call die_abort
         end if
      end if
c
c             when nonlocal analysis in effect, allocate space
c             for block nonlocal state values. we need to store span x 
c             num local values (fixed global in param_def). Space
c             will be passed to material model code to insert values
c             for elements in the block at the integration point
c             being processed.
c             allocate dummy 1x1 for std. local analysis to
c             simplify calls later.
c
      if( local_work%block_has_nonlocal_solids ) then
         allocate( local_work%nonlocal_state_blk(mxvl,
     &             nonlocal_shared_state_size), stat=error ) 
         if( error .ne. 0 ) then
           write(out,9000) 17
           call die_abort
         end if
      else 
         allocate( local_work%nonlocal_state_blk(1,1), stat=error ) 
         if( error .ne. 0 ) then
           write(out,9000) 18
           call die_abort
         end if
      end if
c
      allocate( local_work%weights(mxgp), stat=error ) 
      if( error .ne. 0 ) then
         write(out,9000) 19
         call die_abort
      end if      
c
      allocate( local_work%sv(3), local_work%lv(3), 
     &          local_work%tv(3), stat=error ) 
      if( error .ne. 0 ) then
         write(out,9000) 20
         call die_abort
      end if      
c
      return
c
 9000 format('>> FATAL ERROR: recstr_allocate'
     &  /,   '                failure status= ',i5,
     &  /,   '                job terminated' )
c
      end

c     ****************************************************************
c     *                                                              *
c     *                   subroutine recstr_deallocate               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified :  3/15/2017 rhd             *
c     *                                                              *
c     *     release data structure in local_work for updating        *
c     *     strains-stresses-internal forces.                        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine recstr_deallocate( local_work )
      implicit none 
      include 'common.main'
      include 'include_sig_up'
c
      integer :: local_mt, error
      logical :: local_debug
c
      local_debug = .false.
      if( local_debug ) write(out,*) "..recstr_dell @ 1"
      local_mt = local_work%mat_type
c
      
      deallocate(
     & local_work%ce_0,
     & local_work%ce_n,
     & local_work%ce_mid,
     & local_work%ce_n1, stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 1
         call die_abort
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 2"
      deallocate( local_work%trnmte, stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 100
         call die_abort
      end if
c      
      if( local_debug ) write(out,*) "..recstr_dell @ 3"
      deallocate(
     1 local_work%det_j,
     2 local_work%det_j_mid,
     3 local_work%nxi,
     4 local_work%neta,
     5 local_work%nzeta,
     6 local_work%gama,
     7 local_work%gama_mid, stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 2
         call die_abort
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 4"
      if( local_work%geo_non_flg ) then
        deallocate( local_work%fn,
     &              local_work%fn1,
     &              local_work%dfn1, stat=error )
        if( error .ne. 0 ) then
         write(out,9000) 102
         call die_abort
        end if
      end if

      if( local_debug ) write(out,*) "..recstr_dell @ 5"
      deallocate( 
     &  local_work%vol_block,
     &  local_work%volume_block, local_work%volume_block_0,
     &  local_work%volume_block_n, local_work%volume_block_n1,
     &  local_work%jac,
     &  local_work%b,
     &  local_work%ue,
     &  local_work%due,
     &  local_work%uenh,
     &  local_work%uen1,
     &  local_work%urcs_blk_n,
     &  local_work%urcs_blk_n1,
     &  local_work%rot_blk_n1,
     &  local_work%rtse, stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 3
         call die_abort
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 6"
      deallocate( local_work%ddtse,
     1   local_work%strain_n,
     2   local_work%dtemps_node_blk,
     3   local_work%temps_ref_node_blk,
     4   local_work%temps_node_blk,
     5   local_work%temps_node_ref_blk,
     6   local_work%nu_vec,
     7   local_work%beta_vec,
     8   local_work%h_vec,
     9   local_work%tan_e_vec,     
     a   local_work%e_vec, stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 4
         call die_abort
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 7"
      deallocate( local_work%sigyld_vec,
     1   local_work%alpha_vec,
     2   local_work%e_vec_n,
     3   local_work%nu_vec_n,
     4   local_work%gp_sig_0_vec,
     5   local_work%gp_sig_0_vec_n,
     6   local_work%gp_h_u_vec,
     7   local_work%gp_h_u_vec_n,
     8   local_work%gp_beta_u_vec,
     9   local_work%gp_beta_u_vec_n, stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 5
         call die_abort
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 8"
      deallocate( local_work%gp_delta_u_vec,
     1   local_work%gp_delta_u_vec_n,
     2   local_work%alpha_vec_n,
     3   local_work%h_vec_n,
     4   local_work%n_power_vec, stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 104
         call die_abort
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 9"
      if( local_mt .eq. 3 ) then
        deallocate(      
     1   local_work%f0_vec,
     2   local_work%eps_ref_vec,
     3   local_work%m_power_vec,
     4   local_work%q1_vec,
     5   local_work%q2_vec,
     6   local_work%q3_vec,
     7   local_work%nuc_s_n_vec,
     8   local_work%nuc_e_n_vec,
     9   local_work%nuc_f_n_vec, stat=error )
        if( error .ne. 0 ) then
         write(out,9000) 6
         call die_abort
        end if
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 10"
      deallocate( local_work%eps_curve,
     1    local_work%shape,
     5    local_work%enode_mat_props,
     6    local_work%fgm_flags, stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 7
         call die_abort
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 11"
      if( local_mt .eq. 5 ) then
        deallocate( local_work%mm05_props, stat=error )
        if( error .ne. 0 ) then
         write(out,9000) 106
         call die_abort
        end if
      end if  
c
      if( local_debug ) write(out,*) "..recstr_dell @ 12"
      if( local_mt .eq. 6 ) then
        deallocate( local_work%mm06_props, stat=error )
        if( error .ne. 0 ) then
         write(out,9000) 108
         call die_abort
        end if
      end if  
c
      if( local_debug ) write(out,*) "..recstr_dell @ 13"
      if( local_mt .eq. 7 ) then
        deallocate( local_work%mm07_props, stat=error )
        if( error .ne. 0 ) then
         write(out,9000) 110
         call die_abort
        end if
      end if  
c
      if( local_debug ) write(out,*) "..recstr_dell @ 14"
      if( local_work%is_umat ) then
        deallocate( local_work%umat_props,
     &              local_work%characteristic_length, stat=error )
        if( error .ne. 0 ) then
         write(out,9000) 112
         call die_abort
        end if
      end if  
c
      if( local_debug ) write(out,*) "..recstr_dell @ 15"
      deallocate( local_work%trne, stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 9
         call die_abort
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 16"
      deallocate( local_work%elem_hist, stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 10
         call die_abort
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 16.5"
      deallocate( local_work%elem_hist1, stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 10
         call die_abort
      end if

c
      if( local_debug ) write(out,*) "..recstr_dell @ 17"
      if( local_mt .eq. 10 .or. local_mt .eq. 11 ) then
        deallocate( local_work%debug_flag,
     1    local_work%local_tol,
     2    local_work%ncrystals,
     3    local_work%angle_type,
     4    local_work%angle_convention,
     5    local_work%c_props,
     6    local_work%nstacks,
     7    local_work%nper, stat=error)
        if( error .ne. 0 ) then
         write(out,9000) 114
         call die_abort
        end if
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 18"
      if( local_work%is_cohes_nonlocal ) then
         deallocate( local_work%top_surf_solid_stresses_n,
     1      local_work%bott_surf_solid_stresses_n,
     2      local_work%top_surf_solid_eps_n,
     3      local_work%bott_surf_solid_eps_n,
     4      local_work%top_surf_solid_elements,
     5      local_work%bott_surf_solid_elements,
     6      local_work%nonlocal_stvals_bott_n,
     7      local_work%nonlocal_stvals_top_n,
     8      stat=error )
         if( error .ne. 0 ) then
           write(out,9000) 12
           call die_abort
         end if
         if( local_debug ) write(out,*) "..recstr_dell @ 19"
         deallocate( local_work%top_solid_matl,
     &      local_work%bott_solid_matl, stat=error )
         if( error .ne. 0 ) then
           write(out,9000) 13
           call die_abort
         end if
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 20"
      deallocate( local_work%cohes_rot_block, stat=error )
      if( error .ne. 0 ) then
           write(out,9000) 14
           call die_abort
      end if
      if( local_work%is_cohes_elem ) then
         deallocate( local_work%cohes_temp_ref,
     1      local_work%cohes_dtemp,
     2      local_work%cohes_temp_n, 
     3      local_work%intf_prp_block, stat=error )
         if( error .ne. 0 ) then
           write(out,9000) 14
           call die_abort
         end if
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 21"
      deallocate( local_work%nonlocal_state_blk, stat=error ) 
         if( error .ne. 0 ) then
           write(out,9000) 15
           call die_abort
         end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 22"
      deallocate( local_work%weights, stat=error ) 
         if( error .ne. 0 ) then
           write(out,9000) 16
           call die_abort
         end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 23"
      deallocate( local_work%sv, local_work%lv, local_work%tv,
     &            stat=error ) 
         if( error .ne. 0 ) then
           write(out,9000) 17
           call die_abort
         end if
      return
c
 9000 format('>> FATAL ERROR: recstr_deallocate'
     &  /,   '                failure status= ',i5,
     &  /,   '                job terminated' )
c
      end

