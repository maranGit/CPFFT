c     ****************************************************************
c     *                                                              *
c     *                   subroutine dupstr_blocked                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/9/2016 rhd               *
c     *                                                              *
c     *     creates a separate copy of element                       *
c     *     blocked data necessary for global stress vector recovery *
c     *     each element in  block                                   *
c     *                                                              *
c     ****************************************************************
c
c
c     subroutine dupstr_blocked(
c    & blk, span, felem, ngp, nnode, ndof, totdof, mat_type, 
c    & geonl, step, iter, belinc, bcdst, ce_0, ce_n, ce_mid, ce_n1,
c    & ue, due, local_work )
      subroutine dupstr_blocked( blk, span, felem, mat_type, 
     & geonl, step, iter, local_work )
c
      use elem_block_data, only:  history_blocks,
     &                            eps_n_blocks, urcs_n_blocks,
     &                            history_blk_list
c     use main_data,       only:  dtemp_nodes, dtemp_elems,
c    &                            temper_nodes, temper_elems,
c    &                            temper_nodes_ref, temperatures_ref
c    &                            fgm_node_values,
c    &                            fgm_node_values_defined,
c    &                            fgm_node_values_cols, matprp,
c    &                            lmtprp
c
c     use segmental_curves, only : max_seg_points, max_seg_curves
c
      implicit none
      include 'common.main'
      include 'include_sig_up'
c
c           parameter declarations
c
      integer :: blk, span, felem, ngp, nnode, ndof, totdof, mat_type,
     &           step, iter, hist_size
      logical :: geonl
c     integer :: belinc(nnode,*), bcdst(totdof,*)
c     double precision ::
c    & ce_0(mxvl,*), ce_mid(mxvl,*), ue(mxvl,*), due(mxvl,*),
c    & ce_n(mxvl,*), ce_n1(mxvl,*)
c
c           local declarations
c
      integer :: elem_type, surf, k, j, i, matl_no
      logical ::local_debug, update, update_coords, middle_surface
      double precision ::
     &   half, zero, one, mag, mags(3), djcoh(mxvl)
      data local_debug, half, zero, one
     &  / .false., 0.5d00, 0.0d00, 1.0d00 /
c!DIR$ ASSUME_ALIGNED ce_0:64, ce_N:64, ce_mid:64, ce_n1:64  

c
      if( local_debug ) write(out,9100)
c
      elem_type      = local_work%elem_type
      surf           = local_work%surface
      middle_surface = surf .eq. 2
      djcoh(1:span)  = zero
c

c           pull coordinates at t=0 from global input vector.
c
c     k = 1
c     do j = 1, nnode
c !DIR$ LOOP COUNT MAX=128  
c !DIR$ IVDEP
c        do i = 1, span
c           ce_0(i,k)   = c(bcdst(k,i))
c           ce_0(i,k+1) = c(bcdst(k+1,i))
c           ce_0(i,k+2) = c(bcdst(k+2,i))
c        end do
c        k = k + 3
c     end do
c
c           for geonl, create a set of nodal coordinates
c           at mid-step and at the end of step. for step 1,
c           iter=0 we don't update for imposed displacements
c           since we are using linear stiffness.
c           if any elements are killed, their displacements
c           were set to zero above so we are just
c           setting the initial node coords. for iter=0 in other
c           steps, just make the mid and n1 coordinates the
c           coordinates at start of step. iter=0 means we are
c           computing equiv nodal forces for imposed/extrapolated
c           displacements to start a step.
c
      update = .true.
      if( step .eq. 1 .and. iter .eq. 0 ) update = .false.
      update_coords = geonl .and. update
c
c     if( update_coords ) then
c      do  j = 1, totdof
c!DIR$ LOOP COUNT MAX=128  
c!DIR$ IVDEP
c         do i = 1, span
c           ce_n(i,j)   = ce_0(i,j) + ue(i,j)
c           ce_mid(i,j) = ce_0(i,j) + ue(i,j) + half*due(i,j)
c           ce_n1(i,j)  = ce_0(i,j) + ue(i,j) + due(i,j)
c         end do
c      end do
c
c        if( local_work%is_cohes_elem ) then ! geonl interface elem
c         if( middle_surface ) call chk_cohes_penetrate( span, mxvl,
c    &               felem, mxndel, nnode, elem_type, ce_n1, djcoh )
c         call cohes_ref_surface( span, mxvl, mxecor,
c    &                            local_work%surface, nnode,
c    &                            totdof, ce_0, ce_n1, djcoh )
c         end if
c     end if   !  update_coords
c
c     if( .not. update_coords ) then
c         do  j = 1, totdof
c!DIR$ LOOP COUNT MAX=128  
c!DIR$ IVDEP
c           do i = 1, span
c             ce_n(i,j)   = ce_0(i,j)
c             ce_mid(i,j) = ce_0(i,j)
c             ce_n1(i,j)  = ce_0(i,j)
c           end do
c         end do
c     end if
c
c           gather nodal and element temperature change over load step
c           (if they are defined). construct a set of incremental nodal
c           temperatures for each element in block.
c
c     if( temperatures ) then
c       if( local_debug )  write(out,9610)
c       call gadtemps( dtemp_nodes, dtemp_elems(felem), belinc,
c    &                 nnode, span, felem, local_work%dtemps_node_blk,
c    &                 mxvl )
c     else
c       local_work%dtemps_node_blk(1:span,1:nnode) = zero
c     end if
c
c           gather reference temperatures for element nodes from the
c           global vector of reference values at(if they are defined).
c           construct a set of reference nodal temperatures for each
c           element in block.
c
c     if( temperatures_ref ) then
c       if( local_debug )  write(out,9610)
c       call gartemps( temper_nodes_ref, belinc, nnode, span,
c    &                 felem, local_work%temps_ref_node_blk, mxvl )
c     else
c       local_work%temps_ref_node_blk(1:span,1:nnode) = zero
c     end if
c
c           build nodal temperatures for elements in the block
c           at end of step (includes both imposed nodal and element
c           temperatures)
c
      if( local_debug )  write(out,9620)
c     call gatemps( temper_nodes, temper_elems(felem), belinc,
c    &              nnode, span, felem, local_work%temps_node_blk,
c    &              mxvl, local_work%dtemps_node_blk,
c    &              local_work%temps_node_to_process )

c
c           if the model has fgm properties at the model nodes, build a
c           table of values for nodes of elements in the block
c
c     if( fgm_node_values_defined ) then
c       do j = 1,  fgm_node_values_cols
c         do i = 1, nnode
c!DIR$ LOOP COUNT MAX=128 
c!DIR$ IVDEP 
c           do k = 1, span
c             local_work%enode_mat_props(i,k,j) =
c    &                     fgm_node_values(belinc(i,k),j)
c           end do
c         end do
c       end do
c     end if
c
c
c           gather material specific data for elements
c           in the block. we split operations based on
c           the material model associated with block
c
      if( local_debug )  write(out,9600)
c
c           gather element data at n from global blocks:
c            a) stresses -  unrotated cauchy stresses for geonl
c            b) strain_n -  strains at n
c            b) ddtse -     strains at start of step n, subsequently
c                           updated to strains at n+1 during strain-
c                           stress updating
c            c) elem_hist - integration point history data for material
c                           models
c           History data:
c            o The global blocks are sized(hist_size,ngp,span)
c            o The local block is sized (span,hist_size,ngp).
c              This makes it possible to pass a 2-D array slice for
c              all elements of the block for a single integration pt.
c
      hist_size = local_work%hist_size_for_blk
c
      call dptstf_copy_history(
     &  local_work%elem_hist(1,1,1), history_blocks(blk)%ptr(1),
     &            ngp, hist_size, span )
c
      call recstr_gastr( local_work%ddtse, eps_n_blocks(blk)%ptr(1),
     &                   ngp, nstr, span )
      call recstr_gastr( local_work%strain_n, eps_n_blocks(blk)%ptr(1),
     &                   ngp, nstr, span )
c
      call recstr_gastr( local_work%urcs_blk_n,
     &                   urcs_n_blocks(blk)%ptr(1),
     &                   ngp, nstrs, span )
c
c
      select case ( mat_type )
c
      case ( 1 )
c
c           vectorized mises plasticty model.
c
c!DIR$ LOOP COUNT MAX=128 
c!DIR$ IVDEP 
c       do i = 1, span
c          matl_no = iprops(38,felem+i-1)
c          local_work%tan_e_vec(i) = matprp(4,matl_no)
c       end do
c
      case ( 2 )
c
c           nonlinear elastic material model (deformation plasticity).
c
         if ( local_debug ) write(out,9950)
c
c
      case ( 3 )
c
c           general mises/gurson model.
c
        if ( local_debug ) write(out,9950)
c!DIR$ LOOP COUNT MAX=128  
c!DIR$ IVDEP
c       do i = 1, span
c          matl_no = iprops(38,felem+i-1)
c          local_work%tan_e_vec(i) = matprp(4,matl_no)
c       end do
c
      case ( 4 )
c
c           cohesive zone models
c
c        if ( local_debug ) write(out,9950)
c        if( local_work%is_cohes_nonlocal ) then
c          call recstr_build_cohes_nonlocal( local_work, iprops )
c        end if
c
      case ( 5 )
c
c           advanced cyclic plasticity model
c
        if ( local_debug ) write(out,9950)
c
      case ( 6 )
c
c           creep model
c
        if ( local_debug ) write(out,9960)

      case ( 7 )
c
c           mises model + hydrogen
c
        if ( local_debug ) write(out,9970)

      case ( 8 )
c
c           Abaqus compatible UMAT
c
        if ( local_debug ) write(out,9980)
c
      case ( 10 )
c
c           CP model
c
        if ( local_debug) write(out,9990)
c
c
      case default
          write(out,*) '>>> invalid material model number'
          write(out,*) '    in dupstr_blocked'
          call die_abort
          stop
c
      end select
c
      if ( local_debug ) write(out,9150)
c
      return
c
 9100 format(8x,'>> entered dupstr_blocked...' )
 9150 format(8x,'>> leaving dupstr_blocked...' )
 9600 format(12x,'>> gather element stresses/strains at step n...' )
 9610 format(12x,'>> gather element incremental temperatures...' )
 9620 format(12x,'>> gather element total temperatures...' )
 9800 format(12x,'>> gather material data for type 1...' )
 9900 format(15x,'>> gather plast. parms, back stress, state...' )
 9950 format(12x,'>> gather data for model type 5...' )
 9960 format(12x,'>> gather data for model type creep...' )
 9970 format(12x,'>> gather data for model type 7...' )
 9980 format(12x,'>> gather data for model type 8...' )
 9990 format(12x,'>> gather data for model type 10...' )
c
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine recstr_gastr                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 01/19/2013 rhd             *
c     *                                                              *
c     *     gathers element stresses from the global stress data     *
c     *     structure to local block array                           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine recstr_gastr( mlocal, mglobal, ngp, nprm, span )
      implicit none 
      include 'param_def'
c
c               parameter declarations
c
      integer :: ngp, nprm, span
      double precision ::
     & mlocal(mxvl,nprm,*), mglobal(nprm,ngp,*)
c
      integer :: k, j, i      
c!DIR$ ASSUME_ALIGNED mlocal:64, mglobal:64  
c
c           unroll inner loop for most common number of integration
c           points (ngp).
c
c
      if( ngp .ne. 8 ) then
        do k = 1, ngp
         do  j = 1, nprm
!DIR$ LOOP COUNT MAX=128  
!DIR$ VECTOR ALIGNED
            do  i = 1, span
               mlocal(i,j,k) = mglobal(j,k,i)
            end do
         end do
        end do
        return
      end if
c
c                number of integration points = 8, unroll.
c
      do  j = 1, nprm
!DIR$ LOOP COUNT MAX=128  
!DIR$ VECTOR ALIGNED
        do  i = 1, span
            mlocal(i,j,1) = mglobal(j,1,i)
            mlocal(i,j,2) = mglobal(j,2,i)
            mlocal(i,j,3) = mglobal(j,3,i)
            mlocal(i,j,4) = mglobal(j,4,i)
            mlocal(i,j,5) = mglobal(j,5,i)
            mlocal(i,j,6) = mglobal(j,6,i)
            mlocal(i,j,7) = mglobal(j,7,i)
            mlocal(i,j,8) = mglobal(j,8,i)
        end do
      end do
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *               subroutine dptstf_copy_history                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 09/27/2015 rhd             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine dptstf_copy_history( local_hist, global_hist,
     &                                ngp, hist_size, span )
      implicit none
c
c               parameter declarations
c
      integer :: ngp, hist_size, span
      double precision ::
     & local_hist(span,hist_size,ngp),
     & global_hist(hist_size,ngp,span)
c
      integer :: k, j, i     
c      if( ngp .ne. 8 ) then
        do k = 1, ngp
         do  j = 1, hist_size
!DIR$ LOOP COUNT MAX=128  
!DIR$ VECTOR ALIGNED
            do  i = 1, span
               local_hist(i,j,k) = global_hist(j,k,i)
            end do
         end do
        end do
        return
c      end if
c
c                number of gauss points = 8, unroll.
c
c     do  j = 1, hist_size
c@!DIR$ LOOP COUNT MAX=128  
c@!DIR$ IVDEP
c        do  i = 1, span
c            local_hist(i,j,1) = global_hist(j,1,i)
c            local_hist(i,j,2) = global_hist(j,2,i)
c            local_hist(i,j,3) = global_hist(j,3,i)
c            local_hist(i,j,4) = global_hist(j,4,i)
c            local_hist(i,j,5) = global_hist(j,5,i)
c            local_hist(i,j,6) = global_hist(j,6,i)
c            local_hist(i,j,7) = global_hist(j,7,i)
c            local_hist(i,j,8) = global_hist(j,8,i)
c        end do
c      end do
c
      return
      end
