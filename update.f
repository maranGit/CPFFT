c     ****************************************************************
c     *                                                              *
c     *                      subroutine update                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 2/17/2017 rhd              *
c     *                                                              *
c     *     various updates of vectors required after the iterative  *
c     *     solution procedure for a step has been completed.        *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine update
c
c      use main_data, only :  temper_nodes, temper_elems,
c     &                       dtemp_nodes, dtemp_elems, mdiag
      use elem_block_data, only : history_blocks, history1_blocks,
     &                            eps_n_blocks, eps_n1_blocks,
     &                            urcs_n_blocks, urcs_n1_blocks,
     &                            history_blk_list, eps_blk_list,
     &                            urcs_blk_list
c     use stiffness_data, only : total_lagrange_forces, 
c    &                           d_lagrange_forces
c
      implicit integer (a-z)
c      
      include 'common.main'

      double precision ::
     &    dlf, tlf
c     
      integer :: blk, felem, span, ngp, hblock_size, eblock_size, 
     &           ublock_size, k, n, i, dof
      logical :: chk, update_lag_forces
      logical, parameter :: local_debug = .false.
c
      if( local_debug ) write(out,9300)
c
c                       tell those slaves to get in here.
c
c     call wmpi_alert_slaves (17)
c
c                       update displacements
c
      u(1:nodof) = u(1:nodof) + du(1:nodof)
c
c                       update the velocities and accelerations
c                       to n+1 using Newmark beta method, velocities
c                       and accelerations at n and the converged 
c                       displacement increment (du) over the step.
c                       skip if we are on worker ranks
c
c     if( root_processor ) call newmrk( nodof, nbeta, dt, du, v, a )
c
c                       compute the equivalent loads for the con-
c                       strained dof.
c
c      dof = csthed
c      do while ( dof .ne. -1 )
c         load(dof) = mdiag(dof)*a(dof) + ifv(dof)
c         dof       = cstmap(dof)
c      end do
c
c                      update the structural vectors so that states
c                      at n and n+1 are identical
c                        1) element histories
c                        2) element strains
c                        3) element stresses
c
c$OMP PARALLEL DO PRIVATE( blk, felem, span, ngp, hblock_size,
c$OMP&                     eblock_size, ublock_size )
      do blk = 1, nelblk
         if( myid .ne. elblks(2,blk) ) cycle
         felem        = elblks(1,blk)
         span         = elblks(0,blk)
c        ngp          = iprops(6,felem)
         ngp          = fftngp
         hblock_size  = span * ngp * history_blk_list(blk)
         eblock_size  = span * ngp * nstr
         ublock_size  = span * ngp * nstrs
c
         if( hblock_size .gt. 0 )
     &      call update_copy( history_blocks(blk)%ptr(1),
     &                  history1_blocks(blk)%ptr(1), hblock_size )
         if( eps_blk_list(blk) .eq. 1 )
     &      call update_copy( eps_n_blocks(blk)%ptr(1),
     &                  eps_n1_blocks(blk)%ptr(1), eblock_size)
         if( urcs_blk_list(blk) .eq. 1 )
     &      call update_copy( urcs_n_blocks(blk)%ptr(1),
     &                  urcs_n1_blocks(blk)%ptr(1), ublock_size )
         if( local_debug ) then
           write(*,*) "span: ", span
           write(*,*) "ngp: ", ngp
           write(*,*) "nstrs: ", nstrs
           write(*,*) "span: ", span
           write(*,*) "span: ", span
           write(*,*) "Ran is here!!"
           write(*,*) "current block is: ", blk
           write(*,*) hblock_size, eps_blk_list(blk), urcs_blk_list(blk)
           write(*,*) urcs_n_blocks(blk)%ptr(1)
           write(*,*) urcs_n1_blocks(blk)%ptr(1)
         endif
      end do ! on blk
c$OMP END PARALLEL DO
c
c                      update nonlocal shared state values 
c                      if they exist
c
c      if( nonlocal_analysis ) then
c       n = nonlocal_shared_state_size
c       do i = 1, noelem
c         if( .not. nonlocal_flags(i) ) cycle
c         chk = allocated( nonlocal_data_n(i)%state_values ) .and.
c     &         allocated( nonlocal_data_n1(i)%state_values )
c        if( chk ) then
c!DIR$ VECTOR ALIGNED        
c           nonlocal_data_n(i)%state_values(1:n) =
c     &              nonlocal_data_n1(i)%state_values(1:n)
c        else
c           write(out,9100) i
c           call die_abort
c        end if
c       end do
c      end if
c
c                       update the nodal and element temperatures
c
c      if( temperatures ) then
c!DIR$ VECTOR ALIGNED      
c        temper_nodes(1:nonode) = temper_nodes(1:nonode) + 
c     &                           dtemp_nodes(1:nonode)
c!DIR$ VECTOR ALIGNED
c        temper_elems(1:noelem) = temper_elems(1:noelem) +
c     &                           dtemp_elems(1:noelem)
c      end if
c
c                       update contact geometry terms if needed
c
c     call updt_contact ! does not change contact forces
c
c                       update lagrange nodal forces for MPCs.
c                       not supported on MPI
c
c      update_lag_forces = allocated(total_lagrange_forces) .and.
c     &                    allocated(d_lagrange_forces) 
c      if( update_lag_forces ) then
c!DIR$ VECTOR ALIGNED      
c       total_lagrange_forces(1:nodof) = total_lagrange_forces(1:nodof)
c     &                                 + d_lagrange_forces(1:nodof)
c       if( local_debug ) then
c          write(out,9205); write(out,9210)
c          do i = 1, nodof
c            dlf = d_lagrange_forces(i)
c            tlf = total_lagrange_forces(i)
c            write(out,9220) i, tlf, tlf-dlf, dlf
c          end do
c       end if
c      end if  
c     
      if( local_debug ) write(out,9305)
c
      return
c           
 9200 format(10x,i5,3f15.5)     
 9190 format(15x,"pbar",10x,"bob_lag",10x,"dlag")
      return
c
 9000 format(3x,i5, 2e14.6)
 9100 format(">>>>> FATAL ERROR. update. nonlocal. elem: ",i8,
     &      /,"      Job terminated." )
 9205 format(5x,"... updating Lagrange forces to n+1 ...")
 9210 format(t5,"sdof", t16, "tlf (new)", t31, "tlf(old)", 
     &  t50, "dlf" )
 9220 format(2x,i6,3f16.6)
 9300 format(1x,"--- entering update ---" )
 9305 format(1x,"--- leaving update ---" )
     
c
      end
c     ****************************************************************
c     *                                                              *
c     *                    subroutine update_copy                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 04/20/2015 rhd             *
c     *                                                              *
c     *     copy vector a = b. using this routine hides the blocked  *
c     *     structures of input arrays from pointers. should get     *
c     *     inlined                                                  *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine update_copy( a, b, n )
      implicit none
      integer :: n
      double precision :: a(n), b(n)
c
!DIR$ VECTOR ALIGNED
      a = b
      return
      end
