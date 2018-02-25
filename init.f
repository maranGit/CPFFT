      subroutine history_cep_init()
c
      use elem_block_data, only : history_blocks, history1_blocks,
     &                            history_blk_list,
     &                            cep_blocks, cep_blk_list
c
      implicit integer (a-z)
      include 'common.main'
      double precision
     &    dummy(1)
      integer matl_info(10)
      logical myblk
c
c
c            proc_type :  = 1, create data arrays only and zero
c                         = 2, create and read from save file
c            fileno:      = save file number. needed only for
c                           proc_type = 2
c
c            the data structure is a vector for each element
c            block (dynamically allocated). the vectors are
c            hung from a dynamically allocated pointer vector
c            of length nelblk.
c
      proc_type = 1
c
      allocate( history_blocks(nelblk),stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 1
           call die_abort
           stop
      end if
c
      allocate( history1_blocks(nelblk),stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 2
           call die_abort
           stop
      end if
c
      allocate( history_blk_list(nelblk),stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 5
           call die_abort
           stop
      end if
c
      allocate( cep_blocks(nelblk), stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 7
           call die_abort
           stop
      end if
      allocate( cep_blk_list(nelblk),stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 8
           call die_abort
           stop
      end if
c
      history_blk_list(1:nelblk)  = 0
      cep_blk_list(1:nelblk)      = 0
c
c            loop over all element blocks. allocate blocks and
c            zero or read from restart file.
c
      do blk = 1, nelblk
c
c                       if we are using MPI:
c                         elblks(2,blk) holds which processor owns the
c                         block.  For slave processors, if we don't own
c                         the block, don't allocate it.  For root,
c                         allocate everything.
c                       if we are using the serial version:
c                         the serial process is root, so allocate all blocks.
c
c        myblk = myid .eq. elblks(2,blk)
c        if ( root_processor ) myblk = .true.
c        if ( .not. myblk ) cycle
c
        felem      = elblks(1,blk)
        span       = elblks(0,blk)
c        ngp        = iprops(6,felem)
        ngp = 1
        call material_model_info( felem, 0, 1, hist_size )
        call material_model_info( felem, 0, 2, cep_size )
c
c                      process history blocks
c                      ----------------------
c
        block_size = span * ngp * hist_size
c
        allocate( history_blocks(blk)%ptr(block_size),
     &            stat = alloc_stat )
        if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 3
           call die_abort
           stop
        end if
c
        allocate( history1_blocks(blk)%ptr(block_size),
     &            stat = alloc_stat )
        if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 4
           call die_abort
           stop
        end if
c
        if ( proc_type .eq. 1 ) then
          call ro_zero_vec( history_blocks(blk)%ptr(1), block_size )
          call ro_zero_vec( history1_blocks(blk)%ptr(1), block_size )
        end if
c
        if ( proc_type .eq. 2 ) then
c               read(fileno) history_blocks(blk)%ptr(1:block_size)
c               call chk_data_key( fileno, 200, blk )
c               call vec_ops( history1_blocks(blk)%ptr(1),
c     &                       history_blocks(blk)%ptr(1),
c     &                       dummy, block_size, 5 )
        end if
c
        history_blk_list(blk)  = hist_size
c
c
c                      process [D] (i.e. cep) blocks
c                      -----------------------------
c
        block_size = span * ngp * cep_size
        allocate( cep_blocks(blk)%vector(block_size),
     &             stat = alloc_stat )
          if ( alloc_stat .ne. 0 ) then
             write(out,9900)
             write(out,9910) 5
             call die_abort
             stop
          end if
c
        if ( proc_type .eq. 1 )
     &         call ro_zero_vec( cep_blocks(blk)%vector(1), block_size )
c
        if ( proc_type .eq. 2 ) then
c          read(fileno) cep_blocks(blk)%vector(1:block_size)
c          call chk_data_key( fileno, 200, blk )
        end if
c
        cep_blk_list(blk)  = cep_size
c
      end do
      return
c
 9900 format('>>> FATAL ERROR: memory allocate failure...')
 9910 format('                 history_cep_init: @',i2,/,
     &       '>>> Job terminated....')
c
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine strains_init                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/8/98                    *
c     *                                                              *
c     *     create the blocked data structure for storage of element *
c     *     strain data at step n and n+1. if requested, read the    *
c     *     blocks from the save file.                               *
c     *                                                              *
c     ****************************************************************
c
      subroutine strains_init()
c
      use elem_block_data, only : eps_n_blocks, eps_n1_blocks,
     &                            eps_blk_list
c
      implicit integer (a-z)
      include 'common.main'
      double precision
     &    dummy(1)
      logical myblk
c
c
c            proc_type :  = 1, create data arrays only and zero
c                         = 2, create and read from save file
c            fileno:      = save file number. needed only for
c                           proc_type = 2
c
c            the data structure is a vector for each element
c            block (dynamically allocated). the vectors are
c            hung from a dynamically allocated pointer vector
c            of length nelblk.
c
c     Ran hard code some values here
c     new input, geometric nonlinear, only one gauss point
      proc_type = 1
      ngp = 1
c
      allocate( eps_n_blocks(nelblk),stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 1
           call die_abort
           stop
      end if
c
      allocate( eps_n1_blocks(nelblk),stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 2
           call die_abort
           stop
      end if
c
      allocate( eps_blk_list(nelblk),stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 5
           call die_abort
           stop
      end if
      eps_blk_list(1:nelblk) = 0
c
c            loop over all element blocks. allocate blocks and
c            zero or read from restart file.
c
      do blk = 1, nelblk
c
c                       if we are using MPI:
c                         elblks(2,blk) holds which processor owns the
c                         block.  For slave processors, if we don't own
c                         the block, don't allocate it.  For root,
c                         allocate everything.
c                       if we are using the serial version:
c                         the serial process is root, so allocate all blocks.
c
c        myblk = myid .eq. elblks(2,blk)
c        if (root_processor) myblk = .true.
c        if (.not. myblk) cycle
c
        felem      = elblks(1,blk)
        span       = elblks(0,blk)
c        ngp        = iprops(6,felem)
        block_size = span * ngp * nstr
c
        allocate( eps_n_blocks(blk)%ptr(block_size),
     &            stat = alloc_stat )
        if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 3
           call die_abort
           stop
        endif
c
        allocate( eps_n1_blocks(blk)%ptr(block_size),
     &            stat = alloc_stat )
        if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 4
           call die_abort
           stop
        endif
c
        if ( proc_type .eq. 1 ) then
          call zero_vector( eps_n_blocks(blk)%ptr(1), block_size )
          call zero_vector( eps_n1_blocks(blk)%ptr(1), block_size )
        end if
        if ( proc_type .eq. 2 ) then
c               read(fileno) eps_n_blocks(blk)%ptr(1:block_size)
c               call chk_data_key( fileno, 360, blk )
c               call vec_ops( eps_n1_blocks(blk)%ptr(1),
c     &                       eps_n_blocks(blk)%ptr(1), dummy,
c     &                       block_size, 5 )
        end if
c
        eps_blk_list(blk) = 1
c
      end do
      return
c
 9900 format('>>> FATAL ERROR: memory allocate failure...')
 9910 format('                 strains_init: @',i2,/,
     &       '>>> Job terminated....')
c
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine stresses_init                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/8/98                    *
c     *                                                              *
c     *     create the blocked data structure for storage of element *
c     *     stress data at step n and n+1. if requested, read the    *
c     *     blocks from the save file.                               *
c     *                                                              *
c     ****************************************************************
c
      subroutine stresses_init()
c
      use elem_block_data, only : urcs_n_blocks, urcs_n1_blocks,
     &                            urcs_blk_list
c
      implicit integer (a-z)
      include 'common.main'
      double precision
     &    dummy(1)
      logical myblk
c
c
c            proc_type :  = 1, create data arrays only and zero
c                         = 2, create and read from save file
c            fileno:      = save file number. needed only for
c                           proc_type = 2
c
c            the data structure is a vector for each element
c            block (dynamically allocated). the vectors are
c            hung from a dynamically allocated pointer vector
c            of length nelblk.
c
c     Ran hard code some values here
c     new input, geometric nonlinear, only one gauss point
      proc_type = 1
      ngp = 1
c
      allocate( urcs_n_blocks(nelblk),stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 1
           call die_abort
           stop
      end if
c
      allocate( urcs_n1_blocks(nelblk),stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 2
           call die_abort
           stop
      end if
c
      allocate( urcs_blk_list(nelblk),stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 5
           call die_abort
           stop
      end if
      urcs_blk_list(1:nelblk) = 0
c
c            loop over all element blocks. allocate blocks and
c            zero or read from restart file.
c
      do blk = 1, nelblk
c
c                       if we are using MPI:
c                         elblks(2,blk) holds which processor owns the
c                         block.  For slave processors, if we don't own
c                         the block, don't allocate it.  For root,
c                         allocate everything.
c                       if we are using the serial version:
c                         the serial process is root, so allocate all blocks.
c
c        myblk = myid .eq. elblks(2,blk)
c        if (root_processor) myblk = .true.
c        if (.not. myblk) cycle
c
        felem      = elblks(1,blk)
        span       = elblks(0,blk)
c        ngp        = iprops(6,felem)
        block_size = span * ngp * nstrs
c
        allocate( urcs_n_blocks(blk)%ptr(block_size),
     &            stat = alloc_stat )
        if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 3
           call die_abort
           stop
        endif
c
        allocate( urcs_n1_blocks(blk)%ptr(block_size),
     &            stat = alloc_stat )
        if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 4
           call die_abort
           stop
        endif
c
        if ( proc_type .eq. 1 ) then
         call zero_vector( urcs_n_blocks(blk)%ptr(1), block_size )
          call zero_vector( urcs_n1_blocks(blk)%ptr(1), block_size )
        end if
        if ( proc_type .eq. 2 ) then
c               read(fileno) urcs_n_blocks(blk)%ptr(1:block_size)
c               call chk_data_key( fileno, 370, blk )
c               call vec_ops( urcs_n1_blocks(blk)%ptr(1),
c     &                       urcs_n_blocks(blk)%ptr(1), dummy,
c     &                       block_size, 5 )
        end if
c
        urcs_blk_list(blk) = 1
c
      end do
      return
c
 9900 format('>>> FATAL ERROR: memory allocate failure...')
 9910 format('                 stresses_init: @',i2,/,
     &       '>>> Job terminated....')
c
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine rotation_init                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/1/98                    *
c     *                                                              *
c     *     create the blocked data structure for storage of element *
c     *     rotation matrices at gauss points to support finite      *
c     *     strains/large displacements. allocate space for those    *
c     *     blocks which require space. initialize or read from      *
c     *     restart file                                             *
c     *                                                              *
c     ****************************************************************
c
      subroutine rotation_init()
c
      use elem_block_data, only : rot_n_blocks, rot_n1_blocks,
     &                            rot_blk_list
      implicit integer (a-z)
      include 'common.main'
      double precision
     &    zero, rot_init(9), dummy(1)
      logical geo_non_flg, local_debug
      data local_debug, zero, one / .false., 0.0, 1.0 /
      logical myblk
c
c
c            proc_type :  = 1, create data arrays only and zero
c                         = 2, create and read from save file
c            fileno:      = save file number. needed only for
c                           proc_type = 2
c
c            the data structure is a vector for each element
c            block (dynamically allocated). the vectors are
c            hung from a dynamically allocated pointer vector
c            of length nelblk.
c
c            determine if any blocks of elements have finite strains/
c            large displacements. if so, we have to create the
c            blocked data structures for material point rotation
c            matrices. for a completely small strain analysis,
c            nothing gets done here.
c
c     Ran hard code some values here
c     new input, geometric nonlinear, only one gauss point
      proc_type = 1
      ngp = 1
      geo_non_flg = .true.
c
      if ( .not. allocated(rot_blk_list) ) then
         allocate( rot_blk_list(nelblk),stat=alloc_stat )
         if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 7
           call die_abort
           stop
          end if
          rot_blk_list(1:nelblk) = 0
      end if
c
      do blk = 1, nelblk
        felem       = elblks(1,blk)
c        geo_non_flg = lprops(18,felem)
        if ( geo_non_flg ) go to 100
      end do
c
c            no blocks with geometric nonlinear elements found.
c            this is a small displacement solution. just leave.
c
      return
c
 100  continue
      allocate( rot_n_blocks(nelblk), stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 1
           call die_abort
           stop
      end if
      allocate( rot_n1_blocks(nelblk), stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 2
           call die_abort
           stop
      end if
c
c            loop over all element blocks. allocate blocks which
c            contain geometric nonlinear elements. initialize with
c            unit matrices or read from restart file.
c
      rot_init(1:9) = zero
      rot_init(1)   = one
      rot_init(5)   = one
      rot_init(9)   = one
c
c
      do blk = 1, nelblk
c
c                       if we are using MPI:
c                         elblks(2,blk) holds which processor owns the
c                         block.  For slave processors, if we don't own
c                         the block, don't allocate it.  For root,
c                         allocate everything.
c                       if we are using the serial version:
c                         the serial process is root, so allocate all blocks.
c
c        myblk = myid .eq. elblks(2,blk)
c        if (root_processor) myblk = .true.
c        if (.not. myblk) cycle
c
        felem       = elblks(1,blk)
c        geo_non_flg = lprops(18,felem)
        if ( .not. geo_non_flg ) cycle
        span       = elblks(0,blk)
c        ngp        = iprops(6,felem)
        block_size = span * ngp * 9
c
        allocate( rot_n_blocks(blk)%ptr(block_size),
     &            stat = alloc_stat )
        if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 4
           call die_abort
           stop
        endif
c
        allocate( rot_n1_blocks(blk)%ptr(block_size),
     &            stat = alloc_stat )
        if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 5
           call die_abort
           stop
        endif
c
        rot_blk_list(blk) = 1
c
        if ( proc_type .eq. 1 ) then
          low = 1
          up  = 9
          do i = 1, span
             do j = 1, ngp
                rot_n_blocks(blk)%ptr(low:up) = rot_init(1:9)
                low = low + 9
                up  = up + 9
             end do
          end do
          call vec_ops( rot_n1_blocks(blk)%ptr(1),
     &                  rot_n_blocks(blk)%ptr(1), dummy, block_size, 5 )
        end if
c
        if ( proc_type .eq. 2 ) then
c               read(fileno) rot_n1_blocks(blk)%ptr
c               call chk_data_key( fileno, 300, blk )
c               rot_n_blocks(blk)%ptr = rot_n1_blocks(blk)%ptr
        end if
c
      end do
      return
c
 9900 format('>>> FATAL ERROR: memory allocate failure...')
 9910 format('                 rotation_init: @',i2,/,
     &       '>>> Job terminated....')
c
      end
c     ****************************************************************
c     *                                                              *
c     *               subroutine ro_copy_vec                         *
c     *                                                              *
c     *                    written by : rhd                          *
c     *                                                              *
c     *                last modified : 02/27/13                      *
c     *                                                              *
c     *        copy one vector into another. use integer type        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ro_zero_vec( vec, isize )
      implicit none
c
      double precision
     &  vec(*), zero
c
      integer isize
      data zero / 0.0d00 /
c
      vec(1:isize) = zero
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine zero_vector                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 07/28/94                   *
c     *                                                              *
c     *     zero a vector of specified length w/ floating zero       *
c     *     signle or double based on this port                      *
c     *                                                              *
c     ****************************************************************
c
      subroutine zero_vector( vec, n )
      double precision
     &  vec(*), zero
      data zero / 0.0d00 /
c
!DIR$ IVDEP
      vec(1:n) = zero
c
      return
      end
     
c     ****************************************************************
c     *                                                              *
c     *                      subroutine vec_ops                      *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/27/96                   *
c     *                                                              *
c     *     perform simple vector-vector operations on single or     *
c     *     double precision vectors (based on the port)             *
c     *                                                              *
c     ****************************************************************
c
      subroutine vec_ops( veca, vecb, vecc, n, opcode )
      double precision
     &  veca(*), vecb(*), vecc(*), zero, const
      integer opcode
      data zero / 0.0d0/
c!DIR$ ASSUME_ALIGNED veca:64, vecb:64, vecc:64
c      
      go to ( 100, 200, 300, 400, 500, 600, 700 ), opcode
c
c            opcode 1:   c = a * b
c
 100  continue
!DIR$ IVDEP
      vecc(1:n) = veca(1:n) * vecb(1:n)
      return
c
c            opcode 2:   b = b * a
c
 200  continue
!DIR$ IVDEP
      vecb(1:n) = vecb(1:n) * veca(1:n)
      return
c
c            opcode 3:   c = a / b
c
 300  continue
!DIR$ IVDEP
      vecc(1:n) = veca(1:n) / vecb(1:n)
      return
c
c            opcode v:   c = zero
c
 400  continue
!DIR$ IVDEP
      vecc(1:n) = zero
      return
c
c            opcode v:   a = b
c
 500  continue
!DIR$ IVDEP
      veca(1:n) = vecb(1:n)
      return
c
c            opcode v:   a = a + b
c
 600  continue
!DIR$ IVDEP
      veca(1:n) = veca(1:n) + vecb(1:n)
      return
c
c            opcode v:   a = const * b ; const = vecc(1)
c
 700  continue
      const = vecc(1)
!DIR$ IVDEP
      veca(1:n) = const * vecb(1:n)
      return
c
      end
