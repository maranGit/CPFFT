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
        myblk = myid .eq. elblks(2,blk)
        if ( root_processor ) myblk = .true.
        if ( .not. myblk ) cycle
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
