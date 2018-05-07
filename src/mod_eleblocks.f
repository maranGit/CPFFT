c     ****************************************************************
c     *                                                              *
c     *                    f-90 data modules:                        *
c     *                     elem_block_data                          *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/9/2016 rhd               *
c     *                                                              *
c     *     define the data structures for element data stored       *
c     *     in blocks                                                *
c     *                                                              *
c     ****************************************************************
c
c
c
      module elem_block_data
c     ----------------------
c
c
c      strategies:
c
c         for older data structures
c         because some f90 compilers cannot properly build code to
c         use the "allocated" intrinsic functions on all types of
c         data types, we keep separate allocatable integer vectors
c         to indicate which blocks are actually allocated for this
c         analysis.
c
c         compilers can now handle allocatable types as rows of
c         user-defeined types
c
        type :: blocks_ptr_type
          double precision, dimension(:), pointer :: ptr
        end type
c
        type :: blocks_allocatable_type
          double precision, dimension(:), allocatable :: vector
        end type
c
c        type :: array_blocks_ptr_type
c          double precision, dimension(:,:), pointer :: ptr
c        end type
c
c        type :: int_blocks_ptr_type
c              integer, dimension(:,:), pointer :: ptr
c        end type
c
c                cdest, edest indexes for blocks
c               --------------------------------
c
c        type (int_blocks_ptr_type), save, dimension(:),
c     &         allocatable :: cdest_blocks, edest_blocks
c        integer, dimension (:), allocatable, save :: cdest_blk_list,
c     &                                               edest_blk_list
c
c                global vectors of element stress data
c               --------------------------------------
c
        type (blocks_ptr_type), save, dimension(:),
     &         allocatable :: urcs_n_blocks, urcs_n1_blocks
        integer, dimension (:), allocatable, save :: urcs_blk_list
c
c                global vectors of element strain data
c               --------------------------------------
c
        type (blocks_ptr_type), save, dimension(:),
     &         allocatable :: eps_n_blocks, eps_n1_blocks
        integer, dimension (:), allocatable, save :: eps_blk_list
c
c               element rotation matrices at Gauss points for finite
c               strains
c               ----------------------------------------------------
c
        type (blocks_ptr_type), save, dimension(:),
     &         allocatable :: rot_n_blocks, rot_n1_blocks
        integer, dimension (:), allocatable, save :: rot_blk_list
c
c               material dependent element history data in blocks
c               -------------------------------------------------
c
        type (blocks_ptr_type), save, dimension(:),
     &         allocatable :: history_blocks, history1_blocks
        integer, dimension (:), allocatable, save :: history_blk_list
c
c               material dependent [D] at each integration point in element
c               block form. internalley we call [D] "cep". 
c               store 21 terms of symmetric part. for cohesive, store
c               symmetric part of 3x3
c               ----------------------------------------------------------
c
        type (blocks_allocatable_type), save, dimension(:),
     &         allocatable :: cep_blocks
        integer, dimension (:), allocatable, save :: cep_blk_list
c
c               number of gauss points for elements in blocks
c               ---------------------------------------------
c
c        integer, dimension (:), allocatable, save :: gausspts_blk_list
c
c                global vectors of element volumes
c               ----------------------------------
c
c        type (blocks_ptr_type), save, dimension(:),
c     &         allocatable :: element_vol_blocks
c
c               element mass matrices
c               ---------------------
c
c        type (array_blocks_ptr_type), save, dimension(:),
c     &         allocatable :: mass_blocks
c
c               element stiffness
c               -----------------
c
c        type (array_blocks_ptr_type), save, dimension(:),
c     &         allocatable :: estiff_blocks
c
c               element internal forces
c               -----------------------
c
c        type (array_blocks_ptr_type), save, dimension(:),
c     &         allocatable :: einfvec_blocks
c
c               solid-interface connections for nonlocal
c               ----------------------------------------
c
c      type :: interface_solid_connections
c        integer,  allocatable, dimension(:,:) :: list
c      end type
c      type(interface_solid_connections), save, allocatable,
c     &           dimension(:) :: solid_interface_lists
c
c
c               nonlocal material state variables.
c               ---------------------------------
c
c               created for analyses with interface-cohesive elements.
c               values computed by material models for solids to be
c               shared with cohesive material models
c
c
c      logical, dimension (:), allocatable, save :: nonlocal_flags
c      type :: vec_nonlocal
c      double precision,
c     &    allocatable, dimension(:) ::state_values
c      end type
c      type(vec_nonlocal), save, allocatable,
c     &     dimension(:) :: nonlocal_data_n, nonlocal_data_n1
c
c        intrinsic allocated
c
      end module
