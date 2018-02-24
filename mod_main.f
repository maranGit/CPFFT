c     ****************************************************************
c     *                                                              *
c     *                    f-90 module main_data                     *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *              last modified : 1/12/2017 rhd                   *
c     *                                                              *
c     *     define the data structures for main, large arrays        *
c     *     used in warp3d solutions. also other variables as we     *
c     *     gradually reduce dependence on common.main               *
c     *                                                              *
c     ****************************************************************
c
c
c
      module main_data

c      logical :: windows_os, linux_os, osx_os
c
c
c                element incidencs and incmap.
c
c
c      integer, dimension(:), save, allocatable ::  incmap, incid
c
c
c                inverse of element incidences, i.e., the list
c                of elements connected to each model node. also
c                the inverse dof maps.
c
c
c      type :: elements_on_node
c        integer element_count 
c        integer, allocatable, dimension(:) :: element_list
c      end type
c      type(elements_on_node), save, allocatable,
c     &                            dimension(:) :: inverse_incidences
c
c      type :: inv_dof_table
c        integer, allocatable, dimension(:,:) :: edof_table
c      end type
c      type(inv_dof_table), save, allocatable,
c     &                            dimension(:) :: inverse_dof_map
c
c                temporary storage of element definitions until
c                start of computations.
c
c
c      integer, dimension(:,:), save, allocatable ::  elstor
c
c
c                nodal load definitions (flags, dofs, etc.)
c
c
c      type :: node_loading_cond
c        integer node_count
c        integer how_defined
c        character(len=80) :: user_file_name
c        integer, pointer, dimension(:,:) :: nodal_loads
c      end type
c      type(node_loading_cond), save, allocatable,
c     &                         dimension(:) :: node_load_defs
c      integer, save, allocatable, dimension(:) :: temp_nodmap
c      integer, save, allocatable, dimension(:,:) :: temp_nodlod
c      integer num_loddat_blks, sizeof_loddat_blks, next_loddat_col,
c     &        max_loddat_blks
c      type :: node_load_block
c        real, pointer, dimension(:,:) :: block
c      end type
c      type(node_load_block), save, allocatable,
c     &                       dimension(:) :: loddat_blocks
c
c
c                nodal and element temperatures. current
c                totals and step increments.
c
c
c      double precision,
c     &       dimension(:), save, allocatable ::
c     &        temper_nodes, dtemp_nodes, temper_elems,
c     &        dtemp_elems, temper_nodes_ref
c      logical temperatures_ref
c
c
c                nodal constraints definitions. release constraints
c
c
c      double precision,
c     &       dimension(:), save, allocatable :: cnstrn, cnstrn_in
c     
c      type :: release_con_info
c        integer num_release_steps, remaining_steps_for_release
c      double precision
c     &   reaction_force
c      end type
c      type(release_con_info), save, allocatable,
c     &                 dimension(:,:) :: release_cons_table
c      integer :: release_cons_steps
c
c
c                nodal displacement vector for adaptive save/
c                restart
c
c
c      double precision,
c     &       dimension(:), save, allocatable :: du_nm1
c
c
c                global element-to-block mapping vectors
c                vectors.
c
c
      integer, save, allocatable, dimension(:,:) :: elems_to_blocks
c
c
c                global coordinate mapping
c
c
c      integer, save, allocatable, dimension(:) :: crdmap
c
c
c                global mapping vectors vectors.
c
c
c      integer, save, allocatable, dimension(:) ::     invdst
c      logical, save, allocatable, dimension(:) ::     repeat_incid
c
c
c                global data to store non-global constraint
c                transformations
c
c
c      type :: trn_ptr_type
c      double precision, dimension(:,:), allocatable :: mat
c      end type
c
c      type (trn_ptr_type), save, allocatable, dimension(:) :: trnmat
c      logical, save, allocatable, dimension(:) ::     trn
c
c
c                global data for diagonal mass, pbar
c
c
c      double precision,
c     &      save, allocatable, dimension(:) :: mdiag, pbar
c
c
c                global data for rloads, dloads
c
c
c      double precision,
c     &      save, allocatable, dimension(:) :: rload, dload, rload_nm1,
c     &                                         total_user_nodal_forces
c      double precision,
c     &      save, allocatable, dimension(:,:) :: load_pattern_factors
c
c                global data to store the load pattern numbers and
c                multipliers for each nonlinear load step
c
c      type :: load_data_for_a_step
c        integer :: num_load_patterns
c        integer, allocatable, dimension(:) :: load_patt_num
c        double precision,allocatable,dimension(:)::load_patt_factor
c      end type
c
c      type(load_data_for_a_step), save, allocatable,
c     &                            dimension(:) :: step_load_data
c
c
c                storage for element equivalent force vectors for
c                current load step (applied pressures, tractions,
c                body forces).
c
c
c      logical elem_equiv_loads_now
c      type :: elem_forces
c        integer :: ncols
c        double precision, pointer, dimension(:,:) :: forces
c      end type
c
c      type(elem_forces), save, allocatable,
c     &                         dimension(:) :: elem_eq_loads
c      double precision,
c     &      save, allocatable, dimension(:) :: eq_node_forces
c      integer, save, allocatable, dimension(:) :: eq_node_force_indexes
c      integer eq_node_force_len
c
c
c                 global variables being gradually moved from
c                 common.main
c
c
c      logical :: nonlocal_analysis, modified_mpcs,
c     &           divergence_check, diverge_check_strict
c
c
c                 information for output packets
c
c
c      integer packet_file_no, ascii_packet_file_no
c      logical output_packets
c      character(len=50) :: packet_file_name
c      character(len=80) :: ascii_packet_file_name
c      character(len=80) :: batch_mess_fname
c
c
c                 material properties specified at nodes to support
c                 functionally graded materials
c
c
c      real, save, allocatable, dimension(:,:) :: fgm_node_values
c      logical fgm_node_values_defined
c      integer fgm_node_values_cols
c
c
c                 logical vectors indicating element types with
c                 specific characteristics. initialized in initst.f
c
c
c      logical :: cohesive_ele_types(50),
c     &           linear_displ_ele_types(50),
c     &           adjust_constants_ele_types(50),
c     &           axisymm_ele_types(50),
c     &           implemented_ele_types(50)
c
c
c                 material properties array. these sizes correspond to
c                 mxmtpr x  mxmat in param_def and must always
c                 be consistent !
c
c      integer imatprp(300,500)
c      real  matprp(300,500)
c      logical lmtprp(300,500)
c      double precision dmatprp(300,500)
c      equivalence (matprp,lmtprp)
c      character(len=24), dimension(300,500) :: smatprp
c
c
c
c                 general tables of input for use throughout code.
c                 only "piston" input loading supported at present
c
c
c      type :: table_entry
c       character(len=24) :: table_name
c       character(len=8) ::  table_type
c       integer num_rows
c       integer num_cols
c       real, dimension(:,:), allocatable :: table_values_sgl
c       double precision, dimension(:,:),
c     &                   allocatable :: table_values_dbl
c      end type 
c
c      type (table_entry), save, allocatable,
c     &                          dimension(:) :: tables      
c
c
c                 used defined, named lists of integers.
c                 usually node numbers or element numbers
c
c
c      type :: ulist
c        character(len=24) :: name
c        integer :: length_list
c        integer, allocatable, dimension(:) :: list
c      end type
c
c      type (ulist), dimension(100) :: user_lists ! 100 is set in param_def  
c
c               UMAT model used ? Force serialization of umats?
c
c      logical :: umat_serial, umat_used
c
c
c               creep material appears in solution. will force
c               iter=0 computations
c
c      logical :: creep_model_used
c
c                 convergence information for last few load steps.
c                 used for user_solution_paramters
c
c      type :: step_convergence_data
c        logical :: step_converged  
c        logical :: adaptive_used 
c        integer :: iterations_for_convergence
c        integer :: adapt_substeps
c      end type
c      type( step_convergence_data ), dimension(5) :: 
c     &   convergence_history
c      logical run_user_solution_routine
c
c                 A CP flag, stick here b/c it's a solution parameter
c
c      logical :: cp_unloading
c
c                 Another solution parameter telling us whether 
c                 or not to use asymmetric assembly
c
c      logical :: asymmetric_assembly
c
c          file name for "output commands file ... after steps <list>'
c          bit map to store expanded list of steps
c
c      character(len=80) :: output_command_file
c      integer, save, allocatable, dimension(:) ::
c     &         output_step_bitmap_list
c
c          string names for WARP3D material models
c
c      character(len=20), save, 
c     &   allocatable, dimension(:) :: material_model_names
c
c                 do we have extrapolated displacement increments
c                 for the step and/or non-zero imposed
c                 (user) displacement increments for the step
c                 global extrapolate flag and no extrapolate next
c                 step only
c
c      logical :: extrapolated_du, non_zero_imposed_du,
c     &           extrapolate, extrap_off_next_step
c
c                 line search parameters
c
c      logical :: line_search, ls_details
c        double precision ::
c     & ls_min_step_length, ls_max_step_length, ls_rho,
c     & ls_slack_tol
c
c                 does model have crystal plasticity materials.
c                 no need to savein restart file
c                 see init.f
c
c      integer :: cp_matls_present
c
c                 options for states results on usual output cmd
c
c      integer :: output_states_type_opt1, output_states_type_opt2
c

      end module     
      




