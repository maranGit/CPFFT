      module fft
      
c     ========================== problem size ==========================
      integer, save :: ndim1, ndim2, ndim3, ndim4
      integer, save :: N, N3, Nhalf
      integer, save :: veclen
      integer, save :: dims(3)
      real(8), save :: l_x, l_y, l_z
      integer, dimension(:), save, allocatable ::  incmap, incid

c     ====================== material parameters ======================
c
c                 material properties array. these sizes correspond to
c                 mxmtpr x  mxmat in param_def and must always
c                 be consistent !
c
c     type :: mat_property
c       logical :: assigned
c       integer :: matnum
c       character(len=24) :: matnam
c       real    :: matprp(300)
c       integer :: imatprp(300)
c       logical :: lmtprp(300)
c       equivalence (matprp,lmtprp)
c       equivalence (matprp,imatprp)
c       double precision :: dmatprp(300)
c       character(len=24), dimension(300) :: smatprp
c     end type

c     type (mat_property), save :: mat_props(500)

      character(len=20), save, 
     &   allocatable, dimension(:) :: material_model_names

c
c
c                 material properties array. these sizes correspond to
c                 mxmtpr x  mxmat in param_def and must always
c                 be consistent !
c
      logical          :: mat_assigned(500) ! if assigned
      integer          :: imatprp(300,500)
      real             :: matprp(300,500)
      logical          :: lmtprp(300,500)
      equivalence (matprp,lmtprp)
      equivalence (matprp,imatprp)
      double precision :: dmatprp(300,500)
      character(len=24), dimension(300,500) :: smatprp
      integer          :: cp_matls_present
c
c
c     ================== Newton-Raphson loop control ==================
      real(8), save :: straininc, tolPCG, tolNR
      integer, save :: maxIter, nstep, tstep
      real(8), save :: F_total(9), mults(100000) ! mxstep from param_def
      logical, save :: out_step(100000)
c
c                 Another solution parameter telling us whether
c                 or not to use asymmetric assembly
c
      logical, save :: asymmetric_assembly

c     ===================== global state variables =====================
      real(8), dimension(:,:), save, allocatable :: Ghat4, K4
      real(8), dimension(:,:), save, allocatable :: DbarF, b, dFm
      real(8), dimension(:,:), save, allocatable :: barF, barF_t
      real(8), dimension(:,:), save, allocatable :: Fn1, Fn, Pn1, Pn
      integer, dimension(:), save, allocatable :: matList

c     ================ fft related temporary variables ================
      real(8), allocatable :: real1(:)
      complex(8), allocatable :: cplx3half(:,:,:), cplx1half(:)

c     ================ pcg related temporary variables ================
      real(8), allocatable :: tmpPcg( :, : )

c     =============== G_K_dF related temporary variables ===============
      real(8), dimension(:,:), allocatable :: tmpReal
      complex(8), dimension(:,:), allocatable :: tmpCplx

      end module
