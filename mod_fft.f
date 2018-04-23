      module fft
      
c     ========================== problem size ==========================
      integer, save :: ndim1, ndim2, ndim3, ndim4
      integer, save :: N, N3, Nhalf
      integer, save :: veclen
      integer, save :: dims(3)
      real(8), save :: l_x, l_y, l_z

c     ====================== material parameters ======================
c
c                 material properties array. these sizes correspond to
c                 mxmtpr x  mxmat in param_def and must always
c                 be consistent !
c
      type :: mat_property
        logical :: assigned
        integer :: matnum
        character(len=24) :: matnam
        integer :: imatprp(300)
        logical :: lmatprp(300)
        double precision :: dmatprp(300)
        character(len=24), dimension(300) :: smatprp
      end type

      type (mat_property), save :: mat_props(500)

      character(len=20), save, 
     &   allocatable, dimension(:) :: material_model_names

c     ================== Newton-Raphson loop control ==================
      real(8), save :: straininc, tolPCG, tolNR
      integer, save :: maxIter, nstep

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
