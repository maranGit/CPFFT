      module fft
      
c     ========================== problem size ==========================
      integer, save :: ndim1, ndim2, ndim3, ndim4
      integer, save :: N, N3, Nhalf
      integer, save :: veclen
      integer, save :: dims(3)
      real(8), save :: l_x, l_y, l_z
      integer, dimension(:), save, allocatable ::  incmap, incid
c
c     ================== material properties array ====================
c             these sizes correspond to mxmtpr x mxmat in param_def,
c             and must always be consistent!
c
      character(len=20), save, 
     &   allocatable, dimension(:) :: material_model_names

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
c     ===================== boundary condition ========================
c
      real(8), save :: BC_all(9,100000)
      real(8), save :: FP_max(9), mults(100000) ! mxstep from param_def
      logical, save :: isNBC(9)
c
c     ================== Arrays to convert F to u =====================
c
      real(8), allocatable :: BTB(:)
      integer, allocatable :: ia_btb(:), ja_btb(:)
c
c     ================== Newton-Raphson loop control ==================
      real(8), save :: straininc, tolPCG, tolNR
      integer, save :: maxIter, nstep
      real(8), save :: tstep
      logical, save :: out_step(100000)
c
c                 Another solution parameter telling us whether
c                 or not to use asymmetric assembly
c
      logical, save :: asymmetric_assembly

c     ===================== global state variables =====================
      real(8), dimension(:,:), save, allocatable :: Ghat4, K4
      real(8), dimension(:,:), save, allocatable :: b, dFm
      real(8), dimension(:,:), save, allocatable :: Fn1, Fn, Pn1, Pn
      integer, dimension(:), save, allocatable :: matList

c     ================ fft related temporary variables ================
c     real(8), allocatable :: real1(:)
c     complex(8), allocatable :: cplx3half(:,:,:), cplx1half(:)
c
c            form coefficients for fftshift and ifftshift
c
c          coeffs1 is the real part of e^(-i*x*h/2)
c          coeffs2 is the imaginary part of e^(-i*x*h/2)
c          for forward transformation, use (coeffs1, coeffs2)
c          for backward transformation, use (coeffs1, -coeffs2)
c
      real(8), allocatable :: coeffs1(:), coeffs2(:)

c     ================ pcg related temporary variables ================
      real(8), allocatable :: tmpPcg( :, : )

c     =============== G_K_dF related temporary variables ===============
      real(8), dimension(:,:), allocatable :: real1, real2, real3

      end module
