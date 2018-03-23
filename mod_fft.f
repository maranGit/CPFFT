      module fft
c     problem size
      integer, save :: ndim1, ndim2, ndim3, ndim4
      integer, save :: N, N3, Nhalf
      integer, save :: veclen
      integer, save :: dims(3)

c     global state variables
      real(8), dimension(:,:), save, allocatable :: Ghat4, K4
      real(8), dimension(:,:), save, allocatable :: DbarF, b, dFm
      real(8), dimension(:,:), save, allocatable :: barF, barF_t
      real(8), dimension(:,:), save, allocatable :: Fn1, Fn, Pn1, Pn
      integer, dimension(:), save, allocatable :: matList

c     fft related temporary variables
      real(8), allocatable :: real1(:)
      complex(8), allocatable :: cplx3half(:,:,:), cplx1half(:)

c     pcg related temporary variables
      real(8), allocatable :: tmpPcg( :, : )

c     G_K_dF related temporary variables
      real(8), dimension(:,:), allocatable :: tmpReal
      complex(8), dimension(:,:), allocatable :: tmpCplx

      end module
