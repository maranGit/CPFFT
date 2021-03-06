c     ****************************************************************
c     *                                                              *
c     *                  subroutine  FFT_init                        *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified : 4/5/2018 RM                *
c     *                                                              *
c     *      Initialize the whole project                            *
c     *      Including:                                              *
c     *          arrays in fft.mod, elem_block_data.mod              *
c     *          variables in common.main                            *
c     *                                                              *
c     ****************************************************************
c
      subroutine FFT_init()
      use fft
      use file_info
      implicit none
      include 'common.main'

c                         local
      logical :: debug
      real(8), parameter :: zero = 0.0D0, one = 1.0D0
      integer :: nblank, reclen, endchr
      logical promsw,echosw,comsw,atrdsw,eolsw,eofsw,menusw,ptsw,signsw
      
      integer :: i
c
c                       initialize common.main
c
      use_mpi = .false.
      myid = 0
      ltmstp = 0
c
c                       initialize the file input and output parameters
c
      inlun(1) = 5
      inlun(2) = 80
      inlun(3) = 81
      inlun(4) = 82
      inlun(5) = 83
      inlun(6) = 84
      inlun(7) = 85
      inlun(8) = 86
      inlun(9) = 87
      inlun(10) = 88
c
      outlun(1) = 6
      outlun(2) = 89
c
      filcnt   = 1
      in       = inlun(filcnt)
      out      = outlun(1)
      outing   = .false.
c
c      output_packets = .false.
c      packet_file_no = 97
c      ascii_packet_file_no = 96
c      packet_file_name(1:) = ' '
c      ascii_packet_file_name(1:) = ' '
c      batch_mess_fname(1:) = ' '
c

c                       summary of fortran file numbers used in warp3d
c                       add new ones here....
c
c                       code updates are gradually using the function
c                       warp3d_get_device_number() to obtain an available
c                       file number for use during well-defined block
c                       of execution
c
c         terminal input:               5  (Unix standard input)
c         terminal output:              6  (Unix standard output)
c         ascii energy file:            11
c         other ascii input files
c           for *input from file:       80-88
c         other ascii output files
c           for *output to file:        89
c         ascii data packet file:       96
c         binary data packet file:      97
c         binary results file:          98
c         formatted results file:       99

c
c
c                       initialize scan
c
      call setin(in)
      call setout(out)
      nblank= 80
      reclen= 80
      endchr= 1h$
      promsw= .false.
      echosw= .true.
      comsw= .false.
      atrdsw= .false.
      eolsw= .true.
      eofsw= .true.
      menusw= .false.
      ptsw= .false.
      signsw= .false.
      call check_to_prompt( promsw )
      call scinit(nblank,reclen,endchr,promsw,echosw,comsw,atrdsw,
     &            eolsw,eofsw,menusw,ptsw,signsw)
c
c                       Asymmetric assembly
c
      asymmetric_assembly = .false.
c
c
c                    build-in material library
c
      allocate( material_model_names(mxmat) )
      material_model_names(1:mxmat)(1:) = "not_used"
      material_model_names(1)(1:) = "bilinear"     
      material_model_names(2)(1:) = "deformation"
      material_model_names(3)(1:) = "mises_gurson"     
      material_model_names(4)(1:) = "cohesive"     
      material_model_names(5)(1:) = "cyclic"     
      material_model_names(6)(1:) = "creep"     
      material_model_names(7)(1:) = "mises_hydrogen"     
      material_model_names(8)(1:) = "umat"     
      material_model_names(9)(1:) = "not_used"     
      material_model_names(10)(1:) = "crystal_plasticity"     
      material_model_names(11)(1:) = "interface_damage"
c
      cp_matls_present = -1
      nummat = 0
c
c                       initialize timing parameters.
c
      times(1:mxtim,1) = zero
      times(1:mxtim,2) = zero
c
c     hard code NR loop parameters in fft.mod
c
      debug = .false.

      ndim1 = 3
      ndim2 = ndim1 * ndim1
      ndim3 = ndim2 * ndim1
      ndim4 = ndim3 * ndim1
      N3 = N * N * N
      noelem = N3
      if ( mod(N, 2) .eq. 1 ) Nhalf = (N + 1) / 2
      if ( mod(N, 2) .eq. 0 ) Nhalf = N / 2 + 1
      veclen = N3 * 9
      dims = [N, N, N]
      noelem = N3
c
c     allocate variables in module
c
      call fftAllocate( 1 )
c
c     initialize deformation gradient
      Fn = zero
      Fn(:,[1,5,9]) = one
      call dcopy(N3*nstrs, Fn(1,1), 1, Fn1(1,1), 1)
      Pn = zero
      Pn1 = zero
c
c     initialize BC
c
      isNBC(1:9) = .false.
      FP_max = zero
      BC_all = zero
c
c     form G_hat_4 matrix and store in Ghat4
c
      call formG()
c
c     form coefficients for fftshift and ifftshift
c
      call formfftshift(coeffs1, coeffs2)
c
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine      fftAllocate                  *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 2/3/18                      *
c     *                                                              *
c     *                  allocate array in mod_fft.f                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine fftAllocate( action )
      use fft
      use elem_block_data
      implicit none
      include 'common.main'
      
c     temporary variables
      integer :: blk
      integer :: action, iok, N3ndim2
      
      N3ndim2 = N3 * ndim2
      
      select case ( action )
      case( 1 )
c
c     allocate variables in module
c
        allocate ( Ghat4(N3, ndim4), stat=iok )
        allocate ( K4(N3, ndim4), stat=iok )
        allocate ( Fn(N3, ndim2), stat=iok )
        allocate ( Fn1(N3, ndim2), stat=iok )
        allocate ( Pn(N3, ndim2) )
        allocate ( Pn1(N3, ndim2) )
        allocate ( b(N3, ndim2) )
        allocate ( dFm(N3, ndim2) )
        allocate ( matList(N3) )
  
c     allocate FFT related variables
        allocate ( coeffs1(N3), coeffs2(N3) )
  
c     allocate pcg related variables
        allocate ( tmpPcg(N3ndim2, 4) )
  
c     allocate internal variables
        allocate( real1(N3, ndim2) )
        allocate( real2(N3, ndim2) )
        allocate( real3(N3, ndim2) )

      case ( 3 )
c     allocate global arrays in Warp3d
        call history_cep_init()
        call stresses_init()
        call rotation_init()
        call strains_init()
c
      case( 2 )
c
c     deallocate variables in module
c
        deallocate( Ghat4, K4, Fn, Fn1 )
        deallocate( Pn, Pn1, b, dFm, matList )
        deallocate( coeffs1, coeffs2 )
        deallocate( tmpPcg )
        deallocate( real1, real2, real3 )
        deallocate( history_blk_list )
        do blk = 1, nelblk
          deallocate( cep_blocks(blk)%vector )
        end do
        deallocate( cep_blocks )
        deallocate( incid, incmap )
        if(allocated(BTB)) deallocate(BTB,ia_btb,ja_btb)

      case default
        write(*,*) ">>Error: invalid option. Job terminated"
      end select

      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine       formG                       *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 1/22/18                     *
c     *                                                              *
c     *                  form G_hat_4 matrix                         *
c     *                                                              *
c     ****************************************************************
c
      subroutine formG()
      use fft, only: N, Nhalf, Ghat4, N3, ndim1, ndim2
      implicit none

c     internal variables
      integer :: ii, jj, kk, ll, tmp
      real(8), dimension(:,:), allocatable :: q
      real(8), dimension(:), allocatable :: q_dot_q
      integer :: indices(4,81)

c                         indices array
      indices = 0
      indices(1,1:27) = 1
      indices(1,28:54) = 2
      indices(1,55:81) = 3
c
      indices(2,1:9) = 1
      indices(2,10:18) = 2
      indices(2,19:27) = 3
      indices(2,28:54) = indices(2,1:27)
      indices(2,55:81) = indices(2,1:27)
c
      indices(3,1:81) = 1
      indices(3,4:81:9) = 2
      indices(3,5:81:9) = 2
      indices(3,6:81:9) = 2
      indices(3,7:81:9) = 3
      indices(3,8:81:9) = 3
      indices(3,9:81:9) = 3
c
      indices(4,1:81:3) = 1
      indices(4,2:81:3) = 2
      indices(4,3:81:3) = 3

c     allocate internal variables
      allocate ( q(N3, ndim1) )
      allocate ( q_dot_q(N3) )

      tmp = 1
      do ii = 1, N
        do jj = 1, N
          do kk = 1, N
            q(tmp, :) = dble( [ ii-Nhalf, jj-Nhalf, kk-Nhalf ] )
            tmp = tmp + 1
          enddo
        enddo
      enddo
      q_dot_q = q(:,1)*q(:,1) + q(:,2)*q(:,2) + q(:,3)*q(:,3)

      do tmp = 1, 81
        ii = indices(1,tmp)
        jj = indices(2,tmp)
        kk = indices(3,tmp)
        ll = indices(4,tmp)
        if (ii .eq. kk) Ghat4(:,tmp) = q(:,jj) * q(:,ll)
      end do

      do ii = 1, N3
        if ( abs( q_dot_q(ii) ) .le. 1D-10 ) then
          Ghat4(ii, :) = 0.0D0
        else
          Ghat4(ii, :) = Ghat4(ii, :) / q_dot_q(ii)
        endif
      enddo

c     deallocate internal variables
      deallocate( q, q_dot_q )
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine   formfftshift                    *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 7/7/18                      *
c     *                                                              *
c     *               form coefficient for fftshift                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine formfftshift( a, b )
      use fft, only: N, Nhalf
      implicit none
c
c                         global
c
      real(8) :: a(N,N,*), b(N,N,*)
c
c                         local
c
      integer :: i, j, k
      real(8) :: temp, NhN, cosx, sinx
      real(8), parameter :: zero = 0.0D0, one = 1.0D0
      real(8), parameter :: two = 2.0D0, four = 4.0D0
      real(8), parameter :: PI = four * atan (one)
c
      NhN = dble(Nhalf) * two * PI / dble(N)
c
      do k = 1, N
        do j = 1, N
          do i = 1, N
            temp = NhN * dble( i + j + k - 3 )
            cosx = dcos( temp )
            sinx = dsin( temp )
            a(i,j,k) = cosx
            b(i,j,k) = -sinx
          end do
        end do
      end do
c
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine      inelbk                       *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 1/22/18                     *
c     *                                                              *
c     *     automatic generation of element blocking: simple version *
c     *     for threaded w/o vectorized blocking. optional           *
c     *     assignment of blocks to domains w/ simple algorithm      *
c     *                                                              *
c     ****************************************************************
c
      subroutine inelbk( matList )
      implicit none
      include 'common.main'
c
      integer :: auto_size
      logical :: display
      integer :: felem, current_size, element, param, i
      logical :: newblk, compatible
      integer :: blk_matmodel, ele_matmodel
      integer, intent(in) :: matList(noelem)
c
c     Ran: hard coded variables for autoblock
c
      auto_size    = mxvl  !  from param_def
      display      = .false.
c
c                     first generation of automatic assignment of
c                     elements to blocks.
c
c                     a) sequential pass thru elements
c                     b) assign to current block or open new
c                        block if full or element/material
c                        combinations are not compatible with
c                        current block.
c                     c) no element renumbering or vectorized
c                        blocking in this first set of
c                        features.
c
c
      nelblk       = 1   ! in common main
      current_size = 1
      felem        = 1
      elblks(1,1)  = 1  ! first element in block
      elblks(0,1)  = 1  ! number elements in block
      if( noelem == 1 ) return
c
      blk_matmodel = matList( 1 )
c
c                     1. element fits in current block?
c                     2. if yes, is it compatible with elements now
c                        in the block?
c                     3. if yes, update current number of elements
c                        in the block, next element
c                     4. otherwise start a new block. set first
c                        element in the block, init block size,
c                        load props for first element in block for
c                        subsequent comparisons
c
      do element = 2, noelem
         newblk       = .false.
         current_size = current_size + 1
         if( current_size .gt. auto_size ) newblk = .true.
         if( .not. newblk ) then
           ele_matmodel = matList( element )
           compatible = ( blk_matmodel .eq. ele_matmodel )
           if( .not. compatible ) newblk = .true.
         end if
         if( .not. newblk ) then
           elblks(0,nelblk)  = current_size
           cycle
         endif
         nelblk = nelblk + 1
         if( nelblk .gt. mxnmbl ) then
            param = nelblk
c      call errmsg(74,param,dums,dumr,dumd)
c      call die_abort
            write(out,*) "too many element blocks required"
            stop
         end if
         felem             = element
         elblks(1,nelblk)  = felem
         elblks(0,nelblk)  = 1
         current_size       = 1
         blk_matmodel = matList( felem )
      end do ! on element
c
c                     display blocking table if requested
c
      if( .not. display ) return
      write(out,9000) nelblk, auto_size
      write(out,9010)
      do i = 1, nelblk
        write(out,9020) i, elblks(1,i),  elblks(0,i), elblks(2,i)
      end do
      write(out,*) ' '
c
      return
c
 9000 format(/,'>> Generated element blocking table:',
     & /,      1x,'  number of blocks, target size: ',i7,i5)
 9010 format(/,5x,
     &'block        1st element in blk        # elements in block',
     &5x,'assigned domain' )
 9020 format(1x,i8, 10x,i9,25x,i4,10x,i6)

      end subroutine
