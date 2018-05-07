c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine ModelOut                     *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   write patran input file                    *
c     *                                                              *
c     ****************************************************************
      subroutine ModelOut()
      use fft, only: N, N3, l_x, l_y, l_z
      implicit none
      include 'common.main'
c
c                dummy variables
c
      integer :: dum, dummy
      character :: dums*8
      real :: dumr
      real(8) :: dumd
c
c                parameters for mesh generation part
c
      logical :: debug
      real(8), parameter :: zero = 0.0D0, one = 1.0D0
      integer :: ndofs, ndofs_1
      real(8), allocatable :: x(:,:)
      integer, allocatable :: ix(:,:)
      integer :: rinc,sinc,tinc
      real(8) :: xl(3,8)
      integer :: ixl(8),numnp,numel,nen,nen1
      integer :: i, ii, ll
c
c                  write patran flat file
c
      logical, external :: label
      integer :: nc
      character :: flat_name*80, tmpname*80
      logical :: text_file, stream_file, compressed
      logical :: patran_convention, warp3d_convention
c
c                  global model size
c
      data ixl /1,2,4,3,5,6,8,7/

      xl(1:3, 1:8) = zero

      xl(1,2) = l_x
      xl(1,4) = l_x
      xl(1,6) = l_x
      xl(1,8) = l_x

      xl(2,3) = l_y
      xl(2,4) = l_y
      xl(2,7) = l_y
      xl(2,8) = l_y

      xl(3,5) = l_z
      xl(3,6) = l_z
      xl(3,7) = l_z
      xl(3,8) = l_z
c
c                 read model file name
c
      if(label(dummy)) then
        tmpname= ' '
        call entits(tmpname,nc)
        if(nc.gt.20) nc=20
        flat_name(1:nc)= tmpname(1:nc)
      else
        call errmsg(19,dum,dums,dumr,dumd)
        flat_name(1:14) = 'FFT_model_flat'
      endif
c
c                 model size (reference configuration)
c
      debug = .true.
      ndofs = ( N + 1 ) * ( N + 1 ) * ( N + 1 )
      ndofs_1 = ndofs - 1
      numnp = ndofs
      numel = N3
      nen = 8
      nen1 = 9
      allocate( x(3,numnp), ix(nen1,numel) )
c
c                   get coordinates of initial configuration
c
      rinc = N
      sinc = N
      tinc = N
      call blkgen(rinc,sinc,tinc,xl,ixl,x,ix,numnp,numel)
c
c                      store in common.main
c
      call CommonStore(numnp,numel,nen,x(1,1),ix(1,1))
c
c                      write patran flat file
c
      text_file         = .false.
      stream_file       = .false.
      compressed        = .false.
      patran_convention = .false.
      warp3d_convention = .false.

      text_file         = .true.
      patran_convention = .true.

      call oumodel_flat( flat_name, text_file, stream_file, compressed, 
     &                   warp3d_convention, patran_convention  )

      return
      end subroutine  
c
c     ****************************************************************
c     *                                                              *
c     *                    subroutine    CommonStore                 *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 5/2/18                      *
c     *                                                              *
c     ****************************************************************
c
      subroutine CommonStore(numnp,numel,nen,x,ix)
      use fft, only: incmap, incid
      implicit none
      include 'common.main'

c                              global
      integer :: numnp, numel, nen
      integer :: ix(nen+1,*)
      real(8) :: x(*)

c                              local
      integer :: i, j, nen1, elemptr

      nonode = numnp
      noelem = numel
      nen1 = nen + 1

      allocate( incmap(numel+1), incid(numel*nen) )

      do i = 1, 3*numnp
        c(i) = x(i)
      end do

      elemptr = 1
      do i = 1, numel
        incmap(i) = elemptr
        do j = 1, nen
          incid(elemptr) = ix(j,i)
          elemptr = elemptr + 1
        end do
      end do

      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                    subroutine    blkgen                      *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 5/2/18                      *
c     *                                                              *
c     ****************************************************************
c
      subroutine blkgen(rinc,sinc,tinc,xl,ixl, ! input
     &                  x,ix,numnp,numel)      ! output

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Delete label 100                                 13/12/2006
c       2. Correct 11 node tet node numbers                 07/02/2007
c       3. Add error on too few incrments for higher order  05/03/2007
c          quadratic and cubic elements.
c       4. Add option for 14 and 15 node tet elements       30/08/2007
c       5. Add 'prt' and 'prth' to 'vblke' call             06/09/2007
c       6. Add 'cap' option for spherical cap               11/11/2008
c          and 'xcyl' and 'ycyl' for cylindrical sectors
c       7. Increase number of reads from 4 to 5             19/11/2008
c       8. Add counts for 64-node brick type elements       06/02/2009
c       9. Correct check on 'ntyp' for 16 node quads        28/05/2009
c      10. Correct call order to ck3dblk                    06/07/2010
c      11. Add 'mate' to set of material number             30/11/2010
c      12. Add set of last element number to last_elm       29/01/2012
c      13. Add surface option for generation in 3-d         31/01/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Generate a block of elements and nodes for mesh
c               descriptions.

c      Inputs:
c         ndm        - Dimension of 'x' array
c         nen1       - Dimension for 'ix' array
c         prt        - Print generated data if true
c         prth       - Print headers if true

c      Outputs:
c         x(ndm,*)   - Block of nodal coordinates
c         ix(nen1,*) - Block of elements
c         rben(numel) - Rigid body number associated with element
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

c                hard-coded constant
      integer, parameter :: ndm = 3, ndim1 = 4, nen = 8, nen1 = 9
      integer, parameter :: ntyp = 10, numGlobalNode = 8
      character, parameter :: ctype = 'cart'
      integer, parameter :: node1 = 1, elmt1 = 1, mat = 1

c                     global
      integer :: rinc, sinc, tinc ,numnp, numel
      real(8) :: x(ndm,*), xl(ndm,*)
      integer :: ix(nen1,*), ixl(*)

c                     local
      integer :: x0(3)
      real(8), parameter :: zero = 0.0D0, one = 1.0D0
      integer :: nr, ns, nt, ni, ne, ma, nf, ng, i, ind
      real(8) :: dr, ds, dt
      real(8) :: xll(3,27)
      integer :: ixll(27)
      logical :: local_debug
c     logical   prt,prth,pcomp,errck,pinput,tinput,palloc,bernfl,err
c     logical   eltype,pblktyp
c     character xh*6,layer*15
c     integer   i,j,k,l,n,nm,nn,nr,ns,nt,nf,ng,ni,ntyp,nodinc
c     integer   ndm,nen1,ne,ma,mab, dlayer,nlay
c     integer   ixl(27),rben(*)
c     real*8    dr,ds,dt
c     real*8    shp(3,16),tb(7),tc(4),td(16)

c     save

c     data      xh/' coord'/
      data      x0 / 0, 0, 0 /
      data      local_debug / .false. /

c     Set parameters for block

      do i = 1, numGlobalNode
        ind = ixl(i)
        ixll(ind) = ind
        xll(1,ind) = xl(1,i)
        xll(2,ind) = xl(2,i)
        xll(3,ind) = xl(3,i)
      end do

      nr     = rinc
      ns     = sinc
      nt     = tinc

c     3-d generations

      ni     = node1
      ne     = elmt1
      ma     = mat
c     nodinc = 0

c     Set generation increments of natural coordinates

      dr = 2.d0/nr
      ds = 2.d0/ns
      dt = 2.d0/nt

      if(ntyp.eq.10) then                     ! 8-node hexahedron
        nf = ne + nr*ns*nt - 1
      else
c       write(ilg,4003) ntyp
c       write(iow,4003) ntyp
      endif

      numel = nf
      nr = nr + 1
      ns = ns + 1
      nt = nt + 1
      ng = nr*ns*nt + ni - 1
      numnp = ng

c       Compute node locations
      if ( local_debug ) then
        write(*,*) "xll"
        do i = 1, 8
          write(*,*) xll(1:3,i)
        end do
      end if

      call vblkn(nr,ns,nt,xll,x,ixll,dr,ds,dt,
     &           ni,ndm,ctype,numnp,x0)

      if ( local_debug ) then
        write(*,*) "coordinate"
        do i = 1, 27
          write(*,*) x(1:3,i)
        end do
      end if

c       Compute element connections

      if(ne.gt.0) then
        call vblke(nr,ns,nt,ix,ni,ne,nen1,ma,ntyp,numel)
      endif

      return

 4003 format(1x,'>>>>> Error: Invalid element type: ',i5,
     &          ' in blkgen',/)

      end
c
c     ****************************************************************
c     *                                                              *
c     *                    subroutine    bcor3d                      *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 5/2/18                      *
c     *                                                              *
c     ****************************************************************
c
      subroutine bcor3d(ixl,xl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Revise block numbering to be same as element     10/11/2008
c          numbers for 27-node Lagrange type.
c       2. Restore option for old numbering on flag oldfl   02/02/2009
c       3. Correct numbering error in emid and assignment   03/02/2009
c          of nonzero locations in ixl array
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute missing coordinate values of 27-node element

c      Inputs:
c         ixl(*)    - Nodal connection list
c         xl(3,*)   - Unadjusted coordinate array

c      Outputs:
c         xl(3,*)   - Adjusted coordinate array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

c     include   'corset.h'

      integer    ixl(27), imid(12),amid(12),bmid(12),cmid(6)
c     integer    dmi1(6),dmi2(6),dmi3(6),dmi4(6)
c     integer    dmi5(6),dmi6(6),dmi7(6),dmi8(6)
c     integer    fmi5(6),fmi6(6),fmi7(6),fmi8(6), emid(6),jmid(12)
      real*8     xl(3,27)

      integer    i,j

      save

c     New

      data       imid/ 9,10,11,12, 13,14,15,16, 18,19,20,21/
      data       amid/ 1, 2, 3, 4,  1, 2, 3, 4,  5, 6, 7, 8/
      data       bmid/ 5, 6, 7, 8,  2, 3, 4, 1,  6, 7, 8, 5/
c     data       cmid/21,22,23,24,25,26/
c     data       dmi1/ 1, 2, 3, 1, 1, 5/
c     data       dmi2/ 4, 3, 4, 2, 2, 6/
c     data       dmi3/ 8, 7, 8, 6, 3, 7/
c     data       dmi4/ 5, 6, 7, 5, 4, 8/
c     data       dmi5/12,10,11, 9, 9,13/
c     data       dmi6/16,14,15,13,10,14/
c     data       dmi7/17,18,19,17,11,15/
c     data       dmi8/20,19,20,18,12,16/

c     Old

c     data       jmid/13,14,15,16, 18,19,20,21,  9,10,11,12/
c     data       emid/26,24,25,23,17,22/
c     data       fmi5/16,14,15,13,13,18/
c     data       fmi6/21,19,20,18,14,19/
c     data       fmi7/ 9,10,11, 9,15,20/
c     data       fmi8/12,11,12,10,16,21/

c     Numbering in 'old' form

c     if(oldfl) then

c       Mid edge coordinates

c       do i = 1,12
c         if(ixl(jmid(i)).eq.0) then
c           do j = 1,3
c             xl(j,jmid(i)) = 0.5d0*(xl(j,amid(i)) + xl(j,bmid(i)))
c           end do ! j
c           ixl(jmid(i)) = jmid(i)
c         endif
c       end do ! i

c       Center face nodes

c       do i = 1,6
c         if(ixl(emid(i)).eq.0) then
c           do j = 1,3
c             xl(j,emid(i)) = 0.50d0*(xl(j,fmi5(i)) + xl(j,fmi6(i))
c    &                              + xl(j,fmi7(i)) + xl(j,fmi8(i)))
c    &                      - 0.25d0*(xl(j,dmi1(i)) + xl(j,dmi2(i))
c    &                              + xl(j,dmi3(i)) + xl(j,dmi4(i)))
c           end do ! j
c           ixl(emid(i)) = emid(i)
c         endif
c       end do ! i

c       Center node

c       if(ixl(27).eq.0) then
c         do j = 1,3
c           xl(j,27) = 0.125d0*(xl(j, 1) +xl(j, 2) +xl(j, 3) +xl(j, 4)
c    &                        + xl(j, 5) +xl(j, 6) +xl(j, 7) +xl(j, 8))
c    &               - 0.250d0*(xl(j, 9) +xl(j,10) +xl(j,11) +xl(j,12)
c    &                        + xl(j,13) +xl(j,14) +xl(j,15) +xl(j,16)
c    &                        + xl(j,18) +xl(j,19) +xl(j,20) +xl(j,21))
c    &               + 0.500d0*(xl(j,17) +xl(j,22) +xl(j,23) +xl(j,24)
c    &                        + xl(j,25) +xl(j,26))
c         end do ! j
c         ixl(27) = 27
c       endif

c     Current numbering order on elements

c     else

c       Mid edge coordinates

        do i = 1,12
          if(ixl(imid(i)).eq.0) then
            do j = 1,3
              xl(j,imid(i)) = 0.5d0*(xl(j,amid(i)) + xl(j,bmid(i)))
            end do ! j
            ixl(imid(i)) = imid(i)
          endif
        end do ! i

c       Center face nodes

c       do i = 1,6
c         if(ixl(cmid(i)).eq.0) then
c           do j = 1,3
c             xl(j,cmid(i)) = 0.50d0*(xl(j,dmi5(i)) + xl(j,dmi6(i))
c    &                              + xl(j,dmi7(i)) + xl(j,dmi8(i)))
c    &                      - 0.25d0*(xl(j,dmi1(i)) + xl(j,dmi2(i))
c    &                              + xl(j,dmi3(i)) + xl(j,dmi4(i)))
c           end do ! j
c           ixl(cmid(i)) = cmid(i)
c         endif
c       end do ! i

c     Bottom and top

      if(ixl(17) .eq. 0) then
        do j = 1, 3
          xl(j,17) = 0.50d0*(xl(j,13) + xl(j,14) + xl(j,15) + xl(j,16))
     &             - 0.25d0*(xl(j, 1) + xl(j, 2) + xl(j, 3) + xl(j, 4))
        end do
        ixl(17) = 17
      end if

      if(ixl(22) .eq. 0) then
        do j = 1, 3
          xl(j,22) = 0.50d0*(xl(j,18) + xl(j,19) + xl(j,20) + xl(j,21))
     &             - 0.25d0*(xl(j, 5) + xl(j, 6) + xl(j, 7) + xl(j, 8))
        end do
        ixl(22) = 22
      end if

c     Mid-face

      if(ixl(23) .eq. 0) then
        do j = 1, 3
          xl(j,23) = 0.50d0*(xl(j,13) + xl(j, 9) + xl(j,10) + xl(j,18))
     &             - 0.25d0*(xl(j, 1) + xl(j, 2) + xl(j, 5) + xl(j, 6))
        end do
        ixl(23) = 23
      end if

      if(ixl(24) .eq. 0) then
        do j = 1, 3
          xl(j,24) = 0.50d0*(xl(j,14) + xl(j,10) + xl(j,11) + xl(j,19))
     &             - 0.25d0*(xl(j, 2) + xl(j, 3) + xl(j, 6) + xl(j, 7))
        end do
        ixl(24) = 24
      end if

      if(ixl(25) .eq. 0) then
        do j = 1, 3
          xl(j,25) = 0.50d0*(xl(j,15) + xl(j,11) + xl(j,12) + xl(j,20))
     &             - 0.25d0*(xl(j, 3) + xl(j, 4) + xl(j, 7) + xl(j, 8))
        end do
        ixl(25) = 25
      end if

      if(ixl(26) .eq. 0) then
        do j = 1, 3
          xl(j,26) = 0.50d0*(xl(j,16) + xl(j,12) + xl(j, 9) + xl(j,21))
     &             - 0.25d0*(xl(j, 4) + xl(j, 1) + xl(j, 8) + xl(j, 5))
        end do
        ixl(26) = 26
      end if

c       Center node

        if(ixl(27).eq.0) then
          do j = 1,3
c           xl(j,27) = 0.125d0*(xl(j, 1) +xl(j, 2) +xl(j, 3) +xl(j, 4)
c    &                        + xl(j, 5) +xl(j, 6) +xl(j, 7) +xl(j, 8))
c    &               - 0.250d0*(xl(j, 9) +xl(j,10) +xl(j,11) +xl(j,12)
c    &                        + xl(j,13) +xl(j,14) +xl(j,15) +xl(j,16)
c    &                        + xl(j,17) +xl(j,18) +xl(j,19) +xl(j,20))
c    &               + 0.500d0*(xl(j,21) +xl(j,22) +xl(j,23) +xl(j,24)
c    &                        + xl(j,25) +xl(j,26))
           xl(j,27) = 0.25d0*(xl(j,13) + xl(j,14) + xl(j,15) + xl(j,16)
     &                      + xl(j,18) + xl(j,19) + xl(j,20) + xl(j,21)
     &                      + xl(j,23) + xl(j,24) + xl(j,25) + xl(j,26)
     &                      - xl(j, 1) - xl(j, 2) - xl(j, 3) - xl(j, 4)
     &                      - xl(j, 5) - xl(j, 6) - xl(j, 7) - xl(j, 8))
          end do ! j
          ixl(27) = 27
        endif

c     endif

      end
c
c     ****************************************************************
c     *                                                              *
c     *               subroutine set_element_type                    *
c     *               (requires element type no.)                    *
c     *                                                              *
c     *                       written by: rhd                        *
c     *                   last modified : 02/1/02                    *
c     *                                                              *
c     ****************************************************************
c
      subroutine set_element_type( element_type, threed_solid_elem,
     &                             hex_elem, wedge_elem, tet_elem,
     &                             twod_elem, quad_elem,
     &                             triangle_elem, axisymm_elem,
     &                             cohesive_elem )
      implicit none
      integer element_type
      logical threed_solid_elem, hex_elem, wedge_elem, tet_elem,
     &        twod_elem, quad_elem, triangle_elem,
     &        axisymm_elem, cohesive_elem
c
      hex_elem      = element_type .ge. 1
     &                .and. element_type .le. 5
c
      wedge_elem    = element_type .eq. 7
c
      tet_elem      = element_type .eq. 6
     &                .or. element_type .eq. 13
c
      quad_elem     = element_type .eq. 9
     &                .or. element_type .eq. 10
c
      triangle_elem = element_type .eq. 8
     &                .or. element_type .eq. 11
c
      axisymm_elem  = element_type .eq. 10
     &                .or. element_type .eq. 11
c
      cohesive_elem = element_type .eq. 12
     &                .or. element_type .eq. 14
     &                .or. element_type .eq. 15
c
      threed_solid_elem = hex_elem .or. wedge_elem .or. tet_elem
c
      twod_elem = quad_elem .or. triangle_elem .or. axisymm_elem
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                    subroutine    vblke                       *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 5/2/18                      *
c     *                                                              *
c     ****************************************************************
c
      subroutine vblke(nr,ns,nt,ix,ni,ne,nen1,mat,ntyp,nf)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Set ninc = 2 for ntyp.le.15                      02/07/2007
c       2. Set ninc = 2 for ntyp.le.18                      04/09/2007
c       3. Replace itl by itq                               05/09/2007
c       4. Add 14/15 node generation                        06/09/2007
c       5. Correct allocation of active node for 14/15 node 20/10/2007
c       6. Add sets for 64-node brick type elements         06/02/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Generate a block of 3-d 8-node brick elements

c      Inputs:
c         nr        - Number elements in 1-local coordinate dir.
c         ns        - Number elements in 2-local coordinate dir.
c         nt        - Number elements in 3-local coordinate dir.
c         ni        - Initial node number for block
c         ne        - Initial element number for block
c         ndm       - Spatial dimension of coordinate array
c         nen1      - Dimension of ix array
c         mat       - Material set number for block
c         ntyp      - Element type for generations
c                     10:  8-node hexahedron  elements
c                     11:  4-node tetrahedron elements
c                     12: 27-node hexahedron  elements
c                     13: 10-node tetrahedron elements
c                     14: 20-node hexahedron  elements
c                     15: 11-node tetrahedron elements
c                     17: 14-node tetrahedron elements
c                     18: 15-node tetrahedron elements
c         prt       - Output generatied coordinates if true
c         prth      - Output title/header data if true

c      Outputs:
c         x(ndm,*)  - Nodal coordinates
c         ix(*)     - Element nodal connection list for block
c         nf        - Final   element number for block
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

c     include  'cdata.h'
c     include  'cdat2.h'
c     include  'compac.h'
c     include  'iofile.h'
c     include  'pconstant.h'
c     include  'trdata.h'

c     include  'pointer.h'
c     include  'comblk.h'

      logical   prt,prth, phd
      integer   ni,nf,ne,nn,ndm,nen1,mat,ma,ntyp, mct
      integer   nr,ns,nt,nrs,i,j,k,l,m,n,dlayer, ninc, nodi
      integer   ix(nen1,*),iq(27),it(10),nd(27),itq(10,6)
      integer :: temp
c     real*8    x(ndm,*)

      save

      data      itq / 1,2,4,5,  9, 21, 12, 17, 25, 23,
     &                2,3,4,8, 10, 11, 21, 27, 26, 20,
     &                2,4,5,8, 21, 23, 25, 27, 20, 16,
     &                2,6,3,8, 18, 24, 10, 27, 22, 26,
     &                3,6,7,8, 24, 14, 19, 26, 22, 15,
     &                5,6,2,8, 13, 18, 25, 16, 22, 27 /

c     Check generation order

      do i = 1,10
        it(i) = i
      end do ! i
      do i = 1,27
        iq(i) = i
      end do ! i
c     if    (ntyp.eq.12) then
c       nn = 27
c     elseif(ntyp.eq.14) then
c       nn = 20
c     elseif(ntyp.eq.19) then
c       nn = 64
c     else
c       nn = 8
c     endif
c     if(trdet.lt.0.0d0) then
c       do i = 1,4
c         iq(i+4) = i
c         iq(i  ) = i+4
c       end do ! i
c       if(ntyp.eq.12 .or. ntyp.eq.14) then
c         do i = 9,12
c           iq(i+4) = i
c           iq(i  ) = i+4
c         end do ! i
c         iq(21) = 22
c         iq(22) = 21
c       endif
c       i     = it(2)
c       it(2) = it(3)
c       it(3) = i
c       if(ntyp.eq.13 .or. ntyp.eq.15 .or.
c    &     ntyp.eq.17 .or. ntyp.eq.18) then
c         i      = it( 5)
c         it( 5) = it( 7)
c         it( 7) = i
c         i      = it( 9)
c         it( 9) = it(10)
c         it(10) = i
c       endif
c     endif

      nrs = nr*ns
      if(ntyp.le.11) then
        ninc  =  1
        nd(1) = -1
        nd(2) =  0
        nd(3) =      nr
        nd(4) = -1 + nr
        nd(5) = -1      + nrs
        nd(6) =           nrs
        nd(7) =      nr + nrs
        nd(8) = -1 + nr + nrs
      elseif(ntyp.le.13) then
        ninc  =  2
        nd( 1) = -1
        nd( 2) =  1
        nd( 3) =  1 + nr*2
        nd( 4) = -1 + nr*2
        nd( 5) = -1        + nrs*2
        nd( 6) =  1        + nrs*2
        nd( 7) =  1 + nr*2 + nrs*2
        nd( 8) = -1 + nr*2 + nrs*2
        nd( 9) =  0
        nd(10) =  1 + nr
        nd(11) =      nr*2
        nd(12) = -1 + nr
        nd(13) =             nrs*2
        nd(14) =  1 + nr   + nrs*2
        nd(15) =      nr*2 + nrs*2
        nd(16) = -1 + nr   + nrs*2
        nd(17) = -1        + nrs
        nd(18) =  1        + nrs
        nd(19) =  1 + nr*2 + nrs
        nd(20) = -1 + nr*2 + nrs
        nd(21) =      nr
        nd(22) =      nr   + nrs*2
        nd(23) = -1 + nr   + nrs
        nd(24) =  1 + nr   + nrs
        nd(25) =             nrs
        nd(26) =      nr*2 + nrs
        nd(27) =      nr   + nrs
      endif

c     Compute element connections

c     if(dlayer.ge.0) then
      ma = mat
c     endif
      nf = ne - 1
      temp = nr * ns * ninc
      do i = 1,nr-1,ninc
c       if(dlayer.eq.3) then
c         ma = ilr(k)
c       endif
        do j = 1,ns-1,ninc
c         if(dlayer.eq.2) then
c           ma = ilr(j)
c         endif
c         n = nr*(j-1 + ns*(k-1)) + ni
          n = nr * (j-1) + ni + i - temp
          do k = 1,nt-1,ninc
c           if(dlayer.eq.1) then
c             ma = ilr(i)
c           endif
c           n = n + 1
            n = n + temp

c           8-node hexahedral elements

            if(ntyp.eq.10) then
              nf = nf + 1
              ix(nen1,nf) = ma
              do m = 1,8
                ix(iq(m),nf) = n + nd(m)
              end do ! m

c           4-node tetrahedral elements

            elseif(ntyp.eq.11) then
              do l = 1,6
                nf = nf + 1
                ix(nen1,nf) = ma
                do m = 1,4
                  ix(it(m),nf) = n + nd(itq(m,l))
                end do ! m
              end do ! l

c           20 and 27-node hexahedral elements

            elseif(ntyp.eq.12) then
              nf = nf + 1
              ix(nen1,nf) = ma
              do m = 1,nn
                ix(iq(m),nf) = n + nd(m)
              end do ! m

c           10 and 14-node tetrahedral elements

            elseif(ntyp.eq.13) then
              do l = 1,6
                nf = nf + 1
                ix(nen1,nf) = ma
                do m = 1,10
                  ix(it(m),nf) = n + nd(itq(m,l))
                end do ! m
              end do ! l

            endif

            n = n + ninc - 1

          end do ! i
        end do ! j
      end do ! k

      end
c
c     ****************************************************************
c     *                                                              *
c     *                    subroutine    vblkn                       *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 5/2/18                      *
c     *                                                              *
c     ****************************************************************
c
      subroutine vblkn(nr,ns,nt,xl,x,ixl,dr,ds,dt,
     &                 ni,ndm,ctype,numnp,x0)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add 'cap' option for spherical cap               11/11/2008
c          and 'xcyl' and 'ycyl' for cylindrical sectors
c       2. Add option for generating sperical cap with      13/11/2008
c          specified thickness and constant radius.
c       3. Change print to 1p,5e13.4                        16/04/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Generate a block of 3-d 8-node brick elements

c      Inputs:
c         nr        - Number elements in 1-local coordinate dir.
c         ns        - Number elements in 2-local coordinate dir.
c         nt        - Number elements in 3-local coordinate dir.
c         xl(ndm,*) - Block nodal coordinate array
c         ixl(*)    - Block nodal connection list
c         dr        - 1-local coordinate increment
c         ds        - 2-local coordinate increment
c         dt        - 3-local coordinate increment
c         ni        - Initial node number for block
c         ndm       - Spatial dimension of mesh
c         ctype     - Type of block coordinates
c         prt       - Output generated data if true
c         prth      - Output title/header data if true

c      Outputs:
c         x(ndm,*)  - Nodal coordinates for block
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

c     include  'cdata.h'
c     include  'cdat2.h'
c     include  'iofile.h'
c     include  'trdata.h'

c                        global
      integer :: ndm, nr, ns, nt, ni, numnp
      integer :: x0(*), ixl(*)
      real(8) :: dr, ds, dt
      real(8) :: xl(3,*),x(ndm,*)
      character :: ctype*15

c                        local
      logical   prt,prth,phd, pcomp
      character xh*6
      integer   i,j,k,l,m,n,mct
      real*8    rr,sn2,cn2,sn3,cn3,afac
      real*8    ss(3),xx(3)

      save

      data      xh/' coord'/

c     Check that all corners of brick are defined

c     do k = 1,3
c       xx(k) = 0.0d0
c     end do ! k

      do k = 1,8
        if(ixl(k).ne.k) go to 900
      end do ! k
      call bcor3d(ixl,xl)
      n = ni
      mct = 0
      ss(3) = -1.0d0
      do k = 1,nt
        ss(2) = -1.0d0
        do j = 1,ns
          ss(1) = -1.0d0
          do i = 1,nr

c           Compute coordinates of node

            call xbcor3d(ss,xl, xx)

c           Transform to global coordinates

            do m = 1,ndm
c             x(m,n) = xr(m)+tr(m,1)*xx(1)+tr(m,2)*xx(2)+tr(m,3)*xx(3)
              x(m,n) = xx(m)
            end do ! m

            n = n + 1
            ss(1) = ss(1) + dr
          end do ! i
          ss(2) = ss(2) + ds
        end do ! j
        ss(3) = ss(3) + dt
      end do ! k

 900  return

      end
c
c     ****************************************************************
c     *                                                              *
c     *                    subroutine    xbcor3d                     *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified: 5/2/18                      *
c     *                                                              *
c     ****************************************************************
c
      subroutine xbcor3d(ss,xl,x)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Revise shape function numbering to standare      10/11/2008
c          27 node Lagrange element
c       2. Restore option for old numbering on flag oldfl   28/01/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute shape functions and coordinates for each point

c      Inputs:
c         ss(3)   - Natural coordinates for point
c         xl(3,*) - Nodal coordinates for brick

c      Outputs:
c         x(3)    - Cartesian coordinates for point
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

c     include  'corset.h'

      integer   j,l, ix(27),iy(27),iz(27)
      real*8    ss(3),xl(3,27),x(3), lshp(3,3),shp

c     save

      data      ix/1,3,3,1, 1,3,3,1, 1,3,3,1, 2,3,2,1,2,
     &             2,3,2,1,2, 2,3,2,1,2/
      data      iy/1,1,3,3, 1,1,3,3, 1,1,3,3, 1,2,3,2,2,
     &             1,2,3,2,2, 1,2,3,2,2/
      data      iz/1,1,1,1, 3,3,3,3, 2,2,2,2, 1,1,1,1,1,
     &             3,3,3,3,3, 2,2,2,2,2/

c     data      ox/1,3,3,1, 1,3,3,1, 1,3,3,1, 2,3,2,1,2, 2,3,2,1,2,
c    &             2,3,2,1,2/
c     data      oy/1,1,3,3, 1,1,3,3, 1,1,3,3, 1,2,3,2,2, 1,2,3,2,2,
c    &             1,2,3,2,2/
c     data      oz/1,1,1,1, 3,3,3,3, 2,2,2,2, 1,1,1,1,1, 3,3,3,3,3,
c    &             2,2,2,2,2/


      do j = 1,3
        lshp(1,j) = 0.5d0*ss(j)*(ss(j) - 1.d0)
        lshp(2,j) = (1.d0 - ss(j)*ss(j))
        lshp(3,j) = 0.5d0*ss(j)*(ss(j) + 1.d0)
      end do ! j

      do j = 1,3
        x(j) = 0.0d0
      end do ! j

      do l = 1,27
        shp = lshp(ix(l),1)*lshp(iy(l),2)*lshp(iz(l),3)
        do j = 1,3
          x(j) = x(j) + shp*xl(j,l)
        end do ! j
      end do ! l

      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine intpntb                      *
c     *                      subroutine shlb                         *
c     *                      subroutine shgb                         *
c     *                                                              *
c     *                   last modified: 5/2/18  RM                  *
c     *                                                              *
c     *      generate shape function, integration point, ...         *
c     *                                                              *
c     ****************************************************************
C-----------------------------------------------------------------------
      subroutine intpntb(l,lint,ib,w,ss)
C-----------------------------------------------------------------------
C      Purpose: Gauss quadrature for 3-d Brick element
c         Revised by Timothy Truster, 02/20/2012
C      Inputs:
C         ll     - Number of points/direction
C      Outputs:
C         lint   - Total number of quadrature points
C         s(4,*) - Gauss points (1-3) and weights (4)
C-----------------------------------------------------------------------
      implicit none

      integer   l1,l2,lint, i,j,k,l, ib
      real*8    s(4,125),gauss(10),weight(10)
      real*8    ss(3),w,i1,j1,one,two,three,g
      real*8    ls(9),lt(9)

      integer iin,iout,irsin,irsout
      common /iounit/ iin,iout,irsin,irsout

      data ls /-1.d0,1.d0,1.d0,-1.d0,0.d0,1.d0,0.d0,-1.d0,0.d0/
      data lt /-1.d0,-1.d0,1.d0,1.d0,-1.d0,0.d0,1.d0,0.d0,0.d0/

        one    = 1.0d0
        two    = 2.0d0
        three  = 3.0d0

      if(ib.eq.0) then

c     1 pt. quadrature
      if(lint.eq.1) then

        do i = 1,3
          ss(i) = 0.d0
        end do ! i
        w = 8.0d0

c     2 x 2 x 2 pt. quadrature
      elseif(lint.eq.8)         then
        l1 = 2
        gauss(1)  =-0.577350269189626D0
        gauss(2)  = 0.577350269189626D0
        weight(1) = 1.D0
        weight(2) = 1.D0
        i1 = dble(l)/dble(l1**2)
        i = ceiling(i1)
        l2 = (l-(i-1)*l1**2)
        j1 = dble(l2)/dble(l1)
        j = ceiling(j1)
        k = (l2-(j-1)*l1)
        ss(1)=gauss(k)
        ss(2)=gauss(j)
        ss(3)=gauss(i)
        w=weight(i)*weight(j)*weight(k)

c        for 3*3*3 pts quadrature
      elseif(lint.eq.27)         then
        l1        = 3
        gauss(1)  =-0.774596669241483d0
        gauss(2)  = 0.d0
        gauss(3)  = 0.774596669241483d0
        weight(1) = 0.555555555555556d0
        weight(2) = 0.888888888888889d0
        weight(3) = 0.555555555555556d0
        i1 = dble(l)/dble(l1**2)
        i = ceiling(i1)
        l2 = (l-(i-1)*l1**2)
        j1 = dble(l2)/dble(l1)
        j = ceiling(j1)
        k = (l2-(j-1)*l1)
        ss(1)=gauss(k)
        ss(2)=gauss(j)
        ss(3)=gauss(i)
        w=weight(i)*weight(j)*weight(k)

c        for 4*4*4 pts quadrature
      elseif(lint.eq.64)         then
        l1        = 4
        gauss(1)  =-0.861136311594953d0
        gauss(2)  =-0.339981043584856d0
        gauss(3)  = 0.339981043584856d0
        gauss(4)  = 0.861136311594953d0
        weight(1) = 0.347854845137454d0
        weight(2) = 0.652145154862546d0
        weight(3) = 0.652145154862546d0
        weight(4) = 0.347854845137454d0
        i1 = dble(l)/dble(l1**2)
        i = ceiling(i1)
        l2 = (l-(i-1)*l1**2)
        j1 = dble(l2)/dble(l1)
        j = ceiling(j1)
        k = (l2-(j-1)*l1)
        ss(1)=gauss(k)
        ss(2)=gauss(j)
        ss(3)=gauss(i)
        w=weight(i)*weight(j)*weight(k)

c      For 5*5*5 pts quadrature
       elseif(lint.eq.125)         then
        l1        = 5
        gauss(1)  =-0.906179845938664d0
        gauss(2)  =-0.538469310105683d0
        gauss(3)  = 0.d0
        gauss(4)  =+0.538469310105683d0
        gauss(5)  =+0.906179845938664d0
        weight(1) = 0.236926885056189d0
        weight(2) = 0.478628670449366d0
        weight(3) = 0.568888888888889d0
        weight(4) = 0.478628670449366d0
        weight(5) = 0.236926885056189d0
        i1 = dble(l)/dble(l1**2)
        i = ceiling(i1)
        l2 = (l-(i-1)*l1**2)
        j1 = dble(l2)/dble(l1)
        j = ceiling(j1)
        k = (l2-(j-1)*l1)
        ss(1)=gauss(k)
        ss(2)=gauss(j)
        ss(3)=gauss(i)
        w=weight(i)*weight(j)*weight(k)

c      For 10*10*10 pts quadrature
       elseif(lint.eq.1000)         then
        l1        = 1000

        gauss(1)  = -0.973906528517172
        gauss(2)  = -0.865063366688985
        gauss(3)  = -0.679409568299024
        gauss(4)  = -0.433395394129247
        gauss(5)  = -0.148874338981631
        gauss(6)  =  0.148874338981631
        gauss(7)  =  0.433395394129247
        gauss(8)  =  0.679409568299024
        gauss(9)  =  0.865063366688985
        gauss(10) =  0.973906528517172

        weight(1) = 0.0666713443086881
        weight(2) = 0.149451349150581
        weight(3) = 0.219086362515982
        weight(4) = 0.269266719309996
        weight(5) = 0.295524224714753
        weight(6) = 0.295524224714753
        weight(7) = 0.269266719309996
        weight(8) = 0.219086362515982
        weight(9) = 0.149451349150581
        weight(10)= 0.0666713443086881

        i1 = dble(l)/dble(l1**2)
        i = ceiling(i1)
        l2 = (l-(i-1)*l1**2)
        j1 = dble(l2)/dble(l1)
        j = ceiling(j1)
        k = (l2-(j-1)*l1)
        ss(1)=gauss(k)
        ss(2)=gauss(j)
        ss(3)=gauss(i)
        w=weight(i)*weight(j)*weight(k)

      endif

      else

c
c     2x2 integration.
c     ----------------
      if (lint .eq. 4) then
        l1 = 2
        gauss(1)  =-0.577350269189626D0
        gauss(2)  = 0.577350269189626D0
        weight(1) = 1.D0
        weight(2) = 1.D0
        j1 = dble(l)/dble(l1)
        j = ceiling(j1)
        k = (l-(j-1)*l1)
        ss(1)=gauss(k)
        ss(2)=gauss(j)
        ss(3)=-1.d0
        w=weight(j)*weight(k)
      endif

      endif

      return
      end

c----------------------------------------------------------------------
c
      subroutine shlb(ss,nel,nen,der,bf,shl,shld,shls,bubble)
c
c      3 dimensional 8 noded element + bubble point
c      By Kaiming Xia 09/20/2001
c         Second derivatives of shape functions for cartesian coordinates
c         By Kaiming Xia 11/03/2001
c      * * F E A P * * A Finite Element Analysis Program
c
c      3 dimensional 27 noded element
c      By Kaiming Xia 09/20/2001
c         Second derivatives of shape functions for cartesian coordinates
c         By Kaiming Xia 11/03/2001
c
c      Revised by Ramon Calderer 03/07/2008
c      Revised by Timothy Truster, 02/20/2012
c
c----------------------------------------------------------------------
c
c      Purpose: Compute 3-d isoparametric shape
c               functions and their derivatives w/r x,y,z

c      Inputs:
c         ss(3)     - Natural coordinates of point
c         xl(ndm,*) - Nodal coordinates for element
c         ndm       - Spatial dimension of mesh
c             ideriv          - Flag for activation of second derivatives (default=0
c                     don't compute second derivative
c      Outputs:
c         xsj       - Jacobian determinant at point
c         shp(4,*)  - Shape functions and derivatives at local point
c         shp2(6,*) - Shape functions and derivatives at local point
c                     cartd(1,i) = dN_i/dx
c                     cartd(2,i) = dN_i/dy
c                     cartd(3,i) = dN_i/dz
c                     cartd(4,i) =  N_i
c                     cartd2(1,i) = d^2N_i/dx^2
c                     cartd2(2,i) = d^2N_i/dy^2
c                     cartd2(3,i) = d^2N_i/dz^2
c                     cartd2(4,i) = d^2N_i/dxy
c                     cartd2(5,i) = d^2N_i/dyz
c                     cartd2(6,i) = d^2N_i/dzx
c
c----------------------------------------------------------------------

      implicit none

!     Input Variables
      integer nel,nen
      logical der,bf
      real*8  ss(3)

!     Output Variables
      real*8 shl(1,nen),shld(3,nen),shls(6,nen),bubble(4)

!     Local Variables
      integer   i , j , k
      real*8    onPLz8,onMIz8,onPLy8,onMIy8,x,y,z

      x=ss(1)
      y=ss(2)
      z=ss(3)
      onPLz8=(1.d0+z)/8.d0
      onMIz8=(1.d0-z)/8.d0
      onPLy8=(1.d0+y)/8.d0
      onMIy8=(1.d0-y)/8.d0

      if(nel.eq.8) then
c for 8 noded 3 dimensional brick element
        shld(1,1)=-(1.d0-y)*onMIz8
        shld(2,1)=-(1.d0-x)*onMIz8
        shld(3,1)=-(1.d0-x)*onMIy8
        shl(1,1)=(1.d0-x)*(1.d0-y)*onMIz8
        shld(1,2)=(1.d0-y)*onMIz8
        shld(2,2)=-(1.d0+x)*onMIz8
        shld(3,2)=-(1.d0+x)*onMIy8
        shl(1,2)=(1.d0+x)*(1.d0-y)*onMIz8
        shld(1,3)=(1.d0+y)*onMIz8
        shld(2,3)=(1.d0+x)*onMIz8
        shld(3,3)=-(1.d0+x)*onPLy8
        shl(1,3)=(1.d0+x)*(1.d0+y)*onMIz8
        shld(1,4)=-(1.d0+y)*onMIz8
        shld(2,4)=(1.d0-x)*onMIz8
        shld(3,4)=-(1.d0-x)*onPLy8
        shl(1,4)=(1.d0-x)*(1.d0+y)*onMIz8
        shld(1,5)=-(1.d0-y)*onPLz8
        shld(2,5)=-(1.d0-x)*onPLz8
        shld(3,5)=(1.d0-x)*onMIy8
        shl(1,5)=(1.d0-x)*(1.d0-y)*onPLz8
        shld(1,6)=(1.d0-y)*onPLz8
        shld(2,6)=-(1.d0+x)*onPLz8
        shld(3,6)=(1.d0+x)*onMIy8
        shl(1,6)=(1.d0+x)*(1.d0-y)*onPLz8
        shld(1,7)=(1.d0+y)*onPLz8
        shld(2,7)=(1.d0+x)*onPLz8
        shld(3,7)=(1.d0+x)*onPLy8
        shl(1,7)=(1.d0+x)*(1.d0+y)*onPLz8
        shld(1,8)=-(1.d0+y)*onPLz8
        shld(2,8)=(1.d0-x)*onPLz8
        shld(3,8)=(1.d0-x)*onPLy8
        shl(1,8)=(1.d0-x)*(1.d0+y)*onPLz8

c        Bubble point
        if(bf) then
          bubble(1)=(-2.d0*x)*(1.d0-y*y)*(1.d0-z*z)
          bubble(2)=(1.d0-x*x)*(-2.d0*y)*(1.d0-z*z)
          bubble(3)=(1.d0-x*x)*(1.d0-y*y)*(-2.d0*z)
          bubble(4)=(1.d0-x*x)*(1.d0-y*y)*(1.d0-z*z)
        endif

c    calculate second derivatives for shape funtions in local coordinate

        if(der) then
        shls(1,1)=0.d0
        shls(2,1)=0.d0
        shls(3,1)=0.d0
        shls(4,1)=onMIz8
        shls(5,1)=(1.d0-x)/8.d0
        shls(6,1)=onMIy8
        shls(1,2)=0.d0
        shls(2,2)=0.d0
        shls(3,2)=0.d0
        shls(4,2)=-onMIz8
        shls(5,2)=(1.d0+x)/8.d0
        shls(6,2)=-onMIy8
        shls(1,3)=0.d0
        shls(2,3)=0.d0
        shls(3,3)=0.d0
        shls(4,3)=onMIz8
        shls(5,3)=-(1.d0+x)/8.d0
        shls(6,3)=-onPLy8
        shls(1,4)=0.d0
        shls(2,4)=0.d0
        shls(3,4)=0.d0
        shls(4,4)=-onMIz8
        shls(5,4)=-(1.d0-x)/8.d0
        shls(6,4)=onPLy8
        shls(1,5)=0.d0
        shls(2,5)=0.d0
        shls(3,5)=0.d0
        shls(4,5)=onPLz8
        shls(5,5)=-(1.d0-x)/8.d0
        shls(6,5)=-onMIy8
        shls(1,6)=0.d0
        shls(2,6)=0.d0
        shls(3,6)=0.d0
        shls(4,6)=-onPLz8
        shls(5,6)=-(1.d0+x)/8.d0
        shls(6,6)=onMIy8
        shls(1,7)=0.d0
        shls(2,7)=0.d0
        shls(3,7)=0.d0
        shls(4,7)=onPLz8
        shls(5,7)=(1.d0+x)/8.d0
        shls(6,7)=onPLy8
        shls(1,8)=0.d0
        shls(2,8)=0.d0
        shls(3,8)=0.d0
        shls(4,8)=-onPLz8
        shls(5,8)=(1.d0-x)/8.d0
        shls(6,8)=-onPLy8
        end if

      endif
c
c for 27 noded 3 dimensional brick element
c
      if(nel.eq.27) then

        shld(1,1)=(2.d0*x-1.d0)*y*(y-1.d0)*z*(z-1.d0)/8.d0
        shld(2,1)=x*(x-1.d0)*(2.d0*y-1.d0)*z*(z-1.d0)/8.d0
        shld(3,1)=x*(x-1.d0)*y*(y-1.d0)*(2.d0*z-1.d0)/8.d0
        shl(1,1)=x*(x-1.d0)*y*(y-1.d0)*z*(z-1.d0)/8.d0

        shld(1,2)=(2.d0*x+1.d0)*y*(y-1.d0)*z*(z-1.d0)/8.d0
        shld(2,2)=x*(x+1.d0)*(2.d0*y-1.d0)*z*(z-1.d0)/8.d0
        shld(3,2)=x*(x+1.d0)*y*(y-1.d0)*(2.d0*z-1.d0)/8.d0
        shl(1,2)=x*(x+1.d0)*y*(y-1.d0)*z*(z-1.d0)/8.d0

        shld(1,3)=(2.d0*x+1.d0)*y*(y+1.d0)*z*(z-1.d0)/8.d0
        shld(2,3)=x*(x+1.d0)*(2.d0*y+1.d0)*z*(z-1.d0)/8.d0
        shld(3,3)=x*(x+1.d0)*y*(y+1.d0)*(2.d0*z-1.d0)/8.d0
        shl(1,3)=x*(x+1.d0)*y*(y+1.d0)*z*(z-1.d0)/8.d0

        shld(1,4)=(2.d0*x-1.d0)*y*(y+1.d0)*z*(z-1.d0)/8.d0
        shld(2,4)=x*(x-1.d0)*(2.d0*y+1.d0)*z*(z-1.d0)/8.d0
        shld(3,4)=x*(x-1.d0)*y*(y+1.d0)*(2.d0*z-1.d0)/8.d0
        shl(1,4)=x*(x-1.d0)*y*(y+1.d0)*z*(z-1.d0)/8.d0

        shld(1,5)=(2.d0*x-1.d0)*y*(y-1.d0)*z*(z+1.d0)/8.d0
        shld(2,5)=x*(x-1.d0)*(2.d0*y-1.d0)*z*(z+1.d0)/8.d0
        shld(3,5)=x*(x-1.d0)*y*(y-1.d0)*(2.d0*z+1.d0)/8.d0
        shl(1,5)=x*(x-1.d0)*y*(y-1.d0)*z*(z+1.d0)/8.d0

        shld(1,6)=(2.d0*x+1.d0)*y*(y-1.d0)*z*(z+1.d0)/8.d0
        shld(2,6)=x*(x+1.d0)*(2.d0*y-1.d0)*z*(z+1.d0)/8.d0
        shld(3,6)=x*(x+1.d0)*y*(y-1.d0)*(2.d0*z+1.d0)/8.d0
        shl(1,6)=x*(x+1.d0)*y*(y-1.d0)*z*(z+1.d0)/8.d0

        shld(1,7)=(2.d0*x+1.d0)*y*(y+1.d0)*z*(z+1.d0)/8.d0
        shld(2,7)=x*(x+1.d0)*(2.d0*y+1.d0)*z*(z+1.d0)/8.d0
        shld(3,7)=x*(x+1.d0)*y*(y+1.d0)*(2.d0*z+1.d0)/8.d0
        shl(1,7)=x*(x+1.d0)*y*(y+1.d0)*z*(z+1.d0)/8.d0

        shld(1,8)=(2.d0*x-1.d0)*y*(y+1.d0)*z*(z+1.d0)/8.d0
        shld(2,8)=x*(x-1.d0)*(2.d0*y+1.d0)*z*(z+1.d0)/8.d0
        shld(3,8)=x*(x-1.d0)*y*(y+1.d0)*(2.d0*z+1.d0)/8.d0
        shl(1,8)=x*(x-1.d0)*y*(y+1.d0)*z*(z+1.d0)/8.d0

        shld(1,17)=(2.d0*x-1.d0)*y*(y-1.d0)*(1.d0-z*z)/4.d0
        shld(2,17)=x*(x-1.d0)*(2.d0*y-1.d0)*(1.d0-z*z)/4.d0
        shld(3,17)=x*(x-1.d0)*y*(y-1.d0)*(-2.d0*z)/4.d0
        shl(1,17)=x*(x-1.d0)*y*(y-1.d0)*(1.d0-z*z)/4.d0

        shld(1,18)=(2.d0*x+1.d0)*y*(y-1.d0)*(1.d0-z*z)/4.d0
        shld(2,18)=x*(x+1.d0)*(2.d0*y-1.d0)*(1.d0-z*z)/4.d0
        shld(3,18)=x*(x+1.d0)*y*(y-1.d0)*(-2.d0*z)/4.d0
        shl(1,18)=x*(x+1.d0)*y*(y-1.d0)*(1.d0-z*z)/4.d0

        shld(1,19)=(2.d0*x+1.d0)*y*(y+1.d0)*(1.d0-z*z)/4.d0
        shld(2,19)=x*(x+1.d0)*(2.d0*y+1.d0)*(1.d0-z*z)/4.d0
        shld(3,19)=x*(x+1.d0)*y*(y+1.d0)*(-2.d0*z)/4.d0
        shl(1,19)=x*(x+1.d0)*y*(y+1.d0)*(1.d0-z*z)/4.d0

        shld(1,20)=(2.d0*x-1.d0)*y*(y+1.d0)*(1.d0-z*z)/4.d0
        shld(2,20)=x*(x-1.d0)*(2.d0*y+1.d0)*(1.d0-z*z)/4.d0
        shld(3,20)=x*(x-1.d0)*y*(y+1.d0)*(-2.d0*z)/4.d0
        shl(1,20)=x*(x-1.d0)*y*(y+1.d0)*(1.d0-z*z)/4.d0

        shld(1,9)=(-2.d0*x)*y*(y-1.d0)*z*(z-1.d0)/4.d0
        shld(2,9)=(1.d0-x*x)*(2.d0*y-1.d0)*z*(z-1.d0)/4.d0
        shld(3,9)=(1.d0-x*x)*y*(y-1.d0)*(2.d0*z-1.d0)/4.d0
        shl(1,9)=(1.d0-x*x)*y*(y-1.d0)*z*(z-1.d0)/4.d0

        shld(1,10)=(2.d0*x+1.d0)*(1.d0-y*y)*z*(z-1.d0)/4.d0
        shld(2,10)=x*(x+1.d0)*(-2.d0*y)*z*(z-1.d0)/4.d0
        shld(3,10)=x*(x+1.d0)*(1.d0-y*y)*(2.d0*z-1.d0)/4.d0
        shl(1,10)=x*(x+1.d0)*(1.d0-y*y)*z*(z-1.d0)/4.d0

        shld(1,11)=(-2.d0*x)*y*(y+1.d0)*z*(z-1.d0)/4.d0
        shld(2,11)=(1.d0-x*x)*(2.d0*y+1.d0)*z*(z-1.d0)/4.d0
        shld(3,11)=(1.d0-x*x)*y*(y+1.d0)*(2.d0*z-1.d0)/4.d0
        shl(1,11)=(1.d0-x*x)*y*(y+1.d0)*z*(z-1.d0)/4.d0

        shld(1,12)=(2.d0*x-1.d0)*(1.d0-y*y)*z*(z-1.d0)/4.d0
        shld(2,12)=x*(x-1.d0)*(-2.d0*y)*z*(z-1.d0)/4.d0
        shld(3,12)=x*(x-1.d0)*(1.d0-y*y)*(2.d0*z-1.d0)/4.d0
        shl(1,12)=x*(x-1.d0)*(1.d0-y*y)*z*(z-1.d0)/4.d0

        shld(1,21)=(-2.d0*x)*(1.d0-y*y)*z*(z-1.d0)/2.d0
        shld(2,21)=(1.d0-x*x)*(-2.d0*y)*z*(z-1.d0)/2.d0
        shld(3,21)=(1.d0-x*x)*(1.d0-y*y)*(2.d0*z-1.d0)/2.d0
        shl(1,21)=(1.d0-x*x)*(1.d0-y*y)*z*(z-1.d0)/2.d0

        shld(1,13)=(-2.d0*x)*y*(y-1.d0)*z*(z+1.d0)/4.d0
        shld(2,13)=(1.d0-x*x)*(2.d0*y-1.d0)*z*(z+1.d0)/4.d0
        shld(3,13)=(1.d0-x*x)*y*(y-1.d0)*(2.d0*z+1.d0)/4.d0
        shl(1,13)=(1.d0-x*x)*y*(y-1.d0)*z*(z+1.d0)/4.d0

        shld(1,14)=(2.d0*x+1.d0)*(1.d0-y*y)*z*(z+1.d0)/4.d0
        shld(2,14)=x*(x+1.d0)*(-2.d0*y)*z*(z+1.d0)/4.d0
        shld(3,14)=x*(x+1.d0)*(1.d0-y*y)*(2.d0*z+1.d0)/4.d0
        shl(1,14)=x*(x+1.d0)*(1.d0-y*y)*z*(z+1.d0)/4.d0

        shld(1,15)=(-2.d0*x)*y*(y+1.d0)*z*(z+1.d0)/4.d0
        shld(2,15)=(1.d0-x*x)*(2.d0*y+1.d0)*z*(z+1.d0)/4.d0
        shld(3,15)=(1.d0-x*x)*y*(y+1.d0)*(2.d0*z+1.d0)/4.d0
        shl(1,15)=(1.d0-x*x)*y*(y+1.d0)*z*(z+1.d0)/4.d0

        shld(1,16)=(2.d0*x-1.d0)*(1.d0-y*y)*z*(z+1.d0)/4.d0
        shld(2,16)=x*(x-1.d0)*(-2.0*y)*z*(z+1.d0)/4.d0
        shld(3,16)=x*(x-1.d0)*(1.d0-y*y)*(2.d0*z+1.d0)/4.d0
        shl(1,16)=x*(x-1.d0)*(1.d0-y*y)*z*(z+1.d0)/4.d0

        shld(1,22)=(-2.d0*x)*(1.d0-y*y)*z*(z+1.d0)/2.d0
        shld(2,22)=(1.d0-x*x)*(-2.d0*y)*z*(z+1.d0)/2.d0
        shld(3,22)=(1.d0-x*x)*(1.d0-y*y)*(2.d0*z+1.d0)/2.d0
        shl(1,22)=(1.d0-x*x)*(1.d0-y*y)*z*(z+1.d0)/2.d0

        shld(1,25)=(-2.d0*x)*y*(y-1.d0)*(1.d0-z*z)/2.d0
        shld(2,25)=(1.d0-x*x)*(2.d0*y-1.d0)*(1.d0-z*z)/2.d0
        shld(3,25)=(1.d0-x*x)*y*(y-1.d0)*(-2.d0*z)/2.d0
        shl(1,25)=(1.d0-x*x)*y*(y-1.d0)*(1.d0-z*z)/2.d0

        shld(1,24)=(2.d0*x+1.d0)*(1.d0-y*y)*(1.d0-z*z)/2.d0
        shld(2,24)=x*(x+1.d0)*(-2.d0*y)*(1.d0-z*z)/2.d0
        shld(3,24)=x*(x+1.d0)*(1.d0-y*y)*(-2.d0*z)/2.d0
        shl(1,24)=x*(x+1.d0)*(1.d0-y*y)*(1.d0-z*z)/2.d0

        shld(1,26)=(-2.d0*x)*y*(y+1.d0)*(1.d0-z*z)/2.d0
        shld(2,26)=(1.d0-x*x)*(2.d0*y+1.d0)*(1.d0-z*z)/2.d0
        shld(3,26)=(1.d0-x*x)*y*(y+1.d0)*(-2.d0*z)/2.d0
        shl(1,26)=(1.d0-x*x)*y*(y+1.d0)*(1.d0-z*z)/2.d0

        shld(1,23)=(2.d0*x-1.d0)*(1.d0-y*y)*(1.d0-z*z)/2.d0
        shld(2,23)=x*(x-1.d0)*(-2.d0*y)*(1.d0-z*z)/2.d0
        shld(3,23)=x*(x-1.d0)*(1.d0-y*y)*(-2.d0*z)/2.d0
        shl(1,23)=x*(x-1.d0)*(1.d0-y*y)*(1.d0-z*z)/2.d0

        shld(1,27)=(-2.d0*x)*(1.d0-y*y)*(1.d0-z*z)
        shld(2,27)=(1.d0-x*x)*(-2.d0*y)*(1.d0-z*z)
        shld(3,27)=(1.d0-x*x)*(1.d0-y*y)*(-2.d0*z)
        shl(1,27)=(1.d0-x*x)*(1.d0-y*y)*(1.d0-z*z)

c       Bubble node
        if(bf) then
          bubble(1)=(-4.d0*x*x*x)*(1.d0-y*y*y*y)*(1.d0-z*z*z*z)
          bubble(2)=(1.d0-x*x*x*x)*(-4.d0*y*y*y)*(1.d0-z*z*z*z)
          bubble(3)=(1.d0-x*x*x*x)*(1.d0-y*y*y*y)*(-4.d0*z*z*z)
          bubble(4)=(1.d0-x*x*x*x)*(1.d0-y*y*y*y)*(1.d0-z*z*z*z)
        endif

c     calculate second derivatives for shape funtions in local coordinate

      if(der) then

        shls(1,1)=2.d0*y*(y-1.d0)*z*(z-1.d0)/8.d0
        shls(2,1)=x*(x-1.d0)*2.d0*z*(z-1.d0)/8.d0
        shls(3,1)=x*(x-1.d0)*y*(y-1.d0)*2.d0/8.d0
        shls(4,1)=(2.d0*x-1.d0)*(2.d0*y-1.d0)*z*(z-1.d0)/8.d0
        shls(5,1)=x*(x-1.d0)*(2.d0*y-1.d0)*(2.d0*z-1.d0)/8.d0
        shls(6,1)=(2.d0*x-1.d0)*y*(y-1.d0)*(2.d0*z-1.d0)/8.d0

        shls(1,2)=2.d0*y*(y-1.d0)*z*(z-1.d0)/8.d0
        shls(2,2)=x*(x+1.d0)*2.d0*z*(z-1.d0)/8.d0
        shls(3,2)=x*(x+1.d0)*y*(y-1.d0)*2.d0/8.d0
        shls(4,2)=(2.d0*x+1.d0)*(2.d0*y-1.d0)*z*(z-1.d0)/8.d0
        shls(5,2)=x*(x+1.d0)*(2.d0*y-1.d0)*(2.d0*z-1.d0)/8.d0
        shls(6,2)=(2.d0*x+1.d0)*y*(y-1.d0)*(2.d0*z-1.d0)/8.d0

        shls(1,3)=2.d0*y*(y+1.d0)*z*(z-1.d0)/8.d0
        shls(2,3)=x*(x+1.d0)*2.d0*z*(z-1.d0)/8.d0
        shls(3,3)=x*(x+1.d0)*y*(y+1.d0)*2.d0/8.d0
        shls(4,3)=(2.d0*x+1.d0)*(2.d0*y+1.d0)*z*(z-1.d0)/8.d0
        shls(5,3)=x*(x+1.d0)*(2.d0*y+1.d0)*(2.d0*z-1.d0)/8.d0
        shls(6,3)=(2.d0*x+1.d0)*y*(y+1.d0)*(2.d0*z-1.d0)/8.d0

        shls(1,4)=2.d0*y*(y+1.d0)*z*(z-1.d0)/8.d0
        shls(2,4)=x*(x-1.d0)*2.d0*z*(z-1.d0)/8.d0
        shls(3,4)=x*(x-1.d0)*y*(y+1.d0)*2.d0/8.d0
        shls(4,4)=(2.d0*x-1.d0)*(2.d0*y+1.d0)*z*(z-1.d0)/8.d0
        shls(5,4)=x*(x-1.d0)*(2.d0*y+1.d0)*(2.d0*z-1.d0)/8.d0
        shls(6,4)=(2.d0*x-1.d0)*y*(y+1.d0)*(2.d0*z-1.d0)/8.d0

        shls(1,5)=2.d0*y*(y-1.d0)*z*(z+1.d0)/8.d0
        shls(2,5)=x*(x-1.d0)*2.d0*z*(z+1.d0)/8.d0
        shls(3,5)=x*(x-1.d0)*y*(y-1.d0)*2.d0/8.d0
        shls(4,5)=(2.d0*x-1.d0)*(2.d0*y-1.d0)*z*(z+1.d0)/8.d0
        shls(5,5)=x*(x-1.d0)*(2.d0*y-1.d0)*(2.d0*z+1.d0)/8.d0
        shls(6,5)=(2.d0*x-1.d0)*y*(y-1.d0)*(2.d0*z+1.d0)/8.d0

        shls(1,6)=2.d0*y*(y-1.d0)*z*(z+1.d0)/8.d0
        shls(2,6)=x*(x+1.d0)*2.d0*z*(z+1.d0)/8.d0
        shls(3,6)=x*(x+1.d0)*y*(y-1.d0)*2.d0/8.d0
        shls(4,6)=(2.d0*x+1.d0)*(2.d0*y-1.d0)*z*(z+1.d0)/8.d0
        shls(5,6)=x*(x+1.d0)*(2.d0*y-1.d0)*(2.d0*z+1.d0)/8.d0
        shls(6,6)=(2.d0*x+1.d0)*y*(y-1.d0)*(2.d0*z+1.d0)/8.d0

        shls(1,7)=2.d0*y*(y+1.d0)*z*(z+1.d0)/8.d0
        shls(2,7)=x*(x+1.d0)*2.d0*z*(z+1.d0)/8.d0
        shls(3,7)=x*(x+1.d0)*y*(y+1.d0)*2.d0/8.d0
        shls(4,7)=(2.d0*x+1.d0)*(2.d0*y+1.d0)*z*(z+1.d0)/8.d0
        shls(5,7)=x*(x+1.d0)*(2.d0*y+1.d0)*(2.d0*z+1.d0)/8.d0
        shls(6,7)=(2.d0*x+1.d0)*y*(y+1.d0)*(2.d0*z+1.d0)/8.d0

        shls(1,8)=2.d0*y*(y+1.d0)*z*(z+1.d0)/8.d0
        shls(2,8)=x*(x-1.d0)*2.d0*z*(z+1.d0)/8.d0
        shls(3,8)=x*(x-1.d0)*y*(y+1.d0)*2.d0/8.d0
        shls(4,8)=(2.d0*x-1.d0)*(2.d0*y+1.d0)*z*(z+1.d0)/8.d0
        shls(5,8)=x*(x-1.d0)*(2.d0*y+1.d0)*(2.d0*z+1.d0)/8.d0
        shls(6,8)=(2.d0*x-1.d0)*y*(y+1.d0)*(2.d0*z+1.d0)/8.d0

        shls(1,17)=2.d0*y*(y-1.d0)*(1.d0-z*z)/4.d0
        shls(2,17)=x*(x-1.d0)*2.d0*(1.d0-z*z)/4.d0
        shls(3,17)=x*(x-1.d0)*y*(y-1.d0)*(-2.d0)/4.d0
        shls(4,17)=(2.d0*x-1.d0)*(2.d0*y-1.d0)*(1.d0-z*z)/4.d0
        shls(5,17)=x*(x-1.d0)*(2.d0*y-1.d0)*(-2.d0*z)/4.d0
        shls(6,17)=(2.d0*x-1.d0)*y*(y-1.d0)*(-2.d0*z)/4.d0

        shls(1,18)=2.d0*y*(y-1.d0)*(1.d0-z*z)/4.d0
        shls(2,18)=x*(x+1.d0)*2.d0*(1.d0-z*z)/4.d0
        shls(3,18)=x*(x+1.d0)*y*(y-1.d0)*(-2.d0)/4.d0
        shls(4,18)=(2.d0*x+1.d0)*(2.d0*y-1.d0)*(1.d0-z*z)/4.d0
        shls(5,18)=x*(x+1.d0)*(2.d0*y-1.d0)*(-2.d0*z)/4.d0
        shls(6,18)=(2.d0*x+1.d0)*y*(y-1.d0)*(-2.d0*z)/4.d0

        shls(1,19)=2.d0*y*(y+1.d0)*(1.d0-z*z)/4.d0
        shls(2,19)=x*(x+1.d0)*2.d0*(1.d0-z*z)/4.d0
        shls(3,19)=x*(x+1.d0)*y*(y+1.d0)*(-2.d0)/4.d0
        shls(4,19)=(2.d0*x+1.d0)*(2.d0*y+1.d0)*(1.d0-z*z)/4.d0
        shls(5,19)=x*(x+1.d0)*(2.d0*y+1.d0)*(-2.d0*z)/4.d0
        shls(6,19)=(2.d0*x+1.d0)*y*(y+1.d0)*(-2.d0*z)/4.d0

        shls(1,20)=2.d0*y*(y+1.d0)*(1.d0-z*z)/4.d0
        shls(2,20)=x*(x-1.d0)*2.d0*(1.d0-z*z)/4.d0
        shls(3,20)=x*(x-1.d0)*y*(y+1.d0)*(-2.d0)/4.d0
        shls(4,20)=(2.d0*x-1.d0)*(2.d0*y+1.d0)*(1.d0-z*z)/4.d0
        shls(5,20)=x*(x-1.d0)*(2.d0*y+1.d0)*(-2.d0*z)/4.d0
        shls(6,20)=(2.d0*x-1.d0)*y*(y+1.d0)*(-2.d0*z)/4.d0

        shls(1,9)=(-2.d0)*y*(y-1.d0)*z*(z-1.d0)/4.d0
        shls(2,9)=(1.d0-x*x)*2.d0*z*(z-1.d0)/4.d0
        shls(3,9)=(1.d0-x*x)*y*(y-1.d0)*2.d0/4.d0
        shls(4,9)=(-2.d0*x)*(2.d0*y-1.d0)*z*(z-1.d0)/4.d0
        shls(5,9)=(1.d0-x*x)*(2.d0*y-1.d0)*(2.d0*z-1.d0)/4.d0
        shls(6,9)=(-2.d0*x)*y*(y-1.d0)*(2.d0*z-1.d0)/4.d0

        shls(1,10)=2.d0*(1.d0-y*y)*z*(z-1.d0)/4.d0
        shls(2,10)=x*(x+1.d0)*(-2.d0)*z*(z-1.d0)/4.d0
        shls(3,10)=x*(x+1.d0)*(1.d0-y*y)*2.d0/4.d0
        shls(4,10)=(2.d0*x+1.d0)*(-2.d0*y)*z*(z-1.d0)/4.d0
        shls(5,10)=x*(x+1.d0)*(-2.d0*y)*(2.d0*z-1.d0)/4.d0
        shls(6,10)=(2.D0*x+1.d0)*(1.d0-y*y)*(2.d0*z-1.d0)/4.d0

        shls(1,11)=(-2.d0)*y*(y+1.d0)*z*(z-1.d0)/4.d0
        shls(2,11)=(1.d0-x*x)*2.d0*z*(z-1.d0)/4.d0
        shls(3,11)=(1.d0-x*x)*y*(y+1.d0)*2.d0/4.d0
        shls(4,11)=(-2.d0*x)*(2.D0*y+1.d0)*z*(z-1.d0)/4.d0
        shls(5,11)=(1.d0-x*x)*(2.d0*y+1.d0)*(2.D0*z-1.d0)/4.d0
        shls(6,11)=(-2.D0*x)*y*(y+1.d0)*(2.d0*z-1.d0)/4.d0

        shls(1,12)=2.d0*(1.d0-y*y)*z*(z-1.d0)/4.d0
        shls(2,12)=x*(x-1.d0)*(-2.d0)*z*(z-1.d0)/4.d0
        shls(3,12)=x*(x-1.d0)*(1.d0-y*y)*2.d0/4.d0
        shls(4,12)=(2.d0*x-1.d0)*(-2.D0*y)*z*(z-1.d0)/4.d0
        shls(5,12)=x*(x-1.d0)*(-2.d0*y)*(2.D0*z-1.d0)/4.d0
        shls(6,12)=(2.D0*x-1.d0)*(1.d0-y*y)*(2.d0*z-1.d0)/4.d0

        shls(1,21)=(-2.d0)*(1.d0-y*y)*z*(z-1.d0)/2.d0
        shls(2,21)=(1.d0-x*x)*(-2.d0)*z*(z-1.d0)/2.d0
        shls(3,21)=(1.d0-x*x)*(1.d0-y*y)*2.d0/2.d0
        shls(4,21)=(-2.d0*x)*(-2.D0*y)*z*(z-1.d0)/2.d0
        shls(5,21)=(1.d0-x*x)*(-2.d0*y)*(2.D0*z-1.d0)/2.d0
        shls(6,21)=(-2.D0*x)*(1.d0-y*y)*(2.d0*z-1.d0)/2.d0

        shls(1,13)=(-2.d0)*y*(y-1.d0)*z*(z+1.d0)/4.d0
        shls(2,13)=(1.d0-x*x)*2.d0*z*(z+1.d0)/4.d0
        shls(3,13)=(1.d0-x*x)*y*(y-1.d0)*2.d0/4.d0
        shls(4,13)=(-2.d0*x)*(2.D0*y-1.d0)*z*(z+1.d0)/4.d0
        shls(5,13)=(1.d0-x*x)*(2.d0*y-1.d0)*(2.D0*z+1.d0)/4.d0
        shls(6,13)=(-2.D0*x)*y*(y-1.d0)*(2.d0*z+1.d0)/4.d0

        shls(1,14)=2.d0*(1.d0-y*y)*z*(z+1.d0)/4.d0
        shls(2,14)=x*(x+1.d0)*(-2.d0)*z*(z+1.d0)/4.d0
        shls(3,14)=x*(x+1.d0)*(1.d0-y*y)*2.d0/4.d0
        shls(4,14)=(2.d0*x+1.d0)*(-2.D0*y)*z*(z+1.d0)/4.d0
        shls(5,14)=x*(x+1.d0)*(-2.d0*y)*(2.D0*z+1.d0)/4.d0
        shls(6,14)=(2.D0*x+1.d0)*(1.d0-y*y)*(2.d0*z+1.d0)/4.d0

        shls(1,15)=(-2.d0)*y*(y+1.d0)*z*(z+1.d0)/4.d0
        shls(2,15)=(1.d0-x*x)*2.d0*z*(z+1.d0)/4.d0
        shls(3,15)=(1.d0-x*x)*y*(y+1.d0)*2.d0/4.d0
        shls(4,15)=(-2.d0*x)*(2.D0*y+1.d0)*z*(z+1.d0)/4.d0
        shls(5,15)=(1.d0-x*x)*(2.d0*y+1.d0)*(2.D0*z+1.d0)/4.d0
        shls(6,15)=(-2.D0*x)*y*(y+1.d0)*(2.d0*z+1.d0)/4.d0

        shls(1,16)=2.d0*(1.d0-y*y)*z*(z+1.d0)/4.d0
        shls(2,16)=x*(x-1.d0)*(-2.0)*z*(z+1.d0)/4.d0
        shls(3,16)=x*(x-1.d0)*(1.d0-y*y)*2.d0/4.d0
        shls(4,16)=(2.d0*x-1.d0)*(-2.d0*y)*z*(z+1.d0)/4.d0
        shls(5,16)=x*(x-1.d0)*(-2.0*y)*(2.d0*z+1.d0)/4.d0
        shls(6,16)=(2.d0*x-1.d0)*(1.d0-y*y)*(2.d0*z+1.d0)/4.d0

        shls(1,22)=(-2.d0)*(1.d0-y*y)*z*(z+1.d0)/2.d0
        shls(2,22)=(1.d0-x*x)*(-2.d0)*z*(z+1.d0)/2.d0
        shls(3,22)=(1.d0-x*x)*(1.d0-y*y)*2.d0/2.d0
        shls(4,22)=(-2.d0*x)*(-2.d0*y)*z*(z+1.d0)/2.d0
        shls(5,22)=(1.d0-x*x)*(-2.d0*y)*(2.d0*z+1.d0)/2.d0
        shls(6,22)=(-2.d0*x)*(1.d0-y*y)*(2.d0*z+1.d0)/2.d0

        shls(1,25)=(-2.d0)*y*(y-1.d0)*(1.d0-z*z)/2.d0
        shls(2,25)=(1.d0-x*x)*2.d0*(1.d0-z*z)/2.d0
        shls(3,25)=(1.d0-x*x)*y*(y-1.d0)*(-2.d0)/2.d0
        shls(4,25)=(-2.d0*x)*(2.d0*y-1.d0)*(1.d0-z*z)/2.d0
        shls(5,25)=(1.d0-x*x)*(2.d0*y-1.d0)*(-2.d0*z)/2.d0
        shls(6,25)=(-2.d0*x)*y*(y-1.d0)*(-2.d0*z)/2.d0

        shls(1,24)=2.d0*(1.d0-y*y)*(1.d0-z*z)/2.d0
        shls(2,24)=x*(x+1.d0)*(-2.d0)*(1.d0-z*z)/2.d0
        shls(3,24)=x*(x+1.d0)*(1.d0-y*y)*(-2.d0)/2.d0
        shls(4,24)=(2.d0*x+1.d0)*(-2.d0*y)*(1.d0-z*z)/2.d0
        shls(5,24)=x*(x+1.d0)*(-2.d0*y)*(-2.d0*z)/2.d0
        shls(6,24)=(2.d0*x+1.d0)*(1.d0-y*y)*(-2.d0*z)/2.d0

        shls(1,26)=(-2.d0)*y*(y+1.d0)*(1.d0-z*z)/2.d0
        shls(2,26)=(1.d0-x*x)*2.d0*(1.d0-z*z)/2.d0
        shls(3,26)=(1.d0-x*x)*y*(y+1.d0)*(-2.d0)/2.d0
        shls(4,26)=(-2.d0*x)*(2.d0*y+1.d0)*(1.d0-z*z)/2.d0
        shls(5,26)=(1.d0-x*x)*(2.d0*y+1.d0)*(-2.d0*z)/2.d0
        shls(6,26)=(-2.d0*x)*y*(y+1.d0)*(-2.d0*z)/2.d0

        shls(1,23)=2.d0*(1.d0-y*y)*(1.d0-z*z)/2.d0
        shls(2,23)=x*(x-1.d0)*(-2.d0)*(1.d0-z*z)/2.d0
        shls(3,23)=x*(x-1.d0)*(1.d0-y*y)*(-2.d0)/2.d0
        shls(4,23)=(2.d0*x-1.d0)*(-2.d0*y)*(1.d0-z*z)/2.d0
        shls(5,23)=x*(x-1.d0)*(-2.d0*y)*(-2.d0*z)/2.d0
        shls(6,23)=(2.d0*x-1.d0)*(1.d0-y*y)*(-2.d0*z)/2.d0

        shls(1,27)=(-2.d0)*(1.d0-y*y)*(1.d0-z*z)
        shls(2,27)=(1.d0-x*x)*(-2.d0)*(1.d0-z*z)
        shls(3,27)=(1.d0-x*x)*(1.d0-y*y)*(-2.d0)
        shls(4,27)=(-2.d0*x)*(-2.d0*y)*(1.d0-z*z)
        shls(5,27)=(1.d0-x*x)*(-2.d0*y)*(-2.d0*z)
        shls(6,27)=(-2.d0*x)*(1.d0-y*y)*(-2.d0*z)

      endif

      endif

      return
      end

c----------------------------------------------------------------------
c
      subroutine shgb(xl,nel,shp,shp2,neg,nen,bf,der,xsj,cartd1,
     &                cartd2,bubble,sx)
c
c      8 noded, 27 noded brick element + bubble point
c
c      Written by Kaiming Xia, 09/20/2001
c         Revised by Ramon Calderer, 03/07/2008
c         Revised by Timothy Truster, 02/20/2012
c
c----------------------------------------------------------------------
c
c      Purpose: Compute 3-d isoparametric shape
c               functions and their derivatives w/r x,y,z

c      Inputs:
c         ss(3)     - Natural coordinates of point
c         xl(ndm,*) - Nodal coordinates for element
c         ndm       - Spatial dimension of mesh
c             ideriv          - Flag for activation of second derivatives (default=0
c                     don't compute second derivative
c      Outputs:
c         xsj       - Jacobian determinant at point
c         shp(4,*)  - Shape functions and derivatives at local point
c         shp2(6,*) - Shape functions and derivatives at local point
c                     cartd(1,i) = dN_i/dx
c                     cartd(2,i) = dN_i/dy
c                     cartd(3,i) = dN_i/dz
c                     cartd(4,i) =  N_i
c                     cartd2(1,i) = d^2N_i/dx^2
c                     cartd2(2,i) = d^2N_i/dy^2
c                     cartd2(3,i) = d^2N_i/dz^2
c                     cartd2(4,i) = d^2N_i/dxy
c                     cartd2(5,i) = d^2N_i/dyz
c                     cartd2(6,i) = d^2N_i/dzx
c
c----------------------------------------------------------------------

      implicit none

!     Input Variables
      integer nel,nen,neg
      logical der,bf
      real*8 xl(3,nel),shp(3,nen),shp2(6,nen),bubble(4)

!     Output Variables
      real*8 xsj, cartd1(3,nen),cartd2(6,nen),sx(3,3)

!     Local Variables
      integer   i , j , k, inode, n
      real*8    rxsj,c1(6,3),tc(6,3),tcj(6,3),at1(6),at2(6)
      real*8    t2(6,6),temp(3)
      real*8    xs(3,3),ad(3,3),u

      integer iin,iout,irsin,irsout
      common /iounit/ iin,iout,irsin,irsout

c     Compute jacobian transformation
c        Calculate c1 matrix simultaneously

        do 30 i=1,3
        do 30 j=1,3
        xs(i,j)=0.d0
        c1(i,j)=0.d0
        do 30 inode=1,nel
        xs(i,j)=xs(i,j)+shp(i,inode)*xl(j,inode)
        c1(i,j)=c1(i,j)+shp2(i,inode)*xl(j,inode)
30        continue
        do 40 i=4,6
        do 40 j=1,3
        c1(i,j)=0.d0
        do 40 inode=1,nel
        c1(i,j)=c1(i,j)+shp2(i,inode)*xl(j,inode)
40        continue

! Output transpose which agrees with customary array definition
      do i = 1,3
        do j = 1,3
          sx(i,j) = xs(j,i)
        enddo
      enddo

c below is to calculate determinant and inverswe of jacobian matrix
c
c     Compute adjoint to jacobian

      ad(1,1) = xs(2,2)*xs(3,3) - xs(2,3)*xs(3,2)
      ad(1,2) = xs(3,2)*xs(1,3) - xs(3,3)*xs(1,2)
      ad(1,3) = xs(1,2)*xs(2,3) - xs(1,3)*xs(2,2)
      ad(2,1) = xs(2,3)*xs(3,1) - xs(2,1)*xs(3,3)
      ad(2,2) = xs(3,3)*xs(1,1) - xs(3,1)*xs(1,3)
      ad(2,3) = xs(1,3)*xs(2,1) - xs(1,1)*xs(2,3)
      ad(3,1) = xs(2,1)*xs(3,2) - xs(2,2)*xs(3,1)
      ad(3,2) = xs(3,1)*xs(1,2) - xs(3,2)*xs(1,1)
      ad(3,3) = xs(1,1)*xs(2,2) - xs(1,2)*xs(2,1)

c     Compute determinant of jacobian

      xsj  = xs(1,1)*ad(1,1) + xs(1,2)*ad(2,1) + xs(1,3)*ad(3,1)
        if (xsj.le.(0.d0)) then
c         write(iout,1000) nen,neg
          write(*,1000) nen,neg
          stop
        endif

      rxsj = 1.d0/xsj
c     Compute jacobian inverse
      do j = 1,3
        do i = 1,3
          xs(i,j) = ad(i,j)*rxsj
        end do
      end do

        do 55 n=1,nel
        do 55 i=1,3
        cartd1(i,n)=0.d0
        do 55 k=1,3
55        cartd1(i,n)=cartd1(i,n)+xs(i,k)*shp(k,n)


c     Initialize second derivatives with respect to global coords.
        do i=1,6
          do n=1,nel
            cartd2(i,n)=0.0d0
          end do
        end do

      if (bf) then

          temp(1) = bubble(1)
          temp(2) = bubble(2)
          temp(3) = bubble(3)
          do i = 1,3
            bubble(i) = temp(1)*xs(i,1) + temp(2)*xs(i,2)
     &                + temp(3)*xs(i,3)
          enddo

      endif


        IF(der) then
C---------------------------------------------------------------
C        Below is to calculate second derivative of shape functions

C        Calculate T2 matrix
        t2(1,1)=xs(1,1)**2
        t2(1,2)=xs(1,2)**2
        t2(1,3)=xs(1,3)**2
        t2(1,4)=2.d0*xs(1,1)*xs(1,2)
        t2(1,5)=2.d0*xs(1,2)*xs(1,3)
        t2(1,6)=2.d0*xs(1,3)*xs(1,1)

        t2(2,1)=xs(2,1)**2
        t2(2,2)=xs(2,2)**2
        t2(2,3)=xs(2,3)**2
        t2(2,4)=2.d0*xs(2,1)*xs(2,2)
        t2(2,5)=2.d0*xs(2,2)*xs(2,3)
        t2(2,6)=2.d0*xs(2,3)*xs(2,1)

        t2(3,1)=xs(3,1)**2
        t2(3,2)=xs(3,2)**2
        t2(3,3)=xs(3,3)**2
        t2(3,4)=2.d0*xs(3,1)*xs(3,2)
        t2(3,5)=2.d0*xs(3,2)*xs(3,3)
        t2(3,6)=2.d0*xs(3,3)*xs(3,1)

        t2(4,1)=xs(1,1)*xs(2,1)
        t2(4,2)=xs(1,2)*xs(2,2)
        t2(4,3)=xs(1,3)*xs(2,3)
        t2(4,4)=xs(1,1)*xs(2,2)+xs(1,2)*xs(2,1)
        t2(4,5)=xs(1,2)*xs(2,3)+xs(1,3)*xs(2,2)
        t2(4,6)=xs(1,1)*xs(2,3)+xs(1,3)*xs(2,1)

        t2(5,1)=xs(2,1)*xs(3,1)
        t2(5,2)=xs(2,2)*xs(3,2)
        t2(5,3)=xs(2,3)*xs(3,3)
        t2(5,4)=xs(2,1)*xs(3,2)+xs(2,2)*xs(3,1)
        t2(5,5)=xs(2,2)*xs(3,3)+xs(2,3)*xs(3,2)
        t2(5,6)=xs(2,1)*xs(3,3)+xs(2,3)*xs(3,1)

        t2(6,1)=xs(1,1)*xs(3,1)
        t2(6,2)=xs(1,2)*xs(3,2)
        t2(6,3)=xs(1,3)*xs(3,3)
        t2(6,4)=xs(1,1)*xs(3,2)+xs(1,2)*xs(3,1)
        t2(6,5)=xs(1,2)*xs(3,3)+xs(1,3)*xs(3,2)
        t2(6,6)=xs(1,1)*xs(3,3)+xs(1,3)*xs(3,1)

C        CALCULATE  T1 MATRIX
        do 60 i=1,6
        do 60 j=1,3
        tc(i,j)=0.d0
        do 60 k=1,6
        tc(i,j)=tc(i,j)-t2(i,k)*c1(k,j)
60        continue
        do 70 i=1,6
        do 70 j=1,3
        tcj(i,j)=0.d0
        do 70 k=1,3
        tcj(i,j)=tcj(i,j)+tc(i,k)*xs(k,j)
70        continue

        do 80 n=1,nel
         do 85 i=1,6
         at1(i)=0.d0
         do 85 k=1,3
         at1(i)=at1(i)+tcj(i,k)*shp(k,n)
85     continue
         do 90 i=1,6
         at2(i)=0.d0
         do 90 k=1,6
90         at2(i)=at2(i)+t2(i,k)*shp2(k,n)
         do 100 i=1,6
       cartd2(i,n)=at1(i)+at2(i)
100        continue
80        continue

       endif

 1000 format('1','non-positive determinant in element number  ',i5,
     &          ' in element material group  ',i5)

        return
      end
