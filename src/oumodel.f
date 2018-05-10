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
      logical, external :: label, isstring, string
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
      flat_name = ' '
      if ( label(dummy) ) then
        call entits( flat_name, nc )
      else if( string(dummy) ) then
        call entits( flat_name, nc )
      else
        call errmsg(19,dum,dums,dumr,dumd)
        flat_name(1:14) = 'FFT_model_flat'
      end if
c
c                 model size (reference configuration)
c
      debug = .false.
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

      if(debug) then
        do ii = 1,10
          write(out,*) x(1:3,ii)
        end do
      end if
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

      deallocate( x, ix )

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

      xll(1:3,1:27) = zero
      ixll(1:27) = 0

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

c     save

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
c           ixl(imid(i)) = imid(i)
            ixl(i) = i
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

c     save

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
c     *                   last modified by : RM                      *
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

c     save

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
