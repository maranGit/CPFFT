
c     ****************************************************************
c     *                                                              *
c     *                      subroutine ctran1                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 4/18/2016 rhd              *
c     *                                                              *
c     *     transform [Dt] from a form relating the unrotated stress *
c     *     rate and the unrotated rate of deformation tensor to     *
c     *     one relating the cauchy stress rate and the rate of      *
c     *     (spatial) deformation strain for a block of elements.    *
c     *                                                              *
c     ****************************************************************
c
      subroutine ctran1( span, felem, blk, cep, qn1, cs, qbar, dj, w,
     &                   is_umat, umat_stress_type, is_crys_pls, debug )
      implicit none
      include 'param_def'

c                 global
      integer, intent(in) :: span, felem, blk, umat_stress_type
      double precision :: cep(mxvl,nstr,*), qn1(mxvl,nstr,*), 
     &     cs(mxvl,*), dj(*), w
      logical :: qbar, is_umat, is_crys_pls, debug

c                 local
      double precision :: tc(mxvl,nstr,nstr), half, two, halfw, wf
      integer, parameter :: gpn = fftngp
      integer :: i, j
      logical :: do_transform
      data half, two / 0.5d00, 2.0d00 /
c
c             extract cep(mxvl,nstr,nstr) from mod_elebloks
      if(debug) write(*,*) "Entering extract_D_symmetric()"
      call extract_D_symmetric( gpn, span, blk, felem, cep )
      if(debug) write(*,*) "Leaving extract_D_symmetric()"
c      
c             [cep] (mxvl x 6 x 6) relates increments
c             of unrotated cauchy stress to increments
c             of the unrotated deformation. transform [cep]
c             so it relates increments of cauchy stress to
c             increments of the deformation, both on the
c             spatial coordinates.
c
c             [cep*] = [qn1] * [cep] * trans([qn1])
c
c             [qn1] is a rotation matrix constructed from the
c             [R] obtained by polar decomposition of the deformation
c             gradient, [F] =[R][U].
c
c             for UMATs the UMAT computed [cep] may already refer to
c             the Cauchy stress. No rotation to be done.
c
c             For crystal plasticity model,
c
      do_transform = .true.
      if( is_crys_pls ) do_transform = .true.
      if( is_umat ) then
        if( umat_stress_type .eq. 1 ) do_transform = .false.
      end if
      if( do_transform ) then
c
c             perform multiplication of [tc] = [qn1] * [cep]
c             use 6 as number of stress components to expose
c             value to compiler
c
      do j = 1, 6
         do i = 1, span
c
            tc(i,j,1)= (qn1(i,j,1)*cep(i,1,1)+
     &                  qn1(i,j,2)*cep(i,2,1)+
     &                  qn1(i,j,3)*cep(i,3,1)+
     &                  qn1(i,j,4)*cep(i,4,1)+
     &                  qn1(i,j,5)*cep(i,5,1)+
     &                  qn1(i,j,6)*cep(i,6,1))
c
            tc(i,j,2)= (qn1(i,j,1)*cep(i,1,2)+
     &                  qn1(i,j,2)*cep(i,2,2)+
     &                  qn1(i,j,3)*cep(i,3,2)+
     &                  qn1(i,j,4)*cep(i,4,2)+
     &                  qn1(i,j,5)*cep(i,5,2)+
     &                  qn1(i,j,6)*cep(i,6,2))
c
            tc(i,j,3)= (qn1(i,j,1)*cep(i,1,3)+
     &                  qn1(i,j,2)*cep(i,2,3)+
     &                  qn1(i,j,3)*cep(i,3,3)+
     &                  qn1(i,j,4)*cep(i,4,3)+
     &                  qn1(i,j,5)*cep(i,5,3)+
     &                  qn1(i,j,6)*cep(i,6,3))
c
            tc(i,j,4)= (qn1(i,j,1)*cep(i,1,4)+
     &                  qn1(i,j,2)*cep(i,2,4)+
     &                  qn1(i,j,3)*cep(i,3,4)+
     &                  qn1(i,j,4)*cep(i,4,4)+
     &                  qn1(i,j,5)*cep(i,5,4)+
     &                  qn1(i,j,6)*cep(i,6,4))
c
            tc(i,j,5)= (qn1(i,j,1)*cep(i,1,5)+
     &                  qn1(i,j,2)*cep(i,2,5)+
     &                  qn1(i,j,3)*cep(i,3,5)+
     &                  qn1(i,j,4)*cep(i,4,5)+
     &                  qn1(i,j,5)*cep(i,5,5)+
     &                  qn1(i,j,6)*cep(i,6,5))
c
            tc(i,j,6)= (qn1(i,j,1)*cep(i,1,6)+
     &                  qn1(i,j,2)*cep(i,2,6)+
     &                  qn1(i,j,3)*cep(i,3,6)+
     &                  qn1(i,j,4)*cep(i,4,6)+
     &                  qn1(i,j,5)*cep(i,5,6)+
     &                  qn1(i,j,6)*cep(i,6,6))
c
         end do
      end do
c
c                       perform multiplication of
c                       [cep*] =  [tc] * transpose([qn1])
c
      do j = 1, 6
         do i = 1, span
c
            cep(i,j,1)= tc(i,j,1)*qn1(i,1,1)+
     &                  tc(i,j,2)*qn1(i,1,2)+
     &                  tc(i,j,3)*qn1(i,1,3)+
     &                  tc(i,j,4)*qn1(i,1,4)+
     &                  tc(i,j,5)*qn1(i,1,5)+
     &                  tc(i,j,6)*qn1(i,1,6)
c
            cep(i,j,2)= tc(i,j,1)*qn1(i,2,1)+
     &                  tc(i,j,2)*qn1(i,2,2)+
     &                  tc(i,j,3)*qn1(i,2,3)+
     &                  tc(i,j,4)*qn1(i,2,4)+
     &                  tc(i,j,5)*qn1(i,2,5)+
     &                  tc(i,j,6)*qn1(i,2,6)
c
            cep(i,j,3)= tc(i,j,1)*qn1(i,3,1)+
     &                  tc(i,j,2)*qn1(i,3,2)+
     &                  tc(i,j,3)*qn1(i,3,3)+
     &                  tc(i,j,4)*qn1(i,3,4)+
     &                  tc(i,j,5)*qn1(i,3,5)+
     &                  tc(i,j,6)*qn1(i,3,6)
c
            cep(i,j,4)= tc(i,j,1)*qn1(i,4,1)+
     &                  tc(i,j,2)*qn1(i,4,2)+
     &                  tc(i,j,3)*qn1(i,4,3)+
     &                  tc(i,j,4)*qn1(i,4,4)+
     &                  tc(i,j,5)*qn1(i,4,5)+
     &                  tc(i,j,6)*qn1(i,4,6)
c
            cep(i,j,5)= tc(i,j,1)*qn1(i,5,1)+
     &                  tc(i,j,2)*qn1(i,5,2)+
     &                  tc(i,j,3)*qn1(i,5,3)+
     &                  tc(i,j,4)*qn1(i,5,4)+
     &                  tc(i,j,5)*qn1(i,5,5)+
     &                  tc(i,j,6)*qn1(i,5,6)
c
            cep(i,j,6)= tc(i,j,1)*qn1(i,6,1)+
     &                  tc(i,j,2)*qn1(i,6,2)+
     &                  tc(i,j,3)*qn1(i,6,3)+
     &                  tc(i,j,4)*qn1(i,6,4)+
     &                  tc(i,j,5)*qn1(i,6,5)+
     &                  tc(i,j,6)*qn1(i,6,6)
c
         end do
      end do
      end if ! on do_transform
c
c            Ran: no need to substract Q-bar for FFT program
c                 becaure I transfer directly from Green-Naghdi
c                 rate to dP/dF
c
c            subtract the [Q-bar] matrix from the transformed
c            [cep]. this is the "initial stress" at the material
c            point level. this remains an option indicated by qbar.
c            note: we must multiply in the gauss weight factor and
c            gauss point det[J] for the subtracted terms. the [cep]
c            passed in had these factors included by the cnst...
c            routines. the [Q-bar] 6x6 comes from the tensor
c            expression -2 (de.De):s, where, s is the stress tensor,
c            de is the rate of deformation tensor and De is the virtual
c            rate of deformation tensor. this expression in matrix form
c            is: - trans([B]) * [Q-bar] * [B]. this modification of 
c            [cep] is essential for convergence of nearly homogeneous
c            deformation problems.
c
c     if( qbar ) then
c       do i = 1, span
c        wf    = dj(i) * w
c        halfw = half * wf
c        cep(i,1,1) = cep(i,1,1) - two * cs(i,1) * wf
c        cep(i,2,2) = cep(i,2,2) - two * cs(i,2) * wf
c        cep(i,3,3) = cep(i,3,3) - two * cs(i,3) * wf
c        cep(i,4,1) = cep(i,4,1) - cs(i,4) * wf
c        cep(i,6,1) = cep(i,6,1) - cs(i,6) * wf
c        cep(i,4,2) = cep(i,4,2) - cs(i,4) * wf
c        cep(i,5,2) = cep(i,5,2) - cs(i,5) * wf
c        cep(i,5,3) = cep(i,5,3) - cs(i,5) * wf
c        cep(i,6,3) = cep(i,6,3) - cs(i,6) * wf
c        cep(i,4,4) = cep(i,4,4) - halfw * ( cs(i,1)+cs(i,2) )
c        cep(i,5,5) = cep(i,5,5) - halfw * ( cs(i,2)+cs(i,3) )
c        cep(i,6,6) = cep(i,6,6) - halfw * ( cs(i,1)+cs(i,3) )
c        cep(i,5,4) = cep(i,5,4) - halfw * cs(i,6)
c        cep(i,6,4) = cep(i,6,4) - halfw * cs(i,5)
c        cep(i,6,5) = cep(i,6,5) - halfw * cs(i,4)
c        cep(i,1,4) = cep(i,4,1)
c        cep(i,1,6) = cep(i,6,1)
c        cep(i,2,4) = cep(i,4,2)
c        cep(i,2,5) = cep(i,5,2)
c        cep(i,3,5) = cep(i,5,3)
c        cep(i,3,6) = cep(i,6,3)
c        cep(i,4,5) = cep(i,5,4)
c        cep(i,4,6) = cep(i,6,4)
c        cep(i,5,6) = cep(i,6,5)
c
c                      experiment with symmetrized version of the 
c                      nonsymmetric term. see Crisfield vol. 2, pg55.
c
c      z(1:6,1:6) = 0.0d0      
C      z(1,1:3) = cs(i,1)      
C      z(2,1:3) = cs(i,2)      
C      z(3,1:3) = cs(i,3)      
C      z(4,1:3) = cs(i,4)      
C      z(5,1:3) = cs(i,5)      
C      z(6,1:3) = cs(i,6)      
c      zt = transpose( z)
c      do k = 1, 6
c      do l = 1, 6
c        z(k,l) = 0.5d0*(z(k,l)+zt(l,k))
c        cep(i,k,l) = cep(i,k,l) + z(k,l)*wf
c      end do
c      end do
c         
c       end do
c     end if
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *              subroutine extract_D_symmetric                  *
c     *           - service routine. should be inlined -             *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified :  1/10/2016 rhd             *
c     *                                                              *
c     ****************************************************************
c
c
c     subroutine extract_D_symmetric( gpn, local_work )
      subroutine extract_D_symmetric( gpn, span, now_blk, felem, cep )
c
      use elem_block_data, only : gbl_cep_blocks => cep_blocks
      implicit none
      include 'param_def'
c     include 'include_tan_ek'
c
c                     global variables
c
      integer, intent(in) :: gpn, span, now_blk, felem
      real(8), intent(out) :: cep(mxvl,nstr,*)
c      
c                     local variables
c
      double precision :: weight, symm_part_cep(mxvl,21), f
      integer :: ielem, sloc, k 
      double precision, parameter :: one = 1.0D0
c      
c     span    = local_work%span
c     weight  = local_work%weights(gpn)
c     now_blk = local_work%blk
c     felem   = local_work%felem
c      
c              pull 21 terms of lower-triangle for this element
c              global cep block is 21 x span x num integration
c              points
c
c              expand to 6 x 6 symmetric [D] and scale by
c              integration weight factor
c
c      do ielem = 1, span
c        sloc = ( 21 * span * (gpn-1) ) + 21 * (ielem-1)
c        f = weight * local_work%det_jac_block(ielem,gpn)
c        do k = 1, 21
c         symm_part_cep(ielem,k) = f *
c     &        gbl_cep_blocks(now_blk)%vector(sloc+k)
c        end do
c      end do 
       
      do k = 1, 21
        do ielem = 1, span
        sloc = ( 21 * span * (gpn-1) ) + 21 * (ielem-1)
c       f = weight * local_work%det_jac_block(ielem,gpn)
        f = one
         symm_part_cep(ielem,k) = f *
     &        gbl_cep_blocks(now_blk)%vector(sloc+k)
        end do
      end do  

       do ielem = 1, span
        cep(ielem,1,1) = symm_part_cep(ielem,1)
        cep(ielem,2,1) = symm_part_cep(ielem,2)
        cep(ielem,2,2) = symm_part_cep(ielem,3)
        cep(ielem,3,1) = symm_part_cep(ielem,4)
        cep(ielem,3,2) = symm_part_cep(ielem,5)
        cep(ielem,3,3) = symm_part_cep(ielem,6)
        cep(ielem,4,1) = symm_part_cep(ielem,7)
        cep(ielem,4,2) = symm_part_cep(ielem,8)
        cep(ielem,4,3) = symm_part_cep(ielem,9)
        cep(ielem,4,4) = symm_part_cep(ielem,10)
        cep(ielem,5,1) = symm_part_cep(ielem,11)
        cep(ielem,5,2) = symm_part_cep(ielem,12)
        cep(ielem,5,3) = symm_part_cep(ielem,13)
        cep(ielem,5,4) = symm_part_cep(ielem,14)
        cep(ielem,5,5) = symm_part_cep(ielem,15)
        cep(ielem,6,1) = symm_part_cep(ielem,16)
        cep(ielem,6,2) = symm_part_cep(ielem,17)
        cep(ielem,6,3) = symm_part_cep(ielem,18)
        cep(ielem,6,4) = symm_part_cep(ielem,19)
        cep(ielem,6,5) = symm_part_cep(ielem,20)
        cep(ielem,6,6) = symm_part_cep(ielem,21)
        cep(ielem,1,2) = symm_part_cep(ielem,2)
        cep(ielem,1,3) = symm_part_cep(ielem,4)
        cep(ielem,1,4) = symm_part_cep(ielem,7)
        cep(ielem,1,5) = symm_part_cep(ielem,11)
        cep(ielem,1,6) = symm_part_cep(ielem,16)
        cep(ielem,2,3) = symm_part_cep(ielem,5)
        cep(ielem,2,4) = symm_part_cep(ielem,8)
        cep(ielem,2,5) = symm_part_cep(ielem,12)
        cep(ielem,2,6) = symm_part_cep(ielem,17)
        cep(ielem,3,4) = symm_part_cep(ielem,9)
        cep(ielem,3,5) = symm_part_cep(ielem,13)
        cep(ielem,3,6) = symm_part_cep(ielem,18)
        cep(ielem,4,5) = symm_part_cep(ielem,14)
        cep(ielem,4,6) = symm_part_cep(ielem,19)
        cep(ielem,5,6) = symm_part_cep(ielem,20)
      end do
c
      return
      end       
