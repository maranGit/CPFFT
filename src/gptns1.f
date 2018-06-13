c     ****************************************************************
c     *                                                              *
c     *                      subroutine gptns1                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *               last modified : 8/20/2017 rhd                  *
c     *                                                              *
c     *     computes the contributon to the tangent                  *
c     *     stiffnes matrices for a block of similar elements in     *
c     *     uniform global coordinates for a single gauss point.     *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine gptns1( local_work, cep, qn1 )
c     use main_data, only: asymmetric_assembly
c     use elem_block_data, only : global_cep_blocks => cep_blocks
      implicit none
      include 'param_def'
      include 'include_sig_up'
c
c                         global
      real(8) :: cep(mxvl,nstr,*), qn1(mxvl,nstr,*)
c
c                         local
      integer :: gpn, span, local_iout, drive_cnst
      logical :: local_debug
c
      gpn         = local_work%gpn
      local_iout  = local_work%iout  
      drive_cnst  = local_work%mat_type
      span        = local_work%span
      local_debug = .false.
c
c
c                 branch on material type:
c                      1 = simple mises model- linear hardening                 
c                          (isotropic or mixed kinematic), has geonl            
c                          option                                               
c                          rate effects on flow properties. temperature         
c                          dependent flow properties, modulus, nu,              
c                          alpha.                                               
c                      2 = nonlinear elastic model (deformation                 
c                          plasticity                                           
c                          with linear + power-law stress-strain                
c                          relation. no geonl option.                           
c                      3 = general mises and gurson model -                     
c                          rate dependent/independent, linear (iso)             
c                          hardening or power-law hardening. gurson             
c                          model with w/o nucleation (strain                    
c                          controlled).                                         
c                          temperature dependent flow props, modulus,           
c                          nu, alpha                                            
c                      4 = interface constitutive models                        
c                          supports several cohesive zone models.               
c                          no geometric stiffness matrix. set geonl             
c                          false to bypass those additional                     
c                          computations                                         
c                      5 = adv. cyclic plasticity model                         
c                      6 = creep                                                
c                      7 = mises + hydrogen                                     
c                      8 = Abaqus UMAT                                          
c                     10 = CP model                                             
c               
c                    Always see funtion: warp3d_matl_num                        
c                    for the latest updates.                                    
c               
c                    [Dt] computed in unrotated configuration.                  
c               
c     if( local_work%is_link_elem ) drive_cnst = -1
c 
      select case( drive_cnst ) 
c     case( -1 )
c       continue   ! link elements                                                  
      case ( 1 )
        call drive_01_cnst( gpn, local_iout, local_work, cep )                       
c     case ( 2 )
c       call drive_02_cnst( gpn, local_iout, local_work )                       
c     case ( 3 )
c       call drive_03_cnst( gpn, local_iout, local_work )                       
c     case ( 4 )
c       call drive_04_cnst( gpn, local_iout, local_work )                       
c       geonl = .false.                                                         
c     case ( 5 )
c       call drive_05_cnst( gpn, local_iout, local_work )                       
c     case ( 6 )
c       call drive_06_cnst( gpn, local_iout, local_work )                       
c     case ( 7 )
c       call drive_07_cnst( gpn, local_iout, local_work )                       
c     case ( 8 )
c       call drive_umat_cnst( gpn, local_iout, local_work )                     
      case (10 )
        call drive_10_cnst( gpn, local_iout, local_work, cep )                       
      case default
          write(local_iout,9500)
          call die_abort
      end select
c
c     call ctran1( span, cep, qn1, local_work%is_umat, 
c    &             local_work%umat_stress_type,
c    &             local_work%is_crys_pls, local_debug )
c
      return
 9500 format(1x,'>> Fatal Error: gptns1. invalid material type..',              
     &    /, 1x,'                job terminated' )                              
      end subroutine
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
c     subroutine ctran1( span, felem, blk, cep, qn1, cs, qbar, dj, w,
c    &                   is_umat, umat_stress_type, is_crys_pls, debug )
      subroutine ctran1( span, cep, qn1, is_umat, umat_stress_type,
     &                    is_crys_pls, debug )
      implicit none
      include 'param_def'

c                 global
      integer, intent(in) :: span, umat_stress_type
      double precision :: cep(mxvl,nstr,*), qn1(mxvl,nstr,*)
c    &     cs(mxvl,*), dj(*), w
      logical :: qbar, is_umat, is_crys_pls, debug

c                 local
      double precision :: tc(mxvl,nstr,nstr), half, two, halfw, wf
      integer, parameter :: gpn = fftngp
      integer :: i, j
      logical :: do_transform
      data half, two / 0.5d00, 2.0d00 /
c
c             extract cep(mxvl,nstr,nstr) from mod_elebloks
c     if(debug) write(*,*) "Entering extract_D_symmetric()"
c     call extract_D_symmetric( gpn, span, blk, felem, cep )
c     if(debug) write(*,*) "Leaving extract_D_symmetric()"
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
c                                                                               
c                                                          
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine drive_01_cnst                *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 1/10/2016 rhd              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine drive_01_cnst( gpn, iout, local_work, cep )                         
c                                                                               
      implicit none                                                             
      include 'param_def'                                                       
c     include 'include_tan_ek'                                                  
      include 'include_sig_up'
c                                                                               
c                     parameter declarations                                    
c                                                                               
      integer :: gpn, iout                                                      
      real(8) :: cep(mxvl,nstr,*)
c                                                                               
c                     local variables                                           
c                                                                               
      integer :: span, now_blk, ielem, k, felem                                 
      logical :: local_debug                                                    
c                                                                               
      span    = local_work%span                                                 
      now_blk = local_work%blk                                                  
      felem   = local_work%felem                                                
c                                                                               
      local_debug = .false. ! felem .eq. 1 .and. gpn .eq. 3                     
c                                                                               
c     call extract_D_symmetric( gpn, local_work )                               
      call extract_D_symmetric( gpn, span, now_blk, felem, cep )
c                                                                               
      if( local_debug ) then                                                    
        write(iout,9000) now_blk, felem, gpn, span                              
        do ielem = 1, span                                                      
          write(iout,9010) felem + ielem - 1                                    
          do k = 1, 6                                                           
            write(iout,9020) cep(ielem,k,1:6)                        
          end do                                                                
        end do                                                                  
      end if                                                                    
      return                                                                    
c                                                                               
 9000 format(1x,'.... debug cnst1. now_blk, felem, gpn, span: ', 4i8)           
 9010 format(10x,'[D] for element: ',i7)                                        
 9020 format(15x,6es14.6)                                                       
                                                                                
      end                                                                       
c
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine drive_10_cnst                *          
c     *                                                              *          
c     *                       written by : mcm                       *          
c     *                                                              *          
c     *                   last modified : 12/8/2015 rhd              *          
c     *                                                              *          
c     *              drive [D] consistent update for CP model        *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine drive_10_cnst( gpn, iout, local_work, cepglb )                         
c                                                                               
c     use main_data, only : matprp, lmtprp                                      
c     use elem_block_data, only : gbl_cep_blocks => cep_blocks                  
      use mm10_defs, only : indexes_common, index_crys_hist                     
c                                                                               
      implicit none                                                             
      include 'param_def'                                                       
c                                                                               
c                     parameter declarations                                    
c                                                                               
c     include 'include_tan_ek'                                                  
      include 'include_sig_up'
      integer :: gpn, iout                                                      
      real(8) :: cepglb(mxvl,nstr,*)
c                                                                               
c                     local variables                                           
c                                                                               
      integer :: span, felem, iter, now_blk,                                    
     &           start_loc, k, i, j, eh, sh                                     
      double precision ::                                                       
     & weight, f, cep(6,6), cep_vec(36), tol                                    
      double precision, parameter :: zero = 0.0D0, one = 1.0D0
      logical :: local_debug                                                    
      equivalence ( cep, cep_vec )                                              
c                                                                               
      span             = local_work%span                                        
      felem            = local_work%felem                                       
c     weight           = local_work%weights(gpn)                                
      iter             = local_work%iter                                        
      now_blk          = local_work%blk                                         
      sh  = indexes_common(1,1) ! first index of cep tangent                    
      eh  = indexes_common(1,2) ! last index of cep tangent                     
c                                                                               
      local_debug =  .false.                                                    
      if( local_debug ) then                                                    
        write(iout,9100) now_blk, span, felem,                                  
     &                   local_work%hist_size_for_blk                           
      end if                                                                    
c                                                                               
c                     the consistent [D] for each integration point is          
c                     stored in the 1st 36 positions of history                 
c                                                                               
      do i = 1, span                                                            
c        f = weight * local_work%det_jac_block(i,gpn)                           
         f = one
         cep_vec(1:36) = local_work%elem_hist1(i,sh:eh,gpn)                     
         cepglb(i,1:6,1:6) = cep(1:6,1:6) * f                           
      end do                                                                    
c                                                                               
c                     code to optionally check symmetry of the                  
c                     [D] linear                                                
c                                                                               
      if( local_debug ) then                                                    
      tol = 0.01d00                                                             
      do i = 1, span                                                            
         cep_vec(1:36) = local_work%elem_hist1(i,sh:eh,gpn)                     
         do j = 1, 6                                                            
            if( cep(j,j) .lt. tol ) then                                        
               write(iout,*) ' .. fatal @ 1'                                    
               call die_abort                                                   
            end if                                                              
        end do                                                                  
        do j = 1, 6                                                             
          do k = 1, 6                                                           
            if( abs( cep(j,k) -cep(k,j) )                                       
     &           .gt. 1.0d-8 ) then                                             
               write(iout,*) ' .. fatal @ 2'                                    
               call die_abort                                                   
            end if                                                              
          end do                                                                
       end do                                                                   
      end do                                                                    
      end if                                                                    
c                                                                               
      if( local_debug ) then                                                    
        write(iout,*) ".... linear elastic [D] for CP"                          
        do k = 1, span                                                          
           write(iout,*) '    ... element: ', felem+k-1                         
           do i = 1, 6                                                          
             write(iout,9000) cep(i,1:6)                                        
           end do                                                               
        end do                                                                  
      end if                                                                    
      return                                                                    
c                                                                               
 9000 format(10x,6e14.5)                                                        
 9100 format(2x,".. debug for routine drive_10_cnst",                           
     & /,10x,"now_blk, span, felem, history size:", 4i7)                        
c                                                                               
      end                                                                       
