      subroutine inmat()
      use fft, only: matprp, lmtprp, imatprp, dmatprp, smatprp, 
     &               mat_assigned
      implicit none
      include 'common.main'

c                          local
      integer :: dum, nc, currmat, mattype
      real :: dumr
      real(8) :: dumd
      character lname*80, mname*24, dums*10, unknown*24
      logical, external :: matchs, matchs_exact, endfil, endcrd, numr
      logical, external :: label, warp3d_matl_num, scanms
c
c
c          read and store material name. If not exist, return
c
      if( .not. label(dum) ) then
        call errmsg(5,dum,dums,dumr,dumd)
        return
      end if
      
      lname  = ' '
      mname = ' '
      call entits( lname, nc )
      if( nc .gt. 24 ) nc = 24
      mname(1:nc) = lname(1:nc)

c                                                                               
c                       search for the specified material name in the           
c                       existing material library. if it is found,              
c                       print error message and return to the driver            
c                       subroutine because an existing material can-            
c                       not be overridden, only deleted and then re-            
c                       defined.                                                
c                                                                               
c                                                                               
         currmat = 1
         do while ( mat_assigned(currmat) )
           if( scanms(matnam(currmat),mname,24 ) ) then
             call errmsg(3,dum,mname,dumr,dumd)
             go to 666
           end if
           currmat = currmat + 1
         end do
c                                                                               
c                       not in library. check to make sure library is           
c                       not full.                                               
c                                                                               
         nummat = nummat + 1
         if( nummat.gt.mxmat ) then
            nummat = mxmat
            call errmsg(6,dum,dums,dumr,dumd)
            go to 666
         end if
c                                                                               
c                       find the first open slot in the material lib-           
c                       rary vector matnam and assign it to the spec-           
c                       ified material.                                         
c                                                                               
      mat_assigned(currmat) = .true.
      matnam(currmat) = mname
c
c
c          set default material parameters
c
      call mat_default(currmat)
c
c          read and store material properties. Material block
c          must start with PROPERTIES followed by material type
c
      call readsc()

      if ( .not. matchs('properties',6) ) then
        call errmsg(8,dum,dums,dumr,dumd)
        return
      end if
      if ( .not. label(dum) ) then
        call errmsg(8,dum,dums,dumr,dumd)
        return
      else
        lname = ' '
        call entits(lname,nc)
        mattype = warp3d_matl_num(lname,nc)
        matprp(9,currmat)  = mattype
      end if
c
c               cp material has special reading procedure
c
      if ( mattype .eq. 10 ) then
        call inmat_cp(currmat)
        go to 666
      end if
c
c               read materiap parameters
c
      do while ( .not. endcrd(dum) )
        if ( matchs_exact(',') ) then
          if ( endcrd(dum) ) call readsc()
        elseif ( matchs_exact('e') ) then
          if ( .not. numr( matprp(1,currmat) ) ) then
            call errmsg(7,dum,'e',dumr,dumd)
          end if
        elseif ( matchs_exact('nu') ) then
          if ( .not. numr( matprp(2,currmat) ) ) then
            call errmsg(7,dum,'nu',dumr,dumd)
          end if
        elseif ( matchs_exact('beta') ) then
          if ( .not. numr( matprp(3,currmat) ) ) then
            call errmsg(7,dum,'beta',dumr,dumd)
          end if
        elseif ( matchs_exact('tan_e') ) then
          if ( .not. numr( matprp(4,currmat) ) ) then
            call errmsg(7,dum,'tan_e',dumr,dumd)
          end if
        elseif ( matchs_exact('yld_pt') ) then
          if ( .not. numr( matprp(5,currmat) ) ) then
            call errmsg(7,dum,'yld_pt',dumr,dumd)
          end if
        elseif ( matchs_exact('alpha') ) then
          if ( .not. numr( matprp(6,currmat) ) ) then
            call errmsg(7,dum,'alpha',dumr,dumd)
          end if
        elseif ( matchs_exact('rho') ) then
          if ( .not. numr( matprp(7,currmat) ) ) then
            call errmsg(7,dum,'rho',dumr,dumd)
          end if
        else
          call entits(unknown,nc)
          call errmsg(7,dum,unknown(1:nc),dumr,dumd)
          call scan()
        endif
      end do
        
c             allocate space for material parameters
c             read and store material parameters
c

  666 return
      end subroutine
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine inmat_cp                     *          
c     *                                                              *          
c     *                       written by : mcm                       *          
c     *                                                              *          
c     *                   last modified : 3/21/12                    *          
c     *                                                              *          
c     *     input properties for the crystal plasticity model (#10)  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine inmat_cp(matnum)                                               
      use fft, only : matprp, lmtprp, imatprp, dmatprp, smatprp           
      implicit integer (a-z)                                                    
      integer, intent(in) :: matnum                                             
      integer :: dumi, nc                                                       
      real :: dumr                                                              
      double precision :: dumd                                                  
      character :: dums, lab*24, filen*24                                       
      logical :: reading                                                        
      logical, external :: matchs_exact, isstring, label, numr, numd,           
     &                     matchs, numi, endcrd                                 
                                                                                
c           See above for a detailed summary of each material option, in        
c           general, handle input for 6,7,9,13,22,25,100-113                    
c                                                                               
c           Set defaults                                                        
      lmtprp(13,matnum)=.false.                                                 
      lmtprp(22,matnum)=.true.                                                  
      lmtprp(25,matnum)=.false.                                                 
c           Read in properties                                                  
      reading = .true.                                                          
      do while (reading)                                                        
            if ( matchs_exact('angle_convention')) then                         
                  if (.not. label(dumi)) then                                   
                        call errmsg(356,dumi,dums,dumr,dumd)                    
                  else                                                          
                        lab = ' '                                               
                        call entits(lab,nc)                                     
                  end if                                                        
c                                                                               
                  if (lab(1:nc) .eq. 'kocks') then                              
                        imatprp(102,matnum) = 1                                 
                  else                                                          
                        call errmsg(357,dumi,lab(1:nc),dumr,dumd)               
                  end if                                                        
            elseif ( matchs_exact('alpha')) then                                
                  if (.not. numr(matprp(6,matnum))) then                        
                        call errmsg(23,dumi,'alpha',dumr,dumd)                   
                  end if                                                        
            elseif ( matchs_exact('rho')) then                                  
                  if (.not. numr(matprp(7,matnum))) then                        
                        call errmsg(23,dumi,'rho',dumr,dumd)                     
                  end if                                                        
            elseif ( matchs_exact('tolerance')) then                            
                  if (.not. numd(dmatprp(100,matnum))) then                     
                        call errmsg(23,dumi,'tolerance',dumr,dumd)               
                  end if                                                        
            elseif ( matchs_exact('n_crystals')) then                           
                  if (.not. numi(imatprp(101,matnum))) then                     
                        call errmsg(23,dumi,'n_crystals',dumr,dumd)              
                  end if                                                        
            elseif ( matchs_exact('angle_type')) then                           
                  if (.not. label(dumi)) then                                   
                        call errmsg(358,dumi,dums,dumr,dumd)                    
                  else                                                          
                        lab = ' '                                               
                        call entits(lab,nc)                                     
                  end if                                                        
c                                                                               
                  if (lab(1:nc) .eq. 'degrees') then                            
                        imatprp(103,matnum) = 1                                 
                  elseif (lab(1:nc) .eq. 'radians') then                        
                        imatprp(103,matnum) = 2                                 
                  else                                                          
                        call errmsg(359,dumi,lab(1:nc),dumr,dumd)               
                  end if                                                        
            elseif ( matchs_exact('crystal_input')) then                        
                  if (.not. label(dumi)) then                                   
                        call errmsg(360,dumi,dums,dumr,dumd)                    
                  else                                                          
                        lab = ' '                                               
                        call entits(lab,nc)                                     
                  end if                                                        
c                                                                               
                  if (lab(1:nc) .eq. 'single') then                             
                        imatprp(104,matnum) = 1                                 
                  elseif (lab(1:nc) .eq. 'file') then                           
                        imatprp(104,matnum) = 2                                 
                  else                                                          
                        call errmsg(361,dumi,lab(1:nc),dumr,dumd)               
                  end if                                                        
            elseif ( matchs_exact('crystal_type')) then                         
                  if (.not. numi(imatprp(105,matnum))) then                     
                        call errmsg(23,dumi,'crystal_type',dumr,dumd)            
                  end if                                                        
            elseif ( matchs_exact('orientation_input')) then                    
                  if (.not. label(dumi)) then                                   
                        call errmsg(360,dumi,dums,dumr,dumd)                    
                  else                                                          
                        lab = ' '                                               
                        call entits(lab,nc)                                     
                  end if                                                        
c                                                                               
                  if (lab(1:nc) .eq. 'single') then                             
                        imatprp(107,matnum) = 1                                 
                  elseif (lab(1:nc) .eq. 'file') then                           
                        imatprp(107,matnum) = 2                                 
                  else                                                          
                        call errmsg(361,dumi,lab(1:nc),dumr,dumd)               
                  end if                                                        
            elseif ( matchs_exact('angles') ) then                              
                  if (.not. (numd(dmatprp(108,matnum)) .and.                    
     &                       numd(dmatprp(109,matnum)) .and.                    
     &                       numd(dmatprp(110,matnum)))) then                   
                        call errmsg(362,dumi,dums,dumr,dumd)                    
                  end if                                                        
            elseif ( matchs_exact('filename')) then                             
                  call doscan()                                                 
                  if (.not. isstring()) then                                    
                        call errmsg(363,dumi,dums,dumr,dumd)                    
                  else                                                          
                        filen = ' '                                             
                        call entits(filen,nc)                                   
                        if (nc .gt. 24) then                                    
                              call errmsg(365,dumi,dums,dumr,dumd)              
                        end if                                                  
                        call scan()                                             
                        smatprp(112,matnum) = filen                             
                  end if                                                        
            elseif ( matchs_exact('debug')) then                                
                  if (.not. label(dumi)) then                                   
                        call errmsg(364,dumi,dums,dumr,dumd)                    
                  else                                                          
                        lab = ' '                                               
                        call entits(lab,nc)                                     
                        if (lab(1:nc) .eq. 'on') then                           
                              lmtprp(13,matnum) = .true.                        
                        elseif (lab(1:nc) .eq. 'off') then                      
                              lmtprp(13,matnum) = .false.                       
                        else                                                    
                              call errmsg(364,dumi,dums,dumr,dumd)              
                        end if                                                  
                  end if                                                        
            elseif ( endcrd(dum) ) then                                         
                  reading = .false.                                             
                  cycle                                                         
            elseif ( matchs(',',1) ) then                                       
                  call readsc()                                                 
            else                                                                
                  call entits(lab,nc)                                           
                  call errmsg(355,dumi,lab(1:nc),dumr,dumd)                     
                  call scan()                                                   
                  cycle                                                         
                  end if                                                        
            end do                                                              
                                                                                
      return                                                                    
      end subroutine                                                            
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine mat_default                  *          
c     *                                                              *          
c     *                       written by : RM                        *          
c     *                                                              *          
c     *                   last modified : 5/16/18                    *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mat_default(matnum)
      use fft, only: matprp, lmtprp, imatprp, dmatprp, smatprp
      implicit none
      include 'param_def'

      integer :: matnum

c                                                                               
c                       assign default material values. the current             
c                       ordering of material values is:                         
c                                                                               
c                  (*)  1  -- young's modulus                                   
c                  (*)  2  -- poisson's ratio                                   
c                       3  -- kinematic hardening ratio (beta)                  
c                       4  -- tangent modulus for bilinear strain               
c                             hardening (et)                                    
c                       5  -- inviscid yield stress (sigma-o)                   
c                  (*)  6  -- thermal expansion coefficient (isotropic)         
c                  (*)  7  -- mass density                                      
c                       8  -- linear elastic material flag                      
c                       9  -- material model type                               
c                              = 1 vectorized linear elastic                    
c                                  and invisicid plasticity model.              
c                                  isotropic/kinematic hardening                
c                              = 2 nonlinear elastic with linear +              
c                                  power law. small strains only.               
c                                  rate independent                             
c                              = 3 general gurson/mises model including         
c                                  nucleation, linear hardening,                
c                                  power-law hardening, matrix                  
c                                  viscoplasticity                              
c                              = 4 Cohesive zone models: linear elastic,        
c                                  bilinear, ramp, exponential_1 and            
c                                  exponential_2, ppr, cavitation               
c                              = 5 advanced cyclic plasticity model             
c                                  developed by kristine cochran                
c                              = 6 creep model                                  
c                              = 7 advanced mises model with hydrogen           
c                                  effects developed by yuiming liang           
c                              = 8 general UMAT for warp3d.                     
c                              = 10 (matching back up with file #s)             
c                                   CP model by mark messner                    
c                              = 11 interface damage model                      
c                      10  -- viscoplastic m power                              
c                      11  -- power-law hardening n power                       
c                      12  -- viscoplastic reference strain rate                
c                      13  -- debug material model computations                 
c                      14  -- initial porosity (f sub 0)                        
c                      15  -- gurson model parameter q1                         
c                      16  -- gurson model parameter q2                         
c                      17  -- gurson model parameter q3                         
c                      18  -- nucleation flag (true or false)                   
c                      19  -- nucleation parameter sn                           
c                      20  -- nucleation parameter en                           
c                      21  -- nucleation parameter fn                           
c                      22  -- allow material model to cut step due to           
c                      23  -- flag to allow crack growth element killing        
c                             excessive reversed plasticity                     
c                      24  -- segmental stress-strain curve logical flag        
c                      25  -- flag to indicate anisotropic thermal              
c                             coefficients are defined                          
c                      26-31  anisotropic thermal expansion coefficients        
c                      32  -- interface stiffness in the longitudinal direction 
c                      33  -- interface stiffness in the transverse direction   
c                      34  -- interface stiffness in the normal direction       
c                      35  -- critical normal stress of the interface           
c                      36  -- critical shear stress of the interface            
c                      37  -- shape parameter for bilinear and ramp             
c                      38  -- second ( additional) shape parameter for ramp     
c                      39  -- critical separation distance in sliding           
c                      40  -- critical separation distance in opening           
c                      41  -- equivalent critical separation distance           
c                      42  -- a ratio to determine the equivalent separation    
c                             under mixed mode loading ( =0 => mode I )         
c                      43  -- a flag for identifying the interface element      
c                      44  -- type of interface models                          
c                             1 - linear elastic; 2- bilinear                   
c                             3 - ramp; 4 - exponential_1; 5 - exponential_2    
c                      45  -- for segmental curve model, the segmental curve    
c                             set number                                        
c                  (*) 46  -- ductile material volume fracture for fgm          
c                      47  -- = 0 homogeneous cohesive material                 
c                             = 1 functionally graded cohesive material         
c                      48  -- critical separation distance ductile (fgm)        
c                      49  -- critical separation distance brittle (fgm)        
c                      50  -- critical stress ductile (fgm)                     
c                      51  -- critical stress brittle (fgm)                     
c                      52  -- beta_ductile (fgm)                                
c                      53  -- beta_brittle (fgm)                                
c                      54  -- compression stiffness multiplier for              
c                             cohesive materials                                
c                      55  -- start of props for cyclic plasticity model        
c                              see inmat_cyclic                                 
c                      70 -- start of props for yuemin liang's adv.             
c                             mises model that includes effects of              
c                             staturated hydrogen                               
c                              70 - yl_1                                        
c                              71 - yl_2                                        
c                              72 - yl_3                                        
c                              73 - yl_4                                        
c                              74 - yl_5                                        
c                              75 - yl_6                                        
c                              76 - yl_7                                        
c                              77 - yl_8                                        
c                              78 - yl_9                                        
c                              79 - yl_10                                       
c                                                                               
c                      80-89  creep model                                       
c                                                                               
c                      90-- PPR cohesive model                                  
c                              90-98 currently used. see comments               
c                              in inmat_cohesive routine                        
c                      100-- local tolerance for CP NR loop                     
c                      101-- number of crystals at e. material point (or max n) 
c                      102-- angle convention (1=Kocks)                         
c                      103-- angle type (1=degree, 2=radians)                   
c                      104-- crystal input (1=single, 2=file)                   
c                      105-- crystal number (for single)                        
c                      106-- crystal offset (for list)                          
c                      107-- orientation input (1=single, 2=file)               
c                      108-110-- psi, theta, phi (for single)                   
c                      111-- orientation offset (for list)                      
c                      112-- STRING crystal list (for offset/list)              
c                      113-- STRING orientation list (for offset/list)          
c                                                                               
c                      115-- macroscale material model number                   
c                      116-118-- s vector                                       
c                      119-121-- l vector                                       
c                      122-125-- t vector                                       
c                      126-- l_s                                                
c                      127-- l_l                                                
c                      128-- l_t                                                
c                      129-- alpha_dmg                                          
c                      130-- nstacks (temp, should calculate from element sz)   
c                      131-- nfail ("")                                         
c                       132-- macro_sz                                          
c                       133-- cp_sz                                             
c                                                                               
c                      148-- link2 x-stiffness
c                      149-- link2 y-stiffness
c                      150-- link2 z-stiffness
c                                                                              
c                      151-200 -- Abaqus compatible UMAT                        
c                              151 - um_1                                       
c                              151 - um_2                                       
c                              ...                                              
c                              ...                                              
c                              200 - um_50                                      
c                                                                               
c                      201-230 cavity option for cohesive material              
c                                                                               
c                      the beta_fact is used to assist in construction          
c                      of planar models containing one layer of elements.       
c                      the stiffness and internal forces are mutiplied          
c                      by the beta_fact to simulate the effect of a             
c                      non-unit thickness.                                      
c                                                                               
c                  (*) some material values maybe specified as having the       
c                      value 'fgm' (a string). In such cases the user           
c                      supplied nodal values for the model are interpolated     
c                      at the element gauss points during execution.            
c                                                                               
c                                                                               
c        beta_fact = one             
         matprp(1,matnum)  = 30000.0 
         matprp(2,matnum)  = 0.3     
         matprp(3,matnum)  = 1.0     
         matprp(4,matnum)  = 0.0     
         matprp(5,matnum)  = 0.0     
         matprp(6,matnum)  = 0.0     
         matprp(7,matnum)  = 0.0     
         lmtprp(8,matnum)  = .false. 
         matprp(9,matnum)  = 1       
         matprp(10,matnum) = 0.0     
         matprp(11,matnum) = 0.0     
         matprp(12,matnum) = 0.0     
         lmtprp(13,matnum) = .false. 
         matprp(14,matnum) = 0.0     
         matprp(15,matnum) = 1.5     
         matprp(16,matnum) = 1.0     
         matprp(17,matnum) = 2.25    
         lmtprp(18,matnum) = .false. 
         matprp(19,matnum) = 0.1     
         matprp(20,matnum) = 0.3     
         matprp(21,matnum) = 0.04    
         lmtprp(22,matnum) = .true.  
         lmtprp(23,matnum) = .false. 
         lmtprp(24,matnum) = .false. 
         lmtprp(25,matnum) = .false. 
         matprp(26,matnum) = 0.0     
         matprp(27,matnum) = 0.0     
         matprp(28,matnum) = 0.0     
         matprp(29,matnum) = 0.0     
         matprp(30,matnum) = 0.0     
         matprp(31,matnum) = 0.0     
         matprp(32,matnum) = 300000.0
         matprp(33,matnum) = 300000.0
         matprp(34,matnum) = 300000.0
         matprp(35,matnum) = 30000.0 
         matprp(36,matnum) = 11538.0 
         matprp(37,matnum) = .15     
         matprp(38,matnum) = .5      
         matprp(39,matnum) = 0.01    
         matprp(40,matnum) = 0.01    
         matprp(41,matnum) = 0.01    
         matprp(42,matnum) = 0.0     
         lmtprp(43,matnum) = .false. 
         matprp(46,matnum) = 0.0     
         matprp(47,matnum) = 0.0     
         matprp(48,matnum) = 0.0     
         matprp(49,matnum) = 0.0     
         matprp(50,matnum) = 0.0     
         matprp(51,matnum) = 0.0     
         matprp(52,matnum) = 0.0     
         matprp(53,matnum) = 0.0     
         matprp(54,matnum) = 10.0    
         matprp(55:99,matnum) = 0.0  
         dmatprp(100,matnum) = 1D-06 
         imatprp(101,matnum) = 1     
         imatprp(102,matnum) = 1     
         imatprp(103,matnum) = 1     
         imatprp(104,matnum) = 3     
         imatprp(105,matnum) = 0     
         imatprp(106,matnum) = 0     
         imatprp(107,matnum) = 3     
         dmatprp(108:110,matnum) = 0.0d00
         imatprp(111,matnum) = 0     
         smatprp(112,matnum) = 'ni'  
         smatprp(113,matnum) = 'ni'  
         matprp(114:mxmtpr,matnum) = 0.0          
         dmatprp(151:200,matnum) = 0.0d00   ! UMAT properties                   
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *              function  warp3d_matl_num                       *
c     *                                                              *
c     *       return the internal coed number for a material model   *
c     *               written by: rhd                                *
c     *                   last modified : 5/28/2016 rhd              *
c     *                                                              *
c     ****************************************************************
c
      integer function warp3d_matl_num( material_model_id, nc )
c      
c      use main_data, only : material_model_names
      implicit none
      include 'param_def'
c   
      integer :: nc, m
      character(len=nc) :: material_model_id      
c
      select case( material_model_id )
      
      case( "bilinear" )
       m = 1
      case( "deformation", "deform" )
       m = 2
      case( "mises", "gurson", "mises-gurson", "mises_gurson" )
       m = 3
      case( "cohesive" )
       m = 4
      case( "cyclic", "adv cyclic", "adv. cyclic" )
       m = 5
      case( "creep", "Norton", "norton" )
       m = 6
      case( "mises_hydrogen", "hydrogen" )
       m = 7
      case( "umat", "UMAT", "um" )
       m = 8
      case( "crystal", "CP", "cp", "crystal_plasticity" )
       m = 10
      case default
         write(*,9000) material_model_id
        call die_abort
      end select 
c
      warp3d_matl_num = m
      return
 9000 format('>> FATAL ERROR: routine warp3d_matl_num ',
     & /,    '                invalide material model name: ',a,
     & /,    '                job aborted',//)
      end      
