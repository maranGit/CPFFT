      subroutine inmat()
      use fft, only: mat_props
      implicit none
      include 'param_def'

c                          local
      integer :: dum, nc, currmat
      real :: dumr
      real(8) :: dumd
      character lname*80, mname*24, dums*10, unknown*24
      logical, external :: matchs, matchs_exact, endfil, endcrd, numd
      logical, external :: label, warp3d_matl_num
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

      currmat = 1
      do while ( mat_props(currmat)%assigned )
        currmat = currmat + 1
      end do
        
      if ( currmat .gt. mxmat ) call errmsg(6,dum,dums,dumr,dumd)
      
      mat_props(currmat)%assigned = .true.
      mat_props(currmat)%matnam = mname
c
c
c          read and store material properties
c
      call readsc()
      do while ( .not. endcrd(dum) )
        if ( matchs_exact(',') ) then
          if ( endcrd(dum) ) call readsc()
        elseif ( matchs('properties',6) ) then
          if (.not. label(dum)) then
            call errmsg(8,dum,dums,dumr,dumd)
            return
          else
            lname = ' '
            call entits(lname,nc)
            mat_props(currmat)%matnum = warp3d_matl_num(lname,nc)
          end if
        elseif ( matchs_exact('e') ) then
          if ( .not. numd( mat_props(currmat)%dmatprp(1) ) ) then
            call errmsg(7,dum,'e',dumr,dumd)
          end if
        elseif ( matchs_exact('nu') ) then
          if ( .not. numd( mat_props(currmat)%dmatprp(2) ) ) then
            call errmsg(7,dum,'nu',dumr,dumd)
          end if
        elseif ( matchs_exact('yld_pt') ) then
          if ( .not. numd( mat_props(currmat)%dmatprp(3) ) ) then
            call errmsg(7,dum,'yld_pt',dumr,dumd)
          end if
        elseif ( matchs_exact('tan_e') ) then
          if ( .not. numd( mat_props(currmat)%dmatprp(4) ) ) then
            call errmsg(7,dum,'tan_e',dumr,dumd)
          end if
        elseif ( matchs_exact('rho') ) then
          if ( .not. numd( mat_props(currmat)%dmatprp(5) ) ) then
            call errmsg(7,dum,'rho',dumr,dumd)
          end if
        elseif ( matchs_exact('alphax') ) then
          if ( .not. numd( mat_props(currmat)%dmatprp(6) ) ) then
            call errmsg(7,dum,'alphax',dumr,dumd)
          end if
        elseif ( matchs_exact('alphay') ) then
          if ( .not. numd( mat_props(currmat)%dmatprp(7) ) ) then
            call errmsg(7,dum,'alphay',dumr,dumd)
          end if
        elseif ( matchs_exact('alphaz') ) then
          if ( .not. numd( mat_props(currmat)%dmatprp(8) ) ) then
            call errmsg(7,dum,'alphaz',dumr,dumd)
          end if
        elseif ( matchs_exact('beta') ) then
          if ( .not. numd( mat_props(currmat)%dmatprp(9) ) ) then
            call errmsg(7,dum,'beta',dumr,dumd)
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
                        call errmsg(5,dumi,'alpha',dumr,dumd)                   
                  end if                                                        
            elseif ( matchs_exact('rho')) then                                  
                  if (.not. numr(matprp(7,matnum))) then                        
                        call errmsg(5,dumi,'rho',dumr,dumd)                     
                  end if                                                        
            elseif ( matchs_exact('tolerance')) then                            
                  if (.not. numd(dmatprp(100,matnum))) then                     
                        call errmsg(5,dumi,'tolerance',dumr,dumd)               
                  end if                                                        
            elseif ( matchs_exact('n_crystals')) then                           
                  if (.not. numi(imatprp(101,matnum))) then                     
                        call errmsg(5,dumi,'n_crystals',dumr,dumd)              
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
                        call errmsg(5,dumi,'crystal_type',dumr,dumd)            
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
