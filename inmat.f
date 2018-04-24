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
