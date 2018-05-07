c     ****************************************************************
c     *                                                              *
c     *                      subroutine oumodel_flat                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/15/2014 rhd              *
c     *                                                              *
c     *    output a flat file describing the model with minimum      *
c     *    info for use primarily by pat2exii to decrease runtime.   *
c     *    model size, coords, element info (type, material number   *
c     *    in warp3d input, incidences). element info can be Patran  *
c     *    convention or WARP3D convention. use text, stream or text *
c     *    compressed file formats                                   *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine oumodel_flat( user_file_name, text, stream, 
     &                         compressed, warp3d_convention, 
     &                         patran_convention  )
      use fft, only : incmap, incid
c     implicit integer (a-z)
      implicit none
      include 'common.main'
c
c         parameter declarations
c
      character (len=*) :: user_file_name
      logical :: text, stream, compressed, warp3d_convention,
     &           patran_convention 
c
c         local declarations
c
      character :: file_name*200, tmp_name*200,
     &             sdate_time_tmp*24, command*100, 
     &             convention_str*40
      character :: title*80, date_str*12, time_str*8
      logical   :: there, debug, connected, got_unit, ok  
      integer   :: next_char      
      integer, external  :: warp3d_get_device_number
c
c         implicit integer (a-z)
c
      integer :: out_file, start_extension, length_title, node, k
c
      debug = .false.
c
c        check to see that number of nodes & elements is > 0
c
      if( nonode .eq. 0 .or. noelem .eq. 0 ) then
          write(out,1020)
          return
      end if
c
c         get current time/date stamp for text file
c
      call fdate( sdate_time_tmp )
c
c        1.  find an available device number
c
      out_file = warp3d_get_device_number()
      got_unit = out_file .gt. 0
      if ( .not. got_unit ) then
        write(out,9400)
        call die_abort
      end if
c
c        2. strip extension if recognized. add .text or .str extension.
c           set warp3d or patran convention. open file
c         
      if( text .and. stream ) then
         write(out,1015)
         text   = .true.
         stream = .false.
      end if
c      
      if( warp3d_convention .and. patran_convention ) then
         write(out,1017)
         warp3d_convention = .true.
         patran_convention = .false.
      end if
c     
      file_name = " "    
      file_name(1:) = adjustl( user_file_name )
      start_extension = index( file_name, ".text", .true. ) ! at end
      next_char = index( file_name, " " ) 
      if ( start_extension .eq. 0 ) ! or .str extension
     &   start_extension = index( file_name, ".str", .true. ) ! at end
      if( start_extension .eq. 0 ) then ! no extension found, add one
        if( text )   file_name(next_char:) = ".text"
        if( stream ) file_name(next_char:) = ".str"
      end if
c
c        3. replace ~/ in file name with full path name to home
c           directory
c        
c     tmp_name = " "
c     tmp_name(1:) = file_name(1:)
c     file_name = " "
c     call tilde( tmp_name, file_name, ok ) ! tmp_name -> file_name
c     if( .not. ok ) then
c        write(out,1030) tmp_name(1:80)
c        call die_abort
c     end if
c
c        4. get the file open. 
c        
      inquire(file = file_name, exist = there )
      if( there ) write(out,1010)
      write(out,1000) file_name(1:len_trim(file_name))
      if( text ) open(unit=out_file,file=file_name,status='unknown')
      if( stream ) open(unit=out_file,file=file_name,status='unknown',
     &     access="stream", form="unformatted" )
c
c        5.  output some preface lines if a text file as comments
c            and model sizes
c
      title(1:) = "Structure: " // adjustl( stname )
      length_title = len_trim( title )
      date_str = sdate_time_tmp(5:11) // sdate_time_tmp(21:24)
      time_str = sdate_time_tmp(12:19)
      convention_str(1:) = "Patran element type and node ordering"
      if( warp3d_convention ) convention_str(1:) =
     &     "WARP3D element type and node ordering"  
c        
      if( text ) then 
        write(out_file,9000) title(1:length_title)
        write(out_file,9020) date_str, time_str
        write(out_file,9022) convention_str(1:)
        write(out_file,9010) nonode, noelem
      else
        write(out_file) nonode, noelem
      end if
c
c        6.  output coordinates of nodes
c
c              The global coordinates are stored in vector c such that
c              the x,y,z coordinates of a node are consecutive in the vector.
c
      if( text ) then 
        do node = 1, nonode
           k = 3 * node - 2
           write(out_file,9030) c(k), c(k+1), c(k+2)
        end do
      end if
      if( stream  ) then 
        do node = 1, nonode
           k = 3 * node - 2
           write(out_file) c(k), c(k+1), c(k+2)
        end do
      end if
c
c
c        7.  output element type, warp3d material number in the input
c            and incidences
c
c            we have 2 options: output element type and incidence in
c            warp3d convention or Patran convention
c
      if( patran_convention ) call oumodel_flat_patran
      if( warp3d_convention ) call oumodel_flat_warp3d
c
c        8. close unit
c
      close( unit= out_file )
      if( compressed .and. text ) then 
             command(1:) = ' '
             command(1:) = 'gzip ' // file_name
!win             result = system( command )
      end if
      write(out,1005)
c      
      return
c
 1000 format(1x,'>>>>>> Model description file: ',a)
 1005 format(1x,'>>>>>> Model description file completed')
 1010 format(//1x,'>>>>>> WARNING:  Specified model description file',
     &       ' already exists. '/ 18x,'File overwritten.')
 1015 format(//1x,'>>>>>> WARNING:  .text and .str both set.',
     &       /,1x,'                  Defaulting to .text ' )
 1017 format(//1x,'>>>>>> WARNING:  warp3d and patran conventions ',
     &            'specified',
     &       /,1x,'                  Defaulting to warp3d ',
     &          ' convention ' )
 1020 format(//,1x,'>>>>>> ERROR:  Either the number of nodes or',
     &        ' number of elements in the model is zero.'/
     &       18x, ' Must define model before using this option.' )
 1030 format(//,1x,'>>>>>> ERROR: cannot resolve/process file name:',
     &   /,18x, a80,/,18x,'Job aborted.'/) 
 9000 format("#",/,"#  ",a)
 9010 format( 2i9 )
 9020 format("#",/,"#  Created: ",a12,2x,a8,/,"#")
 9022 format("#  Convention: ",a40,/,"#")
 9030 format(3e25.13)
 9040 format(i4,i4,27i8)
 9400 format('>> FATAL: could not find a unit number to write',
     & /     '          model description file. job aborted...',/)
     
      contains
c     ****************************************************************
c     *                                                              *
c     *                subroutine oumodel_flat_patran                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/16/2016 rhd             *
c     *                                                              *
c     *    output element information to flat file using the Patran  *
c     *    convention for element type, node ordering.               *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine oumodel_flat_patran
      implicit none
c
c         local declarations
c
      integer, parameter :: emap_size=50
      integer :: warp_to_pat_20(20), pat_incid(30),
     &           warp_to_wed15(15), inter_tri12_to_wedge15(12),
     &           elem_types_to_exo_blk_map(emap_size)
      integer :: elem, incptr, etype, num_enodes, matnum,
     &           pat_etype, node, id, eleblk, span, ispan, felem,
     &           blk_knt_exo, blk_for_exo
      logical :: hex, tet, wedge,
     &           interface_hex,  interface_tri,
     &           interface_elem, ldum1, ldum2, ldum3, ldum4, ldum5
      data warp_to_pat_20 / 5,1,4,8,6,2,3,7,13,12,16,20,14,10,15,
     &                     18,17,9,11,19 /
      data warp_to_wed15  / 1,2,3,4,5,6,7,8,9,13,14,15,10,11,12/
      data inter_tri12_to_wedge15 / 1,2,3,7,8,9,4,5,6,13,14,15 /
c
c              for each element, we write the element type number
c              following the Patran scheme, block number for ParaView
c              (.exo) and the element incidences to follow the Patran 
c              numbering scheme.
c
c              EXO block number. Starts with 1 with max value = max
c              number of different element types appearing in the 
c              model. All elements of the same etype are assigned to 
c              the same EXO block. The EXO is not connected to the
c              WARP3D block number that groups elements for processing.
c
c              the code here builds a local vector as WARP3D blocks
c              are processed to store the EXO block numbers as they
c              are assigned.
c
c              at some point, we may want to make more EXO blocks
c              by assigning all element of the same etyp and material
c              to a block. However, this has the potential to create
c              a large number of EXO blocks. In ParaView, all those
c              blocks would need to checked on/off for visibilty.
c
      blk_knt_exo = 0
      elem_types_to_exo_blk_map = 0
      if( nelblk <= 0 ) then 
         write(out,9110)
         call die_abort
      end if
c   
      do eleblk = 1, nelblk
       span       = elblks(0,eleblk)
       felem      = elblks(1,eleblk)
c      etype      = iprops(1,felem)
       etype      = 2   ! l3disop, 8 nodes hex
       if( etype .gt. emap_size ) then
         write(out,9100)
         call die_abort
       end if   
       if( elem_types_to_exo_blk_map(etype). eq. 0 ) then
         blk_knt_exo = blk_knt_exo + 1
         elem_types_to_exo_blk_map(etype) = blk_knt_exo
       end if 
       blk_for_exo =  elem_types_to_exo_blk_map(etype) 
c      num_enodes  = iprops(2,felem)
       num_enodes  = 8    ! 13disop, 8 nodes hex
c      matnum      = iprops(38,felem) ! not used at present
       matnum      = 1 ! not used at present
       call set_element_type( etype, ldum1, hex, wedge, tet, ldum2,
     &                        ldum3, ldum4, ldum5, interface_elem )
       interface_hex = interface_elem .and. num_enodes .eq. 8
       interface_tri = interface_elem .and. 
     &                 (num_enodes .eq. 6 .or. num_enodes .eq. 12) 
c   
       do ispan = 1, span
         elem      = felem + ispan - 1
         pat_incid = 0
         incptr    = incmap(elem) - 1
c
         if( interface_hex ) then
            pat_etype = 8
            do node = 1, 8
               pat_incid(node) = incid(incptr + node)
            end do
            num_enodes = 8
         end if
c
         if( interface_tri ) then
            pat_etype = 7
            if( num_enodes .eq. 6 ) then  ! trint6
              do node = 1, 6
               pat_incid(node) = incid(incptr + node)
              end do
            else ! trint12, write as wedge15 with patran nodes
c                  10, 11, 12 missing
              num_enodes = 15
              pat_incid(1:15) = 0
              do node = 1, 12
                pat_incid(inter_tri12_to_wedge15(node)) =
     &                incid(incptr + node) 
              end do             
            end if
         end if
c
         if( hex ) then
            pat_etype = 8
            if( etype .eq. 1 ) then
              do node = 1, 20
                pat_incid(warp_to_pat_20(node)) = incid(incptr + node)
              end do
            else
              do node = 1, 8
               pat_incid(node) = incid(incptr + node)
              end do
              num_enodes = 8
            end if
         end if
c
         if( tet ) then  ! node scheme same
             pat_etype = 5
             do node = 1, num_enodes
                pat_incid(node) = incid(incptr + node)
             end do
         end if
c
         if( wedge ) then
             pat_etype = 7
             if( num_enodes .eq. 15 ) then
               do node = 1, num_enodes
                 pat_incid(warp_to_wed15(node)) = incid(incptr + node)
               end do
             else
               do node = 1, num_enodes
                 pat_incid(node) = incid(incptr + node)
               end do
             end if
         end if
c
         if( debug ) write(out,5000) elem,
     &        (incid(id), id=incptr+1,incptr+num_enodes),
     &        (pat_incid(id), id=1,num_enodes)
c
         if( text ) write(out_file,9040) 
     &              pat_etype, blk_for_exo, pat_incid(1:27)
         if( stream ) write(out_file)
     &              pat_etype, blk_for_exo, pat_incid(1:27)
c     
       end do ! over block loop for Patran elem type option
      end do ! over blocks 
c
      return    
c
 5000 format(4x,'>>>  Incidences for element : ',i6 //
     &        4x,'     WARP3D  : ',20i6 /
     &        4x,'     PATRAN  : ',20i6)
 9100 format('>> FATAL ERROR: vector overflow: oumodel_flat_patran',
     &  /,   '                job terminated' )
 9110 format('>> FATAL ERROR: element, coordinates, incidences, and',
     &  /,   '                blocking must be defined before this',
     &  /,   '                command. **job terminated**' )
 9040 format(i4,i4,27i8)
c        
      end subroutine oumodel_flat_patran
      
c     ****************************************************************
c     *                                                              *
c     *                      subroutine oumodel_flat_warp3d          *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/16/2016 rhd             *
c     *                                                              *
c     *    output element information to flat file using the WARP3D  *
c     *    convention for element type, node ordering.               *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine oumodel_flat_warp3d
      implicit none
c
c         local declarations
c
      integer, parameter :: emap_size=50
      integer :: wrp_incid(30), elem_types_to_exo_blk_map(emap_size)
      integer :: elem, incptr,etype, num_enodes, matnum, id, node,
     &           eleblk, span, ispan, felem,
     &           blk_knt_exo, blk_for_exo      
c 
c
c              for each element, we write the element type number
c              following the WARP3D scheme, block number for ParaView
c              (.exo) and the element incidences to follow the WARP3D 
c              numbering scheme.
c
c              EXO block number. Starts with 1 with max value = max
c              number of different element types appearing in the 
c              model. All elements of the same etype are assigned to 
c              the same EXO block. The EXO is not connected to the
c              WARP3D block number that groups elements for processing.
c
c              the code here builds a local vector as WARP3D blocks
c              are processed to store the EXO block numbers as they
c              are assigned.
c
c              at some point, we may want to make more EXO blocks
c              by assigning all element of the same etyp and material
c              to a block. However, this has the potential to create
c              a large number of EXO blocks. In ParaView, all those
c              blocks would need to checked on/off for visibilty.
c
      blk_knt_exo = 0
      elem_types_to_exo_blk_map = 0
      if( nelblk <= 0 ) then 
         write(out,9110)
         call die_abort
      end if
c      
      do eleblk = 1, nelblk
       span       = elblks(0,eleblk)
       felem      = elblks(1,eleblk)
c      etype      = iprops(1,felem)
       etype      = 2   ! l3disop, 8 nodes hex
       if( etype .gt. emap_size ) then
         write(out,9100)
         call die_abort
       end if   
       if( elem_types_to_exo_blk_map(etype). eq. 0 ) then
         blk_knt_exo = blk_knt_exo + 1
         elem_types_to_exo_blk_map(etype) = blk_knt_exo
       end if 
       blk_for_exo =  elem_types_to_exo_blk_map(etype) 
c      num_enodes  = iprops(2,felem)
       num_enodes  = 8    ! 13disop, 8 nodes hex
c      matnum      = iprops(38,felem) ! not used at present
       matnum      = 1 ! not used at present
c    
       do ispan = 1, span
         elem = felem + ispan - 1
         wrp_incid(1:30) = 0
         incptr        = incmap(elem) - 1
         do node = 1, num_enodes
            wrp_incid(node) = incid(incptr + node)
         end do
c
         if( debug ) write(out,5000) elem,
     &        (wrp_incid(id), id=1,num_enodes)
c
         if( text ) write(out_file,9040) 
     &              etype, blk_for_exo, wrp_incid(1:27)
         if( stream ) write(out_file)
     &               etype, blk_for_exo, wrp_incid(1:27)
       end do ! on ispan
c     
      end do ! element blocks loop 
c
      return  

 5000 format(4x,'>>>  Incidences for element : ',i6 //
     &        4x,'     WARP3D  : ',20i6 /
     &        4x,'     PATRAN  : ',20i6)
 9100 format('>> FATAL ERROR: vector overflow: oumodel_flat_warp3d',
     &  /,   '                job terminated' )
 9110 format('>> FATAL ERROR: element, coordinates, incidences, and',
     &  /,   '                blocking must be defined before this',
     &  /,   '                command. **job terminated**' )
 9040 format(i4,i4,27i8)
c
      end subroutine oumodel_flat_warp3d
      
      end subroutine oumodel_flat
c
c     ****************************************************************
c     *                                                              *
c     *         integer function warp3d_get_device_number()          *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 7/26/12                    *
c     *                                                              *
c     *     find a non-connected Fortran device number               *
c     *                                                              *
c     ****************************************************************
c
      integer function warp3d_get_device_number()
      implicit none
      integer iunit
      logical connected
c
c        1.  find an available device number, open the neutral
c            file
c
      warp3d_get_device_number = -1
      do iunit = 11, 5000
        inquire( unit=iunit, opened=connected )
        if ( connected ) cycle
        warp3d_get_device_number = iunit
        return
      end do
c
      return
      end
