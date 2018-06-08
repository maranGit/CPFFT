c     ****************************************************************
c     *                                                              *
c     *                      subroutine ouresult                     *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified : 5/07/2018 RM               *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ouresult(step)
      implicit none
c
c                 global
c
      integer :: step
c
c                 print nodal displacement
c
      call oudisp( step )
c
c                 print element stress/strain
c
      call ouelestr( step, 1 ) ! strain
      call ouelestr( step, 2 ) ! stress
c
c                 print element state variables
c
      call oustates()

      return

      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine oudisp                       *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified : 5/11/2018 RM               *
c     *                                                              *
c     *                drive nodal displacement output               *
c     *                                                              *
c     ****************************************************************
c
      subroutine oudisp(step)
      use fft, only: N
      implicit none
      include 'common.main'
c
c                 global
c
      integer :: step
c
c                 local
c
      integer :: dva, bnfile, fmfile, opt, nod, currdof
      integer :: flat_file_number, outype, quantity
      logical :: oubin, ouasc, flat_file, stream_file
      logical :: text_file, compressed
      real(8), parameter :: small_tol=1.0d-30, zero=0.0d0
c
c                 compute nodal displacement from F
c                 store in common.main -> u
c
      call f2disp(step)
c
c                 hard code output options
c
      opt         = 1
      dva         = 1
      bnfile      = 0
      fmfile      = 0
      myid        = 0
      oubin       = .false.
      ouasc       = .false.
      use_mpi     = .false.
      flat_file   = .true.
      text_file   = .true.
      stream_file = .false.
      compressed  = .false.
c
c                 open output file
c
      call ouocdd( dva, step, oubin, ouasc, bnfile, fmfile,
     &             opt, use_mpi, myid, flat_file, stream_file,
     &             text_file, compressed, flat_file_number )

c
c                 print nodal displacement
c                 copied from warp3d ouddpa.f
c
      outype      = 1
      quantity    = 1
      call ouddpa_flat_header(outype, quantity, flat_file_number)

      where( abs( u ) .lt. small_tol ) u = zero

      do nod = 1, nonode
        currdof = 3*nod - 2
        if( flat_file .and. text_file ) 
     &    write(flat_file_number,930) u(currdof:(currdof+2))
 
      end do  ! on all nodes
c
c                 close output file
c
      call ouocdd( dva, step, oubin, ouasc, bnfile, fmfile, 2,
     &             .false., myid, flat_file, stream_file, text_file, 
     &             compressed, flat_file_number )

      return

 930  format(3e15.6)

      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine ouelestr                     *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified : 5/11/2018 RM               *
c     *                                                              *
c     *                 drive element stress/strain output           *
c     *                                                              *
c     ****************************************************************
c
      subroutine ouelestr(step, data_type)
      use elem_block_data, only: urcs_n1_blocks, eps_n1_blocks
      implicit none
      include 'common.main'
c
c                 global
c
      integer :: step, data_type
c
c                 local
c
      integer :: flat_file_number, fileno, quantity
      logical :: flat_file, text_file, stream_file
      logical :: oubin, ouasc, compressed
      character(len=1) dummy_char
      real(8), parameter :: small_tol=1.0d-30, zero=0.0d0
      logical :: stress
      double precision :: elem_results(mxvl,mxstmp)
      integer, external :: warp3d_get_device_number
      integer :: num_short_strain
      integer :: num_short_stress
      integer :: blk, span, felem, num_vals
      integer :: i
c
c                 hard code output options
c
c     ltmstp      = step
      myid        = 0
      oubin       = .false.
      ouasc       = .false.
      flat_file   = .true.
      text_file   = .true.
      stream_file = .false.
      compressed  = .false.
      stress      = .false.
      num_short_stress = 11
      num_short_strain = 7
      if(data_type .eq. 2) stress = .true.
      num_vals = num_short_strain + 15
      if( stress ) num_vals = num_short_stress + 15
c
c                        open output file
c
      call ouocst_elem( data_type, ltmstp, oubin, ouasc, fileno, 1,
     &                  use_mpi, myid, flat_file,
     &                  stream_file, text_file, compressed,
     &                  flat_file_number, dummy_char )
c
c                        header info files
c
      if( flat_file .and. text_file ) then
        quantity = 1
        if( stress ) quantity = 2
        call ouddpa_flat_header( 3, quantity, flat_file_number )
      else
        call errmsg(22)
        call die_abort()
      end if
c
c                        write element data records into Patran or
c                        flat file. For Patran file(s), the element 
c                        type is hard coded into this line as a 
c                        hex element
c
c                        zero small values to prevent 3-digit exponents
c                        in formatted files
c
c     copy stress/strain from pointer array to local stack
      if( stress ) then ! output stress

        do blk = 1, nelblk

          elem_results = zero
          span = elblks(0, blk)
          felem = elblks(1, blk)

          ! ngp = 1, so easier
c         format of urcs_n1_blocks: nstrs * ngp * mxvl
          do i = 1, mxvl
            elem_results(i,1:nstr) 
     &       = urcs_n1_blocks(blk)%ptr((i*nstrs-nstrs+1):(i*nstrs-3))
          end do

          call oupestr(span, num_vals, elem_results(1,1), 
     &                 flat_file_number)

        end do ! loop over block

      else ! output strain
        do blk = 1, nelblk

          elem_results = zero
          span = elblks(0, blk)
          felem = elblks(1, blk)

          ! ngp = 1, so easier
c         format of eps_n1_blocks: nstr * ngp * mxvl
          do i = 1, mxvl
            elem_results(i,1:nstr)
     &       = eps_n1_blocks(blk)%ptr((i*nstr-nstr+1):(i*nstr))
          end do

          call oupestr(span, num_vals, elem_results(1,1),
     &                 flat_file_number)

        end do ! loop over block

      end if
c
c                        close output file
c
      call ouocst_elem( data_type, ltmstp, oubin, ouasc, fileno, 2,
     &                  use_mpi, myid, flat_file,
     &                  stream_file, text_file, compressed,
     &                  flat_file_number, dummy_char )
      return

      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine oupestr                      *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified : 5/12/2018 RM               *
c     *                                                              *
c     *            print element stress/strain to result file        *
c     *                  format is flat text file                    *
c     *                                                              *
c     ****************************************************************
c
      subroutine oupestr(span, n, str, fid)
      implicit none
      include 'common.main'
c                        global
      integer :: span, n, fid
      real(8) :: str(mxvl, mxstmp)
c                        local
      integer :: elem
      real(8), parameter :: small_tol=1.0d-30, zero=0.0d0

c     zero small values
      where( abs(str) .lt. small_tol ) str = zero

c     write to output file ( flat text )
         do elem = 1, span
           if( use_mpi ) then
             write(fid,9100) elem, str(elem,1:n)
           else
             write(fid,9200) str(elem,1:n)
           end if
         end do

      return
 9100 format(i8,30d15.6)
 9200 format(30e15.6)
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine oustates                     *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified : 5/11/2018 RM               *
c     *                                                              *
c     *                 print element state results                  *
c     *                 copied from warp3d oustates.f                *
c     *                 current configuration                        *
c     *                                                              *
c     ****************************************************************
c
      subroutine oustates
      implicit none
c
c                 print element state results
c                 copied from warp3d oustates.f
c                 current configuration
c
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *             subroutine ouddpa_flat_header                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 4/16/2014 (rhd)            *
c     *                                                              *
c     *     Write header lines for flat text files                   *
c     *                                                              *
c     ****************************************************************
      
      subroutine ouddpa_flat_header( type, quantity, 
     &                               flat_file_number )
      implicit none
      include 'common.main'

      integer :: type, quantity, flat_file_number
c
c                       local declarations
c
      integer :: step_num
      character(len=15) :: vector_types(5), tensor_types(2),
     &                     scalar_types(1) 
      character(len=24) :: sdate_time_tmp
c      
      data vector_types / 'displacements', 'velocities',
     &                    'accelerations', 'reactions',
     &                    'temperatures' / 
      data tensor_types / 'strains', 'stresses' /
      data scalar_types / 'states' /
c
c                       type:
c                         = 1 a nodal quantity of vector type
c                         = 2 a nodal tensor type
c                         = 3 an element tensor type
c                         = 4 scalar type
c                       quantity for type = 1, see entries in
c                         vector_types above
c                       quantity for type = 2,3 see entries in
c                         tensor_types
c                       quantity for type = see scalar_types
c
      step_num = ltmstp  ! from common.main
c      
      write(flat_file_number,9000)
      if( type .eq. 1 ) then
        write(flat_file_number,9010) vector_types(quantity)
      end if
      if( type .eq. 2 ) then
        write(flat_file_number,9010) tensor_types(quantity)
      end if
      if( type .eq. 3 ) then
        write(flat_file_number,9012) tensor_types(quantity)
      end if
      if( type .eq. 4 ) then
        write(flat_file_number,9014) scalar_types(quantity)
      end if
      write(flat_file_number,9020) stname
      write(flat_file_number,9030) nonode, noelem
      call fdate( sdate_time_tmp )
      write(flat_file_number,9040) sdate_time_tmp
      write(flat_file_number,9050) step_num
      write(flat_file_number,9000)
c      
      return 
c            
 9000 format('#')
 9010 format('#  WARP3D nodal results: ',a15)
 9012 format('#  WARP3D element results: ',a15)
 9014 format('#  WARP3D element results: ',a15)
 9020 format('#  Structure name: ',a8 )
 9030 format('#  Model nodes, elements: ',2i8)
 9040 format('#  ',a24)
 9050 format('#  Load(time) step: ',i8 )
c 
      end 

c     ****************************************************************
c     *                                                              *
c     *                      subroutine ouocdd                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 1/21/2017 rhd              *
c     *                                                              *
c     *     handle file open/close for Patran and flat file results  *
c     *     make file name based on load step, results type, etc.    *
c     *     this routine handles result types:                       *
c     *     displacements, velocities, accelerations, reactions,     *
c     *     temperatures                                             *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine ouocdd( dva, ltmstp, oubin, ouasc, bnfile, fmfile,
     &                   opt, use_mpi, myid, flat_file, stream_file,
     &                   text_file, compressed, flat_file_number )
      implicit none
c
      integer :: dva, ltmstp, opt, myid, bnfile, fmfile, 
     &           flat_file_number
      logical :: oubin, ouasc, use_mpi, flat_file, 
     &           stream_file, text_file, compressed 
c
      integer :: k, step_number
      integer, external :: warp3d_get_device_number
      logical :: patran_file, ok
      character(len=14) :: bflnam, fflnam
      character(len=4)  :: strtnm, slist(5)
      character(len=3)  :: flat_list(5)
      character(len=30) :: command
      character(len=30), save :: flat_name
!win      external system       
      
      data slist / 'wn?d','wn?v', 'wn?a', 'wn?r', 'wn?t' /
      data flat_list /  'wnd','wnv', 'wna', 'wnr', 'wnt' / 
c
c                       close patran or flat result file
c                       compress for flat text file if option
c                       requested.
c
      patran_file = oubin .or. ouasc
      ok = patran_file .or. flat_file
c      
      if( .not. ok ) then
        k = 0
        call errmsg( 20 )
        call die_abort
      end if
c                     
      if( opt .eq. 2 .and. patran_file ) then  ! close file
        if( oubin ) close(bnfile,status='keep')
        if( ouasc ) close(fmfile,status='keep')
        return
      end if
c      
      if( opt .eq. 2 .and. flat_file ) then
         close(unit=flat_file_number,status='keep')
         if( .not. text_file ) return
          if( compressed ) then
             command(1:) = ' '
             command(1:) = 'gzip ' // flat_name
!win             result = system( command )
          end if
          return
      end if
c
c                       opt =1, create file name and open file
c                       valid request for here?
c
      if( dva .lt. 1 .or. dva .gt. 5 ) then
        k = 0
        call errmsg( 21 )
        call die_abort
      end if
c
c                       patran result files. ltmstp is step number
c
c                         wn(b)(f) +    X     + step no + MPI rank
c                        char*3      char*1      i5.5     i4.4
c
c                        X = d, v, a, r, t
c
      step_number = ltmstp
      if ( patran_file ) then      
        strtnm = slist(dva)
        if( oubin ) then
          strtnm(3:3) = "b"
          call ouflnm( strtnm, bflnam, step_number, use_mpi, myid )
          bnfile = 98
          open(bnfile,file=bflnam,status='unknown',
     &         access='sequential',form='unformatted',recl=350 )
        end if
        if( ouasc)  then
          strtnm(3:3) = "f"
          call ouflnm( strtnm, fflnam, step_number, use_mpi, myid )
          fmfile = 99
          open(fmfile,file=fflnam,status='unknown',
     &         access='sequential',form='formatted',recl=350 )
         end if
         return
      end if         
c
c                       flat result files. name structure
c                        wnX + step # + _text   + .MPI rank
c                        wnX + step # + _stream + .MPI rank
c                              i5.5                i4.4
c                        where X is d, v, a, r, t
c
      flat_file_number = warp3d_get_device_number()
      flat_name(1:30) = ' '
      flat_name(1:3)  = flat_list(dva)
      if( step_number .gt. 99999 ) step_number = step_number - 99999
      write(flat_name(4:),9000) step_number
c      
      if( stream_file ) then
        flat_name(9:) = '_stream'
        if( use_mpi ) then
          flat_name(16:16) = "."
          write(flat_name(17:),9100) myid  
        end if    
        open( unit=flat_file_number, file=flat_name, status='unknown',
     &        access='stream', form='unformatted' )
        return  
      end if
c      
      if( text_file ) then
        flat_name(9:) = '_text'
        if( use_mpi ) then
          flat_name(14:14) = "."
          write(flat_name(15:),9100) myid  
        end if  
        open( unit=flat_file_number, file=flat_name, status='unknown',
     &        access='sequential', form='formatted' )
        return
      end if         
c
 9000 format(i5.5)
 9100 format(i4.4)
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine ouocst_elem                  *
c     *                                                              *
c     *                       written by : kck                       *
c     *                                                              *
c     *                   last modified : 2/2/2017 rhd               *
c     *                                                              *
c     *     this subroutine opens or closes files for (1) Patran     *
c     *     binary or formatted output, or (2) flat file with        *
c     *     stream or text form                                      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ouocst_elem( data_type, stepno, oubin, ouasc, fileno, 
     &                        opt, use_mpi, myid, flat_file,
     &                        stream_file, text_file, compressed,
     &                        flat_file_number, matl_name_id  )
      implicit none
c   
c                       parameters
c
      integer :: data_type, stepno, fileno, opt,  myid,
     &           flat_file_number, dot_pos 
      logical :: oubin, ouasc, use_mpi,
     &           flat_file, stream_file, text_file, compressed
      character(len=*) :: matl_name_id
c   
c                       locals
c
      integer :: now_len
      logical :: ok, patran_file
      character(len=80) :: command
      character(len=80), save :: flat_name
c
c                       branch on whether files are to be opened or  
c                       closed.
c
c                       data_type: 1 = strains, 2 = stresses, 
c                                  3 = material states
c
c                       file_no is only for Patran results. calling 
c                       routine must pass value.
c
c                       flat_file_number. routine here sets value
c                       on file open. calling routine must send
c                       value on close
c
      patran_file = oubin .or. ouasc
c      
      select case( opt )
c
      case( 1 )   ! open files
        if( patran_file ) call ouocst_elem_pat_file
        if( flat_file )   call ouocst_elem_flat_file
        ok = patran_file .or. flat_file
        if( .not. ok ) then
          write(*,9000) opt
          call die_abort
        end if  
c
      case( 2 )  !  close files
        ok = patran_file .or. flat_file
        if( .not. ok ) then
          write(*,9000) opt
          call die_abort
        end if  
c
        if( patran_file ) then
          close( unit=fileno, status='keep' )
          return
        end if
c      
        close( unit=flat_file_number, status='keep' )
        if( stream_file ) return
        if( compressed ) then
          command(1:) = ' '
          now_len = len_trim( flat_name )
          command(1:) = 'gzip ' // flat_name(1:now_len)
!win          result = system( command )
        end if
c          
      end select
      
      return
 9000 format('>> FATAL ERROR: ouocst_elem. opt: ',i2,
     & /,    '                job aborted',// )      
c
      contains
c     ========  
c
c     ****************************************************************
c     *  (in contains)  subroutine ouocst_elem_pat_file              *
c     ****************************************************************
c
      subroutine ouocst_elem_pat_file
      implicit none
c
      integer :: now_len
      character :: patran_file_name*80, strtnm*4, form_type*20     
c
c                       attach patran binary or ascii file
c                       
      if( oubin ) then
         strtnm = 'webe'
         if( data_type .eq. 2 ) strtnm = 'webs'
         if( data_type .eq. 3 ) strtnm = 'webm'
         form_type = "unformatted"
      else ! ascii
         strtnm = 'wefe'
         if( data_type .eq. 2 ) strtnm = 'wefs'
         if( data_type .eq. 3 ) strtnm = 'wefm'
         form_type = "formatted"
      end if
c      
      patran_file_name = " "
      patran_file_name(1:4) = strtnm
      write(patran_file_name(5:),9000) stepno
c      
      if( data_type .eq. 3 .and. matl_name_id(1:1) .ne. " " )
     &     patran_file_name(10:) = "_" // matl_name_id(1:)
c     
      if( use_mpi ) then
         dot_pos = len_trim( patran_file_name ) + 1
         patran_file_name(dot_pos:dot_pos) = "."
         write(patran_file_name(dot_pos+1:),9010) myid
      end if  
c
      open(unit=fileno,file=patran_file_name,status='unknown',
     &     access='sequential', form=form_type )
c
      return
c      
 9000 format( i5.5 )
 9010 format( i4.4 )
c
      end subroutine ouocst_elem_pat_file
c
c
c     ****************************************************************
c     *  (in contains)  subroutine ouocst_elem_flat_file             *
c     ****************************************************************
c
      subroutine ouocst_elem_flat_file
      implicit none
c
      integer :: now_len, dot_pos
      integer, external :: warp3d_get_device_number
!win      integer, external  :: system
      character :: strtnm*4, form_type*20, access_type*20
c      
      flat_file_number = warp3d_get_device_number()
      if( flat_file_number .eq. -1 ) then
        write(*,9110)
        call die_gracefully
      end if
c
      flat_name(1:) = ' '
      flat_name(1:3) = 'wee' 
      if( data_type .eq. 2 ) flat_name(1:3) ='wes'
      if( data_type .eq. 3 ) flat_name(1:3) ='wem'
      write(flat_name(4:),9000) stepno
c      
      if( flat_file .and. text_file ) then
         flat_name(9:) = '_text'
         if( data_type .eq. 3 .and. matl_name_id(1:1) .ne. " " ) then
            now_len = len_trim( flat_name )
            flat_name(now_len+1:) = "_" // matl_name_id(1:)
         end if  
         if( use_mpi ) then
           dot_pos = len_trim( flat_name ) + 1
           flat_name(dot_pos:dot_pos) = "."
           write(flat_name(dot_pos+1:),9010) myid
         end if  
         access_type = 'sequential'
         form_type   = 'formatted'
      end if
c
      if( flat_file .and. stream_file ) then
         flat_name(9:) = '_stream'
         if( data_type .eq. 3 .and. matl_name_id(1:1) .ne. " " ) then
            now_len = len_trim( flat_name )
            flat_name(now_len+1:) = "_" // matl_name_id(1:)
         end if 
         if( use_mpi ) then
           dot_pos = len_trim( flat_name ) + 1
           flat_name(dot_pos:dot_pos) = "."
           write(flat_name(dot_pos+1:),9010) myid
         end if  
         access_type = 'stream'
         form_type   = 'unformatted'
      end if

      open( unit=flat_file_number,file=flat_name, status='unknown',
     &     access=access_type, form=form_type ) 
c
      return      
c
 9000 format( i5.5 )
 9010 format( i4.4 )
 9110 format(/1x,
     &'>>>>> FATAL ERROR: routine ouocst_elem_flat_file',
     & /,16x,'Job terminated....'/)
c
      end subroutine ouocst_elem_flat_file
      end subroutine ouocst_elem
c     ****************************************************************
c     *                                                              *
c     *                      subroutine ouflnm                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 04/13/2014 rhd             *
c     *                                                              *
c     *     creates a file name for patran output                    *
c     *     based on the output type and the step number.            *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ouflnm( strtnm, flnam, stepno, use_mpi, myid )
      implicit integer (a-z)
c      
      character(len=4) :: strtnm
      character(len=*) :: flnam
      logical use_mpi
c
      step_number = stepno 
c
c                check to make sure the step number is not
c                greater than 5 digits. 
c
      if( step_number .gt. 99999 ) step_number = step_number - 99999
c
      flnam(1:) = ' '             
      flnam(1:4) = strtnm  ! first part of name
      write(flnam(5:),9100) step_number
c
c                for mpi, we add the process id as a suffix
c                to the file name.
c
      if ( .not. use_mpi ) return
      flnam(10:10) = '.'
      write(flnam(11:14),9000) myid
      return
c
 9000 format(i4.4)
 9100 format(i5.5)
c
      end
