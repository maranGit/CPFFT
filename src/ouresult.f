c     ****************************************************************
c     *                                                              *
c     *                      subroutine ouresult                     *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                   last modified : 5/07/2018 RM               *
c     *                                                              *
c     *           drive nodal and elemental results output           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ouresult(step)
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
c
c                 compute nodal displacement from F
c                 store in common.main -> u
c
      call f2disp(step)
c
c                 open result file
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

      call ouocdd( dva, step, oubin, ouasc, bnfile, fmfile,
     &             opt, use_mpi, myid, flat_file, stream_file,
     &             text_file, compressed, flat_file_number )

c
c                 write result file
c
      outype      = 1
      quantity    = 1
      call ouddpa_flat_header(outype, quantity, flat_file_number)
      do nod = 1, nonode
        currdof = 3*nod - 2
        if( flat_file .and. text_file ) 
     &    write(flat_file_number,930) u((currdof-2):currdof)
 
      end do  ! on all nodes
c
c                 close patran or flat file. done.
c
      call ouocdd( dva, step, oubin, ouasc, bnfile, fmfile, 2,
     &             .false., myid, flat_file, stream_file, text_file, 
     &             compressed, flat_file_number )
c     

      return
 930  format(3e15.6)
      end subroutine

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
