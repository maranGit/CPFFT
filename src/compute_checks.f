      subroutine compute_checks
      use  fft, only : cp_matls_present, matprp, imatprp
      implicit none
      include 'common.main'
c
      integer :: matnum
      logical :: found_cp

      cp_matls_present = -1 ! check every time this subroutine is called
c
c
c              check for presence of material(s) using crystal
c              plasticity.
c
c              cp_matls_present = -1 we've not yet checked
c                               =  0 we checked and no CP materials
c                               =  1 we check. there are CP materials
c                                  and we need to call CP setup
c                                  routine if not already done.
c                               for MPI, workers need sizes from root
c
      if( cp_matls_present == -1 ) then
         found_cp = .false.
         do matnum = 1, nummat  ! nummat in common.main
           if(  matprp(9,matnum) /= 10 ) cycle
           found_cp = .true.
           exit
         end do
         cp_matls_present = 0
         if( found_cp ) then
            cp_matls_present = 1
            call mm10_set_history_locs
            call wmpi_compute_set_history_locs                                  
         end if                                                                 
      end if                                                                    
c                                                                               
c              add more setups here ... may need to add another                 
c              wmpi_... routine as well                                         
c                       
      return
      end subroutine
