c     ****************************************************************
c     *                                                              *
c     *                subroutine outime                             *
c     *  output cpu times for each major part of solution at         *
c     *  job termination                                             *
c     *                                                              *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 01/20/2017 rhd             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine outime 
c     use global_data ! old common.main
      implicit none   
      include 'common.main'  
c
      integer :: calc   
      character(len=30) :: clctyp
      real :: t1
      real, external :: wcputime
c
      write(out,9000)
c
      t1 = wcputime( 1 )
c
      do calc = 1, mxtim
         if( times(calc,2) .lt. 0.01 ) cycle
         if(calc.eq.1) then
            clctyp= 'pcg solution vector update:   '
         else if(calc.eq.2) then
            clctyp= 'sig-eps & internal force:     '
         end if
         write(out,9010) clctyp
         write(out,9020) times(calc,1),
     &        100.0*times(calc,1)/t1, int(times(calc,2))
c
      end do
c
 9000 format(////1x,'>>>>>  solution timings   <<<<<')
c
 9010 format(//2x,'calculations for ',a30)
c
 9020 format(/3x,'wall time (secs): ',f10.4,
     &       f5.1,' (%) ', 'no. calls: ',i8)
c
      return
      end
