c     ****************************************************************
c     *                      suboutine errmsg                        *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 6/20/2016 rhd              *
c     *                                                              *
c     *     this subroutine prints assorted error messages in re-    *
c     *     ponse to calls from all over the program. virtually all  *
c     *     error messages in the program are generated here.        *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine errmsg(errnum,param,sparam,rparam,dparam)
      implicit none
      character(len=*) :: sparam
c
      integer :: errnum, param
      double precision :: dparam
      real :: erprmr, rparam
c
      include 'common.main'

      return
      end

