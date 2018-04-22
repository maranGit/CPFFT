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
      include 'common.main'

c                         global
      character(len=*) :: sparam
      integer :: errnum, param
      double precision :: dparam
      real :: rparam
c                        local

      select case (errnum)
      case (1)
      case (2)
      case (3)
      case default
      end select

      return
      end

