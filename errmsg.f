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
        write(out,9001) 
 9001   format(/1x,'>>>>> warning: the name of the structure has not ',
     &             'been input. '/)
      case (2)
        write(out,9002)
 9002   format(/1x,'>>>>> Error: number of grids must be integer. '/)
      case (3)
        write(out,9003)
 9003   format(/1x,'>>>>> Error: number of material must be integer. '/)
      case (4)
        write(out,9004)
 9004   format(/1x,'>>>>> Error: unrecogonized term. '/)
      case (5)
        write(out,9005)
 9005   format(/,1x,'>>>>> Error: the name of the material is ',
     &           'expected. abort the input'/16x,'of this material ',
     &          'and scan for another high level',/,16x,'command.',/)
      case (6)
        write(out,9006)
 9006   format(/,1x,'>>>>> Error: maximum materials is 500.',/)
      case (7)
        write(out,9007) sparam
 9007   format(/,1x,'>>>>> Error: Invalid ',a5,' in material ',
     &          'parameters, not a number.',/)
      case (8)
        write(out,9008)
 9008   format(/,1x,'>>>>> Error: properties must have a type. ',
     &          'Returning to high level command not.',/)
      case (9)
        write(out,9009)
 9009   format(/,1x,'>>>>> Error: invalid element list',/)
      case (10)
        write(out,9010) sparam
 9010   format(/,1x,'>>>>> Error: unknown term ',a24,' in element',/)
      case (11)
        write(out,9011)
 9011   format(/,1x,'>>>>> Error: material must be followed by '
     &         'a defined type in element definition',/)
      case (12)
        write(out,9012) sparam
 9012   format(/,1x,'>>>>> Error: in element definition, '
     &         'material ',a8,' is not defined',/)
      case default
        write(out,9999)
 9999   format(/,1x,'>>>>> Error: Unrecogonized option in errmsg.f',/)
      end select

      return
      end

