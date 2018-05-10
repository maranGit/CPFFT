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
      case (13)
        write(out,9013)
 9013   format(/,1x,'>>>>> Error: Invalid in STRAINS, not a number.',/)
      case (14)
        write(out,9014)
 9014   format(/,1x,'>>>>> Error: Invalid step list.',/)
      case (15)
        write(out,9015)
 9015   format(/,1x,'>>>>> Error: Invalid constraints. Not a number.',/)
      case (16)
        write(out,9016) param
 9016   format(/,1x,'>>>>> Error: Invalid step: ',i6,' in loading.',/)
      case (17)
        write(out,9017) sparam
 9017   format(/,1x,'>>>>> Error: ',a24,' in nonlinear analysis '
     &    'parameters must be a number',/)
      case (18)
        write(out,9018)
 9018   format(/,1x,'>>>>> Error: Unknown command in convergence
     &    tolerance.',/)
      case (19)
        write(out,9019)
 9019   format(/,1x,'>>>>> Warning: patran file name is not provided, ',
     &              'default name is used',/)
      case (20)
        write(out,9020)
 9020   format(/,1x,'>>>>> Error: Invalid file option. Should be ',
     &              'either flat file or patran file.',/)
      case (21)
        write(out,9021)
 9021   format(/,1x,'>>>>> Error: Invalid output option. Should be ',
     &              'd, v, a, r or T.',/)
      case default
        write(out,9999)
 9999   format(/,1x,'>>>>> Error: Unrecogonized option in errmsg.f',/)
      end select

      return
      end
