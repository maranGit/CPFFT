c     ****************************************************************
c     *                                                              *
c     *                      subroutine cep2A                        *
c     *                                                              *
c     *                       written by : RM                        *
c     *                                                              *
c     *                 last modified : 3/20/2018 RM                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine cep2A(span, sigma, cep, F_inv, detF, A4, iout)
      implicit none
c     include 'param_def'
      integer, parameter :: mxvl = 128, nstr = 6, nstrs = 9
c
c                   global
      integer, intent(in) :: span, iout
      real(8), intent(in) :: sigma(mxvl,nstr), cep(mxvl,nstr,*)
      real(8), intent(in) :: detF(*), F_inv(mxvl,*)
      real(8), intent(out):: A4(mxvl,*)
c
c                   local
      integer :: ii
      logical :: local_debug
      real(8) :: matgeo(mxvl,nstrs*nstrs)
c     real(8), allocatable :: matgeo(:,:)
c
      local_debug = .false.
c     allocate(matgeo(mxvl,nstrs*nstrs))
c
c     add geometric stiffness to total stiffness
      do ii = 1, span
        matgeo(ii,1) = cep(ii,1,1) + sigma(ii,1) * detF(ii)
        matgeo(ii,2) = cep(ii,1,4)
        matgeo(ii,3) = cep(ii,1,6)
        matgeo(ii,4) = cep(ii,1,4) + sigma(ii,4) * detF(ii)
        matgeo(ii,5) = cep(ii,1,2)
        matgeo(ii,6) = cep(ii,1,5)
        matgeo(ii,7) = cep(ii,1,6) + sigma(ii,6) * detF(ii)
        matgeo(ii,8) = cep(ii,1,5)
        matgeo(ii,9) = cep(ii,1,3)
        matgeo(ii,10) = cep(ii,4,1)
      end do
      do ii = 1, span
        matgeo(ii,11) = cep(ii,4,4) + sigma(ii,1) * detF(ii)
        matgeo(ii,12) = cep(ii,4,6)
        matgeo(ii,13) = cep(ii,4,4)
        matgeo(ii,14) = cep(ii,4,2) + sigma(ii,4) * detF(ii)
        matgeo(ii,15) = cep(ii,4,5)
        matgeo(ii,16) = cep(ii,4,6)
        matgeo(ii,17) = cep(ii,4,5) + sigma(ii,6) * detF(ii)
        matgeo(ii,18) = cep(ii,4,3)
        matgeo(ii,19) = cep(ii,6,1)
        matgeo(ii,20) = cep(ii,6,4)
      end do
      do ii = 1, span
        matgeo(ii,21) = cep(ii,6,6) + sigma(ii,1) * detF(ii)
        matgeo(ii,22) = cep(ii,6,4)
        matgeo(ii,23) = cep(ii,6,2)
        matgeo(ii,24) = cep(ii,6,5) + sigma(ii,4) * detF(ii)
        matgeo(ii,25) = cep(ii,6,6)
        matgeo(ii,26) = cep(ii,6,5)
        matgeo(ii,27) = cep(ii,6,3) + sigma(ii,6) * detF(ii)
        matgeo(ii,28) = cep(ii,4,1) + sigma(ii,4) * detF(ii)
        matgeo(ii,29) = cep(ii,4,4)
        matgeo(ii,30) = cep(ii,4,6)
      end do
      do ii = 1, span
        matgeo(ii,31) = cep(ii,4,4) + sigma(ii,2) * detF(ii)
        matgeo(ii,32) = cep(ii,4,2)
        matgeo(ii,33) = cep(ii,4,5)
        matgeo(ii,34) = cep(ii,4,6) + sigma(ii,5) * detF(ii)
        matgeo(ii,35) = cep(ii,4,5)
        matgeo(ii,36) = cep(ii,4,3)
        matgeo(ii,37) = cep(ii,2,1)
        matgeo(ii,38) = cep(ii,2,4) + sigma(ii,4) * detF(ii)
        matgeo(ii,39) = cep(ii,2,6)
        matgeo(ii,40) = cep(ii,2,4)
      end do
      do ii = 1, span
        matgeo(ii,41) = cep(ii,2,2) + sigma(ii,2) * detF(ii)
        matgeo(ii,42) = cep(ii,2,5)
        matgeo(ii,43) = cep(ii,2,6)
        matgeo(ii,44) = cep(ii,2,5) + sigma(ii,5) * detF(ii)
        matgeo(ii,45) = cep(ii,2,3)
        matgeo(ii,46) = cep(ii,5,1)
        matgeo(ii,47) = cep(ii,5,4)
        matgeo(ii,48) = cep(ii,5,6) + sigma(ii,4) * detF(ii)
        matgeo(ii,49) = cep(ii,5,4)
        matgeo(ii,50) = cep(ii,5,2)
      end do
      do ii = 1, span
        matgeo(ii,51) = cep(ii,5,5) + sigma(ii,2) * detF(ii)
        matgeo(ii,52) = cep(ii,5,6)
        matgeo(ii,53) = cep(ii,5,5)
        matgeo(ii,54) = cep(ii,5,3) + sigma(ii,5) * detF(ii)
        matgeo(ii,55) = cep(ii,6,1) + sigma(ii,6) * detF(ii)
        matgeo(ii,56) = cep(ii,6,4)
        matgeo(ii,57) = cep(ii,6,6)
        matgeo(ii,58) = cep(ii,6,4) + sigma(ii,5) * detF(ii)
        matgeo(ii,59) = cep(ii,6,2)
        matgeo(ii,60) = cep(ii,6,5)
      end do
      do ii = 1, span
        matgeo(ii,61) = cep(ii,6,6) + sigma(ii,3) * detF(ii)
        matgeo(ii,62) = cep(ii,6,5)
        matgeo(ii,63) = cep(ii,6,3)
        matgeo(ii,64) = cep(ii,5,1)
        matgeo(ii,65) = cep(ii,5,4) + sigma(ii,6) * detF(ii)
        matgeo(ii,66) = cep(ii,5,6)
        matgeo(ii,67) = cep(ii,5,4)
        matgeo(ii,68) = cep(ii,5,2) + sigma(ii,5) * detF(ii)
        matgeo(ii,69) = cep(ii,5,5)
        matgeo(ii,70) = cep(ii,5,6)
      end do
      do ii = 1, span
        matgeo(ii,71) = cep(ii,5,5) + sigma(ii,3) * detF(ii)
        matgeo(ii,72) = cep(ii,5,3)
        matgeo(ii,73) = cep(ii,3,1)
        matgeo(ii,74) = cep(ii,3,4)
        matgeo(ii,75) = cep(ii,3,6) + sigma(ii,6) * detF(ii)
        matgeo(ii,76) = cep(ii,3,4)
        matgeo(ii,77) = cep(ii,3,2)
        matgeo(ii,78) = cep(ii,3,5) + sigma(ii,5) * detF(ii)
        matgeo(ii,79) = cep(ii,3,6)
        matgeo(ii,80) = cep(ii,3,5)
        matgeo(ii,81) = cep(ii,3,3) + sigma(ii,3) * detF(ii)
      end do
c
c     check material + geometry
      if (local_debug) then
        write(iout, 1000)
        do ii = 1, span
          write(iout, 1001) ii
          write(iout, 1002) sigma(ii,:)
        enddo
      endif
c     pull back to reference configuration
      do ii = 1, span
        A4(ii,1) = F_inv(ii,1)*matgeo(ii,1)*F_inv(ii,1)
     &           + F_inv(ii,1)*matgeo(ii,2)*F_inv(ii,2)
     &           + F_inv(ii,1)*matgeo(ii,3)*F_inv(ii,3)
     &           + F_inv(ii,2)*matgeo(ii,10)*F_inv(ii,1)
     &           + F_inv(ii,2)*matgeo(ii,11)*F_inv(ii,2)
     &           + F_inv(ii,2)*matgeo(ii,12)*F_inv(ii,3)
     &           + F_inv(ii,3)*matgeo(ii,19)*F_inv(ii,1)
     &           + F_inv(ii,3)*matgeo(ii,20)*F_inv(ii,2)
     &           + F_inv(ii,3)*matgeo(ii,21)*F_inv(ii,3)
      end do
      do ii = 1, span
        A4(ii,2) = F_inv(ii,1)*matgeo(ii,1)*F_inv(ii,4)
     &           + F_inv(ii,1)*matgeo(ii,2)*F_inv(ii,5)
     &           + F_inv(ii,1)*matgeo(ii,3)*F_inv(ii,6)
     &           + F_inv(ii,2)*matgeo(ii,10)*F_inv(ii,4)
     &           + F_inv(ii,2)*matgeo(ii,11)*F_inv(ii,5)
     &           + F_inv(ii,2)*matgeo(ii,12)*F_inv(ii,6)
     &           + F_inv(ii,3)*matgeo(ii,19)*F_inv(ii,4)
     &           + F_inv(ii,3)*matgeo(ii,20)*F_inv(ii,5)
     &           + F_inv(ii,3)*matgeo(ii,21)*F_inv(ii,6)
      end do
      do ii = 1, span
        A4(ii,3) = F_inv(ii,1)*matgeo(ii,1)*F_inv(ii,7)
     &           + F_inv(ii,1)*matgeo(ii,2)*F_inv(ii,8)
     &           + F_inv(ii,1)*matgeo(ii,3)*F_inv(ii,9)
     &           + F_inv(ii,2)*matgeo(ii,10)*F_inv(ii,7)
     &           + F_inv(ii,2)*matgeo(ii,11)*F_inv(ii,8)
     &           + F_inv(ii,2)*matgeo(ii,12)*F_inv(ii,9)
     &           + F_inv(ii,3)*matgeo(ii,19)*F_inv(ii,7)
     &           + F_inv(ii,3)*matgeo(ii,20)*F_inv(ii,8)
     &           + F_inv(ii,3)*matgeo(ii,21)*F_inv(ii,9)
      end do
      do ii = 1, span
        A4(ii,4) = F_inv(ii,1)*matgeo(ii,4)*F_inv(ii,1)
     &           + F_inv(ii,1)*matgeo(ii,5)*F_inv(ii,2)
     &           + F_inv(ii,1)*matgeo(ii,6)*F_inv(ii,3)
     &           + F_inv(ii,2)*matgeo(ii,13)*F_inv(ii,1)
     &           + F_inv(ii,2)*matgeo(ii,14)*F_inv(ii,2)
     &           + F_inv(ii,2)*matgeo(ii,15)*F_inv(ii,3)
     &           + F_inv(ii,3)*matgeo(ii,22)*F_inv(ii,1)
     &           + F_inv(ii,3)*matgeo(ii,23)*F_inv(ii,2)
     &           + F_inv(ii,3)*matgeo(ii,24)*F_inv(ii,3)
      end do
      do ii = 1, span
        A4(ii,5) = F_inv(ii,1)*matgeo(ii,4)*F_inv(ii,4)
     &           + F_inv(ii,1)*matgeo(ii,5)*F_inv(ii,5)
     &           + F_inv(ii,1)*matgeo(ii,6)*F_inv(ii,6)
     &           + F_inv(ii,2)*matgeo(ii,13)*F_inv(ii,4)
     &           + F_inv(ii,2)*matgeo(ii,14)*F_inv(ii,5)
     &           + F_inv(ii,2)*matgeo(ii,15)*F_inv(ii,6)
     &           + F_inv(ii,3)*matgeo(ii,22)*F_inv(ii,4)
     &           + F_inv(ii,3)*matgeo(ii,23)*F_inv(ii,5)
     &           + F_inv(ii,3)*matgeo(ii,24)*F_inv(ii,6)
      end do
      do ii = 1, span
        A4(ii,6) = F_inv(ii,1)*matgeo(ii,4)*F_inv(ii,7)
     &           + F_inv(ii,1)*matgeo(ii,5)*F_inv(ii,8)
     &           + F_inv(ii,1)*matgeo(ii,6)*F_inv(ii,9)
     &           + F_inv(ii,2)*matgeo(ii,13)*F_inv(ii,7)
     &           + F_inv(ii,2)*matgeo(ii,14)*F_inv(ii,8)
     &           + F_inv(ii,2)*matgeo(ii,15)*F_inv(ii,9)
     &           + F_inv(ii,3)*matgeo(ii,22)*F_inv(ii,7)
     &           + F_inv(ii,3)*matgeo(ii,23)*F_inv(ii,8)
     &           + F_inv(ii,3)*matgeo(ii,24)*F_inv(ii,9)
      end do
      do ii = 1, span
        A4(ii,7) = F_inv(ii,1)*matgeo(ii,7)*F_inv(ii,1)
     &           + F_inv(ii,1)*matgeo(ii,8)*F_inv(ii,2)
     &           + F_inv(ii,1)*matgeo(ii,9)*F_inv(ii,3)
     &           + F_inv(ii,2)*matgeo(ii,16)*F_inv(ii,1)
     &           + F_inv(ii,2)*matgeo(ii,17)*F_inv(ii,2)
     &           + F_inv(ii,2)*matgeo(ii,18)*F_inv(ii,3)
     &           + F_inv(ii,3)*matgeo(ii,25)*F_inv(ii,1)
     &           + F_inv(ii,3)*matgeo(ii,26)*F_inv(ii,2)
     &           + F_inv(ii,3)*matgeo(ii,27)*F_inv(ii,3)
      end do
      do ii = 1, span
        A4(ii,8) = F_inv(ii,1)*matgeo(ii,7)*F_inv(ii,4)
     &           + F_inv(ii,1)*matgeo(ii,8)*F_inv(ii,5)
     &           + F_inv(ii,1)*matgeo(ii,9)*F_inv(ii,6)
     &           + F_inv(ii,2)*matgeo(ii,16)*F_inv(ii,4)
     &           + F_inv(ii,2)*matgeo(ii,17)*F_inv(ii,5)
     &           + F_inv(ii,2)*matgeo(ii,18)*F_inv(ii,6)
     &           + F_inv(ii,3)*matgeo(ii,25)*F_inv(ii,4)
     &           + F_inv(ii,3)*matgeo(ii,26)*F_inv(ii,5)
     &           + F_inv(ii,3)*matgeo(ii,27)*F_inv(ii,6)
      end do
      do ii = 1, span
        A4(ii,9) = F_inv(ii,1)*matgeo(ii,7)*F_inv(ii,7)
     &           + F_inv(ii,1)*matgeo(ii,8)*F_inv(ii,8)
     &           + F_inv(ii,1)*matgeo(ii,9)*F_inv(ii,9)
     &           + F_inv(ii,2)*matgeo(ii,16)*F_inv(ii,7)
     &           + F_inv(ii,2)*matgeo(ii,17)*F_inv(ii,8)
     &           + F_inv(ii,2)*matgeo(ii,18)*F_inv(ii,9)
     &           + F_inv(ii,3)*matgeo(ii,25)*F_inv(ii,7)
     &           + F_inv(ii,3)*matgeo(ii,26)*F_inv(ii,8)
     &           + F_inv(ii,3)*matgeo(ii,27)*F_inv(ii,9)
      end do
      do ii = 1, span
        A4(ii,10) = F_inv(ii,4)*matgeo(ii,1)*F_inv(ii,1)
     &           + F_inv(ii,4)*matgeo(ii,2)*F_inv(ii,2)
     &           + F_inv(ii,4)*matgeo(ii,3)*F_inv(ii,3)
     &           + F_inv(ii,5)*matgeo(ii,10)*F_inv(ii,1)
     &           + F_inv(ii,5)*matgeo(ii,11)*F_inv(ii,2)
     &           + F_inv(ii,5)*matgeo(ii,12)*F_inv(ii,3)
     &           + F_inv(ii,6)*matgeo(ii,19)*F_inv(ii,1)
     &           + F_inv(ii,6)*matgeo(ii,20)*F_inv(ii,2)
     &           + F_inv(ii,6)*matgeo(ii,21)*F_inv(ii,3)
      end do
      do ii = 1, span
        A4(ii,11) = F_inv(ii,4)*matgeo(ii,1)*F_inv(ii,4)
     &           + F_inv(ii,4)*matgeo(ii,2)*F_inv(ii,5)
     &           + F_inv(ii,4)*matgeo(ii,3)*F_inv(ii,6)
     &           + F_inv(ii,5)*matgeo(ii,10)*F_inv(ii,4)
     &           + F_inv(ii,5)*matgeo(ii,11)*F_inv(ii,5)
     &           + F_inv(ii,5)*matgeo(ii,12)*F_inv(ii,6)
     &           + F_inv(ii,6)*matgeo(ii,19)*F_inv(ii,4)
     &           + F_inv(ii,6)*matgeo(ii,20)*F_inv(ii,5)
     &           + F_inv(ii,6)*matgeo(ii,21)*F_inv(ii,6)
      end do
      do ii = 1, span
        A4(ii,12) = F_inv(ii,4)*matgeo(ii,1)*F_inv(ii,7)
     &           + F_inv(ii,4)*matgeo(ii,2)*F_inv(ii,8)
     &           + F_inv(ii,4)*matgeo(ii,3)*F_inv(ii,9)
     &           + F_inv(ii,5)*matgeo(ii,10)*F_inv(ii,7)
     &           + F_inv(ii,5)*matgeo(ii,11)*F_inv(ii,8)
     &           + F_inv(ii,5)*matgeo(ii,12)*F_inv(ii,9)
     &           + F_inv(ii,6)*matgeo(ii,19)*F_inv(ii,7)
     &           + F_inv(ii,6)*matgeo(ii,20)*F_inv(ii,8)
     &           + F_inv(ii,6)*matgeo(ii,21)*F_inv(ii,9)
      end do
      do ii = 1, span
        A4(ii,13) = F_inv(ii,4)*matgeo(ii,4)*F_inv(ii,1)
     &           + F_inv(ii,4)*matgeo(ii,5)*F_inv(ii,2)
     &           + F_inv(ii,4)*matgeo(ii,6)*F_inv(ii,3)
     &           + F_inv(ii,5)*matgeo(ii,13)*F_inv(ii,1)
     &           + F_inv(ii,5)*matgeo(ii,14)*F_inv(ii,2)
     &           + F_inv(ii,5)*matgeo(ii,15)*F_inv(ii,3)
     &           + F_inv(ii,6)*matgeo(ii,22)*F_inv(ii,1)
     &           + F_inv(ii,6)*matgeo(ii,23)*F_inv(ii,2)
     &           + F_inv(ii,6)*matgeo(ii,24)*F_inv(ii,3)
      end do
      do ii = 1, span
        A4(ii,14) = F_inv(ii,4)*matgeo(ii,4)*F_inv(ii,4)
     &           + F_inv(ii,4)*matgeo(ii,5)*F_inv(ii,5)
     &           + F_inv(ii,4)*matgeo(ii,6)*F_inv(ii,6)
     &           + F_inv(ii,5)*matgeo(ii,13)*F_inv(ii,4)
     &           + F_inv(ii,5)*matgeo(ii,14)*F_inv(ii,5)
     &           + F_inv(ii,5)*matgeo(ii,15)*F_inv(ii,6)
     &           + F_inv(ii,6)*matgeo(ii,22)*F_inv(ii,4)
     &           + F_inv(ii,6)*matgeo(ii,23)*F_inv(ii,5)
     &           + F_inv(ii,6)*matgeo(ii,24)*F_inv(ii,6)
      end do
      do ii = 1, span
        A4(ii,15) = F_inv(ii,4)*matgeo(ii,4)*F_inv(ii,7)
     &           + F_inv(ii,4)*matgeo(ii,5)*F_inv(ii,8)
     &           + F_inv(ii,4)*matgeo(ii,6)*F_inv(ii,9)
     &           + F_inv(ii,5)*matgeo(ii,13)*F_inv(ii,7)
     &           + F_inv(ii,5)*matgeo(ii,14)*F_inv(ii,8)
     &           + F_inv(ii,5)*matgeo(ii,15)*F_inv(ii,9)
     &           + F_inv(ii,6)*matgeo(ii,22)*F_inv(ii,7)
     &           + F_inv(ii,6)*matgeo(ii,23)*F_inv(ii,8)
     &           + F_inv(ii,6)*matgeo(ii,24)*F_inv(ii,9)
      end do
      do ii = 1, span
        A4(ii,16) = F_inv(ii,4)*matgeo(ii,7)*F_inv(ii,1)
     &           + F_inv(ii,4)*matgeo(ii,8)*F_inv(ii,2)
     &           + F_inv(ii,4)*matgeo(ii,9)*F_inv(ii,3)
     &           + F_inv(ii,5)*matgeo(ii,16)*F_inv(ii,1)
     &           + F_inv(ii,5)*matgeo(ii,17)*F_inv(ii,2)
     &           + F_inv(ii,5)*matgeo(ii,18)*F_inv(ii,3)
     &           + F_inv(ii,6)*matgeo(ii,25)*F_inv(ii,1)
     &           + F_inv(ii,6)*matgeo(ii,26)*F_inv(ii,2)
     &           + F_inv(ii,6)*matgeo(ii,27)*F_inv(ii,3)
      end do
      do ii = 1, span
        A4(ii,17) = F_inv(ii,4)*matgeo(ii,7)*F_inv(ii,4)
     &           + F_inv(ii,4)*matgeo(ii,8)*F_inv(ii,5)
     &           + F_inv(ii,4)*matgeo(ii,9)*F_inv(ii,6)
     &           + F_inv(ii,5)*matgeo(ii,16)*F_inv(ii,4)
     &           + F_inv(ii,5)*matgeo(ii,17)*F_inv(ii,5)
     &           + F_inv(ii,5)*matgeo(ii,18)*F_inv(ii,6)
     &           + F_inv(ii,6)*matgeo(ii,25)*F_inv(ii,4)
     &           + F_inv(ii,6)*matgeo(ii,26)*F_inv(ii,5)
     &           + F_inv(ii,6)*matgeo(ii,27)*F_inv(ii,6)
      end do
      do ii = 1, span
        A4(ii,18) = F_inv(ii,4)*matgeo(ii,7)*F_inv(ii,7)
     &           + F_inv(ii,4)*matgeo(ii,8)*F_inv(ii,8)
     &           + F_inv(ii,4)*matgeo(ii,9)*F_inv(ii,9)
     &           + F_inv(ii,5)*matgeo(ii,16)*F_inv(ii,7)
     &           + F_inv(ii,5)*matgeo(ii,17)*F_inv(ii,8)
     &           + F_inv(ii,5)*matgeo(ii,18)*F_inv(ii,9)
     &           + F_inv(ii,6)*matgeo(ii,25)*F_inv(ii,7)
     &           + F_inv(ii,6)*matgeo(ii,26)*F_inv(ii,8)
     &           + F_inv(ii,6)*matgeo(ii,27)*F_inv(ii,9)
      end do
      do ii = 1, span
        A4(ii,19) = F_inv(ii,7)*matgeo(ii,1)*F_inv(ii,1)
     &           + F_inv(ii,7)*matgeo(ii,2)*F_inv(ii,2)
     &           + F_inv(ii,7)*matgeo(ii,3)*F_inv(ii,3)
     &           + F_inv(ii,8)*matgeo(ii,10)*F_inv(ii,1)
     &           + F_inv(ii,8)*matgeo(ii,11)*F_inv(ii,2)
     &           + F_inv(ii,8)*matgeo(ii,12)*F_inv(ii,3)
     &           + F_inv(ii,9)*matgeo(ii,19)*F_inv(ii,1)
     &           + F_inv(ii,9)*matgeo(ii,20)*F_inv(ii,2)
     &           + F_inv(ii,9)*matgeo(ii,21)*F_inv(ii,3)
      end do
      do ii = 1, span
        A4(ii,20) = F_inv(ii,7)*matgeo(ii,1)*F_inv(ii,4)
     &           + F_inv(ii,7)*matgeo(ii,2)*F_inv(ii,5)
     &           + F_inv(ii,7)*matgeo(ii,3)*F_inv(ii,6)
     &           + F_inv(ii,8)*matgeo(ii,10)*F_inv(ii,4)
     &           + F_inv(ii,8)*matgeo(ii,11)*F_inv(ii,5)
     &           + F_inv(ii,8)*matgeo(ii,12)*F_inv(ii,6)
     &           + F_inv(ii,9)*matgeo(ii,19)*F_inv(ii,4)
     &           + F_inv(ii,9)*matgeo(ii,20)*F_inv(ii,5)
     &           + F_inv(ii,9)*matgeo(ii,21)*F_inv(ii,6)
      end do
      do ii = 1, span
        A4(ii,21) = F_inv(ii,7)*matgeo(ii,1)*F_inv(ii,7)
     &           + F_inv(ii,7)*matgeo(ii,2)*F_inv(ii,8)
     &           + F_inv(ii,7)*matgeo(ii,3)*F_inv(ii,9)
     &           + F_inv(ii,8)*matgeo(ii,10)*F_inv(ii,7)
     &           + F_inv(ii,8)*matgeo(ii,11)*F_inv(ii,8)
     &           + F_inv(ii,8)*matgeo(ii,12)*F_inv(ii,9)
     &           + F_inv(ii,9)*matgeo(ii,19)*F_inv(ii,7)
     &           + F_inv(ii,9)*matgeo(ii,20)*F_inv(ii,8)
     &           + F_inv(ii,9)*matgeo(ii,21)*F_inv(ii,9)
      end do
      do ii = 1, span
        A4(ii,22) = F_inv(ii,7)*matgeo(ii,4)*F_inv(ii,1)
     &           + F_inv(ii,7)*matgeo(ii,5)*F_inv(ii,2)
     &           + F_inv(ii,7)*matgeo(ii,6)*F_inv(ii,3)
     &           + F_inv(ii,8)*matgeo(ii,13)*F_inv(ii,1)
     &           + F_inv(ii,8)*matgeo(ii,14)*F_inv(ii,2)
     &           + F_inv(ii,8)*matgeo(ii,15)*F_inv(ii,3)
     &           + F_inv(ii,9)*matgeo(ii,22)*F_inv(ii,1)
     &           + F_inv(ii,9)*matgeo(ii,23)*F_inv(ii,2)
     &           + F_inv(ii,9)*matgeo(ii,24)*F_inv(ii,3)
      end do
      do ii = 1, span
        A4(ii,23) = F_inv(ii,7)*matgeo(ii,4)*F_inv(ii,4)
     &           + F_inv(ii,7)*matgeo(ii,5)*F_inv(ii,5)
     &           + F_inv(ii,7)*matgeo(ii,6)*F_inv(ii,6)
     &           + F_inv(ii,8)*matgeo(ii,13)*F_inv(ii,4)
     &           + F_inv(ii,8)*matgeo(ii,14)*F_inv(ii,5)
     &           + F_inv(ii,8)*matgeo(ii,15)*F_inv(ii,6)
     &           + F_inv(ii,9)*matgeo(ii,22)*F_inv(ii,4)
     &           + F_inv(ii,9)*matgeo(ii,23)*F_inv(ii,5)
     &           + F_inv(ii,9)*matgeo(ii,24)*F_inv(ii,6)
      end do
      do ii = 1, span
        A4(ii,24) = F_inv(ii,7)*matgeo(ii,4)*F_inv(ii,7)
     &           + F_inv(ii,7)*matgeo(ii,5)*F_inv(ii,8)
     &           + F_inv(ii,7)*matgeo(ii,6)*F_inv(ii,9)
     &           + F_inv(ii,8)*matgeo(ii,13)*F_inv(ii,7)
     &           + F_inv(ii,8)*matgeo(ii,14)*F_inv(ii,8)
     &           + F_inv(ii,8)*matgeo(ii,15)*F_inv(ii,9)
     &           + F_inv(ii,9)*matgeo(ii,22)*F_inv(ii,7)
     &           + F_inv(ii,9)*matgeo(ii,23)*F_inv(ii,8)
     &           + F_inv(ii,9)*matgeo(ii,24)*F_inv(ii,9)
      end do
      do ii = 1, span
        A4(ii,25) = F_inv(ii,7)*matgeo(ii,7)*F_inv(ii,1)
     &           + F_inv(ii,7)*matgeo(ii,8)*F_inv(ii,2)
     &           + F_inv(ii,7)*matgeo(ii,9)*F_inv(ii,3)
     &           + F_inv(ii,8)*matgeo(ii,16)*F_inv(ii,1)
     &           + F_inv(ii,8)*matgeo(ii,17)*F_inv(ii,2)
     &           + F_inv(ii,8)*matgeo(ii,18)*F_inv(ii,3)
     &           + F_inv(ii,9)*matgeo(ii,25)*F_inv(ii,1)
     &           + F_inv(ii,9)*matgeo(ii,26)*F_inv(ii,2)
     &           + F_inv(ii,9)*matgeo(ii,27)*F_inv(ii,3)
      end do
      do ii = 1, span
        A4(ii,26) = F_inv(ii,7)*matgeo(ii,7)*F_inv(ii,4)
     &           + F_inv(ii,7)*matgeo(ii,8)*F_inv(ii,5)
     &           + F_inv(ii,7)*matgeo(ii,9)*F_inv(ii,6)
     &           + F_inv(ii,8)*matgeo(ii,16)*F_inv(ii,4)
     &           + F_inv(ii,8)*matgeo(ii,17)*F_inv(ii,5)
     &           + F_inv(ii,8)*matgeo(ii,18)*F_inv(ii,6)
     &           + F_inv(ii,9)*matgeo(ii,25)*F_inv(ii,4)
     &           + F_inv(ii,9)*matgeo(ii,26)*F_inv(ii,5)
     &           + F_inv(ii,9)*matgeo(ii,27)*F_inv(ii,6)
      end do
      do ii = 1, span
        A4(ii,27) = F_inv(ii,7)*matgeo(ii,7)*F_inv(ii,7)
     &           + F_inv(ii,7)*matgeo(ii,8)*F_inv(ii,8)
     &           + F_inv(ii,7)*matgeo(ii,9)*F_inv(ii,9)
     &           + F_inv(ii,8)*matgeo(ii,16)*F_inv(ii,7)
     &           + F_inv(ii,8)*matgeo(ii,17)*F_inv(ii,8)
     &           + F_inv(ii,8)*matgeo(ii,18)*F_inv(ii,9)
     &           + F_inv(ii,9)*matgeo(ii,25)*F_inv(ii,7)
     &           + F_inv(ii,9)*matgeo(ii,26)*F_inv(ii,8)
     &           + F_inv(ii,9)*matgeo(ii,27)*F_inv(ii,9)
      end do
      do ii = 1, span
        A4(ii,28) = F_inv(ii,1)*matgeo(ii,28)*F_inv(ii,1)
     &           + F_inv(ii,1)*matgeo(ii,29)*F_inv(ii,2)
     &           + F_inv(ii,1)*matgeo(ii,30)*F_inv(ii,3)
     &           + F_inv(ii,2)*matgeo(ii,37)*F_inv(ii,1)
     &           + F_inv(ii,2)*matgeo(ii,38)*F_inv(ii,2)
     &           + F_inv(ii,2)*matgeo(ii,39)*F_inv(ii,3)
     &           + F_inv(ii,3)*matgeo(ii,46)*F_inv(ii,1)
     &           + F_inv(ii,3)*matgeo(ii,47)*F_inv(ii,2)
     &           + F_inv(ii,3)*matgeo(ii,48)*F_inv(ii,3)
      end do
      do ii = 1, span
        A4(ii,29) = F_inv(ii,1)*matgeo(ii,28)*F_inv(ii,4)
     &           + F_inv(ii,1)*matgeo(ii,29)*F_inv(ii,5)
     &           + F_inv(ii,1)*matgeo(ii,30)*F_inv(ii,6)
     &           + F_inv(ii,2)*matgeo(ii,37)*F_inv(ii,4)
     &           + F_inv(ii,2)*matgeo(ii,38)*F_inv(ii,5)
     &           + F_inv(ii,2)*matgeo(ii,39)*F_inv(ii,6)
     &           + F_inv(ii,3)*matgeo(ii,46)*F_inv(ii,4)
     &           + F_inv(ii,3)*matgeo(ii,47)*F_inv(ii,5)
     &           + F_inv(ii,3)*matgeo(ii,48)*F_inv(ii,6)
      end do
      do ii = 1, span
        A4(ii,30) = F_inv(ii,1)*matgeo(ii,28)*F_inv(ii,7)
     &           + F_inv(ii,1)*matgeo(ii,29)*F_inv(ii,8)
     &           + F_inv(ii,1)*matgeo(ii,30)*F_inv(ii,9)
     &           + F_inv(ii,2)*matgeo(ii,37)*F_inv(ii,7)
     &           + F_inv(ii,2)*matgeo(ii,38)*F_inv(ii,8)
     &           + F_inv(ii,2)*matgeo(ii,39)*F_inv(ii,9)
     &           + F_inv(ii,3)*matgeo(ii,46)*F_inv(ii,7)
     &           + F_inv(ii,3)*matgeo(ii,47)*F_inv(ii,8)
     &           + F_inv(ii,3)*matgeo(ii,48)*F_inv(ii,9)
      end do
      do ii = 1, span
        A4(ii,31) = F_inv(ii,1)*matgeo(ii,31)*F_inv(ii,1)
     &           + F_inv(ii,1)*matgeo(ii,32)*F_inv(ii,2)
     &           + F_inv(ii,1)*matgeo(ii,33)*F_inv(ii,3)
     &           + F_inv(ii,2)*matgeo(ii,40)*F_inv(ii,1)
     &           + F_inv(ii,2)*matgeo(ii,41)*F_inv(ii,2)
     &           + F_inv(ii,2)*matgeo(ii,42)*F_inv(ii,3)
     &           + F_inv(ii,3)*matgeo(ii,49)*F_inv(ii,1)
     &           + F_inv(ii,3)*matgeo(ii,50)*F_inv(ii,2)
     &           + F_inv(ii,3)*matgeo(ii,51)*F_inv(ii,3)
      end do
      do ii = 1, span
        A4(ii,32) = F_inv(ii,1)*matgeo(ii,31)*F_inv(ii,4)
     &           + F_inv(ii,1)*matgeo(ii,32)*F_inv(ii,5)
     &           + F_inv(ii,1)*matgeo(ii,33)*F_inv(ii,6)
     &           + F_inv(ii,2)*matgeo(ii,40)*F_inv(ii,4)
     &           + F_inv(ii,2)*matgeo(ii,41)*F_inv(ii,5)
     &           + F_inv(ii,2)*matgeo(ii,42)*F_inv(ii,6)
     &           + F_inv(ii,3)*matgeo(ii,49)*F_inv(ii,4)
     &           + F_inv(ii,3)*matgeo(ii,50)*F_inv(ii,5)
     &           + F_inv(ii,3)*matgeo(ii,51)*F_inv(ii,6)
      end do
      do ii = 1, span
        A4(ii,33) = F_inv(ii,1)*matgeo(ii,31)*F_inv(ii,7)
     &           + F_inv(ii,1)*matgeo(ii,32)*F_inv(ii,8)
     &           + F_inv(ii,1)*matgeo(ii,33)*F_inv(ii,9)
     &           + F_inv(ii,2)*matgeo(ii,40)*F_inv(ii,7)
     &           + F_inv(ii,2)*matgeo(ii,41)*F_inv(ii,8)
     &           + F_inv(ii,2)*matgeo(ii,42)*F_inv(ii,9)
     &           + F_inv(ii,3)*matgeo(ii,49)*F_inv(ii,7)
     &           + F_inv(ii,3)*matgeo(ii,50)*F_inv(ii,8)
     &           + F_inv(ii,3)*matgeo(ii,51)*F_inv(ii,9)
      end do
      do ii = 1, span
        A4(ii,34) = F_inv(ii,1)*matgeo(ii,34)*F_inv(ii,1)
     &           + F_inv(ii,1)*matgeo(ii,35)*F_inv(ii,2)
     &           + F_inv(ii,1)*matgeo(ii,36)*F_inv(ii,3)
     &           + F_inv(ii,2)*matgeo(ii,43)*F_inv(ii,1)
     &           + F_inv(ii,2)*matgeo(ii,44)*F_inv(ii,2)
     &           + F_inv(ii,2)*matgeo(ii,45)*F_inv(ii,3)
     &           + F_inv(ii,3)*matgeo(ii,52)*F_inv(ii,1)
     &           + F_inv(ii,3)*matgeo(ii,53)*F_inv(ii,2)
     &           + F_inv(ii,3)*matgeo(ii,54)*F_inv(ii,3)
      end do
      do ii = 1, span
        A4(ii,35) = F_inv(ii,1)*matgeo(ii,34)*F_inv(ii,4)
     &           + F_inv(ii,1)*matgeo(ii,35)*F_inv(ii,5)
     &           + F_inv(ii,1)*matgeo(ii,36)*F_inv(ii,6)
     &           + F_inv(ii,2)*matgeo(ii,43)*F_inv(ii,4)
     &           + F_inv(ii,2)*matgeo(ii,44)*F_inv(ii,5)
     &           + F_inv(ii,2)*matgeo(ii,45)*F_inv(ii,6)
     &           + F_inv(ii,3)*matgeo(ii,52)*F_inv(ii,4)
     &           + F_inv(ii,3)*matgeo(ii,53)*F_inv(ii,5)
     &           + F_inv(ii,3)*matgeo(ii,54)*F_inv(ii,6)
      end do
      do ii = 1, span
        A4(ii,36) = F_inv(ii,1)*matgeo(ii,34)*F_inv(ii,7)
     &           + F_inv(ii,1)*matgeo(ii,35)*F_inv(ii,8)
     &           + F_inv(ii,1)*matgeo(ii,36)*F_inv(ii,9)
     &           + F_inv(ii,2)*matgeo(ii,43)*F_inv(ii,7)
     &           + F_inv(ii,2)*matgeo(ii,44)*F_inv(ii,8)
     &           + F_inv(ii,2)*matgeo(ii,45)*F_inv(ii,9)
     &           + F_inv(ii,3)*matgeo(ii,52)*F_inv(ii,7)
     &           + F_inv(ii,3)*matgeo(ii,53)*F_inv(ii,8)
     &           + F_inv(ii,3)*matgeo(ii,54)*F_inv(ii,9)
      end do
      do ii = 1, span
        A4(ii,37) = F_inv(ii,4)*matgeo(ii,28)*F_inv(ii,1)
     &           + F_inv(ii,4)*matgeo(ii,29)*F_inv(ii,2)
     &           + F_inv(ii,4)*matgeo(ii,30)*F_inv(ii,3)
     &           + F_inv(ii,5)*matgeo(ii,37)*F_inv(ii,1)
     &           + F_inv(ii,5)*matgeo(ii,38)*F_inv(ii,2)
     &           + F_inv(ii,5)*matgeo(ii,39)*F_inv(ii,3)
     &           + F_inv(ii,6)*matgeo(ii,46)*F_inv(ii,1)
     &           + F_inv(ii,6)*matgeo(ii,47)*F_inv(ii,2)
     &           + F_inv(ii,6)*matgeo(ii,48)*F_inv(ii,3)
      end do
      do ii = 1, span
        A4(ii,38) = F_inv(ii,4)*matgeo(ii,28)*F_inv(ii,4)
     &           + F_inv(ii,4)*matgeo(ii,29)*F_inv(ii,5)
     &           + F_inv(ii,4)*matgeo(ii,30)*F_inv(ii,6)
     &           + F_inv(ii,5)*matgeo(ii,37)*F_inv(ii,4)
     &           + F_inv(ii,5)*matgeo(ii,38)*F_inv(ii,5)
     &           + F_inv(ii,5)*matgeo(ii,39)*F_inv(ii,6)
     &           + F_inv(ii,6)*matgeo(ii,46)*F_inv(ii,4)
     &           + F_inv(ii,6)*matgeo(ii,47)*F_inv(ii,5)
     &           + F_inv(ii,6)*matgeo(ii,48)*F_inv(ii,6)
      end do
      do ii = 1, span
        A4(ii,39) = F_inv(ii,4)*matgeo(ii,28)*F_inv(ii,7)
     &           + F_inv(ii,4)*matgeo(ii,29)*F_inv(ii,8)
     &           + F_inv(ii,4)*matgeo(ii,30)*F_inv(ii,9)
     &           + F_inv(ii,5)*matgeo(ii,37)*F_inv(ii,7)
     &           + F_inv(ii,5)*matgeo(ii,38)*F_inv(ii,8)
     &           + F_inv(ii,5)*matgeo(ii,39)*F_inv(ii,9)
     &           + F_inv(ii,6)*matgeo(ii,46)*F_inv(ii,7)
     &           + F_inv(ii,6)*matgeo(ii,47)*F_inv(ii,8)
     &           + F_inv(ii,6)*matgeo(ii,48)*F_inv(ii,9)
      end do
      do ii = 1, span
        A4(ii,40) = F_inv(ii,4)*matgeo(ii,31)*F_inv(ii,1)
     &           + F_inv(ii,4)*matgeo(ii,32)*F_inv(ii,2)
     &           + F_inv(ii,4)*matgeo(ii,33)*F_inv(ii,3)
     &           + F_inv(ii,5)*matgeo(ii,40)*F_inv(ii,1)
     &           + F_inv(ii,5)*matgeo(ii,41)*F_inv(ii,2)
     &           + F_inv(ii,5)*matgeo(ii,42)*F_inv(ii,3)
     &           + F_inv(ii,6)*matgeo(ii,49)*F_inv(ii,1)
     &           + F_inv(ii,6)*matgeo(ii,50)*F_inv(ii,2)
     &           + F_inv(ii,6)*matgeo(ii,51)*F_inv(ii,3)
      end do
      do ii = 1, span
        A4(ii,41) = F_inv(ii,4)*matgeo(ii,31)*F_inv(ii,4)
     &           + F_inv(ii,4)*matgeo(ii,32)*F_inv(ii,5)
     &           + F_inv(ii,4)*matgeo(ii,33)*F_inv(ii,6)
     &           + F_inv(ii,5)*matgeo(ii,40)*F_inv(ii,4)
     &           + F_inv(ii,5)*matgeo(ii,41)*F_inv(ii,5)
     &           + F_inv(ii,5)*matgeo(ii,42)*F_inv(ii,6)
     &           + F_inv(ii,6)*matgeo(ii,49)*F_inv(ii,4)
     &           + F_inv(ii,6)*matgeo(ii,50)*F_inv(ii,5)
     &           + F_inv(ii,6)*matgeo(ii,51)*F_inv(ii,6)
      end do
      do ii = 1, span
        A4(ii,42) = F_inv(ii,4)*matgeo(ii,31)*F_inv(ii,7)
     &           + F_inv(ii,4)*matgeo(ii,32)*F_inv(ii,8)
     &           + F_inv(ii,4)*matgeo(ii,33)*F_inv(ii,9)
     &           + F_inv(ii,5)*matgeo(ii,40)*F_inv(ii,7)
     &           + F_inv(ii,5)*matgeo(ii,41)*F_inv(ii,8)
     &           + F_inv(ii,5)*matgeo(ii,42)*F_inv(ii,9)
     &           + F_inv(ii,6)*matgeo(ii,49)*F_inv(ii,7)
     &           + F_inv(ii,6)*matgeo(ii,50)*F_inv(ii,8)
     &           + F_inv(ii,6)*matgeo(ii,51)*F_inv(ii,9)
      end do
      do ii = 1, span
        A4(ii,43) = F_inv(ii,4)*matgeo(ii,34)*F_inv(ii,1)
     &           + F_inv(ii,4)*matgeo(ii,35)*F_inv(ii,2)
     &           + F_inv(ii,4)*matgeo(ii,36)*F_inv(ii,3)
     &           + F_inv(ii,5)*matgeo(ii,43)*F_inv(ii,1)
     &           + F_inv(ii,5)*matgeo(ii,44)*F_inv(ii,2)
     &           + F_inv(ii,5)*matgeo(ii,45)*F_inv(ii,3)
     &           + F_inv(ii,6)*matgeo(ii,52)*F_inv(ii,1)
     &           + F_inv(ii,6)*matgeo(ii,53)*F_inv(ii,2)
     &           + F_inv(ii,6)*matgeo(ii,54)*F_inv(ii,3)
      end do
      do ii = 1, span
        A4(ii,44) = F_inv(ii,4)*matgeo(ii,34)*F_inv(ii,4)
     &           + F_inv(ii,4)*matgeo(ii,35)*F_inv(ii,5)
     &           + F_inv(ii,4)*matgeo(ii,36)*F_inv(ii,6)
     &           + F_inv(ii,5)*matgeo(ii,43)*F_inv(ii,4)
     &           + F_inv(ii,5)*matgeo(ii,44)*F_inv(ii,5)
     &           + F_inv(ii,5)*matgeo(ii,45)*F_inv(ii,6)
     &           + F_inv(ii,6)*matgeo(ii,52)*F_inv(ii,4)
     &           + F_inv(ii,6)*matgeo(ii,53)*F_inv(ii,5)
     &           + F_inv(ii,6)*matgeo(ii,54)*F_inv(ii,6)
      end do
      do ii = 1, span
        A4(ii,45) = F_inv(ii,4)*matgeo(ii,34)*F_inv(ii,7)
     &           + F_inv(ii,4)*matgeo(ii,35)*F_inv(ii,8)
     &           + F_inv(ii,4)*matgeo(ii,36)*F_inv(ii,9)
     &           + F_inv(ii,5)*matgeo(ii,43)*F_inv(ii,7)
     &           + F_inv(ii,5)*matgeo(ii,44)*F_inv(ii,8)
     &           + F_inv(ii,5)*matgeo(ii,45)*F_inv(ii,9)
     &           + F_inv(ii,6)*matgeo(ii,52)*F_inv(ii,7)
     &           + F_inv(ii,6)*matgeo(ii,53)*F_inv(ii,8)
     &           + F_inv(ii,6)*matgeo(ii,54)*F_inv(ii,9)
      end do
      do ii = 1, span
        A4(ii,46) = F_inv(ii,7)*matgeo(ii,28)*F_inv(ii,1)
     &           + F_inv(ii,7)*matgeo(ii,29)*F_inv(ii,2)
     &           + F_inv(ii,7)*matgeo(ii,30)*F_inv(ii,3)
     &           + F_inv(ii,8)*matgeo(ii,37)*F_inv(ii,1)
     &           + F_inv(ii,8)*matgeo(ii,38)*F_inv(ii,2)
     &           + F_inv(ii,8)*matgeo(ii,39)*F_inv(ii,3)
     &           + F_inv(ii,9)*matgeo(ii,46)*F_inv(ii,1)
     &           + F_inv(ii,9)*matgeo(ii,47)*F_inv(ii,2)
     &           + F_inv(ii,9)*matgeo(ii,48)*F_inv(ii,3)
      end do
      do ii = 1, span
        A4(ii,47) = F_inv(ii,7)*matgeo(ii,28)*F_inv(ii,4)
     &           + F_inv(ii,7)*matgeo(ii,29)*F_inv(ii,5)
     &           + F_inv(ii,7)*matgeo(ii,30)*F_inv(ii,6)
     &           + F_inv(ii,8)*matgeo(ii,37)*F_inv(ii,4)
     &           + F_inv(ii,8)*matgeo(ii,38)*F_inv(ii,5)
     &           + F_inv(ii,8)*matgeo(ii,39)*F_inv(ii,6)
     &           + F_inv(ii,9)*matgeo(ii,46)*F_inv(ii,4)
     &           + F_inv(ii,9)*matgeo(ii,47)*F_inv(ii,5)
     &           + F_inv(ii,9)*matgeo(ii,48)*F_inv(ii,6)
      end do
      do ii = 1, span
        A4(ii,48) = F_inv(ii,7)*matgeo(ii,28)*F_inv(ii,7)
     &           + F_inv(ii,7)*matgeo(ii,29)*F_inv(ii,8)
     &           + F_inv(ii,7)*matgeo(ii,30)*F_inv(ii,9)
     &           + F_inv(ii,8)*matgeo(ii,37)*F_inv(ii,7)
     &           + F_inv(ii,8)*matgeo(ii,38)*F_inv(ii,8)
     &           + F_inv(ii,8)*matgeo(ii,39)*F_inv(ii,9)
     &           + F_inv(ii,9)*matgeo(ii,46)*F_inv(ii,7)
     &           + F_inv(ii,9)*matgeo(ii,47)*F_inv(ii,8)
     &           + F_inv(ii,9)*matgeo(ii,48)*F_inv(ii,9)
      end do
      do ii = 1, span
        A4(ii,49) = F_inv(ii,7)*matgeo(ii,31)*F_inv(ii,1)
     &           + F_inv(ii,7)*matgeo(ii,32)*F_inv(ii,2)
     &           + F_inv(ii,7)*matgeo(ii,33)*F_inv(ii,3)
     &           + F_inv(ii,8)*matgeo(ii,40)*F_inv(ii,1)
     &           + F_inv(ii,8)*matgeo(ii,41)*F_inv(ii,2)
     &           + F_inv(ii,8)*matgeo(ii,42)*F_inv(ii,3)
     &           + F_inv(ii,9)*matgeo(ii,49)*F_inv(ii,1)
     &           + F_inv(ii,9)*matgeo(ii,50)*F_inv(ii,2)
     &           + F_inv(ii,9)*matgeo(ii,51)*F_inv(ii,3)
      end do
      do ii = 1, span
        A4(ii,50) = F_inv(ii,7)*matgeo(ii,31)*F_inv(ii,4)
     &           + F_inv(ii,7)*matgeo(ii,32)*F_inv(ii,5)
     &           + F_inv(ii,7)*matgeo(ii,33)*F_inv(ii,6)
     &           + F_inv(ii,8)*matgeo(ii,40)*F_inv(ii,4)
     &           + F_inv(ii,8)*matgeo(ii,41)*F_inv(ii,5)
     &           + F_inv(ii,8)*matgeo(ii,42)*F_inv(ii,6)
     &           + F_inv(ii,9)*matgeo(ii,49)*F_inv(ii,4)
     &           + F_inv(ii,9)*matgeo(ii,50)*F_inv(ii,5)
     &           + F_inv(ii,9)*matgeo(ii,51)*F_inv(ii,6)
      end do
      do ii = 1, span
        A4(ii,51) = F_inv(ii,7)*matgeo(ii,31)*F_inv(ii,7)
     &           + F_inv(ii,7)*matgeo(ii,32)*F_inv(ii,8)
     &           + F_inv(ii,7)*matgeo(ii,33)*F_inv(ii,9)
     &           + F_inv(ii,8)*matgeo(ii,40)*F_inv(ii,7)
     &           + F_inv(ii,8)*matgeo(ii,41)*F_inv(ii,8)
     &           + F_inv(ii,8)*matgeo(ii,42)*F_inv(ii,9)
     &           + F_inv(ii,9)*matgeo(ii,49)*F_inv(ii,7)
     &           + F_inv(ii,9)*matgeo(ii,50)*F_inv(ii,8)
     &           + F_inv(ii,9)*matgeo(ii,51)*F_inv(ii,9)
      end do
      do ii = 1, span
        A4(ii,52) = F_inv(ii,7)*matgeo(ii,34)*F_inv(ii,1)
     &           + F_inv(ii,7)*matgeo(ii,35)*F_inv(ii,2)
     &           + F_inv(ii,7)*matgeo(ii,36)*F_inv(ii,3)
     &           + F_inv(ii,8)*matgeo(ii,43)*F_inv(ii,1)
     &           + F_inv(ii,8)*matgeo(ii,44)*F_inv(ii,2)
     &           + F_inv(ii,8)*matgeo(ii,45)*F_inv(ii,3)
     &           + F_inv(ii,9)*matgeo(ii,52)*F_inv(ii,1)
     &           + F_inv(ii,9)*matgeo(ii,53)*F_inv(ii,2)
     &           + F_inv(ii,9)*matgeo(ii,54)*F_inv(ii,3)
      end do
      do ii = 1, span
        A4(ii,53) = F_inv(ii,7)*matgeo(ii,34)*F_inv(ii,4)
     &           + F_inv(ii,7)*matgeo(ii,35)*F_inv(ii,5)
     &           + F_inv(ii,7)*matgeo(ii,36)*F_inv(ii,6)
     &           + F_inv(ii,8)*matgeo(ii,43)*F_inv(ii,4)
     &           + F_inv(ii,8)*matgeo(ii,44)*F_inv(ii,5)
     &           + F_inv(ii,8)*matgeo(ii,45)*F_inv(ii,6)
     &           + F_inv(ii,9)*matgeo(ii,52)*F_inv(ii,4)
     &           + F_inv(ii,9)*matgeo(ii,53)*F_inv(ii,5)
     &           + F_inv(ii,9)*matgeo(ii,54)*F_inv(ii,6)
      end do
      do ii = 1, span
        A4(ii,54) = F_inv(ii,7)*matgeo(ii,34)*F_inv(ii,7)
     &           + F_inv(ii,7)*matgeo(ii,35)*F_inv(ii,8)
     &           + F_inv(ii,7)*matgeo(ii,36)*F_inv(ii,9)
     &           + F_inv(ii,8)*matgeo(ii,43)*F_inv(ii,7)
     &           + F_inv(ii,8)*matgeo(ii,44)*F_inv(ii,8)
     &           + F_inv(ii,8)*matgeo(ii,45)*F_inv(ii,9)
     &           + F_inv(ii,9)*matgeo(ii,52)*F_inv(ii,7)
     &           + F_inv(ii,9)*matgeo(ii,53)*F_inv(ii,8)
     &           + F_inv(ii,9)*matgeo(ii,54)*F_inv(ii,9)
      end do
      do ii = 1, span
        A4(ii,55) = F_inv(ii,1)*matgeo(ii,55)*F_inv(ii,1)
     &           + F_inv(ii,1)*matgeo(ii,56)*F_inv(ii,2)
     &           + F_inv(ii,1)*matgeo(ii,57)*F_inv(ii,3)
     &           + F_inv(ii,2)*matgeo(ii,64)*F_inv(ii,1)
     &           + F_inv(ii,2)*matgeo(ii,65)*F_inv(ii,2)
     &           + F_inv(ii,2)*matgeo(ii,66)*F_inv(ii,3)
     &           + F_inv(ii,3)*matgeo(ii,73)*F_inv(ii,1)
     &           + F_inv(ii,3)*matgeo(ii,74)*F_inv(ii,2)
     &           + F_inv(ii,3)*matgeo(ii,75)*F_inv(ii,3)
      end do
      do ii = 1, span
        A4(ii,56) = F_inv(ii,1)*matgeo(ii,55)*F_inv(ii,4)
     &           + F_inv(ii,1)*matgeo(ii,56)*F_inv(ii,5)
     &           + F_inv(ii,1)*matgeo(ii,57)*F_inv(ii,6)
     &           + F_inv(ii,2)*matgeo(ii,64)*F_inv(ii,4)
     &           + F_inv(ii,2)*matgeo(ii,65)*F_inv(ii,5)
     &           + F_inv(ii,2)*matgeo(ii,66)*F_inv(ii,6)
     &           + F_inv(ii,3)*matgeo(ii,73)*F_inv(ii,4)
     &           + F_inv(ii,3)*matgeo(ii,74)*F_inv(ii,5)
     &           + F_inv(ii,3)*matgeo(ii,75)*F_inv(ii,6)
      end do
      do ii = 1, span
        A4(ii,57) = F_inv(ii,1)*matgeo(ii,55)*F_inv(ii,7)
     &           + F_inv(ii,1)*matgeo(ii,56)*F_inv(ii,8)
     &           + F_inv(ii,1)*matgeo(ii,57)*F_inv(ii,9)
     &           + F_inv(ii,2)*matgeo(ii,64)*F_inv(ii,7)
     &           + F_inv(ii,2)*matgeo(ii,65)*F_inv(ii,8)
     &           + F_inv(ii,2)*matgeo(ii,66)*F_inv(ii,9)
     &           + F_inv(ii,3)*matgeo(ii,73)*F_inv(ii,7)
     &           + F_inv(ii,3)*matgeo(ii,74)*F_inv(ii,8)
     &           + F_inv(ii,3)*matgeo(ii,75)*F_inv(ii,9)
      end do
      do ii = 1, span
        A4(ii,58) = F_inv(ii,1)*matgeo(ii,58)*F_inv(ii,1)
     &           + F_inv(ii,1)*matgeo(ii,59)*F_inv(ii,2)
     &           + F_inv(ii,1)*matgeo(ii,60)*F_inv(ii,3)
     &           + F_inv(ii,2)*matgeo(ii,67)*F_inv(ii,1)
     &           + F_inv(ii,2)*matgeo(ii,68)*F_inv(ii,2)
     &           + F_inv(ii,2)*matgeo(ii,69)*F_inv(ii,3)
     &           + F_inv(ii,3)*matgeo(ii,76)*F_inv(ii,1)
     &           + F_inv(ii,3)*matgeo(ii,77)*F_inv(ii,2)
     &           + F_inv(ii,3)*matgeo(ii,78)*F_inv(ii,3)
      end do
      do ii = 1, span
        A4(ii,59) = F_inv(ii,1)*matgeo(ii,58)*F_inv(ii,4)
     &           + F_inv(ii,1)*matgeo(ii,59)*F_inv(ii,5)
     &           + F_inv(ii,1)*matgeo(ii,60)*F_inv(ii,6)
     &           + F_inv(ii,2)*matgeo(ii,67)*F_inv(ii,4)
     &           + F_inv(ii,2)*matgeo(ii,68)*F_inv(ii,5)
     &           + F_inv(ii,2)*matgeo(ii,69)*F_inv(ii,6)
     &           + F_inv(ii,3)*matgeo(ii,76)*F_inv(ii,4)
     &           + F_inv(ii,3)*matgeo(ii,77)*F_inv(ii,5)
     &           + F_inv(ii,3)*matgeo(ii,78)*F_inv(ii,6)
      end do
      do ii = 1, span
        A4(ii,60) = F_inv(ii,1)*matgeo(ii,58)*F_inv(ii,7)
     &           + F_inv(ii,1)*matgeo(ii,59)*F_inv(ii,8)
     &           + F_inv(ii,1)*matgeo(ii,60)*F_inv(ii,9)
     &           + F_inv(ii,2)*matgeo(ii,67)*F_inv(ii,7)
     &           + F_inv(ii,2)*matgeo(ii,68)*F_inv(ii,8)
     &           + F_inv(ii,2)*matgeo(ii,69)*F_inv(ii,9)
     &           + F_inv(ii,3)*matgeo(ii,76)*F_inv(ii,7)
     &           + F_inv(ii,3)*matgeo(ii,77)*F_inv(ii,8)
     &           + F_inv(ii,3)*matgeo(ii,78)*F_inv(ii,9)
      end do
      do ii = 1, span
        A4(ii,61) = F_inv(ii,1)*matgeo(ii,61)*F_inv(ii,1)
     &           + F_inv(ii,1)*matgeo(ii,62)*F_inv(ii,2)
     &           + F_inv(ii,1)*matgeo(ii,63)*F_inv(ii,3)
     &           + F_inv(ii,2)*matgeo(ii,70)*F_inv(ii,1)
     &           + F_inv(ii,2)*matgeo(ii,71)*F_inv(ii,2)
     &           + F_inv(ii,2)*matgeo(ii,72)*F_inv(ii,3)
     &           + F_inv(ii,3)*matgeo(ii,79)*F_inv(ii,1)
     &           + F_inv(ii,3)*matgeo(ii,80)*F_inv(ii,2)
     &           + F_inv(ii,3)*matgeo(ii,81)*F_inv(ii,3)
      end do
      do ii = 1, span
        A4(ii,62) = F_inv(ii,1)*matgeo(ii,61)*F_inv(ii,4)
     &           + F_inv(ii,1)*matgeo(ii,62)*F_inv(ii,5)
     &           + F_inv(ii,1)*matgeo(ii,63)*F_inv(ii,6)
     &           + F_inv(ii,2)*matgeo(ii,70)*F_inv(ii,4)
     &           + F_inv(ii,2)*matgeo(ii,71)*F_inv(ii,5)
     &           + F_inv(ii,2)*matgeo(ii,72)*F_inv(ii,6)
     &           + F_inv(ii,3)*matgeo(ii,79)*F_inv(ii,4)
     &           + F_inv(ii,3)*matgeo(ii,80)*F_inv(ii,5)
     &           + F_inv(ii,3)*matgeo(ii,81)*F_inv(ii,6)
      end do
      do ii = 1, span
        A4(ii,63) = F_inv(ii,1)*matgeo(ii,61)*F_inv(ii,7)
     &           + F_inv(ii,1)*matgeo(ii,62)*F_inv(ii,8)
     &           + F_inv(ii,1)*matgeo(ii,63)*F_inv(ii,9)
     &           + F_inv(ii,2)*matgeo(ii,70)*F_inv(ii,7)
     &           + F_inv(ii,2)*matgeo(ii,71)*F_inv(ii,8)
     &           + F_inv(ii,2)*matgeo(ii,72)*F_inv(ii,9)
     &           + F_inv(ii,3)*matgeo(ii,79)*F_inv(ii,7)
     &           + F_inv(ii,3)*matgeo(ii,80)*F_inv(ii,8)
     &           + F_inv(ii,3)*matgeo(ii,81)*F_inv(ii,9)
      end do
      do ii = 1, span
        A4(ii,64) = F_inv(ii,4)*matgeo(ii,55)*F_inv(ii,1)
     &           + F_inv(ii,4)*matgeo(ii,56)*F_inv(ii,2)
     &           + F_inv(ii,4)*matgeo(ii,57)*F_inv(ii,3)
     &           + F_inv(ii,5)*matgeo(ii,64)*F_inv(ii,1)
     &           + F_inv(ii,5)*matgeo(ii,65)*F_inv(ii,2)
     &           + F_inv(ii,5)*matgeo(ii,66)*F_inv(ii,3)
     &           + F_inv(ii,6)*matgeo(ii,73)*F_inv(ii,1)
     &           + F_inv(ii,6)*matgeo(ii,74)*F_inv(ii,2)
     &           + F_inv(ii,6)*matgeo(ii,75)*F_inv(ii,3)
      end do
      do ii = 1, span
        A4(ii,65) = F_inv(ii,4)*matgeo(ii,55)*F_inv(ii,4)
     &           + F_inv(ii,4)*matgeo(ii,56)*F_inv(ii,5)
     &           + F_inv(ii,4)*matgeo(ii,57)*F_inv(ii,6)
     &           + F_inv(ii,5)*matgeo(ii,64)*F_inv(ii,4)
     &           + F_inv(ii,5)*matgeo(ii,65)*F_inv(ii,5)
     &           + F_inv(ii,5)*matgeo(ii,66)*F_inv(ii,6)
     &           + F_inv(ii,6)*matgeo(ii,73)*F_inv(ii,4)
     &           + F_inv(ii,6)*matgeo(ii,74)*F_inv(ii,5)
     &           + F_inv(ii,6)*matgeo(ii,75)*F_inv(ii,6)
      end do
      do ii = 1, span
        A4(ii,66) = F_inv(ii,4)*matgeo(ii,55)*F_inv(ii,7)
     &           + F_inv(ii,4)*matgeo(ii,56)*F_inv(ii,8)
     &           + F_inv(ii,4)*matgeo(ii,57)*F_inv(ii,9)
     &           + F_inv(ii,5)*matgeo(ii,64)*F_inv(ii,7)
     &           + F_inv(ii,5)*matgeo(ii,65)*F_inv(ii,8)
     &           + F_inv(ii,5)*matgeo(ii,66)*F_inv(ii,9)
     &           + F_inv(ii,6)*matgeo(ii,73)*F_inv(ii,7)
     &           + F_inv(ii,6)*matgeo(ii,74)*F_inv(ii,8)
     &           + F_inv(ii,6)*matgeo(ii,75)*F_inv(ii,9)
      end do
      do ii = 1, span
        A4(ii,67) = F_inv(ii,4)*matgeo(ii,58)*F_inv(ii,1)
     &           + F_inv(ii,4)*matgeo(ii,59)*F_inv(ii,2)
     &           + F_inv(ii,4)*matgeo(ii,60)*F_inv(ii,3)
     &           + F_inv(ii,5)*matgeo(ii,67)*F_inv(ii,1)
     &           + F_inv(ii,5)*matgeo(ii,68)*F_inv(ii,2)
     &           + F_inv(ii,5)*matgeo(ii,69)*F_inv(ii,3)
     &           + F_inv(ii,6)*matgeo(ii,76)*F_inv(ii,1)
     &           + F_inv(ii,6)*matgeo(ii,77)*F_inv(ii,2)
     &           + F_inv(ii,6)*matgeo(ii,78)*F_inv(ii,3)
      end do
      do ii = 1, span
        A4(ii,68) = F_inv(ii,4)*matgeo(ii,58)*F_inv(ii,4)
     &           + F_inv(ii,4)*matgeo(ii,59)*F_inv(ii,5)
     &           + F_inv(ii,4)*matgeo(ii,60)*F_inv(ii,6)
     &           + F_inv(ii,5)*matgeo(ii,67)*F_inv(ii,4)
     &           + F_inv(ii,5)*matgeo(ii,68)*F_inv(ii,5)
     &           + F_inv(ii,5)*matgeo(ii,69)*F_inv(ii,6)
     &           + F_inv(ii,6)*matgeo(ii,76)*F_inv(ii,4)
     &           + F_inv(ii,6)*matgeo(ii,77)*F_inv(ii,5)
     &           + F_inv(ii,6)*matgeo(ii,78)*F_inv(ii,6)
      end do
      do ii = 1, span
        A4(ii,69) = F_inv(ii,4)*matgeo(ii,58)*F_inv(ii,7)
     &           + F_inv(ii,4)*matgeo(ii,59)*F_inv(ii,8)
     &           + F_inv(ii,4)*matgeo(ii,60)*F_inv(ii,9)
     &           + F_inv(ii,5)*matgeo(ii,67)*F_inv(ii,7)
     &           + F_inv(ii,5)*matgeo(ii,68)*F_inv(ii,8)
     &           + F_inv(ii,5)*matgeo(ii,69)*F_inv(ii,9)
     &           + F_inv(ii,6)*matgeo(ii,76)*F_inv(ii,7)
     &           + F_inv(ii,6)*matgeo(ii,77)*F_inv(ii,8)
     &           + F_inv(ii,6)*matgeo(ii,78)*F_inv(ii,9)
      end do
      do ii = 1, span
        A4(ii,70) = F_inv(ii,4)*matgeo(ii,61)*F_inv(ii,1)
     &           + F_inv(ii,4)*matgeo(ii,62)*F_inv(ii,2)
     &           + F_inv(ii,4)*matgeo(ii,63)*F_inv(ii,3)
     &           + F_inv(ii,5)*matgeo(ii,70)*F_inv(ii,1)
     &           + F_inv(ii,5)*matgeo(ii,71)*F_inv(ii,2)
     &           + F_inv(ii,5)*matgeo(ii,72)*F_inv(ii,3)
     &           + F_inv(ii,6)*matgeo(ii,79)*F_inv(ii,1)
     &           + F_inv(ii,6)*matgeo(ii,80)*F_inv(ii,2)
     &           + F_inv(ii,6)*matgeo(ii,81)*F_inv(ii,3)
      end do
      do ii = 1, span
        A4(ii,71) = F_inv(ii,4)*matgeo(ii,61)*F_inv(ii,4)
     &           + F_inv(ii,4)*matgeo(ii,62)*F_inv(ii,5)
     &           + F_inv(ii,4)*matgeo(ii,63)*F_inv(ii,6)
     &           + F_inv(ii,5)*matgeo(ii,70)*F_inv(ii,4)
     &           + F_inv(ii,5)*matgeo(ii,71)*F_inv(ii,5)
     &           + F_inv(ii,5)*matgeo(ii,72)*F_inv(ii,6)
     &           + F_inv(ii,6)*matgeo(ii,79)*F_inv(ii,4)
     &           + F_inv(ii,6)*matgeo(ii,80)*F_inv(ii,5)
     &           + F_inv(ii,6)*matgeo(ii,81)*F_inv(ii,6)
      end do
      do ii = 1, span
        A4(ii,72) = F_inv(ii,4)*matgeo(ii,61)*F_inv(ii,7)
     &           + F_inv(ii,4)*matgeo(ii,62)*F_inv(ii,8)
     &           + F_inv(ii,4)*matgeo(ii,63)*F_inv(ii,9)
     &           + F_inv(ii,5)*matgeo(ii,70)*F_inv(ii,7)
     &           + F_inv(ii,5)*matgeo(ii,71)*F_inv(ii,8)
     &           + F_inv(ii,5)*matgeo(ii,72)*F_inv(ii,9)
     &           + F_inv(ii,6)*matgeo(ii,79)*F_inv(ii,7)
     &           + F_inv(ii,6)*matgeo(ii,80)*F_inv(ii,8)
     &           + F_inv(ii,6)*matgeo(ii,81)*F_inv(ii,9)
      end do
      do ii = 1, span
        A4(ii,73) = F_inv(ii,7)*matgeo(ii,55)*F_inv(ii,1)
     &           + F_inv(ii,7)*matgeo(ii,56)*F_inv(ii,2)
     &           + F_inv(ii,7)*matgeo(ii,57)*F_inv(ii,3)
     &           + F_inv(ii,8)*matgeo(ii,64)*F_inv(ii,1)
     &           + F_inv(ii,8)*matgeo(ii,65)*F_inv(ii,2)
     &           + F_inv(ii,8)*matgeo(ii,66)*F_inv(ii,3)
     &           + F_inv(ii,9)*matgeo(ii,73)*F_inv(ii,1)
     &           + F_inv(ii,9)*matgeo(ii,74)*F_inv(ii,2)
     &           + F_inv(ii,9)*matgeo(ii,75)*F_inv(ii,3)
      end do
      do ii = 1, span
        A4(ii,74) = F_inv(ii,7)*matgeo(ii,55)*F_inv(ii,4)
     &           + F_inv(ii,7)*matgeo(ii,56)*F_inv(ii,5)
     &           + F_inv(ii,7)*matgeo(ii,57)*F_inv(ii,6)
     &           + F_inv(ii,8)*matgeo(ii,64)*F_inv(ii,4)
     &           + F_inv(ii,8)*matgeo(ii,65)*F_inv(ii,5)
     &           + F_inv(ii,8)*matgeo(ii,66)*F_inv(ii,6)
     &           + F_inv(ii,9)*matgeo(ii,73)*F_inv(ii,4)
     &           + F_inv(ii,9)*matgeo(ii,74)*F_inv(ii,5)
     &           + F_inv(ii,9)*matgeo(ii,75)*F_inv(ii,6)
      end do
      do ii = 1, span
        A4(ii,75) = F_inv(ii,7)*matgeo(ii,55)*F_inv(ii,7)
     &           + F_inv(ii,7)*matgeo(ii,56)*F_inv(ii,8)
     &           + F_inv(ii,7)*matgeo(ii,57)*F_inv(ii,9)
     &           + F_inv(ii,8)*matgeo(ii,64)*F_inv(ii,7)
     &           + F_inv(ii,8)*matgeo(ii,65)*F_inv(ii,8)
     &           + F_inv(ii,8)*matgeo(ii,66)*F_inv(ii,9)
     &           + F_inv(ii,9)*matgeo(ii,73)*F_inv(ii,7)
     &           + F_inv(ii,9)*matgeo(ii,74)*F_inv(ii,8)
     &           + F_inv(ii,9)*matgeo(ii,75)*F_inv(ii,9)
      end do
      do ii = 1, span
        A4(ii,76) = F_inv(ii,7)*matgeo(ii,58)*F_inv(ii,1)
     &           + F_inv(ii,7)*matgeo(ii,59)*F_inv(ii,2)
     &           + F_inv(ii,7)*matgeo(ii,60)*F_inv(ii,3)
     &           + F_inv(ii,8)*matgeo(ii,67)*F_inv(ii,1)
     &           + F_inv(ii,8)*matgeo(ii,68)*F_inv(ii,2)
     &           + F_inv(ii,8)*matgeo(ii,69)*F_inv(ii,3)
     &           + F_inv(ii,9)*matgeo(ii,76)*F_inv(ii,1)
     &           + F_inv(ii,9)*matgeo(ii,77)*F_inv(ii,2)
     &           + F_inv(ii,9)*matgeo(ii,78)*F_inv(ii,3)
      end do
      do ii = 1, span
        A4(ii,77) = F_inv(ii,7)*matgeo(ii,58)*F_inv(ii,4)
     &           + F_inv(ii,7)*matgeo(ii,59)*F_inv(ii,5)
     &           + F_inv(ii,7)*matgeo(ii,60)*F_inv(ii,6)
     &           + F_inv(ii,8)*matgeo(ii,67)*F_inv(ii,4)
     &           + F_inv(ii,8)*matgeo(ii,68)*F_inv(ii,5)
     &           + F_inv(ii,8)*matgeo(ii,69)*F_inv(ii,6)
     &           + F_inv(ii,9)*matgeo(ii,76)*F_inv(ii,4)
     &           + F_inv(ii,9)*matgeo(ii,77)*F_inv(ii,5)
     &           + F_inv(ii,9)*matgeo(ii,78)*F_inv(ii,6)
      end do
      do ii = 1, span
        A4(ii,78) = F_inv(ii,7)*matgeo(ii,58)*F_inv(ii,7)
     &           + F_inv(ii,7)*matgeo(ii,59)*F_inv(ii,8)
     &           + F_inv(ii,7)*matgeo(ii,60)*F_inv(ii,9)
     &           + F_inv(ii,8)*matgeo(ii,67)*F_inv(ii,7)
     &           + F_inv(ii,8)*matgeo(ii,68)*F_inv(ii,8)
     &           + F_inv(ii,8)*matgeo(ii,69)*F_inv(ii,9)
     &           + F_inv(ii,9)*matgeo(ii,76)*F_inv(ii,7)
     &           + F_inv(ii,9)*matgeo(ii,77)*F_inv(ii,8)
     &           + F_inv(ii,9)*matgeo(ii,78)*F_inv(ii,9)
      end do
      do ii = 1, span
        A4(ii,79) = F_inv(ii,7)*matgeo(ii,61)*F_inv(ii,1)
     &           + F_inv(ii,7)*matgeo(ii,62)*F_inv(ii,2)
     &           + F_inv(ii,7)*matgeo(ii,63)*F_inv(ii,3)
     &           + F_inv(ii,8)*matgeo(ii,70)*F_inv(ii,1)
     &           + F_inv(ii,8)*matgeo(ii,71)*F_inv(ii,2)
     &           + F_inv(ii,8)*matgeo(ii,72)*F_inv(ii,3)
     &           + F_inv(ii,9)*matgeo(ii,79)*F_inv(ii,1)
     &           + F_inv(ii,9)*matgeo(ii,80)*F_inv(ii,2)
     &           + F_inv(ii,9)*matgeo(ii,81)*F_inv(ii,3)
      end do
      do ii = 1, span
        A4(ii,80) = F_inv(ii,7)*matgeo(ii,61)*F_inv(ii,4)
     &           + F_inv(ii,7)*matgeo(ii,62)*F_inv(ii,5)
     &           + F_inv(ii,7)*matgeo(ii,63)*F_inv(ii,6)
     &           + F_inv(ii,8)*matgeo(ii,70)*F_inv(ii,4)
     &           + F_inv(ii,8)*matgeo(ii,71)*F_inv(ii,5)
     &           + F_inv(ii,8)*matgeo(ii,72)*F_inv(ii,6)
     &           + F_inv(ii,9)*matgeo(ii,79)*F_inv(ii,4)
     &           + F_inv(ii,9)*matgeo(ii,80)*F_inv(ii,5)
     &           + F_inv(ii,9)*matgeo(ii,81)*F_inv(ii,6)
      end do
      do ii = 1, span
        A4(ii,81) = F_inv(ii,7)*matgeo(ii,61)*F_inv(ii,7)
     &           + F_inv(ii,7)*matgeo(ii,62)*F_inv(ii,8)
     &           + F_inv(ii,7)*matgeo(ii,63)*F_inv(ii,9)
     &           + F_inv(ii,8)*matgeo(ii,70)*F_inv(ii,7)
     &           + F_inv(ii,8)*matgeo(ii,71)*F_inv(ii,8)
     &           + F_inv(ii,8)*matgeo(ii,72)*F_inv(ii,9)
     &           + F_inv(ii,9)*matgeo(ii,79)*F_inv(ii,7)
     &           + F_inv(ii,9)*matgeo(ii,80)*F_inv(ii,8)
     &           + F_inv(ii,9)*matgeo(ii,81)*F_inv(ii,9)
      end do
c     deallocate( matgeo )
      return
 1000 format(5x,'>>> cep2A: Now checking mat+geo')
 1001 format(5x,'    local_element: ',i3)
 1002 format( 6D9.2 )
      end subroutine
