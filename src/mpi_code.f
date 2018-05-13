c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine die_gracefully               *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 03/02/98                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine die_gracefully
c
      stop
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine die_abort                    *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 03/02/98                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine die_abort
c
      stop
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_send_simple_angles      *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 3/12/14                    *
c     *                                                              *
c     *           Send the simple angle properties, if required      *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_send_simple_angles
      implicit none
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_send_crystals           *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 8/16/12                    *
c     *                                                              *
c     *           Send all the crystal properties to the workers     *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_send_crystals
      implicit none
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_dealloc_crystals        *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 8/16/12                    *
c     *                                                              *
c     *           Dealloc crystal data structures                    *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_dealloc_crystals
      implicit none
      return
      end subroutine
