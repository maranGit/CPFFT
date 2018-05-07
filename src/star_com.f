
c     ****************************************************************
c     *                                                              *
c     *                      subroutine star_com                     *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 06/27/2014 rhd             *
c     *                                                              *
c     *     interprets the special star commands:                    *
c     *     commands that are preceeded by an astrick and are        *
c     *     trapped directly by scan.                                *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine star_com
      implicit integer (a-z)             
      include 'common.main'
      double precision
     &  dumd
      real t1, wcputime, dumr
c      external wcputime
      character(len=1) :: dums
      logical promsw,echosw,comsw,atrdsw,eolsw,eofsw,menusw,ptsw,
     &     signsw
      logical matchs, debug
      data debug /.false./
      if (debug) write (out,*) '>>>>>>> inside star_com'
c                       
c               command: input -- getting input from file
c
c      if(matchs('input',5)) then
c         call infile
c
c               command: output -- writing output to a file
c
c      else if (matchs('output',6)) then
c         call outfil
c
c               command: time -- writes out total wall time
c
c      else if(matchs('time',4)) then
c         t1 = wcputime(1)
c         call errmsg (182,dum,dums,t1,dumd)
c
c               command: reset -- resets after a fatal error
c
c      else if(matchs('reset',5)) then
c         input_ok = .true.
c         call errmsg(184,dum,dums,dumr,dumd)
c
c               command: echo -- sets the scan echo on or off
c
c      else if(matchs('echo',4)) then
c         nblank= 20
c         reclen= 80
c         endchr= 1h$
c         promsw= .false.
c         comsw= .false.
c         atrdsw= .false.
c         eolsw= .true.
c         eofsw= .true.
c         menusw= .false.
c         ptsw= .false.
c         signsw= .false.
c         call scinit(nblank,reclen,endchr,promsw,echosw,comsw,atrdsw,
c     &        eolsw,eofsw,menusw,ptsw,signsw)
c     
c         if (matchs('off',3)) then
c            echosw= .false.
c            call scinit(nblank,reclen,endchr,promsw,echosw,comsw,
c     &           atrdsw,eolsw,eofsw,menusw,ptsw,signsw)
c         else 
c            echosw= .true.
c            call scinit(nblank,reclen,endchr,promsw,echosw,comsw,
c     &           atrdsw,eolsw,eofsw,menusw,ptsw,signsw)
c         endif
c
c      else
c         call errmsg (206,dum,dums,dumr,dumd)
c      endif
c
 9999 continue
      if (debug) write (out,*) '<<<<<< leaving star_com'
      return
      end



