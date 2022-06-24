c...+....|....+....|....+....|....+....|....+....|....+....|....+....|..
c  Program "NARC" (Narrow-Angle Rectangular Coordinates) performs 
c  wavefield extrapolation using the finite-difference narrow-angle 
c  one-way wave equation.  Theory: "The gap between seismic ray theory 
c  and full wavefield extrapolation", C. J. Thomson, GJI, Vol. 137, 
c  364-380.  This is a double-precision implementation and stripped 
c  down version of the code used in Angus et al. (2004), Angus (2005) 
c  and Angus & Thomson (2006).
c
c  Copyright (c) 2007 D.A. Angus.
c  All rights reserved by the author(s).
c  Last modified: July 16, 2007.
c
c  Uses a frequency-domain finite-difference scheme.
c***********************************************************************
      program narc
c***********************************************************************

      implicit none
      include '../Input/narc_dp.par'

      integer ix1,jx1,nt,nx1,nx1st,nx2,nx3,iref,mwid,nord,i,j

      real*8 dt,dx1,dx2,dx3,omwin,rpeak,beta
      real*8 wvf0(ntmx),rfilt(ntmx)
      real*8 tinc(nx2mx,nx3mx)
      real*8 p1inc(3,nx2mx,nx3mx),p2inc(3,nx2mx,nx3mx),
     +     p3inc(3,nx2mx,nx3mx)
      real*8 a6(6,6,nx3mx,nx2mx,nx1sto),a3(3,3,3,3),a(6,6)
      real*8 sigvec1(3,nx3mx)

      complex*16 uvec(3,ntmx,nx3mx,nx2mx,3)

c     Variables for simulation timing
      real*4 timediff,timearry(2)

c     Program starts

      timediff=dtime(timearry)

      write(*,*)
      write(*,*)'**************************************************'
      write(*,*)'****           Program Narc                   ****'
      write(*,*)'****                                          ****'
      write(*,*)'****    (c) Original code:                    ****'
      write(*,*)'****          D.A. Angus & C.J. Thomson       ****'
      write(*,*)'****          Queens University, 2004.        ****'
      write(*,*)'****    (c) Modifications: D.A. Angus,        ****'
      write(*,*)'****         University of Bristol, 2006.     ****'
      write(*,*)'**************************************************'
      write(*,*)
      write(*,*)'Stage:'

c     Input/Output files

      write(*,*)'       (1/6) Opening I/O files.'
      open(10,file='./Input/narc_dp.inp')
      open(20,file='./Log/incwave.out')
      open(11,file='./Log/narc.log')

      write(*,*)'       (2/6) Reading extrapolation input parameters.'
      call ingrid(nt,dt,nx1,dx1,nx2,dx2,nx3,dx3,iref,omwin,mwid,
     +     nord,rpeak,beta)
      call inoutfiles

      write(*,*)'       (3/6) Reading elastic file information.'
      ix1=1
      jx1=-1                    !Tag to indicate IC call
      nx1st=nx1sto              !Only need to store initial plane ec's
      call edo(ix1,jx1,nx1st,nx1,nx2,nx3,nt,dt,dx1,dx2,dx3,a6)

c      do i=1,6
c         do j=i,6
c            write(*,*)i,j,a6(i,j,1,1,1)
c         enddo
c      enddo
c      write(*,*)i,j,den
c      read(*,*)

      write(*,*)'       (4/6) Reading incident wavefield information.'
      call setincid(tinc,nx1st,nx2,nx3,dx1,dx2,dx3,a6,uvec,wvf0,
     +     nt,dt,p1inc,p2inc,p3inc,mwid,nord,rpeak,sigvec1,beta)
      close(10)
      close(20)

c      write(*,*)'       (5/6) Evaluating FD stability and dispersion.'
c      do i=1,6
c         do j=1,6
c            a(i,j)=a6(i,j,nx3/2+1,nx2/2+1,1)
c         enddo
c      enddo
c      call fd_ana(nx1,dx1,nx2,dx2,nx3,dx3,nt,dt,omwin,rpeak,
c     +     p1inc,p2inc,p3inc,a,iref)

      write(*,*)'       (6/6) Performing wavefield extrapolation.'
      call extrcart(nt,dt,omwin,nx1st,nx1,dx1,nx2,dx2,nx3,dx3,a6,
     +     uvec,iref,p2inc,p3inc)

      timediff=dtime(timearry)

      write(*,*)
      write(*,*)'       Cartesian extrapolation completed'
      write(*,*)
      write(*,*)'       Waveform time parameters (nt,dt):',nt,dt,'.'
      write(*,*)'       Number of x1 planes=',nx1,' at dx1=',dx1,'.'
      write(*,*)'       Grid dimensions (nx2,nx3):',nx2,',',nx3,'.'
      write(*,*)'       Lateral increments (dx2,dx3):',dx2,',',dx3,'.'
      write(*,*)
      write(*,*)'       CPU time: ',int(timediff),' seconds.'

      close(11)

      stop
      end
