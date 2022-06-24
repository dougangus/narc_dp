c........+.........+.........+.........+.........+.........+.........+..
c  Program Conversion: 
c  The program converts oneway waveform output from binary to ascii or 
c  ascii to binary.
c
c  Copyright (c) 2007 D.A. Angus.
c  All rights reserved by the author(s).
c  Last modified: July 25, 2007.
c***********************************************************************
      program conversion
c***********************************************************************

      implicit none

      integer ntmx,nx2mx,nx3mx
      parameter (ntmx=257,nx2mx=201,nx3mx=201)

      integer iconvert,inumber,istart,increment
      integer ip,jx1,nx1,ix2,jx2,nx2,ix3,jx3,nx3,it,nt,nom,nom2
      integer nta(nx2mx,nx3mx)
      real*8 dx1,x1o,dt
      real*8 tbg(nx2mx,nx3mx),tnd(nx2mx,nx3mx)
      complex*16 uvc(3,ntmx,nx3mx,nx2mx)

      write(*,*)
      write(*,*)'Binary to ascii conversion (yes=1 and no=0)'
      write(*,*)
      read(*,*)iconvert
      write(*,*)

c     Decide which plane to view
      open(24,file='./Output/wavefld.out',form='unformatted')
      open(25,file='./Output/wavefld.asc',form='formatted')
      if(iconvert.eq.1)then
         read(24)inumber,istart,increment,nx1
         read(24)dx1,x1o
	 write(25,*)inumber,istart,increment,nx1
	 write(25,*)dx1,x1o
      elseif(iconvert.eq.0)then
	 read(25,*)inumber,istart,increment,nx1
	 read(25,*)dx1,x1o
         write(24)inumber,istart,increment,nx1
         write(24)dx1,x1o
      endif

c     Read frequency-domain wavefield
      do ip=1,inumber+1
         if(iconvert.eq.1)then
            read(24)jx1,nx1,nx2,nx3
            read(24)nt,dt,nom,nom2
            read(25,*)jx1,nx1,nx2,nx3
            read(25,*)nt,dt,nom,nom2
	 elseif(iconvert.eq.0)then
            read(25,*)jx1,nx1,nx2,nx3
            read(25,*)nt,dt,nom,nom2
            write(24)jx1,nx1,nx2,nx3
            write(24)nt,dt,nom,nom2
	 endif
         do ix2=1,nx2
            do ix3=1,nx3
	       if(iconvert.eq.1)then
                  read(24)jx2,jx3,nta(ix2,ix3),tbg(ix2,ix3),tnd(ix2,ix3)
		  write(25,*)jx2,jx3,nta(ix2,ix3),tbg(ix2,ix3),
     +                 tnd(ix2,ix3)
               elseif(iconvert.eq.0)then
                  read(25,*)jx2,jx3,nta(ix2,ix3),tbg(ix2,ix3),
     +                 tnd(ix2,ix3)
		  write(24)jx2,jx3,nta(ix2,ix3),tbg(ix2,ix3),
     +                 tnd(ix2,ix3)
	       endif
               do it=1,nt
	          if(iconvert.eq.1)then
                     read(24)uvc(1,it,ix3,ix2),uvc(2,it,ix3,ix2),
     +                    uvc(3,it,ix3,ix2)
                     write(25,*)uvc(1,it,ix3,ix2),uvc(2,it,ix3,ix2),
     +                    uvc(3,it,ix3,ix2)
                  elseif(iconvert.eq.0)then
                     read(25,*)uvc(1,it,ix3,ix2),uvc(2,it,ix3,ix2),
     +                    uvc(3,it,ix3,ix2)
                     write(24)uvc(1,it,ix3,ix2),uvc(2,it,ix3,ix2),
     +                    uvc(3,it,ix3,ix2)
		  endif
               enddo
            enddo
         enddo
      enddo
      close(24)
      close(25)

      stop
      end
