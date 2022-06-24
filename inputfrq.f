c...+....|....+....|....+....|....+....|....+....|....+....|....+....|..
c  File Inputfrq.f contains the following 3 subroutines: Ingrid, 
c  Incwft and InOutFiles.
c
c  Last modified: May 23, 2012.
c***********************************************************************
      subroutine ingrid(nt,dt,nx1,dx1,nx2,dx2,nx3,dx3,iref,omwin,mwid,
     +     nord,rpeak,cfreq,beta,ifdstab,isacout)
c     Input various basic parameters of the spatial and temporal grids,
c     the source and the medium.
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer i,nt,nx1,nx2,nx3,iref,mwid,nord,ifdstab,isacout
      real*8 dt,dx1,dx2,dx3,beta,omwin,rpeak,cfreq
      character btemp*66

c     Grid parameters
      read(10,*)
      read(10,*)                !Grid input parameters
      read(10,100)btemp,nt,dt
      read(10,101)btemp,nx1,nx2,nx3
      read(10,102)btemp,dx1,dx2,dx3
      read(10,103)btemp,beta
      read(10,105)btemp,ifdstab
      read(10,*)                !Spatial filtering margins
      read(10,104)btemp,imargin,jmargin
      read(10,*)                !Homogeneous model input format
      read(10,101)btemp,itropic,inhomog,iuse
      read(10,*)                !Simple isotropic, homogeneous parameters
      read(10,102)btemp,vfast,vslow,den
      read(10,*)                !Anisotropic model parameters
      read(10,104)btemp,ielas,ieff
      read(10,102)btemp,alpha1,beta1,gamma1
      read(10,102)btemp,alpha2,beta2,gamma2
      read(10,*)                !Heterogeneous model parameters
      read(10,105)btemp,ianalyt
      read(10,104)btemp,imodel,imtrick
      read(10,102)btemp,x1o,x2o,x3o
      read(10,*)                !Waveform parameters
      read(10,103)btemp,omwin 
      read(10,105)btemp,iricker
      read(10,105)btemp,mwid 
      read(10,105)btemp,nord
      read(10,106)btemp,rpeak,cfreq
      read(10,105)btemp,iref
      do i=1,13
         read(10,*)
      enddo
      read(10,105)btemp,ipol
      rewind(10)
      do i=1,28
         read(10,*)
      enddo

c      if(inhomog.eq.1.and.ianalyt.eq.1.and.imodel.eq.6)dx1=-dx1

c     Output parameters

      read(10,*)            !Output parameters
      read(10,101)btemp,inumber,istart,increment
      read(10,105)btemp,isacout

c     First rotation angles
      alpha1=alpha1*4.d0*pi4/180.d0
      beta1=beta1*4.d0*pi4/180.d0
      gamma1=gamma1*4.d0*pi4/180.d0
c     Second rotation angles
      alpha2=alpha2*4.d0*pi4/180.d0
      beta2=beta2*4.d0*pi4/180.d0
      gamma2=gamma2*4.d0*pi4/180.d0

 100  format(a66,1x,i4,1x,f16.8)
 101  format(a66,1x,3(i5,1x))
 102  format(a66,1x,3(f16.8))
 103  format(a66,1x,f16.8)
 104  format(a66,1x,2(i4,1x))
 105  format(a66,1x,i4)
 106  format(a66,1x,2(f16.8))

      return
      end

c***********************************************************************
      subroutine incwft(iform,temp1,temp2,temp3,temp4)
c     Input incident wavefront data: 
c     Iform=0 ... ray azimuth, plunge, wavefront curvature (in terms of 
c     theta and phi range) and wavetype (1-P, 2-S1, 3-S2).
c     Iform=1 ... ray horizontal slowness and slowness range and 
c     wavetype (1-P, 2-S1, 3-S2).
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer iform
      real*8 temp1,temp2,temp3,temp4
      character btemp*66

      read(10,*)                !Header line
      read(10,105)btemp,iform

      if(iform.eq.0)then        !Read in angles
         read(10,106)btemp,temp1,temp2
         read(10,106)btemp,temp3,temp4
         read(10,*)
         read(10,*)
      elseif(iform.eq.1)then    !Read in slownesses
         read(10,*)
         read(10,*)
         read(10,106)btemp,temp1,temp2
         read(10,106)btemp,temp3,temp4
      endif

      read(10,105)btemp,iwave
      read(10,105)btemp,ipol    !body-waves
      read(10,103)btemp,s1phi

      s1phi=s1phi*4.d0*pi4/180.d0
      nplane=0
      if(temp3.eq.0.d0.and.temp4.eq.0.d0)nplane=1

 103  format(a66,1x,f16.8)
 105  format(a66,1x,i4)
 106  format(a66,1x,2(f16.8,1x))

      return
      end

c***********************************************************************
      subroutine inoutfiles
c     Inputs relevant Atrakd parameters. 
c***********************************************************************

      implicit none
      include 'narc_dp.par'
      character btemp*66

      read(10,*)                !Header line
      read(10,108)btemp,waveout
      read(10,108)btemp,modelcart_e

 108  format(a66,1x,a100)

      return
      end
