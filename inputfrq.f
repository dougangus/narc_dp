c...+....|....+....|....+....|....+....|....+....|....+....|....+....|..
c  File Inputfrq.f contains the following 3 subroutines: Ingrid, 
c  Incwft and InOutFiles.
c
c  Last modified: July 26, 2007.
c***********************************************************************
      subroutine ingrid(nt,dt,nx1,dx1,nx2,dx2,nx3,dx3,iref,omwin,mwid,
     +     nord,rpeak,beta)
c     Input various basic parameters of the spatial and temporal grids,
c     the source and the medium.
c***********************************************************************

      implicit none
      include '../Input/narc_dp.par'

      integer i,nt,nx1,nx2,nx3,iref,mwid,nord
      real*8 dt,dx1,dx2,dx3,beta,omwin,rpeak

c     Grid parameters

      read(10,*)            !Grid input parameters
      read(10,*) nt,dt
      read(10,*) nx1,nx2,nx3
      read(10,*) dx1,dx2,dx3
      read(10,*) beta
      read(10,*)            !Spatial filtering margins
      read(10,*) imargin,jmargin
      read(10,*)            !Homogeneous model input format
      read(10,*) itropic,inhomog,iuse
      read(10,*)            !Simple isotropic, homogeneous parameters
      read(10,*) vfast,vslow,den
      read(10,*)            !Anisotropic model parameters
      read(10,*) ielas,ieff
      read(10,*) alpha1,beta1,gamma1
      read(10,*) alpha2,beta2,gamma2
      read(10,*)            !Heterogeneous model parameters
      read(10,*) ianalyt
      read(10,*) imodel,imtrick
      read(10,*) x1o,x2o,x3o
      read(10,*)            !Waveform parameters
      read(10,*) omwin 
      read(10,*) iricker
      read(10,*) mwid 
      read(10,*) nord
      read(10,*) rpeak
      read(10,*) iref

      do i=1,13
         read(10,*)
      enddo
      read(10,*)ipol
      rewind(10)
      do i=1,26
         read(10,*)
      enddo

c      if(inhomog.eq.1.and.ianalyt.eq.1.and.imodel.eq.6)dx1=-dx1

c     Output parameters

      read(10,*)            !Output parameters
      read(10,*) inumber,istart,increment

c     First rotation angles
      alpha1=alpha1*4.d0*pi4/180.d0
      beta1=beta1*4.d0*pi4/180.d0
      gamma1=gamma1*4.d0*pi4/180.d0
c     Second rotation angles
      alpha2=alpha2*4.d0*pi4/180.d0
      beta2=beta2*4.d0*pi4/180.d0
      gamma2=gamma2*4.d0*pi4/180.d0

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
      include '../Input/narc_dp.par'

      integer iform
      real*8 temp1,temp2,temp3,temp4

      read(10,*)                        !Header line
      read(10,*)iform

      if(iform.eq.0)then                !Read in angles
         read(10,*) temp1,temp2
         read(10,*) temp3,temp4
         read(10,*)
         read(10,*)
      elseif(iform.eq.1)then            !Read in slownesses
         read(10,*)
         read(10,*)
         read(10,*)temp1,temp2
         read(10,*)temp3,temp4
      endif

      read(10,*)iwave

      read(10,*)                        !Initial polarization of
      read(10,*)ipol                    !body-waves
      read(10,*)s1phi

      s1phi=s1phi*4.d0*pi4/180.d0
      nplane=0
      if(temp3.eq.0.d0.and.temp4.eq.0.d0)nplane=1

      return
      end

c***********************************************************************
      subroutine inoutfiles
c     Inputs relevant Atrakd parameters. 
c***********************************************************************

      implicit none
      include '../Input/narc_dp.par'

      read(10,*)                        !Header line
      read(10,*)waveout
      read(10,*)modelcart_e

      return
      end
