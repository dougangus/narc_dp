c...+....|....+....|....+....|....+....|....+....|....+....|....+....|..
c  File Inhomg_ana.f contains the following 9 subroutines: Elast_axsh, 
c  Halite_lat, Sphere, Anelast, Spetzler, Mod_1D, Recfnc, Afar 
c  and Afar_Alt.
c
c  Last modified: April 17, 2012.
c***********************************************************************
      subroutine elast_axsh(ix1,dx1,nx1,nx1st,nx2,nx3,a6)
c     Input elastic constants (divided by density) in Voigt notation.
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer ix1,nx1,nx2,nx3,nx1st,jx1lo,jx1hi,kx1,kx2,kx3,i,j,k,l
      real*8 xm1,x1,dx1,fac
      real*8 a6(6,6,nx3mx,nx2mx,nx1sto),a(6,6)
      real*8 a3(3,3,3,3),a3r(3,3,3,3)

      if(ix1.le.nx1st)then
         jx1lo=1
         jx1hi=nx1st
      else
         jx1lo=mod(ix1,nx1st)
         if(jx1lo.ne.1)then
            write(11,*)
            write(11,*)'Error in elastc_axsh.  Calling for elastic '
            write(11,*)'matrix when jx1lo ne to 1.'
            write(11,*)
            stop
         endif
         jx1lo=1
         jx1hi=nx1st
      endif

c     Maximum depth of model
      xm1=250.d0
      if(xm1.lt.dble(nx1-1)*dx1)then
         write(11,*)
         write(11,*)'Maximum depth exceeds xm1 for linear composite'
         write(11,*)'elastic matrix in elast_axsh.'
         write(11,*)'Check linear scheme.'
         write(11,*)
         stop
      endif

c     Get halite under 200% axial extension
      call halite_axi(a3)
c     Get halite under 600% (10) simple shear
      call halite_shr(a3r)

      do kx1=jx1lo,jx1hi
         do kx2=1,nx2
            do kx3=1,nx3
               x1=dble(ix1-1+kx1-1)*dx1
c     Determine linear factor (pure axial at x1=0 and pure shear at 
c     x1=xm1)
               fac=x1/xm1
c     Linear combination of extension and shear over depth range to form
c     composite elastic matrix
               do i=1,3
                  do j=1,3
                     do k=1,3
                        do l=1,3
                           a3(i,j,k,l)=(1.d0-fac)*a3(i,j,k,l)+
     +                          fac*a3r(i,j,k,l)
                        enddo
                     enddo
                  enddo
               enddo
               call averps(a3,vfast,vslow)
               call rotate(alpha1,beta1,gamma1,a3,a3r)
               call rotate(alpha2,beta2,gamma2,a3r,a3)
               if(itropic.eq.0)then
                  call isotens(a)
               elseif(itropic.eq.1)then
                  call c66(a3,a)
               endif
               do i=1,6
                  do j=1,6
                     a6(i,j,kx3,kx2,kx1)=a(i,j)
                  enddo
               enddo
            enddo
         enddo
      enddo

      return
      end

c***********************************************************************
      subroutine halite_lat(ix1,dx1,nx1,dx2,nx2,dx3,nx3,nx1st,a6)
c     Input elastic constants (divided by density) in Voigt notation.
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer ix1,nx1,nx2,nx3,nx1st,jx1lo,jx1hi,kx1,kx2,kx3,i,j,k,l
      real*8 xm1,x1,dx1,x2,dx2,x3,dx3,fac,rmax,r
      real*8 a6(6,6,nx3mx,nx2mx,nx1sto),a(6,6)
      real*8 a3(3,3,3,3),a3r(3,3,3,3)
      real*8 ax3(3,3,3,3),sh3(3,3,3,3)

      if(ix1.le.nx1st)then
         jx1lo=1
         jx1hi=nx1st
      else
         jx1lo=mod(ix1,nx1st)
         if(jx1lo.ne.1)then
            write(11,*)
            write(11,*)'Error in elastc_axsh.  Calling for elastic '
            write(11,*)'matrix when jx1lo ne to 1.'
            write(11,*)
            stop
         endif
         jx1lo=1
         jx1hi=nx1st
      endif

c     Maximum depth of model
      xm1=250.d0
      if(xm1.lt.dble(nx1-1)*dx1)then
         write(11,*)
         write(11,*)'Maximum depth exceeds xm1 for linear composite'
         write(11,*)'elastic matrix in elast_axsh.'
         write(11,*)'Check linear scheme.'
         write(11,*)
         stop
      endif

c     Get halite under 200% axial extension
      call halite_axi(ax3)
      call rotate(alpha1,beta1,gamma1,ax3,a3r)
      call rotate(alpha2,beta2,gamma2,a3r,ax3)
c     Get halite under 600% (10) simple shear
      call halite_shr(sh3)
      call rotate(alpha1,beta1,gamma1,sh3,a3r)
      call rotate(alpha2,beta2,gamma2,a3r,sh3)

      do kx1=jx1lo,jx1hi
         do kx2=1,nx2
            do kx3=1,nx3
               x1=dble(ix1-1+kx1-1)*dx1
               x2=dble(kx2-1)*dx2-dble(nx2-1)*dx2/2.d0
               x3=dble(kx3-1)*dx3-dble(nx3-1)*dx3/2.d0
               r=dsqrt(x2**2.d0+x3**2.d0)
               if(kx2.eq.1.and.kx3.eq.1)rmax=r
c     Determine linear factor (pure axial at x1=0 and pure shear at 
c     x1=xm1)
               fac=r/rmax
c     Linear combination of extension and shear over lateral range to 
c     form composite elastic matrix (cylindrically symmetric)
               do i=1,3
                  do j=1,3
                     do k=1,3
                        do l=1,3
                           a3(i,j,k,l)=(1.d0-fac)*ax3(i,j,k,l)+
     +                          fac*sh3(i,j,k,l)
                        enddo
                     enddo
                  enddo
               enddo
               call averps(a3,vfast,vslow)
               if(itropic.eq.0)then
                  call isotens(a)
               elseif(itropic.eq.1)then
                  call c66(a3,a6(1,1,kx3,kx2,kx1))
               endif
               do i=1,6
                  do j=1,6
                     a6(i,j,kx3,kx2,kx1)=a(i,j)
                  enddo
               enddo
            enddo
         enddo
      enddo

      return
      end

c***********************************************************************
      subroutine sphere(ix1,dx1,nx1,dx2,nx2,dx3,nx3,nx1st,a6)
c     Input elastic constants (divided by density) in Voigt notation.
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer ix1,nx1,nx2,nx3,nx1st,jx1lo,jx1hi,kx1,kx2,kx3,i,j
      integer icase,i3d,jref
      real*8 xm1,x1,dx1,x2,dx2,x3,dx3,stemp
      real*8 v0p,vsphp,vp,v0s,vsphs,vs,x1s,x2s,x3s,sc,rs
      real*8 a6(6,6,nx3mx,nx2mx,nx1sto),a(6,6)
      real*8 a3(3,3,3,3)

      if(ix1.le.nx1st)then
         jx1lo=1
         jx1hi=nx1st
      else
         jx1lo=mod(ix1,nx1st)
         if(jx1lo.ne.1)then
            write(11,*)
            write(11,*)'Error in sphere.  Calling for elastic '
            write(11,*)'matrix when jx1lo ne to 1.'
            write(11,*)
            stop
         endif
         jx1lo=1
         jx1hi=nx1st
      endif

c     Maximum depth of model
      xm1=1000.d0
      if(xm1.lt.dble(nx1-1)*dx1)then
         write(11,*)
         write(11,*)'Maximum depth exceeds xm1 for linear composite'
         write(11,*)'elastic matrix in sphere.'
         write(11,*)'Check linear scheme.'
         write(11,*)
         stop
      endif

      call halite(a3)
      call averps(a3,vp,vs)

c     Type of spherical inclusion
c      icase=0  !(high velocity sphere)
      icase=1  !(low velocity sphere)
c      i3d=0    !(2D - cylinder)
      i3d=1    !(3D - sphere)
      jref=0

      if(icase.eq.0)then
         vsphp=vp*1.20d0
         vsphs=vs*1.20d0
         v0p=vp
         v0s=vs
      elseif(icase.eq.1)then
         vsphp=vp
         vsphs=vs
         v0p=vp*1.20d0
         v0s=vs*1.20d0
      endif

      x2s=+500.d0
      x3s=+500.d0
      x1s=+300.d0
      sc=9.d0

      do kx1=jx1lo,jx1hi
         x1=x1o+dble(ix1-1+kx1-1)*dx1
         do kx2=1,nx2
            x2=x2o+dble(kx2-1)*dx2
            do kx3=1,nx3
               x3=x3o+dble(kx3-1)*dx3
               if(i3d.eq.0)then
                  rs=dsqrt((x2-x2s)**2.d0+(x1-x1s)**2.d0)
               elseif(i3d.eq.1)then
                  rs=dsqrt((x2-x2s)**2.d0+(x3-x3s)**2.d0+(x1-x1s)**2.d0)
               endif
               stemp=(4.d0*pi4*(rs**3.d0))/(720.d0*(sc**3.d0))
               if(jref.eq.0)then
                  vfast=(vsphp-v0p)*(2.d0/
     +                 (dexp(stemp)+dexp(-stemp)))+v0p
                  vslow=(vsphs-v0s)*(2.d0/
     +                 (dexp(stemp)+dexp(-stemp)))+v0s
               elseif(jref.eq.1)then
                  vfast=v0p
                  vslow=v0s
               endif
               call isotens(a)
               do i=1,6
                  do j=1,6
                     a6(i,j,kx3,kx2,kx1)=a(i,j)
                  enddo
               enddo
            enddo
         enddo
      enddo

      return
      end

c***********************************************************************
      subroutine anelast(ix1,dx1,nx1,nx2,nx3,nx1st,a6)
c     Input elastic constants (divided by density) in Voigt notation.
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer ix1,nx1,nx2,nx3,nx1st,jx1lo,jx1hi,kx1,kx2,kx3,i,j
      integer icase,iwlay
      real*8 xm1,x1,dx1,lmodel,kmodel
      real*8 rho,z,z1,z2,dz,pfac,sfac,fac,vp,vs
      real*8 a6(6,6,nx3mx,nx2mx,nx1sto),a(6,6)
      real*8 a3(3,3,3,3)

      if(ix1.le.nx1st)then
         jx1lo=1
         jx1hi=nx1st
      else
         jx1lo=mod(ix1,nx1st)
         if(jx1lo.ne.1)then
            write(11,*)
            write(11,*)'Error in anelast.  Calling for elastic '
            write(11,*)'matrix when jx1lo ne to 1.'
            write(11,*)
            stop
         endif
         jx1lo=1
         jx1hi=nx1st
      endif

      call halite(a3)
      call averps(a3,vp,vs)
      rho=den

      lmodel=1
      kmodel=1

      do kx1=jx1lo,jx1hi
         x1=x1o+dble(ix1-1+kx1-1)*dx1
         if(lmodel.eq.1)then    !Isotropic homogeneous
            vfast=vp
            vslow=vs
         elseif(lmodel.eq.2)then !Isotropic linear increase
            xm1=2000.d0         !Maximum depth
            if(xm1.lt.dble(nx1-1)*dx1)then
               write(11,*)
               write(11,*)'Maximum depth > xm1 for linear composite'
               write(11,*)'elastic matrix in sphere.'
               write(11,*)'Check linear scheme.'
               write(11,*)
               stop
            endif
            if(kmodel.eq.1)then
               fac=0.05d0
            elseif(kmodel.eq.2)then
               fac=0.50d0
            elseif(kmodel.eq.3)then
               fac=0.05d0
            elseif(kmodel.eq.4)then
               fac=0.50d0
            endif
            if(kmodel.lt.3)then
               vfast=(1.d0+fac*(x1/xm1))*vp
               vslow=(1.d0+fac*(x1/xm1))*vs
            elseif(kmodel.gt.2)then
               vfast=vp
               vslow=vs
               den=(1.d0+fac*(x1/xm1))*rho
            endif   
         elseif(lmodel.ge.3)then !Transition model
            z=x1
            if(lmodel.eq.3)iwlay=3     !Width of transition (sharp)
            if(lmodel.eq.4)iwlay=8     !Width of transition (medium)
            if(lmodel.eq.5)iwlay=32    !Width of transition (wide)
            if(lmodel.eq.6)iwlay=1     !Width of transition (discrete)
c     Transition zones
            if(iwlay.eq.3)then
               z1=1000.d0
               z2=1030.d0
               dz=30.d0
            elseif(iwlay.eq.8)then
               z1=1000.d0
               z2=1080.d0
               dz=80.d0
            elseif(iwlay.eq.32)then
               z1=1000.d0
               z2=1320.d0
            elseif(iwlay.eq.1)then
               z1=1000.d0
               z2=1001.d0
            endif
c     Transition type
            icase=0             !(low to high velocity)
c            icase=1             !(high to low velocity)
            if(kmodel.eq.1)then
               pfac=0.05d0      !5% difference for P-wave
               sfac=0.05d0      !5% difference for S-wave
            elseif(kmodel.eq.2)then
               pfac=0.50d0      !50% difference for P-wave
               sfac=0.50d0      !50% difference for S-wave
            endif
            if(z.lt.z1)then
               fac=0.d0
               if(icase.eq.0)then
                  vfast=vp
                  vslow=vs
               elseif(icase.eq.1)then
                  vfast=vp*(1.d0+pfac)
                  vslow=vs*(1.d0+sfac)
               endif
            elseif(z.ge.z1.and.z.lt.z2)then
               fac=(z-z1)/(z2-z1)
               if(icase.eq.0)then
                  vfast=vp*(1.d0+pfac*fac)
                  vslow=vs*(1.d0+sfac*fac)
               elseif(icase.eq.1)then
                  vfast=vp*(1.d0+pfac*(1.d0-fac))
                  vslow=vs*(1.d0+sfac*(1.d0-fac))
               endif
            elseif(z.ge.z2)then
               fac=0.d0
               if(icase.eq.0)then
                  vfast=vp*(1.d0+pfac)
                  vslow=vs*(1.d0+sfac)
               elseif(icase.eq.1)then
                  vfast=vp
                  vslow=vs
               endif
            endif
         endif

         do kx2=1,nx2
            do kx3=1,nx3
               call isotens(a)
               do i=1,6
                  do j=1,6
                     a6(i,j,kx3,kx2,kx1)=a(i,j)
                  enddo
               enddo
            enddo
         enddo
      enddo

      den=rho

      return
      end

c***********************************************************************
      subroutine spetzler(ix1,dx1,dx2,nx2,dx3,nx3,nx1st,a6,ijflag)
c     1-D isotropic medium to generate caustics (from Spetzler and 
c     Snieder, 2003).  Input elastic constants (divided by density) in 
c     Voigt notation.
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer ix1,nx2,nx3,nx1st,jx1lo,jx1hi,ijflag
      integer kx1,kx2,kx3,i,j
      real*8 pi,vp,vs,rkx,epsi,xno,fac,x1,dx1,x2,dx2,x3,dx3,theta
      real*8 ctemp1,ctemp2,ctemp3,ctemp4,ctemp5,caustic
      real*8 a6(6,6,nx3mx,nx2mx,nx1sto),a(6,6)
      real*8 a3(3,3,3,3)

      pi=4.d0*pi4

      rkx=0.025d0*((dble(nx2-1)*dx2)**4.d0)
      epsi=0.035d0
      xno=0.0d0*(dble(nx2-1)*dx2)
      fac=1.d0

      if(ix1.le.nx1st)then
         jx1lo=1
         jx1hi=nx1st
      else
         jx1lo=mod(ix1,nx1st)
         if(jx1lo.ne.1)then
            write(11,*)
            write(11,*)'Error in Spetzler.  Calling for elastic '
            write(11,*)'matrix when jx1lo ne to 1.'
            write(11,*)
            stop
         endif
         jx1lo=1
         jx1hi=nx1st
      endif

      call halite(a3)
      call averps(a3,vp,vs)

      if(ix1.eq.1)open(unit=39,file='./Output/caustic.out')

      do kx1=jx1lo,jx1hi
         x1=x1o+dble(ix1-1+kx1-1)*dx1
         do kx2=1,nx2
            x2=dble(kx2-1)*dx2
            theta=(x2**4.d0)/rkx
            do kx3=1,nx3
               x3=dble(kx3-1)*dx3
               if(ijflag.eq.0)then
                  vfast=vp
                  vslow=vs
               elseif(ijflag.eq.1)then
c     Spetzler inverse slowness 1-D model
                  vfast=vp*(1.d0+dsqrt(2.d0)*epsi*fac*
     +                 dsin(((x2+xno)**4.d0)/rkx))
                  vslow=vs*(1.d0+dsqrt(2.d0)*epsi*fac*
     +                 dsin(((x2+xno)**4.d0)/rkx))
c     Caustic prediction for Spetzler inverse slowness 1-D model
                  ctemp1=vp+dsqrt(2.d0)*vp*epsi*dsin(theta)
                  ctemp2=3.d0*(x2**2.d0)*dcos(theta)*(ctemp1**(-2.d0))
                  ctemp3=(4.d0*(x2**6.d0)/rkx)*dsin(theta)*
     +                 (ctemp1**(-2.d0))
                  ctemp4=(dsqrt(128.d0)*vp*epsi*(x2**6.d0)/rkx)*
     +                 dcos(theta)*dcos(theta)*(ctemp1**(-3.d0))
                  ctemp5=(-dsqrt(32.d0)*vp*epsi/rkx)*(
     +                 ctemp2-ctemp3-ctemp4)
                  caustic=dsqrt(-2.d0/(vp*ctemp5))
               endif
               call isotens(a)
               do i=1,6
                  do j=1,6
                     a6(i,j,kx3,kx2,kx1)=a(i,j)
                  enddo
               enddo
            enddo
            if(kx1.eq.1.and.ix1.eq.1.and.ijflag.eq.1)then
               write(39,666)x2,caustic
            endif
         enddo
      enddo
 666  format(f7.1,2x,e12.5)

      vfast=vp
      vslow=vs

      if(ix1.eq.1)close(39)

      return
      end

c***********************************************************************
      subroutine mod_1d(zd,vpref,vsref,rhoref,i1d)
c     1D reference model.
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer i1d
      real*8 zd,vpref,vsref,rhoref,z1d,p2sscale

      z1d=zd/1.d3
      p2sscale=0.7d0

      if(i1d.eq.1)then       !IASPEI91
         if(z1d.lt.20.)then
            vpref=5.80
            vsref=3.36
            rhoref=2.72
         elseif((z1d.ge.20.).and.(z1d.lt.36.))then
            vpref=6.5
            vsref=3.75
            rhoref=2.92
         elseif((z1d.ge.36.).and.(z1d.le.78.))then
            vpref=8.04+((z1d-36.)/(78.-36.))*0.005
            vsref=4.47+((z1d-36.)/(78.-36.))*0.015
            rhoref=3.3198+((z1d-36.)/(78.-36.))*0.0257
         elseif((z1d.gt.78.).and.(z1d.le.120.))then
            vpref=8.045+((z1d-78.)/(120.-78.))*0.005
            vsref=4.485+((z1d-78.)/(120.-78.))*0.015
            rhoref=3.3455+((z1d-78.)/(120.-78.))*0.0258
         elseif((z1d.gt.120.).and.(z1d.lt.166.))then
            vpref=8.05+((z1d-120.)/(166.-120.))*0.125
            vsref=4.5+((z1d-120.)/(166.-120.))*0.009
            rhoref=3.3713+((z1d-120.)/(166.-120.))*0.0272
         elseif((z1d.ge.166.).and.(z1d.le.210.))then
            vpref=8.175+((z1d-166.)/(210.-166.))*0.125
            vsref=4.509+((z1d-166.)/(210.-166.))*0.009
            rhoref=3.3985+((z1d-166.)/(210.-166.))*0.0273
         elseif((z1d.gt.210.).and.(z1d.le.260.))then
            vpref=8.3+((z1d-210.)/(260.-210.))*0.1825
            vsref=4.522+((z1d-210.)/(260.-210.))*0.087
            rhoref=3.4258+((z1d-210.)/(260.-210.))*0.0303
         elseif((z1d.gt.260.).and.(z1d.le.310.))then
            vpref=8.4825+((z1d-260.)/(310.-260.))*0.1825
            vsref=4.609+((z1d-260.)/(310.-260.))*0.087
            rhoref=3.4561+((z1d-260.)/(310.-260.))*0.0303
         elseif((z1d.gt.310.).and.(z1d.le.360.))then
            vpref=8.665+((z1d-310.)/(360.-310.))*0.1825
            vsref=4.696+((z1d-310.)/(360.-310.))*0.087
            rhoref=3.4864+((z1d-310.)/(360.-310.))*0.0303
         elseif((z1d.gt.360.).and.(z1d.le.410.))then
            vpref=8.8475+((z1d-360.)/(410.-360.))*0.1825
            vsref=4.783+((z1d-360.)/(410.-360.))*0.087
            rhoref=3.5167+((z1d-360.)/(410.-360.))*0.0303
         elseif((z1d.gt.410.).and.(z1d.le.460.))then
            vpref=9.36+((z1d-410.)/(460.-410.))*0.168
            vsref=5.07+((z1d-410.)/(460.-410.))*0.106
            rhoref=3.7557+((z1d-410.)/(460.-410.))*0.0618
         elseif((z1d.gt.460.).and.(z1d.le.510.))then
            vpref=9.528+((z1d-460.)/(510.-460.))*0.168
            vsref=5.176+((z1d-460.)/(510.-460.))*0.106
            rhoref=3.8175+((z1d-460.)/(510.-460.))*0.0618
         else
            vpref=9.696
            vsref=5.282
            rhoref=3.8793
         endif
      elseif(i1d.eq.2)then   !Tilmann & Ni, 2003, Suppl. Table S3
         if(z1d.lt.6.)then
            vpref=5.80
         elseif((z1d.ge.6.).and.(z1d.lt.36.))then
            vpref=6.0+((z1d-6.)/(36.-6.))*0.1
         elseif((z1d.ge.36.).and.(z1d.le.68.))then
            vpref=6.1+((z1d-36.)/(68.-36.))*0.1
         elseif((z1d.gt.68.).and.(z1d.le.100.))then
            vpref=8.1+((z1d-68.)/(100.-68.))*0.05
         elseif((z1d.gt.100.).and.(z1d.lt.188.))then
            vpref=8.15+((z1d-102.)/(188.-102.))*0.1
         elseif((z1d.ge.188.).and.(z1d.le.210.))then
            vpref=8.25+((z1d-188.)/(210.-188.))*0.05
         elseif((z1d.gt.210.).and.(z1d.le.230.))then
            vpref=8.3+((z1d-210.)/(230.-210.))*0.07
         elseif((z1d.gt.230.).and.(z1d.le.250.))then
            vpref=8.37+((z1d-230.)/(250.-230.))*0.08
         elseif((z1d.gt.250.).and.(z1d.le.270.))then
            vpref=8.45+((z1d-250.)/(270.-250.))*0.07
         elseif((z1d.gt.270.).and.(z1d.le.290.))then
            vpref=8.52+((z1d-270.)/(290.-270.))*0.07
         elseif((z1d.gt.290.).and.(z1d.le.310.))then
            vpref=8.59+((z1d-290.)/(310.-290.))*0.07
         elseif((z1d.gt.310.).and.(z1d.le.330.))then
            vpref=8.66+((z1d-310.)/(330.-310.))*0.08
         elseif((z1d.gt.330.).and.(z1d.le.350.))then
            vpref=8.74+((z1d-330.)/(350.-330.))*0.07
         elseif((z1d.gt.350.).and.(z1d.le.370.))then
            vpref=8.81+((z1d-350.)/(370.-350.))*0.07
         elseif((z1d.gt.370.).and.(z1d.le.390.))then
            vpref=8.88+((z1d-370.)/(390.-370.))*0.08
         elseif((z1d.gt.390.).and.(z1d.le.410.))then
            vpref=8.96+((z1d-390.)/(410.-390.))*0.07
         else
            vpref=9.03
         endif
         vsref=p2sscale*vpref
      endif

      return
      end

c***********************************************************************
      subroutine recfnc(ix1,dx1,nx1,dx2,nx2,dx3,nx3,nx1st,a6)
c     Input elastic constants (divided by density) in Voigt notation.
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer ix1,nx1,nx2,nx3,nx1st,jx1lo,jx1hi,kx1,kx2,kx3,i,j,layer
      real*8 dx1,dx2,dx3,xm1,x1,x2,x3
      real*8 a6(6,6,nx3mx,nx2mx,nx1sto),a(6,6)
      real*8 a3(3,3,3,3)

      if(ix1.le.nx1st)then
         jx1lo=1
         jx1hi=nx1st
      else
         jx1lo=mod(ix1,nx1st)
         if(jx1lo.ne.1.and.nx1st.gt.1)then
            write(11,*)
            write(11,*)'Error in recfnc.  Calling for elastic '
            write(11,*)'matrix when jx1lo ne to 1.'
            write(11,*)
            stop
         endif
         jx1lo=1
         jx1hi=nx1st
      endif

c     Maximum depth of model
      xm1=450000.0
      if(xm1.lt.real(nx1-1)*dx1)then
         write(11,*)
         write(11,*)'Maximum depth exceeds xm1 for linear composite'
         write(11,*)'elastic matrix in recfnc.'
         write(11,*)'Check linear scheme.'
         write(11,*)
         stop
      endif

      do kx1=jx1lo,jx1hi

         x1=x1o+real(ix1-1+kx1-1)*dx1
         if(x1.lt.50000)layer=1
         if(x1.ge.50000.and.x1.lt.150000)layer=2
         if(x1.ge.150000.and.x1.lt.200000)layer=3
         if(x1.ge.200000.and.x1.lt.400000)layer=4
         if(x1.ge.400000.and.x1.le.450000)layer=5
         call rflayer(a3,layer)

         do kx2=1,nx2
            do kx3=1,nx3
               x2=x2o+real(kx2-1)*dx2
               x3=x3o+real(kx3-1)*dx3
               call c66(a3,a)
               do i=1,6
                  do j=1,6
                     a6(i,j,kx3,kx2,kx1)=a(i,j)
                  enddo
               enddo
            enddo
         enddo

      enddo

      return
      end

c***********************************************************************
      subroutine afar(jx1,ix1,dx2,nx2,dx3,nx3,nx1st,a6)
c     Input elastic constants (divided by density) in Voigt notation.
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer jx1,ix1,nx2,nx3,nx1st,jx1lo,jx1hi,kx1,kx2,kx3,i,j
      integer ii,jj
      real*8 dx2,dx3,x2,x3,pi,fac,ltrans,l2trans
      real*8 a6(6,6,nx3mx,nx2mx,nx1sto)
      real*8 a1(6,6),a2(6,6),a_alt(6,6),den1,den2
      real*8 a3r(3,3,3,3)
      character*100 modname2

      pi=4.d0*pi4
      ltrans=10.d3
      l2trans=ltrans/2.d0
      if(ix1.eq.1)then
         open(30,file=modelcart_e)
         do i=1,6
            do j=i,6
               read(30,*)ii,jj,a1(i,j)
               if(dabs(a1(i,j)).lt.1.d0) a1(i,j)=0.d0
               if(i.ne.j)a1(j,i)=a1(i,j)
            enddo
         enddo
         read(30,*)ii,jj,den1
         read(30,*)modname2
         close(30)
         open(30,file=modname2)
         do i=1,6
            do j=i,6
               read(30,*)ii,jj,a2(i,j)
               if(dabs(a2(i,j)).lt.1.d0) a2(i,j)=0.d0
               if(i.ne.j)a2(j,i)=a2(i,j)
            enddo
         enddo
         read(30,*)ii,jj,den2
         close(30)
         den=(den1+den2)/2.d0
c     Setup initial elasticity for incident wavefield calculation
         if(jx1.eq.-1)then
            do i=1,6
               do j=1,6
                  a_alt(i,j)=(a1(i,j)+a2(i,j))/2.d0
               enddo
            enddo
            call aijkl(a_alt,a3r)
            call averps(a3r,vfast,vslow)
            if(ipol.eq.1)call isotens(a_alt)
         endif
      endif

      if(ix1.le.nx1st)then
         jx1lo=1
         jx1hi=nx1st
      else
         jx1lo=mod(ix1,nx1st)
         if(jx1lo.ne.1.and.nx1st.gt.1)then
            write(11,*)
            write(11,*)'Error in afar.  Calling for elastic '
            write(11,*)'matrix when jx1lo ne to 1.'
            write(11,*)
            stop
         endif
         jx1lo=1
         jx1hi=nx1st
      endif

      if(ix1.eq.1)then
         do kx1=jx1lo,jx1hi
            do kx2=1,nx2
               do kx3=1,nx3
                  x2=x2o+dble(kx2-1)*dx2
                  x3=x3o+dble(kx3-1)*dx3
                  if(jx1.eq.-1)then
                     do i=1,6
                        do j=1,6
                           a6(i,j,kx3,kx2,kx1)=a_alt(i,j)
                        enddo
                     enddo
                  elseif(jx1.ne.-1)then
                     if(x2.lt.(28.d3-l2trans))then
                        do i=1,6
                           do j=1,6
                              a6(i,j,kx3,kx2,kx1)=a2(i,j)
                           enddo
                        enddo
                     elseif(x2.ge.(28.d3-l2trans).and.
     .                       x2.le.(28.d3+l2trans))then
                        fac=(1.d0+dcos((28.d3+l2trans-x2)*pi/
     .                       ltrans))/2.d0
                        do i=1,6
                           do j=1,6
                              a6(i,j,kx3,kx2,kx1)=
     .                             a2(i,j)*(1.d0-fac)+a1(i,j)*fac
                           enddo
                        enddo
                     elseif(x2.gt.(28.d3+l2trans).and.
     .                       x2.lt.(68.d3-l2trans))then
                        do i=1,6
                           do j=1,6
                              a6(i,j,kx3,kx2,kx1)=a1(i,j)
                           enddo
                        enddo
                     elseif(x2.ge.(68.d3-l2trans).and.
     .                       x2.le.(68.d3+l2trans))then
                        fac=(1.d0+dcos((68.d3+l2trans-x2)*pi/
     .                       ltrans))/2.d0
                        do i=1,6
                           do j=1,6
                              a6(i,j,kx3,kx2,kx1)=
     .                             a1(i,j)*(1.d0-fac)+a2(i,j)*fac
                           enddo
                        enddo
                     elseif(x2.gt.68.d3)then
                        do i=1,6
                           do j=1,6
                              a6(i,j,kx3,kx2,kx1)=a2(i,j)
                           enddo
                        enddo
                     endif
                  endif
               enddo
            enddo
         enddo
      endif

      return
      end

c***********************************************************************
      subroutine afar_alt(jx1,ix1,dx2,nx2,dx3,nx3,nx1st,a6)
c     Input elastic constants (divided by density) in Voigt notation.
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer jx1,ix1,nx2,nx3,nx1st,jx1lo,jx1hi,kx1,kx2,kx3,i,j
      integer ii,jj
      real*8 dx2,dx3,x2,x3,pi,ltrans,l2trans,rangle
      real*8 a6(6,6,nx3mx,nx2mx,nx1sto)
      real*8 a_0(6,6),a_30(6,6),a_alt(6,6)
      real*8 a3(3,3,3,3),a3r(3,3,3,3),angle1,angle2,angle3
      character*100 modname

      pi=4.d0*pi4

c  JH - changd l2trans, was 4d03
      ltrans=1.d0
      l2trans=ltrans/2.d0
      rangle=30.0d0

      if(ix1.le.10)then
         open(30,file=modelcart_e)
         do i=1,6
            do j=i,6
               read(30,*)ii,jj,a_0(i,j)
               if(dabs(a_0(i,j)).lt.1.d0) a_0(i,j)=0.d0
               if(i.ne.j)a_0(j,i)=a_0(i,j)
            enddo
         enddo
         read(30,*)ii,jj,den
         read(30,*)modname
         close(30)
         open(30,file=modname)
         do i=1,6
            do j=i,6
               read(30,*)ii,jj,a_30(i,j)
               if(dabs(a_30(i,j)).lt.1.d0) a_30(i,j)=0.d0
               if(i.ne.j)a_30(j,i)=a_30(i,j)
            enddo
         enddo
         read(30,*)ii,jj,den
         close(30)
		
c JH - rotating, to fit data

c      call aijkl(a_0,a3)
c      angle1=0.d0
c      angle2=pi/2.d0
c      angle3=0.d0
c      call rotate(angle1,angle2,angle3,a3,a3r)
c      call c66(a3r,a_0)
			
         call aijkl(a_30,a3)
         angle1=0.d0
         angle2=60.*(pi/180.)
         angle3=0.d0
         call rotate(angle1,angle2,angle3,a3,a3r)
         call c66(a3r,a_30)

c JH - rotating into one-way co-ordinates

         call aijkl(a_0,a3)
         angle1=pi/2.d0
         angle2=0.d0
         angle3=0.d0
         call rotate(angle1,angle2,angle3,a3,a3r)
         call c66(a3r,a_0)

         call aijkl(a_30,a3)
         angle1=pi/2.d0
         angle2=0.d0
         angle3=0.d0
         call rotate(angle1,angle2,angle3,a3,a3r)
         call c66(a3r,a_30)

c     Setup initial elasticity for incident wavefield calculation
         if(jx1.eq.-1)then
            do i=1,6
               do j=1,6
                  a_alt(i,j)=(a_0(i,j)+a_30(i,j))/2.d0
               enddo
            enddo
            call aijkl(a_alt,a3r)
            call averps(a3r,vfast,vslow)
            if(ipol.eq.1) call isotens(a_alt)
         endif
      endif

      if(ix1.le.nx1st)then
         jx1lo=1
         jx1hi=nx1st
      else
         jx1lo=mod(ix1,nx1st)
         if(jx1lo.ne.1.and.nx1st.gt.1)then
            write(11,*)
            write(11,*)'Error in afar.  Calling for elastic '
            write(11,*)'matrix when jx1lo ne to 1.'
            write(11,*)
            stop
         endif
         jx1lo=1
         jx1hi=nx1st
      endif

      if(ix1.le.10)then
         do kx1=jx1lo,jx1hi
            do kx2=1,nx2
               do kx3=1,nx3
                  x2=x2o+dble(kx2-1)*dx2
                  x3=x3o+dble(kx3-1)*dx3
                  if(jx1.eq.-1)then !Pre extrapoloation - average Cij
                     do i=1,6
                        do j=1,6
                           a6(i,j,kx3,kx2,kx1)=a_alt(i,j)
                        enddo
                     enddo
                  elseif(jx1.ne.-1)then
                     do i=1,6
                        do j=1,6
                           a6(i,j,kx3,kx2,kx1)=a_0(i,j)
                        enddo
                     enddo
                  endif
               enddo
            enddo
         enddo
      endif

      return
      end

c***********************************************************************
      subroutine karakoram(jx1,ix1,dx2,nx2,dx3,nx3,nx1st,a6)
c     Input elastic constants (divided by density) in Voigt notation.
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer jx1,ix1,nx2,nx3,nx1st,jx1lo,jx1hi,kx1,kx2,kx3,i,j
      integer ii,jj
      real*8 dx2,dx3,x2,x3,pi,rfac,ltrans,l2trans,angle,rangle
      real*8 a6(6,6,nx3mx,nx2mx,nx1sto)
      real*8 a_0(6,6),a_30(6,6),a_alt(6,6),ao(6,6)
      real*8 a3(3,3,3,3),a3r(3,3,3,3),angle1,angle2,angle3
      character*100 modname

      pi=4.d0*pi4

c  JH - changd l2trans, was 4d03
      ltrans=1.d0
      l2trans=ltrans/2.d0
      rangle=30.0d0

      if(ix1.le.10)then
         open(30,file=modelcart_e)
         do i=1,6
            do j=i,6
               read(30,*)ii,jj,a_0(i,j)
               if(dabs(a_0(i,j)).lt.1.d0) a_0(i,j)=0.d0
               if(i.ne.j)a_0(j,i)=a_0(i,j)
            enddo
         enddo
         read(30,*)ii,jj,den
         read(30,*)modname
         close(30)
         do i=1,6
            do j=1,6
               a_0(i,j)=a_0(i,j)/den
            enddo
         enddo
         open(30,file=modname)
         do i=1,6
            do j=i,6
               read(30,*)ii,jj,a_30(i,j)
               if(dabs(a_30(i,j)).lt.1.d0) a_30(i,j)=0.d0
               if(i.ne.j)a_30(j,i)=a_30(i,j)
            enddo
         enddo
         read(30,*)ii,jj,den
         close(30)
         do i=1,6
            do j=1,6
               a_30(i,j)=a_30(i,j)/den
            enddo
         enddo

c Rotating from global into one-way co-ordinates

         call aijkl(a_0,a3)
         angle1=pi/2.d0
         angle2=0.d0
         angle3=0.d0
         call rotate(angle1,angle2,angle3,a3,a3r)
         call c66(a3r,a_0)
         call aijkl(a_30,a3)
         angle1=pi/2.d0
         angle2=0.d0
         angle3=0.d0
         call rotate(angle1,angle2,angle3,a3,a3r)
         call c66(a3r,a_30)

c     Setup initial elasticity for incident wavefield calculation
         if(jx1.eq.-1)then
            do i=1,6
               do j=1,6
                  a_alt(i,j)=(a_0(i,j)+a_30(i,j))/2.d0
               enddo
            enddo
            call aijkl(a_alt,a3r)
            call averps(a3r,vfast,vslow)
            if(ipol.eq.1) call isotens(a_alt)
         endif
      endif

      if(ix1.le.nx1st)then
         jx1lo=1
         jx1hi=nx1st
      else
         jx1lo=mod(ix1,nx1st)
         if(jx1lo.ne.1.and.nx1st.gt.1)then
            write(11,*)
            write(11,*)'Error in afar.  Calling for elastic '
            write(11,*)'matrix when jx1lo ne to 1.'
            write(11,*)
            stop
         endif
         jx1lo=1
         jx1hi=nx1st
      endif

      if(ix1.le.10)then
         do kx1=jx1lo,jx1hi
            do kx2=1,nx2
               do kx3=1,nx3
                  x2=x2o+dble(kx2-1)*dx2
                  x3=x3o+dble(kx3-1)*dx3
                  if(jx1.eq.-1)then
                     do i=1,6
                        do j=1,6
                           a6(i,j,kx3,kx2,kx1)=a_alt(i,j)
                        enddo
                     enddo
                  elseif(jx1.ne.-1)then
                     do i=1,6
                        do j=1,6
                           a6(i,j,kx3,kx2,kx1)=a_30(i,j)
                        enddo
                     enddo
                  endif
               enddo
            enddo
         enddo
      endif

      return
      end
