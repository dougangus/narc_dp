c........+.........+.........+.........+.........+.........+.........+..
c  File wfc_models.f contains the following 7 subroutines: sphere_wfc, 
c  halite_wfc, aij_wfc, aijkl_wfc, averps_wfc, isotens_wfc and rotate.
c
c  Last modified: March 29, 2011.
c***********************************************************************
      subroutine sphere_wfc(ix1,dx1,nx1,dx2,nx2,dx3,nx3,nx1st,a6)
c     Input elastic constants (divided by density) in Voigt notation.
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer ix1,nx1,nx2,nx3,nx1st,jx1lo,jx1hi,kx1,kx2,kx3,i,j,k,l
      integer icase,i3d,jref
      real*8 xm1,x1,dx1,x2,dx2,x3,dx3,stemp
      real*8 v0p,vsphp,vp,v0s,vsphs,vs,x1s,x2s,x3s,sc,rs
      real*8 a6(6,6,nx3mx,nx2mx,nx1sto),a_b(6,6),a_s(6,6)
      real*8 a3(3,3,3,3),a3r(3,3,3,3)

c     This bit of code is for oneway we elastic tensor storage scheme
      if(ix1.le.nx1st)then
         jx1lo=1
         jx1hi=nx1st
      else
         jx1lo=mod(ix1,nx1st)
         if(jx1lo.ne.1)then
            write(11,*)
            write(11,*)'Error in sphere_wfc.  Calling for elastic '
            write(11,*)'matrix when jx1lo ne to 1.'
            write(11,*)
            stop
         endif
         jx1lo=1
         jx1hi=nx1st
      endif

c     Maximum depth of model
      xm1=1000.d0

c     This bit of code is for oneway we checking
      if(xm1.lt.dble(nx1-1)*dx1)then
         write(11,*)
         write(11,*)'Maximum depth exceeds xm1 for linear composite'
         write(11,*)'elastic matrix in sphere.'
         write(11,*)'Check linear scheme.'
         write(11,*)
         stop
      endif

c     Get elastic constants
      call halite_wfc(a_b)
c     Rotate elastic tensor into orientation of interest
      call aijkl_wfc(a_b,a3)
      call rotate(alpha1,beta1,gamma1,a3,a3r)
      call rotate(alpha2,beta2,gamma2,a3r,a3)
      call aij_wfc(a3,a_b)
c      call averps_wfc(a3,vp,vs) !Get isotropic average P and S velocity

c     Type of inclusion
      icase=0                   !(high velocity=0, low velocity=1)
      i3d=1                     !(cylinder=0, sphere=1)
      do i=1,6
         do j=1,6
            if(icase.eq.0)a_s(i,j)=a_b(i,j)*1.05d0
            if(icase.eq.1)a_s(i,j)=a_b(i,j)*0.95d0
         enddo
      enddo

c     Set location of spherical inclusion center in global coordinates
      x2s=+500.d0
      x3s=+500.d0
      x1s=+300.d0
      sc=9.d0                   !This defines size of sphere

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
               do i=1,6
                  do j=1,6
                     a6(i,j,kx3,kx2,kx1)=(a_s(i,j)-a_b(i,j))*(2.d0/
     +                    (dexp(stemp)+dexp(-stemp)))+a_b(i,j)
                  enddo
               enddo
            enddo
         enddo
      enddo

      return
      end

c***********************************************************************
      subroutine halite_wfc(a6)
c     Elasticity of cubic symmetry (halite) after Raymer et al. (2000).
c     Elasticity is devided by density (SI units - Pa and kg/m^3)
c***********************************************************************

      implicit none

      integer i,j,k,l
      real*8 dens,a6(6,6)

      dens=2160.d0
      do i=1,6
         do j=1,6
            a6(i,j)=0.d0
         enddo
      enddo
      a6(1,1)=49.1d9/dens
      a6(1,2)=14.0d9/dens
      a6(4,4)=12.7d9/dens
      a6(2,2)=a6(1,1)
      a6(3,3)=a6(1,1)
      a6(5,5)=a6(4,4)
      a6(6,6)=a6(4,4)
      a6(1,3)=a6(1,2)
      a6(2,3)=a6(1,2)
      do i=1,6
         do j=i,6
            if(i.ne.j)a6(j,i)=a6(i,j)
         enddo
      enddo

      return
      end

c***********************************************************************
      subroutine aij_wfc(a3,a6)
c     Returns the a_ij (Voigt) given the a_ijkl.  
c***********************************************************************

      implicit none
      integer i,j
      real*8 a6(6,6),a3(3,3,3,3)

      a6(1,1)=a3(1,1,1,1)
      a6(2,2)=a3(2,2,2,2)
      a6(3,3)=a3(3,3,3,3)
      a6(1,2)=a3(1,1,2,2)
      a6(1,3)=a3(1,1,3,3)
      a6(2,3)=a3(2,2,3,3)
      a6(6,6)=a3(1,2,1,2)
      a6(5,5)=a3(1,3,1,3)
      a6(4,4)=a3(2,3,2,3)
      a6(1,4)=a3(1,1,2,3)
      a6(1,5)=a3(1,1,1,3)
      a6(1,6)=a3(1,1,1,2)
      a6(2,4)=a3(2,2,2,3)
      a6(2,5)=a3(2,2,1,3)
      a6(2,6)=a3(2,2,1,2)
      a6(3,4)=a3(3,3,2,3)
      a6(3,5)=a3(3,3,1,3)
      a6(3,6)=a3(3,3,1,2)
      a6(4,5)=a3(2,3,1,3)
      a6(4,6)=a3(2,3,1,2)
      a6(5,6)=a3(1,3,1,2)
c     Impose symmetry 
      do i=1,6
         do j=i,6
            if(i.ne.j)a6(j,i)=a6(i,j)
         enddo
      enddo

      return
      end

c***********************************************************************
      subroutine aijkl_wfc(a6,a3)
c     Returns the a_ijkl given the a_mn Voigt notation.  First, the 
c     upper triangle of a_mn is used then imposes the general 
c     symmetries that apply.  Does not assume particular crystal 
c     symmetry (e.g. hexagonal).
c***********************************************************************

      implicit none

      integer i,j,k,l
      real*8 a6(6,6),a3(3,3,3,3)

      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  a3(i,j,k,l)=0.d0
               enddo
            enddo
         enddo
      enddo
      a3(1,1,1,1)=a6(1,1)
      a3(2,2,2,2)=a6(2,2)
      a3(3,3,3,3)=a6(3,3)
      a3(1,1,2,2)=a6(1,2)
      a3(1,1,3,3)=a6(1,3)
      a3(2,2,3,3)=a6(2,3)
      a3(1,2,1,2)=a6(6,6)
      a3(1,3,1,3)=a6(5,5)
      a3(2,3,2,3)=a6(4,4)
      a3(1,1,2,3)=a6(1,4)
      a3(1,1,1,3)=a6(1,5)
      a3(1,1,1,2)=a6(1,6)
      a3(2,2,2,3)=a6(2,4)
      a3(2,2,1,3)=a6(2,5)
      a3(2,2,1,2)=a6(2,6)
      a3(3,3,2,3)=a6(3,4)
      a3(3,3,1,3)=a6(3,5)
      a3(3,3,1,2)=a6(3,6)
      a3(2,3,1,3)=a6(4,5)
      a3(2,3,1,2)=a6(4,6)
      a3(1,3,1,2)=a6(5,6)
c     Impose symmetry -- general form
      a3(2,2,1,1)=a3(1,1,2,2)
      a3(3,3,1,1)=a3(1,1,3,3)
      a3(3,3,2,2)=a3(2,2,3,3)
      a3(2,1,1,2)=a3(1,2,1,2)
      a3(1,2,2,1)=a3(1,2,1,2)
      a3(2,1,2,1)=a3(1,2,1,2)
      a3(3,1,1,3)=a3(1,3,1,3)
      a3(1,3,3,1)=a3(1,3,1,3)
      a3(3,1,3,1)=a3(1,3,1,3)
      a3(3,2,2,3)=a3(2,3,2,3)
      a3(2,3,3,2)=a3(2,3,2,3)
      a3(3,2,3,2)=a3(2,3,2,3)
      a3(1,1,3,2)=a3(1,1,2,3)
      a3(2,3,1,1)=a3(1,1,2,3)
      a3(3,2,1,1)=a3(1,1,2,3)
      a3(1,1,3,1)=a3(1,1,1,3)
      a3(1,3,1,1)=a3(1,1,1,3)
      a3(3,1,1,1)=a3(1,1,1,3)
      a3(1,1,2,1)=a3(1,1,1,2)
      a3(1,2,1,1)=a3(1,1,1,2)
      a3(2,1,1,1)=a3(1,1,1,2)
      a3(2,2,3,2)=a3(2,2,2,3)
      a3(2,3,2,2)=a3(2,2,2,3)
      a3(3,2,2,2)=a3(2,2,2,3)
      a3(2,2,3,1)=a3(2,2,1,3)
      a3(1,3,2,2)=a3(2,2,1,3)
      a3(3,1,2,2)=a3(2,2,1,3)
      a3(2,2,2,1)=a3(2,2,1,2)
      a3(1,2,2,2)=a3(2,2,1,2)
      a3(2,1,2,2)=a3(2,2,1,2)
      a3(3,3,3,2)=a3(3,3,2,3)
      a3(2,3,3,3)=a3(3,3,2,3)
      a3(3,2,3,3)=a3(3,3,2,3)
      a3(3,3,3,1)=a3(3,3,1,3)
      a3(1,3,3,3)=a3(3,3,1,3)
      a3(3,1,3,3)=a3(3,3,1,3)
      a3(3,3,2,1)=a3(3,3,1,2)
      a3(1,2,3,3)=a3(3,3,1,2)
      a3(2,1,3,3)=a3(3,3,1,2)
      a3(2,3,3,1)=a3(2,3,1,3)
      a3(3,2,1,3)=a3(2,3,1,3)
      a3(3,2,3,1)=a3(2,3,1,3)
      a3(1,3,2,3)=a3(2,3,1,3)
      a3(1,3,3,2)=a3(2,3,1,3)
      a3(3,1,2,3)=a3(2,3,1,3)
      a3(3,1,3,2)=a3(2,3,1,3)
      a3(2,3,2,1)=a3(2,3,1,2)
      a3(3,2,1,2)=a3(2,3,1,2)
      a3(3,2,2,1)=a3(2,3,1,2)
      a3(1,2,2,3)=a3(2,3,1,2)
      a3(1,2,3,2)=a3(2,3,1,2)
      a3(2,1,2,3)=a3(2,3,1,2)
      a3(2,1,3,2)=a3(2,3,1,2)
      a3(3,1,1,2)=a3(1,3,1,2)
      a3(1,3,2,1)=a3(1,3,1,2)
      a3(3,1,2,1)=a3(1,3,1,2)
      a3(1,2,1,3)=a3(1,3,1,2)
      a3(2,1,1,3)=a3(1,3,1,2)
      a3(1,2,3,1)=a3(1,3,1,2)
      a3(2,1,3,1)=a3(1,3,1,2)

      return
      end

c***********************************************************************
      subroutine averps_wfc(cc,vpccav,vsccav) 
c     From a given anisotropic elastic tensor the average isotropic P 
c     and S wave velocities are calculated.  Elasticity divided by 
c     density.
c***********************************************************************

      implicit none
      real*8 vpccav,vsccav,tr1,tr2,tr3,ciikk,cikik
      real*8 test,test1,test2,test3
      real*8 cc(3,3,3,3)

      tr1=cc(1,1,1,1)+cc(2,2,2,2)+cc(3,3,3,3)
      tr2=cc(1,1,2,2)+cc(1,1,3,3)+cc(2,2,3,3)
      tr3=cc(2,3,2,3)+cc(1,3,1,3)+cc(1,2,1,2)
      ciikk=tr1+2.d0*tr2
      cikik=tr1+2.d0*tr3
      vsccav=(3.d0*cikik-ciikk)/30.d0
      vpccav=(3.d0*ciikk+cikik)/30.d0+vsccav
      vsccav=dsqrt(vsccav)
      vpccav=dsqrt(vpccav)
c     Test energy constraints
      if(cc(1,1,1,1).lt.0.d0)
     .     print *,'error energy 1 in averps_wfc'
      if(cc(2,2,2,2).lt.0.d0)
     .     print *,'error energy 2 in averps_wfc'
      if(cc(3,3,3,3).lt.0.d0)
     .     print *,'error energy 3 in averps_wfc'
      if(cc(2,3,2,3).lt.0.d0)
     .     print *,'error energy 4 in averps_wfc'
      if(cc(1,2,1,2).lt.0.d0)
     .     print *,'error energy 5 in averps_wfc'
      if(cc(1,3,1,3).lt.0.d0)
     .     print *,'error energy 6 in averps_wfc'
      test1=cc(1,1,1,1)*cc(2,2,2,2)-cc(1,1,2,2)**2
      if(test1.lt.0.d0)
     .     print *,'error energy 7 in averps_wfc'
      test2=cc(1,1,2,2)*cc(2,2,3,3)
     .     -cc(2,2,2,2)*cc(1,1,3,3)
      test3=cc(1,1,1,1)*cc(2,2,3,3)
     .     -cc(1,1,2,2)*cc(1,1,3,3)
      test=test1*cc(3,3,3,3) + test2*cc(1,1,3,3)
     .     - test3*cc(2,2,3,3)
      if(test.lt.0.d0)
     .     print *,'error energy 8 in averps_wfc'

      return
      end

c***********************************************************************
      subroutine isotens_wfc(fast,slow,as6)
c     Elasticity of an isotropic medium.  Elasticity is divided by 
c     density.
c***********************************************************************

      implicit none
      integer i,j
      real*8 fast,slow,a(6,6),as6(6,6)

      a(1,1)=fast**2.d0
      a(4,4)=slow**2.d0
      a(1,2)=a(1,1)-2.d0*a(4,4)
      do i=1,6
         do j=1,6
            as6(i,j)=0.d0
         enddo
      enddo
      as6(1,1)=a(1,1)
      as6(2,2)=a(1,1)
      as6(3,3)=a(1,1)
      as6(1,2)=a(1,2)
      as6(1,3)=a(1,2)
      as6(2,3)=a(1,2)
      as6(4,4)=a(4,4)
      as6(6,6)=a(4,4)
      as6(5,5)=a(4,4)
      as6(2,1)=as6(1,2)          ! impose symmetry
      as6(3,1)=as6(1,3)
      as6(3,2)=as6(2,3)

      return
      end

c***********************************************************************
      subroutine rotate(alpha,beta,gam,cc,cct)
c     This subroutine performs a transformation of the tensor c_ijkl of 
c     elastic constants by rotation about three angles (alpha, beta and
c     gamma)
c     alpha : rotation about x_2
c     beta  : rotation about x_3
c     gamma : rotation about x_1
c     Note that the sequence of the rotation is important (AB.ne.BA)
c     in this case we rotate about x_2 first, than x_3 and finally x_1.
c     (Note all rotations are in a clockwise sense).
c***********************************************************************

      implicit none

      integer i,j,k,l,m,n,r,s
      real*8 alpha,beta,gam,csum
      real*8 cc(3,3,3,3),cct(3,3,3,3),a(3,3)

      a(1,1) = dcos(alpha)*dcos(beta)
      a(1,2) = dsin(beta)
      a(1,3) =-dsin(alpha)*dcos(beta)
      a(2,1) =-dcos(gam)*dsin(beta)*dcos(alpha)+dsin(gam)*dsin(alpha)
      a(2,2) = dcos(gam)*dcos(beta)
      a(2,3) = dcos(gam)*dsin(beta)*dsin(alpha)+dsin(gam)*dcos(alpha)
      a(3,1) = dsin(gam)*dsin(beta)*dcos(alpha)+dcos(gam)*dsin(alpha)
      a(3,2) =-dsin(gam)*dcos(beta)
      a(3,3) =-dsin(gam)*dsin(beta)*dsin(alpha)+dcos(gam)*dcos(alpha)

c     c_ijkl ---> c_mnrs

      do m=1,3
        do n=1,3
          do r=1,3
            do s=1,3
              csum = 0.d0
              do i=1,3
                do j=1,3
                  do k=1,3
                    do l=1,3
                      csum = csum + 
     +                      a(m,i)*a(n,j)*a(r,k)*a(s,l)*cc(i,j,k,l)
                   enddo
                 enddo
               enddo
             enddo
              if(dabs(csum).lt.1.d0) csum = 0.d0
              cct(m,n,r,s) = csum
           enddo
         enddo
       enddo
      enddo

      return
      end
