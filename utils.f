c........+.........+.........+.........+.........+.........+.........+..
c     File Utils.f contains the following 26 subroutines: Roots, 
c     Coeff, Cofac, Xrdvec, Dpoly, Rotate, Rotvec, Rotatelc, Rinv3, 
c     Cinv3, Rinv3a, Crinv3, Trans, Ctrans, Rotmat1, Rmult3a, Cmult3a, 
c     Rmult3b, Cmult3b, Rmult3c, Cmult3c, Rmult3d, Cmult3d, Uvecshft, 
c     Rtaper, and Detm.
c
c     Last modified: May 23, 2012.
c***********************************************************************
      subroutine roots(a1,a2,a3,x1,x2r,x2i,x3r,x3i)
c     Subroutine to evaluate the roots of a third order polynomial.
c***********************************************************************

      implicit none

      real*8 pi,a1,a2,a3
      real*8 q,r,q3,r2,discr,rtdisc,rin,rrin,s,t
      real*8 x1,x2r,x3r,x2i,x3i,cosphi,phi,a,b,theta

      pi=4.d0*datan(1.d0)

      q = (3.d0*a2 - a1*a1)/9.d0
      r = (9.d0*a1*a2 - 27.d0*a3 - 2.d0*a1*a1*a1)/54.d0
      q3 = q*q*q
      r2 = r*r
      discr = q3+r2
      if(discr.lt.0.d0) goto 99
      rtdisc = dsqrt(discr)
      rin = r + rtdisc
      rrin = r - rtdisc
      if(rin.ge.0.d0)  s = rin**(1.d0/3.d0)
      if(rin.lt.0.d0)  s = -((-rin)**(1.d0/3.d0))
      if(rrin.ge.0.d0)  t = rrin**(1.d0/3.d0)
      if(rrin.lt.0.d0)  t = -((-rrin)**(1.d0/3.d0))
      x1 = s+t - a1/3.d0
      x2r = -0.5d0*(s+t) - a1/3.d0
      x3r = x2r
      x2i = dsqrt(3.d0)*(s-t)/2.d0
      x3i = -x2i
      goto 999
99    cosphi = r/dsqrt(-q3)
      phi = dacos(cosphi)
      a = 2.d0*dsqrt(-q)
      b = -a1/3.d0
      theta = phi/3.d0
      x1 = a*dcos(theta) +b
      x2r = a*dcos(theta+2.d0*pi/3.d0) +b
      x2i = 0.d0
      x3r = a*dcos(theta+4.d0*pi/3.d0) +b
      x3i = 0.d0
999   return
      end

c***********************************************************************
      subroutine coeff(c,a1,a2,a3)
c     Subroutine to evaluate the coefficients of the eigenvalue
c     equation for a 3x3 matrix. x3 + a1x2 + a2x + a3 = 0
c***********************************************************************

      implicit none

      real*8 c11,c12,c21,c13,c31,c22,c23,c32,c33,a1,a2,a3
      real*8 c(3,3)

      c11 = c(1,1)
      c12 = c(1,2)
      c21 = c12
      c13 = c(1,3)
      c31 = c13
      c22 = c(2,2)
      c23 = c(2,3)
      c32 = c23
      c33 = c(3,3)

      a1 =  -(c11+c22+c33)
      a2 =   c11*c22+c22*c33+c11*c33-c12*c21-c13*c31-c23*c32
      a3 =  -c11*c22*c33+c11*c23*c32+c12*c21*c33+c13*c22*c31
      a3 = a3 - c12*c23*c31-c13*c21*c32

      return
      end

c***********************************************************************
      subroutine cofac(a,p,cotr,cc)
c     Generates cofactors of matrix a_jk in georg's calculations for 
c     point source problem (jmk uses different definition ?)
c                  a_jk=(a_ijkl p_i p_l - delta_jk)
c     Input: slowness vector and elasticity.
c     Output: transposed cofactor matrix with trace.
c***********************************************************************

      implicit none

      integer i,j,k,l
      real*8 cc
      real*8 a(3,3,3,3),p(3),g(3,3),co(3,3),cotr(3,3)

c     Forming gamma
c     Generating a_jk

      do 65 j=1,3
        do 55 k=1,3
          g(j,k) = 0.0d0
          do 45 i=1,3
            do 35 l=1,3
              g(j,k) = g(j,k) + p(l)*p(i)*a(i,j,k,l)
35          continue
45        continue
55      continue
65    continue

      g(1,1)=g(1,1)-1.d0
      g(2,2)=g(2,2)-1.d0
      g(3,3)=g(3,3)-1.d0

c     Generating cofactors
c     Generating b_kj and b_qq

      co(1,1)=(g(2,2)*g(3,3)-g(2,3)*g(2,3))
      co(1,2)=-(g(2,1)*g(3,3)-g(2,3)*g(3,1))
      co(1,3)=(g(2,1)*g(3,2)-g(2,2)*g(3,1))

      co(2,1)=-(g(1,2)*g(3,3)-g(1,3)*g(3,2))
      co(2,2)=(g(1,1)*g(3,3)-g(1,3)*g(3,1))
      co(2,3)=-(g(1,1)*g(3,2)-g(1,2)*g(3,1))

      co(3,1)=(g(1,2)*g(2,3)-g(1,3)*g(2,2))
      co(3,2)=-(g(1,1)*g(2,3)-g(1,3)*g(2,1))
      co(3,3)=(g(1,1)*g(2,2)-g(1,2)*g(2,1))

      cc=co(1,1)+co(2,2)+co(3,3)

      cotr(1,1)=co(1,1)
      cotr(2,2)=co(2,2)
      cotr(3,3)=co(3,3)
      cotr(1,2)=co(2,1)
      cotr(1,3)=co(3,1)
      cotr(2,1)=co(1,2)
      cotr(2,3)=co(3,2)
      cotr(3,1)=co(1,3)
      cotr(3,2)=co(2,3)

      return
      end

c***********************************************************************
      subroutine xrdvec(fsm,ivec,cdv)
c***********************************************************************

      implicit none

      integer i,j,i1,j1,i2,j2,imx,jmx,ibranch,ivec
      real*8 tst0,test,scale
      real*8 fsm(3,3),cofsm(3,3),cdv(3)

      ibranch=ivec
      tst0=0.d0

      do 100 i=1,3
        do 200 j=1,3
          i1=i+1
          if(i1.gt.3)i1=i1-3
          j1=j+1
          if(j1.gt.3)j1=j1-3
          i2=i+2
          if(i2.gt.3)i2=i2-3
          j2=j+2
          if(j2.gt.3)j2=j2-3
          cofsm(i,j)=fsm(i1,j1)*fsm(i2,j2)-fsm(i1,j2)*fsm(i2,j1)
          test=dabs(cofsm(i,j))
          if(test.gt.tst0)then
            imx=i
            jmx=j
            tst0=test
          endif
200     continue
100   continue

      scale=dabs(cofsm(imx,jmx))
c     scale=scale*scale
      cdv(jmx)=cofsm(imx,jmx)/scale
      jmx=jmx+1
      if(jmx.gt.3)jmx=jmx-3
      cdv(jmx)=cofsm(imx,jmx)/scale
      jmx=jmx+1
      if(jmx.gt.3)jmx=jmx-3
      cdv(jmx)=cofsm(imx,jmx)/scale

      return
      end

c***********************************************************************
      subroutine dpoly(xcof,cof,m,rootr,rooti,ier)
c     Extract roots of a polynomial (double precision).
c     xcof: coeffizienten des polynoms, maximal 37
c           dcoef(1) ist das konstante Glied
c     cof:hilfsvector, der nicht benoetigt wird, aber deklarieren
c     m:   grad des polynoms, maximal 36
c     rootr: vector mit realteilen der wurzeln
c     rooti: vector mit imaginaerteilen der Wurzeln
c     ier:   fehlerangabe
c***********************************************************************

      implicit none

      integer ifit,n,m,ier,nx,nxx,n2,kj1,l,mt,in,ict,i,itemp
      real*8 xo,yo,x,y,dx,dy,xpr,ypr,ux,uy,v,yt,xt,u,temp,xt2,yt2,fi
      real*8 sumsq,alpha
      real*8 cof(37),rooti(36),rootr(36),xcof(37)

c     Check for input data errors
      ifit=0
      n=m
      ier=0
      if(xcof(n+1).eq.0.) then
        ier=4
      elseif(n.le.0) then
        ier=1
      elseif(n.gt.36) then
        ier=2
      else

c     No input data errors, so begin
      nx=n
      nxx=n+1
      n2=1
      kj1=n+1
      do 10 l=1,kj1
        mt=kj1-l+1
        cof(mt)=xcof(l)
10    continue

c     Set initial values
20    xo=-.00500101d0
      yo=.01000101d0

c     Zero initial value counter
      in=0
30    x=xo

c     Increment initial values and counter
      xo=-10.d0*yo
      yo=-10.d0*x

c     Set x and y to current value
      x=xo
      y=yo
      in=in+1
      go to 50
40    ifit=1
      xpr=x
      ypr=y

c     Evaluate polynomial and derivatives
50    ict=0
60    ux=0.d0
      uy=0.d0
      v =0.d0
      yt=0.d0
      xt=1.d0
      u=cof(n+1)
      if(u.eq.0.d0) go to 120
      do 70 i=1,n
        l=n-i+1
        temp=cof(l)
        xt2=x*xt-y*yt
        yt2=x*yt+y*xt
        u=u+temp*xt2
        v=v+temp*yt2
        fi=i
        ux=ux+fi*xt*temp
        uy=uy-fi*yt*temp
        xt=xt2
        yt=yt2
70    continue
      sumsq=ux*ux+uy*uy
      if(sumsq.eq.0.d0) go to 100
      dx=(v*uy-u*ux)/sumsq
      x=x+dx
      dy=-(u*uy+v*ux)/sumsq
      y=y+dy
      if((dabs(dy)+dabs(dx)).lt.1.d-12) go to 80

c     Step iteration counter
      ict=ict+1
      if(ict.lt.500) go to 60
      if(ifit.ne.0) go to 80
      if(in.eq.5) then
        ier=3
        go to 170
      else
        go to 30
      endif
80    do 90 l=1,nxx
        mt=kj1-l+1
        temp=xcof(mt)
        xcof(mt)=cof(l)
        cof(l)=temp
90    continue
      itemp=n
      n=nx
      nx=itemp
      if(ifit.eq.0) then
        go to 40
      else
        go to 110
      endif
100   if(ifit.eq.0) go to 30
      x=xpr
      y=ypr
110   ifit=0
      if(dabs(y).lt.1.d-10*dabs(x)) go to 130
      alpha=x+x
      sumsq=x*x+y*y
      n=n-2
      go to 140
120   x=0.d0
      nx=nx-1
      nxx=nxx-1
130   y=0.d0
      sumsq=0.d0
      alpha=x
      n=n-1
140   cof(2)=cof(2)+alpha*cof(1)
      do 150 l=2,n
        cof(l+1)=cof(l+1)+alpha*cof(l)-sumsq*cof(l-1)
150   continue
160   rooti(n2)=y
      rootr(n2)=x
      n2=n2+1
      if(sumsq.ne.0.d0) then
        y=-y
        sumsq=0.d0
        go to 160
      endif
      if(n.gt.0) go to 20
      endif

170   return
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

c***********************************************************************
      subroutine crotate(alpha,beta,gam,cc,cct)
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
      real*8 alpha,beta,gam
      real*4 a(3,3)
      complex*8 csum
      complex*8 cc(3,3,3,3),cct(3,3,3,3)

      a(1,1) = sngl( dcos(alpha)*dcos(beta))
      a(1,2) = sngl( dsin(beta))
      a(1,3) = sngl(-dsin(alpha)*dcos(beta))
      a(2,1) = sngl(-dcos(gam)*dsin(beta)*dcos(alpha)+
     +     dsin(gam)*dsin(alpha))
      a(2,2) = sngl( dcos(gam)*dcos(beta))
      a(2,3) = sngl( dcos(gam)*dsin(beta)*dsin(alpha)+
     +     dsin(gam)*dcos(alpha))
      a(3,1) = sngl( dsin(gam)*dsin(beta)*dcos(alpha)+
     +     dcos(gam)*dsin(alpha))
      a(3,2) = sngl(-dsin(gam)*dcos(beta))
      a(3,3) = sngl(-dsin(gam)*dsin(beta)*dsin(alpha)+
     +     dcos(gam)*dcos(alpha))

c     c_ijkl ---> c_mnrs

      do m=1,3
        do n=1,3
          do r=1,3
            do s=1,3
              csum = cmplx(0.0,0.0)
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
              if(cabs(csum).lt.1.0) csum = cmplx(0.0,0.0)
              cct(m,n,r,s) = csum
           enddo
         enddo
       enddo
      enddo

      return
      end

c***********************************************************************
      subroutine rotvec(a,b)
c     This subroutine performs a transformation of the vector from
c     one coordinate system to another defined by the transformation
c     matrix Aij.  Aij consists of the local coordinate frame written
c     in terms of the global coordinates (each row consists of a
c     basis vector - first row is wavefront normal .. normalized 
c     slowness).  Note: transformation is from global to local frame.
c***********************************************************************

      implicit none

      integer i,j
      real*8 a(3,3),b(3),c(3)

c     A(X) --> a(x)  (global to local frame)

      do i=1,3
         c(i)=0.d0
         do j=1,3
            c(i)=c(i)+a(i,j)*b(j)
         enddo
      enddo

      do i=1,3
         b(i)=c(i)
      enddo

      return
      end

c***********************************************************************
      subroutine rotatelc(a,cc,cct)
c     This subroutine performs a transformation of the tensor c_ijkl of 
c     elastic constants based on the transformation matrix Aij.  Aij
c     consists of the local coordinate frame written in terms of the 
c     global coordinates (each row consists of a basis vector - first
c     row is wavefront normal .. normalized slowness).
c     Note: transformation is from global to local frame.
c***********************************************************************

      implicit none

      integer i,j,k,l,m,n,r,s
      real*8 csum
      real*8 cc(3,3,3,3),cct(3,3,3,3),a(3,3)

c     c(X) (c_ijkl) --> c(x) (c_mnrs)

      do 65 m=1,3
        do 60 n=1,3
          do 55 r=1,3
            do 50 s=1,3
              csum = 0.0
              do 45 i=1,3
                do 40 j=1,3
                  do 35 k=1,3
                    do 30 l=1,3
                      csum = csum + a(m,i)*a(n,j)*a(r,k)*a(s,l)*
     +                      cc(i,j,k,l)
30                  continue
35                continue
40              continue
45            continue
              cct(m,n,r,s) = csum
50          continue
55        continue
60      continue
65    continue

      return
      end

c***********************************************************************
      subroutine rinv3(a,ainv)
c     The subroutine inverts the 3x3 matrix a.
c***********************************************************************

      implicit none

      integer i,j
      real*8 deta,detb,detc,rnull
      real*8 a(3,3),ainv(3,3),co(3,3)

      rnull=0.d0

      co(1,1)=(a(2,2)*a(3,3)-a(2,3)*a(3,2))
      co(1,2)=-(a(2,1)*a(3,3)-a(2,3)*a(3,1))
      co(1,3)=(a(2,1)*a(3,2)-a(2,2)*a(3,1))

      co(2,1)=-(a(1,2)*a(3,3)-a(1,3)*a(3,2))
      co(2,2)=(a(1,1)*a(3,3)-a(1,3)*a(3,1))
      co(2,3)=-(a(1,1)*a(3,2)-a(1,2)*a(3,1))

      co(3,1)=(a(1,2)*a(2,3)-a(1,3)*a(2,2))
      co(3,2)=-(a(1,1)*a(2,3)-a(1,3)*a(2,1))
      co(3,3)=(a(1,1)*a(2,2)-a(1,2)*a(2,1))

      deta=a(1,1)*co(1,1)+a(1,2)*co(1,2)+a(1,3)*co(1,3)
      detb=a(2,1)*co(2,1)+a(2,2)*co(2,2)+a(2,3)*co(2,3)
      detc=a(3,1)*co(3,1)+a(3,2)*co(3,2)+a(3,3)*co(3,3)

c     Check to see if a is invertable

      if(deta.eq.rnull.or.detb.eq.rnull.or.detc.eq.rnull)then
         write(*,*)'Error in rinv3.  Determinant equal to zero.'
         write(*,*)'Program stopped.'
         stop
      endif

      do 10 i=1,3
        do 11 j=1,3
          ainv(i,j)=co(j,i)/deta
11      continue
10    continue

      return
      end

c***********************************************************************
      subroutine cinv3(a,ainv)
c     The subroutine inverts the complex 3x3 matrix a.
c***********************************************************************

      implicit none

      integer i,j
      complex*8 deta,detb,detc,cnull
      complex*8 a(3,3),ainv(3,3),co(3,3)

      cnull=cmplx(0.0,0.0)

      co(1,1)=(a(2,2)*a(3,3)-a(2,3)*a(3,2))
      co(1,2)=-(a(2,1)*a(3,3)-a(2,3)*a(3,1))
      co(1,3)=(a(2,1)*a(3,2)-a(2,2)*a(3,1))

      co(2,1)=-(a(1,2)*a(3,3)-a(1,3)*a(3,2))
      co(2,2)=(a(1,1)*a(3,3)-a(1,3)*a(3,1))
      co(2,3)=-(a(1,1)*a(3,2)-a(1,2)*a(3,1))

      co(3,1)=(a(1,2)*a(2,3)-a(1,3)*a(2,2))
      co(3,2)=-(a(1,1)*a(2,3)-a(1,3)*a(2,1))
      co(3,3)=(a(1,1)*a(2,2)-a(1,2)*a(2,1))

      deta=a(1,1)*co(1,1)+a(1,2)*co(1,2)+a(1,3)*co(1,3)
      detb=a(2,1)*co(2,1)+a(2,2)*co(2,2)+a(2,3)*co(2,3)
      detc=a(3,1)*co(3,1)+a(3,2)*co(3,2)+a(3,3)*co(3,3)

c     Check to see if a is invertable

      if(deta.eq.cnull.or.detb.eq.cnull.or.detc.eq.cnull)then
         write(*,*)'Error in rinv3.  Determinant equal to zero.'
         write(*,*)'Program stopped.'
         stop
      endif

      do 10 i=1,3
        do 11 j=1,3
          ainv(i,j)=co(j,i)/deta
11      continue
10    continue

      return
      end

c***********************************************************************
      subroutine rinv3a(a,ainv,itag)
c     The subroutine inverts the 3x3 matrix a.
c***********************************************************************

      implicit none

      integer i,j,itag
      real*8 deta,detb,detc
      real*8 a(3,3),ainv(3,3),co(3,3)

      co(1,1)=(a(2,2)*a(3,3)-a(2,3)*a(3,2))
      co(1,2)=-(a(2,1)*a(3,3)-a(2,3)*a(3,1))
      co(1,3)=(a(2,1)*a(3,2)-a(2,2)*a(3,1))

      co(2,1)=-(a(1,2)*a(3,3)-a(1,3)*a(3,2))
      co(2,2)=(a(1,1)*a(3,3)-a(1,3)*a(3,1))
      co(2,3)=-(a(1,1)*a(3,2)-a(1,2)*a(3,1))

      co(3,1)=(a(1,2)*a(2,3)-a(1,3)*a(2,2))
      co(3,2)=-(a(1,1)*a(2,3)-a(1,3)*a(2,1))
      co(3,3)=(a(1,1)*a(2,2)-a(1,2)*a(2,1))

      deta=a(1,1)*co(1,1)+a(1,2)*co(1,2)+a(1,3)*co(1,3)
      detb=a(2,1)*co(2,1)+a(2,2)*co(2,2)+a(2,3)*co(2,3)
      detc=a(3,1)*co(3,1)+a(3,2)*co(3,2)+a(3,3)*co(3,3)

c     Check to see if a is invertable

      if(deta.eq.0.d0.or.detb.eq.0.d0.or.detc.eq.0.d0)then
c         write(*,*)'Error in rinv3.  Eigenvectors not independent.'
         itag=1
         return
      endif

      do 10 i=1,3
        do 11 j=1,3
          ainv(i,j)=co(j,i)/deta
11      continue
10    continue

      return
      end

c***********************************************************************
      subroutine crinv3(a,ainv)
c     The subroutine inverts the 3x3 complex double precision matrix a.
c***********************************************************************

      implicit none

      integer i,j
      real*8 temp1,temp2,temp3
      complex*16 deta,detb,detc
      complex*16 a(3,3),ainv(3,3),co(3,3)

      co(1,1)=(a(2,2)*a(3,3)-a(2,3)*a(3,2))
      co(1,2)=-(a(2,1)*a(3,3)-a(2,3)*a(3,1))
      co(1,3)=(a(2,1)*a(3,2)-a(2,2)*a(3,1))

      co(2,1)=-(a(1,2)*a(3,3)-a(1,3)*a(3,2))
      co(2,2)=(a(1,1)*a(3,3)-a(1,3)*a(3,1))
      co(2,3)=-(a(1,1)*a(3,2)-a(1,2)*a(3,1))

      co(3,1)=(a(1,2)*a(2,3)-a(1,3)*a(2,2))
      co(3,2)=-(a(1,1)*a(2,3)-a(1,3)*a(2,1))
      co(3,3)=(a(1,1)*a(2,2)-a(1,2)*a(2,1))

      deta=a(1,1)*co(1,1)+a(1,2)*co(1,2)+a(1,3)*co(1,3)
      detb=a(2,1)*co(2,1)+a(2,2)*co(2,2)+a(2,3)*co(2,3)
      detc=a(3,1)*co(3,1)+a(3,2)*co(3,2)+a(3,3)*co(3,3)

      temp1=abs(realpart(deta))
      temp2=abs(realpart(detb))
      temp3=abs(realpart(detc))

c     Check to see if a is invertable.

      if(temp1.eq.0.d0.or.temp2.eq.0.d0.or.temp3.eq.0.d0)then
         write(*,*)'Error in cinv3.  Determinant equal to zero.'
         write(*,*)'Program stopped.'
         stop
      endif

      do 10 i=1,3
        do 11 j=1,3
          ainv(i,j)=co(j,i)/deta
11      continue
10    continue

      return
      end

c***********************************************************************
      subroutine trans3(a,atrans)
c     The subroutine transposes the 3x3 matrix a.
c***********************************************************************

      implicit none

      integer i,j
      real*8 a(3,3),atrans(3,3)

      do 10 i=1,3
        do 11 j=1,3
          atrans(i,j)=a(j,i)
11      continue
10    continue

      return
      end

c***********************************************************************
      subroutine ctrans3(a,atrans)
c     The subroutine transposes the complex 3x3 matrix a.
c***********************************************************************

      implicit none

      integer i,j
      complex*8 a(3,3),atrans(3,3)

      do 10 i=1,3
        do 11 j=1,3
          atrans(i,j)=a(j,i)
11      continue
10    continue

      return
      end

c***********************************************************************
      subroutine rotmat1(a,b)
c     This subroutine evaluates b(j,i)a(k,l)b(m,n) and returns the 
c     matrix product as b.  This is done to evaluate the propagator 
c     matrices RT.P.R (see notes).
c***********************************************************************

      implicit none

      integer i,j,m,n
      real*8 a(3,3),b(3,3),c(3,3)

      do i=1,3
         do j=1,3
            c(i,j)=0.d0
            do m=1,3
               do n=1,3
                  c(i,j)=c(i,j)+a(m,i)*a(n,j)*b(m,n)
               enddo
            enddo
         enddo
      enddo

      do i=1,3
         do j=1,3
            b(i,j)=c(i,j)
         enddo
      enddo

      return
      end

c***********************************************************************
      subroutine rmult3a(a,b,c)
c     This subroutine multiplies two real 3x3 matrices.
c***********************************************************************

      implicit none

      integer i,j,k
      real*8 a(3,3),b(3,3),c(3,3)

      do 10 i=1,3
        do 11 j=1,3
          c(i,j) = 0.d0
          do 12 k=1,3
            c(i,j) = c(i,j) + a(i,k)*b(k,j)
12        continue
11      continue
10    continue

      return
      end

c***********************************************************************
      subroutine rmult3b(a,b,c,d)
c     This subroutine multiplies three real 3x3 matrices.
c***********************************************************************

      implicit none

      integer i,j,k
      real*8 a(3,3),b(3,3),c(3,3),d(3,3),work1(3,3)

      do 10 i=1,3
        do 11 j=1,3
          work1(i,j) = 0.d0
          do 12 k=1,3
            work1(i,j) = work1(i,j) + a(i,k)*b(k,j)
12        continue
11      continue
10    continue

      do 15 i=1,3
        do 16 j=1,3
          d(i,j) = 0.d0
          do 17 k=1,3
            d(i,j) = d(i,j) + work1(i,k)*c(k,j)
17        continue
16      continue
15    continue

      return
      end

c***********************************************************************
      subroutine cmult3b(a,b,c,d)
c     This subroutine multiplies three complex 3x3 matrices.
c***********************************************************************

      implicit none

      integer i,j,k
      complex*8 a(3,3),b(3,3),c(3,3),d(3,3),work1(3,3)

      do 10 i=1,3
        do 11 j=1,3
          work1(i,j) = cmplx(0.0,0.0)
          do 12 k=1,3
            work1(i,j) = work1(i,j) + a(i,k)*b(k,j)
12        continue
11      continue
10    continue

      do 15 i=1,3
        do 16 j=1,3
          d(i,j) = cmplx(0.0,0.0)
          do 17 k=1,3
            d(i,j) = d(i,j) + work1(i,k)*c(k,j)
17        continue
16      continue
15    continue

      return
      end

c***********************************************************************
      subroutine rmult3c(a,b,c,d,e)
c     This subroutine multiplies four real 3x3 matrices.
c***********************************************************************

      implicit none

      integer i,j,k
      real*8 a(3,3),b(3,3),c(3,3),d(3,3),e(3,3),work1(3,3),
     +       work2(3,3)

      do 10 i=1,3
        do 11 j=1,3
          work1(i,j) = 0.d0
          do 12 k=1,3
            work1(i,j) = work1(i,j) + a(i,k)*b(k,j)
12        continue
11      continue
10    continue

      do 13 i=1,3
        do 14 j=1,3
          work2(i,j) = 0.d0
          do 15 k=1,3
            work2(i,j) = work2(i,j) + c(i,k)*d(k,j)
15        continue
14      continue
13    continue

      do 16 i=1,3
        do 17 j=1,3
          e(i,j) = 0.d0
          do 18 k=1,3
            e(i,j) = e(i,j) + work1(i,k)*work2(k,j)
18        continue
17      continue
16    continue

      return
      end

c***********************************************************************
      subroutine cmult3c(a,b,c,d,e)
c     This subroutine multiplies four complex 3x3 matrices.
c***********************************************************************

      implicit none

      integer i,j,k
      complex*8 a(3,3),b(3,3),c(3,3),d(3,3),e(3,3),work1(3,3),
     +       work2(3,3)

      do 10 i=1,3
        do 11 j=1,3
          work1(i,j) = (0.0,0.0)
          do 12 k=1,3
            work1(i,j) = work1(i,j) + a(i,k)*b(k,j)
12        continue
11      continue
10    continue

      do 13 i=1,3
        do 14 j=1,3
          work2(i,j) = cmplx(0.0,0.0)
          do 15 k=1,3
            work2(i,j) = work2(i,j) + c(i,k)*d(k,j)
15        continue
14      continue
13    continue

      do 16 i=1,3
        do 17 j=1,3
          e(i,j) = cmplx(0.0,0.0)
          do 18 k=1,3
            e(i,j) = e(i,j) + work1(i,k)*work2(k,j)
18        continue
17      continue
16    continue

      return
      end

c***********************************************************************
      subroutine rmult3d(a,b,c,d,e,f)
c     This subroutine multiplies five real 3x3 matrices.
c***********************************************************************

      implicit none

      integer i,j,k
      real*8 a(3,3),b(3,3),c(3,3),d(3,3),e(3,3),f(3,3),work1(3,3),
     +       work2(3,3),work3(3,3)

      do 10 i=1,3
        do 11 j=1,3
          work1(i,j) = 0.d0
          do 12 k=1,3
            work1(i,j) = work1(i,j) + a(i,k)*b(k,j)
12        continue
11      continue
10    continue

      do 13 i=1,3
        do 14 j=1,3
          work2(i,j) = 0.d0
          do 15 k=1,3
            work2(i,j) = work2(i,j) + c(i,k)*d(k,j)
15        continue
14      continue
13    continue

      do 16 i=1,3
        do 17 j=1,3
          work3(i,j) = 0.d0
          do 18 k=1,3
            work3(i,j) = work3(i,j) + work1(i,k)*work2(k,j)
18        continue
17      continue
16    continue

      do 19 i=1,3
        do 20 j=1,3
          f(i,j) = 0.d0
          do 21 k=1,3
            f(i,j) = f(i,j) + work3(i,k)*e(k,j)
21        continue
20      continue
19    continue

      return
      end

c***********************************************************************
      subroutine cmult3d(a,b,c,d,e,f)
c     This subroutine multiplies five complex 3x3 matrices.
c***********************************************************************

      implicit none

      integer i,j,k
      complex*8 a(3,3),b(3,3),c(3,3),d(3,3),e(3,3),f(3,3),work1(3,3),
     +       work2(3,3),work3(3,3)

      do 10 i=1,3
        do 11 j=1,3
          work1(i,j) = cmplx(0.0,0.0)
          do 12 k=1,3
            work1(i,j) = work1(i,j) + a(i,k)*b(k,j)
12        continue
11      continue
10    continue

      do 13 i=1,3
        do 14 j=1,3
          work2(i,j) = cmplx(0.0,0.0)
          do 15 k=1,3
            work2(i,j) = work2(i,j) + c(i,k)*d(k,j)
15        continue
14      continue
13    continue

      do 16 i=1,3
        do 17 j=1,3
          work3(i,j) =cmplx(0.0,0.0)
          do 18 k=1,3
            work3(i,j) = work3(i,j) + work1(i,k)*work2(k,j)
18        continue
17      continue
16    continue

      do 19 i=1,3
        do 20 j=1,3
          f(i,j) = cmplx(0.0,0.0)
          do 21 k=1,3
            f(i,j) = f(i,j) + work3(i,k)*e(k,j)
21        continue
20      continue
19    continue

      return
      end

c***********************************************************************
      subroutine cmult3a(a,b,c)
c     This subroutine multiplies two complex 3x3 matrices.
c***********************************************************************

      implicit none

      integer i,j,k
      complex*8 a(3,3),b(3,3),c(3,3)

      do 10 i=1,3
        do 11 j=1,3
          c(i,j) = cmplx(0.0,0.0)
          do 12 k=1,3
            c(i,j) = c(i,j) + a(i,k)*b(k,j)
12        continue
11      continue
10    continue

      return
      end

c***********************************************************************
      subroutine uvecshft(uvec,nom,nx3,nx2)
c     Save the vector field of the newest x1 plane into the so-called
c     z-dz buffer and moves the z-dz buffer into the z-2dz buffer. 
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer nom,nx2,nx3,ix2,ix3,it,i,nantest,istop
      complex*16 uvec(3,ntmx,nx3mx,nx2mx,3)

      istop=0
      nantest=1
      do 10 ix2=1,nx2
        do 20 ix3=1,nx3
          do 30 it=1,nom
            do 40 i=1,3
              uvec(i,it,ix3,ix2,1)=uvec(i,it,ix3,ix2,2)
              uvec(i,it,ix3,ix2,2)=uvec(i,it,ix3,ix2,3)
              uvec(i,it,ix3,ix2,3)=dcmplx(0.d0,0.d0)
              if(nantest.eq.1)then
                 if(isnan(realpart(uvec(i,it,ix3,ix2,1))).or.
     +                isnan(imagpart(uvec(i,it,ix3,ix2,1))))then
                    istop=1
                 elseif(isnan(realpart(uvec(i,it,ix3,ix2,2))).or.
     +                   isnan(imagpart(uvec(i,it,ix3,ix2,2))))then
                    istop=1
                 elseif(isnan(realpart(uvec(i,it,ix3,ix2,3))).or.
     +                   isnan(imagpart(uvec(i,it,ix3,ix2,3))))then
                    istop=1
                 endif
              endif
              if(istop.eq.1)then
                 write(*,*)'Problem: NaN in uvec wavefield'
                 write(*,*)'i,it,ix2,ix3: ',i,it,ix2,ix3
                 write(*,*)'uvec(i,1): ',uvec(i,it,ix3,ix2,1)
                 write(*,*)'uvec(i,2): ',uvec(i,it,ix3,ix2,2)
                 write(*,*)'uvec(i,3): ',uvec(i,it,ix3,ix2,3)
                 stop
              endif
40          continue
30        continue
20      continue
10    continue

      return
      end

c***********************************************************************
      SUBROUTINE FOUR1(DATA,NN,ISIGN)
c     Replaces DATA by its discrete Fourier transform, if ISIGN is
c     input as 1; or replaces DATA by NN times its inverse discrete
c     Fourier transform, if ISIGN is input as -1.  DATA is a complex
c     array of length NN or, equivalently, a real array of length 2*NN.
c     NN must be an integer power of 2 (this is not checked for!).
c***********************************************************************
 
      implicit none

      INTEGER NN,ISIGN,I,J,N,M,MMAX,ISTEP
      REAL*8 TEMPR,TEMPI,THETA,WPR,WPI,WR,WI,WTEMP
      REAL*8 DATA(2*NN)

      N=2*NN
      J=1
      DO 11 I=1,N,2
        IF(J.GT.I)THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        ENDIF
        M=N/2
1       IF ((M.GE.2).AND.(J.GT.M)) THEN
          J=J-M
          M=M/2
        GO TO 1
        ENDIF
        J=J+M
11    CONTINUE
      MMAX=2
2     IF (N.GT.MMAX) THEN
        ISTEP=2*MMAX
        THETA=6.28318530717959D0/(ISIGN*MMAX)
        WPR=-2.D0*DSIN(0.5D0*THETA)**2
        WPI=DSIN(THETA)
        WR=1.D0
        WI=0.D0
        DO 13 M=1,MMAX,2
          DO 12 I=M,N,ISTEP
            J=I+MMAX
            TEMPR=WR*DATA(J)-WI*DATA(J+1)
            TEMPI=WR*DATA(J+1)+WI*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
12        CONTINUE
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
13      CONTINUE
        MMAX=ISTEP
      GO TO 2
      ENDIF
      RETURN
      END

c***********************************************************************
      subroutine rtaper(nom2,iomc1,iomc2,rfilt)
c     This subroutine applies a high pass filter to the frequency data
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer iom,nom2,iomc1,iomc2
      real*8 pi,rat,rat1
      real*8 rfilt(ntmx)

      if(iomc2.gt.nom2)write(11,*)'ERROR in rtaper'

      pi=4.d0*datan(1.d0)
  
      rat1=dble(iomc2-iomc1)/pi

      do 100 iom=1,nom2
        if(iomc1.eq.iomc2)then
          rfilt(iom)=1.d0
        else
          if(iom.lt.iomc1)then
            rfilt(iom)=0.d0
          elseif(iom.ge.iomc1.and.iom.le.iomc2)then
            rat=dble(iom-iomc1)/rat1
            rfilt(iom)=(1.d0-dcos(rat))/2.d0
          else
            rfilt(iom)=1.d0
          endif
        endif
100   continue

      return
      end

c***********************************************************************
      subroutine detm(a,d)
c     The subroutine inverts the 3x3 matrix a.
c***********************************************************************

      implicit none

      real*8 deta,detb,detc,d
      real*8 a(3,3),co(3,3)

      co(1,1)=(a(2,2)*a(3,3)-a(2,3)*a(3,2))
      co(1,2)=-(a(2,1)*a(3,3)-a(2,3)*a(3,1))
      co(1,3)=(a(2,1)*a(3,2)-a(2,2)*a(3,1))

      co(2,1)=-(a(1,2)*a(3,3)-a(1,3)*a(3,2))
      co(2,2)=(a(1,1)*a(3,3)-a(1,3)*a(3,1))
      co(2,3)=-(a(1,1)*a(3,2)-a(1,2)*a(3,1))

      co(3,1)=(a(1,2)*a(2,3)-a(1,3)*a(2,2))
      co(3,2)=-(a(1,1)*a(2,3)-a(1,3)*a(2,1))
      co(3,3)=(a(1,1)*a(2,2)-a(1,2)*a(2,1))

      deta=a(1,1)*co(1,1)+a(1,2)*co(1,2)+a(1,3)*co(1,3)
      detb=a(2,1)*co(2,1)+a(2,2)*co(2,2)+a(2,3)*co(2,3)
      detc=a(3,1)*co(3,1)+a(3,2)*co(3,2)+a(3,3)*co(3,3)

c     Check to see if a is invertable

      if(deta.eq.0.d0.or.detb.eq.0.d0.or.detc.eq.0.d0)then
         write(*,*)'Error in detm.  Determinant equal to zero.'
         write(*,*)'Program stopped.'
         stop
      endif

      d=(deta+detb+detc)/3.d0

      return
      end
