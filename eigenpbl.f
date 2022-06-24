c........+.........+.........+.........+.........+.........+.........+..
c  File Eigenpbl.f contains the following 12 subroutines:  
c  Slowness, Anisoev, Isoev, Vertslw, Eig3, Polarr, Velg, Sortp,
c  Preturn, Tql2 and Tred2.
c
c  Last modified: April 17, 2012.
c***********************************************************************
      subroutine slowness(cc,phi,theta,vn,p,v)
c     Input: cc        - elastic constants devided by density
c            phi/theta - take off angles
c            vn        - normal, phase velocity
c            p         - slowness
c            v         - group velocity
c     Numbering of slowness components different to Georg's.
c     wn(1) <---> z
c     wn(2) <---> x
c     wn(3) <---> y
c     This subroutine is called once to set up the "gradients" of the
c     slownesses over the grid.  Based on Gslow but modified for 
c     Oneway coordinate system.
c***********************************************************************

      implicit none

      integer i,j,k,l
      real*8 phi,theta,csum
      real*8 a1,a2,a3,x1,x2r,x2i,x3r,x3i
      real*8 bqq1,bqq2,bqq3

      real*8 cc(3,3,3,3),c(3,3)
      real*8 v(3,3),p(3,3),vn(3),wn(3)
      real*8 e1(3,3,3),e2(3,3,3),e3(3,3,3)
      real*8 f1(3),f2(3),f3(3)
      real*8 p1(3),p2(3),p3(3)
      real*8 b1(3,3),b2(3,3),b3(3,3)

      wn(1) = dcos(theta)
      wn(2) = dsin(theta)*dcos(phi)
      wn(3) = dsin(theta)*dsin(phi)

      do 13 i=1,3
         do 12 l=1,3
            csum = 0.0d0
            do 11 j=1,3
               do 10 k=1,3
                  csum = csum + cc(i,j,k,l)*wn(j)*wn(k)
 10            continue
 11            continue
               c(i,l) = csum
               c(l,i) = csum
 12         continue
 13      continue

      call coeff(c,a1,a2,a3)

      call roots(a1,a2,a3,x1,x2r,x2i,x3r,x3i)

c     Normal velocity for three waves

      vn(1) = dsqrt(x1)
      vn(2) = dsqrt(x2r)
      vn(3) = dsqrt(x3r)

c     Slowness components for three quasi ?-waves

      do 15 i=1,3
         p(1,i) =wn(i)/vn(1)
         p1(i)=p(1,i)
 15   continue
      do 16 i=1,3
         p(2,i) =wn(i)/vn(2)
         p2(i)=p(2,i)
 16   continue
      do 17 i=1,3
         p(3,i) =wn(i)/vn(3)
         p3(i)=p(3,i)
 17   continue

c     Calculate group velocity

      do 20 i=1,3
         do 18 j=1,3
            do 19 k=1,3
               e1(i,j,k)=0.d0
               e2(i,j,k)=0.d0
               e3(i,j,k)=0.d0
               do 21 l=1,3
                  e1(i,j,k) = e1(i,j,k)+cc(i,j,k,l)*p(1,l)
                  e2(i,j,k) = e2(i,j,k)+cc(i,j,k,l)*p(2,l)
                  e3(i,j,k) = e3(i,j,k)+cc(i,j,k,l)*p(3,l)
 21            continue
 19         continue
 18      continue
 20   continue

c     Get transposed cofactors

      call cofac(cc,p1,b1,bqq1)
      call cofac(cc,p2,b2,bqq2)
      call cofac(cc,p3,b3,bqq3)

      do 23 i=1,3
         f1(i)=0.d0
         f2(i)=0.d0
         f3(i)=0.d0
         do 24 j=1,3
            do 25 k=1,3
               f1(i) = f1(i)+e1(i,j,k)*b1(j,k)
               f2(i) = f2(i)+e2(i,j,k)*b2(j,k)
               f3(i) = f3(i)+e3(i,j,k)*b3(j,k)
 25         continue
 24      continue
 23   continue
      
      do 28 j=1,3
         v(1,j)=f1(j)/bqq1
 28   continue
      do 29 j=1,3
         v(2,j)=f2(j)/bqq2
 29   continue
      do 31 j=1,3
         v(3,j)=f3(j)/bqq3
 31   continue

      return
      end

c***********************************************************************
      subroutine anisoev(a,p2,p3,eval,evec,iflag,nreal,k2vec)
c     This subroutine is similar to routine Evprob from Georg.
c     Only real eigenvalues and displacements are of interest and 
c     stresses are not required (although they should be easy to 
c     include).  Eval and evec are real here, and p1r and pol are used 
c     as temporary storage for dealing with lower level subroutines 
c     verstlw and polarr.
c     Only the forward-going (positive) slownesses and vectors are 
c     needed, and they are sorted in order of magnitude (so P is first).
c     Note: When P is evanescent only S1 and S2 will be used to form the
c     propagator and only S2 when S1 goes evanescent.
c     Modifications for Oneway coordinates.
c***********************************************************************

      implicit none

      integer nreal,imax,imid,imin,i,j,ii,kk,jflag
      integer k2vec(2),iflag(6)
      real*8 p2,p3,temp,temp2,temp3
      real*8 a(3,3,3,3)
      real*8 eval(6),evec(6,6)
      real*8 pvec(3),pol(3)
      complex*16 p1(6)

      call vertslw(a,p2,p3,p1,iflag,nreal)

      if(nreal.eq.3)then

c     Three real slowness will contribute to the propagator --
c     first sort them in order of size.

        imax=1
        imin=1
        temp=dble(p1(1))
        temp2=dble(p1(1))
        do 99 ii=2,3
           temp3=dble(p1(ii))
           if(temp3.gt.temp)then
              imax=ii
              temp=temp3
           endif
           if(temp3.lt.temp2)then
              imin=ii
              temp2=temp3
           endif
 99     continue

        imid=6-imax-imin

        if(imin.ne.1.or.imid.ne.2.or.imax.ne.3)then
           write(*,*)'In anisoev,imin,imid,imax=',imin,imid,imax
        elseif(imin.eq.imax.or.imin.eq.imid.or.imid.eq.imax)then
           write(*,*)'** Error-1 in anisoev **'
           stop
        endif
        
        do 100 ii=1,3
           if(ii.eq.1)i=imin
           if(ii.eq.2)i=imid
           if(ii.eq.3)i=imax
           eval(ii)=dble(p1(i))
           pvec(2)=p2
           pvec(3)=p3
           pvec(1)=eval(ii)
           call polarr(a,pvec,pol,i)
           do 200 j=1,3
              evec(j,ii)=pol(j)
 200       continue
 100    continue
        
      elseif(nreal.eq.2)then
         
c     Two real slownesses and their vectors contribute -- these are 
c     stored in the first two elements of eval and first two columns of 
c     evec.
         
         kk=0
         do 300 jflag=1,3
            if(iflag(jflag).eq.1)then
               kk=kk+1
               if(kk.eq.3) write(*,*)'** Error-2 in anisoev, kk=3 **',
     +              ' nreal=',nreal,'  iflag(1-3)=',iflag(1),iflag(2),
     +              iflag(3)
               k2vec(kk)=jflag
            endif
 300     continue
         
         do 400 ii=1,2
            eval(ii)=dble(p1(k2vec(ii)))
            pvec(2)=p2
            pvec(3)=p3
            pvec(1)=eval(ii)  
            call polarr(a,pvec,pol,k2vec(ii))
            
            do 500 j=1,3
               evec(j,ii)=pol(j)
 500        continue
 400     continue
         
      elseif(nreal.eq.1)then

c     One real slowness and vector contribute.

         k2vec(1)=0
         do 600 jflag=1,3
            if(iflag(jflag).eq.1)then
               if(k2vec(1).ne.0)then
                  write(*,*)'** Error-3 in ansisoev **'
                  stop
               endif
               k2vec(1)=jflag
            endif
 600     continue
         eval(1)=dble(p1(k2vec(1)))
         
         pvec(2)=p2
         pvec(3)=p3
         pvec(1)=eval(1)       
         call polarr(a,pvec,pol,k2vec(1))
         
         do 700 j=1,3
            evec(j,1)=pol(j)
 700     continue
         
      endif
      
      return
      end

c***********************************************************************
      subroutine isoev(a,p2,p3,eval,evec,nreal)
c     Modification of isoeig for Oneway coordinate system.
c     NOTE: c11 = lambda + 2mu  and  c66 = mu   in Lame notation.
c     Sean's comments follow:  This subroutine calculates the 
c     displacement-stress eigenvectors and the eigenvalues in the case 
c     of an isotropic medium. The coding follows Fryer and Frazer (1987)
c     with the obvious errors corrected.  The vectors are choosen to be 
c     in a right handed system. P is along the slowness vector, SV is 
c     pointed upwards, and SH is horizontal but its direction depends on
c     the value of p(x) (It goes to the left of the p(x) direction). All
c     of the vectors are normalized so that the displacements are of 
c     unit magnitude.  Orientations rechecked, 07/92 by Sean Guest.
c     Input:  matrix of elastic constants (devided by density) and two 
c             slowness components (px and py).
c     Output: Displacement and stress (all zero for stress) 
c             eigensolutions (real componnets).
c     Orientations check Aug 17, 1999 by Doug Angus (RHR).
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer nreal,itest,i,j
      real*8 pa,p2,p3,fac1,fac2,acomp,bcomp,ccomp,temp
      real*8 c(6,6),a(3,3,3,3),eval(6),evec(6,6)

      c(1,1)=a(1,1,1,1)*den
      c(6,6)=a(2,1,2,1)*den
      pa=p2**2.d0+p3**2.d0
      fac1=den/c(1,1)-pa
      fac2=den/c(6,6)-pa
      nreal=3

c     Define the eigenvalues

      eval(1)=dsqrt(fac1)
      eval(2)=dsqrt(fac2)
      eval(3)=eval(2)
      eval(4)=0.d0
      eval(5)=0.d0
      eval(6)=0.d0

c     Evaluation of the Eigenvectors
      
c     Set the stress components to zero.

      evec(1,4)=0.d0
      evec(2,4)=0.d0
      evec(3,4)=0.d0
      evec(4,4)=0.d0
      evec(5,4)=0.d0
      evec(6,4)=0.d0

      evec(1,5)=0.d0
      evec(2,5)=0.d0
      evec(3,5)=0.d0
      evec(4,5)=0.d0
      evec(5,5)=0.d0
      evec(6,5)=0.d0

      evec(1,6)=0.d0
      evec(2,6)=0.d0
      evec(3,6)=0.d0
      evec(4,6)=0.d0
      evec(5,6)=0.d0
      evec(6,6)=0.d0

c  Solution of the Eigenvectors based on the p2 and p3 values.

      if(p2.ne.0.and.p3.ne.0)then

         evec(1,1)=eval(1)
         evec(2,1)=p2
         evec(3,1)=p3
         evec(4,1)=0.d0
         evec(5,1)=0.d0
         evec(6,1)=0.d0
         
         evec(1,2)=-pa/eval(2)
         evec(2,2)=p2
         evec(3,2)=p3
         evec(4,2)=0.d0
         evec(5,2)=0.d0
         evec(6,2)=0.d0
         
         evec(1,3)=0.d0
         evec(2,3)=-p3
         evec(3,3)=p2
         evec(4,3)=0.d0
         evec(5,3)=0.d0
         evec(6,3)=0.d0
         
      elseif(p2.eq.0.and.p3.ne.0)then
         
         evec(1,1)=eval(1)
         evec(2,1)=0.d0
         evec(3,1)=p3
         evec(4,1)=0.d0
         evec(5,1)=0.d0
         evec(6,1)=0.d0

         evec(1,3)=0.d0
         evec(2,3)=-p3
         evec(3,3)=0.d0
         evec(4,3)=0.d0
         evec(5,3)=0.d0
         evec(6,3)=0.d0

         evec(1,2)=-evec(2,1)*evec(3,3)+evec(3,1)*evec(2,3)
         evec(2,2)=-evec(3,1)*evec(1,3)+evec(1,1)*evec(3,3)
         evec(3,2)=-evec(1,1)*evec(2,3)+evec(2,1)*evec(1,3)
         evec(4,2)=0.d0
         evec(5,2)=0.d0
         evec(6,2)=0.d0
         
      elseif(p2.ne.0.and.p3.eq.0)then
         
         evec(1,1)=eval(1)
         evec(2,1)=p2
         evec(3,1)=0.d0
         evec(4,1)=0.d0
         evec(5,1)=0.d0
         evec(6,1)=0.d0
       
         evec(1,3)=0.d0
         evec(2,3)=0.d0
         evec(3,3)=p2
         evec(4,3)=0.d0
         evec(5,3)=0.d0
         evec(6,3)=0.d0
         
         evec(1,2)=-evec(2,1)*evec(3,3)+evec(3,1)*evec(2,3)
         evec(2,2)=-evec(3,1)*evec(1,3)+evec(1,1)*evec(3,3)
         evec(3,2)=-evec(1,1)*evec(2,3)+evec(2,1)*evec(1,3)
         evec(4,2)=0.d0
         evec(5,2)=0.d0
         evec(6,2)=0.d0
         
      elseif(p2.eq.0.and.p3.eq.0)then
         
         evec(1,1)=eval(1)
         evec(2,1)=0.d0
         evec(3,1)=0.d0
         evec(4,1)=0.d0
         evec(5,1)=0.d0
         evec(6,1)=0.d0
         
         evec(1,2)=0.d0
         evec(2,2)=eval(2)
         evec(3,2)=0.d0
         evec(4,2)=0.d0
         evec(5,2)=0.d0
         evec(6,2)=0.d0
         
         evec(1,3)=0.d0
         evec(2,3)=0.d0
         evec(3,3)=eval(3)
         evec(4,3)=0.d0
         evec(5,3)=0.d0
         evec(6,3)=0.d0
         
      endif

c     Test for orthogonality
      itest=0
      if(itest.eq.1)then
         write(*,*)evec(1,1)*evec(1,2)+evec(2,1)*evec(2,2)+
     +        evec(3,1)*evec(3,2)
         write(*,*)evec(1,2)*evec(1,3)+evec(2,2)*evec(2,3)+
     +        evec(3,2)*evec(3,3)
         write(*,*)evec(1,1)*evec(1,3)+evec(2,1)*evec(2,3)+
     +        evec(3,1)*evec(3,3)
         read(*,*)
      endif

c     Normilise the eigendisplacements.

      do 50 j=1,3
         acomp=evec(1,j)*evec(1,j)
         bcomp=evec(2,j)*evec(2,j)
         ccomp=evec(3,j)*evec(3,j)
         temp=dsqrt(acomp+bcomp+ccomp)
         do 60 i=1,3
            evec(i,j)=evec(i,j)/temp
 60      continue
 50   continue
      
      return
      end

c***********************************************************************
      subroutine vertslw(a,p2,p3,p1,iflag,nreal)
c     Finds vertical slowness for given horizontal slowness.
c     This version sorts real and imaginary parts into positive and
c     negative groups.  It also uses the group velocity to sort up/down 
c     waves when the slowness is real, and it flags the real and complex
c     slownesses for later selection (iflag=+/-1 for real, +/-2 for 
c     complex, where +/- indicates c up or down).
c     Uses Georg's routines to solve the relevant linear systems.
c***********************************************************************
     
      implicit none

      integer ipos,ineg,nreal,l
      integer iflag(6)

      real*8 p2,p3,acc,test,cc
      real*8 a(3,3,3,3),p1r(6),p1i(6)
      real*8 cotr(3,3),p(3),v(3)
      complex*16 p1(6)

      call eigp3(a,p2,p3,p1r,p1i)

c     sort according to size of p3

      ipos=0
      ineg=0
      nreal=0
      acc=1.d-12
      do 100 l=1,6
         test=0.d0
         v(1)=0.d0
         if(p1r(l).eq.0.d0)then
            if(p1i(l).eq.0.d0)then
               write(*,*)'**Error-1 in vertslw'
               stop
            endif
            test=2.d0*acc
         else
            if(p1i(l).ne.0.d0)then
               test=dabs(p1i(l)/p1r(l))
            else
               p(2)=p2
               p(3)=p3
               p(1)=p1r(l)
               call cofac(a,p,cotr,cc)
               call velg(a,p,cotr,v)
            endif
         endif
         if(test.gt.acc)then
            if(p1i(l).gt.0.d0)then
               ipos=ipos+1
               if(ipos.gt.3) write(*,*)'**Error-ipos-1 in vertslw**'
               p1(ipos)=dcmplx(p1r(l),p1i(l))
               iflag(ipos)=2
            else
               ineg=ineg+1
               if(ineg.gt.3)write(*,*)'**Error-ineg-1 in vertslw**'
               p1(ineg+3)=dcmplx(p1r(l),p1i(l))
               iflag(ineg+3)=-2
            endif
         else
            if(v(1).gt.0.d0)then
               ipos=ipos+1
               if(ipos.gt.3)write(*,*)'**Error-ipos-2 in vertslw**'
               p1(ipos)=dcmplx(p1r(l),p1i(l))
               iflag(ipos)=1
               nreal=nreal+1
            else
               ineg=ineg+1
               if(ineg.gt.3)write(*,*)'**Error-ineg-2 in vertslw**'
               p1(ineg+3)=dcmplx(p1r(l),p1i(l))
               iflag(ineg+1)=-1
            endif
         endif
 100  continue
      
      if(ipos.ne.3.or.ineg.ne.3)then
         write(*,*)'**Error-2 in vertslw'
         read(*,*)
      endif
      if(nreal.gt.3)then
         write(*,*)'**Error-3 in vertslw, nreal.gt.3'
      endif
      
      return
      end

c***********************************************************************
      subroutine eigp3(a,p2,p3,p1r,p1i)
c     Input:  matrix of elastic constants (devided by density) and two 
c             slowness components (p1 and p2).
c     Output: six solutions for the slowness component p3 split into 
c             real and imaginary parts.
c     Experience shows that as soon as the eigenvalues are not elements 
c     of the slowness surface the imaginary part becomes unequal 0.
c     Modified from eigp3 to include the isotropic case calling the 
c     isoeig subroutine.
c***********************************************************************

      implicit none

      integer i,j,k,j1,k1,j2,k2,j3,k3,ier
      real*8 p2,p3,ax,bx,cx,dx,ex
      real*8 a(3,3,3,3)
      real*8 alpha(3,3), beta(3,3), gamma(3,3)
      real*8 axx(6),bxx(6),cxx(6),dxx(6),exx(6),fxx(6),gxx(6)
      real*8 coe(7),p1r(6),p1i(6),xcoe(7)

      do 1 j=1,3
         do 2 k=1,3
            alpha(j,k) = a(2,j,k,2)*p2*p2 + a(3,j,k,2)*p2*p3 +
     +           a(2,j,k,3)*p3*p2 + a(3,j,k,3)*p3*p3
            if(j.eq.k) alpha(j,k)=alpha(j,k)-1.d0
            beta(j,k)  = a(1,j,k,2)*p2 + a(1,j,k,3)*p3 +
     +           a(2,j,k,1)*p2 + a(3,j,k,1)*p3
            gamma(j,k) = a(1,j,k,1)
 2       continue
 1    continue
      
      do 5 i=1,6
         if(i.eq.1) then
            j1=3
            k1=3
            j2=1
            k2=1
            j3=2
            k3=2
         else if(i.eq.2) then
            j1=3
            k1=1
            j2=1
            k2=3
            j3=2
            k3=2
         else if(i.eq.3) then
            j1=3
            k1=1
            j2=1
            k2=2
            j3=2
            k3=3
         else if(i.eq.4) then
            j1=3
            k1=2
            j2=1
            k2=1
            j3=2
            k3=3
         else if(i.eq.5) then
            j1=3
            k1=2
            j2=1
            k2=3
            j3=2
            k3=1
         else
            j1=3
            k1=3
            j2=1
            k2=2
            j3=2
            k3=1
         end if

         ax = alpha(j1,k1)*alpha(j2,k2)
         bx = alpha(j1,k1)*beta(j2,k2)+beta(j1,k1)*alpha(j2,k2)
         cx = alpha(j1,k1)*gamma(j2,k2)+
     +        beta(j1,k1)*beta(j2,k2)+
     +        gamma(j1,k1)*alpha(j2,k2)
         dx = beta(j1,k1)*gamma(j2,k2)+gamma(j1,k1)*beta(j2,k2)
         ex = gamma(j1,k1)*gamma(j2,k2)

         axx(i) = alpha(j3,k3)*ax
         bxx(i) = alpha(j3,k3)*bx+beta(j3,k3)*ax
         cxx(i) = alpha(j3,k3)*cx+beta(j3,k3)*bx+gamma(j3,k3)*ax
         dxx(i) = alpha(j3,k3)*dx+beta(j3,k3)*cx+gamma(j3,k3)*bx
         exx(i) = alpha(j3,k3)*ex+beta(j3,k3)*dx+gamma(j3,k3)*cx
         fxx(i) = beta(j3,k3)*ex+gamma(j3,k3)*dx
         gxx(i) = gamma(j3,k3)*ex
         
 5    continue
      
      coe(1)=axx(1)-axx(2)+axx(3)-axx(4)+axx(5)-axx(6)
      coe(2)=bxx(1)-bxx(2)+bxx(3)-bxx(4)+bxx(5)-bxx(6)
      coe(3)=cxx(1)-cxx(2)+cxx(3)-cxx(4)+cxx(5)-cxx(6)
      coe(4)=dxx(1)-dxx(2)+dxx(3)-dxx(4)+dxx(5)-dxx(6)
      coe(5)=exx(1)-exx(2)+exx(3)-exx(4)+exx(5)-exx(6)
      coe(6)=fxx(1)-fxx(2)+fxx(3)-fxx(4)+fxx(5)-fxx(6)
      coe(7)=gxx(1)-gxx(2)+gxx(3)-gxx(4)+gxx(5)-gxx(6)

      coe(1)=coe(1)/coe(7)
      coe(2)=coe(2)/coe(7)
      coe(3)=coe(3)/coe(7)
      coe(4)=coe(4)/coe(7)
      coe(5)=coe(5)/coe(7)
      coe(6)=coe(6)/coe(7)
      coe(7)=1.d0

      call dpoly(coe,xcoe,6,p1r,p1i,ier)

      return
      end

c***********************************************************************
      subroutine polarr(a,p,pol,ivec)
c     Generates polarisation vectors of
c       a_jk=(a_ijkl p_i p_l - delta_jk) pol_k = 0.
c     Input: slowness vector and elasticity
c     Output: polarisation vectors
c     This version for REAL evectors and slownesses only.
c***********************************************************************

      implicit none

      integer ivec,i,j,k,l
      real*8 anorm
      real*8 a(3,3,3,3)
      real*8 p(3),g(3,3),pol(3)

c     Forming gamma
c     Generating g_jk

      do 65 j=1,3
         do 55 k=1,3
            g(j,k) = (0.0d0,0.d0)
            do 45 i=1,3
               do 35 l=1,3
                  g(j,k) = g(j,k) + p(l)*p(i)*a(i,j,k,l)
 35            continue
 45         continue
 55      continue
 65   continue

      g(1,1)=g(1,1)-1.d0
      g(2,2)=g(2,2)-1.d0
      g(3,3)=g(3,3)-1.d0

      call xrdvec(g,ivec,pol)

c     normalise

      anorm=dsqrt(pol(1)**2+pol(2)**2+pol(3)**2) 
      pol(1)=pol(1)/anorm
      pol(2)=pol(2)/anorm
      pol(3)=pol(3)/anorm

      return
      end

c***********************************************************************
      subroutine velg(cc,p,cotr,v)
c     Calculates group velocity for given slowness and transposed 
c     cofactor and elstic constants.
c     p: slowness
c     v: group velocity
c     v(1) --->  x_1
c     v(2) --->  x_2
c     v(3) --->  x_3
c***********************************************************************

      implicit none

      integer i,j,k,l
      real*8 tra
      real*8 v(3),p(3),cc(3,3,3,3)
      real*8 e1(3,3,3),f1(3),cotr(3,3)

      tra=cotr(1,1)+cotr(2,2)+cotr(3,3)

c     calculate group velocity

      do 20 i=1,3
         do 18 j=1,3
            do 19 k=1,3
               e1(i,j,k)=0.d0
               do 21 l=1,3
                  e1(i,j,k) = e1(i,j,k)+cc(i,j,k,l)*p(l)
 21            continue
 19         continue
 18      continue
 20   continue

      do 23 i=1,3
         f1(i)=0.d0
         do 24 j=1,3
            do 25 k=1,3
               f1(i) = f1(i)+e1(i,j,k)*cotr(j,k)
 25         continue
 24      continue
 23   continue

      do 28 j=1,3
         v(j)=f1(j)/tra
 28   continue

      return
      end

c***********************************************************************
      subroutine sortp(p3,p1,iflag,nreal)
c     This subroutine sorts the forward-going (positive) slownesses in
c     order of magnitude (P first).
c***********************************************************************

      implicit none

      integer imax,imid,imin,nreal,i,ii,kk,jflag
      integer iflag(6),k2vec(2)
      real*8 temp,temp2,temp3
      real*8 p1(3)
      complex*16 p3(6)

      if(nreal.eq.3)then        !Three real slowness components
        imax=1
        imin=1
        temp=dble(p3(1))
        temp2=dble(p3(1))
        do 99 ii=2,3
           temp3=dble(p3(ii))
           if(temp3.gt.temp)then
              imax=ii
              temp=temp3
           endif
           if(temp3.lt.temp2)then
              imin=ii
              temp2=temp3
           endif
 99     continue

        imid=6-imax-imin

        if(imin.ne.1.or.imid.ne.2.or.imax.ne.3)then
           write(*,*)'in evprob,imin,imid,imax=',imin,imid,imax
        elseif(imin.eq.imax.or.imin.eq.imid.or.imid.eq.imax)then
           write(*,*)'**Error-1 in evprob**'
           stop
        endif
        
        do 100 ii=1,3
           if(ii.eq.1)i=imin
           if(ii.eq.2)i=imid
           if(ii.eq.3)i=imax
           p1(ii)=dble(p3(i))
 100    continue
        
      elseif(nreal.eq.2)then    !Two real slowness components

         kk=0
         do 300 jflag=1,3
            if(iflag(jflag).eq.1)then
               kk=kk+1
               if(kk.eq.3) write(*,*)'**Error-2 in evprob, kk=3**',
     +              ' nreal=',nreal,'  iflag(1-3)=',iflag(1),iflag(2),
     +              iflag(3)
               k2vec(kk)=jflag
            endif
 300     continue
         
         do 400 ii=1,2
            p1(ii)=dble(p3(k2vec(ii)))
 400     continue
         
      elseif(nreal.eq.1)then    !One real slowness component

         k2vec(1)=0
         do 600 jflag=1,3
            if(iflag(jflag).eq.1)then
               if(k2vec(1).ne.0)then
                  write(*,*)'**Error-3 in evprob**'
                  stop
               endif
               k2vec(1)=jflag
            endif
 600     continue
         p1(1)=dble(p3(k2vec(1)))         
      endif
      
      return
      end

c***********************************************************************
      subroutine preturn(p1,p,nreal)
c     This subroutine returns the appropriate slownesses.
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer nreal
      real*8 p1(3),p(3,3)

      if(iwave.eq.1)then
         if(nreal.eq.3)then
            p(iwave,1)=dble(p1(1))
         else
            write(*,*)'Input problems:parameter iwave=',iwave
            write(*,*)'but nreal=',nreal,' (p-wave evanescent).'
            write(*,*)'Program stopped within Preturn.'
            stop
         endif
      elseif(iwave.eq.2)then
         if(nreal.eq.3)then
            p(iwave,1)=dble(p1(2))
         elseif(nreal.eq.2)then
            p(iwave,1)=dble(p1(1))
         else
            if(itropic.eq.0)then
               p(iwave,1)=dble(p1(1))
            else
               write(*,*)'Input problems:parameter iwave=',iwave
               write(*,*)'but nreal=',nreal,' (p-, s1-wave 
     +              evanescent)'
               write(*,*)'Program stopped within Preturn.'
               stop
            endif
         endif
      elseif(iwave.eq.3)then
         if(nreal.eq.3)then
            p(iwave,1)=dble(p1(3))
         elseif(nreal.eq.2)then
            p(iwave,1)=dble(p1(2))
         elseif(nreal.eq.1)then
            p(iwave,1)=dble(p1(1))
         else
            write(*,*)'Wave type :',iwave,' evanescent.'
            write(*,*)'Only ',nreal,' non-evanescent wave 
     +           type(s).'
            write(*,*)'Program stopped within Preturn.'
            stop
         endif
      endif

      return
      end

c***********************************************************************
      subroutine tred2(nm,n,a,d,e,z)
c     This subroutine is a translation of the algol procedure Tred2,
c     Num. Math. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
c     Handbook for Auto. Comp., Vol.ii-Linear Algebra, 212-226(1971).
c     This subroutine reduces a real symmetric matrix to a
c     symmetric tridiagonal matrix using and accumulating
c     orthogonal similarity transformations.
c     On input:
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c        n is the order of the matrix.
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c     On output:
c        d contains the diagonal elements of the tridiagonal matrix.
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c        z contains the orthogonal transformation matrix
c          produced in the reduction.
c        a and z may coincide.  if distinct, a is unaltered.
c     Questions and comments should be directed to Burton S. Garbow,
c     Mathematics and Computer Science Div, Argonne National Laboratory
c     This version dated august 1983.
c***********************************************************************

      integer i,j,k,l,n,ii,nm,jp1
      real*8 a(nm,n),d(n),e(n),z(nm,n)
      real*8 f,g,h,hh,scale

      do 100 i = 1, n

         do 80 j = i, n
   80    z(j,i) = a(j,i)

         d(i) = a(n,i)
  100 continue

      if (n .eq. 1) go to 510
c     .......... for i=n step -1 until 2 do -- ..........
      do 300 ii = 2, n
         i = n + 2 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 2) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(d(k))

         if (scale .ne. 0.0d0) go to 140
  130    e(i) = d(l)

         do 135 j = 1, l
            d(j) = z(l,j)
            z(i,j) = 0.0d0
            z(j,i) = 0.0d0
  135    continue

         go to 290

  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue

         f = d(l)
         g = -dsign(dsqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
c     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0d0

         do 240 j = 1, l
            f = d(j)
            z(j,i) = f
            g = e(j) + z(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220

            do 200 k = jp1, l
               g = g + z(k,j) * d(k)
               e(k) = e(k) + z(k,j) * f
  200       continue

  220       e(j) = g
  240    continue
c     .......... form p ..........
         f = 0.0d0

         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue

         hh = f / (h + h)
c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - hh * d(j)
c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)

            do 260 k = j, l
  260       z(k,j) = z(k,j) - f * e(k) - g * d(k)

            d(j) = z(l,j)
            z(i,j) = 0.0d0
  280    continue

  290    d(i) = h
  300 continue
c     .......... accumulation of transformation matrices ..........
      do 500 i = 2, n
         l = i - 1
         z(n,l) = z(l,l)
         z(l,l) = 1.0d0
         h = d(i)
         if (h .eq. 0.0d0) go to 380

         do 330 k = 1, l
  330    d(k) = z(k,i) / h

         do 360 j = 1, l
            g = 0.0d0

            do 340 k = 1, l
  340       g = g + z(k,i) * z(k,j)

            do 360 k = 1, l
               z(k,j) = z(k,j) - g * d(k)
  360    continue

  380    do 400 k = 1, l
  400    z(k,i) = 0.0d0

  500 continue

  510 do 520 i = 1, n
         d(i) = z(n,i)
         z(n,i) = 0.0d0
  520 continue

      z(n,n) = 1.0d0
      e(1) = 0.0d0
      return
      end

c***********************************************************************
      subroutine tql2(nm,n,d,e,z,ierr)
c     This subroutine is a translation of the algol procedure tql2,
c     Num. Math. 11, 293-306(1968) by Bowdler, Martin, Reinsch, and
c     Wilkinson.
c     Handbook for Auto. Comp., Vol.ii-Linear Algebra, 227-240(1971).
c     This subroutine finds the eigenvalues and eigenvectors
c     of a symmetric tridiagonal matrix by the ql method.
c     The eigenvectors of a full symmetric matrix can also
c     be found if  tred2  has been used to reduce this
c     full matrix to tridiagonal form.
c     On input:
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c        n is the order of the matrix.
c        d contains the diagonal elements of the input matrix.
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c        z contains the transformation matrix produced in the
c          reduction by  tred2, if performed.  if the eigenvectors
c          of the tridiagonal matrix are desired, z must contain
c          the identity matrix.
c      On output:
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1,2,...,ierr-1.
c        e has been destroyed.
c        z contains orthonormal eigenvectors of the symmetric
c          tridiagonal (or full) matrix.  if an error exit is made,
c          z contains the eigenvectors associated with the stored
c          eigenvalues.
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c     Calls pythag for  dsqrt(a*a + b*b).
c     Questions and comments should be directed to Burton S. Garbow,
c     Mathematics and Computer Science Div, Argonne National Laboratory.
c     This version dated august 1983.
c***********************************************************************

      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      real*8 d(n),e(n),z(nm,n)
      real*8 c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2
      real*4 pythag

      ierr = 0
      if (n .eq. 1) go to 1001

      do 100 i = 2, n
  100 e(i-1) = e(i)

      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0

      do 240 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue

  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = dble(pythag(p,1.0))
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145

         do 140 i = l2, n
  140    d(i) = d(i) - h

  145    f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = dble(pythag(sngl(p),sngl(e(i))))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
c     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue

  200    continue

         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
c     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)

         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue

         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p

         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue

  300 continue

      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end

c***********************************************************************
      real function pythag(a,b)
c     Function call calculates dsqrt(a*a+b*b).
c***********************************************************************

      real*4 a,b

      pythag = sqrt(a*a + b*b)

      return
      end

