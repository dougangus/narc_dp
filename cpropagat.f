c...+....|....+....|....+....|....+....|....+....|....+....|....+....|..
c  File Cpropagat.f contains the following 4 subroutines: cpropagat, 
c  csolvepa, csolvepab and csymtry.
c
c  Last modified: May 23, 2012.
c***********************************************************************
      subroutine cpropagat(a3,p0,p2,p3,p22,p23,p32,p33)
c     Calculates the variable coefficients of the one-way wave
c     equation.
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer i,j,ierr
      real*4 densqrt,sden,cmodul
      complex*8 ci,a3(3,3,3,3),c11(3,3),c11inv(3,3)
      complex*8 c12(3,3),c21(3,3),c13(3,3),c31(3,3)
      complex*8 c23(3,3),c32(3,3),c22(3,3),c33(3,3)
      complex*8 cval(3),cvals(3,3),cvalsqrt(3,3),cvalsqin(3,3)
      complex*8 cvecs(3,3),cvecstrn(3,3)
      complex*8 p0(3,3),p2(3,3),p3(3,3),p22(3,3),p23(3,3),p32(3,3),
     +     p33(3,3)
      real*8 ar(3,3),ai(3,3)
      real*8 wr(3),wi(3),zr(3,3),zi(3,3),fv1(3),fv2(3),
     +     fv3(3)

      sden=sngl(den)
      densqrt=sqrt(sden)
      do i=1,3
         do j=1,3
            c11(i,j) = a3(i,1,1,j)*sden
            c22(i,j) = a3(i,2,2,j)*sden
            c33(i,j) = a3(i,3,3,j)*sden
            c12(i,j) = a3(i,1,2,j)*sden
            c21(i,j) = a3(i,2,1,j)*sden
            c13(i,j) = a3(i,1,3,j)*sden
            c31(i,j) = a3(i,3,1,j)*sden
            c23(i,j) = a3(i,2,3,j)*sden
            c32(i,j) = a3(i,3,2,j)*sden
         enddo
      enddo

c     Evaluate complex eigensolution
      ci=cmplx(0.0,1.0)
      do i=1,3
         do j=1,3
            ar(i,j)=dble(realpart(c11(i,j)))
            ai(i,j)=dble(imagpart(c11(i,j)))
         enddo
      enddo
      call cg(3,3,ar,ai,wr,wi,1,zr,zi,fv1,fv2,fv3,ierr)
      do i=1,3
         cval(i)=cmplx(wr(i))+ci*cmplx(wi(i))
c     Normalize eigenvectors (they are unnormalized)
         cmodul=sqrt(sngl(
     +        zr(i,1)**2.d0+zr(i,2)**2.d0+zr(i,3)**2.d0+
     +        zi(i,1)**2.d0+zi(i,2)**2.d0+zi(i,3)**2.d0))
         do j=1,3
            cvecs(i,j)=cmplx(zr(i,j))/cmodul+ci*cmplx(zi(i,j))/cmodul
         enddo
      enddo

      if(ierr.ne.0) then
         write(11,*)'Eigenvalue/eigenvector error not equal to 0.'
         stop
      endif
      do i=1,3
         do j=1,3
            if(i.eq.j) cvals(i,j) = cval(i)
            if(i.eq.j) cvalsqrt(i,j) = csqrt(cvals(i,j))
            if(i.ne.j) cvals(i,j) = cmplx(0.0,0.0)
            if(i.ne.j) cvalsqrt(i,j) = cmplx(0.0,0.0)
         enddo
      enddo

c     Need to convert to complex value subroutines
      call cinv3(c11,c11inv)
      call cinv3(cvalsqrt,cvalsqin)
      call ctrans3(cvecs,cvecstrn)
      call cmult3b(cvecs,cvalsqin,cvecstrn,p0)

      do i=1,3
         do j=1,3
            p0(i,j) = p0(i,j)*densqrt
         enddo
c         write(*,*)(cvecs(i,j),j=1,3)
c         write(*,*)(p0(i,j),j=1,3)
c         read(*,*)
      enddo
      call csolvepa(cvecstrn,c11inv,c12,c21,cvecs,cvalsqin,p2)
      call csolvepa(cvecstrn,c11inv,c13,c31,cvecs,cvalsqin,p3)
      call csolvepab(cvecstrn,c11inv,c22,c22,c12,c21,c12,c21,
     +     cvecs,cvalsqin,densqrt,p2,p2,p22)
      call csolvepab(cvecstrn,c11inv,c23,c32,c12,c21,c13,c31,
     +     cvecs,cvalsqin,densqrt,p2,p3,p23)
      call csolvepab(cvecstrn,c11inv,c33,c33,c13,c31,c13,c31,
     +     cvecs,cvalsqin,densqrt,p3,p3,p33)
      do i=1,3
         do j=1,3
            p32(i,j)=p23(i,j)
         enddo
      enddo

      return
      end

c***********************************************************************
      subroutine csolvepa(cvecstrn,c11inv,c1a,ca1,cvecs,cvalsqin,pa)
c     Evaluate the local P_alpha coefficient matrix.
c***********************************************************************

      implicit none

      integer i,j
      complex*8 cvecs(3,3),cvecstrn(3,3),cvalsqin(3,3),c11inv(3,3)
      complex*8 c1a(3,3),ca1(3,3)
      complex*8 pa(3,3)
      complex*8 qtemp(3,3),temp1(3,3),temp2(3,3)

      call cmult3d(cvecstrn,c11inv,c1a,cvecs,cvalsqin,temp1)
      call cmult3d(cvecstrn,c11inv,ca1,cvecs,cvalsqin,temp2)
      do i=1,3
         do j=1,3
            qtemp(i,j)=-temp1(i,j)-temp2(i,j)
         enddo
      enddo
      call csymtry(qtemp,cvalsqin,pa)
      call cmult3b(cvecs,pa,cvecstrn,pa)

      return
      end

c***********************************************************************
      subroutine csolvepab(cvecstrn,c11inv,cab,cba,c1a,ca1,c1b,cb1,
     +     cvecs,cvalsqin,rhosqrt,pa,pb,pab)
c     Evaluates the local P_alpha_beta coefficient matrix.
c***********************************************************************

      implicit none

      integer i,j
      real*4 rhosqrt
      complex*8 cvecs(3,3),cvecstrn(3,3),cvalsqin(3,3),c11inv(3,3)
      complex*8 cab(3,3),cba(3,3),c1a(3,3),ca1(3,3),c1b(3,3),cb1(3,3)
      complex*8 ct1(3,3),ct2(3,3),ct3(3,3)
      complex*8 pa(3,3),pb(3,3),pab(3,3)
      complex*8 qtemp(3,3),temp1(3,3),temp2(3,3),temp3(3,3),
     +     temp4(3,3),temp5(3,3)

      do i=1,3
         do j=1,3
            ct1(i,j)=cab(i,j)+cba(i,j)
            ct2(i,j)=c1a(i,j)+ca1(i,j)
            ct3(i,j)=c1b(i,j)+cb1(i,j)
         enddo
      enddo
      call cmult3c(cvecstrn,c11inv,ct1,cvecs,temp1)
      call cmult3d(cvecstrn,c11inv,ct2,pb,cvecs,temp2) 
      call cmult3d(cvecstrn,c11inv,ct3,pa,cvecs,temp3) 
      call cmult3c(cvecstrn,pa,pb,cvecs,temp4)
      call cmult3c(cvecstrn,pb,pa,cvecs,temp5)
      do i=1,3
         do j=1,3
            qtemp(i,j)=(-1.0/rhosqrt)*(0.5*temp1(i,j)+
     +           0.5*temp2(i,j)+0.5*temp3(i,j)+
     +           0.5*temp4(i,j)+0.5*temp5(i,j))
         enddo
      enddo
      call csymtry(qtemp,cvalsqin,pab)
      call cmult3b(cvecs,pab,cvecstrn,pab)

      return
      end

c***********************************************************************
      subroutine csymtry(qtemp,cvalsqin,p) 
c     Evaluates the P-matrix values in two stages using the symmetric
c     and antisymmetric parts of the defining equation.
c***********************************************************************

      implicit none

      integer i,j
      complex*8 cvalsqin(3,3)
      complex*8 p(3,3),panti(3,3),psym(3,3)
      complex*8 qtemp(3,3),qanti(3,3),qsym(3,3)

      do i=1,3
         do j=1,3
            qanti(i,j) = (qtemp(i,j) - qtemp(j,i))/2.0
            qsym(i,j) = (qtemp(i,j) + qtemp(j,i))/2.0
         enddo
      enddo

      panti(1,1) = cmplx(0.0,0.0)
      panti(2,2) = cmplx(0.0,0.0)
      panti(3,3) = cmplx(0.0,0.0)
      panti(1,2) = qanti(1,2)/(cvalsqin(1,1)+cvalsqin(2,2))
      panti(1,3) = qanti(1,3)/(cvalsqin(1,1)+cvalsqin(3,3))
      panti(2,3) = qanti(2,3)/(cvalsqin(2,2)+cvalsqin(3,3))
      panti(2,1) = -panti(1,2)
      panti(3,1) = -panti(1,3)
      panti(3,2) = -panti(2,3)

      psym(1,1) = qsym(1,1)/(2.0*cvalsqin(1,1))
      psym(2,2) = qsym(2,2)/(2.0*cvalsqin(2,2))
      psym(3,3) = qsym(3,3)/(2.0*cvalsqin(3,3))
      psym(1,2) = qsym(1,2)/(cvalsqin(1,1)+cvalsqin(2,2))
      psym(1,3) = qsym(1,3)/(cvalsqin(1,1)+cvalsqin(3,3))
      psym(2,3) = qsym(2,3)/(cvalsqin(2,2)+cvalsqin(3,3))
      psym(2,1) = psym(1,2)
      psym(3,1) = psym(1,3)
      psym(3,2) = psym(2,3)

      do i=1,3
         do j=1,3
            p(i,j) = psym(i,j) + panti(i,j)
         enddo
      enddo

      return
      end
