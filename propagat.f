c...+....|....+....|....+....|....+....|....+....|....+....|....+....|..
c  File Propagat.f contains the following 4 subroutines: Propagat, 
c  Solvepa, Solvepab and Symtry.
c
c  Last modified: October 31, 2006.
c***********************************************************************
      subroutine propagat(a3,c11inv,p0,p2,p3,p22,p23,p32,p33)
c     Calculates the variable coefficients of the one-way wave
c     equation.
c***********************************************************************

      implicit none
      include '../Input/narc_dp.par'

      integer i,j,ierr
      real*8 densqrt
      real*8 a3(3,3,3,3),c11(3,3),c11inv(3,3)
      real*8 c12(3,3),c21(3,3),c13(3,3),c31(3,3)
      real*8 c23(3,3),c32(3,3),c22(3,3),c33(3,3)
      real*8 cval(3),cvals(3,3),cvalsqrt(3,3),cvalsqin(3,3)
      real*8 cvecs(3,3),cvecstrn(3,3),ctemp(3)
      real*8 p0(3,3),p2(3,3),p3(3,3),p22(3,3),p23(3,3),p32(3,3),
     +     p33(3,3)

      densqrt=dsqrt(den)

      do i=1,3
         do j=1,3
            c11(i,j) = a3(i,1,1,j)*den     !According to Thomson 1999
            c22(i,j) = a3(i,2,2,j)*den
            c33(i,j) = a3(i,3,3,j)*den
            c12(i,j) = a3(i,1,2,j)*den
            c21(i,j) = a3(i,2,1,j)*den
            c13(i,j) = a3(i,1,3,j)*den
            c31(i,j) = a3(i,3,1,j)*den
            c23(i,j) = a3(i,2,3,j)*den
            c32(i,j) = a3(i,3,2,j)*den
         enddo
      enddo

      if(itropic.eq.0)then

         if(inhomog.eq.1)call averps(a3,vfast,vslow)

         do i=1,3
            do j=1,3
               if(i.ne.j) cvals(i,j) = 0.d0
               if(i.ne.j) cvalsqrt(i,j) = 0.d0
               if(i.eq.1.and.j.eq.1) cvals(i,j) = vfast*vfast*den
               if(i.eq.2.and.j.eq.2) cvals(i,j) = vslow*vslow*den
               if(i.eq.3.and.j.eq.3) cvals(i,j) = vslow*vslow*den
               if(i.eq.1.and.j.eq.1) cvalsqrt(i,j) = 
     +              dsqrt(cvals(i,j)) 
               if(i.eq.2.and.j.eq.2) cvalsqrt(i,j) = 
     +              dsqrt(cvals(i,j))
               if(i.eq.3.and.j.eq.3) cvalsqrt(i,j) = 
     +              dsqrt(cvals(i,j))
               cvecs(1,1) = 1.d0
               cvecs(2,1) = 0.d0
               cvecs(3,1) = 0.d0
               cvecs(1,2) = 0.d0
               cvecs(2,2) = 1.d0
               cvecs(3,2) = 0.d0
               cvecs(1,3) = 0.d0
               cvecs(2,3) = 0.d0
               cvecs(3,3) = 1.d0
            enddo
         enddo

      elseif(itropic.eq.1)then !Using eispack routines

         call tred2(3,3,c11,cval,ctemp,cvecs)
         call tql2(3,3,cval,ctemp,cvecs,ierr)
         if(ierr.ne.0) then
            write(11,*)'Eigenvalue/eigenvector error not equal to 0.'
            stop
         endif
         do i=1,3
            do j=1,3
               if(i.eq.j) cvals(i,j) = cval(i)
               if(i.eq.j) cvalsqrt(i,j) = dsqrt(cvals(i,j))
               if(i.ne.j) cvals(i,j) = 0.d0
               if(i.ne.j) cvalsqrt(i,j) = 0.d0
            enddo
         enddo

      endif

      call rinv3(c11,c11inv)
      call rinv3(cvalsqrt,cvalsqin)
      call trans3(cvecs,cvecstrn)
      call rmult3b(cvecs,cvalsqin,cvecstrn,p0)

      do i=1,3
         do j=1,3
            p0(i,j) = p0(i,j)*densqrt
         enddo
      enddo
      call solvepa(cvecstrn,c11inv,c12,c21,cvecs,cvalsqin,p2)
      call solvepa(cvecstrn,c11inv,c13,c31,cvecs,cvalsqin,p3)
      call solvepab(cvecstrn,c11inv,c22,c22,c12,c21,c12,c21,
     +     cvecs,cvalsqin,densqrt,p2,p2,p22)
      call solvepab(cvecstrn,c11inv,c23,c32,c12,c21,c13,c31,
     +     cvecs,cvalsqin,densqrt,p2,p3,p23)
      call solvepab(cvecstrn,c11inv,c33,c33,c13,c31,c13,c31,
     +     cvecs,cvalsqin,densqrt,p3,p3,p33)
      do i=1,3
         do j=1,3
            p32(i,j)=p23(i,j)
         enddo
      enddo

      return
      end

c***********************************************************************
      subroutine solvepa(cvecstrn,c11inv,c1a,ca1,cvecs,cvalsqin,pa)
c     Evaluate the local P_alpha coefficient matrix.
c***********************************************************************

      implicit none

      integer i,j
      real*8 cvecs(3,3),cvecstrn(3,3),cvalsqin(3,3),c11inv(3,3)
      real*8 c1a(3,3),ca1(3,3)
      real*8 pa(3,3)
      real*8 qtemp(3,3),temp1(3,3),temp2(3,3)

      call rmult3d(cvecstrn,c11inv,c1a,cvecs,cvalsqin,temp1)
      call rmult3d(cvecstrn,c11inv,ca1,cvecs,cvalsqin,temp2)
      do i=1,3
         do j=1,3
            qtemp(i,j)=-temp1(i,j)-temp2(i,j)
         enddo
      enddo
      call symtry(qtemp,cvalsqin,pa)
      call rmult3b(cvecs,pa,cvecstrn,pa)

      return
      end

c***********************************************************************
      subroutine solvepab(cvecstrn,c11inv,cab,cba,c1a,ca1,c1b,cb1,
     +     cvecs,cvalsqin,rhosqrt,pa,pb,pab)
c     Evaluates the local P_alpha_beta coefficient matrix.
c***********************************************************************

      implicit none

      integer i,j
      real*8 rhosqrt
      real*8 cvecs(3,3),cvecstrn(3,3),cvalsqin(3,3),c11inv(3,3)
      real*8 cab(3,3),cba(3,3),c1a(3,3),ca1(3,3),c1b(3,3),cb1(3,3)
      real*8 ct1(3,3),ct2(3,3),ct3(3,3)
      real*8 pa(3,3),pb(3,3),pab(3,3)
      real*8 qtemp(3,3),temp1(3,3),temp2(3,3),temp3(3,3),
     +     temp4(3,3),temp5(3,3)

      do i=1,3
         do j=1,3
            ct1(i,j)=cab(i,j)+cba(i,j)
            ct2(i,j)=c1a(i,j)+ca1(i,j)
            ct3(i,j)=c1b(i,j)+cb1(i,j)
         enddo
      enddo
      call rmult3c(cvecstrn,c11inv,ct1,cvecs,temp1)
      call rmult3d(cvecstrn,c11inv,ct2,pb,cvecs,temp2) 
      call rmult3d(cvecstrn,c11inv,ct3,pa,cvecs,temp3) 
      call rmult3c(cvecstrn,pa,pb,cvecs,temp4)
      call rmult3c(cvecstrn,pb,pa,cvecs,temp5)
      do i=1,3
         do j=1,3
            qtemp(i,j)=(-1.d0/rhosqrt)*(0.5d0*temp1(i,j)+
     +           0.5d0*temp2(i,j)+0.5d0*temp3(i,j)+
     +           0.5d0*temp4(i,j)+0.5d0*temp5(i,j))
         enddo
      enddo
      call symtry(qtemp,cvalsqin,pab)
      call rmult3b(cvecs,pab,cvecstrn,pab)

      return
      end

c***********************************************************************
      subroutine symtry(qtemp,cvalsqin,p) 
c     Evaluates the P-matrix values in two stages using the symmetric
c     and antisymmetric parts of the defining equation.
c***********************************************************************

      implicit real*8(a-h,o-z)

      dimension cvalsqin(3,3)
      dimension p(3,3),panti(3,3),psym(3,3)
      dimension qtemp(3,3),qanti(3,3),qsym(3,3)

      do i=1,3
         do j=1,3
            qanti(i,j) = (qtemp(i,j) - qtemp(j,i))/2.d0
            qsym(i,j) = (qtemp(i,j) + qtemp(j,i))/2.d0
         enddo
      enddo

      panti(1,1) = 0.d0
      panti(2,2) = 0.d0
      panti(3,3) = 0.d0
      panti(1,2) = qanti(1,2)/(cvalsqin(1,1)+cvalsqin(2,2))
      panti(1,3) = qanti(1,3)/(cvalsqin(1,1)+cvalsqin(3,3))
      panti(2,3) = qanti(2,3)/(cvalsqin(2,2)+cvalsqin(3,3))
      panti(2,1) = -panti(1,2)
      panti(3,1) = -panti(1,3)
      panti(3,2) = -panti(2,3)

      psym(1,1) = qsym(1,1)/(2*cvalsqin(1,1))
      psym(2,2) = qsym(2,2)/(2*cvalsqin(2,2))
      psym(3,3) = qsym(3,3)/(2*cvalsqin(3,3))
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
