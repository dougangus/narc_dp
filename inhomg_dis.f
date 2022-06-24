c...+....|....+....|....+....|....+....|....+....|....+....|....+....|..
c  File Inhomg_dis.f contains the following 2 subroutines: EcCart and 
c  AngCart.
c
c  Last modified: November 3, 2006.
c***********************************************************************
      subroutine eccart(ix1,nx1st,nx2,nx3,a6,numec)
c     Input elastic constants (divided by density) in Voigt notation
c     of the curvilinear grid.
c***********************************************************************

      implicit none
      include '../Input/narc_dp.par'

      integer ix1,nx1st,nx2,nx3,numec,jx1lo,jx1hi,kx1,kx2,kx3,i,j
      real*8 density
      real*8 a6(6,6,nx3mx,nx2mx,nx1sto)

      jx1lo=1
      jx1hi=nx1st

      do 10 kx1=jx1lo,jx1hi

         do 15 kx2=1,nx2
            do 20 kx3=1,nx3
               do i=1,6
                  do j=1,6
                     a6(i,j,kx3,kx2,kx1)=0.d0
                  enddo
               enddo

               if(numec.gt.3)then

                  read(26)a6(1,1,kx3,kx2,kx1),a6(1,2,kx3,kx2,kx1),
     +                 a6(1,3,kx3,kx2,kx1),a6(2,2,kx3,kx2,kx1),
     +                 a6(2,3,kx3,kx2,kx1),a6(3,3,kx3,kx2,kx1),
     +                 a6(4,4,kx3,kx2,kx1),a6(5,5,kx3,kx2,kx1),
     +                 a6(6,6,kx3,kx2,kx1),a6(1,6,kx3,kx2,kx1),
     +                 a6(2,6,kx3,kx2,kx1),a6(3,6,kx3,kx2,kx1),
     +                 a6(4,5,kx3,kx2,kx1),a6(1,4,kx3,kx2,kx1),
     +                 a6(1,5,kx3,kx2,kx1),a6(2,4,kx3,kx2,kx1),
     +                 a6(2,5,kx3,kx2,kx1),a6(3,4,kx3,kx2,kx1),
     +                 a6(3,5,kx3,kx2,kx1),a6(5,6,kx3,kx2,kx1),
     +                 a6(4,6,kx3,kx2,kx1),density

               elseif(numec.eq.2)then

                  read(26)a6(1,1,kx3,kx2,kx1),a6(5,5,kx3,kx2,kx1)

                  a6(1,2,kx3,kx2,kx1)=a6(1,1,kx3,kx2,kx1)-2.d0*
     +                 a6(5,5,kx3,kx2,kx1)

                  a6(2,2,kx3,kx2,kx1)=a6(1,1,kx3,kx2,kx1)
                  a6(3,3,kx3,kx2,kx1)=a6(1,1,kx3,kx2,kx1)
                  a6(1,3,kx3,kx2,kx1)=a6(1,2,kx3,kx2,kx1)
                  a6(2,3,kx3,kx2,kx1)=a6(1,2,kx3,kx2,kx1)
                  a6(4,4,kx3,kx2,kx1)=a6(5,5,kx3,kx2,kx1)
                  a6(6,6,kx3,kx2,kx1)=a6(5,5,kx3,kx2,kx1)

                  a6(2,1,kx3,kx2,kx1)=a6(1,2,kx3,kx2,kx1) !symmetry
                  a6(3,1,kx3,kx2,kx1)=a6(1,3,kx3,kx2,kx1)
                  a6(3,2,kx3,kx2,kx1)=a6(2,3,kx3,kx2,kx1)

               elseif(numec.eq.3)then

                  read(26)a6(1,1,kx3,kx2,kx1),a6(5,5,kx3,kx2,kx1),
     +                 density

                  a6(1,2,kx3,kx2,kx1)=a6(1,1,kx3,kx2,kx1)-2.d0*
     +                 a6(5,5,kx3,kx2,kx1)

                  a6(2,2,kx3,kx2,kx1)=a6(1,1,kx3,kx2,kx1)
                  a6(3,3,kx3,kx2,kx1)=a6(1,1,kx3,kx2,kx1)
                  a6(1,3,kx3,kx2,kx1)=a6(1,2,kx3,kx2,kx1)
                  a6(2,3,kx3,kx2,kx1)=a6(1,2,kx3,kx2,kx1)
                  a6(4,4,kx3,kx2,kx1)=a6(5,5,kx3,kx2,kx1)
                  a6(6,6,kx3,kx2,kx1)=a6(5,5,kx3,kx2,kx1)

                  a6(2,1,kx3,kx2,kx1)=a6(1,2,kx3,kx2,kx1) !symmetry
                  a6(3,1,kx3,kx2,kx1)=a6(1,3,kx3,kx2,kx1)
                  a6(3,2,kx3,kx2,kx1)=a6(2,3,kx3,kx2,kx1)

               endif

               if((numec.ne.2).and.(kx2.eq.1.and.kx3.eq.1))den=density

 20         continue
 15      continue
 10   continue

      return
      end

c***********************************************************************
      subroutine angcart(ix1,nx1st,nx2,nx3,a6,a3,numec)
c     Input elastic constants (divided by density) in Voigt notation
c     of the curvilinear grid.
c***********************************************************************

      implicit none
      include '../Input/narc_dp.par'

      integer ix1,nx1st,nx2,nx3,numec,jx1lo,jx1hi,kx1,kx2,kx3
      integer i,j
      real*8 alpha,beta,gamma
      real*8 a3(3,3,3,3),a3r1(3,3,3,3),a3r2(3,3,3,3)
      real*8 a6(6,6,nx3mx,nx2mx,nx1sto),a(6,6)

      jx1lo=1
      jx1hi=nx1st

      do 10 kx1=jx1lo,jx1hi
         do 15 kx2=1,nx2
            do 20 kx3=1,nx3
               read(26)alpha,beta
               gamma=0.d0
               call rotate(alpha,0.d0,gamma,a3,a3r1)
               call rotate(0.d0,beta,gamma,a3r1,a3r2)
               call c66(a3r2,a)
               do i=1,6
                  do j=1,6
                     a6(i,j,kx3,kx2,kx1)=a(i,j)
                  enddo
               enddo
 20         continue
 15      continue
 10   continue

      return
      end
