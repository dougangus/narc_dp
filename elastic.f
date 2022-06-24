c........+.........+.........+.........+.........+.........+.........+..
c  File Elastic.f contains the following 8 subroutines: Edo, Emodel, 
c  C66, Aijkl, Caijkl, Isotens, Averps, and Symmetry.
c
c  Last modified: October 8, 2012.
c***********************************************************************
      subroutine edo(ix1,jx1,nx1st,nx1,nx2,nx3,nt,dt,dx1,dx2,dx3,a6)
c     This subroutine interfaces with the subroutine emodel: 
c     1) homogeneous,
c     2) inhomogeneous and analytic, or
c     3) inhomogeneous and discrete (input from model file).
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer ix1,jx1,nx1st,nx1,nx2,nx3,itmp,nt
      real*8 dx1,dx2,dx3,dt
      real*8 a6(6,6,nx3mx,nx2mx,nx1sto)

      if(jx1.ne.-1)then
         if(imtrick.eq.0)then
            if(inhomog.eq.0.or.iuse.eq.1)then
               if(ix1.eq.1)then
                  jx1=1
               endif
            elseif(inhomog.eq.1)then
               if(ix1.le.nx1st)then
                  jx1=ix1
               else
                  jx1=mod(ix1,nx1st)
                  if(jx1.eq.0)then
                     jx1=jx1+nx1st
                  endif
               endif
            endif
         elseif(imtrick.eq.1)then
            jx1=1
         endif
      endif

      if(inhomog.eq.0.or.iuse.eq.1)then
         if(ix1.eq.1)then
            call emodel(ix1,jx1,nx1st,nx1,nx2,nx3,nt,dt,dx1,dx2,dx3,a6)
         endif
      elseif(inhomog.eq.2)then
         if(ix1.eq.1)then
            call emodel(ix1,jx1,nx1st,nx1,nx2,nx3,nt,dt,dx1,dx2,dx3,a6)
         endif
      elseif(inhomog.eq.1)then
         if(imtrick.eq.1.and.ianalyt.eq.0)then
          write(*,*)'itmp in edo',itmp
          stop
          if(ix1.eq.1.or.mod(ix1,itmp).eq.0.and.ix1.ne.nx1)then
            call emodel(ix1,jx1,nx1st,nx1,nx2,nx3,nt,dt,dx1,dx2,dx3,a6)
          endif
         elseif((jx1.eq.1.and.inhomog.ge.1).or.(jx1.eq.-1))then
            call emodel(ix1,jx1,nx1st,nx1,nx2,nx3,nt,dt,dx1,dx2,dx3,a6)
         endif
      endif

      return
      end

c***********************************************************************
      subroutine emodel(ix1,jx1,nx1st,nx1,nx2,nx3,nt,dt,dx1,dx2,dx3,a6)
c     This subroutine interfaces with the various elastic models: 
c     1) homogeneous,
c     2) inhomogeneous and analytic, or
c     3) inhomogeneous and discrete (input from model file).
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer ix1,jx1,nx1st,nx1,nx2,nx3,i,j,numec,ijflag,nt
      real*8 dx1,dx2,dx3,tmp1,tmp2,dt
      real*8 a6(6,6,nx3mx,nx2mx,nx1sto),a3(3,3,3,3),a(6,6)

      if(inhomog.eq.0.or.iuse.eq.1)then
         if(ix1.eq.1)
     +        write(*,*)'             -> using homog or iuse option'
         call inelastc(jx1,nx2,nx3,a6)
      elseif(inhomog.eq.2)then
         if(ix1.eq.1)then
            write(*,*)'             -> using homogeneous frequency'
            write(*,*)'             -> dependent elasticity option'
            call freqelastc(jx1,nt,dt)
         endif
      elseif(inhomog.eq.1)then
         if(ianalyt.eq.0)then   !Discrete model (from input file)
            if(ix1.eq.1)then
               open(unit=26,file=modelcart_e,form='unformatted')
               read(26)numec,ncx1,dcx1,ncx2,dcx2,ncx3,dcx3
            endif
            if(numec.eq.2)read(26)den
            if(imtrick.eq.0.and.ix1.eq.1)then
               if(ncx1.lt.nx1.or.ncx2.ne.nx2.or.ncx3.ne.nx3)then
                  write(11,*)
                  write(11,*)'Problem with Cartesian grid dimensions'
                  write(11,*)'extrapolation only.'
                  write(11,*)' nx1, nx2, nx3:',nx1,nx2,nx3
                  write(11,*)'ncx1,ncx2,ncx3:',ncx1,ncx2,ncx3
                  write(11,*)
                  stop
               elseif(dcx1.ne.dx1.or.dcx2.ne.dx2.or.dcx3.ne.dx3)then
                  write(11,*)
                  write(11,*)'Problem with Cartesian grid dimensions'
                  write(11,*)'extrapolation only.'
                  write(11,*)' dx1, dx2, dx3:',dx1,dx2,dx3
                  write(11,*)'dcx1,dcx2,dcx3:',dcx1,dcx2,dcx3
                  write(11,*)
                  stop
               endif
            elseif(imtrick.eq.1.and.ix1.eq.1)then
               if(ncx2.ne.nx2.or.ncx3.ne.nx3)then
                  write(11,*)
                  write(11,*)'Problem with Cartesian grid dimensions'
                  write(11,*)'extrapolation only.'
                  write(11,*)'  nx2, nx3:',nx2,nx3
                  write(11,*)' ncx2,ncx3:',ncx2,ncx3
                  write(11,*)
                  stop
               elseif(dcx2.ne.dx2.or.dcx3.ne.dx3)then
                  write(11,*)
                  write(11,*)'Problem with Cartesian grid dimensions'
                  write(11,*)'extrapolation only.'
                  write(11,*)'  dx2, dx3:',dx2,dx3
                  write(11,*)' dcx2,dcx3:',dcx2,dcx3
                  write(11,*)
                  stop
               endif
               tmp1=dcx1/dx1
               if(mod(tmp1,1.0).ne.0.or.tmp1.lt.1.0)then
                  write(11,*)
                  write(11,*)'dcx1,dx1,dcx1/dx1:',dcx1,dx1,tmp1
                  write(11,*)'Problem with grid increments between'
                  write(11,*)'input file and input model: dx1 not '
                  write(11,*)'an exact factor of dcx1.'
                  write(11,*)
                  stop
               else
                  tmp1=real(ncx1)*dcx1
                  tmp2=real(nx1)*dx1
                  if(tmp1.ne.tmp2)then
                     write(11,*)
                     write(11,*)'Problem with # of x1 nodes in input'
                     write(11,*)'and model files: total depth not an'
                     write(11,*)'exact match.'
                     write(11,*)'x1 for model file:',tmp1
                     write(11,*)'x1 for input file:',tmp2
                     write(11,*)
                  endif
               endif
            endif
            if(jx1.eq.-1)then
               write(*,*)'             -> calling discrete emodel from'
               write(*,*)'                input file'
            endif
            call eccart(nx1st,nx2,nx3,a6,numec)
c     Next command to close model input file after initial read for 
c     elastic tensor used in evaluating wavefront initial conditions
            if(jx1.eq.-1)then
               close(26)
               do i=1,6
                  do j=1,6
                     a(i,j)=a6(i,j,1,1,1)
                  enddo
               enddo
               call aijkl(a,a3)
               call averps(a3,vfast,vslow)
            endif
         elseif(ianalyt.eq.1)then
            if(imodel.eq.1)then
               if(jx1.eq.-1)
     +              write(*,*)'             -> using elast_axsh emodel'
               call elast_axsh(ix1,dx1,nx1,nx1st,nx2,nx3,a6)
            elseif(imodel.eq.2)then
               if(jx1.eq.-1)
     +              write(*,*)'             -> using halite_lat emodel'
               call halite_lat(ix1,dx1,nx1,dx2,nx2,dx3,nx3,nx1st,a6)
            elseif(imodel.eq.3)then
               if(jx1.eq.-1)
     +              write(*,*)'             -> using sphere emodel'
               call sphere(ix1,dx1,nx1,dx2,nx2,dx3,nx3,nx1st,a6)
            elseif(imodel.eq.4)then
               if(jx1.eq.-1)
     +              write(*,*)'             -> using anelast emodel'
               call anelast(ix1,dx1,nx1,nx2,nx3,nx1st,a6)
            elseif(imodel.eq.5)then
               ijflag=0
               if(jx1.eq.-1)
     +              write(*,*)'             -> using spetzler emodel'
               call spetzler(ix1,dx1,dx2,nx2,dx3,nx3,nx1st,a6,
     +              ijflag)
            elseif(imodel.eq.7)then
               if(jx1.eq.-1)
     +              write(*,*)'             -> using recfnc emodel'
               call recfnc(ix1,dx1,nx1,dx2,nx2,dx3,nx3,nx1st,a6)
            elseif(imodel.eq.8)then
               if(jx1.eq.-1)
     +              write(*,*)'             -> using afar emodel'
c               call afar(jx1,ix1,dx2,nx2,dx3,nx3,nx1st,a6)
               call afar_alt(jx1,ix1,dx2,nx2,dx3,nx3,nx1st,a6)
            elseif(imodel.eq.9)then
               if(jx1.eq.-1)
     +              write(*,*)'             -> using karakoram emodel'
               call karakoram(jx1,ix1,dx2,nx2,dx3,nx3,nx1st,a6)
            endif
         endif
         if(jx1.eq.-1)write(*,*)
      endif

      return
      end

c***********************************************************************
      subroutine c66(a3,as6)
c     Returns the a_ij (Voigt) given the a_ijkl.  
c***********************************************************************

      implicit none

      integer i,j
      real*8 as6(6,6),a3(3,3,3,3)

      as6(1,1)=a3(1,1,1,1)
      as6(2,2)=a3(2,2,2,2)
      as6(3,3)=a3(3,3,3,3)
      as6(1,2)=a3(1,1,2,2)
      as6(1,3)=a3(1,1,3,3)
      as6(2,3)=a3(2,2,3,3)
  
      as6(6,6)=a3(1,2,1,2)
      as6(5,5)=a3(1,3,1,3)
      as6(4,4)=a3(2,3,2,3)
 
      as6(1,4)=a3(1,1,2,3)
      as6(1,5)=a3(1,1,1,3)
      as6(1,6)=a3(1,1,1,2)
 
      as6(2,4)=a3(2,2,2,3)
      as6(2,5)=a3(2,2,1,3)
      as6(2,6)=a3(2,2,1,2)
 
      as6(3,4)=a3(3,3,2,3)
      as6(3,5)=a3(3,3,1,3)
      as6(3,6)=a3(3,3,1,2)
 
      as6(4,5)=a3(2,3,1,3)
      as6(4,6)=a3(2,3,1,2)

      as6(5,6)=a3(1,3,1,2)

c     Impose symmetry 

      do i=2,6
         do j=1,i-1
            as6(i,j)=as6(j,i)
         enddo
      enddo

      return
      end


c***********************************************************************
      subroutine ac66(a3,as6)
c     Returns the a_ij (Voigt) given the a_ijkl.  
c***********************************************************************

      implicit none

      integer i,j
      complex*8 as6(6,6),a3(3,3,3,3)

      as6(1,1)=a3(1,1,1,1)
      as6(2,2)=a3(2,2,2,2)
      as6(3,3)=a3(3,3,3,3)
      as6(1,2)=a3(1,1,2,2)
      as6(1,3)=a3(1,1,3,3)
      as6(2,3)=a3(2,2,3,3)
  
      as6(6,6)=a3(1,2,1,2)
      as6(5,5)=a3(1,3,1,3)
      as6(4,4)=a3(2,3,2,3)
 
      as6(1,4)=a3(1,1,2,3)
      as6(1,5)=a3(1,1,1,3)
      as6(1,6)=a3(1,1,1,2)
 
      as6(2,4)=a3(2,2,2,3)
      as6(2,5)=a3(2,2,1,3)
      as6(2,6)=a3(2,2,1,2)
 
      as6(3,4)=a3(3,3,2,3)
      as6(3,5)=a3(3,3,1,3)
      as6(3,6)=a3(3,3,1,2)
 
      as6(4,5)=a3(2,3,1,3)
      as6(4,6)=a3(2,3,1,2)

      as6(5,6)=a3(1,3,1,2)

c     Impose symmetry 

      do i=2,6
         do j=1,i-1
            as6(i,j)=as6(j,i)
         enddo
      enddo

      return
      end

c***********************************************************************
      subroutine aijkl(as6,a3)
c     Returns the a_ijkl given the a_mn Voigt notation.  First, the 
c     upper triangle of a_mn is used then imposes the general 
c     symmetries that apply.  Does not assume particular crystal 
c     symmetry (e.g. hexagonal).
c***********************************************************************

      implicit none

      integer i,j,k,l
      real*8 as6(6,6),a3(3,3,3,3)

      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  a3(i,j,k,l)=0.d0
               enddo
            enddo
         enddo
      enddo

      a3(1,1,1,1)=as6(1,1)
      a3(2,2,2,2)=as6(2,2)
      a3(3,3,3,3)=as6(3,3)
      a3(1,1,2,2)=as6(1,2)
      a3(1,1,3,3)=as6(1,3)
      a3(2,2,3,3)=as6(2,3)
  
      a3(1,2,1,2)=as6(6,6)
      a3(1,3,1,3)=as6(5,5)
      a3(2,3,2,3)=as6(4,4)
 
      a3(1,1,2,3)=as6(1,4)
      a3(1,1,1,3)=as6(1,5)
      a3(1,1,1,2)=as6(1,6)
 
      a3(2,2,2,3)=as6(2,4)
      a3(2,2,1,3)=as6(2,5)
      a3(2,2,1,2)=as6(2,6)
 
      a3(3,3,2,3)=as6(3,4)
      a3(3,3,1,3)=as6(3,5)
      a3(3,3,1,2)=as6(3,6)
 
      a3(2,3,1,3)=as6(4,5)
      a3(2,3,1,2)=as6(4,6)

      a3(1,3,1,2)=as6(5,6)

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
      subroutine caijkl(as6,a3r,a3)
c     Returns the complex a_ijkl given the complex a_mn Voigt notation. 
c     First, the upper triangle of a_mn is used then imposes the general 
c     symmetries that apply.  Does not assume particular crystal 
c     symmetry (e.g. hexagonal).
c***********************************************************************

      implicit none

      integer i,j,k,l
      real*8 a3r(3,3,3,3)
      complex*8 as6(6,6),a3(3,3,3,3)

      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  a3(i,j,k,l)=cmplx(0.d0,0.d0)
               enddo
            enddo
         enddo
      enddo

      a3(1,1,1,1)=as6(1,1)
      a3(2,2,2,2)=as6(2,2)
      a3(3,3,3,3)=as6(3,3)
      a3(1,1,2,2)=as6(1,2)
      a3(1,1,3,3)=as6(1,3)
      a3(2,2,3,3)=as6(2,3)
  
      a3(1,2,1,2)=as6(6,6)
      a3(1,3,1,3)=as6(5,5)
      a3(2,3,2,3)=as6(4,4)
 
      a3(1,1,2,3)=as6(1,4)
      a3(1,1,1,3)=as6(1,5)
      a3(1,1,1,2)=as6(1,6)
 
      a3(2,2,2,3)=as6(2,4)
      a3(2,2,1,3)=as6(2,5)
      a3(2,2,1,2)=as6(2,6)
 
      a3(3,3,2,3)=as6(3,4)
      a3(3,3,1,3)=as6(3,5)
      a3(3,3,1,2)=as6(3,6)
 
      a3(2,3,1,3)=as6(4,5)
      a3(2,3,1,2)=as6(4,6)

      a3(1,3,1,2)=as6(5,6)

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

      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  a3r(i,j,k,l)=dble(realpart(a3(i,j,k,l)))
               enddo
            enddo
         enddo
      enddo

      return
      end

c***********************************************************************
      subroutine isotens(as6)
c     Elasticity of an isotropic medium.  Elasticity is divided by 
c     density.
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer i,j
      real*8 a(6,6),as6(6,6)

      a(1,1)=vfast**2.d0
      a(4,4)=vslow**2.d0
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
      subroutine averps(cc,vpccav,vsccav) 
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
     .     print *,'error energy 1 in averps'
      if(cc(2,2,2,2).lt.0.d0)
     .     print *,'error energy 2 in averps'
      if(cc(3,3,3,3).lt.0.d0)
     .     print *,'error energy 3 in averps'
      if(cc(2,3,2,3).lt.0.d0)
     .     print *,'error energy 4 in averps'
      if(cc(1,2,1,2).lt.0.d0)
     .     print *,'error energy 5 in averps'
      if(cc(1,3,1,3).lt.0.d0)
     .     print *,'error energy 6 in averps'
      test1=cc(1,1,1,1)*cc(2,2,2,2)-cc(1,1,2,2)**2
      if(test1.lt.0.d0)
     .     print *,'error energy 7 in averps'
         
      test2=cc(1,1,2,2)*cc(2,2,3,3)
     .     -cc(2,2,2,2)*cc(1,1,3,3)
      test3=cc(1,1,1,1)*cc(2,2,3,3)
     .     -cc(1,1,2,2)*cc(1,1,3,3)
      test=test1*cc(3,3,3,3) + test2*cc(1,1,3,3)
     .     - test3*cc(2,2,3,3)
      if(test.lt.0.d0)
     .     print *,'error energy 8 in averps'

      return
      end

c***********************************************************************
      subroutine symmetry(isym,as6,a)
c     Returns the a_ijkl given the a_mn Voigt notation.  First, the 
c     upper triangle of a_mn is used then imposes the general 
c     symmetries that apply.  Does not assume particular crystal 
c     symmetry (e.g. hexagonal).
c***********************************************************************

      implicit none

      integer isym,i,j
      real*8 as6(6,6),a(3,3,3,3)

      if(isym.eq.0)then         !For a3(3,3,3,3)
         a(2,2,1,1)=a(1,1,2,2)
         a(3,3,1,1)=a(1,1,3,3)
         a(3,3,2,2)=a(2,2,3,3)
         a(2,1,1,2)=a(1,2,1,2)
         a(1,2,2,1)=a(1,2,1,2)
         a(2,1,2,1)=a(1,2,1,2)
         a(3,1,1,3)=a(1,3,1,3)
         a(1,3,3,1)=a(1,3,1,3)
         a(3,1,3,1)=a(1,3,1,3)
         a(3,2,2,3)=a(2,3,2,3)
         a(2,3,3,2)=a(2,3,2,3)
         a(3,2,3,2)=a(2,3,2,3)
         a(1,1,3,2)=a(1,1,2,3)
         a(2,3,1,1)=a(1,1,2,3)
         a(3,2,1,1)=a(1,1,2,3)
         a(1,1,3,1)=a(1,1,1,3)
         a(1,3,1,1)=a(1,1,1,3)
         a(3,1,1,1)=a(1,1,1,3)
         a(1,1,2,1)=a(1,1,1,2)
         a(1,2,1,1)=a(1,1,1,2)
         a(2,1,1,1)=a(1,1,1,2)
         a(2,2,3,2)=a(2,2,2,3)
         a(2,3,2,2)=a(2,2,2,3)
         a(3,2,2,2)=a(2,2,2,3)
         a(2,2,3,1)=a(2,2,1,3)
         a(1,3,2,2)=a(2,2,1,3)
         a(3,1,2,2)=a(2,2,1,3)
         a(2,2,2,1)=a(2,2,1,2)
         a(1,2,2,2)=a(2,2,1,2)
         a(2,1,2,2)=a(2,2,1,2)
         a(3,3,3,2)=a(3,3,2,3)
         a(2,3,3,3)=a(3,3,2,3)
         a(3,2,3,3)=a(3,3,2,3)
         a(3,3,3,1)=a(3,3,1,3)
         a(1,3,3,3)=a(3,3,1,3)
         a(3,1,3,3)=a(3,3,1,3)
         a(3,3,2,1)=a(3,3,1,2)
         a(1,2,3,3)=a(3,3,1,2)
         a(2,1,3,3)=a(3,3,1,2)
         a(2,3,3,1)=a(2,3,1,3)
         a(3,2,1,3)=a(2,3,1,3)
         a(3,2,3,1)=a(2,3,1,3)
         a(1,3,2,3)=a(2,3,1,3)
         a(1,3,3,2)=a(2,3,1,3)
         a(3,1,2,3)=a(2,3,1,3)
         a(3,1,3,2)=a(2,3,1,3)
         a(2,3,2,1)=a(2,3,1,2)
         a(3,2,1,2)=a(2,3,1,2)
         a(3,2,2,1)=a(2,3,1,2)
         a(1,2,2,3)=a(2,3,1,2)
         a(1,2,3,2)=a(2,3,1,2)
         a(2,1,2,3)=a(2,3,1,2)
         a(2,1,3,2)=a(2,3,1,2)
         a(3,1,1,2)=a(1,3,1,2)
         a(1,3,2,1)=a(1,3,1,2)
         a(3,1,2,1)=a(1,3,1,2)
         a(1,2,1,3)=a(1,3,1,2)
         a(2,1,1,3)=a(1,3,1,2)
         a(1,2,3,1)=a(1,3,1,2)
         a(2,1,3,1)=a(1,3,1,2)
      else                      !For a6(6,6)
         do i=2,6
            do j=1,i-1
               as6(i,j)=as6(j,i)
            enddo
         enddo
      endif

      return
      end
