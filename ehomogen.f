c...+....|....+....|....+....|....+....|....+....|....+....|....+....|..
c  File Ehomogen.f contains the following 2 subroutines: Inelastc and 
c  Effectiv.
c
c  Last modified: July 16, 2007.
c***********************************************************************
      subroutine inelastc(ix1,jx1,nx2,nx3,a6)
c     Input elastic constants (divided by density) in Voigt notation.
c***********************************************************************

      implicit none
      include '../Input/narc_dp.par'

      integer ix1,jx1,nx2,nx3,kx2,kx3,i,j
      real*8 a6(6,6,nx3mx,nx2mx,nx1sto)
      real*8 a3(3,3,3,3),a3r(3,3,3,3),a(6,6)

      if(jx1.eq.-1)then         !Choice of elasticity tensor
      write(*,*)
      write(*,*)'             Choosen elasticity:'
      if(iuse.eq.1)then
         write(*,*)'             Isotropic, homogeneous model using'
         write(*,*)'             values set in input file.'
      else
         if(itropic.eq.0)then
            write(*,*)'             Isotropic, homogeneous model using'
            write(*,*)'             isotropic average of the following'
            write(*,*)'             anisotropic elasticity:'
         elseif(itropic.eq.1)then
            write(*,*)'             Anisotropic, homogeneous model'
            write(*,*)'             using the following elasticity'
         endif
         if(ielas.eq.-1)
     +        write(*,*)'             Elasticity from input file'
         if(ielas.eq.1)
     +        write(*,*)'             1 Crack: Shearer & Chapman'
         if(ielas.eq.2)
     +        write(*,*)'             2 Cubic nickel: Miller & Musgrave'
         if(ielas.eq.3)
     +        write(*,*)'             3 Halite: Raymer et al., 2000'
         if(ielas.eq.4)
     +        write(*,*)'             4 Halite shear: Raymer et al.'
         if(ielas.eq.5)
     +        write(*,*)'             5 Halite axial: Raymer et al.'
         if(ielas.eq.6)
     +        write(*,*)'             6 Mylonite: Lloyd & Kendall'
         if(ielas.eq.7)
     +        write(*,*)'             7 Quartz: Lloyd & Kendall'
         if(ielas.eq.8)
     +        write(*,*)'             8 Quartz shear: Lloyd & Kendall'
         if(ielas.eq.9)
     +        write(*,*)'             9 Olivine: Babuska & Cara'
         if(ielas.eq.10)
     +        write(*,*)'             10 Olivine Atype: Jung & Karato'
         if(ielas.eq.11)
     +        write(*,*)'             11 Olivine Btype: Jung & Karato'
         if(ielas.eq.12)then
            write(*,*)'             12 Elasticity from file:'
            write(*,*)modelcart_e
         endif
      endif
      endif
      write(*,*)

      do kx2=1,nx2
         do kx3=1,nx3
            if(iuse.eq.1)then                       !Iso/hom
               itropic=0
               call isotens(a)
            else
               if(ielas.eq.-1.and.iuse.eq.0)then     !Crack
                  call ec_file(a)
                  call aijkl(a,a3)
               elseif(ielas.eq.1.and.iuse.eq.0)then !Crack
                  call crack_sc(a)
                  call aijkl(a,a3)
               elseif(ielas.eq.2.and.iuse.eq.0)then !Cubic nickel
                  call nickel(a3)
               elseif(ielas.eq.3.and.iuse.eq.0)then !Halite
                  call halite(a3)
               elseif(ielas.eq.4.and.iuse.eq.0)then !Halite shear
                  call halite_shr(a3)
               elseif(ielas.eq.5.and.iuse.eq.0)then !Halite axial
                  call halite_axi(a3)
               elseif(ielas.eq.6.and.iuse.eq.0)then !Mylonite
                  call mylonite(a3)
               elseif(ielas.eq.7.and.iuse.eq.0)then !Quartz
                  call quartz(a3)
               elseif(ielas.eq.8.and.iuse.eq.0)then !Quartz shear
                  call quartz_shr(a3)
               elseif(ielas.eq.9.and.iuse.eq.0)then !Olivine
                  call olivine(a)
                  call aijkl(a,a3)
               elseif(ielas.eq.10.and.iuse.eq.0)then !Olivine
                  call olivine_a(a)
                  call aijkl(a,a3)
               elseif(ielas.eq.11.and.iuse.eq.0)then !Olivine
                  call olivine_b(a)
                  call aijkl(a,a3)
               elseif(ielas.eq.12.and.iuse.eq.0)then !Afar
                  call afar_homoec(a)
                  call aijkl(a,a3)
               endif
               call averps(a3,vfast,vslow)
               call rotate(alpha1,beta1,gamma1,a3,a3r)
               call rotate(alpha2,beta2,gamma2,a3r,a3)
               if(itropic.eq.0)then
                  call isotens(a)
               elseif(itropic.eq.1)then
                  call c66(a3,a)
               endif
            endif
            do i=1,6
               do j=1,6
                  a6(i,j,kx3,kx2,1)=a(i,j)
               enddo
            enddo
         enddo
      enddo

      return
      end

c***********************************************************************
      subroutine freqelastc(ix1,jx1,nt,dt,aw6)
c     Input elastic constants (divided by density) in Voigt notation.
c***********************************************************************

      implicit none
      include '../Input/narc_dp.par'

      integer ix1,jx1,i,j,nt,nom,iom
      real*8 omn,delom,om,dt,pi,ttot
      real*8 a3(3,3,3,3),a3r(3,3,3,3),a(6,6)
      complex*16 ai,aw6(6,6,ntmx)

      pi=4.d0*pi4

      ttot=dt*dble(nt)
      omn=pi/dt
      nom=nt/2 + 1
      delom=2.d0*pi/ttot
      ai=dcmplx(0.d0,1.d0)

      write(*,*)nt,dt
      write(*,*)pi,ttot
      write(*,*)omn,nom,delom
      write(*,*)ai
      read(*,*)

      return
      end
