c........+.........+.........+.........+.........+.........+.........+..
c  File Crystal.f contains the following 14 subroutines: Crack_sc, 
c  Nickel6, Nickel, Halite, Halite_axi, Halite_shr, Quartz, Quartz_shr, 
c  Mylonite, Olivine_A, Olivine_B, RFlayer, Afar_homoec and Ec_file.
c
c  Last modified: October 8, 2012.
c***********************************************************************
      subroutine crack_sc(as6)
c     Elasticity of hexagonally symmetric cracked medium based on 
c     Shearer & Chapman (model 4).  Note that x1-direction is the 
c     axis of rotational symmetry.  Elasticity is divided by density.
c     ichse:
c      1=thin water-filled cracks'
c      2=thick water-filled cracks'
c      3=thin dry cracks'
c      4=thin water-filled cracks- high crack density.'
c***********************************************************************
 
      implicit none
      include 'narc_dp.par'

      integer ichse,i,j
      real*8 a(6,6),as6(6,6)

c     SI units (Pa and kg/m^3)

      den=2800.d0
      ichse=4
      if(ichse.eq.1)then        ! thin water-filled cracks 
         a(1,1)=20.04d+6
         a(2,2)=20.22d+6
         a(6,6)=5.10d+6
         a(4,4)=6.38d+6
         a(1,2)=7.41d+6
      elseif(ichse.eq.2)then    ! thick water-filled cracks 
         a(1,1)=14.02d+6
         a(2,2)=19.40d+6
         a(6,6)=5.10d+6
         a(4,4)=6.38d+6
         a(1,2)=5.18d+6
      elseif(ichse.eq.3)then    ! thin dry cracks 
         a(1,1)=11.91d+6
         a(2,2)=19.11d+6
         a(6,6)=5.10d+6
         a(4,4)=6.38d+6
         a(1,2)=4.40d+6
      elseif(ichse.eq.4)then    ! thin water-filled cracks 
         a(1,1)=19.63d+6           ! (high crack density)
         a(2,2)=20.16d+6
         a(6,6)=3.48d+6
         a(4,4)=6.38d+6
         a(1,2)=7.26d+6
      endif

      do i=1,6
         do j=1,6
            as6(i,j)=0.d0
         enddo
      enddo

      as6(1,1)=a(1,1)
      as6(2,2)=a(2,2)
      as6(3,3)=a(2,2)
      as6(1,2)=a(1,2)
      as6(1,3)=a(1,2)
      as6(4,4)=a(4,4)
      as6(6,6)=a(6,6)
      as6(5,5)=a(6,6)
      as6(2,3)=as6(3,3) - 2.d0*as6(4,4)
 
      as6(2,1)=as6(1,2)           ! impose symmetry
      as6(3,1)=as6(1,3)
      as6(3,2)=as6(2,3)
 
      return
      end

c***********************************************************************
      subroutine nickel6(as6)
c     Elasticity of cubic symmetry (nickel) after Miller and Musgrave 
c     (page 358) - elasticity is divided by density.  Orientation based
c     on the (0,0,1) axis or , according to Musgrave (page 108), with 
c     reference axes coincident with the edges of the cubic cell.
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer i,j
      real*8 c(6,6),as6(6,6)

c     SI units (Pa and kg/m^3)

      den=8950.d0
      c(1,1)=25.2d10/den
      c(1,2)=15.4d10/den
      c(4,4)=12.2d10/den

      do i=1,6
         do j=1,6
            as6(i,j)=0.d0
         enddo
      enddo

      as6(1,1)=c(1,1)
      as6(2,2)=c(1,1)
      as6(3,3)=c(1,1)
      as6(1,2)=c(1,2)
      as6(1,3)=c(1,2)
      as6(2,3)=c(1,2)
      as6(4,4)=c(4,4)
      as6(5,5)=c(4,4)
      as6(6,6)=c(4,4)

c     Impose symmetry

      as6(2,1)=as6(1,2)
      as6(3,1)=as6(1,3)
      as6(3,2)=as6(2,3)

      return
      end

c***********************************************************************
      subroutine nickel(a)
c     Elasticity of cubic symmetry (nickel) after Miller/Musgrave p.358.
c     Elasticity is devided by density (SI units - Pa and kg/m^3)
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer i,j,k,l
      real*8 c(6,6),a(3,3,3,3)

      den=8950.d0
      c(1,1)=25.2d10/den
      c(1,2)=15.4d10/den
      c(4,4)=12.2d10/den

      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  a(i,j,k,l)=0.d0
               enddo
            enddo
         enddo
      enddo

      a(1,1,1,1)=c(1,1)
      a(2,2,2,2)=c(1,1)
      a(3,3,3,3)=c(1,1)
      a(1,1,2,2)=c(1,2)
      a(1,1,3,3)=c(1,2)
      a(2,2,3,3)=c(1,2)
      a(1,2,1,2)=c(4,4)
      a(1,3,1,3)=c(4,4)
      a(2,3,2,3)=c(4,4)

c     Impose symmetry

      a(2,2,1,1)=c(1,2)
      a(3,3,1,1)=c(1,2)
      a(3,3,2,2)=c(1,2)
      a(2,1,1,2)=c(4,4)
      a(1,2,2,1)=c(4,4)
      a(2,1,2,1)=c(4,4)
      a(3,1,1,3)=c(4,4)
      a(1,3,3,1)=c(4,4)
      a(3,1,3,1)=c(4,4)
      a(3,2,2,3)=c(4,4)
      a(2,3,3,2)=c(4,4)
      a(3,2,3,2)=c(4,4)

      return
      end

c***********************************************************************
      subroutine halite(a)
c     Elasticity of cubic symmetry (halite) after Raymer et al. (2000).
c     Elasticity is devided by density (SI units - Pa and kg/m^3)
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer i,j,k,l
      real*8 c(6,6),a(3,3,3,3)

      den=2160.d0
      c(1,1)=49.1d9/den
      c(1,2)=14.0d9/den
      c(4,4)=12.7d9/den

      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  a(i,j,k,l)=0.d0
               enddo
            enddo
         enddo
      enddo

      a(1,1,1,1)=c(1,1)
      a(2,2,2,2)=c(1,1)
      a(3,3,3,3)=c(1,1)
      a(1,1,2,2)=c(1,2)
      a(1,1,3,3)=c(1,2)
      a(2,2,3,3)=c(1,2)
      a(1,2,1,2)=c(4,4)
      a(1,3,1,3)=c(4,4)
      a(2,3,2,3)=c(4,4)

c     Impose symmetry

      a(2,2,1,1)=c(1,2)
      a(3,3,1,1)=c(1,2)
      a(3,3,2,2)=c(1,2)
      a(2,1,1,2)=c(4,4)
      a(1,2,2,1)=c(4,4)
      a(2,1,2,1)=c(4,4)
      a(3,1,1,3)=c(4,4)
      a(1,3,3,1)=c(4,4)
      a(3,1,3,1)=c(4,4)
      a(3,2,2,3)=c(4,4)
      a(2,3,3,2)=c(4,4)
      a(3,2,3,2)=c(4,4)

      return
      end

c***********************************************************************
      subroutine halite_axi(a)
c     Elasticity of halite polycrystal under 200% axial extension using 
c     parameters of Raymer, Tommasi & Kendall (2000, Geophysics, 65, 
c     1272-1280).  Elasticity is divided by density.
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer i,j,k,l
      real*8 c(6,6),a(3,3,3,3)

      den=2160.d0
      c(1,1)=47.3d+09
      c(1,2)=14.9d+09
      c(1,3)=14.9d+09
      c(1,4)=0.d0
      c(1,5)=-0.1d+09
      c(1,6)=0.1d+09
      c(2,2)=45.6d+09
      c(2,3)=16.6d+09
      c(2,4)=0.1d+09
      c(2,5)=0.d0
      c(2,6)=0.d0
      c(3,3)=45.7d+09
      c(3,4)=0.d0
      c(3,5)=0.1d+09
      c(3,6)=0.d0
      c(4,4)=15.0d+09
      c(4,5)=0.d0
      c(4,6)=0.d0
      c(5,5)=13.4d+09
      c(5,6)=0.d0
      c(6,6)=13.4d+09

      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  a(i,j,k,l)=0.d0
               enddo
            enddo
         enddo
      enddo

      a(1,1,1,1)=c(1,1)/den
      a(2,2,2,2)=c(2,2)/den
      a(3,3,3,3)=c(3,3)/den
      a(1,1,2,2)=c(1,2)/den
      a(1,1,3,3)=c(1,3)/den
      a(2,2,3,3)=c(2,3)/den
      a(1,2,1,2)=c(6,6)/den
      a(1,3,1,3)=c(5,5)/den
      a(2,3,2,3)=c(4,4)/den
      a(1,1,2,3)=c(1,4)/den
      a(1,1,1,3)=c(1,5)/den
      a(1,1,1,2)=c(1,6)/den
      a(2,2,2,3)=c(2,4)/den
      a(2,2,1,3)=c(2,5)/den
      a(2,2,1,2)=c(2,6)/den
      a(3,3,2,3)=c(3,4)/den
      a(3,3,1,3)=c(3,5)/den
      a(3,3,1,2)=c(3,6)/den
      a(2,3,1,3)=c(4,5)/den
      a(2,3,1,2)=c(4,6)/den
      a(1,3,1,2)=c(5,6)/den

c     Impose symmetry -- general form

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

      return
      end

c***********************************************************************
      subroutine halite_shr(a)
c     Elasticity of halite polycrystal under 600% simple shear using 
c     parameters of Raymer, Tommasi & Kendall (2000, Geophysics, 65, 
c     1272-1280).  (Shear strain factor 10.)  Elasticity is divided by
c     density.
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer i,j,k,l
      real*8 c(6,6),a(3,3,3,3)

      den=2160.d0
      c(1,1)=45.7d+09
      c(1,2)=14.5d+09
      c(1,3)=16.9d+09
      c(1,4)=0.d0
      c(1,5)=-0.1d+09
      c(1,6)=0.1d+09
      c(2,2)=48.0d+09
      c(2,3)=14.6d+09
      c(2,4)=0.d0
      c(2,5)=0.d0
      c(2,6)=-0.2d+09
      c(3,3)=45.6d+09
      c(3,4)=0.d0
      c(3,5)=0.d0
      c(3,6)=0.d0
      c(4,4)=13.1d+09
      c(4,5)=0.d0
      c(4,6)=0.d0
      c(5,5)=15.3d+09
      c(5,6)=0.d0
      c(6,6)=13.1d+09

      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  a(i,j,k,l)=0.d0
               enddo
            enddo
         enddo
      enddo

      a(1,1,1,1)=c(1,1)/den
      a(2,2,2,2)=c(2,2)/den
      a(3,3,3,3)=c(3,3)/den
      a(1,1,2,2)=c(1,2)/den
      a(1,1,3,3)=c(1,3)/den
      a(2,2,3,3)=c(2,3)/den
      a(1,2,1,2)=c(6,6)/den
      a(1,3,1,3)=c(5,5)/den
      a(2,3,2,3)=c(4,4)/den
      a(1,1,2,3)=c(1,4)/den
      a(1,1,1,3)=c(1,5)/den
      a(1,1,1,2)=c(1,6)/den
      a(2,2,2,3)=c(2,4)/den
      a(2,2,1,3)=c(2,5)/den
      a(2,2,1,2)=c(2,6)/den
      a(3,3,2,3)=c(3,4)/den
      a(3,3,1,3)=c(3,5)/den
      a(3,3,1,2)=c(3,6)/den
      a(2,3,1,3)=c(4,5)/den
      a(2,3,1,2)=c(4,6)/den
      a(1,3,1,2)=c(5,6)/den

c     Impose symmetry -- general form

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

      return
      end

c***********************************************************************
      subroutine quartz(a)
c     Elasticity of quartz.  SI units (pa and kg/m^3).
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer i,j,k,l
      real*8 c(6,6),a(3,3,3,3)

      c(1,1)=3.26241d+07
      c(1,2)=2.64662d+06
      c(1,3)=4.47368d+06
      c(1,4)=-6.78196d+06
      c(1,5)=-3759.40d0
      c(1,6)=-3759.40d0
      c(2,2)=3.26391d+07
      c(2,3)=4.48120d+06
      c(2,4)=6.77820d+06
      c(2,5)=0.d0
      c(2,6)=0.d0
      c(3,3)=3.97556d+07
      c(3,4)=-3759.40d0
      c(3,5)=0.d0
      c(3,6)=0.d0
      c(4,4)=2.18835d+07
      c(4,5)=-3759.40d0
      c(4,6)=0.d0
      c(5,5)=2.18759d+07
      c(5,6)=-6.78196d+06
      c(6,6)=1.49925d+07
      den=2660.d0

      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  a(i,j,k,l)=0.d0
               enddo
            enddo
         enddo
      enddo

      a(1,1,1,1)=c(1,1)
      a(2,2,2,2)=c(2,2)
      a(3,3,3,3)=c(3,3)
      a(1,1,2,2)=c(1,2)
      a(1,1,3,3)=c(1,3)
      a(2,2,3,3)=c(2,3)
      a(1,2,1,2)=c(6,6)
      a(1,3,1,3)=c(5,5)
      a(2,3,2,3)=c(4,4)
      a(1,1,2,3)=c(1,4)
      a(1,1,1,3)=c(1,5)
      a(1,1,1,2)=c(1,6)
      a(2,2,2,3)=c(2,4)
      a(2,2,1,3)=c(2,5)
      a(2,2,1,2)=c(2,6)
      a(3,3,2,3)=c(3,4)
      a(3,3,1,3)=c(3,5)
      a(3,3,1,2)=c(3,6)
      a(2,3,1,3)=c(4,5)
      a(2,3,1,2)=c(4,6)
      a(1,3,1,2)=c(5,6)

c     Impose symmetry -- general form

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

      return
      end

c***********************************************************************
      subroutine quartz_shr(a)
c     Elasticity of quartz shear zone.  SI units (pa and kg/m^3).
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer i,j,k,l
      real*8 c(6,6),a(3,3,3,3)

      c(1,1)=3.62445d+07
      c(1,2)=2.97284d+06
      c(1,3)=4.52914d+06
      c(1,4)=-1.15212d+06
      c(1,5)=872587.d0
      c(1,6)=698825.d0
      c(2,2)=3.59007d+07
      c(2,3)=3.30525d+06
      c(2,4)=1.74140d+06
      c(2,5)=1.33343d+06
      c(2,6)=-804593.d0
      c(3,3)=3.44124d+07
      c(3,4)=1.25411d+06
      c(3,5)=-1.59030d+06
      c(3,6)=967023.d0
      c(4,4)=1.65640d+07
      c(4,5)=1.48831d+06
      c(4,6)=1.52608d+06
      c(5,5)=1.90383d+07
      c(5,6)=-3777.43d0
      c(6,6)=1.71949d+07
      den=2647.30d0

      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  a(i,j,k,l)=0.d0
               enddo
            enddo
         enddo
      enddo

      a(1,1,1,1)=c(1,1)
      a(2,2,2,2)=c(2,2)
      a(3,3,3,3)=c(3,3)
      a(1,1,2,2)=c(1,2)
      a(1,1,3,3)=c(1,3)
      a(2,2,3,3)=c(2,3)
      a(1,2,1,2)=c(6,6)
      a(1,3,1,3)=c(5,5)
      a(2,3,2,3)=c(4,4)
      a(1,1,2,3)=c(1,4)
      a(1,1,1,3)=c(1,5)
      a(1,1,1,2)=c(1,6)
      a(2,2,2,3)=c(2,4)
      a(2,2,1,3)=c(2,5)
      a(2,2,1,2)=c(2,6)
      a(3,3,2,3)=c(3,4)
      a(3,3,1,3)=c(3,5)
      a(3,3,1,2)=c(3,6)
      a(2,3,1,3)=c(4,5)
      a(2,3,1,2)=c(4,6)
      a(1,3,1,2)=c(5,6)

c     Impose symmetry -- general form

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

      return
      end

c***********************************************************************
      subroutine mylonite(a)
c     Elasticity of quartz mylonite zone.  SI units (pa and kg/m^3).
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer i,j,k,l
      real*8 c(6,6),a(3,3,3,3)

      c(1,1)=3.97386d+07
      c(1,2)=3.81898d+06
      c(1,3)=3.04083d+06
      c(1,4)=-400408.d0
      c(1,5)=192649.d0
      c(1,6)=79326.1d0
      c(2,2)=3.21611d+07
      c(2,3)=4.30627d+06
      c(2,4)=309750.d0
      c(2,5)=389076.d0
      c(2,6)=86881.0d0
      c(3,3)=3.41858d+07
      c(3,4)=30219.5d0
      c(3,5)=-1.95293d+06
      c(3,6)=169984.d0
      c(4,4)=1.58350d+07
      c(4,5)=355079.d0
      c(4,6)=-430627.d0
      c(5,5)=1.85812d+07
      c(5,6)=-441960.d0
      c(6,6)=1.84188d+07
      den=2647.30d0

      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  a(i,j,k,l)=0.d0
               enddo
            enddo
         enddo
      enddo

      a(1,1,1,1)=c(1,1)
      a(2,2,2,2)=c(2,2)
      a(3,3,3,3)=c(3,3)
      a(1,1,2,2)=c(1,2)
      a(1,1,3,3)=c(1,3)
      a(2,2,3,3)=c(2,3)
      a(1,2,1,2)=c(6,6)
      a(1,3,1,3)=c(5,5)
      a(2,3,2,3)=c(4,4)
      a(1,1,2,3)=c(1,4)
      a(1,1,1,3)=c(1,5)
      a(1,1,1,2)=c(1,6)
      a(2,2,2,3)=c(2,4)
      a(2,2,1,3)=c(2,5)
      a(2,2,1,2)=c(2,6)
      a(3,3,2,3)=c(3,4)
      a(3,3,1,3)=c(3,5)
      a(3,3,1,2)=c(3,6)
      a(2,3,1,3)=c(4,5)
      a(2,3,1,2)=c(4,6)
      a(1,3,1,2)=c(5,6)

c     Impose symmetry -- general form

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

      return
      end

c***********************************************************************
      subroutine olivine(as6)
c     Elasticity of orthorombic symmetry (olivine) after Babuska and 
c     Cara (page 49) - elasticity is divided by density.  Units are 
c     in pa and kg/m^3.  Orientation ??
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer i,j
      real*8 as6(6,6)

c     SI units (Pa and kg/m^3)

      den=3311.d0
      do i=1,6
         do j=1,6
            as6(i,j)=0.d0
         enddo
      enddo

      as6(1,1)=323.7d9/den
      as6(2,2)=197.6d9/den
      as6(3,3)=235.1d9/den
      as6(1,2)=66.4d9/den
      as6(1,3)=71.6d9/den
      as6(2,3)=75.6d9/den
      as6(4,4)=64.6d9/den
      as6(5,5)=78.7d9/den
      as6(6,6)=79.d9/den

c     Impose symmetry

      as6(2,1)=as6(1,2)
      as6(3,1)=as6(1,3)
      as6(3,2)=as6(2,3)

      return
      end

c***********************************************************************
      subroutine olivine_a(as6)
c     Elasticity of type A olivine from Jung & Karato (see Science, 
c     August 24, 2001) - elasticity is divided by density.  Units are 
c     in pa and kg/m^3.  
c     Orientation: 1: shear direction, 2: normal to shear plane and 
c     3:orthogonal to 1 and 2
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer i,j
      real*8 as6(6,6)

c     SI units (Pa and kg/m^3)

      den=3355.d0

      do i=1,6
         do j=1,6
            as6(i,j)=0.d0
         enddo
      enddo

      as6(1,1)= 2.5057d11/den
      as6(1,2)= 0.7590d11/den
      as6(1,3)= 0.7879d11/den
      as6(1,4)= 0.0032d11/den
      as6(1,5)=-0.0030d11/den
      as6(1,6)=-0.0351d11/den
      as6(2,2)= 2.2103d11/den
      as6(2,3)= 0.7703d11/den
      as6(2,4)=-0.0129d11/den
      as6(2,5)= 0.0004d11/den
      as6(2,6)=-0.0639d11/den
      as6(3,3)= 2.3213d11/den
      as6(3,4)=-0.0194d11/den
      as6(3,5)=-0.0096d11/den
      as6(3,6)=-0.0127d11/den
      as6(4,4)= 0.7740d11/den
      as6(4,5)=-0.0229d11/den
      as6(4,6)= 0.0011d11/den
      as6(5,5)= 0.8049d11/den
      as6(5,6)=-0.0066d11/den
      as6(6,6)= 0.7839d11/den

c     Impose symmetry

      as6(2,1)=as6(1,2)
      as6(3,1)=as6(1,3)
      as6(4,1)=as6(1,4)
      as6(5,1)=as6(1,5)
      as6(6,1)=as6(1,6)
      as6(3,2)=as6(2,3)
      as6(4,2)=as6(2,4)
      as6(5,2)=as6(2,5)
      as6(6,2)=as6(2,6)
      as6(4,3)=as6(3,4)
      as6(5,3)=as6(3,5)
      as6(6,3)=as6(3,6)
      as6(5,4)=as6(4,5)
      as6(6,4)=as6(4,6)
      as6(6,5)=as6(5,6)

      return
      end

c***********************************************************************
      subroutine olivine_b(as6)
c     Elasticity of type B olivine from Jung & Karato (see Science, 
c     August 24, 2001) - elasticity is divided by density.  Units are 
c     in pa and kg/m^3.  
c     Orientation: 1: shear direction, 2: normal to shear plane and 
c     3:orthogonal to 1 and 2
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer i,j
      real*8 as6(6,6)

c     SI units (Pa and kg/m^3)

      den=3355.d0

      do i=1,6
         do j=1,6
            as6(i,j)=0.d0
         enddo
      enddo

      as6(1,1)= 2.3339d11/den
      as6(1,2)= 0.7780d11/den
      as6(1,3)= 0.7600d11/den
      as6(1,4)=-0.0084d11/den
      as6(1,5)= 0.0151d11/den
      as6(1,6)= 0.0223d11/den
      as6(2,2)= 2.1521d11/den
      as6(2,3)= 0.7478d11/den
      as6(2,4)= 0.0138d11/den
      as6(2,5)=-0.0114d11/den
      as6(2,6)= 0.0169d11/den
      as6(3,3)= 2.6226d11/den
      as6(3,4)= 0.0300d11/den
      as6(3,5)= 0.0561d11/den
      as6(3,6)= 0.0111d11/den
      as6(4,4)= 0.7898d11/den
      as6(4,5)= 0.0048d11/den
      as6(4,6)= 0.0078d11/den
      as6(5,5)= 0.7928d11/den
      as6(5,6)=-0.0032d11/den
      as6(6,6)= 0.7224d11/den

c     Impose symmetry

      as6(2,1)=as6(1,2)
      as6(3,1)=as6(1,3)
      as6(4,1)=as6(1,4)
      as6(5,1)=as6(1,5)
      as6(6,1)=as6(1,6)
      as6(3,2)=as6(2,3)
      as6(4,2)=as6(2,4)
      as6(5,2)=as6(2,5)
      as6(6,2)=as6(2,6)
      as6(4,3)=as6(3,4)
      as6(5,3)=as6(3,5)
      as6(6,3)=as6(3,6)
      as6(5,4)=as6(4,5)
      as6(6,4)=as6(4,6)
      as6(6,5)=as6(5,6)

      return
      end

c***********************************************************************
      subroutine rflayer(a,layer)
c     Elasticity of cubic symmetry (nickel) after Miller/Musgrave p.358.
c     Elasticity is devided by density (SI units - Pa and kg/m^3)
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer i,j,layer
      real*8 as6(6,6),a(3,3,3,3)

      do i=1,6
         do j=1,6
            as6(i,j)=0.d0
         enddo
      enddo

      if(layer.eq.1.or.layer.eq.3.or.layer.eq.5)then
         itropic=0
         den=2700.d0
         vfast=8000.d0
         vslow=4500.d0
         call isotens(as6)
      elseif(layer.eq.2)then
         itropic=1
         den=2700.d0
         as6(1,1)=6.784000e+07
         as6(1,2)=2.506450e+07
         as6(1,3)=2.522512e+07
         as6(1,4)=2.206507e+05
         as6(1,5)=0.000000e+00
         as6(1,6)=0.000000e+00
         as6(2,2)=6.655885e+07
         as6(2,3)=2.539745e+07
         as6(2,4)=-9.085934e+05
         as6(2,5)=0.000000e+00
         as6(2,6)=0.000000e+00
         as6(3,3)=6.524549e+07
         as6(3,4)=-8.956165e+05
         as6(3,5)=0.000000e+00
         as6(3,6)=0.000000e+00
         as6(4,4)=2.026783e+07
         as6(4,5)=0.000000e+00
         as6(4,6)=0.000000e+00
         as6(5,5)=2.064972e+07
         as6(5,6)=-5.708633e+05
         as6(6,6)=2.106528e+07
      elseif(layer.eq.4)then
         itropic=1
         den=2700.d0
         as6(1,1)=6.784000e+07
         as6(1,2)=2.532469e+07
         as6(1,3)=2.496494e+07
         as6(1,4)=1.509340e+05
         as6(1,5)=0.000000e+00
         as6(1,6)=0.000000e+00
         as6(2,2)=6.444085e+07
         as6(2,3)=2.538796e+07
         as6(2,4)=-6.071352e+05
         as6(2,5)=0.000000e+00
         as6(2,6)=0.000000e+00
         as6(3,3)=6.738246e+07
         as6(3,4)=-6.270170e+05
         as6(3,5)=0.000000e+00
         as6(3,6)=0.000000e+00
         as6(4,4)=2.025834e+07
         as6(4,5)=0.000000e+00
         as6(4,6)=0.000000e+00
         as6(5,5)=2.132287e+07
         as6(5,6)=-3.904935e+05
         as6(6,6)=2.039213e+07
      endif

      call aijkl(as6,a)

      return
      end

c***********************************************************************
      subroutine afar_homoec(a)
c     Elasticity from files provided by J-M-K for the Afar splitting 
c     study.  
c     Elasticity is devided by density (SI units - Pa and kg/m^3)
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer i,ii,j,jj
      real*8 den1,den2
      real*8 a1(6,6),a2(6,6),a(6,6)
      character*100 modname2

      open(30,file=modelcart_e)
      do i=1,6
         do j=i,6
            read(30,*)ii,jj,a1(i,j)
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
            if(i.ne.j)a2(j,i)=a2(i,j)
         enddo
      enddo
      read(30,*)ii,jj,den2
      close(30)
      den=(den1+den2)/2.d0
      do i=1,6
         do j=1,6
            a(i,j)=(a1(i,j)+a2(i,j))/2.d0
         enddo
      enddo

      return
      end

c***********************************************************************
      subroutine karakoram_homoec(a)
c     Elasticity from files provided by Lloyd for the Karakoram 
c     splitting study.  
c     Elasticity is devided by density (SI units - Pa and kg/m^3)
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer i,ii,j,jj
      real*8 den1,den2
      real*8 a1(6,6),a2(6,6),a(6,6)
      character*100 modname2

      open(30,file=modelcart_e)
      do i=1,6
         do j=i,6
            read(30,*)ii,jj,a1(i,j)
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
            if(i.ne.j)a2(j,i)=a2(i,j)
         enddo
      enddo
      read(30,*)ii,jj,den2
      close(30)
      den=(den1+den2)/2.d0
      do i=1,6
         do j=1,6
            a1(i,j)=a1(i,j)/den1
            a2(i,j)=a2(i,j)/den2
            a(i,j)=(a1(i,j)+a2(i,j))/2.d0
         enddo
      enddo

      return
      end

c***********************************************************************
      subroutine ec_file(a6)
c     Elasticity taken from input file specified in narc.inp
c     Elasticity is devided by density (SI units - Pa and kg/m^3)
c***********************************************************************

      implicit none
      include 'narc_dp.par'

      integer i,ii,j,jj
      real*8 angle1,angle2,angle3,pi
      real*8 a6(6,6),a3(3,3,3,3),a3r(3,3,3,3)

      pi=4.d0*pi4

      open(30,file=modelcart_e)
      do i=1,6
         do j=i,6
            read(30,*)ii,jj,a6(i,j)
            if(i.ne.j)a6(j,i)=a6(i,j)
         enddo
      enddo
      read(30,*)ii,jj,den
      close(30)

c     Rotate from Atrak coorsdinates to Oneway
      call aijkl(a6,a3)
      angle1=pi/2.d0
      angle2=0.d0
      angle3=0.d0
      call rotate(angle1,angle2,angle3,a3,a3r)
      call c66(a3r,a6)

      return
      end
