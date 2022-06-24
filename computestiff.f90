! File computestiff.f90 contains the following subroutine: computestiffness.
!
! Created by Alan Baird.
! Last modified: April 18, 2012.
!***************************************************************************************************
subroutine computestiffness(vp,vs,rho,poro,tau,bulk,fracden,fracrad,vtieps,vtigam,vtidel,freq,c6)
! Computes the frequency dependent complex stiffness tensor (c6) at a given frequency
! using the Chapman (2003) squirt-flow theory with a VTI background elasticity.
! Fractures are assumed to be vertical and oriented normal to the x1 axis.

  double precision, intent(in) :: vp,vs,rho,poro,tau,bulk,fracden,fracrad,vtieps,vtigam,vtidel,freq
  complex, intent(out) :: c6(6,6)
  double precision :: zero,ft,po,gam

  !Definitions
  zero=0.d0
  ft=4.d0/3.d0
  po=0.5d0*(vp*vp-2.d0*vs*vs)/(vp*vp-vs*vs)
  gam = (1.18d0*(1.d0+((ft*vs*vs*rho)/bulk)))/(1.d0-po)  

  !Compute VTI background elasticity
  call cvtic(vp,vs,rho,vtieps,vtigam,vtidel,c6)
  !Include poroelasticity and squirtflow
  call computeFDstiff(c6,rho,fracden,fracrad,zero,poro,gam,tau,freq)
	
end subroutine computestiffness
