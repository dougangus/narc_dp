! File VTI_lib.f90 contains the following 3 subroutines: cvtic, zeros2d and symmetry2d.
!
! Created by Alan Baird.
! Last modified: April 17, 2012.
! From here - functions for generating stiffness tensors based on RP models
!***************************************************************************************************
subroutine cvtic(vp,vs,rho,eps,gam,del,c)
!---------------------------------------------------------------------------------------------------
!
! Function CVTI generates a vti stiffness tensor based on input P and S velocities and
! Thomsen's parameters
!
! Input parameters are:
!    vp = P-wave velocity
!    vs = S-wave velocity
!    eps = Thomsen's epsilon
!    gam = Thomsen's gamma
!    del = Thomsen's delta
!
! Calls to other subroutines:
!    zeros2d - initialise matrix full of zeros
!    symmetry2d - to fill matrix
!
! Written by J.P. Verdon, University of Bristol, 2008-2011
!
!***************************************************************************************************

implicit none

DOUBLE PRECISION, INTENT(IN)           :: vp,vs,rho,eps,gam,del
DOUBLE PRECISION           :: cvti(6,6)
COMPLEX, INTENT(OUT)      :: c(6,6)

call zeros2d(cvti,6,6)
cvti(3,3)=vp*vp*rho
cvti(4,4)=vs*vs*rho
cvti(1,1)=cvti(3,3)*(2*eps+1)
cvti(6,6)=cvti(4,4)*(2*gam+1)     
cvti(2,2)=cvti(1,1)
cvti(5,5)=cvti(4,4)
cvti(1,2)=cvti(1,1)-2*cvti(6,6) 
cvti(1,3)=sqrt(2*del*cvti(3,3)*(cvti(3,3)-cvti(4,4)) + (cvti(3,3)-cvti(4,4))*(cvti(3,3)-cvti(4,4))) - cvti(4,4)     
if (2*del*cvti(3,3)*(cvti(3,3)-cvti(4,4)) + (cvti(3,3)-cvti(4,4))*(cvti(3,3)-cvti(4,4)).lt.0)then
    write(*,*)'stopping program, c(1,3) is not viable'
    read(*,*)
    stop
endif
cvti(2,3)=cvti(1,3)
call symmetry2d(cvti,6)   


 c=cmplx(cvti)   

return
end subroutine cvtic
!***************************************************************************************************


!***************************************************************************************************
  subroutine zeros2d(a,n1,n2)
!***************************************************************************************************
implicit none
integer             :: i,j,n1,n2
DOUBLE PRECISION                :: a(n1,n2)
      
do i=1,n1
    do j=1,n2
        a(i,j)=0.0
    enddo
enddo
return
end subroutine zeros2d
!***************************************************************************************************


!***************************************************************************************************
  subroutine symmetry2d(a,n)
!***************************************************************************************************
implicit none
integer         :: n,i,j
DOUBLE PRECISION           :: a(n,n)

do i=1,n
    do j=i,n
        if(i.ne.j)a(j,i)=a(i,j)
    enddo
enddo
return
end subroutine symmetry2d
!***************************************************************************************************




