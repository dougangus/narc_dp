c...+....|....+....|....+....|....+....|....+....|....+....|....+....|..
c  File Squirt_lib.f contains the following 2 subroutines: 
c     computeFDstiff and fdcv01_background.
c
c  Last modified: April 17, 2012.
c***********************************************************************
      subroutine computeFDstiff(c,rho,fd,af,cd,ppor,gam,tau,w)
c***********************************************************************

      implicit none
      real*8 lam,mu,rho,fd
      real*8 af,cd,ppor,gam,tau,w,zero
      complex c(6,6)
      
      zero=0
      mu=(c(4,4)+c(5,5)+c(6,6))/3.0
      lam=(c(1,1)+c(2,2)+c(3,3))/3.0-2.0*mu
      call fdcv01_background(lam,mu,rho,fd,af,cd,ppor,gam,tau,w,c) 

      return
      end

c***********************************************************************
      subroutine fdcv01_background(lam,mu,rho,fd,af,cd,ppor,gam,tau,w,
     +     cin) 
c     Takes real values for the lamee constants (lam & mu), the density 
c     (rho) fracture density (fd), fracture length in meters (af), the 
c     microcrack density (cd), the porosity (ppor) the parameter gamma 
c     (gam) related to fluid compressibility, the time constant (tau) 
c     related to the inverse of the squirt flow frequency and the wave 
c     frequency (w). The output is a complex elastic tensor, with the 
c     x_3 axis assumed to be the axis of symmetry. This elastic tensor 
c     is characterised by the constants c11, c33, c12, c13 and c44.
c***********************************************************************

      implicit none
      real*8 lam,mu,rho,ppor,cpor,fpor,beta,io
      real*8 tau,tauf,gam,gamp,kc,ar,w,v,pi,crit,l2,l3,l4
      real*8 af
      real*8 cd,fd
      complex i,c11,c22,c33,c13,c12,c23,c44,c55,c66
      complex kap,cin(6,6)
      complex a,b,c,z1,z2,z3,z4,z5
      complex d1,d2,f1,f2,g1,g2,g3,vel
      complex aa1,aa2,aa3,aa4,aa5,aa6,bb2,bb3,bb4,bb5
      complex jj1,jj2,jj3,jj4,kk1,kk2,kk3,kk4,kk5,kk6
      complex mm2,mm4,mm5,mm6
      complex x1,x2,x3,v1,v2,v3

      i=(0,1)
c     parameter values
      ar=1E-3
      gamp=1
      tauf=5000*af*tau
      kc=0
      cpor=4.187*cd*ar
      fpor=4.187*fd*ar
c     some definitions                                                                                          
      pi=3.1415926535
      l2=(lam*lam) + (1.33333333*lam*mu) +(0.8*mu*mu)
      l3=4*(lam*lam + (1.33333*lam*mu) + (0.5333333*mu*mu) )
      l4=lam*lam+(4.0/3.0)*lam*mu+(4.0/15.0)*mu*mu !Eqn. 58
      v=lam/(2*(lam+mu))
      kap=lam+ (0.66666666*mu)
      crit=(pi*mu*ar)/(2-2*v)
      io=(cpor/ar)/( (cpor/ar) + ppor)
      beta=(fpor/ar)/((cpor/ar) + ppor)
c     calculate d1 and d2, f1 and f2
      z1=(io/(3 + 3*kc)) + (gamp - gamp*io)
      z2=((i*w*tau)*( 1/(3 + 3*kc) -gamp))/(1+i*w*tau)
      z3=io +( (io*beta)/(1+i*w*tauf))
      a=z1 - (z2*z3)
      b=beta/( (1+kc)*(1+i*w*tauf) )
      z4=((1-io)*gam) + (( (1-io)*beta)/(1+i*w*tauf))
      z5=(io + ((io*beta)/(1 + i*w*tauf)))*
     +     ((1+i*w*gam*tau)/(1+i*w*tau))       
      c= z4 + z5
      d1=a/c
      d2=b/c
      x1=io*d1*(1+i*w*gam*tau)/(1+i*w*tau)
      x2= d1*(1-io)
      x3=(((1/(3 + 3*kc)) -gamp)*io*i*w*tau)/(1+i*w*tau)
      f1=(x1 + x2 + x3)/(1 + i*w*tauf)
      v1=(i*w*tauf)/(1 + kc)
      v2=(io*d2*(1+i*w*gam*tau))/(1+i*w*tau)
      v3=d2*(1-io)
      f2=(v1+v2+v3)/(1+i*w*tauf)
c     end definiton of d1, d2, f1, f2
c     define g1, g2, g3      
      g1=(i*w*tau)/( (1+i*w*tau)*(1+kc))
      g2=(d1 +(d1*i*w*tau*gam) - (gamp*i*w*tau))
     +     /(1+i*w*tau)
      g3=(d2*(1+i*w*gam*tau))/(1+i*w*tau)
c     end of definition of g1, g2, g3
c     calculation of c_11
      aa1=(l2/crit) + ((32*(1-v)*mu)/(15*(2-v)*pi*ar))
      aa2=(((l2/crit)+kap)*g1)+(( 
     +     ((3*kap*kap)/crit) +3*kap)*g2)+
     +     ( (((lam*kap)/crit)+lam)*g3)
c     crack correction cpor*(aa1-aa2)
      aa3=( (3*lam*lam + 4*lam*mu+ 
     +     (mu*mu*(36+20*v)/(7-5*v)) )*3*(1-v))/(4*(1+v)*mu)
      aa4=(1+(3*kap)/(4*mu))*(3*kap*d1 + lam*d2)
c     porosity correction ppor*(aa3-aa4)
      aa5=lam*lam/crit
      aa6=((3*lam*kap/crit)+(3*kap))*f1 +
     +     ( ( (lam*lam/crit) + lam)*f2)

c     fracture correction fpor*(aa5-aa6)
c      c11= (lam+2*mu) + ((-1)*cpor*(aa1-aa2)) + 
c     +     ((-1)*ppor*(aa3-aa4)) + ((-1)*fpor*(aa5-aa6))
      c33= cin(3,3) + ((-1)*cpor*(aa1-aa2)) + 
     + ((-1)*ppor*(aa3-aa4)) + ((-1)*fpor*(aa5-aa6))
      c22= cin(2,2) + ((-1)*cpor*(aa1-aa2)) + 
     +     ((-1)*ppor*(aa3-aa4)) + ((-1)*fpor*(aa5-aa6))
c     calculation of c33
      bb2=(((l2/crit)+kap)*g1)+
     +     (( ((3*kap*kap)/crit) +3*kap)*g2)+(
     +     ((((lam+2*mu)*kap)/crit)+(lam+2*mu))*g3)
c     crack correction cpor*(aa1-bb2)
      bb3=(1+((3*kap)/(4*mu)))*
     +     ( (3*kap*d1) + ((lam+2*mu)*d2) )
c     porosity correction ppor*(aa3-bb3)
      bb4=(lam +2*mu)*(lam+2*mu)/crit
      bb5=(( (3*kap*(lam+2*mu)/crit) +3*kap)*f1) +
     +     f2*( ((lam+2*mu)*
     +     (lam+2*mu)/crit) + lam + 2*mu)
c     fracture correction fpor*(bb4-bb5)
c      c33=(lam+2*mu) + ((-1)*cpor*(aa1-bb2)) + 
c     +     ((-1)*ppor*(aa3-bb3))  
c     +     +((-1)*fpor*(bb4-bb5))
      c11=cin(1,1)+ ((-1)*cpor*(aa1-bb2)) + 
     +     ((-1)*ppor*(aa3-bb3))  
     +     +((-1)*fpor*(bb4-bb5))
c     calculation of c44
      jj1=(0.2666666666667)*mu*mu*(1-g1)/crit
      jj2=((1.6)*(1-v)*mu)/(pi*ar*(2-v))
c     crack correction cpor*(jj1 + jj2)
      jj3=(15*(1-v)*mu)/(7-5*v)
c     porosity correction ppor*jj3
      vel=ppor*jj3
      jj4=(4*(1-v)*mu)/(pi*(2-v)*ar)
c     fracture correction fpor*jj4
      vel=fpor*jj4
c      c44= mu + ((-1)*cpor*(jj1+jj2)) + 
c     +     ((-1)*ppor*jj3) + ((-1)*fpor*jj4)
      c66= cin(6,6) + ((-1)*cpor*(jj1+jj2)) + 
     +     ((-1)*ppor*jj3) + ((-1)*fpor*jj4)
      c55= cin(5,5) + ((-1)*cpor*(jj1+jj2)) + 
     +     ((-1)*ppor*jj3) + ((-1)*fpor*jj4)
c      c66= mu + ((-1)*cpor*(jj1+jj2)) + 
c     +     ((-1)*ppor*jj3) + ((-1)*fpor*jj4)
c     Calculation of C12 (CHECKED AND FINE) - leads to Eqn 56 (modified)
      kk1=l4/crit-(16.0*(1.0-v)*mu)/(15.0*(2.0-v)*pi*ar)
      kk2=(l4/crit+kap)*g1+(3.0*kap*kap/crit+3.0*kap)*g2+
     +     (lam*kap/crit+lam)*g3
c     Crack correction cpor*(kk1-kk2)
      kk3=(3.0*(1.0-v)/(4.0*mu*(1.0+v)))*
     +     (3.0*lam*lam+4.0*lam*mu-mu*mu*4.0*(1.0+5.0*v)/(7.0-5.0*v))
      kk4=(1.0+3.0*kap/(4.0*mu))*(3.0*kap*d1+lam*d2)
c     Porosity correction ppor*(kk3-kk4)
      kk5=lam*lam/crit
      kk6=f1*3.0*kap*(1.0+lam/crit)+f2*lam*(1.0+lam/crit)
c     Fracture correction fpor*(kk5-kk6)
      c23=cin(2,3)-cpor*(kk1-kk2)-ppor*(kk3-kk4)-fpor*(kk5-kk6)
c     Calculation of C13 (CHECKED AND FINE) - leads to Eqn 60 (modified)
      mm2=(l4/crit+kap)*g1+3.0*kap*(1.0+kap/crit)*g2+
     +     (lam+mu)*(1.0+kap/crit)*g3
c     Crack correction cpor*(kk1-mm2)
      mm4=(1.0+3.0*kap/(4.0*mu))*(3.0*kap*d1+(lam+mu)*d2)
c     Porosity correction ppor*(kk3-mm4)
      mm5=lam*(lam+2.0*mu)/crit
      mm6=f1*3.0*kap*(1.0+(lam+mu)/crit)+f2*(lam+mu+lam*(lam+mu)/crit)
c     Fracture correction fpor*(mm5-mm6)
      c13=cin(1,3)-cpor*(kk1-mm2)-ppor*(kk3-mm4)-fpor*(mm5-mm6)
      c12=cin(1,2)-cpor*(kk1-mm2)-ppor*(kk3-mm4)-fpor*(mm5-mm6)
      c44= cin(4,4) + 0.5*(( ((-1)*cpor*(aa1-aa2)) +
     +     ((-1)*ppor*(aa3-aa4)) + ((-1)*fpor*(aa5-aa6))) - 
     +     ((-1)*cpor*(kk1-kk2)-ppor*(kk3-kk4)-fpor*(kk5-kk6)))
      cin(1,1)=c11
      cin(2,2)=c22
      cin(3,3)=c33
      cin(4,4)=c44
      cin(5,5)=c55
      cin(6,6)=c66
      cin(1,2)=c12
      cin(1,3)=c13
      cin(2,3)=c23
      cin(2,1)=c12
      cin(3,1)=c13
      cin(3,2)=c23

      return
      end 
