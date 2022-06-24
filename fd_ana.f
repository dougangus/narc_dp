c...+....|....+....|....+....|....+....|....+....|....+....|....+....|..
c  File FD_Ana.f contains the following subroutine: Fd_ana, Sortslwc and
c  Trnsfv6.  This subroutine evaluates the stability and dispersion 
c  characteristics of the FD extrapolator given the input initial 
c  conditions and grid paramenters (see Angus & Thomson, 2006).
c
c  Last modified: April 16, 2012.
c***********************************************************************
      subroutine fd_ana(nx1,dx1,nx2,dx2,nx3,dx3,nt,dt,omwin,rpeak,
     +     p1inc,p2inc,p3inc,a6,iref)
c***********************************************************************

      implicit none
      include '../Input/narc_dp.par'

      integer nx1,nx2,nx3,nt,iref,iom,nom,nom2,icomp,ix2,ix3
      integer i,j,ic,jc,kc,ierr,ievaltag,kx2,kx3,ik,jk
      real*8 dx1,dx2,dx3,dt,omwin,rpeak,pi,om,omn,tott,delom,dx
      real*8 Dctrace,Ddtrace,rkctemp,rkdtemp,dkc,rtemp,rktmp
      real*8 p1qP,p1qS1,p1qS2,rk1,rk2,rk3,rhomax,rmax,Al2norm
      real*8 a3(3,3,3,3),c11inv(3,3),a6(6,6),rkval(3)
      real*8 p0(3,3),p2(3,3),p3(3,3),p22(3,3),p23(3,3),p32(3,3),p33(3,3)
      real*8 p1inc(3,nx2mx,nx3mx),p2inc(3,nx2mx,nx3mx),
     +     p3inc(3,nx2mx,nx3mx)
      complex*16 cp0(3,3),cp2(3,3),cp3(3,3),cp22(3,3),cp23(3,3),
     +     cp32(3,3),cp33(3,3)
      complex*16 diag(3,3),dinv(3,3),A(6,6),At(6,6),Atrace
      complex*16 c1dent(3,3),c2dent(3,3),c3dent(3,3)
      complex*16 d1dent(3,3),d2dent(3,3),d3dent(3,3),d22dent(3,3),
     +     d33dent(3,3),d23dent(3,3),d32dent(3,3)
      complex*16 Dc(3,3),Dd(3,3),Dct(3,3),Ddt(3,3)
      complex*16 ai

      real*8 amr(6,6),ami(6,6)
      real*8 zr(6,6),zi(6,6),wr(6),wi(6),Atemp(6,6)
      real*8 zr3(3,3),zi3(3,3),wr3(3),wi3(3),Dctemp(3,3),Ddtemp(3,3)
      real*8 dcr(3,3),dci(3,3),ddr(3,3),ddi(3,3)
      real*8 ec(3,nx3mx,nx2mx,ntmx),ed(3,nx3mx,nx2mx,ntmx)

      pi=4.d0*pi4
      ai=dcmplx(0.d0,1.d0)
      omn=pi/dt
      tott=dble(nt)*dt
      delom=2.d0*pi/tott
      nom=nt/2+1
      nom2=int(dble(nom-1)*dble(omwin))+1

      if(itropic.eq.0)then
         call aijkl(a6,a3)
      elseif(itropic.eq.1)then
         if(ipol.eq.0)then
            call aijkl(a6,a3)
         elseif(ipol.eq.1)then
            call aijkl(a6,a3)
            vfast=dsqrt(a3(1,1,1,1))
            vslow=dsqrt(a3(2,3,2,3))
            call isotens(a6)
            call aijkl(a6,a3)
         endif
      endif
c      call aijkl(a6,a3)

c***********************************************************************
c     This section evaluates the stability characteristics
c***********************************************************************
c     Assume homogeneous.
      call propagat(a3,c11inv,p0,p2,p3,p22,p23,p32,p33)
      do ix2=1,nx2
         do ix3=1,nx3
c     Setup diagonal reference phase matrix
            do i=1,3
               do j=1,3
                  if(i.eq.j)then
                     diag(i,j)=dcmplx(1.d0,0.d0)
                     dinv(i,j)=dcmplx(1.d0,0.d0)
                  else
                     diag(i,j)=dcmplx(0.d0,0.d0)
                     dinv(i,j)=dcmplx(0.d0,0.d0)
                  endif
               enddo
            enddo
            if(iref.eq.1)then
               p1qP=p0(1,1)
               p1qS1=p0(2,2)
               p1qS2=p0(3,3)
               p0(1,1)=0.d0
               p0(2,2)=0.d0
               p0(3,3)=0.d0
            endif
            do ic=1,3
               do jc=1,3
                  cp0(ic,jc)=dcmplx(p0(ic,jc),0.d0)
                  cp2(ic,jc)=dcmplx(p2(ic,jc),0.d0)
                  cp3(ic,jc)=dcmplx(p3(ic,jc),0.d0)
                  cp22(ic,jc)=dcmplx(p22(ic,jc),0.d0)
                  cp33(ic,jc)=dcmplx(p33(ic,jc),0.d0)
                  cp23(ic,jc)=dcmplx(p23(ic,jc),0.d0)
                  cp32(ic,jc)=dcmplx(p32(ic,jc),0.d0)
               enddo
            enddo
            
            do iom=1,nom2
               om=dble(iom-1)*delom
               rk1=om*p1inc(iwave,ix2,ix3)
               rk2=om*p2inc(iwave,ix2,ix3)
               rk3=om*p3inc(iwave,ix2,ix3)
               if(iref.eq.1)then
                  diag(1,1)=cdexp(ai*om*p1qP*dx1)
                  diag(2,2)=cdexp(ai*om*p1qS1*dx1)
                  diag(3,3)=cdexp(ai*om*p1qS2*dx1)
               endif

               call cinv3(diag,dinv)
               rhomax=-1.d10
               rmax=-1.d10

c     Evaluate amplification matrix
               Atrace=dcmplx(0.d0,0.d0) !Atrace = trace of A
               do ic=1,6
                  do jc=1,6
                     if(iom.eq.1)then
                        if(ic.eq.jc)then
                           A(ic,jc)=dcmplx(1.d0,0.d0)
                        else
                           A(ic,jc)=dcmplx(0.d0,0.d0)
                        endif
                     elseif(iom.gt.1)then
                        if(ic.le.3.and.jc.le.3)then
                           A(ic,jc)=2.d0*dx1*(ai*om*cp0(ic,jc)+
     +                          ai*cp2(ic,jc)*dsin(rk2*dx2)/dx2+
     +                          ai*cp3(ic,jc)*dsin(rk3*dx3)/dx3+
     +                          (2.d0/(ai*om))*
     +                          (cp22(ic,jc)*(dcos(rk2*dx2)-1.d0)/
     +                          (dx2*dx2)+
     +                          cp33(ic,jc)*(dcos(rk3*dx3)-1.d0)/
     +                          (dx3*dx3)-
     +                          ((cp23(ic,jc)+cp32(ic,jc))/2.d0)*
     +                          dsin(rk2*dx2)*dsin(rk3*dx3)/(dx2*dx3)))
                        elseif(ic.le.3.and.jc.gt.3)then
                           if(ic.eq.jc-3)then
                              A(ic,jc)=dinv(ic,jc-3)*dcmplx(1.d0,0.d0)
                           else
                              A(ic,jc)=dcmplx(0.d0,0.d0)
                           endif
                        elseif(ic.gt.3.and.jc.le.3)then
                           if(ic-3.eq.jc)then
                              A(ic,jc)=dcmplx(1.d0,0.d0)
                           else
                              A(ic,jc)=dcmplx(0.d0,0.d0)
                           endif
                        else
                           A(ic,jc)=dcmplx(0.d0,0.d0)
                        endif
                     endif
                     if(ic.eq.jc)Atrace=Atrace+A(ic,jc)

                     amr(ic,jc)=realpart(A(ic,jc))
                     ami(ic,jc)=imagpart(A(ic,jc))

                  enddo
               enddo

c     Evaluate Hermitian of amplification matrix
               do ic=1,6
                  do jc=1,6
                     At(ic,jc)=dconjg(A(jc,ic))
                  enddo
               enddo

c     Evaluate (A*)(A) where (A*) is transpose of conjg(A)
               do ic=1,6
                  do jc=1,6
                     Atemp(ic,jc)=0.d0
                     do kc=1,6
                        Atemp(ic,jc)=Atemp(ic,jc)+
     +                       realpart(A(ic,kc)*At(kc,jc))
                     enddo
                  enddo
               enddo

c     Evaluate eigensolution to amplification matrix using EISPACK

               call tred2(6,6,Atemp,wr,wi,zr)
               call tql2(6,6,wr,wi,zr,ierr)
               if(ierr.ne.0)then
                  write(*,*)
                  write(*,*)'Ierr.ne.0 from CG call.'
                  write(*,*)'iom,ix3,ix2',iom,ix3,ix2
c                  pause
               endif
               ievaltag=0
               do i=1,6
                  if(wr(i).gt.1.0)ievaltag=ievaltag+1
                  rhomax=max(rhomax,wr(i))
               enddo

               if(iom.eq.1)then
                  rhomax=1.d0
               endif

c     Evaluate L2 norm of amplification matrix
               do ic=1,6
                  do jc=1,6
                     Al2norm=Al2norm+Atemp(ic,jc)*Atemp(ic,jc)
                  enddo
               enddo
               Al2norm=dsqrt(Al2norm)

               rmax=max(rhomax,rmax)

            enddo

         enddo
      enddo

c***********************************************************************
c     This section evaluates the dispersion characteristics
c***********************************************************************
c     Central grid point
      kx2=nx2/2+1
      kx3=nx3/2+1
      dx=dx2                    !Assuming square lateral grid
      call propagat(a3,c11inv,p0,p2,p3,p22,p23,p32,p33)

      do ix2=1,nx2
         do ix3=1,nx3
c     Setup diagonal reference phase matrix
            if(iref.eq.1)then
               p0(1,1)=0.d0
               p0(2,2)=0.d0
               p0(3,3)=0.d0
            endif
            do ic=1,3
               do jc=1,3
                  cp0(ic,jc)=dcmplx(p0(ic,jc),0.d0)
                  cp2(ic,jc)=dcmplx(p2(ic,jc),0.d0)
                  cp3(ic,jc)=dcmplx(p3(ic,jc),0.d0)
                  cp22(ic,jc)=dcmplx(p22(ic,jc),0.d0)
                  cp33(ic,jc)=dcmplx(p33(ic,jc),0.d0)
                  cp23(ic,jc)=dcmplx(p23(ic,jc),0.d0)
                  cp32(ic,jc)=dcmplx(p32(ic,jc),0.d0)
               enddo
            enddo
            do iom=1,nom2
               om=dble(iom-1)*delom
c     Setup diagonal wavenumber arrays
               do i=1,3
                  c1dent(i,i)=dcmplx(om*p1inc(i,ix2,ix3),0.d0)
                  c2dent(i,i)=dcmplx(om*p2inc(i,ix2,ix3),0.d0)
                  c3dent(i,i)=dcmplx(om*p3inc(i,ix2,ix3),0.d0)
                  d1dent(i,i)=dcmplx(dsin(om*p1inc(i,ix2,ix3)*dx1)/
     +                 dx1,0.d0)
                  d2dent(i,i)=dcmplx(dsin(om*p2inc(i,ix2,ix3)*dx)/
     +                 dx,0.d0)
                  d3dent(i,i)=dcmplx(dsin(om*p3inc(i,ix2,ix3)*dx)/dx,
     +                 0.d0)
                  d22dent(i,i)=dcmplx((dcos(om*p2inc(i,ix2,ix3)*dx)-
     +                 1.d0)/(dx*dx),0.d0)
                  d33dent(i,i)=dcmplx((dcos(om*p3inc(i,ix2,ix3)*dx)-
     +                 1.d0)/(dx*dx),0.d0)
                  d23dent(i,i)=dcmplx((dsin(om*p2inc(i,ix2,ix3)*dx)*
     +                 dsin(om*p3inc(i,ix2,ix3)*dx)/(dx*dx)),0.d0)
                  d32dent(i,i)=d23dent(i,i)
                  do j=1,3
                     if(i.ne.j)then
                        c1dent(i,j)=dcmplx(0.d0,0.d0)
                        c2dent(i,j)=dcmplx(0.d0,0.d0)
                        c3dent(i,j)=dcmplx(0.d0,0.d0)
                        d1dent(i,j)=dcmplx(0.d0,0.d0)
                        d2dent(i,j)=dcmplx(0.d0,0.d0)
                        d3dent(i,j)=dcmplx(0.d0,0.d0)
                        d22dent(i,j)=dcmplx(0.d0,0.d0)
                        d33dent(i,j)=dcmplx(0.d0,0.d0)
                        d23dent(i,j)=dcmplx(0.d0,0.d0)
                        d32dent(i,j)=dcmplx(0.d0,0.d0)
                     endif
                  enddo

                  do icomp=1,3
                     do ic=1,3
                        do jc=1,3
                           ik=icomp
                           jk=icomp
c     Continuous wave-equation dispersion matrix
                           Dc(ic,jc)=om*om*cp0(ic,jc)+om*(
     +                          c2dent(ik,jk)*cp2(ic,jc)+
     +                          c3dent(ik,jk)*cp3(ic,jc)-
     +                          c1dent(ik,jk))+
     +                          (c2dent(ik,jk)*c2dent(ik,jk)*
     +                          cp22(ic,jc)+
     +                          c3dent(ik,jk)*c3dent(ik,jk)*
     +                          cp33(ic,jc)+
     +                          c2dent(ik,jk)*c3dent(ik,jk)*
     +                          (cp23(ic,jc)+cp32(ic,jc)))
c     Discrete wave-equation dispersion matrix
                           Dd(ic,jc)=om*om*cp0(ic,jc)+om*(
     +                          d2dent(ik,jk)*cp2(ic,jc)+
     +                          d3dent(ik,jk)*cp3(ic,jc)-
     +                          d1dent(ik,jk))+(d23dent(ik,jk)*
     +                          (cp23(ic,jc)+cp32(ic,jc))-
     +                          2.d0*d22dent(ik,jk)*cp22(ic,jc)-
     +                          2.d0*d33dent(ik,jk)*cp33(ic,jc))
                           dcr(ic,jc)=realpart(Dc(ic,jc))
                           dci(ic,jc)=imagpart(Dc(ic,jc))
                           ddr(ic,jc)=realpart(Dd(ic,jc))
                           ddi(ic,jc)=imagpart(Dd(ic,jc))
                        enddo
                     enddo
c     Evaluate Hermitian of dispersion matrix
                     if(iref.eq.0)then
                        if(iom.eq.1)then
                           Dctrace=1.d0
                           Ddtrace=1.d0
                        else
                           Dctrace=realpart(Dc(1,1)+Dc(2,2)+Dc(3,3))
                           Ddtrace=realpart(Dd(1,1)+Dd(2,2)+Dd(3,3))
                        endif
                        do ic=1,3
                           do jc=1,3
                              Dc(ic,jc)=Dc(ic,jc)/Dctrace
                              Dd(ic,jc)=Dd(ic,jc)/Ddtrace
                           enddo
                        enddo
                     endif
                     do ic=1,3
                        do jc=1,3
                           Dct(ic,jc)=dconjg(Dc(jc,ic))
                           Ddt(ic,jc)=dconjg(Dd(jc,ic))
                        enddo
                     enddo
c     Evaluate (A*)(A) where (A*) is transpose of conjg(A)
                     do ic=1,3
                        do jc=1,3
                           Dctemp(ic,jc)=0.d0
                           Ddtemp(ic,jc)=0.d0
                           do kc=1,3
                              Dctemp(ic,jc)=Dctemp(ic,jc)+
     +                             realpart(Dc(ic,kc)*Dct(kc,jc))
                              Ddtemp(ic,jc)=Ddtemp(ic,jc)+
     +                             realpart(Dd(ic,kc)*Ddt(kc,jc))
                           enddo
                        enddo
                     enddo
c     Evaluate eigensolution to dispersion matrix using EISPACK
c     Continuous wave-equation
                     call tred2(3,3,Dctemp,wr3,wi3,zr3)
                     call tql2(3,3,wr3,wi3,zr3,ierr)
                     if(ierr.ne.0)then
                        write(*,*)
                        write(*,*)'Ierr.ne.0 from CG call.'
                        write(*,*)'iom,ix3,ix2',iom,ix3,ix2
c                        pause
                     endif
                     if(iref.eq.0)then
                        ec(icomp,ix3,ix2,iom)=dsqrt(wr3(1))
                     elseif(iref.eq.1)then
                        ec(icomp,ix3,ix2,iom)=dsqrt(wr3(3))
                     endif
c     Discrete wave-equation
                     call tred2(3,3,Ddtemp,wr3,wi3,zr3)
                     call tql2(3,3,wr3,wi3,zr3,ierr)
                     if(ierr.ne.0)then
                        write(*,*)
                        write(*,*)'Ierr.ne.0 from CG call.'
                        write(*,*)'iom,ix3,ix2',iom,ix3,ix2
c                        pause
                     endif
                     if(iref.eq.0)then
                        ed(icomp,ix3,ix2,iom)=dsqrt(wr3(1))
                     elseif(iref.eq.1)then
                        ed(icomp,ix3,ix2,iom)=dsqrt(wr3(3))
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
c     Loop to determine maximum and minimum eigenvalues
      do i=1,3
         rkval(i)=-1.d10
         if(iref.eq.0)ik=i
         if(iref.eq.1)ik=i+3
         do ix2=1,nx2
            do ix3=1,nx3
               do iom=1,nom2
                  om=dble(iom-1)*delom
                  if(iref.eq.0)then
                     rkctemp=dsqrt(
     +                    ec(i,ix3,ix2,iom)**2.d0+
     +                    (om*p2inc(i,ix2,ix3))**2.d0+
     +                    (om*p3inc(i,ix2,ix3))**2.d0)
                     rkdtemp=dsqrt(
     +                    ed(i,ix3,ix2,iom)**2.d0+
     +                    (om*p2inc(i,ix2,ix3))**2.d0+
     +                    (om*p3inc(i,ix2,ix3))**2.d0)
                     dkc=dble(om/omn)
                  elseif(iref.eq.1)then
                     rkctemp=dsqrt(
     +                    ec(i,ix3,ix2,iom)**2.d0+
     +                    (om*p2inc(i,ix2,ix3))**2.d0+
     +                    (om*p3inc(i,ix2,ix3))**2.d0)
                     rkdtemp=dsqrt(
     +                    ed(i,ix3,ix2,iom)**2.d0+
     +                    (om*p2inc(i,ix2,ix3))**2.d0+
     +                    (om*p3inc(i,ix2,ix3))**2.d0)
                     dkc=dble(om/omn)
                  endif
                  rtemp=dabs(rkdtemp/rkctemp-1.d0)
                  rkval(i)=max(rkval(i),rtemp)
c                  write(*,*)i,ec(i,ix3,ix2,iom),p2inc(i,ix2,ix3),
c     +                 p3inc(i,ix2,ix3)
c                  read(*,*)
               enddo
            enddo
         enddo
      enddo

c     Write stability and dispersion results to screen
c     Stability results in terms of maximum eigenvalue or spectral 
c     radius (stable if rmax<1.0005)
      write(*,*)'             (i) Stability'
      if(rmax.gt.5.d0)then
      write(*,*)'                 Warning: spectral amplitude for FD'
      write(*,*)'                 scheme greater than 5.0 for at least'
      write(*,*)'                 of the frequency components:'
      write(*,*)'                 max(spectral amplitude)=',rmax
c      pause
      endif
      if(rmax.gt.1.0005d0)then
      write(*,*)'                 Warning: spectral amplitude for FD'
      write(*,*)'                 scheme greater than 1.0005 for at'
      write(*,*)'                 at least one of the frequency'
      write(*,*)'                 components.'
      write(*,*)'                 max(spectral amplitude)=',rmax
c      pause
      endif
      if(rmax.le.1.0005d0)then
      write(*,*)'                 Propagator stable for the FD scheme'
      write(*,*)'                 for all frequencies (i.e., maximum'
      write(*,*)'                 eigenvalue less than or equal to'
      write(*,*)'                 1.0005).'
      write(*,*)'                 -> maximum eigenvalue:',rmax
      endif

c     Dispersion results:
      write(*,*)'             (i) Dispersion'
      do icomp=1,3
      rktmp=rkval(icomp)
      if(icomp.eq.1.and.iwave.eq.icomp)then
      if(rktmp.le.0.1d0)then
      write(*,*)'                 Pwave: optimal dispersion range,'
      elseif(rktmp.gt.0.1d0.and.rktmp.le.0.2d0)then
      write(*,*)'                 Pwave might be dispersive,'
      else
      write(*,*)'                 Pwave in dispersive range'
      write(*,*)'                 Continue with extrapolation?'
      write(*,*)
c      pause
      endif
      elseif(icomp.ge.2.and.iwave.eq.icomp)then
      if(rktmp.le.0.1d0)then
      write(*,*)'                 Swave: optimal dispersion range,'
      elseif(rktmp.gt.0.1d0.and.rktmp.le.0.2d0)then
      write(*,*)'                 Swave might be dispersive'
      else
      write(*,*)'                 Swave in dispersive range'
      write(*,*)'                 Continue with extrapolation?'
      write(*,*)
c      pause
      endif
      endif
      if(icomp.eq.iwave)then
      write(*,*)'                 -> absolute value of difference'
      write(*,*)'                 between one and dispersion ratio,'
      write(*,*)'                 phase velocity (cont/disc):',rktmp
      write(*,*)
      endif
      enddo

      return
      end

c***********************************************************************
      subroutine sortslwc(wr,wi,zr,zi,wc,zc)
c     Sorts eigenvalues/vectors into separate up/down sets.
c     This version sorts real and imaginary parts into positive and
c     negative groups. Does not sort within a group (unecessary?).
c     C. Thomson, February 1995.
c***********************************************************************

      implicit none

      integer ipos,ineg,i,l
      real*8 acc,test
      real*8 wr(6),wi(6),zr(6,6),zi(6,6)
      real*8 wwrkr(6),wwrki(6),zwrkr(6,6),zwrki(6,6)
      complex*16 wc(6),zc(6,6)

c     Sort according to size of p3 (vertical slowness)

      ipos=0
      ineg=0
      acc=1.d-12
      do 100 l=1,6
         test=0.d0
         if(wr(l).eq.0.d0)then
            if(wi(l).eq.0.d0)then
               write(*,*)'Error 1 in Sortslwc'
               stop
            endif
            test=2.d0*acc
         else
            if(wi(l).ne.0.d0)then
               test=dabs(wi(l)/wr(l))
            endif
         endif
         if(test.gt.acc)then
            if(wi(l).gt.0.d0)then
               ipos=ipos+1
               wwrkr(ipos)=wr(l)
               wwrki(ipos)=wi(l)
               call trnsfv6(zr(1,l),zi(1,l),
     .              zwrkr(1,ipos),zwrki(1,ipos))
            else
               ineg=ineg+1
               wwrkr(ineg+3)=wr(l)
               wwrki(ineg+3)=wi(l)
               call trnsfv6(zr(1,l),zi(1,l),
     .              zwrkr(1,ineg+3),zwrki(1,ineg+3))
            endif
         else
            if(wr(l).gt.0.d0)then
               ipos=ipos+1
               wwrkr(ipos)=wr(l)
               wwrki(ipos)=wi(l)
               call trnsfv6(zr(1,l),zi(1,l),
     .              zwrkr(1,ipos),zwrki(1,ipos))
            else
               ineg=ineg+1
               wwrkr(ineg+3)=wr(l)
               wwrki(ineg+3)=wi(l)
               call trnsfv6(zr(1,l),zi(1,l),
     .              zwrkr(1,ineg+3),zwrki(1,ineg+3))
            endif
         endif
 100  continue

      if(ipos.ne.3.or.ineg.ne.3)write(*,*)'Error 2 in Sortslwc'

      do 101 i=1,6
         wc(i)=dcmplx(wwrkr(i),wwrki(i))
         do 102 l=1,6
            zc(l,i)=dcmplx(zwrkr(l,i),zwrki(l,i))
 102     continue
 101  continue

      return
      end

c***********************************************************************
      subroutine trnsfv6(vin1,vin2,vout1,vout2)
c***********************************************************************

      implicit none

      integer l
      real*8 vin1(6),vin2(6),vout1(6),vout2(6)

      do 100 l=1,6
         vout1(l)=vin1(l)
         vout2(l)=vin2(l)
 100  continue

      return
      end
