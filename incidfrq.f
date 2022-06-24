c...+....|....+....|....+....|....+....|....+....|....+....|....+....|..
c  File Incidfrq.f contains the following 5 subroutines: Setincid, Vgrid
c  Scsinc, Rickerw and Dot.
c
c  Last modified: April 16, 2012.
c***********************************************************************
      subroutine setincid(tinc,nx1st,nx2,nx3,dx1,dx2,dx3,a6,uvec,wvf0,
     +     nt,dt,p1inc,p2inc,p3inc,mwid,nord,rpeak,sigvec1,beta)
c     Set up the initial conditions for the one-way operator.  
c     This version sets up a smoothly curved (convex) wavefront. The 
c     curved wavefront is described in terms of incidence angles at the
c     centre of the grid and deviations from those angles at four 
c     cross points.
c     Intended for use in the conical-point problem.
c     Displacements are set using true eigenvectors at the resulting
c     grid slownesses.
c     Waveform: One version sets the time-domain waveform as a boxcar 
c               smoothed (connvolved) several times with itself 
c               -- i.e. a sinc function in the frequency domain.
c               Parameters mwid and nord define the boxcar.
c               An alternative choice uses a cosine-like time pulse.
c***********************************************************************

      implicit none
      include '../Input/narc_dp.par'

      integer nx1st,nx2,nx3,nt,nom,iom,mwid,nord,nreal
      integer icomp,iform,ishft,it,ix1,ix2c,ix3c,ix2,ix3,i,j
      integer iflag(6),k2vec(2)

      real*8 dx,dx1,dx2,dx3,x2c,x3c,x2,x3,rnx,half,h2
      real*8 dt,tc,omn,om,delom,ttot,time,tshift
      real*8 pi,pi2,rpeak,beta,vnor,cphi,sphi
      real*8 walpha,palpha,temp1,temp2,temp3,temp4,temp7,temp8
      real*8 theta,dtheta,phi,dphi,dpx,dpy,pxwid,pywid
      real*8 thran1,thran2,phran1,phran2
      real*8 p1c,p2c,p3c,p2,p3,p22,p23,p32,p33
      real*8 phimx,phimn,thetamx,thetamn,thetmx,thetmn,thet2
      real*8 del1,del2,del3,del4
      real*8 p2mxp,p2mnp,p3mxp,p3mnp,p2mxt,p2mnt,p3mxt,p3mnt
      real*8 sign2,sign3,sdot,atmp,btmp,ctmp,disp1,disp2,disp3

      real*8 tinc(nx2mx,nx3mx),p1inc(3,nx2mx,nx3mx),
     +     p2inc(3,nx2mx,nx3mx),p3inc(3,nx2mx,nx3mx)
      real*8 a6(6,6,nx3mx,nx2mx,nx1sto),wvf0(ntmx)
      real*8 a(6,6),a3r(3,3,3,3),a3iso(3,3,3,3)
      real*8 vn(3),p(3,3),vg(3,3),px1(3),wn(3)
      real*8 p2ran(2),p3ran(2),eval(6),evec(6,6)
      real*8 sigvec(3),sigvec1(3,nx3mx)

      complex*16 ai,ctemp
      complex*16 uvec(3,ntmx,nx3mx,nx2mx,3)
      complex*16 px1c(6)
      complex*8 cfilt(ntmx) 

      real*4 singpi,omegai,sdelom,sdt,srpeak

      if(nx2.gt.nxmx.or.nx3.gt.nxmx)then
         write(11,*)
         write(11,*)'Error: nx2,nx3 = ',nx2,nx3,' but nxmx =',nxmx
         write(11,*)'Stopping program in setincid.'
         stop
      endif

      pi=4.d0*pi4
      pi2=pi/2.d0
      singpi=sngl(pi)

c     Time and elasticity at grid center
      tc=0.d0
      tmin=1.d5
      tmax=-1.d5
      ix2c=(nx2/2)+1
      ix3c=(nx3/2)+1
      ix1=1
      x2c=dble(ix2c-1)*dx2
      x3c=dble(ix3c-1)*dx3

      do i=1,6
         do j=1,6
            a(i,j)=a6(i,j,ix3c,ix2c,ix1)
         enddo
      enddo

      if(itropic.eq.0)then
         call aijkl(a,a3r)
      elseif(itropic.eq.1)then
         if(ipol.eq.0)then
            call aijkl(a,a3r)
         elseif(ipol.eq.1)then
            call aijkl(a,a3r)
            vfast=dsqrt(a3r(1,1,1,1))
            vslow=dsqrt(a3r(2,3,2,3))
            write(*,*)'             -> Note: option choosen to manually'
            write(*,*)'                specify initial wavefront'
            write(*,*)'                body-wave polarizations.'
            write(*,*)
            call isotens(a)
            call aijkl(a,a3r)
         endif
      endif

c     Section to determine max allowable slownesses for narrow-angle
c     propagation
      if(iwave.eq.1)then
         walpha=dsin(15.d0*pi/180.d0)
         palpha=walpha/vfast
      else
         walpha=dsin(15.d0*pi/180.d0)
         palpha=walpha/vslow
      endif

c     Get incident wave type here (qP=1, qS1=2 and qS2=3)
      call incwft(iform,temp1,temp2,temp3,temp4)

c     Set variable grid parameters
      half=dble(nx2-1)*dx2/2.d0
      h2=2.d0*half
      call vgrid(h2,nx2,beta)

      do icomp=1,3

c     Find slownesses at grid centre as well as slowness gradients
         if(iform.eq.0)then     !Incident wvfrnt defined by angles
            theta=temp1*pi/180.d0
            phi=temp2*pi/180.d0
            dtheta=temp3*pi/180.d0
            dphi=temp4*pi/180.d0
            rnx=dble(nx2)      !Assuming square grid
            dx=dx2
            thran1=-dtheta*(rnx/2.d0)*dx
            thran2=dtheta*(rnx/2.d0)*dx
            phran1=-dphi*(rnx/2.d0)*dx
            phran2=dphi*(rnx/2.d0)*dx
            if(itropic.eq.0)then
               wn(1)=dcos(theta)
               wn(2)=dsin(theta)*dcos(phi)
               wn(3)=dsin(theta)*dsin(phi)
               if(icomp.eq.1)then
                  vnor=vfast
               else
                  vnor=vslow
               endif
               p(icomp,1)=wn(1)/vnor
               p(icomp,2)=wn(2)/vnor
               p(icomp,3)=wn(3)/vnor
            elseif(itropic.eq.1)then
               if(ipol.eq.0)then
                  call slowness(a3r,phi,theta,vn,p,vg)
               elseif(ipol.eq.1)then
                  wn(1)=dcos(theta)
                  wn(2)=dsin(theta)*dcos(phi)
                  wn(3)=dsin(theta)*dsin(phi)
                  if(icomp.eq.1)then
                     vnor=vfast
                  else
                     vnor=vslow
                  endif
                  p(icomp,1)=wn(1)/vnor
                  p(icomp,2)=wn(2)/vnor
                  p(icomp,3)=wn(3)/vnor
               endif
            endif
            p1c=p(icomp,1)
            p2c=p(icomp,2)
            p3c=p(icomp,3)
         elseif(iform.eq.1)then !Incident wvfrnt defined by slownesses
            if(itropic.eq.0)then
               if(icomp.eq.1)then
                  p(icomp,1)=dsqrt(1.d0/(vfast**2.d0)-temp1**2.d0-
     +                 temp2**2.d0)
               else
                  p(icomp,1)=dsqrt(1.d0/(vslow**2.d0)-temp1**2.d0-
     +                 temp2**2.d0)
               endif
            elseif(itropic.eq.1)then
               if(ipol.eq.0)then
                  call vertslw(a3r,temp1,temp2,px1c,iflag,nreal)
                  call sortp(px1c,px1,iflag,nreal)
                  if(icomp.eq.1)then
                     if(nreal.eq.3)then
                        p(icomp,1)=dble(px1(1))
                     else
                        write(11,*)'Input problems: iwave=',iwave
                        write(11,*)'but nreal=',nreal
                        write(11,*)'(Pwave evanescent)'
                        write(11,*)'Program stopped within setincid.'
                        stop
                     endif
                  elseif(icomp.eq.2)then
                     if(nreal.eq.3)then
                        p(icomp,1)=dble(px1(2))
                     elseif(nreal.eq.2)then
                        p(icomp,1)=dble(px1(1))
                     else
                        write(11,*)'Input problems: iwave=',icomp
                        write(11,*)'but nreal=',nreal,' (p-, s1-wave 
     +                       evanescent)'
                        write(11,*)'Program stopped within setincid.'
                        stop
                     endif
                  elseif(icomp.eq.3)then
                     if(nreal.eq.3)then
                        p(icomp,1)=dble(px1(3))
                     elseif(nreal.eq.2)then
                        p(icomp,1)=dble(px1(2))
                    elseif(nreal.eq.1)then
                        p(icomp,1)=dble(px1(1))
                     else
                        write(11,*)'Wave type :',icomp,' evanescent.'
                        write(11,*)'Only ',nreal,' non-evanescent wave 
     +                       type(s).'
                        write(11,*)'Program stopped.'
                        stop
                     endif
                  endif
               elseif(ipol.eq.1)then
                  if(icomp.eq.1)then
                     p(icomp,1)=dsqrt(1.d0/(vfast**2.d0)-temp1**2.d0-
     +                    temp2**2.d0)
                  else
                     p(icomp,1)=dsqrt(1.d0/(vslow**2.d0)-temp1**2.d0-
     +                    temp2**2.d0)
                  endif
               endif
            endif

c     Define incident centre slowness vector and range.
            p(icomp,1)=dble(p(icomp,1))
            p(icomp,2)=temp1
            p(icomp,3)=temp2
            dpx=temp3
            dpy=temp4
            rnx=dble(nx2)      !Assuming square grid
            dx=dx2
            pxwid=dpx*(rnx/2.d0)
            pywid=dpy*(rnx/2.d0)
            p1c=p(icomp,1)
            p2c=p(icomp,2)
            p3c=p(icomp,3)
c     Find theta and phi corresponding to this slowness
            phi=datan2(p3c,p2c)
            theta=datan2(p2c/dcos(phi),p1c)

         endif

         if((itropic.eq.1).and.(nreal.ne.3))then
            write(11,*)'Warning: nreal of centre grid < 3 in setinc.'
            write(11,*)'nreal=',nreal
         endif

         cphi=dcos(phi)
         sphi=dsin(phi)

         if(iform.eq.0)then
            if(icomp.eq.1)then
               vnor=vfast
            else
               vnor=vslow
            endif
c     First examine the phi limits at constant theta.
            phimx=phi+datan2(dtan(phran2),dsin(theta))   
            thet2=dacos(dcos(theta)*dcos(phran2))
            if(thet2.gt.pi2) write(11,*)'Error: thet2 in setincid'
            if(itropic.eq.0)then
               wn(1)=dcos(thet2)
               wn(2)=dsin(thet2)*dcos(phimx)
               wn(3)=dsin(thet2)*dsin(phimx)
               p(icomp,1)=wn(1)/vnor
               p(icomp,2)=wn(2)/vnor
               p(icomp,3)=wn(3)/vnor
            elseif(itropic.eq.1)then
               if(ipol.eq.0)then
                  call slowness(a3r,phimx,thet2,vn,p,vg)
               elseif(ipol.eq.1)then
                  wn(1)=dcos(thet2)
                  wn(2)=dsin(thet2)*dcos(phimx)
                  wn(3)=dsin(thet2)*dsin(phimx)
                  p(icomp,1)=wn(1)/vnor
                  p(icomp,2)=wn(2)/vnor
                  p(icomp,3)=wn(3)/vnor
               endif
            endif
            p2mxp=p(icomp,2)
            p3mxp=p(icomp,3)

            phimn=phi-datan2(dtan(dabs(phran1)),dsin(theta)) 
            thet2=dacos(dcos(theta)*dcos(phran1))
            if(thet2.gt.pi2) write(11,*)'Error: thet2 in setincid'
            if(itropic.eq.0)then
               wn(1)=dcos(thet2)
               wn(2)=dsin(thet2)*dcos(phimn)
               wn(3)=dsin(thet2)*dsin(phimn)
               p(icomp,1)=wn(1)/vnor
               p(icomp,2)=wn(2)/vnor
               p(icomp,3)=wn(3)/vnor
            elseif(itropic.eq.1)then
               if(ipol.eq.0)then
                  call slowness(a3r,phimn,thet2,vn,p,vg)
               elseif(ipol.eq.1)then
                  wn(1)=dcos(thet2)
                  wn(2)=dsin(thet2)*dcos(phimn)
                  wn(3)=dsin(thet2)*dsin(phimn)
                  p(icomp,1)=wn(1)/vnor
                  p(icomp,2)=wn(2)/vnor
                  p(icomp,3)=wn(3)/vnor
               endif
            endif
            p2mnp=p(icomp,2)
            p3mnp=p(icomp,3)
         
c     Next look at theta limits at constant phi (almost!) ...
            thetmx=theta+thran2
            if(thetmx.gt.pi2)write(11,*)'***thetmx.gt.pi/2 in setincid'
            if(thetmx.gt.pi)then
               write(11,*)'***Error: thetmx.gt.pi in setincid'
               stop
            endif
            if(itropic.eq.0)then
               wn(1)=dcos(thetmx)
               wn(2)=dsin(thetmx)*dcos(phi)
               wn(3)=dsin(thetmx)*dsin(phi)
               p(icomp,1)=wn(1)/vnor
               p(icomp,2)=wn(2)/vnor
               p(icomp,3)=wn(3)/vnor
            elseif(itropic.eq.1)then
               if(ipol.eq.0)then
                  call slowness(a3r,phi,thetmx,vn,p,vg)
               elseif(ipol.eq.1)then
                  wn(1)=dcos(thetmx)
                  wn(2)=dsin(thetmx)*dcos(phi)
                  wn(3)=dsin(thetmx)*dsin(phi)
                  p(icomp,1)=wn(1)/vnor
                  p(icomp,2)=wn(2)/vnor
                  p(icomp,3)=wn(3)/vnor
               endif
            endif
            p2mxt=p(icomp,2)
            p3mxt=p(icomp,3)

            thetmn=theta+thran1
            if(thetmx.lt.0.d0)then
               thetmn=-thetmn
               phi=phi+pi
            elseif(thetmn.lt.-pi2)then
               write(11,*)'***Error: thetmx.lt.-pi in setincid'
               stop
            endif
            if(itropic.eq.0)then
               wn(1)=dcos(thetmn)
               wn(2)=dsin(thetmn)*dcos(phi)
               wn(3)=dsin(thetmn)*dsin(phi)
               p(icomp,1)=wn(1)/vnor
               p(icomp,2)=wn(2)/vnor
               p(icomp,3)=wn(3)/vnor
            elseif(itropic.eq.1)then
               if(ipol.eq.0)then
                  call slowness(a3r,phi,thetmn,vn,p,vg)
               elseif(ipol.eq.1)then
                  wn(1)=dcos(thetmn)
                  wn(2)=dsin(thetmn)*dcos(phi)
                  wn(3)=dsin(thetmn)*dsin(phi)
                  p(icomp,1)=wn(1)/vnor
                  p(icomp,2)=wn(2)/vnor
                  p(icomp,3)=wn(3)/vnor
               endif
            endif
            p2mnt=p(icomp,2)
            p3mnt=p(icomp,3)

c     Check whether max/min slownesses greater than allowable slownesses
            if(dabs(p2mxp).gt.palpha.or.dabs(p3mxp).gt.palpha.or.
     +           dabs(p2mnp).gt.palpha.or.dabs(p3mnp).gt.palpha.or.
     +           dabs(p2mxt).gt.palpha.or.dabs(p3mxt).gt.palpha.or.
     +           dabs(p2mnt).gt.palpha.or.dabs(p3mnt).gt.palpha)then
               write(*,*)'             Problem with initial conditions:'
               write(*,*)'             -> maximum lateral slowness is:'
               write(*,*)'               ',palpha
               write(*,*)'             -> from input file:'
               write(*,*)'                 p2mxp,p3mxp:',p2mxp,p3mxp
               write(*,*)'                 p2mnp,p3mnp:',p2mnp,p3mnp
               write(*,*)'                 p2mxt,p3mxt:',p2mxt,p3mxt
               write(*,*)'                 p2mnt,p3mnt:',p2mnt,p3mnt
               write(*,*)'             Values exceed narrow-angle'
               write(*,*)'             approximation (+/- 15 degrees).'
               write(*,*)
               stop
            endif

c     Infer transverse slowness gradients
            del1=p2mxt-p2mnt
            del2=p3mxt-p3mnt
            del3=p2mxp-p2mnp
            del4=p3mxp-p3mnp
            p22=(del1*cphi-del3*sphi)/h2
            p23=(del1*sphi+del3*cphi)/h2
            p32=(del2*cphi-del4*sphi)/h2
            p33=(del2*sphi+del4*cphi)/h2
            p23=(p23+p32)/2.d0
            p32=p23

         elseif(iform.eq.1)then
c     Wavefront evaluated using slownesses:
c     Transverse slowness gradients
            p22=dpx
            p33=dpy
            p23=0.d0    !Seems to work well with p23=p32=0
            p32=0.d0    !This is a small loss in the quadratic approx.
         endif

c     Horizontal slowness lateral derivatives now known.
c     Looping over grid to evaluate incident times, slownesses, 
c     polarizations
         p2ran(1)=0.d0
         p3ran(1)=0.d0
         p2ran(2)=0.d0
         p3ran(2)=0.d0

c     Now set grid values of time, etc.

         do 10 ix2=1,nx2
            x2=xwor(ix2)
            do 20 ix3=1,nx3
               x3=xwor(ix3)
               tinc(ix2,ix3)=tc+p2c*(x2-x2c)+p3c*(x3-x3c)+
     +              0.5d0*p22*(x2-x2c)**2.d0+
     +              0.5d0*p33*(x3-x3c)**2.d0+
     +              p23*(x2-x2c)*(x3-x3c)
               tmin=min(tinc(ix2,ix3),tmin)
               tmax=max(tinc(ix2,ix3),tmax)

c     In plane slowness at each point, and its range
               p2=p2c+p22*(x2-x2c)+p23*(x3-x3c)
               p3=p3c+p33*(x3-x3c)+p32*(x2-x2c)

               if(itropic.eq.0)then
                  if(iwave.eq.1)then
                     p(icomp,1)=
     +                    dsqrt(1.d0/(vfast**2.d0)-p2**2.d0-p3**2.d0)
                  else
                     p(icomp,1)=
     +                    dsqrt(1.d0/(vslow**2.d0)-p2**2.d0-p3**2.d0)
                  endif
               elseif(itropic.eq.1)then
                  if(ipol.eq.0)then
                     call vertslw(a3r,p2,p3,px1c,iflag,nreal)
                     call sortp(px1c,px1,iflag,nreal)
                     if(icomp.eq.1)then
                        if(nreal.eq.3)then
                           p(icomp,1)=dble(px1(1))
                        else
                           write(11,*)'Input problems: iwave=',iwave
                           write(11,*)'but nreal=',nreal
                           write(11,*)'(Pwave evanescent)'
                           write(11,*)'Program stopped within setincid.'
                           stop
                        endif
                     elseif(icomp.eq.2)then
                        if(nreal.eq.3)then
                           p(icomp,1)=dble(px1(2))
                        elseif(nreal.eq.2)then
                           p(icomp,1)=dble(px1(1))
                        else
                           write(11,*)'Input problems: iwave=',icomp
                           write(11,*)'but nreal=',nreal,' (p-, s1-wave 
     +                          evanescent)'
                           write(11,*)'Program stopped within setincid.'
                           stop
                        endif
                     elseif(icomp.eq.3)then
                        if(nreal.eq.3)then
                           p(icomp,1)=dble(px1(3))
                        elseif(nreal.eq.2)then
                           p(icomp,1)=dble(px1(2))
                        elseif(nreal.eq.1)then
                           p(icomp,1)=dble(px1(1))
                        else
                           write(11,*)'Wave type :',icomp,' evanescent.'
                           write(11,*)'Only ',nreal,' non-evanescent
     +                          wave type(s).'
                           write(11,*)'Program stopped.'
                           stop
                        endif
                     endif
                  elseif(ipol.eq.1)then
                     if(iwave.eq.1)then
                        p(icomp,1)=
     +                       dsqrt(1.d0/(vfast**2.d0)-p2**2.d0-p3**2.d0)
                     else
                        p(icomp,1)=
     +                       dsqrt(1.d0/(vslow**2.d0)-p2**2.d0-p3**2.d0)
                     endif
                  endif
               endif

               p1inc(icomp,ix2,ix3)=p(icomp,1)
               p2inc(icomp,ix2,ix3)=p2
               p3inc(icomp,ix2,ix3)=p3
               p2ran(1)=min(p2ran(1),p2inc(icomp,ix2,ix3))
               p3ran(1)=min(p3ran(1),p3inc(icomp,ix2,ix3))
               p2ran(2)=max(p2ran(2),p2inc(icomp,ix2,ix3))
               p3ran(2)=max(p3ran(2),p3inc(icomp,ix2,ix3))
               temp7=p2inc(icomp,ix2,ix3)
               temp8=p3inc(icomp,ix2,ix3)

c     True normal slowness, polarization and stress for this 
c     transverse slowness
               if(itropic.eq.0)then
                  call isoev(a3r,p2,p3,eval,evec,nreal)
               elseif(itropic.eq.1)then
                  if(ipol.eq.0)then
                     call anisoev(a3r,p2,p3,eval,evec,iflag,nreal,k2vec)
                  elseif(ipol.eq.1)then
                     call isoev(a3r,p2,p3,eval,evec,nreal)
                  endif
               endif

               if(icomp.eq.iwave)then
                  atmp=0.d0
                  btmp=0.d0
                  ctmp=0.d0

                  if(iwave.eq.1)then
                     if(nreal.eq.3)then
                        atmp=evec(1,iwave)
                        if(atmp.ne.0.d0)then
c     Displacement set for a P-wave with a positive 1 component (i.e.
c     in the direction of propagation)
                           atmp=atmp/dabs(atmp)
                           disp1=atmp*evec(1,iwave)
                           disp2=atmp*evec(2,iwave)
                           disp3=atmp*evec(3,iwave)
                        else
                           write(11,*)'Setincid error 1: nreal,atmp=',
     +                          nreal,atmp,' (stopping)'
                           write(*,*)'ix2,ix3,nreal',ix2,ix3,nreal
c                           pause
                           stop
                        endif
                     else
                        write(11,*)'Input problem. iwave=1 but nreal'
                        write(11,*)'lt 3 at (ix2,ix3)=',ix2,ix3
                        write(11,*)' (stopping)'
                        write(*,*)'ix2,ix3,nreal',ix2,ix3,nreal
c                        pause
                        stop
                     endif

                  elseif((iwave.eq.2.or.iwave.eq.3).and.nreal.ne.1)then
c     Setting for an S-wave needs a more complicated scheme to impose
c     smoothness (except for a possible line discontinuity due to a
c     conical point)
                     ishft=0
                     if(nreal.eq.2)ishft=1
                     if(ipol.eq.0)then
                        disp1=0.0
                        disp2=1.0
                        disp3=0.0
                     elseif(ipol.eq.1)then
                        disp1=0.d0
                        disp2=dcos(s1phi)
                        disp3=dsin(s1phi)
                     endif
                     btmp=disp1*evec(1,2-ishft)+disp2*evec(2,2-ishft)
     +                    +disp3*evec(3,2-ishft)
                     ctmp=disp1*evec(1,3-ishft)+disp2*evec(2,3-ishft)
     +                    +disp3*evec(3,3-ishft)
                     sign2=1.d0
                     if(ix3.eq.1)then
                        sigvec(1)=sign2*evec(1,iwave-ishft)
                        sigvec(2)=sign2*evec(2,iwave-ishft)
                        sigvec(3)=sign2*evec(3,iwave-ishft)
                        if(ix2.eq.1)then
                           sigvec1(1,ix3)=sigvec(1)
                           sigvec1(2,ix3)=sigvec(2)
                           sigvec1(3,ix3)=sigvec(3)
                        else
                           call dot(sigvec,sigvec1(1,ix3),sdot)
                           if(sdot.gt.0.d0)then
                              sign2=1.d0
                           elseif(sdot.lt.0.d0)then
                              sign2=-1.d0
                           else
                              write(11,*)'Setincid error 2, sdot=',
     +                             sdot,'   ix2,ix3=',ix2,ix3
                           endif
                           sigvec(1)=sign2*evec(1,iwave-ishft)
                           sigvec(2)=sign2*evec(2,iwave-ishft)
                           sigvec(3)=sign2*evec(3,iwave-ishft)
                           sigvec1(1,ix3)=sigvec(1)
                           sigvec1(2,ix3)=sigvec(2)
                           sigvec1(3,ix3)=sigvec(3)
                        endif
                        disp1=sigvec(1)
                        disp2=sigvec(2)
                        disp3=sigvec(3)
                     else
                        call dot(sigvec,evec(1,iwave-ishft),sdot)
                        if(sdot.gt.0.d0)then
                           sign3=1.d0
                        elseif(sdot.lt.0.d0)then
                           sign3=-1.d0
                        else
                           write(11,*)'Setincid error 3, sdot=',sdot,
     +                          '   ix2,ix3=',ix2,ix3
                        endif
                        sigvec(1)=sign3*evec(1,iwave-ishft)
                        sigvec(2)=sign3*evec(2,iwave-ishft)
                        sigvec(3)=sign3*evec(3,iwave-ishft)
                        if(ix2.gt.1)then
                           call dot(sigvec,sigvec1(1,ix3),sdot)
                           if((itropic.eq.1).and.(sdot.lt.0.d0))then
                              write(11,*)'Conical point ix2,ix3=',
     +                             ix2,ix3
                           elseif(sdot.eq.0.d0)then
                              write(11,*)'Setincid error 4, sdot=',
     +                             sdot,'   ix2,ix3=',ix2,ix3
                           endif
                        endif
                        disp1=sigvec(1)
                        disp2=sigvec(2)
                        disp3=sigvec(3)
                        sigvec1(1,ix3)=sigvec(1)
                        sigvec1(2,ix3)=sigvec(2)
                        sigvec1(3,ix3)=sigvec(3)
                     endif

                     disp1=btmp*evec(1,2-ishft)+ctmp*evec(1,3-ishft)
                     disp2=btmp*evec(2,2-ishft)+ctmp*evec(2,3-ishft)
                     disp3=btmp*evec(3,2-ishft)+ctmp*evec(3,3-ishft)

c     Moved from earlier ipol.eq.0 statement in nreal=2 if statement
c                     elseif(ipol.eq.1)then
c                        disp1=0.0
c                        disp2=cos(s1phi)
c                        disp3=sin(s1phi)
c                  endif

c     Displacement now set for S1 or S2 wave.
                  elseif(nreal.eq.1)then
                     write(11,*)'Danger: in setincid nreal=',nreal,' at 
     +                    (ix2,ix3)=',ix2,ix3,' (stopping).'
                     stop
                  endif

c     Displacement vector now set - output and initialisation.

                  write(20,*)'ix2,ix3,tinc,nreal'
                  write(20,1000)ix2,ix3,tinc(ix2,ix3),nreal
                  write(20,*)'eval,p2inc,p3inc,disp1,disp2,disp3'
                  write(20,1001)eval(iwave),p2inc(icomp,ix2,ix3),
     +                 p3inc(icomp,ix2,ix3),disp1,disp2,disp3
c                  write(20,1001)eval(1),
c     +                 p2inc(icomp,ix2,ix3),p3inc(icomp,ix2,ix3),
c     +                 evec(1,1),evec(2,1),evec(3,1)
c                  write(20,1001)eval(2),
c     +                 p2inc(icomp,ix2,ix3),p3inc(icomp,ix2,ix3),
c     +                 evec(1,2),evec(2,2),evec(3,2)
c                  write(20,1001)eval(3),
c     +                 p2inc(icomp,ix2,ix3),p3inc(icomp,ix2,ix3),
c     +                 evec(1,3),evec(2,3),evec(3,3)
 1000             format(2(i3,1x),d11.4,1x,i3)
 1001             format(6(1x,d11.4))

c     Now apply the time shift to the incoming wavelet, in
c     the frequency domain for ease ...

                  if(ix2.eq.1.and.ix3.eq.1)then 
                     ttot=dt*dble(nt)
                     omn=pi/dt
                     nom=nt/2 + 1
                     delom=2.d0*pi/ttot
                     sdelom=sngl(delom)
                     sdt=sngl(dt)
                     srpeak=sngl(rpeak)
                     ai=dcmplx(0.d0,1.d0)
c     Set up sinc**nord filter in omega domain (i.e. boxcar
c     convolved with itself nord times in the time domain,
c     the width of the boxcar is 2*mwid time samples)
                     omegai=0.0 ! frequency is real in this code 
                     if(iricker.eq.0)then
                        call scsinc(singpi,nt,nom,mwid,nord,omegai,
     +                       cfilt)
                     elseif(iricker.eq.1)then
                        call rickerw(nt,sdt,nom,sdelom,singpi,cfilt,
     +                       srpeak)
                     endif
                  endif

c     Also apply a constant overall time shift, as the smoothed 
c     boxcar is zero-phase and we prefer a causal pulse even on the 
c     first x1 plane -- an appropriate time delay is greater than
c     half the pulse width -- will actually use full pulse width ...

                  tshift=1.d0*dt*dble(mwid)*dble(nord)
                  time=tshift+tinc(ix2,ix3)
                  om=-delom

                  do iom=1,nom
                     om=om+delom
                     ctemp=cfilt(iom)
                     if(iom.gt.1)ctemp=ctemp*cdexp(ai*om*time)
c     Note the vector component mappings here 
                     uvec(1,iom,ix3,ix2,1)=disp1*ctemp 
                     uvec(2,iom,ix3,ix2,1)=disp2*ctemp
                     uvec(3,iom,ix3,ix2,1)=disp3*ctemp
c     Conjugate the upper half of the omega line
                     if(iom.gt.1.and.iom.lt.nom)then
                        uvec(1,nt+2-iom,ix3,ix2,1)=
     +                       dconjg(uvec(1,iom,ix3,ix2,1))
                        uvec(2,nt+2-iom,ix3,ix2,1)=
     +                       dconjg(uvec(2,iom,ix3,ix2,1))
                        uvec(3,nt+2-iom,ix3,ix2,1)=
     +                       dconjg(uvec(3,iom,ix3,ix2,1))
                     endif
                  enddo
               endif

 20         continue
 10      continue

      enddo

      return
      end

c***********************************************************************
      subroutine vgrid(alen,npts,beta)
c     This subroutine contructs the variable grid based on Anderson 
c     et. al.
c***********************************************************************

      implicit none
      include '../Input/narc_dp.par'

      integer npts,igrid,ipt
      real*8 alen,alen2,xbar,dxbr,ybar,fac2,fac3,beta,pi
      real*8 dxwor(nxmx)

      pi=4.d0*pi4

      if(npts.gt.nxmx)then
         write(11,*)
         write(11,*)'Error: npts = ',npts,' but nxmx =',nxmx
         write(11,*)'Stopping program within vgrid.'
         stop
      endif

      if(((npts+1)/2)*2.ne.npts+1)then
         write(11,*)
         write(11,*)'Error: npts must be odd.'
         write(11,*)'Stopping program within vgrid.'
         stop
      endif

      open(unit=28,file='./Log/vgrid.out')

      igrid=1

      if(igrid.eq.1)then          !As discussed by Anderson et al.

         alen2=alen/2.d0
         dxbr=2.d0/real(npts-1)
         fac2=(beta+1.d0)/(beta-1.d0)

c     Fill midpoint to upper end ... lower range by reflection ...

         do ipt=(npts+1)/2,npts

            xbar=dxbr*real(ipt-(npts+1)/2)

c     Note that coordinate ybar here is measured down from upper end 

            ybar=1.d0-xbar
            fac3=(1.d0-ybar)*dlog(fac2)
            fac3=dexp(fac3)

            xwor(ipt)=alen-alen2*(beta+1.d0-(beta-1.d0)*fac3)/
     +           (fac3+1.d0)

            if(ipt.eq.1)xwor(ipt)=0.d0

            if(ipt.gt.(npts+1)/2)dxwor(ipt)=xwor(ipt)-xwor(ipt-1)

            if(ipt.gt.(npts+1)/2)then

               xwor(npts+1-ipt)=alen2-(xwor(ipt)-alen2)
               dxwor(npts+1-ipt)=dxwor(ipt)

            endif

         enddo

         write(28,*)'beta=',beta
         write(28,*)
         write(28,*)'ipt, xwor, dxwor'

      elseif(igrid.eq.2)then    !Chebyshev distribution (Fornberg)

         do ipt=1,npts
            xwor(ipt)=alen*(1.d0-
     +           dcos(dble(ipt-1)*pi/(2.d0*dble((npts)/2))))/2.d0

            if(ipt.gt.1)dxwor(ipt)=xwor(ipt)-xwor(ipt-1)
         enddo

         write(28,*)'ipt, xwor, dxwor'

      endif

      do ipt=1,npts
         write(28,1000)ipt,xwor(ipt),dxwor(ipt)
      enddo
 1000 format(i3,2(e12.4))

      close(28)

      return
      end

c***********************************************************************
      subroutine scsinc(spi,nt,nomega,mwid,nord,somegai,cfilt)
c     This subroutine creates the waveform pulse:
c     nomega = nt/2 + 1 - redundant information but good as a check
c     mwid   = 0.5*(width of time-domain boxcar) (i.e. mwid=1 means a 
c              boxcar of width 2*dt)
c     nor d   = order of filter (i.e. number of times applied)
c     Copyright (c) 2000 C. J. Thomson - all rights reserved by author.
c***********************************************************************

      implicit none
      include '../Input/narc_dp.par'

      integer nt,nomega,mwid,nord,iom
      real*4 spi,somegai,stemp
      complex*8 cfilt(ntmx),arg

c     Sinc filter for complex frequency
      
      if(nomega.ne.(nt/2+1))write(11,*)'***ERROR in scsinc**'

      cfilt(1)=cmplx(1.0,0.0)

      do iom=2,nomega
         stemp=2.0*spi*real(iom-1)*real(mwid)/real(nt)
         arg=cmplx(0.0,1.0)*cmplx(stemp,somegai)
         if(stemp.ge.spi)then
            cfilt(iom)=cmplx(0.0,0.0)
         else
            cfilt(iom)=(cexp(arg)-cexp(-arg))/2.0/arg
            cfilt(iom)=cfilt(iom)**nord 
         endif
      enddo

      return
      end

c***********************************************************************
      subroutine rickerw(nt,sdt,nom,sdelom,spi,cfilt,srpeak)
c     This subroutine creates a Ricker waveform in the frequency domain.
c     nom = nt/2 + 1 - redundant information but good as a check.
c***********************************************************************

      implicit none
      include '../Input/narc_dp.par'

      integer ideriv,nt,nom,iom,iwp
      real*4 sdt,sdelom,spi,srpeak,rhz,fn,omn,fp,omp,td,tr,om
      complex*8 cfilt(ntmx),arg,ai
      
      ai=cmplx(0.0,1.0)
      ideriv=0

      rhz=2.0*spi                !Conversion from radians to Hz
      fn=1.0/(2.0*sdt)
      omn=2.0*spi*fn

c     Set peak frequency to approximately rpeak (%) of Nyquist frequency
      if(srpeak.eq.0.0)then
         fp=0.025*fn            !2.5% of Nyquist frequency (fn)
         omp=fp*rhz
      else
         fp=srpeak
         omp=fp*rhz
      endif

      iwp=int(omp/sdelom)
      if(iwp.le.2)then
         write(*,*)'Choosen wavelet peak frequency too low for IC'
         stop
      elseif(iwp.ge.nom-1)then
         write(*,*)'Choosen wavelet peak frequency too high for IC'
         stop
      endif
      td=(sqrt(6.0)/spi)/omp
      tr=td/sqrt(3.0)

      write(*,*)
      write(*,*)'             Time domain Ricker wavelet:'
      write(*,*)'             Peak freq.:',omp,' (rads) or ',fp,' Hz.'
      write(*,*)'             Nyq. freq.:',omn,' (rads) or',fn,' Hz.'
      write(*,*)'             -> fp=',(fp/fn)*100.0,' % of fn.'
      write(*,*)'             Max. freq.:',real(nom)*sdelom/rhz,' Hz.'
      write(*,*)'             Pulse width (s):',td
      write(*,*)

      do iom=1,nom
         om=real(iom-1)*sdelom
         if(iom.eq.1)then
            cfilt(iom)=cmplx(0.0,0.0)
         else
            arg=cmplx((om/omp)**2.0,0.0)
            if(ideriv.eq.0)then
               cfilt(iom)=(2.0/sqrt(spi))*arg*cexp(-arg)
            elseif(ideriv.eq.1)then
               cfilt(iom)=ai*om*(2.0/sqrt(spi))*arg*cexp(-arg)
            endif
         endif
      enddo

      return
      end

c***********************************************************************
      subroutine dot(dvec1,dvec2,sdot)
c     This subroutine calculates the scalar product.
c***********************************************************************

      implicit none
      
      real*8 sdot,dvec1(3),dvec2(3)

      sdot=dvec1(1)*dvec2(1)+dvec1(2)*dvec2(2)+dvec1(3)*dvec2(3)

      return
      end

