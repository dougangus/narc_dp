c...+....|....+....|....+....|....+....|....+....|....+....|....+....|..
c  File Extrcart.f contains the following subroutine: extrcart.
c  This version is the modified equi-spaced fd algorithm ... modified to
c  improve the boundary differences using incident lateral slownesses 
c  in a phase shifting scheme.
c
c  Last modified: November 6, 2006.
c***********************************************************************
      subroutine extrcart(nt,dt,omwin,nx1st,nx1,dx1,nx2,dx2,nx3,dx3,a6,
     +     uvec,iref,p2inc,p3inc)
c***********************************************************************

      implicit none
      include '../Input/narc_dp.par'

      integer nx1st,nx1,nx2,nx3,ix1,jx1,jx2,jx3,kx1a,kx1b,kx1c
      integer nt,it,nom,nom2,iom,kom,iomc1,iomc2,ic,jc,kc
      integer itmp,ik,iref,iout,im
      character*17 a1text
      character*4 a2text

      real*8 dt,tott,om,omk,omn,omwin,delom
      real*8 pi,dx,dx1,dx2,dx3,tx1,x1step,root2,vel
      real*8 phas,phas1,phas2,pphas1,pphas2,pphas,dpphas,slp,sls1,sls2
      real*8 a3(3,3,3,3),a3r(3,3,3,3),c11inv(3,3)
      real*8 p0(3,3),p2(3,3),p3(3,3),p22(3,3),p23(3,3),p32(3,3),p33(3,3)
      real*8 a6(6,6,nx3mx,nx2mx,nx1sto)
      real*8 slonor(3,nx3mx,nx2mx),rfilt(ntmx)
      real*8 p2inc(3,nx2mx,nx3mx),p3inc(3,nx2mx,nx3mx)

      complex*16 utmp1,utmp2,utmp3,utmp4,utmp5,utmp6
      complex*16 ai,egtmp1,egtmp2,egtmp3,u3bar,duvec
      complex*16 cp0(3,3),cp2(3,3),cp3(3,3),cp22(3,3),cp23(3,3),
     +     cp32(3,3),cp33(3,3),diag(3,3)
      complex*16 uvec(3,ntmx,nx3mx,nx2mx,3)
      complex*16 ctmp1(4*nx2mx),ctmp2(4*nx2mx),ctmp3(4*nx2mx)

c  ** added by JW
      character tfname*128 ! temporary filename
      integer ifindex ! index count into direct access file

      if(nx2.gt.nxmx.or.nx3.gt.nxmx)then
         write(*,*)
         write(*,*)'Error: nx2,nx3 = ',nx2,nx3,' but nxmx =',nxmx
         write(*,*)'Stopping program in extrcart.'
         stop
      endif

      if(iwave.eq.1) vel=vfast
      if(iwave.ne.1) vel=vslow

      a1text='             ix1:'
      a2text=' x1:'

      pi=4.d0*pi4
      ai=dcmplx(0.d0,1.d0)
      omn=4.d0*pi/dt
      tott=dreal(nt)*dt
      delom=2.d0*pi/tott
      nom=nt/2+1
      nom2=int(dreal(nom-1)*omwin)+1

      iomc1=1
      iomc2=1
      call rtaper(nom2,iomc1,iomc2,rfilt)

      if(imtrick.eq.1)then
         itmp=int(nx1/ncx1)
      endif

c  ** intialise the IFINDEX parameter 
      ifindex = 1


      do 10 ix1=1,nx1
         if(ix1.eq.1)then
            kx1a=1              !first-order extrapolation for 1 to 2
            kx1b=1              !plane
            kx1c=2
            x1step=dx1    
         elseif(ix1.ge.2)then
            kx1a=1              !second-order extrapolation for the rest
            kx1b=2
            kx1c=3
            x1step=2.d0*dx1
         endif
         call edo(ix1,jx1,nx1st,nx1,nx2,nx3,dx1,dx2,dx3,a6)
c     Loop over inner part of grid (not outer grid points)
         do 20 jx2=1,nx2
            do 30 jx3=1,nx3
               if(inhomog.eq.0.or.iuse.eq.1)then
                  if(ix1.eq.1.and.jx2.eq.1.and.jx3.eq.1)then
                     call aijkl(a6(1,1,jx3,jx2,jx1),a3r)
                     call propagat(a3r,c11inv,p0,p2,p3,p22,p23,p32,p33)
                     do ic=1,3
                        slonor(ic,jx3,jx2)=p0(ic,ic)
                     enddo
                  else
                     do ic=1,3
                        slonor(ic,jx3,jx2)=slonor(ic,1,1)
                     enddo
                  endif
               elseif(inhomog.eq.1)then
                  call aijkl(a6(1,1,jx3,jx2,jx1),a3r)
                  call propagat(a3r,c11inv,p0,p2,p3,p22,p23,p32,p33)
                  do ic=1,3
                     slonor(ic,jx3,jx2)=p0(ic,ic)
                  enddo
               endif
               do ic=1,3
                  do jc=1,3
                     if(iref.eq.1)then
                        if(ic.eq.jc)then
                           p0(ic,jc)=0.d0
                        else
                           diag(ic,jc)=dcmplx(0.d0,0.d0)
                        endif
                     endif
                     if(jx2.ne.1.and.jx2.ne.nx2.and.
     +                    jx3.ne.1.and.jx3.ne.nx3)then
                        cp0(ic,jc)=dcmplx(p0(ic,jc),0.d0)
                        cp2(ic,jc)=dcmplx(p2(ic,jc),0.d0)
                        cp3(ic,jc)=dcmplx(p3(ic,jc),0.d0)
                        cp22(ic,jc)=dcmplx(p22(ic,jc),0.d0)
                        cp33(ic,jc)=dcmplx(p33(ic,jc),0.d0)
                        cp23(ic,jc)=dcmplx(p23(ic,jc),0.d0)
                        cp32(ic,jc)=dcmplx(p32(ic,jc),0.d0)
                     endif
                  enddo
               enddo
               if(jx2.eq.1.or.jx2.eq.nx2.or.jx3.eq.1.or.jx3.eq.nx3) 
     +              goto 30
               do 40 iom=1,nom2
                  om=dreal(iom-1)*delom
                  do 50 ic=1,3
                     if(iom.eq.1)then
                        uvec(ic,iom,jx3,jx2,kx1c)=
     +                       uvec(ic,iom,jx3,jx2,kx1a)
                        go to 50
                     endif
                     utmp1=dcmplx(0.d0,0.d0)
                     utmp2=dcmplx(0.d0,0.d0)
                     utmp3=dcmplx(0.d0,0.d0)
                     utmp4=dcmplx(0.d0,0.d0)
                     utmp5=dcmplx(0.d0,0.d0)
                     utmp6=dcmplx(0.d0,0.d0)
                     do 60 jc=1,3
                        utmp1=utmp1+cp0(ic,jc)*uvec(jc,iom,jx3,jx2,kx1b)
                        utmp2=utmp2+cp2(ic,jc)*
     +                       (uvec(jc,iom,jx3,jx2+1,kx1b)-
     +                       uvec(jc,iom,jx3,jx2-1,kx1b))
                        utmp3=utmp3+cp3(ic,jc)*
     +                       (uvec(jc,iom,jx3+1,jx2,kx1b)-
     +                       uvec(jc,iom,jx3-1,jx2,kx1b))
                        utmp4=utmp4+cp22(ic,jc)*
     +                       (uvec(jc,iom,jx3,jx2+1,kx1b)-2.d0*
     +                       uvec(jc,iom,jx3,jx2,kx1b)+
     +                       uvec(jc,iom,jx3,jx2-1,kx1b))
                        utmp5=utmp5+cp33(ic,jc)*
     +                       (uvec(jc,iom,jx3+1,jx2,kx1b)-2.d0*
     +                       uvec(jc,iom,jx3,jx2,kx1b)+
     +                       uvec(jc,iom,jx3-1,jx2,kx1b))
                        utmp6=utmp6+(cp23(ic,jc)+cp32(ic,jc))*
     +                       (uvec(jc,iom,jx3+1,jx2+1,kx1b)-
     +                       uvec(jc,iom,jx3-1,jx2+1,kx1b)-
     +                       uvec(jc,iom,jx3+1,jx2-1,kx1b)+
     +                       uvec(jc,iom,jx3-1,jx2-1,kx1b))
 60                  continue 
c     Step in x1
                     duvec=x1step*(ai*om*utmp1+
     +                    utmp2/(2.d0*dx2)+utmp3/(2.d0*dx3)+
     +                    (-ai/om)*(utmp4/(dx2*dx2)+utmp5/(dx3*dx3)+
     +                    utmp6/(4.d0*dx2*dx3)))
                     if(iref.eq.1)then
                        uvec(ic,iom,jx3,jx2,kx1c)=duvec
                     elseif(iref.eq.0)then
                        uvec(ic,iom,jx3,jx2,kx1c)=
     +                       uvec(ic,iom,jx3,jx2,kx1a)+duvec
                     endif
 50               continue
 40            continue
 30         continue
 20      continue

c     Now pass over edge and corner grid points by linear extrapolation 
c     from the inner grid points 
         kom=10
         omk=delom*dreal(kom-1)
         kc=iwave
         root2=dsqrt(2.d0)
         do 61 jx2=1,nx2
            do 62 jx3=1,nx3
               do 63 iom=1,nom2
                  om=dreal(iom-1)*delom
                  do 64 ic=1,3
                     if(jx2.eq.1.and.jx3.eq.1)then
                        ik=1
                        dx=dx2
                        if(ic.eq.1.and.iom.eq.1.and.ix1.eq.1)then
                           ctmp1(ik)=uvec(kc,kom,jx3+1,jx2+1,kx1a)
                           ctmp2(ik)=uvec(kc,kom,jx3+2,jx2+2,kx1a)
                           ctmp3(ik)=uvec(kc,kom,jx3+3,jx2+3,kx1a)
                        endif
                        if(ctmp1(ik).eq.ctmp2(ik).or.ctmp1(ik).eq.
     +                       dcmplx(0.d0,0.d0))then
                           if(ic.eq.1.and.iom.eq.1)then
                              egtmp3=dcmplx(1.d0,0.d0)
                           endif
                           uvec(ic,iom,jx3,jx2,kx1c)=
     +                          uvec(ic,iom,jx3+1,jx2+1,kx1c)*egtmp3
                        else
                           if(ic.eq.1.and.iom.eq.1)then
                              egtmp1=ctmp2(ik)/ctmp1(ik)
                              egtmp2=ctmp3(ik)/ctmp2(ik)
                              phas1=datan2(dimag(egtmp1),dreal(egtmp1))
                              phas2=datan2(dimag(egtmp2),dreal(egtmp2))
                           endif
                           if(iom.eq.1)then
                              egtmp3=2.d0*egtmp1-egtmp2
                           else
                              pphas1=phas1/(omk*dx)
                              pphas2=phas2/(omk*dx)
                              dpphas=pphas1-pphas2
                              pphas=pphas1+dpphas
                              phas=pphas*om*dx
                              phas=(p2inc(iwave,2,2)+
     +                             p3inc(iwave,2,2))*dx*om
                              egtmp3=dcmplx(dcos(phas),-dsin(phas))
                           endif
                           u3bar=uvec(ic,iom,jx3+2,jx2+2,kx1c)*egtmp3
                           egtmp2=u3bar*
     +                          dconjg(uvec(ic,iom,jx3+1,jx2+1,kx1c))
                           phas=datan2(imagpart(egtmp2),
     +                          realpart(egtmp2))
                           egtmp2=dcmplx(dcos(phas),-dsin(phas))
                           u3bar=u3bar*egtmp2
                           u3bar=2.d0*uvec(ic,iom,jx3+1,jx2+1,kx1c)-
     +                          u3bar
                           uvec(ic,iom,jx3,jx2,kx1c)=u3bar*egtmp3
     +                          *egtmp2
                        endif
                     elseif(jx2.eq.1.and.jx3.eq.nx3)then
                        ik=nx3
                        dx=dx2
                        if(ic.eq.1.and.iom.eq.1.and.ix1.eq.1)then
                           ctmp1(ik)=uvec(kc,kom,jx3-1,jx2+1,kx1a)
                           ctmp2(ik)=uvec(kc,kom,jx3-2,jx2+2,kx1a)
                           ctmp3(ik)=uvec(kc,kom,jx3-3,jx2+3,kx1a)
                        endif
                        if(ctmp1(ik).eq.ctmp2(ik).or.ctmp1(ik).eq.
     +                       dcmplx(0.d0,0.d0))then
                           if(ic.eq.1.and.iom.eq.1)then
                              egtmp3=dcmplx(1.d0,0.d0)
                           endif
                           uvec(ic,iom,jx3,jx2,kx1c)=
     +                          uvec(ic,iom,jx3-1,jx2+1,kx1c)*egtmp3
                        else
                           if(ic.eq.1.and.iom.eq.1)then
                              egtmp1=ctmp2(ik)/ctmp1(ik)
                              egtmp2=ctmp3(ik)/ctmp2(ik)
                              phas1=datan2(dimag(egtmp1),dreal(egtmp1))
                              phas2=datan2(dimag(egtmp2),dreal(egtmp2))
                           endif
                           if(iom.eq.1)then
                              egtmp3=2.d0*egtmp1-egtmp2
                           else
                              pphas1=phas1/(omk*dx)
                              pphas2=phas2/(omk*dx)
                              dpphas=pphas1-pphas2
                              pphas=pphas1+dpphas
                              phas=pphas*om*dx
                              phas=(p2inc(iwave,2,nx3-1)-
     +                             p3inc(iwave,2,nx3-1))*dx*om
                              egtmp3=dcmplx(dcos(phas),-dsin(phas))
                           endif
                           u3bar=uvec(ic,iom,jx3-2,jx2+2,kx1c)*egtmp3
                           egtmp2=u3bar*
     +                            dconjg(uvec(ic,iom,jx3-1,jx2+1,kx1c))
                           phas=datan2(imagpart(egtmp2),
     +                          realpart(egtmp2))
                           egtmp2=dcmplx(dcos(phas),-dsin(phas))
                           u3bar=u3bar*egtmp2
                           u3bar=2.d0*uvec(ic,iom,jx3-1,jx2+1,kx1c)-
     +                          u3bar
                           uvec(ic,iom,jx3,jx2,kx1c)=u3bar*egtmp3
     +                          *egtmp2
                        endif
                     elseif(jx2.eq.nx2.and.jx3.eq.1)then
                        ik=nx3+nx2
                        dx=dx2
                        if(ic.eq.1.and.iom.eq.1.and.ix1.eq.1)then
                            ctmp1(ik)=uvec(kc,kom,jx3+1,jx2-1,kx1a)
                            ctmp2(ik)=uvec(kc,kom,jx3+2,jx2-2,kx1a)
                            ctmp3(ik)=uvec(kc,kom,jx3+3,jx2-3,kx1a)
                        endif
                        if(ctmp1(ik).eq.ctmp2(ik).or.ctmp1(ik).eq.
     +                       dcmplx(0.d0,0.d0))then
                           if(ic.eq.1.and.iom.eq.1)then
                              egtmp3=dcmplx(1.d0,0.d0)
                           endif
                           uvec(ic,iom,jx3,jx2,kx1c)=
     +                          uvec(ic,iom,jx3+1,jx2-1,kx1c)*egtmp3
                        else
                           if(ic.eq.1.and.iom.eq.1)then
                              egtmp1=ctmp2(ik)/ctmp1(ik)
                              egtmp2=ctmp3(ik)/ctmp2(ik)
                              phas1=datan2(dimag(egtmp1),dreal(egtmp1))
                              phas2=datan2(dimag(egtmp2),dreal(egtmp2))
                           endif
                           if(iom.eq.1)then
                              egtmp3=2.d0*egtmp1-egtmp2
                           else
                              pphas1=phas1/(omk*dx)
                              pphas2=phas2/(omk*dx)
                              dpphas=pphas1-pphas2
                              pphas=pphas1+dpphas
                              phas=pphas*om*dx
                              phas=(-p2inc(iwave,nx2-1,2)+
     +                             p3inc(iwave,nx2-1,2))*dx*om
                              egtmp3=dcmplx(dcos(phas),-dsin(phas))
                           endif
                           u3bar=uvec(ic,iom,jx3+2,jx2-2,kx1c)*egtmp3
                           egtmp2=u3bar*
     +                          dconjg(uvec(ic,iom,jx3+1,jx2-1,kx1c))
                           phas=datan2(imagpart(egtmp2),
     +                          realpart(egtmp2))
                           egtmp2=dcmplx(dcos(phas),-dsin(phas))
                           u3bar=u3bar*egtmp2
                           u3bar=2.d0*uvec(ic,iom,jx3+1,jx2-1,kx1c)-
     +                          u3bar 
                           uvec(ic,iom,jx3,jx2,kx1c)=u3bar*egtmp3
     +                          *egtmp2
                        endif
                     elseif(jx2.eq.nx2.and.jx3.eq.nx3)then
                        ik=2*nx3+nx2
                        dx=dx2
                        if(ic.eq.1.and.iom.eq.1.and.ix1.eq.1)then
                            ctmp1(ik)=uvec(kc,kom,jx3-1,jx2-1,kx1a)
                            ctmp2(ik)=uvec(kc,kom,jx3-2,jx2-2,kx1a)
                            ctmp3(ik)=uvec(kc,kom,jx3-3,jx2-3,kx1a)
                        endif
                        if(ctmp1(ik).eq.ctmp2(ik).or.ctmp1(ik).eq.
     +                       dcmplx(0.d0,0.d0))then
                           if(ic.eq.1.and.iom.eq.1)then
                              egtmp3=dcmplx(1.d0,0.d0)
                           endif
                           uvec(ic,iom,jx3,jx2,kx1c)=
     +                          uvec(ic,iom,jx3-1,jx2-1,kx1c)*egtmp3
                        else
                           if(ic.eq.1.and.iom.eq.1)then
                              egtmp1=ctmp2(ik)/ctmp1(ik)
                              egtmp2=ctmp3(ik)/ctmp2(ik)
                              phas1=datan2(dimag(egtmp1),dreal(egtmp1))
                              phas2=datan2(dimag(egtmp2),dreal(egtmp2))
                           endif
                           if(iom.eq.1)then
                              egtmp3=2.d0*egtmp1-egtmp2
                           else
                              pphas1=phas1/(omk*dx)
                              pphas2=phas2/(omk*dx)
                              dpphas=pphas1-pphas2
                              pphas=pphas1+dpphas
                              phas=pphas*om*dx
                              phas=-(p2inc(iwave,nx2-1,nx3-1)
     +                             +p3inc(iwave,nx2-1,nx3-1))*dx*om
                              egtmp3=dcmplx(dcos(phas),-dsin(phas))
                           endif
                           u3bar=uvec(ic,iom,jx3-2,jx2-2,kx1c)*egtmp3
                           egtmp2=u3bar*
     +                          dconjg(uvec(ic,iom,jx3-1,jx2-1,kx1c))
                           phas=datan2(imagpart(egtmp2),
     +                          realpart(egtmp2))
                           egtmp2=dcmplx(dcos(phas),-dsin(phas))
                           u3bar=u3bar*egtmp2
                           u3bar=2.d0*uvec(ic,iom,jx3-1,jx2-1,kx1c)-
     +                          u3bar
                           uvec(ic,iom,jx3,jx2,kx1c)=u3bar*egtmp3
     +                          *egtmp2
                        endif
                     elseif(jx2.eq.1.and.(jx3.gt.1.and.jx3.lt.nx3))then
                        ik=jx3
                        dx=dx2
                        if(ic.eq.1.and.iom.eq.1.and.ix1.eq.1)then
                           ctmp1(ik)=uvec(kc,kom,jx3,jx2+1,kx1a)
                           ctmp2(ik)=uvec(kc,kom,jx3,jx2+2,kx1a)
                           ctmp3(ik)=uvec(kc,kom,jx3,jx2+3,kx1a)
                        endif
                        if(ctmp1(ik).eq.ctmp2(ik).or.ctmp1(ik).eq.
     +                       dcmplx(0.d0,0.d0))then
                           if(ic.eq.1.and.iom.eq.1)then
                              egtmp3=dcmplx(1.d0,0.d0)
                           endif
                           uvec(ic,iom,jx3,jx2,kx1c)=
     +                          uvec(ic,iom,jx3,jx2+1,kx1c)*egtmp3
                        else
                           if(ic.eq.1.and.iom.eq.1)then
                              egtmp1=ctmp2(ik)/ctmp1(ik)
                              egtmp2=ctmp3(ik)/ctmp2(ik)
                              phas1=datan2(dimag(egtmp1),dreal(egtmp1))
                              phas2=datan2(dimag(egtmp2),dreal(egtmp2))
                           endif
                           if(iom.eq.1)then
                              egtmp3=2.d0*egtmp1-egtmp2
                           else
                              pphas1=phas1/(omk*dx)
                              pphas2=phas2/(omk*dx)
                              dpphas=pphas1-pphas2
                              pphas=pphas1+dpphas
                              phas=pphas*om*dx
                              phas=p2inc(iwave,2,jx3)*dx*om
                              egtmp3=dcmplx(dcos(phas),-dsin(phas))
                           endif
                           u3bar=uvec(ic,iom,jx3,jx2+2,kx1c)*egtmp3
                           egtmp2=u3bar*
     +                          dconjg(uvec(ic,iom,jx3,jx2+1,kx1c))
                           phas=datan2(imagpart(egtmp2),
     +                          realpart(egtmp2))
                           egtmp2=dcmplx(dcos(phas),-dsin(phas))
                           u3bar=u3bar*egtmp2
                           u3bar=2.d0*uvec(ic,iom,jx3,jx2+1,kx1c)-
     +                          u3bar 
                           uvec(ic,iom,jx3,jx2,kx1c)=u3bar*egtmp3
     +                          *egtmp2
                        endif
                     elseif((jx2.gt.1.and.jx2.lt.nx2).and.jx3.eq.1)then
                        ik=nx3+jx2-1
                        dx=dx2
                        if(ic.eq.1.and.iom.eq.1.and.ix1.eq.1)then
                           ctmp1(ik)=uvec(kc,kom,jx3+1,jx2,kx1a)
                           ctmp2(ik)=uvec(kc,kom,jx3+2,jx2,kx1a)
                           ctmp3(ik)=uvec(kc,kom,jx3+3,jx2,kx1a)
                        endif
                        if(ctmp1(ik).eq.ctmp2(ik).or.ctmp1(ik).eq.
     +                       dcmplx(0.d0,0.d0))then
                           if(ic.eq.1.and.iom.eq.1)then
                              egtmp3=dcmplx(1.d0,0.d0)
                           endif
                           uvec(ic,iom,jx3,jx2,kx1c)=
     +                          uvec(ic,iom,jx3+1,jx2,kx1c)*egtmp3
                        else
                           if(ic.eq.1.and.iom.eq.1)then
                              egtmp1=ctmp2(ik)/ctmp1(ik)
                              egtmp2=ctmp3(ik)/ctmp2(ik)
                              phas1=datan2(dimag(egtmp1),dreal(egtmp1))
                              phas2=datan2(dimag(egtmp2),dreal(egtmp2))
                           endif
                           if(iom.eq.1)then
                              egtmp3=2.d0*egtmp1-egtmp2
                           else
                              pphas1=phas1/(omk*dx)
                              pphas2=phas2/(omk*dx)
                              dpphas=pphas1-pphas2
                              pphas=pphas1+dpphas
                              phas=pphas*om*dx
                              phas=p3inc(iwave,jx2,2)*dx*om
                              egtmp3=dcmplx(dcos(phas),-dsin(phas))
                           endif
                           u3bar=uvec(ic,iom,jx3+2,jx2,kx1c)*egtmp3
                           egtmp2=u3bar*
     +                          dconjg(uvec(ic,iom,jx3+1,jx2,kx1c))
                           phas=datan2(imagpart(egtmp2),
     +                          realpart(egtmp2))
                           egtmp2=dcmplx(dcos(phas),-dsin(phas))
                           u3bar=u3bar*egtmp2
                           u3bar=2.d0*uvec(ic,iom,jx3+1,jx2,kx1c)-
     +                          u3bar 
                           uvec(ic,iom,jx3,jx2,kx1c)=u3bar*egtmp3
     +                          *egtmp2
                        endif
                     elseif(jx2.eq.nx2.and.(jx3.gt.1.and.jx3.lt.nx3))
     +                       then
                        ik=nx3+nx2+jx3-1
                        dx=dx2
                        if(ic.eq.1.and.iom.eq.1.and.ix1.eq.1)then
                           ctmp1(ik)=uvec(kc,kom,jx3,jx2-1,kx1a)
                           ctmp2(ik)=uvec(kc,kom,jx3,jx2-2,kx1a)
                           ctmp3(ik)=uvec(kc,kom,jx3,jx2-3,kx1a)
                        endif
                        if(ctmp1(ik).eq.ctmp2(ik).or.ctmp1(ik).eq.
     +                       dcmplx(0.d0,0.d0))then
                           if(ic.eq.1.and.iom.eq.1)then
                              egtmp3=dcmplx(1.d0,0.d0)
                           endif
                           uvec(ic,iom,jx3,jx2,kx1c)=
     +                          uvec(ic,iom,jx3,jx2-1,kx1c)*egtmp3
                        else
                           if(ic.eq.1.and.iom.eq.1)then
                              egtmp1=ctmp2(ik)/ctmp1(ik)
                              egtmp2=ctmp3(ik)/ctmp2(ik)
                              phas1=datan2(dimag(egtmp1),dreal(egtmp1))
                              phas2=datan2(dimag(egtmp2),dreal(egtmp2))
                           endif
                           if(iom.eq.1)then
                              egtmp3=2.d0*egtmp1-egtmp2
                           else
                              pphas1=phas1/(omk*dx)
                              pphas2=phas2/(omk*dx)
                              dpphas=pphas1-pphas2
                              pphas=pphas1+dpphas
                              phas=pphas*om*dx
                              phas=-p2inc(iwave,nx2-1,jx3)*dx*om
                              egtmp3=dcmplx(dcos(phas),-dsin(phas))
                           endif
                           u3bar=uvec(ic,iom,jx3,jx2-2,kx1c)*egtmp3
                           egtmp2=u3bar*
     +                          dconjg(uvec(ic,iom,jx3,jx2-1,kx1c))
                           phas=datan2(imagpart(egtmp2),
     +                          realpart(egtmp2))
                           egtmp2=dcmplx(dcos(phas),-dsin(phas))
                           u3bar=u3bar*egtmp2
                           u3bar=2.d0*uvec(ic,iom,jx3,jx2-1,kx1c)-u3bar 
                           uvec(ic,iom,jx3,jx2,kx1c)=u3bar*egtmp3
     +                          *egtmp2
                        endif
                     elseif((jx2.gt.1.and.jx2.lt.nx2).and.jx3.eq.nx3)
     +                       then
                        ik=2*nx3+nx2+jx2-1
                        dx=dx2
                        if(ic.eq.1.and.iom.eq.1.and.ix1.eq.1)then
                           ctmp1(ik)=uvec(kc,kom,jx3-1,jx2,kx1a)
                           ctmp2(ik)=uvec(kc,kom,jx3-2,jx2,kx1a)
                           ctmp3(ik)=uvec(kc,kom,jx3-3,jx2,kx1a)
                        endif
                        if(ctmp1(ik).eq.ctmp2(ik).or.ctmp1(ik).eq.
     +                       dcmplx(0.d0,0.d0))then
                           if(ic.eq.1.and.iom.eq.1)then
                              egtmp3=dcmplx(1.d0,0.d0)
                           endif
                           uvec(ic,iom,jx3,jx2,kx1c)=
     +                          uvec(ic,iom,jx3-1,jx2,kx1c)*egtmp3
                        else
                           if(ic.eq.1.and.iom.eq.1)then
                              egtmp1=ctmp2(ik)/ctmp1(ik)
                              egtmp2=ctmp3(ik)/ctmp2(ik)
                              phas1=datan2(dimag(egtmp1),dreal(egtmp1))
                              phas2=datan2(dimag(egtmp2),dreal(egtmp2))
                           endif
                           if(iom.eq.1)then
                              egtmp3=2.d0*egtmp1-egtmp2
                           else
                              pphas1=phas1/(omk*dx)
                              pphas2=phas2/(omk*dx)
                              dpphas=pphas1-pphas2
                              pphas=pphas1+dpphas
                              phas=pphas*om*dx
                              phas=-p3inc(iwave,jx2,nx3-1)*dx*om
                              egtmp3=dcmplx(dcos(phas),-dsin(phas))
                           endif
                           u3bar=uvec(ic,iom,jx3-2,jx2,kx1c)*egtmp3
                           egtmp2=u3bar*
     +                          dconjg(uvec(ic,iom,jx3-1,jx2,kx1c))
                           phas=datan2(imagpart(egtmp2),
     +                          realpart(egtmp2))
                           egtmp2=dcmplx(dcos(phas),-dsin(phas))
                           u3bar=u3bar*egtmp2
                           u3bar=2.d0*uvec(ic,iom,jx3-1,jx2,kx1c)-u3bar 
                           uvec(ic,iom,jx3,jx2,kx1c)=u3bar*egtmp3
     +                          *egtmp2
                        endif
                     endif
 64               continue
 63            continue
 62         continue
 61      continue

c     Now step forward into next x1 plane if reference phase is removed

         if(iref.eq.1)then
            do 67 jx2=1,nx2
               do 68 jx3=1,nx3
                  slp=slonor(1,jx3,jx2) 
                  sls1=slonor(2,jx3,jx2)
                  sls2=slonor(3,jx3,jx2)
                  do 69 iom=1,nom2
                     om=dreal(iom-1)*delom
                     diag(1,1)=cdexp(ai*om*slp*dx1) 
                     diag(2,2)=cdexp(ai*om*sls1*dx1) 
                     diag(3,3)=cdexp(ai*om*sls2*dx1)
                     do 70 ic=1,3
c     Remove reference phase from uvec in back plane
                        if(ix1.gt.1)then
                           uvec(ic,iom,jx3,jx2,kx1a)=
     +                          uvec(ic,iom,jx3,jx2,kx1a)
     +                          *diag(ic,ic)
                        endif
                        if(iom.eq.1)then
                           uvec(ic,iom,jx3,jx2,kx1c)=
     +                          uvec(ic,iom,jx3,jx2,kx1a)
                           go to 70
                        endif
c     Step forward  
                        uvec(ic,iom,jx3,jx2,kx1c)=
     +                       uvec(ic,iom,jx3,jx2,kx1a)
     +                       +uvec(ic,iom,jx3,jx2,kx1c)
c     Re-introduce reference phase into new plane
                           uvec(ic,iom,jx3,jx2,kx1c)=
     +                          uvec(ic,iom,jx3,jx2,kx1c)
     +                          *diag(ic,ic)
 70                  continue
 69               continue
 68            continue
 67         continue
         endif
 
         iout=0
         if(ix1.eq.1)then
            iout=1
            kx1c=1
         elseif(ix1.gt.1.and.ix1.lt.nx1)then
            do im=1,inumber
               if(ix1.eq.istart+(im-1)*increment)then
                  iout=1
               endif
            enddo
         elseif(ix1.eq.nx1)then
            iout=1
         endif
         tx1=x1o+dreal(ix1-1)*dx1

c-----------------------------------------------------------------------
c     ** rewritten output by James Wookey: Use direct access files for
c     ** quicker ways into the output
c-----------------------------------------------------------------------
         if(iout.eq.1)then
c        ** screen messages
            if(ix1.eq.1)write(*,666)a1text,ix1,a2text,tx1
            if(ix1.ne.1)write(*,666)'',ix1,'',tx1

            if(ix1.eq.1) then
c           ** first time around, open files and write grid and 
c           ** time parameters to header file
               write(tfname,'(a,a)') 
     +            waveout(1:index(waveout,' ')-1),'.bin'

               open(22,file=tfname,form='unformatted',
     +                        access='direct',recl=16)
               
               write(tfname,'(a,a)') 
     +            waveout(1:index(waveout,' ')-1),'.hdr'
               open(23,file=tfname)

               write(23,*) inumber,istart,increment,nx1
               write(23,*) dx1,x1o
               write(23,*) nx1,nx2,nx3
               write(23,*) nt,dt,nom,nom2
            endif

c        ** write out the output current plane and starting index
            write(23,*) ix1,ifindex
            
            do jx2=1,nx2
               do jx3=1,nx3
                  do it=1,nt
                     write(22,rec=ifindex) uvec(1,it,jx3,jx2,kx1c)
                     write(22,rec=ifindex+1) uvec(2,it,jx3,jx2,kx1c)
                     write(22,rec=ifindex+2) uvec(3,it,jx3,jx2,kx1c)
                     ifindex = ifindex + 3
                  enddo
               enddo
            enddo

         endif
c-----------------------------------------------------------------------

         if(ix1.gt.1.and.ix1.lt.nx1)then
            call uvecshft(uvec,nt,nom,nx3,nx2)
         endif

 10   continue

 666  format(1x,a17,1x,i6,1x,a4,1x,f9.2,' m.')

      close(26)

      return
      end
