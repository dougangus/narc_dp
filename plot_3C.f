c........+.........+.........+.........+.........+.........+.........+..
c  Program Plot_3C: 
c  The program plots waveform results from the narrow angle 
c  approximation program Oneway_Omega.  It contains the following 9
c  subroutines: Four1, Taper, Gksinit, Dfclrs, Assorted, Bndary2, 
c  Mytext1, Textl, Mytext2, and Grclse2.
c
c  Copyright (c) 2006 D.A. Angus.
c  All rights reserved by the author(s).
c  Last modified: April 19, 2012.
c***********************************************************************
      program plot_3C
c***********************************************************************

      implicit real*4 (a-h,o-z)

      parameter (ntmx=257,ntrmx=257,nx1mx=5)
      parameter (nx2mx=201,nx3mx=201)
      parameter (ntracwmx=201)

      integer unit1,unit2
      character*5 a1text,a2text,a3text,a4text

      real*8 uvec(3,ntmx,ntrmx),rfilt(ntmx)
      real*8 dt,pi,tott,om,delom,dx1,x1o,omc,xshift
      real*8 tbg(nx2mx,nx3mx),tnd(nx2mx,nx3mx)

      dimension xfun(ntmx),yfun(ntmx)
      dimension tbeg(ntrmx),tend(ntrmx)

      complex*16 uvc(3,ntmx,nx3mx,nx2mx),work1(ntmx),work2(ntmx),
     +     work3(ntmx)
      complex*16 cxshift,ai

      character*80 xtext,ytext,ctext
      character*6 chfmx,chfmy
      character*1 ansm,answer2,answer3

      dimension nta(nx2mx,nx3mx)
      dimension ix2pt(ntracwmx),ix3pt(ntracwmx)
      dimension iplanes(100)

      common/wcord/wxmax,wxmin,wymax,wymin,wdx,wdy
      common/vcord/xmax,xmin,ymax,ymin

      data unit1/14/,unit2/15/

      a1text='Tag:'
      a2text=' ix1:'
      a3text=' x1:'
      a4text=' m.'

      pi=4.d0*datan(1.d0)
      ai=dcmplx(0.d0,1.d0)

      ihidden=0
      iprompt=1

c     Decide which plane to view
      open(24,file='./Output/wavefld.out',form='unformatted')
c      open(24,file='./FP_Study_Results/wavefld.out',form='unformatted')
      read(24)inumber,istart,increment,nx1
      read(24)dx1,x1o

      do im=1,inumber+2
         if(im.eq.1)then
            iplanes(im)=1
         elseif(im.eq.inumber+2)then
            iplanes(im)=nx1
         else
            iplanes(im)=istart+(im-2)*increment
         endif
         if(im.gt.1)then
            if(iplanes(im).eq.iplanes(im-1))then
               ilast=im-1
               goto 11
            elseif(iplanes(im).lt.iplanes(im-1))then
               ilast=im-1
               goto 11
            elseif(iplanes(im).gt.nx1)then
               ilast=im-1
               goto 11
            elseif(iplanes(im).eq.nx1)then
               ilast=im
               goto 11
            else
               ilast=im
            endif
         endif
      enddo

 11   write(*,*)
      write(*,*)'The following planes are stored within input file:'
      do im=1,ilast
c         write(*,*)'Tag:',im,' ix1:',iplanes(im),' x1:',
c     +        x1o+real(iplanes(im)-1)*dx1,' m.'
         write(*,666)a1text,im,a2text,iplanes(im),a3text,
     +        x1o+dble(iplanes(im)-1)*dx1,a4text
      enddo
 666  format(1x,a,i6,a,i6,a,f9.2,a)
      write(*,*)
      write(*,*)'Select plane to view based on tag integer value.'
      read(*,*)iplane

c     Read frequency-domain wavefield

      do ip=1,iplane
         read(24)jx1,nx1,nx2,nx3
         read(24)nt,dt,nom,nom2
         if(ip.eq.iplane.and.iprompt.eq.1)then
            write(*,*)'Plane being viewed is:',jx1,' .'
         endif
         do ix2=1,nx2
            do ix3=1,nx3
               read(24)jx2,jx3,nta(ix2,ix3),tbg(ix2,ix3),tnd(ix2,ix3)
               do it=1,nt
                  read(24)uvc(1,it,ix3,ix2),uvc(2,it,ix3,ix2),
     +                 uvc(3,it,ix3,ix2)
               enddo
            enddo
         enddo
      enddo
      write(*,*)'Finished reading wavefield array.'
      close(24)

c     Determine frequency content information

      omc=1.d0/(2.d0*dt)              !Nyquist frequency
      tott=dble(nt)*(dt)            !time series record length
      delom=1.d0/tott                !frequency step length

      if(iprompt.eq.0)then
         mnt=nt
         scalea=3.0
         iedge=0
         ix3=25
         ntracw=nx2-2*iedge
         do i=1,nx2-2*iedge
            ix2pt(i)=i+iedge
            ix3pt(i)=ix3
         enddo
         ifilt=0
         ilow=1
         ihigh=nom
         call taper(nom2,iomc1,ilow,ihigh,iomc4,ifilttype,rfilt)
         write(*,*)dble(ilow-1)*delom,dble(ihigh-1)*delom,' Hz.'
         xshift=0.d0
c         xshift=(tnd(1,1)-tbg(1,1))*xshift/100.d0
         iorder=0
         istart=1
         iend=ntracw
         iinc=1
         iskip=1
      elseif(iprompt.eq.1)then
c         write(*,*)'Select number of time samples to be viewed'
c         write(*,*)'(must be less than or equal to:',nt,').'
c         write(*,*)
c         read(*,*)mnt
c         write(*,*)
         iskip=1    !Doesn't work as setup . due to FFT 2n requirement
         mnt=nt
         write(*,*)'Select amplitude scaling for traces.'
         read(*,*)scalea
c     Read in the user selected grid points to be viewed.
         write(*,*)'Would you like to view a row of grid points?'
         write(*,*)'(yes=y,no=n)'
         read(*,*)answer2
         if(answer2.eq.'y'.or.answer2.eq.'Y')then
            write(*,*)'To view x2 from 1 to nx2 enter 2.'
            write(*,*)'Or to view x3 from 1 to nx3 enter 3.'
            read(*,*)ians2
            write(*,*)'Avoid plotting outer traces - enter margin value'
            write(*,*)'(for all traces enter 0).'
            read(*,*)iedge
            if(ians2.eq.2)then
               ntracw=nx2-2*iedge
               write(*,*)'Enter x3 coordinate.'
               read(*,*)ix3
               do i=1,nx2-2*iedge
                  ix2pt(i)=i+iedge
                  ix3pt(i)=ix3
               enddo
            elseif(ians2.eq.3)then
               ntracw=nx3-2*iedge
               write(*,*)'Enter x2 coordinate.'
               read(*,*)ix2
               do i=1,nx3-2*iedge
                  ix3pt(i)=i+iedge
                  ix2pt(i)=ix2
               enddo
            else
               write(*,*)
               write(*,*)'Error on input.'
               stop
            endif
         elseif(answer2.eq.'n'.or.answer2.eq.'N')then
            write(*,*)'Input number grid point waveforms to be viewed.'
            write(*,*)'(maximum is ',ntracwmx,' .)'
            read(*,*)ntracw
            write(*,*)'Grid waveforms to be viewed are: ',ntracw,' .'
            write(*,*)'Enter the specific grid points to be viewed.'
            write(*,*)'(ix2,ix3 enter).'
            do 70 i=1,ntracw
               if(i.lt.ntracw)write(*,*)'Point: ',i
               if(i.eq.ntracw)write(*,*)'Last point: ',i
               read(*,*)ix2pt(i),ix3pt(i)
               write(*,*)'Point: ',i,' = (',ix2pt(i),',',ix3pt(i),' )'
 70         continue
         else
            write(*,*)
            write(*,*)'Not a valid input.  Retry.'
            write(*,*)
         endif 

 222     write(*,*)'Do you want to filter the data?'
         read(*,*)answer3

         if(answer3.eq.'y'.or.answer3.eq.'Y')then
            ifilt=1
         elseif(answer3.eq.'n'.or.answer3.eq.'N')then
            ifilt=0
            ilow=1
            ihigh=nom
         else
            write(*,*)
            write(*,*)'Not a valid input.  Retry.'
            write(*,*)
            goto 222 
         endif

         imargin=3

         if(ifilt.eq.1)then
            write(*,*)'The number of frequency samples is:',int(rt),' .'
            write(*,*)'The frequency increment is:',delom,' Hz.'
            write(*,*)'Nyquist frequency (maximum) is:',omc,' Hz.'
            write(*,*)
 333        write(*,*)'Type of filtering: '
            write(*,*)'(1-low-pass,2-bandpass,3-highpass,4-notch)?'
            read(*,*)ifilttype
            if(ifilttype.eq.1)then
               ilow=1
               iomc1=ilow
               write(*,*)'Enter integer value of highest frequency.'
               write(*,*)'Note: must between',imargin,' and ',nom2,' .'
               read(*,*)ihigh
               iomc4=ihigh+imargin
            elseif(ifilttype.eq.2)then
               write(*,*)'Enter integer value of lowest and highest'
               write(*,*)'frequency: between',imargin,' and ',nom2,' .'
               read(*,*)ilow,ihigh
               iomc1=ilow-imargin
               iomc4=ihigh+imargin
            elseif(ifilttype.eq.3)then
               ihigh=nom2
               iomc4=ihigh
               write(*,*)'Enter integer value of lowest frequency.'
               write(*,*)'Note: must between',imargin,' and ',nom2,' .'
               read(*,*)ilow
               iomc1=ilow-imargin
            elseif(ifilttype.eq.4)then
               iomc1=1
               iomc4=nom2
               write(*,*)'Enter integer val of lowest and highest notch'
               write(*,*)'in frequency.  Note: must be between',imargin
               write(*,*)' and',nom2,' .'
               read(*,*)ilow,ihigh
            else
               write(*,*)
               write(*,*)'Invalid input.  Retry.'
               write(*,*)
               goto 333
            endif
         endif

         call taper(nom2,iomc1,ilow,ihigh,iomc4,ifilttype,rfilt)

         write(*,*)'Frequency band is (low,high):'
         write(*,*)dble(ilow-1)*delom,dble(ihigh-1)*delom,' Hz.'

         write(*,*)'Enter phase shift (in %) for waveform display'
         write(*,*)'(right-left: + and left-right: -1).'
         read(*,*)xshift
         xshift=(tnd(1,1)-tbg(1,1))*xshift/100.0

         write(*,*)'Choose trace ordering (0: normal and 1: reverse)'
         read(*,*)iorder
         if(iorder.eq.0)then
            istart=1
            iend=ntracw
            iinc=iskip
         elseif(iorder.eq.1)then
            istart=ntracw
            iend=1
            iinc=-iskip
         endif
      endif

      delom=2.d0*pi*delom

c     Selection of the wanted wavefields.

      ntrac=1
      do 90 i=istart,iend,iinc
         ix2=ix2pt(i)
         ix3=ix3pt(i)
         tbeg(ntrac)=sngl(tbg(ix2,ix3))
         tend(ntrac)=sngl(tnd(ix2,ix3))
         om=-delom
         do iom=1,nom
            om=om+delom
            cxshift=-ai*om*xshift
            if(ifilt.eq.0)then
               rfilt(iom)=1.d0
            endif
            if(iom.le.nom2)then
               work1(iom)=uvc(1,iom,ix3,ix2)*rfilt(iom)*cdexp(cxshift)
               work2(iom)=uvc(2,iom,ix3,ix2)*rfilt(iom)*cdexp(cxshift)
               work3(iom)=uvc(3,iom,ix3,ix2)*rfilt(iom)*cdexp(cxshift)
            else
               work1(iom)=dcmplx(0.d0,0.d0)
               work2(iom)=dcmplx(0.d0,0.d0)
               work3(iom)=dcmplx(0.d0,0.d0)
            endif
            if(ihidden.eq.1.and.iom.eq.2)then
               work1(iom)=dcmplx(0.d0,0.d0)
               work2(iom)=dcmplx(0.d0,0.d0)
               work3(iom)=dcmplx(0.d0,0.d0)
            endif
            if(iom.gt.1.and.iom.lt.nom)then
               work1(nt+2-iom)=dconjg(work1(iom))
               work2(nt+2-iom)=dconjg(work2(iom))
               work3(nt+2-iom)=dconjg(work3(iom))
            endif
         enddo

c     Call to inverse fft to bring data into time-domain

         isign=-1
         call FOUR1(work1,nt,isign)
         call FOUR1(work2,nt,isign)
         call FOUR1(work3,nt,isign)

         do 100 it=1,nt

            uvec(1,it,ntrac)=realpart(work1(it))/dble(nt)
            uvec(2,it,ntrac)=realpart(work2(it))/dble(nt)
            uvec(3,it,ntrac)=realpart(work3(it))/dble(nt)

 100     continue 

         if(ntrac.lt.ntracw)ntrac=ntrac+1

 90   continue

c      open(10,file='./Output/wavefld.time',form='formatted')
c      write(10,*)'nt,dt:',nt,dt
c      do it=1,nt
c         write(10,*)uvec(1,it,1),uvec(2,it,1),uvec(3,it,1)
c      enddo
c      close(10)

c     Find scaling for plotting from 1st level amplitudes

      if(ntrac.eq.ntracw)then 

         smaxx=0.0
         smaxy=0.0
         smaxz=0.0
         smaxa=0.0
            
c         do 2001 irec=1,ntracw
         do 2001 irec=istart,iend,iinc
               
            smaxx0=0.0
            smaxy0=0.0
            smaxz0=0.0
               
            do 2002 ipts=1,mnt
               smaxx0=amax1(smaxx0,sngl(dabs(uvec(1,ipts,irec))))
               smaxy0=amax1(smaxy0,sngl(dabs(uvec(2,ipts,irec))))
               smaxz0=amax1(smaxz0,sngl(dabs(uvec(3,ipts,irec))))
 2002       continue

c            write(*,*)
c            write(*,*)'irec=',irec,' smaxx,y,z: ',smaxx0,smaxy0,
c     +           smaxz0

            smaxx=amax1(smaxx,smaxx0)
            smaxy=amax1(smaxy,smaxy0)
            smaxz=amax1(smaxz,smaxz0)
            smaxa=amax1(smaxa,smaxx)
            smaxa=amax1(smaxa,smaxy)
            smaxa=amax1(smaxa,smaxz)
               
 2001    continue

         write(*,*)
         write(*,*)' smaxx,y,z, smaxa: ',smaxx,smaxy,smaxz,smaxa
         write(*,*)
            
      endif

c***********************************************************************
c     NCAR/GKS Plotting
c***********************************************************************

c     Open gks

      call gksinit(unit1,unit2,0)

      call dfclrs(2)
      call dfclrs(3)

      call assorted

      call bndary2

      call gstxci(0)

c     Create segment (x_1 direction)

      nseg=1
      call gcrsg(nseg)
      
      call gstxal(1,5)
      ctext='1-displacement'
      call mytext1(ctext,0.03,0.800,0.017)
      ctext='2-displacement'
      call mytext1(ctext,0.36,0.800,0.017)
      ctext='3-displacement'
      call mytext1(ctext,0.69,0.800,0.017)
      call gstxal(0,0)

c     Set window and viewport

      chfmx='(f7.2)'
      chfmy='(e9.2)'
c      xtext='ms'
      xtext='s'
      ytext='receiver'

c     Using start/stop times of final trace as world frame

c      scalex=1000.0
      scalex=1.0
      wxmin=tbeg(ntracw/2)*scalex
      wxmax=tend(ntracw)*scalex

      write(*,*)'wxmin,wxmax are ',wxmin,wxmax
      write(*,*)

      wymin=0.0 
      wymax=real(ntracw+1)

      dt=4.0
      wdx=(wxmax-wxmin)/dt
      wdy=1.0
      xmin=0.03
      xmax=0.33
      ymin=0.10
      ymax=0.98

      call gswn(1,wxmin,wxmax,wymin,wymax)
      call gsvp(1,xmin,xmax,ymin,ymax)

      ianotx=1
      ianoty=0
      call mytext2(xtext,ytext,chfmx,chfmy,ianotx,ianoty,nseg)

c      do 199 irec=1,ntracw,iskip
      do 199 irec=istart,iend,iinc

         do 200 ipts=1,mnt
            xfun(ipts)=(tbeg(irec)
     +           + real(ipts-1)*(tend(irec)-tbeg(irec))
     +           /real(mnt-1))*scalex
            yfun(ipts)=scalea*sngl(uvec(1,ipts,irec))/smaxa
     +           +wdy*real(irec)
 200     continue

         call gslwsc(2.5)
         call gsln(1)
         if(ix1.eq.2)call gsln(2)
         if(ix1.ge.3)call gsln(3)
         call gpl(mnt,xfun,yfun)
         call gsln(1)
         call gslwsc(1.0)

 199  continue

      call gclsg

c     Create segment (x_2 direction)

      nseg=nseg+1
      call gcrsg(nseg)

      xmin=0.36
      xmax=0.66
      ymin=0.10
      ymax=0.98

      call gswn(1,wxmin,wxmax,wymin,wymax)
      call gsvp(1,xmin,xmax,ymin,ymax)

      ianotx=1
      ianoty=0
      call mytext2(xtext,ytext,chfmx,chfmy,ianotx,ianoty,nseg)

c      do 299 irec=1,ntracw,iskip
      do 299 irec=istart,iend,iinc

         do 300 ipts=1,mnt
            xfun(ipts)=(tbeg(irec)
     +           + real(ipts-1)*(tend(irec)-tbeg(irec))
     +           /real(mnt-1))*scalex
            yfun(ipts)=scalea*sngl(uvec(2,ipts,irec))/smaxa
     +           +wdy*real(irec)
 300     continue

         call gslwsc(2.5)
         call gsln(1)
         if(ix1.eq.2)call gsln(2)
         if(ix1.ge.3)call gsln(3)
         call gpl(mnt,xfun,yfun)
         call gsln(1)
         call gslwsc(1.0)

 299  continue

      call gclsg

c     Create segment (x_3 direction)

      nseg=nseg+1
      call gcrsg(nseg)

      xmin=0.69
      xmax=0.99
      ymin=0.10
      ymax=0.98

      call gswn(1,wxmin,wxmax,wymin,wymax)
      call gsvp(1,xmin,xmax,ymin,ymax)

      ianotx=1
      ianoty=0
      call mytext2(xtext,ytext,chfmx,chfmy,ianotx,ianoty,nseg)

c      do 399 irec=1,ntracw,iskip
      do 399 irec=istart,iend,iinc

         do 400 ipts=1,mnt
            xfun(ipts)=(tbeg(irec)
     +           + real(ipts-1)*(tend(irec)-tbeg(irec))
     +           /real(mnt-1))*scalex
            yfun(ipts)=scalea*sngl(uvec(3,ipts,irec))/smaxa
     +           +wdy*real(irec)
 400     continue

         call gslwsc(2.5)
         call gsln(1)
         if(ix1.eq.2)call gsln(2)
         if(ix1.ge.3)call gsln(3)
         call gpl(mnt,xfun,yfun)
         call gsln(1)
         call gslwsc(1.0)

 399  continue

      call gclsg

      write(*,*)'Finished plotting.'
      write(*,*)' '
      call frame
      call grclse2(ansm,unit1,unit2)
      close(7)

      stop
      end

c***********************************************************************
      SUBROUTINE FOUR1(DATA,NN,ISIGN)
c     Replaces DATA by its discrete Fourier transform, if ISIGN is
c     input as 1; or replaces DATA by NN times its inverse discrete
c     Fourier transform, if ISIGN is input as -1.  DATA is a complex
c     array of length NN or, equivalently, a real array of length 2*NN.
c     NN must be an integer power of 2 (this is not checked for!).
c***********************************************************************
 
      implicit real*8(a-h,o-z)

      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
      DIMENSION DATA(*)
      N=2*NN
      J=1
      DO 11 I=1,N,2
        IF(J.GT.I)THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        ENDIF
        M=N/2
1       IF ((M.GE.2).AND.(J.GT.M)) THEN
          J=J-M
          M=M/2
        GO TO 1
        ENDIF
        J=J+M
11    CONTINUE
      MMAX=2
2     IF (N.GT.MMAX) THEN
        ISTEP=2*MMAX
        THETA=6.28318530717959/(ISIGN*MMAX)
        WPR=-2.0*SIN(0.5*THETA)**2
        WPI=SIN(THETA)
        WR=1.0
        WI=0.0
        DO 13 M=1,MMAX,2
          DO 12 I=M,N,ISTEP
            J=I+MMAX
            TEMPR=WR*DATA(J)-WI*DATA(J+1)
            TEMPI=WR*DATA(J+1)+WI*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
12        CONTINUE
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
13      CONTINUE
        MMAX=ISTEP
      GO TO 2
      ENDIF
      RETURN
      END

c***********************************************************************
      subroutine taper(nom2,iomc1,iomc2,iomc3,iomc4,ifilt,rfilt)
c     This subroutine applies a high pass filter to the frequency data
c***********************************************************************

      implicit real*8(a-h,o-z)

      dimension rfilt(nom2)

      pi=4.0*atan(1.0)
  
      if(ifilt.eq.1)then                   !Low-pass filter

         rat1=real(iomc4-iomc3)/pi

         do iom=1,nom2

            if(iomc4.eq.iomc3)then

               rfilt(iom)=1.0

            else

               if(iom.lt.iomc3)then
                  rfilt(iom)=1.0
               elseif(iom.ge.iomc3.and.iom.le.iomc4)then
                  rat=real(iom-iomc3)/rat1
                  rfilt(iom)=1.0-(1.0-cos(rat))/2.0
               else
                  rfilt(iom)=0.0
               endif

            endif

         enddo

      elseif(ifilt.eq.2)then               !Band-pass filter

         rat1=real(iomc2-iomc1)/pi
         rat2=real(iomc4-iomc3)/pi

         do iom=1,nom2

            if(iomc1.eq.iomc2)then
               rfilt(iom)=1.0
            elseif(iomc3.eq.iomc4)then
               rfilt(iom)=1.0
            else

               if(iom.lt.iomc1)then
                  rfilt(iom)=0.0
               elseif(iom.ge.iomc1.and.iom.le.iomc2)then
                  rat=real(iom-iomc1)/rat1
                  rfilt(iom)=(1.0-cos(rat))/2.0
               elseif(iom.gt.iomc2.and.iom.lt.iomc3)then
                  rfilt(iom)=1.0
               elseif(iom.ge.iomc3.and.iom.le.iomc4)then
                  rat=real(iom-iomc3)/rat1
                  rfilt(iom)=1.0-(1.0-cos(rat))/2.0
               else
                  rfilt(iom)=0.0
               endif

            endif

         enddo

      elseif(ifilt.eq.3)then               !High-pass filter

         rat1=real(iomc2-iomc1)/pi

         do iom=1,nom2

            if(iomc2.eq.iomc1)then

               rfilt(iom)=1.0

            else

               if(iom.lt.iomc1)then
                  rfilt(iom)=0.0
               elseif(iom.ge.iomc1.and.iom.le.iomc2)then
                  rat=real(iom-iomc2)/rat1
                  rfilt(iom)=1.0-(1.0-cos(rat))/2.0
               else
                  rfilt(iom)=1.0
               endif

            endif

         enddo

      elseif(ifilt.eq.4)then

         do iom=1,nom2

            if(iomc1.eq.iomc2)then
               rfilt(iom)=1.0
            elseif(iomc3.eq.iomc4)then
               rfilt(iom)=1.0
            else

               if(iom.lt.iomc1)then
                  rfilt(iom)=1.0
               elseif(iom.ge.iomc1.and.iom.le.iomc2)then
                  rat=1.0-real(iom-iomc1)/rat1
                  rfilt(iom)=(1.0-cos(rat))/2.0
               elseif(iom.gt.iomc2.and.iom.lt.iomc3)then
                  rfilt(iom)=0.0
               elseif(iom.ge.iomc3.and.iom.le.iomc4)then
                  rat=1.0-real(iom-iomc3)/rat1
                  rfilt(iom)=1.0-(1.0-cos(rat))/2.0
               else
                  rfilt(iom)=1.0
               endif

            endif

         enddo

      endif

      return
      end

c***********************************************************************
      SUBROUTINE GKSINIT(UNIT1,UNIT2,ifile)
c***********************************************************************

      INTEGER UNIT1,UNIT2,ifile

      OPEN(UNIT1,FILE='./Log/ERRORS')
      OPEN(UNIT2,FILE='./Log/METAFILE')

      CALL GOPKS(UNIT1,1024)

C     Open and activate workstation for PS 

c      if(ifile.eq.0) CALL NGSETC('ME','./Output/plot_3C.ps')
c      CALL GOPWK(2,2,20)
      if(ifile.eq.0) CALL NGSETC('ME','./Output/plot_3C.eps')
      CALL GOPWK(2,2,21)
      CALL GACWK(2)

C     Open and activate workstation for X11

      CALL GOPWK(3,0,8)
      CALL GACWK(3)
      CALL NGSETI('SC',2)

C     Define colous for X11 on a continuous spectrum.

      NITER=201
      DO 11 K = 1,NITER 
         H = REAL(K)/REAL(NITER)*360.
C         CALL HLSRGB(H,50.,100.,RV,GV,BV)
         CALL HLSRGB(H,50.,100.,RV,GV,BV)
         CALL GSCR(2,K,RV,GV,BV)
 11   CONTINUE

      CALL GSCR(2,NITER+1,1.,0.,0.)

C      CALL GSCR(2,2,1.,0.,0.)
C      CALL GSCR(2,3,0.,1.,0.)
C      CALL GSCR(2,4,0.,0.,1.)
C      CALL GSCR(2,5,1.,0.5,0.)
C      CALL GSCR(2,6,1.,0.,0.5)
C      CALL GSCR(2,7,0.,0.5,0.2)
C      CALL GSCR(2,8,0.5,0.5,0.5)
C      CALL GSCR(2,9,0.3,0.5,0.1)

C     Open and activate workstation for WISS

      CALL GOPWK(4,UNIT2,3)
      CALL GACWK(4)

      RETURN
      END

c***********************************************************************
      SUBROUTINE DFCLRS(IWKID)
c***********************************************************************

C     Define a set of RGB color triples for colors 1 through 15.

      DIMENSION RGBV(3,15)

C     Define the RGB color triples needed below.

      DATA RGBV / 1.00 , 1.00 , 1.00 ,
     +     0.70 , 0.70 , 0.70 ,
     +     0.75 , 0.50 , 1.00 ,
     +     0.50 , 0.00 , 1.00 ,
     +     0.00 , 0.00 , 1.00 ,
     +     0.00 , 0.50 , 1.00 ,
     +     0.00 , 1.00 , 1.00 ,
     +     0.00 , 1.00 , 0.60 ,
     +     0.00 , 1.00 , 0.00 ,
     +     0.70 , 1.00 , 0.00 ,
     +     1.00 , 1.00 , 0.00 ,
     +     1.00 , 0.75 , 0.00 ,
     +     1.00 , 0.38 , 0.38 ,
     +     1.00 , 0.00 , 0.38 ,
     +     1.00 , 0.00 , 0.00 /
     
C     Define 16 different color indices, for indices 0 through 15.  The
C     color corresponding to index 0 is black and the color 
C     corresponding to index 1 is white.

      CALL GSCR (IWKID,0,0.,0.,0.)

      DO 101 I=1,15
         CALL GSCR (IWKID,I,RGBV(1,I),RGBV(2,I),RGBV(3,I))
 101  CONTINUE

C     Done.

      RETURN
      END

c***********************************************************************
      SUBROUTINE ASSORTED
c***********************************************************************

      INTEGER GSTRP,GSTRKP,LASF(13),GBUNDL,GINDIV
      DATA GBUNDL,GINDIV/0,1/,GSTRP,GSTRKP/0,2/

c     Ideas taken from example gstxfp (text font and precision)   
c     define own style with bundle to set style for each specific  
c     workstation (set asf to bundel)

c     linetype

      LASF(1)=GINDIV

c     linewidth

      LASF(2)=GINDIV

c     polyline colour index

      LASF(3)=GINDIV

c     marker type

      LASF(4)=GINDIV

c     marker size

      LASF(5)=GINDIV

c     polymarker colour index

      LASF(6)=GINDIV
         
c     text font and precision

      LASF(7)=GINDIV

c     character expansion factor

      LASF(8)=GINDIV

c     character spacing

      LASF(9)=GINDIV

c     text colour index

      LASF(10)=GINDIV

c     fill area interior style (set on bundl for fill)

      LASF(11)=GINDIV

c     fill area style (set on bundl for fill)

      LASF(12)=GINDIV

c     fill area colour (set on bundl for fill)

      LASF(13)=GINDIV

      CALL GSASF(LASF)

C     DEFINE FILL AREA COLOUR INDICES

c      DO 1010 ICOL=1,15

c     Workstation, fill area index, inter style, style index, col ind

c         CALL GSFAR(1,ICOL,1,1,ICOL)

c     next version works on geola, because of black background

c         CALL GSFAR(1,ICOL,1,1,16-ICOL)

c         CALL GSFAR(2,ICOL,1,1,1)

c 1010 CONTINUE

C     Workstation-text index-font-precision-expansion-space-colour

c      CALL GSTXR(4,1,-5,GSTRKP,1.,0.,1)

c      CALL GSTXR(1,1,4,GSTRP,1.,0.,15)

        RETURN
        END

c***********************************************************************
      SUBROUTINE BNDARY2
c     Routine to draw the plotter frame and fill the background
c***********************************************************************

      REAL PX(5),PY(5)

C     SET FILL AREA STYLE TO SOLID (GSS*GKS MANUAL V.1 P. 7-14
C     AND V.2, FORTRAN SECTION, P. 13)

      CALL GSFAIS(1)

C     SET THE COLOUR OF THE FILL TO LIGHT GREY (FOR IBMVGA12 DRIVER)
C     CALL GSFACI(9)
C     OR CYAN

      CALL GSFACI(1)

C     SET CORNERS OF THE BOX TO FILL IN WORLD COORDINATES
C     (JUST GUESSING THIS IS THE SIZE OF THE ENTIRE WORLD SPACE)

      PX(1)=0.0
      PX(2)=1.
      PX(3)=PX(2)
      PX(4)=PX(1)
      PX(5)=PX(1)
      PY(1)=0.0
      PY(2)=0.0
      PY(3)=1.
      PY(4)=PY(3)
      PY(5)=PY(1)

C     FILL IT

      CALL GFA(5,PX,PY)

C     PLOT PERIMETER OF BOX IN ANOTHER COLOUR
C     (RED ON IBMVGA12)

      CALL GSLWSC(3.)
      CALL GSPLCI(0)

c     CALL GPL(5,PX,PY)

      RETURN
      END

c***********************************************************************
      subroutine mytext1(ctext,xpos,ypos,xtexh)
c***********************************************************************

      character*80 ctext
      common/vcord/xmax,xmin,ymax,ymin

      call gselnt(0)

c     Plot text

c     Set character height for text

      call gschh(xtexh)

c     Set character font

      call gstxfp(-13,2)

c     Set character up vector

      call gschup(-1.,0.)

c     Set text alignment ... before calling...
c      call gstxal(2,5)

c     Find string length

      call textl(ctext,lenb,length)

      call gtx(xpos,ypos,ctext(1:length))

      call gselnt(1)

      return
      end

c***********************************************************************
      subroutine textl(string,lenb,lene)
c***********************************************************************

      character*(*) string, subst*1

      do 19 i=len(string),1,-1
         subst=string(i:i)
         if (subst.ne.' ') goto 102
 19   continue
 102  lene=i

c zero length string gets length one

      if (i.eq.1) goto 103

      do 21 i=1,len(string),1
         subst=string(i:i)
         if (subst.ne.' ') goto 103
 21   continue
 103  lenb=i
      
      return
      end

c***********************************************************************
      subroutine mytext2(xtext,ytext,chfm1,chfm2,ianotx,ianoty,nseg)
c***********************************************************************

      integer inxt,inyt
      real xmin,xmax,ymin,ymax
      real ndx,ndy,nxt,nyt,tick

      dimension yf(5),xf(5),xt(2),yt(2),ytt(2),ytb(2),xtr(2),xtl(2)
      
      character*80 xtext, ytext
      character chnum*9,numdat*9,chfm1*6,chfm2*6,chfmi1*4,
     +     chfmi2*4,chfm22*8

      common/wcord/wxmax,wxmin,wymax,wymin,wdx,wdy
      common/vcord/xmax,xmin,ymax,ymin

c     Plot with repect to viewport

      call gselnt(0)

c     Set parameter

c     Set ticklength

      tick=0.008

c     Set character font

      call gstxfp(-13,2)

c     Set character spacing

      call gschsp(0.1)

c     Set character height for labels

      xlabh=0.015

c     Set character height for text

      xtexh=xlabh*1.25

c     Set line width for frame

      call gslwsc(1.5)

c     Simple algorithm to determine formats for labels
c     maximum length of a numerical label is 9 characters
c     (seems adequate) f or e or i format should work

      if (chfm2(3:3).eq.'0') then

c     No annotation on y axis

         jleny=0
         chfmi1='(a'//chfm1(3:3)//')'
         write (numdat(1:1),'(a1)') chfm1(3:3)
         read (numdat(1:1),'(i1)') jlenx
      
      else
      
         chfmi1='(a'//chfm1(3:3)//')'
         write (numdat(1:1),'(a1)') chfm1(3:3)
         read (numdat(1:1),'(i1)') jlenx
         chfmi2='(a'//chfm2(3:3)//')'
         write (numdat(1:1),'(a1)') chfm2(3:3)
         read (numdat(1:1),'(i1)') jleny
      
      endif
      
      if (chfm2(2:2).eq.'e') then
         chfm22='(1p'//chfm2(2:6)
      else
         chfm22=chfm2
      endif
      
c     End label setup

      nxt = (wxmax - wxmin)/wdx
      nyt = (wymax - wymin)/wdy
      ndx = (xmax - xmin)/nxt
      ndy = (ymax - ymin)/nyt

c     Plot frame in viewport coordinates

      xf(1)= xmin
      xf(2)= xmax
      xf(3)= xmax
      xf(4)= xmin
      xf(5)= xf(1)
      yf(1)= ymin
      yf(2)= ymin
      yf(3)= ymax
      yf(4)= ymax
      yf(5)= yf(1)

c     Polyline frame

      call gpl(5,xf,yf)

      call gschh(xlabh)

c     Plot y-axis

      inyt=int(nyt)
      do i=0,inyt

         yt(1)=float(i)*ndy+ymin
         yt(2)=yt(1)
         xtl(1)=xmin
         xtl(2)=xmin+tick
         xtr(1)=xmax
         xtr(2)=xmax-tick
         xtx=xmin-2.*tick

         anum = real(i)*wdy+wymin

         if (jleny.ne.0) then
            write (numdat(1:9),chfm22) anum
            read (numdat(1:9),chfmi2) chnum
            
            call gpl(2,xtl,yt)
            call gpl(2,xtr,yt)

c     Set text alignment

            call gstxal(3,3)
            if(ianoty.ne.0)call gtx(xtx,yt(1),chnum(1:jleny))
         endif
      enddo

c     Plot x-axis

      inxt=int(nxt)
      do i=0,inxt
         xt(1)=float(i)*ndx+xmin
         xt(2)=xt(1)
         ytt(1)=ymax
         ytt(2)=ymax-tick
         ytb(1)=ymin
         ytb(2)=ymin+tick
         ytx = ymin - 2.*tick
         anum = real(i)*wdx + wxmin

         if (chfm1(2:2).ne.'i') then
            write (numdat(1:9),chfm1) anum
         else
            write (numdat(1:9),chfm1) int(anum)
         endif
         read (numdat(1:9),chfmi1) chnum

         call gpl(2,xt,ytt)
         call gpl(2,xt,ytb)

c     Set text alignment
         call gstxal(3,3)
         call gschup(-1.,0.)
         if(nseg.eq.2)then
            if(ianotx.ne.0)call gtx(xt(1),ytx,chnum(1:jlenx))
         endif
      enddo
c     Write axis text

c     Set character height for axis text

      call gschh(xtexh)

c     Y-axis

c     Set character up vector

      call gschup(-1.,0.)

c     Set text alignment

      call gstxal(2,5)

c     Find string length

      call textl(ytext,lenb,length)

      xpos=xmin-(float(jleny)+1.)*xlabh
      ypos=ymin+(ymax-ymin)/2.
      if(ianoty.ne.0)call gtx(xpos,ypos,ytext(1:length))

c     X-axis

c     Set character up vector

      call gschup(0.,1.)

c     Set text alignment

      call gstxal(1,1)

c     Find string length

      call textl(xtext,lenb,length)

      if(nseg.eq.3)then
c      xpos=xmin+(xmax-xmin)/2.
c      ypos=ymin-3.0*xlabh
         xpos=xmin
         ypos=ymin-xlabh*.9
         if(ianotx.ne.0)call gtx(xpos,ypos,xtext(1:length))
      endif

c     Set text alignment

      call gstxal(0,0)

      call gselnt(1)

      return
      end

c***********************************************************************
      SUBROUTINE GRCLSE2(ANSM,UNIT1,UNIT2)
c***********************************************************************

      CHARACTER*1 ANSM
      INTEGER UNIT1,UNIT2

      CALL GDAWK(2)
      CALL GCLWK(2)
      CALL GDAWK(3)
      CALL GCLWK(3)
      CALL GDAWK(4)
      CALL GCLWK(4)
C      CALL GDAWK(0)
C      CALL GCLWK(0)
      CALL GCLKS

C      IF(ANSM.EQ.'Y'.OR.ANSM.EQ.'y')THEN
C         CALL GERCLS(UNIT1)
C      ELSEIF(ANSM.EQ.'N'.OR.ANSM.EQ.'n')THEN
          CLOSE(UNIT1)
          CLOSE(UNIT2)
C      ENDIF

      RETURN
      END

c***********************************************************************
      subroutine rotvec(a,b)
c     This subroutine performs a transformation of the vector from
c     one coordinate system to another defined by the transformation
c     matrix Aij.  Aij consists of the local coordinate frame written
c     in terms of the global coordinates (each row consists of a
c     basis vector - first row is wavefront normal .. normalized 
c     slowness).  Note: transformation is from global to local frame.
c***********************************************************************

      implicit real*4 (a-h,o-z)
      dimension a(3,3),b(3),c(3)

c     A(X) --> a(x)  (global to local frame)

      do i=1,3
         c(i)=0.0
         do j=1,3
            c(i)=c(i)+a(i,j)*b(j)
         enddo
      enddo

      do i=1,3
         b(i)=c(i)
      enddo

      return
      end

c***********************************************************************
      subroutine rinv3(a,ainv)
c     The subroutine inverts the 3x3 matrix a.
c***********************************************************************

      implicit real*4 (a-h,o-z)
      dimension a(3,3),ainv(3,3),co(3,3)

      co(1,1)=(a(2,2)*a(3,3)-a(2,3)*a(3,2))
      co(1,2)=-(a(2,1)*a(3,3)-a(2,3)*a(3,1))
      co(1,3)=(a(2,1)*a(3,2)-a(2,2)*a(3,1))

      co(2,1)=-(a(1,2)*a(3,3)-a(1,3)*a(3,2))
      co(2,2)=(a(1,1)*a(3,3)-a(1,3)*a(3,1))
      co(2,3)=-(a(1,1)*a(3,2)-a(1,2)*a(3,1))

      co(3,1)=(a(1,2)*a(2,3)-a(1,3)*a(2,2))
      co(3,2)=-(a(1,1)*a(2,3)-a(1,3)*a(2,1))
      co(3,3)=(a(1,1)*a(2,2)-a(1,2)*a(2,1))

      deta=a(1,1)*co(1,1)+a(1,2)*co(1,2)+a(1,3)*co(1,3)
      detb=a(2,1)*co(2,1)+a(2,2)*co(2,2)+a(2,3)*co(2,3)
      detc=a(3,1)*co(3,1)+a(3,2)*co(3,2)+a(3,3)*co(3,3)

c     Check to see if a is invertable

      if(deta.eq.0.0.or.detb.eq.0.0.or.detc.eq.0.0)then
         write(*,*)'Error in rinv3.  Determinant equal to zero.'
         write(*,*)'Program stopped.'
         stop
      endif

      do 10 i=1,3
        do 11 j=1,3
          ainv(i,j)=co(j,i)/deta
11      continue
10    continue

      return
      end
