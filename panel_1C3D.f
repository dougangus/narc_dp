c........+.........+.........+.........+.........+.........+.........+..
c  Program Panel_1C3D
c
c  This program plots waveform results in a 3D array or panel.  It 
c  contains the following 6 subroutines: Grinit, Grclse, Four1, Taper,
c  Dfclrs, and Bndary2.
c
c  Copyright (c) 2005 D. A. Angus & C. J. Thomson.
c  All rights reserved by the authors.
c  Last modified: December 5, 2006.
c***********************************************************************

c     Variable initialisation.

      implicit real*4 (a-h,o-z)

      parameter (ntmx=257,ntrmx=129,nx1mx=5)
      parameter (nx2mx=201,nx3mx=201)
      parameter (ntracwmx=201)

      integer unit1,unit2

      real*8 uvec(3,ntmx,ntrmx),rfilt(ntmx)
      real*8 dt,pi,tott,om,delom,dx1,x1o,omc,xshift
      real*8 tbg(nx2mx,nx3mx),tnd(nx2mx,nx3mx)

      real*4 xfun(ntmx),yfun1(ntmx)
      real*4 tbeg(ntrmx),tend(ntrmx)
      real*4 wxmax,wxmin,wymax,wymin,wdx,wdy
      real*4 xmax,xmin,ymax,ymin

      complex*16 uvc(3,ntmx,nx3mx,nx2mx),work1(ntmx),work2(ntmx),
     +     work3(ntmx),cxshift,ai

      character*80 xtext,ytext,ctext
      character*6 chfmx,chfmy
      character*1 ansm,answer1,answer2,answer3

      dimension nta(nx2mx,nx3mx),nta2(ntrmx)
      dimension ix2pt(ntracwmx),ix3pt(ntracwmx)
      dimension ix2pt0(ntracwmx), ix2pt1(ntracwmx), ix2pt2(ntracwmx),
     .          ix2pt3(ntracwmx), ix3pt0(ntracwmx), ix3pt1(ntracwmx), 
     .          ix3pt2(ntracwmx), ix3pt3(ntracwmx)
      dimension iplanes(100)

      common/wcord/wxmax,wxmin,wymax,wymin,wdx,wdy
      common/vcord/xmax,xmin,ymax,ymin

      data unit1/14/,unit2/15/

c     Srfplt variables

      parameter (nymx=500,ntot=ntmx*nymx,ntot2=2*ntot)

      COMMON /SRFIP1/ IFR,ISTP,IROTS,IDRX,IDRY,IDRZ,IUPPER,ISKIRT,
     +     NCLA,THETA,HSKIRT,CHI,CLO,CINC,ISPVAL
      COMMON /PZSINT/ IPZSMJ,IPZSMN,IPZSTX

      REAL XX(nymx),YY(ntmx),Z1(ntot),S(6)
      dimension mscrat(ntot2)

      DATA S(1),S(2),S(3),S(4),S(5),S(6)/
     +     12.0,-08.0,08.0,0.0,0.0,0.0/

c Specify coordinates for plot titles.  The values CX and CY
c define the center of the title string in a 0. to 1. range.

      DATA CX/0.5/,CY/0.72/

c     Open input/output files.

      open(unit=24,file='./Output/wavefld.out',form='unformatted')

      pi=4.d0*datan(1.d0)
      spi=sngl(pi)
      ai=dcmplx(0.d0,1.d0)

c     Decide which plane to view

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
      write(*,*)
      do im=1,ilast
         write(*,*)'Tag:',im,' ix1:',iplanes(im)
      enddo
      write(*,*)
      write(*,*)'Select three planes to view using integer tags.'
      write(*,*)
      read(*,*)iplane1,iplane2,iplane3

c     Read frequency-domain wavefield

      jplane=1
      kplane=iplane1

1112  continue

      do ip=jplane,kplane

         read(24)jx1,nx1,nx2,nx3
         read(24)nt,dt,nom,nom2

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

      write(*,*)'Finished reading wavefield frequency-domain data'
      write(*,*)'from the first plane'
      write(*,*)

      if(kplane.eq.iplane1)then

c     Determine frequency content information

      omc=1.d0/(2.d0*dt)            !Nyquist frequency
      tott=dreal(nt)*(dt)           !time series record length
      delom=1.d0/tott               !frequency step length

c     Read in the user selected grid points to be viewed.

111      write(*,*)
         write(*,*)'Primary receiver coordinate definition:'
         write(*,*)'To view x2 from 1 to nx2 enter 2.'
         write(*,*)'Or to view x3 from 1 to nx3 enter 3.'
         read(*,*)ians2
         write(*,*)
         write(*,*)'To avoid plotting outer traces enter margin value'
         write(*,*)'(for all traces enter 0).'
         read(*,*)iedge
         if(ians2.eq.2)then
            ntracw=nx2-2*iedge
            write(*,*)
            write(*,*)'Enter three desired x3 coordinates.'
            read(*,*)ix31,ix32,ix33
            do 50 i=1,nx2-2*iedge
               ix2pt0(i)=i+iedge
               ix3pt1(i)=ix31
               ix3pt2(i)=ix32
               ix3pt3(i)=ix33
 50         continue
         elseif(ians2.eq.3)then
            ntracw=nx3-2*iedge
            write(*,*)
            write(*,*)'Enter three desired x2 coordinates.'
            read(*,*)ix21,ix22,ix23
            do 60 i=1,nx3-2*iedge
               ix3pt0(i)=i+iedge
               ix2pt1(i)=ix21
               ix2pt2(i)=ix22
               ix2pt3(i)=ix23
 60         continue
         else
            write(*,*)
            write(*,*)'Error on input.'
            stop
         endif

      write(*,*)
 222  write(*,*)'Do you want to filter the data?'
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
         write(*,*)
         write(*,*)'The number of frequency samples is:',int(rt),' .'
         write(*,*)'The frequency increment is:',delom,' Hz.'
         write(*,*)'Nyquist frequency (maximum) is:',omc,' Hz.'
         write(*,*)
 333     write(*,*)'Type of filtering: '
         write(*,*)'(1-low-pass,2-bandpass,3-highpass,4-notch)?'
         read(*,*)ifilttype
         if(ifilttype.eq.1)then
            ilow=1
            iomc1=ilow
            write(*,*)
            write(*,*)'Enter integer value of highest frequency.'
            write(*,*)'Note: must be between',imargin,' and ',nom2,' .'
            read(*,*)ihigh
            iomc4=ihigh+imargin
         elseif(ifilttype.eq.2)then
            write(*,*)
            write(*,*)'Enter integer value of lowest and highest'
            write(*,*)'frequency: between',imargin,' and ',nom2,' .'
            read(*,*)ilow,ihigh
            iomc1=ilow-imargin
            iomc4=ihigh+imargin
         elseif(ifilttype.eq.3)then
            ihigh=nom2
            iomc4=ihigh
            write(*,*)
            write(*,*)'Enter integer value of lowest frequency.'
            write(*,*)'Note: must be between',imargin,' and ',nom2,' .'
            read(*,*)ilow
            iomc1=ilow-imargin
         elseif(ifilttype.eq.4)then
            iomc1=1
            iomc4=nom2
            write(*,*)
            write(*,*)'Enter integer value of lowest and highest notch'
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

      write(*,*)
      write(*,*)'Frequency band is (low,high):'
      write(*,*)dreal(ilow-1)*delom,dreal(ihigh-1)*delom,' Hz.'

      write(*,*)
      write(*,*)'Enter phase shift (in %) for waveform display .'
      read(*,*)xshift
      xshift=(tnd(1,1)-tbg(1,1))*xshift/100.d0

      write(*,*)
      write(*,*)'Select the displacement component to plot:'
      read(*,*)iwcomp
      write(*,*)'Plotting component ',iwcomp

      write(*,*)
      write(*,*)'Reverse the time axis?'
      write(*,*)'  no=0, time increases left to right'
      write(*,*)' yes=1, time increases right to left'
      read(*,*)itswitch
      write(*,*)

      delom=2.d0*pi*delom

      endif

c     Selection of the wanted wavefields.

c     Loop on the three `columns'

      do 1111 icol=1,3

         ntrac=1
         do 90 i=1,ntracw

            if(ians2.eq.2)then
               ix2=ix2pt0(i)
               if(icol.eq.1)then
                  ix3=ix3pt1(i)
               elseif(icol.eq.2)then
                  ix3=ix3pt2(i)
               elseif(icol.eq.3)then
                  ix3=ix3pt3(i)
               endif
            endif

            if(ians2.eq.3)then
               ix3=ix3pt0(i)
               if(icol.eq.1)then
                  ix2=ix2pt1(i)
               elseif(icol.eq.2)then
                  ix2=ix2pt2(i)
               elseif(icol.eq.3)then
                  ix2=ix2pt3(i)
               endif
            endif

            tbeg(ntrac)=tbg(ix2,ix3)
            tend(ntrac)=tnd(ix2,ix3)
            nta2(ntrac)=nta(ix2,ix3)

            om=-delom

            do iom=1,nom

               om=om+delom

               cxshift=-ai*om*xshift

               if(ifilt.eq.0)then
                  rfilt(iom)=1.d0
               endif

               if(iom.le.nom2)then
                  work1(iom)=
     +                 uvc(iwcomp,iom,ix3,ix2)*rfilt(iom)*cdexp(cxshift)
               else
                  work1(iom)=dcmplx(0.d0,0.d0)
               endif

               if(ihidden.eq.1.and.iom.eq.2)then
                  work1(iom)=dcmplx(0.d0,0.d0)
               endif

               if(iom.gt.1.and.iom.lt.nom)then
                  work1(nt+2-iom)=dconjg(work1(iom))
               endif

            enddo

c     Call to inverse fft to bring data into time-domain

            isign=-1
            call FOUR1(work1,nt,isign)

            do 303 it=1,nt

c     Note storing in uvec(1)

               if(itswitch.eq.0)then
                  uvec(1,it,ntrac)=realpart(work1(it))/dble(nt)
               elseif(itswitch.eq.1)then
                  uvec(1,it,ntrac)=realpart(work1(nt+1-it))/dble(nt)
               endif

 303        continue 

            if(ntrac.lt.ntracw)ntrac=ntrac+1

 90      continue
         
c     Find scaling for plotting from 1st level amplitudes

 699     if(ntrac.eq.ntracw)then 

            smaxx=0.0
            smaxa=0.0
            
            do 2001 irec=1,ntracw
               
               smaxx0=0.0
               
               do 2002 ipts=1,nt
                  smaxx0=amax1(smaxx0,sngl(dabs(uvec(1,ipts,irec))))
 2002          continue

c               write(*,*)
c               write(*,*)'irec=',irec,' smaxx: ',smaxx0

               smaxx=amax1(smaxx,smaxx0)
               smaxa=amax1(smaxa,smaxx)
               
 2001       continue

c            write(*,*)
c            write(*,*)' smaxx, smaxa: ',smaxx,smaxa
c            write(*,*)
            
         endif

         IFR=0
         IDRY=0
         IUPPER=1

c     SRFACE SECTION  

c         write(*,*)'S(1,2,3) =',S(1),S(2),S(3)

c     Read data

         zmax=0.0
         do 99 it=1,nt
            YY(it)=real(it)/real(nt)
            do 100 jt=1,ntracw
               XX(jt)=real(jt)/real(ntracw)
               Z1(jt+ntracw*(it-1))=uvec(1,it,jt)
               zmax=amax1(zmax,abs(Z1(jt+ntracw*(it-1))))
 100        continue
 99      continue

c         write(*,*)'zmax=',zmax

         zmax=zmax*8.0

         do 101 it=1,nt
            do 102 jt=1,ntracw
               Z1(jt+ntracw*(it-1))=Z1(jt+ntracw*(it-1))/zmax
 102        continue
 101     continue

c         write(*,*)'XX(ntracw)=',XX(ntracw)
c         write(*,*)'YY(nt)=',YY(nt)

         if(kplane.eq.iplane1.and.icol.eq.1)then

c     OPEN GKS FOR MACHINE TYPE

            ansm='N'
            CALL GRINIT(ANSM,UNIT1,UNIT2)

            call dfclrs(2)
            call dfclrs(3)

            call bndary2

            IERROR = 0

c     Select the normalization transformation 0.

            CALL GSELNT(0)

c     Frame 2 -- The SRFACE entry.

            CALL GSTXAL(2,3)

c     Set the character height.

            CALL GSCHH(.016)

c     Write the text.

c            CALL GTX(CX,CY,'X-COMPONENT')

         endif

c     NOTE: SRFACE2 IS DIFFERENT FROM ORIGINAL SRFACE ONLY IN THAT
c     IT DOES NOT ADVANCE THE FRAME INTERNALLY (FRAME2 CALLED,
c     NOT FRAME)
c     ALSO: THE CALL TO SET IN SRFACEs2 IS SUPPRESSED IN ORDER TO
c     CREATE AN ARRAY OF WAVEFORMS

         CALL GSLWSC(1.0)
         WOSW=0.5
         CALL GETSET (XVPL,XVPR,YVPB,YVPT,XWDL,XWDR,YWDB,YWDT,LNLG)

c         WRITE(*,*)'XVPL,XVPR,YVPB,YVPT=',XVPL,XVPR,YVPB,YVPT
c         WRITE(*,*)
c         WRITE(*,*)'XWDL,XWDR,YWDB,YWDT=',XWDL,XWDR,YWDB,YWDT
c         WRITE(*,*)

c     Following values ascertained from SRFACEs[2].

         XWDL=1.0
         XWDR=1024.0
         YWDB=1.0
         YWDT=1024.0
         LNLG=1

         if(icol.eq.1)then

            VPL1=0.04
            VPR1=0.44
            if(kplane.eq.iplane1)then
               VPB1=0.47
               VPT1=0.87
            elseif(kplane.eq.iplane2)then
               VPB1=0.23
               VPT1=0.63
            elseif(kplane.eq.iplane3)then
               VPB1=0.00
               VPT1=0.39
            endif

            CALL SET  (VPL1,VPR1,VPB1,VPT1,XWDL,XWDR,YWDB,YWDT,LNLG)

            CALL SRFACE2 (XX,YY,Z1,mscrat,ntracw,ntracw,nt,S,0.0)

         elseif(icol.eq.2)then

            VPL1=0.32
            VPR1=0.72
            if(kplane.eq.iplane1)then
               VPB1=0.54
               VPT1=0.94
            elseif(kplane.eq.iplane2)then
               VPB1=0.31
               VPT1=0.71
            elseif(kplane.eq.iplane3)then
               VPB1=0.07
               VPT1=0.47
            endif

            CALL SET  (VPL1,VPR1,VPB1,VPT1,XWDL,XWDR,YWDB,YWDT,LNLG)

            CALL SRFACE2 (XX,YY,Z1,mscrat,ntracw,ntracw,nt,S,0.)

            ISIZE=3

         elseif(icol.eq.3)then

            VPL1=0.60
            VPR1=1.00
            if(kplane.eq.iplane1)then
               VPB1=0.61
               VPT1=1.00
            elseif(kplane.eq.iplane2)then
               VPB1=0.39
               VPT1=0.79
            elseif(kplane.eq.iplane3)then
               VPB1=0.15
               VPT1=0.55
            endif

            CALL SET  (VPL1,VPR1,VPB1,VPT1,XWDL,XWDR,YWDB,YWDT,LNLG)

            CALL SRFACE2 (XX,YY,Z1,mscrat,ntracw,ntracw,nt,S,0.)

            ISIZE=3

         endif

 1111 continue

c     Now loop on rows

      if(kplane.eq.iplane1)then
         jplane=iplane1+1
         kplane=iplane2
         go to 1112
      elseif(kplane.eq.iplane2)then
         jplane=iplane2+1
         kplane=iplane3
         go to 1112
      endif

c     Advance the frame.

      write(*,*)'Finished.'

      CALL FRAME

      close(24)

c     Close GKS.

      CALL GRCLSE(ANSM,UNIT1,UNIT2)

 69   stop
      end

c***********************************************************************
      SUBROUTINE GRINIT(ANSM,UNIT1,UNIT2)
c***********************************************************************

      CHARACTER*1 ANSM

      INTEGER UNIT1,UNIT2

      IF(ANSM.EQ.'Y'.OR.ANSM.EQ.'y')THEN

c     OPEN GKS, OPEN WORKSTATION, ACTIVATE WORKSTATION

         CALL GOPKS(UNIT1,1024)

         CALL GOPWK(1,0,1)
         CALL GACWK(1)

c     OPEN AND ACTIVATE METAFILE

         CALL GOPWK(4,0,4)
         CALL GACWK(4)

c     AND WISS FOR SEGMENTS

         CALL GOPWK(0,0,0)
         CALL GACWK(0)

      ELSEIF(ANSM.EQ.'N'.OR.ANSM.EQ.'n')THEN

c     SUN4/60 VERSION

         OPEN(UNIT1,FILE='ERRORS')
         OPEN(UNIT2,FILE='METAFILE')

c     SET WORKSTATION SIZE ETC.

         CALL GOPKS(UNIT1,1024)

c     Open and activate workstation for X11 

         CALL GOPWK(3,0,8)
         CALL GACWK(3)
         CALL NGSETI('SC',2)

c     OPEN AND ACTIVATE METAFILE

c         CALL GACWK(5)

c     OPEN AND ACTIVATE POSTSCRIPT DEVICE

c         CALL GOPWK(4,UNIT2,10)
c         CALL GACWK(4)

         call NGSETC('ME','./Output/wavepanel.ps')
         CALL GOPWK(2,2,20)
         CALL GACWK(2)

c     AND WISS FOR SEGMENTS

c         CALL GOPWK(0,0,1)
c         CALL GACWK(0)

      ENDIF

      RETURN
      END

c***********************************************************************
      SUBROUTINE GRCLSE(ANSM,UNIT1,UNIT2)
c***********************************************************************

      CHARACTER*1 ANSM
      INTEGER UNIT1,UNIT2

      CALL GDAWK(2)
      CALL GCLWK(2)
      CALL GDAWK(3)
      CALL GCLWK(3)
      CALL GCLKS

      IF(ANSM.EQ.'Y'.OR.ANSM.EQ.'y')THEN
c         CALL GERCLS(UNIT1)
      ELSEIF(ANSM.EQ.'N'.OR.ANSM.EQ.'n')THEN
         CLOSE(UNIT1)
         CLOSE(UNIT2)
      ENDIF

      RETURN
      END

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
        THETA=6.28318530717959D0/(ISIGN*MMAX)
        WPR=-2.D0*DSIN(0.5D0*THETA)**2
        WPI=DSIN(THETA)
        WR=1.D0
        WI=0.D0
        DO 13 M=1,MMAX,2
          DO 12 I=M,N,ISTEP
            J=I+MMAX
c           TEMPR=SNGL(WR)*DATA(J)-SNGL(WI)*DATA(J+1)
c           TEMPI=SNGL(WR)*DATA(J+1)+SNGL(WI)*DATA(J)
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

      pi=4.d0*datan(1.d0)
  
      if(ifilt.eq.1)then                   !Low-pass filter

         rat1=dble(iomc4-iomc3)/pi

         do iom=1,nom2

            if(iomc4.eq.iomc3)then

               rfilt(iom)=1.d0

            else

               if(iom.lt.iomc3)then
                  rfilt(iom)=1.d0
               elseif(iom.ge.iomc3.and.iom.le.iomc4)then
                  rat=dble(iom-iomc3)/rat1
                  rfilt(iom)=1.d0-(1.d0-dcos(rat))/2.d0
               else
                  rfilt(iom)=0.d0
               endif

            endif

         enddo

      elseif(ifilt.eq.2)then               !Band-pass filter

         rat1=dble(iomc2-iomc1)/pi
         rat2=dble(iomc4-iomc3)/pi

         do iom=1,nom2

            if(iomc1.eq.iomc2)then
               rfilt(iom)=1.d0
            elseif(iomc3.eq.iomc4)then
               rfilt(iom)=1.d0
            else

               if(iom.lt.iomc1)then
                  rfilt(iom)=0.d0
               elseif(iom.ge.iomc1.and.iom.le.iomc2)then
                  rat=dble(iom-iomc1)/rat1
                  rfilt(iom)=(1.d0-dcos(rat))/2.d0
               elseif(iom.gt.iomc2.and.iom.lt.iomc3)then
                  rfilt(iom)=1.d0
               elseif(iom.ge.iomc3.and.iom.le.iomc4)then
                  rat=dble(iom-iomc3)/rat1
                  rfilt(iom)=1.d0-(1.d0-dcos(rat))/2.d0
               else
                  rfilt(iom)=0.d0
               endif

            endif

         enddo

      elseif(ifilt.eq.3)then               !High-pass filter

         rat1=dble(iomc2-iomc1)/pi

         do iom=1,nom2

            if(iomc2.eq.iomc1)then

               rfilt(iom)=1.d0

            else

               if(iom.lt.iomc1)then
                  rfilt(iom)=0.d0
               elseif(iom.ge.iomc1.and.iom.le.iomc2)then
                  rat=dble(iom-iomc2)/rat1
                  rfilt(iom)=1.d0-(1.d0-dcos(rat))/2.d0
               else
                  rfilt(iom)=1.d0
               endif

            endif

         enddo

      elseif(ifilt.eq.4)then

         do iom=1,nom2

            if(iomc1.eq.iomc2)then
               rfilt(iom)=1.d0
            elseif(iomc3.eq.iomc4)then
               rfilt(iom)=1.d0
            else

               if(iom.lt.iomc1)then
                  rfilt(iom)=1.d0
               elseif(iom.ge.iomc1.and.iom.le.iomc2)then
                  rat=1.d0-dble(iom-iomc1)/rat1
                  rfilt(iom)=(1.d0-dcos(rat))/2.d0
               elseif(iom.gt.iomc2.and.iom.lt.iomc3)then
                  rfilt(iom)=0.d0
               elseif(iom.ge.iomc3.and.iom.le.iomc4)then
                  rat=1.d0-dble(iom-iomc3)/rat1
                  rfilt(iom)=1.d0-(1.d0-dcos(rat))/2.d0
               else
                  rfilt(iom)=1.d0
               endif

            endif

         enddo

      endif

      return
      end

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


C REVERSED FROM PREVIOiUS STATEMENT !!!

      CALL GSCR (IWKID,1,0.,0.,0.)
      CALL GSCR (IWKID,0,RGBV(1,1),RGBV(2,1),RGBV(3,1))

      DO 101 I=2,15
         CALL GSCR (IWKID,I,RGBV(1,I),RGBV(2,I),RGBV(3,I))
 101  CONTINUE

C     Done.

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

      CALL GSFACI(0)

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
