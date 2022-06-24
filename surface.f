c........+.........+.........+.........+.........+.........+.........+..
c  File Surface contains the following 2 subroutines: Srface2 and 
c  Srfgk2.
c
c  Last modified: May 31, 2002.
c***********************************************************************
      SUBROUTINE SRFACE2 (X,Y,Z,M,MX,NX,NY,S,STEREO)
c***********************************************************************
C +-----------------------------------------------------------------+
C |                                                                 |
C |                Copyright (C) 1989 by UCAR                       |
C |        University Corporation for Atmospheric Research          |
C |                    All Rights Reserved                          |
C |                                                                 |
C |                 NCARGRAPHICS  Version 3.00                      |
C |                                                                 |
C +-----------------------------------------------------------------+
C
C DIMENSION OF           X(NX),Y(NY),Z(MX,NY),M(2,NX,NY),S(6)
C ARGUMENTS
C
C PURPOSE                SRFACE draws a perspective picture of a
C                        function of two variables with hidden lines
C                        removed.  The function is approximated by a
C                        two-dimensional array of heights.
C
C USAGE                  If the following assumptions are met, use
C
C                             CALL EZSRFC (Z,M,N,ANGH,ANGV,WORK)
C
C                          assumptions:
C                              .the entire array is to be drawn,
C                              .the data are equally spaced (in the
C                               X-Y plane),
C                              .no stereo pairs,
C                              .scaling is chosen internally.
C
C                        If these assumptions are not met use
C
C                             CALL SRFACE (X,Y,Z,M,MX,NX,NY,S,
C                                          STEREO)
C
C ARGUMENTS
C
C ON INPUT               Z
C FOR EZSRFC               The M by N array to be drawn.
C
C                        M
C                          The first dimension of Z.
C
C                        N
C                          The second dimension of Z.
C
C                        ANGH
C                          Angle in degrees in the X-Y plane to the
C                          line of sight (counter-clockwise from
C                          the plus-X axis).
C
C                        ANGV
C                          Angle in degrees from the X-Y plane to
C                          the line of sight (positive angles are
C                          above the middle Z, negative below).
C
C                        WORK
C                          A scratch storage dimensioned at least
C                          2*M*N+M+N.
C
C ON OUTPUT              Z, M, N, ANGH, ANGV are unchanged.  WORK
C FOR EZSRFC             has been written in.
C
C
C ARGUMENTS
C
C ON INPUT               X
C FOR SRFACE               A linear array NX long containing the X
C                          coordinates of the points in the surface
C                          approximation, in increasing numerical order.
C                          See note, below.
C
C                        Y
C                          A linear array NY long containing the Y
C                          coordinates of the points in the surface
C                          approximation, in increasing numerical order.
C                          See note, below.
C
C                        Z
C                          An array MX by NY containing the surface
C                          to be drawn in NX by NY cells.
C                          Z(I,J) = F(X(I),Y(J)).  See note, below.
C
C                        M
C                          Scratch array at least 2*NX*NY words
C                          long.
C
C                        MX
C                          First dimension of Z.
C
C                        NX
C                          Number of data values in the X direction
C                          (the first subscript direction) in Z.
C                          to be plotted. When plotting an entire
C                          array, MX=NX.
C
C                        NY
C                          Number of data values in the Y direction
C                          (the second subscript direction) to be
C                          plotted.
C
C                        S
C                          S defines the line of sight.  The viewer's
C                          eye is at (S(1), S(2), S(3)) and the
C                          point looked at is at (S(4), S(5), S(6)).
C                          The eye should be outside the block with
C                          opposite corners (X(1), Y(1), ZMIN) and
C                          (X(NX), Y(NY), ZMAX) and the point looked
C                          at should be inside it.  For a nice
C                          perspective effect, the distance between
C                          the eye and the point looked at should be
C                          5 to 10 times the size of the block.  See
C                          note, below.
C
C                        STEREO
C                          Flag to indicate if stereo pairs are to
C                          be drawn.  0.0 means no stereo pair (one
C                          picture).  Non-zero means put out two
C                          pictures.  The value of STEREO is the
C                          relative angle between the eyes.  A value
C                          of 1.0 produces standard separation.
C                          Negative STEREO reverses the left and
C                          right figures.  See the documentation below
C                          for internal variable ISTP for additional
C                          information.
C
C ON OUTPUT              X, Y, Z, MX, NX, NY, S, STEREO are
C FOR SRFACE             unchanged.  M has been written in.
C
C NOTES                  . The range of Z compared with the range
C                          of X and Y determines the shape of the
C                          picture.  They are assumed to be in the
C                          same units and not wildly different in
C                          magnitude.  S is assumed to be in the
C                          same units as X, Y, and Z.
C                        . Picture size can be made relative to
C                          distance.  See comments in SETR.
C                        . TRN32S can be used to translate from 3
C                          space to 2 space.  See comments there.
C                        . Data with extreme discontinuities may
C                          cause visibility errors.  If this problem
C                          occurs, use a distant eye position
C                          away from the +Z axis.
C                        . The default line color is set to
C                          color index 1.  If the user wishes to
C                          change the line color, he can do so by
C                          defining color index 1 before calling
C                          SRFACE, or by putting the common block
C                          SRFINT in his calling program and
C                          defining and using color index ISRFMJ
C                          (defaulted to 1 in BLOCKDATA.)
C
C ENTRY POINTS           SRFACE, SRFGK, EZSRFC, SETR, DRAWS, TRN32S,
C                        CLSET, CTCELL, SRFABD
C
C COMMON BLOCKS          SRFACE, SRFINT, SRFBLK, PWRZIS, SRFIP1
C
C REQUIRED LIBRARY       The SPPS.
C ROUTINES
C
C REQUIRED GKS LEVEL     0A
C
C I/O                    Plots
C
C PRECISION              Single
C
C LANGUAGE               FORTRAN
C
C HISTORY                Converted to FORTRAN 77 and GKS in March 1984.
C
C                        Prepared for SIGGRAPH, August 1976.
C
C                        Standardized in January 1973.
C
C                        Written in December 1971.  Replaced K.S.+G.
C                        algorithm called SOLIDS at NCAR.
C
C
C ALGORITHM              The data are processed from the near side of
C                        of the surface to the far side.  Visibility
C                        information is stored (see reference.)
C                        Highest so far is visible from above.
C
C REFERENCE              Wright, T.J., A Two Space Solution to the
C                        Hidden Line Problem for Plotting a Function
C                        of Two Variables.  IEEE Trans. Comp.,
C                        pp 28-33, January 1973.
C
C ACCURACY               If the ends of a line segment are visible,
C                        the middle is assumed visible.
C
C TIMING                 Proportional to NX*NY.
C
C
C INTERNAL PARAMETERS    name   default  function
C                        ----   -------  --------
C                        IFR        1    -1  Call FRAME first.
C                                         0  Do not call FRAME.
C                                        +1  Call FRAME when done.
C
C                        ISTP       0    STEREO type if STEREO
C                                        non-zero.
C                                        -1  Alternating frames,
C                                            slightly offset (for
C                                            movies,  IROTS = 0).
C                                         0  Blank frame between
C                                            for stereo slide.
C                                            IROTS = 1).
C                                        +1  Both on same frame.
C                                            (left picture to left
C                                            side.  IROTS = 0).
C
C                        IROTS      0     0  +Z in vertical plotting
C                                            direction (CINE mode).
C                                        +1  +Z in horizontal
C                                            plotting direction
C                                            (COMIC mode).
C
C                        IDRX       1    +1  Draw lines of constant
C                                            X.
C                                         0  Do not.
C
C                        IDRY       1    +1  Draw lines of constant
C                                            Y.
C                                         0  Do not.
C
C                        IDRZ       0    +1  Draw lines of constant
C                                            Z (contour lines).
C                                         0  Do not.
C
C                        IUPPER     0    +1  Draw upper side of
C                                            surface.
C                                         0  Draw both sides.
C                                        -1  Draw lower side.
C
C                        ISKIRT     0    +1  Draw a skirt around the
C                                            surface.
C                                            BOTTOM = HSKIRT.
C                                         0  Do not.
C
C                        NCLA       6    Approximate number of
C                                        levels of constant Z that
C                                        are drawn if levels are not
C                                        specified.  40 levels
C                                        maximum.
C
C                        THETA    .02    Angle, in radians, between
C                                        eyes for stereo pairs.
C
C                        HSKIRT    0.    Height of skirt
C                                            (if ISKIRT = 1).
C
C                        CHI       0.    Highest level of constant
C                                        Z.
C
C                        CLO       0.    Lowest level of constant Z.
C
C                        CINC      0.    Increment between levels.
C
C                          [If CHI, CLO, or CINC is zero, a nice
C                           value is generated automatically.]
C
C                        IOFFP     0     Flag to control use of special
C                                        value feature.  Do not have
C                                        both IOFFP=1 and ISKIRT=1.
C                                         0 Feature not in use
C                                        +1  Feature in use.  No lines
C                                            drawn to data points in Z
C                                            that are equal to SPVAL.
C
C                        SPVAL     0.    Special value used to mark un-
C                                        known data when IOFFP=1.
C
C
C
c***********************************************************************

      EXTERNAL        SRFABD
C
      DIMENSION       X(NX)      ,Y(NY)      ,Z(MX,NY)   ,M(2,NX,NY) ,
     1                S(6)
      DIMENSION       WIN1(4)    ,VP1(4)     ,LASF(13)
      COMMON /SRFINT/ ISRFMJ     ,ISRFMN     ,ISRFTX
      CALL Q8QST4 ('GRAPHX','SRFACE','SRFACE','VERSION 01')
C
C     THIS DRIVER SAVES THE CURRENT NORMALIZATION TRANSFORMATION
C     INFORMATION, DEFINES THE NORMALIZATION TRANSFORMATION
C     APPROPRIATE FOR SRFGK, CALLS SRFGK, AND RESTORES THE ORIGINAL
C     NORMALIZATION TRANSFORMATION.
C
C     GET CURRENT NORMALIZATION TRANSFORMATION NUMBER
C
      CALL GQCNTN (IER,NTORIG)
C
C     STORE WINDOW AND VIEWPORT OF NORMALIZATION TRANSFORMATION 1
C
      CALL GQNT (NTORIG,IER,WIN1,VP1)
      CALL GETUSV('LS',IOLLS)
C
C     SET WINDOW AND VIEWPORT FOR SRFGK
C
c     CALL SET(0.,1.,0.,1.,1.,1024.,1.,1024.,1)
C
C     SET LINE COLOR TO INDIVIDUAL (SAVE CURRENT SETTING)
C
      CALL GQASF (IER,LASF)
      LASFSV  = LASF(3)
      LASF(3) = 1
      CALL GSASF(LASF)
C
C     SET LINE COLOR INDEX TO COMMON VARIABLE ISRFMJ (SAVE
C     CURRENT SETTING)
C
      CALL GQPLCI (IER,LCISV)
      CALL GSPLCI (ISRFMJ)
C
C     DRAW PLOT
C
      CALL SRFGK2 (X,Y,Z,M,MX,NX,NY,S,STEREO)
C
C     RESTORE INITIAL LINE COLOR SETTINGS
C
      LASF(3) = LASFSV
      CALL GSASF(LASF)
      CALL GSPLCI (LCISV)
C
C     RESTORE ORIGINAL NORMALIZATION TRANSFORMATION
C
      CALL SET(VP1(1),VP1(2),VP1(3),VP1(4),WIN1(1),WIN1(2),
     -         WIN1(3),WIN1(4),IOLLS)
      CALL GSELNT (NTORIG)
C
      RETURN
      END

c***********************************************************************
      SUBROUTINE SRFGK2 (X,Y,Z,M,MX,NX,NY,S,STEREO)
c***********************************************************************

      DIMENSION       X(NX)      ,Y(NY)      ,Z(MX,NY)   ,M(2,NX,NY) ,
     1                S(6)
      DIMENSION       MXS(2)     ,MXF(2)     ,MXJ(2)     ,MYS(2)     ,
     1                MYF(2)     ,MYJ(2)
      COMMON /SRFBLK/ LIMU(1024) ,LIML(1024) ,CL(41)     ,NCL        ,
     1                LL         ,FACT       ,IROT       ,NDRZ       ,
     2                NUPPER     ,NRSWT      ,BIGD       ,UMIN       ,
     3                UMAX       ,VMIN       ,VMAX       ,RZERO      ,
     4                IOFFP      ,NSPVAL     ,SPVAL      ,BIGEST
      COMMON /PWRZ1S/ XXMIN      ,XXMAX      ,YYMIN      ,YYMAX      ,
     1                ZZMIN      ,ZZMAX      ,DELCRT     ,EYEX       ,
     2                EYEY       ,EYEZ
      COMMON /SRFIP1/ IFR        ,ISTP       ,IROTS      ,IDRX       ,
     1                IDRY       ,IDRZ       ,IUPPER     ,ISKIRT     ,
     2                NCLA       ,THETA      ,HSKIRT     ,CHI        ,
     3                CLO        ,CINC       ,ISPVAL
      DATA JF,IF,LY,LX,ICNST/1,1,2,2,0/
      CALL Q8QST4 ('GRAPHX','SRFACE','SRFGK','VERSION 01')
      BIGEST = R1MACH(2)
      MMXX = MX
      NNXX = NX
      NNYY = NY
      STER = STEREO
      NXP1 = NNXX+1
      NYP1 = NNYY+1
      NLA = NCLA
      NSPVAL = ISPVAL
      NDRZ = IDRZ
      IF (IDRZ .NE. 0)
     1    CALL CLSET (Z,MMXX,NNXX,NNYY,CHI,CLO,CINC,NLA,40,CL,NCL,
     2                ICNST,IOFFP,SPVAL,BIGEST)
      IF (IDRZ .NE. 0) NDRZ = 1-ICNST
      STHETA = SIN(STER*THETA)
      CTHETA = COS(STER*THETA)
      RX = S(1)-S(4)
      RY = S(2)-S(5)
      RZ = S(3)-S(6)
      D1 = SQRT(RX*RX+RY*RY+RZ*RZ)
      D2 = SQRT(RX*RX+RY*RY)
      DX = 0.
      DY = 0.
      IF (STEREO .EQ. 0.) GO TO  20
      D1 = D1*STEREO*THETA
      IF (D2 .GT. 0.) GO TO  10
      DX = D1
      GO TO  20
   10 AGL = ATAN2(RX,-RY)
      DX = D1*COS(AGL)
      DY = D1*SIN(AGL)
   20 IROT = IROTS
      NPIC = 1
      IF (STER .NE. 0.) NPIC = 2
      FACT = 1.
      IF (NRSWT .NE. 0) FACT = RZERO/D1
      IF (ISTP.EQ.0 .AND. STER.NE.0.) IROT = 1
      DO 570 IPIC=1,NPIC
         NUPPER = IUPPER
         IF (IFR .LT. 0) CALL FRAME
C
C SET UP MAPING FROM FLOATING POINT 3-SPACE TO CRT SPACE.
C
         SIGN1 = IPIC*2-3
         EYEX = S(1)+SIGN1*DX
         POIX = S(4)+SIGN1*DX
         EYEY = S(2)+SIGN1*DY
         POIY = S(5)+SIGN1*DY
         EYEZ = S(3)
         POIZ = S(6)
         LL = 0
         XEYE = EYEX
         YEYE = EYEY
         ZEYE = EYEZ
         CALL TRN32S (POIX,POIY,POIZ,XEYE,YEYE,ZEYE,0)
         LL = IPIC+2*ISTP+3
         IF (STER .EQ. 0.) LL = 1
         IF (NRSWT .NE. 0) GO TO 100
         XXMIN = X(1)
         XXMAX = X(NNXX)
         YYMIN = Y(1)
         YYMAX = Y(NNYY)
         UMIN = BIGEST
         VMIN = BIGEST
         ZZMIN = BIGEST
         UMAX = -UMIN
         VMAX = -VMIN
         ZZMAX = -ZZMIN
         DO  40 J=1,NNYY
            DO  30 I=1,NNXX
               ZZ = Z(I,J)
               IF (IOFFP.EQ.1 .AND. ZZ.EQ.SPVAL) GO TO  30
               ZZMAX = AMAX1(ZZMAX,ZZ)
               ZZMIN = AMIN1(ZZMIN,ZZ)
               CALL TRN32S (X(I),Y(J),Z(I,J),UT,VT,DUMMY,1)
               UMAX = AMAX1(UMAX,UT)
               UMIN = AMIN1(UMIN,UT)
               VMAX = AMAX1(VMAX,VT)
               VMIN = AMIN1(VMIN,VT)
   30       CONTINUE
   40    CONTINUE
         IF (ISKIRT .NE. 1) GO TO  70
         NXSTP = NNXX-1
         NYSTP = NNYY-1
         DO  60 J=1,NNYY,NYSTP
            DO  50 I=1,NNXX,NXSTP
               CALL TRN32S (X(I),Y(J),HSKIRT,UT,VT,DUMMY,1)
               UMAX = AMAX1(UMAX,UT)
               UMIN = AMIN1(UMIN,UT)
               VMAX = AMAX1(VMAX,VT)
               VMIN = AMIN1(VMIN,VT)
   50       CONTINUE
   60    CONTINUE
   70    CONTINUE
         WIDTH = UMAX-UMIN
         HIGHT = VMAX-VMIN
         DIF = .5*(WIDTH-HIGHT)
         IF (DIF)  80,100, 90
   80    UMIN = UMIN+DIF
         UMAX = UMAX-DIF
         GO TO 100
   90    VMIN = VMIN-DIF
         VMAX = VMAX+DIF
  100    XEYE = EYEX
         YEYE = EYEY
         ZEYE = EYEZ
         CALL TRN32S (POIX,POIY,POIZ,XEYE,YEYE,ZEYE,0)
         DO 120 J=1,NNYY
            DO 110 I=1,NNXX
               CALL TRN32S (X(I),Y(J),Z(I,J),UT,VT,DUMMY,1)
               M(1,I,J) = UT
               M(2,I,J) = VT
  110       CONTINUE
  120    CONTINUE
C
C INITIALIZE UPPER AND LOWER VISIBILITY ARRAYS
C
         DO 130 K=1,1024
            LIMU(K) = 0
            LIML(K) = 1024
  130    CONTINUE
C
C FIND ORDER TO DRAW LINES
C
         NXPASS = 1
         IF (S(1) .GE. X(NNXX)) GO TO 160
         IF (S(1) .LE. X(1)) GO TO 170
         DO 140 I=2,NNXX
            LX = I
            IF (S(1) .LE. X(I)) GO TO 150
  140    CONTINUE
  150    MXS(1) = LX-1
         MXJ(1) = -1
         MXF(1) = 1
         MXS(2) = LX
         MXJ(2) = 1
         MXF(2) = NNXX
         NXPASS = 2
         GO TO 180
  160    MXS(1) = NNXX
         MXJ(1) = -1
         MXF(1) = 1
         GO TO 180
  170    MXS(1) = 1
         MXJ(1) = 1
         MXF(1) = NNXX
  180    NYPASS = 1
         IF (S(2) .GE. Y(NNYY)) GO TO 210
         IF (S(2) .LE. Y(1)) GO TO 220
         DO 190 J=2,NNYY
            LY = J
            IF (S(2) .LE. Y(J)) GO TO 200
  190    CONTINUE
  200    MYS(1) = LY-1
         MYJ(1) = -1
         MYF(1) = 1
         MYS(2) = LY
         MYJ(2) = 1
         MYF(2) = NNYY
         NYPASS = 2
         GO TO 230
  210    MYS(1) = NNYY
         MYJ(1) = -1
         MYF(1) = 1
         GO TO 230
  220    MYS(1) = 1
         MYJ(1) = 1
         MYF(1) = NNYY
C
C PUT ON SKIRT ON FRONT SIDE IF WANTED
C
  230    IF (NXPASS.EQ.2 .AND. NYPASS.EQ.2) GO TO 490
         IF (ISKIRT .EQ. 0) GO TO 290
         IN = MXS(1)
         IF = MXF(1)
         JN = MYS(1)
         JF = MYF(1)
         IF (NYPASS .NE. 1) GO TO 260
         CALL TRN32S (X(1),Y(JN),HSKIRT,UX1,VX1,DUMMY,1)
         CALL TRN32S (X(NNXX),Y(JN),HSKIRT,UX2,VX2,DUMMY,1)
         QU = (UX2-UX1)/(X(NNXX)-X(1))
         QV = (VX2-VX1)/(X(NNXX)-X(1))
         YNOW = Y(JN)
         DO 240 I=1,NNXX
            CALL TRN32S (X(I),YNOW,HSKIRT,RU,RV,DUMMY,1)
            CALL DRAWS (IFIX(RU),IFIX(RV),M(1,I,JN),M(2,I,JN),1,1)
  240    CONTINUE
         CALL DRAWS (IFIX(UX1),IFIX(VX1),IFIX(UX2),IFIX(VX2),1,1)
         IF (IDRY .NE. 0) GO TO 260
         DO 250 I=2,NNXX
            CALL DRAWS (M(1,I-1,JN),M(2,I-1,JN),M(1,I,JN),M(2,I,JN),1,1)
  250    CONTINUE
  260    IF (NXPASS .NE. 1) GO TO 290
         CALL TRN32S (X(IN),Y(1),HSKIRT,UY1,VY1,DUMMY,1)
         CALL TRN32S (X(IN),Y(NNYY),HSKIRT,UY2,VY2,DUMMY,1)
         QU = (UY2-UY1)/(Y(NNYY)-Y(1))
         QV = (VY2-VY1)/(Y(NNYY)-Y(1))
         XNOW = X(IN)
         DO 270 J=1,NNYY
            CALL TRN32S (XNOW,Y(J),HSKIRT,RU,RV,DUMMY,1)
            CALL DRAWS (IFIX(RU),IFIX(RV),M(1,IN,J),M(2,IN,J),1,1)
  270    CONTINUE
         CALL DRAWS (IFIX(UY1),IFIX(VY1),IFIX(UY2),IFIX(VY2),1,1)
         IF (IDRX .NE. 0) GO TO 290
         DO 280 J=2,NNYY
            CALL DRAWS (M(1,IN,J-1),M(2,IN,J-1),M(1,IN,J),M(2,IN,J),1,1)
  280    CONTINUE
C
C PICK PROPER ALGORITHM
C
  290    LI = MXJ(1)
         MI = MXS(1)-LI
         NI = IABS(MI-MXF(1))
         LJ = MYJ(1)
         MJ = MYS(1)-LJ
         NJ = IABS(MJ-MYF(1))
C
C IF NXPASS IS 1, PUT THE I LOOP OUTERMOST.  OTHERWISE, PUT THE J LOOP
C OUTERMOST.
C
         IF (NXPASS.EQ.2) GO TO 360
         IF (ISKIRT.NE.0 .OR. NYPASS.NE.1) GO TO 310
         I = MXS(1)
         DO 300 J=2,NNYY
            CALL DRAWS (M(1,I,J-1),M(2,I,J-1),M(1,I,J),M(2,I,J),0,1)
  300    CONTINUE
  310    DO 350 II=1,NNXX
            I = MI+II*LI
            IPLI = I+LI
            IF (NYPASS .EQ. 1) GO TO 320
            K = MYS(1)
            L = MYS(2)
            IF (IDRX .NE. 0)
     1          CALL DRAWS (M(1,I,K),M(2,I,K),M(1,I,L),M(2,I,L),1,1)
            IF (NDRZ.NE.0 .AND. II.NE.NI)
     1          CALL CTCELL (Z,MMXX,NNXX,NNYY,M,MIN0(I,I+LI),K)
  320       DO 340 JPASS=1,NYPASS
               LJ = MYJ(JPASS)
               MJ = MYS(JPASS)-LJ
               NJ = IABS(MJ-MYF(JPASS))
               DO 330 JJ=1,NJ
                  J = MJ+JJ*LJ
                  JPLJ = J+LJ
                  IF (IDRX.NE.0 .AND. JJ.NE.NJ)
     1                CALL DRAWS (M(1,I,J),M(2,I,J),M(1,I,JPLJ),
     2                            M(2,I,JPLJ),1,1)
                  IF (I.NE.MXF(1) .AND. IDRY.NE.0)
     1                CALL DRAWS (M(1,IPLI,J),M(2,IPLI,J),M(1,I,J),
     2                            M(2,I,J),1,1)
                  IF (NDRZ.NE.0 .AND. JJ.NE.NJ .AND. II.NE.NNXX)
     1                CALL CTCELL (Z,MMXX,NNXX,NNYY,M,MIN0(I,I+LI),
     2                             MIN0(J,J+LJ))
  330          CONTINUE
  340       CONTINUE
  350    CONTINUE
         GO TO 430
  360    IF (ISKIRT.NE.0 .OR. NXPASS.NE.1) GO TO 380
         J = MYS(1)
         DO 370 I=2,NNXX
            CALL DRAWS (M(1,I-1,J),M(2,I-1,J),M(1,I,J),M(2,I,J),0,1)
  370    CONTINUE
  380    DO 420 JJ=1,NNYY
            J = MJ+JJ*LJ
            JPLJ = J+LJ
            IF (NXPASS .EQ. 1) GO TO 390
            K = MXS(1)
            L = MXS(2)
            IF (IDRY .NE. 0)
     1          CALL DRAWS (M(1,K,J),M(2,K,J),M(1,L,J),M(2,L,J),1,1)
            IF (NDRZ.NE.0 .AND. JJ.NE.NJ)
     1          CALL CTCELL (Z,MMXX,NNXX,NNYY,M,K,MIN0(J,J+LJ))
  390       DO 410 IPASS=1,NXPASS
               LI = MXJ(IPASS)
               MI = MXS(IPASS)-LI
               NI = IABS(MI-MXF(IPASS))
               DO 400 II=1,NI
                  I = MI+II*LI
                  IPLI = I+LI
                  IF (IDRY.NE.0 .AND. II.NE.NI)
     1                CALL DRAWS (M(1,I,J),M(2,I,J),M(1,IPLI,J),
     2                            M(2,IPLI,J),1,1)
                  IF (J.NE.MYF(1) .AND. IDRX.NE.0)
     1                CALL DRAWS (M(1,I,JPLJ),M(2,I,JPLJ),M(1,I,J),
     2                            M(2,I,J),1,1)
                  IF (NDRZ.NE.0 .AND. II.NE.NI .AND. JJ.NE.NNYY)
     1                CALL CTCELL (Z,MMXX,NNXX,NNYY,M,MIN0(I,I+LI),
     2                             MIN0(J,J+LJ))
  400          CONTINUE
  410       CONTINUE
  420    CONTINUE
  430    IF (ISKIRT .EQ. 0) GO TO 520
C
C FIX UP IF SKIRT IS USED WITH LINES ONE WAY.
C
         IF (IDRX .NE. 0) GO TO 460
         DO 450 IPASS=1,NXPASS
            IF (NXPASS .EQ. 2) IF = 1+(IPASS-1)*(NNXX-1)
            DO 440 J=2,NNYY
               CALL DRAWS (M(1,IF,J-1),M(2,IF,J-1),M(1,IF,J),M(2,IF,J),
     1                     1,0)
  440       CONTINUE
  450    CONTINUE
  460    IF (IDRY .NE. 0) GO TO 520
         DO 480 JPASS=1,NYPASS
            IF (NYPASS .EQ. 2) JF = 1+(JPASS-1)*(NNYY-1)
            DO 470 I=2,NNXX
               CALL DRAWS (M(1,I-1,JF),M(2,I-1,JF),M(1,I,JF),M(2,I,JF),
     1                     1,0)
  470       CONTINUE
  480    CONTINUE
         GO TO 520
C
C ALL VISIBLE IF VIEWED FROM DIRECTLY ABOVE OR BELOW.
C
  490    IF (NUPPER.GT.0 .AND. S(3).LT.S(6)) GO TO 520
         IF (NUPPER.LT.0 .AND. S(3).GT.S(6)) GO TO 520
         NUPPER = 1
         IF (S(3) .LT. S(6)) NUPPER = -1
         DO 510 I=1,NNXX
            DO 500 J=1,NNYY
               IF (IDRX.NE.0 .AND. J.NE.NNYY)
     1             CALL DRAWS (M(1,I,J),M(2,I,J),M(1,I,J+1),M(2,I,J+1),
     2                         1,0)
               IF (IDRY.NE.0 .AND. I.NE.NNXX)
     1             CALL DRAWS (M(1,I,J),M(2,I,J),M(1,I+1,J),M(2,I+1,J),
     2                         1,0)
               IF (IDRZ.NE.0 .AND. I.NE.NNXX .AND. J.NE.NNYY)
     1             CALL CTCELL (Z,MMXX,NNXX,NNYY,M,I,J)
  500       CONTINUE
  510    CONTINUE
  520    IF (STER .EQ. 0.) GO TO 560
         IF (ISTP) 540,530,550
  530    CALL FRAME
  540    CALL FRAME
         GO TO 570
  550    IF (IPIC .NE. 2) GO TO 570
  560    IF (IFR .GT. 0) CALL FRAME
  570 CONTINUE
      RETURN
      END
