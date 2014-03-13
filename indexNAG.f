      SUBROUTINE E02ADF(M,KPLUS1,NROWS,X,Y,W,WORK1,WORK2,A,S,IFAIL)
C
C     NAG LIBRARY SUBROUTINE  E02ADF
C
C     E02ADF  COMPUTES WEIGHTED LEAST-SQUARES POLYNOMIAL
C     APPROXIMATIONS TO AN ARBITRARY SET OF DATA POINTS.
C
C     FORSYTHE-CLENSHAW METHOD WITH MODIFICATIONS DUE TO
C     REINSCH AND GENTLEMAN.
C
C     USES NAG LIBRARY ROUTINE  P01AAF.
C     USES BASIC EXTERNAL FUNCTION  SQRT.
C
C     STARTED - 1973.
C     COMPLETED - 1976.
C     AUTHOR - MGC AND JGH.
C
C     WORK1  AND  WORK2  ARE WORKSPACE AREAS.
C     WORK1(1, R)  CONTAINS THE VALUE OF THE  R TH  WEIGHTED
C     RESIDUAL FOR THE CURRENT DEGREE  I.
C     WORK1(2, R)  CONTAINS THE VALUE OF  X(R)  TRANSFORMED
C     TO THE RANGE  -1  TO  +1.
C     WORK1(3, R)  CONTAINS THE WEIGHTED VALUE OF THE CURRENT
C     ORTHOGONAL POLYNOMIAL (OF DEGREE  I)  AT THE  R TH
C     DATA POINT.
C     WORK2(1, J)  CONTAINS THE COEFFICIENT OF THE CHEBYSHEV
C     POLYNOMIAL OF DEGREE  J - 1  IN THE CHEBYSHEV-SERIES
C     REPRESENTATION OF THE CURRENT ORTHOGONAL POLYNOMIAL
C     (OF DEGREE  I).
C     WORK2(2, J)  CONTAINS THE COEFFICIENT OF THE CHEBYSHEV
C     POLYNOMIAL OF DEGREE  J - 1  IN THE CHEBYSHEV-SERIES
C     REPRESENTATION OF THE PREVIOUS ORTHOGONAL POLYNOMIAL
C     (OF DEGREE  I - 1).
C
C     NAG COPYRIGHT 1975
C     MARK 5 RELEASE
C     MARK 6 REVISED  IER-84
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     CHECK THAT THE VALUES OF  M  AND  KPLUS1  ARE REASONABLE
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E02ADF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, KPLUS1, M, NROWS
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NROWS,KPLUS1), S(KPLUS1), W(M), WORK1(3,M),
     *                  WORK2(2,KPLUS1), X(M), Y(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPIP1, BETAI, BJ, BJP1, BJP2, CI, D, DF, DI,
     *                  DIM1, DJ, EPSR, FACTOR, PIJ, SIGMAI, WRPR,
     *                  WRPRSQ, X1, XCAPR, XM
      INTEGER           I, IERROR, IPLUS1, IPLUS2, J, JPLUS1, JPLUS2,
     *                  JREV, K, MDIST, R
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      IERROR = 4
      IF (KPLUS1.LT.1 .OR. M.LT.KPLUS1) GO TO 380
      K = KPLUS1 - 1
C
C     TEST THE VALIDITY OF THE DATA.
C
C     CHECK THAT THE WEIGHTS ARE STRICTLY POSITIVE.
C
      IERROR = 1
      DO 20 R = 1, M
         IF (W(R).LE.0.0D0) GO TO 380
   20 CONTINUE
C
C     CHECK THAT THE VALUES OF  X(R)  ARE NON-DECREASING AND
C     DETERMINE
C     THE NUMBER  (MDIST)  OF DISTINCT VALUES OF  X(R).
C
      IERROR = 2
      MDIST = 1
      DO 40 R = 2, M
         IF (X(R).LT.X(R-1)) GO TO 380
         IF (X(R).EQ.X(R-1)) GO TO 40
         MDIST = MDIST + 1
   40 CONTINUE
C
C     IF THE  X(R)  ALL HAVE THE SAME VALUE, I.E.  MDIST = 1,
C     THE NORMALIZATION OF THE INDEPENDENT VARIABLE IS NOT
C     POSSIBLE.
C
      IERROR = 3
      IF (MDIST.EQ.1) GO TO 380
C
C     IF THE NUMBER OF DISTINCT VALUES OF  X(R)  FAILS TO EXCEED
C     THE MAXIMUM DEGREE  K  THERE IS NO UNIQUE POLYNOMIAL
C     APPROXIMATION OF THAT DEGREE.
C
      IERROR = 4
      IF (MDIST.LE.K) GO TO 380
C
C     CHECK THAT  NROWS  HAS BEEN SET SUFFICIENTLY LARGE.
C
      IERROR = 5
      IF (NROWS.LT.KPLUS1) GO TO 380
      IERROR = 0
C
      X1 = X(1)
      XM = X(M)
      D = XM - X1
C
C     THE INITIAL VALUES  WORK1(1, R)  (R = 1, 2, ..., M)  OF THE
C     WEIGHTED RESIDUALS AND THE VALUES  WORK1(2, R)  (R = 1, 2,
C     ..., M)
C     OF THE NORMALIZED INDEPENDENT VARIABLE ARE COMPUTED.  NOTE
C     THAT
C     WORK1(2, R)  IS COMPUTED FROM THE EXPRESSION BELOW RATHER
C     THAN THE
C     MORE NATURAL FORM
C
C     (2.0*X(R) - X1 - XM)/D,
C
C     SINCE THE FORMER GUARANTEES THE COMPUTED VALUE TO DIFFER FROM
C     THE TRUE VALUE BY AT MOST  4.0*MACHINE ACCURACY,  WHEREAS THE
C     LATTER HAS NO SUCH GUARANTEE.
C
      DO 60 R = 1, M
         WORK1(1,R) = W(R)*Y(R)
         WORK1(2,R) = ((X(R)-X1)-(XM-X(R)))/D
   60 CONTINUE
      I = 1
      BETAI = 0.0D0
      DO 360 IPLUS1 = 1, KPLUS1
C
C        SET STARTING VALUES FOR DEGREE  I.
C
         IPLUS2 = IPLUS1 + 1
         IF (IPLUS1.EQ.KPLUS1) GO TO 100
         DO 80 JPLUS1 = IPLUS2, KPLUS1
            A(IPLUS1,JPLUS1) = 0.0D0
   80    CONTINUE
         WORK2(1,IPLUS2) = 0.0D0
         WORK2(2,IPLUS2) = 0.0D0
  100    ALPIP1 = 0.0D0
         CI = 0.0D0
         DI = 0.0D0
         A(I,IPLUS1) = 0.0D0
         WORK2(1,IPLUS1) = 1.0D0
         IF (KPLUS1.GT.1) WORK2(2,1) = WORK2(1,2)
         DO 260 R = 1, M
            XCAPR = WORK1(2,R)
C
C           THE WEIGHTED VALUE  WORK1(3, R)  OF THE ORTHOGONAL
C           POLYNOMIAL OF
C           DEGREE  I  AT  X = X(R)  IS COMPUTED BY RECURRENCE FROM ITS
C           CHEBYSHEV-SERIES REPRESENTATION.
C
            IF (IPLUS1.GT.1) GO TO 120
            WRPR = W(R)*0.5D0*WORK2(1,1)
            WORK1(3,R) = WRPR
            GO TO 240
  120       J = IPLUS2
            IF (XCAPR.GT.0.5D0) GO TO 200
            IF (XCAPR.GE.-0.5D0) GO TO 160
C
C           GENTLEMAN*S MODIFIED RECURRENCE.
C
            FACTOR = 2.0D0*(1.0D0+XCAPR)
            DJ = 0.0D0
            BJ = 0.0D0
            DO 140 JREV = 1, I
               J = J - 1
               DJ = WORK2(1,J) - DJ + FACTOR*BJ
               BJ = DJ - BJ
  140       CONTINUE
            WRPR = W(R)*(0.5D0*WORK2(1,1)-DJ+0.5D0*FACTOR*BJ)
            WORK1(3,R) = WRPR
            GO TO 240
C
C           CLENSHAW*S ORIGINAL RECURRENCE.
C
  160       FACTOR = 2.0D0*XCAPR
            BJP1 = 0.0D0
            BJ = 0.0D0
            DO 180 JREV = 1, I
               J = J - 1
               BJP2 = BJP1
               BJP1 = BJ
               BJ = WORK2(1,J) - BJP2 + FACTOR*BJP1
  180       CONTINUE
            WRPR = W(R)*(0.5D0*WORK2(1,1)-BJP1+0.5D0*FACTOR*BJ)
            WORK1(3,R) = WRPR
            GO TO 240
C
C           REINSCH*S MODIFIED RECURRENCE.
C
  200       FACTOR = 2.0D0*(1.0D0-XCAPR)
            DJ = 0.0D0
            BJ = 0.0D0
            DO 220 JREV = 1, I
               J = J - 1
               DJ = WORK2(1,J) + DJ - FACTOR*BJ
               BJ = BJ + DJ
  220       CONTINUE
            WRPR = W(R)*(0.5D0*WORK2(1,1)+DJ-0.5D0*FACTOR*BJ)
            WORK1(3,R) = WRPR
C
C           THE COEFFICIENT  CI  OF THE  I TH  ORTHOGONAL POLYNOMIAL AND
C           THE
C           COEFFICIENTS  ALPIP1  AND  BETA I  IN THE
C           THREE-TERM RECURRENCE RELATION FOR THE ORTHOGONAL
C           POLYNOMIALS ARE COMPUTED.
C
  240       WRPRSQ = WRPR**2
            DI = DI + WRPRSQ
            CI = CI + WRPR*WORK1(1,R)
            ALPIP1 = ALPIP1 + WRPRSQ*XCAPR
  260    CONTINUE
         CI = CI/DI
         IF (IPLUS1.NE.1) BETAI = DI/DIM1
         ALPIP1 = 2.0D0*ALPIP1/DI
C
C        THE WEIGHTED RESIDUALS  WORK1(1, R)  (R = 1, 2, ..., M)  FOR
C        DEGREE  I  ARE COMPUTED, TOGETHER WITH THEIR SUM OF SQUARES,
C        SIGMAI.
C
         SIGMAI = 0.0D0
         DO 280 R = 1, M
            EPSR = WORK1(1,R) - CI*WORK1(3,R)
            WORK1(1,R) = EPSR
            SIGMAI = SIGMAI + EPSR**2
  280    CONTINUE
C
C        THE ROOT MEAN SQUARE RESIDUAL  S(I + 1)  FOR DEGREE  I  IS
C        THEORETICALLY UNDEFINED IF  M = I + 1  (THE CONDITION FOR THE
C        POLYNOMIAL TO PASS EXACTLY THROUGH THE DATA POINTS).  SHOULD
C        THIS
C        CASE ARISE THE R.M.S. RESIDUAL IS SET TO ZERO.
C
         IF (IPLUS1.GE.M) GO TO 300
         DF = M - IPLUS1
         S(IPLUS1) = SQRT(SIGMAI/DF)
         GO TO 320
  300    S(IPLUS1) = 0.0D0
C
C        THE CHEBYSHEV COEFFICIENTS  A(I + 1, 1), A(I + 1, 2), ...,
C        A(I + 1, I + 1)  TOGETHER WITH THE COEFFICIENTS
C        WORK2(1, 1), WORK2(1, 2), ..., WORK2(1, I + 1),   IN THE
C        CHEBYSHEV-SERIES REPRESENTATION OF THE  I TH  ORTHOGONAL
C        POLYNOMIAL ARE COMPUTED.
C
  320    DO 340 JPLUS1 = 1, IPLUS1
            JPLUS2 = JPLUS1 + 1
            PIJ = WORK2(1,JPLUS1)
            A(IPLUS1,JPLUS1) = A(I,JPLUS1) + CI*PIJ
            IF (JPLUS2.GT.KPLUS1) GO TO 380
            WORK2(1,JPLUS1) = WORK2(1,JPLUS2) + WORK2(2,JPLUS1) -
     *                        ALPIP1*PIJ - BETAI*WORK2(2,JPLUS2)
            WORK2(2,JPLUS2) = PIJ
  340    CONTINUE
         DIM1 = DI
         I = IPLUS1
  360 CONTINUE
  380 IF (IERROR) 400, 420, 400
  400 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
  420 IFAIL = 0
      RETURN
      END
      INTEGER FUNCTION P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-621 (APR 1988).
C     MARK 13B REVISED. IER-668 (AUG 1988).
C
C     P01ABF is the error-handling routine for the NAG Library.
C
C     P01ABF either returns the value of IERROR through the routine
C     name (soft failure), or terminates execution of the program
C     (hard failure). Diagnostic messages may be output.
C
C     If IERROR = 0 (successful exit from the calling routine),
C     the value 0 is returned through the routine name, and no
C     message is output
C
C     If IERROR is non-zero (abnormal exit from the calling routine),
C     the action taken depends on the value of IFAIL.
C
C     IFAIL =  1: soft failure, silent exit (i.e. no messages are
C                 output)
C     IFAIL = -1: soft failure, noisy exit (i.e. messages are output)
C     IFAIL =-13: soft failure, noisy exit but standard messages from
C                 P01ABF are suppressed
C     IFAIL =  0: hard failure, noisy exit
C
C     For compatibility with certain routines included before Mark 12
C     P01ABF also allows an alternative specification of IFAIL in which
C     it is regarded as a decimal integer with least significant digits
C     cba. Then
C
C     a = 0: hard failure  a = 1: soft failure
C     b = 0: silent exit   b = 1: noisy exit
C
C     except that hard failure now always implies a noisy exit.
C
C     S.Hammarling, M.P.Hooper and J.J.du Croz, NAG Central Office.
C
C     .. Scalar Arguments ..
      INTEGER                 IERROR, IFAIL, NREC
      CHARACTER*(*)           SRNAME
C     .. Array Arguments ..
      CHARACTER*(*)           REC(*)
C     .. Local Scalars ..
      INTEGER                 I, NERR
      CHARACTER*72            MESS
C     .. External Subroutines ..
      EXTERNAL                P01ABZ, X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC               ABS, MOD
C     .. Executable Statements ..
      IF (IERROR.NE.0) THEN
C        Abnormal exit from calling routine
         IF (IFAIL.EQ.-1 .OR. IFAIL.EQ.0 .OR. IFAIL.EQ.-13 .OR.
     *       (IFAIL.GT.0 .AND. MOD(IFAIL/10,10).NE.0)) THEN
C           Noisy exit
            CALL X04AAF(0,NERR)
            DO 20 I = 1, NREC
               CALL X04BAF(NERR,REC(I))
   20       CONTINUE
            IF (IFAIL.NE.-13) THEN
               WRITE (MESS,FMT=99999) SRNAME, IERROR
               CALL X04BAF(NERR,MESS)
               IF (ABS(MOD(IFAIL,10)).NE.1) THEN
C                 Hard failure
                  CALL X04BAF(NERR,
     *                     ' ** NAG hard failure - execution terminated'
     *                        )
                  CALL P01ABZ
               ELSE
C                 Soft failure
                  CALL X04BAF(NERR,
     *                        ' ** NAG soft failure - control returned')
               END IF
            END IF
         END IF
      END IF
      P01ABF = IERROR
      RETURN
C
99999 FORMAT (' ** ABNORMAL EXIT from NAG Library routine ',A,': IFAIL',
     *  ' =',I6)
      END
      SUBROUTINE P01ABZ
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     Terminates execution when a hard failure occurs.
C
C     ******************** IMPLEMENTATION NOTE ********************
C     The following STOP statement may be replaced by a call to an
C     implementation-dependent routine to display a message and/or
C     to abort the program.
C     *************************************************************
C     .. Executable Statements ..
      STOP
      END
      SUBROUTINE X04AAF(I,NERR)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 7C REVISED IER-190 (MAY 1979)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-829 (DEC 1989).
C     IF I = 0, SETS NERR TO CURRENT ERROR MESSAGE UNIT NUMBER
C     (STORED IN NERR1).
C     IF I = 1, CHANGES CURRENT ERROR MESSAGE UNIT NUMBER TO
C     VALUE SPECIFIED BY NERR.
C
C     .. Scalar Arguments ..
      INTEGER           I, NERR
C     .. Local Scalars ..
      INTEGER           NERR1
C     .. Save statement ..
      SAVE              NERR1
C     .. Data statements ..
      DATA              NERR1/6/
C     .. Executable Statements ..
      IF (I.EQ.0) NERR = NERR1
      IF (I.EQ.1) NERR1 = NERR
      RETURN
      END
      SUBROUTINE X04BAF(NOUT,REC)
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     X04BAF writes the contents of REC to the unit defined by NOUT.
C
C     Trailing blanks are not output, except that if REC is entirely
C     blank, a single blank character is output.
C     If NOUT.lt.0, i.e. if NOUT is not a valid Fortran unit identifier,
C     then no output occurs.
C
C     .. Scalar Arguments ..
      INTEGER           NOUT
      CHARACTER*(*)     REC
C     .. Local Scalars ..
      INTEGER           I
C     .. Intrinsic Functions ..
      INTRINSIC         LEN
C     .. Executable Statements ..
      IF (NOUT.GE.0) THEN
C        Remove trailing blanks
         DO 20 I = LEN(REC), 2, -1
            IF (REC(I:I).NE.' ') GO TO 40
   20    CONTINUE
C        Write record to external file
   40    WRITE (NOUT,FMT=99999) REC(1:I)
      END IF
      RETURN
C
99999 FORMAT (A)
      END
      SUBROUTINE E02AEF(NPLUS1,A,XCAP,P,IFAIL)
C     NAG LIBRARY SUBROUTINE  E02AEF
C
C     E02AEF  EVALUATES A POLYNOMIAL FROM ITS CHEBYSHEV-
C     SERIES REPRESENTATION.
C
C     CLENSHAW METHOD WITH MODIFICATIONS DUE TO REINSCH
C     AND GENTLEMAN.
C
C     USES NAG LIBRARY ROUTINES  P01ABF  AND  X02AJF.
C     USES INTRINSIC FUNCTION  ABS.
C
C     STARTED - 1973.
C     COMPLETED - 1976.
C     AUTHOR - MGC AND JGH.
C
C     NAG COPYRIGHT 1975
C     MARK 5 RELEASE
C     MARK 7 REVISED IER-140 (DEC 1978)
C     MARK 9 REVISED. IER-352 (SEP 1981)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E02AEF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  P, XCAP
      INTEGER           IFAIL, NPLUS1
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NPLUS1)
C     .. Local Scalars ..
      DOUBLE PRECISION  BK, BKP1, BKP2, DK, ETA, FACTOR
      INTEGER           IERROR, K, KREV, N, NPLUS2
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
      IERROR = 0
      ETA = X02AJF()
C     INSERT CALL TO X02AJF
C
C     ETA  IS THE SMALLEST POSITIVE NUMBER SUCH THAT
C     THE COMPUTED VALUE OF  1.0 + ETA  EXCEEDS UNITY.
C
      IF (NPLUS1.GE.1) GO TO 20
      IERROR = 2
      GO TO 180
   20 IF (ABS(XCAP).LE.1.0D0+4.0D0*ETA) GO TO 40
      IERROR = 1
      P = 0.0D0
      GO TO 180
   40 IF (NPLUS1.GT.1) GO TO 60
      P = 0.5D0*A(1)
      GO TO 180
   60 N = NPLUS1 - 1
      NPLUS2 = N + 2
      K = NPLUS2
      IF (XCAP.GT.0.5D0) GO TO 140
      IF (XCAP.GE.-0.5D0) GO TO 100
C
C     GENTLEMAN*S MODIFIED RECURRENCE.
C
      FACTOR = 2.0D0*(1.0D0+XCAP)
      DK = 0.0D0
      BK = 0.0D0
      DO 80 KREV = 1, N
         K = K - 1
         DK = A(K) - DK + FACTOR*BK
         BK = DK - BK
   80 CONTINUE
      P = 0.5D0*A(1) - DK + 0.5D0*FACTOR*BK
      GO TO 180
C
C     CLENSHAW*S ORIGINAL RECURRENCE.
C
  100 FACTOR = 2.0D0*XCAP
      BKP1 = 0.0D0
      BK = 0.0D0
      DO 120 KREV = 1, N
         K = K - 1
         BKP2 = BKP1
         BKP1 = BK
         BK = A(K) - BKP2 + FACTOR*BKP1
  120 CONTINUE
      P = 0.5D0*A(1) - BKP1 + 0.5D0*FACTOR*BK
      GO TO 180
C
C     REINSCH*S MODIFIED RECURRENCE.
C
  140 FACTOR = 2.0D0*(1.0D0-XCAP)
      DK = 0.0D0
      BK = 0.0D0
      DO 160 KREV = 1, N
         K = K - 1
         DK = A(K) + DK - FACTOR*BK
         BK = BK + DK
  160 CONTINUE
      P = 0.5D0*A(1) + DK - 0.5D0*FACTOR*BK
  180 IF (IERROR) 200, 220, 200
  200 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
  220 IFAIL = 0
      RETURN
      END
      DOUBLE PRECISION FUNCTION X02AJF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS  (1/2)*B**(1-P)  IF ROUNDS IS .TRUE.
C     RETURNS  B**(1-P)  OTHERWISE
C
      DOUBLE PRECISION X02CON
      DATA X02CON /Z'3CA0000000000001' /
C     .. Executable Statements ..
      X02AJF = X02CON
      RETURN
      END
      DOUBLE PRECISION FUNCTION X02AKF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS  B**(EMIN-1)  (THE SMALLEST POSITIVE MODEL NUMBER)
C
      DOUBLE PRECISION X02CON
      DATA X02CON /Z'0010000000000000' /
C     .. Executable Statements ..
      X02AKF = X02CON
      RETURN
      END
      DOUBLE PRECISION FUNCTION S14ABF(X,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C        LNGAMMA(X) FUNCTION
C        ABRAMOWITZ AND STEGUN  CH.6
C
C     **************************************************************
C
C     TO EXTRACT THE CORRECT CODE FOR A PARTICULAR MACHINE-RANGE,
C     ACTIVATE THE STATEMENTS CONTAINED IN COMMENTS BEGINNING  CDD ,
C     WHERE  DD  IS THE APPROXIMATE NUMBER OF SIGNIFICANT DECIMAL
C     DIGITS REPRESENTED BY THE MACHINE
C     DELETE THE ILLEGAL DUMMY STATEMENTS OF THE FORM
C     * EXPANSION (NNNN) *
C
C     ALSO INSERT APPROPRIATE DATA STATEMENTS TO DEFINE CONSTANTS
C     WHICH DEPEND ON THE RANGE OF NUMBERS REPRESENTED BY THE
C     MACHINE, RATHER THAN THE PRECISION (SUITABLE STATEMENTS FOR
C     SOME MACHINES ARE CONTAINED IN COMMENTS BEGINNING CRD WHERE
C     D IS A DIGIT WHICH SIMPLY DISTINGUISHES A GROUP OF MACHINES).
C     DELETE THE ILLEGAL DUMMY DATA STATEMENTS WITH VALUES WRITTEN
C     *VALUE*
C
C     **************************************************************
C
C        IMPLEMENTATION DEPENDENT CONSTANTS
C
C        IF(X.LT.XSMALL)GAMMA(X)=1/X
C             I.E.   XSMALL*EULGAM.LE.XRELPR
C        LNGAM(XVBIG)=GBIG.LE.XOVFLO
C        LNR2PI=LN(SQRT(2*PI))
C        IF(X.GT.XBIG)LNGAM(X)=(X-0.5)LN(X)-X+LNR2PI
C
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='S14ABF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
      INTEGER                          IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 G, GBIG, LNR2PI, T, XBIG, XSMALL,
     *                                 XVBIG, Y
      INTEGER                          I, M
C     .. Local Arrays ..
      CHARACTER*1                      P01REC(1)
C     .. External Functions ..
      INTEGER                          P01ABF
      EXTERNAL                         P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC                        LOG, DBLE
C     .. Data statements ..
C08   DATA XSMALL,XBIG,LNR2PI/
C08  *1.0D-8,1.2D+3,9.18938533D-1/
C09   DATA XSMALL,XBIG,LNR2PI/
C09  *1.0D-9,4.8D+3,9.189385332D-1/
C12   DATA XSMALL,XBIG,LNR2PI/
C12  *1.0D-12,3.7D+5,9.189385332047D-1/
C15   DATA XSMALL,XBIG,LNR2PI/
C15  *1.0D-15,6.8D+6,9.189385332046727D-1/
      DATA XSMALL,XBIG,LNR2PI/
     *1.0D-17,7.7D+7,9.18938533204672742D-1/
C19   DATA XSMALL,XBIG,LNR2PI/
C19  *1.0D-19,3.1D+8,9.189385332046727418D-1/
C
C     RANGE DEPENDENT CONSTANTS
      DATA XVBIG,GBIG/2.55D+305,1.79D+308/
C     FOR IEEE SINGLE PRECISION
CR0   DATA XVBIG,GBIG/4.08E+36,3.40E+38/
C     FOR IBM 360/370 AND SIMILAR MACHINES
CR1   DATA XVBIG,GBIG/4.29D+73,7.231D+75/
C     FOR DEC10, HONEYWELL, UNIVAC 1100 (S.P.)
CR2   DATA XVBIG,GBIG/2.05D36,1.69D38/
C     FOR ICL 1900
CR3   DATA XVBIG,GBIG/3.39D+74,5.784D+76/
C     FOR CDC 7600/CYBER
CR4   DATA XVBIG,GBIG/1.72D+319,1.26D+322/
C     FOR UNIVAC 1100 (D.P.)
CR5   DATA XVBIG,GBIG/1.28D305,8.98D+307/
C     FOR IEEE DOUBLE PRECISION
CR7   DATA XVBIG,GBIG/2.54D+305,1.79D+308/
C     .. Executable Statements ..
      IF (X.GT.XSMALL) GO TO 20
C        VERY SMALL RANGE
      IF (X.LE.0.0D0) GO TO 160
      IFAIL = 0
      S14ABF = -LOG(X)
      GO TO 200
C
   20 IF (X.GT.15.0D0) GO TO 120
C        MAIN SMALL X RANGE
      M = X
      T = X - DBLE(M)
      M = M - 1
      G = 1.0D0
      IF (M) 40, 100, 60
   40 G = G/X
      GO TO 100
   60 DO 80 I = 1, M
         G = (X-DBLE(I))*G
   80 CONTINUE
  100 T = 2.0D0*T - 1.0D0
C
C      * EXPANSION (0026) *
C
C     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 08E.09
C08   Y = (((((((((((+1.88278283D-6*T-5.48272091D-6)*T+1.03144033D-5)
C08  *    *T-3.13088821D-5)*T+1.01593694D-4)*T-2.98340924D-4)
C08  *    *T+9.15547391D-4)*T-2.42216251D-3)*T+9.04037536D-3)
C08  *    *T-1.34119055D-2)*T+1.03703361D-1)*T+1.61692007D-2)*T +
C08  *    8.86226925D-1
C
C     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 09E.10
C09   Y = ((((((((((((-6.463247484D-7*T+1.882782826D-6)
C09  *    *T-3.382165478D-6)*T+1.031440334D-5)*T-3.393457634D-5)
C09  *    *T+1.015936944D-4)*T-2.967655076D-4)*T+9.155473906D-4)
C09  *    *T-2.422622002D-3)*T+9.040375355D-3)*T-1.341184808D-2)
C09  *    *T+1.037033609D-1)*T+1.616919866D-2)*T + 8.862269255D-1
C
C     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 12E.13
C12   Y = ((((((((((((((((-8.965837291520D-9*T+2.612707393536D-8)
C12  *    *T-3.802866827264D-8)*T+1.173294768947D-7)
C12  *    *T-4.275076254106D-7)*T+1.276176602829D-6)
C12  *    *T-3.748495971011D-6)*T+1.123829871408D-5)
C12  *    *T-3.364018663166D-5)*T+1.009331480887D-4)
C12  *    *T-2.968895120407D-4)*T+9.157850115110D-4)
C12  *    *T-2.422595461409D-3)*T+9.040335037321D-3)
C12  *    *T-1.341185056618D-2)*T+1.037033634184D-1)
C12  *    *T+1.616919872437D-2)*T + 8.862269254528D-1
C
C     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 15E.16
C15   Y = (((((((((((((((-1.243191705600000D-10*T+
C15  *    3.622882508800000D-10)*T-4.030909644800000D-10)
C15  *    *T+1.265236705280000D-9)*T-5.419466096640000D-9)
C15  *    *T+1.613133578240000D-8)*T-4.620920340480000D-8)
C15  *    *T+1.387603440435200D-7)*T-4.179652784537600D-7)
C15  *    *T+1.253148247777280D-6)*T-3.754930502328320D-6)
C15  *    *T+1.125234962812416D-5)*T-3.363759801664768D-5)
C15  *    *T+1.009281733953869D-4)*T-2.968901194293069D-4)
C15  *    *T+9.157859942174304D-4)*T-2.422595384546340D-3
C15   Y = ((((Y*T+9.040334940477911D-3)*T-1.341185057058971D-2)
C15  *    *T+1.037033634220705D-1)*T+1.616919872444243D-2)*T +
C15  *    8.862269254527580D-1
C
C     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 17E.18
      Y = (((((((((((((((-1.46381209600000000D-11*T+
     *    4.26560716800000000D-11)*T-4.01499750400000000D-11)
     *    *T+1.27679856640000000D-10)*T-6.13513953280000000D-10)
     *    *T+1.82243164160000000D-9)*T-5.11961333760000000D-9)
     *    *T+1.53835215257600000D-8)*T-4.64774927155200000D-8)
     *    *T+1.39383522590720000D-7)*T-4.17808776355840000D-7)
     *    *T+1.25281466396672000D-6)*T-3.75499034136576000D-6)
     *    *T+1.12524642975590400D-5)*T-3.36375833240268800D-5)
     *    *T+1.00928148823365120D-4)*T-2.96890121633200000D-4
      Y = ((((((Y*T+9.15785997288933120D-4)*T-2.42259538436268176D-3)
     *    *T+9.04033494028101968D-3)*T-1.34118505705967765D-2)
     *    *T+1.03703363422075456D-1)*T+1.61691987244425092D-2)*T +
     *    8.86226925452758013D-1
C
C     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 19E.19
C19   Y = (((((((((((((((+6.710886400000000000D-13*T-
C19  *    1.677721600000000000D-12)*T+6.710886400000000000D-13)
C19  *    *T-4.152360960000000000D-12)*T+2.499805184000000000D-11)
C19  *    *T-6.898581504000000000D-11)*T+1.859597107200000000D-10)
C19  *    *T-5.676387532800000000D-10)*T+1.725556326400000000D-9)
C19  *    *T-5.166307737600000000D-9)*T+1.548131827712000000D-8)
C19  *    *T-4.644574052352000000D-8)*T+1.393195837030400000D-7)
C19  *    *T-4.178233990758400000D-7)*T+1.252842254950400000D-6)
C19  *    *T-3.754985815285760000D-6)*T+1.125245651030528000D-5
C19   Y = (((((((((Y*T-3.363758423922688000D-5)
C19  *    *T+1.009281502108083200D-4)
C19  *    *T-2.968901215188000000D-4)*T+9.157859971435078400D-4)
C19  *    *T-2.422595384370689760D-3)*T+9.040334940288877920D-3)
C19  *    *T-1.341185057059651648D-2)*T+1.037033634220752902D-1)
C19  *    *T+1.616919872444250674D-2)*T + 8.862269254527580137D-1
C
      S14ABF = LOG(Y*G)
      IFAIL = 0
      GO TO 200
C
  120 IF (X.GT.XBIG) GO TO 140
C        MAIN LARGE X RANGE
      T = 450.0D0/(X*X) - 1.0D0
C
C      * EXPANSION (0059) *
C
C     EXPANSION (0059) EVALUATED AS Y(T)  --PRECISION 08E.09
C08   Y = (+3.89980902D-9*T-6.16502533D-6)*T + 8.33271644D-2
C
C     EXPANSION (0059) EVALUATED AS Y(T)  --PRECISION 09E.10
C09   Y = (+3.899809019D-9*T-6.165025333D-6)*T + 8.332716441D-2
C
C     EXPANSION (0059) EVALUATED AS Y(T)  --PRECISION 12E.13
C12   Y = ((-6.451144077930D-12*T+3.899809018958D-9)
C12  *    *T-6.165020494506D-6)*T + 8.332716440658D-2
C
C     EXPANSION (0059) EVALUATED AS Y(T)  --PRECISION 15E.16
C15   Y = (((+2.002019273379824D-14*T-6.451144077929628D-12)
C15  *    *T+3.899788998764847D-9)*T-6.165020494506090D-6)*T +
C15  *    8.332716440657866D-2
C
C     EXPANSION (0059) EVALUATED AS Y(T)  --PRECISION 17E.18
      Y = ((((-9.94561064728159347D-17*T+2.00201927337982364D-14)
     *    *T-6.45101975779653651D-12)*T+3.89978899876484712D-9)
     *    *T-6.16502049453716986D-6)*T + 8.33271644065786580D-2
C
C     EXPANSION (0059) EVALUATED AS Y(T)  --PRECISION 19E.19
C19   Y = (((((+7.196406678180202240D-19*T-9.945610647281593472D-17)
C19  *    *T+2.001911327279650935D-14)*T-6.451019757796536510D-12)
C19  *    *T+3.899788999169644998D-9)*T-6.165020494537169862D-6)*T +
C19  *    8.332716440657865795D-2
C
      S14ABF = (X-0.5D0)*LOG(X) - X + LNR2PI + Y/X
      IFAIL = 0
      GO TO 200
C
  140 IF (X.GT.XVBIG) GO TO 180
C        ASYMPTOTIC LARGE X RANGE
      S14ABF = (X-0.5D0)*LOG(X) - X + LNR2PI
      IFAIL = 0
      GO TO 200
C
C        FAILURE EXITS
  160 IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      S14ABF = 0.0D0
      GO TO 200
  180 IFAIL = P01ABF(IFAIL,2,SRNAME,0,P01REC)
      S14ABF = GBIG
C
  200 RETURN
C
      END
      DOUBLE PRECISION FUNCTION X02AMF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS THE 'SAFE RANGE' PARAMETER
C     I.E. THE SMALLEST POSITIVE MODEL NUMBER Z SUCH THAT
C     FOR ANY X WHICH SATISFIES X.GE.Z AND X.LE.1/Z
C     THE FOLLOWING CAN BE COMPUTED WITHOUT OVERFLOW, UNDERFLOW OR OTHER
C     ERROR
C
C        -X
C        1.0/X
C        SQRT(X)
C        LOG(X)
C        EXP(LOG(X))
C        Y**(LOG(X)/LOG(Y)) FOR ANY Y
C
      DOUBLE PRECISION X02CON
      DATA X02CON /Z'0010000000000000' /
C     .. Executable Statements ..
      X02AMF = X02CON
      RETURN
      END
      SUBROUTINE S14BAF(A,X,TOL,P,Q,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, HALF, QUARTR, THRTWO, TWO, THREE,
     *                  FOUR, EIGHT, TWENTY, PSEVEN, ONEP4
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0,HALF=0.5D0,QUARTR=0.25D0,
     *                  THRTWO=1.5D0,TWO=2.0D0,THREE=3.0D0,FOUR=4.0D0,
     *                  EIGHT=8.0D0,TWENTY=20.0D0,PSEVEN=0.7D0,
     *                  ONEP4=1.4D0)
      INTEGER           MAXIT, NTERMS
      PARAMETER         (MAXIT=600,NTERMS=22)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='S14BAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, P, Q, TOL, X
      INTEGER           IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION  ALG, ALGP1, ALGS, ALPHA, ALX, AP1, BOT, EPS,
     *                  EPS1, GA, PP, QQ, RHO, RR, SAFE, SS, SUM, TERM,
     *                  TT, U, UNDFL, V, XMA, XPA, Y
      INTEGER           IERR, K, NREC, TIFAIL
C     .. Local Arrays ..
      DOUBLE PRECISION  C(NTERMS)
      CHARACTER*80      REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  S14AAF, S14ABF, S14BAZ, X02AJF, X02AMF
      INTEGER           P01ABF
      EXTERNAL          S14AAF, S14ABF, S14BAZ, X02AJF, X02AMF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, EXP, LOG, MAX, SQRT
C     .. Data statements ..
      DATA              C/5.7721566490153286060651209008D-1,
     *                  -6.5587807152025388107701951515D-1,
     *                  -4.200263503409523552900393488D-2,
     *                  1.665386113822914895017007951D-1,
     *                  -4.21977345555443367482083013D-2,
     *                  -9.6219715278769735621149217D-3,
     *                  7.21894324666309954230950103D-3,
     *                  -1.16516759185906511211971D-3,
     *                  -2.15241674114950972815730D-4,
     *                  1.2805028238811618615320D-4,
     *                  -2.013485478078823865569D-5,
     *                  -1.2504934821426706573D-6,
     *                  1.1330272319816958824D-6,
     *                  -2.056338416977607103D-7, 6.1160951044814158D-9,
     *                  5.0020076444692229D-9, -1.181274570487020D-9,
     *                  1.04342671169110D-10, 7.782263439905D-12,
     *                  -3.696805618642D-12, 5.1003702875D-13,
     *                  -2.058326054D-14/
C     .. Executable Statements ..
C
C     Let  GAMMA(A)  denote the gamma function and  GAM(A,X)  the
C     (complementary) incomplete gamma function,
C
C     GAM(A,X)=integral from T=X to T=infinity of EXP(-T)*T**(A-1).
C
C     Let  GAMSTAR(A,X)  denote Tricomi's form of the incomplete gamma
C     function, which for A.gt.0. is defined by
C
C     GAMSTAR(A,X)=(X**(-A)/GAMMA(A))*integral from T=0 to T=X of
C                EXP(-T)*T**(A-1).
C
C     For the purpose of this subroutine, these functions are normalized
C     as follows:
C
C     P(A,X)  =   (X**A)*GAMSTAR(A,X).
C
C     Q(A,X)  =   GAM(A,X)/GAMMA(A),
C
C     The program below attempts to evaluate  P(A,X)  and  Q(A,X),
C     both to an accuracy of TOL significant decimal digits, for
C     positive A and nonnegative X. There are (rare) instances in
C     which the accuracy attained is somewhat less than the accuracy
C     specified. The discrepancy, however, should never exceed one or
C     two (decimal) orders of accuracy.
C     The functions are evaluated by means of the Maclaurin expansion,
C     Taylor expansion, or Legendre's continued fraction, depending on
C     the values of A and X. However, when A .ge. 20.0 and
C     0.7*A .le. X .le. 1.4*A, the asymptotic expansion of Temme is
C     used for greater efficiency.
C
C     Arguments:
C         A - The first argument of P and Q.
C         X - The second argument of P and Q.
C       TOL - The number of correct significant decimal digits
C             desired in the results.
C         P - An output variable returning the value of P(A,X).
C         Q - An output variable returning the value of Q(A,X).
C     IFAIL - Error flag. The values of IFAIL have these meanings:
C            0 - Successful exit.
C            1 - Illegal negative or zero argument A. The routine
C                exits with the value zero for P and Q.
C            2 - Illegal negative argument X. The routine exits
C                with the value zero for P and Q.
C            3 - Convergence fails within MAXIT (=600) iterations,
C                either in Taylor's series or in Legendre's continued
C                fraction. Reason unknown. The computation is
C                aborted and the routine exits with the value zero
C                for P and Q.
C
C     The data declaration contains the successive coefficients in the
C     Maclaurin expansion of (1/GAMMA(A+1))-1.
C     Values of these coefficients (to 31 decimal places) can be found
C     in Table 5 of J.W.Wrench,jr.:
C     Concerning Two Series for the Gamma Function  , Math. Comput.
C     22, 1968, 617-626.
C
C     This routine is derived from ACM Algorithm 542.
C
C     References:
C     W. Gautschi,   'A Computational Procedure for Incomplete Gamma
C     Functions', ACM Trans. Math. Software, Vol. 5, No. 4, 1979,
C     466-481.
C     N. M. Temme,   'On the computation of the incomplete gamma
C     functions for large values of the parameters', proceedings of
C     Shrivenham conference 'Algorithms for Approximation',
C     ed. J.C.Mason and M.G.Cox, Oxford University Press 1987.
C
      P = ZERO
      Q = ZERO
      NREC = 0
      IERR = 0
      IF (A.LE.ZERO) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99999) A
      ELSE IF (X.LT.ZERO) THEN
         IERR = 2
         NREC = 1
         WRITE (REC,FMT=99998) X
      ELSE IF (X.EQ.ZERO) THEN
         P = ZERO
         Q = ONE
      ELSE
C        1.0E-18 is the accuracy of constants in this routine and
C        auxiliaries.
         EPS = MAX(X02AJF(),1.0D-18)
         IF (TOL.GT.EPS .AND. TOL.LE.ONE) EPS = TOL
         SAFE = X02AMF()
         IF (A.GE.TWENTY .AND. PSEVEN*A.LE.X .AND. X/ONEP4.LE.A) THEN
C           Use the asymptotic expansion of Temme.
            UNDFL = LOG(SAFE)
            IF (A.LE.X) THEN
               Q = S14BAZ(A,X,.FALSE.,EPS,UNDFL)
               P = ONE - Q
            ELSE
               P = S14BAZ(A,X,.TRUE.,EPS,UNDFL)
               Q = ONE - P
            END IF
         ELSE
            ALX = LOG(X)
            IF (X.LT.QUARTR) THEN
               ALPHA = LOG(HALF)/ALX
            ELSE
               ALPHA = X + QUARTR
            END IF
            BOT = LOG(SAFE)
            EPS1 = EPS/100
            AP1 = A + ONE
C
C           Evaluation of the logarithm of GAMMA(A+1.0).
C
            TIFAIL = 1
            ALGP1 = S14ABF(AP1,TIFAIL)
            IF (TIFAIL.EQ.2) THEN
C              ln gamma(a+1) overflows. P and Q are 1.0 and 0.0,
C              or vice versa, to machine precision.
               IF (A.GT.X) THEN
                  P = ZERO
                  Q = ONE
               ELSE
                  P = ONE
                  Q = ZERO
               END IF
            ELSE IF (A.GT.ALPHA) THEN
C
C              Evaluation of P(A,X) for A.gt.ALPHA(X) by Taylor
C              expansion.
C
               TERM = ONE
               SUM = ONE
               K = 0
   20          K = K + 1
               IF (K.GT.MAXIT) GO TO 120
               TERM = X*TERM/(A+K)
               SUM = SUM + TERM
               IF (ABS(TERM).GT.EPS*SUM) GO TO 20
               ALGS = A*ALX - X + LOG(SUM) - ALGP1
               IF (ALGS.LE.BOT) THEN
                  P = ZERO
               ELSE
                  P = EXP(ALGS)
               END IF
               Q = ONE - P
            ELSE IF (X.GT.THRTWO) THEN
C
C              Evaluation of Q(A,X) for X.gt.1.5 and A.le.ALPHA(X) by
C              means of the Legendre continued fraction.
C
               XPA = X + ONE - A
               XMA = X - ONE - A
               IF (XPA.GT.ONE/SQRT(SAFE)) THEN
C                 XPA*XMA would overflow, but P = 1.0 to machine
C                 precision.
                  P = ONE
                  Q = ZERO
               ELSE
                  PP = ZERO
                  QQ = XPA*XMA
                  RR = FOUR*XPA
                  SS = -A + ONE
                  TERM = ONE
                  SUM = ONE
                  RHO = ZERO
                  K = 1
   40             K = K + 1
                  IF (K.GT.MAXIT) GO TO 120
                  PP = PP + SS
                  QQ = QQ + RR
                  RR = RR + EIGHT
                  SS = SS + TWO
                  TT = PP*(ONE+RHO)
                  RHO = TT/(QQ-TT)
                  TERM = RHO*TERM
                  SUM = SUM + TERM
                  IF (ABS(TERM).GT.EPS*SUM) GO TO 40
                  ALG = A*ALX - X + LOG(A*SUM/XPA) - ALGP1
                  IF (ALG.LE.BOT) THEN
                     Q = ZERO
                  ELSE
                     Q = EXP(ALG)
                  END IF
                  P = ONE - Q
               END IF
            ELSE
C
C              Direct evaluation of Q(A,X) and P(A,X) for X.le.1.5
C              and A.le.ALPHA(X).
C
               IF (A.GT.HALF) THEN
                  TIFAIL = -1
                  U = S14AAF(A,TIFAIL) - (X**A)/A
               ELSE
C                 NTERMS = 22 is sufficient for approximately 18
C                 decimal place accuracy in the summation below.
                  SUM = C(NTERMS)
                  DO 60 K = NTERMS - 1, 1, -1
                     SUM = A*SUM + C(K)
   60             CONTINUE
                  GA = -SUM/(ONE+A*SUM)
                  Y = A*ALX
                  SUM = ONE
                  TERM = ONE
                  K = 1
   80             K = K + 1
                  IF (K.GT.MAXIT) GO TO 120
                  TERM = Y*TERM/K
                  SUM = SUM + TERM
                  IF (ABS(TERM).GT.EPS1*SUM) GO TO 80
                  U = GA - SUM*ALX
               END IF
               PP = A*X
               QQ = AP1
               RR = A + THREE
               TERM = ONE
               SUM = ONE
               K = 1
  100          K = K + 1
               IF (K.GT.MAXIT) GO TO 120
               PP = PP + X
               QQ = QQ + RR
               RR = RR + TWO
               TERM = -PP*TERM/QQ
               SUM = SUM + TERM
               IF (ABS(TERM).GT.EPS1*SUM) GO TO 100
               V = (X**AP1)*SUM/AP1
               Q = U + V
               Q = A*Q*EXP(-ALGP1)
               P = ONE - Q
            END IF
         END IF
      END IF
      GO TO 140
  120 IERR = 3
      NREC = 1
      WRITE (REC,FMT=99997) MAXIT
  140 IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, A.le.0.0 : A =',1P,D13.5)
99998 FORMAT (1X,'** On entry, X.lt.0.0 : X =',1P,D13.5)
99997 FORMAT (1X,'** Algorithm fails to terminate in ',I4,' iterations')
      END
      SUBROUTINE G01AEF(N,K2,X,ICLASS,CINT,IFREQ,XMIN,XMAX,IFAIL)
C     MARK 4 RELEASE NAG COPYRIGHT 1974.
C     MARK 4.5 REVISED IER-49/42
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12A REVISED. IER-510 (AUG 1986).
C
C     G01AEF - THE ROUTINE CONSTRUCTS A FREQUENCY DISTRIBUTION
C     FOR VALUES IN AN ARRAY X WITH EITHER ROUTINE
C     CALCULATED OR USER SUPPLIED CLASS INTERVALS
C
C     USES NAG LIBRARY ROUTINE P01AAF
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G01AEF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  XMAX, XMIN
      INTEGER           ICLASS, IFAIL, K2, N
C     .. Array Arguments ..
      DOUBLE PRECISION  CINT(K2), X(N)
      INTEGER           IFREQ(K2)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHA, CR, STEP, XTEMP
      INTEGER           I, IERR, J, JJ, K0, LDN, LUP, NDN, NOB, NUP
      LOGICAL           LK0
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE
C     .. Executable Statements ..
      IF (K2-2) 480, 20, 20
   20 IF (N-1) 500, 40, 40
C     ZERO  FREQUENCIES
   40 DO 60 I = 1, K2
         IFREQ(I) = 0
   60 CONTINUE
      NUP = 0
      NDN = 0
      K0 = K2 - 2
      LK0 = K0 .EQ. 0
C     DETERMINE  MIN  AND  MAX
      XMIN = X(1)
      XMAX = XMIN
      IF (N.EQ.1) GO TO 100
      DO 80 I = 2, N
         XTEMP = X(I)
         IF (XTEMP.LT.XMIN) XMIN = XTEMP
         IF (XTEMP.GT.XMAX) XMAX = XTEMP
   80 CONTINUE
  100 IF (ICLASS.EQ.1) GO TO 200
      IF (LK0) GO TO 420
C     CALCULATE CLASS INTERVALS WITH EQUAL SPACING
      CR = XMAX - XMIN
      ALPHA = 0.001D0
      IF (CR.GT.0.0D0) GO TO 160
      IF (XMIN) 120, 140, 120
  120 CR = ABS(XMIN)*ALPHA
      GO TO 160
  140 CR = ALPHA
  160 STEP = CR*(1.0D0+ALPHA)/DBLE(K0)
      CR = XMIN - 0.5D0*ALPHA*CR
      CINT(1) = CR
      DO 180 I = 1, K0
         XTEMP = CR
         CR = XTEMP + STEP
         CINT(I+1) = CR
  180 CONTINUE
      GO TO 260
C     CLASS INTERVALS SUPPLIED--CHECK IN ASCENDING ORDER
  200 CR = CINT(1)
      IF (LK0) GO TO 440
      JJ = K0 + 1
      DO 240 I = 2, JJ
         XTEMP = CINT(I)
         IF (CR.LT.XTEMP) GO TO 220
         IERR = 3
         GO TO 520
  220    CR = XTEMP
  240 CONTINUE
C     DETERMINE CLASS TO WHICH EACH CASE BELONGS
  260 NOB = K2/2
      CR = CINT(NOB)
      LDN = NOB + 1
      LUP = K2 - 1
      K0 = NOB - 1
      LK0 = K0 .EQ. 0
      DO 380 I = 1, N
         XTEMP = X(I)
         IF (XTEMP-CR) 320, 280, 280
C        VARIATE VALUE ABOVE MIDDLE  INTERVAL
  280    DO 300 J = LDN, LUP
            IF (XTEMP.GE.CINT(J)) GO TO 300
            IFREQ(J) = IFREQ(J) + 1
            GO TO 380
  300    CONTINUE
         NUP = NUP + 1
         GO TO 380
C        VARIATE VALUE BELOW MIDDLE  INTERVAL
  320    IF (LK0) GO TO 360
         DO 340 J = 1, K0
            JJ = NOB - J
            IF (XTEMP.LT.CINT(JJ)) GO TO 340
            JJ = JJ + 1
            IFREQ(JJ) = IFREQ(JJ) + 1
            GO TO 380
  340    CONTINUE
  360    NDN = NDN + 1
  380 CONTINUE
  400 IFREQ(1) = NDN
      IFREQ(K2) = NUP
      IFAIL = 0
      RETURN
  420 CR = (XMAX+XMIN)*0.5D0
      CINT(1) = CR
  440 DO 460 I = 1, N
         IF (X(I).LT.CR) NDN = NDN + 1
  460 CONTINUE
      NUP = N - NDN
      GO TO 400
  480 IERR = 1
      GO TO 520
  500 IERR = 2
  520 IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,P01REC)
      RETURN
      END
      SUBROUTINE G08CGF(K2,IFREQ,CINT,DIST,PAR,IPARAM,PROB,CHISQ,P,NDF,
     *                  EVAL,CHISQI,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G08CGF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CHISQ, P
      INTEGER           IFAIL, IPARAM, K2, NDF
      CHARACTER         DIST
C     .. Array Arguments ..
      DOUBLE PRECISION  CHISQI(K2), CINT(K2-1), EVAL(K2), PAR(2),
     *                  PROB(K2)
      INTEGER           IFREQ(K2)
C     .. Local Scalars ..
      DOUBLE PRECISION  DF, STORP, STORP1, SUM, ULEN, XN
      INTEGER           I, IERR, IF2, IFAULT, N, NREC, NSMALL
      LOGICAL           CHI, EXPON, GAMMA, NORM, UNIF, USERDF
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  G01ECF, G08CGZ, X02AJF
      INTEGER           P01ABF
      EXTERNAL          G01ECF, G08CGZ, X02AJF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE
C     .. Executable Statements ..
C
      CHISQ = 0.0D0
      IERR = 0
      NREC = 1
      N = 0
      NORM = .FALSE.
      UNIF = .FALSE.
      EXPON = .FALSE.
      CHI = .FALSE.
      GAMMA = .FALSE.
      USERDF = .FALSE.
C
C     Check input parameters.
C
      IF (K2.LT.2) THEN
         IERR = 1
         WRITE (P01REC,FMT=99999) K2
      ELSE IF (DIST.NE.'N' .AND. DIST.NE.'n' .AND. DIST.NE.'U' .AND.
     *         DIST.NE.'u' .AND. DIST.NE.'E' .AND. DIST.NE.'e' .AND.
     *         DIST.NE.'C' .AND. DIST.NE.'c' .AND. DIST.NE.'G' .AND.
     *         DIST.NE.'g' .AND. DIST.NE.'A' .AND. DIST.NE.'a') THEN
         IERR = 2
         WRITE (P01REC,FMT=99998) DIST
      ELSE IF (IPARAM.LT.0 .OR. IPARAM.GE.K2-1) THEN
         NREC = 1
         IERR = 3
         WRITE (P01REC,FMT=99997) IPARAM
C
C        Check for elements of IFREQ.lt.0, and that CINT is in ascending
C          order, (and calculate N).
C
      ELSE
C
         IF (DIST.EQ.'N' .OR. DIST.EQ.'n') NORM = .TRUE.
         IF (DIST.EQ.'U' .OR. DIST.EQ.'u') UNIF = .TRUE.
         IF (DIST.EQ.'E' .OR. DIST.EQ.'e') EXPON = .TRUE.
         IF (DIST.EQ.'C' .OR. DIST.EQ.'c') CHI = .TRUE.
         IF (DIST.EQ.'G' .OR. DIST.EQ.'g') GAMMA = .TRUE.
         IF (DIST.EQ.'A' .OR. DIST.EQ.'a') USERDF = .TRUE.
C
         IF (IFREQ(1).LT.0) THEN
            IERR = 4
            WRITE (P01REC,FMT=99996)
            GO TO 140
         END IF
         N = N + IFREQ(1)
         DO 20 I = 2, K2 - 1
            N = N + IFREQ(I)
            IF (IFREQ(I).LT.0) THEN
               IERR = 4
               WRITE (P01REC,FMT=99996)
               GO TO 140
            END IF
            IF (CINT(I).LE.CINT(I-1)) THEN
               IERR = 5
               WRITE (P01REC,FMT=99995)
               GO TO 140
            END IF
   20    CONTINUE
         IF (IFREQ(K2).LT.0) THEN
            IERR = 4
            WRITE (P01REC,FMT=99996)
            GO TO 140
         END IF
         N = N + IFREQ(K2)
         IF ((EXPON .OR. CHI .OR. GAMMA) .AND. CINT(1).LT.0.0D0) THEN
            IERR = 6
            WRITE (P01REC,FMT=99994)
            GO TO 140
         END IF
C
C        Check the parameters of the CDF's
C
         NREC = 2
         IF (NORM .AND. PAR(2).LE.0.0D0) THEN
            IERR = 7
            WRITE (P01REC,FMT=99993) PAR(2)
         ELSE IF (UNIF) THEN
            IF (PAR(1).GE.PAR(2) .OR. PAR(1).GT.CINT(1) .OR. PAR(2)
     *          .LT.CINT(K2-1)) THEN
               IERR = 7
               WRITE (P01REC,FMT=99992) PAR(1), PAR(2)
            END IF
         ELSE IF (EXPON .AND. PAR(1).LE.0.0D0) THEN
            IERR = 7
            WRITE (P01REC,FMT=99991) PAR(1)
         ELSE IF (CHI .AND. PAR(1).LE.0.0D0) THEN
            IERR = 7
            WRITE (P01REC,FMT=99990) PAR(1)
         ELSE IF (GAMMA .AND. (PAR(1).LE.0.0D0 .OR. PAR(2).LE.0.0D0))
     *            THEN
            IERR = 7
            WRITE (P01REC,FMT=99989) PAR(1), PAR(2)
         ELSE IF (USERDF) THEN
            SUM = 0.0D0
            DO 40 I = 1, K2
               IF (PROB(I).LE.0.0D0) THEN
                  IERR = 8
                  WRITE (P01REC,FMT=99988) I, PROB(I)
                  GO TO 140
               END IF
               SUM = SUM + PROB(I)
   40       CONTINUE
            IF (ABS(SUM-1.0D0).GT.X02AJF()) THEN
               IERR = 8
               WRITE (P01REC,FMT=99987) SUM
               GO TO 140
            END IF
         END IF
      END IF
      IF (IERR.EQ.0) THEN
         NREC = 0
         XN = DBLE(N)
C
C        Determine which distribution is to be used and calculate
C          probabilities, expected frequencies.
C
         IF (USERDF) THEN
C
C           User supplied probabilities.
C
            DO 60 I = 1, K2
               EVAL(I) = PROB(I)*XN
   60       CONTINUE
         ELSE IF (UNIF) THEN
C
C           Uniform distribution.
C
            ULEN = PAR(2) - PAR(1)
            EVAL(1) = XN*(CINT(1)-PAR(1))/ULEN
            DO 80 I = 2, K2 - 1
               EVAL(I) = XN*(CINT(I)-CINT(I-1))/ULEN
   80       CONTINUE
            EVAL(K2) = XN*(PAR(2)-CINT(K2-1))/ULEN
         ELSE
C
C           One of the other standard distributions.
C
            IFAULT = 0
            STORP1 = 0.0D0
            DO 100 I = 1, K2 - 1
               STORP = G08CGZ(DIST,CINT(I),PAR,IFAULT)
               IF (IFAULT.NE.0) THEN
                  NREC = 2
                  IERR = 11
                  WRITE (P01REC,FMT=99986)
               END IF
               EVAL(I) = (STORP-STORP1)*XN
               STORP1 = STORP
  100       CONTINUE
            STORP = 1.0D0 - STORP
            EVAL(K2) = STORP*XN
         END IF
C
C        Calculate contributions to the statistic and hence
C        the statistic.
C
         NSMALL = 0
         DO 120 I = 1, K2
            IF (EVAL(I).LT.1.0D0) NSMALL = NSMALL + 1
            IF (EVAL(I).NE.0.0D0) THEN
               CHISQI(I) = (DBLE(IFREQ(I))-EVAL(I))*(DBLE(IFREQ(I))
     *                     -EVAL(I))/EVAL(I)
            ELSE IF (EVAL(I).EQ.0.0D0 .AND. IFREQ(I).EQ.0) THEN
               CHISQI(I) = 0.0D0
            ELSE
               IERR = 9
               NREC = 1
               WRITE (P01REC,FMT=99985)
               GO TO 140
            END IF
            CHISQ = CHISQ + CHISQI(I)
  120    CONTINUE
         IF (NSMALL.NE.0) THEN
            IERR = 10
            NREC = 1
            WRITE (P01REC,FMT=99984) NSMALL
         END IF
C
C        Calculate the degrees of freedom associated with the test.
C        No. intervals - 1 - No. estimated parameters.
C
         NDF = K2 - 1 - IPARAM
         DF = DBLE(NDF)
C
C        Determine the probability associated with the statistic.
C
         IF2 = 1
         P = G01ECF('UPPER',CHISQ,DF,IF2)
      END IF
  140 IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, NCLASS.lt.2 : NCLASS = ',I16)
99998 FORMAT (1X,'** On entry, DIST is not valid : DIST = ',A1)
99997 FORMAT (1X,'** On entry, NPEST.lt.0 or NPEST.ge.NCLASS-1 : NPEST',
     *       ' = ',I16)
99996 FORMAT (1X,'** On entry, at least one class has frequency less t',
     *       'han zero.')
99995 FORMAT (1X,'** On entry, the elements of CINT are not in ascendi',
     *       'ng order.')
99994 FORMAT (1X,'** On entry, the intervals contained in CINT are inv',
     *       'alid.')
99993 FORMAT (1X,'** On entry, the variance of the normal distribution',
     *       ',PAR(2) is invalid:',/3X,'PAR(2) = ',1P,D13.5)
99992 FORMAT (1X,'** On entry, the parameters of the uniform  distribu',
     *       'tion are invalid:',/3X,'PAR(1) = ',1P,D13.5,' PAR(2) = ',
     *       D13.5)
99991 FORMAT (1X,'** On entry, the parameter of the exponential distri',
     *       'bution is invalid:',/4X,'PAR(1) = ',1P,D13.5)
99990 FORMAT (1X,'** On entry, the parameter of the Chi-squared distri',
     *       'bution is invalid:',/4X,'PAR(1) = ',1P,D13.5)
99989 FORMAT (1X,'** On entry, the parameter(s) of the gamma distribut',
     *       'ion are invalid: ',/4X,'PAR(1) = ',1P,D13.5,' PAR(2) = ',
     *       D13.5)
99988 FORMAT (1X,'** On entry, with DIST.eq.''A'' or ''a'' at least on',
     *       'e elements of PROB.le.0.0',/4X,'For the ',I16,' th class',
     *       ' PROB(i) = ',D13.5)
99987 FORMAT (1X,'** On entry with DIST.eq.''A'' or ''a'' the sum of t',
     *       'he elements of PROB',/4X,'does .ne. 1.0.  SUM of PROB(i)',
     *       '  = ',D13.5)
99986 FORMAT (1X,'** The solution has failed to converge whilst comput',
     *       'ing expected values for the',/4X,'gamma orchi-squared di',
     *       'st. The result returned should be an adequate approx.')
99985 FORMAT (1X,'** An expected frequency equals zero, when the obser',
     *       'ved frequency was not.')
99984 FORMAT (1X,'** ',I16,' classes have expected frequency less than',
     *       ' one.')
      END
      DOUBLE PRECISION FUNCTION G01ECF(TAIL,X,DF,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Computes the lower tail probability of X for a
C     CHI-squared distribution with DF degrees of freedom.
C
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G01ECF')
      DOUBLE PRECISION                 ZERO, HALF
      PARAMETER                        (ZERO=0.0D0,HALF=0.5D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 DF, X
      INTEGER                          IFAIL
      CHARACTER*1                      TAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 P, Q, TOL
      INTEGER                          IERR, IFA
C     .. Local Arrays ..
      CHARACTER*80                     REC(1)
C     .. External Functions ..
      INTEGER                          P01ABF
      EXTERNAL                         P01ABF
C     .. External Subroutines ..
      EXTERNAL                         S14BAF
C     .. Executable Statements ..
      G01ECF = ZERO
      IERR = 0
      IF (X.LT.ZERO) THEN
         IERR = 2
         WRITE (REC,FMT=99999) X
      ELSE IF (DF.LE.ZERO) THEN
         IERR = 3
         WRITE (REC,FMT=99998) DF
      ELSE IF (TAIL.NE.'L' .AND. TAIL.NE.'U' .AND. TAIL.NE.'l' .AND.
     *         TAIL.NE.'u') THEN
         IERR = 1
         WRITE (REC,FMT=99997) TAIL
      END IF
      IF (IERR.EQ.0) THEN
C
C           Use transformation of a GAMMA.
C
         IFA = 1
         TOL = 0.5D-5
         CALL S14BAF(HALF*DF,HALF*X,TOL,P,Q,IFA)
         IF (TAIL.EQ.'L' .OR. TAIL.EQ.'l') THEN
            G01ECF = P
         ELSE
            G01ECF = Q
         END IF
         IF (IFA.EQ.3) THEN
            IERR = 4
            WRITE (REC,FMT=99996)
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,1,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, X.lt.0.0: X = ',1P,D13.5)
99998 FORMAT (1X,'** On entry, DF.le.0.0: DF = ',1P,D13.5)
99997 FORMAT (1X,'** On entry, TAIL is not valid: TAIL = ',A1)
99996 FORMAT (1X,'** Full accuracy has not been achieved.')
      END
      DOUBLE PRECISION FUNCTION G08CGZ(CDIST,RINT,PAR,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Provides the probabilities from the appropriate CDF.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 RINT
      INTEGER                          IFAIL
      CHARACTER                        CDIST
C     .. Array Arguments ..
      DOUBLE PRECISION                 PAR(2)
C     .. Local Scalars ..
      DOUBLE PRECISION                 P, PEXP, Q, SR, TOL, Z
      INTEGER                          IFA
C     .. External Functions ..
      DOUBLE PRECISION                 S15ABF, X02AJF, X02AMF
      EXTERNAL                         S15ABF, X02AJF, X02AMF
C     .. External Subroutines ..
      EXTERNAL                         S14BAF
C     .. Intrinsic Functions ..
      INTRINSIC                        EXP, LOG, SQRT
C     .. Executable Statements ..
C
      G08CGZ = 0.0D0
      IF (CDIST.EQ.'N' .OR. CDIST.EQ.'n') THEN
         Z = (RINT-PAR(1))/SQRT(PAR(2))
         IFA = 1
         G08CGZ = S15ABF(Z,IFA)
      ELSE IF (CDIST.EQ.'E' .OR. CDIST.EQ.'e') THEN
         SR = LOG(X02AMF())
         PEXP = PAR(1)*RINT
         IF (-PEXP.GE.SR) THEN
            G08CGZ = 1.0D0 - EXP(-PEXP)
         ELSE
            G08CGZ = 1.0D0
         END IF
      ELSE IF (CDIST.EQ.'C' .OR. CDIST.EQ.'c') THEN
         TOL = X02AJF()*10.0D0
         IFA = 1
         CALL S14BAF(PAR(1)/2.0D0,RINT/2.0D0,TOL,P,Q,IFA)
         G08CGZ = P
      ELSE IF (CDIST.EQ.'G' .OR. CDIST.EQ.'g') THEN
         Z = RINT/PAR(2)
         TOL = X02AJF()*10.0D0
         IFA = 1
         CALL S14BAF(PAR(1),Z,TOL,P,Q,IFA)
         G08CGZ = P
      END IF
      IF (IFA.EQ.3) IFAIL = 1
      RETURN
      END
      DOUBLE PRECISION FUNCTION X01AAF(X)
C     MARK 8 RE-ISSUE. NAG COPYRIGHT 1980.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     RETURNS THE VALUE OF THE MATHEMATICAL CONSTANT PI.
C
C     X IS A DUMMY ARGUMENT
C
C     IT MAY BE NECESSARY TO ROUND THE REAL CONSTANT IN THE
C     ASSIGNMENT STATEMENT TO A SMALLER NUMBER OF SIGNIFICANT
C     DIGITS IN ORDER TO AVOID COMPILATION PROBLEMS.  IF SO, THEN
C     THE NUMBER OF DIGITS RETAINED SHOULD NOT BE LESS THAN
C     .     2 + INT(FLOAT(IT)*ALOG10(IB))
C     WHERE  IB  IS THE BASE FOR THE REPRESENTATION OF FLOATING-
C     .             POINT NUMBERS
C     . AND  IT  IS THE NUMBER OF IB-ARY DIGITS IN THE MANTISSA OF
C     .             A FLOATING-POINT NUMBER.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
C     .. Executable Statements ..
      X01AAF = 3.14159265358979323846264338328D0
      RETURN
      END
      DOUBLE PRECISION FUNCTION S14BAZ(A,X,GETP,EPS,UNDFL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Parameters ..
      DOUBLE PRECISION                 TWO, ONE, HALF, ZERO
      PARAMETER                        (TWO=2.0D0,ONE=1.0D0,HALF=0.5D0,
     *                                 ZERO=0.0D0)
      DOUBLE PRECISION                 RT2PI
      PARAMETER                        (RT2PI=2.5066282746310005024D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 A, EPS, UNDFL, X
      LOGICAL                          GETP
C     .. Local Scalars ..
      DOUBLE PRECISION                 DIF, ETA, U, V, Y
      INTEGER                          IFAIL, S
C     .. External Functions ..
      DOUBLE PRECISION                 S14BAX, S14BAY, S15ADF
      EXTERNAL                         S14BAX, S14BAY, S15ADF
C     .. Intrinsic Functions ..
      INTRINSIC                        EXP, SQRT
C     .. Executable Statements ..
C
      IF (GETP) THEN
         S = -ONE
      ELSE
         S = ONE
      END IF
      DIF = (X-A)/A
      Y = A*S14BAX(DIF)
      IF (Y.LT.ZERO) Y = ZERO
      ETA = SQRT(TWO*Y/A)
      V = SQRT(Y)
      IF (X.LT.A) THEN
         ETA = -ETA
         V = -V
      END IF
      IFAIL = 0
      U = HALF*S15ADF(S*V,IFAIL)
      IF (-Y.GE.UNDFL) THEN
         V = S*EXP(-Y)*S14BAY(ETA,A,EPS)/(RT2PI*SQRT(A))
      ELSE
C        exp(-Y) underflows.
         V = ZERO
      END IF
      S14BAZ = U + V
      RETURN
      END
      DOUBLE PRECISION FUNCTION S14AAF(X,IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 7C REVISED IER-184 (MAY 1979)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     GAMMA FUNCTION
C
C     **************************************************************
C
C     TO EXTRACT THE CORRECT CODE FOR A PARTICULAR MACHINE-RANGE,
C     ACTIVATE THE STATEMENTS CONTAINED IN COMMENTS BEGINNING  CDD ,
C     WHERE  DD  IS THE APPROXIMATE NUMBER OF SIGNIFICANT DECIMAL
C     DIGITS REPRESENTED BY THE MACHINE
C     DELETE THE ILLEGAL DUMMY STATEMENTS OF THE FORM
C     * EXPANSION (NNNN) *
C
C     ALSO INSERT APPROPRIATE DATA STATEMENTS TO DEFINE CONSTANTS
C     WHICH DEPEND ON THE RANGE OF NUMBERS REPRESENTED BY THE
C     MACHINE, RATHER THAN THE PRECISION (SUITABLE STATEMENTS FOR
C     SOME MACHINES ARE CONTAINED IN COMMENTS BEGINNING CRD WHERE
C     D IS A DIGIT WHICH SIMPLY DISTINGUISHES A GROUP OF MACHINES).
C     DELETE THE ILLEGAL DUMMY DATA STATEMENTS WITH VALUES WRITTEN
C     *VALUE*
C
C     **************************************************************
C
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='S14AAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
      INTEGER                          IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 G, GBIG, T, XBIG, XMINV, XSMALL,
     *                                 Y
      INTEGER                          I, M
C     .. Local Arrays ..
      CHARACTER*1                      P01REC(1)
C     .. External Functions ..
      INTEGER                          P01ABF
      EXTERNAL                         P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, SIGN, DBLE
C     .. Data statements ..
C08   DATA XSMALL/1.0D-8/
C09   DATA XSMALL/3.0D-9/
C12   DATA XSMALL/1.0D-12/
C15   DATA XSMALL/3.0D-15/
      DATA XSMALL/1.0D-17/
C19   DATA XSMALL/1.7D-18/
C
      DATA XBIG,GBIG,XMINV/ 1.70D+2,4.3D+304,2.23D-308 /
C     XBIG = LARGEST X SUCH THAT  GAMMA(X) .LT. MAXREAL
C                            AND  1.0/GAMMA(X+1.0) .GT. MINREAL
C             (ROUNDED DOWN TO AN INTEGER)
C     GBIG = GAMMA(XBIG)
C     XMINV = MAX(1.0/MAXREAL,MINREAL)  (ROUNDED UP)
C     FOR IEEE SINGLE PRECISION
CR0   DATA XBIG,GBIG,XMINV /33.0E0,2.6E+35,1.2E-38/
C     FOR IBM 360/370 AND SIMILAR MACHINES
CR1   DATA XBIG,GBIG,XMINV /57.0D0,7.1D+74,1.4D-76/
C     FOR DEC-10, HONEYWELL, UNIVAC 1100 (S.P.)
CR2   DATA XBIG,GBIG,XMINV /34.0D0,8.7D+36,5.9D-39/
C     FOR ICL 1900
CR3   DATA XBIG,GBIG,XMINV /58.0D0,4.0D+76,1.8D-77/
C     FOR CDC 7600/CYBER
CR4   DATA XBIG,GBIG,XMINV /164.0D0,2.0D+291,3.2D-294/
C     FOR UNIVAC 1100 (D.P.)
CR5   DATA XBIG,GBIG,XMINV /171.0D0,7.3D+306,1.2D-308/
C     FOR IEEE DOUBLE PRECISION
CR7   DATA XBIG,GBIG,XMINV /170.0D0,4.3D+304,2.3D-308/
C     .. Executable Statements ..
C
C     ERROR 1 AND 2 TEST
      T = ABS(X)
      IF (T.GT.XBIG) GO TO 160
C     SMALL RANGE TEST
      IF (T.LE.XSMALL) GO TO 140
C     MAIN RANGE REDUCTION
      M = X
      IF (X.LT.0.0D0) GO TO 80
      T = X - DBLE(M)
      M = M - 1
      G = 1.0D0
      IF (M) 20, 120, 40
   20 G = G/X
      GO TO 120
   40 DO 60 I = 1, M
         G = (X-DBLE(I))*G
   60 CONTINUE
      GO TO 120
   80 T = X - DBLE(M-1)
C     ERROR 4 TEST
      IF (T.EQ.1.0D0) GO TO 220
      M = 1 - M
      G = X
      DO 100 I = 1, M
         G = (DBLE(I)+X)*G
  100 CONTINUE
      G = 1.0D0/G
  120 T = 2.0D0*T - 1.0D0
C
C      * EXPANSION (0026) *
C
C     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 08E.09
C08   Y = ((((((((((((+1.88278283D-6)*T-5.48272091D-6)
C08  *    *T+1.03144033D-5)*T-3.13088821D-5)*T+1.01593694D-4)
C08  *    *T-2.98340924D-4)*T+9.15547391D-4)*T-2.42216251D-3)
C08  *    *T+9.04037536D-3)*T-1.34119055D-2)*T+1.03703361D-1)
C08  *    *T+1.61692007D-2)*T + 8.86226925D-1
C
C     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 09E.10
C09   Y = (((((((((((((-6.463247484D-7)*T+1.882782826D-6)
C09  *    *T-3.382165478D-6)*T+1.031440334D-5)*T-3.393457634D-5)
C09  *    *T+1.015936944D-4)*T-2.967655076D-4)*T+9.155473906D-4)
C09  *    *T-2.422622002D-3)*T+9.040375355D-3)*T-1.341184808D-2)
C09  *    *T+1.037033609D-1)*T+1.616919866D-2)*T + 8.862269255D-1
C
C     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 12E.13
C12   Y = (((((((((((((((-7.613347676160D-8)*T+2.218377726362D-7)
C12  *    *T-3.608242105549D-7)*T+1.106350622249D-6)
C12  *    *T-3.810416284805D-6)*T+1.138199762073D-5)
C12  *    *T-3.360744031186D-5)*T+1.008657892262D-4)
C12  *    *T-2.968993359366D-4)*T+9.158021574033D-4)
C12  *    *T-2.422593898516D-3)*T+9.040332894085D-3)
C12  *    *T-1.341185067782D-2)*T+1.037033635205D-1)
C12  *    *T+1.616919872669D-2)*T + 8.862269254520D-1
C
C     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 15E.16
C15   Y = (((((((((((((((-1.243191705600000D-10
C15  *    *T+3.622882508800000D-10)*T-4.030909644800000D-10)
C15  *    *T+1.265236705280000D-9)*T-5.419466096640000D-9)
C15  *    *T+1.613133578240000D-8)*T-4.620920340480000D-8)
C15  *    *T+1.387603440435200D-7)*T-4.179652784537600D-7)
C15  *    *T+1.253148247777280D-6)*T-3.754930502328320D-6)
C15  *    *T+1.125234962812416D-5)*T-3.363759801664768D-5)
C15  *    *T+1.009281733953869D-4)*T-2.968901194293069D-4)
C15  *    *T+9.157859942174304D-4)*T-2.422595384546340D-3
C15   Y = ((((Y*T+9.040334940477911D-3)*T-1.341185057058971D-2)
C15  *    *T+1.037033634220705D-1)*T+1.616919872444243D-2)*T +
C15  *     8.862269254527580D-1
C
C     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 17E.18
      Y = (((((((((((((((-1.46381209600000000D-11
     *    *T+4.26560716800000000D-11)*T-4.01499750400000000D-11)
     *    *T+1.27679856640000000D-10)*T-6.13513953280000000D-10)
     *    *T+1.82243164160000000D-9)*T-5.11961333760000000D-9)
     *    *T+1.53835215257600000D-8)*T-4.64774927155200000D-8)
     *    *T+1.39383522590720000D-7)*T-4.17808776355840000D-7)
     *    *T+1.25281466396672000D-6)*T-3.75499034136576000D-6)
     *    *T+1.12524642975590400D-5)*T-3.36375833240268800D-5)
     *    *T+1.00928148823365120D-4)*T-2.96890121633200000D-4
      Y = ((((((Y*T+9.15785997288933120D-4)*T-2.42259538436268176D-3)
     *    *T+9.04033494028101968D-3)*T-1.34118505705967765D-2)
     *    *T+1.03703363422075456D-1)*T+1.61691987244425092D-2)*T +
     *     8.86226925452758013D-1
C
C     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 19E.20
C19   Y = (((((((((((((((+6.7108864000000000000D-13
C19  *    *T-1.6777216000000000000D-12)*T+6.7108864000000000000D-13)
C19  *    *T-4.1523609600000000000D-12)*T+2.4998051840000000000D-11)
C19  *    *T-6.8985815040000000000D-11)*T+1.8595971072000000000D-10)
C19  *    *T-5.6763875328000000000D-10)*T+1.7255563264000000000D-9)
C19  *    *T-5.1663077376000000000D-9)*T+1.5481318277120000000D-8)
C19  *    *T-4.6445740523520000000D-8)*T+1.3931958370304000000D-7)
C19  *    *T-4.1782339907584000000D-7)*T+1.2528422549504000000D-6)
C19  *    *T-3.7549858152857600000D-6)*T+1.1252456510305280000D-5
C19   Y = (((((((((Y*T-3.3637584239226880000D-5)
C19  *    *T+1.0092815021080832000D-4)
C19  *    *T-2.9689012151880000000D-4)*T+9.1578599714350784000D-4)
C19  *    *T-2.4225953843706897600D-3)*T+9.0403349402888779200D-3)
C19  *    *T-1.3411850570596516480D-2)*T+1.0370336342207529018D-1)
C19  *    *T+1.6169198724442506740D-2)*T + 8.8622692545275801366D-1
C
      S14AAF = Y*G
      IFAIL = 0
      GO TO 240
C
C     ERROR 3 TEST
  140 IF (T.LT.XMINV) GO TO 200
      S14AAF = 1.0D0/X
      IFAIL = 0
      GO TO 240
C
C     ERROR EXITS
  160 IF (X.LT.0.0D0) GO TO 180
      IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      S14AAF = GBIG
      GO TO 240
C
  180 IFAIL = P01ABF(IFAIL,2,SRNAME,0,P01REC)
      S14AAF = 0.0D0
      GO TO 240
C
  200 IFAIL = P01ABF(IFAIL,3,SRNAME,0,P01REC)
      T = X
      IF (X.EQ.0.0D0) T = 1.0D0
      S14AAF = SIGN(1.0D0/XMINV,T)
      GO TO 240
C
  220 IFAIL = P01ABF(IFAIL,4,SRNAME,0,P01REC)
      S14AAF = GBIG
C
  240 RETURN
      END
      DOUBLE PRECISION FUNCTION S15ABF(X,IFAIL)
C     MARK 5A REVISED - NAG COPYRIGHT 1976
C     MARK 5C REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-753 (DEC 1989).
C     CUMULATIVE NORMAL DISTRIBUTION P(X)
C     ******************************************************************
C     TO SELECT THE CORRECT VALUE FOR A PARTICULAR MACHINE-RANGE,
C     ACTIVATE THE STATEMENT CONTAINED IN A COMMENT BEGINNING  CDD ,
C     WHERE  DD  IS THE APPROXIMATE NUMBER OF SIGNIFICANT DECIMAL
C     DIGITS REPRESENTED BY THE MACHINE.
C     DELETE THE ILLEGAL DUMMY STATEMENT OF THE FORM
C     * EXPANSION (DATA) *
C     ******************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
      INTEGER                          IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 RRTWO
C     .. External Functions ..
      DOUBLE PRECISION                 S15ADF
      EXTERNAL                         S15ADF
C     .. Data statements ..
C     PRECISION DEPENDENT CONSTANTS
C      * EXPANSION (DATA) *
C08   DATA RRTWO/7.0710678D-1/
C12   DATA RRTWO/7.07106781187D-1/
C14   DATA RRTWO/7.0710678118655D-1/
      DATA RRTWO/7.071067811865475D-1/
C18   DATA RRTWO/7.07106781186547524D-1/
C     .. Executable Statements ..
C
C     NO FAILURE EXITS - IFAIL SET TO ZERO IN S15ADF
      S15ABF = 0.5D0*S15ADF(-X*RRTWO,IFAIL)
      RETURN
C
      END
      DOUBLE PRECISION FUNCTION S15ADF(X,IFAIL)
C     MARK 5A REVISED - NAG COPYRIGHT 1976
C     MARK 5C REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     COMPLEMENT OF ERROR FUNCTION ERFC(X)
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
      INTEGER                          IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 T, XHI, XLO, Y
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, EXP
C     .. Data statements ..
C     PRECISION DEPENDENT CONSTANTS
C08   DATA XLO/-4.5D0/
C12   DATA XLO/-5.25D0/
C14   DATA XLO/-5.75D0/
      DATA XLO/-6.25D0/
C18   DATA XLO/-6.5D0/
C
C     RANGE DEPENDENT CONSTANTS
      DATA XHI/ 2.66D+1 /
C     XHI = LARGEST X SUCH THAT EXP(-X*X) .GT. MINREAL (ROUNDED DOWN)
CR1   DATA XHI/13.0D0/
CR2   DATA XHI/9.5D0/
CR3   DATA XHI/13.0D0/
CR4   DATA XHI/25.0D0/
CR5   DATA XHI/26.0D0/
C     .. Executable Statements ..
C
C     NO FAILURE EXITS
      IFAIL = 0
C     TEST EXTREME EXITS
      IF (X.GE.XHI) GO TO 20
      IF (X.LE.XLO) GO TO 40
C
C     EXPANSION ARGUMENT
      T = 1.0D0 - 7.5D0/(ABS(X)+3.75D0)
C
C      * EXPANSION (0021) *
C
C     EXPANSION (0021) EVALUATED AS Y(T)  --PRECISION 08E
C08   Y = (((((((((((+3.1475326D-5)*T-1.3874589D-4)*T-6.4127909D-6)
C08  *    *T+1.7866301D-3)*T-8.2316935D-3)*T+2.4151896D-2)
C08  *    *T-5.4799165D-2)*T+1.0260225D-1)*T-1.6357229D-1)
C08  *    *T+2.2600824D-1)*T-2.7342192D-1)*T + 1.4558972D-1
C
C     EXPANSION (0021) EVALUATED AS Y(T)  --PRECISION 12E
C12   Y = ((((((((((((((((-4.21661579602D-8*T-8.63384346353D-8)
C12  *    *T+6.06038693567D-7)*T+5.90655413508D-7)
C12  *    *T-6.12872971594D-6)*T+3.73223486059D-6)
C12  *    *T+4.78645837248D-5)*T-1.52546487034D-4)
C12  *    *T-2.55222360474D-5)*T+1.80299061562D-3)
C12  *    *T-8.22062412199D-3)*T+2.41432185990D-2)
C12  *    *T-5.48023263289D-2)*T+1.02604312548D-1)
C12  *    *T-1.63571895545D-1)*T+2.26008066898D-1)
C12  *    *T-2.73421931495D-1)*T + 1.45589721275D-1
C
C     EXPANSION (0021) EVALUATED AS Y(T)  --PRECISION 14E
C14   Y = (((((((((((((((-2.2356173494379D-9
C14  *    *T+4.5302502889845D-9)*T+2.5918103316137D-8)
C14  *    *T-6.3684846832869D-8)*T-1.7642194353331D-7)
C14  *    *T+6.4907607131235D-7)*T+7.4296952017617D-7)
C14  *    *T-6.1758018478516D-6)*T+3.5866167916231D-6)
C14  *    *T+4.7895180610590D-5)*T-1.5246364229106D-4)
C14  *    *T-2.5534256252531D-5)*T+1.8029626230333D-3)
C14  *    *T-8.2206213481002D-3)*T+2.4143223946968D-2)
C14  *    *T-5.4802326675661D-2)*T+1.0260431203382D-1
C14   Y = (((Y*T-1.6357189552481D-1)*T+2.2600806691658D-1)
C14  *    *T-2.7342193149541D-1)*T + 1.4558972127504D-1
C
C     EXPANSION (0021) EVALUATED AS Y(T)  --PRECISION 16E
      Y = (((((((((((((((+3.328130055126039D-10
     *    *T-5.718639670776992D-10)*T-4.066088879757269D-9)
     *    *T+7.532536116142436D-9)*T+3.026547320064576D-8)
     *    *T-7.043998994397452D-8)*T-1.822565715362025D-7)
     *    *T+6.575825478226343D-7)*T+7.478317101785790D-7)
     *    *T-6.182369348098529D-6)*T+3.584014089915968D-6)
     *    *T+4.789838226695987D-5)*T-1.524627476123466D-4)
     *    *T-2.553523453642242D-5)*T+1.802962431316418D-3)
     *    *T-8.220621168415435D-3)*T+2.414322397093253D-2
      Y = (((((Y*T-5.480232669380236D-2)*T+1.026043120322792D-1)
     *    *T-1.635718955239687D-1)*T+2.260080669166197D-1)
     *    *T-2.734219314954260D-1)*T + 1.455897212750385D-1
C
C     EXPANSION (0021) EVALUATED AS Y(T)  --PRECISION 18E
C18   Y = (((((((((((((((-1.58023488119651697D-11
C18  *    *T-4.94972069009392927D-11)*T+1.86424953544623784D-10)
C18  *    *T+6.29796246918239617D-10)*T-1.34751340973493898D-9)
C18  *    *T-4.84566988844706300D-9)*T+9.22474802259858004D-9)
C18  *    *T+3.14410318645430670D-8)*T-7.26754673242913196D-8)
C18  *    *T-1.83380699508554268D-7)*T+6.59488268069175234D-7)
C18  *    *T+7.48541685740064308D-7)*T-6.18344429012694168D-6)
C18  *    *T+3.58371497984145357D-6)*T+4.78987832434182054D-5)
C18  *    *T-1.52462664665855354D-4)*T-2.55353311432760448D-5
C18   Y = ((((((((Y*T+1.80296241673597993D-3)
C18  *    *T-8.22062115413991215D-3)
C18  *    *T+2.41432239724445769D-2)*T-5.48023266949776152D-2)
C18  *    *T+1.02604312032198239D-1)*T-1.63571895523923969D-1)
C18  *    *T+2.26008066916621431D-1)*T-2.73421931495426482D-1)*T +
C18  *     1.45589721275038539D-1
C
      S15ADF = EXP(-X*X)*Y
      IF (X.LT.0.0D0) S15ADF = 2.0D0 - S15ADF
      RETURN
C
   20 S15ADF = 0.0D0
      RETURN
   40 S15ADF = 2.0D0
      RETURN
C
      END
      DOUBLE PRECISION FUNCTION S14BAX(X)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     Computes the value X - log(1+X), retaining full relative
C     precision. Uses two Chebyshev expansions, the first valid in the
C     range -0.5 .lt. X .lt. 0.0, the second valid in the range
C     0.0 .lt. X .lt. 0.5.
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
C     .. Local Scalars ..
      DOUBLE PRECISION                 Z
C     .. Executable Statements ..
      IF (X.GT.0.0D0) THEN
         Z = 4*X - 1
C        Precision           20 sig. figs.
         S14BAX = ((((((((((((((-2.967444749673129063172D-15+
     *            (5.727889008167279942775D-16)*Z)
     *            *Z+1.254316797763179284544D-14)
     *            *Z+(-6.609208145248494083774D-14))
     *            *Z+3.552093631141621449360D-13)
     *            *Z+(-1.877575086584465734489D-12))
     *            *Z+9.950448268789644316317D-12)
     *            *Z+(-5.297680484461551503860D-11))
     *            *Z+2.832511945923442439052D-10)
     *            *Z+(-1.521770786533239703319D-09))
     *            *Z+8.221622611294801319744D-09)
     *            *Z+(-4.471056587610118468800D-08))
     *            *Z+2.450395094674988857885D-07)
     *            *Z+(-1.355590675105630409968D-06))
     *            *Z+7.586141840685986975605D-06)*Z
         S14BAX = X*X*((((((S14BAX+(-4.307383586344297282484D-05))
     *            *Z+2.492281965528721274994D-04)
     *            *Z+(-1.479382557242297784497D-03))
     *            *Z+9.109536917931723215386D-03)
     *            *Z+(-5.940635794528781547636D-02))
     *            *Z+4.297031789726439077393D-01)
      ELSE
         Z = 4*X + 1
C        Precision           20 sig. figs.
         S14BAX = ((((((((((((((-9.222006652207022772023D-14+
     *            (3.049209823713897715108D-14)*Z)
     *            *Z+8.110971610764310768061D-14)
     *            *Z+(-2.708761967476698009188D-13))
     *            *Z+1.468315445990726325182D-12)
     *            *Z+(-4.549612295389177687275D-12))
     *            *Z+1.318723069570494295267D-11)
     *            *Z+(-4.157371211079606768770D-11))
     *            *Z+1.323612522929628013599D-10)
     *            *Z+(-4.184113387833978015279D-10))
     *            *Z+1.325817212562882710866D-09)
     *            *Z+(-4.217964594365521564431D-09))
     *            *Z+1.346846941394504845741D-08)
     *            *Z+(-4.318289149822655693289D-08))
     *            *Z+1.391089264405774451800D-07)*Z
         S14BAX = X*X*((((((((((((S14BAX+(-4.505690978103018400967D-07))
     *            *Z+1.468654737538004800256D-06)
     *            *Z+(-4.823073082193340091105D-06))
     *            *Z+1.598133959833499951509D-05)
     *            *Z+(-5.353471603473257684307D-05))
     *            *Z+1.817808088831846419015D-04)
     *            *Z+(-6.280405138046454551685D-04))
     *            *Z+2.220117130128515512154D-03)
     *            *Z+(-8.100449505773729833165D-03))
     *            *Z+3.096169990770673931421D-02)
     *            *Z+(-1.275070148763436552823D-01))
     *            *Z+6.029131592284948390275D-01)
      END IF
      RETURN
      END
      DOUBLE PRECISION FUNCTION S14BAY(ETA,A,EPS)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Parameters ..
      DOUBLE PRECISION                 ONE
      PARAMETER                        (ONE=1.0D0)
      INTEGER                          NTERMS
      PARAMETER                        (NTERMS=26)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 A, EPS, ETA
C     .. Local Scalars ..
      DOUBLE PRECISION                 S, T, Y
      INTEGER                          I, M
C     .. Local Arrays ..
      DOUBLE PRECISION                 BM(0:NTERMS-1), FM(0:NTERMS)
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS
C     .. Data statements ..
      DATA                             (FM(I),I=0,18)
     *                                 /1.0000000000000000000D+00,
     *                                 -3.3333333333333333333D-01,
     *                                 8.3333333333333333333D-02,
     *                                 -1.4814814814814814815D-02,
     *                                 1.1574074074074074074D-03,
     *                                 3.5273368606701940035D-04,
     *                                 -1.7875514403292181070D-04,
     *                                 3.9192631785224377817D-05,
     *                                 -2.1854485106799921615D-06,
     *                                 -1.8540622107151599607D-06,
     *                                 8.2967113409530860050D-07,
     *                                 -1.7665952736826079304D-07,
     *                                 6.7078535434014985804D-09,
     *                                 1.0261809784240308043D-08,
     *                                 -4.3820360184533531866D-09,
     *                                 9.1476995822367902342D-10,
     *                                 -2.5514193994946249767D-11,
     *                                 -5.8307721325504250675D-11,
     *                                 2.4361948020667416244D-11/
      DATA                             (FM(I),I=19,26)
     *                                 /-5.0276692801141755891D-12,
     *                                 1.1004392031956134771D-13,
     *                                 3.3717632624009853788D-13,
     *                                 -1.3923887224181620659D-13,
     *                                 2.8534893807047443204D-14,
     *                                 -5.1391118342425726190D-16,
     *                                 -1.9752288294349442835D-15,
     *                                 8.0995211567045613341D-16/
C     .. Executable Statements ..
C     When A .ge. 20.0, NTERMS = 26 is sufficient for approximately
C     18 decimal place accuracy.
      BM(NTERMS-1) = FM(NTERMS)
      BM(NTERMS-2) = FM(NTERMS-1)
      DO 20 M = NTERMS - 1, 2, -1
         BM(M-2) = FM(M-1) + M*BM(M)/A
   20 CONTINUE
C
      S = BM(0)
      Y = ETA
      M = 1
C
   40 T = BM(M)*Y
      S = S + T
      M = M + 1
      Y = Y*ETA
      IF (ABS(T/S).GE.EPS .AND. M.LT.NTERMS) GO TO 40
C
      S14BAY = S/(ONE+BM(1)/A)
      RETURN
      END
      SUBROUTINE E02CBF(MFIRST,MLAST,K,L,X,XMIN,XMAX,Y,YMIN,YMAX,FF,A,
     *                  NA,WORK,NWORK,IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     THIS SUBROUTINE EVALUATES A POLYNOMIAL OF DEGREE K AND L
C     RESPECTIVELY IN THE INDEPENDENT VARIABLES X AND Y.  THE
C     POLYNOMIAL IS GIVEN IN DOUBLE CHEBYSHEV SERIES FORM
C     A(I,J) * TI(XCAP) * TJ(YCAP),
C     SUMMED OVER I = 0,1,...K AND J = 0,1,...L WITH THE CONVENTION
C     THAT TERMS WITH EITHER I OR J ZERO ARE HALVED AND THE TERM
C     WITH BOTH I AND J ZERO IS MULTIPLIED BY 0.25. HERE TI(XCAP)
C     IS THE CHEBYSHEV POLYNOMIAL OF THE FIRST KIND OF DEGREE I
C     WITH ARGUMENT XCAP=((X - XMIN) - (XMAX - X))/(XMAX - XMIN).
C     TJ(YCAP) IS DEFINED SIMILARLY. THE COEFFICIENT A(I,J)
C     SHOULD BE STORED IN ELEMENT (L + 1)*I + J + 1 OF THE SINGLE
C     DIMENSION ARRAY A. THE EVALUATION IS PERFORMED FOR A SINGLE
C     GIVEN VALUE OF Y WITH EACH X VALUE GIVEN IN X(R), FOR R =
C     MFIRST, MFIRST+1,....,MLAST.
C
C     STARTED - 1978.
C     COMPLETED - 1978.
C     AUTHOR - GTA.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E02CBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  XMAX, XMIN, Y, YMAX, YMIN
      INTEGER           IFAIL, K, L, MFIRST, MLAST, NA, NWORK
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NA), FF(MLAST), WORK(NWORK), X(MLAST)
C     .. Local Scalars ..
      DOUBLE PRECISION  D, XCAP, YCAP
      INTEGER           I, IERROR, KP1, LP1, M, R
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          E02AEF
C     .. Executable Statements ..
      KP1 = K + 1
      LP1 = L + 1
      M = MLAST - MFIRST + 1
C
C     CHECK THAT THE INTEGER INPUT PARAMETERS HAVE REASONABLE
C     VALUES
C
      IERROR = 1
      IF (M.LE.0 .OR. K.LT.0 .OR. L.LT.0 .OR. NA.LT.KP1*LP1 .OR.
     *    NWORK.LT.KP1) GO TO 80
C
C     CHECK THAT THE Y RANGE IS REASONABLE AND THAT THE GIVEN
C     VALUE OF Y IS NOT OUTSIDE IT
C
      IERROR = 2
      IF (YMIN.GE.YMAX .OR. Y.LT.YMIN .OR. Y.GT.YMAX) GO TO 80
      D = XMAX - XMIN
C
C     CHECK THAT THE X RANGE IS REASONABLE AND THAT NONE OF
C     THE GIVEN VALUES OF X IS OUTSIDE IT
C
      IERROR = 3
      IF (D.LE.0.0D+0) GO TO 80
      DO 20 R = MFIRST, MLAST
         IF (X(R).LT.XMIN .OR. X(R).GT.XMAX) GO TO 80
   20 CONTINUE
C
C     CALCULATE YCAP, THE NORMALIZED VALUE OF Y
C
      YCAP = ((Y-YMIN)-(YMAX-Y))/(YMAX-YMIN)
      IERROR = 1
      R = -L
C
C     EVALUATE THE COEFFICIENTS OF THE POLYNOMIAL FOR THE GIVEN Y
C
      DO 40 I = 1, KP1
         R = R + LP1
         CALL E02AEF(LP1,A(R),YCAP,WORK(I),IERROR)
         IERROR = IERROR + 1
         IF (IERROR.NE.1) GO TO 80
   40 CONTINUE
C
C     EVALUATE THE POLYNOMAL AT THE GIVEN X VALUES
C
      DO 60 R = MFIRST, MLAST
         XCAP = ((X(R)-XMIN)-(XMAX-X(R)))/D
         IERROR = 1
         CALL E02AEF(KP1,WORK,XCAP,FF(R),IERROR)
         IF (IERROR.EQ.0) GO TO 60
         IERROR = 3
         GO TO 80
   60 CONTINUE
   80 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
      SUBROUTINE E04ABF(FUN,EPS,T,A,B,MAXCAL,X,F,IFAIL)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 8 REVISED. IER-231 (MAR 1980).
C     MARK 8D REVISED. IER-272 (DEC 1980).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C
C     **************************************************************
C
C     E04ABF ATTEMPTS TO FIND A MINIMUM IN AN INTERVAL A .LE. X .LE.
C     B OF A FUNCTION F(X) OF THE SCALAR X, USING FUNCTION VALUES
C     ONLY.
C
C     IT IS BASED ON THE SUBROUTINE UNIFUN IN THE NPL ALGORITHMS
C     LIBRARY (REF. NO. E4/13/F). THE FUNCTION F(X) IS DEFINED BY
C     THE USER-SUPPLIED SUBROUTINE FUN. T AND EPS DEFINE A TOLERANCE
C     TOL = EPS * ABS(X) + T, AND FUN IS NEVER EVALUATED AT TWO
C     POINTS CLOSER THAN TOL. IF FUN IS DELTA-UNIMODAL, FOR SOME
C     DELTA LESS THAN TOL, THEN X APPROXIMATES THE GLOBAL MINIMUM OF
C     FUN WITH AN ERROR LESS THAN 3*TOL. IF FUN IS NOT DELTA-
C     UNIMODAL ON (A, B), THEN X MAY APPROXIMATE A LOCAL, BUT NON
C     GLOBAL, MINIMUM. EPS SHOULD BE NO SMALLER THAN 2*EPSMCH, AND
C     PREFERABLY NOT MUCH LESS THAN SQRT(EPSMCH), WHERE EPSMCH IS
C     THE RELATIVE MACHINE PRECISION. T SHOULD BE POSITIVE. NOTE
C     THAT, FOR CONSISTENCY WITH OTHER E04 DOCUMENTATION, THE NAME
C     FUNCT IS USED INSTEAD OF FUN IN THE WRITE-UP.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN, HAZEL M.
C     BARBER AND MARGARET H. WRIGHT, D.N.A.C., NATIONAL PHYSICAL
C     LABORATORY, ENGLAND
C
C     **************************************************************
C
C     FUN
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E04ABF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, B, EPS, F, T, X
      INTEGER           IFAIL, MAXCAL
C     .. Subroutine Arguments ..
      EXTERNAL          FUN
C     .. Local Scalars ..
      DOUBLE PRECISION  B1, D, E, EPSMCH, F1, F2, FA, FU, FV, FW,
     *                  GTEST1, GTEST2, GU, OLDF, PT2, PT4, PT6, RR,
     *                  RTEPS, SCXBD, SFTBND, SS, TOL, U, X1, X2,
     *                  XLAMDA, XV, XW
      INTEGER           IFLAG, ILOC, ISAVE, NUMF
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          E04ABZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      ISAVE = IFAIL
C
C     A MACHINE-DEPENDENT CONSTANT IS SET HERE. EPSMCH IS THE
C     SMALLEST POSITIVE REAL NUMBER SUCH THAT 1.0 + EPSMCH .GT. 1.0
C
      EPSMCH = X02AJF()
      RTEPS = SQRT(EPSMCH)
      IF (EPS.LT.EPSMCH) EPS = RTEPS
      IF (T.LT.EPSMCH) T = RTEPS
C
C     ERROR IN INPUT PARAMETERS
C
      IFAIL = 1
      IF (A+T.GE.B .OR. MAXCAL.LT.3) GO TO 140
      SFTBND = A
      PT2 = (B-A)*2.0D-1
      PT4 = PT2 + PT2
      PT6 = PT2 + PT4
      X1 = A + PT4
      CALL FUN(X1,F1)
      X2 = B - PT4
      CALL FUN(X2,F2)
      XLAMDA = B
      IF (F1.GT.F2) GO TO 20
      X = X1
      A = -PT4
      B = PT2
      XW = PT2
      B1 = PT2
      RR = 1.0D+0
      D = -PT2
      FW = F2
      FV = F1
      F = F1
C
C     SET STEP TO NEW POINT
C
      U = -PT2
      GO TO 40
   20 X = X2
      A = -PT2
      B = PT4 + EPS*ABS(XLAMDA) + T
      XW = -PT2
      B1 = B
      RR = -1.0D+0
      D = PT2
      FW = F1
      FV = F2
      F = F2
C
C     SET STEP TO NEW POINT
C
      U = PT2
   40 XV = 0.0D+0
      SCXBD = PT4
      E = PT6
      SS = 0.0D+0
      FA = FW + T
      OLDF = FA
      GTEST1 = 0.0D+0
      GTEST2 = 0.0D+0
      TOL = EPS*ABS(X) + T
      CALL FUN(X+U,FU)
      GU = 0.0D+0
      NUMF = 3
C
C     SET ILOC TO 3 SO THAT THE MAIN SECTION OF E04ABZ IS EXECUTED
C     AS THE INITIAL 3 POINTS HAVE ALREADY BEEN SET UP
C
      ILOC = 3
   60 CALL E04ABZ(EPS,T,0.0D+0,SFTBND,XLAMDA,U,FU,GU,X,F,XW,FW,XV,FV,A,
     *            FA,B,OLDF,B1,SCXBD,E,D,RR,SS,GTEST1,GTEST2,TOL,ILOC,
     *            IFLAG)
      IF (IFLAG.NE.1) GO TO 100
      IF (NUMF.GE.MAXCAL) GO TO 80
      CALL FUN(X+U,FU)
      NUMF = NUMF + 1
      GO TO 60
   80 IFAIL = 2
      GO TO 120
  100 IFAIL = 0
  120 MAXCAL = NUMF
      A = A + X
      B = B + X
  140 CONTINUE
      IF (IFAIL.EQ.0) RETURN
      IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
      RETURN
C
C     END OF E04ABF (UNIFUN)
C
      END
      SUBROUTINE E04ABZ(EPS,T,ETA,SFTBND,XLAMDA,U,FU,GU,XMIN,FMIN,XW,FW,
     *                  XV,FV,A,FA,B,OLDF,B1,SCXBD,E,D,RR,SS,GTEST1,
     *                  GTEST2,TOL,ILOC,ITEST)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 7 REISSUE
C     MARK 8 REVISED. IER-239 (APR 1980).
C     MARK 8 REVISED. IER-244 (MAY 1980).
C     MARK 9 REVISED. IER-317 (SEP 1981).
C     MARK 11B REVISED. IER-456 (SEP 1984).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     E04ABZ (NEWPTQ), AN ALGORITHM FOR FINDING A STEPLENGTH, CALLED
C     REPEATEDLY BY NPL OPTIMIZATION ROUTINES WHICH REQUIRE A STEP-
C     LENGTH TO BE COMPUTED USING QUADRATIC INTERPOLATION.
C     THE PARAMETERS SET UP BEFORE THE CALL OF E04ABZ CONTAIN
C     INFORMATION ABOUT THE INTERVAL IN WHICH A LOWER POINT IS TO BE
C     FOUND AND FROM THIS E04ABZ PRODUCES A POINT AT WHICH THE
C     FUNCTION CAN BE EVALUATED OUTSIDE THIS SUBROUTINE.
C     THE VALUE OF THE INTEGER PARAMETER ILOC DETERMINES THE PATH
C     TAKEN THROUGH THE CODE. FOR A FURTHER DESCRIPTION OF ILOC
C     AND THE OTHER PARAMETERS SEE NPL ALGORITHMS LIBRARY REF. NO.
C     E4/15/F.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN,
C     MARGARET H. WRIGHT AND ENID M. LONG
C     D.N.A.C. NATIONAL PHYSICAL LABORATORY, ENGLAND
C
C     **************************************************************
C
C
C     BRANCH TO APPROPRIATE SECTION OF CODE DEPENDING ON THE
C     VALUE OF THE FLAG ILOC
C     THE SIGNIFICANCE OF THE FLAGS ILOC AND ITEST ARE DESCRIBED IN
C     NPL ALGORITHMS LIBRARY DOCUMENT REF. NO. E4/15/F.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, B, B1, D, E, EPS, ETA, FA, FMIN, FU, FV, FW,
     *                  GTEST1, GTEST2, GU, OLDF, RR, SCXBD, SFTBND, SS,
     *                  T, TOL, U, XLAMDA, XMIN, XV, XW
      INTEGER           ILOC, ITEST
C     .. Local Scalars ..
      DOUBLE PRECISION  A1, D1, D2, Q, R, S, XM
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      GO TO (20,40,40,480,460) ILOC
C
C     ILOC = 1
C
C     CHECK INPUT PARAMETERS
C
   20 ITEST = 2
      TOL = T
      IF (U.LE.0.0D+0 .OR. 0.5D0*XLAMDA.LE.T .OR. GU.GT.0.0D+0) RETURN
      ITEST = 1
C
C     A AND B DEFINE THE INTERVAL OF UNCERTAINTY. XMIN DENOTES
C     THE LOWEST POINT OBTAINED SO FAR, XW THE LAST VALUE
C     OF XMIN AND XV THE SCALED VALUE OF ALPHA CORRESPONDING TO
C     THE HIGHEST FUNCTION VALUE OF THE THREE POINTS THROUGH
C     WHICH A PARABOLA MAY BE FITTED. INITIALIZE A, XV, XW, XMIN
C     AT ORIGIN AND CORRESPONDING FUNCTION VALUES AT LATEST
C     ESTIMATE OF MINIMUM.
C
      XMIN = 0.0D+0
      XW = 0.0D+0
      XV = 0.0D+0
      A = 0.0D+0
      OLDF = FU
      FMIN = FU
      FW = FU
      FV = FU
      FA = FU
      D = U
C
C     THE PARAMETER RR HAS TWO USES DURING THE EXECUTION OF THIS
C     SUBROUTINE.  INITIALLY THE SIGN OF RR INDICATES WHETHER OR NOT
C     THE MINIMUM HAS BEEN BRACKETED. LATER, WHEN A POINT SATISFYING
C     THE GTEST2 CRITERION HAS BEEN FOUND, RR IS USED TO COMPUTE A
C     STEPLENGTH WHICH SATISFIES THE SECOND CRITERION INVOLVING
C     GTEST1.
C
      RR = -1.0D+0
C
C     SET UP XBND AS A BOUND ON THE STEP TO BE TAKEN. (XBND IS NOT
C     COMPUTED EXPLICITLY BUT SCXBD IS ITS SCALED VALUE.) SET THE
C     UPPER BOUND ON THE INTERVAL OF UNCERTAINTY INITIALLY TO
C     XLAMDA + TOL(XLAMDA).
C
      SCXBD = XLAMDA
      B = SCXBD + EPS*ABS(SCXBD) + T
      E = 2.0D+0*B
      B1 = B
C
C     COMPUTE THE CONSTANTS REQUIRED FOR THE TWO CONVERGENCE
C     CRITERIA.
C
      GTEST1 = -1.0D-4*GU
      GTEST2 = -ETA*GU
C
C     SET ILOC TO INDICATE THAT ONLY TWO POINTS ARE AVAILABLE
C
      ILOC = 2
      GO TO 400
C
C     ILOC = 2 OR 3
C
C     UPDATE A, B, XV, XW, AND XMIN.
C
   40 IF (FU.GT.FMIN) GO TO 100
C
C     IF FUNCTION VALUE NOT INCREASED, NEW POINT BECOMES
C     NEXT ORIGIN AND OTHER POINTS ARE SCALED ACCORDINGLY.
C
      IF (U.LT.0.0D+0) GO TO 60
      A = 0.0D+0
      FA = FMIN
      IF (XW.EQ.XV .AND. FMIN.EQ.FU) RR = 1.0D+0
      GO TO 80
   60 B = 0.0D+0
      RR = 1.0D+0
   80 XV = XW
      FV = FW
      FW = FMIN
      FMIN = FU
      XMIN = XMIN + U
      A = A - U
      B = B - U
      XV = XV - U
      XW = 0.0D+0 - U
C
C     THIS MAY BE CHANGED TO XW = - U IF THE COMPUTER IS
C     SUCH THAT  - U AND 0.0 - U ARE IDENTICAL.
C
      SCXBD = SCXBD - U
      TOL = EPS*ABS(XMIN) + T
      GO TO 180
C
C     IF FUNCTION VALUE INCREASED, ORIGIN REMAINS UNCHANGED
C     BUT OTHER POINTS MAY BE INTERCHANGED.
C
  100 IF (U.GE.0.0D+0) GO TO 120
      A = U
      FA = FU
      GO TO 140
  120 B = U
      RR = 1.0D+0
  140 IF (FU.GT.FW .AND. XW.NE.0.0D+0) GO TO 160
      XV = XW
      FV = FW
      XW = U
      FW = FU
      GO TO 180
  160 XV = U
      FV = FU
  180 XM = 5.0D-1*(A+B)
C
C     CHECK TERMINATION CRITERIA.
C
      IF (0.5D0*ABS(XM).LE.TOL-0.25D0*(B-A)
     *    .OR. XMIN+B.LE.SFTBND .OR. FA-FMIN.LE.ABS(A)
     *    *GTEST2 .AND. FMIN.LT.OLDF .AND. (ABS(XMIN-XLAMDA)
     *    .GT.TOL .OR. RR.LT.0.0D+0)) GO TO 440
      R = 0.0D+0
      Q = 0.0D+0
      S = 0.0D+0
      IF (ABS(E).LE.TOL) GO TO 240
C
C     FIT PARABOLA THROUGH XMIN, XV, XW.
C
      IF (ILOC.NE.2) GO TO 200
C
C     SPECIAL CASE. ONLY TWO POINTS ARE AVAILABLE FOR
C     QUADRATIC INTERPOLATION
C
      Q = 2.0D+0*(FW-FMIN-XW*GU)
      S = GU*XW*XW
      IF (XMIN.NE.0.0D+0) S = (2.0D+0*(FMIN-FW)+XW*GU)*XW
      GO TO 220
  200 R = XW*(FV-FMIN)
      Q = XV*(FW-FMIN)
      S = R*XW - Q*XV
      Q = 2.0D+0*(Q-R)
  220 IF (Q.GT.0.0D+0) S = -S
      IF (Q.LE.0.0D+0) Q = -Q
      R = E
C
C     IF THE LAST STEP EXPANDED THE INTERVAL OR THE MINIMUM HAS
C     ALREADY BEEN BRACKETED SET E AS THE LAST STEP TAKEN
C
      IF (D.NE.B1 .OR. RR.GT.0.0D+0) E = D
C
C     CONSTRUCT AN ARTIFICIAL BOUND ON THE ESTIMATED STEPLENGTH.
C
  240 A1 = A
      B1 = B
      IF (XMIN.NE.A) GO TO 260
      D = XM
      GO TO 340
  260 IF (RR.GT.0.0D+0) GO TO 280
      D = -4.0D+0*A
      IF (D.GE.SCXBD) D = SCXBD
      GO TO 320
C
C     DETERMINE INTERVAL OF LENGTH D2 IN WHICH TO SET
C     ARTIFICIAL BOUND.
C
  280 D1 = A
      D2 = B
      IF (ABS(D2).GT.TOL .AND. (XW.LE.0.0D+0 .OR. ABS(D1).LE.TOL))
     *    GO TO 300
      U = D1
      D1 = D2
      D2 = U
  300 U = -D1/D2
      IF (U.GE.1.0D+0) D = 5.0D+0*D2*(1.0D-1+1.0D+0/U)/1.1D+1
      IF (U.LT.1.0D+0) D = 5.0D-1*D2*SQRT(U)
C
C     IF THE MINIMUM IS BRACKETED BY XV AND XW THE STEP MUST LIE
C     WITHIN (A, B).
C
  320 IF (XW.LT.0.0D+0 .AND. XV.GT.0.0D+0 .OR. XW.GT.0.0D+0 .AND. XV.LT.
     *    0.0D+0) GO TO 340
C
C     IF THE MINIMUM IS NOT BRACKETED BY XV AND XW THE STEP MUST LIE
C     WITHIN (A1, B1).
C
      IF (D.LE.0.0D+0) A1 = D
      IF (D.GT.0.0D+0) B1 = D
C
C     REJECT THE STEP OBTAINED BY INTERPOLATION IF IT LIES OUTSIDE
C     THE REQUIRED INTERVAL OR IT IS GREATER THAN HALF THAT
C     OBTAINED DURING THE LAST-BUT-ONE ITERATION.
C
  340 IF (ABS(S).GE.ABS(5.0D-1*Q*R) .OR. S.LE.Q*A1 .OR. S.GE.Q*B1)
     *    GO TO 360
C
C     A PARABOLIC INTERPOLATION STEP.
C
      D = S/Q
C
C     F MUST NOT BE EVALUATED TOO CLOSE TO A OR B.
C
      IF (0.5D0*(D-A).GE.TOL .AND. 0.5D0*(B-D).GE.TOL) GO TO 380
      D = TOL
      IF (XM.LE.0.0D+0) D = -TOL
      GO TO 380
C
C     A NON-INTERPOLATION STEP.
C
  360 E = B
      IF (XM.LE.0.0D+0) E = A
  380 ILOC = 3
C
C     CHECK THAT THE NEW STEP LENGTH WILL NOT BE GREATER THAN
C     XLAMDA.
C
  400 IF (D.LT.SCXBD) GO TO 420
C
C     REPLACE THE STEP LENGTH BY THE SCALED BOUND (SO AS TO COMPUTE
C     THE NEW POINT ON THE BOUNDARY.
C
      D = SCXBD
C
C     MOVE SCXBD TO THE LEFT SO THAT NEWXBND + TOL(NEWXBND) = XBND.
C
      SCXBD = (SCXBD-TOL)/(1.0D+0+EPS)
  420 U = D
      IF (ABS(D).LT.TOL .AND. D.LE.0.0D+0) U = -TOL
      IF (ABS(D).LT.TOL .AND. D.GT.0.0D+0) U = TOL
      ITEST = 1
      RETURN
C
C     THE FIRST CONVERGENCE CRITERION HAS BEEN SATISFIED. NOW CHECK
C     THAT THE FUNCTION VALUE HAS BEEN REDUCED SUFFICIENTLY. THE
C     VARIABLE RR IS NOW USED TO REDUCE THE STEP LENGTH.
C
  440 D = RR
      RR = XMIN
      SS = 5.0D-1
      FU = FMIN
      IF (XMIN.EQ.0.0D+0) XMIN = T
C
C     IF XMIN LIES WITHIN TOL OF THE BOUNDARY AND THE MINIMUM HAS
C     BEEN BRACKETED, THEN RECOMPUTE THE POINT ON THE BOUNDARY.
C
  460 IF (ABS(XMIN-XLAMDA).GE.TOL .OR. XMIN.EQ.T) GO TO 480
      IF (SCXBD.LT.0.0D+0 .AND. XW.LT.0.0D+0 .AND. XV.LT.0.0D+0)
     *    XMIN = XLAMDA
      IF (D.LT.0.0D+0) GO TO 480
      U = 0.0D+0
      ILOC = 4
      ITEST = 1
      RETURN
C
C     CHECK THAT THE NEW POINT SATISFIES SAFEGUARD CONDITIONS.
C     IF NECESSARY ATTEMPT TO FIND A SUFFICIENTLY LOWER POINT
C     BY SUCCESSIVELY DECREASING THE STEPLENGTH.
C
  480 IF (XMIN+B.GT.SFTBND) GO TO 500
      ITEST = 4
C
C     ITEST = 4 IMPLIES THAT THE REQUIRED STEP LENGTH IS SMALLER
C     THAN SFTBND.
C
      RETURN
  500 IF (OLDF-FU.LE.GTEST1*XMIN) GO TO 520
      FMIN = FU
      ITEST = 0
C
C     THE ALGORITHM HAS SUCCESSFULLY FOUND A SUFFICIENTLY LOWER
C     POINT.
C
      RETURN
  520 IF (XMIN.NE.T) GO TO 540
      ITEST = 3
C
C     DESPITE REPEATED REDUCTIONS IN THE STEP SIZE, A LOWER POINT
C     COULD NOT BE FOUND.
C
      RETURN
C
C     A SUFFICIENT REDUCTION IN THE FUNCTION VALUE HAS NOT YET BEEN
C     FOUND, TRY A FURTHER REDUCTION IN THE STEP LENGTH.
C
  540 XMIN = RR*SS
      SS = SS*SS
      IF (XMIN.LT.T) XMIN = T
      ITEST = 1
      U = 0.0D+0
      ILOC = 5
      RETURN
C
C     END OF E04ABZ (NEWPTQ)
C
      END
      SUBROUTINE G01AAF(N,X,IWT,WT,XMEAN,S2,S3,S4,XMIN,XMAX,WSUM,IFAIL)
C     MARK 4 RELEASE NAG COPYRIGHT 1974.
C     MARK 4.5 REVISED
C     MARK 10D REVISED. IER-430 (FEB 1984)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C
C     G01AAF - THE ROUTINE RETURNS THE MEAN, THE STANDARD
C     DEVIATION, THE COEFFICIENTS OF SKEWNESS AND
C     KURTOSIS, AND THE LARGEST AND SMALLEST VALUES
C     FROM THE VALUES IN THE ARRAY X - OPTIONALLY
C     USING WEIGHTS IN THE ARRAY WT
C
C     USES NAG LIBRARY ROUTINE P01AAF
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G01AAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  S2, S3, S4, WSUM, XMAX, XMEAN, XMIN
      INTEGER           IFAIL, IWT, N
C     .. Array Arguments ..
      DOUBLE PRECISION  WT(N), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  A2, ATEMP, BTEMP, SMALL, WTEMP, XTEMP
      INTEGER           I, IERR, K, LOW
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AKF
      INTEGER           P01ABF
      EXTERNAL          X02AKF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      IF (N) 20, 20, 80
   20 IERR = 1
   40 IWT = 0
   60 IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,P01REC)
      RETURN
   80 IF (IWT-1) 100, 140, 100
C     INITIALISE  WEIGHTS IF NOT GIVEN
  100 DO 120 I = 1, N
         WT(I) = 1.0D0
  120 CONTINUE
      LOW = 1
      WTEMP = 1.0D0
      IWT = N
      GO TO 300
C     FIND FIRST POSITIVE WEIGHT--CHECK OTHERS NON-NEG
  140 IWT = 0
      DO 260 I = 1, N
         WTEMP = WT(I)
         IF (WTEMP) 280, 260, 160
  160    IF (I-N) 180, 240, 240
  180    LOW = I + 1
         DO 220 K = LOW, N
            IF (WT(K)) 280, 220, 200
  200       IWT = IWT + 1
  220    CONTINUE
  240    LOW = I
         IWT = IWT + 1
         GO TO 300
  260 CONTINUE
  280 IERR = 3
      GO TO 40
C     INITIALISE MEANS ETC
  300 S2 = 0.0D0
      S3 = 0.0D0
      S4 = 0.0D0
      XMEAN = X(LOW)
      WSUM = WTEMP
      XMIN = XMEAN
      XMAX = XMEAN
      IF (IWT-1) 280, 440, 320
  320 A2 = 0.0D0
      SMALL = SQRT(SQRT(X02AKF()))
C     LOOP FOR CALCULATIONS
      LOW = LOW + 1
      DO 360 I = LOW, N
         WTEMP = WT(I)
         IF (WTEMP) 280, 360, 340
  340    XTEMP = X(I)
         IF (XTEMP.LT.XMIN) XMIN = XTEMP
         IF (XTEMP.GT.XMAX) XMAX = XTEMP
         A2 = 2.0D0*WSUM*WTEMP + A2
         WSUM = WTEMP + WSUM
         XTEMP = XTEMP - XMEAN
         ATEMP = WTEMP*XTEMP
         IF (ABS(ATEMP).LE.SMALL) ATEMP = 0.0D0
         BTEMP = ATEMP/WSUM
         ATEMP = ATEMP*(XTEMP-BTEMP)
         WTEMP = 3.0D0*BTEMP*S2
         S4 = BTEMP*(2.0D0*WTEMP-4.0D0*S3) +
     *        ATEMP*(-3.0D0*ATEMP/WSUM+XTEMP*XTEMP) + S4
         S3 = -WTEMP + ATEMP*(-2.0D0*BTEMP+XTEMP) + S3
         S2 = ATEMP + S2
         XMEAN = BTEMP + XMEAN
  360 CONTINUE
      IF (S2) 400, 400, 380
  380 ATEMP = WSUM/A2
      A2 = S2*ATEMP
      S2 = SQRT(A2)
      S3 = ATEMP*S3/(A2*S2)
      S4 = ATEMP*S4/(A2*A2) - 3.0D0
      GO TO 420
  400 S3 = 0.0D0
      S4 = 0.0D0
      S2 = 0.0D0
  420 IFAIL = 0
      RETURN
  440 IERR = 2
      GO TO 60
      END
      DOUBLE PRECISION FUNCTION G01FFF(P,A,B,TOL,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15 REVISED. IER-933 (APR 1991).
C     MARK 18 REVISED. IER-1873 (MAY 1997).
C
C     Based on:
C     Algorithm as 91 Appl. Statist. (1975) Vol.24, P.385.
C     Computes the deviate associated with the lower tail
C     probability P, of the Gamma distribution with shape
C     parameter A and scale parameter B.
C
C     .. Parameters ..
      DOUBLE PRECISION                 ZERO, HALF, ONE, TWO, THREE, SIX,
     *                                 FIFTY, BIG
      PARAMETER                        (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,
     *                                 TWO=2.0D0,THREE=3.0D0,SIX=6.0D0,
     *                                 FIFTY=50.0D0,BIG=1.0D6)
      DOUBLE PRECISION                 C1, C2, C3, C4, C5, C6, C7, C8,
     *                                 C9, C10, C11, C12, C13, C14, C15,
     *                                 C16, C17, C18, C19, C20, C21,
     *                                 C22, C23, C24, C25, C26, C27,
     *                                 C28, C29, C30, C31, C32, C33,
     *                                 C34, C35, C36, C37, C38
      PARAMETER                        (C1=0.01D0,C2=0.111111D0,
     *                                 C3=0.16D0,C4=0.4D0,C5=0.62D0,
     *                                 C6=4.4D0,C7=4.67D0,C8=6.66D0,
     *                                 C9=6.73D0,C10=13.32D0,C11=60.0D0,
     *                                 C12=70.0D0,C13=84.0D0,
     *                                 C14=105.0D0,C15=120.0D0,
     *                                 C16=127.0D0,C17=140.0D0,
     *                                 C18=175.0D0,C19=210.0D0,
     *                                 C20=252.0D0,C21=264.0D0,
     *                                 C22=294.0D0,C23=346.0D0,
     *                                 C24=420.0D0,C25=462.0D0,
     *                                 C26=606.0D0,C27=672.0D0,
     *                                 C28=707.0D0,C29=735.0D0,
     *                                 C30=889.0D0,C31=932.0D0,
     *                                 C32=966.0D0,C33=1141.0D0,
     *                                 C34=1182.0D0,C35=1278.0D0,
     *                                 C36=1740.0D0,C37=2520.0D0,
     *                                 C38=5040.0D0)
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G01FFF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 A, B, P, TOL
      INTEGER                          IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 BB, C, D, DD, EPS, G, GAM,
     *                                 GAMLOG, P1, P2, Q, Q2, S1, S2,
     *                                 S3, S4, S5, S6, T, UFLOW, X
      INTEGER                          IERROR, IFAIL2, IFAULT, ITER,
     *                                 MAXIT
C     .. Local Arrays ..
      CHARACTER*80                     REC(1)
C     .. External Functions ..
      DOUBLE PRECISION                 G01FAZ, S14ABF, X02AJF, X02AMF
      INTEGER                          P01ABF
      EXTERNAL                         G01FAZ, S14ABF, X02AJF, X02AMF,
     *                                 P01ABF
C     .. External Subroutines ..
      EXTERNAL                         S14BAF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, EXP, LOG, MAX, SQRT
C     .. Executable Statements ..
C
C     Test arguments and initialize.
C
      G01FFF = ZERO
      IF (P.LT.ZERO .OR. P.GE.ONE) THEN
         WRITE (REC,FMT=99999) P
         IERROR = 1
      ELSE IF (A.LE.ZERO) THEN
         WRITE (REC,FMT=99998) A
         IERROR = 2
      ELSE IF (A.GT.BIG) THEN
         WRITE (REC,FMT=99993) A
         IERROR = 2
      ELSE IF (B.LE.ZERO) THEN
         WRITE (REC,FMT=99994) B
         IERROR = 2
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.EQ.0 .AND. P.NE.0.0D0) THEN
         EPS = MAX(X02AJF(),1.0D-18)*FIFTY
         IF (TOL.GT.EPS .AND. TOL.LT.ONE) EPS = TOL
         MAXIT = 100
         UFLOW = LOG(X02AMF())
         DD = LOG(TWO)
         IFAULT = 1
         G = S14ABF(A,IFAULT)
         C = A - ONE
         IF (A.LT.-C5*LOG(P)) THEN
            GAMLOG = ONE/A*(LOG(P)+LOG(A)+G+A*DD)
            IF (GAMLOG.LT.UFLOW) THEN
               GAM = ZERO
               IERROR = 3
               WRITE (REC,FMT=99996)
               GO TO 100
            ELSE
               GAM = EXP(GAMLOG)
            END IF
            IF (GAM.LT.EPS) GO TO 80
            GO TO 40
         END IF
C
C        Find starting values
C
         IF (A.LE.C3) THEN
            GAM = C4
            D = LOG(ONE-P)
   20       Q = GAM
            P1 = ONE + GAM*(C7+GAM)
            P2 = GAM*(C9+GAM*(C8+GAM))
            T = -HALF + (C7+TWO*GAM)/P1 - (C9+GAM*(C10+THREE*GAM))/P2
            GAM = GAM - (ONE-EXP(D+G+HALF*GAM+C*DD)*P2/P1)/T
            IF (ABS(Q/GAM-ONE).GT.C1) GO TO 20
            GO TO 40
         END IF
C
C        Call G01FAZ (P HAS BEEN CHECKED)
C
         IFAULT = 1
         X = G01FAZ(P,IFAULT)
C
C        Starting approximation using Wilson and Hilferty estimate.
C
         P1 = C2/A
         GAM = TWO*A*(X*SQRT(P1)+ONE-P1)**3
C
C        Starting approximation for P tending to 1.0.
C
         IF (GAM.GT.C6*A+SIX) GAM = -TWO*(LOG(ONE-P)-C*LOG(HALF*GAM)+G)
C
C        Call S14BAF and calculate seven term taylor series.
C
   40    ITER = 0
   60    CONTINUE
         ITER = ITER + 1
         Q = GAM
         P1 = HALF*GAM
         IFAIL2 = 1
         CALL S14BAF(A,P1,EPS,P2,Q2,IFAIL2)
         IF (IFAIL2.NE.0) THEN
            IERROR = 5
            WRITE (REC,FMT=99995)
            GO TO 100
         END IF
         P2 = P - P2
         T = P2*EXP(A*DD+G+P1-C*LOG(GAM))
         BB = T/GAM
         D = HALF*T - BB*C
         S1 = (C19+D*(C17+D*(C14+D*(C13+D*(C12+C11*D)))))/C24
         S2 = (C24+D*(C29+D*(C32+D*(C33+C35*D))))/C37
         S3 = (C19+D*(C25+D*(C28+C31*D)))/C37
         S4 = (C20+D*(C27+C34*D)+C*(C22+D*(C30+C36*D)))/C38
         S5 = (C13+C21*D+C*(C18+C26*D))/C37
         S6 = (C15+C*(C23+C16*C))/C38
         GAM = GAM + T*(ONE+HALF*T*S1-BB*C*
     *         (S1-BB*(S2-BB*(S3-BB*(S4-BB*(S5-BB*S6))))))
         IF (ABS(Q/GAM-ONE).GT.EPS) THEN
            IF (ITER.GE.MAXIT) THEN
               IERROR = 4
               WRITE (REC,FMT=99997)
            ELSE
               GO TO 60
            END IF
         END IF
   80    G01FFF = HALF*GAM*B
      END IF
  100 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,1,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, P.lt.0.0 .OR. P.ge.1.0: P = ',D13.5)
99998 FORMAT (1X,'** On entry, A.le.0.0 : A = ',D13.5)
99997 FORMAT (1X,'** Solution fails to converge.')
99996 FORMAT (1X,'** P is too close to 0.0 or 1.0')
99995 FORMAT (1X,'** Convergence failure in calculating gamma integral.'
     *       )
99994 FORMAT (1X,'** On entry, B .le. 0.0 : B = ',D13.5)
99993 FORMAT (1X,'** On entry, A is too large: A = ',D13.5)
      END
      DOUBLE PRECISION FUNCTION G01FAZ(P,IFAIL)
C     MARK 18 RELEASE. NAG COPYRIGHT 1997.
C
C     G01FAZ RETURNS THE DEVIATE ASSOCIATED WITH THE LOWER
C     TAIL PROBABILITY P FROM THE STANDARD NORMAL
C     DISTRIBUTION.
C
C     WRITTEN BY N.M.MACLAREN AND COLLEAGUES
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C
C     RELATIVE ACCURACY IS 5.0E-13
C
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G01FAZ')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 P
      INTEGER                          IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 P1, P2, P3, P4, P5, P6, P7, P8,
     *                                 Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8,
     *                                 R2PI, S, T2B1, T2B2, T2B3, T2B4,
     *                                 T2B5, T2B6, T2T1, T2T2, T2T3,
     *                                 T2T4, T2T5, T2T6, T3B1, T3B2,
     *                                 T3B3, T3B4, T3B5, T3B6, T3B7,
     *                                 T3T1, T3T2, T3T3, T3T4, T3T5,
     *                                 T3T6, T3T7, X, Y, YD
C     .. Local Arrays ..
      CHARACTER*1                      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION                 X01AAF
      INTEGER                          P01ABF
      EXTERNAL                         X01AAF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, MIN, SIGN, LOG, SQRT
C     .. Data statements ..
      DATA                             P1, P2, P3, P4, P5, P6, P7,
     *                                 P8/0.1000000000000000D+1,
     *                                 -0.4163362762616374D+2,
     *                                 0.6494128979404664D+2,
     *                                 -0.3689355386573687D+2,
     *                                 0.8984374887949291D+1,
     *                                 -0.8486205099916682D+0,
     *                                 0.1799227452322391D-1,
     *                                 0.1243484778425483D-3/
      DATA                             Q1, Q2, Q3, Q4, Q5, Q6, Q7,
     *                                 Q8/0.1000000000000000D+1,
     *                                 -0.4180029429283041D+2,
     *                                 0.7184967217618504D+2,
     *                                 -0.4645534714071778D+2,
     *                                 0.1357696314152835D+2,
     *                                 -0.1714887182301168D+1,
     *                                 0.6929108687616727D-1, 0.15D-8/
      DATA                             T2T1, T2T2, T2T3, T2T4, T2T5,
     *                                 T2T6/8.416212335733377D-1,
     *                                 -6.002063030711376D0,
     *                                 -4.722617278101355D-1,
     *                                 68.25030402764662D0,
     *                                 -96.23966784668191D0,
     *                                 -10.47655255702662D0/
      DATA                             T2B1, T2B2, T2B3, T2B4, T2B5,
     *                                 T2B6/1.0D0, -11.37563668160269D0,
     *                                 41.33878042998575D0,
     *                                 -43.59278882675467D0,
     *                                 -21.25589975172773D0,
     *                                 25.08610603077810D0/
      DATA                             T3T1, T3T2, T3T3, T3T4, T3T5,
     *                                 T3T6, T3T7/-3.012472607913056D0,
     *                                 2.982419254309647D0,
     *                                 22.70676420727861D0,
     *                                 7.670763790822552D0,
     *                                 5.519186629029667D-1,
     *                                 7.985443046076538D-3,
     *                                 4.150184151574655D-6/
      DATA                             T3B1, T3B2, T3B3, T3B4, T3B5,
     *                                 T3B6, T3B7/1.0D0,
     *                                 1.239758125817922D0,
     *                                 -1.205350980555889D1,
     *                                 -12.02359058219926D0,
     *                                 -2.373035796200643D0,
     *                                 -1.193440282031508D-1,
     *                                 -1.216250189896074D-3/
C     .. Executable Statements ..
      IF (P.LE.0.0D0 .OR. P.GE.1.0D0) GO TO 60
      IFAIL = 0
      X = P - 0.5D0
      IF (ABS(X).GT.0.3D0) GO TO 20
C
C     BREAK UP RANGE
C     FOR 0.2 LE P LE 0.8 WE USE A RATIONAL TCHEBYCHEV (P/Q)
C
      R2PI = SQRT(2.0D0*X01AAF(0.0D0))
      S = X*R2PI
      X = S*S
      Y = P1 + X*(P2+X*(P3+X*(P4+X*(P5+X*(P6+X*(P7+X*P8))))))
      YD = Q1 + X*(Q2+X*(Q3+X*(Q4+X*(Q5+X*(Q6+X*(Q7+X*Q8))))))
      Y = Y/YD
      G01FAZ = Y*S
      RETURN
C
C     FOR 0.08 LE P LE 0.2 OR 0.8 LE P LE 0.92
C     WE USE RATIONAL TCHEBYCHEV (T2T/T2B)
C
   20 S = SIGN(1.0D0,X)
      IF (ABS(X).GT.0.42D0) GO TO 40
C
C     BREAK UP RANGE SOME MORE
C
      X = ABS(X) - 0.3D0
      Y = T2T1 + X*(T2T2+X*(T2T3+X*(T2T4+X*(T2T5+X*T2T6))))
      YD = T2B1 + X*(T2B2+X*(T2B3+X*(T2B4+X*(T2B5+X*T2B6))))
      Y = Y/YD
      G01FAZ = Y*S
      RETURN
C
C     THE FOLLOWING CASE HANDLES ASYMPTOTIC BEHAVIOUR.
C     UNFORTUNATELY WE MUST MAKE A TRANSFORMATION.
C
   40 X = SQRT((-2.0D0)*LOG(MIN(P,1.0D0-P)))
      Y = T3T1 + X*(T3T2+X*(T3T3+X*(T3T4+X*(T3T5+X*(T3T6+X*T3T7)))))
      YD = T3B1 + X*(T3B2+X*(T3B3+X*(T3B4+X*(T3B5+X*(T3B6+X*T3B7)))))
      Y = Y/YD + X
      G01FAZ = S*Y
      RETURN
C
C     ERROR EXITS
C
   60 IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      G01FAZ = 0.0D0
      RETURN
      END
      SUBROUTINE M01DAF(RV,M1,M2,ORDER,IRANK,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     M01DAF RANKS A VECTOR OF REAL NUMBERS IN ASCENDING
C     OR DESCENDING ORDER.
C
C     M01DAF USES A VARIANT OF LIST-MERGING, AS DESCRIBED
C     BY KNUTH. THE ROUTINE TAKES ADVANTAGE OF NATURAL
C     ORDERING IN THE DATA, AND USES A SIMPLE LIST INSERTION
C     IN A PREPARATORY PASS TO GENERATE ORDERED LISTS OF
C     LENGTH AT LEAST 10. THE RANKING IS STABLE: EQUAL ELEMENTS
C     PRESERVE THEIR ORDERING IN THE INPUT DATA.
C
C     THE MINIMUM LENGTH OF THE LISTS AT THE END OF THE
C     PREPARATORY PASS IS DEFINED BY THE VARIABLE MAXINS.
C
C     WRITTEN BY N.M.MACLAREN, UNIVERSITY OF CAMBRIDGE.
C     REVISED BY NAG CENTRAL OFFICE.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='M01DAF')
      INTEGER           MAXINS
      PARAMETER         (MAXINS=10)
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M1, M2
      CHARACTER*1       ORDER
C     .. Array Arguments ..
      DOUBLE PRECISION  RV(M2)
      INTEGER           IRANK(M2)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, B, C
      INTEGER           I, I1, I2, IERR, ILIST, J, J1, J2, K, K1, K2, L,
     *                  LIST1, LIST2, NLAST, NPREV, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Executable Statements ..
C
C       CHECK THE ARGUMENTS AND DEAL WITH THE TRIVIAL CASE.
C
      IF (M2.LT.1 .OR. M1.LT.1 .OR. M1.GT.M2) THEN
         IERR = 1
         WRITE (P01REC,FMT=99999) M1, M2
         NREC = 2
      ELSE IF (ORDER.NE.'A' .AND. ORDER.NE.'a' .AND. ORDER.NE.'D' .AND.
     *         ORDER.NE.'d') THEN
         IERR = 2
         WRITE (P01REC,FMT=99998) ORDER
         NREC = 1
      ELSE IF (M1.EQ.M2) THEN
         IRANK(M2) = M2
         IERR = 0
      ELSE
         IERR = 0
C
C        INITIALISE, USING NATURAL RUNS IN BOTH DIRECTIONS AND
C        STRAIGHT LIST INSERTION FOR SMALL LISTS.
C
C        I  POINTS TO THE SMALLEST ELEMENT IN THE CURRENT LIST
C        J  POINTS TO THE LARGEST  ELEMENT IN THE CURRENT LIST
C        B  IS THE VALUE OF THE SMALLEST ELEMENT IN CURRENT LIST
C        C  IS THE VALUE OF THE LARGEST  ELEMENT IN CURRENT LIST
C
         ILIST = -1
         K = M1
         I = K
         J = K
         L = K + MAXINS
         B = RV(K)
         C = B
         DO 40 K = M1 + 1, M2
C
C           DEAL WITH ADDITIONS AT EITHER END.
C
            A = RV(K)
            IF (A.GE.C) THEN
               IRANK(J) = K
               J = K
               C = A
            ELSE IF (A.LT.B) THEN
               IRANK(K) = I
               I = K
               B = A
            ELSE
C
C              DO AN ASCENDING LIST INSERTION.
C
               IF (K.LT.L) THEN
                  I2 = I
   20             I1 = I2
                  I2 = IRANK(I1)
                  IF (A.GE.RV(I2)) GO TO 20
                  IRANK(I1) = K
                  IRANK(K) = I2
               ELSE
C
C                 ADD THE CURRENT LIST ON TO THE OTHERS.
C
                  IF (ILIST.LT.0) THEN
                     LIST1 = -I
                     ILIST = 0
                  ELSE IF (ILIST.EQ.0) THEN
                     LIST2 = -I
                     ILIST = 1
                     NPREV = NLAST
                  ELSE
                     IRANK(NPREV) = -I
                     NPREV = NLAST
                  END IF
C
                  NLAST = J
                  I = K
                  J = K
                  L = K + MAXINS
                  B = RV(K)
                  C = B
               END IF
            END IF
   40    CONTINUE
C
C        TIDY UP AT THE END.
C
         IRANK(J) = 0
         IF (ILIST.LT.0) THEN
            LIST1 = -I
            GO TO 280
         ELSE IF (ILIST.EQ.0) THEN
            LIST2 = -I
         ELSE
            IRANK(NPREV) = -I
         END IF
         IRANK(NLAST) = 0
C
C        AT THIS POINT:
C        LIST1 = - (INDEX OF LEAST ELEMENT IN THE FIRST LIST)
C        LIST2 = - (INDEX OF LEAST ELEMENT IN THE SECOND LIST)
C        FOR EACH K, IRANK(K) IS THE INDEX OF THE NEXT ELEMENT IN THE
C        CURRENT LIST, EXCEPT THAT, IF THERE IS NO SUCH ELEMENT,
C        IRANK(K) IS - (INDEX OF THE LEAST ELEMENT IN THE NEXT LIST
C        BUT 1)  OR 0 IF THERE IS NO SUCH LIST.
C
C        START MERGING LISTS BY PAIRS.
C
   60    ILIST = -1
         I = -LIST1
         J = -LIST2
   80    K = I
         IF (RV(I).GT.RV(J)) K = J
         IF (ILIST.LT.0) THEN
            LIST1 = -K
            ILIST = 0
         ELSE IF (ILIST.EQ.0) THEN
            LIST2 = -K
            ILIST = 1
            NLAST = L
         ELSE
            IRANK(NLAST) = -K
            NLAST = L
         END IF
C
C        MOVE ALONG THE LISTS UNTIL ONE FINISHES.
C
C        NEW VARIABLES I2, J2 AND K2 ARE USED INSTEAD OF I, J AND K
C        WITHIN THE INNERMOST BLOCK TO ENCOURAGE OPTIMISING COMPILERS TO
C        STORE THEM IN REGISTERS.
C         I2 POINTS TO THE CURRENT ELEMENT IN THE FIRST LIST
C         J2 POINTS TO THE CURRENT ELEMENT IN THE SECOND LIST
C         K2 POINTS TO THE CURRENT ELEMENT IN THE MERGED LIST
C
         I2 = I
         J2 = J
         IF (K.NE.I2) GO TO 140
  100    A = RV(J2)
         K2 = I2
  120    I2 = K2
         K2 = IRANK(I2)
         IF (K2.LE.0) GO TO 180
         IF (A.GE.RV(K2)) GO TO 120
         IRANK(I2) = J2
         I2 = K2
  140    A = RV(I2)
         K2 = J2
  160    J2 = K2
         K2 = IRANK(J2)
         IF (K2.LE.0) GO TO 200
         IF (A.GT.RV(K2)) GO TO 160
         IRANK(J2) = I2
         J2 = K2
         GO TO 100
C
C        ADD THE REMAINS OF ONE LIST TO THE OTHER.
C
  180    K = 1
         I1 = K2
         GO TO 220
  200    K = 2
         J1 = K2
  220    I = I2
         J = J2
         IF (K.EQ.1) THEN
C
C           FIRST LIST IS EXHAUSTED
C
            IRANK(I) = J
            I = -I1
            J1 = J
  240       J = J1
            J1 = IRANK(J)
            IF (J1.GT.0) GO TO 240
            L = J
            J = -J1
         ELSE
C
C           SECOND LIST IS EXHAUSTED
C
            IRANK(J) = I
            J = -J1
            I1 = I
  260       I = I1
            I1 = IRANK(I)
            IF (I1.GT.0) GO TO 260
            L = I
            I = -I1
         END IF
C
C        TIDY UP AND CARRY ON IF NOT FINISHED.
C
         IF ((I.NE.0) .AND. (J.NE.0)) GO TO 80
         IRANK(L) = 0
         K = I + J
         IF (ILIST.GT.0) THEN
            IRANK(NLAST) = -K
            GO TO 60
         ELSE IF (K.NE.0) THEN
            LIST2 = -K
            GO TO 60
         END IF
C
C        IF DESCENDING, REVERSE ALL POINTERS BETWEEN EQUALITY
C        BLOCKS.
C
  280    IF ((ORDER.EQ.'D') .OR. (ORDER.EQ.'d')) THEN
            I = 0
            J = -LIST1
  300       K = J
            K1 = K
            A = RV(K)
  320       K = K1
            K1 = IRANK(K)
            IF (K1.NE.0) THEN
               IF (A.EQ.RV(K1)) GO TO 320
            END IF
            IRANK(K) = I
            I = J
            J = K1
            IF (J.NE.0) GO TO 300
            LIST1 = -I
         END IF
C
C        CONVERT THE LIST FORM TO RANKS AND RETURN.
C
         K = M1
         I = -LIST1
  340    I1 = IRANK(I)
         IRANK(I) = K
         K = K + 1
         I = I1
         IF (I.GT.0) GO TO 340
C
      END IF
C
      IF (IERR.NE.0) THEN
         IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,P01REC)
      ELSE
         IFAIL = 0
      END IF
      RETURN
C
99999 FORMAT (' ** On entry, one or more of the following parameter va',
     *  'lues is illegal',/'    M1 =',I16,'  M2 =',I16)
99998 FORMAT (' ** On entry, ORDER has an illegal value: ORDER = ',A1)
      END
