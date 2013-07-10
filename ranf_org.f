*###[ ranf: random number generator:
	DOUBLE PRECISION function ranf(dummy)
***#[*comment:***********************************************************
*									*
*	Switchyard routine which is compatible with the interface used	*
*	by VEGAS and calls one of the random number generators		*
*	described below.  The numbers are the same as in Kankaala's	*
*	thesis (the missing ones are commercial).			*
*									*
*	Input:	dummy	integer		  to fool optimizers which want	*
*					  to take ranf out of loops.	*
*	Parameter:							*
*		ifunc	integer		  1: GGL			*
*					  5: R250			*
*					  6: RAN3			*
*					  7: RANMAR			*
*					  8: RCARRY			*
*					  9: RANLUX			*
*					  10: RICHTMEYER		*
*									*
*	Output:	ranf	double precision  random number in (0,1)	*
*									*
***#]*comment:***********************************************************
*  #[ declarations:
	implicit none
*
*	arguments
*
	integer dummy
*
*	local variables
*
	integer ifunc,init
	DOUBLE PRECISION GGL,ran3,R250,ds,x1(1)
	save ds,init
*
*	parameter,data
*
	parameter(ifunc=5)
	data init /0/
	data ds /667790.d0/
*
*  #] declarations:
*  #[ switchyard:
*
	if ( ifunc.eq.1 ) then
	    if ( init.eq.0 ) then
		write(*,'(a)')'# ranf: using GGL generator - '
		write(*,'(a)')'#       watch the Marsiglia planes...'
	    endif
	    init = init+1
	    if ( init.eq.2147483646 ) print *,'ranf: error: ',
     +		'cycle length exhausted!!'
	    ranf = GGL(ds)
	elseif ( ifunc.eq.5 ) then
	    if ( init.eq.0 ) then
***		write(*,'(a)')'# ranf: using R250 generator'
		init = 1
	    endif
	    ranf = R250(dummy)
	elseif ( ifunc.eq.6 ) then
	    if ( init.eq.0 ) then
		write(*,'(a)')'# ranf: using RAN3 - '
		write(*,'(a)')'#       not too good, this one'
		init = 1
	    endif
	    ranf = ran3(dummy)
	elseif ( ifunc.eq.7 ) then
	    if ( init.eq.0 ) then
		write(*,'(a)')'# ranf: using RANMAR - '
		write(*,'(a)')'#       not too good, this one'
		init = 1
	    endif
	    call RANMAR(x1,1)
	    ranf = x1(1)
	elseif ( ifunc.eq.8 ) then
	    if ( init.eq.0 ) then
		write(*,'(a)')'# ranf: using RCARRY - '
		write(*,'(a)')'#       not too good, this one'
		init = 1
	    endif
	    call RCARRY(x1,1)
	    ranf = x1(1)
	elseif ( ifunc.eq.9 ) then
	    if ( init.eq.0 ) then
		write(*,'(a)')'# ranf: using RANLUX - very very nice'
		init = 1
	    endif
	    call RANLUX(x1,1)
	    ranf = x1(1)
	elseif ( ifunc.eq.10 ) then
	    if ( init.eq.0 ) then
		write(*,'(a)')'# ranf: using richtmeyer quasi-random '//
     +			'I hope you know what you are doing ''cause I'//
     +			'''m not'
		init = 1
	    endif
	    call RANLUX(x1,1)
	    ranf = x1(1)
	else
	    print *,'ranf: error: unknown value for ifunc: ',ifunc
	    stop
	endif
*
*  #] switchyard:
*###] ranf:
	end
*###[ resetranf: reset random number generator:
	subroutine resetranf(dummy)
***#[*comment:***********************************************************
*									*
*	Switchyard routine which is compatible with the interface used	*
*	by VEGAS and calls one of the random number generators		*
*	described below.  The numbers are the same as in Kankaala's	*
*	thesis (the missing ones are commercial).			*
*									*
*	Input:	dummy	integer		  to fool optimizers which want	*
*					  to take ranf out of loops.	*
*	Parameter:							*
*		ifunc	integer		  1: GGL			*
*					  5: R250			*
*					  6: RAN3			*
*					  7: RANMAR			*
*					  8: RCARRY			*
*					  9: RANLUX			*
*					  10: RICHTMEYER		*
*									*
*	Output:	ranf	double precision  random number in (0,1)	*
*									*
***#]*comment:***********************************************************
*  #[ declarations:
	implicit none
*
*	arguments
*
	integer dummy
*
*	local variables
*
	double precision xdum,resetr250
	integer ifunc,init,resetggli
	external resetggli,resetr250
*
*	parameter,data
*
	parameter(ifunc=5)
	data init /0/
*
*  #] declarations:
*  #[ switchyard:
*
	if ( ifunc.eq.1 ) then
	    xdum = resetggli(dummy)
	elseif ( ifunc.eq.5 ) then
	    xdum = resetr250(dummy)
	else
	    print *,'resetranf: error: unknown value for ifunc: ',ifunc
	    stop
	endif
*
*  #] switchyard:
*###] resetranf:
	end
*###[ comment:
*\\
*FORTRAN SOURCE CODES FOR SOME RANDOM NUMBER GENERATORS
*submitted by I. Vattulainen, K. Kankaala, J. Saarinen and
*T. Ala-Nissila.
*\\
*We have tested extensively several random number generators.
*The test results are available for example from the high energy
*physics preprint server at Florida State University:
*hep-lat@ftp.scri.fsu.edu preprint number 9304008.
*Further tests hasve been performed on various versions
*of RCARRY. These results are available also from
*hep-lat@ftp.scri.fsu.edu as a preprint number 9306008.
*The source codes are given below for those random number generator
*codes that are public domain ones.
*These are the exact implementations tested by us.
*\\
*###] comment:
*###[ GGL:
C ================================================================
C
C     RANDOM NUMBER GENERATOR    GGL
C     original reference for algorithm
C     P. Lewis, A. Goodman, J. Miller,IBM Sys. J., 2, 136 (1969)

      DOUBLE PRECISION   FUNCTION GGL (DS)

      DOUBLE PRECISION   DS, D2

      DATA               D2/2147483647.D0/

      DS = DMOD(16807.D0*DS,D2)
      GGL = DS/D2

      RETURN
*###] GGL:
      END
*###[ RAN3:
C ========================================================================
C
C RAN3:
C -----
C This is the RAN 3 from Numerical Recipes as proposed
C by D. Knuth and implemented by Press et al.

      double precision function ran3(idum)
c         implicit REAL*8*4(m)
c         parameter (mbig=4000000.,mseed=1618033.,mz=0.,fac=2.5e-7)
      parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1.e-9)
      dimension ma(55)
      data iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=mseed-iabs(idum)
        mj=mod(mj,mbig)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.mz)mk=mk+mbig
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.mz)ma(i)=ma(i)+mbig
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz)mj=mj+mbig
      ma(inext)=mj
      ran3=mj*fac
      return
*###] RAN3:
      end
*###[ RANMAR:
C =====================================================================
      SUBROUTINE RANMAR(RVEC,LENV)
C
C This version is identical to that in CPC software library
C
C Universal random number generator proposed by Marsaglia and Zaman
C in report FSU-SCRI-87-50
C        modified by F. James, 1988 and 1989, to generate a vector
C        of pseudorandom numbers RVEC of length LENV, and to put in
C        the COMMON block everything needed to specify currrent state,
C        and to add input and output entry points RMARIN, RMARUT.
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C!!!  Calling sequences for RANMAR:                                  ++
C!!!      CALL RANMAR (RVEC, LEN)   returns a vector RVEC of LEN     ++
C!!!                   32-bit random floating point numbers between  ++
C!!!                   zero and one.                                 ++
C!!!      CALL RMARIN(I1,N1,N2)   initializes the generator from one ++
C!!!                   32-bit integer I1, and number counts N1,N2    ++
C!!!                  (for initializing, set N1=N2=0, but to restart ++
C!!!                    a previously generated sequence, use values  ++
C!!!                    output by RMARUT)                            ++
C!!!      CALL RMARUT(I1,N1,N2)   outputs the value of the original  ++
C!!!                  seed and the two number counts, to be used     ++
C!!!                  for restarting by initializing to I1 and       ++
C!!!                  skipping N1*100000000+N2 numbers.              ++
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      double precision RVEC(*)
      COMMON/RASET1/U(97),C,I97,J97
      PARAMETER (MODCNS=1000000000)
      SAVE CD, CM, TWOM24, NTOT, NTOT2, IJKL
      DATA NTOT,NTOT2,IJKL/-1,0,0/
C
      IF (NTOT .GE. 0)  GO TO 50
C
C        Default initialization. User has called RANMAR without RMARIN.
      IJKL = 54217137
      NTOT = 0
      NTOT2 = 0
      KALLED = 0
      GO TO 1
C
      ENTRY      RMARIN(IJKLIN, NTOTIN,NTOT2N)
C         Initializing routine for RANMAR, may be called before
C         generating pseudorandom numbers with RANMAR. The input
C         values should be in the ranges:  0<=IJKLIN<=900 OOO OOO
C                                          0<=NTOTIN<=999 999 999
C                                          0<=NTOT2N<<999 999 999!
C To get the standard values in Marsaglia's paper, IJKLIN=54217137
C                                            NTOTIN,NTOT2N=0
      IJKL = IJKLIN
      NTOT = MAX(NTOTIN,0)
      NTOT2= MAX(NTOT2N,0)
      KALLED = 1
C          always come here to initialize
    1 CONTINUE
      IJ = IJKL/30082
      KL = IJKL - 30082*IJ
      I = MOD(IJ/177, 177) + 2
      J = MOD(IJ, 177)     + 2
      K = MOD(KL/169, 178) + 1
      L = MOD(KL, 169)
      WRITE(6,'(A,I10,2X,2I10)') ' RANMAR INITIALIZED:',IJKL,NTOT,NTOT2
CCC      PRINT '(A,4I10)', '   I,J,K,L= ',I,J,K,L
      DO 2 II= 1, 97
      S = 0.
      T = .5
      DO 3 JJ= 1, 24
         M = MOD(MOD(I*J,179)*K, 179)
         I = J
         J = K
         K = M
         L = MOD(53*L+1, 169)
         IF (MOD(L*M,64) .GE. 32)  S = S+T
    3    T = 0.5*T
    2 U(II) = S
      TWOM24 = 1.0
      DO 4 I24= 1, 24
    4 TWOM24 = 0.5*TWOM24
      C  =   362436.*TWOM24
      CD =  7654321.*TWOM24
      CM = 16777213.*TWOM24
      I97 = 97
      J97 = 33
C       Complete initialization by skipping
C            (NTOT2*MODCNS + NTOT) random numbers
      DO 45 LOOP2= 1, NTOT2+1
      NOW = MODCNS
      IF (LOOP2 .EQ. NTOT2+1)  NOW=NTOT
      IF (NOW .GT. 0)  THEN
        WRITE(6,'(A,I15)') ' RMARIN SKIPPING OVER ',NOW
       DO 40 IDUM = 1, NTOT
       UNI = U(I97)-U(J97)
       IF (UNI .LT. 0.)  UNI=UNI+1.
       U(I97) = UNI
       I97 = I97-1
       IF (I97 .EQ. 0)  I97=97
       J97 = J97-1
       IF (J97 .EQ. 0)  J97=97
       C = C - CD
       IF (C .LT. 0.)  C=C+CM
   40  CONTINUE
      ENDIF
   45 CONTINUE
      IF (KALLED .EQ. 1)  RETURN
C
C          Normal entry to generate LENV random numbers
   50 CONTINUE
      DO 100 IVEC= 1, LENV
      UNI = U(I97)-U(J97)
      IF (UNI .LT. 0.)  UNI=UNI+1.
      U(I97) = UNI
      I97 = I97-1
      IF (I97 .EQ. 0)  I97=97
      J97 = J97-1
      IF (J97 .EQ. 0)  J97=97
      C = C - CD
      IF (C .LT. 0.)  C=C+CM
      UNI = UNI-C
      IF (UNI .LT. 0.) UNI=UNI+1.
      RVEC(IVEC) = UNI
C             Replace exact zeros by uniform distr. *2**-24
         IF (UNI .EQ. 0.)  THEN
         ZUNI = TWOM24*U(2)
C             An exact zero here is very unlikely, but let's be safe.
         IF (ZUNI .EQ. 0.) ZUNI= TWOM24*TWOM24
         RVEC(IVEC) = ZUNI
         ENDIF
  100 CONTINUE
      NTOT = NTOT + LENV
         IF (NTOT .GE. MODCNS)  THEN
         NTOT2 = NTOT2 + 1
         NTOT = NTOT - MODCNS
         ENDIF
      RETURN
C           Entry to output current status
      ENTRY RMARUT(IJKLUT,NTOTUT,NTOT2T)
      IJKLUT = IJKL
      NTOTUT = NTOT
      NTOT2T = NTOT2
      RETURN
*###] RANMAR:
      END
*###[ RCARRY:
C =====================================================================
      SUBROUTINE RCARRY(RVEC,LENV)
C
C     This version is identical to that in CPC software library
C
C         Add-and-carry random number generator proposed by
C         Marsaglia and Zaman in SIAM J. Scientific and Statistical
C             Computing, to appear probably 1990.
C         modified with enhanced initialization by F. James, 1990
C
C
C  NB! Recently, F. James informed us that there is a slight mistake
C      in this implementation on line 13 and suggested the following
C      change:
C
C      Original line:
C
C      UNI = SEEDS(I24) - SEEDS(J24) - CARRY
C
C      Suggested modification:
C
C      UNI = SEEDS(J24) - SEEDS(I24) - CARRY
C
C      We have also tested this new versionand the test results
C      are available in the hep-lat preprint 9306008.
C
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C!!!  Calling sequences for RCARRY:                                  ++
C!!!      CALL RCARRY (RVEC, LEN)   returns a vector RVEC of LEN     ++
C!!!                   32-bit random floating point numbers between  ++
C!!!                   zero and one.                                 ++
C!!!      CALL RCARGO(INT)     initializes the generator from one    ++
C!!!                   32-bit integer INT                            ++
C!!!      CALL RCARIN(IVEC)    restarts the generator from vector    ++
C!!!                   IVEC of 25 32-bit integers (see RCARUT)       ++
C!!!      CALL RCARUT(IVEC)    outputs the current values of the 25  ++
C!!!                 32-bit integer seeds, to be used for restarting ++
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      double precision RVEC(LENV)
      DIMENSION SEEDS(24), ISEEDS(24), ISDEXT(25)
      PARAMETER (TWOP12=4096.)
      PARAMETER (ITWO24=2**24, ICONS=2147483563)
      SAVE NOTYET, I24, J24, CARRY, SEEDS, TWOM24
      LOGICAL NOTYET
      DATA NOTYET/.TRUE./
      DATA I24,J24,CARRY/24,10,0./
C
C              Default Initialization by Multiplicative Congruential
      IF (NOTYET) THEN
         NOTYET = .FALSE.
         JSEED = 314159265
         WRITE(6,'(A,I12)') ' RCARRY DEFAULT INITIALIZATION: ',JSEED
            TWOM24 = 1.
         DO 25 I= 1, 24
            TWOM24 = TWOM24 * 0.5
         K = JSEED/53668
         JSEED = 40014*(JSEED-K*53668) -K*12211
         IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
         ISEEDS(I) = MOD(JSEED,ITWO24)
   25    CONTINUE
         DO 50 I= 1,24
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24
   50    CONTINUE
         I24 = 24
         J24 = 10
         CARRY = 0.
         IF (SEEDS(24) .LT. SEEDS(14)) CARRY = TWOM24
      ENDIF
C
C          The Generator proper: "Subtract-with-borrow",
C          as proposed by Marsaglia and Zaman,
C          Florida State University, March, 1989
C
      DO 100 IVEC= 1, LENV
      UNI = SEEDS(I24) - SEEDS(J24) - CARRY
      IF (UNI .LT. 0.)  THEN
         UNI = UNI + 1.0
         CARRY = TWOM24
      ELSE
         CARRY = 0.
      ENDIF
      SEEDS(I24) = UNI
      I24 = I24 - 1
      IF (I24 .EQ. 0)  I24 = 24
      J24 = J24 - 1
      IF (J24 .EQ. 0)  J24 = 24
      RVEC(IVEC) = UNI
  100 CONTINUE
      RETURN
C           Entry to input and float integer seeds from previous run
      ENTRY RCARIN(ISDEXT)
         TWOM24 = 1.
         DO 195 I= 1, 24
  195    TWOM24 = TWOM24 * 0.5
      WRITE(6,'(A)') ' FULL INITIALIZATION OF RCARRY WITH 25 INTEGERS:'
      WRITE(6,'(5X,5I12)') ISDEXT
      DO 200 I= 1, 24
      SEEDS(I) = REAL(ISDEXT(I))*TWOM24
  200 CONTINUE
      CARRY = REAL(MOD(ISDEXT(25),10))*TWOM24
      ISD = ISDEXT(25)/10
      I24 = MOD(ISD,100)
      ISD = ISD/100
      J24 = ISD
      RETURN
C                    Entry to ouput seeds as integers
      ENTRY RCARUT(ISDEXT)
      DO 300 I= 1, 24
         ISDEXT(I) = INT(SEEDS(I)*TWOP12*TWOP12)
  300 CONTINUE
      ICARRY = 0
      IF (CARRY .GT. 0.)  ICARRY = 1
      ISDEXT(25) = 1000*J24 + 10*I24 + ICARRY
      RETURN
C                    Entry to initialize from one integer
      ENTRY RCARGO(INSEED)
      JSEED = INSEED
      WRITE(6,'(A,I12)') ' RCARRY INITIALIZED FROM SEED ',INSEED
C      TWOM24 = 1.
         DO 325 I= 1, 24
           TWOM24 = TWOM24 * 0.5
         K = JSEED/53668
         JSEED = 40014*(JSEED-K*53668) -K*12211
         IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
         ISEEDS(I) = MOD(JSEED,ITWO24)
  325    CONTINUE
         DO 350 I= 1,24
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24
  350    CONTINUE
         I24 = 24
         J24 = 10
         CARRY = 0.
         IF (SEEDS(24) .LT. SEEDS(14)) CARRY = TWOM24
      RETURN
*###] RCARRY:
      END
*###[ r250: random number generator R250
	double precision function R250(dummy)
***#[*comment:***********************************************************
*									*
*	Random number generator according to				*
*	S.Kirkpatrick & E.P.Stoll,,J.Comp.Phys.40(1981)517		*
*	initialized with a congruent generator.				*
*	See also K.Kankaala,CSC Res.Rep.R03/93 (PhD thesis)		*
*	I hope I made no errors.   GJvO 11-jan-94			*
*									*
***#]*comment:***********************************************************
*  #[ declarations:
	implicit none
*
*	argument
*
	integer dummy
*
*	local variables
*
	integer*4 ix(0:249),ip,init,i,ggli
	double precision range,resetr250,xdum
	save ix,ip,init,range
	integer resetggli
	external resetggli
*
*	data
*
	data init /0/
	data ip /0/
*
*  #] declarations:
*  #[ reset:
*	
	goto 1
	entry resetr250(dummy)
	init = 0
	ip = 0
	xdum = resetggli(dummy)
	return
    1	continue
*
*  #] reset:
*  #[ init:
*
	if ( init.eq.0 ) then
	    init = 1
	    do 10 i=0,249
		ix(i) = ggli(i)
   10	    continue
	    range = 1/dble(2)**31
	endif
*
*  #] init:
*  #[ work:
*
	ix(ip) = ieor( ix(ip), ix(mod(ip+147,250)) )
	R250 = ix(ip)*range
	ip = ip + 1
	if ( ip.gt.249 ) ip = 0
*
*  #] work:
*###] r250:
	end
*###[ ggli: random number generator GGL
	integer*4 function ggli(dummy)
***#[*comment:***********************************************************
*									*
*	The GGL linear congruent MC generator, use only to initialize	*
*	R250 as the Marsiglia planes are quite bad and the period short	*
*									*
***#]*comment:***********************************************************
*  #[ declarations:
	implicit none
*
*	arguments
*
	integer dummy
*
*	local variables
*
	double precision x,m
	save x,m
	integer resetggli
*
*	data
*
	data x/667790.d0/
	data m/2147483647.d0/
*
*  #] declarations:
*  #[ reset:
*	
	goto 1
	entry resetggli(dummy)
	x = 667790.d0
	resetggli = 0
	return
    1	continue
*
*  #] reset:
*  #[ work:
*
	x = mod(16807.d0*x,m)
	ggli = x
*
*  #] work:
*###] ggli:
	end
*###[ ranlux:
      SUBROUTINE RANLUX(RVEC,LENV)
C         Subtract-and-borrow random number generator proposed by
C         Marsaglia and Zaman, implemented by F. James with the name
C         RCARRY in 1991, and later improved by Martin Luescher
C         in 1993 to produce "Luxury Pseudorandom Numbers".
C     Fortran 77 coded by F. James, 1993
C
C       references:
C  M. Luscher, Computer Physics Communications  79 (1994) 100
C  F. James, Computer Physics Communications 79 (1994) 111
C
C   LUXURY LEVELS.
C   ------ ------      The available luxury levels are:
C
C  level 0  (p=24): equivalent to the original RCARRY of Marsaglia
C           and Zaman, very long period, but fails many tests.
C  level 1  (p=48): considerable improvement in quality over level 0,
C           now passes the gap test, but still fails spectral test.
C  level 2  (p=97): passes all known tests, but theoretically still
C           defective.
C  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible
C           correlations have very small chance of being observed.
C  level 4  (p=389): highest possible luxury, all 24 bits chaotic.
C
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C!!!  Calling sequences for RANLUX:                                  ++
C!!!      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++
C!!!                   32-bit random floating point numbers between  ++
C!!!                   zero (not included) and one (also not incl.). ++
C!!!      CALL RLUXGO(LUX,INT,K1,K2) initializes the generator from  ++
C!!!               one 32-bit integer INT and sets Luxury Level LUX  ++
C!!!               which is integer between zero and MAXLEV, or if   ++
C!!!               LUX .GT. 24, it sets p=LUX directly.  K1 and K2   ++
C!!!               should be set to zero unless restarting at a break++
C!!!               point given by output of RLUXAT (see RLUXAT).     ++
C!!!      CALL RLUXAT(LUX,INT,K1,K2) gets the values of four integers++
C!!!               which can be used to restart the RANLUX generator ++
C!!!               at the current point by calling RLUXGO.  K1 and K2++
C!!!               specify how many numbers were generated since the ++
C!!!               initialization with LUX and INT.  The restarting  ++
C!!!               skips over  K1+K2*E9   numbers, so it can be long.++
C!!!   A more efficient but less convenient way of restarting is by: ++
C!!!      CALL RLUXIN(ISVEC)    restarts the generator from vector   ++
C!!!                   ISVEC of 25 32-bit integers (see RLUXUT)      ++
C!!!      CALL RLUXUT(ISVEC)    outputs the current values of the 25 ++
C!!!                 32-bit integer seeds, to be used for restarting ++
C!!!      ISVEC must be dimensioned 25 in the calling program        ++
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DIMENSION RVEC(LENV)
      DIMENSION SEEDS(24), ISEEDS(24), ISDEXT(25)
      PARAMETER (MAXLEV=4, LXDFLT=3)
      DIMENSION NDSKIP(0:MAXLEV)
      DIMENSION NEXT(24)
      PARAMETER (TWOP12=4096., IGIGA=1000000000,JSDFLT=314159265)
      PARAMETER (ITWO24=2**24, ICONS=2147483563)
      SAVE NOTYET, I24, J24, CARRY, SEEDS, TWOM24, TWOM12, LUXLEV
      SAVE NSKIP, NDSKIP, IN24, NEXT, KOUNT, MKOUNT, INSEED
      INTEGER LUXLEV
      LOGICAL NOTYET
      DATA NOTYET, LUXLEV, IN24, KOUNT, MKOUNT /.TRUE., LXDFLT, 0,0,0/
      DATA I24,J24,CARRY/24,10,0./
C                               default
C  Luxury Level   0     1     2   *3*    4
      DATA NDSKIP/0,   24,   73,  199,  365 /
Corresponds to p=24    48    97   223   389
C     time factor 1     2     3     6    10   on slow workstation
C                 1    1.5    2     3     5   on fast mainframe
C
C  NOTYET is .TRUE. if no initialization has been performed yet.
C              Default Initialization by Multiplicative Congruential
      IF (NOTYET) THEN
         NOTYET = .FALSE.
         JSEED = JSDFLT
         INSEED = JSEED
         WRITE(6,'(A,I12)') ' RANLUX DEFAULT INITIALIZATION: ',JSEED
         LUXLEV = LXDFLT
         NSKIP = NDSKIP(LUXLEV)
         LP = NSKIP + 24
         IN24 = 0
         KOUNT = 0
         MKOUNT = 0
         WRITE(6,'(A,I2,A,I4)')  ' RANLUX DEFAULT LUXURY LEVEL =  ',
     +        LUXLEV,'      p =',LP
            TWOM24 = 1.
         DO 25 I= 1, 24
            TWOM24 = TWOM24 * 0.5
         K = JSEED/53668
         JSEED = 40014*(JSEED-K*53668) -K*12211
         IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
         ISEEDS(I) = MOD(JSEED,ITWO24)
   25    CONTINUE
         TWOM12 = TWOM24 * 4096.
         DO 50 I= 1,24
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24
         NEXT(I) = I-1
   50    CONTINUE
         NEXT(1) = 24
         I24 = 24
         J24 = 10
         CARRY = 0.
         IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24
      ENDIF
C
C          The Generator proper: "Subtract-with-borrow",
C          as proposed by Marsaglia and Zaman,
C          Florida State University, March, 1989
C
      DO 100 IVEC= 1, LENV
      UNI = SEEDS(J24) - SEEDS(I24) - CARRY
      IF (UNI .LT. 0.)  THEN
         UNI = UNI + 1.0
         CARRY = TWOM24
      ELSE
         CARRY = 0.
      ENDIF
      SEEDS(I24) = UNI
      I24 = NEXT(I24)
      J24 = NEXT(J24)
      RVEC(IVEC) = UNI
C  small numbers (with less than 12 "significant" bits) are "padded".
      IF (UNI .LT. TWOM12)  THEN
         RVEC(IVEC) = RVEC(IVEC) + TWOM24*SEEDS(J24)
C        and zero is forbidden in case someone takes a logarithm
         IF (RVEC(IVEC) .EQ. 0.)  RVEC(IVEC) = TWOM24*TWOM24
      ENDIF
C        Skipping to luxury.  As proposed by Martin Luscher.
      IN24 = IN24 + 1
      IF (IN24 .EQ. 24)  THEN
         IN24 = 0
         KOUNT = KOUNT + NSKIP
         DO 90 ISK= 1, NSKIP
         UNI = SEEDS(J24) - SEEDS(I24) - CARRY
         IF (UNI .LT. 0.)  THEN
            UNI = UNI + 1.0
            CARRY = TWOM24
         ELSE
            CARRY = 0.
         ENDIF
         SEEDS(I24) = UNI
         I24 = NEXT(I24)
         J24 = NEXT(J24)
   90    CONTINUE
      ENDIF
  100 CONTINUE
      KOUNT = KOUNT + LENV
      IF (KOUNT .GE. IGIGA)  THEN
         MKOUNT = MKOUNT + 1
         KOUNT = KOUNT - IGIGA
      ENDIF
      RETURN
C
C           Entry to input and float integer seeds from previous run
      ENTRY RLUXIN(ISDEXT)
*     IF block added by Phillip Helbig after correpondence with James
      IF (NOTYET) THEN
         WRITE(6,'(A)')  ' PROPER RESULTS ONLY WITH INITIALISATION FROM
     $25 INTEGERS OBTAINED WITH RLUXUT'
         NOTYET = .FALSE.
      ENDIF
         TWOM24 = 1.
         DO 195 I= 1, 24
         NEXT(I) = I-1
  195    TWOM24 = TWOM24 * 0.5
         NEXT(1) = 24
         TWOM12 = TWOM24 * 4096.
      WRITE(6,'(A)') ' FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:'
      WRITE(6,'(5X,5I12)') ISDEXT
      DO 200 I= 1, 24
      SEEDS(I) = REAL(ISDEXT(I))*TWOM24
  200 CONTINUE
      CARRY = 0.
      IF (ISDEXT(25) .LT. 0)  CARRY = TWOM24
      ISD = IABS(ISDEXT(25))
      I24 = MOD(ISD,100)
      ISD = ISD/100
      J24 = MOD(ISD,100)
      ISD = ISD/100
      IN24 = MOD(ISD,100)
      ISD = ISD/100
      LUXLEV = ISD
        IF (LUXLEV .LE. MAXLEV) THEN
          NSKIP = NDSKIP(LUXLEV)
          WRITE (6,'(A,I2)') ' RANLUX LUXURY LEVEL SET BY RLUXIN TO: ',
     +                         LUXLEV
        ELSE  IF (LUXLEV .GE. 24) THEN
          NSKIP = LUXLEV - 24
          WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXIN TO:',LUXLEV
        ELSE
          NSKIP = NDSKIP(MAXLEV)
          WRITE (6,'(A,I5)') ' RANLUX ILLEGAL LUXURY RLUXIN: ',LUXLEV
          LUXLEV = MAXLEV
        ENDIF
      INSEED = -1
      RETURN
C
C                    Entry to ouput seeds as integers
      ENTRY RLUXUT(ISDEXT)
      DO 300 I= 1, 24
         ISDEXT(I) = INT(SEEDS(I)*TWOP12*TWOP12)
  300 CONTINUE
      ISDEXT(25) = I24 + 100*J24 + 10000*IN24 + 1000000*LUXLEV
      IF (CARRY .GT. 0.)  ISDEXT(25) = -ISDEXT(25)
      RETURN
C
C                    Entry to output the "convenient" restart point
      ENTRY RLUXAT(LOUT,INOUT,K1,K2)
      LOUT = LUXLEV
      INOUT = INSEED
      K1 = KOUNT
      K2 = MKOUNT
      RETURN
C
C                    Entry to initialize from one or three integers
      ENTRY RLUXGO(LUX,INS,K1,K2)
         IF (LUX .LT. 0) THEN
            LUXLEV = LXDFLT
         ELSE IF (LUX .LE. MAXLEV) THEN
            LUXLEV = LUX
         ELSE IF (LUX .LT. 24 .OR. LUX .GT. 2000) THEN
            LUXLEV = MAXLEV
            WRITE (6,'(A,I7)') ' RANLUX ILLEGAL LUXURY RLUXGO: ',LUX
         ELSE
            LUXLEV = LUX
            DO 310 ILX= 0, MAXLEV
              IF (LUX .EQ. NDSKIP(ILX)+24)  LUXLEV = ILX
  310       CONTINUE
         ENDIF
      IF (LUXLEV .LE. MAXLEV)  THEN
         NSKIP = NDSKIP(LUXLEV)
         WRITE(6,'(A,I2,A,I4)') ' RANLUX LUXURY LEVEL SET BY RLUXGO :',
     +        LUXLEV,'     P=', NSKIP+24
      ELSE
          NSKIP = LUXLEV - 24
          WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXGO TO:',LUXLEV
      ENDIF
      IN24 = 0
      IF (INS .LT. 0)  WRITE (6,'(A)')
     +   ' Illegal initialization by RLUXGO, negative input seed'
      IF (INS .GT. 0)  THEN
        JSEED = INS
        WRITE(6,'(A,3I12)') ' RANLUX INITIALIZED BY RLUXGO FROM SEEDS',
     +      JSEED, K1,K2
      ELSE
        JSEED = JSDFLT
        WRITE(6,'(A)')' RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED'
      ENDIF
      INSEED = JSEED
      NOTYET = .FALSE.
      TWOM24 = 1.
         DO 325 I= 1, 24
           TWOM24 = TWOM24 * 0.5
         K = JSEED/53668
         JSEED = 40014*(JSEED-K*53668) -K*12211
         IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
         ISEEDS(I) = MOD(JSEED,ITWO24)
  325    CONTINUE
      TWOM12 = TWOM24 * 4096.
         DO 350 I= 1,24
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24
         NEXT(I) = I-1
  350    CONTINUE
      NEXT(1) = 24
      I24 = 24
      J24 = 10
      CARRY = 0.
      IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24
C        If restarting at a break point, skip K1 + IGIGA*K2
C        Note that this is the number of numbers delivered to
C        the user PLUS the number skipped (if luxury .GT. 0).
      KOUNT = K1
      MKOUNT = K2
      IF (K1+K2 .NE. 0)  THEN
        DO 500 IOUTER= 1, K2+1
          INNER = IGIGA
          IF (IOUTER .EQ. K2+1)  INNER = K1
          DO 450 ISK= 1, INNER
            UNI = SEEDS(J24) - SEEDS(I24) - CARRY
            IF (UNI .LT. 0.)  THEN
               UNI = UNI + 1.0
               CARRY = TWOM24
            ELSE
               CARRY = 0.
            ENDIF
            SEEDS(I24) = UNI
            I24 = NEXT(I24)
            J24 = NEXT(J24)
  450     CONTINUE
  500   CONTINUE
C         Get the right value of IN24 by direct calculation
        IN24 = MOD(KOUNT, NSKIP+24)
        IF (MKOUNT .GT. 0)  THEN
           IZIP = MOD(IGIGA, NSKIP+24)
           IZIP2 = MKOUNT*IZIP + IN24
           IN24 = MOD(IZIP2, NSKIP+24)
        ENDIF
C       Now IN24 had better be between zero and 23 inclusive
        IF (IN24 .GT. 23) THEN
           WRITE (6,'(A/A,3I11,A,I5)')
     +    '  Error in RESTARTING with RLUXGO:','  The values', INS,
     +     K1, K2, ' cannot occur at luxury level', LUXLEV
           IN24 = 0
        ENDIF
      ENDIF
      RETURN
      END
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*      PROGRAM LUXTST
*C         Exercise for the RANLUX Pseudorandom number generator.
*C
*      DIMENSION RVEC(1000)
*      DIMENSION ISDEXT(25)
*C
*C         check that we get the right numbers (machine-indep.)
*      WRITE (6,'(/A)')  '  CALL RANLUX(RVEC,100)'
*      CALL RANLUX(RVEC,100)
*      WRITE (6,'(A/9X,5F12.8)') ' RANLUX default numbers   1-  5:',
*     +    (RVEC(L),L=1,5)
*      CALL RANLUX(RVEC,100)
*      WRITE (6,'(A/9X,5F12.8)') ' RANLUX default numbers 101-105:',
*     +    (RVEC(L),L=1,5)
*C
*      WRITE (6,'(/A)')  ' CALL RLUXGO(0,0,0,0)'
*      CALL RLUXGO(0,0,0,0)
*      CALL RANLUX(RVEC,100)
*      WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury level 0,   1-  5:',
*     +    (RVEC(L),L=1,5)
*      CALL RANLUX(RVEC,100)
*      WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury level 0, 101-105:',
*     +    (RVEC(L),L=1,5)
*C
*      WRITE (6,'(/A)')  '   CALL RLUXGO(389,1,0,0)'
*      CALL RLUXGO(389,1,0,0)
*      CALL RANLUX(RVEC,100)
*      WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury p=389,   1-  5:',
*     +    (RVEC(L),L=1,5)
*      CALL RANLUX(RVEC,100)
*      WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury p=389, 101-105:',
*     +    (RVEC(L),L=1,5)
*C
*      WRITE (6,'(/A)')  '  CALL RLUXGO(75,0,0,0)'
*      CALL RLUXGO(75,0,0,0)
*      CALL RANLUX(RVEC,100)
*      WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury p= 75,   1-  5:',
*     +    (RVEC(L),L=1,5)
*      CALL RANLUX(RVEC,100)
*      WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury p= 75, 101-105:',
*     +    (RVEC(L),L=1,5)
*C
*      WRITE (6,'(/A)')  '  test restarting from the full vector'
*      CALL RLUXUT(ISDEXT)
*      WRITE (6,'(/A/(1X,5I14))') '  current RANLUX status saved:',ISDEXT
*      CALL RANLUX(RVEC,100)
*      WRITE (6,'(A/9X,5F12.8)') ' RANLUX numbers 1- 5:',
*     +    (RVEC(L),L=1,5)
*      CALL RANLUX(RVEC,100)
*      WRITE (6,'(A/9X,5F12.8)') ' RANLUX numbers 101-105:',
*     +    (RVEC(L),L=1,5)
*C
*      WRITE (6,'(/A)')   '   previous RANLUX status will be restored'
*      CALL RLUXIN(ISDEXT)
*      CALL RANLUX(RVEC,100)
*      WRITE (6,'(A/9X,5F12.8)') ' RANLUX numbers 1- 5:',
*     +    (RVEC(L),L=1,5)
*      CALL RANLUX(RVEC,100)
*      WRITE (6,'(A/9X,5F12.8)') ' RANLUX numbers 101-105:',
*     +    (RVEC(L),L=1,5)
*C
*      WRITE (6,'(/A)')  '     test the restarting by skipping'
*      CALL RLUXGO(4,7674985,0,0)
*      CALL RLUXAT(I1,I2,I3,I4)
*      WRITE (6,'(A,4I10)')  '  RLUXAT values =',I1,I2,I3,I4
*      DO 150 LI= 1, 10
*  150 CALL RANLUX(RVEC,1000)
*      CALL RLUXAT(I1,I2,I3,I4)
*      WRITE (6,'(A,4I10)')  '  RLUXAT values =',I1,I2,I3,I4
*      CALL RANLUX(RVEC,200)
*      WRITE (6,'(A,2F10.6)')  '  Next and 200th numbers are:',
*     +                             RVEC(1), RVEC(200)
*      CALL RLUXGO(I1,I2,I3,I4)
*      CALL RANLUX(RVEC,200)
*      WRITE (6,'(A,2F10.6)')  '  Next and 200th numbers are:',
*     +                             RVEC(1), RVEC(200)
*C
*      WRITE (6,'(/A)') ' The following should provoke an error message'
*      CALL RLUXGO(4,11111,31,0)
*      STOP
*C
*C   OUTPUT FROM THE ABOVE TEST PROGRAM SHOULD BE:
*C   --------------------------------------------
*C  CALL RANLUX(RVEC,100)
*C RANLUX DEFAULT INITIALIZATION:    314159265
*C RANLUX DEFAULT LUXURY LEVEL =   3      p = 223
*C RANLUX default numbers   1-  5:
*C           0.53981817  0.76155043  0.06029940  0.79600263  0.30631220
*C RANLUX default numbers 101-105:
*C           0.43156743  0.03774416  0.24897110  0.00147784  0.90274453
*C
*C  CALL RLUXGO(0,0,0,0)
*C RANLUX LUXURY LEVEL SET BY RLUXGO : 0     P=  24
*C RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED
*C RANLUX luxury level 0,   1-  5:
*C           0.53981817  0.76155043  0.06029940  0.79600263  0.30631220
*C RANLUX luxury level 0, 101-105:
*C           0.41538775  0.05330932  0.58195311  0.91397446  0.67034441
*C
*C   CALL RLUXGO(389,1,0,0)
*C RANLUX LUXURY LEVEL SET BY RLUXGO : 4     P= 389
*C RANLUX INITIALIZED BY RLUXGO FROM SEEDS           1           0           0
*C RANLUX luxury p=389,   1-  5:
*C           0.94589490  0.47347850  0.95152789  0.42971975  0.09127384
*C RANLUX luxury p=389, 101-105:
*C           0.02618265  0.03775346  0.97274780  0.13302165  0.43126065
*C
*C  CALL RLUXGO(75,0,0,0)
*C RANLUX P-VALUE SET BY RLUXGO TO:   75
*C RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED
*C RANLUX luxury p= 75,   1-  5:
*C           0.53981817  0.76155043  0.06029940  0.79600263  0.30631220
*C RANLUX luxury p= 75, 101-105:
*C           0.25600731  0.23443210  0.59164381  0.59035838  0.07011414
*C
*C  test restarting from the full vector
*C
*C  current RANLUX status saved:
*C       16156027      16534309      15243811       2751687       6002207
*C        7979506       1301976       4567313       4305996       5872599
*C       12003090       2146823      12606367       4111505       5979640
*C       12739666      10489318      14036909      11729352       8061448
*C        7832659       6069758       3197719       1832730      75080216
*C RANLUX numbers 1- 5:
*C           0.22617835  0.60655993  0.86417443  0.43920082  0.23382509
*C RANLUX numbers 101-105:
*C           0.08107197  0.21466845  0.84856731  0.94078046  0.85626233
*C
*C   previous RANLUX status will be restored
*C FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:
*C         16156027    16534309    15243811     2751687     6002207
*C          7979506     1301976     4567313     4305996     5872599
*C         12003090     2146823    12606367     4111505     5979640
*C         12739666    10489318    14036909    11729352     8061448
*C          7832659     6069758     3197719     1832730    75080216
*C RANLUX P-VALUE SET BY RLUXIN TO:   75
*C RANLUX numbers 1- 5:
*C           0.22617835  0.60655993  0.86417443  0.43920082  0.23382509
*C RANLUX numbers 101-105:
*C           0.08107197  0.21466845  0.84856731  0.94078046  0.85626233
*C
*C     test the restarting by skipping
*C RANLUX LUXURY LEVEL SET BY RLUXGO : 4     P= 389
*C RANLUX INITIALIZED BY RLUXGO FROM SEEDS     7674985           0           0
*C  RLUXAT values =         4   7674985         0         0
*C  RLUXAT values =         4   7674985    161840         0
*C  Next and 200th numbers are:  0.019648  0.590586
*C RANLUX LUXURY LEVEL SET BY RLUXGO : 4     P= 389
*C RANLUX INITIALIZED BY RLUXGO FROM SEEDS     7674985      161840           0
*C  Next and 200th numbers are:  0.019648  0.590586
*C
*C The following should provoke an error message
*C RANLUX LUXURY LEVEL SET BY RLUXGO : 4     P= 389
*C RANLUX INITIALIZED BY RLUXGO FROM SEEDS       11111          31           0
*C  Error in RESTARTING with RLUXGO:
*C  The values      11111         31          0 cannot occur at luxury level    4
*      END
*###] ranlux:
*###[ richtmeyer quasi-random:
	DOUBLE PRECISION function richtm(idummy)
	implicit none
	integer idummy,init,j,k,n
	parameter(n=6)
	save init
	DOUBLE PRECISION s(n),sum
	data init/0/
	if ( init.eq.0 ) then
	    init = 1
	    s(1) = sqrt(2d0) - 1
	    s(2) = sqrt(3d0) - 1
	    s(3) = sqrt(5d0) - 2
	    s(4) = sqrt(7d0) - 2
	    s(5) = sqrt(11d0)- 3
	    s(6) = sqrt(13d0)- 3
	endif
	do j = 1,n
	    sum = 0
	    do k = 1,n
		sum = sum + s(k)
	    enddo
	    s(j) = mod(sum,1d0)
	enddo
	richtm = s(1)
*###] richtmeyer quasi-random:
      end


