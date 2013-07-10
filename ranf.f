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
	    ranf = R250(dummy)
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
	xdum = resetr250(dummy)
*
*  #] switchyard:
*###] resetranf:
	end
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
	integer ix(0:249),ip,init,i,ggli
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
	integer function ggli(dummy)
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
