*###[ ffsort:
	subroutine ffsort(a,ii,nn)
***#[*comment:***********************************************************
*									*
*	Sort the array a(nn): give the position of the smallest element	*
*	in ii(1), ..., largest in ii(nn).  I use a fancy merge-sort	*
*	algorithm which is probably not the samrtest thing to do with	*
*	the small arrays for which it is used, but it was fun to program*
*	To extend to larger arrays: just change 1024 to some power of 2	*
*									*
*	Input:	a	real(nn)		array			*
*		nn	integer						*
*	Output:	ii	integer(nn)	a(ii(1))<=a(ii(2))<=.<=a(ii(nn))*
*									*
***#]*comment:***********************************************************
*  #[ declarations:
	implicit none
*
*	arguments
*
	integer nn,ii(nn)
	real a(nn)
*
*	local variables
*
	integer nmax
	parameter (nmax=65536)
	integer i,j,k,jj(nmax,2),h,j12,j21,l,m,n,o
*
*	common
*
***	include 'ff.h'
*
*  #] declarations:
*  #[ work:
	if ( nn.gt.nmax ) then
	    print  *,'ffsort: can only sort up to ',nmax
     +		  ,' elements, not ',nn
	    stop
	endif
	do 10 i=1,nn
	    jj(i,1) = i
   10	continue
	j12 = 1
	j21 = 2
*
*	do the first sweep faster
*
	do 15 i=1,nn-1,2
	    if ( a(jj(i,j12)) .le. a(jj(i+1,j12)) ) then
		jj(i,j21) = jj(i,j12)
		jj(i+1,j21) = jj(i+1,j12)
	    else
		jj(i,j21) = jj(i+1,j12)
		jj(i+1,j21) = jj(i,j12)
	    endif
   15	continue
	if ( mod(nn,2).ne.0 ) jj(nn,j21) = jj(nn,j12)
	o   = j12
	j12 = j21
	j21 = o
*
*	and do the other sweeps (works also for k=1,...)
*
	do 100 k=2,nint(log(dble(nmax))/log(dble(2)))
	    h = 2**k
	    do 90 j=1,nn,h
		l = j
		n = j
		m = j+h/2
		if ( m.gt.nn ) then
		    do 17 o=j,nn
			jj(o,j21) = jj(o,j12)
   17		    continue
		    goto 90
		endif
		do 20 i=1,2*nmax
		    if ( a(jj(l,j12)) .le. a(jj(m,j12)) ) then
			jj(n,j21) = jj(l,j12)
			l = l+1
			n = n+1
			if ( l.ge.j+h/2 ) then
			    do 18 o=m,min(j+h-1,nn)
				jj(n,j21) = jj(o,j12)
				n = n+1
   18			    continue
			    goto 21
			endif
		    else
			jj(n,j21) = jj(m,j12)
			m = m+1
			n = n+1
			if ( m.ge.j+h .or. m.gt.nn ) then
			    do 19 o=l,j+h/2-1
				jj(n,j21) = jj(o,j12)
				n = n+1
   19			    continue
			    goto 21
			endif
		    endif
   20		continue
   21		continue
		if ( n.ne.j+h .and. n.ne.nn+1 ) print *,'n wrong: ',n
   90	    continue
	    o   = j12
	    j12 = j21
	    j21 = o
	    if ( h.ge.nn ) goto 900
  100	continue
  900	continue
	do 901 i=1,nn
	    ii(i) = jj(i,j12)
  901	continue
*  #] work:
*  #[ debug output:
*	if ( lwrite ) then
*	    print *,'This should be sorted:'
*	    do 910 i=1,nn
*		print '(i5,f20.8)',ii(i),a(ii(i))
*  910	    continue
*	endif
*  #] debug output:
*###] ffsort:
	end
