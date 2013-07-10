	program invalidfield
*
*	make a field full of 3e33 to fill up the prcp database
*
	implicit none
	integer i,n,m
#include "recfac.h"
	character *80 string
	integer iargc
*
	call getarg(1,string)
	read(string,*) n
	call getarg(2,string)
	read(string,*) m
*
	open(1,file='invalid.dat',form='unformatted',access='direct',
     +		recl=recfa4*n*m)
	write(1,rec=1) (3e33,i=1,n*m)		
	close(1)
	end
