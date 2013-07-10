	program grads2extra
*  #[ comment:
*
*	Quick hack to convert a GrADS data file to an extra file
*	with NO time, field or level information
*
*  #] comment:
*  #[ declarations:
	implicit none
#include "recfac.h"
	integer NMAX
	parameter (NMAX=20000)
	integer i,j,nsize,date,type
	real field(NMAX)
	character size*12, infile*128, outfile*128
	integer iargc
*  #] declarations:
*  #[ arguments:
	if ( iargc().lt.3 ) then
	    print *,'usage: grads2extra nsize infile outfile '//
     +		  '[date [type]]'
	    stop
	endif
	date = 0
	type = 0
	call getarg(1,size)
	read(size,'(i12)') nsize
	call getarg(2,infile)
	call getarg(3,outfile)
	if ( iargc().ge.4 ) then
	    call getarg(4,size)
	    read(size,*) date
	endif
	if ( iargc().ge.5 ) then
	    call getarg(5,size)
	    read(size,*) type
	endif
*  #] arguments:
*  #[ convert:
*
*	basically, all we do is to convert from an unformatted file to 
*	a formatted file...
*
	open(1,file=infile,status='old',form='unformatted',
     +		access='direct',recl=recfac*nsize)
	open(2,file=outfile,status='new',form='unformatted')
	do i=1,10000
	    read(1,rec=i,end=200,err=200) (field(j),j=1,nsize)
	    write(2) date,type,1,nsize
	    write (*,*) date,type,1,nsize
	    write(2) (field(j),j=1,nsize)
	    date = date + 1
	    if ( date.gt.150000 .and. date.lt.202000 ) then
		if ( mod(date,100).ge.13 ) then
		    date = 100*(date/100 + 1) + mod(date,100) - 12
		endif
	    endif
	enddo
	print *,'more than 10000 records???'
  200	continue
*  #] convert:
	end
