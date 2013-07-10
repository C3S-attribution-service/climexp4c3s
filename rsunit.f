*  #[ rsunit:
	subroutine rsunit(irsunit)
*
*       find a free unit number below 100
*
	implicit none
	integer irsunit
	logical lopen
	do irsunit=99,10,-1
	    inquire(irsunit,opened=lopen)
	    if ( .not.lopen ) goto 20
	enddo
	print '(a)','rsunit: error: no free units under 100!'
	call abort
   20	continue
*  #] rsunit:
	end
