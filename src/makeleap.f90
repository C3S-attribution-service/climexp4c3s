program makepleap
	!	
	! make a time series that denotes how many years have elapsed since the last leap year
	!
	implicit none
	integer yrbeg,yrend
	parameter (yrbeg=1600,yrend=2300)
	integer mo,yr
	real delta,odelta
	integer,external :: leap

	print '(a)','# delta [dy] difference between tropical calendar and civiv calendar'
	print '(a)','# to investigate spurious trends in spring and autumn'
	delta = 0
	do yr=yrbeg,yrend
		delta = delta + 0.25 - 0.0075
		if ( leap(yr).eq.2 ) then
			odelta = delta
			delta = delta - 1
		end if
		print '(i4,12f8.4)',yr,(-odelta,mo=1,2),(-delta,mo=3,12)
		odelta = delta
	end do
end program