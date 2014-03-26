	integer function leap(yr)
!
!	returns 1 if the year yr is not a leap year
!	returns 2 if it is a leap year
!	Gregorian calendar assumed throughout
!	
	implicit none
	integer yr
	if ( mod(yr,4).ne.0 ) then
	    leap = 1
	elseif ( mod(yr,100).ne.0 ) then
	    leap = 2
	elseif ( mod(yr,400).ne.0 ) then
	    leap = 1
	else	
	    leap = 2
	endif
	end


