	integer function leap(yr)
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


