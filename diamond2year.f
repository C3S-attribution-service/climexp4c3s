	program diamond2year
*
*	replace the strings '\Diamond' by the year
*
	implicit none
	integer i,year,month,oldmonth,l,init,n,m
	real x(10)
	character line*200,line2*180,file*1000
	logical lmonth
	integer iargc,llen,getnumwords
*
	if ( iargc().ne.2 ) then
	    print *,'usage: diamond2year psfile datfile > outfile'
	    stop
	endif
	call getarg(1,file)
	open(1,file=file,status='old')
	call getarg(2,file)
	open(2,file=file,status='old')
*
	init = 0
	oldmonth = 0
	lmonth = .FALSE.
  100	continue
	read(1,'(a)',end=800,err=900) line
	i = index(line,' Pls    ')
	l = llen(line)
	if ( i.ne.0 ) then
	    if ( init.eq.0 ) then
		init = 1
		write(*,'(a)')
     +		      '(Helvetica) findfont 100 scalefont setfont'
	    endif
  110	    continue
	    read(2,'(a)') line2
	    if ( line2(1:1).eq.'#' .or. line2.eq.' ' ) goto 110
	    n = getnumwords(line2)
	    read(line2,*,err=901) (x(m),m=1,n-4),year,month
	    if ( lmonth .or. oldmonth.eq.0 ) then
		if ( month.lt.10 ) then
		    write(*,'(2a,i1,a,i4,a)') line(1:i-1),
     +			  ' M (',month,'.',year,') Cshow'
		elseif ( month.lt.100 ) then
		    write(*,'(2a,i2,a,i4,a)') line(1:i-1),
     +			  ' M (',month,'.',year,') Cshow'
		else
		    write(*,'(2a,i3,a,i4,a)') line(1:i-1),
     +			  ' M (',month,'.',year,') Cshow'
		endif
	    else
		write(*,'(2a,i4,a)') line(1:i-1),' M (',year,') Cshow'
	    endif
	    if ( month.ne.oldmonth ) then
		if ( oldmonth.gt.0 ) lmonth = .TRUE.
		oldmonth = month
	    endif
	else
	    if ( init.eq.1 ) then
		init = 0
		write(*,'(a)')
     +		      '(Helvetica) findfont 140 scalefont setfont'
	    endif
	    write(*,'(a)') line(:l)
	endif
	goto 100
*
  800	continue
	goto 999
  900	print *,'error reading input'
	print *,line
	call abort
  901	print *,'error reading x,year,month from '
	print *,line2
	call abort
  999	continue
	end
	function llen(a)
	character*(*) a
	do 10 i=len(a),1,-1
	    if(a(i:i).ne.'?' .and. a(i:i).ne.' ')goto 20
   10	continue
	llen=len(a)
   20	continue
	llen = i
	end

