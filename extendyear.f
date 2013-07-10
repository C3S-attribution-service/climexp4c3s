        program extendyear
*
*       proglet to repeat a year with months number 12-24
*
        implicit none
        integer month,i,j,oldmonth
        character line*200,c*160,rest(10,0:12)*160
        integer llen
        external llen
*
        oldmonth = -1
        do month=0,12
            do i=1,10
                rest(i,month) = ' '
            enddo
        enddo
  100   continue
        read(*,'(a)',end=800,err=800) line
        if ( line.eq.' ' ) then
            goto 100
        endif
        read(line,'(i3,a)') month,c
        if ( month.ne.oldmonth ) then
            i = 0
            oldmonth = month
        endif
        if ( month.ge.0 .and. month.le.12 ) then
            i = i+1
            if ( i.gt.10 ) then
                print *,'error: more than 10 copies of month ',month
                call abort
            endif
            rest(i,month) = c
        endif
        goto 100
*
  800   continue
        do j=1,i
	    if ( llen(rest(j,0)).gt.1 ) then
            	print '(i3,a)',-1,rest(j,0)(1:llen(rest(j,0)))
            	print '(i3,a)', 0,rest(j,0)(1:llen(rest(j,0)))
            endif
            do month=1,12
		if ( llen(rest(j,month)).gt.1 ) then
                    print '(i3,a)',
     +                month,rest(j,month)(1:llen(rest(j,month)))
		else
		    print '(a)'
		endif
            enddo
            do month=1,12
		if ( llen(rest(j,month)).gt.1 ) then
                    print '(i3,a)',12+
     +                month,rest(j,month)(1:llen(rest(j,month)))
		else
		    print '(a)'
		endif
            enddo
            print '(a)'
            print '(a)'
        enddo
        end
