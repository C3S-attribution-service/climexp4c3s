        program plotsoi
*
*       convert the SOI data from http://nic.fb4.noaa.gov:80/data/cddb/
*	or the KNMI "Phil Jones" series from Marleen Kaltofen.
*       into a format suitable for gnuplot
*
        implicit none
        integer yrbeg,yrend
        parameter (yrbeg=1800,yrend=2020)
        integer i,j,k,l,year,loop,type,mean
        real soi(12,yrbeg:yrend),soimean(12,yrbeg:yrend),
     +        soisd(12,yrbeg:yrend),s
        character line*80
	integer iargc
*
*       init
*
	type = 3
	if ( iargc().ge.1 ) then
	    call getarg(1,line)
	    if ( line(1:4).eq.'knmi' ) then
		type = 1
            elseif ( line(1:4).eq.'jone' ) then
                type = 3
	    elseif ( line(1:4).eq.'noaa' ) then
		type = 2
	    else
		print '(2a)','Please specify knmi, jones or noaa, not '
     +                ,line
		stop
	    endif
	endif
        if ( iargc().ge.2 ) then
            call getarg(2,line)
            read(line,*) mean
        else
            mean = 12
        endif
        call makeabsent(soi,12,yrbeg,yrend)
*
	if ( type.eq.1 ) then
*
*           read KNMI data
*
            open(1,file='soi.knmi',status='old')
            do i=1,5
                read(1,'(a)') line
                print '(2a)','# ',line
            enddo
            do i=1866,1997
                read(1,*) year,(soi(j,i),j=1,12)
                if ( year.ne.i ) print *,'error: year != i: ',year,i
            enddo
            close(1)
	elseif ( type.eq.3 ) then
*
*       read Phil Jones data
*
            open(1,file='soi.dat',status='old')
            do i=1,5
                read(1,'(a)') line
                print '(2a)','# ',line
            enddo
            do i=1866,1997
                read(1,*) year,(soi(j,i),j=1,12)
                if ( year.ne.i ) print *,'error: year != i: ',year,i
            enddo
            close(1)
	elseif ( type.eq.2 ) then
*
*           read historical data
*
            open(1,file='soi.his',status='old')
            do i=1,3
                read(1,'(a)') line
                print '(2a)','# ',line
            enddo
            do i=1882,1950
                read(1,'(i4,x,12f6.1)') year,(soi(j,i),j=1,12)
                if ( year.ne.i ) print *,'error: year != i: ',year,i
            enddo
            close(1)
*           
*           read modern data
*           
            open(1,file='soi',status='old')
            do loop = 1,2
                do i=1,4
                    read(1,'(a)') line
                    if ( loop.eq.2 ) print '(2a)','# ',line
                enddo
                do i=1951,yrend
                    read(1,'(i4,x,12f6.1)') year,(soi(j,i),j=1,12)
                    if ( year.ne.i ) print *,'error: year != i: ',year,i
                enddo
            enddo
            close(1)
	else
	    print *,'error: unknown type ',type
	    stop
	endif
*
*       get n-month running mean
*
        call rmean(soi(1,yrbeg),soimean(1,yrbeg),soisd(1,yrbeg),yrbeg
     +        ,yrend,mean)
*
*       print
*
        do i=yrbeg,yrend
            do j=1,12
                if ( abs(soi(j,i)).lt.1e33 ) then
                    if ( abs(soimean(j,i)).lt.1e33 ) then
                        print '(i4,i3,2f8.2)',i,j,soi(j,i),soimean(j,i)
                    else
                        print '(i4,i3,2f8.2)',i,j,soi(j,i)
                    endif
                else
                    print '(a)'
                endif
            enddo
        enddo
*        
        end
