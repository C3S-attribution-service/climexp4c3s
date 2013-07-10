        program plotsoi
*
*       convert the NAO data from Wilco and Phil Jones.
*       into a format suitable for gnuplot
*
        implicit none
        integer yrbeg,yrend
        parameter (yrbeg=1800,yrend=2020)
        integer i,j,k,l,year,loop,type,mean
        real nao(12,yrbeg:yrend,2),naomean(yrbeg:yrend,2),s(2)
        character line*80
	integer iargc
*
*       init
*
	type = 1
        mean = 2                ! 5-yr running mean
	if ( iargc().ge.1 ) then
	    call getarg(1,line)
	    if ( line(1:4).eq.'DJF ' ) then
		type = 3
            elseif ( line(1:4).eq.'DJFM' ) then
                type = 4
	    else
		print '(2a)','Please specify DJFor DJFM, not '
     +                ,line
		stop
	    endif
	endif
        call makeabsent(nao(1,yrbeg,1),12,yrbeg,yrend)
        call makeabsent(nao(1,yrbeg,2),12,yrbeg,yrend)
*
        call readdat(nao(1,yrbeg,1),yrbeg,yrend,'../CRUData/nao.dat')
        call readdat(nao(1,yrbeg,2),yrbeg,yrend
     +        ,'../CRUData/nao_ijs_azo.dat')
	if ( type.eq.1 ) then
*
*           print
*
            do i=yrbeg,yrend
                do j=1,12
                    if ( abs(nao(j,i,2)).lt.1e33 ) then
                        print '(i4,i3,2f8.2)',i,j,nao(j,i,1),nao(j,i,2)
                    elseif ( abs(nao(j,i,1)).lt.1e33 ) then
                        print '(i4,i3,2f8.2)',i,j,nao(j,i,1)
                    else
                        print '(a)'
                    endif
                enddo
            enddo
        else
            do i=yrbeg+1,yrend
                do k=1,2
                    naomean(i,k) = nao(12,i-1,k)
                    do j=1,type-1
                        naomean(i,k) = naomean(i,k) + nao(j,i,k)
                    enddo
                    if ( naomean(i,k).lt.1e33 ) then
                        naomean(i,k) = naomean(i,k)/type
                    else
                        naomean(i,k) = 3e33
                    endif
                enddo
            enddo
            do i=yrbeg,yrend
                do k=1,2
                    s(k) = 0
                    do j=i-mean,i+mean
                        if ( j.ge.yrbeg .and. j.le.yrend ) then
                            s(k) = s(k) + naomean(j,k)
                        else
                            s(k) = 3e33
                        endif
                    enddo
                    if ( s(k).lt.1e33 ) then
                        s(k) = s(k)/(2*mean+1)
                    else
                        s(k) = 3e33
                    endif
                enddo
 1000           format(i5,4f8.2)
                if ( abs(s(2)).lt.1e33 ) then
                    print 1000,i+1,naomean(i,1),s(1),naomean(i,2),s(2)
                elseif ( abs(s(1)).gt.1e33 .and. abs(naomean(i,2)).lt
     +                    .1e33 )then
                    print 1000,i+1,naomean(i,1),0.,naomean(i,2)
                elseif ( abs(s(1)).lt.1e33 .and. abs(naomean(i,2)).lt
     +                    .1e33 )then
                    print 1000,i+1,naomean(i,1),s(1),naomean(i,2)
                elseif ( abs(s(1)).lt.1e33 ) then
                    print 1000,i+1,naomean(i,1),s(1)
                elseif ( abs(naomean(i,1)).lt.1e33 ) then
                    print 1000,i+1,naomean(i,1)
                else
                    print '(a)'
                endif
            enddo
        endif
*
        end
