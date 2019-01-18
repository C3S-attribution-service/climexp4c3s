subroutine getseriesmoment(moment,data,npermax,yrbeg,yrend, &
    nperyear,yr1,yr2,result)

!   compute the moment of the time series in data over yr1:yr2
!   I do not use the fancy ep-trick of numerical recipes.

    implicit none
    integer :: moment,npermax,yrbeg,yrend,nperyear,yr1,yr2
    real :: data(npermax,yrbeg:yrend),result
    integer :: i,j,n
    real :: s1,s2,s3,s4
    logical,parameter :: lwrite=.false.

    if ( lwrite ) then
        print *,'getseriesmoment: input: moment = ',moment
        print *,'                        yr1,yr2= ',yr1,yr2
    endif

    s1 = 3e33
    s2 = 3e33
    s3 = 3e33
    s4 = 3e33
    if ( moment >= 1 ) then
        s1 = 0
        n = 0
        do i=yr1,yr2
            do j=1,nperyear
                if ( data(j,i) < 1e30 ) then
                    n = n + 1
                    s1 = s1 + data(j,i)
                endif
            enddo
        enddo
        if ( n > 0 ) then
            s1 = s1/n
        else
            goto 800
        endif
        if ( lwrite ) print *,'s1 = ',s1
    endif
    if ( moment >= 2 ) then
        s2 = 0
        do i=yr1,yr2
            do j=1,nperyear
                if ( data(j,i) < 1e30 ) then
                    s2 = s2 + (data(j,i)-s1)**2
                endif
            enddo
        enddo
        if ( n > 1 ) then
            s2 = s2/(n-1)
        else
            goto 800
        endif
        if ( s2 < 0 ) then
            write(0,*) 'getmoment: error: s2<0 ',s2
            goto 800
        else
            s2 = sqrt(s2)
        endif
        if ( lwrite ) print *,'s2 = ',s2
    endif
    if ( moment >= 3 ) then
        s3 = 0
        do i=yr1,yr2
            do j=1,nperyear
                if ( data(j,i) < 1e30 ) then
                    s3 = s3 + (data(j,i)-s1)**3
                endif
            enddo
        enddo
        if ( s2**3 == 0 ) then
        !**             write(0,*) 'getmoment: error: s2 ~ 0',s2
            goto 800
        else
            s3 = s3/(n*s2**3)
        endif
        if ( lwrite ) print *,'s3 = ',s3
    endif
    if ( moment >= 4 ) then
        s4 = 0
        do i=yr1,yr2
            do j=1,nperyear
                if ( data(j,i) < 1e30 ) then
                    s4 = s4 + (data(j,i)-s1)**4
                endif
            enddo
        enddo
        if ( s2**4 == 0 ) then
            goto 800
        else
            s4 = s4/(n*s2**4) - 3
        endif
        if ( lwrite ) print *,'s4 = ',s4
    endif

800 continue
    if ( moment == 1 ) then
        result = s1
    elseif ( moment == 2 ) then
        result = s2
    elseif ( moment == 3 ) then
        result = s3
    elseif ( moment == 4 ) then
        result = s4
    else
        result = 3e33
    endif

end subroutine getseriesmoment
