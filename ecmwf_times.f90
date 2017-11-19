program ecmwf_times
!
!   compute the parameters needed for the new ECMWF web site
!
    implicit none
    integer year,mo,lead
    integer i,j,y,m,n,yrend,moend,dyend,dpm(12,2)
    integer leap,iargc
    character string*100,version*10,nformat*2
    data dpm /31,28,31,30,31,30,31,31,30,31,30,31, &
    &         31,29,31,30,31,30,31,31,30,31,30,31/
    
    if ( iargc().lt.3 ) then
        write(0,*) 'usage: ecmwf_time year mo lead'
        write(0,*) 'gives: analysis time, hours to end of averaging period, date of end of average period'
        call exit(-1)
    end if
    call getarg(1,string)
    read(string,*) year
    call getarg(2,string)
    read(string,*) mo
    call getarg(3,string)
    read(string,*) lead
    call getarg(4,version)
    
    if ( version == ' ' ) then
        yrend = year
        moend = mo + lead + 3 - 1
        n = 0
        do i=mo,moend
            j = i
            y = year
            if ( j > 12  ) then
                j = j - 12
                y = y + 1
            end if
            n = n + dpm(j,leap(y))
        end do
        if ( moend > 12 ) then
            moend = moend - 12
            yrend = yrend + 1
        end if
        n = n - 1
        dyend = dpm(moend,leap(yrend))
    else if ( version == 'seas5alt' ) then
        yrend = year
        moend = mo + lead
        n = 0
        do i=mo,moend-1
            j = i
            y = year
            if ( j > 12  ) then
                j = j - 12
                y = y + 1
            end if
            n = n + dpm(j,leap(y))
        end do
        if ( moend > 12 ) then
            moend = moend - 12
            yrend = yrend + 1
        end if
        dyend = 1
    else if ( version == 'seas5' ) then
        yrend = year
        moend = mo + lead
        if ( lead == 1 ) then
            n = 31 
        else if ( lead == 2 ) then
            n = 62
        else if ( lead == 3 ) then
            n = 92
        else if ( lead == 4 ) then
            n = 123
        else
            write(0,*) 'ecmwf_times: error: cannot handle lead >1 yet:',lead
            call exit(-1)
        end if
        m = 0
        do i=mo,moend-1
            j = i
            y = year
            if ( j > 12  ) then
                j = j - 12
                y = y + 1
            end if
            m = m + dpm(j,leap(y))
        end do
        if ( moend > 12 ) then
            moend = moend - 12
            yrend = yrend + 1
        end if
        dyend = n - m + 1
    else
        write(0,*) 'ecmwf_time: unknown version ',trim(version)
        call exit(-1)
    end if
    n = 24*n
    if ( n < 1000 ) then
        nformat = 'i3'
    else if ( n < 10000 ) then
        nformat = 'i4'
    else
        nformat = 'i5'
    end if
    write(*,'(i4,3i2.2,a,'//nformat//',a,i4,3i2.2)') year,mo,1,0,',',n,',', &
&       yrend,moend,dyend,0

end program