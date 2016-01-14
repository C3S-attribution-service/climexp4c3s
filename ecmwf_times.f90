program ecmwf_times
!
!   compute the parameters needed for the new ECMWF web site
!
    implicit none
    integer year,mo,lead
    integer i,j,y,n,yrend,moend,dpm(12,2)
    integer leap,iargc
    character string*100
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
   
   write(*,'(i4,3i2.2,a,i4,a,i4,3i2.2)') year,mo,1,0,',',24*(n-1),',', &
&     yrend,moend,dpm(moend,leap(yrend)),0

end program