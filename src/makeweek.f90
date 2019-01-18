program makeweek

!	make a time series with 1 for Mon-Fri, 0 for Saturday and Sunday
!	and one with 0 for Sunday only

    implicit none
    integer,parameter :: yrbeg=1600,yrend=2100
    integer :: yr,mo,dy,jul,jul0,jul1,ikind,weekday
    real :: week
    logical :: lwrite
    integer,external :: julday,leap

    do ikind=1,2
        if ( ikind == 1 ) then
            open(1,file='week.dat')
            write(1,'(a)') '# Mon-Fri = 1, Sat-Sun= = 0'
            write(1,'(a)') '# weekday [1] working days'
        else
            open(1,file='oldweek.dat')
            write(1,'(a)') '# Mon-Sat = 1, SSun= = 0'
            write(1,'(a)') '# weekday [1] old working days'
        end if
        jul0 = julday( 1, 1,yrbeg)
        jul1 = julday(12,31,yrend)
        do jul=jul0,jul1
            call caldat(jul,mo,dy,yr)
            weekday = 1 + mod(jul,7)
            if ( weekday <= 5 .or. &
            weekday == 6 .and. ikind == 2 ) then
                week = 1
            else
                week = 0
            end if
            write(1,'(i4,i2.2,i2.2,f4.0)') yr,mo,dy,week
        end do
    end do
end program makeweek