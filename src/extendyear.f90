program extendyear

!   proglet to repeat a year with months number 12-24

    implicit none
    integer :: month,i,j,oldmonth
    character :: line*200,c*160,rest(10,0:12)*160

    oldmonth = -1
    do month=0,12
        do i=1,10
            rest(i,month) = ' '
        end do
    end do
100 continue
    read(*,'(a)',end=800,err=800) line
    if ( line == ' ' ) then
        goto 100
    end if
    read(line,'(i3,a)') month,c
    if ( month /= oldmonth ) then
        i = 0
        oldmonth = month
    end if
    if ( month >= 0 .and. month <= 12 ) then
        i = i+1
        if ( i > 10 ) then
            print *,'error: more than 10 copies of month ',month
            call exit(-1)
        end if
        rest(i,month) = c
    end if
    goto 100

800 continue
    do j=1,i
        if ( len_trim(rest(j,0)) > 1 ) then
            print '(i3,a)',-1,trim(rest(j,0))
            print '(i3,a)', 0,trim(rest(j,0))
        end if
        do month=1,12
            if ( len_trim(rest(j,month)) > 1 ) then
                print '(i3,a)', &
                month,trim(rest(j,month))
            else
                print '(a)'
            end if
        end do
        do month=1,12
            if ( len_trim(rest(j,month)) > 1 ) then
                print '(i3,a)',12+month,trim(rest(j,month))
            else
                print '(a)'
            end if
        end do
        print '(a)'
        print '(a)'
    end do
end program extendyear
