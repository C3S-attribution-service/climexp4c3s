program diamond2year

!   replace the strings '\Diamond' by the year

    implicit none
    integer :: i,year,month,oldmonth,l,init,n,m
    real :: x(10),xmin,xmax,ymin,ymax
    character line*200,line2*180,file*1000
    logical :: lmonth
    integer,external :: getnumwords

    if ( command_argument_count() < 2 ) then
        print *,'usage: diamond2year psfile datfile [xmin xmax ymin ymax] > outfile'
        call exit(-1)
    endif
    call get_command_argument(1,file)
    open(1,file=file,status='old')
    call get_command_argument(2,file)
    open(2,file=file,status='old')
    if ( command_argument_count() > 2 ) then
        call get_command_argument(3,line)
        read(line,*) xmin
        call get_command_argument(4,line)
        read(line,*) xmax
        call get_command_argument(5,line)
        read(line,*) ymin
        call get_command_argument(6,line)
        read(line,*) ymax
    else
        xmin = -3e33
        xmax = +3e33
        ymin = -3e33
        ymax = +3e33
    end if

    init = 0
    oldmonth = 0
    lmonth = .false. 
    100 continue
    read(1,'(a)',end=800,err=900) line
    i = index(line,' Pls    ')
    l = len_trim(line)
    if ( i /= 0 ) then
        if ( init == 0 ) then
            init = 1
            write(*,'(a)') 'LTa'
            write(*,'(a)') '(Helvetica) findfont 100 scalefont setfont'
        endif
    110 continue
        read(2,'(a)') line2
        if ( line2(1:1) == '#' .or. line2 == ' ' ) goto 110
        n = getnumwords(line2)
        read(line2,*,err=901) (x(m),m=1,2),year,month
        if ( n == 6 ) then
            if ( x(1) < xmin .or. x(1) > xmax ) goto 110
            if ( x(2) < ymin .or. x(2) > ymax ) goto 110
        end if
        if ( lmonth .or. oldmonth == 0 ) then
            if ( month < 10 ) then
                write(*,'(2a,i1,a,i4,a)') line(1:i-1), &
                ' M (',month,'.',year,') Cshow'
            elseif ( month < 100 ) then
                write(*,'(2a,i2,a,i4,a)') line(1:i-1), &
                ' M (',month,'.',year,') Cshow'
            else
                write(*,'(2a,i3,a,i4,a)') line(1:i-1), &
                ' M (',month,'.',year,') Cshow'
            endif
        else
            write(*,'(2a,i4,a)') line(1:i-1),' M (',year,') Cshow'
        endif
        if ( month /= oldmonth ) then
            if ( oldmonth > 0 ) lmonth = .true. 
            oldmonth = month
        endif
    else
        if ( init == 1 ) then
            init = 0
            write(*,'(a)') 'LTa'
            write(*,'(a)') '(Helvetica) findfont 140 scalefont setfont'
        endif
        write(*,'(a)') line(:l)
    endif
    goto 100

800 continue
    goto 999
900 print *,'error reading input'
    print *,line
    call exit(-1)
901 print *,'error reading x,year,month from '
    print *,line2
    call exit(-1)
999 continue
end program diamond2year
