program untransform

!   remove the transformation in a gumbel or log plot

    implicit none
    integer :: i,n,itrans,yyyymm,nwords
    real :: xx(4)
    character :: line*256,flag*10
    integer,external :: getnumwords

    if ( command_argument_count() /= 3 ) then
        write(0,*) &
        'usage: untransform true|false true|false true|false'
        write(0,*) '       transforms a log|sqrt|square back'
        stop
    endif

    itrans = 0
    do i=1,3
        call get_command_argument(i,flag)
        if ( flag /= 'false' .and. flag /= 'off' ) then
            if ( itrans /= 0 ) then
                write(0,*) 'untransform: error: can only handle'// &
                ' one transformation at a time'
                call abort
            endif
            itrans = i
        endif
    enddo
100 continue
    read(*,'(a)',end=800) line
    if ( line(1:1) == '<' ) then
        print '(2a)','# ',trim(line)
        goto 100
    endif
    if ( line(1:1) == '#' .or. line == ' ' ) then
        print '(a)',trim(line)
        goto 100
    endif
    if ( index(line,'amoeba') /= 0 ) goto 100
    nwords = getnumwords(line)
    if ( nwords == 5 ) then
        read(line,*) n,xx
    elseif ( nwords == 6 ) then
        read(line,*) n,xx,yyyymm
    else
        write(0,'(2a)') 'untransform: funny line: ',trim(line)
        write(*,'(2a)') '# untransform: funny line: ',trim(line)
    endif
    if ( itrans == 1 ) then
        if ( xx(2) /= -999.9 ) xx(2) = 10.**xx(2)
        if ( xx(3) /= -999.9 ) xx(3) = 10.**xx(3)
    elseif ( itrans == 2 ) then
        if ( xx(2) /= -999.9 ) xx(2) = xx(2)**2
        if ( xx(3) /= -999.9 ) xx(3) = xx(3)**2
    elseif ( itrans == 3 ) then
        if ( xx(2) /= -999.9 ) then
            if ( xx(2) >= 0 ) then
                xx(2) = sqrt(xx(2))
            else
                write(0,*) 'untransform: error: negative number: ' &
                ,xx(2)
            endif
        endif
        if ( xx(3) /= -999.9 ) then
            if ( xx(3) >= 0 ) then
                xx(3) = sqrt(xx(3))
            else
                write(0,*) 'untransform: error: negative number: ' &
                ,xx(3)
            endif
        endif
    endif
    if ( nwords == 5 ) then
        write(*,'(i8,4g22.6,i10)') i,xx
    elseif ( nwords == 6 ) then
        write(*,'(i8,4g22.6,i10)') i,xx,yyyymm
    endif
    goto 100
800 continue
end program untransform
