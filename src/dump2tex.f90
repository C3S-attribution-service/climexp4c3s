program dump2tex

!   hak

    implicit none
    integer :: yr,mo
    real :: indx,var,scale
    character string*80

    scale = 1
    if ( command_argument_count() == 1 ) then
        call get_command_argument(1,string)
        read(string,*) scale
    endif

    100 continue
    read (*,'(a)',end=200,err=900) string
    if ( string(1:1) == '#' .or. string == ' ' ) goto 100
    read(string,*) indx,var,yr,mo
    var = var*scale
!**        if ( yr.gt.1800 ) yr = yr-1900
    1000 format(a,f8.4,a,f8.4,a,i4.2,a)
    if ( yr > 2000 ) then
        print *,'donot know how to format year ',yr
        stop
    elseif ( yr > 1900 ) then
        print 1000,'\put(',indx,',',var &
        ,'){\makebox(0,0){$',yr-1900,'$}}'
    elseif ( yr > 1800 ) then
        print 1000,'\put(',indx,',',var &
        ,'){\makebox(0,0){$\underline{',yr-1800,'}$}}'
    elseif ( yr > 1700 ) then
        print 1000,'\put(',indx,',',var &
        ,'){\makebox(0,0){$\underline{\underline{',yr-1700 &
        ,'}}$}}'
    else
        print *,'donot know how to format year ',yr
        stop
    endif
    goto 100
200 continue
    stop
    900 print *,'error reading data'
end program dump2tex
