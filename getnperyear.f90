 program getnperyear

!       get nperyear

    implicit none
    include 'param.inc'
    integer :: nperyear
    real :: data(npermax,yrbeg:yrend)
    character line*255,var*40,units*20
    logical :: lwrite

!       init

    lwrite = .false. 
    if ( command_argument_count() == 0 ) then
        print *,'usage: getnperyear datafile'
        call exit(-1)
    endif

    call get_command_argument(1,line)
    call readseries(line,data(1,yrbeg),npermax,yrbeg,yrend,nperyear &
    ,var,units, .false. ,lwrite)
    if ( nperyear < 10 ) then
        print '(i1)',nperyear
    elseif ( nperyear < 100 ) then
        print '(i2)',nperyear
    elseif ( nperyear < 1000 ) then
        print '(i3)',nperyear
    elseif ( nperyear < 10000 ) then
        print '(i4)',nperyear
    endif
end program getnperyear
