    program makesnow

!       crude estimate of the amount of snow

    implicit none
    integer :: yrbeg,yrend,npermax
    parameter(yrbeg=1700,yrend=2020,npermax=366)
    integer :: m,yr,type,n,nperyear
    real :: temp(npermax,yrbeg:yrend),prcp(npermax,yrbeg:yrend), &
    snow(npermax,yrbeg:yrend),sd
    character string*256,var*80,units*80
    logical :: lwrite
    real :: erfcc
    external erfcc

    lwrite = .false. 
    type = 2
    sd = 5
    if ( command_argument_count() /= 2 ) then
        print *,'usage: makesnow temp_file prcp_file'
        call exit(-1)
    endif

    call get_command_argument(1,string)
    call readseries(string,temp,npermax,yrbeg,yrend,nperyear, &
    var,units, .true. ,lwrite)
    print '(2a)','# temperature file ',trim(string)
    call get_command_argument(2,string)
    call readseries(string,prcp,npermax,yrbeg,yrend,n, &
    var,units, .true. ,lwrite)
    if ( n /= nperyear ) then
        write(0,*) 'error: unequal time scales ',n,nperyear
        call exit(-1)
    end if
    print '(2a)','# precipitation file ',trim(string)
    if ( type == 1 .or. nperyear > 12 ) then
        print '(3a)','# snow [',trim(units), &
        '] ''Snow'' computed as P when T<=0'
    elseif ( type == 2 ) then
        print '(3a,f5.0,a)','# snow [',units, &
        '] ''Snow'' computed as P when T<0 with ' &
        ,sd,'o Gauss'
    elseif ( type == 3 ) then
        print '(a)','# ''Snow'' computed as T*P when T<0'
    else
        print '(a)','wrong type of snow'
        stop
    endif

    do yr=yrbeg,yrend
        do m=1,nperyear
            if ( temp(m,yr) < 1e33 .and. prcp(m,yr) < 1e33 ) then
                if ( type == 1 .or. nperyear > 12 ) then
                !                   zeroth approximation
                    if ( temp(m,yr) <= 0 ) then
                        snow(m,yr) = prcp(m,yr)
                    else
                        snow(m,yr) = 0
                    endif
                elseif ( type == 2 ) then
                !                       same but smeared with a Gaussian with s.d. sd
                    snow(m,yr) = prcp(m,yr)*erfcc(temp(m,yr)/sd)/2
                elseif ( type == 3 ) then
                !                       give more weight to T: T*P when T<0.
                    if ( temp(m,yr) < 0 ) then
                        snow(m,yr) = prcp(m,yr)*temp(m,yr)
                    else
                        snow(m,yr) = 0
                    endif
                else
                    print '(a)','wrong type of snow'
                    stop
                endif
            else
                snow(m,yr) = 3e33
            endif
        enddo
    enddo
    call printdatfile(6,snow,npermax,nperyear,yrbeg,yrend)
end program makesnow
