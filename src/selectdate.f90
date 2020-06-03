program selectdate

!   scale a time series and optionally apply an offset

    implicit none
    include 'param.inc'
    integer :: i,j,nperyear,yr,month,day
    real :: data(npermax,yrbeg:yrend),newdata(yrbeg:yrend)
    character :: string*1024,var*40,units*40,lvar*120,svar*120,metadata(2,100)*2000, &
        history*50000,title*500
    logical :: lwrite
    integer,save :: dpm(12)
    character(3),save :: months(12),seasons(4)
    character(7),save :: halfyears(2)
    data months &
    /'Jan','Feb','Mar','Apr','May','Jun' &
    ,'Jul','Aug','Sep','Oct','Nov','Dec'/
    data seasons /'DJF','MAM','JJA','SON'/
    data halfyears /'Oct-Mar','Apr-Sep'/
    data dpm /31,29,31,30,31,30,31,31,30,31,30,31/
    logical,external :: isnumchar

    lwrite = .false. 
    if ( command_argument_count() < 2 ) then
        print *,'usage: selectdate file month [date]'
        call exit(-1)
    endif
    call get_command_argument(1,string)
    call readseriesmeta(string,data,npermax,yrbeg,yrend,nperyear,var,units,lvar,svar,history,metadata, &
        .false.,lwrite)
    call get_command_argument(2,string)
    read(string,*) month
    day = 1
    if ( command_argument_count() > 2 ) then
        call get_command_argument(3,string)
        if ( isnumchar(string(1:1)) ) then
            read(string,*) day
        end if
    end if
    title = ' '
    if ( month <= 0 .or. month > min(nperyear,12) ) then
        write(0,*) 'selectdate: invalid value for month ',month, &
            ', should be less than or equal to ',min(nperyear,12)
        call exit(-1)
    end if
    if ( nperyear > 1 ) then
        j = month
        if ( nperyear == 2 ) then
            lvar = lvar//' in '//halfyears(month)
        else if ( nperyear == 4 ) then
            lvar = lvar//' in '//seasons(month)
        else if ( nperyear == 12 ) then
            lvar = lvar//' in '//months(month)
        else
            write(lvar,'(2a,i2,2a)') trim(lvar),' on ',day,' ',months(month)
            j = day
            do i=1,month-1
                j = j + dpm(month)
            end do
        end if
    end if
    call printvar(6,var,units,lvar)
    call get_command_argument(1,string)
    call printmetadata(6,string,' ',title,history,metadata)
    do yr=yrbeg,yrend
        newdata(yr) = data(j,yr)
    end do
    call printdatfile(6,newdata,1,1,yrbeg,yrend)
end program selectdate

