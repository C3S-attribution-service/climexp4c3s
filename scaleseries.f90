program scaletimeseries

!   scale a time series

    implicit none
    include 'param.inc'
    integer :: i,j,nperyear,ndpm
    real :: data(npermax,yrbeg:yrend),factor,offset
    character string*1024,var*40,units*20
    logical :: lwrite
    integer :: iargc

    lwrite = .false. 
    if ( iargc() < 2 ) then
        print *,'usage: scaleseries factor[:offset] file [ndpm]'
        stop
    endif
    call getarg(1,string)
    i = index(string,':')
    if ( i == 0 ) then
        read(string,*) factor
        offset = 0
    else
        read(string(:i-1),*) factor
        read(string(i+1:),*) offset
    end if
    if ( iargc() > 2 ) then
        call getarg(3,string)
        read(string,*,err=10) ndpm
        goto 20
        10 continue
        ndpm = 0
        20 continue
    else
        ndpm = 0
    endif
    call getarg(2,string)
    call readseries(string,data,npermax,yrbeg,yrend,nperyear,var &
    ,units, .false. ,lwrite)
    if ( ndpm /= 0 .and. nperyear /= 12 ) then
        write(0,*) 'scaleseries: error: can only scale by dpm if nperyear = 12, not ',nperyear
        write(*,*) 'scaleseries: error: can only scale by dpm if nperyear = 12, not ',nperyear
        call abort
    endif
    call copyheader(string,6)
    write(6,'(a,g16.6)') '# scaled with a factor ',factor
    if ( ndpm == 1 ) then
        write(6,'(a)')'# and multiplied by the number of days in a month'
    elseif ( ndpm == 2 ) then
        write(6,'(a)') '# and multiplied by the number of days in a month squared'
    elseif ( ndpm == -1 ) then
        write(6,'(a)')'# and divided by the number of days in a month'
    elseif ( ndpm == -2 ) then
        write(6,'(a)') '# and divided by the number of days in a month squared'
    elseif ( ndpm /= 0 ) then
        write(6,'(a,i3)') '# and multiplied by dpm**',ndpm
    endif
    if ( offset /= 0 ) then
        write(6,'(a,g16.6)') '# added offset ',offset
    end if
    call scaleseries(data,npermax,nperyear,yrbeg,yrend,factor,offset,ndpm)
    call printdatfile(6,data,npermax,nperyear,yrbeg,yrend)
end program scaletimeseries
