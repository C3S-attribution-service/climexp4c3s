program operate

!    Apply a logarithm, sqrt, ... to a timeseries in standard format and return the resulting series.

    implicit none
    include 'param.inc'
    integer :: i,j,yr,mo,nyr,nperyear
    real :: data(npermax,yrbeg:yrend)
    character :: file*1024,line*128,oper*3,var*100,units*100,lvar*200,svar*100, &
        history*10000,metadata(2,100)*2000,title*100
    logical :: lwrite,lstandardunits

    lwrite = .false. 
    lstandardunits = .false. 
    if ( command_argument_count() < 2 ) then
        write(0,*) &
        'usage: operate log|sqrt|exp file'
        stop
    endif
    call get_command_argument(1,oper)
    call get_command_argument(2,file)
    call readseriesmeta(file,data,npermax,yrbeg,yrend,nperyear,var,units, &
        lvar,svar,history,metadata,lstandardunits,lwrite)

    if ( oper == 'log' ) then
        call takelog(data,npermax,nperyear,yrbeg,yrend)
        if ( units /= '1' ) units = 'log('//trim(units)//')'
        var = 'log('//trim(var)//')'
        lvar = 'logarithm of '//lvar
    elseif ( oper == 'sqrt' ) then
        call takesqrt(data,npermax,nperyear,yrbeg,yrend)
        if ( units /= '1' ) units = 'sqrt('//trim(units)//')'
        var = 'sqrt('//trim(var)//')'
        lvar = 'sqrt of '//lvar
    elseif ( oper == 'exp' ) then
        call takeexp(data,npermax,nperyear,yrbeg,yrend)
        if ( units /= '1' ) units = 'exp('//trim(units)//')'
        var = 'exp('//trim(var)//')'
        lvar = 'exp of '//lvar
    elseif ( oper == 'inv' ) then
        call takeinv(data,npermax,nperyear,yrbeg,yrend)
        if ( units /= '1' ) units = '1/('//trim(units)//')'
        var = '1/('//trim(var)//')'
        lvar = 'inverse of '//lvar
    else
        write(0,'(2a)') 'operate: unknown operation ',oper
        call abort
    endif

    call printvar(6,var,units,lvar)
    title = trim(oper)//' of '//trim(file)
    call copyheadermeta(file,6,title,history,metadata)
    call printdatfile(6,data,npermax,nperyear,yrbeg,yrend)

end program operate