program diffdat
!
!   compute first derivative of a time series
!
    implicit none
    include 'param.inc'
    integer :: nperyear,n
    real :: data(npermax,yrbeg:yrend),diff(npermax,yrbeg:yrend),minfac
    character :: file*256,var*80,units*40,string*10,lvar*100,svar*100,history*10000, &
        metadata(2,100)*2000
    logical :: lwrite

    minfac = 0.7            ! arbitrary
    call get_command_argument(1,file)
    if ( file == ' ' ) then
        write(0,*) 'usage: diffdat file [n [debug]]'
        write(0,*) '       computes first derivative of file using n consecutive values'
        call exit(-1)
    end if
    call get_command_argument(2,string)
    lwrite = .false. 
    if ( string /= ' ' ) then
        read(string,*) n
        if ( n < 2 ) n = 2
        if ( n >= 3 .and. mod(n,2) /= 1 ) n = n + 1
        call get_command_argument(3,string)
        if ( string == 'lwrite' .or. string == 'debug' ) then
            lwrite = .true. 
        end if
    else
        n = 3
    end if
    call readseriesmeta(file,data,npermax,yrbeg,yrend,nperyear,var,units,lvar,svar, &
        history,metadata,.false.,lwrite)
    call derivative(n,data,diff,npermax,yrbeg,yrend,nperyear,minfac,lwrite)
    if ( n == 2 ) then
        print '(3a)','# difference with previous data point of ',trim(file)
    else
        print '(2a)','# centered derivative of ',trim(file)
        if ( n > 3 ) print '(a,i3,a)','# using linear regression over ',n,' data points'
    end if
    var = 'd_'//trim(var)//'_dt'
    lvar = 'time derivative of '//lvar
    call nperyear2units(nperyear,string)
    units = trim(units)//'/'//string
    call printvar(6,var,units,lvar)
    call copyheadermeta(file,6,' ',history,metadata)
    call printdatfile(6,diff,npermax,nperyear,yrbeg,yrend)
end program diffdat

                
