program diffdat
!
!   compute first derivative of a time series
!
    implicit none
    include 'param.inc'
    integer :: nperyear,n
    real :: data(npermax,yrbeg:yrend),diff(npermax,yrbeg:yrend),minfac
    character file*256,var*80,units*40,string*10
    logical :: lwrite

    minfac = 0.7            ! arbitrary
    call getarg(1,file)
    if ( file == ' ' ) then
        write(0,*) 'usage: diffdat file [n [debug]]'
        write(0,*) '       computes first derivative of file using n consecutive values'
        call exit(-1)
    end if
    call getarg(2,string)
    lwrite = .false. 
    if ( string /= ' ' ) then
        read(string,*) n
        if ( n < 2 ) n = 2
        if ( n >= 3 .and. mod(n,2) /= 1 ) n = n + 1
        call getarg(3,string)
        if ( string == 'lwrite' .or. string == 'debug' ) then
            lwrite = .true. 
        end if
    else
        n = 3
    end if
    call readseries(file,data,npermax,yrbeg,yrend,nperyear,var,units, .false. ,lwrite)
    if ( n == 2 ) then
        call ndiffit(data,npermax,nperyear,yrbeg,yrend,1,minfac)
        diff = data
    else
        call derivative(n,data,diff,npermax,yrbeg,yrend,nperyear,minfac,lwrite)
    end if
    if ( n == 2 ) then
        print '(3a)','# difference with previous year of ',trim(file)
    else
        print '(3a,i3,a)','# centered derivative of ',trim(file) &
            ,' using linear regression over ',n,' data points'
    end if
    call copyheader(file,6)
    call printdatfile(6,diff,npermax,nperyear,yrbeg,yrend)
end program diffdat

                
