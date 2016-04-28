program timeshift
!
!   shift a time series in time
!
    implicit none
    include 'param.inc'
    integer nshift,nperyear
    real data(npermax,yrbeg:yrend)
    character file*256,string*10,var*40,units*80
    logical lwrite
    integer iargc
    lwrite = .false.

    if ( iargc().lt.2 ) then
        write(0,*) 'usage: timeshift series nmonth'
        call exit(-1)
    endif
    call getarg(1,file)
    call getarg(2,string)
    read(string,*) nshift

    call readseries(file,data,npermax,yrbeg,yrend,nperyear,var,units,.false.,lwrite)
    call shiftseries(data,npermax,nperyear,yrbeg,yrend,nshift)
    call copyheader(file,6)
    call nperyear2units(nperyear,units)
    if ( abs(nshift).eq.1 ) then
        write(6,'(a,i4,a)') '# time shifted by ',nshift,' '//trim(units)
    elseif ( abs(nshift).gt.1 ) then
        write(6,'(a,i4,a)') '# time shifted by ',nshift,' '//trim(units)//'s'
    endif
    call printdatfile(6,data,npermax,nperyear,yrbeg,yrend)
end program
