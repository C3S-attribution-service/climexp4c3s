program timeshift
!
!   shift a time series in time
!
    implicit none
    include 'param.inc'
    integer :: nshift,nperyear
    real :: data(npermax,yrbeg:yrend)
    character :: file*1023,string*10,var*80,units*120,svar*120,lvar*120, &
        history*10000,metadata(2,100)*2000
    logical :: lwrite
    lwrite = .false.

    if ( command_argument_count().lt.2 ) then
        write(0,*) 'usage: timeshift series nmonth'
        call exit(-1)
    endif
    call get_command_argument(1,file)
    call get_command_argument(2,string)
    read(string,*) nshift

    call readseriesmeta(file,data,npermax,yrbeg,yrend,nperyear,var,units,lvar,svar, &
        history,metadata,.false.,lwrite)
    call shiftseries(data,npermax,nperyear,yrbeg,yrend,nshift)
    call copyheader(file,6)
    call printmetadata(6,file,' ',' ',history,metadata)
    if ( svar /= ' ' ) then
        print '(2a)','# variable_standard_name :: ',trim(svar)
    end if
    call nperyear2units(nperyear,units)
    if ( abs(nshift).eq.1 ) then
        write(6,'(a,i4,a)') '# time shifted by ',nshift,' '//trim(units)
    elseif ( abs(nshift).gt.1 ) then
        write(6,'(a,i4,a)') '# time shifted by ',nshift,' '//trim(units)//'s'
    endif
    call printdatfile(6,data,npermax,nperyear,yrbeg,yrend)
end program
