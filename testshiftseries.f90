program testshiftseries
!
!   test shiftseries_yr
!
    implicit none
    include 'param.inc'
    integer dyr,nperyear
    real data(366,yrbeg:yrend)
    logical lwrite
    character file*1024,var*20,units*40
    
    print *,'dyr?'
    read *,dyr
    call get_command_argument(1,file)
    call readseries(trim(file),data,366,yrbeg,yrend,nperyear,var,units,.false.,.false.)
    lwrite = .true.
    call shiftseriesyear(data,366,nperyear,yrbeg,yrend,dyr,lwrite)
    !!!call printdatfile(6,data,366,nperyear,yrbeg,yrend)
end program