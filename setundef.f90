program setundef

!   set the values of series 1 that are undefined in series 2 to undefined.
!   ((like cdo)

    implicit none
    include 'param.inc'
    integer :: yr,mo,nperyear,n
    real :: series1(npermax,yrbeg:yrend),series2(npermax,yrbeg:yrend)
    logical :: lstandardunits,lwrite
    character :: file*1024,var*40,lvar*80,units*40

    if ( command_argument_count() /= 2 ) then
        write(0,*) 'usage: setundef series1 series2'
        call exit(-1)
    end if
    lwrite = .false.
    lstandardunits = .false.
    call get_command_argument(2,file)
    call readseries(file,series2,npermax,yrbeg,yrend,n,var,units,lstandardunits,lwrite)
    call get_command_argument(1,file)
    call readseries(file,series1,npermax,yrbeg,yrend,nperyear,var,units,lstandardunits,lwrite)
    if ( n /= nperyear ) then
        write(0,*) 'setundef: error: unequal time scales ',n,nperyear
        call exit(-1)
    end if
    do yr=yrbeg,yrend
        do mo=1,nperyear
            if ( series2(mo,yr).gt.1e33 ) then
                series1(mo,yr) = 3e33
            end if
        end do
    end do
    call copyheader(file,6)
    call printdatfile(6,series1,npermax,nperyear,yrbeg,yrend)
end program