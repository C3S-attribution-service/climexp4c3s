program runningmean

!       program to get the running mean from a data file

    implicit none
    integer :: yrbeg,yrend
    parameter (yrbeg=1700,yrend=2020)
    integer :: i,j,mean,year,nperyear
    real :: data(12,yrbeg:yrend),rdata(12,yrbeg:yrend), &
    sdata(12,yrbeg:yrend)
    character file*128,var*40,units*20
    logical :: lwrite

    lwrite = .false. 
    if ( command_argument_count() < 2 .or. command_argument_count() > 3 ) then
        print *,'usage: runningmean mean infile [sd]'
        call exit(-1)
    endif
    call get_command_argument(1,file)
    read(file,*) mean
    call get_command_argument(2,file)
    call readseries(file,data,12,yrbeg,yrend,nperyear,var,units,lstandardunits,lwrite)
    call copyheader(file,6)
    print '(a,i3,2a)','# with a ',mean,'-month running mean '

    call rmean(data,rdata,sdata,yrbeg,yrend,mean)

    if ( command_argument_count() == 2 ) then
        call printdatfile(6,rdata,12,12,yrbeg,yrend)
    else
        call printdatfile(6,sdata,12,12,yrbeg,yrend)
    endif

end program runningmean