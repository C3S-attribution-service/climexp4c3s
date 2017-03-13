program netcdf2dat

!   convert a netCDF time series to a standard .dat file
!
    implicit none
    include 'netcdf.inc'
    include 'param.inc'
    integer :: npermx
    parameter(npermx=24*366)
    character file*255,var*60,units*80
    integer :: i,j,nperyear,status,ncid
    real :: data(npermx,yrbeg:yrend)
    logical :: lwrite,lstandardunits
    integer :: iargc

    lwrite = .false. 
    lstandardunits = .false. 
    if ( iargc() < 1 ) then
        print *,'usage: netcdf2dat infile'
        call exit(-1)
    endif
    call getarg(1,file)
    call readseries(file,data,npermx,yrbeg,yrend,nperyear,var,units,lstandardunits,lwrite)
    call copyheader(file,6)
    j = 0
    do i=len_trim(file),1,-1
        if ( file(i:i) == '/' ) j = j + 1
        if ( j == 2 ) exit
    end do
    print '(2a)','# converted from ',trim(file(i+1:))
    call printdatfile(6,data,npermx,nperyear,yrbeg,yrend)
end program netcdf2dat