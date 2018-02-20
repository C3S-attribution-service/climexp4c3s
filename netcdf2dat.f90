program netcdf2dat

!   convert a netCDF time series to a standard .dat file
!
    implicit none
    include 'netcdf.inc'
    include 'param.inc'
    integer :: npermx
    parameter(npermx=24*366)
    character :: file*1023,var*80,units*80,lvar*120,svar*120,history*10000,metadata(2,100)*2000
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
    call readseriesmeta(file,data,npermx,yrbeg,yrend,nperyear,var,units,lvar,svar, &
        history,metadata,lstandardunits,lwrite)
    call copyheader(file,6)
    call printmetadata(6,file,' ',' ',history,metadata)
    if ( svar /= ' ' ) then
        print '(2a)','# variable_standard_name :: ',trim(svar)
    end if
    
    j = 0
    call printdatfile(6,data,npermx,nperyear,yrbeg,yrend)
end program netcdf2dat