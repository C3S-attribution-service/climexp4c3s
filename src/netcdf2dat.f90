program netcdf2dat

!   convert a netCDF time series to a standard .dat file
!
    implicit none
    include 'netcdf.inc'
    include 'param.inc'
    integer :: npermx
    parameter(npermx=366)
    character :: file*1023,var*80,units*80,lvar*120,svar*120,history*50000,metadata(2,100)*2000
    integer :: i,j,nperyear,status,ncid,mens1,mens,iens
    real,allocatable :: data(:,:,:)
    logical :: lwrite,lstandardunits

    allocate(data(npermx,yrbeg:yrend,0:nensmax))
    lwrite = .false. 
    lstandardunits = .false. 
    if ( command_argument_count() < 1 ) then
        print *,'usage: netcdf2dat infile'
        call exit(-1)
    endif
    call get_command_argument(1,file)
    if ( lwrite ) print *,'calling readensseriesmeta'
    call readensseriesmeta(file,data,npermx,yrbeg,yrend,nensmax, &
        nperyear,mens1,mens,var,units,lvar,svar,history,metadata,lstandardunits,lwrite)
    call printmetadata(6,file,' ',' ',history,metadata)
    if ( svar /= ' ' ) then
        print '(2a)','# variable_standard_name :: ',trim(svar)
    end if
    
    j = 0
    do iens=mens1,mens
        if ( mens1 /= mens ) then
            print '(a,i4)','# ensemble member ',iens
            print '(a)'
        end if
        call printdatfile(6,data(1,yrbeg,iens),npermx,nperyear,yrbeg,yrend)
    end do
    deallocate(data)
end program netcdf2dat