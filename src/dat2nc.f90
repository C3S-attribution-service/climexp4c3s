program dat2nc
!
!   program to convert a .dat-style time series into a netCDF file 
!   using the climexp libraries
!
    implicit none
    include "param.inc"
    integer             :: yr,mo,nperyear,i,j,n,mens,mens1
    real*4              :: x,y,x0,y0
    real*4,allocatable  :: data(:,:,:)
    character           :: file*256,ncfile*256,since*32,var*40,lvar*120,svar*120,description*2000
    character           :: type*1,comment*10000,line*80,units*20,metadata(2,100)*2000
    character           :: history*20000,title*200
    logical             :: lwrite,lstandardunits,lexist
    integer,external    :: leap
    lwrite = .false.
!
!   read data
!
    if ( command_argument_count().ne.4 ) then
        print *,'usage: dat2nc infile type name outfile ',command_argument_count()
        stop
    endif
    call get_command_argument(1,file)
    if ( lwrite ) print *,'reading from ',trim(file)
    ! make sure we do not accidentally read from an old netcdf file :-(
    i = index(file,'.dat')
    if ( i /= 0 ) then
        ncfile = file(:i)//'nc'
        inquire(file=trim(ncfile),exist=lexist)
        if ( lexist ) then
            open(1,file=trim(ncfile))
            close(1,status='delete')
        end if
    end if
    ! remove output file if it exists
    call get_command_argument(4,ncfile)
    inquire(file=trim(ncfile),exist=lexist)
    if ( lexist ) then
        open(1,file=trim(ncfile))
        close(1,status='delete')
    end if
    allocate(data(npermax,yrbeg:yrend,0:nensmax))
    lstandardunits = .false.
    call readensseriesmeta(file,data,npermax,yrbeg,yrend,nensmax,nperyear,mens1,mens, &
        var,units,lvar,svar,history,metadata,lstandardunits,lwrite)
    if ( var.eq.' ' ) then
        call get_command_argument(2,file)
        if ( file.eq.'i') then
            var = 'index'
        elseif ( file.eq.'t' ) then
            var = 'temperature'
        elseif ( file.eq.'p' ) then
            var = 'precipitation'
        elseif ( file.eq.'s' ) then
            var = 'pressure'
        elseif ( file.eq.'l' ) then
            var = 'sealevel'
        else
            var = 'unknown'
        endif
    endif
    call get_command_argument(3,description)
    do i=1,len(description)
        if ( description(i:i).eq.'_') description(i:i) = ' '
    enddo
!
!   and write out
!
    comment = ' '
    if ( mens1 == mens ) then
        call writencseries(ncfile,data,npermax,yrbeg,yrend,nperyear, &
            title,description,comment,metadata,history,var,lvar,units)
    else
        call writencensseries(ncfile,data,npermax,yrbeg,yrend,nperyear, &
            title,description,comment,metadata,history,var,lvar,units,mens1,mens)
    end if
end program
