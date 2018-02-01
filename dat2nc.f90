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
    character           :: file*256,ncfile*256,since*32,var*40,lvar*80,description*2000
    character           :: type*1,comment*10000,line*80,units*20,metadata(2,100)*2000
    character           :: history*20000,title*200
    logical             :: lwrite,lstandardunits,lexist
    integer             :: iargc
    integer,external    :: llen,leap
    lwrite = .FALSE.
!
!   read data
!
    if ( iargc().ne.4 ) then
        print *,'usage: dat2nc infile type name outfile ',iargc()
        stop
    endif
    call getarg(1,file)
    call getarg(4,ncfile)
    allocate(data(npermax,yrbeg:yrend,0:nensmax))
    if ( lwrite ) print *,'reading from ',trim(file)
    lstandardunits = .false.
    ! make sure we do not accidentally read from an old netcdf file :-(
    inquire(file=trim(ncfile),exist=lexist)
    if ( lexist ) then
        open(1,file=trim(ncfile))
        close(1,status='delete')
    end if
    call readensseries(file,data,npermax,yrbeg,yrend,nensmax,nperyear,mens1,mens, &
    &   var,units,lstandardunits,lwrite)
    open(1,file=file)
    comment = ' '
    lvar = ' '
    metadata = ' '
    history = ' '
    n = 0
    do i=1,100
        read(1,'(a)') description
        if ( description(1:1) /= '#' ) exit
        j = index(description,' :: ')
        if ( j /= 0 ) then
            if ( description(3:j-1) == 'title' .or. description(3:j-1) == 'Title' ) then
                title = description(j+4:)
            else if ( description(3:j-1) == 'history' .or. description(3:j-1) == 'History' ) then
                history = description(j+4:)
            else
                n = n + 1
                metadata(1,n) = description(3:j-1)
                metadata(2,n) = description(j+4:)
            end if
        else
            j = index(description,']')
            if ( j.eq.0 ) then
                comment(1+len_trim(comment):) = ', '//trim(description(3:))
            else
                lvar = description(j+2:)
                j = index(description,'[')
                var = description(3:j-1)
            end if
        end if
    end do
    close(1)
    comment = comment(3:) ! get rid of the leading ", "
    comment(1+llen(comment):) = ', via the KNMI Climate Explorer (http://climexp.knmi.nl)'
    if ( var.eq.' ' ) then
        call getarg(2,file)
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
    call getarg(3,description)
    do i=1,len(description)
        if ( description(i:i).eq.'_') description(i:i) = ' '
    enddo
!
!   and write out
!
    if ( mens1 == mens ) then
        call writencseries(ncfile,data,npermax,yrbeg,yrend,nperyear &
        &   ,title,description,comment,metadata,history,var,lvar,units)
    else
        write(0,*) 'dat2nc: error: writing netcdf ensembles not yet ready'
        write(*,*) 'dat2nc: error: writing netcdf ensembles not yet ready'
    end if
end program
