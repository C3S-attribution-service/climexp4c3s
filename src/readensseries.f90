subroutine readensseries(file,data,npermax,yrbeg,yrend,nensmax &
        ,nperyear,mens1,mens,var,units,lstandardunits,lwrite)
    implicit none
    integer :: npermax,yrbeg,yrend,nperyear,nensmax,mens1,mens
    real :: data(npermax,yrbeg:yrend,0:nensmax)
    logical :: lstandardunits,lwrite
    character :: file*(*),var*(*),units*(*)
    character ::  lvar*120,svar*120,history*50000,metadata(2,100)*1000

    call readensseriesmeta(file,data,npermax,yrbeg,yrend,nensmax &
        ,nperyear,mens1,mens,var,units,lvar,svar,history,metadata,lstandardunits,lwrite)

end subroutine readensseries

subroutine readensseriesmeta(file,data,npermax,yrbeg,yrend,nensmax &
    ,nperyear,mens1,mens,var,units,lvar,svar,history,metadata,lstandardunits,lwrite)

!   read a data file or an ensemble of data files in the array data
!   also returns the number of ensemble members in mens, 0=not an ensemble

    implicit none
    include 'netcdf.inc'
    integer :: npermax,yrbeg,yrend,nperyear,nensmax,mens1,mens
    real :: data(npermax,yrbeg:yrend,0:nensmax)
    logical :: lstandardunits,lwrite
    character :: file*(*),var*(*),units*(*)
    character ::  lvar*(*),svar*(*),history*(*),metadata(2,100)*(*)
    integer :: i,m,iens,status,ncid,unit
    logical :: lexist,lfirst
    character :: ensfile*1023,line*1023,saveunits*100

    if ( lwrite ) then
        print *,'readensseries: file = ',trim(file)
        print *,'               npermax,yrbeg,yrend = ',npermax,yrbeg,yrend
    endif

!   ensembles are denoted by a '%' or '++' in the file name for separate file, @@ if the
!   ensemble dimension is in the netcdf file (no other formats are supported)

    if ( index(file,'@@') /= 0 ) then
        ! it is a file with an ensemble dimension
        i = index(file,'.dat')
        if ( i > 0 ) then
            file(i:) = '.nc'
        end if
        inquire(file=trim(file),exist=lexist)
        if ( lexist ) then ! netcdf file
            if ( lwrite ) print *,'calling readncseriesensmeta ',trim(file)
            call readncseriesensmeta(file,data,npermax,nperyear,yrbeg,yrend,nensmax,mens1,mens &
                ,ncid,var,units,lvar,svar,history,metadata,lwrite)
        else ! only an ascii file
            call rsunit(unit)
            file(i:) = '.dat'
            if ( lwrite ) print *,'readensseriesmeta: opening ',trim(file), &
                ' as ascii file on unit ',unit
            open(unit,file=trim(file),status='old',err=902)
            call readseriesheader(var,units,lvar,svar,history,metadata,line,unit,mens1,lwrite)
            if ( lwrite) print *,'readensseriesmeta: read metadata, mens1 = ',mens1
            m = mens1
            data = 3e33
            do while ( m >= 0 )
                ! reads one ensemble member and returns the number of the next one
                call readoneseries(unit,file,line,data(1,yrbeg,mens),npermax,yrbeg,yrend,nperyear,m,lwrite)
                ! make sure line contains the next numerical data
            101 continue
                if ( line == ' ' .or. line(1:1) == '#' .or. line(2:2) == '#' ) then
                    read(unit,'(a)',err=102,end=102) line
                    goto 101
                end if
                go to 103
            102 continue
                m = -1
            103 continue
                if ( lwrite) print *,'readensseriesmeta: read ensemble member ',mens,', next one is ',m
                if ( m >= 0 ) mens = m
            end do 
        end if 
        if ( lstandardunits ) then
            saveunits = units
            do iens=mens1,mens
                units = saveunits
                call makestandardseries(data(1,yrbeg,iens),npermax,yrbeg,yrend, &
                    nperyear,var,units,lwrite)
            end do
        end if
    else if ( index(file,'%') == 0 .and. index(file,'++') == 0 ) then
        mens = 0
        mens1 = 0
        if ( lwrite ) print *,'not an ensemble, calling readseriesmeta'
        call readseriesmeta(file,data(1,yrbeg,0),npermax,yrbeg,yrend &
            ,nperyear,var,units,lvar,svar,history,metadata,lstandardunits,lwrite)
    else
        mens = -1
        mens1 = 0
        lfirst = .true.
        do iens=0,nensmax
            ensfile = file
            call filloutens(ensfile,iens)
            inquire(file=trim(ensfile),exist=lexist)
            if ( .not. lexist ) then
                if ( mens == -1 ) then
                    mens1 = iens + 1
                end if
                if ( iens > 1 ) exit
            else
                mens = iens
            end if
        end do
        if ( mens < 0 ) then
            write(0,*) 'readensseries: error: could not find ensemble ' &
                ,trim(file),' with member ',trim(ensfile)
            call exit(-1)
        end if
        do iens=mens1,mens
            ensfile = file
            call filloutens(ensfile,iens)
            inquire(file=trim(ensfile),exist=lexist)
            if ( .not. lexist ) then
                if ( lwrite ) print *,'skipping file ',trim(ensfile)
                data(:,:,iens) = 3e33
            else
                call keepalive1('Reading ensemble member',iens,mens)
                if ( lfirst ) then
                    lfirst = .false.
                    if ( lwrite ) print *,'reading file ',trim(ensfile),' with metadata'
                    call readseriesmeta(ensfile,data(1,yrbeg,iens),npermax &
                        ,yrbeg,yrend,nperyear,var,units,lvar,svar,history,metadata &
                        ,lstandardunits,lwrite)
                else
                    if ( lwrite ) print *,'reading file ',trim(ensfile)
                    call readseries(ensfile,data(1,yrbeg,iens),npermax &
                        ,yrbeg,yrend,nperyear,var,units,lstandardunits,lwrite)
                end if
            end if
        end do
    end if
    return
902 print *,'readensseries: error opening ',trim(file)
    call exit(-1)
end subroutine readensseriesmeta
