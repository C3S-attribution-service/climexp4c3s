program compute_wetbulb
!
!   compute wet bulb temperatures at a station
!
    implicit none
    integer npermax,yrbeg,yrend
    parameter(npermax=366*24,yrbeg=1850,yrend=2100)
    integer hh,dd,mm,yyyy,i,j,nperyear,nperyear2
    real,allocatable :: temp(:,:),dewtemp(:,:),pres(:,:),wetbulb(:,:)
    logical lwrite,lstandardunits,lsimplified
    character tempfile*255,dewtempfile*255,presfile*255,var*40,units*80,longname*80,line*255
    integer iargc
    real,external :: wetbulbdew,wetbulbsim

    lwrite = .false.
    if ( iargc() /= 3 ) then
        write(0,*) 'usage: compute_wetbulb tempfile dewtempfile pressurefile|simplified'
        call exit(-1)
    end if
    allocate(temp(npermax,yrbeg:yrend),dewtemp(npermax,yrbeg:yrend))
    allocate(pres(npermax,yrbeg:yrend),wetbulb(npermax,yrbeg:yrend))
    call getarg(1,tempfile)
    call getarg(2,dewtempfile)
    call getarg(3,presfile)
    open(1,file=trim(tempfile),status='old')
    do
        read(1,'(a)') line
        if ( line(1:1) /= '#' .and. line(2:2) /= '#' ) exit
        if ( index(line,'href') /= 0 .or. index(line,'DATA') /= 0 ) print '(a)',trim(line)
    end do
    close(1)
    lstandardunits = .false.
    write(0,*) 'reading ',trim(tempfile)
    call readseries(tempfile,temp,npermax,yrbeg,yrend,nperyear,var,units,lstandardunits,lwrite)
    if ( units == 'K' ) then
        temp = temp - 273.15
        units = 'Celsius'
    end if
    if ( units /= 'Celsius' ) then
        write(0,*) 'compute_wetbulb: error: expecting Celsius or K for temperature, not ',trim(units)
        call exit(-1)
    end if
    write(0,*) 'reading ',trim(dewtempfile)
    call readseries(dewtempfile,dewtemp,npermax,yrbeg,yrend,nperyear2,var,units,lstandardunits,lwrite)
    if ( nperyear /= nperyear2 ) then
        write(0,*) 'compute_wetbulb: error: unequal time steps ',nperyear,nperyear2
        call exit(-1)
    end if
    if ( units == 'K' ) then
        dewtemp = dewtemp - 273.15
        units = 'Celsius'
    end if
    if ( units /= 'Celsius' ) then
        write(0,*) 'compute_wetbulb: error: expecting Celsius or K for dewpoint temperature, not ',trim(units)
        call exit(-1)
    end if
    if ( presfile.eq.'simplified' ) then
        lsimplified = .true.
        do i=yrbeg,yrend
            do j=1,nperyear
                wetbulb(j,i) = wetbulbsim(temp(j,i),dewtemp(j,i),lwrite)
            end do
        end do
        print '(a)','# simplified wet bulb temperature computed from temperature and dew point temperature'
        print '(a)','# Twet [Celsius] simplified wet bulb temperature'        
    else
        lsimplified = .false.
    write(0,*) 'reading ',trim(presfile)
    call readseries(presfile,pres,npermax,yrbeg,yrend,nperyear2,var,units,lstandardunits,lwrite)
        if ( units == 'Pa' ) then
            pres = pres/100
            units = 'hPa'
        end if
        if ( units /= 'hPa' .and. units /= 'mb' ) then
            write(0,*) 'compute_wetbulb: error: expecting Pa, hPa or mb for dewpoint temperature, not ',trim(units)
            call exit(-1)
        end if
        if ( nperyear /= nperyear2 ) then
            write(0,*) 'compute_wetbulb: error: unequal time steps ',nperyear,nperyear2
            call exit(-1)
        end if
        write(0,*) 'Computing'    
        do i=yrbeg,yrend
            do j=1,nperyear
                wetbulb(j,i) = wetbulbdew(temp(j,i),pres(j,i),dewtemp(j,i),lwrite)
            end do
        end do
        print '(a)','# wet bulb temperature computed from temperature, dew point temperature and sea-level pressure'
        print '(a)','# Twet [Celsius] wet bulb temperature'
    end if
    call printdatfile(6,wetbulb,npermax,nperyear,yrbeg,yrend)
end program