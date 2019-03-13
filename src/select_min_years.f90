program select_min_years
!
!   select stations with a minimum number of years from a stationlist
!   or min/max elevastion, or min/max lngitude/latitude
!
    implicit none
    integer :: minyrs,i,j,n
    real :: xval,elev,lon,lat,elevmin,elevmax,lon1,lon2,lat1,lat2
    logical :: lok
    character :: file*1024,lines(10)*256,string*80,value*80

    call get_command_argument(2,string)
    if ( string == ' ' ) then
        write(0,*) 'usage: select_min_years inlist minyrs [elevmin N] [elevmax M] [lon{12} degrees] [lat{12} degrees]'
        call exit(-1)
    end if
    read(string,*) minyrs
    elevmin = 3e33
    elevmax = 3e33
    lon1 = 3e33
    lon2 = 3e33
    lat1 = 3e33
    lat2 = 3e33
    do i=3,command_argument_count(),2
        call get_command_argument(i,string)
        call get_command_argument(i+1,value)
        if ( value == ' ' ) exit
        write(0,*) 'string,value = ',trim(string),',',trim(value)
        read(value,*) xval
        if ( string(1:7) == 'elevmin' ) then
            elevmin = xval
        else if ( string(1:7) == 'elevmax' ) then
            elevmax = xval
        else if ( string(1:4) == 'lon1' ) then
            lon1 = xval
        else if ( string(1:4) == 'lon2' ) then
            lon2 = xval
        else if ( string(1:4) == 'lat1' ) then
            lat1 = xval
        else if ( string(1:4) == 'lat2' ) then
            lat2 = xval
        else
            write(0,*) 'select_min_years: error: unknown option ',trim(string),' ',trim(value)
        end if
    end do
    call get_command_argument(1,file)
    open(1,file=trim(file),status='old')
    read(1,'(a)') lines(1) ! located...
    i = index(lines(1),'ocated ') + 7
    if ( lines(1)(i:i) /= 's' .and. lines(1)(i:i) /= 'S' ) then
        ! delete number of stations
        j = i
        do while ( lines(1)(j:j) == ' ' )
            j = j + 1
        end do
        do while ( lines(1)(j:j) /= ' ' )
            j = j + 1
        end do
        lines(1)(i:) = lines(1)(j:)
    end if
    write(*,'(a)') trim(lines(1)) 
    do
        lok = .true.
        do i=1,10
            read(1,'(a)',end=800) lines(i)
            j = index(lines(i),'oordinates')
            if ( j > 0 ) then
                call parsecoord(lines(i),lat,lon,elev)
                if ( elevmin < 1e33 .and. elev < elevmin ) lok = .false.
                if ( elevmax < 1e33 .and. elev > elevmax ) lok = .false.
                if ( lon1 < 1e33 .and. lon < lon1 ) lok = .false.
                if ( lon2 < 1e33 .and. lon > lon2 ) lok = .false.
                if ( lat1 < 1e33 .and. lat < lat1 ) lok = .false.
                if ( lat2 < 1e33 .and. lat > lat2 ) lok = .false.
            end if
            j = index(lines(i),'Found ') + index(lines(i),'found ')
            if ( j > 0 ) then
                read(lines(i)(j+6:),*) n
                if ( n >= minyrs .and. lok ) then
                    do j=1,i
                        write(*,'(a)') trim(lines(j))
                    end do
                end if
                exit
            end if
        end do
    end do
800 continue
    do j=1,i-1
        write(*,'(a)') trim(lines(j))
    end do
end program select_min_years