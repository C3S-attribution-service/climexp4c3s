program coordinates2kml

!   convert the Climate Explorer coordinates to KML coordinates

    implicit none
    integer :: i
    real :: lon1,lon2,lat1,lat2,maplon1,maplon2,maplat1,maplat2
    real :: lookat_lon,lookat_lat,lookat_range
    logical :: lwrite,lwrap
    character :: string*80

    lwrite = .false. 
    if ( command_argument_count() < 4 ) then
        print *,'usage: coordinates2kml lon1|@ lon2|@ lat1|@ lat2|@ [mapfile]'
        call exit(-1)
    endif
    call getcoordinate(1,lon1)
    call getcoordinate(2,lon2)
    call getcoordinate(3,lat1)
    call getcoordinate(4,lat2)
    if ( lwrite ) print *,'found arguments ',lon1,lon2,lat1,lat2
    if ( command_argument_count() >= 5 ) then
!       I expect in this file the grads copmmands to restrict the
!       plot area, i.e., 'set lon l1 l2' and/or 'set lat l1 l2'
        call get_command_argument(5,string)
        open(1,file=string,status='old',err=902)
        maplon1 = 3e33
        maplon2 = 3e33
        maplat1 = 3e33
        maplat2 = 3e33
        200 continue
        read(1,'(a)',end=210) string
        i = index(string,' lon ')
        if ( i /= 0 ) then
            read(string(i+5:),*) maplon1,maplon2
            if ( lwrite ) print *,'found maplons ',maplon1,maplon2
        endif
        i = index(string,' lat ')
        if ( i /= 0 ) then
            read(string(i+5:),*) maplat1,maplat2
            if ( lwrite ) print *,'found maplats ',maplat1,maplat2
        endif
        goto 200
        210 continue
        close(1)
!           the ordering in getfieldopts.cgi and grads.cgi is
!           1) user-defined, if not there
!           2) from map, if not there
!           3) -180,180 -90,90
        call overwriteifundefined(lon1,maplon1)
        call overwriteifundefined(lon2,maplon2)
        call overwriteifundefined(lat1,maplat1)
        call overwriteifundefined(lat2,maplat2)
    endif
    call overwriteifundefined(lon1,-180.)
    call overwriteifundefined(lon2,+180.)
    call overwriteifundefined(lat1,-90.)
    call overwriteifundefined(lat2,+90.)

!   whole world?

    if ( lon2-lon1 == 360 ) then
        lwrap = .true. 
    else
        lwrap = .false. 
    endif

!       view point

    lookat_lat = (lat1+lat2)/2
    lookat_lon = (lon1+lon2)/2
    lookat_range = max(10000.,min(10000000.,100000*max(abs(lon2-lon1),abs(lat2-lat1))))

!       and output for grads.cgi

    if ( lwrap ) print '(a)','lwrap=TRUE'
    print '(a,f8.3,a)','lookat_lon="',lookat_lon,'"'
    print '(a,f8.3,a)','lookat_lat="',lookat_lat,'"'
    print '(a,f9.0,a)','lookat_range="',lookat_range,'"'
    print '(a,f8.3,a)','latlonbox_north="',lat2,'"'
    print '(a,f8.3,a)','latlonbox_south="',lat1,'"'
    print '(a,f8.3,a)','latlonbox_east="',lon2,'"'
    print '(a,f8.3,a)','latlonbox_west="',lon1,'"'

!       error messages

    goto 999
902 write(0,*) 'coordinates2kml: error opening mapfile ',trim(string)
    call exit(-1)
999 continue
end program coordinates2kml

subroutine getcoordinate(i,x)
    implicit none
    integer :: i
    real :: x
    character string*80
    call get_command_argument(i,string)
    if ( string(1:1) /= '@' ) then
        read(string,*,err=901) x
    else
        x = 3e33
    endif
    return
901 write(0,*) 'coordinates2kml: error reading coordinate from ',trim(string)
    call exit(-1)
end subroutine getcoordinate

subroutine overwriteifundefined(x,y)
    implicit none
    real :: x,y
    if ( x > 1e33 ) x = y
end subroutine overwriteifundefined
