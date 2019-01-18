    program list2kml

!       convert a simple list produced by stationlist list to kml for
!       Google Earth

    implicit none
    integer :: i
    real :: lat,lon,height,lat1,lat2,lon1,lon2,dum
    character file*255,email*80,prog*40,extraargs*255,line*255, &
    code*20,name*80,nameu*80,climate*40

    if ( command_argument_count() < 4 ) then
        print *,'usage: list2kml plotlist var email prog [extraargs]'
        call exit(-1)
    endif

    call get_command_argument(1,file)
    call get_command_argument(2,climate)
    do i=1,len_trim(climate)
        if ( climate(i:i) == '_' ) climate(i:i) = ' '
    enddo
    call get_command_argument(3,email)
    call get_command_argument(4,prog)
    if ( command_argument_count() == 5 ) then
        call get_command_argument(5,extraargs)
    else
        extraargs = ' '
    end if
    open(1,file=file,status='old')
    print '(a)','<?xml version="1.0" encoding="UTF-8"?>'
    print '(a)','<kml xmlns="http://earth.google.com/kml/2.0">'
    print '(a)','<Document>'
    print '(3a)','<name>',trim(file),'</name>'
    10 continue
    read(1,'(a)') line
    if ( line(1:1) /= '#' ) goto 10
    read(line(3:),*) lon1,lon2,lat1,lat2
    print '(a)','<LookAt>'
    if ( abs(lon2-lon1-360) < 1 ) then
    !           whole earth, let's be eurocentric
        print '(a)','<longitude>0</longitude>'
        print '(a)','<latitude>0</latitude>'
        print '(a)','<range>10000000</range>'
    else
        print '(a,f8.3,a)','<longitude>',(lon1+lon2)/2,'</longitude>'
        print '(a,f8.3,a)','<latitude>',(lat1+lat2)/2,'</latitude>'
        print '(a,f14.3,a)','<range>',max(10000., &
        &            100000*max(abs(lon2-lon1),abs(lat2-lat1))),'</range>'
    endif
    print '(a)','</LookAt>'
    100 continue
    read(1,'(a)',end=800) line
    if ( line(1:1) == '#' .or. line == ' ' ) goto 100
    i = index(line,' ')
    code = line(:i-1)
    read(line(i+1:),'(2f10.4,g11.4,f8.5,a)') lon,lat,height,dum,name
    nameu = name
    do i=1,len_trim(name)
        if ( name(i:i) == '_' ) name(i:i) = ' '
    enddo
    print '(a)','<Placemark>'
    print '(3a)','<name>',trim(name),'</name>'
    print '(a)','<Point><coordinates>'
    print '(3(f10.4,a))',lon,',',lat,',',0.
    print '(a)','</coordinates></Point>'
    print '(a)','<description><![CDATA['
!**        print '(3a)','<h2>',trim(name),'</h2>'
    if ( extraargs /= ' ' ) then
        print '(13a)','<a href="http://climexp.knmi.nl/',trim(prog) &
            ,'.cgi?id=',trim(email),'&WMO=',trim(code),'&STATION=', &
            trim(nameu),'&extraargs=',trim(extraargs),'">', &
            trim(climate),'</a>'
    else
        print '(13a)','<a href="http://climexp.knmi.nl/',trim(prog) &
            ,'.cgi?id=',trim(email),'&WMO=',trim(code),'&STATION=', &
            trim(nameu),'">', &
            trim(climate),'</a>'        
    end if
    print '(a)',']]></description>'
    print '(a)','</Placemark>'
    goto 100
800 continue
    print '(a)','</Document>'
    print '(a)','</kml>'
end program list2kml
