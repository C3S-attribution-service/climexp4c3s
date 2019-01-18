real function in_polygon(polygon,npol,x,y,pole,lwrite)
    !
    ! checks whether the point (x,y) is inside the polygon
    ! if pole == 'sp' then the South Pole is inside the polygon, otherwise outside
    ! returns 0 is outside, 1 if inside, 0.5 if on the edge and 3e33 if it cannot be determined
    ! (the overlying routine will try again perturbing (x,y) slightly to avoid the ambiguity)
    !
    implicit none
    integer npol
    double precision polygon(2,npol)
    double precision x,y
    character pole*(*)
    logical lwrite
    integer ipol
    double precision ynull
    real result,partial
    real segments_crossed
    logical foundtouch,llwrite

    foundtouch = .false.
    ynull = -90
    if ( pole == 'sp') then
        result = 1
    else
        result = 0
    end if
    if ( y == ynull ) then
        in_polygon = result
        return
    end if
    !
    if ( .false. .and. lwrite ) print *,'started with ',result
    do ipol=1,npol-1
        ! skip "segments" that include the polygon-separating marker 3e33,3e33
        if ( polygon(1,ipol) > 1e33 ) cycle
        if ( polygon(1,ipol+1) > 1e33 ) cycle
        llwrite = .false.
        !!!if ( x == 14.75 .and. y > 48.5 .and. &
        !!!& (polygon(1,ipol)-x)*(polygon(1,ipol+1)-x) < 0 ) llwrite = .true.
        partial = segments_crossed(polygon(1:2,ipol),polygon(1:2,ipol+1),x,y,ynull,llwrite)
        if ( lwrite .and. partial /= 0 ) print '(a,i6,f3.0,4f10.5)','in_polygon: segments_crossed returns ', &
            ipol,partial,polygon(1:2,ipol),polygon(1:2,ipol+1)
        if ( partial == 0.5 ) then
            result = partial
            if ( lwrite ) print *,'end point on line ',result
            exit
        end if
        if ( partial == 1 ) then
            result = 1-result
            if ( .false. .and. lwrite ) print *,'flipped from ',1-result,' to ',result
        else if ( partial == 3e33 ) then
            if ( lwrite ) then
                print *,'undefined ',result
                partial = segments_crossed(polygon(1:2,ipol),polygon(1:2,ipol+1),x,y,ynull,lwrite)
            end if
            foundtouch = .true.
        end if
    end do
    if ( foundtouch .and. result /= 0.5 ) then
        in_polygon = 3e33
    else
        in_polygon = result
    end if
end function in_polygon

real function segments_crossed(p1,p2,xin,yin,ynull,lwrite)
    !
    ! determines whther the two line segments p1 - p2 and (x,ynull) - (x,y) cross or not
    ! 0: no cross
    ! 1: crossed inside line segment
    ! 0.5: touched at (x,y)
    ! 3e33: touched at another point
    !
    implicit none
    double precision p1(2),p2(2)
    double precision xin,yin,ynull
    logical lwrite
    double precision x,y,x1,y1,x2,y2,lam,mu
    !
    if ( lwrite ) then
        print *,'segments_crossed: input:'
        print *,'p1     = ',p1
        print *,'p2     = ',p2
        print *,'x,y,y0 = ',xin,y,ynull
    end if
    if ( abs(ynull) /= 90 ) then
        write(0,*) 'segments_crossed: error: ynull should be +/-90, not ',ynull
        call exit(-1)
    end if
    !
    ! first reduce the longitude to the shortest possible form assuming only one wrap off
    !
    x1 = p1(1)
    x2 = p2(1)
    x = xin
    y = yin
    if ( abs(x2-x1) > abs(x2-360-x1) .and. abs(x2-360-x1) > 0.01 ) then
        x2 = x2 - 360
    end if
    if ( abs(x2-x1) > abs(x2+360-x1) .and. abs(x2+360-x1) > 0.01 ) then
        x2 = x2 + 360
    end if
    if ( abs(x-(x1+x2)/2) > abs(x-360-(x1+x2)/2) ) then
        x = x - 360
    end if
    if ( abs(x-(x1+x2)/2) > abs(x+360-(x1+x2)/2) ) then
        x = x + 360
    end if
    y1 = p1(2)
    y2 = p2(2)
    !
    ! a few special cases that will cause grief later on
    !
    if ( abs(y1) == 90 .and. abs(y2) == 90 ) then
        if ( lwrite ) print *,'segment is one point on pole'
        segments_crossed = 0
        return
    end if
    if ( lwrite ) then
        print *,'special case? x1,x2,x = ',x1,x2,x
    end if
    if ( abs(x1-x2) < 0.00000001 ) then
        if ( x /= x1 ) then
            segments_crossed = 0
        else if ( ynull == -90 .and. min(y1,y2) > y ) then
            segments_crossed = 0 ! north of the grid point
        else if ( ynull == -90 .and. max(y1,y2) > y ) then
            segments_crossed = 0.5 ! on the line
        else if ( ynull == +90 .and. max(y1,y2) < y ) then
            segments_crossed = 0 ! south of the grid point
        else if ( ynull == +90 .and. min(y1,y2) < y ) then
            segments_crossed = 0.5 ! on the line
        else
            segments_crossed = 3e33 ! no idea whether it crossed or not
        end if
        if ( lwrite ) then
            print *,'special case x1=x2, segments_crossed = ',segments_crossed
        end if
        return
    end if
    if ( y == ynull ) then
        segments_crossed = 0.5
        return
    end if
    ! position along p1 - p2
    lam = (x-x1)/(x2-x1)
    if ( lwrite ) print *,'lam = ',lam
    if ( lam < 0 .or. lam > 1 ) then
        segments_crossed = 0
    else
        ! position along (x,ynull) - (x,y)
        mu = (y1-ynull + lam*(y2-y1))/(y-ynull)
        if ( mu > 1 ) then
            segments_crossed = 0
        else if ( mu == 1 ) then
            segments_crossed = 0.5
        else
            if ( lam == 0 .or. lam == 1 ) then
                segments_crossed = 3e33
            else
                segments_crossed = 1
            end if
        end if      
    end if
end function segments_crossed

subroutine read_polygon(datfile,npol,npolmax,polygon,lwrite)
    ! reads ghe polygons from file
    ! these are assumed to be in the format
    ! # comment
    ! x11 y11
    ! x12 y12
    ! ...
    ! x1n y1n
    !
    ! x21 y21
    ! ...
    ! x2m y2m
    ! ...
    ! where the first index denotes the number of the polygon and the second the 
    ! vertex of the polygon. If the polygon is not closed this is done explicitly.
    ! In the (linear) output a line with 3e33,3e33 denotes the start of a new polygon
    implicit none
    integer npol,npolmax
    double precision polygon(2,npolmax)
    character datfile*(*)
    integer ipol
    character string*100
    logical lwrite
    integer npol1,npol2,i,n
    logical leof
!
    npol = 0
    npol1 = npol + 1
    polygon = 3e33
    if ( lwrite ) print *,'opening ',trim(datfile)
    open(1,file=trim(datfile),status='old')
10  continue
    do
        leof = .true.
        read(1,'(a)',end=100) string
        if ( string == ' ' ) then
            ! signifies the beginning of a new polygon
            leof = .false.
            if ( npol == 0 ) then
                goto 10 ! skip empty first line
            else
                goto 100 ! new polygon
            end if
        end if
        if ( string(1:1) == '#' .or. string(2:2) == '#' ) cycle
        npol = npol + 1
        ! excel often uses commas rather than decimal points
        if ( index(string,'.') == 0 ) then
            n = 0
            do i=1,len(string)
                if ( string(i:i) == ',' ) n = n + 1
            end do
            if ( n == 2 ) then
                do i=1,len(string)
                    if ( string(i:i) == ',' ) string(i:i) = '.'
                end do
            end if
        end if
        read(string,*,err=901,end=901) polygon(1,npol),polygon(2,npol)
        if ( polygon(2,npol) < -90 .or. polygon(2,npol) > 90 ) then
            write(0,*) 'polygon2mask: error: latitude ',npol,' is ',polygon(2,npol)
            call exit(-1)
        end if
    end do
100 continue
    if ( polygon(1,npol1) /= polygon(1,npol) .or. &
         polygon(2,npol1) /= polygon(2,npol) ) then
        ! close loop
        npol = npol + 1
        polygon(1,npol) = polygon(1,npol1)
        polygon(2,npol) = polygon(2,npol1)
    end if
    call keepalive1('Reading shapefile',npol,0)
    if ( lwrite ) then
        print *,'read_polygon: found ',npol-npol1,'-edged polygon'
        do ipol=npol1,npol
            print *,ipol,polygon(1,ipol),polygon(2,ipol)
        end do
    end if
    if ( .not.leof ) then
        ! put marker in list and read another polygon
        npol = npol + 1
        polygon(1,npol) = 3e33
        polygon(2,npol) = 3e33
        ! mark start of new polygon
        npol1 = npol + 1
        ! continue reading
        goto 10
    end if
    close(1)
    return
901 write(0,*) 'read_polygon: error reading two reals from line ',trim(string)
    call exit(-1)
end subroutine read_polygon
