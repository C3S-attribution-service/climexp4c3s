subroutine getlatlonwindow(lat1,lat2,lon1,lon2, &
    xx,nx,xwrap,avex,yy,ny,avey,x1,x2,y1,y2,lwrite)

!   translate the requested window lon1,lon2,lat1,lat2 to grid
!   coordinates (we do not do partial boxes (yet))

    implicit none
    integer :: nx,avex,ny,avey,x1,x2,y1,y2
    real :: lat1,lat2,lon1,lon2,xx(nx),yy(ny)
    logical :: xwrap,lwrite
    integer :: i,j
    real :: lon1c,lon2c,lat1c,lat2c

    if ( lwrite ) then
        print *,'getlatlonwindow: input'
        print *,'lat1,lat2,lon1,lon2    = ',lat1,lat2,lon1,lon2
        print *,'xx(1),xx(2),...,xx(',nx,') = ',xx(1),xx(2),'...',xx(nx),xwrap
        print *,'yy(1),yy(2),...,yy(',ny,') = ',yy(1),yy(2),'...',yy(ny)
        print *,'avex,avey = ',avex,avey
    endif

    if ( lon1 /= 0 .or. lon2 /= 360 ) then
        call getlonwindow(lon1,lon2,x1,x2,xx,nx,lon1c,lon2c,lwrite)
        if ( lwrite ) print *,'getlonwindow returned x1,x2 = ',x1,x2
!       safety first
        if ( xwrap ) then
            if ( x1 < 1 ) then
                x1 = x1 + nx
                x2 = x2 + nx
            endif
            if ( x2 < x1 ) then
                x2 = x2 + nx
            endif
        else
            if ( x1 > nx ) then
                x1 = x1 - nx
                x2 = x2 - nx
                do i=1,nx
                    xx(i) = xx(i) + 360
                enddo
            endif
            x1 = max(1,x1)
            x2 = min(nx,x2)
        endif
    else
        x1 = 1
        x2 = nx
    endif
    i = mod(x2-x1+1,avex)
    if ( i /= 0 .and. .not. (1+x2-x1 == nx .and. xwrap) ) then
!       make it an integer number of intervals
        if ( lwrite ) print *,'getlatlonwindow: adjusting x interval from ',x1,x2
        if ( i <= avex/2 ) then
            j  = x1+i/2
            x2 = x2-i+j-x1
            x1 = j
        else
            i = avex - i
            j = max(1,x1-i/2)
            x2 = min(nx,x2+i-x1+j)
            x1 = j
        endif
        if ( lwrite ) print *,'getlatlonwondow: adjusting x interval to   ',x1,x2
    endif
    if ( x2 < x1 ) then
        write(0,*) 'getlatlonwindow: error: request to average over' &
            //' more X grid points ','(',avex &
            ,') than are available between ',lon1c,' and &
            ',lon2c,'degrees east'
        call exit(-1)
    endif
    if ( lat1 /= -90 .or. lat2 /= +90 ) then
        call getlatwindow(lat1,lat2,y1,y2,yy,ny,lat1c,lat2c,lwrite)
    else
        y1 = 1
        y2 = ny
    endif
    j = mod(y2-y1+1,avey)
    if ( j /= 0 ) then
!       make it an integer number of intervals
        if ( lwrite ) print *,'getlatlonwindow: adjusting y interval from ',y1,y2
        if ( j <= avey/2 ) then
            i  = y1+j/2
            y2 = y2-j+i-y1
            y1 = i
        else
            j = avey - j
            i = max(1,y1-j/2)
            y2 = min(ny,y2+j-y1+i)
            y1 = i
        endif
        if ( lwrite ) print *,'getlatlonwindow: adjusting y interval to   ',y1,y2
    endif
    if ( y2 < y1 ) then
        write(0,*) 'getlatlonwindow: error: request to average over' &
            //' more Y grid points ','(',avey &
            ,') than are available between ',lat1c,' and &
            ',lat2c,' degrees north'
        call exit(-1)
    endif
    if ( lwrite ) print '(a,2i4,a,2i4)','getlatlonwindow: cutting out window x ',x1,x2,' and y ',y1,y2
end subroutine getlatlonwindow

subroutine getlonwindow(lon1,lon2,x1,x2,xx,nx,lon1c,lon2c,lwrite)

!	translate a window in longitude lon1...lon2 into grid points x1...x2
!	and return the exact cut-offs in lon1c,lon2c.  The axis is stored in xx

    implicit none
    integer :: x1,x2,nx
    real :: lon1,lon2,xx(nx),lon1c,lon2c,lonc
    logical :: lwrite
    integer :: xx1,xx2

!   first the trivial case

    if ( nx == 1 ) then
        if ( lwrite ) print * &
        ,'getlonwindow: X-axis only has one point'
        x1 = 1
        x2 = 1
        lon1c = xx(1)
        lon2c = xx(1)
        return
    end if

!   longitude is a mess - topologically not trivial

!   first reduce to 0..720
    lon1c = mod(lon1,360.)
    if ( lon1c < 0 ) lon1c = lon1c + 360
    lon2c = mod(lon2,360.)
    if ( lon2c < 0 ) lon2c = lon2c + 360
    if ( lon2c < lon1c ) lon2c = lon2c + 360
!   look out for whole-world averages that have been squashed to zero...
    if ( abs(lon2-lon1-360) < 0.001 ) lon2c = lon2c + 360 - 0.0001
!   next warn if the section comprises more than 180 degrees
    if ( .false. .and. lon2c-lon1c > 180 ) write(0,'(a)') 'Warning: more than 180 degrees longitude<br>'
    if ( lwrite ) then
        print *,'getlonwindow: lon1,lon2   = ',lon1,lon2
        print *,'              lon1c,lon2c = ',lon1c,lon2c
        print *,'              xx(1),xx(nx)= ',xx(1),xx(nx)
    endif

!   for non-wrapping grids, shift to maximize overlap

    if ( abs(xx(nx)-xx(1)+xx(2)-xx(1)-360) > 5 ) then
        if ( lon1c < xx(1) ) then
            if ( xx(nx) - (lon1c+360) > lon2c-xx(1) ) then
                lon1c = lon1c + 360
                lon2c = lon2c + 360
            endif
        endif
        if ( lon2c > xx(nx) ) then
            if ( (lon2c-360) - xx(1) > xx(nx)-lon1c ) then
                lon1c = lon1c - 360
                lon2c = lon2c - 360
            endif
        endif
        if ( lon1c < 2*xx(1)-xx(2) .or. lon2c > 2*xx(nx)-xx(nx-1) ) then
            write(0,*) 'getlonwindow: error: cannot extraplate to ', &
                lon1c,lon2c,' with an X-axis ',xx(1),xx(2),'...',xx(nx-1),xx(nx)
            lon1c = 3e33
            lon2c = 3e33
            return
        end if
        if ( lwrite ) then
            print *,'getlonwindow: input:     ',lon1,lon2
            print *,'              reduced to ',lon1c,lon2c
        endif
    endif

!	get grid points

    if ( xx(2) > xx(1) ) then
        do x1=1,nx
            if ( xx(x1) >= lon1c ) goto 110
        enddo
        do x1=nx+1,2*nx
            if ( xx(x1-nx)+360 >= lon1c ) goto 110
        enddo
        write(0,*) 'getlonwindow: error: cannot locate lon1c ',lon1c
        write(0,*) '                     in axis ',xx,' (also +360)'
        lon1c = 3e33
        lon2c = 3e33
        return
    110 continue
        do x2=1,nx
            if ( xx(x2) > lon2c ) goto 120
        enddo
        do x2=nx+1,2*nx
            if ( xx(x2-nx)+360 > lon2c ) goto 120
        enddo
        do x2=2*nx+1,3*nx
            if ( xx(x2-2*nx)+720 > lon2c ) goto 120
        enddo
        write(0,*) 'getlonwindow: error: cannot locate lon2c ',lon2c
        lon1c = 3e33
        lon2c = 3e33
        return
    120 continue
        x2 = x2 - 1
        if ( x2 < x1 ) then
!           the requested interval was between two grid points
            xx1 = 1 + mod(x1-1,nx)
            xx2 = 1 + mod(x2+nx-1,nx) ! x2 can be 0
            lonc = (lon1c + lon2c)/2
            if (  min(abs(xx(xx1)-lonc),abs(xx(xx1)+360-lon1c)) < &
                  min(abs(xx(xx2)-lonc),abs(xx(xx2)-360-lon1c),abs(xx(xx2)-720-lon1c)) ) then
!               round up
                x2 = x2 + 1
            else
!               round down
                x1 = x1 - 1
                if ( x1 <= 0 .or. x2 <= 0 ) then
                    x1 = x1 + nx
                    x2 = x2 + nx
                end if
            endif
        endif
!       hack
        if ( x1 <= nx ) then
            call getlon(x1,xx,nx,lon1c)
        else
            call getlon(x1-nx,xx,nx,lon1c)
        endif
        call getlon(x2+1,xx,nx,lon2c)
        if ( x2 > x1 .and. lon2c <= lon1c ) lon1c = lon1c - 360
        if ( lwrite ) then
            print *,'getlonwindow: output     ',x1,x2
            print *,'              corresponds',lon1c,lon2c
        endif
    else                    ! decreasing longitudes
        if ( lon2c > 720 ) then
            do x1=1-2*nx,-nx
                if ( xx(x1+2*nx)+720 <= lon2c ) goto 130
            enddo
        elseif ( lon2c > 360 ) then
            do x1=1-nx,0
                if ( xx(x1+nx)+360 <= lon2c ) goto 130
            enddo
        elseif ( lon2c > 0 ) then
            do x1=1,nx
                if ( xx(x1) <= lon2c ) goto 130
            enddo
        endif
        write(0,*) 'getlonwindow: error: cannot locate lon2c ',lon2c
        lon1c = 3e33
        lon2c = 3e33
        return
    130 continue
        if ( lon1c > 360 ) then
            do x2=1-nx,0
                if ( xx(x2+nx)+360 <= lon1c ) goto 140
            enddo
        elseif ( lon1 > 0 ) then
            do x2=1,nx
                if ( xx(x2) <= lon1c ) goto 140
            enddo
        endif
        write(0,*) 'getlonwindow: error: cannot locate lon1c ',lon1c
        write(0,*) '                     in axis ',xx
        lon1c = 3e33
        lon2c = 3e33
        return
    140 continue
        x2 = x2 - 1
        if ( x2 < x1 ) then
!           the requested interval was between two grid points
            lonc = (lon1c + lon2c)/2
            write(0,*) 'getlonwindow: untested case, please check'
            if (  min(abs(xx(x2)-lonc),abs(xx(x2)+360-lon1c)) < &
                  min(abs(xx(x1)-lonc),abs(xx(x1)-360-lon1c),abs(xx(x1)-720-lon1c)) ) then
!               round up
                x2 = x2 + 1
            else
!               round down
                x1 = x1 - 1
            endif
        endif
        call getlon(x1,xx,nx,lon2c)
        call getlon(x2+1,xx,nx,lon1c)
    endif
end subroutine getlonwindow

subroutine getlatwindow(lat1,lat2,y1,y2,yy,ny,lat1c,lat2c,lwrite)

!	translate a window in latitude lat1...lat2 into grid points y1...y2
!	and return the exact cut-offs in lat1c,lat2c.  The axis is stored in yy

    implicit none
    integer :: y1,y2,ny
    real :: lat1,lat2,yy(ny),lat1c,lat2c
    logical :: lwrite

!   special case...
    if ( ny == 1 ) then
        if ( lwrite ) print * &
        ,'getlatwindow: Y-axis only has one point'
        y1 = 1
        y2 = 1
        lat1c = yy(1)
        lat2c = yy(1)
        return
    endif

!   latitude: accept both directions
    if ( lat1 > lat2 ) then
        lat1c = lat2
        lat2c = lat1
    else
        lat1c = lat1
        lat2c = lat2
    endif
    if ( yy(2) > yy(1) ) then
        if ( lat1c < 2*yy(1)-yy(2) ) then
            write(0,*) 'getlatwindow: error: cannot extraplate to ', &
                lat1c,' with a Y-axis ',yy(1),yy(2),'...'
            lat1c = 3e33
            lat2c = 3e33
            return
        end if
        if ( lat2c > 2*yy(ny)-yy(ny-1) ) then
            write(0,*) 'getlatwindow: error: cannot extraplate to ', &
                lat2c,' with a Y-axis ',yy(ny),yy(ny-1),'...'
            lat1c = 3e33
            lat2c = 3e33
            return
        end if
        do y1=1,ny
            if ( yy(y1) >= lat1c ) goto 210
        enddo
        write(0,*) 'getlatwindow: error: cannot locate lat1c ',lat1c
        lat1c = 3e33
        lat2c = 3e33
        return
    210 continue
        do y2=1,ny
            if ( yy(y2) > lat2c ) goto 220
        enddo
    220 continue
        y2 = y2 - 1
        if ( y2 < y1 ) then
            if (  abs(yy(y1)-(lat1+lat2)/2) < &
            abs(yy(y2)-(lat1+lat2)/2) ) then
                y2 = y2 + 1
            else
                y1 = y1 - 1
            endif
        endif
        call getlon(y1,yy,ny,lat1c)
        call getlon(y2+1,yy,ny,lat2c)
    else
        if ( lat1c < 2*yy(ny)-yy(ny-1) ) then
            write(0,*) 'getlatwindow: error: cannot extraplate to ', &
                lat1c,' with a Y-axis ',yy(ny),yy(ny-1),'...'
            lat1c = 3e33
            lat2c = 3e33
            return
        end if
        if ( lat2c > 2*yy(1)-yy(2) ) then
            write(0,*) 'getlatwindow: error: cannot extraplate to ', &
                lat2c,' with a Y-axis ',yy(1),yy(2),'...'
            lat1c = 3e33
            lat2c = 3e33
            return
        end if
        do y1=1,ny
            if ( yy(y1) <= lat2c ) goto 230
        enddo
        write(0,*) 'getlatwindow: error: cannot locate lat2c ',lat2c
        lat1c = 3e33
        lat2c = 3e33
        return
    230 continue
        do y2=1,ny
            if ( yy(y2) < lat1c ) goto 240
        enddo
    240 continue
        y2 = y2 - 1
        if ( y2 < y1 ) then
            if (  abs(yy(y1)-(lat1+lat2)/2) < abs(yy(y2)-(lat1+lat2)/2) ) then
                y2 = y2 + 1
            else
                y1 = y1 - 1
            endif
        endif
        call getlon(y1,yy,ny,lat2c)
        call getlon(y2+1,yy,ny,lat1c)
    endif
    lat1c = max(-90.,lat1c)
    lat2c = min(+90.,lat2c)
    y1 = max(1,y1)
    y2 = min(ny,y2)
end subroutine getlatwindow

subroutine getlevwindow(lev1,lev2,zz,nz,z1,z2,lwrite)

!	translate a window in levels lev1...lev2 into grid points z1...z2

    implicit none
    integer :: z1,z2,nz
    real :: lev1,lev2,zz(nz)
    logical :: lwrite
    real :: lev1c,lev2c

!   special case...
    if ( nz == 1 ) then
        if ( lwrite ) print *,'getlevwindow: Z-axis only has one point'
        z1 = 1
        z2 = 1
        return
    endif

!   levels: accept both directions
    if ( lev1 > lev2 ) then
        lev1c = lev2
        lev2c = lev1
    else
        lev1c = lev1
        lev2c = lev2
    endif
    if ( zz(2) > zz(1) ) then
        do z1=1,nz
            if ( zz(z1) >= lev1c ) goto 210
        enddo
        write(0,*) 'getlevwindow: error: cannot locate lev1c ',lev1c
        lev1c = 3e33
        lev2c = 3e33
        return
        210 continue
        do z2=1,nz
            if ( zz(z2) > lev2c ) goto 220
        enddo
        220 continue
        z2 = z2 - 1
        if ( z2 < z1 ) then
            if (  abs(zz(z1)-(lev1+lev2)/2) < &
                  abs(zz(z2)-(lev1+lev2)/2) ) then
                z2 = z2 + 1
            else
                z1 = z1 - 1
            endif
        endif
        call getlon(z1,zz,nz,lev1c)
        call getlon(z2+1,zz,nz,lev2c)
    else
        do z1=1,nz
            if ( zz(z1) <= lev2c ) goto 230
        enddo
        write(0,*) 'getlevwindow: error: cannot locate lev2c ',lev2c
        lev1c = 3e33
        lev2c = 3e33
        return
    230 continue
        do z2=1,nz
            if ( zz(z2) < lev1c ) goto 240
        enddo
        240 continue
        z2 = z2 - 1
        if ( z2 < z1 ) then
            if (  abs(zz(z1)-(lev1+lev2)/2) < &
                  abs(zz(z2)-(lev1+lev2)/2) ) then
                z2 = z2 + 1
            else
                z1 = z1 - 1
            endif
        endif
        call getlon(z1,zz,nz,lev2c)
        call getlon(z2+1,zz,nz,lev1c)
    endif
    z1 = max(1,z1)
    z2 = min(nz,z2)
end subroutine getlevwindow

subroutine getlon(x1,xx,nx,lon1)

!   get the longitude corresponding to starting from x1
!   or stopping at x1-1

    implicit none
    integer :: x1,nx
    real :: xx(nx),lon1
    logical :: lwrite
    parameter (lwrite= .false. )

    if ( lwrite ) write(0,*) 'getlon: x1,nx = ',x1,nx
    if ( mod(x1,nx) == 1 ) then
        if ( x1 == 1 ) then
            lon1 = (3*xx(1) - xx(2))/2
            if ( lwrite) write(0,*) 'getlon: first point, lon1 = ',lon1
        else
            lon1 = (3*xx(nx) - xx(nx-1))/2
            if ( lwrite) write(0,*) 'getlon: last point, lon1 = ',lon1
        endif
    else
        lon1 = (xx(mod(x1-1,nx)) + xx(mod(x1-1,nx)+1))/2
        if ( lwrite) write(0,*) 'getlon: normal point, lon1 = ',lon1
    endif
end subroutine getlon

subroutine enscutoutwindow(x1,x2,y1,y2,xx,nx,xwrap,xrev,avex,yy &
    ,ny,avey,wx,wy,field,nxf,nyf,nens1,nens2,nperyear,firstyr &
    ,lastyr,yr1,yr2,lwrite)

!	cut out the region x1,x2,y1,y2 from a field and save into the
!   same field, adjusting the axes.

    implicit none
    integer :: x1,x2,y1,y2,nx,avex,ny,avey,nxf,nyf,nens1,nens2,nperyear &
        ,firstyr,lastyr,yr1,yr2
    real :: field(nxf,nyf,nperyear,firstyr:lastyr,0:nens2), &
        xx(nx),yy(ny),sum,wgt,fac
    logical :: xwrap,xrev,lwrite
    integer :: i,ii,i1,j,jj,j1,yr,mo,nnx,iens
    real :: wx(nx),wy(ny)

    if ( lwrite ) then
        print *,'enscutoutwindow: called with'
        print *,'x1,x2,y1,y2   = ',x1,x2,y1,y2
        print *,'nx,xwrap,xrev = ',nx,xwrap,xrev,xx(1),xx(2),'..',xx(nx-1),xx(nx)
        print *,'ny,yy         = ',ny,yy(1),yy(2),'..',yy(ny-1),yy(ny)
    endif

!   in order to avoid overlaps, move the data first in case the X
!   coordinate wraps over the selected area
    if ( xwrap .and. x2 > nx ) then
        if ( lwrite ) print *,'First shifting, overlap ',x1,x2
        if ( 1+x2-x1 > nx ) then
            write(0,*) 'cutoutwindow: internal error ',1+x2-x1,nx
            call exit(-1)
        endif
        nnx = 1+x2-x1
        if ( nnx < 1 ) nnx = nnx + nx
        if ( lwrite ) print *,'nnx = ',nnx
        do iens=nens1,nens2
            do yr=yr1,yr2
                do mo=1,nperyear
                    do j=1,1+y2-y1
                        do i=1,nnx
                            ii = 1 + mod(i+x1-2,nx)
                            wx(i) = field(ii,j+y1-1,mo,yr,iens)
                        enddo
                        do i=1,nnx
                            field(i,j,mo,yr,iens) = wx(i)
                            if ( abs(field(i,j,mo,yr,iens)) > 1e15 &
                                .and. field(i,j,mo,yr,iens) < 1e33 ) then
                                print *,'cutoutwindow: internal error:',i,j,mo,yr &
                                    ,field(i,j,mo,yr,iens)
                                call exit(-1)
                            endif
                        enddo
                    enddo
                enddo
            enddo
        enddo
        if ( lwrite ) then
            print *,'Shifted x grid from '
            print '(10f8.2)',(xx(i),i=1,nx)
        endif
        do i=1,nnx
            ii = 1 + mod(i+x1-2,nx)
            wx(i) = xx(ii)
            if ( lwrite ) print *,'shifting ',xx(ii),' from ',ii,' to ',i
        enddo
        do i=1,nnx
            xx(i) = wx(i)
            if ( i > 1 ) then
                if ( xx(i) > xx(i-1) .eqv. xrev ) then
                    if ( xrev ) then
                        xx(i) = xx(i) - 360
                    else
                        xx(i) = xx(i) + 360
                    endif
                endif
            endif
        enddo
        x1 = 1
        x2 = nnx
        if ( nnx < nx ) xwrap = .false. 
        nx = nnx
        if ( y1 > 1 ) then
            do j=1,y2-y1+1
                yy(j) = yy(j+y1-1)
            enddo
            ny = y2-y1+1
            y1 = 1
            y2 = ny
        endif
        if ( lwrite ) then
            print *,'into '
            print '(10f8.2)',(xx(i),i=1,nx)
        endif
    endif

!   get weights along old axes
    call getweights('x',xx,wx,nx,xwrap,lwrite)
    call getweights('y',yy,wy,ny,.false.,lwrite)

!   real work
    write(0,*) '@@@ x1,y1 = ',x1,y1
    if ( y1 > 1 .or. x1 > 1 .or. avex > 1 .or. avey > 1 ) then
        if ( avex > 1 .or. avey > 1 ) then
            if ( lwrite ) print *,'Averaging fields'
!           area-weighted average of field
            do iens=nens1,nens2
                do yr=yr1,yr2
                    do mo=1,nperyear
                        do j=0,(y2-y1)/avey
                            do i=0,(x2-x1)/avex
                                sum = 0
                                wgt = 0
                                fac = 0
                                do j1=0,min(avey-1,y2-y1-avey*j)
                                    do i1=0,min(avex-1,x2-x1-avex*i)
                                        ii = x1+avex*i+i1
                                        jj = y1+avey*j+j1
                                        if ( field(ii,jj,mo,yr,iens) < 1e33) then
                                            sum = sum +wx(ii)*wy(jj)*field(ii,jj,mo,yr,iens)
                                            wgt = wgt +wx(ii)*wy(jj)
                                        endif
                                        fac = fac + wx(ii)*wy(jj)
                                    enddo
                                enddo
                                if ( wgt >= fac/2 .and. wgt > 0 ) then
                                    field(1+i,1+j,mo,yr,iens) = sum/wgt
                                else
                                    field(1+i,1+j,mo,yr,iens) = 3e33
                                endif
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        else                ! avex>1 || avey>1
            do iens=nens1,nens2
                do yr=yr1,yr2
                    do mo=1,nperyear
                        do j=0,y2-y1
                            do i=0,x2-x1
                                field(1+i,1+j,mo,yr,iens) = field(x1+i,y1+j,mo,yr,iens)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        endif               ! avex>1 || avey>1
!       adjust grids as well
        if ( x2-x1+1 < nx ) then
            xwrap = .false. 
        endif
        if ( avex > 1 ) then
            if ( lwrite ) then
                print *,'transformed x grid from '
                print '(10f8.2)',(xx(i),i=1,nx)
                print '(10f8.2)',(wx(i),i=1,nx)
            endif
            do i=0,(x2-x1)/avex
                sum = 0
                wgt = 0
                do i1=0,min(avex-1,x2-x1-avex*i)
                    ii = 1+mod(x1+avex*i+i1-1,nx)
                    wgt = wgt + wx(ii)
                    if ( i1 > 0 .and. ( &
                    xx(ii) > xx(1+mod(x1+avex*i-1,nx)) .eqv. xrev ) ) then
                        if ( xrev ) then
                            sum = sum + (xx(ii)-360)*wx(ii)
                        else
                            sum = sum + (xx(ii)+360)*wx(ii)
                        endif
                    else
                        sum = sum + xx(ii)*wx(ii)
                    endif
                enddo
                if ( wgt == 0 ) then
                    write(0,*) 'cutoutwindow: error: weight=0'
                    call exit(-1)
                endif
                xx(i+1) = sum/wgt
                wx(i+1) = wgt
                if ( i > 0 ) then
                    if ( xx(i+1) > xx(i) .eqv. xrev ) then
                        if ( xrev ) then
                            xx(i+1) = xx(i+1) - 360
                        else
                            xx(i+1) = xx(i+1) + 360
                        endif
                    endif
                endif
            enddo
            nx = 1 + (x2-x1)/avex
            if ( lwrite ) then
                print *,'into '
                print '(10f8.2)',(xx(i),i=1,nx)
                print '(10f8.2)',(wx(i),i=1,nx)
            endif
        elseif ( x1 > 1 ) then
            if ( lwrite ) then
                print *,'shifted x grid from '
                print '(10f8.2)',(xx(i),i=1,nx)
                print '(10f8.2)',(wx(i),i=1,nx)
            endif
            do i=0,x2-x1
                xx(i+1) = xx(i+x1)
            enddo
            do i=0,x2-x1
                wx(i+1) = wx(i+x1)
            enddo
            nx = 1 + x2 - x1
            if ( lwrite ) then
                print *,'to'
                print '(10f8.2)',(xx(i),i=1,nx)
                print '(10f8.2)',(wx(i),i=1,nx)
            endif
        endif
        if ( avey > 1 ) then
            if ( lwrite ) then
                print *,'transformed y grid from '
                print '(10f8.2)',(yy(i),i=1,ny)
                print '(10f8.2)',(wy(i),i=1,ny)
            endif
            do j=0,(y2-y1)/avey
                sum = 0
                wgt = 0
                do j1=0,min(avey-1,y2-y1-avey*j)
                    jj = y1+avey*j+j1
                    wgt = wgt + wy(jj)
                    sum = sum + wy(jj)*yy(jj)
                enddo
                yy(j+1) = sum/wgt
                wy(j+1) = wgt
            enddo
            ny = 1 + (y2-y1)/avey
            if ( lwrite ) then
                print *,'into'
                print '(10f8.2)',(yy(i),i=1,ny)
                print '(10f8.2)',(wy(i),i=1,ny)
            endif
        elseif ( y1 > 1 ) then
            if ( lwrite ) then
                print *,'shifted y grid from '
                print '(10f8.2)',(yy(i),i=1,ny)
                print '(10f8.2)',(wy(i),i=1,ny)
            endif
            do j=0,y2-y1
                yy(j+1) = yy(j+y1)
            enddo
            do j=0,y2-y1
                wy(j+1) = wy(j+y1)
            enddo
            ny = 1 + y2 - y1
            if ( lwrite ) then
                print *,'to'
                print '(10f8.2)',(yy(i),i=1,ny)
                print '(10f8.2)',(wy(i),i=1,ny)
            endif
        else
            if ( lwrite ) print *,'just truncated the y grid'
            ny = y2
        endif
    else                    ! adjust upper boundaries
        if ( lwrite ) print *,'nothing to do'
        nx = x2
        ny = y2
    endif                   ! anything to do at all?
end subroutine enscutoutwindow

subroutine biggerwindow(xx,nx,yy,ny,xwrap, &
    lon1,lon2,lat1,lat2,x1,x2,y1,y2)

!   increase the area x1,x2,y1,y2 so that lon1,lon2,lat1,lat2 is
!   inside the grid points (and not closest to the halfway points)

    implicit none
    integer :: nx,ny,x1,x2,y1,y2
    real :: xx(nx),yy(ny),lon1,lon2,lat1,lat2
    logical :: xwrap
    integer :: n,xx1,xx2
    real :: l1,l2
    logical :: lwrite
    lwrite = .false. 

!   latitude is easy

    if ( lwrite ) then
        write(0,*) 'biggerwindow: requested: ',lat1,lat2,'<br>'
        write(0,*) 'biggerwindow: original:  ',y1,y2,yy(y1),yy(y2),'<br>'
    endif
    if ( yy(1) < yy(2) ) then
        if ( yy(y1) > min(lat1,lat2) .and. y1 > 1 )  y1 = y1 - 1
        if ( yy(y2) < max(lat1,lat2) .and. y2 < ny ) y2 = y2 + 1
    else
        if ( yy(y1) < max(lat1,lat2) .and. y1 > 1 )  y1 = y1 - 1
        if ( yy(y2) > min(lat1,lat2) .and. y2 < ny ) y2 = y2 + 1
    endif
    if ( lwrite ) then
        write(0,*) 'biggerwindow: adjusted:  ',y1,y2,yy(y1),yy(y2),'<br>'
    endif

!   longitude is not easy

    if ( lwrite ) then
        write(0,*) 'biggerwindow: requested: ',lon1,lon2,'<br>'
        write(0,*) 'biggerwindow: original:  ',x1,x2, &
        xx(1+mod(x1-1,nx)),xx(1+mod(x2-1,nx)),xwrap,'<br>'
    endif
    n = x2 - x1 + 1
    if ( xwrap .and. n <= 0 ) n = n + nx
    if ( n /= nx ) then
        xx1 = x1
        if ( xwrap ) then
            if ( xx1 < 1 )  xx1 = xx1 + nx
            if ( xx1 > nx ) xx1 = xx1 - nx
        endif
        xx2 = x2
        if ( xwrap ) then
            if ( xx2 < 1 )  xx2 = xx2 + nx
            if ( xx2 > nx ) xx2 = xx2 - nx
        endif
        if ( abs(xx(xx1)-lon1+360) < &
        min(abs(xx(xx1)-lon1),abs(xx(xx1)-lon1-360)) ) then
            l1 = lon1 - 360
        elseif ( abs(xx(xx1)-lon1-360) < &
            min(abs(xx(xx1)-lon1),abs(xx(xx1)-lon1+360)) ) then
            l1 = lon1 + 360
        else
            l1 = lon1
        endif
        if ( abs(xx(xx2)-lon2+360) < &
        min(abs(xx(xx2)-lon2),abs(xx(xx2)-lon2-360)) ) then
            l2 = lon2 - 360
        elseif ( abs(xx(xx2)-lon2-360) < &
            min(abs(xx(xx2)-lon2),abs(xx(xx2)-lon2+360)) ) then
            l2 = lon2 + 360
        else
            l2 = lon2
        endif
        if ( xwrap ) then
            if ( xx(xx1) > l1 ) x1 = x1 - 1
            if ( xx(xx2) < l2 ) x2 = x2 + 1
        else
            if ( xx(xx1) > l1 .and. xx1 > 1 )  x1 = x1 - 1
            if ( xx(xx2) < l2 .and. xx2 < nx ) x2 = x2 + 1
        endif
    endif
    if ( lwrite ) then
        write(0,*) 'biggerwindow: adjusted:  ',x1,x2, &
        xx(1+mod(x1-1,nx)),xx(1+mod(x2-1,nx)),'<br>'
    endif
end subroutine biggerwindow
