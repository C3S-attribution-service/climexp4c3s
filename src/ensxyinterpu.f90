subroutine xyinterpu(field1,xx1,nx1,yy1,ny1,field2,xx2,nx2,yy2 &
    ,ny2,xx,nx,yy,ny,yr11,yr21,yr12,yr22,nxf,nyf,nzf,nz &
    ,nperyear,intertype,lwrite)
    implicit none

!   single file version of xyinterou belopw

    include "params.h"
    integer :: nx1,ny1,nx2,ny2,nx,ny,yr11,yr21,yr12,yr22,nxf,nyf,nzf,nz,nperyear,intertype
    real :: xx(nxmax),yy(nymax), &
        field1(nxf,nyf,nzf,nperyear,yr11:yr21),xx1(nx1),yy1(ny1), &
        field2(nxf,nyf,nzf,nperyear,yr12:yr22),xx2(nx2),yy2(ny2)
    logical :: lwrite
    call ensxyinterpu(field1,xx1,nx1,yy1,ny1,0,0 &
        ,field2,xx2,nx2,yy2,ny2,0,0,xx,nx,yy,ny,yr11,yr21 &
        ,yr12,yr22,nxf,nyf,nzf,nz,nperyear,intertype,lwrite)
end subroutine xyinterpu

subroutine ensxyinterpu(field1,xx1,nx1,yy1,ny1,mens11,mens1 &
    ,field2,xx2,nx2,yy2,ny2,mens12,mens2,xx,nx,yy,ny,yr11,yr21 &
    ,yr12,yr22,nxf,nyf,nzf,nz,nperyear,intertype,lwrite)

!   interpolate and truncate the fields so that they are
!   on a common grid.  The fields are modified in situ
!   to save on memory.  Undefined is assumed 3e33.
!   intertype:
!   0: highest resolution grid
!      -1: lowest resolution grid
!   1: grid 1
!   2: grid 2

    implicit none
    include 'params.h'
    integer :: nlevmax
    parameter(nlevmax=2)
    integer :: nx1,ny1,nx2,ny2,nx,ny,yr11,yr21,yr12,yr22,nxf,nyf,nzf,nz &
        ,nperyear,intertype,mens11,mens1,mens12,mens2
    real :: xx(nxmax),yy(nymax), &
        field1(nxf,nyf,nzf,nperyear,yr11:yr21,0:mens1),xx1(nx1),yy1(ny1), &
        field2(nxf,nyf,nzf,nperyear,yr12:yr22,0:mens2),xx2(nx2),yy2(ny2)
    logical :: lwrite
    integer :: i,j,i1x1,i2x1,i1x2,i2x2,i1y1,i2y1,i1y2,i2y2, &
        ix(nxmax),iy(nymax),yr,mo,jz,iens
    real :: ax(nxmax),ay(nymax),f(nxmax,nymax)
    logical :: x1rev,x2rev,y1rev,y2rev,lintx1,lintx2,linty1,linty2

!   debug output

    if ( lwrite ) then
        print *,'ensxyinterpu: input'
        print *,'field1 ',nx1,ny1,nz,mens11,mens1,yr11,yr21
        print *,'field2 ',nx2,ny2,nz,mens12,mens2,yr12,yr22
    endif

!   sanity check (<= I was insane once)

    if ( nx1 > nxf .or. nx2 > nxf .or. &
    ny1 > nyf .or. ny2 > nyf ) then
        write(0,*) 'ensxyinterpu: error: grid larger than physical', &
            ' dimensions: ',1,nx1,ny1,2,nx2,ny2,'f',nxf,nyf
        write(*,*) 'ensxyinterpu: error: grid larger than physical', &
            ' dimensions: ',1,nx1,ny1,2,nx2,ny2,'f',nxf,nyf
        call exit(-1)
    endif

!   all the grid-based work

    call gridprops(xx1,nx1,yy1,ny1,xx2,nx2,yy2,ny2,xx,nx,yy,ny,intertype,lwrite, &
        lintx1,lintx2,linty1,linty2,x1rev,x2rev,y1rev,y2rev,ax,ix,nxmax,ay,iy,nymax, &
        i1x1,i2x1,i1x2,i2x2,i1y1,i2y1,i1y2,i2y2)

!   anything to do?

    if ( .not. lintx1 .and. .not. lintx2 .and. &
         .not. linty1 .and. .not. linty2 .and. &
        (x1rev.eqv.x2rev) .and. (y1rev.eqv.y2rev) ) then
        if ( lwrite ) print *,'no interpolation necessary'
!       reverse the axes back to what they were
        if ( x1rev ) then
            call revit(xx1,nx1)
            call revit(xx2,nx2)
            call revit(xx,nx)
        endif
        if ( y1rev ) then
            call revit(yy1,ny1)
            call revit(yy2,ny2)
            call revit(yy,ny)
        endif
        return
    endif

!   perform the interpolations

    do iens=mens11,mens1
        if ( lwrite ) then
            if ( lintx1 ) then
                print *,'interpolating f1 x ',iens,yr11,yr21,nperyear,nz
            elseif ( linty1 .or. i1x1 /= 1 .or. i1y1 /= 1 .or. x1rev .or. y1rev ) then
                print *,'shifting f1 x ',iens,yr11,yr21,nperyear,nz
            endif
            if ( linty1 ) then
                print *,'interpolating f1 y ',iens,yr11,yr21,nperyear,nz
            elseif ( lintx1 .or. i1x1 /= 1 .or. i1y1 /= 1 .or. x1rev .or. y1rev ) then
                print *,'shifting f1 y ',iens,yr11,yr21,nperyear,nz
            endif
        endif
        do yr=yr11,yr21
            do mo=1,nperyear
                do jz=1,nz
!                   interpolate/copy from field1 to f
                    if ( lintx1 ) then
                        call doxint(f,nxmax,nymax,field1(1,1,jz,mo &
                            ,yr,iens),nxf,nyf,nx1,ny1,ax,ix,nx,x1rev)
                    elseif ( linty1 .or. i1x1 /= 1 .or. i1y1 /= 1 &
                         .or. x1rev .or. y1rev ) then
                        call doxnint(f,nxmax,nymax,field1(1,1,jz,mo &
                            ,yr,iens),nxf,nyf,nx1,ny1,i1x1,nx,x1rev)
                    endif
!                   interpolate/copy from f back to field1
                    if ( linty1 ) then
                        call doyint(field1(1,1,jz,mo,yr,iens),nxf &
                            ,nyf,f,nxmax,nymax,nx1,ny1,ay,iy,ny,nx,y1rev)
                    elseif ( lintx1 .or. i1x1 /= 1 .or. i1y1 /= 1 .or. x1rev .or. y1rev ) then
                        call doynint(field1(1,1,jz,mo,yr,iens),nxf &
                            ,nyf,f,nxmax,nymax,nx1,ny1,i1y1,nx,ny,y1rev)
                    endif
                enddo       ! jz
            enddo           ! mo
        enddo               ! yr1
    enddo                   ! iens
    do iens=mens12,mens2
        if ( lwrite ) then
            if ( lintx2 ) then
                print *,'interpolating f2 x ',iens,yr12,yr22,nperyear,nz
            elseif ( linty2 .or. i1x2 /= 1 .or. i1y2 /= 1 .or. x2rev .or. y2rev ) then
                print *,'shifting      f2 x ',iens,yr12,yr22,nperyear,nz
            endif
            if ( linty2 ) then
                print *,'interpolating f2 y ',iens,yr12,yr22,nperyear,nz
            elseif ( lintx2 .or. i1x2 /= 1 .or. i1y2 /= 1 .or. x2rev .or. y2rev ) then
                print *,'shifting      f2 y ',iens,yr12,yr22,nperyear,nz
            endif
        endif
        do yr=yr12,yr22
            do mo=1,nperyear
                do jz=1,nz
!                   interpolate/copy from field2 to f
                    if ( lintx2 ) then
                        call doxint(f,nxmax,nymax,field2(1,1,jz,mo &
                            ,yr,iens),nxf,nyf,nx2,ny2,ax,ix,nx,x2rev)
                    elseif ( linty2 .or. i1x2 /= 1 .or. i1y2 /= 1 .or. x2rev .or. y2rev ) then
                        call doxnint(f,nxmax,nymax,field2(1,1,jz,mo &
                            ,yr,iens),nxf,nyf,nx2,ny2,i1x2,nx,x2rev)
                    endif
!                   interpolate/copy from f back to field2
                    if ( linty2 ) then
                        call doyint(field2(1,1,jz,mo,yr,iens),nxf &
                            ,nyf,f,nxmax,nymax,nx2,ny2,ay,iy,ny,nx,y2rev)
                    elseif ( lintx2 .or. i1x2 /= 1 .or. i1y2 /= 1 .or. x2rev .or. y2rev ) then
                        call doynint(field2(1,1,jz,mo,yr,iens),nxf &
                            ,nyf,f,nxmax,nymax,nx2,ny2,i1y2,nx,ny,y2rev)
                    endif
                enddo       ! jz
            enddo           ! mo
        enddo               ! yr2
    enddo                   ! iens

!   that's it

end subroutine ensxyinterpu

subroutine gridprops(xx1,nx1,yy1,ny1,xx2,nx2,yy2,ny2, &
    xx,nx,yy,ny,intertype,lwrite, &
    lintx1,lintx2,linty1,linty2,x1rev,x2rev,y1rev,y2rev, &
    ax,ix,nxmax,ay,iy,nymax, &
    i1x1,i2x1,i1x2,i2x2,i1y1,i2y1,i1y2,i2y2)

!   establish some properties of the grids

    implicit none
    integer :: nx1,ny1,nx2,ny2,nx,ny,intertype,nxmax,nymax
    real :: xx1(nx1),yy1(ny1),xx2(nx2),yy2(ny2),xx(nxmax),yy(nxmax)
    logical :: lwrite
    logical :: lintx1,lintx2,linty1,linty2,x1rev,x2rev,y1rev,y2rev
    integer :: ix(nxmax),iy(nymax), &
    i1x1,i2x1,i1x2,i2x2,i1y1,i2y1,i1y2,i2y2
    real :: ax(nxmax),ay(nymax)
    integer :: i
    real :: lat1,lat2,lon1,lon2,dxmin1,dxmin2,dymin1,dymin2,d1,d2
    logical :: x1wrap,x2wrap

!   reversed grids?

    if ( nx1 > 1 ) then
        x1rev = xx1(2) < xx1(1)
    else
        x1rev = .false. 
    endif
    if ( nx2 > 1 ) then
        x2rev = xx2(2) < xx2(1)
    else
        x2rev = .false. 
    endif
    if ( ny1 > 1 ) then
        y1rev = yy1(2) < yy1(1)
    else
        y1rev = .false. 
    endif
    if ( ny2 > 1 ) then
        y2rev = yy2(2) < yy2(1)
    else
        y2rev = .false. 
    endif
!   save myself a lot of headaches later on
    if ( x1rev ) call revit(xx1,nx1)
    if ( x2rev ) call revit(xx2,nx2)
    if ( y1rev ) call revit(yy1,ny1)
    if ( y2rev ) call revit(yy2,ny2)
!   grid wraps around the earth?  (I assume it is in degrees)
    if ( nx1 > 1 ) then
        x1wrap = abs(2*xx1(nx1)-xx1(nx1-1)-360-xx1(1)) < 1e-3
    else
        x1wrap = .false. 
    end if
    if ( nx2 > 1 ) then
        x2wrap = abs(2*xx2(nx2)-xx2(nx2-1)-360-xx2(1)) < 1e-3
    else
        x2wrap = .false. 
    end if
    if ( lwrite ) then
        print *,'interpu grid properties'
        if ( x1rev ) print *,'grid1 x reversed'
        if ( x2rev ) print *,'grid2 x reversed'
        if ( y1rev ) print *,'grid1 y reversed'
        if ( y2rev ) print *,'grid2 y reversed'
        if ( x1wrap ) print *,'grid1 x wraps'
        if ( x2wrap ) print *,'grid2 x wraps'
    endif

!   cut-offs: longitude

    if ( x1wrap ) then
        i1x2 = 1
        i2x2 = nx2
        if ( x2wrap ) then
!           both wrap - no cutout region
            i1x1 = 1
            i2x1 = nx1
        else                ! x2wrap
!           cut out the region corresponding with field2
            call getlonwindow(xx2(1),xx2(nx2),i1x1,i2x1,xx1,nx1,d1,d2,lwrite)
            if ( d1 > 1e33 .or. d2 > 1e33 ) then
                write(0,*) 'ensxyinterpu: cannot find lon ',lon1,lon2
                call exit(-1)
            endif
            if ( i1x1 > 1   .and. abs(xx1(i1x1)-xx2(1)) > 1e-3 ) i1x1 = i1x1-1
            if ( i2x1 < nx1 .and. abs(xx1(1+mod(i2x1-1,nx1))-xx2(nx2)) > 1e-3 ) i2x1 = i2x1+1
        endif
    else                    ! x1wrap
        if ( x2wrap ) then
!           cut out the region corresponding with field1
            i1x1 = 1
            i2x1 = nx1
            lon1 = xx1(1)
            lon2 = xx1(nx1)
            call getlonwindow(xx1(1),xx1(nx1),i1x2,i2x2,xx2,nx2,d1,d2,lwrite)
            if ( d1 > 1e33 .or. d2 > 1e33 ) then
                write(0,*) 'ensxyinterpu: cannot find lon ',lon1,lon2
                call exit(-1)
            endif
            if ( i1x2 > 1   .and. abs(xx2(i1x2)-xx1(1)) > 1e-3 ) i1x2 = i1x2-1
            if ( i2x2 < nx2 .and. abs(xx2(i2x1)-xx1(nx1)) > 1e-3) i2x2 = i2x2+1
        else                ! x2wrap
!           cut out the region corresponding with the intersection of field1,field2
            lon1 = max(xx1(1),xx2(1))
            lon2 = min(xx1(nx1),xx2(nx2))
            if ( lon1 > lon2 ) then
                if ( xx1(1) > xx2(nx2) ) then
                    xx1(1:nx1) = xx1(1:nx1) - 360
                elseif ( xx2(1) < xx1(nx1) ) then
                    xx2(1:nx2) = xx2(1:nx2) - 360
                endif
                lon1 = max(xx1(1),xx2(1))
                lon2 = min(xx1(nx1),xx2(nx2))
                if ( lon1 > lon2 ) then
                    write(0,*) 'ensxyinterpu: error: no overlap in longitude' &
                        ,xx1(1),xx1(nx1),xx2(1),xx2(nx2)
                    call exit(-1)
                endif
            endif
            call getlonwindow(lon1,lon2,i1x1,i2x1,xx1,nx1,d1,d2,lwrite)
            if ( d1 > 1e33 .or. d2 > 1e33 ) then
                write(0,*) 'ensxyinterpu: cannot find lon ',lon1,lon2
                call exit(-1)
            endif
            if ( i1x1 > 1   .and. abs(xx1(i1x1)-lon1) > 1e-3 ) i1x1 = i1x1-1
            if ( i2x1 < nx1 .and. abs(xx1(i2x1)-lon2) > 1e-3 ) i2x1 = i2x1+1
            call getlonwindow(lon1,lon2,i1x2,i2x2,xx2,nx2,d1,d2,lwrite)
            if ( d1 > 1e33 .or. d2 > 1e33 ) then
                write(0,*) 'ensxyinterpu: cannot find lon ',lon1,lon2
                call exit(-1)
            endif
            if ( i1x2 > 1   .and. abs(xx2(i1x2)-lon1) > 1e-3 ) i1x2 = i1x2-1
            if ( i2x2 < nx2 .and. abs(xx2(i2x1)-lon2) > 1e-3 ) i2x2 = i2x2+1
        endif
    endif
    if ( lwrite ) then
        print '(a)',' ensxyinterpu: cutout region x'
        print '(a,2i7,a,2i7)',' grid1: cut out ',i1x1,i2x1,' of ',1,nx1
        print '(a,2f7.2,a,2f7.2)',' grid1: cut out ', &
            xx1(i1x1),xx1(1+mod(i2x1-1,nx1)),' of ',xx1(1),xx1(nx1)
        print '(a,2i7,a,2i7)',' grid2: cut out ',i1x2,i2x2,' of ',1,nx2
        print '(a,2f7.2,a,2f7.2)',' grid2: cut out ', &
            xx2(i1x2),xx2(1+mod(i2x2-1,nx2)),' of ',xx2(1),xx2(nx2)
    endif

!   cut-offs latitude: simpler (no wrapping)

    lat1 = max(yy1(1),yy2(1))
    lat2 = min(yy1(ny1),yy2(ny2))
    if ( lat1 > lat2 ) then
        write(0,*) 'ensxyinterpu: error: no overlap in latitude' &
            ,yy1(1),yy1(ny1),yy2(1),yy2(ny2)
        call exit(-1)
    endif
    call getlatwindow(lat1,lat2,i1y1,i2y1,yy1,ny1,d1,d2,lwrite)
    if ( i1y1 > 1   .and. abs(yy1(i1y1)-lat1) > 1e-3 ) i1y1 = i1y1-1
    if ( i2y1 < ny1 .and. abs(yy1(i2y1)-lat2) > 1e-3 ) i2y1 = i2y1+1
    call getlatwindow(lat1,lat2,i1y2,i2y2,yy2,ny2,d1,d2,lwrite)
    if ( i1y2 > 1   .and. abs(yy2(i1y2)-lat1) > 1e-3 ) i1y2 = i1y2-1
    if ( i2y2 < ny2 .and. abs(yy2(i2y2)-lat2) > 1e-3 ) i2y2 = i2y2+1
    if ( lwrite ) then
        print '(a)','ensxyinterpu: cutout region y'
        print '(a,2i7,a,2i7)',' grid1: cut out ',i1y1,i2y1,' of ',1,ny1
        print '(a,2f7.2,a,2f7.2)',' grid1: cut out ', &
            yy1(i1y1),yy1(i2y1),' of ',yy1(1),yy1(ny1)
        print '(a,2i7,a,2i7)',' grid2: cut out ',i1y2,i2y2,' of ',1,ny2
        print '(a,2f7.2,a,2f7.2)',' grid2: cut out ', &
            yy2(i1y2),yy2(i2y2),' of ',yy2(1),yy2(ny2)
    endif

!   minimum grid distances

    if ( x1wrap .and. i1x1 == 1 .and. i2x1 == nx1 ) then
        dxmin1 = xx1(1) - xx1(nx1) + 360
    else
        dxmin1 = 360
    endif
    do i=i1x1,i2x1-1
        dxmin1 = min(dxmin1, &
        abs(xx1(1+mod(i,nx1))-xx1(1+mod(i-1,nx1))), &
        abs(xx1(1+mod(i,nx1))-xx1(1+mod(i-1,nx1))+360), &
        abs(xx1(1+mod(i,nx1))-xx1(1+mod(i-1,nx1))-360))
    enddo
    if ( x2wrap .and. i1x2 == 1 .and. i2x2 == nx2 ) then
        dxmin2 = xx2(1) - xx2(nx2) + 360
    else
        dxmin2 = 360
    endif
    do i=i1x2,i2x2-1
        dxmin2 = min(dxmin2, &
        abs(xx2(1+mod(i,nx2))-xx2(1+mod(i-1,nx2))), &
        abs(xx2(1+mod(i,nx2))-xx2(1+mod(i-1,nx2))+360), &
        abs(xx2(1+mod(i,nx2))-xx2(1+mod(i-1,nx2))-360))
    enddo
    dymin1 = 180
    do i=i1y1,i2y1-1
        dymin1 = min(dymin1,yy1(i+1)-yy1(i))
    enddo
    dymin2 = 180
    do i=i1y2,i2y2-1
        dymin2 = min(dymin2,yy2(i+1)-yy2(i))
    enddo
    call chooseint(lintx1,lintx2,dxmin1,dxmin2,i1x1,i2x1,i1x2,i2x2 &
        ,xx1,xx2,intertype,lwrite)
    call chooseint(linty1,linty2,dymin1,dymin2,i1y1,i2y1,i1y2,i2y2 &
        ,yy1,yy2,intertype,lwrite)
    if ( lwrite ) then
        print *,'grid1 x min ',dxmin1
        print *,'grid2 x min ',dxmin2
        print *,'grid1 y min ',dymin1
        print *,'grid2 y min ',dymin2
        if ( lintx1 ) print *,'interpolating field1 x'
        if ( lintx2 ) print *,'interpolating field2 x'
        if ( linty1 ) print *,'interpolating field1 y'
        if ( linty2 ) print *,'interpolating field2 y'
    endif

!   get new grid, compute interpolation coefficients
    if ( lintx1 ) then
        if ( lintx2 ) then
            write(0,*) 'ensxyinterpu: cannot interpolate both x fields'
            call exit(-1)
        else
            call setupint(xx,ax,ix,nx,xx1,nx1,x1wrap,xx2,nx2,i1x2 &
                ,i2x2,x2wrap,lwrite)
        endif
    else
        if ( lintx2 ) then
            call setupint(xx,ax,ix,nx,xx2,nx2,x2wrap,xx1,nx1,i1x1 &
                ,i2x1,x1wrap,lwrite)
        else
!           both grids are already equal, copy grid info
            call setupint(xx,ax,ix,nx,xx2,nx2,x2wrap,xx1,nx1,i1x1 &
                ,i2x1,x1wrap,lwrite)
        endif
    endif
    if ( linty1 ) then
        if ( linty2 ) then
            write(0,*) 'ensxyinterpu: cannot interpolate both y fields'
            call exit(-1)
        else
            call setupint(yy,ay,iy,ny,yy1,ny1, .false. ,yy2,ny2,i1y2 &
                ,i2y2,.false.,lwrite)
        endif
    else
        if ( linty2 ) then
            call setupint(yy,ay,iy,ny,yy2,ny2, .false. ,yy1,ny1,i1y1 &
                ,i2y1,.false.,lwrite)
        else
!           both grids are already equal, copy grid info
            call setupint(yy,ay,iy,ny,yy2,ny2, .false. ,yy1,ny1,i1y1 &
                ,i2y1,.false.,lwrite)
        endif
    endif

end subroutine gridprops
