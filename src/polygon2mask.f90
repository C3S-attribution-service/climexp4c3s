program polygon2mask
    !
    ! converts a polygon to a mask, given a grid
    !
    use lsdata
    implicit none
    include 'params.h'
    include 'netcdf.inc'
    integer,parameter :: npolmax=500000,nvarmax=10
    integer i,j,npol,mens,mens1,ncid,nx,ny,nz,nt,nperyear,firstyr,firstmo,endian,nvars,jvars(6,nvarmax),ntvarid
    integer nens1,nens2,ntmax,itimeaxis(1),iarg,ncoords
    real xx(nxmax),yy(nymax),zz(nzmax),undef,xxls(nxmax),yyls(nymax)
    double precision polygon(2,npolmax),coords(3,npolmax)
    real mask(nxmax,nymax)
    character infile*255,datfile*255,maskfile*255,lz(3)*20,ltime*120,title*1047,history*4095, &
        & string*80,lsmasktype*4
    character vars(nvarmax)*50,lvars(nvarmax)*100,svars(nvarmax)*100,units(nvarmax)*100 &
        & ,cell_methods(nvarmax)*100,metadata(2,100)*2000
    character pole*2,clwrite*1
    logical lwrite,xrev,yrev,xwrap

    lwrite = .false.
    lsmasktype = 'all '
    if ( command_argument_count() < 3 ) then
        write(0,*) 'usage: polygon2mask grid.nc polygon.txt [sp] '// &
            & '[lsmask file.nc all|land|sea|notland|notsea|all] [debug] mask.nc'
        write(0,*) 'with grid.nc a netcdf file that defines the grid'
        write(0,*) 'with polygon.txt a text file with xi yi on each ine'
        write(0,*) 'sp denotes that the South Pole is inside the polygon'
        write(0,*) 'mask.nc is 1 when the grid point is inside the polygon, '
        write(0,*) '0 when outside (to be exact odd/even no of crossings).'
        stop
    end if
    do i=3,command_argument_count()-1
        call get_command_argument(i,string)
        call tolower(string)
        if ( string == 'lwrite' .or. string == 'debug') then
            lwrite = .true.
        end if
    end do
    !
    ! read grid xx,yy from netcdf (or grads) file
    !
    call get_command_argument(1,infile)
    call getmetadata(infile,mens1,mens,ncid,datfile,nxmax,nx &
        ,xx,nymax,ny,yy,nzmax,nz,zz,lz,nt,nperyear,firstyr,firstmo &
        ,ltime,undef,endian,title,history,nvarmax,nvars,vars,jvars &
        ,lvars,svars,units,cell_methods,metadata,lwrite)
    i = nf_close(ncid)
    if ( ny <= 1 ) then
        write(0,*) 'polygon2mask: error: expecting ta least two grid points in latitude'
        call exit(-1)
    end if
    if ( nx <= 1 ) then
        write(0,*) 'polygon2mask: error: expecting ta least two grid points in longitude'
        call exit(-1)
    end if
    call getxyprop(xx,nx,yy,ny,xrev,yrev,xwrap)
    !
    ! read polygon from text file
    !
    call get_command_argument(2,datfile)
    call read_polygon(datfile,npol,npolmax,polygon,lwrite)
    !
    if ( command_argument_count() > 3 ) then
        call get_command_argument(3,pole)
        if ( pole /= 'sp' ) pole='np'
        do iarg=3,command_argument_count()-1
            call get_command_argument(iarg,string)
            if ( string(1:7) == 'lsmask ' ) then
                j = iarg ! getlsmask overwrites its first argument :-(
                call getlsmask(j,lsmasktype,nxmax,xxls,nymax,yyls,lwrite)
                if ( lsmasktype /= 'all' ) call checkgridequal(nx,ny,xx,yy,nxls,nyls,xxls,yyls)
            end if
        end do
    end if
    !
    call get_command_argument(command_argument_count(),maskfile)
    !
    ! fill mask array
    !
    call keepalive1('Generating mask file ',0,0)
    call fillmask(polygon,npol,pole,xx,nx,yy,ny,mask,nxmax,nymax,lwrite)
    !
    ! apply land/sea mask if requested
    !
    if ( lsmasktype == 'sea' .or. lsmasktype == 'nots' ) then
        do j=1,ny
            do i=1,nx
                if ( abs(lsmask(i,j)) > 1e-4 .eqv. lsmasktype == 'sea' ) then
                    mask(i,j) = 0
                end if
            end do
        end do
    else if ( lsmasktype == 'land' .or. lsmasktype == 'notl' ) then
        do j=1,ny
            do i=1,nx
                if ( abs(lsmask(i,j)-1) > 1e-4 .eqv. lsmasktype == 'land' ) then
                    mask(i,j) = 0
                end if
            end do
        end do
    else if ( lsmasktype == '5lan' .or. lsmasktype == '5sea' ) then
        do j=1,ny
            do i=1,nx
                if ( lsmask(i,j) < 0.5 .eqv. lsmasktype == '5lan' ) then
                    mask(i,j) = 0
                end if
            end do
        end do
    end if
    !
    ! make sure there is at least one grid box non-zero inside *or close to* each polygon
    ! this gives more equal weight to the polygons (eg Hawaii)
    !
    call getcoordspolygon(polygon,npol,coords,ncoords,lwrite)
    call setmasknonzero(coords,ncoords,xx,nx,yy,ny,xwrap,mask,nxmax,nymax,lwrite)
    !   
    !
    ! and write to file
    !
    nvars = 1
    units(1) = '1'
    vars(1) = 'mask'
    lvars(1) = 'mask from polygon in '//trim(datfile)
    nt = 1
    undef = 3e33
    title = 'mask derived from '//trim(datfile)
    nens1 = 0
    nens2 = 0
    ntmax = 1
    nz = 0
    jvars(1,1) = 0
    ! the convention in the Climate Explorer is that fields that apply to all months are in Dec0000
    firstyr = 0
    firstmo = 12
    call writenc(maskfile,ncid,ntvarid,itimeaxis,ntmax,nx,xx,ny,yy &
        ,nz,zz,nt,nperyear,firstyr,firstmo,undef,title,nvars &
        ,vars,jvars,lvars,units,nens1,nens2)
    call writencslice(ncid,ntvarid,itimeaxis,ntmax,jvars,mask &
        ,nxmax,nymax,nzmax,nx,ny,nz,1,0)
    i = nf_close(ncid) ! otherwise the cache is not fluched to disk :-(

end program polygon2mask

subroutine fillmask(polygon,npol,pole,xx,nx,yy,ny,mask,nxmax,nymax,lwrite)
!
! fills the mask(i,j) with a lat/lon grid of reals:
!
! 0 if the grid point is outside the polygon
! 1 if it is inside
! 0.5 if it is on the edge
! a very simplistic crossing-counting algorithm is used that
! gives funny results if the polygon crosses itself.
!
    implicit none
    integer npol,nx,ny,nxmax,nymax
    double precision polygon(2,npol)
    real xx(nx),yy(ny),mask(nxmax,nymax)
    character pole*2
    logical lwrite
    integer ipol,ix,iy,ixmin,ixmax,iymin,iymax,count,countmax
    double precision xmin,xmax,ymin,ymax,res(2),epsilon
    real in_polygon

    mask = 0
    if ( .false. ) then ! too many exceptions... (eg Arctic)
        ! THIS CODE DOES NOT WORK AND IS NOT EXECUTED
        ! search for lowest/highest latitude to save time in the subsequent loops
        !
        ymin = max(yy(1),yy(ny))
        ymax = min(yy(1),yy(ny))
        do ipol=1,npol - 1 ! last point repeats the first one
            if ( polygon(2,ipol) < 1e33 ) then
                ymin = min(ymin,polygon(2,ipol))
                ymax = max(ymax,polygon(2,ipol))
            end if
        end do
        if ( yy(ny) > yy(1) ) then ! bloody reversed grids...
            do iy = 1,ny
                if ( yy(iy) >= ymin ) exit
            end do
            iymin = iy
            do iy = ny,1,-1
                if ( yy(iy) <= ymax ) exit
            end do
            iymax = iy
        else
            do iy = ny,1,-1
                if ( yy(iy) >= ymin ) exit
            end do
            iymax = iy
            do iy = 1,ny
                if ( yy(iy) <= ymax ) exit
            end do
            iymin = iy
        end if
        if ( lwrite ) then
            print *,'fillmask: restricting computation to latitude range'
            print *,iymin,yy(iymin)
            print *,iymax,yy(iymax)
        end if
        !
        ! do not search for lowest/highest longitude to save time in the subsequent loops
        !
        xmin = max(xx(1),xx(nx))
        xmax = min(xx(1),xx(nx))
        do ipol=1,npol - 1 ! last point repeats the first one
            if ( polygon(1,ipol) < 3e33 ) then
                xmin = min(xmin,polygon(1,ipol))
                xmax = max(xmax,polygon(1,ipol))
            end if
        end do
        if ( xmin < xx(1) .or. xmax > xx(nx) ) then
            ! probably wraps, take all points
            ixmin = 1
            ixmax = nx
        else
            do ix = 1,nx
                if ( xx(ix) >= xmin ) exit
            end do
            ixmin = ix
            do ix = nx,1,-1
                if ( xx(ix) <= xmax ) exit
            end do
            ixmax = ix
        end if
        if ( lwrite ) then
             print *,'fillmask: restricting computation to longitude range'
             print *,ixmin,xx(ixmin)
             print *,ixmax,xx(ixmax)
        end if
    else
        ! KEEP IT SIMPLE FOR THE TIME BEING...
        ixmin = 1
        ixmax = nx
        iymin = 1
        iymax = ny
    end if
    !
    ! loop over grid
    !
    do ix=ixmin,ixmax
        call keepalive1('Grid point',ix-ixmin+1,ixmax-ixmin+1)
        do iy=iymin,iymax
            mask(ix,iy) = in_polygon(polygon,npol,dble(xx(ix)),dble(yy(iy)),pole,lwrite)
            if ( lwrite ) print '(a,i5,f9.3,i5,f8.3,f4.1)', &
               'Found point ',ix,xx(ix),iy,yy(iy),mask(ix,iy)
            if ( mask(ix,iy) == 3e33 ) then
                ! ambiguous - shift point by epsilon
                ! I expect this to happen a lot as both polygons and grids tend to be on integer numbers
                epsilon = (xx(2) - xx(1))/100
                res(1) = -3e33
                res(2) = +3e33
                do while ( res(1) /= res(2) )
                    if ( lwrite ) print *,'Found undefined point, trying with offset ',epsilon
                    res(1) = in_polygon(polygon,npol,dble(xx(ix))-epsilon,dble(yy(iy)),pole,lwrite)
                    if ( lwrite ) print *,'   x- point ',ix,xx(ix)-epsilon,iy,yy(iy),res(1)
                    if ( res(1) == 3e33 ) then
                        res(1) = in_polygon(polygon,npol,dble(xx(ix))-epsilon,dble(yy(iy))-epsilon/sqrt(3.),pole,lwrite)
                        if ( lwrite ) print *,'  xy- point ',ix,xx(ix)-epsilon,iy,yy(iy)-epsilon/sqrt(3.),res(1)
                        if ( res(1) == 3e33 ) then
                            res(1) = in_polygon(polygon,npol,dble(xx(ix)),dble(yy(iy))-epsilon/sqrt(3.),pole,lwrite)
                            if ( lwrite ) print *,'   y- point ',ix,xx(ix),iy,yy(iy)-epsilon/sqrt(3.),res(1)
                        end if
                    end if
                    res(2) = in_polygon(polygon,npol,dble(xx(ix))+epsilon,dble(yy(iy)),pole,lwrite)
                    if ( lwrite ) print *,'   x+ point ',ix,xx(ix)+epsilon,iy,yy(iy),res(2)
                    if ( res(2) == 3e33 ) then
                        res(2) = in_polygon(polygon,npol,dble(xx(ix))+epsilon,dble(yy(iy))+epsilon/sqrt(3.),pole,lwrite)
                        if ( lwrite ) print *,'  xy+ point ',ix,xx(ix)+epsilon,iy,yy(iy)+epsilon/sqrt(3.),res(2)
                        if ( res(2) == 3e33 ) then
                            res(2) = in_polygon(polygon,npol,dble(xx(ix)),dble(yy(iy))+epsilon/sqrt(3.),pole,lwrite)
                            if ( lwrite ) print *,'   y+ point ',ix,xx(ix),iy,yy(iy)+epsilon/sqrt(3.),res(2)
                        end if
                    end if
                    if ( (res(1) == 3e33 .or. res(2) == 3e33) .and. epsilon >= 1e-4 ) then
                        ! try again, we were unlucky
                        res(1) = -1
                        res(2) = +1
                    end if 
                    epsilon = epsilon/3.14 ! not a nice number
                    if ( epsilon < 1e-4 ) exit
                end do
                if ( res(1) == 3e33 .or. res(2) == 3e33 .or. res(1) /= res(2) ) then
                    write(0,*) 'fillmask: internal error ',res,epsilon
                    write(0,*) '          for point ',ix,iy,xx(ix),yy(iy)
                    call exit(-1)
                end if
                mask(ix,iy) = res(1)
            end if
        end do
    end do
end subroutine fillmask

subroutine getcoordspolygon(polygon,npol,coords,ncoords,lwrite)
    !
    ! get central coordinates of polygons in list
    ! for the time being just the centre of the boundingbox
    !
    implicit none
    integer npol,nmax,ncoords
    double precision polygon(2,npol),coords(3,npol)
    logical lwrite
    integer ipol,i
    double precision xmin,xmax,ymin,ymax,x,xold,area

    xmin = 3e33
    xmax = -3e33
    xold = 3e33
    ymin = 3e33
    ymax = -3e33
    area = 0
    ncoords = 0
    do ipol=1,npol
        if ( polygon(1,ipol) > 1e33 ) then
            ncoords = ncoords + 1
            coords(1,ncoords) = (xmin + xmax)/2
            coords(2,ncoords) = (ymin + ymax)/2
            coords(3,ncoords) = abs(area)/2
            if ( lwrite ) then
                print *,'getcoordspolygon: n,x,y,area = ',ncoords,(coords(i,ncoords),i=1,3)
            end if
            xmin = 3e33
            xmax = -3e33
            xold = 3e33
            ymin = 3e33
            ymax = -3e33
            area = 0
        else
            x = polygon(1,ipol)
            if ( xold < 1e33 ) then
                if ( abs(x-xold) > abs(x-xold+360) ) then
                    x = x + 360
                else if ( abs(x-xold) > abs(x-xold-360) ) then
                    x = x - 360
                end if
            end if
            xmin = min(xmin,x)
            xmax = max(xmax,x)
            ymin = min(ymin,polygon(2,ipol))
            ymax = max(ymax,polygon(2,ipol))
            if ( xold < 1e33 ) then
                ! http://www.mathopenref.com/coordpolygonarea.html
                area = area + xold*polygon(2,ipol) - x*polygon(2,ipol-1)
            end if
            xold = x
        end if
    end do
    ! and the last polygon
    if ( xmin < 1e33 ) then
        ncoords = ncoords + 1
        coords(1,ncoords) = (xmin + xmax)/2
        coords(2,ncoords) = (ymin + ymax)/2
        coords(3,ncoords) = abs(area)/2
        if ( lwrite ) then
            print *,'getcoordspolygon: n,x,y,area = ',ncoords,(coords(i,ncoords),i=1,3)
        end if
    end if

end subroutine getcoordspolygon

subroutine setmasknonzero(coords,ncoords,xx,nx,yy,ny,xwrap,mask,nxmax,nymax,lwrite)
    !
    !   set mask to nonzero for the central point of the polygon if there is nothing there
    !
    implicit none
    integer ncoords,nx,ny,nxmax,nymax
    double precision coords(3,ncoords)
    real xx(nx),yy(ny),mask(nxmax,nymax)
    logical xwrap,lwrite
    integer icoords,ix,iy
    real x,dist,pi,dx,dy
    integer,allocatable :: ijmin(:,:)
    real,allocatable :: distmin(:)

    pi = 4*atan(1.)
    allocate(ijmin(2,ncoords))
    allocate(distmin(ncoords))
    ijmin = -999
    distmin = 3e33
    do iy=1,ny
        do ix=1,nx
            do icoords=1,ncoords
                x = xx(ix)
                if ( abs(x-coords(1,icoords)+360) < abs(x-coords(1,icoords)) ) then
                    x = x+360
                else if ( abs(x-coords(1,icoords)-360) < abs(x-coords(1,icoords)) ) then
                    x = x-360
                end if
                dist = sqrt((cos(yy(iy)/180*pi)*(x-coords(1,icoords)))**2 + &
                    & (yy(iy)-coords(2,icoords))**2) ! small distances so approx flat
                if ( dist < distmin(icoords) ) then
                    ijmin(1,icoords) = ix
                    ijmin(2,icoords) = iy
                    distmin(icoords) = dist
                end if
            end do
        end do
    end do
    do icoords=1,ncoords
        ix = ijmin(1,icoords)
        iy = ijmin(2,icoords)
        if ( mask(ix,iy) < 0.5 ) then
            if ( xwrap ) then
                if ( ix == 1 ) then
                    dx = (xx(2) - xx(nx) + 360)/2
                else if ( ix == nx ) then
                    dx = (xx(1) - xx(nx-1) + 360)/2
                else
                    dx = (xx(ix+1) - xx(ix-1))/2
                end if        
            else
                if ( ix == 1 ) then
                    dx = xx(2) - xx(1)
                else if ( ix == nx ) then
                    dx = xx(nx) - xx(nx-1)
                else
                    dx = (xx(ix+1) - xx(ix-1))/2
                end if        
            end if
            if ( iy == 1 ) then
                dy = yy(2) - yy(1)
            else if ( iy == ny ) then
                dy = yy(ny) - yy(ny-1)
            else
                dy = (yy(iy+1) - yy(iy-1))/2
            end if
            ! only if the distance is smaller than half the diagonal of the grid box -
            ! necessary for regional grids
            if ( distmin(icoords) < sqrt((cos(yy(iy)/180*pi)*dx)**2 + dy**2)/2 ) then
                mask(ix,iy) = mask(ix,iy) + coords(3,icoords)/abs(dx*dy)
                ! area of polygon divided by grid box area, no factors cos(lat)
                mask(ix,iy) = min(1.,mask(ix,iy)) ! Sulawesi is bigger than one grid box
                if ( lwrite ) then
                    print *,'setmasknonzero: ix,iy,x,y,dx,dy = ',ix,iy,xx(ix),yy(iy),dx,dy
                    print *,'                icoords,distmin = ',icoords,distmin(icoords)
                    print *,'                area,value      = ',coords(3,icoords),coords(3,icoords)/abs(dx*dy)
                    print *,'                mask            = ',mask(ix,iy)
                end if
            end if
        end if
    end do
end subroutine setmasknonzero
