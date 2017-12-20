program polygon2mask
    !
    ! converts a polygon to a mask, given a grid
    !
    use lsdata
    implicit none
    include 'params.h'
    include 'netcdf.inc'
    integer,parameter :: npolmax=50000,nvarmax=10
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
    integer iargc

    lwrite = .false.
    lsmasktype = 'all '
    if ( iargc().lt.3 ) then
        write(0,*) 'usage: polygon2mask grid.nc polygon.txt [sp] '// &
            & '[lsmask file.nc all|land|sea|notland|notsea|all] [debug] mask.nc'
        write(0,*) 'with grid.nc a netcdf file that defines the grid'
        write(0,*) 'with polygon.txt a text file with xi yi on each ine'
        write(0,*) 'sp denotes that the South Pole is inside the polygon'
        write(0,*) 'mask.nc is 1 when the grid point is inside the polygon, '
        write(0,*) '0 when outside (to be exact odd/even no of crossings).'
        stop
    end if
    do i=3,iargc()-1
     call getarg(i,string)
     call tolower(string)
     if ( string.eq.'lwrite' .or. string.eq.'debug') then
        lwrite = .true.
     end if
  end do
  !
  ! read grid xx,yy from netcdf (or grads) file
  !
  call getarg(1,infile)
  call getmetadata(infile,mens1,mens,ncid,datfile,nxmax,nx &
       & ,xx,nymax,ny,yy,nzmax,nz,zz,lz,nt,nperyear,firstyr,firstmo &
       & ,ltime,undef,endian,title,history,nvarmax,nvars,vars,jvars &
       & ,lvars,svars,units,cell_methods,metadata,lwrite)
  i = nf_close(ncid)
  if ( ny.le.1 ) then
     write(0,*) 'polygon2mask: error: expecting ta least two grid points in latitude'
     call abort
  end if
  if ( nx.le.1 ) then
     write(0,*) 'polygon2mask: error: expecting ta least two grid points in longitude'
     call abort
  end if
  call getxyprop(xx,nx,yy,ny,xrev,yrev,xwrap)
  !
  ! read polygon from text file
  !
  call getarg(2,datfile)
  call read_polygon(datfile,npol,npolmax,polygon,lwrite)
  !
  if ( iargc().gt.3 ) then
     call getarg(3,pole)
     if ( pole.ne.'sp' ) pole='np'
     do iarg=3,iargc()-1
        call getarg(iarg,string)
        if ( string(1:6).eq.'lsmask' ) then
           j = iarg ! getlsmask overwrites its first argument :-(
           call getlsmask(j,lsmasktype,nxmax,xxls,nymax,yyls,lwrite)
           if ( lsmasktype.ne.'all' ) call checkgridequal(nx,ny,xx,yy,nxls,nyls,xxls,yyls)
        end if
     end do
  end if
  !
  call getarg(iargc(),maskfile)
  !
  ! fill mask array
  !
  call keepalive1('Generating mask file ',0,0)
  call fillmask(polygon,npol,pole,xx,nx,yy,ny,mask,nxmax,nymax,lwrite)
  !
  ! apply land/sea mask if requested
  !
  if ( lsmasktype.eq.'sea' .or. lsmasktype.eq.'nots' ) then
     do j=1,ny
        do i=1,nx
           if ( abs(lsmask(i,j)).gt.1e-4 .eqv. lsmasktype.eq.'sea' ) then
              mask(i,j) = 0
           end if
        end do
     end do
  else if ( lsmasktype.eq.'land' .or. lsmasktype.eq.'notl' ) then
     do j=1,ny
        do i=1,nx
           if ( abs(lsmask(i,j)-1).gt.1e-4 .eqv. lsmasktype.eq.'land' ) then
              mask(i,j) = 0
           end if
        end do
     end do
  else if ( lsmasktype.eq.'5lan' .or. lsmasktype.eq.'5sea' ) then
     do j=1,ny
        do i=1,nx
           if ( lsmask(i,j).lt.0.5 .eqv. lsmasktype.eq.'5lan' ) then
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
       &  ,nz,zz,nt,nperyear,firstyr,firstmo,undef,title,nvars &
       & ,vars,jvars,lvars,units,nens1,nens2)
  call writencslice(ncid,ntvarid,itimeaxis,ntmax,jvars,mask &
       & ,nxmax,nymax,nzmax,nx,ny,nz,1,0)
  i = nf_close(ncid) ! otherwise the cache is not fluched to disk :-(

end program polygon2mask

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
    integer npol1,npol2,i
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
        if ( string.eq.' ' ) then
            ! signifies the beginning of a new polygon
            leof = .false.
            if ( npol.eq.0 ) then
                goto 10 ! skip empty first line
            else
                goto 100 ! new polygon
            end if
        end if
        if ( string(1:1).eq.'#' .or. string(2:2).eq.'#' ) cycle
        npol = npol + 1
        ! excel often uses commas rather than decimal points
        if ( index(string,'.') == 0 ) then
            do i=1,len(string)
                if ( string(i:i) == ',' ) string(i:i) = '.'
            end do
        end if
        read(string,*) polygon(1,npol),polygon(2,npol)
        if ( polygon(2,npol).lt.-90 .or. polygon(2,npol).gt.90 ) then
            write(0,*) 'polygon2mask: error: latitude ',npol,' is ',polygon(2,npol)
            call abort
        end if
    end do
100 continue
    if ( polygon(1,npol1).ne.polygon(1,npol) .or. &
        & polygon(2,npol1).ne.polygon(2,npol) ) then
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
end subroutine read_polygon

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
         if ( polygon(2,ipol).lt.1e33 ) then
            ymin = min(ymin,polygon(2,ipol))
            ymax = max(ymax,polygon(2,ipol))
         end if
      end do
      if ( yy(ny).gt.yy(1) ) then ! bloody reversed grids...
         do iy = 1,ny
            if ( yy(iy).ge.ymin ) exit
         end do
         iymin = iy
         do iy = ny,1,-1
            if ( yy(iy).le.ymax ) exit
         end do
         iymax = iy
      else
         do iy = ny,1,-1
            if ( yy(iy).ge.ymin ) exit
         end do
         iymax = iy
         do iy = 1,ny
            if ( yy(iy).le.ymax ) exit
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
         if ( polygon(1,ipol).lt.3e33 ) then
            xmin = min(xmin,polygon(1,ipol))
            xmax = max(xmax,polygon(1,ipol))
         end if
      end do
      if ( xmin.lt.xx(1) .or. xmax.gt.xx(nx) ) then
         ! probably wraps, take all points
         ixmin = 1
         ixmax = nx
      else
         do ix = 1,nx
            if ( xx(ix).ge.xmin ) exit
         end do
         ixmin = ix
         do ix = nx,1,-1
            if ( xx(ix).le.xmax ) exit
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
        &   'Found point ',ix,xx(ix),iy,yy(iy),mask(ix,iy)
        if ( mask(ix,iy).eq.3e33 ) then
           ! ambiguous - shift point by epsilon
           ! I expect this to happen a lot as both polygons and grids tend to be on integer numbers
           epsilon = (xx(2) - xx(1))/100
           res(1) = -3e33
           res(2) = +3e33
           do while ( res(1).ne.res(2) )
              if ( lwrite ) print *,'Found undefined point, trying with offset ',epsilon
              res(1) = in_polygon(polygon,npol,dble(xx(ix))-epsilon,dble(yy(iy)),pole,lwrite)
              if ( lwrite ) print *,'   x- point ',ix,xx(ix)-epsilon,iy,yy(iy),res(1)
              if ( res(1).eq.3e33 ) then
                 res(1) = in_polygon(polygon,npol,dble(xx(ix))-epsilon,dble(yy(iy))-epsilon/sqrt(3.),pole,lwrite)
                 if ( lwrite ) print *,'  xy- point ',ix,xx(ix)-epsilon,iy,yy(iy)-epsilon/sqrt(3.),res(1)
                 if ( res(1).eq.3e33 ) then
                    res(1) = in_polygon(polygon,npol,dble(xx(ix)),dble(yy(iy))-epsilon/sqrt(3.),pole,lwrite)
                    if ( lwrite ) print *,'   y- point ',ix,xx(ix),iy,yy(iy)-epsilon/sqrt(3.),res(1)
                 end if
              end if
              res(2) = in_polygon(polygon,npol,dble(xx(ix))+epsilon,dble(yy(iy)),pole,lwrite)
              if ( lwrite ) print *,'   x+ point ',ix,xx(ix)+epsilon,iy,yy(iy),res(2)
              if ( res(2).eq.3e33 ) then
                 res(2) = in_polygon(polygon,npol,dble(xx(ix))+epsilon,dble(yy(iy))+epsilon/sqrt(3.),pole,lwrite)
                 if ( lwrite ) print *,'  xy+ point ',ix,xx(ix)+epsilon,iy,yy(iy)+epsilon/sqrt(3.),res(2)
                 if ( res(2).eq.3e33 ) then
                    res(2) = in_polygon(polygon,npol,dble(xx(ix)),dble(yy(iy))+epsilon/sqrt(3.),pole,lwrite)
                    if ( lwrite ) print *,'   y+ point ',ix,xx(ix),iy,yy(iy)+epsilon/sqrt(3.),res(2)
                 end if
              end if
              if ( (res(1).eq.3e33 .or. res(2).eq.3e33) .and. epsilon.ge.1e-4 ) then
                  ! try again, we were unlucky
                  res(1) = -1
                  res(2) = +1
              end if 
              epsilon = epsilon/3.14 ! not a nice number
              if ( epsilon < 1e-4 ) exit
           end do
           if ( res(1).eq.3e33 .or. res(2).eq.3e33 .or. res(1).ne.res(2) ) then
              write(0,*) 'fillmask: internal error ',res,epsilon
              write(0,*) '          for point ',ix,iy,xx(ix),yy(iy)
              call abort
           end if
           mask(ix,iy) = res(1)
        end if
     end do
  end do
end subroutine fillmask

real function in_polygon(polygon,npol,x,y,pole,lwrite)
  !
  ! checks whether the point (x,y) is inside the polygon
  ! if pole.eq.'sp' then the South Pole is inside the polygon, otherwise outside
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
  if ( pole.eq.'sp') then
     result = 1
  else
     result = 0
  end if
  if ( y.eq.ynull ) then
     in_polygon = result
     return
  end if
!
  if ( .false. .and. lwrite ) print *,'started with ',result
  do ipol=1,npol-1
     ! skip "segments" that include the polygon-separating marker 3e33,3e33
     if ( polygon(1,ipol).gt.1e33 ) cycle
     if ( polygon(1,ipol+1).gt.1e33 ) cycle
     llwrite = .false.
     !!!if ( x.eq.14.75 .and. y.gt.48.5 .and. &
     !!!& (polygon(1,ipol)-x)*(polygon(1,ipol+1)-x).lt.0 ) llwrite = .true.
     partial = segments_crossed(polygon(1:2,ipol),polygon(1:2,ipol+1),x,y,ynull,llwrite)
     if ( lwrite .and. partial.ne.0 ) print '(a,i6,f3.0,4f10.5)','in_polygon: segments_crossed returns ', &
     & ipol,partial,polygon(1:2,ipol),polygon(1:2,ipol+1)
     if ( partial.eq.0.5 ) then
        result = partial
        if ( lwrite ) print *,'end point on line ',result
        exit
     end if
     if ( partial.eq.1 ) then
        result = 1-result
        if ( .false. .and. lwrite ) print *,'flipped from ',1-result,' to ',result
     else if ( partial.eq.3e33 ) then
        if ( lwrite ) then
        	print *,'undefined ',result
			partial = segments_crossed(polygon(1:2,ipol),polygon(1:2,ipol+1),x,y,ynull,lwrite)
		end if
        foundtouch = .true.
     end if
  end do
  if ( foundtouch .and. result.ne.0.5 ) then
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
  if ( abs(ynull).ne.90 ) then
     write(0,*) 'segments_crossed: error: ynull should be +/-90, not ',ynull
     call abort
  end if
  !
  ! first reduce the longitude to the shortest possible form assuming only one wrap off
  !
  x1 = p1(1)
  x2 = p2(1)
  x = xin
  y = yin
  if ( abs(x2-x1).gt.abs(x2-360-x1) .and. abs(x2-360-x1).gt.0.01 ) then
     x2 = x2 - 360
  end if
  if ( abs(x2-x1).gt.abs(x2+360-x1) .and. abs(x2+360-x1).gt.0.01 ) then
     x2 = x2 + 360
  end if
  if ( abs(x-(x1+x2)/2).gt.abs(x-360-(x1+x2)/2) ) then
     x = x - 360
  end if
  if ( abs(x-(x1+x2)/2).gt.abs(x+360-(x1+x2)/2) ) then
     x = x + 360
  end if
  y1 = p1(2)
  y2 = p2(2)
  !
  ! a few special cases that will cause grief later on
  !
  if ( abs(y1).eq.90 .and. abs(y2).eq.90 ) then
     if ( lwrite ) print *,'segment is one point on pole'
     segments_crossed = 0
     return
  end if
  if ( lwrite ) then
     print *,'special case? x1,x2,x = ',x1,x2,x
  end if
  if ( abs(x1-x2).lt.0.00000001 ) then
     if ( x.ne.x1 ) then
        segments_crossed = 0
     else if ( ynull.eq.-90 .and. min(y1,y2).gt.y ) then
        segments_crossed = 0 ! north of the grid point
     else if ( ynull.eq.-90 .and. max(y1,y2).gt.y ) then
        segments_crossed = 0.5 ! on the line
     else if ( ynull.eq.+90 .and. max(y1,y2).lt.y ) then
        segments_crossed = 0 ! south of the grid point
     else if ( ynull.eq.+90 .and. min(y1,y2).lt.y ) then
        segments_crossed = 0.5 ! on the line
     else
        segments_crossed = 3e33 ! no idea whether it crossed or not
     end if
     if ( lwrite ) then
        print *,'special case x1=x2, segments_crossed = ',segments_crossed
     end if
     return
  end if
  if ( y.eq.ynull ) then
     segments_crossed = 0.5
     return
  end if
  ! position along p1 - p2
  lam = (x-x1)/(x2-x1)
  if ( lwrite ) print *,'lam = ',lam
  if ( lam.lt.0 .or. lam.gt.1 ) then
     segments_crossed = 0
  else
     ! position along (x,ynull) - (x,y)
     mu = (y1-ynull + lam*(y2-y1))/(y-ynull)
     if ( mu.gt.1 ) then
        segments_crossed = 0
     else if ( mu.eq.1 ) then
        segments_crossed = 0.5
     else
        if ( lam.eq.0 .or. lam.eq.1 ) then
           segments_crossed = 3e33
        else
           segments_crossed = 1
        end if
     end if      
  end if
end function segments_crossed

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
     if ( polygon(1,ipol).gt.1e33 ) then
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
        if ( xold.lt.1e33 ) then
            if ( abs(x-xold).gt.abs(x-xold+360) ) then
                x = x + 360
            else if ( abs(x-xold).gt.abs(x-xold-360) ) then
                x = x - 360
            end if
        end if
        xmin = min(xmin,x)
        xmax = max(xmax,x)
        ymin = min(ymin,polygon(2,ipol))
        ymax = max(ymax,polygon(2,ipol))
        if ( xold.lt.1e33 ) then
            ! http://www.mathopenref.com/coordpolygonarea.html
            area = area + xold*polygon(2,ipol) - x*polygon(2,ipol-1)
        end if
        xold = x
     end if
  end do
  ! and the last polygon
  if ( xmin.lt.1e33 ) then
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
                if ( abs(x-coords(1,icoords)+360).lt.abs(x-coords(1,icoords)) ) then
                    x = x+360
                else if ( abs(x-coords(1,icoords)-360).lt.abs(x-coords(1,icoords)) ) then
                    x = x-360
                end if
                dist = sqrt((cos(yy(iy)/180*pi)*(x-coords(1,icoords)))**2 + &
                    & (yy(iy)-coords(2,icoords))**2) ! small distances so approx flat
                if ( dist.lt.distmin(icoords) ) then
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
        if ( mask(ix,iy).lt.0.5 ) then
            if ( xwrap ) then
                if ( ix.eq.1 ) then
                    dx = (xx(2) - xx(nx) + 360)/2
                else if ( ix.eq.nx ) then
                    dx = (xx(1) - xx(nx-1) + 360)/2
                else
                    dx = (xx(ix+1) - xx(ix-1))/2
                end if        
            else
                if ( ix.eq.1 ) then
                    dx = xx(2) - xx(1)
                else if ( ix.eq.nx ) then
                    dx = xx(nx) - xx(nx-1)
                else
                    dx = (xx(ix+1) - xx(ix-1))/2
                end if        
            end if
            if ( iy.eq.1 ) then
                dy = yy(2) - yy(1)
            else if ( iy.eq.ny ) then
                dy = yy(ny) - yy(ny-1)
            else
                dy = (yy(iy+1) - yy(iy-1))/2
            end if
            ! only if the distance is smaller than half the diagonal of the grid box -
            ! necessary for regional grids
            if ( distmin(icoords).lt.sqrt((cos(yy(iy)/180*pi)*dx)**2 + dy**2)/2 ) then
                mask(ix,iy) = mask(ix,iy) + &
                    & coords(3,icoords)/abs(dx*dy)
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
