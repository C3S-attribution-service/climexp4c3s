program polygon2mask
  !
  ! converts a polygon to a mask, given a grid
  !
  use lsdata
  implicit none
  include 'params.h'
  include 'netcdf.inc'
  integer,parameter :: npolmax=10000,nvarmax=10
  integer i,j,npol,mens,mens1,ncid,nx,ny,nz,nt,nperyear,firstyr,firstmo,endian,nvars,jvars(6,nvarmax),ntvarid
  integer nens1,nens2,ntmax,itimeaxis(1),iarg
  real xx(nxmax),yy(nymax),zz(nzmax),undef,xxls(nxmax),yyls(nymax)
  real polygon(2,npolmax),mask(nxmax,nymax)
  character infile*255,datfile*255,maskfile*255,lz(3)*20,ltime*120,title*1047,history*4095, &
       & string*80,lsmasktype*4
  character vars(nvarmax)*50,lvars(nvarmax)*100,svars(nvarmax)*100,units(nvarmax)*100 &
       & ,cell_methods(nvarmax)*100
  character pole*2,clwrite*1
  logical lwrite
  integer iargc

  lwrite = .false.
  call getenv('POLYGON2MAsK_LWRITE',clwrite)
  if ( clwrite.eq.'T' .or. clwrite.eq.'t' ) lwrite = .true.
  lsmasktype = 'all '
  if ( iargc().lt.3 ) then
     write(0,*) 'usage: polygon2mask grid.nc polygon.txt [sp] '// &
          & '[lsmask file.nc all|land|sea|notland|notsea|all] mask.nc'
     write(0,*) 'with grid.nc a netcdf file that defines the grid'
     write(0,*) 'with polygon.txt a text file with xi yi on each ine'
     write(0,*) 'sp denotes that the South Pole is inside the polygon'
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
       & ,lvars,svars,units,cell_methods,lwrite)
  i = nf_close(ncid)
  if ( ny.le.1 ) then
     write(0,*) 'polygon2mask: error: expecting ta least two grid points in latitude'
     call abort
  end if
  if ( nx.le.1 ) then
     write(0,*) 'polygon2mask: error: expecting ta least two grid points in longitude'
     call abort
  end if
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
  implicit none
  integer npol,npolmax
  real polygon(2,npolmax)
  character datfile*(*)
  integer ipol
  character string*100
  logical lwrite
!
  npol = 0
  polygon = 3e33
  open(1,file=trim(datfile),status='old')
  do
     read(1,'(a)',end=100) string
     if ( string.eq.' ' ) cycle
     if ( string(1:1).eq.'#' ) cycle
     npol = npol + 1
     read(string,*) polygon(1,npol),polygon(2,npol)
     if ( polygon(2,npol).lt.-90 .or. polygon(2,npol).gt.90 ) then
        write(0,*) 'polygon2mask: error: latitude ',npol,' is ',polygon(2,npol)
        call abort
     end if
  end do
100 continue
  if ( polygon(1,1).ne.polygon(1,npol) .or. &
       & polygon(2,1).ne.polygon(2,npol) ) then
     ! close loop
     npol = npol + 1
     polygon(1,npol) = polygon(1,1)
     polygon(2,npol) = polygon(2,1)
  end if
  if ( lwrite ) then
     print *,'read_polygon: found ',npol-1,'-edged polygon'
     do ipol=1,npol
        print *,ipol,polygon(1,ipol),polygon(2,ipol)
     end do
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
  real polygon(2,npol),xx(nx),yy(ny),mask(nxmax,nymax)
  character pole*2
  logical lwrite
  integer ipol,ix,iy,ixmin,ixmax,iymin,iymax
  real xmin,xmax,ymin,ymax,res(2),epsilon
  real in_polygon

  mask = 0
  if ( .false. ) then ! too many exceptions... (eg Arctic)
  !
  ! search for lowest/highest latitude to save time in the subsequent loops
  !
  ymin = max(yy(1),yy(ny))
  ymax = min(yy(1),yy(ny))
  do ipol=1,npol - 1 ! last point repeats the first one
     ymin = min(ymin,polygon(2,ipol))
     ymax = max(ymax,polygon(2,ipol))
  end do
  if ( yy(ny).gt.yy(1) ) then ! bloody rersved grids...
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
     xmin = min(xmin,polygon(1,ipol))
     xmax = max(xmax,polygon(1,ipol))
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
     ixmin = 1
     ixmax = nx
     iymin = 1
     iymax = ny
  end if
  !
  ! loop over grid
  !
  do ix=ixmin,ixmax
     do iy=iymin,iymax
        mask(ix,iy) = in_polygon(polygon,npol,xx(ix),yy(iy),pole,lwrite)
        if ( lwrite ) print *,'Found point ',ix,xx(ix),iy,yy(iy),mask(ix,iy)
        if ( mask(ix,iy).eq.3e33 ) then
           ! ambiguous - shift point by epsilon
           ! I expect this to happen a lot as both polygons and grids tend to be on integer numbers
           epsilon = (xx(2) - xx(1))/100
           res(1) = -3e33
           res(2) = +3e33
           do while ( res(1).ne.res(2) )
              if ( lwrite ) print *,'Found undefined point, trying with offset ',epsilon
              res(1) = in_polygon(polygon,npol,xx(ix)-epsilon,yy(iy),pole,lwrite)
              if ( lwrite ) print *,'   x- point ',ix,xx(ix)-epsilon,iy,yy(iy),res(1)
              if ( res(1).eq.3e33 ) then
                 res(1) = in_polygon(polygon,npol,xx(ix)-epsilon,yy(iy)-epsilon/sqrt(3.),pole,lwrite)
                 if ( lwrite ) print *,'  xy- point ',ix,xx(ix)-epsilon,iy,yy(iy)-epsilon/sqrt(3.),res(1)
                 if ( res(1).eq.3e33 ) then
                    res(1) = in_polygon(polygon,npol,xx(ix),yy(iy)-epsilon/sqrt(3.),pole,lwrite)
                    if ( lwrite ) print *,'   y- point ',ix,xx(ix),iy,yy(iy)-epsilon/sqrt(3.),res(1)
                 end if
              end if
              res(2) = in_polygon(polygon,npol,xx(ix)+epsilon,yy(iy),pole,lwrite)
              if ( lwrite ) print *,'   x+ point ',ix,xx(ix)+epsilon,iy,yy(iy),res(2)
              if ( res(2).eq.3e33 ) then
                 res(2) = in_polygon(polygon,npol,xx(ix)+epsilon,yy(iy)+epsilon/sqrt(3.),pole,lwrite)
                 if ( lwrite ) print *,'  xy+ point ',ix,xx(ix)+epsilon,iy,yy(iy)+epsilon/sqrt(3.),res(2)
                 if ( res(2).eq.3e33 ) then
                    res(2) = in_polygon(polygon,npol,xx(ix),yy(iy)+epsilon/sqrt(3.),pole,lwrite)
                    if ( lwrite ) print *,'   y+ point ',ix,xx(ix),iy,yy(iy)+epsilon/sqrt(3.),res(2)
                 end if
              end if
              epsilon = epsilon/10
           end do
           if ( res(1).eq.3e33 ) then
              write(0,*) 'fillmask: internal error ',res
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
  real polygon(2,npol),x,y
  character pole*(*)
  logical lwrite
  integer ipol
  real ynull,result,partial
  real segments_crossed
  logical foundtouch

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
  if ( lwrite ) print *,'started with ',result
  do ipol=1,npol-1
     partial = segments_crossed(polygon(1:2,ipol),polygon(1:2,ipol+1),x,y,ynull,.false.)
     if ( lwrite .and. partial.ne.0 ) print *,'in_polygon: segments_crossed returns ',ipol,partial
     if ( partial.eq.0.5 ) then
        result = partial
        if ( lwrite ) print *,'end point on line ',result
        exit
     end if
     if ( partial.eq.1 ) then
        result = 1-result
        if ( lwrite ) print *,'flipped from ',1-result,' to ',result
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

real function segments_crossed(p1,p2,xin,y,ynull,lwrite)
  !
  ! determines whther the two line segments p1 - p2 and (x,ynull) - (x,y) cross or not
  ! 0: no cross
  ! 1: crossed inside line segment
  ! 0.5: touched at (x,y)
  ! 3e33: touched at another point
  !
  implicit none
  real p1(2),p2(2),xin,y,ynull
  logical lwrite
  real x,x1,y1,x2,y2,lam,mu
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
  if ( abs(x1-x2).lt.0.001 ) then
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
