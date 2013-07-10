        subroutine interpu(field1,xx1in,yy1in,nx1,ny1,field2,xx2in,yy2in
     +       ,nx2,ny2,xx,nx,yy,ny,yr11,yr21,yr12,yr22,nxf,nyf,nperyear
     +       ,intertype,lwrite)
*
*       interpolate and truncate the fields so that they are
*       on a common grid.  The fields are modified in situ
*       to save on memory.  Undefined is assumed 3e33.
*       intertype: 
*       0: highest resolution grid
*      -1: lowest resolution grid
*       1: grid 1
*       2: grid 2
*
        implicit none
#include "params.h"
        integer nx1,ny1,nx2,ny2,nx,ny,yr11,yr21,yr12,yr22,nxf,nyf
     +        ,nperyear,intertype
        real xx(nxmax),yy(nymax),
     +       field1(nxf,nyf,nperyear,yr11:yr21),xx1in(nx1),yy1in(ny1),
     +       field2(nxf,nyf,nperyear,yr12:yr22),xx2in(nx2),yy2in(ny2)
        logical lwrite
        integer i,j,i1x1,i2x1,i1x2,i2x2,i1y1,i2y1,i1y2,i2y2,
     +        ix(nxmax),iy(2,nymax),yr,mo
        real xx1(nxmax),yy1(nymax),xx2(nxmax),yy2(nymax),
     +       dxmin1,dxmin2,dymin1,dymin2,lon1,lon2,lat1,lat2,d1,d2,d3,
     +       ax(nxmax),ay(nymax),f(nxmax,nymax)
        logical x1rev,x2rev,y1rev,y2rev,x1wrap,x2wrap,
     +        lintx1,lintx2,linty1,linty2
*
        if ( lwrite ) then
            print *,'interpu: input'
            if ( nx1.gt.1 ) print *,'xx1 :',xx1in(1),xx1in(2),'...'
     +           ,xx1in(ny1)
            if ( ny1.gt.1 ) print *,'yy1 :',yy1in(1),yy1in(2),'...'
     +           ,yy1in(ny1)
            if ( nx2.gt.1 ) print *,'xx2 :',xx2in(1),xx2in(2),'...'
     +           ,xx2in(ny2)
            if ( ny2.gt.1 ) print *,'yy2 :',yy2in(1),yy2in(2),'...'
     +           ,yy2in(ny2)
        endif
*
*       anything to do?
*
        if ( nx1.eq.nx2 .and. ny1.eq.ny2 ) then
            do i=1,nx1
                if ( abs(xx1in(i)-xx2in(i)).gt.1e-3 ) goto 10
            enddo
            do i=1,ny1
                if ( abs(yy1in(i)-yy2in(i)).gt.1e-3 ) goto 10
            enddo
            nx = nx1
            ny = ny1
            xx(1:nx) = xx1in(1:nx)
            yy(1:ny) = yy1in(1:ny)
            if ( lwrite ) print *,'no interpolation necessary'
            return
        endif
 10     continue
        xx1(1:nx1) = xx1in(1:nx1)
        yy1(1:ny1) = yy1in(1:ny1)
        xx2(1:nx2) = xx2in(1:nx2)
        yy2(1:ny2) = yy2in(1:ny2)
*
*       establish some properties of the grids
*
*       reversed grids?
        if ( nx1.gt.1 ) then
            x1rev = xx1(2).lt.xx1(1)
        else
            x1rev = .FALSE.
        endif
        if ( nx2.gt.1 ) then
            x2rev = xx2(2).lt.xx2(1)
        else
            x2rev = .FALSE.
        endif
        if ( ny1.gt.1 ) then
            y1rev = yy1(2).lt.yy1(1)
        else
            y1rev = .FALSE.
        endif
        if ( ny2.gt.1 ) then
            y2rev = yy2(2).lt.yy2(1)
        else
            y2rev = .FALSE.
        endif
*       save myself a lot of headaches later on
        if ( x1rev ) call revit(xx1,nx1)
        if ( x2rev ) call revit(xx2,nx2)
        if ( y1rev ) call revit(yy1,ny1)
        if ( y2rev ) call revit(yy2,ny2)
*       grid wraps around the earth?  (I assume it is in degrees)
        x1wrap = abs(2*xx1(nx1)-xx1(nx1-1)-360-xx1(1)).lt.1e-3
        x2wrap = abs(2*xx2(nx2)-xx2(nx2-1)-360-xx2(1)).lt.1e-3
        if ( lwrite ) then
            print *,'interpu grid properties'
            if ( x1rev ) print *,'grid1 x reversed'
            if ( x2rev ) print *,'grid2 x reversed'
            if ( y1rev ) print *,'grid1 y reversed'
            if ( y2rev ) print *,'grid2 y reversed'
            if ( x1wrap ) print *,'grid1 x wraps'
            if ( x2wrap ) print *,'grid2 x wraps'
        endif
*
*       cut-offs: longitude
*
        if ( x1wrap ) then
            i1x2 = 1
            i2x2 = nx2
            if ( x2wrap ) then
*               both wrap - no cutout region
                i1x1 = 1
                i2x1 = nx1
            else                ! x2wrap
*               cut out the region corresponding with field2
                call getlonwindow(xx2(1),xx2(nx2),i1x1,i2x1,xx1,nx1,
     +                d1,d2,lwrite)
                if ( d1.gt.1e33 .or. d2.gt.1e33 ) then
                    write(0,*) 'interpu: cannot find lon ',lon1,lon2
                    write(*,*) 'interpu: cannot find lon ',lon1,lon2
                    call abort
                endif
                if ( i1x1.gt.1   .and. abs(xx1(i1x1)-xx2(1)).gt.1e-3 )
     +                i1x1 = i1x1-1
                if ( i2x1.lt.nx1 .and. abs(xx1(1+mod(i2x1-1,nx1))
     +                -xx2(nx2)).gt.1e-3 ) i2x1 = i2x1+1
            endif
        else                    ! x1wrap
            if ( x2wrap ) then
*               cut out the region corresponding with field1
                i1x1 = 1
                i2x1 = nx1
                lon1 = xx1(1)
                lon2 = xx1(nx1)
                call getlonwindow(xx1(1),xx1(nx1),i1x2,i2x2,xx2,nx2,
     +                d1,d2,lwrite)
                if ( d1.gt.1e33 .or. d2.gt.1e33 ) then
                    write(0,*) 'interpu: cannot find lon ',lon1,lon2
                    write(*,*) 'interpu: cannot find lon ',lon1,lon2
                    call abort
                endif
                if ( i1x2.gt.1   .and. abs(xx2(i1x2)-xx(1)).gt.1e-3 )
     +                i1x2 = i1x2-1
                if ( i2x2.lt.nx2 .and. abs(xx2(i2x1)-xx(nx1)).gt.1e-3 )
     +                i2x2 = i2x2+1
            else                ! x2wrap
*               cut out the region corresponding with the intersection
*               of field1,field2
                d1 = min(xx1(nx1),xx2(nx2)) - max(xx1(1),xx2(1))
                d2 = min(xx1(nx1)-360,xx2(nx2)) - max(xx1(1)-360,xx2(1))
                d3 = min(xx1(nx1)+360,xx2(nx2)) - max(xx1(1)+360,xx2(1))
                if ( d2.gt.d1 .and. d2.gt.d3 ) then
                    do i=1,nx1
                        xx1(i) = xx1(i) - 360
                    enddo
                elseif ( d3.gt.d1 .and. d3.gt. d2 ) then
                    do i=1,nx1
                        xx1(i) = xx1(i) + 360
                    enddo
                endif
                lon1 = max(xx1(1),xx2(1))
                lon2 = min(xx1(nx1),xx2(nx2))
                if ( lon1.gt.lon2 ) then
                    write(0,*) 'interpu: error: no overlap in longitude'
     +                    ,xx1(1),xx1(nx1),xx2(1),xx2(nx2)
                    write(*,*) 'interpu: error: no overlap in longitude'
     +                    ,xx1(1),xx1(nx1),xx2(1),xx2(nx2)
                    call abort
                endif
                call getlonwindow(lon1,lon2,i1x1,i2x1,xx1,nx1,d1,d2
     +               ,lwrite)
                if ( d1.gt.1e33 .or. d2.gt.1e33 ) then
                    write(0,*) 'interpu: error: cannot find lon ',lon1
     +                   ,lon2
                    write(*,*) 'interpu: error: cannot find lon ',lon1
     +                   ,lon2
                    call abort
                endif
                if ( i1x1.gt.1   .and. abs(xx1(i1x1)-lon1).gt.1e-3 )
     +                i1x1 = i1x1-1
                if ( i2x1.lt.nx1 .and. abs(xx1(i2x1)-lon2).gt.1e-3 )
     +                i2x1 = i2x1+1
                call getlonwindow(lon1,lon2,i1x2,i2x2,xx2,nx2,d1,d2
     +               ,lwrite)
                if ( d1.gt.1e33 .or. d2.gt.1e33 ) then
                    write(0,*) 'interpu: cannot find lon ',lon1,lon2
                    call abort
                endif
                if ( i1x2.gt.1   .and. abs(xx2(i1x2)-lon1).gt.1e-3 )
     +                i1x2 = i1x2-1
                if ( i2x2.lt.nx2 .and. abs(xx2(i2x2)-lon2).gt.1e-3 )
     +                i2x2 = i2x2+1
            endif
        endif
        if ( lwrite ) then
            print '(a)',' interpu: cutout region x'
            print '(a,2i7,a,2i7)',' grid1: cut out ',i1x1,i2x1,' of ',
     +            1,nx1
            print '(a,2f7.2,a,2f7.2)',' grid1: cut out ',
     +            xx1(i1x1),xx1(1+mod(i2x1-1,nx1)),' of ',
     +            xx1(1),xx1(nx1)
            print '(a,2i7,a,2i7)',' grid2: cut out ',i1x2,i2x2,' of ',
     +            1,nx2
            print '(a,2f7.2,a,2f7.2)',' grid2: cut out ',
     +            xx2(i1x2),xx2(1+mod(i2x2-1,nx2)),' of ',
     +            xx2(1),xx2(nx2)
        endif
*
*       cut-offs latitude: simpler (no wrapping)
*
        lat1 = max(yy1(1),yy2(1))
        lat2 = min(yy1(ny1),yy2(ny2))
        if ( lat1.gt.lat2 ) then
            write(0,*) 'interpu: error: no overlap in latitude'
     +            ,yy1(1),yy1(ny1),yy2(1),yy2(ny2)
            write(*,*) 'interpu: error: no overlap in latitude'
     +            ,yy1(1),yy1(ny1),yy2(1),yy2(ny2)
            call abort
        endif
        call getlatwindow(lat1,lat2,i1y1,i2y1,yy1,ny1,d1,d2,lwrite)
        if ( i1y1.gt.1   .and. abs(yy1(i1y1)-lat1).gt.1e-3 )
     +        i1y1 = i1y1-1
        if ( i2y1.lt.ny1 .and. abs(yy1(i2y1)-lat2).gt.1e-3 )
     +        i2y1 = i2y1+1
        call getlatwindow(lat1,lat2,i1y2,i2y2,yy2,ny2,d1,d2,lwrite)
        if ( i1y2.gt.1   .and. abs(yy2(i1y2)-lat1).gt.1e-3 )
     +        i1y2 = i1y2-1
        if ( i2y2.lt.ny2 .and. abs(yy2(i2y2)-lat2).gt.1e-3 )
     +        i2y2 = i2y2+1
        if ( lwrite ) then
            print '(a)',' interpu: cutout region y'
            print '(a,2i7,a,2i7)',' grid1: cut out ',i1y1,i2y1,' of ',
     +            1,ny1
            print '(a,2f7.2,a,2f7.2)',' grid1: cut out ',
     +            yy1(i1y1),yy1(i2y1),' of ',yy1(1),yy1(ny1)
            print '(a,2i7,a,2i7)',' grid2: cut out ',i1y2,i2y2,' of ',
     +            1,ny2
            print '(a,2f7.2,a,2f7.2)',' grid2: cut out ',
     +            yy2(i1y2),yy2(i2y2),' of ',yy2(1),yy2(ny2)
        endif
*
*       minimum grid distances
*
        if ( x1wrap .and. i1x1.eq.1 .and. i2x1.eq.nx1 ) then
            dxmin1 = xx1(1) - xx1(nx1) + 360
        else
            dxmin1 = 360
        endif
        do i=i1x1,i2x1-1
            dxmin1 = min(dxmin1,xx1(1+mod(i,nx1))-xx1(1+mod(i-1,nx1)))
        enddo
        if ( x2wrap .and. i1x2.eq.1 .and. i2x2.eq.nx2 ) then
            dxmin2 = xx2(1) - xx2(nx2) + 360
        else
            dxmin2 = 360
        endif
        do i=i1x2,i2x2-1
            dxmin2 = min(dxmin2,xx2(1+mod(i,nx2))-xx2(1+mod(i-1,nx2)))
        enddo
        dymin1 = 180
        do i=i1y1,i2y1-1
            dymin1 = min(dymin1,yy1(i+1)-yy1(i))
        enddo
        dymin2 = 180
        do i=i1y2,i2y2-1
            dymin2 = min(dymin2,yy2(i+1)-yy2(i))
        enddo
        call chooseint(lintx1,lintx2,dxmin1,dxmin2,i1x1,i2x1,i1x2,i2x2
     +        ,xx1,xx2,intertype,lwrite)
        call chooseint(linty1,linty2,dymin1,dymin2,i1y1,i2y1,i1y2,i2y2
     +        ,yy1,yy2,intertype,lwrite)
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
*       
*       get new grid, compute interpolation coefficients
        if ( lintx1 ) then
            if ( lintx2 ) then
                write(0,*) 'interpu: cannot interpolate both x fields'
                call abort
            else
                call setupint(xx,ax,ix,nx,xx1,nx1,x1wrap,xx2,nx2,i1x2
     +                ,i2x2,x2wrap,lwrite)
            endif
        else
            if ( lintx2 ) then
                call setupint(xx,ax,ix,nx,xx2,nx2,x2wrap,xx1,nx1,i1x1
     +                ,i2x1,x1wrap,lwrite)
            else
*               both grids are already equal, copy grid info
                call setupint(xx,ax,ix,nx,xx2,nx2,x2wrap,xx1,nx1,i1x1
     +                ,i2x1,x1wrap,lwrite)
            endif
        endif
        if ( linty1 ) then
            if ( linty2 ) then
                write(0,*) 'interpu: cannot interpolate both y fields'
                call abort
            else
                call setupint(yy,ay,iy,ny,yy1,ny1,.FALSE.,yy2,ny2,i1y2
     +                ,i2y2,.FALSE.,lwrite)
            endif
        else
            if ( linty2 ) then
                call setupint(yy,ay,iy,ny,yy2,ny2,.FALSE.,yy1,ny1,i1y1
     +                ,i2y1,.FALSE.,lwrite)
            else
*               both grids are already equal, copy grid info
                call setupint(yy,ay,iy,ny,yy2,ny2,.FALSE.,yy1,ny1,i1y1
     +                ,i2y1,.FALSE.,lwrite)
            endif
        endif
*
*       anything to do?
*
        if (  .not.lintx1 .and. .not.lintx2 .and.
     +        .not.linty1 .and. .not.linty2 .and.
     +        (x1rev.eqv.x2rev) .and. (y1rev.eqv.y2rev) ) then
            if ( lwrite ) print *,'no interpolation necessary'
            return
        endif
*
*       perform the interpolations
*
        do yr=yr11,yr21
            do mo=1,nperyear
*               interpolate/copy from field1 to f
                if ( lintx1 ) then
                    if ( lwrite ) print *,'interpolating f1->f x ',mo,yr
                   call doxint(f,nxmax,nymax,field1(1,1,mo,yr),nxf,nyf
     +                    ,nx1,ny1,ax,ix,nx,x1rev)
                elseif ( linty1 .or. i1x1.ne.1 .or. i1y1.ne.1 .or. x1rev
     +                    .or.y1rev ) then
                    if ( lwrite ) print *,'shifting      f1->f x ',mo,yr
                    call doxnint(f,nxmax,nymax,field1(1,1,mo,yr),nxf,nyf
     +                    ,nx1,ny1,i1x1,nx,x1rev)
                endif
                if ( lwrite .and. yr.eq.yr11 ) then
                    print *,'f(1,1,',mo,yr,') = ',f(1,1)
                endif
*               interpolate/copy from f back to field1
                if ( linty1 ) then
                    if ( lwrite ) print *,'interpolating f->f1 y ',mo,yr
                    call doyint(field1(1,1,mo,yr),nxf,nyf,f,nxmax,nymax
     +                    ,nx1,ny1,ay,iy,ny,nx,y1rev)
                elseif ( lintx1 .or. i1x1.ne.1 .or. i1y1.ne.1 .or. x1rev
     +                    .or. y1rev ) then
                    if ( lwrite ) print *,'shifting      f->f1 y ',mo,yr
                    call doynint(field1(1,1,mo,yr),nxf,nyf,f,nxmax,nymax
     +                    ,nx1,ny1,i1y1,nx,ny,y1rev)
                endif
                if ( lwrite .and. yr.eq.yr11 ) then
                    print *,'field1(1,1,',mo,yr,') = ',field1(1,1,mo,yr)
                endif
            enddo               ! mo
        enddo                   ! yr1
        do yr=yr12,yr22
            do mo=1,nperyear        
*               interpolate/copy from field2 to f
                if ( lintx2 ) then
                    if ( lwrite ) print *,'interpolating f2->f x ',mo,yr
                    call doxint(f,nxmax,nymax,field2(1,1,mo,yr),nxf,nyf
     +                    ,nx2,ny2,ax,ix,nx,x2rev)
                elseif ( linty2 .or. i1x2.ne.1 .or. i1y2.ne.1 .or. x2rev
     +                    .or.y2rev ) then
                    if ( lwrite ) print *,'shifting      f2->f x ',mo,yr
                    call doxnint(f,nxmax,nymax,field2(1,1,mo,yr),nxf,nyf
     +                    ,nx2,ny2,i1x2,nx,x2rev)
                endif
*               interpolate/copy from f back to field2
                if ( linty2 ) then
                    if ( lwrite ) print *,'interpolating f->f2 y ',mo,yr
                    call doyint(field2(1,1,mo,yr),nxf,nyf,f,nxmax,nymax
     +                    ,nx2,ny2,ay,iy,ny,nx,y2rev)
                elseif ( lintx2 .or. i1x2.ne.1 .or. i1y2.ne.1 .or. x2rev
     +                    .or.y2rev ) then
                    if ( lwrite ) print *,'shifting      f->f2 y ',mo,yr
                    call doynint(field2(1,1,mo,yr),nxf,nyf,f,nxmax,nymax
     +                    ,nx2,ny2,i1y2,nx,ny,y2rev)
                endif
            enddo               ! mo
        enddo                   ! yr2
*
*       that's it
*
        end

        subroutine chooseint(lintx1,lintx2,dxmin1,dxmin2,
     +        i1x1,i2x1,i1x2,i2x2,xx1,xx2,intertype,lwrite)
*
*       choose the grid with the highest resolution (intertype=0),
*       lowest resolution (intertype=-1), grid 1 (1), grid (2).
*       in case of a tie choose the grid with the largest gaps (least points)
*       also check whether the endpoints agree
*       in case of a tie choose field 2 if the grid points do not agree
*       otherwise interpolation is not necessary
*       
        implicit none
        logical lintx1,lintx2
        integer i1x1,i2x1,i1x2,i2x2,intertype
        real dxmin1,dxmin2,xx1(i2x1),xx2(i2x2)
        logical lwrite
        integer i
*
        lintx1 = .FALSE.
        lintx2 = .FALSE.
        if ( intertype.eq.1 ) then
            if ( i1x2.ne.i1x1 .or. i2x2.ne.i2x1 ) then
                if ( lwrite ) print *,'chooseint: grid2 differs ',
     +                i1x2,i2x2,i1x1,i2x1
                lintx2 = .TRUE.
            else
                do i=i1x1,i2x1
                    if ( abs(xx1(i)-xx2(i)).gt.1e-3 ) then
                        if ( lwrite ) print *,'chooseint: grid2(',i
     +                        ,') differs ',xx1(i),xx2(i)
                        lintx2 = .TRUE.
                        goto 10
                    endif
                enddo
   10           continue
            endif
        elseif ( intertype.eq.2 ) then
            if ( i1x2.ne.i1x1 .or. i2x2.ne.i2x1 ) then
                if ( lwrite ) print *,'chooseint: grid1 differs ',
     +                i1x2,i2x2,i1x1,i2x1
                lintx1 = .TRUE.
            else
                do i=i1x1,i2x1
                    if ( abs(xx1(i)-xx2(i)).gt.1e-3 ) then
                        if ( lwrite ) print *,'chooseint: grid1(',i
     +                        ,') differs ',xx1(i),xx2(i)
                        lintx1 = .TRUE.
                        goto 20
                    endif
                enddo
   20           continue
            endif
        elseif ( dxmin1.lt.dxmin2 ) then
            if ( lwrite ) print *,'chooseint: grid 1 is finer: ',dxmin1,
     +           dxmin2
            if ( intertype.eq.0 ) then
                lintx2 = .TRUE.
            else
                lintx1 = .TRUE.
            endif
        elseif ( dxmin2.lt.dxmin1 ) then
            if ( lwrite ) print *,'chooseint: grid 2 is finer: ',dxmin2,
     +           dxmin1
            if ( intertype.eq.0 ) then
                lintx1 = .TRUE.
            else
                lintx2 = .TRUE.
            endif
        elseif ( i2x1-i1x1.lt.i2x2-i1x2 ) then
            if ( lwrite ) print *,'chooseint: grid 1 has fewer points'
            lintx2 = .TRUE.
        elseif ( i2x2-i1x2.lt.i2x1-i1x1 ) then
            if ( lwrite ) print *,'chooseint: grid 2 has fewer points'
            lintx1 = .TRUE.
        elseif ( xx1(i2x1)-xx1(i1x1).lt.xx2(i2x2)-xx2(i1x2) ) then
            if ( lwrite ) print *,'chooseint: grid 1 is smaller'
            lintx2 = .TRUE.
        elseif ( xx2(i2x2)-xx2(i1x2).lt.xx1(i2x1)-xx1(i1x1) ) then
            if ( lwrite ) print *,'chooseint: grid 2 is smaller'
            lintx1 = .TRUE.
        elseif ( i1x1.ne.i1x2 ) then
            if ( lwrite ) print *,'chooseint: grid 2 is shifted'
            lintx2 = .TRUE.     ! but probably a trivial shift
        else
            do i=i1x1,i2x1
                if ( abs(xx1(i)-xx2(i)).gt.1e-3 ) then
                    if ( lwrite ) print *,'chooseint: grid2(',i
     +                    ,') differs ',xx1(i),xx2(i)
                    lintx2 = .TRUE.
                    goto 100
                endif
            enddo
  100       continue
        endif
        end

        subroutine setupint(xx,ax,ix,nx,xx1,nx1,x1wrap,xx2,nx2,i1x2,i2x2
     +        ,x2wrap,lwrite)
*
*       For interpolation of field1 to the grid of field2 
*       sets up the interpolation coefficients ax, the indices of the
*       points between which the interpolation takes place ix, and fill
*       xx and nx.
*
        implicit none
        integer nx,nx1,nx2,i1x2,i2x2,ix(i2x2)
        real ax(i2x2),xx(i2x2),xx1(nx1),xx2(nx2)
        logical x1wrap,x2wrap,lwrite
        integer i,j,i1,i2
        real x1,x2
*
*       define grid
*
        nx = 1+i2x2-i1x2
        do i=1,nx
            j = i+i1x2-1
            if ( x2wrap .and. j.lt.1 ) then
                xx(i) = xx2(j+nx2) - 360
            elseif ( j.le.nx2 ) then
                xx(i) = xx2(j)
            elseif ( x2wrap .and. j.le.2*nx2 ) then
                xx(i) = xx2(j-nx2) + 360
            else
                write(0,*) 'error: more than 720o'
                call abort
            endif
        enddo
*
*       yet more wrapping misery
*
        if ( x2wrap ) then
            if ( xx(1).gt.xx1(nx1) ) then
                if ( lwrite ) print *,'shifting x(i) 360o down'
                do i=1,nx
                    xx(i) = xx(i) - 360
                enddo
            elseif ( xx(nx).lt.xx1(1) ) then
                if ( lwrite ) print *,'shifting x(i) 360o up'
                do i=1,nx
                    xx(i) = xx(i) + 360
                enddo
            endif
        endif
        if ( xx(1).gt.xx1(nx1) .or. xx(nx).lt.xx1(1) ) then
            write(0,*) 'setupint: error: no overlapping coordinates'
            write(0,*) (xx(i),i=1,nx)
            write(0,*) (xx1(i),i=1,nx1)
        endif
*
*       special case
*
        if ( nx.eq.1 .and. nx1.eq.1 ) then
            ix(1) = 1
            ax(1) = 1
            if ( lwrite ) print '(a,i5)'
     +            ,'setupint: interpolation set-up special case'
            return
        endif
*
*       set interpolation coefficients
*
        do i=1,nx
*           search for largest neighbours under xx(i) in grid1
            ix(i) = 1
            if ( x1wrap ) then
                do j=1-nx1,0
                    if ( xx1(j+nx1)-360.lt.xx(i) ) ix(i) = j
                enddo
            endif
            do j=1,nx1-1
                if ( xx1(j).lt.xx(i) ) ix(i) = j
            enddo
            if ( x1wrap ) then
                if ( xx1(nx1).lt.xx(i) ) ix(i) = nx1
                do j=nx1+1,2*nx1
                    if ( xx1(j-nx1)+360.lt.xx(i) ) ix(i) = j
                enddo
            endif
*           wrapping problems
            i1 = ix(i)
            i2 = ix(i) + 1
            if ( x1wrap ) then
                if ( i1.le.0 )   i1 = i1 + nx1
                if ( i1.gt.nx1 ) i1 = i1 - nx1
                if ( i2.le.0 )   i2 = i2 + nx1
                if ( i2.gt.nx1 ) i2 = i2 - nx1
            endif
            if (  i1.le.0 .or. i1.gt.nx1 .or.
     +            i2.le.0 .or. i2.gt.nx1 ) then
                write(0,*) 'setupint: error: out of boundaries',i1,i2
     +                ,nx1
                write(*,*) 'setupint: error: out of boundaries',i1,i2
     +                ,nx1
                call abort
            endif
            x1 = xx1(i2) - xx(i)
            if ( x1.lt.-180 ) x1 = x1 + 360
            if ( x1.gt.+180 ) x1 = x1 - 360
            x2 = xx1(i2) - xx1(i1)
            if ( x2.lt.-180 ) x2 = x2 + 360
            if ( x2.gt.+180 ) x2 = x2 - 360
            ax(i) = x1/x2
        enddo
*
*       debug output
*
        if ( lwrite ) then
            print '(a,i5)','setupint: interpolation set-up ',nx
            do i=1,nx
                i1 = ix(i)
                if ( i1.le.0 )   i1 = i1 + nx1
                if ( i1.gt.nx1 ) i1 = i1 - nx1
                i2 = ix(i) + 1
                if ( i2.le.0 )   i2 = i2 + nx1
                if ( i2.gt.nx1 ) i2 = i2 - nx1
                print '(i4,f8.2,i4,f8.2,f8.4,i4,f8.2,f8.4)',i,xx(i),
     +                i1,xx1(i1),ax(i),i2,xx1(i2),1-ax(i)
            enddo
        endif
        end

        subroutine doxint(f,nxmax,nymax,f1,nxf,nyf,nx1,ny1,ax,ix,nx,
     +        xrev)
*
*       perform the interpolation in x from f1 into f
*
        implicit none
        integer nxmax,nymax,nxfmnyf,nxf,nyf,nx1,ny1,nx
        integer ix(nx)
        real f(nxmax,nymax),f1(nxf,nyf),ax(nx)
        logical xrev
        integer i,j,i1,i2
        real y
*
        do j=1,ny1
            if ( xrev ) then
                do i=1,nx1/2
                    y = f1(i,j)
                    f1(i,j) = f1(nx1-i+1,j)
                    f1(nx1-i+1,j) = y
                enddo
            endif
            do i=1,nx
                i1 = ix(i)
                if ( i1.lt.1 )   i1 = i1 + nx1
                if ( i1.gt.nx1 ) i1 = i1 - nx1
                if ( ax(i).eq.1 ) then
                    f(i,j) = f1(i1,j)
                else
                    i2 = ix(i) + 1
                    if ( i2.lt.1 )   i2 = i2 + nx1
                    if ( i2.gt.nx1 ) i2 = i2 - nx1
                    if ( ax(i).eq.0 ) then
                        f(i,j) = f1(i2,j)
                    elseif ( f1(i1,j).lt.1e33 ) then
                        if ( f1(i2,j).lt.1e33 ) then
                            f(i,j) = ax(i)*f1(i1,j) + (1-ax(i))*f1(i2,j)
                        elseif ( ax(i).gt.0.75 ) then
                            f(i,j) = f1(i1,j)
                        else
                            f(i,j) = 3e33
                        endif
                    else
                        if ( f1(i2,j).lt.1e33 .and. 1-ax(i).gt.0.75 )
     +                        then
                            f(i,j) = f1(i2,j)
                        else
                            f(i,j) = 3e33
                        endif
                    endif
                    if ( .false. ) then
                        print *,'ax(i) = ',ax(i)
                        print *,'f1(',i1,j,') = ',f1(i1,j)
                        print *,'f1(',i2,j,') = ',f1(i2,j)
                        print *,'f(',i,j,') = ',f(i,j)
                    endif
                endif
            enddo
        enddo
        end
        subroutine doxnint(f,nxmax,nymax,f1,nxf,nyf,nx1,ny1,i1x1,nx,
     +        xrev)
*       
*       special case: just a shift plus maybe a reversal
*       
        implicit none
        integer nxmax,nymax,nxf,nyf,nx1,ny1,i1x1,nx
        real f(nxmax,nymax),f1(nxf,nyf)
        logical xrev
        integer i,j,i1
        real y
*
        do j=1,ny1
            if ( xrev ) then
                do i=1,nx1/2
                    y = f1(i,j)
                    f1(i,j) = f1(nx1-i+1,j)
                    f1(nx1-i+1,j) = y
                enddo
            endif
            do i=1,nx
                i1 = i+i1x1-1
                if ( i1.lt.1 )   i1 = i1 + nx1
                if ( i1.gt.nx1 ) i1 = i1 - nx1
                f(i,j) = f1(i1,j)
            enddo
        enddo
        end

        subroutine doyint(f1,nxf,nyf,f,nxmax,nymax,nx1,ny1,ay,iy,ny,nx
     +        ,yrev)
*
*       perform the interpolation in y from f into f1
*
        implicit none
        integer nxf,nyf,nx1,ny1,nxmax,nymax,ny,nx
        integer iy(ny)
        real f1(nxf,nyf),f(nxmax,nymax),ay(ny)
        logical yrev
        integer i,j,j1
        real y
*        
        if ( yrev ) then
            do i=1,nx
                do j=1,ny1/2
                    y = f(i,j)
                    f(i,j) = f(i,ny1-j+1)
                    f(i,ny1-j+1) = y
                enddo
            enddo
        endif
        do j=1,ny
            j1 = iy(j)
            if ( ay(j).eq.1 ) then
                do i=1,nx
                    f1(i,j) = f(i,j1)
                enddo
            else
                do i=1,nx
                    if ( f(i,j1).lt.1e33 ) then
                        if ( f(i,j1+1).lt.1e33 ) then
                            f1(i,j) = ay(j)*f(i,j1) +(1-ay(j))*f(i,j1+1)
                        elseif ( ay(j).gt.0.75 ) then
                            f1(i,j) = f(i,j1)
                        else
                            f1(i,j) = 3e33
                        endif
                    else
                        if ( f(i,j1+1).lt.1e33 .and. 1-ay(j).gt.0.75 )
     +                        then
                            f1(i,j) = f(i,j1+1)
                        else
                            f1(i,j) = 3e33
                        endif
                    endif
                enddo
            endif
        enddo
        end
        subroutine doynint(f1,nxf,nyf,f,nxmax,nymax,nx1,ny1,i1y1,nx,ny,
     +        yrev)
*       
*       special case: just a shift and maybe reverse
*
        implicit none
        integer nxf,nyf,nx1,ny1,nxmax,nymax,i1y1,nx,ny
        real f1(nxf,nyf),f(nxmax,nymax)
        logical yrev
        integer i,j
        real y
*
        if ( yrev ) then
            do i=1,nx
                do j=1,ny1/2
                    y = f(i,j)
                    f(i,j) = f(i,ny1-j+1)
                    f(i,ny1-j+1) = y
                enddo
            enddo
        endif
        do j=1,ny
            do i=1,nx
                f1(i,j) = f(i,j+i1y1-1)
            enddo
        enddo
        end

