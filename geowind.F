        program geowind
!
!       given a sea-level pressure field, compute geostrophic wind components
!
!       (c) Geert Jan van Oldenborgh, KNMI, Nov 2007
!
        implicit none
#include "netcdf.inc"
#include "recfac.h"
        integer nxmax,nymax,nzmax,nvmax,yrbeg,yrend,npermax,ndatmax
     +       ,indxmax
        parameter(nxmax=360,nymax=180,nzmax=1,nvmax=1,
     +       yrbeg=1800,yrend=2300,npermax=12,ndatmax=1000,indxmax=6)
        integer i,j,k,l,n,mo,yr,i1,i2,j1,j2
     +       ,status,irec,jul1,jul2
        integer ncid,nx,ny,nz,nt,nperyear,firstyr,firstmo,endian
     +       ,nvars,ivars(2,2*indxmax+2),jvars(6,nvmax),lastyr
        real xx(nxmax),yy(nymax),zz(nzmax),undef
        real lon,lat,xsize,ysize,f,deg2m,mb2pa,deg2rad,rho,dx,dy
        real,allocatable :: slp(:,:,:,:),ugeo(:,:,:,:,:)
        logical lwrite,xrev,yrev,xwrap
        character slpfile*255,datfile*255,outfile(3)*255,title*255
        character vars(nvmax)*10,lvars(nvmax)*80,units(nvmax)*20
        character ovars(2*indxmax+2)*10,olvars(2*indxmax+2)*80,
     +       ounits(2*indxmax+2)*20,otitle*255,octlfile*255,odatfile*255
        character string*80
!
        integer iargc,julday,leap
!
!       usage
!
        if ( iargc().lt.4 ) then
            print *,'usage: ugeo slp.{ctl|nc} '//
     +           'gwest.ctl gsouth.ctl vort.ctl'
            call abort
        endif
!
!       init
!
        lwrite = .false.
        call getarg(5,string)
        if ( string(1:5).eq.'debug' .or.
     +       string(1:5).eq.'lwrite' ) lwrite = .true.
        xsize = 20              ! degrees
        ysize = 10              ! degrees
        deg2m = 1e7/90
        deg2rad = atan(1.)/45
        mb2pa = 100
        rho = 1
        f = 8*atan(1.)/(24*60*60)
        call getarg(1,slpfile)
        call getarg(2,outfile(1))
        call getarg(3,outfile(2))
        call getarg(4,outfile(3))
!
!       read field
!
        print *,'reading field'
        status = nf_open(slpfile,nf_nowrite,ncid)
        if ( status.ne.nf_noerr ) then
            call parsectl(slpfile,datfile,nxmax,nx,xx,nymax,ny,yy
     +           ,nzmax,nz,zz,nt,nperyear,firstyr,firstmo,undef
     +           ,endian,title,1,nvars,vars,ivars,lvars,units)
            ncid = -1
        else
            call parsenc(slpfile,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax
     +           ,nz,zz,nt,nperyear,firstyr,firstmo,undef,title,1
     +           ,nvars,vars,jvars,lvars,units)
        endif
        lastyr = firstyr + (nt+firstmo-2)/nperyear
        if ( lwrite) print *,'allocating fields(',nx,ny,nperyear,firstyr
     +       ,lastyr,')'
        allocate(slp(nx,ny,nperyear,firstyr:lastyr))

        if ( ncid.eq.-1 ) then
            call readdatfile(datfile,slp,nx,ny,nx,ny,
     +           nperyear,firstyr,lastyr,firstyr,firstmo,nt,undef,endian
     +           ,lwrite,firstyr,lastyr,1,1)
        else
            call readncfile(ncid,slp,nx,ny,nx,ny,
     +           nperyear,firstyr,lastyr,firstyr,firstmo,nt,undef,
     +           lwrite,firstyr,lastyr,jvars)
        endif
        call makestandardfield(slp,nx,ny,1,nperyear,firstyr,lastyr
     +       ,nx,ny,1,nperyear,firstyr,lastyr,vars(1),units(1),lwrite)
        call getxyprop(xx,nx,yy,ny,xrev,yrev,xwrap)
        if ( lwrite ) print *,'xwrap = ',xwrap
!
!       compute geostrophic wind components
!
        print *,'computing geostrophic wind'
        if ( lwrite ) print *,'allocating ',
     +       nx,ny,nperyear,firstyr,lastyr,3
        allocate(ugeo(nx,ny,nperyear,firstyr:lastyr,3))
        do j=1,ny
            call keepalive(j,ny)
            lat = yy(j)
            if ( yy(j)+ysize/2.gt.max(yy(1),yy(ny)) .or.
     +           yy(j)-ysize/2.lt.min(yy(1),yy(ny)) ) then
                if ( lwrite ) print *,'skipping lat point ',j,lat
                ugeo(1:nx,j,1:nperyear,firstyr:lastyr,1:3) = 3e33
                cycle
            endif
            do i=1,nx
                lon = xx(i)
                if ( .not.xwrap .and. (
     +               xx(i)+xsize/2.gt.max(xx(1),xx(nx)) .or.
     +               xx(i)-xsize/2.lt.min(xx(1),xx(nx)) ) ) then
                    if ( lwrite ) print *,'skipping lon point ',i,lon
     +                   ,xwrap
                    ugeo(i,j,1:nperyear,firstyr:lastyr,1:3) = 3e33
                    cycle
                endif
                if ( lwrite ) print *,'calling getlatlonwindow'
                call getlatlonwindow(lat-ysize/2,lat+ysize/2,lon-xsize/2
     +               ,lon+xsize/2,xx,nx,xwrap,1,yy,ny,1,i1,i2,j1,j2
     +               ,.false.)
                if ( yy(j1)*yy(j2).le.0 .or. j1.eq.j2 .or.
     +               .not.xwrap .and. i1.eq.i2 ) then
!                   does not work at edges or across the equator
                    if ( lwrite ) print *,'skipping lat point ',j,lat
                    ugeo(i,j,1:nperyear,firstyr:lastyr,1:3) = 3e33
                    cycle
                endif
                if ( xwrap ) then
                    if ( i1.lt.1 ) i1 = i1 + nx
                    if ( i2.lt.1 ) i2 = i2 + nx
                    if ( i1.gt.nx ) i1 = i1 - nx
                    if ( i2.gt.nx ) i2 = i2 - nx
                    if ( i2.gt.nx ) i2 = i2 - nx
                endif
                if ( i1.eq.i2 ) then
                    write(0,*) 'vsm: error: i1=i2: ',i1,i2,j1,j2,lon,lat
                    call abort
                endif
                dx = xx(i2)-xx(i1)
                if ( xwrap ) then
                    if ( abs(dx-360).lt.abs(dx) ) then
                        dx = dx - 360
                    elseif ( abs(dx+360).lt.abs(dx) ) then
                        dx = dx + 360
                    endif
                endif
                dx = dx*deg2m*cos(yy(j)*deg2rad)
                dy = (yy(j2) - yy(j1))*deg2m
                do yr=firstyr,lastyr
                    do mo=1,nperyear
                        if ( slp(i1,j1,mo,yr).gt.1e33 .or.
     +                       slp(i1,j2,mo,yr).gt.1e33 .or.
     +                       slp(i2,j1,mo,yr).gt.1e33 .or.
     +                       slp(i2,j2,mo,yr).gt.1e33 .or.
     +                       slp(i,j,mo,yr).gt.1e33 ) then
                            if ( lwrite ) print *,'invalid SLP ',i1,i,i2
     +                           ,j1,j,j2,mo,yr
                            ugeo(i,j,mo,yr,1:3) = 3e33
                            cycle
                        endif
                        ugeo(i,j,mo,yr,1) = mb2pa/rho*(
     +                       (slp(i1,j1,mo,yr) + slp(i2,j1,mo,yr))
     +                       - (slp(i1,j2,mo,yr) + slp(i2,j2,mo,yr)))
     +                       /2/dy/(2*f*sin(yy(j)*deg2rad))
                        ugeo(i,j,mo,yr,2) = -mb2pa/rho*(
     +                       (slp(i1,j1,mo,yr) + slp(i1,j2,mo,yr))
     +                       - (slp(i2,j1,mo,yr) + slp(i2,j2,mo,yr)))
     +                       /2/dx/(2*f*sin(yy(j)*deg2rad))
                        ugeo(i,j,mo,yr,3) = 
     +                       ( slp(i1,j1,mo,yr) + slp(i2,j1,mo,yr)
     +                       + slp(i1,j2,mo,yr) + slp(i2,j2,mo,yr) )/4
     +                       - slp(i,j,mo,yr)
                        if ( .false. .and. lwrite ) print *,'ugeo(',i,j
     +                       ,mo,yr,'1-3) = ',(ugeo(i,j,mo,yr,k),k=1,3)
                    enddo       ! mo
                enddo           ! yr
            enddo
        enddo
!
!       write output: geostrophic wind
!
        print *,'writing geo.ctl/grd'
        otitle = 'geostrophic wind components of '//trim(slpfile)
        undef = 3e33
        do k=1,3
            i = index(outfile(k),'.ctl') - 1
            if ( i.le.0 ) i = 1 + len_trim(outfile(i))
            odatfile = outfile(k)(1:i)//'.grd'
            if ( k.eq.1 ) then
                ovars(1) = 'ugeo'
                olvars(1) = 'zonal geostrophic wind'
                ounits(1) = 'm/s'
            elseif ( k.eq.2 ) then
                ovars(1) = 'vgeo'
                olvars(1) = 'meridional geostrophic wind'
                ounits(1) = 'm/s'
            else
                ovars(1) = 'vort'
                olvars(1) = 'geostrophic vorticity'
                ounits(1) = units(1)
            endif
            ivars(1,1) = 0
            ivars(2,1) = 99
            if ( 366*(nperyear/366).ne.nperyear ) then
                nt = (lastyr-firstyr+1)*nperyear
            else
                jul1 = julday( 1, 1,firstyr)
                jul2 = julday(12,31,lastyr)
                nt = (nperyear/366)*(jul2 - jul1 + 1)
            endif
            call writectl(outfile(k),odatfile,nx,xx,ny,yy,nz,zz,
     +           nt,nperyear,firstyr,1,undef,otitle,1,ovars,ivars,olvars
     +           ,ounits)
            open(1,file=odatfile,form='unformatted',access='direct',
     +           recl=recfa4*nx*ny)
            irec = 0
            do yr=firstyr,lastyr
                do mo=1,nperyear
                    if ( nperyear.ne.366 .or.
     +                   mo.ne.60 .or. leap(yr).eq.2 ) then
                        irec = irec + 1
                        write(1,rec=irec)
     +                       ((ugeo(i,j,mo,yr,k),i=1,nx),j=1,ny)
                    endif
                enddo
            enddo
            if ( irec.ne.nt ) then
                write(0,*) 'geowind: error in daily data ',nt,irec
            endif
            close(1)
        enddo
!
!       finito
!
        end
