        program extractfield
*
*       program to compute a netcdf file containing observations
*       for all grid points in a specified region, to be processed
*       further by the R routines
*
*       (c) Geert Jan van Oldenborgh, 31-may-2005, KNMI
*
*       This file may be copied and modified  freely as long as the
*       above copyright notice is kept intact
*
        implicit none
#include "params.h"
#include "recfac.h"
#include "netcdf.inc"
#include "getopts.inc"
        integer ncid,nx,ny,nz,nt,firstyr,f,firstmo,nvars
     +       ,jvars(6,nvmax),endian,status,nperyear
     +       ,lastyr,mens1,mens,iens,nn,ntmax,month,yr1s,yr2s,
     +       yrstart,yrstop
        integer, allocatable :: itimeaxis(:)
        integer i,j,j1,j2,k,n,jx,jy,mo,m,yr,y,nrec,ivars(2,6)
     +       ,ldir,iout,x1,x2,y1,y2,mobegin,ntvarid,ndup(npermax)
     +       ,validens(0:nensmax)
        real xx(nxmax),yy(nymax),zz(nzmax),
     +       wx(nxmax),wy(nymax),u
        real,allocatable :: field(:,:,:,:,:),fcst(:,:,:)
     +       ,newfcst(:,:,:,:,:)
        logical lexist,xrev,xwrap,yrev,lvalid
        character infile*256,datfile*256,title*256,vars(1)*20,
     +       lvars(1)*80,units(1)*20,file*256,tmpunits*20
        character line*80
        integer iargc,llen
*
*       check arguments
*
        lwrite = .FALSE.
        n = iargc()
        if ( n.lt.3 ) then
            print *,'usage: extractfield '//
     +           'field.[ctl|nc] '//
     +           ' [sum|ave|max|min|sel n] '//
     +           '[begin yr] [end yr] plot file.nc'
            stop
        endif
        call killfile(infile,datfile,file,0)
        call getarg(1,infile)
        mens1 = 0
        if ( index(infile,'%').gt.0 .or. index(infile,'++').gt.0 ) then
            ensemble = .true.
            call filloutens(infile,0)
            if ( lwrite ) print *,'check whether ',trim(infile)
     +           ,' exists'
            inquire(file=infile,exist=lexist)
            if ( .not.lexist ) then
                if ( lwrite ) print *,'no, start at 1'
                call getarg(1,infile)
                mens1 = 1
                call filloutens(infile,1)
            endif
            do mens=mens1+1,nensmax
                call getarg(1,file)
                if ( lwrite ) print *,'check whether ',trim(infile)
     +               ,' exists'
                call filloutens(file,mens)
                inquire(file=file,exist=lexist)
                if ( .not.lexist ) exit
            enddo
            mens = mens - 1
            print '(a,i2,a,i2)','# found ensemble members ',mens1,' to '
     +           ,mens
        else
            ensemble = .false.
            mens1 = 0
            mens = 0
        endif
        if ( lwrite ) print *,'extractfield: nf_opening file '
     +        ,infile(1:llen(infile))
        status = nf_open(infile,nf_nowrite,ncid)
        if ( status.ne.nf_noerr ) then
            call parsectl(infile,datfile,nxmax,nx,xx,nymax,ny,yy
     +           ,nzmax,nz,zz,nt,nperyear,firstyr,firstmo,u,endian
     +           ,title,1,nvars,vars,ivars,lvars,units)
            ncid = -1
        else
            datfile = infile
            call parsenc(infile,ncid,nxmax,nx,xx,nymax,ny,yy
     +            ,nzmax,nz,zz,nt,nperyear,firstyr,firstmo,u,title,1
     +            ,nvars,vars,jvars,lvars,units)
        endif
*
        lastyr = firstyr + (nt+firstmo-2)/nperyear
        call getopts(2,iargc(),nperyear,firstyr,lastyr,.true.,
     +       mens1,mens)
        if ( lag1.ne.0 .or. lag2.ne.0 ) print *,'ignoring lags'
        if ( .not.plot ) write(0,*)
     +       'extractfield: specify output file with ''plot'''
        if ( lks ) write(0,*) 'extractfield: K-S not supported'
        if ( mdiff.ne.0 .or. mdiff2.ne.0 ) write(0,*) 
     +      'extractfield: monthly anomalies not supported'
        if ( lconting ) write(0,*) 'extractfield: contingency '//
     +        'tables not supported'
        do i=1,indxuse
            if ( lincl(i) ) write(0,*) 'extractfield: what do ',
     +            'you mean with ',strindx(i),'?'
        enddo
        yr1 = max(yr1,firstyr)
        yr2 = min(yr2,lastyr)
        f = firstyr
        firstyr = yr1
        lastyr  = yr2
        if ( lwrite ) print *,'yr1,yr2 = ',yr1,yr2
        if ( lwrite ) print *,'nens1,nens2 = ',nens1,nens2
*
*       pgf90 on RHEL happily allocates >2GB :-(
        if ( nx*ny*nperyear*(lastyr-firstyr+1)*(nens2+1).gt.2**29 )
     +       then
            write(0,*) 'error: array too large, reduce resolution, '
     +           ,'number of years or number of ensemble members ',nx
     +           *ny*nperyear*(lastyr-firstyr+1)*(nens2+1)
            write(*,*) 'error: array too large, reduce resolution, '
     +           ,'number of years or number of ensemble members ',nx
     +           *ny*nperyear*(lastyr-firstyr+1)*(nens2+1)
            call abort
        endif
*
        allocate(field(nx,ny,nperyear,firstyr:lastyr,0:nens2))
        allocate(fcst(nperyear,firstyr:lastyr,0:nens2))
*
*       read fields
*
        if ( ncid.eq.-1 ) then
            if ( ensemble ) then
                do iens=nens1,nens2
                    call getarg(1,infile)
                    call filloutens(infile,iens)
                    call parsectl(infile,datfile,nxmax,nx,xx,nymax
     +                   ,ny,yy,nzmax,nz,zz,nt,nperyear,f
     +                   ,firstmo,u,endian,title,1,nvars,vars
     +                   ,ivars,lvars,units)
                    call keepalive(iens,nens2-nens1+2)
                    call readdatfile(datfile
     +                   ,field(1,1,1,firstyr,iens)
     +                   ,nx,ny,nx,ny,nperyear,firstyr,lastyr
     +                   ,f,firstmo,nt,u,endian,lwrite
     +                   ,max(firstyr,yr1),min(lastyr,yr2),1,1)
                 enddo
            else
                call keepalive(1,2)
                call readdatfile(datfile,field,nx,ny,nx,ny
     +               ,nperyear,firstyr,lastyr,f,firstmo,nt,u
     +               ,endian,lwrite,max(firstyr,yr1),min(lastyr,yr2)
     +               ,1,1)
            endif
        else
            if ( ensemble ) then
                do iens=nens1,nens2
                    call getarg(1,infile)
                    call filloutens(infile,iens)
                    call parsenc(infile,ncid,nxmax,nx,xx,nymax,ny
     +                   ,yy,nzmax,nz,zz,nt,nperyear,f,firstmo
     +                   ,u,title,1,nvars,vars,jvars,lvars,units)
                    call keepalive(iens,nens2-nens1+2)
                    call readncfile(ncid,field(1,1,1,firstyr,iens),nx
     +                   ,ny,nx,ny,nperyear,firstyr,lastyr,f
     +                   ,firstmo,nt,u,lwrite
     +                   ,max(firstyr,yr1),min(lastyr,yr2),jvars)
                enddo
            else
                call keepalive(1,2)
                call readncfile(ncid,field,nx,ny,nx,ny,nperyear
     +               ,firstyr,lastyr,firstyr,firstmo,nt,u,lwrite
     +               ,max(firstyr,yr1),min(lastyr,yr2),jvars)
            endif
        endif
*       convert to standard units
        do iens=nens1,nens2
            tmpunits = units(1) ! they are otherwise adjusted
            call makestandardfield(field(1,1,1,firstyr,iens),nx,ny,1
     +           ,nperyear,firstyr,lastyr,nx,ny,1,nperyear
     +           ,max(firstyr,yr1),min(lastyr,yr2),vars(1),tmpunits
     +           ,lwrite)
        enddo
        units(1) = tmpunits
        call keepalive(nens2-nens1+2,nens2-nens1+2)
*
*       cut out region of interest
*       well, slightly larger to avoid white region in plot later on
*
        call getxyprop(xx,nx,yy,ny,xrev,yrev,xwrap)
        call getlatlonwindow(lat1,lat2,lon1,lon2,xx,nx,xwrap,1,yy
     +        ,ny,1,x1,x2,y1,y2,lwrite)
        call biggerwindow(xx,nx,yy,ny,xwrap,lon1,lon2,lat1,lat2,
     +       x1,x2,y1,y2)
        call enscutoutwindow(x1,x2,y1,y2,xx,nx,xwrap,xrev,1,yy,ny
     +       ,1,wx,wy,field,nx,ny,nens1,nens2,nperyear,firstyr
     +       ,lastyr,max(firstyr,yr1),min(lastyr,yr2),lwrite)
*       my routines may have shifted the longitude axis - shift back
        if ( abs(xx(1)-lon1).gt.abs(xx(1)-lon1-360) ) then
            if ( abs(xx(nx)-lon2).gt.abs(xx(nx)-lon1-360) ) then
                xx = xx - 360
            endif
        endif
*
***        print '(a,f7.3,a,f7.3,a,f8.3,a,f8.3)','# lat = ',yy(1),' to '
***  +       ,yy(ny),', lon = ',xx(1),' to ',xx(nx)
        nn = 0
        if ( .not.dump ) then
            if ( m1.eq.m2 .and. nperyear.le.12 ) then
                allocate(newfcst(nx,ny,1,yr1:yr2,nens1:nens2))
            else
                allocate(newfcst(nx,ny,nperyear,yr1:yr2,nens1:nens2))
            endif
        endif
        yr1s = yr1
        yr2s = yr2
        yrstart = yr2
        yrstop  = yr1
        do jy=1,ny
            call keepalive(jy,ny)
            do jx=1,nx
*
*               create timeseries from fields
*
                fcst(1:nperyear,firstyr:max(yr1-1,firstyr-1),0:nens2) =
     +               3e33
                fcst(1:nperyear,min(yr2+1,lastyr+1):lastyr,0:nens2) =
     +               3e33
                n = 0
                do yr=yr1,yr2
                    do mo=1,nperyear
                        do iens=nens1,nens2
                            if ( yr.lt.firstyr .or. yr.gt.lastyr )
     +                           then
                                fcst(mo,yr,iens) = 3e33
                            else
                                fcst(mo,yr,iens) =
     +                               field(jx,jy,mo,yr,iens)
                            endif
                        enddo
                    enddo       ! mo
                enddo           ! yr
*
*               anomalies
*
                if ( anom ) then
                    if ( lwrite ) print '(a)','# Taking anomalies '
                    do iens=nens1,nens2
                        call anomal(fcst(1,firstyr,iens),nperyear
     +                       ,nperyear,firstyr,lastyr,yr1,yr2)
                    enddo
                endif
*
*               sum
*
                if ( lsum.gt.1 ) then
                    if ( lwrite .and.jx.eq.1 .and. jy.eq.1 ) 
     +                   print '(a,i3)','# Summing series ',lsum
                    do iens=nens1,nens2
                        call sumit(fcst(1,yr1,iens),nperyear,nperyear
     +                       ,yr1,yr2,lsum2,'v')
                    enddo
                endif
*
*               detrending
*
                if ( ldetrend ) then
                    do iens=nens1,nens2
                        if ( lwrite ) print *,'Detrending ens ',iens
                        call detrend(fcst(1,firstyr,iens),nperyear,
     +                       nperyear,firstyr,lastyr,yr1,yr2,m1,m2,lsel)
                    enddo
                endif
*
*               copy ensemble members so that there is the same
*               number of valid ones at every time step
*
                ndup = 0
                if ( nens2.gt.nens1 .and. lmakeensfull ) then
                    call makeensfull(ndup,nperyear,fcst,nperyear,firstyr
     +                   ,lastyr,nens1,nens2,validens,lwrite)
                endif
*               
*               copy to newfcst
*
                do yr=yr1,yr2
                    do month=m1,m2
                        call getj1j2(j1,j2,month,nperyear,.false.)
                        do j=j1,j2
                            if ( m1.eq.m2 .and. nperyear.le.12 ) then
                                mo = 1
                            else
                                mo = j
                            endif
                            do iens=nens1,nens2
                                newfcst(jx,jy,mo,yr,iens) = 
     +                               fcst(j,yr,iens)
                                if ( fcst(j,yr,iens).lt.1e33 ) then
                                    yrstart = min(yrstart,yr)
                                    yrstop  = max(yrstop,yr)
                                endif
                            enddo ! j
                        enddo   ! iens
                    enddo       ! month
                enddo           ! yr
            enddo               ! jx
        enddo                   ! jy
*
*       write the field to netcdf
*
        if ( lwrite ) print *,'yr1,yr2 = ',yr1,yr2
        nt = (yr2-yr1+1)
        if ( m1.ne.m2 .or. nperyear.gt.12 ) nt = nt*nperyear
        allocate(itimeaxis(nt))
        if ( nz.gt.1 ) then
            write(0,*) 'extractfield: error: cannot handle '//
     +           'nz>1 yet: ',nz
            call abort
        else
            nz = 1
            ivars(1,1) = 0    
        endif
        if ( m1.ne.m2 .or.nperyear.gt.12 ) then
            mobegin = 1
        else
            mobegin = m1
            nperyear = 1
        endif
*
*       forecasts
*
        call writenc(plotfile,ncid,ntvarid,itimeaxis,nt,nx,xx,ny
     +           ,yy,nz,zz,nt,nperyear,yr1,mobegin,3e33,title,nvars
     +       ,vars,ivars,lvars,units,nens1,nens2)
        do iens=nens1,nens2
            do yr=yr1,yr2
                if ( m1.eq.m2 .and. nperyear.le.12 ) then
                    call writencslice(ncid,ntvarid,itimeaxis,nt
     +                   ,ivars,newfcst(1,1,1,yr,iens),nx,ny,nz,nx
     +                   ,ny,nz,yr-yr1+1,iens-nens1+1)
                else
                    do mo=1,nperyear
                        call writencslice(ncid,ntvarid,itimeaxis,nt
     +                       ,ivars,newfcst(1,1,mo,yr,iens),nx,ny
     +                       ,nz,nx,ny,nz,mo+nperyear*(yr-yr1),iens
     +                       -nens1+1)
                    enddo
                endif
            enddo
        enddo
        i = nf_close(ncid)
*
*       and export yrstart, yrstop with interval with valid data
*
        call savestartstop(yrstart,yrstop)
        end
