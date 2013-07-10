        program findmax
*
*       compute the maximum of a field along lattude or longitude
*       averaging over a band inthe other direction.  Output is a
*       time series of the position of the max or min.
*
        implicit none
#include "params.h"
#include "netcdf.inc"
#include "getopts.inc"
        integer nyrmax,nlevmax,mensmax
        parameter (nyrmax=121,nlevmax=1,mensmax=1)
        integer n,nx,ny,nz,nt,firstyr,lastyr,firstmo,nvars,
     +          ivars(2,nvmax),jvars(6,nvmax),ncid,endian,
     +          status,nperyear,mens1,mens
        logical lexist
        real xx(nxmax),yy(nymax),zz(nzmax),undef,
     +        field(nxmax,nymax,12,nyrmax,mensmax)
        character infile*255,datfile*255,outfile*255,line*255
     +        ,vars(nvmax)*10,lvars(nvmax)*40,title*255,
     +        units(nvmax)*10
        integer iargc,llen
*
*       for f2c
        common /ccorrelatefield/ field
*
*       process command line
*
        n = iargc()
        if ( n.lt.3 ) then
            write(0,*) 'usage: findmax infile.[ctl|nc] min|max|nul'//
     +           'lat|lon [lon1 x1] [lon2 x2] [lat1 y1] [lat2 y2]'
            stop
        endif
        call getarg(1,infile)
        if ( index(infile,'%').gt.0 .or. index(infile,'++').gt.0 ) then
            ensemble = .true.
            call filloutens(infile,0)
            write(0,*) 'findmax: error: ensembles not yet supported'
            call abort
        else
            mens1 = 0
            mens = 0
        endif
        if ( lwrite ) print *,'findmax: nf_opening file '
     +        ,infile(1:llen(infile))
        status = nf_open(infile,nf_nowrite,ncid)
        if ( status.ne.nf_noerr ) then
            call parsectl(infile,datfile,nxmax,nx,xx,nymax,ny,yy,nzmax
     +            ,nz,zz,nt,nperyear,firstyr,firstmo,undef,endian,title
     +            ,1,nvars,vars,ivars,lvars,units)
            if ( nz.ne.1 ) then
                write(0,*) 'findmax: error: cannot handle 3D fields yet'
     +               ,nz
                call abort
            endif
            ncid = -1
            if ( ensemble ) then
                do mens=1,nensmax
                    call getarg(1,line)
                    call filloutens(line,mens)
                    inquire(file=line,exist=lexist)
                    if ( .not.lexist ) goto 100
                enddo
  100           continue
                mens = mens - 1
                write(0,*) 'located ',mens+1,' ensemble members<br>'
            endif
        else
            call parsenc(infile,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax
     +            ,nz,zz,nt,nperyear,firstyr,firstmo,undef,title,1,nvars
     +            ,vars,jvars,lvars,units)
            if ( ensemble ) then
                do mens=1,nensmax
                    call getarg(1,line)
                    call filloutens(line,mens)
                    status = nf_open(line,nf_nowrite,ncid)
                    if ( status.ne.nf_noerr ) goto 200
                enddo
  200           continue
                mens = mens - 1
                write(0,*) 'located ',mens+1,' ensemble members<br>'
            endif
        endif
        lastyr = firstyr + (firstmo+nt-2)/nperyear
*       process arguments
        call getopts(4,n,nperyear,yrbeg,yrend,.true.,mens1,mens)
        if ( ensemble ) write(0,*) 'Using ensemble members ',nens1
     +        ,' to ',nens2,'<br>'
        call getarg(1,infile)
*
*       check dimensions
*
        if ( nx*ny*nz*nperyear*(lastyr-firstyr+1)*(nens2-nens1+1)
     +        .gt.nxmax*nymax*12*nyrmax*mensmax ) then
            write(0,*) 'correlatefield: error: field too large '
            write(0,*) '  nx       = ',nx
            write(0,*) '  ny       = ',ny
            write(0,*) '  nz       = ',nz
            write(0,*) '  nperyear = ',nperyear
            write(0,*) '  years    = ',firstyr,lastyr
            write(0,*) '  ensemble = ',nens1,nens2
            write(0,*) 'total request',
     +            nx*ny*nz*nperyear*(lastyr-firstyr+1)*(nens2-nens1+1)
            write(0,*) 'available    ',
     +            nxmax*nymax*12*nyrmax*mensmax
            write(*,*) 'correlatefield: error: field too large '
            write(*,*) '  nx       = ',nx
            write(*,*) '  ny       = ',ny
            write(*,*) '  nz       = ',nz
            write(*,*) '  nperyear = ',nperyear
            write(*,*) '  years    = ',firstyr,lastyr
            write(*,*) '  ensemble = ',nens1,nens2
            write(*,*) 'total request',
     +            nx*ny*nz*nperyear*(lastyr-firstyr+1)*(nens2-nens1+1)
            write(*,*) 'available    ',
     +            nxmax*nymax*12*nyrmax*mensmax
            call abort
        endif
        if ( nx*ny*nz .gt. nxmax*nymax*nlevmax ) then
            write(0,*) 'correlatefield: error: fields too large: '
     +            ,nx,ny,nz,nxmax,nymax,nlevmax
            write(*,*) 'correlatefield: error: fields too large: '
     +            ,nx,ny,nz,nxmax,nymax,nlevmax
            call abort
        endif
*
*       to save on RAM usage
*
        call gfield(datfile,ncid,field,nx,xx,ny,yy,nz,zz,nt
     +        ,nperyear,firstyr,lastyr,firstmo,undef,endian,jvars)
        end

        subroutine gfield(datfile,ncid,field,nx,xx,ny,yy,nz,zz,nt
     +        ,nperyear,firstyr,lastyr,firstmo,undef,endian,jvars)
*
*       break to use field() compactly
*
        implicit none
#include "recfac.h"
#include "params.h"
#include "getopts.inc"
        real absent
        parameter (absent=3e33)
*
        integer ncid,endian,nx,ny,nz,nt,nperyear,firstyr,lastyr,
     +          firstmo,jvars(6,nvmax),x1,x2,y1,y2
        real field(nx,ny,nz,nperyear,firstyr:lastyr,nens1:nens2),
     +          undef,xx(nx),yy(ny),zz(nz),wx(nx),wy(ny),wz(nz)
        character datfile*(*)
*
        integer i,j,k,mo,yr,minmax,latlon,negpos,iens
        real fxy(npermax,yrbeg:yrend),vals(npermax),ff(nxmax),sum,wgt
     +       ,xmax,ymax
        logical xrev,yrev,xwrap,valid,lexist
        character line*256,dir*256,string*7
        integer iargc,llen,rindex
        external iargc,getarg,llen,rindex
*
        call getarg(2,line)
        if ( line(1:3).eq.'min' ) then
            minmax = -1
        elseif ( line(1:3).eq.'max' ) then
            minmax = +1
        elseif ( line(1:3).eq.'nul' ) then
            minmax = 0
            if ( line(4:6).eq.'neg' ) then
                negpos = -1
            elseif ( line(4:6).eq.'pos' ) then
                negpos = +1
            else
                negpos = 0
            endif
        else
            goto 901
        endif
        call getarg(3,line)
        if ( line(1:3).eq.'lat' ) then
            latlon = 1
        elseif ( line(1:3).eq.'lon' ) then
            latlon = 2
        else
            goto 902
        endif
*       range of years
        yr1 = max(yr1,firstyr)
        yr2 = min(yr2,firstyr + (firstmo+nt-2)/nperyear)
*
*       read field, change absent values to our convention
*
        if ( ensemble ) then
*           put the %% back in datfile...
            if ( nens2.lt.10 ) then
                i = 1
            elseif ( nens2.lt.100 ) then
                i = 2
            elseif ( nens2.lt.1000 ) then
                i = 3
            else
                write(0,*)
     +                'correlatefield: cannot handle ensembles up to '
     +                ,nens2,' yet'
                call abort
            endif
            string = '0000000'
            j = rindex(datfile,string(1:i))
            if ( j.eq.0 ) then
                write(0,*) 'correlatefield: error: cannot find '
     +                ,string(1:i),' in ',datfile(1:llen(datfile))
                call abort
            endif
            do k=j,j+i-1
                datfile(k:k) = '%'
            enddo
        endif
        do iens=nens1,nens2
            call keepalive(iens-nens1+1,nens2-nens1+1)
            if ( ncid.eq.-1 ) then
                dir=datfile
                if ( ensemble ) call filloutens(dir,iens)
                print *,'looking for '//dir(1:llen(dir))
                inquire(file=dir,exist=lexist)
                if ( .not.lexist ) then
                    print *,'looking for '//dir(1:llen(dir))//'.gz'
                    inquire(file=dir(1:llen(dir))//'.gz',exist=lexist)
                    if ( .not.lexist ) then
                        nens2 = iens-1
                        if ( nens2.ge.nens1 ) then
                            write(0,*) 'Found ensemble 0 to ',nens2
     +                            ,'<br>'
                            goto 5
                        else
                            write(0,*) 'Cannot locate file '
     +                            ,dir(1:llen(dir))
                            call abort
                        endif
                    endif
                endif
                if ( lwrite ) then
                    print *,'opening file ',dir(1:llen(dir))
                endif
                call zreaddatfile(dir,field(1,1,1,1,firstyr,iens),
     +                nx,ny,nz,nx,ny,nz,nperyear,firstyr,lastyr,
     +                firstyr,firstmo,nt,undef,endian,lwrite,yr1,yr2,1,1
     +                )
            else
                if ( nz.ne.1 ) then
                    write(0,*) 'cannot read 3D netCDF files yet'
                    call abort
                endif
                if ( ensemble ) then
                    write(0,*)
     +                    'cannot handle ensembles of netcdf files yet'
                    call abort
                endif
                call readncfile(ncid,field,nx,ny,nx,ny,nperyear,firstyr
     +                ,lastyr,firstyr,firstmo,nt,undef,lwrite,yr1,yr2
     +                ,jvars)
            endif
        enddo
    5   continue
*        
*       time series options
*
        if ( anom .or. lsum.gt.1 .or. ldetrend ) then
            do iens=nens1,nens2
                do k=1,nz
                    do j=1,ny
                        do i=1,nx
                            do yr=yr1,yr2
                                do mo=1,nperyear
                                    fxy(mo,yr) = field(i,j,k,mo,yr,iens)
                                enddo
                            enddo
*
*                           anomalies
*
                            if ( anom ) then
                                call anomal(fxy,npermax,nperyear,yrbeg
     +                               ,yrend,yr1,yr2)
                            endif
*       
*                           sum
*                           
                            if ( lsum.gt.1 ) then
                                call sumit(fxy,npermax,nperyear,yrbeg
     +                               ,yrend,lsum,oper)
                            endif
*
*                           detrend
*
                            if ( ldetrend ) then
                                if ( lwrite ) print *,'Detrending field'
                                call detrend(fxy,npermax,nperyear
     +                               ,yrbeg,yrend,yr1,yr2,1,12,1)
                            endif
                            do yr=yr1,yr2
                                do mo=1,nperyear
                                    field(i,j,k,mo,yr,iens) = fxy(mo,yr)
                                enddo
                            enddo
                        enddo   ! i
                    enddo       ! j
                enddo           ! k
            enddo               ! iens
        endif
*
*	get boundaries in grid points
*
        call getxyprop(xx,nx,yy,ny,xrev,yrev,xwrap)
        call getlatlonwindow(lat1,lat2,lon1,lon2,xx,nx,xwrap,avex,yy,ny
     +        ,avey,x1,x2,y1,y2,lwrite)
*
*       average, cut out window
*       
        call enscutoutwindow(x1,x2,y1,y2,xx,nx,xwrap,xrev,avex,yy,ny
     +       ,avey,wx,wy,field,nx,ny,nens1,nens2,nperyear,firstyr,lastyr
     +       ,yr1,yr2,lwrite)
        enddo
*
*       output header
*
        do iens=nens1,nens2
            if ( minmax.eq.-1 ) then
                write(0,'(a)')'# position of maximum along '
            elseif ( minmax.eq.+1 ) then
                write(0,'(a)')'# position of minimum along '
            elseif ( minmax.eq.0 ) then
                write(0,'(a)')'# position of zero along '
            else
                write(0,*) 'findmax: error: minmax = ',minmax
                call abort
            endif
            if ( latlon.eq.1 ) then
                write(0,'(a,f8.2,a,f8.2,a)')
     +               '# latitude (longitude averaged from ',lon1,' to '
     +               ,lon2,')'
            elseif ( latlon.eq.2 ) then
                write(0,'(a,f8.2,a,f8.2,a)')
     +               '# longitude (latitude averaged from ',lat1,' to '
     +               ,lat2,')'
            else
                write(0,*) 'findmax: error: latlon = ',latlon
                call abort
            endif
            write(0,'(2a)')'# of field in ',datfile(1:llen(datfile))
            write(0,'(a)')'#'
            write(0,'(a)')'#'
            do yr=yr1,yr2
                valid = .false.
                do mo=1,nperyear
*                   
*                   fill linear array
*                   
                    if ( latlon.eq.1 ) then
                        do j=1,ny
                            sum = 0
                            wgt = 0
                            do i=1,nx
                                if ( sum.lt.1e33 .and.
     +                               field(i,j,1,mo,yr,iens).lt.1e33 )
     +                               then
                                    sum = sum +wx(i)*field(i,j,1,mo,yr
     +                                   ,iens)
                                    wgt = wgt + wx(i)
                                else
                                    sum = 3e33
                                endif
                            enddo
                            if ( sum.lt.1e33 .and. wgt.ne.0 ) then
                                ff(j) = sign(1,minmax)*sum/wgt
                            else
                                ff(j) = 3e33
                            endif
                        enddo
                        if ( minmax.eq.0 ) then
                            call zerolin(xmax,negpos,xx,ff,nx)
                        else
                            call maxquad(xmax,ymax,yy,ff,ny)
                        endif
                    elseif ( latlon.eq.2 ) then
                        do i=1,nx
                            sum = 0
                            wgt = 0
                            do j=1,ny
                                if ( sum.lt.1e33 .and.
     +                               field(i,j,1,mo,yr,iens).lt.1e33 )
     +                               then
                                    sum = sum +wy(j)*field(i,j,1,mo,yr
     +                                   ,iens)
                                    wgt = wgt + wy(j)
                                else
                                    sum = 3e33
                                endif
                            enddo
                            if ( sum.lt.1e33 .and. wgt.ne.0 ) then
                                ff(i) = sign(1,minmax)*sum/wgt
                            else
                                ff(i) = 3e33
                            endif
                        enddo
                        if ( minmax.eq.0 ) then
                            call zerolin(xmax,negpos,xx,ff,nx)
                        else
                            call maxquad(xmax,ymax,xx,ff,nx)
                        endif
                        if ( lwrite .and. xmax.lt.1e33 ) then
                            print *,'findmax: found maximum at ',xmax
     +                           ,ymax,'in ',mo,yr
                            print *,'xx = '
                            print '(10f8.2)',(xx(i),i=1,nx)
                            print *,'ff = '
                            print '(10f8.2)',(ff(i),i=1,nx)
                        endif
                    else
                        write(0,*) 'findmax: error: latlon = ',latlon
                        call abort
                    endif
                    if ( xmax.lt.1e33 ) then
                        vals(mo) = xmax
                        valid = .true.
                    else
                        vals(mo) = -999.9
                    endif
                enddo
                if ( valid ) then
                    write(0,'(i4,366f8.2)')yr,(vals(mo),mo=1,nperyear)
                endif
            enddo
        enddo
*
*       error messages
*
        goto 999
  901   print *,'findmax: error reading min|max from '
     +        ,line(1:llen(line))
        call abort
  902   print *,'findmax: error reading lat|lon from '
     +        ,line(1:llen(line))
        call abort
  903   print *,'error reading date from file ',line(1:index(line,' ')-1
     +        ),' at record ',k
        call abort
  904   print *,'error cannot locate field file file ',line(1:index(line
     +        ,' ')-1)
        call abort
  920   print *,'error cannot open new correlations file '
     +        ,datfile(1:index(datfile,' ')-1)
        call abort
  999   continue
        end
