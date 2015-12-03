        program get_index
!
!       get an index from a gridded (GrADS) data file
!       14-jan-2000: fixed bug that I cannot cut out a region outside 
!                    0..360 longitude
!       7-jul-2000: rewritten to enable getting the data from the .ctl file,
!                   .nc file not yet implemented
!       13-jul-2000 better treatment of missing points: compute
!                   anomalies first.
!       8-sep-2000 bring in synch with the other programs, conserve RAM
!       17-jun-005 add L/S mask option
!       1-nov-2005 adapted for daily data
!       2-feb-2010 added field output
!       22-jun-2012 added mask option
!
!       This file is part of the KNMI Climate Explorer and may be modified
!       freely as long as the below copyright notoice is kept intact.
!       (c) 1999-2012 KNMI
!
        implicit none
#include "params.h"
#include "netcdf.inc"
        integer nyrmax
        parameter (nyrmax=19)
        integer nx,ny,nz,nt,nperyear,firstyr,firstmo,lastyr,nvars,
     +        ivars(2,nvmax),endian,status,ncid,jvars(6,nvmax)
        integer,allocatable :: nn(:,:,:)
        real xx(nxmax),yy(nymax),zz(nzmax),undef
        real,allocatable :: field(:,:,:,:),mean(:,:,:)
        logical lwrite
        character file*255,datfile*255,title*1023,vars(nvmax)*40
     +        ,lvars(nvmax)*100,units(nvmax)*60
        integer iargc
*
        lwrite = .false.
        call getarg(iargc(),file)
        if ( file.eq.'debug' .or. file.eq.'lwrite' ) lwrite = .true.
        if ( iargc().lt.3 ) then
            print *,'usage: get_index file.[nc|ctl] '
     +            //'lon1 lon2 lat1 lat2|file listfile|mask file '
     +            //'[minfac r] [dipole tt|bb|ll|rr|tl|tr|bl|br] '
     +            //'[noisemodel 1|2] [lsmask file all|land|sea] '
     +            //'[interp|nearest] [nomissing] '
     +            //'[gridpoints fieldname] [outfield outfile.nc]'
     +           //' [debug]'
            stop
        endif
        call killfile(file,datfile,title,1)
*
        call getarg(1,file)
        if ( lwrite ) print *,'get_index: nf_opening file '
     +        ,trim(file)
        status = nf_open(file,nf_nowrite,ncid)
        if ( status.ne.nf_noerr ) then
            call parsectl(file,datfile,nxmax,nx,xx,nymax,ny,yy,nzmax,nz
     +            ,zz,nt,nperyear,firstyr,firstmo,undef,endian,title,1
     +            ,nvars,vars,ivars,lvars,units)
            ncid = -1
        else
            call parsenc(file,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax
     +            ,nz,zz,nt,nperyear,firstyr,firstmo,undef,title,1,nvars
     +            ,vars,jvars,lvars,units)
            datfile = file
        endif
*       range of years
        if ( nperyear.eq.366 ) then
            lastyr = firstyr + int((firstmo+nt-2)/365.24)
            if ( lwrite ) print *,'lastyr: ',lastyr,' = ', firstyr,
     +           ' + (',firstmo,'+',nt,'-2)/',365.24
        else
            lastyr = firstyr + (firstmo+nt-2)/nperyear
            if ( lwrite ) print *,'lastyr: ',lastyr,' = ', firstyr,
     +           ' + (',firstmo,'+',nt,'-2)/',nperyear
        end if
*
*       allocate arrays
*
        allocate(nn(nx,ny,nperyear))
        allocate(mean(nx,ny,nperyear))
        if ( lwrite ) print *,'allocating field(',nx,ny,nperyear,
     +       firstyr,':',lastyr,'), size ',4*nx*ny*nperyear*
     +       (lastyr-firstyr+1)
        allocate(field(nx,ny,nperyear,firstyr:lastyr),stat=status)
        if ( status.ne.0 ) then
            write(0,*) 'could not allocate field'
            call abort
        end if
*
*       now superfluous...
*
        call gindx(file,datfile,ncid,field,mean,nn,undef,endian
     +       ,nx,xx,ny,yy,nz,zz,nt,nperyear,firstyr,firstmo,lastyr,jvars
     +       ,title,vars,lvars,units,lwrite)
        end

        subroutine gindx(file,datfile,ncid,field,mean,nn,undef
     +       ,endian,nx,xx,ny,yy,nz,zz,nt,nperyear,firstyr,firstmo
     +       ,lastyr,jvars,title,vars,lvars,units,lwrite)
*
*       use field more efficiently
*
        use lsdata
        implicit none
#include "recfac.h"
#include "params.h"
#include "netcdf.inc"
        integer mypermax
        parameter(mypermax=4*366)
        real absent
        parameter (absent=3e33)
*
        integer ncid,endian,nx,ny,nz,nt,nperyear,firstyr,firstmo,lastyr,
     +        jvars(6,nvmax),nn(nx,ny,nperyear)
        real field(nx,ny,nperyear,firstyr:lastyr),undef,xx(nx),yy(ny)
     +       ,zz(nz),mean(nx,ny,nperyear)
        character file*(*),datfile*(*),title*(*),vars(1)*(*),
     +       lvars(1)*100,units(1)*(*)
        logical lwrite
*
        integer i,j,k,k1,n,x1,x2,y1,y2,year,month,iskip,sgn,ii,
     +       npoints,ipoints,in,mopts,out,noisemodel,status,nvars,
     +       ivars(2,1),lsncid,ndpm,yr1,yr2,nyr,ntvarid,nxf,nyf,
     +       nzf,it,mo,yr,irec
        integer ncidmask,nxmask,nymask,nzmask,ntmask,fyr,fmo,
     +       jvarsmask(6,1)
        integer,allocatable :: itimeaxis(:)
        integer xlist(nxmax*nymax),ylist(nxmax*nymax)
        real minfac,ave(mypermax),avemean(mypermax),lon1,lon2,lat1,lat2,
     +       wx(nxmax),wy(nymax),f,w,sw,sw0,pi,lon1c,lon2c,lat1c,
     +       lat2c,dx,dy,dl,xxls(nxmax),yyls(nymax),zzls(1),meanish
     +       ,offset,slope,undefmask
        real xxmask(nxmax),yymask(nymax),zzmask(nzmax)
        real,allocatable :: mask(:,:)
        logical xrev,yrev,xwrap,interp,missing,lstandardunits
     +       ,gridpoints,anom,outfield,lmask
        character string*255,dipole*2,shortfilename*255,datadir*255
     +       ,letter*1,lsmaskfile*255,lsmasktype*4,lstitle*255,
     +       fieldname*100,outfile*255
        character varsmask*40,lvarsmask(1)*80,titlemask*1000,unitsmask*40
        character lsvars(1)*10,lslvars(1)*40,lsunits(1)*20,newunits*20
        integer iargc
        integer,external :: normx,get_endian,leap
*
*       init
*
        pi = 4*atan(1d0)
        minfac = 40
        interp = .false.
        missing = .true.
        anom = .false.
        lstandardunits = .false.
        gridpoints = .false.
        outfield = .false.
        datadir = './data'
        lmask = .false.
        allocate(mask(nx,ny))
        mask = 1 ! default
*
*       the area to be cut out
*
        call getarg(2,string)
        if ( string(1:4).eq.'file' ) then
            npoints = 1000000
            call getarg(3,string)
            call rsunit(in)
            open(in,file=string,status='old')
            read(in,'(a)') datadir
            mopts = 4
        else if ( string(1:4).eq.'mask' ) then
            npoints = 1 ! for gthe time being only one mask can be requested at a time
            lmask = .true.
            call getarg(3,string)
            write(0,'(2a)') 'cutting out region defined by mask ',
     +           trim(string)
            write(*,'(2a)') '# cutting out region defined by mask ',
     +           trim(string)
            ncidmask = 0
            call parsenc(string,ncidmask,nxmax,nxmask,xxmask,
     +           nymax,nymask,yymask,nzmax,nzmask,zzmask,ntmask,n,
     +           fyr,fmo,undefmask,titlemask,1,nvars,varsmask,jvarsmask,
     +           lvarsmask,unitsmask)
            if ( nzmask.gt.1 ) then
                write(0,*) 'get_index: error: can only handle X-Y ',
     +               'masks, not with ',nzmask,' levels'
                call abort
            end if
            call checkgridequal(nx,ny,xx,yy,nxmask,nymask,xxmask,yymask)
            call readonencfield(ncidmask,jvarsmask,mask,nxmask,nymask,
     +           lwrite)
            call checklsmask(mask,nx,ny,lwrite)
            mopts = 4
        else
            npoints = 1
            read(string,*,err=901) lon1
            call getarg(3,string)
            read(string,*,err=902) lon2
            call getarg(4,string)
            read(string,*,err=903) lat1
            call getarg(5,string)
            read(string,*,err=904) lat2
            if ( lat2.lt.lat1 ) then
                w = lat2
                lat2 = lat1
                lat1 = w
            endif
            mopts = 6
        endif
*
*       rest of options
*
        iskip = 0
        dipole = 'no'
        lsmasktype = 'all '
        do i=mopts,iargc()
            if ( iskip.gt.0 ) then
                iskip = iskip - 1
            else
                call getarg(i,string)
                if ( string(1:6).eq.'minfac' ) then
                    call getarg(i+1,string)
                    read(string,*) minfac
                    iskip = 1
                    if ( missing ) then
                        print '(a,f6.2)',
     +                       '# using minimal fraction of valid points'
     +                       ,minfac
                    endif
                endif
                if ( string(1:6).eq.'interp' ) then
                    interp = .TRUE.
                endif
                if ( string(1:7).eq.'nearest' ) then
                    interp = .FALSE.
                endif
                if ( string(1:4).eq.'anom' ) then
                    anom = .TRUE.
                    write(*,'(a)') '# anomalies'
                endif
                if ( string(1:6).eq.'dipole' ) then
                    call getarg(i+1,dipole)
                    iskip = 1
                    if ( dipole.ne.'no' ) then
                        print '(2a)','# using dipole ',dipole
                    endif
                endif
                if ( string(1:5).eq.'noise' ) then
                    call getarg(i+1,string)
                    iskip = 1
                    read(string,*) noisemodel
                    if (noisemodel.ne.1 ) then
                        write(0,*) 'noisemodel ',noisemodel
     +                        ,' not yet ready'
                        call abort
                    endif
                endif
                if ( string(1:6).eq.'lsmask' ) then
                    iskip = 2
                    j = i ! getlsmask overwrites its first argument :-(
                    call getlsmask(j,lsmasktype,nxmax,xxls,nymax
     +                   ,yyls,lwrite)
                    if ( lsmasktype.ne.'all' ) then
                        call checkgridequal(nx,ny,xx,yy,nxls,nyls,
     +                       xxls,yyls)
                    end if
                endif
                if ( string(1:6).eq.'lwrite' .or. string(1:5).eq.'debug'
     +                ) then
                    lwrite = .TRUE.
                    if ( lwrite ) print *
     +                    ,'Debugging information requested'
                endif
                if ( string(1:5).eq.'nomis' ) then
                    missing = .false.
                    if ( lwrite ) print *,'No missing data'
                endif
                if ( string(1:4).eq.'grid' ) then
                    npoints = 100000
                    gridpoints = .true.
                    iskip = 1
                    call getarg(i+1,fieldname)
                    if ( fieldname.eq.' ' ) then
                        write(0,*)
     +                       'get_index: error: expecting fieldname'
                        write(*,*)
     +                       'get_index: error: expecting fieldname'
                        call abort
                    endif
!                   if the field is home-constructed, it may have slashes in it
                    fieldname=fieldname(1+index(fieldname,'/',.true.):)
                    if ( lwrite ) print *,'make grid points of '
     +                   ,trim(fieldname),' in region'
*                   no possibility to skip points yet
                endif
                if ( string(1:4).eq.'outf' ) then
                    outfield = .true.
                    iskip = 1
                    call getarg(i+1,outfile)
                    if ( outfile.eq.' ' ) then
                        write(0,*)
     +                       'get_index: error: expecting outfile'
                        write(*,*)
     +                       'get_index: error: expecting outfile'
                        call abort
                    end if
                end if
                if ( string(1:13).eq.'standardunits' ) then
                    lstandardunits = .true.
                    if ( lwrite ) print *
     +                    ,'Converting to standard units'
                endif
            endif
        enddo
        if ( interp ) then
            letter = 'i'
        else
            letter = 'n'
        endif
*
*       read file
*
        call keepalive1('Reading file',0,5)
        do j=len(file),1,-1
            if ( file(j:j).eq.'/' ) goto 101
        enddo
  101   continue
        shortfilename = file(j+1:index(file,' ')-1)
        if ( ncid.eq.-1 ) then
            call readdatfile(datfile,field,nx,ny,nx,ny,nperyear,firstyr
     +            ,lastyr,firstyr,firstmo,nt,undef,endian,lwrite,firstyr
     +            ,lastyr,1,1)
        else
            call readncfile(ncid,field,nx,ny,nx,ny,nperyear,firstyr
     +            ,lastyr,firstyr,firstmo,nt,undef,lwrite,firstyr,lastyr
     +            ,jvars)
        endif
        call keepalive1('Reading file',1,5)
*
*       manage units - I convert them at the end to save time
*
        if ( lstandardunits ) then
            call estimatemean(field,nx,ny,1,nperyear,firstyr,lastyr,
     +           nx,ny,1,nperyear,firstyr,lastyr,meanish,lwrite)
            call makestandardunits(meanish,nperyear,vars(1),units(1)
     +           ,newunits,offset,slope,ndpm,lwrite)
        else
            newunits = units(1)
        endif
        call keepalive1('Converted units',2,5)
        do yr=firstyr,lastyr
            do mo=1,nperyear
                do j=1,ny
                    do i=1,nx
                        if ( field(i,j,mo,yr).lt.-999 .and.
     +                       field(i,j,mo,yr).gt.-1000 ) then
                            write(0,*) 'get_index: suspicious field(',
     +                           i,j,mo,yr,') = ', field(i,j,mo,yr)
                        end if
                    end do
                end do
            end do
        end do
        call getxyprop(xx,nx,yy,ny,xrev,yrev,xwrap)
!
!       field output
!
        if ( outfield ) then
            call getlatlonwindow(lat1,lat2,lon1,lon2,xx,nx,xwrap,1,yy,ny
     +           ,1,x1,x2,y1,y2,lwrite)
            if ( (y2-y1)*(x2-x1).gt.0.2*nx*ny ) then
                write(0,*) 'Region too large (',(y2-y1)*(x2-x1),
     +           	'&gt;0.2*',nx*ny,').  Please download '//
     +               	'the whole dataset.'
                call abort
            end if
            nxf = nx
            nyf = ny
            nzf = nz
            call enscutoutwindow(x1,x2,y1,y2,xx,nx,xwrap,xrev,1,yy
     +           ,ny,1,wx,wy,field,nx,ny,0,0,nperyear,firstyr,lastyr
     +           ,firstyr,lastyr,lwrite)
            if ( lstandardunits) then
                call makestandardfield(field,nxf,nyf,1,nperyear,firstyr
     +               ,lastyr,nx,ny,1,nperyear,firstyr,lastyr,vars(1)
     +               ,units(1),lwrite)
            end if
            if ( index(outfile,'.ctl').eq.0 ) then
                allocate(itimeaxis(nt))
                ivars(1,1) = 0
                if ( nperyear.eq.366 ) then
                    yr=firstyr
                    mo=firstmo
                    irec = 0
                    do it=1,nt
                        if ( nperyear.eq.366 .and. mo.eq.60 .and.
     +                       leap(yr).eq.1 ) then
                            mo = mo + 1
                            cycle
                        end if
                        irec = irec + 1
                        mo = mo + 1
                        if ( mo.gt.nperyear ) then
                            mo = mo - nperyear
                            yr = yr + 1
                        end if
                    end do
                else
                    irec = nt
                end if
                call writenc(outfile,ncid,ntvarid,itimeaxis,nt,nx,xx
     +               ,ny,yy,nz,zz,irec,nperyear,firstyr,firstmo,3e33
     +               ,title,1,vars,ivars,lvars,units,0,0)
                yr=firstyr
                mo=firstmo
                irec = 0
                do it=1,nt
                    if ( nperyear.eq.366 .and. mo.eq.60 .and.
     +                   leap(yr).eq.1 ) then
                        mo = mo + 1
                        cycle
                    end if
                    irec = irec + 1
                    call writencslice(ncid,ntvarid,itimeaxis,nt,ivars
     +                   ,field(1,1,mo,yr),nxf,nyf,nzf,nx,ny,nz,irec,1)
                    mo = mo + 1
                    if ( mo.gt.nperyear ) then
                        mo = mo - nperyear
                        yr = yr + 1
                    end if
                end do
                status = nf_close(ncid)
            else
                write(0,*) 'GrADS .ctl/.grd file output not yet ready'
                call abort
            end if
            goto 999
        end if
*
*       get mean of whole field if multiple points are requested
*
        if ( npoints.gt.1 .and. (missing.or.anom) ) then
            call getmean(mean,nn,nx,ny,nperyear,field,nx,ny
     +            ,nperyear,firstyr,lastyr,nx,ny,firstyr,firstmo,nt
     +            ,lwrite)
        endif
        call keepalive1('Computed climatology',3,5)
*
*       compute weights
*
        call getweights('x',xx,wx,nx,xwrap,lwrite)
        call getweights('y',yy,wy,ny,.FALSE.,lwrite)
*
*       manage collection of grid points
*
        if ( gridpoints ) then
            write(*,'(2a,4(f10.3,a))') trim(fieldname),
     +           ' grid points in ',
     +           lat1,'N:',lat2,'N,',lon1,'E:',lon2,'E'
            call getlonwindow(lon1,lon2,x1,x2,xx,nx,lon1c,lon2c,lwrite)
            if ( lon1c.gt.1e33 .or. lon2c.gt.1e33 ) then
                write(0,*)
     +               'get_index: something went wrong in getlonwindow'
     +               ,lon1c,lon2c
                call abort
            endif
            call getlatwindow(lat1,lat2,y1,y2,yy,ny,lat1c,lat2c,lwrite)
            if ( lat1c.gt.1e33 .or. lat2c.gt.1e33 ) then
                write(0,*)
     +               'get_index: something went wrong in getlatwindow'
     +               ,lat1c,lat2c
                call abort
            endif

            npoints = 0
            do j=y1,y2
                do i=x1,x2
                    ii = normx(i,nx)
                    if ( lsmasktype.eq.'land' .or. 
     +                   lsmasktype.eq.'notl' ) then
                        if ( abs(lsmask(ii,j)-1).gt.1e-4 .eqv.
     +                       lsmasktype.eq.'land' ) then
                            if ( lwrite ) print *,'not ',lsmasktype,
     +                           ' point ',ii,j,lsmask(ii,j)
                            cycle
                        endif
                    elseif ( lsmasktype.eq.'sea ' .or. 
     +                       lsmasktype.eq.'nots' ) then
                        if ( abs(lsmask(ii,j)-0).gt.1e-4 .eqv.
     +                       lsmasktype.eq.'sea' ) then
                            if ( lwrite ) print *,'not ',lsmasktype,
     +                           ' point ',ii,j,lsmask(ii,j)
                            cycle
                        endif
                    elseif ( lsmasktype.eq.'5lan' .or. 
     +                       lsmasktype.eq.'5sea' ) then
                        if ( lsmask(ii,j).lt.0.5 .eqv.
     +                       lsmasktype.eq.'5lan' ) then
                            if ( lwrite ) print *,'not ',lsmasktype,
     +                           ' point ',ii,j,lsmask(ii,j)
                            cycle
                        endif
                    endif
                    npoints = npoints + 1
                    xlist(npoints) = ii
                    ylist(npoints) = j
                enddo
            enddo
            print '(a,i8,a)','found ',npoints,' grid points'
            print '(a)','========='
        endif
*
*       loop over all points
*
        if ( npoints.eq.1 .and. .not. gridpoints) then
            out = 6
        else
            call rsunit(out)
        endif
        do ipoints=1,npoints
            yr1 = +10000
            yr2 = -10000
            nyr = 0
            if ( npoints.gt.1 .and. .not.gridpoints ) then
                fieldname = shortfilename
                read(in,*,end=999) lon1,lat1
                lon2 = lon1
                lat2 = lat1
            endif
            if ( gridpoints ) then
                lon1 = xx(xlist(ipoints))
                lon2 = lon1
                lat1 = yy(ylist(ipoints))
                lat2 = lat1
                write(*,'(a,f6.2,a,f7.2,a)') 'coordinates: ',lat1,'N, '
     +               ,lon1,'E'
!               stationlist expects the pattern ':(.*)N.*,(.*)E'
                write(string,'(a,f7.2,a,f6.2,2a)') 'grid point: _',lon1,
     +               '_',lat1,'_',letter
                do i=14,len_trim(string)
                    if ( string(i:i).eq.' ' ) string(i:i) = '0'
                enddo
                write(*,'(3a)') trim(string),' ',trim(fieldname)
!               this should be identical to the file name opened below
            endif
            if ( npoints.gt.1 .or. gridpoints ) then
                write(string,'(4a,f7.2,a,f6.2,3a)')
     +                trim(datadir),'/grid'
     +               ,trim(fieldname),'_',lon1,'_',lat1,'_',letter
     +               ,'.dat'
                do i=1,len_trim(string)
                    if ( string(i:i).eq.' ' ) string(i:i) = '0'
                enddo
                if ( lwrite ) print *,'opening ',trim(string)
                open(out,file=string,status='unknown',err=900)
            endif
*
            if ( title.eq.' ' ) title = shortfilename
            write(out,'(6a)') '# ',trim(vars(1)),' [',trim(newunits)
     +           ,'] from ',trim(title)
*
*           compute indices of region to be cut out
*
            if ( lwrite ) then
                write(0,'(a,i4,1000f7.1)') 'get_index: found X axis ',nx
     +                ,(xx(i),i=1,nx)
                write(0,'(a,i4,1000f7.1)') 'get_index: found Y axis ',ny
     +                ,(yy(i),i=1,ny)
            endif
            if ( gridpoints ) then
                x1 = xlist(ipoints)
                x2 = x1
                y1 = ylist(ipoints)
                y2 = y1
                write(out,1000) '# grid point lon,lat =',lon1,lat1
            else if ( lmask ) then
                ! save time by computing the bounding box of the mask
                call getfirstnonzero(mask,nx,ny,'y',1,y1)
                call getfirstnonzero(mask,nx,ny,'y',-1,y2)
                call getfirstnonzero(mask,nx,ny,'x',1,x1)
                call getfirstnonzero(mask,nx,ny,'x',-1,x2)
                if ( lwrite ) print *,'mask is only non-zero over x ',
     +               x1,x2,' and y ',y1,y2
            else
                call getlonwindow(lon1,lon2,x1,x2,xx,nx,lon1c,lon2c
     +               ,lwrite)
                if ( lon1c.gt.1e33 .or. lon2c.gt.1e33 ) goto 900
                call getlatwindow(lat1,lat2,y1,y2,yy,ny,lat1c,lat2c
     +               ,lwrite)
                if ( lat1c.gt.1e33 .or. lat2c.gt.1e33 ) goto 900
                if ( interp ) then
                    if ( lon1.ne.lon2 .or. lat1.ne.lat2 ) then
                        write(0,*)
     +                       'get_index: error: cannot interpolate '/
     +                       /'area yet'
                        call abort
                    else
*                       find other points around the requested point
                        dx = lon1 - xx(normx(x1,nx))
                        if ( lwrite ) write(0,*) 'x1,x2,dx = ',x1,x2,dx
                        if ( abs(dx).gt.180 ) dx = dx - 360*nint(dx/360)
                        if ( dx.gt.0 .neqv. xrev) then
                            x2 = x1 + 1
                        else
                            x1 = x2 - 1
                            dx = lon1-xx(normx(x1,nx))
                            if ( abs(dx).gt.180 ) dx = dx - 360*nint(dx
     +                           /360)
                        endif
                        if ( lwrite ) write(0,*) 'x1,x2,dx = ',x1,x2,dx
     +                       ,xx(x1),xx(x2)
                        dl = xx(normx(x2,nx)) - xx(normx(x1,nx))
                        if ( abs(dl).gt.180 ) dl = dl - 360*nint(dl/360)
                        if ( .not.xrev .and. (dx.lt.0 .or. dx.gt.dl) .or
     +                       .xrev .and. (dx.gt.0 .or. dx.lt.dl) ) then
                            write(0,*) 'get_index: error: dx,dl = ',dx
     +                           ,dl,xrev
                        endif
                        ii = normx(x1,nx)
                        wx(ii) = wx(ii)*(1 - dx/dl)
                        lon1c = xx(ii)
                        ii = normx(x2,nx)
                        wx(ii) = wx(ii)*dx/dl
                        lon2c = xx(ii)
                        if ( lwrite ) write(0,*) 'w1,w2 = ',
     +                       wx(normx(x1,nx)),normx(x1,nx),
     +                       wx(normx(x2,nx)),normx(x2,nx)
*                       
                        dy = lat1 - yy(y1)
                        if ( lwrite ) write(0,*) 'y1,y2,dy = ',y1,y2,dy
                        if ( dy.gt.0 .neqv. yrev ) then
                            y2 = y1 + 1
                        else
                            y1 = y2 - 1
                            dy = lat1 - yy(y1)
                        endif
                        if ( lwrite ) write(0,*) 'y1,y2,dy = ',y1,y2,dy
     +                       ,yy(y1),yy(y2)
                        dl = yy(y2) - yy(y1)
                        if ( .not.yrev .and.(dy.lt.0 .or. dy.gt.dl) .or.
     +                       yrev .and. (dy.gt.0 .or. dy.lt.dl) ) then
                            write(0,*) 'get_index: error: dy,dl = ',dy
     +                           ,dl,yrev
                        endif
                        if ( lwrite ) write(0,*) 'w1,w2 = ',1 - dy/dl,dy
     +                       /dl
                        wy(y1) = wy(y1)*(1 - dy/dl)
                        lat1c = yy(y1)
                        wy(y2) = wy(y2)*dy/dl
                        lat2c = yy(y2)
                    endif
                    if ( npoints.eq.1 .and. .not.gridpoints ) then
                        write(0,1000) 'interpolating points lon=',lon1c
     +                       ,lon2c,', lat=',lat1c,lat2c,'<br>'
                    endif
                    write(out,1000) '# interpolating points lon=',lon1c
     +                   ,lon2c,', lat=',lat1c,lat2c
                    if ( lwrite ) write(0,'(a,2i4,a,2i4)')
     +                   'This corresponds to grid points x=',x1,x2
     +                   ,',y=',y1,y2
                else
 1000               format(a,2f9.3,a,2f9.3,a)
                    if ( npoints.eq.1 ) then
                        write(0,1000) 'cutting out region lon=',lon1c
     +                       ,lon2c,', lat=',lat1c,lat2c,'<br>'
                    endif
                    write(out,1000) '# cutting out region lon=',lon1c
     +                   ,lon2c,', lat=',lat1c,lat2c
                    if ( lwrite ) write(0,'(a,2i4,a,2i4)')
     +                'This corrsponds to grid points x=',x1,x2,',y=',y1
     +                   ,y2
                endif
            endif
            call keepalive1('Computed weights',4,5)
*
*           get mean of just the requested area if only one point is
*           requested
*
            if ( npoints.eq.1 .and. (missing.or.anom) ) then
                call getwinmean(mean,nn,nx,ny,nperyear,field,nx,ny
     +                ,nperyear,firstyr,lastyr,nx,ny,firstyr,firstmo,nt
     +                ,x1,x2,y1,y2,lwrite)
            endif
            call keepalive1('Computed area mean',5,5)
*
*           cut out region
*
            if ( missing.or.anom ) then
                month = -nperyear ! the first round (month<1) compute avemean
                k1 = -nperyear+1
            else
                avemean = 0     ! for safety, it should not be used
                month = 0
                k1 = 1
            endif
            ave(1:nperyear) = 3e33
            year = firstyr
            do k=k1,nt+firstmo-1
                month = month + 1
                if ( month.gt.nperyear ) then
                    call outputave(out,year,ave,nperyear,lstandardunits,
     +                   offset,slope,ndpm,year,yr1,yr2,nyr)
                    month = month - nperyear
                    year = year + 1
                    do i=1,nperyear
                        ave(i) = 3e33
                    enddo
                endif
                if ( month.gt.0 ) then
                    ave(month) = 0
                else
                    avemean(month+nperyear) = 0
                endif
                sw = 0
                do j=y1,y2
                    do i=x1,x2
                        ii = normx(i,nx)
                        if ( month.gt.0 ) then
                            f = field(ii,j,month,year)
                        else
                            f = mean(ii,j,month+nperyear)
                        endif
                        if ( lsmasktype.eq.'land' .or.
     +                       lsmasktype.eq.'notl' ) then
                            if ( abs(lsmask(ii,j)-1).gt.1e-4 .eqv.
     +                           lsmasktype.eq.'land' ) then
                                f = 3e33
                                if ( lwrite ) print *,'not ',lsmasktype,
     +                               ' point ',ii,j,lsmask(ii,j)
                            endif
                        elseif ( lsmasktype.eq.'sea ' .or.
     +                           lsmasktype.eq.'nots' ) then
                            if ( abs(lsmask(ii,j)-0).gt.1e-4 .eqv.
     +                           lsmasktype.eq.'sea' ) then
                                f = 3e33
                                if ( lwrite ) print *,'not ',lsmasktype,
     +                               ' point ',ii,j,lsmask(ii,j)
                            endif
                        elseif ( lsmasktype.eq.'5sea' .or.
     +                           lsmasktype.eq.'5lan' ) then
                            if ( lsmask(ii,j).lt.0.5 .eqv.
     +                           lsmasktype.eq.'5lan' ) then
                                f = 3e33
                                if ( lwrite ) print *,'not ',lsmasktype,
     +                               ' point ',ii,j,lsmask(ii,j)
                            endif
                        endif
                        if ( f.lt.-999 .and. f.gt.-1000 ) then
                            write(0,*)'get_index: warning: suspicious ',
     +                           'point field(',ii,j,month,year,') = ',f
                        end if
                        if ( f.lt.1e33 ) then
                            if (  dipole.eq.'no' .or.
     +                        dipole.eq.'ll' .and. i.lt.(x1+x2)/2. .or.
     +                        dipole.eq.'rr' .and. i.gt.(x1+x2)/2. .or.
     +                        dipole.eq.'tt' .and. j.lt.(y1+y2)/2. .or.
     +                        dipole.eq.'bb' .and. j.gt.(y1+y2)/2. .or.
     +                        (dipole.eq.'bl' .or. dipole.eq.'lb').and.
     +                        (i-x1)*(y1-y2).lt.(x2-x1)*(j-y2) .or.
     +                        (dipole.eq.'tr' .or. dipole.eq.'rt').and.
     +                        (i-x1)*(y1-y2).gt.(x2-x1)*(j-y2) .or.
     +                        (dipole.eq.'br' .or. dipole.eq.'rb').and.
     +                        (i-x1)*(y2-y1).gt.(x2-x1)*(j-y1) .or.
     +                        (dipole.eq.'tl' .or. dipole.eq.'lt').and.
     +                        (i-x1)*(y2-y1).lt.(x2-x1)*(j-y1) ) then
                                sgn = +1
                            else
                                sgn = -1
                            endif
                            w = wx(ii)*wy(j)
                            if ( lmask ) w = w*mask(ii,j)
                            if ( lwrite .and. .false. ) then
                                if ( month.gt.0 ) then
                                    write(0,*) 'adding field(',ii,j
     +                                    ,month,year,') = ',field(ii,j
     +                                    ,month,year),mean(ii,j,month)
     +                                    ,sgn,w
                                else
                                    write(0,*) 'adding mean(',ii,j,
     +                                    month+nperyear,') = ',
     +                                    mean(ii,j,month+nperyear),sgn
     +                                    ,w
                                endif
                            endif
                            if ( month.gt.0 ) then
                                ave(month) = ave(month) + sgn*w*
     +                                field(ii,j,month,year) 
                                if ( missing.or.anom ) ave(month) =
     +                               ave(month) - sgn*w*mean(ii,j,month)
                            else
                                avemean(month+nperyear) = 
     +                               avemean(month+nperyear) + 
     +                               sgn*w*mean(ii,j,month+nperyear)
                            endif
                            sw = sw + w
                        else
                            if ( lwrite .and. .false. ) print *,
     +                           'get_index: invalid point'
                        endif
                    enddo
                enddo
                if ( month.gt.0 ) then
                    if ( sw.lt.minfac/100*sw0 .or. sw.eq.0 ) then
                        ave(month) = 3e33
                        if ( lwrite ) write(0,*) 'sw<minfac*area 1: ',sw
     +                        ,minfac/100,sw0,month
                    elseif ( missing .and. .not.anom .and.
     +                       avemean(month).gt.1e33 ) then
                        ave(month) = 3e33
                        if ( lwrite ) write(0,*) 'mean missing: ',
     +                       avemean(month)
                    else
                        ave(month) = ave(month)/sw
                        if ( missing .and. .not.anom ) ave(month) =
     +                       ave(month) + avemean(month)
                        if ( lwrite ) write(*,*) 'ave(',month,') = '
     +                        ,ave(month)
                    endif
                else
                    if ( sw.eq.0 ) then
                        avemean(month+nperyear) = 3e33
                        sw0 = 0
                        if ( lwrite ) write(*,*) 'sw<minfac*area 2: ',sw
     +                        ,minfac/100,month+nperyear
                    else
                        sw0 = sw
                        avemean(month+nperyear) =
     +                       avemean(month+nperyear)/sw
                        if ( lwrite ) write(*,*) 'avemean(',month
     +                        +nperyear,') = ',avemean(month+nperyear),
     +                       sw
                    endif
                endif
                call keepalive1('Time step',k-k1+1,nt+firstmo-k1)
            enddo               ! loop over all months in file
*           print last (possibly incomplete) record
            call outputave(out,year,ave,nperyear,lstandardunits,
     +           offset,slope,ndpm,year,yr1,yr2,nyr)
            if ( nyr.eq.0 ) then
                write(0,*) 'get_index: error: cannot find any data'
                call abort
            end if
            if ( gridpoints ) then
                if ( nyr.gt.0 ) then
                    print '(a,i4,a,i4,a,i4)','found ',nyr
     +                   ,' years with data in ',yr1,'-',yr2
                else
                    print '(a)','could not locate any data'
                endif
                print '(a)','========='
            endif
*           finito
            if ( npoints.gt.1 ) then
                close(out)
            endif
  900       continue
            call keepalive1('Writing grid points',ipoints,npoints)
        enddo
        goto 999
 901    write(0,*) 'get_index: error reading lon1 from ',trim(string)
        call abort
 902    write(0,*) 'get_index: error reading lon2 from ',trim(string)
        call abort
 903    write(0,*) 'get_index: error reading lat1 from ',trim(string)
        call abort
 904    write(0,*) 'get_index: error reading lat2 from ',trim(string)
        call abort
  999   continue
        if ( allocated(lsmask) ) deallocate(lsmask)
        end

        subroutine outputave(out,year,ave,nperyear,lstandardunits,offset
     +       ,slope,ndpm,yr,yr1,yr2,nyr)
        implicit none
        integer out,year,nperyear,ndpm,yr,yr1,yr2,nyr
        logical lstandardunits
        real ave(nperyear),offset,slope
        integer i,dy,mo,dpm365(12),dpm(12),init
        save init,dpm,dpm365
        integer leap
        data init /0/
        data dpm/ 31,29,31,30,31,30,31,31,30,31,30,31/
        data dpm365/ 31,28,31,30,31,30,31,31,30,31,30,31/
*
        if ( lstandardunits ) then
            do mo=1,nperyear
                if ( ave(mo).lt.1e33 ) then
                    if ( ndpm.ne.0 ) then
                        if ( nperyear.eq.366 .and. mo.eq.2 .and. 
     +                       leap(yr).eq.2 ) then
                            ave(mo) = ave(mo)*29.**ndpm
                        else
                            ave(mo) = ave(mo)*real(dpm(mo))
     +                           **ndpm
                        endif
                    endif
                    ave(mo) = ave(mo)*slope + offset
                endif
            enddo
        endif
*       
        do mo=1,nperyear
            if ( ave(mo).lt.1e33 ) goto 201
        enddo
        return
 201    continue
        yr1 = min(yr1,year)
        nyr = nyr + 1
        yr2 = max(yr2,year)
*
        call printdatfile(out,ave,nperyear,nperyear,yr,yr)
        end

        subroutine getfirstnonzero(mask,nx,ny,axis,sign,y1)
        implicit none
        integer nx,ny,sign,y1
        real mask(nx,ny)
        character axis
        integer i,j,i1,i2,di,j1,j2,dj
        logical foundit

        foundit = .false.
        if ( axis.eq.'y' ) then
            if ( sign.eq.1 ) then
                j1 = 1
                j2 = ny
                dj = +1
            else if ( sign.eq.-1 ) then
                j1 = ny
                j2 = 1
                dj = -1
            else
                write(0,*) 'getfirstnonzero: error: sign should be +/-1'
     +               ,', not ',sign
                write(*,*) 'getfirstnonzero: error: sign should be +/-1'
     +               ,', not ',sign
                call abort
            end if
            do j=j1,j2,dj
                do i=1,nx
                    if ( mask(i,j).gt.0 ) then
                        foundit = .true.
                        exit
                    end if
                end do
                if ( foundit ) exit
            end do
            y1 = j
        else if ( axis.eq.'x') then
           if ( sign.eq.1 ) then
                i1 = 1
                i2 = nx
                di = +1
            else if ( sign.eq.-1 ) then
                i1 = nx
                i2 = 1
                di = -1
            else
                write(0,*) 'getfirstnonzero: error: sign should be +/-1'
     +               ,', not ',sign,axis
                write(*,*) 'getfirstnonzero: error: sign should be +/-1'
     +               ,', not ',sign,axis
                call abort
            end if
            do i=i1,i2,di
                do j=1,ny
                    if ( mask(i,j).gt.0 ) then
                        foundit = .true.
                        exit
                    end if
                end do
                if ( foundit ) exit
            end do
            y1 = i
        else
            write(0,*) 'getfirstnonzero: error: axis should be x or y'
     +               ,', not ',axis
            write(*,*) 'getfirstnonzero: error: axis should x or y'
     +               ,', not ',axis
            call abort
        end if
        if ( .not.foundit ) then
            write(0,*) 'getfirstnonzero: error: mask is all zero ',axis
     +           ,sign
            write(*,*) 'getfirstnonzero: error: mask is all zero ',axis
     +           ,sign
            if ( axis.eq.'x' ) then
                print *,'i1,i2,di = ',i1,i2,di
            else
                print *,'j1,j2,dj = ',j1,j2,dj
            end if
            call abort
        end if
 
        end
        
