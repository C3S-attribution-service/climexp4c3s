program get_index

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

!       This file is part of the KNMI Climate Explorer and may be modified
!       freely as long as the below copyright notoice is kept intact.
!       (c) 1999-2012 KNMI

    use lsdata
    implicit none
    include 'params.h'
    include 'netcdf.inc'
    integer,parameter :: nyrmax=19
    integer,parameter :: mypermax=4*366,ntmax=5000000
    real,parameter :: absent=3e33
    integer :: nx,ny,nz,nt,nperyear,firstyr,firstmo,lastyr,nvars, &
        ivars(2,nvmax),endian,status,ncid,jvars(6,nvmax),nens1,nens2
    integer :: i,j,k,k1,n,x1,x2,y1,y2,year,month,iskip,sgn,ii, &
        npoints,ipoints,in,mopts,out,noisemodel, &
        lsncid,ndpm,yr1,yr2,nyr,ntvarid,nxf,nyf, &
        nzf,it,mo,yr,irec,ix,iy
    integer :: ncidmask,nxmask,nymask,nzmask,ntmask,fyr,fmo, &
        jvarsmask(6,1),iarray1(13),iarray2(13),iret1,iret2,iret
    integer :: xlist(nxmax*nymax),ylist(nxmax*nymax)
    integer,allocatable :: itimeaxis(:)
    integer,allocatable :: nn(:,:,:)
    real :: xx(nxmax),yy(nymax),zz(nzmax),undef
    real :: minfac,ave(mypermax),avemean(mypermax),lon1,lon2,lat1,lat2, &
        wx(nxmax),wy(nymax),f,w,sw,sw0,pi,lon1c,lon2c,lat1c, &
        lat2c,dx,dy,dl,xxls(nxmax),yyls(nymax),zzls(1),meanish &
        ,offset,slope,undefmask
    real :: xxmask(nxmax),yymask(nymax),zzmask(nzmax)
    real,allocatable :: mask(:,:)
    real,allocatable :: field(:,:,:,:),mean(:,:,:)
    logical :: lwrite
    logical :: tdefined(ntmax)
    logical :: xrev,yrev,xwrap,interp,missing,lstandardunits &
        ,gridpoints,anom,outfield,lmask,lexist
    character file*255,datfile*255,title*1023,vars(nvmax)*40 &
        ,lvars(nvmax)*100,units(nvmax)*60,FORM_field*100, &
        lz(3)*20,svars(100)*100,ltime*120,history*50000, &
        cell_methods(100)*100,longoper*100,metaoper*20,metadata(2,100)*2000
    character :: string*255,dipole*2,shortfilename*255,datadir*255 &
        ,letter*1,lsmaskfile*255,lsmasktype*4,lstitle*255, &
        fieldname*100,outfile*255
    character :: varsmask*40,lvarsmask(1)*80,titlemask*1000, &
        unitsmask*40
    character :: lsvars(1)*10,lslvars(1)*40,lsunits(1)*20,newunits*20
    character :: oper*1,masktype*4
    integer,external :: normx,get_endian,leap

    lwrite = .false.
    call get_command_argument(command_argument_count(),file)
    if ( file == 'debug' .or. file == 'lwrite' ) lwrite = .true.
    if ( command_argument_count() < 3 ) then
        print *,'usage: get_index file.[nc|ctl] ' &
        //'lon1 lon2 lat1 lat2|file listfile|mask file ' &
        //'[minfac r] [dipole tt|bb|ll|rr|tl|tr|bl|br] ' &
        //'[noisemodel 1|2] [lsmask file all|land|sea] ' &
        //'[interp|nearest] [nomissing] ' &
        //'[gridpoints fieldname] [outfield outfile.nc]' &
        //' [debug]'
        stop
    endif
    call killfile(file,datfile,title,1)

    call get_command_argument(1,file)
    if ( lwrite ) print *,'get_index: nf_opening file ',trim(file)
    status = nf_open(trim(file),nf_nowrite,ncid)
    if ( lwrite ) print *,'           returns ',status
    if ( status /= nf_noerr ) then
        history = ' '
        svars =  ' '
        cell_methods = ' '
        metadata = ' '
        call parsectl(file,datfile,nxmax,nx,xx,nymax,ny,yy,nzmax,nz &
            ,zz,nt,nperyear,firstyr,firstmo,undef,endian,title,1 &
            ,nvars,vars,ivars,lvars,units)
        ncid = -1
    else
        nens1 = 0
        nens2 = 0
        call ensparsenc(file,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax &
            ,nz,zz,lz,nt,nperyear,firstyr,firstmo,ltime,tdefined,ntmax &
            ,nens1,nens2,undef,title,history,1,nvars &
            ,vars,jvars,lvars,svars,units,cell_methods,metadata)
        datfile = file
    endif
    ! until I get a new server...
!   range of years
    if ( nperyear == 366 ) then
        lastyr = firstyr + int((firstmo+nt-2)/365.24)
        if ( lwrite ) print *,'lastyr: ',lastyr,' = ', firstyr, &
        ' + (',firstmo,'+',nt,'-2)/',365.24
    else
        lastyr = firstyr + (firstmo+nt-2)/nperyear
        if ( lwrite ) print *,'lastyr: ',lastyr,' = ', firstyr, &
        ' + (',firstmo,'+',nt,'-2)/',nperyear
    end if

!   allocate arrays

    allocate(nn(nx,ny,nperyear))
    allocate(mean(nx,ny,nperyear))
    if ( lwrite ) print *,'allocating field(',nx,ny,nperyear, &
        firstyr,':',lastyr,'), size ',4d0*nx*ny*nperyear*(lastyr-firstyr+1)
    allocate(field(nx,ny,nperyear,firstyr:lastyr),stat=status)
    if ( status /= 0 ) then
        write(0,*) 'get_index: error: could not allocate field, choose a lower-resolution version'
        call exit(-1)
    end if

!   init

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

!   the area to be cut out

    call get_command_argument(2,string)
    if ( string(1:4) == 'file' ) then
        npoints = 1000000
        call get_command_argument(3,string)
        call rsunit(in)
        open(in,file=string,status='old')
        read(in,'(a)') datadir
        mopts = 4
    else if ( string(1:4) == 'mask' ) then
        npoints = 1 ! for the time being only one mask can be requested at a time
        lmask = .true.
        call get_command_argument(3,string)
        write(0,'(3a)') 'cutting out region defined by mask ',trim(string),'<br>'
        write(*,'(2a)') '# cutting out region defined by mask ',trim(string)
        ncidmask = 0
        call parsenc(string,ncidmask,nxmax,nxmask,xxmask, &
            nymax,nymask,yymask,nzmax,nzmask,zzmask,ntmask,n, &
            fyr,fmo,undefmask,titlemask,1,nvars,varsmask,jvarsmask, &
            lvarsmask,unitsmask)
        if ( nzmask > 1 ) then
            write(0,*) 'get_index: error: can only handle X-Y masks, not with ',nzmask,' levels'
            call exit(-1)
        end if
        call checkgridequal(nx,ny,xx,yy,nxmask,nymask,xxmask,yymask)
        call readonencfield(ncidmask,jvarsmask,mask,nxmask,nymask,lwrite)
        call checklsmask(mask,nx,ny,lwrite)
        do i=1,100
            if ( metadata(1,i) == ' ' ) exit
        end do
        metadata(1,i) = 'mask'
        metadata(2,i) = string
        mopts = 4
    else
        npoints = 1
        read(string,*,err=901) lon1
        call get_command_argument(3,string)
        read(string,*,err=902) lon2
        call get_command_argument(4,string)
        read(string,*,err=903) lat1
        call get_command_argument(5,string)
        read(string,*,err=904) lat2
        if ( lat2 < lat1 ) then
            w = lat2
            lat2 = lat1
            lat1 = w
        endif
        mopts = 6
    endif

!   rest of options

    iskip = 0
    dipole = 'no'
    oper = 'v'
    lsmasktype = 'all '
    do i=mopts,command_argument_count()
        if ( iskip > 0 ) then
            iskip = iskip - 1
        else
            call get_command_argument(i,string)
            if ( string(1:6) == 'minfac' ) then
                call get_command_argument(i+1,string)
                read(string,*) minfac
                iskip = 1
                if ( missing ) then
                    print '(a,f6.2)','# minimal_valid_fraction :: ',minfac
                endif
            endif
            if ( string(1:6) == 'interp' ) then
                interp = .true.
            endif
            if ( string(1:7) == 'nearest' ) then
                interp = .false.
            endif
            if ( string(1:4) == 'anom' ) then
                anom = .true.
                write(*,'(a)') '# anomalies'
            endif
            if ( string(1:6) == 'dipole' ) then
                call get_command_argument(i+1,dipole)
                iskip = 1
                if ( dipole /= 'no' ) then
                    print '(2a)','# using dipole ',dipole
                endif
            endif
            if ( string(1:5) == 'noise' ) then
                call get_command_argument(i+1,string)
                iskip = 1
                read(string,*) noisemodel
                if (noisemodel /= 1 ) then
                    write(0,*) 'noisemodel ',noisemodel,' not yet ready'
                    call exit(-1)
                endif
            endif
            if ( string(1:6) == 'lsmask' ) then
                iskip = 2
                j = i ! getlsmask overwrites its first argument :-(
                call getlsmask(j,lsmasktype,nxmax,xxls,nymax,yyls,lwrite)
                if ( lsmasktype /= 'all' ) then
                    call checkgridequal(nx,ny,xx,yy,nxls,nyls,xxls,yyls)
                end if
            endif
            if ( string(1:6) == 'lwrite' .or. string(1:5) == 'debug' ) then
                lwrite = .true.
                if ( lwrite ) print *,'Debugging information requested'
            endif
            if ( string(1:5) == 'nomis' ) then
                missing = .false.
                if ( lwrite ) print *,'No missing data'
            endif
            if ( string(1:4) == 'grid' ) then
                npoints = 100000
                gridpoints = .true.
                iskip = 1
                call get_command_argument(i+1,fieldname)
                if ( fieldname == ' ' ) then
                    write(0,*) 'get_index: error: expecting fieldname'
                    write(*,*) 'get_index: error: expecting fieldname'
                    call exit(-1)
                endif
!               if the field is home-constructed, it may have slashes in it
                fieldname=fieldname(1+index(fieldname,'/', .true.):)
                if ( lwrite ) print *,'make grid points of ',trim(fieldname),' in region'
!               no possibility to skip points yet
            endif
            if ( string(1:4) == 'outf' ) then
                outfield = .true.
                iskip = 1
                call get_command_argument(i+1,outfile)
                if ( outfile == ' ' ) then
                    write(0,*) &
                    'get_index: error: expecting outfile'
                    write(*,*) &
                    'get_index: error: expecting outfile'
                    call exit(-1)
                end if
            end if
            if ( string(1:13) == 'standardunits' ) then
                lstandardunits = .true.
                if ( lwrite ) print * &
                ,'Converting to standard units'
            endif
            if ( string(1:3) == 'max' ) then
                oper = 'x'
                if ( lwrite ) print '(a)','# taking max over region'
            endif
            if ( string(1:3) == 'min' ) then
                oper = 'n'
                if ( lwrite ) print '(a)','# taking min over region'
            endif
        endif
    enddo
    if ( interp ) then
        letter = 'i'
    else
        letter = 'n'
    endif

!   read file

    call keepalive1('Reading file',0,5)
    do j=len(file),1,-1
        if ( file(j:j) == '/' ) goto 101
    enddo
101 continue
    shortfilename = file(j+1:index(file,' ')-1)
    if ( ncid == -1 ) then
        call readdatfile(datfile,field,nx,ny,nx,ny,nperyear,firstyr &
            ,lastyr,firstyr,firstmo,nt,undef,endian,lwrite,firstyr &
            ,lastyr,1,1)
    else
        call readncfile(ncid,field,nx,ny,nx,ny,nperyear,firstyr &
            ,lastyr,firstyr,firstmo,nt,undef,lwrite,firstyr,lastyr &
            ,jvars)
    endif
    call keepalive1('Reading file',1,5)

!   manage units - I convert them at the end to save time

    if ( lstandardunits ) then
        call estimatemean(field,nx,ny,1,nperyear,firstyr,lastyr &
            ,nx,ny,1,nperyear,firstyr,lastyr,meanish,lwrite)
        call makestandardunits(meanish,nperyear,vars(1),units(1) &
            ,newunits,offset,slope,ndpm,lwrite)
    else
        newunits = units(1)
    endif
    call keepalive1('Converted units',2,5)
    do yr=firstyr,lastyr
        do mo=1,nperyear
            do j=1,ny
                do i=1,nx
                    if ( field(i,j,mo,yr) < -999 .and. &
                         field(i,j,mo,yr) > -1000 ) then
                        write(0,*) 'get_index: suspicious field(', &
                            i,j,mo,yr,') = ', field(i,j,mo,yr)
                    end if
                end do
            end do
        end do
    end do
    call getxyprop(xx,nx,yy,ny,xrev,yrev,xwrap)

!   field output

    if ( outfield ) then
        if ( lmask ) then
            ! apply mask to field
            call getmaskbox(mask,nx,ny,xwrap,0.5,x1,x2,y1,y2,lwrite)
            masktype = '5lan'
            if ( x1 == 0 .and. x2 == 0 ) then
                ! only islands, try again
                call getmaskbox(mask,nx,ny,xwrap,0.1,x1,x2,y1,y2,lwrite)
                masktype = 'land'
            end if
            if ( x1 == 0 .and. x2 == 0 ) then
                ! only small islands, try again
                call getmaskbox(mask,nx,ny,xwrap,0.01,x1,x2,y1,y2,lwrite)
            end if
            if ( x1 == 0 .and. x2 == 0 ) then
                ! only tiny islands, give up
                write(0,*) 'get_index: error: cannot find points with weight > 0.01, giving up'
                call exit(-1)
            end if
            call applylsmask(field,mask,nx,ny,1,nperyear,firstyr,lastyr,0,0,masktype,lwrite)
            if ( lwrite ) then
                do iy=1,ny
                    do ix=1,nx
                        if ( field(ix,iy,1,firstyr) < 1e33 ) print *,'field(',ix,iy,'1,',firstyr,') = ',field(ix,iy,1,firstyr)
                    end do
                end do
            end if
            lon1 = xx(x1)
            if ( x2 > nx ) then
                lon2 = xx(x2-nx)
            else
                lon2 = xx(x2)
            end if
            lat1 = yy(y1)
            lat2 = yy(y2)
        end if
        call getlatlonwindow(lat1,lat2,lon1,lon2,xx,nx,xwrap,1,yy,ny,1,x1,x2,y1,y2,lwrite)
        if ( (y2-y1)*(x2-x1) > 0.2*nx*ny ) then
            write(0,*) 'Region too large (',(y2-y1)*(x2-x1),'&gt;0.2*',nx*ny, &
                ').  Please download the whole dataset.'
            call exit(-1)
        end if
        nxf = nx
        nyf = ny
        nzf = nz
        call keepalive1('Cutting out',0,1)
        call enscutoutwindow(x1,x2,y1,y2,xx,nx,xwrap,xrev,1,yy &
            ,ny,1,wx,wy,field,nx,ny,0,0,nperyear,firstyr,lastyr &
            ,firstyr,lastyr,lwrite)
        call keepalive1('Cutting out',1,1)
        if ( lstandardunits) then
            call makestandardfield(field,nxf,nyf,1,nperyear,firstyr &
                ,lastyr,nx,ny,1,nperyear,firstyr,lastyr,vars(1),units(1),lwrite)
        end if
        if ( index(outfile,'.ctl') == 0 ) then
            allocate(itimeaxis(nt))
            ivars(1,1) = 0
            call subtractleapyears(nt,firstyr,firstmo,nperyear,irec)
            title = 'subset of '//trim(title)
            call trim_geospatial_metadata(metadata,lon1,lon2,lat1,lat2)
            call enswritenc(outfile,ncid,ntvarid,itimeaxis,nt,nx,xx &
                ,ny,yy,nz,zz,lz,irec,nperyear,firstyr,firstmo,ltime,3e33 &
                ,title,history,1,vars,ivars,lvars,svars,units,cell_methods &
                ,metadata,0,0)
            yr=firstyr
            mo=firstmo
            irec = 0
            do it=1,nt
                if ( nperyear == 366 .and. mo == 60 .and. &
                leap(yr) == 1 ) then
                    mo = mo + 1
                    cycle
                end if
                irec = irec + 1
                call keepalive1('Writing day',it,nt)
                call writencslice(ncid,ntvarid,itimeaxis,nt,ivars &
                    ,field(1,1,mo,yr),nxf,nyf,nzf,nx,ny,nz,irec,1)
                mo = mo + 1
                if ( mo > nperyear ) then
                    mo = mo - nperyear
                    yr = yr + 1
                end if
            end do
            status = nf_close(ncid)
        else
            write(0,*) 'GrADS .ctl/.grd file output not yet ready'
            call exit(-1)
        end if
        goto 999
    end if

!   series metadata output

    print '(5a)','# operating on ',trim(title)
    call getenv('FORM_field',FORM_field)
    title = 'spatial statistic of '//title
    call delete_geospatial_metadata(metadata)
    call printmetadata(6,file,FORM_field,title,history,metadata)

!   get mean of whole field if multiple points are requested

    if ( npoints > 1 .and. (missing .or. anom) .and. oper == 'v' ) then
        call getmean(mean,nn,nx,ny,nperyear,field,nx,ny &
            ,nperyear,firstyr,lastyr,nx,ny,firstyr,firstmo,nt &
            ,lwrite)
    endif
    call keepalive1('Computed climatology',3,5)

!   compute weights

    call getweights('x',xx,wx,nx,xwrap,lwrite)
    call getweights('y',yy,wy,ny, .false.,lwrite)

!   manage collection of grid points

    if ( gridpoints ) then
        if ( lmask ) then
            !!!call get_command_argument(3,string)
            !!!write(*,'(3a)') trim(fieldname),' grid points within mask ',trim(string)
            ! compute the bounding box of the mask
            call getfirstnonzero(mask,nx,ny,'y',1,y1)
            call getfirstnonzero(mask,nx,ny,'y',-1,y2)
            call getfirstnonzero(mask,nx,ny,'x',1,x1)
            call getfirstnonzero(mask,nx,ny,'x',-1,x2)
            if ( lwrite ) print *,'mask is only non-zero over x ',x1,x2,' and y ',y1,y2
            lon1 = xx(x1)
            lon2 = xx(x2)
            lat1 = yy(y1)
            lat2 = yy(y2)
            write(*,'(2a,4(f10.3,a))') trim(fieldname),' grid points in ', &
                lat1,'N:',lat2,'N,',lon1,'E:',lon2,'E'
        else
            write(*,'(2a,4(f10.3,a))') trim(fieldname),' grid points in ', &
                lat1,'N:',lat2,'N,',lon1,'E:',lon2,'E'
            call getlonwindow(lon1,lon2,x1,x2,xx,nx,lon1c,lon2c,lwrite)
            if ( lon1c > 1e33 .or. lon2c > 1e33 ) then
                write(0,*) 'get_index: something went wrong in getlonwindow',lon1c,lon2c
                call exit(-1)
            endif
            call getlatwindow(lat1,lat2,y1,y2,yy,ny,lat1c,lat2c,lwrite)
            if ( lat1c > 1e33 .or. lat2c > 1e33 ) then
                write(0,*) 'get_index: something went wrong in getlatwindow',lat1c,lat2c
                call exit(-1)
            endif
        end if

        npoints = 0
        do j=y1,y2
            do i=x1,x2
                ii = normx(i,nx)
                if ( lsmasktype == 'land' .or. lsmasktype == 'notl' ) then
                    if ( abs(lsmask(ii,j)-1) > 1e-4 .eqv. lsmasktype == 'land' ) then
                        if ( lwrite ) print *,'not ',lsmasktype,' point ',ii,j,lsmask(ii,j)
                        cycle
                    endif
                elseif ( lsmasktype == 'sea ' .or. lsmasktype == 'nots' ) then
                    if ( abs(lsmask(ii,j)-0) > 1e-4 .eqv. lsmasktype == 'sea' ) then
                        if ( lwrite ) print *,'not ',lsmasktype,' point ',ii,j,lsmask(ii,j)
                        cycle
                    endif
                elseif ( lsmasktype == '5lan' .or. lsmasktype == '5sea' ) then
                    if ( lsmask(ii,j) < 0.5 .eqv. lsmasktype == '5lan' ) then
                        if ( lwrite ) print *,'not ',lsmasktype,' point ',ii,j,lsmask(ii,j)
                        cycle
                    endif
                endif
                if ( lmask ) then
                    if ( mask(ii,j) == 0 ) cycle
                end if
                npoints = npoints + 1
                xlist(npoints) = ii
                ylist(npoints) = j
            enddo
        enddo
        print '(a,i8,a)','found ',npoints,' grid points'
        print '(a)','========='
    endif

!   loop over all points

    if ( npoints == 1 .and. .NOT. gridpoints) then
        out = 6
    else
        call rsunit(out)
    endif
    do ipoints=1,npoints
        yr1 = +10000
        yr2 = -10000
        nyr = 0
        if ( npoints > 1 .and. .NOT. gridpoints ) then
            fieldname = shortfilename
            read(in,*,end=999) lon1,lat1
            lon2 = lon1
            lat2 = lat1
        endif
        if ( gridpoints ) then
            ! stationlist expects a country these days
            write(*,'(i8,a)') ipoints,' (grid)'
            lon1 = xx(xlist(ipoints))
            lon2 = lon1
            lat1 = yy(ylist(ipoints))
            lat2 = lat1
            write(*,'(a,f6.2,a,f7.2,a)') 'coordinates: ',lat1,'N, ',lon1,'E'
        !   stationlist expects the pattern ':(.*)N.*,(.*)E'
            write(string,'(a,f7.2,a,f6.2,2a)') 'grid point: _',lon1,'_',lat1,'_',letter
            do i=14,len_trim(string)
                if ( string(i:i) == ' ' ) string(i:i) = '0'
            enddo
            write(*,'(3a)') trim(string),' ',trim(fieldname)
        !   this should be identical to the file name opened below
        endif
        if ( npoints > 1 .or. gridpoints ) then
            write(string,'(4a,f7.2,a,f6.2,3a)') trim(datadir),'/grid' &
                ,trim(fieldname),'_',lon1,'_',lat1,'_',letter,'.dat'
            do i=1,len_trim(string)
                if ( string(i:i) == ' ' ) string(i:i) = '0'
            enddo
            inquire(file=trim(string),exist=lexist)
            if ( lexist ) then
            ! check whether it is up-to-date
                call mystat(trim(fieldname),iarray1,iret1)
                call mystat(trim(string),iarray2,iret2)
                if ( iret1 == 0 .and. iret2 == 0 ) then
                    if ( iarray2(10) >= iarray1(10) ) then
                        ! file already there and seems valid, next point
                        if ( lwrite ) print *,'file ',trim(string), &
                            'already there and seems valid', &
                            iarray2(10) >= iarray1(10)
                        cycle
                    end if
                end if
                if ( lwrite ) print *,'file ',trim(string), &
                    ' invalid or outdated ',iret1,iarray1(10), &
                    iret2,iarray2(10),', deleting'
                open(out,file=trim(string),status='old')
                close(out,status='delete')
            end if
            if ( lwrite ) print *,'opening ',trim(string)
            open(out,file=trim(string),status='unknown',err=900)
        endif
    
!       compute indices of region to be cut out
    
        if ( lwrite ) then
            write(0,'(a,i4,10000f7.1)') 'get_index: found X axis ',nx,(xx(i),i=1,nx)
            write(0,'(a,i4,10000f7.1)') 'get_index: found Y axis ',ny,(yy(i),i=1,ny)
        endif
        if ( gridpoints ) then
            x1 = xlist(ipoints)
            x2 = x1
            y1 = ylist(ipoints)
            y2 = y1
            write(out,'(a,2f9.3)') '# grid point lon,lat =',lon1,lat1
        else if ( lmask ) then
            ! save time by computing the bounding box of the mask
            call getfirstnonzero(mask,nx,ny,'y',1,y1)
            call getfirstnonzero(mask,nx,ny,'y',-1,y2)
            call getfirstnonzero(mask,nx,ny,'x',1,x1)
            call getfirstnonzero(mask,nx,ny,'x',-1,x2)
            if ( lwrite ) print *,'mask is only non-zero over x ', &
            x1,x2,' and y ',y1,y2
        else
            call getlonwindow(lon1,lon2,x1,x2,xx,nx,lon1c,lon2c,lwrite)
            if ( lon1c > 1e33 .or. lon2c > 1e33 ) goto 900
            call getlatwindow(lat1,lat2,y1,y2,yy,ny,lat1c,lat2c,lwrite)
            if ( lat1c > 1e33 .or. lat2c > 1e33 ) goto 900
            if ( interp ) then
                if ( lon1 /= lon2 .or. lat1 /= lat2 ) then
                    write(0,*) 'get_index: error: cannot interpolate area yet'
                    call exit(-1)
                else
!                   find other points around the requested point
                    dx = lon1 - xx(normx(x1,nx))
                    if ( lwrite ) write(0,*) 'x1,x2,dx = ',x1,x2,dx
                    if ( abs(dx) > 180 ) dx = dx - 360*nint(dx/360)
                    if ( dx > 0 .neqv. xrev) then
                        x2 = x1 + 1
                    else
                        x1 = x2 - 1
                        dx = lon1-xx(normx(x1,nx))
                        if ( abs(dx) > 180 ) dx = dx - 360*nint(dx/360)
                    endif
                    if ( lwrite ) write(0,*) 'x1,x2,dx = ',x1,x2,dx,xx(x1),xx(x2)
                    dl = xx(normx(x2,nx)) - xx(normx(x1,nx))
                    if ( abs(dl) > 180 ) dl = dl - 360*nint(dl/360)
                    if ( .NOT. xrev .and. (dx < 0 .or. dx > dl) .or. &
                               xrev .and. (dx > 0 .or. dx < dl) ) then
                        write(0,*) 'get_index: error: dx,dl = ',dx,dl,xrev
                    endif
                    ii = normx(x1,nx)
                    wx(ii) = wx(ii)*(1 - dx/dl)
                    lon1c = xx(ii)
                    ii = normx(x2,nx)
                    wx(ii) = wx(ii)*dx/dl
                    lon2c = xx(ii)
                    if ( lwrite ) write(0,*) 'w1,w2 = ', &
                        wx(normx(x1,nx)),normx(x1,nx), &
                        wx(normx(x2,nx)),normx(x2,nx)
                
                    dy = lat1 - yy(y1)
                    if ( lwrite ) write(0,*) 'y1,y2,dy = ',y1,y2,dy
                    if ( dy > 0 .neqv. yrev ) then
                        y2 = y1 + 1
                    else
                        y1 = y2 - 1
                        dy = lat1 - yy(y1)
                    endif
                    if ( lwrite ) write(0,*) 'y1,y2,dy = ',y1,y2,dy,yy(y1),yy(y2)
                    dl = yy(y2) - yy(y1)
                    if ( .NOT. yrev .and. (dy < 0 .or. dy > dl) .or. &
                               yrev .and. (dy > 0 .or. dy < dl) ) then
                        write(0,*) 'get_index: error: dy,dl = ',dy ,dl,yrev
                    endif
                    if ( lwrite ) write(0,*) 'w1,w2 = ',1 - dy/dl,dy/dl
                    wy(y1) = wy(y1)*(1 - dy/dl)
                    lat1c = yy(y1)
                    wy(y2) = wy(y2)*dy/dl
                    lat2c = yy(y2)
                endif
                if ( npoints == 1 .and. .NOT. gridpoints ) then
                    write(0,'(a,2f9.3,a,2f9.3,a)') 'interpolating points lon=',lon1c &
                        ,lon2c,', lat=',lat1c,lat2c,'<br>'
                endif
                write(out,'(a,2f9.3,a,2f9.3,a)') '# interpolating points lon=',lon1c &
                    ,lon2c,', lat=',lat1c,lat2c
                if ( lwrite ) write(0,'(a,2i4,a,2i4)') &
                    'This corresponds to grid points x=',x1,x2,',y=',y1,y2
            else
                if ( oper == 'v' ) then
                    if ( lon1 == lon2 .and. lat1 == lat2 ) then
                        if ( interp ) then
                            longoper = 'interpolating over'
                            metaoper = 'interp_region'
                        else
                            longoper = 'taking grid box'
                            metaoper = 'ave_region'
                        end if
                    else
                        longoper = 'averaging anomalies over'
                        metaoper = 'ave_region'
                    end if
                else if ( oper == 'x' ) then
                    longoper = 'taking maximum of'
                    metaoper = 'max_region'
                else if ( oper == 'n' ) then
                    longoper = 'taking minimum of'
                    metaoper = 'min_region'
                else
                    longoper = 'cutting out'
                    metaoper = 'region'
                end if                    
                if ( npoints == 1 ) then
                    write(0,'(2a,2f9.3,a,2f9.3,a)') trim(longoper),' region lon=',lon1c &
                        ,lon2c,', lat=',lat1c,lat2c,'<br>'
                endif
                write(out,'(3a,2f9.3,a,2f9.3,a)') '# ',trim(metaoper),' :: lon=',lon1c &
                    ,lon2c,', lat=',lat1c,lat2c
                write(out,'(3a,2f9.3,a,2f9.3,a)') '# ',trim(longoper),' region lon=',lon1c &
                    ,lon2c,', lat=',lat1c,lat2c
                if ( lwrite ) write(0,'(a,2i4,a,2i4)') &
                    'This corrsponds to grid points x=',x1,x2,',y=',y1,y2
            endif
        endif
        call keepalive1('Computed weights',4,5)
        write(out,'(6a)') '# ',trim(vars(1)),' [',trim(newunits),'] ',trim(lvars(1))

!       get mean of just the requested area if only one point is requested
    
        if ( npoints == 1 .and. (missing .or. anom) .and. oper == 'v' ) then
            call getwinmean(mean,nn,nx,ny,nperyear,field,nx,ny &
                ,nperyear,firstyr,lastyr,nx,ny,firstyr,firstmo,nt &
                ,x1,x2,y1,y2,lwrite)
        endif
        call keepalive1('Computed area mean',5,5)
    
!       cut out region
    
        if ( (missing .or. anom) .and. oper == 'v' ) then
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
            if ( month > nperyear ) then
                call outputave(out,year,ave,nperyear,lstandardunits, &
                    offset,slope,ndpm,year,yr1,yr2,nyr)
                month = month - nperyear
                year = year + 1
                do i=1,nperyear
                    ave(i) = 3e33
                enddo
            endif
            if ( month > 0 ) then
                if ( oper == 'v' ) then
                    ave(month) = 0
                else if ( oper == 'n' ) then
                    ave(month) = 3e33
                else if ( oper == 'x' ) then
                    ave(month) = -3e33
                end if
            else
                avemean(month+nperyear) = 0
            endif
            sw = 0
            do j=y1,y2
                do i=x1,x2
                    ii = normx(i,nx)
                    if ( month > 0 ) then
                        f = field(ii,j,month,year)
                    else
                        f = mean(ii,j,month+nperyear)
                    endif
                    if ( lsmasktype == 'land' .or. &
                    lsmasktype == 'notl' ) then
                        if ( abs(lsmask(ii,j)-1) > 1e-4 .eqv. &
                        lsmasktype == 'land' ) then
                            f = 3e33
                            if ( lwrite ) print *,'not ',lsmasktype, &
                            ' point ',ii,j,lsmask(ii,j)
                        endif
                    elseif ( lsmasktype == 'sea ' .or. &
                        lsmasktype == 'nots' ) then
                        if ( abs(lsmask(ii,j)-0) > 1e-4 .eqv. &
                        lsmasktype == 'sea' ) then
                            f = 3e33
                            if ( lwrite ) print *,'not ',lsmasktype, &
                            ' point ',ii,j,lsmask(ii,j)
                        endif
                    elseif ( lsmasktype == '5sea' .or. &
                        lsmasktype == '5lan' ) then
                        if ( lsmask(ii,j) < 0.5 .eqv. &
                        lsmasktype == '5lan' ) then
                            f = 3e33
                            if ( lwrite ) print *,'not ',lsmasktype, &
                            ' point ',ii,j,lsmask(ii,j)
                        endif
                    endif
                    if ( f < -999 .and. f > -1000 ) then
                        write(0,*)'get_index: warning: suspicious ', &
                            'point field(',ii,j,month,year,') = ',f
                    end if
                    if ( f < 1e33 ) then
                        if (  dipole == 'no' .or. &
                        dipole == 'll' .and. i < (x1+x2)/2. .or. &
                        dipole == 'rr' .and. i > (x1+x2)/2. .or. &
                        dipole == 'tt' .and. j < (y1+y2)/2. .or. &
                        dipole == 'bb' .and. j > (y1+y2)/2. .or. &
                        (dipole == 'bl' .or. dipole == 'lb') .and. &
                        (i-x1)*(y1-y2) < (x2-x1)*(j-y2) .or. &
                        (dipole == 'tr' .or. dipole == 'rt') .and. &
                        (i-x1)*(y1-y2) > (x2-x1)*(j-y2) .or. &
                        (dipole == 'br' .or. dipole == 'rb') .and. &
                        (i-x1)*(y2-y1) > (x2-x1)*(j-y1) .or. &
                        (dipole == 'tl' .or. dipole == 'lt') .and. &
                        (i-x1)*(y2-y1) < (x2-x1)*(j-y1) ) then
                            sgn = +1
                        else
                            sgn = -1
                        endif
                        w = wx(ii)*wy(j)
                        if ( lmask ) w = w*mask(ii,j)
                        if ( lwrite .and. .false.) then
                            if ( month > 0 ) then
                                write(0,*) 'adding field(',ii,j &
                                ,month,year,') = ',field(ii,j &
                                ,month,year),mean(ii,j,month) &
                                ,sgn,w
                            else
                                write(0,*) 'adding mean(',ii,j,month+nperyear,') = ', &
                                    mean(ii,j,month+nperyear),sgn,w
                            endif
                        endif
                        if ( month > 0 ) then
                            if ( oper == 'v' ) then
                                ave(month) = ave(month) + sgn*w*field(ii,j,month,year)
                                if ( missing .or. anom ) &
                                    ave(month) = ave(month) - sgn*w*mean(ii,j,month)
                            else if ( oper == 'x' ) then
                                if ( lmask ) then
                                    if ( mask(ii,j).gt.0.5 ) ave(month) = max(ave(month),field(ii,j,month,year))
                                else
                                    ave(month) = max(ave(month),field(ii,j,month,year))
                                end if
                            else if ( oper == 'n' ) then
                                if ( lmask ) then
                                    if ( mask(ii,j).gt.0.5 ) ave(month) = min(ave(month),field(ii,j,month,year))
                                else
                                    ave(month) = min(ave(month),field(ii,j,month,year))
                                end if
                            else
                                write(0,*) 'get_index: error: unrecognised oper ',oper
                                call exit(-1)
                            end if
                        else
                            avemean(month+nperyear) = avemean(month+nperyear) + &
                                sgn*w*mean(ii,j,month+nperyear)
                        endif
                        sw = sw + w
                    else
                        if ( lwrite .and. .false.) print *,'get_index: invalid point'
                    endif
                enddo
            enddo
            if ( month > 0 ) then
                if ( sw < minfac/100*sw0 .or. sw == 0 ) then
                    ave(month) = 3e33
                    if ( lwrite ) write(0,*) 'sw<minfac*area 1: ',sw,minfac/100,sw0,month
                elseif ( missing .and. .NOT. anom .and. oper == 'v' .and. avemean(month) > 1e33 ) then
                    ave(month) = 3e33
                    if ( lwrite ) write(0,*) 'mean missing: ',avemean(month)
                else if ( oper == 'v' ) then
                    ave(month) = ave(month)/sw
                    if ( missing .and. .NOT. anom ) ave(month) = ave(month) + avemean(month)
                    if ( lwrite ) write(*,*) 'ave(',month,') = ',ave(month)
                else if ( oper == 'n' .or. oper == 'x' ) then
                    if ( lwrite ) write(*,*) 'ave(',month,') = ',ave(month)
                endif
            else
                if ( sw == 0 ) then
                    avemean(month+nperyear) = 3e33
                    sw0 = 0
                    if ( lwrite ) write(*,*) 'sw<minfac*area 2: ',sw,minfac/100,month+nperyear
                else
                    sw0 = sw
                    avemean(month+nperyear) = avemean(month+nperyear)/sw
                    if ( lwrite ) write(*,*) 'avemean(',month+nperyear,') = ',avemean(month+nperyear),sw
                endif
            endif
            call keepalive1('Time step',k-k1+1,nt+firstmo-k1)
        enddo               ! loop over all months in file
!       print last (possibly incomplete) record
        call outputave(out,year,ave,nperyear,lstandardunits, &
            offset,slope,ndpm,year,yr1,yr2,nyr)
        if ( nyr == 0 ) then
        !!!write(0,*) 'get_index: error: cannot find any data'
            if ( out /= 6 ) then
                close(out,status='delete')
            end if
        end if
        if ( gridpoints ) then
            if ( nyr > 0 ) then
                print '(a,i4,a,i4,a,i4)','found ',nyr &
                ,' years with data in ',yr1,'-',yr2
            else
                print '(a)','could not locate any data'
            endif
            print '(a)','========='
        endif
!       finito
        if ( npoints > 1 ) then
            close(out)
            ! make a netcdf copy for faster access by next program
            call mysystem('./bin/dat2nc '//trim(string)// &
                ' i get_index '//trim(string(:index(string,'.dat') &
                ))//'nc',iret)
        endif
    900 continue
        call keepalive1('Writing grid points',ipoints,npoints)
    enddo
    goto 999
901 write(0,*) 'get_index: error reading lon1 from ',trim(string)
    call exit(-1)
902 write(0,*) 'get_index: error reading lon2 from ',trim(string)
    call exit(-1)
903 write(0,*) 'get_index: error reading lat1 from ',trim(string)
    call exit(-1)
904 write(0,*) 'get_index: error reading lat2 from ',trim(string)
    call exit(-1)
999 continue
    if ( allocated(lsmask) ) deallocate(lsmask)
end program get_index

subroutine outputave(out,year,ave,nperyear,lstandardunits,offset,slope,ndpm,yr,yr1,yr2,nyr)
    implicit none
    integer :: out,year,nperyear,ndpm,yr,yr1,yr2,nyr
    logical :: lstandardunits
    real :: ave(nperyear),offset,slope
    integer :: i,dy,mo
    integer,save :: dpm(12),init
    integer,external :: leap
    data init /0/
    data dpm / 31,28,31,30,31,30,31,31,30,31,30,31/

    if ( lstandardunits ) then
        do mo=1,nperyear
            if ( ave(mo) < 1e33 ) then
                if ( ndpm /= 0 ) then
                    if ( nperyear /= 12 ) then
                        write(0,*) 'get_index: can only convert for monthly data'
                        call exit(-1)
                    end if
                    if ( mo == 2 .and. leap(yr) == 2 ) then
                        ave(mo) = ave(mo)*29.**ndpm
                    else
                        ave(mo) = ave(mo)*real(dpm(mo))**ndpm
                    endif
                endif
                ave(mo) = ave(mo)*slope + offset
            endif
        enddo
    endif

    do mo=1,nperyear
        if ( ave(mo) < 1e33 ) goto 201
    enddo
    return
201 continue
    yr1 = min(yr1,year)
    nyr = nyr + 1
    yr2 = max(yr2,year)

    call printdatfile(out,ave,nperyear,nperyear,yr,yr)
end subroutine outputave

subroutine getfirstnonzero(mask,nx,ny,axis,sign,y1)
    implicit none
    integer :: nx,ny,sign,y1
    real :: mask(nx,ny)
    character axis
    integer :: i,j,i1,i2,di,j1,j2,dj
    logical :: foundit

    foundit = .false.
    if ( axis == 'y' ) then
        if ( sign == 1 ) then
            j1 = 1
            j2 = ny
            dj = +1
        else if ( sign == -1 ) then
            j1 = ny
            j2 = 1
            dj = -1
        else
            write(0,*) 'getfirstnonzero: error: sign should be +/-1, not ',sign
            write(*,*) 'getfirstnonzero: error: sign should be +/-1, not ',sign
            call exit(-1)
        end if
        do j=j1,j2,dj
            do i=1,nx
                if ( mask(i,j) > 0 ) then
                    foundit = .true.
                    exit
                end if
            end do
            if ( foundit ) exit
        end do
        y1 = j
    else if ( axis == 'x') then
        if ( sign == 1 ) then
            i1 = 1
            i2 = nx
            di = +1
        else if ( sign == -1 ) then
            i1 = nx
            i2 = 1
            di = -1
        else
            write(0,*) 'getfirstnonzero: error: sign should be +/-1, not ',sign,axis
            write(*,*) 'getfirstnonzero: error: sign should be +/-1, not ',sign,axis
            call exit(-1)
        end if
        do i=i1,i2,di
            do j=1,ny
                if ( mask(i,j) > 0 ) then
                    foundit = .true.
                    exit
                end if
            end do
            if ( foundit ) exit
        end do
        y1 = i
    else
        write(0,*) 'getfirstnonzero: error: axis should be x or y, not ',axis
        write(*,*) 'getfirstnonzero: error: axis should be x or y, not ',axis
        call exit(-1)
    end if
    if ( .not. foundit ) then
        write(0,*) 'getfirstnonzero: error: mask is all zero ',axis,sign
        write(*,*) 'getfirstnonzero: error: mask is all zero ',axis,sign
        if ( axis == 'x' ) then
            print *,'i1,i2,di = ',i1,i2,di
        else
            print *,'j1,j2,dj = ',j1,j2,dj
        end if
        call exit(-1)
    end if
     
end subroutine getfirstnonzero
            
