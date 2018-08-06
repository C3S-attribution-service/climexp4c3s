program regionverification

!   program to compute a table containing observations and forecasts
!   for all grid points in a specified region, to be processed
!   further by the R routines

!   (c) Geert Jan van Oldenborgh, 7-nov-2005, KNMI

!   This file may be copied and modified freely as long as the
!   above copyright notice is kept intact

    implicit none
    include 'params.h'
    include 'netcdf.inc'
    include 'getopts.inc'
    integer, parameter :: recfa4=4
    integer :: ncid1,ncid2,nx1,ny1,nz1,nt1,nper1,firstyr1,f1,firstmo1 &
        ,fm1,nx2,ny2,nz2,nt2,nper2,firstyr2,firstmo2,nvars &
        ,jvars1(6,nvmax),jvars2(6,nvmax),endian1,endian2,status &
        ,nxf,nyf,nperyear,firstyr,lastyr,lastyr1,mens1,mens,nmodel1 &
        ,nmodel2,iens,jens,kens,nn,ntmax,month,yr1s,yr2s,it,mineen &
        ,yrstart,yrstop,multimodel(0:nmodelmax),nmodel
    integer, allocatable :: itimeaxis(:)
    integer :: i,j,j1,j2,k,n,jx,jy,mo,m,yr,y,nt,nrec, &
        ivars1(2,6),ivars2(2,6) &
        ,ldir,nx,ny,nz,ncid,iout,x1,x2,y1,y2,mobegin,ntvarid
    real :: xx1(nxmax),yy1(nymax),zz1(nzmax), &
        xx2(nxmax),yy2(nymax),zz2(nzmax), &
        xx(nxmax),yy(nymax),zz(nzmax), &
        wx1(nxmax),wy1(nymax),wx2(nxmax),wy2(nymax), &
        u1,u2
    real,allocatable :: field1(:,:,:,:,:),field2(:,:,:,:), &
        fcst(:,:,:),obs(:,:),newfcst(:,:,:,:,:),newobs(:,:,:,:)
    logical :: lexist,xrev1,xwrap1,yrev1,xrev2,xwrap2,yrev2,lvalid
    character infile*256,datfile1*256,datfile2*256, &
        title1*1024,oldtitle1*1024,title2*1024,vars1(1)*20, &
        lvars1(1)*200,units1(1)*40,vars2(1)*20,lvars2(1)*200, &
        units2(1)*40,file*256,tmpunits*20,oldmodel*20,newmodel*20
    character line*80

!   check arguments

    lwrite = .FALSE. 
    n = command_argument_count()
    if ( n < 4 ) then
        print *,'usage: regionverification '// &
        'field1.[ctl|nc] field2.[ctl|nc] '// &
        ' [sum|ave|max|min|sel n] '// &
        '[begin yr] [end yr] dump table|plot file.nc'
        stop
    endif
    call killfile(infile,datfile1,datfile2,0)
    call get_command_argument(1,infile)
    mens1 = 0
    if ( index(infile,'%') > 0 .OR. index(infile,'++') > 0 ) then
        ensemble = .TRUE. 
        call filloutens(infile,0)
        if ( lwrite ) print *,'check whether ',trim(infile) &
        ,' exists'
        inquire(file=infile,exist=lexist)
        if ( .NOT. lexist ) then
            if ( lwrite ) print *,'no, start at 1'
            call get_command_argument(1,infile)
            mens1 = 1
            call filloutens(infile,1)
        endif
        do mens=mens1+1,nensmax
            call get_command_argument(1,file)
            if ( lwrite ) print *,'check whether ',trim(infile) &
            ,' exists'
            call filloutens(file,mens)
            inquire(file=file,exist=lexist)
            if ( .NOT. lexist ) exit
        enddo
        mens = mens - 1
        print '(a,i2,a,i2)','# found ensemble members ',mens1,' to ' &
        ,mens
    else
        ensemble = .FALSE. 
        mens1 = 0
        mens = 0
    endif
    if ( lwrite ) print *,'regionverification: nf_opening file ',trim(infile)
    status = nf_open(infile,nf_nowrite,ncid1)
    if ( status /= nf_noerr ) then
        call parsectl(infile,datfile1,nxmax,nx1,xx1,nymax,ny1,yy1 &
        ,nzmax,nz1,zz1,nt1,nper1,firstyr1,firstmo1,u1,endian1 &
        ,title1,1,nvars,vars1,ivars1,lvars1,units1)
        ncid1 = -1
    else
        datfile1 = infile
        call parsenc(infile,ncid1,nxmax,nx1,xx1,nymax,ny1,yy1 &
        ,nzmax,nz1,zz1,nt1,nper1,firstyr1,firstmo1,u1,title1,1 &
        ,nvars,vars1,jvars1,lvars1,units1)
    endif

    call get_command_argument(2,infile)
    if ( infile == 'perfectmodel' ) then
        nmodel1 = mens1
        nmodel2 = mens
        vars2   = vars1
        ivars2  = ivars1
        lvars2  = lvars2
        units2  = units1
    else
        nmodel1 = -1
        nmodel2 = -1
        if ( lwrite ) print *,'regionverification: nf_opening file ',trim(infile)
        status = nf_open(infile,nf_nowrite,ncid2)
        if ( status /= nf_noerr ) then
            call parsectl(infile,datfile2,nxmax,nx2,xx2,nymax,ny2 &
            ,yy2,nzmax,nz2,zz2,nt2,nper2,firstyr2,firstmo2,u2 &
            ,endian2,title2,1,nvars,vars2,ivars2,lvars2,units2)
            ncid2 = -1
        else
            datfile2 = infile
            call parsenc(infile,ncid2,nxmax,nx2,xx2,nymax,ny2,yy2 &
            ,nzmax,nz2,zz2,nt2,nper2,firstyr2,firstmo2,u2 &
            ,title2,1,nvars,vars2,jvars2,lvars2,units2)
        endif
    endif

    if ( nmodel1 == -1 ) then
        nxf = max(nx1,nx2)
        nyf = max(ny1,ny2)
        firstyr = max(firstyr1,firstyr2)
        nperyear = max(nper1,nper2)
        lastyr1 = firstyr1 + (nt1+firstmo1-2)/nperyear
        lastyr = min(lastyr1,firstyr2 + (nt2+firstmo2-2)/nperyear)
    else
        nxf = nx1
        nyf = ny1
        firstyr = firstyr1
        nperyear = nper1
        lastyr1 = firstyr1 + (nt1+firstmo1-2)/nperyear
        lastyr = lastyr1
    endif
    nt = nperyear*(lastyr-firstyr+1)
    call getopts(3,command_argument_count(),nperyear,firstyr,lastyr, .TRUE. , &
    mens1,mens)
    if ( lag1 /= 0 .OR. lag2 /= 0 ) print *,'ignoring lags'
    if ( .NOT. dump .AND. .NOT. plot ) write(0,*) &
    'regionverification: specify table'
    if ( lks ) write(0,*) 'regionverification: K-S not supported'
    if ( mdiff /= 0 .OR. mdiff2 /= 0 ) write(0,*) &
    'regionverification: monthly anomalies not supported'
    if ( lconting ) write(0,*) 'regionverification: contingency '// &
    'tables not supported'
    do i=1,indxuse
        if ( lincl(i) ) write(0,*) 'regionverification: what do ', &
        'you mean with ',strindx(i),'?'
    enddo
    yr1 = max(yr1,firstyr)
    yr2 = min(yr2,lastyr)
    firstyr = yr1
    lastyr  = yr2
    if ( lwrite ) print *,'yr1,yr2 = ',yr1,yr2
    if ( nmodel1 >= 0 ) then
        nmodel1 = nens1
        nmodel2 = nens2
    endif
    if ( lwrite ) print *,'nens1,nens2 = ',nens1,nens2

!   pgf90 on RHEL happily allocates >2GB :-(
    if ( nxf*nyf*nperyear*(lastyr1-firstyr1+1)*(nens2+1) > 2**30 ) &
    then
        write(0,*) 'error: array too large, reduce resolution, ' &
        ,'number of years or number of ensemble members ',nxf &
        *nyf*nperyear*(lastyr1-firstyr1+1)*(nens2+1),'>',2**30
        write(*,*) 'error: array too large, reduce resolution, ' &
        ,'number of years or number of ensemble members ',nxf &
        *nyf*nperyear*(lastyr1-firstyr1+1)*(nens2+1),'>',2**30
        print *,'nxf      = ',nxf
        print *,'nyf      = ',nyf
        print *,'nperyear = ',nperyear
        print *,'firstyr1 = ',firstyr1
        print *,'lastyr1  = ',lastyr1
        print *,'nens2    = ',nens2
        call abort
    endif

    allocate(field1(nxf,nyf,nperyear,firstyr1:lastyr1,0:nens2))
    allocate(fcst(nperyear,firstyr:lastyr,0:nens2))
    if ( nmodel1 == -1 ) then
        allocate(field2(nxf,nyf,nperyear,firstyr:lastyr))
        allocate(obs(nperyear,firstyr:lastyr))
    endif

!   read fields

    multimodel = 0
    nmodel = 0
    oldmodel = ' '
    oldtitle1 = ' '
    if ( ncid1 == -1 ) then
        if ( ensemble ) then
            do iens=nens1,nens2
                call get_command_argument(1,infile)
                call filloutens(infile,iens)
                call parsectl(infile,datfile1,nxmax,nx1,xx1,nymax &
                ,ny1,yy1,nzmax,nz1,zz1,nt1,nper1,f1,fm1 &
                ,u1,endian1,title1,1,nvars,vars1 &
                ,ivars1,lvars1,units1)
                call datfile2modelname(datfile1,newmodel)
                if ( newmodel /= oldmodel ) then
                    oldmodel = newmodel
                    multimodel(nmodel) = iens
                    nmodel = nmodel + 1
                endif
                call keepalive2('Reading ensemble member ',iens, &
                nens2-nens1+1, .TRUE. )
                call readdatfile(datfile1 &
                ,field1(1,1,1,firstyr1,iens) &
                ,nxf,nyf,nx1,ny1,nperyear,firstyr1,lastyr1 &
                ,f1,fm1,nt1,u1,endian1,lwrite,yr1,yr2,1,1)
            enddo
        else
            call keepalive(1,2)
            call readdatfile(datfile1,field1,nxf,nyf,nx1,ny1 &
            ,nperyear,firstyr1,lastyr1,firstyr1,firstmo1,nt1,u1 &
            ,endian1,lwrite,max(firstyr1,yr1),min(lastyr1,yr2) &
            ,1,1)
        endif
    else
        if ( ensemble ) then
            do iens=nens1,nens2
                call get_command_argument(1,infile)
                call filloutens(infile,iens)
                status = nf_open(infile,nf_nowrite,ncid1)
                call parsenc(infile,ncid1,nxmax,nx1,xx1,nymax,ny1 &
                ,yy1,nzmax,nz1,zz1,nt1,nper1,f1,fm1 &
                ,u1,title1,1,nvars,vars1,jvars1,lvars1,units1)
                if ( title1 /= oldtitle1 ) then
                    oldtitle1 = title1
                    multimodel(nmodel) = iens
                    nmodel = nmodel + 1
                endif
                call keepalive2('Reading ensemble member ',iens, &
                nens2-nens1+1, .TRUE. )
                call readncfile(ncid1,field1(1,1,1,firstyr,iens),nxf &
                ,nyf,nx1,ny1,nperyear,firstyr1,lastyr1 &
                ,f1,fm1,nt1,u1,lwrite &
                ,max(f1,yr1),min(lastyr1,yr2),jvars1)
            enddo
        else
            call keepalive(1,2)
            call readncfile(ncid1,field1,nxf,nyf,nx1,ny1,nperyear &
            ,firstyr1,lastyr1,firstyr1,firstmo1,nt1,u1,lwrite &
            ,max(firstyr1,yr1),min(lastyr1,yr2),jvars1)
        endif
    endif
    multimodel(nmodel) = nens2+1
    print '(a,100i3)','# nmodel, multimodel = ',nmodel, &
    multimodel(0:nmodel)
!       convert to standard units
    do iens=nens1,nens2
        tmpunits = units1(1) ! they are otherwise adjusted
        call makestandardfield(field1(1,1,1,firstyr1,iens),nxf,nyf,1 &
        ,nperyear,firstyr1,lastyr1,nx1,ny1,1,nperyear &
        ,firstyr1,lastyr1,vars1(1),tmpunits,lwrite)
    enddo
    units1(1) = tmpunits
    call keepalive(nens2-nens1+2,nens2-nens1+2)
    if ( nmodel1 == -1 ) then
        if ( ncid2 == -1 ) then
            call readdatfile(datfile2,field2,nxf,nyf,nx2,ny2 &
            ,nperyear,firstyr,lastyr,firstyr2,firstmo2,nt2,u2 &
            ,endian2,lwrite,yr1,yr2,1,1)
        else
            call readncfile(ncid2,field2,nxf,nyf,nx2,ny2,nperyear &
            ,firstyr,lastyr,firstyr2,firstmo2,nt2,u2,lwrite,yr1 &
            ,yr2,jvars2)
        endif
    !           convert to standard units
        call makestandardfield(field2,nxf,nyf,1,nperyear,firstyr &
        ,lastyr,nx2,ny2,1,nperyear,yr1,yr2,vars2(1) &
        ,units2(1),lwrite)
        call keepalive(0,0)
    endif

!       cut out region of interest

    call getxyprop(xx1,nx1,yy1,ny1,xrev1,yrev1,xwrap1)
    call getlatlonwindow(lat1,lat2,lon1,lon2,xx1,nx1,xwrap1,1,yy1 &
    ,ny1,1,x1,x2,y1,y2,lwrite)
    if ( .NOT. dump ) then
        call biggerwindow(xx1,nx1,yy1,ny1,xwrap1,lon1,lon2,lat1,lat2 &
        ,x1,x2,y1,y2)
    endif
    call enscutoutwindow(x1,x2,y1,y2,xx1,nx1,xwrap1,xrev1,1,yy1,ny1 &
    ,1,wx1,wy1,field1,nxf,nyf,nens1,nens2,nperyear,firstyr1 &
    ,lastyr1,max(firstyr1,yr1),min(lastyr1,yr2),lwrite)
    if ( nmodel1 == -1 ) then
        call getxyprop(xx2,nx2,yy2,ny2,xrev2,yrev2,xwrap2)
        call getlatlonwindow(lat1,lat2,lon1,lon2,xx2,nx2,xwrap2,1 &
        ,yy2,ny2,1,x1,x2,y1,y2,lwrite)
        if ( .NOT. dump ) then
            call biggerwindow(xx2,nx2,yy2,ny2,xwrap2,lon1,lon2,lat1 &
            ,lat2,x1,x2,y1,y2)
        endif
        call enscutoutwindow(x1,x2,y1,y2,xx2,nx2,xwrap2,xrev2,1,yy2 &
        ,ny2,1,wx2,wy2,field2,nxf,nyf,0,0,nperyear,firstyr &
        ,yrend,yr1,yr2,lwrite)
    
    !           interpolate fields to common grid
    
        call ensxyinterpu( &
        field1,xx1,nx1,yy1,ny1,nens1,nens2, &
        field2,xx2,nx2,yy2,ny2,0,0, &
        xx,nx,yy,ny,firstyr1,lastyr1,firstyr,lastyr &
        ,nxf,nyf,1,nz1,nperyear,intertype,lwrite)
    else
        nx = nx1
        ny = ny1
        xx(1:nx) = xx1(1:nx)
        yy(1:ny) = yy1(1:ny)
    endif
!       my routines may have shifted the longitude axis - shift back
    if ( abs(xx(1)-lon1) > abs(xx(1)-lon1-360) ) then
        if ( abs(xx(nx)-lon2) > abs(xx(nx)-lon1-360) ) then
            xx = xx - 360
        endif
    endif
    if ( dump .AND. plot ) then
        write(11,'(a,4f8.2)') '#',lon1,lon2,lat1,lat2
        do j=1,ny
            do i=1,nx
                write(11,'(2i4.4,4f9.3,a,i4.4,a,i4.4)') &
                i,j,xx(i),yy(j),0.2,-1.,' point_',i,'_',j
            enddo
        enddo
    endif

!**        print '(a,f7.3,a,f7.3,a,f8.3,a,f8.3)','# lat = ',yy(1),' to '
!**  +       ,yy(ny),', lon = ',xx(1),' to ',xx(nx)
    nn = 0
    if ( .NOT. dump ) then
        if ( m1 == m2 .AND. nperyear <= 12 ) then
            if ( nmodel1 == -1 ) then
                allocate(newobs(nx,ny,1,yr1:yr2))
            endif
            allocate(newfcst(nx,ny,1,yr1:yr2,nens1:nens2))
        else
            if ( nmodel1 == -1 ) then
                allocate(newobs(nx,ny,nperyear,yr1:yr2))
            endif
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
        
            ! create timeseries from fields
        
            fcst(1:nperyear,firstyr:max(yr1-1,firstyr1-1),0:nens2) = &
            &                3e33
            fcst(1:nperyear,min(yr2+1,lastyr1+1):lastyr,0:nens2) = &
            &                3e33
            obs(1:nperyear,firstyr:yr1-1) = 3e33
            obs(1:nperyear,yr2+1:lastyr) = 3e33
            n = 0
            if ( nmodel1 == -1 ) then
                do yr=yr1,yr2
                    do mo=1,nperyear
                        obs(mo,yr) = field2(jx,jy,mo,yr)
                    enddo
                enddo
            endif
            do yr=yr1,yr2
                do mo=1,nperyear
                    do iens=nens1,nens2
                        if ( yr < firstyr1 .OR. yr > lastyr1 ) &
                        then
                            fcst(mo,yr,iens) = 3e33
                        else
                            fcst(mo,yr,iens) = &
                            field1(jx,jy,mo,yr,iens)
                        endif
                    enddo
                enddo       ! mo
            enddo           ! yr
        
        !               transform time series
            if ( lwrite ) then
                print *,'@@@ ',jx,xx(jx),jy,yy(jy)
                write(*,'(a)') '@@@ fcst'
                do yr=yr1,yr2
                    write(*,'(i4,12f8.1)') yr,(fcst(mo,yr,nens1),mo &
                    =1,12)
                end do
                write(*,'(a)') '@@@ obs'
                do yr=yr1,yr2
                    write(*,'(i4,12f8.1)') yr,(obs(mo,yr),mo=1,12)
                end do
            end if
            call manipulatetimeseries(fcst,obs,nperyear,firstyr &
            ,lastyr,nperyear,nmodel1,vars2(1),jx,jy, &
            multimodel,nmodel)
        
            ! print table if requested
        
            if ( dump ) then
                call printtable(fcst,obs,nperyear,firstyr,lastyr, &
                nmodel1,nmodel2,nens1,nens2,yr1,yr2,m1,m2 &
                ,nperyear,nn,yrstart,yrstop)
            else
                ! were "adjusted" by manipulatetimeseries :-(
                yr1 = yr1s
                yr2 = yr2s
                do yr=yr1,yr2
                    do month=m1,m2
                        call getj1j2(j1,j2,month,nperyear, .FALSE. )
                        do j=j1,j2
                            if ( m1 == m2 .AND. nperyear <= 12 ) &
                            then
                                mo = 1
                            else
                                mo = j
                            endif
                            if ( nmodel1 == -1 ) then
                                newobs(jx,jy,mo,yr) = obs(j,yr)
                            endif
                            do iens=nens1,nens2
                                newfcst(jx,jy,mo,yr,iens) = &
                                fcst(j,yr,iens)
                                if ( fcst(j,yr,iens) < 1e33 ) then
                                    yrstart = min(yrstart,yr)
                                    yrstop  = max(yrstop,yr)
                                endif
                            enddo ! j
                        enddo ! iens
                    enddo   ! month
                enddo       ! yr
            endif           ! dump
        enddo               ! jx
    enddo                   ! jy
    if ( dump ) then
        if ( nn == 0 ) then
            write(*,*) &
            'regionverification: did not find any valid data'
            close(10,status='delete')
        else
            close(10)
        endif
    else                    ! .NOT. dump
    
        ! write the field to netcdf
    
        if ( lwrite ) print *,'yr1,yr2 = ',yr1,yr2
        nt = (yr2-yr1+1)
        if ( m1 /= m2 .OR. nperyear > 12 ) nt = nt*nperyear
        if ( nmodel1 /= -1 ) then
            nt = nt*(nens2-nens1+1)
            mineen = -1
        else
            mineen = 0
        endif
        allocate(itimeaxis(nt))
        if ( nz1 > 1 .OR. nz2 > 1 ) then
            write(0,*) 'regionverification: error: cannot handle '// &
            'nz>1 yet: ',nz1,nz2
            call abort
        else
            nz = 1
            ivars1(1,1) = 0
            ivars2(1,1) = 0
        endif
        if ( m1 /= m2 .OR. nperyear > 12 ) then
            mobegin = 1
        else
            mobegin = m1
            nperyear = 1
        endif
    
        ! forecasts
    
        call writenc(plotfile,ncid1,ntvarid,itimeaxis,nt,nx,xx,ny &
        ,yy,nz,zz,nt,nperyear,yr1,mobegin,3e33,title1,nvars &
        ,vars1,ivars1,lvars1,units1,nens1,nens2+mineen)
        if ( nmodel1 /= -1 ) then
        
            ! perfect model
        
            do iens=nens1,nens2-1
                do yr=yr1,yr2
                    do jens=nens1,nens2
                        if ( iens < jens ) then
                            kens = iens
                        else
                            kens = iens + 1
                        endif
                        if ( yr == yr1 ) print '(a,i3,a,i3,a,i4)' &
                        ,'# writing newfcst(',kens &
                        ,') to forecast file as member ' &
                        ,iens-nens1+1,' at yr ',yr+(yr2-yr1+1) &
                        *(jens-nens1)
                        if ( m1 == m2 .AND. nperyear <= 12 ) then
                            call writencslice(ncid1,ntvarid &
                            ,itimeaxis,nt,ivars1 &
                            ,newfcst(1,1,1,yr,kens) &
                            ,nx,ny,nz,nx,ny,nz &
                            ,yr-yr1+1+(yr2-yr1+1)*(jens-nens1) &
                            ,iens-nens1+1)
                        else
                            do mo=1,nperyear
                                call writencslice(ncid1,ntvarid &
                                ,itimeaxis,nt,ivars1 &
                                ,newfcst(1,1,mo,yr,kens) &
                                ,nx,ny,nz,nx,ny,nz &
                                ,mo+nperyear* &
                                (yr-yr1+(yr2-yr1+1)*(jens-nens1 &
                                )),iens-nens1+1)
                            enddo
                        endif
                    enddo   ! jens
                enddo       ! yr
            enddo           ! iens
        else
        
            ! normal case
        
            do iens=nens1,nens2
                do yr=yr1,yr2
                    if ( m1 == m2 .AND. nperyear <= 12 ) then
                        call writencslice(ncid1,ntvarid,itimeaxis,nt &
                        ,ivars1,newfcst(1,1,1,yr,iens),nx,ny,nz &
                        ,nx,ny,nz,yr-yr1+1,iens-nens1+1)
                    else
                        do mo=1,nperyear
                            call writencslice(ncid1,ntvarid &
                            ,itimeaxis,nt,ivars1,newfcst(1,1,mo &
                            ,yr,iens),nx,ny,nz,nx,ny,nz,mo &
                            +nperyear*(yr-yr1),iens-nens1+1)
                        enddo
                    endif
                enddo
            enddo
        endif
    
        ! observations
    
        i = index(plotfile,'.nc')
        if ( i == 0 ) i = len_trim(plotfile) + 1
        plotfile(i:) = '_obs.nc'
        call writenc(plotfile,ncid2,ntvarid,itimeaxis,nt,nx,xx,ny &
        ,yy,nz,zz,nt,nperyear,yr1,mobegin,3e33,title2,nvars &
        ,vars2,ivars2,lvars2,units2,0,0)
        if ( nmodel1 /= -1 ) then
        
            ! perfect model
        
            do yr=yr1,yr2
                do iens=nens1,nens2
                    if ( m1 == m2 .AND. nperyear <= 12 ) then
                        if ( yr == yr1 ) print '(a,i3,a,i4)' &
                        ,'# writing newfcst(',iens &
                        ,') to observation file at yr ', &
                        yr+(yr2-yr1+1)*(iens-nens1)
                        call writencslice(ncid2,ntvarid,itimeaxis,nt &
                        ,ivars2,newfcst(1,1,1,yr,iens) &
                        ,nx,ny,nz,nx,ny,nz &
                        ,yr-yr1+1+(yr2-yr1+1)*(iens-nens1),1)
                    else
                        do mo=1,nperyear
                            call writencslice(ncid2,ntvarid &
                            ,itimeaxis,nt,ivars2 &
                            ,newfcst(1,1,mo,yr,iens) &
                            ,nx,ny,nz,nx,ny,nz &
                            ,mo+nperyear &
                            *(yr-yr1+(yr2-yr1+1)*(iens-nens1)) &
                            ,1)
                        enddo
                    endif
                enddo
            enddo
        else
        
            ! normal case
        
            do yr=yr1,yr2
                if ( m1 == m2 .AND. nperyear <= 12 ) then
                    call writencslice(ncid2,ntvarid,itimeaxis,nt &
                    ,ivars2,newobs(1,1,1,yr),nx,ny,nz,nx,ny,nz &
                    ,yr-yr1+1,1)
                else
                    do mo=1,nperyear
                        call writencslice(ncid2,ntvarid,itimeaxis,nt &
                        ,ivars2,newobs(1,1,mo,yr),nx,ny,nz,nx &
                        ,ny,nz,mo+nperyear*(yr-yr1),1)
                    enddo
                endif
            enddo
        endif
        ! essential!
        i = nf_close(ncid1)
        i = nf_close(ncid2)
    endif

!   and export yrstart, yrstop with interval with valid data

    call savestartstop(yrstart,yrstop)
end program
