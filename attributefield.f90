    program attributefield
!
!   do an empirical attribution study by fitting a time series to a
!   GEV or GPD with the position parameter linearly dependent on a covariate
!   and studying the difference in return time in the current climate and a
!   previous climate.
!
    use lsdata
    implicit none
    include 'params.h'
    include 'getopts.inc'
    integer nvarmax
    parameter(nvarmax=1)
    integer nperyear,nperyear1,mens1,mens,ncid,iens,firstyr,firstmo,endian
    integer i,yr,mo,n,j1,j2,off,nx,ny,nz,nt,nvars,fyr,lyr,nxf,nyf,nzf
        integer x1,x2,y1,y1
    real xx(nxmax),yy(nymax),zz(nzmax),undef,wx(nxmax),wy(nymax)
    real,allocatable :: field(:,:,:,:,:,:),covariate(:,:,:),series(:,:,:)
    logical xrev,yrev,xwrap
    character file*1024,datfile*1024,covariatefile*1024,distribution*6,assume*5,string*80
    character var*40,units*80,var1*40,units1*80,lz(3)*20,ltime*120,title*255,history*2048
    character vars(nvarmax)*20,svars(nvarmax)*20,lvars(nvarmax)*120,units(nvarmax)*40
    character cell_methods(nvarmax)*100
    integer iargc
!
    if ( iargc().lt.8 ) then
        write(0,*) 'usage: attributefield field covariate_series ', &
        & 'GEV|Gumbel|GPD|Gauss assume shift|scale', &
        & 'mon n [sel m] [ave N] [log|sqrt] ', &
        & 'begin2 past_climate_year end2 year_under_study', &
        & 'plot FAR_plot_file [dgt threshold%]'
        write(0,*) 'note that n and m are in months even if the series is daily.'
        write(0,*) 'N is always in the same units as the series.'
        write(0,*) 'the covariate series is averaged to the same time scale.'
        stop
    end if
!
!   initialisation
!
    call attribute_init(file,distribution,assume,off,nperyear,yrbeg,yrend,nensmax,lwrite)
    call getmetadata(file,mens1,mens,ncid,datfile,nxmax,nx &
     &       ,xx,nymax,ny,yy,nzmax,nz,zz,lz,nt,nperyear,firstyr,firstmo &
     &       ,ltime,undef,endian,title,history,nvarmax,nvars,vars,jvars &
     &       ,lvars,svars,units,cell_methods,lwrite)
!
    nxf = nx
    nyf = ny
    nzf = nz
    fyr = firstyr
    lyr = firstyr + (nt-1)/nperyear
    allocate(field(nxf,nyf,nzf,nperyear,fyr:lyr,0:mens))
    call readfield(ncid,file,datfile,field,nxf,nyf,nzf &
     &       ,nperyear,fyr,lyr,nens1,nens2,nx,ny,nz,nperyear,yr1,yr2 &
     &       ,firstyr,firstmo,nt,undef,endian,vars,units,lstandardunits &
     &       ,lwrite)
    
    call getarg(2+off,covariatefile)
    allocate(covariate(npermax,yrbeg:yrend,0:nensmax))
    if ( index(covariatefile,'%%') == 0 .and. &
    &    index(covariatefile,'++') == 0 ) then
        call readseries(covariatefile,covariate,npermax,yrbeg,yrend &
        & ,nperyear1,var1,units1,lstandardunits,lwrite)
        do iens=mens1+1,mens
            covariate(:,:,iens) = covariate(:,:,mens1)
        end do
    else
        call readensseries(covariatefile,covariate,npermax,yrbeg,yrend &
        & ,nperyear1,var1,units1,lstandardunits,lwrite)
    end if
    
    call getopts(6+off,iargc(),nperyear,yrbeg,yrend,.true.,mens1,mens)
    if ( yr1a.lt.yr1 .or. yr1a.lt.yrbeg ) then
        write(0,*) 'attribute: error: reference year should be after start of series ',yr1,yr1a
        call abort
    end if
    if ( yr2a.gt.yr2 .or. yr2a.gt.yrend ) then
        write(0,*) 'attribute: error: current year should be before end of series ',yr2,yr2a
        call abort
    end if
    call getj1j2(j1,j2,m1,nperyear,lwrite)
    call print_bootstrap_message(max(1,nint(decor)),j1,j2)
    if ( distribution.eq.'gpd' ) print '(a)','# after declustering.'
!
!	get boundaries in grid points
!
    call getxyprop(xx,nx,yy,ny,xrev,yrev,xwrap)
    call getlatlonwindow(lat1,lat2,lon1,lon2,xx,nx,xwrap,avex,yy,ny &
 &        ,avey,x1,x2,y1,y2,lwrite)
!
!       average, cut out window - everything to make the arrays smaller
!
    call keepalive1('Cutting out region',0,-1)
    call enscutoutwindow(x1,x2,y1,y2,xx,nx,xwrap,xrev,avex,yy,ny, &
 &       avey,wx,wy,field,nxf,nyf,nens1,nens2,nperyear, &
 &       fyr,lyr,yr1,yr2,lwrite)
!
!   process data
!
    do iy=y1,y2
        do ixx=x1,x2
            ix = ixx
            if ( ix.gt.nx ) ix = ix - nx
            if ( ix.lt.1 ) ix = ix + nx
            do iens=nens1,nens2
                do i=yr1,yr2
                    do j=1,nperyear
                        series(j,i,iens) = field(ix,iy,j,i,iens)
                    end do
                end do
            end do
            
    if ( lchangesign ) then
        do iens=mens1,mens
            call changesign(series(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend)
        end do
    endif
    if ( ldetrend ) then
        do iens=mens1,mens
            call detrend(series(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend,yr1,yr2,m1,m2,lsel)
            call detrend(covariate(1,yrbeg,iens),npermax,nperyear1,yrbeg,yrend,yr1,yr2,m1,m2,lsel)
        end do
    end if
    if ( anom ) then
        do iens=mens1,mens
            call anomal(series(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend,yr1,yr2)
        end do
    end if
    if ( lsum.gt.1 ) then
        do iens=mens1,mens
            call sumit(series(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend,lsum,oper)
        end do
    endif
    if ( logscale ) then
        print '(a)','# taking logarithm'
        do iens=mens1,mens
            call takelog(series(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend)
        end do
        if ( xyear.lt.1e33 ) xyear = log(xyear)
    endif
    if ( sqrtscale ) then
        print '(a)','# taking sqrt'
        do iens=mens1,mens
            call takesqrt(series(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend)
        end do
        if ( xyear.lt.1e33 ) xyear = sqrt(xyear)
    endif
    if ( squarescale ) then
        print '(a)','# taking square'
        do iens=mens1,mens
            call takesquare(series(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend)
        end do
        if ( xyear.lt.1e33 ) xyear = xyear**2
    endif

    call attribute_dist(series,nperyear,covariate,nperyear1,npermax,yrbeg,yrend,&
    &   mens1,mens,assume,distribution)

end program attributefield

