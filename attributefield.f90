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
    include 'netcdf.inc'
    integer nvarmax,nresmax,ntmax
    parameter(nvarmax=200,nresmax=100,ntmax=1)
    integer nperyear,nperyear1,mens1,mens,ncid,iens,firstyr,firstmo,endian
    integer i,j,k,l,yr,mo,n,j1,j2,off,nx,ny,nz,nt,nvars,fyr,lyr,nxf,nyf,nzf,nresults
    integer x1,x2,y1,y2,ix,ixx,iy,iz,iyrs(10),ntvarid,itimeaxis(ntmax),ivars(2,nvarmax), &
    &   jvars(6,nvarmax)
    real xx(nxmax),yy(nymax),zz(nzmax),undef,wx(nxmax),wy(nymax),results(3,nresmax)
    real,allocatable :: field(:,:,:,:,:,:),covariate(:,:,:),series(:,:,:),res(:,:,:,:,:)
    logical xrev,yrev,xwrap,lfirst,lprint
    character file*1024,datfile*1024,covariatefile*1024,distribution*6,assume*5,string*80
    character var*40,var1*40,units1*80,lz(3)*20,ltime*120,title*255,history*2048
    character vars(nvarmax)*20,svars(nvarmax)*20,lvars(nvarmax)*120,units(nvarmax)*40
    character cell_methods(nvarmax)*100,orgunits*40,outfile*1024,format*20,seriesids(0:nensmax)*10
    integer iargc
    data iyrs /10,20,50,100,200,500,1000,2000,5000,10000/
!
    if ( iargc().lt.8 ) then
        write(0,*) 'usage: attributefield field covariate_series ', &
        & 'GEV|Gumbel|GPD|Gauss assume shift|scale ', &
        & 'mon n [sel m] [ave N] [log|sqrt] ', &
        & 'begin2 past_climate_year end2 year_under_study ', &
        & '[dgt threshold%] outfile.nc'
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
    yr1 = max(yr1,fyr)
    yr2 = min(lyr,lyr)
    nens1 = max(nens1,mens1)
    nens2 = min(nens2,mens)
    allocate(field(nxf,nyf,nzf,nperyear,fyr:lyr,0:mens))
    call readfield(ncid,file,datfile,field,nxf,nyf,nzf &
     &       ,nperyear,fyr,lyr,nens1,nens2,nx,ny,nz,nperyear,yr1,yr2 &
     &       ,firstyr,firstmo,nt,undef,endian,vars,orgunits,lstandardunits &
     &       ,lwrite)
    
    call getarg(2+off,covariatefile)
    allocate(covariate(npermax,fyr:lyr,0:nensmax))
    if ( index(covariatefile,'%%') == 0 .and. &
    &    index(covariatefile,'++') == 0 ) then
        call readseries(covariatefile,covariate,npermax,fyr,lyr &
        & ,nperyear1,var1,units1,lstandardunits,lwrite)
        do iens=mens1+1,mens
            covariate(:,:,iens) = covariate(:,:,mens1)
        end do
    else
        call readensseries(covariatefile,covariate,npermax,fyr,lyr &
        & ,nperyear1,var1,units1,lstandardunits,lwrite)
    end if
    
    call getopts(6+off,iargc()-1,nperyear,yrbeg,yrend,.true.,mens1,mens)
    yr1 = max(yr1,fyr) ! messed up by getopts :-(
    yr2 = min(yr2,lyr)
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
    lfirst = .true.
    iz = 1
    do iy=y1,y2
        do ixx=x1,x2
            ix = ixx
            if ( ix.gt.nx ) ix = ix - nx
            if ( ix.lt.1 ) ix = ix + nx
            if ( lfirst ) allocate(series(npermax,fyr:lyr,0:mens))
            do iens=nens1,nens2
                do i=yr1,yr2
                    do j=1,nperyear
                        series(j,i,iens) = field(ix,iy,iz,j,i,iens)
                    end do
                end do
            end do
            
            if ( lchangesign ) then
                do iens=mens1,mens
                    call changesign(series(1,fyr,iens),npermax,nperyear,fyr,lyr)
                end do
            end if
            if ( ldetrend ) then
                do iens=mens1,mens
                    call detrend(series(1,fyr,iens),npermax,nperyear,fyr,lyr,yr1,yr2,m1,m2,lsel)
                    if ( lfirst) call detrend(covariate(1,fyr,iens),npermax,nperyear1,fyr,lyr,yr1,yr2,m1,m2,lsel)
                end do
            end if
            if ( anom ) then
                do iens=mens1,mens
                    call anomal(series(1,fyr,iens),npermax,nperyear,fyr,lyr,yr1,yr2)
                end do
            end if
            if ( lsum.gt.1 ) then
                do iens=mens1,mens
                    call sumit(series(1,fyr,iens),npermax,nperyear,fyr,lyr,lsum,oper)
                end do
            endif
            if ( logscale ) then
                print '(a)','# taking logarithm'
                do iens=mens1,mens
                    call takelog(series(1,fyr,iens),npermax,nperyear,fyr,lyr)
                end do
                if ( xyear.lt.1e33 ) xyear = log(xyear)
            endif
            if ( sqrtscale ) then
                print '(a)','# taking sqrt'
                do iens=mens1,mens
                    call takesqrt(series(1,fyr,iens),npermax,nperyear,fyr,lyr)
                end do
                if ( xyear.lt.1e33 ) xyear = sqrt(xyear)
            endif
            if ( squarescale ) then
                print '(a)','# taking square'
                do iens=mens1,mens
                    call takesquare(series(1,fyr,iens),npermax,nperyear,fyr,lyr)
                end do
                if ( xyear.lt.1e33 ) xyear = xyear**2
            endif

            if ( lwrite ) print *,'attributefield: calling attribute_dist for ix,iy,iz = ', &
            & ix,iy,iz,xx(ix),yy(iy),zz(iz)
            call keepalive1('Grid point',ix+(iy-1)*ny+(nz-1)*nx*ny,nx*ny*nz)
            lprint = .false.
            do i=mens1,mens
                write(seriesids(i),'(i3.3)') i
            enddo
            call attribute_dist(series,nperyear,covariate,nperyear1,npermax,fyr,lyr,&
            &   mens1,mens,assume,distribution,seriesids,results,nresmax,nresults,lprint)
            if ( lwrite ) then
                print *,'results = '
                do i=1,nresults
                    print *,(results(j,i),j=1,3)
                    if ( results(1,i) == 3e33 ) exit
                end do
            end if
            if ( lfirst ) allocate(res(nxf,nyf,nzf,3,nresmax))
            do i=1,nresults
                do j=1,3
                    res(ix,iy,iz,j,i) = results(j,i)
                end do
            end do

            lfirst = .false.
        end do  ! x
    end do      ! y
    if ( lwrite ) print *,'deallocating field'
    deallocate(field)
    if ( lwrite ) print *,'deallocating covariate'
    deallocate(covariate)
    if ( lwrite ) print *,'deallocating series'
    deallocate(series)
    
    ! write output
    
    if ( nresults /= 38 ) then
        write(0,*) 'attributefield: internal error: expecting nresults = 37, not ',nresults
        call exit(-1)
    end if
!
!   metadata
!  
    vars(1) = 'mu'
    units(1) = orgunits
    lvars(1) = 'position parameter'
    vars(2) = 'sigma'
    units(2) = orgunits
    lvars(2) = 'scale parameter'
    vars(3) = 'xi'
    units(3) = '1'
    lvars(3) = 'shape parameter'
    vars(4) = 'alpha'
    units(4) = trim(orgunits)//'/'//trim(units1)
    lvars(4) = 'trend parameter'
    vars(5) = 'beta'
    units(5) = trim(orgunits)//'/'//trim(units1)
    lvars(5) = 'second trend parameter'
    write(vars(6),'(a,i4.4)') 'tx',yr1a
    units(6) = 'yr'
    write(lvars(6),'(a,i4.4)') 'return time in climate of ',yr1a
    write(vars(7),'(a,i4.4)') 'tx',yr2a
    units(7) = 'yr'
    write(lvars(7),'(a,i4.4)') 'return time in climate of ',yr2a
    vars(8) = 'ratio'
    units(8) = '1'
    write(lvars(5),'(a,i4.4)') 'ratio of return times'
    k = 8
    do i=1,10
        l = 2 + (i-1)/3
        write(format,'(a,i1,a)') '(a,i',l,',a,i4.4)'
        k = k + 1
        write(vars(k),format) 't',iyrs(i),'_',yr1a
        units(k) = orgunits
        write(lvars(k),format) ' ',iyrs(i),'-yr-return value in climate of ',yr1a
    end do
    do i=1,10
        l = 2 + (i-1)/3
        write(format,'(a,i1,a)') '(a,i',l,',a,i4.4)'
        k = k + 1
        write(vars(k),format) 't',iyrs(i),'_',yr2a
        units(k) = orgunits
        write(lvars(k),format) ' ',iyrs(i),'-yr-return value in climate of ',yr2a
    end do
!   error bars
    do i=1,k
        vars(k+i) = 'lo'//vars(i)
        units(k+i) = units(i)
        lvars(k+i) = 'lower bound of 95% CI of '//lvars(i)
        vars(2*k+i) = 'hi'//vars(i)
        units(2*k+i) = units(i)
        lvars(2*k+i) = 'upper bound of 95% CI of '//lvars(i)
    end do
    nvars = 3*k
    if ( nz.le.1 ) then
        ivars(1,1:nvars) = 0
    else
        ivars(1,1:nvars) = 1
    end if

    call getarg(iargc(),outfile)
    title = 'Fit of '//trim(file)//' with '//trim(distribution)//' dependent on'// &
    &   trim(covariatefile)
    call writenc(outfile,ncid,ntvarid,itimeaxis,ntmax,nx,xx,ny,yy &
     &       ,nz,zz,1,1,2000,1,3e33,title,nvars &
     &       ,vars,ivars,lvars,units,0,0)
!
!   data
!
    do i=1,3 ! central value; lower, upper bounds 95% CI
        do j=1,nvars/3
            call writencslice(ncid,ntvarid,itimeaxis,ntmax,ivars,res(1,1,1,i,j) &
     &        ,nx,ny,nz,nx,ny,nz,1,1)
        end do
    end do
    i = nf_close(ncid)
    if ( lwrite ) print *,'deallocating res'
    deallocate(res)

end program attributefield

