program attribute
!
!   do an empirical attribution study by fitting a time series to a
!   Gaussian, Gumbell, GEV or GPD with the position parameter and/or scale parameter
!   linearly dependent on a covariate and studying the difference in return time
!   in the current climate and a previous climate.
!
    implicit none
    include 'param.inc'
    include 'getopts.inc'
    integer nresmax
    parameter(nresmax=100)
    integer nperyear,nperyear1,mens1,mens,mmens1,mmens,iens,nresults
    integer i,yr,mo,n,off
    real results(3,nresmax)
    real,allocatable :: series(:,:,:),covariate(:,:,:)
    character seriesfile*1024,covariatefile*1024,distribution*6,assume*5,string*80
    character var*40,units*80,var1*40,units1*80,seriesids(0:nensmax)*30
    logical lprint
    integer iargc
    real scalingpower
    common /c_scalingpower/ scalingpower

    if ( iargc().lt.8 ) then
        write(0,*) 'usage: attribute series covariate_series|none ', &
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
    call attribute_init(seriesfile,distribution,assume,off,nperyear,yrbeg,yrend,nensmax,lwrite)

    allocate(series(npermax,yrbeg:yrend,0:nensmax))
    if ( seriesfile == 'file' ) then
        ! set of stations
        call readsetseries(series,seriesids,npermax,yrbeg,yrend &
        & ,nensmax,nperyear,mens1,mens,var,units,lstandardunits,lwrite)
    else if ( seriesfile == 'gridpoints' ) then
        ! netcdf file with gridpoints
        call readgridpoints(series,seriesids,npermax,yrbeg,yrend &
        & ,nensmax,nperyear,mens1,mens,var,units,lstandardunits,lwrite)
    else
        ! simple data
        call readensseries(seriesfile,series,npermax,yrbeg,yrend &
        & ,nensmax,nperyear,mens1,mens,var,units,lstandardunits,lwrite)
        if ( mens.gt.mens1 ) then
            do i=mens1,mens
                write(seriesids(i),'(i3.3)') i
            enddo
        else
            seriesids(mens) = ' '
        end if
    end if
    
    call getarg(2+off,covariatefile)
    allocate(covariate(npermax,yrbeg:yrend,0:nensmax))
    if ( covariatefile == 'none' ) then
        assume = 'none'
        covariate = 0
        nperyear1 = 1
    else if ( index(covariatefile,'%%') == 0 .and. &
    &    index(covariatefile,'++') == 0 ) then
        call readseries(covariatefile,covariate,npermax,yrbeg,yrend &
        &   ,nperyear1,var1,units1,lstandardunits,lwrite)
        do iens=mens1,mens
            if ( iens /= 0 ) then
                covariate(:,:,iens) = covariate(:,:,0)
            end if
        end do
    else
        call readensseries(covariatefile,covariate,npermax,yrbeg,yrend &
        & ,nensmax,nperyear1,mmens1,mmens,var1,units1,lstandardunits,lwrite)
        if ( mmens1 /= mens1 .or. mmens /= mens ) then
            write(0,*) 'attribute: error: number of ensemble members should be the same ' &
            & ,'found covariate: ',mmens1,'-',mmens,', series ',mens1,'-',mens
            call exit(-1)
        end if
    end if
    
    call getopts(6+off,iargc(),nperyear,yrbeg,yrend,.true.,mens1,mens)
    if ( yr1a.lt.yr1 .or. yr1a.lt.yrbeg ) then
        write(0,*) 'attribute: error: reference year should be after start of series ',yr1,yr1a
        call exit(-1)
    end if
    if ( yr2a.gt.yr2 .or. yr2a.gt.yrend ) then
        write(0,*) 'attribute: error: current year should be before end of series ',yr2,yr2a
        call exit(-1)
    end if
!
!   process data
!
    if ( biasmul /= 1 .or. biasadd /= 0 ) then
        write(0,*) 'Applying bias correction of scale ',biasmul,' and offset ',biasadd
        print '(a,g20.4,a,g20.4,a)','# Applying bias correction of scale ',biasmul,' and offset ',biasadd
        do iens=mens1,mens
            call scaleseries(series(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend,biasmul,biasadd,0)
        end do
    end if
    if ( ldetrend ) then
        do iens=mens1,mens
            call detrend(series(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend,yr1,yr2,m1,m2,lsel)
            if ( assume /= 'none' ) then
                call detrend(covariate(1,yrbeg,iens),npermax,nperyear1,yrbeg,yrend,yr1,yr2,m1,m2,lsel)
            end if
        end do
    end if
    if ( anom ) then
        if ( assume == 'scale' ) then
            write(0,*) 'error: it makes no sense to use "scale" on anomalies'
            write(*,*) 'error: it makes no sense to use "scale" on anomalies'
            call exit(-1)
        end if
        do iens=mens1,mens
            call anomal(series(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend,yr1,yr2)
        end do
    end if
    if ( lsum.gt.1 ) then
        do iens=mens1,mens
            call sumit(series(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend,lsum,oper)
        end do
    endif
    scalingpower = 1
    if ( logscale ) then
        print '(a)','# taking logarithm'
        do iens=mens1,mens
            call takelog(series(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend)
        end do
        if ( xyear.lt.1e33 ) xyear = log(xyear)
        scalingpower = 0
    endif
    if ( sqrtscale ) then
        print '(a)','# taking sqrt'
        do iens=mens1,mens
            call takesqrt(series(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend)
        end do
        if ( xyear.lt.1e33 ) xyear = sqrt(xyear)
        scalingpower = scalingpower*0.5
    endif
    if ( squarescale ) then
        print '(a)','# taking square'
        do iens=mens1,mens
            call takesquare(series(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend)
        end do
        if ( xyear.lt.1e33 ) xyear = xyear**2
        scalingpower = scalingpower*2
    endif
    if ( cubescale ) then
        print '(a)','# taking square'
        do iens=mens1,mens
            call takecube(series(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend)
        end do
        if ( xyear.lt.1e33 ) xyear = xyear**3
        scalingpower = scalingpower*3
    endif
    if ( twothirdscale ) then
        print '(a)','# taking power two-third'
        do iens=mens1,mens
            call taketwothird(series(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend)
        end do
        if ( xyear.lt.1e33 .and. xyear.ge.0 ) xyear = xyear**(2/3.)
        scalingpower = scalingpower*2./3.
    endif
    if ( lchangesign ) then
        do iens=mens1,mens
            call changesign(series(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend)
        end do
        if ( xyear.lt.1e30 ) then
            xyear = -xyear
        end if
    endif

    lprint = .true.
    if ( lwrite ) print *,'attribute: calling attribute_dist'
    call attribute_dist(series,nperyear,covariate,nperyear1,npermax,yrbeg,yrend,&
    &   mens1,mens,assume,distribution,seriesids,results,nresmax,nresults,lprint)

end program attribute

