program attribute
!
!   do an empirical attribution study by fitting a time series to a
!   GEV or GPD with the position parameter linearly dependent on a covariate
!   and studying the difference in return time in the current climate and a
!   previous climate.
!
    implicit none
    include 'param.inc'
    include 'getopts.inc'
    integer nperyear,nperyear1,mens1,mens
    integer i,yr,mo,n,j1,j2
    real,allocatable :: series(:,:),covariate(:,:)
    character seriesfile*1024,covariatefile*1024,distribution*6,assume*5,string*80
    character var*40,units*80,var1*40,units1*80
    integer iargc

    if ( iargc().lt.8 ) then
        write(0,*) 'usage: attribute series covariate_series ', &
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
!   getopts comes way too late for these options...
    do i=6,iargc()
        call getarg(i,string)
        if ( string.eq.'debug' .or. string.eq.'lwrite' ) then
            lwrite = .true.
        end if
        if ( string(1:9).eq.'standardu' ) then
            lstandardunits = .true.
        end if
    end do

    call getarg(3,distribution)
    call tolower(distribution)
    if ( distribution.ne.'gev' .and. distribution.ne.'gumbel' .and. &
    &    distribution.ne.'gpd' .and. distribution.ne.'gauss' ) then
        write(0,*) 'attribute: error: only GEV, GPD or Gauss supported, not ',distribution
        call abort
    end if

    call getarg(4,string)
    if ( string(1:4).ne.'assu' ) then
        write(0,*) 'attribute: error: expecting "assume", not ',trim(string)
        call abort
    end if
    call getarg(5,assume)
    call tolower(assume)
    if ( assume.ne.'shift' .and. assume.ne.'scale' .and. assume.ne.'both' ) then
        write(0,*) 'attribute: error: only shift, scale or both supported, not ',assume
        call abort
    end if

    call getarg(1,seriesfile)
    allocate(series(npermax,yrbeg:yrend))
    call readseries(seriesfile,series,npermax,yrbeg,yrend &
    & ,nperyear,var,units,lstandardunits,lwrite)
    
    call getarg(2,covariatefile)
    allocate(covariate(npermax,yrbeg:yrend))
    call readseries(covariatefile,covariate,npermax,yrbeg,yrend &
    & ,nperyear1,var1,units1,lstandardunits,lwrite)
    
    mens1 = 0
    mens = 0
    call getopts(6,iargc(),nperyear,yrbeg,yrend,.true.,mens1,mens)
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
!
!   process data
!
    if ( lchangesign ) then
        call changesign(series,npermax,nperyear,yrbeg,yrend)
    endif
    if ( ldetrend ) then
        call detrend(series,npermax,nperyear,yrbeg,yrend,yr1,yr2,m1,m2,lsel)
    end if
    if ( anom ) then
        call anomal(series,npermax,nperyear,yrbeg,yrend,yr1,yr2)
    end if
    if ( lsum.gt.1 ) then
        call sumit(series,npermax,nperyear,yrbeg,yrend,lsum,oper)
    endif
    if ( logscale ) then
        print '(a)','# taking logarithm'
        call takelog(series,npermax,nperyear,yrbeg,yrend)
    endif
    if ( sqrtscale ) then
        print '(a)','# taking sqrt'
        call takesqrt(series,npermax,nperyear,yrbeg,yrend)
    endif
    if ( squarescale ) then
        print '(a)','# taking square'
        call takesquare(series,npermax,nperyear,yrbeg,yrend)
    endif

    call attribute_dist(series,nperyear,covariate,nperyear1,npermax,yrbeg,yrend,assume,distribution)

end program attribute

