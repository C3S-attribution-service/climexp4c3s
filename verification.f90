    program verification
!  #[ comment:

!       verify forecasts against observations.  Based on correlations.F
!       this is an empty shell, which does all the time series
!       manipulation station parameters from www.ncdc.noaa.gov or
!       other(seasons etc), and writes a table for R to do the actual work.

!       (c) Geert Jan van Oldenborgh, 7-nov-2005, KNMI

!       This file may be copied and modified  freely as long as the
!       above copyright notice is kept intact

!  #] comment:
!  #[ declarations:

    implicit none
    include 'param.inc'
    include 'getopts.inc'
    integer :: nperyear,mens1,mens,mo,yr,j,iu,iens,nperyear2,yrmin &
        ,yrmax,n,nmodel1,nmodel2,nmodel,multimodel(0:nmodelmax) &
        ,yrstart,yrstop
    real,allocatable :: obs(:,:),fcst(:,:,:)
    real :: s1,s2
    character line*255,varobs*40,unitobs*20,varfcst*4,unitfcst*20
    integer :: iargc

!  #] declarations:
!  #[ check arguments:

!       check arguments

    lwrite = .false. 
    lstandardunits = .true. 
    if ( iargc() < 4 ) then
        print *,'usage: verification forecastfile '// &
        'file obsfile dump outtable '// &
        '[month m[:n] [sum|ave|max|min|sel m] '// &
        '[begin yr] [end yr] [detrend]'// &
        '[debias none|mean|var|all] '
        call exit(-1)
    endif
!  #] check arguments:
!  #[ read forecast:

    call getarg(1,line)
    allocate(obs(npermax,yrbeg:yrend),fcst(npermax,yrbeg:yrend,0:nensmax))
    call readensseries(line,fcst,npermax,yrbeg,yrend,nensmax &
        ,nperyear2,mens1,mens,varfcst,unitfcst,lstandardunits,lwrite)
    print '(a,i2,a,i2)','# located ensemble members ',mens1,' to ',mens

!  #] read forecast:
!  #[ options:

    call getopts(2,iargc(),nperyear2,yrbeg,yrend, .true. ,mens1,mens)
    if ( .not. dump ) then
        write(0,*) 'verification: error: cannot find ''dump'' in argument list'
        write(*,*) 'verification: error: cannot find ''dump'' in argument list'
        call exit(-1)
    endif
    if ( lag1 /= 0 .or. lag2 /= 0 ) then
        write(0,*) 'verification: error: cannot handle lags ',lag1,lag2
        write(*,*) 'verification: error: cannot handle lags ',lag1,lag2
        call exit(-1)
    endif

!       read observations file (indicated by "file obsfile" in the args)

    if ( indxuse < 10 .or. .not. lincl(10) ) then
        write(0,*) 'verification: error: no obsfile specified with ''file obsfile'''
        write(0,*) 'verification: error: no obsfile specified with ''file obsfile'''
        call exit(-1)
    endif
    if ( indexfiles(10) == 'perfectmodel' ) then
        print '(a)','# computing perfect model statistic'
        nmodel1 = nens1
        nmodel2 = nens2
        nperyear = nperyear2
    else
        nmodel1 = -1
        nmodel2 = -1
        call readseries(indexfiles(10),obs,npermax,yrbeg,yrend &
        ,nperyear,varobs,unitobs,lstandardunits,lwrite)
        if ( nperyear /= nperyear2 ) then
            write(0,*) 'verification: error: unequal time intervals',nperyear,nperyear2
            write(*,*) 'verification: error: unequal time intervals',nperyear,nperyear2
            call exit(-1)
        endif
    endif

!  #] options:
!  #[ manipulate time series:

!  for the time being no multimodel information
    multimodel(0) = nens1
    multimodel(1) = nens2+1
    nmodel = 1
    call manipulatetimeseries(fcst,obs,npermax,yrbeg,yrend &
        ,nperyear,nmodel1,varobs,0,0,multimodel,nmodel)

!  #] manipulate time series:
!  #[ make verification table:

    n = 0
    call printtable(fcst,obs,npermax,yrbeg,yrend,nmodel1,nmodel2 &
        ,nens1,nens2,yr1,yr2,m1,m2,nperyear,n,yrstart,yrstop)
    if ( n == 0 ) then
        write(*,*) 'verification: did not find any valid data'
        close(10,status='delete')
    else
        close(10)
    endif
!  #] make verification table:
end program verification