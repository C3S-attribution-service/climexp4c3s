program spectrum

!   compute a primitive spectrum of a climate series

    implicit none
    include 'param.inc'
    include 'getopts.inc'
    integer :: i,j,j1,j2,n,nout,nperyear,mens1,mens &
        ,iens,nn,yr,mo,yrstart,yrstop,imc,nmc,nprob,nimax
    real :: epx(4*ndata),epy(4*ndata),ofac,fac,xmax,ymax,alpha,s2,sd,probmax,mean,mean1
    real,allocatable :: epx1(:),epy1(:),epyall(:,:),prob(:),data(:,:,:)
    integer,allocatable :: imax(:)
    character line*255,file*1024,var*80,units*60,lvar*120,svar*120,history*50000,metadata(2,100)*2000

!   init

    lwrite = .false. 
    if ( command_argument_count() == 0 ) then
        print *,'usage: spectrum datafile [month n] [sum m] '// &
            '[detrend] [anom] [ensanom] [diff [n]] [begin yyyy] [end yyyy] [xave n]'
        call exit(-1)
    end if

    n = 0
    call get_command_argument(1,file)
    allocate(data(npermax,yrbeg:yrend,0:nensmax))
    call readensseriesmeta(file,data,npermax,yrbeg,yrend,nensmax &
        ,nperyear,mens1,mens,var,units,lvar,svar,history,metadata,.false.,lwrite)
    call getopts(2,command_argument_count(),nperyear,yrbeg,yrend,.true.,mens1,mens)
!   I am abusing yr1a,yr2a here...
    if ( yr1a == yr1 ) yr1a = yrbeg
    if ( yr2a == yr2 ) yr2a = yrend + 1
    yrstart = yr2
    yrstop  = yr1
    if ( mens > 0 ) then
        write(0,*) 'The spectra are computed for ensemble members ',nens1,' to ',nens2
        if ( nens2-nens1 > 0 ) write(0,*) ' and then averaged'
        write(0,*) '<br>'
        print '(a,i4,a,i4)','# taking ensemble members ',nens1,' to ',nens2
    end if
    if ( m1 == 0 ) then
        j1 = 1
        j2 = nperyear
    else
        j1 = m1
        j2 = m1
        if ( nperyear /= 12 ) then
            call month2period(j1,nperyear,1)
            call month2period(j2,nperyear,0)
            if ( j2 < j1 ) j2 = j2 + nperyear
        end if
    end if
    if ( lwrite ) then
        if ( nperyear == 12  ) then
            print '(a,i2,a,i2)','# spectrum: using monhs ',j1,' to ',j2
        else
            print '(a,i2,a,i2)','# spectrum: using periods ',j1,' to ',j2
        end if
    end if

!   sum

    if ( lsum > 1 ) then
        if ( lwrite ) print '(a)','# Summing data'
        do iens=nens1,nens2
            call sumit(data(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend,lsum,oper)
        end do
    end if

!   logscale

    if ( logscale ) then
        if ( lwrite ) print '(a,2i3)','# Taking log of series '
        do iens=nens1,nens2
            call takelog(data(1,yrbeg,iens),nperyear,nperyear,yrbeg,yrend)
        end do
    end if

!   sqrtscale

    if ( sqrtscale ) then
        if ( lwrite ) print '(a,2i3)','# Taking sqrt of series '
        do iens=nens1,nens2
            call takesqrt(data(1,yrbeg,iens),nperyear,nperyear,yrbeg,yrend)
        end do
    end if

!   squarescale

    if ( squarescale ) then
        if ( lwrite ) print '(a,2i3)','# Taking sqrt of series '
        do iens=nens1,nens2
            call takesquare(data(1,yrbeg,iens),nperyear,nperyear,yrbeg,yrend)
        end do
    end if

!   detrending

    if ( ldetrend ) then
        if ( lwrite ) print *,'Detrending data'
        do iens=nens1,nens2
            call detrend(data(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend,yr1,yr2,m1,m1,1)
        end do
    end if
    if ( ndiff /= 0 ) then
        if ( lwrite ) print *,'Taking anomalies - data'
        do iens=nens1,nens2
            call diffit(data(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend,ndiff)
        end do
    end if

!   anomalies w.r.t. seasonal cycle

    if ( anom ) then
        do iens=nens1,nens2
            call anomal(data(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend,yrbeg,yrend)
        end do
    end if

!   anomalies wrt ensemble mean

    if ( lensanom .and. min(nens2,mens) > max(nens1,mens1) ) then
        if ( lwrite ) print *,'Taking anomalies wrt ensemble mean'
        call anomalensemble(data,npermax,nperyear,yrbeg,yrend, &
            yr1,yr2,max(nens1,mens1),min(nens2,mens))
    end if

    call getspectrum(data,npermax,nperyear,yrbeg,yrend,yr1,yr2,mens1 &
        ,mens,nens1,nens2,j1,j2,avex,epx,epy,ndata,nout,ofac,mean &
        ,yrstart,yrstop,lwrite,yr1a,yr2a)

!       get error bars by comparing it to an AR(1) process

    i = (yrstop-yrstart+1)*(j2-j1+1)
    nmc = 1000000000/i**2
    nmc = min(5000,nmc)
    nmc = max(200,nmc)
    !!!print '(a,3i6)','# i,nmc = ',nout,nmc
    allocate(epx1(4*ndata))
    allocate(epy1(4*ndata))
    allocate(epyall(nmc,nout))
    allocate(prob(nout))
    allocate(imax(nout))
    call getred1(alpha,s2,j1,j2,nperyear,data,npermax,yrbeg,yrend,nensmax)
    sd = sqrt(s2)
    if ( lwrite ) print *,'alpha = ',alpha
    do imc=1,nmc
        call keepalive1('Bootstrapping',imc,nmc)
!!!         if ( mod(imc,10).eq.0 ) print '(a,i6)','# ',imc
        call make1mcseries(alpha,sd,data,npermax,nperyear,yrbeg &
            ,yrend,nens1,nens2,j1,j2,yrstart,yrstop,lwrite)
        call getspectrum(data,npermax,nperyear,yrbeg,yrend,yr1,yr2 &
            ,mens1,mens,nens1,nens2,j1,j2,avex,epx1,epy1,ndata,nout &
            ,ofac,mean1,yrstart,yrstop,lwrite,yr1a,yr2a)
        do i=nint(ofac),nout
            if ( epx(i) /= epx1(i) ) then
                write(0,*) 'warning: x-coordinates not identical ',i,epx(i),epx1(i)
            end if
            epyall(imc,i) = epy1(i)
        end do
    end do
    nprob = 0
    probmax = 1
    nimax = 0
    do i=nint(ofac),nout
        call nrsort(nmc,epyall(1,i))
        call getsign1(prob(i),epy(i),epyall(1,i),nmc,lwrite)
        prob(i) = 1-prob(i)
        if ( prob(i) < 0.05 ) nprob = nprob + 1
        if ( prob(i) < probmax ) then
            nimax = 1
            imax(nimax) = i
            probmax = prob(i)
        else if ( prob(i) == probmax ) then
            nimax = nimax + 1
            imax(nimax) = i
        end if
    end do

!   output; do not yet know how to treat dependent data.

    call printmetadata(6,trim(file),' ','spectrum of ',history,metadata)
    print '(a)','# frequency power AR(1)  prob'
    print '(3a)','# power in [',trim(units),'^2]'
    if ( nperyear >= 53 ) then
        print '(a)','# frequency in [mo^-1]'
        fac = 12
    else
        print '(a)','# frequency in [yr^-1]'
        fac = 1
    end if
    do i=nint(ofac),nout
        if ( epx(i) > 0 ) then
            write(*,'(4f12.6)') epx(i)/fac,epy(i)*fac, &
            epyall(nint(0.95*nmc),i)*fac,prob(i)
        else
            write(0,*) 'spectrum: something is wrong',i,epx(i)
            call exit(-1)
        end if
    end do
    if ( avex > 1 ) then
        write(0,'(a)') 'The spectrum has been obtained by averaging'
        write(0,'(i3,a)') avex,' bins of the periodogram<br>'
    end if
    write(0,'(a,i6,a,f6.3,a)') 'The second line denotes the 95% highest spectrum of ',nmc &
        ,' AR(1) processes with the same autocorrelation (',alpha,').'
    if ( nprob <= 1 ) then
        write(0,'(a,i6)') 'There is ',nprob
    else
        write(0,'(a,i6)') 'There are ',nprob
    end if
    write(0,'(a,f10.1,a)') ' bins with p&lt;0.05, by chance one would expect ', &
        0.05*(nout-nint(ofac)+1),'.'
    if ( nprob > 1.2*0.05*(nout-nint(ofac)+1) ) then
        if ( nimax == 1 ) then
            write(0,'(a,f8.2,a,f6.3,a)') 'The most significant peak is at ',1/epx(imax(1)) &
                ,' years (p=',probmax,').'
        else
            write(0,'(a)') 'The most significant peaks are at '
            do i=nimax,2,-1
                write(0,'(f8.2,a)') 1/epx(imax(i)),', '
            end do
            write(0,'(a,f8.2,a,f6.3,a)') ' and ',1/epx(imax(1)), ' years (p=',probmax,').'
        end if
    end if
    if ( yr1a > yrbeg .or. yr2a < yrend ) then
        write(0,'(2a,i4,a,i4,a,f16.3)') &
            'The mean of the spectrum exp(&lang;log(1/f)&rang;) ', &
            ' band-pass filtered between ',yr1a-yrbeg,' and ', &
        yr2a-yrbeg,' yrs is at period ',mean, ' yr '
    else
        write(0,'(2a,f16.3)') &
            'The mean of the spectrum exp(&lang;log(1/f)&rang;) ', &
            'is at period ',mean, ' yr '
    end if
    call savestartstop(yrstart,yrstop)
    goto 999
900 print *,'error: cannot read integer from ',file
    call exit(-1)
999 continue
end program spectrum
