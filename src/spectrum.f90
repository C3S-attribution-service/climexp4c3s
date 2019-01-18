program spectrum

!   compute a primitive spectrum of a climate series

    implicit none
    include 'param.inc'
    include 'getopts.inc'
    integer :: i,j,j1,j2,n,nout,nperyear,mens1,mens &
        ,iens,nn,yr,mo,yrstart,yrstop,imc,nmc,nprob,nimax
    real :: ofac,fac,xmax,ymax,alpha,s2,sd,probmax,mean,mean1
    real,allocatable :: epx(:),epy(:),epx1(:),epy1(:),epyall(:,:),prob(:),data(:,:,:)
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

    allocate(epx(4*ndata))
    allocate(epy(4*ndata))
    call getspectrum(data,npermax,nperyear,yrbeg,yrend,yr1,yr2,mens1 &
        ,mens,nens1,nens2,j1,j2,avex,epx,epy,ndata,nout,ofac,mean &
        ,yrstart,yrstop,lwrite,yr1a,yr2a)

!   get error bars by comparing it to an AR(1) process

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

subroutine getspectrum(data,npermax,nperyear,yrbeg,yrend,yr1,yr2 &
    ,mens1,mens,nens1,nens2,j1,j2,avex,epx,epy,ndata,nout &
    ,ofac,mean,yrstart,yrstop,lwrite,yr1a,yr2a)

!   compute a spectrum from the time series in data

    implicit none
    integer :: npermax,nperyear,yrbeg,yrend,yr1,yr2,mens1,mens,nens1 &
        ,nens2,j1,j2,avex,ndata,nout,yrstart,yrstop,yr1a,yr2a
    real :: data(npermax,yrbeg:yrend,0:nens2),mean
    real :: epx(4*ndata),epy(4*ndata),ofac,prob
    logical :: lwrite
    integer :: iens,n,n1,n2,nn,yr,mo,i,j,jmax
    integer,allocatable :: nepx(:)
    real :: hifac,sx,sy,somx,somy
    real,allocatable :: x(:),y(:),px(:),py(:)
    integer,external :: leap

    allocate(x(ndata))
    allocate(y(ndata))
    allocate(px(4*ndata))
    allocate(py(4*ndata))
    allocate(nepx(4*ndata))

!   loop over ensemble members

    somx = 0
    somy = 0
    do iens=nens1,nens2
    
!       fill arrays
    
        n = 0
        n1 = 0
        do yr=yr1-1,yr2
            do mo=j1,j2
                j = mo
                call normon(j,yr,i,nperyear)
                if ( i < yr1 .or. yr2 > yr2 ) cycle
                if ( abs(data(j,i,iens)) < 1e33 ) then
                    n = n + 1
                    yrstart = min(yrstart,i)
                    yrstop  = max(yrstop,i)
                    if ( nperyear /= 366 ) then
                        x(n) = i + (j-0.5)/nperyear
                    else
                        if ( leap(i) == 1 ) then
                            if ( j < 60 ) then
                                x(n) =  i + (j-0.5)/365
                            else
                                x(n) =  i + (j-1.5)/365
                            end if
                        else
                            x(n) = i + (j-0.5)/366
                        end if
                    end if
                    y(n) = data(j,i,iens)
                    if ( n1 == 0 ) then
                        n1 = nperyear*i+j
                    end if
                    n2 = nperyear*i+j
                end if
            end do
        end do
        if ( n == 0 ) then
            write(0,*) 'Spectrum: error: no valid data found'
            call exit(-1)
        end if
    
!       compute/estimate other parameters
    
        if ( j1 /= j2 ) then
            hifac = real(n)/real(n2-n1+1)
        else
            hifac = real(n)/real((n2-n1)/nperyear+1)
        end if
        if ( hifac == 1 ) then
            ofac = 1
        else
            ofac = 4        ! see how it works
        end if
    
!       call period (Numerical recipes p 572)
!       take care of dependent data!

        if ( lwrite ) then
            print *,'call Numerical Recipes period'
            print *,'x = ',(x(i),i=1,min(n,5))
            print *,'y = ',(y(i),i=1,min(n,5))
            print *,'ofac = ',ofac
            print *,'hifac = ',hifac
            print *,'4*ndata = ',4*ndata
        end if
        call period(x,y,n,ofac,hifac,px,py,4*ndata,nout,jmax,prob)
        if ( lwrite ) then
            print *,'back from period'
            print *,'px = ',(px(i),i=1,min(nout,5))
            print *,'py = ',(py(i),i=1,min(nout,5))
        end if
    
!       compute mean
    
        do i=2,nout-1
            if ( 1/px(i) > yr1a-yrbeg .and. 1/px(i) < yr2a-yrbeg ) &
            then
                somx = somx + py(i)*(px(i+1)-px(i-1))/2
                somy = somy + py(i)*log(px(i))*(px(i+1)-px(i-1))/2
            else if ( lwrite ) then
                write(*,*) 'disregarding point at T=',1/px(i)
            end if
        end do
    
!       average
    
        if ( avex > 1 ) then
            do i=nint(ofac),nout-avex,avex
                sx = px(i)
                sy = py(i)
                do j=1,avex-1
                    sx = sx + px(i+j)
                    sy = sy + py(i+j)
                end do
                n = 1+i/avex
                px(n) = sx/avex
                py(n) = sy/avex
            end do
            ofac = 1+nint(ofac)/avex
            nout = n
        end if
    
!       collect ensemble informnation
    
        if ( iens == nens1 ) then
            nn = nout
            do i=1,nn
                nepx(i) = 1
                epx(i) = px(i)
                epy(i) = py(i)
                if ( lwrite ) print *,i,epx(i),epy(i),nepx(i)
            end do
        else
            do i=1,nout
                if ( epx(i) == px(i) ) then
                    nepx(i) = nepx(i) + 1
                    epy(i) = epy(i) + py(i)
                    if ( lwrite ) print *,i,epx(i),epy(i),nepx(i)
                else        ! unequal array sizes - choose nearest
                    do j=1,nn-1
                        if ( px(i) < (epx(j)+epx(j+1))/2 ) then
                            goto 100
                        end if
                    end do
                    100 continue
                    epx(j) = (nepx(j)*epx(j) + px(i))/(nepx(j) + 1)
                    nepx(j) = nepx(j) + 1
                    epy(j) = epy(j) + py(i)
                    if ( lwrite ) print *,j,epx(j),epy(j),nepx(j)
                end if
            end do
        end if
    end do
    do i=1,nn
        epy(i) = epy(i)/nepx(i)
        if ( lwrite ) print *,i,epx(i),epy(i),nepx(i)
    end do
    mean = exp(-somy/somx)
    if ( lwrite ) write(*,'(a,3f16.4)') '# getspectrum: mean = ',mean,somy,somx
end subroutine getspectrum
