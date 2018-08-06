program lomb

!   compute the Lomb periodogram of a climate series

    implicit none
    include 'param.inc'
    include 'getopts.inc'
    integer :: i,j,j1,j2,n,n1,n2,nout,jmax,nperyear,mens1,mens &
        ,iens,nepx(4*ndata),nn,yr,mo,yrstart,yrstop
    real :: x(ndata),y(ndata),px(4*ndata),py(4*ndata), &
        epx(4*ndata),epy(4*ndata),prob,hifac,ofac,fac,sx,sy, &
        xmax,ymax,timescale
    real,allocatable :: data(:,:,:)
    character :: line*255,var*40,units*20

!   init

    lwrite = .false. 
    if ( command_argument_count() == 0 ) then
        print *,'usage: lomb datafile [month n] [sum m] [detrend]'// &
            ' [anom] [ensanom] [diff [n]] [begin yyyy] [end yyyy] [xave n]'
        call exit(-1)
    endif

    n = 0
    call get_command_argument(1,line)
    allocate(data(npermax,yrbeg:yrend,0:nensmax))
    call readensseries(line,data,npermax,yrbeg,yrend,nensmax &
        ,nperyear,mens1,mens,var,units, .false. ,lwrite)
    call getopts(2,command_argument_count(),nperyear,yrbeg,yrend, .true. ,mens1,mens)
    yrstart = yr2
    yrstop  = yr1
    if ( mens > 0 ) then
        write(0,*) 'The spectra are computed for ensemble members ',nens1,' to ',nens2
        if ( nens2-nens1 > 0 ) write(0,*) ' and then averaged'
        write(0,*) '<br>'
        print '(a,i4,a,i4)','# taking ensemble members ',nens1,' to ',nens2
    endif
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
        endif
    endif
    if ( lwrite ) then
        if ( nperyear == 12  ) then
            print '(a,i2,a,i2)','# lomb: using monhs ',j1,' to ',j2
        else
            print '(a,i2,a,i2)','# lomb: using periods ',j1,' to ',j2
        endif
    endif

!   sum

    if ( lsum > 1 ) then
        if ( lwrite ) print '(a)','# Summing data'
        do iens=nens1,nens2
            call sumit(data(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend,lsum,oper)
        enddo
    endif

!   logscale

    if ( logscale ) then
        if ( lwrite ) print '(a,2i3)','# Taking log of series '
        do iens=nens1,nens2
            call takelog(data(1,yrbeg,iens),nperyear,nperyear,yrbeg,yrend)
        enddo
    endif

!   sqrtscale

    if ( sqrtscale ) then
        if ( lwrite ) print '(a,2i3)','# Taking sqrt of series '
        do iens=nens1,nens2
            call takesqrt(data(1,yrbeg,iens),nperyear,nperyear,yrbeg,yrend)
        enddo
    endif

!   squarescale

    if ( squarescale ) then
        if ( lwrite ) print '(a,2i3)','# Taking sqrt of series '
        do iens=nens1,nens2
            call takesquare(data(1,yrbeg,iens),nperyear,nperyear,yrbeg,yrend)
        enddo
    endif

!   detrending

    if ( ldetrend ) then
        if ( lwrite ) print *,'Detrending data'
        do iens=nens1,nens2
            call detrend(data(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend,yr1,yr2,m1,m1,1)
        enddo
    endif
    if ( ndiff /= 0 ) then
        if ( lwrite ) print *,'Taking anomalies - data'
        do iens=nens1,nens2
            call diffit(data(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend,ndiff)
        enddo
    endif

!       anomalies w.r.t. seasonal cycle

    if ( anom ) then
        do iens=nens1,nens2
            call anomal(data(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend,yrbeg,yrend)
        enddo
    endif

!   anomalies wrt ensemble mean

    if ( lensanom .and. min(nens2,mens) > max(nens1,mens1) ) then
        if ( lwrite ) print *,'Taking anomalies wrt ensemble mean'
        call anomalensemble(data,npermax,nperyear,yrbeg,yrend, &
            yr1,yr2,max(nens1,mens1),min(nens2,mens))
    endif

!   loop over ensemble members

    do iens=nens1,nens2
    
!       fill arrays
    
        n = 0
        n1 = 0
        do yr=yr1-1,yr2
            do mo=j1,j2
                j = mo
                call normon(j,yr,i,nperyear)
                if ( i < yr1 .or. i > yr2 ) cycle
                if ( abs(data(j,i,iens)) < 1e33 ) then
                    n = n + 1
                    yrstart = min(yrstart,i)
                    yrstop  = max(yrstop,i)
                    x(n) = i + (j-0.5)/nperyear
                    if ( nperyear > 12 ) then
                        x(n) = x(n)*365.24
                    endif
                    y(n) = data(j,i,iens)
                    if ( n1 == 0 ) then
                        n1 = nperyear*i+j
                    endif
                    n2 = nperyear*i+j
                endif
            enddo
        enddo
    
!       compute/estimate other parameters
    
        if ( m1 == 0 ) then
            hifac = real(n)/real(n2-n1+1)
        else
            hifac = real(n)/real((n2-n1)/nperyear+1)
        endif
        if ( hifac == 1 ) then
            ofac = 1
        else
            ofac = 4        ! see how it works
        endif
    
!       call period (Numerical recipes p 572)
!       take care of dependent data!
    
        if ( lwrite ) then
            print *,'call Numerical Recipes period'
            print *,'x = ',(x(i),i=1,min(n,5))
            print *,'y = ',(y(i),i=1,min(n,5))
            print *,'ofac = ',ofac
            print *,'hifac = ',hifac
            print *,'4*ndata = ',4*ndata
        endif
        call period(x,y,n,ofac,hifac,px,py,4*ndata,nout,jmax,prob)
    
!       average
    
        if ( avex > 1 ) then
            do i=nint(ofac),nout-avex,avex
                sx = px(i)
                sy = py(i)
                do j=1,avex-1
                    sx = sx + px(i+j)
                    sy = sy + py(i+j)
                enddo
                n = 1+i/avex
                px(n) = sx/avex
                py(n) = sy/avex
            enddo
            ofac = 1+nint(ofac)/avex
            nout = n
        endif
    
!       collect ensemble informnation
    
        if ( iens == nens1 ) then
            nn = nout
            do i=1,nn
                nepx(i) = 1
                epx(i) = px(i)
                epy(i) = py(i)
                if ( lwrite ) print *,i,epx(i),epy(i),nepx(i)
            enddo
        else
            do i=1,nout
                if ( epx(i) == px(i) ) then
                    nepx(i) = nepx(i) + 1
                    epy(i) = epy(i) + py(i)
                else        ! unequal array sizes - choose nearest
                    do j=1,nn-1
                        if ( px(j) < (epx(j)+epx(j+1))/2 ) then
                            goto 100
                        endif
                    enddo
                    100 continue
                    epx(j) = (nepx(j)*epx(j) + px(i))/(nepx(j) + 1)
                    nepx(j) = nepx(j) + 1
                    epy(j) = epy(j) + py(i)
                endif
            enddo
        endif
    enddo
    do i=1,nn
        epy(i) = epy(i)/nepx(i)
    enddo

!   output; do not yet know how to treat dependent data.

    print '(a)','# frequency power'
    do i=nint(ofac),nout
        if ( epx(i) > 0 ) then
            write(*,'(2f12.6)') epx(i),epy(i)
        else
            write(0,*) 'lomb: something is wrong',i,epx(i)
            call abort
        endif
    enddo

!   the 95% significance line - correct for dependent data

    write(0,'(a)') '<font color="#ff0000">Changed normalisation '// &
        ' to equal-area with log(T) x-axis</font></br>'
    fac = 1
    if ( m1 == 0 ) then
        if ( lsum > 1 ) fac = lsum
    else
        if ( lsum > nperyear ) fac = lsum/real(nperyear)
    endif
    call maxquad(xmax,ymax,epx(nint(ofac)),epy(nint(ofac)), &
    nout-nint(ofac)+1)
    if ( fac == 1 .and. avex == 1 .and. mens == 0 ) then
        write(0,'(a)') 'The horizontal line denotes the 95% significance limit.'
        write(0,'(a,f8.2,a,f6.2,a)') 'The peak at ',1/px(jmax) &
            ,'yr has a significance of ',100*(1-prob),'%.'
        write(0,'(a)') 'These numbers refer to the chance of '// &
            'any peak not exceeding this value in a white-noise spectrum.<br>'
        write(*,'(a)')
        do i=1,nout
            sy = -log(0.05/(2*nout/ofac))*fac
            if ( epx(i)*sy < xmax*ymax ) then
                write(*,'(2f12.6)') epx(i),sy
            endif
        enddo
    endif
    if ( avex > 1 ) then
        write(0,'(a)') 'The spectrum has been obtained by averaging'
        write(0,'(i3,a)') avex,' bins of the periodogram<br>'
    endif
    write(0,'(a,2f12.6,a)') 'The highest peak is at ',1/xmax,xmax*ymax,' (in a log(T) plot)<br>'

    call savestartstop(yrstart,yrstop)
    goto 999
900 print *,'error: cannot read integer from ',line
    call exit(-1)
999 continue
end program lomb
