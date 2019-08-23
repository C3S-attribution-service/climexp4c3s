program runningmoments

!   analyse running mean, s.d., skew, ... of timeseries

    implicit none
    include 'param.inc'
    include 'getopts.inc'
    integer,parameter :: nboot=799
    integer :: i,ii,j,jj,nperyear,imom,iens,mens1,mens,n,yr,j1,j2,mboot &
        ,yrmin,yrmax,iboot,iran,yr1s,yr2s,ndatens
    real,allocatable :: data(:,:,:),xx(:),yy(:)
    real :: result(0:nboot),rmin(0:nboot),rmax(0:nboot),rdif(0:nboot) &
        ,ave,adev,sdev,sd2,skew,curt,ave1,adev1,sdev1,sd21,skew1 &
        ,curt1,signmin,signmax,signdif,noise,alpha
    logical :: first
    character :: file*1024,string*100,momentvar*4,var*80,units*80,lvar*120,svar*120,history*50000, &
        metadata(2,100)*2000
    real :: gasdev
    external gasdev
    data iran /0/

    n = command_argument_count()
    if ( n < 2 ) then
        print *,'usage: runningmoments timeseries ', &
            '1|2|3|ave|sd|skew runw n outfile [options]'
        call exit(-1)
    end if

    call get_command_argument(1,file)
    print '(2a)','# file ',file(1:index(file,' ')-1)
    allocate(data(npermax,yrbeg:yrend,0:nensmax))
    call readensseriesmeta(file,data,npermax,yrbeg,yrend,nensmax &
        ,nperyear,mens1,mens,var,units,lvar,svar,history,metadata,lstandardunits,lwrite)
    call get_command_argument(2,string)
    if ( string(1:3) == 'mea' .or. string(1:3) == 'ave' .or. string(1:1) == '1' ) then
        imom = 1
    else if ( string(1:2) == 'sd' .or. string(1:4) == 's.d.' .or. &
              string(1:3) == 'momentvar' .or. string(1:1) == '2' ) then
        imom = 2
    else if ( string(1:3) == 'ske' .or. string(1:1) == '3' ) then
        imom = 3
    else if ( string(1:4) == 'curt' .or. string(1:1) == '4' ) then
        imom = 4
    else
        write(0,*) 'runningmoments: unknown moment ',string
        call exit(-1)
    end if
    call getopts(3,n,nperyear,yrbeg,yrend, .true. ,mens1,mens)
    if ( mens > 0 ) then
        print '(a,i3,a,i3)','# Using ensemble members ',nens1,' to ',nens2
    end if
    if ( nyrwindow <= 0 ) then
        write(0,*) 'runningmoments: please specify running window length >0, not ',nyrwindow
        call exit(-1)
    end if

!   sum series

    do iens=nens1,nens2
        if ( lsum > 1 ) then
            if ( iens == nens1 ) print '(a,i4,a)','# summing over ',lsum,' periods'
            call sumit(data(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend,lsum,oper)
            decor = max(decor,real(lsum-1))
        end if

!       apply minindx, maxindx cuts

        if ( maxindx < 1e33 .or. minindx > -1e33 ) then
            do i=yr1,yr2
                do j=1,nperyear
                    if (  data(j,i,iens) > maxindx .or. data(j,i,iens) < minindx ) then
                        data(j,i,iens) = 3e33
                    end if
                end do
            end do
        end if

!       log, sqrt
    
        if ( logscale ) then
            if ( iens == nens1 ) print '(a)','# taking logarithm'
            call takelog(data(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend)
        end if
        if ( sqrtscale ) then
            if ( iens == nens1 ) print '(a)','# taking sqrt'
            call takesqrt(data(1,yrbeg,iens),npermax ,nperyear,yrbeg,yrend)
        end if
        if ( squarescale ) then
            if ( iens == nens1 ) print '(a)','# taking square'
            call takesquare(data(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend)
        end if

!       differentiate data

        if ( ndiff /= 0 ) then
            if ( iens == nens1 ) print '(a,i4)','# taking differences/averaging ',ndiff
            call diffit(data(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend,ndiff)
        end if

!       anomalies - necessary if we consider more than one month

        if ( anom ) then
            if ( iens == nens1 ) print '(a)','# taking anomalies'
            call anomal(data(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend,yr1,yr2)
        end if
    end do

!   anomalies wrt ensemble mean

    if ( lensanom .and. nens1 /= nens2 ) then
        if ( lwrite ) print '(a)','# taking anomalies wrt ensemble mean'
        call anomalensemble(data,npermax,nperyear,yrbeg,yrend,yr1,yr2,nens1,nens2)
    end if

!   overall moment

    call getj1j2(j1,j2,m1,nperyear,.true.)
    if ( minnum <= 0 ) then
        minnum = nint(0.8*nyrwindow)
        write(14,'(a,i6,a,i6,a)') '# Demanding at least ',minnum &
            ,' years in a sliding window of ',nyrwindow,' years'
    end if
    n=0
    first = .true. 
    ndatens = nperyear*(yr2-yr1+1)*(nens2-nens1+1)
    allocate(xx(ndatens))
    allocate(yy(ndatens))
    do iens=nens1,nens2
        do i=yr1-1,yr2+1
            do jj=j1,j2
                j = jj
                call normon(j,i,ii,nperyear)
                if ( ii < yr1 .or. ii > yr2 ) goto 600
                if ( data(j,ii,iens) < 1e33 ) then
                    if ( first ) then
                        first = .false. 
                        yrmin = i
                    end if
                    yrmax = i
                    n = n + 1
                    if ( n > ndatens ) goto 901
                    xx(n) = data(j,ii,iens)
                end if
            600 continue
            end do
        end do
    end do
    yr1 = max(yrmin,yr1)
    yr2 = min(yrmax,yr2)
    if ( n < minnum*(j2-j1+1)*(nens2-nens1+1) ) then
        if ( lwrite ) print *,'not enough points: ',n,minnum
        goto 900
    end if
    if ( lwrite ) print *,'analysing ',n,' points'
    call moment(xx,n,ave,adev,sdev,sd2,skew,curt)
    call bootmoment(imom,xx,yy,n,nboot,decor,result)
    if ( imom == 1 ) then
        momentvar = 'mean'
    else if ( imom == 2 ) then
        momentvar = 's.d.'
    else if ( imom == 3 ) then
        momentvar = 'skew'
    else if ( imom == 4 ) then
        momentvar = 'curt'
    else
        momentvar = '????'
    end if
    call printmetadata(14,file,' ','running '//momentvar//' of',history,metadata)
    if ( lweb ) then
        write(14,'(a)') '# <table class="realtable" '// &
            'border=0 cellpadding=0 cellspacing=0>'// &
            '<tr><th>moment</th><th>value</th>'// &
            '<th>95% CI</th>'
        write(14,'(3a,g20.4,a,g20.4,a,g20.4,a,g20.4,a)') &
            '# <tr><td>',momentvar,'</td><td>',result(0), &
            '</td><td>', &
        result(nint(real(nboot+1)*0.025)),'...', &
        result(nint(real(nboot+1)*0.975)),'</td></tr></table>'
    else
        write(14,'(7a)') '# moment ', &
            '            mean    ', &
            '            2.5%    ', &
            '             17%    ', &
            '             50%    ', &
            '             83%    ', &
            '           97.5%    '
        write(14,'(3a,6g20.6)') '# ',momentvar,' = ',result(0), &
        result(nint(real(nboot+1)*0.025)), &
        result(nint(real(nboot+1)*0.17)), &
        result(nint(real(nboot+1)*0.5)), &
        result(nint(real(nboot+1)*0.83)), &
        result(nint(real(nboot+1)*0.975))
    end if

!   running window analysis

    rmin(0) = +3e33
    rmax(0) = -3e33
    do yr=yr1-1,yr2-nyrwindow+1
        call keepalive(yr-yr1+1,yr2-yr1-nyrwindow+2)
        n = 0
        do iens=nens1,nens2
            yr1s = max(yr1,yr)
            yr2s = min(yr2,yr+nyrwindow-1)
            do i=yr1s-1,yr2s
                do jj=j1,j2
                    j = jj
                    call normon(j,i,ii,nperyear)
                    if ( ii < yr1s .or. ii > yr2s ) goto 700
                    if ( data(j,ii,iens) < 1e33 ) then
                        n = n + 1
                        xx(n) = data(j,ii,iens)
                    end if
                700 continue
                end do
            end do
        end do
        if ( n < minnum*(j2-j1+1)*(nens2-nens1+1) ) then
            if ( lwrite ) print *,'not enough points: ',n,minnum
            goto 800
        end if
        mboot = (nboot+1)/2-1
        call bootmoment(imom,xx,yy,n,mboot,decor,result)
        if ( lwrite ) then
            print *,'with ',n,' points in ',yr1s,yr2s,' found ',result(0)
        end if
        rmin(0) = min(rmin(0),result(0))
        rmax(0) = max(rmax(0),result(0))
        write(14,'(i4,6g20.6)') yr+nyrwindow/2,result(0), &
        result(nint(real(mboot+1)*0.025)), &
        result(nint(real(mboot+1)*0.17)), &
        result(nint(real(mboot+1)*0.5)), &
        result(nint(real(mboot+1)*0.83)), &
        result(nint(real(mboot+1)*0.975))
    800 continue
    end do
    rdif(0) = rmax(0)-rmin(0)

!   MC test on significance of variations

    print '(a)','#<br>For overall significance disregarding nonzero '// &
        'skewness and curtosis for the moment.'
    if ( decor == 0 ) then
        alpha = 0
    else
        alpha = exp(-1/decor)
    end if
    do iboot=1,nboot
        call keepalive(iboot,nboot)
        do iens=nens1,nens2
            noise = gasdev(iran)
            do i=yr1,yr2
                do j=1,nperyear
                    noise = alpha*noise + sqrt(1-alpha**2)*gasdev(iran)
                    if ( data(j,i,iens) < 1e33 ) then
                        data(j,i,iens) = ave + sdev*noise
                    end if
                end do
            end do
        end do
        rmin(iboot) = +3e33
        rmax(iboot) = -3e33
        do yr=yr1,yr2-nyrwindow+1
            n = 0
            do iens=nens1,nens2
                yr1s = max(yr1,yr)
                yr2s = min(yr2,yr+nyrwindow-1)
                do i=yr1s-1,yr2s
                    do jj=j1,j2
                        j = jj
                        call normon(j,i,ii,nperyear)
                        if ( ii < yr1s .or. ii > yr2s ) goto 850
                        if ( data(j,ii,iens) < 1e33 ) then
                            n = n + 1
!                           should include higher moments RSN
                            xx(n) = data(j,ii,iens)
                        end if
                    850 continue
                    end do
                end do
            end do
            if ( n < minnum*(j2-j1+1)*(nens2-nens1+1) ) then
                if ( lwrite ) print *,'not enough points: ',n,minnum
                goto 890
            end if
            call moment(xx,n,ave1,adev1,sdev1,sd21,skew1,curt1)
            call getmoment(imom,xx,n,result(0))
            rmin(iboot) = min(rmin(iboot),result(0))
            rmax(iboot) = max(rmax(iboot),result(0))
        890 continue
        end do
        rdif(iboot) = rmax(iboot)-rmin(iboot)
    end do
    if ( decor == 0 ) then
        print '(a)','# Assuming all data points are independent'
    else
        print '(a,f5.1)','# Assuming a decorrelation length of ',decor
    end if
    print '(2a,i6,a)','# Significances are computed' &
        ,' against a ',nboot,' sample Monte Carlo.'
    if ( lweb) print '(a)','<table class="realtable" '// &
        'border=0 cellpadding=0 cellspacing=0>'// &
        '<tr><th colspan=3>Probability that the distribution is '// &
        '<br>a chance fluctuation around a constant'// &
        '<tr><th>statistic</th><th>value</th><th>'// &
        'p-value</th></tr>'
    call getsign('minimum',rmin(0),rmin(1),nboot,-1,signmin, .true. )
    call getsign('maximum',rmax(0),rmax(1),nboot,1,signmax, .true. )
    call getsign('difference',rdif(0),rdif(1),nboot,1,signdif,.true.)
    if ( lweb) print '(a)','</table>'

900 continue
    goto 999
901 write(0,*) 'runningmoments: error: too many points, increase ndatens (currently',ndatens,')'
    call exit(-1)
999 continue
end program runningmoments
