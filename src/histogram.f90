program histogram

!   make a histogram of timeseries
!   and do a maximum-likelihood fit if requested

    implicit none
    include 'param.inc'
    include 'getopts.inc'
    integer :: maxbin,maxldat,nboot
    parameter (maxbin=1000,maxldat=10000000,nboot=800)

    integer :: i,ii,j,jj,j1,j2,n,m,nperyear,nbin,month,yr,nn(maxbin) &
        ,Nless,nfit,ntot,ntype,iens,mens1,mens,l,nnn &
        ,iboot,mboot,ndecor,nfitted,it,tmax,off,yrstart,yrstop &
        ,dy,mo,iret,jmax,iensmax,nzero
    integer,allocatable :: yrs(:)
    logical :: lastvalid,lprinted,lexist
    logical,allocatable :: lf(:)
    real,allocatable :: data(:,:,:),xx(:),xs(:)
    real :: mindat,maxdat,frac &
        ,d,a,b,c,da,db,dc,xi,dxi,yy(maxbin) &
        ,x,x1,x2,s,sqrt2,snorm,f,ff,offset,xn(maxbin),xy(maxbin),df &
        ,chsq,prob,lon,lat,elev,threshold,z,tx,tx25 &
        ,tx975,t(10),t25(10),t975(10),xmax
    real :: mean(0:nboot),sd(0:nboot),skew(0:nboot),ranf
    character :: file*255,string*100,var*20,units*20,prog*100, &
        number*100,line*255,command*255,extraargs*100, &
        extraargs_*100,ids(0:nensmax)*30,assume*5
    character :: lvar*120,svar*120,history*50000,metadata(2,100)*2000
    real,external :: gammln,erf,erfcc,gammq,gammp,invcumpois,invcumgamm,serfi

    real :: scalingpower
    common /c_scalingpower/ scalingpower

    call killfile(string,line,file,0)
    allocate(yrs(0:maxldat))
    allocate(lf(maxldat))
    allocate(xx(maxldat))
    allocate(xs(maxldat))
    lwrite = .false. 
    n = command_argument_count()
    if ( n < 2 ) then
        print *,'usage: histogram {timeseries|file list prog} nbins' &
            //' hist hist|qq|gumbel|log|sqrtlog' &
            //' [fit poisson|gauss|gamma|gumbel|gev|gpd] [options]'
        call exit(-1)
    end if

!    record the exact command used in the output file
    command = '#'
    do i=0,command_argument_count()
        call get_command_argument(i,line)
        command = trim(command)//' '//trim(line)
    end do
    print '(a)',trim(command)

!   getopts comes way too late for these options...
    do i=4,command_argument_count()
        call get_command_argument(i,string)
        if ( string == 'debug' .or. string == 'lwrite' ) then
            lwrite = .true. 
        end if
        if ( string(1:9) == 'standardu' ) then
            lstandardunits = .true. 
        end if
    end do
    call get_command_argument(1,file)
    allocate(data(npermax,yrbeg:yrend,0:nensmax))
    if ( file /= 'file' ) then
!       simple data file (possibly an ensemble)
        off = 0
        call readensseriesmeta(file,data,npermax,yrbeg,yrend,nensmax, &
            nperyear,mens1,mens,var,units,lvar,svar,history,metadata, &
            lstandardunits,lwrite)
    else
        off = 2
        call readsetseries(data,ids,npermax,yrbeg,yrend,nensmax, &
            nperyear,mens1,mens,var,units,lvar,svar,history,metadata, &
            lstandardunits,lwrite)
    end if
    call get_command_argument(2+off,string)
    read(string,*) nbin
    if ( nbin > maxbin ) then
        write(0,*) 'histogram: error: increase maxbin'
        write(*,*) '# histogram error: increase maxbin'
        call exit(-1)
    end if
    i = 3 + off
    nfit = 0
    ntype = 0
    yrs(0) = 0
    assume = 'shift'
    if ( n >= i ) then
     10 continue
        call get_command_argument(i,string)
        if ( string(1:3) == 'fit' ) then
            call get_command_argument(i+1,string)
            i = i + 2
            if ( string(1:2) == 'no' ) then
                nfit = 0
            else if ( string(1:3) == 'poi' ) then
                nfit = 1
            else if ( string(1:3) == 'gau' ) then
                nfit = 2
            else if ( string(1:3) == 'gam' ) then
                nfit = 3
            else if ( string(1:3) == 'gum' ) then
                nfit = 4
            else if ( string(1:3) == 'gev' ) then
                nfit = 5
            else if ( string(1:3) == 'gpd' ) then
                nfit = 6
            else
                write(0,*) 'histogram: error: unknown distribution ',string
                write(*,*) '# histogram error: unknown distribution ',string
                call exit(-1)
            end if
            goto 10
        end if
        if ( string(1:4) == 'hist' ) then
            call get_command_argument(i+1,string)
            i = i + 2
            call tolower(string)
            if ( string(1:4) == 'hist' ) then
                ntype = 0
            else if ( string(1:3) == 'qq' ) then
                ntype = 1
            else if ( string(1:4) == 'gumb' ) then
                ntype = 2
            else if ( string(1:3) == 'log' ) then
                ntype = 3
            else if ( string(1:7) == 'sqrtlog' ) then
                ntype = 4
            else
                write(0,*) 'histogram: error: unknown plot ',string
                write(*,*) '# histogram error: unknown plot ',string
                call exit(-1)
            end if
            goto 10
        end if
        if ( string(1:6) == 'assume' ) then
            call get_command_argument(i+1,assume)
            i = i + 2
            goto 10
        end if
        call getopts(i,n,nperyear,yrbeg,yrend, .true. ,mens1,mens)
        if ( mens > 0 ) then
            print '(a,i3,a,i3)','# using ensemble members ',nens1,' to ',nens2
        end if
    end if
    if ( lchangesign ) then
        if ( ntype < 2 ) then ! no effect on histogram or QQ plot
            lchangesign = .false.
        end if
    end if
    if ( plot .and. .not. lbootstrap ) then
        mboot = 0
    else
        mboot = nboot
    end if
    yrstart = yr2
    yrstop  = yr1

!   sum series

    do iens=nens1,nens2
        if ( lsum > 1 ) then
            call sumit(data(1,yrbeg,iens),npermax,nperyear &
                ,yrbeg,yrend,lsum,oper)
        end if
        scalingpower = 1
        if ( logscale ) then
            if ( iens == nens1 ) print '(a)','# taking logarithm'
            call takelog(data(1,yrbeg,iens),npermax &
                ,nperyear,yrbeg,yrend)
            scalingpower = 0
        end if
        if ( sqrtscale ) then
            if ( iens == nens1 ) print '(a)','# taking sqrt'
            call takesqrt(data(1,yrbeg,iens),npermax &
                ,nperyear,yrbeg,yrend)
            scalingpower = scalingpower*0.5
        end if
        if ( squarescale ) then
            if ( iens == nens1 ) print '(a)','# taking square'
            call takesquare(data(1,yrbeg,iens),npermax &
                ,nperyear,yrbeg,yrend)
            scalingpower = scalingpower*2
        end if
        if ( cubescale ) then
            if ( iens == nens1 ) print '(a)','# taking square'
            call takecube(data(1,yrbeg,iens),npermax &
                ,nperyear,yrbeg,yrend)
            scalingpower = scalingpower*3
        end if
        if ( twothirdscale ) then
            if ( iens == nens1 ) print '(a)','# taking power twothird'
            call taketwothird(data(1,yrbeg,iens),npermax &
                ,nperyear,yrbeg,yrend)
            scalingpower = scalingpower*2./3.
        end if
    end do
    if ( lchangesign ) xyear = -xyear

!       apply minindx, maxindx cuts;
!       in case of a Poisson distribution round to unambiguous
!       numbers

    if ( m1 /= m2 ) then
        write(0,*) 'please specify only one starting month'
        call exit(-1)
    end if
    call getj1j2(j1,j2,m1,nperyear, .false. )
    if ( pminindx >= 0 ) then
        call getenscutoff(minindx,pminindx,data,npermax,nperyear &
            ,yrbeg,yrend,nensmax,nens1,nens2,yr1,yr2,j1,j2,0)
        if ( lwrite ) print *,'pminindx,minindx = ',pminindx,minindx
    end if
    if ( pmaxindx >= 0 ) then
        call getenscutoff(maxindx,pmaxindx,data,npermax,nperyear &
            ,yrbeg,yrend,nensmax,nens1,nens2,yr1,yr2,j1,j2,0)
        if ( lwrite ) print *,'pmaxindx,maxindx = ',pmaxindx,maxindx
    end if
    if ( nfit == 1 ) then
        if ( minindx >= 0 ) minindx = 0.5 + int(minindx+1e-5)
        if ( maxindx < 2.**31 ) maxindx = 0.5 + int(maxindx - 1e-5)
    end if
    do iens=nens1,nens2
        if ( maxindx < 1e33 .or. minindx > -1e33 ) then
            do i=yr1,yr2
                do j=1,nperyear
                    if (  data(j,i,iens) > maxindx .or. &
                        data(j,i,iens) < minindx ) then
                        data(j,i,iens) = 3e33
                    end if
                end do
            end do
        end if
    
!       minus sign for lower tails
    
        if ( lchangesign ) then
            if ( iens == nens1 ) print '(a)','# changing sign'
            call changesign(data(1,yrbeg,iens),npermax &
                ,nperyear,yrbeg,yrend)
        end if
    
!       detrend data
    
        if ( ldetrend ) then
            if ( iens == nens1 ) print '(a,i4)','# detrending data'
            call detrend(data(1,yrbeg,iens),npermax &
                ,nperyear,yrbeg,yrend,yr1,yr2,m1,m2,lsel)
        end if
    
!       differentiate data
    
        if ( ndiff /= 0 ) then
            if ( iens == nens1 ) print '(a,i4)','# taking differences/averaging ',ndiff
            call diffit(data(1,yrbeg,iens),npermax,nperyear &
                ,yrbeg,yrend,ndiff)
        end if
    
!       anomalies
    
        if ( anom ) then
            if ( iens == nens1 ) print '(a)','# taking anomalies'
            call anomal(data(1,yrbeg,iens),npermax,nperyear &
                ,yrbeg,yrend,yr1,yr2)
        end if
    end do

!   anomalies wrt ensemble mean

    if ( lensanom .and. nens1 /= nens2 ) then
        if ( lwrite ) print '(a)','# taking anomalies wrt ensemble mean'
        call anomalensemble(data,npermax,nperyear,yrbeg, &
            yrend,yr1,yr2,nens1,nens2)
    end if

!   compute cut-off for GPD - has to be done after all the operations on the series
!   to get a consistent computation (eg otherwise it would vary by season...)

    if ( nfit == 6 .and. pmindata >= 0 ) then
        if ( lchangesign ) then
            call getenscutoff(mindata,100-pmindata,data,npermax &
                ,nperyear,yrbeg,yrend,nensmax,nens1,nens2,yr1 &
                ,yr2,j1,j2,0)
        else
            call getenscutoff(mindata,pmindata,data,npermax &
                ,nperyear,yrbeg,yrend,nensmax,nens1,nens2,yr1 &
                ,yr2,j1,j2,0)
        end if
        if ( lwrite ) print *,'pmindata,mindata = ',pmindata,mindata
    end if

!   compute min, max

    do month=m1,m2
        print '(a,i2)','# month ',month
        if ( month == 0 ) then
            j1 = 1
            j2 = nperyear
        else
            j1 = month
            j2 = month + lsel - 1
            if ( nperyear > 12 ) then
                call month2period(j1,nperyear,1)
                call month2period(j2,nperyear,0)
            end if
        end if
        ntot = 0
        do iboot=0,mboot
            mean(iboot) = 0
            sd(iboot) = 0
            skew(iboot) = 0
        end do
        mindat = +3e33
        maxdat = -3e33
        offset = 0
        tx = 3e33
        tx25 = 3e33
        tx975 = 3e33
        t = 3e33
        t25 = 3e33
        t975 = 3e33
        ! first find the largest value in the year under study
        xmax = -3e33
        do iens=nens1,nens2
            do yr=yr1-1,yr2
                do jj=j1,j2
                    j = jj
                    call normon(j,yr,i,nperyear)
                    if ( i >= yr1 .and. i <= yr2) then
                        if ( i == yr2a .and. data(j,i,iens) < 1e33 ) then
                            if ( xmax < data(j,i,iens) ) then
                                jmax = j
                                iensmax = iens
                                xmax = data(j,i,iens)
                            end if ! largest value?
                        end if ! year to be taken out?
                    end if ! in range?
                end do  ! jj
            end do ! yr
        end do ! iens
        do iens=nens1,nens2
            lastvalid = .false. 
            do yr=yr1-1,yr2
                if ( j1 /= j2 .and. j2-j1+1 /= nperyear ) then
                    lastvalid = .false. 
                end if
                do jj=j1,j2
                    j = jj
                    call normon(j,yr,i,nperyear)
                    if ( i >= yr1 .and. i <= yr2 ) then
                        if ( i == yr2a .and. j == jmax ) then
                            ! set xyear to the largest value
                            if ( iens == iensmax ) then
                                xyear = data(j,i,iens)
                            end if
                            ! and set all ensemble members to undef
                            ! to avoid possibly dependent data
                            data(j,i,iens) = 3e33
                        end if
                        if ( data(j,i,iens) < 1e33 ) then
                            ntot = ntot + 1
                            if ( ntot > maxldat ) then
                                print *,'# histogram error increase maxldat'
                                write(0,*) 'histogram: error: increase maxldat'
                                call exit(-1)
                            end if
                            xx(ntot) = data(j,i,iens)
                            call getdymo(dy,mo,j,nperyear)
                            yrs(ntot) = 10000*i + 100*mo + dy
                            if ( nperyear > 1000 ) then
                                yrs(ntot) = 100*yrs(ntot) + mod(ntot,nint(nperyear/366.))
                            end if
                            if ( nens1 /= nens2 ) then
                                yrs(ntot) = 100*yrs(ntot) + iens
                            end if
                            offset = offset + xx(ntot)
                            yrstart = min(yrstart,i)
                            yrstop  = max(yrstop,i)
                            lf(ntot) = .not. lastvalid
                            lastvalid = .true. 
                            if ( .false. .and. lwrite ) print *,'xx(',ntot,') = ',xx(ntot)
                            if ( lwrite .and. lf(ntot) ) print *,'boundary at ',ntot,j,i,iens
                        else
                            lastvalid = .false. 
                        end if
                    else
                        lastvalid = .false. 
                    end if
                end do
            end do
        end do
        if ( ntot == 0 ) then
            write(*,'(a)') '# histogram: no valid points'
            stop
        end if
        offset = offset/ntot
        if ( lwrite ) print '(a,g20.10)','# offset = ',offset
        do i=1,ntot
            x = xx(i) - offset
            mean(0) = mean(0) + x
            sd(0) = sd(0) + x**2
            skew(0) = skew(0) + x**3
            mindat = min(mindat,xx(i))
            maxdat = max(maxdat,xx(i))
        end do
    
!       enough points?
    
        if ( ntot < minnum ) then
            print *,'not enough points'
            goto 999
        end if
    
!       compute errors with a bootstrap
    
        decor = max(decor,real(lsum)-1)
        if ( j1 == j2 ) then
            ndecor = 1+int(decor/nperyear)
        else
            ndecor = 1+int(decor)
        end if
        call print_bootstrap_message(max(1,ndecor),j1,j2)
        do iboot=1,mboot
            call keepalive1('Bootstrapping',iboot,mboot)
            n = ntot/ndecor
            do i=1,n
                ii = 0
                200 continue
                ii = ii + 1
                if ( ii > 100 ) then
                    write(*,'(a,i4,a,i6)') '# histogram: error: '// &
                        'cannot find blocks of size ',ndecor &
                        ,' for moving block bootstrap ',iboot
                    goto 300
                end if
                call random_number(ranf)
                j = 1+min(ntot-ndecor,int((ntot-ndecor)*ranf))
                ! check whether there was a break in this segment
                do jj=1,ndecor-1
                    if ( lf(j+jj) ) then
                        if ( lwrite ) print *,'crossing boundary at ',j+jj,ii,iboot
                        goto 200
                    end if
                end do
                do jj=0,ndecor-1
                    x = xx(j+jj) - offset
                    mean(iboot) = mean(iboot) + x
                    sd(iboot) = sd(iboot) + x**2
                    skew(iboot) = skew(iboot) + x**3
                end do
            end do
        end do
    300 continue
        mboot = iboot-1
        do iboot=0,mboot
            if ( iboot == 0 ) then
                n = ntot
            else
                n = ndecor*(ntot/ndecor)
            end if
            mean(iboot) = mean(iboot)/n
            sd(iboot) = sd(iboot)/n - mean(iboot)**2
            if ( sd(iboot) >= 0 ) then
                sd(iboot) = sqrt(sd(iboot))
            else
                print '(a,2g20.10,i6)','# error: variance<0 ' &
                    ,mean(iboot),sd(iboot),iboot
                do i=1,n
                    print '(a,i6,g20.10)','# ',i,xx(i)
                end do
                sd(iboot) = 0
            end if
            skew(iboot) = skew(iboot)/n - 3*sd(iboot)**2 &
                *mean(iboot) - mean(iboot)**3
            if ( sd(iboot) > 0 ) then
                skew(iboot) = skew(iboot)/sd(iboot)**3
            else
                skew(iboot) = 3e33
            end if
            mean(iboot) = mean(iboot) + offset
        end do
        call nrsort(mboot,mean(1))
        call nrsort(mboot,sd(1))
        call nrsort(mboot,skew(1))
        print '(a,i16)','# N:         ',ntot
        call print_bootstrap_message(ndecor,j1,j2)
        if ( lweb ) then
            if ( namestring /= ' ' ) then ! currently not used
                print '(4a)','# <tr><th colspan=3>',trim(namestring),'</th></tr>'
            end if
            if ( abs(nint(confidenceinterval)-confidenceinterval) < 0.001 ) then
                print '(2a,i2,a)','# <tr><th>parameter</th><th>', &
                    'value&plusmn;2&sigma;</th><th>', &
                    nint(confidenceinterval),'% CI</th></tr>'
            else
                print '(2a,f6.3,a)','# <tr><th>parameter</th><th>', &
                    'value&plusmn;2&sigma;</th><th>', &
                    confidenceinterval,'% CI</th></tr>'
            end if
            print '(a,i9,a)','# <tr><td>N:</td><td>',ntot,'</td><td>&nbsp;</td></tr>'
        else
            print '(8a)',       '#        ', &
                '                ', &
                '            2.5%', &
                '             16%', &
                '             50%', &
                '             84%', &
                '           97.5%', &
                '          +/-95%'
            print '(a,i9)','# N:         ',ntot
        end if
        call printvalerr('# mean:      ',mean,mboot,plot,11,lweb,lchangesign)
        call printuntransf(mean(0))
        call printvalerr('# s.d.(n):   ',sd,mboot,plot,11,lweb,.false. )
        call printval   ('# s.d.(n-1): ',sd(0)*sqrt(real(ntot)/real(ntot-1)),plot,11,lweb, .false. )
        call printvalerr('# skew:      ',skew,mboot,plot,11,lweb,lchangesign)
        call printval   ('# min:       ',mindat,plot,11,lweb,lchangesign)
        call printuntransf(mindat)
        call printval   ('# max:       ',maxdat,plot,11,lweb,lchangesign)
        call printuntransf(maxdat)
    
!       fit to distribution
    
        xs(1:ntot) = xx(1:ntot)
        if ( nfit == 0 ) then
!           no fit requested
            snorm = 1
        else if ( nfit == 1 ) then
!           Poisson distribution
            call fitpoi(xx,ntot,mean(0),a,j1,j2,lweb,ntype &
                ,lchangesign,yr2a,xyear,t,t25,t975,tx &
                ,tx25,tx975,confidenceinterval,.true.,.true.,lwrite)
            call poisnorm(a,snorm)
            nfitted = 2
        else if ( nfit == 2 ) then
!           Gaussian distribution
            call fitgau(xx,ntot,mean(0),sd(0),a,b,minindx,maxindx, &
                ntype,j1,j2,yr2a,xyear,t,t25,t975,tx,tx25,tx975, &
                confidenceinterval,.true.,.true.,lweb,lchangesign, &
                lwrite)
            call gausnorm(a,b,snorm)
            nfitted = 3
        else if ( nfit == 3 ) then
!           Gamma distribution
            if ( units(1:4) == 'mm/d' ) then
                ! impose the ETCCDI cut-off of 1 mm/dy for wet days.
                do i=1,ntot
                    if ( xx(i) > 0 .and. xx(i) < 1 ) xx(i) = 0
                end do
            end if
            nzero = 0
            do i=1,ntot
                if ( xx(i) == 0 ) nzero = nzero + 1
            end do
            call fitgam(xx,ntot,mean(0),sd(0),a,b,j1,j2,lweb,ntype &
                ,lchangesign,yr2a,xyear,t,t25,t975,tx &
                ,tx25,tx975,confidenceinterval,.true.,.true.,lwrite)
            if ( minindx > -1e33 ) then
                if ( maxindx < 1e33 ) then
                    snorm = gammp(a,maxindx/b) - gammp(a,minindx/b)
                else
                    snorm = gammq(a,minindx/b)
                end if
            else
                if ( maxindx < 1e33 ) then
                    snorm = gammp(a,maxindx/b)
                else
                    snorm = 1
                end if
            end if
            nfitted = 3
        else if ( nfit == 4 ) then
!           Gumbel distribution
            if ( nperyear > 12 .and. oper /= 'a' .and. oper /= 'i' ) then
                write(0,*) 'Warning: data do not appear to be block maxima or minima.<p>'
            end if
            call fitgum(xx,ntot,mean(0),sd(0),a,b,j1,j2,lweb,ntype, &
                lchangesign,yr2a,xyear,t,t25,t975,tx, &
                tx25,tx975,confidenceinterval,.true.,.true.,lwrite)
            call gumbnorm(a,b,snorm)
            nfitted = 3
            if ( lchangesign ) a = -a
        else if ( nfit == 5 ) then
!           GEV distribution
            if ( nperyear > 12 .and. oper /= 'a' .and. oper /= 'i' ) then
                write(0,*) 'Warning: data do not appear to be block maxima or minima.<p>'
            end if
            call fitgev(xx,ntot,mean(0),sd(0),a,b,xi,j1,j2,lweb &
                ,ntype,lchangesign,yr2a,xyear,t,t25,t975,tx &
                ,tx25,tx975,restrain,confidenceinterval,.true. &
                ,.true.,lwrite)
            call gevnorm(a,b,xi,snorm)
            nfitted = 4
        else if ( nfit == 6 ) then
!           GPD distribution
            call fitgpd(xx,ntot,mean(0),sd(0),b,xi,j1,j2,lweb, &
                ntype,lchangesign,pmindata,threshold,yr2a,xyear, &
                t,t25,t975,tx,tx25,tx975,restrain,assume, &
                confidenceinterval,.true.,.true.,lwrite)
            snorm = 1
            nfitted = 3
        else
            write(0,*) 'histogram: error: unknown distribution ',nfit
            call exit(-1)
        end if
        if ( ntype == 0 ) then
        
!           adjust binsize to nice numbers (heuristics)
        
!           if the minimum is close to zero make exactly zero
            if ( mindat > 0 .and. mindat < maxdat/3 ) then
                if ( lwrite ) print *,'adjusted mindat from ',mindat,' to 0'
                mindat = 0
            end if
            if ( maxdat < 0 .and. maxdat > mindat/3 ) then
                if ( lwrite ) print *,'adjusted maxdat from ',maxdat,' to 0'
                maxdat = 0
            end if
!           if the minimum and maxmimum are similar in absolute value
!		    make them equal
            if ( mindat < 0 .and. maxdat > 0 .and. abs(mindat+maxdat)/2 < maxdat/3 ) then
                if ( mindat < -maxdat ) then
                    if ( lwrite ) print *,'adjusted maxdat from ',maxdat,' to -mindat ',-mindat
                    maxdat = -mindat
                else
                    if ( lwrite ) print *,'adjusted mindat from ',mindat,' to -maxdat ',-maxdat
                    mindat = -maxdat
                end if
            end if
!           give the user a chnce to specify the end of the scale
            if ( maxdata < 1e33 ) then
                if ( lwrite ) print *,'adjusted maxdat from ',maxdat &
                    ,' to user-defined value ',maxdata
                maxdat = maxdata
                mindat = 0
                d = maxdat/nbin
            else
                d = (maxdat-mindat)/(nbin-1)
                s = d
                if ( d >= 1 ) then
                    i = int(log10(d))
                else if ( d > 0 ) then
                    i = int(log10(d)) - 1
                else
                    i = 1
                end if
                d = d/10.**i ! now d is between 1 and 10
                if ( d > 5 ) then
                    d = 10
                else if ( d > 2 ) then
                    d = 5
                else if ( d > 1 ) then
                    d = 2
                else
                    d = 1
                end if
                d = d*10.**i
                if ( lwrite ) then
                    print '(a,f16.4,a,f16.4)','# binsize rounded from ',s,' to ',d
                end if
                s = mindat
                if ( mindat >= 0 ) then
                    mindat = d*(int(mindat/d))
                else if ( mindat /= -maxdat ) then
                    mindat = d*(-1+int(mindat/d))
                else
                    mindat = -nbin*d/2
                end if
                if ( lwrite ) then
                    print '(a,f16.4,a,f16.4)','# minimum rounded from ',s,' to ',mindat
                end if
            end if
        
!           fill histogram
        
            do i=1,nbin
                nn(i) = 0
            end do
            nless = 0
            do yr=yr1-1,yr2
                do jj=j1,j2
                    j = jj
                    call normon(j,yr,i,nperyear)
                    if ( i >= yr1 .and. i <= yr2 ) then
                        do iens=nens1,nens2
                            if ( data(j,i,iens) < 1e33 ) then
                                n = 1 + int((data(j,i,iens)-mindat)/d)
                                if ( n < 0 .or. n > nbin ) then
                                    write(*,*)'# histogram warning: n>nbin:',n,nbin
                                else
                                    nn(n) = nn(n) + 1
                                    if ( data(j,i,iens) < mean(0) ) nless= nless + 1
                                end if
                            end if
                        end do
                    end if
                end do
            end do
        
!           fill fitted curve array
        
            if ( nfit == 0 ) then
                ! no fit requested
                do i=1,nbin
                    yy(i) = 0
                end do
            else if ( nfit == 1 ) then
                ! Poisson distribution
                do i=1,nbin
                    x1 = mindat + d*(i-1)
                    yy(i) = 0
                    do j=1,nint(d)
                        x2 = x1 + (j-1)
                        yy(i) = yy(i) + ntot*exp(x2*log(a)-a-gammln(1+x2))
                    end do
                    yy(i) = yy(i)/snorm
                end do
            else if ( nfit == 2 ) then
                ! Gaussian distribution
                sqrt2 = sqrt(2.)
                do i=1,nbin
                    x1 = (mindat + d*(i-1) - a)/b
                    x2 = (mindat + d*i     - a)/b
                    if ( x2 > 0 ) then
                        yy(i) = ntot*(erfcc(x1/sqrt2) - erfcc(x2/sqrt2))/2
                    else
                        yy(i) = ntot*(erfcc(-x2/sqrt2) - erfcc(-x1/sqrt2))/2
                    end if
                    yy(i) = yy(i)/snorm
                end do
            else if ( nfit == 3 ) then
                ! Gamma distribution
                do i=1,nbin
                    x1 = (mindat + d*(i-1))/abs(b)
                    x2 = (mindat + d*i    )/abs(b)
                    if ( abs(x2) < a-1 ) then
                        s = gammp(a,abs(x2)) - gammp(a,abs(x1))
                    else
                        s = gammq(a,abs(x1)) - gammq(a,abs(x2))
                    end if
                    yy(i) = (ntot-nzero)*abs(s)/snorm
                end do
            else if ( nfit == 4 ) then
                ! Gumbel distribution
                do i=1,nbin
                    x1 = (mindat + d*(i-1) - a)/b
                    x2 = (mindat + d*i     - a)/b
                    yy(i) = ntot*(exp(-exp(-x2))-exp(-exp(-x1)))/snorm
                end do
            else if ( nfit == 5 ) then
                ! GEV distribution
                do i=1,nbin
                    x1 = (mindat + d*(i-1) - a)/b
                    x2 = (mindat + d*i     - a)/b
                    if ( xi == 0 ) then
                        yy(i) = ntot*(exp(-exp(-x2))-exp(-exp(-x1)))/snorm
                    else
                        if ( 1+xi*x1 <= 0 .or. 1+xi*x2 <= 0 ) then
                            yy(i) = -999.9
                        else
                            yy(i) = ntot*(exp(-(1+xi*x2)**(-1/xi)) - exp(-(1+xi*x1)**(-1/xi)))
                        end if
                    end if
                end do
            else if ( nfit == 6 ) then
                ! GPD distribution
                if ( lchangesign ) then
                    nnn = 0
                    do i=1,nbin
                        x2 = -(mindat + d*i + mindata)/b
                        print *,i,x2
                        if ( x2 > 0 ) nnn = nnn + nn(i)
                    end do
                else
                    nnn = 0
                    do i=1,nbin
                        x1 = (mindat + d*(i-1) - mindata)/b
                        if ( x1 > 0 ) nnn = nnn + nn(i)
                    end do
                end if
                do i=1,nbin
                    if ( lchangesign ) then
                        x1 = -(mindat + d*(i-1) + mindata)/b
                        x2 = -(mindat + d*i     + mindata)/b
                    else
                        x1 = (mindat + d*(i-1) - mindata)/b
                        x2 = (mindat + d*i     - mindata)/b
                    end if
                    if ( ( .not. lchangesign .and. x1 < 0) .or. (lchangesign .and. x2 < 0) ) then
                        yy(i) = -999.9
                    else if ( abs(xi) <= 1e-4 ) then
                        yy(i) = nnn*(exp(-x1+0.5*x1**2*xi) - exp(-x2+0.5*x2**2*xi))/snorm
                    else
                        if ( 1+xi*x1 <= 0 .or. 1+xi*x2 <= 0 ) then
                            yy(i) = -999.9
                        else
                            yy(i) = nnn*((1+xi*x1)**(-1/xi) - (1+xi*x2)**(-1/xi))/snorm
                        end if
                    end if
                    if ( lwrite ) write(*,*) x1,x2,yy(i),snorm,nnn
                end do
            else
                write(0,*) 'histogram: error: unknown distribution ',nfit
                call exit(-1)
            end if
        
!           compute \chi^2 and P of fit
        
            if ( nfit /= 0 ) then
                ! leave out trivial points with N=0, n~0 in the tails
                n = 0
                nnn = 0
                do i=1,nbin
                    if ( nfit == 6 ) then
                        ! only count the points above the threshold...
                        nnn = nnn + nn(i)
                        if ( lchangesign ) then
                            x2 = mindat + d*(i-1)
                            if ( lwrite ) then
                                print *,'i,x2,nn(i),nnn,yy(i) = ',i,x2,nn(i),nnn,yy(i)
                                print *,'x2,mindata = ',x2,mindata
                            end if
                            if ( x2 < -mindata ) cycle
                        else
                            x1 = mindat + d*i
                            if ( lwrite ) then
                                print *,'i,x1,nn(i),nnn,yy(i) = ',i,x1,nn(i),nnn,yy(i)
                                print *,'x1,mindata = ',x1,mindata
                            end if
                            if ( x1 < mindata ) cycle
                        end if
                    end if
                    if ( yy(i) /= -999.9 .and. ( yy(i) > 0.1 .or. nn(i) > 0 .or. n > 0 ) ) then
                        if ( yy(i) /= 0 .or. nn(i) /= 0 ) then
                            if ( nn(i) > 0 .and. yy(i) <= 0 ) then
                                write(0,*) 'histogram: error: events with p=0 occurred: ',i,nn(i)
                            end if
                            n = n + 1
                            xn(n) = nn(i)
                            xy(n) = yy(i)
                        end if
                    else if ( n == 0 ) then
                        if ( lwrite ) print *,'n = 0!'
                    end if
                end do
                do i=n,1,-1
                    if ( lwrite ) print *,'looking for non-zero bin: ',i,xy(i),xn(i)
                    if ( xy(i) > 0.1 .or. xn(i) > 0 ) goto 410
                end do
            410 continue
                n = i
                if ( lwrite ) then
                    do i=1,n
                        print '(a,i4,2f9.2,g10.2)','# ',i,xn(i),xy(i),(xn(i)-xy(i))**2/xy(i)
                    end do
                end if
                ! Numerical Recipes, p615
                if ( n > nfitted ) then
                    call chsone(xn,xy,n,nfitted,df,chsq,prob)
                    if ( plot ) write(11,'(g16.4)') chsq/df
                    if ( lweb ) then
                        print '(a,g9.2,a,i4,a,g9.3,a,f6.4)', '# <tr><td>&chi;^2/df</td><td>' &
                            ,chsq,'/',nint(df),' =',chsq/df,'</td><td>p=',prob
                    else
                        print '(a,g9.2,a,i4,a,g9.3)','# chi^2/df:',chsq,'/',nint(df),' =',chsq/df
                        print '(a,f8.4,a)','# probability:   ',100*prob, &
                            '% (this is the probability that the observed distribution was drawn ' &
                            //'from the fitted one)'
                    end if
                    if ( plot ) write(11,'(g16.4)') prob
                else
                    if ( lwrite ) print *,'# not enough non-zero bins to determine \chi^2: ',n,nfitted
                end if
            end if
        
!           print out
        
            n = 0
            s = 0
            do i=1,nbin
                n = n + nn(i)
                if ( yy(i) /= -999.9 ) then
                    s = s + yy(i)
                else
                    s = n
                end if
                x = mindat+ d*(i-0.5)
                if ( lchangesign ) then
                    x = -x
                end if
                print '(i4,g12.4,2i9,2f12.2)',i,x,nn(i),n,yy(i),s
            end do

        else if ( ntype == 1 ) then
        
!           QQ plot
        
            call nrsort(ntot,xx)
            do i=1,ntot
                f = real(i)/real(ntot+1)
                if ( nfit == 0 ) then
                    ! no fit requested
                    write(0,*) 'histogram: error: can only make QQ plot when fitting to a distribution'
                    write(*,*) '# histogram error: can only make QQ plot when fitting to a distribution'
                    call exit(-1)
                else if ( nfit == 1 ) then
                    ! Poisson distribution - only the last point of a bin makes sense
                    if ( xx(i) /= xx(i+1) ) then
                        f = snorm*f
                        if ( minindx > 0 ) then
                            f = f + gammq(minindx+0.5,a)
                        end if
                        s = invcumpois(f,1-f,a)
                    else
                        s = 3e33
                    end if
                else if ( nfit == 2 ) then
                    ! Gaussian distribution
                    sqrt2 = sqrt(2.)
                    ff = 2*snorm*f
                    if ( minindx > -1e33 ) then
                        ff = ff + erf((minindx-a)/(sqrt2*b))
                    else
                        ff = ff - 1
                    end if
                    ! netlib routine
                    z = serfi(ff)
                    s = a + sqrt2*b*z
                else if ( nfit == 3 ) then
                    ! Gamma distribution
                    if ( i <= nzero ) then
                        s = 0
                    else
                        f = (ntot*f - nzero)/(ntot - nzero)
                        f = snorm*f
                        if ( minindx > 0 ) then
                            f = f + gammp(a,minindx/b)
                        end if
                        s = invcumgamm(f,1-f,a,abs(b))
                    end if
                else if ( nfit == 4 ) then
                    ! Gumbel distribution
                    s = snorm*f
                    if ( minindx > -1e33 ) then
                        s = s + exp(-exp(-(minindx-a)/b))
                    end if
                    s = a - b*log(-log(s))
                else if ( nfit == 5 ) then
                    ! GEV distribution
                    s = snorm*f
                    if ( minindx > -1e33 ) then
                        s = s + exp(-(1+(minindx-a)/b)**(-1/xi))
                    end if
                    if ( xi == 0 ) then
                        s = a - b*log(-log(s))
                    else
                        s = a + b/xi*((-log(s))**(-xi)-1)
                    end if
                else if ( nfit == 6 ) then
                    ! GPD distribution
                    write(0,*) 'QQ plot for GPD not yet ready'
                    call exit(-1)
                else
                    write(0,*) 'histogram: error: unknown distribution ',nfit
                    call exit(-1)
                end if
                if ( s < 1e33 ) then
                    x = xx(i)
                    if ( lchangesign ) then
                        if ( x /= -999.9 ) x = -x
                    end if
                    print '(i8,2g22.6)',i,x,s
                end if
            end do
        else if ( ntype == 2 .or. ntype == 3 .or. ntype == 4 ) then

!           CDF

            if ( plot ) then
                do i=1,10
                    write(11,'(3g16.4)') t(i),t25(i),t975(i)
                end do
                write(11,'(3g16.4)') tx,tx25,tx975
            end if
            if ( lchangesign ) xyear = -xyear
            if ( lchangesign .and. nfit == 6 ) then ! heuristic
                a = threshold
                mindata = threshold
                b = abs(b)
            end if
            if ( lchangesign .and. nfit == 5 ) then ! heuristic
                a = -a
                b = -b
            end if
            if ( lchangesign .and. nfit == 3 ) then ! heuristic
                b = -b
            end if
            if ( nfit == 3 .and. units(1:4) == 'mm/d' ) then
                ! impose the ETCCDI cut-off of 1 mm/dy for wet days.
                do i=1,ntot
                    if ( xx(i) > 0 .and. xx(i) < 1 ) xx(i) = 0
                    if ( xs(i) > 0 .and. xs(i) < 1 ) xs(i) = 0
                end do                
            end if
            frac = 1
            call plot_ordered_points(xx,xs,yrs,ntot,ntype,nfit, &
                frac,a,b,xi,j1,j2,1,1,minindx,mindata,pmindata, &
                yr2a,xyear,snorm,lchangesign,lwrite, .true. )
        else
            write(0,*) 'Unknown plot type ',ntype
            write(*,*) 'Unknown plot type ',ntype
            call exit(-1)
        end if
    end do
999 continue
    call printmetadata(6,file,' ','histogram of',history,metadata)
    call savestartstop(yrstart,yrstop)
end program

