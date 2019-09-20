subroutine getopts(iarg1,iarg2,nperyear,yrbeg,yrend,loutin,mens1,mens)

!   process options common to correlat*.[fF]
!   the switches are in common in getopts.inc

    implicit none
    integer :: iarg1,iarg2,nperyear,yrbeg,yrend,mens1,mens
    logical :: loutin
    include 'getopts.inc'
    integer :: i,j,k,l,lskip,idiff,iret
    real :: s,minfacinit
    logical :: lout,lexist,lopen
    character line*511,string*511

!   init

    lout = loutin
    logscale = .FALSE. 
    sqrtscale = .FALSE. 
    lchangesign = .FALSE. 
    lstandardunits = .FALSE. 
    lallobs = .FALSE. 
    logfield = .FALSE. 
    sqrtfield = .FALSE. 
    squarescale = .FALSE. 
    twothirdscale = .FALSE. 
    lnormsd = .FALSE. 
    lnomissing = .FALSE. 
    lrank = .FALSE. 
    lfitnoise = .FALSE. 
    ldetrend = .FALSE. 
    debias = 0
    normalization = 1
    biasmul = 1
    biasadd = 0
    lincludelast = .false.
    biasrt = -3e33
    add_option = 0
    lks = .FALSE. 
    lconting = .FALSE. 
    lbootstrap = .FALSE. 
    lrmse = .FALSE. 
    lmae = .FALSE. 
    irunvar = 0
    noisetype = 0
    anom = .FALSE. 
    lensanom = .FALSE. 
    composite = .FALSE. 
    lrandom = .FALSE. 
    idiff = 0
    ndiff = 0
    ndiff2 = 0
    lnooverlap = .FALSE. 
    ncrossvalidate = 0
    bbfile = ' '
    nfittime = 0
    fitfunc = 1
    fix2 = .FALSE. 
    lag1 = 0
    lag2 = 0
    lsum = 1
    lsum2 = -1
    mdiff = 0
    mdiff2 = -1
    nyrwindow = 0
    lsel = 1
    decor = 0
    m1 = 0
    if ( nperyear == 1 ) then
        m2 = 0
    else
        m2 = min(12,nperyear)
    endif
    yr1 = yrbeg
    yr2 = yrend
    yr1a = yrbeg
    yr2a = yrend + 1
    yr2b = -9999
    lstartstop = .FALSE. 
    minindx = -2e33
    maxindx = +2e33
    mindata = -2e33
    maxdata = +2e33
    pminindx = -1
    pmaxindx = -1
    pmindata = -1
    pmaxdata = -1
    xyear = 3e33
    minfacinit = -1
    minfac = minfacinit
    minfacsum = minfacinit
    minnum = -1
    lon1 = 0
    lon2 = 360
    lat1 = -90
    lat2 = +90
    lev1 = -3e33
    lev2 = 3e33
    altlon1 = 0
    altlon2 = 360
    altlat1 = -90
    altlat2 = +90
    altlev1 = -3e33
    altlev2 = 3e33
    avex = 1
    avey = 1
    altavex = 1
    altavey = 1
    intertype = 0
    confidenceinterval = 95
    nens1 = mens1
    nens2 = mens
    lmakeensfull = .FALSE. 
    restrain = 0
    lead1 = 1
    lead2 = 1
    leads = .FALSE. 
    lsubtract = .FALSE. 
    lwrite = .FALSE. 
    dump = .FALSE. 
    plot = .FALSE. 
    lweb = .FALSE. 
    namestring = ' '
    do i=1,19
        lincl(i) = .FALSE. 
    enddo
    indxuse = 9
    strindx(1) = 'SOI'
    strindx(2) = 'NINO12'
    strindx(3) = 'NINO3'
    strindx(4) = 'NINO4'
    strindx(5) = 'NINO3.4'
    strindx(6) = 'NAO'
    strindx(7) = 'sunspot'
    strindx(8) = 'sslength'
    strindx(9) = 'time'
    call getenv('REMOTE_ADDR',line)
    if ( line /= ' ' ) then
        lweb = .TRUE. 
        lout = .FALSE. 
    endif

!       loop over arguments

    lskip = 0
    do i=iarg1,iarg2
        if ( lskip > 0 ) then
            lskip = lskip - 1
            goto 50
        endif
        call get_command_argument(i,line)
        if ( lwrite ) print *,'getopts: parsing ',i,trim(line)
        if ( line(1:8) == 'logfield' ) then
            logfield = .TRUE. 
            if ( lout ) print '(a)','# set logscale on for field'
        elseif ( line(1:3) == 'log' ) then
            logscale = .TRUE. 
            if ( lout ) print '(a)','# set logscale on for timeseries'
        elseif ( line(1:9) == 'sqrtfield' ) then
            sqrtfield = .TRUE. 
            if ( lout ) print '(a)','# set sqrt on for field'
        elseif ( line(1:4) == 'sqrt' ) then
            sqrtscale = .TRUE. 
            if ( lout ) print '(a)','# set sqrt on for time series'
        elseif ( line(1:6) == 'square' ) then
            squarescale = .TRUE. 
            if ( lout ) print '(a)','# set square on'
        elseif ( line(1:6) == 'cube' ) then
            cubescale = .TRUE. 
            if ( lout ) print '(a)','# set square on'
        elseif ( line(1:8) == 'twothird' ) then
            twothirdscale = .TRUE. 
            if ( lout ) print '(a)','# set two-third power on'
        elseif ( line(1:10) == 'changesign' ) then
            lchangesign = .TRUE. 
            if ( lout ) print '(a)','# change sign of data'
        elseif ( line(1:6) == 'normsd' ) then
            lnormsd = .TRUE. 
            if ( lout ) print '(a)','# normalize to standard deviation'
        elseif ( line(1:6) == 'nomiss' ) then
            lnomissing = .TRUE. 
            if ( lout ) print '(a)','# no missing data'
        elseif ( line(1:4) == 'rank' ) then
            lrank = .TRUE. 
            if ( lout ) print '(a)','# set rank correllations on'
        elseif ( line(1:8) == 'fitnoise' ) then
            lfitnoise = .TRUE. 
            if ( lout ) print '(a)','# fit to noise model'
        elseif ( line(1:7) == 'fitfunc' ) then
            call get_command_argument(i+1,line)
            if ( line(1:3) == 'lin' ) then
                fitfunc = 1
            elseif ( line(1:4) == 'quad' .OR. line(1:5) == 'parab' ) then
                fitfunc = 2
            elseif ( line(1:3) == 'cub' ) then
                fitfunc = 3
            else
                read(line,*,err=927) fitfunc
            endif
            if ( lout ) print '(a,i2)','# fit to polynomial of degree ',fitfunc
            lskip = 1
        elseif ( line(1:7) == 'fittime' ) then
            if ( lout ) print '(a)','# fitting time derivative with relaxation'
            call get_command_argument(i+1,line)
            if (  ichar(line(1:1)) >= ichar('0') .AND. &
                  ichar(line(1:1)) <= ichar('9') ) then
                read(line,*,err=919) nfittime
                if ( nfittime <= 0 ) goto 919
                if ( lout ) print '(a,i3,a)','# computing time derivative with ',nfittime,' points'
                lskip = 1
            else
                nfittime = 1
            endif
        elseif ( line(1:3) == 'dec' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=913) decor
            decor = abs(decor)
            if ( lout ) print '(a,f6.2)','# computing significances with decorrelation ',decor
            lskip = 1
        else if ( line(1:4) == 'add_' ) then
            if ( line(1:8) == 'add_anom' ) then
                add_option = 0
                if ( lwrite ) print '(a)','# filling missing data with anomalies'
            else if ( line(1:8) == 'add_clim' ) then
                add_option = 1
                if ( lwrite ) print '(a)','# filling missing data with climatology'
            else if ( line(1:9) == 'add_trend' ) then
                add_option = 2
                if ( lwrite ) print '(a)','# filling missing data with climatology plus trend'
            else if ( line(1:8) == 'add_pers' ) then
                add_option = 3
                if ( lwrite ) print '(a)','# filling missing data with persistence'
            else if ( line(1:8) == 'add_damp' ) then
                add_option = 4
                if ( lwrite ) print '(a)','# filling missing data with damped persistence'
            else
                write(0,*) 'getopts: unknown option ',trim(line)
            end if
        elseif ( line(1:2) == 'ks' ) then
            lks = .TRUE. 
            call get_command_argument(i+1,line)
            if ( line(1:4) == 'plot' ) then
                kscut = -3e33
                if ( lout ) print *,'Performing KS test, plot in dump.dat'
                dump = .TRUE. 
                open(10,file='dump.dat')
            else
                read(line,*,err=911) kscut
                if ( lout ) print *,'Performing KS test at cut-off ',kscut
            endif
            lskip = 1
        elseif ( line(1:6) == 'contin' ) then
            lconting = .TRUE. 
            if ( lout ) print '(2a)','# computing 3x3 contingency table with cuts min/max index/data'
        elseif ( line(1:3) == 'rms' ) then
            lrmse = .TRUE. 
            if ( lout ) print '(2a)','# computing rms error, bias, ratio variance'
        elseif ( line(1:3) == 'mae' ) then
            lmae = .TRUE. 
            if ( lout ) print '(2a)','# computing MAE, bias, ratio variance'
        elseif ( line(1:4) == 'boot' ) then
            lbootstrap = .TRUE. 
            if ( lout ) print '(a)','# Bootstrapping to get 1,2sigma bounds on r'
        elseif ( line(1:6) == 'fix1' ) then
            fix2 = .FALSE. 
            if ( lout ) print '(a)','# starting months refer to field1/series'
        elseif ( line(1:6) == 'fix2' ) then
            fix2 = .TRUE. 
            if ( lout ) print '(a)','# starting months refer to field2/field/index'
        elseif ( line(1:3) == 'lag' ) then
            call get_command_argument(i+1,line)
            call deletebackslash(line)
            j = index(line,':')
            call numbersonly(line)
            if ( j /= 0 ) then
                read(line(:j-1),*,err=905) lag1
                read(line(j+1:),*,err=905) lag2
                if ( lout ) print '(a,2i6)','# using lags ',lag1,lag2
            else
                read(line,*,err=905) lag1
                if ( lout ) print '(a,i6)','# using lag ',lag1
                lag2 = lag1
            endif
            lskip = 1
        elseif ( line(1:3) == 'day' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=928) day0
            if ( day0 < 0 .OR. day0 > 31 ) goto 928
            if ( lout ) print '(a,i2)','# only using day ',day0
            lskip = 1
        elseif ( line(1:3) == 'mon' ) then
            call get_command_argument(i+1,line)
            call deletebackslash(line)
            j = index(line,':')
            if ( j /= 0 ) then
                read(line(:j-1),*,err=906) m1
                read(line(j+1:),*,err=906) m2
            else
                read(line,*,err=906) m1
                m2 = m1
            endif
            if ( m1 < 0 .OR. m1 > min(nperyear,12) ) goto 912
            if ( m2 < 0 .OR. m2 > min(nperyear,12) ) goto 912
            if ( m1 == 0 ) then
                if ( lout ) print '(a)','# using all months of the year'
            endif
            if ( m1 > 0 .OR. m1 /= m2 ) then
                if ( lout ) print '(a,i2,a,i2)','# only using starting months ',m1,':',m2
            endif
            lskip = 1
        elseif ( line(1:6) == 'begin2' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=907) yr1a
            yr1a = max(yr1a,yrbeg)
            if ( lout ) print '(a,i4)','# starting 2 at year ',yr1a
            lskip = 1
        elseif ( line(1:4) == 'end2' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=908) yr2a
            yr2a = min(yr2a,yrend)
            if ( lout ) print '(a,i4)','# ending 2 at year ',yr2a
            lskip = 1
        elseif ( line(1:4) == 'end3' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=908) yr2b
            if ( lout ) print '(a,i4)','# ending 3 at year ',yr2b
            lskip = 1
        elseif ( line(1:3) == 'beg' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=907) yr1
            yr1 = max(yr1,yrbeg)
            if ( yr1a == yrbeg ) yr1a = yr1
            if ( lout ) print '(a,i4)','# starting at year ',yr1
            lskip = 1
        elseif ( line(1:3) == 'end' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=908) yr2
            yr2 = min(yr2,yrend)
            if ( yr2a == yrend ) yr2a = yr2
            if ( lout ) print '(a,i4)','# ending at year ',yr2
            lskip = 1
        elseif ( line(1:9) == 'startstop' ) then
            lstartstop = .TRUE. 
            inquire(unit=12,opened=lopen)
            if ( .NOT. lopen ) then ! maybe getopts is called twice...
                call get_command_argument(i+1,line)
                inquire(file=trim(line),exist=lexist)
                l = len_trim(line)
                if ( lexist ) then
                    do j=0,9999
                        write(line(l+1:),'(i4.4)') j
                        inquire(file=trim(line),exist=lexist)
                        if ( .NOT. lexist ) exit
                    end do
                    ! move the old one out of the way
                    call mysystem('mv '//line(1:l)//' '//trim(line),iret)
                end if
                open(12,file=trim(line(1:l)),status='new')
                if ( lout ) print '(2a)','# writing start and stop years on file ',trim(line)
            end if
            lskip = 1
        elseif ( line(1:2) == 'lt' ) then
            call get_command_argument(i+1,line)
            j = index(line,'%')
            if ( line == 'n' ) then
                pmaxindx = 19712000
            else if ( j == 0 ) then
                read(line,*,err=909) maxindx
                if ( lout ) print '(a,g20.4)','# using maximum index ',maxindx
            else
                read(line(1:j-1),*,err=909) pmaxindx
                if ( lout ) print '(a,f8.2,a)','# using maximum index',pmaxindx,'%'
            endif
            lskip = 1
        elseif ( line(1:2) == 'gt' ) then
            call get_command_argument(i+1,line)
            j = index(line,'%')
            if ( line == 'n' ) then
                pminindx=19712000
            else if ( j == 0 ) then
                read(line,*,err=910) minindx
                if ( lout ) print '(a,g20.4)','# using minimum index ',minindx
            else
                read(line(1:j-1),*,err=910) pminindx
                if ( lout ) print '(a,f8.2,a)','# using minimum index',pminindx,'%'
            endif
            lskip = 1
        elseif ( line(1:3) == 'dlt' ) then
            call get_command_argument(i+1,line)
            j = index(line,'%')
            if ( j == 0 ) then
                read(line,*,err=909) maxdata
                if ( lout ) print '(a,g20.4)','# using maximum data ',maxdata
            else
                read(line(1:j-1),*,err=909) pmaxdata
                if ( lout ) print '(a,f8.2,a)','# using maximum data ',pmaxdata,'%'
            endif
            lskip = 1
        elseif ( line(1:3) == 'dgt' ) then
            call get_command_argument(i+1,line)
            j = index(line,'%')
            if ( j == 0 ) then
                read(line,*,err=910) mindata
                if ( lout ) print '(a,g20.4)','# using minimum data ',mindata
            else
                read(line(1:j-1),*,err=910) pmindata
                if ( lout ) print '(a,f8.2,a)','# using minimum data ',pmindata,'%'
            endif
            lskip = 1
        elseif ( line(1:5) == 'xyear' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=908) xyear
            if ( lout ) print '(a,g20.4)','# computing return value for ',xyear
            yr2a = 9999
            lskip = 1
        elseif ( line(1:3) == 'sum' ) then
            oper = '+'
            if ( line(4:4) == ' ' .OR. line(4:4) == 'm' .OR. line(4:4) == '1' ) then
                call get_command_argument(i+1,line)
                read(line,*,err=906) lsum
                if ( lsum <= 0 ) then
                    write(0,*) 'error: expecting sum>0,not ',lsum
                    call exit(-1)
                endif
                if ( lout ) print '(a,i4,a)','# summing ',lsum,' months/periods'
            elseif ( line(4:4) == '2' ) then
                call get_command_argument(i+1,line)
                read(line,*,err=906) lsum2
                if ( lsum <= 0 ) then
                    write(0,*) 'error: expecting sum2>0,not ',lsum2
                    call exit(-1)
                endif
                if ( lout ) print '(a,i4,a)','# summing ',lsum2,' months/periods of index/field/field2'
            else
                goto 906
            endif
            lskip = 1
        elseif ( line(1:3) == 'ave' ) then
            oper = 'v'
            if ( line(4:4) == ' ' .OR. line(4:4) == '1' ) then
                call get_command_argument(i+1,line)
                call numbersonly(line)
                read(line,*,err=906) lsum
                if ( lout ) print '(a,i4,a)','# averaging ',lsum,' months/periods'
            elseif ( line(4:4) == '2' ) then
                call get_command_argument(i+1,line)
                call numbersonly(line)
                read(line,*,err=906) lsum2
                if ( lout ) print '(a,i4,a)','# averaging ',lsum2,' months/periods of index/field/field2'
            else
                goto 906
            endif
            lskip = 1
        elseif ( line(1:3) == 'max' ) then
            oper = 'a'
            if ( line(4:4) == ' ' .OR. line(4:4) == '1' ) then
                call get_command_argument(i+1,line)
                call numbersonly(line)
                read(line,*,err=906) lsum
                if ( lout ) print '(a,i4,a)','# Max-ing ',lsum,' months/periods'
            elseif ( line(4:4) == '2' ) then
                call get_command_argument(i+1,line)
                call numbersonly(line)
                read(line,*,err=906) lsum2
                if ( lout ) print '(a,i4,a)','# Max-ing ',lsum2,' months/periods of index/field/field2'
            else
                goto 906
            endif
            lskip = 1
        elseif (  line(1:4) == 'min ' .OR. &
            line(1:4) == 'min1' .OR. &
            line(1:4) == 'min2' ) then
            oper = 'i'
            if ( line(4:4) == ' ' .OR. line(4:4) == '1' ) then
                call get_command_argument(i+1,line)
                call numbersonly(line)
                read(line,*,err=906) lsum
                if ( lout ) print '(a,i4,a)','Min-ing ',lsum,' months/periods'
            else
                call get_command_argument(i+1,line)
                call numbersonly(line)
                read(line,*,err=906) lsum2
                if ( lout ) print '(a,i4,a)','Min-ing ',lsum2,' months/periods  of index/field/field2'
            endif
            lskip = 1
        elseif ( line(1:5) == 'mdiff' ) then
            if ( line(6:6) == ' ' .OR. line(6:6) == '1' ) then
                call get_command_argument(i+1,line)
                read(line,*,err=906) mdiff
                if ( lout ) print '(a,i4,a)' &
                ,'# Taking anomalies wrt ',mdiff &
                ,' previous months/periods'
            elseif ( line(4:4) == '2' ) then
                call get_command_argument(i+1,line)
                read(line,*,err=906) mdiff2
                if ( lout ) print '(a,i4,a)','# Taking anomalies wrt ',mdiff2, &
                    ' previous months/periods of index/field/field2'
            else
                goto 906
            endif
            lskip = 1
        elseif ( line(1:4) == 'runc' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=923) nyrwindow
            if ( nyrwindow <= 0 ) goto 923
            if ( lout ) print '(a,i4,a)' &
            ,'# doing a running correlation analysis of ',nyrwindow,' years'
            irunvar = 1
            call get_command_argument(i+2,plotrunfile)
            if ( plotrunfile == ' ' ) goto 924
            open(14,file=trim(plotrunfile))
            lskip = 2
        elseif ( line(1:4) == 'runr' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=923) nyrwindow
            if ( nyrwindow <= 0 ) goto 923
            if ( lout ) print '(a,i4,a)','# doing a running regression analysis of ' , &
                nyrwindow,' years'
            irunvar = 2
            call get_command_argument(i+2,plotrunfile)
            if ( plotrunfile == ' ' ) goto 924
            open(14,file=trim(plotrunfile))
            lskip = 2
        elseif ( line(1:9) == 'runstartr' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=923) nyrwindow ! not used
            nyrwindow = 1
            if ( lout ) print '(a,i4,a)','# doing a regression with running start point'
            irunvar = -2
            call get_command_argument(i+2,plotrunfile)
            if ( plotrunfile == ' ' ) goto 924
            open(14,file=trim(plotrunfile))
            lskip = 2
        elseif ( line(1:8) == 'runstart' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=923) nyrwindow
            nyrwindow = 1
            if ( lout ) print '(a,i4,a)','# doing a correlation with running start point'
            irunvar = -1
            call get_command_argument(i+2,plotrunfile)
            if ( plotrunfile == ' ' ) goto 924
            open(14,file=trim(plotrunfile))
            lskip = 2
        elseif ( line(1:3) == 'run' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=923) nyrwindow
            if ( nyrwindow <= 0 ) goto 923
            if ( lout ) print '(a,i4,a)','# doing a running analysis of ' &
                ,nyrwindow,' years'
            irunvar = 0
            call get_command_argument(i+2,plotrunfile)
            if ( plotrunfile == ' ' ) goto 924
            open(14,file=trim(plotrunfile))
            lskip = 2
        elseif ( line(1:7) == 'obsplot' ) then
            call get_command_argument(i+1,plotrunfile)
            if ( plotrunfile == ' ' ) goto 924
            open(15,file=trim(plotrunfile))
            lskip = 1
        elseif ( line(1:6) == 'random' ) then
            call get_command_argument(i+1,line)
            if ( line(1:6) == 'series' ) then
                lrandom = .FALSE. 
                if ( lout ) print '(3a)','# significance running', &
                    ' correlations assesed by substituting Gauss for series'
            elseif ( line(1:5) == 'index' .OR. line(1:5) == 'field' ) then
                lrandom = .TRUE. 
                if ( lout ) print '(4a)','# significance running', &
                    ' correlations assesed by substiting Gauss for ',line(1:5)
            else
                print *,'getopts: error: unknown option for random'
                print *,trim(line)
                call exit(-1)
            endif
            lskip = 1
        elseif ( line(1:5) == 'noise' ) then
            call get_command_argument(i+1,line)
            if ( line(1:5) == 'white' ) then
                noisetype = 0
                if ( lout ) print '(3a)','# assuming white noise'
            elseif ( line(1:3) == 'red' ) then
                noisetype = 1
                if ( lout ) print '(3a)','# assuming red noise'
            else
                print *,'getopts: error: unknown option for noise'
                print *,trim(line)
                call exit(-1)
            endif
            lskip = 1
        elseif ( line(1:4) == 'conf' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=926) confidenceinterval
            if ( confidenceinterval < 50 .OR. &
            confidenceinterval >= 100 ) then
                write(0,*) 'getopts: error: confidence interval ', &
                    confidenceinterval,'% does not make sense, ', &
                    'use a number between 50 and 100'
                call exit(-1)
            end if
            if ( lout ) print '(a,g12.4,a)','# confidence interval ' &
                ,confidenceinterval,'%'
            lskip = 1
        elseif ( line(1:3) == 'sel' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=914) lsel
            if ( lsel < 1 .OR. lsel > nperyear ) goto 914
            if ( lout ) print '(a,i4,a)','# selecting ',lsel,' months'
            lskip = 1
        elseif ( line(1:9) == 'minfacsum' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=915) minfacsum
            if ( minfacsum >= 0.995 ) minfacsum = minfacsum/100
            if ( minfacsum < 0 .OR. minfacsum > 1 ) goto 915
            if ( lout ) print '(a,f4.2,a)','# requiring at least ' &
                ,minfacsum,' fraction valid points in sums'
            lskip = 1
        elseif ( line(1:6) == 'minfac' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=915) minfac
            if ( minfac >= 0.995 ) minfac = minfac/100
            if ( minfac < 0 .OR. minfac > 1 ) goto 915
            if ( lout ) print '(a,f4.2,a)','# requiring at least ' &
                ,minfac,' fraction valid points'
            lskip = 1
        elseif ( line(1:6) == 'minnum' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=920) minnum
            if ( lout ) print '(a,i6,a)','# requiring at least ' &
                ,minnum,' valid points'
            lskip = 1
        elseif ( line(1:4) == 'lon1' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=916) lon1
            if ( lout ) print '(a,f7.2)','# starting near longitude ',lon1
            lskip = 1
        elseif ( line(1:4) == 'lon2' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=916) lon2
            if ( lout ) print '(a,f7.2)' &
            ,'# stopping near longitude ',lon2
            lskip = 1
        elseif ( line(1:4) == 'lat1' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=916) lat1
            if ( lat1 < -90 .OR. lat1 > 90 ) goto 916
            if ( lout ) print '(a,f7.2)','# starting near latitude '
            lskip = 1
        elseif ( line(1:4) == 'lat2' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=916) lat2
            if ( lat2 < -90 .OR. lat2 > 90 ) goto 916
            if ( lout ) print '(a,f7.2)','# stopping near latitude ',lat2
            lskip = 1
        elseif ( line(1:4) == 'lev1' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=916) lev1
            if ( lout ) print '(a,f7.2)','# starting near level ',lev1
            lskip = 1
        elseif ( line(1:4) == 'lev2' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=916) lev2
            if ( lout ) print '(a,f7.2)','# stopping near level ',lev2
            lskip = 1
        elseif ( line(1:4) == 'xave' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=917) avex
            if ( avex < 1 ) goto 917
            if ( lout ) print '(a,i5,a)','# averaging over ',avex,' x grid points'
            lskip = 1
        elseif ( line(1:4) == 'yave' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=918) avey
            if ( avey < 1 ) goto 918
            if ( lout ) print '(a,i5,a)','# averaging over ',avey,' y grid points'
            lskip = 1
        elseif ( line(1:7) == 'altlon1' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=916) altlon1
            if ( lout ) print '(a,f7.2)','# field 2 starting near longitude ',altlon1
            lskip = 1
        elseif ( line(1:7) == 'altlon2' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=916) altlon2
            if ( lout ) print '(a,f7.2)','# field 2 stopping near longitude ',altlon2
            lskip = 1
        elseif ( line(1:7) == 'altlat1' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=916) altlat1
            if ( altlat1 < -90 .OR. altlat1 > 90 ) goto 916
            if ( lout ) print '(a,f7.2)' &
                ,'# field 2 starting near latitude ',altlat1
            lskip = 1
        elseif ( line(1:7) == 'altlat2' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=916) altlat2
            if ( altlat2 < -90 .OR. altlat2 > 90 ) goto 916
            if ( lout ) print '(a,f7.2)' &
                ,'# field 2 stopping near latitude ',altlat2
            lskip = 1
        elseif ( line(1:7) == 'altlev1' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=916) altlev1
            if ( lout ) print '(a,f7.2)' &
                ,'# field 2 starting near level ',altlev1
            lskip = 1
        elseif ( line(1:7) == 'altlev2' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=916) altlev2
            if ( lout ) print '(a,f7.2)' &
                ,'# field 2 stopping near level ',altlev2
            lskip = 1
        elseif ( line(1:7) == 'altxave' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=917) altavex
            if ( altavex < 1 ) goto 917
            if ( lout ) print '(a,i5,a)','# averaging field 2 over ' &
                ,altavex,' x grid points'
            lskip = 1
        elseif ( line(1:7) == 'altyave' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=918) altavey
            if ( altavey < 1 ) goto 918
            if ( lout ) print '(a,i5,a)','# averaging field 2 over ' &
                ,altavey,' y grid points'
            lskip = 1
        elseif ( line(1:7) == 'ensanom' ) then
            lensanom = .TRUE. 
            if ( lout ) print '(a)' &
            ,'# taking anomalies wrt ensemble'
        elseif ( line(1:3) == 'ens' ) then
            call get_command_argument(i+2,line)
            if ( line /= ' ' ) then
                read(line,*,err=922) nens2
                call get_command_argument(i+1,line)
                read(line,*,err=922) nens1
                if ( lout ) print '(a,i3,a,i3)','# requested ensemble members ',nens1,' to ',nens2
                lskip = 2
            end if
        elseif ( line(1:11) == 'makeensfull' ) then
            lmakeensfull = .TRUE. 
            if ( lout ) print '(a)','# augmenting ensembles with duplicates to constant size'
        elseif ( line(1:13) == 'standardunits' ) then
            lstandardunits = .TRUE. 
            if ( lout ) print '(a)' &
            ,'# converting to standard units'
        elseif ( line(1:6) == 'allobs' ) then
            lallobs = .TRUE. 
            if ( lout ) print '(a)','# using all observations to define the climatology'
        elseif ( line(1:4) == 'lead' ) then
            leads = .TRUE. 
            call get_command_argument(i+1,line)
            read(line,*,err=925) lead1
            call get_command_argument(i+2,line)
            read(line,*,err=925) lead2
            if ( lead1 /= 0 .OR. lead2 /= 0 ) then
                if ( lout ) print '(a,i3,a,i3)','# using lead times ',lead1,' to ',lead2
            endif
            lskip = 2
        elseif ( line(1:8) == 'restrain' ) then
            call get_command_argument(i+1,line)
            read(line,*,err=925) restrain
            if ( restrain /= 0 ) then
                if ( lout ) print '(a,f4.2)' ,'# restraining shape parameter to ',restrain
            endif
            lskip = 1
        elseif ( line(1:5) == 'inter' ) then
            call get_command_argument(i+1,line)
            if ( line(1:1) == '1' ) then
                intertype = 1
            elseif ( line(1:1) == '2' ) then
                intertype = 2
            elseif ( line(1:2) == 'hi' ) then
                intertype = 0
            elseif ( line(1:2) == 'lo' ) then
                intertype = -1
            else
                goto 921
            endif
            if ( lout ) print '(2a)','# interpolating to grid ',trim(line)
            lskip = 1
        elseif ( line(1:4) == 'dump' ) then
            if ( i < iarg2 ) then
                call get_command_argument(i+1,plotfile)
                lskip = 1
            else
                plotfile = 'dump.dat'
            endif
            if ( .NOT. dump ) then
                dump = .TRUE. 
                if ( lwrite ) print '(2a)','# writing dumpfile on ',trim(plotfile)
                open(10,file=trim(plotfile))
                call get_command_argument(0,line)
                k = index(line,'/', .TRUE. )
                line = line(k+1:)
                do j=1,command_argument_count()
                    call get_command_argument(j,string)
                    k = index(string,'/', .TRUE. )
                    line = trim(line)//' '//string(k+1:)
                end do
                write(10,'(2a)') '# ',trim(line)
            end if
        elseif ( line(1:4) == 'plot' ) then
            if ( .NOT. plot ) then
                plot = .TRUE. 
                call get_command_argument(i+1,plotfile)
                if ( lwrite ) print *,'Writing plotfile on ',trim(plotfile)
                open(11,file=trim(plotfile))
            end if
            lskip = 1
        elseif ( line(1:6) == 'lwrite' .OR. line(1:5) == 'debug' ) then
            lwrite = .TRUE. 
            if ( lout ) print *,'turned on debug printing'
        elseif ( line(1:4) == 'diff' ) then
            if ( i <= iarg2 ) then
                call get_command_argument(i+1,line)
                if ( line(1:1) == '-' .OR. &
                    ichar(line(1:1)) >= ichar('0') .AND. &
                    ichar(line(1:1)) <= ichar('9') ) then
                    read(line,*) j
                    lskip = 1
                else
                    j = 1
                endif
            else
                j = 1
            endif
            if ( j > 0 ) then
                if ( lout ) print '(a,i3,a)','# taking anomaly to ',j,' previous years'
            elseif ( j < 0 ) then
                if ( lout ) print '(a,i3,a)','# summing with ',-j,' previous years'
            endif
            if ( j /= 0 ) then
                idiff = idiff + 1
                if ( idiff == 1 ) then
                    ndiff = j
                elseif ( idiff == 2 ) then
                    ndiff2 = j
                else
                    write(0,*) 'getopts: error: more than 2 diffs: ',idiff
                endif
            endif
        elseif ( line(1:6) == 'noover' ) then
            if ( lout ) print '(a)','# using non-overlapping intervals'
            lnooverlap = .TRUE. 
        elseif ( line(1:6) == 'crossv' ) then
            lskip = 2
            call get_command_argument(i+1,line)
            read(line,*) ncrossvalidate
            if ( lout ) print '(a,i4,a)' &
            ,'# computing cross-validated correlations and fits leaving out ', &
                ncrossvalidate,' time steps'
            call get_command_argument(i+2,bbfile)
            if ( bbfile == 'none' ) then
                bbfile = ' '
            end if
            if ( lout .AND. bbfile /= ' ' ) then
                print '(3a)','# writing cross-validated regression in file ',trim(bbfile)
            end if
        elseif ( line(1:4) == 'detr' ) then
            ldetrend = .TRUE. 
            if ( lout ) print '(a)','# detrending both fields'
        elseif ( line(1:6) == 'debias' ) then
            lskip = 1
            call get_command_argument(i+1,line)
            if ( line(1:4) == 'none' ) then
                debias = 0
                if ( lout ) print '(a)','# no bias correction'
            elseif ( line(1:4) == 'mean' ) then
                debias = 1
                if ( lout ) print '(a)' &
                ,'# removing bias in mean (with jackknife)'
            elseif ( line(1:3) == 'var' ) then
                debias = 2
                if ( lout ) print '(a)','# removing bias in mean and variance (with jackknife)'
            elseif ( line(1:3) == 'all' ) then
                debias = 3
                if ( lout ) print '(a)','# removing bias in whole PDF'
            else
                lskip = 0
            endif
        elseif ( line(1:7) == 'biasmul' ) then
            lskip = 1
            call get_command_argument(i+1,line)
            j = index(line,'%')
            if ( j == 0 ) then ! factor
                read(line,*,err=929) biasmul
            else ! percentage
                read(line(:j-1),*) biasmul
                biasmul = 1 + biasmul/100 ! convert to factor
            end if
            if ( lout ) print '(a)','# used multiplicative bias correction ',biasmul
        elseif ( line(1:7) == 'biasadd' ) then
            lskip = 1
            call get_command_argument(i+1,line)
            read(line,*,err=929) biasadd
            if ( lout ) print '(a)','# used additive bias correction ',biasadd
        elseif ( line(1:7) == 'biasrt' ) then
            lskip = 1
            call get_command_argument(i+1,line)
            read(line,*,err=929) biasrt
            if ( lout ) print '(a,f10.1)','# evaluate model for return time ',biasrt
        elseif ( line(1:11) == 'includelast' ) then
            lincludelast = .true.
            if ( lout ) print '(a)','# in,cude event itself in fit'
        elseif ( line(1:6) == 'normal' ) then
            lskip = 1
            call get_command_argument(i+1,line)
            if ( line(1:8) == 'maxspace' ) then
                normalization = 1
                if ( lout ) print '(a)','# normalize spatial pattern to maximum 1'
            elseif ( line(1:8) == 'varspace' ) then
                normalization = 2
                if ( lout ) print '(a)','# normalize spatial pattern to variance 1'
            elseif ( line(1:7) == 'maxtime' ) then
                normalization = -1
                if ( lout ) print '(a)','# normalize time series maximum to 1 '
            elseif ( line(1:7) == 'vartime' ) then
                normalization = -2
                if ( lout ) print '(a)','# normalize time series variance to 1'
            else
                write(0,*) 'getopts: unrecognized suboption: ',trim(line)
                lskip = 0
            endif
        elseif ( line(1:4) == 'anom' ) then
            anom = .TRUE. 
            if ( lout ) print '(a)','# subtracting seasonal cycle'
        elseif ( line(1:4) == 'comp' ) then
            composite = .TRUE. 
            if ( lout ) print '(a)','# making composites'
        elseif ( line(1:3) == 'sub' ) then
            lsubtract = .TRUE. 
            if ( lout ) print '(a)' &
            ,'subtracting best fit to produce new field'
        elseif ( line(1:4) == 'name' ) then
            call get_command_argument(i+1,namestring)
            do j=1,len(namestring)
                if ( namestring(j:j) == '_' ) namestring(j:j) = ' '
            enddo
            lskip = 1
            if ( lout ) print '(2a)','# name = ',namestring
        elseif ( line(1:3) == 'soi' ) then
            lincl(1) = .TRUE. 
        elseif ( line(1:6) == 'nino12' ) then
            lincl(2) = .TRUE. 
        elseif ( line(1:6) == 'nino3 ' ) then
            lincl(3) = .TRUE. 
        elseif ( line(1:6) == 'nino4 ' ) then
            lincl(4) = .TRUE. 
        elseif ( line(1:7) == 'nino3.4' .OR. line(1:7) == 'nino34' ) then
            lincl(5) = .TRUE. 
        elseif ( line(1:3) == 'nao' ) then
            lincl(6) = .TRUE. 
        elseif ( line(1:7) == 'sunspot' ) then
            lincl(7) = .TRUE. 
        elseif ( line(1:9) == 'sunlength' ) then
            lincl(8) = .TRUE. 
        elseif ( line(1:4) == 'time' ) then
            lincl(9) = .TRUE. 
        elseif ( line(1:4) == 'file' ) then
            indxuse = indxuse + 1
            lincl(indxuse) = .TRUE. 
            call get_command_argument(i+1,indexfiles(indxuse))
            do j=len(indexfiles(indxuse)),1,-1
                if ( indexfiles(indxuse)(j:j) == '/' ) goto 10
            enddo
        10  continue
            strindx(indxuse) = indexfiles(indxuse)(j+1:j+len(strindx(indxuse)))
            lskip = 1
        elseif ( line /= ' ' ) then
            if ( lout ) print '(a)','getopts: warning: unrecognized option'
            if ( lout ) print '(a)',trim(line)
        endif
        50 continue
    enddo
    if ( lconting ) then
        if ( minindx < -1e30 .AND. pminindx < 0 ) then
            pminindx = 100/3.
            if ( lout ) print *,'assumed lower cut on index to be ',pminindx
        endif
        if ( maxindx > +1e30 .AND. pmaxindx < 0 ) then
            pmaxindx = 200/3.
            if ( lout ) print *,'assumed upper cut on index to be ',pmaxindx
        endif
        if ( mindata < -1e30 .AND. pmindata < 0 ) then
            pmindata = 100/3.
            if ( lout ) print *,'assumed lower cut on data to be ',pmindata
        endif
        if ( maxdata > +1e30 .AND. pmaxdata < 0 ) then
            pmaxdata = 200/3.
            if ( lout ) print *,'assumed upper cut on data to be ' ,pmaxdata
        endif
    endif
    if ( lsum2 < 0 ) lsum2 = lsum
    if ( mdiff2 < 0 ) mdiff2 = mdiff
    if ( minfac < 0 ) minfac = 0.5 ! the original algorithm does not work well
    if ( minfacsum < 0 ) minfacsum = 0.999
    if ( lrank .and. ncrossvalidate > 0 ) then
        write(0,*) 'getopts: error: cannot use cross-validation with rank correlations yet'
        call exit(-1)
    end if
    return

!   error messages

905 print *,'getopts: error reading lag value ',trim(line)
    call exit(-1)
906 print *,'getopts: error reading sum value ',trim(line)
    call exit(-1)
907 print *,'getopts: error reading begin year value ',trim(line)
    call exit(-1)
908 print *,'getopts: error reading end year value ',trim(line)
    call exit(-1)
909 print *,'getopts: error reading maximum data value ',trim(line)
    call exit(-1)
910 print *,'getopts: error reading minimum data value ',trim(line)
    call exit(-1)
911 print *,'getopts: error reading minimum ks cut-off or "plot"',trim(line)
    call exit(-1)
912 print *,'getopts: error reading month value ',trim(line), &
        ', should be integer between 1 and ',nperyear
    call exit(-1)
913 print *,'getopts: error reading ldecorrelation length ',trim(line)
    call exit(-1)
914 print *,'getopts: error reading selection value ',trim(line),nperyear
    call exit(-1)
915 print *,'getopts: error reading minfac value ',trim(line)
    call exit(-1)
916 print *,'getopts: error reading lon,lat value ',trim(line)
    call exit(-1)
917 print *,'getopts: error: invalid value for xave: ',trim(line)
    call exit(-1)
918 print *,'getopts: error: invalid value for yave: ',trim(line)
    call exit(-1)
919 print *,'getopts: error: invalid value for fittime: ',trim(line)
    call exit(-1)
920 print *,'getopts: error reading minnum value ',trim(line)
    call exit(-1)
921 print *,'getopts: error reading interpolation type ',trim(line)
    call exit(-1)
922 print *,'getopts: error reading ensemble numbers from ',trim(line)
    call exit(-1)
923 print *,'getopts: error: expecting a window length of more than 0: ',trim(line)
    call exit(-1)
924 print *,'getopts: error: expecting a file name for running correlations'
    call exit(-1)
925 print *,'getopts: error reading leads from ',trim(line)
    call exit(-1)
926 print *,'getopts: error reading factor from ',trim(line)
    call exit(-1)
927 print *,'getopts: error reading degree of fit polynomial from ',trim(line)
    call exit(-1)
928 print *,'getopts: error: invalid value for day: ',trim(line)
    call exit(-1)
929 print *,'getopts: error: cannot read bias correction from: ',trim(line)
    call exit(-1)
end subroutine getopts

subroutine month2period(m,nperyear,offset)

!   offset = 1: beginning of month
!   offset = 0: end of month

    implicit none
    integer :: m,nperyear,offset,y
    integer :: j,k,dpm(12)
    logical :: lwrite
    parameter(lwrite= .FALSE. )
    data dpm /31,29,31,30,31,30,31,31,30,31,30,31/

!   normalize m
    if ( m < 1 ) then
        print *,'month2period: error: m<1: ',m
        call exit(-1)
    endif
    if ( lwrite ) print *,'y,m,nperyear=',y,m,nperyear
    y = (m-1)/12
    m = 1 + mod(m-1,12)
    if ( lwrite ) print *,'y,m         =',y,m

!   a few easy cases

    if ( nperyear <= 12 ) then
        return
    elseif ( nperyear == 36 ) then
        if ( offset == 0 ) then
            m = 3*m
        else
            m = 3*m-2
        endif
        return
    elseif ( nperyear == 360 ) then
        if ( offset == 0 ) then
            m = 30*m
        else
            m = 30*m-29
        endif
        return
    endif

!   compute day of year corresponding to the first day of m
    j=offset
    do k=1,m-offset
        j = j+dpm(k)
    enddo
    if ( lwrite ) print *,'k,offset,j = ',k,offset,j
!   and the corresponding period
    if ( nperyear == 366 ) then
        m = j
    elseif ( nperyear == 365 ) then
        if ( m <= 2 ) then
            m = j
        else
            m = j-1
        endif
    else
        if ( m <= 2 ) then
            m = nint(j/(365./nperyear))
        else
            m = nint((j-1)/(365./nperyear))
        endif
        m = m + offset
    endif
    if ( lwrite ) print *,'m = ',m
    m = m + nperyear*y
    if ( lwrite ) print *,'m = ',m
end subroutine month2period

subroutine deletebackslash(line)
    implicit none
    character*(*) line
    integer :: j
    j = index(line,'\\')
    !!!print *,'@@@ j,line = ',j,trim(line)
    do while ( j /= 0 )
        line = line(:j-1) // line(j+1:)
        j = index(line,'\\')
    !!!print *,'@@@ j,line = ',j,trim(line)
    end do
end subroutine

subroutine numbersonly(line)
!
!   delete everything that is not part of a number
!
    implicit none
    character :: line*(*)
    integer i,ic
    character :: c

    do i=1,len_trim(line)
        c = line(i:i)
        ic = ichar(c)
        if ( ic >= ichar('0') .and. ic <= ichar('9') ) cycle
        if ( c == '.' ) cycle
        if ( c == '+' .or. c== '-' ) cycle
        if ( c == 'e' .or. c == 'E' ) cycle ! accept scientific notation
        if ( c == 'd' .or. c == 'D' ) cycle ! accept old Fortran double precision scientific notation
        line(i:i) = ' ' ! everything else is thrown away
    end do
end subroutine