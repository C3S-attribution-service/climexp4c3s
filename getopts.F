        subroutine getopts(iarg1,iarg2,nperyear,yrbeg,yrend,loutin,mens1
     +       ,mens)
*
*       process options common to correlat*.[fF]
*       the switches are in common in getopts.inc
*
        implicit none
        integer iarg1,iarg2,nperyear,yrbeg,yrend,mens1,mens
        logical loutin
#include "getopts.inc"
        integer i,j,k,l,lskip,idiff
        real s,minfacinit
        logical lout,lexist
        character line*511,string*511
        integer iargc,llen
*
*       init
*
        lout = loutin
        logscale = .FALSE.
        sqrtscale = .FALSE.
        lchangesign = .FALSE.
        lstandardunits = .false.
        lallobs = .false.
        logfield = .FALSE.
        sqrtfield = .FALSE.
        squarescale = .FALSE.
        lnormsd = .FALSE.
        lnomissing = .FALSE.
        lrank = .FALSE.
        lfitnoise = .FALSE.
        ldetrend = .FALSE.
        debias = 0
        normalization = 1
	lks = .FALSE.
	lconting = .FALSE.
        lbootstrap = .FALSE.
        lrmse = .FALSE.
        irunvar = 0
        noisetype = 0
        anom = .FALSE.
        lensanom = .FALSE.
        composite = .FALSE.
        lrandom = .FALSE.
        idiff = 0
        ndiff = 0
        ndiff2 = 0
        lnooverlap = .false.
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
        if ( nperyear.eq.1 ) then
            m2 = 0
        else
            m2 = min(12,nperyear)
        endif
        yr1 = yrbeg
        yr2 = yrend
        yr1a = yrbeg
        yr2a = yrend + 1
        lstartstop = .FALSE.
        minindx = -2e33
        maxindx = +2e33
        mindata = -2e33
        maxdata = +2e33
        pminindx = -1
        pmaxindx = -1
        pmindata = -1
        pmaxdata = -1
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
        factor = 1
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
        if ( line.ne.' ' ) then
            lweb = .true.
            lout = .false.
        endif
*
*       loop over arguments
*
	lskip = 0
        do i=iarg1,iarg2
            if ( lskip.gt.0 ) then
                lskip = lskip - 1
                goto 50
            endif
            call getarg(i,line)
            if ( lwrite ) print *,'getopts: parsing ',i,trim(line)
            if ( line(1:8).eq.'logfield' ) then
                logfield = .TRUE.
                if ( lout ) print '(a)','# set logscale on for field'
            elseif ( line(1:3).eq.'log' ) then
                logscale = .TRUE.
                if ( lout ) print '(a)',
     +               '# set logscale on for timeseries'
            elseif ( line(1:9).eq.'sqrtfield' ) then
                sqrtfield = .TRUE.
                if ( lout ) print '(a)','# set sqrt on for field'
            elseif ( line(1:4).eq.'sqrt' ) then
                sqrtscale = .TRUE.
                if ( lout ) print '(a)',
     +               '# set sqrt on for time series'
            elseif ( line(1:6).eq.'square' ) then
                squarescale = .TRUE.
                if ( lout ) print '(a)','# set square on'
            elseif ( line(1:10).eq.'changesign' ) then
                lchangesign = .TRUE.
                if ( lout ) print '(a)','# change sign of data'
            elseif ( line(1:6).eq.'normsd' ) then
                lnormsd = .TRUE.
                if ( lout ) print '(a)'
     +                ,'# normalize to standard deviation'
            elseif ( line(1:6).eq.'nomiss' ) then
                lnomissing = .TRUE.
                if ( lout ) print '(a)','# no missing data'
            elseif ( line(1:4).eq.'rank' ) then
                lrank = .TRUE.
                if ( lout ) print '(a)','# set rank correllations on'
            elseif ( line(1:8).eq.'fitnoise' ) then
                lfitnoise = .TRUE.
                if ( lout ) print '(a)','# fit to noise model'
            elseif ( line(1:7).eq.'fitfunc' ) then
                call getarg(i+1,line)
                if ( line(1:3).eq.'lin' ) then
                    fitfunc = 1
                elseif ( line(1:4).eq.'quad' .or. line(1:5).eq.'parab' )
     +               then
                    fitfunc = 2
                elseif ( line(1:3).eq.'cub' ) then
                    fitfunc = 3
                else
                    read(line,*,err=927) fitfunc
                endif
                if ( lout ) print '(a,i2)'
     +               ,'# fit to polynomial of degree ',fitfunc
                lskip = 1
            elseif ( line(1:7).eq.'fittime' ) then
                if ( lout ) print '(a)'
     +                ,'# fitting time derivative with relaxation'
                call getarg(i+1,line)
                if (  ichar(line(1:1)).ge.ichar('0') .and.
     +                ichar(line(1:1)).le.ichar('9') ) then
                    read(line,*,err=919) nfittime
                    if ( nfittime.le.0 ) goto 919
                    if ( lout ) print '(a,i3,a)'
     +                    ,'# computing time derivative with ',nfittime,
     +                    ' points'
                    lskip = 1
                else
                    nfittime = 1
                endif
            elseif ( line(1:3).eq.'dec' ) then
                call getarg(i+1,line)
                read(line,*,err=913) decor
                decor = abs(decor)
                if ( lout ) print '(a,f6.2)'
     +                ,'# computing significances with decorrelation '
     +                ,decor
                lskip = 1
            elseif ( line(1:2).eq.'ks' ) then
                lks = .TRUE.
                call getarg(i+1,line)
                if ( line(1:4).eq.'plot' ) then
                    kscut = -3e33
                    if ( lout ) print *
     +                    ,'Performing KS test, plot in dump.dat'
                    dump = .TRUE.
                    open(10,file='dump.dat')
                else
                    read(line,*,err=911) kscut
                    if ( lout ) print *,'Performing KS test at cut-off '
     +                    ,kscut
                endif
                lskip = 1
            elseif ( line(1:6).eq.'contin' ) then
                lconting = .TRUE.
                if ( lout ) print '(2a)'
     +                ,'# computing 3x3 contingency table with cuts '
     +                ,'min/max index/data'
            elseif ( line(1:3).eq.'rms' ) then
                lrmse = .TRUE.
                if ( lout ) print '(2a)'
     +                ,'# computing rms error, bias, ratio variance'
            elseif ( line(1:4).eq.'boot' ) then
                lbootstrap = .TRUE.
                if ( lout ) print '(a)'
     +                ,'# Bootstrapping to get 1,2sigma bounds on r'
            elseif ( line(1:6).eq.'fix1' ) then
                fix2 = .FALSE.
                if ( lout ) print '(a)'
     +                ,'# starting months refer to field1/series'
            elseif ( line(1:6).eq.'fix2' ) then
                fix2 = .TRUE.
                if ( lout ) print '(a)'
     +                ,'# starting months refer to field2/field/index'
            elseif ( line(1:3).eq.'lag' ) then
                call getarg(i+1,line)
	        j = index(line,':')
	        if ( j.ne.0 ) then
	            read(line(:j-1),*,err=905) lag1
	            read(line(j+1:),*,err=905) lag2
                    if ( lout ) print '(a,2i6)','# using lags ',lag1
     +                    ,lag2
                else
                    read(line,*,err=905) lag1
                    if ( lout ) print '(a,i6)','# using lag ',lag1
                    lag2 = lag1
                endif
                lskip = 1
            elseif ( line(1:3).eq.'day' ) then
                call getarg(i+1,line)
                read(line,*,err=928) day0
                if ( day0.lt.0 .or. day0.gt.31 ) goto 928
                if ( lout ) print '(a,i2)','# only using day ',day0
                lskip = 1
            elseif ( line(1:3).eq.'mon' ) then
                call getarg(i+1,line)
	        j = index(line,':')
	        if ( j.ne.0 ) then
	            read(line(:j-1),*,err=906) m1
	            read(line(j+1:),*,err=906) m2
	        else
                    read(line,*,err=906) m1
                    m2 = m1
	        endif
                if ( m1.lt.0 .or. m1.gt.min(nperyear,12) ) goto 912
                if ( m2.lt.0 .or. m2.gt.min(nperyear,12) ) goto 912
                if ( m1.eq.0 ) then
                    if ( lout ) print '(a)'
     +                    ,'# using all months of the year'
                endif
                if ( m1.gt.0 .or. m1.ne.m2 ) then
                    if ( lout ) print '(a,i2,a,i2)'
     +                    ,'# only using starting months ',m1,':',m2
                endif
                lskip = 1
            elseif ( line(1:6).eq.'begin2' ) then
                call getarg(i+1,line)
                read(line,*,err=907) yr1a
                yr1a = max(yr1a,yrbeg)
                if ( lout ) print '(a,i4)','# starting 2 at year ',yr1a
                lskip = 1
            elseif ( line(1:4).eq.'end2' ) then
                call getarg(i+1,line)
                read(line,*,err=908) yr2a
                yr2a = min(yr2a,yrend)
                if ( lout ) print '(a,i4)','# ending 2 at year ',yr2a
                lskip = 1
            elseif ( line(1:3).eq.'beg' ) then
                call getarg(i+1,line)
                read(line,*,err=907) yr1
                yr1 = max(yr1,yrbeg)
                if ( yr1a.eq.yrbeg ) yr1a = yr1
                if ( lout ) print '(a,i4)','# starting at year ',yr1
                lskip = 1
            elseif ( line(1:3).eq.'end' ) then
                call getarg(i+1,line)
                read(line,*,err=908) yr2
                yr2 = min(yr2,yrend)
                if ( yr2a.eq.yrend ) yr2a = yr2
                if ( lout ) print '(a,i4)','# ending at year ',yr2
                lskip = 1
            elseif ( line(1:9).eq.'startstop' ) then
                call getarg(i+1,line)
                inquire(file=trim(line),exist=lexist)
                if ( lexist ) then
                    l = len_trim(line)
                    do j=0,9999
                        write(line(l:),'(i4.4)') j
                        inquire(file=trim(line),exist=lexist)
                        if ( .not.lexist ) exit
                    end do
                end if
                open(12,file=trim(line),status='new')
                if ( lout ) print '(2a)'
     +               ,'# writing start and stop years on file '
     +               ,trim(line)
                lskip = 1
            elseif ( line(1:2).eq.'lt' ) then
                call getarg(i+1,line)
                j = index(line,'%')
                if ( line.eq.'n' ) then
                    pmaxindx = 19712000
                else if ( j.eq.0 ) then
                    read(line,*,err=909) maxindx
                    if ( lout ) print '(a,g20.4)'
     +                    ,'# using maximum index ',maxindx
                else
                    read(line(1:j-1),*,err=909) pmaxindx
                    if ( lout ) print '(a,f8.2,a)'
     +                    ,'# using maximum index',pmaxindx,'%'
                endif
                lskip = 1
            elseif ( line(1:2).eq.'gt' ) then
                call getarg(i+1,line)
                j = index(line,'%')
                if ( line.eq.'n' ) then
                    pminindx=19712000
                else if ( j.eq.0 ) then
                    read(line,*,err=910) minindx
                    if ( lout ) print '(a,g20.4)'
     +                    ,'# using minimum index ',minindx
                else
                    read(line(1:j-1),*,err=910) pminindx
                    if ( lout ) print '(a,f8.2,a)'
     +                    ,'# using minimum index',pminindx,'%'
                endif
                lskip = 1
            elseif ( line(1:3).eq.'dlt' ) then
                call getarg(i+1,line)
                j = index(line,'%')
                if ( j.eq.0 ) then
                    read(line,*,err=909) maxdata
                    if ( lout ) print '(a,g20.4)'
     +                    ,'# using maximum data ',maxdata
                else
                    read(line(1:j-1),*,err=909) pmaxdata
                    if ( lout ) print '(a,f8.2,a)'
     +                    ,'# using maximum data ',pmaxdata,'%'
                endif
                lskip = 1
            elseif ( line(1:3).eq.'dgt' ) then
                call getarg(i+1,line)
                j = index(line,'%')
                if ( j.eq.0 ) then
                    read(line,*,err=910) mindata
                    if ( lout ) print '(a,g20.4)'
     +                    ,'# using minimum data ',mindata
                else
                    read(line(1:j-1),*,err=910) pmindata
                    if ( lout ) print '(a,f8.2,a)'
     +                    ,'# using minimum data ',pmindata,'%'
                endif
                lskip = 1
            elseif ( line(1:3).eq.'sum' ) then
                oper = '+'
                if ( line(4:4).eq.' ' .or. line(4:4).eq.'m' .or.
     +               line(4:4).eq.'1' ) then
                    call getarg(i+1,line)
                    read(line,*,err=906) lsum
                    if ( lsum.le.0 ) then
                        write(0,*) 'error: expecting sum>0,not ',lsum
                        call abort
                    endif
                    if ( lout ) print '(a,i4,a)','# summing ',lsum
     +                    ,' months/periods'
                elseif ( line(4:4).eq.'2' ) then
                    call getarg(i+1,line)
                    read(line,*,err=906) lsum2
                    if ( lsum.le.0 ) then
                        write(0,*) 'error: expecting sum2>0,not ',lsum2
                        call abort
                    endif
                    if ( lout ) print '(a,i4,a)','# summing ',lsum2
     +                    ,' months/periods of index/field/field2'
                else
                    goto 906
                endif
                lskip = 1
            elseif ( line(1:3).eq.'ave' ) then
                oper = 'v'
                if ( line(4:4).eq.' ' .or. line(4:4).eq.'1' ) then
                    call getarg(i+1,line)
                    read(line,*,err=906) lsum
                    if ( lout ) print '(a,i4,a)','# averaging ',lsum
     +                    ,' months/periods'
                elseif ( line(4:4).eq.'2' ) then
                    call getarg(i+1,line)
                    read(line,*,err=906) lsum2
                    if ( lout ) print '(a,i4,a)','# averaging ',lsum2
     +                    ,' months/periods of index/field/field2'
                else
                    goto 906
                endif
                lskip = 1
            elseif ( line(1:3).eq.'max' ) then
                oper = 'a'
                if ( line(4:4).eq.' ' .or. line(4:4).eq.'1' ) then
                    call getarg(i+1,line)
                    read(line,*,err=906) lsum
                    if ( lout ) print '(a,i4,a)','# Max-ing ',lsum
     +                    ,' months/periods'
                elseif ( line(4:4).eq.'2' ) then
                    call getarg(i+1,line)
                    read(line,*,err=906) lsum2
                    if ( lout ) print '(a,i4,a)','# Max-ing ',lsum2
     +                    ,' months/periods of index/field/field2'
                else
                    goto 906
                endif
                lskip = 1
            elseif (  line(1:4).eq.'min ' .or. 
     +                line(1:4).eq.'min1' .or. 
     +                line(1:4).eq.'min2' ) then
                oper = 'i'
                if ( line(4:4).eq.' ' .or. line(4:4).eq.'1' ) then
                    call getarg(i+1,line)
                    read(line,*,err=906) lsum
                    if ( lout ) print '(a,i4,a)','Min-ing ',lsum
     +                    ,' months/periods'
                else
                    call getarg(i+1,line)
                    read(line,*,err=906) lsum2
                    if ( lout ) print '(a,i4,a)','Min-ing ',lsum2
     +                    ,' months/periods  of index/field/field2'
                endif
                lskip = 1
            elseif ( line(1:5).eq.'mdiff' ) then
                if ( line(6:6).eq.' ' .or. line(6:6).eq.'1' ) then
                    call getarg(i+1,line)
                    read(line,*,err=906) mdiff
                    if ( lout ) print '(a,i4,a)'
     +                    ,'# Taking anomalies wrt ',mdiff
     +                    ,' previous months/periods'
                elseif ( line(4:4).eq.'2' ) then
                    call getarg(i+1,line)
                    read(line,*,err=906) mdiff2
                    if ( lout ) print '(a,i4,a)'
     +                    ,'# Taking anomalies wrt ',mdiff2,' previous'
     +                    //' months/periods of index/field/field2'
                else
                    goto 906
                endif
                lskip = 1
            elseif ( line(1:4).eq.'runc' ) then
                call getarg(i+1,line)
                read(line,*,err=923) nyrwindow
                if ( nyrwindow.le.0 ) goto 923
                if ( lout ) print '(a,i4,a)'
     +               ,'# doing a running correlation analysis of '
     +               ,nyrwindow,' years'
                irunvar = 1
                call getarg(i+2,plotrunfile)
                if ( plotrunfile.eq.' ' ) goto 924
                open(14,file=trim(plotrunfile))
                lskip = 2
            elseif ( line(1:4).eq.'runr' ) then
                call getarg(i+1,line)
                read(line,*,err=923) nyrwindow
                if ( nyrwindow.le.0 ) goto 923
                if ( lout ) print '(a,i4,a)'
     +               ,'# doing a running regression analysis of '
     +               ,nyrwindow,' years'
                irunvar = 2
                call getarg(i+2,plotrunfile)
                if ( plotrunfile.eq.' ' ) goto 924
                open(14,file=trim(plotrunfile))
                lskip = 2
            elseif ( line(1:3).eq.'run' ) then
                call getarg(i+1,line)
                read(line,*,err=923) nyrwindow
                if ( nyrwindow.le.0 ) goto 923
                if ( lout ) print '(a,i4,a)'
     +               ,'# doing a running analysis of '
     +               ,nyrwindow,' years'
                irunvar = 0
                call getarg(i+2,plotrunfile)
                if ( plotrunfile.eq.' ' ) goto 924
                open(14,file=trim(plotrunfile))
                lskip = 2
            elseif ( line(1:6).eq.'random' ) then
                call getarg(i+1,line)
                if ( line(1:6).eq.'series' ) then
                    lrandom = .false.
                    if ( lout ) print '(3a)','# significance running',
     +                   ' correlations assesed by substituting',
     +                   ' Gauss for series'
                elseif ( line(1:5).eq.'index' .or. line(1:5).eq.'field'
     +                   ) then
                    lrandom = .true.
                    if ( lout ) print '(4a)','# significance running',
     +                   ' correlations assesed by substiting',
     +                   ' Gauss for ',line(1:5)
                else
                    print *,'getopts: error: unknown option for random'
                    print *,line(1:llen(line))
                    call abort
                endif
                lskip = 1
            elseif ( line(1:5).eq.'noise' ) then
                call getarg(i+1,line)
                if ( line(1:5).eq.'white' ) then
                    noisetype = 0
                    if ( lout ) print '(3a)','# assuming white noise'
                elseif ( line(1:3).eq.'red' ) then
                    noisetype = 1
                    if ( lout ) print '(3a)','# assuming red noise'
                else
                    print *,'getopts: error: unknown option for noise'
                    print *,line(1:llen(line))
                    call abort
                endif
                lskip = 1
            elseif ( line(1:4).eq.'fact' ) then
                call getarg(i+1,line)
                read(line,*,err=926) factor
                if ( lout ) print '(a,g12.4)','# multiplying field1 by '
     +               ,factor
                lskip = 1
            elseif ( line(1:3).eq.'sel' ) then
                call getarg(i+1,line)
                read(line,*,err=914) lsel
                if ( lsel.lt.1 .or. lsel.gt.nperyear ) goto 914
                if ( lout ) print '(a,i4,a)','# selecting ',lsel
     +                ,' months'
                lskip = 1
            elseif ( line(1:9).eq.'minfacsum' ) then
                call getarg(i+1,line)
                read(line,*,err=915) minfacsum
                if ( minfacsum.ge.0.995 ) minfacsum = minfacsum/100
		if ( minfacsum.lt.0 .or. minfacsum.gt.1 ) goto 915
                if ( lout ) print '(a,f4.2,a)','# requiring at least '
     +               ,minfacsum,' fraction valid points in sums'
                lskip = 1
            elseif ( line(1:6).eq.'minfac' ) then
                call getarg(i+1,line)
                read(line,*,err=915) minfac
                if ( minfac.ge.0.995 ) minfac = minfac/100
		if ( minfac.lt.0 .or. minfac.gt.1 ) goto 915
                if ( lout ) print '(a,f4.2,a)','# requiring at least '
     +                ,minfac,' fraction valid points'
                lskip = 1
            elseif ( line(1:6).eq.'minnum' ) then
                call getarg(i+1,line)
                read(line,*,err=920) minnum
                if ( lout ) print '(a,i6,a)','# requiring at least '
     +                ,minnum,' valid points'
                lskip = 1
            elseif ( line(1:4).eq.'lon1' ) then
                call getarg(i+1,line)
                read(line,*,err=916) lon1
                if ( lout ) print '(a,f7.2)'
     +                ,'# starting near longitude ',lon1
                lskip = 1
            elseif ( line(1:4).eq.'lon2' ) then
                call getarg(i+1,line)
                read(line,*,err=916) lon2
                if ( lout ) print '(a,f7.2)'
     +                ,'# stopping near longitude ',lon2
                lskip = 1
            elseif ( line(1:4).eq.'lat1' ) then
                call getarg(i+1,line)
                read(line,*,err=916) lat1
                if ( lat1.lt.-90 .or. lat1.gt.90 ) goto 916
                if ( lout ) print '(a,f7.2)','# starting near latitude '
     +                ,lat1
                lskip = 1
            elseif ( line(1:4).eq.'lat2' ) then
                call getarg(i+1,line)
                read(line,*,err=916) lat2
                if ( lat2.lt.-90 .or. lat2.gt.90 ) goto 916
                if ( lout ) print '(a,f7.2)','# stopping near latitude '
     +                ,lat2
                lskip = 1
            elseif ( line(1:4).eq.'lev1' ) then
                call getarg(i+1,line)
                read(line,*,err=916) lev1
                if ( lout ) print '(a,f7.2)','# starting near level '
     +                ,lev1
                lskip = 1
            elseif ( line(1:4).eq.'lev2' ) then
                call getarg(i+1,line)
                read(line,*,err=916) lev2
                if ( lout ) print '(a,f7.2)','# stopping near level '
     +                ,lev2
                lskip = 1
            elseif ( line(1:4).eq.'xave' ) then
                call getarg(i+1,line)
                read(line,*,err=917) avex
                if ( avex.lt.1 ) goto 917
                if ( lout ) print '(a,i5,a)','# averaging over ',avex,
     +			' x grid points'
                lskip = 1
            elseif ( line(1:4).eq.'yave' ) then
                call getarg(i+1,line)
                read(line,*,err=918) avey
                if ( avey.lt.1 ) goto 918
                if ( lout ) print '(a,i5,a)','# averaging over ',avey,
     +			' y grid points'
                lskip = 1
            elseif ( line(1:7).eq.'altlon1' ) then
                call getarg(i+1,line)
                read(line,*,err=916) altlon1
                if ( lout ) print '(a,f7.2)'
     +                ,'# field 2 starting near longitude ',altlon1
                lskip = 1
            elseif ( line(1:7).eq.'altlon2' ) then
                call getarg(i+1,line)
                read(line,*,err=916) altlon2
                if ( lout ) print '(a,f7.2)'
     +                ,'# field 2 stopping near longitude ',altlon2
                lskip = 1
            elseif ( line(1:7).eq.'altlat1' ) then
                call getarg(i+1,line)
                read(line,*,err=916) altlat1
                if ( altlat1.lt.-90 .or. altlat1.gt.90 ) goto 916
                if ( lout ) print '(a,f7.2)'
     +               ,'# field 2 starting near latitude ',altlat1
                lskip = 1
            elseif ( line(1:7).eq.'altlat2' ) then
                call getarg(i+1,line)
                read(line,*,err=916) altlat2
                if ( altlat2.lt.-90 .or. altlat2.gt.90 ) goto 916
                if ( lout ) print '(a,f7.2)'
     +               ,'# field 2 stopping near latitude ',altlat2
                lskip = 1
            elseif ( line(1:7).eq.'altlev1' ) then
                call getarg(i+1,line)
                read(line,*,err=916) altlev1
                if ( lout ) print '(a,f7.2)'
     +               ,'# field 2 starting near level ',altlev1
                lskip = 1
            elseif ( line(1:7).eq.'altlev2' ) then
                call getarg(i+1,line)
                read(line,*,err=916) altlev2
                if ( lout ) print '(a,f7.2)'
     +               ,'# field 2 stopping near level ',altlev2
                lskip = 1
            elseif ( line(1:7).eq.'altxave' ) then
                call getarg(i+1,line)
                read(line,*,err=917) altavex
                if ( altavex.lt.1 ) goto 917
                if ( lout ) print '(a,i5,a)','# averaging field 2 over '
     +               ,altavex,' x grid points'
                lskip = 1
            elseif ( line(1:7).eq.'altyave' ) then
                call getarg(i+1,line)
                read(line,*,err=918) altavey
                if ( altavey.lt.1 ) goto 918
                if ( lout ) print '(a,i5,a)','# averaging field 2 over '
     +               ,altavey,' y grid points'
                lskip = 1
            elseif ( line(1:7).eq.'ensanom' ) then
                lensanom = .TRUE.
                if ( lout ) print '(a)'
     +               ,'# taking anomalies wrt ensemble'
            elseif ( line(1:3).eq.'ens' ) then
                call getarg(i+2,line)
                if ( line.ne.' ' ) then
                    read(line,*,err=922) nens2
                    call getarg(i+1,line)
                    read(line,*,err=922) nens1
                    if ( lout ) print '(a,i3,a,i3)'
     +                   ,'# requested ensemble members ',nens1,' to '
     +                   ,nens2
                    lskip = 2
                end if
            elseif ( line(1:11).eq.'makeensfull' ) then
                lmakeensfull = .TRUE.
                if ( lout ) print '(a)'
     +                ,'# augmenting ensembles with duplicates to '//
     +               'constant size'
            elseif ( line(1:13).eq.'standardunits' ) then
                lstandardunits = .TRUE.
                if ( lout ) print '(a)'
     +                ,'# converting to standard units'
            elseif ( line(1:6).eq.'allobs' ) then
                lallobs = .TRUE.
                if ( lout ) print '(a)','# using all observations '//
     +               'to define the climatology'
            elseif ( line(1:4).eq.'lead' ) then
                leads = .TRUE.
                call getarg(i+1,line)
                read(line,*,err=925) lead1
                call getarg(i+2,line)
                read(line,*,err=925) lead2
                if ( lead1.ne.0 .or. lead2.ne.0 ) then
                    if ( lout ) print '(a,i3,a,i3)'
     +                    ,'# using lead times ',lead1,' to ',lead2
                endif
                lskip = 2
            elseif ( line(1:8).eq.'restrain' ) then
                call getarg(i+1,line)
                read(line,*,err=925) restrain
                if ( restrain.ne.0 ) then
                    if ( lout ) print '(a,f4.2)'
     +                   ,'# restraining shape parameter to ',restrain
                endif
                lskip = 1
            elseif ( line(1:5).eq.'inter' ) then
                call getarg(i+1,line)
                if ( line(1:1).eq.'1' ) then
                    intertype = 1
                elseif ( line(1:1).eq.'2' ) then
                    intertype = 2
                elseif ( line(1:2).eq.'hi' ) then
                    intertype = 0
                elseif ( line(1:2).eq.'lo' ) then
                    intertype = -1
                else
                    goto 921
                endif
                if ( lout ) print '(2a)','# interpolating to grid ',
     +                line(1:llen(line))
                lskip = 1
            elseif ( line(1:4).eq.'dump' ) then
                dump = .TRUE.
                if ( i.lt.iarg2 ) then
                    call getarg(i+1,plotfile)
                    lskip = 1
                else
                    plotfile = 'dump.dat'
                endif
                if ( lwrite ) print '(2a)','# writing dumpfile on ',
     +                trim(plotfile)
                open(10,file=trim(plotfile))
                call getarg(0,line)
                k = index(line,'/',.true.)
                line = line(k+1:)
                do j=1,iargc()
                    call getarg(j,string)
                    k = index(string,'/',.true.)
                    line = trim(line)//' '//string(k+1:)
                end do
                write(10,'(2a)') '# ',trim(line)
            elseif ( line(1:4).eq.'plot' ) then
                plot = .TRUE.
                call getarg(i+1,plotfile)
                lskip = 1
                if ( lwrite ) print *,'Writing plotfile on ',
     +                plotfile(1:index(plotfile,' ')-1)
                open(11,file=trim(plotfile))
            elseif ( line(1:6).eq.'lwrite' .or. line(1:5).eq.'debug' )
     +                then
                lwrite = .TRUE.
                if ( lout ) print *,'turned on debug printing'
            elseif ( line(1:4).eq.'diff' ) then
                if ( i.le.iarg2 ) then
                    call getarg(i+1,line)
                    if ( line(1:1).eq.'-' .or.
     +                    ichar(line(1:1)).ge.ichar('0') .and.
     +                    ichar(line(1:1)).le.ichar('9') ) then
                        read(line,*) j
                        lskip = 1
                    else
                        j = 1
                    endif
                else
                    j = 1
                endif
                if ( j.gt.0 ) then
                    if ( lout ) print '(a,i3,a)','# taking anomaly to '
     +                    ,j,' previous years'
                elseif ( j.lt.0 ) then
                    if ( lout ) print '(a,i3,a)','# summing with ',-j
     +                    ,' previous years'
                endif
                if ( j.ne.0 ) then
                    idiff = idiff + 1
                    if ( idiff.eq.1 ) then
                        ndiff = j
                    elseif ( idiff.eq.2 ) then
                        ndiff2 = j
                    else
                        write(0,*) 'getopts: error: more than '//
     +                            '2 diffs: ',idiff
                    endif
                endif
            elseif ( line(1:6).eq.'noover' ) then
                if ( lout ) print '(a)','# using non-overlapping '//
     +               'intervals'
                lnooverlap = .true.
            elseif ( line(1:6).eq.'crossv' ) then
                lskip = 2
                call getarg(i+1,line)
                read(line,*) ncrossvalidate
                if ( lout ) print '(a,i4,a)'
     +               ,'# computing cross-validated correlations '//
     +               'and fits leaving out ',ncrossvalidate,
     +               ' time steps'
                call getarg(i+2,bbfile)
                if ( bbfile.eq.'none' ) then
                    bbfile = ' '
                end if
                if ( lout .and. bbfile.ne.' ' ) then
                    print '(3a)','# writing cross-validated ',
     +                   'regression in file ',trim(bbfile)
                end if
            elseif ( line(1:4).eq.'detr' ) then
                ldetrend = .TRUE.
                if ( lout ) print '(a)','# detrending both fields'
            elseif ( line(1:6).eq.'debias' ) then
                lskip = 1
                call getarg(i+1,line)
                if ( line(1:4).eq.'none' ) then
                    debias = 0
                    if ( lout ) print '(a)','# no bias correction'
                elseif ( line(1:4).eq.'mean' ) then
                    debias = 1
                    if ( lout ) print '(a)'
     +                  ,'# removing bias in mean (with jackknife)'
                elseif ( line(1:3).eq.'var' ) then
                    debias = 2
                    if ( lout ) print '(a)'
     +               ,'# removing bias in mean and variance '//
     +                '(with jackknife)'
                elseif ( line(1:3).eq.'all' ) then
                    debias = 3
                    if ( lout ) print '(a)'
     +               ,'# removing bias in whole PDF'
                else
                    lskip = 0
                endif
            elseif ( line(1:6).eq.'normal' ) then
                lskip = 1
                call getarg(i+1,line)
                if ( line(1:8).eq.'maxspace' ) then
                    normalization = 1
                    if ( lout ) print '(a)'
     +                   ,'# normalize spatial pattern to maximum 1'
                elseif ( line(1:8).eq.'varspace' ) then
                    normalization = 2
                    if ( lout ) print '(a)'
     +                   ,'# normalize spatial pattern to variance 1'
                elseif ( line(1:7).eq.'maxtime' ) then
                    normalization = -1
                    if ( lout ) print '(a)'
     +                   ,'# normalize time series maximum to 1 '
                elseif ( line(1:7).eq.'vartime' ) then
                    normalization = -2
                    if ( lout ) print '(a)'
     +                   ,'# normalize time series variance to 1'
                else
                    write(0,*) 'getopts: unrecognized suboption: '
     +                   ,trim(line)
                    lskip = 0
                endif
            elseif ( line(1:4).eq.'anom' ) then
                anom = .TRUE.
                if ( lout ) print '(a)','# subtracting seasonal cycle'
            elseif ( line(1:4).eq.'comp' ) then
                composite = .TRUE.
                if ( lout ) print '(a)','# making composites'
            elseif ( line(1:3).eq.'sub' ) then
                lsubtract = .TRUE.
                if ( lout ) print '(a)'
     +                ,'subtracting best fit to produce new field'
            elseif ( line(1:4).eq.'name' ) then
                call getarg(i+1,namestring)
                do j=1,len(namestring)
                    if ( namestring(j:j).eq.'_' ) namestring(j:j) = ' '
                enddo
                lskip = 1
                if ( lout ) print '(2a)','# name = ',namestring
            elseif ( line(1:3).eq.'soi' ) then
                lincl(1) = .TRUE.
            elseif ( line(1:6).eq.'nino12' ) then
                lincl(2) = .TRUE.
            elseif ( line(1:6).eq.'nino3 ' ) then
                lincl(3) = .TRUE.
            elseif ( line(1:6).eq.'nino4 ' ) then
                lincl(4) = .TRUE.
            elseif ( line(1:7).eq.'nino3.4' .or. line(1:7).eq.'nino34'
     +                ) then
                lincl(5) = .TRUE.
            elseif ( line(1:3).eq.'nao' ) then
                lincl(6) = .TRUE.
            elseif ( line(1:7).eq.'sunspot' ) then
                lincl(7) = .TRUE.
            elseif ( line(1:9).eq.'sunlength' ) then
                lincl(8) = .TRUE.
            elseif ( line(1:4).eq.'time' ) then
                lincl(9) = .TRUE.
            elseif ( line(1:4).eq.'file' ) then
                indxuse = indxuse + 1
                lincl(indxuse) = .TRUE.
                call getarg(i+1,indexfiles(indxuse))
                do j=len(indexfiles(indxuse)),1,-1
                    if ( indexfiles(indxuse)(j:j).eq.'/' ) goto 10
                enddo
   10           continue
                strindx(indxuse) = indexfiles(indxuse)
     +               (j+1:j+len(strindx(indxuse)))
                lskip = 1
            elseif ( line.ne.' ' ) then
                if ( lout ) print '(a)'
     +                ,'getopts: warning: unrecognized option'
                if ( lout ) print '(a)',trim(line)
            endif
   50       continue
        enddo
        if ( lconting ) then
            if ( minindx.lt.-1e30 .and. pminindx.lt.0 ) then
                pminindx = 100/3.
                if ( lout ) print *,'assumed lower cut on index to be '
     +                ,pminindx
            endif
            if ( maxindx.gt.+1e30 .and. pmaxindx.lt.0 ) then
                pmaxindx = 200/3.
                if ( lout ) print *,'assumed upper cut on index to be '
     +                ,pmaxindx
            endif
            if ( mindata.lt.-1e30 .and. pmindata.lt.0 ) then
                pmindata = 100/3.
                if ( lout ) print *,'assumed lower cut on data to be '
     +                ,pmindata
            endif
            if ( maxdata.gt.+1e30 .and. pmaxdata.lt.0 ) then
                pmaxdata = 200/3.
                if ( lout ) print *,'assumed upper cut on data to be '
     +                ,pmaxdata
            endif
        endif
        if ( lsum2.lt.0 ) lsum2 = lsum
        if ( mdiff2.lt.0 ) mdiff2 = mdiff
        if ( minfac.lt.0 ) minfac = 0.5 ! the original algorithm does not work well
        if ( minfacsum.lt.0 ) minfacsum = 0.76
        if ( lrank .and. ncrossvalidate.gt.0 ) then
            write(0,*) 'getopts: error: cannot use cross-validation ',
     +           'with rank correlations yet'
            call abort
        end if
        return
*
*       error messages
*
 905    print *,'getopts: error reading lag value ',trim(line)
        call abort
 906    print *,'getopts: error reading sum value ',trim(line)
        call abort
 907    print *,'getopts: error reading begin year value ',trim(line)
        call abort
 908    print *,'getopts: error reading end year value ',trim(line)
        call abort
 909    print *,'getopts: error reading maximum data value ',trim(line)
        call abort
 910    print *,'getopts: error reading minimum data value ',trim(line)
        call abort
 911    print *,'getopts: error reading minimum ks cut-off or "plot" '
     +       ,trim(line)
        call abort
 912    print *,'getopts: error reading month value ',trim(line),
     +       ', should be integer between 1 and ',nperyear
        call abort
 913    print *,'getopts: error reading ldecorrelation length '
     +       ,trim(line)
        call abort
 914    print *,'getopts: error reading selection value ',trim(line)
     +       ,nperyear
        call abort
 915    print *,'getopts: error reading minfac value ',trim(line)
        call abort
 916    print *,'getopts: error reading lon,lat value ',trim(line)
        call abort
 917    print *,'getopts: error: invalid value for xave: ',trim(line)
        call abort
 918    print *,'getopts: error: invalid value for yave: ',trim(line)
        call abort
 919    print *,'getopts: error: invalid value for fittime: ',trim(line)
        call abort
 920    print *,'getopts: error reading minnum value ',trim(line)
        call abort
 921    print *,'getopts: error reading interpolation type ',trim(line)
        call abort
 922    print *,'getopts: error reading ensemble numbers from '
     +       ,trim(line)
        call abort
 923    print *,'getopts: error: expecting a window length of more '//
     +       'than 0: ',trim(line)
        call abort
 924    print *,'getopts: error: expecting a file name for running '//
     +       'correlations'
        call abort
 925    print *,'getopts: error reading leads from ',trim(line)
        call abort
 926    print *,'getopts: error reading factor from ',trim(line)
        call abort
 927    print *,'getopts: error reading degree of fit polynomial from '
     +       ,trim(line)
        call abort
 928    print *,'getopts: error: invalid value for day: ',trim(line)
        call abort
        end

        subroutine month2period(m,nperyear,offset)
*
*       offset = 1: beginning of month
*       offset = 0: end of month
*
        implicit none
        integer m,nperyear,offset,y
        integer j,k,dpm(12)
        logical lwrite
        parameter(lwrite=.FALSE.)
        data dpm /31,29,31,30,31,30,31,31,30,31,30,31/
*
*       normalize m
        if ( m.lt.1 ) then
            print *,'month2period: error: m<1: ',m
            call abort
        endif
        if ( lwrite ) print *,'y,m,nperyear=',y,m,nperyear
        y = (m-1)/12
        m = 1 + mod(m-1,12)
        if ( lwrite ) print *,'y,m         =',y,m
*
*       a few easy cases
*
        if ( nperyear.le.12 ) then
            return
        elseif ( nperyear.eq.36 ) then
            if ( offset.eq.0 ) then
                m = 3*m
            else
                m = 3*m-2
            endif
            return
        elseif ( nperyear.eq.360 ) then
            if ( offset.eq.0 ) then
                m = 30*m
            else
                m = 30*m-29
            endif
            return
        endif
*
*       compute day of year corresponding to the first day of m
        j=offset
        do k=1,m-offset
            j = j+dpm(k)
        enddo
        if ( lwrite ) print *,'k,offset,j = ',k,offset,j
*       and the corresponding period
        if ( nperyear.eq.366 ) then
            m = j
        elseif ( nperyear.eq.365 ) then
            if ( m.le.2 ) then
                m = j
            else
                m = j-1
            endif
        else
            if ( m.le.2 ) then
                m = nint(j/(365./nperyear))
            else
                m = nint((j-1)/(365./nperyear))
            endif
	    m = m + offset
        endif
        if ( lwrite ) print *,'m = ',m
        m = m + nperyear*y
        if ( lwrite ) print *,'m = ',m
        end
