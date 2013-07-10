        program lomb
*
*       compute the Lomb periodogram of a climate series
*
        implicit none
        include 'param.inc'
        integer i,j,j1,j2,n,n1,n2,mon,sum,nout,jmax,ndiff,nperyear,yr1
     +       ,yr2,ave,mens1,mens,iens,nepx(4*ndata),nn,nens1,nens2,yr,mo
        real data(npermax,yrbeg:yrend,0:nensmax),
     +        x(ndata),y(ndata),px(4*ndata),py(4*ndata),
     +        epx(4*ndata),epy(4*ndata),prob,hifac,ofac,fac,sx,sy,
     +        xmax,ymax,timescale
        character line*255,oper*1,var*40,units*20
        logical lwrite,ldetrend,anom
	integer iargc
*
*       init
*
        lwrite = .false.
        yr1 = yrbeg
        yr2 = yrend
        ave = 1
        if ( iargc().eq.0 ) then
            print *,'usage: lomb datafile [month n] [sum m] [detrend]'//
     +            ' [anom] [diff [n]] [begin yyyy] [end yyyy] [ave n]'
            stop
        endif
*
        n = 0
        mon = 0
        sum = 1
        ldetrend = .FALSE.
        anom = .FALSE.
	ndiff = 0
        nens1 = 0
        nens2 = 0
        call getarg(1,line)
        call readensseries(line,data,npermax,yrbeg,yrend,nensmax
     +       ,nperyear,mens1,mens,var,units,.false.,lwrite)
        if ( mens1.gt.0 ) nens1 = mens1
        if ( mens.gt.0 )  nens2 = mens
        do i=2,iargc()
            if ( n.gt.0 ) then
                n = n-1
            else
                call getarg(i,line)
                if ( line(1:3).eq.'mon' ) then
                    call getarg(i+1,line)
                    n = n + 1
                    read(line,*,err=900) mon
                    if ( mon.eq.0 ) then
                        write(*,'(a,i2)') '# using monthly data'
                    else
                        write(*,'(a,i2)')
     +                    '# only using yearly data starting in month '
     +                        ,mon
                    endif
                elseif ( line(1:3).eq.'sum' ) then
                    call getarg(i+1,line)
                    n = n + 1
                    read(line,*,err=900) sum
                    oper = '+'
                    if ( nperyear.eq.12 ) then
                        write(*,'(a,i2,a)') '# summed over ',sum
     +                        ,' months'
                    else
                        write(*,'(a,i2,a)') '# summed over ',sum
     +                        ,' periods'
                    endif
                elseif ( line(1:3).eq.'det' ) then
                    ldetrend = .TRUE.
                    write(*,'(a,i2,a)') '# detrended'
                elseif ( line(1:3).eq.'beg' ) then
                    call getarg(i+1,line)
                    n = n + 1
                    read(line,*,err=900) yr1
                    write(*,'(a,i4)') '# first year ',yr1
                elseif ( line(1:3).eq.'end' ) then
                    call getarg(i+1,line)
                    n = n + 1
                    read(line,*,err=900) yr2
                    write(*,'(a,i4)') '# last year ',yr2
                elseif ( line(1:3).eq.'ave' ) then
                    call getarg(i+1,line)
                    n = n + 1
                    read(line,*,err=900) ave
                    write(*,'(a,i4,a)') '# averaging over ',ave,' bins'
                elseif ( line(1:3).eq.'ano' ) then
                    anom = .TRUE.
                    write(*,'(a,i2,a)') '# cut out yearly cycle'
                elseif ( line(1:4).eq.'diff' ) then
                    if ( i.lt.iargc() ) then
                        call getarg(i+1,line)
                        if ( ichar(line(1:1)).ge.ichar('0') .and.
     +                        ichar(line(1:1)).le.ichar('9') ) then
                            read(line,*) ndiff
                            n = n + 1
                        else
                            ndiff = 1
                        endif
                    else
                        ndiff = 1
                    endif
                    print '(a,i3,a)','# taking anomaly to ',ndiff
     +                    ,' previous years'
                elseif ( line(1:3).eq.'ens' ) then
                    call getarg(i+1,line)
                    read(line,*,err=900) nens1
                    nens1 = max(nens1,mens1)
                    call getarg(i+2,line)
                    read(line,*,err=900) nens2
                    nens2 = min(nens2,mens)
                    if ( nens1.ne.0 .or. nens2.ne.0 ) then
                        if ( lwrite )print '(a,i3,a,i3)'
     +                        ,'# using ensemble members ',nens1,' to '
     +                        ,nens2
                    endif
                    n = n + 2
                elseif ( line(1:5).eq.'debug' ) then
                    lwrite = .true.
                    print *,'Turned on debug information'
                else
                    print '(2a)','lomb: error: unrecognized option '
     +                   ,trim(line)
                endif
            endif
        enddo
        if ( mens.gt.0 ) then
            write(0,*) 'The spectra are computed for ensemble members '
     +            ,nens1,' to ',nens2
            if ( nens2-nens1.gt.0 ) write(0,*) ' and then averaged'
            write(0,*) '<br>'
            print '(a,i4,a,i4)','# taking ensemble members ',nens1
     +            ,' to ',nens2
        endif
        if ( mon.eq.0 ) then
            j1 = 1
            j2 = nperyear
        else
            j1 = mon
            j2 = mon
	    if ( nperyear.ne.12 ) then
		call month2period(j1,nperyear,1)
		call month2period(j2,nperyear,0)
                if ( j2.lt.j1 ) j2 = j2 + nperyear
	    endif
        endif
        if ( lwrite ) then
            if ( nperyear.eq.12  ) then
                print *,'lomb: using monhs ',j1,' to ',j2
            else
                print *,'lomb: using periods ',j1,' to ',j2
            endif
        endif
*
*       sum
*
        if ( sum.gt.1 ) then
            if ( lwrite ) print *,'Summing data'
            do iens=nens1,nens2
                call sumit(data(1,yrbeg,iens),npermax,nperyear,yrbeg
     +                ,yrend,sum,oper)
            enddo
        endif        
*
*       detrending
*       
        if ( ldetrend ) then
            if ( lwrite ) print *,'Detrending data'
            do iens=nens1,nens2
                call detrend(data(1,yrbeg,iens),npermax,nperyear,
     +               yrbeg,yrend,yr1,yr2,mon,mon,1)
            enddo
        endif
        if ( ndiff.gt.0 ) then
            if ( lwrite ) print *,'Taking anomalies - data'
            do iens=nens1,nens2
                call diffit(data(1,yrbeg,iens),npermax,nperyear,yrbeg
     +                ,yrend,ndiff)
            enddo
        endif
*
*       anomalies w.r.t. seasonal cycle
*
        if ( anom ) then
            do iens=nens1,nens2
                call anomal(data(1,yrbeg,iens),npermax,nperyear,yrbeg
     +                ,yrend,yrbeg,yrend)
            enddo
        endif
*       
*       loop over ensemble members
*       
        do iens=nens1,nens2
*
*       fill arrays
*       
            n = 0
            n1 = 0
            do yr=yr1,yr2
                do mo=j1,j2
                    j = mo
                    call normon(j,yr,i,nperyear)
                    if ( abs(data(j,i,iens)).lt.1e33 ) then
                        n = n + 1
                        x(n) = i + (j-0.5)/nperyear
                        if ( nperyear.gt.12 ) then
                            x(n) = x(n)*365.24
                        endif
                        y(n) = data(j,i,iens)
                        if ( n1.eq.0 ) then
                            n1 = nperyear*i+j
                        endif
                        n2 = nperyear*i+j
                    endif
                enddo
            enddo
*       
*       compute/estimate other parameters
*       
            if ( mon.eq.0 ) then
                hifac = real(n)/real(n2-n1+1)
            else
                hifac = real(n)/real((n2-n1)/nperyear+1)
            endif
            if ( hifac.eq.1 ) then
                ofac = 1
            else
                ofac = 4        ! see how it works
            endif
*       
*       call period (Numerical recipes p 572)
*       take care of dependent data!
*       
            call period(x,y,n,ofac,hifac,px,py,4*ndata,nout,jmax,prob)
*       
*       average
*       
            if ( ave.gt.1 ) then
                do i=nint(ofac),nout-ave,ave
                    sx = px(i)
                    sy = py(i)
                    do j=1,ave-1
                        sx = sx + px(i+j)
                        sy = sy + py(i+j)
                    enddo
                    n = 1+i/ave
                    px(n) = sx/ave
                    py(n) = sy/ave
                enddo
                ofac = 1+nint(ofac)/ave
                nout = n
            endif
*       
*       collect ensemble informnation
*            
            if ( iens.eq.nens1 ) then
                nn = nout
                do i=1,nn
                    nepx(i) = 1
                    epx(i) = px(i)
                    epy(i) = py(i)
                    if ( lwrite ) print *,i,epx(i),epy(i),nepx(i)
                enddo
            else
                do i=1,nout
                    if ( epx(i).eq.px(i) ) then
                        nepx(i) = nepx(i) + 1
                        epy(i) = epy(i) + py(i)
                    else        ! unequal array sizes - choose nearest
                        do j=1,nn-1
                            if ( px(j).lt.(epx(j)+epx(j+1))/2 ) then
                                goto 100
                            endif
                        enddo
  100                   continue
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
*
*       output; do not yet know how to treat dependent data.
*
        print '(a)','# frequency power'
        do i=nint(ofac),nout
            if ( epx(i).gt.0 ) then
                write(*,'(2f12.6)') epx(i),epy(i)
            else
                write(0,*) 'lomb: something is wrong',i,epx(i)
                call abort
            endif
        enddo
*
*       the 95% significance line - correct for dependent data
*
        write(0,'(a)') '<font color="#ff0000">Changed normalisation '//
     $       ' to equal-area with log(T) x-axis</font></br>'
        fac = 1
        if ( mon.eq.0 ) then
            if ( sum.gt.1 ) fac = sum
        else
            if ( sum.gt.nperyear ) fac = sum/real(nperyear)
        endif
        call maxquad(xmax,ymax,epx(nint(ofac)),epy(nint(ofac)),
     +        nout-nint(ofac)+1)
        if ( fac.eq.1 .and. ave.eq.1 .and. mens.eq.0 ) then
            write(0,'(a)') 'The horizontal line denotes the 95% '//
     +            'significance limit.'
            write(0,'(a,f8.2,a,f6.2,a)') 'The peak at ',1/px(jmax)
     +            ,'yr has a significance of ',100*(1-prob),'%.'
            write(0,'(a)') 'These numbers refer to the chance of '//
     +            'any peak not exceeding this value in a '//
     +            'white-noise spectrum.<br>'
            write(*,'(a)')
            do i=1,nout
                sy = -log(0.05/(2*nout/ofac))*fac
                if ( epx(i)*sy.lt.xmax*ymax ) then
                    write(*,'(2f12.6)') epx(i),sy
                endif
            enddo
        endif
        if ( ave.gt.1 ) then
            write(0,'(a)') 'The spectrum has been obtained by averaging'
            write(0,'(i3,a)') ave,' bins of the periodogram<br>'
        endif
        write(0,'(a,2f12.6,a)') 'The highest peak is at ',1/xmax,
     +        xmax*ymax,' (in a log(T) plot)<br>'
*
        goto 999
  900   print *,'error: cannot read integer from ',line
        call abort
  999   continue
        end
