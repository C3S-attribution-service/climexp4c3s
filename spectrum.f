        program spectrum
*
*       compute a primitive spectrum of a climate series
*
        implicit none
        include 'param.inc'
        include 'getopts.inc'
        integer i,j,j1,j2,n,nout,nperyear,mens1,mens
     +       ,iens,nn,yr,mo,yrstart,yrstop,imc,nmc,nprob,nimax
        real data(npermax,yrbeg:yrend,0:nensmax),epx(4*ndata),
     +       epy(4*ndata),ofac,fac,xmax,ymax,alpha,s2,sd,probmax,mean,
     +       mean1
        real,allocatable :: epx1(:),epy1(:),epyall(:,:),prob(:)
        integer,allocatable :: imax(:)
        character line*255,var*40,units*20
	integer iargc
*
*       init
*
        lwrite = .false.
        if ( iargc().eq.0 ) then
            print *,'usage: spectrum datafile [month n] [sum m] '//
     +           '[detrend] [anom] [ensanom] [diff [n]] '//
     +           '[begin yyyy] [end yyyy] [xave n]'
            stop
        endif
*
        n = 0
        call getarg(1,line)
        call readensseries(line,data,npermax,yrbeg,yrend,nensmax
     +       ,nperyear,mens1,mens,var,units,.false.,lwrite)
        call getopts(2,iargc(),nperyear,yrbeg,yrend,.true.,mens1,mens)
!       I am abusing yr1a,yr2a here...
        if ( yr1a.eq.yr1 ) yr1a=yrbeg
        if ( yr2a.eq.yr2 ) yr2a = yrend + 1
        yrstart = yr2
        yrstop  = yr1
        if ( mens.gt.0 ) then
            write(0,*) 'The spectra are computed for ensemble members '
     +            ,nens1,' to ',nens2
            if ( nens2-nens1.gt.0 ) write(0,*) ' and then averaged'
            write(0,*) '<br>'
            print '(a,i4,a,i4)','# taking ensemble members ',nens1
     +            ,' to ',nens2
        endif
        if ( m1.eq.0 ) then
            j1 = 1
            j2 = nperyear
        else
            j1 = m1
            j2 = m1
	    if ( nperyear.ne.12 ) then
		call month2period(j1,nperyear,1)
		call month2period(j2,nperyear,0)
                if ( j2.lt.j1 ) j2 = j2 + nperyear
	    endif
        endif
        if ( lwrite ) then
            if ( nperyear.eq.12  ) then
                print '(a,i2,a,i2)','# spectrum: using monhs ',j1,' to '
     +               ,j2
            else
                print '(a,i2,a,i2)','# spectrum: using periods ',j1
     +               ,' to ',j2
            endif
        endif
*
*       sum
*
        if ( lsum.gt.1 ) then
            if ( lwrite ) print '(a)','# Summing data'
            do iens=nens1,nens2
                call sumit(data(1,yrbeg,iens),npermax,nperyear,yrbeg
     +                ,yrend,lsum,oper)
            enddo
        endif
*
*       logscale
*
        if ( logscale ) then
            if ( lwrite ) print '(a,2i3)','# Taking log of series '
            do iens=nens1,nens2
                call takelog(data(1,yrbeg,iens),nperyear,nperyear,yrbeg
     +                ,yrend)
            enddo
        endif
*
*       sqrtscale
*
        if ( sqrtscale ) then
            if ( lwrite ) print '(a,2i3)','# Taking sqrt of series '
            do iens=nens1,nens2
                call takesqrt(data(1,yrbeg,iens),nperyear,nperyear,yrbeg
     +               ,yrend)
            enddo
        endif
*
*       squarescale
*
        if ( squarescale ) then
            if ( lwrite ) print '(a,2i3)','# Taking sqrt of series '
            do iens=nens1,nens2
                call takesquare(data(1,yrbeg,iens),nperyear,nperyear
     +               ,yrbeg,yrend)
            enddo
        endif
*
*       detrending
*       
        if ( ldetrend ) then
            if ( lwrite ) print *,'Detrending data'
            do iens=nens1,nens2
                call detrend(data(1,yrbeg,iens),npermax,nperyear,
     +               yrbeg,yrend,yr1,yr2,m1,m1,1)
            enddo
        endif
        if ( ndiff.ne.0 ) then
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
!
!       anomalies wrt ensemble mean
!
        if ( lensanom .and. min(nens2,mens).gt.max(nens1,mens1) ) then
            if ( lwrite ) print *,'Taking anomalies wrt ensemble mean'
            call anomalensemble(data,npermax,nperyear,yrbeg,yrend,
     +           yr1,yr2,max(nens1,mens1),min(nens2,mens))
        endif

        call getspectrum(data,npermax,nperyear,yrbeg,yrend,yr1,yr2,mens1
     +       ,mens,nens1,nens2,j1,j2,avex,epx,epy,ndata,nout,ofac,mean
     +       ,yrstart,yrstop,lwrite,yr1a,yr2a)
*
*       get error bars by comparing it to an AR(1) process
*
        i = (yrstop-yrstart+1)*(j2-j1+1)
        nmc = 1000000000/i**2
        nmc = min(5000,nmc)
        nmc = max(200,nmc)
        print '(a,3i6)','# i,nmc = ',nout,nmc
        allocate(epx1(4*ndata))
        allocate(epy1(4*ndata))
        allocate(epyall(nmc,nout))
        allocate(prob(nout))
        allocate(imax(nout))
        call getred1(alpha,s2,j1,j2,nperyear,data,npermax,yrbeg,yrend
     +       ,nensmax)
        sd = sqrt(s2)
        if ( lwrite ) print *,'alpha = ',alpha
        do imc=1,nmc
            call keepalive(imc,nmc)
!!!            if ( mod(imc,10).eq.0 ) print '(a,i6)','# ',imc
            call make1mcseries(alpha,sd,data,npermax,nperyear,yrbeg
     +           ,yrend,nens1,nens2,j1,j2,yrstart,yrstop,lwrite)
            call getspectrum(data,npermax,nperyear,yrbeg,yrend,yr1,yr2
     +           ,mens1,mens,nens1,nens2,j1,j2,avex,epx1,epy1,ndata,nout
     +           ,ofac,mean1,yrstart,yrstop,lwrite,yr1a,yr2a)
            do i=nint(ofac),nout
                if ( epx(i).ne.epx1(i) ) then
                    write(0,*) 'warning: x-coordinates not identical ',i
     +                   ,epx(i),epx1(i)
                endif
                epyall(imc,i) = epy1(i)
            enddo
        enddo
        nprob = 0
        probmax = 1
        nimax = 0
        do i=nint(ofac),nout
            call nrsort(nmc,epyall(1,i))
            call getsign1(prob(i),epy(i),epyall(1,i),nmc,lwrite)
            prob(i) = 1-prob(i)
            if ( prob(i).lt.0.05 ) nprob = nprob + 1
            if ( prob(i).lt.probmax ) then
                nimax = 1
                imax(nimax) = i
                probmax = prob(i)
            elseif ( prob(i).eq.probmax ) then
                nimax = nimax + 1
                imax(nimax) = i
            endif
        enddo        
*
*       output; do not yet know how to treat dependent data.
*
        print '(a)','# frequency power AR(1)  prob'
        print '(3a)','# power in [',trim(units),'^2]'
        if ( nperyear.ge.53 ) then
            print '(a)','# frequency in [mo^-1]'
            fac = 12
        else
            print '(a)','# frequency in [yr^-1]'
            fac = 1
        endif
        do i=nint(ofac),nout
            if ( epx(i).gt.0 ) then
                write(*,'(4f12.6)') epx(i)/fac,epy(i)*fac,
     +               epyall(nint(0.95*nmc),i)*fac,prob(i)
            else
                write(0,*) 'spectrum: something is wrong',i,epx(i)
                call abort
            endif
        enddo
        if ( avex.gt.1 ) then
            write(0,'(a)') 'The spectrum has been obtained by averaging'
            write(0,'(i3,a)') avex,' bins of the periodogram<br>'
        endif
        write(0,'(a,i6,a,f6.3,a)')
     +       'The second line denotes the 95% highest spectrum of ',nmc
     +       ,' AR(1) processes with the same autocorrelation (',alpha,
     +       ').'        
        if ( nprob.le.1 ) then
            write(0,'(a,i6)') 'There is ',nprob
        else
            write(0,'(a,i6)') 'There are ',nprob
        endif
        write(0,'(a,f10.1,a)') ' bins with p&lt;0.05, '//
     +       'by chance one would expect ',0.05*(nout-nint(ofac)+1),'.'
        if ( nprob.gt.1.2*0.05*(nout-nint(ofac)+1) ) then
            if ( nimax.eq.1 ) then
                write(0,'(a,f8.2,a,f6.3,a)')
     +               'The most significant peak is at ',1/epx(imax(1))
     +               ,' years (p=',probmax,').'
            else
                write(0,'(a)') 'The most significant peaks are at '
                do i=nimax,2,-1
                    write(0,'(f8.2,a)') ,1/epx(imax(i)),', '
                enddo
                write(0,'(a,f8.2,a,f6.3,a)') ' and ',1/epx(imax(1)),
     +               ' years (p=',probmax,').'
            endif
        endif
        if ( yr1a.gt.yrbeg .or. yr2a.lt.yrend ) then
            write(0,'(2a,i4,a,i4,a,f16.3)')
     +           'The mean of the spectrum exp(&lang;log(1/f)&rang;) ',
     +           ' band-pass filtered between ',yr1a-yrbeg,' and ',
     +           yr2a-yrbeg,' yrs is at period ',mean, ' yr '

        else
            write(0,'(2a,f16.3)')
     +           'The mean of the spectrum exp(&lang;log(1/f)&rang;) ',
     +           'is at period ',mean, ' yr '
        end if
        call savestartstop(yrstart,yrstop)
        goto 999
  900   print *,'error: cannot read integer from ',line
        call abort
  999   continue
        end
