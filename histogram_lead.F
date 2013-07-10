*  #[ histogram:
        program histogram
*
*       make a histogram of timeseries
*       and do a maximum-likelihood fit if requested
*
        implicit none
#include "param.inc"
#include "getopts.inc"
        integer maxbin,maxldat,nboot
        parameter (maxbin=1000,maxldat=1000000,nboot=800)

        integer i,ii,j,jj,j1,j2,n,nperyear,nbin,month,yr,nn(maxbin)
     +       ,nless,nfit,ntot,ntype,ier,iens,mens1,mens,lead,mop,yrp,l
     +       ,iboot,mboot,ndecor,nfitted,t,tmax
        logical lf(maxldat),lastvalid
        real data(npermax,yrbeg:yrend,0:nensmax,leadmax),mindat,maxdat
     +       ,d,a,b,c,da,db,dc,xi,dxi,xx(maxldat),yy(maxbin),x,x1,x2,s
     +       ,sqrt2,snorm,f,offset,xn(maxbin),xy(maxbin),df,chsq,prob
        real mean(0:nboot),sd(0:nboot),skew(0:nboot)
        character file*255,string*100

        integer iargc
        real gammln,erf,erfcc,gammq,gammp,invcumpois,invcumgamm
        real *8 ranf
        external gammln,erf,erfcc,gammp,gammq,invcumpois
     +        ,invcumgamm,ranf
*
        n = iargc()
        if ( n.lt.2 ) then
            print *,'usage: histogram timeseries nbins'//
     +           ' hist hist|qq|gumbel'//
     +           ' [fit poisson|gauss|gamma|gumbel|gev] [options]'
            stop
        endif
*
        call getarg(1,file)
        print '(2a)','# file ',file(1:index(file,' ')-1)
        call readleadensdat(file,data,npermax,yrbeg,yrend,nensmax
     +        ,leadmax,nperyear,mens1,mens)
        call getarg(2,string)
        read(string,*) nbin
        if ( nbin.gt.maxbin ) then
            write(0,*) 'histogram: error: increase maxbin'
            write(*,*) 'histogram: error: increase maxbin'
            call abort
        endif
        i = 3
        nfit = 0
        ntype = 0
        if ( n.gt.2 ) then
 10         continue
            call getarg(i,string)
            if ( string(1:3).eq.'fit' ) then
                call getarg(i+1,string)
                i = i + 2
                if ( string(1:2).eq.'no' ) then
                    nfit = 0
                elseif ( string(1:3).eq.'poi' ) then
                    nfit = 1
                elseif ( string(1:3).eq.'gau' ) then
                    nfit = 2
                elseif ( string(1:3).eq.'gam' ) then
                    nfit = 3
                elseif ( string(1:3).eq.'gum' ) then
                    nfit = 4
                elseif ( string(1:3).eq.'gev' ) then
                    nfit = 5
                else
                    write(0,*) 'histogram: error: unknown distribution '
     +                    ,string
                    write(*,*) 'histogram: error: unknown distribution '
     +                    ,string
                    call abort
                endif
                goto 10
            endif
            if ( string(1:4).eq.'hist' ) then
                call getarg(i+1,string)
                i = i + 2
                call tolower(string)
                if ( string(1:4).eq.'hist' ) then
                    ntype = 0
                elseif ( string(1:3).eq.'qq' ) then
                    ntype = 1
                elseif ( string(1:4).eq.'gumb' ) then
                    ntype = 2
                else
                    write(0,*) 'histogram: error: unknown plot ',string
                    write(*,*) 'histogram: error: unknown plot ',string
                    call abort
                endif
                goto 10
            endif
            call getopts(i,n,nperyear,yrbeg,yrend,.true.,mens1,mens)
            if ( mens.gt.0 ) then
                print '(a,i3,a,i3)','# using ensemble members ',nens1,
     +                ' to ',nens2
            endif
            if ( lead2.gt.leadmax ) then
                write(0,*) 'error: lead2 > leadmax: ',lead2,leadmax
                call abort
            endif
        endif
        if ( plot .and. .not. lbootstrap ) then
            mboot = 0
        else
            mboot = nboot
        endif
*       
*       take max
*       
        if ( lsum.gt.1 .and. oper.eq.'a' ) then
            print '(a,i4,a)','# taking max of ',lsum,' periods'
            do iens=nens1,nens2
                do lead=lead1,lead2
                    if ( m1.ne.m2 ) then
                        write(0,*) ' histogram: error: can only take '//
     +                        'max for one month'
                        call abort
                    endif
                    do yr=yr1,yr2-1
                        do month=m1+1,m1+lsum-1
                            mop = month
                            call normon(mop,yr,yrp,nperyear)
                            if ( data(mop,yrp,iens,lead).lt.1e33 ) then
                                data(month,yr,iens,lead) = max(
     +                                data(month,yr,iens,lead),
     +                                data(mop,yrp,iens,lead))
                                data(mop,yrp,iens,lead) = 3e33 ! safety
                            else ! search at other lead times
                                do l=lead+1,lead2
                                    if ( data(mop,yrp,iens,l).lt.1e33 )
     +                                    then
                                        data(month,yr,iens,lead) = max(
     +                                        data(month,yr,iens,lead),
     +                                        data(mop,yrp,iens,l))
                                        data(mop,yrp,iens,l) = 3e33 ! not twice!
                                        goto 100
                                    endif
                                enddo
                                data(month,yr,iens,lead) = 3e33
  100                           continue
                            endif ! possible?
                        enddo   ! lsum
                    enddo       ! yr
                enddo           ! lead
            enddo               ! iens
            lsum = 1            ! otherwise the sum kicks in...
        endif
*
*       sum series
*       
        do lead=lead1,lead2
            do iens=nens1,nens2
                if ( lsum.gt.1 ) then
                    if ( iens.eq.nens1 ) print '(a,i4,a)'
     +                    ,'# summing over ',lsum,' periods'
                    call sumit(data(1,yrbeg,iens,lead),npermax,nperyear
     +                    ,yrbeg,yrend,lsum,oper)
                endif
*       
*               apply minindx, maxindx cuts
*               In case of a Poisson distribution round to unambiguous
*               numbers
*       
                if ( nfit.eq.1 ) then
                    if ( minindx.ge.0 ) minindx = 0.5 + int(minindx+1e-5
     +                    )
                    if ( maxindx.lt.2.**31 ) maxindx = 0.5 + int(maxindx
     +                    -1e-5)
                endif
                if ( maxindx.lt.1e33 .or. minindx.gt.-1e33 ) then
                    do i=yr1,yr2
                        do j=1,nperyear
                            if (  data(j,i,iens,lead).gt.maxindx .or.
     +                            data(j,i,iens,lead).lt.minindx ) then
                                data(j,i,iens,lead) = 3e33
                            endif
                        enddo
                    enddo
                endif
*       
*               log, sqrt
*
                if ( logscale ) then
                    if ( iens.eq.nens1 ) print '(a)'
     +                    ,'# taking logarithm'
                    call takelog(data(1,yrbeg,iens,lead),npermax
     +                    ,nperyear,yrbeg,yrend)
                endif
                if ( sqrtscale ) then
                    if ( iens.eq.nens1 ) print '(a)','# taking sqrt'
                    call takesqrt(data(1,yrbeg,iens,lead),npermax
     +                    ,nperyear,yrbeg,yrend)
                endif
                if ( squarescale ) then
                    if ( iens.eq.nens1 ) print '(a)','# taking square'
                    call takesquare(data(1,yrbeg,iens,lead),npermax
     +                    ,nperyear,yrbeg,yrend)
                endif
                if ( lchangesign ) then
                    if ( iens.eq.nens1 ) print '(a)','# changing sign'
                    call changesign(data(1,yrbeg,iens,lead),npermax
     +                    ,nperyear,yrbeg,yrend)
                endif
*
*               detrend data
*
                if ( ldetrend ) then
                    if ( iens.eq.nens1 ) print '(a,i4)'
     +                    ,'# detrending data'
                    call detrend(data(1,yrbeg,iens,lead),npermax
     +                   ,nperyear,yrbeg,yrend,yr1,yr2,m1,m2,lsel)
                endif
*
*               differentiate data
*
                if ( ndiff.ne.0 ) then
                    if ( iens.eq.nens1 ) print '(a,i4)'
     +                    ,'# taking differences/averaging ',ndiff
                    call diffit(data(1,yrbeg,iens,lead),npermax,nperyear
     +                    ,yrbeg,yrend,ndiff)
                endif
*
*               anomalies - necessary if we consider more than one month
*
                if ( anom ) then
                    if ( iens.eq.nens1 ) print '(a)'
     +                    ,'# taking anomalies'
                    call anomal(data(1,yrbeg,iens,lead),npermax,nperyear
     +                    ,yrbeg,yrend,yr1,yr2)
                endif
            enddo
        enddo
*
*       compute min, max
*
        do month=m1,m2
            print '(a,i2)','# month ',month
            if ( month.eq.0 ) then
                j1 = 1
                j2 = nperyear
            else
                j1 = month
                j2 = month + lsel - 1
                if ( nperyear.gt.12 ) then
                    call month2period(j1,nperyear,1)
                    call month2period(j2,nperyear,0)
                endif
            endif
            ntot = 0
            do iboot=0,mboot
                mean(iboot) = 0
                sd(iboot) = 0
                skew(iboot) = 0
            enddo
            mindat = +3e33
            maxdat = -3e33
            offset = 0
            do iens=nens1,nens2
                do lead=lead1,lead2
                    lastvalid = .false.
                    do yr=yr1-1,yr2
                        if ( j1.ne.j2 .and. j2-j1+1.ne.nperyear ) then
                            lastvalid = .false.
                        endif
                        do jj=j1,j2
                            j = jj
                            call normon(j,yr,i,nperyear)
                            if ( i.ge.yr1 .and. i.le.yr2) then
                                if (  data(j,i,iens,lead).lt.1e33 ) then
                                    ntot = ntot + 1
                                    if ( ntot.gt.maxldat ) then
                                        print *,'histogram: error '//
     +                                        'increase maxldat'
                                        write(0,*) 'histogram: error: '
     +                                        //'increase maxldat'
                                        call abort
                                    endif
                                    xx(ntot) = data(j,i,iens,lead)
                                    offset = offset + xx(ntot)
                                    lf(ntot) = .not.lastvalid
                                    lastvalid = .true.
                                    if ( lwrite .and. lf(ntot) ) print *
     +                                   ,'boundary at ',ntot,j,i,iens
     +                                   ,lead
                                else
                                    lastvalid = .false.
                                endif
                            else
                                lastvalid = .false.
                            endif
                        enddo
                    enddo
                enddo
            enddo
            if ( ntot.eq.0 ) then
                write(*,'(a)') '# histogram: no valid points'
                stop
            endif
            offset = offset/ntot
            if ( lwrite ) print '(a,g20.10)','# offset = ',offset
            do i=1,ntot
                x = xx(i) - offset
                mean(0) = mean(0) + x
                sd(0) = sd(0) + x**2
                skew(0) = skew(0) + x**3
                mindat = min(mindat,xx(i))
                maxdat = max(maxdat,xx(i))
            enddo
*
*           enough points?
*
            if ( ntot.lt.minnum ) then
                print *,'not enough points'
                goto 999
            endif
*
*           compute errors with a bootstrap
*
            if ( j1.eq.j2 ) then
                ndecor = 1+nint(decor/nperyear)
            else
                ndecor = 1+nint(decor)
            endif
            do iboot=1,mboot
                n = ntot/ndecor
                do i=1,n
                    ii = 0
 200                continue
                    ii = ii + 1
                    if ( ii.gt.100 ) then
                        write(*,'(a,i4,a,i6)') '# histogram: error: '//
     +                       'cannot find blocks of size ',ndecor
     +                       ,' for moving block bootstrap ',iboot
                        goto 300
                    endif
                    j = 1+min(ntot-ndecor,
     +                   int((ntot-ndecor)*ranf(i*iboot)))
*                   check whether there was a break in this segment
                    do jj=1,ndecor-1
                        if ( lf(j+jj) ) then
                            if ( lwrite ) print *
     +                           ,'crossing boundary at ',j+jj,ii,iboot
                            goto 200
                        endif
                    enddo
                    do jj=0,ndecor-1
                        x = xx(j+jj) - offset
                        mean(iboot) = mean(iboot) + x
                        sd(iboot) = sd(iboot) + x**2
                        skew(iboot) = skew(iboot) + x**3
                    enddo
                enddo
            enddo
 300        continue
            mboot = iboot-1
            do iboot=0,mboot
                if ( iboot.eq.0 ) then
                    n = ntot
                else
                    n = ndecor*(ntot/ndecor)
                endif
                mean(iboot) = mean(iboot)/n
                sd(iboot) = sd(iboot)/n - mean(iboot)**2
                if ( sd(iboot).gt.0 ) then
                    sd(iboot) = sqrt(sd(iboot))
                else
                    print '(a,2g20.10,i6)','# error: variance<0 '
     +                   ,mean(iboot),sd(iboot),iboot
                    do i=1,n
                        print '(a,i6,g20.10)','# ',i,xx(i)
                    enddo
                    sd(iboot) = 0
                endif
                skew(iboot) = skew(iboot)/n - 3*sd(iboot)**2
     +               *mean(iboot) - mean(iboot)**3
                if ( sd(iboot).gt.0 ) then
                    skew(iboot) = skew(iboot)/sd(iboot)**3
                else
                    skew(iboot) = 3e33
                endif
                mean(iboot) = mean(iboot) + offset
            enddo
            call nrsort(mboot,mean(1))
            call nrsort(mboot,sd(1))
            call nrsort(mboot,skew(1))
            print '(a,i16)','# N:         ',ntot
            if ( ndecor.eq.1 ) then
                print '(3a)'
     +               ,'# The error margins were computed with a '
     +               ,'bootstrap method that assumes all points are '
     +               ,'independent'
            elseif ( j1.ne.j2 ) then
                print '(2a,i4,a)'
     +               ,'# The error margins were computed with a '
     +               ,'moving block bootstrap with block size '
     +               ,ndecor,' months'
            else
                print '(2a,i4)'
     +               ,'# The error margins were computed with a '
     +               ,'moving block bootstrap with block size '
     +               ,ndecor,' years'
            endif
            if ( lweb ) then
                print '(8a)'
     +               ,'# <tr><th>parameter</th><th>value</th><th>'
     +               ,'95% CI</th></tr>'
            else
                print '(8a)',       '#        ',
     +               '                ',
     +               '            2.5%',
     +               '             16%',
     +               '             50%',
     +               '             84%',
     +               '           97.5%',
     +               '          +/-95%'
            endif
            call printvalerr('# mean:      ',mean,mboot,plot,11,lweb)
            call printuntransf(mean(0))
            call printvalerr('# s.d.(n):   ',sd,mboot,plot,11,lweb)
            call printval   ('# s.d.(n-1): ',sd(0)*sqrt(real(ntot)
     +            /real(ntot-1)),plot,11,lweb)
            call printvalerr('# skew:      ',skew,mboot,plot,11,lweb)
            call printval   ('# min:       ',mindat,plot,11,lweb)
            call printuntransf(mindat)
            call printval   ('# max:       ',maxdat,plot,11,lweb)
            call printuntransf(maxdat)
*
*           fit to distribution
*       
            if ( nfit.eq.0 ) then
*               no fit requested
                snorm = 1
            elseif ( nfit.eq.1 ) then
*               Poisson distribution
                call fitpoi(xx,ntot,mean(0),a)
                call poisnorm(a,snorm)
                nfitted = 2
            elseif ( nfit.eq.2 ) then
*               Gaussian distribution
                if ( minindx.lt.-1e33 .and. maxindx.gt.1e33 ) then
                    a = mean(0)
                    b = sd(0)
                    print '(a)','# for Gaussian used mean, sd.'
                else
                    call fitgau(xx,ntot,mean(0),sd(0),a,b)
                endif
                call gausnorm(a,b,snorm)
                nfitted = 3
            elseif ( nfit.eq.3 ) then
*               Gamma distribution
                call fitgam(xx,ntot,mean(0),sd(0),a,b)
                if ( minindx.gt.-1e33 ) then
                    if ( maxindx.lt.1e33 ) then
                        snorm = gammp(a,maxindx/b) - gammp(a,minindx/b)
                    else
                        snorm = gammq(a,minindx/b)
                    endif
                else
                    if ( maxindx.lt.1e33 ) then
                        snorm = gammp(a,maxindx/b)
                    else
                        snorm = 1
                    endif
                endif
                nfitted = 3
            elseif ( nfit.eq.4 ) then
*               Gumbel distribution
                call fitgum(xx,ntot,mean(0),sd(0),a,b,j1,j2,lweb,ntype
     +               ,lchangesign)
                call gumbnorm(a,b,snorm)
                nfitted = 3
            elseif ( nfit.eq.5 ) then
*               GEV distribution
                call fitgev(xx,ntot,mean(0),sd(0),a,b,xi,j1,j2,lweb
     +               ,ntype,lchangesign)
                call gevnorm(a,b,xi,snorm)
                nfitted = 4
            else
                write(0,*) 'histogram: error: unknow distribution ',nfit
                call abort
            endif
            if ( ntype.eq.0 ) then
*
*               adjust binsize to nice numbers (heuristics)
*
*               if the minimum is close to zero make exactly zero
                if ( mindat.gt.0 .and. mindat.lt.maxdat/3 ) mindat = 0
                if ( maxdat.lt.0 .and. maxdat.gt.mindat/3 ) maxdat = 0
*               if the minimum and maxmimum are similar in absolute value 
*		make them equal
		if ( mindat.lt.0 .and. maxdat.gt.0 .and. 
     +			abs(mindat+maxdat)/2.lt.maxdat/3 ) then
		    if ( mindat.lt.-maxdat ) then
			maxdat = -mindat
		    else
			mindat = -maxdat
		    endif
		endif
                d = (maxdat-mindat)/(nbin-1)
                s = d
                if ( d.ge.1 ) then
                    i = int(log10(d))
                else
                    i = int(log10(d)) - 1
                endif
                d = d/10.**i    ! now d is between 1 and 10
                if ( d.gt.5 ) then
                    d = 10
                elseif ( d.gt.2 ) then
                    d = 5
                elseif ( d.gt.1 ) then
                    d = 2
                else
                    d = 1
                endif
                d = d*10.**i
                if ( lwrite ) then
                    print '(a,f16.4,a,f16.4)','# binsize rounded from '
     +                    ,s,' to ',d
                endif
                s = mindat
                if ( mindat.ge.0 ) then
                    mindat = d*(int(mindat/d))
                elseif ( mindat.ne.-maxdat ) then
                    mindat = d*(-1+int(mindat/d))
                else
                    mindat = -nbin*d/2
                endif
                if ( lwrite ) then
                    print '(a,f16.4,a,f16.4)','# minimum rounded from '
     +                    ,s,' to ',mindat
                endif
*
*               fill histogram
*
                do i=1,nbin
                    nn(i) = 0
                enddo
                nless = 0
                do yr=yr1-1,yr2
                    do jj=j1,j2
                        j = jj
                        call normon(j,yr,i,nperyear)
                        if ( i.gt.yr1 .and. i.lt.yr2 ) then
                            do lead=lead1,lead2
                                do iens=nens1,nens2
                                    if ( data(j,i,iens,lead).lt.1e33 )
     +                                    then
                                        n = 1 + int((data(j,i,iens,lead)
     +                                        -mindat)/d)
                                        if ( n.lt.0 .or. n.gt.nbin )
     +                                        then
                                            write(0,*)
     +                                            'histogram: error
     +                                            :'//' n is ',n,nbin
                                            write(*,*)
     +                                            'histogram: error
     +                                            :'//' n is ',n,nbin
                                            write(*,*) j,i,iens,lead,
     +                                           data(j,i,iens,lead),
     +                                           mindat,maxdat,d
                                            call abort
                                        endif
                                        nn(n) = nn(n) + 1
                                        if ( data(j,i,iens,lead)
     +                                        .lt. mean(0)) 
     +                                        nless= nless + 1
                                    endif
                                enddo
                            enddo
                        endif
                    enddo
                enddo
*       
*               fill fitted curve array
*
                if ( nfit.eq.0 ) then
*                   no fit requested
                    do i=1,nbin
                        yy(i) = 0
                    enddo
                elseif ( nfit.eq.1 ) then
*                   Poisson distribution
                    do i=1,nbin
                        x1 = mindat + d*(i-1)
                        yy(i) = 0
                        do j=1,nint(d)
                            x2 = x1 + (j-1)
                            yy(i) = yy(i) + ntot*
     +                            exp(x2*log(a)-a-gammln(1+x2))
                        enddo
                        yy(i) = yy(i)/snorm
                    enddo
                elseif ( nfit.eq.2 ) then
*                   Gaussian distribution
                    sqrt2 = sqrt(2.)
                    do i=1,nbin
                        x1 = (mindat + d*(i-1) - a)/b
                        x2 = (mindat + d*i     - a)/b
                        if ( x2.gt.0 ) then
                            yy(i) = ntot*(erfcc(x1/sqrt2) - 
     +                            erfcc(x2/sqrt2))/2
                        else
                            yy(i) = ntot*(erfcc(-x2/sqrt2) - 
     +                            erfcc(-x1/sqrt2))/2
                        endif
                        yy(i) = yy(i)/snorm
                    enddo
                elseif ( nfit.eq.3 ) then
*                   Gamma distribution
                    do i=1,nbin
                        x1 = (mindat + d*(i-1))/b
                        x2 = (mindat + d*i    )/b
                        if ( x2.lt.a-1 ) then
                            s = gammp(a,x2) - gammp(a,x1)
                        else
                            s = gammq(a,x1) - gammq(a,x2)
                        endif
                        yy(i) = ntot*s/snorm
                    enddo                
                elseif ( nfit.eq.4 ) then
*                   Gumbel distribution
                    do i=1,nbin
                        x1 = (mindat + d*(i-1) - a)/b
                        x2 = (mindat + d*i     - a)/b
                        yy(i) = ntot*(exp(-exp(-x2))-exp(-exp(-x1)))
     +                        /snorm
                    enddo
                elseif ( nfit.eq.5 ) then
*                   GEV distribution
                    do i=1,nbin
                        x1 = (mindat + d*(i-1) - a)/b
                        x2 = (mindat + d*i     - a)/b
                        if ( xi.eq.0 ) then
                            yy(i) = ntot*(exp(-exp(-x2))-exp(-exp(-x1)))
     +                            /snorm
                        else
                            if ( 1+xi*x1.le.0 .or. 1+xi*x2.le.0 ) then
                                yy(i) = 0
                            else
***                                write(0,*) x1,x2,xi,snorm
                                yy(i) = ntot*(exp(-(1+xi*x2)**(-1/xi)) -
     +                                exp(-(1+xi*x1)**(-1/xi)))/snorm
                            endif
                        endif
                    enddo
                else
                    write(0,*) 'histogram: error: unknown distribution '
     +                    ,nfit
                    call abort
                endif
*
*               compute \chi^2 and P of fit
*
                if ( nfit.ne.0 ) then
*                   leave out trivial points with N=0, n~0 in the tails
                    n = 0
                    do i=1,nbin
                        if ( yy(i).gt.0.1 .or. nn(i).gt.0 .or. n.gt.0 )
     +                       then
                            n = n + 1
                            xn(n) = nn(i)
                            xy(n) = yy(i)
                        endif                        
                    enddo
                    do i=n,1,-1
                        if ( xy(i).gt.0.1 .or. xn(i).gt.0 ) goto 410
                    enddo
 410                continue
                    n = i
                    if ( lwrite ) then
                        do i=1,n
                            print '(a,i4,2f8.2,g9.2)','# ',i,xn(i),xy(i)
     +                           ,(xn(i)-xy(i))**2/xy(i)
                        enddo
                    endif
*                   Numerical Recipes, p615
                    call chsone(xn,xy,n,nfitted,df,chsq,prob)
                    if ( plot ) write(11,'(g16.4)') chsq/df
                    if ( lweb ) then
                        print '(a,g9.2,a,i4,a,g9.3,a,f6.4)',
     +                       '# <tr><td>&chi;^2/df</td><td>'
     +                       ,chsq,'/',nint(df),' =',chsq/df
     +                       ,'</td><td>p=',prob
                    else
                        print '(a,g9.2,a,i4,a,g9.3)','# chi^2/df:',
     +                       chsq,'/',nint(df),' =',chsq/df
                        print '(a,f8.4,a)','# probability:   ',100*prob,
     +                   '% (this is the probability that the observed '
     +                   //'distribution was drawn from the fitted one)'
                    endif
                    if ( plot ) write(11,'(g16.4)') prob
                endif
*
*               print out
*
                n = 0
                s = 0
                do i=1,nbin
                    n = n + nn(i)
                    s = s + yy(i)
                    print '(i4,f12.4,2i9,2f12.2)',i,mindat+ d*(i-0.5)
     +                    ,nn(i),n,yy(i),s
                enddo
***                print '(a,i8,a,f6.2a)','# there are ',nless,' (',100
***     +                *nless/real(n),'%) points less than the mean'

            elseif ( ntype.eq.1 ) then
*       
*               QQ plot
*       
                call nrsort(ntot,xx)
                do i=1,ntot
                    f = real(i)/real(ntot+1)
                    if ( nfit.eq.0 ) then
*                       no fit requested
                        write(0,*) 'histogram: error: can only make ',
     +                        'QQ plot when fitting to a distribution'
                        write(*,*) 'histogram: error: can only make ',
     +                        'QQ plot when fitting to a distribution'
                        call abort
                    elseif ( nfit.eq.1 ) then
*                       Poisson distribution - only the last point of a bin
*			makes sense
			if ( xx(i).ne.xx(i+1) ) then
                            f = snorm*f
                            if ( minindx.gt.0 ) then
                                f = f + gammq(minindx+0.5,a)
			    endif
                            s = invcumpois(f,a)
			else
			    s = 3e33
                        endif
                    elseif ( nfit.eq.2 ) then
*                       Gaussian distribution
                        sqrt2 = sqrt(2.)
                        f = 2*snorm*f
                        if ( minindx.gt.-1e33 ) then
                            f = f + erf((minindx-a)/(sqrt2*b))
                        else
                            f = f - 1
                        endif
*                       IMSL routine
                        call merfi(f,x,ier)
                        s = a + sqrt2*b*x
                    elseif ( nfit.eq.3 ) then
*                       Gamma distribution
                        f = snorm*f
                        if ( minindx.gt.0 ) then
                            f = f + gammp(a,minindx/b)
                        endif
                        s = invcumgamm(f,a,b)
                    elseif ( nfit.eq.4 ) then
*                       Gumbel distribution
                        s = snorm*f
                        if ( minindx.gt.-1e33 ) then
                            s = s + exp(-exp(-(minindx-a)/b))
                        endif
                        s = a - b*log(-log(s))
                    elseif ( nfit.eq.5 ) then
*                       GEV distribution
                        s = snorm*f
                        if ( minindx.gt.-1e33 ) then
                            s = s + exp(-(1+(minindx-a)/b)**(-1/xi))
                        endif
                        if ( xi.eq.0 ) then
                            s = a - b*log(-log(s))
                        else
                            s = a + b/xi*((-log(s))**(-xi)-1)
                        endif
                    else
                        write(0,*)
     +                        'histogram: error: unknown distribution '
     +                        ,nfit
                        call abort
                    endif
                    if ( s.lt.1e33 ) then
                        x = xx(i)
                        if ( lchangesign ) then
                            if ( x.ne.-999.9 ) x = -x
                        endif
                        print '(i8,2g22.6)',i,x,s
                    endif
                enddo
            elseif ( ntype.eq.2 ) then
*       
*               Gumbel plot
*       
                call nrsort(ntot,xx)
                do i=1,ntot+100
                    if ( i.le.ntot ) then
                        f = real(i)/real(ntot+1)
                        x = xx(i)
                    else
                        f = 1 - 1/real(ntot+1)*0.9**(i-ntot)
                        if ( abs(f-1).lt.1e-6 ) goto 800
                        x = -999.9
                    endif
                    if ( nfit.eq.0 ) then
                        s = -999.9
                    elseif ( nfit.eq.1 ) then
*                       Poisson distribution - only the last point
*                       of a bin makes sense
                        if ( i.gt.ntot .or. 
     +                       xx(min(ntot,i)).ne.xx(min(ntot,i+1)) ) then
                            f = snorm*f
                            if ( minindx.gt.0 ) then
                                f = f + gammq(minindx+0.5,a)
                            endif
                            s = invcumpois(f,a)
                        else
                            s = -999.9
                        endif
                    elseif ( nfit.eq.2 ) then
*                       Gaussian distribution
                        sqrt2 = sqrt(2.)
                        f = 2*snorm*f
                        if ( minindx.gt.-1e33 ) then
                            f = f + erf((minindx-a)/(sqrt2*b))
                        else
                            f = f - 1
                        endif
*                       IMSL routine
                        call merfi(f,x,ier)
                        s = a + sqrt2*b*x
                    elseif ( nfit.eq.3 ) then
*                       Gamma distribution
                        f = snorm*f
                        if ( minindx.gt.0 ) then
                            f = f + gammp(a,minindx/b)
                        endif
                        s = invcumgamm(f,a,b)
                    elseif ( nfit.eq.4 ) then
*                       Gumbel distribution
                        s = snorm*f
                        if ( minindx.gt.-1e33 ) then
                            s = s + exp(-exp(-(minindx-a)/b))
                        endif
                        s = a - b*log(-log(s))
                    elseif ( nfit.eq.5 ) then
*                       GEV distribution
                        s = snorm*f
                        if ( minindx.gt.-1e33 ) then
                            s = s + exp(-(1+(minindx-a)/b)**(-1/xi))
                        endif
                        if ( xi.eq.0 ) then
                            s = a - b*log(-log(s))
                        else
                            s = a + b/xi*((-log(s))**(-xi)-1)
                        endif
                    else
                        write(0,*) 'histogram: error: '//
     +                       'unknown distribution ',nfit
                        call abort
                    endif
                    if ( lchangesign ) then
                        if ( x.ne.-999.9 ) x = -x
                    endif
                    print '(i8,3g22.6)',i,-log(-log(f)),x,s
                    if ( i.gt.ntot .and. (1-f)*(j2-j1+1).lt.0.0001 ) 
     +                   goto 800
                enddo
 800            continue
*
*               convert Gumbel variates to return periods
                tmax = max(10000,int((ntot+1)/real(j2-j1+1))) ! in years
                print '(a,f8.4,a)','#@ set xrange [:',
     +               -log(-log(1-1/(tmax*real(j2-j1+1)+1))),']'
                print '(a)','#@ set xtics (\\'
                t = 1
                do i=1,100
                    f = 1-1/(t*real(j2-j1+1)+1)
                    write(string,'(a,i1,a)') '(a,i',
     +                   1+int(log10(real(t))),',a,f8.4,a)'
                    if ( t.gt.10 .and. mod(i,3).ne.1 ) then
                        print '(a,f8.4,a)','#@ "" ',-log(-log(f)),
     +                   '\\'
                    else
                        print string,'#@ "',t,'" ',-log(-log(f)),
     +                   '\\'
                    endif
                    if ( mod(i,3).eq.2 ) then
                        t = nint(2.5*t)
                    else
                        t = 2*t
                    endif
                    if ( t.gt.tmax ) goto 801
                    print '(a)','#@ ,\\'
                enddo
 801            continue
                print '(a)','#@ )'
            else
                write(0,*) 'Unknown plot type ',ntype
                write(*,*) 'Unknown plot type ',ntype
                call abort
            endif
        enddo
 999    continue
        end
*  #] histogram:
*  #[ fitpoi:
        subroutine fitpoi(xx,ntot,mean,a)
*       
*       a fit a Poisson distribution to the data
*       
        implicit none
*        
        integer ntot
        real xx(ntot),mean,a
*
        integer i,nx
        real ax,bx,cx,fa,fb,fc,tol
*
        integer nmax,ncur
        parameter(nmax=100000)
        real data(nmax),restrain
        common /fitdata1/ data
        common /fitdata2/ restrain,ncur
*
        real llpoisson
        external llpoisson
*
*       check input
*       
        do i=1,ntot
            if ( abs(xx(i)-nint(xx(i))).gt.0.01 ) then
                write(0,*) 'histogram: error: cannot fit Poisson to '//
     +                'non-integer data: ',xx(i)
                call abort
            endif
        enddo
*
*       copy to common for routine llpoisson
*       
        ncur = ntot
        do i=1,ncur
            data(i) = xx(i)
        enddo
*
*       fit, using Numerical Recipes routines
*
        ax = mean/2
        bx = mean
        call mnbrak(ax,bx,cx,fa,fb,fc,llpoisson)
        tol = 1e-6
        call brent(ax,bx,cx,llpoisson,tol,a)
*
*       output
*
        print '(a)','# Fitted to Poisson distribution'
        print '(a)','# p(n) = mu^n*exp(-mu)/n! with'
        print '(a,f16.3,a)','# mu = ',a
        end
*  #] fitpoi:
*  #[ llpoisson:
        real function llpoisson(a)
*
*       computes the log-likelihood function for a Poisson distribution
*       with parameter a and data in common.
*       
        implicit none
*       
        real a
*       
        integer i
        real s
*
        integer nmax,ncur
        parameter(nmax=100000)
        real data(nmax),restrain
        common /fitdata1/ data
        common /fitdata2/ restrain,ncur
*
        real gammln
        external gammln
*
        llpoisson = 0
        do i=1,ncur
            llpoisson = llpoisson + 
     +            data(i)*log(a) - a - gammln(1+data(i))
        enddo
*       normalization is not 1 in case of cut-offs
        call poisnorm(a,s)
        llpoisson = llpoisson - ncur*log(s)
*       minimum, not maximum
        llpoisson = -llpoisson
*        
        end
*  #] llpoisson:
*  #[ poisnorm:
        subroutine poisnorm(a,s)
*
*       get the normalization of the poisson function on the interval 
*       minindx,maxindx
*
        implicit none
#include "getopts.inc"
        real a,s
        integer i
        real xk1,xk2
        real poisson,gammp,gammq
        external poisson,gammp,gammq
*
        if ( minindx.gt.0 ) then
            xk1 = nint(minindx+0.5)
            if ( maxindx.lt.2.**31 ) then
                xk2 = nint(maxindx-0.5) + 1
                s = gammq(xk2,a) - gammq(xk1,a)
            else
                s = gammp(xk1,a)
            endif
        else
            if ( maxindx.lt.1e33 ) then
                xk2 = nint(maxindx-0.5) + 1
                s = gammq(xk2,a)
            else
                s = 1
            endif
        endif
***        print *,'poisnorm: norm = ',s
        end
*  #] poisnorm:
*  #[ poisson:
        real function poisson(a,n)
        implicit none
        real a
        integer n
        real factrl
        external factrl
*
        poisson = a**n*exp(-a)/factrl(n)
***        print *,'poisson(',a,n,') = ',poisson
        end
*  #] poisson:
*  #[ cumpois:
        real function cumpois(x)
*
*       compute the cumulative Poisson probability minus some
*       requested frequency.  Parameters are passed in common
*
        implicit none
        real x
        real pc,muc
        common /ccumpois/ pc,muc
        real gammq
        external gammq
*
        if ( x.lt.0 ) then
            cumpois = -pc
        else
            cumpois = gammq(x+1,muc) - pc
        endif
*
        end
*  #] cumpois:
*  #[ invcumpois:
        real function invcumpois(p,mu)
*
*       compute the inverse of the cumulative Poisson distribution Q(k,mu)
*       as long as I do not find or make an explicit function just solve
*       the equation.
*
        implicit none
        real p,mu
        integer i
        real x,x1,x2,tol
        real pc,muc
        common /ccumpois/ pc,muc
        real cumpois,zbrent
        external cumpois,zbrent
*
*       check argument
*
        if ( p.eq.0 ) then
            x = 0
            goto 999
        elseif ( p.lt.0 .or. p.ge.1 ) then
            write(0,*) 'invcumpois: illegal argument ',p
            x = 3e33
            goto 999
        endif
*
*       parameters for function cumpois
*
        muc = mu
        pc = p
*
*       bracket zero
*
        x1 = 0
        x2 = mu
        i = 0
  100   continue
        if ( cumpois(x2).lt.0 ) then
            x2 = 1.6*x2
            i = i + 1
            if ( i.gt.100 ) then
                write(0,*) 'invcumpois: error: cannot find root'
                x = 3e33
                goto 999
            endif
            goto 100
        endif
*
*       get root
*
        tol = 1e-5*x2
        x = zbrent(cumpois,x1,x2,tol)
*
*       finito
*
  999   continue
        invcumpois = x
        return
        end
*  #] invcumpois:
*  #[ fitgau:
        subroutine fitgau(xx,ntot,mean,sd,a,b)
*
*       a fit a gaussian distribution to the data
*
        implicit none
*
        integer ntot
        real xx(ntot),mean,sd,a,b
*
        integer i,nx,iter
        real tol,p(3,2),q(2),xmin,y(3)
*
        integer nmax,ncur
        parameter(nmax=100000)
        real data(nmax),restrain
        common /fitdata1/ data
        common /fitdata2/ restrain,ncur
*
        real llgauss
        external llgauss
*
*       copy to common for routine llgauss
*       
        ncur = ntot
        do i=1,ncur
            data(i) = xx(i)
        enddo
*       
*       fit, using Numerical Recipes routines
*
        p(1,1) = mean *0.9
        p(1,2) = sd   *0.9
        p(2,1) = p(1,1) *1.2
        p(2,2) = p(1,2)
        p(3,1) = p(1,1)
        p(3,2) = p(1,2) *1.2
        do i=1,3
            q(1) = p(i,1)
            q(2) = p(i,2)
            y(i) = llgauss(q)
        enddo
        tol = 1e-4
        call amoeba(p,y,3,2,2,tol,llgauss,iter)
*       maybe add restart later
        a = p(1,1)
        b = p(1,2)
*
*       output
*
        print '(a,i5,a)','# Fitted to Gaussian distribution in ',iter
     +        ,' iterations'
        print '(a)','# p(x) = exp(-(x-a)^2/(2*b^2))/(b*sqrt(2*pi)) with'
        print '(a,f16.3)','# a = ',a
        print '(a,f16.3)','# b = ',b
        end
*  #] fitgau:
*  #[ llgauss:
        real function llgauss(p)
*
*       computes the log-likelihood function for a normal distribution
*       with parameters alpha,beta=p(1),p(2) and data in common.
*
        implicit none
*       
        real p(2)
*
        integer i
        real s
*
        integer nmax,ncur
        parameter(nmax=100000)
        real data(nmax),restrain
        common /fitdata1/ data
        common /fitdata2/ restrain,ncur
*       
        llgauss = 0
        do i=1,ncur
            llgauss = llgauss - ((data(i)-p(1))/p(2))**2/2
        enddo
        llgauss = llgauss - ncur*log(p(2))
*       normalization is not 1 in case of cut-offs
        call gausnorm(p(1),p(2),s)
        llgauss = llgauss - ncur*log(s)
*       minimum, not maximum
        llgauss = -llgauss
***        print *,'a,b,llgauss = ',p(1),p(2),llgauss
*
        end
*  #] llgauss:
*  #[ gausnorm:
        subroutine gausnorm(a,b,s)
        implicit none
#include "getopts.inc"
        real a,b,s
        real z1,z2,sqrt2
        real erfcc
        external erfcc
        sqrt2 = sqrt(2.)
        if ( minindx.gt.-1e33 ) then
            z1 = (minindx-a)/b
            if ( maxindx.lt.1e33 ) then
                z2 = (maxindx-a)/b
                s = (erfcc(z1/sqrt2) - erfcc(z2/sqrt2))/2
            else
                s = erfcc(z1/sqrt2)/2
            endif
        else
            if ( maxindx.lt.1e33 ) then
                z2 = (maxindx-a)/b
                s = erfcc(-z2/sqrt2)/2
            else
                s = 1
            endif
        endif
***        print *,'gausnorm: norm = ',a,b,s
        end
*  #] gausnorm:
*  #[ fitgam:
        subroutine fitgam(xx,ntot,mean,sd,a,b)
*
*       a fit a gamma distribution to the data
*
        implicit none
*
        integer ntot
        real xx(ntot),mean,sd,a,b
*
        integer i,nx,iter,init
        real tol,p(3,2),q(2),xmin,fret,y(3)
        save init
*
        integer nmax,ncur
        parameter(nmax=100000)
        real data(nmax),restrain
        common /fitdata1/ data
        common /fitdata2/ restrain,ncur
*
        real llgamma
        external llgamma,dllgamma
*       
        data init /0/
*
*       check input
*       
        xmin = 3e33
        do i=1,ntot
            if ( xx(i).lt.0 ) then
                write(0,*) 'histogram: error: cannot fit Gamma to '//
     +                'negative data: ',xx(i)
                call abort
            endif
            if ( xx(i).ne.0 ) then
                xmin = min(xmin,xx(i))
            endif
        enddo
        xmin = xmin/5
        do i=1,ntot
            if ( xx(i).eq.0 ) then
                if ( init.eq.0 ) then
                    init = 1
                    print '(a,f16.6)','# changed zero to ',xmin
                endif
                xx(i) = xmin
            endif
        enddo
*
*       copy to common for routines llgamma, dllgamma
*       
        ncur = ntot
        do i=1,ncur
            data(i) = xx(i)
        enddo
*       
*       fit, using Numerical Recipes routines
*
***     This frequently crashes
***        q(1) = mean**2/sd**2
***        q(2) = sd**2/mean
***        tol = 1e-4
***        call dfpmin(q,2,tol,iter,fret,llgamma,dllgamma)
***        a = q(1)
***        b = q(2)
***     so try amoeba - slow but sure
        p(1,1) = mean**2/sd**2-0.05
        p(1,2) = sd**2/mean    *0.9
        p(2,1) = p(1,1) +0.1
        p(2,2) = p(1,2)
        p(3,1) = p(1,1)
        p(3,2) = p(1,2) *1.2
        do i=1,3
            q(1) = p(i,1)
            q(2) = p(i,2)
            y(i) = llgamma(q)
        enddo
        tol = 1e-4
        call amoeba(p,y,3,2,2,tol,llgamma,iter)
*       maybe add restart later
        a = p(1,1)
        b = p(1,2)
*
*       output
*
        print '(a,i5,a)','# Fitted to gamma distribution in ',iter
     +        ,' iterations'
        print '(a)','# p(x) = (x/b)^(a-1)*exp(-x/b)/(b*Gamma(a)) with'
        print '(a,f16.3)','# a = ',a
        print '(a,f16.3)','# b = ',b
        end
*  #] fitgam:
*  #[ llgamma:
        real function llgamma(p)
*
*       computes the log-likelihood function for a Gamma distribution
*       with parameters alpha,beta=p(1),p(2) and data in common.
*       
        implicit none
*       
        real p(2)
*
        integer i
*
        integer nmax,ncur
        parameter(nmax=100000)
        real data(nmax),restrain
        common /fitdata1/ data
        common /fitdata2/ restrain,ncur
#include "getopts.inc"
*       
        real gammln,gammp,gammq
        external gammln,gammp,gammq
*
        llgamma = 0
        do i=1,ncur
            llgamma = llgamma + 
     +            (p(1)-1)*log(data(i)/p(2)) - data(i)/p(2) - log(p(2))
     +            - gammln(p(1))
        enddo
*       normalization is not 1 in case of cut-offs
        if ( minindx.gt.-1e33 ) then
            if ( maxindx.lt.1e33 ) then
                llgamma = llgamma - ncur*log(gammp(p(1),maxindx/p(2))
     +                - gammp(p(1),minindx/p(2)))
            else
                llgamma = llgamma - ncur*log(gammq(p(1),minindx/p(2)))
            endif
        else
            if ( maxindx.lt.1e33 ) then
                llgamma = llgamma - ncur*log(gammp(p(1),maxindx/p(2)))
            endif
        endif
*       minimum, not maximum
        llgamma = -llgamma
***        print *,'a,b,llgamma = ',p(1),p(2),llgamma
*
        end
*  #] llgamma:
*  #[ dllgamma:
        subroutine dllgamma(p,dp)
*
*       computes the derivatives of the log-likelihood function for a
*       Gamma distribution with parameters p(1),p(2) and data in common.
*       currently unused.
*       
        implicit none
*       
        real p(2),dp(2)
*
        integer i
        real p1(2),p2(2),d
*
        integer nmax,ncur
        parameter(nmax=100000)
        real data(nmax),restrain
        common /fitdata1/ data
        common /fitdata2/ restrain,ncur
*
        real dgammln,llgamma
        external dgammln,llgamma
*
        dp(1) = 0
        dp(2) = 0
        do i=1,ncur
            dp(1) = dp(1)
     +            + log(data(i)/p(2)) - dgammln(p(1))
            dp(2) = dp(2)
     +            - p(1)/p(2) + data(i)/p(2)**2
        enddo
*       minimum, not maximum
        dp(1) = -dp(1)
        dp(2) = -dp(2)
*
        print *,'dp(1) = ',dp(1)
        p1(1) = p(1) + 1e-3
        p1(2) = p(2)
        p2(1) = p(1) - 1e-3
        p2(2) = p(2)
        d = (llgamma(p1)-llgamma(p2))/2e-3
        print *,'    cmp ',d
        print *,'dp(2) = ',dp(2)
        p1(1) = p(1)
        p1(2) = p(2) + 10
        p2(1) = p(1)
        p2(2) = p(2) - 10
        d = (llgamma(p1)-llgamma(p2))/20
        print *,'    cmp ',d
*
        end
*  #] dllgamma:
*  #[ gamdist:
        real function gamdist(x)
*
*       Gamma distribution, parameters passed in common
*       currently unused
*
        implicit none
        real x
        real alpha,beta
        common /gamcom/ alpha,beta
        real z,y
        real gammln
        external gammln
        if ( x.lt.0 ) then
            write(0,*) 'gamdist: error: cannot evaluate for x&lt;0: ',x
            call abort
        elseif ( x.eq.0 ) then
            if ( alpha.eq.1 ) then
                gamdist = 1/beta
                return
            elseif ( alpha.gt.1 ) then
                gamdist = 0
                return
            else
                write(0,*) 'gamdist: error: infinite for x=0, alpha= '
     +                ,alpha
                call abort
            endif
        endif
        z = x/beta
        y = (alpha-1)*log(z) - z - log(beta) - gammln(alpha)
        gamdist = exp(y)
        end
*  #] gamdist:
*  #[ cumgamm:
        real function cumgamm(x)
*
*       compute the cumulative Gamma probability minus some
*       requested frequency.  Parameters are passed in common
*
        implicit none
        real x
        real pc,ac,bc
        common /ccumgamm/ pc,ac,bc
        real gammp
        external gammp
*
        if ( x.le.0 ) then
            cumgamm = -pc
        else
            cumgamm = gammp(ac,x/bc) - pc
        endif
*
        end
*  #] cumgamm:
*  #[ invcumgamm:
        real function invcumgamm(p,a,b)
*
*       compute the inverse of the cumulative Gamma distribution P(a,x/b)
*       as long as I do not find or make an explicit function just solve
*       the equation.
*
        implicit none
        real p,a,b
        integer i
        real x,x1,x2,tol
        real pc,ac,bc
        common /ccumgamm/ pc,ac,bc
        real cumgamm,zbrent
        external cumgamm,zbrent
*
*       check argument
*
        if ( p.eq.0 ) then
            x = 0
            goto 999
        elseif ( p.lt.0 .or. p.ge.1 ) then
            write(0,*) 'invcumgamm: illegal argument ',p
            x = 3e33
            goto 999
        endif
*
*       parameters for function cumgamm
*
        ac = a
        bc = b
        pc = p
*
*       bracket zero
*
        x1 = 0
        x2 = b
        i = 0
  100   continue
        if ( cumgamm(x2).lt.0 ) then
            x2 = 1.6*x2
            i = i + 1
            if ( i.gt.100 ) then
                write(0,*) 'invcumgamm: error: cannot find root'
                x = 3e33
                goto 999
            endif
            goto 100
        endif
*
*       get root
*
        tol = 1e-5*x2
        x = zbrent(cumgamm,x1,x2,tol)
*
*       finito
*
  999   continue
        invcumgamm = x
        return
        end
*  #] invcumpois:
*  #[ fitgum:
        subroutine fitgum(xx,ntot,mean,sd,a,b,j1,j2,lweb,ntype
     +       ,lchangesign)
*
*       a fit a Gumbel distribution to the data
*
        implicit none
*
        integer nmc
        parameter (nmc=500)
        integer ntot,j1,j2,ntype
        real xx(ntot),mean,sd,a,b
        logical lweb,lchangesign
*
        integer i,j,nx,iter,iens
        real x,t(4),aa(nmc),bb(nmc),tt(nmc,4),a5,b5,a1,b1,t5(4),t1(4)
*
        integer nmax,ncur
        parameter(nmax=100000)
        real data(nmax),restrain
        common /fitdata1/ data
        common /fitdata2/ restrain,ncur
*
        real llgumbel
        real *8 ranf
        external llgumbel,ranf
*
*       copy to common for routine llgumbel
*
        ncur = ntot
        do i=1,ncur
            data(i) = xx(i)
        enddo
        b = sd*sqrt(6.)/(4*atan(1.))
        a = mean - 0.57721*b
        call fit1gum(a,b,iter)
        do i=1,4
            x = i + log10(real(j2-j1+1))
            if ( x.le.8 ) then
                t(i) = a - b*log(-log(1-dble(10)**(dble(-x))))
            else
                t(i) = a + b*x*log(10.)
            endif
            if ( lchangesign ) then
                t(i) = -t(i)
            endif
        enddo
*       
*       bootstrap to find error estimates
*       
        if ( .not.lweb ) print '(a,i6,a)','# Doing a ',nmc
     +        ,'-member bootstrap to obtain error estimates'
        do iens=1,nmc
            do i=1,ncur
                j = 1+int(ncur*ranf(i))
                if ( j.lt.1 .or. j.gt.ncur ) then
                    write(0,*) 'fitgum: error: j = ',j
                    call abort
                endif
                data(i) = xx(j)
            enddo
            aa(iens) = a
            bb(iens) = b
            call fit1gum(aa(iens),bb(iens),iter)
            if ( lchangesign ) then
                aa(iens) = -aa(iens)
                bb(iens) = -bb(iens)
            endif
            do i=1,4
                x = i + log10(real(j2-j1+1))
                if ( x.le.8 ) then
                    tt(iens,i) = aa(iens)
     +                    - bb(iens)*log(-log(1-dble(10)**(dble(-x))))
                else
                    tt(iens,i) = aa(iens) + bb(iens)*x*log(10.)
                endif
            enddo
        enddo
        call getcut(a5,95.,nmc,aa)
        call getcut(b5,95.,nmc,bb)
        call getcut(a1,99.,nmc,aa)
        call getcut(b1,99.,nmc,bb)
        do i=1,4
            call getcut(t5(i),95.,nmc,tt(1,i))
            call getcut(t1(i),99.,nmc,tt(1,i))
        enddo
        if ( lchangesign ) then
            a = -a
            b = -b
        endif
*
*       output
*
        if ( lweb ) then
            print '(a)','# <tr><td colspan="3">Fitted to '//
     +           'Gumbel distribution P(x) = exp(-exp(-(x-a)/b))'//
     +           '</td></tr>'
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td>a:</td><td>'
     +           ,a,'</td><td> 5%&gt;',a5,', 1%&gt;',a1,'</td></tr>'
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td>b:</td><td>'
     +           ,b,'</td><td> 5%&gt;',b5,', 1%&gt;',b1,'</td></tr>'
            do i=1,4
                print '(a,i5,a,f16.3,a,f16.3,a,f16.3,a)'
     +               ,'# <tr><td>return period ',10**i,' yr</td><td>'
     +               ,t(i),'</td><td> 5%&gt;',t5(i),', 1%&gt;',t1(i)
     +               ,'</td></tr>'
                call print3untransf(t(i),t5(i),t1(i))
            enddo
        else
            print '(a,i5,a)','# Fitted to Gumbel distribution in ',iter
     +        ,' iterations'
            print '(a)','# P(x) = exp(-exp(-(x-a)/b)) with'
            print '(a,f16.3,a,f16.3,a,f16.3)','# a = ',a,' 5%>',a5
     +           ,' 1%>',a1
            print '(a,f16.3,a,f16.3,a,f16.3)','# b = ',b,' 5%>',b5
     +           ,' 1%>',b1
            do i=1,4
                print '(a,i5,a,f16.3,a,f16.3,a,f16.3)'
     +               ,'# value for return period ',10**i,' year: ',t(i)
     +               ,' 5%>',t5(i),' 1%>',t1(i)
                call print3untransf(t(i),t5(i),t1(i))
            enddo
        endif
        end
*  #] fitgum:
*  #[ fit1gum:
        subroutine fit1gum(a,b,iter)
        implicit none
        integer iter
        real a,b
        integer i
        real q(2),p(3,2),y(3),tol
        real llgumbel
        external llgumbel
*
*       fit, using Numerical Recipes routines
*       
        q(1) = a
        q(2) = b
        p(1,1) = q(1) *0.9
        p(1,2) = q(2) *0.9
        p(2,1) = p(1,1) *1.2
        p(2,2) = p(1,2)
        p(3,1) = p(1,1)
        p(3,2) = p(1,2) *1.2
        do i=1,3
            q(1) = p(i,1)
            q(2) = p(i,2)
            y(i) = llgumbel(q)
        enddo
        tol = 1e-4
        call amoeba(p,y,3,2,2,tol,llgumbel,iter)
*       maybe add restart later
        a = p(1,1)
        b = abs(p(1,2))
        end
*  #] fit1gum:
*  #[ llgumbel:
        real function llgumbel(p)
*
*       computes the log-likelihood function for a Gumbel distribution
*       with parameters a,b=p(1),p(2) and data in common.
*
        implicit none
*       
        real p(2)
*
        integer i
        real z,s
*
        integer nmax,ncur
        parameter(nmax=100000)
        real data(nmax),restrain
        common /fitdata1/ data
        common /fitdata2/ restrain,ncur
*
        llgumbel = 0
        do i=1,ncur
            z = (data(i)-p(1))/abs(p(2))
            llgumbel = llgumbel - exp(-z) - z
        enddo
        llgumbel = llgumbel - ncur*log(abs(p(2)))
*       normalization is not 1 in case of cut-offs
        call gumbnorm(p(1),abs(p(2)),s)
        llgumbel = llgumbel - ncur*log(s)
*       minimum, not maximum
        llgumbel = -llgumbel
***        print *,'a,b,llgumbel = ',p(1),p(2),llgumbel
*
        end
*  #] llgumbel:
*  #[ gumbnorm:
        subroutine gumbnorm(a,b,s)
        implicit none
#include "getopts.inc"
        real a,b,s
        real z1,z2
        if ( minindx.gt.-1e33 ) then
            z1 = (minindx-a)/b
            if ( maxindx.lt.1e33 ) then
                z2 = (maxindx-a)/b
                s = exp(-exp(-z2)) - exp(-exp(-z1))
            else
                s = 1 - exp(-exp(-z1))
            endif
        else
            if ( maxindx.lt.1e33 ) then
                z2 = (maxindx-a)/b
                s = exp(-exp(-z2))
            else
                s = 1
            endif
        endif
***        print *,'gumbnorm: norm = ',a,b,s
        end
*  #] gumbnorm:
*  #[ fitgev:
        subroutine fitgev(xx,ntot,mean,sd,a,b,xi,j1,j2,lweb,ntype
     +       ,lchangesign)
*
*       a fit a GEV distribution to the data
*
        implicit none
*       
        integer nmc
        parameter(nmc=500)
        integer ntot,j1,j2,ntype
        real xx(ntot),mean,sd,a,b
        logical lweb,lchangesign
*
        integer i,j,nx,iter,iens
        real x,xi,t(4),aa(nmc),bb(nmc),xixi(nmc),tt(nmc,4),a25,a975,b25
     +        ,b975,xi25,xi975,t5(4),t1(4),t25(4),t975(4),da,db,dxi,f
        character lgt*4
*
        integer nmax,ncur
        parameter(nmax=100000)
        real data(nmax),restrain
        common /fitdata1/ data
        common /fitdata2/ restrain,ncur
*
        real llgev
        real *8 ranf
        external llgev,ranf
*
*       copy to common for routine llgev
*
        ncur = ntot
        do i=1,ncur
            data(i) = xx(i)
        enddo
        b = sd*sqrt(6.)/(4*atan(1.))
        a = mean - 0.57721*b
        xi = 0
        call fit1gev(a,b,xi,iter)
        if ( xi.eq.0 ) then
            do i=1,4
                x = i + log10(real(j2-j1+1))
                if ( x.le.8 ) then
                    t(i) = a - b*log(-log(1-dble(10)**(dble(-x))))
                else
                    t(i) = a + b*x*log(10.)
                endif
            enddo
        else
            do i=1,4
                x = i + log10(real(j2-j1+1))
                if ( x.le.8 ) then
                    t(i) = a - b/xi*(1-(-log(1-dble(10)**(dble(-x))))
     +                    **(-xi))
                else
                    t(i) = a - b/xi*(1-10.**(x*xi))
                endif
                if ( lchangesign ) t(i) = -t(i)
            enddo
        endif
*       
*       bootstrap to find error estimates
*       
        if ( .not.lweb ) print '(a,i6,a)','# doing a ',nmc
     +        ,'-member bootstrap to obtain error estimates'
        do iens=1,nmc
            if ( mod(iens,100).eq.0 ) print '(a,i6)','# ',iens
            do i=1,ncur
                j = 1+int(ncur*ranf(i))
                if ( j.lt.1 .or. j.gt.ncur ) then
                    write(0,*) 'fitgev: error: j = ',j
                    call abort
                endif
                data(i) = xx(j)
            enddo
            aa(iens) = a
            bb(iens) = b
            call fit1gev(aa(iens),bb(iens),xixi(iens),iter)
            if ( lchangesign ) then
                aa(iens) = -aa(iens)
                bb(iens) = -bb(iens)
            endif
            do i=1,4
                x = i + log10(real(j2-j1+1))
                if ( xi.eq.0 ) then
                    if ( x.le.8 ) then
                        tt(iens,i) = aa(iens) - bb(iens)*
     +                        log(-log(1-dble(10)**(dble(-x))))
                    else
                        tt(iens,i) = aa(iens) + bb(iens)*x*log(10.)
                    endif
                else
                    if ( x.le.8 ) then
                        tt(iens,i) = aa(iens) - bb(iens)/xi*
     +                        (1-(-log(1-dble(10)**(dble(-x))))**(-xi))
                    else
                        tt(iens,i) = aa(iens)
     +                        - bb(iens)/xi*(1-10.**(x*xi))
                    endif
                endif
            enddo
        enddo
        if ( lchangesign ) then
            a = -a
            b = -b
        endif
        call getcut( a25, 2.5,nmc,aa)
        call getcut(a975,97.5,nmc,aa)
        call getcut( b25, 2.5,nmc,bb)
        call getcut(b975,97.5,nmc,bb)
        call getcut( xi25, 2.5,nmc,xixi)
        call getcut(xi975,97.5,nmc,xixi)
        do i=1,4
            if ( lchangesign ) then
                lgt = '&lt;'
                call getcut(t5(i),5.,nmc,tt(1,i))
                call getcut(t1(i),1.,nmc,tt(1,i))
            else
                lgt = '&gt;'
                call getcut(t5(i),95.,nmc,tt(1,i))
                call getcut(t1(i),99.,nmc,tt(1,i))
            endif
            call getcut(t25(i),2.5,nmc,tt(1,i))
            call getcut(t975(i),97.5,nmc,tt(1,i))
        enddo
*
*       output
*
        if ( lweb ) then
            print '(a)','# <tr><td colspan="3">Fitted to GEV '//
     +           'distribution P(x) = exp(-(1+&xi;*(x-a)/b)**(-1/&xi;))'
     +           //'</td></tr>'
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td>a:</td><td>'
     +           ,a,'</td><td>',a25,'...',a975,'</td></tr>'
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td>b:</td><td>'
     +           ,b,'</td><td>',b25,'...',b975,'</td></tr>'
            print '(a,f16.3,a,f16.3,a,f16.3,a)'
     +           ,'# <tr><td>&xi;:</td><td>',xi,'</td><td>',xi25,'...'
     +           ,xi975,'</td></tr>'
            do i=1,4
                print '(a,i5,a,f16.3,2a,f16.3,2a,f16.3,a)'
     +               ,'# <tr><td>return period ',10**i,' yr</td><td>'
     +               ,t(i),'</td><td> 5%',lgt,t5(i),', 1%',lgt,t1(i)
     +               ,'</td></tr>'
                call print3untransf(t(i),t5(i),t1(i))
            enddo
            if ( ntype.eq.2 ) then ! Gumbel plot
                do i=1,4
                    f = 1 - 10.**(-i)
                    print '(i8,3g22.6)',0,-log(-log(f)),-999.9,t25(i)
                enddo
                print '(a)'
                do i=1,4
                    f = 1 - 10.**(-i)
                    print '(i8,3g22.6)',0,-log(-log(f)),-999.9,t975(i)
                enddo
                print '(a)'
            endif
        else
            print '(a,i5,a)','# Fitted to GEV distribution in ',iter
     +           ,' iterations'
            print '(a)','# P(x) = exp(-(1+xi*(x-a)/b)**(-1/xi)) with'
***         print '(a,f16.3,a,f16.3,a,f16.3)','# a = ',a,' 5%&gt;',a5
***         +        ,' 1%&gt;',a1
***         print '(a,f16.3,a,f16.3,a,f16.3)','# b = ',b,' 5%&gt;',b5
***     +        ,' 1%&gt;',b1
***         print '(a,f16.3,a,f16.3,a,f16.3)','# xi= ',xi,' 95%',xi25
***         +        ,'&lt;xi&lt;',xi975
            print '(a,f16.3,a,f16.3,a,f16.3)','# a = ',a,' \\pm ',a975
     +           -a25
            print '(a,f16.3,a,f16.3,a,f16.3)','# b = ',b,' \\pm ',b975
     +           -b25
            print '(a,f16.3,a,f16.3,a,f16.3)','# xi= ',xi,' \\pm ',xi975
     +           -xi25
            do i=1,4
                print '(a,i5,a,f16.3,a,f16.3,a,f16.3)'
     +               ,'# value for return period ',10**i,' year: ',t(i)
     +               ,' 5%>',t5(i),' 1%> ',t1(i)
                call print3untransf(t(i),t5(i),t1(i))
            enddo
        endif
        end
*  #] fitgev:
*  #[ fit1gev:
        subroutine fit1gev(a,b,xi,iter)
        implicit none
        integer iter
        real a,b,xi
        integer i
        real q(3),p(4,3),y(4),tol
        real llgev
        external llgev
*
*       fit, using Numerical Recipes routines
*       
        q(1) = a
        q(2) = b
        q(3) = xi
        p(1,1) = q(1) *0.9
        p(1,2) = q(2) *0.9
        p(1,3) = q(3) *0.9
        p(2,1) = p(1,1) *1.2
        p(2,2) = p(1,2)
        p(2,3) = p(1,3)
        p(3,1) = p(1,1)
        p(3,2) = p(1,2) *1.2
        p(3,3) = p(1,3)
        p(4,1) = p(1,1)
        p(4,2) = p(1,2)
        p(4,3) = p(1,3) *1.2 + 0.1
        do i=1,4
            q(1) = p(i,1)
            q(2) = p(i,2)
            q(3) = p(i,3)
            y(i) = llgev(q)
        enddo
        tol = 1e-4
        call amoeba(p,y,4,3,3,tol,llgev,iter)
*       maybe add restart later
        a = p(1,1)
        b = abs(p(1,2))
        xi = p(1,3)
        end
*  #] fit1gev:
*  #[ llgev:
        real function llgev(p)
*
*       computes the log-likelihood function for a GEV distribution
*       with parameters a,b=p(1),p(2) and data in common.
*
        implicit none
*       
        real p(3)
*
        integer i
        real z,xi,s
*
        integer nmax,ncur
        parameter(nmax=100000)
        real data(nmax),restrain
        common /fitdata1/ data
        common /fitdata2/ restrain,ncur
*
        llgev = 0
        if ( p(2).eq.0 ) then
            write(0,*) 'llgev: cannot handle b=0'
            call abort
        endif
        do i=1,ncur
            z = (data(i)-p(1))/abs(p(2))
            xi = p(3)
            if ( abs(xi).lt.1e-4 ) then
                llgev = llgev - exp(-z+xi*z**2/2) - z*(1+xi-xi*z/2)
            else
                if ( 1+xi*z.le.0 ) then
***                    write(0,*) 'GEV undefined',(1+xi*z)
                    llgev = llgev - 3e33
                else
                    llgev = llgev - (1+1/xi)*log(1+xi*z)
     +                    - (1+xi*z)**(-1/xi)
                endif
            endif
        enddo
        llgev = llgev - ncur*log(abs(p(2)))
*       normalization is not 1 in case of cut-offs
        call gevnorm(p(1),abs(p(2)),p(3),s)
        llgev = llgev - ncur*log(s)
*       minimum, not maximum
        llgev = -llgev
***        print *,'ncur,a,b,xi,llgev = ',ncur,p(1),p(2),p(3),llgev
*
        end
*  #] llgev:
*  #[ gevnorm:
        subroutine gevnorm(a,b,xi,s)
        implicit none
#include "getopts.inc"
        real a,b,xi,s
        real z1,z2
        if ( minindx.gt.-1e33 ) then
            z1 = (minindx-a)/b
            if ( maxindx.lt.1e33 ) then
                z2 = (maxindx-a)/b
                if ( xi.eq.0 ) then
                    s = exp(-exp(-z2)) - exp(-exp(-z1))
                else
                    s = exp(-(1+xi*z2)**(-1/xi)) - 
     +                    exp(-(1+xi*z1)**(-1/xi))
                endif
            else
                if ( xi.eq.0 ) then
                    s = 1 - exp(-exp(-z1))
                else
                    s = 1 - exp(-(1+xi*z1)**(-1/xi))
                endif
            endif
        else
            if ( maxindx.lt.1e33 ) then
                z2 = (maxindx-a)/b
                if ( xi.eq.0 ) then
                    s = exp(-exp(-z2))
                else
                    s = exp(-(1+xi*z2)**(-1/xi))
                endif
            else
                s = 1
            endif
        endif
***        print *,'gevnorm: norm = ',a,b,s
        end
*  #] gevnorm:
*  #[ printval:
      subroutine printval(string,val,plot,iuplot,lweb)
      implicit none
      integer iuplot
      real val
      logical plot,lweb
      character string*(*)
      integer i
      if ( lweb ) then
          print '(3a,g16.4,10a)','# <tr><td>',string(3:),'</td><td>',val
     +         ,'</td><td>&nbsp;</td></tr>'
      else
          print '(a,6g16.4)',string,val
      endif
      if ( plot ) write(iuplot,'(6g16.4)') val
      end
*  #] printval:
*  #[ printvalerr:
      subroutine printvalerr(string,array,mboot,plot,iuplot,lweb)
      implicit none
      integer iuplot,mboot
      real array(0:mboot)
      logical plot,lweb
      character string*(*)
      if ( mboot.eq.0 ) then
          call printval(string,array,plot,iuplot,lweb)
      else
          if ( lweb ) then
              print '(3a,7(g16.6,a))','# <tr><td>',string(3:),
     +             '</td><td>',array(0),'&plusmn;',
     +             (array(nint(0.975*mboot))-array(nint(0.025*mboot)))/2
     +             ,'</td><td>',
     +             array(nint(0.025*mboot)),'...',
     +             array(nint(0.975*mboot)),'</td></tr>'
          else
              print '(a,7g16.6)',string,array(0),
     +             array(nint(0.025*mboot)),
     +             array(nint(0.16*mboot)),
     +             array(nint(0.50*mboot)),
     +             array(nint(0.84*mboot)),
     +             array(nint(0.975*mboot)),
     +             (array(nint(0.975*mboot))-array(nint(0.025*mboot)))/2
          endif
          if ( plot ) write(iuplot,'(7g16.6)') array(0),
     +         array(nint(0.025*mboot)),
     +         array(nint(0.16*mboot)),
     +         array(nint(0.50*mboot)),
     +         array(nint(0.84*mboot)),
     +         array(nint(0.975*mboot))
      endif
      end
*  #] printvalerr:
*  #[ printuntransf:
        subroutine printuntransf(xin)
        implicit none
#include "getopts.inc"
        real xin
        real x
        x = xin
        if ( logscale .or. sqrtscale .or. squarescale ) then
            if ( logscale ) x = 10.**(x)
            if ( sqrtscale ) x = x**2
            if ( squarescale ) x = sqrt(x)
            if ( lweb ) then
                print '(a,f16.2)','# w/o transf.',x
            else
                print '(a,f16.2,a)','# <tr><td>w/o transf.</td><td>',x
     +               ,'</td><td>&nbsp;</td></tr>'
            endif
        endif
        end
*  #] printuntransf:
*  #[ print3untransf:
        subroutine print3untransf(x1in,x2in,x3in)
        implicit none
#include "getopts.inc"
        real x1in,x2in,x3in
        real x(3)
        integer i
        if ( logscale .or. sqrtscale .or. squarescale ) then
            do i=1,3
                if ( i.eq.1 ) then
                    x(i) = x1in
                elseif ( i.eq.2 ) then
                    x(i) = x2in
                else
                    x(i) = x3in
                endif
                if ( logscale ) x(i) = 10.**(x(i))
                if ( sqrtscale ) x(i) = x(i)**2
                if ( squarescale ) x(i) = sqrt(x(i))
            enddo
            print '(a,f16.2,a,f16.3,a,f16.3)','# w/o transf.',x(1)
     +            ,' 5%&gt;',x(2),' 1%&gt;',x(3)
        endif
        end
*  #] print3untransf:
