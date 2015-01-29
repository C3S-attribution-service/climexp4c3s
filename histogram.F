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
        parameter (maxbin=1000,maxldat=10000000,nboot=800)

        integer i,ii,j,jj,j1,j2,n,m,nperyear,nbin,month,yr,nn(maxbin)
     +       ,Nless,nfit,ntot,ntype,ier,iens,mens1,mens,l,nnn
     +       ,iboot,mboot,ndecor,nfitted,it,tmax,off,yrstart,yrstop
     +       ,dy,mo,iret,jmax,iensmax
        integer,allocatable :: yrs(:)
        logical lastvalid,lprinted,lexist
        logical,allocatable :: lf(:)
        real,allocatable :: data(:,:,:),xx(:),xs(:)
        real mindat,maxdat,frac
     +       ,d,a,b,c,da,db,dc,xi,dxi,yy(maxbin)
     +       ,x,x1,x2,s,sqrt2,snorm,f,ff,offset,xn(maxbin),xy(maxbin),df
     +       ,chsq,prob,lon,lat,elev,threshold,z,tx,tx25
     +       ,tx975,t(10),t25(10),t975(10),xmax
        real mean(0:nboot),sd(0:nboot),skew(0:nboot)
        character file*255,string*100,var*20,units*20,prog*100,
     +       number*100,line*255,command*255,extraargs*100,
     +       extraargs_*100,ids(0:nensmax)*30

        integer iargc
        real gammln,erf,erfcc,gammq,gammp,invcumpois,invcumgamm
        real *8 ranf
        external gammln,erf,erfcc,gammp,gammq,invcumpois
     +        ,invcumgamm,ranf
*
        call killfile(string,line,file,0)
        allocate(yrs(0:maxldat))
        allocate(lf(maxldat))
        allocate(xx(maxldat))
        allocate(xs(maxldat))
        lwrite = .false.
        n = iargc()
        if ( n.lt.2 ) then
            print *,'usage: histogram {timeseries|file list prog} nbins'
     +           //' hist hist|qq|gumbel|log|sqrtlog'
     +           //' [fit poisson|gauss|gamma|gumbel|gev|gpd] [options]'
            stop
        endif
!
!       getopts comes way too late for these options...
        do i=4,iargc()
            call getarg(i,string)
            if ( string.eq.'debug' .or. string.eq.'lwrite' ) then
                lwrite = .true.
            end if
            if ( string(1:9).eq.'standardu' ) then
                lstandardunits = .true.
            end if
        end do
        call getarg(1,file)
        allocate(data(npermax,yrbeg:yrend,0:nensmax))
        if ( file.ne.'file' ) then
*           simple data file (possibly an ensemble)
            off = 0
            print '(2a)','# file ',trim(file)
            call readensseries(file,data,npermax,yrbeg,yrend,nensmax
     +           ,nperyear,mens1,mens,var,units,lstandardunits,lwrite)
        else
            off = 2
            call readsetseries(data,ids,npermax,yrbeg,yrend,nensmax,
     +           nperyear,mens1,mens,var,units,lstandardunits,lwrite)
        endif
        call getarg(2+off,string)
        read(string,*) nbin
        if ( nbin.gt.maxbin ) then
            write(0,*) 'histogram: error: increase maxbin'
            write(*,*) '# histogram error: increase maxbin'
            call abort
        endif
        i = 3 + off
        nfit = 0
        ntype = 0
        yrs(0) = 0
        if ( n.ge.i ) then
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
                elseif ( string(1:3).eq.'gpd' ) then
                    nfit = 6
                else
                    write(0,*) 'histogram: error: unknown distribution '
     +                    ,string
                    write(*,*)
     +                   '# histogram error: unknown distribution '
     +                   ,string
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
                elseif ( string(1:3).eq.'log' ) then
                    ntype = 3
                elseif ( string(1:7).eq.'sqrtlog' ) then
                    ntype = 4
                else
                    write(0,*) 'histogram: error: unknown plot ',string
                    write(*,*) '# histogram error: unknown plot ',string
                    call abort
                endif
                goto 10
            endif
            call getopts(i,n,nperyear,yrbeg,yrend,.true.,mens1,mens)
            if ( mens.gt.0 ) then
                print '(a,i3,a,i3)','# using ensemble members ',nens1,
     +                ' to ',nens2
            endif
        endif
        if ( plot .and. .not. lbootstrap ) then
            mboot = 0
        else
            mboot = nboot
        endif
        yrstart = yr2
        yrstop  = yr1
*
*       sum series
*
        do iens=nens1,nens2
            if ( lsum.gt.1 ) then
                call sumit(data(1,yrbeg,iens),npermax,nperyear
     +               ,yrbeg,yrend,lsum,oper)
            endif
            if ( logscale ) then
                if ( iens.eq.nens1 ) print '(a)'
     +               ,'# taking logarithm'
                call takelog(data(1,yrbeg,iens),npermax
     +               ,nperyear,yrbeg,yrend)
            endif
            if ( sqrtscale ) then
                if ( iens.eq.nens1 ) print '(a)','# taking sqrt'
                call takesqrt(data(1,yrbeg,iens),npermax
     +               ,nperyear,yrbeg,yrend)
            endif
            if ( squarescale ) then
                if ( iens.eq.nens1 ) print '(a)','# taking square'
                call takesquare(data(1,yrbeg,iens),npermax
     +               ,nperyear,yrbeg,yrend)
            endif
            if ( twothirdscale ) then
                if ( iens.eq.nens1 ) print '(a)',
     +               '# taking power twothird'
                call taketwothird(data(1,yrbeg,iens),npermax
     +               ,nperyear,yrbeg,yrend)
            endif
        enddo
*
*       apply minindx, maxindx cuts;
*       in case of a Poisson distribution round to unambiguous
*       numbers
*
        if ( m1.ne.m2 ) then
            write(0,*) 'please specify only one starting month'
            call abort
        endif
        call getj1j2(j1,j2,m1,nperyear,.false.)
        if ( pminindx.ge.0 ) then
            call getenscutoff(minindx,pminindx,data,npermax,nperyear
     +           ,yrbeg,yrend,nensmax,nens1,nens2,yr1,yr2,j1,j2,0)
            if ( lwrite ) print *,'pminindx,minindx = ',pminindx
     +           ,minindx
        endif
        if ( pmaxindx.ge.0 ) then
            call getenscutoff(maxindx,pmaxindx,data,npermax,nperyear
     +           ,yrbeg,yrend,nensmax,nens1,nens2,yr1,yr2,j1,j2,0)
            if ( lwrite ) print *,'pmaxindx,maxindx = ',pmaxindx
     +           ,maxindx
        endif
        if ( nfit.eq.6 .and. pmindata.ge.0 ) then
            if ( lchangesign ) then
                call getenscutoff(mindata,100-pmindata,data,npermax
     +               ,nperyear,yrbeg,yrend,nensmax,nens1,nens2,yr1
     +               ,yr2,j1,j2,0)
            else
                call getenscutoff(mindata,pmindata,data,npermax
     +               ,nperyear,yrbeg,yrend,nensmax,nens1,nens2,yr1
     +               ,yr2,j1,j2,0)
            endif
            if ( lwrite ) print *,'pmindata,mindata = ',pmindata
     +           ,mindata
        endif
        if ( nfit.eq.1 ) then
            if ( minindx.ge.0 ) minindx = 0.5 + int(minindx+1e-5)
            if ( maxindx.lt.2.**31 ) maxindx = 0.5 + int(maxindx -
     +           1e-5)
        endif
        do iens=nens1,nens2
            if ( maxindx.lt.1e33 .or. minindx.gt.-1e33 ) then
                do i=yr1,yr2
                    do j=1,nperyear
                        if (  data(j,i,iens).gt.maxindx .or.
     +                       data(j,i,iens).lt.minindx ) then
                            data(j,i,iens) = 3e33
                        endif
                    enddo
                enddo
            endif
*
*           minus sign for lower tails
*
            if ( lchangesign ) then
                if ( iens.eq.nens1 ) print '(a)','# changing sign'
                call changesign(data(1,yrbeg,iens),npermax
     +               ,nperyear,yrbeg,yrend)
            endif
*           
*           detrend data
*           
            if ( ldetrend ) then
                if ( iens.eq.nens1 ) print '(a,i4)'
     +               ,'# detrending data'
                call detrend(data(1,yrbeg,iens),npermax
     +               ,nperyear,yrbeg,yrend,yr1,yr2,m1,m2,lsel)
            endif
*           
*           differentiate data
*           
            if ( ndiff.ne.0 ) then
                if ( iens.eq.nens1 ) print '(a,i4)'
     +               ,'# taking differences/averaging ',ndiff
                call diffit(data(1,yrbeg,iens),npermax,nperyear
     +               ,yrbeg,yrend,ndiff)
            endif
*
*           anomalies
*
            if ( anom ) then
                if ( iens.eq.nens1 ) print '(a)'
     +               ,'# taking anomalies'
                call anomal(data(1,yrbeg,iens),npermax,nperyear
     +               ,yrbeg,yrend,yr1,yr2)
            endif
        enddo
!
!       anomalies wrt ensemble mean
!
        if ( lensanom .and. nens1.ne.nens2 ) then
            if ( lwrite ) print '(a)'
     +           ,'# taking anomalies wrt ensemble mean'
            call anomalensemble(data,npermax,nperyear,yrbeg,
     +           yrend,yr1,yr2,nens1,nens2)
        endif
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
                        if ( i.ge.yr1 .and. i.le.yr2) then
                            if ( i.eq.yr2a .and. data(j,i,iens).lt.1e33
     +                           ) then
                                if ( xmax.lt.data(j,i,iens) ) then
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
                    if ( j1.ne.j2 .and. j2-j1+1.ne.nperyear ) then
                        lastvalid = .false.
                    endif
                    do jj=j1,j2
                        j = jj
                        call normon(j,yr,i,nperyear)
                        if ( i.ge.yr1 .and. i.le.yr2) then
                            if ( i.eq.yr2a .and. j.eq.jmax ) then
                                ! set xyear to the largest value
                                if ( iens.eq.iensmax ) then
                                    xyear = data(j,i,iens)
                                end if
                                ! and set all ensemble members to undef
                                ! to avoid possibly dependent data
                                data(j,i,iens) = 3e33
                            endif
                            if ( data(j,i,iens).lt.1e33 ) then
                                ntot = ntot + 1
                                if ( ntot.gt.maxldat ) then
                                    print *,'# histogram error '//
     +                                   'increase maxldat'
                                    write(0,*) 'histogram: error: '
     +                                   //'increase maxldat'
                                    call abort
                                endif
                                xx(ntot) = data(j,i,iens)
                                call getdymo(dy,mo,j,nperyear)
                                yrs(ntot) = 10000*i + 100*mo + dy
                                if ( nperyear.gt.1000 ) then
                                    yrs(ntot) = 100*yrs(ntot)
     +                                   + mod(ntot,nint(nperyear/366.))
                                end if
                                if ( nens1.ne.nens2 ) then
                                    yrs(ntot) = 100*yrs(ntot) + iens
                                endif
                                offset = offset + xx(ntot)
                                yrstart = min(yrstart,i)
                                yrstop  = max(yrstop,i)
                                lf(ntot) = .not.lastvalid
                                lastvalid = .true.
                                if ( lwrite .and. lf(ntot) ) print *
     +                               ,'boundary at ',ntot,j,i,iens
                            else
                                lastvalid = .false.
                            endif
                        else
                            lastvalid = .false.
                        endif
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
                call keepalive1('Bootstrapping',iboot,mboot)
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
                if ( sd(iboot).ge.0 ) then
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
            call print_bootstrap_message(ndecor,j1,j2)
            if ( lweb ) then
                if ( namestring.ne.' ' ) then ! currently not used
                    print '(4a)','# <tr><th colspan=3>',trim(namestring)
     +                   ,'</th></tr>'
                endif
                print '(9a)'
     +               ,'# <tr><th>parameter</th><th>value&plusmn;2&sigma'
     +               ,'</th><th>95% CI</th></tr>'
                print '(a,i9,a)','# <tr><td>N:</td><td>',ntot,
     +               '</td><td>&nbsp;</td></tr>'
            else
                print '(8a)',       '#        ',
     +               '                ',
     +               '            2.5%',
     +               '             16%',
     +               '             50%',
     +               '             84%',
     +               '           97.5%',
     +               '          +/-95%'
                print '(a,i9)','# N:         ',ntot
            endif
            call printvalerr('# mean:      ',mean,mboot,plot,11,lweb
     +           ,lchangesign)
            call printuntransf(mean(0))
            call printvalerr('# s.d.(n):   ',sd,mboot,plot,11,lweb,
     +           .false.)
            call printval   ('# s.d.(n-1): ',sd(0)*sqrt(real(ntot)
     +            /real(ntot-1)),plot,11,lweb,.false.)
            call printvalerr('# skew:      ',skew,mboot,plot,11,lweb
     +           ,lchangesign)
            call printval   ('# min:       ',mindat,plot,11,lweb
     +           ,lchangesign)
            call printuntransf(mindat)
            call printval   ('# max:       ',maxdat,plot,11,lweb
     +           ,lchangesign)
            call printuntransf(maxdat)
*
*           fit to distribution
*       
            xs(1:ntot) = xx(1:ntot)
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
                call fitgau(xx,ntot,mean(0),sd(0),a,b,minindx,maxindx,
     +               ntype,j1,j2,yr2a,xyear,t,t25,t975,tx,tx25,tx975,
     +               .true.,.true.,lweb,lchangesign,lwrite)
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
                if ( nperyear.gt.12 .and. oper.ne.'a' .and. oper.ne.'i'
     +               ) then
                    write(0,*) 'Warning: data do not appear to be block'
     +                   //' maxima or minima.<p>'
                endif
                call fitgum(xx,ntot,mean(0),sd(0),a,b,j1,j2,lweb,ntype
     +               ,lchangesign,yr2a,xyear,t,t25,t975,tx
     +               ,tx25,tx975,.true.,.true.,lwrite)
                call gumbnorm(a,b,snorm)
                nfitted = 3
            elseif ( nfit.eq.5 ) then
*               GEV distribution
                if ( nperyear.gt.12 .and. oper.ne.'a' .and. oper.ne.'i'
     +               ) then
                    write(0,*) 'Warning: data do not appear to be block'
     +                   //' maxima or minima.<p>'
                endif
                call fitgev(xx,ntot,mean(0),sd(0),a,b,xi,j1,j2,lweb
     +               ,ntype,lchangesign,yr2a,xyear,t,t25,t975,tx
     +               ,tx25,tx975,restrain,.true.,.true.,lwrite)
                call gevnorm(a,b,xi,snorm)
                nfitted = 4
            elseif ( nfit.eq.6 ) then
*               GPD distribution
                call fitgpd(xx,ntot,mean(0),sd(0),b,xi,j1,j2,lweb,
     +               ntype,lchangesign,pmindata,threshold,yr2a,xyear,
     +               t,t25,t975,tx,tx25,tx975,restrain,.true.,.true.,
     +               lwrite)
***                print *,'mindata,threshold = ',mindata,threshold
                snorm = 1
                nfitted = 3
            else
                write(0,*) 'histogram: error: unknown distribution '
     +               ,nfit
                call abort
            endif
            if ( ntype.eq.0 ) then
*
*               adjust binsize to nice numbers (heuristics)
*
*               if the minimum is close to zero make exactly zero
                if ( mindat.gt.0 .and. mindat.lt.maxdat/3 ) then
                    if ( lwrite ) print *,'adjusted mindat from ',mindat
     +                   ,' to 0'
                    mindat = 0
                end if
                if ( maxdat.lt.0 .and. maxdat.gt.mindat/3 ) then
                    if ( lwrite ) print *,'adjusted maxdat from ',maxdat
     +                   ,' to 0'
                    maxdat = 0
                end if
*               if the minimum and maxmimum are similar in absolute value 
*		make them equal
		if ( mindat.lt.0 .and. maxdat.gt.0 .and. 
     +			abs(mindat+maxdat)/2.lt.maxdat/3 ) then
		    if ( mindat.lt.-maxdat ) then
                        if ( lwrite ) print *,'adjusted maxdat from '
     +                       ,maxdat,' to -mindat ',-mindat
			maxdat = -mindat
		    else
                        if ( lwrite ) print *,'adjusted mindat from '
     +                       ,mindat,' to -maxdat ',-maxdat
			mindat = -maxdat
		    endif
		endif
!               give the user a chnce to specify the end of the scale
                if ( maxdata.lt.1e33 ) then
                    if ( lwrite ) print *,'adjusted maxdat from ',maxdat
     +                   ,' to user-defined value ',maxdata
                    maxdat = maxdata
                    mindat = 0
                    d = maxdat/nbin
                else
                    d = (maxdat-mindat)/(nbin-1)
                    s = d
                    if ( d.ge.1 ) then
                        i = int(log10(d))
                    elseif ( d.gt.0 ) then
                        i = int(log10(d)) - 1
                    else
                        i = 1
                    endif
                    d = d/10.**i ! now d is between 1 and 10
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
                        print '(a,f16.4,a,f16.4)'
     +                       ,'# binsize rounded from ',s,' to ',d
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
                        print '(a,f16.4,a,f16.4)'
     +                       ,'# minimum rounded from ',s,' to ',mindat
                    endif
                end if
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
                        if ( i.ge.yr1 .and. i.le.yr2 ) then
                            do iens=nens1,nens2
                                if ( data(j,i,iens).lt.1e33 )
     +                               then
                                    n = 1 + int((data(j,i,iens)
     +                                   -mindat)/d)
                                    if ( n.lt.0 .or. n.gt.nbin )
     +                                   then
                                        write(*,*)
     +                                       '# histogram warning
     +                                       :'//' n>nbin:',n,nbin
                                    else
                                        nn(n) = nn(n) + 1
                                        if ( data(j,i,iens).lt. mean(0)
     +                                       ) nless= nless + 1
                                    end if
                                endif
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
                                yy(i) = -999.9
                            else
                                yy(i) = ntot*(exp(-(1+xi*x2)**(-1/xi)) -
     +                                exp(-(1+xi*x1)**(-1/xi)))
                            endif
                        endif
                    enddo
                elseif ( nfit.eq.6 ) then
*                   GPD distribution
                    if ( lchangesign ) then
                        nnn = 0
                        do i=1,nbin
                            x2 = -(mindat + d*i + mindata)/b
                            print *,i,x2
                            if ( x2.gt.0 ) nnn = nnn + nn(i)
                        enddo
                    else
                        nnn = 0
                        do i=1,nbin
                            x1 = (mindat + d*(i-1) - mindata)/b
                            if ( x1.gt.0 ) nnn = nnn + nn(i)
                        enddo
                    endif
                    do i=1,nbin
                        if ( lchangesign ) then
                            x1 = -(mindat + d*(i-1) + mindata)/b
                            x2 = -(mindat + d*i     + mindata)/b
                        else
                            x1 = (mindat + d*(i-1) - mindata)/b
                            x2 = (mindat + d*i     - mindata)/b
                        endif
                        if ( (.not.lchangesign .and. x1.lt.0) .or. 
     +                       (lchangesign .and. x2.lt.0) ) then
                            yy(i) = -999.9
                        elseif ( abs(xi).le.1e-4 ) then
                            yy(i) = nnn*(exp(-x1+0.5*x1**2*xi)
     +                           - exp(-x2+0.5*x2**2*xi))/snorm
                        else
                            if ( 1+xi*x1.le.0 .or. 1+xi*x2.le.0 )
     +                           then
                                yy(i) = -999.9
                            else
                                yy(i) = nnn*((1+xi*x1)**(-1/xi)
     +                               - (1+xi*x2)**(-1/xi))/snorm
                            endif
                        endif
                        if ( lwrite ) write(*,*) x1,x2,yy(i),snorm,nnn
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
                    nnn = 0
                    do i=1,nbin
                        if ( nfit.eq.6 ) then
!                           only count the points above the threshold...
                            nnn = nnn + nn(i)
                            if ( lchangesign ) then
                                x2 = mindat + d*(i-1)
                                if ( lwrite ) then
                                    print *,'i,x2,nn(i),nnn,yy(i) = ',i
     +                                   ,x2,nn(i),nnn,yy(i)
                                    print *,'x2,mindata = ',x2,mindata
                                endif
                                if ( x2.lt.-mindata ) cycle
                            else
                                x1 = mindat + d*i
                                if ( lwrite ) then
                                    print *,'i,x1,nn(i),nnn,yy(i) = ',i
     +                                   ,x1,nn(i),nnn,yy(i)
                                    print *,'x1,mindata = ',x1,mindata
                                endif
                                if ( x1.lt.mindata ) cycle
                            endif
                        endif
                        if ( yy(i).ne.-999.9 .and. (
     +                       yy(i).gt.0.1 .or. nn(i).gt.0 .or. n.gt.0 )
     +                       ) then
                            if ( yy(i).ne.0 .or. nn(i).ne.0 ) then
                                if ( nn(i).gt.0 .and. yy(i).le.0 ) then
                                    write(0,*)
     +                                   'histogram: error: events '//
     +                                   'with p=0 occurred: ',i,nn(i)
                                endif
                                n = n + 1
                                xn(n) = nn(i)
                                xy(n) = yy(i)
                            endif
                        elseif ( n.eq.0 ) then
                            if ( lwrite ) print *,'n = 0!'
                        endif
                    enddo
                    do i=n,1,-1
                        if ( lwrite ) print *
     +                       ,'looking for non-zero bin: ',i,xy(i),xn(i)
                        if ( xy(i).gt.0.1 .or. xn(i).gt.0 ) goto 410
                    enddo
 410                continue
                    n = i
                    if ( lwrite ) then
                        do i=1,n
                            print '(a,i4,2f9.2,g10.2)','# ',
     +                           i,xn(i),xy(i),(xn(i)-xy(i))**2/xy(i)
                        enddo
                    endif
*                   Numerical Recipes, p615
                    if ( n.gt.nfitted ) then
                        call chsone(xn,xy,n,nfitted,df,chsq,prob)
                        if ( plot ) write(11,'(g16.4)') chsq/df
                        if ( lweb ) then
                            print '(a,g9.2,a,i4,a,g9.3,a,f6.4)',
     +                           '# <tr><td>&chi;^2/df</td><td>'
     +                           ,chsq,'/',nint(df),' =',chsq/df
     +                           ,'</td><td>p=',prob
                        else
                            print '(a,g9.2,a,i4,a,g9.3)','# chi^2/df:',
     +                           chsq,'/',nint(df),' =',chsq/df
                            print '(a,f8.4,a)','# probability:   ',
     +                           100*prob,
     +                           '% (this is the probability that the '
     +                           //'observed distribution was drawn '
     +                           //'from the fitted one)'
                        endif
                        if ( plot ) write(11,'(g16.4)') prob
                    else
                        if ( lwrite ) print *,'# not enough non-zero '
     +                       //'bins to determine \chi^2: ',n,nfitted
                    endif
                endif
*
*               print out
*
                n = 0
                s = 0
                do i=1,nbin
                    n = n + nn(i)
                    if ( yy(i).ne.-999.9 ) then
                        s = s + yy(i)
                    else
                        s = n
                    endif
                    x = mindat+ d*(i-0.5)
                    if ( lchangesign ) then
                        x = -x
                    endif
                    print '(i4,g12.4,2i9,2f12.2)',i,x,nn(i),n,yy(i),s
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
                        write(*,*) '# histogram error: can only make ',
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
                        ff = 2*snorm*f
                        if ( minindx.gt.-1e33 ) then
                            ff = ff + erf((minindx-a)/(sqrt2*b))
                        else
                            ff = ff - 1
                        endif
*                       IMSL routine
                        call merfi(ff,z,ier)
                        s = a + sqrt2*b*z
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
                    elseif ( nfit.eq.6 ) then
*                       GPD distribution
                        write(0,*) 'QQ plot for GPD not yet ready'
                        call abort
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
            elseif ( ntype.eq.2 .or. ntype.eq.3 .or. ntype.eq.4 ) then
                if ( plot ) then
                    do i=1,10
                        write(11,'(3g16.4)') t(i),t25(i),t975(i)
                    enddo
                    write(11,'(3g16.4)') tx,tx25,tx975
                endif
                if ( lchangesign ) xyear = -xyear
                frac = 1
                call plot_ordered_points(xx,xs,yrs,ntot,ntype,nfit,
     +               frac,a,b,xi,j1,j2,minindx,mindata,pmindata,
     +               yr2a,xyear,snorm,lchangesign,lwrite,.true.)
            else
                write(0,*) 'Unknown plot type ',ntype
                write(*,*) 'Unknown plot type ',ntype
                call abort
            endif
        enddo
 999    continue
        call savestartstop(yrstart,yrstop)
        end
*  #] histogram:

