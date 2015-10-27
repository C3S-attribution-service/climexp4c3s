subroutine attribute_dist(series,nperyear,covariate,nperyear1,npermax,yrbeg,yrend,&
&   mens1,mens,assume,distribution,seriesids,results,nresmax,nresults,lprint)
!
!   take block maxima, convert to linear arrays and call fitgevcov / fitgumcov or
!   take average, convert to linear arrays and call fitgaucov / fitgpdcov
!   if ( lprint ) print the output to screen, if not return in results(1:nresults)
!
    implicit none
    include "getopts.inc"
    integer nperyear,nperyear1,npermax,yrbeg,yrend,mens1,mens,nresmax,nresults
    real series(npermax,yrbeg:yrend,0:mens),covariate(npermax,yrbeg:yrend,0:mens)
    real results(3,nresmax)
    character seriesids(0:mens)*(*),assume*(*),distribution*(*)
    integer fyr,lyr,ntot,i,j,k,ntype,nmax,npernew,j1,j2,iens,jens,ensmax,init,ndecor,n
    integer,allocatable :: yrs(:)
    real a(3),b(3),xi(3),alpha(3),beta(3),cov1,cov2,offset,t(3,10,3),tx(3,3)
    real,allocatable :: xx(:,:),yrseries(:,:,:),yrcovariate(:,:,:),yy(:),crosscorr(:,:)
    logical lboot,lprint,subtract_offset
    character operation*4,file*1024,idmax*30
    data init /0/

    results = 3e33
    nresults = 0
    fyr = yr1
    lyr = yr2
    nmax = nperyear*(yr2-yr1+1)*(1+mens-mens1)
    !!!write(0,*) 'assume,distribution = ',assume,distribution
    allocate(xx(2,nmax))
    allocate(yrs(0:nmax))
    allocate(yy(nmax))
    allocate(crosscorr(0:mens,0:mens))
    xx = 3e33
    if ( lwrite ) print *,'attribute_dist: series(:,2000,0) = ',(series(j,2000,0),j=1,min(15,nperyear))
    
    if ( distribution.eq.'gev' .or. distribution.eq.'gumbel' ) then
        allocate(yrseries(1,fyr:lyr,0:mens))
        allocate(yrcovariate(1,fyr:lyr,0:mens))
        yrseries = 3e33
        yrcovariate = 3e33
        if ( nperyear.eq.nperyear1 .and. nperyear.gt.1 ) then
            if ( lwrite ) print *,'attribute_dist: calling make_two_annual_values with max'
            call make_two_annual_values(series,covariate,nperyear,npermax,yrbeg,yrend,mens1,mens, &
            &   yrseries,yrcovariate,fyr,lyr,'max')
            if ( xyear.lt.1e33 ) then
                call get_covariate_extrayear(covariate,nperyear1,npermax,yrbeg,yrend,mens1,mens, &
                & yrcovariate,fyr,lyr)
            end if
        else
            if ( lwrite ) print *,'attribute_dist: calling make_annual_values for series with max'
            call make_annual_values(series,nperyear,npermax,yrbeg,yrend,mens1,mens, &
            &   yrseries,fyr,lyr,'max')
            if ( lwrite ) print *,'attribute_dist: calling make_annual_values for covariate with mean'
            call make_annual_values(covariate,nperyear1,npermax,yrbeg,yrend,mens1,mens, &
            &   yrcovariate,fyr,lyr,'mean')
        end if
        npernew = 1
        j1 = 1
        j2 = 1
        m1 = 1
        lsel = 1
        ndecor = 1+nint(decor) ! not used if all goes well
        if ( init == 0 ) then
            init = 1
            call print_bootstrap_message(max(1,ndecor),j1,j2)
        end if
    else if ( distribution.eq.'gpd' .or. distribution.eq.'gauss' ) then
        ! in the others lsum inicates teh block maxima length...
        call getj1j2(j1,j2,m1,nperyear,lwrite)
        decor = max(decor,real(lsum)-1)
        if ( j1.eq.j2 ) then
            ndecor = 1+int(decor/nperyear)
        else
            ndecor = 1+int(decor)
        endif
        if ( init == 0 ) then
            init = 1
            call print_bootstrap_message(max(1,ndecor),j1,j2)
        end if
        allocate(yrseries(nperyear,fyr:lyr,0:mens))
        allocate(yrcovariate(nperyear,fyr:lyr,0:mens))
        yrseries = 3e33
        yrcovariate = 3e33
        ! copy series to keep the code easier to handle GPD and Gauss at the same time :-(
        yrseries(1:nperyear,fyr:lyr,mens1:mens) = series(1:nperyear,fyr:lyr,mens1:mens)
        ! change covariate to the same time resolution as series
        if ( nperyear1.lt.nperyear ) then
            if ( lwrite ) print *,'repeating covariate from nperyear = ',nperyear1,' to ',nperyear
           do iens=mens1,mens
                call annual2shorter(covariate(1,yrbeg,iens),npermax,yrbeg,yrend,nperyear1, &
                &   yrcovariate(1,fyr,iens),nperyear,fyr,lyr,nperyear,1,nperyear,1,lwrite)
            end do
        else if ( nperyear1.gt.nperyear ) then
            ! this should not occur in the web interface
            write(0,*) 'atribute_dist: error: covariate should not have higher time resolution than series: ', &
            &   nperyear1,nperyear
            write(*,*) 'atribute_dist: error: covariate should not have higher time resolution than series: ', &
            &   nperyear1,nperyear
            call abort
        else ! equal already
            if ( lwrite ) print *,'copying covariate ',nperyear
            yrcovariate(1:nperyear,fyr:lyr,mens1:mens) = covariate(1:nperyear,fyr:lyr, &
            &   mens1:mens)
        end if
        npernew = nperyear
        call getj1j2(j1,j2,m1,npernew,.false.)
    else
        write(*,*) 'attribute_dist: error: unknown distribution ',trim(distribution)
        write(0,*) 'attribute_dist: error: unknown distribution ',trim(distribution)
        call abort
    end if
    
    if ( assume == 'scale' ) then
        call checknonegative(yrseries,npernew,fyr,lyr,mens1,mens,j1,j2,assume, &
        &   lchangesign,lwrite)
    end if
    if ( lnormsd .and. mens.gt.mens1 ) then
        if ( lwrite ) print *,'attrbute_dist: normalising all series'
        call normaliseseries(yrseries,npernew,fyr,lyr,mens1,mens,j1,j2,assume,lwrite)
    end if
    
    if ( lwrite ) print *,'attribute_dist: calling handle_then_now'
    call handle_then_now(yrseries,yrcovariate,npernew,fyr,lyr,j1,j2,yr1a,yr2a,mens1,mens, &
    & xyear,ensmax,cov1,cov2,lprint,lwrite)
    if ( cov1.gt.1e33 .or. cov2.gt.1e33 ) then
        if ( lwrite ) print *,'giving up, cov1,cov2 = ',cov1,cov2
        return
    end if
    if ( ensmax.ge.mens1 ) then
        idmax = seriesids(ensmax)
    else
        idmax = ' ' ! undefined
    end if
    if ( lprint ) then
        if ( namestring.ne.' ' ) then
            print '(4a)','# <tr><th colspan=4>',trim(namestring),'</th></tr>'
        endif
        if ( abs(nint(confidenceinterval)-confidenceinterval).lt.0.0001 ) then
            print '(a,i2,a)' ,'# <tr><th>parameter</th><th>year</th><th>value</th><th>' &
        &       ,nint(confidenceinterval),'% CI</th></tr>'
        else
            print '(a,f6.3,a)' ,'# <tr><th>parameter</th><th>year</th><th>value</th><th>' &
        &       ,confidenceinterval,'% CI</th></tr>'
        end if
        print '(a,i4,a,g16.5,a)','# <tr><td>covariate:</td><td>',yr1a,'</td><td>',cov1, &
        &   '</td><td>&nbsp;</td></tr>'
        print '(a,i4,a,g19.5,a)','# <tr><td>&nbsp;</td><td>',yr2a,'</td><td>',cov2, &
        &   '</td><td>&nbsp;</td></tr>'
    end if
    subtract_offset = .true.
    if ( subtract_offset ) then
        if ( lwrite ) print *,'attribute_dist: calling subtract_constant'
        call subtract_constant(yrcovariate,yrseries,npernew,fyr,lyr,mens1,mens, &
        &   cov1,cov2,offset,lwrite)
    else
        offset = 0
    end if
    if ( mens.gt.mens1 ) then
        n = (mens-mens1+1)
        n = n*(n-1)/2
        i = 0
        do iens = mens1,mens
            crosscorr(iens,iens) = 1
            do jens=iens+1,mens
                call keepalive1('Computing correlations',i,n)
                i = i + 1
                call getcorr(yrseries(1,fyr,iens),npernew,fyr,lyr, &
                &            yrseries(1,fyr,jens),npernew,fyr,lyr, &
                &            j1,j2,crosscorr(iens,jens),lwrite)
                crosscorr(jens,iens) = crosscorr(iens,jens)
            end do
        end do
    else
        crosscorr(mens1,mens1) = 1
    end if

    lboot = .true.
    if ( distribution.eq.'gev' ) then
        ntype = 2 ! Gumbel plot
        if ( lwrite ) print *,'attribute_dist: calling fitgevcov',j1,j2
        call fitgevcov(yrseries,yrcovariate,npernew,fyr,lyr,mens1,mens & 
    &       ,crosscorr,a,b,xi,alpha,beta,j1,j2 &
    &       ,lweb,ntype,lchangesign,yr1a,yr2a,xyear,idmax,cov1,cov2,offset &
    &       ,t,tx,restrain,assume,confidenceinterval,ndecor,lboot,lprint,dump,plot,lwrite)
    else if ( distribution.eq.'gpd' ) then
        ntype = 3 ! log plot
        !!!print *,'DEBUG'
        !!!lboot = .false.
        !!!lwrite = .true.
        if ( lwrite ) print *,'attribute_dist: calling fitgpdcov'
        call fitgpdcov(yrseries,yrcovariate,npernew,fyr,lyr,mens1,mens & 
    &       ,crosscorr,a,b,xi,alpha,beta,j1,j2 &
    &       ,lweb,ntype,lchangesign,yr1a,yr2a,xyear,idmax,cov1,cov2,offset &
    &       ,t,tx,pmindata,restrain,assume,confidenceinterval,ndecor,lboot &
    &       ,lprint,dump,plot,lwrite)
    else if  ( distribution.eq.'gumbel' ) then
        ntype = 2 ! Gumbel plot
        xi = 0
        if ( lwrite ) print *,'attribute_dist: calling fitgumcov'
        call fitgumcov(yrseries,yrcovariate,npernew,fyr,lyr,mens1,mens & 
    &       ,crosscorr,a,b,alpha,beta,j1,j2 &
    &       ,lweb,ntype,lchangesign,yr1a,yr2a,xyear,idmax,cov1,cov2,offset &
    &       ,t,tx,assume,confidenceinterval,ndecor,lboot,lprint,dump,plot,lwrite)
    else if  ( distribution.eq.'gauss' ) then
        ntype = 4 ! sqrtlog plot
        xi = 0
        if ( lwrite ) print *,'attribute_dist: calling fitgaucov'
        call fitgaucov(yrseries,yrcovariate,npernew,fyr,lyr,mens1,mens & 
    &       ,crosscorr,a,b,alpha,beta,j1,j2 &
    &       ,lweb,ntype,lchangesign,yr1a,yr2a,xyear,idmax,cov1,cov2,offset &
    &       ,t,tx,assume,confidenceinterval,ndecor,lboot,lprint,dump,plot,lwrite)
    else
        write(0,*) 'attribute_dist: error: unknown distribution ',trim(distribution)
    end if
    if ( .not.lprint ) then
        results(:,1) = a
        results(:,2) = b
        results(:,3) = xi
        results(:,4) = alpha
        results(:,5) = beta
        do i=1,3
            results(:,5+i) = tx(:,i)
        end do
        k = 8
        do j=1,3
            do i=1,10
                k = k + 1
                if ( k.gt.nresmax ) then
                    write(0,*) 'attrubute_dist: internal error: results array to small: ',nresmax
                    call abort
                end if
                results(:,k) = t(:,i,j)
            end do
        end do
        nresults = k
    end if

end subroutine attribute_dist

subroutine attribute_init(file,distribution,assume,off,nperyear,yrbeg,yrend,nensmax,lwrite)
!
!   initialisation common to time series and field programs
!
    implicit none
    integer off,nperyear,yrbeg,yrend,nensmax
    character file*(*),distribution*(*),assume*(*)
    logical lwrite
    character string*255,string1*255,string2*200   
    integer iargc
!
    call killfile(string,string1,string2,0) ! random strings
!   the usual call to getopts comes way too late for these options...
    nperyear = 12
    call getarg(1,file)
    if ( file /= 'file' ) then
        off = 0
    else
        off = 2
    end if
    call getopts(6+off,iargc(),nperyear,yrbeg,yrend,.false.,0,nensmax)

    call getarg(3+off,distribution)
    call tolower(distribution)
    if ( distribution.ne.'gev' .and. distribution.ne.'gumbel' .and. &
    &    distribution.ne.'gpd' .and. distribution.ne.'gauss' ) then
        write(0,*) 'attribute: error: only GEV, GPD or Gauss supported, not ',distribution
        call abort
    end if

    call getarg(4+off,string)
    if ( string(1:4).ne.'assu' ) then
        write(0,*) 'attribute: error: expecting "assume", not ',trim(string)
        call abort
    end if
    call getarg(5+off,assume)
    call tolower(assume)
    if ( assume.ne.'shift' .and. assume.ne.'scale' .and. assume.ne.'both' ) then
        write(0,*) 'attribute: error: only shift, scale or both supported, not ',assume
        call abort
    end if
end subroutine

subroutine getdpm(dpm,nperyear)
    ! make a list of number of days per month
    implicit none
    integer dpm(12),nperyear
    integer dpm366(12)
    data dpm366 /31,29,31,30,31,30,31,31,30,31,30,31/
    if ( nperyear.eq.366 ) then
        dpm = dpm366
    else if ( nperyear.eq.365 ) then
        dpm = dpm366
        dpm(2) = 28
    else if ( nperyear.eq.360 ) then
        dpm = 30
    else if ( nperyear.lt.360 .and. nperyear.ge.12 ) then
        dpm = nperyear/12
    else if ( nperyear.le.4 .and. nperyear.ge.1 ) then
        dpm = 1
    else
        write(0,*) 'getdpm: error: unknown nperyear = ',nperyear
        call abort
    end if
end subroutine getdpm

subroutine make_annual_values(series,nperyear,npermax,yrbeg,yrend,mens1,mens, &
&   yrseries,fyr,lyr,operation)
    
    ! construct an annual time series with the maxima

    implicit none
    include 'getopts.inc'
    integer nperyear,npermax,yrbeg,yrend,mens1,mens,fyr,lyr
    real series(npermax,yrbeg:yrend,0:mens),yrseries(1,fyr:lyr,0:mens)
    character operation*(*)
    integer j1,j2,yy,yr,mm,mo,dd,dy,k,m,mtot,n,dpm(12),iens
    real s

    if ( lwrite ) then
        print *,'make_annual_values: taking ',trim(operation)
        print *,'nperyear,npermax = ',nperyear,npermax
    end if
    if ( nperyear == 1 ) then
        do iens=mens1,mens
            do yy=yr1,yr2
                yrseries(1,yy,iens) = series(1,yy,iens)
            end do
        end do
    else if ( nperyear < 12 ) then
        do iens=mens1,mens
            do yr=yr1,yr2
                m = 0
                if ( operation.eq.'max' ) then
                    s = -3e33
                else if ( operation.eq.'min' ) then
                    s = 3e33
                else
                    s = 0
                end if
                do dy=1,nperyear
                    if ( series(dy,yr,iens).lt.1e33 ) then
                        m = m + 1
                        if ( operation.eq.'max' ) then
                            s = max(s,series(dy,yr,iens))
                        else if ( operation.eq.'min' ) then
                            s = min(s,series(dy,yr,iens))
                        else if ( operation.eq.'mean' ) then
                            s = s + series(dy,yr,iens)
                        else
                            write(0,*) 'make_annual_values: unknown operation ',trim(operation)
                            call abort
                        end if
                    end if
                end do
                if ( m.gt.minfac*nperyear ) then
                    if ( operation.eq.'min' .or. operation.eq.'max' ) then
                        yrseries(1,yr,iens) = s
                    else if ( operation.eq.'mean' ) then
                        yrseries(1,yr,iens) = s/m
                    else
                        write(0,*) 'make_annual_values: unknown operation ',trim(operation)
                        call abort
                    end if
                end if
            end do ! yr
        end do ! iens                
    else if ( nperyear >= 12 ) then
        call getj1j2(j1,j2,m1,12,lwrite)
        if ( lwrite ) print *,'make_annual_values: j1,j2,nperyear = ',j1,j2,nperyear
        call getdpm(dpm,nperyear)
        if ( lwrite ) print *,'                    dpm = ',dpm
        do iens=mens1,mens
            do yy=yr1,yr2
                if ( operation.eq.'max' ) then
                    s = -3e33
                else if ( operation.eq.'min' ) then
                    s = 3e33
                else
                    s = 0
                end if
                m = 0
                mtot = 0
                dd = 0
                do mm=1,j1-1
                    dd = dd + dpm(mm)
                end do
                do mm=j1,j2
                    mo = mm
                    call normon(mo,yy,yr,min(12,nperyear))
                    do k=dd+1,dd+dpm(mo)
                        dy = k
                        call normon(dy,yy,yr,nperyear)
                        if ( yr.le.yr2 ) then
                            mtot = mtot + 1
                            if ( series(dy,yr,iens).lt.1e33 ) then
                                m = m + 1
                                if ( operation.eq.'max' ) then
                                    s = max(s,series(dy,yr,iens))
                                else if ( operation.eq.'min' ) then
                                    s = min(s,series(dy,yr,iens))
                                else if ( operation.eq.'mean' ) then
                                    s = s + series(dy,yr,iens)
                                else
                                    write(0,*) 'make_annual_values: unknown operation ',trim(operation)
                                    call abort
                                end if
                            end if
                        end if
                    end do
                    dd = dd + dpm(mo)
                end do
                if ( m.gt.minfac*mtot ) then
                    if ( operation.eq.'min' .or. operation.eq.'max' ) then
                        yrseries(1,yy,iens) = s
                    else if ( operation.eq.'mean' ) then
                        yrseries(1,yy,iens) = s/m
                    else
                        write(0,*) 'make_annual_values: unknown operation ',trim(operation)
                        call abort
                    end if
                end if
            end do ! yy
        end do ! iens
    end if ! nperyear >= 12
end subroutine make_annual_values

subroutine make_two_annual_values(series,covariate,nperyear,npermax,yrbeg,yrend,mens1,mens, &
&   yrseries,yrcovariate,fyr,lyr,operation)
    
    ! construct two annual time series with the maxima of series and the corresponding values of covariate

    implicit none
    include 'getopts.inc'
    integer nperyear,npermax,yrbeg,yrend,mens1,mens,fyr,lyr
    real series(npermax,yrbeg:yrend,0:mens),covariate(npermax,yrbeg:yrend,0:mens), &
    & yrseries(1,fyr:lyr,0:mens),yrcovariate(1,fyr:lyr,0:mens)
    character operation*(*)
    integer j1,j2,yy,yr,mm,mo,dd,dy,k,m,mtot,n,dpm(12),iens
    real s

    if ( lwrite ) then
        print *,'make_two_annual_values: taking ',trim(operation)
        print *,'nperyear,npermax = ',nperyear,npermax
    end if
    if ( nperyear == 1 ) then
        do iens=mens1,mens
            do yy=yr1,yr2
                yrseries(1,yy,iens) = series(1,yy,iens)
                yrcovariate(1,yy,iens) = covariate(1,yy,iens)
            end do
        end do
    else if ( nperyear < 12 ) then
        do iens=mens1,mens
            do yr=yr1,yr2
                m = 0
                if ( operation.eq.'max' ) then
                    s = -3e33
                else if ( operation.eq.'min' ) then
                    s = 3e33
                else
                    s = 0
                    yrcovariate(1,yr,iens) = 0
                end if
                do dy=1,nperyear
                    if ( series(dy,yr,iens).lt.1e33 .and. covariate(dy,yr,iens).lt.1e33 ) then
                        m = m + 1
                        if ( operation.eq.'max' ) then
                            if ( series(dy,yr,iens).gt.s ) then
                                s = series(dy,yr,iens)
                                yrcovariate(1,yr,iens) = covariate(dy,yr,iens)
                            end if
                        else if ( operation.eq.'min' ) then
                            if ( series(dy,yr,iens).lt.s ) then
                                s = series(dy,yr,iens)
                                yrcovariate(1,yr,iens) = covariate(dy,yr,iens)
                            end if
                        else if ( operation.eq.'mean' ) then
                            s = s + series(dy,yr,iens)
                            yrcovariate(1,yr,iens) = yrcovariate(1,yr,iens) + covariate(dy,yr,iens)
                        else
                            write(0,*) 'make_annual_values: unknown operation ',trim(operation)
                            call abort
                        end if
                    end if
                end do
                if ( m.gt.minfac*nperyear ) then
                    if ( operation.eq.'min' .or. operation.eq.'max' ) then
                        yrseries(1,yr,iens) = s
                    else if ( operation.eq.'mean' ) then
                        yrseries(1,yr,iens) = s/m
                        yrcovariate(1,yr,iens) = yrcovariate(1,yr,iens)/m
                    else
                        write(0,*) 'make_annual_values: unknown operation ',trim(operation)
                        call abort
                    end if
                    if ( lwrite ) print *,yr,yrcovariate(1,yr,iens),yrseries(1,yr,iens)
                end if
            end do ! yr
        end do ! iens                
    else if ( nperyear >= 12 ) then
        call getj1j2(j1,j2,m1,12,lwrite)
        if ( lwrite ) print *,'make_annual_values: j1,j2,nperyear = ',j1,j2,nperyear
        call getdpm(dpm,nperyear)
        if ( lwrite ) print *,'                    dpm = ',dpm
        do iens=mens1,mens
            do yy=yr1,yr2
                if ( operation.eq.'max' ) then
                    s = -3e33
                else if ( operation.eq.'min' ) then
                    s = 3e33
                else
                    s = 0
                    yrcovariate(1,yr,iens) = 0
                end if
                m = 0
                mtot = 0
                dd = 0
                do mm=1,j1-1
                    dd = dd + dpm(mm)
                end do
                do mm=j1,j2
                    mo = mm
                    call normon(mo,yy,yr,min(12,nperyear))
                    do k=dd+1,dd+dpm(mo)
                        dy = k
                        call normon(dy,yy,yr,nperyear)
                        if ( yr.le.yr2 ) then
                            mtot = mtot + 1
                            if ( series(dy,yr,iens).lt.1e33 ) then
                                m = m + 1
                                if ( operation.eq.'max' ) then
                                    s = max(s,series(dy,yr,iens))
                                    yrcovariate(1,yr,iens) = covariate(dy,yr,iens)
                                else if ( operation.eq.'min' ) then
                                    s = min(s,series(dy,yr,iens))
                                    yrcovariate(1,yr,iens) = covariate(dy,yr,iens)
                                else if ( operation.eq.'mean' ) then
                                    s = s + series(dy,yr,iens)
                                    yrcovariate(1,yr,iens) = yrcovariate(1,yr,iens) + covariate(dy,yr,iens)
                                else
                                    write(0,*) 'make_annual_values: unknown operation ',trim(operation)
                                    call abort
                                end if
                            end if
                        end if
                    end do
                    dd = dd + dpm(mo)
                end do
                if ( m.gt.minfac*mtot ) then
                    if ( operation.eq.'min' .or. operation.eq.'max' ) then
                        yrseries(1,yy,iens) = s
                    else if ( operation.eq.'mean' ) then
                        yrseries(1,yy,iens) = s/m
                        yrcovariate(1,yr,iens) = yrcovariate(1,yr,iens)/m
                    else
                        write(0,*) 'make_annual_values: unknown operation ',trim(operation)
                        call abort
                    end if
                end if
            end do ! yy
        end do ! iens
    end if ! nperyear >= 12
end subroutine make_two_annual_values

subroutine get_covariate_extrayear(covariate,nperyear,npermax,yrbeg,yrend,mens1,mens,yrcovariate,fyr,lyr)
    
    ! construct two annual time series with the maxima of series and the corresponding values of covariate

    implicit none
    include 'getopts.inc'
    integer nperyear,npermax,yrbeg,yrend,mens1,mens,fyr,lyr
    real covariate(npermax,yrbeg:yrend,0:mens),yrcovariate(1,fyr:lyr,0:mens)
    integer n,j1,j2,iens,i,j,k
    real s
    
    if ( yrcovariate(1,yr2a,0).gt.1e33 ) then
        ! it may have been left out because there is no data in series
        ! take the mean value over the months (we do not know which time of year the
        ! user-supplied value refers to...)
        s = 0
        n = 0
        call getj1j2(j1,j2,m1,nperyear,.false.)
        do iens = mens1,mens
            do j=j1,j2
                k = j
                call normon(k,yr2a,i,nperyear)
                if ( covariate(k,i,iens).lt.1e33 ) then
                    n = n + 1
                    s = s + covariate(k,i,iens)
                end if
            end do
        end do
        if ( n.gt.0 ) then
            yrcovariate(1,yr2a,mens1:mens) = s/n
        end if
    end if
end subroutine

subroutine fill_linear_array(series,covariate,nperyear,j1,j2,fyr,lyr,mens1,mens,&
&   xx,yrs,nmax,ntot,lwrite)

    ! transfer the valid pairs in series, covariate to xx(1:2,1:ntot)
    
    implicit none
    integer nperyear,j1,j2,fyr,lyr,mens1,mens,nmax,ntot
    integer yrs(0:nmax)
    real series(nperyear,fyr:lyr,0:mens),covariate(nperyear,fyr:lyr,0:mens),xx(2,nmax)
    integer yy,yr,mm,mo,day,month,yrstart,yrstop,iens
    logical lwrite

    if ( lwrite ) then
        print *,'fill_linear_array: nperyear,j1,j2,fyr,lyr,mens1,mens = ', &
        &   nperyear,j1,j2,fyr,lyr,mens1,mens
        if ( .true. ) then
            do yr=fyr,lyr
                if ( series(j1,yr,mens1).lt.1e33 .and. &
     &               covariate(j1,yr,mens1).lt.1e33 ) then
                    print *,yr,series(j1,yr,mens1),covariate(j1,yr,mens1)
                end if
            end do
        end if
    end if
    yrstart = lyr
    yrstop = fyr
    ntot = 0
    do iens=mens1,mens
        do yy=fyr,lyr
            do mm=j1,j2
                mo = mm
                call normon(mo,yy,yr,nperyear)
                if ( yr.ge.fyr .and. yr.le.lyr ) then
                    if ( series(mo,yr,iens).lt.1e33 .and. covariate(mo,yr,iens).lt.1e33 ) then
                        ntot = ntot + 1
                        xx(1,ntot) = series(mo,yr,iens)
                        xx(2,ntot) = covariate(mo,yr,iens)
                        call getdymo(day,month,mo,nperyear)
                        yrs(ntot) = 10000*yr + 100*month + day
                        if ( nperyear.gt.1000 ) then
                            yrs(ntot) = 100*yrs(ntot) + mod(ntot,nint(nperyear/366.))
                            end if
                        yrstart = min(yrstart,yr)
                        yrstop = max(yrstop,yr)
                    end if
                end if
            end do
        end do
    end do
    call savestartstop(yrstart,yrstop)
end subroutine fill_linear_array

subroutine sample_bootstrap(series,covariate,nperyear,j1,j2,fyr,lyr,mens1,mens,&
&   crosscorr,ndecor,xx,nmax,ntot,sdecor,lwrite)

    ! transfer random samples of the valid pairs in series, covariate to xx(1:2,1:ntot)
    ! serial autocrrelations are taken into account with a moving block of length ndecor
    ! for j1==j2 this is counted in years, for j1/=j2 in months/days
    
    implicit none
    integer nperyear,j1,j2,fyr,lyr,mens1,mens,ndecor,nmax,ntot
    real series(nperyear,fyr:lyr,0:mens),covariate(nperyear,fyr:lyr,0:mens),sdecor
    real crosscorr(0:mens,0:mens),xx(2,nmax)
    logical lwrite
    integer yy,yr,mm,mmm,mo,day,month,yrstart,yrstop,iens,i,j,ntry,nleft,jens,nens
    real ranf,cutoff
    logical ok

    cutoff = exp(-0.5)  ! based on the moving block going to 1/e, 
                        ! this is half the distance as we go in all directions.
    cutoff = exp(-1.)
    if ( lwrite ) then
        print *,'sample_bootstrap: nperyear,j1,j2,fyr,lyr,mens1,mens = ', &
        &   nperyear,j1,j2,fyr,lyr,mens1,mens
        print *,'                  ndecor = ',ndecor
    end if
!
!   as this routine will often be called with far too large ranges in years,
!   let's first narrow that down to get an acceptable hit rate later on
!
    yrstart = lyr
    yrstop = fyr
    ntot = 0
    do iens=mens1,mens
        do yy=fyr,lyr
            do mm=j1,j2
                mo = mm
                call normon(mo,yy,yr,nperyear)
                if ( yr.ge.fyr .and. yr.le.lyr ) then
                    if ( series(mo,yr,iens).lt.1e33 .and. covariate(mo,yr,iens).lt.1e33 ) then
                        ntot = ntot + 1
                        yrstart = min(yrstart,yr)
                        yrstop = max(yrstop,yr)
                    end if
                end if
            end do
        end do
    end do
    if ( lwrite ) print *,'sample_bootstrap: yrstart,yrstop = ',yrstart,yrstop
 !
 !  the sampling itself.
 !
    j = 0
    nens = 0
    do while ( j < ntot )
        nens = nens + 1
        ok = .false.
        ntry = 0
        do while ( .not.ok )
            ntry = ntry + 1
            if ( ntry.gt.10000 ) then
                write(0,*) 'sample_bootstrap: too many tries ',ntry
                call exit(-1)
            end if
            if ( mens1 == mens ) then
                iens = mens1
            else
                call random_number(ranf)
                iens = 0 + int((mens-mens1+1)*ranf)
            end if
            call random_number(ranf)
            yy = yrstart + int((yrstop-yrstart+1)*ranf)
            if ( j1 == j2 ) then
                mo = j1
                yr = yy
            else
                call random_number(ranf)
                mo = j1 + int((j2-j1+1)*ranf)
                call normon(mo,yy,yr,nperyear)
            end if
            if ( lwrite ) print *,'reference ',mo,yr,iens
            do jens = mens1,mens
                if ( crosscorr(jens,iens).lt.1e33 .and. &
                &    crosscorr(jens,iens).gt.cutoff ) then ! always true for the diagonal
                    ! include in spatial moving block
                    if ( series(mo,yr,jens).lt.1e33 .and. covariate(mo,yr,jens).lt.1e33 ) then
                        ok = .true.
                        j = j + 1
                        if ( lwrite ) print *,'crosscorr(',jens,iens, &
                        &   ') = ',crosscorr(jens,iens),'>',cutoff
                        xx(1,j) = series(mo,yr,jens)
                        xx(2,j) = covariate(mo,yr,jens)
                        if ( lwrite ) print *,'xx(:,',j,') = ',xx(:,j),mo,yr,jens
                        if ( ndecor.gt.1 ) then
                            ! handle serial autocorrelations with a moving block
                            if ( j1 == j2 ) then
                                ! annual values, ndecor refers to year-on-year autocorrelations
                                ! just skip undefineds, we do the same with the end of the series.
                                do yy=yr+1,yr+ndecor-1
                                    if ( yy <= yrstop .and. j < ntot ) then
                                        if ( series(mo,yy,jens).lt.1e33 .and. &
                                        &    covariate(mo,yy,jens).lt.1e33 ) then
                                            j = j + 1
                                            if ( j.gt.ntot ) goto 800
                                            if ( lwrite ) print *,'also taking ',mo,yy,jens
                                            xx(1,j) = series(mo,yy,jens)
                                            xx(2,j) = covariate(mo,yy,jens)
                                            if ( lwrite ) print *,'xx(:,',j,') = ',xx(:,j)
                                        end if
                                    end if
                                end do
                            else
                                ! monthly (or daily) values, ndecor refers to the 
                                ! autocorrelation of these.
                                nleft = ndecor - 2
                                do mmm=mo+1,min(j2,mo+ndecor-1)
                                    mm = mmm
                                    call normon(mm,yr,yy,nperyear)
                                    if ( yy <= yrstop .and. j < ntot ) then
                                        if ( series(mm,yy,jens).lt.1e33 .and. &
                                        &    covariate(mm,yy,jens).lt.1e33 ) then
                                            j = j + 1
                                            if ( j.gt.ntot ) goto 800
                                            nleft = nleft - 1
                                            if ( lwrite ) print *,'also taking ',mm,yy,jens
                                            xx(1,j) = series(mm,yy,jens)
                                            xx(2,j) = covariate(mm,yy,jens)
                                            if ( lwrite ) print *,'xx(:,',j,') = ',xx(:,j)
                                        end if
                                    end if
                                end do
                                ! if there are months/days left after we hit the end of the season,
                                ! try going backwards
                                do mmm=mo-1,max(j1,mo-nleft)
                                    mm = mmm
                                    call normon(mm,yr,yy,nperyear)
                                    if ( yy >= yrstart .and. j < ntot ) then
                                        if ( series(mm,yy,jens).lt.1e33 .and. &
                                        &    covariate(mm,yy,jens).lt.1e33 ) then
                                            j = j + 1
                                            if ( j.gt.ntot ) goto 800
                                            nleft = nleft - 1
                                            if ( lwrite ) print *,'also taking ',mm,yy,jens
                                            xx(1,j) = series(mm,yy,jens)
                                            xx(2,j) = covariate(mm,yy,jens)
                                            if ( lwrite ) print *,'xx(:,',j,') = ',xx(:,j)
                                        end if ! valid data
                                    end if ! valid year
                                end do ! loop over moving block
                            end if ! j1 == j2
                        end if ! ndecor > 0 
                    end if ! valid data
                end if ! crosscorr < cutoff
            end do ! ensemble members for spatial moving block
        end do ! while ( .not.ok )
    end do ! ntot, length of array
    800 continue
    do i=1,ntot
        if ( xx(1,i).gt.1e10.or. xx(2,i).gt.1e10 ) then
            write(0,*) 'sample_bootstrap: error: xx:,',i,') = ',xx(:,i)
        end if
    end do
    sdecor = real(ntot)/real(nens)
end subroutine sample_bootstrap

subroutine print_spatial_scale(scross)
    implicit none
    real scross
    if ( scross.gt.1.05 ) then
        print '(a,f6.1,a)','# Used spatial decorrelation blocks of ', & 
        &   scross,' stations on average in the bootstrap.'
    else
        print '(a)','# All series are considered independent in the bootstrap.'
    end if
end subroutine print_spatial_scale

subroutine handle_then_now(series,covariate,nperyear,fyr,lyr,j1,j2,yr1a,yr2a, &
    &   mens1,mens,xyear,ensmax,cov1,cov2,lprint,lwrite)

    ! handle the conversion from "then" (yr1a) and "now" (yr2a) to the variables
    ! fitgevcov expects (xyear,idmax,cov1,cov2), sets series to undef at "now".
    
    implicit none
    integer nperyear,fyr,lyr,j1,j2,yr1a,yr2a,mens1,mens,ensmax
    real series(nperyear,fyr:lyr,0:mens),covariate(nperyear,fyr:lyr,0:mens)
    real xyear,cov1,cov2
    logical lprint,lwrite
    real dummy

    call find_cov(series,covariate,nperyear,fyr,lyr,mens1,mens,j1,j2,yr1a,cov1,dummy,ensmax,1,lprint,lwrite)
    call find_cov(series,covariate,nperyear,fyr,lyr,mens1,mens,j1,j2,yr2a,cov2,xyear,ensmax,2,lprint,lwrite)

end subroutine handle_then_now

subroutine find_cov(series,covariate,nperyear,fyr,lyr,mens1,mens,j1,j2,yr,cov,xyear,ensmax,i12,lprint,lwrite)
    implicit none
    integer nperyear,fyr,lyr,mens1,mens,j1,j2,yr,ensmax,i12
    real series(nperyear,fyr:lyr,0:mens),covariate(nperyear,fyr:lyr,0:mens),cov,xyear
    logical lchangesign,lprint,lwrite
    integer i,j,mo,momax,yrmax,iens
    real s

    cov = 3e33
    ensmax = -1
    if ( yr.ne.9999 ) then
        s = -3e33
        do iens=mens1,mens
            do mo=j1,j2
                j = mo
                call normon(j,yr,i,nperyear)
                if ( i.ge.fyr .and. i.le.lyr ) then
                    if ( series(j,i,iens).lt.1e33 .and. series(j,i,iens).gt.s ) then
                        s = series(j,i,iens)
                        momax = j
                        yrmax = i
                        ensmax = iens
                    end if
                end if
            end do
        end do
    end if
    if ( abs(s).gt.1e33 ) then
        momax = 1
        yrmax = yr
        ensmax = mens1
        if ( yrmax.lt.fyr .or. yrmax.gt.lyr ) then
            if ( lprint ) then
                write(0,*) 'find_cov: error: yr ',yr,' outside range of data<br>'
            end if
        end if
    end if
    cov = covariate(momax,yrmax,ensmax)
    if ( cov.gt.1e33 ) then
        if ( lprint ) then
            write(0,*) '<p>find_cov: error: no valid value in covariate(', &
            &   momax,yrmax,ensmax,') = ',cov,'<br>'
        end if
        if ( lwrite ) then
            print *,'find_cov: error: no valid value in covariate(', &
            &   momax,yrmax,ensmax,') = ',cov
        end if
    end if
    if ( i12.eq.2 ) then
        if ( abs(xyear).gt.1e33 ) then
            xyear = series(momax,yrmax,ensmax)
        else
            ensmax = -1 ! xyear was given by user
        end if
        series(momax,yrmax,mens1:mens) = 3e33 ! for GPD we should also make a few values to the sides undef
        if ( lwrite ) print *,'find_cov: xyear = ',xyear,momax,yrmax,ensmax
        if ( xyear.gt.1e33 ) then
            if ( lprint ) then
                write(0,*) 'find_cov: error: cannot find valid data in ',yr,', periods ',j1,j2, &
                & ', ensemble members ',mens1,mens
            end if
            cov = 3e33
        end if
    end if
end subroutine find_cov

subroutine subtract_constant(covariate,series,nperyear,fyr,lyr,mens1,mens,cov1,cov2,offset,lwrite)

    ! subtract a constant from the covariate series and the reference points
    ! to keep numbers small
    
    implicit none
    integer nperyear,fyr,lyr,mens1,mens
    real covariate(nperyear,fyr:lyr,0:mens),series(nperyear,fyr:lyr,0:mens)
    real cov1,cov2,offset
    logical lwrite
    integer yr,mo,n,iens
    real s

    s = 0
    n = 0
    do iens=mens1,mens
        do yr=fyr,lyr
            do mo=1,nperyear
                if ( covariate(mo,yr,iens).lt.1e33 .and. series(mo,yr,iens).lt.1e33 ) then
                    n = n + 1
                    s = s + covariate(mo,yr,iens)
                end if
            end do
        end do
    end do
    if ( n.eq.0 ) return
    s = s/n
    do iens=mens1,mens
        do yr=fyr,lyr
            do mo=1,nperyear
                if ( covariate(mo,yr,iens).lt.1e33 ) then
                    covariate(mo,yr,iens) = covariate(mo,yr,iens) - s
                end if
            end do
        end do
    end do
    if ( lwrite ) print *,'subtract_constant: cov1,cov2 were ',cov1,cov2
    cov1 = cov1 - s
    cov2 = cov2 - s
    offset = s
    if ( lwrite ) print *,'subtract_constant: cov1,cov2 are  ',cov1,cov2,offset

end subroutine subtract_constant

subroutine decluster(xx,yrs,ntot,threshold,lwrite)
!
!   set all but the local maximum in a clustered maximum (t-tsep,t+tsep)
!   equal to a low value. tsep is determined as in Roth et al, 2014.
!
!   input:  xx(2,ntot)  values of time series (1,1:ntot) and covariate (2,1:ntot)
!           yrs(0:ntot) 10000*yr + 100*mo + dy
!   output: xx          with values adjusted so that only the maximum of each cluster remains
!
    implicit none
    integer nmax
    parameter(nmax=25)
    integer ntot,yrs(0:ntot)
    real xx(2,ntot),threshold
    logical lwrite
    integer i,j,m,n,nn,yr,mo,dy,jul1,jul2,tsep,jmax
    real p95,pcut,xmin,fracn(2:nmax),cutoff,s
    real,allocatable :: yy(:)
    integer,external :: julday
    
    cutoff = 0.05*0.002 ! number from Martin Roth, not in paper.
!
!   first obtain the 95th percentile
!
    allocate(yy(ntot))
    do i=1,ntot
        yy(i) = xx(1,i)
    end do
    pcut = 95
    call getcut(p95,pcut,ntot,yy)
    if ( lwrite ) then
        print *,'decluster: p95 = ',p95
    end if
    xmin = yy(1)
!
!   next the fraction of clusters with length >= n for which val>=p95
!
    fracn = 0
    do i=1,ntot
        if ( xx(1,i).gt.p95 ) then
            yr = yrs(i)/10000
            mo = mod(yrs(i),10000)/100
            dy = mod(yrs(i),100)
            jul1 = julday(mo,dy,yr)
            do n=2,nmax
                m = i + n - 1
                if ( m.le.ntot ) then
                    if ( xx(1,m).gt.p95 ) then
                        yr = yrs(m)/10000
                        mo = mod(yrs(m),10000)/100
                        dy = mod(yrs(m),100)
                        jul2 = julday(mo,dy,yr)
                        nn = jul2 - jul1 + 1
                        if ( nn.eq.n ) then ! no break in the time series...
                            fracn(n) = fracn(n) + 1
                        else
                            exit
                        end if
                    else
                        exit ! end of run of points exceeding p95
                    end if ! data(offset)>p95
                end if ! in range
            end do ! n
        end if ! data>p95
    end do ! i
    fracn = fracn/ntot
    if ( lwrite ) then
        do i=2,nmax
            print *,'decluster: fracn(',i,') = ',fracn(i)
        end do
    end if
!
!   the n for which the fraction is low enough defines tsep
!
    do n=2,nmax
        if ( fracn(n).lt.cutoff ) then
            tsep = n-2
            exit
        end if
    end do
    if ( lwrite ) then
        print *,'decluster: tsep,cutoff = ',tsep,cutoff
    end if
 !
 !  set xx(1,i) to xmin when it is not the maximum value in a cluster
 !
    if ( tsep.gt.0 ) then
        yy = xmin
        do i=1+tsep,ntot-tsep
            s = xx(1,i-tsep)
            jmax = -tsep
            do j=-tsep+1,tsep
                if ( xx(1,i+j).gt.s ) then
                    s = xx(1,i+j)
                    jmax = j
                end if
            end do
            if ( jmax.eq.0 ) then ! local maxmimum
                yy(i) = xx(1,i)
            end if
        end do
        do i=1,ntot
            xx(1,i) = yy(i)
        end do
!
!       adjust threshold if necessary
!
        call nrsort(ntot,yy)
        do i=1,ntot
            if ( yy(i).gt.xmin ) exit
        end do
        s = 100*(i+1)/real(ntot+1)
        if ( lwrite ) then
            print *,'last invalid value is yy(',max(1,i-1),') = ',yy(max(1,i-1))
            print *,'first valid value is yy(',i,') =  ',yy(i)
            print *,'corresponding to threshold =      ',s,'% of ',ntot,' points'
            print *,'compare to user threshold =       ',threshold,'%'
        end if
        if ( s.gt.threshold ) then
            write(0,*) 'decluster: adjusting threshold from ',threshold,' to ',s,'<br>'
            threshold = s
        end if
    end if
!
end
    
subroutine copyab3etc(a3,b3,xi3,alpha3,beta3,t3,tx3, &
     &           a,a25,a975,b,b975,xi,xi25,xi975,alpha,alpha25,alpha975, &
     &           beta,beta25,beta975,t,t25,t975,tx,tx25,tx975)
!
!   copy separate best fit and 95% lower and upper bounds to first dimension
!   of *3 variables
!
    implicit none
    real :: a3(3),b3(3),xi3(3),alpha3(3),beta3(3),t3(3,10,3),tx3(3,3)
    real :: a,a25,a975,b,b25,b975,xi,xi25,xi975,alpha,alpha25,alpha975, &
    &   beta,beta25,beta975,t(10,3),t25(10,3),t975(10,3),tx(3),tx25(3),tx975(3)
    integer i,j
    
    call copy3scalar(a3,a,a25,a975)
    if ( b.ge.0 ) then
        call copy3scalar(b3,b,b25,b975)
    else
        call copy3scalar(b3,-b,-b975,-b25)
    end if
    call copy3scalar(xi3,xi,xi25,xi975)
    call copy3scalar(alpha3,alpha,alpha25,alpha975)
    call copy3scalar(beta3,beta,beta25,beta975)
    do j=1,10
        do i=1,3
            call copy3scalar(t3(1,j,i),t(j,i),t25(j,i),t975(j,i))
        end do
    end do
    do i=1,3
        call copy3scalar(tx3(1,i),tx(i),tx25(i),tx975(i))
    end do
end subroutine

subroutine copy3scalar(a3,a,a25,a975)
    implicit none
    real a3(3),a,a25,a975
    a3(1) = a
    a3(2) = a25
    a3(3) = a975
end subroutine

subroutine normaliseseries(series,nperyear,fyr,lyr,mens1,mens,j1,j2,assume,lwrite)
!
!   normalise all ensemble time series in series to have the same mean
!
    implicit none
    integer nperyear,fyr,lyr,mens1,mens,j1,j2
    real series(nperyear,fyr:lyr,0:mens)
    logical lwrite
    character assume*(*)
    integer iens
    real mean,ensmean

    if ( mens1.eq.mens ) return

    if ( assume.eq.'shift' ) then
        ! compute overall mean
        call getmeanseries(series,nperyear,mens1,mens,fyr,lyr,j1,j2,.false.,ensmean)
        do iens=mens1,mens
            ! compute mean of this member
            call getmeanseries(series,nperyear,iens,iens,fyr,lyr,j1,j2,.false.,mean)
            ! shift series so that it has the same mean as the ensemble mean
            series(:,:,iens) = series(:,:,iens) - mean + ensmean
        end do
    else if ( assume.eq.'scale' ) then
        ! compute overall multiplicative mean
        call getmeanseries(series,nperyear,mens1,mens,fyr,lyr,j1,j2,.true.,ensmean)
        do iens=mens1,mens
            ! compute mean of this member
            call getmeanseries(series,nperyear,iens,iens,fyr,lyr,j1,j2,.true.,mean)
            ! scale series so that it has the same mean as the ensemble mean
            series(:,:,iens) = series(:,:,iens) / mean * ensmean
        end do 
    else if ( assume.eq.'both'  ) then
        write(0,*) 'normaliseseries: error: cannot normalise with assumption "both"'
        write(*,*) 'normaliseseries: error: cannot normalise with assumption "both"'
        call exit(-1)
    else
        write(0,*) 'normaliseseries: error: unknown assumption ',trim(assume)
        write(*,*) 'normaliseseries: error: unknown assumption ',trim(assume)
        call exit(-1)
    end if
end subroutine

subroutine checknonegative(series,nperyear,fyr,lyr,mens1,mens,j1,j2,assume, &
&   lchangesign,lwrite)
!
!   check that series does not contain negative values
!
    implicit none
    integer nperyear,fyr,lyr,mens1,mens,j1,j2
    real series(nperyear,fyr:lyr,0:mens)
    logical lchangesign,lwrite
    character assume*(*)
    integer iens,yy,yr,mm,mo

    if ( assume /= 'scale' ) return

    do iens=mens1,mens
        do yy=fyr,lyr
            do mm=j1,j2
                mo = mm
                call normon(mo,yy,yr,nperyear)
                if ( yr.ge.fyr .and. yr.le.lyr ) then
                    if ( .not.lchangesign .and. series(mo,yr,iens).lt.0 .or. &
                        & lchangesign .and. series(mo,yr,iens).gt.0 .and. &
                        &       series(mo,yr,iens).lt.1e33 ) then
                        write(0,*) 'error: option "scale" is not compatible '// &
                        &   'with negative values in the time series', &
                        &   series(mo,yr,iens)
                        write(*,*) 'error: option "scale" is not compatible '// &
                        &   'with negative values in the time series', &
                        &   series(mo,yr,iens)
                        call exit(-1)
                    end if
                end if
            end do
        end do
    end do
end subroutine

subroutine getmeanseries(series,nperyear,mens1,mens,fyr,lyr,j1,j2,lmult,mean)
    implicit none
    integer nperyear,fyr,lyr,mens1,mens,j1,j2
    real series(nperyear,fyr:lyr,0:mens),mean
    logical lmult
    integer iens,n,yr,yy,mo,mm
    real s

    n = 0
    s = 0
    do iens=mens1,mens
        do yy=fyr,lyr
            do mm=j1,j2
                mo = mm
                call normon(mo,yy,yr,nperyear)
                if ( yr.ge.fyr .and. yr.le.lyr ) then
                    if ( lmult ) then
                        if ( series(mo,yr,iens).lt.1e33 .and. series(mo,yr,iens).gt.0 ) then
                            n = n + 1
                            s = s + log(series(mo,yr,iens))
                        end if
                    else
                        if ( series(mo,yr,iens).lt.1e33 ) then
                            n = n + 1
                            s = s + series(mo,yr,iens)
                        end if
                    end if
                end if
            end do
        end do
    end do
    if ( lmult ) then
        mean = exp(s/n)
    else
        mean = s/n
    end if    
end subroutine
