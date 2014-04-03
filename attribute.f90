program attribute
!
!   do an empirical attribution study by fitting a time series to a
!   GEV or GPD with the position parameter linearly dependent on a covariate
!   and studying the difference in return time in the current climate and a
!   previous climate.
!
    implicit none
    include 'param.inc'
    include 'getopts.inc'
    integer nperyear,nperyear1,mens1,mens
    integer i,yr,mo,n,j1,j2
    real,allocatable :: series(:,:),covariate(:,:)
    character seriesfile*1024,covariatefile*1024,distribution*6,assume*5,string*80
    character var*40,units*80,var1*40,units1*80
    integer iargc

    if ( iargc().lt.8 ) then
        write(0,*) 'usage: attribute series covariate_series ', &
        & '(GEV)|(Gumbel)|(GPD gt threshold)|(Gauss) assume shift|scale', &
        & 'mon n [sel m] [ave N] [log|sqrt] ', &
        & 'begin2 past_climate_year end2 year_under_study', &
        & 'plot FAR_plot_file'
        write(0,*) 'note that n and m are in months even if the series is daily.'
        write(0,*) 'N is always in the same units as the series.'
        write(0,*) 'the covariate series is averaged to the same time scale.'
        stop
    end if
!
!   initialisation
!
!   getopts comes way too late for these options...
    do i=6,iargc()
        call getarg(i,string)
        if ( string.eq.'debug' .or. string.eq.'lwrite' ) then
            lwrite = .true.
        end if
        if ( string(1:9).eq.'standardu' ) then
            lstandardunits = .true.
        end if
    end do

    call getarg(3,distribution)
    call tolower(distribution)
    if ( distribution.ne.'gev' .and. distribution.ne.'gumbel' .and. &
    &    distribution.ne.'gpd' .and. distribution.ne.'gauss' ) then
        write(0,*) 'attribute: error: only GEV, GPD or Gauss supported, not ',distribution
        call abort
    end if

    call getarg(4,string)
    if ( string(1:4).ne.'assu' ) then
        write(0,*) 'attribute: error: expecting "assume", not ',trim(string)
        call abort
    end if
    call getarg(5,assume)
    call tolower(assume)
    if ( assume.ne.'shift' .and. assume.ne.'scale' .and. assume.ne.'both' ) then
        write(0,*) 'attribute: error: only shift, scale or both supported, not ',assume
        call abort
    end if

    call getarg(1,seriesfile)
    allocate(series(npermax,yrbeg:yrend))
    call readseries(seriesfile,series,npermax,yrbeg,yrend &
    & ,nperyear,var,units,lstandardunits,lwrite)
    
    call getarg(2,covariatefile)
    allocate(covariate(npermax,yrbeg:yrend))
    call readseries(covariatefile,covariate,npermax,yrbeg,yrend &
    & ,nperyear1,var1,units1,lstandardunits,lwrite)
    
    mens1 = 0
    mens = 0
    call getopts(6,iargc(),nperyear,yrbeg,yrend,.true.,mens1,mens)
    if ( yr1a.lt.yr1 .or. yr1a.lt.yrbeg ) then
        write(0,*) 'attribute: error: reference year should be after start of series ',yr1,yr1a
        call abort
    end if
    if ( yr2a.gt.yr2 .or. yr2a.gt.yrend ) then
        write(0,*) 'attribute: error: current year should be before end of series ',yr2,yr2a
        call abort
    end if
    call getj1j2(j1,j2,m1,nperyear,lwrite)
    call print_bootstrap_message(max(1,nint(decor)),j1,j2)
!
!   process data
!
    if ( lchangesign ) then
        call changesign(series,npermax,nperyear,yrbeg,yrend)
    endif
    if ( ldetrend ) then
        call detrend(series,npermax,nperyear,yrbeg,yrend,yr1,yr2,m1,m2,lsel)
    end if
    if ( anom ) then
        call anomal(series,npermax,nperyear,yrbeg,yrend,yr1,yr2)
    end if
    if ( lsum.gt.1 ) then
        call sumit(series,npermax,nperyear,yrbeg,yrend,lsum,oper)
    endif
    if ( logscale ) then
        print '(a)','# taking logarithm'
        call takelog(series,npermax,nperyear,yrbeg,yrend)
    endif
    if ( sqrtscale ) then
        print '(a)','# taking sqrt'
        call takesqrt(series,npermax,nperyear,yrbeg,yrend)
    endif
    if ( squarescale ) then
        print '(a)','# taking square'
        call takesquare(series,npermax,nperyear,yrbeg,yrend)
    endif

    call attribute_dist(series,nperyear,covariate,nperyear1,npermax,yrbeg,yrend,assume,distribution)

end program attribute

subroutine attribute_dist(series,nperyear,covariate,nperyear1,npermax,yrbeg,yrend,assume,distribution)
!
!   take block maxima, convert to linear arrays and call fitgevcov / fitgumcov or
!   take average, convert to linear arrays and call fitgaucov
!
    implicit none
    include "getopts.inc"
    integer nperyear,nperyear1,npermax,yrbeg,yrend
    real series(npermax,yrbeg:yrend),covariate(npermax,yrbeg:yrend)
    character assume*(*),distribution*(*)
    integer fyr,lyr,ntot,i,ntype,nmax,npernew,j1,j2
    integer,allocatable :: yrs(:)
    real a,b,xi,alpha,beta,xyear,cov1,cov2,offset
    real t(10,3),t25(10,3),t975(10,3),tx(3),tx25(3),tx975(3)
    real,allocatable :: xx(:,:),yrseries(:,:),yrcovariate(:,:),yy(:)
    logical lboot,lprint
    character operation*4,file*1024

    fyr = yr1
    lyr = yr2
    nmax = nperyear*(yr2-yr1+1)
    allocate(xx(2,nmax))
    allocate(yrs(0:nmax))
    allocate(yy(nmax))
    xx = 3e33
    
    if ( distribution.eq.'gev' .or. distribution.eq.'gumbel' ) then
        allocate(yrseries(1,fyr:lyr))
        allocate(yrcovariate(1,fyr:lyr))
        yrseries = 3e33
        yrcovariate = 3e33
        if ( lwrite ) print *,'attribute_dist: calling make_annual_series for series with max'
        call make_annual_values(series,nperyear,npermax,yrbeg,yrend,yrseries,fyr,lyr,'max')
        if ( lwrite ) print *,'attribute_dist: calling make_annual_series for covariate with mean'
        call make_annual_values(covariate,nperyear1,npermax,yrbeg,yrend,yrcovariate,fyr,lyr,'mean')
        npernew = 1
        j1 = 1
        j2 = 1
        m1 = 1
        lsel = 1
    else if ( distribution.eq.'gpd' .or. distribution.eq.'gauss' ) then
        allocate(yrseries(nperyear,fyr:lyr))
        allocate(yrcovariate(nperyear,fyr:lyr))
        yrseries = 3e33
        yrcovariate = 3e33
        ! copy series to keep the code easier to handle GEV and Gauss at the same time :-(
        yrseries(1:nperyear,fyr:lyr) = series(1:nperyear,fyr:lyr)
        ! change covariate to the same time resolution as series
        if ( nperyear1.le.nperyear ) then
            call annual2shorter(covariate,npermax,yrbeg,yrend,nperyear1, &
            &   yrcovariate,nperyear,fyr,lyr,nperyear,1,nperyear,1,lwrite)
        else if ( nperyear1.gt.nperyear ) then
            ! this should not occur in the web interface
            write(0,*) 'atribute_gev: error: covariate should not have higher time resolution than series: ', &
            &   nperyear1,nperyear
            write(*,*) 'atribute_gev: error: covariate should not have higher time resolution than series: ', &
            &   nperyear1,nperyear
            call abort
        else ! equal already
            yrcovariate(1:nperyear,fyr:lyr) = covariate(1:nperyear,fyr:lyr)
        end if
        npernew = nperyear
        call getj1j2(j1,j2,m1,npernew,.false.)
    else
        write(*,*) 'attribute_dist: error: unknown distribution ',trim(distribution)
        write(0,*) 'attribute_dist: error: unknown distribution ',trim(distribution)
        call abort
    end if    
    if ( lwrite ) print *,'attribute_dist: calling handle_then_now'
    call handle_then_now(yrseries,yrcovariate,npernew,fyr,lyr,j1,j2,yr1a,yr2a,xyear,cov1,cov2,lwrite)
    print '(a,i4,a,g16.5,a)','# <tr><td>covariate:</td><td>',yr1a,'</td><td>',cov1, &
    &   '</td><td>&nbsp;</td></tr>'
    print '(a,i4,a,g19.5,a)','# <tr><td>&nbsp;</td><td>',yr2a,'</td><td>',cov2, &
    &   '</td><td>&nbsp;</td></tr>'
    if ( lwrite ) print *,'attribute_dist: calling subtract_constant'
    call subtract_constant(yrcovariate,yrseries,npernew,fyr,lyr,cov1,cov2,offset,lwrite)
    if ( lwrite ) print *,'attribute_dist: calling fill_linear_array'
    call fill_linear_array(yrseries,yrcovariate,npernew,j1,j2,fyr,lyr,xx,yrs,nmax,ntot)
    if ( distribution.eq.'gpd' .and. nperyear.ge.360 ) then
        call decluster(xx,yrs,ntot,lwrite)
    end if

    if ( lweb ) then
        print '(8a)' ,'# <tr><th>parameter</th><th>year</th><th>value</th><th>' &
        &   ,'95% CI</th></tr>'
        print '(a,i9,a)','# <tr><td>N:</td><td>&nbsp;</td><td>',ntot, &
        &   '</td><td>&nbsp;</td></tr>'
    end if

    lboot = .true.
    lprint = .true.
    if ( distribution.eq.'gev' ) then
        ntype = 2 ! Gumbel plot
        if ( lwrite ) print *,'attribute_dist: calling fitgevcov'
        call fitgevcov(xx,yrs,ntot,a,b,xi,alpha,beta,j1,j2 &
    &       ,lweb,ntype,lchangesign,yr1a,yr2a,xyear,cov1,cov2,offset &
    &       ,t,t25,t975,tx,tx25,tx975,restrain,assume,lboot,lprint,plot,lwrite)
    else if ( distribution.eq.'gpd' ) then
        ntype = 3 ! log plot
        print *,'DEBUG'
        lboot = .false.
        lwrite = .true.
        if ( lwrite ) print *,'attribute_dist: calling fitgpdcov'
        call fitgpdcov(xx,yrs,ntot,a,b,xi,alpha,beta,j1,j2 &
    &       ,lweb,ntype,lchangesign,yr1a,yr2a,xyear,cov1,cov2,offset &
    &       ,t,t25,t975,tx,tx25,tx975,pmindata,restrain,assume,lboot,lprint,plot,lwrite)
    else if  ( distribution.eq.'gumbel' ) then
        ntype = 2 ! Gumbel plot
        if ( lwrite ) print *,'attribute_dist: calling fitgumcov'
        call fitgumcov(xx,yrs,ntot,a,b,alpha,beta,j1,j2 &
    &       ,lweb,ntype,lchangesign,yr1a,yr2a,xyear,cov1,cov2,offset &
    &       ,t,t25,t975,tx,tx25,tx975,assume,lboot,lprint,plot,lwrite)
    else if  ( distribution.eq.'gauss' ) then
        ntype = 4 ! sqrtlog plot
        if ( lwrite ) print *,'attribute_dist: calling fitgaucov'
        call fitgaucov(xx,yrs,ntot,a,b,alpha,beta,j1,j2 &
    &       ,lweb,ntype,lchangesign,yr1a,yr2a,xyear,cov1,cov2,offset &
    &       ,t,t25,t975,tx,tx25,tx975,assume,lboot,lprint,plot,lwrite)
    else
        write(0,*) 'attribute_dist: error: unknown distribution ',trim(distribution)
    end if

end subroutine attribute_dist

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

subroutine make_annual_values(series,nperyear,npermax,yrbeg,yrend,yrseries,fyr,lyr, &
& operation)
    
    ! construct an annual time series with the maxima

    implicit none
    include 'getopts.inc'
    integer nperyear,npermax,yrbeg,yrend,fyr,lyr
    real series(npermax,yrbeg:yrend),yrseries(1,fyr:lyr)
    character operation*(*)
    integer j1,j2,yy,yr,mm,mo,dd,dy,k,m,mtot,n,dpm(12)
    real s

    call getj1j2(j1,j2,m1,min(12,nperyear),lwrite)
    if ( lwrite ) print *,'make_annual_values: j1,j2,nperyear = ',j1,j2,nperyear
    call getdpm(dpm,nperyear)
    if ( lwrite ) print *,'                    dpm = ',dpm
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
                    if ( series(dy,yr).lt.1e33 ) then
                        m = m + 1
                        if ( operation.eq.'max' ) then
                            s = max(s,series(dy,yr))
                        else if ( operation.eq.'min' ) then
                            s = min(s,series(dy,yr))
                        else if ( operation.eq.'mean' ) then
                            s = s + series(dy,yr)
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
                yrseries(1,yy) = s
            else if ( operation.eq.'mean' ) then
                yrseries(1,yy) = s/m
            else
                write(0,*) 'make_annual_values: unknown operation ',trim(operation)
                call abort
            end if    
        end if
    end do

end subroutine make_annual_values

subroutine fill_linear_array(series,covariate,nperyear,j1,j2,fyr,lyr,xx,yrs,nmax,ntot)

    ! transfer the valid pairs in series, covariate to xx(1:2,1:ntot)
    
    implicit none
    integer nperyear,j1,j2,fyr,lyr,nmax,ntot
    integer yrs(0:nmax)
    real series(nperyear,fyr:lyr),covariate(nperyear,fyr:lyr),xx(2,nmax)
    integer yy,yr,mm,mo,day,month,yrstart,yrstop
    
    yrstart = lyr
    yrstop = fyr
    ntot = 0
    do yy=fyr,lyr
        do mm=j1,j2
            mo = mm
            call normon(mo,yy,yr,nperyear)
            if ( series(mo,yr).lt.1e33 .and. covariate(mo,yr).lt.1e33 ) then
                ntot = ntot + 1
                xx(1,ntot) = series(mo,yr)
                xx(2,ntot) = covariate(mo,yr)
                call getdymo(day,month,mo,nperyear)
                yrs(ntot) = 10000*yr + 100*month + day
                if ( nperyear.gt.1000 ) then
                    yrs(ntot) = 100*yrs(ntot) + mod(ntot,nint(nperyear/366.))
                end if
                yrstart = min(yrstart,yr)
                yrstop = max(yrstop,yr)
            end if
        end do
    end do

    call savestartstop(yrstart,yrstop)
end subroutine fill_linear_array

subroutine handle_then_now(series,covariate,nperyear,fyr,lyr,j1,j2,yr1a,yr2a, &
    &   xyear,cov1,cov2,lwrite)

    ! handle the conversion from "then" (yr1a) and "now" (yr2a) to the variables
    ! fitgevcov expects (xyear,cov1,cov2), sets series to undef at "now".
    
    implicit none
    integer nperyear,fyr,lyr,j1,j2,yr1a,yr2a
    real series(nperyear,fyr:lyr),covariate(nperyear,fyr:lyr),xyear,cov1,cov2
    logical lwrite
    real dummy

    call find_cov(series,covariate,nperyear,fyr,lyr,j1,j2,yr1a,cov1,dummy,1,lwrite)
    call find_cov(series,covariate,nperyear,fyr,lyr,j1,j2,yr2a,cov2,xyear,2,lwrite)

end subroutine handle_then_now

subroutine find_cov(series,covariate,nperyear,fyr,lyr,j1,j2,yr,cov,xyear,i12,lwrite)
    implicit none
    integer nperyear,fyr,lyr,j1,j2,yr,i12
    real series(nperyear,fyr:lyr),covariate(nperyear,fyr:lyr),cov,xyear
    logical lchangesign,lwrite
    integer mo,momax
    real s

    s = -3e33
    do mo=j1,j2
        if ( series(mo,yr).lt.1e33 .and. series(mo,yr).gt.s ) then
            s = series(mo,yr)
            momax = mo
        end if
    end do
    if ( abs(s).lt.1e33 ) then
        cov = covariate(momax,yr)
        if ( cov.gt.1e33 ) then
            write(0,*) 'find_cov: error: no valid value in cavariate(',momax,yr,') = ',cov
            call abort
        end if
        if ( i12.eq.2 ) then
            xyear = series(momax,yr)
            series(momax,yr) = 3e33 ! for GPD we should also make a few values to the sides undef
            if ( lwrite ) print *,'find_cov: xyear = ',xyear,momax,yr
        end if
    else
        write(0,*) 'find_cov: error: cannot find valid data in ',yr
        call abort
    end if

end subroutine find_cov

subroutine subtract_constant(covariate,series,nperyear,fyr,lyr,cov1,cov2,offset,lwrite)

    ! subtract a constant from the covariate series and the reference points
    ! to keep numbers small
    
    implicit none
    integer nperyear,fyr,lyr
    real covariate(nperyear,fyr:lyr),series(nperyear,fyr:lyr),cov1,cov2,offset
    logical lwrite
    integer yr,mo,n
    real s

    s = 0
    n = 0
    do yr=fyr,lyr
        do mo=1,nperyear
            if ( covariate(mo,yr).lt.1e33 .and. series(mo,yr).lt.1e33 ) then
                n = n + 1
                s = s + covariate(mo,yr)
            end if
        end do
    end do
    if ( n.eq.0 ) return
    s = s/n
    do yr=fyr,lyr
        do mo=1,nperyear
            if ( covariate(mo,yr).lt.1e33 ) then
                covariate(mo,yr) = covariate(mo,yr) - s
            end if
        end do
    end do
    if ( lwrite ) print *,'subtract_constant: cov1,cov2 were ',cov1,cov2
    cov1 = cov1 - s
    cov2 = cov2 - s
    offset = s
    if ( lwrite ) print *,'subtract_constant: cov1,cov2 are  ',cov1,cov2,offset

end subroutine subtract_constant

subroutine decluster(xx,yrs,ntot,lwrite)
!
!   set all but the local maximum in a clustered maximum (t-tsep,t_+tsep)
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
    real xx(2,ntot)
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
    xmin = min(yy(1),yy(ntot)) ! not sure which way yy has been sorted :-)
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
            if ( j.eq.0 ) then ! local maxmimum
                yy(i) = xx(1,i)
            end if
        end do
        do i=1,ntot
            xx(1,i) = yy(i)
        end do
    end if
!
    end
