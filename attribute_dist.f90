subroutine attribute_dist(series,nperyear,covariate,nperyear1,npermax,yrbeg,yrend,&
&   mens1,mens,assume,distribution)
!
!   take block maxima, convert to linear arrays and call fitgevcov / fitgumcov or
!   take average, convert to linear arrays and call fitgaucov / fitgpdcov
!
    implicit none
    include "getopts.inc"
    integer nperyear,nperyear1,npermax,yrbeg,yrend,mens1,mens
    real series(npermax,yrbeg:yrend,0:mens),covariate(npermax,yrbeg:yrend,0:mens)
    character assume*(*),distribution*(*)
    integer fyr,lyr,ntot,i,ntype,nmax,npernew,j1,j2,iens
    integer,allocatable :: yrs(:)
    real a,b,xi,alpha,beta,cov1,cov2,offset
    real t(10,3),t25(10,3),t975(10,3),tx(3),tx25(3),tx975(3)
    real,allocatable :: xx(:,:),yrseries(:,:,:),yrcovariate(:,:,:),yy(:)
    logical lboot,lprint,subtract_offset
    character operation*4,file*1024

    fyr = yr1
    lyr = yr2
    nmax = nperyear*(yr2-yr1+1)*(1+mens-mens1)
    allocate(xx(2,nmax))
    allocate(yrs(0:nmax))
    allocate(yy(nmax))
    xx = 3e33
    
    if ( distribution.eq.'gev' .or. distribution.eq.'gumbel' ) then
        allocate(yrseries(1,fyr:lyr,0:mens))
        allocate(yrcovariate(1,fyr:lyr,0:mens))
        yrseries = 3e33
        yrcovariate = 3e33
        if ( lwrite ) print *,'attribute_dist: calling make_annual_series for series with max'
        call make_annual_values(series,nperyear,npermax,yrbeg,yrend,mens1,mens, &
        &   yrseries,fyr,lyr,'max')
        if ( lwrite ) print *,'attribute_dist: calling make_annual_series for covariate with mean'
        call make_annual_values(covariate,nperyear1,npermax,yrbeg,yrend,mens1,mens, &
        &   yrcovariate,fyr,lyr,'mean')
        npernew = 1
        j1 = 1
        j2 = 1
        m1 = 1
        lsel = 1
    else if ( distribution.eq.'gpd' .or. distribution.eq.'gauss' ) then
        allocate(yrseries(nperyear,fyr:lyr,0:mens))
        allocate(yrcovariate(nperyear,fyr:lyr,0:mens))
        yrseries = 3e33
        yrcovariate = 3e33
        ! copy series to keep the code easier to handle GEV and Gauss at the same time :-(
        yrseries(1:nperyear,fyr:lyr,mens1:mens) = series(1:nperyear,fyr:lyr,mens1:mens)
        ! change covariate to the same time resolution as series
        if ( nperyear1.le.nperyear ) then
            do iens=mens1,mens
                call annual2shorter(covariate(1,yrbeg,iens),npermax,yrbeg,yrend,nperyear1, &
                &   yrcovariate(1,fyr,iens),nperyear,fyr,lyr,nperyear,1,nperyear,1,lwrite)
            end do
        else if ( nperyear1.gt.nperyear ) then
            ! this should not occur in the web interface
            write(0,*) 'atribute_gev: error: covariate should not have higher time resolution than series: ', &
            &   nperyear1,nperyear
            write(*,*) 'atribute_gev: error: covariate should not have higher time resolution than series: ', &
            &   nperyear1,nperyear
            call abort
        else ! equal already
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
    if ( lwrite ) print *,'attribute_dist: calling handle_then_now'
    call handle_then_now(yrseries,yrcovariate,npernew,fyr,lyr,j1,j2,yr1a,yr2a,mens1,mens,xyear,cov1,cov2,lwrite)
    print '(a,i4,a,g16.5,a)','# <tr><td>covariate:</td><td>',yr1a,'</td><td>',cov1, &
    &   '</td><td>&nbsp;</td></tr>'
    print '(a,i4,a,g19.5,a)','# <tr><td>&nbsp;</td><td>',yr2a,'</td><td>',cov2, &
    &   '</td><td>&nbsp;</td></tr>'
    subtract_offset = .true.
    if ( subtract_offset ) then
        if ( lwrite ) print *,'attribute_dist: calling subtract_constant'
        call subtract_constant(yrcovariate,yrseries,npernew,fyr,lyr,mens1,mens, &
        &   cov1,cov2,offset,lwrite)
    else
        offset = 0
    end if
    if ( lwrite ) print *,'attribute_dist: calling fill_linear_array'
    call fill_linear_array(yrseries,yrcovariate,npernew,j1,j2,fyr,lyr,mens1,mens, &
    &   xx,yrs,nmax,ntot)
    if ( distribution.eq.'gpd' .and. nperyear.ge.360 ) then
        call decluster(xx,yrs,ntot,pmindata,lwrite)
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
        !!!print *,'DEBUG'
        !!!lboot = .false.
        !!!lwrite = .true.
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

subroutine fill_linear_array(series,covariate,nperyear,j1,j2,fyr,lyr,mens1,mens,&
&   xx,yrs,nmax,ntot)

    ! transfer the valid pairs in series, covariate to xx(1:2,1:ntot)
    
    implicit none
    integer nperyear,j1,j2,fyr,lyr,mens1,mens,nmax,ntot
    integer yrs(0:nmax)
    real series(nperyear,fyr:lyr,0:mens),covariate(nperyear,fyr:lyr,0:mens),xx(2,nmax)
    integer yy,yr,mm,mo,day,month,yrstart,yrstop,iens

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

subroutine handle_then_now(series,covariate,nperyear,fyr,lyr,j1,j2,yr1a,yr2a, &
    &   mens1,mens,xyear,cov1,cov2,lwrite)

    ! handle the conversion from "then" (yr1a) and "now" (yr2a) to the variables
    ! fitgevcov expects (xyear,cov1,cov2), sets series to undef at "now".
    
    implicit none
    integer nperyear,fyr,lyr,j1,j2,yr1a,yr2a,mens1,mens
    real series(nperyear,fyr:lyr,0:mens),covariate(nperyear,fyr:lyr,0:mens)
    real xyear,cov1,cov2
    logical lwrite
    real dummy

    call find_cov(series,covariate,nperyear,fyr,lyr,mens1,mens,j1,j2,yr1a,cov1,dummy,1,lwrite)
    call find_cov(series,covariate,nperyear,fyr,lyr,mens1,mens,j1,j2,yr2a,cov2,xyear,2,lwrite)

end subroutine handle_then_now

subroutine find_cov(series,covariate,nperyear,fyr,lyr,mens1,mens,j1,j2,yr,cov,xyear,i12,lwrite)
    implicit none
    integer nperyear,fyr,lyr,mens1,mens,j1,j2,yr,i12
    real series(nperyear,fyr:lyr,0:mens),covariate(nperyear,fyr:lyr,0:mens),cov,xyear
    logical lchangesign,lwrite
    integer i,j,mo,momax,yrmax,ensmax,iens
    real s

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
        if ( yrmax.lt.fyr .or. yrmax.gt.lyr ) then
            write(0,*) 'find_cov: error: yr ',yr,' outside range of data'
                call abort
        end if
        ensmax = mens1
    end if
    cov = covariate(momax,yrmax,ensmax)
    if ( cov.gt.1e33 ) then
        write(0,*) 'find_cov: error: no valid value in cavariate(', &
        &   momax,yrmax,iens,') = ',cov
        call abort
    end if
    if ( i12.eq.2 ) then
        if ( xyear.gt.1e33 ) xyear = series(momax,yrmax,ensmax)
        series(momax,yrmax,mens1:mens) = 3e33 ! for GPD we should also make a few values to the sides undef
        if ( lwrite ) print *,'find_cov: xyear = ',xyear,momax,yrmax
    end if
    if ( xyear.gt.1e33 ) then
        write(0,*) 'find_cov: error: cannot find valid data in ',yr,', periods ',j1,j2, &
        & ', ensemble members ',mens1,mens
        call abort
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
