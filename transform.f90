program transform
!
!   tranformation program belonging to the KNMI'14 scenarios
!   Based on the R versions by Alexander Bakker
!
    implicit none
    integer npermax,yrbeg,yrend
    parameter(npermax=366,yrbeg=1800,yrend=2200)
    integer horizon,nperyear
    real,allocatable :: series(:,:),newseries(:,:)
    real deltas(5,12)
    character variable*10,scenario*2,region*3,infile*255,var*40,units*80,scaling*5
    character string*100
    logical lwrite
    integer iargc
    lwrite = .false.

    if ( iargc().lt.5 ) then
        write(0,*) 'usage: transform infile variable horizon scenario[.subscenario] region'
        write(0,*) '       gives transformed time series following the KNMI''14 scenario'
        stop
    end if

    call get_parameters(infile,variable,horizon,scenario,scaling,region)
    allocate(series(npermax,yrbeg:yrend),newseries(npermax,yrbeg:yrend))
    call readseries(infile,series,npermax,yrbeg,yrend,nperyear, &
&       var,units,.false.,lwrite)
    if ( nperyear /= 366 ) then
        write(0,*) 'transform: error: only for daily data'
        write(*,*) 'transform: error: only for daily data'
        call abort
    end if
    call tolower(variable)
    call tolower(units)
    if ( variable == 'rr' ) then
        if ( units /= 'mm/dy' .and. units /= 'mm/day' .and. units /= 'mm/dag' ) then
            write(0,*) 'transform: error: rr units should be mm/dy, not ',trim(units)
            write(*,*) 'transform: error: rr units should be mm/dy, not ',trim(units)
            call abort
        end if
    else if ( variable == 'tg' .or. variable == 'tn' .or. variable == 'tx' ) then
        if ( units /= 'celsius' .and. units /= 'c' .and. units /= 'degrees_celsius' ) then
            write(0,*) 'transform: error: rr units should be celsius, not ',trim(units)
            write(*,*) 'transform: error: rr units should be celsius, not ',trim(units)
            call abort
        end if
    else
        write(0,*) 'transform: error: can only handle variables rr,tx,tg,tn at he moment, not ',trim(variable)
        write(*,*) 'transform: error: can only handle variables rr,tx,tg,tn at he moment, not ',trim(variable)
        call abort
    end if
    ! only the 30 years 1981-2010
    series(:,yrbeg:1980) = 3e33
    series(:,2011:yrend) = 3e33
    if ( variable == 'tg' .or. variable == 'tn' .or. variable == 'tx' ) then
        call read_deltas(variable,horizon,scenario,region,deltas,scaling,lwrite)
        call lintransform(series,npermax,yrbeg,yrend,horizon,deltas,newseries,lwrite)
    else if ( variable == 'rr' ) then
        call read_deltas(variable,horizon,scenario,region,deltas,scaling,lwrite)
        call rrtransform(series,npermax,yrbeg,yrend,horizon,deltas,newseries,lwrite)
    else
        write(0,*) 'transform: error: unknow value for varibale, expecting ', &
&           'tg,tn,tx or rr, found ',variable
        call abort
    end if
    call shiftseriesyear(newseries,npermax,nperyear,yrbeg,yrend,horizon-1995,lwrite)
    string = scenario
    if ( trim(variable) == 'rr' ) then
        string = trim(scenario)//'('//trim(scaling)//')'
    end if
    print '(3a)','# Transformed following the KNMI''14 scenario ',trim(string), ', v3.0'
    if ( trim(variable) /= 'rr' ) then
        print '(4a)','# for the region ',trim(region)
    end if
    call copyheader(infile,6)
    call transform_printdatfile(6,newseries,npermax,nperyear,yrbeg,yrend)

end program transform

subroutine get_parameters(infile,variable,horizon,scenario,scaling,region)
!
!   read parameters from the command line
!
    implicit none
    integer horizon
    character infile*(*), variable*(*), scenario*(*), region*(*), scaling*(*)
    character string*80

    call getarg(1,infile)
    call getarg(2,variable)
    call getarg(3,string)
    read(string,*,err=901) horizon
    call getarg(4,string)
    scenario = string(1:2)
    if ( string(3:3) == '.' ) then
        scaling = string(4:)
    else
        scaling = ' '
    end if
    if ( scenario.ne.'__' .and. &
    &   ( index('WG',scenario(1:1)).eq.0 .or. index('LH',scenario(2:2)).eq.0 ) ) then
        write(*,*) 'transform: error: expecting scenario __, GL, GH, WL or WH, not ',scenario
        write(0,*) 'transform: error: expecting scenario __, GL, GH, WL or WH, not ',scenario
        call abort
    end if
    call getarg(5,region)
    return
    
901 write(0,*) 'transform: error reading horizon from ',trim(string)
    call abort
end subroutine get_parameters

subroutine read_deltas(variable,horizon,scenario,region,deltas,scaling,lwrite)
!
!   read the deltas from file for the years 2030, 2050 and 2085
!   and interpolate to the requested horizon
!
    implicit none
    integer nregions,nyears
    parameter(nregions=7,nyears=3)
    integer horizon,yrs(nyears)
    real deltas(5,12)
    character variable*(*),scenario*(*),region*(*),scaling*(*)
    logical lwrite
    integer i,j,iyr,yr,iregion,mo,ndeltas,ip99
    character file*255,line*255,scen*2
    real delta0(5,12,nyears),a,aa(7)
    logical foundregion

    if ( scenario == '__' ) then
        yrs = (/1995,2030,2030/)
    else
        yrs = (/1995,2050,2085/)
    end if

    if ( variable == 'rr' ) then
        ndeltas = 7
        if ( scaling == 'lower' ) then
            ip99 = 5
        else if ( scaling == 'centr' ) then
            ip99 = 6
        else if ( scaling == 'upper' ) then
            ip99 = 7
        else
            write(0,*) 'read_deltas: error: expecting scaling lower|centr|upper, not ',scaling
            write(0,*) 'transform: error: expecting scaling lower|centr|upper, not ',scaling
            call abort
        end if
    else
        ndeltas = 5
        ip99 = 5
    end if

    delta0 = 3e33
    do iyr=1,nyears
        yr = yrs(iyr)
        if ( yr == 1995 ) then
            delta0(1:5,1:12,iyr) = 0
            cycle
        end if
        scen = scenario
        write(file,'(5a,i4,a)') 'KNMI14/deltas-KNMI14__',trim(variable),'___', &
&           trim(scen),'__',yr,'.txt'
        open(1,file=trim(file),status='old',err=900)
        read(1,'(a)') line
        foundregion = .false.
        do
            read(1,'(a)',end=100,err=901) line
            if ( variable == 'rr' ) then
                foundregion = .true.
                if ( trim(scen).eq.'__' ) then !^%$#$%^&*
                    mo = index('JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC',line(1:3))
                    if ( mo.eq.0 ) then
                        write(0,*) 'transform: error: cannot recognise month ',line(1:3)
                        call abort
                    end if
                    mo = (mo+2)/3
                    read(line(4:),*) (aa(j),j=1,ndeltas)
                else
                    read(line,*) mo,(aa(j),j=1,ndeltas)
                end if
                do j=1,5
                    delta0(j,mo,iyr) = aa(j)
                end do
                if ( ip99 /= 5 ) then
                    delta0(5,mo,iyr) = aa(ip99)
                end if
                if ( lwrite ) then
                    print *,'delta0(:,',mo,yrs(iyr),') = ',delta0(:,mo,iyr)
                end if
            else
                if ( line(1:3) == region ) then
                    foundregion = .true.
                    read(line(4:),*) mo,(aa(j),j=1,ndeltas)
                    do j=1,5
                        delta0(j,mo,iyr) = aa(j)
                    end do
                    if ( lwrite ) then
                        print *,'delta0(:,',mo,yrs(iyr),') = ',delta0(:,mo,iyr)
                    end if
                end if
            end if
        end do
100     continue
        if ( .not. foundregion ) then
            write(*,*) 'transform: error: cannot find region ',region
            write(0,*) 'transform: error: cannot find region ',region
            call abort
        end if
    end do ! iyr

    if ( horizon.lt.yrs(1) ) then
        write(0,*) 'get_deltas: error: horizon should be equal to or larger than ',yrs(1)
        write(*,*) 'transform: error: horizon should be equal to or larger than ',yrs(1)
        call abort
    end if
    if ( horizon.gt.yrs(nyears) ) then
        write(0,*) 'get_deltas: error: horizon should be equal to or smaller than ',yrs(nyears)
        write(*,*) 'transform: error: horizon should be equal to or smaller than ',yrs(nyears)
        call abort
    end if
    do iyr=1,nyears-1
        if ( horizon >= yrs(iyr) .and. horizon < yrs(iyr+1) .or. &
        &    iyr == nyears-1 .and. horizon == yrs(iyr+1) ) then  ! 2085 is special case
            ! linear interpolation
            a = real(horizon-yrs(iyr))/real(yrs(iyr+1)-yrs(iyr))
            if ( lwrite ) then
                print *,'linear interpolation ',(1-a),yrs(iyr),' and ',a,yrs(iyr+1)
            end if
            do mo=1,12
                do j=1,5
                    deltas(j,mo) = (1-a)*delta0(j,mo,iyr) + a*delta0(j,mo,iyr+1)
                end do ! j
                if ( lwrite ) then
                    print *,'deltas(:,',mo,') = ',deltas(:,mo)
                end if
            end do ! mo
        end if
    end do ! iyr
    return

900 write(0,*) 'read_deltas: error: cannot locate file ',trim(file)
    call abort
    
901 write(0,*) 'read_deltas: error reading frome file ',trim(file)
    write(0,*) 'read_deltas: around line ',trim(line)
    call abort
    
end subroutine read_deltas

subroutine lintransform(series,npermax,yrbeg,yrend,horizon,deltas,newseries,lwrite)
!
!   transform series over 1981-2010 into newseries over 30 years centered
!   on horizon using the deltas of the 5,50,95 percentiles
!   Based on transform_temperatuur.R by Alexander Bakker
!
    implicit none
    integer npermax,yrbeg,yrend,horizon
    real series(npermax,yrbeg:yrend),newseries(npermax,yrbeg:yrend),deltas(5,12)
    logical lwrite
    integer mo,d,dy,yr,ip
    real Xp(5,12),Yp(5,12),a,b
    
    newseries = 3e33
    do mo=1,12
        call get_quantiles(series,npermax,yrbeg,yrend,mo,Xp(1,mo))
    end do ! mo
    Yp = Xp + deltas
    do yr=max(1981,yrbeg),min(2010,yrend)
        do d=1,366
            call getdymo(dy,mo,d,366)
            do ip=2,4-1
                if ( series(d,yr).lt.Xp(ip,mo) ) exit
            end do
            a = (Yp(ip,mo)-Yp(ip-1,mo))/(Xp(ip,mo)-Xp(ip-1,mo))
            b = Yp(ip,mo) - a*Xp(ip,mo)
            newseries(d,yr+horizon-1995) = a*series(d,yr) + b
        end do ! dy
    end do ! yr
    
end subroutine lintransform

subroutine rrtransform(oldseries,npermax,yrbeg,yrend,horizon,deltas,series,lwrite)
    implicit none
    integer :: npermax,yrbeg,yrend,horizon
    real :: oldseries(npermax,yrbeg:yrend),series(npermax,yrbeg:yrend),deltas(5,12)
    logical :: lwrite
    integer :: i,j,n,mm,mo,nn(12),nntot(12),ndry,dwet,id,d,dy,yr,d1,yr1,im,idry,imin,nadd
    integer :: ntvalues,nprec,ntot,nmax
    real :: th,qq1,qq2,x1,x2,tol,a,b,c,diff,diffmin,stap
    real :: wdf_obs(12),mean_obs(12),mwet_obs(12),q2_obs(12),q1_obs(12),ratio(12)
    real :: wdf_fut(12),mean_fut(12),mwet_fut(12),q2_fut(12),q1_fut(12)
    integer,allocatable :: dyr(:,:,:),tmonths(:),ii(:),jj(:),itemp(:),add(:)
    real*8,allocatable :: X(:,:),Xw(:,:),tvalues(:),temp(:)
    real,allocatable :: Xm(:),X1m(:),precwet(:),aa(:)
    logical :: leftdry,rightdry,available
    character :: version*4
    real,external :: rr_transform_f,zbrent
    ! to the fit routine to determine b
    integer nXmax
    parameter(nXmax=31*30) ! 30 years, 31 days per month
    integer nXm
    real qfut,mfut,qobs,mobs,XXm(nXmax)
    common /crr_transform_f/ qfut,mfut,qobs,mobs,nXm,XXm

    th = 0.1
    qq1 = 0.99
    qq2 = 0.90
    version = 'v1.1'
    if ( lwrite ) print *,'rrtransform: th,qq1,qq1 = ',th,qq1,qq1
    
    open(1,file='KNMI14/ratio_Q99_Q90.txt',status='old')
    do mo=1,12
        read(1,*) mm,ratio(mo)
        if ( mm /= mo ) then
            write(0,*) 'rrtrans: error: mo /= mm ',mo,mm
            call abort
        end if
    end do
    
    series = oldseries ! transform in-place

    nmax = 31*30
    allocate(X(nmax,12),Xw(nmax,12),Xm(nmax),X1m(nmax),precwet(nmax),aa(nmax))
    allocate(dyr(2,nmax,12))
    allocate(tvalues(nmax))
    allocate(temp(nmax))
    allocate(tmonths(nmax))
    allocate(ii(nmax),itemp(nmax),add(nmax),jj(nmax))
    nn = 0
    nntot = 0
    mean_obs = 0
    mwet_obs = 0
    do yr=max(1981,yrbeg),min(2010,yrend)
        do d=1,366
            call getdymo(dy,mo,d,366)
            if ( series(d,yr).lt.1e33 ) then
                nntot(mo) = nntot(mo) + 1
                mean_obs(mo) = mean_obs(mo) + series(d,yr)
                if ( series(d,yr) >= th ) then
                    nn(mo) = nn(mo) + 1
                    X(nn(mo),mo) = series(d,yr)
                    mwet_obs(mo) = mwet_obs(mo) + series(d,yr)
                end if ! wet?
            end if ! valid?
        end do ! d        
    end do ! yr
    
    wdf_obs = real(nn)/real(nntot)
    mean_obs = mean_obs/nntot
    mwet_obs = mwet_obs/nn
    do mo=1,12
        ! I cannot use the climexp routine because it has different conventions from R
        Xm(1:nn(mo)) = X(1:nn(mo),mo)
        call nrsort(nn(mo),Xm)
        call quantile(q2_obs(mo),0.90,nn(mo),Xm,.false.)
    end do
    q1_obs = q2_obs*ratio
    if ( lwrite ) then
        print *,'wdf_obs = '
        print '(i2,f10.6)',(mo,wdf_obs(mo),mo=1,12)
        print *,'mean_obs = '
        print '(i2,f10.6)',(mo,mean_obs(mo),mo=1,12)
        print *,'mwet_obs = '
        print '(i2,f10.6)',(mo,mwet_obs(mo),mo=1,12)
        print *,'q2_obs = '
        print '(i2,f12.4)',(mo,q2_obs(mo),mo=1,12)
        print *,'q1_obs = '
        print '(i2,f12.4)',(mo,q1_obs(mo),mo=1,12)
    end if

    do mo=1,12
        wdf_fut(mo) = wdf_obs(mo)*(1+deltas(1,mo)/100)
        mean_fut(mo) = mean_obs(mo)*(1+deltas(2,mo)/100)
        mwet_fut(mo) = mean_fut(mo)/wdf_fut(mo)
        q1_fut(mo) = q1_obs(mo)*(1+deltas(5,mo)/100)
    end do ! mo
    if ( lwrite ) then
        print *,'wdf_fut = '
        print '(i2,f10.6)',(mo,wdf_fut(mo),mo=1,12)
        print *,'mean_fut = '
        print '(i2,f10.6)',(mo,mean_fut(mo),mo=1,12)
        print *,'mwet_fut = '
        print '(i2,f10.6)',(mo,mwet_fut(mo),mo=1,12)
        print *,'q1_fut = '
        print '(i2,f12.4)',(mo,q1_fut(mo),mo=1,12)
    end if

    ! drying wet days

    n = 0
    do im=1,12
        if ( lwrite ) print *,'deltas(1,',im,') = ',deltas(1,im)
        if ( deltas(1,im) < 0 ) n = n + 1
    end do
    if ( n.gt.0 ) then
        if ( version == 'v1.1' ) then
            ! select target values/days
            nn = 0
            nntot = 0
            do yr=max(1981,yrbeg),min(2010,yrend)
                do d=1,366
                    if ( series(d,yr) > 1e33 ) cycle
                    call getdymo(dy,mo,d,366)
                    nntot(mo) = nntot(mo) + 1
                    if ( nntot(mo) > nmax ) then
                        write(0,*) 'rrtransform: error: incraese nmax ',nmax,nntot(mo),mo
                        call abort
                    end if
                    dyr(1,nntot(mo),mo) = d
                    dyr(2,nntot(mo),mo) = yr
                    if ( series(d,yr) >= th ) then
                        ! make unique, the series are given to 1 decimal place 
                        X(nntot(mo),mo) = dble(0.1d0*nint(10*series(d,yr))) + 1d-10*(10000d0*yr + 100d0*mo + 1d0*dy)
                    else
                        X(nntot(mo),mo) = dble(series(d,yr))
                    end if
                    if ( deltas(1,mo) >= 0 ) cycle
                    if ( series(d,yr) >= th ) then
                        nn(mo) = nn(mo) + 1
                        if ( nn(mo) > nmax ) then
                            write(0,*) 'rrtransform: error: incraese nmax ',nmax,nn(mo),mo
                            call abort
                        end if
                        Xw(nn(mo),mo) = X(nntot(mo),mo)
                    end if
                end do ! d
            end do ! yr
            ntvalues = 0
            do mo=1,12
                if ( deltas(1,mo) >= 0 ) cycle
                call nrsort_d(nn(mo),Xw(1,mo))
                ndry = nint(-nn(mo)*deltas(1,mo)/100)
                if ( ndry > 0 ) then
                    stap = real(nn(mo))/ndry
                    if ( lwrite ) print *,'mo,ndry,stap = ',mo,ndry,stap
                    do i=1,ndry
                        ntvalues = ntvalues + 1
                        ! imitate the "round to even" algorithm of R
                        j = nint((i-0.5+0.0001*(mod(i,2)-0.5))*stap)
                        if ( j.gt.nn(mo) ) then
                            write(0,*) 'transform: error: j>nn(mo) ',j,nn(mo),mo
                            call abort
                        end if
                        tvalues(ntvalues) = Xw(j,mo)
                        tmonths(ntvalues) = mo
                    end do
                end if
            end do ! mo
            call ffsort_d(tvalues,ii,ntvalues)
            do i=1,ntvalues
                itemp(i) = tmonths(ii(i))
            end do
            tmonths = itemp
            do i=1,ntvalues
                temp(i) = tvalues(ii(i))
            end do
            tvalues = temp
            ! actual drying
            do idry=1,ntvalues ! select days for drying
                diffmin = 3e33
                mo = tmonths(idry)
                do i=1,nntot(mo)
                    d = dyr(1,i,mo)
                    yr = dyr(2,i,mo)

                    leftdry = .true.
                    d1 = d - 1
                    yr1 = yr
10                  continue
                    if ( d1 <= 0 ) then
                        d1 = d1 + 366
                        yr1 = yr1 - 1
                    endif
                    if ( yr1 >= max(1981,yrbeg) ) then
                        if ( series(d1,yr1).gt.1e33 ) then
                            d1 = d1 - 1
                            goto 10
                        end if
                        if ( series(d1,yr1) >= th ) leftdry = .false.
                    end if

                    rightdry = .true.
                    d1 = d + 1
                    yr1 = yr
20                  continue
                    if ( d1 > 366 ) then
                        d1 = d1 - 366
                        yr1 = yr1 + 1
                    endif
                    if ( yr1 <= min(2010,yrend) ) then
                        if ( series(d1,yr1).gt.1e33 ) then
                            d1 = d1 + 1
                            goto 20
                        end if
                        if ( series(d1,yr1) >= th ) rightdry = .false.
                    end if

                    if ( X(i,mo) > th .and. (leftdry .or. rightdry) ) then
                        available = .true.
                    else
                        available = .false.
                    end if
                    if ( available ) then
                        diff = abs(X(i,mo)-tvalues(idry))
                        if ( diff.lt.diffmin ) then
                            diffmin = diff
                            imin = i
                        end if
                    end if
                end do
                if ( lwrite ) print *,'series(',dyr(1,imin,mo),dyr(2,imin,mo),') = ', &
                &   X(imin,mo),' set to zero'
                X(imin,mo) = 0
                series(dyr(1,imin,mo),dyr(2,imin,mo)) = 0
            end do
        else if ( version == 'v1.2' ) then
            write(0,*) 'rrtransform: error: not yet ready ',version
            write(*,*) 'transform: error: not yet ready ',version
            call abort    
        else
            write(0,*) 'rrtransform: error: unknown version ',version
            write(*,*) 'transform: error: unknown version ',version
            call abort
        end if
    end if ! there are months that become drier
    
    ! wetting dry days
    
    do im=1,12
        if ( deltas(1,im) > 0 ) then
            n = 0
            ntot = 0
            do yr=max(1981,yrbeg),min(2010,yrend)
                do d=1,366
                    if ( series(d,yr).gt.1e33 ) cycle
                    call getdymo(dy,mo,d,366)
                    if ( mo == im ) then
                        ntot = ntot + 1
                        dyr(1,ntot,mo) = d
                        dyr(2,ntot,mo) = yr
                        Xm(ntot) = oldseries(d,yr)
                        d1 = d - 1
                        yr1 = yr
100                     continue
                        if ( d1 < 1 ) then
                            d1 = d1 + 366
                            yr1 = yr1 - 1
                        end if
                        if ( yr1 < max(1981,yrbeg) ) then
                            X1m(ntot) = 1
                        else
                            if ( oldseries(d1,yr1).gt.1e33 ) then
                                ! Feb 29 on non-leap days, or just undef in a user-supplied series
                                ! try again
                                d1 = d1 - 1
                                goto 100
                            end if
                            X1m(ntot) = oldseries(d1,yr1)
                        end if
                        if ( Xm(ntot) >= th ) then
                            n = n + 1
                            Xw(n,1) = Xm(ntot)
                        end if
                    end if ! mo == im
                end do ! d
            end do ! yr
            call nrsort_d(n,Xw)                    
            dwet = nint(deltas(1,im)/100*n) ! not yet round to even
            if ( lwrite ) print *,'im,deltas,n,dwet = ',im,deltas(1,im),n,dwet
            if ( dwet > 0 ) then
                ! select target values
                stap = real(n)/dwet ! I think it is real. Checked, indeed.
                if ( lwrite ) print *,'stap = ',stap
                do i=1,dwet
                    tvalues(i) = Xw(nint((i-0.5+0.0001*(mod(i,2)-0.5))*stap),1)
                end do
                if ( lwrite ) then
                    print *,'tvalues = '
                    print *,tvalues(1:dwet)
                end if
                ! select days to wet
                nprec = 0
                do i=1,ntot
                    if ( Xm(i) >= 0.1 ) then
                        nprec = nprec + 1
                    end if
                    precwet(i) = nprec
                end do
                precwet(1:ntot) = precwet(1:ntot) + stap/2
                nadd = 0
                do id=1,dwet
                    do i=1,ntot
                        if ( Xm(i) < th .and. X1m(i) >= 0.1 .and. precwet(i) >= stap ) then
                            if ( Xm(i).lt.1e33 ) then
                                nadd = nadd + 1
                                add(nadd) = i
                                precwet = precwet - stap
                                precwet(1:i) = 0
                                !!!exit ! break out of the i-loop (not needed in Fortran)
                            end if
                        end if
                    end do ! i=1,ntot
                end do ! id=1,dwet
                do i=1,nadd
                    aa(i) = X1m(add(i))
                end do
                call ffsort(aa,ii,nadd)
                do i=1,nadd
                    jj(ii(i)) = i
                end do
                if ( lwrite ) then
                    print *,'X1m(add(i)) = '
                    print *,aa(1:nadd)
                    print *,'rank = '
                    print *,jj(1:nadd)
                end if
                do i=1,nadd
                    series(dyr(1,add(i),im),dyr(2,add(i),im)) = tvalues(jj(i))
                end do
                if ( lwrite ) then
                    print *,'was',(oldseries(dyr(1,add(i),im),dyr(2,add(i),im)),i=1,nadd)
                    print *,' is',tvalues(jj(i:nadd))
                end if
            end if ! dfwet > 0
        end if ! days need to be added
    end do ! calendar month im
    
    ! transforming wet days
    ! determine coefficients
    do im=1,12
        nXm = 0
        do yr=max(1981,yrbeg),min(2010,yrend)
            do d=1,366
                if ( series(d,yr).gt.1e33 ) cycle
                call getdymo(dy,mo,d,366)
                if ( mo.ne.im ) cycle
                if ( series(d,yr) < th ) cycle
                nXm = nXm + 1
                XXm(nXm) = series(d,yr)
            end do  
        end do
        mobs = mwet_obs(im)
        qobs = q1_obs(im)
        mfut = mwet_fut(im)
        qfut = q1_fut(im)
        if ( lwrite ) print *,'mobs,qobs,mfut,qfut = ',mobs,qobs,mfut,qfut

        x1 = 0.1
        x2 = 3.
        tol = 0.0001
        if ( rr_transform_f(x1)*rr_transform_f(x2) < 0 ) then
            b = zbrent(rr_transform_f,x1,x2,tol)
            if ( lwrite ) print *,'b via zbrent = ',b
        else
            b = 3e33
            do i=1,300
                x1 = real(i)/100
                b = min(b,rr_transform_f(x1))
            end do
            if ( lwrite ) print *,'b via loop = ',b
        end if
        a = qfut/(qobs**b)
        c = a*(qobs**b)/qobs ! factor voor waarden groter dan q99
        if ( lwrite ) print *,'a,c = ',a,c

        ! transformeer
    
        do yr=max(1981,yrbeg),min(2010,yrend)
            do d=1,366
                if ( series(d,yr).gt.1e33 ) cycle
                call getdymo(dy,mo,d,366)
                if ( mo.ne.im ) cycle
                if ( series(d,yr) < th ) cycle
                
                if ( series(d,yr) < qobs ) then
                    series(d,yr) = a*series(d,yr)**b
                else
                    series(d,yr) = c*series(d,yr)
                end if
                if ( series(d,yr).lt.th ) series(d,yr) = th ! prevent days being dried by the transformation
            end do ! d
        end do ! yr
    
    end do ! im
    
end subroutine rrtransform

real function rr_transform_f(b)
    implicit none
    real b
    integer i
    real mean

    integer nXmax
    parameter(nXmax=31*30) ! 30 years, 31 days per month
    integer nXm
    real qfut,mfut,qobs,mobs,XXm(nXmax)
    common /crr_transform_f/ qfut,mfut,qobs,mobs,nXm,XXm
    
    mean = 0
    do i=1,nXm
        if ( XXm(i) < qobs ) then
            mean = mean + XXm(i)**b
        else
            mean = mean + XXm(i)*qobs**b/qobs
        end if
    end do
    mean = mean/nXm
    rr_transform_f = qfut/mfut - (qobs**b)/mean
end function

subroutine get_quantiles(series,npermax,yrbeg,yrend,mo,Xp)
!
!   compute the 1%,5%,50%,95%,99% quantiles of series in month mo in 1981-2010
!
    implicit none
    integer npermax,yrbeg,yrend,mo
    real series(npermax,yrbeg:yrend),Xp(5)
    integer :: yr,dy,i,m,n,dpm,offset
    real :: pcut(5)
    real,allocatable :: a(:)
    logical lwrite
    integer,external :: get_dpm
    lwrite = .false.

    pcut = (/ 1., 5., 50., 95., 99. /)
    allocate(a(30*31))
    offset = 0
    do m=1,mo-1
        offset = offset + get_dpm(m,2000)
    end do ! m
    n = 0
    do yr=max(1981,yrbeg),min(2010,yrend)
        dpm = get_dpm(mo,yr)
        do dy=offset+1,offset+dpm
            if ( series(dy,yr).lt.1e33 ) then
                n = n + 1
                a(n) = series(dy,yr)
            end if
        end do ! dy
    end do ! yr
    if ( n.lt.200 ) then ! totally arbitrary cut-off
        write(0,*) 'transform: error: need at least 200 days with valid data in month ', &
        & mo,' to compute quantiles'
        call abort
    end if
    call nrsort(n,a)
    do i=1,5
        call quantile(Xp(i),pcut(i)/100,n,a,lwrite)
    end do

end subroutine get_quantiles

integer function get_dpm(mo,yr)
    implicit none
    integer :: mo,yr
    integer,external :: leap
    if ( mo == 1 .or. mo == 3 .or. mo == 5 .or. mo == 7 .or. &
    &   mo == 8 .or. mo == 10 .or. mo == 12 ) then
        get_dpm = 31
    else if ( mo == 4 .or. mo == 6 .or. mo == 9 .or. mo == 11 ) then
        get_dpm = 30
    else
        if ( leap(yr) == 1 ) then
            get_dpm = 28
        else
            get_dpm = 29
        end if
    end if
end function get_dpm

subroutine quantile(q,cut,n,a,lwrite)
!
!   compute quantiles the same way as R, see http://stat.ethz.ch/R-manual/R-patched/library/stats/html/quantile.html
!   0<q<1
!   assume a has been sorted
!
    implicit none
    integer n
    real q,cut,a(n)
    logical lwrite
    integer i,type
    real x,eps
    eps = 1e-5
    type = 7 ! 7 is the R default, Geert used 5.

    if ( n.le.1 ) then
        q = 3e33
        return
    endif
    ! find elements around cut
    if ( type.eq.5 ) then
        ! matlab; hydrologists love this
        x = 0.5 + cut*n
    else if ( type.eq.6 ) then
        x = cut*(n+1)
    else if ( type.eq.7 ) then
        ! weird S/R definition...
        x = 1 + cut*(n-1)
    else
        write(0,*) 'quantile: error: type not yet supported: ',type
        call abort
    end if
    i = nint(x)
	if ( lwrite ) print *,'x,i,a(i),a(i+1) = ',x,i,a(max(i,n)),a(min(i-1,1))
    if ( abs(x-i).lt.eps ) then
        ! exact hit, demand unambiguous results
        if ( cut < 0.5 ) then
            q = a(i)*(1+eps)
        else
            q = a(i)*(1-eps)
        endif
    else
        ! interpolate
        i = int(x)
        x = x - i
        if ( i.lt.1 ) then
            q = (2-x)*a(1) + (x-1)*a(2)
        elseif ( i.gt.n-1 ) then
            q = -x*a(n-1) + (1+x)*a(n)
        else
            q = (1-x)*a(i) + x*a(i+1)
        endif
    endif
	if ( lwrite ) print *,'quantile: q = ',q
end subroutine

subroutine transform_printdatfile(unit,data,npermax,nperyear,yrbeg,yrend)
    implicit none
    integer unit,npermax,nperyear,yrbeg,yrend
    real data(npermax,yrbeg:yrend)
    integer year,i,dy,mo,dpm(12),ical
    data dpm /31,29,31,30,31,30,31,31,30,31,30,31/
    double precision val(360),offset
    call flush(unit)
    do year=yrbeg,yrend
        i = 0
        do mo=1,12
            do dy=1,dpm(mo)
                i = i + 1
                if ( data(i,year).lt.1e33 ) then
                    !!!write(unit,'(i6,2i2.2,f11.1)') year,mo,dy,data(i,year)*(1.00002)
                    write(unit,'(i6,2i2.2,f11.4)') year,mo,dy,data(i,year)
                end if
            end do
        end do
    end do
end subroutine

