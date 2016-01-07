program autocor
!
!   correlate autocorrelation function of a time series
!
!   GJvO hack, jan-1998
!
    implicit none
    include "getopts.inc"
    include "param.inc"
    integer lagmax
    parameter (lagmax=100)
    integer i,ii,j,jj,lag,n,nperyear,iens,mens1,mens,yr,mo,yrstart,yrstop,lagunit
    real data(npermax,yrbeg:yrend,0:nensmax),s,adata,cc(0:lagmax),absent,sig(0:lagmax),norlag,lagfac
    parameter (absent=3e33)
    character line*128,var*40,units*20
    integer iargc
!
!   check arguments
!
    lwrite = .false.
    n = iargc()
    if ( n.lt.1 ) then
        print *,'usage: autocor datafile [log] [month n] [sum n] [anom] [decor n]'
        stop
    endif
!
!       read data
!
    lstandardunits = .false.
    call getarg(1,line)
    call readensseries(line,data,npermax,yrbeg,yrend,nensmax &
 &       ,nperyear,mens1,mens,var,units,lstandardunits,lwrite)
!
!   arguments
!
    call getopts(2,iargc(),nperyear,yrbeg,yrend,.true.,mens1,mens)
    yrstart = yr2
    yrstop  = yr1
    if ( mens.gt.0 ) then
        write(0,*) 'Concatenating ensemble members',nens1,' to ',nens2,'<br>'
        print '(a,i4,a,i4)','# taking ensemble members ',nens1,' to ',nens2
    endif
!
!   sum
!
    if ( lsum.gt.1 ) then
        do iens=nens1,nens2
            call sumit(data(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend,lsum,'+')
        enddo
    endif
!
!   log
!
    if ( logscale ) then
        do iens=nens1,nens2
            call takelog(data(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend)
        enddo
    endif
!
!   sqrt
!
    if ( sqrtscale ) then
        do iens=nens1,nens2
            call takesqrt(data(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend)
        enddo
    endif
!
!   square
!
    if ( squarescale ) then
        do iens=nens1,nens2
            call takesquare(data(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend)
        enddo
    endif
!
!       detrend
!
    if ( ldetrend ) then
        do iens=nens1,nens2
            call detrend(data(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend,yr1,yr2,m1,m2,1)
        enddo
    endif
!
!   anomalies
!
    if ( anom ) then
        do iens=nens1,nens2
            call anomal(data(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend,yr1,yr2,yrend)
        enddo
    endif
    if ( lensanom .and. min(nens2,mens).gt.max(nens1,mens1) ) then
        if ( lwrite ) print *,'Taking anomalies wrt ensemble mean'
        call anomalensemble(data,npermax,nperyear,yrbeg,yrend, &
 &           yr1,yr2,max(nens1,mens1),min(nens2,mens))
    endif
!
!   autocor!
!
    call autocov(data,npermax,yrbeg,yrend,nperyear,lagmax,cc,sig,yrstart,yrstop,lagunit)

    print '(a)','# lag     autocor     autocov   2/sqrt(N)'
    if ( nperyear/lagunit.eq.1 ) then
        print '(a)','# lag in yr'
        lagfac = 1
    elseif ( nperyear/lagunit.eq.4 ) then 
        print '(a)','# lag in yr'
        lagfac = 1/4.
    elseif ( nperyear/lagunit.eq.12 ) then
        print '(a)','# lag in mo'
        lagfac = 1
    elseif ( nperyear/lagunit.eq.36 ) then
        print '(a)','# lag in mo'
        lagfac = 1/3.
    elseif ( nperyear/lagunit.eq.48 ) then
        print '(a)','# lag in mo'
        lagfac = 1/4.
    elseif ( nperyear/lagunit.eq.52 ) then
        print '(a)','# lag in dy'
        lagfac = 7
    elseif ( nperyear/lagunit.eq.73 ) then
        print '(a)','# lag in dy'
        lagfac = 5
    elseif ( nperyear/lagunit.ge.360 .and. &
 &           nperyear/lagunit.le.366 ) then
        print '(a)','# lag in dy'
        lagfac = 1
    else
        write(0,*) 'autocor: error: cannot handle ',nperyear,' data points per year yet'
        call abort
    endif
    do lag=0,lagmax
        norlag = lagfac*lag
        if ( cc(lag).lt.1e33 ) then
            if ( lag.ge.decor .and. sig(lag).le.1 ) then
                print '(f7.1,f12.4,g16.4,f12.4)',norlag,cc(lag)/cc(0),cc(lag),sig(lag)
            else
                print '(f7.1,f12.4,g16.4,f12.4)',norlag,cc(lag)/cc(0),cc(lag)
            endif
        endif
    enddo
    call savestartstop(yrstart,yrstop)
end program
