program autocor
!
!   correlate autocorrelation function of a time series
!
!   GJvO hack, jan-1998
!
    implicit none
    include 'getopts.inc'
    include 'param.inc'
    integer,parameter :: lagmax=100
    real,parameter :: absent=3e33
    integer :: i,ii,j,jj,lag,n,nperyear,iens,mens1,mens,yr,mo,yrstart,yrstop,lagunit
    real :: s,adata,cc(0:lagmax),sig(0:lagmax),norlag,lagfac
    real,allocatable :: data(:,:,:)
    character :: line*128,file*1024,var*80,units*80,lvar*120,svar*120,history*50000, &
        metadata(2,100)*2000
!
!   check arguments
!
    lwrite = .false.
    n = command_argument_count()
    if ( n.lt.1 ) then
        print *,'usage: autocor datafile [log] [month n] [sum n] [anom] [decor n]'
        stop
    end if
!
!   read data
!
    lstandardunits = .false.
    allocate(data(npermax,yrbeg:yrend,0:nensmax))
    call get_command_argument(1,file)
    call readensseriesmeta(file,data,npermax,yrbeg,yrend,nensmax &
        ,nperyear,mens1,mens,var,units,lvar,svar,history,metadata,lstandardunits,lwrite)
!
!   arguments
!
    call getopts(2,command_argument_count(),nperyear,yrbeg,yrend,.true.,mens1,mens)
    yrstart = yr2
    yrstop  = yr1
    if ( mens.gt.0 ) then
        write(0,*) 'Concatenating ensemble members',nens1,' to ',nens2,'<br>'
        print '(a,i4,a,i4)','# taking ensemble members ',nens1,' to ',nens2
    end if
!
!   sum
!
    if ( lsum.gt.1 ) then
        do iens=nens1,nens2
            call sumit(data(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend,lsum,'+')
        end do
    end if
!
!   log
!
    if ( logscale ) then
        do iens=nens1,nens2
            call takelog(data(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend)
        end do
    end if
!
!   sqrt
!
    if ( sqrtscale ) then
        do iens=nens1,nens2
            call takesqrt(data(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend)
        end do
    end if
!
!   square
!
    if ( squarescale ) then
        do iens=nens1,nens2
            call takesquare(data(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend)
        end do
    end if
!
!       detrend
!
    if ( ldetrend ) then
        do iens=nens1,nens2
            call detrend(data(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend,yr1,yr2,m1,m2,1)
        end do
    end if
!
!   anomalies
!
    if ( anom ) then
        do iens=nens1,nens2
            call anomal(data(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend,yr1,yr2,yrend)
        end do
    end if
    if ( lensanom .and. min(nens2,mens).gt.max(nens1,mens1) ) then
        if ( lwrite ) print *,'Taking anomalies wrt ensemble mean'
        call anomalensemble(data,npermax,nperyear,yrbeg,yrend, &
            yr1,yr2,max(nens1,mens1),min(nens2,mens))
    end if
!
!   autocor!
!
    call autocov(data,npermax,yrbeg,yrend,nperyear,lagmax,cc,sig,yrstart,yrstop,lagunit)


    call printmetadata(6,trim(file),' ','autocorrelation of ',history,metadata)
    print '(a)','# lag     autocor     autocov   2/sqrt(N)'
    if ( nperyear/lagunit.eq.1 ) then
        print '(a)','# lag in yr'
        lagfac = 1
    else if ( nperyear/lagunit.eq.4 ) then 
        print '(a)','# lag in yr'
        lagfac = 1/4.
    else if ( nperyear/lagunit.eq.12 ) then
        print '(a)','# lag in mo'
        lagfac = 1
    else if ( nperyear/lagunit.eq.36 ) then
        print '(a)','# lag in mo'
        lagfac = 1/3.
    else if ( nperyear/lagunit.eq.48 ) then
        print '(a)','# lag in mo'
        lagfac = 1/4.
    else if ( nperyear/lagunit.eq.52 ) then
        print '(a)','# lag in dy'
        lagfac = 7
    else if ( nperyear/lagunit.eq.73 ) then
        print '(a)','# lag in dy'
        lagfac = 5
    else if ( nperyear/lagunit.ge.360 .and. &
             nperyear/lagunit.le.366 ) then
        print '(a)','# lag in dy'
        lagfac = 1
    else
        write(0,*) 'autocor: error: cannot handle ',nperyear,' data points per year yet'
        call exit(-1)
    end if
    do lag=0,lagmax
        norlag = lagfac*lag
        if ( cc(lag).lt.1e33 ) then
            if ( lag.ge.decor .and. sig(lag).le.1 ) then
                print '(f7.1,f12.4,g16.4,f12.4)',norlag,cc(lag)/cc(0),cc(lag),sig(lag)
            else
                print '(f7.1,f12.4,g16.4,f12.4)',norlag,cc(lag)/cc(0),cc(lag)
            end if
        end if
    end do
    call savestartstop(yrstart,yrstop)
end program
