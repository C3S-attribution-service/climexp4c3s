program extractseries
!
!       Convert a time series in my standard format to a table for R
!       (c) Geert Jan van Oldenborgh, 7-jun-2006, KNMI
!       This file may be copied and modified  freely as long as the
!       above copyright notice is kept intact
!
    implicit none
    include 'param.inc'
    include 'getopts.inc'
    integer :: nperyear,mens1,mens,mo,yr,i,j,j1,j2,mo1,month,iu,iens &
        ,nperyear2,yrmin,yrmax,n,nmodel1,nmodel2,ndup(npermax) &
        ,validens(npermax)
    real,allocatable :: fcst(:,:,:)
    real :: s1,s2
    character :: line*255,varfcst*20,unitfcst*40

!   check arguments

    lwrite = .false. 
    lstandardunits = .true. 
    if ( command_argument_count() < 4 ) then
        print *,'usage: extractseries infile '// &
            '[month m[:n] [sum|ave|max|min|sel m] '// &
            '[begin yr] [end yr] [detrend]'// &
            '[debias none|mean|var|all] dump outtable'
        call exit(-1)
    endif

    call get_command_argument(1,line)
    allocate(fcst(npermax,yrbeg:yrend,0:nensmax))
    call readensseries(line,fcst,npermax,yrbeg,yrend,nensmax &
        ,nperyear,mens1,mens,varfcst,unitfcst,lstandardunits &
        ,lwrite)
    print '(a,i2,a,i2)','# located ensemble members ',mens1,' to ',mens

    call getopts(2,command_argument_count(),nperyear,yrbeg,yrend, .true. ,mens1,mens)
    if ( .not. dump ) then
        write(0,*) 'extractseries: error: cannot find ''dump'' in argument list'
        write(*,*) 'extractseries: error: cannot find ''dump'' in argument list'
        call exit(-1)
    endif
    if ( lag1 /= 0 .or. lag2 /= 0 ) then
        write(0,*) 'extractseries: error: cannot handle lags ',lag1,lag2
        write(*,*) 'extractseries: error: cannot handle lags ',lag1,lag2
        call exit(-1)
    endif

!    nomalies

    if ( anom ) then
        if ( lwrite ) print '(a)','# Taking anomalies '
        do iens=nens1,nens2
            call anomal(fcst(1,yrbeg,iens),nperyear,nperyear,yrbeg,yrend,yr1,yr2)
        enddo
    endif

!       sum

    if ( lsum > 1 ) then
        if ( lwrite ) print '(a,i3)','# Summing series ',lsum
        do iens=nens1,nens2
            call sumit(fcst(1,yr1,iens),nperyear,nperyear,yr1,yr2,lsum2,'v')
        enddo
    endif

!       detrending

    if ( ldetrend ) then
        do iens=nens1,nens2
            if ( lwrite ) print *,'Detrending ens ',iens
            call detrend(fcst(1,yrbeg,iens),nperyear,nperyear,yrbeg,yrend,yr1,yr2,m1,m2,lsel)
        enddo
    endif
    if ( ndiff /= 0 ) then
        if ( lwrite ) print *,'Taking differences - series'
        do iens=nens1,nens2
            call diffit(fcst(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend,ndiff)
        enddo
    endif

!   copy ensemble members so that there is the same
!   number of valid ones at every time step

    ndup = 0
    if ( nens2 > nens1 .and. lmakeensfull ) then
        call makeensfull(ndup,nperyear,fcst,nperyear,yrbeg,yrend,nens1,nens2,validens,lwrite)
    endif

!   make table

    do i=1,len_trim(varfcst)
        if ( varfcst(i:i) == ' ' ) varfcst(i:i) = '_'
    enddo
    write(10,'(a)') varfcst
    do yr=yr1,yr2
        do month=m1,m2
            call getj1j2(j1,j2,month,nperyear, .false. )
            do mo1=j1,j2
                mo = mo1
                if ( mo > 12 ) mo = mo - 12
                do iens=nens1,nens2
                    if ( fcst(mo,yr,iens) == 3e33 ) fcst(mo,yr,iens) = -999.9
                enddo       ! iens
                !!!write(*,'(2i4,100g14.6)') yr,mo,(fcst(mo,yr,iens),iens=nens1,nens2)
                write(10,'(100g14.6)') (fcst(mo,yr,iens),iens=nens1,nens2)
            enddo           ! mo!
        enddo               ! month
    enddo                   ! yr
    close(10)

end program extractseries
