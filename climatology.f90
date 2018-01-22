program climatology

!   compute a yearly climatology of a timeseries
!   output:
!   2000mmdd mean 2.5% 17% 50% 83% 97.5%

    implicit none
    include 'param.inc'
    integer,parameter :: ncut=5
    integer :: i,j,k,l,yr,mo,dy,nperyear,yr1,yr2,smooth,times, &
        nn(npermax),dpm(12),mens1,mens,nens1,nens2,iens,lsum, &
        firstyr,lastyr
    real :: data(npermax,yrbeg:yrend,0:nensmax),mean(npermax,0:ncut), &
        x(yrend-yrbeg+1),pcut(ncut)
    logical :: lexist,lstandardunits,lwrite
    character :: file*256,ensfile*255,string*10,oper*1,var*40,units*20
    character :: lvar*120,svar*120,history*50000,metadata(2,100)*1000
    integer :: iargc
    data pcut /.025,.17,.5,.83,.975/
    data dpm /31,29,31,30,31,30,31,31,30,31,30,31/

    lwrite = .FALSE. 
    if ( iargc() < 1 ) then
        print *,'usage: climatology file [begin yr] [end yr]'// &
            ' [ave|sum n ] [smooth n [times m]]'
        stop
    endif

    yr1 = yrbeg
    yr2 = yrend
    firstyr = 9999
    lastyr = -9999
    lsum = 1
    smooth = 1
    times = 1
    do i=2,iargc()-1,2
        call getarg(i,string)
        call getarg(i+1,file)
        if ( string(1:9) /= 'startstop' ) read(file,*) j
        if ( string(1:3) == 'beg' ) then
            yr1 = max(j,yrbeg)
        elseif ( string(1:3) == 'end' ) then
            yr2 = min(j,yrend)
        elseif ( string(1:3) == 'ave' ) then
            lsum = j
            oper = 'v'
        elseif ( string(1:3) == 'sum' ) then
            lsum = j
            oper = '+'
        elseif ( string(1:3) == 'smo' ) then
            smooth = j
        elseif ( string(1:3) == 'tim' ) then
            times = j
        elseif ( string(1:9) == 'startstop' ) then
            call getarg(i+1,file)
            inquire(file=trim(file),exist=lexist)
            if ( lexist ) then
                l = len_trim(file)
                do j=0,9999
                    write(file(l:),'(i4.4)') j
                    inquire(file=trim(file),exist=lexist)
                    if ( .NOT. lexist ) exit
                end do
            end if
            open(12,file=trim(file),status='new')
            if ( lwrite ) print '(2a)','# writing start and stop years on file ',trim(file)
        else
            print '(2a)','# climatology: disregarding unknown option ',trim(string)
        endif
    enddo

!   read data

    call getarg(1,file)
    lstandardunits = .FALSE. 
    call readensseriesmeta(file,data,npermax,yrbeg,yrend,nensmax &
        ,nperyear,mens1,mens,var,units,lvar,svar,history,metadata,lstandardunits,lwrite)
    nens1 = mens1
    nens2 = mens

!   average

    if ( lsum /= 1 ) then
        do iens=nens1,nens2
            call sumit(data(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend,lsum,'v')
        enddo
    endif

!   compute mean

    do j=1,nperyear
        mean(j,0) = 0
        nn(j) = 0
    enddo
    do iens=nens1,nens2
        do i=yr1,yr2
            do j=1,nperyear
                if ( data(j,i,iens) < 1e33 ) then
                    firstyr = min(firstyr,i)
                    lastyr = max(lastyr,i)
                    nn(j) = nn(j) + 1
                    mean(j,0) = mean(j,0) + data(j,i,iens)
                endif
            enddo
        enddo
    enddo
    do j=1,nperyear
        if ( nn(j) > 0 ) then
            mean(j,0) = mean(j,0)/nn(j)
        else
            mean(j,0) = 3e33
        endif
    enddo

!   compute percentiles

    do j=1,nperyear
        call keepalive1('Computing percentiles',j,nperyear)
        do i=1,ncut
            if ( (i >= 2 .OR. i <= ncut-1) .AND. nn(j) > 2 .OR. &
            (i == 1 .OR. i == ncut) .AND. nn(j) > 20 ) then
                call getenscutoff(mean(j,i),100*pcut(i),data,npermax &
                    ,nperyear,yrbeg,yrend,nensmax,nens1,nens2,yr1 &
                    ,yr2,j,j,0)
            else
                mean(j,i) = 3e33
            endif
        enddo
    enddo

!   smooth

    if ( smooth > 1 ) then
        write(0,'(a)') '# smoothing not yet ready'
        call abort
    endif

!   print results - twice for two year

    call savestartstop(firstyr,lastyr)

    print '(2a)','# Climatology of ',trim(file)
    call copyheadermeta(file,6,' ',history,metadata)
    print '(a)','#2000 + date        mean'// &
        '                2.5%'// &
        '                 17%'// &
        '                 50%'// &
        '                 83%'// &
        '               97.5%'
    do yr=2000,2001
        1000 format(i4.4,2i2.2,100g20.6)
        if ( nperyear == 2 ) then
            do j=1,2
                if ( mean(j,0) < 1e33 ) then
                    if ( j == 1 ) then
                        print 1000,yr-1,6*j+4,1,(mean(j,i),i=0,ncut)
                    else
                        print 1000,yr,6*j-8,1,(mean(j,i),i=0,ncut)
                    endif
                else
                    print '(a)'
                endif
            enddo
            if ( yr == 2001 .AND. mean(1,0) < 1e33 ) then
                print 1000,2001,10,1,(mean(1,i),i=0,ncut)
            endif
        else if ( nperyear == 4 ) then
            do j=1,4
                if ( mean(j,0) < 1e33 ) then
                    if ( j == 1 ) then
                        print 1000,yr-1,3*j+9,1,(mean(j,i),i=0,ncut)
                    else
                        print 1000,yr,3*j-3,1,(mean(j,i),i=0,ncut)
                    endif
                else
                    print '(a)'
                endif
            enddo
            if ( yr == 2001 .AND. mean(1,0) < 1e33 ) then
                print 1000,2001,12,1,(mean(1,i),i=0,ncut)
            endif
        else if ( nperyear == 12 ) then
            do j=1,12
                if ( mean(j,0) < 1e33 ) then
                    print 1000,yr,j,1,(mean(j,i),i=0,ncut)
                else
                    print '(a)'
                endif
            enddo
        elseif ( nperyear == 36 ) then
            do j=1,12
                do k=1,3
                    if ( mean(k+3*(j-1),0) < 1e33 ) then
                        print 1000,yr,j,1+10*(k-1),(mean(k+3*(j-1),i),i=0,ncut)
                    else
                        print '(a)'
                    endif
                enddo
            enddo
        elseif ( nperyear == 360 ) then
            do j=1,12
                do k=1,30
                    if ( mean(k+30*(j-1),0) < 1e33 ) then
                        print 1000,yr,j,k,(mean(k+30*(j-1),i),i=0,ncut)
                    else
                        print '(a)'
                    endif
                enddo
            enddo
        else
            k = nint(365.24/nperyear)
            dy = 0
            mo = 1
            do j=1,nperyear
                dy = dy + k
                if ( k > 1 .AND. mo == 2 .AND. dy >= 29 ) &
                dy = dy + 1
                if ( dy > dpm(mo) ) then
                    dy = dy - dpm(mo)
                    mo = mo + 1
                endif
                if ( mean(j,0) < 1e33 ) then
                    print 1000,yr,mo,dy,(mean(j,i),i=0,ncut)
                elseif ( mo /= 2 .OR. dy /= 29 ) then
                    print '(a)'
                endif
            enddo
        endif
    enddo
!   circumvent bug in gnuplot
    if ( nperyear > 4 .AND. mean(1,0) < 1e33 ) then
        print 1000,2002,1,1,(mean(1,i),i=0,ncut)
    endif

end program
