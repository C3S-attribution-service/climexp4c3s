program series

!   ancient program to convert a standard .dat file into a text file suitable for
!   gnuplot with columns for
!        1:   year
!        2:13 each month
!       14:17 each season (DJF,MAM,JJA,SON)
!       18:19 each half year (Oct-Mar, Apr-Sep)
!       20    full year Jan-Dec
!       21    full year Jul-Jun
!       22:41 same, but with 10-yr running mean (added 14-jul-2001)

    implicit none
    integer,parameter :: yrbeg=1600,yrend=2300,npermax=366,nensmax=750
    integer :: i,j,k,nperyear,itype,iens,mens1,mens,n,m,ilast
    real :: s(41,yrbeg:yrend),absent,ss,a1,sx,sy,w,t(41)
    real,allocatable :: data(:,:,:)
    parameter (absent=3e33)
    character line*256,var*40,units*20
    logical :: lprec,lwrite,lastvalid,lplot,linvalid
    integer :: iargc

!       init

    lwrite = .false. 
    if ( iargc() /= 1 .and. iargc() /= 3 ) then
        write(0,*) 'usage: series datafile [plot plotfile]'
        stop
    endif

    lplot = .false. 
    if ( iargc() == 3 ) then
        lplot = .true. 
        call getarg(2,line)
        if ( line(1:4) /= 'plot') then
            write(0,*) 'series: error: expecting "plot", found ',trim(line)
            call exit(-1)
        end if
        call getarg(3,line)
        open(1,file=trim(line),status='unknown') ! allow overwriting...
    end if
    call getarg(1,line)
    allocate(data(npermax,yrbeg:yrend,0:nensmax))
    call readensseries(line,data,npermax,yrbeg,yrend,nensmax &
        ,nperyear,mens1,mens,var,units, .false. ,lwrite)
    print '(5a)','# ',trim(var),' [',trim(units),']'
    if ( index(line,'/dd') > 0 .or. index(line,'dd') == 1 ) then
        write(0,*) 'Hmm. This looks like a (wind) direction to me.'
        write(0,*) 'Averaging over a unit circle.<p>'
        itype = 1
    else
        itype = 0
    endif
    a1 = atan(1.)/45

!   determine type of quantity - precipitation or not?

    do j=80,1,-1
        if ( line(j:j) == '/' ) goto 100
    enddo
    100 continue
    j = j+1
    if ( line(j:j) == 'p' ) then
        lprec = .true. 
        print '(2a)','# series: summed file ',line(1:index(line,' ') &
        -1)
    else
        lprec = .false. 
        print '(2a)','# series: averaged file ',line(1:index(line &
        ,' ')-1)
    endif
    print '(5a)','#year     Jan     Feb     Mar     Apr     May', &
        '     Jun     Jul     Aug     Sep     Oct     Nov', &
        '     Dec     DJF     MAM     JJA     SON Oct-Mar', &
        ' Apr-Sep Jan-Dec Jul-Jun; repeated with 10-yr running average'
    print '(a,101x,a)','#','(the year is the year of the last month)'

!   sum

    if ( nperyear /= 12 ) then
        write(0,*) 'series: unable to handle nperyear != 12 yet ',nperyear
        write(0,*) 'Sorry'
        call exit(-1)
    endif
    do i=yrbeg,yrend
    
!       months
    
        do j=1,12
            s(j+1,i) = 0
            do iens=mens1,mens
                s(j+1,i) = s(j+1,i) + data(j,i,iens)
            enddo
            s(j+1,i) = s(j+1,i)/(mens-mens1+1)
        enddo
    
!       seasons
    
        do j=1,4
            n = 0
            s(j+13,i) = 0
            do iens=mens1,mens
                if ( j == 1 ) then
                    if ( i == yrbeg .or. &
                         data(12,max(yrbeg,i-1),iens) > 1e33 .or. &
                         data(1,i,iens) > 1e33 .or. &
                         data(2,i,iens) > 1e33 ) then
                        s(j+13,i) = s(j+13,i)
                    elseif ( itype == 0 ) then
                        s(j+13,i) = s(j+13,i) + data(12,i-1,iens) + data(1,i,iens) + data(2,i,iens)
                        n = n + 1
                    elseif ( itype == 1 ) then
                        s(j+13,i) = s(j+13,i) + atan2( &
                            sin(a1*data(12,i-1,iens)) + &
                            sin(a1*data(1,i,iens)) + &
                            sin(a1*data(2,i,iens)), &
                            cos(a1*data(12,i-1,iens)) + &
                            cos(a1*data(1,i,iens)) + &
                            cos(a1*data(2,i,iens)))/a1
                        n = n + 1
                    endif
                else
                    if ( data(3*j-3,i,iens) > 1e33 .or. &
                         data(3*j-2,i,iens) > 1e33 .or. &
                         data(3*j-1,i,iens) > 1e33 ) then
                        s(j+13,i) = s(j+13,i)
                    elseif ( itype == 0 ) then
                        s(j+13,i) = s(j+13,i) + data(3*j-3,i,iens) + data(3*j-2,i,iens) + data(3*j-1,i,iens)
                        n = n + 1
                    elseif ( itype == 1 ) then
                        s(j+13,i) = s(j+13,i) + atan2( &
                            sin(a1*data(3*j-3,i,iens)) + &
                            sin(a1*data(3*j-2,i,iens)) + &
                            sin(a1*data(3*j-1,i,iens)), &
                            cos(a1*data(3*j-3,i,iens)) + &
                            cos(a1*data(3*j-2,i,iens)) + &
                            cos(a1*data(3*j-1,i,iens)))/a1
                        n = n + 1
                    endif
                endif
            enddo            ! iens
            if ( n == 0 ) then
                s(j+13,i) = 3e33
            elseif ( .not. lprec .and. itype == 0 ) then
                s(j+13,i) = s(j+13,i)/(3*n)
            else
                s(j+13,i) = s(j+13,i)/n
            endif
        enddo
    
!       half years
    
        n = 0
        s(18,i) = 0
        do iens=mens1,mens
            if ( i == yrbeg .or. &
                 data(10,max(yrbeg,i-1),iens) > 1e33 .or. &
                 data(11,max(yrbeg,i-1),iens) > 1e33 .or. &
                 data(12,max(yrbeg,i-1),iens) > 1e33 .or. &
                 data(1,i,iens) > 1e33 .or. &
                 data(2,i,iens) > 1e33 .or. &
                 data(3,i,iens) > 1e33 ) then
                s(18,i) = s(18,i)
            else
                if ( itype == 0 ) then
                    s(18,i) = s(18,i) + &
                        data(10,i-1,iens) + data(11,i-1,iens) + &
                        data(12,i-1,iens) + data(1,i,iens) + &
                        data(2,i,iens) + data(3,i,iens)
                    n = n + 1
                else
                    s(18,i) = s(18,i) + atan2( &
                        sin(a1*data(10,i-1,iens)) + &
                        sin(a1*data(11,i-1,iens)) + &
                        sin(a1*data(12,i-1,iens)) + &
                        sin(a1*data(1,i,iens)) + &
                        sin(a1*data(2,i,iens)) + &
                        sin(a1*data(3,i,iens)), &
                        cos(a1*data(10,i-1,iens)) + &
                        cos(a1*data(11,i-1,iens)) + &
                        cos(a1*data(12,i-1,iens)) + &
                        cos(a1*data(1,i,iens)) + &
                        cos(a1*data(2,i,iens)) + &
                        cos(a1*data(3,i,iens)))/a1
                    n = n + 1
                endif
            endif
        enddo
        if ( n > 0 ) then
            s(18,i) = s(18,i)/n
        else
            s(18,i) = 3e33
        endif
        n = 0
        s(19,i) = 0
        do iens=mens1,mens
            if ( data(4,i,iens) > 1e33 .or. &
                 data(5,i,iens) > 1e33 .or. &
                 data(6,i,iens) > 1e33 .or. &
                 data(7,i,iens) > 1e33 .or. &
                 data(8,i,iens) > 1e33 .or. &
                 data(9,i,iens) > 1e33 ) then
                s(19,i) = s(19,i)
            elseif ( itype == 0 ) then
                s(19,i) = s(19,i) + &
                    data(4,i,iens) + data(5,i,iens) + &
                    data(6,i,iens) + data(7,i,iens) + &
                    data(8,i,iens) + data(9,i,iens)
                n = n + 1
            else
                sx = 0
                sy = 0
                do j=4,9
                    sx = sx + cos(a1*data(j,i,iens))
                    sy = sy + sin(a1*data(j,i,iens))
                enddo
                s(19,i) = s(19,i) + atan2(sy,sx)/a1
                n = n + 1
            endif
        enddo
        if ( n == 0 ) then
            s(19,i) = 3e33
        else
            s(19,i) = s(19,i)/n
        endif
        if ( .not. lprec .and. itype == 0 ) then
            if ( s(18,i) < 1e30 ) s(18,i) = s(18,i)/6
            if ( s(19,i) < 1e30 ) s(19,i) = s(19,i)/6
        endif
    
!       years
    
        n = 0
        s(20,i) = 0
        do iens=mens1,mens
            if ( itype == 0 ) then
                sx = 0
                m = 0
                do j=1,12
                    if ( data(j,i,iens) < 1e30 ) then
                        sx = sx + data(j,i,iens)
                        m = m + 1
                    endif
                enddo
                if ( m == 12 ) then
                    s(20,i) = s(20,i) + sx
                    n = n + 1
                endif
            elseif ( itype == 1 ) then
                sx = 0
                sy = 0
                m = 0
                do j=1,12
                    if ( data(j,i,iens) < 1e33 ) then
                        sx = sx + cos(a1*data(j,i,iens))
                        sy = sy + sin(a1*data(j,i,iens))
                        m = m + 1
                    endif
                enddo
                if ( m == 12 ) then
                    s(20,i) = s(20,i) + atan2(sy,sx)/a1
                    n = n + 1
                endif
            endif
        enddo
        if ( n > 0 ) then
            s(20,i) = s(20,i)/n
        else
            s(20,i) = 3e33
        endif
        if ( i == yrbeg ) then
            s(21,i) = 3e33
        else
            n = 0
            s(21,i) = 0
            do iens=mens1,mens
                if ( itype == 0 ) then
                    sx = 0
                    m = 0
                    do j=1,6
                        if ( data(j,i,iens) < 1e33 ) then
                            sx = sx + data(j,i,iens)
                            m = m + 1
                        endif
                        if ( data(j+6,i-1,iens) < 1e33 ) then
                            sx = sx + data(j+6,i-1,iens)
                            m = m + 1
                        endif
                    enddo
                    if ( m == 12 ) then
                        s(21,i) = s(21,i) + sx
                        n = n + 1
                    endif
                elseif ( itype == 1 ) then
                    sx = 0
                    sy = 0
                    m = 0
                    do j=1,6
                        if ( data(j,i,iens) < 1e33 .and. &
                             data(j+6,i-1,iens) < 1e33 ) then
                            sx = sx + cos(a1*data(j,i,iens)) + &
                                cos(a1*data(j+6,i-1,iens))
                            sy = sy + sin(a1*data(j,i,iens)) + &
                                sin(a1*data(j+6,i-1,iens))
                            m = m + 2
                        endif
                    enddo
                    if ( m == 12 ) then
                        s(21,i) = s(21,i) + atan2(sy,sx)/a1
                        n = n + 1
                    endif
                endif
            enddo
            if ( n > 0 ) then
                s(21,i) = s(21,i)/n
            else
                s(21,i) = 3e33
            endif
        endif
        if ( .not. lprec .and. itype == 0 ) then
            if ( s(20,i) < 1e33 ) s(20,i) = s(20,i)/12
            if ( s(21,i) < 1e33 ) s(21,i) = s(21,i)/12
        endif
    
!       normalize
    
        if ( itype == 1 ) then
            do j=14,21
                if ( s(j,i) < 0 ) s(j,i) = s(j,i) + 360
            enddo
        endif
    enddo

!   compute 10-yr running average

    do i=yrbeg,yrend
        if ( i-5 < yrbeg .or. i+5 > yrend ) then
            do j=22,41
                s(j,i) = 3e33
            enddo
        else
            if ( itype == 0 ) then
                do j=2,21
                    ss = 0
                    w = 0
                    do k=-5,5
                        if ( s(j,i-k) < 1e33 ) then
                            if ( abs(k) == 5 ) then
                                w = w + 0.5
                                ss = ss + s(j,i-k)/2
                            else
                                w = w + 1
                                ss = ss + s(j,i-k)
                            endif
                        endif
                    enddo
                    if ( w > 7.9 ) then
                        s(j+20,i) = ss/w
                    else
                        s(j+20,i) = 3e33
                    endif
                enddo
            elseif ( itype == 1 ) then
                do j=2,21
                    sx = (cos(a1*s(j,i-5)) + cos(a1*s(j,i+5)))/2
                    sy = (sin(a1*s(j,i-5)) + sin(a1*s(j,i+5)))/2
                    if ( sx > 1e30 ) sx = 3e33
                    do k=-4,4
                        if ( s(j,i-k) > 1e33 ) then
                            sx = 3e33
                        else
                            sx = sx + cos(a1*s(j,i-k))
                            sy = sy + sin(a1*s(j,i-k))
                        endif
                    enddo
                    if ( sx > 1e33 ) then
                        s(j+20,i) = 3e33
                    else
                        s(j+20,i) = atan2(sy,sx)/a1
                    endif
                enddo
            
!               normalize
            
                do j=22,41
                    if ( s(j,i) < 0 ) s(j,i) = s(j,i) + 360
                enddo
            endif
        endif
    enddo

!   print data output

    lastvalid = .false. 
    do i=yrbeg,yrend
        s(1,i) = 3e33 ! it seems s(1,i) is not used and never set,
        ! we use it here as a flag whether there is at least one valid value
        do j=2,41
            if ( s(j,i) >= absent/3 ) then
                s(j,i) = -999.9
            else
                s(1,i) = 0
            endif
        enddo
        if ( s(1,i) >= absent/3 ) then
            if ( lastvalid ) print '(a)'
            lastvalid = .false. 
        else
            print '(i5,40g14.6)',i,(s(j,i),j=2,41)
            lastvalid = .true. 
        endif
    enddo

!   print plotting output, defining a few extra points for gnuplot

    if ( lplot ) then
        write(1,'(2a)') '# Warning: some fictitious data points ', &
            'were added to let gnuplot make nice plots!'
        do i=yrbeg,yrend
            if ( s(1,i) >= absent/3 ) then
                if ( lastvalid ) then
                    ! one line with the next year but the previous year's data
                    write(1,'(i5,40g14.6)') i,(s(j,i-1),j=2,41)
                    write(1,'(a)')
                end if
                lastvalid = .false. 
            else
                linvalid = .false. 
                do j=2,21
                    if ( s(j,i) == -999.9 ) then
                        linvalid = .true. 
                    end if
                end do
                if ( linvalid .and. i > yrbeg ) then ! at least one invalid value - trick
!                   one line with the undefs replaced by the previous line
                    t = s(:,i)
                    do j=2,21
                        if ( t(j) == -999.9 ) t(j) = s(j,i-1)
                    end do
                    write(1,'(i5,40g14.6)') i,(t(j),j=2,41)
!                   one line with the next year and the undef
!                   also if the new value is undef, replace with the
!                   previous one (i.e., repeat the first step already)
                    t = s(:,i+1)
                    do j=2,21
                        if ( s(j,i) == -999.9 ) then
                            t(j) = -999.9
                        else if ( t(j) == -999.9 ) then
                            t(j) = s(j,i)
                        end if
                    end do
                    write(1,'(i5,40g14.6)') i+1,(t(j),j=2,41)
                    ! and a blank line
                    write(1,'(a)')
                else        ! normal print
                    write(1,'(i5,40g14.6)') i,(s(j,i),j=2,41)
                end if
                lastvalid = .true. 
            endif
        end do
    end if
end program series
