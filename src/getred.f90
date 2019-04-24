subroutine getred(alpha,j1,j2,lag,k,nperyear,imens,indxmx,indx &
    ,data,npermax,yrbeg,yrend,nensmax,a,b)

!       compute the lag-1 autocorrelation that is needed to make red
!       noise

    implicit none
    include 'getopts.inc'
    integer :: j1,j2,lag,k,nperyear,indxmx,npermax,imens(0:indxmx),yrbeg,yrend,nensmax
    real :: data(npermax,yrbeg:yrend,0:nensmax),indx(npermax,yrbeg:yrend,0:nensmax,indxmx)
    real :: a,b,alpha
    integer :: iiens,iens,jens,yr,jj,j,i,m,ii,n
    integer :: yrfirst,yrlast,yrmin,dyr,jjfirst,jjlast,jjmin,djj
    real :: res,res0,s1a,s1b
    double precision :: s2a,s2b,s12
    logical :: llwrite

    llwrite = lwrite
!!!        llwrite = .true.
    if ( llwrite ) then
        print *,'getred: input'
        print *,'j1,j2 = ',j1,j2
        print *,'lag,k = ',lag,k
        print *,'nperyear,npermax = ',nperyear,npermax
        print *,'a,b   = ',a,b
        print *,'imens = ',imens
    endif
    if ( noisetype == 0 ) then
    ! (I cannot remember what this means...)
        if ( lwrite ) print *,'getred: noisetype = 0 => alpha = 0'
        alpha = 0
        return
    elseif ( noisetype /= 1 ) then
        write(0,*) 'getred: error: only know about noisetype=0,1, ','not ',noisetype
        write(*,*) 'getred: error: only know about noisetype=0,1, ','not ',noisetype
        call exit(-1)
    endif

!	first compute minimum time interval

    yrmin = 9999
    jjmin = 9999
    yrfirst = 9999
    jjfirst = 9999
    yrlast = 9999
    jjlast = 9999
    do iiens=nens1,nens2
        if ( imens(0) > 0 ) then
            iens = iiens
        else
            iens = 0
        endif
        if ( imens(k) > 0 ) then
            jens = iiens
        else
            jens = 0
        endif
        do yr=yr1-1,yr2
            do jj=j1,j2
                if ( jj == j1 .and. j1 /= j2 .and. ( j1 /= 1 .or. j2 /= nperyear ) ) then
                    res0 = 3e33
                endif
                if ( .not. fix2 ) then
                    j = jj-lag
                else
                    j = jj
                endif
                call normon(j,yr,i,nperyear)
                if ( i < yr1 .or. i > yr2 ) cycle
                m = j+lag
                call normon(m,i,ii,nperyear)
                if ( ii < yr1 .or. ii > yr2 ) cycle
                if (  data(j,i,iens) < 1e33 .and. indx(m,ii,jens,k) < 1e33 .and. &
                    ( lconting .or. ( &
                        (data(j,i,iens) <= maxdata) .eqv. &
                        (data(j,i,iens) >= mindata) .eqv. &
                        (maxdata >= mindata) ) .and. ( &
                        (indx(m,ii,jens,k) <= maxindx) .eqv. &
                        (indx(m,ii,jens,k) >= minindx) .eqv. &
                        (maxindx >= minindx) ) ) ) then
                    yrfirst = min(yr,yrfirst)
                    jjfirst = min(jj,jjfirst)
                    if ( yrlast < 9999 .and. yr /= yrlast ) then
                        yrmin = min(yrmin,yr-yrlast)
                    end if
                    yrlast = yr
                    			
                    if ( jjlast < 9999 .and. jj /= jjlast ) then
                        jjmin = min(jjmin,jj-jjlast)
                    end if
                    jjlast = jj
                end if
            end do
        end do
    end do
    if ( yrmin < 9999 ) then
        dyr = yrmin
    else
        dyr = 1
    end if
    if ( jjmin < 9999 ) then
        djj = jjmin
    else
        djj = 1
    end if
    if ( lwrite ) then
        print *,'getred: found yrfirst,jjfirst = ',yrfirst,jjfirst
        print *,'        found dyr,djj         = ',dyr,djj
    end if

!	real work

    n   = 0
    s1a = 0
    s1b = 0
    s2a = 0
    s2b = 0
    s12 = 0
    res0 = 3e33
    do iiens=nens1,nens2
        if ( imens(0) > 0 ) then
            iens = iiens
        else
            iens = 0
        endif
        if ( imens(k) > 0 ) then
            jens = iiens
        else
            jens = 0
        endif
        do yr=yrfirst,yr2,dyr
            do jj=jjfirst,j2,djj
                if ( jj == j1 .and. j1 /= j2 .and. &
                ( j1 /= 1 .or. j2 /= nperyear ) ) then
                    res0 = 3e33
                endif
                if ( .not. fix2 ) then
                    j = jj-lag
                else
                    j = jj
                endif
                call normon(j,yr,i,nperyear)
                if ( i < yr1 .or. i > yr2 ) go to 710
                m = j+lag
                call normon(m,i,ii,nperyear)
                if ( ii < yr1 .or. ii > yr2 ) go to 710
                if (  data(j,i,iens) < 1e33 .and. indx(m,ii,jens,k) < 1e33 .and. &
                    ( lconting .or. ( &
                        (data(j,i,iens) <= maxdata) .eqv. &
                        (data(j,i,iens) >= mindata) .eqv. &
                        (maxdata >= mindata) ) .and. ( &
                        (indx(m,ii,jens,k) <= maxindx) .eqv. &
                        (indx(m,ii,jens,k) >= minindx) .eqv. &
                        (maxindx >= minindx) ) ) ) then
                    res = data(j,i,iens) - a - b*indx(m,ii,jens,k)
                    if ( .false. .and. llwrite ) then
                        print '(3i5,f9.3,4i5,f9.3,f16.3)',j,i,iens &
                        ,data(j,i,iens),m,ii,jens,k,indx(m,ii,jens,k),res
                    endif
                else
                    res = 3e33
                endif
                if ( abs(res) < 1e33 .and. abs(res0) < 1e33 ) then
                    n = n+1
                    s1a = s1a + res0
                    s1b = s1b + res
                    s2a = s2a + res0**2
                    s2b = s2b + res**2
                    s12 = s12 + res0*res
                endif
                goto 720
            710 continue
                res = 3e33
            720 continue
                res0 = res
            enddo           ! month jj
        enddo               ! year yr
    enddo                   ! iens
!   do not assume that the mean is zero - I do not know why,
!   but I get only 0.1 times the standard deviation.  Bug in the fit
!   routine?
    if ( n > 6 ) then
        s1a = s1a/n
        s1b = s1b/n
        s12 = s12/n - s1a*s1b
        s2a = s2a/n - s1a**2
        s2b = s2b/n - s1b**2
        if ( llwrite ) then
            write(*,*) 'mean^2, var : ',s1a**2,s2a,n
            write(*,*) 'mean^2, var : ',s1b**2,s2b,n
        endif
        if ( s2a > 1e-10 .and. s2b > 1e-10 ) then
            alpha = s12/sqrt(s2a*s2b)
        else
            alpha = 0
        endif
    else
        alpha = 0
    endif
    if ( llwrite ) then
        print *,'getred: alpha = ',alpha
    endif
end subroutine getred

subroutine getred1(alpha,s2,j1,j2,nperyear,data,npermax,yrbeg,yrend,nensmax)

!   compute the lag-1 autocorrelation that is needed to make red noise

    implicit none
    include 'getopts.inc'
    integer :: j1,j2,nperyear,npermax,yrbeg,yrend,nensmax
    real :: data(npermax,yrbeg:yrend,0:nensmax)
    real :: alpha,s2
    integer :: iens,yr,jj,j,i,m,ii,n
    real :: res,res0,s2a,s2b,s12,s
    logical :: llwrite

    llwrite = lwrite
!**        llwrite = .true.
    if ( llwrite ) then
        print *,'getred1: input'
        print *,'j1,j2 = ',j1,j2
    endif
    n = 0
    s = 0
    do iens=nens1,nens2
        do yr=yr1-1,yr2
            do jj=j1,j2
                j = jj
                call normon(jj,yr,i,nperyear)
                if ( i < yr1 .or. i > yr2 ) cycle
                if ( data(j,i,iens) < 1e33 ) then
                    n = n + 1
                    s = s + data(j,i,iens)
                endif
            enddo
        enddo
    enddo
    if ( n < 6 ) then
        alpha = 3e33
        s2 = 3e33
        return
    endif
    s = s/n

    n   = 0
    s2a = 0
    s2b = 0
    s12 = 0
    res0 = 3e33
    do iens=nens1,nens2
        do yr=yr1-1,yr2
            do jj=j1,j2
                j = jj
                call normon(j,yr,i,nperyear)
                if ( i < yr1 .or. i > yr2 ) cycle
                if ( j == j1 .and. j1 /= j2 .and. ( j1 /= 1 .or. j2 /= nperyear ) ) then
                    res0 = 3e33
                endif
                if ( i < yrbeg .or. i > yrend ) then
                    res = 3e33
                elseif ( data(j,i,iens) < 1e33 ) then
                    res = data(j,i,iens) - s
                    if ( res < 1e33 .and. res0 < 1e33 ) then
                        n = n+1
                        s2a = s2a + res0**2
                        s2b = s2b + res**2
                        s12 = s12 + res0*res
                    endif
                else
                    res = 3e33
                endif
                res0 = res
            enddo           ! month jj
        enddo               ! year yr
    enddo                   ! iens
    if ( n > 10 ) then
        if ( s2a > 0 .and. s2b > 0 ) then
            alpha = s12/sqrt(s2a*s2b)
        else
            alpha = 0
        endif
    else
        alpha = 0
    endif
    s2 = sqrt(s2a*s2b)
    if ( llwrite ) then
        print *,'getred1: alpha,s2 = ',alpha,s2
    endif
end subroutine getred1
