 subroutine autocov(data,npermax,yrbeg,yrend,nperyear,lagmax,cc,sig,yrstart,yrstop,lagunit)

!   conpute autocovariance of data

    implicit none
    include 'getopts.inc'
    integer,intent(in) :: npermax,yrbeg,yrend,nperyear,lagmax
    integer,intent(inout) :: yrstart,yrstop
    integer,intent(out) :: lagunit
    real,intent(in) :: data(npermax,yrbeg:yrend,0:nens2)
    real,intent(out) :: cc(0:lagmax),sig(0:lagmax)
    integer :: n,i,j,ii,jj,j1,j2,m,iens,lag,yr,mo
    real :: adata,s

    if ( m1 == 0 ) then
        j1 = 1
        j2 = nperyear
        lagunit = 1
    else
        j1 = m1
        if ( nperyear <= 12 ) then
!           consider year-on-year autocorrelations
            lagunit = nperyear
            j2 = m1
            lsum = lsum
        else
!           consider a seasonal selection
            lagunit = 1
            j2 = m1 + lsum - 1
            lsum = 1
        endif
        if ( nperyear /= 12 ) then
            call month2period(j1,nperyear,1)
            call month2period(j2,nperyear,0)
            if ( j2 < j1 ) j2 = j2 + nperyear
        !**             print *,'j1,j2 = ',j1,j2
        endif
    endif
    n = 0
    adata = 0
    do iens=nens1,nens2
        do ii=yr1,yr2
            do jj=j1,j2
                j = jj
                call normon(j,ii,i,nperyear)
                if ( data(j,i,iens) < 1e33 ) then
                    n = n + 1
                    adata = adata + data(j,i,iens)
                endif
            enddo
        enddo
    enddo
    if ( n == 0 ) then
        cc = 3e33
        sig = 0
        return
    endif
    adata = adata/n
    do lag=0,lagmax
        s = 0
        n = 0
        do iens=nens1,nens2
            do yr=yr1-1,yr2
                do mo=j1,j2
                    j = mo
                    call normon(j,yr,i,nperyear)
                    m = j + lagunit*lag
                    call normon(m,i,ii,nperyear)
                    if ( ii <= yr2 .and. ii >= yr1 .and. i <= yr2 .and. i >= yr1 ) then
                        if ( data(j,i,iens) < 1e33 .and. data(m,ii,iens) < 1e33 ) then
                            n = n+1
                            s = s + (data(j,i,iens)-adata)*(data(m,ii,iens)-adata)
                            yrstart = min(yrstart,i,ii)
                            yrstop  = max(yrstop,i,ii)
                        endif
                    endif
                enddo
            enddo
        enddo
        if ( n > 0 ) then
            cc(lag) = s/n
            if ( cc(lag) > +cc(0) ) cc(lag) = +cc(0)
            if ( cc(lag) < -cc(0) ) cc(lag) = -cc(0)
            sig(lag) = 2/sqrt(real(n)/max(real(1+(lsum-1)/lagunit),real(decor/lagunit)))
        else
            cc(lag) = 3e33
            sig(lag) = 0
        endif
    enddo
end subroutine autocov
