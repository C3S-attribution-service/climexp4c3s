subroutine detrend(data,npermax,nperyear,yrbeg,yrend,yr1,yr2,m1,m2,lsel)

!   fit a line a + b*x to the data (one point per year) and subtract it

    implicit none
    integer,intent(in) :: npermax,nperyear,yrbeg,yrend,yr1,yr2,m1,m2,lsel
    real,intent(inout) :: data(npermax,yrbeg:yrend)

    integer,parameter :: nmax=60000
    integer :: i,n,yr,mo,i1,i2,month,mm2
    real :: x(nmax),y(nmax),sig(1),a,b,da,db,chi2,q
    integer,save :: init=0

    do month=m1,m2
        if ( month == 0 ) then
            i1 = 1
            i2 = nperyear
        else
            call getj1j2(i1,i2,month,nperyear, .false. )
        endif
        n = 0
        do yr=yr1,yr2
            do i=i1,i2
                mo = i
                if ( mo > nperyear ) mo = mo - nperyear
                if ( data(mo,yr) < 1e33 ) then
                    n = n + 1
                    if ( n > nmax ) then
                        print *,'detrend: compile with nmax larger'
                        call exit(-1)
                    endif
                    x(n) = yr + (mo-0.5)/nperyear - 2000
                    y(n) = data(mo,yr)
                endif
            enddo
        enddo
        if ( n > 1 ) then
            call fit(x,y,n,sig,0,a,b,da,db,chi2,q)
            if ( init < 10 ) then
                init = init + 1
                print '(a,f8.4,a,f8.4,a,i2)' &
                    ,'# detrend: subtracting a trend of ',b,' +/- ' &
                    ,db,' for month ',month
            !**                print *,'(a = ',a,' +/- ',da,')'
            endif
            do yr=yrbeg,yrend
                do i=i1,i2
                    mo = i
                    if ( mo > nperyear ) mo = mo - nperyear
                    if ( data(mo,yr) < 1e33 ) then
                        data(mo,yr) = data(mo,yr) - b*(yr + (mo-0.5)/nperyear - 2000)
                    endif
                enddo
            enddo
        endif
    enddo
    if ( .false. ) then
        do yr=yrbeg,yrend
            print '(i4,73g10.2)',yr,(data(mo,yr),mo=1,nperyear)
        enddo
    endif

end subroutine detrend
