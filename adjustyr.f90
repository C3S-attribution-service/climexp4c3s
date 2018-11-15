subroutine adjustyr(year1,year2,data,npermax,nperyear,yrbeg,yrend)

!   adjust the range yr1,yr2 to be no longer than the valid data

    implicit none
    include 'getopts.inc'
    integer :: npermax,nperyear,year1,year2,yrbeg,yrend
    real :: data(npermax,yrbeg:yrend)
    integer :: yr,month,i,j,jj,j1,j2,lag
    integer :: x1,x2,incr

    call getj1j2(j1,j2,m1,nperyear, .false. )
    if ( lwrite ) print *,'adjustyr: input  ',year1,year2,m1,m2
    x1 = year1
    x2 = year2
    if ( year1 < year2 ) then
        incr = 1
        x1 = max(yrbeg,x1)
        x2 = min(yrend,x2)
    else
        incr = -1
        x1 = min(yrend,x1)
        x2 = max(yrbeg,x2)
    endif
    do yr=x1,x2,incr
        do month=m1,m2
            if ( month == 0 ) then
                j1 = 1
                j2 = nperyear
            else
                call getj1j2(j1,j2,month,nperyear,.false.)
            endif
            do jj=j1,j2
                do lag=lag1,lag2
                    if ( fix2 ) then
                        j = jj+lag
                    else
                        j = jj
                    endif
                    call normon(j,yr,i,nperyear)
                    if ( i >= yrbeg .and. i <= yrend ) then
                        if ( data(j,i) < 1e33 ) go to 110
                    endif
                enddo
            enddo
        enddo
        if ( incr == 1 ) then
            yr1 = yr
        else
            yr2 = yr
        endif
    enddo
110 continue
    if ( lwrite ) print *,'adjustyr: output ',yr1,yr2
end subroutine adjustyr
