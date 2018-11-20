subroutine mhipass(data,npermax,nperyear,yrbeg,yrend,mdiff,minfac)

!   Take anomalies wrt the mdiff months around the value

    implicit none
    integer,intent(in) :: npermax,nperyear,yrbeg,yrend,mdiff
    real,intent(in) :: minfac
    real,intent(inout) :: data(npermax,yrbeg:yrend)
    integer :: i,j,k,m,ii,iold,jold,n
    real :: s,sold
    real,parameter :: absent=3e33
    integer,external :: leap

    if ( mdiff == 1 ) then
        ! take simple differences
        sold = 3e33
        jold = -999
        do i=yrbeg,yrend
            do j=1,nperyear
                if ( nperyear == 366 .and. j == 31+29 .and. leap(i) == 1 ) cycle
                if ( sold < 1e33 .and. data(j,i) < 1e33 ) then
                    s = data(j,i) - sold
                else
                    s = 3e33
                end if
                if ( jold >= 1 ) data(jold,iold) = s
                jold = j
                iold = i
                sold = data(j,i)
            end do
        end do
    else if ( mdiff > 1 ) then
        do i=yrbeg,yrend
            do j=1,nperyear
                s = 0
                n = 0
                do k=-mdiff/2,mdiff/2
                    m = j+k
                    call normon(m,i,ii,nperyear)
                    if ( ii < yrbeg .or. ii > yrend ) cycle
                    if ( data(m,ii) < 0.9*absent ) then
                        if ( abs(k) == mdiff/2 .and. mod(mdiff,2) == 0 ) then
                            s = s + data(m,ii)/2
                            n = n + 0.5
                        else
                            s = s + data(m,ii)
                            n = n + 1
                        endif
                    endif
                enddo
                m = j - mdiff/2
                call normon(m,i,ii,nperyear)
                if ( ii >= yrbeg .and. ii <= yrend ) then
                    if ( n > minfac*mdiff ) then
                        data(m,ii) = data(j,i) - s/n
                    else
                        data(m,ii) = absent
                    end if
                end if
            enddo
        enddo
        call shiftseries(data,npermax,nperyear,yrbeg,yrend,mdiff/2)
    endif

end subroutine mhipass

