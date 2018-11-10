subroutine mdiffit(data,npermax,nperyear,yrbeg,yrend,mdiff)

!   Take anomalies wrt the mdiff previous months

    implicit none
    integer,intent(in) :: npermax,nperyear,yrbeg,yrend,mdiff
    real,intent(inout) :: data(npermax,yrbeg:yrend)
    integer :: i,j,m,n,ii
    real :: s
    real,parameter :: absent=3e33

    if ( mdiff > 0 ) then
        do i=yrend,yrbeg,-1
            do j=nperyear,1,-1
                if ( data(j,i) < 0.9*absent ) then
                    s = 0
                    do n=1,mdiff
                        m = j-n
                        call normon(m,i,ii,nperyear)
                        if ( ii < yrbeg ) then
                            s = absent
                        else if ( s < 0.9*absent .and. &
                            data(m,ii) < 0.9*absent ) then
                            s = s + data(m,ii)
                        else
                            s = absent
                        end if
                    end do
                    if ( s < 0.9*absent ) then
                        data(j,i) = data(j,i) - s/mdiff
                    else
                        data(j,i) = absent
                    end if
                end if
            end do
        end do
    end if

end subroutine mdiffit
