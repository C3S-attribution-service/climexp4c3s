subroutine takerelanom(data,mean,npermax,yrbeg,yrend,nens1,nens2,nperyear,mon,lsum,lwrite)
!   take relative anomalies.
!   the variable lsum indocates the length of the seaon over
!   which the climatology should be averaged: annual, bi-annual, seasonal, monthly.
    implicit none
    integer,intent(in) :: npermax,yrbeg,yrend,nens1,nens2,nperyear,mon,lsum
    real,intent(inout) :: data(npermax,yrbeg:yrend,0:nens2),mean(npermax)
    logical,intent(in) :: lwrite
    integer :: i,j,iens,n
    real :: minfac,s

    minfac = 0.6 ! require at least this many months with data to give a value for a season...
    if ( mon == 0 .and. nperyear > 1 ) then
        ! special case for the IPCC WG1 AR5 Annex I Atlas - should never have made it's way into
        ! the general Climate Explorer programs, but it was a quick solution at a deadline :-(
        if ( lsum == 12 ) then
            s = sum(mean(1:nperyear))
            mean = s/nperyear
        else if ( lsum == nperyear/2 ) then
            if ( nperyear == 12 ) then
                s = sum(mean(4:9))
                mean(4:9) = s/6
                s = sum(mean(1:3)) + sum(mean(10:12))
                mean(1:3) = s/6
                mean(10:12) = s/6
            else
                goto 903
            end if
        else if ( lsum == 3 ) then
            if ( nperyear == 12 ) then
                do i=2,4
                    s = sum(mean(3*(i-1):3*i-1))
                    mean(3*(i-1):3*i-1) = s/3
                end do
                s = mean(12) + mean(1) + mean(2)
                mean(12) = s/3
                mean(1) = s/3
                mean(2) = s/3
            else if ( nperyear == 4 ) then
                if ( lwrite ) print *,'do nothing'
            else
                goto 903
            end if
        else if ( lsum == 12 ) then
            if ( nperyear == 12 ) then
                if ( lwrite ) print *,'do nothing'
            else
                goto 903
            end if
        else
            goto 903
        end if
    else if ( nperyear == 1 ) then
        if ( lwrite ) print *,'do nothing'
    else
    ! more general solution for the KNMI Atlas
        s = 0
        n = 0
        do i=1,lsum
            j = 1 + mod(mon+i-2,12)
            if ( mean(j) < 1e33 ) then
                n = n + 1
                s = s + mean(j)
            end if
        end do
        if ( n > minfac*lsum ) then
            mean = s/n
        else
            mean = 3e33
        end if
    end if
    do iens=nens1,nens2
        do i=yrbeg,yrend
            do j=1,nperyear
                if ( data(j,i,iens) < 1e30 .and. mean(j) /= 0 &
                 .and. mean(j) < 1e33 ) then
                    data(j,i,iens) = data(j,i,iens)/mean(j)
                else
                    data(j,i,iens) = 3e33
                end if
            end do
        end do
    end do
    goto 999
903 write(0,*) 'takerelanom: error: cannot take relative anomalies ' &
        ,'with nperyear = ',nperyear,' and lsum = ',lsum,' yet'
    call exit(-1)
999 continue
end subroutine takerelanom
