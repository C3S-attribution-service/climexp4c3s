subroutine diffit(data,npermax,nperyear,yrbeg,yrend,ndiff)

!   old entry point

    implicit none
    integer,intent(in) :: npermax,nperyear,yrbeg,yrend,ndiff
    real,intent(inout) :: data(npermax,yrbeg:yrend)
    real :: minfacsum
    minfacsum = 0.75
    call ndiffit(data,npermax,nperyear,yrbeg,yrend,ndiff,minfacsum)
end subroutine diffit

    subroutine ndiffit(data,npermax,nperyear,yrbeg,yrend,ndiff,minfacsum)

!   Take anomaly with respect to ndiff previous years, less when
!   undefined or near beginning of series (high-pass filter)
!   When ndiff < 0, take the sum of the n previous years
!   (poor man's low-pass filer).
!   7-dec-2001: only accept the point when at least 3/4 is valid data
!   29 november 2006: made this user-chosable.

    implicit none
    integer,intent(in) :: npermax,nperyear,yrbeg,yrend,ndiff
    real,intent(inout) :: data(npermax,yrbeg:yrend),minfacsum
    integer :: i,j,k,n,m
    real :: s
    real,parameter :: absent=3e33
    logical,parameter :: lwrite=.false.

    if ( lwrite) then
        print *,'ndiffit: ndiff,minfacsum = ',ndiff,minfacsum
        print *,'         yrbeg,yrend  = ',yrbeg,yrend
    end if
    if ( minfacsum < 0 ) minfacsum = 0.75
    if ( ndiff > 0 ) then
        do i=yrend,yrbeg,-1
            do j=1,nperyear
                if ( data(j,i) < 0.9*absent ) then
                    n = 0
                    s = 0
                    do k=max(i-abs(ndiff),yrbeg),i-1
                        if ( data(j,k) < 0.9*absent ) then
                            s = s + data(j,k)
                            n = n + 1
                        endif
                    enddo
                    if ( n >= minfacsum*abs(ndiff) ) then
                        data(j,i) = data(j,i) - s/n
                    else
                        data(j,i) = absent
                    endif
                else
                    data(j,i) = absent
                endif
            enddo
        enddo
    elseif ( ndiff < 0 ) then
        do i=yrend,yrbeg,-1
            m = 0
            do j=1,nperyear
                n = 0
                s = 0
                do k=max(i-abs(ndiff),yrbeg),i
                    if ( data(j,k) < 0.9*absent ) then
                        s = s + data(j,k)
                        n = n + 1
                    endif
                enddo
            !!!print *,'ndiff: ',i,j,n
                if ( n >= minfacsum*(1-ndiff) ) then
                    data(j,i) = s/n
                    m = m + 1
                else
                    data(j,i) = absent
                endif
            enddo
            if ( lwrite .and. m > 0 ) then
                print *,i,(data(j,i),j=1,nperyear)
            end if
        enddo
    endif

end subroutine ndiffit

subroutine dooverlap(fcst,npermax,yrbeg,yrend,nens1,nens2,nperyear,ndiff)
    implicit none
    integer,intent(in) :: npermax,yrbeg,yrend,nens1,nens2,nperyear,ndiff
    real,intent(inout) :: fcst(npermax,yrbeg:yrend,0:nens2)
    integer :: yr,mo,i,j,iens
    logical :: lvalid

!   make sure that the summed regions do not overlap

    do mo=1,nperyear
        yr = yrbeg
    101 continue
        do i=yr,yrend
!           search for a mo,yr with at least one valid point
            lvalid = .false. 
            do iens=nens1,nens2
                if ( fcst(mo,i,iens) < 1e30 ) then
                    lvalid = .true. 
                    exit
                end if
            end do
            if ( lvalid .and. i < yrend ) then
!               put the next -ndiff years to undefined
                do j=i+1,min(yrend,i-ndiff)
                    do iens=nens1,nens2
                        fcst(mo,j,iens) = 3e33
                    end do
                end do
!               and continue after that block
                yr = j
                goto 101
            end if
        end do
    end do
end subroutine dooverlap
