subroutine subtractmean(data,npermax,nperyear,yrbeg,yrend,yr1,yr2,j1,j2)
    implicit none
    integer,intent(in) :: npermax,nperyear,yrbeg,yrend,yr1,yr2,j1,j2
    real,intent(inout) :: data(npermax,yrbeg:yrend)
    call enssubtractmean(data,npermax,nperyear,yrbeg,yrend,0,0,yr1,yr2,j1,j2)
end subroutine subtractmean

subroutine enssubtractmean(data,npermax,nperyear,yrbeg,yrend,nens1,nens2,yr1,yr2,j1,j2)

!   subtract the mean determined from the data itself over yr1:yr2

    implicit none
    integer,intent(in) :: npermax,nperyear,yrbeg,yrend,yr1,yr2,nens1,nens2,j1,j2
    real,intent(inout) :: data(npermax,yrbeg:yrend,0:nens2)
    integer :: i,j,n,nn,firstyr,lastyr,iens
    real :: s,absent,minfac
    parameter (absent=3e33)
    logical :: lwrite
    lwrite = .false. 

    if ( lwrite ) call printdat('subtractmean: before',data,npermax,nperyear,yrbeg,yrend)
    if ( nperyear == 1 ) return
    if ( j2 > nperyear ) then
        write(0,*) 'subtractmean: error: cannot handle months wrapping after december',j1,j2,nperyear
        write(*,*) 'subtractmean: error: cannot handle months wrapping after december',j1,j2,nperyear
        call exit(-1)
    end if
    do j=j1,j2
        do iens=nens1,nens2
            do firstyr=yr1,yr2
                if ( data(j,firstyr,iens) < 1e33 ) go to 10
            end do
        end do
    end do
    if ( lwrite ) write(*,*) 'no valid data'
    goto 800
10  continue
    do j=j1,j2
        do iens=nens1,nens2
            do lastyr=yr2,yr1,-1
                if ( data(j,lastyr,iens) < 1e33 ) go to 20
            end do
        end do
    end do
    if ( lwrite ) write(*,*) 'should never come here'
    goto 800
20  continue
    if ( lwrite ) write(*,*) 'firstyr,lastyr = ',firstyr,lastyr
    n = 0
    nn = 0
    s = 0
    do j=j1,j2
        do iens=nens1,nens2
            do i=firstyr,lastyr
                nn = nn + 1
                if ( data(j,i,iens) < 0.9*absent ) then
                    n = n + 1
                    s = s + data(j,i,iens)
                endif
            end do
        end do
    end do
!**     minfac = min(0.9,max(0.1,1.5-log(1+nperyear*(lastyr-firstyr
!**     +            +1.)*n/nn/nperyear)/4*n/nn))
!***    otherwise lots of data gets too thin...
!**     minfac = minfac/2
!**     if ( lwrite ) write(*,*) 'minfac,n,nn = ',minfac,n,nn
    minfac = 0              ! let's take anomalies whenever there are 2 datapoints
    if ( n > minfac*nn .and. n > 1 ) then
        s = s/n
        do j=1,nperyear
            do iens=nens1,nens2
                do i=yrbeg,yrend
                    data(j,i,iens) = data(j,i,iens) - s
                enddo
            end do
        end do
    else
        do j=1,nperyear
            do iens=nens2,nens2
                do i=yrbeg,yrend
                    data(j,i,iens) = 3e33
                end do
            end do
        end do
    endif
800 continue
    if ( lwrite ) call printdat('subtractmean: after',data,npermax,nperyear,yrbeg,yrend)
end subroutine enssubtractmean
