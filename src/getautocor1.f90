subroutine getautocor1(x,y,n,a1,a2,lwrite)

!   get the lag-1 and lag-2 autocorrelation of the residuals of a straight-line fit through (x,y)

    implicit none
    integer,intent(in) :: n
    real,intent(in) :: x(n),y(n)
    real,intent(out) :: a1,a2
    logical,intent(in) :: lwrite
    integer :: i
    real :: sig(1),mean1,var1,a,b,siga,sigb,chi2,q
    real,allocatable :: res(:)

    a1 = 3e33
    a2 = 3e33
    if ( n <= 2 ) then
        return
    end if
!   first compute residuals to the fit (which means we'll be fitting yet another time -
!   computers are fast these days)
    allocate(res(n))
    call fit(x,y,n,sig,0,a,b,siga,sigb,chi2,q)
    if ( a > 1e33 .or. b > 1e33 ) return
    do i=1,n
        res(i) = y(i) - (a+b*x(i))
    end do
!   the mean of residuals should be close to zero...
    mean1 = 0
    do i=1,n
        mean1 = mean1 + res(i)
    enddo
    mean1 = mean1/n
!   variance
    var1 = 0
    do i=1,n
        var1 = var1 + (res(i)-mean1)**2
    enddo
    var1 = var1/(n-1)
    if ( var1 == 0 ) then
        a1 = 1
        a2 = 1
    else
!       and autocorrelations the old-fashioned way.
        a1 = 0
        do i=1,n-1
            a1 = a1 + (res(i)-mean1)*(res(i+1)-mean1)
        enddo
        a1 = a1/(n-1)
        a1 = a1/var1
        a2 = 0
        do i=1,n-2
            a2 = a2 + (res(i)-mean1)*(res(i+2)-mean1)
        enddo
        a2 = a2/(n-2)
        a2 = a2/var1
    end if
    deallocate(res)
    if ( lwrite ) then
        print *,'getautocor1: found'
        print *,'   mean = ',mean1
        print *,'   s.d. = ',sqrt(var1)
        print *,'   a1,a2= ',a1,a2
    end if
end subroutine getautocor1
