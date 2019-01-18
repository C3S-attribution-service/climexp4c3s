subroutine getsign(text,x,a,n,isgn,sign,lprint)

!   compute the significance of x given a null-hypothesis array a(n)

    implicit none
    include 'getopts.inc'
    integer,intent(in) :: n,isgn
    real,intent(in) :: x,a(n)
    real,intent(out) :: sign
    character,intent(in) :: text*(*)
    integer :: i
    logical :: lprint

!   watch out for nonsense input

    if ( abs(x) > 1e33 ) then
        if ( lprint ) print '(3a)','# ',text,': no valid data'
        return
    endif

!   otherwise, sort and find percentile

    call nrsort(n,a)
    call getsign1(sign,x,a,n,lwrite)
    if ( isgn < 0 ) sign = 1-sign
    if ( lprint ) then
        if ( lweb ) then
            print '(3a,f6.2,a,f6.4,a)','<tr><td> ',text,' </td><td>' &
                ,x,' </td><td> ',(1-sign),' </td></tr>'
        else
            print '(3a,f6.2,a,f6.4,a)','# ',text,' = ',x,', p = ',(1-sign)
        endif
    endif
end subroutine getsign
            
subroutine getsign1(sign,x,a,n,lwrite)

!   get probability of x given a sorted array of null-hypotheses a(n)

    implicit none
    integer,intent(in) :: n
    real,intent(in) :: x,a(n)
    real,intent(out) :: sign
    logical,intent(in) :: lwrite
    integer :: i

    if ( lwrite ) then
        print *,'x = ',x
        print *,'a = '
        do i=1,n
            print *,a(i)
        enddo
    endif
!   how many are larger/smaller than x?
!   should be done by bisection... or call to GSL...
    do i=1,n
        if ( a(i) > x ) go to 10
    enddo
    sign = 1 - 1./(n+1)
    goto 20
10  continue
    if ( i == 1 ) then
        sign = 1./(n+1)
        goto 20
    endif
    sign = (i*(x-a(i-1)) + (i-1)*(a(i)-x))/((n+1)*(a(i)-a(i-1)))
20  continue
end subroutine getsign1
