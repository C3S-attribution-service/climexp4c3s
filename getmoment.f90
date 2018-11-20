subroutine getmoment(moment,xx,n,result)

!   compute the moment of xx(1:n) which does not contain missing data
!   I do not use the fancy ep-trick of numerical recipes.

    implicit none
    integer,intent(in) :: moment,n
    real,intent(in) :: xx(n)
    real,intent(out) :: result
    integer :: i
    real :: s1,s2,s3,s4
    logical,parameter :: lwrite=.false.

    if ( lwrite ) then
        print *,'getmoment: input: moment = ',moment
        print *,'                  n      = ',n
        if ( n < 20 ) then
            do i=1,n
                print *,i,xx(i)
            enddo
        endif
    endif

    if ( n < moment ) then
        result = 3e33
        return
    end if
    if ( moment >= 1 ) then
        s1 = 0
        do i=1,n
            s1 = s1 + xx(i)
        enddo
        s1 = s1/n
        if ( lwrite ) print *,'s1 = ',s1
    endif
    if ( moment >= 2 ) then
        s2 = 0
        do i=1,n
            s2 = s2 + (xx(i)-s1)**2
        enddo
        s2 = s2/(n-1)
        if ( s2 < 0 ) then
            write(0,*) 'getmoment: error: s2<0 ',s2
            s2 = 3e33
        else
            s2 = sqrt(s2)
        endif
        if ( lwrite ) print *,'s2 = ',s2
    endif
    if ( moment >= 3 ) then
        s3 = 0
        do i=1,n
            s3 = s3 + (xx(i)-s1)**3
        enddo
        if ( s2**3 == 0 ) then
!**         write(0,*) 'getmoment: error: s2 ~ 0',s2
            s3 = 3e33
        else
            s3 = s3/(n*s2**3)
        endif
        if ( lwrite ) print *,'s3 = ',s3
    endif
    if ( moment >= 4 ) then
        s4 = 0
        do i=1,n
            s4 = s4 + (xx(i)-s1)**4
        enddo
        if ( s2**4 == 0 ) then
            s4 = 3e33
        else
            s4 = s4/(n*s2**4) - 3
        endif
        if ( lwrite ) print *,'s4 = ',s4
    endif

    if ( moment == 1 ) then
        result = s1
    elseif ( moment == 2 ) then
        result = s2
    elseif ( moment == 3 ) then
        result = s3
    elseif ( moment == 4 ) then
        result = s4
    else
        result = 3e33
    endif

end subroutine getmoment
