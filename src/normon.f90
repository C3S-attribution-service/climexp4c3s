subroutine normon(m,i,ii,nperyear)

!   normalize the month m to (1:nperyear), adjusting the year i to ii

    implicit none
    integer,intent(in) :: i,nperyear
    integer,intent(out) :: ii
    integer,intent(inout) :: m
    logical,parameter :: lwrite= .false.

    if ( lwrite ) then
        print *,'normon: input: m,i,nperyear = ',m,i,nperyear
    end if
    if ( m <= 0 ) then
        ii = i + m/nperyear - 1
        m = mod(m,nperyear) + nperyear
    elseif ( m > nperyear ) then
        ii = i + (m-1)/nperyear
        m = 1 + mod(m-1,nperyear)
    else
        ii = i
    endif
    if ( lwrite ) then
        print *,'normon: outnput: m,ii = ',m,ii
    end if
end subroutine normon
