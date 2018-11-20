subroutine makeabsent(data,nperyear,yrbeg,yrend)

!   fill data with absent, should be inline in f90...

    implicit none
    integer :: nperyear,yrbeg,yrend
    real :: data(nperyear,yrbeg:yrend)
    real :: absent
    parameter (absent=3e33)

    data = absent

end subroutine makeabsent
