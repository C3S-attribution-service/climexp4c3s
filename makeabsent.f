        subroutine makeabsent(data,nperyear,yrbeg,yrend)
*
*       fill data with absent
*
        implicit none
        integer nperyear,yrbeg,yrend
        real data(nperyear,yrbeg:yrend)
        integer i,j
        real absent
        parameter (absent=3e33)
*
        do j=yrbeg,yrend
            do i=1,nperyear
                data(i,j) = absent
            enddo
        enddo
*
        end
