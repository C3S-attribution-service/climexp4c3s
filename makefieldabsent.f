        subroutine makefieldabsent(field,nx,ny,yrbeg,yrend,yrfirst
     +        ,yrlast)
*
*       fill field with absent
*
        implicit none
        integer nx,ny,yrbeg,yrend,yrfirst,yrlast
        real field(nx,ny,12,yrbeg:yrend)
        integer i,j,jx,jy
        real absent
        parameter (absent=3e33)
*
        do i=yrfirst,yrlast
            do j=1,12
                do jy=1,ny
                    do jx=1,nx
                        field(jx,jy,j,i) = absent
                    enddo
                enddo
            enddo
        enddo
*
        end
