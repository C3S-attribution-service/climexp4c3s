subroutine fieldanomal(field,nxf,nyf,nperyear,yr1,yr2,nx,ny)

!   take anomalies wrt the seasonal cycle

    implicit none
    integer,intent(in) :: nxf,nyf,nperyear,yr1,yr2,nx,ny
    real,intent(inout) :: field(nxf,nyf,nperyear,yr1:yr2)
    integer :: mo,yr,i,j,n
    real :: s

    do j=1,ny
        do i=1,nx
            do mo=1,nperyear
                s = 0
                n = 0
                do yr=yr1,yr2
                    if ( field(i,j,mo,yr) < 1e30 ) then
                        n = n + 1
                        s = s + field(i,j,mo,yr)
                    end if
                end do
                if ( n > 4 ) then
                    s = s/n
                else
                    s = 3e33
                end if
                do yr=yr1,yr2
                    if ( field(i,j,mo,yr) < 1e30 .and. s < 1e30 ) then
                        field(i,j,mo,yr) = field(i,j,mo,yr) - s
                    else
                        field(i,j,mo,yr) = 3e33
                    end if
                end do
            end do
        end do
    end do
end subroutine fieldanomal

