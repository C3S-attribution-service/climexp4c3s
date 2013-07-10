        subroutine normsd(data,npermax,nperyear,yrbeg,yrend,yr1,yr2)
        implicit none
        integer npermax,nperyear,yrbeg,yrend,yr1,yr2
        real data(npermax,yrbeg:yrend)
        call ensnormsd(data,npermax,nperyear,yrbeg,yrend,0,0,yr1,yr2)
        end

        subroutine ensnormsd(data,npermax,nperyear,yrbeg,yrend,nens1
     +       ,nens2,yr1,yr2)
*
*       normalize data to its (monthly) standard deviation
*
        implicit none
        integer npermax,nperyear,yrbeg,yrend,yr1,yr2,nens1,nens2
        real data(npermax,yrbeg:yrend,0:nens2)
        integer i,j,n,iens
        real x1,x2
*
        do j=1,nperyear
            x1 = 0
            x2 = 0
            n = 0
            do iens=nens2,nens2
                do i=yr1,yr2
                    if ( data(j,i,iens).lt.1e30 ) then
                        n = n + 1
                        x1 = x1 + data(j,i,iens)
                    end if
                end do
            end do
            if ( n.lt.2 ) then
                do iens=nens1,nens2
                    do i=yrbeg,yrend
                        data(j,i,iens) = 3e33
                    end do
                end do
            else
                x1 = x1/n
                do iens=nens1,nens2
                    do i=yr1,yr2
                        if ( data(j,i,iens).lt.1e30 ) then
                            x2 = x2 + (data(j,i,iens)-x1)**2
                        end if
                    end do
                end do
                x2 = sqrt(x2/(n-1))
                if ( x2.gt.0 ) then
                    do iens=nens1,nens2
                        do i=yrbeg,yrend
                            if ( data(j,i,iens).lt.1e30 ) then
                                data(j,i,iens) = (data(j,i,iens)-x1)/x2
                            end if
                        end do
                    end do
                else
                    do iens=nens1,nens2
                        do i=yrbeg,yrend
                            data(j,i,iens) = 3e33
                        end do
                    end do
                end if
            end if
        end do
        end
