subroutine seriesmoment(data,npermax,yrbeg,yrend,nperyear,yr1,yr2,ave,adev,sdev,var,skew,curt)
!
!   compute some moments of a time series
!
    implicit none
    integer npermax,yrbeg,yrend,nperyear,yr1,yr2
    real data(npermax,yrbeg:yrend)
    real ave(nperyear),adev(nperyear),sdev(nperyear),var(nperyear),skew(nperyear),curt(nperyear)
    integer yr,mo,nmax,i,j,m,n,i1,i2
    real,allocatable :: xx(:)
    
    nmax = 31*(yr2-yr1+1)
    allocate(xx(nmax))
    do mo=1,nperyear
        n = 0
        if ( nperyear <= 36 ) then ! just that season / month / dekade
            i1 = mo
            i2 = mo
        else if ( nperyear >= 360 ) then ! take 31-day running period
            i1 = mo-15
            i2 = mo+15
        else ! pentads
            i1 = mo - 3
            i2 = mo + 3
        end if
        do yr=yr1,yr2
            do j=i1,i2
                m = j
                call normon(m,yr,i,nperyear)
                if ( i >= yr1 .and. i <= yr2 ) then
                    if ( data(m,i) < 1e33 ) then
                        n = n + 1
                        xx(n) = data(m,i)
                    end if
                end if
            end do
        end do
        call moment(xx,n,ave(mo),adev(mo),sdev(mo),var(mo),skew(mo),curt(mo))
    end do
    deallocate(xx)
end subroutine