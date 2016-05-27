subroutine seriesautocor(data,npermax,yrbeg,yrend,nperyear,yr1,yr2,ave,var,aa1,aa2)
!
!   compute day-to-day (month-to-month, aa1) and year-to-year (aa2) autocorrelations 
!   of a time series as a function of the season, given ave and var
!
    implicit none
    integer npermax,yrbeg,yrend,nperyear,yr1,yr2
    real data(npermax,yrbeg:yrend),ave(npermax),var(npermax)
    real s1,s2,aa1(nperyear),aa2(nperyear)
    integer yr,mo,nmax,i,ii1,j,m,m1,n1,n2,i1,i2
    
    do mo=1,nperyear
        n1 = 0
        n2 = 0
        s1 = 0
        s2 = 0
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
                m1 = j-1
                call normon(m1,yr,ii1,nperyear)
                if ( ii1 >= yr1 .and. i <= yr2 ) then
                    if ( data(m,i) < 1e33 .and. data(m1,ii1) < 1e33 ) then
                        n1 = n1 + 1
                        s1 = s1 + (data(m,i)-ave(m))*(data(m1,ii1)-ave(m1))
                    end if
                end if
                ii1 = i-1
                if ( ii1 >= yr1 .and. i <= yr2 ) then
                    if ( data(m,i) < 1e33 .and. data(m,ii1) < 1e33 ) then
                        n2 = n2 + 1
                        s2 = s2 + (data(m,i)-ave(m))*(data(m,ii1)-ave(m))
                    end if
                end if
            end do
        end do
        if ( n1 > 0 ) then
            s1 = s1/n1
            m1 = mo-1
            if ( m1 < 1 ) m1 = m1 + nperyear
            aa1(mo) = s1/sqrt(var(m1)*var(mo))
        else
            aa1(mo) = 3e33
        end if
        if ( n1 > 0 ) then
            s2 = s2/n2
            aa2(mo) = s2/var(mo)
        else    
            aa1(mO) = 3e33
        end if
    end do

end subroutine