subroutine rmean(data,rdata,sdata,yrbeg,yrend,mean)

!   stupid subroutine to compute the mean-month running mean and standard deviation of data

    implicit none
    integer,intent(in) :: yrbeg,yrend,mean
    real,intent(in) :: data(12,yrbeg:yrend)
    real,intent(out) :: rdata(12,yrbeg:yrend),sdata(12,yrbeg:yrend)
    integer :: year,month,y1,m1,i
    integer :: yrmax
    real :: s

    do year=yrbeg,yrend
        do month=1,12
            rdata(month,year) = 0
            sdata(month,year) = 0
            do i=-mean/2,mean/2
                m1 = month + i
                call normon(m1,year,y1,12)
                if ( y1 < yrbeg .or. y1 > yrend ) then
                    rdata(month,year) = 3e33
                    sdata(month,year) = 3e33
                else
                    if ( rdata(month,year) < 1e33 .and. data(m1,y1) < 1e33 ) then
                        if ( mod(mean,2) == 0 .and. abs(i) == mean/2 ) then
                            rdata(month,year) = rdata(month,year) + data(m1,y1)/2
                            sdata(month,year) = sdata(month,year) + data(m1,y1)**2/4
                            if ( i == -mean/2 ) then
                                s = data(m1,y1)
                            elseif ( i == mean/2 ) then
                                sdata(month,year) = sdata(month,year) + data(m1,y1)*s/2
                            endif
                        else
                            rdata(month,year) = rdata(month,year) + data(m1,y1)
                            sdata(month,year) = sdata(month,year) + data(m1,y1)**2
                        endif
                    else
                        rdata(month,year) = 3e33
                        sdata(month,year) = 3e33
                    endif
                endif
            enddo
            if ( rdata(month,year) < 1e33 ) then
                rdata(month,year) = rdata(month,year)/mean
                sdata(month,year) = sdata(month,year)/mean - rdata(month,year)**2
                if ( sdata(month,year) >= 0 ) then
                    sdata(month,year) = sqrt(sdata(month,year))
                else
                    write(0,*) 'rmean: error: variance < 0 ',sdata(month,year)
                    sdata(month,year) = 3e33
                endif
            endif
        !**                print *,'rdata(',month,year,') = ',rdata(month,year)
        enddo
    enddo

    return
end subroutine rmean

subroutine runmean(xx,yy,nperyear,k)

!   even more simple subroutine to take the k-day running mean of xx into yy

    implicit none
    integer,intent(in) :: nperyear,k
    real,intent(in) :: xx(nperyear)
    real,intent(out) :: yy(nperyear)
    integer :: i,ii,j,n
    real :: s
            
    do ii=1,nperyear
        s = 0
        n = 0
        do j=-k/2,k/2
            i = ii + j
            if ( i <= 0 ) i = i + nperyear
            if ( i > nperyear) i = i - nperyear
            if ( xx(i) < 1e33 ) then
                n = n + 1
                s = s + xx(i)
            end if
        end do
        if ( n >= k/2 ) then
            yy(ii) = s/n
        else
            yy(ii) = 3e33
        end if
    end do
end subroutine runmean
