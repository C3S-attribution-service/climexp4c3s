        subroutine rmean(data,rdata,sdata,yrbeg,yrend,mean)
*       
*       stupid subroutine to compute the mean-month running mean 
*       and standard deviation of data
*
        implicit none
        integer yrbeg,yrend,mean
        real data(12,yrbeg:yrend),rdata(12,yrbeg:yrend),
     +        sdata(12,yrbeg:yrend)
        integer year,month,y1,m1,i
        integer yrmax
        real s
*
        do year=yrbeg,yrend
            do month=1,12
                rdata(month,year) = 0
                sdata(month,year) = 0
                do i=-mean/2,mean/2
                    m1 = month + i
                    call normon(m1,year,y1,12)
                    if ( y1.lt.yrbeg .or. y1.gt.yrend ) then
                        rdata(month,year) = 3e33
                        sdata(month,year) = 3e33
                    else
                        if ( rdata(month,year).lt.1e33 .and. data(m1,y1)
     +                        .lt.1e33 ) then
                            if ( mod(mean,2).eq.0 .and. abs(i).eq.mean/2
     +                            ) then
                                rdata(month,year) = rdata(month,year) +
     +                                data(m1,y1)/2
                                sdata(month,year) = sdata(month,year) +
     +                                data(m1,y1)**2/4
                                if ( i.eq.-mean/2 ) then
                                    s = data(m1,y1)
                                elseif ( i.eq.mean/2 ) then
                                    sdata(month,year) = sdata(month,year
     +                                    ) + data(m1,y1)*s/2
                                endif
***                                print *,'adding 1/2 data(',m1,y1,
***     +                                ') to rdata(',month,year,')'
                            else
                                rdata(month,year) = rdata(month,year) +
     +                                data(m1,y1)
                                sdata(month,year) = sdata(month,year) +
     +                                data(m1,y1)**2
***                                print *,'adding data(',m1,y1,
***     +                                ') to rdata(',month,year,')'
                            endif
                        else
                            rdata(month,year) = 3e33
                            sdata(month,year) = 3e33
                        endif
                    endif
                enddo
                if ( rdata(month,year).lt.1e33 ) then
                    rdata(month,year) = rdata(month,year)/mean
                    sdata(month,year) = sdata(month,year)/mean - 
     +                    rdata(month,year)**2
                    if ( sdata(month,year).ge.0 ) then
                        sdata(month,year) = sqrt(sdata(month,year))
                    else
                        write(0,*) 'rmean: error: variance < 0 '
     +                        ,sdata(month,year)
                        sdata(month,year) = 3e33
                    endif
                endif
***                print *,'rdata(',month,year,') = ',rdata(month,year)
            enddo
        enddo
*
        return
        end

