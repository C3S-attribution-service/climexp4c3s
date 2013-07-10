        subroutine printdatfile(unit,data,npermax,nperyear,yrbeg,yrend)
        implicit none
        integer unit,npermax,nperyear,yrbeg,yrend
        real data(npermax,yrbeg:yrend)
        integer year,i,dy,mo,dpm(12,3),ical
        data dpm /
     +       30,30,30,30,30,30,30,30,30,30,30,30,
     +       31,28,31,30,31,30,31,31,30,31,30,31,
     +       31,29,31,30,31,30,31,31,30,31,30,31/
        double precision val(360),offset
*
        if ( nperyear.lt.360 ) then
            ical = 0
        elseif ( nperyear.eq.360 ) then
            ical = 1
        elseif ( nperyear.eq.365 ) then
            ical = 2
        elseif ( nperyear.eq.366 ) then
            ical = 3
        else
            ical = 4
        endif
        call flush(unit)
        do year=yrbeg,yrend
            if ( ical.eq.0 ) then
                do i=1,nperyear
                    if ( data(i,year).lt.1e33 ) goto 200
                enddo
*               no valid points
                goto 210
*               there are valid points - print out
 200            continue
                do i=1,nperyear
                    if ( data(i,year).lt.1e33 ) then
                        val(i) = data(i,year)
                    else
                        val(i) = -999.9d0
                    endif
                enddo
                write(unit,'(i5,2000g15.7)') year,(val(i),i=1,nperyear)
 210            continue
            elseif ( ical.le.3 ) then
                i = 0
                do mo=1,12
                    do dy=1,dpm(mo,ical)
                        i = i + 1
                        if ( data(i,year).lt.1e33 ) then
                            write(unit,'(i4,2i3,g15.7)') year,mo,dy
     +                           ,data(i,year)
                        endif
                    enddo
                enddo
            else
                if ( nperyear.le.366 ) then
                    offset = 0.5d0 ! most likely averages
                else
                    offset = 1.0d0  ! most likely point measurements
                end if
                do i=1,nperyear
                    if ( data(i,year).lt.1e33 ) then
                        write(unit,'(2g15.7)') year+(i-offset)/nperyear
     +                       ,data(i,year)
                    endif
                enddo
            endif
        enddo
        return
        end
