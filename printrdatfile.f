        subroutine printrdatfile(unit,data,npermax,nperyear,yrbeg,yrend
     +       ,yr1,yr2)
        implicit none
        integer unit,npermax,nperyear,yrbeg,yrend,yr1,yr2
        real data(npermax,yrbeg:yrend)
        integer year,i,dy,mo,dpm(12,3),ical
        data dpm /
     +       30,30,30,30,30,30,30,30,30,30,30,30,
     +       31,28,31,30,31,30,31,31,30,31,30,31,
     +       31,29,31,30,31,30,31,31,30,31,30,31/
        double precision val,offset
        integer,external ::leap
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
        end if
        call flush(unit)
        do year=yr1,yr2
            if ( ical.eq.0 ) then
                if ( nperyear.eq.12 ) then
                    do mo=1,12
                        if ( data(mo,year).lt.1e33 ) then
                            val = data(mo,year)
                        else
                            val = -999.9
                        end if
                        write(unit,'(i4,i3,a,g14.6)') year,mo,' 00',val
                    end do
                elseif ( nperyear.eq.1 ) then
                    if ( data(1,year).lt.1e33 ) then
                        val = data(i,year)
                    else
                        val =-999.9
                    end if
                    write(unit,'(i4,a,g14.6)') year,' 00 00',val
                else
                    write(0,*) 'printrdatfile: cannot handle nperyear ='
     +                   ,nperyear
                    call abort
                end if
            else if ( ical.le.3 ) then
                i = 0
                do mo=1,12
                    do dy=1,dpm(mo,ical)
                        i = i + 1
                        if ( data(i,year).lt.1e33 ) then
                            val = data(i,year)
                        else
                            val =-999.9
                        end if
                        if ( ical.eq.3 .and. mo.eq.2 .and. dy.eq.29 .and
     +                       .leap(year).eq.1 ) cycle
                        write(unit,'(i4,2i3,g14.6)') year,mo,dy,val
                    end do
                end do
            else
                write(0,*) 'printrdatfile: cannot handle nperyear ='
     +               ,nperyear
                call abort
            end if
        end do
        return
        end
