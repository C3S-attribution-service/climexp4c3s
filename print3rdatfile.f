        subroutine print3rdatfile(unit,data,npermax,nperyear,yrbeg,yrend
     +       )
        implicit none
        integer unit,npermax,nperyear,yrbeg,yrend
        real data(npermax,yrbeg:yrend,3)
        integer year,i,dy,mo,dpm(12,3),ical,ivar
        double precision val(3),offset
        logical todo
        data dpm /
     +       30,30,30,30,30,30,30,30,30,30,30,30,
     +       31,28,31,30,31,30,31,31,30,31,30,31,
     +       31,29,31,30,31,30,31,31,30,31,30,31/
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
        if ( ical.ne.3 ) then
            write(0,*) 'error: only for daily data, bot ',nperyear
            call abort
        end if
        call flush(unit)
        do year=yrbeg,yrend
            i = 0
            do mo=1,12
                do dy=1,dpm(mo,ical)
                    i = i + 1
                    val = -99.9
                    todo = .false.
                    do ivar=1,3
                        if ( data(i,year,ivar).lt.1e33 ) then
                            val(ivar) = data(i,year,ivar)
                            todo = .true.
                        end if
                    end do
                    if ( todo ) then
                        write(unit,'(i4,2i3,3f8.1)') year,mo,dy
     +                       ,(val(ivar),ivar=1,3)
                    end if
                end do
            end do
        end do
        return
        end
