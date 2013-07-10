        program dat2sound
*
*       convert a standard .,dat file into a raw data file 
*       that sox can read
*
        implicit none
        integer firstyear,lastyear,rate
        parameter(firstyear=1,lastyear=2020,rate=440*12)
        integer yr,mo,offset,repeat,last
        real xmax,xmin,data(12,firstyear:lastyear)
        character file*256
        integer iargc
*
        if ( iargc().ne.1 ) then
            print *,'usage: dat2sound infile.dat'
            stop
        endif
        call getarg(1,file)
        call makeabsent(data,12,firstyear,lastyear)
        call readdat(data,firstyear,lastyear,file)
        print '(a,i8)','; Sample Rate ',rate
        xmax = -3e33
        xmin = 3e33
        offset = 0
        do yr=firstyear,lastyear
            do mo=1,12
                if ( data(mo,yr).lt.1e33 ) then
                    if ( offset.eq.0 ) offset = 12*yr+mo
                    last = 12*yr+mo
                    xmax = max(xmax,data(mo,yr))
                    xmin = min(xmin,data(mo,yr))
                endif
            enddo
        enddo
        do repeat=1,10
            do yr=firstyear,lastyear
                do mo=1,12
                    if ( data(mo,yr).lt.1e33 ) then
                        print *,(12*yr+mo - offset + (repeat-1)*
     +                        (last-offset+1))/real(rate),
     +                        2*(data(mo,yr)-xmin)/(xmax-xmin) - 1
                    endif
                enddo
            enddo
        enddo
        end
