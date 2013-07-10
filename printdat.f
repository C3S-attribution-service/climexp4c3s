        subroutine printdat(unit,data,yrbeg,yrend)
        implicit none
        integer unit,yrbeg,yrend
        real data(12,yrbeg:yrend)
        integer year,i
        real val(12)
*
        do year=yrbeg,yrend
            do i=1,12
                if ( data(i,year).ne.3e33 ) goto 200
            enddo
*           no valid points
            goto 210
*           there are valid points - print out
  200       continue
            do i=1,12
                if ( data(i,year).lt.1e33 ) then
                    val(i) = data(i,year)
                else
                    val(i) = -999.9
                endif
            enddo
            write(unit,'(i5,12f8.2)') year,val
  210       continue
        enddo
        return
        end
