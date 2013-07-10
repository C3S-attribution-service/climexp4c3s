        subroutine mdiffit(data,npermax,nperyear,yrbeg,yrend,mdiff)
*
*       Take anomalies wrt the mdiff previous months
*
        implicit none
        integer npermax,nperyear,yrbeg,yrend,mdiff
        real data(npermax,yrbeg:yrend)
        integer i,j,m,n,ii
        real s,absent
        parameter (absent=3e33)
*
        if ( mdiff.gt.0 ) then
            do i=yrend,yrbeg,-1
                do j=nperyear,1,-1
                    if ( data(j,i).lt.0.9*absent ) then
                        s = 0
                        do n=1,mdiff
                            m = j-n
                            call normon(m,i,ii,nperyear)
                            if ( ii.lt.yrbeg ) then
                                s = absent
                            elseif ( s.lt.0.9*absent .and. 
     +                                data(m,ii).lt.0.9*absent ) then
                                s = s + data(m,ii)
                            else
                                s = absent
                            endif
                        enddo
                        if ( s.lt.0.9*absent ) then
                            data(j,i) = data(j,i) - s/mdiff
                        else
                            data(j,i) = absent
                        endif
                    endif
                enddo
            enddo
        endif
*
        end
