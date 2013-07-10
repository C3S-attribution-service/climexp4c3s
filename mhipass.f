        subroutine mhipass(data,npermax,nperyear,yrbeg,yrend,mdiff
     +       ,minfac)
*
*       Take anomalies wrt the mdiff months around the value
*
        implicit none
        integer npermax,nperyear,yrbeg,yrend,mdiff
        real data(npermax,yrbeg:yrend),minfac,n
        integer i,j,k,m,ii
        real s,absent
        parameter (absent=3e33)
*
        if ( mdiff.gt.0 ) then
            do i=yrbeg,yrend
                do j=1,nperyear
                    s = 0
                    n = 0
                    do k=-mdiff/2,mdiff/2
                        m = j+k
                        call normon(m,i,ii,nperyear)
                        if ( ii.lt.yrbeg .or. ii.gt.yrend ) cycle
                        if ( data(m,ii).lt.0.9*absent ) then
                            if ( abs(k).eq.mdiff/2 .and. 
     +                           mod(mdiff,2).eq.0 ) then
                                s = s + data(m,ii)/2
                                n = n + 0.5
                            else
                                s = s + data(m,ii)
                                n = n + 1
                            endif
                        endif
                    enddo
                    m = j - mdiff/2
                    call normon(m,i,ii,nperyear)
                    if ( ii.ge.yrbeg .and. ii.le.yrend ) then
                        if ( n.gt.minfac*mdiff ) then
                            data(m,ii) = data(j,i) - s/n
                        else
                            data(m,ii) = absent
                        end if
                    end if
                enddo
            enddo
            call shiftseries(data,npermax,nperyear,yrbeg,yrend,mdiff/2)
        endif
*
        end

