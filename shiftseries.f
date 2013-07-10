        subroutine shiftseries(data,npermax,nperyear,yrbeg,yrend,nshift)
*       
*       shift the time series over nshift units (normally months)
*
        implicit none
        integer npermax,nperyear,yrbeg,yrend,nshift
        real data(npermax,yrbeg:yrend)
        integer i,j,ii,k
*       
        if ( nshift.gt.0 ) then
            do i=yrend,yrbeg,-1
                do j=nperyear,1,-1
                    k = j + nshift
                    call normon(k,i,ii,nperyear)
                    if ( ii.le.yrend ) then
                        data(k,ii) = data(j,i)
                        data(j,i) = 3e33
                    endif
                enddo
            enddo
        elseif ( nshift.lt.0 ) then
            do i=yrbeg,yrend
                do j=1,nperyear
                    k = j + nshift
                    call normon(k,i,ii,nperyear)
                    if ( ii.ge.yrbeg ) then
                        data(k,ii) = data(j,i)
                        data(j,i) = 3e33
                    endif
                enddo
            enddo
        endif
        end
