        subroutine hipass(data,npermax,nperyear,yrbeg,yrend,ndiff,minfac
     +       )
*
*       Take anomalies wrt the mdiff years around the value
*
        implicit none
        integer npermax,nperyear,yrbeg,yrend,ndiff
        real data(npermax,yrbeg:yrend),minfac
        integer i,j,k,ii
        real s,absent,n
        parameter (absent=3e33)
*
        if ( ndiff.gt.0 ) then
            do i=yrbeg,yrend
                do j=1,nperyear
                    s = 0
                    n = 0
                    do k=-ndiff/2,ndiff/2
                        ii = i + k
                        if ( ii.lt.yrbeg .or. ii.gt.yrend ) cycle
                        if ( data(j,ii).lt.0.9*absent ) then
                            if ( abs(k).eq.ndiff/2 .and. 
     +                           mod(ndiff,2).eq.0 ) then
                                n = n + 0.5
                                s = s + data(j,ii)/2
                            else
                                n = n + 1
                                s = s + data(j,ii)
                            end if
                        end if
                    end do
                    ii = i - ndiff/2
                    if ( ii.ge.yrbeg ) then
                        if ( n.gt.minfac*ndiff ) then
                            data(j,ii) = data(j,i) - s/ndiff
                        else
                            data(j,ii) = absent
                        endif
                    end if
                end do
            end do
            call shiftseries(data,npermax,nperyear,yrbeg,yrend,
     +            nperyear*(ndiff/2))
        endif
*
        end

