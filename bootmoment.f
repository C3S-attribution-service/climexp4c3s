        subroutine bootmoment(moment,xx,yy,n,nboot,decor,result)
*
*       compute a moment of the array xx(1:n) and a bootstrap error
*       estimate by a sorted result(1:nboot) array.  The array yy(1:n) 
*       is working space. 
*
        implicit none
        integer moment,n,nboot
        real xx(n),yy(n),decor,result(0:nboot)
        integer i,ii,j,jj,iboot,nblock
        double precision ranf
        external ranf
        
        call getmoment(moment,xx,n,result(0))
        nblock = nint(decor+1)
        if ( nboot.gt.0 ) then
            do iboot=1,nboot
                do i=1,n
                    if ( nblock.eq.1 .or. mod(i,nblock).eq.1 ) then
                        j = int(1+n*ranf(i))
                        do ii=i,min(n,i+nblock-1)
                            jj = j + (ii-i)
                            if ( jj.gt.n ) jj = jj - n
                            yy(ii) = xx(jj)
                        enddo
                    endif
                enddo
                call getmoment(moment,yy,n,result(iboot))
            enddo
            call nrsort(nboot,result(1))
        endif
        end
