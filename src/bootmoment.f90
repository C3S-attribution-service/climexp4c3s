subroutine bootmoment(moment,xx,yy,n,nboot,decor,result)

!   compute a moment of the array xx(1:n) and a bootstrap error
!   estimate by a sorted result(1:nboot) array.  The array yy(1:n)
!   is working space.

    implicit none
    integer :: moment,n,nboot
    real :: xx(n),yy(n),decor,result(0:nboot),mean
    integer :: i,ii,j,jj,iboot,nblock
    real :: ranf

    if ( moment > 0 ) then
        call getmoment(moment,xx,n,result(0))
    else if ( moment == -2 ) then
        call getmoment(1,xx,n,mean)
        call getmoment(2,xx,n,result(0))
        result(0) = result(0)/mean
    else
        write(0,*) 'bootmoment: error: cannot handle moment = ',moment
        call exit(-1)
    end if
    nblock = nint(decor+1)
    if ( nboot > 0 ) then
        do iboot=1,nboot
            do i=1,n
                if ( nblock == 1 .or. mod(i,nblock) == 1 ) then
                    call random_number(ranf)
                    j = int(1+n*ranf)
                    do ii=i,min(n,i+nblock-1)
                        jj = j + (ii-i)
                        if ( jj > n ) jj = jj - n
                        yy(ii) = xx(jj)
                    enddo
                endif
            enddo
            if ( moment > 0 ) then
                call getmoment(moment,yy,n,result(iboot))
            else if ( moment == -2 ) then
                call getmoment(1,yy,n,mean)
                call getmoment(2,yy,n,result(iboot))
                result(iboot) = result(iboot)/mean
            end if
        enddo
        call nrsort(nboot,result(1))
    endif
end subroutine bootmoment
