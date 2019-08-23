subroutine getenscutoff(cut,pcut,dat,npermax,nperyear,yrbeg &
    ,yrend,nensmax,nens1,nens2,yr1,yr2,j1,j2,lag)

!       compute the cut-off in absolute units given a cut-off as a
!       percentage (pcut), data and the range over which the data
!       should be considered
!       14-jan-2002
!       If the data is daily, a 5-day window is considered.

    implicit none
    integer,intent(in) :: npermax,nperyear,yrbeg,yrend,yr1,yr2,nensmax,nens1,nens2,j1,j2,lag
    real,intent(in) :: pcut,dat(npermax,yrbeg:yrend,0:nensmax)
    real,intent(out) :: cut

    integer :: yr,jj,m,j,i,ii,n,k,kdif,iens,nmax
    real,allocatable :: a(:)
    logical :: lwrite
    parameter (lwrite= .false. )

!   get linear array with data
    nmax = (j2-j1+1)*(yrend-yrbeg+1)*(nens2-nens1+1)
10  continue
    allocate(a(nmax))
    n = 0
    if ( j2 > nperyear ) then
        write(0,*) 'getenscutoff: error: cannot handle years wrapping: ',1,j2,nperyear
        write(*,*) 'getenscutoff: error: cannot handle years wrapping: ',j1,j2,nperyear
        call exit(-1)
    end if
    do iens=nens1,nens2
        do i=yr1,yr2
            do j=j1,j2
                if ( dat(j,i,iens) < 1e33 ) then
                    n = n+1
                    if ( n > nmax ) then
                        deallocate(a)
                        nmax = nmax*2
                        goto 10
                    endif
                    a(n) = dat(j,i,iens)
                endif
            enddo
        enddo
    enddo

    call getcut(cut,pcut,n,a)

    if ( lwrite ) then
        print *,'getcutoff: pcut          = ',pcut
        print *,'           yr1,yr2,j1,j2 = ',yr1,yr2,j1,j2
        print *,'           i,n           = ',i,n
        print *,'           cut           = ',cut
    endif
    deallocate(a)
end subroutine getenscutoff
