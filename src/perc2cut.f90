subroutine perc2cut(lag,k,j1,j2,nperyear,imens1,imens,indxmx, &
    indx,data,npermax,yrbeg,yrend,nensmax)

!   if the boundaries are given as a percentage, convert to real numbers

    implicit none
    include 'getopts.inc'
    integer :: lag,k,j1,j2,nperyear,indxmx,npermax,yrbeg,yrend,nensmax
    integer :: imens1(0:indxmx),imens(0:indxmx)
    real :: indx(npermax,yrbeg:yrend,0:nensmax,indxmx),data(npermax,yrbeg:yrend,0:nensmax)
    integer :: m,n

    if ( fix2 ) then
        m = 0
        n = -lag
    else
        m = lag
        n = 0
    endif
    if ( pminindx > 0 .and. pminindx < 100 ) then
        call getenscutoff(minindx,pminindx,indx(1,yrbeg,0,k), &
            npermax,nperyear,yrbeg,yrend,nensmax,imens1(k),imens(k) &
            ,yr1,yr2,j1,j2,m)
    endif
    if ( pmaxindx > 0 .and. pmaxindx < 100 ) then
        call getenscutoff(maxindx,pmaxindx,indx(1,yrbeg,0,k), &
            npermax,nperyear,yrbeg,yrend,nensmax,imens1(k),imens(k) &
            ,yr1,yr2,j1,j2,m)
    endif
    if ( pmindata > 0 .and. pmindata < 100 ) then
        call getenscutoff(mindata,pmindata,data,npermax,nperyear, &
            yrbeg,yrend,nensmax,imens1(0),imens(0),yr1,yr2,j1,j2,n)
    endif
    if ( pmaxdata > 0 .and. pmaxdata < 100 ) then
        call getenscutoff(maxdata,pmaxdata,data,npermax,nperyear, &
            yrbeg,yrend,nensmax,imens1(0),imens(0),yr1,yr2,j1,j2,n)
    endif
end subroutine perc2cut
