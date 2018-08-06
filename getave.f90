program getave

!   compute the average of variable 2 as a function of a bin in
!   variable one, given a file with (x,y) pairs

    implicit none
    integer,parameter :: nmax=10000
    integer :: i,j,n,m,ii(nmax),imin,imax
    real :: bin,x(nmax),y(nmax),xmin,xmax,s1,s2,sx
    character(80) :: string

    if ( command_argument_count() /= 1 ) then
        print *,'usage: getave binsize < infile'
        stop
    endif
    call get_command_argument(1,string)
    read(string,*) bin

!       read data

    n = 0
    100 continue
    read(*,'(a)',err=900,end=200) string
    if ( string(1:1) == '#' .or. string == ' ' ) goto 100
    n = n + 1
    if ( n > nmax ) goto 901
    read(string,*) x(n),y(n)
    goto 100
    200 continue

!       process data

    call ffsort(x,ii,n)
    imin = x(ii(1))/bin
    if ( imin < 0 ) imin = imin - 1
    imax = x(ii(n))/bin
    if ( imax < 0 ) imax = imax - 1
    j = 0
    print '(a)','#   x      ymean      +/-    no points'
    do i=imin,imax
        m = 0
        s1 = 0
        s2 = 0
        sx = 0
        xmin = i*bin
        xmax = (i+1)*bin
    
    !           collect bin
    
        300 continue
        j = j + 1
        if ( j > n ) goto 400
        if ( x(ii(j)) > xmax ) goto 400
        m = m + 1
        s1 = s1 + y(ii(j))
        s2 = s2 + y(ii(j))**2
        sx = sx + x(ii(j))
        goto 300
        400 continue
        if ( m > 1 ) then
            print '(f7.3,2f10.3,i6)',sx/m,s1/m,sqrt(s2-s1**2/m)/(m-1 &
            ),m
        else
            print '(a)','# less than two points in bin'
        endif
    enddo

!       error messages

    stop
900 print *,'error reading input '
    print *,string
    stop
901 print *,'too many points, increase nmax = ',nmax
end program getave
