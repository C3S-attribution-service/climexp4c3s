program getchance

!       read a file in 'dump.dat' format (see correlations.f) and
!       given the independent variable x computes chance intervals
!       of the dependent variable y

    implicit none
    integer :: nmax,npmax
    parameter (nmax=10000,npmax=100)
    integer :: i,j,ndata,npart,ii(nmax),ix,nx
    real :: x(nmax),y(nmax),xmin,xmax,xgiven,perc(npmax),frac,sig(1),a &
        ,b,siga,sigb,chi2,q,xx(npmax),yy,x0
    logical :: logscale,lwrite,sqrtscale
    character(1024) :: string

!       process arguments

    if ( command_argument_count() < 1 ) then
        print *,'usage: getchance [xgiven|:] [npart [file]]'
        stop
    endif
    lwrite = .false. 
    npart = 3
    call get_command_argument(1,string)
    if ( string == ':' ) then
        xmin = -3e33
        xmax = +3e33
    else
        read(string,*) xmin
        xmax = xmin
    endif
    if ( command_argument_count() >= 2 ) then
        call get_command_argument(2,string)
        read(string,*) npart
        if ( npart > npmax ) then
            print *,'error: increase npmax = ',npmax,' to ',npart
            call abort
        endif
    endif
    if ( command_argument_count() >= 3 ) then
        call get_command_argument(3,string)
        open(1,file=string)
        if ( lwrite ) print *,'using file &
        ',string(1:index(string,' ')-1)
    else
        open(1,file='dump.dat')
        if ( lwrite ) print *,'using file dump.dat'
    endif
    if ( xmin /= xmax ) then
        open(3,file=string(1:index(string,' ')-1)//'.trc')
    endif
!**	open(2,file='dump_cor.dat')

!       read data

    ndata = 0
    logscale = .false. 
    sqrtscale = .false. 
    100 continue
    read(1,'(a)',err=900,end=200) string
    if ( string(1:1) == '#' ) then
        if ( lwrite ) print '(a)',string
        if ( string == '# logarithmic plot' ) logscale = .true. 
        if ( string == '# sqrt plot' ) sqrtscale = .true. 
    !**            write(2,'(a)') string
        goto 100
    endif
    if ( string == ' ' ) goto 100
    ndata = ndata + 1
    read(string,*) x(ndata),y(ndata)
    goto 100
    200 continue
    if ( ndata < npart ) then
        print *,'error: less data (',ndata,') than partitions (' &
        ,npart,')'
        call exit(-1)
    endif

!       sort data and get n-ciles

    call ffsort(y,ii,ndata)
    print *,'original percentiles'
    do i=1,npart-1
        frac = (i*ndata/real(npart) + 0.5)
        j = int(frac)
        frac = frac - j
        perc(i) = y(ii(j))*(1-frac) + y(ii(j+1))*frac
        if ( logscale ) then
            print '(f6.2,f12.2)',100*i/real(npart),exp(perc(i))
        else
            print '(f6.2,f12.2)',100*i/real(npart),perc(i)
        endif
    enddo
!       get the median
    if ( mod(npart,2) == 0 ) then
        yy = perc(npart/2)
    else
        frac = (ndata/2 + 0.5)
        j = int(frac)
        frac = frac - j
        yy = y(ii(j))*(1-frac) + y(ii(j+1))*frac
    endif

!       fit data to straight line

    call fit(x,y,ndata,sig,0,b,a,sigb,siga,chi2,q)
    x0 = (yy-a)/b
    if ( lwrite ) then
        print *,'best fit: a = ',a,' +/- ',siga
        if ( logscale ) then
            print *,'          b = ',exp(b),' +/- ',sigb*exp(b)
            print *,'    median at ',x0,exp(yy)
        elseif ( sqrtscale ) then
            print *,'best fit: b = ',b**2,' +/- ',2*b*sigb
            print *,'    median at ',x0,yy**2
        else
            print *,'best fit: b = ',b,' +/- ',sigb
            print *,'    median at ',x0,yy
        endif
    endif

!	compute new distribution y -> y - a*(x-x0)

    do i=1,ndata
        y(i) = y(i) - a*(x(i)-x0)
    enddo
    call ffsort(y,ii,ndata)
!**	do i=1,ndata
!**	    write(2,'(i5,2f12.4)') i,x(ii(i)),y(ii(i))
!**	enddo

!       compute range of x-values

    if ( xmin < -1e33 ) then
        xmin = 1e30
        do i=1,ndata
            xmin = min(xmin,x(i))
        enddo
    !            if ( xmin.lt.0 ) then
    !                xmin = 1.2*xmin
    !            else
    !                xmin = 0.8*xmin
    !            endif
    endif
    if ( xmax > 1e33 ) then
        xmax = -1e30
        do i=1,ndata
            xmax = max(xmax,x(i))
        enddo
    !            if ( xmax.gt.0 ) then
    !                xmax = 1.2*xmax
    !            else
    !                xmax = 0.8*xmax
    !            endif
    endif
    if ( xmin == xmax ) then
        nx = 0
    else
        nx = 50
    endif
    do ix=0,nx
        if ( nx == 0 ) then
            xgiven = xmin
        else
            xgiven = xmin + ix*(xmax-xmin)/nx
        endif
    
    !           compute new percentiles given x=xgiven
    
        if ( nx == 0 ) print *,'percentiles given x = ',xgiven
        do i=1,npart-1
            yy = perc(i) - a*(xgiven-x0)
        !               (should be binary search)
            do j=1,ndata
                if ( y(ii(j)) > yy ) goto 300
            enddo
            300 continue
            if ( j == 1 ) then
                xx(i) = 0
            elseif ( j > ndata ) then
                xx(i) = 100
            else
                frac = (yy-y(ii(j-1)))/(y(ii(j))-y(ii(j-1)))
                xx(i) = 100*((j-1) + frac)/dble(ndata)
            endif
            if ( nx == 0 ) then
                if ( logscale ) then
                    print '(f6.2,2f12.2)',xx(i),exp(perc(i)),exp(yy)
                elseif ( sqrtscale ) then
                    print '(f6.2,2f12.2)',xx(i),perc(i)**2,yy**2
                else
                    print '(f6.2,2f12.2)',xx(i),perc(i),yy
                endif
            endif
        enddo
        if ( nx /= 0 ) then
            write(3,'(f16.4,99f8.2)') xgiven,(xx(i),i=1,npart-1)
        endif
    enddo
    goto 999

!       error messages

900 print *,'error reading from file'
    print '(a)',string
    call exit(-1)
999 continue
end program getchance