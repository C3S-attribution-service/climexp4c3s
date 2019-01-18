program multifit

!   read a correlate-dumpfile and do a simultaneous fit to data with n predictors
!   and write the resulting optimal compound index on stdout in dump format
!
!       Uses Numerical Recipes (I am lazy)

    implicit none
    integer,parameter :: indxmx=12,ndatmx=12000
    integer :: i,j,k,n,year(ndatmx),month(ndatmx)
    real :: data(ndatmx),indx(ndatmx,indxmx),chisq,sig(ndatmx), &
        x(ndatmx),u(ndatmx,indxmx),v(indxmx,indxmx),w(indxmx), &
        a(indxmx),da(indxmx,indxmx),nindx(ndatmx)
    common /cindx/ indx
    logical :: lwrite
    character string*255,indname(indxmx)*20
    external findx

!   process arguments

    lwrite = .false. 
    n = command_argument_count()
    if ( n < 2 ) then
        print *,'usage: multifit dumpfile index1 index2 ...'
        call exit(-1)
    end if
    call get_command_argument(1,string)
    i = index(string,' ')
    if ( lwrite ) print *,'opening old dumpfile ',string(1:i-1)
    open(1,file=string,status='old')
    string(i:i) = '1'
    if ( lwrite ) print *,'opening new dumpfile ',string(1:i)
    open(2,file=string)
    indname(1) = 'constant'
    do i=2,n
        call get_command_argument(i,indname(i))
    end do

!   read data

    i = 0
100 continue
    read (1,'(a)',end=200,err=900) string
    if ( string(1:1) == '#' ) then
        if ( lwrite ) print '(a)',trim(string)
        write(2,'(a)') trim(string)
        goto 100
    end if
    if ( string == ' ' ) goto 100
    if ( index(string,'999.9') /= 0 ) goto 100
    i = i + 1
    read(string,*,err=901)(indx(i,j),j=1,n-1),data(i),year(i) &
    ,month(i)
    goto 100
200 continue
    close(1)

!   call fit routine

    if ( i < n ) goto 902
    do j=1,i
        x(j) = j
    end do
    do j=1,i
        sig(j) = 1
    end do
    if ( lwrite ) print *,'calling svdfit'
    call svdfit(x,data,sig,i,a,n,u,v,w,ndatmx,indxmx,chisq,findx)
    print '(a,f12.2,a,i5,a,f12.4)','chisq/ndf = ',chisq,'/',i-n,' = ',chisq/(i-n)

!   get covariances of fit parameters

    if ( lwrite ) print *,'calling svdvar'
    call svdvar(v,n,indxmx,w,da,indxmx)

!   print output

    print '(a)','Best fit:'
    do j=1,n
        print '(a,g16.6,a,g16.6)',indname(j),a(j),' +/- ',sqrt(da(j,j))
    end do
    print '(a)','Correlation matrix:'
    do j=1,n
        print '(a,20f7.3)',indname(j),(da(k,j)/sqrt(da(j,j)*da(k,k)),k=1,j)
    end do

!   construct new index (without constant term)

    do j=1,i
        nindx(j) = 0
        do k=2,n
            nindx(j) = nindx(j) + a(k)*indx(j,k-1)
        end do
    end do

!   write new dumpfile

    do j=1,i
        write(2,'(2f12.4,2i5)') nindx(j),data(j),year(j),month(j)
    end do
    write(2,'(a,i2.2,a,f10.6)') ('# a',k,'=',a(k),k=1,n)
    close(2)

!   error messages

    goto 999
900 print *,'error reading dumpfile at ',string
    call exit(-1)
901 print *,'error reading data from string ',string
    call exit(-1)
902 print *,'not enough data points to fit ',i,n
    call exit(-1)
999 continue
end program multifit

subroutine findx(xi,f,n)
    implicit none
    integer :: n
    real :: xi,f(n)
    integer :: indxmx,ndatmx
    parameter (indxmx=12,ndatmx=12000)
    real :: indx(ndatmx,indxmx)
    common /cindx/ indx
    integer :: i,j

    i = nint(xi)
    if ( abs(xi-i) > 0.01 ) print *,'findx: wrong input! ',xi
    if ( n > indxmx ) print *,'findx: wrong input! ',n,indxmx
    f(1) = 1
    do j=2,n
        f(j) = indx(i,j-1)
    end do
end subroutine findx
