subroutine fit2(x,y,n,a,b,c,da,db,dc,r,prob,df,ncrossvalidate,lwrite)

!       compute a quadratic fit
!       y = a + b*x + c*x**2
!       with 1-sigma error bounds, chi2 and sqrt(explained variance)

    implicit none
    integer,parameter :: nmax=366*500
    integer,intent(in) :: n,ncrossvalidate
    real,intent(in) :: x(n),y(n)
    real,intent(out) :: a,b,c,da,db,dc,r,prob,df
    logical,intent(in) :: lwrite
    integer :: i
    real :: aa(3),daa(3,3),u(nmax,3),v(3,3),w(3),chi2
    real :: z,adata,sxx,aindx,syy,sxy
    real :: xi(nmax),sig(nmax)
    integer :: nn
    real :: xx(nmax)
    common /c_findx2/ nn,xx
    external findx2

    if ( lwrite ) then
        print *,'fit2: input: df = ',df
        print *,'x,y = '
        do i=1,n,n-1
            print *,i,x(i),y(i)
        enddo
    endif

!   copy to common
    nn = n
    xx(1:n) = x(1:n)
!   prepare arguments
    do i=1,n
        xi(i) = i
    enddo
    do i=1,n
        sig(i) = 1
    enddo
!   call Numerical Recipes routine
    call svdfit(xi,y,sig,n,aa,3,u,v,w,nmax,3,chi2,findx2)
    call svdvar(v,3,3,w,daa,3)
    a = aa(1)
    b = aa(2)
    c = aa(3)
    da = sqrt(daa(1,1))
    db = sqrt(daa(2,2))
    dc = sqrt(daa(3,3))
    if ( lwrite ) then
        print *,'a,b,c    = ',a,b,c
        print *,'da,db,dc = ',da,db,dc
    endif

!   compute explained variance

    do i=1,n
        xx(i) = a + b*x(i) + c*x(i)**2
    enddo
    call pearsnxx(y,xx,n,r,prob,z,adata,sxx,aindx,syy,sxy,df,ncrossvalidate)

end subroutine fit2

subroutine findx2(xi,f,n)

!   used by quadratic fitting routine fit2

    implicit none
    integer,parameter :: nmax=366*500
    integer,intent(in) :: n
    real,intent(in) :: xi
    real,intent(out) :: f(n)
    integer :: nn
    real :: xx(nmax)
    common /c_findx2/ nn,xx
    integer :: i,j

    if ( n /= 3 ) go to 901
    i = nint(xi)
    if ( abs(xi-i) > 0.01 .or. i < 1 .or. i > nn ) go to 902
    f(1) = 1
    f(2) = xx(i)
    f(3) = xx(i)**2
    return
901 print *,'findx2: should be called with n=3, not ',n
    call exit(-1)
902 print *,'findx2: wrong input! ',xi
    call exit(-1)
end subroutine findx2
