subroutine fitgam(xx,ntot,mean,sd,a,b)

!   a fit a gamma distribution to the data

    implicit none
    integer :: ntot
    real :: xx(ntot),mean,sd,a,b
    integer :: i,nx,iter,init
    real :: tol,p(3,2),q(2),xmin,fret,y(3)
    integer,parameter :: nmax=10000000
    integer :: ncur
    real :: data(nmax),xrestrain
    logical :: llwrite,llchangesign
    common /fitdata1/ data
    common /fitdata2/ xrestrain,ncur,llwrite,llchangesign
    real,external :: llgamma,dllgamma

    data init /0/

!   check input

    xmin = 3e33
    do i=1,ntot
        if ( xx(i) < 0 ) then
            write(0,*) 'histogram: error: cannot fit Gamma to negative data: ',xx(i)
            call exit(-1)
        end if
        if ( xx(i) /= 0 ) then
            xmin = min(xmin,xx(i))
        end if
    end do
    xmin = xmin/5
    do i=1,ntot
        if ( xx(i) == 0 ) then
            if ( init == 0 ) then
                init = 1
                print '(a,f16.6)','# changed zero to ',xmin
            end if
            xx(i) = xmin
        end if
    end do

!   copy to common for routines llgamma, dllgamma

    ncur = ntot
    do i=1,ncur
        data(i) = xx(i)
    end do

!       fit, using Numerical Recipes routines

!**     This frequently crashes
!**        q(1) = mean**2/sd**2
!**        q(2) = sd**2/mean
!**        tol = 1e-4
!**        call dfpmin(q,2,tol,iter,fret,llgamma,dllgamma)
!**        a = q(1)
!**        b = q(2)
!**     so try amoeba - slow but sure
    p(1,1) = mean**2/sd**2-0.05
    p(1,2) = sd**2/mean    *0.9
    p(2,1) = p(1,1) +0.1
    p(2,2) = p(1,2)
    p(3,1) = p(1,1)
    p(3,2) = p(1,2) *1.2
    do i=1,3
        q(1) = p(i,1)
        q(2) = p(i,2)
        y(i) = llgamma(q)
    end do
    tol = 1e-4
    call amoeba(p,y,3,2,2,tol,llgamma,iter)
!   maybe add restart later
    a = p(1,1)
    b = p(1,2)

!   output

    print '(a,i5,a)','# Fitted to gamma distribution in ',iter,' iterations'
    print '(a)','# p(x) = (x/b)^(a-1)*exp(-x/b)/(b*Gamma(a)) with'
    print '(a,f16.3)','# a = ',a
    print '(a,f16.3)','# b = ',b
end subroutine fitgam

real function llgamma(p)

!   computes the log-likelihood function for a Gamma distribution
!   with parameters alpha,beta=p(1),p(2) and data in common.

    implicit none
    include 'getopts.inc'
    real :: p(2)
    integer :: i
    integer,parameter :: nmax=10000000
    integer :: ncur
    real :: data(nmax),xrestrain
    logical :: llwrite,llchangesign
    common /fitdata1/ data
    common /fitdata2/ xrestrain,ncur,llwrite,llchangesign

    real,external :: gammln,gammp,gammq

    llgamma = 0
    do i=1,ncur
        llgamma = llgamma + &
        (p(1)-1)*log(data(i)/p(2)) - data(i)/p(2) - log(p(2)) - gammln(p(1))
    end do
!   normalization is not 1 in case of cut-offs
    if ( minindx > -1e33 ) then
        if ( maxindx < 1e33 ) then
            llgamma = llgamma - ncur*log(gammp(p(1),maxindx/p(2)) - gammp(p(1),minindx/p(2)))
        else
            llgamma = llgamma - ncur*log(gammq(p(1),minindx/p(2)))
        end if
    else
        if ( maxindx < 1e33 ) then
            llgamma = llgamma - ncur*log(gammp(p(1),maxindx/p(2)))
        end if
    end if
!   minimum, not maximum
    llgamma = -llgamma
!**        print *,'a,b,llgamma = ',p(1),p(2),llgamma

end function llgamma

subroutine dllgamma(p,dp)

!   computes the derivatives of the log-likelihood function for a
!   Gamma distribution with parameters p(1),p(2) and data in common.
!   currently unused.

    implicit none
    real :: p(2),dp(2)
    integer :: i
    real :: p1(2),p2(2),d
    integer,parameter :: nmax=10000000
    integer :: ncur
    real :: data(nmax),xrestrain
    logical :: llwrite,llchangesign
    common /fitdata1/ data
    common /fitdata2/ xrestrain,ncur,llwrite,llchangesign
    real,external :: dgammln,llgamma

    dp(1) = 0
    dp(2) = 0
    do i=1,ncur
        dp(1) = dp(1) &
        + log(data(i)/p(2)) - dgammln(p(1))
        dp(2) = dp(2) &
        - p(1)/p(2) + data(i)/p(2)**2
    end do
!   minimum, not maximum
    dp(1) = -dp(1)
    dp(2) = -dp(2)

    print *,'dp(1) = ',dp(1)
    p1(1) = p(1) + 1e-3
    p1(2) = p(2)
    p2(1) = p(1) - 1e-3
    p2(2) = p(2)
    d = (llgamma(p1)-llgamma(p2))/2e-3
    print *,'    cmp ',d
    print *,'dp(2) = ',dp(2)
    p1(1) = p(1)
    p1(2) = p(2) + 10
    p2(1) = p(1)
    p2(2) = p(2) - 10
    d = (llgamma(p1)-llgamma(p2))/20
    print *,'    cmp ',d

end subroutine dllgamma

real function gamdist(x)

!   Gamma distribution, parameters passed in common
!   currently unused

    implicit none
    real :: x
    real :: alpha,beta
    common /gamcom/ alpha,beta
    real :: z,y
    real,external :: gammln
    if ( x < 0 ) then
        write(0,*) 'gamdist: error: cannot evaluate for x&lt;0: ',x
        call exit(-1)
    else if ( x == 0 ) then
        if ( alpha == 1 ) then
            gamdist = 1/beta
            return
        else if ( alpha > 1 ) then
            gamdist = 0
            return
        else
            write(0,*) 'gamdist: error: infinite for x=0, alpha= ' &
            ,alpha
            call exit(-1)
        end if
    end if
    z = x/beta
    y = (alpha-1)*log(z) - z - log(beta) - gammln(alpha)
    gamdist = exp(y)
end function gamdist

real function cumgamm(x)

!   compute the cumulative Gamma probability minus some
!   requested frequency.  Parameters are passed in common

    implicit none
    real :: x
    real :: pc,ac,bc
    common /ccumgamm/ pc,ac,bc
    real,external :: gammp

    if ( x <= 0 ) then
        cumgamm = -pc
    else
        cumgamm = gammp(ac,x/bc) - pc
    end if

end function cumgamm

real function invcumgamm(p,a,b)

!       compute the inverse of the cumulative Gamma distribution P(a,x/b)
!       as long as I do not find or make an explicit function just solve
!       the equation.

    implicit none
    real :: p,a,b
    integer :: i
    real :: x,x1,x2,tol
    real :: pc,ac,bc
    common /ccumgamm/ pc,ac,bc
    real,external :: cumgamm,zbrent

!   check argument

    if ( p == 0 ) then
        x = 0
        goto 999
    else if ( p < 0 .or. p >= 1 ) then
        write(0,*) 'invcumgamm: illegal argument ',p
        x = 3e33
        goto 999
    end if

!   parameters for function cumgamm

    ac = a
    bc = b
    pc = p

!   bracket zero

    x1 = 0
    x2 = b
    i = 0
100 continue
    if ( cumgamm(x2) < 0 ) then
        x2 = 1.6*x2
        i = i + 1
        if ( i > 100 ) then
            write(0,*) 'invcumgamm: error: cannot find root'
            x = 3e33
            goto 999
        end if
        goto 100
    end if

!   get root

    tol = 1e-5*x2
    x = zbrent(cumgamm,x1,x2,tol)

!   finito

999 continue
    invcumgamm = x
    return
end function invcumgamm
