subroutine fitgam(xx,ntot,mean,sd,a,b,j1,j2,lweb,ntype &
    ,lchangesign,year,xyear,t,t25,t975,tx &
    ,tx25,tx975,confidenceinterval,lboot,lprint,lwrite)

!   a fit a gamma distribution to the data

    implicit none

    integer,intent(in) :: ntot,j1,j2,ntype,year
    real,intent(in) :: xx(ntot),mean,sd,xyear,confidenceinterval
    real,intent(out) :: a,b,t(10),t25(10),t975(10),tx,tx25,tx975
    logical,intent(in) :: lweb,lchangesign,lboot,lprint,lwrite

    integer :: nmc,i,j,nx,iter,init,iens,nzero
    real,allocatable :: aa(:),bb(:),xnxn(:),tt(:,:),txtx(:)
    real :: x,t5(10),t1(10),db,f,threshold,thens,z,ll,ll1,a25,a975,b25,b975,xn25,xn975, &
        ranf,xn
    character :: lgt*4

    integer,parameter :: nmax=10000000
    integer :: ncur
    real :: data(nmax),xrestrain
    logical :: llwrite,llchangesign
    common /fitdata1/ data
    common /fitdata2/ xrestrain,ncur,llwrite,llchangesign

    real,external :: llgamma,gamreturnlevel,gamreturnyear

    if ( lwrite ) then
        print *,'fitgam: input:'
        print *,'j1,j2      = ',j1,j2
        print *,'mean, sd   = ',mean,sd
        print *,'year,xyear = ',year,xyear
        if ( .false. ) then
            do i=1,ntot
                print *,i,xx(i)
            end do
        end if
    end if

!   check input

    nzero = 0
    do i=1,ntot
        if ( xx(i) < 0 .neqv. lchangesign ) then
            write(0,*) 'fitgam: error: cannot fit Gamma to negative data: ',xx(i),lchangesign
            call exit(-1)
        end if
        if ( xx(i) == 0 ) nzero = nzero + 1
    end do
    if ( nzero == ntot ) then
        write(0,*) 'fitgam: error: only zero points'
        call exit(-1)
    end if
    xn = 1 - nzero/real(ntot) ! wet day fraction

!   determine number of bootstrap samples needed, demand at least 25 samples above the threshold

    nmc = max(1000,nint(25*2/(1-confidenceinterval/100)))
    allocate(aa(nmc),bb(nmc),xnxn(nmc),tt(nmc,10),txtx(nmc))

!   copy to common for routines llgamma, dllgamma

    ncur = ntot - nzero
    j = 0
    do i=1,ntot
        if ( xx(i) /= 0 ) then
            j = j + 1
            data(j) = xx(i)
        end if
    end do
    if ( j /= ncur ) then
        write(0,*) 'fitgam: internal error: ',j,ncur
    end if
    llwrite = lwrite

    a = mean**2/sd**2-0.05
    b = sd**2/mean    *0.9
    call fit1gam(a,b,iter)
    do i=1,10
        if ( mod(i,3) == 1 ) then
            x = 1 + i/3
        else if ( mod(i,3) == 2 ) then
            x = 1 + log10(2.) + i/3
        else
            x = log10(5.) + i/3
        end if
        x = x + log10(real(j2-j1+1))
        t(i) = gamreturnlevel(a,b,xn,x)
    end do

    if ( xyear < 1e33 ) then
        tx = gamreturnyear(a,b,xyear)/xn
        if ( lwrite ) then
            print *,'return time = ',tx
        end if
    end if

!   bootstrap to find error estimates

    if ( .not. lboot ) then
        if ( lchangesign ) then
            b = -b
        end if
        return
    end if
    if ( lprint ) then
        if ( lweb ) then
            write(0,*) 'Doing a ',nmc,'-member bootstrap to obtain error estimates'
        else
            print '(a,i6,a)','# Doing a ',nmc,'-member bootstrap to obtain error estimates'
        end if
    end if
    
    do iens=1,nmc
        call keepalive1('Computing bootstrap sample ',iens,nmc)
        ncur = 0
        do i=1,ntot
            call random_number(ranf)
            j = 1+int(ntot*ranf)
            if ( j < 1 .or. j > ntot ) then
                write(0,*) 'fitgam: error: j = ',j
                call exit(-1)
            end if
            if ( xx(j) /= 0 ) then
                ncur = ncur + 1
                data(ncur) = xx(j)
            end if
        end do
        aa(iens) = a
        bb(iens) = b
        xnxn(iens) = ncur/real(ntot)
        call fit1gam(aa(iens),bb(iens),iter)
        do i=1,10
            if ( mod(i,3) == 1 ) then
                x = 1 + i/3
            else if ( mod(i,3) == 2 ) then
                x = 1 + log10(2.) + i/3
            else
                x = log10(5.) + i/3
            end if
            x = x + log10(real(j2-j1+1))
            tt(iens,i) = gamreturnlevel(aa(iens),bb(iens),xnxn(iens),x)
            if ( .false. .and. tt(iens,i) < 10 ) then
                print *,'gamreturnlevel(',aa(iens),xnxn(iens),bb(iens),x/xn,' = ',tt(iens,i),iens,i
            end if
        end do
        if ( xyear < 1e33 ) then
            txtx(iens) = gamreturnyear(aa(iens),bb(iens),xyear)/xn
            if ( txtx(iens) < 1 ) then
                print *,'gamreturnyear(',aa(iens),bb(iens),xyear,') = ',txtx(iens),iens
            end if
        end if
    end do
    if ( lchangesign ) then
        b = -b
        bb = -bb
    end if
    call getcut( a25,(100-confidenceinterval)/2,nmc,aa)
    call getcut(a975,(100+confidenceinterval)/2,nmc,aa)
    call getcut( b25,(100-confidenceinterval)/2,nmc,bb)
    call getcut(b975,(100+confidenceinterval)/2,nmc,bb)
    call getcut( xn25,(100-confidenceinterval)/2,nmc,xnxn)
    call getcut(xn975,(100+confidenceinterval)/2,nmc,xnxn)
    do i=1,10
        if ( lchangesign ) then
            lgt = '&lt;'
            call getcut(t5(i),5.,nmc,tt(1,i))
            call getcut(t1(i),1.,nmc,tt(1,i))
        else
            lgt = '&gt;'
            call getcut(t5(i),95.,nmc,tt(1,i))
            call getcut(t1(i),99.,nmc,tt(1,i))
        end if
        call getcut( t25(i),(100-confidenceinterval)/2,nmc,tt(1,i))
        call getcut(t975(i),(100+confidenceinterval)/2,nmc,tt(1,i))
    end do
    if ( xyear < 1e33 ) then
        call getcut( tx25,(100-confidenceinterval)/2,nmc,txtx)
        call getcut(tx975,(100+confidenceinterval)/2,nmc,txtx)
    end if

!   output

    if ( .not. lprint .and. .not. lwrite ) return
    if ( lweb ) then
        if ( nzero == 0 ) then
            print '(a)','# <tr><td colspan="4">Fitted to gamma distribution '// &
                'G(x) = (x/&sigma;)^(&xi;-1)*exp(-x/&sigma;)/(&sigma;&Gamma;(&xi;))</td></tr>'
        else
            print '(a)','# <tr><td colspan="4">Fitted to gamma distribution '// &
                'G''(x) = (1-N)&delta;(x) + N(x/&sigma;)^(&xi;-1)*exp(-x/&sigma;)/(&sigma;&Gamma;(&xi;))</td></tr>'
            print '(a,f16.3,a,f16.3,a,f16.3,a)','#<tr><td> N (wet day fraction)</td><td>',xn, &
                '</td><td>',xn25,'...',xn975,'</td?<tr>'
        end if
        print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td>&sigma; (scale parameter)</td><td>',b,'</td><td>',b25,'...',b975,'</td?<tr>'
        print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td>&xi; (shape parameter)</td><td>',a,'</td><td>',a25,'...',a975,'</td?<tr>'
    else
        print '(a,i5,a)','# Fitted to gamma distribution in ',iter,' iterations'
        if ( nzero == 0 ) then
            print '(a)','# G(x) = (x/sigma)^(xi-1)*exp(-x/sigma)/(sigma*Gamma(xi)) with'
        else
            print '(a)','# G(x) = (1-N)delta(x) + N*(x/sigma)^(xi-1)*exp(-x/sigma)/(sigma*Gamma(xi)) with'    
            print '(a,f16.3,a,f16.3,a)','# N     = ',xn,' \\pm ',xn975-xn25,' fraction wet days'
        end if
        print '(a,f16.3,a,f16.3,a)','# sigma = ',b,' \\pm ',b975-b25,' scale parameter'
        print '(a,f16.3,a,f16.3,a)','# xi    = ',a,' \\pm ',a975-a25,' shape parameter'
    end if
    call printreturnvalue(ntype,t,t25,t975,lweb)
    if ( lchangesign .and. xyear < 1e33 ) then
        call printreturntime(year,-xyear,tx,tx25,tx975,lweb)
    else
        call printreturntime(year,xyear,tx,tx25,tx975,lweb)
    end if
    call plotreturnvalue(ntype,t25,t975,j2-j1+1)
end subroutine fitgam

subroutine fit1gam(a,b,iter)

!   fit, using Numerical Recipes routines

    implicit none
    real,intent(inout) :: a,b
    integer,intent(out) :: iter
    integer :: i
    real :: tol,p(3,2),q(2),y(3)
    real,external :: llgamma

    p(1,1) = a
    p(1,2) = b
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
end subroutine fit1gam

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
        if ( data(i)/p(2) <= 0 ) then
            write(0,*) 'llgamma: error: data(',i,')/p(2) <= 0 ',data(i),p(2)
            call exit(-1)
        end if
        llgamma = llgamma + (p(1)-1)*log(data(i)/p(2)) - data(i)/p(2) - log(abs(p(2))) - gammln(p(1))
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
    !!!print *,'a,b,llgamma,ncur = ',p(1),p(2),llgamma,ncur

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
            write(0,*) 'gamdist: error: infinite for x=0, alpha= ',alpha
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

!   compute the inverse of the cumulative Gamma distribution P(a,x/b)
!   as long as I do not find or make an explicit function just solve
!   the equation.

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

real function gamreturnlevel(a,b,xn,x)

!   compute return times given the gamma distribution parameters a=xi,b=sigma and x = log10(returntime)
!   Should use a complement function and/or a few Taylor series approximations but does not yet

    implicit none
    real,intent(in) :: a,b,xn,x
    integer,save :: init=0
    real :: p
    real,external :: invcumgamm

    if ( b > 0 ) then ! upper tail
        p = 1 - 10**(-x)
    else
        p = 10**(-x)
    end if
    if ( p < 1-xn ) then
        gamreturnlevel = 0
    else
        p = (p+xn-1)/xn
        if ( p > 0.99999 .and. init == 0 ) then
            init = 1
            write(0,*) 'gamreturnlevel: loss of precision, please upgrade algorithms ',p
        end if
        gamreturnlevel = invcumgamm(p,a,abs(b))
    end if
end function gamreturnlevel

real function gamreturnyear(a,b,xyear)

!   compute the return time of the value xyear with the fitted values

    implicit none
    real,intent(in) :: a,b,xyear
    integer,save :: init=0
    real :: p
    real,external :: cumgamm

    real :: pc,ac,bc
    common /ccumgamm/ pc,ac,bc

    pc = 0
    ac = a
    bc = b
    p = cumgamm(xyear)
    if ( p < 0.00001 .and. init == 0 ) then
        init = 1
        write(0,*) 'gamreturnyear: loss of precision, please upgrade algorithms ',p
    end if
    gamreturnyear = 1/(1-p)
end function gamreturnyear

