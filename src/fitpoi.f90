subroutine fitpoi(xx,ntot,mean,a,j1,j2,lweb,ntype &
    ,lchangesign,year,xyear,t,t25,t975,tx &
    ,tx25,tx975,confidenceinterval,lboot,lprint,lwrite)

!   a fit a Poisson distribution to the data

    use BrentToGSL
    implicit none

    integer,intent(in) :: ntot,j1,j2,ntype,year
    real,intent(in)    :: xx(ntot),mean,xyear,confidenceinterval
    real,intent(out)   :: a,t(10),t25(10),t975(10),tx,tx25,tx975
    logical,intent(in) :: lweb,lchangesign,lboot,lprint,lwrite

    integer :: nmc,i,j,nx,iter,init,iens,nzero
    real :: ax,bx,cx,fa,fb,fc,tol
    real,allocatable :: aa(:),tt(:,:),txtx(:)
    real :: x,t5(10),t1(10),db,f,threshold,thens,z,ll,ll1,a25,a975,ranf,xn
    character :: lgt*4

    integer :: ncur
    integer,parameter :: nmax=100000
    real :: data(nmax),restrain
    logical :: llwrite,llchangesign
    common /fitdata1/ data
    common /fitdata2/ restrain,ncur,llwrite,llchangesign

    real,external :: llpoisson,poireturnlevel,poireturnyear

    if ( lwrite ) then
        print *,'fitpoi: input:'
        print *,'j1,j2      = ',j1,j2
        print *,'mean       = ',mean
        print *,'year,xyear = ',year,xyear
        if ( .false. ) then
            do i=1,ntot
                print *,i,xx(i)
            end do
        end if
    end if

!   check input

    do i=1,ntot
        if ( abs(xx(i)-nint(xx(i))) > 0.01 ) then
            write(0,*) 'histogram: error: cannot fit Poisson to non-integer data: ',xx(i)
            call exit(-1)
        endif
    enddo

!   determine number of bootstrap samples needed, demand at least 25 samples above the threshold

    nmc = max(1000,nint(25*2/(1-confidenceinterval/100)))
    allocate(aa(nmc),tt(nmc,10),txtx(nmc))

!   copy to common for routine llpoisson

    ncur = ntot
    do i=1,ncur
        data(i) = xx(i)
    enddo
    llwrite = lwrite

!   fit, using Numerical Recipes routines

    ax = mean/2
    bx = 1.5*mean
    call mnbrak(ax,bx,cx,fa,fb,fc,llpoisson)
    tol = 1e-5
    fa = brent(ax,bx,cx,llpoisson,tol,a)
    do i=1,10
        if ( mod(i,3) == 1 ) then
            x = 1 + i/3
        else if ( mod(i,3) == 2 ) then
            x = 1 + log10(2.) + i/3
        else
            x = log10(5.) + i/3
        end if
        x = x + log10(real(j2-j1+1))
        t(i) = poireturnlevel(a,x)
    end do

    if ( xyear < 1e33 ) then
        tx = poireturnyear(a,xyear)/(j2-j1+1)
        if ( lwrite ) then
            print *,'return time = ',tx
        end if
    end if

!   bootstrap to find error estimates

    if ( .not. lboot ) then
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
        do i=1,ntot
            call random_number(ranf)
            j = 1+int(ntot*ranf)
            if ( j < 1 .or. j > ntot ) then
                write(0,*) 'fitpoi: error: j = ',j
                call exit(-1)
            end if
            data(i) = xx(j)
        end do
        ax = a/2
        bx = 2*a
        call mnbrak(ax,bx,cx,fa,fb,fc,llpoisson)
        tol = 1e-5
        fa = brent(ax,bx,cx,llpoisson,tol,a)
        aa(iens) = a
        do i=1,10
            if ( mod(i,3) == 1 ) then
                x = 1 + i/3
            else if ( mod(i,3) == 2 ) then
                x = 1 + log10(2.) + i/3
            else
                x = log10(5.) + i/3
            end if
            x = x + log10(real(j2-j1+1))
            tt(iens,i) = poireturnlevel(aa(iens),x)
        end do
        if ( xyear < 1e33 ) then
            txtx(iens) = poireturnyear(aa(iens),xyear)
            txtx(iens) = txtx(iens)/(j2-j1+1)
        end if
    end do
    call getcut( a25,(100-confidenceinterval)/2,nmc,aa)
    call getcut(a975,(100+confidenceinterval)/2,nmc,aa)
    do i=1,10
        lgt = '&gt;'
        call getcut(t5(i),95.,nmc,tt(1,i))
        call getcut(t1(i),99.,nmc,tt(1,i))
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
        print '(a)','# <tr><td colspan="4">Fitted to Poisson distribution '// &
                'p(n) = &mu;^n exp(-&mu;)/n!</td></tr>'
        print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td>&mu;</td><td>',a,'</td><td>',a25,'...',a975,'</td?<tr>'
    else
        print '(a)','# Fitted to Poisson distribution'
        print '(a)','# p(n) = mu^n*exp(-mu)/n! with'
        print '(a,f16.3,a,f16.3,a)','# mu = ',a,' \\pm ',a975-a25
    end if
    call printreturnvalue(ntype,t,t25,t975,lweb)
    call printreturntime(year,xyear,tx,tx25,tx975,lweb)
    call plotreturnvalue(ntype,t25,t975,j2-j1+1)
end subroutine fitpoi

real function llpoisson(a)

!   computes the log-likelihood function for a Poisson distribution
!   with parameter a and data in common.

    implicit none

    real :: a

    integer :: i
    real :: s

    integer :: ncur
    integer,parameter :: nmax=100000
    real :: data(nmax),restrain
    logical :: llwrite,llchangesign
    common /fitdata1/ data
    common /fitdata2/ restrain,ncur,llwrite,llchangesign

    real,external :: gammln

    llpoisson = 0
    if ( a <= 0 ) then
        llpoisson = 3e33
        return
    end if
    do i=1,ncur
        llpoisson = llpoisson + data(i)*log(a) - a - gammln(1+data(i))
    enddo
!   normalization is not 1 in case of cut-offs
    call poisnorm(a,s)
    llpoisson = llpoisson - ncur*log(s)
!   minimum, not maximum
    llpoisson = -llpoisson
    if ( llwrite .and. .false. ) print *,'a,llpoisson = ',a,llpoisson
end function llpoisson

subroutine poisnorm(a,s)

!   get the normalization of the poisson function on the interval, minindx,maxindx

    implicit none
    include 'getopts.inc'
    real :: a,s
    integer :: i
    real :: xk1,xk2
    real,external :: poisson,gammp,gammq

    if ( minindx > 0 ) then
        xk1 = nint(minindx+0.5)
        if ( maxindx < 2.**31 ) then
            xk2 = nint(maxindx-0.5) + 1
            s = gammq(xk2,a) - gammq(xk1,a)
        else
            s = gammp(xk1,a)
        endif
    else
        if ( maxindx < 1e33 ) then
            xk2 = nint(maxindx-0.5) + 1
            s = gammq(xk2,a)
        else
            s = 1
        endif
    endif
!**        print *,'poisnorm: norm = ',s
end subroutine poisnorm

real function poisson(a,n)
    implicit none
    real :: a
    integer :: n
    real :: factrl
    external factrl

    poisson = a**n*exp(-a)/factrl(n)
!**        print *,'poisson(',a,n,') = ',poisson
end function poisson

real function cumpois(x)

!   compute the cumulative Poisson probability minus some
!   requested frequency.  Parameters are passed in common

    implicit none
    real :: x
    real :: pc,muc,p1c
    common /ccumpois/ pc,muc,p1c
    real,external :: gammp,gammq

    if ( x < 0 ) then
        cumpois = -pc
    else if ( pc < 0.5 ) then
        cumpois = gammq(x+1,muc) - pc
    else
        cumpois = p1c - gammp(x+1,muc)
    endif
    
end function cumpois

real function invcumpois(p,p1,mu)

!   compute the inverse of the cumulative Poisson distribution Q(k,mu)
!   as long as I do not find or make an explicit function just solve
!   the equation.

    use ZbrentToGSL
    implicit none
    real :: p,p1,mu
    integer :: i
    real :: x,x1,x2,tol
    real :: pc,muc,p1c
    common /ccumpois/ pc,muc,p1c
    real,external :: cumpois

!   check argument

    if ( p == 0 ) then
        x = 0
        goto 999
    elseif ( p < 0 .or. p1 < 0 ) then
        write(0,*) 'invcumpois: illegal argument ',p
        x = 3e33
        goto 999
    endif
    if ( abs(p+p1-1) > 0.0001 ) then
        write(0,*) 'invcumpois: error: p1 /= 1-p ',p1,1-p
        call exit(-1)
    end if

!   parameters for function cumpois

    muc = mu
    pc = p
    p1c = p1

!   bracket zero

    x1 = 0
    x2 = mu
    i = 0
100 continue
    if ( cumpois(x2) < 0 ) then
        x2 = 1.6*x2
        i = i + 1
        if ( i > 100 ) then
            write(0,*) 'invcumpois: error: cannot find root'
            x = 3e33
            goto 999
        endif
        goto 100
    endif

!   get root

    tol = 1e-5*x2
    x = zbrent(cumpois,x1,x2,tol)

!   finito

999 continue
    invcumpois = x
    return
end function invcumpois

real function poireturnlevel(a,x)

!   compute return times given the poisson distribution parameters a and x = log10(returntime)

    implicit none
    real,intent(in) :: a,x
    real :: p,p1
    real,external :: invcumpois

    p1 = 10**(-x)
    p = 1 - p1
    poireturnlevel = invcumpois(p,p1,a)
    
end function poireturnlevel

real function poireturnyear(a,xyear)

!   compute the return time of the value xyear with the fitted value

    implicit none
    real,intent(in) :: a,xyear
    real,external :: gammp

    poireturnyear = 1/gammp(xyear+1,a)

end function poireturnyear

