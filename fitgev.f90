subroutine fitgev(xx,ntot,mean,sd,a,b,xi,j1,j2,lweb,ntype &
    ,lchangesign,year,xyear,t,t25,t975,tx &
    ,tx25,tx975,inrestrain,confidenceinterval,lboot,lprint &
    ,lwrite)

!       a fit a GEV distribution to the data, which is already assumed to be block max
!       input:
!       xx(ntot) data
!       mean,sd  for first guess parameters
!       j1,j2    use days/months/... j1 to j2
!       year     leave out this year from the fit and compute return time for it
!       xyear    value for year, has been set to undef in the series
!       inrestrain restrain xi parameter by adding a normal distribution of width 0.5*inrestrain to the cost function
!       output
!       a,b,xi     parameters of fit
!       t(10)    return values for 10, 20, 50, ..., 10000 years
!       t25,t975   2.5%, 97.5% quantiles of these return values
!       tx      return time of the value of year (xyear) in the context of the other values
!       tx25,tx975 2.5%, 97.5% quantiles of these return times

    implicit none

    integer :: nmc
    parameter(nmc=1000)
    integer :: ntot,j1,j2,ntype,year
    real :: xx(ntot),mean,sd,a,b,xi,xyear,inrestrain, &
        t(10),t25(10),t975(10),tx,tx25,tx975,confidenceinterval
    logical :: lweb,lchangesign,lboot,lprint,lwrite

    integer :: i,j,k,n,nx,iter,iens
    real :: x,bb(nmc),xixi(nmc),tt(nmc,10),b25 &
        ,b975,xi25,xi975,t5(10),t1(10),db,dxi,f &
        ,threshold,thens,z,ll,ll1,txtx(nmc),aa(nmc) &
        ,a25,a975,ranf
    character :: lgt*4

    integer,parameter :: nmax=100000
    integer :: ncur
    real :: data(nmax),restrain
    logical :: llwrite,llchangesign
    common /fitdata1/ data
    common /fitdata2/ restrain,ncur,llwrite,llchangesign

    real,external :: llgev,gevreturnlevel,gevreturnyear

    if ( lwrite ) then
        print *,'fitgev: input:'
        print *,'j1,j2      = ',j1,j2
        print *,'mean, sd   = ',mean,sd
        print *,'year,xyear = ',year,xyear
        if ( .false. ) then
            do i=1,ntot
                print *,i,xx(i)
            end do
        end if
    end if

!   ill-defined case

    if ( sd == 0 ) then
        xi = 3e33
        a = 3e33
        b = 3e33
        t = 3e33
        t25 = 3e33
        t975 = 3e33
        tx = 3e33
        tx25 = 3e33
        tx975 = 3e33
        return
    end if

!   copy to common for routine llgev

    ncur = ntot
    do i=1,ncur
        data(i) = xx(i)
    end do
    restrain = inrestrain
    llwrite = lwrite
            
    b = sd*sqrt(6.)/(4*atan(1.))
    a = mean - 0.57721*b
    xi = 0
    call fit1gev(a,b,xi,iter)

    do i=1,10
        if ( mod(i,3) == 1 ) then
            x = 1 + i/3
        else if ( mod(i,3) == 2 ) then
            x = 1 + log10(2.) + i/3
        else
            x = log10(5.) + i/3
        end if
        x = x + log10(real(j2-j1+1))
        t(i) = gevreturnlevel(a,b,xi,x)
    end do

    if ( xyear < 1e33 ) then
        tx = gevreturnyear(a,b,xi,xyear)
        if ( lwrite ) then
            print *,'return time = ',tx
        end if
    end if

!   bootstrap to find error estimates

    if ( .not. lboot ) then
        if ( lchangesign ) then
            b = -b
            t = -t
        end if
        return
    end if
    if ( .not. lweb ) print '(a,i6,a)','# doing a ',nmc,'-member bootstrap to obtain error estimates'
    do iens=1,nmc
        call keepalive1('Bootstrapping',iens,nmc)
        if ( .not. lweb .and. mod(iens,100) == 0 ) print '(a,i6)','# ',iens
        do i=1,ntot
            call random_number(ranf)
            j = 1+int(ntot*ranf)
            if ( j < 1 .or. j > ntot ) then
                write(0,*) 'fitgev: error: j = ',j
                call exit(-1)
            end if
            data(i) = xx(j)
        end do
        aa(iens) = a
        bb(iens) = b
        xixi(iens) = xi
        llwrite = .false. 
        call fit1gev(aa(iens),bb(iens),xixi(iens),iter)
        do i=1,10
            if ( mod(i,3) == 1 ) then
                x = 1 + i/3
            else if ( mod(i,3) == 2 ) then
                x = 1 + log10(2.) + i/3
            else
                x = log10(5.) + i/3
            end if
            x = x + log10(real(j2-j1+1))
            tt(iens,i) = gevreturnlevel(aa(iens),bb(iens),xixi(iens),x)
            if ( .false. .and. tt(iens,i) < 10 ) then
                print *,'gevreturnlevel(',aa(iens),bb(iens),xixi(iens),x,' = ',tt(iens,i),iens,i
            end if
        end do
        if ( xyear < 1e33 ) then
            txtx(iens) = gevreturnyear(aa(iens),bb(iens),xixi(iens),xyear)
            if ( .false. .and. txtx(iens) < 1 ) then
                print *,'gevreturnyear(',aa(iens),bb(iens),xixi(iens),xyear,') = ',txtx(iens),iens
            end if
        end if
    end do
    if ( lchangesign ) then
        a = -a
        aa = -aa
        b = -b
        t = -t
        tt = -tt
    end if
    call getcut( a25,(100-confidenceinterval)/2,nmc,aa)
    call getcut(a975,(100+confidenceinterval)/2,nmc,aa)
    call getcut( b25,(100-confidenceinterval)/2,nmc,bb)
    call getcut(b975,(100+confidenceinterval)/2,nmc,bb)
    call getcut( xi25,(100-confidenceinterval)/2,nmc,xixi)
    call getcut(xi975,(100+confidenceinterval)/2,nmc,xixi)
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

!       output

    if ( .not. lprint .and. .not. lwrite ) return
    if ( lweb ) then
        print '(a)','# <tr><td colspan="3">Fitted to GEV '// &
            'distribution P(x) = exp(-(1+&xi;*(x-a)/b)**(-1/&xi;))</td></tr>'
        print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td>a:</td><td>' &
            ,a,'</td><td>',a25,'...',a975,'</td></tr>'
        print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td>b:</td><td>' &
            ,abs(b),'</td><td>',b25,'...',b975,'</td></tr>'
        print '(a,f16.3,a,f16.3,a,f16.3,a)' &
            ,'# <tr><td>&xi;:</td><td>',xi,'</td><td>',xi25,'...',xi975,'</td></tr>'
    else
        print '(a,i5,a)','# Fitted to GEV distribution in ',iter,' iterations'
        print '(a)','# P(x) = exp(-(1+xi*(x-a)/b)**(-1/xi)) with'
        print '(a,f16.3,a,f16.3,a,f16.3)','# a = ',a,' \\pm ',a975-a25
        print '(a,f16.3,a,f16.3,a,f16.3)','# b = ',b,' \\pm ',b975-b25
        print '(a,f16.3,a,f16.3,a,f16.3)','# xi= ',xi,' \\pm ',xi975-xi25
    end if
    call printreturnvalue(ntype,t,t25,t975,lweb)
    if ( lchangesign .and. xyear < 1e33 ) xyear = -xyear
    call printreturntime(year,xyear,tx,tx25,tx975,lweb)
    if ( lchangesign .and. xyear < 1e33 ) xyear = -xyear
    call plotreturnvalue(ntype,t25,t975,j2-j1+1)
end subroutine fitgev

subroutine fit1gev(a,b,xi,iter)
    use AmoebaToGSL
    implicit none
    integer :: iter
    real :: a,b,xi
    integer :: i
    real :: q(3),p(4,3),y(4),tol
    real,external :: llgev

    q(1) = a
    q(2) = b
    q(3) = xi
    p(1,1) = q(1) *0.9
    p(1,2) = q(2) *0.9
    p(1,3) = q(3) *0.9
    p(2,1) = p(1,1) *1.2
    p(2,2) = p(1,2)
    p(2,3) = p(1,3)
    p(3,1) = p(1,1)
    p(3,2) = p(1,2) *1.2
    p(3,3) = p(1,3)
    p(4,1) = p(1,1)
    p(4,2) = p(1,2)
    p(4,3) = p(1,3) *1.2 + 0.1
    do i=1,4
        q(1) = p(i,1)
        q(2) = p(i,2)
        q(3) = p(i,3)
        y(i) = llgev(q)
    end do
    tol = 1e-4
    call amoeba(p,y,4,3,3,tol,llgev,iter)
!   maybe add restart later
    a = p(1,1)
    b = abs(p(1,2))
    xi = p(1,3)
end subroutine fit1gev

real function llgev(p)

!   computes the log-likelihood function for a GEV distribution
!   with parameters a,b,xi=p(1-3) and data in common.

    implicit none

    real,intent(in) :: p(3)
    integer :: i
    real :: z,xi,s

    integer,parameter :: nmax=100000
    integer :: ncur
    real :: data(nmax),restrain
    logical :: llwrite,llchangesign
    common /fitdata1/ data
    common /fitdata2/ restrain,ncur,llwrite,llchangesign

    llgev = 0
    if ( p(2) <= 0 ) then
        llgev = 3e33
        goto 999
    end if
    if ( abs(p(3)) > 10 ) then
        llgev = 3e33
        goto 999
    end if
    if ( restrain < 0 ) then
        write(0,*) 'llgev: restrain<0 ',restrain
        call exit(-1)
    end if
    do i=1,ncur
        if ( abs(p(2)) < 1e-30 ) then
            llgev = 3e33
            goto 999
        end if
        z = (data(i)-p(1))/abs(p(2))
        xi = p(3)
        if ( abs(xi) < 1e-4 ) then
            if ( -z+xi*z**2/2 > log(3e33) ) then
                llgev = 3e33
                goto 999
            end if
            llgev = llgev - exp(-z+xi*z**2/2) - z*(1+xi-xi*z/2)
        else
            if ( 1+xi*z <= 0 ) then
                llgev = 3e33
                goto 999
            else if ( -log(1+xi*z)/xi > log(3e33) ) then
            ! too large...
                llgev = 3e33
                goto 999
            else
                llgev = llgev - (1+1/xi)*log(1+xi*z) &
                - (1+xi*z)**(-1/xi)
            end if
        end if
    end do
    llgev = llgev - ncur*log(abs(p(2)))
!   normalization is not 1 in case of cut-offs
    call gevnorm(p(1),abs(p(2)),p(3),s)
    if ( s < 1e33 ) then
        llgev = llgev - ncur*log(s)
    else
        llgev = 3e33
        goto 999
    end if
    if ( restrain /= 0 ) then
!       preconditioning on xi with gaussian of width restrain/2 around 0
        llgev = llgev - (xi/(restrain/2))**2/2
    end if
!   minimum, not maximum
    llgev = -llgev

999 continue
    if ( llwrite ) print *,'ncur,a,b,xi,llgev = ',ncur,p(1),p(2),p(3),llgev,z
end function llgev

subroutine gevnorm(a,b,xi,s)
    implicit none
    include 'getopts.inc'
    real :: a,b,xi,s
    real :: z1,z2

    if ( minindx > -1e33 ) then
        z1 = (minindx-a)/b
        if ( maxindx < 1e33 ) then
            z2 = (maxindx-a)/b
            if ( abs(xi) < 1e-4 ) then
                s = exp(-exp(-z2)) - exp(-exp(-z1))
            else
                s = exp(-(1+xi*z2)**(-1/xi)) - exp(-(1+xi*z1)**(-1/xi))
            end if
        else
            if ( abs(xi) < 1e-4 ) then
                s = 1 - exp(-exp(-z1))
            else
                z2 = (1+xi*z1)
                if ( z2 > 0 ) then
                    s = 1 - exp(-z2**(-1/xi))
                    if ( s < 1e-4 ) s = 3e33
                else
                    s = 3e33
                end if
            end if
        end if
    else
        if ( maxindx < 1e33 ) then
            z2 = (maxindx-a)/b
            if ( xi == 0 ) then
                s = exp(-exp(-z2))
            else
                s = exp(-(1+xi*z2)**(-1/xi))
            end if
        else
            s = 1
        end if
    end if
!**        print *,'gevnorm: norm = ',a,b,s
end subroutine gevnorm

real function gevreturnlevel(a,b,xi,x)

!   compute return times given the GEV distribution parameters a,b,xi and x = log10(returntime)
!   Uses a few Taylor series approximation for xi small and/or return time large

    implicit none
    real,intent(in) :: a,b,xi,x
    real :: y,t
    if ( abs(xi) > 10 ) then
        gevreturnlevel = 3e33
    else if ( abs(xi) < 1e-4 ) then
        if ( x <= 8 ) then
            y = log(-log(1-dble(10)**(dble(-x))))
        else
            y = -x*log(10.)
        end if
        t = a - b*y + b*xi/2*y**2
    else
        if ( x <= 8 ) then
            t = a - b/xi*(1-(-log(1-dble(10)**(dble(-x))))**(-xi))
        else
            t = a - b/xi*(1-10.**(x*xi))
        end if
    end if
    gevreturnlevel = t
end function gevreturnlevel

real function gevreturnyear(a,b,xi,xyear)

!   compute the return time of the value xyear with the fitted values

    implicit none
    real,intent(in) :: a,b,xi,xyear
    real :: y,z,tx

    z = (1 + xi*(xyear-a)/b)
    if ( z < 0 ) then
        y = 1e20
    else if ( abs(xi) > 1e-3 ) then
        y = -z**(-1/xi)
    else
        y = -exp(-(xyear-a)/b + xi/2*((xyear-a)/b)**2)
    end if
    if ( y >= 1e20 ) then
        tx = 1e20
    else if ( abs(y) > 1e-3 ) then
        tx = 1/(1 - exp(y))
    else if ( abs(y) > 1e-19 ) then
        tx = -1/(y + 0.5*y**2)
    else
        tx = 1e20
    end if
    if ( .false. .and. tx > 1e20 ) then
        write(0,*) 'gevreturnyear: tx > 1e20: ',tx
        write(0,*) 'b,xi,xyear = ',b,xi,xyear
    end if
    gevreturnyear = tx
end function gevreturnyear
