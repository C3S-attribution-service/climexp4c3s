subroutine fitgum(xx,ntot,mean,sd,a,b,j1,j2,lweb,ntype &
    ,lchangesign,year,xyear,t,t25,t975,tx &
    ,tx25,tx975,confidenceinterval,lboot,lprint,lwrite)

!   a fit a Gumbel distribution to the data, which is already assumed to be block max
!   input:
!   xx(ntot) data
!   mean,sd  for first guess parameters
!   j1,j2    use days/months/... j1 to j2
!   output
!   a,b      parameters of fit
!   t(10)    return values for 10, 20, 50, ..., 10000 years

    implicit none

    integer,intent(in) :: ntot,j1,j2,ntype,year
    real,intent(in) :: xx(ntot),mean,sd,xyear,confidenceinterval
    real,intent(out) :: a,b,t(10),t25(10),t975(10),tx,tx25,tx975
    logical,intent(in) :: lweb,lchangesign,lboot,lprint,lwrite

    integer :: nmc,i,j,nx,iter,iens
    real,allocatable :: aa(:),bb(:),xixi(:),tt(:,:),txtx(:)
    real :: x,xi,b25,b975,t5(10),t1(10),db,f,threshold,thens,z,ll,ll1,a25,a975,ranf
    character :: lgt*4

    integer :: ncur
    integer,parameter :: nmax=100000
    real :: data(nmax),restrain
    logical :: llwrite,llchangesign
    common /fitdata1/ data
    common /fitdata2/ restrain,ncur,llwrite,llchangesign

    real,external :: llgumbel,gevreturnlevel,gevreturnyear

    if ( lwrite ) then
        print *,'fitgum: input:'
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
        b = 3e33
        t = 3e33
        t25 = 3e33
        t975 = 3e33
        tx = 3e33
        tx25 = 3e33
        tx975 = 3e33
        return
    end if

!   determine number of bootstrap samples needed, demand at least 25 samples above the threshold

    nmc = max(1000,nint(25*2/(1-confidenceinterval/100)))
    allocate(aa(nmc),bb(nmc),xixi(nmc),tt(nmc,10),txtx(nmc))

!   copy to common for routine llgumbel

    ncur = ntot
    do i=1,ncur
        data(i) = xx(i)
    end do
    llwrite = lwrite

    b = sd*sqrt(6.)/(4*atan(1.))
    a = mean - 0.57721*b
    call fit1gum(a,b,iter)
    do i=1,10
        if ( mod(i,3) == 1 ) then
            x = 1 + i/3
        else if ( mod(i,3) == 2 ) then
            x = 1 + log10(2.) + i/3
        else
            x = log10(5.) + i/3
        end if
        x = x + log10(real(j2-j1+1))
        xi = 0
        t(i) = gevreturnlevel(a,b,xi,x)
    end do

    if ( xyear < 1e33 ) then
        xi = 0
        tx = gevreturnyear(a,b,xi,xyear)
        if ( lwrite ) then
            print *,'return time = ',tx
        end if
    end if

!   bootstrap to find error estimates

    if ( .not. lboot ) then
        if ( lchangesign ) then
            a = -a
            b = -b
            t = -t
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
        do i=1,ncur
            call random_number(ranf)
            j = 1+int(ntot*ranf)
            if ( j < 1 .or. j > ncur ) then
                write(0,*) 'fitgum: error: j = ',j
                call exit(-1)
            end if
            data(i) = xx(j)
        end do
        aa(iens) = a
        bb(iens) = b
        call fit1gum(aa(iens),bb(iens),iter)
        do i=1,10
            if ( mod(i,3) == 1 ) then
                x = 1 + i/3
            else if ( mod(i,3) == 2 ) then
                x = 1 + log10(2.) + i/3
            else
                x = log10(5.) + i/3
            end if
            x = x + log10(real(j2-j1+1))
            xixi(iens) = 0
            tt(iens,i) = gevreturnlevel(aa(iens),bb(iens),xixi(iens),x)
            if ( .false. .and. tt(iens,i) < 10 ) then
                print *,'gevreturnlevel(',aa(iens),bb(iens),xixi(iens),x,' = ',tt(iens,i),iens,i
            end if
        end do
        if ( xyear < 1e33 ) then
            txtx(iens) = gevreturnyear(aa(iens),bb(iens),xixi(iens),xyear)
            if ( txtx(iens) < 1 ) then
                print *,'gevreturnyear(',aa(iens),bb(iens),xixi(iens),xyear,') = ',txtx(iens),iens
            end if
        end if
    end do
    if ( lchangesign ) then
        a = -a
        aa = -aa
        t = -t
        tt = -tt
    end if
    call getcut( a25,(100-confidenceinterval)/2,nmc,aa)
    call getcut(a975,(100+confidenceinterval)/2,nmc,aa)
    call getcut( b25,(100-confidenceinterval)/2,nmc,bb)
    call getcut(b975,(100+confidenceinterval)/2,nmc,bb)
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
        print '(a)','# <tr><td colspan="3">Fitted to '// &
            'Gumbel distribution P(x) = exp(-exp(-(x-&mu;)/'// &
            '&sigma;))</td></tr>'
        print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td>&mu;:'// &
            '</td><td>',a,'</td><td>',a25,'...',a975,'</td></tr>'
        print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td>&sigma;:'// &
            '</td><td>',b,'</td><td>',b25,'...',b975,'</td></tr>'
    else
        print '(a,i5,a)','# Fitted to Gumbel distribution in ',iter,' iterations'
        print '(a)','# P(x) = exp(-exp(-(x-a)/b)) with'
        print '(a,f16.3,a,f16.3,a,f16.3)','# a = ',a,' \\pm ',a975-a25
        print '(a,f16.3,a,f16.3,a,f16.3)','# b = ',b,' \\pm ',b975-b25
    end if
    call printreturnvalue(ntype,t,t25,t975,lweb)
    call printreturntime(year,xyear,tx,tx25,tx975,lweb)
    call plotreturnvalue(ntype,t25,t975,j2-j1+1)
end subroutine fitgum

subroutine fit1gum(a,b,iter)

!   fit, using Numerical Recipes routines

    use AmoebaToGSL
    implicit none
    real,intent(inout) :: a,b
    integer,intent(out) :: iter
    integer :: i
    real :: q(2),p(3,2),y(3),tol
    real,external :: llgumbel
    q(1) = a
    q(2) = b
    p(1,1) = q(1) *0.9
    p(1,2) = q(2) *0.9
    p(2,1) = p(1,1) *1.2
    p(2,2) = p(1,2)
    p(3,1) = p(1,1)
    p(3,2) = p(1,2) *1.2
    do i=1,3
        q(1) = p(i,1)
        q(2) = p(i,2)
        y(i) = llgumbel(q)
    end do
    tol = 1e-4
    call amoeba(p,y,3,2,2,tol,llgumbel,iter)
!   maybe add restart later
    a = p(1,1)
    b = abs(p(1,2))
end subroutine fit1gum

real function llgumbel(p)

!   computes the log-likelihood function for a Gumbel distribution
!   with parameters a,b=p(1),p(2) and data in common.

    implicit none
    real,intent(in) :: p(2)
    integer :: i
    real :: z,s

    integer :: ncur
    integer,parameter :: nmax=100000
    real :: data(nmax),restrain
    logical :: llwrite,llchangesign
    common /fitdata1/ data
    common /fitdata2/ restrain,ncur,llwrite,llchangesign

    llgumbel = 0
    do i=1,ncur
        z = (data(i)-p(1))/abs(p(2))
        llgumbel = llgumbel - exp(-z) - z
    end do
    llgumbel = llgumbel - ncur*log(abs(p(2)))
!   normalization is not 1 in case of cut-offs
    call gumbnorm(p(1),abs(p(2)),s)
    llgumbel = llgumbel - ncur*log(s)
!   minimum, not maximum
    llgumbel = -llgumbel
!**        print *,'a,b,llgumbel = ',p(1),p(2),llgumbel

end function llgumbel

subroutine gumbnorm(a,b,s)
    implicit none
    include 'getopts.inc'
    real,intent(in) :: a,b
    real,intent(out) :: s
    real :: z1,z2
    if ( minindx > -1e33 ) then
        z1 = (minindx-a)/b
        if ( maxindx < 1e33 ) then
            z2 = (maxindx-a)/b
            s = exp(-exp(-z2)) - exp(-exp(-z1))
        else
            s = 1 - exp(-exp(-z1))
        end if
    else
        if ( maxindx < 1e33 ) then
            z2 = (maxindx-a)/b
            s = exp(-exp(-z2))
        else
            s = 1
        end if
    end if
!**        print *,'gumbnorm: norm = ',a,b,s
end subroutine gumbnorm