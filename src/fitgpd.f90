subroutine fitgpd(xx,ntot,mean,sd,b,xi,j1,j2,lweb,ntype &
    ,lchangesign,pthreshold,threshold,year,xyear,t,t25,t975,tx &
    ,tx25,tx975,inrestrain,assume,confidenceinterval,lboot &
    ,lprint,lwrite)

!   a fit a GPD distribution to the data
!   input:
!       xx(ntot) data
!       mean,sd  for first guess parameters
!       j1,j2    use days/months/... j1 to j2
!   output
!       b,xi     parameters of fit
!       t(10)    return values for 10, 20, 50, ..., 10000 years

    implicit none

    integer :: nmc
    parameter(nmc=1000)
    integer :: ntot,j1,j2,ntype,year
    real :: xx(ntot),mean,sd,b,pthreshold,inrestrain,xyear, &
        t(10),t25(10),t975(10),tx,tx25,tx975,confidenceinterval
    logical :: lweb,lchangesign,lboot,lprint,lwrite
    character assume*5

    integer :: i,j,k,n,nx,iter,iens
    real :: x,xi,bb(nmc),xixi(nmc),tt(nmc,10),b25 &
        ,b975,xi25,xi975,t5(10),t1(10),db,dxi,f &
        ,threshold,thens,z,ll,ll1,txtx(nmc),ranf
    character :: lgt*4
    integer,save :: iweird=0,nweird=1

    integer :: ncur
    integer,parameter :: nmax=100000
    real :: data(nmax),restrain,cthreshold
    logical :: llwrite,llchangesign
    common /fitdata1/ data
    common /fitdata2/ restrain,ncur,llwrite,llchangesign,cthreshold
    character cassume*5
    common /fitdata4/ cassume

    real,external :: llgpd

    if ( ntot > nmax ) then
        write(0,*) 'fitgpd: error: recompile with nmax at least ',ntot
        call exit(-1)
    end if
    llwrite = lwrite
    if ( lwrite ) then
        print *,'fitgpd: input:'
        print *,'j1,j2      = ',j1,j2
        print *,'pthreshold = ',pthreshold,'%'
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
        b = 3e33
        t = 3e33
        t25 = 3e33
        t975 = 3e33
        tx = 3e33
        tx25 = 3e33
        tx975 = 3e33
        return
    end if

!   copy to common for routine llgpd

    if ( pthreshold > 1e33 .or. pthreshold == 0 ) then
        write(*,'(a)') '# histogram: please specify threshold when fitting GPD to tail'
        call exit(-1)
    end if
    call getcut(threshold,pthreshold,ntot,xx)
    ncur = 0
    do i=1,ntot
        if ( xx(i) >= threshold ) then
            ncur = ncur + 1
            data(ncur) = xx(i) - threshold
        end if
    end do
    restrain = inrestrain
    cthreshold = threshold
    cassume = assume
    llchangesign = lchangesign
    if ( lwrite ) print *,'fitgpd: found ',ncur,' points above threshold'
!   make sure that points with a lot of equal values at the
!   threshold are not included...
    if ( ncur < 3 .or. ncur > 2*(1-pthreshold/100)*ntot ) then
        xi = 3e33
        b = 3e33
        t = 3e33
        t25 = 3e33
        t975 = 3e33
        tx = 3e33
        tx25 = 3e33
        tx975 = 3e33
        return
    end if
    b = sd
    if ( lchangesign .and. cassume == 'scale' ) then
        xi = -0.1
    else
        xi = 0.1
    end if
    call fit1gpd(b,xi,iter)
    if ( b < 1e-6*sd ) then
    !           something went wrong, throw away results
        iweird = iweird + 1
        if ( iweird >= nweird ) then
            nweird = 2*nweird
            print '(a,2g16.4,2i8)','# fitgpd: weird result for b = ',b,xi,ncur,iweird
        end if
        xi = 3e33
        b = 3e33
        t = 3e33
        t25 = 3e33
        t975 = 3e33
        tx = 3e33
        tx25 = 3e33
        tx975 = 3e33
        return
    end if
    do i=1,10
        if ( mod(i,3) == 1 ) then
            x = 1 + i/3
        else if ( mod(i,3) == 2 ) then
            x = 1 + log10(2.) + i/3
        else
            x = log10(5.) + i/3
        end if
        x = x + log10(real(j2-j1+1)) + log10(1-pthreshold/100)
        if ( abs(xi) < 1e-4 ) then
            t(i) = b*x*log(10.) + 0.5*xi*(x*log(10.))**2
        else
            t(i) = b/xi*(-1 + exp(xi*x*log(10.)))
        end if
        t(i) = t(i) + threshold
    end do

    if ( xyear < 1e33 ) then
        if ( xyear > threshold ) then
            x =  xyear - threshold
            z = (1 + xi*x/b)
            if ( z > 1e10 ) then
                tx = 3e33
            else if ( z <= 0 ) then
                tx = 1e20 ! infinity
            else if ( abs(xi) < 1e-3 ) then
                tx = exp((x/b) - 0.5*xi*(x/b)**2)/(1-pthreshold/100)
            else if ( log(z)/xi > 45 ) then
                tx = 1e20
            else
                tx = z**(1/xi)/(1-pthreshold/100)
            end if
            if ( .false. .and. tx > 1e20 ) then
                write(0,*) 'fitgpd: tx > 1e20: ',tx
                write(0,*) 'z,xi = ',z,xi
            end if
        else
            n = 0
            do i=1,ntot
                if ( xx(i) > xyear ) n = n + 1
            end do
            tx = real(ntot+1)/real(n)
        end if
        ! convert to years
        tx = tx/(j2-j1+1)
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
        if ( .not. lweb .and. mod(iens,100) == 0 ) &
        print '(a,i6)','# ',iens
        do i=1,ntot
            call random_number(ranf)
            j = 1+int(ntot*ranf)
            if ( j < 1 .or. j > ntot ) then
                write(0,*) 'fitgpd: error: j = ',j
                call exit(-1)
            end if
            data(i) = xx(j)
        end do
        call getcut(thens,pthreshold,ntot,data)
!       as a side effect, data is now sorted
        do i=ntot,1,-1
            data(i) = data(i) - thens
            if ( data(i) < 0 ) exit
        end do
        ncur = ntot - i
        do i=1,ncur
            data(i) = data(i+ntot-ncur)
        end do
        bb(iens) = b
        xixi(iens) = 0
        call fit1gpd(bb(iens),xixi(iens),iter)
        do i=1,10
            if ( mod(i,3) == 1 ) then
                x = 1 + i/3
            else if ( mod(i,3) == 2 ) then
                x = 1 + log10(2.) + i/3
            else
                x = log10(5.) + i/3
            end if
            x = x + log10(real(j2-j1+1)) + log10(1-pthreshold/100)
            if ( abs(xixi(iens)) < 1e-4 ) then
                tt(iens,i) = bb(iens)*x*log(10.)
            else
                tt(iens,i) = bb(iens)/xixi(iens)*(-1 + exp(xixi(iens)*x*log(10.)))
            end if
            tt(iens,i) = tt(iens,i) + thens
        end do
        if ( xyear < 1e33 ) then
            if ( xyear > thens ) then
                x =  xyear - thens
                if ( abs(xixi(iens)) < 1e-4 ) then
                    z = x/bb(iens) - 0.5*(x/bb(iens))**2
                    txtx(iens) = exp(z)/(1-pthreshold/100)
                else
                    if ( bb(iens) > 1e-20 ) then
                        z = (1 + xixi(iens)*x/bb(iens))
                    else
                        z = 1e20
                    end if
                    if ( z > 0 ) then
                        z = log(z)/xixi(iens)
                        if ( z < 65 ) then
                            txtx(iens) = exp(z)/(1-pthreshold/100)
                        else
                            txtx(iens) = 1e20
                        end if
                    else
                        txtx(iens) = 1e20
                    end if
                end if
            else
                n = 0
                do i=1,ntot
                    if ( data(i) > xyear-thens ) n = n + 1
                end do
                txtx(iens) = real(ntot+1)/real(n)
            end if
            ! convert to years
            txtx(iens) = txtx(iens)/(j2-j1+1)
            if ( lwrite ) print *,'return time ',iens,txtx(iens)
        end if
    end do
    if ( lchangesign ) then
        t = -t
        tt = -tt
    end if
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

!   output

    if ( .not. lprint .and. .not. lwrite ) return
    if ( lweb ) then
        print '(a)','# <tr><td colspan="3">Fitted to GPD '// &
            'distribution H(x) = 1 - (1+&xi;(x-&mu;)/&sigma;)^(-1/&xi;)</td></tr>'
        if ( lchangesign ) then
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td>&mu;:'// &
                '</td><td>',-threshold,'</td><td>(',100-pthreshold,'%)</td></tr>'
        else
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td>&mu;:'// &
                '</td><td>',threshold,'</td><td>(',pthreshold,'%)</td></tr>'
        end if
        print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td>&sigma;:'// &
            '</td><td>',b,'</td><td>',b25,'...',b975,'</td></tr>'
        print '(a,f16.3,a,f16.3,a,f16.3,a)' &
            ,'# <tr><td>&xi;:</td><td>',xi,'</td><td>',xi25,'...' ,xi975,'</td></tr>'
    else
        print '(a,i5,a)','# Fitted to GPD distribution in ',iter,' iterations'
        print '(a)','# H(x) = 1-(1+xi*(x-a)/b)**(-1/xi) with'
        print '(a,f16.3,a,f16.3,a,f16.3)','# a = ',threshold,' (',pthreshold,'%)'
        print '(a,f16.3,a,f16.3,a,f16.3)','# b = ',b,' \\pm ',b975-b25
        print '(a,f16.3,a,f16.3,a,f16.3)','# xi= ',xi,' \\pm ',xi975-xi25
    end if
    call printreturnvalue(ntype,t,t25,t975,lweb)
    if ( lchangesign .and. xyear < 1e33 ) xyear = -xyear
    call printreturntime(year,xyear,tx,tx25,tx975,lweb)
    if ( lchangesign .and. xyear < 1e33 ) xyear = -xyear
    call plotreturnvalue(ntype,t25,t975,j2-j1+1)
end subroutine fitgpd

subroutine fit1gpd(b,xi,iter)
    use AmoebaToGSL
    implicit none
    integer :: iter
    real :: b,xi
    integer :: i
    real :: q(2),p(3,2),y(3),tol
    logical :: lok

    integer :: ncur
    integer,parameter :: nmax=100000
    real :: data(nmax),restrain,cthreshold
    logical :: llwrite,llchangesign
    common /fitdata1/ data
    common /fitdata2/ restrain,ncur,llwrite,llchangesign,cthreshold
    character cassume*5
    common /fitdata4/ cassume
    real,external :: llgpd

!   fit, using Numerical Recipes routines

10  continue
    q(1) = b
    q(2) = xi
    p(1,1) = q(1) *0.9
    p(1,2) = q(2) *0.9
    p(2,1) = p(1,1) *1.2
    p(2,2) = p(1,2)
    p(3,1) = p(1,1)
    p(3,2) = p(1,2) *1.2 + 0.1
    lok = .false. 
    do i=1,3
        q(1) = p(i,1)
        q(2) = p(i,2)
        y(i) = llgpd(q)
        if ( y(i) < 1e33 ) lok = .true. 
    end do
    if ( .not. lok ) then
        if ( xi /= 0 ) then
            xi = xi + 0.1*xi/abs(xi)
        else if ( cassume == 'scale' .and. llchangesign ) then
            xi = -0.1
        else
            xi = 0.1
        end if
        if ( abs(xi) < 2 ) then
            goto 10
        else
            write(0,*) 'fit1gpd: error: cannot find initial values'
            b = 3e33
            xi = 3e33
            return
        end if
    end if
    tol = 1e-4
    call amoeba(p,y,3,2,2,tol,llgpd,iter)
!   maybe add restart later
    b = p(1,1)
    xi = p(1,2)
    if ( abs(xi) > 10 ) then
        write(0,*) 'fit1gpd: error: shape parameter xi = ',xi
    end if
    if ( abs(b) > 1e15 ) then
        write(0,*) 'fit1gpd: error: position parameter b = ',b
    end if
end subroutine fit1gpd

real function llgpd(p)

!   computes the log-likelihood function for a GPD distribution
!   with parameters a,b=p(1),p(2) and data in common.

    implicit none

    real :: p(2)

    integer :: i
    real :: b,xi,s,z,llold
    integer,save :: init
    data init /0/

    integer :: nmax,ncur
    parameter(nmax=100000)
    real :: data(nmax),restrain,cthreshold
    logical :: llwrite,llchangesign
    common /fitdata1/ data
    common /fitdata2/ restrain,ncur,llwrite,llchangesign,cthreshold
    character cassume*5
    common /fitdata4/ cassume

    llgpd = 0
    b = p(1)
    xi = p(2)
!!!print *,'llgpd:b,xi = ',b,xi
    if ( b <= 1e-15 ) then
        llgpd = 3e33
        goto 999
    end if
    if ( abs(xi) > 10 ) then
        llgpd = 3e33
        goto 999
    end if
    if ( restrain < 0 ) then
        write(0,*) 'llgpd: restrain<0 ',restrain
        call exit(-1)
    end if
    if ( llchangesign .and. cassume == 'scale' ) then
        if ( init == 0 ) then
            init = 1
            write(0,*) 'Enforcing a hard lower bound of zero.'
        end if
        if ( xi >= 0 ) then
            ! scaling implies (for climate) that the distribution cannot cross zero,
            ! so xi must be < 0
            llgpd = 3e33
            goto 999
        end if
        if ( b > abs(cthreshold*xi) ) then
            ! scaling implies (for climate) that the distribution cannot cross zero,
            ! so the upper limit must be <0 (for flipped signs)
            llgpd = 3e33
            goto 999
        end if
    end if
                
    do i=1,ncur
        z = data(i)
        if ( z < 0 ) then
            write(0,*) 'llgpd: z<0 ',z,i,ncur
            call exit(-1)
        end if
        if ( 1+xi*z/b <= 0 ) then
            llgpd = 3e33
            goto 999
        end if
        llold = llgpd
        if ( abs(xi) < 1e-4 ) then
            llgpd = llgpd - z/b + (z/b)**2*xi/2
        !!!                print *,i,z, - z/b + (z/b)**2*xi/2 - log(b)
        else
            llgpd = llgpd - (1+1/xi)*log(1+xi*z/b)
        !!!                print *,i,z, - (1+1/xi)*log(1+xi*z/b) - log(b)
        end if
        if ( llwrite ) print *,i,z,b,llgpd-llold
    end do
    if ( restrain == 0 ) then
        llgpd = llgpd - ncur*log(b)
    else
        !   the -2 gives visually much better results;
        !   I am unsure of the mathematical derivation
        llgpd = llgpd - (ncur-2)*log(b)
        !   preconditioning on xi with gaussian of width restrain/2
        !   around 0
        llgpd = llgpd - (xi/(restrain/2))**2/2
    end if
    !       minimum, not maximum
    llgpd = -llgpd
999 continue
!!!        print '(a,i3,2f6.2,f10.4)','ncur,b,xi,llgpd = ',ncur,p(1),p(2)
!!!     +       ,llgpd

end function llgpd
