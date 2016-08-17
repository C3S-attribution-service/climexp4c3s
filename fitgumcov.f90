!  #[ fitgumcov:
subroutine fitgumcov(yrseries,yrcovariate,npernew,fyr,lyr        &
 &       ,mens1,mens,crosscorr,a3,b3,alpha3,beta3,j1,j2          &
 &       ,lweb,ntype,lchangesign,yr1a,yr2a,xyear,idmax,cov1,cov2 &
 &       ,offset,t3,tx3,assume,confidenceinterval,ndecor         &
 &       ,lboot,lprint,dump,plot,lwrite)
!
!   a fit a Gumbel distribution to the data, which is already assumed to be block max
!   input:
!   xx(ntot) data
!   j1,j2    use days/months/... j1 to j2
!   output
!   a,b,alpha,beta  parameters of fit
!   t(10,3)    return values for 10, 20, 50, ..., 10000 years for cov=cov1,cov2 and the difference
!   t25,t975   (100-confidenceinterval)/2%, (100+confidenceinterval)/2% quantiles of these return values
!   tx(3)      return time of the value of year (xyear) in the context of the other values and difference
!   tx25,tx975 (100-confidenceinterval)/2%, (100-confidenceinterval)/2% quantiles of these return times and their differences
!
    implicit none
!
    integer npernew,fyr,lyr,mens1,mens,nmc,ntot,j1,j2,ntype,yr1a,yr2a,ndecor
    real yrseries(npernew,fyr:lyr,0:mens),yrcovariate(npernew,fyr:lyr,0:mens), &
 &       crosscorr(0:mens,0:mens),a3(3),b3(3),alpha3(3),beta3(3),xyear,        &
 &       cov1,cov2,offset,t3(3,10,3),tx3(3,3),confidenceinterval
    character assume*(*),idmax*(*)
    logical lweb,lchangesign,lboot,lprint,dump,plot,lwrite
!
    integer i,j,nx,iter,iens,nfit,year
    integer,allocatable :: yrs(:)
    real,allocatable :: aa(:),bb(:),xixi(:),tt(:,:,:),         &
 &       txtx(:,:),alphaalpha(:),betabeta(:),aacov(:,:)
    real x,a,b,xi,alpha,beta,t5(10,3),t1(10,3),db,f            &
 &       ,threshold,thens,z,ll,ll1                             &
 &       ,a25,a975,b25,b975,alpha25,alpha975,beta25,beta975    &
 &       ,ranf,mean,sd,dalpha,dbeta,mindata,minindx,pmindata   &
 &       ,snorm,s,frac,t(10,3),t25(10,3),t975(10,3)            &
 &       ,tx(3),tx25(3),tx975(3),ttt(10,3),txtxtx(3),xi3(3)    &
 &       ,acov(3,2),cmin,cmax,plo,phi,scross,sdecor
    real adev,var,skew,curt,aaa,bbb,siga,chi2,q
    real,allocatable :: xx(:,:),yy(:),ys(:),zz(:),sig(:)
    character lgt*4,method*3
!
    integer nmax,ncur
    parameter(nmax=100000)
    real data(2,nmax),restrain
    logical llwrite,llchangesign
    common /fitdata3/ data
    common /fitdata2/ restrain,ncur,llwrite,llchangesign
    character cassume*5
    common /fitdata4/ cassume
!
    real llgumbelcov,gevcovreturnlevel,gevcovreturnyear
    external llgumbelcov,gevcovreturnlevel,gevcovreturnyear
!
!       estimate number of bootstrap samples needed, demand at least 25 above threshold
!
    nmc = max(1000,nint(25*2/(1-confidenceinterval/100)))
    allocate(yrs(0:nmax))
    allocate(xx(2,nmax),aa(nmc),bb(nmc),xixi(nmc),tt(nmc,10,3),  &
 &       txtx(nmc,3),alphaalpha(nmc),betabeta(nmc),aacov(nmc,2))
    year = yr2a

    if ( lwrite ) print *,'fitgumcov: calling fill_linear_array'
    call fill_linear_array(yrseries,yrcovariate,npernew,j1,j2,   &
 &       fyr,lyr,mens1,mens,xx,yrs,nmax,ntot,lwrite)
    if ( lprint .and. lweb ) then
        print '(a,i9,a)','# <tr><td>N:</td><td>&nbsp;</td><td>', &
 &           ntot,'</td><td>&nbsp;</td></tr>'
    end if


    if ( lwrite ) then
        print *,'fitgumcov: input:'
        print *,'j1,j2      = ',j1,j2
        print *,'year,xyear = ',year,xyear
        print *,'cov1,cov2,offset ',cov1,cov2,offset
        if ( .false. ) then
            do i=1,ntot
                print *,i,xx(:,i)
            end do
        end if
    end if
!
!   compute first-guess parameters
!
    allocate(yy(ntot))
    allocate(ys(ntot))
    allocate(zz(ntot))
    allocate(sig(ntot))
    cmin = 3e33
    cmax = -3e33
    do i=1,ntot
        yy(i) = xx(1,i)
        zz(i) = xx(2,i)
        cmin = min(cmin,xx(2,i))
        cmax = max(cmax,xx(2,i))
    end do
    call write_obscov(xx,yrs,ntot,-3e33,cov2,xyear,year,offset,lchangesign)

    sig = 0
    call moment(yy,ntot,mean,adev,sd,var,skew,curt)
    call fit(zz,yy,ntot,sig,0,aaa,alpha,siga,dalpha,chi2,q)
    if ( lwrite ) then
        print *,'fitgevcov: computed initialisation values:'
        print *,'mean,sd,alpha,dalpha = ',mean,sd,alpha,dalpha
    end if
!
!   ill-defined case
!
    if ( sd == 0 ) then
        a3 = 3e33
        b3 = 3e33
        alpha3 = 3e33
        beta3 = 3e33
        t3 = 3e33
        tx3 = 3e33
        return
    end if
!
!   copy to common for routine llgumbel
!
    ncur = ntot
    do i=1,ncur
        data(:,i) = xx(:,i)
    end do
    llwrite = lwrite
    cassume = assume
!
    b = sd*sqrt(6.)/(4*atan(1.))
    a = mean - 0.57721*b
    if ( assume == 'shift' .or. assume == 'scale' ) then
        beta = 3e33
        call fit1gumcov(a,b,alpha,dalpha,iter)
    else if ( assume == 'both' ) then
        beta = alpha
        dbeta = dalpha
        call fit2gumcov(a,b,alpha,beta,dalpha,dbeta,iter)
    else
        write(0,*) 'fitgumcov: error: unknown value for assume ', &
 &           assume
    end if
    xi = 0
    call getreturnlevels(a,b,xi,alpha,beta,cov1,cov2,             &
 &       gevcovreturnlevel,j1,j2,assume,t)
    if ( xyear < 1e33 ) then
        call getreturnyears(a,b,xi,alpha,beta,xyear,cov1,cov2,    &
 &           gevcovreturnyear,j1,j2,tx,lchangesign,lwrite)
    end if
    call getabfromcov(a,b,alpha,beta,cov1,aaa,bbb)
    acov(1,1) = aaa
    call getabfromcov(a,b,alpha,beta,cov2,aaa,bbb)
    acov(1,2) = aaa
    call write_threshold(cmin,cmax,a,b,alpha,beta,offset,         &
 &       lchangesign)
!
!   bootstrap to find error estimates
!
    if ( .not.lboot ) then
        if ( lchangesign ) then
            a = -a
            t = -t
            alpha = -alpha
            beta = -beta
        end if
        a3(1) = a
        a3(2:3) = 3e33
        b3(1) = b
        b3(2:3) = 3e33
        xi3 = 0
        alpha3(1) = alpha
        alpha3(2:3) = 3e33
        beta3(1) = beta
        beta3(2:3) = 3e33
        t3(1,:,:) = t(:,:)
        t3(2:3,:,:) = 3e33
        tx3(1,:) = tx(:)
        tx3(2:3,:) = 3e33            
        return
    end if
    if ( lprint .and. .not.lweb ) print '(a,i6,a)','# Doing a ',nmc  &
 &        ,'-member bootstrap to obtain error estimates'
    scross = 0
    do iens=1,nmc
        call keepalive1('Computing bootstrap sample ',iens,nmc)
        method = 'new'
        if ( method == 'old' ) then
            do i=1,ncur
                call random_number(ranf)
                j = 1+int(ntot*ranf)
                if ( j < 1 .or. j > ncur ) then
                    write(0,*) 'fitgumcov: error: j = ',j
                    call abort
                end if
                data(:,i) = xx(:,j)
            end do
        else
            call sample_bootstrap(yrseries,yrcovariate,             &
 &               npernew,j1,j2,fyr,lyr,mens1,mens,crosscorr,        &
 &               ndecor,data,nmax,ntot,sdecor,lwrite)
            scross = scross + sdecor
        end if
        aa(iens) = a
        bb(iens) = b
        alphaalpha(iens) = alpha
        if ( assume == 'shift' .or. assume == 'scale' ) then
            betabeta(iens) = 3e33
            call fit1gumcov(aa(iens),bb(iens),                       &
 &               alphaalpha(iens),dalpha,iter)
        else if ( assume == 'both' ) then
            betabeta(iens) = beta
            call fit2gumcov(aa(iens),bb(iens),                       &
 &               alphaalpha(iens),betabeta(iens),dalpha,dbeta,iter)
        else
            write(0,*) 'fitgumcov: cannot handle assume = ',assume
            call abort
        end if
        call getabfromcov(aa(iens),bb(iens),                         &
 &           alphaalpha(iens),betabeta(iens),cov1,aaa,bbb)
        aacov(iens,1) = aaa
        call getabfromcov(aa(iens),bb(iens),                         &
 &           alphaalpha(iens),betabeta(iens),cov2,aaa,bbb)
        aacov(iens,2) = aaa
        xi = 0
        call getreturnlevels(aa(iens),bb(iens),xi,alphaalpha(iens),  &
 &           betabeta(iens),cov1,cov2,gevcovreturnlevel,j1,j2,assume,ttt)
        do i=1,10
            do j=1,3
                tt(iens,i,j) = ttt(i,j)
            end do
        end do
        if ( xyear < 1e33 ) then
            call getreturnyears(aa(iens),bb(iens),xi,                 &
 &               alphaalpha(iens),betabeta(iens),xyear,cov1,cov2,     &
 &               gevcovreturnyear,j1,j2,txtxtx,lchangesign,lwrite)
            do j=1,3
                txtx(iens,j) = txtxtx(j)
            end do
        end if
    end do
    if ( mens > mens1 ) call print_spatial_scale(scross/nmc)
    iens = nmc
    if ( lchangesign ) then
        a = -a
        acov = -acov
        aa = -aa
        aacov = -aacov
        alpha = -alpha
        alphaalpha = -alphaalpha
        if ( assume == 'both' ) then
            beta = -beta
            betabeta = -betabeta
        end if
        t = -t
        tt = -tt
    end if
    plo = (100-confidenceinterval)/2
    phi = (100+confidenceinterval)/2
    call getcut( a25,plo,nmc,aa)
    call getcut(a975,phi,nmc,aa)
    call getcut( b25,plo,nmc,bb)
    call getcut(b975,phi,nmc,bb)
    call getcut( alpha25,plo,nmc,alphaalpha)
    call getcut(alpha975,phi,nmc,alphaalpha)
    if ( assume == 'both' ) then
        call getcut( beta25,plo,nmc,betabeta)
        call getcut(beta975,phi,nmc,betabeta)
    else
        beta25 = 3e33
        beta975 = 3e33
    end if
    do i=1,10
        do j=1,3
            if ( lchangesign ) then
                lgt = '&lt;'
                call getcut(t5(i,j),5.,nmc,tt(1,i,j))
                call getcut(t1(i,j),1.,nmc,tt(1,i,j))
            else
                lgt = '&gt;'
                call getcut(t5(i,j),95.,nmc,tt(1,i,j))
                call getcut(t1(i,j),99.,nmc,tt(1,i,j))
            end if
            call getcut( t25(i,j),plo,nmc,tt(1,i,j))
            call getcut(t975(i,j),phi,nmc,tt(1,i,j))
        end do
    end do
    do j=1,3
        if ( xyear < 1e33 ) then
            call getcut( tx25(j),plo,nmc,txtx(1,j))
            call getcut(tx975(j),phi,nmc,txtx(1,j))
            if ( lchangesign ) xyear = -xyear
        end if
    end do
    call getcut(acov(2,1),plo,iens,aacov(1,1))
    call getcut(acov(3,1),phi,iens,aacov(1,1))
    call getcut(acov(2,2),plo,iens,aacov(1,2))
    call getcut(acov(3,2),phi,iens,aacov(1,2))
    call write_dthreshold(cov1,cov2,acov,offset,lchangesign)
!
!   output
!
    if ( .not.lprint ) then
        call copyab3etc(a3,b3,xi3,alpha3,beta3,t3,tx3,                &
 &           a,a25,a975,b,b975,xi,xi,xi,alpha,alpha25,alpha975,       &
 &           beta,beta25,beta975,t,t25,t975,tx,tx25,tx975)
        if ( .not.lwrite ) return
    end if
    if ( lweb ) then
        print '(a)','# <tr><td colspan="4">Fitted to Gumbel '//       &
 &           'distribution P(x) = exp(-exp(-(x-&mu;'')/&sigma;''))'   &
 &           //'</td></tr>'
        call printab(lweb)
        print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'//  &
 &           '&mu;:</td><td>',a,'</td><td>',a25,'...',a975,           &
 &           '</td></tr>'
        print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'//  &
 &           '&sigma;:</td><td>',b,'</td><td>',b25,'...',b975,        &
 &           '</td></tr>'
        print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'//  &
 &           '&alpha;:</td><td>',alpha,'</td><td>',alpha25,'...',     &
 &           alpha975,'</td></tr>'
        if ( assume == 'both' ) then
            print '(a,f16.3,a,f16.3,a,f16.3,a)',                       &
 &               '# <tr><td colspan=2>&beta;:</td><td>',beta,          &
 &               '</td><td>',beta25,'...',beta975,'</td></tr>'
        end if
    else
        print '(a,i5,a)','# Fitted to Gumbel distribution in ',iter,' iterations'
        print '(a)','# P(x) = exp(-exp(-(x-a'')/b''))'
        call printab(lweb)
        print '(a,f16.3,a,f16.3,a,f16.3)','# a = ',a,' \\pm ',(a975-a25)/2
        print '(a,f16.3,a,f16.3,a,f16.3)','# b = ',b,' \\pm ',(b975-b25)/2
        print '(a,f16.3,a,f16.3,a,f16.3)','# alpha ',alpha,' \\pm ',  &
 &           (alpha975-alpha25)/2
        if ( assume == 'both' ) then
            print '(a,f16.3,a,f16.3,a,f16.3)','# beta  ',beta,        &
 &               ' \\pm ',(beta975-beta25)/2
        end if
    end if
    call printcovreturnvalue(ntype,t,t25,t975,yr1a,yr2a,lweb,plot,assume)
    call printcovreturntime(year,xyear,idmax,tx,tx25,tx975,yr1a,yr2a  &
 &       ,lweb,plot)
    call printcovpvalue(txtx,nmc,nmc,lweb)

    if ( dump ) then
        call plot_tx_cdfs(txtx,nmc,nmc,ntype,j1,j2)
    end if
    if ( plot ) write(11,'(3g20.4,a)') alpha,alpha25,alpha975,' alpha'

    ! no cuts
    mindata = -2e33
    minindx = -2e33
    pmindata = -1
    snorm = 1
    frac = 1
    ! gumbel fit
    nfit = 4
    if ( lchangesign ) then
        a = -a
        alpha = -alpha
    end if

    ! compute distribution at past year and plot it
    call adjustyy(ntot,xx,assume,a,b,alpha,beta,cov1,      &
 &       yy,zz,aaa,bbb,lchangesign,lwrite)
    ys(1:ntot) = yy(1:ntot)
    print '(a,i5)','# distribution in year ',yr1a
    call plotreturnvalue(ntype,t25(1,1),t975(1,1),j2-j1+1)
    call plot_ordered_points(yy,ys,yrs,ntot,ntype,nfit,    &
 &       frac,aaa,bbb,xi,j1,j2,minindx,mindata,pmindata,   &
 &       year,xyear,snorm,lchangesign,lwrite,.false.)

    ! compute distribution at past year and plot it
    call adjustyy(ntot,xx,assume,a,b,alpha,beta,cov2,      &
 &       yy,zz,aaa,bbb,lchangesign,lwrite)
    ys(1:ntot) = yy(1:ntot)
    print '(a)'
    print '(a)'
    print '(a,i5)','# distribution in year ',yr2a
    call plotreturnvalue(ntype,t25(1,2),t975(1,2),j2-j1+1)
    call plot_ordered_points(yy,ys,yrs,ntot,ntype,nfit,    &
 &       frac,aaa,bbb,xi,j1,j2,minindx,mindata,pmindata,   &
 &       year,xyear,snorm,lchangesign,lwrite,.true.)

end subroutine
!  #] fitgumcov:
!  #[ fit1gumcov:
subroutine fit1gumcov(a,b,alpha,dalpha,iter)
    implicit none
    integer iter
    real a,b,alpha,dalpha
    integer i
    real q(4),p(4,3),y(4),tol
    real llgumbelcov
    external llgumbelcov
!
!   fit, using Numerical Recipes routines
!
    q(1) = a
    q(2) = b
    q(3) = alpha
    q(4) = 3e33
    p(1,1) = q(1) *0.9
    p(1,2) = q(2) *0.9
    p(1,3) = q(3) - dalpha
    p(2,1) = p(1,1) *1.2
    p(2,2) = p(1,2)
    p(2,3) = p(1,3)
    p(3,1) = p(1,1)
    p(3,2) = p(1,2) *1.2
    p(3,3) = p(1,3)
    p(4,1) = p(1,1)
    p(4,2) = p(1,2)
    p(4,3) = p(1,3) + 2*dalpha
    do i=1,4
        q(1) = p(i,1)
        q(2) = p(i,2)
        q(3) = p(i,3)
        y(i) = llgumbelcov(q)
    end do
    tol = 1e-4
    call amoeba(p,y,4,3,3,tol,llgumbelcov,iter)
!   maybe add restart later
    a = p(1,1)
    b = p(1,2)
    alpha = p(1,3)
end subroutine
!  #] fit1gumcov:
!  #[ fit2gumcov:
subroutine fit2gumcov(a,b,alpha,beta,dalpha,dbeta,iter)
    implicit none
    integer iter
    real a,b,alpha,beta,dalpha,dbeta
    integer i
    real q(4),p(5,4),y(5),tol
    real llgumbelcov
    external llgumbelcov
!
!   fit, using Numerical Recipes routines
!
    q(1) = a
    q(2) = b
    q(3) = alpha
    q(4) = beta
    q(4) = beta
    p(1,1) = q(1) *0.9
    p(1,2) = q(2) *0.9
    p(1,3) = q(3) - dalpha
    p(1,4) = q(4) - dbeta
    p(2,1) = p(1,1) *1.2
    p(2,2) = p(1,2)
    p(2,3) = p(1,3)
    p(2,4) = p(1,4)
    p(3,1) = p(1,1)
    p(3,2) = p(1,2) *1.2
    p(3,3) = p(1,3)
    p(3,4) = p(1,4)
    p(4,1) = p(1,1)
    p(4,2) = p(1,2)
    p(4,3) = p(1,3) + 2*dalpha
    p(4,4) = p(1,4)
    p(5,1) = p(1,1)
    p(5,2) = p(1,2)
    p(5,3) = p(1,3)
    p(5,4) = p(1,4) + 2*dbeta
    do i=1,5
        q(1) = p(i,1)
        q(2) = p(i,2)
        q(3) = p(i,3)
        q(4) = p(i,4)
        y(i) = llgumbelcov(q)
    end do
    tol = 1e-4
    call amoeba(p,y,5,4,4,tol,llgumbelcov,iter)
!   maybe add restart later
    a = p(1,1)
    b = p(1,2)
    alpha = p(1,3)
    beta = p(1,4)
end subroutine
!  #] fit2gumcov:
!  #[ llgumbelcov:
real function llgumbelcov(p)
!
!   computes the log-likelihood function for a Gumbel distribution
!   with parameters a,b,alpha=p(1),p(2),p(3) and data in common.
!
    implicit none
!   
    real p(4)
!
    integer i
    real z,s,aa,bb
!
    integer nmax,ncur
    parameter(nmax=100000)
    real data(2,nmax),restrain
    logical llwrite,llchangesign
    common /fitdata3/ data
    common /fitdata2/ restrain,ncur,llwrite,llchangesign
!
    llgumbelcov = 0
    do i=1,ncur
        call getabfromcov(p(1),p(2),p(3),p(4),data(2,i),aa,bb)
        if ( abs(p(2)) < 1e-30 ) then
            llgumbelcov = 3e33
            goto 999
        end if
        z = (data(1,i) - aa)/bb
        if ( z < -log(3e33) ) then
            llgumbelcov = 3e33
            goto 999
        else
            llgumbelcov = llgumbelcov - exp(-z) - z - log(abs(bb))
        end if
    end do
!   normalization is not 1 in case of cut-offs
    call gumbcovnorm(p(1),p(2),p(3),p(4),s)
    if ( s < 1e33 ) then
        llgumbelcov = llgumbelcov - ncur*log(s)
    else
        llgumbelcov = 3e33
        goto 999
    end if
!   minimum, not maximum
    llgumbelcov = -llgumbelcov
    if ( llwrite ) print *,'a,b,alpha,llgumcov = ',    &
 &       p(1),p(2),p(3),llgumbelcov
!
999 continue
end function
!  #] llgumbelcov:
!  #[ gumbcovnorm:
subroutine gumbcovnorm(a,b,alpha,beta,s)
    implicit none
	include "getopts.inc"
    real a,b,alpha,beta,s
    real z1,z2

    if ( minindx > -1e33 .or. maxindx < 1e33 ) then
        write(0,*) 'gumbcovnorm: boundaries not yet available for fit of gumbel(t)'
        call abort
    else
        s = 1
    end if
!!!    print *,'gumbcovnorm: norm = ',a,b,alpha,s
end subroutine
!  #] gumbcovnorm:
