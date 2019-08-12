subroutine fitgumcov(yrseries,yrcovariate,npernew,fyr,lyr             &
        ,mens1,mens,crosscorr,a3,b3,alpha3,beta3,j1,j2,nens1,nens2   &
        ,lweb,ntype,lchangesign,yr1a,yr2a,yr2b,xyear,idmax,cov1,cov2 &
        ,cov3,offset,t3,tx3,assume,confidenceinterval,ndecor              &
        ,lboot,lprint,dump,plot,lwrite)
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
    integer npernew,fyr,lyr,mens1,mens,nmc,ntot,j1,j2,nens1,nens2, &
        ntype,yr1a,yr2a,yr2b,ndecor
    real yrseries(npernew,fyr:lyr,0:mens),yrcovariate(npernew,fyr:lyr,0:mens), &
        crosscorr(0:mens,0:mens),a3(3),b3(3),alpha3(3),beta3(3),xyear,        &
        cov1,cov2,cov3,offset,t3(3,10,3),tx3(3,3),confidenceinterval
    character assume*(*),idmax*(*)
    logical lweb,lchangesign,lboot,lprint,dump,plot,lwrite
!
    integer i,j,nx,iter,iens,nfit,year,nj
    integer,allocatable :: yrs(:)
    real,allocatable :: aa(:),bb(:),baba(:),xixi(:),           &
        tt(:,:,:),txtx(:,:),alphaalpha(:),betabeta(:),        &
        aacov(:,:)
    real x,a,b,ba,xi,alpha,beta,t5(10,4),t1(10,4),db,f         &
        ,threshold,thens,z,ll,ll1                             &
        ,a25,a975,b25,b975,alpha25,alpha975,beta25,beta975    &
        ,aa25,aa975,bb25,bb975,ba25,ba975                     &
        ,ranf,mean,sd,dalpha,dbeta,mindata,minindx,pmindata   &
        ,snorm,s,frac,t(10,4),t25(10,4),t975(10,4)            &
        ,tx(4),tx25(4),tx975(4),ttt(10,4),txtxtx(4),xi3(4)    &
        ,acov(3,3),cmin,cmax,plo,phi,scross,sdecor
    real adev,var,skew,curt,aaa,bbb,siga,chi2,q
    real,allocatable :: xx(:,:),yy(:),ys(:),zz(:),sig(:)
    logical lnone,last
    character lgt*4,method*3
!
    integer nmax,ncur
    parameter(nmax=1000000)
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
!   estimate number of bootstrap samples needed, demand at least 25 above threshold
!
    nmc = max(1000,nint(25*2/(1-confidenceinterval/100)))
    allocate(yrs(0:nmax))
    allocate(xx(2,nmax),aa(nmc),bb(nmc),baba(nmc),xixi(nmc),tt(nmc,10,4),  &
        txtx(nmc,4),alphaalpha(nmc),betabeta(nmc),aacov(nmc,3))
    year = yr2a
    if ( cov1 == 0 .and. cov2 == 0 ) then
        lnone = .true.
    else
        lnone = .false.
    end if
    restrain = 0 ! no shape parameter

    if ( lwrite ) print *,'fitgumcov: calling fill_linear_array'
    call fill_linear_array(yrseries,yrcovariate,npernew,j1,j2,   &
        fyr,lyr,mens1,mens,xx,yrs,nmax,ntot,lwrite)
    if ( lprint .and. lweb ) then
        print '(a,i9,a)','# <tr><td>N:</td><td>&nbsp;</td><td>', &
            ntot,'</td><td>&nbsp;</td></tr>'
    end if
    if ( ntot < 5 ) then
        print '(a)','</table>'
        print '(a)','There are not enough points to fit a Gumbel distribution. Please use a (much) longer time series.'
        print '(a)','<p>'
        call exit(-1)
    end if

    if ( lwrite ) then
        print *,'fitgumcov: input:'
        print *,'j1,j2      = ',j1,j2
        print *,'year,xyear = ',year,xyear
        print *,'cov1,cov2,offset ',cov1,cov2,offset
        if ( cov3 < 1e33 ) then
            print *,'cov3           = ',cov3
        end if
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
    if ( dump ) call write_obscov(xx,yrs,ntot,-3e33,cov2,xyear,year,offset,lchangesign)

    sig = 0
    call moment(yy,ntot,mean,adev,sd,var,skew,curt)
    call fit(zz,yy,ntot,sig,0,aaa,alpha,siga,dalpha,chi2,q)
    if ( lwrite ) then
        print *,'fitgumcov: computed initialisation values:'
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
    if ( lnone ) then
        alpha = 3e33
        beta = 3e33
        call fit0gumcov(a,b,iter)
    else if ( assume == 'shift' .or. assume == 'scale' ) then
        beta = 3e33
        call fit1gumcov(a,b,alpha,dalpha,iter)
    else if ( assume == 'both' ) then
        beta = alpha
        dbeta = dalpha
        call fit2gumcov(a,b,alpha,beta,dalpha,dbeta,iter)
    else
        write(0,*) 'fitgumcov: error: unknown value for assume ',assume
    end if
    if ( assume == 'scale' ) ba = b/a
    xi = 0
    call getreturnlevels(a,b,xi,alpha,beta,cov1,cov2,cov3,gevcovreturnlevel,j1,j2,assume,t)
    if ( xyear < 1e33 ) then
        call getreturnyears(a,b,xi,alpha,beta,xyear,cov1,cov2,cov3, &
            gevcovreturnyear,j1,j2,tx,lchangesign,lwrite)
    end if
    call getabfromcov(a,b,alpha,beta,cov1,aaa,bbb)
    acov(1,1) = aaa
    call getabfromcov(a,b,alpha,beta,cov2,aaa,bbb)
    acov(1,2) = aaa
    if ( cov3 < 1e33 ) then
        call getabfromcov(a,b,alpha,beta,cov3,aaa,bbb)
        acov(1,3) = aaa
    end if
    if ( dump ) call write_threshold(cmin,cmax,a,b,xi,alpha,beta,offset,lchangesign,gevcovreturnlevel)
!
!   bootstrap to find error estimates
!
    if ( .not.lboot ) then
        if ( lchangesign ) then
            a = -a
            b = -b
            t = -t
            if ( .not. lnone ) then
                alpha = -alpha
                if ( cassume == 'both' ) then
                    beta = -beta
                end if
            end if
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
        t3(1,:,1:3) = t(:,1:3)
        t3(2:3,:,:) = 3e33
        tx3(1,:) = tx(1:3)
        tx3(2:3,:) = 3e33            
        return
    end if
    if ( lprint .and. .not.lweb ) print '(a,i6,a)','# Doing a ',nmc  &
         ,'-member bootstrap to obtain error estimates'
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
                    call exit(-1)
                end if
                data(:,i) = xx(:,j)
            end do
        else
            call sample_bootstrap(yrseries,yrcovariate,             &
                npernew,j1,j2,fyr,lyr,nens1,nens2,crosscorr,       &
                ndecor,data,nmax,ntot,sdecor,lwrite)
            scross = scross + sdecor
        end if
        aa(iens) = a
        bb(iens) = b
        alphaalpha(iens) = alpha
        if ( lnone ) then
            alphaalpha(iens) = 3e33
            betabeta(iens) = 3e33
            call fit0gumcov(aa(iens),bb(iens),iter)
        elseif ( assume == 'shift' .or. assume == 'scale' ) then
            betabeta(iens) = 3e33
            call fit1gumcov(aa(iens),bb(iens),alphaalpha(iens),dalpha,iter)
        else if ( assume == 'both' ) then
            betabeta(iens) = beta
            call fit2gumcov(aa(iens),bb(iens),                       &
                alphaalpha(iens),betabeta(iens),dalpha,dbeta,iter)
        else
            write(0,*) 'fitgumcov: cannot handle assume = ',assume
            call exit(-1)
        end if
        if ( assume == 'scale' ) baba(iens) = bb(iens)/aa(iens)
        call getabfromcov(aa(iens),bb(iens),                         &
            alphaalpha(iens),betabeta(iens),cov1,aaa,bbb)
        aacov(iens,1) = aaa
        call getabfromcov(aa(iens),bb(iens),                         &
            alphaalpha(iens),betabeta(iens),cov2,aaa,bbb)
        aacov(iens,2) = aaa
        if ( cov3 < 1e33 ) then
            call getabfromcov(aa(iens),bb(iens),                     &
                alphaalpha(iens),betabeta(iens),cov3,aaa,bbb)
            aacov(iens,3) = aaa
        end if
        xi = 0
        call getreturnlevels(aa(iens),bb(iens),xi,alphaalpha(iens),  &
            betabeta(iens),cov1,cov2,cov3,gevcovreturnlevel,j1,j2,assume,ttt)
        do i=1,10
            do j=1,4
                tt(iens,i,j) = ttt(i,j)
            end do
        end do
        if ( xyear < 1e33 ) then
            call getreturnyears(aa(iens),bb(iens),xi,                 &
                alphaalpha(iens),betabeta(iens),xyear,cov1,cov2,cov3, &
                gevcovreturnyear,j1,j2,txtxtx,lchangesign,lwrite)
            do j=1,4
                txtx(iens,j) = txtxtx(j)
            end do
        end if
    end do
    if ( mens > mens1 ) call print_spatial_scale(scross/(nmc*ndecor),ntot/real(mens-mens1+1)/real(ndecor))
    iens = nmc
    if ( lchangesign ) then
        a = -a
        acov = -acov
        aa = -aa
        aacov = -aacov
        b = -b
        bb = -bb
        if ( .not. lnone ) then
            alpha = -alpha
            alphaalpha = -alphaalpha
            if ( assume == 'both' ) then
                beta = -beta
                betabeta = -betabeta
            end if
        end if
        do j=1,10
            do i=1,4
                if ( assume == 'scale' .and. i == 3 ) cycle ! a relative change does not change sign...
                if ( t(j,i) < 1e30 ) then
                    t(j,i) = -t(j,i)
                end if
                do iens = 1,nmc
                    if ( tt(iens,j,i) < 1e33 ) then
                        tt(iens,j,i) = -tt(iens,j,i)
                    end if
                end do
            end do
        end do        
    end if
    plo = (100-confidenceinterval)/2
    phi = (100+confidenceinterval)/2
    call getcut( a25,plo,nmc,aa)
    call getcut(a975,phi,nmc,aa)
    call getcut( b25,plo,nmc,bb)
    call getcut(b975,phi,nmc,bb)
    if ( assume == 'scale' ) then
        call getcut( ba25,plo,iens,baba)
        call getcut(ba975,phi,iens,baba)
    end if
    if ( .not. lnone ) then
        call getcut( alpha25,plo,nmc,alphaalpha)
        call getcut(alpha975,phi,nmc,alphaalpha)
        if ( assume == 'both' ) then
            call getcut( beta25,plo,nmc,betabeta)
            call getcut(beta975,phi,nmc,betabeta)
        end if
    end if
    if ( cov3 < 1e33 ) then
        nj = 4
    else
        nj = 3
    end if
    do j=1,nj-1
        call getcut(acov(2,j),plo,iens,aacov(1,j))
        call getcut(acov(3,j),phi,iens,aacov(1,j))
    end do
    do i=1,10
        do j=1,nj
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
    do j=1,nj
        if ( xyear < 1e33 ) then
            call getcut( tx25(j),plo,nmc,txtx(1,j))
            call getcut(tx975(j),phi,nmc,txtx(1,j))
        end if
    end do
    if ( xyear < 1e33 ) then
        if ( lchangesign ) xyear = -xyear
    end if
!
!   output
!
    if ( .not.lprint ) then
        call copyab3etc(a3,b3,xi3,alpha3,beta3,t3,tx3,                &
            a,a25,a975,b,b975,xi,xi,xi,alpha,alpha25,alpha975,       &
            beta,beta25,beta975,t,t25,t975,tx,tx25,tx975)
        if ( .not.lwrite ) return
    end if
    if ( lweb ) then
        if ( lnone ) then
            print '(a)','# <tr><td colspan="4">Fitted to Gumbel '//            &
                'distribution P(x) = exp(-exp(-(x-&mu;)/&sigma;))</td></tr>'
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'//         &
                '&mu;:</td><td>',a,'</td><td>',a25,'...',a975,'</td></tr>'
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'//         &
                '&sigma;:</td><td>',b,'</td><td>',b25,'...',b975,'</td></tr>'
            if ( assume == 'scale' ) then
                print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'// &
                       '&sigma;/&mu;:</td><td>',ba,'</td><td>',ba25,'...',ba975,'</td></tr>'
            end if
        else
            print '(a)','# <tr><td colspan="4">Fitted to Gumbel '//       &
                'distribution P(x) = exp(-exp(-(x-&mu;'')/&sigma;''))'   &
                //'</td></tr>'
            call printab(restrain,lnone,lweb)
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'//  &
                '&mu;:</td><td>',a,'</td><td>',a25,'...',a975,'</td></tr>'
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'//  &
                '&sigma;:</td><td>',b,'</td><td>',b25,'...',b975,'</td></tr>'
            if ( assume == 'scale' ) then
                print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'// &
                       '&sigma;/&mu;:</td><td>',ba,'</td><td>',ba25,'...',ba975,'</td></tr>'
            end if
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'//  &
                '&alpha;:</td><td>',alpha,'</td><td>',alpha25,'...',     &
                alpha975,'</td></tr>'
            if ( assume == 'both' ) then
                print '(a,f16.3,a,f16.3,a,f16.3,a)',                       &
                    '# <tr><td colspan=2>&beta;:</td><td>',beta,          &
                    '</td><td>',beta25,'...',beta975,'</td></tr>'
            end if
        end if
    else
        print '(a,i5,a)','# Fitted to Gumbel distribution in ',iter,' iterations'
        if ( lnone ) then
            print '(a)','# P(x) = exp(-exp(-(x-a)/b))'
            print '(a,f16.3,a,f16.3,a,f16.3)','# a = ',a,' \\pm ',(a975-a25)/2
            print '(a,f16.3,a,f16.3,a,f16.3)','# b = ',b,' \\pm ',(b975-b25)/2
            if ( assume == 'scale' ) print '(a,f16.3,a,f16.3,a,f16.3)', &
                '# b/a = ',ba,' \\pm ',(ba975-ba25)/2
        else
            print '(a)','# P(x) = exp(-exp(-(x-a'')/b''))'
            call printab(restrain,lnone,lweb)
            call getabfromcov(a,b,alpha,beta,cov1,aaa,bbb)
            call getabfromcov(a25,b25,alpha,beta,cov1,aa25,bb25)
            call getabfromcov(a975,b975,alpha,beta,cov1,aa975,bb975)
            print '(a,i4,a,f16.3,a,f16.3,a,f16.3)','# ',yr1a,': a = ',aaa,' \\pm ',(aa975-aa25)/2
            print '(a,i4,a,f16.3,a,f16.3,a,f16.3)','# ',yr1a,': b = ',bbb,' \\pm ',(bb975-bb25)/2
            if ( assume == 'scale' ) print '(a,f16.3,a,f16.3,a,f16.3)', &
                '# b/a = ',ba,' \\pm ',(ba975-ba25)/2
            call getabfromcov(a,b,alpha,beta,cov2,aaa,bbb)
            call getabfromcov(a25,b25,alpha,beta,cov2,aa25,bb25)
            call getabfromcov(a975,b975,alpha,beta,cov2,aa975,bb975)
            print '(a,i4,a,f16.3,a,f16.3,a,f16.3)','# ',yr2a,': a = ',aaa,' \\pm ',(aa975-aa25)/2
            print '(a,i4,a,f16.3,a,f16.3,a,f16.3)','# ',yr2a,': b = ',bbb,' \\pm ',(bb975-bb25)/2
            print '(a,f16.3,a,f16.3,a,f16.3)','# alpha ',alpha,' \\pm ',(alpha975-alpha25)/2
            if ( assume == 'both' ) then
                print '(a,f16.3,a,f16.3,a,f16.3)','# beta  ',beta,' \\pm ',(beta975-beta25)/2
            end if
        end if
    end if
    call printcovreturntime(year,xyear,idmax,tx,tx25,tx975,yr1a,yr2a,yr2b,lweb,plot,assume,lnone,1)
    if ( .not. lnone ) call printcovpvalue(txtx,nmc,nmc,lweb,plot)
    call printcovreturnvalue(ntype,t,t25,t975,yr1a,yr2a,lweb,plot,assume,lnone)
    if ( .not. lnone ) call printcovreturntime(year,xyear,idmax,tx,tx25,tx975,yr1a,yr2a,yr2b,lweb,plot,assume,lnone,2)

    if ( dump ) then
        call plot_tx_cdfs(txtx,nmc,nmc,ntype,j1,j2)
    end if
    if ( plot ) write(11,'(3g20.4,a)') alpha,alpha25,alpha975,' alpha'
    call write_dthreshold(cov1,cov2,cov3,acov,offset,lchangesign)

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
        b = -b
        if ( .not. lnone ) then
            alpha = -alpha
            if ( assume == 'both' ) beta = -beta
        end if
    end if

    if ( lnone ) then
        call plotreturnvalue(ntype,t25(1,1),t975(1,1),j2-j1+1)
        ys(1:ntot) = yy(1:ntot)
        call plot_ordered_points(yy,ys,yrs,ntot,ntype,nfit,   &
            frac,a,b,xi,j1,j2,minindx,mindata,pmindata,      &
            year,xyear,snorm,lchangesign,lwrite,.true.)
    else
        ! compute distribution at past year and plot it
        call adjustyy(ntot,xx,assume,a,b,alpha,beta,cov1,      &
            yy,zz,aaa,bbb,lchangesign,lwrite)
        ys(1:ntot) = yy(1:ntot)
        print '(a,i5)','# distribution in year ',yr1a
        call plotreturnvalue(ntype,t25(1,1),t975(1,1),j2-j1+1)
        call plot_ordered_points(yy,ys,yrs,ntot,ntype,nfit,    &
            frac,aaa,bbb,xi,j1,j2,minindx,mindata,pmindata,   &
            year,xyear,snorm,lchangesign,lwrite,.false.)

        ! compute distribution at current year and plot it
        call adjustyy(ntot,xx,assume,a,b,alpha,beta,cov2,      &
            yy,zz,aaa,bbb,lchangesign,lwrite)
        ys(1:ntot) = yy(1:ntot)
        print '(a)'
        print '(a)'
        print '(a,i5)','# distribution in year ',yr2a
        if ( cov3 < 1e33 ) then
            last = .false.
        else
            last = .true.
        end if
        call plotreturnvalue(ntype,t25(1,2),t975(1,2),j2-j1+1)
        call plot_ordered_points(yy,ys,yrs,ntot,ntype,nfit,    &
            frac,aaa,bbb,xi,j1,j2,minindx,mindata,pmindata,   &
            year,xyear,snorm,lchangesign,lwrite,last)
        if ( cov3 < 1e33 ) then
            ! compute distribution at optional future year and plot it
            ! only plot the data points if they (almost) go up to that year
            ! (eg model data), but not if it is an extrapolation (eg obs)
            call adjustyy(ntot,xx,assume,a,b,alpha,beta,cov3,                   &
                yy,zz,aaa,bbb,lchangesign,lwrite)
            if ( cov3 > cmax + 0.1*(cmax-cmin) .or. cov3 < cmin - 0.1*(cmax-cmin) ) then
                 yy(1:ntot) = 3e33
            end if
            ys(1:ntot) = yy(1:ntot)
            print '(a)'
            print '(a)'
            print '(a,i5)','# distribution in year ',yr2b
            call plotreturnvalue(ntype,t25(1,4),t975(1,4),j2-j1+1)
            call plot_ordered_points(yy,ys,yrs,ntot,ntype,nfit,                 &
                frac,aaa,bbb,xi,j1,j2,minindx,mindata,pmindata,                &
                year,xyear,snorm,lchangesign,lwrite,.true.)
        end if
    end if

end subroutine

subroutine fit0gumcov(a,b,iter)
    use AmoebaToGSL
    implicit none
    integer iter
    real a,b,alpha,dalpha
    integer i
    real q(4),p(3,2),y(3),tol
    real llgumbelcov
    external llgumbelcov
!
!   fit, using Numerical Recipes routines
!
    q(1) = a
    q(2) = b
    q(3) = 0
    q(4) = 0
    p(1,1) = q(1) *0.9
    p(1,2) = q(2) *0.9
    p(2,1) = p(1,1) *1.2
    p(2,2) = p(1,2)
    p(3,1) = p(1,1)
    p(3,2) = p(1,2) *1.2
    do i=1,3
        q(1) = p(i,1)
        q(2) = p(i,2)
        q(3) = 0
        q(4) = 0
        y(i) = llgumbelcov(q)
    end do
    tol = 1e-4
    call amoeba(p,y,3,2,2,tol,llgumbelcov,iter)
!   maybe add restart later
    a = p(1,1)
    b = p(1,2)
end subroutine

subroutine fit1gumcov(a,b,alpha,dalpha,iter)
    use AmoebaToGSL
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

subroutine fit2gumcov(a,b,alpha,beta,dalpha,dbeta,iter)
    use AmoebaToGSL
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
    parameter(nmax=1000000)
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
    if ( llwrite ) print *,'a,b,alpha,llgumcov = ',p(1),p(2),p(3),llgumbelcov
!
999 continue
end function

subroutine gumbcovnorm(a,b,alpha,beta,s)
    implicit none
	include "getopts.inc"
    real a,b,alpha,beta,s
    real z1,z2

    if ( minindx > -1e33 .or. maxindx < 1e33 ) then
        write(0,*) 'gumbcovnorm: boundaries not yet available for fit of gumbel(t)'
        call exit(-1)
    else
        s = 1
    end if
!!!    print *,'gumbcovnorm: norm = ',a,b,alpha,s
end subroutine