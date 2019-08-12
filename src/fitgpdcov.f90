subroutine fitgpdcov(yrseries,yrcovariate,npernew,fyr,lyr &
             ,mens1,mens,crosscorr,a3,b3,xi3,alpha3,beta3,j1,j2,nens1,nens2 &
             ,lweb,ntype,lchangesign,yr1a,yr2a,yr2b,xyearin,idmax,cov1,cov2 &
             ,cov3,offset,t3,tx3,threshold,inrestrain,assume &
             ,confidenceinterval,ndecor,lboot,lprint,dump,plot,lwrite)
!
!       fit a GPD distribution to the data, which is already assumed to be declustered
!       input:
!       xx(2,ntot) data,covariate
!       j1,j2    use days/months/... j1 to j2
!       year     leave out this year from teh fit and compute return time for it
!       xyear    value for year, has been set to undef in the series
!       inrestrain restrain xi parameter by adding a normal distribution of width 0.5*inrestrain to the cost function
!       threshold in percent
!       assume   shift: only vary threshold, scale: vary threshold & b in unison, both: independently
!       output
!       a,b,xi,alpha,beta     parameters of fit
!       assume   shift: alpha modifies the position parameter a(cov) = a + alpha*cov
!                scale: alpha modifies both the position and shape parameters:
!                       a(cov) = a*exp(alpha*cov/a), b(cov) = b*exp(alpha*cov/a)
!                both:  a(cov) = a + alpha*cov, b(cov) = b + beta*cov
!       t(10,3)    return values for 10, 20, 50, ..., 10000 years for cov=cov1,cov2 and the difference
!       t25,t975   (100-confidenceinterval)/2%, (100+confidenceinterval)/25% quantiles of these return values
!       tx(3)      return time of the value of year (xyear) in the context of the other values and difference
!       t25,t975   (100-confidenceinterval)/2%, (100+confidenceinterval)/2% quantiles of these return values
!
    implicit none
!
    integer nmc
    parameter(nmc=1000)
    integer npernew,fyr,lyr,mens1,mens,ntot,j1,j2,nens1,nens2,ntype,yr1a,yr2a,yr2b,ndecor
    real yrseries(npernew,fyr:lyr,0:mens), &
             yrcovariate(npernew,fyr:lyr,0:mens),crosscorr(0:mens,0:mens), &
             a3(3),b3(3),xi3(3),alpha3(3),beta3(3),xyearin, &
             cov1,cov2,cov3,offset,inrestrain,t3(3,10,3),tx3(3,3),threshold, &
             confidenceinterval
    character assume*(*),idmax*(*)
    logical lweb,lchangesign,lboot,lprint,dump,plot,lwrite
!
    integer i,j,jj,k,l,n,nx,iter,iter1,iens,iiens,nfit,year,tsep,nj
    integer,allocatable :: yrs(:),bootyrs(:)
    real x,a,b,ba,xi,alpha,beta,t(10,4),t25(10,4),t975(10,4), &
             tx(4),tx25(4),tx975(4),aa(nmc),bb(nmc),xixi(nmc),baba(nmc), &
             alphaalpha(nmc),betabeta(nmc),tt(nmc,10,4),a25,a975, &
             b25,b975,ba25,ba975,xi25,xi975,alpha25,alpha975,beta25,beta975, &
             t5(10,4),t1(10,4),db,dxi,f,z,ll,ll1,txtx(nmc,4), &
             ranf,mean,sd,dalpha,dbeta,mindata,minindx,pmindata,snorm,s, &
             xmin,cmin,cmax,c,xxyear,frac,ttt(10,4),txtxtx(4), &
             acov(3,3),aacov(nmc,3),plo,phi,xyear,scross,sdecor, &
             aa25,aa975,bb25,bb975
    real adev,var,skew,curt,aaa,bbb,siga,chi2,q,p(4)
    integer,allocatable :: ii(:),yyrs(:)
    real,allocatable :: xx(:,:),yy(:),ys(:),zz(:),sig(:)
    logical lllwrite,lnone,last
    character lgt*4,string*1000,arg*250,method*3
!
    integer nmax,ncur
    parameter(nmax=1000000)
    real data(2,nmax),restrain
    logical llwrite,llchangesign
    common /fitdata3/ data
    common /fitdata2/ restrain,ncur,llwrite,llchangesign
    character cassume*5
    common /fitdata4/ cassume
    integer nthreshold
    real athreshold,pthreshold
    common /fitdata5/ nthreshold,athreshold,pthreshold
!
    real llgpdcov,gpdcovreturnlevel,gpdcovreturnyear
    external llgpdcov,gpdcovreturnlevel,gpdcovreturnyear
!
    allocate(yrs(0:nmax),bootyrs(0:nmax))
    allocate(xx(2,nmax))
    if ( cov1 == 0 .and. cov2 == 0 ) then
        lnone = .true.
    else
        lnone = .false.
    end if
    if ( lwrite ) print *,'fitgpdcov: calling fill_linear_array'
    call fill_linear_array(yrseries,yrcovariate,npernew,j1,j2, &
             fyr,lyr,mens1,mens,xx,yrs,nmax,ntot,lwrite)
    if ( lprint .and. lweb ) then
        print '(a,i9,a)','# <tr><td>N:</td><td>&nbsp;</td><td>', &
                 ntot,'</td><td>&nbsp;</td></tr>'
    end if
    if ( ntot < 6 ) then
        print '(a)','</table>'
        print '(a)','There are not enough points to fit a Generalised Pareto distribution. Please use a (much) longer time series.'
        print '(a)','<p>'
        call exit(-1)
    end if
    
    tsep = ndecor - 1
    if ( npernew >= 360 ) then
        tsep = -9999
        call decluster(xx,yrs,ntot,threshold,tsep,lwrite)
    end if

    year = yr2a
    pthreshold = threshold
    xyear = xyearin
    if ( lwrite ) then
        print *,'fitgpdcov: input:'
        print *,'assume         = ',assume
        print *,'j1,j2          = ',j1,j2
        print *,'threshold      = ',threshold,'%'
        print *,'year,xyear     = ',year,xyear
        print *,'cov1,cov2,offset ',cov1,cov2,offset
        if ( cov3 < 1e33 ) then
            print *,'cov3           = ',cov3
        end if
        print *,'ntot           = ',ntot
        print *,'ndecor         = ',ndecor
        if ( .false. ) then
            do i=1,ntot
                print *,i,(xx(j,i),j=1,2)
            end do
        end if
    end if
    a3 = 3e33
    b3 = 3e33
    xi3 = 3e33
    alpha3 = 3e33
    beta3 = 3e33
    t3 = 3e33
    tx3 = 3e33
    if ( ntot < 5 ) then
        return
    end if
!
!   thin data to exclude all the values equal to xmin (result of declustering)
!
    xmin = 3e33
    cmin = 3e33
    cmax = -3e33
    do i=1,ntot
        xmin = min(xmin,xx(1,i))
        cmin = min(cmin,xx(2,i))
        cmax = max(cmax,xx(2,i))
    end do
    ncur = 0
    do i=1,ntot
        if ( xx(1,i) > xmin ) then
            ncur = ncur + 1
        end if
    end do
    if ( ncur < 10 ) then
        return
    end if
!
!   dump (covariate,observation) pairs to plotfile on unit 15
!
    if ( dump ) call write_obscov(xx,yrs,ntot,xmin,cov2,xyear,year,offset,lchangesign)
!       
!   compute first-guess parameters
!
    allocate(ii(ncur))
    allocate(yyrs(0:ncur))
    allocate(yy(ncur))
    allocate(ys(ncur))
    allocate(zz(2*ncur))
    allocate(sig(ncur))
    j = 0
    yyrs(0) = yrs(0)
    do i=1,ntot
        if ( xx(1,i) > xmin ) then
            j = j + 1
            ii(j) = i
            yy(j) = xx(1,i)
            zz(j) = xx(2,i)
            yyrs(j) = yrs(i)
        end if
    end do
    if ( j /= ncur ) then
        write(0,*) 'fitgpdcov: error: j != ncur ',j,ncur
        write(*,*) 'fitgpdcov: error: j != ncur ',j,ncur
        call exit(-1)
    end if
    sig = 0
    call moment(yy,ncur,mean,adev,sd,var,skew,curt)
    call fit(zz,yy,ncur,sig,0,aaa,alpha,siga,dalpha,chi2,q)
    if ( lwrite ) then
        print *,'fitgpdcov: computed initialisation values:'
        print *,'mean,sd,alpha,dalpha = ',mean,sd,alpha,dalpha
    end if
!
!   ill-defined case
!
    if ( sd == 0 ) then
        goto 801 ! deallocate and return
    end if
!
!   copy to common for routine llgpdcov
!
!   number of points above threshold, threshold is relative to the full set (ntot)
    nthreshold = nint(ntot*(1-threshold/100))
    if ( nthreshold < 10 ) then
        if ( lprint ) then
            write(0,*) 'fitgpdcov: error: not enough points above threshold: ',nthreshold
            write(*,*) 'fitgpdcov: error: not enough points above threshold: ',nthreshold
        end if
        goto 801 ! deallocate and return
    end if
    if ( nthreshold >= ncur ) then
        write(0,*) 'fitgpdcov: error: nthreshold &gt ncur ',nthreshold,ncur
        write(*,*) 'fitgpdcov: error: nthreshold &gt ncur ',nthreshold,ncur
        nthreshold = ncur -1
    end if
    do i=1,ncur
        data(:,i) = xx(:,ii(i))
    end do
    restrain = inrestrain
    llwrite = .false. ! lwrite
    llchangesign = lchangesign
    cassume = assume

    ! first guess
    ys(1:ncur) = yy(1:ncur)
    call nrsort(ncur,yy)
    athreshold = (yy(ncur-nthreshold) + yy(ncur-nthreshold+1))/2
    if ( lwrite ) then
        print *,'fitgpdcov: nthreshold,ncur,ntot,athreshold = ',nthreshold,ncur,ntot,athreshold
    end if
    ! needed later on...
    yy(1:ncur) = ys(1:ncur)
    b = sd/3 ! should set the scale roughly right...
    if ( cassume == 'scale' .and. lchangesign ) then
        xi = -0.1
    else
        xi = 0.1
    end if
    if ( lnone ) then
        alpha = 3e33
        beta = 3e33
        call fit0gpdcov(a,b,xi,iter)
    else if ( assume == 'shift' .or. assume == 'scale' ) then
        beta = 3e33
        call fit1gpdcov(a,b,xi,alpha,dalpha,iter)
    else if ( assume == 'both' ) then
        beta = alpha
        dbeta = dalpha
        call fit2gpdcov(a,b,xi,alpha,beta,dalpha,dbeta,iter)
    else
        write(0,*) 'fitgpdcov: error: unknown value for assume ',assume
    end if
    if ( assume == 'scale' ) ba = b/a
    call getreturnlevels(a,b,xi,alpha,beta,cov1,cov2,cov3,gpdcovreturnlevel,j1,j2,assume,t)
    if ( xyear < 1e33 ) then
        call getreturnyears(a,b,xi,alpha,beta,xyear,cov1,cov2,cov3, &
             gpdcovreturnyear,j1,j2,tx,lchangesign,lwrite)
        ! convert to years
        do i=1,2
            if ( tx(i).lt.1e33 ) then
                tx(i) = tx(i)/(j2-j1+1)
            end if
        end do
    end if
    call getabfromcov(a,b,alpha,beta,cov1,aaa,bbb)
    if ( lchangesign) then
        acov(1,1) = -aaa
    else
        acov(1,1) = aaa
    end if
    call getabfromcov(a,b,alpha,beta,cov2,aaa,bbb)
    if ( lchangesign) then
        acov(1,2) = -aaa
    else
        acov(1,2) = aaa
    end if
    if ( cov3 < 1e33 ) then
        call getabfromcov(a,b,alpha,beta,cov3,aaa,bbb)
        if ( lchangesign) then
            acov(1,3) = -aaa
        else
            acov(1,3) = aaa
        end if
    end if
    if ( dump ) call write_threshold(cmin,cmax,a,b,xi,alpha,beta,offset,lchangesign,gpdcovreturnlevel)
!
!   bootstrap to find error estimates
!
    if ( .not.lboot ) then
        if ( lchangesign ) then
            a = -a
            t = -t
            if ( alpha < 1e33 ) then
                alpha = -alpha
                if ( beta < 1e33 ) then
                    beta = -beta
                end if
            end if
        end if
        a3(1) = a
        a3(2:3) = 3e33
        b3(1) = abs(b)
        b3(2:3) = 3e33
        xi3(1) = xi
        xi3(2:3) = 3e33
        alpha3(1) = alpha
        alpha3(2:3) = 3e33
        beta3(1) = beta
        beta3(2:3) = 3e33
        t3(1,:,1:3) = t(:,1:3)
        t3(2:3,:,:) = 3e33
        tx3(1,1:3) = tx(1:3)
        tx3(2:3,:) = 3e33            
        goto 801 ! deallocate and return
    end if
    if ( lprint .and. .not.lweb ) print '(a,i6,a)','# doing a ',nmc &
              ,'-member bootstrap to obtain error estimates'
    iens = 0
    scross = 0
    do iiens=1,nmc
        iens = iens + 1
        if ( lprint ) call keepalive1('Bootstrapping',iiens,nmc)
        if ( lprint .and. .not.lweb .and. mod(iiens,100) == 0 ) print '(a,i6)','# ',iiens
        n = 1 + (ncur-1)/ndecor
        method = 'new'
        if ( method == 'old' ) then
            do i=1,n
                ! we do not have the information here to check whether the
                ! data points were contiguous in the original series...
                ! TODO: propagate that information              
                call random_number(ranf)
                j = 1 + min(ncur-ndecor,int((ncur-ndecor)*ranf))
                if ( j < 1 .or. j > ncur ) then
                    write(0,*) 'fitgpdcov: error: j = ',j
                    call exit(-1)
                end if
                if ( i < n ) then ! the blocks that fit in whole
                    do jj=0,ndecor-1
                        data(:,1+(i-1)*ndecor+jj) = xx(:,ii(j+jj))
                    end do
                else
                    do jj=0,ndecor-1 ! one more block to the end, the previous block is shortened
                        data(:,1+ncur-ndecor+jj) = xx(:,ii(j+jj))
                    end do
                end if
            end do
        else
            ndecor = tsep + 1
            lllwrite = .false. ! lwrite
            call sample_bootstrap(yrseries,yrcovariate, &
                     npernew,j1,j2,fyr,lyr,nens1,nens2,crosscorr, &
                     ndecor,data,nmax,ntot,sdecor,lllwrite)
            if ( npernew >= 360 ) then
                bootyrs = -9999 ! cannot yet keep track of discontinuities, just hope they are not too bad
                ! or use the original ones
                ! use the same tsep as the first call, otherwise the answer is *very* different
                call decluster(data,yrs,ntot,threshold,tsep,lwrite)
            end if
            scross = scross + sdecor
        end if
        aa(iens) = a
        bb(iens) = b
        xixi(iens) = xi
        alphaalpha(iens) = alpha
        llwrite = .false.
        if ( lnone ) then
            alphaalpha(iens) = 3e33
            betabeta(iens) = 3e33
            call fit0gpdcov(aa(iens),bb(iens),xixi(iens),iter)
        else if ( assume == 'shift' .or. assume == 'scale' ) then
            betabeta(iens) = 3e33
            call fit1gpdcov(aa(iens),bb(iens),xixi(iens),alphaalpha(iens),dalpha,iter1)
        else if ( assume == 'both' ) then
            betabeta(iens) = beta
            call fit2gpdcov(aa(iens),bb(iens),xixi(iens), &
                     alphaalpha(iens),betabeta(iens),dalpha,dbeta,iter1)
        else
            write(0,*) 'fitgpdcov: error: unknown value for assume ',assume
        end if
        if ( assume == 'scale' ) baba(iens) = bb(iens)/aa(iens)
        call getabfromcov(aa(iens),bb(iens),alphaalpha(iens),betabeta(iens),cov1,aaa,bbb)
        aacov(iens,1) = aaa
        call getabfromcov(aa(iens),bb(iens),alphaalpha(iens),betabeta(iens),cov2,aaa,bbb)
        aacov(iens,2) = aaa
        if ( cov3 < 1e33 ) then
            call getabfromcov(aa(iens),bb(iens),alphaalpha(iens),betabeta(iens),cov3,aaa,bbb)
            aacov(iens,3) = aaa
        end if
        call getreturnlevels(aa(iens),bb(iens),xixi(iens), &
                 alphaalpha(iens),betabeta(iens), &
                 cov1,cov2,cov3,gpdcovreturnlevel,j1,j2,assume,ttt)
        do i=1,10
            do j=1,4
                tt(iens,i,j) = ttt(i,j)
            end do
        end do
        if ( xyear < 1e33 ) then
            call getreturnyears(aa(iens),bb(iens),xixi(iens), &
                     alphaalpha(iens),betabeta(iens),xyear,cov1,cov2,cov3, &
                     gpdcovreturnyear,j1,j2,txtxtx,lchangesign,lwrite)
            do j=1,4
                if ( j == 3 ) cycle
                ! convert to years
                txtx(iens,j) = txtxtx(j)/(j2-j1+1)
            end do
            txtx(iens,3) = txtxtx(3)
            ! if the event would have been impossible in the current climate
            ! disregard this bootstrap sample
            if ( txtxtx(2) > 1e19 ) then
                iens = iens - 1
            end if
            ! if the return time is not infinite but something went completely wrong
            ! also disregard it
            if ( txtxtx(1) > 1e30 ) then
                iens = iens - 1
            end if
        end if
    end do
    if ( mens > mens1 ) call print_spatial_scale(scross/(nmc*ndecor),ntot/real(mens-mens1+1)/real(ndecor))
    if ( lchangesign ) then
        if ( a < 1e33 ) a = -a
        if ( ba < 1e33 ) ba = -ba
        do iiens=1,iens
            if ( aa(iiens) < 1e33 ) aa(iiens) = -aa(iiens)
            if ( baba(iiens) < 1e33 ) baba(iiens) = -baba(iiens)
            do j=1,3
                if ( aacov(iiens,j) < 1e33 ) aacov(iiens,j) = -aacov(iiens,j)
            end do
        end do
        if ( alpha < 1e33 ) then
            alpha = -alpha
            alphaalpha = -alphaalpha
            if ( beta < 1e33 ) then
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
            end do
        end do
        do iiens=1,iens
            do j=1,10
                do i=1,4
                    if ( assume == 'scale' .and. i == 3 ) cycle ! a relative change does not change sign...
                    if ( tt(iiens,j,i) < 1e33 ) then
                        tt(iiens,j,i) = -tt(iiens,j,i)
                    end if
                end do
            end do
        end do
    end if
    plo = (100-confidenceinterval)/2
    phi = (100+confidenceinterval)/2
    call getcut( a25,plo,iens,aa)
    call getcut(a975,phi,iens,aa)
    call getcut( b25,plo,iens,bb)
    call getcut(b975,phi,iens,bb)
    if ( assume == 'scale' ) then
        call getcut( ba25,plo,iens,baba)
        call getcut(ba975,phi,iens,baba)
    end if
    call getcut( xi25,plo,iens,xixi)
    call getcut(xi975,phi,iens,xixi)
    if ( alpha < 1e33 ) then
        call getcut( alpha25,plo,iens,alphaalpha)
        call getcut(alpha975,phi,iens,alphaalpha)
        if ( beta < 1e33 ) then
            call getcut( beta25,plo,iens,betabeta)
            call getcut(beta975,phi,iens,betabeta)
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
                call getcut(t5(i,j),5.,iens,tt(1,i,j))
                call getcut(t1(i,j),1.,iens,tt(1,i,j))
            else
                lgt = '&gt;'
                call getcut(t5(i,j),95.,iens,tt(1,i,j))
                call getcut(t1(i,j),99.,iens,tt(1,i,j))
            end if
            call getcut( t25(i,j),plo,iens,tt(1,i,j))
            call getcut(t975(i,j),phi,iens,tt(1,i,j))
        end do
    end do
    do j=1,nj
        if ( xyear < 1e33 ) then
            call getcut(tx25(j),plo,iens,txtx(1,j))
            call getcut(tx975(j),phi,iens,txtx(1,j))
        end if
    end do
    if ( xyear < 1e33 ) then
        if ( lchangesign ) xyear = -xyear
    end if
!
!   output
!
    if ( .not.lprint ) then
        call copyab3etc(a3,b3,xi3,alpha3,beta3,t3,tx3, &
                 a,a25,a975,b,b975,xi,xi,xi,alpha,alpha25,alpha975, &
                 beta,beta25,beta975,t,t25,t975,tx,tx25,tx975)
        goto 801 ! deallocate and return
        !!!if ( .not.lwrite ) goto 801 ! deallocate and return
    end if
    if ( lweb ) then
        if ( lnone ) then
            print '(a)','# <tr><td colspan="4">Fitted to GPD '// &
                     'distribution H(x+&mu;) = 1 - (1+&xi;*x/&sigma;)^(-1/&xi;)</td></tr>'
            call printab(restrain,lnone,lweb)
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'// &
                     '&mu;:</td><td>',a,'</td><td>',a25,'...',a975,'</td></tr>'
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'// &
                     '&sigma;:</td><td>',b,'</td><td>',b25,'...',b975,'</td></tr>'
            if ( assume == 'scale' ) then
                print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'// &
                    '&sigma;/&mu;:</td><td>',ba,'</td><td>',ba25,'...',ba975,'</td></tr>'
            end if
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'// &
                     '&xi;:</td><td>',xi,'</td><td>',xi25,'...',xi975,'</td></tr>'
        else
            print '(a)','# <tr><td colspan="4">Fitted to GPD '// &
                     'distribution H(x+&mu;'') = 1 - (1+&xi;*x/&sigma;'')^(-1/&xi;)</td></tr>'
            call printab(restrain,lnone,lweb)
            call getabfromcov(a,b,alpha,beta,cov1,aaa,bbb)
            call getabfromcov(a25,b25,alpha,beta,cov1,aa25,bb25)
            call getabfromcov(a975,b975,alpha,beta,cov1,aa975,bb975)
            print '(a,i5,a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td>'//         &
                 '&mu;'':</td><td>',yr1a,'</td><td>',aaa,'</td><td>',       &
                 aa25,'...',aa975,'</td></tr>'
            print '(a,i5,a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td>'//         &
                 '&sigma;'':</td><td>',yr1a,'</td><td>',bbb,'</td><td>',    &
                 bb25,'...',bb975,'</td></tr>'
            call getabfromcov(a,b,alpha,beta,cov2,aaa,bbb)
            call getabfromcov(a25,b25,alpha,beta,cov2,aa25,bb25)
            call getabfromcov(a975,b975,alpha,beta,cov2,aa975,bb975)
            print '(a,i5,a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td>'//         &
                 '&mu;'':</td><td>',yr2a,'</td><td>',aaa,'</td><td>',       &
                 aa25,'...',aa975,'</td></tr>'
            print '(a,i5,a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td>'//         &
                '&sigma;'':</td><td>',yr2a,'</td><td>',bbb,'</td><td>',    &
                bb25,'...',bb975,'</td></tr>'
            if ( assume == 'scale' ) then
                print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'// &
                    '&sigma;/&mu;:</td><td>',ba,'</td><td>',ba25,'...',ba975,'</td></tr>'
            end if
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'// &
                    '&xi;:</td><td>',xi,'</td><td>',xi25,'...',xi975,'</td></tr>'
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'// &
                    '&alpha;:</td><td>',alpha,'</td><td>',alpha25,'...',alpha975,'</td></tr>'
            if ( assume == 'both' ) then
                print '(a,f16.3,a,f16.3,a,f16.3,a)', &
                        '# <tr><td colspan=2>&beta;:</td><td>',beta, &
                        '</td><td>',beta25,'...',beta975,'</td></tr>'
            end if
        end if
    else
        print '(a,i5,a)','# Fitted to GPD distribution in ',iter,' iterations'
        if ( lnone ) then
            print '(a)','# H(x+a) = 1-(1+xi*x/b)**(-1/xi)'
            call printab(restrain,lnone,lweb)
            print '(a,f16.3,a,f16.3,a,f16.3)','# a = ',a,' \\pm ',(a975-a25)/2
            print '(a,f16.3,a,f16.3,a,f16.3)','# b = ',b,' \\pm ',(b975-b25)/2
            if ( assume == 'scale' ) print '(a,f16.3,a,f16.3,a,f16.3)', &
                '# b/a = ',ba,' \\pm ',(ba975-ba25)/2
            print '(a,f16.3,a,f16.3,a,f16.3)','# xi  = ',xi,' \\pm ',(xi975-xi25)/2
        else
            print '(a)','# H(x+a'') = 1-(1+xi*x/b'')**(-1/xi)'
            call printab(restrain,lnone,lweb)
            call getabfromcov(a,b,alpha,beta,cov1,aaa,bbb)
            call getabfromcov(a25,b25,alpha,beta,cov1,aa25,bb25)
            call getabfromcov(a975,b975,alpha,beta,cov1,aa975,bb975)
            print '(a,i4,a,f16.3,a,f16.3,a,f16.3)','# ',yr1a,': a = ',aaa,' \\pm ',(aa975-aa25)/2
            print '(a,i4,a,f16.3,a,f16.3,a,f16.3)','# ',yr1a,': b = ',bbb,' \\pm ',(bb975-bb25)/2
            call getabfromcov(a,b,alpha,beta,cov2,aaa,bbb)
            call getabfromcov(a25,b25,alpha,beta,cov2,aa25,bb25)
            call getabfromcov(a975,b975,alpha,beta,cov2,aa975,bb975)
            print '(a,i4,a,f16.3,a,f16.3,a,f16.3)','# ',yr2a,': a = ',aaa,' \\pm ',(aa975-aa25)/2
            print '(a,i4,a,f16.3,a,f16.3,a,f16.3)','# ',yr2a,': b = ',bbb,' \\pm ',(bb975-bb25)/2
            if ( assume == 'scale' ) print '(a,f16.3,a,f16.3,a,f16.3)', &
                '# b/a  = ',ba,' \\pm ',(ba975-ba25)/2
            print '(a,f16.3,a,f16.3,a,f16.3)','# xi  = ',xi,' \\pm ',(xi975-xi25)/2
            print '(a,f16.3,a,f16.3,a,f16.3)','# alpha ',alpha,' \\pm ',(alpha975-alpha25)/2
            if ( assume == 'both' ) then
                print '(a,f16.3,a,f16.3,a,f16.3)','# beta  ',beta,' \\pm ',(beta975-beta25)/2
            end if
        end if
    end if
    call printcovreturntime(year,xyear,idmax,tx,tx25,tx975,yr1a,yr2a,yr2b,lweb,plot,assume,lnone,1)
    if ( .not. lnone ) call printcovpvalue(txtx,nmc,iens,lweb,plot)
    call printcovreturnvalue(ntype,t,t25,t975,yr1a,yr2a,lweb,plot,assume,lnone)
    if ( .not. lnone ) call printcovreturntime(year,xyear,idmax,tx,tx25,tx975,yr1a,yr2a,yr2b,lweb,plot,assume,lnone,2)

    if ( dump ) then
        call plot_tx_cdfs(txtx,nmc,iens,ntype,j1,j2)
    end if
    if ( plot ) write(11,'(3g20.4,a)') alpha,alpha25,alpha975,' alpha'
    if ( dump ) call write_dthreshold(cov1,cov2,cov3,acov,offset,lchangesign)

    ! no cuts
    mindata = -2e33
    minindx = -2e33
    pmindata = threshold
    snorm = 1
    frac = ncur/real(ntot)
    ! GPD fit
    nfit = 6

    do i=1,ncur
        data(:,i) = xx(:,ii(i))
    end do
    if ( lchangesign ) then
        ! we had flipped the sign on a,alpha but not yet on the rest, flip back...
        a = -a
        if ( alpha < 1e33 ) then
            alpha = -alpha
            if ( assume == 'both' ) beta = -beta
        end if
    end if
    
    if ( lnone ) then
        call plotreturnvalue(ntype,t25(1,1),t975(1,1),j2-j1+1)
        ys(1:ncur) = yy(1:ncur)
        mindata = a
        call plot_ordered_points(yy,ys,yrs,ncur,ntype,nfit,                 &
            frac,a,b,xi,j1,j2,minindx,mindata,pmindata,                &
            year,xyear,snorm,lchangesign,lwrite,.true.)
    else
    ! compute distribution at past year and plot it
        call adjustyy(ncur,data,assume,a,b,alpha,beta,cov1,yy,zz,aaa,bbb,lchangesign,lwrite)
        ys(1:ncur) = yy(1:ncur)
        print '(a,i5)','# distribution in year ',yr1a
        call plotreturnvalue(ntype,t25(1,1),t975(1,1),j2-j1+1)
        mindata = aaa
        call plot_ordered_points(yy,ys,yyrs,ncur,ntype,nfit, &
            frac,aaa,bbb,xi,j1,j2,minindx,mindata,pmindata, &
            year,xyear,snorm,lchangesign,lwrite,.false.)

        ! compute distribution at present year and plot it
        call adjustyy(ncur,data,assume,a,b,alpha,beta,cov2,yy,zz,aaa,bbb,lchangesign,lwrite)
        ys(1:ncur) = yy(1:ncur)
        print '(a)'
        print '(a)'
        print '(a,i5)','# distribution in year ',yr2a
        if ( cov3 < 1e33 ) then
            last = .false.
        else
            last = .true.
        end if
        call plotreturnvalue(ntype,t25(1,2),t975(1,2),j2-j1+1)
        mindata = aaa
        call plot_ordered_points(yy,ys,yyrs,ncur,ntype,nfit, &
            frac,aaa,bbb,xi,j1,j2,minindx,mindata,pmindata, &
            year,xyear,snorm,lchangesign,lwrite,last)
        if ( cov3 < 1e33 ) then
            ! compute distribution at present year and plot it
            call adjustyy(ncur,data,assume,a,b,alpha,beta,cov3,yy,zz,aaa,bbb,lchangesign,lwrite)
            if ( cov3 > cmax + 0.1*(cmax-cmin) .or. cov3 < cmin - 0.1*(cmax-cmin) ) then
                 yy(1:ncur) = 3e33
            end if
            ys(1:ncur) = yy(1:ncur)
            print '(a)'
            print '(a)'
            print '(a,i5)','# distribution in year ',yr2b
            call plotreturnvalue(ntype,t25(1,4),t975(1,4),j2-j1+1)
            mindata = aaa
            call plot_ordered_points(yy,ys,yyrs,ncur,ntype,nfit, &
                frac,aaa,bbb,xi,j1,j2,minindx,mindata,pmindata, &
                year,xyear,snorm,lchangesign,lwrite,.true.)
        end if
    end if
801 continue
    deallocate(ii)
    deallocate(yyrs)
    deallocate(yy)
    deallocate(ys)
    deallocate(zz)
    deallocate(sig)
    return
end subroutine

subroutine fit0gpdcov(a,b,xi,iter)
    use AmoebaToGSL
    implicit none
    integer iter
    real a,b,xi
    integer i,j,n
    real q(4),p(3,2),y(3),tol
    logical lok
    real,external :: llgpdcov
    integer ncur
    real restrain
    logical llwrite,llchangesign
    common /fitdata2/ restrain,ncur,llwrite,llchangesign
    character cassume*5
    common /fitdata4/ cassume
    integer nthreshold
    real athreshold,pthreshold
    common /fitdata5/ nthreshold,athreshold,pthreshold
!
    n = 0
10  continue
    q(1) = b
    q(2) = xi
    q(3) = 0
    q(4) = 0
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
        y(i) = llgpdcov(q)
        if ( y(i) < 1e33 ) lok = .true.
    end do
    if ( .not.lok ) then
        n = n + 1
        if ( llwrite ) then
            print *,n,'searching for good initial values, these all give 3e33'
            do i=1,3
                print *,(p(i,j),j=1,2),y(i)
            enddo
        end if
        if ( xi /= 0 ) then
            if ( cassume == 'scale' ) then
                xi = xi + 0.03*xi/abs(xi)
            else
                xi = xi + 0.03
            end if
        else if ( cassume == 'scale' ) then
            if ( athreshold < 0 ) then
                xi = -0.1
            else
                xi = +0.1 ! most precip has positive \xi
            end if
        else
            xi = -0.1 ! most temperature extremes negative
        end if
        if ( abs(xi).lt.2 ) then
            goto 10
        else
            write(0,*) 'fit1gpdcov: error: cannot find initial fit values'
            a = 3e33
            b = 3e33
            xi = 3e33
        end if
    end if
    tol = 1e-4
    call amoeba(p,y,3,2,2,tol,llgpdcov,iter)
!   maybe add restart later
    a = athreshold
    b = p(1,1)
    xi = p(1,2)
end subroutine

subroutine fit1gpdcov(a,b,xi,alpha,dalpha,iter)
    use AmoebaToGSL
    implicit none
    integer iter
    real a,b,xi,alpha,dalpha
    integer i,j,n
    real q(4),p(4,3),y(4),tol
    logical lok
    real,external :: llgpdcov
    integer ncur
    real restrain
    logical llwrite,llchangesign
    common /fitdata2/ restrain,ncur,llwrite,llchangesign
    character cassume*5
    common /fitdata4/ cassume
    integer nthreshold
    real athreshold,pthreshold
    common /fitdata5/ nthreshold,athreshold,pthreshold
!
    n = 0
10  continue
    q(1) = b
    q(2) = xi
    q(3) = alpha
    q(4) = 3e33
    p(1,1) = q(1) *0.9
    p(1,2) = q(2) *0.9
    p(1,3) = q(3) - dalpha
    p(2,1) = p(1,1) *1.2
    p(2,2) = p(1,2)
    p(2,3) = p(1,3)
    p(3,1) = p(1,1)
    p(3,2) = p(1,2) *1.2 + 0.1
    p(3,3) = p(1,3)
    p(4,1) = p(1,1)
    p(4,2) = p(1,2)
    p(4,3) = p(1,3) + 2*dalpha
    lok = .false.
    do i=1,4
        q(1) = p(i,1)
        q(2) = p(i,2)
        q(3) = p(i,3)
        y(i) = llgpdcov(q)
        if ( y(i) < 1e33 ) lok = .true.
    end do
    if ( .not.lok ) then
        n = n + 1
        if ( llwrite ) then
            print *,n,'searching for good initial values, these all give 3e33'
            do i=1,4
                print *,(p(i,j),j=1,3),y(i)
            enddo
        end if
        if ( xi /= 0 ) then
            if ( cassume == 'scale' ) then
                xi = xi + 0.03*xi/abs(xi)
            else
                xi = xi + 0.03
            end if
        else if ( cassume == 'scale' ) then
            if ( athreshold < 0 ) then
                xi = -0.1
            else
                xi = +0.1 ! most precip has positive \xi
            end if
        else
            xi = -0.1 ! most temperature extremes negative
        end if
        if ( abs(xi).lt.2 ) then
            if ( llwrite ) print *,'trying again with xi = ',xi
            goto 10
        else
            write(0,*) 'fit1gpdcov: error: cannot find initial fit values'
            a = 3e33
            b = 3e33
            xi = 3e33
            alpha = 3e33
        end if
    end if
    tol = 1e-4
    call amoeba(p,y,4,3,3,tol,llgpdcov,iter)
!   maybe add restart later
    a = athreshold
    b = p(1,1)
    xi = p(1,2)
    alpha = p(1,3)
end subroutine

subroutine fit2gpdcov(a,b,xi,alpha,beta,dalpha,dbeta,iter)
    use AmoebaToGSL
    implicit none
    integer iter
    real a,b,xi,alpha,beta,dalpha,dbeta
    integer i
    real q(4),p(5,4),y(5),tol
    logical lok
    real llgpdcov
    external llgpdcov
    character cassume*5
    common /fitdata4/ cassume
    integer nthreshold
    real athreshold,pthreshold
    common /fitdata5/ nthreshold,athreshold,pthreshold
!
10  continue
    q(1) = b
    q(2) = xi
    q(3) = alpha
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
    p(3,2) = p(1,2) *1.2 + 0.1
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
    lok = .false.
    do i=1,5
        q(1) = p(i,1)
        q(2) = p(i,2)
        q(3) = p(i,3)
        q(4) = p(i,4)
        y(i) = llgpdcov(q)
        if ( y(i) < 1e33 ) lok = .true.
    end do
    if ( .not.lok ) then
        if ( xi /= 0 ) then
            xi = xi + 0.03*xi/abs(xi)
        else if ( cassume == 'scale' .and. athreshold < 0 ) then
            xi = -0.1
        else 
            xi = 0.1
        end if
        if ( abs(xi).lt.1 ) then
            goto 10
        else
            write(0,*) 'fit2gpdcov: error: cannot find initial fit values'
            a = 3e33
            b = 3e33
            xi = 3e33
            alpha = 3e33
            beta = 3e33
        end if
    end if
    tol = 1e-4
    call amoeba(p,y,5,4,4,tol,llgpdcov,iter)
!   maybe add restart later
    a = athreshold
    b = p(1,1)
    xi = p(1,2)
    alpha = p(1,3)
    beta = p(1,4)
end subroutine

real function llgpdcov(p)
!
!       computes the log-likelihood function for a covariant-dependent GPD distribution
!       with parameters a,b,xi,alpha=p(1-4) and data in common.
!
    implicit none
    integer maxloop
    parameter(maxloop=50)
!
    real p(4)
!
    integer i,j,n,nlo,nhi,iloop,ilohi
    real x,z,xi,s,aa,bb,llold,alpha,alo,ahi,delta,aathreshold(2),arg
    real,allocatable :: xx(:),xxunsorted(:)
!
    integer nmax,ncur
    parameter(nmax=1000000)
    real data(2,nmax),restrain
    logical llwrite,llchangesign
    common /fitdata3/ data
    common /fitdata2/ restrain,ncur,llwrite,llchangesign
    character cassume*5
    common /fitdata4/ cassume
    integer nthreshold
    real athreshold,pthreshold
    common /fitdata5/ nthreshold,athreshold,pthreshold
    integer,save :: init=0
!
    llgpdcov = 0
    allocate(xx(ncur))
    allocate(xxunsorted(ncur))
    xi = p(2)
    if ( abs(xi) > 10 ) then
        if ( llwrite ) print *,'llgpdcov: |xi|>10: ',p(2)
        llgpdcov = 3e33
        goto 999 
    end if
    if ( restrain < 0 ) then
        write(0,*) 'llgpdcov: error: restrain<0 ',restrain
        call exit(-1)
    end if
!
!   get threshold
!
    alpha = p(3)
    if ( cassume /= 'scale' ) then
        do i=1,ncur
            xx(i) = data(1,i) - alpha*data(2,i)
            !!!print *,'xx(',i,') = ',xx(i),data(1,i),alpha*data(2,i)
        end do
        if ( llwrite ) xxunsorted = xx
        call nrsort(ncur,xx)
        athreshold = (xx(ncur-nthreshold) + xx(ncur-nthreshold+1))/2
        if ( llwrite ) then
            print *,'llgpdcov: nthreshold,athreshold = ',nthreshold,athreshold
            print *,'          xx(',ncur-nthreshold,') = ',xx(ncur-nthreshold)
            print *,'          xx(',ncur-nthreshold+1,') = ',xx(ncur-nthreshold+1)
        end if
        xx = xx - athreshold
        if ( llwrite ) xxunsorted = xxunsorted - athreshold
    else
        if ( llchangesign .and. init == 0 ) then
            init = 1
            write(0,*) 'Enforcing a hard lower bound of zero.'
        end if
        if ( llchangesign .and. xi >= 0 ) then
            ! scaling implies (for climate) that the distribution cannot cross zero, 
            ! so xi must be < 0
            llgpdcov = 3e33
            goto 999
        end if
        ! iterative procedure I am afraid...
        ! assume athreshold has been set to a reasonable value higher up
        if ( llwrite ) then
            print *,'llgpdcov: iterative threshold determination'
            print *,'athreshold,alpha = ',athreshold,alpha
        end if
        do ilohi=1,2 ! find threshold values for which n jumps to nthreshold and away from it
            nlo = -1
            nhi = -1
            alo = -3e33
            ahi = +3e33
            delta = 1.2
            aathreshold(ilohi) = athreshold
            do iloop=1,maxloop
                n = 0
                do i=1,ncur
                    arg = alpha*data(2,i)/aathreshold(ilohi)
                    if ( arg < 70 ) then
                        xx(i) = data(1,i) - aathreshold(ilohi)*exp(arg)
                    else
                        xx(i) = -3e33
                    end if
                    if ( xx(i) > 0 ) then
                        n = n + 1
                    end if
                end do
                if ( nhi == -1 .or. nlo == -1 ) then
                    ! bracket solution
                    if ( n <= nthreshold+1-ilohi .eqv. athreshold > 0 ) then
                        ahi = aathreshold(ilohi)
                        nhi = n
                        aathreshold(ilohi) = aathreshold(ilohi) / delta
                    else
                        alo = aathreshold(ilohi)
                        nlo = n
                        aathreshold(ilohi) = aathreshold(ilohi) * delta
                    end if
                    if ( .false. .and. llwrite ) print *,'bracketing ',ilohi,nthreshold,alo,nlo,ahi,nhi
                else
                    ! bisect
                    if ( n <= nthreshold+1-ilohi .eqv. athreshold > 0 ) then
                        ahi = aathreshold(ilohi)
                        nhi = n
                    else
                        alo = aathreshold(ilohi)
                        nlo = n
                    end if
                    aathreshold(ilohi) = (alo+ahi)/2
                    if ( abs((alo-ahi)/(alo+ahi)).lt.2e-7 ) exit
                    if ( .false. .and. llwrite ) print *,'bisecting  ',ilohi,nthreshold,alo,nlo,ahi,nhi
                end if
            end do
            if ( iloop > maxloop ) then
                if ( llwrite ) then
                    write(*,'(a,i3,2i5,f10.2)') '# llgpdcov: warning:  ' &
         &                   //'threshold computation did not converge ', &
         &                   iloop,n,nthreshold,aathreshold(ilohi)
                end if
            end if
            if ( llwrite ) print *,'aathreshold(ilohi) = ',aathreshold(ilohi)
        end do ! ilohoi
        athreshold = (aathreshold(1) + aathreshold(2))/2
        if ( llchangesign .and. p(1) > athreshold*xi ) then
            ! scaling implies (for climate) that the distribution cannot cross zero, 
            ! so the upper limit must be <0 (for flipped signs)
            if ( llwrite ) print *,'GPD limit below zero, invalid: ',p(1),' > ',athreshold,'*',xi
            llgpdcov = 3e33
            goto 999
        else
            if ( llwrite ) print *,'GPD limit above zero, OK: ',p(1),' <= ',athreshold,'*',xi
        end if
        n = 0
        do i=1,ncur
            arg = alpha*data(2,i)/athreshold
            if ( arg < 70 ) then
                xx(i) = data(1,i) - athreshold*exp(arg)
            else if ( athreshold < 0 ) then
                xx(i) = 3e33
            else
                xx(i) = -3e33
            end if
            if ( xx(i) >= 0 ) then
                n = n + 1
            end if
        end do
        if ( llwrite ) print *,'athreshold = ',athreshold,n
        call nrsort(ncur,xx)
    end if
    if ( xx(ncur) > 1e33 ) then
        if ( llwrite ) print *,'llgpdcov: x'' undefined '
        llgpdcov = 3e33
        goto 999
    end if
!
!   and compute cost function
!
    do i=ncur-nthreshold+1,ncur
        call getabfromcov(athreshold,p(1),p(3),p(4),data(2,i),aa,bb)
        if ( abs(aa).gt.1e33 ) then
            if ( llwrite ) print *,'llgpdcov: aa undefined ',aa
            llgpdcov = 3e33
            goto 999            
        end if 
        if ( abs(bb) < 1e-15 ) then
            if ( llwrite ) print *,'llgpdcov: bb too small ',bb
            llgpdcov = 3e33
            goto 999
        end if
        if ( bb < 0 ) then
            if ( llwrite ) print *,'llgpdcov: bb negative ',bb
            llgpdcov = 3e33
            goto 999
        end if
        if ( cassume == 'scale' .and. ( llchangesign .eqv. athreshold > 0 ) ) then
            if ( llwrite ) print *,'llgpdcov: threshold has wrong sign ',athreshold
            llgpdcov = 3e33
            goto 999
        end if

        z = xx(i)
        if ( z < 0 ) then
            if ( z < -5e-5*abs(athreshold) ) then
                write(0,*) 'llgpdcov: error: z<0 ',z,n,nthreshold
                write(*,*) '# llgpdcov: error: z<0 ',z,n,nthreshold
                if ( llwrite ) then
                    do j=1,ncur
                        print *,i,j,xx(j)
                    end do
                end if
            end if
            z = 0
        end if
        s = 1+xi*z/bb
        if ( s <= 0 ) then
            if ( llwrite ) print *,'llgpdcov: 1+xi*z/bb < 0 ',i,s,xi,z,bb
            llgpdcov = 3e33
            goto 999
        end if
        llold = llgpdcov
        if ( abs(xi) < 1e-4 ) then
            llgpdcov = llgpdcov - z/bb + (z/bb)**2*xi/2
!!!                print *,i,z, - z/b + (z/b)**2*xi/2 - log(b)
        else
            llgpdcov = llgpdcov - (1+1/xi)*log(1+xi*z/bb)
!!!                print *,i,z, - (1+1/xi)*log(1+xi*z/b) - log(b)
        end if
        if ( llwrite ) then
            do j=1,ncur
                if ( xxunsorted(j) == xx(i) ) exit
            end do
            print *,i,data(1,j),data(2,j),z,bb,llgpdcov-llold
        end if
        llgpdcov = llgpdcov - log(bb)
    end do
!   normalization is not 1 in case of cut-offs
    call gpdcovnorm(athreshold,abs(p(1)),p(2),p(3),p(4),s)
    if ( s < 1e33 ) then
        llgpdcov = llgpdcov - ncur*log(s)
    else
        if ( llwrite ) print *,'cannot handle cuts yet'
        llgpdcov = 3e33
        goto 999
    end if
    if ( restrain /= 0 ) then
!       preconditioning on xi with gaussian of width restrain/2 around 0
        llgpdcov = llgpdcov - (xi/(restrain/2))**2/2
    end if
!   minimum, not maximum
    llgpdcov = -llgpdcov
!
999 continue
    deallocate(xx)
    if ( llwrite ) print *,'a,b,xi,alpha,llgpdcov = ',athreshold,p(1),p(2),p(3),llgpdcov
end function

subroutine gpdcovnorm(a,b,xi,alpha,beta,s)
implicit none
include 'getopts.inc'
    real a,b,xi,alpha,beta,s
    real z1,z2

    if ( minindx > -1e33 .or. maxindx < 1e33 ) then
        write(0,*) 'gpdcovnorm: boundaries not yet avaiable for fit of GPD(t)'
        call exit(-1)
    else
        s = 1
    end if
!!!    print *,'gpdcovnorm: norm = ',a,b,s
end subroutine

real function gpdcovreturnlevel(a,b,xi,alpha,beta,xx,cov)
!
!       compute return times given the GPD distribution parameters a,b,xi and 
!       x = log10(returntime) for covariant cov and fit parameter alpha
!       Uses a few Taylor series approximation for xi small and/or return time large
!
    implicit none
    real a,b,xi,alpha,beta,xx,cov
    real aa,bb,x,y,t
    logical lwrite
    integer ncur
    real restrain
    logical llwrite,llchangesign
    common /fitdata2/ restrain,ncur,llwrite,llchangesign
    character cassume*5
    common /fitdata4/ cassume
    integer nthreshold
    real athreshold,pthreshold
    common /fitdata5/ nthreshold,athreshold,pthreshold

    call getabfromcov(a,b,alpha,beta,cov,aa,bb)
    if ( abs(bb) < 1e-10 .or. bb.gt.1e33 ) then
        gpdcovreturnlevel = 3e33
        return
    end if
    x = xx + log10(1-pthreshold/100)
    if ( abs(xi) > 10 ) then
        gpdcovreturnlevel = 3e33
    else if ( abs(xi) < 1e-4 ) then
        t = bb*x*log(10.) + 0.5*xi*(x*log(10.))**2
    else
        y = xi*x*log(10.)
        if ( y < 46 ) then
            t = bb/xi*(-1 + exp(y))
        else
            t = 1e20
        end if
    end if
    t = t + aa ! threshold
    if ( llwrite ) then
        print *,'gpdcovreturnlevel:'
        print *,'a,b,xi,alpha,beta = ',a,b,xi,alpha,beta
        print *,'t,cov             = ',exp(x*log(10.)),cov
        print *,'aa,bb             = ',aa,bb
        print *,'gpdcovreturnlevel = ',t
    end if
    if ( cassume == 'scale' .and. ( llchangesign .eqv. t > 0 ) ) then
        write(0,*) 'gpdcovreturnlevel: error: wrong sign return level ',t
        t = 3e33
    end if
    gpdcovreturnlevel = t
end function

real function gpdcovreturnyear(a,b,xi,alpha,beta,xyear,cov,lchangesign)
!
!   compute the return time of the value xyear with the fitted values
!
    implicit none
    real a,b,xi,alpha,beta,xyear,cov
    integer i,n,ntot
    real x,y,z,tx,aa,bb,aaa,bbb,p
    real,allocatable :: yy(:),zz(:)
    logical lchangesign
    integer nmax,ncur
    parameter(nmax=1000000)
    real data(2,nmax),restrain
    logical llwrite,llchangesign
    common /fitdata3/ data
    common /fitdata2/ restrain,ncur,llwrite,llchangesign
    character cassume*5
    common /fitdata4/ cassume
    integer nthreshold
    real athreshold,pthreshold
    common /fitdata5/ nthreshold,athreshold,pthreshold

    llwrite = .false.
    if ( llwrite ) then
        print *,'gpdcovreturnyear: input'
        print *,'a,b,xi     = ',a,b,xi
        print *,'alpha,beta = ',alpha,beta
        print *,'xyear,cov  = ',xyear,cov,lchangesign
    end if
    x = xyear
    call getabfromcov(a,b,alpha,beta,cov,aa,bb)
    if ( llwrite ) print *,'aa,bb = ',aa,bb
    if ( xyear > aa ) then
        x = xyear - aa
        z = (1 + xi*x/bb)
        if ( z > 1e10 ) then
            !!!write(0,*) 'gpdcovreturnyear: error: z = ',z,xi,x,bb
            tx = 3e33
        else if ( z <= 0 ) then
            tx = 1e20
        else if ( abs(xi) < 1e-3 ) then
            tx = exp(x/bb - 0.5*xi*(x/bb)**2)/(1-pthreshold/100)
        else if ( log(z)/xi > 45 ) then
            tx = 1e20
        else
            tx = z**(1/xi)/(1-pthreshold/100)
        end if
        if ( .false. .and. tx > 1e20 ) then
            write(0,*) 'fitgpd: tx > 1e20: ',tx
            write(0,*) 'z,xi = ',z,xi
        end if
        if ( llwrite ) then
            print *,'gpdcovreturnyear: tx from GPD fit = ',tx
            print *,'  aa,pthreshold = ',aa,pthreshold
            print *,'  x,z = ',x,z
        end if
    else
        n = 0
        allocate(yy(ncur),zz(ncur))
        call adjustyy(ncur,data,cassume,a,b,alpha,beta,cov,yy,zz,aaa,bbb,lchangesign,llwrite)
        if ( .true. ) then
            do i=1,ncur
                if ( yy(i) > xyear ) n = n + 1
            end do
        else
            if ( cassume == 'shift' ) then
                do i=1,ncur
                    if ( data(1,i) - alpha*(data(2,i)-cov) > xyear ) n = n + 1
                end do
            else if ( cassume == 'scale' ) then
                do i=1,ncur
                    if ( data(1,i)*exp(alpha/athreshold*(data(2,i)-cov)) > xyear ) n = n + 1
                end do
            else if ( cassume == 'both' ) then
                do i=1,ncur
                    ! not sure whether this is correct
                    if ( (data(1,i)-xyear - alpha*(data(2,i)-cov)) > 0 ) n = n + 1 
                end do
            else
                write(0,*) 'gpdcovreturnyear: error: unknown value for assume: ',cassume
                call exit(-1)
            end if
        end if
        ! approximately... I do not have this information here.
        ntot = nint(nthreshold/(1-pthreshold/100))
        if ( n < nthreshold/2 ) then
            !!!write(0,*) 'gpdcovreturnyear: error: n<nthreshold/2 but xyear < aa ',n,nthreshold/2,xyear,aa
            if ( llwrite ) write(*,*) 'gpdcovreturnyear: error: n<nhreshold/2 but xyear < aa ',n,nthreshold/2,xyear,aa
            tx = 3e33
        else
            tx = real(ntot+1)/real(n)
        end if
        if ( llwrite ) then
            print *,'gpdcovreturnyear: interpolated tx = ',tx
            print *,'  n,ntot = ',n,ntot
        end if
        deallocate(yy,zz)
    end if
    if ( .true. .and. tx < 1e-10 ) then
        write(0,*) 'gpdcovreturnyear: tx < 1e-10: ',tx
        write(0,*) 'a,b,xi,alpha,xyear = ',a,b,alpha,xi,xyear
        tx = 3e33
    end if
    gpdcovreturnyear = tx
end function
