subroutine fitgevcov(yrseries,yrcovariate,npernew,fyr,lyr                   &
             ,mens1,mens,crosscorr,a3,b3,xi3,alpha3,beta3,j1,j2,nens1,nens2 &
             ,lweb,ntype,lchangesign,yr1a,yr2a,yr2b,xyear,idmax,cov1,cov2   &
             ,cov3,offset,t3,tx3,inrestrain,assume,confidenceinterval       &
             ,ndecor,lboot,lprint,dump,plot,lwrite)
!
!       fit a GEV distribution to the data, which is already assumed to be block max
!       input:
!       xx(2,ntot) data,covariate
!       j1,j2    use days/months/... j1 to j2
!       year     leave out this year from the fit and compute return time for it
!       xyear    value for year, has been set to undef in the series
!       inrestrain restrain xi parameter by adding a normal distribution of width 0.5*inrestrain to the cost function
!       assume   shift: only vary a, scale: vary a & b in unison, both: independently
!       output
!       a,b,xi,alpha,beta  parameters of fit
!       assume   shift: alpha modifies the position parameter a(cov) = a + alpha*cov
!                scale: alpha modifies both the position and shape parameters:
!                       a(cov) = a*exp(alpha*cov), b(cov) = b*exp(alpha*cov)
!                both:  a(cov) = a + alpha*cov, b(cov) = b + beta*cov
!       t(10,4)    return values for 10, 20, 50, ..., 10000 years for cov=cov1,cov2,the difference
!                  and optionally cov3
!       t25,t975   (100-confidenceinterval)/2%, (100+confidenceinterval)/2% quantiles of these return values
!       tx(3)      return time of the value of year (xyear) in the context of the other values and difference
!       tx25,tx975 (100-confidenceinterval)/2%, (100+confidenceinterval)/2% quantiles of these return times and their differences
!
    implicit none
    integer nmc
    parameter(nmc=1000)
    integer npernew,fyr,lyr,mens1,mens,ntot,j1,j2,nens1,nens2,ntype,yr1a,yr2a,yr2b,ndecor
    real yrseries(npernew,fyr:lyr,0:mens),yrcovariate(npernew,fyr:lyr,0:mens),   &
 &       crosscorr(0:mens,0:mens),a3(3),b3(3),xi3(3),alpha3(3),beta3(3),xyear,   &
 &       cov1,cov2,cov3,offset,inrestrain,t3(3,10,3),tx3(3,3),confidenceinterval
    character assume*(*),idmax*(*)
    logical lweb,lchangesign,lboot,lprint,dump,plot,lwrite
!
    integer i,j,k,l,n,nx,iter,iens,jens,iiens,nfit,year,yr,nj
    integer,allocatable :: yrs(:)
    real x,a,b,ba,xi,alpha,beta,t(10,4),t25(10,4),t975(10,4),           &
 &       tx(4),tx25(4),tx975(4),aa(nmc),bb(nmc),baba(nmc),xixi(nmc),    &
 &       alphaalpha(nmc),betabeta(nmc),tt(nmc,10,4),b25,b975,ba25,ba975 &
 &       ,xi25,xi975,alpha25,alpha975,beta25,beta975,t5(10,4),t1(10,4)  &
 &       ,db,dxi,f,threshold,thens,z,ll,ll1,txtx(nmc,4)                 &
 &       ,a25,a975,ranf,mean,sd,dalpha,dbeta                            &
 &       ,mindata,minindx,pmindata,snorm,s,xxyear,frac,                 &
 &       acov(3,3),aacov(nmc,3),plo,phi,cmin,cmax,scross,sdecor
    real ttt(10,4),txtxtx(4)
    real adev,var,skew,curt,aaa,bbb,aa25,aa975,bb25,bb975,siga,chi2,q
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
    real llgevcov,gevcovreturnlevel,gevcovreturnyear
    external llgevcov,gevcovreturnlevel,gevcovreturnyear
!
    allocate(yrs(0:nmax))
    allocate(xx(2,nmax))
    if ( cov1 == 0 .and. cov2 == 0 ) then
        lnone = .true.
    else
        lnone = .false.
    end if
    year = yr2a

    if ( lwrite ) print *,'fitgevcov: calling fill_linear_array'
    call fill_linear_array(yrseries,yrcovariate,npernew,j1,j2,          &
         fyr,lyr,mens1,mens,xx,yrs,nmax,ntot,lwrite)
    if ( lprint .and. lweb ) then
        print '(a,i9,a)','# <tr><td>N:</td><td>&nbsp;</td><td>',        &
             ntot,'</td><td>&nbsp;</td></tr>'
    end if
    if ( ntot < 6 ) then
        print '(a)','</table>'
        print '(a)','There are not enough points to fit a GEV distribution. Please use a (much) longer time series.'
        print '(a)','<p>'
        call exit(-1)
    end if

    if ( lwrite ) then
        print *,'fitgevcov: input:'
        print *,'assume         = ',assume
        print *,'j1,j2          = ',j1,j2
        print *,'year,xyear     = ',year,xyear
        print *,'cov1,cov2,offset ',cov1,cov2,offset
        if ( cov3 < 1e33 ) then
            print *,'cov3           = ',cov3
        end if
        if ( .true. ) then
            do i=1,ntot
                print *,i,(xx(j,i),j=1,2)
            enddo
        endif
    endif
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
    sig = 0
    call moment(yy,ntot,mean,adev,sd,var,skew,curt)
!
!   ill-defined case
!
    if ( sd == 0 ) then
        a3 = 3e33
        b3 = 3e33
        xi3 = 3e33
        alpha3 = 3e33
        beta3 = 3e33
        t3 = 3e33
        tx3 = 3e33
        return
    endif

!   initialisation values

    if ( .not.lnone ) then
        call fit(zz,yy,ntot,sig,0,aaa,alpha,siga,dalpha,chi2,q)
        if ( lwrite ) then
            print *,'fitgevcov: computed initialisation values:'
            print *,'mean,sd,alpha,dalpha = ',mean,sd,alpha,dalpha
        end if
    end if
!
!   copy to common for routine llgevcov
!
    ncur = ntot
    do i=1,ncur
        data(:,i) = xx(:,i)
    enddo
    restrain = inrestrain
    llwrite = lwrite
    llchangesign = lchangesign
    cassume = assume
!
!   dump (covariate,observation) pairs to plotfile on unit 15
!
    if ( dump ) call write_obscov(xx,yrs,ntot,-3e33,cov2,xyear,year,offset,lchangesign)
!
!   first-guess estaimtes of the parameters
!
    b = sd*sqrt(6.)/(4*atan(1.))
    a = mean - 0.57721*b
    if ( assume == 'scale' .and. lchangesign ) then
        xi = -0.1
    else
        xi = 0.1
    end if
    if ( lnone ) then
        alpha = 3e33
        beta = 3e33
        call fit0gevcov(a,b,xi,iter)
    else if ( assume == 'shift' .or. assume == 'scale' ) then
        beta = 3e33
        call fit1gevcov(a,b,xi,alpha,dalpha,iter)
    else if ( assume == 'both' ) then
        beta = alpha
        dbeta = dalpha
        call fit2gevcov(a,b,xi,alpha,beta,dalpha,dbeta,iter)
    else
        write(0,*) 'fitgevcov: error: unknown value for assume ',assume
    end if
    if ( assume == 'scale' ) ba = b/a
    call getreturnlevels(a,b,xi,alpha,beta,cov1,cov2,cov3,gevcovreturnlevel, &
        j1,j2,assume,t)
    if ( xyear < 1e33 ) then
        call getreturnyears(a,b,xi,alpha,beta,xyear,cov1,cov2,cov3          &
            ,gevcovreturnyear,j1,j2,tx,lchangesign,lwrite)
    endif
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
        endif
        a3(1) = a
        a3(2:3) = 3e33
        b3(1) = b
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
        return
    endif
    if ( lprint .and. .not.lweb ) print '(a,i6,a)','# doing a ',nmc     &
          ,'-member bootstrap to obtain error estimates'
    iens = 0
    scross = 0
    do iiens=1,nmc
        iens = iens + 1
        call keepalive1('Bootstrapping',iiens,nmc)
        if ( lprint .and. .not.lweb .and. mod(iiens,100) == 0 )         &
             print '(a,i6)','# ',iiens
        method = 'new'
        if ( method == 'old' ) then
            do i=1,ntot
                call random_number(ranf)
                j = 1+int(ntot*ranf)
                if ( j < 1 .or. j > ntot ) then
                    write(0,*) 'fitgev: error: j = ',j
                    call exit(-1)
                endif
                data(:,i) = xx(:,j)
            enddo
        else
            call sample_bootstrap(yrseries,yrcovariate,                 &
                 npernew,j1,j2,fyr,lyr,nens1,nens2,crosscorr,           &
                 ndecor,data,nmax,ntot,sdecor,lwrite)
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
            call fit0gevcov(aa(iens),bb(iens),xixi(iens),iter)
        else if ( assume == 'shift' .or. assume == 'scale' ) then
            betabeta(iens) = 3e33
            call fit1gevcov(aa(iens),bb(iens),xixi(iens),alphaalpha(iens),dalpha,iter)
        else if ( assume == 'both' ) then
            betabeta(iens) = beta
            call fit2gevcov(aa(iens),bb(iens),xixi(iens),               &
                 alphaalpha(iens),betabeta(iens),dalpha,dbeta,iter)
        else
            write(0,*) 'fitgevcov: error: unknown value for assume ',assume
        end if
        if ( assume == 'scale' ) baba(iens) = bb(iens)/aa(iens)
        if ( lwrite ) print *,'a,b,xi,alpha = ',aa(iens),bb(iens),xixi(iens),alphaalpha(iens)
        call getabfromcov(aa(iens),bb(iens),alphaalpha(iens),betabeta(iens),cov1,aaa,bbb)
        aacov(iens,1) = aaa
        call getabfromcov(aa(iens),bb(iens),alphaalpha(iens),betabeta(iens),cov2,aaa,bbb)
        aacov(iens,2) = aaa
        if ( cov3 < 1e33 ) then
            call getabfromcov(aa(iens),bb(iens),alphaalpha(iens),betabeta(iens),cov3,aaa,bbb)
            aacov(iens,3) = aaa
        end if
        call getreturnlevels(aa(iens),bb(iens),xixi(iens),              &
            alphaalpha(iens),betabeta(iens),cov1,cov2,cov3,gevcovreturnlevel, &
            j1,j2,assume,ttt)
        do i=1,10
            do j=1,4
                tt(iens,i,j) = ttt(i,j)
            end do
        end do
        if ( xyear < 1e33 ) then
            call getreturnyears(aa(iens),bb(iens),xixi(iens),           &
                 alphaalpha(iens),betabeta(iens),xyear,cov1,cov2,cov3,  &
                 gevcovreturnyear,j1,j2,txtxtx,lchangesign,lwrite)
            do j=1,4
                txtx(iens,j) = txtxtx(j)
            end do
            ! if the event would have been impossible in the current climate
            ! disregard this bootstrap sample
            if ( txtxtx(2) > 1e19 ) then
                if ( lwrite ) print*,'disregarding bootstrap sample'
                iens = iens - 1
            else if ( txtxtx(3) > 1e33 ) then
                write(0,*) 'fitgevcov: something went wrong',txtxtx
                iens = iens - 1
            end if
            if ( txtxtx(2) < 0.5 ) then
                write(0,*) 'fitgevcov: suspicious point ',txtxtx
            end if
        endif
    enddo
    if ( mens > mens1 ) call print_spatial_scale(scross/(nmc*ndecor),ntot/real(mens-mens1+1)/real(ndecor))
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
                do jens = 1,nmc
                    if ( tt(jens,j,i) < 1e33 ) then
                        tt(jens,j,i) = -tt(jens,j,i)
                    end if
                end do
            end do
        end do        
    endif
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
    if ( .not. lnone ) then
        call getcut( alpha25,plo,iens,alphaalpha)
        call getcut(alpha975,phi,iens,alphaalpha)
        if ( assume == 'both' ) then
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
            endif
            call getcut( t25(i,j),plo,iens,tt(1,i,j))
            call getcut(t975(i,j),phi,iens,tt(1,i,j))
        enddo
    end do
    do j=1,nj
        if ( xyear < 1e33 ) then
            call getcut( tx25(j),plo,iens,txtx(1,j))
            call getcut(tx975(j),phi,iens,txtx(1,j))
        endif
    end do
    if ( xyear < 1e33 ) then
        if ( lchangesign ) xyear = -xyear
    end if
!
!       output
!
    if ( .not.lprint ) then
        call copyab3etc(a3,b3,xi3,alpha3,beta3,t3,tx3,                  &
             a,a25,a975,b,b975,xi,xi,xi,alpha,alpha25,alpha975,         &
             beta,beta25,beta975,t,t25,t975,tx,tx25,tx975)
        if ( .not.lwrite ) return
    end if
    if ( lweb ) then
        if ( lnone ) then
            print '(a)','# <tr><td colspan="4">Fitted to GEV '//            &
                 'distribution P(x) = exp(-(1+&xi;(x-&mu;)'//             &
                     '/&sigma;)^(-1/&xi;))</td></tr>'
            call printab(restrain,lnone,lweb)
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'//         &
                 '&mu;:</td><td>',a,'</td><td>',       &
                 a25,'...',a975,'</td></tr>'
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'//         &
                 '&sigma;:</td><td>',b,'</td><td>',    &
                 b25,'...',b975,'</td></tr>'
            if ( assume == 'scale' ) then
                print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'// &
                 '&sigma;/&mu;:</td><td>',ba,'</td><td>',ba25,'...',ba975,'</td></tr>'
            end if
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'//    &
                 '&xi;:</td><td>',xi,'</td><td>',xi25,'...',xi975,          &
                 '</td></tr>'
        else
            print '(a)','# <tr><td colspan="4">Fitted to GEV '//            &
                 'distribution P(x) = exp(-(1+&xi;(x-&mu;'')'//             &
                     '/&sigma;'')^(-1/&xi;))</td></tr>'
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
             &           '&sigma;/&mu;:</td><td>',ba,'</td><td>',ba25,'...',ba975,'</td></tr>'
            end if
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'//    &
                 '&xi;:</td><td>',xi,'</td><td>',xi25,'...',xi975,          &
                 '</td></tr>'
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'//    &
                 '&alpha;:</td><td>',alpha,'</td><td>',alpha25,'...',       &
                 alpha975,'</td></tr>'
            if ( assume == 'both' ) then
                print '(a,f16.3,a,f16.3,a,f16.3,a)',                        &
                     '# <tr><td colspan=2>&beta;:</td><td>',beta,           &
                     '</td><td>',beta25,'...',beta975,'</td></tr>'
            end if
        end if
    else
        print '(a,i5,a)','# Fitted to GEV distribution in ',iter,' iterations'
        if ( lnone ) then
            print '(a)','# P(x) = exp(-(1+xi*(x-a/b)**(-1/xi))'
            call printab(restrain,lnone,lweb)
            print '(a,f16.3,a,f16.3,a,f16.3)','# a = ',a,' \\pm ',(a975-a25)/2
            print '(a,f16.3,a,f16.3,a,f16.3)','# b = ',b,' \\pm ',(b975-b25)/2
        else
            print '(a)','# P(x) = exp(-(1+xi*(x-a''/b'')**(-1/xi)) with'
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
        end if
        if ( assume == 'scale' ) print '(a,f16.3,a,f16.3,a,f16.3)', &
            '# b/a  = ',ba,' \\pm ',(ba975-ba25)/2
        print '(a,f16.3,a,f16.3,a,f16.3)','# xi  = ',xi,' \\pm ',(xi975-xi25)/2
        if ( .not. lnone ) then
            print '(a,f16.3,a,f16.3,a,f16.3)','# alpha ',alpha,' \\pm ',(alpha975-alpha25)/2
            if ( assume == 'both' ) then
                print '(a,f16.3,a,f16.3,a,f16.3)','# beta  ',beta,' \\pm ',(beta975-beta25)/2
            end if
        end if
    end if
    call printcovreturntime(year,xyear,idmax,tx,tx25,tx975,yr1a,yr2a,yr2b,lweb,plot,assume,lnone,1)
    if ( .not. lnone ) call printcovpvalue(txtx,nmc,iens,lweb,plot)
    call printcovreturnvalue(ntype,t,t25,t975,yr1a,yr2a,lweb,plot,assume,lnone)
    if ( .not. lnone ) call printcovreturntime(year,xyear,idmax,tx,tx25,tx975,yr1a,yr2a,yr2b,lweb,.false.,assume,lnone,2)

    if ( dump ) then
        call plot_tx_cdfs(txtx,nmc,iens,ntype,j1,j2)
    end if
    if ( plot ) write(11,'(3g20.4,a)') alpha,alpha25,alpha975,' alpha'
    if ( dump ) call write_dthreshold(cov1,cov2,cov3,acov,offset,lchangesign)

    ! no cuts
    mindata = -2e33
    minindx = -2e33
    pmindata = -1
    snorm = 1
    frac = 1
    ! GEV fit
    nfit = 5

    if ( lchangesign ) then
        ! we had flipped the sign on a,b,alpha but not yet on the rest, flip back...
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
        call plot_ordered_points(yy,ys,yrs,ntot,ntype,nfit,                 &
             frac,a,b,xi,j1,j2,minindx,mindata,pmindata,                &
             year,xyear,snorm,lchangesign,lwrite,.true.)
    else
    ! compute distribution at past year and plot it
        call adjustyy(ntot,xx,assume,a,b,alpha,beta,cov1,yy,zz,aaa,bbb,lchangesign,lwrite)
        ys(1:ntot) = yy(1:ntot)
        print '(a,i5)','# distribution in year ',yr1a
        call plotreturnvalue(ntype,t25(1,1),t975(1,1),j2-j1+1)
        call plot_ordered_points(yy,ys,yrs,ntot,ntype,nfit,                 &
             frac,aaa,bbb,xi,j1,j2,minindx,mindata,pmindata,                &
             year,xyear,snorm,lchangesign,lwrite,.false.)

        ! compute distribution at present year and plot it
        call adjustyy(ntot,xx,assume,a,b,alpha,beta,cov2,                   &
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
        call plot_ordered_points(yy,ys,yrs,ntot,ntype,nfit,                 &
             frac,aaa,bbb,xi,j1,j2,minindx,mindata,pmindata,                &
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

end subroutine fitgevcov

subroutine fit0gevcov(a,b,xi,iter)
    use AmoebaToGSL
    implicit none
    integer iter
    real a,b,xi
    integer i
    real q(5),p(4,3),y(4),tol
    logical lok
    character cassume*5
    common /fitdata4/ cassume
    real llgevcov
    external llgevcov
!
10  continue
    q(1) = a
    q(2) = b
    q(3) = xi
    q(4) = 0
    q(5) = 0
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
    lok = .false.
    do i=1,4
        q(1) = p(i,1)
        q(2) = p(i,2)
        q(3) = p(i,3)
        q(4) = 0
        q(5) = 0
        y(i) = llgevcov(q)
        if ( y(i) < 1e33 ) lok = .true.
    enddo
    if ( .not.lok ) then
        if ( xi /= 0 ) then
            xi = xi + 0.1*xi/abs(xi)
        else if ( cassume == 'scale' .and. a < 0 ) then
            xi = -0.1
        else 
            xi = 0.1
        end if
        if ( abs(xi) < 1 ) then
            goto 10
        else
            write(0,*) 'fit1gevcov: error: cannot find initial fit values',a,b,xi
            a = 3e33
            b = 3e33
            xi = 3e33
        end if
    end if
    tol = 1e-4
    call amoeba(p,y,4,3,3,tol,llgevcov,iter)
!       maybe add restart later
    a = p(1,1)
    b = abs(p(1,2))
    xi = p(1,3)
end subroutine fit0gevcov

subroutine fit1gevcov(a,b,xi,alpha,dalpha,iter)
    use AmoebaToGSL
    implicit none
    integer iter
    real a,b,xi,alpha,dalpha
    integer i
    real q(5),p(5,4),y(5),tol
    logical lok
    character cassume*5
    common /fitdata4/ cassume
    real llgevcov
    external llgevcov
!
10  continue
    q(1) = a
    q(2) = b
    q(3) = xi
    q(4) = alpha
    q(5) = 3e33
    p(1,1) = q(1) *0.9
    p(1,2) = q(2) *0.9
    p(1,3) = q(3) *0.9
    p(1,4) = q(4) - dalpha
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
    p(4,3) = p(1,3) *1.2 + 0.1
    p(4,4) = p(1,4)
    p(5,1) = p(1,1)
    p(5,2) = p(1,2)
    p(5,3) = p(1,3)
    p(5,4) = p(1,4) + 2*dalpha
    lok = .false.
    do i=1,5
        q(1) = p(i,1)
        q(2) = p(i,2)
        q(3) = p(i,3)
        q(4) = p(i,4)
        y(i) = llgevcov(q)
        if ( y(i) < 1e33 ) lok = .true.
    enddo
    if ( .not.lok ) then
        if ( xi /= 0 ) then
            xi = xi + 0.1*xi/abs(xi)
        else if ( cassume == 'scale' .and. a < 0 ) then
            xi = -0.1
        else 
            xi = 0.1
        end if
        if ( abs(xi) < 1 ) then
            goto 10
        else
            write(0,*) 'fit1gevcov: error: cannot find initial fit values',a,b,xi,alpha
            a = 3e33
            b = 3e33
            xi = 3e33
            alpha = 3e33
        end if
    end if
    tol = 1e-4
    call amoeba(p,y,5,4,4,tol,llgevcov,iter)
!       maybe add restart later
    a = p(1,1)
    b = abs(p(1,2))
    xi = p(1,3)
    alpha = p(1,4)
end subroutine fit1gevcov

subroutine fit2gevcov(a,b,xi,alpha,beta,dalpha,dbeta,iter)
    use AmoebaToGSL
    implicit none
    integer iter
    real a,b,xi,alpha,beta,dalpha,dbeta
    integer i
    real q(5),p(6,5),y(6),tol
    logical lok
    character cassume*5
    common /fitdata4/ cassume
    real llgevcov
    external llgevcov
!
10  continue
    q(1) = a
    q(2) = b
    q(3) = xi
    q(4) = alpha
    q(5) = beta
    p(1,1) = q(1) *0.9
    p(1,2) = q(2) *0.9
    p(1,3) = q(3) *0.9
    p(1,4) = q(4) - dalpha
    p(1,5) = q(5) - dbeta
    p(2,1) = p(1,1) *1.2
    p(2,2) = p(1,2)
    p(2,3) = p(1,3)
    p(2,4) = p(1,4)
    p(2,5) = p(1,5)
    p(3,1) = p(1,1)
    p(3,2) = p(1,2) *1.2
    p(3,3) = p(1,3)
    p(3,4) = p(1,4)
    p(3,5) = p(1,5)
    p(4,1) = p(1,1)
    p(4,2) = p(1,2)
    p(4,3) = p(1,3) *1.2 + 0.1
    p(4,4) = p(1,4)
    p(4,5) = p(1,5)
    p(5,1) = p(1,1)
    p(5,2) = p(1,2)
    p(5,3) = p(1,3)
    p(5,4) = p(1,4) + 2*dalpha
    p(5,5) = p(1,5)
    p(6,1) = p(1,1)
    p(6,2) = p(1,2)
    p(6,3) = p(1,3)
    p(6,4) = p(1,4)
    p(6,5) = p(1,5) + 2*dbeta
    lok = .false.
    do i=1,6
        q(1) = p(i,1)
        q(2) = p(i,2)
        q(3) = p(i,3)
        q(4) = p(i,4)
        q(5) = p(i,5)
        y(i) = llgevcov(q)
        if ( y(i) < 1e33 ) lok = .true.
    enddo
    if ( .not.lok ) then
        if ( xi /= 0 ) then
            xi = xi + 0.1*xi/abs(xi)
        else if ( cassume == 'scale' .and. a < 0 ) then
            xi = -0.1
        else 
            xi = 0.1
        end if
        if ( abs(xi) < 1 ) then
            goto 10
        else
            write(0,*) 'fit1gevcov: error: cannot find initial fit values'
            a = 3e33
            b = 3e33
            xi = 3e33
            alpha = 3e33
            beta = 3e33
        end if
    end if
    tol = 1e-4
    call amoeba(p,y,6,5,5,tol,llgevcov,iter)
!       maybe add restart later
    a = p(1,1)
    b = abs(p(1,2))
    xi = p(1,3)
    alpha = p(1,4)
    beta = p(1,5)
end subroutine fit2gevcov

real function llgevcov(p)
!
!   computes the log-likelihood function for a covariant-dependent GEV distribution
!   with parameters a,b,xi,alpha,beta=p(1-5) and data in common.
!
    implicit none
!
    real p(5)
!
    integer i
    real x,z,xi,s,aa,bb
!
    integer nmax,ncur
    parameter(nmax=1000000)
    real data(2,nmax),restrain
    logical llwrite,llchangesign
    common /fitdata3/ data
    common /fitdata2/ restrain,ncur,llwrite,llchangesign
    character cassume*5
    common /fitdata4/ cassume
    integer,save :: init=0
!
    llgevcov = 0
    if ( abs(p(3)) > 10 ) then
        llgevcov = 3e33
        goto 999
    endif
    if ( restrain < 0 ) then
        write(0,*) 'llgevcov: restrain<0 ',restrain
        call exit(-1)
    end if
    if ( cassume == 'scale' .and. llchangesign ) then
        if ( init == 0 ) then
            init = 1
            write(0,*) 'Enforcing a hard lower bound of zero.'
        end if
        if ( p(3) >= 0 ) then
            ! scaling implies (for climate) that the distribution cannot cross zero, 
            ! so xi must be < 0
            llgevcov = 3e33
            goto 999
        end if
        if ( p(2) > p(1)*p(3) ) then
            ! scaling implies (for climate) that the distribution cannot cross zero, 
            ! so the upper limit must be <0 (for flipped signs)
            llgevcov = 3e33
            goto 999
        end if
    end if
    do i=1,ncur
        call getabfromcov(p(1),p(2),p(4),p(5),data(2,i),aa,bb)
        if ( abs(bb) < 1e-30 ) then
            llgevcov = 3e33
            goto 999
        end if
        z = (data(1,i)-aa)/bb
        xi = p(3)
        if ( abs(xi) < 1e-4 ) then
            if ( -z+xi*z**2/2 > log(3e33) ) then
                llgevcov = 3e33
                goto 999
            end if
            s = - exp(-z+xi*z**2/2) - z*(1+xi-xi*z/2)
        else
            if ( 1+xi*z <= 0 ) then
!!!                 write(0,*) 'GEV undefined',(1+xi*z)
                llgevcov = 3e33
                goto 999
            else if ( -log(1+xi*z)/xi > log(3e33) ) then
                ! too large...
                llgevcov = 3e33
                goto 999
            else
                s = - (1+1/xi)*log(1+xi*z) - (1+xi*z)**(-1/xi)
            endif
        endif
        s = s - log(abs(bb))
        llgevcov = llgevcov + s            
        if ( .false. .and. llwrite ) then
            print *,i,data(1,i),aa,bb,xi,z,-s,-llgevcov
        end if
    enddo
!   normalization is not 1 in case of cut-offs
    call gevcovnorm(p(1),p(2),p(3),p(4),p(5),s)
    if ( s < 1e33 ) then
        llgevcov = llgevcov - ncur*log(s)
    else
        llgevcov = 3e33
        goto 999
    end if
    if ( restrain /= 0 ) then
!           preconditioning on xi with gaussian of width restrain/2
!           around 0
        llgevcov = llgevcov - (xi/(restrain/2))**2/2
    endif
!       minimum, not maximum
    llgevcov = -llgevcov
!
999 continue
    if ( llwrite .and. .false. ) print *,'a,b,xi,alpha,llgevcov = ',p(1),p(2),p(3),p(4),llgevcov
end function llgevcov

subroutine gevcovnorm(a,b,xi,alpha,beta,s)
    implicit none
	include 'getopts.inc'
    real a,b,xi,alpha,beta,s
    real z1,z2

    if ( minindx > -1e33 .or. maxindx < 1e33 ) then
        write(0,*) 'gevcovnorm: boundaries not yet avaiable for fit of GEV(t)'
        call exit(-1)
    else
        s = 1
    endif
!!!        print *,'gevcovnorm: norm = ',a,b,s
end subroutine gevcovnorm

real function gevcovreturnlevel(a,b,xi,alpha,beta,x,cov)
!
!       compute return times given the GEV distribution parameters a,b,xi and 
!       x = log10(returntime) for covariant cov and fit parameter alpha
!       Uses a few Taylor series approximation for xi small and/or return time large
!
    implicit none
    real a,b,xi,alpha,beta,x,cov
    real aa,bb,y,t
    call getabfromcov(a,b,alpha,beta,cov,aa,bb)
    if ( abs(xi) > 10 ) then
        gevcovreturnlevel = 3e33
    else if ( abs(xi) < 1e-4 ) then
        if ( x <= 8 ) then
            y = log(-log(1-dble(10)**(dble(-x))))
        else
            y = -x*log(10.)
        end if
        t = aa - bb*y + bb*xi/2*y**2
    else
        if ( x <= 8 ) then
            t = aa - bb/xi*(1-(-log(1-dble(10)**(dble(-x))))**(-xi))
        else
            t = aa - bb/xi*(1-10.**(x*xi))
        end if
    end if
    gevcovreturnlevel = t
end function gevcovreturnlevel

real function gevcovreturnyear(a,b,xi,alpha,beta,xyear,cov,lchangesign)
!
!   compute the return time of the value xyear with the fitted values
!
    implicit none
    real a,b,xi,alpha,beta,xyear,cov
    logical lchangesign
    real x,y,z,tx,arg,aa,bb

    x = xyear
    call getabfromcov(a,b,alpha,beta,cov,aa,bb)
    z = (1 + xi*(x-aa)/bb)
    if ( z < 0 ) then
        y = 1e20
    else if ( abs(xi) > 1e-3 ) then
        !!!y = -z**(-1/xi) sometimes does not fit
        y = log(z)*(-1/xi)
        if ( y < 20*log(10.) ) then
            y = -exp(y)
        else
            y = -1e20
        end if
    else
        if ( xi == 0 ) then
            arg = -(x-aa)/bb
        else
            arg = -(x-aa)/bb + xi/2*((x-aa)/bb)**2
        end if
        if ( arg > log(3e33) ) then
            y = -1e20
        else
            y = -exp(arg)
        end if
    end if
    if ( y > 1e19 ) then
        tx = 1e20 ! infinity, not undefined!
    else if ( y > log(3e33) ) then
        tx = 0
    else if ( abs(y) > 1e-3 ) then
        tx = 1/(1 - exp(y))
    else if ( abs(y) > 1e-19 ) then
        tx = -1/(y + 0.5*y**2)
    else
        tx = 1e20
    end if
    if ( .false. .and. tx > 1e20 ) then
        write(0,*) 'gevcovreturnyear: tx > 1e20: ',tx
        write(0,*) 'a,b,xi,alpha,xyear = ',a,b,alpha,xi,xyear
    endif
    gevcovreturnyear = tx
end function gevcovreturnyear
