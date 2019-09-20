subroutine fitgaucov(yrseries,yrcovariate,npernew,fyr,lyr &
            ,mens1,mens,crosscorr,a3,b3,alpha3,beta3,j1,j2,nens1,nens2 &
            ,lweb,ntype,lchangesign,yr1a,yr2a,yr2b,xyear,idmax,cov1,cov2 &
            ,cov3,offset,t3,tx3,assume,confidenceinterval,ndecor,lboot &
            ,lprint,dump,plot,lwrite)
!
!   a fit a gaussian distribution with mean linearly dependent on a covariate 
!   to the data
!
    implicit none
!
    integer npernew,fyr,lyr,mens1,mens,ntot,ntype,j1,j2,nens1,nens2,yr1a,yr2a,yr2b,ndecor
    real yrseries(npernew,fyr:lyr,0:mens), &
            yrcovariate(npernew,fyr:lyr,0:mens),crosscorr(0:mens,0:mens), &
            a3(3),b3(3),alpha3(3),beta3(3),xyear,cov1,cov2,cov3, &
            offset,t3(3,10,3),tx3(3,3),confidenceinterval
    logical lweb,lchangesign,lboot,lprint,dump,plot,lwrite
    character assume*(*),idmax*(*)
!
    integer nmc
    parameter(nmc=1000)
    integer i,j,jj,n,nx,iter,iens,nfit,imc,ier,year,nj
    integer,allocatable :: yrs(:)
    real a,b,ba,t(10,4),t25(10,4),t975(10,4),tx(4),tx25(4),tx975(4), &
            aa(nmc),bb(nmc),baba(nmc),tt(nmc,10,4),xi,alpha,beta,dalpha,dbeta, &
            xmin,z,x,f,txtx(nmc,4),alphaalpha(nmc),betabeta(nmc), &
            mean,sd,ranf,mindata,minindx,pmindata,snorm,s,frac,scross,sdecor
    real a25,a975,b25,b975,ba25,ba975,alpha25,alpha975,beta25,beta975,aa25,aa975,bb25,bb975
    real adev,var,skew,curt,aaa,bbb,siga,chi2,q,cmin,cmax,plo,phi
    real ttt(10,4),txtxtx(4),dum,xi3(3),acov(3,3),aacov(nmc,3)
    real,allocatable :: xx(:,:),yy(:),ys(:),zz(:),sig(:)
    logical last
    character lgt*4,method*3
    external gaucovreturnlevel,gaucovreturnyear
!
    integer nmax
    parameter(nmax=1000000)
    integer ncur
    real data(2,nmax),restrain
    logical llwrite,llchangesign,lnone
    common /fitdata3/ data
    common /fitdata2/ restrain,ncur,llwrite,llchangesign
    character cassume*5
    common /fitdata4/ cassume
!
    year = yr2a
    allocate(yrs(0:nmax))
    allocate(xx(2,nmax))
    if ( cov1 == 0 .and. cov2 == 0 ) then
        lnone = .true.
    else
        lnone = .false.
    end if
    restrain = 0 ! no shape parameter

    if ( lwrite ) print *,'fitgaucov: calling fill_linear_array'
    call fill_linear_array(yrseries,yrcovariate,npernew,j1,j2, &
            fyr,lyr,mens1,mens,xx,yrs,nmax,ntot,lwrite)
    if ( lprint .and. lweb ) then
        print '(a,i9,a)','# <tr><td>N:</td><td>&nbsp;</td><td>', &
                ntot,'</td><td>&nbsp;</td></tr>'
    end if
    if ( ntot < 5 ) then
        print '(a)','</table>'
        print '(a)','There are not enough points to fit a normal distribution. Please use a (much) longer time series.'
        print '(a)','<p>'
        call exit(-1)
    end if

    if ( lwrite ) then
        print *,'fitgaucov: input'
        print *,'year,xyear  = ',year,xyear
        print *,'cov1,cov2,offset ',cov1,cov2,offset
        if ( cov3 < 1e33 ) then
            print *,'cov3           = ',cov3
        end if
    end if
    if ( ntot.eq.0 ) then
        if ( lwrite ) print *,'fitgaucov: ntot=0'
        a3 = 3e33
        b3 = 3e33
        xi3 = 0
        alpha3 = 3e33
        beta3 = 3e33
        t3 = 3e33
        tx3 = 3e33
        return
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
    if ( dump ) call write_obscov(xx,yrs,ntot,-3e33,cov2,xyear,year,offset,lchangesign)
!
    sig = 0
    call moment(yy,ntot,mean,adev,sd,var,skew,curt)
    call fit(zz,yy,ntot,sig,0,aaa,alpha,siga,dalpha,chi2,q)
    if ( lwrite ) then
        print *,'fitgaucov: computed initialisation values:'
        print *,'mean,sd,alpha,dalpha = ',mean,sd,alpha,dalpha
    end if
!
!   a trivial case which causes no end of trouble
!
    if ( sd.eq.0 ) then
        if ( lwrite ) print *,'fitgaucov: sd=0, everything undfined'
        a3 = 3e33
        b3 = 3e33
        xi3 = 0
        alpha3 = 3e33
        beta3 = 3e33
        t3 = 3e33
        tx3 = 3e33
        return
    endif
!
!   copy to common for routine llgausscov
!
    ncur = ntot
    do i=1,ncur
        data(:,i) = xx(:,i)
    enddo
    cassume = assume
!   
!   fit, using Numerical Recipes routines
!
    a = mean
    b = sd
    if ( lnone ) then
        alpha = 3e33
        beta = 3e33
        call fit0gaucov(a,b,iter)
    else if ( assume.eq.'shift' .or. assume.eq.'scale' ) then
        beta = 3e33
        call fit1gaucov(a,b,alpha,dalpha,iter)
    else if ( assume.eq.'both' ) then
        call fit2gaucov(a,b,alpha,beta,dalpha,dbeta,iter)
    else
        write(0,*) 'fitgaucov: cannot handle assume = ',assume
        call exit(-1)
    end if
    if ( assume == 'scale' ) ba = b/a
    dum = 0
    call getreturnlevels(a,b,dum,alpha,beta,cov1,cov2,cov3,gaucovreturnlevel,j1,j1,assume,t)
    if ( xyear.lt.1e33 ) then
        call getreturnyears(a,b,dum,alpha,beta,xyear,cov1,cov2,cov3 &
            ,gaucovreturnyear,j1,j1,tx,lchangesign,lwrite)
    endif
    call getabfromcov(a,b,alpha,beta,cov1,aaa,bbb)
    acov(1,1) = aaa
    call getabfromcov(a,b,alpha,beta,cov2,aaa,bbb)
    acov(1,2) = aaa
    if ( cov3 < 1e33 ) then
        call getabfromcov(a,b,alpha,beta,cov3,aaa,bbb)
        acov(1,3) = aaa
    end if
    if ( dump ) call write_threshold(cmin,cmax,a,b,xi,alpha,beta,offset,lchangesign,gaucovreturnlevel)
!
!   Bootstrap for error estimate
!
    if ( .not.lboot ) then
        if ( lchangesign ) then
            a = -a
            acov = -acov
            aacov = -aacov
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
        xi3 = 0
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
    if ( lprint .and. .not.lweb ) print '(a,i6,a)','# doing a ',nmc &
             ,'-member bootstrap to obtain error estimates'
    scross = 0
    do iens=1,nmc
        if ( lprint .and. .not.lweb .and. mod(iens,100).eq.0 ) print '(a,i6)','# ',iens
        method = 'new'
        if ( method == 'old' ) then
            n = 1 + (ntot-1)/ndecor
            do i=1,n
                ! we do not have the information here to check whether the
                ! data points were contiguous in the original series...
                ! TODO: propagate that information              
                call random_number(ranf)
                j = 1 + min(ntot-ndecor,int((ntot-ndecor)*ranf))
                if ( j.lt.1 .or. j.gt.ntot ) then
                    write(0,*) 'fitgaucov: error: j = ',j
                    call exit(-1)
                endif
                if ( i.lt.n ) then ! the blocks that fit in whole
                    do jj=0,ndecor-1
                        data(:,1+(i-1)*ndecor+jj) = xx(:,j+jj)
                    end do
                else
                    do jj=0,ndecor-1 ! one more block to the end, the previous block is shortened
                        data(:,1+ntot-ndecor+jj) = xx(:,j+jj)
                    end do
                end if
            enddo
        else
            call sample_bootstrap(yrseries,yrcovariate, &
                    npernew,j1,j2,fyr,lyr,nens1,nens2,crosscorr, &
                    ndecor,data,nmax,ntot,sdecor,lwrite)
            scross = scross + sdecor
        end if
        aa(iens) = a
        bb(iens) = b
        alphaalpha(iens) = alpha
        llwrite = .false.
        if ( lnone ) then
            alphaalpha(iens) = 3e33
            betabeta(iens) = 3e33
            call fit0gaucov(aa(iens),bb(iens),iter)
        else if ( assume.eq.'shift' .or. assume.eq.'scale' ) then
            betabeta(iens) = 3e33
            call fit1gaucov(aa(iens),bb(iens),alphaalpha(iens),dalpha,iter)
        else if ( assume.eq.'both' ) then
            betabeta(iens) = beta
            call fit2gaucov(aa(iens),bb(iens),alphaalpha(iens),betabeta(iens),dalpha,dbeta,iter)
        else
            write(0,*) 'fitgaucov: cannot handle assume = ',assume
            call exit(-1)
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
        call getreturnlevels(aa(iens),bb(iens),dum, &
                alphaalpha(iens),betabeta(iens), &
                cov1,cov2,cov3,gaucovreturnlevel,j1,j1,assume,ttt)
        do i=1,10
            do j=1,4
                tt(iens,i,j) = ttt(i,j)
            end do
        end do
        if ( xyear.lt.1e33 ) then
            call getreturnyears(aa(iens),bb(iens),dum, &
                    alphaalpha(iens),betabeta(iens),xyear,cov1,cov2,cov3, &
                    gaucovreturnyear,j1,j1,txtxtx,lchangesign,lwrite)
            do j=1,4
                txtx(iens,j) = txtxtx(j)
            end do
        endif
    enddo
    if ( j1 /= j2 ) then
        write(0,'(2a,i3,a)') 'The return times are for a given month, for any month ', &
            'in the season they are a factor ',j2-j1+1,' lower.'
    end if
    if ( mens > mens1 ) call print_spatial_scale(scross/(nmc*ndecor),ntot/real(mens-mens1+1)/real(ndecor))
    iens = nmc
    if ( lchangesign ) then
        a = -a
        acov = -acov
        aa = -aa
        aacov = -aacov
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
    endif
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
    call getcut( alpha25,plo,nmc,alphaalpha)
    call getcut(alpha975,phi,nmc,alphaalpha)
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
            call getcut( t25(i,j),plo,nmc,tt(1,i,j))
            call getcut(t975(i,j),phi,nmc,tt(1,i,j))
        enddo
    end do
    do j=1,nj
        if ( xyear.lt.1e33 ) then
            call getcut( tx25(j),plo,nmc,txtx(1,j))
            call getcut(tx975(j),phi,nmc,txtx(1,j))
        endif
    end do
    if ( xyear < 1e33 ) then
        if ( lchangesign ) xyear = -xyear
    end if
!
!   output
!
    if ( .not.lprint ) then
        xi = 0
        call copyab3etc(a3,b3,xi3,alpha3,beta3,t3,tx3, &
                a,a25,a975,b,b975,xi,xi,xi,alpha,alpha25,alpha975, &
                beta,beta25,beta975,t,t25,t975,tx,tx25,tx975)
        if ( .not.lwrite ) return
    end if
    if ( lweb ) then
        if ( lnone ) then
            print '(a)','# <tr><td colspan="4">Fitted to normal '// &
                'distribution P(x) = exp(-(x-&mu;)&sup2;'// &
                '/(2&sigma;&sup2;))/(&sigma;&radic;(2&pi;))</td></tr>'
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'// &
                '&mu;:</td><td>',a,'</td><td>',a25,'...',a975,'</td></tr>'
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'// &
                '&sigma;:</td><td>',b,'</td><td>',b25,'...',b975,'</td></tr>'
            if ( assume == 'scale' ) then
                print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'// &
                    '&sigma;/&mu;:</td><td>',ba,'</td><td>',ba25,'...',ba975,'</td></tr>'
            end if
        else
            print '(a)','# <tr><td colspan="4">Fitted to normal '// &
                'distribution P(x) = exp(-(x-&mu;'')&sup2;'// &
                '/(2&sigma;''&sup2;))/(&sigma;''&radic;(2&pi;))</td></tr>'
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
                    '&alpha;:</td><td>',alpha,'</td><td>',alpha25,'...', &
                    alpha975,'</td></tr>'
            if ( assume.eq.'both' ) then
                print '(a,f16.3,a,f16.3,a,f16.3,a)', &
                        '# <tr><td colspan=2>&beta;:</td><td>',beta, &
                        '</td><td>',beta25,'...',beta975,'</td></tr>'
            end if
        end if
    else
        print '(a,i5,a)','# Fitted to Gaussian distribution in ',iter,' iterations'
        if ( lnone ) then
            print '(a)','# p(x) = exp(-(x-a)^2/(2*b^2))/(b''*sqrt(2*pi))'
            print '(a,f16.3,a,f16.3)','# a = ',a,' \\pm ',(a975-a25)/2
            print '(a,f16.3,a,f16.3)','# b = ',b,' \\pm ',(b975-b25)/2
            if ( assume == 'scale' ) print '(a,f16.3,a,f16.3,a,f16.3)', &
                '# b/a = ',ba,' \\pm ',(ba975-ba25)/2
        else
            print '(a)','# p(x) = exp(-(x-a'')^2/(2*b''^2))/(b''*sqrt(2*pi)) with'
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
                '# b/a = ',ba,' \\pm ',(ba975-ba25)/2
            print '(a,f16.3,a,f16.3)','# alpha = ',alpha,' \\pm ',(alpha975-alpha25)/2
            if ( assume.eq.'both' ) then
                print '(a,f16.3,a,f16.3)','# beta  ',beta,' \\pm ',(beta975-beta25)/2
            end if
        end if
    endif
    call printcovreturntime(year,xyear,idmax,tx,tx25,tx975,yr1a,yr2a,yr2b,lweb,plot,assume,lnone,1)
    if ( .not. lnone ) call printcovpvalue(txtx,nmc,nmc,lweb,plot)
    call printcovreturnvalue(ntype,t,t25,t975,yr1a,yr2a,lweb,plot,assume,lnone)
    if ( .not. lnone ) call printcovreturntime(year,xyear,idmax,tx,tx25,tx975,yr1a,yr2a,yr2b,lweb,plot,assume,lnone,2)

    if ( dump ) then
        call plot_tx_cdfs(txtx,nmc,nmc,ntype,j1,j1)
    end if
    if ( plot ) write(11,'(3g20.4,a)') alpha,alpha25,alpha975,' alpha'
    if ( dump ) call write_dthreshold(cov1,cov2,cov3,acov,offset,lchangesign)

    ! no cuts
    mindata = -2e33
    minindx = -2e33
    pmindata = -1
    snorm = 1
    frac = 1
    ! fit to gauss (normal distribution)
    nfit = 2
    if ( lchangesign ) then
        a = -a
        if ( .not. lnone ) then
            alpha = -alpha
            if ( assume == 'both' ) beta = -beta
        end if
    end if

    if ( lnone ) then
        ! we give the probability per month, not per N months selection...
        call plotreturnvalue(ntype,t25(1,1),t975(1,1),1)
        ys(1:ntot) = yy(1:ntot)
        call plot_ordered_points(yy,ys,yrs,ntot,ntype,nfit,                 &
            frac,a,b,xi,j1,j1,minindx,mindata,pmindata,                &
            year,xyear,snorm,lchangesign,lwrite,.true.)
    else
        ! compute distribution at past year and plot it
        call adjustyy(ntot,xx,assume,a,b,alpha,beta,cov1,yy,zz,aaa,bbb,lchangesign,lwrite)
        ys(1:ntot) = yy(1:ntot)
        print '(a,i5)','# distribution in year ',yr1a
        call plotreturnvalue(ntype,t25(1,1),t975(1,1),1)
        call plot_ordered_points(yy,ys,yrs,ntot,ntype,nfit, &
            frac,aaa,bbb,dum,j1,j1,minindx,mindata,pmindata, &
            year,xyear,snorm,lchangesign,lwrite,.false.)

        ! compute distribution at present year and plot it
        call adjustyy(ntot,xx,assume,a,b,alpha,beta,cov2,yy,zz,aaa,bbb,lchangesign,lwrite)
        ys(1:ntot) = yy(1:ntot)
        print '(a)'
        print '(a)'
        print '(a,i5)','# distribution in year ',yr2a
        if ( cov3 < 1e33 ) then
            last = .false.
        else
            last = .true.
        end if
        call plotreturnvalue(ntype,t25(1,2),t975(1,2),1)
        call plot_ordered_points(yy,ys,yrs,ntot,ntype,nfit, &
            frac,aaa,bbb,dum,j1,j1,minindx,mindata,pmindata, &
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
end subroutine fitgaucov

subroutine fit0gaucov(a,b,iter)
    use AmoebaToGSL
    implicit none
    integer iter
    real a,b
    integer i
    real q(4),p(3,2),y(3),tol
    real llgausscov
    external llgausscov
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
        y(i) = llgausscov(q)
    enddo
    tol = 1e-4
    call amoeba(p,y,3,2,2,tol,llgausscov,iter)
!   maybe add restart later
    a = p(1,1)
    b = p(1,2)
end subroutine fit0gaucov

subroutine fit1gaucov(a,b,alpha,dalpha,iter)
    use AmoebaToGSL
    implicit none
    integer iter
    real a,b,alpha,dalpha
    integer i
    real q(4),p(4,3),y(4),tol
    real llgausscov
    external llgausscov
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
        y(i) = llgausscov(q)
    enddo
    tol = 1e-4
    call amoeba(p,y,4,3,3,tol,llgausscov,iter)
!   maybe add restart later
    a = p(1,1)
    b = p(1,2)
    alpha = p(1,3)
end subroutine fit1gaucov

subroutine fit2gaucov(a,b,alpha,beta,dalpha,dbeta,iter)
    use AmoebaToGSL
    implicit none
    integer iter
    real a,b,alpha,beta,dalpha,dbeta
    integer i
    real q(4),p(5,4),y(5),tol
    real llgausscov
    external llgausscov
!
!   fit, using Numerical Recipes routines
!   
    q(1) = a
    q(2) = b
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
        y(i) = llgausscov(q)
    enddo
    tol = 1e-4
    call amoeba(p,y,5,4,4,tol,llgausscov,iter)
!   maybe add restart later
    a = p(1,1)
    b = p(1,2)
    alpha = p(1,3)
    beta = p(1,4)
end subroutine fit2gaucov

real function llgausscov(p)
!
!   computes the log-likelihood function for a normal distribution
!   with parameters alpha,beta=p(1),p(2) and data in common.
!
    implicit none
!   
    real p(4)
!
    integer i
    real z,s,aa,bb
!
    integer nmax
    parameter(nmax=1000000)
    integer ncur
    real data(2,nmax),restrain
    logical llwrite,llchangesign
    common /fitdata3/ data
    common /fitdata2/ restrain,ncur,llwrite,llchangesign
    character cassume*5
    common /fitdata4/ cassume
!   
    llgausscov = 0
    do i=1,ncur
        call getabfromcov(p(1),p(2),p(3),p(4),data(2,i),aa,bb)
        z = (data(1,i) - aa)/bb
        llgausscov = llgausscov - z**2/2 - log(abs(bb))
    enddo
!   normalization is not 1 in case of cut-offs
    call gauscovnorm(aa,bb,s)
    llgausscov = llgausscov - ncur*log(s)
!   minimum, not maximum
    llgausscov = -llgausscov
!!!        print *,'a,b,llgausscov = ',p(1),p(2),llgausscov
!
end function llgausscov

subroutine gauscovnorm(a,b,s)
    implicit none
    include 'getopts.inc'
    real a,b,s
    real z1,z2,sqrt2
    real erfcc
    external erfcc
    if ( minindx.gt.-1e33 .or. maxindx.lt.1e33 ) then
        write(0,*) 'gauscovnorm: boundaries not yet available for fit of gauss(t)'
        call exit(-1)
    else
        s = 1
    endif
!!!        print *,'gauscovnorm: norm = ',a,b,s
end subroutine gauscovnorm

real function gaucovreturnlevel(a,b,xi,alpha,beta,x,cov)
!
!   compute return times given the normal distribution parameters a,b and 
!   x = log10(returntime) for covariant cov and fit parameter alpha
!
    implicit none
    real :: a,b,xi,alpha,beta,x,cov
    real :: aa,bb,f,z,t
    real,external :: serfci
!
    !!!print *,'# gauscovreturnlevel: input: a,b,xi,alpha,beta,x,cov = ',a,b,xi,alpha,beta,x,cov
    call getabfromcov(a,b,alpha,beta,cov,aa,bb)
    f = 10.**x
    f = 2/f
    z = serfci(f)
    t = aa + sqrt(2.)*bb*z
    gaucovreturnlevel = t
    !!!!print *,'# gauscovreturnlevel: output: ',t
end function gaucovreturnlevel

real function gaucovreturnyear(a,b,xi,alpha,beta,xyear,cov,lchangesign)
!
!   compute the return time of the value xyear with the fitted values
!
    implicit none
    real a,b,xi,alpha,beta,xyear,cov
    logical lchangesign
    real z,tx,aa,bb
    real erfc

    call getabfromcov(a,b,alpha,beta,cov,aa,bb)        
    z = (xyear - aa)/bb
    if ( z.gt.12 ) then
        tx = 3e33
    else
        tx = 2/erfc(z/sqrt(2.))
    end if
    gaucovreturnyear = tx
end function gaucovreturnyear
