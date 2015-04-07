*  #[ fitgpdcov:
        subroutine fitgpdcov(xx,yrs,ntot,a3,b3,xi3,alpha3,beta3,j1,j2
     +       ,lweb,ntype,lchangesign,yr1a,yr2a,xyearin,idmax,cov1,cov2
     +       ,offset,t3,tx3,threshold,inrestrain,assume
     +       ,confidenceinterval,ndecor,lboot,lprint,dump,plot,lwrite)
*
*       fit a GPD distribution to the data, which is already assumed to be declustered
*       input:
*       xx(2,ntot) data,covariate
*       j1,j2    use days/months/... j1 to j2
*       year     leave out this year from teh fit and compute return time for it
*       xyear    value for year, has been set to undef in the series
*       inrestrain restrain xi parameter by adding a normal distribution of width 0.5*inrestrain to the cost function
*       threshold in percent
*       assume   shift: only vary threshold, scale: vary threshold & b in unison, both: independently
*       output
*       a,b,xi,alpha,beta     parameters of fit
*       assume   shift: alpha modifies the position parameter a(cov) = a + alpha*cov
*                scale: alpha modifies both the position and shape parameters:
*                       a(cov) = a*exp(alpha*cov/a), b(cov) = b*exp(alpha*cov/a)
*                both:  a(cov) = a + alpha*cov, b(cov) = b + beta*cov
*       t(10,3)    return values for 10, 20, 50, ..., 10000 years for cov=cov1,cov2 and the difference
*       t25,t975   (100-confidenceinterval)/2%, (100+confidenceinterval)/25% quantiles of these return values
*       tx(3)      return time of the value of year (xyear) in the context of the other values and difference
*       t25,t975   (100-confidenceinterval)/2%, (100+confidenceinterval)/2% quantiles of these return values
*
        implicit none
*
        integer nmc
        parameter(nmc=1000)
        integer ntot,j1,j2,ntype,yr1a,yr2a,ndecor
        integer yrs(0:ntot)
        real xx(2,ntot),a3(3),b3(3),xi3(3),alpha3(3),beta3(3),xyearin,
     +       cov1,cov2,offset,inrestrain,t3(3,10,3),tx3(3,3),threshold,
     +       confidenceinterval
        character assume*(*),idmax*(*)
        logical lweb,lchangesign,lboot,lprint,dump,plot,lwrite
*
        integer i,j,jj,k,l,n,nx,iter,iter1,iens,iiens,nfit,year
        real x,a,b,xi,alpha,beta,t(10,3),t25(10,3),t975(10,3),
     +       tx(3),tx25(3),tx975(3),aa(nmc),bb(nmc),xixi(nmc),
     +       alphaalpha(nmc),betabeta(nmc),tt(nmc,10,3),
     +       b25,b975,xi25,xi975,alpha25,alpha975,t5(10,3),t1(10,3),
     +       db,dxi,f,z,ll,ll1,txtx(nmc,3),a25,a975,beta25,beta975,
     +       ranf,mean,sd,dalpha,dbeta,mindata,minindx,pmindata,snorm,s,
     +       xmin,cmin,cmax,c,xxyear,frac,ttt(10,3),txtxtx(3),
     +       acov(3,2),aacov(nmc,2),plo,phi,xyear
        real adev,var,skew,curt,aaa,bbb,siga,chi2,q,p(4)
        integer,allocatable :: ii(:),yyrs(:)
        real,allocatable :: yy(:),ys(:),zz(:),sig(:)
        logical lopen
        character lgt*4,string*1000,arg*250
        integer iargc
*
        integer nmax,ncur
        parameter(nmax=100000)
        real data(2,nmax),restrain
        logical llwrite
        common /fitdata3/ data
        common /fitdata2/ restrain,ncur,llwrite
        character cassume*5
        common /fitdata4/ cassume
        integer nthreshold
        real athreshold,pthreshold
        common /fitdata5/ nthreshold,athreshold,pthreshold
*
        real llgpdcov,gpdcovreturnlevel,gpdcovreturnyear
        external llgpdcov,gpdcovreturnlevel,gpdcovreturnyear
*
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
            print *,'ntot           = ',ntot
            if ( .false. ) then
                do i=1,ntot
                    print *,i,(xx(j,i),j=1,2)
                enddo
            endif
        endif
        a3 = 3e33
        b3 = 3e33
        xi3 = 3e33
        alpha3 = 3e33
        beta3 = 3e33
        t3 = 3e33
        tx3 = 3e33
        if ( ntot.lt.5 ) then
            return
        end if
!
!       thin data to exclude all the values equal to xmin (result of declustering)
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
            if ( xx(1,i).gt.xmin ) then
                ncur = ncur + 1
            end if
        end do
        if ( ncur.lt.10 ) then
            return
        end if
!
!       dump (covariate,observation) pairs to plotfile on unit 15
!
        call write_obscov(xx,yrs,ntot,xmin,cov2,xyear,year,offset,
     +       lchangesign)
*       
*       compute first-guess parameters
*
        allocate(ii(ncur))
        allocate(yyrs(0:ncur))
        allocate(yy(ncur))
        allocate(ys(ncur))
        allocate(zz(2*ncur))
        allocate(sig(ncur))
        j = 0
        yyrs(0) = yrs(0)
        do i=1,ntot
            if ( xx(1,i).gt.xmin ) then
                j = j + 1
                ii(j) = i
                yy(j) = xx(1,i)
                zz(j) = xx(2,i)
                yyrs(j) = yrs(i)
            end if
        end do
        if ( j.ne.ncur ) then
            write(0,*) 'fitgpdcov: error: j != ncur ',j,ncur
            write(*,*) 'fitgpdcov: error: j != ncur ',j,ncur
            call abort
        end if
        sig = 0
        call moment(yy,ncur,mean,adev,sd,var,skew,curt)
        call fit(zz,yy,ncur,sig,0,aaa,alpha,siga,dalpha,chi2,q)
        if ( lwrite ) then
            print *,'fitgpdcov: computed initialisation values:'
            print *,'mean,sd,alpha,dalpha = ',mean,sd,alpha,dalpha
        end if
*
*       ill-defined case
*
        if ( sd.eq.0 ) then
            goto 801 ! deallocate and return
        endif
*
*       copy to common for routine llgpdcov
*
*       number of points above threshold, threshold is relative to the full set (ntot)
        nthreshold = nint(ntot*(1-threshold/100))
        if ( nthreshold.lt.10 ) then
            if ( lprint ) then
                write(0,*) 'fitgpdcov: error: not enough points above'//
     +               ' threshold: ',nthreshold
                write(*,*) 'fitgpdcov: error: not enough points above'//
     +               ' threshold: ',nthreshold
            end if
            goto 801 ! deallocate and return
        end if
        if ( nthreshold.ge.ncur ) then
            write(0,*) 'fitgpdcov: error: nthreshold &gt ncur ',
     +           nthreshold,ncur
            write(0,*) 'fitgpdcov: error: nthreshold &gt ncur ',
     +           nthreshold,ncur
            nthreshold = ncur -1
        end if
        do i=1,ncur
            data(:,i) = xx(:,ii(i))
        enddo
        restrain = inrestrain
        llwrite = lwrite
        cassume = assume

        ! first guess
        ys(1:ncur) = yy(1:ncur)
        call nrsort(ncur,yy)
        athreshold = (yy(ncur-nthreshold) + yy(ncur-nthreshold+1))/2
        if ( lwrite ) then
            print *,'fitgpdcov: nthreshold,ncur,ntot,athreshold = ',
     +           nthreshold,ncur,ntot,athreshold
        end if
        ! needed later on...
        yy(1:ncur) = ys(1:ncur)
        b = sd/3 ! should set the scale roughly right...
        xi = 0
        if ( assume.eq.'shift' .or. assume.eq.'scale' ) then
            beta = 3e33
            call fit1gpdcov(a,b,xi,alpha,dalpha,iter)
        else if ( assume.eq.'both' ) then
            beta = alpha
            dbeta = dalpha
            call fit2gpdcov(a,b,xi,alpha,beta,dalpha,dbeta,iter)
        else
            write(0,*) 'fitgpdcov: error: unknown value for assume ',
     +           assume
        end if
        call getreturnlevels(a,b,xi,alpha,beta,cov1,cov2,
     +       gpdcovreturnlevel,j1,j2,t)
        if ( xyear.lt.1e33 ) then
            call getreturnyears(a,b,xi,alpha,beta,xyear,cov1,cov2,
     +           gpdcovreturnyear,j1,j2,tx,lwrite)
            ! convert to years
            tx(1:2) = tx(1:2)/(j2-j1+1)
        endif
        call getabfromcov(a,b,alpha,beta,cov1,aaa,bbb)
        acov(1,1) = aaa
        call getabfromcov(a,b,alpha,beta,cov2,aaa,bbb)
        acov(1,2) = aaa
        call write_threshold(cmin,cmax,a,b,alpha,beta,offset,
     +       lchangesign)
*
*       bootstrap to find error estimates
*
        if ( .not.lboot ) then
            if ( lchangesign ) then
                a = -a
                b = -b
                t = -t
                alpha = -alpha
                if ( cassume.eq.'both' ) then
                    beta = -beta
                end if
            endif
            a3(1) = a
            b3(1) = abs(b)
            xi3(1) = xi
            alpha3(1) = alpha
            beta3(1) = beta
            t3(1,:,:) = t(:,:)
            tx3(1,:) = tx(:)
            goto 801 ! deallocate and return
        endif
        if ( lprint .and. .not.lweb ) print '(a,i6,a)','# doing a ',nmc
     +        ,'-member bootstrap to obtain error estimates'
        iens = 0
        do iiens=1,nmc
            iens = iens + 1
            if ( lprint ) call keepalive1('Bootstrapping',iiens,nmc)
            if ( lprint .and. .not.lweb .and. mod(iiens,100).eq.0 )
     +           print '(a,i6)','# ',iiens
            n = 1 + (ncur-1)/ndecor
            do i=1,n
                ! we do not have the information here to check whether the
                ! data points were contiguous in the original series...
                ! TODO: propagate that information              
                call random_number(ranf)
                j = 1 + min(ncur-ndecor,int((ncur-ndecor)*ranf))
                if ( j.lt.1 .or. j.gt.ncur ) then
                    write(0,*) 'fitgpdcov: error: j = ',j
                    call abort
                endif
                if ( i.lt.n ) then ! the blocks that fit in whole
                    do jj=0,ndecor-1
                        data(:,1+(i-1)*ndecor+jj) = xx(:,ii(j+jj))
                    end do
                else
                    do jj=0,ndecor-1 ! one more block to the end, the previous block is shortened
                        data(:,1+ncur-ndecor+jj) = xx(:,ii(j+jj))
                    end do
                end if
            enddo
            aa(iens) = a
            bb(iens) = b
            xixi(iens) = xi
            alphaalpha(iens) = alpha
            llwrite = .false.
            if ( assume.eq.'shift' .or. assume.eq.'scale' ) then
                betabeta(iens) = 3e33
                call fit1gpdcov(aa(iens),bb(iens),xixi(iens),
     +               alphaalpha(iens),dalpha,iter1)
            else if ( assume.eq.'both' ) then
                betabeta(iens) = beta
                call fit2gpdcov(aa(iens),bb(iens),xixi(iens),
     +               alphaalpha(iens),betabeta(iens),dalpha,dbeta,iter1)
            else
                write(0,*) 'fitgpdcov: error: unknown value for assume '
     +               ,assume
            end if
            call getabfromcov(aa(iens),bb(iens),
     +           alphaalpha(iens),betabeta(iens),cov1,aaa,bbb)
            aacov(iens,1) = aaa
            call getabfromcov(aa(iens),bb(iens),
     +           alphaalpha(iens),betabeta(iens),cov2,aaa,bbb)
            aacov(iens,2) = aaa
            call getreturnlevels(aa(iens),bb(iens),xixi(iens),
     +           alphaalpha(iens),betabeta(iens),
     +           cov1,cov2,gpdcovreturnlevel,j1,j2,ttt)
            do i=1,10
                do j=1,3
                    tt(iens,i,j) = ttt(i,j)
                end do
            end do
            if ( xyear.lt.1e33 ) then
                call getreturnyears(aa(iens),bb(iens),xixi(iens),
     +               alphaalpha(iens),betabeta(iens),xyear,cov1,cov2,
     +               gpdcovreturnyear,j1,j2,txtxtx,lwrite)
                do j=1,2
                    ! convert to years
                    txtx(iens,j) = txtxtx(j)/(j2-j1+1)
                end do
                txtx(iens,3) = txtxtx(3)
                ! if the event would have been impossible in the current climate
                ! disregard this bootstrap sample
                if ( txtxtx(3).gt.1e19 ) then
                    iens = iens - 1
                end if
            endif
        enddo
        if ( lchangesign ) then
            a = -a
            acov = -acov
            aa = -aa
            aacov = -aacov
            b = -b
            bb = -bb
            alpha = -alpha
            alphaalpha = -alphaalpha
            if ( assume.eq.'both' ) then
                beta = -beta
                betabeta = -betabeta
            end if
            t = -t
            tt = -tt
        endif
        plo = (100-confidenceinterval)/2
        phi = (100+confidenceinterval)/2
        call getcut( a25,plo,iens,aa)
        call getcut(a975,phi,iens,aa)
        call getcut( b25,plo,iens,bb)
        call getcut(b975,phi,iens,bb)
        call getcut( xi25,plo,iens,xixi)
        call getcut(xi975,phi,iens,xixi)
        call getcut( alpha25,plo,iens,alphaalpha)
        call getcut(alpha975,phi,iens,alphaalpha)
        if ( assume.eq.'both' ) then
            call getcut( beta25,plo,iens,betabeta)
            call getcut(beta975,phi,iens,betabeta)
        end if
        call getcut(acov(2,1),plo,iens,aacov(1,1))
        call getcut(acov(3,1),phi,iens,aacov(1,1))
        call getcut(acov(2,2),plo,iens,aacov(1,2))
        call getcut(acov(3,2),phi,iens,aacov(1,2))
        do i=1,10
            do j=1,3
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
        do j=1,3
            if ( xyear.lt.1e33 ) then
                call getcut(tx25(j),plo,iens,txtx(1,j))
                call getcut(tx975(j),phi,iens,txtx(1,j))
                if ( lchangesign ) xyear = -xyear
            endif
        end do
*
*       output
*
        if ( .not.lprint ) then
            call copyab3etc(a3,b3,xi3,alpha3,beta3,t3,tx3,
     +           a,a25,a975,b,b975,xi,xi,xi,alpha,alpha25,alpha975,
     +           beta,beta25,beta975,t,t25,t975,tx,tx25,tx975)
            if ( .not.lwrite ) goto 801 ! deallocate and return
        end if
        if ( lweb ) then
            print '(a)','# <tr><td colspan="4">Fitted to GPD '//
     +           'distribution H(x+a'') = 1 - (1+&xi;*x/b'')^(-1/&xi;)'
     +           //'</td></tr>'
            call printab(lweb)
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'//
     +           'a:</td><td>',a,'</td><td>',a25,'...',a975,'</td></tr>'
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'//
     +           'b:</td><td>',b,'</td><td>',b25,'...',
     +           b975,'</td></tr>'
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'//
     +           '&xi;:</td><td>',xi,'</td><td>',xi25,'...',xi975,
     +           '</td></tr>'
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'//
     +           '&alpha;:</td><td>',alpha,'</td><td>',alpha25,'...',
     +           alpha975,'</td></tr>'
            if ( assume.eq.'both' ) then
                print '(a,f16.3,a,f16.3,a,f16.3,a)',
     +               '# <tr><td colspan=2>&beta;:</td><td>',beta,
     +               '</td><td>',beta25,'...',beta975,'</td></tr>'
            end if
        else
            print '(a,i5,a)','# Fitted to GPD distribution in ',iter
     +           ,' iterations'
            print '(a)','# H(x+a) = 1-(1+xi*x/b)**(-1/xi) with'
            print '(a,f16.3,a,f16.3,a,f16.3)','# a = ',a,' \\pm ',a975
     +           -a25
            print '(a,f16.3,a,f16.3,a,f16.3)','# b = ',b,' \\pm ',b975
     +           -b25
            print '(a,f16.3,a,f16.3,a,f16.3)','# xi  = ',xi,' \\pm ',
     +           xi975-xi25
            print '(a,f16.3,a,f16.3,a,f16.3)','# alpha ',alpha,' \\pm ',
     +           alpha975-alpha25
        end if
        call printcovreturnvalue(ntype,t,t25,t975,yr1a,yr2a,lweb,plot)
        call printcovreturntime(year,xyear,idmax,tx,tx25,tx975,yr1a,yr2a
     +       ,lweb,plot)
!       plot fit for present-day climate
        call plotreturnvalue(ntype,t25(1,2),t975(1,2),j2-j1+1)

        if ( dump ) then
            call plot_tx_cdfs(txtx,nmc,iens,ntype,j1,j2)
        end if
        if ( plot ) write(11,'(3g20.4,a)') alpha,alpha25,alpha975,
     +       ' alpha'
        call write_dthreshold(cov1,cov2,acov,offset,lchangesign)

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
        enddo
        ! compute distribution at past year and plot it
        call adjustyy(ncur,data,assume,a,b,alpha,beta,cov1,
     +       yy,zz,aaa,bbb,lchangesign,lwrite)
        ys(1:ncur) = yy(1:ncur)
        mindata = aaa
        print '(a,i5)','# distribution in year ',yr1a
        call plotreturnvalue(ntype,t25(1,1),t975(1,1),j2-j1+1)
        call plot_ordered_points(yy,ys,yyrs,ncur,ntype,nfit,
     +       frac,aaa,bbb,xi,j1,j2,minindx,mindata,pmindata,
     +       year,xyear,snorm,lchangesign,lwrite,.false.)

        ! compute distribution at present year and plot it
        call adjustyy(ncur,data,assume,a,b,alpha,beta,cov2,
     +       yy,zz,aaa,bbb,lchangesign,lwrite)
        ys(1:ncur) = yy(1:ncur)
        mindata = aaa
        print '(a)'
        print '(a)'
        print '(a,i5)','# distribution in year ',yr2a
        call plotreturnvalue(ntype,t25(1,2),t975(1,2),j2-j1+1)
        call plot_ordered_points(yy,ys,yyrs,ncur,ntype,nfit,
     +       frac,aaa,bbb,xi,j1,j2,minindx,mindata,pmindata,
     +       year,xyear,snorm,lchangesign,lwrite,.true.)
 801    continue
        deallocate(ii)
        deallocate(yyrs)
        deallocate(yy)
        deallocate(ys)
        deallocate(zz)
        deallocate(sig)
        return
        end
*  #] fitgpdcov:
*  #[ fit1gpdcov:
        subroutine fit1gpdcov(a,b,xi,alpha,dalpha,iter)
        implicit none
        integer iter
        real a,b,xi,alpha,dalpha
        integer i
        real q(4),p(4,3),y(4),tol
        real llgpdcov
        external llgpdcov
        integer nthreshold
        real athreshold,pthreshold
        common /fitdata5/ nthreshold,athreshold,pthreshold
*
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
        do i=1,4
            q(1) = p(i,1)
            q(2) = p(i,2)
            q(3) = p(i,3)
            y(i) = llgpdcov(q)
        enddo
        tol = 1e-4
        call amoeba(p,y,4,3,3,tol,llgpdcov,iter)
*       maybe add restart later
        a = athreshold
        b = p(1,1)
        xi = p(1,2)
        alpha = p(1,3)
        end
*  #] fit1gpdcov:
*  #[ fit2gpdcov:
        subroutine fit2gpdcov(a,b,xi,alpha,beta,dalpha,dbeta,iter)
        implicit none
        integer iter
        real a,b,xi,alpha,beta,dalpha,dbeta
        integer i
        real q(4),p(5,4),y(5),tol
        real llgpdcov
        external llgpdcov
        integer nthreshold
        real athreshold,pthreshold
        common /fitdata5/ nthreshold,athreshold,pthreshold
*
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
        do i=1,5
            q(1) = p(i,1)
            q(2) = p(i,2)
            q(3) = p(i,3)
            q(4) = p(i,4)
            y(i) = llgpdcov(q)
        enddo
        tol = 1e-4
        call amoeba(p,y,5,4,4,tol,llgpdcov,iter)
*       maybe add restart later
        a = athreshold
        b = p(1,1)
        xi = p(1,2)
        alpha = p(1,3)
        beta = p(1,4)
        end
*  #] fit2gpdcov:
*  #[ llgpdcov:
        real function llgpdcov(p)
*
*       computes the log-likelihood function for a covariant-dependent GPD distribution
*       with parameters a,b,xi,alpha=p(1-4) and data in common.
*
        implicit none
        integer maxloop
        parameter(maxloop=20)
*
        real p(4)
*
        integer i,j,n,nlo,nhi,iloop
        real x,z,xi,s,aa,bb,llold,alpha,alo,ahi,delta
        real,allocatable :: xx(:),xxunsorted(:)
*
        integer nmax,ncur
        parameter(nmax=100000)
        real data(2,nmax),restrain
        logical llwrite
        common /fitdata3/ data
        common /fitdata2/ restrain,ncur,llwrite
        character cassume*5
        common /fitdata4/ cassume
        integer nthreshold
        real athreshold,pthreshold
        common /fitdata5/ nthreshold,athreshold,pthreshold
*
        llgpdcov = 0
        allocate(xx(ncur))
        allocate(xxunsorted(ncur))
        xi = p(2)
        if ( abs(xi).gt.10 ) then
            if ( llwrite ) print *,'llgpdcov: |xi|>10: ',p(2)
            llgpdcov = 3e33
            goto 999
        endif
        if ( restrain.lt.0 ) then
            write(0,*) 'llgpdcov: restrain<0 ',restrain
            call abort
        end if
!
!       get threshold
!
        alpha = p(3)
        if ( cassume.ne.'scale' ) then
            do i=1,ncur
                xx(i) = data(1,i) - alpha*data(2,i)
                !!!print *,'xx(',i,') = ',xx(i),data(1,i),alpha*data(2,i)
            end do
            if ( llwrite ) xxunsorted = xx
            call nrsort(ncur,xx)
            athreshold = (xx(ncur-nthreshold) + xx(ncur-nthreshold+1))/2
            if ( llwrite ) then
                print *,'llgpdcov: nthreshold,athreshold = ',nthreshold,
     +               athreshold
                print *,'          xx(',ncur-nthreshold,') = ',
     +               xx(ncur-nthreshold)
                print *,'          xx(',ncur-nthreshold+1,') = ',
     +               xx(ncur-nthreshold+1)
            end if
            xx = xx - athreshold
            if ( llwrite ) xxunsorted = xxunsorted - athreshold
        else
            ! iterative procedure I am afraid...
            ! assume athreshold has been set to a reasonable value higher up
            if ( llwrite ) then
                print *,'llgpdcov: iterative threshold determination'
                print *,'athreshold,alpha = ',athreshold,alpha
            end if
            nlo = -1
            nhi = -1
            alo = -3e33
            ahi = +3e33
            delta = abs(p(1))
            do iloop=1,maxloop
                n = 0
                do i=1,ncur
                    xx(i) = data(1,i)
     +                   - athreshold*exp(alpha*data(2,i)/athreshold)
                    if ( xx(i).gt.0 ) then
                        n = n + 1
                    end if
                end do
                if ( n.eq.nthreshold ) then
                    exit
                end if
                if ( nhi.eq.-1 .or. nlo.eq.-1 ) then
                    ! bracket solution
                    if ( n.lt.nthreshold ) then
                        ahi = athreshold
                        nhi = n
                        athreshold = athreshold - delta
                    else
                        alo = athreshold
                        nlo = n
                        athreshold = athreshold + delta
                    end if
                    delta = 2*delta
                else
                    ! bisect
                    if ( n.lt.nthreshold ) then
                        ahi = athreshold
                        nhi = n
                    else
                        alo = athreshold
                        nlo = n
                    end if
                    athreshold = (alo+ahi)/2
                end if
                if ( llwrite ) print *,nthreshold,alo,nlo,ahi,nhi
            end do
            if ( iloop.gt.maxloop ) then
                athreshold = alo ! better a bit too low than too high, which gives z<0
                if ( llwrite ) then
                    write(*,'(a,i3,2i5,f10.2)') '# llgpdcov: warning:  '
     +                   //'threshold computation did not converge ',
     +                   iloop,n,nthreshold,athreshold
                end if
            end if
            if ( llwrite ) print *,'athreshold = ',athreshold
            call nrsort(ncur,xx)
        end if
!
!       and compute cost function
!
        do i=ncur-nthreshold+1,ncur
            call getabfromcov(athreshold,p(1),p(3),p(4),data(2,i),aa,bb)
            if ( abs(bb).lt.1e-30 ) then
                if ( llwrite ) print *,'llgpdcov: bb too small ',bb
                llgpdcov = 3e33
                goto 999
            end if
            z = xx(i)
            if ( z.lt.0 ) then
                if ( z.lt.-5e-5*abs(athreshold) ) then
                    write(0,*) 'llgpdcov: error: z<0 ',z,i,ncur
                end if
                z = 0
            endif
            s = 1+xi*z/bb
            if ( s.le.0 ) then
                if ( llwrite ) print *,'llgpdcov: 1+xi*z/bb < 0 ',s
                llgpdcov = 3e33
                goto 999
            endif
            llold = llgpdcov
            if ( abs(xi).lt.1e-4 ) then
                llgpdcov = llgpdcov - z/bb + (z/bb)**2*xi/2
***                print *,i,z, - z/b + (z/b)**2*xi/2 - log(b)
            else
                llgpdcov = llgpdcov - (1+1/xi)*log(1+xi*z/bb)
***                print *,i,z, - (1+1/xi)*log(1+xi*z/b) - log(b)
            endif
            if ( llwrite ) then
                do j=1,ncur
                    if ( xxunsorted(j).eq.xx(i) ) exit
                end do
                print *,i,data(1,j),data(2,j),z,bb,llgpdcov-llold
            end if
            llgpdcov = llgpdcov - log(bb)
        enddo
*       normalization is not 1 in case of cut-offs
        call gpdcovnorm(athreshold,abs(p(1)),p(2),p(3),p(4),s)
        if ( s.lt.1e33 ) then
            llgpdcov = llgpdcov - ncur*log(s)
        else
            if ( llwrite ) print *,'cannot handle cuts yet'
            llgpdcov = 3e33
            goto 999
        end if
        if ( restrain.ne.0 ) then
*           preconditioning on xi with gaussian of width restrain/2
*           around 0
            llgpdcov = llgpdcov - (xi/(restrain/2))**2/2
        endif
*       minimum, not maximum
        llgpdcov = -llgpdcov
*
  999   continue
        deallocate(xx)
        if ( llwrite ) print *,'a,b,xi,alpha,llgpdcov = ',
     +       athreshold,p(1),p(2),p(3),llgpdcov
        end
*  #] llgpdcov:
*  #[ gpdcovnorm:
        subroutine gpdcovnorm(a,b,xi,alpha,beta,s)
        implicit none
	include "getopts.inc"
        real a,b,xi,alpha,beta,s
        real z1,z2

        if ( minindx.gt.-1e33 .or. maxindx.lt.1e33 ) then
            write(0,*) 'gpdcovnorm: boundaries not yet avaiable for '//
     +           'fit of GPD(t)'
            call abort
        else
            s = 1
        endif
***        print *,'gpdcovnorm: norm = ',a,b,s
        end
*  #] gpdcovnorm:
*  #[ gpdcovreturnlevel:
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
        integer nthreshold
        real athreshold,pthreshold
        common /fitdata5/ nthreshold,athreshold,pthreshold

        lwrite = .false. ! local variable
        call getabfromcov(a,b,alpha,beta,cov,aa,bb)
        x = xx + log10(1-pthreshold/100)
        if ( abs(xi).gt.10 ) then
            gpdcovreturnlevel = 3e33
        else if ( abs(xi).lt.1e-4 ) then
            t = bb*x*log(10.) + 0.5*xi*(x*log(10.))**2
        else
            y = xi*x*log(10.)
            if ( y.lt.46 ) then
                t = bb/xi*(-1 + exp(y))
            else
                t = 1e20
            end if
        end if
        t = t + aa ! threshold
        if ( lwrite ) then
            print *,'gpdcovreturnlevel:'
            print *,'a,b,xi,alpha,beta = ',a,b,xi,alpha,beta
            print *,'x,cov             = ',x,cov
            print *,'aa,bb             = ',aa,bb
            print *,'gpdcovreturnlevel = ',t
        end if
        gpdcovreturnlevel = t
        end
*  #] gpdcovreturnlevel:
*  #[ gpdcovreturnyear:
        real function gpdcovreturnyear(a,b,xi,alpha,beta,xyear,cov)
!
!       compute the return time of the value xyear with the fitted values
!
        implicit none
        real a,b,xi,alpha,beta,xyear,cov
        integer i,n,ntot
        real x,y,z,tx,aa,bb
        integer nmax,ncur
        parameter(nmax=100000)
        real data(2,nmax),restrain
        logical llwrite
        common /fitdata3/ data
        common /fitdata2/ restrain,ncur,llwrite
        character cassume*5
        common /fitdata4/ cassume
        integer nthreshold
        real athreshold,pthreshold
        common /fitdata5/ nthreshold,athreshold,pthreshold

        if ( llwrite ) then
            print *,'gpdcovreturnyear: input'
            print *,'a,b,xi     = ',a,b,xi
            print *,'alpha,beta = ',alpha,beta
            print *,'xyear,cov  = ',xyear,cov
        end if
        x = xyear
        call getabfromcov(a,b,alpha,beta,cov,aa,bb)
        if ( xyear.gt.aa ) then
            x = xyear - aa
            z = (1 + xi*x/bb)
            if ( z.gt.0 .and. abs(xi).gt.1e-3 ) then
                tx = z**(1/xi)/(1-pthreshold/100)
            else if ( z.gt.0 ) then
                tx = exp(z - 0.5*xi*z**2)/(1-pthreshold/100)
            else
                tx = 1e20
            end if
            if ( .false. .and. tx.gt.1e20 ) then
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
            if ( cassume.eq.'shift' ) then
                do i=1,ncur
                    if ( data(1,i) - alpha*(data(2,i)-cov)
     +                   .gt.xyear ) n = n + 1
                enddo
            else if ( cassume.eq.'scale' ) then
                do i=1,ncur
                    if ( data(1,i)*exp(alpha/athreshold*(data(2,i)-cov))
     +                   .gt.xyear ) n = n + 1
                enddo
            else if ( cassume.eq.'both' ) then
                do i=1,ncur
                    ! not sure whether this is correct
                    if ( (data(1,i)-xyear - alpha*(data(2,i)-cov))
     +                   .gt.0 ) n = n + 1 
                enddo
            else
                write(0,*) 'gpdcovreturnyear: error: unknown value '//
     +               'for assume: ',cassume
                call abort
            end if
            ! approximately... I do not have this information here.
            ntot = nint(nthreshold/(1-pthreshold/100))
            tx = real(ntot+1)/real(n)
            if ( llwrite ) then
                print *,'gpdcovreturnyear: interpolated tx = ',tx
                print *,'  n,ntot = ',n,ntot
            end if
        endif
        if ( .true. .and. tx.lt.1e-10 ) then
            write(0,*) 'gpdcovreturnyear: tx < 1e-10: ',tx
            write(0,*) 'a,b,xi,alpha,xyear = ',a,b,alpha,xi,xyear
        endif
        gpdcovreturnyear = tx
        end
*  #] gpdcovreturnyear:
