*  #[ fitgau:
        subroutine fitgaucov(xx,yrs,ntot,a,b,alpha,beta,j1,j2
     +       ,lweb,ntype,lchangesign,yr1a,yr2a,xyear,cov1,cov2,offset
     +       ,t,t25,t975,tx,tx25,tx975,assume,lboot,lprint,plot,lwrite)
*
*       a fit a gaussian distribution with mean linearly dependent on a covariate 
*       to the data
*
        implicit none
*
        integer ntot,ntype,j1,j2,yr1a,yr2a
        integer yrs(0:ntot)
        real xx(2,ntot),a,b,alpha,beta,xyear,cov1,cov2,offset,
     +       t(10,3),t25(10,3),t975(10,3),tx(3),tx25(3),tx975(3)
        real ttt(10,3),txtxtx(3),dum
        logical lweb,lchangesign,lboot,lprint,plot,lwrite
        character assume*(*)
*
        integer nmc
        parameter(nmc=1000)
        integer i,j,nx,iter,iens,nfit,imc,ier,year
        real xmin,aa(nmc),bb(nmc),tt(nmc,10,3),dalpha,dbeta
     +       ,z,x,f,txtx(nmc,3),alphaalpha(nmc),betabeta(nmc)
     +       ,mean,sd,ranf,mindata,minindx,pmindata,snorm,s
        real a25,a975,b25,b975,alpha25,alpha975,beta25,beta975
        real adev,var,skew,curt,aaa,siga,chi2,q
        real,allocatable :: yy(:),ys(:),zz(:),sig(:)
        character lgt*4
        external gaucovreturnlevel,gaucovreturnyear
*
        integer nmax
        parameter(nmax=100000)
        integer ncur
        real data(2,nmax),restrain
        logical llwrite
        common /fitdata3/ data
        common /fitdata2/ restrain,ncur,llwrite
        character cassume*5
        common /fitdata4/ cassume
*
        year = yr2a
        if ( lwrite ) then
            print *,'fitgaucov: input'
            print *,'year,xyear  = ',year,xyear
            print *,'cov1,cov2,offset ',cov1,cov2,offset
        end if
*
*       compute first-guess parameters
*
        allocate(yy(ntot))
        allocate(ys(ntot))
        allocate(zz(ntot))
        allocate(sig(ntot))
        do i=1,ntot
            yy(i) = xx(1,i)
            zz(i) = xx(2,i)
        end do
        sig = 0
        call moment(yy,ntot,mean,adev,sd,var,skew,curt)
        call fit(zz,yy,ntot,sig,0,aaa,alpha,siga,dalpha,chi2,q)
        if ( lwrite ) then
            print *,'fitgaucov: computed initialisation values:'
            print *,'mean,sd,alpha,dalpha = ',mean,sd,alpha,dalpha
        end if
*
*       a trivial case which causes no end of trouble
*
        if ( sd.eq.0 ) then
            if ( lwrite ) print *,'fitgaucov: sd=0, everything undfined'
            a = 3e33
            b = 3e33
            alpha = 3e33
            beta = 3e33
            t = 3e33
            t25 = 3e33
            t975 = 3e33
            tx = 3e33
            tx25 = 3e33
            tx975 = 3e33
            return
        endif
*
*       copy to common for routine llgausscov
*
        ncur = ntot
        do i=1,ncur
            data(:,i) = xx(:,i)
        enddo
        cassume = assume
*       
*       fit, using Numerical Recipes routines
*
        a = mean
        b = sd
        if ( assume.eq.'shift' .or. assume.eq.'scale' ) then
            beta = 3e33
            call fit1gaucov(a,b,alpha,dalpha,iter)
        else if ( assume.eq.'both' ) then
            call fit2gaucov(a,b,alpha,beta,dalpha,dbeta,iter)
        else
            write(0,*) 'fitgaucov: cannot handle assume = ',assume
            call abort
        end if
        dum = 0
        call getreturnlevels(a,b,dum,alpha,beta,cov1,cov2,
     +       gaucovreturnlevel,j1,j2,t)
        if ( xyear.lt.1e33 ) then
            call getreturnyears(a,b,dum,alpha,beta,xyear,cov1,cov2,
     +           gaucovreturnyear,j1,j2,tx,lwrite)
        endif
*
*       Bootstrap for error estimate
*
        if ( .not.lboot ) then
            if ( lchangesign ) then
                a = -a
                t = -t
                alpha = -alpha
            endif
            return
        endif
        if ( .not.lweb ) print '(a,i6,a)','# doing a ',nmc
     +        ,'-member bootstrap to obtain error estimates'
        do iens=1,nmc
            if ( .not.lweb .and. mod(iens,100).eq.0 )
     +           print '(a,i6)','# ',iens
            do i=1,ntot
                call random_number(ranf)
                j = 1+int(ntot*ranf)
                if ( j.lt.1 .or. j.gt.ntot ) then
                    write(0,*) 'fitgaucov: error: j = ',j
                    call abort
                endif
                data(:,i) = xx(:,j)
            enddo
            aa(iens) = a
            bb(iens) = b
            alphaalpha(iens) = alpha
            llwrite = .false.
            if ( assume.eq.'shift' .or. assume.eq.'scale' ) then
                betabeta(iens) = 3e33
                call fit1gaucov(aa(iens),bb(iens),
     +               alphaalpha(iens),dalpha,iter)
            else if ( assume.eq.'both' ) then
                betabeta(iens) = beta
                call fit2gaucov(aa(iens),bb(iens),
     +               alphaalpha(iens),betabeta(iens),dalpha,dbeta,iter)
            else
                write(0,*) 'fitgaucov: cannot handle assume = ',assume
                call abort
            end if
            call getreturnlevels(aa(iens),bb(iens),dum,
     +           alphaalpha(iens),betabeta(iens),
     +           cov1,cov2,gaucovreturnlevel,j1,j2,ttt)
            do i=1,10
                do j=1,3
                    tt(iens,i,j) = ttt(i,j)
                end do
            end do
            if ( xyear.lt.1e33 ) then
                call getreturnyears(aa(iens),bb(iens),dum,
     +               alphaalpha(iens),betabeta(iens),xyear,cov1,cov2,
     +               gaucovreturnyear,j1,j2,txtxtx,lwrite)
                do j=1,3
                    txtx(iens,j) = txtxtx(j)
                end do
            endif
        enddo
        if ( lchangesign ) then
            a = -a
            aa = -aa
            alpha = -alpha
            alphaalpha = -alphaalpha
            t = -t
            tt = -tt
        endif
        call getcut( a25, 2.5,nmc,aa)
        call getcut(a975,97.5,nmc,aa)
        call getcut( b25, 2.5,nmc,bb)
        call getcut(b975,97.5,nmc,bb)
        call getcut( alpha25, 2.5,nmc,alphaalpha)
        call getcut(alpha975,97.5,nmc,alphaalpha)
        if ( assume.eq.'both' ) then
            call getcut( beta25, 2.5,nmc,betabeta)
            call getcut(beta975,97.5,nmc,betabeta)
        end if
        do i=1,10
            do j=1,3
                call getcut(t25(i,j),2.5,nmc,tt(1,i,j))
                call getcut(t975(i,j),97.5,nmc,tt(1,i,j))
            enddo
        end do
        do j=1,3
            if ( xyear.lt.1e33 ) then
                call getcut(tx25(j), 2.5,nmc,txtx(1,j))
                call getcut(tx975(j),97.5,nmc,txtx(1,j))
                if ( lchangesign ) xyear = -xyear
            endif
        end do
*
*       output
*
        if ( lweb ) then
            print '(a)','# <tr><td colspan="4">Fitted to normal '//
     +           'distribution P(x) = exp(-(x-a'')&sup2;'//
     +           '/(2b''&sup2;))/(b''&radic;(2&pi;))</td></tr>'
            call printab(lweb)
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'//
     +           'a:</td><td>',a,'</td><td>',a25,'...',a975,'</td></tr>'
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'//
     +           'b:</td><td>',b,'</td><td>',b25,'...',b975,'</td></tr>'
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td colspan=2>'//
     +           '&alpha;:</td><td>',alpha,'</td><td>',alpha25,'...',
     +           alpha975,'</td></tr>'
            if ( assume.eq.'both' ) then
                print '(a,f16.3,a,f16.3,a,f16.3,a)',
     +               '# <tr><td colspan=2>&beta;:</td><td>',beta,
     +               '</td><td>',beta25,'...',beta975,'</td></tr>'
            end if
        else
            print '(a,i5,a)','# Fitted to Gaussian distribution in '
     +           ,iter,' iterations'
            print '(a)','# p(x) = exp(-(x-a'')^2/(2*b''^2))/'//
     +           '(b''*sqrt(2*pi)) with'
            call printab(lweb)
            print '(a,f16.3)','# a = ',a
            print '(a,f16.3)','# b = ',b
            print '(a,f16.3)','# alpha = ',alpha
            if ( assume.eq.'both' ) then
                print '(a,f16.3,a,f16.3,a,f16.3)','# beta  ',beta,
     +               ' \\pm ',beta975-beta25
            end if
        endif
        call printcovreturnvalue(ntype,t,t25,t975,yr1a,yr2a,lweb)
        call printcovreturntime(year,xyear,tx,tx25,tx975,yr1a,yr2a,lweb)
        call plotreturnvalue(ntype,t25(1,2),t975(1,2),j2-j1+1)

        if ( plot ) then
            call plot_tx_cdfs(txtx,nmc,ntype)
        end if

        if ( cassume.eq.'both' ) then
            write(0,*) 'fitgaucov: error: cannot handle plottng yet for'
     +           //'assume = both'
            write(*,*) 'fitgevcov: error: cannot handle plottng yet for'
     +           //'assume = both'
            stop
        end if
        ! compute distribution at year and plot it
        if ( lchangesign ) then
            do i=1,ntot
                yy(i) = yy(i) + alpha*(zz(i)-cov2)
            end do
        else
            do i=1,ntot
                yy(i) = yy(i) - alpha*(zz(i)-cov2)
            end do
        end if
        ys(1:ntot) = yy(1:ntot)
        ! no cuts
        mindata = -2e33
        minindx = -2e33
        pmindata = -1
        snorm = 1
        ! fit to gauss (normal distribution)
        nfit = 2
        aa = a+cov2*alpha
        if ( lchangesign ) b = -b
        call plot_ordered_points(yy,ys,yrs,ntot,ntype,nfit,
     +       aa,b,dum,j1,j2,minindx,mindata,pmindata,
     +       year,xyear,snorm,lchangesign,lwrite)
        ! and print xyear if the distribution would be at yr1a
        print '(a)'
        print '(a)'
        f = 1 - 1/real(ntot+1)*0.9**100
        call printpoint(0,1/real(ntot+1),ntype,-999.9,
     +       xyear+alpha*(cov2-cov1),100*yr1a)
        call printpoint(0,f,ntype,-999.9,
     +       xyear+alpha*(cov2-cov1),100*yr1a)

        end
*  #] fitgaucov:
*  #[ fit1gaucov:
        subroutine fit1gaucov(a,b,alpha,dalpha,iter)
        implicit none
        integer iter
        real a,b,alpha,dalpha
        integer i
        real q(4),p(4,3),y(4),tol
        real llgausscov
        external llgausscov
*
*       fit, using Numerical Recipes routines
*       
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
*       maybe add restart later
        a = p(1,1)
        b = p(1,2)
        alpha = p(1,3)
        end
*  #] fit1gaucov:
*  #[ fit2gaucov:
        subroutine fit2gaucov(a,b,alpha,beta,dalpha,dbeta,iter)
        implicit none
        integer iter
        real a,b,alpha,beta,dalpha,dbeta
        integer i
        real q(4),p(5,4),y(5),tol
        real llgausscov
        external llgausscov
*
*       fit, using Numerical Recipes routines
*       
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
*       maybe add restart later
        a = p(1,1)
        b = p(1,2)
        alpha = p(1,3)
        beta = p(1,4)
        end
*  #] fit2gaucov:
*  #[ llgausscov:
        real function llgausscov(p)
*
*       computes the log-likelihood function for a normal distribution
*       with parameters alpha,beta=p(1),p(2) and data in common.
*
        implicit none
*       
        real p(4)
*
        integer i
        real z,s,aa,bb
*
        integer nmax
        parameter(nmax=100000)
        integer ncur
        real data(2,nmax),restrain
        logical llwrite
        common /fitdata3/ data
        common /fitdata2/ restrain,ncur,llwrite
        character cassume*5
        common /fitdata4/ cassume
*       
        llgausscov = 0
        do i=1,ncur
            call getabfromcov(p(1),p(2),p(3),p(4),data(2,i),aa,bb)
            z = (data(1,i) - aa)/bb
            llgausscov = llgausscov - z**2/2 - log(abs(bb))
        enddo
*       normalization is not 1 in case of cut-offs
        call gauscovnorm(aa,bb,s)
        llgausscov = llgausscov - ncur*log(s)
*       minimum, not maximum
        llgausscov = -llgausscov
***        print *,'a,b,llgausscov = ',p(1),p(2),llgausscov
*
        end
*  #] llgausscov:
*  #[ gauscovnorm:
        subroutine gauscovnorm(a,b,s)
        implicit none
#include "getopts.inc"
        real a,b,s
        real z1,z2,sqrt2
        real erfcc
        external erfcc
        if ( minindx.gt.-1e33 .or. maxindx.lt.1e33 ) then
            write(0,*) 'gauscovnorm: boundaries not yet available for'//
     +           ' fit of gauss(t)'
            call abort
        else
            s = 1
        endif
***        print *,'gauscovnorm: norm = ',a,b,s
        end
*  #] gauscovnorm:
*  #[ gaucovreturnlevel:
        real function gaucovreturnlevel(a,b,xi,alpha,beta,x,cov)
!
!       compute return times given the normal distribution parameters a,b and 
!       x = log10(returntime) for covariant cov and fit parameter alpha
!
        implicit none
        real a,b,xi,alpha,beta,x,cov
        integer ier
        real aa,bb,f,z,t
!
        call getabfromcov(a,b,alpha,beta,cov,aa,bb)
        f = 10.**x
        f = 1-2/f
        call merfi(f,z,ier)
        t = aa + sqrt(2.)*bb*z
        gaucovreturnlevel = t
        end
*  #] gaucovreturnlevel:
*  #[ gaucovreturnyear:
        real function gaucovreturnyear(a,b,xi,alpha,beta,xyear,cov)
!
!       compute the return time of the value xyear with the fitted values
!
        implicit none
        real a,b,xi,alpha,beta,xyear,cov
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

        end
*  #] gaucovreturnyear:
