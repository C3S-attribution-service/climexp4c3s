*  #[ fitgpd:
        subroutine fitgpd(xx,ntot,mean,sd,b,xi,j1,j2,lweb,ntype
     +       ,lchangesign,pthreshold,threshold,year,xyear,t,t25,t975,tx
     +       ,tx25,tx975,inrestrain,lboot,lprint,lwrite)
*
*       a fit a GPD distribution to the data
*       input:
*       xx(ntot) data
*       mean,sd  for first guess parameters
*       j1,j2    use days/months/... j1 to j2
*       output
*       b,xi     parameters of fit
*       t(10)    return values for 10, 20, 50, ..., 10000 years
*
        implicit none
*
        integer nmc
        parameter(nmc=1000)
        integer ntot,j1,j2,ntype,year
        real xx(ntot),mean,sd,b,pthreshold,inrestrain,xyear,
     +       t(10),t25(10),t975(10),tx,tx25,tx975
        logical lweb,lchangesign,lboot,lprint,lwrite
*
        integer i,j,k,n,nx,iter,iens,iweird,nweird
        real x,xi,bb(nmc),xixi(nmc),tt(nmc,10),b25
     +       ,b975,xi25,xi975,t5(10),t1(10),db,dxi,f
     +       ,threshold,thens,z,ll,ll1,txtx(nmc),ranf
        character lgt*4
        save iweird,nweird
        data iweird,nweird /0,1/
*
        integer nmax,ncur
        parameter(nmax=100000)
        real data(nmax),restrain
        logical llwrite
        common /fitdata1/ data
        common /fitdata2/ restrain,ncur,llwrite
*
        real llgpd
        external llgpd
*
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
                enddo
            endif
        endif
*
*       ill-defined case
*
        if ( sd.eq.0 ) then
            xi = 3e33
            b = 3e33
            t = 3e33
            t25 = 3e33
            t975 = 3e33
            tx = 3e33
            tx25 = 3e33
            tx975 = 3e33
            return
        endif
*
*       copy to common for routine llgpd
*
        if ( pthreshold.gt.1e33 .or. pthreshold.eq.0 ) then
            write(*,'(a)') '# histogram: please specify threshold when'
     +           //' fitting GPD to tail'
            call abort
        endif
        call getcut(threshold,pthreshold,ntot,xx)
        ncur = 0
        do i=1,ntot
            if ( xx(i).ge.threshold ) then
                ncur = ncur + 1
                data(ncur) = xx(i) - threshold
            endif
        enddo
        restrain = inrestrain
        if ( lwrite ) print *,'fitgpd: found ',ncur
     +       ,' points above threshold'
!       make sure that points with a lot of equal values at the
!       threshold are not included...
        if ( ncur.lt.3 .or. ncur.gt.2*(1-pthreshold/100)*ntot ) then
            xi = 3e33
            b = 3e33
            t = 3e33
            t25 = 3e33
            t975 = 3e33
            tx = 3e33
            tx25 = 3e33
            tx975 = 3e33
            return
        endif            
        b = sd
        xi = 0
        call fit1gpd(b,xi,iter)
        if ( b.lt.1e-6*sd ) then
*           something went wrong, throw away results
            iweird = iweird + 1
            if ( iweird.ge.nweird ) then
                nweird = 2*nweird
                print '(a,2g16.4,2i8)','# fitgpd: weird result for b = '
     +               ,b,xi,ncur,iweird
            endif
            xi = 3e33
            b = 3e33
            t = 3e33
            t25 = 3e33
            t975 = 3e33
            tx = 3e33
            tx25 = 3e33
            tx975 = 3e33
            return
        endif
        do i=1,10
            if ( mod(i,3).eq.1 ) then
                x = 1 + i/3
            elseif ( mod(i,3).eq.2 ) then
                x = 1 + log10(2.) + i/3
            else
                x = log10(5.) + i/3
            endif
            x = x + log10(real(j2-j1+1)) + log10(1-pthreshold/100)
            if ( abs(xi).lt.1e-4 ) then
                t(i) = b*x*log(10.) + 0.5*xi*(x*log(10.))**2
            else
                t(i) = b/xi*(-1 + exp(xi*x*log(10.)))
            endif
            t(i) = t(i) + threshold
        enddo

        if ( xyear.lt.1e33 ) then
            if ( xyear.gt.threshold ) then
                x =  xyear - threshold
                z = (1 + xi*x/b)
                if ( z.gt.0 .and. abs(xi).gt.1e-3 ) then
                    tx = z**(1/xi)/(1-pthreshold/100)
                else if ( z.gt.0 ) then
                    tx = exp(z - 0.5*xi*z**2)/(1-pthreshold/100)
                else
                    tx = 1e20
                end if
                if ( tx.gt.1e20 ) then
                    write(0,*) 'fitgpd: tx > 1e20: ',tx
                    write(0,*) 'z,xi = ',z,xi
                end if
            else
                n = 0
                do i=1,ntot
                    if ( xx(i).gt.xyear ) n = n + 1
                enddo
                tx = real(ntot+1)/real(n)
            endif
            if ( lwrite ) then
                print *,'return time = ',tx
            end if
        endif
*
*       bootstrap to find error estimates
*
        if ( .not.lboot ) then
            if ( lchangesign ) then
                b = -b
                t = -t
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
                    write(0,*) 'fitgpd: error: j = ',j
                    call abort
                endif
                data(i) = xx(j)
            enddo
            call getcut(thens,pthreshold,ntot,data)
*           as a side effect, data is now sorted
            do i=ntot,1,-1
                data(i) = data(i) - thens
                if ( data(i).lt.0 ) exit
            enddo
            ncur = ntot - i
            do i=1,ncur
                data(i) = data(i+ntot-ncur)
            enddo
            bb(iens) = b
            xixi(iens) = 0
            call fit1gpd(bb(iens),xixi(iens),iter)
            do i=1,10
                if ( mod(i,3).eq.1 ) then
                    x = 1 + i/3
                elseif ( mod(i,3).eq.2 ) then
                    x = 1 + log10(2.) + i/3
                else
                    x = log10(5.) + i/3
                endif
                x = x + log10(real(j2-j1+1)) + log10(1-pthreshold/100)
                if ( abs(xixi(iens)).lt.1e-4 ) then
                    tt(iens,i) = bb(iens)*x*log(10.)
                else
                    tt(iens,i) = bb(iens)/xixi(iens)*(-1 + 
     +                   exp(xixi(iens)*x*log(10.)))
                endif
                tt(iens,i) = tt(iens,i) + thens
            enddo
            if ( xyear.lt.1e33 ) then
                if ( xyear.gt.thens ) then
                    x =  xyear - thens
                    if ( abs(xixi(iens)).lt.1e-4 ) then
                        z = x/bb(iens) - 0.5*(x/bb(iens))**2
                        txtx(iens) = exp(z)/(1-pthreshold/100)
                    else
                        if ( bb(iens).gt.1e-20 ) then
                            z = (1 + xixi(iens)*x/bb(iens))
                        else
                            z = 1e20
                        endif
                        if ( z.gt.0 ) then
                            z = log(z)/xixi(iens)
                            if ( z.lt.65 ) then
                                txtx(iens) = exp(z)/(1-pthreshold/100)
                            else
                                txtx(iens) = 1e20
                            endif
                        else
                            txtx(iens) = 1e20
                        endif
                    endif
                else
                    n = 0
                    do i=1,ntot
                        if ( data(i).gt.xyear-thens ) n = n + 1
                    enddo
                    txtx(iens) = real(ntot+1)/real(n)
                endif
                if ( lwrite ) print *,'return time ',iens,txtx(iens)
            endif
        enddo
        if ( lchangesign ) then
            b = -b
            bb = -bb
            t = -t
            tt = -tt
        endif
        call getcut( b25, 2.5,nmc,bb)
        call getcut(b975,97.5,nmc,bb)
        call getcut( xi25, 2.5,nmc,xixi)
        call getcut(xi975,97.5,nmc,xixi)
        do i=1,10
            if ( lchangesign ) then
                lgt = '&lt;'
                call getcut(t5(i),5.,nmc,tt(1,i))
                call getcut(t1(i),1.,nmc,tt(1,i))
            else
                lgt = '&gt;'
                call getcut(t5(i),95.,nmc,tt(1,i))
                call getcut(t1(i),99.,nmc,tt(1,i))
            endif
            call getcut(t25(i),2.5,nmc,tt(1,i))
            call getcut(t975(i),97.5,nmc,tt(1,i))
        enddo
        if ( xyear.lt.1e33 ) then
            call getcut(tx25, 2.5,nmc,txtx)
            call getcut(tx975,97.5,nmc,txtx)
        endif
*
*       output
*
        if ( .not.lprint .and. .not.lwrite ) return
        if ( lweb ) then
            print '(a)','# <tr><td colspan="3">Fitted to GPD '//
     +           'distribution H(x) = 1 - (1+&xi;*x/b)^(-1/&xi;)'
     +           //'</td></tr>'
            print '(a,f16.3,a,f16.3,a,f16.3,a)','# <tr><td>b:</td><td>'
     +           ,b,'</td><td>',b25,'...',b975,'</td></tr>'
            print '(a,f16.3,a,f16.3,a,f16.3,a)'
     +           ,'# <tr><td>&xi;:</td><td>',xi,'</td><td>',xi25,'...'
     +           ,xi975,'</td></tr>'
        else
            print '(a,i5,a)','# Fitted to GPD distribution in ',iter
     +           ,' iterations'
            print '(a)','# H(x) = 1-(1+xi*x/b)**(-1/xi) with'
            print '(a,f16.3,a,f16.3,a,f16.3)','# b = ',b,' \\pm ',b975
     +           -b25
            print '(a,f16.3,a,f16.3,a,f16.3)','# xi= ',xi,' \\pm ',xi975
     +           -xi25
        endif
        call printreturnvalue(ntype,t,t25,t975,lweb)
        call printreturntime(year,xyear,tx,tx25,tx975,lweb)
        call plotreturnvalue(ntype,t25,t975,j2-j1+1)
        end
*  #] fitgpd:
*  #[ fit1gpd:
        subroutine fit1gpd(b,xi,iter)
        implicit none
        integer iter
        real b,xi
        integer i
        real q(2),p(3,2),y(3),tol
        real llgpd
        external llgpd
*
*       fit, using Numerical Recipes routines
*
        !!!print *,'fit1gpd: b,xi = ',b,xi
        q(1) = b
        q(2) = xi
        p(1,1) = q(1) *0.9
        p(1,2) = q(2) *0.9
        p(2,1) = p(1,1) *1.2
        p(2,2) = p(1,2)
        p(3,1) = p(1,1)
        p(3,2) = p(1,2) *1.2 + 0.1
        do i=1,3
            q(1) = p(i,1)
            q(2) = p(i,2)
            y(i) = llgpd(q)
        enddo
        tol = 1e-4
        call amoeba(p,y,3,2,2,tol,llgpd,iter)
*       maybe add restart later
        b = p(1,1)
        xi = p(1,2)
        if ( abs(xi).gt.10 ) then
            write(0,*) 'fit1gpd: error: shape parameter xi = ',xi
        end if
        if ( abs(b).gt.1e15 ) then
            write(0,*) 'fit1gpd: error: position parameter b = ',b
        end if
        end
*  #] fit1gpd:
*  #[ llgpd:
        real function llgpd(p)
*
*       computes the log-likelihood function for a GPD distribution
*       with parameters a,b=p(1),p(2) and data in common.
*
        implicit none
*       
        real p(2)
*
        integer i
        real b,xi,s,z,llold
*
        integer nmax,ncur
        parameter(nmax=100000)
        real data(nmax),restrain
        logical llwrite
        common /fitdata1/ data
        common /fitdata2/ restrain,ncur,llwrite
*
        llgpd = 0
        b = p(1)
        xi = p(2)
        !!!print *,'llgpd:b,xi = ',b,xi
        if ( b.le.1e-15 ) then
            llgpd = 3e33
            goto 999
        endif
        if ( abs(xi).gt.10 ) then
            llgpd = 3e33
            goto 999
       endif
        if ( restrain.lt.0 ) then
            write(0,*) 'llgpd: restrain<0 ',restrain
            call abort
        end if
        do i=1,ncur
            z = data(i)
            if ( z.lt.0 ) then
                write(0,*) 'llgpd: z<0 ',z,i,ncur
                call abort
            endif
            if ( 1+xi*z/b.le.0 ) then
                llgpd = 3e33
                goto 999
            endif
            llold = llgpd
            if ( abs(xi).lt.1e-4 ) then
                llgpd = llgpd - z/b + (z/b)**2*xi/2
***                print *,i,z, - z/b + (z/b)**2*xi/2 - log(b)
            else
                llgpd = llgpd - (1+1/xi)*log(1+xi*z/b)
***                print *,i,z, - (1+1/xi)*log(1+xi*z/b) - log(b)
            endif
            if ( llwrite ) print *,i,z,b,llgpd-llold
        enddo
        if ( restrain.eq.0 ) then
            llgpd = llgpd - ncur*log(b)
        else
*           the -2 gives visually much better results;
*           I am unsure of the mathematical derivation
            llgpd = llgpd - (ncur-2)*log(b)
*           preconditioning on xi with gaussian of width restrain/2
*           around 0
            llgpd = llgpd - (xi/(restrain/2))**2/2
        endif
*       minimum, not maximum
        llgpd = -llgpd
 999    continue
!!!        print '(a,i3,2f6.2,f10.4)','ncur,b,xi,llgpd = ',ncur,p(1),p(2)
!!!     +       ,llgpd
*
        end
*  #] llgpd:
