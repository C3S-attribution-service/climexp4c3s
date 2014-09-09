*  #[ fitgam:
        subroutine fitgam(xx,ntot,mean,sd,a,b)
*
*       a fit a gamma distribution to the data
*
        implicit none
*
        integer ntot
        real xx(ntot),mean,sd,a,b
*
        integer i,nx,iter,init
        real tol,p(3,2),q(2),xmin,fret,y(3)
        save init
*
        integer nmax,ncur
        parameter(nmax=100000)
        real data(nmax),restrain
*
        real llgamma
        external llgamma,dllgamma
*       
        data init /0/
*
*       check input
*       
        xmin = 3e33
        do i=1,ntot
            if ( xx(i).lt.0 ) then
                write(0,*) 'histogram: error: cannot fit Gamma to '//
     +                'negative data: ',xx(i)
                call abort
            endif
            if ( xx(i).ne.0 ) then
                xmin = min(xmin,xx(i))
            endif
        enddo
        xmin = xmin/5
        do i=1,ntot
            if ( xx(i).eq.0 ) then
                if ( init.eq.0 ) then
                    init = 1
                    print '(a,f16.6)','# changed zero to ',xmin
                endif
                xx(i) = xmin
            endif
        enddo
*
*       copy to common for routines llgamma, dllgamma
*       
        ncur = ntot
        do i=1,ncur
            data(i) = xx(i)
        enddo
*       
*       fit, using Numerical Recipes routines
*
***     This frequently crashes
***        q(1) = mean**2/sd**2
***        q(2) = sd**2/mean
***        tol = 1e-4
***        call dfpmin(q,2,tol,iter,fret,llgamma,dllgamma)
***        a = q(1)
***        b = q(2)
***     so try amoeba - slow but sure
        p(1,1) = mean**2/sd**2-0.05
        p(1,2) = sd**2/mean    *0.9
        p(2,1) = p(1,1) +0.1
        p(2,2) = p(1,2)
        p(3,1) = p(1,1)
        p(3,2) = p(1,2) *1.2
        do i=1,3
            q(1) = p(i,1)
            q(2) = p(i,2)
            y(i) = llgamma(q)
        enddo
        tol = 1e-4
        call amoeba(p,y,3,2,2,tol,llgamma,iter)
*       maybe add restart later
        a = p(1,1)
        b = p(1,2)
*
*       output
*
        print '(a,i5,a)','# Fitted to gamma distribution in ',iter
     +        ,' iterations'
        print '(a)','# p(x) = (x/b)^(a-1)*exp(-x/b)/(b*Gamma(a)) with'
        print '(a,f16.3)','# a = ',a
        print '(a,f16.3)','# b = ',b
        end
*  #] fitgam:
*  #[ llgamma:
        real function llgamma(p)
*
*       computes the log-likelihood function for a Gamma distribution
*       with parameters alpha,beta=p(1),p(2) and data in common.
*       
        implicit none
*       
        real p(2)
*
        integer i
*
        integer nmax,ncur
        parameter(nmax=10000000)
        real data(nmax),xrestrain
        logical llwrite
        common /fitdata1/ data
        common /fitdata2/ xrestrain,ncur,llwrite
#include "getopts.inc"
*       
        real gammln,gammp,gammq
        external gammln,gammp,gammq
*
        llgamma = 0
        do i=1,ncur
            llgamma = llgamma + 
     +            (p(1)-1)*log(data(i)/p(2)) - data(i)/p(2) - log(p(2))
     +            - gammln(p(1))
        enddo
*       normalization is not 1 in case of cut-offs
        if ( minindx.gt.-1e33 ) then
            if ( maxindx.lt.1e33 ) then
                llgamma = llgamma - ncur*log(gammp(p(1),maxindx/p(2))
     +                - gammp(p(1),minindx/p(2)))
            else
                llgamma = llgamma - ncur*log(gammq(p(1),minindx/p(2)))
            endif
        else
            if ( maxindx.lt.1e33 ) then
                llgamma = llgamma - ncur*log(gammp(p(1),maxindx/p(2)))
            endif
        endif
*       minimum, not maximum
        llgamma = -llgamma
***        print *,'a,b,llgamma = ',p(1),p(2),llgamma
*
        end
*  #] llgamma:
*  #[ dllgamma:
        subroutine dllgamma(p,dp)
*
*       computes the derivatives of the log-likelihood function for a
*       Gamma distribution with parameters p(1),p(2) and data in common.
*       currently unused.
*       
        implicit none
*       
        real p(2),dp(2)
*
        integer i
        real p1(2),p2(2),d
*
        integer nmax,ncur
        parameter(nmax=10000000)
        real data(nmax),restrain
        logical llwrite
        common /fitdata1/ data
        common /fitdata2/ restrain,ncur,llwrite
*
        real dgammln,llgamma
        external dgammln,llgamma
*
        dp(1) = 0
        dp(2) = 0
        do i=1,ncur
            dp(1) = dp(1)
     +            + log(data(i)/p(2)) - dgammln(p(1))
            dp(2) = dp(2)
     +            - p(1)/p(2) + data(i)/p(2)**2
        enddo
*       minimum, not maximum
        dp(1) = -dp(1)
        dp(2) = -dp(2)
*
        print *,'dp(1) = ',dp(1)
        p1(1) = p(1) + 1e-3
        p1(2) = p(2)
        p2(1) = p(1) - 1e-3
        p2(2) = p(2)
        d = (llgamma(p1)-llgamma(p2))/2e-3
        print *,'    cmp ',d
        print *,'dp(2) = ',dp(2)
        p1(1) = p(1)
        p1(2) = p(2) + 10
        p2(1) = p(1)
        p2(2) = p(2) - 10
        d = (llgamma(p1)-llgamma(p2))/20
        print *,'    cmp ',d
*
        end
*  #] dllgamma:
*  #[ gamdist:
        real function gamdist(x)
*
*       Gamma distribution, parameters passed in common
*       currently unused
*
        implicit none
        real x
        real alpha,beta
        common /gamcom/ alpha,beta
        real z,y
        real gammln
        external gammln
        if ( x.lt.0 ) then
            write(0,*) 'gamdist: error: cannot evaluate for x&lt;0: ',x
            call abort
        elseif ( x.eq.0 ) then
            if ( alpha.eq.1 ) then
                gamdist = 1/beta
                return
            elseif ( alpha.gt.1 ) then
                gamdist = 0
                return
            else
                write(0,*) 'gamdist: error: infinite for x=0, alpha= '
     +                ,alpha
                call abort
            endif
        endif
        z = x/beta
        y = (alpha-1)*log(z) - z - log(beta) - gammln(alpha)
        gamdist = exp(y)
        end
*  #] gamdist:
*  #[ cumgamm:
        real function cumgamm(x)
*
*       compute the cumulative Gamma probability minus some
*       requested frequency.  Parameters are passed in common
*
        implicit none
        real x
        real pc,ac,bc
        common /ccumgamm/ pc,ac,bc
        real gammp
        external gammp
*
        if ( x.le.0 ) then
            cumgamm = -pc
        else
            cumgamm = gammp(ac,x/bc) - pc
        endif
*
        end
*  #] cumgamm:
*  #[ invcumgamm:
        real function invcumgamm(p,a,b)
*
*       compute the inverse of the cumulative Gamma distribution P(a,x/b)
*       as long as I do not find or make an explicit function just solve
*       the equation.
*
        implicit none
        real p,a,b
        integer i
        real x,x1,x2,tol
        real pc,ac,bc
        common /ccumgamm/ pc,ac,bc
        real cumgamm,zbrent
        external cumgamm,zbrent
*
*       check argument
*
        if ( p.eq.0 ) then
            x = 0
            goto 999
        elseif ( p.lt.0 .or. p.ge.1 ) then
            write(0,*) 'invcumgamm: illegal argument ',p
            x = 3e33
            goto 999
        endif
*
*       parameters for function cumgamm
*
        ac = a
        bc = b
        pc = p
*
*       bracket zero
*
        x1 = 0
        x2 = b
        i = 0
  100   continue
        if ( cumgamm(x2).lt.0 ) then
            x2 = 1.6*x2
            i = i + 1
            if ( i.gt.100 ) then
                write(0,*) 'invcumgamm: error: cannot find root'
                x = 3e33
                goto 999
            endif
            goto 100
        endif
*
*       get root
*
        tol = 1e-5*x2
        x = zbrent(cumgamm,x1,x2,tol)
*
*       finito
*
  999   continue
        invcumgamm = x
        return
        end
*  #] invcumpois:
