*  #[ fitpoi:
        subroutine fitpoi(xx,ntot,mean,a)
*       
*       a fit a Poisson distribution to the data
*       
        implicit none
*        
        integer ntot
        real xx(ntot),mean,a
*
        integer i,nx
        real ax,bx,cx,fa,fb,fc,tol
*
        integer nmax,ncur
        parameter(nmax=10000000)
        real data(nmax),restrain
        logical llwrite
        common /fitdata1/ data
        common /fitdata2/ restrain,ncur,llwrite
*
        real llpoisson
        external llpoisson
*
*       check input
*       
        do i=1,ntot
            if ( abs(xx(i)-nint(xx(i))).gt.0.01 ) then
                write(0,*) 'histogram: error: cannot fit Poisson to '//
     +                'non-integer data: ',xx(i)
                call abort
            endif
        enddo
*
*       copy to common for routine llpoisson
*       
        ncur = ntot
        do i=1,ncur
            data(i) = xx(i)
        enddo
*
*       fit, using Numerical Recipes routines
*
        ax = mean/2
        bx = mean
        call mnbrak(ax,bx,cx,fa,fb,fc,llpoisson)
        tol = 1e-6
        call brent(ax,bx,cx,llpoisson,tol,a)
*
*       output
*
        print '(a)','# Fitted to Poisson distribution'
        print '(a)','# p(n) = mu^n*exp(-mu)/n! with'
        print '(a,f16.3,a)','# mu = ',a
        end
*  #] fitpoi:
*  #[ llpoisson:
        real function llpoisson(a)
*
*       computes the log-likelihood function for a Poisson distribution
*       with parameter a and data in common.
*       
        implicit none
*       
        real a
*       
        integer i
        real s
*
        integer nmax,ncur
        parameter(nmax=10000000)
        real data(nmax),restrain
        logical llwrite
        common /fitdata1/ data
        common /fitdata2/ restrain,ncur,llwrite
*
        real gammln
        external gammln
*
        llpoisson = 0
        do i=1,ncur
            llpoisson = llpoisson + 
     +            data(i)*log(a) - a - gammln(1+data(i))
        enddo
*       normalization is not 1 in case of cut-offs
        call poisnorm(a,s)
        llpoisson = llpoisson - ncur*log(s)
*       minimum, not maximum
        llpoisson = -llpoisson
*        
        end
*  #] llpoisson:
*  #[ poisnorm:
        subroutine poisnorm(a,s)
*
*       get the normalization of the poisson function on the interval 
*       minindx,maxindx
*
        implicit none
#include "getopts.inc"
        real a,s
        integer i
        real xk1,xk2
        real poisson,gammp,gammq
        external poisson,gammp,gammq
*
        if ( minindx.gt.0 ) then
            xk1 = nint(minindx+0.5)
            if ( maxindx.lt.2.**31 ) then
                xk2 = nint(maxindx-0.5) + 1
                s = gammq(xk2,a) - gammq(xk1,a)
            else
                s = gammp(xk1,a)
            endif
        else
            if ( maxindx.lt.1e33 ) then
                xk2 = nint(maxindx-0.5) + 1
                s = gammq(xk2,a)
            else
                s = 1
            endif
        endif
***        print *,'poisnorm: norm = ',s
        end
*  #] poisnorm:
*  #[ poisson:
        real function poisson(a,n)
        implicit none
        real a
        integer n
        real factrl
        external factrl
*
        poisson = a**n*exp(-a)/factrl(n)
***        print *,'poisson(',a,n,') = ',poisson
        end
*  #] poisson:
*  #[ cumpois:
        real function cumpois(x)
*
*       compute the cumulative Poisson probability minus some
*       requested frequency.  Parameters are passed in common
*
        implicit none
        real x
        real pc,muc
        common /ccumpois/ pc,muc
        real gammq
        external gammq
*
        if ( x.lt.0 ) then
            cumpois = -pc
        else
            cumpois = gammq(x+1,muc) - pc
        endif
*
        end
*  #] cumpois:
*  #[ invcumpois:
        real function invcumpois(p,mu)
*
*       compute the inverse of the cumulative Poisson distribution Q(k,mu)
*       as long as I do not find or make an explicit function just solve
*       the equation.
*
        implicit none
        real p,mu
        integer i
        real x,x1,x2,tol
        real pc,muc
        common /ccumpois/ pc,muc
        real cumpois,zbrent
        external cumpois,zbrent
*
*       check argument
*
        if ( p.eq.0 ) then
            x = 0
            goto 999
        elseif ( p.lt.0 .or. p.ge.1 ) then
            write(0,*) 'invcumpois: illegal argument ',p
            x = 3e33
            goto 999
        endif
*
*       parameters for function cumpois
*
        muc = mu
        pc = p
*
*       bracket zero
*
        x1 = 0
        x2 = mu
        i = 0
  100   continue
        if ( cumpois(x2).lt.0 ) then
            x2 = 1.6*x2
            i = i + 1
            if ( i.gt.100 ) then
                write(0,*) 'invcumpois: error: cannot find root'
                x = 3e33
                goto 999
            endif
            goto 100
        endif
*
*       get root
*
        tol = 1e-5*x2
        x = zbrent(cumpois,x1,x2,tol)
*
*       finito
*
  999   continue
        invcumpois = x
        return
        end
*  #] invcumpois:
