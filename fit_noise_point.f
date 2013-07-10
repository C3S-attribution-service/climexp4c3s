        subroutine fit_noise_point(x,y,n,p,iter,e2)
*
*       fit to full noise model
*
        implicit none
        integer n,iter
        real x(n),y(n),p(5,4),e2
        integer i,j
        real q(4),f(5),tol,ff
        real fitfun
        external fitfun
*
***        print *,'fit_noise_point: n     = ',n
***        print *,'fit_noise_point: p(in) = ',(p(1,j),j=1,4)
        do i=2,5
            do j=1,4
                p(i,j) = p(1,j)
            enddo
        enddo
        p(2,1) = p(2,1)*1.2 + 0.1
        p(3,2) = p(3,2)*1.2 + 0.1
        p(4,3) = p(4,3)*1.2 + 0.1
        p(5,4) = p(5,4)*1.2 + 0.1
        ff = 0
        do i=1,5
            do j=1,4
                q(j) = p(i,j)
            enddo
            f(i) = fitfun(q)
            ff = max(ff,abs(f(i)))
        enddo
        tol = 1e-4
        call amoeba(p,f,5,4,4,tol,fitfun,iter)
        e2 = f(1)
        end

        real function fitfun(q)
        implicit none
        real q(4)
        integer nmax
        parameter (nmax=105000)
        integer n
        real x(nmax),y(nmax)
        common /cfitfun/ x,y,n
        integer i
	real e,eps,s
        eps=1e-6

***        print *,'n = ',n
        fitfun = 0
        s = 0
        do i=1,n
            e = (y(i) - q(1) - q(2)*x(i)*(1 + q(3)*x(i)))/
     +            sqrt(abs(1 + q(4)*x(i))+eps)
	    fitfun = fitfun + e**2
            s = s + x(i)**2
	enddo
        s = s/n
        fitfun = fitfun/n
        s = q(4)*sqrt(s)
        if ( s.lt.-0.5 ) fitfun = fitfun + (s+0.5)**2
        if ( s.gt. 1.0 ) fitfun = fitfun + (s-1.0)**2
        fitfun = fitfun
***        print '(5f16.8)',fitfun,q(1),q(2),q(3),q(4)
        end

        real function getskew(q)
        implicit none
        real q(4)
        integer nmax
        parameter (nmax=105000)
        integer n
        real x(nmax),y(nmax)
        common /cfitfun/ x,y,n
        integer i
	real e,e1,e2,e3,eps,emax,yy
        eps=1e-6
        emax = 0

***        print *,'getskew: n = ',n
        e1 = 0
        do i=1,n
            e = (y(i) - q(1) - q(2)*x(i)*(1 + q(3)*x(i)))/
     +            sqrt(abs(1 + q(4)*x(i))+eps)
	    e1 = e1 + e
	enddo
	e1 = e1/n
        e2 = 0
	e3 = 0
        do i=1,n
            yy = (q(1) + q(2)*x(i)*(1 + q(3)*x(i)))/
     +            sqrt(abs(1 + q(4)*x(i))+eps)
            e = y(i) - yy
***            if ( abs(e).gt.abs(emax) ) then
***                emax = e
***                print *,i,e,x(i),y(i),yy
***            endif
	    e2 = e2 + (e-e1)**2
	    e3 = e3 + (e-e1)**3
	enddo
	e2 = e2/n
	e3 = e3/n
	getskew = e3/(sqrt(e2))**3
        end


