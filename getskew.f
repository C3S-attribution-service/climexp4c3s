        real function getskew(q)
        implicit none
        real q(4)
        integer nmax
        parameter (nmax=105000)
        integer n
        real x(nmax),y(nmax)
        common /cfitfun/ x,y,n
        integer i
	real e,e1,e2,e3

        print *,'getskew: n = ',n
        e1 = 0
        do i=1,n
            e = (y(i) - q(1) - q(2)*x(i)*(1 + q(3)*x(i)))/
     +            sqrt(1 + q(4)*x(i))
	    e1 = e1 + e
	enddo
	e1 = e1/n
        e2 = 0
	e3 = 0
        do i=1,n
            e = (y(i) - q(1) - q(2)*x(i)*(1 + q(3)*x(i)))/
     +            sqrt(1 + q(4)*x(i))
	    e2 = e2 + (e-e1)**2
	    e3 = e3 + (e-e1)**3
	enddo
	e2 = e2/n
	e3 = e3/n
	getskew = e3/(sqrt(e2))**3
        end
