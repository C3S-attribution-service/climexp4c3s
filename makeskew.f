        subroutine makeskew(x,s,iran)
!
!       generate a random variable with zero mean, unit variance and
!       skewness s by transforming a normal distrbution
!
        implicit none
        real x,s
        integer iran
        integer i
        real x1,x2,tol
        real,save :: sold,a
        real skew
        common /cskew/ skew
        real,external :: gasdev,zbrent,zeroskew
        data sold /3e33/
!
!       generate a random variable from a N(0,1)
        x = gasdev(iran)
        if ( s.eq.0 ) then
            return
        end if
!
!       compute A in transformation y1 = x + A*x
        if ( s.ne.sold) then
            sold = s
            skew  = s
            x1 = -10
            x2 = 10
            tol = 1e-4
            a = zbrent(zeroskew,x1,x2,tol)
        end if
!
!       transform
        x = (x + a*(x**2-1))/sqrt(1+2*a**2)
!
        end

        real function zeroskew(x)
        real x
        real skew
        common /cskew/ skew
        zeroskew = (6*x + 8*x**3)/sqrt((1 + 2*x**2)**3) - skew
!!!        print *,x,skew,zeroskew
        end

