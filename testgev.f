        program testgev
        implicit none
        integer i
        real a,b,xi,x,y,z,xyear
        a = 0

 2      continue
        print *,'b,xi,xyear'
        read (*,*,end=1) b,xi,xyear
        z = (1 + xi*(xyear-a)/b)
        print *,'full'
        y = -z**(-1/xi)
        print *,1/(1 - exp(y))
        print *,'xyear large'
        print *,-1/y
        print *,-1/(y + 0.5*y**2)
        print *,'xi~0'
        y = -exp(-(xyear-a)/b + xi/2*((xyear-a)/b)**2)
        print *,1/(1 - exp(y))
        print *,'xi~0, xyear large'
        print *,-1/y
        print *,-1/(y + 0.5*y**2)
        goto 2

 1      continue
        print *,'b,xi = '
        read *,b,xi
        do i=1,10
            if ( mod(i,3).eq.1 ) then
                x = 1 + i/3
            elseif ( mod(i,3).eq.2 ) then
                x = 1 + log10(2.) + i/3
            else
                x = log10(5.) + i/3
            endif
            print *,'x~0, P~0'
            y = -x*log(10.)
            print *,a - b*y + b*xi/2*y**2
            print *,'x~0'
            y = log(-log(1-dble(10)**(dble(-x))))
            print *,a - b*y + b*xi/2*y**2
            print *,'P~0'
            print *,a - b/xi*(1-10.**(x*xi))
            print *,'full'
            print *,a - b/xi*(1-(-log(1-dble(10)**(dble(-x))))
     +                    **(-xi))
        end do
        goto 1
        end
