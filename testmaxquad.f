        program testmaxquad
        real x(10),y(10)
        double precision ranf
        external ranf
    1   do i=1,10
            x(i) = i
            y(i) = ranf(i)
        enddo
        call maxquad(xmax,ymax,x,y,10)
        read *,i
        goto 1
        
        end
