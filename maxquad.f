        subroutine maxquad(xmax,ymax,x,y,n)
*
*       compute the maximum of a peak by fitting a parabola through the 
*       3 points around the maximum.  Does not assume equidistant points.
*
        implicit none
        integer n
        real xmax,ymax,x(n),y(n)
        integer i,imax
        real xm,xp,ym,yp,a,b,xx,yy
        logical lwrite
        parameter (lwrite=.false.)
        if ( lwrite ) then
            print *,'maxquad: input n=',n
            print *,'x = ',x
            print *,'y = ',y
        endif
*
*       safety first
*       
        xmax = 3e33
        ymax = 3e33
        if ( n.lt.3 ) then
            write(0,*) 'maxquad: error: array too short',n,'<br>'
            return
        endif
*
*       find max not on edge
*
        ymax = -3e33
        imax = 0
        do i=2,n-1
            if ( y(i).lt.1e33 .and. y(i).gt.ymax .and. 
     +           y(i-1).lt.1e33 .and. y(i).gt.y(i-1) .and.
     +           y(i+1).lt.1e33 .and. y(i).gt.y(i+1) ) then
                imax = i
                ymax = y(i)
            endif
        enddo
        if ( imax.eq.0 ) then
*            write(*,*) 'maxquad: no peak found<br>'
*            write(0,*) (x(i),i=1,n),'<br>'
*            write(0,*) (y(i),i=1,n),'<br>'
            return
        endif
*
*       reduced variables
*
        xm = x(imax-1) - x(imax)
        xp = x(imax+1) - x(imax)
        ym = y(imax-1) - y(imax)
        yp = y(imax+1) - y(imax)
*       
*       solve 2 equations with 2 unknowns
*       
        a = ( xp   *ym - xm*   yp)/((xm-xp)*xm*xp)
        b = (-xp**2*ym + xm**2*yp)/((xm-xp)*xm*xp)
        if ( a.ge.0 ) then
            write(0,*) 'maxquad: mysterious error: a >= 0: ',a,imax
     +            ,'<br>'
            write(0,*) (x(i),i=1,n),'<br>'
            write(0,*) (y(i),i=1,n),'<br>'
        endif
        xmax = x(imax) - b/(2*a)
        ymax = y(imax) - b**2/(4*a)
*       
*       debug
*       
        if ( lwrite ) then
            print *,xmax,ymax
*            print *
*            do i=1,101
*                xx = x(imax-1)+(i-1)*(x(imax+1)-x(imax-1))/100
*                yy = y(imax) + a*(xx-x(imax))**2 + b*(xx-x(imax))
*                print *,xx,yy
*            enddo
*            print *
            do i=-1,1
                print *,x(imax+i),y(imax+i)
            enddo
        endif
*
        end
