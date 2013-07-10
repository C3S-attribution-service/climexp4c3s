        subroutine zerolin(xnul,negpos,x,y,n)
*
*       compute the first zero fitting a straight line through the 
*       2 points around the maximum.  Does not assume equidistant points.
*       if negpos<0 look for descending slope, if >0 ascending, 
*       if 0 don't care.
*
        implicit none
        integer negpos,n
        real xnul,x(n),y(n)
        integer i,imax
        real xm,xp,ym,yp,a,b,xx,yy
        logical lwrite
        parameter (lwrite=.false.)
        if ( lwrite ) then
            print *,'zerolin: input n=',n
            print *,'x = ',x
            print *,'y = ',y
        endif
*
*       safety first
*       
        xnul = 3e33
        if ( n.lt.2 ) then
            write(0,*) 'zerolin: error: array too short',n,'<br>'
            return
        endif
*
*       find first zero
*
        do i=2,n
            if ( y(i).lt.1e33 .and. y(i-1).lt.1e33 .and. (
     +           negpos.eq.0 .and. y(i)*y(i-1).lt.0 .or.
     +           negpos.lt.0 .and. y(i-1).gt.0 .and. y(i).lt.0 .or.
     +           negpos.gt.0 .and. y(i-1).lt.0 .and. y(i).gt.0 ) ) then
                xnul = (x(i-1)*y(i) - x(i)*y(i-1))/(y(i)-y(i-1))
                if ( lwrite ) print *,'zerolin: found zero at ',xnul
                return
            endif
        enddo
        if ( lwrite ) then
            write(*,*) 'zerolin: no zero crossing found'
        endif
        return
        end
