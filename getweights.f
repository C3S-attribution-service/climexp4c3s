        subroutine getweights(xyz,xx,wx,nx,xwrap,lwrite)
*
*       compute wights (length between halfway points) of axis
*       for the points between x1 and x2
*
        implicit none
        integer nx
        real xx(nx),wx(nx)
        logical xwrap,lwrite
        character*1 xyz
*
        integer i,ip1,im1
        real pi
*
        if ( lwrite ) then
            print *,'getweights: input ',xyz
            print *,'xwrap = ',xwrap
            if ( nx.ge.2 ) then
                print *,'xx    = ',xx(1),xx(2),'...',xx(nx-1),xx(nx)
            endif
        endif
        pi = 4*atan(1d0)
        if ( xyz.ne.'x' .and. xwrap ) then
            write(0,*) 'getweights: error: only x can wrap! ',xyz
            call abort
        endif
        if ( nx.eq.1 ) then
            wx(1) = 1
            return
        endif
        do i=1,nx
            if ( xwrap ) then
                im1 = 1+mod(i+nx-2,nx)
                ip1 = 1+mod(i     ,nx)
                wx(i) = min(
     +                abs(xx(ip1)-xx(im1)),
     +                abs(xx(ip1)-xx(im1)-360),
     +                abs(xx(ip1)-xx(im1)+360))/2
            elseif ( i.eq.1 ) then
                wx(i) = abs(xx(2)-xx(1))
            elseif ( i.eq.nx ) then
                wx(i) = abs(xx(nx)-xx(nx-1))
            else
                wx(i) = abs(xx(i+1)-xx(i-1))/2
            endif
        enddo
        if ( xyz.eq.'y' ) then
*       due to round-off errors this can become -\epsilon at the poles.
            do i=1,nx
                wx(i) = max(0.,wx(i)*cos(xx(i)*pi/180))
            enddo
        endif
        if ( lwrite ) then
            print *,'w',xyz,' = '
            print '(10f8.2)',(wx(i),i=1,nx)
        endif
        end

        subroutine getxyprop(xx,nx,yy,ny,xrev,yrev,xwrap)
        implicit none
        integer nx,ny
        real xx(nx),yy(ny)
        logical xrev,yrev,xwrap
*
        if ( nx.gt.1 ) then
            xrev  = xx(2).lt.xx(1)
        else
            xrev = .FALSE.
        endif
        if ( ny.gt.1 ) then
            yrev  = yy(2).lt.yy(1)
        else
            yrev = .FALSE.
        endif
        if ( nx.eq.1 ) then
            xwrap = .FALSE.
        elseif ( xrev ) then
            xwrap = abs(2*xx(nx)-xx(nx-1)+360-xx(1)).lt.1e-3
        else
            xwrap = abs(2*xx(nx)-xx(nx-1)-360-xx(1)).lt.1e-3
        endif
        end
