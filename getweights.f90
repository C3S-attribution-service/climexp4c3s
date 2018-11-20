subroutine getweights(xyz,xx,wx,nx,xwrap,lwrite)

!   compute wights (length between halfway points) of axis
!   for the points between x1 and x2

    implicit none
    integer,intent(in) :: nx
    real,intent(in) :: xx(nx)
    real,intent(out) :: wx(nx)
    logical,intent(in) :: xwrap,lwrite
    character,intent(in) :: xyz

    integer :: i,ip1,im1
    real :: pi

    if ( lwrite ) then
        print *,'getweights: input ',xyz
        print *,'xwrap = ',xwrap
        if ( nx >= 2 ) then
            print *,'xx    = ',xx(1),xx(2),'...',xx(nx-1),xx(nx)
        endif
    endif
    pi = 4*atan(1d0)
    if ( xyz /= 'x' .and. xwrap ) then
        write(0,*) 'getweights: error: only x can wrap! ',xyz
        call abort
    endif
    if ( nx == 1 ) then
        wx(1) = 1
        return
    endif
    do i=1,nx
        if ( xwrap ) then
            im1 = 1+mod(i+nx-2,nx)
            ip1 = 1+mod(i     ,nx)
            wx(i) = min( &
                abs(xx(ip1)-xx(im1)), &
                abs(xx(ip1)-xx(im1)-360), &
                abs(xx(ip1)-xx(im1)+360))/2
        elseif ( i == 1 ) then
            wx(i) = abs(xx(2)-xx(1))
        elseif ( i == nx ) then
            wx(i) = abs(xx(nx)-xx(nx-1))
        else
            wx(i) = abs(xx(i+1)-xx(i-1))/2
        endif
    enddo
    if ( xyz == 'y' ) then
!       due to round-off errors this can become -\epsilon at the poles.
        do i=1,nx
            wx(i) = max(0.,wx(i)*cos(xx(i)*pi/180))
        enddo
    endif
    if ( lwrite ) then
        print *,'w',xyz,' = '
        print '(10f8.2)',(wx(i),i=1,nx)
    endif
end subroutine getweights

subroutine getxyprop(xx,nx,yy,ny,xrev,yrev,xwrap)
    implicit none
    integer,intent(in) :: nx,ny
    real,intent(in) :: xx(nx),yy(ny)
    logical,intent(out) :: xrev,yrev,xwrap

    if ( nx > 1 ) then
        xrev  = xx(2) < xx(1)
    else
        xrev = .false. 
    endif
    if ( ny > 1 ) then
        yrev  = yy(2) < yy(1)
    else
        yrev = .false. 
    endif
    if ( nx == 1 ) then
        xwrap = .false. 
    elseif ( xrev ) then
        xwrap = abs(2*xx(nx)-xx(nx-1)+360-xx(1)) < 1e-3
    else
        xwrap = abs(2*xx(nx)-xx(nx-1)-360-xx(1)) < 1e-3
    endif
end subroutine getxyprop
