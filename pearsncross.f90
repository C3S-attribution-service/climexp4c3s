subroutine pearsncross(x,y,n,r,prob,z,ax,sxx,ay,syy,sxy,df,ncrossvalidate)

!   compute Pearson's r with cross-validation if requested
!   df = 3e33 - forget about significance...

    implicit none
    integer :: n,ncrossvalidate
    real :: x(n),y(n),r,prob,z,ax,sxx,ay,syy,sxy,df
    integer :: i,j,m
    real :: s,t,xt,yt,betai,arg3
    integer,save :: ierr,nerr
    data ierr /0/
    data nerr /1/

    ax = sum(x)/n
    ay = sum(y)/n
    sxx = 0
    syy = 0
    sxy = 0
    do j=1,n
        if ( ncrossvalidate <= 1 ) then
            xt = x(j) - ax
            yt = y(j) - ay
        else
            s = ax*n
            t = ay*n
            m = n
            do i=max(1,j-ncrossvalidate/2), &
                min(n,j+(ncrossvalidate-1)/2)
                s = s - x(i)
                t = t - y(i)
                m = m - 1
            end do
            if ( m <= 1 ) then
                r = 3e33
                prob = 3e33
                return
            end if
            xt = x(j) - s/m
            yt = y(j) - t/m
        end if
        if ( abs(xt) < 1e-6*abs(ax) ) xt = 0
        if ( abs(yt) < 1e-6*abs(ay) ) yt = 0
        sxx = sxx + xt**2
        syy = syy + yt**2
        sxy = sxy + xt*yt
    end do
    if ( ncrossvalidate == 1 ) then
        sxx = n**2*sxx/(n-1)**2
        syy = n**2*syy/(n-1)**2
        sxy = n**2*sxy/(n-1)**2
    end if
    if ( sxx*syy == 0 ) then
        r = 0
    else
        r = sxy/sqrt(sxx*syy)
    end if
    if ( r > 1 ) r = 1
    if ( r < -1 ) r = -1
    if ( abs(r-1) < 1e-6 ) then
        z = 3e33
        prob = 0
        return
    end if
    z = (1+r)/(1-r)
    if ( z > 0 ) then
        z = 0.5*log(z)
    else
        write(0,*) 'pearsnx: error: arg log<=0: ',z
        write(0,*) '                r = ',r
        z = 3e33
    end if
    if ( df > 1e33 ) then
        prob = 3e33
        return
    end if
    t = (1-r)*(1+r)
    if ( .not. df > 0 ) then
        ierr = ierr + 1
        if ( ierr >= nerr ) then
            nerr = 2*nerr
            write(*,*) 'pearsncross: df<=0: ',df,ierr
            write(0,*) 'pearsncross: df<=0: ',df,ierr
        end if
        prob = 3e33
        return
    end if
    if ( t > 0 ) then
        t = r*sqrt(df/t)
    else
        write(0,*) 'pearsnx: error: arg sqrt<=0: ',t
        write(0,*) '                r = ',r
        t = 3e33
    endif
    arg3 = df/(df+t**2)
    if ( .not. (arg3 >= 0. .and. arg3 <= 1.) ) then
        print *,'pearsnx: error: arg betai = ',df/(df+t**2)
        print *,'                df = ',df
        print *,'                r  = ',r
        print *,'                t  = ',t
        do j=1,n
            print *,j,x(j),x(j)-ax,y(j),y(j)-ay
        enddo
        return
    end if
    prob = betai(df/2,0.5,arg3)
end subroutine pearsncross
