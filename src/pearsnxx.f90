subroutine pearsnxx(x,y,n,r,prob,z,ax,sxx,ay,syy,sxy,df)
    implicit none
    integer,intent(in) :: n
    real,intent(in) :: x(n),y(n),df
    real,intent(out) :: r,prob,z,ax,sxx,ay,syy,sxy
    call pearsnxxgen(x,y,n,r,prob,z,ax,sxx,ay,syy,sxy,df,.true.)
end subroutine pearsnxx
subroutine apearsnxx(x,y,n,r,prob,z,ax,sxx,ay,syy,sxy,df)
    implicit none
    integer,intent(in) :: n
    real,intent(in) :: x(n),y(n),df
    real,intent(out) :: r,prob,z,ax,sxx,ay,syy,sxy
    call pearsnxxgen(x,y,n,r,prob,z,ax,sxx,ay,syy,sxy,df,.true.)
end subroutine apearsnxx

subroutine pearsnxxgen(x,y,n,r,prob,z,ax,sxx,ay,syy,sxy,df,lsubmean)
!
!   re-implementation of the Numerical Recipes subroutine with a few extra parameters
!   using the GSL routines
!
    use fgsl
    implicit none
    integer,intent(in) :: n
    real,intent(in) :: x(n),y(n),df
    real,intent(out) :: r,prob,z,ax,sxx,ay,syy,sxy
    logical,intent(in) :: lsubmean
    integer :: i
    integer(fgsl_size_t) :: dim
    real(fgsl_double) :: dax,day,dsxx,dsxy,dsyy,dr,dz,dt,arg1,arg2,arg3
    real(fgsl_double),allocatable :: xx(:),yy(:)
!
    dim = n
    allocate(xx(n),yy(n))
    xx = x
    yy = y
    if ( lsubmean ) then
        dax = fgsl_stats_mean(xx,1_fgsl_size_t,dim)
        day = fgsl_stats_mean(yy,1_fgsl_size_t,dim)
    else
        dax = 0
        day = 0
    end if
    dsxx = 0
    dsxy = 0
    dsyy = 0
    do i=1,n
        dsxx = dsxx + (xx(i) - dax)**2
        dsxy = dsxy + (xx(i) - dax)*(yy(i) - day)
        dsyy = dsyy + (yy(i) - day)**2
    end do
    deallocate(xx,yy)
    if ( dsxx*dsyy == 0 ) then
	    dr = 0 ! arbitrary choice, least surprise.
    else
        dr = dsxy/sqrt(dsxx*dsyy)
    end if
!
!   I have seen lots of rounding problems...
!
    if ( dr > 1 ) dr = 1
    if ( dr < -1 ) dr = -1
    ! single precision return values
    r = dr
    sxx = sxx
    sxy = dsxy
    syy = dsyy
    if ( abs(r-1).lt.1e-6 ) then
        z = 3e33
        prob = 0
        return
    end if
    dz = (1+dr)/(1-dr)
    if ( dz > 0 ) then
        dz = log(dz)/2
        z = dz
    else
        write(0,*) 'pearsnxx: error: arg log<=0: ',z
        write(0,*) '                 r = ',r
        z = 3e33
    end if
    if ( df > 100 ) then ! normal distribution
        dt = abs(dz)*sqrt(df/2)
        prob = fgsl_sf_erfc(dt)
    else ! Student's t
        dt = (1-dr)*(1+dr)
        if ( df <= 0 ) then
            write(0,*) 'pearsnxx: error: df<=0: ',df
            prob = 3e33
            return
        end if
        dt = dr*sqrt(df/dt)
        arg1 = df/2
        arg2 = 0.5d0
        arg3 = df/(df+dt**2)
        prob = fgsl_sf_beta_inc(arg1,arg2,arg3) ! Wikipedia
    end if
end subroutine pearsnxxgen
