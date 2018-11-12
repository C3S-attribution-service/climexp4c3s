subroutine polint(xx,yy,n,x,y,erry)
!
!   wrapper with Numerical Recipes calling conventions to GSL interpolation routine
!   it would be better to rewrite the loop calling this to take advantage of the acceleration...
!
    use fgsl
    integer,intent(in) :: n
    real,intent(in)    :: xx(n),yy(n),x
    real,intent(out)   :: y,erry
    
    real(fgsl_double),allocatable :: dxx(:),dyy(:)
    real(fgsl_double)       :: dx,dy,derry
    type(fgsl_interp)       :: finterp
    type(fgsl_interp_accel) :: faccel
    integer(fgsl_size_t)    :: ndim,status
    integer                 :: i
    logical                 :: lwrite=.false.
    character               :: name*100
!
!   first a few special cases
!   in contrast to the Numerical Recipes routine we do not do extrapolation
!
    if ( x < xx(1) ) then
        y = yy(1)
        dy = 3e33
        return
    else if ( x > xx(n) ) then
        y = yy(n)
        dy = 3e33
        return
    end if

    allocate(dxx(n),dyy(n))
    ndim = n
    dxx = xx
    dyy = yy
    dx = x
    if ( n < 5 ) then
        finterp = fgsl_interp_alloc(fgsl_interp_polynomial,ndim)
    else
        finterp = fgsl_interp_alloc(fgsl_interp_cspline,ndim)
    end if
    faccel = fgsl_interp_accel_alloc()
    name = fgsl_interp_name(finterp)
    if ( lwrite ) print *,'polint: using ',trim(name)
    status = fgsl_interp_init(finterp,dxx,dyy)
    dy = fgsl_interp_eval(finterp,dxx,dyy,dx,faccel)
    derry = 3e33 ! not yet used in the calling routine
    call fgsl_interp_free(finterp)
    call fgsl_interp_accel_free(faccel)

    y = dy
    erry = derry
    deallocate(dxx,dyy)

end subroutine polint
