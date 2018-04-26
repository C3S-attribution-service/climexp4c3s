subroutine pearsnxx(x,y,n,r,prob,z,ax,sxx,ay,syy,sxy,df)
    implicit none
    integer,intent(in) :: n
    real,intent(in) :: x(n),y(n),df
    real,intent(out) :: r,prob,z,ax,sxx,ay,syy,sxy
    call pearsnxx_gen(x,y,n,r,prob,z,ax,sxx,ay,syy,sxy,df,.true.,0)
end subroutine pearsnxx
subroutine apearsnxx(x,y,n,r,prob,z,ax,sxx,ay,syy,sxy,df)
    implicit none
    integer,intent(in) :: n
    real,intent(in) :: x(n),y(n),df
    real,intent(out) :: r,prob,z,ax,sxx,ay,syy,sxy
    call pearsnxx_gen(x,y,n,r,prob,z,ax,sxx,ay,syy,sxy,df,.false.,0)
end subroutine apearsnxx
subroutine pearsncross(x,y,n,r,prob,z,ax,sxx,ay,syy,sxy,df)
    implicit none
    integer,intent(in) :: n
    real,intent(in) :: x(n),y(n),df,ncrossvalidate
    real,intent(out) :: r,prob,z,ax,sxx,ay,syy,sxy
    call pearsnxx_gen(x,y,n,r,prob,z,ax,sxx,ay,syy,sxy,df,.true.,ncrossvalidate)
end subroutine apearsnxx

subroutine pearsnxx_gen(x,y,n,r,prob,z,ax,sxx,ay,syy,sxy,df,lsubmean,ncrossvalidate)
!
!   re-implementation of the Numerical Recipes subroutine with a few extra parameters
!   using the GSL routines. The extra flags allow for not subtracting the mean or 
!   computing cross-validated results.
!
    use fgsl
    implicit none
    integer,intent(in) :: n,ncrossvalidate
    real,intent(in) :: x(n),y(n),df
    real,intent(out) :: r,prob,z,ax,sxx,ay,syy,sxy
    logical,intent(in) :: lsubmean
    integer :: i
    integer(fgsl_size_t) :: dim
    real,type(fgsl_double) :: dax,day,dsxx,dsxy,dsyy,dr,dz,dt,arg3
    real,type(fgsl_double),allocatable :: xx(:),yy(:)
    
    if ( ncrossvalidate > 0 ) then
        write(0,*) 'presnxx_gen: cross-validation not yet copied from pearsncross.f'
        call exit(-1)
    end if
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
        dsxy = dsxx + (xx(i) - dax)*(yy(i) - day)
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
    ! Fisher z
    dz = (1+dr)/(1-dr)
    z = dz
    if ( z > 0 ) then
        z = log(z)/2
    else
        write(0,*) 'pearsnxx: error: arg log<=0: ',z
        write(0,*) '                 r = ',r
        z = 3e33
    end if
    ! Student t
    dt = (1-dr)*(1+dr)
    if ( df <= 0 ) then
        write(0,*) 'pearsnxx: error: df<=0: ',df
        prob = 3e33
        return
    end if
    dt = dr*sqrt(df/dt)
    arg3 = df/(df+dt**2)
    prob = fgsl_sf_beta_inc(df/2,0.5d0,arg3) ! Wikipedia formula
end subroutine pearsnxx
