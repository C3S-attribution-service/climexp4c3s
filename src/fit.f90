subroutine fit(xx,yy,ndata,sig,mwt,a,b,siga,sigb,chi2,q)
!
!   wrapper that uses GSL routines instead of Numerical recipes code
!   fits a+b*x to (x,y)
!
    use fgsl
    implicit none
    integer,intent(in) :: ndata,mwt
    real,intent(in) :: xx(ndata),yy(ndata),sig(ndata)
    real,intent(out) :: a,b,siga,sigb,chi2,q
    integer :: i
    integer(fgsl_size_t) :: dim
    integer(fgsl_int) :: iret
    real(fgsl_double) :: dc0,dc1,dcov00,dcov01,dcov11,dsumsq,da,db,dq
    real(fgsl_double), allocatable :: dxx(:),dyy(:),dww(:)
!
!   special cases cause crashes later on... (Numerical recipes never checks) (nor does GSL)
!
    a = 3e33
    b = 3e33
    siga = 3e33
    sigb = 3e33
    chi2 = 0
    q = 3e33
    if ( ndata == 0 ) then
        return
    end if
    if ( ndata == 1 ) then
        a = yy(1)
        b = 0
        chi2 = 0
        return
    end if
    if ( ndata == 2 ) then
        if ( xx(1) == xx(2) ) then
            if ( yy(1) == yy(2) ) then
                a = yy(1)
                b = 0
            else
                a = 3e33
                b = 3e33
                chi2 = 3e33
            end if
            return
        end if
        a = (yy(2) - yy(1))/(xx(2) - xx(1))
        b = yy(1) - a*xx(1)
        chi2 = 0
        return
    end if
    if ( all(xx == xx(1)) ) then ! GSL crashes on this
        a = 3e33
        b = sum(yy)/ndata
        return
    end if
    
    dim = ndata
    ! GSL uses double precision
    if ( mwt == 0 ) then
        ! unweighted
        allocate(dxx(ndata),dyy(ndata))
        dxx = xx
        dyy = yy
        iret = fgsl_fit_linear(dxx, 1_fgsl_size_t, dyy, 1_fgsl_size_t, &
            dim, dc0, dc1, dcov00, dcov01, dcov11, dsumsq)
        deallocate(dxx,dyy)
    else
        ! weighted
        allocate(dxx(ndata),dyy(ndata),dww(ndata))
        dxx = xx
        dyy = yy
        dww = 1/sig**2
        iret = fgsl_fit_wlinear(dxx, 1_fgsl_size_t, dww, 1_fgsl_size_t, dyy, 1_fgsl_size_t, &
            dim, dc0, dc1, dcov00, dcov01, dcov11, dsumsq)
        deallocate(dxx,dyy,dww)
    end if 
    a = dc0
    b = dc1
    siga = sqrt(dcov00)
    sigb = sqrt(dcov11)
    chi2 = dsumsq
!
!   only the p-value left to compute, Eq 15.2.12
!
    da = (ndata - 2)/dble(2)
    db = dsumsq/2
    dq = fgsl_sf_gamma_inc_Q(da,db)
    q = dq
end subroutine fit