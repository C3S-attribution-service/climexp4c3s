subroutine hunt(xx,n,x,j)

!   Numerical recipes compatible wrapper of GSL routine
!   Finds j such that x is between xx(j) and xx(j+1)

    use fgsl
    implicit none
    integer,intent(in)  :: n
    integer,intent(out) :: j
    real,intent(in)     :: xx(n),x

    integer(fgsl_size_t)          :: index_lo,index_hi
    real(fgsl_double)             :: xd
    real(fgsl_double),allocatable :: xa(:)

!   first special cases

    if ( x < xx(1) ) then
        j = 0
        return
    else if ( x > xx(n) ) then
        j = n
        return
    end if

    allocate(xa(n))
    xa = xx
    xd = x
    index_lo = 1
    index_hi = n
    j = fgsl_interp_bsearch(xa, xd, index_lo, index_hi)
    deallocate(xa)
end subroutine hunt
