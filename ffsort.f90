subroutine ffsort(a,ii,n)
!
!   call GSL routine, do not try to do everything myself anymore
!
    use fgsl
    implicit none
    integer :: n,ii(n)
    real :: a(n)
    integer :: i
    integer(fgsl_size_t) :: dim
    integer(fgsl_size_t), allocatable :: p(:)
    real(fgsl_double), allocatable :: aa(:)

    dim = n
    allocate(aa(n),p(n))
    aa = a
    call fgsl_sort_index(p,aa,1_fgsl_size_t,dim)
    ii = p+1 ! the routine returns C offset-0 indices rather than Fortran offset-1
    deallocate(aa,p)
end subroutine ffsort
subroutine ffsort_d(a,ii,n)
!
!   call GSL routine, do not try to do everything myself anymore
!
    use fgsl
    implicit none
    integer :: n,ii(n)
    real(kind=8) :: a(n)
    integer(fgsl_size_t) :: dim
    integer(fgsl_size_t), allocatable :: p(:)
    real(fgsl_double), allocatable :: aa(:)

    dim = n
    allocate(aa(n),p(n))
    aa = a
    call fgsl_sort_index(p,aa,1_fgsl_size_t,dim)
    ii = p+1 ! the routine returns C offset-0 indices rather than Fortran offset-1
    deallocate(aa,p)
end subroutine ffsort_d
