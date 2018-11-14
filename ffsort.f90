subroutine ffsort(a,ii,n)
!
!   call GSL routine, do not try to do everything myself anymore
!
    use fgsl
    implicit none
    integer :: n,ii(n)
    real :: a(n)
    integer(fgsl_size_t) :: dim
    integer(fgsl_size_t), allocatable :: p(:)
    real(fgsl_double), allocatable :: aa(:)

    dim = n
    allocate(aa(n),p(n))
    aa = a
    call fgsl_sort_index(p,aa,1_fgsl_size_t,dim)
    ii = p
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
    ii = p
    deallocate(aa,p)
end subroutine ffsort_d
