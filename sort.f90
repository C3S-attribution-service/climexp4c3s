subroutine nrsort(n,a)
!
!   wrapper for the equivalent (F)GSL routine
!
    use fgsl
    implicit none
    integer,intent(in) :: n
    real,intent(inout) :: a(n)
    integer(fgsl_size_t) :: dim
    real(fgsl_double), allocatable :: aa(:)
    allocate(aa(n))
    aa = a
    dim = n
    call fgsl_sort(aa,1_fgsl_size_t,dim)
    a = aa
    deallocate(aa)
end subroutine nrsort

subroutine nrsort_d(n,a)
!
!   wrapper for the equivalent (F)GSL routine
!
    use fgsl
    implicit none
    integer,intent(in) :: n
    real(kind=8),intent(inout) :: a(n)
    integer(fgsl_size_t) :: dim
    real(fgsl_double), allocatable :: aa(:)
    allocate(aa(n))
    aa = a
    dim = n
    call fgsl_sort(aa,1_fgsl_size_t,dim)
    a = aa
    deallocate(aa)
end subroutine nrsort_d
