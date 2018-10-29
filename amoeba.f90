module AmoebaToGSL
!
!   call (F)GSL simplex minimisation routine, provides an interface similar to Numerical recipes ameba
!   with two big differences: 
!   - initialisation is included
!   - parameters are passed as arguments and not in common
!   Adapted from FGSL sample code by Camiel, attempting to understand by GJvO
!
    use, intrinsic :: iso_c_binding
    use fgsl
    implicit none
    private

    public :: amoeba

    interface
        function funk_interface(x,p)
            use fgsl
            real(fgsl_double), intent(in) :: x(:) ! data
            real(fgsl_double), intent(in) :: p(:) ! parameters
            real(fgsl_double) :: funk_interface
        end function funk_interface
    end interface

    integer                            :: ndimensions  =  0
     procedure(funk_interface), pointer :: funk_pointer => Null()

contains

    function funk_wrapper(v, p) bind(c)
        type(c_ptr), value :: v, p
        real(c_double)     :: funk_wrapper
!
        type(fgsl_vector)          :: vec
        real(fgsl_double), pointer :: par(:), pvec(:)
        integer(fgsl_int)          :: status
        call fgsl_obj_c_ptr(vec, v)
        call c_f_pointer(p, par, (/ ndimensions /))
        status = fgsl_vector_align(pvec, vec)
        funk_wrapper = funk_pointer(par,pvec)
        return
    end function funk_wrapper


    subroutine amoeba(p, y, mp, np, ndim, ftol, funk, iter)
        integer(fgsl_size_t), intent(in)    :: mp, np, ndim
        integer(fgsl_int),    intent(out)   :: iter
        real(fgsl_double),    intent(inout) :: p(mp,np), y(mp)
        real(fgsl_double),    intent(in)    :: ftol
        interface
            function funk(x, p)
                use fgsl
                real(fgsl_double), intent(in) :: x(:) ! data
                real(fgsl_double), intent(in) :: p(:) ! parameters
                real(fgsl_double) :: funk
            end function funk
        end interface

        integer(fgsl_int), parameter :: ITMAX = 5000

        type(fgsl_multimin_fminimizer) :: mmin_fmin
        type(fgsl_multimin_function)   :: mmin_f
        type(fgsl_vector)              :: xvec, step

        real(fgsl_double), target      :: fpar(ndim), xv(ndim), stepv(ndim)
        real(fgsl_double)              :: dmx
        type(c_ptr)                    :: ptr
        integer(fgsl_int)              :: status, i
        character(kind=fgsl_char,len=fgsl_strmax) :: name
        real(fgsl_double), pointer :: xptr(:)

        funk_pointer => funk
        ndimensions  =  ndim
    
        mmin_fmin = fgsl_multimin_fminimizer_alloc(fgsl_multimin_fminimizer_nmsimplex,ndim)
        name      = fgsl_multimin_fminimizer_name(mmin_fmin)
        call unit_assert_equal('fgsl_multiroot_fsolver_name','nmsimplex',trim(name))
  
        ! initialize function object
        fpar(:) = p(1,1:ndim)
        ptr     = c_loc(fpar)
        mmin_f  = fgsl_multimin_function_init(funk_wrapper,ndim,ptr)
!
        xv(:)  = p(1,1:ndim)
        xvec   = fgsl_vector_init(1.0_fgsl_double)
        status = fgsl_vector_align(xv,ndim,xvec,ndim,0_fgsl_size_t,1_fgsl_size_t)
        if (status /= fgsl_success) then
            print*, 'Failed to construct the starting vector: status=',status
        end if

        do i=1,ndim
            stepv(i) = p(1+i,i) - p(1,i)
        end do
        step   = fgsl_vector_init(1.0_fgsl_double)
        status = fgsl_vector_align(stepv,ndim,step,ndim,0_fgsl_size_t,1_fgsl_size_t)
        if (status /= fgsl_success) then
            print*, 'Failed to construct the step size vector: status=',status
        end if
    
        status = fgsl_multimin_fminimizer_set(mmin_fmin, mmin_f, xvec, step)
        if (status /= fgsl_success) then
            print *, 'Failed to configure the minimizer: status=',status
        end if

        call fgsl_vector_free(xvec)
        call fgsl_vector_free(step)
    
        iter = 0
        do
            iter = iter + 1
            status = fgsl_multimin_fminimizer_iterate(mmin_fmin)
            if (status /= fgsl_success .or. iter > ITMAX) exit

            dmx    = fgsl_multimin_fminimizer_size(mmin_fmin)
            status = fgsl_multimin_test_size(dmx, ftol)
            if (status == fgsl_success) exit
        end do

        xvec   = fgsl_multimin_fminimizer_x(mmin_fmin)
        status = fgsl_vector_align(xptr, xvec)
        if (status /= fgsl_success) then
            print *, 'Failed to retrieve the current best estimate: status=',status
        end if
        p(1,1:ndim) = xptr(:)

        call fgsl_multimin_fminimizer_free(mmin_fmin)
        call fgsl_multimin_function_free(mmin_f)

        funk_pointer => Null()
        ndimensions  =  0
    
        return
    end subroutine amoeba

end module AmoebaToGSL
