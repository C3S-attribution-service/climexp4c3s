module AmoebaToGSL
!
!   call (F)GSL simplex minimisation routine, provides an interface identical to Numerical recipes amoeba
!   Adapted from FGSL sample code by Camiel, attempting to understand and adjusted to single precision by GJvO
!
    use, intrinsic :: iso_c_binding
    use fgsl
    implicit none
    private

    public :: amoeba

    interface
        function funk_interface(x)
            real, intent(in) :: x(10) ! old-school Fortran array, if I specify x(:) I get an f90-array
            real :: funk_interface ! because amoeba gets a pointer to a function that returns a real
        end function funk_interface
    end interface

    integer                            :: ndimensions  =  0
    procedure(funk_interface), pointer :: funk_pointer => Null()

contains

    function funk_wrapper(x, params) bind(c)
        
!       this is what GSL expects: two pointers to fgsl_vector arguments, fgsl_double/c_double return value
!       The function funk passed to amoeba, which funk_pointer points to, expects one real array 
!       as arguments and returns a real.
!       Note that params is not used, as the amoeba interface does not allow for these and hence
!       these are passed in common. As soon as backwards compatibility is no longer required this
!       should be cleaned up.

        type(c_ptr), value :: x, params
        real(c_double)     :: funk_wrapper
!
        type(fgsl_vector)          :: vec
        real(fgsl_double), pointer :: pvec(:)
        real,allocatable           :: rr(:)
        integer(fgsl_int)          :: status
        logical,parameter          :: lwrite=.false.
        
        call fgsl_obj_c_ptr(vec, x)
        status = fgsl_vector_align(pvec, vec)
        if ( lwrite ) print *,'funk_wrapper: pvec = ',pvec(1:ndimensions)
        allocate(rr(ndimensions))
        rr(1:ndimensions) = pvec(1:ndimensions)
        funk_wrapper = dble(funk_pointer(rr))
        if ( lwrite) print *,'funk_wrapper: funk_pointer returned ',funk_wrapper
        deallocate(rr)
        return
    end function funk_wrapper

    subroutine amoeba(p, y, mpin, npin, ndimin, ftolin, funk, iterout)
        integer,intent(in)  :: mpin,npin,ndimin
        integer,intent(out) :: iterout
        real,intent(inout)  :: p(mpin,npin)
        real,intent(out)    :: y(mpin) ! input values are not used
        real,intent(in)     :: ftolin

        integer(fgsl_size_t) :: mp, np, ndim
        integer(fgsl_int)    :: iter
        real(fgsl_double),allocatable,target :: pp(:,:), yy(:)
        real(fgsl_double)    :: ftol
        real,external        :: funk
        interface
            function funk(x)
                real,intent(in) :: x(10) ! data
            end function funk
        end interface

        integer(fgsl_int), parameter :: ITMAX = 5000

        type(fgsl_multimin_fminimizer) :: mmin_fmin
        type(fgsl_multimin_function)   :: mmin_f
        type(fgsl_vector)              :: xvec, step

        real(fgsl_double), target      :: fpar(ndimin), xv(ndimin), stepv(ndimin)
        real(fgsl_double)              :: dmx
        type(c_ptr)                    :: yptr,ptr
        integer(fgsl_int)              :: status, i
        logical,parameter              :: lwrite=.false.
        character(kind=fgsl_char,len=fgsl_strmax) :: name
        real(fgsl_double), pointer     :: xptr(:)
        real                           :: xret,z(10)        
!
!       convert to fgsl_size_t/fgsl_double
!
        if ( lwrite ) then
            print *,'amoeba: mpin,npin,ndimin,ftolin = ',mpin, npin, ndimin, ftolin
            print *,'        p(1,1:2) = ',p(1,1),p(1,2)
        end if
        mp = mpin
        np = npin
        ndim = ndimin
        ftol = ftolin
        allocate(pp(mp,np),yy(mp))
        pp = p

        funk_pointer => funk
        ndimensions  =  ndim
        if ( ndim > 10 ) then
            write(0,*) 'amoeba: error: arrays too small, please increase size to ',ndim
            call exit(-1)
        end if

        mmin_fmin = fgsl_multimin_fminimizer_alloc(fgsl_multimin_fminimizer_nmsimplex,ndim)
        name      = fgsl_multimin_fminimizer_name(mmin_fmin)
        if ( lwrite ) print *,'amoeba: name = ',trim(name)
  
        ! initialize function object
        fpar(:) = pp(1,1:ndim)
        ptr     = c_loc(fpar)
        mmin_f  = fgsl_multimin_function_init(funk_wrapper,ndim,ptr)
!
        xv(:)  = pp(1,1:ndim)
        xvec   = fgsl_vector_init(1.0_fgsl_double)
        status = fgsl_vector_align(xv,ndim,xvec,ndim,0_fgsl_size_t,1_fgsl_size_t)
        if (status /= fgsl_success) then
            print*, 'Failed to construct the starting vector: status=',status
        end if

        do i=1,ndim
            stepv(i) = pp(1+i,i) - pp(1,i)
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
            if (status /= fgsl_success .or. iter >= ITMAX) exit

            dmx    = fgsl_multimin_fminimizer_size(mmin_fmin)
            status = fgsl_multimin_test_size(dmx, ftol)
            if (status == fgsl_success) then
                if ( lwrite ) print *,'amoeba: converged after ',iter,' iterations'
                exit
            end if
        end do

        xvec   = fgsl_multimin_fminimizer_x(mmin_fmin)
        status = fgsl_vector_align(xptr, xvec)
        if (status /= fgsl_success) then
            print *, 'Failed to retrieve the current best estimate: status=',status
        end if
        pp(1,1:ndim) = xptr(:)

        call fgsl_multimin_fminimizer_free(mmin_fmin)
        call fgsl_multimin_function_free(mmin_f)

        funk_pointer => Null()
        ndimensions  =  0

        iterout = iter
        
        p(1,1:ndim) = pp(1,1:ndim)
        do i=2,mpin
            p(i,1:ndim) = p(1,1:ndim)
        end do
        do i=1,ndim
            z(i) = p(1,i)
        end do
        y(1) = funk(z)
        do i=2,mpin
            y(i) = y(1)
        end do
        deallocate(pp,yy)
        if ( lwrite ) then
            print *,'amoeba: p = ',p(1,1:ndimin)
            print *,'amoeba: y = ',y(1)
        end if
        return
    end subroutine amoeba

end module AmoebaToGSL
