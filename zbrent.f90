module ZbrentToGSL
!
!   call (F)GSL simplex minimisation routine, provides an interface identical to Numerical recipes ameba
!   Adapted from FGSL sample code by Camiel, attempting to understand and adjusted to single precision by GJvO
!
    use, intrinsic :: iso_c_binding
    use fgsl
    implicit none
    private

    public :: zbrent

    interface
        function funk_interface(x)
            real, intent(in) :: x
            real :: funk_interface ! because zbrent gets a pointer to a function that returns a real
        end function funk_interface
    end interface

    procedure(funk_interface), pointer :: funk_pointer => Null()

contains

    function funk_wrapper_z(x, params) bind(c)

!       this is what GSL expects: two pointers to fgsl_vector arguments, fgsl_double/c_double return value
!       The function funk passed to zbrent, which funk_pointer points to, expects one real  
!       as arguments and returns a real.
!       Note that params is not used, as the zbrent interface does not allow for these and hence
!       these are passed in common. As soon as backwards compatibility is no longer required this
!       should be cleaned up.

        real(c_double), value      :: x
        type(c_ptr), value         :: params
        real(c_double)             :: funk_wrapper_z
!
        real                       :: r
        integer(fgsl_int)          :: status
        logical,parameter          :: lwrite=.false.

        r = x
        if ( lwrite) print *,'funk_wrapper: calling funk_pointer with arg ',r
        funk_wrapper_z = dble(funk_pointer(r))
        if ( lwrite) print *,'funk_wrapper: funk_pointer returned ',funk_wrapper_z
        return
    end function funk_wrapper_z

    function zbrent(funk,x1in,x2in,ftolin)
        real,intent(in)     :: x1in,x2in,ftolin
        real                :: zbrent

        integer(fgsl_int)    :: iter
        real(fgsl_double)    :: x1,x2,ftol
        real(fgsl_double), target :: fpar(1) ! increase as soon as it replaces the common block
        real,external        :: funk
        interface
            function funk(x)
                real,intent(in) :: x ! data
            end function funk
        end interface

        integer(fgsl_int), parameter   :: ITMAX = 5000

        type(fgsl_root_fsolver)        :: root_fslv
        type(fgsl_function)            :: stdfunc

        real(fgsl_double), parameter   :: eps5 = 1.0d-5
        real(fgsl_double)              :: dmx, ra, xlo, xhi
        type(c_ptr)                    :: yptr,ptr
        integer(fgsl_int)              :: status, i
        logical,parameter              :: lwrite=.false.
        character(kind=fgsl_char,len=fgsl_strmax) :: name
        real(fgsl_double), pointer     :: xptr(:)
        real                           :: xret,z(4)        
!
!       convert to fgsl_double
!
        if ( lwrite ) then
            print *,'zbrent: x1in,x2in,ftolin = ',x1in,x2in,ftolin
        end if
        x1 = x1in
        x2 = x2in
        ftol = ftolin

        funk_pointer => funk

        root_fslv = fgsl_root_fsolver_alloc(fgsl_root_fsolver_brent)
        name      = fgsl_root_fsolver_name(root_fslv)
        if ( lwrite ) print *,'zbrent: name = ',trim(name)
  
        ! initialize function object
        fpar = (/ 0. /) ! currently unused, parameters are still passed in common for backward compatibility
        ptr     = c_loc(fpar)
        if ( lwrite ) print *,'zbrent: calling fgsl_function_init'
        stdfunc  = fgsl_function_init(funk_wrapper_z,ptr)
!
        if ( lwrite ) print *,'zbrent: calling fgsl_root_fsolver_set'
        status = fgsl_root_fsolver_set(root_fslv, stdfunc, x1, x2)
        if (status /= fgsl_success) then
            print *, 'Failed to configure the solver: status=',status
        end if
        name = fgsl_root_fsolver_name(root_fslv)
        if ( lwrite ) print *,'zbrent: name = ',trim(name)

        iter = 0
        do
            iter = iter + 1
            if ( lwrite ) print *,'zbrent: calling fgsl_root_fsolver_iterate ',iter
            status = fgsl_root_fsolver_iterate(root_fslv)
            if ( status /= fgsl_success .or. iter >= ITMAX ) then
                write(0,*) 'zbrent: error or no convergence ',status,iter
                exit
            end if
            ra = fgsl_root_fsolver_root(root_fslv)
            xlo = fgsl_root_fsolver_x_lower(root_fslv)
            xhi = fgsl_root_fsolver_x_upper(root_fslv)
            status = fgsl_root_test_interval(xlo,xhi,0.0_fgsl_double,eps5)
            if ( status == fgsl_success ) then
                if ( lwrite ) print *,'zbrent: converged after ',iter,' iterations'
                exit
            end if
        end do
        call fgsl_root_fsolver_free(root_fslv)
        call fgsl_function_free(stdfunc)
        if ( lwrite ) print *,'zbrent: copying to output'
        zbrent = ra
        if ( lwrite ) then
            print *,'zbrent = ',zbrent
        end if
        return
    end function zbrent

end module ZbrentToGSL
