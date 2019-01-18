module BrentToGSL
!
!   call (F)GSL simplex minimisation routine, provides an interface identical to Numerical recipes ameba
!   Adapted from FGSL sample code by Camiel, attempting to understand and adjusted to single precision by GJvO
!
    use, intrinsic :: iso_c_binding
    use fgsl
    implicit none
    private

    public :: brent,mnbrak

    interface
        function funk_interface(x)
            real, intent(in) :: x
            real :: funk_interface ! because brent gets a pointer to a function that returns a real
        end function funk_interface
    end interface

    procedure(funk_interface), pointer :: funk_pointer => Null()

contains

    function funk_wrapper_1(x, params) bind(c)

!       this is what GSL expects: two pointers to fgsl_vector arguments, fgsl_double/c_double return value
!       The function funk passed to brent, which funk_pointer points to, expects one real  
!       as arguments and returns a real.
!       Note that params is not used, as the brent interface does not allow for these and hence
!       these are passed in common. As soon as backwards compatibility is no longer required this
!       should be cleaned up.

        real(c_double), value      :: x
        type(c_ptr), value         :: params
        real(c_double)             :: funk_wrapper_1
!
        real                       :: r
        integer(fgsl_int)          :: status
        logical,parameter          :: lwrite=.false.

        r = x
        if ( lwrite) print *,'funk_wrapper: calling funk_pointer with arg ',r
        funk_wrapper_1 = dble(funk_pointer(r))
        if ( lwrite) print *,'funk_wrapper: funk_pointer returned ',funk_wrapper_1
        return
    end function funk_wrapper_1

    function brent(a,b,c,funk,ftolin,xmin)
        real,intent(in)     :: a,b,c,ftolin
        real,intent(out)    :: xmin
        real                :: brent

        integer(fgsl_int)    :: iter
        real(fgsl_double)    :: xlo,xmid,xhi,ftol
        real(fgsl_double), target :: fpar(1) ! increase as soon as it replaces the common block
        real,external        :: funk
        interface
            function funk(x)
                real,intent(in) :: x ! data
            end function funk
        end interface

        integer(fgsl_int), parameter   :: ITMAX = 5000

        type(fgsl_min_fminimizer)      :: min_fslv
        type(fgsl_function)            :: stdfunc

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
            print *,'brent: a,b,c,ftolin = ',a,b,c,ftolin
        end if
        xlo = a
        xmid = b
        xhi = c
        ftol = ftolin

        funk_pointer => funk

        min_fslv = fgsl_min_fminimizer_alloc(fgsl_min_fminimizer_brent)
        name     = fgsl_min_fminimizer_name(min_fslv)
        if ( lwrite ) print *,'brent: name = ',trim(name)
  
        ! initialize function object
        fpar = (/ 0. /) ! currently unused, parameters are still passed in common for backward compatibility
        ptr     = c_loc(fpar)
        if ( lwrite ) print *,'brent: calling fgsl_function_init'
        stdfunc  = fgsl_function_init(funk_wrapper_1,ptr)
!
        if ( lwrite ) print *,'brent: calling fgsl_min_fminimizer_set'
        status = fgsl_min_fminimizer_set(min_fslv, stdfunc, xmid, xlo, xhi)
        if (status /= fgsl_success) then
            print *, 'Failed to configure the solver: status=',status
        end if
        name = fgsl_min_fminimizer_name(min_fslv)
        if ( lwrite ) print *,'brent: name = ',trim(name)

        iter = 0
        do
            iter = iter + 1
            if ( lwrite ) print *,'brent: calling fgsl_min_fminimizer_iterate ',iter
            status = fgsl_min_fminimizer_iterate(min_fslv)
            if ( status /= fgsl_success .or. iter >= ITMAX ) then
                write(0,*) 'brent: error or no convergence ',status,iter
                exit
            end if
            xmid = fgsl_min_fminimizer_x_minimum(min_fslv)
            xlo = fgsl_min_fminimizer_x_lower(min_fslv)
            xhi = fgsl_min_fminimizer_x_upper(min_fslv)
            status = fgsl_root_test_interval(xlo,xhi,0.0_fgsl_double,ftol)
            if ( status == fgsl_success ) then
                if ( lwrite ) print *,'brent: converged after ',iter,' iterations'
                exit
            end if
        end do
        brent = fgsl_min_fminimizer_x_minimum(min_fslv)
        call fgsl_min_fminimizer_free(min_fslv)
        call fgsl_function_free(stdfunc)
        if ( lwrite ) print *,'brent: copying to output'
        xmin = xmid
        if ( lwrite ) then
            print *,'brent = ',brent
        end if
        return
    end function brent

    subroutine mnbrak(a,b,c,fa,fb,fc,func)

!       Bracket a minimum given two points ax,bx
!       Returns three points a,b,c such that fa > fb < fc.
!       Not in GSL so reinvented the wheel.
 
        implicit none
        real,intent(inout) :: a,b
        real,intent(out)   :: c,fa,fb,fc
        integer            :: iter
        real               :: x,xlo,flo,xhi,fhi
        logical,parameter  :: lwrite=.false.
        real,external      :: func
 
        if ( a > b ) then ! not sure this is needed but just in case
            x = a
            a = b
            b = x
        end if
        fa = func(a)
        c = b
        fc = func(c)
        b = (a+c)/2
        fb = func(b)
    
        iter = 0
        do
            if ( lwrite ) then
                print *,'mnbrak: a,b,c = ',a,b,c
                print *,'     fa,fb,fc = ',fa,fb,fc
            end if
            iter = iter + 1
            if ( fb < fa .and. fb < fc ) exit
            if ( fa < fb .and. fb < fc ) then ! slope to the left, extend to the left
                x = a - (c-a)
                b = a
                fb = fa
                a = x
                fa = func(a)
            else if ( fa > fb .and. fb > fc ) then ! same to the right
                x = c + (c-a)
                b = c
                fb = fc
                c = x
                fc = func(c)
            else ! maximum in range, not very useful
                xlo = a - 10*(c-a)
                flo = func(xlo)
                xhi = c + 10*(c-a)
                fhi = func(xhi)
                if ( fhi < flo ) then ! steeper to the right, let's look to the left
                    x = a - (c-a)
                    b = a
                    fb = fa
                    a = x
                    fa = func(a)
                else
                    x = c + (c-a)
                    b = c
                    fb = fc
                    c = x
                    fc = func(c)
                end if            
            end if
            if ( abs(a) > 1e10 .or. abs(c) > 1e10 ) then
                write(0,*) 'mnbrak: cannot find minimum',a,b,c,iter
                return
            end if
        end do

    end subroutine mnbrak

end module brentToGSL

