real function gammq(a,x)
!
!   wrapper around ghe (F)GSL routine to compute the incomplete Gamma function Q(a,x)
!   Note that GSL uses double precision
!
    use fgsl
    implicit none
    real :: a,x
    real(fgsl_double) :: aa,xx
    
    if ( x < 0 .or. a <= 0 ) then
        write(0,*) 'gammq: error: bad arguments in gammq: ',x,a
        call exit(-1)
    end if
    aa = a
    xx = x
    gammq = fgsl_sf_gamma_inc_Q(aa,xx)
end function gammq

real function gammp(a,x)
!
!   wrapper around ghe (F)GSL routine to compute the incomplete Gamma function P(a,x)
!   Note that GSL uses double precision
!
    use fgsl
    implicit none
    real :: a,x
    real(fgsl_double) :: aa,xx
    
    if ( x < 0 .or. a <= 0 ) then
        write(0,*) 'gammp: error: bad arguments in gammp: ',x,a
        call exit(-1)
    end if
    aa = a
    xx = x
    gammp = fgsl_sf_gamma_inc_P(dble(a),dble(x))
end function gammp

real function gammln(x)
!
!     wrapper fir GSL function
!
    use fgsl
    implicit none
    real :: x
    real(fgsl_double) :: xx
    xx = x
    gammln = fgsl_sf_lngamma(xx)
end function gammln

real function factrl(n)
!
!    wrapper for FGSL function
!
    use fgsl
    implicit none
    integer :: n
    integer(fgsl_int) :: nn
    nn = n
    factrl = fgsl_sf_fact(nn)
end function factrl

real function factln(n)
!
!    wrapper for FGSL function
!
    use fgsl
    implicit none
    integer :: n
    integer(fgsl_int) :: nn
    nn = n
    factln = fgsl_sf_lnfact(nn)
end function factln

