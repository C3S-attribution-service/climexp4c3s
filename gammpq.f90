real function gammq(a,x)
!
!   wrapper around ghe (F)GSL routine to compute the incomplete Gamma function Q(a,x)
!   Note that GSL uses double precision
!
    use fgsl
    implicit none
    real :: a,x
    
    if ( x < 0 .or. a <= 0 ) then
        write(0,*) 'gammq: error: bad arguments in gammq: ',x,a
        call exit(-1)
    end if
    gammq = fgsl_sf_gamma_inc_Q(dble(a),dble(x))
end function gammq

real function gammp(a,x)
!
!   wrapper around ghe (F)GSL routine to compute the incomplete Gamma function P(a,x)
!   Note that GSL uses double precision
!
    use fgsl
    implicit none
    real :: a,x
    
    if ( x < 0 .or. a <= 0 ) then
        write(0,*) 'gammp: error: bad arguments in gammp: ',x,a
        call exit(-1)
    end if
    gammp = fgsl_sf_gamma_inc_P(dble(a),dble(x))
end function gammp

real function gammln(x)
!
!     wrapper fir GSL function
!
    use fgsl
    implicit none
    real :: x
    gammln = fgsl_sf_lngamma(dble(x))
end function gammln

real function factrl(n)
!
!    wrapper for FGSL function
!
    use fgsl
    implicit none
    integer :: n
    factrl = fgsl_sf_fact(n)
end function factrl

