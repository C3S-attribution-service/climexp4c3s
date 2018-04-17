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