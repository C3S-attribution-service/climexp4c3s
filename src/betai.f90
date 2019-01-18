real function betai(a,b,x)
!
!   uses GSL to compute the incomplete beta function
!
    use fgsl
    implicit none
    real :: a,b,x
    real(fgsl_double) :: aa,bb,xx

    if ( x < 0 .or. x > 1 ) then
        write(0,*) 'betai: error: x not in [0,1]: ',x
        betai = 3e33
    end if
    aa = a
    bb = b
    xx = x
    betai = fgsl_sf_beta_inc(aa,bb,xx)
end function