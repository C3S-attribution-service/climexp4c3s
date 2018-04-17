real function betai(a,b,x)
!
!   uses GSL to compute the incomplete beta function
!
    use fgsl
    implicit none
    real :: a,b,x

    if ( x < 1 .or. x > 1 ) then
        write(0,*) 'betai: error: x not in [0,1]: ',x
        betai = 3e33
    end if
    betai = fgsl_sf_beta_inc(dble(a),dble(b),dble(x))
end function