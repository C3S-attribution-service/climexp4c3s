real function erfcc(x)
!
!   wrapper to compute the complementary error function using (F)GSL
!
    use fgsl
    implicit none
    real :: x
    real(fgsl_double) :: xx
    xx = x
    erfcc = fgsl_sf_erfc(xx)
end function erfcc
real function erfc(x)
!
!   wrapper to compute the complementary error function using (F)GSL
!
    use fgsl
    implicit none
    real :: x
    real(fgsl_double) :: xx
    xx = x
    erfc = fgsl_sf_erfc(xx)
end function erfc
