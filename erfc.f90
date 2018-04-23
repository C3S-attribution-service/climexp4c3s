real function erfcc(x)
!
!   wrapper to compute the complementary error function using (F)GSL
!
    use fgsl
    implicit none
    real :: x
    erfcc = fgsl_sf_erfc(dble(x))
end function erfcc
real function erfc(x)
!
!   wrapper to compute the complementary error function using (F)GSL
!
    use fgsl
    implicit none
    real :: x
    erfc = fgsl_sf_erfc(dble(x))
end function erfc
