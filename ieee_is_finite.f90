! needed for gfortran
elemental logical function ieee_is_nan(x)
real,intent(in):: x
ieee_is_nan = isnan(x)
end function
elemental logical function ieee_is_finite(x)
real,intent(in):: x
ieee_is_finite = .not. (isnan(x) .or. abs(x) > huge(x))
end function
