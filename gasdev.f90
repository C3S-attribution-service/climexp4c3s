real function gasdev(iseed)
!
!   generates Gaussian-distributed random numbers following Wikipedia
!
    integer,intent(in) :: iseed ! not used
    logical,save :: lgenerate = .true.
    real :: u(2)
    real,save :: z0,z1
    real,parameter :: two_pi = 2*3.14159265358979323846

    if ( lgenerate ) then
        lgenerate = .false.
        u = 0
        do while ( u(1) == 0 )
            call random_number(u)
        end do
        z0 = sqrt(-2*log(u(1))) * cos(two_pi*u(2))
        z1 = sqrt(-2*log(u(1))) * sin(two_pi*u(2))
        gasdev = z0
    else
        lgenerate = .true.
        gasdev = z1
    end if
end function gasdev
