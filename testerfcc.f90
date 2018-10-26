program test
    implicit none
    integer :: ier
    real :: z,erfc,erfcc,p,pi,p1,serfi,serfci
    pi = 4*atan(1.)
1   continue
    read *,z
    p = erfc(z/sqrt(2.))
    print *,'z,p = ',z,p
    z = serfci(p)
    z = z*sqrt(2.)
    print *,'z,p = ',z,p
    goto 1
end program test
