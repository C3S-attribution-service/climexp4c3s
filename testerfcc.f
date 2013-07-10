        program test
        implicit none
        integer ier
        real z,erfcc,p
 1      continue
        read *,z
        p = erfcc(z/sqrt(2.))
        print *,'z,p = ',z,p
        call merfi(1-p,z,ier)
        z = z*sqrt(2.)
        print *,'z,ier = ',z,ier
        goto 1
        end
