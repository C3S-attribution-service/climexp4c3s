!       get the factor with which one has to multiply the number 
!       of reals (or integrs) in a direct access file when opening it.
!       (recfac), and the same parameter for single-precision reals (recfa4)
!       this used to be all over the place but has been standardised since
        integer recfac,recfa4
        parameter(recfac=8)
        parameter(recfa4=4)
