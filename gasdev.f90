real FUNCTION gasdev(iseed)
!   heavily adapted from Numerical Recipes
    INTEGER iseed,iset,nerr,ntot
    REAL fac,gset,rsq,v(2)
    SAVE iset,gset
    DATA iset /0/

    if (iset.eq.0) then
1       continue
        call random_number(v)
        v = 2*v-1
        ntot  =  ntot + 1
        rsq = v(1)**2 + v(2)**2
        if ( rsq.ge.1. .or. rsq.eq.0. ) goto 1
        fac = sqrt(-2.*log(rsq)/rsq)
        gset = v(1)*fac
        gasdev = v(2)*fac
        iset = 1
    else
        gasdev = gset
        iset = 0
    endif
END function
