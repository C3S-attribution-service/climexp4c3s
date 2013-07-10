*       get the factor with which one has to multiply the number 
*       of reals (or integrs) in a direct access file when opening it.
*       (recfac), and the same parameter for single-precision reals (recfa4)
        integer recfac,recfa4
#if defined(sun) || defined(__sun__) || defined (__NeXT__) || defined(linux)
*       but the Sun (with Sun f77 2.0 || g77) and NeXT/Linux (with f2c) in bytes
#ifdef REAL64
        parameter(recfac=8)
#else
        parameter(recfac=4)
#endif
        parameter(recfa4=4)
#else
#if defined(__alpha) || defined(__sgi)
*       the alpha and indies counts record lengths in 4-byte units
#ifdef REAL64
        parameter(recfac=2)
#else
        parameter(recfac=1)
#endif
        parameter(recfa4=1)
#else
*       the Cray measures direct access files in bytes
        parameter(recfac=8)
        parameter(recfa4=4)
#endif
#endif

