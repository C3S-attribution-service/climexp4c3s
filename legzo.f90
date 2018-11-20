SUBROUTINE LEGZO(n,x,w)
    IMPLICIT NONE
    INTEGER, INTENT(in)  :: n
    REAL,    INTENT(out) :: x(n), w(n)
!
!     =========================================================
!     Purpose : Compute the zeros of Legendre polynomial Pn(x)
!               in the interval [-1,1], and the corresponding
!               weighting coefficients for Gauss-Legendre
!               integration
!     Input :   n    --- Order of the Legendre polynomial
!     Output:   X(n) --- Zeros of the Legendre polynomial
!               W(n) --- Corresponding weighting coefficients
!     =========================================================
!    
    INTEGER :: i, j, n0, nr, k
    REAL*8  :: z, z0, p, pd, pf, q, f0, f1, fd, wp, gd, Pi
    REAL*8, PARAMETER :: Epsilon = 3.0e-16

    Pi = 4.0*ATan(1.0)
    n0 = (n + 1)/2

    x(:) = 0.0; w(:) = 0.0;
    DO nr=1,n0
        z  = Cos(Pi*(nr - 0.25)/n)
        DO
            z0 = z
            p  = 1.0
            DO i=1,nr-1
                p = p*(z - x(i))
            END DO
            f0 = 1.0
            IF (nr == n0 .AND. n /= 2*Int(n/2)) z = 0.0
            f1 = z
            DO k=2,n
                pf = (2.0 - 1.0/k)*z*f1 - (1.0 - 1.0/k)*f0
                pd = k*(f1 - z*pf)/(1.0 - z*z)
                f0 = f1
                f1 = pf
            END DO

            IF (z == 0.0) EXIT

!!!         if ( p.eq.0 ) then
!!!             write(0,*) 'legzo: error: p = 0, pf = ',pf
!!!             call abort
!!!         end if
            fd = pf/p
            q  = 0.0
            DO i=1,nr
                wp = 1.0
                DO j=1,nr
                    IF (j /= i) wp = wp*(z - x(j))
                END DO
                q = q + wp
            END DO
            gd = (pd - q*fd)
!!!         if ( p.eq.0 ) then
!!!             write(0,*) 'legzo: error: p = 0, gd = ',gd
!!!             call abort
!!!         end if
            gd = gd/p
!!!         if ( gd.eq.0 ) then
!!!             write(0,*) 'legzo: error: gd = 0, fd = ',fd
!!!             call abort
!!!         end if
            z  = z - fd/gd
    
            if (Abs(z - z0) <= Abs(z)*Epsilon) EXIT
        END DO

        x(    nr) =  z
        x(n+1-nr) = -z
        w(    nr) = 2.0/((1.0 - z*z)*pd*pd)
        w(n+1-nr) = w(nr)
    END DO
    RETURN
END SUBROUTINE LEGZO
