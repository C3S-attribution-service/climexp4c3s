subroutine period(x,y,n,ofac,hifac,px,py,np,nout,jmax,prob)
!
!   wrapper for Numerical Recipes routine to open source routine below
!
    use fgsl
    implicit none
    integer :: n,nout,np,jmax
    real :: x(n),y(n),ofac,hifac,px(np),py(np),prob

    integer,save :: init=0
    integer :: i
    logical :: lregular

    integer(fgsl_size_t) :: ndim
    integer(fgsl_int) :: status
    real :: dx,xmin,xmax
    real(fgsl_double),allocatable :: data(:)
    type(fgsl_fft_real_workspace) :: work
    type(fgsl_fft_real_wavetable) :: wave

    integer :: nt,nw,m,ncpu
    real(kind=8),allocatable :: tt(:),yy(:),sig(:),ww(:),pp(:)
    real(kind=8) :: xxmin,xxmax,pi2 = 6.2831853071795862d0

    character :: string*20
    
!   add test for regular x

    dx = x(2) - x(1)
    lregular = .true.
    do i=2,n-1
        if ( abs(x(i+1) - x(i) - dx) > 0.01*dx ) then
            lregular = .false.
            exit
        end if
    end do
    
    if ( lregular .and. .false. ) then ! normalisation or something else still wrong
        if ( init == 0 ) then
            init = 1
            write(0,*) 'Usng FFT<br>'
        end if
    
!       regular, use FFT

        nout = n/2
        ndim = n
        work = fgsl_fft_real_workspace_alloc(ndim)
        wave = fgsl_fft_real_wavetable_alloc(ndim)
        allocate(data(n))
        data(1:n) = y(1:n)
        status = fgsl_fft_real_transform(data, 1_fgsl_size_t, ndim, wave, work)
        do i=1,nout
            ! see https://www.gnu.org/software/gsl/doc/html/fft.html
            px(i) = i/(dx*(n-1))
            py(i) = (data(2*i-1)**2 + data(2*i)**2)
        end do
        jmax = -1
        prob = 3e33
        
    else
        if ( init == 0 ) then
            init = 1
            write(0,*) 'Using Lomb-Scargle periodogram<br>'
        end if

!       irregular array x, use Lomb-Scargle method
!       copy, generate input parameters

        nt = n
        allocate(tt(nt),yy(nt),sig(nt))
        tt = x
        yy = y
        sig = 1
        call getenv('OMP_NUM_THREADS',string)
        if ( string /= ' ' ) then
            read(string,*) ncpu
        else
            ncpu = 1
        end if
        nout = ofac*hifac*n/2
        nw = nout
        if ( nout > np ) then
            write(0,*) 'period: error: output arrays too small ',nout,np
            call exit(-1)
        end if
        allocate(ww(nout),pp(nout))
        xxmin = minval(tt)
        xxmax = maxval(tt)
        do i=1,nout
            ww(i) = i*pi2/(dble(ofac)*(xxmax-xxmin))
        end do

!       call
    
        call glombscargle(tt,yy,sig,ww,ncpu,pp,m,nt,nw)

!       copy output parameters

        px(1:nout) = ww(1:nout)/pi2
        py(1:nout) = pp(1:nout)
        jmax = -1
        prob = 3e33
    end if
end subroutine period

SUBROUTINE GLOMBSCARGLE(T, Y, SIG, W, NCPU, P, M, NT, NW)
    !!!use omp_lib
    IMPLICIT NONE

    !*  Input arguments
    INTEGER NT, NW, M
    INTEGER :: NCPU
    REAL (KIND=8) T(NT), Y(NT), SIG(NT)
    REAL (KIND=8) W(NW), P(NW)

    !*  Calculate the Generalized Lomb-Scargle (GLS) periodogram as defined by
    !*  Zechmeister, M. & Kürster, M., A&A 496, 577-584, 2009 (ZK2009) 
    !*
    !*  Arguments
    !*  =========
    !* 
    !*  T   (input) REAL (KIND=8) array, dimension (NT)
    !*      Sample times
    !* 
    !*  Y   (input) REAL (KIND=8) array, dimension (NT)
    !*      Measurement values
    !* 
    !* SIG  (input) REAL (KIND=8) array, dimension (NT)
    !*      Measured uncertainties
    !*
    !*  W   (input) REAL (KIND=8) array, dimension (NW)
    !*      Angular frequencies for output periodogram
    !* 
    !* NCPU (input) INTEGER, OPTIONAL
    !*      Number of CPUs to distribute calculation
    !*
    !*  P   (output) REAL (KIND=8) array, dimension (NW)
    !*      Lomb-Scargle periodogram
    !* 
    !*  NT (input) integer
    !*      Dimension of input arrays
    !* 
    !*  NW (output) integer
    !*      Dimension of frequency and output arrays

    !*  Local variables
    INTEGER I, J
    REAL (KIND=8) sig2(NT), ww(NT), th(NT), yh(NT)
    REAL (KIND=8) upow(NW), a(NW), b(NW), off(NW)
    REAL (KIND=8) pi, pi2, pi4
    REAL (KIND=8) Y_, YY, C, S, YC, YS, CCh, CSh, SSh, CC, SS, CS, D
    REAL (KIND=8) x(NT), cosx(NT), sinx(NT), wcosx(NT), wsinx(NT)

    pi2 = 6.2831853071795862d0
    pi4 = 12.566370614359172d0

    th = T

    sig2 = sig*sig
    ww = (1.d0 / sum(1.d0/sig2)) / sig2  ! w of ZK2009, the normalized weights

    Y_ = sum(ww*Y)              ! Eq. (7)
    yh = Y - Y_                 ! Subtract weighted mean

    YY = sum(ww * yh**2)        ! Eq. (16)

    ! Unnormalized power
    upow = 0.d0
    a = 0.d0
    b = 0.d0
    off = 0.d0

    ! If NCPU is present, distribute the calculation using OpenMP
    !!!if (present(NCPU).and.(.not.(NCPU.eq.-1))) then
    !!!    call OMP_SET_NUM_THREADS(NCPU)
    !!!else
    !!!    call OMP_SET_NUM_THREADS(1)
    !!!endif

    !$OMP PARALLEL DO &
    !$OMP PRIVATE(x,cosx,sinx,wcosx,wsinx,C,S,YC,YS,CCh,CSh,SSh,CC,SS,CS,D)
    do I=1,NW

        x = W(I) * th
        cosx = cos(x)
        sinx = sin(x)
        wcosx = ww*cosx         ! attach weights
        wsinx = ww*sinx         ! attach weights

        C = sum(wcosx)         ! Eq. (8)
        S = sum(wsinx)         ! Eq. (9)

        YC = sum(yh*wcosx)     ! Eq. (17)
        YS = sum(yh*wsinx)     ! Eq. (18)
        CCh = sum(wcosx*cosx)  ! Eq. (13)
        CSh = sum(wcosx*sinx)  ! Eq. (15)
        SSh = 1.d0 - CCh
        CC = CCh - C*C         ! Eq. (13)
        SS = SSh - S*S         ! Eq. (14)
        CS = CSh - C*S         ! Eq. (15)
        D = CC*SS - CS*CS      ! Eq. (6)

        a(I) = (YC*SS-YS*CS) / D
        b(I) = (YS*CC-YC*CS) / D
        off(I) = -a(I)*C - b(I)*S

        upow(I) = (SS*YC*YC + CC*YS*YS - 2.d0*CS*YC*YS) / (YY*D) ! Eq. (5)

    end do
    !$OMP END PARALLEL DO

    ! An ad-hoc estimate of the number of independent frequencies 
    ! see discussion following Eq. (24)
    M = (maxval(W/pi2) - minval(W/pi2)) * (maxval(th) - minval(th))
    M = -1 ! not used
    P = upow
        
END SUBROUTINE GLOMBSCARGLE

