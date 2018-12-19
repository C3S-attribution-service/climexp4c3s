subroutine getweightedquant(point,nnens,nens,nmod,quant,lwrite)

!   Compute some qunatiles of a set of values:
!   min, 5%, 25%, 50%, 75%, 95%, max
!   All models are weighted equally

    implicit none
    integer,intent(in) :: nens,nmod,nnens(nens)
    real,intent(in) :: point(nens)
    real,intent(out) :: quant(-6:6)
    logical,intent(in) :: lwrite
    integer :: iens,imod,i,n,ii(1024),icut
    real :: aa(1024),ww(1024),cuts(-6:6),w,w1,x,x1,x2,sw,wn
    logical :: lskipundefined
    data cuts /0,0.025,0.05,0.10,0.17,0.25,0.50,0.75,0.83,0.90,0.95,0.975,1./

    lskipundefined = .true. 
    if ( lwrite ) print *,'getweightedquant: point,nnens = ',nens
    n = 0
    sw = 0
    do iens=1,nens
        if ( point(iens) < 1e33 .or. .not. lskipundefined ) then
            n = n + 1
            if ( n > 1024 ) then
                write(0,*) 'getweightedquant: error: array too small ',1024
                call exit(-1)
            end if
            aa(n) = point(iens)
            if ( nnens(iens) /= 0 ) then
                ww(n) = 1./nnens(iens)
                sw = sw + ww(n)
            else
                write(0,*) 'getweightedquant: error: found 0 ensemble members of iens = ',iens
                ww(n) = 1
            end if
        end if
        if ( lwrite ) print *,iens,point(iens),nnens(iens),sw
    end do
    if ( abs(sw-nmod) > 0.01 .and. .not. lskipundefined ) then
        write(0,*) 'getweightedquant: error: \sum w_n = ',sw,'but nmod = ',nmod
        call exit(-1)
    end if
    if ( n < 2 ) then
        if ( lwrite ) print *,'gethistquant: not enough valid points ',n
        quant = 3e33
        return
    end if
    call ffsort(aa,ii,n)
    if ( lwrite ) then
        print *,'gethistquant: sorted array:'
        do i=1,n
            print '(i3,2f10.3)',i,aa(ii(i)),ww(ii(i))
        end do
    endif

!   min & max are easy...
    quant(-6) = aa(ii(1))
    quant(+6) = aa(ii(n))

!   quantiles slightly less so
    wn = ww(ii(n))
    do icut=-5,5
        w = 0
        do i=1,n
            w1 = w
            w = w + ww(ii(i))
            if ( w > (nmod+wn)*cuts(icut) ) then
                exit
            end if
        end do
        if ( lwrite ) then
            print '(a,f7.3)','found quantile ',cuts(icut)
            if ( i > 1 ) print '(i4,3f7.3)',i-1,w1/(nmod+wn),aa(ii(i-1)),ww(ii(i-1))
            if ( i <= n ) print '(i4,3f7.3)',i,w/(nmod+wn),aa(ii(i)),ww(ii(i))
        end if
        if ( i == 1 ) then  ! extrapolate - not yet ready
            if ( lwrite ) print '(a)','taking min'
            quant(icut) = aa(ii(1))
        else if ( i <= n ) then ! interpolate
            if ( abs(aa(ii(i-1))) > 1e33 .or. &
            abs(aa(ii(i))) > 1e33 ) then
                ! in case lskipundef is .false.
                quant(icut) = 3e33 ! refuse to interpolate to undefs
            else
                x1 = cuts(icut) - w1/(nmod+wn) ! distance to previous point
                x2 = w/(nmod+wn) - cuts(icut) ! distance to next point
                x = x2/(x1+x2) ! unweighted
                quant(icut) = x*aa(ii(i-1)) + (1-x)*aa(ii(i)) ! interpolate
                if ( lwrite ) then
                    print '(a,2f7.3)','interpolating with weights ',x,1-x
                    print '(a,f7.3)','gives ',quant(icut)
                end if
            end if
        else                ! extrapolate - not yet ready
            if ( lwrite ) print '(a)','taking max'
            quant(icut) = aa(ii(n))
        end if
    end do
end subroutine getweightedquant
