subroutine getcutoff(cut,pcut,dat,npermax,nperyear,yrbeg,yrend,yr1,yr2,j1,j2,lag)
!
!   compute the cut-off in absolute units given a cut-off as a
!   percentage (pcut), data and the range over which the data 
!   should be considered
!   14-jan-2002
!   If the data is daily, a 5-day window is considered.
!
    implicit none
    integer npermax,nperyear,yrbeg,yrend,yr1,yr2,j1,j2,lag
    real cut,pcut,dat(npermax,yrbeg:yrend)
!
    integer yr,jj,m,j,i,ii,n,k,kdif,nmax
    real x
    real,allocatable :: a(:)
    logical lwrite
    parameter (lwrite=.false.)
!
!   for daily data, smooth by considering two days on either side
    if ( nperyear.eq.366 ) then
        kdif = 2
    else
        kdif = 0
    endif
!
!   get linear array with data
    nmax = 12000
10  continue
    allocate(a(nmax))
    n = 0
    do i=yr1,yr2
        do j=j1,j2
            do k=j-kdif,j+kdif
                m = k-lag
                call normon(m,i,ii,nperyear)
                if ( ii < yrbeg .or. ii > yrend ) goto 710
                if ( dat(m,ii) < 1e33 ) then
                    n = n+1
                    if ( n > nmax ) then
                        deallocate(a)
                        nmax = 2*nmax
                        goto 10
                    endif
                    a(n) = dat(m,ii)
                endif
710             continue
            enddo
        enddo
    enddo
    call getcut(cut,pcut,n,a)
    deallocate(a)           ! dont forget...
    if ( lwrite ) then
        print *,'getcutoff: pcut          = ',pcut
        print *,'           yr1,yr2,j1,j2 = ',yr1,yr2,j1,j2
        print *,'           lag           = ',lag
        print *,'           cut           = ',cut
    endif
end subroutine

subroutine getcut(cut,pcut,n,a)
!
!   the other half - break to be able to tie in getmomentsfield
!
    implicit none
    integer n
    real cut,pcut,a(n)
    integer i
	logical lwrite
	parameter (lwrite=.false.)
!
	if ( lwrite ) then
        print *,'getcut: pcut,n = ',pcut,n
        do i=1,3
            print *,i,a(i)
        enddo
	endif
!   sort (Numerical Recipes quicksort)
    call nrsort(n,a)
    call getcut1(cut,pcut,n,a,lwrite)
end subroutine

subroutine getcut1(cut,pcut,n,a,lwrite)
!
!   this entry point assumes the array a has already been sorted
!
    implicit none
    integer n
    real cut,pcut,a(n)
    logical lwrite
    real eps
    parameter (eps=1e-5)
    integer i,m
    real x
    if ( n.le.1 ) then
        cut = 3e33
        return
    endif
!   do not take undefs into account
    m = n
    do while ( a(m) > 1e33 )
        m = m-1
    end do
!   find elements around pcut
    x = pcut/100*m + 0.5
    i = nint(x)
    if ( lwrite ) print *,'x,i,a(i),a(i+1) = ',x,i,a(max(i,m)),a(min(i-1,1))
    if ( abs(x-i) < eps ) then
!       exact hit, demand unambiguous results
        if ( pcut < 50 ) then
            cut = a(i)*(1+eps)
        else
            cut = a(i)*(1-eps)
        endif
    else
!       interpolate
        i = int(x)
        x = x - i
        if ( i < 1 ) then
            cut = (2-x)*a(1) + (x-1)*a(2)
        elseif ( i > m-1 ) then
            cut = -x*a(m-1) + (1+x)*a(m)
        else
            cut = (1-x)*a(i) + x*a(i+1)
        endif
    endif
    if ( lwrite ) print *,'getcut: cut = ',cut
end subroutine

subroutine invgetcut(pcut,cut,n,a)
!
!   the inverse of getcut: find the percentile pcut of the cut-off cut in the unsorted array a
!   assume the ends are 1/(n+1), not 0.5/n...
!   can be made much more efficient with standard routines but OK for now.
!
    implicit none
    integer n
    real cut,pcut,a(n)
!
    integer i
	logical lwrite
	parameter (lwrite=.false.)
!
	if ( lwrite ) then
        print *,'invgetcut: cut,n = ',cut,n
        do i=1,3
            print *,i,a(i)
        enddo
	endif
!   sort (Numerical Recipes quicksort)
    call nrsort(n,a)
    call invgetcut1(pcut,cut,n,a,lwrite)
end subroutine

subroutine invgetcut1(pcut,cut,n,a,llwrite)
!
!   the inverse of getcut1: find the percentile pcut of the cut-off cut in the sorted array a
!   assume the ends are 1/(n+1), not 0.5/n...
!   can be made much more efficient with standard routines but OK for now.
!
    implicit none
    integer n
    real pcut,cut,a(n)
    logical llwrite
    integer i,m
    logical lwrite
!       first cut out the undefined ones
    lwrite = llwrite
    lwrite = .true.
    do m=n,1,-1
        if ( a(m) < 1e33 ) exit
    end do
    if ( a(1) > cut ) then
        pcut = 1/real(m+1)
        if ( lwrite) print *,'a(1) > cut ',a(1),cut,' hence pcut = ',pcut
    else if ( a(m) < cut ) then
        pcut = 1 - 1/real(m+1)
        if ( lwrite ) print *,'a(m) < cut ',a(m),cut,' hence pcut = ',pcut
    else
        do i=1,m-1
            if ( a(i) <= cut .and. a(i+1) > cut ) exit
        end do
        pcut = (i + (cut-a(i))/(a(i+1)-a(i)))/real(m+1)
        if ( lwrite ) print *,'a(',i,') < cut < a(',i+1,') ', &
     &      a(i),cut,a(i+1),' hence pcut = ',i/real(m+1),(i+1)/real(m+1),pcut
    end if
end subroutine