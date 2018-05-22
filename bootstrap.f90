subroutine bootstrap(y,x,f,ys,xs,n,rr,kind,ndecor,w1,w2,ncrossvalidate)
!
!   compute error bounds on r using a bootstrap technique
!   Input: y(n),x(n)    data
!          f(n)         .true. if the first element of a new series
!          ys(n),xs(n)  scratch arrays
!          w1(n),w2(n)  scratch arrays
!          n            size of arrays
!          kind         1: normal correlation
!                       2: rank correlation
!                       3: RMSE
!                       4: MAE
!          ndecor        decorrelation length
!   Output:rr(-2:2)     -2sigma,-1sigma,median,1sigma,2sigma of r
!
    implicit none
    integer :: n,kind,ndecor,ncrossvalidate
    real :: y(n),x(n),xs(n),ys(n),rr(-2:2),a,w1(n),w2(n)
    logical :: f(n)

    integer,parameter :: nboot=1000
    integer :: i,j,k,l,m,ii
    real :: r(nboot),prob,z,ax,ay,sx,sy,sxy,df,d,probd,val(-2:2),bias
    real,allocatable :: wy(:),wx(:)
    logical :: lwrite
    save val,lwrite

    real*8 :: ranf
    external ranf

    data val /-95.,-68.3,0,68.3,95./
    data lwrite / .false. /

!   input

    if ( lwrite ) then
        print *,'bootstrap input: '
        do i=1,n
            print *,i,x(i),y(i),f(i)
        enddo
    endif

!   generate bootstrap ensembles and compute correlation

    if ( kind == 3 .or. kind == 4 ) then
        allocate(wy(n))
        allocate(wx(n))
        wy = 1
        wx = 1
    end if
    m = n/ndecor
    m = m*ndecor
    df = m-2
    do i=1,nboot
        do j=1,m,ndecor
            ii = 0
            100 continue
            ii = ii + 1
            k = 1+min(n-ndecor,int((n-ndecor)*ranf(i+j)))
            do l=1,ndecor-1
                if ( f(k+l) ) then
                    if ( ii < 100 ) then
!**                     print *,'crossing boundary',k+l,ii
                        goto 100
                    endif
                    write(0,*) 'bootstrap: error: cannot find contiguous block of size ',ndecor
                    write(*,*) 'bootstrap: error: cannot find contiguous block of size ',ndecor
                    rr = 3e33
                    return
                endif
            enddo
            do l=0,ndecor-1
                xs(j+l) = x(k+l)
                ys(j+l) = y(k+l)
            enddo
        enddo
        if ( kind == 1 ) then
            call pearsncross(ys,xs,m,r(i),prob,z,ax,sx,ay,sy,sxy,df,ncrossvalidate)
        elseif ( kind == 2 ) then
            call spearx(ys,xs,m,w1,w2,d,z,probd,r(i),prob,m/(df+2),ax,sx,ay,sy)
        elseif ( kind == 3 ) then
            call getrms(ys,wy,xs,wx,m,bias,ax,sx,ay,sy,r(i))
        elseif ( kind == 4 ) then
            call getmae(ys,wy,xs,wx,m,bias,ax,sx,ay,sy,r(i))
        else
            write(0,*) 'bootstrap: error: kind = ',kind,' not yet implemented'
            write(*,*) 'bootstrap: error: kind = ',kind,' not yet implemented'
            call exit(-1)
        endif
    enddo
    if ( kind == 3 .or. kind == 4 ) then
        deallocate(wy)
        deallocate(wx)
    end if

!   sort and return 95.4, 68.3% levels

    call nrsort(nboot,r)
    do i=-2,2
        call getcut1(rr(i),(100+val(i))/2,nboot,r,lwrite)
    enddo
end subroutine bootstrap
