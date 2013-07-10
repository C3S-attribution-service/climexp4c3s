        subroutine bootstrap(y,x,f,ys,xs,n,rr,kind,ndecor,w1,w2
     +       ,ncrossvalidate)
*
*       compute error bounds on r using a bootstrap technique
*       Input: y(n),x(n)    data
*              f(n)         .true. if the first element of a new series
*              ys(n),xs(n)  scratch arrays
*              w1(n),w2(n)  scratch arrays
*              n            size of arrays
*              kind         1: normal correlation
*                           2: rank correlation
*              ndecor        decorrelation length
*       Output:rr(-2:2)     -2sigma,-1sigma,median,1sigma,2sigma of r
*
        implicit none
        integer n,kind,ndecor,ncrossvalidate
        real y(n),x(n),xs(n),ys(n),rr(-2:2),a,w1(n),w2(n)
        logical f(n)
*
        integer nboot
        parameter (nboot=800)
        integer i,j,k,l,m,ii
        real r(nboot),prob,z,ax,ay,sx,sy,sxy,df,d,probd,val(-2:2)
        logical lwrite
        save val,lwrite
*       
        real*8 ranf
        external ranf
*
        data val /-95.,-68.3,0,68.3,95./
        data lwrite /.false./
*
*       input
*
        if ( lwrite ) then
            print *,'bootstrap input: '
            do i=1,n
                print *,i,x(i),y(i),f(i)
            enddo
        endif
*
*       generate bootstrap ensembles and compute correlation
*
        m = n/ndecor
        m = m*ndecor
        df = m-2
        do i=1,nboot
            do j=1,m,ndecor
                ii = 0
 100            continue
                ii = ii + 1
                k = 1+min(n-ndecor,int((n-ndecor)*ranf(i+j)))
                do l=1,ndecor-1
                    if ( f(k+l) ) then
                        if ( ii.lt.100 ) then
***                            print *,'crossing boundary',k+l,ii
                            goto 100
                        endif
                        write(0,*) 'bootstrap: error: cannot find '//
     +                       'contiguous block of size ',ndecor
                        write(*,*) 'bootstrap: error: cannot find '//
     +                       'contiguous block of size ',ndecor
                        rr = 3e33
                        return
                    endif
                enddo
                do l=0,ndecor-1
                    xs(j+l) = x(k+l)
                    ys(j+l) = y(k+l)
                enddo
            enddo
            if ( kind.eq.1 ) then
                call pearsncross(ys,xs,m,r(i),prob,z,ax,sx,ay,sy,sxy,df
     +               ,ncrossvalidate)
            elseif ( kind.eq.2 ) then
                call spearx(ys,xs,m,w1,w2,d,z,probd,r(i),prob,m/(df+2)
     +                ,ax,sx,ay,sy)
            else
                write(0,*) 'bootstrap: error: kind = ',kind
     +                ,' not yet implemented'
                write(*,*) 'bootstrap: error: kind = ',kind
     +                ,' not yet implemented'
                call abort
            endif
        enddo
*
*       sort and return 95.4, 68.3% levels
*
        call nrsort(nboot,r)
        do i=-2,2
            a = 0.5 + nboot*(100+val(i))/200.
            j = int(a)
            a = a - j
            if ( j.lt.1 .or. j.ge.nboot ) then
                write(0,*) 'bootstrap: error: will not extrapolate'
                write(*,*) 'bootstrap: error: will not extrapolate'
                rr(i) = 3e33
            else
                rr(i) = (1-a)*r(j) + a*r(j+1)
            end if
        enddo
        end
