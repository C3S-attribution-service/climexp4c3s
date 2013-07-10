        subroutine fitnoisemodel(xx,yy,nn,a1,da1,a2,da2,bn,dbn,
     +       skeweps,dskeweps,sdeps,lboot,lwrite)
*
*       fit a state-dependent noise term
*
        implicit none
        integer nmax,nboot
        parameter (nmax=105000,nboot=100)
        integer nn
        real xx(nn),yy(nn),a1,da1,a2,da2,bn,dbn,skeweps,dskeweps,sdeps
        integer i,j,n,iter,iboot
        logical lboot,lwrite
        real x(nmax),y(nmax),p(5,4),q(4),f(5),tol,pp(5,4),qq(nboot,6),
     +       qqq(nboot),qcut(2,6),x16,x84,d,sig(1),da,db,chi2,q1,e2
        common /cfitfun/ x,y,n
        real getskew
        real*8 ranf
        external getskew,ranf
*
        if ( lwrite ) print *,'fitnoisemodel: input nn = ',nn
        if ( nn.lt.10 ) then
            a1 = 3e33
            da1 = 3e33
            a2 = 3e33
            da2 = 3e33
            bn = 3e33
            dbn = 3e33
            skeweps = 3e33
            dskeweps = 3e33
            sdeps = 3e33
            return
        endif
        if ( nn.gt.nmax ) then
            write(0,*) 'fitnoisemodel: increase nmax to ',nn
            write(*,*) 'fitnoisemodel: increase nmax to ',nn
            call abort
        endif
        n = nn
        do i=1,n
            x(i) = xx(i)
            y(i) = yy(i)
            if ( lwrite ) print *,i,xx(i),yy(i)
        enddo
*
*       first guess
*
        call fit(x,y,n,sig,0,p(1,1),p(1,2),da,db,chi2,q1)
        p(1,3) = 0
        p(1,4) = 0
        call fit_noise_point(x,y,n,p,iter,e2)
        a1 = p(1,2)
        if ( p(1,2).gt.0 ) then
            a2 = +p(1,3)
        else
            a2 = -p(1,3)
        endif
        bn = p(1,4)
        do j=1,4
            q(j) = p(1,j)
        enddo
	skeweps = getskew(q)
        sdeps = sqrt(e2)
        if ( lboot ) then
            do iboot=1,nboot
                do i=1,n
                    j = 1 + int(ranf(i)*n)
                    if ( j.lt.1 .or. j.gt.n ) then
                        write(0,*) 'error: j = ',j,n
                        call abort
                    endif
                    x(i) = xx(j)
                    y(i) = yy(j)
                enddo
                do j=1,4
                    pp(1,j) = p(1,j)
                enddo
                call fit_noise_point(xx,yy,n,pp,iter,e2)
                do j=1,4
                    qq(iboot,j) = pp(1,j)
                enddo
                do j=1,4
                    q(j) = p(1,j)
                enddo
                qq(iboot,5) = getskew(q)
            enddo
            do i=3,5
                do iboot=1,nboot
                    qqq(iboot) = qq(iboot,i)
                enddo
                call getcut(x16,16.0,nboot,qqq)
                call getcut(x84,84.0,nboot,qqq)
                d = (x84-x16)/2
                if ( i.eq.3 ) then
                    da2 = d
                elseif ( i.eq.4 ) then
                    dbn = d
                elseif ( i.eq.5 ) then
                    dskeweps = d
                else
                    write(0,*) 'error bnc571892'
                    call abort
                endif
            enddo
        else
            da1 = 3e33
            da2 = 3e33
            dbn = 3e33
            dskeweps = 3e33
        endif
        if ( lwrite ) then
            print *,'iter ',iter
            print *,'a1 = ',a1,da1
            print *,'a2 = ',a2,da2
            print *,'bn = ',bn,dbn
            print *,'skew ',skeweps,dskeweps
            print *,'s.d. ',sdeps
        endif
        end
