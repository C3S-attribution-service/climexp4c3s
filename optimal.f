        program testcor
*
* zoekt optimale combinatie predictors gegeven correlaties
*
        implicit none
        integer nmax
        parameter (nmax=50)
        integer i,j,n,indx(nmax)
        real xx(nmax,nmax),xy(nmax),xxi(nmax,nmax),d,a(nmax),var
        character string*80
*
    1   continue
        read (*,'(a)') string
        print '(a)',string
        print '(a)','give correlations <xi y>, end with EOF'
        do i=1,nmax
            read (*,*,end=100,err=100) xy(i)
        enddo
        print '(a,i5)','maximum is ',nmax
        stop
  100   continue
        n = i-1
        if ( n.eq.0 ) stop
        print '(a,i5,a)','Found ',n,' values'
        print '(10f8.2)',(xy(i),i=1,n)
        print '(a)','give correlations <xi xj> per row'
        do i=1,n
            read (*,*) (xx(j,i),j=1,n)
            print '(10f8.2)',(xx(j,i),j=1,n)
            if ( xx(i,i).ne.1 ) then
                print '(a)','error: expecting <xi xi> = 1, fixed'
                xx(i,i) = 1
            endif
        enddo
*
*       invert matrix xx (NumRec)
*
        do i=1,n
            do j=1,n
                xxi(j,i) = 0
            enddo
        enddo
        do i=1,n
            xxi(i,i) = 1
        enddo
        call ludcmp(xx,n,nmax,indx,d)
        do j=1,n
            call lubksb(xx,n,nmax,indx,xxi(1,j))
        enddo
*
*       compute loadings a
*
        do i=1,n
            a(i) = 0
            do j=1,n
                a(i) = a(i) + xy(j)*xxi(j,i)
            enddo
        enddo
        print '(a)','Coefficients xI are (minus factor sd(y)/sd(xj))'
        print '(10f8.2)',(a(i),i=1,n)
*
*       compute explained variance
*
        var = 0
        do i=1,n
            var = var + a(i)*xy(i)
        enddo
        print '(a)','Explained percentage variance is'
        print '(f8.2)',100*var
        goto 1
        end
*
*       Gerrits version
        subroutine aap
      x11 = 1.
      x22 = 1.
      x12 = 0.67
      s1  = 0.21
      s2  = 0.34
      s1  = 0.28
      s2  = 0.35
      x12 = 0.47
      s1  = 0.34
      s2  = 0.25
      s1  = 0.35
      s2  = 0.30
*
      x11i = 1/(1-x12**2)
      x22i = 1/(1-x12**2)
      x12i = -x12/(1-x12**2)
*
      rmax = s1**2*x11i+s2**2*x22i+2*s1*s2*x12i
*
      write(6,*) sqrt(rmax)
*
      x11 = 1.
      x22 = 1.
      x33 = 1.
      x12 = 0.67
      x23 = 0.47
      x13 = 0.50
      s1  = 0.21
      s2  = 0.34
      s3  = 0.25
      s1  = 0.28
      s2  = 0.35
      s3  = 0.30
      det = x11*x22*x33+2.*x12*x23*x13-x11*x23**2-x22*x13**2-x33*x12**2
      x11i= (x22*x33-x23**2)/det
      x22i= (x11*x33-x13**2)/det
      x33i= (x22*x33-x12**2)/det
      x12i= (x13*x23-x12*x33)/det
      x13i= (x12*x23-x13*x22)/det
      x23i= (x12*x13-x23*x11)/det
      write(6,*) 'ch11' , (x11*x11i+x12*x12i+x13*x13i)
      write(6,*) 'ch12' , (x11*x12i+x12*x22i+x13*x23i)
      write(6,*) 'ch13' , (x11*x13i+x12*x23i+x13*x33i)
      write(6,*) 'ch22' , (x12*x12i+x22*x22i+x23*x23i)
      write(6,*) 'ch23' , (x12*x13i+x22*x23i+x23*x33i)
      write(6,*) 'ch33' , (x13*x13i+x23*x23i+x33*x33i)
*
      write(6,*)
      rmax = s1**2*x11i+s2**2*x22i+2*s1*s2*x12i
     .      +s3**2*x33i+2*s1*s3*x13i+2*s2*s3*x23i
      write(6,*) sqrt(rmax)
      end

