        program testsavgol
*
*       print out the Savitzky-Golay filtering coefficients
*
        implicit none
        integer i,j,n,m
        real c(1000)
    1   print *,'Give # of left=right points, degree'
        read *,n,m
        call savgol(c,1000,n,n,1,m)
        print *,'Coefficients for first derivative are'
        do i=n,1,-1
            print '(i3,f8.4)',-i,c(1000-i+1)
        enddo
        do i=1,n+1
            print '(i3,f8.4)',i-1,c(i)
        enddo
        goto 1
        end
