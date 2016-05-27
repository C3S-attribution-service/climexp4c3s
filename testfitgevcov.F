        program testgevcov
        implicit none
#include "param.inc"
#include "getopts.inc"
        integer i,n,j1,j2,year,ntot,ntype
        integer nperyear,mens1,mens
        integer,allocatable :: yrs(:)
        real a,b,xi,alpha,x,y,z,xyear,mean,sd,dalpha,cov1,cov2,offset
        real inrestrain,t(10,3),t25(10,3),t975(10,3),
     +       tx(3),tx25(3),tx975(3)
        real,allocatable :: xx(:,:),yy(:)
        logical lboot,lprint
        real,external :: gasdev

        call getopts(1,0,nperyear,yrbeg,yrend,.false.,mens1
     +       ,mens)

        lweb = .true.
        lchangesign = .false.
        lboot = .true.
        lprint = .true.
        lwrite = .false.
        ntype = 2 ! gumbel plot

        n = 100
        allocate(xx(2,n))
        allocate(yy(n))
        allocate(yrs(0:n))
        print *,'a,b,alpha: '
        read *,a,b,alpha
        ! covariate
        do i=1,n
            xx(2,i) = i
            yrs(i) = i
        enddo
        ! time series
        ntot = n
        do i=1,n
            x = gasdev(i)
            xx(1,i) = (a + b*x + alpha*xx(2,i))
        end do
        mean = a
        sd = b
        dalpha = abs(alpha)
        j1 = 1
        j2 = 1
        cov1 = 1
        cov2 = n
        if ( .false. ) then
            year = yrend + 1
            xyear = 3e33
        else
            year = n
            xyear = xx(1,n)
            xx(1,n) = 3e33
            ntot = ntot - 1 ! do not send this point to the fit routine
        end if
        inrestrain = 0.2

        if ( alpha.eq.0 ) then
            do i=1,n
                yy(i) = xx(1,i)
            end do
            call fitgev(yy,ntot,mean,sd,a,b,xi,j1,j2
     +           ,lweb,ntype,lchangesign,year,xyear,t,t25,t975
     +           ,tx,tx25,tx975,inrestrain,lboot,lprint,lwrite)
        else
            call fitgevcov(xx,yrs,ntot,a,b,xi,alpha,j1,j2
     +           ,lweb,ntype,lchangesign,year,xyear,cov1,cov2,offset,t,t25,t975
     +           ,tx,tx25,tx975,inrestrain,lboot,lprint,lwrite)
        end if
        end
