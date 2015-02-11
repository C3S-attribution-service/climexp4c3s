*  #[ fitgau:
        subroutine fitgau(xx,ntot,mean,sd,a,b,minindx,maxindx,ntype,
     +       j1,j2,year,xyear,t,t25,t975,tx,tx25,tx975,
     +       confidenceinterval,lboot,lprint,lweb,lchangesign,lwrite)
*
*       a fit a gaussian distribution to the data
*
        implicit none
*
        integer ntot,ntype,j1,j2,year
        real xx(ntot),mean,sd,a,b,minindx,maxindx,xyear,tx,tx25,tx975
     +       ,t(10),t25(10),t975(10),confidenceinterval
        logical lboot,lprint,lweb,lchangesign,lwrite
*
        integer nmc
        parameter(nmc=1000)
        integer i,j,nx,iter,imc,ier
        real tol,p(3,2),q(2),xmin,y(3),aa(nmc),bb(nmc),tt(nmc,10)
     +       ,z,zz(10),x,f,txtx(nmc)
*
        integer nmax
        parameter(nmax=10000000)
        integer ncur
        real data(nmax),restrain
        logical llwrite
        common /fitdata1/ data
        common /fitdata2/ restrain,ncur,llwrite
*
        real llgauss,erfc
        real *8 ranf 
        external llgauss,ranf,erfc
!
        if ( lwrite ) then
            print *,'fitgau: input'
            print *,'mean,sd     = ',mean,sd
            print *,'min/maxindx = ',minindx,maxindx
            print *,'year,xyear  = ',year,xyear
        end if
*
*       first a trivial case which causes no end of trouble
*
        if ( sd.eq.0 ) then
            if ( lwrite ) print *,'fitgau: sd = 0, everything undfined'
            t = 3e33
            t25 = 3e33
            t975 = 3e33
            tx = 3e33
            tx25 = 3e33
            tx975 = 3e33
            return
        endif
*
        if ( minindx.lt.-1e33 .and. maxindx.gt.1e33 ) then
*
*           use mean and sd as parameters of the normal distribution
*
            a = mean
            b = sd
            if ( lprint ) print '(a)','# for Gaussian used mean, sd.'
            do i=1,10
                if ( mod(i,3).eq.1 ) then
                    f = 10.**((i+2)/3)
                elseif ( mod(i,3).eq.2 ) then
                    f = 2*10.**((i+1)/3)
                else
                    f = 5*10.**(i/3)
                endif
                f = f*(j2-j1+1)
                f = 1-2/f
                call merfi(f,zz(i),ier)
                t(i) = a + sqrt(2.)*b*zz(i)
            enddo
            if ( xyear.lt.1e33 ) then
                z = (xyear - a)/b
                if ( z.gt.12 ) then
                    tx = 3e33
                else
                    tx = 2/erfc(z/sqrt(2.))
                end if
                if ( lwrite ) print *,'tx = ',tx
            endif
            if ( .not.lboot ) then
                if ( lchangesign ) then
                    a = -a
                    t = -t
                endif
                return
            endif
            do imc=1,nmc
                do i=1,ntot
                    j = 1+int(ntot*ranf(i))
                    if ( j.lt.1 .or. j.gt.ntot ) then
                        write(0,*) 'fitgau: error: j = ',j
                        call abort
                    endif
                    data(i) = xx(j)
                enddo
                aa(imc) = 0
                bb(imc) = 0
                do i=1,ntot
                    aa(imc) = aa(imc) + data(i) - mean
                    bb(imc) = bb(imc) + (data(i) - mean)**2
                enddo
                aa(imc) = aa(imc)/ntot
                bb(imc) = bb(imc)/ntot - aa(imc)**2
                if ( bb(imc).lt.0 ) then
                    write(0,*) 'fitgau: error: var<0: ',bb(imc)
                    bb(imc) = 0
                endif
                bb(imc) = sqrt(bb(imc))
                aa(imc) = aa(imc) + mean
                do i=1,10
                    tt(imc,i) = aa(imc) + sqrt(2.)*bb(imc)*zz(i)
                enddo
                if ( xyear.lt.1e33 ) then
                    z = (xyear - aa(imc))/bb(imc)
                    txtx(imc) = 2/erfc(z/sqrt(2.))
                endif
            enddo               ! imc
            do i=1,10
                call getcut( t25(i),(100-confidenceinterval)/2,nmc,
     +               tt(1,i))
                call getcut(t975(i),(100+confidenceinterval)/2,nmc,
     +               tt(1,i))
            enddo
            if ( lchangesign ) then
                t25 = -t25
                t975 = -t975
            endif
            if ( xyear.lt.1e33 ) then
                call getcut( tx25,(100-confidenceinterval)/2,nmc,txtx)
                call getcut(tx975,(100+confidenceinterval)/2,nmc,txtx)
            endif
            if ( lprint ) then
                call printreturnvalue(ntype,t,t25,t975,lweb)
                call printreturntime(year,xyear,tx,tx25,tx975,lweb)
                call plotreturnvalue(ntype,t25,t975,j2-j1+1)
            endif
            if ( lchangesign ) then
                a = -a
                b = -b
            endif
            return
        else
*
*           copy to common for routine llgauss
*       
            ncur = ntot
            do i=1,ncur
                data(i) = xx(i)
            enddo
*       
*           fit, using Numerical Recipes routines
*
            p(1,1) = mean *0.9
            p(1,2) = sd   *0.9
            p(2,1) = p(1,1) *1.2
            p(2,2) = p(1,2)
            p(3,1) = p(1,1)
            p(3,2) = p(1,2) *1.2
            do i=1,3
                q(1) = p(i,1)
                q(2) = p(i,2)
                y(i) = llgauss(q)
            enddo
            tol = 1e-4
            call amoeba(p,y,3,2,2,tol,llgauss,iter)
*           maybe add restart later
            a = p(1,1)
            b = p(1,2)
*
*           output
*
            print '(a,i5,a)','# Fitted to Gaussian distribution in '
     +           ,iter,' iterations'
            print '(a)'
     +           ,'# p(x) = exp(-(x-a)^2/(2*b^2))/(b*sqrt(2*pi)) with'
            print '(a,f16.3)','# a = ',a
            print '(a,f16.3)','# b = ',b
        endif
        end
*  #] fitgau:
*  #[ llgauss:
        real function llgauss(p)
*
*       computes the log-likelihood function for a normal distribution
*       with parameters alpha,beta=p(1),p(2) and data in common.
*
        implicit none
*       
        real p(2)
*
        integer i
        real s
*
        integer nmax,ncur
        parameter(nmax=10000000)
        real data(nmax),restrain
        logical llwrite
        common /fitdata1/ data
        common /fitdata2/ restrain,ncur,llwrite
*       
        llgauss = 0
        do i=1,ncur
            llgauss = llgauss - ((data(i)-p(1))/p(2))**2/2
        enddo
        llgauss = llgauss - ncur*log(p(2))
*       normalization is not 1 in case of cut-offs
        call gausnorm(p(1),p(2),s)
        llgauss = llgauss - ncur*log(s)
*       minimum, not maximum
        llgauss = -llgauss
***        print *,'a,b,llgauss = ',p(1),p(2),llgauss
*
        end
*  #] llgauss:
*  #[ gausnorm:
        subroutine gausnorm(a,b,s)
        implicit none
#include "getopts.inc"
        real a,b,s
        real z1,z2,sqrt2
        real erfcc
        external erfcc
        sqrt2 = sqrt(2.)
        if ( minindx.gt.-1e33 ) then
            z1 = (minindx-a)/b
            if ( maxindx.lt.1e33 ) then
                z2 = (maxindx-a)/b
                s = (erfcc(z1/sqrt2) - erfcc(z2/sqrt2))/2
            else
                s = erfcc(z1/sqrt2)/2
            endif
        else
            if ( maxindx.lt.1e33 ) then
                z2 = (maxindx-a)/b
                s = erfcc(-z2/sqrt2)/2
            else
                s = 1
            endif
        endif
***        print *,'gausnorm: norm = ',a,b,s
        end
*  #] gausnorm:
