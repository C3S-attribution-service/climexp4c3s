        subroutine derivative(ndiff,data,diff,npermax,yrbeg,yrend
     +       ,nperyear,minfac,lwrite)
!
!       take a centered derivative of the data in data and return it in
!       diff
!       ndff: number of points in the derivative (must be odd)
!       
        implicit none
        integer nmax
        parameter (nmax=101)
        integer ndiff,npermax,yrbeg,yrend,nperyear
        real data(npermax,yrbeg:yrend),diff(npermax,yrbeg:yrend),minfac
        logical lwrite
        integer mo,yr,mm,ym,yp,mp,k,m,nn
        real xx(nmax),yy(nmax),sig(nmax),a,b,siga,sigb,chi2,q

        if ( ndiff.eq.3 ) then      ! old two-point algorithm
            do yr=yrbeg,yrend
                do mo=1,nperyear
                    mm = mo-1
                    call normon(mm,yr,ym,nperyear)
                    mp = mo + 1
                    call normon(mp,yr,yp,nperyear)
                    if ( ym.ge.yrbeg .and. yp.le.yrend ) then
                        if ( data(mm,ym).lt.1e33 .and. 
     +                       data(mp,yp).lt.1e33 ) then
                            diff(mo,yr) = (data(mp,yp) - data(mm,ym))/2
                        else
                            diff(mo,yr) = 3e33
                        end if
                    else
                        diff(mo,yr) = 3e33
                    end if
                end do
            end do
        else                    ! new regression algorithm
            do yr=yrbeg,yrend
                do mo=1,nperyear
                    nn = 0
                    do k=-ndiff/2,(ndiff-1)/2 ! rpund the center point up
                        m = mo + k
                        call normon(m,yr,yp,nperyear)
                        if ( yp.ge.yrbeg .and. yp.le.yrend ) then
                            if ( data(m,yp).lt.1e33 ) then
                                nn =nn + 1
                                if ( nn.gt.nmax ) then
                                    write(0,*) 'derivative: error: ',
     +                                   'can only handle ',nmax,
     +                                   '-point derivatives'
                                    call abort
                                end if
                                xx(nn) = k
                                yy(nn) = data(m,yp)
                            end if
                        end if
                        if ( nn.ge.minfac*ndiff ) then
                            call fit(xx,yy,nn,sig,0,a,diff(mo,yr),
     +                           siga,sigb,chi2,q)
                        else
                            diff(mo,yr) = 3e33
                        end if
                    end do
                end do
            end do
        end if
        end
