subroutine derivative(ndiff,data,diff,npermax,yrbeg,yrend,nperyear,minfac,lwrite)

!   take a centered derivative of the data in data and return it in diff
!   ndff: number of points in the derivative

    implicit none
    integer :: nmax
    parameter (nmax=101)
    integer :: ndiff,npermax,yrbeg,yrend,nperyear
    real :: data(npermax,yrbeg:yrend),diff(npermax,yrbeg:yrend),minfac
    logical :: lwrite
    integer :: mo,yr,mm,ym,yp,mp,k,m,nn
    real :: xx(nmax),yy(nmax),sig(nmax),a,b,siga,sigb,chi2,q

    if ( ndiff == 2 ) then
        call mdiffit(data,npermax,nperyear,yrbeg,yrend,1,minfac)
        diff = data
    else if ( ndiff == 3 ) then      ! old two-point algorithm
        do yr=yrbeg,yrend
            do mo=1,nperyear
                mm = mo-1
                call normon(mm,yr,ym,nperyear)
                mp = mo + 1
                call normon(mp,yr,yp,nperyear)
                if ( ym >= yrbeg .and. yp <= yrend ) then
                    if ( data(mm,ym) < 1e33 .and. &
                    data(mp,yp) < 1e33 ) then
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
                do k=-ndiff/2,(ndiff-1)/2 ! round the center point up
                    m = mo + k
                    call normon(m,yr,yp,nperyear)
                    if ( yp >= yrbeg .and. yp <= yrend ) then
                        if ( data(m,yp) < 1e33 ) then
                            nn =nn + 1
                            if ( nn > nmax ) then
                                write(0,*) 'derivative: error: ', &
                                'can only handle ',nmax, &
                                '-point derivatives'
                                call exit(-1)
                            end if
                            xx(nn) = k
                            yy(nn) = data(m,yp)
                        end if
                    end if
                    if ( nn >= minfac*ndiff ) then
                        call fit(xx,yy,nn,sig,0,a,diff(mo,yr), &
                        siga,sigb,chi2,q)
                    else
                        diff(mo,yr) = 3e33
                    end if
                end do
            end do
        end do
    end if
end subroutine derivative
