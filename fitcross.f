        subroutine fitcross(x,y,ndata,sig,mwt,a,b,siga,sigb,chi2,q
     +       ,ncrossvalidate,aa,bb,laabb)
!
!       based on fit out of Numerical Recipes, but with cross-validation
!       (and many programming improvements, especially error checking)
!       26-oct-2009: added cross-validation
!       when ncrossvalidate>0, compute a,b with xbar with ncrossvalidate
!       points withheld
!       when laabb.eq..true. also compute bb(ndata) with b's with
!       ncrossvalidate points withheld.
!
        implicit none
        integer mwt,ndata,ncrossvalidate
        real sig(ndata),x(ndata),y(ndata)
        real a,b,chi2,q,siga,sigb,aa(ndata),bb(ndata)
        logical laabb
!
        integer i,j
        real sigdat,ss,st2,sx,sxoss,sy,t,wt
        real,allocatable :: sst2(:),ax(:),ssx(:),ssy(:),sss(:)
        logical lwrite
        real,external :: gammq

        lwrite = .false.
*       special case causes crashes later on...
        if ( ndata.eq.1 ) then
            a = y(1)
            b = 0
            siga = 3e33
            sigb = 3e33
            chi2 = 0
            q = 3e33
            if ( ncrossvalidate.gt.0 .and. laabb ) bb = 3e33
            return
        endif
        if ( ncrossvalidate.ge.ndata ) then
            a = 3e33
            b = 3e33
            siga = 3e33
            sigb = 3e33
            chi2 = 0
            q = 3e33
            if ( ncrossvalidate.gt.0 .and. laabb ) bb = 3e33
            return
        end if
        sx = 0
        sy = 0
        st2 = 0
        b = 0
        if ( ncrossvalidate.gt.0 ) then
            allocate(sst2(ndata))
            allocate(ax(ndata))
            allocate(ssx(ndata))
            allocate(ssy(ndata))
            allocate(sss(ndata))
            bb = 0
            sst2 = 0
        end if
        if ( mwt.ne.0 ) then
            ss = 0
            do i=1,ndata
                wt = 1/(sig(i)**2)
                ss = ss + wt
                sx = sx + x(i)*wt
                sy = sy + y(i)*wt
            end do
            do i=1,ndata
                if ( ncrossvalidate.eq.0 ) then
                    t = (x(i) - sx/ss)/sig(i)
                else
                    sss(i) = ss
                    ssx(i) = sx
                    ssy(i) = sy
                    do j=max(1,i-(ncrossvalidate-1)/2),
     +                   min(ndata,i+ncrossvalidate/2)
                        wt = 1/sig(j)**2
                        sss(i) = sss(i) - wt
                        ssx(i) = ssx(i) - x(j)*wt
                        ssy(i) = ssy(i) - y(j)*wt
                    end do
                    ax(i) = ssx(i)/sss(i)
                    t = (x(i) - ax(i))/sig(i)
                end if
                st2 = st2 + t**2
                b = b + t*y(i)/sig(i)
            end do
            if ( ncrossvalidate.gt.0 .and. laabb ) then
                do j=1,ndata
                    do i=1,ndata
                        if ( j.lt.i-(ncrossvalidate)/2 .or.
     +                       j.gt.i+(ncrossvalidate-1)/2 ) then
                            t = (x(i) - ax(j))/sig(i)
                            sst2(j) = sst2(j) + t**2
                            bb(j) = bb(j) + t*y(i)/sig(i)
                        end if
                    end do
                end do
            end if
        else                    ! mwt.eq.0, no weights
            do i=1,ndata
                sx = sx + x(i)
                sy = sy + y(i)
            end do
            ss = ndata
            do i=1,ndata
                if ( ncrossvalidate.eq.0 ) then
                    t = x(i) - sx/ss
                else
                    sss(i) = ss
                    ssx(i) = sx
                    ssy(i) = sy
                    do j=max(1,i-(ncrossvalidate-1)/2),
     +                   min(ndata,i+ncrossvalidate/2)
                        sss(i) = sss(i) - 1
                        ssx(i) = ssx(i) - x(j)
                        ssy(i) = ssy(i) - y(j)
                    end do
                    ax(i) = ssx(i)/sss(i)
                    t = x(i) - ax(i)
                end if
                st2 = st2 + t**2
                b = b + t*y(i)
            end do
            if ( ncrossvalidate.gt.0 .and. laabb ) then
                do j=1,ndata
                    do i=1,ndata
                        if ( j.lt.i-(ncrossvalidate)/2 .or.
     +                       j.gt.i+(ncrossvalidate-1)/2 ) then
                            t = x(i) - ax(j)
                            sst2(j) = sst2(j) + t**2
                            bb(j) = bb(j) + t*y(i)
                        end if
                    end do
                end do
            end if
        end if
        if ( st2.eq.0 ) then    ! all x-values equal, put b arbitrarily to zero
            b = 0
            sigb = 0
            a = sy/ss
            if ( ncrossvalidate.gt.0 .and. laabb ) then
                bb = 0
                do i=1,ndata
                    aa(i) = ssy(i)/sss(i)
                end do
            end if
            siga = 1/sqrt(ss)
        else
            b = b/st2
            a = (sy-sx*b)/ss
            if ( ncrossvalidate.gt.0 .and. laabb ) then
                do i=1,ndata
                    if ( sst2(i).gt.0 ) then
                        bb(i) = bb(i)/sst2(i)
                    else
                        bb(i) = 0
                    end if
                    aa(i) = (ssy(i) - bb(i)*ssx(i))/sss(i)
                end do
            end if
            siga = sqrt((1 + sx*sx/(ss*st2))/ss)
            if ( st2.lt.1e-10 ) then
                sigb = 3e33
            else
                sigb = sqrt(1/st2)
            end if
        end if
        chi2=0
        if ( ndata.eq.2 ) then
            siga = 3e33
            sigb = 3e33
            q = 1
        else
            if ( mwt.ne.0 ) then
                do i=1,ndata
                    chi2 = chi2 + ((y(i) - a - b*x(i))/sig(i))**2
                end do
                q = gammq(0.5*(ndata-2),0.5*chi2)
            else
                do i=1,ndata
                    chi2 = chi2 + (y(i) - a - b*x(i))**2
                end do
                q = 1
                sigdat = sqrt(chi2/(ndata-2))
                siga = siga*sigdat
                sigb = sigb*sigdat
            end if
        end if
        if ( lwrite ) print *,'sx,sy,ss,st2,a,b = ',sx,sy,ss,st2,a,b
        if ( ncrossvalidate.gt.0 .and. laabb ) then
            if ( lwrite ) then
                print *,'i,ssx,ssy,sss,sst2,aa,bb'
                do i=1,ndata
                    print *,i,ssx(i),ssy(i),sss(i),sst2(i),aa(i),bb(i)
                end do
            end if
        end if
        if ( ncrossvalidate.gt.0 ) then
            deallocate(sst2)
            deallocate(ax)
            deallocate(ssx)
            deallocate(ssy)
            deallocate(sss)
        end if
        end
