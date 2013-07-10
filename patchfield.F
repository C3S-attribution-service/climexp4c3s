        program patchfield
*
*       fill in missing values in a field or extend it in time using
*       data from a second field. A linear regression on the overlap is 
*       used to adjust the second field to the first one at each grid point.
*
        implicit none
#include "params.h"
        integer i,j,dy,mo,yr,nperyear,nperyear2,n,nperday,dpm(12)
        real xx(30*(yrend-yrbeg+1)),yy(30*(yrend-yrbeg+1))
        real scale(12),offset(12),siga,sigb,chi2,q,sig(1)
        real,allocatable :: mainfield(:,:,:,:,:),auxdata(:,:,:,:,:)
        character string*1024,var*40,units*20,var2*40,units2*20,
     +       mainfile*255,auxfile*255
        logical lwrite
        integer iargc
        data dpm /31,29,31,20,31,30,31,31,30,31,30,31/

        lwrite = .false.
        if ( iargc().lt.2 ) then
            write(0,*) 'usage: patchfield mainfield auxfield'
            write(0,*) 'patches holes in mainfield using data from '//
     +           'auxfield linearly regressed on mainfield at each '//
     +           'grid point'
            stop
        endif
        call getarg(1,mainfile)
        i = index(string,':')
        call readseries(mainfile,maindata,npermax,yrbeg,yrend,nperyear
     +       ,var,units,.false.,lwrite)
        call getarg(2,auxfile)
        call readseries(auxfile,auxdata,npermax,yrbeg,yrend,nperyear2
     +       ,var2,units2,.false.,lwrite)
        if ( nperyear.ne.nperyear2 ) then
            write(0,*) 'patchseries: error: can only handle series '//
     +           'with the same time resolution, not ',nperyear
     +           ,nperyear2
            write(*,*) 'patchseries: error: can only handle series '//
     +           'with the same time resolution, not ',nperyear
     +           ,nperyear2
            call abort
        end if
!
!       determine regression coefficents per month or less
!
        if ( nperyear.le.12 ) then
            do mo=1,nperyear
                n = 0
                do yr=yrbeg,yrend
                    if ( maindata(mo,yr).lt.1e33 .and.
     +                   auxdata(mo,yr).lt.1e33 ) then
                        n = n + 1
                        xx(n) = auxdata(mo,yr)
                        yy(n) = maindata(mo,yr)
                    end if
                end do
                if ( n.lt.7 ) then ! arbitrary
                    scale(mo) = 3e33
                    offset(mo) = 3e33
                else
                    call fit(xx,yy,n,sig,0,offset(mo),scale(mo),siga
     +                   ,sigb,chi2,q)
                end if
            end do
            do yr=yrbeg,yrend
                do mo=1,nperyear
                    if ( maindata(mo,yr).gt.1e33 .and.
     +                   auxdata(mo,yr).lt.1e33 ) then
                        maindata(mo,yr) =
     +                       offset(mo) + scale(mo)*auxdata(mo,yr)
                    end if
                end do
            end do
        else if ( nperyear.ge.260 ) then ! daily or 6-hourly frequency
            if ( nperyear.eq.360 .or. nperyear.eq.4*360 ) then
                dpm = 30
                nperday = nperyear/360
            else if ( nperyear.eq.365 .or. nperyear.eq.4*365 ) then
                dpm(2) = 28
                nperday = nperyear/365
            else
                nperday = nint(nperyear/366.)
            end if
            do mo=1,12
                n = 0
                do yr=yrbeg,yrend
                    do dy=1,nperday*dpm(mo)
                        call invgetdymo(dy,mo,j,nperyear)
                        if ( maindata(j,yr).lt.1e33 .and.
     +                       auxdata(j,yr).lt.1e33 ) then
                            n = n + 1
                            xx(n) = auxdata(j,yr)
                            yy(n) = maindata(j,yr)
                        end if
                    end do
                end do
                if ( n.lt.7 ) then ! arbitrary
                    scale(mo) = 3e33
                    offset(mo) = 3e33
                else
                    call fit(xx,yy,n,sig,0,offset(mo),scale(mo),siga
     +                   ,sigb,chi2,q)
                end if
            end do
            do yr=yrbeg,yrend
                do mo=1,12
                    do dy=1,nperday*dpm(mo)
                        call invgetdymo(dy,mo,j,nperyear)
                        if ( maindata(j,yr).gt.1e33 .and.
     +                       auxdata(j,yr).lt.1e33 ) then
                            maindata(j,yr) =
     +                           offset(mo) + scale(mo)*auxdata(j,yr)
                        end if
                    end do
                end do
            end do
        else
           write(0,*) 'merging pentad or decadal time series not ready'
           call abort 
        end if
!
!       output
!
        call copyheader(mainfile,6)
        write(6,'(a,g16.6)') '# patched with data from ',trim(auxfile)
        write(6,'(a,g16.6)') '# using offset,scale '
        do mo=1,min(12,nperyear)
            write(6,'(a,2f8.3)') '# ',offset(mo),scale(mo)
        end do
        call printdatfile(6,maindata,npermax,nperyear,yrbeg,yrend)
        end
