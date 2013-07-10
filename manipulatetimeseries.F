        subroutine manipulatetimeseries(fcst,obs,npermax,yrbeg,yrend
     +       ,nperyear,nmodel1,var,jx,jy,multimodel,nmodel)
        implicit none
#include "getopts.inc"
        integer npermax,yrbeg,yrend,nperyear,nmodel1,jx,jy,nmodel
        integer multimodel(0:nmodel)
        real fcst(npermax,yrbeg:yrend,0:nens2),
     +       obs(npermax,yrbeg:yrend)
        character*(*) var
        integer yrmin,yrmax,i,j,n,yr,mo,iens,imodel,iens1,iens2
        integer,allocatable :: validens(:),ndup(:)
        real s1,s2,xmin
        logical lvalid
*
        allocate(validens(nens2+1))
        allocate(ndup(nperyear))
*
*       get mean and begin- and end years
*
        yrmin = yrend
        yrmax = yrbeg
        if ( nmodel1.eq.-1 ) then
            s1 = 0
            n = 0
            do yr=yr1,yr2
                do mo=1,nperyear
                    if ( obs(mo,yr).lt.1e33 ) then
                        if ( obs(mo,yr).gt.1e12 ) then
                            write(0,*) 'manupulatetimeseries: '
     +                           ,'suspect value of obs(',mo,yr,') = '
     +                           ,obs(mo,yr)
                        endif
                        yrmin = min(yrmin,yr)
                        yrmax = max(yrmax,yr)
                        s1 = s1 + obs(mo,yr)
                        n = n + 1
                    endif
                enddo
            enddo
            if ( n.eq.0 ) then
                fcst = 3e33
                goto 999
            endif
            s1 = s1/n
            yr1 = max(yr1,yrmin)
            yr2 = min(yr2,yrmax)
        else
            yr1 = max(yr1,yrbeg)
            yr2 = min(yr2,yrend)
        endif
        s2 = 0
        n = 0
        yrmin = yrend
        yrmax = yrbeg
*       we'd like to keep all observations, even if there are no
*       forecasts, to better define the climatology
        do iens=nens1,nens2
            do yr=yr1,yr2
                do mo=1,nperyear
                    if ( fcst(mo,yr,iens).lt.1e33 ) then
                        yrmin = min(yrmin,yr)
                        yrmax = max(yrmax,yr)
                        s2 = s2 + fcst(mo,yr,iens)
                        n = n + 1
                    endif
                enddo
            enddo
        enddo
        if ( n.eq.0 ) then
            if ( lwrite ) then
                write(*,*) 'manipulatetimeseries: no valid forecasts at'
     +               ,jx,jy
            end if
            goto 999
        endif
        s2 = s2/n
#if 0
        yr1 = max(yr1,yrmin)
        yr2 = min(yr2,yrmax)
#endif
*
*       sum
*
        if ( nmodel1.eq.-1 .and. lsum.gt.1 ) then
            if ( lwrite ) print '(a,i3)','# Summing series ',lsum
            call sumit(obs(1,yrbeg),npermax,nperyear,yrbeg,yrend,lsum
     +           ,oper)
        endif
        if ( lsum2.gt.1 ) then
            do iens=nens1,nens2
                call sumit(fcst(1,yrbeg,iens),npermax,nperyear
     +                        ,yrbeg,yrend,lsum2,'v')
            enddo
        endif
*
*       logscale
*
        if ( logscale ) then
            if ( lwrite ) print '(a,2i3)','# Taking log of series '
            call takelog(obs(1,yrbeg),npermax,nperyear,yrbeg,yrend)
            do iens=nens1,nens2
                call takelog(fcst(1,yrbeg,iens),npermax,nperyear,yrbeg
     +               ,yrend)
            enddo
        endif
*
*       sqrtscale
*
        if ( sqrtscale ) then
            if ( lwrite ) print '(a,2i3)','# Taking sqrt of series '
            call takesqrt(obs(1,yrbeg),nperyear,nperyear,yrbeg,yrend)
            do iens=nens1,nens2
                call takesqrt(fcst(1,yrbeg,iens),npermax,nperyear,yrbeg
     +               ,yrend)
            enddo
        endif
!
!       simple high-pas and low-pass filters
!
        if ( ndiff.ne.0 ) then
            if ( nmodel1.eq.-1 ) then
                if ( lwrite ) print *,'Taking differences of obs'
                call diffit(obs,npermax,nperyear,yrbeg,yrend,ndiff)
            end if
            do iens=nens1,nens2
                if ( lwrite ) print *,'Taking difference of ens ',iens
                minfacsum = 0.99
                call ndiffit(fcst(1,yrbeg,iens),npermax,nperyear
     +               ,yrbeg,yrend,ndiff,minfacsum)
            end do
            if ( lnooverlap .and. ndiff.lt.0 ) then
                call dooverlap(fcst,npermax,yrbeg,yrend,nens1,nens2
     +               ,nperyear,ndiff)
            end if
        end if
*
*       detrending
*
        if ( ldetrend ) then
            if ( nmodel1.eq.-1 .) then
                if ( lwrite ) print *,'Detrending obs'
                call detrend(obs,npermax,nperyear,yrbeg,yrend,yr1,yr2,m1
     +               ,m2,lsel)
            endif
            do iens=nens1,nens2
                if ( lwrite ) print *,'Detrending ens ',iens
                call detrend(fcst(1,yrbeg,iens),npermax,nperyear
     +               ,yrbeg,yrend,yr1,yr2,m1,m2,lsel)
            enddo
        endif
*
*       copy ensemble members so that there is the same
*       number of valid ones at every time step
*
        ndup = 0
        if ( nens2.gt.nens1 .and. lmakeensfull ) then
            call makeensfull(ndup,nperyear,fcst,npermax,yrbeg,yrend
     +           ,nens1,nens2,validens,lwrite)
        endif
*
*       bias corrections
*
        if ( nmodel1.eq.-1 ) then
            do imodel=1,nmodel
                iens1 = multimodel(imodel-1)
                iens2 = multimodel(imodel) - 1
                if ( debias.eq.1 ) then
                    call debiasmean(obs,fcst,npermax,nperyear,yrbeg
     +                   ,yrend,yr1,yr2,iens1,iens2,var,lwrite)
                elseif ( debias.eq.2 ) then
                    call debiasvar(obs,fcst,npermax,nperyear,yrbeg,yrend
     +                   ,yr1,yr2,iens1,iens2,var,lwrite)
                elseif ( debias.eq.3 ) then
                    call debiasall(obs,fcst,npermax,nperyear,yrbeg,yrend
     +                   ,yr1,yr2,iens1,iens2,var,lwrite)
                endif
            enddo
            if ( debias.gt.0 ) then
                ncrossvalidate = 0 ! is already included
            end if
        endif
*
*       anomalies
*
        if ( anom ) then
            if ( lwrite ) print '(a)','# Taking anomalies '
            if ( nmodel1.eq.-1 ) then
                call anomal(obs,npermax,nperyear,yrbeg,yrend,yr1,yr2)
            endif
            do iens=nens1,nens2
                call anomal(fcst(1,yrbeg,iens),npermax,nperyear
     +               ,yrbeg,yrend,yr1,yr2)
            enddo
        endif
*
 999    continue
        deallocate(validens)
        deallocate(ndup)
*
        end
