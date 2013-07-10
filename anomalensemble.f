        subroutine anomalensemble(data,mpermax,nperyear,yrbeg,yrend,
     +           yr1,yr2,nens1,nens2)
!
!       change data to anomalies relative to the ensemble mean
!
        implicit none
        integer mpermax,nperyear,yrbeg,yrend,yr1,yr2,nens1,nens2
        real data(mpermax,yrbeg:yrend,0:nens2)
        integer mo,yr,iens,n
        real,allocatable :: mean(:,:)
!
!       sanity check
!
        if ( nens2.eq.nens1 ) then
            write(0,*)
     +           'anomalensemble: cannot take anomalies, no ensemble'
            write(0,*) nens1,nens2
            return
        endif
!
!       compute ensemble mean
!
        allocate(mean(nperyear,yr1:yr2))
        do yr=yr1,yr2
            do mo=1,nperyear
                mean(mo,yr) = 0
                n = 0
                do iens=nens1,nens2
                    if ( data(mo,yr,iens).lt.1e33 ) then
                        n = n + 1
                        mean(mo,yr) = mean(mo,yr) + data(mo,yr,iens)
                    endif
                enddo
                if ( n.gt.0 ) then
                    mean(mo,yr) = mean(mo,yr)/n
                else
                    mean(mo,yr) = 3e33
                endif
            enddo
        enddo
!
!       subtract from data
!
        do iens=nens1,nens2
            do yr=yr1,yr2
                do mo=1,nperyear
                    if ( data(mo,yr,iens).lt.1e33 .and.
     +                   mean(mo,yr).lt.1e33 ) then
                        data(mo,yr,iens) = data(mo,yr,iens)-mean(mo,yr)
                    else
                        data(mo,yr,iens) = 3e33
                    endif
                enddo
            enddo
        enddo
!
!       set rest of data to undefined
!
        do iens=0,nens1-1
            data(:,:,iens) = 3e33
        enddo
        do yr=yrbeg,yr1-1
            data(:,yr,:) = 3e33
        enddo
        do yr=yr2+1,yrend
            data(:,yr,:) = 3e33
        enddo
!
!       clean up
!
        deallocate(mean)
        end
