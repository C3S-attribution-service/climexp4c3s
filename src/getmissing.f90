subroutine getmissing(data,npermax,yrbeg,yrend,nens1,nens2,nperyear,j1,j2,y1,y2,nmissing,nn,n0missing,n0,nlength)
!
!   actually get the missing data statistics for the adta
!
    implicit none
    integer npermax,yrbeg,yrend,nens1,nens2,nperyear,j1,j2,y1,y2
    integer nmissing(yrbeg:yrend),nn(yrbeg:yrend),n0missing,n0,nlength(10)
    real data(npermax,yrbeg:yrend,0:nens2)
    integer yr1,yr2,yr,m,mo,mo1,mo2,mm1,mm2,nperday,length,iens
    logical doit
    integer,external :: leap

    n0 = 0
    n0missing = 0
    length = 0
    nlength = 0
    nmissing = 0
    nn = 0
    do yr1=y1,y2
        do mo1=1,nperyear
            do iens=nens1,nens2
                if ( data(mo1,yr1,iens) < 1e33 ) goto 10
            end do
        end do
    end do
    write(0,*) 'getmissing: error: no valid data'
    return
10  continue
    do yr2=y2,y1,-1
        do mo2=nperyear,1,-1
            do iens=nens1,nens2
                if ( data(mo2,yr2,iens) < 1e33 ) goto 20
            end do
        end do
    end do
    write(0,*) 'getmissing: error: no valid data'
    return
20  continue
    do yr=yr1,yr2
        if ( yr == yr1 ) then
            mm1 = mo1
        else
            mm1 = 1
        end if
        if ( yr == yr2 ) then
            mm2 = mo2
        else
            mm2 = nperyear
        end if
        do iens=nens1,nens2
            ! any valid data in this year?
            doit = .false.
            do m=max(j1,mm1),j2
                mo = m
                if ( mo > nperyear ) mo = mo - nperyear
                if ( mo > mm2 ) cycle
                if ( data(mo,yr,iens) < 1e33 ) doit = .true.
            end do
            ! if not, do not consider this station-year
            if ( .not.doit ) cycle
            do m=max(j1,mm1),j2
                mo = m
                if ( mo > nperyear ) mo = mo - nperyear
                if ( mo > mm2 ) cycle
                nn(yr) = nn(yr) + 1
                if ( data(mo,yr,iens) > 1e33 ) then
                    if ( mod(nperyear,366) == 0 ) then
                        nperday = nint(nperyear/365.)
                        if ( leap(yr) == 2 .or. mo <= 59*nperday .or. mo > 60*nperday ) then
                            nmissing(yr) = nmissing(yr) + 1
                            length = length + 1
                        end if
                    else
                        nmissing(yr) = nmissing(yr) + 1
                        length = length + 1
                    end if
                else
                    if ( length > 0 ) nlength(min(length,10)) = nlength(min(length,10)) + 1
                    length = 0
                end if
            end do
        end do
        n0 = n0 + nn(yr)
        n0missing = n0missing + nmissing(yr)
    end do
    y1 = yr1 
    y2 = yr2
end subroutine