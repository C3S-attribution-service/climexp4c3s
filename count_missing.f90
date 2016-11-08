program count_missing
!
!   count how many values are missing in a time series
!
    implicit none
    include 'param.inc'
    integer yr,mo,nperyear,nmissing(yrbeg:yrend),nn(yrbeg:yrend)
    integer n0,n0missing,y1,y2,nlength(10),length,nperday
    real data(npermax,yrbeg:yrend),frac(yrbeg:yrend)
    logical lstandardunits,lwrite,nomissing
    character file*1023,var*80,units*40,line*80
    integer iargc
    integer,external :: leap
    lwrite = .false.
    
    if ( iargc() < 1 .or. iargc() > 3 ) then
        print *,'count_missing file [yr1 yr2]'
        call exit(-1)
    end if
    
    call getarg(1,file)
    lstandardunits = .false.
    call readseries(file,data,npermax,yrbeg,yrend,nperyear,var,units,lstandardunits,lwrite)
    if ( iargc() > 1 ) then
        call getarg(2,line)
        read(line,*) y1
    else
        y1 = yrbeg
    end if
    if ( iargc() > 2 ) then
        call getarg(3,line)
        read(line,*) y2
    else
        y2 = yrend
    end if

    call getmissing(data,npermax,yrbeg,yrend,0,0,nperyear,y1,y2,nmissing,nn,n0missing,n0,nlength)
    
    nomissing = .true.
    do yr=y1,y2
        if ( mod(nperyear,366) /= 0 ) then
            frac(yr) = real(nmissing(yr))/real(nn(yr))
        else
            nperday = nint(nperyear/366.)
            if ( leap(yr) == 1 ) then
                frac(yr) = real(nmissing(yr))/(nperday*365.)
            else
                frac(yr) = real(nmissing(yr))/(nperday*366.)
            end if
        end if
        if ( nmissing(yr) > 0 ) nomissing = .false.
    end do
    if ( nomissing ) then
        print '(a)','# no missing data'
    else
        call copyheader(file,6)
        print '(a)','# miss [1] fraction of missing data'
        do yr=y1,y2
            print '(i4,f8.4)',yr,frac(yr)
        end do
        print '(a4,2i6,f7.1)','# all ',n0,n0missing,100*real(n0missing)/real(n0)
        print '(a)','# length of gaps'
        do length = 1,10
            print '(a,i2,i6)','#',length,nlength(length)
        end do
        print '(a)','# (the last line includes all longer gaps)'     
    end if
end program