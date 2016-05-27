program count_missing
!
!   count how many values are missing in a time series
!
    implicit none
    include 'param.inc'
    integer yr,mo,nperyear,n,nmissing,yr1,yr2,mo1,mo2
    integer n0,n0missing,m1,m2,y1,y2,length,nlength(10)
    real data(npermax,yrbeg:yrend)
    logical lstandardunits,lwrite
    character file*1023,var*80,units*40,line*80
    integer iargc
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

    do yr1=y1,y2
        do mo1=1,nperyear
            if ( data(mo1,yr1) < 1e33 ) goto 10
        end do
    end do
    write(0,*) 'count_missing: error: no valid data'
    call exit(-1)
10  continue
    do yr2=y2,y1,-1
        do mo2=nperyear,1,-1
            if ( data(mo2,yr2) < 1e33 ) goto 20
        end do
    end do
    write(0,*) 'count_missing: error: no valid data'
    call exit(-1)
20  continue
    n0 = 0
    n0missing = 0
    length = 0
    nlength = 0
    print '(a)','  yr  ntot nmiss percentage'
    do yr=yr1,yr2
        n = 0
        nmissing = 0
        if ( yr == yr1 ) then
            m1 = mo1
        else
            m1 = 1
        end if
        if ( yr == yr2 ) then
            m2 = mo2
        else
            m2 = nperyear
        end if
        do mo=m1,m2
            n = n + 1
            if ( data(mo,yr) > 1e33 ) then
                nmissing = nmissing + 1
                length = length + 1
            else
                if ( length > 0 ) nlength(min(length,10)) = nlength(min(length,10)) + 1
                length = 0
            end if
        end do
        print '(i4,2i6,f7.1)',yr,n,nmissing,100*real(nmissing)/real(n)
        n0 = n0 + n
        n0missing = n0missing + nmissing
    end do
    print '(a)','============================'
    print '(a4,2i6,f7.1)','all ',n0,n0missing,100*real(n0missing)/real(n0)
    print '(a)'
    print '(a)','length of gaps'
    do length = 1,10
        print '(i2,i6)',length,nlength(length)
    end do
    print '(a)','(the last line includes all longer gaps)'     
    
end program