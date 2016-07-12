program count_missing
!
!   count how many values are missing in a time series
!
    implicit none
    include 'param.inc'
    integer yr,mo,nperyear,nmissing(yrbeg:yrend),nn(yrbeg:yrend)
    integer n0,n0missing,y1,y2,nlength(10),length
    real data(npermax,yrbeg:yrend)
    logical lstandardunits,lwrite
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
    
    print '(a)','  yr  ntot nmiss percentage'
    do yr=y1,y2
        print '(i4,2i6,f7.1)',yr,nn(yr),nmissing(yr),100*real(nmissing(yr))/real(nn(yr))
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