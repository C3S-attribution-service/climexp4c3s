program count_missing
!
!   count how many values are missing in a time series
!
    implicit none
    include 'param.inc'
    include 'getopts.inc'
    integer   :: yr,mo,nperyear,nmissing(yrbeg:yrend),nn(yrbeg:yrend)
    integer   :: n0,n0missing,j1,j2,y1,y2,nlength(10),length,nperday,mens1,mens,mon
    real      :: frac(yrbeg:yrend)
    real,allocatable :: data(:,:,:)
    logical   :: nomissing
    character :: file*1023,var*80,units*40,line*80
    integer   :: iargc
    integer,external :: leap
    lwrite = .false.
    
    if ( iargc() < 1 .or. iargc() > 7 ) then
        print *,'count_missing file [yr1 yr2] [mon m sel n]'
        call exit(-1)
    end if
    
    call getarg(1,file)
    lstandardunits = .false.
    allocate(data(npermax,yrbeg:yrend,0:nensmax))
    call readensseries(file,data,npermax,yrbeg,yrend,nensmax,nperyear,mens1,mens, &
        var,units,lstandardunits,lwrite)
    
    y1 = yrbeg
    y2 = yrend
    mon = 1
    lsel = nperyear
    if ( iargc() > 2 ) then
        call getarg(2,line)
        if ( line(1:3) == 'mon' .and. iargc() > 2 ) then
            call getopts(2,iargc(),nperyear,yrbeg,yrend,.true.,mens1,mens)
        else
            read(line,*) y1
            call getarg(3,line)
            read(line,*) y2
            call getopts(4,iargc(),nperyear,yrbeg,yrend,.true.,mens1,mens)            
        end if
    end if

    call getj1j2(j1,j2,m1,nperyear,.false.)
    !!!print *,'@@@ m1,lsel,j1,j2 = ',m1,lsel,j1,j2
    call getmissing(data,npermax,yrbeg,yrend,mens1,mens,nperyear,j1,j2,y1,y2,nmissing,nn,n0missing,n0,nlength)
    
    nomissing = .true.
    do yr=y1,y2
        if ( nn(yr) > 0 ) then
            frac(yr) = real(nmissing(yr))/real(nn(yr))
            if ( nmissing(yr) > 0 ) nomissing = .false.
        else
            frac(yr) = 3e33
        end if
    end do
    if ( nomissing ) then
        print '(a)','# no missing data'
    else
        call copyheader(file,6)
        print '(a)','# miss [1] fraction of missing data'
        do yr=y1,y2
            if ( frac(yr) < 1e33 ) then
                print '(i4,f8.4)',yr,frac(yr)
            end if
        end do
        print '(a,2i9,f7.1)','# all ',n0,n0missing,100*real(n0missing)/real(n0)
        print '(a)','# length of gaps'
        do length = 1,10
            print '(a,i2,i9)','#',length,nlength(length)
        end do
        print '(a)','# (the last line includes all longer gaps)'     
    end if
end program