program printbigtable
!
!   print a set of time series as a (big_) ascii table
!
    implicit none
    include 'param.inc'
    integer mensmax
    parameter(mensmax=999)
    integer i,mo,yr,iens,mens1,mens,nperyear,iarg,nens1,nens2
    integer fyr,lyr,datum,d,m
    real,allocatable :: data(:,:,:)
    logical lstandardunits,lset,lwrite,printit
    character file*255,var*20,units*20,string*80,command*500,ids(0:mensmax)*30,format*2
    integer,external :: leap

    if ( command_argument_count().lt.1 ) then
        print *,'usage: printtable ensfile|(file listfile prog) [dummy]'
        call exit(-1)
    endif

    allocate(data(1:npermax,yrbeg:yrend,0:mensmax))

    lstandardunits = .false.
    lwrite = .false.
    call get_command_argument(1,file)
    if ( file == 'file' ) then
        call get_command_argument(2,file)
        lset = .true.
        iarg = 4
    else
        lset = .false.
        iarg = 2
    end if
    nens1 = 0
    nens2 = 999
    do while ( iarg.le.command_argument_count() )
        call get_command_argument(iarg,string)
        if ( string.eq.'ens' ) then
            call get_command_argument(iarg+1,string)
            read(string,*) nens1
            call get_command_argument(iarg+2,string)
            read(string,*) nens2
            iarg = iarg + 3
        else if ( string.eq.'debug' .or. string.eq.'lwrite' ) then
            lwrite = .true.
            print *,'turned on debug output'
            iarg = iarg + 1
        else
            iarg = iarg + 1
        end if
    end do
    if ( lset ) then
        call readsetseries(data,ids,npermax,yrbeg,yrend,mensmax    &
 &           ,nperyear,mens1,mens,var,units,lstandardunits,lwrite)
    else
        call readensseries(file,data,npermax,yrbeg,yrend,mensmax   &
 &           ,nperyear,mens1,mens,var,units,lstandardunits,lwrite)
    end if
    nens1 = max(nens1,mens1)
    nens2 = min(nens2,mens)

    command = '#'
    do i=0,command_argument_count()
        call get_command_argument(i,string)
        command = trim(command)//' '//trim(string(1+index(string,'/',.true.):))
    end do
    print '(a)',trim(command)
    if ( .not.lset ) then
        call filloutens(file,mens1)
        call copyheader(file,6)
        print '(a9,1000i30)','# member ',(iens,iens=nens1,nens2)            
    else
        ! this information is lost by the time I get here...
        print '(10a)','# ',trim(var),' [',trim(units),'] ',trim(var)
        print '(a9,1000a30)','# id     ',(ids(iens),iens=nens1,nens2)    
    end if
    do yr=yrbeg,yrend
        do mo=1,npermax
            printit = .false.
            do iens=nens1,nens2
                if ( data(mo,yr,iens).lt.1e33 ) then
                    printit = .true.
                end if
            end do
            if ( printit ) then
                if ( nperyear == 1 ) then
                    datum = yr
                    format='i4'
                else if ( nperyear == 2 ) then
                    datum = 100*yr + 6*mo - 5
                    format='i6'
                else if ( nperyear == 4 ) then
                    datum = 100*yr + 3*mo - 2
                    format='i6'
                else if ( nperyear == 12 ) then
                    datum = 100*yr + mo
                    format='i6'
                else if ( nperyear == 36 ) then
                    datum = 100*(100*yr + 1+(mo-1)/3) + 10*mod(mo-1,3) + 1
                    format='i8'
                else if ( nperyear >= 360 .and. nperyear <= 366 ) then
                    call getdymo(d,m,mo,nperyear)
                    datum = 100*(100*yr + m) + d
                    format='i8'
                else
                    write(0,*) 'printtable: cannot yet handle nperyear = ',nperyear
                    write(*,*) 'printtable: cannot yet handle nperyear = ',nperyear
                    call exit(-1)
                end if
                print '('//format//',1000g30.7)',datum,(data(mo,yr,iens),iens=nens1,nens2)
            end if
        end do
    end do

end program
