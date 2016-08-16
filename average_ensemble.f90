program average_ensemble
!
!   compute the ensemble average of a set of time series
!
    implicit none
    include 'param.inc'
    integer mensmax
    parameter(mensmax=1999)
    integer i,mo,yr,iens,mens1,mens,nperyear,iarg,nens1,nens2
    integer fyr,lyr,nmissing(yrbeg:yrend),nnn(yrbeg:yrend),n0missing,n0,nlength(10)
    integer,allocatable :: nn(:,:)
    real,allocatable :: data(:,:,:),mean(:,:)
    logical lstandardunits,lset,lwrite
    character file*255,var*20,units*20,string*80,oper*3,ids(0:mensmax)*30,command*500
    integer,external :: leap
    integer iargc

    if ( iargc().lt.2 ) then
        print *,'usage: average_ensemble ensfile|(file setfile) [ens n1 n2] '// &
 &           '[mean|min|max|num] [debug] [dummy]'
        call exit(-1)
    endif

    allocate(data(1:npermax,yrbeg:yrend,0:mensmax))
    allocate(mean(1:npermax,yrbeg:yrend))
    allocate(nn(1:npermax,yrbeg:yrend))

    lstandardunits = .false.
    lwrite = .false.
    call getarg(1,file)
    if ( file == 'file' ) then
        call getarg(2,file)
        lset = .true.
        iarg = 4
    else
        lset = .false.
        iarg = 2
    end if
    nens1 = 0
    nens2 = mensmax
    oper = 'mea'
    do while ( iarg.le.iargc() )
        call getarg(iarg,string)
        if ( string.eq.'ens' ) then
            call getarg(iarg+1,string)
            read(string,*) nens1
            call getarg(iarg+2,string)
            read(string,*) nens2
            iarg = iarg + 3
        else if ( string.eq.'mean' ) then
            oper = 'mea'
            iarg = iarg + 1
        else if ( string.eq.'min' ) then
            oper = 'min'
            iarg = iarg + 1
        else if ( string.eq.'max' ) then
            oper = 'max'
            iarg = iarg + 1
        else if ( string.eq.'num' ) then
            oper = 'num'
            iarg = iarg + 1
        else if ( string.eq.'mis' ) then
            oper = 'mis'
            iarg = iarg + 1
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

    nn(1:nperyear,yrbeg:yrend) = 0
    if ( oper.eq.'mea' .or. oper.eq.'num' ) then
        mean(1:nperyear,yrbeg:yrend) = 0
    else if ( oper.eq.'min' ) then
        mean(1:nperyear,yrbeg:yrend) = 3e33
    else if ( oper.eq.'max' ) then
        mean(1:nperyear,yrbeg:yrend) = -3e33
    else if ( oper.eq.'mis' ) then
        fyr = yrbeg
        lyr = yrend
        call getmissing(data,npermax,yrbeg,yrend,nens1,nens2, &
        &   nperyear,fyr,lyr,nmissing,nnn,n0missing,n0,nlength)
    else
        write(0,*) 'average_ensemble: error: unknown operation ',oper
        call exit(-1)
    end if
    if ( oper /= 'mis' ) then
        do iens=nens1,nens2
            do yr=yrbeg,yrend
                do mo=1,nperyear
                    if ( data(mo,yr,iens).lt.1e33 ) then
                        nn(mo,yr) = nn(mo,yr) + 1
                        if ( oper.eq.'mea' ) then
                            mean(mo,yr) = mean(mo,yr) + data(mo,yr,iens)
                        else if ( oper.eq.'min' ) then
                            mean(mo,yr) = min(mean(mo,yr),data(mo,yr,iens))
                        else if ( oper.eq.'max' ) then
                            mean(mo,yr) = max(mean(mo,yr),data(mo,yr,iens))
                        else if ( oper.eq.'num' ) then
                            mean(mo,yr) = mean(mo,yr) + 1
                        else
                            write(0,*) 'average_ensemble: error: '//    &
     &                           'unknown operation 2 ',oper
                            call exit(-1)
                        end if
                    endif
                enddo
            enddo
        enddo
        if ( oper.eq.'num' ) then
            do fyr=yrbeg,yrend
                do mo=1,nperyear
                    if ( nn(mo,fyr).gt.0 ) goto 110
                end do
            end do
110         continue
            do lyr=yrend,yrbeg,-1
                do mo=1,nperyear
                    if ( nn(mo,lyr).gt.0 ) goto 120
                end do
            end do
120         continue
            if ( fyr.gt.yrbeg ) mean(:,yrbeg:fyr-1) = 3e33
            if ( lyr.lt.yrend ) mean(:,lyr+1:yrend) = 3e33
            if ( nperyear.eq.366 ) then
                do yr=fyr,lyr
                    if ( leap(yr).eq.1 ) then
                        if ( mean(60,yr).eq.0 ) then
                            mean(60,yr) = 3e33
                        else
                            write(0,*) 'average_ensemble: suspicious'//    &
     &                           ' value on 29 Feb ',yr,mean(60,yr)
                        end if
                    end if
                end do
            end if
        else
            do yr=yrbeg,yrend
                do mo=1,nperyear
                    if ( nn(mo,yr).gt.1 ) then
                        if ( oper.eq.'mea' ) then
                            mean(mo,yr) = mean(mo,yr)/nn(mo,yr)
                        end if
                    else
                        mean(mo,yr) = 3e33
                    endif
                enddo
            enddo
        end if
    else ! oper == 'mis'
        nperyear = 1
        mean = 3e33
        do yr=fyr,lyr
            if ( nnn(yr) > 0 ) then
                mean(1,yr) = real(nmissing(yr))/real(nnn(yr))
            end if
        end do
    end if
    command = '#'
    do i=0,iargc()
        call getarg(i,string)
        command = trim(command)//' '//trim(string(1+index(string,'/',.true.):))
    end do
    print '(a)',trim(command)
    print '(4a)','# ensemble ',oper,' of file ',trim(file)
    print '(a,i3,a,i3)','# using members ',nens1,' to ',nens2
    if ( oper.eq.'num' ) then
        print '(a)','# N [1] number of series with valid data'
    else if ( oper.eq.'mis' ) then
        print '(a)','# miss [1] fraction with invalid data'
    else if ( .not.lset ) then
        call filloutens(file,mens1)
        call copyheader(file,6)
    else
        ! this information is lost by the time I get here...
        print '(10a)','# ',trim(var),' [',trim(units),'] ',    &
 &           trim(oper),' of ',trim(var)
    end if
    call printdatfile(6,mean,npermax,nperyear,yrbeg,yrend)

end program
