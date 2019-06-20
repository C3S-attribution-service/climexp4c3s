subroutine readsetseries(data,ids,npermax,yrbeg,yrend,nensmax, &
    nperyear,mens1,mens,var,units,lstandardunits,lwrite)
    implicit none
    integer   :: npermax,yrbeg,yrend,nensmax,nperyear,mens1,mens
    real      :: data(npermax,yrbeg:yrend,0:nensmax)
    logical   :: lwrite,lstandardunits
    character :: ids(0:nensmax)*(*),var*(*),units*(*)
    character :: lvar*120,svar*120,history*50000,metadata(2,100)*2000
    call readsetseriesmeta(data,ids,npermax,yrbeg,yrend,nensmax, &
        nperyear,mens1,mens,var,units,lvar,svar,history,metadata,lstandardunits,lwrite)
end subroutine readsetseries

subroutine readsetseriesmeta(data,ids,npermax,yrbeg,yrend,nensmax, &
    nperyear,mens1,mens,var,units,lvar,svar,history,metadata,lstandardunits,lwrite)
!
!   read a set of time series files in a list in the old stupid getstations
!   output format; after a line with coordinates there is a line
!   with after the colon a 'number', which is used by the program
!   (the next argument) to retrieve the data and store it in
!   data/prognumber.dat
!   metadata from the first series is returned.
!
    implicit none
    integer   :: npermax,yrbeg,yrend,nensmax,nperyear,mens1,mens
    real      :: data(npermax,yrbeg:yrend,0:nensmax)
    logical   :: lwrite,lstandardunits
    character :: ids(0:nensmax)*(*),var*(*),units*(*)
    character :: lvar*(*),svar*(*),history*(*),metadata(2,100)*(*)
    integer   :: i,m,iret,filetime,stringtime,ntot
    real      :: lon,lat,elev
    logical   :: lexist,pexist
    character :: prog*100,extraargs*100,extraargs_*100,file*255, &
                 line*255,number*20,string*255,command*1023
    integer,external :: getfiletime

    call get_command_argument(3,prog)
    if ( (prog(1:3) /= 'get' .and. prog(1:3) /= 'eca' .and. &
         prog(1:4) /= 'beca' .and. prog(1:4) /= 'gdcn' .and. &
         prog(1:4) /= 'grid' ) .or. &
         index(prog,'/') /= 0 .or. &
         index(prog,';') /= 0 .or. index(prog,'`') /= 0 .or. &
         index(prog,'&') /= 0 .or. index(prog,'|') /= 0 ) then
        print *,'readsetseries: invalid argument ',prog
        call exit(-1)
    endif
!   one prog has underscores in its name, which conflicts with the convention to put the extra
!   arguments there
    if ( prog(1:4) == 'grid' ) then
        ! the name already contains underscores, so the trick to append the
        ! extraargs separeted by dashes does not work.
        i = 0
    else
        i = index(prog,'_')
    end if
    if ( i /= 0 ) then
        extraargs = prog(i+1:)
        prog = prog(:i-1)
        extraargs_ = extraargs
        do i=1,len(extraargs)
            if ( extraargs(i:i) == '_' ) extraargs(i:i) = ' '
        end do
        write(0,*) 'readsetseries: transforming stations to '// &
            'lower time resolution ',trim(extraargs),'<br>'
    else
        extraargs = ' '
        extraargs_ = ' '
    end if
    call get_command_argument(2,file)
    open(1,file=file,status='old',err=101)
    goto 102
101 continue
    do i=1,len(file)-1
        if ( file(i:i+1) == '\\:' ) then
            file(i:) = file(i+1:)
        end if
    end do
    open(1,file=file,status='old')
102 continue
    mens = 0
    mens1 = 0
    call readntot(1,ntot)
1   continue
    call readcoord(1,lon,lat,elev)
    if ( lon > 1e33 ) goto 2
    read(1,'(a)',err=1,end=1) line
    call readcodename(line,number,string,lwrite)
!   this assumes that it is not an ensemble!
!   first generate original timescale file
    file = './data/'//trim(prog)//trim(number)//'.dat'
    inquire(file=trim(file),exist=lexist)
    inquire(file='./bin/'//trim(prog),exist=pexist)
    if ( lwrite .and. .not.pexist ) print *,'./bin/'//trim(prog)//' does not exist'
!   the getdutch* , gdcn*, eca* and GHCN get* programs
!   check themselves whether the file needs to be regenerated
    if ( pexist .and. &
        prog(1:8) == 'getdutch' .or. &
        prog(1:4) == 'gdcn' .or. &
        prog(1:3) == 'eca' .or. prog(1:4) == 'beca' .or. &
        prog(1:7) == 'getprcp' .or. prog(1:7) == 'gettemp' .or. &
        prog(1:6) == 'getmin' .or. prog(1:6) == 'getmax' .or. &
        prog(1:6) == 'getslp' ) then
        ! generate it
        write(command,'(10a)') './bin/',trim(prog),' ',trim(number),' ',trim(file)
        if ( lwrite ) print *,trim(command),'<br>'
        call mysystem(trim(command),iret)
    else if ( .not. lexist .and. pexist ) then
        ! generate it
        write(command,'(10a)') './bin/',trim(prog), &
            ' ',trim(number),' ',trim(file),' > ',trim(file)
        if ( lwrite ) print *,trim(command),'<br>'
        call mysystem(trim(command),iret)
    end if
    if ( extraargs /= ' ' ) then
    ! generate lower-resolution file
        i = index(file,'.dat')
        string = file(:i-1)//'_'//trim(extraargs_)//'.dat'
        inquire(file=trim(string),exist=lexist)
        if ( lexist ) then
            filetime = getfiletime(trim(file))
            stringtime = getfiletime(trim(string))
        else
            filetime = 0
            stringtime = 0
        end if
        if ( .not. lexist .or. stringtime < filetime ) then
            write(command,'(10a)') './bin/daily2longer ', &
                trim(file),' ',trim(extraargs),' > ',trim(string)
            if ( lwrite ) print *,trim(command),'<br>'
            call mysystem(trim(command),iret)
        end if
        file = string
    end if
    if ( lwrite ) print *,'reading ',trim(file)
    inquire(file=trim(file),exist=lexist)
    if ( .not. lexist ) then
        write(0,'(2a)') 'skipping file ',trim(file)
        goto 1
    end if
    call keepalive1('Reading series',mens,ntot)
    if ( mens == 0 ) then
        call readseriesmeta(file,data(1,yrbeg,mens),npermax,yrbeg,yrend &
            ,m,var,units,lvar,svar,history,metadata,lstandardunits,.false.)
    else
        call readseries(file,data(1,yrbeg,mens),npermax,yrbeg,yrend &
            ,m,var,units,lstandardunits,.false.)
    end if
    ids(mens) = number
    if ( m == 0 ) goto 1
    if ( mens == 0 ) then
        nperyear = m
    else
        if ( m /= nperyear ) then
            write(0,*) 'readsetseries: error: different time scales',m,nperyear, &
                ' in file ',trim(file)
            call exit(-1)
        endif
    endif
    mens = mens + 1
    if ( mens > nensmax ) then
        write(0,*) 'readsetseries: error: increase nensmax'
        call exit(-1)
    endif
    goto 1
2   continue
    mens = mens - 1
end subroutine readsetseriesmeta

