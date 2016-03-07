        subroutine readsetseries(data,ids,npermax,yrbeg,yrend,nensmax,
     +       nperyear,mens1,mens,var,units,lstandardunits,lwrite)
*       read set of data files in a list in the old stupid getstations
*       output format; after a line with coordinates there is a line
*       with after the colon a 'number', which is used by the program
*       (the next argument) to retrieve the data and store it in 
*       data/prognumber.dat
        implicit none
        integer npermax,yrbeg,yrend,nensmax,nperyear,mens1,mens
        real data(npermax,yrbeg:yrend,0:nensmax)
        logical lwrite,lstandardunits
        character ids(0:nensmax)*(*),var*(*),units*(*)
        integer i,m,iret,filetime,stringtime
        real lon,lat,elev
        logical lexist
        character prog*100,extraargs*100,extraargs_*100,file*255,
     +       line*255,number*20,string*255,command*1023
        integer,external :: getfiletime

        call getarg(3,prog)
        if ( (prog(1:3).ne.'get' .and. prog(1:3).ne.'eca' .and. 
     +       prog(1:4).ne.'beca' .and. prog(1:4).ne.'gdcn' .and.
     +       prog(1:4).ne.'grid' ) .or. 
     +        index(prog,'/').ne.0 .or.
     +        index(prog,';').ne.0 .or. index(prog,'`').ne.0 .or.
     +        index(prog,'&').ne.0 .or. index(prog,'|').ne.0 ) then
            print *,'readsetseries: invalid argument ',prog
            call abort
        endif
        ! one prog has underscores in its name, which conflicts with the convention to put the extra 
        ! arguments there
        if ( prog(1:4) == 'grid' ) then
            ! the name already contains underscores, so the trick to append the 
            ! extraargs separeted by dashes does not work.       
            i = 0
        else
            i = index(prog,'_')
        end if
        if ( i.ne.0 ) then
            extraargs = prog(i+1:)
            prog = prog(:i-1)
            extraargs_ = extraargs
            do i=1,len(extraargs)
                if ( extraargs(i:i).eq.'_' ) 
     +               extraargs(i:i) = ' '
            end do
            write(0,*) 'readsetseries: transforming stations to '//
     +           'lower time resolution ',trim(extraargs),'<br>'
        else
            extraargs = ' '
            extraargs_ = ' '
        end if
        call getarg(2,file)
        open(1,file=file,status='old')
        mens = 0
        mens1 = 0
 1      continue
        call readcoord(1,lon,lat,elev)
        if ( lon.gt.1e33 ) goto 2
        read(1,'(a)',err=1,end=1) line
        call readcodename(line,number,string,lwrite)
!       this assumes that it is not an ensemble!
!       first generate original timescale file
        file = './data/'//trim(prog)//trim(number)//'.dat'
        inquire(file=trim(file),exist=lexist)
        ! the getdutch* , gdcn*, eca* and GHCN get* programs 
        ! check themselves whether the file needs to be regenerated
        if ( prog(1:8).eq.'getdutch' .or.
     +       prog(1:4).eq.'gdcn' .or.
     +       prog(1:3).eq.'eca' .or. prog(1:4).eq.'beca' .or.
     +       prog(1:7).eq.'getprcp' .or. prog(1:7).eq.'gettemp' .or.
     +       prog(1:6).eq.'getmin' .or. prog(1:6).eq.'getmax' .or.
     +       prog(1:6).eq.'getslp' ) then
            ! generate it
            write(command,'(10a)') './bin/',trim(prog),
     +           ' ',trim(number),' ',trim(file)
            if ( lwrite ) print *,trim(command),'<br>'
            call mysystem(trim(command),iret)
        else if ( .not.lexist ) then
            ! generate it
            write(command,'(10a)') './bin/',trim(prog),
     +           ' ',trim(number),' ',trim(file),' > ',trim(file)
            if ( lwrite ) print *,trim(command),'<br>'
            call mysystem(trim(command),iret)
        end if
        if ( extraargs.ne.' ' ) then
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
            if ( .not.lexist .or. stringtime.lt.filetime ) then
                write(command,'(10a)') './bin/daily2longer ',
     +               trim(file),' ',trim(extraargs),' > ',trim(string)
                if ( lwrite ) print *,trim(command),'<br>'
                call mysystem(trim(command),iret)
            end if
            file = string
        end if
        if ( lwrite ) print *,'reading ',trim(file)
        call keepalive1('Reading series',mens,-1)
        call readseries(file,data(1,yrbeg,mens),npermax,yrbeg,yrend
     +       ,m,var,units,lstandardunits,.false.)
        ids(mens) = number
        if ( m.eq.0 ) goto 1
        if ( mens.eq.0 ) then
            nperyear = m
        else
            if ( m.ne.nperyear ) then
                write(0,*) 'readsetseries: error: different time scales'
     +               ,m,nperyear
                call abort
            endif
        endif
        mens = mens + 1
        if ( mens.gt.nensmax ) then
            write(0,*) 'readsetseries: error: increase nensmax'
            call abort
        endif
        goto 1
 2      continue
        mens = mens - 1
        
        end subroutine

