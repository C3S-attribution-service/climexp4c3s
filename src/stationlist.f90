    program stationlist
!
!   plough through a statonlist as composed by the get* programs
!   and make a list suitable for feeding to a plotting program:
!
!   # lon1 lon2 lat1 lat2
!   code lon lat value name
!   code lon lat value name
!
!   should do this in python...
!
    implicit none
    integer,parameter :: nmax=250000
    integer :: i,j,k,n,ntot,ldir,nextarg,retval,nsigns,isign,loptions
    real :: lon,lon1,lon2,lat,lat1,lat2,val,sign,regr,x1,s1,x2,s2,z, &
        signs(nmax),psigns(5),mean,sd,sd1,skew,xmin,xmax,chsq,prob, &
        tx,tx25,tx975,t(10),t25(10),t975(10),alpha,alpha25,alpha975 &
        ,ttx(3),ttx25(3),ttx975(3),tt(10,3),tt25(10,3),tt975(10,3)
    character(20) :: code
    real :: elev
    logical :: lexist,first,ltwosided
    character line*255,file*256,prog*120,filterargs*100, &
        filterext*100,oper*4,name*40,options*1024,datfile*256 &
        ,rawfile*256,command*1024,tmpfile*29,aveline*10000,dir*256 &
        ,coordlist*30,field*255,fieldfile*255,letter*1,runfile*256 &
        ,histoptions*80,email*80,attribute_args*1024,country*40 &
        ,pid*20
    logical :: lwrite,interp
    integer :: getpid
    real :: erfc
    data psigns /0.1,0.05,0.01,0.005,0.001/

    lwrite = .false. 
    nsigns = 0
    ltwosided = .false.
    if ( command_argument_count() < 4 ) then
        print *,'usage: stationlist listfile outfile prog operation ...'
        call exit(-1)
    end if
    call getenv('EMAIL',email)
    call get_command_argument(4,line)
    if ( line(1:6) == 'ifield' ) then
        oper = line(7:)
        interp = .true. 
        letter = 'i'
        call get_command_argument(5,field)
        nextarg = 6
    elseif ( line(1:5) == 'field' ) then
        oper = line(6:)
        interp = .false. 
        letter = 'n'
        call get_command_argument(5,field)
        nextarg = 6
    else
        oper = line
        field = ' '
        nextarg = 5
    end if
    if ( lwrite ) print *,'stationlist: oper = ',oper
    call get_command_argument(3,prog)
    if ( (prog(1:3) /= 'get' .and. prog(1:3) /= 'eca' .and. &
          prog(1:4) /= 'beca' .and. prog(1:4) /= 'gdcn' .and. &
          prog(1:4) /= 'grid' ) .or. &
        index(prog,'/') /= 0 .or. &
        index(prog,';') /= 0 .or. index(prog,'`') /= 0 .or. &
        index(prog,'&') /= 0 .or. index(prog,'|') /= 0 ) then
        print *,'stationlist: invalid argument ',prog
        call exit(-1)
    end if
    if ( lwrite ) print *,'stationlist: prog = ',trim(prog)
    i = index(prog,'_')
    if ( i /= 0 ) then
!       the output has to be filtered
        filterext = prog(i:)
        prog(i:) = ' '
        filterargs = filterext
        do i=1,len_trim(filterargs)
            if ( filterargs(i:i) == '_' ) filterargs(i:i) = ' '
        end do
    else
        filterext = ' '
    end if
    call get_command_argument(2,file)
    open(2,file=trim(file),status='unknown',err=900)
    call get_command_argument(1,file)
    open(1,file=trim(file),status='old',err=1)
    goto 2
1   continue
!   try again ith backslashes removed; Mac OS and linux differ in how
!   they propagate arguments to cgi scripts. Some backslashes may have crept in...
    do i=1,len(file)
        if ( file(i:i) == '\\' ) then
            file(I:) = file(i+1:)
        end if
    end do
    open(1,file=trim(file),status='old',err=901)
2   continue
    j=1
    options = ' '
    do i=nextarg,command_argument_count()
        call get_command_argument(i,options(j:))
        if ( options(j:j+3) == 'runc' .or. options(j:j+3) == 'runr' ) then
            call get_command_argument(i+2,runfile)
        end if
        if ( options(j:j+4) == 'debug' ) then
            lwrite = .true. 
        end if
        j = min(len(options)-10,len_trim(options)+2)
    end do
    if ( oper(1:2) == 'hi' ) then
        if ( oper(3:4) == 'gr' .or. oper(3:4) == 'gR' ) then
            histoptions = '20 hist sqrtlog '
        else if ( oper(3:4) == 'pr' .or. oper(3:4) == 'pR' ) then
            histoptions = '20 hist log '
        else if ( oper(3:4) == 'vr' .or. oper(3:4) == 'vR' ) then
            histoptions = '20 hist gumbel '
        else
            histoptions = '20 hist hist '
        end if
        if ( oper(1:3) == 'hig' ) then
            histoptions(len_trim(histoptions)+2:) = 'fit gauss'
        elseif ( oper(1:3) == 'hip' ) then
            histoptions(len_trim(histoptions)+2:) = 'fit gpd'
        elseif ( oper(1:3) == 'hiv' ) then
            histoptions(len_trim(histoptions)+2:) = 'fit gev'
        end if
        if ( lwrite ) print *,'histoptions are ',histoptions(1:len_trim(histoptions))
    else if ( oper(1:2) == 'at' ) then
        call getenv('attribute_args',attribute_args)
    end if
    write(pid,'(i20.20)') getpid()
    tmpfile = '/tmp/corr'//pid
    options(len_trim(options)+2:) = 'plot '//tmpfile
    loptions = len_trim(options)
    call getenv('DIR',dir)
    if ( dir == ' ' ) dir = '.'
    ldir = len_trim(dir)
    if ( lwrite ) print *,'DIR = ',dir(1:max(1,ldir))
    if ( ldir == 0 ) then
        call getenv('HOME',dir)
        ldir = len_trim(dir)
        dir = dir(1:len_trim(dir))//'/climexp/'
        ldir = len_trim(dir)
        if ( lwrite ) print *,'DIR = ',dir(1:max(1,ldir))
    elseif ( dir(ldir:ldir) /= '/' ) then
        ldir = ldir + 1
        dir(ldir:ldir) = '/'
        if ( lwrite ) print *,'DIR = ',dir(1:max(1,ldir))
    end if

!   write own PID in kill metafile

    call killfile(command,line,datfile,0)

!   for field-station correlations generate a list of coordinates

    if ( field /= ' ' ) then
        write(coordlist,'(a,i20.20)') '/tmp/coord',getpid()
        open(3,file=coordlist)
        write(3,'(a)') trim(dir)//'data'
    10  continue
        call readcoord(1,lon,lat,elev)
        if ( lon > 1e33 ) goto 20
        write(3,'(2f8.2)') lon,lat
        goto 10
    20  continue
        close(3)
        rewind(1)
    
!       and generate the required time series files
    
        command = 'bin/get_index '//trim(field)//' file '//coordlist
        if ( interp ) then
            command(len_trim(command)+1:) = ' interpolate'
        end if
        if ( .false. ) then
            command(len_trim(command)+1:) = ' debug'
        end if
        if ( lwrite ) print *,'executing ',trim(command)
        call mysystem(command,retval)
        call checkval(retval,command)
    end if

!   search for box

100 continue
    read(1,'(a)',end=902,err=902) line
    i = index(line,'earching for stations') + &
        index(line,'ocated stations') + index(line,'stations in') + &
        index(line,'regions in') + index(line,'rid points in') + &
        index(line,'Found')*index(line,' stations')
    if ( i == 0 ) goto 100
    if ( index(line,'inside polygon') /= 0 ) goto 100

    i = index(line,'ions in') + index(line,'ints in')
    if ( i /= 0 ) then
!       found box
        i = i + 7
        j = i + index(line(i:),'N') - 2
        read(line(i:j),*) lat1
        i = i + index(line(i:),':')
        j = i + index(line(i:),'N') - 2
        read(line(i:j),*) lat2
        i = i + index(line(i:),',')
        j = i + index(line(i:),'E') - 2
        read(line(i:j),*) lon1
        i = i + index(line(i:),':')
        j = i + index(line(i:),'E') - 2
        read(line(i:j),*) lon2
    else
        i = index(line,'stations near')
        if ( i /= 0 ) then
!           found circle
            call readcoord(1,lon1,lat1,elev)
            if ( lon1 > 1e33 ) goto 180
            lon2 = lon1
            lat2 = lat1
        110 continue
            call readcoord(1,lon,lat,elev)
            if ( lon > 1e33 ) goto 180
            if ( abs(lon1-lon) < abs(lon1-lon-360) .and. &
                 abs(lon1-lon) < abs(lon1-lon+360) ) then
                lon1 = min(lon1,lon)
            elseif ( abs(lon1-lon-360) < abs(lon1-lon+360) ) then
                lon1 = min(lon1,lon+360)
            else
                lon1 = min(lon1,lon-360)
            end if
            if ( abs(lon2-lon) < abs(lon2-lon-360) .and. &
                 abs(lon2-lon) < abs(lon2-lon+360) ) then
                lon2 = max(lon2,lon)
            elseif ( abs(lon2-lon-360) < abs(lon2-lon+360) ) then
                lon1 = max(lon2,lon+360)
            else
                lon2 = max(lon2,lon-360)
            end if
            lat1 = min(lat1,lat)
            lat2 = max(lat2,lat)
            goto 110
        180 continue
                            
!           extra 10% on all sides
            lon1 = lon1 - (lon2-lon1)/10
            lon2 = lon2 + (lon2-lon1)/11
            lat1 = max(-90.,lat1 - (lat2-lat1)/10)
            lat2 = min(+90.,lat2 + (lat2-lat1)/11)
            rewind(1)
        else
!           no co-ordinates - take the whole world
            lon1 = -30
            lon2 = 330
            lat1 = -90
            lat2 =  90
        end if
    end if
    if ( lon2-lon1 == 360 ) then
        lon1 = -30
        lon2 = 330
    end if
    write(2,'(a,4f10.4)') '# ',lon1,lon2,lat1,lat2
!
!   search for number of stations
!
    rewind(1)
    call readntot(1,ntot)
!
!   get station data
!
    first = .true. 
    if ( oper(1:2) == 'hi' .or. oper(1:2) == 'at' ) then
        call mysystem('echo "<table class=realtable width="100%" '// &
            'border=0 cellpadding= cellspacing=0>"',retval)
    end if
    if ( oper == 'aver' ) then
        aveline = trim(dir)//'bin/averageseries'
    end if
    n = 0
200 continue
    call keepalive1('Station ',n,ntot)
    call readcountry(1,country)
    if ( lwrite ) print *,'stationlist: attempting to read station ',n+1
    call readcoord(1,lon,lat,elev)
    if ( lon > 1e33 ) goto 800
    read(1,'(a)',end=800,err=903) line
    call readcodename(line,code,name,lwrite)

!   operate

    sign = -1
    if ( oper == 'plot' ) then
        val = 0.2
        sign = 0
    elseif ( oper == 'list' ) then
        val = elev
        sign = 0
    elseif ( oper == 'corr' .or. oper == 'sign' .or. oper == 'slop' &
         .or. oper == 'nslo' .or. oper == 'regr' .or. &
        oper == 'nreg' .or. oper == 'aver' .or. &
        oper == 'auco' .or. oper == 'ausl' .or. oper == 'aure' &
         .or. oper == 'val' .or. oper == 'anom' .or. &
        oper == 'frac' .or. oper == 'zval' .or. &
        oper == 'runc' .or. oper == 'zdif' .or. &
        oper == 'bdif' .or. &
        oper == 'hime' .or. oper == 'hisd' .or. &
        oper == 'hisk' .or. &
        oper == 'himi' .or. oper == 'hima' .or. &
        oper(1:3) == 'hig' .or. &
        oper(1:3) == 'hip' .or. oper(1:3) == 'hiv' .or. &
        oper(1:2) == 'at' ) then
        if ( email == 'someone@somewhere' .and. n >= 100 ) then
            print *,'anonymous users can only use 100 stations'
            goto 800
        end if
!       retrieve data
        write(datfile,'(4a)') trim(dir)//'data/',trim(prog),trim(code),'.dat'
        do i=1,len_trim(datfile)
            if ( datfile(i:i) == ' ' ) datfile(i:i) = '_'
        end do
        ! check whether the file exists and has length >0
        open(3,file=datfile,err=300)
        read(3,'(a)',err=300,end=300) line
        close(3)
        goto 310
    300 continue
        close(3)
        if ( prog(1:4) == 'grid' ) then
            write(*,*) 'stationlist: error: the grid point ',trim(datfile) &
                ,' should have been generated already'
            write(0,*) 'stationlist: error: the grid point ',trim(datfile) &
                ,' should have been generated already'
            call exit(-1)
        end if
        if ( prog(1:8) == 'getdutch' .or. prog(1:4) == 'gdcn' .or. &
             prog(1:3) == 'eca' .or. prog(1:4) == 'beca' .or. &
             prog(1:7) == 'getprcp' .or. prog(1:7) == 'gettemp' .or. &
             prog(1:6) == 'getmin' .or. prog(1:6) == 'getmax' .or. &
             prog(1:6) == 'getslp' .or. prog(1:9) == 'getfrench' .or. &
             prog(1:10) == 'getaustria' ) then
            write(command,'(6a)') trim(dir)//'bin/',trim(prog), &
                ' ',trim(code),' ',trim(datfile)//'.'//pid
        else
            write(command,'(6a)') trim(dir)//'bin/',trim(prog), &
                ' ',trim(code),' > ',trim(datfile)//'.'//pid
        end if
        if ( lwrite ) print *,trim(command)
        call mysystem(command,retval)
        !!!call checkval(retval,command) do not mind a few going wrong...
        if ( retval /= 0 ) then
            call mysystem('rm -f '//trim(datfile)//'.'//pid,retval)
            goto 200
        end if
        command = 'mv '//trim(datfile)//'.'//pid//' '//trim(datfile)
        if ( lwrite ) print *,trim(command)
        call mysystem(command,retval)
        call checkval(retval,command)
        inquire(file=trim(datfile)//'.'//pid//'.nc',exist=lexist)
        if ( lexist ) then
            command = 'mv '//trim(datfile)//'.'//pid//'.nc '//datfile(1:index(datfile,'.dat')-1)//'.nc'
            if ( lwrite ) print *,trim(command)
            call mysystem(command,retval)
            call checkval(retval,command)            
        end if
    310 continue
        if ( filterext /= ' ' ) then
            ! filter datfile to longer-period file
            rawfile = datfile
            i = index(rawfile,'.dat')
            datfile = rawfile(:i-1)//trim(filterext)//'.dat'
            open(3,file=datfile,err=400)
            read(3,'(a)',err=400,end=400) line
            close(3)
            goto 410
        400 continue
            close(3)
            command = trim(dir)//'bin/daily2longer '// &
                trim(rawfile)//trim(filterargs)//' > ' &
                //trim(datfile)//'.'//pid
            if ( lwrite ) print *, trim(command)
            call mysystem(command,retval)
            call checkval(retval,command)
            command = 'mv '//trim(datfile)//'.'//pid//' '//trim(datfile)
            call mysystem(command,retval)
            call checkval(retval,command)            
        410 continue
        end if
        if ( oper == 'aver' ) then
            i = len_trim(aveline)
            if ( i > len(aveline)-30 ) then
                print *,'stationlist: error: too many stations for average'
                call exit(-1)
            end if
            aveline(i+2:) = datfile
        else
            ! correlate or plot data
            n = n + 1
            write(line,'(i6)') n
            if ( ntot > 0 ) then
                line(7:8) = '_/'
                write(line(9:),'(i6)') ntot
                do i=9,14
                    if ( line(i:i) == ' ' ) line(i:i) = '_'
                end do
            end if
            do i=1,len_trim(name)
                if ( name(i:i) == ' ' ) name(i:i) = '_'
            end do
            options(loptions+2:) = 'name '//trim(line)//'_'//trim(name)
            if ( oper == 'val' .or. oper == 'frac' .or. &
                 oper == 'anom' .or. oper == 'zval' ) then
                ! plot data at a given time
                write(command,'(6a)') trim(dir)//'/bin/getval ',trim(datfile),' ',trim(options)
                if ( lwrite ) then
                    print *,'stationlist: plot data'
                    print *,'             command = ',trim(command),'<p>'
                end if
            elseif ( field /= ' ' ) then
                do i=len_trim(field),1,-1
                    if ( field(i:i) == '/' ) goto 420
                end do
            420 continue
                write(fieldfile,'(3a,f7.2,a,f6.2,3a)') &
                    'data/grid',field(i+1:len_trim(field)),'_',lon,'_', &
                    lat,'_',letter,'.dat'
                do i=1,len_trim(fieldfile)
                    if ( fieldfile(i:i) == ' ' ) fieldfile(i:i)='0'
                end do
                write(command,'(6a)') trim(dir)//'bin/correlate ', &
                    trim(datfile),' file ',trim(fieldfile),' ',trim(options)
                if ( lwrite ) then
                    print *,'stationlist: correlate with field'
                    print *,'             command = ',trim(command)
                end if
            elseif ( oper(1:2) == 'au' ) then
                write(command,'(6a)') trim(dir)//'bin/correlate ', &
                    trim(datfile),' file ',trim(datfile),' ',trim(options)
                if ( lwrite ) then
                    print *,'stationlist: autocorrelation'
                    print *,'             command = ',trim(command)
                end if
            elseif ( oper(1:2) == 'hi' ) then
                write(command,'(6a)') trim(dir)//'bin/histogram ', &
                    trim(datfile),' ',trim(histoptions),' ', &
                    trim(options)
                if ( lwrite ) then
                    print *,'stationlist: histogram'
                    print *,'             command = ',trim(command)
                end if
            elseif ( oper(1:2) == 'at' ) then
                write(command,'(6a)') trim(dir)//'bin/attribute ', &
                    trim(datfile),' ',trim(attribute_args),' ', trim(options)
                if ( lwrite ) then
                    print *,'stationlist: attribute'
                    print *,'             command = ',trim(command)
                end if
            else
                write(command,'(4a)') trim(dir)//'bin/correlate ', &
                    trim(datfile),' ',trim(options)
                if ( lwrite ) then
                    print *,'stationlist: correlate (I think)'
                    print *,'             command = ',trim(command)
                end if
            end if
            command(len_trim(command)+2:)= '|tr -c -d "\\n\\r[:print:]"'
            if ( oper(1:2) == 'hi' .or. oper(1:2) == 'at' ) then
                if ( .not. first ) then
                    if ( oper == 'higa' ) then
                        command(len_trim(command)+2:) = '|egrep  ''(&chi;|th colspan)'''
                    elseif ( oper(1:3) == 'hig' .or. &
                             oper(1:3) == 'hip' .or. &
                             oper(1:3) == 'hiv' ) then
                        command(len_trim(command)+2:) = &
                        '|egrep  ''(return perio|th colspan)'''
                    elseif ( oper == 'hisd' ) then
                        command(len_trim(command)+2:) = '|egrep  '// &
                            '''(^# <tr><td>s.d..n-1|th colspan)'''
                    elseif ( oper(1:2) == 'at' ) then ! print the var name in a HTML comment
                        command(len_trim(command)+2:) = &
                            '|egrep ''('//oper//'|th colspan)'''
                    else
                        command(len_trim(command)+2:) = &
                            '|egrep  ''(^# <tr><td>'//oper(3:4) &
                            //'|th colspan)'''
                    end if
                else
                    command(len_trim(command)+2:) = '|egrep  ''(^# <)'''
                end if
                command(len_trim(command)+2:) = '| sed -e "s/# //"'
            else
                if ( .not. first ) then
                    command(len_trim(command)+2:) = &
                        '|sed -e ''s/# \\([bz]dif\\)/\\1/'''// &
                        '|egrep -v ''(#|===|sign.  no|% pnts|'// &
                        'Month.*:|All year|Requiring)'''
                end if
            end if
            if ( lwrite ) print *,trim(command),'<p>'
            first = .false. 
            call mysystem(command,retval)
            ! do not call checkval(retval,command)
            ! as the return value is the one of grep...
            ! retrieve results
            if ( oper == 'zdif' .or. oper == 'bdif' .or. oper == 'runc' ) then
                if ( lwrite ) print *,'opening file ',trim(runfile)
                open(3,file=runfile,status='old')
            500 continue
                read(3,'(a)',end=510) line
                goto 520
            510 continue
                call mysystem('echo "error computing running correlations"',retval)
                goto 200
            520 continue
                if ( line(1:6) /= '# zdif' .and. &
                    line(1:6) /= '# bdif' ) goto 500
                if ( index(line,'****') /= 0 ) then
                    call mysystem('echo "not enough data"',retval)
                    goto 200
                end if
                read(line,'(9x,f6.2,3x,e14.6)',err=904) val,sign
                if ( lwrite ) print *,'read val,sign = ',val,sign
                close(3,status='delete')
            elseif ( oper(1:2) == 'hi' ) then
                if ( lwrite ) print *,'sttaionlist: opening file ',trim(tmpfile)
                open(3,file=tmpfile,status='old',err=200)
                read(3,*,end=200,err=200) mean
                read(3,*,end=200,err=200) sd
                read(3,*,end=200,err=200) sd1
                read(3,*,end=200,err=200) skew
                read(3,*,end=200,err=200) xmin
                read(3,*,end=200,err=200) xmax
                if ( oper == 'higa' ) then
                    read(3,*,end=200,err=200) chsq
                    read(3,*,end=200,err=200) prob
                elseif ( oper(1:3) == 'hig' .or. &
                    oper(1:3) == 'hip' .or. &
                    oper(1:3) == 'hiv' ) then
                    do i=1,10
                        read(3,*,end=200,err=200)t(i),t25(i),t975(i)
                    end do
                    read(3,*,end=200,err=200) tx,tx25,tx975
                end if
                if ( lwrite ) then
                    print *,'mean = ',mean
                    print *,'sd   = ',sd
                    print *,'sd1  = ',sd1
                    print *,'skew = ',skew
                    print *,'xmin = ',xmin
                    print *,'xmax = ',xmax
                    print *,'chsq = ',chsq
                    print *,'prob = ',prob
                    do i=1,10
                        print *,'t(',i,') = ',t(i),t25(i),t975(i)
                    end do
                    print *,'tx   = ',tx,tx25,tx975
                    close(3)
                else
                    close(3,status='delete')
                end if
                if ( oper == 'hime' ) then
                    val = mean
                elseif ( oper == 'hisd' ) then
                    val = sd1
                elseif ( oper == 'hisk' ) then
                    val = skew
                elseif ( oper == 'himi' ) then
                    val = xmin
                elseif ( oper == 'hima' ) then
                    val = xmax
                elseif ( oper == 'higa' ) then
                    val = chsq
                    sign = prob
                elseif ( oper == 'higr' .or. oper == 'hipr' .or. &
                    oper == 'hivr' ) then
                    val = tx
                    if ( tx < 1e33 ) sign = 1/tx
                elseif ( oper == 'higR' .or. oper == 'hipR' .or. &
                    oper == 'hivR' ) then
                    val = tx25
                elseif ( oper(1:3) == 'hig' .or. oper(1:3) == 'hip' &
                     .or. oper(1:3) == 'hiv' &
                    ) then
                    read(oper(4:4),'(i1)') i
                    val = t(i+1)
                else
                    write(*,*) 'error: unknown operation ',oper
                    write(0,*) 'error: unknown operation ',oper
                    call exit(-1)
                end if
            elseif ( oper(1:2) == 'at' ) then
                if ( lwrite ) print *,'stationlist: opening file ',trim(tmpfile)
                open(3,file=trim(tmpfile),status='old',err=200)
                do j=1,3
                    read(3,*,end=200,err=200) ttx(j),ttx25(j),ttx975(j)
                end do
                read(3,*,end=200,err=200) sign
                do i=1,10
                    do j=1,2
                        read(3,*,end=200,err=200) tt(i,j),tt25(i,j),tt975(i,j)
                    end do
                end do
                read(3,*,end=200,err=200) alpha,alpha25,alpha975
                if ( lwrite ) then
                    do j=1,3
                        print *,'tx(',j,')   = ',ttx(j),ttx25(j),ttx975(j)
                    end do
                    print *,'sign = ',sign
                    do i=1,10
                        do j=1,2
                            print *,'tt(',i,j,') = ',tt(i,j),tt25(i,j),tt975(i,j)
                        end do
                    end do
                    print *,'alpha = ',alpha,alpha25,alpha975
                    close(3)
                else
                    close(3,status='delete')
                end if
                if ( oper == 'atal' ) then
                    val = alpha
                elseif ( oper == 'atr2' ) then
                    val = ttx(1)
                    sign = 1/ttx(1)
                elseif ( oper == 'atr1' ) then
                    val = ttx(2)
                    sign = 1/ttx(2)
                elseif ( oper == 'atra' ) then
                    val = log10(ttx(3))
                    sign = 2*sign ! the output is one-sided, here we use two-sided
                    if ( .false. ) then ! old code
                        ! assuming the PDF of log(ratio) is gaussian the 95% CI corresponds to ±2\sigma
                        sd = (log(ttx975(3)) - log(ttx25(3)))/4
                        z = abs(log(ttx(3))/sd)
                        ! convert to p-value, assuming a 2-sided test, even though we
                        ! often know from theory which way it should go.
                        sign = erfc(z/sqrt(2.))
                    end if
                elseif ( ichar(oper(3:3)) >= ichar('0') .and. &
                         ichar(oper(3:3)) <= ichar('9') ) then
                    read(oper(3:3),'(i1)') i
                    read(oper(4:4),'(i1)') j
                    j = 3-j ! stupid clash of conventions
                    val = tt(i+1,j)
                else
                    write(*,*) 'error: unknown operation ',oper
                    write(0,*) 'error: unknown operation ',oper
                    call exit(-1)
                end if
            else            ! correlation or getval
                open(3,file=tmpfile,status='old',err=200)
                read(3,*,end=700,err=700) i,j,val,sign,k,x1,s1,x2,s2,regr
                close(3,status='delete')
            end if
            if ( val > 1e33 ) goto 200
            if ( oper == 'sign' .or. oper == 'ausi' .or. oper == 'runc' ) then
                if ( sign == 0 ) sign=1e-35
                val = -log10(sign)
            elseif ( oper == 'slop' .or. oper == 'ausl' .or. &
                     oper == 'regr' .or. oper == 'aure' ) then
                val = regr
            elseif ( oper == 'nslo' .or. oper == 'nreg' ) then
                if ( x1 == 0 ) then
                    if ( regr == 0 ) then
                        val = 0
                    else
                        val = 3e33
                    end if
                else
                    val = regr/x1
                end if
            elseif ( oper == 'frac' ) then
                if ( x1 == 0 ) then
                    if ( val == 0 ) then
                        val = 0
                    else
                        val = 3e33
                    end if
                else
                    val = val/x1-1
                end if
                sign = 0
            elseif ( oper == 'val' ) then
                sign = 0
            elseif ( oper == 'anom' ) then
                val = val-x1
                sign = 0
            elseif ( oper == 'zval' ) then
                val = (val-x1)/s1
                sign = 0
            end if
            if ( sign >= 0 ) then
                if ( nsigns < nmax ) then
                    nsigns = nsigns + 1
                    if ( val < 0 ) then
                        ltwosided = .true.
                        signs(nsigns) = 1-sign
                    else
                        signs(nsigns) = sign
                    end if
                else
                    print *,'last stations not in histogram'
                end if
            end if
        end if
    else
        print *,'stationlist: unknown operation code ',oper
        call exit(-1)
    end if
1000 format(a,2f10.4,' ',g10.4,f8.4,1x,4a)
    if ( sign == -1 ) sign = 0
    if ( oper /= 'aver' ) then
        if ( sign >= 0 .and. sign < 1  ) then
            write(2,1000) trim(code),lon,lat,val,sign,trim(name),' (',trim(country),')'
        else
            write(0,*) 'warning: skipping point with sign = ',sign,' ',trim(name)
        end if
    end if

!   next

    goto 200
700 continue
    close(3,status='delete')
    goto 200
800 continue

!   finish up

    if ( coordlist /= ' ' ) then
        open(3,file=coordlist,err=810)
        close(3,status='delete')
    810 continue
    end if
    if ( oper == 'aver' ) then
        call mysystem(aveline,retval)
        call checkval(retval,command)
    end if
    if ( oper(1:2) == 'hi' .or. oper(1:2) == 'at' ) then
        call mysystem('echo "</table>"',retval)
    end if

!   how many points are 'significant'?

    if ( nsigns > 0 .and. (oper == 'higa' .or. oper(1:3) /= 'hig') ) then
        call flush(6)
        call nrsort(nsigns,signs)
        if ( ltwosided ) print '(a)','<b>Positive values</b><br>'
        do isign=1,5
            do i=1,nsigns
                if ( signs(i) > psigns(isign) ) exit
            end do
            if ( i-1 > 0 ) then
                print '(a,i5,a,i5,a,f7.2,a,f6.2,a)','There are ',i-1,'/' &
                    ,nsigns,' (',100*(i-1)/real(n) &
                    ,'%) stations with <i>p</i>&lt;',100*psigns(isign),'%<br>'
                else
                    print '(a,f6.2,a)','No significant trends at <i>p</i>&lt;',100*psigns(isign),'%<br>'
                    exit
            end if
        end do
        if ( ltwosided ) then
            print '(a)','<b>Negative values</b><br>'
            do isign=1,5
                do i=nsigns,1,-1
                    if ( signs(i) < 1-psigns(isign) ) exit
                end do
                i = nsigns - i
                if ( i-1 > 0 ) then
                    print '(a,i5,a,i5,a,f7.2,a,f6.2,a)','There are ',i-1,'/' &
                        ,nsigns,' (',100*(i-1)/real(n) &
                        ,'%) stations with <i>p</i>&lt;',100*psigns(isign),'%<br>'
                else
                    print '(a,f6.2,a)','No significant trends at <i>p</i>&lt;',100*psigns(isign),'%<br>'
                    exit
                end if
            end do
        end if
    end if
    goto 999
900 print *,'stationlist: error: cannot open outfile ',trim(file)
    call exit(-1)
901 print *,'stationlist: error: cannot open file ',trim(file)
    call exit(-1)
902 print *,'stationlist: error: cannot locate ''stations in'' in ',trim(file)
    call exit(-1)
903 print *,'stationlist: error: error reading in ',trim(file)
    call exit(-1)
904 write(0,'(a)') 'stationlist: error reading zdif,sign from '
    write(0,'(a)') trim(line)
999 continue
end program stationlist
