        program stationlist
*
*       plough through a statonlist as composed by the get* programs
*       and make a list suitable for feeding to a plotting program:
*
*       # lon1 lon2 lat1 lat2
*       code lon lat value name
*       code lon lat value name
*
*       should do this in python...
*
        implicit none
        integer nmax
        parameter(nmax=250000)
        integer i,j,k,n,ldir,nextarg,retval,nsigns,isign,loptions
        real lon,lon1,lon2,lat,lat1,lat2,val,sign,regr,x1,s1,x2,s2,z,
     +       signs(nmax),psigns(5),mean,sd,sd1,skew,xmin,xmax,chsq,prob,
     +       tx,tx25,tx975,t(10),t25(10),t975(10),alpha,alpha25,alpha975
     +       ,ttx(3),ttx25(3),ttx975(3),tt(10,3),tt25(10,3),tt975(10,3)
        character*20 code
        real elev
        logical lexist,first
        character line*255,file*256,prog*120,filterargs*100,
     +       filterext*100,oper*4,name*40,options*1024,datfile*256
     +       ,rawfile*256,command*1024,tmpfile*29,aveline*10000,dir*256
     +       ,coordlist*30,field*255,fieldfile*255,letter*1,runfile*256
     +       ,histoptions*80,email*80,attribute_args*1024
        logical lwrite,interp
        integer iargc,llen,getpid
        real erfc
        data psigns /0.1,0.05,0.01,0.005,0.001/
*       
        lwrite = .false.
        nsigns = 0
        if ( iargc().lt.4 ) then
            print *,'usage: stationlist listfile outfile prog '//
     +            'operation ...'
            stop
        endif
        call getenv('EMAIL',email)
        call getarg(4,line)
        if ( line(1:6).eq.'ifield' ) then
            oper = line(7:)
            interp = .TRUE.
            letter = 'i'
            call getarg(5,field)
            nextarg = 6
        elseif ( line(1:5).eq.'field' ) then
            oper = line(6:)
            interp = .FALSE.
            letter = 'n'
            call getarg(5,field)
            nextarg = 6
        else
            oper = line
            field = ' '
            nextarg = 5
        endif
        if ( lwrite ) print *,'stationlist: oper = ',oper
        call getarg(3,prog)
        if ( (prog(1:3).ne.'get' .and. prog(1:3).ne.'eca' .and. 
     +       prog(1:4).ne.'beca' .and. prog(1:4).ne.'gdcn' .and.
     +       prog(1:4).ne.'grid' ) .or. 
     +        index(prog,'/').ne.0 .or.
     +        index(prog,';').ne.0 .or. index(prog,'`').ne.0 .or.
     +        index(prog,'&').ne.0 .or. index(prog,'|').ne.0 ) then
            print *,'stationlist: invalid argument ',prog
            call abort
        endif
        if ( lwrite ) print *,'stationlist: prog = ',prog(1:llen(prog))
        i = index(prog,'_')
        if ( i.ne.0 ) then
*           the output has to be filtered
            filterext = prog(i:)
            prog(i:) = ' '
            filterargs = filterext
            do i=1,llen(filterargs)
                if ( filterargs(i:i).eq.'_' ) filterargs(i:i) = ' '
            enddo
        else
            filterext = ' '
        endif
        call getarg(2,file)
        open(2,file=file,status='unknown',err=900)
        call getarg(1,file)
        open(1,file=file,status='old',err=901)
        j=1
        options = ' '
        do i=nextarg,iargc()
            call getarg(i,options(j:))
            if ( options(j:j+3).eq.'runc' .or. options(j:j+3).eq.'runr'
     +           ) then
                call getarg(i+2,runfile)
            endif
            if ( options(j:j+4).eq.'debug' ) then
                lwrite = .true.
            endif
            j = min(len(options)-10,llen(options)+2)
        enddo
        if ( oper(1:2).eq.'hi' ) then
            if ( oper(3:4).eq.'gr' .or. oper(3:4).eq.'gR' ) then
                histoptions = '20 hist sqrtlog '
            else if ( oper(3:4).eq.'pr' .or. oper(3:4).eq.'pR' ) then
                histoptions = '20 hist log '
            else if ( oper(3:4).eq.'vr' .or. oper(3:4).eq.'vR' ) then
                histoptions = '20 hist gumbel '
            else
                histoptions = '20 hist hist '
            endif
            if ( oper(1:3).eq.'hig' ) then
                histoptions(llen(histoptions)+2:) = 'fit gauss'
            elseif ( oper(1:3).eq.'hip' ) then
                histoptions(llen(histoptions)+2:) = 'fit gpd'
            elseif ( oper(1:3).eq.'hiv' ) then
                histoptions(llen(histoptions)+2:) = 'fit gev'
            endif
            if ( lwrite ) print *,'histoptions are ',
     +           histoptions(1:llen(histoptions))
        else if ( oper(1:2).eq.'at' ) then
            call getenv('attribute_args',attribute_args)
        endif
        write(tmpfile,'(a,i20.20)') '/tmp/corr',getpid()
        options(llen(options)+2:) = 'plot '//tmpfile
        loptions = llen(options)
        call getenv('DIR',dir)
        if ( dir.eq.' ' ) dir = '.'
        ldir = llen(dir)
        if ( lwrite ) print *,'DIR = ',dir(1:max(1,ldir))
        if ( ldir.eq.0 ) then
            call getenv('HOME',dir)
            ldir = llen(dir)
            dir = dir(1:llen(dir))//'/climexp/'
            ldir = llen(dir)
            if ( lwrite ) print *,'DIR = ',dir(1:max(1,ldir))
        elseif ( dir(ldir:ldir).ne.'/' ) then
            ldir = ldir + 1
            dir(ldir:ldir) = '/'
            if ( lwrite ) print *,'DIR = ',dir(1:max(1,ldir))
        endif
*
*       write own PID in kill metafile
*
        call killfile(command,line,datfile,0)
*
*       for field-station correlations generate a list of coordinates
*
        if ( field.ne.' ' ) then
            write(coordlist,'(a,i20.20)') '/tmp/coord',getpid()
            open(3,file=coordlist)
            write(3,'(a)') dir(1:ldir)//'data'
   10       continue
            call readcoord(1,lon,lat,elev)
            if ( lon.gt.1e33 ) goto 20
            write(3,'(2f8.2)') lon,lat
            goto 10
   20       continue
            close(3)
            rewind(1)
*
*           and generate the required time series files
*
            command = 'bin/get_index '//field(1:llen(field))//' file '
     +            //coordlist
            if ( interp ) then
                command(llen(command)+1:) = ' interpolate'
            endif
            if ( .false. ) then
                command(llen(command)+1:) = ' debug'
            endif
            if ( lwrite ) print *,'executing ',command(1:llen(command))
            call mysystem(command,retval)
            call checkval(retval,command)
        endif
*
*       search for box
*
  100   continue
        read(1,'(a)',end=902,err=902) line
        i = index(line,'earching for stations') + 
     +       index(line,'ocated stations') + index(line,'regions in') +
     +       index(line,'rid points in') + 
     +       index(line,'Found')*index(line,' stations')
        if ( i.eq.0 ) goto 100
        i = index(line,'ions in') + index(line,'ints in')
        if ( i.ne.0 ) then
*           found box
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
            if ( i.ne.0 ) then
*               found circle
                call readcoord(1,lon1,lat1,elev)
                if ( lon1.gt.1e33 ) goto 180
                lon2 = lon1
                lat2 = lat1
 110            continue
                call readcoord(1,lon,lat,elev)
                if ( lon.gt.1e33 ) goto 180
                if ( abs(lon1-lon).lt.abs(lon1-lon-360) .and.
     +               abs(lon1-lon).lt.abs(lon1-lon+360) ) then
                    lon1 = min(lon1,lon)
                elseif ( abs(lon1-lon-360).lt.abs(lon1-lon+360) ) then
                    lon1 = min(lon1,lon+360)
                else
                    lon1 = min(lon1,lon-360)
                endif
                if ( abs(lon2-lon).lt.abs(lon2-lon-360) .and.
     +               abs(lon2-lon).lt.abs(lon2-lon+360) ) then
                    lon2 = max(lon2,lon)
                elseif ( abs(lon2-lon-360).lt.abs(lon2-lon+360) ) then
                    lon1 = max(lon2,lon+360)
                else
                    lon2 = max(lon2,lon-360)
                endif
                lat1 = min(lat1,lat)
                lat2 = max(lat2,lat)
                goto 110
 180            continue
                
*               extra 10% on all sides
                lon1 = lon1 - (lon2-lon1)/10
                lon2 = lon2 + (lon2-lon1)/11
                lat1 = max(-90.,lat1 - (lat2-lat1)/10)
                lat2 = min(+90.,lat2 + (lat2-lat1)/11)
                rewind(1)
            else
!               no co-ordinates - take the whole world
                lon1 = -30
                lon2 = 330
                lat1 = -90
                lat2 =  90
            end if
        endif
        if ( lon2-lon1.eq.360 ) then
            lon1 = -30
            lon2 = 330
        endif
        write(2,'(a,4f10.4)') '# ',lon1,lon2,lat1,lat2
*
*       get station data
*       
        first = .TRUE.
        if ( oper(1:2).eq.'hi' .or. oper(1:2).eq.'at' ) then
            call mysystem('echo "<table class=realtable width=451 '//
     +           'border=0 cellpadding= cellspacing=0>"',retval)
        endif
        if ( oper.eq.'aver' ) then
            aveline = dir(1:ldir)//'bin/averageseries'
        endif
        n = 0
  200   continue
        if ( lwrite ) print *,'stationlist: attempting to read station '
     +       ,n+1
        call readcoord(1,lon,lat,elev)
        if ( lon.gt.1e33 ) goto 800
        read(1,'(a)',end=800,err=903) line
        call readcodename(line,code,name,lwrite)
*
*       operate
*
        sign = -1
        if ( oper.eq.'plot' ) then
            val = 0.2
        elseif ( oper.eq.'list' ) then
            val = elev
        elseif ( oper.eq.'corr' .or. oper.eq.'sign' .or. oper.eq.'slop'
     +           .or. oper.eq.'nslo' .or. oper.eq.'regr' .or. 
     +           oper.eq.'nreg' .or. oper.eq.'aver' .or.  
     +           oper.eq.'auco' .or. oper.eq.'ausl' .or. oper.eq.'aure' 
     +           .or. oper.eq.'val' .or. oper.eq.'anom' .or. 
     +           oper.eq.'frac' .or. oper.eq.'zval' .or. 
     +           oper.eq.'runc' .or. oper.eq.'zdif' .or. 
     +           oper.eq.'bdif' .or. 
     +           oper.eq.'hime' .or. oper.eq.'hisd' .or. 
     +           oper.eq.'hisk' .or.
     +           oper.eq.'himi' .or. oper.eq.'hima' .or. 
     +           oper(1:3).eq.'hig' .or.
     +           oper(1:3).eq.'hip' .or. oper(1:3).eq.'hiv' .or.
     +           oper(1:2).eq.'at' ) then
            if ( email.eq.'someone@somewhere' .and. n.ge.100 ) then
                print *,'anonymous users can only use 100 stations'
                goto 800
            endif
*           retrieve data
            write(datfile,'(4a)') dir(1:ldir)//'data/',
     +           prog(1:llen(prog)),code(1:llen(code)),'.dat'
            do i=1,llen(datfile)
                if ( datfile(i:i).eq.' ' ) datfile(i:i) = '_'
            enddo
!           the getdutch* and gdcn* script check themselves
!           whether the files need to be regenerated
            if ( prog(1:8).eq.'getdutch' .or. prog(1:4).eq.'gdcn' .or.
     +           prog(1:3).eq.'eca' .or. prog(1:4).eq.'beca' .or.
     +           prog(1:7).eq.'getprcp' .or. prog(1:7).eq.'gettemp' .or.
     +           prog(1:6).eq.'getmin' .or. prog(1:6).eq.'getmax' .or.
     +           prog(1:6).eq.'getslp' .or. prog(1:9).eq.'getfrench' ) 
     +           goto 300
*           check whether the file exists and has length >0
            open(3,file=datfile,err=300)
            read(3,'(a)',err=300,end=300) line
            close(3)
            goto 310
  300       continue
            close(3)
            if ( prog(1:4).eq.'grid' ) then
                write(*,*) 'stationlist: error: the grid point '
     +               ,trim(datfile)
     +               ,' should have been generated already'
                write(0,*) 'stationlist: error: the grid point '
     +               ,trim(datfile)
     +               ,' should have been generated already'
                call abort
            endif
            if ( prog(1:8).eq.'getdutch' .or. prog(1:4).eq.'gdcn' .or.
     +           prog(1:3).eq.'eca' .or. prog(1:4).eq.'beca' .or.
     +       prog(1:7).eq.'getprcp' .or. prog(1:7).eq.'gettemp' .or.
     +       prog(1:6).eq.'getmin' .or. prog(1:6).eq.'getmax' .or.
     +       prog(1:6).eq.'getslp' .or. prog(1:9).eq.'getfrench' ) then
                write(command,'(6a)') trim(dir)//'bin/',trim(prog),
     +                ' ',trim(code),' ',trim(datfile)
            else
                write(command,'(6a)') trim(dir)//'bin/',trim(prog),
     +                ' ',trim(code),' > ',trim(datfile)
            end if
            if ( lwrite ) print *,trim(command)
            call mysystem(command,retval)
            !!!call checkval(retval,command) do not mind a few going wrong...
            if ( retval.ne.0 ) then
                call mysystem('rm -f '//trim(datfile),retval)
                goto 200
            end if
  310       continue
            if ( filterext.ne.' ' ) then
*               filter datfile to longer-period file
                rawfile = datfile
                i = index(rawfile,'.dat')
                datfile = rawfile(:i-1)//filterext(1:llen(filterext))//
     +                '.dat'
                open(3,file=datfile,err=400)
                read(3,'(a)',err=400,end=400) line
                close(3)
                goto 410
  400           continue
                close(3)
                command = dir(1:ldir)//'bin/daily2longer '//
     +                rawfile(1:llen(rawfile))
     +                //filterargs(1:llen(filterargs))//' > '
     +                //datfile
                if ( lwrite ) print *, command(1:llen(command))
                call mysystem(command,retval)
                call checkval(retval,command)
  410           continue
            endif
            if ( oper.eq.'aver' ) then
                i = llen(aveline)
                if ( i.gt.len(aveline)-30 ) then
                    print *,'stationlist: error:'/
     +                    /'too many stations for average'
                    call abort
                endif
                aveline(i+2:) = datfile
            else
*               correlate or plot data
                n = n + 1
                write(line,'(i6)') n
***                call mysystem('echo "'//line(1:7)//name//'"',retval)
                options(loptions+2:) = 'name '//line(1:6)//'_'
     +               //trim(name)
                if ( oper.eq.'val' .or. oper.eq.'frac' .or. 
     +               oper.eq.'anom' .or. oper.eq.'zval' ) then
*                   plot data at a given time
                    write(command,'(6a)') dir(1:ldir)//'bin/getval ',
     +                    datfile(1:llen(datfile)),' ',
     +                    options(1:llen(options))
                    if ( lwrite ) then
                        print *,'stationlist: plot data'
                        print *,'             command = '
     +                       ,command(1:llen(command))
                    endif
                elseif ( field.ne.' ' ) then
                    do i=llen(field),1,-1
                        if ( field(i:i).eq.'/' ) goto 420
                    enddo
  420               continue
                    write(fieldfile,'(3a,f7.2,a,f6.2,3a)')
     +                   'data/grid',field(i+1:llen(field)),'_',lon,'_',
     +                   lat,'_',letter,'.dat'
                    do i=1,llen(fieldfile)
                        if ( fieldfile(i:i).eq.' ' ) fieldfile(i:i)='0'
                    enddo
                    write(command,'(6a)') dir(1:ldir)//'bin/correlate ',
     +                    datfile(1:llen(datfile)),' file ',
     +                    fieldfile(1:llen(fieldfile)),' ',
     +                    options(1:llen(options))
                    if ( lwrite ) then
                        print *,'stationlist: correlate with field'
                        print *,'             command = '
     +                       ,command(1:llen(command))
                    endif
                elseif ( oper(1:2).eq.'au' ) then
                    write(command,'(6a)') dir(1:ldir)//'bin/correlate ',
     +                    datfile(1:llen(datfile)),' file ',
     +                    datfile(1:llen(datfile)),' ',
     +                    options(1:llen(options))
                    if ( lwrite ) then
                        print *,'stationlist: autocorrelation'
                        print *,'             command = '
     +                       ,command(1:llen(command))
                    endif
                elseif ( oper(1:2).eq.'hi' ) then
                    write(command,'(6a)') dir(1:ldir)//'bin/histogram ',
     +                    datfile(1:llen(datfile)),' ',
     +                    histoptions(1:llen(histoptions)),' ',
     +                    options(1:llen(options))
                    if ( lwrite ) then
                        print *,'stationlist: histogram'
                        print *,'             command = '
     +                       ,command(1:llen(command))
                    endif
                elseif ( oper(1:2).eq.'at' ) then
                    write(command,'(6a)') trim(dir)//'bin/attribute ',
     +                    trim(datfile),' ',trim(attribute_args),' ',
     +                   trim(options)
                    if ( lwrite ) then
                        print *,'stationlist: attribute'
                        print *,'             command = '
     +                       ,command(1:llen(command))
                    endif
                else
                    write(command,'(4a)') dir(1:ldir)//'bin/correlate ',
     +                    datfile(1:llen(datfile)),' ',
     +                    options(1:llen(options))
                    if ( lwrite ) then
                        print *,'stationlist: correlate (I think)'
                        print *,'             command = '
     +                       ,command(1:llen(command))
                    endif
                endif
                command(llen(command)+2:)= '|tr -c -d "\\n\\r[:print:]"'
                if ( oper(1:2).eq.'hi' .or. oper(1:2).eq.'at' ) then
                    if ( .not.first ) then
                        if ( oper.eq.'higa' ) then
                            command(llen(command)+2:) = 
     +                           '|egrep  ''(&chi;|th colspan)'''
                        elseif ( oper(1:3).eq.'hig' .or. 
     +                           oper(1:3).eq.'hip' .or.
     +                           oper(1:3).eq.'hiv' ) then
                            command(llen(command)+2:) = 
     +                           '|egrep  ''(return perio|th colspan)'''
                        elseif ( oper.eq.'hisd' ) then
                            command(llen(command)+2:) = '|egrep  '//
     +                           '''(^# <tr><td>s.d..n-1|th colspan)'''
                        elseif ( oper(1:2).eq.'at' ) then ! print the var name in a HTML comment
                            command(llen(command)+2:) = 
     +                           '|egrep ''('//oper//'|th colspan)'''
                        else
                            command(llen(command)+2:) = 
     +                           '|egrep  ''(^# <tr><td>'//oper(3:4)
     +                           //'|th colspan)'''
                        endif
                    else
                        command(llen(command)+2:) = 
     +                       '|egrep  ''(^# <)'''
                    endif
                    command(llen(command)+2:) = '| sed -e "s/# //"'
                else
                    if ( .not.first ) then
                        command(llen(command)+2:) = 
     +                       '|sed -e ''s/# \\([bz]dif\\)/\\1/'''//
     +                       '|egrep -v ''(#|===|sign.  no|% pnts|'//
     +                       'Month.*:|All year|Requiring)'''
                    endif
                endif
                if ( lwrite ) print *,trim(command)
                first = .FALSE.
                call mysystem(command,retval)
*               do not call checkval(retval,command) 
*               as the return value is the one of grep...
*               retrieve results
                if ( oper.eq.'zdif' .or. oper.eq.'bdif' .or. 
     +               oper.eq.'runc' ) then
                    if ( lwrite ) print *,'opening file ',
     +                   runfile(1:llen(runfile))
                    open(3,file=runfile,status='old')
 500                continue
                    read(3,'(a)',end=510) line
                    goto 520
 510                continue
                    call mysystem(
     +                   'echo "error computing running correlations"',
     +                   retval)
***                    call mysystem('echo "'//command//'"',retval)
                    goto 200
 520                continue
                    if ( line(1:6).ne.'# zdif' .and.
     +                   line(1:6).ne.'# bdif' ) goto 500
                    if ( index(line,'****').ne.0 ) then
                       call mysystem('echo "not enough data"',retval)
***                    call mysystem('echo "'//command//'"',retval)
                       goto 200
                    endif
                    read(line,'(9x,f6.2,3x,e14.6)',err=904) val,sign
                    if ( lwrite ) print *,'read val,sign = ',val,sign
                    close(3,status='delete')
                elseif ( oper(1:2).eq.'hi' ) then
                    if ( lwrite ) print *,'sttaionlist: opening file '
     +                   ,tmpfile(1:llen(tmpfile))
                    open(3,file=tmpfile,status='old',err=200)
                    read(3,*,end=200,err=200) mean
                    read(3,*,end=200,err=200) sd
                    read(3,*,end=200,err=200) sd1
                    read(3,*,end=200,err=200) skew
                    read(3,*,end=200,err=200) xmin
                    read(3,*,end=200,err=200) xmax
                    if ( oper.eq.'higa' ) then
                        read(3,*,end=200,err=200) chsq
                        read(3,*,end=200,err=200) prob
                    elseif ( oper(1:3).eq.'hig' .or. 
     +                       oper(1:3).eq.'hip' .or. 
     +                       oper(1:3).eq.'hiv' ) then
                        do i=1,10                            
                            read(3,*,end=200,err=200)t(i),t25(i),t975(i)
                        enddo
                        read(3,*,end=200,err=200) tx,tx25,tx975
                    endif
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
                        enddo
                        print *,'tx   = ',tx,tx25,tx975
                        close(3)
                    else
                        close(3,status='delete')
                    endif
                    if ( oper.eq.'hime' ) then
                        val = mean
                    elseif ( oper.eq.'hisd' ) then
                        val = sd1
                    elseif ( oper.eq.'hisk' ) then
                        val = skew
                    elseif ( oper.eq.'himi' ) then
                        val = xmin
                    elseif ( oper.eq.'hima' ) then
                        val = xmax
                    elseif ( oper.eq.'higa' ) then
                        val = chsq
                        sign = prob
                    elseif ( oper.eq.'higr' .or. oper.eq.'hipr' .or.
     +                       oper.eq.'hivr' ) then
                        val = tx
                        if ( tx.lt.1e33 ) sign = 1/tx
                    elseif ( oper.eq.'higR' .or. oper.eq.'hipR' .or.
     +                       oper.eq.'hivR' ) then
                        val = tx25
                    elseif ( oper(1:3).eq.'hig' .or. oper(1:3).eq.'hip'
     +                       .or. oper(1:3).eq.'hiv'
     +                       ) then
                        read(oper(4:4),'(i1)') i
                        val = t(i+1)
                    else
                        write(*,*) 'error: unknown operation ',oper
                        write(0,*) 'error: unknown operation ',oper
                        call abort
                    endif
                elseif ( oper(1:2).eq.'at' ) then
                    if ( lwrite ) print *,'sttaionlist: opening file '
     +                   ,trim(tmpfile)
                    open(3,file=trim(tmpfile),status='old',err=200)
                    do i=1,10
                        do j=1,2
                            read(3,*,end=200,err=200)tt(i,j),tt25(i,j),
     +                           tt975(i,j)
                        enddo
                    end do
                    do j=1,3
                        read(3,*,end=200,err=200) ttx(j),ttx25(j),
     +                       ttx975(j)
                    end do
                    read(3,*,end=200,err=200) alpha,alpha25,alpha975
                    if ( lwrite ) then
                        print *,'alpha = ',alpha
                        do i=1,10
                            do j=1,2
                                print *,'tt(',i,j,') = ',tt(i,j),
     +                               tt25(i,j),tt975(i,j)
                            end do
                        end do
                        do j=1,3
                            print *,'tx(',j,')   = ',ttx(j),ttx25(j),
     +                           ttx975(j)
                        end do
                        close(3)
                    else
                        close(3,status='delete')
                    endif
                    if ( oper.eq.'atal' ) then
                        val = alpha
                    elseif ( oper.eq.'atr2' ) then
                        val = ttx(1)
                        sign = 1/ttx(1)
                    elseif ( oper.eq.'atr1' ) then
                        val = ttx(2)
                        sign = 1/ttx(2)
                    elseif ( oper.eq.'atra' ) then
                        val = log10(ttx(3))
                        ! assuming the PDF of log(ratio) is gaussian the 95% CI corresponds to ±2\sigma
                        sd = (log(ttx975(3)) - log(ttx25(3)))/4
                        z = abs(log(ttx(3))/sd)
                        ! convert to p-value, assuming a 2-sided test, even though we 
                        ! often know from theory which way it should go.
                        sign = erfc(z/sqrt(2.))
                    elseif ( ichar(oper(3:3)).ge.ichar('0') .and.
     +                       ichar(oper(3:3)).le.ichar('9') ) then
                        read(oper(3:3),'(i1)') i
                        read(oper(4:4),'(i1)') j
                        j = 3-j ! stupid clash of conventions
                        val = tt(i+1,j)
                    else
                        write(*,*) 'error: unknown operation ',oper
                        write(0,*) 'error: unknown operation ',oper
                        call abort
                    endif
                else            ! correlation or getval
                    open(3,file=tmpfile,status='old',err=200)
                    read(3,*,end=700,err=700) i,j,val,sign,k,x1,s1,x2,s2
     +                   ,regr
                    close(3,status='delete')
                endif
                if ( sign.ge.0 ) then
                    if ( nsigns.lt.nmax ) then
                        nsigns = nsigns + 1
                        signs(nsigns) = sign
                    else
                        print *,'last stations not in histogram'
                    endif
                endif
                if ( val.gt.1e33 ) goto 200
                if ( oper.eq.'sign' .or. oper.eq.'ausi' .or. 
     +               oper.eq.'runc' ) then
                    if ( sign.eq.0 ) sign=1e-35
                    val = -log10(sign)
                elseif ( oper.eq.'slop' .or. oper.eq.'ausl' .or.
     +                   oper.eq.'regr' .or. oper.eq.'aure' ) then
                    val = regr
                elseif ( oper.eq.'nslo' .or. oper.eq.'nreg' ) then
                    if ( x1.eq.0 ) then
                        if ( regr.eq.0 ) then
                            val = 0
                        else
                            val = 3e33
                        endif
                    else
                        val = regr/x1
                    endif
                elseif ( oper.eq.'frac' ) then
                    if ( x1.eq.0 ) then
                        if ( val.eq.0 ) then
                            val = 0
                        else
                            val = 3e33
                        endif
                    else
                        val = val/x1-1
                    endif
                    sign = 0
                elseif ( oper.eq.'anom' ) then
                    val = val-x1
                    sign = 0
                elseif ( oper.eq.'zval' ) then
                    val = (val-x1)/s1
                    sign = 0
                endif
            endif
        else
            print *,'stationlist: unknown operation code ',oper
            call abort
        endif
 1000   format(a,2f10.4,' ',g10.4,f8.4,1x,a)
        if ( oper.ne.'aver' ) then
            write(2,1000) code(1:llen(code)),lon,lat,val,sign,trim(name)
        endif
*       
*       next
*       
        goto 200
  700   continue
        close(3,status='delete')
        goto 200
  800   continue
*
*       finish up
*
        if ( coordlist.ne.' ' ) then
            open(3,file=coordlist,err=810)
            close(3,status='delete')
  810       continue
        endif
        if ( oper.eq.'aver' ) then
            call mysystem(aveline,retval)
            call checkval(retval,command)
        endif
        if ( oper(1:2).eq.'hi' .or. oper(1:2).eq.'at' ) then
            call mysystem('echo "</table>"',retval)
        end if
*
*       how many points are 'significant'?
*
        if ( nsigns.gt.0 .and. (oper.eq.'higa' .or. oper(1:3).ne.'hig')
     +       ) then
            call flush(6)
            call nrsort(nsigns,signs)
            do isign=1,5
                do i=1,nsigns
                    if ( signs(i).gt.psigns(isign) ) goto 820
                enddo
 820            print '(a,i5,a,i5,a,f7.2,a,f6.2,a)','There are ',i-1,'/'
     +               ,nsigns,' (',100*(i-1)/real(n)
     +               ,'%) stations with P < ',100*psigns(isign),'%<br>'
            enddo
        endif
        goto 999
  900   print *,'stationlist: error: cannot open outfile ',
     +        file(1:llen(file))
        call abort
  901   print *,'stationlist: error: cannot open file ',
     +        file(1:llen(file))
        call abort
  902   print *,'stationlist: error: cannot locate ''stations in'' in ',
     +        file(1:llen(file))
        call abort
  903   print *,'stationlist: error: error reading in ',
     +        file(1:llen(file))
        call abort
  904   write(0,'(a)') 'stationlist: error reading zdif,sign from '
        write(0,'(a)') line
  999   continue
        end
