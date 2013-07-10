        subroutine parsectl(file,datfile,nxmax,nx,xx,nymax,ny,yy,nzmax
     +        ,nz,zz,nt,nperyear,yrbegin,mobegin,undef,endian,title
     +        ,nvarmax,nvars,vars,ivars,lvars,units)
*
*       open GrADS ctl file file, and read out grid info
*
        implicit none
*       arguments
        integer nxmax,nymax,nzmax,nx,ny,nz,nt,nperyear,yrbegin,mobegin
     +        ,endian,nvarmax,nvars,ivars(2,nvarmax)
        real xx(nxmax),yy(nymax),zz(nzmax),undef
        character file*(*),datfile*(*),title*(*),vars(nvarmax)*(*)
     +        ,lvars(nvarmax)*(*),units(nvarmax)*(*)
*       local variables
        integer i,j,k,n,nfile,unit,dpm(12),nperyear0
        real s
        logical xrev,yrev,zrev,lwrite,foundit(0:5)
        character string*1024,months(0:12,3)*3,clwrite*10,
     +       foundvar(0:5)*4
*       externals
        integer llen
        external llen
*       data
        data months /
     +        '???','JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG'
     +        ,'SEP','OCT','NOV','DEC','???','jan','feb','mar','apr'
     +       ,'may','jun','jul','aug','sep','oct','nov','dec','???'
     +       ,'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'
     +       ,'Oct','Nov','Dec'/
        data dpm /31,29,31,30,31,30,31,31,30,31,30,31/
        data foundvar /'dset','xdef','ydef','zdef','tdef','vars'/
        lwrite = .false.
        call getenv('PARSECTL_LWRITE',clwrite)
        if ( index(clwrite,'T') + index(clwrite,'t') .gt.0 ) then
            lwrite = .true.
        endif
        nperyear0 = 366
        foundit = .false.
*
*       open file
        call getnfile(file,nfile)
        string = file(1:nfile)//'.ctl'
        if ( lwrite ) print '(2a)','parsectl: opening input file ',
     +        string(1:nfile+4)
        call rsunit(unit)
        open(unit,file=trim(string),status='old',err=903)
        xrev = .FALSE.
        yrev = .FALSE.
        zrev = .FALSE.
        endian = 0
        nx = 0
        ny = 0
        nz = 0
        nt = 0
        title = ' '
  100   continue
        read(unit,'(a)',end=800,err=900) string
        if ( lwrite ) print '(a)','parsectl: parsing string'
     +        ,string(1:llen(string))
*       comment?
        if ( string.eq.' ' .or. string(1:1).eq.'*' ) goto 100
*       find dataset
        if ( string(1:4).eq.'DSET' .or. string(1:4).eq.'dset' ) then
            j = index(string,'^')
            if ( j.ne.0 ) then
                do i=nfile,1,-1
                    if ( file(i:i).eq.'/' ) goto 110
                enddo
  110           continue
                datfile = file(1:i)//string(j+1:)
            else
                datfile = string(6:)
            endif
            if ( lwrite ) print '(2a)','parsectl: found dset '
     +            ,datfile(1:llen(datfile))
            foundit(0) = .true.
*       title
        elseif ( string(1:5).eq.'TITLE' .or. string(1:5).eq.'title' )
     +            then
            title = string(7:)
            if ( lwrite ) print '(2a)','parsectl: found title '
     +            ,title(1:llen(title))
*       undef
        elseif ( string(1:5).eq.'UNDEF' .or. string(1:5).eq.'undef' )
     +            then
            read(string(7:),*) undef
            if ( lwrite ) print '(a,g10.2)','parsectl: found undef '
     +            ,undef
*       options
        elseif ( string(1:6).eq.'OPTION' .or. string(1:6).eq.'option' )
     +            then
            i = index(string,'XREV') + index(string,'xrev')
            if ( i.ne.0 ) xrev = .TRUE.
            i = index(string,'YREV') + index(string,'yrev')
            if ( i.ne.0 ) yrev = .TRUE.
            i = index(string,'ZREV') + index(string,'zrev')
            if ( i.ne.0 ) zrev = .TRUE.
            if ( lwrite ) print '(a,3l1)','parsectl: found rev options '
     +            ,xrev,yrev,zrev
            i = index(string,'BIG_ENDIAN') + index(string,'big_endian')
            if ( i.ne.0 ) endian = +1
            i = index(string,'LITTLE_ENDIAN') + index(string
     +            ,'little_endian')
            if ( i.ne.0 ) endian = -1
            if ( lwrite ) print '(a,i2)','parsectl: found endian: '
     +           ,endian
            i = index(string,'360_DAY') + index(string,'360_day')
            if ( i.ne.0 ) nperyear0 = 360
            i = index(string,'365_DAY') + index(string,'365_day')
            if ( i.ne.0 ) nperyear0 = 365
            if ( lwrite ) print '(a,i5)','parsectl: found nperyear0: '
     +           ,nperyear0
*       XDEF, YDEF, ZDEF
        elseif ( string(1:4).eq.'XDEF' .or. string(1:4).eq.'xdef' ) then
            call getdef(unit,string,xx,nx,nxmax)
            if ( xrev ) call revit(xx,nx)
            if ( lwrite ) print '(a,i4,1000f7.1)'
     +            ,'parsectl: found X axis ',nx,(xx(i),i=1,nx)
            foundit(1) = .true.
        elseif ( string(1:4).eq.'YDEF' .or. string(1:4).eq.'ydef' ) then
            call getdef(unit,string,yy,ny,nymax)
            if ( yrev ) call revit(yy,ny)
            if ( lwrite ) print '(a,i4,1000f7.1)'
     +            ,'parsectl: found  axis ',ny,(yy(i),i=1,ny)
            foundit(2) = .true.
        elseif ( string(1:4).eq.'ZDEF' .or. string(1:4).eq.'zdef' ) then
            call getdef(unit,string,zz,nz,nzmax)
            if ( zrev ) call revit(zz,nz)
            if ( lwrite ) print '(a,i4,1000f7.1)'
     +            ,'parsectl: found Z axis ',nz,(zz(i),i=1,nz)
            foundit(3) = .true.
*       TDEF
        elseif ( string(1:4).eq.'TDEF' .or. string(1:4).eq.'tdef' ) then
            read(string(6:),*) nt
            if ( nt.le.0 ) then
                write(0,*) 'parsectl: error: nt = ',nt
                write(*,*) 'parsectl: error: nt = ',nt
                call abort
            end if
            i = index(string,'linear') + index(string,'LINEAR')
            if ( i.ne.0 ) then
                i = i + 6
*               skip spaces
  120           continue
                i = i + 1
                if ( string(i:i).eq.' ' ) goto 120
*               skip over initial time string - I want to parse it backwards
  130           continue
                i = i + 1
                if ( string(i:i).ne.' ' ) goto 130
*               parse time step
                j = i
 135            continue
                j = j + 1
                if ( string(j:j).eq.' ' .or.
     +               ichar(string(j:j)).ge.ichar('0') .and.
     +               ichar(string(j:j)).le.ichar('9') ) goto 135
                j = j + 1
                read(string(i+1:j-2),*) s
                if ( lwrite) print *,'read timestep ',s,' from '
     +                ,string(i+1:j-2)
                if (  string(j-1:j).eq.'YR' .or. 
     +                string(j-1:j).eq.'yr' ) then
                    if ( s.eq.1 ) then
                        nperyear = 1
                    else
                        write(0,'(a,f6.0,a)')
     +                       'parsectl: cannot handle ',s
     +                       ,'-yearly reacords'
                        write(*,'(a,f6.0,a)')
     +                       'parsectl: cannot handle ',s
     +                       ,'-yearly reacords'
                        call abort
                    endif
                elseif (  string(j-1:j).eq.'MO' .or. 
     +                string(j-1:j).eq.'mo' ) then
                    nperyear = nint(12/s)
                    if ( abs(nperyear-12/s).gt.0.01 ) then
                        write(*,'(a,f6.2,a)')
     +                        'parsectl: error: cannot handle ',s
     +                        ,' monthly records'
                        write(0,'(a,f6.2,a)')
     +                        'parsectl: error: cannot handle ',s
     +                        ,' monthly records'
                        call abort
                    endif
                elseif (  string(j-1:j).eq.'WK' .or. 
     +                    string(j-1:j).eq.'wk' ) then
                    nperyear = nint(52/s)
                elseif (  string(j-1:j).eq.'DY' .or. 
     +                    string(j-1:j).eq.'dy' ) then
                    nperyear = nint(nperyear0/s)
                    if ( nperyear.eq.37 ) nperyear = 36
                elseif (  string(j-1:j).eq.'HR' .or. 
     +                    string(j-1:j).eq.'hr' ) then
                    nperyear = nint(24*nperyear0/s)
                elseif (  string(j-1:j).eq.'MN' .or. 
     +                    string(j-1:j).eq.'mn' ) then
                    nperyear = nint(60*24*nperyear0/s)
                else
                    write(*,'(a)') 'parsectl: error: cannot parse '
     +                    ,string(j-1:j),' as unit of time'
                    call abort
                endif
*               locate year
                i = i - 4
                if (  ichar(string(i:i)).ge.ichar('0') .and.
     +                ichar(string(i:i)).le.ichar('9') ) then
                    read(string(i:i+3),'(i4)') yrbegin
                else
                    i = i+2
                    read(string(i:i+1),'(i2)') yrbegin
                    if ( yrbegin.lt.30 ) then
                        yrbegin = yrbegin + 2000
                    else
                        yrbegin = yrbegin + 1900
                    endif
                endif
                do mobegin=1,12
                    if (  string(i-3:i-1).eq.months(mobegin,1) .or.
     +                    string(i-3:i-1).eq.months(mobegin,2) .or.
     +                    string(i-3:i-1).eq.months(mobegin,3) ) goto
     +                    140
                enddo
                write(0,*) 'parsectl: error: cannot find month in '
     +                ,string(i-3:i-1)
                mobegin = -1
  140           continue
!               convert month of first data to period of first data
!               it is a bit unfortunate that I use the same symbol for both
                if ( nperyear.eq.1 ) then
!                   no months per year
                    mobegin = 1
                elseif ( nperyear.eq.4 ) then
!                   round of to nearest season
                    mobegin = nint((mobegin+2)/3.)
                elseif ( nperyear.gt.12 ) then
                    read(string(i-5:i-2),'(i2)') n
                    do j=1,mobegin-1
                        n = n + dpm(j)
                    enddo
                    mobegin = 1 + (n-1)*nperyear/365
                endif
                if ( lwrite ) print '(a,i6,1x,i3,1x,i4,i5)'
     +               ,'parsectl: found time axis ',nt
     +               ,mobegin,yrbegin,nperyear
            else
                write(0,*) 'parsectl: error: cannot handle ',string
                call abort
            endif
            foundit(4) = .true.
*       VARS
        elseif ( string(1:4).eq.'VARS' .or. string(1:4).eq.'vars' ) then
            read(string(5:),*) nvars
            if ( lwrite ) print '(a,i3,a)','parsectl: found ',nvars
     +            ,' variables'
            if ( nvars.gt.nvarmax ) goto 902
            do i=1,nvars
                read(unit,'(a)',end=800,err=900) string
                j = 0
*               skip over leading blanks and tabs
  210           continue
                j = j + 1
                if ( j.gt.len(string) ) goto 901
                if ( string(j:j).eq.' ' .or. ichar(string(j:j)).eq.9 )
     +                goto 210
*               find name
                k = j
  220           continue
                k = k + 1
                if ( k.gt.len(string) ) goto 901
                if ( string(k:k).ne.' ' .and. ichar(string(k:k)).ne.9 )
     +                goto 220
                vars(i) = string(j:k-1)
                call checkstring(vars(i))
*               find two numbers
                do n=1,2
                    j = k
  230               continue
                    j = j + 1
                    if ( j.gt.len(string) ) goto 901
                    if ( string(j:j).eq.' ' .or. 
     +                    ichar(string(j:j)).eq.9 ) goto 230
                    k = j
  240               continue
                    k = k + 1
                    if ( k.gt.len(string) ) goto 901
                    if ( string(k:k).ne.' ' .and. 
     +                    ichar(string(k:k)).ne.9 ) goto 240
                    read(string(j:k),*) ivars(n,i)
                enddo
*               rest is long name
                j = k
  250           continue
                j = j + 1
                if ( j.gt.len(string) ) goto 901
                if ( string(j:j).eq.' ' .or. 
     +                ichar(string(j:j)).eq.9 ) goto 250
                lvars(i) = string(j:)
                call checkstring(lvars(i))
                j = index(lvars(i),'[')
                k = index(lvars(i),']')
                if ( j.gt.0 .and. k.gt.j+1 ) then
                    units(i) = lvars(i)(j+1:k-1)
                    lvars(i)(j:k) = ' '
                    call checkstring(units(i))
                else
                    units(i) = ' '
                endif
                if ( lwrite ) print '(2a,2i4,1x,a)'
     +                ,'parsectl: found variable ',vars(i),ivars(1,i)
     +                ,ivars(2,i),lvars(i),' in ',units(i)
            enddo
            read(unit,'(a)',end=800,err=900) string
            if ( index(string,'endvars').eq.0 .and. 
     +            index(string,'ENDVARS').eq.0 ) then
                write(0,*) 'parsectl: warning: no ENDVARS found'
            endif
            if ( nvars.eq.1 ) then
*               consider only the part of the Z axis that is used
                nz = max(1,min(nz,ivars(1,1)))
            endif
            foundit(5) = .true.
        else
            write(0,*) 'parsectl: warning: cannot parse'
            write(0,'(a)') string(1:llen(string))
        endif
        goto 100
  800   continue
        close(unit)
        do i=0,5
            if ( .not.foundit(i) ) then
                write(0,*) 'parsectl: error: did not find ',foundvar(i),
     +               ' in file ',trim(file)
                write(*,*) 'parsectl: error: did not find ',foundvar(i),
     +               ' in file ',trim(file)
                call abort
            endif
        enddo
        if ( lwrite ) then
            print '(a)','parsectl: finshed'
            do i=1,nvars
                print *,i,trim(vars(i)),trim(units(i))
            end do
        end if
        return
  900   write(0,*) 'parsectl: error reading input file'
        write(0,'(a)') string
        call abort
  901   write(0,*) 'parsectl: error reading decsription variable ',i
        write(0,'(a)') string
        call abort
  902   write(0,*) 'parsectl: error: more variables (',nvars
     +        ,') than size of arrays (',nvarmax,')'
        call abort
 903    write(0,*) 'parsectl: error: cannot find input file '
     +       ,trim(string)
        call abort
        end

        subroutine getdef(inunit,string,xx,nx,nxmax)
*       
*       read the [X|Y|Z]DEF statement from unit unit
*       the first line is in string
*       
        implicit none
        integer inunit,nx,nxmax
        character string*256
        real xx(nxmax)
*       
        integer i,n,n0
        real x0,dx
*
        integer llen,getnumnum
        external llen,getnumnum
*
        if ( string(6:6).eq.'*' ) then
            write(0,*) 'getdef: cannot read n from ',trim(string)
            write(*,*) 'getdef: cannot read n from ',trim(string)
            call abort
        endif
        read(string(6:),*) nx
        if ( nx.gt.nxmax ) then
            write(0,*) 'getdef: error: nx>nxmax: ',nx,nxmax
            write(0,*) '        ',string
            call abort
        endif
        i = index(string,'LINEAR') + index(string,'linear')
        if ( i.ne.0 ) then
            read(string(i+7:),*) x0,dx
            do i=1,nx
                xx(i) = x0 + (i-1)*dx
            enddo
        else
            i = index(string,'LEVELS') + index(string,'levels')
            n0 = 0
            if ( i.ne.0 ) then
                i = i + 7
*               loop over lines
  100           continue
                n = getnumnum(string(i:))
                read(string(i:),*) (xx(i),i=n0+1,min(n0+n,nx))
                n0 = n0+n
                if ( n0.lt.nx ) then
                    read(inunit,'(a)') string
                    i = 1
                    goto 100
                endif
            else
                write(0,*) 'error: cannot parse ',string(1:llen(string))
                call abort
            endif
        endif
        end

        subroutine getnfile(file,nfile)
        implicit none
        integer nfile
        character*(*) file
        nfile = index(file,'.ctl')
        if ( nfile.eq.0 ) then
            nfile = index(file,'.dat')
            if ( nfile.eq.0 ) then
                nfile = index(file,' ')
                if ( nfile.eq.0 ) then
                    write(0,*) 'error: filename too long ',file
                    call abort
                endif
            endif
        endif
        nfile = nfile - 1
        end

        integer function getnumnum(string)
*       
*       return the number of numbers (i.e., blank-delimited entities)
*       in string
*       
        implicit none
        character*(*) string
        integer i,n
        n = 0
        i = 0
  100   continue
        i = i + 1
        if ( i.gt.len(string) ) then
            getnumnum = n
            return
        endif
        if ( string(i:i).eq.' ' ) then
            goto 100
        endif
        n = n + 1
  200   continue
        i = i + 1
        if ( i.gt.len(string) ) then
            getnumnum = n
            return
        endif
        if ( string(i:i).eq.' ' ) then
            goto 100
        else
            goto 200
        endif
        end
*
        subroutine revit(x,n)
        implicit none
        integer i,n
        real y,x(n)
        do i=1,n/2
            y = x(i)
            x(i) = x(n-i+1)
            x(n-i+1) = y
        enddo
        end

        subroutine datfile2modelname(datfile,modelname)
!
!       find a unique model name from the datfile name
!       only works for Demeter & eurosip data for now - but that is the
!       only place we have multi-model ensembles at the moment.
!
        character datfile*(*),modelname*(*)
        integer k
        integer,external :: rindex
        k = rindex(datfile,'/')
        modelname = datfile(k+1:)
        k = index(modelname,'_')
        if ( k.gt.0 ) modelname(k:) = ' '
        end
