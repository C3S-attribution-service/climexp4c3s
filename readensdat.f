*  #[ readensdat:
        subroutine readensdat(file,data,npermax,yrbeg,yrend,nensmax
     +        ,nperyear,mens1,mens)
*       
*       read a data file or an ensemble of data files in the array data
*       also returns the number of ensemble members in mens, 0=not an
*       ensemble
*       
        implicit none
        integer npermax,yrbeg,yrend,nperyear,nensmax,mens1,mens,lead
        real data(npermax,yrbeg:yrend,0:nensmax)
        character file*(*)
        integer iens
        logical lexist,lwrite
        character ensfile*255,line*10
        integer llen
        external llen
        lwrite = .false.
        call getenv('READENSDAT_LWRITE',line)
        if ( line(1:1).eq.'t' .or. line(1:1).eq.'T' ) then
            print '(a)','# turned on debug output'
            lwrite = .true.
        endif
        if ( lwrite ) print *,'readensdat: file = ',trim(file)
*       
*       ensembles are denoted by a '%' or '++' in the file name
*
        if ( index(file,'%').eq.0 .and. index(file,'++').eq.0 ) then
            mens = 0
            call makeabsent(data(1,yrbeg,0),npermax,yrbeg,yrend)
            call readdat(data(1,yrbeg,0),npermax,nperyear,yrbeg,yrend
     +            ,file)
            if ( lwrite ) then
                print *,'readensdat: read'
                call printdatfile(6,data(1,yrbeg,0),npermax,nperyear
     +               ,yrbeg,yrend)
            endif
        else
            mens = -1
            mens1 = 0
            do iens=0,nensmax
                ensfile = file
                call filloutens(ensfile,iens)
                if ( lwrite) print *,'looking for file ',trim(ensfile)
                inquire(file=ensfile,exist=lexist)
                if ( .not.lexist ) then
                    if ( mens.eq.-1 ) mens1 = iens + 1
                    cycle
                endif
                if ( lwrite) print *,'reading file '
     +                ,ensfile(1:index(ensfile,' ')-1)
                mens = iens
                call makeabsent(data(1,yrbeg,mens),npermax,yrbeg
     +                ,yrend)
                call readdat(data(1,yrbeg,mens),npermax,nperyear
     +                ,yrbeg,yrend,ensfile)
                if ( mens.eq.mens1 .and. lwrite ) then
                    print *,'readensdat: read'
                    call printdatfile(6,data(1,yrbeg,0),npermax,nperyear
     +                   ,yrbeg,yrend)
                endif
            enddo
            if ( mens.lt.0 ) then
                write(0,*)
     +                'readensdat: error: could not find ensemble '
     +                ,file(1:llen(file)),ensfile(1:llen(ensfile))
                call abort
            endif
        endif
        end
*  #] readensdat:
*  #[ readleadensdat:
        subroutine readleadensdat(file,data,npermax,yrbeg,yrend,nensmax
     +        ,leadmax,nperyear,mens1,mens)
*       
*       read a data file or an ensemble of data files in the array data
*       also returns the number of ensemble members in mens, 0=not an
*       ensemble
*       
        implicit none
        integer npermax,yrbeg,yrend,nperyear,nensmax,leadmax,mens1,mens
     +       ,lead
        real data(npermax,yrbeg:yrend,0:nensmax,leadmax)
        character file*(*)
        integer iens,jens
        logical lexist,lwrite
        character ensfile*255
        integer llen
        external llen
        lwrite = .FALSE.
*       
*       ensembles are denoted by a '%%' or '++' in the file name
*
        if ( index(file,'%%').eq.0 .and. index(file,'++').eq.0 ) then
            mens1 = 0
            mens = 0
            call makeabsent(data(1,yrbeg,0,1),npermax,yrbeg,yrend)
            call readdat(data(1,yrbeg,0,1),npermax,nperyear,yrbeg,yrend
     +            ,file)
        else
            mens = -1
            mens1 = 0
            do lead=1,leadmax
                do iens=0,nensmax
                    ensfile = file
                    call filloutleadens(ensfile,lead,iens)
                    inquire(file=ensfile,exist=lexist)
                    if ( lwrite) print *,'looking for file '
     +                    ,ensfile(1:index(ensfile,' ')-1)
                    if ( .not.lexist ) then
                        if ( mens.eq.-1 ) mens1 = iens + 1
                        cycle
                    endif
                    if ( lwrite) print *,'reading file ',trim(ensfile)
                    mens = max(mens,iens)
                    call makeabsent(data(1,yrbeg,iens,lead),npermax
     +                    ,yrbeg,yrend)
                    call readdat(data(1,yrbeg,iens,lead),npermax
     +                    ,nperyear,yrbeg,yrend,ensfile)
   10               continue
                enddo
            enddo
            if ( mens.lt.0 ) then
                write(0,*)
     +                'readensdat: error: could not find ensemble '
     +                ,file(1:llen(file)),ensfile(1:llen(ensfile))
                call abort
            endif
        endif
        end
*  #] readleadensdat:
