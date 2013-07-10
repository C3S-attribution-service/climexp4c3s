        program txt2dat
        implicit none
        integer i,j,yr,iens
        real val
        logical ensemble
        character file*255,line*255,ensfile*255,type*10,wmo*255
*
        call getarg(1,file)
        if ( index(file,'++').ne.0 .or. index(file,'%%').ne.0 ) then
            ensemble = .true.
            iens = -1
            call getenv('TYPE',type)
        endif
        open(1,file=file,status='old')
        if ( .not.ensemble ) print '(2a)','# ',trim(file)
        i = 1
 10     continue
        read(1,'(a)') line
 11     continue
        if ( line(1:1).eq.'#' ) then
            if ( ensemble ) then
                i = index(line,'ensemble member')
                if ( i.ne.0 ) then
                    read(line(i+15:),*) iens
                    ensfile = trim(file)//'.dat'
                    i = index(ensfile,'/')
                    ensfile = ensfile(1:i)//'i'//trim(ensfile(i+1:))
                    call filloutens(ensfile,iens)
                    if ( iens.ge.0 ) close(2)
                    print '(2a)','# ',trim(ensfile)
                    open(2,file=ensfile)
                endif
                write(2,'(a)') trim(line)
            else
                print '(a)',trim(line)
            endif
            goto 10
        endif
 20     continue
        read(line,*) yr,val
        if ( ensemble ) then
            write(2,*) yr,val
        else
            print *,yr,val
        endif
 30     continue
        read(1,'(a)',end=800) line
        if ( line.eq.' ' ) goto 30
        if ( line(1:1).eq.'#' ) then
            if ( ensemble .and. index(line,'ensemble member').ne.0 )
     +           goto 11
            goto 30
        endif
        goto 20
 800    continue
        close(1)
        Close(2)
        end
