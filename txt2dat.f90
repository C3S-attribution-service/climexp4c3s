program txt2dat
!
!   convert the output of a running correlation/regression into a datafile to analyse further
!   or convert an ensemble anomaly file into a set of dat files to analyse further.
!
    implicit none
    integer i,j,iens
    real val,yr
    logical ensemble
    character file*255,line*255,ensfile*255,type*10,wmo*255,headerfile*256
!
    call getarg(1,file)
    if ( index(file,'++').ne.0 .or. index(file,'%%').ne.0 ) then
        ensemble = .true.
        iens = -1
        call getenv('TYPE',type)
    endif
    open(1,file=trim(file),status='old')
    if ( .not.ensemble ) then
        print '(2a)','# ',trim(file)
    else
        i = index(file,'/',.true.) ! rightmost slash
        headerfile = '/tmp/'//file(i+1:)
        open(2,file=trim(headerfile))
    end if
!   loop over all lines in the input file
10  continue
    read(1,'(a)') line
11  continue
    if ( line(1:1).eq.'#' ) then
        if ( .not.ensemble ) then
            ! just copy the comments over
            print '(a)',trim(line)
        else
            i = index(line,'ensemble member')
            if ( i.ne.0 ) then
                ! start a new ensemble member file 
                read(line(i+15:),*) iens
                ensfile = trim(file)//'.dat'
                i = index(ensfile,'.txt.dat')
                if ( i.ne.0 ) ensfile(i:) = '.dat'
                i = index(ensfile,'/')
                if ( ensfile(i+1:i+1).ne.type(1:1) ) then
                    ensfile = ensfile(1:i)//trim(type)//trim(ensfile(i+1:))
                end if
                call filloutens(ensfile,iens)
                close(2)
                i = index(ensfile,'.dat')
                print '(2a)','# ',ensfile(:i-1)
                open(2,file=trim(ensfile),err=901)
                call copyheader(headerfile,2)
            endif
            write(2,'(a)') trim(line)
        endif
        goto 10
    endif
20  continue
    if ( line.eq.' ' ) goto 10
    read(line,*) yr,val
    if ( ensemble ) then
        write(2,*) yr,val
    else
        print *,yr,val
    endif
30  continue
    read(1,'(a)',end=800) line
    if ( line.eq.' ' ) goto 30
    if ( line(1:1).eq.'#' ) then
        if ( ensemble .and. index(line,'ensemble member').ne.0 ) goto 11
        goto 30
    endif
    goto 20
800 continue
    close(1)
    if ( ensemble ) close(2)
    goto 999
901 write(0,*) 'txt2dat: error opening file ',trim(ensfile)
999 continue
end program
