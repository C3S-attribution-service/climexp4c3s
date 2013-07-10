        subroutine killfile(scriptpid,email,file,who)
*
*       write the PID into a kill metafile in the Climate Explorer
*       arguments: three scratch character strings
*       
        implicit none
        character scriptpid*(*),email*(*),file*(*)
        integer iu,who
        integer getpid,getppid,llen
*       
        call getenv('SCRIPTPID',scriptpid)
        if ( scriptpid.eq.' ' ) return
        call getenv('EMAIL',email)
        if ( email.eq.' ' ) then
            call getenv('FORM_EMAIL',email)
            if ( email.eq.' ' ) return
        endif
        call rsunit(iu)
        write(file,'(4a)') 'pid/',scriptpid(1:llen(scriptpid)),'.',
     +       email(1:llen(email))
        open(iu,file=file,status='old',err=10)
        read(iu,'(a)') email
        read(iu,'(a)') email
        if ( who.eq.0 ) then
            write(iu,'(i9)') getpid()
        else
            write(iu,'(i9)') getppid()
        endif
        close(iu)
 10     continue
        end
