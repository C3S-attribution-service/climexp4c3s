        subroutine rkillfile
        implicit none
        character scriptpid*20,email*60,file*85
        integer iu,ls,le
        logical lopen
        integer getpid
*       
        call getenv('SCRIPTPID',scriptpid)
        call getenv('FORM_EMAIL',email)
        do iu=99,10,-1
            inquire(iu,opened=lopen)
            if ( .not.lopen ) goto 20
        enddo
        print '(a)','rsunit: error: no free units under 100!'
        call abort
   20   continue
        do ls=len(scriptpid),1,-1
            if ( scriptpid(ls:ls).ne.' ' ) goto 1
        enddo
 1      continue
        do le=len(email),1,-1
            if ( email(le:le).ne.' ' ) goto 2
        enddo
 2      continue
        write(file,'(4a)') 'pid/',scriptpid(1:ls),'.',
     +       email(1:le)
        open(iu,file=file,status='old',err=10)
        read(iu,'(a)') email
        read(iu,'(a)') email
        write(iu,'(i9)') getpid()
        close(iu)
 10     continue
        end

