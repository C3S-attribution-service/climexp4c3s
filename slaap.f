        program slaap
*
*       test hoe ik asynchrone keepalive boodschappen kan sturen
*
        implicit none
        integer i,n,getpid
        character*10 string,string1

        call getarg(1,string)
        if ( string.eq.' ' ) then
            print *,'usage: slaap n'
            stop
        endif
        read(string,*) n

        write(string1,'(i8)') getpid()
        call putenv('EOFID='//string1)
        call system('./stillsleeping.cgi '//string//' &')
        do i=1,n
            print *,'hard aan het werk ',i,'<p>'
            call sleep(1)
        enddo
        
        end
