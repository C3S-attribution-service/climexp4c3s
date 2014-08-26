        integer function getfiletime(file)
        implicit none
        character file*(*)
        integer sarray(13),iret
        integer stat

        iret = stat(trim(file)//char(0),sarray)
        if ( iret /= 0 ) then
            getfiletime = 0
        else
            getfiletime = sarray(10) ! modification time
        end if

        end function
        
