        subroutine savestartstop(yrstart,yrstop)
        implicit none
        integer yrstart,yrstop
        logical lopen
        inquire(unit=12,opened=lopen)
        if ( lopen ) then
            if ( yrstart.ge.-999 ) then
                write(12,'(i4)') yrstart
            else
                write(12,'(i5)') yrstart
            end if
            write(12,'(i4)') yrstop
        end if
        close(12)
        end
