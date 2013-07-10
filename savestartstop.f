        subroutine savestartstop(yrstart,yrstop)
        implicit none
        integer yrstart,yrstop
        logical lopen
        inquire(unit=12,opened=lopen)
        if ( lopen ) write(12,'(i4)') yrstart,yrstop
        close(12)
        end
