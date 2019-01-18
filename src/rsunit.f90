subroutine rsunit(irsunit)

!       find a free unit number below 100

    implicit none
    integer :: irsunit
    logical :: lopen
    do irsunit=99,10,-1
        inquire(irsunit,opened=lopen)
        if ( .not. lopen ) exit
    end do
    if ( lopen ) then
        print '(a)','rsunit: error: no free units under 100!'
        call exit(-1)
    end if

end subroutine rsunit
