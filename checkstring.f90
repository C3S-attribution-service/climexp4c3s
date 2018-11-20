subroutine checkstring(string)

!   check string for suspicious characters

    implicit none
    character,intent(inout) :: string*(*)
    character :: allowed*78
    integer :: i,j
    allowed = &
        'abcdefghijklmnopqrstuvwxyz'// &
        'ABCDEFGHIJKLMNOPQRSTUVWXYZ'// &
        '0123456789%^*()-_+ {}[],./'

    do i=1,len(string)
        if ( index(allowed,string(i:i)) == 0 ) then
            string(i:i) = ' '
        end if
    end do

end subroutine checkstring
