subroutine toupper(string)
    implicit none
    character*(*) :: string
    integer :: i
    do i=1,len(string)
        if ( ichar(string(i:i)) >= ichar('a') .and. &
             ichar(string(i:i)) <= ichar('z') ) then
            string(i:i) = char(ichar(string(i:i)) - ichar('a') + ichar('A'))
        end if
    end do
end subroutine toupper
