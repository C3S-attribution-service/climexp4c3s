integer function getnumwords(string)

!   count the number of space-delimited words in string

    implicit none
    character,intent(in) :: string*(*)
    integer :: i,n

!   hey, this is Fortran!

    i = 0
    n = 0
100 continue
    if ( i >= len(string) ) go to 800
    i = i + 1
!   make sure tabs and windows-style newlines (\r\n) are not accidentally counted as a separate word
    if ( string(i:i) == ' ' .or. string(i:i) == char(9) .or. string(i:i) == char(13) ) go to 100
!   found the beginning of a word
    n = n + 1
200 continue
    if ( i >= len(string) ) go to 800
    i = i + 1
    if ( string(i:i) /= ' ' .and. string(i:i) /= char(9) .and. string(i:i) /= char(13) ) go to 200
    goto 100
800 continue
    getnumwords = n
end function getnumwords
