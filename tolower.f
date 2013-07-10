*  #[ tolower:
        subroutine tolower(string)
        implicit none
        character*(*) string
        integer i
        do i=1,len(string)
            if ( ichar(string(i:i)).ge.ichar('A') .and. 
     +           ichar(string(i:i)).le.ichar('Z') ) then
                string(i:i) = char(ichar(string(i:i)) - ichar('A') +
     +              ichar('a'))
            endif
        enddo
*  #] tolower:
        end
