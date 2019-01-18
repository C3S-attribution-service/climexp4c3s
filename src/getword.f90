subroutine getword(n,string,delim,word)

!   get the Nth word in string delimited by delim

    integer,intent(in) :: n
    character,intent(in) :: string*(*),delim*1
    character,intent(out) :: word*(*)
    integer :: i,j,m,ii(0:100)

    if ( n > 99 ) then
        write(0,*) 'getword: error: can only handle up to 99 words'
        call exit(-1)
    endif
            
    if ( n <= 0 ) then
        word = ' '
        return
    end if
    m = 0
    ii(0) = 0
    do i=1,len(string)
        if ( string(i:i) == delim ) then
            m = m + 1
            if ( n > 100 ) then
                write(0,*) 'getword: error: can only handle up to 99 words'
                call exit(-1)
            endif
            ii(m) = i
            if ( m >= n ) exit
        end if
    end do

    if ( n == m+1 ) then
        word = string(ii(m)+1:)
    else if ( n > m ) then
        word = ' '
    else if ( ii(n) > ii(n-1)+1 ) then
        word = string(ii(n-1)+1:ii(n)-1)
    else
        word = ' '
    end if
end subroutine getword
