        subroutine getword(n,string,delim,word)
!
!       get the Nth word in string delimited by delim
!
        integer n
        character string*(*),delim*1,word*(*)
        integer i,j,m,ii(0:100)

        if ( n.gt.99 ) then
            write(0,*) 'getword: error: can only handle up to 99 words'
            call abort
        endif
        
        if ( n.le.0 ) then
            word = ' '
            return
        end if
        m = 0
        ii(0) = 0
        do i=1,len(string)
            if ( string(i:i).eq.delim ) then
                m = m + 1
                if ( n.gt.100 ) then
                    write(0,*) 'getword: error: can only handle up to '/
     +                   /'99 words'
                    call abort
                endif
                ii(m) = i
                if ( m.ge.n ) exit
            end if
        end do

        if ( n.eq.m+1 ) then
            word = string(ii(m)+1:)
        else if ( n.gt.m ) then
            word = ' '
        else if ( ii(n).gt.ii(n-1)+1 ) then
            word = string(ii(n-1)+1:ii(n)-1)
        else
            word = ' '
        end if
        end
