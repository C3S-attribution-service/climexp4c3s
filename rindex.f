*  #[ rindex:
        integer function rindex(s1,s2)
        character*(*) s1,s2
        integer i
        do i=len(s1)-len(s2)+1,1,-1
            if ( s1(i:i+len(s2)-1).eq.s2 ) goto 100
        enddo
  100   continue
        rindex = i
        end
*  #] rindex: