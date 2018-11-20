integer function rindex(s1,s2)

!   to be replaced by the f90 intrinsic one day

    character*(*) :: s1,s2
    integer :: i

    do i=len(s1)-len(s2)+1,1,-1
        if ( s1(i:i+len(s2)-1) == s2 ) go to 100
    enddo
100 continue
    rindex = i

end function rindex
