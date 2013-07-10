        subroutine rankit(aa,ii,n,absent)
*
*       replace the values of aa by their rank in a sorted version
*
        implicit none
        integer n,ii(n)
        real aa(n),absent
        integer i
*
*       work
        call ffsort(aa,ii,n)
        do i=1,n
            if ( aa(ii(i)).lt. 0.9*absent ) then
                print *,aa(ii(i)),i
                aa(ii(i)) = i - n/2.
            endif
        enddo
        end
