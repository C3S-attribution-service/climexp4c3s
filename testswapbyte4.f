        program testswapbyte4
        implicit none
        integer field(10,10)
        integer i,j
        do i=1,10
            do j=1,10
                field(i,j) = i + 2**16*j
            enddo
        enddo
        print *,'voor'
        print '(10z8)',field
        call swapbyte4(field,10*10)
        print *,'na'
        print '(10z8)',field
        end
