        subroutine swapbyte4a(in,out,n)
c       does a byteswap on a series of integrer*4 numbers to give a real*4 output
        integer n
        integer*4 in(n)
        real*4 out(n)
        
        integer i
        integer*1 ii(4), jj(4)
        integer*4 s
        real*4 t
        equivalence (s,ii)
        equivalence (t,jj)
        
        do i=1,n
            s = in(i)
            jj(1) = ii(4)
            jj(2) = ii(3)
            jj(3) = ii(2)
            jj(4) = ii(1)
            print *,'t,jj = ',t,jj
            out(i) = t
        end do
        
        end
