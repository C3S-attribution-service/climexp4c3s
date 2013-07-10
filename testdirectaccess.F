        program testdirectaccess
!
!       test whether writing and reading grads binary files works OK
!
        implicit none
#include "recfac.h"
        integer i,j
        real field(100)
        
        print *,'recfac = ',recfac
        print *,'recfa4 = ',recfa4

        do i=1,100
            field(i) = i
        enddo
        open(1,file='test.dat',form='unformatted',access='direct',
     +       recl=recfa4*100)
        write(1,rec=1) field
        close(1)

        field = 0
        open(1,file='test.dat',form='unformatted',access='direct',
     +       recl=recfa4*100)
        read(1,rec=1) field
        print *,'field(1)   = ',field(1)
        print *,'field(100) = ',field(100)

        end

