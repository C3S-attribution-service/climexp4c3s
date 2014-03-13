        subroutine mysystem(string,retval)
        implicit none
        integer retval
        character*(*) string

        !!!write(0,*) trim(string)
#ifdef __linux
*       pgf90, but possibly more
        integer system
        retval = system(string)
#else
*       curiously, g95 does not define __linux...
        call system(string,retval)
#endif
        end
