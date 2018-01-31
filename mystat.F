        subroutine mystat(nm,statb,retval)
        implicit none
        integer retval
        character*(*) nm
        integer statb(13)

        !!!write(0,*) trim(string)
#ifdef __linux
*       pgf90, but possibly more
        integer stat
        retval = stat(nm,statb)
#else
*       curiously, g95 & gfortran do not define __linux...
        call stat(nm,statb,retval)
#endif
        end
        subroutine myctime(itime,string)
        implicit none
        integer itime
        character*(*) string

#ifdef __linux
*       pgf90, but possibly more
        character*30 ctime
        external ctime
        string = ctime(itime)
#else
*       curiously, g95 & gfortran do not define __linux...
        call ctime(itime,string)
#endif
        end
