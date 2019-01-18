subroutine mysystem(string,retval)

!   remannt of the Good Old Days when somecompilers demanded a function
!   and others a subroutine. pgf90 requires a function and gfortran accepts both forms.
    implicit none
    integer :: retval
    character :: string*(*)
    integer :: system
    retval = system(string)
end subroutine mysystem
