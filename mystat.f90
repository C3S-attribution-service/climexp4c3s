subroutine mystat(nm,statb,retval)
    implicit none
    integer :: retval
    character*(*) nm
    integer :: statb(13)
    integer :: stat
    retval = stat(trim(nm)//char(0),statb)
end subroutine mystat
