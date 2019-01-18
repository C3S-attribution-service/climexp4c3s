integer function get_endian()

!   try to figure out whether I/O is big-endian or little-endian

    implicit none
    integer :: grib=1196575042,birg=1112101447,iu
    integer,save :: endian=0
    integer*4 :: i

    if ( endian == 0 ) then
        call rsunit(iu)
        open(iu,file='/tmp/get_endian',form='unformatted')
        write(iu) 'GRIB'
        rewind(iu)
        read(iu) i
        close(iu,status='delete')
        if ( i == grib ) then
            endian = +1
        elseif ( i == birg ) then
            endian = -1
        endif
    endif
    get_endian = endian
end function get_endian
