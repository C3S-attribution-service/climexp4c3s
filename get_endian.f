        integer function get_endian()
*
*       try to figure out whether I/O is big-endian or little-endian
*
        implicit none
        integer endian,grib,birg,iu
        integer*4 i
        save endian
        data endian /0/
        data grib,birg /1196575042,1112101447/

        if ( endian.eq.0 ) then
            call rsunit(iu)
            open(iu,file='/tmp/get_endian',form='unformatted')
            write(iu) 'GRIB'
            rewind(iu)
            read(iu) i
            close(iu,status='delete')
            if ( i.eq.grib ) then
                endian = +1
            elseif ( i.eq.birg ) then
                endian = -1
            endif
        endif
        get_endian = endian
        end
