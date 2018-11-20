program lastvalid

!   print out the last valid time step of a series (used is climate diagnostics script)

    implicit none
    include 'param.inc'
    integer :: yr,mo,nperyear
    real :: data(npermax,yrbeg:yrend),undef
    character file*255,var*60,units*60
            
    if ( command_argument_count() /= 1 ) then
        write(0,*) 'usage: lastvalid series.dat'
        call exit(-1)
    end if

    call get_command_argument(1,file)
    call readseries(file,data,npermax,yrbeg,yrend,nperyear,var,units,.false.,.false.)
    do yr=yrend,yrbeg,-1
        do mo=nperyear,1,-1
            if ( data(mo,yr) < 1e33 ) go to 100
        end do
    end do
    write(0,*) 'no valid data in '//trim(file)
    write(*,*) 'no valid data in '//trim(file)
    call exit(-1)
100 continue
    print '(2i5,g15.7,2a)',yr,mo,data(mo,yr),' ',trim(var)
end program
