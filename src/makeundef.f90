program makeundef

!   generate an undefined field

    implicit none
    integer,parameter :: recfa4=4
    integer :: nx,ny,nt,i,j,k
    real :: absent
    character(512) :: string

    if ( command_argument_count() /= 5 ) then
        print *,'usage: makeundef nx ny nt undef outfile'
        stop
    endif
    call get_command_argument(1,string)
    read(string,*) nx
    call get_command_argument(2,string)
    read(string,*) ny
    call get_command_argument(3,string)
    read(string,*) nt
    call get_command_argument(4,string)
    read(string,*) absent
    call get_command_argument(5,string)
    open(1,file=string,access='direct',form='unformatted',recl=recfa4*nx*ny)
    do k=1,nt
        write(1,rec=k) ((absent,i=1,nx),j=1,ny)
    enddo
    close(1)
end program