program makeundef

!   generate an undefined field

    implicit none
    integer,parameter :: recfa4=4
    integer :: nx,ny,nt,i,j,k
    real :: absent
    character(512) :: string
    integer :: iargc

    if ( iargc() /= 5 ) then
        print *,'usage: makeundef nx ny nt undef outfile'
        stop
    endif
    call getarg(1,string)
    read(string,*) nx
    call getarg(2,string)
    read(string,*) ny
    call getarg(3,string)
    read(string,*) nt
    call getarg(4,string)
    read(string,*) absent
    call getarg(5,string)
    open(1,file=string,access='direct',form='unformatted',recl=recfa4*nx*ny)
    do k=1,nt
        write(1,rec=k) ((absent,i=1,nx),j=1,ny)
    enddo
    close(1)
end program