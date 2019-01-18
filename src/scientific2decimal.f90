program scientific2decimal

!   convert a number from scientific notation to decimal point notation

    implicit none
    real :: x
    integer :: n
    character :: string*20,format*10
    call get_command_argument(1,string)
    if ( string == ' ' ) then
        format='(f20.8)'
    else
        read(string,*) n
        write(format,'(a,i1,a)') '(f10.',n,')'
    end if
    read(*,*,end=999,err=900) x
    print format,x
    goto 999
900 print '(f6.8)',-999.9
999 continue
end program scientific2decimal
