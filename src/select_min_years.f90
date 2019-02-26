program select_min_years
!
!   select stations with a minimum number of years from a stationlist
!
    implicit none
    integer :: minyrs,i,j,n
    character :: file*1024,lines(10)*256

    call get_command_argument(2,file)
    if ( file == ' ' ) then
        write(0,*) 'usage: select_min_years inlist minyrs'
        call exit(-1)
    end if
    read(file,*) minyrs
    call get_command_argument(1,file)
    open(1,file=trim(file),status='old')
    read(1,'(a)') lines(1) ! located...
    i = index(lines(1),'ocated ') + 7
    if ( lines(1)(i:i) /= 's' .and. lines(1)(i:i) /= 'S' ) then
        ! delete number of stations
        j = i
        do while ( lines(1)(j:j) == ' ' )
            j = j + 1
        end do
        do while ( lines(1)(j:j) /= ' ' )
            j = j + 1
        end do
        lines(1)(i:) = lines(1)(j:)
    end if
    write(*,'(a)') trim(lines(1)) 
    do
        do i=1,10
            read(1,'(a)',end=800) lines(i)
            j = index(lines(i),'Found ') + index(lines(i),'found ')
            if ( j > 0 ) then
                read(lines(i)(j+6:),*) n
                if ( n >= minyrs ) then
                    do j=1,i
                        write(*,'(a)') trim(lines(j))
                    end do
                end if
                exit
            end if
        end do
    end do
800 continue
    do j=1,i-1
        write(*,'(a)') trim(lines(j))
    end do
end program select_min_years