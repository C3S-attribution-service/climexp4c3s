subroutine printmetadata(file,FORM_field,title,history,metadata)

    implicit none
    character :: file*(*),FORM_field*(*),title*(*),history*(*),metadata(2,100)*(*)
    integer :: i,j,k
    logical :: linstitution
    
    if ( FORM_field /= ' ' ) then
        print '(5a)','# <a href=http://climexp.knmi.nl/select.cgi?field=',trim(FORM_field), &
            '>climexp.knmi.nl/select.cgi?field=',trim(FORM_field),'</a>'
    end if
    if ( file /= ' ' ) print '(2a)','# file :: ',trim(file)
    if ( title /= ' ' ) print '(2a)','# title :: ',trim(title)
    linstitution = .false.
    do i=1,100
        if ( metadata(1,i) /= ' ' ) then
            if ( metadata(1,i) == 'conventions' .or. metadata(1,i) == 'Conventions' ) cycle
            if ( metadata(1,i) == 'institution' .or. metadata(1,i) == 'Institution' ) then
                if ( metadata(2,i)(1:10) /= 'KNMI Clima' ) then ! avoid recursion...
                    metadata(2,i) = 'KNMI Climate Explorer and '//metadata(2,i)
                    linstitution = .true.
                end if
            end if
            ! gnuplot and others choke on ascii zeros
            do k=1,2
                do j=1,1000
                    if ( metadata(k,i)(j:j) == char(0) ) metadata(k,i)(j:j) = ' '
                end do
            end do
            print '(4a)','# ',trim(metadata(1,i)),' :: ',trim(metadata(2,i))
        end if
    end do
    if ( .not.linstitution ) then
        print '(a)','# institution :: KNMI Climate Explorer'
    end if
    call extend_history(history)
    print '(2a)','# history :: ',trim(history)
    
end subroutine printmetadata

subroutine extend_history(history)
    implicit none
    character history*(*)
    integer :: i,j,ii(8),l
    character :: string*1000,line*1000
    integer :: iargc
    
    string = ' '
    call getenv('USER',string)
    if ( string /= ' ' ) string = 'user '//string 
    call date_and_time(values=ii)
    write(string(len_trim(string)+2:),'(i4,a,i2.2,a,i2.2)') ii(1),'-',ii(2),'-',ii(3)
    write(string(len_trim(string)+2:),'(i2,a,i2.2,a,i2.2)') ii(5),':',ii(6),':',ii(7)
    do i=0,iargc()
        l = min(len_trim(string) + 2,len(string)-2)
        call getarg(i,line)
        if ( line(1:1) == '/' ) then
            ! convert from absolute path to relative one
            j = index(line,'/climexp/')
            if ( j /= 0 ) then
                line = line(j+9:)
            end if
        end if
        string(l:) = line
        if ( index(string(l:),'startstop') /= 0 ) then
            string(l:) = ' '
        endif
    enddo
    history = trim(string)//'\\n'//trim(history)

end subroutine extend_history
