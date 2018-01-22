subroutine printmetadata(lun,file,FORM_field,title,history,metadata)

    implicit none
    integer :: lun
    character :: file*(*),FORM_field*(*),title*(*),history*(*),metadata(2,100)*(*)
    integer :: i,j,k,ititle,ifile
    logical :: linstitution
    
    if ( FORM_field /= ' ' ) then
        write(lun,'(5a)') '# <a href=http://climexp.knmi.nl/select.cgi?field=',trim(FORM_field), &
            '>climexp.knmi.nl/select.cgi?field=',trim(FORM_field),'</a>'
    end if
    ifile = 0
    if ( file /= ' ' .and. index(file,'aap') == 0 .and. index(file,'noot') == 0 ) then
        ifile = 1
        write(lun,'(2a)') '# file :: ',trim(file)
    end if
    ititle = 0
    if ( title /= ' ' ) then
        ititle = 1
        write(lun,'(2a)') '# title :: ',trim(title)
    end if
    linstitution = .false.
    do i=1,100
        if ( metadata(1,i) /= ' ' ) then
            if ( metadata(1,i) == 'conventions' .or. metadata(1,i) == 'Conventions' ) cycle
            if ( metadata(1,i) == 'institution' .or. metadata(1,i) == 'Institution' ) then
                if ( metadata(2,i)(1:10) /= 'KNMI Clima' ) then ! avoid recursion...
                    metadata(2,i) = 'KNMI Climate Explorer and '//metadata(2,i)
                end if
                linstitution = .true.
            end if
            if ( metadata(1,i) == 'file' .and. ifile > 0 ) then
                ifile = ifile + 1
                metadata(1,i) = 'oldfile'
            end if
            if ( metadata(1,i) == 'oldfile' .and. ifile > 1 ) then
                ifile = ifile + 1
                metadata(1,i) = 'olderfile'
            end if
            if ( metadata(1,i) == 'title' .and. ititle > 0 ) then
                ititle = ititle + 1
                metadata(1,i) = 'oldtitle'
            end if
            if ( metadata(1,i) == 'oldtitle' .and. ititle > 1 ) then
                ititle = ititle + 1
                metadata(1,i) = 'oldertitle'
            end if
            ! gnuplot and others choke on ascii zeros
            ! we also do not want newlines...
            do k=1,2
                do j=1,1000
                    if ( metadata(k,i)(j:j) == char(0) ) metadata(k,i)(j:j) = ' '
                    if ( metadata(k,i)(j:j) == '\n' ) metadata(k,i)(j:j) = ' '
                end do
            end do
            write(lun,'(4a)') '# ',trim(metadata(1,i)),' :: ',trim(metadata(2,i))
        end if
    end do
    if ( .not.linstitution ) then
        write(lun,'(a)') '# institution :: KNMI Climate Explorer'
    end if
    call extend_history(history)
    write(lun,'(2a)') '# history :: ',trim(history)
    
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

subroutine merge_metadata(metadata1,nmetadata1,metadata2,title2,history2,prefix)

!   merge metadata from two datasets:
!   - eliminate common items - often the first field was also used to construct the second one
!   - prefix the rest from dataset 2 with prefix

    implicit none
    integer nmetadata1
    character metadata1(2,100)*(*),metadata2(2,100)*(*),title2*(*),history2*(*),prefix*(*)
    integer i,j,k

!   get number of non-blank items in metadata1
    do i=1,100
        if ( metadata1(1,i) == ' ' ) exit
    end do
    nmetadata1 = i-1
!   set items in metadata2 that also occur in metadata1 to blank
    do j=1,100
        if ( metadata2(1,j) == ' ' ) exit
        do k=1,100
            if ( metadata1(1,k) == ' ' ) exit
            if ( metadata2(1,j) == metadata1(1,k) .and. &
                 metadata2(2,j) == metadata1(2,k) ) then
                metadata2(1,j) = ' '
                metadata2(2,j) = ' '
            end if
        end do
    end do
    do k=1,j-1
        if ( metadata2(1,k) /= ' ' ) then
            nmetadata1 = nmetadata1 + 1
            metadata1(1,nmetadata1) = trim(prefix)//metadata2(1,k)
            metadata1(2,nmetadata1) = metadata2(2,k)
        end if
    end do
!   also store the history of field2
    if ( history2 /= ' ' ) then
        nmetadata1 = nmetadata1 + 1
        metadata1(1,nmetadata1) = trim(prefix)//'history'
        metadata1(2,nmetadata1) = history2
    end if
!   also store the title of field2
    if ( title2 /= ' ' ) then
        nmetadata1 = nmetadata1 + 1
        metadata1(1,nmetadata1) = trim(prefix)//'title'
        metadata1(2,nmetadata1) = title2
    end if
end subroutine merge_metadata

subroutine add_varnames_metadata(var,lvar,svar,metadata,variable)

!   add variable names to global metadata (because they are changed in the calling routine)

    implicit none
    character :: var*(*),lvar*(*),svar*(*),metadata(2,100)*(*),variable*(*)
    integer :: i

    do i=1,100
        if ( metadata(1,i) == ' ' ) exit
    end do
    i = i - 1
    if ( var /= ' ' ) then
        i = i+1
        metadata(1,i) = variable
        metadata(2,i) = var
    end if
    if ( lvar /= ' ' ) then
        i = i+1
        metadata(1,i) = trim(variable)//'_long_name'
        metadata(2,i) = lvar
    end if
    if ( svar /= ' ' ) then
        i = i+1
        metadata(1,i) = trim(variable)//'_standard_name'
        metadata(2,i) = svar
    end if
end subroutine add_varnames_metadata