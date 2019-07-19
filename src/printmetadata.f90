subroutine printmetadata(lun,file,FORM_field,title,history,metadata)

    implicit none
    integer :: lun
    character :: file*(*),FORM_field*(*),title*(*),history*(*),metadata(2,100)*(*),scripturl*2000
    integer :: i,j,k,ititle,ifile,ii,iscripturl
    logical :: linstitution
    logical,external :: isnumchar
    
    if ( FORM_field /= ' ' ) then
        write(lun,'(5a)') '# <a href=http://climexp.knmi.nl/select.cgi?field=',trim(FORM_field), &
            '>climexp.knmi.nl/select.cgi?field=',trim(FORM_field),'</a>'
    end if
    ifile = 0
    iscripturl = 0
    if ( file /= ' ' .and. index(file,'aap') == 0 .and. index(file,'noot') == 0 ) then
        ifile = 1
        write(lun,'(2a)') '# file :: ',trim(file)
    end if
    ititle = 0
    if ( title /= ' ' ) then
        i = len_trim(title)
        if ( title(i-2:i) == ' of' ) then
            do j=1,100
                if ( metadata(1,j) == 'title' ) then
                    metadata(2,j) = trim(title)//' '//trim(metadata(2,j))
                    exit
                end if
            end do
        else
            ititle = 1
            write(lun,'(2a)') '# title :: ',trim(title)
        end if
    end if
    linstitution = .false.
    do i=1,100
        if ( metadata(1,i) /= ' ' ) then
            if ( metadata(1,i) == 'conventions' .or. metadata(1,i) == 'Conventions' ) cycle
            if ( metadata(1,i) == 'institution' .or. metadata(1,i) == 'Institution' ) then
                if ( metadata(2,i)(1:10) /= 'KNMI Clima' ) then ! avoid recursion...
                    metadata(2,i) = 'KNMI Climate Explorer and '//trim(metadata(2,i))
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
                do j=1,len(metadata(k,i))
                    if ( metadata(k,i)(j:j) == char(0) ) metadata(k,i)(j:j) = ' '
                    if ( metadata(k,i)(j:j) == '\n' ) metadata(k,i)(j:j) = ' '
                end do
            end do
            write(lun,'(4a)') '# ',trim(metadata(1,i)),' :: ',trim(metadata(2,i))
            if ( metadata(1,i)(1:9) == 'scripturl' .and. &
                 isnumchar(metadata(1,i)(10:10)) .and. isnumchar(metadata(1,i)(11:11)) ) then
                read(metadata(1,i)(10:11),*) ii
                iscripturl = max(iscripturl,ii)
            end if
       end if
    end do
    call getenv('SCRIPTURL',scripturl)
    if ( scripturl /= ' ' ) then
        do i=1,100
            if ( trim(metadata(2,i)) == trim(scripturl) ) exit
        end do
        if ( i > 100 ) then
            iscripturl = iscripturl + 1
            write(lun,'(a,i2.2,2a)') '# scripturl',iscripturl,' :: ',trim(scripturl)
        end if
    end if
    if ( .not.linstitution ) then
        write(lun,'(a)') '# institution :: KNMI Climate Explorer'
    end if
    call extend_history(history)
    write(lun,'(2a)') '# history :: ',trim(history)
    
end subroutine printmetadata

subroutine printvar(unit,var,units,lvar)
    implicit none
    integer :: unit
    character :: var*(*),units*(*),lvar*(*)
    write(unit,'(6a)') '# ',trim(var),' [',trim(units),'] ',trim(lvar)
end subroutine printvar

subroutine extend_history(history)
    implicit none
    character history*(*)
    integer :: i,j,ii(8),l
    character :: string*1000,line*1000
    
    string = ' '
    call getenv('USER',string)
    if ( string /= ' ' ) string = 'user '//string 
    call date_and_time(values=ii)
    write(string(len_trim(string)+2:),'(i4,a,i2.2,a,i2.2)') ii(1),'-',ii(2),'-',ii(3)
    write(string(len_trim(string)+2:),'(i2,a,i2.2,a,i2.2)') ii(5),':',ii(6),':',ii(7)
    do i=0,command_argument_count()
        l = min(len_trim(string) + 2,len(string)-2)
        call get_command_argument(i,line)
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
    integer,intent(out) :: nmetadata1
    character,intent(inout) :: metadata1(2,100)*(*),metadata2(2,100)*(*)
    character,intent(in) :: title2*(*),history2*(*),prefix*(*)
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
!   change occurrences of prefix in metadata1 to another one (hopefully just one level of conflicts...)
    do i=1,nmetadata1
        if ( metadata1(1,i)(1:len_trim(prefix)) == trim(prefix) ) then
            metadata1(1,i) = 'main_'//trim(metadata1(1,i))
        end if
    end do
!   merge
    do k=1,j-1
        if ( metadata2(1,k) /= ' ' ) then
            nmetadata1 = nmetadata1 + 1
            if ( nmetadata1 > 100 ) exit
            metadata1(1,nmetadata1) = trim(prefix)//metadata2(1,k)
            metadata1(2,nmetadata1) = metadata2(2,k)
        end if
    end do
!   also store the history of field2
    if ( history2 /= ' ' .and. nmetadata1 < 100 ) then
        nmetadata1 = nmetadata1 + 1
        metadata1(1,nmetadata1) = trim(prefix)//'history'
        metadata1(2,nmetadata1) = history2
    end if
!   also store the title of field2
    if ( title2 /= ' ' .and. nmetadata1 < 100 ) then
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

subroutine trim_geospatial_metadata(metadata,lon1,lon2,lat1,lat2)

!   adjust the geospatial* metadata to reflect the new box

    implicit none
    real :: lon1,lon2,lat1,lat2
    character :: metadata(2,100)*(*)
    integer :: i

    do i=1,100
        if ( metadata(1,i) == ' ' ) exit
        if ( metadata(1,i) == 'geospatial_lat_min' ) write(metadata(2,i),'(f7.1)') lat1
        if ( metadata(1,i) == 'geospatial_lat_max' ) write(metadata(2,i),'(f7.1)') lat2
        if ( metadata(1,i) == 'geospatial_lon_min' ) write(metadata(2,i),'(f7.1)') lon1
        if ( metadata(1,i) == 'geospatial_lon_max' ) write(metadata(2,i),'(f7.1)') lon2
    enddo
end subroutine trim_geospatial_metadata

subroutine delete_geospatial_metadata(metadata)

!   delete the geospatial* metadata to reflect that a time series has been made from a field.

    implicit none
    character :: metadata(2,100)*(*)
    integer :: i,j

    do i=1,100
        if ( metadata(1,i) == ' ' ) exit
        if ( metadata(1,i)(1:11) == 'geospatial_' ) then
            do j=i+1,100
                metadata(1,j-1) = metadata(1,j)
                metadata(2,j-1) = metadata(2,j)
                if ( metadata(1,j) == ' ' ) exit
            end do
            if ( j > 100 ) then
                metadata(1,j-1) = ' '
                metadata(2,j-1) = ' '
            end if
        end if
    end do
                
end subroutine delete_geospatial_metadata

subroutine coarsen_geospatial_metadata(metadata,avex,avey)

!   adjust the resoloution in the geospatial* metadata to reflect that the field has been coarsened

    implicit none
    character :: metadata(2,100)*(*)
    integer :: avex,avey
    integer :: i
    real :: res

    do i=1,100
        if ( metadata(1,i) == ' ' ) exit
        if ( metadata(1,i) == 'geospatial_lon_resolution' ) then
            read(metadata(2,i),*,err=101,end=101) res
            res = avex*res
            write(metadata(2,i),'(f7.1)') res
        101 continue
        end if
        if ( metadata(1,i) == 'geospatial_lat_resolution' ) then
            read(metadata(2,i),*,err=201,end=201) res
            res = avey*res
            write(metadata(2,i),'(f7.1)') res
        201 continue
        end if
    enddo
                
end subroutine coarsen_geospatial_metadata

subroutine add_metadata(name,value,metadata)

!   add name,value pair to metadata

    implicit none
    character :: name*(*),value*(*),metadata(2,100)*(*)
    integer i
    
    do i=1,100
        if ( metadata(1,i) == ' ' ) exit
    end do
    if ( i > 100 ) return ! not important enough to crash out
    metadata(1,i) = name
    metadata(2,i) = value

end subroutine add_metadata