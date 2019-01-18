subroutine filloutens(file,iens)
    implicit none
    character*(*) :: file
    integer :: iens
    integer :: i,j
    character :: formaat*4,teken*1

    i = index(file,'%')
    if ( i == 0 ) then
        i = index(file,'++')
        if ( i == 0 ) then
            write(0,*) 'filloutens: error: cannot find + or % in file ',trim(file)
            call exit(-1)
        endif
        teken = '+'
    else
        teken = '%'
    endif
    j = i
10  continue
    j = j + 1
    if ( j <= len(file) ) then
        if ( file(j:j) == teken ) go to 10
    endif
    j = j - 1
!**        write(formaat,'(a,i1,a,i1)') 'i',j-i+1,'.',j-i+1
!**        write(file(i:j),formaat) iens
    if ( j == i ) then
        write(file(i:j),'(i1.1)') iens
    elseif ( j == i+1 ) then
        write(file(i:j),'(i2.2)') iens
    elseif ( j == i+2 ) then
        write(file(i:j),'(i3.3)') iens
    elseif ( j == i+3 ) then
        write(file(i:j),'(i4.4)') iens
    else
        write(0,*) 'filloutens: iens too big ',iens,10**(j-i+1)-1
        call exit(-1)
    endif
end subroutine filloutens

subroutine filloutleadens(file,lead,iens)
    implicit none
    character*(*) :: file
    integer :: iens,lead
    integer :: i,j
    character :: formaat*4,teken*1

    i = index(file,'%%')
    if ( i == 0 ) then
        i = index(file,'++')
        if ( i == 0 ) then
            write(0,*) &
            'filloutleadens: error: '// &
            'cannot find ++ or %% in file ',trim(file)
            call exit(-1)
        endif
        teken = '+'
    else
        teken = '%'
    endif
    j = i
10  continue
    j = j + 1
    if ( j <= len(file) ) then
        if ( file(j:j) == teken ) go to 10
    endif
    j = j - 1
!**        write(formaat,'(a,i1,a,i1)') 'i',j-i+1,'.',j-i+1
!**        write(file(i:j),formaat) iens
    if ( iens < 10 .and. j == i ) then
        write(file(i:j),'(i1.1)') iens
    elseif ( iens < 100 .and. j == i+1 ) then
        write(file(i:j),'(i2.2)') iens
    elseif ( iens < 1000 .and. j == i+2 ) then
        write(file(i:j),'(i3.3)') iens
    elseif ( iens < 10000 .and. j == i+3 ) then
        write(file(i:j),'(i4.4)') iens
    else
        write(0,*) 'filloutleadens: iens too large ',iens,10**(j-i)
        call exit(-1)
    endif

    j = i-1
    i = index(file(1:j),'%')
    if ( i == 0 ) then
        i = index(file(1:j),'+_')
        if ( i == 0 ) then
            goto 900
        endif
        teken = '+'
    else
        teken = '%'
    endif
    j = i
20  continue
    j = j + 1
    if ( j <= len(file) ) then
        if ( file(j:j) == teken ) go to 20
    endif
    j = j - 1
!**        write(formaat,'(a,i1,a,i1)') 'i',j-i+1,'.',j-i+1
!**        write(file(i:j),formaat) lead
    if ( lead < 10 .and. j == i ) then
        write(file(i:j),'(i1)') lead
    elseif ( lead < 10 .and. j == i+1 ) then
        write(file(i:j),'(2i1)') 0,lead
    elseif ( lead < 100 .and. j == i+1 ) then
        write(file(i:j),'(i2)') lead
    else
        write(0,*) 'filloutleadens: cannot handle lead>99 yet'
        call exit(-1)
    endif
900 continue
end subroutine filloutleadens
