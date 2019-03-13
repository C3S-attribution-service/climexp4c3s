subroutine readcoord(unit,lon,lat,elev)

!   search file on unit for a line containing 'oordinates'
!   and read the coordinates from it.  Not very robust yet.

    implicit none
    integer :: unit
    real :: lon,lat,elev
    integer :: i
    character(80) :: line

200 continue
    read(unit,'(a)',end=800,err=800) line
    i = index(line,'oordinates:')
    if ( i == 0 ) goto 200
    call parsecoord(line,lat,lon,elev)
    return
800 continue
    lon = 3e33
    lat = 3e33
    return
end subroutine readcoord

subroutine parsecoord(line,lat,lon,elev)
    implicit none
    real,intent(out) :: lat,lon,elev
    character,intent(in) :: line*(*)
    integer i,j
    i = index(line,'oordinates:')
    i = i+11
    j = index(line,'N')-1
    read(line(i:j),*) lat
    i = index(line,',') + 1
    j = index(line,'E')-1
    read(line(i:j),*) lon
    i = i + index(line(i:),',') + 1
    j = index(line,'m')-1
    if ( j < i ) then
        elev = -999
    else
        read(line(i:j),*) elev
    endif
end subroutine parsecoord

subroutine readcodename(line,code,name,lwrite)

!   read the code and name from line in getstation format:
!   ^[^:]*: code name$

    implicit none
    character line*(*),code*(*),name*(*)
    logical :: lwrite
    integer :: i,j

    i = index(line,':')
210 continue
    i = i+1
    if ( line(i:i) == ' ' ) goto 210
    j = i
220 continue
    j = j+1
    if ( line(j:j) /= ' ' ) goto 220
    code = line(i:j-1)
    if ( lwrite ) print *,'readcodename: found code ',trim(code)
    goto 226
225 continue
    print *, 'readcodename: error reading code from ',line(i:j-1)
    print *, line
    226 continue
    i = j
230 continue
    i = i+1
    if ( line(i:i) == ' ' ) goto 230
    name = line(i:)
    do j=1,len(name)
        if ( name(j:j) == ';' .or. name(j:j) == '`' .or. &
        name(j:j) == '"' .or. name(j:j) == '''' .or. &
        name(j:j) == '&' .or. name(j:j) == '|' .or. &
        name(j:j) == '%' .or. name(j:j) == '<' .or. &
        name(j:j) == '>') name(j:j) = ' '
        if ( name(j:j) == '(' .or. name(j:j) == ')' ) name(j:j)='/'
    enddo
end subroutine readcodename

subroutine checkval(retval,command)
    implicit none
    integer :: retval
    character command*(*)
    if ( retval /= 0 ) then
        write(0,*) 'error ',retval,' in '
        write(0,*) trim(command)
        write(*,*) 'error ',retval,' in '
        write(*,*) trim(command)
        call exit(-1)
    endif
end subroutine checkval

subroutine readcountry(unit,country)

!   search file on unit for a line containing '(country)'
!   and read the country from it.  Not very robust yet.

    implicit none
    integer :: unit
    character country*(*)
    integer :: i,j
    character(80) :: line

200 continue
    read(unit,'(a)',end=800,err=800) line
    i = index(line,'(')
    if ( i == 0 ) goto 200
    i = i+1
    j = index(line,')')-1
    country = line(i:j)
    if ( country(1:4) == 'pop.' ) goto 200
    return
800 continue
    country = 'unknown'
    return
end subroutine readcountry

subroutine readntot(unit,ntot)

!   search file on unit for a line with the total number of stations

    implicit none
    integer :: unit,ntot
    integer :: i,j
    character :: line*100
!
    ntot = 0
    do j=1,50
        read(1,'(a)',end=902,err=902) line
        if ( line(3:4) == '==' ) exit
        i = index(line,'ound ')
        if ( i > 0 ) then
            read(line(i+5:),*,err=190,end=190) ntot
            exit
        end if
        i = index(line,'ocated ')
        if ( i > 0 ) then
            read(line(i+7:),*,err=190,end=190) ntot
            exit
        end if
    end do
190 continue
    return
902 continue
    write(0,*) 'readntot: error reading file, last line was ',trim(line)
    call exit(-1)
end subroutine readntot