program dat2grads

!   convert a (set of) time series files into a grads ctl/dat file combo
!   with 1 data point.  Useful to check field versions of programs against
!   point versions :-(

    implicit none
    integer,parameter :: recfa4=4
    include 'param.inc'
    integer :: i,nperyear,yr,mo,yr1,yr2,ix,nx
    real,allocatable :: data(:,:,:)
    real :: lon,lat,alt
    character infile*256,file*256,unit*3,var*40,units*20
    logical :: lwrite
    integer,external :: get_endian

    lwrite = .false. 
    if ( command_argument_count() < 2 ) then
        print *,'usage: dat2grads infile1.dat ... outfile.ctl'
        call exit(-1)
    endif
    nx = command_argument_count()-1
    allocate(data(npermax,yrbeg:yrend,nx))
    do ix=1,nx
        call get_command_argument(ix,infile)
        if ( ix == 1 ) then
            call getcoordinates(infile,lat,lon,alt)
        end if
        call readseries(infile,data(1,yrbeg,ix),npermax,yrbeg,yrend, &
            nperyear,var,units, .false. ,lwrite)
    end do
    yr1 = yrend
    yr2 = yrbeg
    do ix=1,nx
        do yr=yrbeg,yrend
            do mo=1,nperyear
                if ( data(mo,yr,ix) < 1e33 ) then
                    yr1 = min(yr1,yr)
                    yr2 = max(yr2,yr)
                end if
            end do
        end do
    end do
    call get_command_argument(command_argument_count(),file)
    open(2,file=file,status='unknown')
            
    i = index(file,'.ctl')
    if ( i == 0 ) then
        write(0,*) 'dat2grads: error: cannot find ''.ctl'' in ',file
        call exit(-1)
    endif
    file(i:) = '.grd'
    open(3,file=file,status='unknown',form='unformatted', &
    access='direct',recl=nx*recfa4)
    i = 0
    do yr=yr1,yr2
        do mo=1,nperyear
            i = i + 1
            write(3,rec=i) (data(mo,yr,ix),ix=1,nx)
        enddo
    enddo
    close(3)

    write(2,'(2a)') 'DSET ^',trim(file)
    write(2,'(2a)') 'TITLE converted with dat2grads from ',trim(infile)
    write(2,'(a)') 'UNDEF 3e33'
    if ( get_endian() == -1 ) then
        write(2,'(a)') 'OPTIONS LITTLE_ENDIAN'
    else
        write(2,'(a)') 'OPTIONS BIG_ENDIAN'
    endif
    write(2,'(a,i4,a,f8.3,a)') 'XDEF ',nx,' LINEAR ',lon,' 1'
    write(2,'(a,f8.3,a)') 'YDEF 1 LINEAR ',lat,' 1'
    write(2,'(a,f8.1,a)') 'ZDEF 1 LINEAR ',alt,' 1'
    if ( nperyear == 1 ) then
        unit = '1YR'
    elseif ( nperyear == 2 ) then
        unit = '6MO'
    elseif ( nperyear == 4 ) then
        unit = '3MO'
    elseif ( nperyear == 12 ) then
        unit = '1MO'
    elseif ( nperyear == 73 ) then
        unit = '5DY'
    elseif ( nperyear == 366 ) then
        unit = '1DY'
    else
        write(0,*) 'dat2grads: error: cannot handle nperyear = ', &
        nperyear,' yet'
        call exit(-1)
    endif
    write(2,'(a,i6,a,i4.4,2a)') 'TDEF ',nperyear*(yr2-yr1+1), &
    ' LINEAR 1JAN',yr1,' ',unit
    write(2,'(a)') 'VARS 1'
    write(2,'(6a)') trim(var),' 0 99 ',trim(var),' [',trim(units) &
    ,']'
    write(2,'(a)') 'ENDVARS'
    close(2)
    END PROGRAM

    subroutine getcoordinate(line,i,lat)
    implicit none
    integer :: i
    real :: lat
    character line*(*)
    integer :: j
    if ( i /= 0 ) then
        j = i-1
        do while ( line(j:j) == '.' .or. line(j:j) == '-' .or. &
            (ichar(line(j:j)) >= ichar('0') .and. &
            ichar(line(j:j)) <= ichar('9') ) )
            j = j - 1
        end do
        j = j + 1
        if ( j == i ) return
        print *,'getcoordinate: coordinate = ',line(j:i-1)
        read(line(j:i-1),*) lat
        print *,'getcoordinate: coordinate = ',lat
    endif
    end subroutine getcoordinate

    subroutine getcoordinates(infile,lat,lon,alt)
    implicit none
    real :: lat,lon,alt
    character infile*(*)
    integer :: iu,i
    character line*128
    lon = 0
    lat = 0
    alt = 0
    call rsunit(iu)
    open(iu,file=infile,status='old')
    do
    read(iu,'(a)') line
    if ( line(1:1) /= '#' ) exit
    i = index(line,'N, ')
    if ( i /= 0 ) call getcoordinate(line,i,lat)
    i = index(line,'E, ')
    if ( i /= 0 ) call getcoordinate(line,i,lon)
    i = index(line,'m       ')
    if ( i /= 0 ) call getcoordinate(line,i,alt)
end do
    close(iu)
    end subroutine
