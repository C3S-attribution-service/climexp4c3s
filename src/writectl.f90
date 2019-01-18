subroutine writectl(file,datfile,nx,xx,ny,yy,nz,zz, &
    nt,nperyear,yrbegin,mobegin,undef,title,nvars,vars,ivars &
    ,lvars,units)

!       write a GrADS .ctl file

    implicit none
!       arguments
    integer :: nx,ny,nz,nt,nperyear,yrbegin,mobegin,nvars, &
    ivars(2,nvars)
    real :: xx(nx),yy(ny),zz(nz),undef
    character :: file*(*),datfile*(*),title*(*),vars(nvars)*(*) &
        ,lvars(nvars)*(*),units(nvars)*(*)
!       local variables
    integer :: i,j,k,n,nfile,unit,dpm(12),yr1
    logical :: lwrite
    character string*1024,months(0:12,2)*3,clwrite*10
!       externals
    integer,external :: lastslash,get_endian
!       date
    data months / &
    'JAN','JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG' &
    ,'SEP','OCT','NOV','DEC','jan','jan','feb','mar','apr' &
    ,'may','jun','jul','aug','sep','oct','nov','dec'/
    data dpm /31,28,31,30,31,30,31,31,30,31,30,31/

    lwrite = .false. 
    call getenv('WRITECTL_LWRITE',clwrite)
    if ( index(clwrite,'T') + index(clwrite,'t') > 0 ) then
        lwrite = .true. 
    endif
    if ( lwrite ) then
        print '(a)','writectl: input'
        print '(2a)','file:     ',trim(file)
        print '(2a)','datfile:  ',trim(datfile)
        print '(a,i4,1000f9.3)','X axis:   ',nx,(xx(i),i=1,nx)
        print '(a,i4,1000f9.3)','Y axis:   ',ny,(yy(i),i=1,ny)
        print '(a,i4,1000f9.1)','Z axis:   ',nz,(zz(i),i=1,nz)
        print '(a,i4,4i5)','T axis:   ',nt,nperyear,yrbegin,mobegin
        print '(a,g12.3)','undef:    ',undef
        print '(2a)','title:    ',trim(title)
        print '(a,i4)','nvars:    ',nvars
        do i=1,nvars
            if ( units(i) /= ' ' .and. ichar(units(i)(1:1)) /= 0 ) &
            then
                print '(a,2i4,1x,4a)',vars(i),ivars(1,i),ivars(2,i) &
                ,trim(lvars(i)),' [',trim(units(i)),']'
            else
                print '(a,2i4,1x,4a)',vars(i),ivars(1,i),ivars(2,i) &
                ,trim(lvars(i))
            endif
        enddo
    endif

!       open file
    call getnfile(file,nfile)
    string = file(1:nfile)//'.ctl'
    if ( lwrite ) print '(2a)','writectl: opening output file ', &
    string(1:nfile+4)
    call rsunit(unit)
    open(unit,file=trim(string),status='new')

!       write file
    if ( datfile(1:1) == '/' .and. &
    datfile(1:lastslash(datfile)) /= file(1:lastslash(file)) ) &
    then
        write(unit,'(2a)') 'DSET ',datfile(1:len_trim(datfile))
    else
        write(unit,'(2a)') 'DSET ^', &
        datfile(lastslash(datfile)+1:len_trim(datfile))
    endif
    write(unit,'(2a)') 'TITLE ',title(1:min(244,len_trim(title)))
    write(unit,'(a,g12.3)') 'UNDEF ',undef
    string = 'OPTIONS'
    if ( get_endian() == +1 ) then
        string(len_trim(string)+2:) = 'BIG_ENDIAN'
    elseif ( get_endian() == -1 ) then
        string(len_trim(string)+2:) = 'LITTLE_ENDIAN'
    endif
    if ( nx > 1 ) then
        if ( xx(2) < xx(1) ) string(len_trim(string)+2:) = 'XREV'
    endif
    if ( ny > 1 ) then
        if ( yy(2) < yy(1) ) string(len_trim(string)+2:) = 'YREV'
    endif
    if ( nz > 1 ) then
        if ( zz(2) < zz(1) ) string(len_trim(string)+2:) = 'ZREV'
    endif
    if ( len_trim(string) > 7 ) &
    write(unit,'(a)') trim(string)
    call writedef(unit,'X',nx,xx)
    call writedef(unit,'Y',ny,yy)
    call writedef(unit,'Z',nz,zz)
    yr1 = max(1,yrbegin)
    if ( nperyear == 12 ) then
        if ( mobegin < 1 .or. mobegin > 12 ) then
            write(0,*) 'writectl: error: month = ',mobegin
            mobegin = 0
        endif
        write(unit,'(a,i6,2a,i4.4,i3,a)') 'TDEF ',nt,' LINEAR 15' &
        ,months(mobegin,1),yr1,12/nperyear,'MO'
    elseif ( nperyear == 4 ) then
        if ( mobegin < 1 .or. mobegin > 4 ) then
            write(0,*) 'writectl: error: season = ',mobegin
            mobegin = 0
        endif
        write(unit,'(a,i6,2a,i4.4,i3,a)') 'TDEF ',nt,' LINEAR 15' &
        ,months(12/nperyear*mobegin-2,1),yr1,12/nperyear, &
        'MO'
    elseif ( nperyear == 2 ) then
        if ( mobegin < 1 .or. mobegin > 2 ) then
            write(0,*) 'writectl: error: season = ',mobegin
            mobegin = 0
        endif
        write(unit,'(a,i6,2a,i4.4,i3,a)') 'TDEF ',nt,' LINEAR 15' &
        ,months(12/nperyear*mobegin-2,1),yr1,12/nperyear, &
        'MO'
    elseif ( nperyear == 1 ) then
        write(unit,'(a,i6,2a,i4.4,a)') 'TDEF ',nt,' LINEAR 1', &
        months(mod(mobegin+5,12)+1,1),yr1+(mobegin+6)/12, &
        ' 1YR'
    elseif ( nperyear < 12 ) then
        write(0,*) 'writectl: cannot handle nperyear = ',nperyear, &
        ' properly yet'
        call exit(-1)
    else
        if ( mobegin < 1 .or. mobegin > nperyear ) then
            write(0,*) 'writectl: error: period = ',mobegin
            i = 0
            j = 0
        else
            n = nint(366.01/nperyear)
            i = nint(n*(mobegin-0.5))
            j = 1
            100 continue
            if ( i > dpm(j) ) then
                i = i - dpm(j)
                j = j + 1
                goto 100
            endif
        endif
        if ( nperyear == 36 ) then
            write(unit,'(a,i6,a,i2,a,i4.4,a)') 'TDEF ',nt &
            ,' LINEAR ',i,months(j,1),yr1,' 14610MN'
        elseif ( nperyear == 48 ) then
            write(unit,'(a,i6,a,i2,a,i4.4,a)') 'TDEF ',nt &
            ,' LINEAR ',i,months(max(j,1),1),yr1,' 10958MN'
        elseif ( nperyear <= 366 ) then
            write(unit,'(a,i6,a,i2,a,i4.4,i4,a)') 'TDEF ',nt &
            ,' LINEAR ',i,months(j,1),yr1,n,'DY'
        elseif ( nperyear <= 24*366 ) then
            n = nint(24*366.001d0/nperyear)
            k = n*mod(mobegin-1,nint(nperyear/366.))
            i = 1 + mobegin/nint(nperyear/366.)
            j = 1
            110 continue
            if ( i > dpm(j) ) then
                i = i - dpm(j)
                j = j + 1
                goto 110
            endif
            write(unit,'(a,i6,a,i2.2,a,i2.2,a,i4.4,i4,a)') 'TDEF ' &
            ,nt,' LINEAR ',k,'Z',i,months(j,1),yr1,n,'HR'
        else
            write(0,*) 'writectl: error: cannot handle nperyear = ' &
            ,nperyear,' yet'
            call exit(-1)
        endif
    endif
    write(unit,'(a,i4)') 'VARS ',nvars
    do i=1,nvars
        if ( units(i) /= ' ' .and. ichar(units(i)(1:1)) /= 0 ) &
        then
            write(unit,'(a,2i6,1x,4a)')vars(i),(ivars(j,i),j=1,2), &
            trim(lvars(i)),' [',trim(units(i)),']'
        else
            write(unit,'(a,2i6,1x,4a)')vars(i),(ivars(j,i),j=1,2), &
            trim(lvars(i))
        endif
    enddo
    write(unit,'(a)') 'ENDVARS'
    close(unit)
    end subroutine writectl

    subroutine writedef(unit,x,nx,xx)

!       write a [XYZ]DEF record, it is assumed the OPTIONS record
!       already contains an [XYZ]REV statement, so the levels are
!       always printed ordered up

    implicit none
    integer :: unit,nx
    real :: xx(nx)
    character(1) :: x
    integer :: i,j,n
    real :: dx

    if ( nx == 1 ) then
        write(unit,'(2a,g12.4,a)') x,'DEF 1 LINEAR ',xx(1),' 1'
        return
    endif
    dx = xx(2) - xx(1)
    do i=2,nx-1
        if ( abs(xx(i+1)-xx(i)-dx) > 1e-5*abs(xx(i)) ) goto 100
    enddo
!       linear
    if ( dx > 0 ) then
        write(unit,'(2a,i6,a,2g14.6)') x,'DEF ',nx,' LINEAR ',xx(1) &
        ,dx
    else
        write(unit,'(2a,i6,a,2g14.6)') x,'DEF ',nx,' LINEAR ',xx(nx) &
        ,-dx
    endif
    return
    100 continue
!       levels
    if ( dx > 0 ) then
        write(unit,'(2a,i6,a,f14.6)') x,'DEF ',nx,' LEVELS ',xx(1)
        do i=1,(nx+4)/6
            write(unit,'(6f14.6)') &
            (xx(j),j=6*i-4,min(6*i+1,nx))
        enddo
    else
        write(unit,'(2a,i6,a,f14.6)') x,'DEF ',nx,' LEVELS ',xx(nx)
        do i=1,(nx+4)/6
            write(unit,'(6f14.6)') &
            (xx(j),j=nx-6*i+5,max(nx-6*i,1),-1)
        enddo
    endif
    return
end subroutine writedef

integer function lastslash(name)
    implicit none
    character*(*) name
    integer :: i
    do i=len(name),1,-1
        if ( name(i:i) == '/' ) goto 100
    enddo
    100 continue
    lastslash = i
end function lastslash

subroutine args2title(title)

!       construct a title from the unix argument list

    implicit none
    character title*(*)
    integer :: i,j,k
    character line*255

    title = ' '
    j = 1
    do i=0,command_argument_count()
        call get_command_argument(i,line)
        if ( index(line,'startstop') /= 0 ) cycle
        k = index(line,'climexp/')
        if ( k > 0 ) then
            k = k + 7
        else
            k = index(line,'LINUX/')
            if ( k > 0 ) k = k + 5
        endif
        title(j:) = line(k+1:)
        j = min(len_trim(title) + 2,len(title))
    enddo
end subroutine args2title

subroutine divideunits(units,units1,units2)

!       make a units string units = units1/units2
!       for the time being very simple

    implicit none
    character units*(*),units1*(*),units2*(*)
    integer :: i
    logical,external ::  multiunit

    if ( units1 == ' ' .or. units2 == ' ' ) then
        units = ' '         ! no idea
    elseif ( units2 == '1' ) then
        units = units1
    elseif ( multiunit(units2) ) then
        units = trim(units1)//'/('//trim(units2)//')'
    else
        units = trim(units1)//'/'//trim(units2)
    endif
end subroutine divideunits

logical function multiunit(unit)

!       checks whether the unit is a product of more than one factor

    implicit none
    character unit*(*)
    integer :: l

    l = len_trim(unit)
    if ( index(unit(:l),'*') /= 0 .or. &
    index(unit(:l),'/') /= 0 .or. &
    index(unit(:l),' ') /= 0 ) then
        multiunit = .true. 
    else
        multiunit = .false. 
    endif
end function multiunit
