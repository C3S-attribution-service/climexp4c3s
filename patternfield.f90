program patternfield

!   computes the projection of a pattern onto a field, which gives a
!   time series on stdout

    implicit none
    include 'params.h'
    include 'netcdf.inc'
    integer,parameter :: recfa4=4
!	maximum number of variables in the pattern file
    integer,parameter :: nvarmax=80
    integer :: i,j,k,n,ncid1,ncid2,nx1,ny1,nz1,nt1,nper1,firstyr1 &
        ,firstmo1,nx2,ny2,nz2,nt2,nper2,firstyr2,firstmo2,nvars &
        ,ivars1(6,nvarmax),ivars2(6,nvarmax),ivar &
        ,month,endian1,endian2,status,nperyear,nxf,nyf,mens1,mens &
        ,nmetadata1,yrpat1,yrpat2
    integer :: yr1,yr2,nx,ny,nz,yr,mo,iskip
    real ::  xx1(nxmax),yy1(nymax),zz1(nzmax),xx2(nxmax),yy2(nymax),zz2(nzmax),u1,u2
    real :: xx(nxmax),yy(nymax),zz(nymax),minfac
    real,allocatable :: var(:,:),field1(:,:,:,:),field2(:,:,:,:),field2a(:,:,:,:), &
        fxy(:,:),fxynew(:,:)
    logical :: lwrite
    logical :: anom,lstandardunits
    character :: infile*1023,patfile*1023,line*80,datfile1*256,datfile2*256,variable*20,FORM_field*250
    character :: title1*256,vars1(nvarmax)*20,lvars1(nvarmax)*80,units1(nvarmax)*20,svars1(nvarmax)*80
    character :: cell_methods1(nvarmax)*100,history1*50000,ltime1*120,lz1(3)*20,metadata1(2,100)*2000
    character :: title2*256,vars2(nvarmax)*20,lvars2(nvarmax)*80,units2(nvarmax)*20,svars2(nvarmax)*80
    character :: cell_methods2(nvarmax)*100,history2*50000,ltime2*120,lz2(3)*20,metadata2(2,100)*2000

    lwrite = .false. 
    anom = .false. 
    lstandardunits = .false. 

!   check arguments

    n = command_argument_count()
    if ( n < 4 ) then
        print *,'usage: patternfield field.[ctl|nc] '// &
            'pattern.[ctl|nc] variable month [minfac n]'
        call exit(-1)
    endif
    call killfile(title1,datfile1,datfile2,0)
    call get_command_argument(command_argument_count(),line)
    call tolower(line)
    if ( line(1:5) == 'debug' .or. line(1:6) == 'lwrite' ) then
        lwrite = .true. 
        print *,'turned debug output on'
    end if
    call get_command_argument(1,infile)
    call getmetadata(infile,mens1,mens,ncid1,datfile1,nxmax,nx1 &
        ,xx1,nymax,ny1,yy1,nzmax,nz1,zz1,lz1,nt1,nper1,firstyr1,firstmo1 &
        ,ltime1,u1,endian1,title1,history1,1,nvars,vars1,ivars1 &
        ,lvars1,svars1,units1,cell_methods1,metadata1,lwrite)

    call get_command_argument(2,patfile)
    call getmetadata(patfile,mens1,mens,ncid2,datfile2,nxmax,nx2 &
        ,xx2,nymax,ny2,yy2,nzmax,nz2,zz2,lz2,nt2,nper2,firstyr2,firstmo2 &
        ,ltime1,u2,endian2,title2,history2,nvarmax,nvars,vars2,ivars2 &
        ,lvars2,svars2,units2,cell_methods2,metadata2,lwrite)

    if ( nper1 < nper2 ) then
        write(0,*) 'patternfield: error: cannot handle EOF at higher frequency than field ', &
            nper1,nper2
        call exit(-1)
    end if
    nperyear = max(nper1,nper2)

    call get_command_argument(3,variable)
    do ivar=1,nvars
        if ( vars2(ivar) == variable ) then
            goto 100
        endif
    enddo
    write(0,*) 'patternfield: cannot locate ',trim(variable), &
        ' in pattern file ',trim(line)
    write(0,*) '              I only have ',(vars2(ivar),ivar=1,nvars)
    call exit(-1)
100 continue
    if ( ncid2 >= 0 ) then
    !   make sure the variable is the first one in the jvar array
        if ( ivar > 1 ) then
            do i=1,5
                ivars2(i,1) = ivars2(i,ivar)
            enddo
        endif
    endif
    if ( lwrite ) print *,'located ',trim(variable),ivar
    call get_command_argument(4,line)
    read(line,*,err=903) month
    if ( month < 0 .or. month > nperyear ) goto 903
    if ( lwrite ) print *,'picking pattern for month ',month
    do i=len(datfile1),1,-1
        if ( datfile1(i:i) == '/' ) goto 200
    enddo
200 continue
    i = i + 1
    print '(3a,i2,4a)','# patternfield: projecting variable ' &
        ,trim(variable),', month ',month,' of pattern ',trim(title2),' on field ' &
        ,trim(title1)

!   get metadata right...

    call merge_metadata(metadata1,nmetadata1,metadata2,title2,history2,'pattern_')
    nmetadata1 = nmetadata1 + 1
    metadata1(1,nmetadata1) = 'pattern_file'
    metadata1(2,nmetadata1) = patfile
    call getenv('FORM_field',FORM_field)
    call printmetadata(6,infile,FORM_field,title1,history1,metadata1)
    
!   allocate arrays

    nxf = max(nx1,nx2)
    nyf = max(ny1,ny2)
    if ( firstyr2 == 0 .or. firstyr2 == 1 ) then
        yrpat1 = 0
        yrpat2 = 1
    else if ( firstyr2 == 2000 .or. firstyr2 == 2001 ) then
        yrpat1 = 2000
        yrpat2 = 2001
    else
        write(0,*) 'patternfield: warning: unknown convention ',firstyr2
        yrpat1 = firstyr2 - 1
        yrpat2 = yrpat1 + 1
    end if
    allocate(field1(nxf,nyf,nper1,firstyr1:yrend))
    allocate(field2(nxf,nyf,nper2,yrpat1:yrpat2))

    minfac = 0.5
    do i=5,command_argument_count()-1
        if ( iskip > 0 ) then
            iskip = iskip - 1
            cycle
        endif
        call get_command_argument(i,line)
        if ( line(1:6) == 'minfac' ) then
            call get_command_argument(6,line)
            iskip = 1
            read(line,*,err=904) minfac
            if ( minfac > 1 ) minfac=minfac/100
        elseif ( line(1:5) == 'debug' .or. line(1:6) == 'lwrite' ) &
            then
            lwrite = .true. 
            print *,'turned debug output on'
        elseif ( line(1:4) == 'stan' ) then
            lstandardunits = .true. 
            print '(a)','# converting to standard units'
        else
            write(0,*) 'error: do not understand argument ',line
            call abort
        endif
    enddo

!   read fields

    call keepalive(0,3)
    yr1 = firstyr1
    yr2 = min(yrend,firstyr1 + (firstmo1+nt1-1)/nper1)
    write(0,*) 'Reading field...<p>'
    if ( ncid1 == -1 ) then
        call readdatfile(datfile1,field1,nxf,nyf,nx1,ny1,nper1 &
            ,firstyr1,yrend,firstyr1,firstmo1,nt1,u1,endian1 &
            ,lwrite,yr1,yr2,1,1)
    else
        call readncfile(ncid1,field1,nxf,nyf,nx1,ny1,nper2 &
            ,firstyr1,yrend,firstyr1,firstmo1,nt1,u1,lwrite,yr1 &
            ,yr2,ivars1)
    endif
    if ( lstandardunits ) then
    !   convert to standard units
        call makestandardfield(field1,nxf,nyf,1 &
            ,nper1,firstyr1,yrend,nx1,ny1,1,nperyear &
            ,firstyr1,yrend,vars1(1),units1(1),lwrite)
        if ( lwrite ) then
            print *,'patternfield: just after standard units'
            print *,'field1(',(nx1+1)/2,(ny1+1)/2,firstmo1 &
                ,firstyr1,') = ',field1((nx1+1)/2,(ny1+1)/2 &
                ,firstmo1,firstyr1)
        endif
    endif
    call keepalive(1,3)
    write(0,*) 'Reading pattern...<p>'
    if ( ncid2 == -1 ) then
        call readdatfile(datfile2,field2,nxf,nyf,nx2,ny2,nper2,yrpat1 &
            ,yrpat2,firstyr2,firstmo2,nt2,u2,endian2,lwrite,yrpat1,yrpat2,ivar,nvars)
    else
       call readncfile(ncid2,field2,nxf,nyf,nx2,ny2,nper2,yrpat1,yrpat2 &
            ,firstyr2,firstmo2,nt2,u2,lwrite,yrpat1,yrpat2,ivars2)
    endif
    call keepalive(2,3)
    
    if ( nper1 > nper2 ) then
!       extend EOF field to the time scale of field if necessary
        allocate(field2a(nxf,nyf,nperyear,yrpat1:yrpat2))
        allocate(fxy(nper2,0:1))
        allocate(fxynew(nperyear,0:1))
        write(0,*) 'patternfield: error: different time scales not yet ready'
        call exit(-1)
    end if

!   interpolate fields to common grid

    call interpu( &
        field1,xx1,yy1,nx1,ny1, &
        field2,xx2,yy2,nx2,ny2, &
        xx,nx,yy,ny,firstyr1,yr2,yrpat1,yrpat2,nxf,nyf,nperyear,0,lwrite)
    call keepalive(3,3)
    allocate(var(nperyear,yr1:yr2))
    call project(var,nperyear,yr1,yr2,0,0,xx,nx,yy,ny, &
        field1,field2,nxf,nyf,nperyear,firstyr1,yrend, &
        month,minfac,.false.,.false.,anom,lwrite)
    call printdatfile(6,var,nperyear,nperyear,yr1,yr2)
    return

!   error messages

    goto 999
902 write(0,*) 'error: firstyr1,lastyr1 = ',firstyr1,firstyr1 + (nt1-1)/nperyear, &
        ' should be between ',yrbeg,' and ',yrend
    write(0,*) '       recompile if this is too restrictive'
    call exit(-1)
903 write(0,*)'error: month should be between 0 and ',nperyear,', not ',trim(line)
    call exit(-1)
904 write(0,*)'error: cannot read minfac from ',trim(line)
    call exit(-1)
999 continue
end program