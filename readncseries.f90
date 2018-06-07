subroutine readncseries(file,data,npermx,nperyear,yr1,yr2,ncid,var,units,lwrite)
    implicit none
    integer :: npermx,nperyear,yr1,yr2,ncid
    real :: data(npermx,yr1:yr2)
    character :: file*(*),var*(*),units*(*)
    character ::  lvar*120,svar*120,history*50000,metadata(2,100)*1000
    logical :: lwrite

    call readncseriesmeta(file,data,npermx,nperyear,yr1,yr2,ncid,var,units,lvar,svar, &
        history,metadata,lwrite)

end subroutine readncseries

subroutine readncseriesmeta(file,data,npermx,nperyear,yr1,yr2,ncid,var,units,lvar,svar, &
    history,metadata,lwrite)

!   read the data in a 1D netCDF file into data.

    implicit none
    integer :: npermx,nperyear,yr1,yr2,ncid
    real :: data(npermx,yr1:yr2)
    character :: file*(*),var*(*),units*(*)
    character :: lvar*(*),svar*(*),history*(*),metadata(2,100)*(*)
    logical :: lwrite

    include 'param.inc'
    integer,parameter :: mtmax=24*366*200
    integer,parameter :: mxmax=1,mymax=1,mzmax=1,nvmax=1
    integer :: i,j,n,nx,ny,nz,nt,firstyr,firstmo,nvars,ivars(6,nvmax),lastyr
    integer :: nens1,nens2
    real :: xx(mxmax),yy(mymax),zz(mzmax),undef
    real,allocatable :: ddata(:,:,:,:)
    logical :: tdefined(mtmax)
    character :: title*512,vars(nvmax)*20,lvars(nvmax)*80, &
        lz(3)*20,svars(100)*100,ltime*120,cell_methods(100)*100

    ncid = 0
    call ensparsenc(file,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax &
        ,nz,zz,lz,nt,nperyear,firstyr,firstmo,ltime,tdefined,mtmax &
        ,nens1,nens2,undef,title,history,nvmax,nvars,vars,ivars &
        ,lvars,svars,units,cell_methods,metadata)
    if ( lwrite ) print *,'readncseries: nperyear = ',nperyear
    var = vars(1)
    lvar = lvars(1)
    svar = svars(1)
    if ( title /= ' ' ) then
        do n=1,100
            if ( metadata(1,n) == ' ' ) exit
        end do
        if ( n <= 100 ) then
            metadata(1,n) = 'title'
            metadata(2,n) = title
        end if
    end if
    if ( nvars /= 1 ) then
        write(0,*) 'readncseries: error: not just one time-dependent variable, found ',nvars
        call exit(-1)
    end if
    if ( nx > 1 .and. ivars(2,1) /= 0 ) then
        write(0,*) 'readncseries: error: found x-dependent variable, found ',nx
        call exit(-1)
    end if
    if ( ny > 1 .and. ivars(3,1) /= 0 ) then
        write(0,*) 'readncseries: error: found y-dependent variable, found ',ny
        call exit(-1)
    end if
    if ( nz > 1 .and. ivars(4,1) /= 0 ) then
        write(0,*) 'readncseries: error: found z-dependent variable, found ',nz
        call exit(-1)
    end if
    lastyr = firstyr + (nt+firstmo-2)/nperyear
    allocate(ddata(1,1,nperyear,firstyr:lastyr))
    call readncfile(ncid,ddata,1,1,nx,ny,nperyear,firstyr,lastyr, &
        firstyr,firstmo,nt,undef,lwrite, &
        max(yr1,yrbeg,firstyr),min(yr2,yrend,lastyr),ivars)
!!!    ntp = @@@
!!!    call fixholefield(ddata,nx,ny,nz,nperyear,firstyr &
!!!        ,lastyr,firstyr,firstmo,ntp,nt,tdefined,lwrite)
    do i=max(yr1,yrbeg,firstyr),min(yr2,yrend,lastyr)
        do j=1,nperyear
            data(j,i) = ddata(1,1,j,i)
        end do
    end do
    deallocate(ddata)
end subroutine readncseriesmeta

subroutine readncseriesensmeta(file,data,npermx,nperyear,yr1,yr2,nensmx,mens1,mens &
                ,ncid,var,units,lvar,svar,history,metadata,lwrite)

!   read the data in a 2D netCDF file (time,ensemble) into data.

    implicit none
    integer :: npermx,nperyear,yr1,yr2,nensmx,mens1,mens,ncid
    real :: data(npermx,yr1:yr2,0:nensmx)
    character :: file*(*),var*(*),units*(*)
    character :: lvar*(*),svar*(*),history*(*),metadata(2,100)*(*)
    logical :: lwrite

    include 'param.inc'
    integer,parameter :: mtmax=24*366*200
    integer,parameter :: mxmax=1,mymax=1,mzmax=1,nvmax=1
    integer :: i,j,n,nx,ny,nz,nt,firstyr,firstmo,nvars,ivars(6,nvmax),lastyr
    integer :: nens1,nens2,iens
    real :: xx(mxmax),yy(mymax),zz(mzmax),undef
    real,allocatable :: ddata(:,:,:,:,:)
    logical :: tdefined(mtmax)
    character :: title*512,vars(nvmax)*20,lvars(nvmax)*80, &
        lz(3)*20,svars(100)*100,ltime*120,cell_methods(100)*100

    ncid = 0
    if ( lwrite ) print *,'calling ensparsenc'
    call ensparsenc(file,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax &
            ,nz,zz,lz,nt,nperyear,firstyr,firstmo,ltime,tdefined,mtmax &
            ,nens1,nens2,undef,title,history,nvmax,nvars,vars,ivars &
            ,lvars,svars,units,cell_methods,metadata)
    mens1 = 0
    mens = nens2 - nens1
    if ( lwrite ) then
        print *,'readncseriesensmeta: nperyear = ',nperyear
        print *,'                     nens1,nens2 = ',nens1,nens2
    end if
    var = vars(1)
    lvar = lvars(1)
    svar = svars(1)
    if ( title /= ' ' ) then
        do n=1,100
            if ( metadata(1,n) == ' ' ) exit
        end do
        if ( n <= 100 ) then
            metadata(1,n) = 'title'
            metadata(2,n) = title
        end if
    end if
    if ( nvars /= 1 ) then
        write(0,*) 'readncseriesensmeta: error: not just one time-dependent variable, found ',nvars
        call exit(-1)
    end if
    if ( nx > 1 .and. ivars(2,1) /= 0 ) then
        write(0,*) 'readncseriesensmeta: error: found x-dependent variable, found ',nx
        call exit(-1)
    end if
    if ( ny > 1 .and. ivars(3,1) /= 0 ) then
        write(0,*) 'readncseriesensmeta: error: found y-dependent variable, found ',ny
        call exit(-1)
    end if
    if ( nz > 1 .and. ivars(4,1) /= 0 ) then
        write(0,*) 'readncseriesensmeta: error: found z-dependent variable, found ',nz
        call exit(-1)
    end if
    lastyr = firstyr + (nt+firstmo-2)/nperyear
    allocate(ddata(1,1,nperyear,firstyr:lastyr,0:mens))
    if ( lwrite ) print *,'calling readncfileens'
    call readncfileens(ncid,ddata,1,1,1,nx,ny,nz,nperyear,firstyr,lastyr,mens, &
        firstyr,firstmo,nt,mens1,mens,undef,lwrite,max(yr1,yrbeg,firstyr),min(yr2,yrend,lastyr),ivars)
    if ( lwrite ) print *,'copying output to 1D array'
    data = 3e33
    do iens=mens1,mens
        if ( iens > nensmax ) then
            write(0,*) 'readncseriesensmeta: disreagrding ensemble members ',nensmax+1,' to ',mens
            mens = nensmax
            exit
        end if
        do i=max(yr1,yrbeg,firstyr),min(yr2,yrend,lastyr)
            do j=1,nperyear
                data(j,i,iens) = ddata(1,1,j,i,iens)
            end do
        end do
    end do
    deallocate(ddata)
end subroutine readncseriesensmeta