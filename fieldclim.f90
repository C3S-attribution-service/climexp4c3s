program fieldclim
!
!   conmpute field climatology
!
    implicit none
    include 'params.h'
    include 'netcdf.inc'
    include 'getopts.inc'
    integer,parameter :: nvarmax=1
    integer :: nx,ny,nz,nt,nperyear,firstyr,firstmo,lastyr,nvars, &
          ivars(2,nvmax),endian,status,ncid,jvars(6,nvmax)
    integer :: yr,mo,i,j,n,yrbegin,ntmax,ntvarid,mens1,mens
    integer,allocatable :: nn(:,:,:),itimeaxis(:)
    real :: xx(nxmax),yy(nymax),zz(nzmax),undef,lsmask(nxmax,nymax)
    real,allocatable :: field(:,:,:,:),mean(:,:,:),mean2(:,:,:),fxy(:,:), &
        fy(:,:,:)
    character :: file*255,datfile*255,title*255,vars(nvmax)*40, &
        lvars(nvmax)*80,units(nvmax)*20,yesno*1,cell_methods(nvarmax)*100, &
        history*50000,ltime*120,lz(3)*20,svars(nvarmax)*80, &
        metadata(2,100)*2000
    logical :: exist
!
    lwrite = .false.
    if ( command_argument_count() < 2 ) then
        print *,'usage: fieldclim file.[nc|ctl] [begin yr1] [end yr2] [ave n] clim.ctl'
        print *,'computes climatology of field'
        stop
    end if
    call get_command_argument(command_argument_count(),file)
    inquire(file=file,exist=exist)
    if ( exist ) then
        yesno = 'y'
        if ( yesno == 'y' .or. yesno == 'Y' .or. &
             yesno == 'j' .or. yesno == 'J' ) then
            open(1,file=file)
            close(1,status='delete')
            i=index(file,'.ctl')
            datfile=file(1:i)//'.grd'
            open(1,file=datfile)
            close(1,status='delete')
        end if
    end if
    nens1 = 0
    nens2 = 0
    call get_command_argument(1,file)
    call getmetadata(file,mens1,mens,ncid,datfile,nxmax,nx &
        ,xx,nymax,ny,yy,nzmax,nz,zz,lz,nt,nperyear,firstyr,firstmo &
        ,ltime,undef,endian,title,history,nvarmax,nvars,vars,jvars &
        ,lvars,svars,units,cell_methods,metadata,lwrite)
    if ( nz > 1 ) then
        write(0,*) 'fieldclim: error: cannot handle vertical axis yet ',nz
        call exit(-1)
    end if
    yrbegin = firstyr
!   range of years
    lastyr = firstyr + (firstmo+nt-2)/nperyear
!
!   other options
!
    n = command_argument_count()
    call getopts(2,n,nperyear,firstyr,lastyr,.false.,0,0)
    firstyr = max(yr1,firstyr)
    yr1 = firstyr
    lastyr = min(yr2,lastyr)
    yr2 = lastyr
!
!   allocate arrays
!
    allocate(nn(nx,ny,nperyear))
    allocate(mean(nx,ny,nperyear))
    allocate(mean2(nx,ny,nperyear))
    allocate(field(nx,ny,nperyear,firstyr:lastyr))
    allocate(fxy(nperyear,firstyr:lastyr))
    allocate(fy(nperyear,firstyr:lastyr,nx))
!
!   read data
!
    if ( ncid == -1 ) then
        call readdatfile(datfile,field,nx,ny,nx,ny,nperyear,firstyr &
             ,lastyr,yrbegin,firstmo,nt,undef,endian,lwrite,yr1,yr2 &
             ,1,1)
    else
        call readncfile(ncid,field,nx,ny,nx,ny,nperyear,firstyr &
             ,lastyr,yrbegin,firstmo,nt,undef,lwrite,yr1,yr2,jvars)
    end if
!
!   take N-period averages
!
    if ( lsum > 1 ) then
        ! faster
        n = 0
!$omp parallel do private(i,j,yr,mo,fy)
        do j=1,ny
!$omp atomic update
            n = n + 1
            call keepalive1('Summing latitude',n,ny)
            do yr=firstyr,lastyr
                do mo=1,nperyear
                    do i=1,nx
                        fy(mo,yr,i) = field(i,j,mo,yr)
                    end do
                end do
            end do
            do i=1,nx
                call sumit(fy(1,firstyr,i),nperyear,nperyear,firstyr,lastyr,lsum,oper)
            end do
            do yr=firstyr,lastyr
                do mo=1,nperyear
                    do i=1,nx
                        field(i,j,mo,yr) = fy(mo,yr,i)
                    end do
                end do
            end do
        end do
!$omp end parallel do
    end if
!
!   compute climatology
!
    nn = 0
    mean = 0
    do yr=yr1,yr2
        call keepalive1('Processing year',yr-yr1+1,yr2-yr1+1)
        do mo=1,nperyear
            do j=1,ny
                do i=1,nx
                    if ( field(i,j,mo,yr) < 1e33 ) then
                        nn(i,j,mo) = nn(i,j,mo) + 1
                        mean(i,j,mo) = mean(i,j,mo) + field(i,j,mo,yr)
                    end if
                end do
            end do
        end do
    end do
    do mo=1,nperyear
        do j=1,ny
            do i=1,nx
                if ( nn(i,j,mo) > 3 ) then ! arbitrary
                    mean(i,j,mo) = mean(i,j,mo)/nn(i,j,mo)
                else
                    mean(i,j,mo) = 3e33
                end if
            end do
        end do
    end do
!
!   smooth daily climatology with twice a 5-day running mean
!
    if ( nperyear.ge.360 ) then
        call smooth(mean,mean2,nn,nx,ny,nperyear,5)
        call smooth(mean2,mean,nn,nx,ny,nperyear,5)
    end if
!
!   write out
!
    undef = 3e33
    if ( lsum > 1 ) then
        title = 'climatology of running mean of '//title
    else
        title = 'climatology of '//title
    end if
    ivars(1,1) = 0
    ivars(2,1) = 99
    call get_command_argument(command_argument_count(),file)
    i = index(file,'.ctl')
    if ( i == 0 ) then
        ! netcdf output
        ntmax = nperyear
        allocate(itimeaxis(ntmax))
        if ( cell_methods(1) /= ' ' ) then
            cell_methods(1) = trim(cell_methods(1))//', climatology over'
        else
            cell_methods(1) = 'climatology over'
        end if
        write(cell_methods(1),'(2a,i4.4,a,i4.4)') trim(cell_methods(1)),' ',yr1,'-',yr2
        if ( nperyear >= 360 ) cell_methods(1) = trim(cell_methods(1))//', smoothed twice with a 5-day running mean'
        call enswritenc(file,ncid,ntvarid,itimeaxis,ntmax,nx,xx,ny &
            ,yy,nz,zz,lz,nperyear,nperyear,2000,1,ltime,3e33,title &
            ,history,nvars,vars,ivars,lvars,svars,units,cell_methods &
            ,metadata,nens1,nens2)
        do mo=1,nperyear
            call writencslice(ncid,ntvarid,itimeaxis,nt,ivars,mean(1,1,mo),nx,ny,nz,nx,ny,nz,mo,1)
        end do
        status = nf_close(ncid)        
    else
        ! GrADS ctl/dat output - deprecated
        datfile = file
        datfile(i:) = '.grd'
        call writectl(file,datfile,nx,xx,ny,yy,nz,zz, &
            nperyear,nperyear,2000,1,undef,title,nvars,vars,ivars &
            ,lvars,units)
        open(2,file=datfile,form='unformatted',access='direct',recl=4*nx*ny)
        do mo=1,nperyear
            !!!print *,'writing mo ',mo
            write(2,rec=mo) ((mean(i,j,mo),i=1,nx),j=1,ny)
        end do
        close(1)
    end if
end program fieldclim

subroutine smooth(mean,mean2,nn,nx,ny,nperyear,nsmooth)
    implicit none
    integer nx,ny,nperyear,nsmooth
    integer nn(nx,ny,nperyear)
    real mean(nx,ny,nperyear),mean2(nx,ny,nperyear)
    integer mo,i,j,k,mo1
    nn = 0
    mean2 = 0
    do mo=1,nperyear
        do k=-nsmooth/2,nsmooth/2
            mo1 = mo + k
            if ( mo1 < 1 ) mo1 = mo1 + nperyear
            if ( mo1 > nperyear ) mo1 = mo1 - nperyear
            do j=1,ny
                do i=1,nx
                    if ( mean(i,j,mo1) < 1e33 ) then
                        nn(i,j,mo) = nn(i,j,mo) + 1
                        mean2(i,j,mo) = mean2(i,j,mo) + mean(i,j,mo1)
                    end if
                end do
            end do
        end do
    end do
    do mo=1,nperyear
        do j=1,ny
            do i=1,nx
                if ( nn(i,j,mo).ge.2 ) then
                    mean2(i,j,mo) = mean2(i,j,mo)/nn(i,j,mo)
                else
                    mean2(i,j,mo) = 3e33
                end if
            end do
        end do
    end do
end subroutine
