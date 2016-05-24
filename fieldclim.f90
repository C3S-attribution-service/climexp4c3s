program fieldclim
!
!   conmpute field climatology
!
    implicit none
    include 'params.h'
    include 'netcdf.inc'
    include 'getopts.inc'
    include 'recfac.inc'
    integer nx,ny,nz,nt,nperyear,firstyr,firstmo,lastyr,nvars, &
 &        ivars(2,nvmax),endian,status,ncid,jvars(6,nvmax)
    integer yr,mo,i,j,n,yrbegin
    integer,allocatable :: nn(:,:,:)
    real xx(nxmax),yy(nymax),zz(nzmax),undef,lsmask(nxmax,nymax)
    real,allocatable :: field(:,:,:,:),mean(:,:,:),mean2(:,:,:),fxy(:,:), &
 &      fy(:,:,:)
    character file*255,datfile*255,title*255,vars(nvmax)*40 &
 &        ,lvars(nvmax)*80,units(nvmax)*20,yesno*1
    logical exist
    integer iargc
!
    lwrite = .false.
    if ( iargc().lt.2 ) then
        print *,'usage: fieldclim file.[nc|ctl] '// &
 &           ' [begin yr1] [end yr2] [ave n] clim.ctl'
        print *,'computes climatology of field'
        stop
    endif
    call getarg(iargc(),file)
    inquire(file=file,exist=exist)
    if ( exist ) then
!!!            print *,'output file exists, overwrite?'
!!!            read(*,*) yesno
        yesno = 'y'
        if ( yesno.eq.'y' .or. yesno.eq.'Y' .or. &
 &           yesno.eq.'j' .or. yesno.eq.'J' ) then
            open(1,file=file)
            close(1,status='delete')
            i=index(file,'.ctl')
            datfile=file(1:i)//'.grd'
            open(1,file=datfile)
            close(1,status='delete')
        endif
    endif
    call getarg(1,file)
    if ( lwrite ) print *,'fieldclim: nf_opening file ',trim(file)
    status = nf_open(file,nf_nowrite,ncid)
    if ( status.ne.nf_noerr ) then
        call parsectl(file,datfile,nxmax,nx,xx,nymax,ny,yy,nzmax,nz &
 &            ,zz,nt,nperyear,firstyr,firstmo,undef,endian,title,1 &
 &            ,nvars,vars,ivars,lvars,units)
        ncid = -1
    else
        call parsenc(file,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax &
 &            ,nz,zz,nt,nperyear,firstyr,firstmo,undef,title,1,nvars &
 &            ,vars,jvars,lvars,units)
        datfile = file
    endif
    yrbegin = firstyr
!       range of years
    lastyr = firstyr + (firstmo+nt-2)/nperyear
!
!   other options
!
    n = iargc()
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
    if ( ncid.eq.-1 ) then
        call readdatfile(datfile,field,nx,ny,nx,ny,nperyear,firstyr &
 &           ,lastyr,yrbegin,firstmo,nt,undef,endian,lwrite,yr1,yr2 &
 &           ,1,1)
    else
        call readncfile(ncid,field,nx,ny,nx,ny,nperyear,firstyr &
 &           ,lastyr,yrbegin,firstmo,nt,undef,lwrite,yr1,yr2,jvars)
    endif
!
!   take N-period averages
!
    if ( lsum.gt.1 ) then
        ! faster
        do j=1,ny
            call keepalive1('Summing latitude',j,ny)
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
                    if ( field(i,j,mo,yr).lt.1e33 ) then
                        nn(i,j,mo) = nn(i,j,mo) + 1
                        mean(i,j,mo) = mean(i,j,mo)+field(i,j,mo,yr)
                    endif
                enddo
            enddo
        enddo
    enddo
    do mo=1,nperyear
        do j=1,ny
            do i=1,nx
                if ( nn(i,j,mo).gt.5 ) then ! arbitrary
                    mean(i,j,mo) = mean(i,j,mo)/nn(i,j,mo)
               else
                    mean(i,j,mo) = 3e33
                endif
            enddo
        enddo
    enddo
!
!   smooth daily climatology with twice a 5-day running mean
!
    if ( nperyear.ge.360 ) then
        call smooth(mean,mean2,nn,nx,ny,nperyear,5)
        call smooth(mean2,mean,nn,nx,ny,nperyear,5)
    endif
!
!   write out
!
    call getarg(iargc(),file)
    i=index(file,'.ctl')
    if ( i.eq.0 ) then
        write(0,*) 'fieldclim: error: need .ctl in outfile'
        call abort
    endif
    datfile = file
    datfile(i:) = '.grd'
    undef = 3e33
    if ( lsum.gt.1 ) then
        title = 'climatology of running mean of '//title
    else
        title = 'climatology of '//title
    end if
    ivars(1,1) = 0
    ivars(2,1) = 99
    call writectl(file,datfile,nx,xx,ny,yy,nz,zz, &
 &       nperyear,nperyear,2000,1,undef,title,nvars,vars,ivars &
 &        ,lvars,units)
    open(2,file=datfile,form='unformatted',access='direct',recl=recfa4*nx*ny)
    do mo=1,nperyear
        !!!print *,'writing mo ',mo
        write(2,rec=mo) ((mean(i,j,mo),i=1,nx),j=1,ny)
    enddo
    close(1)
end program

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
            if ( mo1.lt.1 ) mo1 = mo1 + nperyear
            if ( mo1.gt.nperyear ) mo1 = mo1 - nperyear
            do j=1,ny
                do i=1,nx
                    if ( mean(i,j,mo1).lt.1e33 ) then
                        nn(i,j,mo) = nn(i,j,mo) + 1
                        mean2(i,j,mo) = mean2(i,j,mo) + mean(i,j,mo1)
                    endif
                enddo
            enddo
        enddo
    enddo
    do mo=1,nperyear
        do j=1,ny
            do i=1,nx
                if ( nn(i,j,mo).ge.2 ) then
                    mean2(i,j,mo) = mean2(i,j,mo)/nn(i,j,mo)
                else
                    mean2(i,j,mo) = 3e33
                endif
            enddo
        enddo
    enddo
end subroutine
