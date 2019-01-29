program averagefieldspace

!   average a few grid boxes together
!   should do this in ferret

    implicit none
    include 'param.inc'
    include 'netcdf.inc'
    integer,parameter :: recfa4=4
    integer :: i,j,ii,jj,iii,jjj,it,n,mo,yr,jx,jy,avex,avey,status,nxf,nyf,nzf
    integer :: mens,mens1,ncid,nx,ny,nz,nt,nperyear,firstyr,firstmo,endian,nvars &
        ,ivars(6,1),lastyr,ntvarid,itimeaxis(ntmax),nnt
    real :: xx(nxmax),yy(nymax),zz(1),undef,minfac,s,ss
    real,allocatable :: field(:,:,:,:)
    character :: file*255,datfile*255,title*1000,vars(1)*40,lvars(1)*200,svars(1)*200, &
        units(1)*40,lz(3)*10,ltime*100,history*20000,cell_methods(1)*100,metadata(2,100)*2000
    logical :: lwrite
    integer :: leap

    if ( command_argument_count() /= 4 ) then
        print *,'usage: averagefieldspace infield avex avey outfield'
        call exit(-1)
    endif

    lwrite = .false. 
    minfac = 0.3
    call get_command_argument(2,file)
    read(file,*) avex
    call get_command_argument(3,file)
    read(file,*) avey
    call get_command_argument(1,file)
    call getmetadata(file,mens1,mens,ncid,datfile,nxmax,nx &
        ,xx,nymax,ny,yy,nzmax,nz,zz,lz,nt,nperyear,firstyr,firstmo &
        ,ltime,undef,endian,title,history,1,nvars,vars,ivars &
        ,lvars,svars,units,cell_methods,metadata,lwrite)
!   range of years
    call getlastyr(firstyr,firstmo,nt,nperyear,lastyr)

!   allocate field

    nxf = nx
    nyf = ny
    nzf = nz
    if ( nz > 1 ) then
        write(0,*) 'averagefieldspace: error: cannot handle fields with vertical extend ',nz,' yet'
        call exit(-1)
    end if
    if ( mens > mens1 ) then
        write(0,*) 'averagefieldspace: error: cannot handle muliple ensemble members ',mens1,mens,' yet'
        call exit(-1)
    end if
    if ( lwrite ) print *,'allocating field ',nx,ny,nperyear,firstyr,lastyr
    allocate(field(nx,ny,nperyear,firstyr:lastyr))

!   read field

    if ( ncid == -1 ) then
        call readdatfile(datfile,field,nx,ny,nx,ny,nperyear,firstyr &
            ,lastyr,firstyr,firstmo,nt,undef,endian,lwrite,firstyr &
            ,lastyr,1,1)
    else
        call readncfile(ncid,field,nx,ny,nx,ny,nperyear,firstyr &
            ,lastyr,firstyr,firstmo,nt,undef,lwrite,firstyr,lastyr &
            ,ivars)
        status = nf_close(ncid)
    end if

!   average

    call spatialaverage(field,xx,yy,nx,ny,nperyear,firstyr,lastyr,avex,avey,lwrite)

!   write output

    nnt = nt
    if ( nperyear == 366 ) then
        nnt = 0
        yr = firstyr
        mo = firstmo - 1
        do it=1,nt
            nnt = nnt + 1
            mo = mo + 1
            if ( mo > nperyear ) then
                mo = mo - nperyear
                yr = yr + 1
            endif
            if ( mo == 60 ) then
                if ( leap(yr) == 1 ) then
                    nnt = nnt - 1
                end if
            endif
        end do
    end if
    undef = 3e33
    write(title(len_trim(title)+2:),'(a,i2,a,i2,a)') ' averaged over ',avex,'x',avey,' grid boxes'
    call coarsen_geospatial_metadata(metadata,avex,avey)
    call get_command_argument(4,file)
    i = index(file,'.ctl')
    if ( i /= 0 ) then
        datfile = file(1:i)
        datfile(i:) = '.grd'
        call writectl(file,datfile,nx,xx,ny,yy,nz,zz, &
            nnt,nperyear,firstyr,firstmo,undef,title,nvars,vars,ivars, &
            lvars,units)
        open(1,file=datfile,access='direct',recl=recfa4*nx*ny)
        yr = firstyr
        mo = firstmo
        do it=1,nt
            write(1,rec=it) ((field(i,j,mo,yr),i=1,nx),j=1,ny)
            mo = mo + 1
            if ( mo > nperyear ) then
                mo = mo - nperyear
                yr = yr + 1
            end if
            if ( nperyear == 366 .and. mo == 60 ) then
                if ( leap(yr) == 1 ) then
                    mo = mo + 1
                end if
            end if
        end do
    else if ( index(file,'.nc') /= 0 ) then
        ivars = 0
        call enswritenc(file,ncid,ntvarid,itimeaxis,ntmax,nx,xx,ny,yy, &
            nz,zz,lz,nnt,nperyear,firstyr,firstmo,ltime,undef,title, &
            history,nvars,vars,ivars,lvars,svars,units,cell_methods, &
            metadata,0,0)
        yr = firstyr
        mo = firstmo
        do it=1,nnt
            if ( lwrite ) print *,'writing ',yr,mo,it,nnt,nt
            call keepalive1('Writing time step ',it,nnt)
            call writencslice(ncid,ntvarid,itimeaxis,ntmax,ivars &
                ,field(1,1,mo,yr),nxf,nyf,nzf,nx,ny,nz,it,1)
            mo = mo + 1
            if ( mo > nperyear ) then
                mo = mo - nperyear
                yr = yr + 1
            endif
            if ( nperyear == 366 .and. mo == 60 ) then
                if ( leap(yr) == 1 ) then
                    mo = mo + 1
                end if
            end if
        end do
        status = nf_close(ncid) ! otherwise the buffer won't be flushed...
    else
        write(0,*) 'need .ctl or .nc in output file name'
        call exit(-1)
    end if
end program averagefieldspace
