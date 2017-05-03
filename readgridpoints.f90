subroutine readgridpoints(data,ids,npermax,yrbeg,yrend,nensmax, &
    nperyear,mens1,mens,nens1,nens2,outvar,outunits,lstandardunits,lwrite)
!
!   read the grid points of a netcdf file into a set of time series
!
    implicit none
    integer,parameter :: nxmax=1000,nymax=1000
    integer   :: npermax,yrbeg,yrend,nensmax,nperyear,mens1,mens,nens1,nens2
    real      :: data(npermax,yrbeg:yrend,0:nensmax)
    logical   :: lwrite,lstandardunits
    character :: ids(0:nensmax)*(*),outvar*(*),outunits*(*)
    integer   :: ncid,nx,ny,nz,nt,firstyr,firstmo,endian,nvars,ivars(6,1)
    integer   :: nens,ix,iy,iz,yr,mo,iens,lastyr
    real      :: xx(nxmax),yy(nymax),zz(1),undef
    real,allocatable :: field(:,:,:,:,:,:)
    logical   :: validdata
    character :: file*1024,datfile*1014,lz(3)*20,ltime*100,title*1024,history*20000, &
        vars(1)*40,svars(1)*80,lvars(1)*80,units(1)*80,cell_methods(1)*128

    call getarg(2,file)
    call getmetadata(file,mens1,mens,ncid,datfile,nxmax,nx      &
        ,xx,nymax,ny,yy,1,nz,zz,lz,nt,nperyear,firstyr,firstmo  &
        ,ltime,undef,endian,title,history,1,nvars,vars,ivars    &
        ,lvars,svars,units,cell_methods,lwrite)
    lastyr = firstyr + (firstmo+nt-2)/nperyear
    mens1 = max(nens1,nens1)
    mens = min(mens,nens2)
    allocate(field(nx,ny,nz,nperyear,firstyr:lastyr,0:mens))
    field = 3e33
    call readfield(ncid,file,datfile,field,nx,ny,nz &
        ,nperyear,firstyr,lastyr,mens1,mens,nx,ny,nz,nperyear,firstyr,lastyr &
        ,firstyr,firstmo,nt,undef,endian,vars,units,lstandardunits,lwrite)

    iz = 1
    nens = 0
    do ix=1,nx
        do iy=1,ny
            do iens=mens1,mens
                validdata = .false.
                do yr=yrbeg,yrend
                    if ( yr < firstyr .or. yr > lastyr ) then
                        data(1:nperyear,yr,nens) = 3e33
                    else
                        do mo=1,nperyear
                            data(mo,yr,nens) = field(ix,iy,iz,mo,yr,iens)
                            if ( data(mo,yr,nens) < 1e33 ) validdata = .true.
                        end do
                    end if
                end do            
                if ( validdata ) then
                    write(ids(nens),'(f9.3,a,f9.3,a)') yy(iy),'N,',xx(ix),'E'
                    if ( nens < nensmax ) then
                        nens = nens + 1
                    else
                        write(0,*) 'readgridpoints: error: too many grid points with data, can only handle ',nensmax
                        call exit(-1)
                    end if
                end if
            end do
        end do
    end do
    ! shift meaning from input to output
    mens1 = 0
    mens = nens-1
    outunits = units(1)
    outvar = vars(1)
end subroutine

