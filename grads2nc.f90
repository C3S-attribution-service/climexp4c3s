program grads2nc

!   convert a GrADS ctl/dat file combo to netcdf

    implicit none
    include 'params.h'
    include 'netcdf.inc'
    integer,parameter :: nvarmax=1000,ntmax=366*1000
    integer :: nx,ny,nz,nt,nperyear,yrbegin,mobegin,endian,nvars &
        ,ivars(2,nvarmax)
    integer :: status,ncid,i,j,k,itimeaxis(ntmax),dpm(12,2), &
        ix,iy,iz,it,iv,localendian,ntvarid
    real :: xx(nxmax),yy(nymax),zz(nzmax),undef
    real,allocatable :: data(:,:,:)
    character :: vars(nvarmax)*10,lvars(nvarmax)*80,units(nvarmax)*40
    character(255) :: ctlfile,datfile,ncfile,title
    logical :: lwrite
    integer,external :: get_endian
    data dpm &
    /31,28,31,30,31,30,31,31,30,31,30,31 &
    ,31,29,31,30,31,30,31,31,30,31,30,31/

    lwrite = .false. 
    if ( command_argument_count() /= 2 ) then
        write(0,*) 'Usage: grads2nc infile.ctl outfile.nc'
        stop
    endif
    localendian = get_endian()

!   metadata

    call get_command_argument(1,ctlfile)
    call parsectl(ctlfile,datfile,nxmax,nx,xx,nymax,ny,yy,nzmax &
        ,nz,zz,nt,nperyear,yrbegin,mobegin,undef,endian,title &
        ,nvarmax,nvars,vars,ivars,lvars,units)
    call get_command_argument(2,ncfile)
    call writenc(ncfile,ncid,ntvarid,itimeaxis,ntmax,nx,xx,ny,yy &
        ,nz,zz,nt,nperyear,yrbegin,mobegin,undef,title,nvars &
        ,vars,ivars,lvars,units,0,0) ! no ensembles yet

!   data

    allocate(data(nx,ny,nz))
    open(1,file=datfile,access='direct',recl=4*nx*ny*nz)
    if ( lwrite ) print *,'endian,localendian = ',endian,localendian
    do it=1,nt
        do iv=1,nvars
            read(1,rec=nvars*(it-1)+iv) &
            (((data(ix,iy,iz),ix=1,nx),iy=1,ny),iz=1,nz)
            if ( localendian*endian == -1 ) then
                if ( lwrite ) print *,'before ',(((data(ix,iy,iz), &
                ix=1,nx),iy=1,ny),iz=1,nz)
                do iz=1,nz
                    do iy=1,ny
                        call swapbyte4(data(1,iy,iz),nx)
                    end do
                end do
                if ( lwrite ) print *,'after  ',(((data(ix,iy,iz) &
                    ,ix=1,nx),iy=1,ny),iz=1,nz)
            end if
            call writencslice(ncid,ntvarid,itimeaxis,ntmax, &
            ivars(1,iv),data,nx,ny,nz,nx,ny,nz,it,1)
        end do
    enddo
    status = nf_close(ncid)
end program grads2nc

