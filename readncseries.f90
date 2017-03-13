subroutine readncseries(file,data,npermx,nperyear,yr1,yr2,ncid,var,units,lwrite)

!   read the data in a 1D netCDF file into data.

    implicit none
    integer :: npermx,nperyear,yr1,yr2,ncid
    real :: data(npermx,yr1:yr2)
    character file*(*),var*(*),units*(*)

    include 'param.inc'
    integer :: mxmax,mymax,mzmax,nvmax
    parameter(mxmax=1,mymax=1,mzmax=1,nvmax=1)
    integer :: i,j,nx,ny,nz,nt,firstyr,firstmo,nvars,ivars(6,nvmax),lastyr
    logical :: lwrite
    real :: xx(mxmax),yy(mymax),zz(mzmax),undef
    real,allocatable :: ddata(:,:,:,:)
    character title*512,vars(nvmax)*20,lvars(nvmax)*80
            
    call parsenc(file,ncid,mxmax,nx,xx,mymax,ny,yy,mzmax &
        ,nz,zz,nt,nperyear,firstyr,firstmo,undef,title,nvmax,nvars &
        ,vars,ivars,lvars,units)
    if ( lwrite ) print *,'readncseries: nperyear = ',nperyear
    var = vars(1)
    if ( nvars /= 1 ) then
        write(0,*) 'readncseries: error: not just one time-dependent variable, found ',nvars
        call exit(-1)
    endif
    if ( nx > 1 .and. ivars(2,1) /= 0 ) then
        write(0,*) 'readncseries: error: found x-dependent variable, found ',nx
        call exit(-1)
    endif
    if ( ny > 1 .and. ivars(3,1) /= 0 ) then
        write(0,*) 'readncseries: error: found y-dependent variable, found ',ny
        call exit(-1)
    endif
    if ( nz > 1 .and. ivars(4,1) /= 0 ) then
        write(0,*) 'readncseries: error: found z-dependent variable, found ',nz
        call exit(-1)
    endif
    lastyr = firstyr + (nt+firstmo-2)/nperyear
    allocate(ddata(1,1,nperyear,firstyr:lastyr))
    call readncfile(ncid,ddata,1,1,nx,ny,nperyear,firstyr,lastyr, &
        firstyr,firstmo,nt,undef,lwrite, &
        max(yr1,yrbeg,firstyr),min(yr2,yrend,lastyr),ivars)
    do i=max(yr1,yrbeg,firstyr),min(yr2,yrend,lastyr)
        do j=1,nperyear
            data(j,i) = ddata(1,1,j,i)
        enddo
    enddo
    deallocate(ddata)
end subroutine readncseries