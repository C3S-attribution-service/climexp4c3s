subroutine readncslice(ncid,ivars,it,field,nx,ny,nz,lwrite)
!
!       write a slice of a netcdf file at time it
!
    implicit none
    integer ncid,ivars(6),it,nx,ny,nz
    real field(nx,ny,nz)
    logical lwrite
    call ensreadncslice(ncid,ivars,it,1,field,nx,ny,nz,lwrite)
end subroutine

subroutine ensreadncslice(ncid,ivars,it,ie,field,nx,ny,nz,lwrite)
!
!   write a slice of a netcdf file at time it and ensemble member ie
!
    implicit none
    include 'netcdf.inc'
    integer ncid,ivars(6),it,ie,nx,ny,nz
    real field(nx,ny,nz)
    logical lwrite
    integer i,k,stride(5),start(5),count(5),status
    k = 0
    stride = 1
    if ( ivars(2).gt.0 ) then
        k = k + 1
        start(k) = 1
        count(k) = max(1,nx)
    endif
    if ( ivars(3).gt.0 ) then
        k = k + 1
        start(k) = 1
        count(k) = max(1,ny)
    endif
    if ( ivars(4).gt.0 ) then
        k = k + 1
        start(k) = 1
        count(k) = max(1,nz)
    endif
    if ( ivars(6).gt.0 ) then
        k = k + 1
        start(k) = ie
        count(k) = 1
    endif
    if ( ivars(5).gt.0 ) then
        k = k + 1
        start(k) = it
        count(k) = 1
    endif
    if ( lwrite ) then
        print '(a,6i6)','readncslice: variable = ',ivars
        print '(a,6i6)','             nx,ny,nz = ',nx,ny,nz
        print '(a,6i6)','             startvec = ',(start(i),i=1,k)
        print '(a,6i6)','             countvec = ',(count(i),i=1,k)
    endif
    status = nf_get_vara_real(ncid,ivars(1),start,count,field)
    if ( status.ne.nf_noerr ) then
        call handle_err(status,'readncslice nf_get_vara_real')
    endif
    if ( lwrite ) print *,'field(',1+nx/2,1+ny/2,1+nz/2,') = ' &
 &       ,field(1+nx/2,1+ny/2,1+nz/2)
end subroutine
