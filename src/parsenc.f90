subroutine parsenc(file,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax &
    ,nz,zz,nt,nperyear,firstyr,firstmo,undef,title,nvarmax &
    ,ntvars,vars,ivars,lvars,units)

!   old entry point

    implicit none
    integer :: ntmax,nensmax
    parameter (ntmax=5000000,nensmax=750)
    integer :: ncid,nxmax,nymax,nzmax,nx,ny,nz,nt,nperyear,firstyr &
        ,firstmo,nvarmax,ntvars,ivars(6,nvarmax)
    real :: xx(nxmax),yy(nymax),zz(nzmax),undef
    character file*(*),title*(*),vars(nvarmax)*(*) &
        ,lvars(nvarmax)*(*),units(nvarmax)*(*)
    integer :: nens1,nens2,i,j
    character lz(3)*20,svars(100)*100,ltime*120,history*50000, &
        cell_methods(100)*100,metadata(2,100)*2000
    logical :: tdefined(ntmax)
    nens1 = 0
    nens2 = 0
    call ensparsenc(file,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax &
        ,nz,zz,lz,nt,nperyear,firstyr,firstmo,ltime,tdefined,ntmax &
        ,nens1,nens2,undef,title,history,nvarmax,ntvars,vars,ivars &
        ,lvars,svars,units,cell_methods,metadata)
    if ( nens1 /= 0 .or. nens2 /= 0 ) then
        write(0,*) 'parsenc: error: found ensemble in file ',trim(file)
        call exit(-1)
    end if
end subroutine parsenc

subroutine ensparsenc(file,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax &
    ,nz,zz,lz,nt,nperyear,firstyr,firstmo,ltime,tdefined,ntmax &
    ,nens1,nens2,undef,title,history,nvarmax,ntvars,vars,ivars &
    ,lvars,svars,units,cell_methods,metadata)

!   extract field metainformation from NetCDF file
!   BUGS: cannot handle multiple latide,longitude axes yet.

    implicit none
    include 'netcdf.inc'
    integer :: nensmax
    parameter (nensmax=750)
!   arguments
    integer :: ncid,nxmax,nymax,nzmax,nx,ny,nz,nt,nperyear,firstyr &
        ,firstmo,ntmax,nens1,nens2,nvarmax,ntvars,ivars(6,nvarmax)
    real :: xx(nxmax),yy(nymax),zz(nzmax),undef
    character file*(*),title*(*),history*(*),lz(3)*(*),ltime*(*), &
        vars(nvarmax)*(*),lvars(nvarmax)*(*),svars(nvarmax)*(*), &
        units(nvarmax)*(*),cell_methods(nvarmax)*(*),metadata(2,100)*(*)
    logical :: tdefined(ntmax)
!   local variables
    integer :: status,ndims,nvars,ngatts,unlimdimid,varid,xtype &
        ,ndimvar,dimids(nf_max_var_dims),natts,dimid,len,ix,iy,iz &
        ,it,ie,i,j,n,l,iperyear
    real*8,allocatable :: tt(:)
    character :: name*(nf_max_name),dimname*(nf_max_name),clwrite*10,axis*2,host*100
    logical :: lwrite,foundtime
    integer :: leap

    allocate(tt(ntmax))
    foundtime = .false. 
    lwrite = .false.
    tdefined(1) = .true.
    call getenv('PARSENC_LWRITE',clwrite)
    if ( index(clwrite,'T') + index(clwrite,'t') > 0 ) then
        lwrite = .true. 
    end if
    xx(1) = 0
    yy(1) = 0
    zz(1) = 0
    lz = ' '
    ltime = ' '
    ivars = 0

!   open file

    if ( ncid == 0 ) then
        if ( lwrite ) print *,'parsenc: opening file ',trim(file)
        status = nf_open(trim(file),nf_nowrite,ncid)
        if ( status /= nf_noerr ) call handle_err(status,trim(file))
        if ( lwrite ) print *,'parsenc: opened with ncid = ',ncid
    else
        if ( lwrite ) print *,'parsenc: already open with ncid = ',ncid
    end if
    call gettitle(ncid,title,lwrite)
    call gettextattopt(ncid,nf_global,'history',history,lwrite)
    if ( history == ' ' ) call gettextattopt(ncid,nf_global,'History',history,lwrite)
    if ( history == ' ' ) call gettextattopt(ncid,nf_global,'HISTORY',history,lwrite) ! Yes, unfortunately...
    call getglobalatts(ncid,metadata,lwrite)
    call getnumbers(ncid,ndims,nvars,ngatts,unlimdimid,lwrite)
    call getdims(ncid,ndims,ix,nx,nxmax,iy,ny,nymax,iz,nz,nzmax,it &
        ,nt,ntmax,ie,nens1,nens2,nensmax,lwrite)

!   loop over variables

    if ( lwrite ) print *,'parsenc: loop over variables 1-',nvars
    ntvars = 0
    nperyear = 0
    do varid=1,nvars
!       get axis information (if any)
        call gettextattopt(ncid,varid,'axis',axis,lwrite)
        call tolower(axis)
!       get dimensions of variable
        status = nf_inq_var(ncid,varid,name,xtype,ndimvar,dimids,natts)
        if ( status /= nf_noerr ) call handle_err(status,'nf_inq_var')
        if ( lwrite ) then
            print *,'parsenc: variable: ',varid
            print *,'         name:     ',trim(name)
            print *,'         dims:     ',ndimvar,':', &
            (dimids(i),i=1,ndimvar)
            print *,'         natts:    ',natts
            print *,'         axis:     ',axis
        end if
        if ( index(name,'_bnd') /= 0 .or. index(name,'_bound') /= 0 ) then
            if ( lwrite ) print *,'parsenc: disregarding boundary ',trim(name)
            cycle
        end if
        if ( name == 'average_T1' .or. name == 'average_T2' .or. name == 'average_DT') then
            if ( lwrite ) print *,'parsenc: disregarding extra time axis information ',trim(name)
            cycle
        end if
        if ( name == 'time_weights' ) then
            if ( lwrite ) print *,'parsenc: disregarding weights ',trim(name)
            cycle
        end if
        if ( index(name,'_FILENAME') /= 0 ) then
            if ( lwrite ) print *,'parsenc: disregarding boundary ',trim(name)
            cycle
        end if
!       what kind of variable do we have?
        if ( ndimvar == 1 .and. dimids(1) == ix ) then
            call getdiminfo('x',ncid,varid,xx,nx,lwrite)
            call makelonreasonable(xx,nx)
        elseif ( ndimvar == 1 .and. dimids(1) == iy ) then
            call getdiminfo('y',ncid,varid,yy,ny,lwrite)
        elseif ( ndimvar == 1 .and. dimids(1) == iz ) then
            call getzdiminfo('z',ncid,varid,zz,nz,lz,lwrite)
        elseif ( ndimvar == 1 .and. dimids(1) == ie ) then
            if ( lwrite ) print *,'renumbering ensemble members'// &
            ' to 0 ... nens-1'
        elseif ( ndimvar == 1 .and. dimids(1) == it .and. &
             .not. foundtime .and. ( axis == 't' .or. &
            name(1:4) == 'time' .or. name(1:4) == 'TIME' .or. &
            name == 'T' .or. name == 't' .or. name == 'T1' .or. &
            (name(1:2) == 't_' .and. name /= 't_ref') .or. &
            name == 'Time'  .or. name == 'Tims' .or. name == 'day' )) then
            foundtime = .true. 
!           (it could have been a timeseries)
            if ( lwrite ) print *,'parsenc: found time axis'
            status = nf_get_var_double(ncid,varid,tt)
            if ( status /= nf_noerr ) call handle_err(status &
                ,'nf_get_var_real(tt)')
            if ( lwrite ) print *,'tt(1-5) = ',(tt(i),i=1,min(nt,5))
            call getperyear(ncid,varid,tt,nt,firstmo,firstyr &
                ,nperyear,iperyear,ltime,tdefined,ntmax,lwrite)
        else
            n = 0
            do i=1,ndimvar
                if ( it /= 0 .and. dimids(i) == it .and. &
                index(name,'bounds') == 0 ) then
                    n = n+1
                    if ( lwrite ) print *,'parsenc: time-varying variable ',varid
                end if
            end do
            if ( lwrite ) then
                print *,'         checking for a lat-lon var'
                print *,'         it,n          = ',it,n
                print *,'         ndimvar,ix,iy = ',ndimvar,ix,iy
                print *,'         dimids(1,2)   = ',dimids(1) &
                ,dimids(2)
            end if
            if ( it <= 0 .and. n == 0 .and. ndimvar == 2 &
             .and. ix /= 0 .and. iy /= 0 .and. &
            ( dimids(1) == ix .and. dimids(2) == iy .or. &
            dimids(2) == ix .and. dimids(1) == iy ) ) then
                n = n+1
                if ( lwrite ) print * &
                ,'parsenc: lat-lon variable ',varid
            end if
            if ( n == 1 ) then
                call addonevariable(ncid,varid,name,ntvars,nvarmax &
                ,ndimvar,dimids,ix,iy,iz,it,ie,vars,ivars,lvars &
                ,svars,units,cell_methods,undef,lwrite)
            end if           ! one time variable?
        end if               ! variable-recognition case
    end do

!   if there is no spatial dimension  we may be able to find it in the metadata

    if ( nx == 1 .and. ix == -1 ) then
        call getlonfrommetadata(ncid,xx(1),lwrite)
    end if
    if ( ny == 1 .and. iy == -1 ) then
        call getlatfrommetadata(ncid,yy(1),lwrite)
    end if

!   the rest cannot handle nz=0

    if ( nz == 0 ) then
        nz = 1
        zz(1) = 0
    end if

!   nor nperyear = 0

    if ( nperyear == 0 ) then
        nperyear = 1
        firstyr = 1
        firstmo = 1
        nt = 1
    end if
    deallocate(tt)

end subroutine ensparsenc

