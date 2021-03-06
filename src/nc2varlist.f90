program nc2varlist

!   reads the metadat of a netcdf file and converts it to an HTML
!   list suitable for use in the Climate Explorer.
!    Replaces a very hairy sed-stateemnt w/o units

    implicit none
    include 'netcdf.inc'
    include 'params.h'
    integer,parameter :: nvarmax=100
    integer :: ncid,nx,ny,nz,nt,nperyear,firstyr &
        ,firstmo,ntvars,ivars(6,nvarmax),ivar
    real :: xx(nxmax),yy(nymax),zz(nzmax),undef
    character file*255,title*10000,vars(nvarmax)*40 &
        ,lvars(nvarmax)*200,units(nvarmax)*60

    if ( command_argument_count() /= 1 ) then
        write(0,*) 'usage: nc2varlist file'
        call exit(-1)
    end if
    call get_command_argument(1,file)
    ncid = 0
    call parsenc(file,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax &
        ,nz,zz,nt,nperyear,firstyr,firstmo,undef,title,nvarmax &
        ,ntvars,vars,ivars,lvars,units)

    do ivar=1,ntvars
        print '(8a)','<input type="radio" class="formradio" ', &
            'name="var" value="',trim(vars(ivar)),'">', &
            trim(lvars(ivar)),' [',trim(units(ivar)),']<br>'
    end do

end program nc2varlist
