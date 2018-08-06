program wetbulb_field
!
!   compute an approximation to the maximum wet bulb temperature of a day
!   by using the max dew-point temperature and the maximum temperature (and surface pressure)
!
    implicit none
    include 'params.h'
    integer nvarmax
    parameter(nvarmax=1)
    integer i,j,yr,mo,idtmax,idtmean,idtdew,idpres,idtwet,status,ntvarid
    integer nx,ny,nz,nx1,ny1,nz1,it,nperday,nvars,lastyr
    integer ivtmax(6,nvarmax),ivtdew(6,nvarmax),ivpres(6,nvarmax),ivars(2,nvarmax)
    integer nperyear,nt,firstyr,firstmo,ntvars,nperyear1,nt1,firstyr1,firstmo1
    integer,allocatable :: itimeaxis(:)
    real xx(nxmax),yy(nymax),zz(nzmax),xx1(nxmax),yy1(nxmax),zz1(nzmax),undef,undef1,undef2
    real,allocatable ::  tmax(:,:,:,:),tdew(:,:,:,:),pres(:,:,:,:),twet(:,:,:,:)
    logical lwrite,lstandardunits
    character tmaxfile*255,tdewfile*255,presfile*255,twetfile*255
    character title*255,vars(nvarmax)*50,lvars(nvarmax)*100
    character tmaxunits(nvarmax)*50,tdewunits(nvarmax)*50,presunits(nvarmax)*50,units(nvarmax)*50
    integer,external :: leap,nf_close
    real,external :: wetbulbdew
    
    if ( command_argument_count() /= 4 ) then
        write(0,*) 'usage: wetbulbfield Tmax Tdew pressure Twetbulb'
        call exit(-1)
    end if
    
    call get_command_argument(1,tmaxfile)
    call get_command_argument(2,tdewfile)
    call get_command_argument(3,presfile)
    call get_command_argument(4,twetfile)
    
    lwrite = .false.
    idtmax = 0
    call parsenc(trim(tmaxfile),idtmax,nxmax,nx,xx,nymax,ny,yy,nzmax &
     &        ,nz,zz,nt,nperyear,firstyr,firstmo,undef,title,nvarmax &
     &        ,ntvars,vars,ivtmax,lvars,tmaxunits)
    if ( tmaxunits(1) /= 'K' .and. tmaxunits(1) /= 'Celsius' .and. tmaxunits(1) /= 'degrees_celsius' ) then
        write(0,*) 'wetbulbfield: error: expecting Celsius or K for Tmax, not ',tmaxunits
        call exit(-1)
    end if
    if ( nz > 1 ) then
        write(0,*) 'wetbulb_field: error: only for horizontal fields, not nz = ',nz
        call exit(-1)
    end if
    idtdew = 0
    call parsenc(trim(tdewfile),idtdew,nxmax,nx1,xx1,nymax,ny1,yy1,nzmax &
     &        ,nz1,zz1,nt1,nperyear1,firstyr1,firstmo1,undef1,title,nvarmax &
     &        ,ntvars,vars,ivtdew,lvars,tdewunits)
    if ( tdewunits(1) /= 'K' .and. tdewunits(1) /= 'Celsius' .and. tdewunits(1) /= 'degrees_celsius' ) then
        write(0,*) 'wetbulbfield: error: expecting Celsius or K for Tdew, not ',tdewunits
        call exit(-1)
    end if
    if ( nz1 > 1 ) then
        write(0,*) 'wetbulb_field: error: only for horizontal fields, not nz = ',nz
        call exit(-1)
    end if
    call checkgridequal(nx,ny,xx,yy, nx1,ny1,xx1,yy1)
    nt = min(nt,nt1)
    if ( nperyear1 /= nperyear .or. firstyr1 /= firstyr .or. firstmo1 /= firstmo ) then
        write(0,*) 'wetbulb_field: error: timescale or start time not equal: ', &
        &   nperyear1,nperyear,firstyr1,firstyr, firstmo1,firstmo
        call exit(-1)
    end if
    idpres = 0
    call parsenc(trim(presfile),idpres,nxmax,nx1,xx1,nymax,ny1,yy1,nzmax &
     &        ,nz1,zz1,nt1,nperyear1,firstyr1,firstmo1,undef2,title,nvarmax &
     &        ,ntvars,vars,ivpres,lvars,presunits)
    if ( presunits(1) /= 'Pa' .and. presunits(1) /= 'hPa' .and. presunits(1) /= 'mb' ) then
        write(0,*) 'wetbulbfield: error: expecting Pa, hPa or mb for pressure, not ',presunits
        call exit(-1)
    end if
    if ( nz1 > 1 ) then
        write(0,*) 'wetbulb_field: error: only for horizontal fields, not nz = ',nz
        call exit(-1)
    end if
    call checkgridequal(nx,ny,xx,yy, nx1,ny1,xx1,yy1)
    nt = min(nt,nt1)
    if ( nperyear1 /= nperyear .or. firstyr1 /= firstyr .or. firstmo1 /= firstmo ) then
        write(0,*) 'wetbulb_field: error: timescale or start time not equal: ', &
        &   nperyear1,nperyear,firstyr1,firstyr, firstmo1,firstmo
        call exit(-1)
    end if
    
    lastyr = firstyr + (firstmo+nt-2)/nperyear
    
    allocate(tmax(nx,ny,nperyear,firstyr:lastyr),pres(nx,ny,nperyear,firstyr:lastyr))
    allocate(tdew(nx,ny,nperyear,firstyr:lastyr),twet(nx,ny,nperyear,firstyr:lastyr))
    
    write(0,*) 'Reading ',trim(tmaxfile)
    nt1 = nt
    call readncfile(idtmax,tmax,nx,ny,nx,ny,nperyear,firstyr,lastyr,firstyr,firstmo,nt1, &
    &   undef,lwrite,firstyr,lastyr,ivtmax)
    if ( tmaxunits(1) == 'K' ) then
        tmaxunits = 'Celsius'
        tmax = tmax - 273.15
    end if
    write(0,*) 'Reading ',trim(tdewfile)
    nt1 = nt
    call readncfile(idtdew,tdew,nx,ny,nx,ny,nperyear,firstyr,lastyr,firstyr,firstmo,nt1, &
    &   undef1,lwrite,firstyr,lastyr,ivtdew)
    if ( tdewunits(1) == 'K' ) then
        tdewunits = 'Celsius'
        tdew = tdew - 273.15
    end if
    write(0,*) 'Reading ',trim(presfile)
    nt1 = nt
    call readncfile(idpres,pres,nx,ny,nx,ny,nperyear,firstyr,lastyr,firstyr,firstmo,nt1, &
    &   undef2,lwrite,firstyr,lastyr,ivpres)
    if ( presunits(1) == 'Pa' ) then
        tmaxunits = 'hPa'
        pres = pres/100
    end if

    write(0,*) 'Computing'
    do yr=firstyr,lastyr
        do mo=1,nperyear
            do j=1,ny
                do i=1,nx
                    twet(i,j,mo,yr) = wetbulbdew(tmax(i,j,mo,yr),pres(i,j,mo,yr),tdew(i,j,mo,yr),lwrite)
                end do  ! i
            end do  ! j
        end do  ! mo
    end do  ! yr
    
    write(0,*) 'Writing ',trim(twetfile)
    nvars = 1
    vars(1) = 'Twetbulb'
    lvars(1) = 'daily maximum wetbulb temperature estimated from Tmax, Tdew and pressure'
    units(1) = 'Celsius'
    ivars(1,1) = 1
    allocate(itimeaxis(nt1))
    call writenc(twetfile,idtwet,ntvarid,itimeaxis,nt1,nx,xx,ny,yy,nz,zz,nt,nperyear,firstyr,firstmo, &
    &   3e33,title,nvars,vars,ivars,lvars,units,0,0)
    it = 0
    nperday = max(1,nint(nperyear/366.24))
    do yr=firstyr,lastyr
        do mo=1,nperyear
            if ( yr == firstyr .and. mo < firstmo ) cycle
            if ( nperyear == 366*nperday .and. leap(yr) == 2 .and. (mo-1)/nperday == 59 ) cycle
            it = it + 1
            if ( it > nt ) cycle
            call writencslice(idtwet,ntvarid,itimeaxis,nt,ivars,twet(1,1,mo,yr),nx,ny,nz,nx,ny,nz,it,0)
        end do ! mo
    end do ! yr
    status = nf_close(idtwet)

end program