    program catnc

!   concatenate netcdf files with the same x,y,z axes in the time direction
!   filling up the gaps with undefineds, preventing holes in the time axis.
!   (this differs from ncrcat and cdo copy).

    implicit none
    include 'params.h'
    include 'netcdf.inc'
    integer :: nvarmax,ntmax,fyr,lyr
    parameter(nvarmax=1,ntmax=366*1000,fyr=1850,lyr=2100)
    integer :: nx,ny,nz,nt,nnt,nperyear,yrbegin,mobegin,endian,nvars,ivars(2,nvarmax), &
        jvars(6,nvarmax),nperday,dy,dy1
    integer :: status,ncid,i,j,k,itimeaxis(ntmax),dpm(12,2), &
        ix,iy,iz,it,iv,localendian,ntvarid,nargs,fyrall,lyrall,f1yr,l1yr,iarg, &
        yr,mo,nens1,nens2
    real :: xx(nxmax),yy(nymax),zz(nzmax),undef
    real,allocatable :: field(:,:,:,:,:),onefield(:,:,:,:,:)
    character :: vars(nvarmax)*10,lvars(nvarmax)*80,svars(nvarmax)*80,units(nvarmax)*40
    character :: ncfile*255,title*1023,lz(3)*20,ltime*120,history*20000, &
        cell_methods(100)*100,metadata(2,100)*2000
    logical :: lwrite,tdefined(ntmax)
    integer,external :: leap
    data dpm &
    /31,28,31,30,31,30,31,31,30,31,30,31 &
    ,31,29,31,30,31,30,31,31,30,31,30,31/

    lwrite = .false. 
    if ( command_argument_count() < 2 ) then
        write(0,*) 'Usage: catnc infile1.nc infiel2.nc ... outfile.nc'
        call exit(-1)
    endif

!   metadata

    nargs = command_argument_count()
    fyrall =  9999
    lyrall = -9999
    do iarg=1,nargs-1
        call keepalive1('Concatenating file ',iarg,nargs-1)
        call get_command_argument(iarg,ncfile)
        if ( lwrite ) print *,'reading metadata ',trim(ncfile)
        ncid = 0
        call ensparsenc(ncfile,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax &
            ,nz,zz,lz,nt,nperyear,yrbegin,mobegin,ltime,tdefined,ntmax &
            ,nens1,nens2,undef,title,history,nvarmax,nvars,vars,jvars,lvars &
            ,svars,units,cell_methods,metadata)
        if ( nens1 /= nens2 ) then
            write(0,*) 'catnc: error: cannot handle ensembles yet'
            call exit(-10)
        end if
        if ( iarg == 1 ) then
            allocate(field(nx,ny,nz,nperyear,fyr:lyr))
            field = 3e33
        end if
        f1yr = yrbegin
        l1yr = yrbegin + (mobegin+nt-2)/nperyear
        if ( f1yr < fyr .or. l1yr > lyr ) then
            write(0,*) 'catnc: error: can only handle data between ',fyr,' and ',lyr, &
                ', not ',f1yr,' to ',l1yr
            call exit(-1)
        end if
        fyrall = min(fyrall,f1yr)
        lyrall = max(lyrall,l1yr)
        allocate(onefield(nx,ny,nz,nperyear,f1yr:l1yr))
        if ( lwrite ) print *,'reading data ',trim(ncfile)
        call zreadncfile(ncid,onefield,nx,ny,nz,nx,ny,nz,nperyear, &
            f1yr,l1yr,yrbegin,mobegin,nt,undef,lwrite,f1yr,l1yr,jvars)
        do yr=f1yr,l1yr
            do mo=1,nperyear
                do k=1,nz
                    do j=1,ny
                        do i=1,nx
                            if ( onefield(i,j,k,mo,yr) < 1e33 ) then
                                field(i,j,k,mo,yr) = onefield(i,j,k,mo,yr)
                            end if
                        end do
                    end do
                end do
            end do
        end do
        deallocate(onefield)
    end do        

    call get_command_argument(nargs,ncfile)
    nt = nperyear*(lyrall - fyrall + 1)
    if ( lwrite ) print *,'writing metadata'
    nnt = nt
    nperday = 1
    if ( nperyear >= 366 ) then
        nperday = nint(nperyear/366.)
        if ( nperyear/nperday == 366 ) then
            nnt = 0
            do yr=fyrall,lyrall
                if ( leap(yr) == 1 ) then
                    nnt = nnt + 365
                else
                    nnt = nnt + 366
                end if
            end do
        end if
    end if
    call enswritenc(ncfile,ncid,ntvarid,itimeaxis,ntmax,nx,xx,ny,yy &
        ,nz,zz,lz,nnt,nperyear,fyrall,1,ltime,undef,title,history,nvars &
        ,vars,ivars,lvars,svars,units,cell_methods,metadata,0,0) ! no ensembles yet
    it = 0
    if ( lwrite ) print *,'writing data'
    do yr=fyrall,lyrall
        do dy=1,nperyear
            if ( leap(yr) == 1 .and. dy/nperday == 60 ) cycle
            it = it + 1
            call keepalive1('Writing slice ',it,(lyrall-fyrall+1)*nperyear)
            call writencslice(ncid,ntvarid,itimeaxis,ntmax, &
                    ivars(1,1),field(1,1,1,dy,yr),nx,ny,nz,nx,ny,nz,it,1)
        end do
    enddo
    status = nf_close(ncid)
end program

