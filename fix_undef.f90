program fix_undef

!   fix wrong undef declarations

    implicit none
    include 'params.h'
    include 'netcdf.inc'
    include 'getopts.inc'
    integer,parameter :: recfa4=4
    real,parameter :: absent=3e33,eps=1e-5
    integer :: nx,ny,nz,nt,firstyr,lastyr,firstmo,nvars, &
        ivars(2,nvmax),jvars(6,nvmax),ncid,endian, &
        status,nperyear,mens,mens1,i1,i2,n1,n2,n,irec,it,ntvarid
    integer :: jx,jy,jz,mo,yr,nunique,i,j,nmax
    integer,allocatable :: nn(:),itimeaxis(:)
    real :: xx(nxmax),yy(nymax),zz(nzmax),undef,u1,u2,uu1,uu2
    real,allocatable :: field(:,:,:,:,:),unique(:)
    character :: infile*255,datfile*255,outfile*255,line*255,lz(3)*20 &
        ,vars(nvmax)*15,lvars(nvmax)*255,svars(nvmax)*80 &
        ,title*255,history*20000,units(nvmax)*40, &
        cell_methods(nvmax)*100,ltime*120,yesno*1, &
        metadata(2,100)*2000
    logical :: lexist
    integer,external :: leap

!   process command line

    if ( command_argument_count() < 2 ) then
        write(0,*) 'usage: fix_undef infile.[ctl|nc] ' &
        //'[lwrite] outfile.[ctl|nc]'
        stop
    end if
    call get_command_argument(1,infile)
    call getmetadata(infile,mens1,mens,ncid,datfile,nxmax,nx &
        ,xx,nymax,ny,yy,nzmax,nz,zz,lz,nt,nperyear,firstyr,firstmo &
        ,ltime,undef,endian,title,history,nvmax,nvars,vars,ivars &
        ,lvars,svars,units,cell_methods,metadata,lwrite)
    if ( mens1 /= 0 .or. mens /= 0 ) then
        write(0,*) 'fix_undef: cannot handle ensembles'
        call exit(-1)
    end if
    lastyr = firstyr + (firstmo+nt-2)/nperyear
    call getopts(2,iargc()-1,nperyear,yrbeg,yrend,.false.,mens1,mens)
    yr1 = firstyr
    yr2 = lastyr
    call get_command_argument(iargc(),outfile)

!   allocate field

    if ( lwrite ) print *,'allocating field(',nx,ny,nz,nperyear,firstyr,lastyr,')'
    allocate(field(nx,ny,nz,nperyear,firstyr:lastyr))
    nmax = nx*ny*nz
    allocate(unique(nmax),nn(nmax))

!   read field

    write(0,*) 'reading field...<p>'
    call readfield(ncid,infile,datfile,field,nx,ny,nz,nperyear &
        ,firstyr,lastyr,nens1,nens2,nx,ny,nz,nperyear,yr1,yr2 &
        ,firstyr,firstmo,nt,undef,endian,vars,units,lstandardunits &
        ,lwrite)

!   fill array with unique values and count them

    write(0,*) 'searching for undef...<p>'
    nunique = 0
    nn = 0
    n = 0
    do yr=firstyr,lastyr
        do mo=1,nperyear
            if ( lwrite ) print *,yr,mo,nunique,n
            do jz=1,nz
                do jy=1,ny
                    do jx=1,nx
                        if ( field(jx,jy,jz,mo,yr) /= absent ) then
                            do i=1,nunique
                                if ( abs(field(jx,jy,jz,mo,yr) - unique(i)) < eps*abs(unique(i)) ) then
!                                   match
                                    nn(i) = nn(i) + 1
                                    goto 100
                                end if
                            end do
!                           no match
                            if ( nunique < nmax ) then
                                nunique = nunique + 1
                                unique(nunique) = &
                                field(jx,jy,jz,mo,yr)
                            end if
                        100 continue
                            n = n + 1
                            if ( n > 6*nmax ) go to 200
                        end if
                    end do
                end do
            end do
        end do
    end do
200 continue

!   find two most frequent values

    if ( lwrite ) print *,'find two most frequently used values'
    u1 = 0
    u2 = 0
    n1 = 0
    n2 = 0
    uu1 = 0
    uu2 = 0
    do i=1,nunique
        if ( nn(i) > n2 ) then
            if ( nn(i) > n1 ) then
                u2 = u1
                n2 = n1
                u1 = unique(i)
                n1 = nn(i)
            else
                u2 = unique(i)
                n2 = nn(i)
            end if
        end if
        if ( abs(unique(i)) > uu2 ) then
            if ( abs(unique(i)) > uu1 ) then
                uu2 = uu1
                uu1 = abs(unique(i))
            else
                uu2 = abs(unique(i))
            end if
        end if
    end do
    if ( lwrite ) then
        print *,'n1,u1 = ',n1,u1
        print *,'n2,u2 = ',n2,u2
        print *,'uu1,uu2 ',uu1,uu2
    end if
!   heuristics - either much bigger or much more frequent
    if ( ( abs(uu1-abs(u1)) < 1e-4*uu1 .and. uu1 > 10*uu2 ) .or. &
         ( n1 > nx*n2 .and. uu1 > uu2 ) ) then
        write(0,*) 'Found wrong undef ',u1,', patching...<p>'
        do yr=firstyr,lastyr
            do mo=1,nperyear
                do jz=1,nz
                    do jy=1,ny
                        do jx=1,nx
                            if ( abs(field(jx,jy,jz,mo,yr) - u1) < eps*abs(u1) ) then
                                field(jx,jy,jz,mo,yr) = absent
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end if

!       write field

    write(0,*) 'writing output...<p>'
    if ( index(outfile,'.nc') == 0 ) then
        i = index(outfile,'.ctl')
        if ( i == 0 ) i = len_trim(outfile)+1
        datfile = outfile(:i-1)//'.grd'
        inquire(file=trim(outfile),exist=lexist)
        if ( lexist ) then
            print *,'outfile ',trim(outfile), &
            ' already exists, overwrite?'
            read(*,*) yesno
            if ( yesno == 'y' .or. yesno == 'Y' ) then
                open(1,file=trim(outfile))
                close(1,status='delete')
                open(1,file=trim(datfile))
                close(1,status='delete')
            end if
        end if
        if ( ncid /= -1 ) then
            ivars(1,1) = nz
            ivars(2,1) = 99
        end if
        call writectl(outfile,datfile,nx,xx,ny,yy,nz,zz, &
            nt,nperyear,firstyr,firstmo,3e33,title,1,vars,ivars &
            ,lvars,units)
        open(1,file=trim(datfile),form='unformatted',access='direct' &
            ,recl=nx*ny*nz*recfa4)
        yr = firstyr
        mo = firstmo
        do it=1,nt
            write(1,rec=it) (((field(jx,jy,jz,mo,yr),jx=1,nx),jy=1,ny),jz=1,nz)
            mo = mo + 1
            if ( mo > nperyear ) then
                mo = 1
                yr = yr + 1
            end if
        end do
        close(1)
    else
        if ( .not. allocated(itimeaxis) ) then
            allocate(itimeaxis(nt))
        end if
        call enswritenc(outfile,ncid,ntvarid,itimeaxis,nt,nx,xx &
            ,ny,yy,nz,zz,lz,nt,nperyear,firstyr,firstmo,ltime,3e33 &
            ,title,history,nvars,vars,ivars,lvars,svars,units &
            ,cell_methods,metadata,nens1,nens2)
        yr = firstyr
        mo = firstmo
        irec = 0
        do it=1,nt
            if ( nperyear == 366 .and. mo == 60 .and. &
            leap(yr) == 1 ) cycle
            irec = irec + 1
            call writencslice(ncid,ntvarid,itimeaxis,nt,ivars &
                ,field(1,1,1,mo,yr),nx,ny,nz,nx,ny,nz,irec,1)
            mo = mo + 1
            if ( mo > nperyear ) then
                mo = 1
                yr = yr + 1
            end if
        end do
        status = nf_close(ncid)
    end if
end program fix_undef
