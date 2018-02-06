program subfieldseries

!   written to subtract global mean temperature from a field

    implicit none
    include 'params.h'
    include 'netcdf.inc'
    integer,parameter :: recfa4=4,ntmax=npermax*(yrend-yrbeg)
    integer :: n,ncid,ncid2,nx,ny,nz,nt,nper1,firstyr,firstmo, &
        nper2,nvars,ivars1(6,nvmax),ivars(6,nvmax),endian,endian2 &
        ,status,nperyear,lastyr,mens1,mens,itimeaxis(ntmax),nens1,nens2 &
        ,ntvarid
    integer :: i,j,jx,jy,irec
    real ::  xx(nxmax),yy(nymax),zz(nzmax),undef,data(npermax,yrbeg:yrend)
    real,allocatable :: field(:,:,:,:)
    logical :: lwrite,lexist
    character :: datfile*1023,outfile*1023,title*500,vars(1)*40, &
        lvars(1)*80,units(1)*80,var2*80,units2*80,lvar2*120, &
        svars(1)*80,history*50000,metadata(2,100)*1000,lz(3)*100, &
        svar2(1)*80,history2*50000,metadata2(2,100)*1000, &
        ltime*100,cell_methods(1)*80,seriesfile*1023,fieldfile*1023
    integer :: iargc

!   check arguments

    lwrite = .false. 
    n = iargc()
    if ( n /= 3 ) then
        print *,'usage: subfieldseries field.[ctl|nc] series.[nc|dat] outfield.cftl'
        call exit(-1)
    end if
    call getarg(1,fieldfile)
    call getmetadata(fieldfile,mens1,mens,ncid,datfile,nxmax,nx &
        ,xx,nymax,ny,yy,nzmax,nz,zz,lz,nt,nper1,firstyr,firstmo &
        ,ltime,undef,endian,title,history,1,nvars,vars,ivars1 &
        ,lvars,svars,units,cell_methods,metadata,lwrite)

    lastyr = firstyr + (firstmo+nt-2)/nper1
    if ( lwrite ) print *,'allocating field ',nx,ny,nper1,firstyr,lastyr
    allocate(field(nx,ny,nper1,firstyr:lastyr))

    call getarg(2,seriesfile)
    call readseriesmeta(seriesfile,data,npermax,yrbeg,yrend,nper2,var2,units2,lvar2,svar2, &
        history2,metadata2,.false.,lwrite)

    nperyear = max(nper1,nper2)
    if ( nper1 /= nper2 ) then
        write(0,*) 'subfieldseries: error: cannot handle different time scales yet',nper1,nper2
        write(*,*) 'subfieldseries: error: cannot handle different time scales yet',nper1,nper2
        call exit(-1)
    end if
    if ( nperyear == 366 ) then
        write(0,*) 'subfieldseries: error: cannot handle leap years yet',nperyear
        write(*,*) 'subfieldseries: error: cannot handle leap years yet',nperyear
        call exit(-1)
    end if
    if ( mens1 /= mens ) then
        write(0,*) 'subfieldseries: error: cannot handle ensemble yet',mens1,mens
        write(*,*) 'subfieldseries: error: cannot handle ensemble yet',mens1,mens
        call exit(-1)
    end if
    nens1 = mens1
    nens2 = mens

!   init

    call getarg(3,outfile)
    inquire(file=outfile,exist=lexist)
    if ( lexist ) then
        open(1,file=outfile)
        close(1,status='delete')
    end if

!   read fields

    if ( ncid == -1 ) then
        call readdatfile(datfile,field,nx,ny,nx,ny,nperyear &
            ,firstyr,lastyr,firstyr,firstmo,nt,undef,endian,lwrite &
            ,firstyr,lastyr,1,1)
    else
        call readncfile(ncid,field,nx,ny,nx,ny,nperyear &
            ,firstyr,lastyr,firstyr,firstmo,nt,undef,lwrite,firstyr &
            ,lastyr,ivars1)
    end if

!   subtract field and series (no regression, just straight subtraction)

    do i=firstyr,lastyr
        do j=1,nperyear
            do jy=1,ny
                do jx=1,nx
                    if ( field(jx,jy,j,i) < 1e33 .and. data(j,i) < 1e33 ) then
                        field(jx,jy,j,i) = field(jx,jy,j,i) - data(j,i)
                    else
                        field(jx,jy,j,i) = 3e33
                    end if
                end do
            end do
        end do
    end do

!   output

    call merge_metadata(metadata,n,metadata2,'global mean temperature',history2,'series_')
    n = n + 1
    metadata(1,n) = 'series_file'
    metadata(2,n) = trim(seriesfile)
    if ( lwrite ) then
        print *,'vars = ',trim(vars(1))
        print *,'lvars = ',trim(lvars(1))
        print *,'svars = ',trim(svars(1))
        print *,'units = ',trim(units(1))
        print *,'cell_methods(1) = ',trim(cell_methods(1))
        print *,'lz(1) = ',trim(lz(1))
        print *,'ltime = ',trim(ltime)
        print *,'title = ',trim(title)
        print *,'history = ',trim(history)
        do i=1,100
        if ( metadata(1,i) == ' ' ) exit
            print *,'metadata(1,',i,') = ',trim(metadata(1,i))
            print *,'metadata(2,',i,') = ',trim(metadata(2,i))
        end do
    end if

    i = index(outfile,'.ctl')
    if ( i == 0 ) then
        nt = nperyear*(lastyr-firstyr+1)
        ivars = 0
        call enswritenc(outfile,ncid,ntvarid,itimeaxis,ntmax,nx,xx,ny &
            ,yy,nz,zz,lz,nt,nperyear,firstyr,1,ltime,3e33,title &
            ,history,nvars,vars,ivars,lvars,svars,units,cell_methods &
            ,metadata,nens1,nens2)
        irec = 0
        do i=firstyr,lastyr
            do j=1,nperyear
                irec = irec + 1
                call writencslice(ncid,ntvarid,itimeaxis,ntmax,ivars, &
                    field(1,1,j,i),nx,ny,nz,nx,ny,nz,irec,1)
            end do
        end do
        status = nf_close(ncid)
    else
        datfile = outfile
        datfile(i:) = '.dat'
        if ( lwrite ) print *,'calling writectl'
        call writectl(outfile,datfile,nx,xx,ny,yy,1,zz &
            ,nperyear*(lastyr-firstyr+1),nperyear,firstyr,1,3e33,title, &
            nvars,vars,ivars,lvars,units)
        if ( lwrite ) print *,'writing data, recl= ',recfa4*nx*ny
        open(unit=2,file=datfile,form='unformatted',access='direct' &
            ,recl=recfa4*nx*ny,err=900)
        irec = 0
        do i=firstyr,lastyr
            do j=1,nperyear
                irec = irec + 1
                write(2,rec=irec) ((field(jx,jy,j,i),jx=1,nx),jy=1,ny)
            end do
        end do
        close(2)
    end if
    goto 999
900 write(0,*) 'error opening data outfile'
    call exit(-1)
999 continue
end program subfieldseries
