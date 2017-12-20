program yearly2shorterfield
!
!   convert a yearly time series to a monthly one
!
    implicit none
    integer,parameter :: nvarmax=1
    include 'params.h'
    include 'getopts.inc'
    integer yr,mo,dy,nperyear,npernew,k,i,nfac,nxf,nyf,nzf,fyr,lyr,iens,ix,iy,iz
    integer mens1,mens,ncid,nx,ny,nz,nt,firstyr,firstmo,endian,yrbegin,mobegin
    integer nvars,ivars(2,nvarmax),jvars(6,nvarmax),ntvarid,ntmax
    integer,allocatable :: itimeaxis(:)
    real xx(nxmax),yy(nymax),zz(nzmax),undef
    real,allocatable :: data(:,:,:,:,:,:),newdata(:,:,:,:,:,:),fxy(:,:),newfxy(:,:)
    character file*1024,datfile*1024,vars(1)*40,svars(1)*100,lvars(1)*120,units(1)*60
    character lz(3)*10,ltime*100,title*1000,history*20000,cell_methods(nvarmax)*100
    character metadata(2,100)*2000
    character string1*50,string2*50
    integer iargc

    call killfile(file,string1,string2,0) ! random strings
    if ( iargc().lt.2 ) then
        print *,'usage: yearly2shorterfield infile.nc nperyearnew [mon n ave|sum m]'
        stop
    endif
    call getarg(1,file)
    call getarg(2,string1)
    read(string1,*,err=901) npernew

    call getmetadata(file,mens1,mens,ncid,datfile,nxmax,nx &
     &       ,xx,nymax,ny,yy,nzmax,nz,zz,lz,nt,nperyear,firstyr,firstmo &
     &       ,ltime,undef,endian,title,history,nvarmax,nvars,vars,jvars &
     &       ,lvars,svars,units,cell_methods,metadata,lwrite)
    call getopts(3,iargc()-1,nperyear,yrbeg,yrend,.false.,0,nensmax)
    if ( iargc().gt.2 ) then
        if ( oper.eq.'v' ) then
            nfac = 1
        else if ( oper.eq.'+' ) then
            nfac = lsum
        else
            goto 904
        endif
    else
        m1 = 1
        lsum = npernew
        nfac = 1
    endif
    
    nxf = nx
    nyf = ny
    nzf = nz
    fyr = firstyr
    lyr = firstyr + (nt-1)/nperyear
    yr1 = max(yr1,fyr)
    yr2 = min(lyr,lyr)
    nens1 = max(nens1,mens1)
    nens2 = min(nens2,mens)
    allocate(data(nxf,nyf,nzf,nperyear,fyr:lyr,0:mens))
    allocate(newdata(nxf,nyf,nzf,npernew,fyr:lyr,0:mens))
    call keepalive1('Reading data',1,1)
    call readfield(ncid,file,datfile,data,nxf,nyf,nzf &
     &       ,nperyear,fyr,lyr,nens1,nens2,nx,ny,nz,nperyear,yr1,yr2 &
     &       ,firstyr,firstmo,nt,undef,endian,vars,units,lstandardunits &
     &       ,lwrite)
    
    allocate(fxy(nperyear,fyr:lyr))
    allocate(newfxy(npernew,fyr:lyr))
    newdata = 3e33 ! safety first
    i = 0 
    do iens=nens1,nens2
        do iz=1,nz
            do iy=1,ny
                do ix=1,nx
                    i = i + 1
                    call keepalive('Transforming',i,nx*ny*nz*(nens2-nens1+1))
                    do mo=1,nperyear
                        do yr=fyr,lyr
                            fxy(mo,yr) = data(ix,iy,iz,mo,yr,iens)
                        end do
                    end do
                    call annual2shorter(fxy,nperyear,fyr,lyr,nperyear, &
                    &       newfxy,npernew,fyr,lyr,npernew,m1,lsum,nfac,lwrite)
                    do mo=1,npernew
                        do yr=fyr,lyr
                            newdata(ix,iy,iz,mo,yr,iens) = newfxy(mo,yr)
                        end do
                    end do
                end do
            end do
        end do
    end do
    call getarg(iargc(),file) ! outfile
    yrbegin = fyr
    mobegin = 1
    undef = 3e33
    nt = npernew*(lyr-fyr+1)
    ntmax = nt
    allocate(itimeaxis(ntmax))
    call enswritenc(file,ncid,ntvarid,itimeaxis,ntmax,nx,xx,ny              &
     &       ,yy,nz,zz,lz,nt,npernew,yrbegin,mobegin,ltime,undef,title      &
     &       ,history,nvars,vars,ivars,lvars,svars,units,cell_methods       &
     &       ,metadata,nens1,nens2)
    do iens=nens1,nens2
        do yr=fyr,lyr
            do mo=1,npernew
                call writencslice(ncid,ntvarid,itimeaxis,nt,ivars,          &
                &   newdata(1,1,1,mo,yr,iens),nx,ny,nz,nx,ny,nz,mo+npernew*(yr-fyr),1+iens)
            end do
        end do
    end do
    call nf_close(ncid) ! do not forget
    goto 999
901 write(0,*) 'yearly2shorterfield: error reading npernew from ',trim(file)
    call exit(-1)
904 write(0,*) 'yearly2shorterfield: error: expecting ''ave|sum'', not ',oper
    call exit(-1)
999 continue
end program
