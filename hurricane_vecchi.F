        program hurricane_vecchi
!
!       compute the expected number of hurricane (\lambda in a Poisson
!       distribution) from MDR and tropical SST
!       From Vecchi et al, 2010, MWR
!
        use lsdata
        implicit none
#include "params.h"
#include "getopts.inc"
        integer nvarmax
        parameter(nvarmax=1)
        integer yr,mo,nperyear,mens1,mens,ncid,nx,ny,nz,nt,firstyr,nvars
     +       ,firstmo,fyr,lyr,iarg,iregion
        integer ix,iy,ixx,x1,x2,y1,y2
        integer endian,jvars(6,nvarmax)
        integer,allocatable :: nn(:,:,:)
        real trop1,mdr1,s
        real xx(nxmax),yy(nymax),zz(nzmax),undef,w,ww,xxls(nxmax)
     +       ,yyls(nymax)
        real,allocatable :: field(:,:,:,:),
     +       mdr(:,:),trop(:,:),lambda(:),mean(:,:,:),wx(:),wy(:)
        logical xwrap
        character file*255,datfile*255,
     +       vars(nvarmax)*40,lvars(nvarmax)*80
     +       ,svars(nvarmax)*80,units(nvarmax)*40,lz(3)*20,ltime*120
     +       ,title*255,history*2048,cell_methods(nvarmax)*100
     +       ,metadata(2,100)*2000
        character lsmasktype*4
        integer iargc
!
        if ( iargc().lt.2 ) then
            write(0,*) 'usage: hurricane_vecchi field lsmask [options]'
            call abort
        end if
        lstandardunits = .true.
        call getarg(1,file)
        nperyear = 12
        iarg = 2
        call getlsmask(iarg,lsmasktype,nxmax,xxls,nymax,yyls,lwrite)
        call getmetadata(file,mens1,mens,ncid,datfile,nxmax,nx
     +       ,xx,nymax,ny,yy,nzmax,nz,zz,lz,nt,nperyear,firstyr,firstmo
     +       ,ltime,undef,endian,title,history,nvarmax,nvars,vars,jvars
     +       ,lvars,svars,units,cell_methods,metadata,lwrite)
        call getopts(iarg,iargc(),nperyear,yrbeg,yrend,.false.,mens1
     +       ,mens)
        if ( nperyear.ne.12 ) then
            write(0,*) 'error: can only handle month;ly data, not '
     +           ,nperyear
            call abort
        end if
        if ( lwrite ) print *,'calling getlsmask'
        if ( lsmasktype.ne.'all' ) then
            call checkgridequal(nx,ny,xx,yy,nxls,nyls,xxls,yyls)
        endif
        lyr = firstyr + (nt+firstmo-1)/nperyear
        lyr = min(yr2,lyr)
        fyr = max(yr1,firstyr)
        yr1 = max(yr1,fyr)
        yr2 = min(yr2,lyr)
        if ( lwrite ) print *,'allocating field(',nx,ny,nperyear,fyr,lyr
     +       ,nens2,')'
!
!       read data
!
        allocate(field(nx,ny,nperyear,fyr:lyr))
        allocate(nn(nx,ny,nperyear))
        allocate(mean(nx,ny,nperyear))
        allocate(wx(nx))
        allocate(wy(ny))
        allocate(mdr(nperyear,fyr:lyr))
        allocate(trop(nperyear,fyr:lyr))
        allocate(lambda(fyr:lyr))
        call readfield(ncid,file,datfile,field,nx,ny,nz
     +       ,nperyear,fyr,lyr,nens1,nens2,nx,ny,nz,nperyear,yr1,yr2
     +       ,firstyr,firstmo,nt,undef,endian,vars,units,lstandardunits
     +       ,lwrite)
*
*       apply land/sea mask
*
        call applylsmask(field,lsmask,nx,ny,nz,nperyear,fyr,lyr,
     +       nens1,nens2,lsmasktype,lwrite)
!
!       make anomalies
!
        call getmean(mean,nn,nx,ny,nperyear,field,nx,ny
     +       ,nperyear,fyr,lyr,nx,ny,1982,1,nperyear*(2005-1982+1)
     +       ,lwrite)
        do yr=fyr,lyr
            do mo=1,nperyear
                do iy=1,ny
                    do ix=1,nx
                        if ( field(ix,iy,mo,yr).lt.1e33 .and.
     +                       mean(ix,iy,mo) .lt.1e33 ) then
                            field(ix,iy,mo,yr) =
     +                           field(ix,iy,mo,yr) - mean(ix,iy,mo)
                        else
                            field(ix,iy,mo,yr) = 3e33
                        end if
                    end do
                end do
            end do
        end do
!
!       cut out MDR and Trop SST averages
!
        call getweights('x',xx,wx,nx,xwrap,lwrite)
        call getweights('y',yy,wy,ny,.FALSE.,lwrite)
        do iregion=1,2
            if ( iregion.eq.1 ) then
                lat1 = 10
                lat2 = 25
                lon1 = -80
                lon2 = -20
            else
                lat1 = -30
                lat2 =  30
                lon1 =   0
                lon2 = 360
            end if
            call getlatlonwindow(lat1,lat2,lon1,lon2,
     +           xx,nx,xwrap,avex,yy,ny,avey,x1,x2,y1,y2,lwrite)
            do yr=fyr,lyr
                do mo=1,nperyear
                    s = 0
                    w = 0
                    ww = 0
                    do iy = y1,y2
                        do ixx=x1,x2
                            ix = ixx
                            if ( ixx.lt.1 ) ix = ix + nx
                            if ( ixx.gt.nx ) ix = ix - nx
                            if ( field(ix,iy,mo,yr).lt.1e33 ) then
                                s = s + wx(ix)*wy(iy)*field(ix,iy,mo,yr)
                                w = w + wx(ix)*wy(iy)
                            end if
                            ww = ww + wx(ix)*wy(iy)
                        end do
                    end do
                    if ( iregion.eq.1 ) then
                        if ( w.gt.ww/2 ) then
                            mdr(mo,yr) = s/w
                        else
                            mdr(mo,yr) = 3e33
                        end if
                    else
                        if ( w.gt.ww/2 ) then
                            trop(mo,yr) = s/w
                        else
                            trop(mo,yr) = 3e33
                        end if
                    end if
                end do          ! yr
            end do              ! mo
        end do                  ! iregion
!       
!       hurricane model
!
        do yr=fyr,lyr
            mdr1 = (mdr(8,yr) + mdr(9,yr) + mdr(10,yr))/3
            trop1 = (trop(8,yr) + trop(9,yr) + trop(10,yr))/3
            if ( mdr1.lt.1e29 .and. trop1.lt.1e29 ) then
                lambda(yr) = exp(1.707 + 1.388*mdr1 - 1.521*trop1)
            else
                lambda(yr) = 3e33
            end if
        end do
!
!       output
!
        print '(a)','# Number of hurricanes over the Atlantic '//
     +       'computed from '//trim(title)
        print '(a)','# using the statistical hurricane model of '//
     +       'Vecchi et al, 2011, Monthly Weather Review'
        print '(a)','# hurr [1] expected seasonal hurricane count'
        call printdatfile(6,lambda,1,1,fyr,lyr)
!
        end
