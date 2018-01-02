program patchfield
!
!       fill in missing values in a field or extend it in time using
!       data from a second field. A linear regression on the overlap is
!       used to adjust the second field to the first one at each grid point.
!
    implicit none
    include 'params.h'
    include 'getopts.inc'
    include 'netcdf.inc'
    integer nvarmax,ntmax
    parameter (nvarmax=1,ntmax=500*12)
    integer mens,mens1,ncid,nx,ny,nz,nt,nperyear,firstyr,firstmo,endian, &
 &      nvars,ivars(6,nvarmax)
    integer mens01,mens11,ncid1,nx1,ny1,nz1,nt1,nperyear1,firstyr1,firstmo1,endian1, &
 &      nvars1,ivars1(6,nvarmax)
    integer ix,iy,iz,i,j,dy,mo,yr,n,nperday,dpm(12),fyr,lyr,status,itimeaxis(ndata),ntvarid
    real xi(30*(yrend-yrbeg+1)),yi(30*(yrend-yrbeg+1))
    real scale(12),offset(12),a,b,siga,sigb,chi2,q,sig(1),sx(12),sy(12)
    real xx(nxmax),yy(nymax),zz(nzmax),undef
    real xx1(nxmax),yy1(nymax),zz1(nzmax),undef1
    real,allocatable :: mainfield(:,:,:,:,:),auxfield(:,:,:,:,:),maindata(:,:),auxdata(:,:)
    logical lreversefit
    character &
        mainfile*255,auxfile*255,outfile*255,datfile*255,datfile1*255,lz(3)*20,ltime*120, &
        title*2000,title1*2000,history*20000,vars(nvarmax)*40,lvars(nvarmax)*120,svars(nvarmax)*80, &
        units(nvarmax)*40,units1(nvarmax)*40,cell_methods(nvarmax)*40,metadata(2,100)*2000, &
        metadata1(2,100)*2000
    character method*4
    integer iargc
    data dpm /31,29,31,20,31,30,31,31,30,31,30,31/
!
!   init
!
    lwrite = .false.
    lreversefit = .false.
    if ( iargc().lt.3 ) then
        write(0,*) 'usage: patchfield mainfield auxfield [regr|bias|none] outfield'
        write(0,*) 'patches holes in mainfield using data from '// &
 &           'auxfield using a linear regression, bias correction or straight copy '// &
 &           'for each grid point'
        call exit(-1)
    endif
!
!   read metadata
!
    call getarg(2,auxfile)
    call getmetadata(auxfile,mens11,mens01,ncid1,datfile1,nxmax,nx1 &
 &       ,xx1,nymax,ny1,yy1,nzmax,nz1,zz1,lz,nt1,nperyear1,firstyr1,firstmo1 &
 &       ,ltime,undef1,endian1,title1,history,nvarmax,nvars,vars,ivars1 &
 &       ,lvars,svars,units1,cell_methods,metadata1,lwrite)
    call getarg(1,mainfile)
    call getmetadata(mainfile,mens1,mens,ncid,datfile,nxmax,nx &
 &       ,xx,nymax,ny,yy,nzmax,nz,zz,lz,nt,nperyear,firstyr,firstmo &
 &       ,ltime,undef,endian,title,history,nvarmax,nvars,vars,ivars &
 &       ,lvars,svars,units,cell_methods,metadata,lwrite)
    if ( iargc() >= 4 ) then
        call getarg(3,method)
        ! "noscale" was recognised by a version that I mistakenly skipped
        if ( method == 'nosc' ) method = 'bias'
    else
        method = 'none'
    end if
!
!   check grids are equal
!
    if ( nperyear /= nperyear1 ) then
        write(0,*) 'error: nperyear unequal: ',nperyear,nperyear1
        call exit(-1)
    end if
    if ( mens /= 0 .or. mens1 /= 0 .or. mens01 /= 0 .or. mens11 /= 0 ) then
        write(0,*) 'error: cannot handle ensembles yet'
        call exit(-1)
    end if
    call checkgridequal3d(nx,ny,nz,xx,yy,zz,nx1,ny1,nz1,xx1,yy1,zz1)
!
!   allocate fields
!
    fyr = min(firstyr1,firstyr)
    lyr = max(firstyr1 + (firstmo1+nt1-2)/nperyear1, firstyr + (firstmo+nt-2)/nperyear)
    yr1 = fyr
    yr2 = lyr
    if ( lwrite ) print *,'allocating arrays ',nx,ny,nz,nperyear,fyr,lyr
    allocate(mainfield(nx,ny,nz,nperyear,fyr:lyr),maindata(nperyear,fyr:lyr))
    allocate(auxfield(nx,ny,nz,nperyear,fyr:lyr),auxdata(nperyear,fyr:lyr))
!
!   read data
!
    write(0,*) 'Reading ',trim(mainfile)
    call readfield(ncid,mainfile,datfile,mainfield,nx,ny &
 &       ,nz,nperyear,fyr,lyr,0,0,nx,ny,nz,nperyear,fyr,lyr &
 &       ,firstyr,firstmo,nt,undef,endian,vars,units,.true. &
 &       ,lwrite)
    write(0,*) 'Reading ',trim(auxfile)
    call readfield(ncid1,auxfile,datfile1,auxfield,nx,ny &
 &       ,nz,nperyear,fyr,lyr,0,0,nx,ny,nz,nperyear,fyr,lyr &
 &       ,firstyr1,firstmo1,nt1,undef1,endian1,vars,units1,.true. &
 &       ,lwrite)
!
!   loop over all grid points
!
    do iz=1,nz
        do iy=1,ny
            do ix=1,nx
                call keepalive1('point ',ix+(iy-1)*nx+(iz-1)*nx*ny,nx*ny*nz)
                do yr=fyr,lyr
                    do mo=1,nperyear
                        maindata(mo,yr) = mainfield(ix,iy,iz,mo,yr)
                        auxdata(mo,yr) = auxfield(ix,iy,iz,mo,yr)
                    end do
                end do
!
!               determine regression coefficents per month or less
!
                if ( nperyear.le.12 ) then
                    sx = 0
                    sy = 0
                    do mo=1,nperyear
                        n = 0
                        do yr=fyr,lyr
                            if ( maindata(mo,yr).lt.1e33 .and. &
             &                   auxdata(mo,yr).lt.1e33 ) then
                                n = n + 1
                                xi(n) = auxdata(mo,yr)
                                sx(mo) = sx(mo) + xi(n)
                                yi(n) = maindata(mo,yr)
                                sy(mo) = sy(mo) + yi(n)
                            end if
                        end do
                        if ( n.gt.0 ) then
                            sx(mo) = sx(mo)/n
                            sy(mo) = sy(mo)/n
                        end if
                        if ( method == 'none' ) then
                            scale(mo) = 1
                            offset(mo) = 0
                        else if ( method == 'bias' ) then
                            if ( n.lt.7 ) then
                                scale(mo) = 3e33
                                offset(mo) = 3e33
                            else
                                scale(mo) = 1
                                offset(mo) = sy(mo) - sx(mo)
                            end if
                        else if ( method == 'regr' ) then
                            if ( n.lt.7 ) then ! arbitrary
                                scale(mo) = 3e33
                                offset(mo) = 3e33
                            else if ( lreversefit ) then
                                call fit(xi,yi,n,sig,0,offset(mo),scale(mo),siga &
                 &                   ,sigb,chi2,q)
                            else
                                call fit(yi,xi,n,sig,0,a,b,siga,sigb,chi2,q)
                                scale(mo) = 1/b
                                offset(mo) = -a/b
                                if ( scale(mo).lt.0.67 .or. scale(mo).gt.1.5 ) then
                                    if ( lwrite ) write(0,*) 'warning: scale(',mo,') = ', &
                 &                      scale(mo),xi(ix),yi(iy),zz(iz)
                                    scale(mo) = 3e33
                                    maindata = 3e33 ! otherwise we get a discontinuity
                                end if
                            end if
                        end if
                    end do
                    do yr=fyr,lyr
                        do mo=1,nperyear
                            if ( maindata(mo,yr).gt.1e33 .and. &
             &                   auxdata(mo,yr).lt.1e33 .and. &
             &                   scale(mo).lt.1e33 ) then
                                maindata(mo,yr) = &
             &                       offset(mo) + scale(mo)*auxdata(mo,yr)
                            end if
                        end do
                    end do
                else if ( nperyear.ge.260 ) then ! daily or 6-hourly frequency
                    if ( nperyear.eq.360 .or. nperyear.eq.4*360 ) then
                        dpm = 30
                        nperday = nperyear/360
                    else if ( nperyear.eq.365 .or. nperyear.eq.4*365 ) then
                        dpm(2) = 28
                        nperday = nperyear/365
                    else
                        nperday = nint(nperyear/366.)
                    end if
                    sx = 0
                    sy = 0
                    do mo=1,12
                        n = 0
                        do yr=fyr,lyr
                            do dy=1,nperday*dpm(mo)
                                call invgetdymo(dy,mo,j,nperyear)
                                if ( maindata(j,yr).lt.1e33 .and. &
             &                       auxdata(j,yr).lt.1e33 ) then
                                    n = n + 1
                                    xi(n) = auxdata(j,yr)
                                    sx(mo) = sx(mo) + xi(n)
                                    yi(n) = maindata(j,yr)
                                    sy(mo) = sy(mo) + yi(n)
                                end if
                            end do
                        end do
                        if ( n.gt.0 ) then
                            sx(mo) = sx(mo)/n
                            sy(mo) = sy(mo)/n
                        end if
                        if ( method == 'none' ) then
                            scale(mo) = 1
                            offset(mo) = 0
                        else if ( method == 'bias' ) then
                            if ( n.lt.7 ) then
                                scale(mo) = 3e33
                                offset(mo) = 3e33
                            else
                                scale(mo) = 1
                                offset(mo) = sy(mo) - sx(mo)
                            end if
                        else if ( method == 'regr' ) then
                            if ( n.lt.7 ) then ! arbitrary
                                scale(mo) = 3e33
                                offset(mo) = 3e33
                            else if ( lreversefit ) then
                                call fit(xi,yi,n,sig,0,offset(mo),scale(mo),siga &
                 &                   ,sigb,chi2,q)
                            else
                                call fit(yi,xi,n,sig,0,a,b,siga,sigb,chi2,q)
                                scale(mo) = 1/b
                                offset(mo) = -a/b
                                if ( scale(mo).lt.0.67 .or. scale(mo).gt.1.5 ) then
                                    if ( lwrite ) write(0,*) 'warning: scale(',mo,') = ', &
                 &                      scale(mo),xi(ix),yi(iy),zz(iz)
                                    scale(mo) = 3e33
                                    maindata = 3e33 ! otherwise we get a discontinuity
                                end if
                            end if
                        end if
                    end do
                    do yr=fyr,lyr
                        do mo=1,12
                            do dy=1,nperday*dpm(mo)
                                call invgetdymo(dy,mo,j,nperyear)
                                if ( maindata(j,yr).gt.1e33 .and. &
             &                       auxdata(j,yr).lt.1e33 .and. &
             &                       scale(mo).lt.1e33 ) then
                                    maindata(j,yr) = &
             &                           offset(mo) + scale(mo)*auxdata(j,yr)
                                end if
                            end do
                        end do
                    end do
                else
                   write(0,*) 'merging pentad or decadal time series not ready'
                   call abort
                end if
!
!               copy back to main field
!
                do yr=fyr,lyr
                    do mo=1,nperyear
                        mainfield(ix,iy,iz,mo,yr) = maindata(mo,yr)
                    end do
                end do
            end do
        end do
    end do
!
!   output
!
    call getarg(iargc(),outfile)
    i = index(outfile,'.ctl')
    if ( i.ne.0 ) then
        write(0,*) 'error: grads output not supported'
        call exit(-1)
    else
!       netcdf output
        undef = 3e33
        title = trim(title)//' extended with '//trim(title1)
        do i=1,1000
            if ( metadata(1,i) == ' ' ) exit
        end do
        do j=1,1000-i
            if ( metadata1(1,j) == ' ' ) exit
            metadata(1,j+i-1) = 'aux_'//trim(metadata1(1,j))
            metadata(1,j+i-1) = metadata1(2,j)
        end do
        call enswritenc(outfile,ncid,ntvarid,itimeaxis,ndata,nx,xx,ny &
 &           ,yy,nz,zz,lz,nperyear*(lyr-fyr+1),nperyear &
 &           ,fyr,1,ltime,undef,title,history,nvars,vars,ivars &
 &           ,lvars,svars,units,cell_methods,metadata,0,0)
        i = 0
        do yr=fyr,lyr
            do j=1,nperyear
                i = i + 1
                call writencslice(ncid,0,0,0,ivars, &
 &                   mainfield(1,1,1,j,yr),nx,ny,nz,nx,ny,nz,i,1)
            end do
        end do
    endif
    status = nf_close(ncid)
end program
