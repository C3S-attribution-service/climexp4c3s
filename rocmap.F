        program rocmap
*
*       compute a map of the area under the ROC curve
*
*       input: netcdf files of obervations and forecasts
*       output: netcdf file of area under the ROC curve
*
*       GJvO 2006
*
        implicit none
#include "netcdf.inc"
#include "params.h"
        integer ntmax
        parameter(ntmax=10000)
        integer i,j,k,iens,ncid,nx,ny,nz,nt,nperyear,firstyr,firstmo
     +       ,nens1,nens2,nvars,jvars(6,1),status,itimeaxis(1),yrbegin
     +       ,ntvarid,mobegin,ivars(2,1),it
        real thobs,pthobs,thmod,pthmod
        real xx(nxmax),yy(nymax),zz(nzmax),undef
        real,allocatable :: obs(:,:,:),fcst(:,:,:,:),
     +       oneobs(:),onefcst(:,:),area(:,:)
        logical lwrite,lsame,tdefined(ntmax)
        character file*255,title*255,history*2000,vars(1)*20,
     +       lvars(1)*120,svars(1)*120,units(1)*60,lz(3)*20,ltime*100,
     +       cell_methods(1)*128,metadata(2,100)*2000,metadata2(2,100)*2000
        integer iargc

        lwrite = .false.
        if ( iargc().lt.5 ) then
            print * ,'Usage: threshold_obs threshold_model '//
     +           'obs.nc fcst.nc out.nc [debug]'
            stop
        endif
        if ( iargc().ge.6 ) then
            call getarg(6,file)
            if ( file.eq.'debug' .or. file.eq.'lwrite' ) then
                lwrite =.true.
            endif
        endif
*
*       read thresholds
*
        call getarg(1,file)
        call readthreshold(file,thobs,pthobs)
        if ( lwrite ) print *,'read obs threshold ',thobs,pthobs
        call getarg(2,file)
        lsame = .false.
        if ( file(1:4).eq.'same' ) then
            if ( thobs.lt.1e30 ) then
                thmod = thobs
                pthmod = 3e33
            else
                lsame = .true.
            endif
        else
            call readthreshold(file,thmod,pthmod)
        endif
        if ( lwrite ) print *,'read mod threshold ',thobs,pthobs
*
*       read observations
*
        call getarg(3,file)
        ncid = 0
        if ( lwrite ) print *,'opening file ',trim(file)
        call ensparsenc(file,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax
     +       ,nz,zz,lz,nt,nperyear,firstyr,firstmo,ltime,tdefined,ntmax
     +       ,nens1,nens2,undef,title,history,1,nvars,vars,jvars,lvars
     +       ,svars,units,cell_methods,metadata)
        if ( nz.gt.1 ) then
            write(0,*) 'rocmap: error: cannot handle 3D fields yet'
            call abort
        endif
        do it=1,nt
            if ( .not.tdefined(it) ) then
                write(0,*) 'rocmap: error: cannot handle holes in time '
     +               ,'axis yet'
                call abort
            end if
        end do
        allocate(obs(nx,ny,nt))
        status = nf_get_var_real(ncid,jvars(1,1),obs)
        if ( status.ne.nf_noerr ) then
            write(0,*) 'rocmap: error reading observations from '
     +           ,trim(file)
            call abort
        endif
*
*       read forecasts, these are assumed to be on the same axes.
*
        call getarg(4,file)
        ncid = 0
        if ( lwrite ) print *,'opening file ',trim(file)
        call ensparsenc(file,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax
     +       ,nz,zz,lz,nt,nperyear,firstyr,firstmo,ltime,tdefined,ntmax
     +       ,nens1,nens2,undef,title,history,1,nvars,vars,jvars,lvars
     +       ,svars,units,cell_methods,metadata2)
        if ( lwrite ) print *,'allocating forecast array'
        allocate(fcst(nx,ny,nt,nens1:nens2))
        status = nf_get_var_real(ncid,jvars(1,1),fcst)
        if ( lwrite ) print *,'reading forecasts'
        if ( status.ne.nf_noerr ) then
            write(0,*) 'rocmap: error reading observations from '
     +           ,trim(file)
            call abort
        endif
*
*       and compute the area under the roc curve
*
        if ( lwrite ) print *,'allocating arrays'
        allocate(oneobs(nt))
        allocate(onefcst(nens2+1,nt))
        allocate(area(nx,ny))
        do j=1,ny
            do i=1,nx
                do k=1,nt
                    oneobs(k) = obs(i,j,k)
                    do iens=nens1,nens2
                        onefcst(iens+1,k) = fcst(i,j,k,iens)
                    enddo
                enddo
                call probroc(oneobs,onefcst,nt,nens2+1,thobs,pthobs,
     +               thmod,pthmod,lsame,.false.,lwrite,area(i,j))
            enddo
        enddo
*
*       and write to netcdf file
*
        call getarg(5,file)
        nt = 1
        ntvarid = 0
        itimeaxis(1) = 0
        yrbegin = 1
        mobegin = 1
        undef = 3e33
        nvars = 1
        vars(1) = 'roc'
        lvars(1) = 'area under the roc curve'
        units(1) = ''
        ivars(1,1) = 0
        call writenc(file,ncid,ntvarid,itimeaxis,ntmax,nx,xx,ny,yy
     +       ,nz,zz,nt,nperyear,yrbegin,mobegin,undef,title,nvars
     +       ,vars,ivars,lvars,units,0,0)
        call writencslice(ncid,ntvarid,itimeaxis,ntmax,ivars,area
     +        ,nx,ny,nz,nx,ny,nz,1,1)
        i = nf_close(ncid)
        end
