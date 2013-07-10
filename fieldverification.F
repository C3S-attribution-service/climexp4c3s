        program fieldverification
*
*       prepares fields for analysis with the R routines.
*       - select months, seasons etc
*       - detrend, debias
*       - regrid
*
        implicit none
#include "params.h"
#include "netcdf.inc"
#include "getopts.inc"
        integer yr,mo,i,jx,jy,mens,mens1
        integer nxobs,nyobs,nzobs,ntobs,nperyearobs,firstmoobs
     +       ,firstyrobs,endianobs,ivarobs(2,1),jvarobs(6,1),ncidobs,
     +       nx,ny,nz,nt,nvars,nperyear,firstmo,firstyr,endian
     +       ,ivars(2,1),jvars(6,1),ncid,status
        real xxobs(nxmax),yyobs(nymax),zzobs(nzmax),undefobs,
     +       xx(nxmax),yy(nymax),zz(nzmax),undef
        real,allocatable :: obs(:,:,:,:),fcst(:,:,:,:,:)
        character obsfile*255,obsdat*255,ensfile*255,ensdat*255
        character titleobs*255,varsobs(1)*20,lvarsobs(1)*126,
     +       unitsobs(1)*40
        character title*255,vars(1)*20,lvars(1)*128,units(1)*40
        integer iargc
*
        if ( iargc().lt.4 ) then
            print *,'usage: fieldverification '//
     +           'obsin obsout fcstin fcstout '//
     +           '[mon n] [sum m] [detrend] '//
     +           '[debias mean|var|all]'
            stop
        endif
*
        call getarg(1,obsfile)
        if ( lwrite ) print *,'fieldverification: nf_opening file '
     +       ,trim(obsfile)
        status = nf_open(infile,nf_nowrite,ncid)
        if ( status.ne.nf_noerr ) then
            call parsectl(obsfile,obsdat,nxmax,nxobs,xxobs,nymax,nyobs
     +           ,yyobs,nzmax,nzobs,zzobs,ntobs,nperyearobs,firstyrobs
     +           ,firstmoobs,undefobs,endianobs,titleobs,1,nvars,varsobs
     +           ,ivarsobs,lvarsobs,unitsobs)
            ncid = -1
        else
            call parsenc(obsfile,ncidobs,nxmax,nxobs,xxobs,nymax,nyobs
     +           ,yyobs,nzmax,nzobs,zzobs,ntobs,nperyearobs,firstyrobs
     +           ,firstmoobs,undefobs,titleobs,1,nvars,varsobs,jvarsobs
     +           ,lvarsobs,unitsobs)
        endif
        call getarg(3,fcstfile)
        if ( index(fcstfile,'%').gt.0 .or. index(fcstfile,'++').gt.0 )
     +       then
            ensemble = .true.
            call filloutens(fcstfile,0)
            inquire(file=infile,exist=lexist)
            if ( .not.lexist ) then
                mens1 = 1
                call filloutens(fcstfile,1)
            else
                mens1 = 0
            endif
        else
            ensemble = .false.
            mens1 = 0
            mens = 0
        endif
        if ( lwrite ) print *,'fieldverification: nf_opening file '
     +       ,trim(fcstfile)
        status = nf_open(infile,nf_nowrite,ncid)
        if ( status.ne.nf_noerr ) then
            call parsectl(fcstfile,fcstdat,nxmax,nx,xx,nymax,ny,yy,nzmax
     +           ,nz,zz,nt,nperyear,firstyr,firstmo,undef,endian,title,1
     +           ,nvars,vars,ivars,lvars,units)
            ncid = -1
        else
            call parsenc(fcstfile,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax,nz
     +           ,zz,nt,nperyear,firstyr,firstmo,undef,title,1,nvars
     +           ,vars,jvars,lvars,units)
        endif
        if ( ensemble ) then
            do mens=mens1+1,nensmax
                call getarg(1,line)
                call filloutens(line,mens)
                inquire(file=line,exist=lexist)
                if ( .not.lexist ) goto 100
            enddo
 100        continue
            mens = mens - 1
            write(0,*) 'located ',mens-mens1+1
     +           ,' ensemble members<br>'
        endif
        
