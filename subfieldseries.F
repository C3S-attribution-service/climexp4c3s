        program addfieldseries
*
*       written to subtract global mean temperature from a field
*       copied from correlatefieldfield
*
        implicit none
#include "params.h"
#include "recfac.h"
#include "netcdf.inc"
        integer n,ncid,ncid2,nx,ny,nz,nt,nper1,firstyr,firstmo,
     +       nper2,nvars,ivars(2,1),jvars1(6,nvmax),endian,endian2
     +       ,status,nperyear,lastyr
        integer i,j,jx,jy,irec
        real  xx(nxmax),yy(nymax),zz(nzmax),
     +        undef,data(12,yrbeg:yrend)
        real,allocatable :: field(:,:,:,:)
        logical lwrite,lexist
*
        character line*256,datfile*256,outfile*256,title*256,vars(1)*40,
     +       lvars(1)*80,units1(1)*40,var2*40,units2*40,yesno*1
        integer iargc,llen
*
*       check arguments
*       
        lwrite = .false.
        n = iargc()
        if ( n.ne.3 ) then
            print *,'usage: subfieldseries field.[ctl|nc] '//
     +           'series.[nc|dat] outfield.cftl'
            call abort
        endif
        call getarg(1,line)
        if ( lwrite ) print *,'correlatefieldfield: nf_opening file '
     +        ,line(1:llen(line))
        status = nf_open(line,nf_nowrite,ncid)
        if ( status.ne.nf_noerr ) then
            call parsectl(line,datfile,nxmax,nx,xx,nymax,ny,yy
     +            ,nzmax,nz,zz,nt,nper1,firstyr,firstmo,undef,endian
     +            ,title,1,nvars,vars,ivars,lvars,units1)
            ncid = -1
        else
            datfile = line
            call parsenc(line,ncid,nxmax,nx,xx,nymax,ny,yy
     +           ,nzmax,nz,zz,nt,nper1,firstyr,firstmo,undef,title,1
     +           ,nvars,vars,jvars1,lvars,units1)
        endif
        lastyr = firstyr + (firstmo+nt-2)/nper1
        if ( lwrite ) print *,'allocating field ',nx,ny,12,firstyr
     +       ,lastyr
        allocate(field(nx,ny,12,firstyr:lastyr))
*
        call getarg(2,line)
        call readseries(line,data,12,yrbeg,yrend,nper2,
     +       var2,units2,.false.,lwrite)
*
        nperyear = max(nper1,nper2)
        if ( nper1.ne.nper2 ) then
            write(0,*) 'correlatefield: error: cannot handle different'
     +            //' time scales yet',nper1,nper2
            write(*,*) 'correlatefield: error: cannot handle different'
     +            //' time scales yet',nper1,nper2
            call abort            
        endif
*
*       init
*
        call getarg(3,outfile)
        inquire(file=outfile,exist=lexist)
        if ( lexist ) then
            print *,'output file ',outfile(1:index(outfile,' ')-1),
     +            ' already exists, overwrite? [y/n]'
            read(*,'(a)') yesno
            if (  yesno.ne.'y' .and. yesno.ne.'Y' .and. 
     +            yesno.ne.'j' .and. yesno.ne.'J' ) then
                stop
            endif
            open(1,file=outfile)
            close(1,status='delete')
        endif
	print *,'init'
*
*       read fields
*
        if ( ncid.eq.-1 ) then
            call readdatfile(datfile,field,nx,ny,nx,ny,nperyear
     +           ,firstyr,lastyr,firstyr,firstmo,nt,undef,endian,lwrite
     +           ,firstyr,lastyr,1,1)
        else
            call readncfile(ncid,field,nx,ny,nx,ny,nperyear
     +           ,firstyr,lastyr,firstyr,firstmo,nt,undef,lwrite,firstyr
     +           ,lastyr,jvars1)
        endif
*
*       subtract field and series
*
        do i=firstyr,lastyr
            do j=1,12
                do jy=1,ny
                    do jx=1,nx
                        if ( field(jx,jy,j,i).lt.1e33 .and.
     +                       data(j,i).lt.1e33 ) then
                            field(jx,jy,j,i) = field(jx,jy,j,i)
     +                           - data(j,i)
                        else
                            field(jx,jy,j,i) = 3e33
                        endif
                    enddo
                enddo
            enddo
        enddo
*       
*       output
*
        datfile = outfile
        i = index(datfile,'.ctl')
        if ( i.eq.0 ) then
            write(0,*) 'error: expecting .ctl in outfile = '
     +           ,trim(outfile)
            call abort
        endif
        datfile(i:) = '.dat'
        if ( lwrite ) print *,'calling writectl'
        call writectl(outfile,datfile,nx,xx,ny,yy,1,zz
     +       ,nperyear*(lastyr-firstyr+1),nperyear,firstyr,1,3e33,title,
     +       nvars,vars,ivars,lvars,units1)
        if ( lwrite ) print *,'writing data, recl= ',recfa4*nx*ny
        open(unit=2,file=datfile,form='unformatted',access='direct'
     +            ,recl=recfa4*nx*ny,err=900)
        irec = 0
        do i=firstyr,lastyr
            do j=1,nperyear
                irec= irec + 1
                write(2,rec=irec) ((field(jx,jy,j,i),jx=1,nx),jy=1,ny)
            enddo
        enddo
        close(2)
        goto 999
 900    write(0,*) 'error opening data outfile'
        call abort
 999    continue
        end
