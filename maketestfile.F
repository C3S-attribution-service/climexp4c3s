        program maketestfile
*
*       create a test .ctl/.dat or .nc file to test the EOF routine
*
        implicit none
#include "recfac.h"
#include "params.h"
        integer i,j,k,nx,ny,nz,nt,ivars(2,1),ncid,istart(4),icount(4)
     +        ,ntvarid,itimeaxis(ntmax)
        real xx(nxmax),yy(nxmax),zz(nzmax),x(nymax),y(nymax),
     +        data(nxmax,nymax)
        character vars(1)*10,lvars(1)*40,units(1)*10,file*256,
     +        datfile*256
        real*8 ranf
        external ranf
        call getarg(1,file)
        if ( file.eq.' ' ) then
            print *,'usage: matetestfile name[.ctl|.nc]'
            call abort
        elseif ( index(file,'.nc').eq.0 .and.  index(file,'.ctl').eq.0 )
     +            then
            print *,'usage: matetestfile name[.ctl|.nc]'
            call abort
        endif
*
        nx = 3
        ny = 3
        nz = 1
        nt = 1200
        do i=1,nx
            xx(i) = i
        enddo
        do i=1,ny
            yy(i) = i
        enddo
        zz(1) = 0
        vars(1) = 'test'
        lvars(1) = 'test data for EOF'
        ivars(1,1) = 1
        ivars(2,1) = 99
        if ( index(file,'.ctl').ne.0 ) then
            datfile = file
            i = index(file,'.ctl')
            datfile(i:) = '.dat'
            call writectl(file,datfile,nx,xx,ny,yy,nz,zz,nt,1900,1,3e33,
     +            'test file for EOF ',1,vars,ivars,lvars,units)
            open(1,file='test.dat',form='unformatted',access='direct'
     +            ,recl=nx*ny*recfa4)
        else
            call writenc(file,ncid,ntvarid,itimeaxis,ntmax,nx,xx,ny,yy
     +            ,nz,zz,nt,12,1900,1,3e33,'test file for EOF ',1,vars
     +            ,ivars,lvars,units,0,0)
        endif
        do i=1,nt
            do j=1,ny
                x(j) = -1+2*ranf(2*i-1+nt*j)
                y(j) = -1+2*ranf(2*i  +nt*j)
                data(1,j) = -x(1)
                data(2,j) = x(1)+y(1)
                data(3,j) = -y(1)
            enddo
            if ( index(file,'.ctl').ne.0 ) then
                write(1,rec=i) ((data(k,j),k=1,3),j=1,ny)
            else
                print *, ((data(k,j),k=1,3),j=1,ny)
                call writencslice(ncid,ivars(1,1),data,nxmax,nymax,1,nx
     +                ,ny,nz,i,1)
            endif
        enddo
        if ( index(file,'.ctl').ne.0 ) then
            close(1)
        else
            i = nf_close(ncid)
        endif
*
        end

