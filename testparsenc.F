        program testparsenc
        implicit none
#include "params.h"
        integer nvarmax
        parameter (nvarmax=10)
        integer nx,ncid,nx,ny,nz,nt,yrbegin,mobegin,nvars,
     +        ivars(6,nvarmax),i,j
        real xx(nxmax),yy(nymax),zz(nzmax),undef
        character title*255,vars(nvarmax)*10,lvars(nvarmax)*40,unit*10
        integer llen
        external llen
*
        call parsenc('aap.nc',ncid,nxmax,nx,xx,nymax,ny,yy,nzmax
     +        ,nz,zz,nt,yrbegin,mobegin,undef,title,nvarmax,nvars,vars
     +        ,ivars,lvars,units)

        print *,'ncid  = ',ncid
        print *,'title = ',title(1:llen(title))
        print *,'undef = ',undef
        print *,'nx,xx = ',nx,(xx(i),i=1,nx)
        print *,'ny,yy = ',ny,(yy(i),i=1,ny)
        print *,'nz,zz = ',nz,(zz(i),i=1,nz)
        print *,'nt,yrbegin.mobegin = ',nt,yrbegin,mobegin
        print *,'nvars = ',nvars
        do i=1,nvars
            print *,vars(i),lvars(i),(ivars(j,i),j=1,5)
        enddo

        end

