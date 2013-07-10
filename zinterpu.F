        subroutine zinterpu(field1,zz1,nx1,ny1,nz1,field2,zz2,nx2,ny2,
     +        nz2,zz,nz,yr11,yr21,yr12,yr22,nxf,nyf,nzf,nperyear,lwrite)
*
*       interpolate the z-coordinates of field1 and field2 to a common
*       grid.
*       NOTE.  Stopgap version that can only do extension of a 2D grid
*       to the other grid.
*
        implicit none
#include "params.h"
        integer nx1,ny1,nz1,nx2,ny2,nz2,yr11,yr21,yr12,yr22,nxf,nyf
     +        ,nzf,nz,nperyear
        real zz(nzmax),
     +        field1(nxf,nyf,nzf,nperyear,yr11:yr21),zz1(nz1),
     +        field2(nxf,nyf,nzf,nperyear,yr12:yr22),zz2(nz2)
        logical lwrite
        integer i,j,k,mo,yr
        if ( lwrite ) then
            print *,'zinterpu: input parameters:'
            print *,'field1: ',nx1,ny1,nz1,yr11,yr21
            print *,'field2: ',nx2,ny2,nz2,yr12,yr22
            print *,'n[xyz]f:',nxf,nyf,nzf
            print *,'nperyear',nperyear
        endif
*       
*       if we just have to extend a 2D field into 3D, do it.
*
        if ( nz1.eq.1 ) then
            nz = nz2
            do k=1,nz2
                zz(k) = zz2(k)
            enddo
            if ( nz2.eq.1 ) then
                if ( lwrite ) print *,'zinterpu: both fields are 2D'
                return
            endif
            if ( lwrite ) print *
     +            ,'zinterpu: extending field1 from 2D to 3D'
            do yr=yr11,yr21
                do mo=1,nperyear        
                    do k=2,nz2
                        do j=1,ny1
                            do i=1,nx1
                                field1(i,j,k,mo,yr)=field1(i,j,1,mo,yr)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            return
        elseif ( nz2.eq.1 ) then
            if ( lwrite ) print *
     +            ,'zinterpu: extending field2 from 2D to 3D'
            nz = nz1
            do k=1,nz1
                zz(k) = zz1(k)
            enddo
            do yr=yr12,yr22
                do mo=1,nperyear        
                    do k=2,nz1
                        do j=1,ny2
                            do i=1,nx2
                                field2(i,j,k,mo,yr)=field2(i,j,1,mo,yr)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            return
        endif
*       
*       check if there is anything to do
*       
        if ( nz1.ne.nz2 ) goto 100
        do k=1,nz1
            if ( zz1(k).ne.zz2(k) ) goto 100
        enddo
*       identical grids - do nothing
        if ( lwrite ) print *,'zinterpu: identical depth grids'
        nz = nz1
        do k=1,nz
            zz(k) = zz1(k)
        enddo
        return
  100   continue
*       
*       the hard work.
*       
        write(0,*) 'zinterpu: cannot interpolate depths yet.'
        write(*,*) 'zinterpu: cannot interpolate depths yet.'
        call abort
        end
