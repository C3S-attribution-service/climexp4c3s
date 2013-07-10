        program convert2rclimexdex
!
!       convert an (almost) arbitrary time series to the format the R
!       RClimDex require
!
        implicit none
#include "param.inc"
        integer yr1,yr2,yr,mo,dy,nperyear
        real data(npermax,yrbeg:yrend,3)
        integer ivar
        character file*255,var*40,units*40
        logical lwrite
        lwrite = .false.

        data = 3e33
        do ivar=1,3
            call getarg(ivar,file)
            if ( file.eq.' ' ) then
                write(0,*) 'convert2r prcp.dat|null '//
     +               'tmax.dat|null tmin.dat|null'
                write(0,*) 'converts 3 files to the RClimDex format'
                call abort
            end if
            if ( file.ne.'null' ) then
                call readseries(file,data(1,yrbeg,ivar),npermax,yrbeg
     +               ,yrend,nperyear,var,units,.false.,lwrite)
            end if
        end do
 
        call print3rdatfile(6,data,npermax,nperyear,yrbeg,yrend)
        end
