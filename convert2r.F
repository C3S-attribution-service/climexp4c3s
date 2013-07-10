        program convert2r
!
!       convert an (almost) arbitrary time series to the format the R
!       routines from CCCMA require
!
        implicit none
#include "param.inc"
        integer yr1,yr2,yr,mo,dy,nperyear
        real data(npermax,yrbeg:yrend)
        character file*255,var*40,units*40
        logical lwrite
        lwrite = .false.

        call getarg(1,file)
        if ( file.eq.' ' ) then
            write(0,*) 'convert2r file'
            write(0,*) 'converts the file to the format required by R'
            call abort
        end if
        call readseries(file,data,npermax,yrbeg,yrend,nperyear,
     +       var,units,.false.,lwrite)
        yr1 = yrend
        yr2 = yrbeg
        do yr=yrbeg,yrend
            do mo=1,nperyear
                if ( data(mo,yr).lt.1e33 ) then
                    yr1 = min(yr1,yr)
                    yr2 = max(yr2,yr)
                end if
            end do
        end do
        call printrdatfile(6,data,npermax,nperyear,yrbeg,yrend
     +       ,yr1,yr2)
        end
