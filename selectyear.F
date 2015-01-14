        program selectyear
*
*       trivial little program to select a given range of years from a
*       standard station data file
*
        implicit none
#include "param.inc"
        integer yr,mo,i,nperyear,yr1,yr2
        real data(npermax,yrbeg:yrend)
        character file*255,line*128,var*40,units*20
        logical lvalid,lwrite
        integer iargc,llen
*
        lwrite = .false.
        if ( iargc().lt.3 .or. iargc().gt.4 ) then
            print *,'usage: selectyear yrbegin yrend datfile [dummy]'
            stop
        endif
        call getarg(1,line)
        read(line,*,err=901) yr1
        call getarg(2,line)
        read(line,*,err=901) yr2
        call getarg(3,file)
        call copyheader(file,6)
        call readseries(file,data,npermax,yrbeg,yrend,nperyear,var,units
     +       ,.false.,lwrite)
        call printdatfile(6,data(1,yr1),npermax,nperyear,yr1,yr2)
        goto 999
  901   write(0,*) 'selectyear: error reading year from ',line
        call abort
  999   continue
        end
