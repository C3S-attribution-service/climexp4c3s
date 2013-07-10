        program operate
*
*       Apply a logarithm, sqrt, ... to a timeseries in standard
*       format and return the resulting series.
*       
        implicit none
#include "param.inc"
        integer i,j,yr,mo,nyr,nperyear
        real data(npermax,yrbeg:yrend)
        character file*256,line*128,oper*3,var*40,units*20
        logical lwrite,lstandardunits
        integer iargc,llen

        lwrite = .false.
        lstandardunits = .false.
        if ( iargc().lt.2 ) then
            write(0,*)
     +            'usage: operate log|sqrt|exp file'
            stop
        endif
        call getarg(1,oper)
        call getarg(2,file)
        write(*,'(4a)') '# ',oper,' of ',trim(file)
        call copyheader(file,6)
        call readseries(file,data,npermax,yrbeg,yrend,nperyear,var,units
     +       ,lstandardunits,lwrite)
*
        if ( oper.eq.'log' ) then
            call takelog(data,npermax,nperyear,yrbeg,yrend)
        elseif ( oper.eq.'sqrt' ) then
            call takesqrt(data,npermax,nperyear,yrbeg,yrend)
        elseif ( oper.eq.'exp' ) then
            call takeexp(data,npermax,nperyear,yrbeg,yrend)
        elseif ( oper.eq.'inv' ) then
            call takeinv(data,npermax,nperyear,yrbeg,yrend)
        else
            write(0,'(2a)') 'operate: unknown operation ',oper
            call abort
        endif
*       
        call printdatfile(6,data,npermax,nperyear,yrbeg,yrend)
*
        end
