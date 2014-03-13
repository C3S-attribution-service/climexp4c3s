        program yearly2shorter
*
*       convert a yearly time series to a monthly one
*
        implicit none
#include "param.inc"
        integer yr,mo,dy,nperyear,npernew,m1,n,k,i,nfac
        real data(npermax,yrbeg:yrend),newdata(npermax,yrbeg:yrend)
        character file*256,var*40,units*20
        logical lwrite,lstandardunits
        integer iargc,llen
        
        lwrite = .false.
        lstandardunits = .false.
        if ( iargc().lt.2 ) then
            print *,'usage: yearly2shorter infile.dat nperyearnew '//
     +           '[mon n ave|sum m]'
            stop
        endif
        call getarg(1,file)
        call readseries(file,data,npermax,yrbeg,yrend,nperyear,var,units
     +       ,lstandardunits,lwrite)
        call copyheader(file,6)
        call getarg(2,file)
        read(file,*,err=901) npernew
        if ( iargc().gt.2 ) then
            call getarg(3,file)
            if ( file(1:3).ne.'mon' ) goto 902
            call getarg(4,file)
            read(file,*,err=903) m1
            call getarg(6,file)
            read(file,*,err=905) n
            call getarg(5,file)
            if ( file(1:3).eq.'ave' ) then
                nfac = 1
            elseif ( file(1:3).eq.'sum' ) then
                nfac = n
            else
                goto 904
            endif
        else
            m1 = 1
            n = npernew
            nfac = 1
        endif
        call annual2shorter(data,npermax,yrbeg,yrend,nperyear,
     +       newdata,npermax,yrbeg,yrend,npernew,m1,n,nfac,lwrite)
*
        call printdatfile(6,newdata,npermax,npernew,yrbeg,yrend)
*
        goto 999
 901    write(0,*) 'yearly2shorter: error reading npernew from '
     +       ,file(1:llen(file))
        call abort
 902    write(0,*) 'yearly2shorter: error: expecting ''month'', not '
     +       ,file(1:llen(file))
        call abort
 903    write(0,*) 'yearly2shorter: error reading first month from '
     +       ,file(1:llen(file))
        call abort
 904    write(0,*) 'yearly2shorter: error: expecting ''ave|sum'', not '
     +       ,file(1:llen(file))
        call abort
 905    write(0,*)
     +       'yearly2shorter: error reading number of months from '
     +       ,file(1:llen(file))
        call abort
 999    continue
        end
