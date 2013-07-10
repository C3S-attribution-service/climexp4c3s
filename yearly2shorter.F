        program yearly2shorter
*
*       convert a yearly time series to a monthly one
*
        implicit none
#include "param.inc"
        integer yr,mo,dy,nperyear,npernew,m1,n,k,i,nfac,dpm(12,2)
        real data(npermax,yrbeg:yrend),newdata(npermax,yrbeg:yrend)
        character file*256,var*40,units*20
        logical lwrite,lstandardunits
        integer iargc,llen,leap
        data dpm /31,28,31,30,31,30,31,31,30,31,30,31,
     +            31,29,31,30,31,30,31,31,30,31,30,31/
        
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
*
        call makeabsent(newdata,npermax,yrbeg,yrend)
        if ( nperyear.eq.1 ) then
            do yr=yrbeg,yrend
                do k=1,n
                    mo = m1 + k - 1
                    call normon(mo,yr,i,npernew)
                    if ( i.ge.yrbeg .and. i.lt.yrend ) then
                        if ( data(1,yr).lt.1e33 ) then
                            newdata(mo,i) = data(1,yr)/nfac
                        endif
                    endif
                enddo
            enddo
        elseif ( mod(npernew,nperyear).eq.0 ) then
            if ( nfac.ne.1 ) nfac = npernew/nperyear
            do yr=yrbeg,yrend
                do mo=1,npernew
                    k = 1 + (mo-1)*nperyear/npernew
                    if ( data(k,yr).lt.1e33 ) then
                        newdata(mo,yr) = data(k,yr)/nfac
                    endif
                enddo
            enddo
        elseif ( npernew.eq.366 .and. nperyear.eq.12 ) then
            do yr=yrbeg,yrend
                k = 0
                do mo=1,12
                    if ( nfac.ne.1 ) nfac = dpm(mo,leap(yr))
                    do dy=1,dpm(mo,2)
                        k = k + 1
                        if ( k.eq.60 .and. leap(yr).eq.1 ) then
                            newdata(k,yr) = 3e33
                        elseif ( data(mo,yr).lt.1e33 ) then
                            newdata(k,yr) = data(mo,yr)/nfac
                        endif
                    enddo
                enddo
            enddo
        else
            write(0,*) 'yearly2longer: error: cannet handle ',
     +           'conversion from ',nperyear,' to ',npernew ,
     +           ' points per year yet'
            write(*,*) 'yearly2longer: error: cannet handle ',
     +           'conversion from ',nperyear,' to ',npernew ,
     +           ' points per year yet'
        endif
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
