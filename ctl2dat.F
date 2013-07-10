        program ctl2dat
*
*       convert a time series in GrADS format to my ASCII format
*       (netcdf2dat cannot handle days yet)
*
        implicit none
#include "recfac.h"
#include "param.inc"
        integer i,yr,mo,nx,ny,nz,nt,nperyear,yrbegin,mobegin,endian
     +       ,nvars,ivars(2,1),month,day
        real xx(1),yy(1),zz(1),undef
        real data(npermax,yrbeg:yrend)
        character file*256,datfile*256,vars(1)*10,title*128,lvars(1)*80
     +       ,units(1)*10
        integer leap,iargc,llen
        
        if ( iargc().ne.1 ) then
            print *,'usage: ctl2dat in.ctl'
            stop
        endif
        call getarg(1,file)
        call parsectl(file,datfile,1,nx,xx,1,nx,yy,1,nz,zz,nt,nperyear
     +       ,yrbegin,mobegin,undef,endian,title,1,nvars,vars,ivars
     +       ,lvars,units)
        print '(2a)','# ',title(1:llen(title))
        print '(2a)','# ',vars(1)(1:llen(vars(1)))
        print '(2a)','# ',lvars(1)(1:llen(lvars(1)))
        print '(a,f8.2,a,f8.2,a,f8.2)','# at lon ',xx(1),', lat ',yy(1)
     +       ,', lev ',zz(1)
        print '(a)','#'
        call makeabsent(data,npermax,yrbeg,yrend)        
        open(1,file=datfile,access='direct',recl=recfa4,status='old')
        yr = yrbegin
        mo = mobegin
        do i=1,nt
            if ( nperyear.eq.366 .and. mo.eq.31+29 .and. leap(yr).eq.1
     +           ) then
                mo = mo + 1
            endif

            read(1,rec=i) data(mo,yr)
            call monthday(yr,mo,nperyear,month,day)
            print '(i4,2i3,g14.6)',yr,month,day,data(mo,yr)

            mo = mo + 1
            if ( mo.gt.nperyear ) then
                mo = mo - nperyear
                yr = yr + 1
            endif
        enddo
**        call printdatfile(6,data,npermax,nperyear,yrbeg,yrend)
        end

        subroutine monthday(yr,mo,nperyear,month,day)
        implicit none
        integer yr,mo,nperyear,month,day
        integer dpm(12)
        data dpm /31,29,31,30,31,30,31,31,30,31,30,31/
        if ( nperyear.eq.12 ) then
            month = mo
            day = 15
        elseif ( nperyear.eq.366 ) then
            day = mo
            do month=1,12
                if ( day - dpm(month).le.0 ) then
                    return
                endif
                day = day - dpm(month)
            enddo
        else
            write(0,*) 'momthday: error : cannot handle nperyear = '
     +           ,nperyear,'yet'
            call abort
        endif
        end
