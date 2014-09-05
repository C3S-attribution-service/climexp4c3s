        program getval
*
*       retrieve a value from a time series, 
*       accompanied by mean and standard deviation
*
        implicit none
#include "param.inc"
#include "getopts.inc"
        integer unit,yr,mo,n,nperyear,i,j,m
        real data(npermax,yrbeg:yrend),xx(yrend-yrbeg+1),x1,s1
        logical lopen
        character file*256,var*20,units*20
        integer iargc,llen
        
        lwrite = .false.
        if ( iargc().lt.5 ) then
            print *,'usage: getval file beg yr mon mo [name name]'
            stop
        endif
        call getarg(1,file)
        call readseries(file,data,npermax,yrbeg,yrend,nperyear,var,units
     +       ,lstandardunits,lwrite)
        if ( nperyear.eq.0 ) goto 999
        call getopts(2,iargc(),nperyear,yrbeg,yrend,.false.,0,0)
        if ( lsum.gt.1 ) then
            call sumit(data,npermax,nperyear,yrbeg,yrend,lsum,oper)
        endif
        n = 0
        do yr=yrbeg,yrend
            if ( data(m1,yr).lt.1e33 ) then
                n = n + 1
                xx(n) = data(m1,yr)
            endif
        enddo
        call getmoment(1,xx,n,x1)
        call getmoment(2,xx,n,s1)
        inquire(unit=11,opened=lopen)
        if ( .not.lopen ) then
            unit = 6
        else
            unit = 11
        endif
        i = index(file,'/',.true.) + 1
        j = index(file,'.') - 1
        if ( yr2a.gt.yrend ) then
            write(*,*) 'getval: error: year undefined'
            write(0,*) 'getval: error: year undefined'
            call abort
        endif
        call invgetdymo(day0,m1,m,nperyear)
        if ( lwrite ) then
            print *,'getval: dy,mo  ',day0,m1
            print *,'getval: data(',m,yr2a,') ',data(m,yr2a)
        end if
        if ( data(m,yr2a).lt.1e30 ) then
            write(*,'(3a,g16.4,3a)')
     +           '<table class=realtable width=451><tr><td width=300>'
     +           ,trim(namestring),'<td>',data(m,yr2a),' '
     +           ,trim(units),'</table>'
        endif
        write(unit,'(2i5,g12.4,f6.1,i2,2g12.4,3f6.1)') 
     +       yr2a,m1,data(m,yr2a),-1.,0,x1,s1,0.,0.,0.
 999    continue
        end
