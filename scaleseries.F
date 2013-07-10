        program scaleseries
*
*       scale a time series
*
        implicit none
#include "param.inc"
        integer i,j,nperyear,ndpm,dpm(12,2)
        real data(npermax,yrbeg:yrend),factor,offset
        character string*1024,var*40,units*20
        logical lwrite
        integer iargc,llen,leap
        data dpm 
     +       /31,28,31,30,31,30,31,31,30,31,30,31
     +       ,31,29,31,30,31,30,31,30,30,31,30,31/

        lwrite = .false.
        if ( iargc().lt.2 ) then
            print *,'usage: scaleseries factor[:offset] file [ndpm]'
            stop
        endif
        call getarg(1,string)
        i = index(string,':')
        if ( i.eq.0 ) then
            read(string,*) factor
            offset = 0
        else
            read(string(:i-1),*) factor
            read(string(i+1:),*) offset
        end if
        if ( iargc().gt.2 ) then
            call getarg(3,string)
            read(string,*,err=10) ndpm
            goto 20
 10         continue
            ndpm = 0
 20         continue
        else
            ndpm = 0
        endif
        call getarg(2,string)
        call readseries(string,data,npermax,yrbeg,yrend,nperyear,var
     +       ,units,.false.,lwrite)
        if ( ndpm.ne.0 .and. nperyear.ne.12 ) then
            write(0,*) 'scaleseries: error: can only scale by dpm if '//
     +           'nperyear = 12, not ',nperyear
            write(*,*) 'scaleseries: error: can only scale by dpm if '//
     +           'nperyear = 12, not ',nperyear
            call abort
        endif
        call copyheader(string,6)
        write(6,'(a,g16.6)') '# scaled with a factor ',factor
        if ( ndpm.eq.1 ) then
            write(6,'(a)')'# and multiplied by '//
     +           'the number of days in a month'
        elseif ( ndpm.eq.2 ) then
            write(6,'(a)') '# and multiplied by '//
     +           'the number of days in a month squared'
        elseif ( ndpm.eq.-1 ) then
            write(6,'(a)')'# and divided by '//
     +           'the number of days in a month'
        elseif ( ndpm.eq.-2 ) then
            write(6,'(a)') '# and divided by '//
     +           'the number of days in a month squared'
        elseif ( ndpm.ne.0 ) then
            write(6,'(a,i3)') '# and multiplied by dpm**',ndpm
        endif
        if ( offset.ne.0 ) then
            write(6,'(a,g16.6)') '# added offset ',offset
        end if
        do i=yrbeg,yrend
            do j=1,nperyear
                if ( data(j,i).lt.1e33 ) then
                    data(j,i) = factor*data(j,i)
                    if ( ndpm.ne.0 ) then
                        data(j,i) = (real(dpm(j,leap(i)))**ndpm)*
     +                       data(j,i)
                    endif
                    data(j,i) = data(j,i) + offset
                endif
            enddo
        enddo
        call printdatfile(6,data,npermax,nperyear,yrbeg,yrend)
        end
