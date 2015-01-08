        program diffdat
*
        implicit none
#include "param.inc"
        integer nperyear,n
        real data(npermax,yrbeg:yrend),diff(npermax,yrbeg:yrend),minfac
        character file*256,var*80,units*40,string*10
        logical lwrite

        minfac = 0.7            ! arbitrary
        call getarg(1,file)
        if ( file.eq.' ' ) then
            write(0,*) 'usage: diffdat file [n [debug]]'
            call abort
        end if
        call getarg(2,string)
        lwrite = .false.
        if ( string.ne.' ' ) then
            read(string,*) n
            if ( n.lt.2 ) n = 2
            if ( n.ge.3 .and. mod(n,2).ne.1 ) n = n + 1
            call getarg(3,string)
            if ( string.eq.'lwrite' .or. string.eq.'debug' ) then
                lwrite = .true.
            end if
        else
            n = 3
        end if
        call readseries(file,data,npermax,yrbeg,yrend,nperyear,var,units
     +       ,.false.,lwrite)
        if ( n.eq.2 ) then
            call ndiffit(data,npermax,nperyear,yrbeg,yrend,1,minfac)
            diff = data
        else
            call derivative(n,data,diff,npermax,yrbeg,yrend,nperyear,
     +           minfac,lwrite)
        end if
        if ( n.eq.2 ) then
            print '(3a)','# difference with previous year of ',
     +           trim(file)
        else
            print '(3a,i3,a)','# centered derivative of ',trim(file)
     +       ,' using linear regression over ',n,' data points'
        end if
        call copyheader(file,6)
        call printdatfile(6,diff,npermax,nperyear,yrbeg,yrend)
        end program

            
