        program diffdat
*
        implicit none
#include "param.inc"
        integer nperyear,n
        real data(npermax,yrbeg:yrend),diff(npermax,yrbeg:yrend),minfac
        character file*256,var*80,units*40,string*10
        logical lwrite

        minfac = 0.5            ! arbitrary
        call getarg(1,file)
        if ( file.eq.' ' ) then
            write(0,*) 'usage: diffdat file [n [debug]]'
            call abort
        end if
        call getarg(2,string)
        lwrite = .false.
        if ( string.ne.' ' ) then
            read(string,*) n
            if ( n.lt.3 ) n = 3
            if ( mod(n,2).ne.1 ) n = n + 1
            call getarg(3,string)
            if ( string.eq.'lwrite' .or. string.eq.'debug' ) then
                lwrite = .true.
            end if
        else
            n = 3
        end if
        call readseries(file,data,npermax,yrbeg,yrend,nperyear,var,units
     +       ,.false.,lwrite)
        call derivative(n,data,diff,npermax,yrbeg,yrend,nperyear,minfac
     +       ,lwrite)
        print '(3a,i3,a)','# centered derivative of ',trim(file)
     +       ,' using ',n,' data points'
        call copyheader(file,6)
        call printdatfile(6,diff,npermax,nperyear,yrbeg,yrend)
        end program

            
