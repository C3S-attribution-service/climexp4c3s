        program filteryearseries
*
*       Apply hi- or lo-pass filters to a time series in standard
*       format and return the resulting series.
*       
        implicit none
#include "param.inc"
#include "getopts.inc"
        integer i,j,yr,mo,nyr,nperyear,n,mens,mens1
        real data(npermax,yrbeg:yrend)
        character file*256,line*128,hilo*2,filtertype*12,var*20,units*20
        integer iargc

        lwrite = .false.
        if ( iargc().lt.4 ) then
            write(0,*) 'usage: filteryearseries hi|lo filtertype nyr '//
     +           ' [minfac N] file garbage'
            stop
        endif
        call getarg(1,hilo)
        if ( hilo.ne.'hi' .and. hilo.ne.'lo' ) then
            write(0,*) 'filterseries: error: say hi or lo, not ',hilo
            call abort
        endif
        call getarg(2,filtertype)
        call getarg(3,line)
        read(line,*,err=901) nyr
        call getarg(4,file)
        call readseries(file,data,npermax,yrbeg,yrend,nperyear,var,
     +       units,.false.,lwrite)
        call copyheader(file,6)
        write(6,'(3a,i4,3a)') '# each calendar month ',hilo,
     +       '-pass filtered with a ',nyr,'-yr '
     +       ,trim(filtertype),' filter'
        call getopts(5,iargc()-1,nperyear,yrbeg,yrend,.true.,mens1,mens
     +       )
        if ( minfac.le.0 ) minfac = 0.75
*
        if ( filtertype.eq.'running-mean' .or. 
     +       filtertype(1:3).eq.'box' )
     +        then
            if ( hilo.eq.'hi' ) then
                call hipass(data,npermax,nperyear,yrbeg,yrend,nyr,minfac
     +               )
            else
                call ndiffit(data,npermax,nperyear,yrbeg,yrend,-nyr+1,
     +               minfacsum)
                call shiftseries(data,npermax,nperyear,yrbeg,yrend,
     +                -nperyear*((nyr-1)/2))
            endif
        else if ( filtertype.eq.'loess1' .or. filtertype.eq.'loess2' )
     +           then
            call myloess(data,npermax,nperyear,yrbeg,yrend,nyr/2
     +           ,minfac,filtertype,hilo,'year','gaussian',lwrite)
        else
            write(0,*) 'filterseries: error: filtertype ',filtertype
     +            ,' not yet implemeneted'
            call abort
        endif
*
        call printdatfile(6,data,npermax,nperyear,yrbeg,yrend)
*
        goto 999
 901    write(0,*) 'filteryearseries: expcting an integer, not '
     +       ,trim(line)
        call abort
 902    write(0,*) 'filteryearseries: expcting a number, not '
     +       ,trim(line)
        call abort
  999   continue
        end
