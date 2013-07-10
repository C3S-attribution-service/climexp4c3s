        program filterseries
*
*       Apply hi- or lo-pass filters to a time series in standard
*       format and return the resulting series.
*       
        implicit none
#include "param.inc"
#include "getopts.inc"
        integer i,j,yr,mo,nmonth,nperyear,mens1,mens
        real data(npermax,yrbeg:yrend)
        character file*256,line*256,hilo*2,filtertype*12,var*40,units*20
        integer iargc

        lwrite = .false.
        if ( iargc().lt.4 ) then
            write(0,*)
     +            'usage: filtermonthseries hi|lo filtertype nmon file'
            stop
        endif
        call getarg(1,hilo)
        if ( hilo.ne.'hi' .and. hilo.ne.'lo' ) then
            write(0,*) 'filterseries: error: say hi or lo, not ',hilo
            call abort
        endif
        call getarg(2,filtertype)
        call getarg(3,line)
        read(line,*,err=901) nmonth
        call getarg(4,file)
        call readseries(file,data,npermax,yrbeg,yrend,nperyear,var,units
     +       ,.false.,lwrite)
        call getopts(5,iargc()-1,nperyear,yrbeg,yrend,.true.,mens1,mens
     +       )
        if ( minfac.lt.0 ) minfac = 0.75
        call copyheader(file,6)
        if ( nperyear.eq.12 ) then
            write(*,'(5a,i4,a)') '# filtered with a ',hilo,
     +           '-pass ',filtertype,' filter with cut-off ',nmonth
     +           ,' months'
        elseif ( nperyear.eq.4 ) then
            write(*,'(5a,i4,2a)') '# filtered with a ',hilo,
     +           '-pass ',trim(filtertype),' filter with cut-off '
     +           ,nmonth,' seasons '
        elseif ( nperyear.eq.360 .or. nperyear.eq.365 .or.
     +           nperyear.eq.366 ) then
            write(*,'(5a,i4,a)') '# filtered with a ',hilo,
     +           '-pass ',trim(filtertype),' filter with cut-off '
     +           ,nmonth,' days'
        else
            write(*,'(5a,i4,a)') '# filtered with a ',hilo,
     +           '-pass ',trim(filtertype),' filter with cut-off '
     +           ,nmonth,' periods'
        endif
*
        if ( filtertype.eq.'running-mean' .or. filtertype.eq.'box' )
     +		then
            if ( hilo.eq.'hi' ) then
                call mhipass(data,npermax,nperyear,yrbeg,yrend,nmonth-1
     +               ,minfac)
            else
                call sumit(data,npermax,nperyear,yrbeg,yrend,nmonth,'v')
                call shiftseries(data,npermax,nperyear,yrbeg,yrend,
     +                nmonth/2)
            endif
        else if ( filtertype.eq.'loess1' .or. filtertype.eq.'loess2' )
     +           then
            call myloess(data,npermax,nperyear,yrbeg,yrend,nmonth/2
     +           ,minfac,filtertype,hilo,'month','gaussian',lwrite)
        else
            write(0,*) 'filterseries: error: filtertype ',filtertype
     +            ,' not yet implemented'
            call abort
        endif
*       
        call printdatfile(6,data,npermax,nperyear,yrbeg,yrend)
*
        goto 999
  901   write(0,*) 'filterseries: expcting an integer, not ',file
        call abort
  999   continue
        end
