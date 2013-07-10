        program extractseries
*  #[ comment:
*
*       Convert a time series in my standard format to a table for R
*
*       (c) Geert Jan van Oldenborgh, 7-jun-2006, KNMI
*
*       This file may be copied and modified  freely as long as the
*       above copyright notice is kept intact
*
*  #] comment:
*  #[ declarations:
*
        implicit none
#include "param.inc"
#include "getopts.inc"
        integer nperyear,mens1,mens,mo,yr,i,j,j1,j2,mo1,month,iu,iens
     +       ,nperyear2,yrmin,yrmax,n,nmodel1,nmodel2,ndup(npermax)
     +       ,validens(npermax)
        real fcst(npermax,yrbeg:yrend,0:nensmax)
        real s1,s2
        character line*255,varfcst*20,unitfcst*40
        integer iargc
*
*  #] declarations:
*  #[ check arguments:
*
*       check arguments
*       
        lwrite = .FALSE.
        lstandardunits = .TRUE.
        if ( iargc().lt.4 ) then
            print *,'usage: extractseries infile '//
     +           '[month m[:n] [sum|ave|max|min|sel m] '//
     +           '[begin yr] [end yr] [detrend]'//
     +           '[debias none|mean|var|all] dump outtable'
            stop
        endif
*  #] check arguments:
*  #[ read forecast:
*
        call getarg(1,line)
        call readensseries(line,fcst,npermax,yrbeg,yrend,nensmax
     +       ,nperyear,mens1,mens,varfcst,unitfcst,lstandardunits
     +       ,lwrite)
        print '(a,i2,a,i2)','# located ensemble members ',mens1,' to '
     +       ,mens
*
*  #] read forecast:
*  #[ options:
*
        call getopts(2,iargc(),nperyear,yrbeg,yrend,.true.,mens1,mens)
        if ( .not.dump ) then
            write(0,*) 'extractseries: error: cannot find ''dump'''//
     +           ' in argument list'
            write(*,*) 'extractseries: error: cannot find ''dump'''//
     +           ' in argument list'
            call abort
        endif
        if ( lag1.ne.0 .or. lag2.ne.0 ) then
            write(0,*) 'extractseries: error: cannot handle lags ',lag1
     +           ,lag2
            write(*,*) 'extractseries: error: cannot handle lags ',lag1
     +           ,lag2
            call abort
        endif
*
*  #] options:
*  #[ manipulate time series:
*
*       anomalies
*       
        if ( anom ) then
            if ( lwrite ) print '(a)','# Taking anomalies '
            do iens=nens1,nens2
                call anomal(fcst(1,yrbeg,iens),nperyear,nperyear
     +               ,yrbeg,yrend,yr1,yr2)
            enddo
        endif
*
*       sum
*
        if ( lsum.gt.1 ) then
            if ( lwrite ) print '(a,i3)','# Summing series ',lsum
            do iens=nens1,nens2
                call sumit(fcst(1,yr1,iens),nperyear,nperyear
     +               ,yr1,yr2,lsum2,'v')
            enddo
        endif
*
*       detrending
*
        if ( ldetrend ) then
            do iens=nens1,nens2
                if ( lwrite ) print *,'Detrending ens ',iens
                call detrend(fcst(1,yrbeg,iens),nperyear
     +               ,nperyear,yrbeg,yrend,yr1,yr2,m1,m2,lsel)
            enddo
        endif
        if ( ndiff.ne.0 ) then
            if ( lwrite ) print *,'Taking differences - series'
            do iens=nens1,nens2
                call diffit(fcst(1,yrbeg,iens),npermax,nperyear,yrbeg
     +                ,yrend,ndiff)
            enddo
        endif
*
*       copy ensemble members so that there is the same
*       number of valid ones at every time step
*       
        ndup = 0
        if ( nens2.gt.nens1 .and. lmakeensfull ) then
            call makeensfull(ndup,nperyear,fcst,nperyear,yrbeg
     +           ,yrend,nens1,nens2,validens,lwrite)
        endif
*
*  #] manipulate time series:
*  #[ make extractseries table:
*
        do i=1,len_trim(varfcst)
            if ( varfcst(i:i).eq.' ' ) varfcst(i:i) = '_'
        enddo
        write(10,'(a)') varfcst
        do yr=yr1,yr2
            do month=m1,m2
                call getj1j2(j1,j2,month,nperyear,.false.)
                do mo1=j1,j2
                    mo = mo1
                    if ( mo.gt.12 ) mo = mo - 12
                    do iens=nens1,nens2
                        if ( fcst(mo,yr,iens).eq.3e33 )
     +                       fcst(mo,yr,iens) = -999.9
                    enddo       ! iens
***                    write(*,'(2i4,100g14.6)') yr,mo,(fcst(mo,yr,iens),
***     +                   iens=nens1,nens2)
                    write(10,'(100g14.6)') (fcst(mo,yr,iens),
     +                   iens=nens1,nens2)
                enddo           ! mo!
            enddo               ! month
        enddo                   ! yr
        close(10)
*
*  #] make extractseries table:
        end
