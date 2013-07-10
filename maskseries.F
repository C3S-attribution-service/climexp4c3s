        program normdiff
*
*       mask out part of a series based on values of another series
*       e.g. all months in the NAO when NINO3>0
*
        implicit none
#include "param.inc"
#include "getopts.inc"
        integer i,j,nperyear,my1,my2,n(npermax)
        real  data(npermax,yrbeg:yrend),mask(npermax,yrbeg:yrend)
        character file*255,var*40,units*20,line*255
        integer iargc,llen,rindex

        if ( iargc().lt.2 ) then
            print *,'usage: maskseries series1 series2 options'
            stop
        endif
        
        lwrite = .false.
        call getarg(1,file)
        call copyheader(file,6)
        call readseries(file,data,npermax,yrbeg,yrend,nperyear,var,units
     +       ,lstandardunits,lwrite)
        call getarg(2,file)
        call readseries(file,mask,npermax,yrbeg,yrend,i,var,units
     +       ,lstandardunits,lwrite)
        if ( i.ne.nperyear ) then
            write(0,*)
     +            'maskseries: error: cannot inperpolate in time (yet)'
            write(*,*)
     +            'maskseries: error: cannot inperpolate in time (yet)'
            call abort
        endif
        call getopts(3,iargc(),nperyear,yrbeg,yrend,.false.,0,0)
        i = rindex(file,'/')
        j = index(file,'.dat')
        if ( j.eq.0 ) then
            j = index(file,'.nc')
            if ( j.eq.0 ) then
                j = len_trim(file) + 1
            endif
        endif
        if ( maxindx.lt.1e33 ) then
            if ( minindx.gt.-1e33 ) then
                print '(a,g10.2,3a,g10.3)','# only when ',minindx
     +               ,' &llt; ',trim(file(i+1:j-1)),' &lt; ',maxindx
            else
                print '(3a,g10.3)','# for ',trim(file(i+1:j-1)),' &lt; '
     +               ,maxindx
            endif
        else
            if ( minindx.gt.-1e33 ) then
                print '(3a,g10.3)','# for ',trim(file(i+1:j-1)),' &gt; '
     +               ,minindx
            else
                write(*,*) 'maskseries: error: expected ''lt val'' '//
     +                ' or ''gt val'''
                write(0,*) 'maskseries: error: expected ''lt val'' '//
     +                ' or ''gt val'''
                call abort
            endif
        endif
        if ( mdiff.gt.0 ) then
            if ( lwrite ) print '(a)','# Taking monthly anomalies'
            call mdiffit(mask,npermax,nperyear,yrbeg,yrend,mdiff)
        endif
        if ( lsum.gt.1 ) then
            if ( lwrite ) print '(a,i3)','# Summing ',lsum
            call sumit(mask,npermax,nperyear,yrbeg,yrend,lsum,oper)
        endif
        if ( ldetrend ) then
            if ( lwrite ) print *,'Detrending'
            call detrend(mask,npermax,nperyear,yrbeg,yrend,yr1,yr2,
     +           m1,m2,lsel)
        endif
        if ( ndiff.ne.0 ) then
            if ( lwrite ) print *,'Taking differences - data'
            call diffit(mask,npermax,nperyear,yrbeg,yrend,ndiff)
        endif
        if ( anom .or. 
     +        (lsel.gt.1 .or. nfittime.gt.0 ) .and. ndiff.le.0 ) then
            if ( lwrite ) print *,'Taking anomalies'
            call anomal(mask,npermax,nperyear,yrbeg,yrend,yr1,yr2)
        endif
       
        do i=yr1,yr2
            do j=1,nperyear
                if ( mask(j,i).lt.maxindx .neqv. mask(j,i).gt.minindx
     +                ) then
                    data(j,i) = 3e33
                endif
            enddo
        enddo
        call printdatfile(6,data,npermax,nperyear,yrbeg,yrend)
        end
