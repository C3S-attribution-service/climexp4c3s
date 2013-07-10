        program daily2longer
*
*       compute aggregate quantities from daily data
*       input: daily time series
*       output: yearly/monthly/10-dy/5-dy time series
*
        implicit none
#include "param.inc"
        integer mpermax
        parameter(mpermax=24*366)
        integer nperyear,nperyearnew,yr,mo,dy,i,j,n,itype
        integer*2 nn(mpermax)
        real olddata(mpermax,yrbeg:yrend),newdata(npermax,yrbeg:yrend),
     +          cut(mpermax),pcut,normdata(npermax),s,oldclim(mpermax),
     +          newclim(npermax)
        character file*512,string*512,lgt*1,oper*3,var*10,units*20
     +       ,climunits*10
        logical normal,lwrite,lstandardunits,lnomissing
        integer iargc,llen
        lwrite = .true.
        lstandardunits = .true.
        lnomissing=.false.
*
        if ( iargc().lt.3 ) then
            print *,'usage: daily2longer infile nperyearnew '//
     +          'mean|sd|sum|abo|bel|min|max|num [<> val[%|p]]'
            print *,'(more options will come as requested)'
            stop
        endif
*
*       read data
*
        call getarg(1,file)
        call readseries(file,olddata,mpermax,yrbeg,yrend,nperyear,var
     +       ,units,lstandardunits,lwrite)
        if ( index(file,'/dd').gt.0 .or. index(file,'dd').eq.1 ) then
            write(0,*) 'Hmm. This looks like a (wind) direction to me.'
            write(0,*) 'Averaging over a unit circle.<p>'
            itype = 1
        else
            itype = 0
        endif
*
*       read operation
*
        call getarg(2,string)
        read(string,*,err=901) nperyearnew
        if ( nperyearnew.gt.npermax ) then
            write(0,*) 'daily2longer: error: nperyearnew = ',nperyearnew
     +           ,' not yet supported'
            write(*,*) 'daily2longer: error: nperyearnew = ',nperyearnew
     +           ,' not yet supported'
            call abort
        endif
        call getarg(3,oper)
        if ( oper.ne.'mea' .and. oper.ne.'sd ' .and.
     +       oper.ne.'min' .and. oper.ne.'max' .and.
     +       oper.ne.'num' .and. oper.ne.'sum' .and.
     +       oper.ne.'bel' .and. oper.ne.'abo' ) then
            write(0,*) 'daily2longer: error: unknown operation ',oper
            call abort
        endif
        normal = .false.
        if ( iargc().eq.4 ) then
            call getarg(4,string)
            if ( string(1:5).eq.'nomis' ) lnomissing = .true.
        endif
        if ( iargc().ge.5 ) then
            call getarg(4,lgt)
            if ( lgt.eq.'l' ) lgt = '<'
            if ( lgt.eq.'g' ) lgt = '>'
            if ( lgt.ne.'<' .and. lgt.ne.'>' ) then
                write(0,*) 'daily2longer: error: unknown comparison '
     +               ,lgt
                call abort
            endif
            call getarg(5,string)
            if ( string.eq.'n' ) then
*               take normals wrt 1971-2000
                do j=1,nperyear
                    normdata(j) = 0
                    n = 0
                    do yr=1971,2000
                        if ( olddata(j,yr).lt.1e33 ) then
                            n = n + 1
                            normdata(j) = normdata(j) + olddata(j,yr)
                        endif
                    enddo
                    if ( n.gt.5 ) then ! arbitrary number
                        normdata(j) = normdata(j)/n
                    else
                        normdata(j) = 3e33
                    endif
                    if ( lwrite ) print *,j,normdata(j)
                enddo
*               no smoothing for the time being
                do yr=yrbeg,yrend
                    do j=1,nperyear
                        if ( olddata(j,yr).lt.1e33 .and. 
     +                       normdata(j).lt.1e33 ) then
                            olddata(j,yr) = olddata(j,yr) - normdata(j)
                        else
                            olddata(j,yr) = 3e33
                        endif
                    enddo
                enddo
                if ( lwrite ) then
                    do j=1,nperyear
                        s = 0
                        n = 0
                        do yr=1971,2000
                            if ( olddata(j,yr).lt.1e33 ) then
                                n = n + 1
                                s = s + olddata(j,yr)
                            endif
                        enddo
                        if ( n.gt.0 ) then
                            print *,j,s/n,n
                        endif
                    enddo
                endif
                lgt = ' '
            else
                read(string,*,err=902) pcut
                if ( index(string,'%').ne.0 .or. index(string,'p').ne.0
     +               ) then
                    do j=1,nperyear
                        call getcutoff(cut(j),pcut,olddata,mpermax
     +                       ,nperyear,yrbeg,yrend,yrbeg,yrend,j,j,0)
                    enddo
                else
                    do j=1,nperyear
                        cut(j) = pcut
                    enddo
                endif
            endif
        else
            lgt = ' '
        endif
*
*       compute climatology and anomalies
*
        if ( oper.eq.'mea' .and. lgt.eq.' ' .and. .not.lnomissing ) then
            do j=1,nperyear
                oldclim(j) = 0
            enddo
            do j=1,nperyear
                nn(j) = 0
            enddo
            do yr=yrbeg,yrend
                do j=1,nperyear
                    if ( olddata(j,yr).lt.1e33 ) then
                        nn(j) = nn(j) + 1
                        oldclim(j) = oldclim(j) + olddata(j,yr)
                    endif
                enddo
            enddo
            do j=1,nperyear
                if ( nn(j).gt.0 ) then
                    oldclim(j) = oldclim(j)/nn(j)
                else
                    oldclim(j) = 3e33
                endif
            enddo
            do yr=yrbeg,yrend
                do j=1,nperyear
                    if ( olddata(j,yr).lt.1e33 ) then
                        olddata(j,yr) = olddata(j,yr) - oldclim(j) 
                    endif
                enddo
            enddo
            climunits = units
            call allday2period(
     +          oldclim,mpermax,nperyear,
     +          newclim,npermax,nperyearnew,
     +          0,0,oper,lgt,cut,itype,var,climunits,lwrite)
        else
            do j=1,nperyearnew
                newclim(j) = 0
            enddo
        endif
*
*       perform operation
*
        call makeabsent(newdata,npermax,yrbeg,yrend)
        call allday2period(
     +     olddata,mpermax,nperyear,
     +     newdata,npermax,nperyearnew,
     +     yrbeg,yrend,oper,lgt,cut,itype,var,units,lwrite)
        do yr=yrbeg,yrend
            do j=1,nperyearnew
                if ( newdata(j,yr).lt.1e33 ) then
                    newdata(j,yr) = newdata(j,yr) + newclim(j)
                endif
            enddo
        enddo
*
*       print output
*
        if ( lgt.eq.' ' ) then
            print '(a,i4,2a)','# daily2longer ',nperyearnew,' ',oper
        else
            print '(a,i4,3a,g20.4)','# daily2longer ',nperyearnew,' ',
     +          oper,lgt,pcut
        endif
        open(1,file=file)
        do
            read(1,'(a)') file
            if ( file(1:1).ne.'#' ) exit
            call adjustunitstring(file,units)
            print '(a)',trim(file)
        enddo
        close(1)
        call printdatfile(6,newdata,npermax,nperyearnew,yrbeg,yrend)
*
*       error messages
*
        goto 999
  901   write(0,*) 'daily2longer: expecting nperyearnew, not ',string
        call abort
  902   write(0,*) 'daily2longer: expecting value[%|p], not ',string
        call abort
  999   continue
        end

        subroutine adjustunitstring(string,units)
        implicit none
        character string*(*),units*(*)
        integer i,j
        i = index(string,'[')
        if ( i.ne.0 ) then
            j = index(string,']')
            if ( j.ne.0 ) then
                string = string(:i)//trim(units)//string(j:)
            endif
        endif
        end
