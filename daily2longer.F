        program daily2longer
*
*       compute aggregate quantities from daily data
*       input: daily time series
*       output: yearly/monthly/10-dy/5-dy time series
*
        implicit none
#include "param.inc"
#include "getopts.inc"
        integer mpermax
        parameter(mpermax=24*366)
        integer nperyear,nperyearnew,yr,mo,dy,i,j,n,itype,nperyear2
        integer*2 nn(mpermax)
        real olddata(mpermax,yrbeg:yrend),newdata(npermax,yrbeg:yrend),
     +       cut(mpermax),pcut,normdata(npermax),s,oldclim(mpermax),
     +       oldclim2(mpermax),newclim(npermax),newclim2(mpermax)
        character file*1023,string*512,lgt*1,moper*3,var*20,units*20,
     +       climunits*10
        integer iargc,leap
        external leap
        lwrite = .false.
        lstandardunits = .true.
        lnomissing=.false.
*
        if ( iargc().lt.3 ) then
            print *,'usage: daily2longer infile nperyearnew'//
     +           ' mean|sd|sum|abo|bel|num|min|max|mintime|maxtime'//
     +           '|firsttime|lasttime|con [<> val[%|p]'//
     +           ' add_anom|add_clim|add_trend|add_damped]'
            print *,'(more options will come as requested)'
            stop
        endif
*
*       read data
*
        call getarg(3,string)
        if ( string.eq.'mintime' ) string = 'nti'
        if ( string.eq.'maxtime' ) string = 'xti'
        if ( string.eq.'firsttime' ) string = 'fti'
        if ( string.eq.'lasttime' ) string = 'lti'
        moper = string
        if ( moper.ne.'mea' .and. moper.ne.'sd ' .and.
     +       moper.ne.'min' .and. moper.ne.'max' .and.
     +       moper.ne.'nti' .and. moper.ne.'xti' .and.
     +       moper.ne.'fti' .and. moper.ne.'lti' .and.
     +       moper.ne.'num' .and. moper.ne.'sum' .and.
     +       moper.ne.'bel' .and. moper.ne.'abo' .and.
     +       moper.ne.'con' ) then
            write(0,*) 'daily2longer: error: unknown operation ',moper
            call abort
        endif
        if ( moper.eq.'max' .or. moper.eq.'min' .or. moper(2:3).eq.'ti'
     +       .or. moper.eq.'num' ) lstandardunits = .false.
        call getarg(1,file)
        call readseries(file,olddata,mpermax,yrbeg,yrend,nperyear,var
     +       ,units,lstandardunits,lwrite)
        if ( index(file,'/dd').gt.0 .or. index(file,'dd').eq.1 ) then
            write(0,*) 'Hmm. This looks like a (wind) direction to me.'
            write(0,*) 'Averaging over a unit circle.<p>'
            itype = 360
!
!           only uncomment if you REALLY need it. It should work,
!        else if ( units.eq.'dy' .or. units.eq.'day' ) then
!            write(0,*) 'Hmm. This looks like a day of the year to me.'
!            write(0,*) 'Averaging over the seasonal cycle.<p>'
!           itype = nperyear
!        else if ( units.eq.'mo' .or. units.eq.'month' ) then
!            write(0,*) 'Hmm. This looks like a month of the year to me.'
!            write(0,*) 'Averaging over the seasonal cycle.<p>'
!            itype = nperyear
        else
            itype = 0
        endif
*
*       read operation
*
        call getarg(2,string)
        read(string,*,err=901) nperyearnew
        if ( abs(nperyearnew).gt.npermax ) then
            write(0,*) 'daily2longer: error: nperyearnew = ',nperyearnew
     +           ,' not yet supported'
            write(*,*) 'daily2longer: error: nperyearnew = ',nperyearnew
     +           ,' not yet supported'
            call abort
        endif
        call getopts(4,iargc(),nperyear,yrbeg,yrend,.true.,0,0)
        if ( minfac.lt.0 ) minfac = 0.5
        if ( minindx.gt.-1e33 .or. pminindx.ge.0 ) then
            if ( maxindx.lt.1e33 .or. pmaxindx.ge.0 ) then
                write(0,*) 'daily2longer: error: unknown comparison '
     +               ,minindx,pminindx,maxindx,pmaxindx
            else
                lgt = '>'
                pcut = pminindx
                cut = minindx
                if ( lwrite ) print *,'using minval
     +               ',minindx,pminindx,'%'
            end if
        else
            if ( maxindx.lt.1e33 .or. pmaxindx.ge.0 ) then
                lgt = '<'
                pcut = pmaxindx
                cut = maxindx
                if ( lwrite ) print *,'using maxval
     +               ',maxindx,pmaxindx,'%'
            else
                lgt = ' '
            endif
        end if
        if ( lsum.gt.1 ) then
            print '(a,i3,a)','# ',lsum,' running average'
            call sumit(olddata,mpermax,nperyear,yrbeg,yrend,lsum,'v')
        end if
        if ( pminindx.eq.19712000 .or. pmaxindx.eq.19712000 ) then
            if ( lwrite ) print *,'take normals wrt 1971-2000'
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
*           no smoothing for the time being
            do yr=yrbeg,yrend
                do j=1,nperyear
                    if ( olddata(j,yr).lt.1e33 .and. 
     +                   normdata(j).lt.1e33 ) then
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
            if ( pcut.ge.0 ) then
                do j=1,nperyear
                    call getcutoff(cut(j),pcut,olddata,mpermax
     +                   ,nperyear,yrbeg,yrend,yrbeg,yrend,j,j,0)
                enddo
            endif
        endif
*
*       compute climatology and anomalies
*
        if ( ( moper.eq.'mea' .or. moper.eq.'sum' )
     +       .and. lgt.eq.' ' .and. .not.lnomissing ) then
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
!           construct climatologies for non-leap years
            if ( 366*(nperyear/366).eq.nperyear ) then
                n = nperyear/366
                do j=1,n*59
                    oldclim2(j) = oldclim(j)
                end do
                do j=n*60-(n-1),n*365
                    oldclim2(j) = oldclim(j+1)
                end do
                nperyear2 = n*365
            end if
            climunits = units
!           compute new climatology
            call allday2period(
     +           oldclim,mpermax,nperyear,
     +           newclim,npermax,nperyearnew,
     +           0,0,moper,lgt,cut,minfac,itype,var,climunits,lwrite)
            if ( 366*(nperyear/366) .eq. nperyear ) then
                call allday2period(
     +               oldclim2,mpermax,nperyear2,
     +               newclim2,npermax,nperyearnew,
     +               0,0,moper,lgt,cut,minfac,itype,var,climunits,lwrite
     +               )
            end if
!
!           fill in missing data when requested
!
            call fillmissingdata(olddata,mpermax,yrbeg,yrend,nperyear,
     +           add_option,.true.,lwrite)
        else
            do j=1,abs(nperyearnew)
                newclim(j) = 0
                newclim2(j) = 0
            enddo
        endif
*
*       perform operation
*
        call makeabsent(newdata,npermax,yrbeg,yrend)
        call allday2period(
     +     olddata,mpermax,nperyear,
     +     newdata,npermax,nperyearnew,
     +     yrbeg,yrend,moper,lgt,cut,minfac,itype,var,units,lwrite)
        do yr=yrbeg,yrend
            do j=1,abs(nperyearnew)
                if ( newdata(j,yr).lt.1e33 ) then
                    if ( (366*(nperyear/366) .eq. nperyear) 
     +                   .and. leap(yr).eq.1 ) then
                        newdata(j,yr) = newdata(j,yr) + newclim2(j)
                    else
                        newdata(j,yr) = newdata(j,yr) + newclim(j)
                    end if
                endif
            enddo
        enddo
*
*       print output
*
        if ( lgt.eq.' ' ) then
            print '(a,i4,2a)','# daily2longer ',nperyearnew,' ',moper
        elseif ( abs(pcut).lt.1 ) then
            print '(a,i4,3a,g20.4,a)','# daily2longer ',nperyearnew,' ',
     +          moper,lgt,pcut,'%'
        else
            print '(a,i4,3a,g20.4)','# daily2longer ',nperyearnew,' ',
     +          moper,lgt,cut(1)
        endif
        open(1,file=file)
        do
            read(1,'(a)') file
            if ( file(1:2).eq.' #' ) file = file(2:)
            if ( file(1:1).ne.'#' ) exit
            call adjustunitstring(file,units,var)
            print '(a)',trim(file)
        enddo
        close(1)
        call printdatfile(6,newdata,npermax,abs(nperyearnew),
     +       yrbeg,yrend)
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

        subroutine adjustunitstring(string,units,var)
        implicit none
        character string*(*),units*(*),var*(*)
        integer i,j
        i = index(string,'[')
        if ( i.ne.0 ) then
            j = index(string,']')
            if ( j.ne.0 ) then
                string = '# '//trim(var)//' ['//trim(units)//string(j:)
            endif
        endif
        end
