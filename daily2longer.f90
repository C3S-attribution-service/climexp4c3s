program daily2longer

!   compute aggregate quantities from daily data
!   input: daily time series
!   output: yearly/monthly/10-dy/5-dy time series

    implicit none
    include 'param.inc'
    include 'getopts.inc'
    integer :: mpermax
    parameter(mpermax=24*366)
    integer :: nperyear,nperyearnew,yr,mo,dy,i,j,n,itype,nperyear2
    integer*2 :: nn(mpermax)
    real :: olddata(mpermax,yrbeg:yrend),newdata(npermax,yrbeg:yrend), &
        cut(mpermax),pcut,normdata(npermax),s,oldclim(mpermax), &
        oldclim2(mpermax),newclim(npermax),newclim2(mpermax), &
        refs(yrbeg:yrend)
    logical :: lvalidclim(mpermax),lvalidclim2(mpermax), &
        lvalid(mpermax,yrbeg:yrend)
    character file*1023,string*512,lgt*1,moper*3,var*20,units*20, &
        climunits*10,lvar*120,svar*120,history*50000,metadata(2,100)*2000
    integer,external :: leap

    lwrite = .false. 
    lstandardunits = .true. 
    lnomissing = .false. 

    if ( command_argument_count() < 3 ) then
        print *,'usage: daily2longer infile nperyearnew'// &
            ' mean|sd|sum|abo|bel|num|min|max|mintime|maxtime'// &
            '|firsttime|lasttime|con [<> val[%|p]'// &
            ' add_anom|add_clim|add_trend|add_persist|add_damped]'
        print *,'(more options will come as requested)'
        call exit(-1)
    endif

!   read data

    call get_command_argument(3,string)
    if ( string == 'mintime' ) string = 'nti'
    if ( string == 'maxtime' ) string = 'xti'
    if ( string == 'firsttime' ) string = 'fti'
    if ( string == 'lasttime' ) string = 'lti'
    moper = string
    if ( moper /= 'mea' .and. moper /= 'sd ' .and. &
         moper /= 'min' .and. moper /= 'max' .and. &
         moper /= 'nti' .and. moper /= 'xti' .and. &
         moper /= 'fti' .and. moper /= 'lti' .and. &
         moper /= 'num' .and. moper /= 'sum' .and. &
         moper /= 'bel' .and. moper /= 'abo' .and. &
         moper /= 'con' ) then
        write(0,*) 'daily2longer: error: unknown operation ',moper
        call exit(-1)
    endif
    if ( moper == 'max' .or. moper == 'min' .or. moper(2:3) == 'ti' &
        .or. moper == 'num' ) lstandardunits = .false. 
    call get_command_argument(1,file)
    call readseriesmeta(file,olddata,mpermax,yrbeg,yrend,nperyear,var &
        ,units,lvar,svar,history,metadata,lstandardunits,lwrite)
    if ( index(file,'/dd') > 0 .or. index(file,'dd') == 1 ) then
        write(0,*) 'Hmm. This looks like a (wind) direction to me.'
        write(0,*) 'Averaging over a unit circle.<p>'
        itype = 360
    
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

!   read operation

    call get_command_argument(2,string)
    read(string,*,err=901) nperyearnew
    if ( abs(nperyearnew) > npermax ) then
        write(0,*) 'daily2longer: error: nperyearnew = ',nperyearnew,' not yet supported'
        write(*,*) 'daily2longer: error: nperyearnew = ',nperyearnew,' not yet supported'
        call exit(-1)
    endif
    call getopts(4,command_argument_count(),nperyear,yrbeg,yrend,.true.,0,0)
    if ( minfac < 0 ) minfac = 0.5
    pcut = -999.9
    if ( minindx > -1e33 .or. pminindx >= 0 ) then
        if ( maxindx < 1e33 .or. pmaxindx >= 0 ) then
            write(0,*) 'daily2longer: error: unknown comparison ' &
                ,minindx,pminindx,maxindx,pmaxindx
        else
            lgt = '>'
            pcut = pminindx
            cut = minindx
            if ( lwrite ) print *,'using minval ',minindx,pminindx,'%'
        end if
    else
        if ( maxindx < 1e33 .or. pmaxindx >= 0 ) then
            lgt = '<'
            pcut = pmaxindx
            cut = maxindx
            if ( lwrite ) print *,'using maxval ',maxindx,pmaxindx,'%'
        else
            lgt = ' '
        endif
    end if
    if ( lsum > 1 ) then
        print '(a,i3,a)','# ',lsum,' running average'
        call sumit(olddata,mpermax,nperyear,yrbeg,yrend,lsum,'v')
    end if
    if ( pminindx == 19712000 .or. pmaxindx == 19712000 ) then
        if ( lwrite ) print *,'take normals wrt 1971-2000'
        do j=1,nperyear
            normdata(j) = 0
            n = 0
            do yr=1971,2000
                if ( olddata(j,yr) < 1e33 ) then
                    n = n + 1
                    normdata(j) = normdata(j) + olddata(j,yr)
                endif
            enddo
            if ( n > 5 ) then ! arbitrary number
                normdata(j) = normdata(j)/n
            else
                normdata(j) = 3e33
            endif
            if ( lwrite ) print *,j,normdata(j)
        enddo
!       no smoothing for the time being
        do yr=yrbeg,yrend
            do j=1,nperyear
                if ( olddata(j,yr) < 1e33 .and. normdata(j) < 1e33 ) then
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
                    if ( olddata(j,yr) < 1e33 ) then
                        n = n + 1
                        s = s + olddata(j,yr)
                    endif
                enddo
                if ( n > 0 ) then
                    print *,j,s/n,n
                endif
            enddo
        endif
        lgt = ' '
    else
        if ( pcut >= 0 ) then
            do j=1,nperyear
                call getcutoff(cut(j),pcut,olddata,mpermax &
                    ,nperyear,yrbeg,yrend,yrbeg,yrend,j,j,0)
            enddo
        endif
    endif

!   compute climatology and anomalies

!   This test should be exactly the same as in day2period
    if ( ( moper == 'mea' .or. moper == 'sum' ) .and. lgt == ' ' ) then
        do j=1,nperyear
            oldclim(j) = 0
        enddo
        do j=1,nperyear
            nn(j) = 0
        enddo
        do yr=yrbeg,yrend
            do j=1,nperyear
                if ( olddata(j,yr) < 1e33 ) then
                    nn(j) = nn(j) + 1
                    oldclim(j) = oldclim(j) + olddata(j,yr)
                endif
            enddo
        enddo
        do j=1,nperyear
            if ( nn(j) > 0 ) then
                oldclim(j) = oldclim(j)/nn(j)
            else
                oldclim(j) = 3e33
            endif
        enddo
        do yr=yrbeg,yrend
            do j=1,nperyear
                if ( olddata(j,yr) < 1e33 ) then
                    olddata(j,yr) = olddata(j,yr) - oldclim(j)
                endif
            enddo
        enddo
!       construct climatologies for non-leap years
        if ( 366*(nperyear/366) == nperyear ) then
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
!       compute new climatology
        lvalidclim = .true. 
        lvalidclim2 = .true. 
        call allday2period( &
        oldclim,mpermax,nperyear,lvalidclim, &
            newclim,npermax,nperyearnew, &
            0,0,moper,lgt,cut,minfac,itype,var,climunits,lwrite)
        if ( 366*(nperyear/366) == nperyear ) then
            call allday2period( &
                oldclim2,mpermax,nperyear2,lvalidclim2, &
                newclim2,npermax,nperyearnew, &
                0,0,moper,lgt,cut,minfac,itype,var,climunits,lwrite)
        end if
    
!       fill in missing data when requested
    
        call fillmissingdata(olddata,lvalid,refs,mpermax,yrbeg,yrend &
            ,nperyear,add_option, .true. ,lwrite)
    else
        do j=1,abs(nperyearnew)
            newclim(j) = 0
            newclim2(j) = 0
        enddo
    endif

!   perform operation

    call makeabsent(newdata,npermax,yrbeg,yrend)
    call allday2period( &
        olddata,mpermax,nperyear,lvalid, &
        newdata,npermax,nperyearnew, &
        yrbeg,yrend,moper,lgt,cut,minfac,itype,var,units,lwrite)
    do yr=yrbeg,yrend
        do j=1,abs(nperyearnew)
            if ( newdata(j,yr) < 1e33 ) then
                if ( (366*(nperyear/366) == nperyear) .and. leap(yr) == 1 ) then
                    newdata(j,yr) = newdata(j,yr) + newclim2(j)
                else
                    newdata(j,yr) = newdata(j,yr) + newclim(j)
                end if
            endif
        enddo
    enddo

!   print output

    call adjustlvar(moper,lvar,nperyearnew,lwrite)
    call printvar(6,var,units,lvar)
    call printmetadata(6,file,' ',' ',history,metadata)
    call printdatfile(6,newdata,npermax,abs(nperyearnew),yrbeg,yrend)

!   error messages

    goto 999
901 write(0,*) 'daily2longer: expecting nperyearnew, not ',string
    call exit(-1)
902 write(0,*) 'daily2longer: expecting value[%|p], not ',string
    call exit(-1)
999 continue
end program daily2longer