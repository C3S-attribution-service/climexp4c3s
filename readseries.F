        subroutine readseries(file,data,npermax,yrbeg,yrend,nperyear,
     +       var,units,lstandardunits,lwrite)
*
*       file (in) filename to read.
*       this routines accepts 3 formats:
*       - a netCDF file with one time series
*       - an ascii file with lines with 
*           'yyyy val_jan val_feb ... val_dec [val_year]'
*       - an ascii file with lines with year.fractionofyear value
*       - an ascii file with lines with year month value
*       - an ascii file with lines with year month day value
*       Comment line start with '#'
*       data (out) the data is returned in data, which is filled with 
*       undef's where no data exists.
*       nperyear,yrbeg,yrend (out) dimensions of data
*       units,var (out) properties of file
*
*       I try to convert the units to my standard ones:
*       temp [C], prcp [mm/period], slp [hPa], ...
*
        implicit none
#include "netcdf.inc"
        integer npermax,nperyear,yrbeg,yrend
        real data(npermax,yrbeg:yrend)
        character file*(*),units*(*),var*(*)
        logical lstandardunits,lwrite
        integer i,j,k,m,year,month,day,unit,status,ncid,dpm(12),
     +       dpm0(12),imonth,year1
        logical lagain
        real*8 x,y,x1,y1,val12(366)
        character line*10000
        integer getnumwords,monthinline
        save dpm,dpm0        
        data dpm0 /31,29,31,30,31,30,31,31,30,31,30,31/
        do i=1,12
            dpm(i) = dpm0(i)
        enddo
*
        if ( lwrite ) then
            print *,'readseries: opening file ',trim(file)
            print *,'npermax,yrbeg,yrend = ',npermax,yrbeg,yrend
        endif
        call makeabsent(data,npermax,yrbeg,yrend)
*
*       also handle netCDF time series
*
        if ( file.eq.'-' ) then
            unit = 5
        else
            status = nf_open(file,nf_nowrite,ncid)
            if ( status.ne.nf_noerr ) then
                call rsunit(unit)
                open(unit,file=trim(file),status='old',err=902)
            else
*               it is a netCDF file
                if ( lwrite ) print *,'calling readncseries ',trim(file)
                call readncseries(file,data,npermax,nperyear,yrbeg,yrend
     +               ,ncid,var,units,lwrite)
                return
            endif
        endif
*       it is an ascii file.
        call readseriesheader(var,units,line,unit,lwrite)
        if ( line.eq.' ' ) then
            nperyear = 0        ! is undefined when there is no data
            if ( unit.ne.5 ) close(unit)
            return
        endif
*       first look for explicit months
        imonth = monthinline(line)
        nperyear = getnumwords(line) - 1
        if ( nperyear.eq.4 ) then
!           it is either seasonal data: yr valDJF valMAM valJJA valSON
!           or hourly data: yr mo dy hr val
            read(line,*) year
            read(unit,*) year1
            if ( year.eq.year1 ) then
                if ( lwrite ) print *,'hourly data'
                nperyear = 366*24
            else
                if ( lwrite ) print *,'seasonal data'
            endif
            backspace(unit)
        endif
*       often there is a yearly mean/sum at the end
        if ( nperyear.eq.13 ) nperyear = 12
        if ( nperyear.gt.npermax ) then
            write(0,*) 'readseries: error: increase npermax ',npermax
     +            ,nperyear
            write(*,*) 'readseries: error: increase npermax ',npermax
     +            ,nperyear
            call abort
        endif
        if ( imonth.gt.0 ) then
            call readdymoyrval(data,npermax,yrbeg,yrend,nperyear,line
     +           ,unit,lwrite)
        elseif ( nperyear.eq.1 ) then
            call readyrfracval(data,npermax,yrbeg,yrend,nperyear
     +           ,line,unit,lwrite)
        elseif ( nperyear.eq.2 ) then
            call readyrmoval(data,npermax,yrbeg,yrend,nperyear,line,unit
     +           ,lwrite)
        elseif ( nperyear.eq.3 ) then
            call readyrmodyval(data,npermax,yrbeg,yrend,nperyear,line
     +           ,unit,lwrite)
        elseif ( nperyear.eq.366*24 ) then
            call readyrmodyhrval(data,npermax,yrbeg,yrend,nperyear,line
     +           ,unit,lwrite)
        else
            call readyrval(data,npermax,yrbeg,yrend,nperyear,line,unit
     +           ,lwrite)
        endif
        if ( unit.ne.5 ) then
            if ( lwrite ) print *,'closing unit ',unit
            close(unit)
        end if
        if ( lstandardunits ) then
            call makestandardseries(data,npermax,yrbeg,yrend,nperyear
     +           ,var,units,lwrite)
        endif
        return
 902    print *,'readseries: error opening ',trim(file)
        call abort
        end

        subroutine readseriesheader(var,units,line,unit,lwrite)
*       read header, try to find the variable name and the units
*       it leaves the last (non-header) line read in line
        implicit none
        integer unit
        logical lwrite
        character var*(*),units*(*),line*(*)
        integer i,j,k,iline
        logical onlynumbers
        units = ' '
        var = ' '
        iline = 0
 100    continue
        read(unit,'(a)',err=903,end=900) line
        if ( lwrite ) print *,'readseriesheader: processing line '
     +       ,trim(line)
*       the Mac sometimes gets confused and puts something else before
*       the '#'
        if ( line(1:1).eq.'#' .or. line(2:2).eq.'#' .or.
     +       line.eq.' ' ) then
*           do not count lines, but look for '#'
            if ( lwrite ) print *,'found new-style #-comments'
            iline = -1
        else if ( iline.eq.-1 ) then
            if ( lwrite ) print *,'end of #-comments - return'
            return
        else if ( onlynumbers(line) ) then
            if ( lwrite ) print *,'line with only numbers - return'
            return
        else
            iline = iline + 1
        endif
        j = index(line,'[')
        if ( j.gt.0 ) then
            k = index(line,']')
            if ( k.gt.j ) then
                units = line(j+1:k-1)
                call checkstring(units)
                j = 1 + index(line,';',.true.)
                ! ugly hack
                if ( j.gt.3 ) then
                    if ( index(line(max(1,j-8):j-1),'&').ne.0 ) then
                        j = 1   ! not a break but an esacped HTML character code...
                    end if
                end if
 110            continue
                if ( line(j:j).eq.' ' ) then
                    j = j + 1
                    goto 110
                endif
                if ( line(j:j).eq.'#' ) then
                    j = j + 2
                endif
                k = j + index(line(j:),' ') - 2
                var = line(j:k)
                if ( var.eq.'GDCN' ) then
*                   hard-coded...
                    j = 13
                    k = j + index(line(j:),' ') - 2
                    var = line(j:k)
                    call tolower(var)
                endif
                call checkstring(var)
            endif
            if ( lwrite ) then
                print *,'found var  =',trim(var)
                print *,'found units=',trim(units)
            endif
        endif
        if ( iline.lt.5 ) goto 100
        read(unit,'(a)') line
        if ( lwrite ) print *,'First line with data ',trim(line)
        return
 900    line = ' '
        return
 903    print *,'readseriesheader: error reading header'
        print *,line
        call abort
        end

        logical function onlynumbers(line)
        implicit none
        character*(*) line
        integer i
        if ( line.eq.' ' ) then
            onlynumbers = .false.
            return
        endif
        onlynumbers = .true.
        do i=1,4
            if ( lge(line(i:i),'0') .and. lle(line(i:i),'9') .or.
     +           line(i:i).eq.' ' .or. line(i:i).eq.'\t' .or.
     +           line(i:i).eq.'.' .or. line(i:i).eq.'-' .or. 
     +           line(i:i).eq.'\n' .or. line(i:i).eq.'\r' ) then
                cycle
            else
                onlynumbers = .false.
                return
            endif
        end do
        do i=5,len(line)
            if ( lge(line(i:i),'0') .and. lle(line(i:i),'9') .or.
     +           line(i:i).eq.' ' .or. line(i:i).eq.'.' .or.
     +           line(i:i).eq.'-' .or. line(i:i).eq.'+' .or.
     +           line(i:i).eq.'e' .or. line(i:i).eq.'E' .or.
     +           line(i:i).eq.'d' .or. line(i:i).eq.'D' .or.
     +           line(i:i).eq.'\t' .or.
     +           line(i:i).eq.'\n' .or. line(i:i).eq.'\r' ) then
                cycle
            else
                onlynumbers = .false.
                return
            endif
        end do
        end

        integer function monthinline(line)
        implicit none
        character line*(*)
        integer imonth
        call tolower(line)
        imonth = index(line,'jan')
        if ( imonth.eq.0 ) imonth = index(line,'feb')
        if ( imonth.eq.0 ) imonth = index(line,'mar')
        if ( imonth.eq.0 ) imonth = index(line,'apr')
        if ( imonth.eq.0 ) imonth = index(line,'may')
        if ( imonth.eq.0 ) imonth = index(line,'jun')
        if ( imonth.eq.0 ) imonth = index(line,'jul')
        if ( imonth.eq.0 ) imonth = index(line,'aug')
        if ( imonth.eq.0 ) imonth = index(line,'sep')
        if ( imonth.eq.0 ) imonth = index(line,'oct')
        if ( imonth.eq.0 ) imonth = index(line,'nov')
        if ( imonth.eq.0 ) imonth = index(line,'dec')
        monthinline = imonth
        end

        subroutine readdymoyrval(data,npermax,yrbeg,yrend,nperyear,line
     +       ,unit,lwrite)
*
*       format is NN-MON-YYYY, possibly with the dashes.
*       Assume daily data in a Gregorian calendar
*       TODO: other calendars, other nperyears.
*
        implicit none
        integer npermax,yrbeg,yrend,nperyear,unit
        real data(npermax,yrbeg:yrend)
        logical lwrite
        character line*(*)
        integer imonth,j,m,day,month,year,dpm(12),dpm0(12)
        real*8 val
        integer monthinline
        save dpm,dpm0        
        data dpm0 /31,29,31,30,31,30,31,31,30,31,30,31/

        if ( lwrite ) print *,'readdymoyrval: reading data as '//
     +       'dd-mon-yyyy or ddmonyyyy'
        do month=1,12
            dpm(month) = dpm0(month)
        enddo
        nperyear = 366
 10     continue
        imonth = monthinline(line)
        if ( imonth.le.0 ) then
            write(0,'(a,i3,a)')
     +           'readseries: cannot recognize month in ',imonth
     +           ,line(imonth:imonth+2)
            write(0,'(a)') trim(line)
            goto 20
        endif
        month = (index
     +       ('jan feb mar apr may jun jul aug sep oct nov dec'
     +       ,line(imonth:imonth+2)) + 3)/4
        j = imonth-1
        if ( line(j:j).eq.'-' ) j = j-1
        if ( j-1.le.0 ) goto 11
        read(line(j-1:j),*,end=1800,err=11) day
        if ( day.lt.1 .or. day.gt.dpm(month) ) goto 11
        goto 12
 11     continue
        write(0,'(a,i8,a,i3,a)') 'readseries: strange day ',day
     +       ,' in ',j,line(j-1:j)
        write(0,'(a)') trim(line)
        goto 20
 12     continue
        j = imonth+3
        if ( line(j:j).eq.'-' ) j = j + 1
        read(line(j:),*,err=13) year
        goto 14
 13     continue
        write(0,'(a,i8,a)') 'readseries: error reading year in'
        write(0,'(a)') trim(line)
        goto 20
 14     continue
        do m=1,month-1
            day = day + dpm(m)
        enddo
        read(line(imonth+8:),*) val
        call checkvalid(year,day,val)
        if ( val.lt.1e33 ) data(day,year) = val
***     print *,day,year,data(day,year)
 20     continue
        read(unit,'(a)',err=904,end=1800) line
        call tolower(line)
        goto 10
 1800   continue
        return
  904   print *,'readdymoyrval: error reading data of line'
        print '(a)',trim(line)
        print '(a,2i3,g14.6)','last data read ',year,day,val
        call abort
        end

        subroutine readyrval(data,npermax,yrbeg,yrend,nperyear,line,unit
     +       ,lwrite)
*
*       format is year val(1) val(2) ...val(nperyear) (nperyear >= 4)
*
        implicit none
        integer npermax,yrbeg,yrend,nperyear,unit
        real data(npermax,yrbeg:yrend)
        logical lwrite
        character line*(*)
        integer year,j
        real*8,allocatable :: val12(:)
        logical lfirst
        save lfirst
        data lfirst /.true./
*
        if ( lwrite ) then
            print *,'readyrval: reading data as '//
     +           'yyyy val1 val2 ... valN (N>3)'
            print *,'readyrval: npermax,yrbeg,yrend,nperyear = ',npermax
     +           ,yrbeg,yrend,nperyear
        endif
        allocate(val12(npermax))
        read(line,*,err=904,end=300) year,(val12(j),j=1,nperyear)
 100    continue
        if ( .false. .and. lwrite ) print *,'read year,val12 = ',year
     +       ,(val12(j),j=1,nperyear)
        if ( year.lt.yrbeg .or. year.gt.yrend ) then
            if ( lfirst ) then
                lfirst = .false.
                write(0,'(a,3i5)') '# disregarding year ',year,yrbeg
     +               ,yrend
            endif
        else
            do j=1,nperyear
                call checkvalid(year,j,val12(j))
                if ( val12(j).lt.1e33 ) data(j,year) = val12(j)
            enddo
        endif
        read(unit,*,err=904,end=300) year,(val12(j),j=1,nperyear)
        goto 100
 300    continue
        deallocate(val12)
        return
 904    print *,'readyrval: error reading data of line'
        print '(a)',trim(line)
        print '(a,i4,1000g14.6)','last data read ',year,(val12(j),j=1
     +       ,nperyear)
        call abort
        end

        subroutine readyrfracval(data,npermax,yrbeg,yrend,nperyear,line
     +       ,unit,lwrite)
*
*       format is year.frac value or yyyy[mm[dd[hh]]] value
*       Figure out nperyear from the data...
*
        implicit none
        integer npermax,yrbeg,yrend,nperyear,unit
        real data(npermax,yrbeg:yrend)
        logical lwrite
        character line*(*)
        integer i,j,k,n,year,month,day,hour,year1,month1,day1,hour1,ymdh
     +       ,dpm(12),nperyear1,jul,jul1
        integer,external :: julday,leap,getj
        save dpm
        real*8 x,y,x1,y1
        logical lfrac,lfirst
        logical onlynumbers
        save lfirst
        data dpm /31,29,31,30,31,30,31,31,30,31,30,31/
        data lfirst /.true./
*
        i = scan(line,'123456789')
        if ( i.eq.0 ) then
            write(0,*)
     +           'readyrfracval: error: cannot find 1 or 2 in line '
     +           ,trim(line)
            call abort
        end if
        j = scan(line(i:),' ')
        if ( index(line(i:i+j),'.').ne.0 ) then
            lfrac = .true.
            if ( lwrite ) print *,'readyrfracval: reading data as '//
     +           'yyyy.frac val'
        else
            lfrac = .false.
            if ( lwrite ) print *,'readyrfracval: reading data as '//
     +           'yyyy[mm[dd[hh]]] val'
        endif
        nperyear = 0
 401    continue
        call tolower(line)
        if ( lfrac ) then
            read(line,*,err=904,end=500) x,y
            year = int(x+1e-4)
        else
            read(line,*,err=904,end=500) ymdh,y
            call yyyymmddhh(ymdh,year,month,day,hour)
            if ( lwrite ) print *,'read and interpreted ',year,month,day
     +           ,hour,y
        endif
        call checkvalid(year,1,y)
        if ( y.gt.1e33 ) then
            if ( lwrite ) print *,'invalid data, try again'
 4015       continue            
            read(unit,'(a)',err=904,end=1800) line            
            if ( .not.onlynumbers(line) ) then
                if ( lwrite ) print *,'found line with non-numbers '
     +               ,trim(line)
                goto 4015
            endif
            if ( lwrite ) print *,'really try again'
            goto 401
        endif
        i = 0
 402    continue
        i = i + 1
 405    continue
        read(unit,'(a)',err=904,end=500) line
        if ( lwrite ) print *,'read line ',trim(line)
        if ( .not.onlynumbers(line) ) then
            if ( lwrite ) print *,'Found non-numbers in ',trim(line)
            goto 405
        end if
        if ( lfrac ) then
            read(line,*,err=904,end=500) x1,y1
            year1 = int(x+1e-4)
        else
            read(line,*,err=904,end=500) ymdh,y1
            call yyyymmddhh(ymdh,year1,month1,day1,hour1)
            if ( lwrite ) print *,'read and interpreted ',year1,month1
     +           ,day1,hour1,y1
            if ( day.gt.dpm(month) .or.
     +           day.gt.28 .and. month.eq.2 .and. leap(year).eq.1 ) then
                write(0,'(a,2i3,i5)')
     +               '# readseries: disregarded invalid date ',day,month
     +               ,year
                goto 402
            end if
        endif
        call checkvalid(year1,2,y1)
        if ( y1.gt.1e33 ) then
            if ( lwrite ) print *,'invalid data, try again'
            goto 402
        endif
        if ( lfrac ) then
            nperyear = max(nperyear,nint(i/(x1-x)))
        else
            jul = julday(month,day,year)
            jul1 = julday(month1,day1,year1)
            if ( jul.eq.jul1 .and. hour.eq.hour1 ) then
                write(0,*) 'readyrfracval: duplicate date 1 ',hour,day
     +               ,month,year,hour1,day1,month1,year1
                goto 402
            end if
            nperyear = max(nperyear,
     +           nint((24*366.)/(24*(jul1-jul) + (hour1-hour))))
        end if
!       round-off errors...
        call adjustnperyear(nperyear,lwrite)
        if ( lwrite ) then
            if ( lfrac ) then
                write(*,*) 'readyrfracval: found nperyear = ',
     +               nperyear,i,x1,x,i/(x1-x)
            else
                write(*,*) 'readyrfracval: found nperyear = ',
     +               nperyear
                write(*,*) year1,year
                write(*,*) month1,month
                write(*,*) day1,day
                write(*,*) hour1,hour
            endif
        endif
        if ( nperyear.gt.npermax ) then
            write(0,*) 'readseries: error: increase npermax ',npermax
     +           ,nperyear
            write(*,*) 'readseries: error: increase npermax ',npermax
     +           ,nperyear
            call abort
        endif
        if ( lfrac ) then
            j = nint(nperyear*(x-year)+0.75001)
        else
            j = getj(month,day,hour,nperyear)
        endif
        call checkvalid(year,j,y)
        if ( year.lt.yrbeg .or. year.gt.yrend ) then
            if ( lfirst ) then
                lfirst = .false.
                write(0,'(a,i5)') '# disregarding year ',year
            endif
        elseif ( y.lt.1e33 ) then
            if ( lwrite ) then
                print *,'data(',year,j,') = ',y
            endif
            data(j,year) = y
        endif
 450    continue
        if ( lfrac ) then
            year = int(x1+1e-4)
            j = nint(nperyear*(x1-year)+0.75001)
        else
            j = getj(month1,day1,hour1,nperyear)
            year = year1
        endif
        call checkvalid(year,j,y1)
        if ( year.lt.yrbeg .or. year.gt.yrend ) then
            if ( lfirst ) then
                lfirst = .false.
                write(0,'(a,i5)') '# disregarding year ',year
            endif
        elseif ( y1.lt.1e33 ) then
            if ( lwrite ) then
                print *,'data(',year,j,') = ',y1
            endif
            data(j,year) = y1
        endif
        x = x1
        y = y1
        year = year1
        month = month1
        day = day1
        hour = hour1
 460    continue
        read(unit,'(a)',err=904,end=500) line
        if ( .not.onlynumbers(line) ) goto 460        
        if ( lfrac ) then
            read(line,*,err=904,end=500) x1,y1
            nperyear1 = nint(1/(x1-x))
            call adjustnperyear(nperyear1,lwrite)
        else
            read(line,*,err=904,end=500) ymdh,y1
            call yyyymmddhh(ymdh,year1,month1,day1,hour1)
            if ( day1.gt.dpm(month1) .or.
     +           day1.gt.28 .and. month1.eq.2 .and. leap(year1).eq.1 )
     +           then
                write(0,'(a,2i3,i5)')
     +               '# readseries: disregarded invalid date ',day1
     +               ,month1,year1
                goto 460
            end if
            jul = julday(month,day,year)
            jul1 = julday(month1,day1,year1)
            if ( jul.eq.jul1 .and. hour.eq.hour1 ) then
                write(0,*) 'readyrfracval: duplicate date 2 ',hour,day
     +               ,month,year,hour1,day1,month1,year1
                goto 460
            end if            
            nperyear1 = nint((24*366.)/(24*(jul1-jul) + (hour1-hour)))
            if ( lwrite .and. nperyear1.ne.nperyear ) then
                print *,'nperyear1 = ',nperyear1
                print *,'date  = ',year,month,day,hour
                print *,'date1 = ',year1,month1,day1,hour1
            end if
            call adjustnperyear(nperyear1,lwrite)
        endif
        if ( nperyear1.gt.nperyear ) then
            if ( lwrite ) print *,'adjusting nperyear from ',nperyear
     +           ,' to ',nperyear1
            nperyear = nperyear1
!           and start all over again
            data = 3e33
            rewind(unit)
 480        continue
            read(unit,'(a)') line
            if ( .not.onlynumbers(line) ) goto 480
            if ( lwrite ) print *,'restarting with line = ',trim(line)
            goto 401
        endif
        goto 450
 500    continue
        return
 1800   continue
        nperyear = 1
        if ( y.lt.1e33 ) data(1,year) = y
        return
 904    print *,'readyrfracval: error reading data of line'
        print '(a)',trim(line)
        if ( lfrac ) then
            print '(a,4g14.6)','last data read ',x1,y1,x,y
        else
            print '(a,4i4,g14.6)','last data read ',year1,month1,day1
     +           ,hour1,y1
        endif
        call abort
        end

        subroutine readyrmoval(data,npermax,yrbeg,yrend,nperyear,line
     +       ,unit,lwrite)
*
*       format is year mon value
*       !!!Figure out nperyear from the data...
*       nperyear is set to 12 always
*
        implicit none
        integer npermax,yrbeg,yrend,nperyear,unit
        real data(npermax,yrbeg:yrend)
        logical lwrite
        character line*(*)
        integer year,j
        real*8 y
        logical lfirst
        save lfirst
        data lfirst /.true./

        if ( lwrite ) then
            print *,'readyrmoval: reading data as yyyy mm val'
            print *,'npermax,yrbeg,yrend = ',npermax,yrbeg,yrend
        endif
        nperyear = 12
        read(line,*,err=904,end=700) year,j,y
 600    continue
        call checkvalid(year,j,y)
        if ( year.lt.yrbeg .or. year.gt.yrend ) then
            if ( lfirst ) then
                lfirst = .false.
                write(0,'(a,i5)') '# disregarding year ',year
            endif
        elseif ( y.lt.1e33 ) then
***         print *,'!! data(',j,year,') = ',y
            data(j,year) = y
        endif
        read(unit,*,err=904,end=700) year,j,y
        goto 600
 700    continue
        return
 904    print *,'readyrmoval: error reading data of line'
        print '(a)',trim(line)
        print '(a,i4,1000g14.6)','last data read ',year,j,y
        call abort
        end

        subroutine readyrmodyval(data,npermax,yrbeg,yrend,nperyear,line
     +       ,unit,lwrite)
*
*       format is year mon day value
*
        implicit none
        integer npermax,yrbeg,yrend,nperyear,unit
        real data(npermax,yrbeg:yrend)
        logical lwrite
        character line*(*)
        integer i,j,year,month,day,dpm(12),dpm0(12)
        real y,xyr,xmo,xdy
        logical lfirst
        save dpm,dpm0,lfirst
        data dpm0 /31,29,31,30,31,30,31,31,30,31,30,31/
        data lfirst /.true./

        if ( lwrite ) print *,'readyrmoval: reading data as '//
     +       'yyyy mm dd val'
        do month=1,12
            dpm(month) = dpm0(month)
        enddo
        nperyear = 366          ! first try
        read(line,*,err=300,end=900) year,month,day,y
        goto 400
 300    continue
        read(line,*,err=905,end=900) xyr,xmo,xdy,y
        year = nint(xyr)
        month = nint(xmo)
        day = nint(xdy)
 400    continue
 800    continue
        if ( y.eq.-999.9 ) goto 899
        if ( month.lt.1 .or. month.gt.12 ) then
            print *,'invalid month ',month
            goto 899
        endif
        if ( day.lt.1 .or. day.gt.dpm(month) ) then
            if ( day.eq.30 .and. month.eq.2 .and. 
     +           data(31,year).gt.1e33) then
                nperyear = 360
                do i=1,12
                    dpm(i) = 30
                enddo
                do i=32,60
                    data(i-1,year) = data(i,year)
                enddo
            else
                write(0,*) 'invalid day ',day
                goto 899
            endif
        endif
        j = day
        do i=1,month-1
            j = j + dpm(i)
        enddo
        if ( year.lt.yrbeg .or. year.gt.yrend ) then
            if ( year.ne.-999 ) then
                if ( lfirst ) then
                    lfirst = .false.
                    write(0,'(a,i5)') '# disregarding year ',year
                endif
            endif
        elseif ( y.lt.1e33 ) then
***         print *,'data(',j,year,') = ',y,day,month
            data(j,year) = y
        endif
 899    continue
8990    continue
        read(unit,'(a)',err=905,end=900) line
        if ( line(1:1).eq.'#' .or. line.eq.' ' ) goto 8990
        read(line,*,err=915,end=900) year,month,day,y
        goto 800
 915    continue
        read(line,*,err=905,end=900) xyr,xmo,xdy,y
        year = nint(xyr)
        month = nint(xmo)
        day = nint(xdy)
        goto 800
 900    continue
        return
 905    print *,'readseries: error reading yr,mo,dy,val data'
        print *,trim(line)
        call abort
        end

        subroutine readyrmodyhrval(data,npermax,yrbeg,yrend,nperyear
     +       ,line,unit,lwrite)
*
*       format is year mon day hour value
*
        implicit none
        integer npermax,yrbeg,yrend,nperyear,unit
        real data(npermax,yrbeg:yrend)
        logical lwrite
        character line*(*)
        integer i,j,year,month,day,hour,dpm(12),dpm0(12)
        real y
        logical lfirst
        save dpm,dpm0,lfirst
        data dpm0 /31,29,31,30,31,30,31,31,30,31,30,31/
        data lfirst /.true./

        if ( lwrite ) print *,'readyrmoval: reading data as '//
     +       'yyyy mm dd hh val'
        do month=1,12
            dpm(month) = dpm0(month)
        enddo
        nperyear = 366*24
        read(line,*,err=905,end=900) year,month,day,hour,y
 800    continue
        if ( month.lt.1 .or. month.gt.12 ) then
            print *,'invalid month ',month
            goto 899
        endif
        if ( day.lt.1 .or. day.gt.dpm(month) ) then
            if ( day.eq.30 .and. month.eq.2 .and. 
     +           data(31,year).gt.1e33) then
                nperyear = 360
                do i=1,12
                    dpm(i) = 30
                enddo
                do i=32,60
                    data(i-1,year) = data(i,year)
                enddo
            else
                write(0,*) 'invalid day ',day
                goto 899
            endif
        endif
        if ( hour.lt.0 .or. hour.gt.24 ) then
            write(0,*) 'invalid hour ',hour
            goto 899
        endif

        j = day
        do i=1,month-1
            j = j + dpm(i)
        enddo
        j = 24*(j-1) + hour + 1
        if ( year.lt.yrbeg .or. year.gt.yrend ) then
            if ( year.ne.-999 ) then
                if ( lfirst ) then
                    lfirst = .false.
                    write(0,'(a,i5)') '# disregarding year ',year
                endif
            endif
        elseif ( y.lt.1e33 ) then
***         print *,'data(',j,year,') = ',y,hour,day,month
            data(j,year) = y
        endif
 899    continue
        read(unit,*,err=905,end=900) year,month,day,hour,y
        goto 800
 900    continue
        return
 905    print *,'readseries: error reading yr,mo,dy,hr,val data'
        print *,year,month,day,hour,y
        call abort
        end

        subroutine copyheader(file,unit)
*
*       copy the header of time series file file to unit unit
*
        implicit none
#include "netcdf.inc"
        integer mxmax,mymax,mzmax,nvmax
        parameter(mxmax=1,mymax=1,mzmax=1,nvmax=1)
        character file*(*)
        integer unit
        integer i,status,nx,ny,nz,nt,firstyr,firstmo,nvars,
     +       ivars(6,nvmax),ncid,nperyear,iunit
        real xx(mxmax),yy(mymax),zz(mzmax),undef
        character line*256,vars(nvmax)*40,lvars(nvmax)*120,
     +       units(nvmax)*80,title*1000
        logical onlynumbers
*
        status = nf_open(file,nf_nowrite,ncid)
        if ( status.ne.0 ) then
            call rsunit(iunit)
            open(iunit,file=file)
            read(iunit,'(a)') line
            if ( onlynumbers(line) ) then
!               no header
                close(iunit)
                return
            elseif ( line(1:1).ne.'#' .and. line(2:2).ne.'#' ) then
                write(unit,'(a)') trim(line)
                do i=2,5
                    read(iunit,'(a)') line
                    if ( onlynumbers(line) ) exit
                    write(unit,'(a)') trim(line)
                enddo
            else
                write(unit,'(a)') trim(line)
                do i=1,1000
                    read(iunit,'(a)') line
                    if ( line(1:1).ne.'#' ) exit
                    write(unit,'(a)') trim(line)
                enddo
            endif
            close(iunit)
        else                    ! netcdf file
            call parsenc(file,ncid,mxmax,nx,xx,mymax,ny,yy,mzmax
     +        ,nz,zz,nt,nperyear,firstyr,firstmo,undef,title,nvmax,nvars
     +        ,vars,ivars,lvars,units)
            status = nf_close(ncid)
            write(unit,'(2a)') '# ',trim(title)
            write(unit,'(6a)') '# ',trim(vars(1)),' [',trim(units(1))
     +           ,'] ',trim(lvars(1))
        end if
        end

        subroutine yyyymmddhh(ymdh,year,month,day,hour)
        implicit none
        integer ymdh,year,month,day,hour
        if ( ymdh.lt.9999 ) then
            year = ymdh
            month = 1
            day = 1
            hour = 0
        elseif ( ymdh.lt.999999 ) then
            year = ymdh/100
            month = mod(ymdh,100)
            day = 1
            hour = 0
        elseif ( ymdh.lt.99999999 ) then
            year = ymdh/10000
            month = mod(ymdh,10000)/100
            day = mod(ymdh,100)
            hour = 0
        else
            year = ymdh/1000000
            month = mod(ymdh,1000000)/10000
            day = mod(ymdh,10000)/100
            hour = mod(ymdh,100)
        endif
        if ( month.lt.1 .or. month.gt.12 ) then
            write(0,*) 'yyyymmddhh: error: found month ',month,' in ',
     +		ymdh
     	    call abort
     	endif
        if ( day.lt.1 .or. day.gt.31 ) then
            write(0,*) 'yyyymmddhh: error: found day ',day,' in ',
     +		ymdh
     	    call abort
     	endif
        if ( hour.lt.0 .or. month.gt.24 ) then
            write(0,*) 'yyyymmddhh: error: found hour ',hour,' in ',
     +		ymdh
     	    call abort
     	endif
        end

        integer function getj(month,day,hour,nperyear)
!
!       compute index in array given a date and nperyear
!
        implicit none
        integer month,day,hour,nperyear
        integer j,k,dpm(12)
        data dpm /31,29,31,30,31,30,31,31,30,31,30,31/
!
        if ( nperyear.eq.1 ) then
            j = 1
        elseif ( nperyear.eq.4 ) then
            j = 1 + month/3
            if ( j.gt.4 ) j = j - 4
        elseif ( nperyear.eq.12 ) then
            j = month                
        else
            j = day
            do k=1,month-1
                j = j + dpm(k)
            enddo
            if ( nperyear.ge.4*360 ) then
                j = 4*(j-1) + mod(hour,6) + 1
            elseif ( nperyear.lt.360 ) then
                j = 1 + (j-1)*nperyear/365
            endif
        endif
        getj = j
        end

        subroutine adjustnperyear(nperyear,lwrite)
!
!       some heuristics to get one of the offically-supported values
!       (1, 4, 12, 36, 73, 360,365,366,1464)
!
        implicit none
        integer nperyear
        logical lwrite
        integer nperyearold
!
        nperyearold = nperyear
        if ( nperyear.eq.0 ) then
            nperyear = 1        ! N-yearly data is yearly
        else if ( nperyear.le.4 ) then
            if ( .false. ) print *,'do nothing to nperyear'
        else if ( nperyear.le.14 ) then
            nperyear = 12       ! monthly data
        else if ( nperyear.le.40 ) then
            nperyear = 36
        else if ( nperyear.le.350 ) then
            nperyear = 73
        else if ( nperyear.lt.450 ) then
            nperyear = 366
        else if ( nperyear.lt.1550 ) then
            nperyear = 1464
        else if ( nperyear.lt.10000 ) then
            nperyear = 8784
        else
            write(0,*) 'Cannot handle more than 1464 values per year'
            write(0,*) '(6-hourly data)'
            write(0,*) 'but attempted to read ',nperyear
            call abort
        end if
        if ( lwrite .and. nperyear.ne.nperyearold ) then
            print *,'adjustnperyear: adjusted ',nperyearold,' to '
     +           ,nperyear
        end if
        end

        subroutine checkvalid(year,mon,val)
*
*       check for missing data and change to 3e33
*       note that this does NOT check for 999.9 - for too many
*       quantities his is a perfectly legal value (e.g., SLP, precip)
*       
        implicit none
        real*8 val
        real v
        integer year,mon
        logical lwrite
        parameter (lwrite = .false.)
*
        v = val
        if ( lwrite) print *,'checkvalid: val = ',v,year,mon
        if ( v.gt.1e28 .or. v.eq.-999.9 .or. v.eq.-999.8 .or. 
     +        v.eq.-999 .or. v.eq.-999.99 .or. v.eq.-9999 .or. 
     +        v.eq.-888.8 ) then
            if ( v.lt.1e33 .and. v.gt.-999 .and. v.ne.-888.8 )
     +            write(0,'(a,2i5,g20.4)')'# disregarding data point',
     +            year,mon,v
            val = 3e33
            if ( lwrite) print *,'checkvalid: val invalid '
        endif
        end
