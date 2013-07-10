        subroutine readdat(data,npermax,nperyear,yrbeg,yrend,file)
*
*       this routines accepts 3 formats:
*       - a netCDF file with one time series
*       - an ascii file with
*         5 lines of comment
*         n lines with 'yyyy val_jan val_feb ... val_dec [val_year]'
*       - an ascii file with
*         5 lines of comment
*         n lines with year.fractionofyear value
*       - an ascii file with
*         5 lines of comment
*         n lines with year month value
*       - an ascii file with
*         5 lines of comment
*         n lines with year month day value
*       the data is returned in data, which is assumed to be filled
*       with undef's before.
*
*       I try to convert the units to my standard ones:
*       temp [C], prcp [mm/period], slp [hPa], ...
*
        implicit none
#include "netcdf.inc"
        integer npermax,nperyear,yrbeg,yrend
        real data(npermax,yrbeg:yrend)
        character file*(*)
        integer i,j,k,m,year,month,day,unit,status,ncid,dpm(12),
     +       dpm0(12),imonth
        logical lagain
        real*8 x,y,x1,y1,val12(366)
        character line*10000,units*10,var*10
        integer getnumwords,llen
        external getnumwords,llen
        save dpm,dpm0        
        data dpm0 /31,29,31,30,31,30,31,31,30,31,30,31/
        do i=1,12
            dpm(i) = dpm0(i)
        enddo
*
***        print *,'readddat: opening file ',file(1:index(file,' ')-1)
*
*       also handle netCDF time series
*
        if ( file.eq.'-' ) then
            unit = 5
        else
            status = nf_open(file,nf_nowrite,ncid)
            if ( status.ne.nf_noerr ) then
                call rsunit(unit)
                open(unit,file=file,status='old',err=902)
            else
*               it is a netCDF file
                call readncseries(file,data,npermax,nperyear,yrbeg,yrend
     +                ,ncid,var,units,.false.)
                return
            endif
        endif
*       it is an ascii file.
*       skip header, but try to find the variable name and the units
        units = ' '
        var = ' '
        do i=1,5
            read(unit,'(a)',err=903) line
***         print *,'read line ',line
        enddo
*       data is in the format
*       yyyy data(1) ... nao(nperyear)
*       first figure out nperyear
        nperyear = 12
        read(unit,'(a)',err=904,end=1800) line
        call tolower(line)
        imonth =  
     +       index(line,'jan')+index(line,'feb')+index(line,'mar')+
     +       index(line,'apr')+index(line,'may')+index(line,'jun')+
     +       index(line,'jul')+index(line,'aug')+index(line,'sep')+
     +       index(line,'oct')+index(line,'nov')+index(line,'dec')
        if ( imonth.gt.0 ) then
*       format is NN-MON-YYYY, possibly with the dashes.  Assume daily
           nperyear = 366
           if ( nperyear.gt.npermax ) then
              write(0,*) 'readdat: error: increase npermax ',npermax
     +             ,nperyear
              write(*,*) 'readdat: error: increase npermax ',npermax
     +             ,nperyear
              call abort
           endif
 10        continue
           lagain = .false.
 101       continue
           month = (index
     +          ('jan feb mar apr may jun jul aug sep oct nov dec'
     +          ,line(imonth:imonth+2)) + 3)/4
           if ( month.le.0 .or. month.gt.12 ) then
              if ( lagain ) then
                 write(0,'(a,i3,a)')
     +                '# readdat: cannot recognize month in ',imonth
     +                ,line(imonth:imonth+2)
                 write(0,'(2a)') '# ',line(1:llen(line))
                 goto 20
              else
                 lagain = .true.
                 imonth =  
     +                index(line,'jan') + index(line,'feb') +
     +                index(line,'mar') + index(line,'apr') +
     +                index(line,'may') + index(line,'jun') +
     +                index(line,'jul') + index(line,'aug') +
     +                index(line,'sep') + index(line,'oct') +
     +                index(line,'nov') + index(line,'dec')
                 goto 101
              endif
           endif
           j = imonth-1
           if ( line(j:j).eq.'-' ) j = j-1
           if ( j-1.le.0 ) goto 11
           read(line(j-1:j),*,end=1800,err=11) day
           if ( day.lt.1 .or. day.gt.31 ) goto 11
           goto 12
 11        continue
           write(0,'(a,i8,a,i3,a)') '# readdat: strange day ',day
     +          ,' in ',j,line(j-1:j)
           write(0,'(2a)') '# ',line(1:llen(line))
           goto 20
 12        continue
           read(line(imonth+3:),*,err=13) year
           goto 14
 13        continue
           write(0,'(a,i8,a)') '# readdat: error reading year in'
           write(0,'(2a)') '# ',line(1:llen(line))
           goto 20
 14        continue
           if ( year.lt.0 ) year = -year
           do m=1,month-1
              day = day + dpm(m)
           enddo
           read(line(imonth+8:),*) data(day,year)
***           print *,day,year,data(day,year)
 20        continue
           read(unit,'(a)',err=904,end=1800) line
           call tolower(line)
           goto 10
        endif
        nperyear = getnumwords(line) - 1
*       often there is a yearly mean/sum at the end
        if ( nperyear.eq.13 ) nperyear = 12
        if ( nperyear.gt.npermax ) then
            write(0,*) 'readdat: error: increase npermax ',npermax
     +            ,nperyear
            write(*,*) 'readdat: error: increase npermax ',npermax
     +            ,nperyear
            call abort
        endif
        if ( nperyear.gt.3 ) then
            read(line,*,err=904,end=300) year,(val12(j),j=1,nperyear)
  100       continue
***         print *,'read year,val12 = ',year,(val12(j),j=1,nperyear)
            if ( year.lt.yrbeg .or. year.gt.yrend ) then
                write(0,'(a,i8)') '# disregarding year ',year
            else
                do j=1,nperyear
                    call checkvalid(year,j,val12(j))
                    data(j,year) = val12(j)
                enddo
            endif
            read(unit,*,err=904,end=300) year,(val12(j),j=1,nperyear)
            goto 100
  300       continue
        elseif ( nperyear.eq.1 ) then
  401       continue
            call tolower(line)
            read(line,*,err=904,end=500) x,y
            year = int(x+1e-4)
            call checkvalid(year,1,y)
            if ( y.gt.1e33 ) then
               read(unit,'(a)',err=904,end=1800) line
               goto 401
            endif
            i = 0
  402       continue
            i = i + 1
            read(unit,*,err=904,end=500) x1,y1
            year = int(x+1e-4)
            j = nint(nperyear*(x-year)+0.75001)
            call checkvalid(year,j,y1)
            if ( y1.gt.1e33 ) goto 402
            nperyear = nint(i/(x1-x))
***            write(0,*) 'plotdat: found nperyear = ',
***     +            nperyear,i,x1,x,i/(x1-x)
            if ( nperyear.gt.npermax ) then
                write(0,*) 'readdat: error: increase npermax ',npermax
     +                ,nperyear
                write(*,*) 'readdat: error: increase npermax ',npermax
     +                ,nperyear
                call abort
            endif
            year = int(x+1e-4)
            j = nint(nperyear*(x-year)+0.75001)
            call checkvalid(year,j,y)
            if ( year.lt.yrbeg .or. year.gt.yrend ) then
                write(0,'(a,i6)') '# disregarding year ',year
            else                
***                print *,'data(',year,j,') = ',y
                data(j,year) = y
            endif
  450       continue
            year = int(x1+1e-4)
            j = nint(nperyear*(x1-year)+0.75001)
            call checkvalid(year,j,y1)
            if ( year.lt.yrbeg .or. year.gt.yrend ) then
                write(0,'(a,i6)') '# disregarding year ',year
            else
***                print *,'data(',year,j,') = ',y1
                data(j,year) = y1
            endif
            read(unit,*,err=904,end=500) x1,y1
            goto 450
  500       continue
        elseif ( nperyear.eq.2 ) then
            nperyear = 12       ! hard-coded
            read(line,*,err=904,end=700) year,j,y
  600       continue
            if ( j.gt.nperyear ) then
                year = year + (j-1)/nperyear
                j = 1 + mod(j-1,nperyear)
            endif
            call checkvalid(year,j,y)
            if ( year.lt.yrbeg .or. year.gt.yrend ) then
                write(0,'(a,i6)') '# disregarding year ',year
            else                
***                print *,'!! data(',j,year,') = ',y
                data(j,year) = y
            endif
            read(unit,*,err=904,end=700) year,j,y
            goto 600
  700       continue
        elseif ( nperyear.eq.3 ) then
            nperyear = 366       ! hard-coded
            if ( nperyear.gt.npermax ) then
                print *,'readdat: recompile main program with npermax'//
     +                ' at least 366, now it is ',npermax
                call abort
            endif
            read(line,*,err=905,end=900) year,month,day,y
  800       continue
            if ( month.lt.1 .or. month.gt.12 ) then
                print *,'invalid month ',month
                goto 899
            endif
            if ( day.lt.1 .or. day.gt.dpm(month) ) then
               if ( day.eq.30 .and. month.eq.2 .and. 
     $              data(31,year).gt.1e33) then
                  nperyear = 360
                  do i=1,12
                     dpm(i) = 30
                  enddo
                  do i=32,60
                     data(i-1,year) = data(i,year)
                  enddo
               else
                  write(0,'(a,i6)') '# invalid day ',day
                  goto 899
               endif
            endif
            j = day
            do i=1,month-1
                j = j + dpm(i)
            enddo
            if ( year.lt.yrbeg .or. year.gt.yrend ) then
                if ( year.ne.-999 ) write(0,'(a,i6)')
     +               '# disregarding year ',year
            else                
***                print *,'data(',j,year,') = ',y,day,month
                data(j,year) = y
            endif
  899       continue
            read(unit,*,err=905,end=900) year,month,day,y
            goto 800
  900       continue
        endif
 1800   continue
        if ( unit.ne.5 ) close(unit)
        return
 902    print *,'readdat: error opening ',trim(file)
        call abort
 903    print *,'readdat: error reading header of ',trim(file)
        print *,line
        call abort
  904   print *,'readdat: error reading data of ',trim(file)
        print '(a)',line(1:llen(line))
        print *,year,(val12(j),j=1,nperyear)
        call abort
  905   print *,'readdat: error reading yr,mo,dy,val data of ',
     +       trim(file)
        print *,year,month,day,y
        call abort
        end

