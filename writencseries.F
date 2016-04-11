        subroutine writencseries(file,data,npermax,yrbeg,yrend,nperyear,
     +       title,comment,var,lvar,units)
*
*       writes a time series to a netCDF file
*       part of the KNMI Climate Explorer, GJvO, 2011
*
        implicit none
#include "netcdf.inc"
*       arguments
        integer npermax,yrbeg,yrend,nperyear
        real data(npermax,yrbeg:yrend)
        character file*(*),title*(*),comment*(*),var*(*),lvar*(*),
     +       units*(*)
*       local variables
        integer status,ncid,ntdimid,idim,i,j,l,ii(8),dy,mo,yr,yr0,yr1
     +       ,yr2,nt,ivar,ntvarid,year,month,n,nperday
        integer dpm(12)
        integer,allocatable :: itimeaxis(:)
        real array(1)
        real,allocatable :: linear(:)
        logical lwrite
        character string*10000,months(0:12,2)*3,clwrite*10
*       externals
        integer iargc,julday,leap
        external julday,leap
*       date
        data months /
     +        '???','JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG'
     +        ,'SEP','OCT','NOV','DEC','???','jan','feb','mar','apr'
     +        ,'may','jun','jul','aug','sep','oct','nov','dec'/
        data dpm /31,28,31,30,31,30,31,31,30,31,30,31/
        lwrite = .FALSE.
        call getenv('WRITENC_LWRITE',clwrite)
        if ( index(clwrite,'T') + index(clwrite,'t') .gt.0 ) then
            lwrite = .true.
        endif
        if ( lwrite ) then
            print *,'writencseries called with'
            print *,'file = ',trim(file)
            print *,'npermax,nperyear,yrbeg,yrend = ',npermax
     +           ,nperyear,yrbeg,yrend
            print *,'var,lvar,units = ',var,lvar,units
        endif
!
!       find beginning, end of data (rounded to the nearest whole year)
!
        yr1 = yrend
        yr2 = yrbeg
        do yr=yrbeg,yrend
            do mo=1,nperyear
                if ( data(mo,yr).lt.1e33 ) then
                    yr1 = min(yr1,yr)
                    yr2 = max(yr2,yr)
                end if
            end do
        end do
        nperday = max(1,nint(nperyear/365.24))
        if ( nperday.gt.1 ) then
            n = nperyear/nperday
        else
            n = nperyear
        end if
        if ( n.ne.366 ) then
            nt = nperyear*(yr2 - yr1 + 1)
        else
            i = julday(1,1,yr1)
            j = julday(12,31,yr2)
            nt = j-i+1
            if ( nperday.gt.1 ) nt = nperday*nt
        end if
        if ( lwrite ) then
            print *,'writencseries: first,last year with data ',yr1,yr2
            print *,'               nt = ',nt
        end if
        allocate(itimeaxis(nt))
        allocate(linear(nt))
*
*       open file, overwriting old file if it exists
*
        if ( lwrite ) print *,'calling nf_create(',trim(file),0,ncid,')'
        status = nf_create(file,0,ncid)
        if ( status.ne.nf_noerr ) call handle_err(status,file)
        if ( lwrite ) print *,'writenc: created ',file(1:len_trim(file))
     +        ,' with ncid = ',ncid
        status = nf_put_att_text(ncid,nf_global,'title',len_trim(title)
     +       ,title)
        if ( status.ne.nf_noerr ) call handle_err(status,'put att title'
     +        )
        status = nf_put_att_text(ncid,nf_global,'comment'
     +       ,len_trim(comment),comment)
        if ( status.ne.nf_noerr ) call handle_err(status
     +       ,'put att comment')
        status = nf_put_att_text(ncid,nf_global,'Conventions',6,'CF-1.0'
     +       )
        if ( status.ne.nf_noerr ) call handle_err(status
     +        ,'put att conventions')
        string = '\nwritten by writencseries (GJvO, KNMI) by '
        l = min(len_trim(string) + 2,len(string)-2)
        call getenv('USER',string(l:))
        l = min(len_trim(string) + 2,len(string)-2)
        call date_and_time(values=ii)
        ii(3) = ii(3)
        write(string(l:),'(i4,a,i2.2,a,i2.2)') ii(1),'-',ii(2),'-',ii(3)
        l = min(len_trim(string) + 2,len(string)-2)
        write(string(l:),'(i2,a,i2.2,a,i2.2)') ii(5),':',ii(6),':',ii(7)
        do i=0,iargc()
            l = min(len_trim(string) + 2,len(string)-2)
            call getarg(i,string(l:))
            if ( index(string(l:),'startstop').ne.0 ) then
                string(l:) = ' '
            endif
        enddo
        if ( lwrite ) then
            print *,'History: ',string(1:len_trim(string))
        endif
        status = nf_put_att_text(ncid,nf_global,'history'
     +       ,len_trim(string),string)
        if ( status.ne.nf_noerr ) call handle_err(status
     +        ,'put att history')
*
*       define dimension
*
        if ( lwrite ) print *,'defining time dimension with length ',nt
        if ( nt.gt.0 ) then
            status = nf_def_dim(ncid,'time',nt,ntdimid)
        else
            status = nf_def_dim(ncid,'time',nf_unlimited,ntdimid)
        endif
        if ( status.ne.nf_noerr ) call handle_err(status,'def time dim')
*
*       define variables: first the axis
*
        if ( lwrite ) print *,'defining time axis'
        status = nf_def_var(ncid,'time',nf_float,1,ntdimid,ntvarid)
        if ( status.ne.nf_noerr ) call handle_err(status,'def time var')
        if ( nperyear.eq.1 ) then
            string = 'years since '
        elseif ( nperyear.le.12 ) then
            string = 'months since '
        elseif ( nperyear.gt.12 .and. nperyear.le.366 ) then
            string = 'days since '
        elseif ( nperyear.gt.366 .and. nperyear.le.366*24 ) then
            string = 'hours since '
        else
            write(0,*) 'writencseries: cannot handle nperyear = ',
     +           nperyear
            write(0,*) '               in defining units string'
            call exit(-1)
        endif
        l = len_trim(string) + 2
        write(string(l:),'(i4,a,i2.2,a,i2.2)') yr1,'-',1,'-',1
        if ( nperyear.gt.366 ) string = trim(string)//' 00:00:00'
        if ( lwrite ) print *,'units = ',trim(string)
        status = nf_put_att_text(ncid,ntvarid,'units',len_trim(string)
     +       ,string)
        if ( status.ne.nf_noerr ) call handle_err(status
     +       ,'put time units')
        status = nf_put_att_text(ncid,ntvarid,'standard_name',4,'time')
        if ( status.ne.nf_noerr ) call handle_err(status
     +       ,'put time standard_name')
        status = nf_put_att_text(ncid,ntvarid,'long_name',4,'time')
        if ( status.ne.nf_noerr ) call handle_err(status
     +       ,'put time long_name')
        status = nf_put_att_text(ncid,ntvarid,'axis',1,'T')
        if ( status.ne.nf_noerr ) call handle_err(status
     +       ,'put time axis')
        if ( nperyear.lt.360 .or. n.eq.366 ) then
            status = nf_put_att_text(ncid,ntvarid,'calendar',9
     +           ,'gregorian')
        elseif ( n.eq.360 ) then
            status = nf_put_att_text(ncid,ntvarid,'calendar',7,
     +           '360_day')
        elseif ( n.eq.365 ) then
            status = nf_put_att_text(ncid,ntvarid,'calendar',7,
     +           '365_day')
        endif
*
*       next the variable itself
*
        if ( lwrite ) print *,'define variable'
        status = nf_def_var(ncid,var,nf_float,1,ntdimid,ivar)
        if ( status.ne.nf_noerr ) then ! concatenation does not work in f2c
            write(0,*) 'netCDF error: arguments were '
            write(0,*) 'ncid = ',ncid
            write(0,*) 'vars = ',var
            write(0,*) 'idim = ',ntdimid
            write(0,*) 'ivar = ',ivar
            string = 'def var '//var
            call handle_err(status,trim(string))
        endif
        status = nf_put_att_text(ncid,ivar,'long_name',len_trim(lvar)
     +       ,lvar)
        if ( status.ne.nf_noerr ) then
            string = 'def long_name '//lvar
            call handle_err(status,trim(string))
        endif
        if ( units.ne.' ' ) then
            status = nf_put_att_text(ncid,ivar,'units',len_trim(units)
     +           ,units)
            if ( status.ne.nf_noerr ) then
                string = 'def units '//lvar
                call handle_err(status,trim(string))
            endif
        endif
        array(1) = 3e33
        status = nf_put_att_real(ncid,ivar,'_FillValue',nf_float,1,array
     +       )
        if ( status.ne.nf_noerr ) then
            string = 'def _FillValue '//lvar
            call handle_err(status,trim(string))
        endif
*
*       end definition mode, put in data mode
*       
        if ( lwrite ) print *,'put in data mode'
        status = nf_enddef(ncid)
        if ( status.ne.nf_noerr ) call handle_err(status,'enddef')
*
*       write axes
*
        if ( n.eq.1 .or. n.eq.12 .or. 
     +       n.eq.366 .or. n.eq.365 .or. n.eq.360 ) then
            if ( nperday.le.1 ) then
                do i=1,nt
                    itimeaxis(i) = i-1
                end do
            else
                do i=1,nt
                    itimeaxis(i) = (24/nperday)*i
                end do
            end if
        else if ( nperyear.lt.12 ) then
            do i=1,nt
                itimeaxis(i) = (i-1)*12/nperyear
            end do
        else if ( nperyear.eq.36 ) then
            month = 1
            year = yr1
            itimeaxis(1) = 0
            do i=2,nt
                if ( mod(i,3).eq.2 .or. mod(i,3).eq.0 ) then
                    itimeaxis(i) = itimeaxis(i-1) + 10
                else
                    itimeaxis(i) = itimeaxis(i-1) + dpm(month) - 20
                    if ( leap(year).eq.2 ) then
                        itimeaxis(i) = itimeaxis(i) + 1
                    end if
                    month = month + 1
                    if ( month.gt.12 ) then
                        month = month - 12
                        year = year + 1
                    end if
                end if
            end do
        else if ( nperyear.lt.360 ) then
            do i=1,nt
                itimeaxis(i) = (i-1)*nint(365./nperyear)
            end do
        else
            write(0,*) 'writencseries: cannot handle nperyear = ',
     +           nperyear
            write(0,*) '               in defining time axis' 
            call exit(-1)
        end if
        if ( lwrite ) print *,'put time axis ',itimeaxis(1:nt)
        status = nf_put_var_int(ncid,ntvarid,itimeaxis)
        if ( status.ne.nf_noerr ) call handle_err(status,'put time')
*
*       write data
*
        i = 0
        do yr=yr1,yr2
            do mo=1,nperyear
                if ( leap(yr).eq.1 .and. 
     +               n.eq.366 .and. 1+(mo-1)/nperday.eq.60 ) cycle
                i = i + 1
                linear(i) = data(mo,yr)
            end do
        end do
        if ( i.ne.nt ) then
            write(0,*) 'writencseries: error: i != nt: ',i,nt
        end if
        status= nf_put_var_real(ncid,ivar,linear)
        if ( status.ne.nf_noerr ) call handle_err(status,'put var')
!
!       end game  - do not forget to cos ehe file or th last bits are
!       not flushed to disk...
!
        status = nf_close(ncid)
        end
