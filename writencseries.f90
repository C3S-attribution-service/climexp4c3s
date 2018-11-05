subroutine writencseries(file,data,npermax,yrbeg,yrend,nperyear, &
    title,description,comment,metadata,history,var,lvar,units)
!
!   writes a time series to a netCDF file
!   part of the KNMI Climate Explorer, GJvO, 2011
!
    implicit none
    include 'netcdf.inc'

!   arguments
    integer :: npermax,yrbeg,yrend,nperyear
    real :: data(npermax,yrbeg:yrend)
    character :: file*(*),title*(*),description*(*),comment*(*),metadata(2,100)*(*), &
        history*(*),var*(*),lvar*(*),units*(*)

!   local variables
    integer :: status,ncid,ntdimid,idim,i,j,l,ii(8),dy,mo,yr,yr0,yr1
    integer :: yr2,nt,ivar,ntvarid,year,month,n,nperday
    integer :: firstmo,firstdy,chunks(5)
    integer,allocatable :: itimeaxis(:)
    real :: array(1)
    real,allocatable :: linear(:)
    logical :: lwrite,lcompress,lcomment
    character :: string*10000,months(0:12,2)*3,clwrite*10

!   externals
    integer,external :: julday,leap

!   date
    data months &
 &        /'???','JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG'  &
 &        ,'SEP','OCT','NOV','DEC','???','jan','feb','mar','apr' &
 &        ,'may','jun','jul','aug','sep','oct','nov','dec'/
    lwrite = .FALSE.
    call getenv('WRITENC_LWRITE',clwrite)
    if ( index(clwrite,'T') + index(clwrite,'t') .gt.0 ) then
        lwrite = .true.
    end if
    if ( lwrite ) then
        print *,'writencseries called with'
        print *,'file = ',trim(file)
        print *,'title = ',trim(title)
        print *,'description = ',trim(description)
        print *,'npermax,nperyear,yrbeg,yrend = ',npermax,nperyear,yrbeg,yrend
        print *,'var,lvar,units = ',var,lvar,units
    end if
!
!   find beginning, end of data (rounded to the nearest whole year)
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
!
!   open file, overwriting old file if it exists
!
    if ( lwrite ) print *,'calling nf_create(',trim(file),NF_NETCDF4+NF_CLASSIC_MODEL,ncid,')'
    status = nf_create(file,NF_NETCDF4+NF_CLASSIC_MODEL,ncid)
    if ( status.ne.nf_noerr ) call handle_err(status,file)
    if ( lwrite ) print *,'writenc: created ',file(1:len_trim(file)),' with ncid = ',ncid
    status = nf_put_att_text(ncid,nf_global,'title',len_trim(title),title)
    if ( status.ne.nf_noerr ) call handle_err(status,'put att title')
    status = nf_put_att_text(ncid,nf_global,'description',len_trim(description),description)
    if ( status.ne.nf_noerr ) call handle_err(status,'put att description')
    lcomment = .false.
    do i=1,100
        if ( metadata(1,i) == ' ' ) exit
        if ( metadata(1,i) == 'comment' .or. metadata(1,i) == 'Comment' ) then
            metadata(2,i) = trim(comment)//'\\n'//trim(metadata(2,i))
            lcomment = .true.
        end if
    end do
    if ( i < 100 .and. .not.lcomment ) then
        metadata(1,i) = 'comment'
        metadata(2,i) = comment
    end if
    call write_metadata(ncid,metadata,lwrite)
    call extend_history(history)
    if ( lwrite ) print *,'History: ',string(1:len_trim(string))
    status = nf_put_att_text(ncid,nf_global,'history',len_trim(history),history)
    if ( status.ne.nf_noerr ) call handle_err(status,'put att history')
    status = nf_put_att_text(ncid,nf_global,'Conventions',6,'CF-1.0')
    if ( status.ne.nf_noerr ) call handle_err(status,'put att conventions')
!
!   define dimension
!
    if ( lwrite ) print *,'defining time dimension with length ',nt
    if ( nt.gt.0 ) then
        status = nf_def_dim(ncid,'time',nt,ntdimid)
    else
        status = nf_def_dim(ncid,'time',nf_unlimited,ntdimid)
    end if
    if ( status.ne.nf_noerr ) call handle_err(status,'def time dim')
!
!   define variables: first the axis
!
    if ( lwrite ) print *,'defining time axis'
    status = nf_def_var(ncid,'time',nf_float,1,ntdimid,ntvarid)
    if ( status.ne.nf_noerr ) call handle_err(status,'def time var')
    if ( nperyear.eq.1 ) then
        string = 'years since '
        firstmo = 7  ! half-way through the year
        firstdy = 1
    elseif ( nperyear.le.12 ) then
        string = 'months since '
        firstmo = 1
        firstdy = 15  ! half-way through the month
    elseif ( nperyear.gt.12 .and. nperyear.le.366 ) then
        string = 'days since '
        firstmo = 1
        firstdy = 1
    elseif ( nperyear.gt.366 .and. nperyear.le.366*24 ) then
        string = 'hours since '
        firstmo = 1
        firstdy = 1
    else
        write(0,*) 'writencseries: cannot handle nperyear = ',nperyear
        write(0,*) '               in defining units string'
        call exit(-1)
    end if
    l = len_trim(string) + 2
    write(string(l:),'(i4,a,i2.2,a,i2.2)') yr1,'-',firstmo,'-',firstdy
    if ( nperyear.gt.366 ) string = trim(string)//' 00:00:00'
    if ( lwrite ) print *,'units = ',trim(string)
    status = nf_put_att_text(ncid,ntvarid,'units',len_trim(string),string)
    if ( status.ne.nf_noerr ) call handle_err(status,'put time units')
    status = nf_put_att_text(ncid,ntvarid,'standard_name',4,'time')
    if ( status.ne.nf_noerr ) call handle_err(status,'put time standard_name')
    status = nf_put_att_text(ncid,ntvarid,'long_name',4,'time')
    if ( status.ne.nf_noerr ) call handle_err(status,'put time long_name')
    status = nf_put_att_text(ncid,ntvarid,'axis',1,'T')
    if ( status.ne.nf_noerr ) call handle_err(status,'put time axis')
    if ( nperyear.lt.360 .or. n.eq.366 ) then
        status = nf_put_att_text(ncid,ntvarid,'calendar',9,'gregorian')
    elseif ( n.eq.360 ) then
        status = nf_put_att_text(ncid,ntvarid,'calendar',7,'360_day')
    elseif ( n.eq.365 ) then
        status = nf_put_att_text(ncid,ntvarid,'calendar',7,'365_day')
    end if
!
!   next the variable itself
!
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
    end if
    status = nf_put_att_text(ncid,ivar,'long_name',len_trim(lvar),lvar)
    if ( status.ne.nf_noerr ) then
        string = 'def long_name '//lvar
        call handle_err(status,trim(string))
    end if
    if ( units.ne.' ' ) then
        status = nf_put_att_text(ncid,ivar,'units',len_trim(units),units)
        if ( status.ne.nf_noerr ) then
            string = 'def units '//lvar
            call handle_err(status,trim(string))
        end if
    end if
    array(1) = 3e33
    status = nf_put_att_real(ncid,ivar,'_FillValue',nf_float,1,array)
    if ( status.ne.nf_noerr ) then
        string = 'def _FillValue '//lvar
        call handle_err(status,trim(string))
    end if
    lcompress = .false. ! increases file size...
    if ( lcompress ) then
        chunks = 1
        chunks(4) = nperyear
        if ( lwrite ) print *,'chunks = ',chunks
        status = nf_def_var_chunking(ncid,ivar,nf_chunked,chunks)
        if ( status /= nf_noerr ) then
            string = 'def chunking '
            string(16:) = lvar
            call handle_err(status,trim(string))
        end if
        status = nf_def_var_deflate(ncid,ivar,0,1,1)
        if ( status /= nf_noerr ) then
            string = 'def deflate '
            string(16:) = lvar
            call handle_err(status,trim(string))
        end if
    end if
!
!   end definition mode, put in data mode
!       
    if ( lwrite ) print *,'put in data mode'
    status = nf_enddef(ncid)
    if ( status.ne.nf_noerr ) call handle_err(status,'enddef')
!
!   write axes
!
    call maketimeaxis(nperyear,nt,yr1,itimeaxis,lwrite)
    if ( lwrite ) print *,'put time axis ',itimeaxis(1:nt)
    status = nf_put_var_int(ncid,ntvarid,itimeaxis)
    if ( status.ne.nf_noerr ) call handle_err(status,'put time')
!
!   write data
!
    i = 0
    do yr=yr1,yr2
        do mo=1,nperyear
            if ( leap(yr).eq.1 .and. n.eq.366 .and. 1+(mo-1)/nperday.eq.60 ) cycle
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
!   end game  - do not forget to close the file or the last bits are
!   not flushed to disk...
!
    status = nf_close(ncid)
end subroutine
