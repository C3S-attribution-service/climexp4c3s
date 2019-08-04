subroutine writenc(file,ncid,ntvarid,itimeaxis,ntmax,nx,xx,ny,yy &
    ,nz,zz,nt,nperyear,yrbegin,mobegin,undef,title,nvars &
    ,vars,ivars,lvars,units,nens1,nens2)

    implicit none
    integer :: ncid,ntvarid,ntmax,itimeaxis(ntmax),nx,ny,nz,nt,nperyear &
        ,yrbegin,mobegin,nvars,ivars(2,nvars),nens1,nens2
    real :: xx(nx),yy(ny),zz(nz),undef
    character file*(*),title*(*),vars(nvars)*(*),lvars(nvars)*(*), &
        units(nvars)*(*),history*100,cell_methods(nvars)*100
    character svars(1000)*1,ltime*10,lz(3)*20,metadata(2,100)*2000
    if ( nvars > 1000 ) then
        write(0,*) 'writenc: can only handle up to 1000 variables'
        call exit(-1)
    endif
    svars(1:nvars) = ' '
    ltime = 'time'
    history = ' '
    cell_methods = ' '
    lz = ' '
    metadata = ' '
    call enswritenc(file,ncid,ntvarid,itimeaxis,ntmax,nx,xx,ny &
        ,yy,nz,zz,lz,nt,nperyear,yrbegin,mobegin,ltime,undef,title &
        ,history,nvars,vars,ivars,lvars,svars,units,cell_methods &
        ,metadata,nens1,nens2)
end subroutine writenc

subroutine enswritenc(file,ncid,ntvarid,itimeaxis,ntmax,nx,xx,ny &
    ,yy,nz,zz,lz,nt,nperyear,yrbegin,mobegin,ltime,undef,title &
    ,history,nvars,vars,ivars,lvars,svars,units,cell_methods &
    ,metadata,nens1,nens2)

!       opens a netCDF file and puts all definitions in it.  The real
!       data is written later.  Note that this routine has been derived from
!       writectl and therefore assume GrADS-style data, i.e., a common grid,
!       only 2D/3D variables etc.

!       Arguments:
!       file:     (in)  name of output file
!       ncid:     (out) handle for subsequent netCDF calls
!       ntvarid:  (out) id of time variable, subsequent cals may need to extend
!                       the time axis
!       itimeaxis:(in/out) integer array of time tics
!       ntmax:    (in)  dimension of itemaxis
!       nx,xx:    (in)  X axis
!       ny,yy:    (in)  Y axis
!       nz,zz:    (in)  Z axis, not used for 2D field (ivars(1,i)==0)
!       lz(3)     (in)  1:units, 2:standard_name, 3:positive
!       nt:       (in)  time axis at this moment
!       nperyear: (in)  number of periods per year (12==monthly, 366==daily)
!       yrbegin:  (in)  year of origin of time axis
!       mobegin:  (in)  month of origin of time axis
!       ltime     (in)  long name of time variable
!       undef:    (in)  undefined
!       title:    (in)  title
!       nvars:    (in)  number of variables
!       vars:     (in)  short names of variables
!       ivars(1,):(in)  if 0, then 2D field, else 3D field
!       ivars(2,):(out) ID of variable for later use in writing data
!       lvars:    (in)  long names of variables
!       svars:    (in)  standard names of variables
!       units:    (in)  units of variables
!       cell_methods: (in)  cell_methods of variables
!       nens1,2:  (in)  ensemble members

    implicit none
    include 'netcdf.inc'
!   arguments
    integer :: ncid,ntvarid,ntmax,itimeaxis(ntmax),nx,ny,nz,nt,nperyear &
        ,yrbegin,mobegin,nvars,ivars(2,nvars),nens1,nens2
    real :: xx(nx),yy(ny),zz(nz),undef
    character file*(*),title*(*),history*(*),vars(nvars)*(*), &
        lvars(nvars)*(*),svars(nvars)*(*),units(nvars)*(*), &
        cell_methods(nvars)*(*),ltime*(*),lz(3)*(*),metadata(2,100)*(*)
!   local variables
    integer :: status,ntdimid,nxdimid,nydimid,nzdimid,nensdimid,ndims &
        ,idims(5),nxvarid,nyvarid,nzvarid,nensvarid,ivar,i,j,l, &
        dy,mo,chunks(5),n,nperday
    integer, allocatable :: iens(:)
    real :: array(1)
    integer :: dpm(12)
    logical :: lwrite,lhasz,lnetcdf4,linstitution
    character string*10000,months(0:12,2)*3,clwrite*10,FORM_field*100
!   date
    data months / &
        '???','JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG' &
        ,'SEP','OCT','NOV','DEC','???','jan','feb','mar','apr' &
        ,'may','jun','jul','aug','sep','oct','nov','dec'/
    data dpm /31,28,31,30,31,30,31,31,30,31,30,31/
    lwrite = .FALSE. 
    if ( nt > ntmax ) then
        write(0,*) 'writenc: error: increase ntmax to ',nt,' and recompile'
        call exit(-1)
    endif
    call getenv('WRITENC_LWRITE',clwrite)
    if ( index(clwrite,'T') + index(clwrite,'t') > 0 ) then
        lwrite = .TRUE. 
    endif
    
    if ( lwrite ) then
        print *,'writenc called with'
        print *,'file = ',trim(file)
        print *,'nx,xx = ',nx,(xx(i),i=1,nx)
        print *,'ny,yy = ',ny,(yy(i),i=1,ny)
        print *,'nz,zz = ',nz,(zz(i),i=1,nz)
        print *,'nt,nperyear,yrbegin,mobegin = ',nt,nperyear,yrbegin,mobegin
        print *,'undef = ',undef
        print *,'title = ',trim(title)
        print *,'nvars,vars,ivars,lvars,units = ',nvars
        do i=1,nvars
            print *,i,vars(i),ivars(1,i),lvars(i),units(i)
        enddo
        print *,'nens1,nens2 = ',nens1,nens2
    endif

!   open file, overwriting old file if it exists

    if ( index(file,'regionverification') /= 0 ) then
        if ( lwrite ) print *,'calling nf_create(',trim(file),'NF_64BIT_OFFSETL,ncid)'
        status = nf_create(file,NF_64BIT_OFFSET,ncid)
        lnetcdf4 = .FALSE. 
    else
        if ( lwrite ) print *,'calling nf_create(',trim(file),'NF_NETCDF4+NF_CLASSIC_MODEL,ncid)'
        status = nf_create(file,NF_NETCDF4+NF_CLASSIC_MODEL,ncid)
        lnetcdf4 = .TRUE. 
    end if
    if ( status /= nf_noerr ) call handle_err(status,file)
    if ( lwrite ) print *,'writenc: created ',trim(file),' with ncid = ',ncid
    if ( lwrite ) print *,'writenc: writing title ',trim(title)
    status = nf_put_att_text(ncid,nf_global,'title',len_trim(title),title)
    if ( status /= nf_noerr ) call handle_err(status,'put att title')
    status = nf_put_att_text(ncid,nf_global,'Conventions',6,'CF-1.0')
    if ( status /= nf_noerr ) call handle_err(status,'put att conventions')
    call getenv('FORM_field',FORM_field)
    if ( FORM_field /= ' ' ) then
        string = 'https://climexp.knmi.nl/select.cgi?field='//trim(FORM_field)
        if ( lwrite ) print *,'writenc: writing source_field ',trim(string)
        status = nf_put_att_text(ncid,nf_global,'source_field',len_trim(string),string)
        if ( status /= nf_noerr ) call handle_err(status,'put att derived_from')
    end if
    call write_metadata(ncid,metadata,lwrite)
    call extend_history(history)
    if ( lwrite ) print *,'writenc: writing history: ',trim(history)
    status = nf_put_att_text(ncid,nf_global,'history',len_trim(history),history)
    if ( status /= nf_noerr ) call handle_err(status,'put att history')

!   define dimensions

    if ( lwrite ) print *,'defining dimensions'
    if ( nt > 0 ) then
        status = nf_def_dim(ncid,'time',nt,ntdimid)
    else
        status = nf_def_dim(ncid,'time',nf_unlimited,ntdimid)
    endif
    if ( status /= nf_noerr ) call handle_err(status,'def time dim')
    status = nf_def_dim(ncid,'lon',nx,nxdimid)
    if ( status /= nf_noerr ) call handle_err(status,'def lon dim')
    status = nf_def_dim(ncid,'lat',ny,nydimid)
    if ( status /= nf_noerr ) call handle_err(status,'def lat dim')
    lhasz = .FALSE. 
    do i=1,nvars
        if ( ivars(1,i) > 0 ) then
            lhasz = .TRUE. 
        endif
    enddo
    if ( lhasz ) then
        status = nf_def_dim(ncid,'level',nz,nzdimid)
        if ( status /= nf_noerr ) call handle_err(status,'def level dim')
    endif
    if ( nens2 > 0 ) then
        status = nf_def_dim(ncid,'ensemble',nens2-nens1+1,nensdimid)
        if ( status /= nf_noerr ) call handle_err(status,'def level dim')
    endif

!   define variables: first the axes

    if ( lwrite ) print *,'defining axes'
    idims(1) = ntdimid
    status = nf_def_var(ncid,'time',nf_float,1,idims,ntvarid)
    if ( status /= nf_noerr ) call handle_err(status,'def time var')
    if ( nperyear == 1 ) then
        string = 'years since '
    elseif ( nperyear <= 12 ) then
        string = 'months since '
    elseif ( nperyear > 12 .and. nperyear <= 366 ) then
        string = 'days since '
    elseif ( nperyear > 366 ) then
        string = 'hours since '
    else
        write(0,*) 'writenc: cannot handle nperyear = ',nperyear
        write(0,*) '         in defining units string'
        call exit(-1)
    endif
    call getdymo(dy,mo,mobegin,nperyear)
    l = len_trim(string) + 2
    write(string(l:),'(i4,a,i2.2,a,i2.2)') yrbegin,'-',mo,'-',dy
    status = nf_put_att_text(ncid,ntvarid,'units',len_trim(string),string)
    if ( status /= nf_noerr ) call handle_err(status,'put time units')
    status = nf_put_att_text(ncid,ntvarid,'standard_name',4,'time')
    if ( status /= nf_noerr ) call handle_err(status,'put time standard_name')
    if ( ltime == ' ' ) ltime = 'time'
    status = nf_put_att_text(ncid,ntvarid,'long_name',len_trim(ltime),ltime)
    if ( status /= nf_noerr ) call handle_err(status,'put time long_name')
    status = nf_put_att_text(ncid,ntvarid,'axis',1,'T')
    if ( status /= nf_noerr ) call handle_err(status,'put time axis')
    if ( nperyear == 1 .or. nperyear == 12 .or. nperyear == 366 ) then
        status = nf_put_att_text(ncid,ntvarid,'calendar',9,'gregorian')
    elseif ( nperyear == 360 ) then
        status = nf_put_att_text(ncid,ntvarid,'calendar',7,'360_day')
    elseif ( nperyear == 365 ) then
        status = nf_put_att_text(ncid,ntvarid,'calendar',7,'365_day')
    endif

    idims(1) = nxdimid
    status = nf_def_var(ncid,'lon',nf_float,1,idims,nxvarid)
    if ( status /= nf_noerr ) call handle_err(status,'def lon var')
    status = nf_put_att_text(ncid,nxvarid,'units',12,'degrees_east')
    if ( status /= nf_noerr ) call handle_err(status,'put lat units')
    status = nf_put_att_text(ncid,nxvarid,'long_name',9,'Longitude')
    if ( status /= nf_noerr ) call handle_err(status,'put lat long_name')
    status = nf_put_att_text(ncid,nxvarid,'standard_name',9,'longitude')
    if ( status /= nf_noerr ) call handle_err(status,'put lat standard_name')
    status = nf_put_att_text(ncid,nxvarid,'axis',1,'X')
    if ( status /= nf_noerr ) call handle_err(status,'put lat axis')

    idims(1) = nydimid
    status = nf_def_var(ncid,'lat',nf_float,1,idims,nyvarid)
    if ( status /= nf_noerr ) call handle_err(status,'def lat var')
    status = nf_put_att_text(ncid,nyvarid,'units',13,'degrees_north')
    if ( status /= nf_noerr ) call handle_err(status,'put lat units')
    status = nf_put_att_text(ncid,nyvarid,'long_name',8,'Latitude')
    if ( status /= nf_noerr ) call handle_err(status,'put lat long_name')
    status = nf_put_att_text(ncid,nyvarid,'standard_name',8,'latitude')
    if ( status /= nf_noerr ) call handle_err(status,'put lat standard_name')
    status = nf_put_att_text(ncid,nyvarid,'axis',1,'Y')
    if ( status /= nf_noerr ) call handle_err(status,'put lat axis')

    if ( lhasz ) then
        idims(1) = nzdimid
        status = nf_def_var(ncid,'level',nf_float,1,idims,nzvarid)
        if ( status /= nf_noerr ) call handle_err(status,'def level var')
        if ( lz(1) /= ' ' ) then
            status = nf_put_att_text(ncid,nzvarid,'units',len_trim(lz(1)),lz(1))
            if ( status /= nf_noerr ) call handle_err(status,'put level units')
        endif
        if ( lz(2) /= ' ' ) then
            status = nf_put_att_text(ncid,nzvarid,'standard_name',len_trim(lz(2)),lz(2))
            if ( status /= nf_noerr ) call handle_err(status,'put level standard_name')
        endif
        if ( lz(3) /= ' ' ) then
            status = nf_put_att_text(ncid,nzvarid,'positive',len_trim(lz(3)),lz(3))
            if ( status /= nf_noerr ) call handle_err(status,'put level positive')
        endif
        status = nf_put_att_text(ncid,nzvarid,'axis',1,'Z')
    endif
    if ( nens2 > 0 ) then
        idims(1) = nensdimid
        status = nf_def_var(ncid,'ensemble',nf_int,1,idims,nensvarid)
        if ( status /= nf_noerr ) call handle_err(status,'def ensemble var')
        status = nf_put_att_text(ncid,nensvarid,'long_name',22,'ensemble member number')
        if ( status /= nf_noerr ) call handle_err(status,'put ensemble long_name')
    endif

!   next the variables themselves

    if ( lwrite ) print *,'define variables'
    do ivar=1,nvars
        if ( lhasz .and. ivars(1,ivar) > 0 ) then
            ndims = 4
            idims(1) = nxdimid
            idims(2) = nydimid
            idims(3) = nzdimid
            idims(4) = ntdimid
        else
            ndims = 3
            idims(1) = nxdimid
            idims(2) = nydimid
            idims(3) = ntdimid
        endif
        if ( nens2 > 0 ) then
            ndims = ndims + 1
            idims(ndims) = nensdimid
        endif
        status = nf_def_var(ncid,vars(ivar),nf_float,ndims,idims,ivars(2,ivar))
        if ( status /= nf_noerr ) then ! concatenation does not work in f2c
            write(0,*) 'netCDF error: arguments were '
            write(0,*) 'ncid = ',ncid
            write(0,*) 'vars = ',vars(ivar)
            write(0,*) 'ndims= ',ndims
            write(0,*) 'idims= ',idims
            write(0,*) 'ivars= ',ivars(2,ivar)
            string = 'def var'
            string(9:) = vars(ivar)
            call handle_err(status,trim(string))
        endif
        if ( lwrite ) print *,ivars(2,ivar),(trim(vars(ivar)))
        status = nf_put_att_text(ncid,ivars(2,ivar),'long_name',len_trim(lvars(ivar)),lvars(ivar))
        if ( status /= nf_noerr ) then
            string = 'def long_name '
            string(15:) = lvars(ivar)
            call handle_err(status,trim(string))
        endif
        if ( svars(ivar) /= ' ' ) then
            status = nf_put_att_text(ncid,ivars(2,ivar) &
                ,'standard_name',len_trim(svars(ivar)),svars(ivar))
            if ( status /= nf_noerr ) then
                string = 'def standard_name '
                string(19:) = svars(ivar)
                call handle_err(status,trim(string))
            endif
        endif
        if ( units(ivar) /= ' ' ) then
            status = nf_put_att_text(ncid,ivars(2,ivar),'units',len_trim(units(ivar)),units(ivar))
            if ( status /= nf_noerr ) then
                string = 'def units'
                string(11:) = lvars(ivar)
                call handle_err(status,trim(string))
            endif
        endif
        if ( cell_methods(ivar) /= ' ' ) then
            status = nf_put_att_text(ncid,ivars(2,ivar),'cell_methods', &
                len_trim(cell_methods(ivar)),cell_methods(ivar))
            if ( status /= nf_noerr ) then
                string = 'def cell_methods'
                string(18:) = lvars(ivar)
                call handle_err(status,trim(string))
            endif
        endif
        array(1) = undef
        status = nf_put_att_real(ncid,ivars(2,ivar),'_FillValue',nf_float,1,array)
        if ( status /= nf_noerr ) then
            string = 'def _FillValue '
            string(16:) = lvars(ivar)
            call handle_err(status,trim(string))
        endif
        if ( lnetcdf4 ) then
            chunks(1) = max(1,nx)
            chunks(2) = max(1,ny)
            chunks(3) = max(1,nz) ! can be 0
            chunks(4) = 1
            chunks(5) = 1
            if ( lwrite ) print *,'chunks = ',chunks
            status = nf_def_var_chunking(ncid,ivars(2,ivar),nf_chunked,chunks)
            if ( status /= nf_noerr ) then
                string = 'def chunking '
                string(16:) = lvars(ivar)
                call handle_err(status,trim(string))
            endif
            status = nf_def_var_deflate(ncid,ivars(2,ivar),0,1,3)
            if ( status /= nf_noerr ) then
                string = 'def deflate '
                string(16:) = lvars(ivar)
                call handle_err(status,trim(string))
            endif
        end if
    enddo

!   end definition mode, put in data mode

    if ( lwrite ) print *,'put in data mode'
    status = nf_enddef(ncid)
    if ( status /= nf_noerr ) call handle_err(status,'enddef')

!   write axes

    if ( nt > 0 ) then
        j = nt
    else
        j = ntmax
    endif
    nperday = max(1,nint(nperyear/365.))
    if ( nperday > 1 ) then
        n = nperyear/nperday
    else
        n = nperyear
    end if
    call maketimeaxis(nperyear,j,yrbegin,itimeaxis,lwrite)
    if ( lwrite ) print *,'put time axis ',itimeaxis(1:j)
    status = nf_put_var_int(ncid,ntvarid,itimeaxis)
    if ( status /= nf_noerr ) call handle_err(status,'put time')

    if ( lwrite ) print *,'put x axis ',xx(1:nx)
    status = nf_put_var_real(ncid,nxvarid,xx)
    if ( status /= nf_noerr ) call handle_err(status,'put lon')

    if ( lwrite ) print *,'put y axis ',yy(1:ny)
    status = nf_put_var_real(ncid,nyvarid,yy)
    if ( status /= nf_noerr ) call handle_err(status,'put lat')

    if ( lhasz ) then
        if ( lwrite ) print *,'put z axis ',zz(1:nz)
        status = nf_put_var_real(ncid,nzvarid,zz)
        if ( status /= nf_noerr ) call handle_err(status,'put level')
    endif
    if ( nens2 > 0 ) then
        allocate(iens(nens2+1))
        do j=1,nens2-nens1+1
            iens(j) = j+nens1-1
        enddo
        status = nf_put_var_int(ncid,nensvarid,iens)
        if ( status /= nf_noerr ) call handle_err(status,'put time')
        deallocate(iens)
    endif

!   that's it for this routine; data will be written later

end subroutine enswritenc
            
subroutine writencslice(ncid,ntvarid,itimeaxis,ntmax,ivars,data &
    ,nxf,nyf,nzf,nx,ny,nz,it,iens)

!       write a timeslice of data in the netCDF file

!       Arguments:
!       ncid:     (in) handle of open netCDF file in data mode
!       ntvarid:  (in) if != 0: ID of unlimited time axis to be updated
!       itimeaxis:(in) values of time axis to be updated
!       ntmax:    (in) dimension if itimeaxis
!       ivars(1): (in) if 0: 2D array
!       ivars(2): (in) ID of variable
!       data:     (in) data to be writte to netCDF file
!       n[xyz]f:  (in) dimensions of data array
!       n[xyz]:   (in) length of axes
!       it:       (in) number of time slice to be written
!       iens:     (in) number of ensemble member to be written (starting
!                      from nens1+1, not 0)

    implicit none
    include 'netcdf.inc'
!   arguments
    integer :: ncid,ntvarid,ntmax,itimeaxis(ntmax),ivars(2),nxf,nyf,nzf &
        ,nx,ny,nz,it,iens
    real :: data(nxf,nyf,nzf)
!   local variables
    integer :: status,icount(5),istart(5),istride(5),imap(5),k,indx(1),i,j
    logical :: lwrite
    character clwrite*10
    lwrite = .FALSE. 
    call getenv('WRITENC_LWRITE',clwrite)
    if ( index(clwrite,'T') + index(clwrite,'t') > 0 ) then
        lwrite = .TRUE. 
    endif
    if ( lwrite ) then
        print *,'writencslice '
        print *,'ncid,ivars ',ncid,ivars
        print *,'n[xyz]  =  ',nx,ny,nz
        print *,'n[xyz]f =  ',nxf,nyf,nzf
        print *,'it,iens =  ',it,iens
    endif
    if ( nx == nxf .and. ny == nyf .and. nz == nzf .and. it == 1 &
        .and. iens == 1 .and. ntvarid == 0 ) then
        if ( lwrite ) then
            print *,'calling nf_put_var_real with ',ncid,ivars(2)
            if ( .FALSE. ) then
                do j=1,ny
                    print '(i3,1000f6.1)',j,(data(i,j,1),i=1,nx)
                enddo
            endif
        endif
        status = nf_put_var_real(ncid,ivars(2),data)
        if ( status /= nf_noerr ) call handle_err(status,'put var')
    else
    
    !   call nf_put_vara_real
    
        icount(1) = nx
        istart(1) = 1
        istride(1)= 1
        imap(1)   = 1
        icount(2) = ny
        istart(2) = 1
        istride(2)= 1
        imap(2)   = nxf
        if ( ivars(1) == 0 ) then
            k = 3
            imap(k)   = nxf*nyf
        else
            icount(3) = nz
            istart(3) = 1
            istride(3)= 1
            imap(3)   = nxf*nyf
            k = 4
            imap(k)   = nxf*nyf*nzf
        endif
        icount(k) = 1
        istart(k) = it
        istride(k)= 1
        k = k + 1
        icount(k) = 1
        istart(k) = iens
        istride(k) = 1
        imap(k) = imap(k-1)
        if ( lwrite ) then
            print *,'writing to netCDF'
            print *,'count = ',icount(1:k)
            print *,'start = ',istart(1:k)
            print *,'stride= ',istride(1:k)
            print *,'map   = ',imap(1:k)
            print *,'data(',1+nx/2,1+ny/2,1+nz/2,') = ',data(1+nx/2,1+ny/2,1+nz/2)
        endif
        status = nf_put_varm_real(ncid,ivars(2),istart,icount,istride,imap,data)
        if ( status /= nf_noerr ) call handle_err(status,'put var')
    endif

!   and update the time axis

    if ( ntvarid /= 0 .and. ntmax > 0 ) then
        indx(1) = it
        status = nf_put_var1_int(ncid,ntvarid,indx,itimeaxis(it))
        if ( status /= nf_noerr ) call handle_err(status,'update T')
    endif

end subroutine writencslice

subroutine subtractleapyears(nt,firstyr,firstmo,nperyear,nnew)

!   subtract the undef leap days from nt so that the length of the real data is left
!   only needed when nperyear == 366 or a multiple thereof, which indicates a Gregorian calendar

    implicit none
    integer :: nt,firstyr,firstmo,nperyear,nnew
    integer :: yr,mo,it
    integer,external :: leap

    if ( mod(nperyear,366 ) /= 0 ) then
        nnew = nt
        return
    end if

    yr=firstyr
    mo=firstmo
    nnew = 0
    do it=1,nt
        if ( nperyear == 366 .and. mo == 60 .and. leap(yr) == 1 ) then
            mo = mo + 1
            cycle
        end if
        nnew = nnew + 1
        mo = mo + 1
        if ( mo > nperyear ) then
            mo = mo - nperyear
            yr = yr + 1
        end if
    end do

end subroutine subtractleapyears

subroutine write_metadata(ncid,metadata,lwrite)
!
!   write metadata to netcdf
!
    implicit none
    include 'netcdf.inc'
    integer :: ncid
    character :: metadata(2,100)*(*)
    logical :: lwrite
    integer :: i,j,status,ii,iscripturl
    real :: array(1)
    logical :: linstitution,illegal
    character :: scripturl*2000,string*20
    logical,external :: isnumchar,isalphachar
!
    linstitution = .false.
    iscripturl = 0
    do i=1,100
        if ( len_trim(metadata(1,i)) > 0 ) then
            illegal = .false.
            do j=1,len_trim(metadata(1,i))
                if ( .not. ( isnumchar(metadata(1,i)(j:j)) .or. isalphachar(metadata(1,i)(j:j)) .or. &
                        metadata(1,i)(j:j) == '_' .or. metadata(1,i)(j:j) == '-' ) ) then
                    illegal = .true.
                    if ( metadata(1,i)(j:j) == char(0) ) then
                        metadata(1,i)(j:j) = ' '
                    else
                        metadata(1,i)(j:j) = '_'
                    end if
                end if
            end do
            if ( illegal ) then
                write(0,*) 'writenc: warning: found illegal character  in metadata(1, ',i,'), char ', &
                    j,': ',trim(metadata(1,i))
            end if
            if ( metadata(1,i) == 'conventions' .or. metadata(1,i) == 'Conventions' ) cycle
            if ( metadata(1,i) == 'institution' .or. metadata(1,i) == 'Institution' ) then
                if ( metadata(2,i)(1:10) /= 'KNMI Clima' ) then ! avoid recursion...
                    metadata(2,i) = 'KNMI Climate Explorer and '//metadata(2,i)
                end if
                linstitution = .true.
            end if
            ! a few should be written as floats...
            if ( metadata(1,i) == 'geospatial_lon_min' .or. metadata(1,i) == 'geospatial_lon_max' .or. &
                 metadata(1,i) == 'geospatial_lat_min' .or. metadata(1,i) == 'geospatial_lat_max' .or. &
                 metadata(1,i) == 'geospatial_lon_resolution' .or. metadata(1,i) == 'geospatial_lat_resolution' ) then
                read(metadata(2,i),*,end=101,err=101) array(1)
                status = nf_put_att_real(ncid,nf_global,trim(metadata(1,i)),nf_float,1,array)
            101 continue
            else
                if ( lwrite ) print *,'writenc: writing ',trim(metadata(1,i)),': ',trim(metadata(2,i))
                status = nf_put_att_text(ncid,nf_global,trim(metadata(1,i)), &
                    len_trim(metadata(2,i)),trim(metadata(2,i)))
            end if
            if ( status /= nf_noerr ) call handle_err(status,'put att '//trim(metadata(1,i)))
            if ( metadata(1,i)(1:9) == 'scripturl' .and. &
                 isnumchar(metadata(1,i)(10:10)) .and. isnumchar(metadata(1,i)(11:11)) ) then
                read(metadata(1,i)(10:11),*) ii
                iscripturl = max(iscripturl,ii)
            end if
        end if
    end do
    if ( .not.linstitution ) then
        status = nf_put_att_text(ncid,nf_global,'institution',21,'KNMI Climate Explorer')
        if ( status /= nf_noerr ) call handle_err(status,'put att institution')
    end if

!   add scripturl of this call

    call getenv('SCRIPTURL',scripturl)
    if ( scripturl /= ' ' ) then
        ! avoid duplicates
        do i=1,100
            if ( trim(metadata(2,i)) == trim(scripturl) ) exit
        end do
        if ( i > 100 ) then
            iscripturl = iscripturl + 1
            write(string,'(a,i2.2)') 'scripturl',iscripturl
            status = nf_put_att_text(ncid,nf_global,trim(string),len_trim(scripturl),trim(scripturl))
            if ( status /= nf_noerr ) call handle_err(status,'put att scripturl')
        end if
    end if
end subroutine

subroutine maketimeaxis(nperyear,nt,yr1,itimeaxis,lwrite)
!
!   make a time axis from our nperyear convention
!
    implicit none
    integer :: nperyear,nt,yr1,itimeaxis(nt)
    logical :: lwrite
    integer :: dpm(12),nperday,n,i,month,year
    integer,external :: leap
    data dpm /31,28,31,30,31,30,31,31,30,31,30,31/

    nperday = max(1,nint(nperyear/365.))
    if ( nperday > 1 ) then
        n = nperyear/nperday
    else
        n = nperyear
    end if
    if ( n == 1 .or. n == 12 .or. n == 366 .or. n == 365 .or. n == 360 ) then
        if ( nperday == 1 .or. nperday == 24 ) then
            do i=1,nt
                itimeaxis(i) = i-1
            end do
        else
            do i=1,nt
                itimeaxis(i) = (24/nperday)*(i-1)
            end do
        end if
    else if ( nperyear < 12 ) then
        do i=1,nt
            itimeaxis(i) = (i-1)*12/nperyear
        end do
    else if ( nperyear == 36 ) then
        month = 1
        year = yr1
        itimeaxis(1) = 0
        do i=2,nt
            if ( mod(i,3) == 2 .or. mod(i,3) == 0 ) then
                itimeaxis(i) = itimeaxis(i-1) + 10
            else
                itimeaxis(i) = itimeaxis(i-1) + dpm(month) - 20
                if ( leap(year) == 2 ) then
                    itimeaxis(i) = itimeaxis(i) + 1
                end if
                month = month + 1
                if ( month.gt.12 ) then
                    month = month - 12
                    year = year + 1
                end if
            end if
        end do
    else if ( nperyear < 360 ) then
        do i=1,nt
            itimeaxis(i) = (i-1)*nint(365./nperyear)
        end do
    else
        write(0,*) 'maketimeaxis: cannot handle nperyear = ',nperyear
        write(0,*) '              in defining time axis' 
        call exit(-1)
    end if
end subroutine maketimeaxis
