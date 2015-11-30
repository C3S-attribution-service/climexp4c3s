        subroutine writenc(file,ncid,ntvarid,itimeaxis,ntmax,nx,xx,ny,yy
     +       ,nz,zz,nt,nperyear,yrbegin,mobegin,undef,title,nvars
     +       ,vars,ivars,lvars,units,nens1,nens2)
        implicit none
        integer ncid,ntvarid,ntmax,itimeaxis(ntmax),nx,ny,nz,nt,nperyear
     +       ,yrbegin,mobegin,nvars,ivars(2,nvars),nens1,nens2
        real xx(nx),yy(ny),zz(nz),undef
        character file*(*),title*(*),vars(nvars)*(*),lvars(nvars)*(*),
     +       units(nvars)*(*),history*100,cell_methods(nvars)*100
        character svars(1000)*1,ltime*10,lz(3)*20
        if ( nvars.gt.1000 ) then
            write(0,*) 'writenc: can only handle up to 1000 variables'
            call abort
        endif
        svars(1:nvars) = ' '
        ltime = 'time'
        history = ' '
        cell_methods = ' '
        lz = ' '
        call enswritenc(file,ncid,ntvarid,itimeaxis,ntmax,nx,xx,ny
     +       ,yy,nz,zz,lz,nt,nperyear,yrbegin,mobegin,ltime,undef,title
     +       ,history,nvars,vars,ivars,lvars,svars,units,cell_methods
     +       ,nens1,nens2)
        end

        subroutine enswritenc(file,ncid,ntvarid,itimeaxis,ntmax,nx,xx,ny
     +       ,yy,nz,zz,lz,nt,nperyear,yrbegin,mobegin,ltime,undef,title
     +       ,history,nvars,vars,ivars,lvars,svars,units,cell_methods
     +       ,nens1,nens2)
*
*       opens a netCDF file and puts all definitions in it.  The real
*       data is written later.  Note that this routine has been derived from 
*       writectl and therefore assume GrADS-style data, i.e., a common grid,
*       only 2D/3D variables etc.
*
*       Arguments:
*       file:     (in)  name of output file
*       ncid:     (out) handle for subsequent netCDF calls
*       ntvarid:  (out) id of time variable, subsequent cals may need to extend
*                       the time axis
*       itimeaxis:(in/out) integer array of time tics
*       ntmax:    (in)  dimension of itemaxis
*       nx,xx:    (in)  X axis
*       ny,yy:    (in)  Y axis
*       nz,zz:    (in)  Z axis, not used for 2D field (ivars(1,i)==0)
*       lz(3)     (in)  1:units, 2:standard_name, 3:positive
*       nt:       (in)  time axis at this moment
*       nperyear: (in)  number of periods per year (12==monthly, 366==daily)
*       yrbegin:  (in)  year of origin of time axis
*       mobegin:  (in)  month of origin of time axis
*       ltime     (in)  long name of time variable
*       undef:    (in)  undefined
*       title:    (in)  title
*       nvars:    (in)  number of variables
*       vars:     (in)  short names of variables
*       ivars(1,):(in)  if 0, then 2D field, else 3D field
*       ivars(2,):(out) ID of variable for later use in writing data
*       lvars:    (in)  long names of variables
*       svars:    (in)  standard names of variables
*       units:    (in)  units of variables
*       cell_methods: (in)  cell_methods of variables
*       nens1,2:  (in)  ensemble members
*
        implicit none
#include "netcdf.inc"
*       arguments
        integer ncid,ntvarid,ntmax,itimeaxis(ntmax),nx,ny,nz,nt,nperyear
     +       ,yrbegin,mobegin,nvars,ivars(2,nvars),nens1,nens2
        real xx(nx),yy(ny),zz(nz),undef
        character file*(*),title*(*),history*(*),vars(nvars)*(*),
     +       lvars(nvars)*(*),svars(nvars)*(*),units(nvars)*(*),
     +       cell_methods(nvars)*(*),ltime*(*),lz(3)*(*)
*       local variables
        integer status,ntdimid,nxdimid,nydimid,nzdimid,nensdimid,ndims
     +       ,idims(5),nxvarid,nyvarid,nzvarid,nensvarid,ivar,i,j,l,
     +       ii(8),dy,mo,chunks(5)
        integer, allocatable :: iens(:)
        real array(1)
        integer dpm(12)
        logical lwrite,lhasz
        character string*10000,months(0:12,2)*3,clwrite*10
*       externals
        integer llen,iargc
*       date
        data months /
     +        '???','JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG'
     +        ,'SEP','OCT','NOV','DEC','???','jan','feb','mar','apr'
     +        ,'may','jun','jul','aug','sep','oct','nov','dec'/
        data dpm /31,28,31,30,31,30,31,31,30,31,30,31/
        lwrite = .FALSE.
        if ( nt.gt.ntmax ) then
            write(0,*) 'writenc: error: increase ntmax to ',nt
     +            ,' and recompile'
            call abort
        endif
        call getenv('WRITENC_LWRITE',clwrite)
        if ( index(clwrite,'T') + index(clwrite,'t') .gt.0 ) then
            lwrite = .true.
        endif
        if ( lwrite ) then
            print *,'writenc called with'
            print *,'file = ',trim(file)
            print *,'nx,xx = ',nx,(xx(i),i=1,nx)
            print *,'ny,yy = ',ny,(yy(i),i=1,ny)
            print *,'nz,zz = ',nz,(zz(i),i=1,nz)
            print *,'nt,nperyear,yrbegin,mobegin = ',nt,nperyear,yrbegin
     +           ,mobegin
            print *,'undef = ',undef
            print *,'title = ',trim(title)
            print *,'nvars,vars,ivars,lvars,units = ',nvars
            do i=1,nvars
                print *,i,vars(i),ivars(1,i),lvars(i),units(i)
            enddo
            print *,'nens1,nens2 = ',nens1,nens2
        endif
*
*       open file, overwriting old file if it exists
*
        if ( lwrite ) print *,'calling nf_create(',trim(file),
     +       'NF_NETCDF4,ncid)'
        status = nf_create(file,NF_NETCDF4,ncid)
        if ( status.ne.nf_noerr ) call handle_err(status,file)
        if ( lwrite ) print *,'writenc: created ',file(1:llen(file))
     +        ,' with ncid = ',ncid
        status = nf_put_att_text(ncid,nf_global,'title',llen(title)
     +       ,title)
        if ( status.ne.nf_noerr ) call handle_err(status,'put att title'
     +        )
        status = nf_put_att_text(ncid,nf_global,'Conventions',6,'CF-1.0'
     +       )
        if ( status.ne.nf_noerr ) call handle_err(status
     +        ,'put att conventions')
        string = trim(history)//'\nwritten by writenc (GJvO, KNMI) by '
        l = min(llen(string) + 2,len(string)-2)
        call getenv('USER',string(l:))
        l = min(llen(string) + 2,len(string)-2)
        call date_and_time(values=ii)
        ii(3) = ii(3)
        write(string(l:),'(i4,a,i2.2,a,i2.2)') ii(1),'-',ii(2),'-',ii(3)
        l = min(llen(string) + 2,len(string)-2)
        write(string(l:),'(i2,a,i2.2,a,i2.2)') ii(5),':',ii(6),':',ii(7)
        do i=0,iargc()
            l = min(llen(string) + 2,len(string)-2)
            call getarg(i,string(l:))
            if ( index(string(l:),'startstop').ne.0 ) then
                string(l:) = ' '
            endif
        enddo
        if ( lwrite ) then
            print *,'History: ',string(1:llen(string))
        endif
        status = nf_put_att_text(ncid,nf_global,'history',llen(string)
     +       ,string)
        if ( status.ne.nf_noerr ) call handle_err(status
     +        ,'put att history')
*
*       define dimensions
*
        if ( lwrite ) print *,'defining dimensions'
        if ( nt.gt.0 ) then
            status = nf_def_dim(ncid,'time',nt,ntdimid)
        else
            status = nf_def_dim(ncid,'time',nf_unlimited,ntdimid)
        endif
        if ( status.ne.nf_noerr ) call handle_err(status,'def time dim')
        status = nf_def_dim(ncid,'lon',nx,nxdimid)
        if ( status.ne.nf_noerr ) call handle_err(status,'def lon dim')
        status = nf_def_dim(ncid,'lat',ny,nydimid)
        if ( status.ne.nf_noerr ) call handle_err(status,'def lat dim')
        lhasz = .false.
        do i=1,nvars
            if ( ivars(1,i).gt.0 ) then
                lhasz = .true.
            endif
        enddo
        if ( lhasz ) then
            status = nf_def_dim(ncid,'level',nz,nzdimid)
            if ( status.ne.nf_noerr ) call handle_err(status
     +           ,'def level dim')
        endif
        if ( nens2.gt.0 ) then
            status = nf_def_dim(ncid,'ensemble',nens2-nens1+1,nensdimid)
            if ( status.ne.nf_noerr ) call handle_err(status
     +           ,'def level dim')
        endif
*
*       define variables: first the axes
*
        if ( lwrite ) print *,'defining axes'
        idims(1) = ntdimid
        status = nf_def_var(ncid,'time',nf_float,1,idims,ntvarid)
        if ( status.ne.nf_noerr ) call handle_err(status,'def time var')
        if ( nperyear.eq.1 ) then
            string = 'years since '
        elseif ( nperyear.le.12 ) then
            string = 'months since '
        elseif ( nperyear.gt.12 .and. nperyear.le.366 ) then
            string = 'days since '
        else
            write(0,*) 'writenc: cannot handle nperyear = ',nperyear
            call abort
        endif
        call getdymo(dy,mo,mobegin,nperyear)
        l = llen(string) + 2
        write(string(l:),'(i4,a,i2.2,a,i2.2)') yrbegin,'-',mo,'-',dy
        status = nf_put_att_text(ncid,ntvarid,'units',llen(string)
     +       ,string)
        if ( status.ne.nf_noerr ) call handle_err(status
     +       ,'put time units')
        status = nf_put_att_text(ncid,ntvarid,'standard_name',4,'time')
        if ( status.ne.nf_noerr ) call handle_err(status
     +       ,'put time standard_name')
        if ( ltime.eq.' ' ) ltime = 'time'
        status = nf_put_att_text(ncid,ntvarid,'long_name',llen(ltime)
     +       ,ltime)
        if ( status.ne.nf_noerr ) call handle_err(status
     +       ,'put time long_name')
        status = nf_put_att_text(ncid,ntvarid,'axis',1,'T')
        if ( status.ne.nf_noerr ) call handle_err(status
     +       ,'put time axis')
        if ( nperyear.eq.1 .or. nperyear.eq.12 .or. nperyear.eq.366 )
     +       then
            status = nf_put_att_text(ncid,ntvarid,'calendar',9
     +           ,'gregorian')
        elseif ( nperyear.eq.360 ) then
            status = nf_put_att_text(ncid,ntvarid,'calendar',7,
     +           '360_day')
        elseif ( nperyear.eq.365 ) then
            status = nf_put_att_text(ncid,ntvarid,'calendar',7,
     +           '365_day')
        endif
*
        idims(1) = nxdimid
        status = nf_def_var(ncid,'lon',nf_float,1,idims,nxvarid)
        if ( status.ne.nf_noerr ) call handle_err(status,'def lon var')
        status = nf_put_att_text(ncid,nxvarid,'units',12,'degrees_east')
        if ( status.ne.nf_noerr ) call handle_err(status,'put lat units'
     +       )
        status = nf_put_att_text(ncid,nxvarid,'long_name',9,'Longitude')
        if ( status.ne.nf_noerr ) call handle_err(status
     +       ,'put lat long_name')
        status = nf_put_att_text(ncid,nxvarid,'standard_name',9
     +       ,'longitude')
        if ( status.ne.nf_noerr ) call handle_err(status
     +       ,'put lat standard_name')
        status = nf_put_att_text(ncid,nxvarid,'axis',1,'X')
        if ( status.ne.nf_noerr ) call handle_err(status
     +       ,'put lat axis')
*
        idims(1) = nydimid
        status = nf_def_var(ncid,'lat',nf_float,1,idims,nyvarid)
        if ( status.ne.nf_noerr ) call handle_err(status,'def lat var')
        status = nf_put_att_text(ncid,nyvarid,'units',13,'degrees_north'
     +       )
        if ( status.ne.nf_noerr ) call handle_err(status,'put lat units'
     +       )
        status = nf_put_att_text(ncid,nyvarid,'long_name',8,'Latitude')
        if ( status.ne.nf_noerr ) call handle_err(status
     +        ,'put lat long_name')
        status = nf_put_att_text(ncid,nyvarid,'standard_name',8
     +       ,'latitude')
        if ( status.ne.nf_noerr ) call handle_err(status
     +        ,'put lat standard_name')
        status = nf_put_att_text(ncid,nyvarid,'axis',1,'Y')
        if ( status.ne.nf_noerr ) call handle_err(status
     +        ,'put lat axis')
*        
        if ( lhasz ) then
            idims(1) = nzdimid
            status = nf_def_var(ncid,'level',nf_float,1,idims,nzvarid)
            if ( status.ne.nf_noerr ) call handle_err(status
     +           ,'def level var')
            if ( lz(1).ne.' ' ) then
                status = nf_put_att_text(ncid,nzvarid,'units'
     +               ,len_trim(lz(1)),lz(1))
                if ( status.ne.nf_noerr ) call handle_err(status
     +               ,'put level units')
            endif
            if ( lz(2).ne.' ' ) then
                status = nf_put_att_text(ncid,nzvarid,'standard_name',
     +               len_trim(lz(2)),lz(2))
                if ( status.ne.nf_noerr ) call handle_err(status
     +               ,'put level standard_name')
            endif
            if ( lz(3).ne.' ' ) then
                status = nf_put_att_text(ncid,nzvarid,'positive',
     +               len_trim(lz(3)),lz(3))
                if ( status.ne.nf_noerr ) call handle_err(status
     +               ,'put level positive')
            endif
            status = nf_put_att_text(ncid,nzvarid,'axis',1,'Z')
        endif
        if ( nens2.gt.0 ) then
            idims(1) = nensdimid
            status = nf_def_var(ncid,'ensemble',nf_int,1,idims,nensvarid
     +           )
            if ( status.ne.nf_noerr ) call handle_err(status
     +           ,'def ensemble var')
            status = nf_put_att_text(ncid,nensvarid,'long_name',22
     +           ,'ensemble member number')
            if ( status.ne.nf_noerr ) call handle_err(status
     +           ,'put ensemble long_name')
        endif
*
*       next the variables themselves
*
        if ( lwrite ) print *,'define variables'
        do ivar=1,nvars
            if ( lhasz .and. ivars(1,ivar).gt.0 ) then
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
            if ( nens2.gt.0 ) then
                ndims = ndims + 1
                idims(ndims) = nensdimid
            endif
            status = nf_def_var(ncid,vars(ivar),nf_float,ndims,idims
     +            ,ivars(2,ivar))
            if ( status.ne.nf_noerr ) then ! concatenation does not work in f2c
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
            status = nf_put_att_text(ncid,ivars(2,ivar),'long_name'
     +           ,llen(lvars(ivar)),lvars(ivar))
            if ( status.ne.nf_noerr ) then
                string = 'def long_name '
                string(15:) = lvars(ivar)
                call handle_err(status,string)
            endif
            if ( svars(ivar).ne.' ' ) then
                status = nf_put_att_text(ncid,ivars(2,ivar)
     +               ,'standard_name',llen(svars(ivar)),svars(ivar))
                if ( status.ne.nf_noerr ) then
                    string = 'def standard_name '
                    string(19:) = svars(ivar)
                    call handle_err(status,string)
                endif
            endif
            if ( units(ivar).ne.' ' ) then
                status = nf_put_att_text(ncid,ivars(2,ivar),'units'
     +               ,llen(units(ivar)),units(ivar))
                if ( status.ne.nf_noerr ) then
                    string = 'def units'
                    string(11:) = lvars(ivar)
                    call handle_err(status,string)
                endif
            endif
            if ( cell_methods(ivar).ne.' ' ) then
                status = nf_put_att_text(ncid,ivars(2,ivar)
     +               ,'cell_methods',llen(cell_methods(ivar))
     +               ,cell_methods(ivar))
                if ( status.ne.nf_noerr ) then
                    string = 'def cell_methods'
                    string(18:) = lvars(ivar)
                    call handle_err(status,string)
                endif
            endif
            array(1) = undef
            status = nf_put_att_real(ncid,ivars(2,ivar),'_FillValue'
     +            ,nf_float,1,array)
            if ( status.ne.nf_noerr ) then
                string = 'def _FillValue '
                string(16:) = lvars(ivar)
                call handle_err(status,string)
            endif
            chunks(1) = max(1,nx)
            chunks(2) = max(1,ny)
            chunks(3) = max(1,nz) ! can be 0
            chunks(4) = 1
            chunks(5) = 1
            if ( lwrite ) print *,'chunks = ',chunks
            status = nf_def_var_chunking(ncid,ivars(2,ivar),nf_chunked,
     +           chunks)
            if ( status.ne.nf_noerr ) then
                string = 'def chunking '
                string(16:) = lvars(ivar)
                call handle_err(status,string)
            endif
            status = nf_def_var_deflate(ncid,ivars(2,ivar),0,1,3)
            if ( status.ne.nf_noerr ) then
                string = 'def deflate '
                string(16:) = lvars(ivar)
                call handle_err(status,string)
            endif
        enddo
*
*       end definition mode, put in data mode
*       
        if ( lwrite ) print *,'put in data mode'
        status = nf_enddef(ncid)
        if ( status.ne.nf_noerr ) call handle_err(status,'enddef')
*       
*       write axes
*       
        if ( nt.gt.0 ) then
            j = nt
        else
            j = ntmax
        endif
        if ( nperyear.eq.1 .or. nperyear.eq.12 .or. nperyear.eq.366 .or.
     +       nperyear.eq.365 .or. nperyear.eq.360 ) then
            do i=1,j
                itimeaxis(i) = i-1
            enddo
        elseif ( nperyear.lt.12 ) then
            do i=1,j
                itimeaxis(i) = (i-1)*12/nperyear
            enddo
        elseif ( nperyear.lt.360 ) then
            do i=1,j
                itimeaxis(i) = (i-1)*nint(365./nperyear)
            enddo
        else
            write(0,*) 'writenc: cannot handle nperyear = ',nperyear
            call abort
        endif
        if ( lwrite ) print *,'put time axis ',itimeaxis(1:j)
        status = nf_put_var_int(ncid,ntvarid,itimeaxis)
        if ( status.ne.nf_noerr ) call handle_err(status,'put time')
*
        if ( lwrite ) print *,'put x axis ',xx(1:nx)
        status = nf_put_var_real(ncid,nxvarid,xx)
        if ( status.ne.nf_noerr ) call handle_err(status,'put lon')
*
        if ( lwrite ) print *,'put y axis ',yy(1:ny)
        status = nf_put_var_real(ncid,nyvarid,yy)
        if ( status.ne.nf_noerr ) call handle_err(status,'put lat')
*
        if ( lhasz ) then
            if ( lwrite ) print *,'put z axis ',zz(1:nz)
            status = nf_put_var_real(ncid,nzvarid,zz)
            if ( status.ne.nf_noerr ) call handle_err(status
     +           ,'put level')
        endif
        if ( nens2.gt.0 ) then
            allocate(iens(nens2+1))
            do j=1,nens2-nens1+1
                iens(j) = j+nens1-1
            enddo
            status = nf_put_var_int(ncid,nensvarid,iens)
            if ( status.ne.nf_noerr ) call handle_err(status,'put time')
            deallocate(iens)
        endif
*
*       that's it for this routine; data will be written later
*
        end
        
        subroutine writencslice(ncid,ntvarid,itimeaxis,ntmax,ivars,data
     +        ,nxf,nyf,nzf,nx,ny,nz,it,iens)
*
*       write a timeslice of data in the netCDF file
*
*       Arguments:
*       ncid:     (in) handle of open netCDF file in data mode
*       ntvarid:  (in) if != 0: ID of unlimited time axis to be updated
*       itimeaxis:(in) values of time axis to be updated
*       ntmax:    (in) dimension if itimeaxis
*       ivars(1): (in) if 0: 2D array
*       ivars(2): (in) ID of variable
*       data:     (in) data to be writte to netCDF file
*       n[xyz]f:  (in) dimensions of data array
*       n[xyz]:   (in) length of axes
*       it:       (in) number of time slice to be written
*       iens:     (in) number of ensemble member to be written (starting
*                      from nens1+1, not 0)
*
        implicit none
#include "netcdf.inc"
*       arguments
        integer ncid,ntvarid,ntmax,itimeaxis(ntmax),ivars(2),nxf,nyf,nzf
     +        ,nx,ny,nz,it,iens
        real data(nxf,nyf,nzf)
*       local variables
        integer status,icount(5),istart(5),istride(5),imap(5),k,indx(1)
     +       ,i,j
        logical lwrite
        character clwrite*10
        lwrite = .FALSE.
        call getenv('WRITENC_LWRITE',clwrite)
        if ( index(clwrite,'T') + index(clwrite,'t') .gt.0 ) then
            lwrite = .true.
        endif
        if ( lwrite ) then
            print *,'writencslice '
            print *,'ncid,ivars ',ncid,ivars
            print *,'n[xyz]  =  ',nx,ny,nz
            print *,'n[xyz]f =  ',nxf,nyf,nzf
            print *,'it,iens =  ',it,iens
        endif
        if ( nx.eq.nxf .and. ny.eq.nyf .anD. nz.eq.nzf .and. it.eq.1
     +       .and. iens.eq.1 .and. ntvarid.eq.0 ) then
            if ( lwrite ) then
                print *,'calling nf_put_var_real with ',ncid,ivars(2)
                if ( .false. ) then
                    do j=1,ny
                        print '(i3,1000f6.1)',j,(data(i,j,1),i=1,nx)
                    enddo
                endif
            endif
            status = nf_put_var_real(ncid,ivars(2),data)
            if ( status.ne.nf_noerr ) call handle_err(status,'put var')
        else
*
*           call nf_put_vara_real
*
            icount(1) = nx
            istart(1) = 1
            istride(1)= 1
            imap(1)   = 1
            icount(2) = ny
            istart(2) = 1
            istride(2)= 1
            imap(2)   = nxf
            if ( ivars(1).eq.0 ) then
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
            if ( lwrite ) then
                print *,'writing to netCDF'
                print *,'count = ',icount(1:k)
                print *,'start = ',istart(1:k)
                print *,'stride= ',istride(1:k)
                print *,'map   = ',imap(1:k)
                print *,'data(',1+nx/2,1+ny/2,1+nz/2,') = ',
     +               data(1+nx/2,1+ny/2,1+nz/2)
            endif
            status = nf_put_varm_real(ncid,ivars(2),istart,icount
     +           ,istride,imap,data)
            if ( status.ne.nf_noerr ) call handle_err(status,'put var')
        endif
*
*       and update the time axis
*
        if ( ntvarid.ne.0 .and. ntmax.gt.0 ) then
            indx(1) = it
            status = nf_put_var1_int(ncid,ntvarid,indx,itimeaxis(it))
            if ( status.ne.nf_noerr ) call handle_err(status,'update T')
        endif
*
        end
