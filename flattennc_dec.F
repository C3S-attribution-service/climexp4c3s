        program flattennc
!
!       flatten a netcdf file with ECMWF conventions to one
!       that the Climate Explorer can handle, with one time 
!       axis without holes.
!
        implicit none
        integer nvarmax
        parameter (nvarmax=20)
#include "params.h"
#include "netcdf.inc"
        integer i,j,k,l,m,n,dy,mo,yr,lead,dyref,moref,yrref,status,varid
     +       ,xtype,ndimvar,dimids(nf_max_var_dims),natts,ifile
     +       ,itype,offset,yr1,yr2,iarray(8),mtmax
        integer nx,ny,nz,nens1,nens2,nt,mt,ntmax
        integer ncid1,ndims1,nvars1,ngatts1,unlimdimid1
        integer ncid2,ntvarid,itimeaxis(ndata)
        integer ix,iy,iz,it,jt,ie,ntvars,firstmo,firstyr,nperyear
     +       ,iperyear,mperyear,jvars(6,nvarmax),ivars(2,nvarmax),lmax
     +       ,dpm(12),jul0,jul1,nfac,firstmo0,firstyr0
        real,allocatable :: data(:,:,:,:,:,:,:)
        real undef,xx(nxmax),yy(nymax),zz(nzmax),dt
        real*8 tl(ndata),tr(ndata),tt(ndata),t
        character infile*255,outfile*255,title*(nf_max_name),
     +       history*2000,name*(nf_max_name),
     +       leadunits*(nf_max_name),source*100,source0*100
        character vars(nvarmax)*40,lvars(nvarmax)*120,svars(nvarmax)*120
     +       ,units(nvarmax)*40,cell_methods*128,ltime*100,lz(3)*20
     +       ,months(12)*3,clwrite*10,lunits*10,metadata(2,100)*2000
        logical lwrite,foundreftime,foundleadtime
        integer iargc,julday
        data months 
     +       /'jan','feb','mar','apr','may','jun'
     +       ,'jul','aug','sep','oct','nov','dec'/
        data dpm /31,29,31,30,31,30,31,31,30,31,30,31/

        tt = 3e33
        lwrite = .false.
        call getenv('FLATTENNC_LWRITE',clwrite)
        if ( index(clwrite,'T') + index(clwrite,'t') .gt.0 ) then
            lwrite = .true.
            print *,'flattennc: debug output requested'
        endif

        if ( iargc().lt.2 ) then
            print *,'usage: flattennc infile1 infile2 ... outfile'
            stop
        endif
        call getarg(iargc(),outfile)
        do ifile=1,iargc()-1
            call getarg(ifile,infile)
!
!           open files
!
            status = nf_open(infile,nf_nowrite,ncid1)
            if ( status.ne.nf_noerr ) call handle_err(status,infile)
!
!           construct a sensible title
!
            if ( ifile.eq.1 ) then
                call gettitle(ncid1,title,lwrite)
                call gettextattopt(ncid1,nf_global,'history',history
     +               ,lwrite)
            end if
            call getnumbers(ncid1,ndims1,nvars1,ngatts1,unlimdimid1
     +           ,lwrite)
            do varid=1,nvars1
*               get name of variable
                status = nf_inq_var(ncid1,varid,name,xtype,ndimvar
     +               ,dimids,natts)
                if ( status.ne.nf_noerr ) call handle_err(status
     +               ,'nf_inq_var')
                if ( name.eq.'source' ) then
!                   get lengths
                    n = 1
                    do i=1,ndimvar
                        status = nf_inq_dim(ncid1,dimids(i),name,m)
                        if ( status.ne.nf_noerr ) call handle_err(status
     +                       ,'nf_inq_dim')
                        n = m*n
                    enddo
                    status = nf_get_var_text(ncid1,varid,source)
                    if ( status.ne.nf_noerr ) call handle_err(status
     +                   ,'nf_get_var_text')
                    do i=1,len(source)
                        if ( source(i:i).eq.char(0) ) source(i:i) = ' '
                    end do
                    if ( ifile.eq.1 ) then
                        title(len_trim(title)+2:) = source
                        source0 = source
                    else
                        if ( source.ne.source0 ) then
                            write(0,*) 'flattennc: error: data from '//
     +                           'different source: ',source,source0
                            call abort
                        end if
                    end if
                endif
            enddo
            if ( lwrite ) print *,'title = ',trim(title)
!       
!           read axes
!
            call getdims(ncid1,ndims1,ix,nx,nxmax,iy,ny,nymax,iz,nz
     +           ,nzmax,it,nt,ndata,ie,nens1,nens2,nensmax,lwrite)
            if ( nens1.ne.nens2 ) then
                write(0,*) 'flattennc: multiple ensemble members does'
     +               ,' not yet work'
                call abort
            end if
            ntvars = 0
            undef = 3e33
            foundreftime = .false.
            foundleadtime = .false.
            xx(1) = 0
            yy(1) = 0
            zz(1) = 0
            do varid=1,nvars1
*               get dimensions of variable
                status = nf_inq_var(ncid1,varid,name,xtype,ndimvar
     +               ,dimids,natts)
                if ( status.ne.nf_noerr ) call handle_err(status
     +               ,'nf_inq_var')
                if ( lwrite ) print *,'investigating variable ',varid
     +               ,trim(name),dimids(1:ndimvar)
                if ( index(name,'_bnd').ne.0 ) then
                    if ( lwrite ) print *
     +                   ,'flattennc: disregarding boundary ',trim(name)
                    cycle
                endif
*               what kind of variable do we have?
                if ( ndimvar.eq.1 .and. dimids(1).eq.ix ) then
                    call getdiminfo('x',ncid1,varid,xx,nx,lwrite)
                    call makelonreasonable(xx,nx)
                elseif ( ndimvar.eq.1 .and. dimids(1).eq.iy ) then
                    call getdiminfo('y',ncid1,varid,yy,ny,lwrite)
                elseif ( ndimvar.eq.1 .and. dimids(1).eq.iz ) then
                    call getzdiminfo('z',ncid1,varid,zz,nz,lz,lwrite)
                elseif ( ndimvar.eq.1 .and. dimids(1).eq.ie ) then
                    if ( lwrite ) print *,'ensemble axis'
                elseif ( name.eq.'reftime' .and. dimids(1).eq.it  ) then
                    foundreftime = .true.
                    call getreftime(ncid1,varid,tr,nt,firstmo,firstyr
     +                   ,nperyear,iperyear,lwrite)
                    if ( ifile.eq.1 ) then
                        firstmo0 = firstmo
                        firstyr0 = firstyr
                        if ( nperyear.eq.366 ) then
                            jul0 = julday(1,firstmo0,firstyr0)
                        else if ( nperyear.eq.12 ) then
                            jul0 = julday(firstmo0,1,firstyr0)
                        else
                            write(0,*) 'flattennc: cannot handle '//
     +                           'nperyear = ',nperyear,' yet'
                            call abort
                        end if
                        if ( lwrite ) then
                            call caldat(jul0,mo,dy,yr)
                            print '(a,i4,a,i2.2,a,i2.2)',
     +                           'Starting date: ',yr,'-',mo,'-',dy
                        end if
                    end if
                elseif ( name.eq.'leadtime' .and. dimids(1).eq.it ) then
                    foundleadtime = .true.
                    call getleadtime(ncid1,varid,tl,nt,leadunits,lwrite)
                else
                    n = 0
                    m = 0
                    do i=1,ndimvar
                        if ( it.ne.0 .and. dimids(i).eq.it ) then
                            n = n+1
                            if ( lwrite ) print *
     +                           ,'flattennc: time-varying variable '
     +                           ,varid
                        elseif ( ix.ne.0 .and. dimids(i).eq.ix ) then
                            m = m+1
                            if ( lwrite ) print *
     +                           ,'flattennc: x-dependent variable '
     +                           ,varid
                        elseif ( iy.ne.0 .and. dimids(i).eq.iy ) then
                            m = m+1
                            if ( lwrite ) print *
     +                           ,'flattennc: y-dependent variable '
     +                           ,varid
                        elseif ( iz.ne.0 .and. dimids(i).eq.iz ) then
                            m = m+1
                            if ( lwrite ) print *
     +                           ,'flattennc: z-dependent variable '
     +                           ,varid
                        endif
                    enddo
                    if ( n.gt.0 .and. m.gt.0 ) then
                        call addonevariable(ncid1,varid,name,ntvars
     +                       ,nvarmax,ndimvar,dimids,ix,iy,iz,it,ie,vars
     +                       ,jvars,lvars,svars,units,cell_methods,undef
     +                       ,lwrite)
                        if ( jvars(4,ntvars).eq.0 ) then
                            ivars(1,ntvars) = 0
                        else
                            ivars(1,ntvars) = nz
                        endif
                    endif
                endif
            enddo
!
!           flatten time axis
!
            if ( .not.foundreftime ) then
                write(0,*) 'flattennc: error: did not find reftime'
                write(*,*) 'flattennc: error: did not find reftime'
                call abort
            endif
            if ( .not.foundleadtime ) then
                write(0,*) 'flattennc: error: did not find leadtime'
                write(*,*) 'flattennc: error: did not find leadtime'
                call abort
            endif
!
            if ( iperyear.eq.366 .and. leadunits.eq.'days' .or.
     +           iperyear.eq.12  .and. leadunits.eq.'months' ) then
                nfac = 1
            elseif ( iperyear.eq.366 .and. leadunits.eq.'hours' ) then
                nfac = 24
            else
                write(0,*) 'flattennc: cannot handle iperyear,leadunits
     +               = ',iperyear,trim(leadunits),' yet'
                write(*,*) 'flattennc: cannot handle iperyear,leadunits
     +               = ',iperyear,trim(leadunits),' yet'
                call abort
            endif
!
!           define or extend T axis
!
            if ( ifile.eq.1 ) then
                do it=1,nt
                    tt(it) = tr(it) + tl(it)/nfac
                    if ( lwrite ) print *,'tt(',it,') = ',tt(it)
                end do
                dt = tt(2) - tt(1)
                ntmax = nt
                offset = 0
            else
                t =  tr(1) + tl(1)/nfac
                if ( lwrite ) then
                    print *,'searching for ',t
                end if
                do it=1,ntmax
                    if ( abs(tt(it)-t).lt.0.2*dt ) exit
                end do
                offset = it-1
                if ( lwrite ) print *,'offset = ',offset
                do it=1,nt
                    ntmax = max(ntmax,it+offset)
                    jt = it + offset
                    if ( tt(jt).ge.1e30 ) then
                        tt(jt) =  tr(it) + tl(it)/nfac
                        if ( lwrite ) print *,'tt(',jt,') = ',tt(jt)
                        if ( abs(tt(jt)-tt(jt-1)-dt).gt.0.2*dt ) then
                            write(0,*) 'flattennc_dec: error: unequal'//
     +                           'steps:',tt(jt-1),tt(jt),tt(jt)-tt(jt-1
     +                           ),dt
                            call abort
                        end if
                    else
                        t = tr(it) + tl(it)/nfac
                        if ( abs(tt(jt)-t).gt.0.2*dt ) then
                            write(0,*) 'flattennc: error: times do not '
     +                           //'agree, tt(',jt,') = ',tt(jt),' != '
     +                           ,t
                            call abort
                        end if
                    end if
                end do
            end if
!
!           figure out real time steps
!
            if ( ifile.eq.1 ) then
                mperyear = nint(iperyear/(tt(2)-tt(1)))
                if ( mperyear.eq.13 ) mperyear = 12 ! round-off errors in February
                if ( lwrite ) then
                    print *,'flattennc: deduced that mperyear = '
     +                   ,mperyear,iperyear,tt(2),tt(1)
                endif
!
!               allocate data array
!
                yr1 = firstyr0
                call date_and_time(values=iarray)
                yr2 = iarray(1)
                mtmax = nt + mperyear*(yr2-yr1+1)
                if ( lwrite ) print *,'yr1,yr2,nt,mtmax,ntvars = ',yr1
     +               ,yr2,nt,mtmax,ntvars
                if ( lwrite ) print *,'allocating data ',nx,ny,nz,mtmax
     +               ,nens1,nens2,4*nx*ny*nz*mtmax*(nens2-nens2+1)
     +               *ntvars
                allocate(data(nx,ny,nz,mtmax,nens1:nens2,4,ntvars))
                data = undef
            end if
!
!           read data
!
            do it=1,nt
                do ie=nens1,nens2
                    i = (it-1)/mperyear
                    if ( i.eq.0 ) then
                        itype = 1 ! first year
                    elseif ( i.lt.5 ) then
                        itype = 3 ! year 2-5
                    elseif ( i.lt.9 ) then
                        itype = 4 ! year 6-9
                    else
                        cycle ! throw away year 10
                    end if
                    do varid=1,ntvars
                        if ( data(1,1,1,it+offset,ie,itype,varid).ne.
     +                       undef ) then
                            write(0,*) 'flattennc: error: overwriting '
     +                           //' data',it+offset,ie,itype,varid,
     +                           data(1,1,1,it+offset,ie,itype,varid)
                            call abort
                        end if
                        if ( it+offset.gt.ntmax ) then
                            write(0,*) 'flattennc: error: array too
     +                           small: ',it+offset,ntmax
                            call abort
                        end if
                        if ( lwrite ) print *,'reading slice at t,e = ',
     +                       it+offset,ie,varid,itype
                        call ensreadncslice(ncid1,jvars(1,varid),it,
     +                       ie-nens1+1,
     +                       data(1,1,1,it+offset,ie,itype,varid),
     +                       nx,ny,nz,lwrite)
                    end do
                end do
            end do
        end do                  ! ifile
!
!       write header
!
!       convert starting dates to new mperyear
        if ( mperyear.ne.nperyear ) then
            if ( mperyear.ne.366 ) then
                firstmo0 = 1 + (firstmo0-1)
     +               *max(mperyear,12)/max(nperyear,12)
            elseif ( nperyear.le.12 ) then
                j = 1
                do i=1,firstmo0-1
                    j = j + dpm(i)
                enddo
                firstmo0 = j
            else
                write(0,*) 'flattennc: cannot transform starting date'
                write(*,*) 'flattennc: cannot transform starting date'
                call abort
            endif
            if ( lwrite ) print *,'transformed firstmo0  to ',firstmo0
        endif
!
!       a few more informative variables
!
        if ( mperyear.le.12 .or. leadunits.eq.'months' ) then
            if ( leadunits.eq.'months' ) then
                k = min(12,1+lmax)
            elseif ( leadunits.eq.'days' ) then
                k = min(12,nint(12*lmax/365.25))
            elseif ( leadunits.eq.'hours' ) then
                k = min(12,nint(12*lmax/(24*365.25)))
            endif
            lunits = 'months'
        else
            k = 1+lmax
            lunits = 'days'
        endif
        ltime = 'verification time (reftime+leadtime)'
        call getdymo(dy,mo,firstmo0,mperyear)
        write(title,'(aa,i2,a)') trim(title),', reftime ',dy
     +       ,months(mo)
        do itype=1,4
            if ( itype.eq.2 ) cycle ! no more yr2
            call getarg(iargc(),outfile)
            i = index(outfile,'.nc')
            if ( i.eq.0 ) then
                i = len_trim(outfile) + 1
            end if
            if ( itype.eq.1 ) then
                outfile(i:) = '_1.nc'
            else if ( itype.eq.2 ) then
                outfile(i:) = '_2.nc'
            else if ( itype.eq.3 ) then
                outfile(i:) = '_2-5.nc'
            else if ( itype.eq.4 ) then
                outfile(i:) = '_6-9.nc'
            else
                write(0,*) 'flattennc: error: itype should be 1-4, not '
     +               ,itype
                call abort
            end if
            call enswritenc(outfile,ncid2,ntvarid,itimeaxis,ndata,nx,xx
     +           ,ny,yy,nz,zz,lz,ntmax,mperyear,firstyr0,firstmo0,ltime
     +           ,undef,title,history,ntvars,vars,ivars,lvars,svars
     +           ,units,cell_methods,metadata,nens1,nens2)
            do it=1,ntmax
                do ie=nens1,nens2
                    do varid=1,ntvars
                        if ( lwrite ) then
                            print *,'write field ',it,ie,varid
                        end if
                        call writencslice(ncid2,0,itimeaxis,ndata,ivars
     +                       ,data(1,1,1,it,ie,itype,varid),nx,ny,nz,nx
     +                       ,ny,nz,it,ie-nens1)
                    end do
                end do
            end do
!           do not forget to close files!  otherwise the last bit is
!           lost
            status = nf_close(ncid1)
            status = nf_close(ncid2)
        end do                  ! itype
!
        end
