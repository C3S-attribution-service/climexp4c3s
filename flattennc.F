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
     +       ,xtype,ndimvar,dimids(nf_max_var_dims),natts
        integer nx,ny,nz,nens1,nens2,nt,mt
        integer ncid1,ndims1,nvars1,ngatts1,unlimdimid1
        integer ncid2,ntvarid,itimeaxis(ndata)
        integer ix,iy,iz,it,ie,ntvars,firstmo,firstyr,nperyear,iperyear,
     +       mperyear,jvars(6,nvarmax),ivars(2,nvarmax),lmax,dpm(12),
     +       jul0,nfac
        real,allocatable :: data(:,:,:)
        real undef,xx(nxmax),yy(nymax),zz(nzmax)
        real*8 tt(ndata),tl(ndata)
        character infile*255,outfile*255,title*(nf_max_name),
     +       history*2000,name*(nf_max_name),
     +       leadunits*(nf_max_name)
        character vars(nvarmax)*40,lvars(nvarmax)*120,svars(nvarmax)*120
     +       ,units(nvarmax)*40,cell_methods*128,ltime*100,lz(3)*20
     +       ,months(12)*3,clwrite*10,lunits*10,metadata(2,100)*2000
        logical lwrite,foundreftime,foundleadtime
        integer iargc,julday
        data months 
     +       /'jan','feb','mar','apr','may','jun'
     +       ,'jul','aug','sep','oct','nov','dec'/
        data dpm /31,29,31,30,31,30,31,31,30,31,30,31/

        lwrite = .false.
        call getenv('FLATTENNC_LWRITE',clwrite)
        if ( index(clwrite,'T') + index(clwrite,'t') .gt.0 ) then
            lwrite = .true.
            print *,'flattennc: debug output requested'
        endif

        if ( iargc().ne.2 ) then
            print *,'usage: flattennc infile outfile'
            stop
        endif
        call getarg(1,infile)
        call getarg(2,outfile)
!
!       open files
!
        status = nf_open(infile,nf_nowrite,ncid1)
        if ( status.ne.nf_noerr ) call handle_err(status,infile)
!
!       construct a sensible title
!
        call gettitle(ncid1,title,lwrite)
        call gettextattopt(ncid1,nf_global,'history',history,lwrite)
        call getnumbers(ncid1,ndims1,nvars1,ngatts1,unlimdimid1,lwrite)
        do varid=1,nvars1
*           get name of variable
            status = nf_inq_var(ncid1,varid,name,xtype,ndimvar,dimids
     +            ,natts)
            if ( status.ne.nf_noerr ) call handle_err(status
     +            ,'nf_inq_var')
            if ( name.eq.'source' .or. name.eq.'experiment_id' ) then
!               get lengths
                n = 1
                do i=1,ndimvar
                    status = nf_inq_dim(ncid1,dimids(i),name,m)
                    if ( status.ne.nf_noerr ) call handle_err(status
     +               ,'nf_inq_dim')
                    n = m*n
                enddo
                l = min(len_trim(title),len(title)-2)
                if ( l+n+1.lt.len(title) ) then
                    status = nf_get_var_text(ncid1,varid,title(l+2:))
                    if ( status.ne.nf_noerr ) call handle_err(status
     +                   ,'nf_get_var_text')
                    title(l+n+2:) = ' '
!                   often there are zeros left
                    do i=l+2,l+n+1
                        if ( title(i:i).eq.char(0) ) title(i:i) = ' '
                    enddo
                endif
            endif
        enddo
        if ( lwrite ) print *,'title = ',trim(title)
!
!       read axes
!
        call getdims(ncid1,ndims1,ix,nx,nxmax,iy,ny,nymax,iz,nz,nzmax,it
     +       ,nt,ndata,ie,nens1,nens2,nensmax,lwrite)
        if ( nens1.ne.nens2 ) then
            write(0,*) 'flattennc: error: cannot handle ensembles yet'
            write(0,*) 'flattennc: error: cannot handle ensembles yet'
            call abort
        endif
        ntvars = 0
        undef = 3e33
        foundreftime = .false.
        foundleadtime = .false.
        xx(1) = 0
        yy(1) = 0
        zz(1) = 0
        do varid=1,nvars1
*           get dimensions of variable
            status = nf_inq_var(ncid1,varid,name,xtype,ndimvar,dimids
     +            ,natts)
            if ( status.ne.nf_noerr ) call handle_err(status
     +            ,'nf_inq_var')
            if ( lwrite ) print *,'investigating variable ',varid
     +           ,trim(name),dimids(1:ndimvar)
            if ( index(name,'_bnd').ne.0 ) then
                if ( lwrite ) print *
     +               ,'flattennc: disregarding boundary ',trim(name)
                cycle
            endif
*           what kind of variable do we have?
            if ( ndimvar.eq.1 .and. dimids(1).eq.ix ) then
                call getdiminfo('x',ncid1,varid,xx,nx,lwrite)
                call makelonreasonable(xx,nx)
            elseif ( ndimvar.eq.1 .and. dimids(1).eq.iy ) then
                call getdiminfo('y',ncid1,varid,yy,ny,lwrite)
            elseif ( ndimvar.eq.1 .and. dimids(1).eq.iz ) then
                call getzdiminfo('z',ncid1,varid,zz,nz,lz,lwrite)
            elseif ( ndimvar.eq.1 .and. dimids(1).eq.ie ) then
                if ( lwrite ) print *,'renumbering ensemble members'//
     +               ' to 0 ... nens-1'
            elseif ( name.eq.'reftime' .and. dimids(1).eq.it  ) then
                foundreftime = .true.
                call getreftime(ncid1,varid,tt,nt,firstmo,firstyr
     +                   ,nperyear,iperyear,lwrite)
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
     +                       ,'flattennc: time-varying variable ',varid
                    elseif ( ix.ne.0 .and. dimids(i).eq.ix ) then
                        m = m+1
                        if ( lwrite ) print *
     +                       ,'flattennc: x-dependent variable ',varid
                    elseif ( iy.ne.0 .and. dimids(i).eq.iy ) then
                        m = m+1
                        if ( lwrite ) print *
     +                       ,'flattennc: y-dependent variable ',varid
                    elseif ( iz.ne.0 .and. dimids(i).eq.iz ) then
                        m = m+1
                        if ( lwrite ) print *
     +                       ,'flattennc: z-dependent variable ',varid
                    endif
                enddo
                if ( n.gt.0 .and. m.gt.0 ) then
                    call addonevariable(ncid1,varid,name,ntvars,nvarmax
     +                   ,ndimvar,dimids,ix,iy,iz,it,ie,vars,jvars,lvars
     +                   ,svars,units,cell_methods,undef,lwrite)
                    if ( jvars(4,ntvars).eq.0 ) then
                        ivars(1,ntvars) = 0
                    else
                        ivars(1,ntvars) = nz
                    endif
                endif
            endif
        enddo
!
!       flatten time axis
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
     +       iperyear.eq.12  .and. leadunits.eq.'months' ) then
            nfac = 1
        elseif ( iperyear.eq.366 .and. leadunits.eq.'hours' ) then
            nfac = 24
        else
            write(0,*) 'flattennc: cannot handle iperyear,leadunits = '
     +           ,iperyear,trim(leadunits),' yet'
            write(*,*) 'flattennc: cannot handle iperyear,leadunits = '
     +           ,iperyear,trim(leadunits),' yet'
            call abort
        endif
        lmax = 0
        do it=1,nt
            if ( tl(it)/nfac.lt.iperyear ) then
                tt(it) = tt(it) + tl(it)/nfac
                lmax = max(lmax,nint(tl(it)))
            else
!               my system cannot handle runs of more than 1 year
                if ( lwrite ) print *,'set tt(',it,') to undef, tl='
     +               ,tl(it)
                tt(it) = 3e33
            endif
        enddo
!       get rid of undefineds at the end
        do it=nt,1,-1
            if ( tt(it).lt.1e33 ) exit
        enddo
        nt = it
        if ( it.le.1 ) then
            write(0,*) 'flattennc: error: time axis too short: ',it
            call abort
        endif
!
        if ( tt(2).eq.tt(1) ) then
            write(0,*) 'flattennc: error: tt(1) = tt(2): ',tt(1),tt(2)
            write(*,*) 'flattennc: error: tt(1) = tt(2): ',tt(1),tt(2)
            call abort
        endif
        mperyear = nint(iperyear/(tt(2)-tt(1)))
        if ( mperyear.eq.13 ) mperyear = 12 ! round-off errors in February
        if ( iperyear.ne.366 ) then
            mt = 1 + nint(mperyear*(tt(nt) - tt(1))/iperyear)
        elseif ( mperyear.eq.366 ) then
            mt = 1 + nint(tt(nt) - tt(1))
        else
            mt = 1 + nint(mperyear*(tt(nt) - tt(1))/365.24)
        endif
!       convert starting dates to new mperyear
        if ( mperyear.ne.nperyear ) then
            if ( mperyear.ne.366 ) then
                firstmo = 1 + (firstmo-1)
     +               *max(mperyear,12)/max(nperyear,12)
            elseif ( nperyear.le.12 ) then
                j = 1
                do i=1,firstmo-1
                    j = j + dpm(i)
                enddo
                firstmo = j
            else
                write(0,*) 'flattennc: cannot transform starting date'
                write(*,*) 'flattennc: cannot transform starting date'
                call abort
            endif
        endif
        if ( lwrite ) then
            print *,'flattennc: deduced that mperyear = ',mperyear
            print *,'                        mt       = ',mt
            print *,'                        firstmo  = ',firstmo
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
        write(ltime,'(a,i3,2a)') 'verification time (reftime+leadtime,'
     +       //'leadtime &lt; ',k,trim(lunits),')'
        call getdymo(dy,mo,firstmo,mperyear)
        write(title,'(aa,i2,a)') trim(title),', reftime ',dy,months(mo)
!
!       write header
!
        call enswritenc(outfile,ncid2,ntvarid,itimeaxis,ndata,nx,xx,ny
     +       ,yy,nz,zz,lz,mt,mperyear,firstyr,firstmo,ltime,undef,title
     +       ,history,ntvars,vars,ivars,lvars,svars,units,cell_methods
     +       ,metadata,nens1,nens2)
!
!       read/write data
!
        allocate(data(nx,ny,nz))
        if ( lwrite ) then
            call getdymo(dy,mo,firstmo,mperyear)
            jul0 = julday(mo,dy,firstyr)
        endif
        j = 1
        do i=1,mt
!           start dates are the same, are the offsets OK?
            do while ( tt(j).gt.1e33 )
                j = j + 1
            enddo
            if ( iperyear.eq.366 ) then
                if ( mperyear.eq.366 ) then
                    n = nint(tt(j)-tt(1))
                elseif ( mperyear.lt.360 ) then
                    n = nint((tt(j)-tt(1))*mperyear/365.24)
                else
                    write(0,*) 'flattennc: cannot convert ',iperyear,
     +                  ' to ',mperyear,' yet'
                    write(*,*) 'flattennc: cannot convert ',iperyear,
     +                  ' to ',mperyear,' yet'
                    call abort
                endif
            else
                n = nint((tt(j)-tt(1))*mperyear/iperyear)
            endif
            if ( n.eq.i-1 ) then
                do varid=1,ntvars
                    if ( lwrite ) then
                        call caldat(jul0+nint(tt(j)-tt(1)),mo,dy,yr)
                        print '(a,2i6,i3,a,i4)','read field ',n,j,dy
     +                       ,months(mo),yr
                    endif
                    call readncslice(ncid1,jvars(1,varid),j,data,nx,ny
     +                   ,nz,lwrite)
                    if ( lwrite ) then
                        print '(a,i6,i3,a,i4)','write field',i
                    endif
                    call writencslice(ncid2,0,itimeaxis,ndata,ivars
     +                   ,data,nx,ny,nz,nx,ny,nz,i,0)
                enddo
                j = j + 1
            elseif ( n.gt.i-1 ) then
                do varid=1,ntvars
                    if ( lwrite ) print *,'write field of undefineds',i
                    data = undef
                    call writencslice(ncid2,0,itimeaxis,ndata,ivars
     +                   ,data,nx,ny,nz,nx,ny,nz,i,0)
                enddo
            else                ! n.lt.i-1
                write(0,*) 'flattennc: error: should never happen'
                write(0,*) 'output time step = ',i
                write(0,*) 'input slice      = ',n,j,tt(j)
            endif
        enddo
!       do not forget to close files!  otherwise the last bit is lost
        status = nf_close(ncid1)
        status = nf_close(ncid2)
!
        end
