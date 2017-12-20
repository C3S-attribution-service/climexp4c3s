        program del_dimension
!
!       flatten a netcdf file with ECMWF conventions to one
!       that the Climate Explorer can handle, with one time 
!       axis without holes.
!
        implicit none
        integer nvarmax,ntmax
        parameter (nvarmax=20,ntmax=60000)
#include "params.h"
#include "netcdf.inc"
        integer i,j,k,l,m,n,dy,mo,yr,lead,dyref,moref,yrref,status,varid
     +       ,xtype,ndimvar,dimids(nf_max_var_dims),natts
        integer nx,ny,nz,nens1,nens2,nt,mt
        integer ncid1,ndims1,nvars1,ngatts1,unlimdimid1
        integer ncid2,ntvarid,itimeaxis(ndata)
        integer ix,iy,iz,it,ie,ntvars,firstmo,firstyr,nperyear,iperyear,
     +       mperyear,ivars(6,nvarmax),lmax,dpm(12),
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
        logical lwrite,tdefined(ntmax)
        integer iargc,julday
        data months 
     +       /'jan','feb','mar','apr','may','jun'
     +       ,'jul','aug','sep','oct','nov','dec'/
        data dpm /31,29,31,30,31,30,31,31,30,31,30,31/

        lwrite = .false.
        call getenv('DEL_DIMENSION_LWRITE',clwrite)
        if ( index(clwrite,'T') + index(clwrite,'t') .gt.0 ) then
            lwrite = .true.
            print *,'del_dimension: debug output requested'
        endif

        if ( iargc().ne.2 ) then
            print *,'usage: del_dimension infile outfile'
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
!       read axes
!
        call ensparsenc(infile,ncid1,nxmax,nx,xx,nymax,ny,yy,nzmax
     +       ,nz,zz,lz,nt,nperyear,firstyr,firstmo,ltime,tdefined,ntmax
     +       ,nens1,nens2,undef,title,history,nvarmax,ntvars,vars,ivars
     +       ,lvars,svars,units,cell_methods,metadata)
!
!       write header
!
        call enswritenc(outfile,ncid2,ntvarid,itimeaxis,ndata,nx,xx,ny
     +       ,yy,nz,zz,lz,nt,nperyear,firstyr,firstmo,ltime,undef,title
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
        do i=1,nt
            do varid=1,ntvars
                call readncslice(ncid1,ivars(1,varid),i,data,nx,ny
     +               ,nz,lwrite)
                if ( lwrite ) then
                    print '(a,i6,i3,a,i4)','read field',i
                endif
                call writencslice(ncid2,0,itimeaxis,ndata,ivars
     +               ,data,nx,ny,nz,nx,ny,nz,i,0)
            enddo
        enddo
!       do not forget to close files!  otherwise the last bit is lost
        status = nf_close(ncid1)
        status = nf_close(ncid2)
!
        end
