        program averagefieldspace
*
*       average a few grid boxes together
*       should do this in ferret
*
        implicit none
#include "param.inc"
#include "netcdf.inc"
#include "recfac.h"
        integer i,j,ii,jj,iii,jjj,it,n,mo,yr,jx,jy,avex,avey,status,nxf
     +       ,nyf,nzf
        integer ncid,nx,ny,nz,nt,nperyear,firstyr,firstmo,endian,nvars
     +       ,ivars(2,1),jvars(6,1),lastyr,ntvarid,itimeaxis(ntmax)
     +       ,nnt
        real xx(nxmax),yy(nymax),zz(1),undef,minfac,s,ss
        real,allocatable :: field(:,:,:,:)
        character file*255,datfile*255,title*255,vars(1)*20,lvars(1)*60
     +       ,units(1)*40
        logical lwrite
        integer iargc,llen,leap
*
        if ( iargc().ne.4 ) then
            print *,
     +           'usage: averagefieldspace infield avex avey outfield'
            call abort
        endif

        lwrite = .false.
        minfac = 0.3
        call getarg(2,file)
        read(file,*) avex
        call getarg(3,file)
        read(file,*) avey
        call getarg(1,file)
        if ( lwrite ) print *,'nf_opening file ',trim(file)
        status = nf_open(file,nf_nowrite,ncid)
        if ( status.ne.nf_noerr ) then
            call parsectl(file,datfile,nxmax,nx,xx,nymax,ny,yy,nzmax,nz
     +            ,zz,nt,nperyear,firstyr,firstmo,undef,endian,title,1
     +            ,nvars,vars,ivars,lvars,units)
            ncid = -1
        else
            call parsenc(file,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax
     +            ,nz,zz,nt,nperyear,firstyr,firstmo,undef,title,1,nvars
     +            ,vars,jvars,lvars,units)
            datfile = file
        endif
*       range of years
        lastyr = firstyr + (firstmo+nt-2)/nperyear
*
*       allocate field
*
        nxf = nx
        nyf = ny
        nzf = nz
        if ( lwrite ) print *,'allocating field ',nx,ny,nperyear,firstyr
     +       ,lastyr
        allocate(field(nx,ny,nperyear,firstyr:lastyr))
*
*       read field
*
        if ( ncid.eq.-1 ) then
            call readdatfile(datfile,field,nx,ny,nx,ny,nperyear,firstyr
     +            ,lastyr,firstyr,firstmo,nt,undef,endian,lwrite,firstyr
     +            ,lastyr,1,1)
        else
            call readncfile(ncid,field,nx,ny,nx,ny,nperyear,firstyr
     +            ,lastyr,firstyr,firstmo,nt,undef,lwrite,firstyr,lastyr
     +            ,jvars)
            status = nf_close(ncid)
        endif
*
*       average
*
        call spatialaverage(field,xx,yy,nx,ny,nperyear,firstyr,lastyr,
     +       avex,avey,lwrite)
*
*       write output
*
        nnt = nt
        if ( nperyear == 366 ) then
            nnt = 0
            yr = firstyr
            mo = firstmo - 1
            do it=1,nt
                nnt = nnt + 1
                mo = mo + 1
                if ( mo.gt.nperyear ) then
                    mo = mo - nperyear
                    yr = yr + 1
                endif
                if ( mo.eq.60 ) then
                    if ( leap(yr).eq.1 ) then
                        nnt = nnt - 1
                    end if
                endif
            end do
        end if
        undef = 3e33
        write(title(llen(title)+2:),'(a,i2,a,i2,a)')
     +       ' averaged over ',avex,'x',avey,' grid boxes'
        call getarg(4,file)
        i = index(file,'.ctl')
        if ( i.ne.0 ) then
            datfile = file(1:i)
            datfile(i:) = '.grd'
            call writectl(file,datfile,nx,xx,ny,yy,nz,zz,
     +        nnt,nperyear,firstyr,firstmo,undef,title,nvars,vars,ivars,
     +        lvars,units)
            open(1,file=datfile,access='direct',recl=recfa4*nx*ny)
            yr = firstyr
            mo = firstmo
            do it=1,nt
                write(1,rec=it) ((field(i,j,mo,yr),i=1,nx),j=1,ny)
                mo = mo + 1
                if ( mo.gt.nperyear ) then
                    mo = mo - nperyear
                    yr = yr + 1
                endif
                if ( nperyear.eq.366 .and. mo.eq.60 ) then
                    if ( leap(yr).eq.1 ) then
                        mo = mo + 1
                    end if
                endif
            enddo
        elseif ( index(file,'.nc').ne.0 ) then
            ivars = 0
            call writenc(file,ncid,ntvarid,itimeaxis,ntmax,nx,xx,ny,yy
     +           ,nz,zz,nnt,nperyear,firstyr,firstmo,undef,title,nvars
     +           ,vars,ivars,lvars,units,0,0)
            open(1,file=datfile,access='direct',recl=recfa4*nx*ny)
            yr = firstyr
            mo = firstmo
            do it=1,nnt
                if ( lwrite ) print *,'writing ',yr,mo,it,nnt,nt
                call keepalive1('Writing time step ',it,nnt)
                call writencslice(ncid,ntvarid,itimeaxis,ntmax,ivars
     +               ,field(1,1,mo,yr),nxf,nyf,nzf,nx,ny,nz,it,0)
                mo = mo + 1
                if ( mo.gt.nperyear ) then
                    mo = mo - nperyear
                    yr = yr + 1
                endif
                if ( nperyear.eq.366 .and. mo.eq.60 ) then
                    if ( leap(yr).eq.1 ) then
                        mo = mo + 1
                    end if
                endif
            enddo
            status = nf_close(ncid) ! otherwise the buffer won't be flushed...
        else
            write(0,*) 'need .ctl or .nc in output file name'
            call abort
        endif
        end
