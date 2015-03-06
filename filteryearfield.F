        program filteryearfield
*
*       Apply hi- or lo-pass filters to a field
*       
        implicit none
#include "params.h"
#include "recfac.h"
#include "netcdf.inc"
#include "getopts.inc"
        integer ntmax
        parameter(ntmax=500*npermax)
        integer i,j,yr,mo,nyr,nperyear,ncid,nx,ny,nz,firstyr,firstmo,nt
     +       ,nvars,ivars(2,nvmax),jvars(6,nvmax),endian,lastyr,jx,jy,jz
     +       ,iu,status,mens1,mens,n,ntvarid,itimeaxis(ntmax)
        real xx(nxmax),yy(nymax),zz(nzmax),undef
        real,allocatable :: data(:,:,:,:,:),fxy(:,:)
        character file*256,line*128,hilo*2,filtertype*12,vars(nvmax)*20
     +       ,lvars(nvmax)*256,units(nvmax)*20,title*255,datfile*256
     +       ,prog*100,yearmonth*5,scale*2
        integer iargc,llen,leap

        call getarg(0,prog)
        if ( iargc().lt.5 ) then
            write(0,*) 'usage: ',trim(prog),' hi|lo filtertype nyr '//
     +           'infile outfile'
            stop
        endif
        if ( index(prog,'filteryearfield').ne.0 ) then
            yearmonth = 'year'
        else if ( index(prog,'filtermonthfield').ne.0 ) then
            yearmonth = 'month'
        else
            write(0,*) trim(prog),': error: cannot identify program'
            call abort
        end if
        call getarg(1,hilo)
        if ( hilo.ne.'hi' .and. hilo.ne.'lo' ) then
            write(0,*) 'filterfield: error: say hi or lo, not ',hilo
            call abort
        endif
        call getarg(2,filtertype)
        call getarg(3,line)
        read(line,*,err=901) nyr
        call getarg(4,file)
        status = nf_open(file,nf_nowrite,ncid)
        if ( status.ne.nf_noerr ) then
            ncid = -1
            call parsectl(file,datfile,nxmax,nx,xx,nymax,ny,yy
     +           ,nzmax,nz,zz,nt,nperyear,firstyr,firstmo,undef,endian
     +           ,title,1,nvars,vars,ivars,lvars,units)
        else
            call parsenc(file,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax
     +           ,nz,zz,nt,nperyear,firstyr,firstmo,undef,title,1,nvars
     +           ,vars,jvars,lvars,units)
        endif
        call getopts(5,iargc()-1,nperyear,yrbeg,yrend,.true.,mens1,mens
     +       )
        if ( nperyear.ne.366 ) then
            lastyr = firstyr + (nt+firstmo-2)/nperyear
        else
            lastyr = firstyr + int((nt+firstmo-2)/365.25)
        end if
        allocate(data(nx,ny,nz,nperyear,firstyr:lastyr))
        allocate(fxy(nperyear,firstyr:lastyr))
        call keepalive(0,1)
        if ( ncid.eq.-1 ) then
            call zreaddatfile(datfile,data,nx,ny,nz,nx,ny,nz
     +           ,nperyear,firstyr,lastyr,firstyr,firstmo,nt,undef
     +           ,endian,lwrite,firstyr,lastyr,1,1)
        else
            if ( nz.ne.1 ) then
                write(*,*) 'filteryearfield: error: '//
     +               'cannot handle 3D netcdf file yet'
                write(0,*) 'filteryearfield: error: '//
     +               'cannot handle 3D netcdf file yet'
                call abort
            endif
            ivars(1,1) = 0
            ivars(2,1) = 99
            call readncfile(ncid,data,nx,ny,nx,ny,nperyear,firstyr
     +           ,lastyr,firstyr,firstmo,nt,undef,lwrite,firstyr,lastyr
     +           ,jvars)
        endif
        call keepalive(1,1)
*
        do jz=1,nz
            do jy=1,ny
                do jx=1,nx
                    do yr=firstyr,lastyr
                        do mo=1,nperyear
                            fxy(mo,yr) = data(jx,jy,jz,mo,yr)
                        end do
                    end do
                    if ( filtertype.eq.'running-mean' .or. 
     +                   filtertype(1:3).eq.'box' ) then
                        if ( hilo.eq.'hi' ) then
                            if ( yearmonth.eq.'year' ) then
                                call hipass(fxy,nperyear,nperyear
     +                               ,firstyr,lastyr,nyr,minfac)
                            else
                                call mhipass(fxy,nperyear,nperyear
     +                               ,firstyr,lastyr,nyr-1,minfac)
                            end if
                        else
                            if ( yearmonth.eq.'year' ) then
                                call ndiffit(fxy,nperyear,nperyear
     +                               ,firstyr,lastyr,-nyr+1,minfacsum)
                                call shiftseries(fxy,nperyear,nperyear
     +                               ,firstyr,lastyr,
     +                               -nperyear*((nyr-1)/2))
                            else
                                call sumit(fxy,nperyear,nperyear
     +                               ,firstyr,lastyr,nyr,'v')
                                call shiftseries(fxy,nperyear,nperyear
     +                               ,firstyr,lastyr,(nyr-1)/2)
                            end if
                        end if
                    else if ( filtertype.eq.'loess1' .or. 
     +                       filtertype.eq.'loess2' )then
                        call myloess(fxy,nperyear,nperyear,firstyr
     +                       ,lastyr,nyr/2,minfac,filtertype,hilo
     +                       ,yearmonth,'gaussian',lwrite)
                    else
                        write(0,*) 'filterseries: error: filtertype '
     +                       ,filtertype,' not yet implemeneted'
                        call abort
                    end if
*       
                    do yr=firstyr,lastyr
                        do mo=1,nperyear
                            data(jx,jy,jz,mo,yr) = fxy(mo,yr)
                        enddo
                    enddo
                enddo
                call keepalive((jz-1)*ny + jy,nz*ny)
            enddo
        enddo
*
        call getarg(iargc(),file)
        if ( yearmonth.eq.'year' ) then
            scale = 'yr'
        else
            select case (nperyear)
            case (4)
                scale='se'
            case(12)
                scale = 'mo'
            case(36)
                scale = 'de'
            case(360:366)
                scale = 'dy'
            case default
                scale='pe'    
            end select
        end if
        write(lvars(1),'(2a,i3,7a)') trim(lvars(1)),' with a ',nyr
     +       ,'-',scale,' ',trim(filtertype),' ',hilo,'-pass filter'
        undef = 3e33
        i = index(file,'.ctl')
        if ( i.ne.0 ) then
            ! grads ct/grd output
            datfile = file
            datfile(i:) = '.grd'
            call writectl(file,datfile,nx,xx,ny,yy,nz,zz,
     +           nt,nperyear,firstyr,firstmo,undef,title,nvars,vars,
     +           ivars,lvars,units)
            call rsunit(iu)
            open(iu,file=trim(datfile),status='new',access='direct',
     +           recl=recfa4*nx*ny*nz)
            yr = firstyr
            mo = firstmo
            n = 0
            do i=1,nt
                if ( nperyear.eq.366 .and. mo.eq.60 .and. leap(yr).eq.1
     +               ) then
                    if ( lwrite ) print *,'skipping Feb 29 in non-leap'
     +                   //' year'
                else
                    n = n + 1
                    write(iu,rec=n) (((data(jx,jy,jz,mo,yr),jx=1,nx),
     +                   jy=1,ny),jz=1,nz)
                end if
                mo = mo + 1
                if ( mo.gt.nperyear ) then
                    mo = mo - nperyear
                    yr = yr + 1
                endif
            enddo
            close(iu)
            nt = n
        else
            ! netcdf output
            call writenc(file,ncid,ntvarid,itimeaxis,ntmax,nx,xx,ny,yy,
     +           nz,zz,nt,nperyear,firstyr,firstmo,undef,title,nvars,
     +           vars,jvars,lvars,units,0,0)
            do i=1,nt
                call writencslice(ncid,ntvarid,itimeaxis,ntmax,jvars,
     +               data,nx,ny,nz,nx,ny,nz,i,0)
            end do
            status = nf_close(ncid)
        end if
*
        goto 999
  901   write(0,*) 'filterseries: expcting an integer, not ',file
        call abort
  999   continue
        end
