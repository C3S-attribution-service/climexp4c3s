        program get_depth
*
*       compute the depth of the x-degree isotherm or thermocline
*       or mixed layer from 3D ocean data
*
        implicit none
#include "params_s.h"
#include "netcdf.inc"
        logical lwrite,ensemble
        integer nens1,nens2
        common /options/ lwrite,ensemble,nens1,nens2
        integer nyrmax,nlevmax,mensmax
        parameter (nyrmax=60,nlevmax=3,mensmax=1)
        integer n,nx,ny,nz,nt,firstyr,lastyr,firstmo,nvars,
     +          ivars(2,nvmax),jvars(6,nvmax),ncid,endian,
     +          status,nperyear,mens
        logical lexist
        real xx(nxmax),yy(nymax),zz(nzmax),undef,
     +        field(nxmax,nymax,12,nyrmax),
     +        res(nxmax,nymax,12,nyrmax)
        character infile*255,datfile*255,outfile*255,line*255
     +        ,vars(nvmax)*10,lvars(nvmax)*40,title*255,
     +        units(nvmax)*10
        integer iargc,llen
*
*       for f2c
        common /ccorrelatefield/ field,res
*
*       process command line
*
        lwrite = .FALSE.
        n = iargc()
        if ( n.lt.2 ) then
            write(0,*) 'usage: get_{depth|heat|thermocline|mld}'//
     +           'infile.[ctl|nc] {temp|heat|delta_t} outfile.[ctl|nc]'
            stop
        endif
        call getarg(1,infile)
        if ( index(infile,'%').gt.0 .or. index(infile,'++').gt.0 ) then
            ensemble = .true.
            call filloutens(infile,0)
        else
            mens = 0
        endif
        if ( lwrite ) print *,'get_depth: nf_opening file '
     +        ,infile(1:llen(infile))
        status = nf_open(infile,nf_nowrite,ncid)
        if ( status.ne.nf_noerr ) then
            call parsectl(infile,datfile,nxmax,nx,xx,nymax,ny,yy,nzmax
     +            ,nz,zz,nt,nperyear,firstyr,firstmo,undef,endian,title
     +            ,1,nvars,vars,ivars,lvars,units)
            ncid = -1
            if ( ensemble ) then
                do mens=1,nensmax
                    call getarg(1,line)
                    call filloutens(line,mens)
                    inquire(file=line,exist=lexist)
                    if ( .not.lexist ) goto 100
                enddo
  100           continue
                mens = mens - 1
                write(0,*) 'located ',mens+1,' ensemble members<br>'
            endif
        else
            call parsenc(infile,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax
     +            ,nz,zz,nt,nperyear,firstyr,firstmo,undef,title,1,nvars
     +            ,vars,jvars,lvars,units)
            if ( ensemble ) then
                do mens=1,nensmax
                    call getarg(1,line)
                    call filloutens(line,mens)
                    status = nf_open(line,nf_nowrite,ncid)
                    if ( status.ne.nf_noerr ) goto 200
                enddo
  200           continue
                mens = mens - 1
                write(0,*) 'located ',mens+1,' ensemble members<br>'
            endif
        endif
        lastyr = firstyr + (firstmo+nt-2)/nperyear
*       process arguments
        if ( ensemble ) write(0,*) 'Using ensemble members ',nens1
     +        ,' to ',nens2,'<br>'
        call getarg(1,infile)
*
*       check dimensions
*
        if ( nx*ny*nz*nperyear*(lastyr-firstyr+1)
     +        .gt.nxmax*nymax*12*nyrmax ) then
            write(0,*) 'get_deoth: error: field too large '
            write(0,*) '  nx       = ',nx
            write(0,*) '  ny       = ',ny
            write(0,*) '  nz       = ',nz
            write(0,*) '  nperyear = ',nperyear
            write(0,*) '  years    = ',firstyr,lastyr
            write(0,*) 'total request',
     +            nx*ny*nz*nperyear*(lastyr-firstyr+1)
            write(0,*) 'available    ',
     +            nxmax*nymax*12*nyrmax
            write(*,*) 'get_depth: error: field too large '
            write(*,*) '  nx       = ',nx
            write(*,*) '  ny       = ',ny
            write(*,*) '  nz       = ',nz
            write(*,*) '  nperyear = ',nperyear
            write(*,*) '  years    = ',firstyr,lastyr
            write(*,*) 'total request',
     +            nx*ny*nz*nperyear*(lastyr-firstyr+1)
            write(*,*) 'available    ',
     +            nxmax*nymax*12*nyrmax
            call abort
        endif
        if ( nx*ny*nz .gt. nxmax*nymax*nlevmax ) then
            write(0,*) 'get_depth: error: fields too large: '
            write(0,*) 'nx,ny,nz            = ',nx,ny,nz,nx*ny*nz
            write(0,*) 'nxmax,nymax,nlevmax = ',nxmax,nymax,nlevmax
     +           ,nxmax*nymax*nlevmax
            write(*,*) 'get_depth: error: fields too large: '
            write(*,*) 'nx,ny,nz            = ',nx,ny,nz,nx*ny*nz
            write(*,*) 'nxmax,nymax,nlevmax = ',nxmax,nymax,nlevmax
     +           ,nxmax*nymax*nlevmax
            call abort
        endif
*
*       to save on RAM usage
*
        call gfield(datfile,ncid,field,res,nx,xx,ny,yy,nz,zz,nt
     +        ,nperyear,firstyr,lastyr,firstmo,undef,endian,jvars)
        end

        subroutine gfield(datfile,ncid,field,res,nx,xx,ny,yy,nz,zz,nt
     +        ,nperyear,firstyr,lastyr,firstmo,undef,endian,jvars)
*
*       break to use field() compactly
*
        implicit none
#include "recfac.h"
#include "params.h"
        logical lwrite,ensemble
        integer nens1,nens2
        common /options/ lwrite,ensemble,nens1,nens2
        real absent
        parameter (absent=3e33)
*
        integer ncid,endian,nx,ny,nz,nt,nperyear,firstyr,lastyr,
     +          firstmo,jvars(6,nvmax)
        real field(nx,ny,nz,nperyear,firstyr:lastyr),
     +          res(nx,ny,nperyear,firstyr:lastyr),
     +          undef,xx(nx),yy(ny),zz(nz)
        character datfile*(*)
*
        integer jx,jy,jz,jz0,i,j,k,k1,k2,n,mo,yr,itype,ldir,
     +          nvars,ivars(2,nvmax),iens
        real val,depth(nzmax),gradients(nzmax),gradmax,temp,tmin,tmax,dt
     +       ,grad,z
        logical lexist
        character outfile*255,line*255,yesno*1,dir*255,string*10
     +       ,vars(nvmax)*10,lvars(nvmax)*40,units(nvmax)*20,title*255
        integer iargc,llen,rindex
*
        if ( nz.le.4 ) then
            write(0,*) 'error: cannot get depth properties with nz = '
     +           ,nz
            write(*,*) 'error: cannot get depth properties with nz = '
     +           ,nz
            call abort
        endif
        call getarg(0,line)
        if ( index(line,'get_depth').ne.0 ) then
            itype = 1
        elseif ( index(line,'get_heat').ne.0 ) then
            itype = 2
        elseif ( index(line,'get_therm').ne.0 ) then
            itype = 0
        elseif ( index(line,'get_mld').ne.0 ) then
            itype = 3
        else
            goto 901
        endif
        if ( itype.gt.0 ) then
            call getarg(2,line)
            read(line,*,err=902) val
        endif
*
*       read field, change absent values to our convention
*
        if ( ensemble ) then
*           put the %% back in datfile...
            if ( nens2.lt.10 ) then
                i = 1
            elseif ( nens2.lt.100 ) then
                i = 2
            elseif ( nens2.lt.1000 ) then
                i = 3
            else
                write(0,*)
     +                'get_depth: cannot handle ensembles up to '
     +                ,nens2,' yet'
                call abort
            endif
            string = '0000000'
            j = rindex(datfile,string(1:i))
            if ( j.eq.0 ) then
                write(0,*) 'get_depth: error: cannot find '
     +                ,string(1:i),' in ',datfile(1:llen(datfile))
                call abort
            endif
            do k=j,j+i-1
                datfile(k:k) = '%'
            enddo
        endif
        do iens=nens1,nens2
            call keepalive(iens-nens1+1,nens2-nens1+1)
            if ( ncid.eq.-1 ) then
                dir=datfile
                if ( ensemble ) call filloutens(dir,iens)
                print *,'looking for '//dir(1:llen(dir))
                inquire(file=dir,exist=lexist)
                if ( .not.lexist ) then
                    print *,'looking for '//dir(1:llen(dir))//'.gz'
                    inquire(file=dir(1:llen(dir))//'.gz',exist=lexist)
                    if ( .not.lexist ) then
                        nens2 = iens-1
                        if ( nens2.ge.nens1 ) then
                            write(0,*) 'Found ensemble 0 to ',nens2
     +                            ,'<br>'
                            goto 5
                        else
                            write(0,*) 'Cannot locate file '
     +                            ,dir(1:llen(dir))
                            call abort
                        endif
                    endif
                endif
                if ( lwrite ) then
                    print *,'opening file ',dir(1:llen(dir))
                endif
                call zreaddatfile(dir,field(1,1,1,1,firstyr),
     +                nx,ny,nz,nx,ny,nz,nperyear,firstyr,lastyr,
     +               firstyr,firstmo,nt,undef,endian,lwrite,firstyr
     +               ,lastyr,1,1)
            else
                if ( nz.ne.1 ) then
                    write(0,*) 'cannot read 3D netCDF files yet'
                    call abort
                endif
                if ( ensemble ) then
                    write(0,*)
     +                    'cannot handle ensembles of netcdf files yet'
                    call abort
                endif
                call readncfile(ncid,field,nx,ny,nx,ny,nperyear,firstyr
     +               ,lastyr,firstyr,firstmo,nt,undef,lwrite,firstyr
     +               ,lastyr,jvars)
            endif
 5          continue
*           open output file
*
            call getarg(iargc(),outfile)
            if ( ensemble ) call filloutens(outfile,iens)
            inquire(file=outfile,exist=lexist)
            if ( lexist ) then
                print *,'output file ',outfile(1:index(outfile,' ')-1),
     +               ' already exists, overwrite? [y/n]'
                read(*,'(a)') yesno
                if (  yesno.ne.'y' .and. yesno.ne.'Y' .and. 
     +               yesno.ne.'j' .and. yesno.ne.'J' ) then
                    stop
                endif
                open(2,file=outfile)
                close(2,status='delete')
            endif
            if ( index(outfile,'.nc').ne.0 ) then
                print *,'netCDF output not yet ready'
                stop
            else
                i = index(outfile,'.ctl')
                if ( i.ne.0 ) then
                    datfile = outfile(:i-1)//'.dat'
                else
                    datfile = outfile
                endif
                open(unit=2,file=datfile,form='unformatted',access
     +               ='direct',recl=recfac*nx*ny,err=920)
            endif
*
*           loop over time, gridpoints
*
            yr = firstyr
            mo = firstmo - 1
            do i=1,nt
                call keepalive(i,nt)
                mo = mo + 1
                if ( mo.gt.nperyear ) then
                    mo = mo - nperyear
                    yr = yr + 1
                endif
                do jy=1,ny
                    do jx=1,nx
                        if ( itype.eq.0 ) then
*                           find thermocline - steepest gradient
                            tmin = 0
                            tmax = 1000
                            dt = 50
                            n = 1
                            grad = -1
                            do jz=1,nz-1
                                if ( field(jx,jy,jz,mo,yr).lt.1e33 )
     +                               then
                                    do k=jz+1,nz
                                        if ( field(jx,jy,k,mo,yr)
     +                                       .lt.1e33 ) then
                                            gradients(n) = 
     +                                           (field(jx,jy,jz,mo,yr)
     +                                           - field(jx,jy,k,mo,yr))
     +                                           /(zz(k) - zz(jz))
                                            if ( abs(gradients(n)-grad)
     +                                           .gt.0.01 ) then
                                                grad = gradients(n)
                                                depth(n) = (zz(jz) +
     +                                               zz(k))/2
                                                n = n + 1
                                                jz0 = jz
                                            else
                                                depth(n-1) = (zz(jz0) +
     +                                               zz(k))/2
                                            endif
                                            goto 110
                                        endif
                                    enddo
 110                                continue
                                endif
                            enddo
                            n = n-1
                            if ( n.gt.2 ) then
                                call maxquad(res(jx,jy,mo,yr),gradmax
     +                               ,depth,gradients,n)
                                if ( res(jx,jy,mo,yr).lt.1e33 ) then
                                    z = res(jx,jy,mo,yr)
                                    do jz=1,n-1
                                        if ( (zz(jz)-z)*(zz(jz+1)-z)
     +                                       .lt.0 ) then
                                            do k1=jz,1,-1
                                                if ( field(jx,jy,k1,mo
     +                                               ,yr).lt.1e33 )
     +                                               goto 120
                                            enddo
                                            goto 190
 120                                        continue
                                            do k2=jz+1,nz
                                                if ( field(jx,jy,k2,mo
     +                                               ,yr).lt.1e33 )
     +                                               goto 130
                                            enddo
                                            goto 190
 130                                        continue
                                            temp = (
     +                                           field(jx,jy,k1,mo,yr)
     +                                           *(zz(k2)-z) + 
     +                                           field(jx,jy,k2,mo,yr)
     +                                           *(z-zz(k1)))/
     +                                           (zz(k2)-zz(k1))
                                            if ( temp.lt.15 .or. 
     +                                           temp.gt.27 ) then
                                                print *,'weird thermo'//
     +                                               'cline at ',xx(jx)
     +                                               ,yy(jy),jz,jz+1,mo
     +                                               ,yr,z,temp,field(jx
     +                                               ,jy,k1,mo,yr),k1
     +                                               ,field(jx,jy,k2,mo
     +                                               ,yr),k2
                                                do k=1,n
                                                    print *,k,depth(k)
     +                                                   ,gradients(k)
                                                enddo
                                            endif
                                            goto 190
                                        endif
                                    enddo
 190                                continue
                                endif
                            else
                                res(jx,jy,mo,yr) = 3e33
                            endif
                        elseif ( itype.eq.1 ) then
*                           find val-degree isotherm
                            tmin=0
                            tmax=6000
                            dt = 10
                            res(jx,jy,mo,yr) = 3e33
                            do jz=1,nz-1
                                if ( field(jx,jy,jz,mo,yr).lt.1e33 .and.
     +                            field(jx,jy,jz+1,mo,yr).lt.1e33 ) then
                                    if ( (field(jx,jy,jz,mo,yr)-val)*
     +                                   (field(jx,jy,jz+1,mo,yr)-val)
     +                                   .le.0 ) then
                                        res(jx,jy,mo,yr) = (
     +                                       (field(jx,jy,jz,mo,yr)-val)
     +                                       *zz(jz+1) + 
     +                                       (val-field(jx,jy,jz+1,mo,yr
     +                                       ))*zz(jz))/
     +                                       (field(jx,jy,jz,mo,yr)
     +                                       -field(jx,jy,jz+1,mo,yr))
                                    endif
                                endif
                            enddo
                        elseif ( itype.eq.2 ) then
*                           find heat content of top val meters
                            tmin=val*3
                            tmax=val*30
                            dt = 50
                            res(jx,jy,mo,yr) = 0
                            do jz=2,nz-1
                                if ( field(jx,jy,jz,mo,yr).lt.1e33.and.
     +                               res(jx,jy,mo,yr).le.1e33 .and.
     +                               zz(jz).lt.val ) then
                                    if ( jz.eq.1 ) then
                                        res(jx,jy,mo,yr) = 
     +                                       field(jx,jy,jz,mo,yr)*
     +                                       (zz(1)+zz(2))/2
                                    else
                                        res(jx,jy,mo,yr) = 
     +                                       res(jx,jy,mo,yr) + 
     +                                       field(jx,jy,jz,mo,yr)*
     +                                       (zz(jz+1)-zz(jz-1))/2
                                    endif
                                elseif ( zz(jz).lt.val ) then
                                    res(jx,jy,mo,yr) = 3e33
                                endif
                            enddo
                        elseif ( itype.eq.3 ) then
                            tmin=0
                            tmax=5000
                            dt = 10
*                           find mixed layer depth with a val-degree
*                           criterium
                            res(jx,jy,mo,yr) = 3e33
                            if ( field(jx,jy,1,mo,yr).lt.1e33 ) then
                                temp = field(jx,jy,1,mo,yr) - val
                                do jz=1,nz-1
                                    if ( field(jx,jy,jz,mo,yr).lt.1e33
     +                                   .and.
     +                                   field(jx,jy,jz+1,mo,yr).lt.1e33
     +                                   ) then
                                        if ( (field(jx,jy,jz,mo,yr)-temp
     +                                       )*(field(jx,jy,jz+1,mo,yr)
     +                                       -temp).le.0 ) then
                                            res(jx,jy,mo,yr) = (
     +                                           (field(jx,jy,jz,mo,yr)
     +                                           -temp)*zz(jz+1) +(temp
     +                                           -field(jx,jy,jz+1,mo,yr
     +                                           ))*zz(jz))/(field(jx,jy
     +                                           ,jz,mo,yr) -field(jx,jy
     +                                           ,jz+1,mo,yr))
                                        endif
                                    endif
                                enddo
                            endif
                            if ( res(jx,jy,mo,yr).le.0 ) then
                                print *,'error: mld(',xx(jx),yy(jy),mo
     +                               ,yr,') = ',res(jx,jy,mo,yr)
                                do jz=1,nz
                                    if ( field(jx,jy,jz,mo,yr).lt.1e33
     +                                   ) then
                                        print *,jz,field(jx,jy,jz,mo,yr)
                                    endif
                                enddo
                            endif
                        else
                            write(0,*) 'error: unknown itype ',itype
                            call abort
                        endif
                        if ( lwrite ) then
                            if ( res(jx,jy,mo,yr).lt.1e33 ) then
                                print *,'res(',jx,jy,mo,yr,') = ',
     +                               res(jx,jy,mo,yr)
                            endif
                        endif
                    enddo       ! nx
                enddo           ! ny
                call latlonint(res(1,1,mo,yr),1,1,ny,nx,yy,xx,1,tmin
     +               ,tmax,dt)
                call latlonint(res(1,1,mo,yr),1,1,ny,nx,yy,xx,1,tmin
     +               ,tmax,dt)
            enddo               ! nt
*
*           write output field in GrADS format
*
            print *,'writing output'
            yr = firstyr
            mo = firstmo - 1
            do i=1,nt
                call keepalive(i,nt)
                mo = mo + 1
                if ( mo.gt.nperyear ) then
                    mo = mo - nperyear
                    yr = yr + 1
                endif
                write(2,rec=i) ((res(jx,jy,mo,yr),jx=1,nx),jy=1,ny)
            enddo
            close(2)
            if ( index(outfile,'.ctl').ne.0 ) then
                call getenv('DIR',dir)
                ldir = llen(dir)
                if ( ldir.eq.0 ) ldir=1
                if ( dir(ldir:ldir).ne.'/' ) then
                    ldir = ldir + 1
                    dir(ldir:ldir) = '/'
                endif
                title = ' '
                n = 1
                do i=0,iargc()-1
                    call getarg(i,line)
                    if ( line(1:ldir).eq.dir(1:ldir) ) then
                        title(n:) = line(ldir+1:)
                    else
                        title(n:) = line
                    endif
                    n = llen(title) + 2
                enddo
                nvars = 1
                if ( itype.eq.0 ) then
                    vars(1) = 'zth'
                    lvars(1) = 'thermocline depth [m]'
                elseif ( itype.eq.1 ) then
                    if ( abs(val-real(nint(val))).lt.0.01 ) then
                        write(vars(1),'(a,i2.2)') 'z',nint(val)
                        write(lvars(1),'(f4.0,a)') val
     +                       ,'-degree isotherm [m]'
                    elseif ( abs(10*val-real(nint(10*val))).lt.0.01 )
     +                       then
                        write(vars(1),'(a,i3.3)') 'z',nint(10*val)
                        write(lvars(1),'(f5.1,a)') val
     +                       ,'-degree isotherm [m]'
                    else
                        write(vars(1),'(a,i3.3)') 'z',nint(100*val)
                        write(lvars(1),'(f8.3,a)') val
     +                       ,'-degree isotherm [m]'
                    endif
                elseif ( itype.eq.2 ) then
                    vars(1) = 'heat'
                    write(lvars(1),'(a,i5,a)')
     +                   'heat contenty of the top ',nint(val)
     +                   ,' meter [Cm]'
                elseif ( itype.eq.3 ) then
                    vars(1) = 'mld'
                    write(lvars(1),'(a,f5.2,a)')
     +                   'mixes layer depth (dT=',val,'C) [m]'
                else
                    write(0,*) 'get_depth: error: itype = ',itype
                    call abort
                endif
                ivars(1,1) = 0
                ivars(2,1) = 99
                call writectl(outfile,datfile,nx,xx,ny,yy,1,zz
     +               ,nt,nperyear,firstyr,firstmo,absent,title,nvars
     +               ,vars,ivars,lvars,units)
            endif
        enddo
*
*       error messages
*
        goto 999
 901    print *,line,' not known'
        call abort
 902    print *,'get_depth: error reading val from ',line(1:llen(line))
        call abort
 903    print *,'error reading date from file ',line(1:index(line,' ')-1
     +        ),' at record ',k
        call abort
 904    print *,'error cannot locate field file file ',line(1:index(line
     +        ,' ')-1)
        call abort
 920    print *,'error cannot open new file '
     +        ,datfile(1:index(datfile,' ')-1)
        call abort
 999    continue
        end
