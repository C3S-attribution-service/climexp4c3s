        program daily2longerfield
*
*       compute aggregate quantities from daily data
*       input: daily field series
*       output: yearly/monthly/10-dy/5-dy field
*
        implicit none
#include "params.h"
#include "recfac.h"
#include "netcdf.inc"
#include "getopts.inc"
        integer ntmax
        parameter(ntmax=60000)
        integer nperyear,nperyearnew,yr,mo,dy,i,j,n,itype,jx,jy
        integer nx,ny,nz,nt,firstyr,firstmo,lastyr,nvars,
     +       ivars(2,nvmax),jvars(6,nvmax),ncid,endian,status,
     +       ntvarid,itimeaxis(ndata),it,ntp,unit
        integer,allocatable :: nn(:,:)
        real pcut,s
        real,allocatable :: cut(:,:,:),normfield(:,:,:),oldseries(:,:)
        real xx(nxmax),yy(nymax),zz(nzmax),undef
        real,allocatable :: oldfield(:,:,:,:),newfield(:,:,:,:)
        character file*256,datfile*256,string*40,lgt*1,opera*3
        character vars(nvmax)*20,lvars(nvmax)*255,title*200,
     +       history*10000,units(nvmax)*20,cell_methods(nvmax)*200,
     +       lz(3)*20,ltime*100,svars(nvmax)*50,punits*40
        logical normal,tdefined(ntmax),lexist
        integer iargc,llen
        lwrite = .false.
        itype = 0               ! for the time being - no vector fields
        ivars(1,:) = 0
        ivars(2,:) = 99
        lgt = ' '
*
        if ( iargc().lt.4 ) then
            print *,'usage: daily2longerfield infield nperyearnew '//
     +           'mean|sd|var|sum|abo|bel|min|max|num|mintime|maxtime'//
     +           '|firsttime|lasttime|con [<> val[%|p]] outfield.ctl'
            print *,'(more options will come as requested)'
            stop
        end if
*
*       read data
*
        call keepalive(0,0)
        call getarg(1,file)
        if ( lwrite ) print *,'daily2longerfield: nf_opening file '
     +       ,trim(file)
        status = nf_open(file,nf_nowrite,ncid)
        if ( status.ne.nf_noerr ) then
            call parsectl(file,datfile,nxmax,nx,xx,nymax,ny,yy,nzmax,nz
     +            ,zz,nt,nperyear,firstyr,firstmo,undef,endian,title,1
     +            ,nvars,vars,ivars,lvars,units)
            lz = ' '
            svars = ' '
            ltime = 'time'
            history = ' '
            ncid = -1
            cell_methods = ' '
        else
            if ( lwrite ) print *,'calling parsenc on ',trim(file)
            call ensparsenc(file,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax
     +           ,nz,zz,lz,nt,nperyear,firstyr,firstmo,ltime,tdefined
     +           ,ntmax,nens1,nens2,undef,title,history,1,nvars,vars
     +           ,jvars,lvars,svars,units,cell_methods)
            ntp = 0
            do it=1,nt
                if ( tdefined(it) ) ntp = ntp + 1
            end do
        endif
        if ( nz.ne.1 ) then
            write(0,*)
     +           'daily2longefield: error: cannot handle 3D fields yet'
            write(*,*)
     +           'daily2longefield: error: cannot handle 3D fields yet'
            call abort
        endif
        if ( nperyear.ne.366 ) then
            lastyr = firstyr + (nt+firstmo-2)/nperyear
        else
            lastyr = firstyr + int((nt+firstmo-2)/365.25)
        end if
        if ( lwrite ) then
            print *,'daily2longerfield: allocating oldfield'
            print *,'nx,ny,nperyear,firstyr,lastyr = ',nx,ny,nperyear
     +           ,firstyr,lastyr
            print *,'size = ',4*nx*ny*nperyear*(lastyr-firstyr+1)
        endif
        allocate(oldfield(nx,ny,nperyear,firstyr:lastyr))
        if ( nperyear.le.366 ) then
            allocate(cut(nx,ny,nperyear))
        else
            allocate(cut(1,1,nperyear))
        endif
        call keepalive(0,0)
        if ( ncid.eq.-1 ) then
            call zreaddatfile(datfile,oldfield,nx,ny,nz,nx,ny,nz
     +           ,nperyear,firstyr,lastyr,firstyr,firstmo,nt,undef
     +           ,endian,lwrite,firstyr,lastyr,1,nvars)
        else
            call readncfile(ncid,oldfield,nx,ny,nx,ny,nperyear
     +           ,firstyr,lastyr,firstyr,firstmo,ntp,undef,lwrite
     +           ,firstyr,lastyr,jvars)
            if ( nt.ne.ntp ) then
                call fixholefield(oldfield,nx,ny,1,nperyear,firstyr
     +               ,lastyr,firstyr,firstmo,ntp,nt,tdefined,lwrite)
            endif
        end if
        call keepalive(0,0)
*
*       read operation
*
        call getarg(2,string)
        read(string,*,err=901) nperyearnew
        call getarg(3,string)
        if ( string.eq.'maxtime' ) string = 'xti'
        if ( string.eq.'mintime' ) string = 'nti'
        if ( string.eq.'firsttime' ) string = 'fti'
        if ( string.eq.'lasttime' ) string = 'lti'
        opera = string
        if ( opera.ne.'mea' .and. opera.ne.'sd ' .and.
     +       opera.ne.'var' .and.
     +       opera.ne.'min' .and. opera.ne.'max' .and.
     +       opera.ne.'nti' .and. opera.ne.'xti' .and.
     +       opera.ne.'fti' .and. opera.ne.'lti' .and.
     +       opera.ne.'num' .and. opera.ne.'sum' .and.
     +       opera.ne.'bel' .and. opera.ne.'abo' .and.
     +       opera.ne.'con' ) then
            write(0,*) 'daily2longerfield: error: unknown operation '
     +           ,opera
            call abort
        end if
        normal = .false.
        if ( iargc().ge.5 ) then
            call getarg(4,string)
            if ( string.eq.'lt' ) lgt = '<'
            if ( string.eq.'gt' ) lgt = '>'
            if ( string.eq.'<' ) lgt = '<'
            if ( string.eq.'>' ) lgt = '>'
            if ( lgt.eq.'<' .or. lgt.eq.'>' ) then
                call getarg(5,string)
                if ( string.eq.'n' ) then
*                   take normals wrt 1971-2000
                    allocate(normfield(nx,ny,nperyear))
                    allocate(nn(nx,ny))
                    do j=1,nperyear
                        normfield(1:nx,1:ny,j) = 0
                        nn(1:nx,1:ny) = 0
                        do yr=1971,2000
                            do jy=1,ny
                                do jx=1,nx
                                    if ( oldfield(jx,jy,j,yr).lt.1e33 )
     +                                   then
                                        nn(jx,jy) = nn(jx,jy) + 1
                                        normfield(jx,jy,j) = 
     +                                       normfield(jx,jy,j)
     +                                       + oldfield(jx,jy,j,yr)
                                    end if
                                end do
                            end do
                        end do
                        do jy=1,ny
                            do jx=1,nx
                                if ( nn(jx,jy).gt.5 ) then ! arbitrary number
                                    normfield(jx,jy,j) =
     +                                   normfield(jx,jy,j)/nn(jx,jy)
                                else
                                    normfield(jx,jy,j) = 3e33
                                end if
                            end do
                        end do
                    end do
*                   no smoothing for the time being
                    do yr=yrbeg,yrend
                        do j=1,nperyear
                            do jy=1,ny
                                do jx=1,nx
                                    if ( oldfield(jx,jy,j,yr).lt.1e33
     +                                   .and.normfield(jx,jy,j).lt.1e33
     +                                   ) then
                                        oldfield(jx,jy,j,yr) = 
     +                                       oldfield(jx,jy,j,yr)
     +                                       - normfield(jx,jy,j)
                                    else
                                        oldfield(jx,jy,j,yr) = 3e33
                                    end if
                                end do
                            end do
                        end do
                    enddo
                    deallocate(normfield)
                    deallocate(nn)
                    lgt = ' '
                else            ! string != 'n'
                    if ( lgt.ne.' ' ) then
                        if ( nperyear.gt.366 ) then
                            write(*,*) 'daily2longerfield: cuts not yet'
     +                           //' implemented for sub-daily data'
                            write(0,*) 'daily2longerfield: cuts not yet'
     +                           //' implemented for sub-daily data'
                            call abort
                        endif
                        i = index(string,'%')
                        if ( i.eq.0 ) i = index(string,'p')
                        if ( i.eq.0 ) i = len(string) + 1
                        read(string(:i-1),*,err=902) pcut
                        if ( i.le.len(string) ) then
                            punits = '%'
                            write(0,*) 'Computing percentage cuts will '
     +                           ,'take a lot of time...<p>'
                            allocate(oldseries(nperyear,firstyr:lastyr))
                            do jy=1,ny
                                do jx=1,nx
                                    do yr=firstyr,lastyr
                                        do j=1,nperyear
                                            oldseries(j,yr) = 
     +                                           oldfield(jx,jy,j,yr)
                                        end do
                                    end do
                                    do j=1,nperyear
                                        call getcutoff(cut(jx,jy,j),pcut
     +                                       ,oldseries,nperyear
     +                                       ,nperyear,firstyr,lastyr
     +                                       ,firstyr,lastyr,j,j,0)
                                    end do
                                end do
                                call keepalive(jy,ny)
                            end do
                            deallocate(oldseries)
                        else
                            cut = pcut
                            punits = units(1)
                        endif
                    end if
                end if
                i = 6
            else
                i = 4
            end if
            call getopts(i,iargc()-1,nperyear,yrbeg,yrend,.true.,0,0)
        end if
        if ( minfac.lt.0 ) minfac = 0.5
*      
*       no climatology -- too expensive, only very limited usefulness.
*
*
*       perform operation
*
        if ( lwrite ) print *,'allocating newfield'
        allocate(newfield(nx,ny,abs(nperyearnew),firstyr:lastyr))
        if ( lwrite ) print *,'calling fieldallday2period'
        call fieldallday2period(
     +       oldfield,nperyear,
     +       newfield,nperyearnew,
     +       nx,ny,firstyr,lastyr,opera,lgt,cut,minfac,
     +       itype,vars(1),units(1),lwrite)
        deallocate(oldfield)
        nperyearnew = abs(nperyearnew)
*
*       standard units
*
        if ( lstandardunits ) then
            call makestandardfield(newfield,nx,ny,1,abs(nperyearnew),
     +           firstyr,lastyr,nx,ny,1,nperyearnew,firstyr,lastyr,
     +           vars(1),units(1),lwrite)
        endif
*
*       adjust names
*
        call adjustnames(opera,nperyear,nperyearnew,lgt,pcut,punits
     +       ,lvars(1),cell_methods(1))
        call adjustvar(opera,vars(1),lwrite)
*
*       output field
*
        call getarg(iargc(),file)
        i = index(file,'.ctl')
        if ( i.ne.0 ) then
!           grads file output (deprecated)
            datfile = file
            datfile(i:) = '.grd'
            undef = 3e33
            inquire(file=trim(file),exist=lexist)
            if ( lexist ) then
                call rsunit(unit)
                open(unit,file=trim(file))
                close(unit,status='delete')
            end if
            inquire(file=trim(datfile),exist=lexist)
            if ( lexist ) then
                call rsunit(unit)
                open(unit,file=trim(datfile))
                close(unit,status='delete')
            end if
            call writectl(file,datfile,nx,xx,ny,yy,nz,zz,
     +           nperyearnew*(lastyr-firstyr+1),nperyearnew,
     +           firstyr,1,undef,title,nvars,vars,ivars,lvars,units)
            open(1,file=trim(datfile),access='direct',form='unformatted'
     +           ,recl=recfa4*nx*ny,status='new')
            i=0
            do yr=firstyr,lastyr
                do j=1,nperyearnew
                    i = i + 1
                    write(1,rec=i) ((newfield(jx,jy,j,yr),jx=1,nx),jy=1
     +                   ,ny)
                end do
            end do
        else
!           netcdf output
            undef = 3e33
            call enswritenc(file,ncid,ntvarid,itimeaxis,ndata,nx,xx,ny
     +           ,yy,nz,zz,lz,nperyearnew*(lastyr-firstyr+1),nperyearnew
     +           ,firstyr,1,ltime,undef,title,history,nvars,vars,ivars
     +           ,lvars,svars,units,cell_methods,nens1,nens2)
            i = 0
            do yr=firstyr,lastyr
                do j=1,nperyearnew
                    i = i + 1
                    call writencslice(ncid,0,0,0,ivars,
     +                   newfield(1,1,j,yr),nx,ny,nz,nx,ny,nz,i,1)
                end do
            end do
          endif
*
*       error messages
*
        goto 999
  901   write(0,*) 'daily2longer: expecting nperyearnew, not ',string
        call abort
  902   write(0,*) 'daily2longer: expecting value[%|p], not ',string
        call abort
  999   continue
        end
