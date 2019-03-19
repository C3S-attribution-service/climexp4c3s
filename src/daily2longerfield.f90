program daily2longerfield

!   compute aggregate quantities from daily data
!   input: daily field series
!   output: yearly/monthly/10-dy/5-dy field

    implicit none
    include 'params.h'
    include 'netcdf.inc'
    include 'getopts.inc'
    integer,parameter :: ntmax=100000,recfa4=4
    integer :: nperyear,nperyearnew,yr,mo,dy,i,j,n,itype,jx,jy
    integer :: nx,ny,nz,nt,firstyr,firstmo,lastyr,nvars, &
        ivars(2,nvmax),jvars(6,nvmax),ncid,endian,status, &
        ntvarid,itimeaxis(ndata),it,ntp,unit,ncount
    integer,allocatable :: nn(:,:)
    real :: pcut,s
    real,allocatable :: cut(:,:,:),normfield(:,:,:),oldseries(:,:)
    real :: xx(nxmax),yy(nymax),zz(nzmax),undef
    real,allocatable :: oldfield(:,:,:,:),newfield(:,:,:,:),fxy(:,:),refs(:)
    character file*256,datfile*256,string*40,lgt*1,opera*3
    character vars(nvmax)*20,lvars(nvmax)*255,title*200, &
        history*20000,units(nvmax)*20,cell_methods(nvmax)*200, &
        lz(3)*20,ltime*100,svars(nvmax)*50,punits*40,metadata(2,100)*2000
    logical :: normal,tdefined(ntmax),lexist
    logical,allocatable :: lvalid(:,:,:,:),lvalid1(:,:)
    integer :: command_argument_count
    lwrite = .false. 
    itype = 0               ! for the time being - no vector fields
    ivars(1,:) = 0
    ivars(2,:) = 99
    lgt = ' '

    if ( command_argument_count() < 4 ) then
        print *,'usage: daily2longerfield infield nperyearnew '// &
            'mean|sd|var|sum|abo|bel|min|max|num|mintime|maxtime'// &
            '|firsttime|lasttime|con [<> val[%|p]] [ave|sum N] outfield.ctl'
        print *,'(more options will come as requested)'
        stop
    end if

!   read data

    call keepalive1('Reading field',0,0)
    call get_command_argument(1,file)
    if ( lwrite ) print *,'daily2longerfield: nf_opening file ',trim(file)
    status = nf_open(file,nf_nowrite,ncid)
    if ( status /= nf_noerr ) then
        call parsectl(file,datfile,nxmax,nx,xx,nymax,ny,yy,nzmax,nz &
            ,zz,nt,nperyear,firstyr,firstmo,undef,endian,title,1 &
            ,nvars,vars,ivars,lvars,units)
        lz = ' '
        svars = ' '
        ltime = 'time'
        history = ' '
        ncid = -1
        cell_methods = ' '
        metadata = ' '
    else
        if ( lwrite ) print *,'calling parsenc on ',trim(file)
        call ensparsenc(file,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax &
            ,nz,zz,lz,nt,nperyear,firstyr,firstmo,ltime,tdefined &
            ,ntmax,nens1,nens2,undef,title,history,1,nvars,vars &
            ,jvars,lvars,svars,units,cell_methods,metadata)
        ntp = 0
        do it=1,nt
            if ( tdefined(it) ) ntp = ntp + 1
        end do
    endif
    if ( nz /= 1 ) then
        write(0,*) 'daily2longefield: error: cannot handle 3D fields yet'
        write(*,*) 'daily2longefield: error: cannot handle 3D fields yet'
        call exit(-1)
    endif
    call getlastyr(firstyr,firstmo,nt,nperyear,lastyr)
    if ( lwrite ) then
        print *,'daily2longerfield: allocating oldfield'
        print *,'nx,ny,nperyear,firstyr,lastyr = ',nx,ny,nperyear,firstyr,lastyr
        print *,'size = ',4*nx*ny*nperyear*(lastyr-firstyr+1)
    endif
    allocate(oldfield(nx,ny,nperyear,firstyr:lastyr))
    allocate(lvalid(nx,ny,nperyear,firstyr:lastyr))
    if ( nperyear <= 366 ) then
        allocate(cut(nx,ny,nperyear))
    else
        allocate(cut(1,1,nperyear))
    endif
    call keepalive1('Reading field',0,0)
    if ( ncid == -1 ) then
        call zreaddatfile(datfile,oldfield,nx,ny,nz,nx,ny,nz &
            ,nperyear,firstyr,lastyr,firstyr,firstmo,nt,undef &
            ,endian,lwrite,firstyr,lastyr,1,nvars)
    else
        call readncfile(ncid,oldfield,nx,ny,nx,ny,nperyear &
            ,firstyr,lastyr,firstyr,firstmo,ntp,undef,lwrite &
            ,firstyr,lastyr,jvars)
        if ( nt /= ntp ) then
            call fixholefield(oldfield,nx,ny,1,nperyear,firstyr &
            ,lastyr,firstyr,firstmo,ntp,nt,tdefined,lwrite)
        endif
    end if
    call keepalive1('Starting computation',0,0)

!   read operation

    call get_command_argument(2,string)
    read(string,*,err=901) nperyearnew
    call get_command_argument(3,string)
    if ( string == 'maxtime' ) string = 'xti'
    if ( string == 'mintime' ) string = 'nti'
    if ( string == 'firsttime' ) string = 'fti'
    if ( string == 'lasttime' ) string = 'lti'
    opera = string
    if ( opera /= 'mea' .and. opera /= 'sd ' .and. &
         opera /= 'var' .and. &
         opera /= 'min' .and. opera /= 'max' .and. &
         opera /= 'nti' .and. opera /= 'xti' .and. &
         opera /= 'fti' .and. opera /= 'lti' .and. &
         opera /= 'num' .and. opera /= 'sum' .and. &
         opera /= 'bel' .and. opera /= 'abo' .and. &
         opera /= 'con' ) then
        write(0,*) 'daily2longerfield: error: unknown operation ',opera
        call exit(-1)
    end if
    normal = .false. 
    if ( command_argument_count() >= 5 ) then
        call get_command_argument(4,string)
        if ( string == 'lt' ) lgt = '<'
        if ( string == 'gt' ) lgt = '>'
        if ( string == '<' ) lgt = '<'
        if ( string == '>' ) lgt = '>'
        if ( lgt == '<' .or. lgt == '>' ) then
            call get_command_argument(5,string)
            if ( string == 'n' ) then
!               take normals wrt 1971-2000
                allocate(normfield(nx,ny,nperyear))
                allocate(nn(nx,ny))
                do j=1,nperyear
                    normfield(1:nx,1:ny,j) = 0
                    nn(1:nx,1:ny) = 0
                    do yr=1971,2000
                        do jy=1,ny
                            do jx=1,nx
                                if ( oldfield(jx,jy,j,yr) < 1e33 ) &
                                then
                                    nn(jx,jy) = nn(jx,jy) + 1
                                    normfield(jx,jy,j) = normfield(jx,jy,j) + oldfield(jx,jy,j,yr)
                                end if
                            end do
                        end do
                    end do
                    do jy=1,ny
                        do jx=1,nx
                            if ( nn(jx,jy) > 5 ) then ! arbitrary number
                                normfield(jx,jy,j) = normfield(jx,jy,j)/nn(jx,jy)
                            else
                                normfield(jx,jy,j) = 3e33
                            end if
                        end do
                    end do
                end do
!               no smoothing for the time being
                do yr=yrbeg,yrend
                    do j=1,nperyear
                        do jy=1,ny
                            do jx=1,nx
                                if ( oldfield(jx,jy,j,yr) < 1e33 &
                                 .and. normfield(jx,jy,j) < 1e33 ) then
                                    oldfield(jx,jy,j,yr) = oldfield(jx,jy,j,yr) - normfield(jx,jy,j)
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
                if ( lgt /= ' ' ) then
                    if ( nperyear > 366 ) then
                        write(*,*) 'daily2longerfield: cuts not yet implemented for sub-daily data'
                        write(0,*) 'daily2longerfield: cuts not yet implemented for sub-daily data'
                        call exit(-1)
                    endif
                    i = index(string,'%')
                    if ( i == 0 ) i = index(string,'p')
                    if ( i == 0 ) i = len(string) + 1
                    read(string(:i-1),*,err=902) pcut
                    if ( i <= len(string) ) then
                        punits = '%'
                        write(0,*) 'Computing percentage cuts will take a lot of time...<p>'
                        allocate(oldseries(nperyear,firstyr:lastyr))
                        do jy=1,ny
                            do jx=1,nx
                                do yr=firstyr,lastyr
                                    do j=1,nperyear
                                        oldseries(j,yr) = oldfield(jx,jy,j,yr)
                                    end do
                                end do
                                if ( lastyr > firstyr ) then
                                    ! compute percentiles over time for each day/month of the year
                                    do j=1,nperyear
                                        call getcutoff(cut(jx,jy,j),pcut,oldseries,nperyear &
                                            ,nperyear,firstyr,lastyr,firstyr,lastyr,j,j,0)
                                    end do
                                else
                                    ! compute percentiles over the seasonal cycle
                                    call getcutoff(cut(jx,jy,1),pcut,oldseries,nperyear &
                                        ,nperyear,firstyr,lastyr,firstyr,lastyr,1,nperyear,0)
                                    cut(jx,jy,2:nperyear) = cut(jx,jy,1)
                                end if
                            end do
                            call keepalive1('Normals for latitude',jy,ny)
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
        call getopts(i,command_argument_count()-1,nperyear,yrbeg,yrend, .true. ,0,0)
    end if
    if ( opera == 'fti' .or. opera == 'lti' ) then
        if ( lgt == ' ' ) then
            write(0,*) 'daily2longerfield: error: need a threshold for first/last time'
            call exit(-1)
        end if
    end if
    if ( minfac < 0 ) minfac = 0.5
    yr1 = max(yr1,firstyr)
    yr2 = min(yr2,lastyr)
!   no climatology -- too expensive, only very limited usefulness.
!   fill in missing data when requested
!   this test should be exactly the same as in fieldday2period
    if ( ( opera == 'mea' .or. opera == 'sum' ) &
            .and. lgt == ' ' .and. add_option > 0 ) then
        allocate(fxy(nperyear,firstyr:lastyr))
        allocate(refs(firstyr:lastyr))
        allocate(lvalid1(nperyear,firstyr:lastyr))
        do j=1,ny
            call keepalive1('Handling missing data',j,ny)
            do i=1,nx
                do mo=1,nperyear
                    do yr=yr1,yr2
                        fxy(mo,yr) = oldfield(i,j,mo,yr)
                    end do
                end do
                call fillmissingdata(fxy,lvalid1,refs,nperyear, &
                    firstyr,lastyr,nperyear,add_option, .false. ,lwrite)
                do mo=1,nperyear
                    do yr=yr1,yr2
                        lvalid(i,j,mo,yr) = lvalid1(mo,yr)
                        if ( oldfield(i,j,mo,yr) > 1e33 .and. fxy(mo,yr) < 1e33 ) then
                            oldfield(i,j,mo,yr) = fxy(mo,yr)
                        end if
                    end do
                end do
            end do
        end do
        deallocate(fxy)
        deallocate(refs)
        deallocate(lvalid1)
    end if

!   take running mean if requested

    if ( lsum > 1 ) then
        allocate(fxy(nperyear,firstyr:lastyr))
        ncount = 0
!$omp parallel do private(j,i,mo,yr,fxy)
        do j=1,ny
!$omp atomic update
            ncount = ncount + 1
            call keepalive1('Summing latitude',ncount,ny)
            do i=1,nx
                do mo=1,nperyear
                    do yr=yr1,yr2
                        fxy(mo,yr) = oldfield(i,j,mo,yr)
                    end do
                end do
                call sumit(fxy,nperyear,nperyear,firstyr,lastyr,lsum,'v')
                do mo=1,nperyear
                    do yr=yr1,yr2
                        oldfield(i,j,mo,yr) = fxy(mo,yr)
                    end do
                end do
            end do
        end do
!$omp end parallel do
        deallocate(fxy)
    end if

!   perform operation

    if ( lwrite ) print *,'allocating newfield'
    allocate(newfield(nx,ny,abs(nperyearnew),firstyr:lastyr))
    if ( lwrite ) print *,'calling fieldallday2period'
    call fieldallday2period(oldfield,nperyear,lvalid,newfield,nperyearnew, &
        nx,ny,firstyr,lastyr,opera,lgt,cut,minfac,add_option, &
        itype,vars(1),units(1),lwrite)
    deallocate(oldfield)
    nperyearnew = abs(nperyearnew)

!   standard units

    if ( lstandardunits ) then
        call makestandardfield(newfield,nx,ny,1,abs(nperyearnew), &
            firstyr,lastyr,nx,ny,1,nperyearnew,firstyr,lastyr, &
            vars(1),units(1),lwrite)
    endif

!   adjust names

    call adjustnames(opera,nperyear,nperyearnew,lgt,pcut,punits &
        ,lvars(1),cell_methods(1))
    call adjustvar(opera,vars(1),lwrite)

!   output field

    call get_command_argument(command_argument_count(),file)
    i = index(file,'.ctl')
    if ( i /= 0 ) then
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
        call writectl(file,datfile,nx,xx,ny,yy,nz,zz, &
        nperyearnew*(lastyr-firstyr+1),nperyearnew, &
        firstyr,1,undef,title,nvars,vars,ivars,lvars,units)
        open(1,file=trim(datfile),access='direct',form='unformatted' &
            ,recl=recfa4*nx*ny,status='new')
        i=0
        do yr=firstyr,lastyr
            do j=1,nperyearnew
                i = i + 1
                write(1,rec=i) ((newfield(jx,jy,j,yr),jx=1,nx),jy=1 &
                ,ny)
            end do
        end do
    else
    !           netcdf output
        undef = 3e33
        call enswritenc(file,ncid,ntvarid,itimeaxis,ndata,nx,xx,ny &
            ,yy,nz,zz,lz,nperyearnew*(lastyr-firstyr+1),nperyearnew &
            ,firstyr,1,ltime,undef,title,history,nvars,vars,ivars &
            ,lvars,svars,units,cell_methods,metadata,nens1,nens2)
        i = 0
        do yr=firstyr,lastyr
            do j=1,nperyearnew
                i = i + 1
                call writencslice(ncid,0,0,0,ivars, &
                newfield(1,1,j,yr),nx,ny,nz,nx,ny,nz,i,1)
            end do
        end do
        status = nf_close(ncid)
    endif

!       error messages

    goto 999
    901 write(0,*) 'daily2longer: expecting nperyearnew, not ',string
    call exit(-1)
    902 write(0,*) 'daily2longer: expecting value[%|p], not ',string
    call exit(-1)
    999 continue
end program daily2longerfield
