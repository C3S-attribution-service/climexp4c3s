program difffield
!
!   compute difference map between two fields
!
    use lsdata
    implicit none
    include 'params.h'
    include 'recfac.inc'
    include 'getopts.inc'
    include 'netcdf.inc'
    integer ifield,nperyear,nperyear1,iens,mens1,mens,i,j,k,l,mo,month,n,jx,jy     &
     &       ,jz,nx1,ny1,nz1,nx2,ny2,nz2,nxf,nyf,nzf,j1,j2,startopts,fyr,lyr
    integer nx,ny,nz,nt,firstyr,firstmo,lastyr,nvars,  &
     &       ivars(2,3),jvars(6,3),ncid,endian,status,pow, &
     &       ntvarid,itimeaxis(npermax),nmetadata1
    integer,allocatable :: nn(:,:,:,:)
    real pcut,s,ss(npermax,2),d,e,z,p
    real xx(nxmax),yy(nymax),zz(nzmax),undef
    real xxls(nxmax),yyls(nymax),zzls(nzmax)
    real xx1(nxmax),yy1(nymax),zz1(nzmax),xx2(nxmax),yy2(nymax),zz2(nzmax)
    real,allocatable :: field(:,:,:,:,:,:),ave1(:,:,:,:,:),ave2(:,:,:,:,:)
    character file*1023,file1*1023,outfile*1023,datfile*1023,string*40
    character vars(3)*20,lvars(3)*120,svars(3)*120,title*4096,           &
             history*50000,lz(3)*20,units(3)*20,title1*5100,title2*500,   &
             cell_methods(3)*100,ltime*20,lsmasktype*4,tmpunits*20,      &
             yesno*3,metadata(2,100)*2000,metadata1(2,100)*2000,         &
             history1*50000
    logical lexist,lequal
    integer get_endian
    real erfcc
    lstandardunits = .true.
    lwrite = .false.
    do i=1,command_argument_count()
        call get_command_argument(i,string)
        if ( string == 'debug' .or. string == 'lwrite' ) then
            lwrite = .true.
        end if
    end do

    if ( command_argument_count() < 3 ) then
        print *,'usage: difffield field1 field2 [lsmask landseamask all|land|sea] '//     &
                 '[mon n] [ave|sum n] [begin yr1] [end yr2] [begin2 yr1a] [end2 yr2a] '// &
                 '[normsd] outfield.nc'
        call exit(-1)
    end if
!
!   read landseamask
!
    startopts = 3
    call getlsmask(startopts,lsmasktype,nxmax,xxls,nymax,yyls,lwrite)
!
!   read data
!
    do ifield=1,2
        call get_command_argument(ifield,file)
        call getmetadata(file,mens1,mens,ncid,datfile,nxmax,nx           &
                 ,xx,nymax,ny,yy,nzmax,nz,zz,lz,nt,nperyear,firstyr      &
                 ,firstmo,ltime,undef,endian,title,history,1,nvars,vars  &
                 ,jvars,lvars,svars,units,cell_methods,metadata,lwrite)
        
        if ( ifield == 1 ) then
            call add_varnames_metadata(vars(1),lvars(1),svars(1),metadata,'variable1')
            nperyear1 = nperyear
            metadata1 = metadata
            file1 = file
            history1 = history
        else
            call add_varnames_metadata(vars(1),lvars(1),svars(1),metadata,'variable2')
            if ( nperyear /= nperyear1 ) then
                write(0,*) 'difffield: error: nperyear unequal ',nperyear1,nperyear
                call exit(-1)
            end if
        end if
        call getopts(startopts,command_argument_count()-1,nperyear,yrbeg,yrend,.false.,mens1,mens)
        if ( m1 == 0 ) then
            if ( m2 == 0 ) then
                m1 = 1
                m2 = nperyear
            else
                m1 = 1
            endif
        endif
        lastyr = firstyr + (nt+firstmo-2)/nperyear
        if ( ifield == 1 ) then
            yr1 = max(yr1,firstyr)
            yr2 = min(yr2,lastyr)
        else
            yr1 = max(yr1a,firstyr)
            yr2 = min(yr2a,lastyr)
        endif
        if ( lwrite ) then
            print *,'firstyr,lastyr,yr1,yr2 = ',firstyr,lastyr,yr1,yr2
            print *,'difffield: allocating field'
            print *,'nx,ny,nz,nperyear,yr1,yr2 = ',nx,ny,nz,nperyear,yr1,yr2,nens1,nens2
            print *,'size = ',4*nx*ny*nz*nperyear*(yr2-yr1+1)*(nens2-nens1+1)
        endif
        fyr = yr1
        lyr = yr2
        allocate(field(nx,ny,nz,nperyear,fyr:lyr,nens1:nens2))
        call readfield(ncid,file,datfile,field,nx,ny                       &
                 ,nz,nperyear,fyr,lyr,nens1,nens2,nx,ny,nz,nperyear,yr1    &
                 ,yr2,firstyr,firstmo,nt,undef,endian,vars,units           &
                 ,lstandardunits,lwrite)
!
!       apply land/sea mask, before the axes are transformed
!
        if ( lwrite ) print *,'selecting ',lsmasktype,' points of land/sea mask'
        if ( lsmasktype /= 'all' ) then
            call aregridsequal(nxls,nyls,xxls,yyls,nx,ny,xx,yy,lequal,lwrite)
            if ( lwrite ) print *,'field ',ifield,' grids equal? ',lequal
            if ( lequal ) then
                call applylsmask(field,lsmask,nx,ny,1,nperyear             &
                         ,yr1,yr2,nens1,nens2,lsmasktype,lwrite)
            endif
        endif
!
!       compute average
!
        if ( ifield == 1 ) then
            allocate(ave1(nx,ny,nz,nperyear,2))
            ave1 = 3e33
            allocate(nn(nx,ny,nz,nperyear))
            nn = 0
            nxf = nx
            nxf = ny
            nzf = nz
            nx1 = nx
            ny1 = ny
            nz1 = nz
            xx1(1:nx) = xx(1:nx)
            yy1(1:ny) = yy(1:ny)
            zz1(1:nz) = zz(1:nz)
            nt = (yr2-yr1+1)*nperyear
            write(title1,'(2a,i4,a,i4)') trim(title),' ',yr1,':',yr2
            call getenswinmean(ave1,nn,nx,ny,nz,nperyear,field,nx,ny        &
                     ,nz,nperyear,fyr,lyr,nens1,nens2,nx,ny,nz,yr1,1,nt     &
                     ,1,nx,1,ny,lwrite)
!           convert variance to error estimate
            do month=1,nperyear
                do jz=1,nz
                    do jy=1,ny
                        do jx=1,nx
                            if ( ave1(jx,jy,jz,month,2) < 1e33 .and.   &
                                 nn(jx,jy,jz,month) > 0 ) then
                                ave1(jx,jy,jz,month,2) = ave1(jx,jy,jz,month,2)/nn(jx,jy,jz,month)
                            else
                                ave1(jx,jy,jz,month,2) = 3e33
                            end if
                        end do
                    end do
                end do
            end do
            if ( lwrite ) print *,'ave1(',nx/2,ny/2,(nz+1)/2,1,        &
                     ') = ',ave1(nx/1,ny/2,(nz+1)/2,1,1),              &
                     ave1(nx/1,ny/2,(nz+1)/2,1,2)
        else if ( ifield == 2 ) then
            nxf = max(nx1,nx)
            nyf = max(ny1,ny)
            nzf = max(nz1,nz)
            allocate(ave2(nxf,nyf,nzf,nperyear,2))
            if ( nxf /= nx1 .or. nyf /= ny1 .or. nzf /= nz1 ) then
!                   make sure ave1 has the same size as ave2
                ave2(1:nx1,1:ny1,1:nz1,:,:) = ave1(1:nx1,1:ny1,1:nz1,:,:)
                if ( lwrite )print *,'resizing ave1 to ',nxf,nyf,nzf
                deallocate(ave1)
                deallocate(nn)
                allocate(ave1(nxf,nyf,nzf,nperyear,2))
                ave1 = 3e33
                allocate(nn(nxf,nyf,nzf,nperyear))
                nn = 0
                ave1(1:nx1,1:ny1,1:nz1,:,:) = ave2(1:nx1,1:ny1,1:nz1,:,:)
            end if
            nx2 = nx
            ny2 = ny
            nz2 = nz
            xx2(1:nx) = xx(1:nx)
            yy2(1:ny) = yy(1:ny)
            zz2(1:nz) = zz(1:nz)
            nt = (yr2-yr1+1)*nperyear
            call getenswinmean(ave2,nn,nxf,nyf,nzf,nperyear,field,nx       &
                     ,ny,nz,nperyear,fyr,lyr,nens1,nens2,nx,ny,nz,yr1,1    &
                     ,nt,1,nx,1,ny,lwrite)
!           convert variance to error estimate
            do month=1,nperyear
                do jz=1,nz
                    do jy=1,ny
                        do jx=1,nx
                            if ( ave2(jx,jy,jz,month,2) < 3e33 .and. nn(jx,jy,jz,month) > 0 ) then
                                ave2(jx,jy,jz,month,2) = ave2(jx,jy,jz,month,2)/nn(jx,jy,jz,month)
                            else
                                ave2(jx,jy,jz,month,2) = 3e33
                            end if
                        end do
                    end do
                end do
            end do
            if ( lwrite ) print *,'ave2(',nx/2,ny/2,(nz+1)/2,1,    &
                     ') = ',ave2(nx/1,ny/2,(nz+1)/2,1,1),          &
                     ave2(nx/1,ny/2,(nz+1)/2,1,2)
        else
            write(0,*) '???'
            call exit(-1)
        end if
        deallocate(field)
!
!       average seasons
!
        if ( ifield == 1 .and. lsum > 1 .or.     &
             ifield == 2 .and. lsum2 > 1 ) then
            do jz=1,nz
                do jy=1,ny
                    do jx=1,nx
                        ss = 0
                        do month=m1,m2
                            n = 0
                            if ( ifield == 1 ) then
                                l = lsum
                            else
                                l = lsum2
                            endif
                            do pow=1,2
                                do k=1,l
                                    mo = month + k - 1
                                    mo = 1 + mod(mo-1,12)
                                    call getj1j2(j1,j2,mo,nperyear,.false.)
                                    do j=j1,j2
                                        if ( ifield == 1 ) then
                                            if ( ave1(jx,jy,jz,j,pow) < 1e30 ) then
                                                ss(month,pow) = ss(month,pow)+ave1(jx,jy,jz,j,pow)
                                                n = n + 1
                                            endif
                                        else
                                            if ( ave2(jx,jy,jz,j,pow) < 1e30 ) then
                                                ss(month,pow) = ss(month,pow)+ave2(jx,jy,jz,j,pow)
                                                n = n + 1
                                            endif
                                        endif
                                    enddo
                                enddo
                                if ( n > 0 ) then
                                    ss(month,pow) = ss(month,pow)/n
                                else
                                    ss(month,pow) = 3e33
                                endif
                            end do
                        end do
                        if ( ifield == 1 ) then
                            do month=m1,m2
                                ave1(jx,jy,jz,month,:) = ss(month,:)
                            end do
                        else
                            do month=m1,m2
                                ave2(jx,jy,jz,month,:) = ss(month,:)
                            end do
                        end if
                    end do
                end do
            end do
        end if
    end do
!
!   interpolate to a common grid
!
    call xyinterpu(ave1,xx1,nx1,yy1,ny1,ave2,xx2,nx2,yy2,ny2,             &
             xx,nx,yy,ny,1,2,1,2,nxf,nyf,nzf,nz,nperyear,intertype,lwrite)
    if ( lwrite ) then
        print *,'ave1(',nx/2,ny/2,(nz+1)/2,m1,') = ',ave1(nx/1,ny/2,(nz+1)/2,m1,1),  &
                     ave1(nx/1,ny/2,(nz+1)/2,m1,2)
        print *,'ave2(',nx/2,ny/2,(nz+1)/2,m1,') = ',ave2(nx/1,ny/2,(nz+1)/2,m1,1),  &
                     ave2(nx/1,ny/2,(nz+1)/2,m1,2)
    end if
!
!   output field
!
    call get_command_argument(command_argument_count(),outfile)
    datfile = outfile
    i = index(outfile,'.ctl')
    if ( i /= 0 ) then
        datfile(i:) = '.dat'
    endif
    title2 = title
    if ( lwrite ) then
        print *,'title1 = ',trim(title1)
        print *,'title2 = ',trim(title2)
    endif
    write(title,'(4a,i4,a,i4)') trim(title1),' minus ',trim(title2),' ',yr1,':',yr2
!   note that we used yr1,yr2 for the second field...
    if ( lwrite ) then
        print *,'title  = ',trim(title)
    endif
    if ( history == history1 ) history = ' '
    call merge_metadata(metadata1,nmetadata1,metadata,title2,history,'field2_')
    undef = 3e33
    nvars = 3
    if ( lnormsd ) then
        vars(1) = 'reldiff'
        lvars(1) = trim(lvars(1))//' relative difference'
        units(1) = '1'
    else
        vars(1) = 'diff'
        lvars(1) = trim(lvars(1))//' difference'
    endif
    ivars(1,1) = 1
    ivars(2,1) = 99
    vars(2) = 'error'
    lvars(2) = 'error on '//trim(lvars(1))
    units(2) = units(1)
    cell_methods(2) = cell_methods(1)
    ivars(1:2,2) = ivars(1:2,1)
    vars(3) = 'prob'
    lvars(3) = 'p-value under a two-sided normal test'
    units(3) = '1'
    cell_methods(3) = cell_methods(1)
    ivars(1:2,3) = ivars(1:2,1)
    svars = ' '
    if ( index(file,'.ctl') /= 0 ) then
        ncid = 0
        call writectl(outfile,datfile,nx,xx,ny,yy,nz,zz,              &
     &           m2-m1+1,nperyear,1,m1,undef,title,nvars,vars,ivars   &
     &           ,lvars,units)
        open(1,file=trim(datfile),access='direct',form='unformatted'  &
     &           ,recl=recfa4*nx*ny*nz,status='new')
    else
        call enswritenc(outfile,ncid,ntvarid,itimeaxis,npermax,nx,xx  &
                ,ny,yy,nz,zz,lz,m2-m1+1,nperyear,1,m1,ltime,undef     &
                ,title,history,nvars,vars,ivars,lvars,svars,units     &
                ,cell_methods,metadata1,0,0)
    end if
    i=0
    if ( nx > nxf .or. ny > nyf ) then
        write(0,*) 'error: ',nxf,nyf,nx,ny
        call exit(-1)
    endif
    do j=m1,m2
        if ( lwrite ) print *,'writing ',nx*ny*nz,' reals at month ',j,' record ',i+1
        do jz=1,nz
            do jy=1,ny
                do jx=1,nx
                    if ( ave1(jx,jy,jz,j,1) < 1e30 .and. ave2(jx,jy,jz,j,1) < 1e30 ) then
                        if ( lnormsd ) then
                            if ( abs(ave1(jx,jy,jz,j,1)) > 1e-10 .and.         &
     &                               abs(ave2(jx,jy,jz,j,1)) > 1e-10 .and.     &
     &                               abs(ave2(jx,jy,jz,j,1)) > mindata ) then
                                s = ave1(jx,jy,jz,j,1)/ave2(jx,jy,jz,j,1)
                                ave1(jx,jy,jz,j,2) = s                                   &
     &                                   *sqrt(ave1(jx,jy,jz,j,2)/ave1(jx,jy,jz,j,1)**2  &
     &                                   + ave2(jx,jy,jz,j,2)/ave2(jx,jy,jz,j,1)**2)
                                ave1(jx,jy,jz,j,1) = s - 1
                            else
                                ave1(jx,jy,jz,j,:) = 3e33
                            end if
                        else
                            if ( abs(ave1(jx,jy,jz,j,1)) > mindata .and.          &
     &                               abs(ave2(jx,jy,jz,j,1)) > mindata ) then
                                ave1(jx,jy,jz,j,1) = ave1(jx,jy,jz,j,1) - ave2(jx,jy,jz,j,1)
                                if ( ave1(jx,jy,jz,j,2) < 0 ) then
                                    write(0,*) 'difffield: error: ave1(',jx,jy,jz,j,',2)<0 ' &
     &                                       ,ave1(jx,jy,jz,j,2)
                                    ave1(jx,jy,jz,j,2) = 3e33
                                else if ( ave2(jx,jy,jz,j,2) < 0 ) then
                                    write(0,*) 'difffield: error: ave2(',jx,jy,jz,j,',2)<0 ' &
     &                                       ,ave2(jx,jy,jz,j,2)
                                    ave1(jx,jy,jz,j,2) = 3e33
                                else if ( ave1(jx,jy,jz,j,2) > 1e30 .or.     &
     &                                    ave2(jx,jy,jz,j,2) > 1e30 ) then
                                    ave1(jx,jy,jz,j,2) = 3e33
                                else
                                    ave1(jx,jy,jz,j,2) = sqrt(ave1(jx,jy,jz,j,2) + ave2(jx,jy,jz,j,2))
                                end if
                            else
                                ave1(jx,jy,jz,j,:) = 3e33
                            end if
                        end if
                    else
                        ave1(jx,jy,jz,j,:) = 3e33
                    end if
                end do
            end do
        end do
        if ( lwrite ) then
            print *,'diff,sd(',nx/2,ny/2,(nz+1)/2,j,') = ',ave1(nx/1,ny/2,(nz+1)/2,j,1),  &
                         ave1(nx/1,ny/2,(nz+1)/2,j,2)
        end if
        if ( ncid == 0 ) then
            do pow=1,2
                i = i + 1
                write(1,rec=i) (((ave1(jx,jy,jz,j,pow),jx=1,nx),jy=1,ny),jz=1,nz)
            end do
        else
            do pow=1,2
                call writencslice(ncid,ntvarid,itimeaxis,npermax                &
     &                   ,ivars(1,pow),ave1(1,1,1,j,pow),nxf,nyf,nzf,nx         &
     &                   ,ny,nz,j-m1+1,0)
            end do
        end if
        do jz=1,nz
            do jy=1,ny
                do jx=1,nx
!                       z-value
                    if ( ave1(jx,jy,jz,j,2) > 0 .and. ave1(jx,jy,jz,j,2) < 1e30 ) then
                        d = ave1(jx,jy,jz,j,1)
                        e = ave1(jx,jy,jz,j,2)
                        ave1(jx,jy,jz,j,1) = d/e
!                           p-value, two-sided, for the time being
!                           normal, not t.
                        z = ave1(jx,jy,jz,j,1)
                        p = erfcc(abs(z)/sqrt(2.))
                        ave1(jx,jy,jz,j,1) = p
                        !!!print *,'d,e,z,p = ',d,e,z,p
                    else
                        ave1(jx,jy,jz,j,1) = 3e33
                    end if
                end do
            end do
        end do
        if ( ncid == 0 ) then
            i = i + 1
            write(1,rec=i) (((ave1(jx,jy,jz,j,1),jx=1,nx),jy=1,ny),jz=1,nz)
        else
            call writencslice(ncid,ntvarid,itimeaxis,npermax,             &
     &               ivars(1,3),ave1(1,1,1,j,1),nxf,nyf,nzf,nx,ny,nz,     &
     &               j-m1+1,0)
        end if
    end do
    if ( ncid /= 0 ) status = nf_close(ncid)
!
    end program
