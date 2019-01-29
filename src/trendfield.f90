program trendfield

!       compute the derivative of a field over N steps

    use lsdata
    implicit none
    include 'params.h'
    include 'netcdf.inc'
    integer,parameter :: recfa4=4
    real :: absent
    parameter (absent=3e33)
    integer :: nvarmax,ntmax
    parameter(nvarmax=13,ntmax=5000)
    integer :: nyrmax,nlevmax,mensmax
    parameter (nyrmax=121,nlevmax=1,mensmax=1)
    integer :: nx,ny,nz,nt,firstyr,lastyr,firstmo,nvars, &
    ivars(2,nvarmax),jvars(6,nvarmax),ncid,endian, &
    status,nperyear,mens,mens1
    integer :: jx,jy,jz,i,mo,yr,it,ndiff,n, &
    iens,yrstart,yrstop,f,irec,year,ntvarid,itimeaxis(ntmax)
    logical :: lwrite,lexist
    real :: xx(nxmax),yy(nymax),zz(nzmax),undef,minfac
    real,allocatable :: field(:,:,:,:,:),fxy(:,:),diff(:,:)
    character infile*255,datfile*255,outfile*255,line*255,lz(3)*20 &
    ,vars(nvarmax)*40,outvars(nvarmax)*40,lvars(nvarmax)*255 &
    ,outlvars(nvarmax)*255,svars(nvarmax)*80,title*255,history &
    *10000,units(nvarmax)*80,outunits(nvarmax)*80 &
    ,cell_methods(nvarmax)*100,ltime*120,metadata(2,100)*2000
    character string*15

!       process command line

    lwrite = .false. 
    minfac = 0.5            ! arbitrary
    n = command_argument_count()
    if ( n < 3 ) then
        write(0,*) 'usage: trendfield infile.[ctl|nc] n [debug] ' &
        //'outfile.[ctl|nc]'
        stop
    endif
!       process arguments
    call get_command_argument(2,line)
    read(line,*) ndiff
    call get_command_argument(3,line)
    if ( line == 'debug' .or. line == 'lwrite' ) lwrite = .true. 
    call get_command_argument(1,infile)
    call getmetadata(infile,mens1,mens,ncid,datfile,nxmax,nx &
        ,xx,nymax,ny,yy,nzmax,nz,zz,lz,nt,nperyear,firstyr,firstmo &
        ,ltime,undef,endian,title,history,nvarmax,nvars,vars,ivars &
        ,lvars,svars,units,cell_methods,metadata,lwrite)
    call getlastyr(firstyr,firstmo,nt,nperyear,lastyr)
    if ( index(infile,'%') > 0 .or. &
    index(infile,'++') > 0 ) then
        write(0,*) 'Using ensemble members ',mens1,' to ',mens &
        ,'<br>'
    endif

!       allocate fields

    allocate(field(nx,ny,nz,nperyear,firstyr:lastyr))
    allocate(fxy(nperyear,firstyr:lastyr))
    allocate(diff(nperyear,firstyr:lastyr))

!       read field, change absent values to our convention

    do iens=mens1,mens
        call keepalive(iens-mens1+1,mens-mens1+1)
        if ( ncid == -1 ) then
            call get_command_argument(1,infile)
            if ( index(infile,'%') > 0 .or. &
            index(infile,'++') > 0 ) then
                call filloutens(infile,iens)
                if ( lwrite ) print *,'calling parsectl on ', &
                trim(infile)
            endif
            call parsectl(infile,datfile,nxmax,nx,xx,nymax,ny,yy &
                ,nzmax,nz,zz,nt,nperyear,f,firstmo,undef &
                ,endian,title,1,nvars,vars,ivars,lvars,units)
            if (lwrite) print '(2a)','# looking for ',trim(datfile)
            inquire(file=datfile,exist=lexist)
            if ( .not. lexist ) then
                if (lwrite) print '(3a)','# looking for ' &
                ,trim(datfile),'.gz'
                inquire(file=trim(datfile)//'.gz',exist=lexist)
                if ( .not. lexist ) then
                    mens = iens-1
                    if ( mens >= mens1 ) then
                        write(0,*) 'Found ensemble 0 to ',mens &
                        ,'<br>'
                        goto 999
                    else
                        write(0,*)'trendfield: error: cannot locate' &
                        ,' file ',trim(datfile)
                        call exit(-1)
                    endif
                endif
            endif
            if ( lwrite ) then
                print *,'opening file ',trim(datfile)
            endif
            call zreaddatfile(datfile,field(1,1,1,1,firstyr), &
                nx,ny,nz,nx,ny,nz,nperyear,firstyr,lastyr, &
                f,firstmo,nt,undef,endian,lwrite,firstyr,lastyr,1,1 &
                )
        else
            call get_command_argument(1,infile)
            if ( index(infile,'%') > 0 .or. &
            index(infile,'++') > 0 ) then
                call filloutens(infile,iens)
                if ( lwrite ) print *,'calling parsenc on ', &
                trim(infile)
                status = nf_open(infile,nf_nowrite,ncid)
            endif
            call parsenc(infile,ncid,nxmax,nx,xx,nymax,ny,yy &
                ,nzmax,nz,zz,nt,nperyear,f,firstmo &
                ,undef,title,1,nvars,vars,jvars,lvars,units)
            call zreadncfile(ncid,field(1,1,1,1,firstyr) &
                ,nx,ny,nz,nx,ny,nz,nperyear,firstyr,lastyr,f &
                ,firstmo,nt,undef,lwrite,firstyr,lastyr,jvars)
        endif
    
    !           compute derivatives
    
        do jz=1,nz
            do jy=1,ny
                do jx=1,nx
                    do yr=firstyr,lastyr
                        do mo=1,nperyear
                            fxy(mo,yr) = field(jx,jy,jz,mo,yr)
                        end do
                    end do
                    call derivative(ndiff,fxy,diff,nperyear,firstyr &
                    ,lastyr,nperyear,minfac,lwrite)
                    do yr=firstyr,lastyr
                        do mo=1,nperyear
                            field(jx,jy,jz,mo,yr) = diff(mo,yr)
                        end do
                    end do
                end do      ! jx
            end do          ! jy
        end do              ! jz
    
    !           write output
    
        call get_command_argument(1,infile)
        call get_command_argument(command_argument_count(),outfile)
        if ( index(infile,'%') > 0 .or. &
        index(infile,'++') > 0 ) then
            call filloutens(outfile,iens)
        end if
        outvars(1) = 'd'//trim(vars(1))//'dt'
        write(outlvars(1),'(i3,2a)') ndiff &
        ,'-point time derivative of ',trim(lvars(1))
        if ( nperyear == 1 ) then
            outunits(1) = '('//trim(units(1))//')/yr'
        else if ( nperyear == 12 ) then
            outunits(1) = '('//trim(units(1))//')/mo'
        else if ( nperyear < 360 ) then
            n = nint(366./nperyear)
            write(outunits(1),'(3a,i2,a)') '(',trim(units(1)),')/',n &
            ,'dy'
        else if  ( nperyear < 366 ) then
            outunits(1) = '('//trim(units(1))//')/dy'
        else if ( nperyear < 360 ) then
            n = nint(24*366./nperyear)
            write(outunits(1),'(3a,i3,a)') '(',trim(units(1)),')/',n &
            ,'hr'
        end if
        i = index(outfile,'.ctl')
        if ( i > 0 ) then
            datfile = outfile(1:i-1)//'.grd'
            call writectl(outfile,datfile,nx,xx,ny,yy,nz,zz, &
                nt,nperyear,firstyr,firstmo,3e33,title,nvars &
                ,outvars,ivars,outlvars,outunits)
            open(1,file=trim(datfile),form='unformatted', &
            access='direct',recl=recfa4*nx*ny*nz)
            yr = firstyr
            mo = firstmo
            do it=1,nt
                write(1,rec=it) (((field(jx,jy,jz,mo,yr),jx=1,nx), &
                jy=1,ny),jz=1,nz)
                mo = mo + 1
                if ( mo > nperyear ) then
                    mo = 1
                    yr = yr + 1
                end if
            end do
        else
            call enswritenc(outfile,ncid,ntvarid,itimeaxis,ntmax,nx &
                ,xx,ny,yy,nz,zz,lz,nt,nperyear,firstyr,firstmo &
                ,ltime,3e33,title,history,nvars,outvars,jvars &
                ,outlvars,svars,outunits,cell_methods,metadata,0,0)
            yr = firstyr
            mo = firstmo
            do it=1,nt
                call writencslice(ncid,ntvarid,itimeaxis,nt,jvars, &
                    field(1,1,1,mo,yr),nx,ny,nz,nx,ny,nz,it,1)
                mo = mo + 1
                if ( mo > nperyear ) then
                    mo = 1
                    yr = yr + 1
                end if
            end do
            status = nf_close(ncid)
        end if
    end do                  ! iens

999 continue
end program trendfield
