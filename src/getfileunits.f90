subroutine getfileunits(file,nx,ny,nz,nt,nperyear,nvarmax,nvars &
    ,var,units,newunits,lvar,svar,xwrap,lwrite)
    implicit none
    include 'netcdf.inc'
    include 'param.inc'
    integer,parameter :: mpermax=24*366
    integer :: nx,ny,nz,nt,nperyear,nvarmax,nvars
    character :: file*(*),var(nvarmax)*(*),units(nvarmax)*(*),newunits(nvarmax)*(*), &
        lvar(nvarmax)*(*),svar(nvarmax)*(*)
    logical :: xwrap,lwrite
    integer :: status,ncid,iu,i,i1,i2
    integer :: firstyr,firstmo,lastyr,lastmo,endian,ivars(2,100),jvars(6,100),ndpm, &
        nens1,nens2
    character :: datfile*255,title*255,history*50000,metadata(2,100)*2000
    character :: lz(3)*20,ltime*120,cell_methods(100)*100
    real :: xx(nxmax),yy(nymax),zz(nzmax),undef
    real :: data(mpermax,yrbeg:yrend),mean,slope,offset
    real,allocatable :: field(:,:,:,:,:)
    logical :: tdefined(ntmax)
    logical :: xrev,yrev,lexist,lcheck
    character line*255

    lcheck = .false. 
    if ( lwrite ) then
        print *,'getfileunits: called with '
        print *,' nvarmax = ',nvarmax
    endif
    units = ' '
    lvar = ' '
    svar = ' '
    if ( index(file,'++') /= 0 .or. index(file,'%') /= 0 ) then
        if ( lwrite ) print *,'getfileunits: ensemble'
        datfile = file
        call filloutens(file,0)
        inquire(file=file,exist=lexist)
        if ( .not. lexist ) then
            file = datfile
            call filloutens(file,1)
        endif
    endif
    if ( lwrite ) print *,'nf_open ',trim(file)
    status = nf_open(file,nf_nowrite,ncid)
    if ( status /= nf_noerr ) then
        if ( index(file,'.ctl ') /= 0 ) then
            if ( lwrite ) print *,'getfileunits: GrADS ctl file'
            call parsectl(file,datfile,nxmax,nx,xx,nymax,ny,yy,nzmax &
                ,nz,zz,nt,nperyear,firstyr,firstmo,undef,endian &
                ,title,nvarmax,nvars,var,ivars,lvar,units)
            if ( lwrite ) then
                do i=1,nvars
                    print *,i,trim(var(i)),trim(units(i))
                end do
            end if
            if ( lcheck ) then
!               read in part of the file to check units
                lastyr = firstyr + min(10,(firstmo+nt -2)/nperyear)
                nt = min(nt,nperyear*(lastyr-firstyr+1) - (firstmo-1))
                allocate(field(nx,ny,nz,nperyear,firstyr:lastyr))
                call zreaddatfile(datfile,field,nx,ny,nz,nx,ny,nz &
                    ,nperyear,firstyr,lastyr,firstyr,firstmo,nt &
                    ,undef,endian,lwrite,firstyr,lastyr,1,1)
                call estimatemean(field,nx,ny,nz,nperyear,firstyr &
                    ,lastyr,nx,ny,nz,nperyear,firstyr,lastyr,mean &
                    ,lwrite)
            else
!               I am not sure I want to read in part of the file
                if ( units(1) == 'K' ) then
                    mean = 273.13
                else
                    mean = 0
                endif
            end if
            call makestandardunits(mean,nperyear,var,units(1) &
                ,newunits(1),offset,slope,ndpm,lwrite)
            call getxyprop(xx,nx,yy,ny,xrev,yrev,xwrap)
        else
            if ( lwrite ) print *,'getfileunits: ascii time series'
            nx = 1
            ny = 1
            nz = 1
            nvars = 1
            xwrap = .false. 
            call rsunit(iu)
            open(iu,file=file,status='old',err=900)
            do i=1,1000
                read(iu,'(a)',end=100) line
                if ( line(1:1) /= '#' .and. line(2:2) /= '#' ) exit
                i1 = index(line,'[')
                i2 = index(line,']')
                if ( i1 > 0 .and. i2 > i1+1 ) then
                    units(1) = line(i1+1:i2-1)
                    call checkstring(units(1))
                    if ( lwrite ) then
                        print *,'found units= ',trim(units(1))
                    endif
                endif
            enddo
        100 continue
            close(iu)
            if ( lwrite ) print *,'calling readseries with '// &
                'mpermax,yrbeg,yrend,lwrite = ',mpermax,yrbeg,yrend,lwrite
            call readseriesmeta(file,data,mpermax,yrbeg,yrend, &
                nperyear,var,newunits(1),lvar,svar,history,metadata,.true.,lwrite)
            ! this is not exact for daily data with leap years
            do firstyr=yrbeg,yrend
                do firstmo=1,nperyear
                    if ( data(firstmo,firstyr) < 1e33 ) goto 110
                enddo
            enddo
            nt = 0
            goto 999
        110 continue
            do lastyr=yrend,yrbeg,-1
                do lastmo=nperyear,1,-1
                    if ( data(lastmo,lastyr) < 1e33 ) goto 120
                enddo
            enddo
        120 continue
            nt = 1 + (lastyr-firstyr)*nperyear - firstmo + lastmo
        endif
    else
        if ( lwrite ) print *,'getfileunits: netcdf file'
        nens1 = 0
        nens2 = 1
        call ensparsenc(file,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax &
            ,nz,zz,lz,nt,nperyear,firstyr,firstmo,ltime,tdefined,ntmax &
            ,nens1,nens2,undef,title,history,nvarmax,nvars,var,jvars,lvar,svar,units &
            ,cell_methods,metadata)
!       I am not sure I want to read in part of the file
        if ( units(1) == 'K' ) then
            mean = 273.13
        else
            mean = 0
        endif
        call makestandardunits(mean,nperyear,var,units(1) &
            ,newunits(1),offset,slope,ndpm,lwrite)
        call getxyprop(xx,nx,yy,ny,xrev,yrev,xwrap)
    endif
    if ( lwrite ) then
        print *,'units    = ',trim(units(1))
        print *,'newunits = ',trim(newunits(1))
    end if
    goto 999
900 print *,'getfileunits: error: cannot open file ',trim(file)
    call exit(-1)
999 continue
end subroutine getfileunits