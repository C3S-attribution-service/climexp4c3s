program netcdf2ascii

!   dump a netcdf file combo to ASCII
!   for the moment assume there is only one level.

    implicit none
    include 'params.h'
    include 'netcdf.inc'
    integer,parameter :: nvarmax=100,ntmax=5000000
    integer :: nx,ny,nz,nt,nperyear,yrbegin,mobegin,endian,nvars,ivars(2,nvarmax),jvars(6,nvarmax)
    integer :: status,ncid,i,j,k,start(4),count(4),nens1,nens2
    integer :: it,iv,ix,iy,yr,mo,dy,dpm(12,2),leap,localendian
    real :: xx(nxmax),yy(nymax),zz(1),undef
    real,allocatable :: data(:,:,:)
    character :: vars(nvarmax)*40,lvars(nvarmax)*120,units(nvarmax)*80
    character(1023) :: infile,datfile,title
    character :: lz(3)*20,svars(100)*100,ltime*120,history*50000, &
        cell_methods(100)*100,metadata(2,100)*2000
    logical :: tdefined(ntmax)
    logical :: lwrite
    integer :: iargc,llen,get_endian
    data dpm &
    /31,28,31,30,31,30,31,31,30,31,30,31 &
    ,31,29,31,30,31,30,31,31,30,31,30,31/

    lwrite = .false. 
    if ( iargc() /= 1 ) then
        write(0,*) 'Usage: netcdf2ascii infile'
        call exit(-1)
    end if
    localendian = get_endian()
    call getarg(1,infile)
    status = nf_open(infile,nf_nowrite,ncid)
    if ( status /= 0 ) then
        ncid = 0
        call parsectl(infile,datfile,nxmax,nx,xx,nymax,ny,yy,nzmax &
            ,nz,zz,nt,nperyear,yrbegin,mobegin,undef,endian,title &
            ,nvarmax,nvars,vars,ivars,lvars,units)
        history = ' '
        metadata = ' '
        cell_methods = ' '
        svars = ' '
        nens1 = 0
        nens2 = 0
        tdefined = .true.
        ltime = ' '
        lz = ' '
        open(1,file=datfile,access='direct',recl=4*nx*ny*nvars)
    else
        call ensparsenc(infile,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax &
            ,nz,zz,lz,nt,nperyear,yrbegin,mobegin,ltime,tdefined,ntmax &
            ,nens1,nens2,undef,title,history,nvarmax,nvars,vars,jvars &
            ,lvars,svars,units,cell_methods,metadata)
!           jvars(1,iv) = no of variable
!           jvars(2,iv) = no of first dimension (should be X)
!           jvars(3,iv) = no of second dimension (should be Y)
!           jvars(4,iv) = no of third dimension (should be Z)
!           jvars(5,iv) = no of fourth dimension (should be T)
!
!           Check that all variables are on the same grid...
        do iv=2,nvars
            do j=2,5
                if ( jvars(j,iv) /= jvars(j,iv-1) ) then
                    write(0,*) 'netcdf2ascii: error: can only handle all variables on the same grid '
                    write(0,*) iv-1,(jvars(k,iv-1),k=2,5)
                    write(0,*) iv,(jvars(k,iv),k=2,5)
                    call exit(-1)
                end if
            end do
        end do
        iv = 1
        k = 0
        if ( jvars(2,iv) > 0 ) then
            k = k + 1
            start(k) = 1
            count(k) = nx
        end if
        if ( jvars(3,iv) > 0 ) then
            k = k + 1
            start(k) = 1
            count(k) = ny
        end if
        if ( jvars(4,iv) > 0 ) then
            k = k + 1
            start(k) = 1
            count(k) = max(nz,1)
        end if
        k = k + 1 ! time is filled out later
    end if
    allocate(data(nx,ny,nvars))
    if ( lwrite ) print *,'endian,localendian = ',endian,localendian
    yr = yrbegin
    if ( nperyear == 12 ) then
        mo = mobegin
        dy = 1
    else if ( nperyear == 366 ) then
        mo = 1
        dy = mobegin
    end if
    call printmetadata(6,infile,' ',title,history,metadata)
    do iv=1,nvars
        if ( svars(iv) == ' ' ) then
            print '(11a)','# ',trim(vars(iv)),' [',trim(units(iv)),'] ',trim(lvars(iv))
        else
            print '(11a)','# ',trim(vars(iv)),' [',trim(units(iv)),'] ',trim(lvars(iv)),' (',trim(svars(iv)),')'
        end if
    end do
    do it=1,nt
        if ( nx > 1 .or. ny > 1 ) then
            print '(a,i4.4,a,i2.2,a,i2.2)','# date: ',yr,'-',mo,'-',dy
        end if
        if ( ncid == 0 ) then
            read(1,rec=it) (((data(ix,iy,iv),ix=1,nx),iy=1,ny),iv=1,nvars)
            if ( localendian*endian == -1 ) then
                if ( lwrite ) print *,'before ',(((data(ix,iy,iv),ix=1,nx),iy=1,ny),iv=1,nvars)
                do iv=1,nvars
                    do iy=1,ny
                        call swapbyte4(data(1,iy,iv),nx)
                    end do
                end do
                if ( lwrite ) print *,'after  ',(((data(ix,iy,iv),ix=1,nx),iy=1,ny),iv=1,nvars)
            end if
        else
            start(k) = it
            count(k) = 1
            do iv=1,nvars
                if ( lwrite ) then
                    print *,'calling nf_get_vara_real with'
                    print *,'ncid,jvars(1,',iv,') = ',ncid,jvars(1,iv)
                    print *,'start = ',(start(j),j=1,k)
                    print *,'count = ',(count(j),j=1,k)
                end if
                status = nf_get_vara_real(ncid,jvars(1,iv),start,count,data(1,1,iv))
            end do
        end if
        do iv=1,nvars
            do iy=1,ny
                do ix=1,nx
                    if ( data(ix,iy,iv) == undef ) data(ix,iy,iv) = 3e33
                end do
            end do
        end do
        if ( nx > 1 .or. ny > 1 ) then
            print '(1000a)','# longitude  latitude      ',(vars(iv),iv=1,nvars)
            do iy=1,ny
                do ix=1,nx
                    do iv=1,nvars
                        if ( data(ix,iy,iv) < 1e33 ) goto 101
                    end do
                    goto 102
                101 continue
                    do iv=1,nvars
                        if ( data(ix,iy,iv) > 1e33 ) data(ix,iy,iv) = -999.9
                    end do
                    print '(2f10.4,1000g20.8)',xx(ix),yy(iy),(data(ix,iy,iv),iv=1,nvars)
                102 continue
                end do
            end do
        else
            do iv=1,nvars
                if ( data(1,1,iv) < 1e33 ) goto 201
            end do
            goto 202
        201 continue
            do iv=1,nvars
                if ( data(1,1,iv) > 1e33 ) &
                data(1,1,iv) = -999.9
            end do
            if ( nperyear == 12 ) then
                print '(i4.4,i3,100g20.8)',yr,mo,(data(1,1,iv),iv=1,nvars)
            else
                print '(i4.4,2i3,g20.8)',yr,mo,dy,(data(1,1,iv),iv=1,nvars)
            end if
        202 continue
        end if
        if ( nperyear == 12 ) then
            mo = mo + 1
        else if ( nperyear == 1 ) then
            yr = yr + 1
        else if ( nperyear == 366 ) then
            dy = dy + 1
        else
            write(0,*) 'Sorry, nperyear = ',nperyear,' not yet implemented'
            call exit(-1)
        end if
        if ( mod(yr,4) /= 0 .or. mod(yr,100) == 0 .and. mod(yr,400) /= 0 ) then
            leap = 1
        else
            leap = 2
        end if
        if ( nperyear > 12 ) then
            if ( dy > dpm(mo,leap) ) then
                dy = dy - dpm(mo,leap)
                mo = mo + 1
            end if
        end if
        if ( mo > 12 ) then
            mo = mo - 12
            yr = yr + 1
        end if
    end do
end program netcdf2ascii

