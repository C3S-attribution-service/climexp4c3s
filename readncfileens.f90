subroutine readncfileens(ncid,field,nxf,nyf,nzf,nx,ny,nz,nperyear, &
    yrbeg,yrend,nensmax,firstyr,firstmo,nt,mens1,mens,undef,lwrite,yr1,yr2,jvars)

!   read the data in a netcdf file into field, replacing undef
!   with 3e33 if necessary, restricting ourselves to the years
!   yr1...yr2.

    implicit none
    include 'netcdf.inc'
    integer :: ncid,nxf,nyf,nzf,nx,ny,nz,nperyear,yrbeg,yrend,nensmax,firstyr &
        ,firstmo,nt,mens1,mens,yr1,yr2,jvars(6)
    real :: field(nxf,nyf,nzf,nperyear,yrbeg:yrend,0:nensmax),undef
    logical :: lwrite
    integer :: ilast,jlast,k1,k2,k3,k_ens,k_time,jul0,jul1,jul2,dy,mo,nn
    integer :: unit,jx,jy,jz,i,j,k,l,iens,start(4),count(4),stride(4), &
        imap(4),status,startmo,ck,noleap,itoobig,ntoobig
    real :: scale,offset
    real,allocatable :: aux2(:,:)
    logical :: lwsave,notfirst
    integer,external :: leap,julday
    logical :: isnan ! ieee_is_nan is more standard but not yet implemented in pgf90 nor gfortran
    itoobig = 0
    ntoobig = 1
    lwsave = lwrite
    !!!lwrite = .true.
    if ( lwrite ) then
        print *,'readncfileens: reading yr1-yr2  ',yr1,yr2
        print *,'            yrbeg,yrend,nper ',yrbeg,yrend,nperyear
        print *,'            field starts at  ',firstyr,firstmo
        print *,'            and has size     ',nx,ny,nz,nt,mens1,mens
        print *,'            varid,vardims    ',jvars
        print *,'            undef            ',undef
    endif
    if ( firstyr > 10000 ) then
        write(0,*) 'readncfileens: error: wrong value for firstyr ',firstyr
        call exit(-1)
    end if
    if ( yr1 < yrbeg ) then
        write(0,*) 'readncfileens: error: begin date before limit',yr1,yrbeg
        write(*,*) 'readncfileens: error: begin date before limit',yr1,yrbeg
        call exit(-1)
    endif
    if ( yr2 > yrend ) then
        write(0,*) 'readncfileens: error: end date after limit',yr2,yrend
        write(*,*) 'readncfileens: error: end date after limit',yr2,yrend
        call exit(-1)
    endif
    if ( yr2 < yr1 ) then
        write(0,*) 'readncfileens: error: end date before begin date',yr2,yr1
        write(*,*) 'readncfileens: error: end date before begin date',yr2,yr1
        call exit(-1)
    endif
    nn = max(1,nint(nperyear/365.24))
    notfirst = .true. 

!   put beginning to undefined

    if ( lwrite .and. yrbeg < max(firstyr,yr1)) print '(a,i4,a,i4)', &
        'readncfileens: zeroing years ',yrbeg,'-',max(firstyr,yr1)-1
    do iens=mens1,mens
        do i=max(firstyr,yrbeg),yr1-1
            do j=1,nperyear
                do jz=1,max(nz,1)
                    do jy=1,max(ny,1)
                        do jx=1,max(nx,1)
                            field(jx,jy,jz,j,i,iens) = 3e33
                        end do
                    end do
                end do
            end do
        end do
        do i=yr1,firstyr-1
            do j=1,nperyear
                do jz=1,max(nz,1)
                    do jy=1,max(ny,1)
                        do jx=1,max(nx,1)
                            field(jx,jy,jz,j,i,iens) = 3e33
                        end do
                    end do
                end do
            end do
        end do
    end do
    if ( firstyr >= yr1 ) then
        if ( lwrite .and. firstmo > 1 ) print '(a,i3)','readncfileens: zeroing months 1-',firstmo-1
        do iens=mens1,mens
            do j=1,firstmo-1
                do jz=1,max(nz,1)
                    do jy=1,max(ny,1)
                        do jx=1,max(nx,1)
                            field(jx,jy,jz,j,firstyr,iens) = 3e33
                        end do
                    end do
                end do
            end do
        end do
    endif

!   read data

    if ( lwrite ) print *,'reading from NetCDF unit ',ncid
    k = 0
    do i=1,4
        stride(i) = 1
    end do
    if ( jvars(2) > 0 ) then
        k = k + 1
        start(k) = 1
        count(k) = max(1,nx)
    endif
    if ( jvars(3) > 0 ) then
        k = k + 1
        start(k) = 1
        count(k) = max(1,ny)
    endif
    if ( jvars(4) > 0 ) then
        k = k + 1
        start(k) = 1
        count(k) = max(1,nz)
    endif
    if ( jvars(6) > 0 ) then
        k = k + 1
        start(k) = 1
        count(k) = mens-mens1+1
        k_ens = k
    endif
    if ( jvars(5) > 0 ) then
        k = k + 1
        k_time = k
        if ( nperyear /= 366 ) then
            if ( yr1 > firstyr ) then
                startmo = 1
                start(k) = 1 + nperyear*(yr1-firstyr) - firstmo + 1
                count(k) = min(nt-start(k)+1, nperyear*(yr2-yr1+1))
            else
                startmo = firstmo
                start(k) = 1
                count(k) = min(nt, nperyear*(yr2-yr1+1) - firstmo + 1)
            end if
        else
            if ( leap(firstyr) == 2 .or. firstmo <= 60 ) then
                call getdymo(dy,mo,firstmo,nperyear)
            else
                call getdymo(dy,mo,firstmo-1,nperyear)
            end if
            jul0 = julday(mo,dy,firstyr)
            jul1 = julday(1,1,yr1)
            jul2 = julday(12,31,yr2)
            if ( yr1 > firstyr ) then
                startmo = 1
                start(k) = nn*(jul1 - jul0) + 1
                count(k) = min(nt,jul2-jul0+1)
            else
                startmo = firstmo
                start(k) = 1
                count(k) = min( nt, nn*(jul2 - jul0 + 1))
            end if
        endif
    else
        startmo = 1
        if ( lwrite ) then
            print *,'readncfile: warning: time undefined in ncid ',ncid
        end if
    endif
    if ( lwrite ) then
        print *,'readncfileens: startvec = ',(start(i),i=1,k)
        print *,'               countvec = ',(count(i),i=1,k)
    endif
!   Is the data contiguous?  This is much easier.
!   teh time and ensemble dimensions are swapped between the netcdf and the array
!   so this does not work if there are multiple ensemble members in the netcdf.
    k1 = 1
    if ( jvars(2) == 0 ) then ! no X axis
        k2 = k1 ! Y is the first dimension in the netcdf file
    else
        k2 = k1 + 1 ! Y axis is second dimension with data
    end if
    if ( jvars(3) == 0 ) then
        k3 = k2
    else
        k3 = k2 + 1 
    end if
    
    if ( (nxf == nx .and. (jvars(2) == k1 .or. jvars(2) == 0) ) .and. &
         (nyf == ny .and. (jvars(3) == k2 .or. jvars(3) == 0) ) .and. &
         (nzf == nz .and. (jvars(4) == k3 .or. jvars(4) == 0) ) ) then
        if ( nperyear /= 366 .and. mens1 == mens ) then
            if ( lwrite ) print '(a,i8,i4,a,2i5,a)' &
                ,'readncfileens: calling nf_get_vara_real(',ncid &
                ,jvars(1),',start,count,field(1,1,1,',startmo &
                ,max(firstyr,yr1),'))'
            status = nf_get_vara_real(ncid,jvars(1),start,count &
                ,field(1,1,1,startmo,max(firstyr,yr1),mens1))
            if ( lwrite ) print '(a)','readncfileens: back from nf_get_vara_real'
            if ( status /= nf_noerr ) call handle_err(status &
                ,'readncfileens: nf_get_vara_real: ')
        else                ! leap years or multiple ensemble members
            if ( lwrite ) then
                if ( nperyear == 366 ) print *,'leap years...'
                if ( mens1 /= mens ) print *,'multiple ensemble members...'
                print *,'ck = ',count(k)
            endif
            ck = count(k)
            count(k) = 1
            i = max(firstyr,yr1)
            j = startmo
            noleap = 0
            do l=1,ck
                if ( nperyear == 366 .and. j == nn*(31+28)+1 .and. leap(i) == 1 ) then
                    if ( lwrite ) print *,'setting ',i,j,' to undefined'
                    do iens=mens1,mens
                        do jz=1,max(nz,1)
                            do jy=1,max(ny,1)
                                do jx=1,max(nx,1)
                                    field(jx,jy,jz,j,i,iens) = undef
                                end do
                            end do
                        end do
                    end do
                    j = j + 1
                    noleap = noleap + 1
                endif
                if ( nperyear == 366 ) then
                    call keepalive1("Reading day",l,ck)
                else
                    call keepalive1("Reading period",l,ck)
                end if
                start(k_ens) = 1
                count(k_ens) = 1
                do iens=mens1,mens
                    if ( .false. .and. lwrite ) then
                        print *,'reading ',i,j,iens
                        print *,'  start  = ',start
                        print *,'  count  = ',count
                        print *,'  stride = ',stride
                        print *,'  imap   = ',imap
                    end if
                    status = nf_get_vara_real(ncid,jvars(1),start,count &
                        ,field(1,1,1,j,i,iens))
                    start(k_ens) = start(k_ens) + 1
                end do
                ilast = i
                jlast = j
                j = j + 1
                if ( j > nperyear ) then
                    j = j - nperyear
                    i = i + 1
                endif
                start(k_time) = start(k_time) + 1
            end do
            count(k) = ck + noleap
            nt = nt + noleap
        endif
    else if ( ( jvars(3) == 0 .or. jvars(3) > jvars(2) ) .and. &
        ( jvars(4) == 0 .or. jvars(4) > jvars(3) &
         .and. jvars(4) > jvars(2) ) .and. &
        ( jvars(5) == 0 .or. jvars(5) > jvars(4) &
         .and. jvars(5) > jvars(3) &
         .and. jvars(5) > jvars(2) ) ) then
!       the co-ordinates are ordered as per COARDS/CF conventions...
!       I have to construct a mapping vector to skip a bit of field
        k = 1
        imap(k) = 1
        if ( jvars(2) > 0 ) then
            k = k + 1
            imap(k) = nxf
        end if
        if ( jvars(3) > 0 ) then
            k = k + 1
            imap(k) = nxf*nyf
        end if
        if ( jvars(4) > 0 ) then
            k = k + 1
            imap(k) = nxf*nyf*nzf
        end if
        if ( jvars(6) > 0 ) then
            k = k + 1
            imap(k) = nxf*nyf*nzf*(nensmax+1)
        end if
        if ( lwrite ) then
            print *,'            imapvec  = ',(imap(i),i=1,k)
        endif
        if ( nperyear /= 366 .and. mens1 == mens ) then
            if ( lwrite ) print '(a)','readncfileens: calling nf_get_varm_real 1'
            status = nf_get_varm_real(ncid,jvars(1),start,count &
                ,stride,imap,field(1,1,1,startmo,max(firstyr,yr1),mens1))
            if ( status /= nf_noerr ) call handle_err(status, &
            '   readncfileens: nf_get_varm_real')
        else
            if ( lwrite ) then
                if ( nperyear == 366 ) print *,'leap years...'
                if ( mens1 /= mens ) print *,'multiple ensemble members...'
            end if
            ck = count(k)
            count(k) = 1
            i = max(firstyr,yr1)
            j = startmo
            noleap = 0
            do l=1,ck
                if ( nperyear == 366 .and. j == 31+29 .and. leap(i) == 1 ) then
                    if ( lwrite ) print *,'setting ',i,j,' to undefined'
                    do iens=mens1,mens
                        do jz=1,max(nz,1)
                            do jy=1,max(ny,1)
                                do jx=1,max(nx,1)
                                    field(jx,jy,jz,j,i,iens) = undef
                                end do
                            end do
                        end do
                    end do
                    j = j + 1
                    noleap = noleap + 1
                endif
                call keepalive1('Reading day',l,ck)
                start(k_ens) = 1
                count(k_ens) = 1
                do iens=mens1,mens
                    if ( .false. .and. lwrite ) then
                        print *,'reading ',i,j,iens
                        print *,'  start  = ',start
                        print *,'  count  = ',count
                        print *,'  stride = ',stride
                        print *,'  imap   = ',imap
                    end if
                    status = nf_get_varm_real(ncid,jvars(1),start &
                        ,count,stride,imap,field(1,1,1,j,i,iens))
                    start(k_ens) = start(k_ens) + 1
                end do
                ilast = i
                jlast = j
                j = j + 1
                if ( j > nperyear ) then
                    j = j - nperyear
                    i = i + 1
                endif
                start(k_time) = start(k_time) + 1
            end do
            count(k) = ck + noleap
            nt = nt + noleap
        endif
    else
        if ( jvars(2) /= 0 ) print *,'X-axis is #',jvars(2)
        if ( jvars(3) /= 0 ) print *,'Y-axis is #',jvars(3)
        if ( jvars(4) /= 0 ) print *,'Z-axis is #',jvars(4)
        if ( jvars(5) /= 0 ) print *,'T-axis is #',jvars(5)
        write(*,*) 'zreadncfield: error: cannot handle '// &
            'co-ordinates in another order than '// &
            'time,lev,lat,lon (CF-convention)'
        write(0,*) 'zreadncfield: error: cannot handle '// &
            'co-ordinates in another order than '// &
            'time,lev,lat,lon (CF-convention)'
        call exit(-1)
    endif

!   for completeness, these loops will almost never be executed

    if ( nperyear /= 366 ) then
        i = max(yr1,firstyr) + count(k)/nperyear
        k = startmo + mod(count(k),nperyear)
        if ( k > nperyear ) then
            k = k-nperyear
            i = i + 1
        endif
    else
        k = jlast + 1
        i = ilast
    endif
    if ( k <= nperyear .and. i <= yrend ) then
        if ( lwrite ) print '(a,i3,a,i3,a,i4)','readncfileens: zeroing months ',k,'-',nperyear &
            ,' of year ',i
        do iens=mens1,mens
            do j=k,nperyear
                do jz=1,max(nz,1)
                    do jy=1,max(ny,1)
                        do jx=1,max(nx,1)
                            field(jx,jy,jz,j,i,iens) = 3e33
                        end do
                    end do
                end do
            end do
        end do
    end if
    if ( lwrite .and. i+1 <= yr2 ) print '(a,i4,a,i4)' &
        ,'readncfileens: zeroing years ',i+1,'-',min(yr2,yrend)
    do iens=mens1,mens
        do k=i+1,min(yr2,yrend)
            do j=1,nperyear
                do jz=1,max(nz,1)
                    do jy=1,max(ny,1)
                        do jx=1,max(nx,1)
                            field(jx,jy,jz,j,k,iens) = 3e33
                        end do
                    end do
                end do
            end do
        end do
    end do

!   change undef convention to ours

    if ( undef /= 3e33 ) then
        if ( lwrite ) print *,'changing undefineds to 3e33'
        if ( abs(undef) > 1000 ) then
            do iens=mens1,mens
                do i=max(yr1,yrbeg),min(yr2,yrend)
                    do j=1,nperyear
                        do jz=1,max(nz,1)
                            do jy=1,max(ny,1)
                                do jx=1,max(nx,1)
                                    if ( abs(field(jx,jy,jz,j,i,iens)-undef) < 1e-6*abs(undef) ) then
                                        field(jx,jy,jz,j,i,iens) = 3e33
                                    end if
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        else
            do iens=mens1,mens
                do i=max(yr1,yrbeg),min(yr2,yrend)
                    do j=1,nperyear
                        do jz=1,max(nz,1)
                            do jy=1,max(ny,1)
                                do jx=1,max(nx,1)
                                    if ( abs(field(jx,jy,jz,j,i,iens)-undef) < 1e-5 ) then
                                        field(jx,jy,jz,j,i,iens) = 3e33
                                    end if
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        endif
    endif

!   there may be a scale and/or offset in the file!

    call applyscaleoffsetens(ncid,jvars(1),field,nxf,nyf,nzf,nperyear &
        ,yrbeg,yrend,nensmax,nx,ny,nz,yr1,yr2,mens1,mens,lwrite)

!   check for values that will cause crashes later on when taking squaares

    do iens=mens1,mens
        do i=max(yr1,yrbeg),min(yr2,yrend)
            do j=1,nperyear
                do jz=1,max(nz,1)
                    do jy=1,max(ny,1)
                        do jx=1,max(nx,1)
                            if ( field(jx,jy,jz,j,i,iens) < 2e33 .and. &
                                    abs(field(jx,jy,jz,j,i,iens)) > 1e16 ) then
                                itoobig = itoobig+1
                                if ( itoobig >= ntoobig ) then
                                    ntoobig = 2*ntoobig
                                    write(0,*) 'readncfileens: warning: setting field(', &
                                        jx,jy,jz,j,i,iens,') = ',field(jx,jy,jz,j,i,iens), &
                                        ' to undefined ',itoobig
                                end if
                                field(jx,jy,jz,j,i,iens) = 3e33
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end do

!   finito

    i = nf_close(ncid)
    ncid = 0                ! otherwise the next parsenc goes wrong.
    if ( lwrite ) then
        print *,'readncfileens:'
        do i=yr1,yr2
            print *,'field(',(nx+1)/2,(ny+1)/2,(nz+1)/2,firstmo &
                ,i,mens1 + (mens-mens1)/2,') = ',field((nx+1)/2,(ny+1)/2,(nz+1) &
                /2,firstmo,i,mens1 + (mens-mens1)/2)
        end do
    end if
    lwrite = lwsave
end subroutine readncfileens
