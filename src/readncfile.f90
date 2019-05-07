subroutine readncfile(ncid,field,nxf,nyf,nx,ny,nperyear &
    ,yrbeg,yrend,firstyr,firstmo,nt,undef,lwrite,yr1,yr2,jvars)

!   read the data in a netcdf file into field, replacing undef
!   with 3e33 if necessary, restricting ourselves to the years
!   yr1...yr2.

    implicit none
    include 'netcdf.inc'
    integer :: ncid,nxf,nyf,nx,ny,nperyear,yrbeg,yrend,firstyr &
        ,firstmo,nt,yr1,yr2,jvars(6)
    real :: field(nxf,nyf,nperyear,yrbeg:yrend),undef
    logical :: lwrite
    integer :: unit,jx,jy,i,j,k,l,start(4),count(4),stride(4),imap(4) &
        ,status,startmo,ck,noleap,ilast,jlast,nn,jj,itoobig, &
        ntoobig,jul0,jul1,jul2,dy,mo,nperday,n
    logical :: notfirst
    integer :: leap,julday
    logical :: isnan

    itoobig = 0
    ntoobig = 1
    if ( lwrite ) then
        print *,'readncfile: reading yr1-yr2  ',yr1,yr2
        print *,'            yrbeg,yrend,nper ',yrbeg,yrend,nperyear
        print *,'            field starts at  ',firstyr,firstmo
        print *,'            and has size     ',nx,ny,nt
        print *,'            varid,vardims    ',jvars
        print *,'            undef            ',undef
    endif
    if ( yr1 < yrbeg ) then
        write(0,*) 'readncfile: error: begin date before limit',yr1,yrbeg
        write(*,*) 'readncfile: error: begin date before limit',yr1,yrbeg
        call exit(-1)
    endif
    if ( yr2 > yrend ) then
        write(0,*) 'readncfile: error: end date after limit' &
        ,yr2,yrend
        write(*,*) 'readncfile: error: end date after limit' &
        ,yr2,yrend
        call exit(-1)
    endif
    if ( yr2 < yr1 ) then
        write(0,*) 'readncfile: error: end date before begin date' &
        ,yr2,yr1
        write(*,*) 'readncfile: error: end date before begin date' &
        ,yr2,yr1
        call exit(-1)
    endif
    notfirst = .true. 
    nn = max(1,nint(nperyear/365.24))

!   for completeness, these loops will almost never be executed --- not true

    if ( lwrite .and. yrbeg < max(firstyr,yr1)) print '(a,i4,a,i4)', &
        'readncfile: zeroing years ',yrbeg,'-',max(firstyr,yr1)-1
    do i=max(firstyr,yrbeg),yr1-1
        do j=1,nperyear
            do jy=1,max(ny,1)
                do jx=1,max(nx,1)
                    field(jx,jy,j,i) = 3e33
                enddo
            enddo
        enddo
    enddo
    do i=yr1,firstyr-1
        do j=1,nperyear
            do jy=1,max(ny,1)
                do jx=1,max(nx,1)
                    field(jx,jy,j,i) = 3e33
                enddo
            enddo
        enddo
    enddo
    if ( firstyr >= yr1 ) then
        if ( lwrite .and. firstmo > 1 ) print '(a,i3)','readncfile: zeroing months 1-',firstmo-1
        do j=1,firstmo-1
            do jy=1,max(ny,1)
                do jx=1,max(nx,1)
                    field(jx,jy,j,firstyr) = 3e33
                enddo
            enddo
        enddo
    endif

!   read data

    if ( lwrite ) print *,'reading from NetCDF unit ',ncid
    k = 0
    do i=1,4
        stride(i) = 1
    enddo
    start = 0
    count = 0
    if ( jvars(2) > 0 ) then
        k = k + 1
        start(k) = 1
        count(k) = nx
    endif
    if ( jvars(3) > 0 ) then
        k = k + 1
        start(k) = 1
        count(k) = ny
    endif
    if ( jvars(4) > 0 ) then
        k = k + 1
        start(k) = 1
        count(k) = 1
    endif
    if ( jvars(6) > 0 ) then
    ! only read one ensemble member
        k = k + 1
        start(k) = 1
        count(k) = 1
    endif
    if ( jvars(5) > 0 ) then
        k = k + 1
        if ( nperyear/nn /= 366 ) then
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
            ! leap years...
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
                count(k) = nn*(jul2 - jul1 + 1)
            else
                startmo = firstmo
                start(k) = 1
                count(k) = min( nt, nn*(jul2 - jul0 + 1))
            endif
        end if
    else
        startmo = 1
        if ( lwrite ) then
            print *,'readncfile: warning: time undefined in ncid ',ncid
        end if
    endif
    if ( lwrite ) then
        print *,'readncfile: startvec = ',(start(i),i=1,k)
        print *,'            countvec = ',(count(i),i=1,k)
    endif
    if ( nxf == nx .and. jvars(2) == 1 .or. jvars(2) == 0 ) then
        if ( nyf == ny .and. jvars(3) == 2 .or. jvars(3) == 0 ) then
            if ( nperyear/nn /= 366 ) then
                if ( lwrite ) then
                    print *,'readncfile: calling nf_get_vara_real(', &
                        ncid,jvars(1),start,count,',field(',1,1, &
                        startmo,max(firstyr,yr1),'))'
                end if
                status = nf_get_vara_real(ncid,jvars(1),start,count &
                    ,field(1,1,startmo,max(firstyr,yr1)))
                if ( status /= nf_noerr ) call handle_err(status &
                    ,'readncfile: nf_get_vara_real: ')
            else            ! leap years...
                if ( lwrite ) then
                    print *,'leap years...'
                    print *,'ck = ',count(k)
                endif
                ck = count(k)
                count(k) = 1
                i = max(firstyr,yr1)
                j = startmo
                noleap = 0
                do l=1,ck
                    if ( i > yrend ) then
                        write(0,*) 'readncfile: error: i>yrend: ',i,yrend,l,ck
                        exit
                    end if
                    if ( j == nn*(31+28)+1 .and. leap(i) == 1 ) then
                        do jj=0,nn-1
                            if ( lwrite ) print *,'setting ',i,j+jj,' to undefined'
                            do jy=1,max(ny,1)
                                do jx=1,max(nx,1)
                                    field(jx,jy,j+jj,i) = undef
                                enddo
                            enddo
                        end do
                        j = j + nn
                        noleap = noleap + nn
                    endif
                    call keepalive1('Reading day ',l,ck)
                    if ( .false. .and. lwrite ) print *, &
                        'reading day ',l,'/',ck,':',start(k),i,j
                    status = nf_get_vara_real(ncid,jvars(1) &
                        ,start,count,field(1,1,j,i))
                    ilast = i
                    jlast = j
                    j = j + 1
                    if ( j > nperyear ) then
                        j = j - nperyear
                        i = i + 1
                        if ( i > min(yrend,yr2) ) exit
                    endif
                    start(k) = start(k) + 1
                enddo
                count(k) = ck + noleap
                nt = nt + noleap
            endif
        elseif ( jvars(3) == 2 ) then
!           I have to construct a mapping vector to skip a bit of field
            imap(1) = 1
            imap(2) = nxf
            imap(3) = nxf*nyf ! z is either missing or 1
            imap(4) = nxf*nyf
            if ( lwrite ) then
                print *,'            imapvec  = ',(imap(i),i=1,k)
            endif
            if ( 366*nn /= nperyear ) then
                if ( lwrite ) print '(a)','readncfile: calling nf_get_varm_real 1'
                status = nf_get_varm_real(ncid,jvars(1),start,count &
                    ,stride,imap,field(1,1,startmo,max(firstyr,yr1)))
                if ( status /= nf_noerr ) call handle_err(status,'readncfile: nf_get_varm_real 1')
            else
                if ( lwrite ) print *,'leap years...'
                ck = count(k)
                count(k) = 1
                i = max(firstyr,yr1)
                j = startmo
                noleap = 0
                do l=1,ck
                    if ( i > yrend ) then
                        write(0,*) 'readncfile: error: i>yrend: ',i,yrend,l,ck
                        exit
                    end if
                    if ( j/nn == 31+29 .and. leap(i) == 1 ) then
                        do jj=0,nn-1
                            if ( lwrite ) print *,'setting ',i,j+jj,' to undefined'
                            do jy=1,max(ny,1)
                                do jx=1,max(nx,1)
                                    field(jx,jy,j+jj,i) = undef
                                enddo
                            enddo
                        end do
                        ilast = i
                        jlast = j + nn - 1
                        j = j + nn
                        noleap = noleap + nn
                    endif
                    call keepalive1('Reading day ',l,ck)
                    if ( lwrite ) print *,'reading ',start(k),i,j
                    status = nf_get_varm_real(ncid,jvars(1),start &
                        ,count,stride,imap,field(1,1,j,i))
                    j = j + 1
                    if ( j > nperyear ) then
                        j = j - nperyear
                        i = i + 1
                        if ( i > min(yrend,yr2) ) exit
                    endif
                    start(k) = start(k) + 1
                enddo
                count(k) = ck + noleap
                nt = nt + noleap
            endif
        else
            write(0,*) 'readncfile: cannot handle jvars(3) != 2 yet'
            call exit(-1)
        endif
    else
        if ( jvars(2) == 1 ) then
!           I have to construct a mapping vector to skip a bit of field
            imap(1) = 1
            imap(2) = nxf
            imap(3) = nxf*nyf ! z is either missing or 1
            imap(4) = nxf*nyf
            if ( lwrite ) then
                print *,'            imapvec  = ',(imap(i),i=1,k)
            endif
            if ( nperyear /= 366 ) then
                if ( lwrite ) print '(a)' &
                ,'readncfile: calling nf_get_varm_real 2'
                status = nf_get_varm_real(ncid,jvars(1),start,count &
                    ,stride,imap,field(1,1,startmo,max(firstyr,yr1)))
                if ( status /= nf_noerr ) call handle_err(status,' readncfile: nf_get_varm_real 2')
            else
                print *,'leap years...'
                ck = count(k)
                count(k) = 1
                i = max(firstyr,yr1)
                j = startmo
                noleap = 0
                do l=1,ck
                    if ( j/nn == 31+29 .and. leap(i) == 1 ) then
                        do jj=1,nn-1
                            if ( lwrite ) print *,'setting ',i,j+jj, &
                            ' to undefined'
                            do jy=1,max(ny,1)
                                do jx=1,max(nx,1)
                                    field(jx,jy,j+jj,i) = undef
                                enddo
                            enddo
                        end do
                        ilast = i
                        jlast = j + nn - 1
                        j = j + nn
                        noleap = noleap + nn
                    endif
                    call keepalive1('Reading day ',l,ck)
                    if ( lwrite ) print *,'reading ',start(k),i,j
                    status = nf_get_varm_real(ncid,jvars(1),start &
                        ,count,stride,imap,field(1,1,j,i))
                    j = j + 1
                    if ( j > nperyear ) then
                        j = j - nperyear
                        i = i + 1
                    endif
                    start(k) = start(k) + 1
                enddo
                count(k) = ck + noleap
                nt = nt + noleap
            endif
        else
            write(0,*) 'readncfile: cannot handle jvars(2) != 0,1 yet: ',jvars(2)
            call exit(-1)
        endif
    endif

!   set the end to undefined

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
        if ( lwrite ) print '(a,i3,a,i3,a,i4)' &
            ,'readncfile: zeroing months ',k,'-',nperyear,' of year ',i
        do j=k,nperyear
            do jy=1,max(ny,1)
                do jx=1,max(nx,1)
                    field(jx,jy,j,i) = 3e33
                enddo
            enddo
        enddo
    endif
    if ( lwrite .and. i+1 <= yr2 ) print '(a,i4,a,i4)' &
        ,'readncfile: zeroing years ',i+1,'-',min(yr2,yrend)
    do k=i+1,min(yr2,yrend)
        do j=1,nperyear
            do jy=1,max(ny,1)
                do jx=1,max(nx,1)
                    field(jx,jy,j,k) = 3e33
                enddo
            enddo
        enddo
    enddo

!   yes, Virginia, there are netcdf files in the wild that encoded undefs by NaNs...
!   Pity, these loops eat a lot of time.
!   Made sure these do not occur in the Climate Explorer to save a lot of processing time.
    if ( .false. ) then
    n = 0
!$omp parallel do private(i,j,jy,jx)
    do i=max(yr1,yrbeg),min(yr2,yrend)
!$omp atomic update
        n = n + 1
        call keepalive1('Checking for weird undefs ',n,min(yr2,yrend)-max(yr1,yrbeg)+1)
        do j=1,nperyear
            do jy=1,max(ny,1)
                do jx=1,max(nx,1)
                    ! isnan does not catch some NaNs under pgf90...
                    if ( isnan(field(jx,jy,j,i)) .or. .not. &
                        field(jx,jy,j,i) > 3e33 .and. &
                        field(jx,jy,j,i) < -3e33 ) then
                            field(jx,jy,j,i) = 3e33
                    end if
                end do
            end do
        end do
    end do
!$omp end parallel do
    end if

!   change other undef conventions to ours

    if ( undef /= 3e33 ) then
        if ( lwrite ) print *,'changing undef from ',undef,' to 3e33'
        if ( abs(undef) > 1000 ) then
            ! relative diff
            do i=max(yr1,yrbeg),min(yr2,yrend)
                call keepalive1('Changing undefs to 3e33', &
                i-max(yr1,yrbeg)+1, &
                min(yr2,yrend)-max(yr1,yrbeg)+1)
                do j=1,nperyear
                    do jy=1,max(ny,1)
                        do jx=1,max(nx,1)
                            if ( abs(field(jx,jy,j,i)-undef) < 1e-6*abs(undef) ) then
                                field(jx,jy,j,i) = 3e33
                            endif
                        enddo
                    enddo
                enddo
            enddo
        else
        ! absolute diff
            do i=max(yr1,yrbeg),min(yr2,yrend)
                call keepalive1('Changing undefs to 3e33', &
                i-max(yr1,yrbeg)+1, &
                min(yr2,yrend)-max(yr1,yrbeg)+1)
                do j=1,nperyear
                    do jy=1,max(ny,1)
                        do jx=1,max(nx,1)
                            if ( abs(field(jx,jy,j,i)-undef) < 1e-5 ) field(jx,jy,j,i) = 3e33
                        enddo
                    enddo
                enddo
            enddo
        endif
    endif

!   there may be a scale and/or offset in the file!

    call applyscaleoffset(ncid,jvars(1),field,nxf,nyf,1,nperyear, &
        yrbeg,yrend,nx,ny,1,yr1,yr2,lwrite)

!   check for values that will cause crashes later on when taking squaares

    n = 0
!$omp parallel do private(i,j,jy,jx)
    do i=max(yr1,yrbeg),min(yr2,yrend)
!$omp atomic update
        n = n + 1
        call keepalive1('Checking a few other things ',n,min(yr2,yrend)-max(yr1,yrbeg)+1)
        do j=1,nperyear
            do jy=1,max(ny,1)
                do jx=1,max(nx,1)
                    if ( field(jx,jy,j,i) < 2e33 .and. abs(field(jx,jy,j,i)) > 1e16 ) then
                        itoobig = itoobig+1
                        if ( itoobig >= ntoobig ) then
                            ntoobig = 2*ntoobig
                            write(0,*) 'readncfile: warning: set ',field(jx,jy,j,i), &
                                ' to undefined ',itoobig
                        end if
                        field(jx,jy,j,i) = 3e33
                    end if
                end do
            end do
        end do
    end do

!   finito

    i = nf_close(ncid)
    ncid = 0                ! otherwise the next parsenc goes wrong.
    if ( lwrite ) then
        print *,'readncfile: field(',1+nx/2,1+ny/2,startmo &
            ,max(firstyr,yr1),') = ',field(1+nx/2,1+ny/2,startmo &
            ,max(firstyr,yr1))
    endif
    return
end subroutine readncfile

subroutine fixholefield(field,nx,ny,nz,nperyear,firstyr &
    ,lastyr,fyr,fmo,ntp,nt,tdefined,lwrite)

!   fixes a field that was read as if there were no holes in it
!   ntp = number of fields read in (sequentially)
!   nt  = differnce between first and last field + 1
!   tdefined = mask which fields are missing
!   (no possibility yet to throw away fields)

    implicit none
    integer :: nx,ny,nz,nperyear,firstyr,lastyr,fyr,fmo,ntp,nt
    real :: field(nx,ny,nz,nperyear,firstyr:lastyr)
    logical :: tdefined(nt),lwrite
    integer :: it,itp,ix,iy,iz,mo,mop,yr,yrp

!!!        lwrite = .true.
    itp = ntp
    do it=nt,1,-1
        mo = fmo + it - 1
        call normon(mo,fyr,yr,nperyear)
        if ( tdefined(it) ) then
            if ( itp < it ) then
                mop = fmo + itp - 1
                call normon(mop,fyr,yrp,nperyear)
                if ( lwrite ) print *,'copying from ',itp,mop,yrp,' to ',it,mo,yr
                field(:,:,:,mo,yr) = field(:,:,:,mop,yrp)
                itp = itp - 1
            end if
        else
            if ( lwrite ) print *,'setting ',it,mo,yr,' to undefined'
            field(:,:,:,mo,yr) = 3e33
        end if
    end do
!!!        lwrite = .false.
end subroutine fixholefield
