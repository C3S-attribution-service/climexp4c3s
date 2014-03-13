        subroutine readncfile(ncid,field,nxf,nyf,nx,ny,nperyear
     +        ,yrbeg,yrend,firstyr,firstmo,nt,undef,lwrite,yr1,yr2,jvars
     +        )
*
*       read the data in a netcdf file into field, replacing undef
*       with 3e33 if necessary, restricting ourselves to the years 
*       yr1...yr2.
*
        implicit none
#include "netcdf.inc"
        integer ncid,nxf,nyf,nx,ny,nperyear,yrbeg,yrend,firstyr
     +        ,firstmo,nt,yr1,yr2,jvars(5)
        real field(nxf,nyf,nperyear,yrbeg:yrend),undef
        logical lwrite
        integer unit,jx,jy,i,j,k,l,start(4),count(4),stride(4),imap(4)
     +        ,status,startmo,ck,noleap,ilast,jlast,nn,jj,itoobig,
     +       ntoobig
        logical notfirst
        integer leap
        logical isnan
!
        itoobig = 0
        ntoobig = 1
        if ( lwrite ) then
            print *,'readncfile: reading yr1-yr2  ',yr1,yr2
            print *,'            yrbeg,yrend,nper ',yrbeg,yrend
     +            ,nperyear
            print *,'            field starts at  ',firstyr,firstmo
            print *,'            and has size     ',nx,ny,nt
            print *,'            varid,vardims    ',jvars
            print *,'            undef            ',undef
        endif
        if ( yr1.lt.yrbeg ) then
            write(0,*) 'readncfile: error: begin date before limit'
     +           ,yr1,yrbeg
            write(*,*) 'readncfile: error: begin date before limit'
     +           ,yr1,yrbeg
            call abort
        endif
        if ( yr2.gt.yrend ) then
            write(0,*) 'readncfile: error: end date after limit'
     +           ,yr2,yrend
            write(*,*) 'readncfile: error: end date after limit'
     +           ,yr2,yrend
            call abort
        endif
        if ( yr2.lt.yr1 ) then
            write(0,*) 'readncfile: error: end date before begin date'
     +           ,yr2,yr1
            write(*,*) 'readncfile: error: end date before begin date'
     +           ,yr2,yr1
            call abort
        endif
        notfirst = .TRUE.
        nn = nperyear/366
*
*       for completeness, these loops will almost never be executed
*
        if ( lwrite .and. yrbeg.lt.max(firstyr,yr1)) print
     +        '(a,i4,a,i4)','readncfile: zeroing years ',yrbeg,'-',
     +        max(firstyr,yr1)-1
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
        if ( firstyr.ge.yr1 ) then
            if ( lwrite.and.firstmo.gt.1 ) print '(a,i3)'
     +            ,'readncfile: zeroing months 1-',firstmo-1
            do j=1,firstmo-1
                do jy=1,max(ny,1)
                    do jx=1,max(nx,1)
                        field(jx,jy,j,firstyr) = 3e33
                    enddo
                enddo
            enddo
        endif
*
*       read data
*
	if ( lwrite ) print *,'reading from NetCDF unit ',ncid
        k = 0
        do i=1,4
            stride(i) = 1
        enddo
        if ( jvars(2).gt.0 ) then
            k = k + 1
            start(k) = 1
            count(k) = nx
        endif
        if ( jvars(3).gt.0 ) then
            k = k + 1
            start(k) = 1
            count(k) = ny
        endif
        if ( jvars(4).gt.0 ) then
            k = k + 1
            start(k) = 1
            count(k) = 1
        endif
        if ( jvars(5).gt.0 ) then
            k = k + 1
            if ( yr1.gt.firstyr ) then
                startmo = 1
                start(k) = 1 + nperyear*(yr1-firstyr) - firstmo + 1
                count(k) = min( nt - nperyear*(yr1-firstyr) + firstmo -
     +                1,nperyear*(yr2-yr1+1))
            else
                startmo = firstmo
                start(k) = 1
                count(k) = min( nt, nperyear*(yr2-yr1+1) - firstmo + 1)
            endif
        else
            startmo = 1
            if ( lwrite ) then
                print *,'readncfile: warning: time undefined in ncid '
     +               ,ncid
            end if
        endif
        if ( lwrite ) then
            print *,'readncfile: startvec = ',(start(i),i=1,k)
            print *,'            countvec = ',(count(i),i=1,k)
        endif
        if ( nxf.eq.nx .and. jvars(2).eq.1 .or. jvars(2).eq.0 ) then
            if ( nyf.eq.ny .and. jvars(3).eq.2 .or. jvars(3).eq.0 ) then
                if ( 366*nn.ne.nperyear ) then
                    if ( lwrite ) then
                        print *,'readncfile: calling ',
     +                       'nf_get_vara_real(',
     +                       ncid,jvars(1),start,count,',field(',1,1,
     +                       startmo,max(firstyr,yr1),')'
                        write(0,*) 'field(1) = ',
     +                       field(1,1,startmo,max(firstyr,yr1))
                        write(0,*) 'field(',
     +                       ((min(yrend,yr2)-max(firstyr,yr1))*nperyear
     +                       +nperyear-startmo+1)*ny*nx,')   = ',
     +                       field(nx,ny,nperyear,min(yrend,yr2))
                    end if
                    status = nf_get_vara_real(ncid,jvars(1),start,count
     +                   ,field(1,1,startmo,max(firstyr,yr1)))
                    if ( status.ne.nf_noerr ) call handle_err(status
     +                   ,'readncfile: nf_get_vara_real')
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
                        if ( i.gt.yrend ) then
                            write(0,*) 'readncfile: error: i>yrend: ',i
     +                           ,yrend,l,ck
                            exit
                        end if
                        if ( j/nn.eq.31+29 .and. leap(i).eq.1 ) then
                            do jj=0,nn-1
                                if ( lwrite ) print *,'setting ',i,j+jj,
     +                               ' to undefined'
                                do jy=1,max(ny,1)
                                    do jx=1,max(nx,1)
                                        field(jx,jy,j+jj,i) = undef
                                    enddo
                                enddo
                            end do
                            j = j + nn
                            noleap = noleap + nn
                        endif
                        call keepalive(l,ck)
                        if ( lwrite ) print *,'reading ',start(k),i,j
                        status = nf_get_vara_real(ncid,jvars(1)
     +                           ,start,count,field(1,1,j,i))
                        ilast = i
                        jlast = j
                        j = j + 1
                        if ( j.gt.nperyear ) then
                            j = j - nperyear
                            i = i + 1
                        endif
                        start(k) = start(k) + 1
                    enddo
                    count(k) = ck + noleap
                    nt = nt + noleap
                endif
            elseif ( jvars(3).eq.2 ) then
*               I have to construct a mapping vector to skip a bit of field
                imap(1) = 1
                imap(2) = nxf
                imap(3) = nxf*nyf ! z is either missing or 1
                imap(4) = nxf*nyf
                if ( lwrite ) then
                    print *,'            imapvec  = ',(imap(i),i=1,k)
                endif
                if ( 366*nn.ne.nperyear ) then
                    if ( lwrite ) print '(a)'
     +                   ,'readncfile: calling nf_get_varm_real 1'
                    status = nf_get_varm_real(ncid,jvars(1),start,count
     +                   ,stride,imap,field(1,1,startmo,max(firstyr,yr1)
     +                   ))
                    if ( status.ne.nf_noerr ) call handle_err(status,
     +                   'readncfile: nf_get_varm_real 1')
                else
                    if ( lwrite ) print *,'leap years...'
                    ck = count(k)
                    count(k) = 1
                    i = max(firstyr,yr1)
                    j = startmo
                    noleap = 0
                    do l=1,ck
                        if ( i.gt.yrend ) then
                            write(0,*) 'readncfile: error: i>yrend: ',i
     +                           ,yrend,l,ck
                            exit
                        end if
                        if ( j/nn.eq.31+29 .and. leap(i).eq.1 ) then
                            do jj=0,nn-1
                                if ( lwrite ) print *,'setting ',i,j+jj,
     +                               ' to undefined'
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
                        call keepalive(l,ck)
                        if ( lwrite ) print *,'reading ',start(k),i,j
                        status = nf_get_varm_real(ncid,jvars(1),start
     +                       ,count,stride,imap,field(1,1,j,i))
                        j = j + 1
                        if ( j.gt.nperyear ) then
                            j = j - nperyear
                            i = i + 1
                        endif
                        start(k) = start(k) + 1
                    enddo
                    count(k) = ck + noleap
                    nt = nt + noleap
                endif
            else
                write(0,*) 'readncfile: cannot handle jvars(3) != 2 yet'
                call abort
            endif
        else
            if ( jvars(2).eq.1 ) then
*               I have to construct a mapping vector to skip a bit of field
                imap(1) = 1
                imap(2) = nxf
                imap(3) = nxf*nyf ! z is either missing or 1
                imap(4) = nxf*nyf
                if ( lwrite ) then
                    print *,'            imapvec  = ',(imap(i),i=1,k)
                endif
                if ( nperyear.ne.366 ) then
                    if ( lwrite ) print '(a)'
     +                   ,'readncfile: calling nf_get_varm_real 2'
                    status = nf_get_varm_real(ncid,jvars(1),start,count
     +                   ,stride,imap,field(1,1,startmo,max(firstyr,yr1)
     +                   ))
                    if ( status.ne.nf_noerr ) call handle_err(status,
     +                   ' readncfile: nf_get_varm_real 2')
                else
                    print *,'leap years...'
                    ck = count(k)
                    count(k) = 1
                    i = max(firstyr,yr1)
                    j = startmo
                    noleap = 0
                    do l=1,ck
                        if ( j/nn.eq.31+29 .and. leap(i).eq.1 ) then
                            do jj=1,nn-1
                                if ( lwrite ) print *,'setting ',i,j+jj,
     +                               ' to undefined'
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
                        call keepalive(l,ck)
                        if ( lwrite ) print *,'reading ',start(k),i,j
                        status = nf_get_varm_real(ncid,jvars(1),start
     +                       ,count,stride,imap,field(1,1,j,i))
                        j = j + 1
                        if ( j.gt.nperyear ) then
                            j = j - nperyear
                            i = i + 1
                        endif
                        start(k) = start(k) + 1
                    enddo
                    count(k) = ck + noleap
                    nt = nt + noleap
                endif
            else
                write(0,*)
     +                'readncfile: cannot handle jvars(2) != 0,1 yet: '
     +                ,jvars(2)
                call abort
            endif
        endif
!       yes, Virginia, there are netcdf files in the wild that encoded undefs by NaNs...
        do i=max(yr1,yrbeg),min(yr2,yrend)
            do j=1,nperyear
                do jy=1,max(ny,1)
                    do jx=1,max(nx,1)
                        ! isnan does not catch some NaNs under pgf90...
                        if ( isnan(field(jx,jy,j,i)) .or. .not.
     +                       field(jx,jy,j,i).gt.3e33 .and.
     +                       field(jx,jy,j,i).lt.-3e33 ) then
                            field(jx,jy,j,i) = 3e33
                        end if
                    end do
                end do
            end do
        end do
!       change other undef conventions to ours
        if ( undef.ne.3e33 ) then
            if ( lwrite ) print *,'changing undef from ',undef,
     +           ' to 3e33'
            if ( abs(undef).gt.1000 ) then
                ! relative diff
                do i=max(yr1,yrbeg),min(yr2,yrend)
                    do j=1,nperyear
                        do jy=1,max(ny,1)
                            do jx=1,max(nx,1)
                                if ( abs(field(jx,jy,j,i)-undef).lt.
     +                               1e-6*abs(undef) ) then
                                    field(jx,jy,j,i) = 3e33
                                endif
                            enddo
                        enddo
                    enddo
                enddo
            else
                ! absolute diff
                do i=max(yr1,yrbeg),min(yr2,yrend)
                    do j=1,nperyear
                        do jy=1,max(ny,1)
                            do jx=1,max(nx,1)
                                if ( abs(field(jx,jy,j,i)-undef).lt.
     +                               1e-5 ) field(jx,jy,j,i) = 3e33
                            enddo
                        enddo
                    enddo
                enddo
            endif
        endif
*       
*       there may be a scale and/or offset in the file!
*
        call applyscaleoffset(ncid,jvars(1),field,nxf,nyf,1,nperyear,
     +       yrbeg,yrend,nx,ny,1,yr1,yr2,lwrite)
!
!       check for values that will cause crashes later on when taking squaares
!
        do i=max(yr1,yrbeg),min(yr2,yrend)
            do j=1,nperyear
                do jy=1,max(ny,1)
                    do jx=1,max(nx,1)
                        if ( field(jx,jy,j,i).lt.2e33 .and.
     +                       abs(field(jx,jy,j,i)).gt.1e16 ) then
                            itoobig = itoobig+1
                            if ( itoobig.ge.ntoobig ) then
                                ntoobig = 2*ntoobig
                                write(0,*) 'readncfile: warning: set ',
     +                               field(jx,jy,j,i),
     +                               ' to undefined ',itoobig
                            end if
                            field(jx,jy,j,i) = 3e33
                        end if
                    end do
                end do
            end do
        end do
*
*       for completeness, these loops will almost never be executed
*
        if ( nperyear.ne.366 ) then
            i = max(yr1,firstyr) + count(k)/nperyear
            k = startmo + mod(count(k),nperyear)
            if ( k.gt.nperyear ) then
                k = k-nperyear
                i = i + 1
            endif
        else
            k = jlast + 1
            i = ilast
        endif 
        if ( k.le.nperyear .and. i.le.yrend ) then
            if ( lwrite ) print '(a,i3,a,i3,a,i4)'
     +            ,'readncfile: zeroing months ',k,'-',nperyear
     +            ,' of year ',i
            do j=k,nperyear
                do jy=1,max(ny,1)
                    do jx=1,max(nx,1)
                        field(jx,jy,j,i) = 3e33
                    enddo
                enddo
            enddo
        endif
        if ( lwrite .and. i+1.le.yr2 ) print '(a,i4,a,i4)'
     +        ,'readncfile: zeroing years ',i+1,'-',min(yr2,yrend)
        do k=i+1,min(yr2,yrend)
            do j=1,nperyear
                do jy=1,max(ny,1)
                    do jx=1,max(nx,1)
                        field(jx,jy,j,k) = 3e33
                    enddo
                enddo
            enddo
        enddo
*
*       finito
*
        i = nf_close(ncid)
        ncid = 0                ! otherwise the next parsenc goes wrong.
        if ( lwrite ) then
            print *,'readncfile: field(',1+nx/2,1+ny/2,startmo
     +           ,max(firstyr,yr1),') = ',field(1+nx/2,1+ny/2,startmo
     +           ,max(firstyr,yr1))
        endif
        return
        end

        subroutine applyscaleoffset(ncid,jvars1,field,nxf,nyf,nzf
     +       ,nperyear,yrbeg,yrend,nx,ny,nz,yr1,yr2,lwrite)
!       
!       there may be a scale and/or offset in the file!
!
        implicit none
        integer ncid,jvars1,nxf,nyf,nzf,nperyear,yrbeg,yrend
        integer nx,ny,nz,yr1,yr2
        real field(nxf,nyf,nzf,nperyear,yrbeg:yrend)
        logical lwrite
        integer i,j,jx,jy,jz
        real scale,offset
!
        call getrealattopt(ncid,jvars1,'scale_factor',scale,lwrite)
        call getrealattopt(ncid,jvars1,'add_offset',offset,lwrite)
*       returns 3e33 if not found...
        if ( offset.gt.1e33 ) offset = 0
        if ( scale.lt.1e33 .and. scale.ne.1 ) then
            if ( offset.ne.0 ) then
!               going through all the data takes an awful lot of time,
!               so do not do it twice...
                if ( lwrite ) print *
     +               ,'readncfile: scaling and offsetting with ',scale
     +               ,offset
                call keepalive(0,0)
                do i=yr1,yr2
                    do j=1,nperyear
                        do jz=1,max(nz,1)
                            do jy=1,max(ny,1)
                                do jx=1,max(nx,1)
                                    if ( field(jx,jy,jz,j,i).lt.2e33 ) 
     +                                   field(jx,jy,jz,j,i) = 
     +                                   field(jx,jy,jz,j,i)*scale +
     +                                   offset
                                end do
                            end do
                        end do
                    end do
                end do
            else
                if ( lwrite ) print *,'readncfile: scaling with ',scale
                call keepalive(0,0)
                do i=yr1,yr2
                    do j=1,nperyear
                        do jz=1,max(nz,1)
                            do jy=1,max(ny,1)
                                do jx=1,max(nx,1)
                                    if ( field(jx,jy,jz,j,i).lt.2e33 ) 
     +                                   field(jx,jy,jz,j,i) = 
     +                                   field(jx,jy,jz,j,i)*scale
                                end do
                            end do
                        end do
                    end do
                end do
            endif
        else if ( offset.ne.0 ) then
            if ( lwrite ) print *,'readncfile: offsetting ',offset
            call keepalive(0,0)
            do i=yr1,yr2
                do j=1,nperyear
                    do jz=1,max(nz,1)
                        do jy=1,max(ny,1)
                            do jx=1,max(nx,1)
                                if ( field(jx,jy,jz,j,i).lt.2e33 ) 
     +                               field(jx,jy,jz,j,i) = 
     +                               field(jx,jy,jz,j,i) + offset
                            end do
                        end do
                    end do
                end do
            end do
        end if
        end

        subroutine fixholefield(field,nx,ny,nz,nperyear,firstyr
     +       ,lastyr,fyr,fmo,ntp,nt,tdefined,lwrite)
!
!       fixes a field that was read as if there were no holes in it
!       ntp = number of fields read in (sequentially)
!       nt  = differnce between first and last field + 1
!       tdefined = mask which fields are missing
!       (no possibility yet to throw away fields)
!
        implicit none
        integer nx,ny,nz,nperyear,firstyr,lastyr,fyr,fmo,ntp,nt
        real field(nx,ny,nz,nperyear,firstyr:lastyr)
        logical tdefined(nt),lwrite
        integer it,itp,ix,iy,iz,mo,mop,yr,yrp
!
!!!        lwrite = .true.
        itp = ntp
        do it=nt,1,-1
            mo = fmo + it - 1
            call normon(mo,fyr,yr,nperyear)
            if ( tdefined(it) ) then
                if ( itp.lt.it ) then
                    mop = fmo + itp - 1
                    call normon(mop,fyr,yrp,nperyear)
                    if ( lwrite ) print *,'copying from ',itp,mop,yrp,
     +                   ' to ',it,mo,yr
                    field(:,:,:,mo,yr) = field(:,:,:,mop,yrp)
                    itp = itp - 1
                end if
            else
                if ( lwrite ) print *,'setting ',it,mo,yr,
     +               ' to undefined'
                field(:,:,:,mo,yr) = 3e33
            end if
        end do
!!!        lwrite = .false.
        end
