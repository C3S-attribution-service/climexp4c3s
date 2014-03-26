        subroutine zreadncfile(ncid,field,nxf,nyf,nzf,nx,ny,nz,nperyear,
     +       yrbeg,yrend,firstyr,firstmo,nt,undef,lwrite,yr1,yr2,jvars)
*
*       read the data in a netcdf file into field, replacing undef
*       with 3e33 if necessary, restricting ourselves to the years 
*       yr1...yr2.
*
        implicit none
#include "netcdf.inc"
        integer ncid,nxf,nyf,nzf,nx,ny,nz,nperyear,yrbeg,yrend,firstyr
     +       ,firstmo,nt,yr1,yr2,jvars(5),ilast,jlast,k1,k2,k3,j1,j2,j3
     +       ,dy,mo
        real field(nxf,nyf,nzf,nperyear,yrbeg:yrend),undef
        logical lwrite
        integer unit,jx,jy,jz,i,j,k,l,start(4),count(4),stride(4),
     +       imap(4),status,startmo,ck,noleap,itoobig,ntoobig
        real scale,offset
        real,allocatable :: aux2(:,:)
        logical lwsave,notfirst
        integer,external :: leap,julday
        logical isnan ! ieee_is_nan is more standard but not yet implemented in pgf90 nor gfortran
        itoobig = 0
        ntoobig = 1
        lwsave = lwrite
***        lwrite = .TRUE.
***        lwrite = .FALSE.
        if ( lwrite ) then
            print *,'zreadncfile: reading yr1-yr2  ',yr1,yr2
            print *,'            yrbeg,yrend,nper ',yrbeg,yrend
     +            ,nperyear
            print *,'            field starts at  ',firstyr,firstmo
            print *,'            and has size     ',nx,ny,nz,nt
            print *,'            varid,vardims    ',jvars
            print *,'            undef            ',undef
        endif
        if ( firstyr.gt.10000 ) then
            write(0,*) 'zreadncfile: error: wrong value for firstyr '
     +           ,firstyr
            call abort
        end if
        if ( yr1.lt.yrbeg ) then
            write(0,*) 'zreadncfile: error: begin date before limit'
     +           ,yr1,yrbeg
            write(*,*) 'zreadncfile: error: begin date before limit'
     +           ,yr1,yrbeg
            call abort
        endif
        if ( yr2.gt.yrend ) then
            write(0,*) 'zreadncfile: error: end date after limit'
     +           ,yr2,yrend
            write(*,*) 'zreadncfile: error: end date after limit'
     +           ,yr2,yrend
            call abort
        endif
        if ( yr2.lt.yr1 ) then
            write(0,*) 'zreadncfile: error: end date before begin date'
     +           ,yr2,yr1
            write(*,*) 'zreadncfile: error: end date before begin date'
     +           ,yr2,yr1
            call abort
        endif
        notfirst = .TRUE.
*
*       for completeness, these loops will almost never be executed
*
        if ( lwrite .and. yrbeg.lt.max(firstyr,yr1)) print
     +        '(a,i4,a,i4)','zreadncfile: zeroing years ',yrbeg,'-',
     +        max(firstyr,yr1)-1
        do i=max(firstyr,yrbeg),yr1-1
            do j=1,nperyear
                do jz=1,max(nz,1)
                    do jy=1,max(ny,1)
                        do jx=1,max(nx,1)
                            field(jx,jy,jz,j,i) = 3e33
                        enddo
                    enddo
                enddo
            enddo
        enddo
        do i=yr1,firstyr-1
            do j=1,nperyear
                do jz=1,max(nz,1)
                    do jy=1,max(ny,1)
                        do jx=1,max(nx,1)
                            field(jx,jy,jz,j,i) = 3e33
                        enddo
                    enddo
                enddo
            enddo
        enddo
        if ( firstyr.ge.yr1 ) then
            if ( lwrite.and.firstmo.gt.1 ) print '(a,i3)'
     +            ,'zreadncfile: zeroing months 1-',firstmo-1
            do j=1,firstmo-1
                do jz=1,max(nz,1)
                    do jy=1,max(ny,1)
                        do jx=1,max(nx,1)
                            field(jx,jy,jz,j,firstyr) = 3e33
                        enddo
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
            count(k) = max(1,nx)
        endif
        if ( jvars(3).gt.0 ) then
            k = k + 1
            start(k) = 1
            count(k) = max(1,ny)
        endif
        if ( jvars(4).gt.0 ) then
            k = k + 1
            start(k) = 1
            count(k) = max(1,nz)
        endif
        if ( jvars(5).gt.0 ) then
            k = k + 1
            if ( yr1.gt.firstyr ) then
                startmo = 1
                if ( nperyear.ne.366 ) then
                    start(k) = 1 + nperyear*(yr1-firstyr) - firstmo + 1
                    count(k) = min(nt-start(k), nperyear*(yr2-yr1+1))
                else ! leap years :-(
                    if ( leap(firstyr).eq.2 .or. firstmo.le.60 ) then
                        call getdymo(dy,mo,firstmo,nperyear)
                    else
                        call getdymo(dy,mo,firstmo-1,nperyear)
                    end if
                    j1 = julday(mo,dy,firstyr)
                    j2 = julday(1,1,yr1)
                    start(k) = 1 + j2-j1
                    j3 = julday(12,31,yr2)
                    count(k) = min(nt-start(k),1+j3-j2)
                end if
            else
                startmo = firstmo
                start(k) = 1
                count(k) = min( nt, nperyear*(yr2-yr1+1) - firstmo + 1)
            endif
        else
            write(0,*) 'zreadncfile: error: time undefined in ncid '
     +           ,ncid
            call abort
        endif
        if ( lwrite ) then
            print *,'zreadncfile: startvec = ',(start(i),i=1,k)
            print *,'             countvec = ',(count(i),i=1,k)
        endif
!       Is the data contiguous?  This is much easier.
        k1 = 1
        if ( jvars(2).eq.0 ) then
            k2 = 1
            if ( jvars(3).eq.0 ) then
                k3 = 1
            else
                k3 = 2
            end if
        else
            k2 = 2
            if ( jvars(3).eq.0 ) then
                k3 = 2
            else
                k3 = 3
            end if
        end if
        if ( (nxf.eq.nx .and. (jvars(2).eq.k1 .or. jvars(2).eq.0) ).and.
     +       (nyf.eq.ny .and. (jvars(3).eq.k2 .or. jvars(3).eq.0) ).and.
     +       (nzf.eq.nz .and. (jvars(4).eq.k3 .or. jvars(4).eq.0) ) 
     +       ) then
            if ( nperyear.ne.366 ) then
                if ( lwrite ) print '(a,i8,i4,a,2i5,a)'
     +               ,'zreadncfile: calling nf_get_vara_real(',ncid
     +               ,jvars(1),',start,count,field(1,1,1,',startmo
     +               ,max(firstyr,yr1),'))'
                status = nf_get_vara_real(ncid,jvars(1),start,count
     +               ,field(1,1,1,startmo,max(firstyr,yr1)))
                if ( lwrite ) print '(a)'
     +               ,'zreadncfile: back from nf_get_vara_real'
                if ( status.ne.nf_noerr ) call handle_err(status
     +               ,'zreadncfile: nf_get_vara_real: ')
            else                ! leap years...
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
                    if ( j.eq.31+29 .and. leap(i).eq.1 ) then
                        if ( lwrite ) print *,'setting ',i,j,
     +                       ' to undefined'
                        do jy=1,max(ny,1)
                            do jx=1,max(nx,1)
                                do jz=1,max(nz,1)
                                    field(jx,jy,jz,j,i) = undef
                                enddo
                            enddo
                        enddo
                        j = j + 1
                        noleap = noleap + 1
                    endif
                    if ( l.gt.ck-noleap ) exit
                    call keepalive(l,ck)
                    if ( lwrite ) print *,'reading ',start(k),i,j
                    status = nf_get_vara_real(ncid,jvars(1),start,count
     +                   ,field(1,1,1,j,i))
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
        else if ( ( jvars(3).eq.0 .or. jvars(3).gt.jvars(2) ) .and.
     +            ( jvars(4).eq.0 .or. jvars(4).gt.jvars(3) 
     +                           .and. jvars(4).gt.jvars(2) ) .and.
     +            ( jvars(5).eq.0 .or. jvars(5).gt.jvars(4)
     +                           .and. jvars(5).gt.jvars(3)
     +                           .and. jvars(5).gt.jvars(2) ) ) then
!           the co-ordinates are ordered as per COARDS/CF conventions...
*           I have to construct a mapping vector to skip a bit of field
            k = 1
            imap(k) = 1
            if ( jvars(2).gt.0 ) then
                k = k + 1
                imap(k) = nxf
            end if
            if ( jvars(3).gt.0 ) then
                k = k + 1
                imap(k) = nxf*nyf
            end if
            if ( jvars(4).gt.0 ) then
                k = k + 1
                imap(k) = nxf*nyf*nzf
            end if
            if ( lwrite ) then
                print *,'            imapvec  = ',(imap(i),i=1,k)
            endif
            if ( nperyear.ne.366 ) then
                if ( lwrite ) print '(a)'
     +               ,'zreadncfile: calling nf_get_varm_real 1'
                status = nf_get_varm_real(ncid,jvars(1),start,count
     +               ,stride,imap,field(1,1,1,startmo,max(firstyr,yr1)))
                if ( status.ne.nf_noerr ) call handle_err(status,
     +               'zreadncfile: nf_get_varm_real')
            else
                if ( lwrite ) print *,'leap years...'
                ck = count(k)
                count(k) = 1
                i = max(firstyr,yr1)
                j = startmo
                noleap = 0
                do l=1,ck
                    if ( j.eq.31+29 .and. leap(i).eq.1 ) then
                        if ( lwrite ) print *,'setting ',i,j,
     +                       ' to undefined'
                        do jz=1,max(nz,1)
                            do jy=1,max(ny,1)
                                do jx=1,max(nx,1)
                                    field(jx,jy,jz,j,i) = undef
                                enddo
                            enddo
                        enddo
                        j = j + 1
                        noleap = noleap + 1
                    endif
                    call keepalive(l,ck)
                    if ( lwrite ) print *,'reading ',start(k),i,j
                    status = nf_get_varm_real(ncid,jvars(1),start
     +                   ,count,stride,imap,field(1,1,1,j,i))
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
        else if ( jvars(2).eq.0 .and. jvars(3).eq.2 .and.
     +            jvars(4).eq.1 .and. jvars(5).eq.3 .and.
     +            nperyear.ne.366 ) then
!           special case for HadCM3 amoc files that are lev,lat,time
!           ordered
            imap = 1
            allocate(aux2(nz,ny))
            i = start(1)
            start(1) = start(2)
            start(2) = i
            i = count(1)
            count(1) = count(2)
            count(2) = i
            ck = count(k)
            count(k) = 1
            i = max(firstyr,yr1)
            j = startmo
            do l=1,ck
                call keepalive(l,ck)
                if ( lwrite ) then
                    print *,'start  = ',start
                    print *,'count  = ',count
                end if
                status = nf_get_vara_real(ncid,jvars(1),start,count
     +                   ,aux2)
                do jz=1,nz
                    do jy = 1,ny
                        if ( .false. .and. aux2(jz,jy).ne.undef )
     +                       print *,'aux2(',jz,jy,') = ',aux2(jz,jy)
                        field(1,jy,jz,j,i) = aux2(jz,jy)
                    end do
                end do
                ilast = i
                jlast = j
                j = j + 1
                if ( j.gt.nperyear ) then
                    j = j - nperyear
                    i = i + 1
                endif
                start(k) = start(k) + 1
            enddo
            count(k) = ck
            deallocate(aux2)
        else
            if ( jvars(2).ne.0 ) print *,'X-axis is #',jvars(2)
            if ( jvars(3).ne.0 ) print *,'Y-axis is #',jvars(3)
            if ( jvars(4).ne.0 ) print *,'Z-axis is #',jvars(4)
            if ( jvars(5).ne.0 ) print *,'T-axis is #',jvars(5)
            write(*,*) 'zreadncfield: error: cannot handle '//
     +           'co-ordinates in another order than '//
     +           'time,lev,lat,lon (CF-convention)'
            write(0,*) 'zreadncfield: error: cannot handle '//
     +           'co-ordinates in another order than '//
     +           'time,lev,lat,lon (CF-convention)'
            call abort
        endif
!       yes, Virginia, there are netcdf files in the wild that encoded undefs by NaNs...
        do i=max(yr1,yrbeg),min(yr2,yrend)
            do j=1,nperyear
                do jz=1,max(nz,1)
                    do jy=1,max(ny,1)
                        do jx=1,max(nx,1)
                            ! isnan does not catch some NaNs under pgf90...
                            if ( isnan(field(jx,jy,jz,j,i)) .or. .not.
     +                           field(jx,jy,jz,j,i).gt.3e33 .and.
     +                           field(jx,jy,jz,j,i).lt.-3e33 ) then
                                field(jx,jy,jz,j,i) = 3e33
                            end if
                        end do
                    end do
                end do
            end do
        end do
!       change undef convention to ours
        if ( undef.ne.3e33 ) then
            if ( lwrite ) print *,'changing undefineds to 3e33'
            if ( abs(undef).gt.1000 ) then
                do i=max(yr1,yrbeg),min(yr2,yrend)
                    do j=1,nperyear
                        do jz=1,max(nz,1)
                            do jy=1,max(ny,1)
                                do jx=1,max(nx,1)
                                    if ( abs(field(jx,jy,jz,j,i)-undef)
     +                                   .lt.1e-6*abs(undef) ) then
                                        field(jx,jy,jz,j,i) = 3e33
                                    endif
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            else
                do i=max(yr1,yrbeg),min(yr2,yrend)
                    do j=1,nperyear
                        do jz=1,max(nz,1)
                            do jy=1,max(ny,1)
                                do jx=1,max(nx,1)
                                    if ( abs(field(jx,jy,jz,j,i)-undef)
     +                                   .lt.1e-5 ) then
                                        field(jx,jy,jz,j,i) = 3e33
                                    endif
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            endif
        endif
*       
*       there may be a scale and/or offset in the file!
*
        call applyscaleoffset(ncid,jvars(1),field,nxf,nyf,nzf,nperyear
     +       ,yrbeg,yrend,nx,ny,nz,yr1,yr2,lwrite)
!
!       check for values that will cause crashes later on when taking squaares
!
        do i=max(yr1,yrbeg),min(yr2,yrend)
            do j=1,nperyear
                do jz=1,max(nz,1)
                    do jy=1,max(ny,1)
                        do jx=1,max(nx,1)
                            if ( field(jx,jy,jz,j,i).lt.2e33 .and.
     +                           abs(field(jx,jy,jz,j,i)).gt.1e16 ) then
                                itoobig = itoobig+1
                                if ( itoobig.ge.ntoobig ) then
                                    ntoobig = 2*ntoobig
                                    write(0,*) 'zreadncfile: warning: ',
     +                                   'setting ',field(jx,jy,jz,j,i),
     +                                   ' to undefined ',itoobig
                                end if
                                field(jx,jy,jz,j,i) = 3e33
                            end if
                        end do
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
     +           ,'zreadncfile: zeroing months ',k,'-',nperyear
     +            ,' of year ',i
            do j=k,nperyear
                do jz=1,max(nz,1)
                    do jy=1,max(ny,1)
                        do jx=1,max(nx,1)
                            field(jx,jy,jz,j,i) = 3e33
                        enddo
                    enddo
                enddo
            enddo
        endif
        if ( lwrite .and. i+1.le.yr2 ) print '(a,i4,a,i4)'
     +       ,'zreadncfile: zeroing years ',i+1,'-',min(yr2,yrend)
        do k=i+1,min(yr2,yrend)
            do j=1,nperyear
                do jz=1,max(nz,1)
                    do jy=1,max(ny,1)
                        do jx=1,max(nx,1)
                            field(jx,jy,jz,j,k) = 3e33
                        enddo
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
            print *,'zreadncfile:'
            do i=yr1,yr2
                print *,'field(',(nx+1)/2,(ny+1)/2,(nz+1)/2,firstmo
     +               ,i,') = ',field((nx+1)/2,(ny+1)/2,(nz+1)
     +               /2,firstmo,i)
            enddo
        endif
        lwrite = lwsave
        return
        end
