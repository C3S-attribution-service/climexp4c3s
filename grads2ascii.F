        program grads2ascii
*
*       dump a GrADS ctl/dat file combo to ASCII
*       for the moment assume there is only one level.
*
        implicit none
#include "params.h"
#include "netcdf.inc"
        integer nvarmax
        parameter(nvarmax=100)
        integer nx,ny,nz,nt,nperyear,yrbegin,mobegin,endian,nvars
     +        ,ivars(2,nvarmax),jvars(6,nvarmax)
        integer status,ncid,i,j,k,start(4),count(4)
        integer it,iv,ix,iy,yr,mo,dy,dpm(12,2),leap,localendian
        real xx(nxmax),yy(nymax),zz(1),undef
        real,allocatable :: data(:,:,:)
        character vars(nvarmax)*10,lvars(nvarmax)*80,units(nvarmax)*40
        character*255 ctlfile,datfile,title
        logical lwrite
        integer iargc,llen,get_endian
        data dpm 
     +        /31,28,31,30,31,30,31,31,30,31,30,31
     +        ,31,29,31,30,31,30,31,31,30,31,30,31/

        lwrite = .false.
        if ( iargc().ne.1 ) then
            write(0,*) 'Usage: grads2ascii infile'
            stop
        endif
        localendian = get_endian()
        call getarg(1,ctlfile)
        status = nf_open(ctlfile,nf_nowrite,ncid)
        if ( status.ne.0 ) then
            ncid = 0
            call parsectl(ctlfile,datfile,nxmax,nx,xx,nymax,ny,yy,nzmax
     +           ,nz,zz,nt,nperyear,yrbegin,mobegin,undef,endian,title
     +           ,nvarmax,nvars,vars,ivars,lvars,units)
            open(1,file=datfile,access='direct',recl=4*nx*ny*nvars)
        else
            call parsenc(ctlfile,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax
     +        ,nz,zz,nt,nperyear,yrbegin,mobegin,undef,title,nvarmax
     +        ,nvars,vars,jvars,lvars,units)
!           jvars(1,iv) = no of variable
!           jvars(2,iv) = no of first dimension (should be X)
!           jvars(3,iv) = no of second dimension (should be X)
!           jvars(4,iv) = no of third dimension (should be Z)
!           jvars(5,iv) = no of fourth dimension (should be T)
!           Check that all variables are on the same grid...
            do iv=2,nvars
                do j=2,5
                    if ( jvars(j,iv).ne.jvars(j,iv-1) ) then
                        write(0,*) 'grads2ascii: error: can only handle'
     +                       ,' all variables on the same grid '
                        write(0,*) iv-1,(jvars(k,iv-1),k=2,5)
                        write(0,*) iv,(jvars(k,iv),k=2,5)
                        call abort
                    end if
                end do
            enddo
            iv = 1
            k = 0
            if ( jvars(2,iv).gt.0 ) then
                k = k + 1
                start(k) = 1
                count(k) = nx
            endif
            if ( jvars(3,iv).gt.0 ) then
                k = k + 1
                start(k) = 1
                count(k) = ny
            endif
            if ( jvars(4,iv).gt.0 ) then
                k = k + 1
                start(k) = 1
                count(k) = max(nz,1)
            endif
            k = k + 1 ! time is filled out later
        end if
        allocate(data(nx,ny,nvars))
        if ( lwrite ) print *,'endian,localendian = ',endian
     +       ,localendian
        yr = yrbegin
        if ( nperyear.eq.12 ) then
            mo = mobegin
            dy = 1
        elseif ( nperyear.eq.366 ) then
            mo = 1
            dy = mobegin
        endif
        print '(2a)','# ',title(1:llen(title))
        do iv=1,nvars
            print '(8a)','# ',trim(vars(iv)),' [',trim(units(iv)),'] '
     +           ,trim(lvars(iv))
        enddo
        do it=1,nt
            if ( nx.gt.1 .or. ny.gt.1 ) then
                print '(a,i4.4,a,i2.2,a,i2.2)','# date: ',yr,'-',mo,'-'
     +               ,dy
            endif
            if ( ncid.eq.0 ) then
                read(1,rec=it) (((data(ix,iy,iv),ix=1,nx),iy=1,ny),iv=1
     +            ,nvars)
                if ( localendian*endian.eq.-1 ) then
                    if ( lwrite ) print *,'before ',(((data(ix,iy,iv),ix
     +                   =1,nx),iy=1,ny),iv=1,nvars)
                    do iv=1,nvars
                        do iy=1,ny
                            call swapbyte4(data(1,iy,iv),nx)
                        end do
                    end do
                    if ( lwrite ) print *,'after  ',(((data(ix,iy,iv),ix
     +                   =1,nx),iy=1,ny),iv=1,nvars)
                end if
            else
                start(k) = it
                count(k) = 1
                do iv=1,nvars
                    if ( lwrite ) then
                        print *,'calling nf_get_vara_real with'
                        print *,'ncid,jvars(1,',iv,') = ',ncid,jvars(1
     +                       ,iv)
                        print *,'start = ',(start(j),j=1,k)
                        print *,'count = ',(count(j),j=1,k)
                    end if
                    status = nf_get_vara_real(ncid,jvars(1,iv),start
     +                   ,count,data(1,1,iv))
                end do
            end if
            do iv=1,nvars
                do iy=1,ny
                    do ix=1,nx
                        if ( data(ix,iy,iv).eq.undef )
     +                       data(ix,iy,iv) = 3e33
                    end do
                end do
            end do
            if ( nx.gt.1 .or. ny.gt.1 ) then
                print '(1000a)','# longitude  latitude      ',
     +               (vars(iv),iv=1,nvars)
                do iy=1,ny
                    do ix=1,nx
                        do iv=1,nvars
                            if ( data(ix,iy,iv).lt.1e33 ) goto 101
                        enddo
                        goto 102
 101                    continue
                        do iv=1,nvars
                            if ( data(ix,iy,iv).gt.1e33 ) data(ix,iy,iv)
     +                           = -999.9
                        end do
                        print '(2f10.4,1000g20.8)',xx(ix),yy(iy)
     +                       ,(data(ix,iy,iv),iv=1,nvars)
 102                    continue
                    enddo
                enddo
            else
                do iv=1,nvars
                    if ( data(1,1,iv).lt.1e33 ) goto 201
                enddo
                goto 202
 201            continue
                do iv=1,nvars
                    if ( data(1,1,iv).gt.1e33 )
     +                   data(1,1,iv) = -999.9
                end do
                if ( nperyear.eq.12 ) then
                    print '(i4.4,i3,100g20.8)',yr,mo,(data(1,1,iv),
     +                   iv=1,nvars)
                else
                    print '(i4.4,2i3,g20.8)',yr,mo,dy,(data(1,1,iv),
     +                   iv=1,nvars)
                endif
 202            continue
            endif
            if ( nperyear.eq.12 ) then
                mo = mo + 1
            elseif ( nperyear.eq.1 ) then
                yr = yr + 1
            elseif ( nperyear.eq.366 ) then
                dy = dy + 1
            else
                write(0,*) 'Sorry, nperyear = ',nperyear
     +                ,' not yet implemented'
                call abort
            endif
            if ( mod(yr,4).ne.0 .or. mod(yr,100).eq.0 .and. mod(yr,400)
     +            .ne.0 ) then
                leap = 1
            else
                leap = 2
            endif
            if ( nperyear.gt.12 ) then
                if ( dy.gt.dpm(mo,leap) ) then
                    dy = dy - dpm(mo,leap)
                    mo = mo + 1
                endif
            endif
            if ( mo.gt.12 ) then
                mo = mo - 12
                yr = yr + 1
            endif
        enddo
        end

