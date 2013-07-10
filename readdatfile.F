        subroutine readdatfile(datfile,field,nxf,nyf,nx,ny,nperyear
     +        ,yrbeg,yrend,firstyr,firstmo,nt,undef,endian,lwrite,yr1
     +        ,yr2,ivar,nvar)
*
*       read the data in a GrADS .dat file into field, replacing undef
*       with 3e33 if necessary, restricting ourselves to the years 
*       yr1...yr2.
*
        implicit none
#include "recfac.h"
        integer nxf,nyf,nx,ny,nperyear,yrbeg,yrend,firstyr,firstmo,nt
     +        ,endian,yr1,yr2,ivar,nvar,n
        real field(nxf,nyf,nperyear,yrbeg:yrend),undef
        character*(*) datfile
        logical lwrite,notfirst,lexist,deletefile
        integer unit,jx,jy,i,ii,j,jj,k,localendian,retval,noleap,ierr
     +       ,nerr,nn
        character*255 tmpfile,command*512
        integer llen,getpid,get_endian,leap
*
        ierr = 0
        nerr = 1
        localendian = get_endian()
        if ( lwrite ) then
            print *,'readdatfile called with parameters'
            print *,'datfile = ',datfile(1:llen(datfile))
            print *,'nxf,nyf = ',nxf,nyf
            print *,'nx,ny   = ',nx,ny
            print *,'yrbeg,end ',yrbeg,yrend
            print *,'firstyr,mo',firstyr,firstmo
            print *,'nt,undef= ',nt,undef
            print *,'endian,loc',endian,localendian
            print *,'yr1,yr2 = ',yr1,yr2
            print *,'ivar,nvar ',ivar,nvar
        endif
        if ( yr1.lt.yrbeg ) then
            write(0,*) 'readdatfile: error: begin date before limit'
     +           ,yr1,yrbeg
            write(*,*) 'readdatfile: error: begin date before limit'
     +           ,yr1,yrbeg
            call abort
        endif
        if ( yr2.gt.yrend ) then
            write(0,*) 'readdatfile: error: end date after limit'
     +           ,yr2,yrend
            write(*,*) 'readdatfile: error: end date after limit'
     +           ,yr2,yrend
            call abort
        endif
        if ( yr2.lt.yr1 ) then
            write(0,*) 'readdatfile: error: end date before begin date'
     +           ,yr2,yr1
            write(*,*) 'readdatfile: error: end date before begin date'
     +           ,yr2,yr1
            call abort
        endif
        notfirst = .TRUE.
*
*       for completeness, these loops will almost never be executed
*
        do i=max(firstyr,yrbeg),yr1-1
            do j=1,nperyear
                do jy=1,ny
                    do jx=1,nx
                        field(jx,jy,j,i) = 3e33
                    enddo
                enddo
            enddo
        enddo
        do i=yr1,firstyr-1
            do j=1,nperyear
                do jy=1,ny
                    do jx=1,nx
                        field(jx,jy,j,i) = 3e33
                    enddo
                enddo
            enddo
        enddo
        if ( firstyr.ge.yr1 ) then
            do j=1,firstmo-1
                do jy=1,ny
                    do jx=1,nx
                        field(jx,jy,j,firstyr) = 3e33
                    enddo
                enddo
            enddo
        endif
*
*       read data
*
        do i=len(datfile),1,-1
            if ( datfile(i:i).eq.'/' ) goto 200
        enddo
  200   continue
        i = i + 1
	if ( lwrite ) print '(3a,i2)','# reading file ',trim(datfile)
     +       ,' with record length ',recfac
        call rsunit(unit)
        open(unit=unit,file=trim(datfile),form='unformatted',
     +       access='direct',recl=recfa4*nx*ny,status='old',err=300)
        deletefile=.FALSE.
        goto 310
  300   continue
*
*       also allow gzipped files
*       
	datfile(index(datfile,' '):) = '.gz'
        if ( lwrite ) print '(2a)','trying ',trim(datfile)
        inquire(file=trim(datfile),exist=lexist)
        if ( .not.lexist ) goto 904
        n = 0
  305   continue
        if ( index(datfile,'SOS').ne.0 ) then
            write(tmpfile,'(a,i16.16,a)') 'data/readdatfile',getpid()+n
     +        ,'.dat'
        else
            write(tmpfile,'(a,i16.16,a)') '/tmp/readdatfile',getpid()+n
     +        ,'.dat'
        endif
        inquire(file=tmpfile,exist=lexist)
        if ( lexist ) then
            n = n + 1
            goto 305
        endif
        write(command,'(4a)') 'gunzip -c ',trim(datfile),' > '
     +       ,trim(tmpfile)
        if ( lwrite ) print '(a)',trim(command)
        call mysystem(trim(command),retval)
        if ( lwrite ) print *,'retval = ',retval
        if ( retval.ne.0 ) then
            print '(2a)','# ',trim(command)
            print '(a,i10)','# readdatfile: retval = ',retval
            call mysystem(trim(command),retval)
        endif
	if ( lwrite ) print '(3a,i2)','reading file ',trim(tmpfile)
     +       ,' with record length ',recfac
        open(unit=unit,file=trim(tmpfile),form='unformatted',
     +       access='direct',recl=recfac*nx*ny,status='old',err=904)
!       to avoid leaving files behind if the process crashes
        write(command,'(2a)') 'rm ',trim(tmpfile)
        call mysystem(command,retval)
        deletefile=.TRUE.
  310   continue
        i = firstyr
        j = firstmo - 1
        noleap = 0
        do k=1,nt
            j = j + 1
            if ( j.gt.nperyear ) then
                j = j - nperyear
                i = i + 1
            endif
            if ( i.ge.yr1 .and. i.le.yr2 ) then
                nn = (nperyear/366)
                if ( 366*nn.eq.nperyear .and. leap(i).eq.1 ) then
                    if ( j/nn.eq.31+29 ) then
                        if ( lwrite ) print *,'setting Feb 29 of ',i,
     +                       ' to undefined'
                        do jj=0,nn-1
                            do jy=1,ny
                                do jx=1,nx
                                    field(jx,jy,j+jj,i) = 3e33
                                enddo
                            enddo
                        end do
                        j = j + nn
                        noleap = noleap + nn
                    endif
                endif
                if ( lwrite ) print *,'reading field ',ivar
     +                +(k-1)*nvar,i,j
*               last field that was read
                ii = i
                jj = j
                read(unit,rec=ivar+(k-1)*nvar,err=903)
     +                ((field(jx,jy,j,i),jx=1,nx),jy=1,ny)
                if ( endian*localendian.eq.-1 ) then
                    call swapbyte4(field(1,1,j,i),nxf*nyf)
                endif
                if ( undef.ne.3e33 ) then
                    do jy=1,ny
                        do jx=1,nx
                            if ( abs((field(jx,jy,j,i)-undef)/undef)
     +                           .lt.1e-6 ) then
                                field(jx,jy,j,i) = 3e33
                            elseif ( field(jx,jy,j,i).gt.1e19 ) then
                                ierr = ierr + 1
                                if ( ierr.ge.nerr ) then
                                    nerr = 2*nerr
                                    write(0,*) 'readdatfile: error: '
     +                                   ,ierr,' field(',jx,jy,j,i
     +                                   ,') > 1e19: ',field(jx,jy,j,i)
     +                                   ,undef
                                end if
                                field(jx,jy,j,i) = 3e33
                            endif
                        enddo
                    enddo
                endif
                if ( lwrite.and.notfirst ) then
                    do jy=1,ny
                        do jx=1,nx
                            if ( notfirst.and. field(jx,jy,j,i).lt.1e33
     +                            ) then
                                notfirst = .FALSE.
                                print *,'readdatfile: first valid point'
     +                                ,jx,jy,j,i,field(jx,jy,j,i)
                            endif
                        enddo
                    enddo
                endif
            endif
        enddo
        if ( deletefile ) then
            close(unit)
***         gives error message with Portland
***         close(unit,status='delete')
***         I now rm it right after opening
        else
            close(unit,status='keep')
        endif
*
*       increase the number of time steps to account for all the Feb 29
*       in non-leap years
*
        nt = nt + noleap
*
*       for completeness, these loops will almost never be executed
*
        k = jj + 1
        i = ii
        if ( i.le.yrend .and. k.le.nperyear ) then
            if ( lwrite ) print *,'readdatfile: making rest of year '//
     +           'absent ',k,nperyear,i
            do j=k,nperyear
                do jy=1,ny
                    do jx=1,nx
                        field(jx,jy,j,i) = 3e33
                    enddo
                enddo
            enddo
        endif
        if ( yr2.ge.i+1 .and. lwrite ) print *
     +       ,'readdatfile: making absent rest of years absent ',i+1,yr2
        do k=i+1,yr2
            do j=1,nperyear
                do jy=1,ny
                    do jx=1,nx
                        field(jx,jy,j,k) = 3e33
                    enddo
                enddo
            enddo
        enddo
        if ( lwrite .and. notfirst ) then
            print *,'Found no valid data in file'
            call abort
        endif
        return
  903   write(0,*) 'readdatfile: error reading from datafile '
     +        ,trim(datfile)
        call abort
  904   write(0,*) 'readdatfile: error cannot locate datafile '
     +        ,trim(datfile)
        call abort
        end
