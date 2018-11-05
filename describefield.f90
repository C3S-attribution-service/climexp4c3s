program describefield

!   write in HTML some information about the field from the
!   GrDAS ctl file or the netCDF file

    implicit none
    include 'netcdf.inc'
    integer,parameter :: nxmax=3600,nymax=1800,nzmax=95, &
        yrbeg=1700,yrend=2300,nvmax=1,npermax=73, &
        ndata=npermax*(yrend-yrbeg+1),nensmax=750,ntmax=2000000
    integer :: i,n,nx,ny,nz,nt,nperyear,firstyr,firstmo,nvars,ivars(2 &
        ,nvmax),jvars(6,nvmax),ncid,endian,status,iens,nens1,nens2, &
        dy1,mo1,dy2,mo2,yr,ntp,ntt,it,iarg,nargs,nt1,firstyr1, &
        firstmo1
    real :: xx(nxmax),yy(nymax),zz(nzmax),undef,dx,dy,dz,yg(nymax) &
        ,wg(nymax),pi
    character :: infile*255,datfile*255,outfile*255,units(nvmax)*60 &
        ,vars(nvmax)*40,lvars(nvmax)*200,svars(nvmax)*80,title*1024 &
        ,history*50000,months(12)*3,seasons(4)*3,lz(3)*20,ltime*100 &
        ,cell_methods(nvmax)*128,ew(2),ns(2),string*100,format*100  &
        ,metadata(2,100)*2000,halfyears(2)*7
    logical :: xwrap,xrev,yrev,lwrite,lexist,ensemble,tdefined(ntmax)
    integer :: leap
    data months /'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/
    data seasons/'DJF','MAM','JJA','SON'/
    data halfyears /'Oct-Mar','Apr-Sep'/ ! not 'winter','summer' because of the Australians
    lwrite = .false. 

!   process command line

    nargs = command_argument_count()
    if ( nargs < 1 ) then
        write(0,*) 'usage: describefield infile.[ctl|nc] [file2 ...]'
        stop
    endif
    call get_command_argument(nargs,string)
    if ( string == 'debug' .or. string == 'lwrite' ) then
        lwrite = .true. 
        nargs = nargs - 1
    end if
    ensemble = .false. 
    iens = 0
    call get_command_argument(1,infile)
    if ( index(infile,'%') /= 0 .or. index(infile,'++') /= 0 ) then
        ensemble = .true. 
        nens1 = 0
        nens2 = 0
        do iarg=1,nargs
            call get_command_argument(iarg,infile)
            do iens=0,nensmax
                outfile = infile
                call filloutens(outfile,iens)
                inquire(file=trim(outfile),exist=lexist)
                if ( iens > 0 .and. .not. lexist ) goto 10
                if ( iens == 0 .and. .not. lexist ) then
                    nens1 = max(nens1,1)
                endif
                nens2 = max(nens2,iens)
            enddo
         10 continue
        end do
        if ( nens2 >= nens1 ) then
            write(0,'(a,i4,a,i4,a)') 'Found ensemble members ',nens1,' to ',nens2,' <br>'
        else
            ensemble = .false.
        end if
    endif
    nt = 0
    ntt = 0
    firstmo = 9999
    firstyr = 9999
    do iarg = 1,nargs
        call keepalive1('Processing chunk',iarg,nargs)
        call get_command_argument(iarg,infile)
        if ( ensemble ) call filloutens(infile,nens1)
        if ( lwrite ) print *,'describefield: nf_opening file ',trim(infile)
        inquire(file=trim(infile),exist=lexist)
        if ( .not.lexist ) then
            write(0,*) 'describefield: error: cannot find ',trim(infile)
            cycle
        end if
        status = nf_open(infile,nf_nowrite,ncid)
        if ( status /= nf_noerr ) then
            if ( lwrite ) print *,'error opening as netcdf file',status
            call parsectl(infile,datfile,nxmax,nx,xx,nymax,ny,yy,nzmax, &
                nz,zz,nt1,nperyear,firstyr1,firstmo1,undef,endian, &
                title,1,nvars,vars,ivars,lvars,units)
            tdefined(1:nt1) = .true. 
            ltime = ' '
            svars = ' '
            cell_methods = ' '
            history = ' '
            lz = ' '
            ncid = -1
        else
            call ensparsenc(infile,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax &
                ,nz,zz,lz,nt1,nperyear,firstyr1,firstmo1,ltime,tdefined &
                ,ntmax,nens1,iens,undef,title,history,1,nvars,vars &
                ,jvars,lvars,svars,units,cell_methods,metadata)
            if ( .not. ensemble .and. iens /= 0 ) then
                ensemble = .true. 
                write(0,'(a,i4,a,i4,a)') 'Found ensemble members ',nens1 &
                    ,' to ',iens,'<br>'
            endif
        endif
!       we assume that if there are multiple files
!       all properties are the same except the time information
        if ( 10000*firstyr1 + firstmo1 < 10000*firstyr + firstmo ) &
        then
            firstyr = firstyr1
            firstmo = firstmo1
        end if
        nt = nt + nt1
        do it=1,nt1
            if ( tdefined(it) ) ntt = ntt + 1
        end do

    end do ! loop over files
    if ( nt == 0 ) call exit(-1)
    write(0,'(2a)') trim(title),'<br>'
    call getxyprop(xx,nx,yy,ny,xrev,yrev,xwrap)

!   X

    if ( nx == 1 ) then
        if ( jvars(2,1) == 0 ) then
            call getlonfrommetadata(ncid,xx,lwrite)
        end if
        if ( xx(1) < 1e33 ) then
            if ( xx(1) < 0 ) then
                write(0,'(a,f8.2,a)') 'X at ',-xx(1),'&deg; W'
            else
                write(0,'(a,f8.2,a)') 'X at ',xx(1),'&deg; E'
            endif
        endif
    elseif ( nx > 1 ) then
        dx = xx(2) - xx(1)
        do i=2,nx-1
            if ( abs(xx(i+1)-xx(i)-dx) > max(0.001,1e-3*abs(xx(i))) ) then
                if ( lwrite ) print *,'found irregular step: ', &
                    xx(i+1),xx(i),dx,xx(i+1)-xx(i)-dx
                goto 100
            endif
        enddo
        if ( xwrap ) then
            write(0,'(a,i4,f8.2,a)') 'X axis: whole world in ',nx,dx,'&deg; steps, '
        else
            write(0,'(a,i5,f8.2,a)') 'X axis: regular grid with  ',nx,dx,'&deg; steps, '
        endif
        if ( xx(1) >= 0 ) then
            ew(1) = 'E'
        else
            ew(1) = 'W'
        end if
        if ( xx(nx) >= 0 ) then
            ew(2) = 'E'
        else
            ew(2) = 'W'
        end if
        write(0,'(a,f8.2,3a,f8.2,2a)') 'first point at ',abs(xx(1)) &
           ,'&deg; ',ew(1),', last point at ',abs(xx(nx)),'&deg; ',ew(2)
        goto 110
    100 continue
        write(0,'(a,i4,a,10000f8.2,a)') 'X axis: irregular grid of ' &
            ,nx,' points at ',(xx(i),i=1,nx)
    110 continue
    endif
    write(0,'(a)') '<br>'

!   Y

    if ( ny == 1 ) then
        if ( jvars(3,1) == 0 ) then
            call getlatfrommetadata(ncid,yy(1),lwrite)
        end if
        if ( yy(1) < 1e33 ) then
            if ( yy(1) < 0 ) then
                write(0,'(a,f8.2,a)') 'Y at ',-yy(1),'&deg; S'
            else
                write(0,'(a,f8.2,a)') 'Y at ',yy(1),'&deg; N'
            endif
        endif
    elseif ( ny > 1 ) then
        if ( yy(1) > 0 ) then
            ns(1) = 'N'
        else
            ns(1) = 'S'
        endif
        if ( yy(ny) > 0 ) then
            ns(2) = 'N'
        else
            ns(2) = 'S'
        endif
        dy = yy(2) - yy(1)
        do i=2,ny-1
            if ( abs(yy(i+1)-yy(i)-dy) > max(0.001,1e-3*abs(yy(i))) ) goto 200
        enddo
        write(0,'(a,i5,f8.2,a)') 'Y axis: regular grid with  ',ny,dy,'&deg; steps, '
        write(0,'(a,f8.2,3a,f8.2,2a)') 'first point at ',abs(yy(1)), &
            '&deg; ',ns(1),', last point at ',abs(yy(ny)),'&deg; ',ns(2)
        goto 210
    200 continue
!       recognise Gaussian grid...
        call legzo(ny,yg,wg)
        pi = 4*atan(1.)
        do i=1,ny
            yg(i) = -90 + 180*acos(yg(i))/pi
        end do
        do i=1,ny/2
            !!!write(0,*) 'i,yy(i),yg(i) = ',i,yy(i),yg(i),yy(i)-yg(i)
            if ( abs(abs(yy(i))-abs(yg(i))) > max(0.01,1e-3*abs(yy(i))) ) goto 205
        end do
        write(0,'(a,i5,a)') 'Y axis: Gaussian grid with  ',ny,' steps, '
        write(0,'(a,f8.2,3a,f8.2,2a)') 'first point at ',abs(yy(1)), &
            '&deg; ',ns(1),', last point at ',abs(yy(ny)),'&deg; ',ns(2)
        goto 210
    205 continue
        do i=1,ny/2
            !!!write(0,*) i,yy(i),yg(i) = ',i,yy(i),yg(i),yy(i)-yg(i)
            if ( abs(abs(yy(i))-abs(yg(i))) > max(0.1,1e-2*abs(yy(i))) ) goto 206
        end do
        write(0,'(a,i5,a)') 'Y axis: Gaussian grid with  ',ny,' steps, '
        write(0,'(a,f8.2,3a,f8.2,2a)') 'first point at ',abs(yy(1)), &
            '&deg; ',ns(1),', last point at ',abs(yy(ny)),'&deg; ',ns(2)
        goto 210
    206 continue
        write(0,'(a,i4,a,10000f8.2,a)') 'Y axis: irregular grid of ' &
            ,ny,' points at ',(yy(i),i=1,ny)
    210 continue
    endif
    write(0,'(a)') '<br>'

!   Z

    if ( nz == 1 ) then
        if ( jvars(4,1) /= 0 ) then
            write(0,'(a,f8.2,4a)') 'Z at ',zz(1),' ',trim(lz(1)),' ',trim(lz(2))
        endif
    elseif ( nz > 1 ) then
        dz = zz(2) - zz(1)
        do i=2,nz-1
            if ( abs(zz(i+1)-zz(i)-dz) > 1e-3*abs(zz(i)) ) goto 300
        enddo
        write(0,'(a,i5,f8.2,a)') 'Z axis: regular grid with  ',nz,dz,' steps, '
        write(0,'(a,f8.4,4a)') 'first point at ',zz(1),' ',trim(lz(1)),' ',trim(lz(2))
        goto 310
    300 continue
        write(0,'(a,1000f10.2)') 'Z axis: irregular grid ',(zz(i),i=1,nz)
        write(0,'(4a)') ' ',trim(lz(1)),' ',trim(lz(2))
    310 continue
    endif
    write(0,'(a)') '<br>'

!   T

    if ( ntt > 0 ) then
        if ( nperyear == 0 ) then
            write(0,'(a)') 'Time axis was not identified correctly'
        elseif ( nperyear == 1 ) then
            write(0,'(a,i4,a,i4.4,a,i4,a)') &
                'Yearly data available from ',firstyr &
                ,' to ',firstyr + nt - 1,' (',ntt,' years)'
        elseif ( nperyear == 2 ) then
            write(0,'(2a,i4.4,2a,i4.4,a,i4,a)') &
                'Biannual data available from ',halfyears(firstmo) &
                ,firstyr,' to ',halfyears(1+mod(firstmo+nt-2,2)) &
                ,firstyr+(firstmo+nt-2)/2,' (',ntt,' half years)'
        elseif ( nperyear == 4 ) then
            write(0,'(2a,i4.4,2a,i4.4,a,i4,a)') &
                'Seasonal data available from ',seasons(firstmo) &
                ,firstyr,' to ',seasons(1+mod(firstmo+nt-2,4)) &
                ,firstyr+(firstmo+nt-2)/4,' (',ntt,' seasons)'
        elseif ( nperyear == 12 ) then
            if ( ntt < 10000 ) then
                format = '(2a,i4.4,2a,i4.4,a,i4,a)'
            else
                format = '(2a,i4.4,2a,i4.4,a,i5,a)'
            end if
            write(0,format) &
                'Monthly data available from ',months(firstmo) &
                ,firstyr,' to ',months(1+mod(firstmo+nt-2,12)) &
                ,firstyr +(firstmo+nt-2)/12,' (',ntt,' months)'
        elseif ( nperyear < 12 ) then
            write(0,'(i2,2a,i4.4,2a,i4.4,a,i4,a)') 12/nperyear, &
                '-monthly data available from ',months(firstmo*12 &
                /nperyear),firstyr,' to ',months(1+mod((firstmo+nt &
                -2)*12/nperyear,12)),firstyr+(firstmo+nt-2) &
                /nperyear,' (',ntt,' times)'
        elseif ( nperyear == 36 ) then
            write(0,'(a,i2,a,i4.4,a,i2,a,i4.4,a,i4,a)') &
                'Decadal data available from ',5+10*mod(firstmo-1,3), &
                months(1+(firstmo-1)/3),firstyr,' to ',5+10*mod(firstmo+ntt-2,3), &
                months(1+mod((firstmo-1)/3+nt-1,12)) &
                ,firstyr+(firstmo+nt-2)/36,' (',ntt,' decades)'
        elseif ( nperyear <= 366 ) then
            call getdymo(dy1,mo1,firstmo,nperyear)
            ntp = nt
            if ( nperyear == 366 ) then
                do yr=firstyr,firstyr+(firstmo+nt-2)/nperyear
                    if ( leap(yr) == 1 .and. (firstmo+nt-2) >= 60) then
                        ntp = ntp + 1 ! we count 366 days per year
                    endif
                enddo
                if ( firstmo >= 60 .and. leap(firstyr) == 1 ) ntp = ntp - 1
            endif
            call getdymo(dy2,mo2,firstmo+ntp-1,nperyear)
            write(0,'(i2,a,i2.2,a,i4.4,a,i2.2,a,i4.4,a,i6,a)') &
                nint(366./nperyear),'-daily data available from ', &
                dy1,months(mo1),firstyr,' to ',dy2,months(mo2), &
                firstyr+(firstmo+ntp-2)/nperyear,' (',ntp,' times)'
            if ( nperyear == 360 ) then
                write(0,'(a)') '(360dy calendar)'
            elseif ( nperyear == 365) then
                write(0,'(a)') '(365dy calendar)'
            endif
        else
            call getdymo(dy1,mo1,firstmo,nperyear)
            ntp = nt
            if ( nperyear == nint(nperyear/366.)*366 ) then
                do yr=firstyr,firstyr+(firstmo+nt-2)/nperyear
                    if ( leap(yr) == 1 ) then
                        ntp = ntp + nint(nperyear/366.) ! we count 366 days per year
                    endif
                enddo
                if ( firstmo >= 60*nint(nperyear/366.) .and. leap(firstyr) == 1 ) then
                    ntp = ntp - nint(nperyear/366.)
                end if
            endif
            if ( lwrite ) print *,'corrected nt from ',nt,' to ',ntp
            call getdymo(dy2,mo2,firstmo+ntp-1,nperyear)
            write(0,'(i2,a,i2.2,a,i4.4,a,i2.2,a,i4.4,a,i6,a)') &
                nint(24*366./nperyear),'-hourly data available from &
                ',dy1,months(mo1),firstyr,' to &
                ',dy2,months(mo2),firstyr+(firstmo+nt-2)/nperyear,' &
                (',ntt,' times)'
        endif
        call tolower(ltime)
        if ( ltime /= 'time' .and. ltime /= ' ' ) then
            write(0,'(3a)') ', time refers to the ',trim(ltime),'<br>'
        end if
    end if
    write(0,'(a)') '<br>'

!   vars

    do i=1,nvars
        if ( .false. .and. cell_methods(i) /= ' ' ) then
            write(0,'(9a)') 'Variable ',vars(i),'(', &
                trim(lvars(i)),') in ',trim(units(i)),' {', &
                trim(cell_methods(i)),'}<br>'
        else
            write(0,'(9a)') 'Variable ',vars(i),'(', &
                trim(lvars(i)),') in ',trim(units(i)),'<br>'
        endif
    enddo
end program
