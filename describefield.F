        program describefield
*
*       write in HTML some information about the field from the 
*       GrDAS ctl file or the netCDF file
*
        implicit none
#include "netcdf.inc"
        integer nxmax,nymax,nzmax,yrbeg,yrend,nvmax,npermax,ndata
     +        ,nensmax,ntmax
        parameter (nxmax=1894,nymax=1974,nzmax=95,
     +       yrbeg=1700,yrend=2300,nvmax=1,npermax=73,
     +       ndata=npermax*(yrend-yrbeg+1),nensmax=230,ntmax=2000000)
*
        integer i,n,nx,ny,nz,nt,nperyear,firstyr,firstmo,nvars,ivars(2
     +       ,nvmax),jvars(6,nvmax),ncid,endian,status,iens,nens1,
     +       dy1,mo1,dy2,mo2,yr,ntp,ntt,it,iarg,nargs,nt1,firstyr1,
     +       firstmo1
        real xx(nxmax),yy(nymax),zz(nzmax),undef,dx,dy,dz,yg(nymax)
     +       ,wg(nymax),pi
        character infile*255,datfile*255,outfile*255,units(nvmax)*60
     +       ,vars(nvmax)*40,lvars(nvmax)*200,svars(nvmax)*80,title*1024
     +       ,history*10000,months(12)*3,seasons(4)*3,lz(3)*20,ltime*100
     +       ,cell_methods(nvmax)*128,ew(2),ns(2),string*20,format*100
        logical xwrap,xrev,yrev,lwrite,lexist,ensemble,tdefined(ntmax)
        integer iargc,llen,leap
        data months /'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug'
     +        ,'Sep','Oct','Nov','Dec'/
        data seasons/'DJF','MAM','JJA','SON'/
        lwrite = .false.
*
*       process command line
*
        n = iargc()
        if ( n.lt.1 ) then
            write(0,*) 'usage: describefield infile.[ctl|nc] '//
     +           '[file2 ...]'
            stop
        endif
        nargs = iargc()
        call getarg(nargs,string)
        if ( string.eq.'debug' .or. string.eq.'lwrite' ) then
            lwrite = .true.
            nargs = nargs - 1
        end if
        ensemble = .false.
        iens = 0
        call getarg(1,infile)
        if ( index(infile,'%').ne.0 .or. index(infile,'++').ne.0 ) then
*           ensemble
            ensemble = .true.
            nens1 = 0
            do iens=0,nensmax
                outfile = infile
                call filloutens(outfile,iens)
                inquire(file=outfile,exist=lexist)
                if ( iens.gt.0 .and. .not.lexist ) goto 10
                if ( iens.eq.0 .and. .not.lexist ) then
                    nens1 = max(nens1,1)
                endif
            enddo
   10       continue
            iens = iens-1
            write(0,'(a,i4,a,i4,a)') 'Found ensemble members ',nens1
     +           ,' to ',iens,' <br>'
            call filloutens(infile,nens1)
        endif
        nt = 0
        ntt = 0
        firstmo = 9999
        firstyr = 9999
        do iarg = 1,nargs
!       no indent for the time being

        call getarg(iarg,infile)
        if ( ensemble ) call filloutens(infile,nens1)
        if ( lwrite ) print *,'describefield: nf_opening file '
     +        ,infile(1:llen(infile))
        status = nf_open(infile,nf_nowrite,ncid)
        if ( status.ne.nf_noerr ) then
            if ( lwrite ) print *,'error opening as netcdf file',status
            call parsectl(infile,datfile,nxmax,nx,xx,nymax,ny,yy,nzmax,
     +           nz,zz,nt1,nperyear,firstyr1,firstmo1,undef,endian,
     +           title,1,nvars,vars,ivars,lvars,units)
            tdefined(1:nt1) = .true.
            ltime = ' '
            svars = ' '
            cell_methods = ' '
            history = ' '
            lz = ' '
            ncid = -1
        else
            call ensparsenc(infile,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax
     +           ,nz,zz,lz,nt1,nperyear,firstyr1,firstmo1,ltime,tdefined
     +           ,ntmax,nens1,iens,undef,title,history,1,nvars,vars
     +           ,jvars,lvars,svars,units,cell_methods)
            if ( .not. ensemble .and. iens.ne.0 ) then
                ensemble = .true.
                write(0,'(a,i4,a,i4,a)') 'Found ensemble members ',nens1
     +           ,' to ',iens,'<br>'
            endif
        endif
!       we assume that if there are multiple files 
!       all properties are the same except the time information
        if ( 10000*firstyr1 + firstmo1 .lt. 10000*firstyr + firstmo )
     +       then
            firstyr = firstyr1
            firstmo = firstmo1
        end if
        nt = nt + nt1
        do it=1,nt1
            if ( tdefined(it) ) ntt = ntt + 1
        end do

        end do ! loop over files
        write(0,'(2a)') title(1:llen(title)),'<br>'
        call getxyprop(xx,nx,yy,ny,xrev,yrev,xwrap)
*
*       X
*
        if ( nx.eq.1 ) then
            if ( jvars(2,1).ne.0 ) then
                if ( xx(1).lt.0 ) then
                    write(0,'(a,f8.2,a)') 'X at ',-xx(1),'&deg; W'
                else
                    write(0,'(a,f8.2,a)') 'X at ',xx(1),'&deg; E'
                endif
            endif
        elseif ( nx.gt.1 ) then
            dx = xx(2) - xx(1)
            do i=2,nx-1
                if ( abs(xx(i+1)-xx(i)-dx).gt.max(0.001,1e-3*abs(xx(i)))
     +               ) then
                    if ( lwrite ) print *,'found irregular step: ',
     +                   xx(i+1),xx(i),dx,xx(i+1)-xx(i)-dx
                    goto 100
                endif
            enddo
            if ( xwrap ) then
                write(0,'(a,i4,f8.2,a)') 'X axis: whole world in ',nx,dx
     +                ,'&deg; steps, '
            else
                write(0,'(a,i5,f8.2,a)') 'X axis: regular grid with  '
     +                ,nx,dx,'&deg; steps, '
            endif
            if ( xx(1).ge.0 ) then
                ew(1) = 'E'
            else
                ew(1) = 'W'
            end if
            if ( xx(nx).ge.0 ) then
                ew(2) = 'E'
            else
                ew(2) = 'W'
            end if
            write(0,'(a,f8.2,3a,f8.2,2a)') 'first point at ',abs(xx(1))
     +           ,'&deg; ',ew(1),', last point at ',abs(xx(nx)),'&deg; '
     +           ,ew(2)
            goto 110
  100       continue
            write(0,'(a,i4,a,10000f8.2,a)') 'X axis: irregular grid of '
     +           ,nx,' points at ',(xx(i),i=1,nx)
  110       continue
        endif
        write(0,'(a)') '<br>'
*
*       Y
*
        if ( ny.eq.1 ) then
            if ( jvars(3,1).ne.0 ) then
                if ( yy(1).lt.0 ) then
                    write(0,'(a,f8.2,a)') 'Y at ',-yy(1),'&deg; S'
                else
                    write(0,'(a,f8.2,a)') 'Y at ',yy(1),'&deg; N'
                endif
            endif
        elseif ( ny.gt.1 ) then
            if ( yy(1).gt.0 ) then
                ns(1) = 'N'
            else
                ns(1) = 'S'
            endif
            if ( yy(ny).gt.0 ) then
                ns(2) = 'N'
            else
                ns(2) = 'S'
            endif
            dy = yy(2) - yy(1)
            do i=2,ny-1
                if ( abs(yy(i+1)-yy(i)-dy).gt.max(0.001,1e-3*abs(yy(i)))
     +               ) goto 200
            enddo
            write(0,'(a,i5,f8.2,a)') 'Y axis: regular grid with  ',ny,dy
     +            ,'&deg; steps, '
            write(0,'(a,f8.2,3a,f8.2,2a)') 'first point at ',abs(yy(1)),
     +           '&deg; ',ns(1),', last point at ',abs(yy(ny)),'&deg; '
     +           ,ns(2)
            goto 210
  200       continue
!           recognise Gaussian grid...
            call legzo(ny,yg,wg)
            pi = 4*atan(1.)
            do i=1,n
                yg(i) = -90 + 180*acos(yg(i))/pi
            end do
            do i=1,n/2
                !!!write(0,*) i,yy(i),yg(i),yy(i)-yg(i)
                if ( abs(abs(yy(i))-abs(yg(i))) .gt.
     +               max(0.01,1e-3*abs(yy(i))) )goto 205
            end do
            write(0,'(a,i5,a)') 'Y axis: Gaussian grid with  ',ny
     +           ,' steps, '
            write(0,'(a,f8.2,3a,f8.2,2a)') 'first point at ',abs(yy(1)),
     +           '&deg; ',ns(1),', last point at ',abs(yy(ny)),'&deg; '
     +           ,ns(2)
            goto 210
 205        continue
            do i=1,n/2
                write(0,*) i,yy(i),yg(i),yy(i)-yg(i)
                if ( abs(abs(yy(i))-abs(yg(i))) .gt.
     +               max(0.1,1e-2*abs(yy(i))) )goto 206
            end do
            write(0,'(a,i5,a)') 'Y axis: almost Gaussian grid with  ',ny
     +           ,' steps, '
            write(0,'(a,f8.2,3a,f8.2,2a)') 'first point at ',abs(yy(1)),
     +           '&deg; ',ns(1),', last point at ',abs(yy(ny)),'&deg; '
     +           ,ns(2)
            goto 210
 206        continue
            write(0,'(a,i4,a,10000f8.2,a)') 'Y axis: irregular grid of '
     +           ,ny,' points at ',(yy(i),i=1,ny)
  210       continue
        endif
        write(0,'(a)') '<br>'
*
*       Z
*
        if ( nz.eq.1 ) then
            if ( jvars(4,1).ne.0 ) then
                write(0,'(a,f8.2,4a)') 'Z at ',zz(1),' ',trim(lz(1)),' '
     +               ,trim(lz(2))
            endif
        elseif ( nz.gt.1 ) then
            dz = zz(2) - zz(1)
            do i=2,nz-1
                if ( abs(zz(i+1)-zz(i)-dz).gt.1e-3*abs(zz(i)) ) goto 300
            enddo
            write(0,'(a,i5,f8.2,a)') 'Z axis: regular grid with  ',nz
     +            ,dz,' steps, '
            write(0,'(a,f8.4,4a)') 'first point at ',zz(1),
     +           ' ',trim(lz(1)),' ',trim(lz(2))
            goto 310
  300       continue
            write(0,'(a,1000f10.2)') 'Z axis: irregular grid ',
     +            (zz(i),i=1,nz)
            write(0,'(4a)') ' ',trim(lz(1)),' ',trim(lz(2))
  310       continue
        endif
        write(0,'(a)') '<br>'
*
*       T
*
        if ( ntt.gt.0 ) then
            if ( nperyear.eq.0 ) then
                write(0,'(a)') 'Time axis was not identified correctly'
            elseif ( nperyear.eq.4 ) then
                write(0,'(2a,i4,2a,i4,a,i4,a)')
     +               'Seasonal data available from ',seasons(firstmo)
     +               ,firstyr,' to ',seasons(1+mod(firstmo+nt-2,4))
     +               ,firstyr+(firstmo+nt-2)/4,' (',ntt,' seasons)'
            elseif ( nperyear.eq.12 ) then
                if ( ntt.lt.10000 ) then
                    format = '(2a,i4,2a,i4,a,i4,a)'
                else
                    format = '(2a,i4,2a,i4,a,i5,a)'
                end if
                write(0,format)
     +               'Monthly data available from ',months(firstmo)
     +               ,firstyr,' to ',months(1+mod(firstmo+nt-2,12))
     +               ,firstyr +(firstmo+nt-2)/12,' (',ntt,' months)'
            elseif ( nperyear.eq.1 ) then
                write(0,'(a,i4,a,i4,a,i4,a)')
     +               'Yearly data available from ',firstyr
     +               ,' to ',firstyr + nt - 1,' (',ntt,' years)'
            elseif ( nperyear.lt.12 ) then
                write(0,'(i2,2a,i4,2a,i4,a,i4,a)') 12/nperyear,
     +               '-monthly data available from ',months(firstmo*12
     +               /nperyear),firstyr,' to ',months(1+mod((firstmo+nt
     +               -2)*12/nperyear,12)),firstyr+(firstmo+nt-2)
     +               /nperyear,' (',ntt,' times)'
            elseif ( nperyear.le.366 ) then
                call getdymo(dy1,mo1,firstmo,nperyear)
                ntp = nt
                if ( nperyear.eq.366 ) then
                    do yr=firstyr,firstyr+(firstmo+nt-2)/nperyear
                        if ( leap(yr).eq.1 .and. (firstmo+nt-2).ge.60)
     +                       then
                            ntp = ntp + 1 ! we count 366 days per year
                        endif
                    enddo
                    if ( firstmo.ge.60 .and. leap(firstyr).eq.1 )
     +                   ntp = ntp - 1
                endif
                call getdymo(dy2,mo2,firstmo+ntp-1,nperyear)
                write(0,'(i2,a,i2.2,a,i4,a,i2.2,a,i4,a,i6,a)') 
     +               nint(366./nperyear),'-daily data available from ',
     +           dy1,months(mo1),firstyr,' to ',dy2,months(mo2),
     +               firstyr+(firstmo+ntp-2)/nperyear,' (',ntp,' times)'
                if ( nperyear.eq.360 ) then
                    write(0,'(a)') '(360dy calendar)'
                elseif ( nperyear.eq.365) then
                    write(0,'(a)') '(365dy calendar)'
                endif
            else
                call getdymo(dy1,mo1,firstmo,nperyear)
                ntp = nt
                if ( nperyear.eq.nint(nperyear/366.)*366 ) then
                    do yr=firstyr,firstyr+(firstmo+nt-2)/nperyear
                        if ( leap(yr).eq.1 ) then
                            ntp = ntp + nint(nperyear/366.) ! we count 366 days per year
                        endif
                    enddo
                    if ( firstmo.ge.60*nint(nperyear/366.)
     +                   .and. leap(firstyr).eq.1 ) then
                        ntp = ntp - nint(nperyear/366.)
                    end if
                endif
                if ( lwrite ) print *,'corrected nt from ',nt,' to ',ntp
                call getdymo(dy2,mo2,firstmo+ntp-1,nperyear)
                write(0,'(i2,a,i2,a,i4,a,i2,a,i4,a,i6,a)') 
     +               nint(24*366./nperyear),'-hourly data available from
     +               ',dy1,months(mo1),firstyr,' to
     +               ',dy2,months(mo2),firstyr+(firstmo+nt-2)/nperyear,'
     +               (',ntt,' times)'
            endif
            call tolower(ltime)
            if ( ltime.ne.'time' .and. ltime.ne.' ' ) then
                write(0,'(3a)') ', time refers to the ',trim(ltime)
     +               ,'<br>'
            end if
        end if
        write(0,'(a)') '<br>'
*
*       vars
*
        do i=1,nvars
            if ( .false. .and. cell_methods(i).ne.' ' ) then
                write(0,'(9a)') 'Variable ',vars(i),'(',
     +               trim(lvars(i)),') in ',trim(units(i)),' {',
     +               trim(cell_methods(i)),'}<br>'
            else
                write(0,'(9a)') 'Variable ',vars(i),'(',
     +               trim(lvars(i)),') in ',trim(units(i)),'<br>'
            endif
        enddo
        end
