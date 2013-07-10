*  #[ getperyear:
        subroutine getperyear(ncid,varid,tt,nt,firstmo,firstyr,nperyear
     +       ,iperyear,ltime,tdefined,ntmax,lwrite)
*
*       convert the time axis to increment 1 starting from 0
*       so that it can be described by firstmo,firstyr and nperyear
*
        implicit none
#include "netcdf.inc"
        integer ncid,varid,nt,firstmo,firstyr,nperyear,ntmax
        real*8 tt(ntmax)
        character ltime*(*)
        logical tdefined(ntmax),lwrite
        integer i,j,k,n,it,status,iperyear,firstdy,dpm(12,2)
        real*8 dtt,tt0,dtt1,dtt2
        character units*(nf_max_name),timeorigin*(nf_max_name),
     +       calendar*(nf_max_name)
        logical all_leap,lwarn
        integer julday
        data dpm /31,28,31,30,31,30,31,31,30,31,30,31,
     +            31,29,31,30,31,30,31,31,30,31,30,31/

*       get long_name attribute
        call gettextattopt(ncid,varid,'long_name',ltime,lwrite)
        if ( ltime.eq.' ' ) ltime = 'time'
*       get units attribute
        call gettextatt(ncid,varid,'units',units,lwrite)
        call tolower(units)
        if ( nt.gt.1 ) then
!           choose second-smallest positive increment
            dtt1 = 3e33
            dtt2 = 3e33
            do it=2,nt
                dtt = tt(it) - tt(it-1)
                if ( dtt.gt.0 ) then
                    if ( dtt.lt.dtt1 ) then
                        dtt2 = dtt1
                        dtt1 = dtt
                    else if ( dtt.lt.dtt2 ) then
                        dtt2 = dtt
                    end if
                end if
            end do
            if ( nt.gt.2 ) then
                dtt = dtt2
            else
                dtt = dtt1
            end if
            if ( dtt.le.0 ) then
                write(0,*) 'getperyear: error: time step is zero'
                if ( lwrite ) then
                    print *,1,tt(1)
                    do it=2,nt
                        print *,it,tt(it),tt(it)-tt(it-1)
                    end do
                end if
                call abort
            end if
            if ( lwrite ) print *,'dtt = ',dtt
        else
            dtt = 1
        endif
        i = index(units,'since')
        if ( i.eq.0 ) then
            call gettextatt(ncid,varid,'time_origin',timeorigin,lwrite)
            units(len_trim(units)+1:) = ' since '//timeorigin(1:i)
            n = n + 7 + i
        endif
        i = index(units,'years since')
        j = index(units,'year since')
        all_leap = .false.
        if ( i.eq.0 .and.j.eq.0 ) then
            i = index(units,'months since')
            k = index(units,'month since')
            j = index(units,'mon_julian since')
            if ( i.eq.0 .and. j.eq.0 .and. k.eq.0 ) then
                i = index(units,'days since')
                if ( i.eq.0 ) then
                    i = index(units,'hours since')
                    if ( i.eq.0 ) then
                        i = index(units,'minutes since')
                        if ( i.eq.0 ) then
                            i = index(units,'seconds since')
                            if ( i.eq.0 ) then
                                write(0,*) 'getperyear: cannot handle '
     +                               ,'unit ''',trim(units),''', only '
     +                               ,'''days/months/hours/minutes/'
     +                               ,'seconds since'''
                                call abort
                            else
                                i = i+13
                                iperyear = nint(365.24*24*3600)
                                nperyear = nint(365.24d0*24*3600/dtt)
                            end if
                        else
                            i = i+13
                            iperyear = nint(365.24*24*60)
                            nperyear = nint(365.24d0*24*60/dtt)
                        end if
                    else
                        i = i+11
                        iperyear = nint(365.24*24)
                        nperyear = nint(365.24d0*24/dtt)
                    endif
                else
                    i = i+10
                    iperyear = 366
                    nperyear = nint(366/dtt)
*                   needed for first month of NCAR CCSM data...
                    if ( nperyear.gt.12 .and. nperyear.le.17 )
     +                   nperyear = 12
                endif
*               calendar type
                call gettextattopt(ncid,varid,'calendar',calendar,lwrite
     +               )
                if ( calendar.eq.' ' ) then
                    iperyear = 366*nint(iperyear/366.) ! gregorian
                else
                    call tolower(calendar)
                    if ( calendar(1:9).eq.'gregorian' .or.
     +                   calendar(1:19).eq.'proleptic_gregorian'
     +                   .or.calendar(1:8).eq.'standard' 
     +                   .or.calendar(1:8).eq.'all_leap' ) then
                        iperyear = 366*nint(iperyear/366.)
                        if ( calendar(1:8).eq.'all_leap' )
     +                       all_leap = .true.
                    elseif ( calendar(1:6).eq.'noleap' .or. 
     +                       calendar(1:7).eq.'365_day' ) then
                        iperyear = 365*nint(iperyear/365.)
                    elseif ( calendar(1:7).eq.'360_day' ) then
                        iperyear = 360*nint(iperyear/360.)
                    else
                        write(0,*)'getperyear: error: unknown calendar '
     +                       ,trim(calendar)
                        write(*,*)'getperyear: error: unknown calendar '
     +                       ,trim(calendar)
                        call abort
                    endif
                endif
*               this assumes nperyear < iperyear
                if ( lwrite ) print *,'iperyear,nperyear = '
     +               ,iperyear,nperyear
                if ( nperyear.gt.0 .and. iperyear.gt.0 .and. 
     +               nperyear.le.iperyear ) then
                    nperyear = nint(real(iperyear)/
     +                   nint(real(iperyear)/nperyear))
                end if
            else
                if ( j.eq.0 .and. k.eq.0 ) then
                    i = i+12
                elseif ( j.eq.0 ) then
                    i = i+11
                else
                    i = j+16
                endif
                iperyear = 12
                nperyear = nint(12/dtt)
            endif
            if ( nperyear.eq.0 ) then
                write(0,*) 'getperyear: error: cannot handle frequency '
     +               //'lower than once per year'
                call abort
            end if
            if ( nperyear.gt.12 .and. nperyear.le.17 ) then 
!               this occurs when Feb is the first month or when
!               there is little data in a month
                nperyear = 12
                if ( lwrite ) print *,'adjusted nperyear to 12'
            endif
            if ( nperyear.gt.24 .and. nperyear.le.26 ) then
                nperyear = 24
                if ( lwrite ) print *,'adjusted nperyear to 24'
            endif
        else
            if ( i.ne.0 ) then
                i = i+11
            else
                i = j+10
            end if
            iperyear = 1
            nperyear = 1
            if ( abs(dtt-1/12.).lt.0.001 ) then
                nperyear = 12
            else if ( abs(dtt-1/4.).lt.0.001 ) then
                nperyear = 4
            else if ( abs(dtt-1).lt.0.001 ) then
                nperyear = 1
            else
                write(0,*) 'getperyear: error: cannot handle '//
     +               'timestep of ',dtt,' years yet'
                call abort
            endif
        endif
        if ( lwrite ) print *,'getperyear: set nperyear to ',nperyear
 110    continue
        i = i+1
        if ( units(i:i).eq.' ' ) goto 110
        j = i
 111    continue
        j = j+1
        if ( ichar(units(j:j)).ge.ichar('0') .and. 
     +       ichar(units(j:j)).le.ichar('9') ) goto 111
        read(units(i:j-1),'(i4)') firstyr
        if ( lwrite ) print *,'read firstyr=',firstyr
        i = j+1
        j = i
 112    continue
        j = j+1
        if ( ichar(units(j:j)).ge.ichar('0') .and. 
     +       ichar(units(j:j)).le.ichar('9') ) goto 112
        read(units(i:j-1),'(i2)') firstmo
        if ( lwrite ) print *,'read firstmo=',firstmo
        i = j+1
        j = i
 113    continue
        j = j+1
        if ( ichar(units(j:j)).ge.ichar('0') .and. 
     +       ichar(units(j:j)).le.ichar('9') ) goto 113
        read(units(i:j-1),'(i2)') firstdy
        if ( lwrite ) print *,'read firstdy=',firstdy
*       NCEP/NCAR reanalysis have hours since 1-1-1; 
*       renormalize to avoid round-off error and calendar
*       problems (should be using udunits...)
*       NCEP ocean reanalysis also have hours since 1-1-1; 
        if ( .false. .and. iperyear.eq.24*366 ) then
            if ( firstyr.eq.1 .and. firstmo.eq.1 ) then
                if ( abs(tt(1)-16208052).lt.abs(tt(1)-17067072)
     +               ) then
                    firstyr = 1850
                    firstmo = 1
                    do i=1,nt
                        tt(i) = tt(i) - 16208052 ! from ncfile...
                    enddo
                elseif ( abs(tt(1)-17067072).lt.
     +                   abs(tt(1)-17298624) ) then
                    firstyr = 1948
                    firstmo = 1
                    do i=1,nt
                        tt(i) = tt(i) - 17067072 ! from ncfile...
                    enddo
                else
                    firstyr = 1974
                    firstmo = 6
                    do i=1,nt
                        tt(i) = tt(i) - 17298624 ! from ncfile...
                    enddo
                endif
            elseif ( firstyr.eq.1900 .and. firstmo.eq.1 ) then
*               ERA40 data have hours since 1900
                firstyr = 1957
                firstmo = 1
                do i=1,nt
                    tt(i) = tt(i) - 505488 + 243*24
*                   from ncfile plus 8 months
                enddo                            
            elseif ( firstyr.eq.1800 .and. firstmo.eq.1 ) then
*               R2 has hours since 1800
                firstyr = 2008
                firstmo = 1
                do i=1,nt
                    tt(i) = tt(i) - 1823280
*                   from ncfile
                enddo                            
            endif
        endif
*       renormalize to avoid round-off error and calendar
*       problems (should be using udunits...)
        if ( nint(iperyear/366.)*366.eq.iperyear .and. .not.all_leap )
     +       then
            j = julday(firstmo,firstdy,firstyr)
            tt0 = tt(1)/nint(iperyear/366.)
            j = j + tt0
            do i=1,nt
                tt(i) = tt(i) - tt0*nint(iperyear/366.)
            enddo
            call caldat(j,firstmo,firstdy,firstyr)
 120        continue
            if ( lwrite ) then
                print *,'transformed to'
                print *,'firstyr,firstmo,firstdy = ',firstyr,firstmo
     +               ,firstdy
                print *,'tt(1-5) = ',(tt(i),i=1,5)
            endif
            if ( firstmo.gt.1 ) then
                firstmo = firstmo - 1
                firstdy = firstdy + dpm(firstmo,2)*nint(iperyear/366.)
                goto 120
            endif
*           ancient nomenclature
            if ( nperyear.eq.1 ) then
                firstmo = 1
            else if ( nperyear.lt.12 ) then
                firstmo = nint(0.5+real(firstdy)/iperyear*12)
            else if ( nperyear.lt.360 ) then
*               12: this makes both 1-1 and 16-1 go to Jan...
                firstmo = nint(0.5+real(firstdy)/iperyear*nperyear)
            else if ( nperyear.gt.366 ) then
                firstmo = 1 + (firstdy-1)*nint(iperyear/366.) + tt(1)
            else
                firstmo = firstdy
            endif
            if ( lwrite ) print *,'finally, firstmo,firstyr = ',
     +           firstmo,firstyr,nperyear
            if ( firstmo.gt.nperyear ) then
                firstmo = firstmo - nperyear
                firstyr = firstyr + 1
                if ( lwrite ) print *,'finallier, firstmo,firstyr = ',
     +               firstmo,firstyr
            end if
        else
            i = nint(tt(1)/iperyear*nperyear-0.25)
            if ( lwrite ) print *,'adding offset ',i,' periods of '
     +           ,nperyear
            firstyr = firstyr + i/nperyear
            firstmo = 1 + (firstmo-1)*nperyear/12 + mod(i,nperyear)
            if ( lwrite ) print *,'getperyear: firstyr,firstmo = '
     +                   ,firstyr,firstmo
            if ( firstmo.lt.1 ) then
                firstyr = firstyr + firstmo/nperyear - 1
                firstmo = firstmo - nperyear*(firstmo/nperyear - 1)
            elseif ( firstmo.gt.nperyear ) then
                firstyr = firstyr + (firstmo-1)/nperyear
                firstmo = firstmo - nperyear*((firstmo-1)/nperyear)
            endif
        endif
        if ( lwrite ) print *,'getperyear: firstyr,firstmo = ',firstyr
     +       ,firstmo
*
        if ( nperyear.le.360 .eqv. iperyear.le.360 ) then
            dtt = iperyear/real(nperyear)
        else
            dtt = 365.24d0*nint(iperyear/365.)/real(nperyear)
        end if
        if ( lwrite ) print *,'adjusted dtt to ',dtt
        it = 1
        tdefined(it) = .true.
        lwarn = .false.
        do i=2,nt
            if ( abs(tt(i)-tt(i-1)-dtt).lt.0.11*dtt ) then
                it = it + 1
                tdefined(it) = .true.
            else
                if ( tt(i)-tt(i-1).eq.0 ) then
                    write(0,*) 'error: time step zero at step ',i-1,i
                    call abort
                end if
                if ( tt(i)-tt(i-1).lt.0 ) then
                    write(0,*) 'error: time step negative at step ',i-1
     +                   ,i,tt(i-1),tt(i)
                    if ( nperyear.ge.12 ) then
                        write(0,*) '       this corresponds to year ',
     +                       firstyr + real(firstmo)/nperyear
     +                       + tt(i-1)/iperyear,' ~ ',firstyr
     +                       + real(firstmo)/nperyear
     +                       + real(i-1)/nperyear
                    end if
                    call abort
                end if
                if ( .not.lwarn ) then
                    write(0,*) 'warning: irregular time axis'
                    lwarn = .true.
                end if
                do j=1,ntmax-i
                    it = it + 1
                    if ( tt(i)-tt(i-1).lt.(j+0.5)*dtt ) then
                        if ( lwrite ) then
                            print *,'getperyear: at step ',i,
     +                           ' an irregular gap was found between ',
     +                           tt(i-1),tt(i)
                            print *,'            expected ',dtt,
     +                           ' found ',tt(i)-tt(i-1)
                            print *,'            marked ',j-1,' time ',
     +                           'steps as undefined'
                        end if
                        tdefined(it) = .true.
                        exit
                    else
                        tdefined(it) = .false.
                    end if
                end do
            endif
        enddo
        if ( it.ne.nt ) then
            nt = it
        end if
!
        end
*  #] getperyear:
*  #[ addonevariable:
        subroutine addonevariable(ncid,varid,name,ntvars,nvarmax,ndimvar
     +       ,dimids,ix,iy,iz,it,ie,vars,ivars,lvars,svars,units
     +       ,cell_methods,undef,lwrite)
        implicit none
#include "netcdf.inc"
        integer ncid,varid,ntvars,nvarmax,ndimvar,dimids(ndimvar)
     +       ,ix,iy,iz,it,ie,ivars(6,nvarmax)
        real undef
        character name*(*),vars(nvarmax)*(*),lvars(nvarmax)*(*),
     +       svars(nvarmax)*(*),units(nvarmax)*(*),
     +       cell_methods(nvarmax)*(*)
        logical lwrite
        integer i,j,n,status
        real x
        character string*500

        if ( ntvars.ge.nvarmax ) then
            write(0,*) 'addonevariable: found more than ',nvarmax
     +           ,' time-varying variables:'
            do i=1,ntvars
                write(0,*) i,trim(vars(i))
            enddo
            write(0,*) i,trim(name)
            call abort
        endif
        ntvars = ntvars + 1
        vars(ntvars) = name
        call checkstring(vars(ntvars))
        ivars(1,ntvars) = varid
        do i=2,6
            ivars(i,ntvars) = 0
        enddo
        if ( lwrite ) print *,'addonevariable ',varid,trim(name)
        do i=1,ndimvar
            if ( dimids(i).eq.it ) then
                ivars(5,ntvars) = i
                if ( lwrite ) print *,'         dimension ',i,' is time'
            elseif ( dimids(i).eq.ix ) then
                if ( ivars(2,ntvars).eq.0 ) then
                    ivars(2,ntvars) = i
                else
                    goto 801
                endif
                if ( lwrite ) print *,'         dimension ',i
     +               ,' is longitude'
            elseif ( dimids(i).eq.iy ) then
                if ( ivars(3,ntvars).eq.0 ) then
                    ivars(3,ntvars) = i
                else
                    goto 801
                endif
                if ( lwrite ) print *,'         dimension '
     +               ,i,' is latitude'
            elseif ( dimids(i).eq.iz ) then
                if ( ivars(4,ntvars).eq.0 ) then
                    ivars(4,ntvars) = i
                else
                    if ( lwrite ) print *,' ivars(4,ntvars) != 0: ',
     +                   ivars(4,ntvars)
                    goto 801
                endif
                if ( lwrite ) print *,'         dimension '
     +               ,i,' is level'
            elseif ( dimids(i).eq.ie ) then
                if ( ivars(6,ntvars).eq.0 ) then
                    ivars(6,ntvars) = i
                else
                    goto 801
                endif
                if ( lwrite ) print *,'         dimension '
     +               ,i,' is ensemble'
            else
                if ( lwrite ) print *,'dimids(',i,') = ',dimids(i),
     +               ' != ',ix,iy,iz,it,ie
                goto 801
            endif
        enddo
        if ( ntvars.gt.1 ) then
*           check that the variable is on the same grid
            if ( ivars(2,ntvars).gt.0 ) then
                if ( dimids(ivars(2,ntvars)).ne.ix ) goto 902
            endif
            if ( ivars(3,ntvars).gt.0 ) then
                if ( dimids(ivars(3,ntvars)).ne.iy ) goto 902
            endif
            if ( ivars(4,ntvars).gt.0 ) then
                if ( dimids(ivars(4,ntvars)).ne.iz ) goto 902
            endif
            if ( ivars(5,ntvars).gt.0 ) then
                if ( dimids(ivars(5,ntvars)).ne.it ) goto 902
            endif
        endif                   ! first variable?
*       get long name of variable
        call gettextattopt(ncid,varid,'long_name',string,lwrite)
        if ( string.eq.' ' )then
            if ( lwrite ) print *,'addonevariable: warning: variable '
     +           ,trim(name),' does not have long_name'
            lvars(ntvars) = name
        else
            lvars(ntvars) = string
        endif
        call checkstring(lvars(ntvars))
*       get standard name of variable
        call gettextattopt(ncid,varid,'standard_name',svars(ntvars)
     +       ,lwrite)
        if ( svars(ntvars).eq.' ' ) then
            if ( lwrite ) print *,'addonevariable: warning: variable '
     +           ,trim(name),' does not have standard_name'
            svars(ntvars) = name
            call checkstring(svars(ntvars))
        endif
*       get units
        call  gettextattopt(ncid,varid,'units',units(ntvars),lwrite)
        if ( units(ntvars).eq.' ' ) then
            if ( lwrite ) print *
     +           ,'addonevariable: cannot find units attribute'
            i = index(lvars(ntvars),'[')
            j = index(lvars(ntvars),']')
            if ( i.ne.0 .and.j.gt.i+1 ) then
                units(ntvars) = lvars(ntvars)(i+1:j-1)
            endif
        endif
        call checkstring(units(ntvars))
*       get cell_methods
        call gettextattopt(ncid,varid,'cell_methods'
     +       ,cell_methods(ntvars),lwrite)
        call checkstring(cell_methods(ntvars))
*       get fill value
        call getrealattopt(ncid,varid,'_FillValue',undef,lwrite)
*       get missing value
        call getrealattopt(ncid,varid,'missing_value',x,lwrite)
        if ( undef.eq.3e33 .and. x.ne.3e33 ) undef = x
        if ( lwrite ) print *,'addonevariable: found undef = ',undef
*       get max range
*       call getrealattopt(ncid,varid,'valid_max',x,lwrite)
*       undef = min(undef,2.01*x) ! I use undef/2 as cut-off
        return
* 
*       errors
*
 801    if ( lwrite ) write(0,*) 'addonevariable: warning: '//
     +       'disregarding variable with strange dimensions ',trim(name)
        ntvars = ntvars - 1
        return
 902    write(0,*) 'addonevariable: error: found variable with '//
     +       'different dimensions ',trim(name)
        call abort
        end
*  #] addonevariable:
*  #[ getdiminfo:
        subroutine getdiminfo(xy,ncid,varid,xx,nx,lwrite)
*
*       get X,Y.Z axis
*
        implicit none
#include "netcdf.inc"
        integer ncid,varid,nx
        real xx(nx)
        character xy*1
        logical lwrite
        integer status,length,i
        character dimname*(nf_max_name)
*
        if ( lwrite ) print *,'getdiminfo: found ',xy,' axis'
*       get units attribute
        call gettextatt(ncid,varid,'units',dimname,lwrite)
        call tolower(dimname)
        if ( ( xy.eq.'x'
     +        .and. dimname(1:9).ne.'degrees_e'
     +        .and. dimname(1:8).ne.'degree_e' 
     +        ) .or. ( xy.eq.'y'
     +        .and. dimname(1:9).ne.'degrees_n'
     +        .and. dimname(1:8).ne.'degree_n' 
     +        ) ) then
            write(0,*) 'getdiminfo: illegal unit "'
     +           ,trim(dimname), '" on the ',xy,' axis'
            if ( dimname(1:3).ne.'deg' ) call abort
        endif
*       get data, we already checked that nx <= nxmax
        status = nf_get_var_real(ncid,varid,xx)
        if ( status.ne.nf_noerr ) call handle_err(status
     +        ,'nf_get_var_real')
        if ( lwrite ) then
            print *,'axis is ',xx(1),xx(min(2,nx)),'...',xx(nx)
        end if
        do i=1,nx
            if ( abs(xx(i)).gt.1e10 ) then
                write(0,*) 'getdiminfo: error: coordinate value is ',
     +               'invalid: ',i,xx(i)
                call abort
            end if
        end do
        return
        end
*  #] getdiminfo:
*  #[ getzdiminfo:
        subroutine getzdiminfo(xy,ncid,varid,xx,nx,lx,lwrite)
*
*       get X,Y.Z axis with full information
*       lx(1); units, lx(2): standard_name, lx(3): positive
*
        implicit none
#include "netcdf.inc"
        integer ncid,varid,nx
        real xx(nx)
        character xy*1,lx(3)*(*)
        logical lwrite
        integer status,length,i
        character dimname*(nf_max_name)
*
        if ( lwrite ) print *,'getzdiminfo: found ',xy,' axis'
*       get units attribute
        call gettextattopt(ncid,varid,'units',dimname,lwrite)
        lx(1) = dimname
        call tolower(dimname)
        if ( ( xy.eq.'x'
     +        .and. dimname(1:9).ne.'degrees_e'
     +        .and. dimname(1:8).ne.'degree_e' 
     +        ) .or. ( xy.eq.'y'
     +        .and. dimname(1:9).ne.'degrees_n'
     +        .and. dimname(1:8).ne.'degree_n' 
     +        ) ) then
            write(0,*) 'getzdiminfo: cannot yet handle unit '
     +            ,dimname(1:length), ' yet'
            call abort
        endif
*       get standard_name attribute
        call gettextattopt(ncid,varid,'standard_name',lx(2),lwrite)
        call tolower(dimname)
        if ( lx(2).ne.' ' ) then
            if ( ( xy.eq.'x' .and. dimname.ne.'longitude' ) .or. 
     +           ( xy.eq.'y' .and. dimname.ne.'latitude' ) ) then
                write(0,*)
     +               'getzdiminfo: cannot yet handle standard_name '
     +               ,dimname(1:length), ' yet'
                call abort
            endif
        endif
*       get positive
        call gettextattopt(ncid,varid,'positive',lx(3),lwrite)
*       get data, we already checked that nx <= nxmax
        status = nf_get_var_real(ncid,varid,xx)
        if ( status.ne.nf_noerr ) call handle_err(status
     +        ,'nf_get_var_real')
        if ( lwrite ) then
            print *,'axis is ',xx(1),xx(min(2,nx)),'...',xx(nx)
        end if
        do i=1,nx
            if ( abs(xx(i)).gt.1e10 ) then
                write(0,*) 'getdiminfo: error: coordinate value is ',
     +               'invalid: ',i,xx(i)
                call abort
            end if
        end do
        return
        end
*  #] getzdiminfo:
*  #[ handle_err:
        subroutine handle_err(status,string)
        implicit none
#include "netcdf.inc"
        integer status
        character*(*) string
        if ( status.ne.nf_noerr ) then
            write(0,*)'netcdf error: ',status,string,nf_strerror(status)
            call abort
        endif
        end
*  #] handle_err:
*  #[ getrealattopt:
        subroutine getrealattopt(ncid,varid,name,value,lwrite)
*       get one real from an optional attribute
        implicit none
#include "netcdf.inc"
        integer ncid,varid
        real value
        character*(*) name
        logical lwrite
        integer n,status
        logical isnan
*
        status = nf_inq_attlen(ncid,varid,name,n)
        if ( status.ne.nf_noerr ) then
            if ( lwrite ) print *,'getrealattopt: warning: cannot find '
     +           ,name
            value = 3e33
        elseif ( n.ne.1 ) then
            write(0,*) 'getrealattopt: error: ',name,' has len>1: ',n
            call abort
        else
            status = nf_get_att_real(ncid,varid,name,value)
            if ( status.ne.nf_noerr ) call handle_err(status,name)
        endif
        if ( isnan(value) ) value = 3e33
*       
        end
*  #] getrealattopt:
*  #[ gettextattopt:
        subroutine gettextatt(ncid,varid,name,value,lwrite)
        implicit none
        integer ncid,varid
        character*(*) name,value
        logical lwrite
        call gettextattall(ncid,varid,name,value,.true.,lwrite)
        end
        subroutine gettextattopt(ncid,varid,name,value,lwrite)
        implicit none
        integer ncid,varid
        character*(*) name,value
        logical lwrite
        call gettextattall(ncid,varid,name,value,.false.,lwrite)
        end
        subroutine gettextattall(ncid,varid,name,value,req,lwrite)
*       get a text string from an (optional) attribute
        implicit none
#include "netcdf.inc"
        integer ncid,varid
        character*(*) name,value
        character string*4095
        logical req,lwrite
        integer i,n,status
*
        value = ' '
        status = nf_inq_attlen(ncid,varid,name,n)
        if ( lwrite ) print *,'length of text string is ',n
        if ( status.ne.nf_noerr ) then
            if ( lwrite ) print *,'gettextattall: warning: cannot find '
     +           ,trim(name)
        elseif ( n.ge.len(value) ) then
            if ( n.lt.len(string)-1 ) then
                status = nf_get_att_text(ncid,varid,name,string)
                string(n+1:) = ' '
                call stripnonprint(string,lwrite)
                if ( status.ne.nf_noerr ) call handle_err(status,name)
                value = string
            else
!               I should allocate a longer string and get as much as I
!               can...
                write(0,*) 'gettextattall: warning: attribute '
     +               ,trim(name),' longer than string ',len(value),n
     +               ,', not read'
            end if
        else if ( n.gt.0 ) then
            status = nf_get_att_text(ncid,varid,name,value)
            if ( status.ne.nf_noerr ) call handle_err(status,name)
            value(n+1:) = ' '
            call stripnonprint(value,lwrite)
        endif
        if ( req .and. value.eq.' ' ) then
            call handle_err(status,name)
        endif
        if ( lwrite ) print *,'attribute ',trim(name),' is ',trim(value)
*       
        end
*  #] gettextattopt:
*  #[ gettitle:
        subroutine gettitle(ncid,title,lwrite)
!
!       retrieve title from netcdf file
!
        implicit none
#include "netcdf.inc"
        integer ncid
        character title*(*)
        logical lwrite
        integer status,r,i,p
        character format*40
        call gettextattopt(ncid,nf_global,'title',title,lwrite)
        if ( title.eq.' ' ) then
            call gettextattopt(ncid,nf_global,'Title',title,lwrite)
        endif
        status = nf_get_att_int(ncid,nf_global,'realization',r)
        if ( lwrite ) print *,'r,status= ',r,status,nf_noerr
        if ( status.eq.nf_noerr ) then
            status = nf_get_att_int(ncid,nf_global,
     +          'initialization_method',i)
            if ( status.ne.nf_noerr ) i = 0
            status = nf_get_att_int(ncid,nf_global,'physics_version',p)
            if ( status.ne.nf_noerr ) p = 0
            if ( r.lt.10 ) then
            	format = '(a,i1,a,i1,a,i1)'
            else
            	format = '(a,i2,a,i1,a,i1)'
            end if
            write(title(min(len(title)-10,len_trim(title)+2):),format)
     +           'r',r,'i',i,'p',p
        end if
        end
*  #] gettitle:
*  #[ getnumbers:
        subroutine getnumbers(ncid,ndims,nvars,ngatts,unlimdimid,lwrite)
*
*       get number of dimension, variables, attributes, unlimited dimension
*
        implicit none
#include "netcdf.inc"        
        integer ncid,ndims,nvars,ngatts,unlimdimid
        logical lwrite
        integer status
        status = nf_inq(ncid,ndims,nvars,ngatts,unlimdimid)
        if ( status.ne.nf_noerr ) call handle_err(status,'nf_inq')
        if ( lwrite ) then
            print *,'getnumbers: found ',ndims,' dimensions'
            print *,'            found ',nvars,' variables'
            print *,'            found ',ngatts,' global attributes'
            if ( unlimdimid.ge.0 ) then
                print *,'            dimension #',unlimdimid
     +               ,' is unlimited'
            endif
        endif
        if ( nvars.lt.2 ) then
            write(0,*) 'getnumbers: error: only found ',nvars
     +           ,' variables, expecting at least 2'
            call abort
        endif
        end
*  #] getnumbers:
*  #[ getdims:
        subroutine getdims(ncid,ndims,ix,nx,nxmax,iy,ny,nymax,iz,nz
     +       ,nzmax,it,nt,ntmax,ie,nens1,nens2,nensmax,lwrite)
!
!       retrieve dimension information from the netcdf file
!       note that we do not follow CF conventions completely:
!       I just look for recognisable names, not name(var)=name(dim)
!
        implicit none
#include "netcdf.inc"
        integer ncid,ndims,ix,nx,nxmax,iy,ny,nymax,iz,nz,nzmax,it,nt
     +       ,ntmax,ie,nens1,nens2,nensmax
        logical lwrite
        integer status,dimid,len,l
        character name*(nf_max_name)
        if ( lwrite ) print *,'getdims: ncid,ndims = ',ncid,ndims
*
*       which dimension is time?  latitude?  longitude?
*
        nx = 1
        ny = 1
        nz = 1
        nt = 1
        nens1 = 0
        nens2 = 0
        ix = -1
        iy = -1
        iz = -1
        it = -1
        ie = -1
        do dimid=1,ndims
            status = nf_inq_dim(ncid,dimid,name,len)
            if ( status.ne.nf_noerr ) call handle_err(status
     +            ,'nf_inq_dim')
            if ( lwrite ) then
                print *,'getdims: dimension:',dimid
                print *,'         name:     ',trim(name)
                print *,'         len:      ',len
            endif
            call tolower(name)
            if ( index(name,'_bnd').ne.0 ) then
                if ( lwrite ) print *
     +               ,'getdims: disregarding boundary '
     +               ,trim(name)
                cycle
            endif
            l = len_trim(name)
            if ( name(1:3).eq.'tim' .or. name.eq.'t' .or. name.eq.'t1'
     +           .or. name(1:2).eq.'t_' ) then
                if ( lwrite ) print *,'getdims: found time'
                nt = len
                it = dimid
                if ( nt.gt.ntmax ) then
                    write(0,*) 'getdims: error: increase ntmax from '
     +                    ,ntmax,' to at least ',nt
                    write(*,*) 'getdims: error: increase ntmax from '
     +                    ,ntmax,' to at least ',nt
                    call abort
                endif
            elseif ( name.eq.'x' .or. name(1:2).eq.'x_' .or. 
     +               name(1:3).eq.'lon'
     +               .and. name(max(1,l-4):l).ne.'edges' .or.
     +               name(1:4).eq.'xdim' ) then
                if ( nx.eq.1 ) then
                    if ( lwrite ) print *,'getdims: found lon'
                    nx = len
                    if ( nx.gt.nxmax ) then
                        write(0,*) 'getdims: too many values on X axis'
     +                        ,nx,nxmax
                        call abort
                    endif
                    ix = dimid
                else
                    write(0,*) 'getdims: ignoring duplicate lon'
                endif
            elseif ( name.eq.'y' .or. name(1:2).eq.'y_' .or. 
     +               name(1:3).eq.'lat'
     +               .and. name(max(1,l-4):l).ne.'edges' .or.
     +               name(1:4).eq.'ydim') then
                if ( ny.eq.1 ) then
                    if ( lwrite ) print *,'getdims: found lat'
                    ny = len
                    if ( ny.gt.nymax ) then
                        write(0,*) 'getdims: too many values on Y axis'
     +                        ,ny,nymax
                        call abort
                    endif
                    iy = dimid
                else
                    write(0,*) 'getdims: ignoring duplicate lat'
                endif
            elseif ( ( name(1:3).eq.'lev' .or. name.eq.'z' .or.
     +               name.eq.'z1' .or. name(1:3).eq.'z1_' .or.
     +               name(1:3).eq.'sig' .or. name(1:4).eq.'plev' .or.
     +               name(1:3).eq.'dep' .or. name(1:4).eq.'heig' .or.
     +               index(name,'press').ne.0 ).and. name(max(1,l-4):l)
     +               .ne.'edges' .or.name.eq.'unspecified' ) then
                if ( nz.eq.1 ) then
                    if ( lwrite ) print *,'getdims: found lev'
                    nz = len
                    if ( nz.gt.nzmax ) then
                        write(0,*) 'getdims: too many values on Z axis'
     +                        ,nz,nzmax
                        call abort
                    endif
                    iz = dimid
                else
                    write(0,*) 'getdims: ignoring duplicate lev'
                endif
            elseif ( name(1:3).eq.'ens' ) then
                if ( ie.eq.-1 ) then
                    if ( lwrite ) print *,'getdims: found ensemble'
                    nens1 = 0
                    nens2 = len-1
                    ie = dimid
                    if ( nens2.gt.nensmax ) then
                        write(0,*)
     +                       'getdims: too many values on ensemble axis'
     +                       ,nens2,nensmax
                        call abort
                    endif
                else
                    write(0,*) 'getdims: ignoring duplicate ens'
                endif
            endif
        enddo
        if ( lwrite ) then
            print *,'getdims: nx,ny,nz,nt,nens2 = ',nx,ny,nz,nt,nens2
            print *,'         ix,iy,iz,it,ie    = ',ix,iy,iz,it,ie
        end if
        end
*  #] getdims:
*  #[ makelonreasonable:
        subroutine makelonreasonable(xx,nx)
!
!       make sure lons start in -180:180
!       and is monotonically incraesing (or decreasing)
!
        implicit none
        integer nx
        real xx(nx)
        integer i,j
        real xmin,dx
*       make a reasonable interval...
        xmin = xx(1)
        do i=2,nx
            if ( xx(i).lt.xmin ) xmin = xx(i)
        enddo
        if ( xmin.gt.180 ) then
            j = (xmin+180)/360
            do i=1,nx
                xx(i) = xx(i) - j*360
            enddo
        endif
*       search for 360-degree jumps
        do i=2,nx
            dx = xx(i) - xx(i-1)
            if ( dx.gt.180 ) then
                do j=i,nx
                    xx(j) = xx(j) - 360
                enddo
            elseif ( dx.lt.-180 ) then
                do j=i,nx
                    xx(j) = xx(j) + 360
                enddo
            endif
        enddo
        end
*  #] makelonreasonable:
*  #[ stripnonprint:
        subroutine stripnonprint(string,lwrite)
        character*(*) string
        logical lwrite
        integer i
        character validchars*95
        validchars = 
     +       'abcdefghijklmnopqrstuvwxyz'//
     +       'ABCDEFGHIJKLMNOPQRSTUVWXYZ'//
     +       '0123456789'//
     +       ' !@#$%^&*()_-+={}[]:;"|\<,>.?/~`'//char(39)
!       strip non-ascii characters...
        if ( lwrite ) print *,'before weeding ',trim(string)
        do while ( verify(string,validchars).ne.0 )
            i = verify(string,validchars)
            if ( string(i:i).eq.char(0) ) then ! C-style string termination
                string(i:) = ' '
            else ! other garbage
                string(i:i) = ' '
            end if
        end do
        if ( lwrite ) print *,'after weeding  ',trim(string)
        end
*  #] stripnonprint:
