subroutine makestandardseries(data,npermax,yrbeg,yrend,nperyear,var,units,lwrite)
!
!       convert a time series in data to standard units:
!
    implicit none
    integer npermax,yrbeg,yrend,nperyear
    real data(npermax,yrbeg:yrend)
    logical lwrite
    character var*(*),units*(*)
    integer yr,mo,ndpm,dpm(12),dpm0(12)
    real offset,slope,mean
    character newunits*60
    integer leap
    save dpm0
    data dpm0 /31,28,31,30,31,30,31,31,30,31,30,31/
    dpm = dpm0
!
    if ( lwrite ) then
        print *,'makestandardseries: npermax,yrbeg,yrend,nperyear = ', &
        &   npermax,yrbeg,yrend,nperyear
        print *,'var = ',trim(var),', units = ',units
    endif
    call getseriesmoment(1,data,npermax,yrbeg,yrend,nperyear,yrbeg,yrend,mean)
    call makestandardunits(mean,nperyear,var,units,newunits,offset,slope,ndpm,lwrite)
    if ( slope.eq.1 .and. offset.eq.0 .and. ndpm.eq.0 ) then
        units = newunits    ! for K => Celsius, hPa => mb, ...
        return
    endif
    if ( lwrite ) then
        print *,'makestandardseries: converting ',trim(var), &
     &           ' from ',trim(units),' to ',trim(newunits)
        print *,'offset = ',offset
        print *,'slope  = ',slope
        print *,'ndpm   = ',ndpm
    endif
    if ( ndpm.ne.0 ) then
        if ( nperyear.eq.360 ) then
            dpm = 30
        endif
    endif
    do yr=yrbeg,yrend
        do mo=1,nperyear
            if ( data(mo,yr).lt.1e33 ) then
                if ( ndpm.ne.0 ) then
                if ( nperyear.eq.366 .and. mo.eq.2 .and. leap(yr).eq.2 ) then
                        data(mo,yr) = data(mo,yr)*29.**ndpm
                    else
                        data(mo,yr) = data(mo,yr)*real(dpm(mo))**ndpm
                    endif
                endif
                data(mo,yr) = data(mo,yr)*slope + offset
            endif
        enddo
    enddo
    units = newunits
end subroutine

subroutine makestandardfield(data,nxf,nyf,nzf,npermax,yrbeg &
     &       ,yrend,nx,ny,nz,nperyear,yr1,yr2,var,units,lwrite)
!
!   convert a field in data to standard units:
!
    implicit none
    integer nxf,nyf,nzf,npermax,yrbeg,yrend,nx,ny,nz,nperyear,yr1,yr2,day,month
    real data(nxf,nyf,nzf,npermax,yrbeg:yrend)
    logical lwrite
    character var*(*),units*(*)
    integer jx,jy,jz,yr,mo,n,ndpm,dpm(12),dpm0(12)
    real offset,slope,mean,factor
    character newunits*60
    integer leap
    save dpm0
    data dpm0 /31,28,31,30,31,30,31,31,30,31,30,31/
    dpm = dpm0
!
    if ( lwrite ) then
        print *,'makestandardfield: nxf,nyf,nzf,npermax,yrbeg,yrend' &
     &           //',nperyear = ',nxf,nyf,nzf,npermax,yrbeg,yrend,nperyear
        print *,'                   nx,ny,nz,nperyear,yr1,yr2 = ',nx &
     &           ,ny,nz,nperyear,yr1,yr2
        print *,'data(',(nx+1)/2,(ny+1)/2,(nz+1)/2,(nperyear+1)/2, &
     &		(yr1+yr2)/2,') = ',data((nx+1)/2,(ny+1)/2,(nz+1)/2, &
     &		(nperyear+1)/2,(yr1+yr2)/2)
        print *,'var = ',trim(var),', units = ',units
    endif
    call estimatemean(data,nxf,nyf,nzf,npermax,yrbeg,yrend, &
     &       nx,ny,nz,nperyear,yr1,yr2,mean,lwrite)
    call makestandardunits(mean,nperyear,var,units,newunits,offset &
     &       ,slope,ndpm,lwrite)
    if ( slope.eq.1 .and. offset.eq.0 .and. ndpm.eq.0 ) then
        units = newunits    ! for K => Celsius, hPa => mb, ...
        return
    endif
    if ( lwrite ) then
        print *,'makestandardfield: converting ',trim(var), &
     &           ' from ',trim(units),' to ',trim(newunits)
        print *,'offset = ',offset
        print *,'slope  = ',slope
        print *,'ndpm   = ',ndpm
    endif
    do yr=yr1,yr2
        do mo=1,nperyear
            if ( ndpm.ne.0 ) then
                if ( nperyear.eq.1 ) then
                    dpm(1) = 365.24/12
                else if ( nperyear.eq.4 ) then
                    if ( leap(yr).eq.1 ) then
                dpm(1) = 90./3
            else
                dpm(1) = 91./3
            end if
            dpm(2) = 92./3
            dpm(3) = 92./3
            dpm(4) = 91./3
            else if ( nperyear.eq.360 ) then
                    dpm = 30
                end if
        call getdymo(day,month,mo,nperyear)
                if ( nperyear.eq.366 .and. month.eq.2 .and. leap(yr).eq.2 ) then
                    factor = slope*29.**ndpm
                else
                    factor = slope*real(dpm(month))**ndpm
                endif
            else
                factor = slope
            endif
            do jz=1,nz
                do jy=1,ny
                    do jx=1,nx
                        if ( data(jx,jy,jz,mo,yr).lt.1e33 ) then
                            data(jx,jy,jz,mo,yr) = offset + data(jx,jy,jz,mo,yr)*factor
                        end if
                    end do
                end do
            end do
        end do
    end do
    units = newunits
    if ( lwrite ) then
        print *,'data(',(nx+1)/2,(ny+1)/2,(nz+1)/2,(nperyear+1)/2, &
     &		(yr1+yr2)/2,') = ',data((nx+1)/2,(ny+1)/2,(nz+1)/2, &
     &		(nperyear+1)/2,(yr1+yr2)/2)
    end if
end subroutine

subroutine makestandardunits(mean,nperyear,invar,units,newunits,offset,slope,ndpm,lwrite)
!
!       convert units to standard units:
!
!       celsius       temperature
!       mm/(interval) precipitation, interval given by nperyear
!       hPa           pressure
!       m/s           wind
!       m             sealevel
!       ??            z500
!
!       newvar = oldvar*slope*dpm^ndpm + offset
!       (dpm=days per month)
!       highly heuristic...
!
    implicit none
    integer nperyear,ndpm
    real offset,slope,mean
    logical lwrite
    character invar*(*),units*(*),newunits*(*)
    integer i
    character var*60,time*10
!
    !!!lwrite = .true.
    if ( lwrite ) then
        print *,'makestandardunits: input'
        print *,'mean,nperyear = ',mean,nperyear
        print *,'invar         = ',trim(invar)
        print *,'units         = ',trim(units)
    endif
    var = invar
    call tolower(var)
    if ( var(1:4) == 'max_' .or. var(1:4) == 'min_' ) then
        var = var(5:)
    end if
    offset = 0
    slope = 1
    ndpm = 0
    newunits = units
!
    if ( units.eq.' ' ) then
!**            write(0,*) 'makestandardunits: empty unit string ',var
        if ( lwrite ) print *,'makestandardunits: empty unit string '
        return
    endif
    if ( units.eq.'1' .or. units.eq.'class' ) then
        slope = 1
        if ( lwrite ) print *,'makestandardunits: dimensionless'
    elseif ( units.eq.'%' .or. units(1:4).eq.'perc' ) then
        slope = 0.01
        newunits = '1'
        if ( lwrite ) print *,'makestandardunits: percentage to fraction: ',slope
!
!       (Wind) speeds
!
    elseif ( (var(1:1).eq.'u' .or. var(1:1).eq.'v' .or. var(1:1).eq.'w' ) .and. &
     &           (units.eq.'m/s' .or. units.eq.'ms-1') ) then
        slope = 1
        newunits = 'm/s'
    elseif ( (var(1:1).eq.'u' .or. var(1:1).eq.'v' .or. var(1:1).eq.'w' ) .and. &
     &           (units.eq.'cm/s' .or. units.eq.'cm s-1') ) then
        slope = 0.01
        newunits = 'm/s'
    elseif ( (var(1:1).eq.'u' .or. var(1:1).eq.'v' .or. var(1:1).eq.'w' ) .and. &
     &           (units.eq.'mm/s' .or. units.eq.'mm s-1') ) then
        slope = 0.001
        newunits = 'm/s'
!
!       Pressure - easy if not wind stress
!
    elseif ( var(1:1).ne.'u' .and. var(1:1).ne.'v' .and. &
     &       var(1:1).ne.'x' .and. var(1:1).ne.'y' .and. &
     &       var(1:4).ne.'taux' .and. var(1:4).ne.'tauy' .and. &
     &       ( units.eq.'hPa' .or. units.eq.'mb' .or. &
     &         units.eq.'millibars' ) ) then
        slope = 1
        newunits = 'mb'
        if ( lwrite ) print *,'makestandardunits: pressure, no conversion'
    elseif ( var(1:1).ne.'u' .and. var(1:1).ne.'v' .and. &
     &       var(1:1).ne.'x' .and. var(1:1).ne.'y' .and. &
     &       var(1:4).ne.'taux' .and. var(1:4).ne.'tauy' .and. &
     &       ( units.eq.'Pa' .or. units.eq.'N/M**2' .or. units.eq.'N m-2' .or. &
     &         units.eq.'N/m2' ) ) then
        slope = 0.01
        newunits = 'mb'
        if ( lwrite ) print *,'makestandardunits: pressure, conversion to mb ',slope
!
!       Geopotential height - relatively easy
!
    elseif ( ( var(1:1).eq.'z' .or. var(1:3).eq.'hgt' ) .and. units.eq.'m' ) then
        slope = 1
        newunits = 'm'
        if ( lwrite ) print *,'makestandardunits: geopotential height, no conversion'
    elseif ( ( var(1:1).eq.'z' .or. var(1:3).eq.'hgt' ) .and. &
     &           ( units.eq.'m**2 s**-2' .or. units.eq.'m2 s-2' ) ) then
        slope = 1/9.81
        newunits = 'm'
        if ( lwrite ) print *,'makestandardunits: pressure, conversion to m ',slope
!
!       Temperature - slightly harder
!
    elseif ( units.eq.'C' .or. units.eq.'Celsius' .or. units.eq.'celsius' .or. &
     &       units(1:8).eq.'degree_C' .or. units(1:8).eq.'degree_c' .or. &
     &       units(1:7).eq.'degreeC' .or. &
     &       units(1:5).eq.'deg_C' .or. units(1:5).eq.'deg_c' .or. &
     &       units(1:4).eq.'degC' .or. units(1:5).eq.'deg C' ) then
        offset = 0
        newunits = 'Celsius'
        if ( lwrite ) print *,'makestandardunits: temperature, no conversion'
        if ( mean.gt.100 ) then
            write(0,*) 'makestandardunits: error: units were ' &
     &               ,trim(units),' but estimated mean is ',mean
        end if
    elseif ( units.eq.'K' .or. units.eq.'Kelvin' .or. units.eq.'kelvin' .or. &
     &       units(1:8).eq.'degree_K' .or. units(1:8).eq.'degree_k' .or. &
     &       units(1:7).eq.'degreeK' .or. units(1:5).eq.'deg_K' .or. &
     &       units(1:5).eq.'deg_k' .or. units(1:4).eq.'degK' ) then
!           it may be an interval
        if ( mean.gt.100 ) then
            offset = -273.15
            if ( lwrite ) print *,'makestandardunits: mean = ',mean, &
 &               'temperature was in Kelvin, add ',offset
        else
            offset = 0
            if ( lwrite ) print *,'makestandardunits: mean = ',mean, &
 &               'temperature interval, no conversion'
        endif
        newunits = 'Celsius'
!
!       Precipitation - the worst
!
    elseif ( var(1:3).eq.'prc' .or. var(1:4).eq.'prec' .or. &
     &       var(1:4).eq.'pre ' .or.  var(1:5).eq.'prate' .or. &
     &       var.eq.'tp' .or. var.eq.'cp'   .or. var.eq.'lsp' .or. &
     &       var.eq.'rr' .or. var.eq.'rh' .or. var.eq.'pr' .or. &
     &       var(1:2).eq.'pp' .or. var.eq.'prlr' .or. &
     &       var(1:3).eq.'pme' .or. var(1:7).eq.'evspsbl' .or. &
     &       var(1:4).eq.'evap' ) then
!       first the weight/height units
        if ( newunits(1:6).eq.'kg m-2' .or. newunits(1:5).eq.'kg/m2' &
     &       .or. newunits(1:5).eq.'KG/M2' &
     &       .or. newunits(1:6).eq.'Kg/m^2' &
     &       .or. newunits(1:6).eq.'kg/m^2' ) then
            if ( lwrite ) print *,'makestandardunits: '// &
     &               'precipitation: converting from kg/m2 to mm'
            slope = 1000    ! density of water [kg/m3], approx
            slope = .001*slope ! m => mm
            newunits = 'mm'//newunits(index(units,'2')+1:)
        endif
        if ( newunits(1:2).eq.'m/' .or. newunits(1:2).eq.'M/' .or. &
     &       newunits(1:2).eq.'m ' .or. newunits(1:2).eq.'M ' ) then
            if ( lwrite ) print *,'makestandardunits: precipitation: converting from m to mm'
            slope = 1000*slope
            newunits = 'mm'//newunits(2:)
        endif
        if ( newunits(1:3).eq.'cm/' .or. newunits(1:3).eq.'CM/'.or. &
     &       newunits(1:3).eq.'cm ' .or. newunits(1:3).eq.'CM ' ) then
            if ( lwrite ) print *,'makestandardunits: precipitation: converting from cm to mm'
            slope = 10*slope
            newunits = 'mm'//newunits(3:)
        endif
        if ( newunits(1:3).eq.'in/' .or. newunits(1:3).eq.'IN/' .or. &
     &       newunits(1:3).eq.'in ' .or. newunits(1:3).eq.'IN ' .or. &
     &       newunits(1:5).eq.'inch/' .or. newunits(1:5).eq.'INCH/' &
     &       .or. newunits(1:5).eq.'inch ' &
     &       .or. newunits(1:5).eq.'INCH ' &
     &       .or. newunits(1:7).eq.'inches/' &
     &       .or. newunits(1:5).eq.'INCHES/' &
     &       .or. newunits(1:7).eq.'inches ' &
     &       .or. newunits(1:5).eq.'INCHES ' ) &
     &           then
            if ( lwrite ) print *,'makestandardunits: precipitation: converting from in to mm'
            slope = 254*slope
            newunits = 'mm'//newunits(3:)
        endif
!           next the time units
        time = newunits(3:)
        if ( time(1:1).eq.' ' ) time = time(2:)
        if ( time(1:1).eq.'/' ) time = time(2:)
        i = index(time,'**') + index(time,'^') + index(time,'-')
        if ( i.ne.0 ) time(i:) = ' '
        if ( lwrite ) print *,'time unit = ',time
!!!            if ( time(1:1).eq.'s' .or. time(1:1).eq.'S' ) then
        if ( .true. .or. nperyear.eq.12 .or. &
     &           nperyear.eq.360 .or. nperyear.eq.365 .or. &
     &           nperyear.eq.366 ) then
!           use mm/day
            if ( time(1:3).eq.'sea' .or. time(1:1).eq.'SEA' ) then
                ndpm = -1
                slope = slope/3
            elseif ( time(1:1).eq.'s' .or. time(1:1).eq.'S' ) then
                slope = slope*60*60*24
            elseif ( time(1:1).eq.'h' .or. time(1:2).eq.'H' ) then
                slope = slope*24
            elseif ( time(1:1).eq.'d' .or. time(1:2).eq.'D' ) then
                slope = slope
            elseif ( time(1:3).eq.'(5d' .or. time(1:3).eq.'(5D' .or. &
     &               time(1:2).eq.'5d' .or. time(1:2).eq.'5D' ) then
                slope = slope/5
            elseif ( time(1:4).eq.'(10d' .or. time(1:4).eq.'(10D' ) then
                slope = slope/10
            elseif ( time(1:1).eq.'m' .or. time(1:2).eq.'M' ) then
                ndpm = -1
            elseif ( time(1:1).eq.'y' .or. time(1:2).eq.'Y' ) then
                slope = slope/365.24
            else
                write(0,*) 'makestandardunits: error: cannot '// &
     &                   ' recognize time unit ',time,' in ',units,var
                write(*,*) 'makestandardunits: error: cannot '// &
     &                   ' recognize time unit ',time,' in ',units,var
                call abort
            endif
            newunits(3:) = '/day'
            if ( lwrite ) print *,'makestandardunits: converted '// &
     &               'precipitation from ',units,' to ',newunits,slope,ndpm
!           does this make sense in mm/dy?
            if ( mean.ne.0 .and. mean < 1e33 .and. &
     &               var(1:3).ne.'pme' .and. var(1:7).ne.'evspsbl' &
     &               .and. ( mean*slope*30.**ndpm.lt.1e-3 .or. &
     &               mean*slope*30.**ndpm.gt.1200 ) ) then
                write(0,*) 'makestandardunits: warning: units were ' &
     &                   ,trim(units),' but this makes about ', &
     &                   mean*slope*30.**ndpm,' mm/day'
                write(0,*) 'mean = ',mean
            end if
        elseif ( nperyear.eq. 12 ) then
!               use mm/month
             if ( time(1:3).eq.'sea' .or. time(1:1).eq.'SEA' ) then
                ndpm = 0
                slope = slope/3
            elseif ( time(1:1).eq.'s' .or. time(1:1).eq.'S' ) then
                slope = slope*60*60*24
                ndpm = 1
            elseif ( time(1:1).eq.'h' .or. time(1:2).eq.'H' ) then
                slope = slope*24
                ndpm = 1
            elseif ( time(1:1).eq.'d' .or. time(1:2).eq.'D' ) then
                ndpm = 1
            elseif ( time(1:1).eq.'m' .or. time(1:2).eq.'M' ) then
                ndpm = 0
            elseif ( time(1:1).eq.'y' .or. time(1:2).eq.'Y' ) then
                slope = slope/365.24
                ndpm = 1
            else
                write(0,*) 'makestandardunits: error: cannot '// &
     &                   ' recognize time unit ',time,' in ',units,var
                write(*,*) 'makestandardunits: error: cannot '// &
     &                   ' recognize time unit ',time,' in ',units,var
                call abort
            endif
            newunits = 'mm/month'
            if ( lwrite ) print *,'makestandardunits: converted '// &
     &               'precipitation from ',units,' to ',newunits,slope,ndpm
        elseif ( nperyear.eq. 1 ) then
!               use mm/year
            if ( time(1:1).eq.'s' .or. time(1:1).eq.'S' ) then
                slope = slope*60*60*24*12
                ndpm = 1
            elseif ( time(1:1).eq.'h' .or. time(1:2).eq.'H' ) then
                slope = slope*24*12
                ndpm = 1
            elseif ( time(1:1).eq.'d' .or. time(1:2).eq.'D' ) then
                slope = slope*12
                ndpm = 1
            elseif ( time(1:1).eq.'m' .or. time(1:2).eq.'M' ) then
                slope = slope*12
            elseif ( time(1:1).eq.'y' .or. time(1:2).eq.'Y' ) then
                slope = slope
            else
                write(0,*) 'makestandardunits: error: cannot '// &
     &                   ' recognize time unit ',time,' in ',units,var
                write(*,*) 'makestandardunits: error: cannot '// &
     &                   ' recognize time unit ',time,' in ',units,var
                call abort
            endif
            newunits = 'mm/month'
            if ( lwrite ) print *,'makestandardunits: converted '// &
     &               'precipitation from ',units,' to ',newunits,slope,ndpm
        endif
!!!            endif
    elseif ( units(1:3).eq.'hr/' ) then
        if ( lwrite ) print *,'makestandardunits: ',units,' seems OK to me'
    else
!
!           nothing fits...
!
!!!         write(0,*) 'makestandardunits: unknown unit string ',units
!!!         write(0,*) '                   nothing done'
        return
    endif
!
end subroutine

subroutine estimatemean(data,nxf,nyf,nzf,npermax,yrbeg,yrend, &
     &       nx,ny,nz,nperyear,yr1,yr2,mean,lwrite)
!
!       make a rough estimate of the mean of the field - principally to 
!       distibguish Celsius from Kelvin.
!
    implicit none
    integer nxf,nyf,nzf,npermax,yrbeg,yrend,nx,ny,nz,nperyear,yr1,yr2
    real data(nxf,nyf,nzf,npermax,yrbeg:yrend),mean
    logical lwrite
    integer n,yr,mo,jx,jy,jz
!   only a few numbers should be enough...
    n = 0
    mean = 0
    do yr=yr1,yr2
        do mo=1,nperyear
            do jz=max(1,nz/2),nz
                do jy=max(1,ny/2),ny ! do not start in Antarctica...
                    do jx=max(1,nx/2),nx
                        if ( data(jx,jy,jz,mo,yr).lt.1e30 .and. &
 &                           data(jx,jy,jz,mo,yr).ne.0 ) then
                            n = n + 1
                            mean = mean + data(jx,jy,jz,mo,yr)
                            if ( n.ge.1000 ) goto 100
                        end if
                    end do
                end do
            end do
        end do
    end do
100 continue
    if ( n.eq.0 ) then
        mean = 3e33
    else
        mean = mean/n
    endif
    if ( lwrite ) print *,'mean = ',mean
end subroutine
!  #] estimatemean:
