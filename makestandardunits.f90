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
            npermax,yrbeg,yrend,nperyear
        print *,'var = ',trim(var),', units = ',units
    end if
    call getseriesmoment(1,data,npermax,yrbeg,yrend,nperyear,yrbeg,yrend,mean)
    call makestandardunits(mean,nperyear,var,units,newunits,offset,slope,ndpm,lwrite)
    if ( slope==1 .and. offset==0 .and. ndpm==0 ) then
        units = newunits    ! for K => Celsius, hPa => mb, ...
        return
    end if
    if ( lwrite ) then
        print *,'makestandardseries: converting ',trim(var), &
                ' from ',trim(units),' to ',trim(newunits)
        print *,'offset = ',offset
        print *,'slope  = ',slope
        print *,'ndpm   = ',ndpm
    end if
    if ( ndpm /= 0 ) then
        if ( nperyear == 360 ) then
            dpm = 30
        end if
    end if
    do yr=yrbeg,yrend
        do mo=1,nperyear
            if ( data(mo,yr) < 1e33 ) then
                if ( ndpm /= 0 ) then
                    if ( nperyear == 366 .and. mo == 2 .and. leap(yr) == 2 ) then
                        data(mo,yr) = data(mo,yr)*29.**ndpm
                    else
                        data(mo,yr) = data(mo,yr)*real(dpm(mo))**ndpm
                    end if
                end if
                data(mo,yr) = data(mo,yr)*slope + offset
            end if
        end do
    end do
    units = newunits
end subroutine makestandardseries

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
                 //',nperyear = ',nxf,nyf,nzf,npermax,yrbeg,yrend,nperyear
        print *,'                   nx,ny,nz,nperyear,yr1,yr2 = ',nx &
                 ,ny,nz,nperyear,yr1,yr2
        print *,'data(',(nx+1)/2,(ny+1)/2,(nz+1)/2,(nperyear+1)/2, &
      		(yr1+yr2)/2,') = ',data((nx+1)/2,(ny+1)/2,(nz+1)/2, &
      		(nperyear+1)/2,(yr1+yr2)/2)
        print *,'var = ',trim(var),', units = ',units
    end if
    call estimatemean(data,nxf,nyf,nzf,npermax,yrbeg,yrend, &
             nx,ny,nz,nperyear,yr1,yr2,mean,lwrite)
    call makestandardunits(mean,nperyear,var,units,newunits,offset &
             ,slope,ndpm,lwrite)
    if ( slope == 1 .and. offset == 0 .and. ndpm == 0 ) then
        units = newunits    ! for K => Celsius, hPa => mb, ...
        return
    end if
    if ( lwrite ) then
        print *,'makestandardfield: converting ',trim(var), &
                ' from ',trim(units),' to ',trim(newunits)
        print *,'offset = ',offset
        print *,'slope  = ',slope
        print *,'ndpm   = ',ndpm
    end if
    do yr=yr1,yr2
        do mo=1,nperyear
            if ( ndpm /= 0 ) then
                if ( nperyear == 1 ) then
                    dpm(1) = 365.24/12
                else if ( nperyear == 4 ) then
                    if ( leap(yr) == 1 ) then
                dpm(1) = 90./3
            else
                dpm(1) = 91./3
            end if
            dpm(2) = 92./3
            dpm(3) = 92./3
            dpm(4) = 91./3
            else if ( nperyear == 360 ) then
                    dpm = 30
                end if
        call getdymo(day,month,mo,nperyear)
                if ( nperyear == 366 .and. month == 2 .and. leap(yr) == 2 ) then
                    factor = slope*29.**ndpm
                else
                    factor = slope*real(dpm(month))**ndpm
                end if
            else
                factor = slope
            end if
            do jz=1,nz
                do jy=1,ny
                    do jx=1,nx
                        if ( data(jx,jy,jz,mo,yr) < 1e33 ) then
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
     	    (yr1+yr2)/2,') = ',data((nx+1)/2,(ny+1)/2,(nz+1)/2, &
     	    (nperyear+1)/2,(yr1+yr2)/2)
    end if
end subroutine makestandardfield

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
    integer i,nday
    character var*60,time*10
    logical,external :: isnumchar
!
    !!!lwrite = .true.
    if ( lwrite ) then
        print *,'makestandardunits: input'
        print *,'mean,nperyear = ',mean,nperyear
        print *,'invar         = ',trim(invar)
        print *,'units         = ',trim(units)
    end if
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
    if ( units == ' ' ) then
!**            write(0,*) 'makestandardunits: empty unit string ',var
        if ( lwrite ) print *,'makestandardunits: empty unit string '
        return
    end if
    if ( units == '1' .or. units == 'class' ) then
        slope = 1
        if ( lwrite ) print *,'makestandardunits: dimensionless'
    else if ( units == '%' .or. units(1:4) == 'perc' ) then
        slope = 0.01
        newunits = '1'
        if ( lwrite ) print *,'makestandardunits: percentage to fraction: ',slope
!
!       (Wind) speeds
!
    else if ( (var(1:1) == 'u' .or. var(1:1) == 'v' .or. var(1:1) == 'w' ) .and. &
     &           (units == 'm/s' .or. units == 'ms-1') ) then
        slope = 1
        newunits = 'm/s'
    else if ( (var(1:1) == 'u' .or. var(1:1) == 'v' .or. var(1:1) == 'w' ) .and. &
     &           (units == 'cm/s' .or. units == 'cm s-1') ) then
        slope = 0.01
        newunits = 'm/s'
    else if ( (var(1:1) == 'u' .or. var(1:1) == 'v' .or. var(1:1) == 'w' ) .and. &
     &           (units == 'mm/s' .or. units == 'mm s-1') ) then
        slope = 0.001
        newunits = 'm/s'
!
!       Pressure - easy if not wind stress
!
    else if ( var(1:1)/='u' .and. var(1:1)/='v' .and. &
     &       var(1:1)/='x' .and. var(1:1)/='y' .and. &
     &       var(1:4)/='taux' .and. var(1:4)/='tauy' .and. &
     &       ( units == 'hPa' .or. units == 'mb' .or. &
     &         units == 'millibars' ) ) then
        slope = 1
        newunits = 'mb'
        if ( lwrite ) print *,'makestandardunits: pressure, no conversion'
    else if ( var(1:1)/='u' .and. var(1:1)/='v' .and. &
     &       var(1:1)/='x' .and. var(1:1)/='y' .and. &
     &       var(1:4)/='taux' .and. var(1:4)/='tauy' .and. &
     &       ( units == 'Pa' .or. units == 'N/M**2' .or. units == 'N m-2' .or. &
     &         units == 'N/m2' ) ) then
        slope = 0.01
        newunits = 'mb'
        if ( lwrite ) print *,'makestandardunits: pressure, conversion to mb ',slope
!
!       Geopotential height - relatively easy
!
    else if ( ( var(1:1) == 'z' .or. var(1:3) == 'hgt' ) .and. units == 'm' ) then
        slope = 1
        newunits = 'm'
        if ( lwrite ) print *,'makestandardunits: geopotential height, no conversion'
    else if ( ( var(1:1) == 'z' .or. var(1:3) == 'hgt' ) .and. &
     &           ( units == 'm**2 s**-2' .or. units == 'm2 s-2' ) ) then
        slope = 1/9.81
        newunits = 'm'
        if ( lwrite ) print *,'makestandardunits: pressure, conversion to m ',slope
!
!       Temperature - slightly harder
!
    else if ( units == 'C' .or. units == 'Celsius' .or. units == 'celsius' .or. &
     &       units == 'degrees Celsius' .or. &
     &       units(1:8) == 'degree_C' .or. units(1:8) == 'degree_c' .or. &
     &       units(1:7) == 'degreeC' .or. units(1:8) == 'degree C' .or. &
     &       units(1:5) == 'deg_C' .or. units(1:5) == 'deg_c' .or. &
     &       units(1:4) == 'degC' .or. units(1:5) == 'deg C' ) then
        offset = 0
        newunits = 'Celsius'
        if ( lwrite ) print *,'makestandardunits: temperature, no conversion'
        if ( mean>100 ) then
            write(0,*) 'makestandardunits: error: units were ' &
     &               ,trim(units),' but estimated mean is ',mean
        end if
    else if ( units == 'K' .or. units == 'Kelvin' .or. units == 'kelvin' .or. &
     &       units(1:8) == 'degree_K' .or. units(1:8) == 'degree_k' .or. &
     &       units(1:7) == 'degreeK' .or. units(1:5) == 'deg_K' .or. &
     &       units(1:5) == 'deg_k' .or. units(1:4) == 'degK' ) then
!           it may be an interval
        if ( mean>100 ) then
            offset = -273.15
            if ( lwrite ) print *,'makestandardunits: mean = ',mean, &
 &               'temperature was in Kelvin, add ',offset
        else
            offset = 0
            if ( lwrite ) print *,'makestandardunits: mean = ',mean, &
 &               'temperature interval, no conversion'
        end if
        newunits = 'Celsius'
!
!       Precipitation - the worst
!
    else if ( var(1:3) == 'prc' .or. var(1:4) == 'prec' .or. &
     &       var(1:4) == 'pre ' .or.  var(1:5) == 'prate' .or. &
     &       var == 'tp' .or. var == 'cp'   .or. var == 'lsp' .or. &
     &       var == 'rr' .or. var == 'rh' .or. var == 'pr' .or. &
     &       var(1:2) == 'pp' .or. var == 'prlr' .or. &
     &       var(1:3) == 'pme' .or. var(1:7) == 'evspsbl' .or. &
     &       var(1:4) == 'evap' ) then
!       first the weight/height units
        if ( newunits(1:6) == 'kg m-2' .or. newunits(1:5) == 'kg/m2' &
     &       .or. newunits(1:5) == 'KG/M2' &
     &       .or. newunits(1:6) == 'Kg/m^2' &
     &       .or. newunits(1:6) == 'kg/m^2' ) then
            if ( lwrite ) print *,'makestandardunits: '// &
     &               'precipitation: converting from kg/m2 to mm'
            slope = 1000    ! density of water [kg/m3], approx
            slope = .001*slope ! m => mm
            newunits = 'mm'//newunits(index(units,'2')+1:)
        end if
        if ( newunits(1:2) == 'm/' .or. newunits(1:2) == 'M/' .or. &
     &       newunits(1:2) == 'm ' .or. newunits(1:2) == 'M ' ) then
            if ( lwrite ) print *,'makestandardunits: precipitation: converting from m to mm'
            slope = 1000*slope
            newunits = 'mm'//newunits(2:)
        end if
        if ( newunits(1:3) == 'cm/' .or. newunits(1:3) == 'CM/'.or. &
     &       newunits(1:3) == 'cm ' .or. newunits(1:3) == 'CM ' ) then
            if ( lwrite ) print *,'makestandardunits: precipitation: converting from cm to mm'
            slope = 10*slope
            newunits = 'mm'//newunits(3:)
        end if
        if ( newunits(1:3) == 'in/' .or. newunits(1:3) == 'IN/' .or. &
     &       newunits(1:3) == 'in ' .or. newunits(1:3) == 'IN ' .or. &
     &       newunits(1:5) == 'inch/' .or. newunits(1:5) == 'INCH/' &
     &       .or. newunits(1:5) == 'inch ' &
     &       .or. newunits(1:5) == 'INCH ' &
     &       .or. newunits(1:7) == 'inches/' &
     &       .or. newunits(1:5) == 'INCHES/' &
     &       .or. newunits(1:7) == 'inches ' &
     &       .or. newunits(1:5) == 'INCHES ' ) &
     &           then
            if ( lwrite ) print *,'makestandardunits: precipitation: converting from in to mm'
            slope = 254*slope
            newunits = 'mm'//newunits(3:)
        end if
!           next the time units
        time = newunits(3:)
        if ( time(1:1) == ' ' ) time = time(2:)
        if ( time(1:1) == '/' ) time = time(2:)
        i = index(time,'**') + index(time,'^') + index(time,'-')
        if ( i/=0 ) time(i:) = ' '
        if ( lwrite ) print *,'time unit = ',time
!!!            if ( time(1:1) == 's' .or. time(1:1) == 'S' ) then
        if ( .true. .or. nperyear == 12 .or. &
     &           nperyear == 360 .or. nperyear == 365 .or. &
     &           nperyear == 366 ) then
!           use mm/day
            if ( isnumchar(time(1:1)) ) then ! 'mm/NNday'
                do i=2,10
                    if ( .not. isnumchar(time(i:i)) ) exit
                end do
                read(time(:i-1),*) nday
                slope = slope/nday
                time = time(i:)
            end if
            if ( time(1:3) == 'sea' .or. time(1:1) == 'SEA' ) then
                ndpm = -1
                slope = slope/3
            else if ( time(1:1) == 's' .or. time(1:1) == 'S' ) then
                slope = slope*60*60*24
            else if ( time(1:1) == 'h' .or. time(1:2) == 'H' ) then
                slope = slope*24
            else if ( time(1:1) == 'd' .or. time(1:2) == 'D' ) then
                slope = slope
            else if ( time(1:3) == '(5d' .or. time(1:3) == '(5D' .or. &
     &               time(1:2) == '5d' .or. time(1:2) == '5D' ) then
                slope = slope/5
            else if ( time(1:4) == '(10d' .or. time(1:4) == '(10D' ) then
                slope = slope/10
            else if ( time(1:1) == 'm' .or. time(1:2) == 'M' ) then
                ndpm = -1
            else if ( time(1:1) == 'y' .or. time(1:2) == 'Y' ) then
                slope = slope/365.24
            else
                write(0,*) 'makestandardunits: error: cannot '// &
                           ' recognise time unit ',time,' in ',units,var
                write(*,*) 'makestandardunits: error: cannot '// &
                           ' recognise time unit ',time,' in ',units,var
                call exit(-1)
            end if
            newunits(3:) = '/day'
            if ( lwrite ) print *,'makestandardunits: converted '// &
                     'precipitation from ',trim(units),' to ',trim(newunits),slope,ndpm
!           does this make sense in mm/dy?
            if ( mean /= 0 .and. mean < 1e33 .and. &
                     var(1:3) /= 'pme' .and. var(1:7) /= 'evspsbl' &
                     .and. ( mean*slope*30.**ndpm < 1e-3 .or. &
                     mean*slope*30.**ndpm > 1200 ) ) then
                write(0,*) 'makestandardunits: warning: units were ' &
                         ,trim(units),' but this makes about ', &
                         mean*slope*30.**ndpm,' mm/day'
                write(0,*) 'mean = ',mean
            end if
        else if ( nperyear ==  12 ) then
!               use mm/month
             if ( time(1:3) == 'sea' .or. time(1:1) == 'SEA' ) then
                ndpm = 0
                slope = slope/3
            else if ( time(1:1) == 's' .or. time(1:1) == 'S' ) then
                slope = slope*60*60*24
                ndpm = 1
            else if ( time(1:1) == 'h' .or. time(1:2) == 'H' ) then
                slope = slope*24
                ndpm = 1
            else if ( time(1:1) == 'd' .or. time(1:2) == 'D' ) then
                ndpm = 1
            else if ( time(1:1) == 'm' .or. time(1:2) == 'M' ) then
                ndpm = 0
            else if ( time(1:1) == 'y' .or. time(1:2) == 'Y' ) then
                slope = slope/365.24
                ndpm = 1
            else
                write(0,*) 'makestandardunits: error: cannot recognize time unit ', &
                    trim(time),' in ',trim(units),' ',trim(var)
                write(*,*) 'makestandardunits: error: cannot recognize time unit ', &
                    trim(time),' in ',trim(units),' ',trim(var)
                call exit(-1)
            end if
            newunits = 'mm/month'
            if ( lwrite ) print *,'makestandardunits: converted '// &
                     'precipitation from ',units,' to ',newunits,slope,ndpm
        else if ( nperyear ==  1 ) then
!               use mm/year
            if ( time(1:1) == 's' .or. time(1:1) == 'S' ) then
                slope = slope*60*60*24*12
                ndpm = 1
            else if ( time(1:1) == 'h' .or. time(1:2) == 'H' ) then
                slope = slope*24*12
                ndpm = 1
            else if ( time(1:1) == 'd' .or. time(1:2) == 'D' ) then
                slope = slope*12
                ndpm = 1
            else if ( time(1:1) == 'm' .or. time(1:2) == 'M' ) then
                slope = slope*12
            else if ( time(1:1) == 'y' .or. time(1:2) == 'Y' ) then
                slope = slope
            else
                write(0,*) 'makestandardunits: error: cannot '// &
                           ' recognize time unit ',time,' in ',units,var
                write(*,*) 'makestandardunits: error: cannot '// &
                           ' recognize time unit ',time,' in ',units,var
                call exit(-1)
            end if
            newunits = 'mm/month'
            if ( lwrite ) print *,'makestandardunits: converted '// &
                     'precipitation from ',units,' to ',newunits,slope,ndpm
        end if
!!!            end if
    else if ( units(1:3) == 'hr/' ) then
        if ( lwrite ) print *,'makestandardunits: ',units,' seems OK to me'
    else
!
!           nothing fits...
!
        if ( lwrite ) then
            print *,'makestandardunits: unknown unit string ',trim(units)
            print *,'                   nothing done'
        end if
        return
    end if
!
end subroutine makestandardunits

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
    if ( lwrite ) then
        print *,'estimatemean: yrbeg,yrend = ',yrbeg,yrend
        print *,'              yr1,yr2     = ',yr1,yr2
    end if
    n = 0
    mean = 0
    do yr=yr1,yr2
        do mo=1,nperyear
            do jz=max(1,nz/2),nz
                do jy=max(1,ny/2),ny ! do not start in Antarctica...
                    do jx=max(1,nx/2),nx
                        if ( data(jx,jy,jz,mo,yr) < 1e30 .and. &
 &                           data(jx,jy,jz,mo,yr) /= 0 ) then
                            n = n + 1
                            mean = mean + data(jx,jy,jz,mo,yr)
                            if ( n>=1000 ) goto 100
                        end if
                    end do
                end do
            end do
        end do
    end do
100 continue
    if ( n == 0 ) then
        mean = 3e33
    else
        mean = mean/n
    end if
    if ( lwrite ) print *,'mean = ',mean
end subroutine estimatemean