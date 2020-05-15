program cumul
!
!   compute cumulative time series
!
    implicit none
    include 'param.inc'
    include 'getopts.inc'
    integer yr,mo,nperyear,mens1,mens,yy,mm,month,lastyr,j1,j2
    real data(npermax,yrbeg:yrend),cum(npermax,yrbeg:yrend),sum
    character file*1024,var*100,units*50,lvar*120,svar*120,history*50000, &
        metadata(2,100)*2000,time_units*20,title*1000
    
    call get_command_argument(1,file)
    if ( file == ' ' ) then
        write(0,*) 'cumul: usage: cumul file [mon m1 sel lsel] [normsd]'
        write(0,*) '       computes cumulative values of series, '// &
            'optionally over months m1-m2, optionally setting negative values to zero'
        call exit(-1)
    end if
    if ( index(file,'%%') /= 0 .or. index(file,'++') /= 0 .or. index(file,'@@') /= 0 ) then
        ensemble = .true.
        write(0,*) 'cumul: cannot handle ensembles yet'
        call exit(-1)
    end if
    call readensseriesmeta(file,data,npermax,yrbeg,yrend,nensmax,nperyear,mens1,mens, &
        var,units,lvar,svar,history,metadata,lstandardunits,lwrite)
    call getopts(2,command_argument_count(),nperyear,yrbeg,yrend,.true.,mens1,mens)
    
    month = m1
    call getj1j2(j1,j2,month,nperyear,.false.)
    lastyr = -9999
    sum = 3e33
    cum = 3e33
    do yy=yrbeg,yrend
        if ( m1 > 0 ) then
            if ( data(j1,yr) < 1e33 ) then
                sum = 0 ! start at zero each year
            else
                sum = 3e33
            end if
        else
            if ( data(j1,yr) < 1e33 ) then
                sum = 0
            end if
        end if
        do mm=j1,j2
            mo = mm
            call normon(mo,yy,yr,nperyear)
            if ( sum < 1e33 .and. data(mo,yr) < 1e33 ) then
                lastyr = max(lastyr,yr)
                sum = sum + data(mo,yr)
                if ( lnormsd ) then
                    if ( sum < 0 ) sum = 0
                end if
                cum(mo,yr) = sum
            else
                cum(mo,yr) = 3e33
            end if
        end do
    end do
    var = 'cum_'//var
    lvar = 'cumulative of '//lvar
    call nperyear2units(nperyear,time_units)
    units = trim(units)//' '//time_units
    call printvar(6,var,units,lvar)
    title = 'cumulative '//var
    if ( m1 > 0 ) then
        write(title,'(2a,i2,a,i2)') trim(title),' over months ',m1,' to ',m2
    end if
    call printmetadata(6,' ',' ',title,history,metadata)
    call printdatfile(6,cum,npermax,nperyear,yrbeg,lastyr)
end program cumul
