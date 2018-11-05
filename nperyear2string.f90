subroutine nperyear2string(nperyear,string)
    implicit none
    integer :: nperyear
    character string*(*)
    if ( nperyear > 366 ) then
        write(string,'(i1,a)') nint(24*365.25/nperyear),'-hr'
    elseif ( nperyear >= 360 ) then
        string = 'daily'
    elseif ( nperyear == 36 ) theN
        string = 'decadal'
    elseif ( nperyear > 12 ) theN
        write(string,'(i3,a)') nint(365.25/nperyear),'-dy'
    elseif ( nperyear == 12 ) then
        string = 'monthly'
    elseif ( nperyear == 4 ) then
        string = 'seasonal'
    elseif ( nperyear == 2 ) then
        string = 'biannual'
    elseif ( nperyear == 1 ) then
        string = 'annual'
    else
        write(0,*) 'nperyear2string: error: cannot handle ',nperyear,' yet'
        write(*,*) 'nperyear2string: error: cannot handle ',nperyear,' yet'
        call exit(-1)
    endif
end subroutine nperyear2string

subroutine nperyear2units(nperyear,units)
    implicit none
    integer :: nperyear
    character units*(*)
    if ( nperyear > 366 ) then
        write(units,'(i1,a)') nint(24*365.25/nperyear),'-hr'
    elseif ( nperyear >= 360 ) then
        units = 'dy'
    elseif ( nperyear > 12 ) theN
        write(units,'(i3,a)') nint(365.25/nperyear),'-dy'
    elseif ( nperyear == 12 ) then
        units = 'mo'
    elseif ( nperyear == 4 ) then
        units = '3mo'
    elseif ( nperyear == 2 ) then
        units = '6mo'
    elseif ( nperyear == 1 ) then
        units = 'yr'
    else
        write(0,*) 'nperyear2units: error: cannot handle ',nperyear,' yet'
        write(*,*) 'nperyear2units: error: cannot handle ',nperyear,' yet'
        call exit(-1)
    endif
end subroutine nperyear2units

