subroutine adjustunits(oper,nperyear,nperyearnew,units,lwrite)

!   adjust units to reflect a new time scale

    implicit none
    integer :: nperyear,nperyearnew
    character oper*3,units*(*)
    logical :: lwrite

    if ( oper == 'mea' .or. oper == 'sd ' .or. &
         oper == 'min' .or. oper == 'max' ) then
        units = units
    else if ( oper == 'var' ) then
        if ( index(trim(units),' ') == 0 .and. &
             index(trim(units),'*') == 0 .and. &
             index(trim(units),'/') == 0 .and. &
             index(trim(units),'^') == 0 ) then
            units = trim(units)//'^2'
        else
            units = '('//trim(units)//')^2'
        end if
    else if ( oper == 'sum' .or. oper == 'bel' .or. oper == 'abo' ) then
        if ( lwrite ) print *,'changed units from ',trim(units)
        call units2longer(units,nperyear,nperyearnew)
        if ( lwrite ) print *,'to ',trim(units)
    else if ( oper == 'num' ) then
        call nperyear2units(nperyearnew,units)
        units = '1/'//units
        if ( lwrite ) print *,'changed units to ',trim(units)
    else if ( oper == 'xti' .or. oper == 'nti' .or. &
        oper == 'fti' .or. oper == 'lti' .or. &
        oper == 'con' ) then
        if ( nperyear == 360 .or. nperyear == 365 .or. &
        nperyear == 366 ) then
            units = 'dy'
        else if ( nperyear == 12 ) then
            units = 'mo'
        else if ( nperyear == 4 ) then
            units = 'season'
        else
            units = 'period'
        end if
    else
        write(0,*) 'adjustunits: error: unknown operation ',oper
        call exit(-1)
    end if
end subroutine

subroutine adjustvar(oper,var,lwrite)

!   adjust variable name

    implicit none
    character oper*3,var*(*)
    logical :: lwrite

    if ( oper == 'num' ) then
        var = 'n'//var
    else if ( oper == 'max' .or. oper == 'min' ) then
        var = oper//'_'//var
    else if ( oper == 'nti' ) then
        var = 'time_min_'//var
    else if ( oper == 'xti' ) then
        var = 'time_max_'//var
    else if ( oper == 'fti' ) then
        var = 'time_first_'//var
    else if ( oper == 'lti' ) then
        var = 'time_last_'//var
    end if

end subroutine

subroutine adjustlvar(oper,lvar,nperyear,lwrite)

!   adjust long lvariable name

    implicit none
    integer :: nperyear
    character :: oper*3,lvar*(*)
    logical :: lwrite
    character :: string*100

    call nperyear2string(nperyear,string)
    if ( oper == 'mea' ) then
        lvar = trim(string)//' mean of '//lvar
    else if ( oper == 'sd ' ) then
        lvar = trim(string)//' standard deviation of '//lvar
    else if ( oper == 'sum' ) then
        lvar = trim(string)//' sum of '//lvar
    else if ( oper == 'num' ) then
        lvar = trim(string)//' number of '//lvar
    else if ( oper == 'max' .or. oper == 'min' ) then
        lvar = trim(string)//' '//oper//'imum of '//lvar
    else if ( oper == 'nti' ) then
        lvar = 'time of '//trim(string)//' minimum of '//lvar
    else if ( oper == 'xti' ) then
        lvar = 'time of '//trim(string)//' maximum  of '//lvar
    else if ( oper == 'fti' ) then
        lvar = 'time of first occurrence of '//lvar
    else if ( oper == 'lti' ) then
        lvar = 'time of last occurrence of '//lvar
    end if

end subroutine
