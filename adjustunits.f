        subroutine adjustunits(oper,nperyear,nperyearnew,units,lwrite)
*
*       adjust units
*
        implicit none
        integer nperyear,nperyearnew
        character oper*3,units*(*)
        logical lwrite
*
        if ( oper.eq.'mea' .or. oper.eq.'sd ' .or.
     +       oper.eq.'min' .or. oper.eq.'max' ) then
            units = units
        else if ( oper.eq.'var' ) then
            if ( index(trim(units),' ').eq.0 .and. 
     +           index(trim(units),'*').eq.0 .and.
     +           index(trim(units),'/').eq.0 .and.
     +           index(trim(units),'^').eq.0 ) then
                units = trim(units)//'^2'
            else
                units = '('//trim(units)//')^2'
            end if
        else if ( oper.eq.'sum' .or. oper.eq.'bel' .or. oper.eq.'abo' )
     +           then
            if ( lwrite ) print *,'changed units from ',trim(units)
            call units2longer(units,nperyear,nperyearnew)
            if ( lwrite ) print *,'to ',trim(units)
        else if ( oper.eq.'num' ) then
            call nperyear2units(nperyearnew,units)
            units = '1/'//units
            if ( lwrite ) print *,'changed units to ',trim(units)
        else if ( oper.eq.'xti' .or. oper.eq.'nti' .or.
     +           oper.eq.'fti' .or. oper.eq.'lti' .or.
     +           oper.eq.'con' ) then
            if ( nperyear.eq.360 .or. nperyear.eq.365 .or. 
     +           nperyear.eq.366 ) then
                units = 'dy'
            else if ( nperyear.eq.12 ) then
                units = 'mo'
            else if ( nperyear.eq.4 ) then
                units = 'season'
            else
                units = 'period'
            end if
        else
            write(0,*) 'adjustunits: error: unknown operation ',oper
            call abort
        endif
        end subroutine

        subroutine adjustvar(oper,var,lwrite)
!
!       adjust variable name
!
        implicit none
        character oper*3,var*(*)
        logical lwrite
!
        if ( oper.eq.'num' ) then
            var = 'n'//var
        else if ( oper.eq.'max' .or. oper.eq.'min' ) then
            var = oper//'_'//var
        else if ( oper.eq.'nti' ) then
            var = 'time_min_'//var
        else if ( oper.eq.'xti' ) then
            var = 'time_max_'//var
        else if ( oper.eq.'fti' ) then
            var = 'time_first_'//var
        else if ( oper.eq.'lti' ) then
            var = 'time_last_'//var
        end if
!
        end subroutine
