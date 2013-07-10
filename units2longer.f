        subroutine units2longer(units,nperyear,nperyearnew)
*
*       adjust units to reflect the fact that we have summed from
*       nperyear to nperyearenw
*
        implicit none
        integer nperyear,nperyearnew
        character units*(*)
        integer i,j
        character*10 timeunits
*
        if ( index(units,'/day')+index(units,'/dy').ne.0 .and.
     +       ( nperyear.eq.360 .or. nperyear.eq.365 .or. 
     +         nperyear.eq.366 ) ) then
*           erase the /day and replace by new unit
            call nperyear2units(nperyearnew,timeunits)
            i = index(units,'/day')
            j = 4
            if ( i.eq.0 ) then
                i = index(units,'/dy')
                j = 3
            endif
            units = units(:i-1)//'/'//trim(timeunits)//units(i+j:)
        elseif ( index(units,'/month').ne.0 .and. nperyear.eq.12 ) then
            i = index(units,'/month')
            call nperyear2units(nperyearnew,timeunits)
            units = units(:i-1)//'/'//trim(timeunits)//units(i+6:)
        elseif ( index(units,'/year').ne.0 .and. nperyear.eq.12 ) then
            i = index(units,'/year')
            call nperyear2units(nperyearnew,timeunits)
            units = units(:i-1)//'/'//trim(timeunits)//units(i+5:)
        else
            call nperyear2units(nperyear,timeunits)
            units = trim(units)//' '//trim(timeunits)
        endif
        end

        subroutine nperyear2units(nperyear,timeunits)
        implicit none
        integer nperyear
        character timeunits*(*)
        integer i
*
        if ( nperyear.eq.1 ) then
            timeunits = 'year'
        elseif ( nperyear.eq.4 ) then
            timeunits = 'season'
        elseif ( nperyear.eq.12 ) then
            timeunits = 'month'
        elseif ( nperyear.lt.360 ) then
            write(timeunits,'(i6)') nint(365.24/nperyear)
            do i=1,6
                if ( timeunits(1:1).eq.' ' ) timeunits = timeunits(2:)
            enddo
            timeunits = '('//trim(timeunits)//'days)'
        elseif ( nperyear.eq.360 .or. nperyear.eq.365 .or.
     +           nperyear.eq.365 ) then
            timeunits = 'day'
        else
            write(timeunits,'(i6)') nint(24*365.24/nperyear)
            do i=1,6
                if ( timeunits(1:1).eq.' ' ) timeunits = timeunits(2:)
            enddo
            timeunits = '('//trim(timeunits)//'hours)'
        endif
        end
