        subroutine getj1j2(j1,j2,month,nperyear,lprintin)
*
*       compute first and last month (day,week, ...) of the main loop
*       based on first and last month of the season
*
        implicit none
        integer j1,j2,month,nperyear
        logical lprintin
#include "getopts.inc"
        integer m,l
        logical lprint
        character*3 months(12),seasons(4),halfyears(2)
        save months
        data months 
     +        /'Jan','Feb','Mar','Apr','May','Jun'
     +        ,'Jul','Aug','Sep','Oct','Nov','Dec'/
        data seasons /'DJF','MAM','JJA','SON'/
        data halfyears /'A-S','O-M'/

        if ( lwrite ) then
            print *,'getj1j2: month,nperyear = ',month,nperyear
            print *,'         lsel,lsum      = ',lsel,lsum
        end if
        lprint = lprintin
        if ( lweb ) lprint = .false.
        if ( nperyear.eq.1 ) then
            j1 = 1
            j2 = 1
        elseif ( month.eq.0 ) then
            j1 = 1
            j2 = nperyear
            if ( lprint ) then
                print '(a)','All year:'
                if ( dump ) write(10,'(a)') '# All year'
            endif
            corrmonths = 'all year'
        else
            if ( lsel.eq.1 ) then
                j1 = month
                j2 = month
                m = 1+mod(month+lsum-2,12)
                if ( m.eq.month ) then
                    if ( nperyear.ge.12 ) then
                        corrmonths = months(month)
                    elseif ( nperyear.eq.4 ) then
                        corrmonths = seasons(month)
                    elseif ( nperyear.eq.2 ) then
                        corrmonths = halfyears(month)
                    else
                        write(0,*)
     +                       'getj1j2: error: cannot handle nperyear = '
     +                       ,nperyear,' yet'
                        call abort
                    endif
                else
                    if ( nperyear.ge.12 ) then
                        write(corrmonths,'(3a)') months(month),'-'
     +                       ,months(1+mod(month+lsum-2,12))
                    else
                        write(corrmonths,'(3a)') seasons(month),'-'
     +                       ,seasons(1+mod(month+lsum-2,12))
                    endif
                endif
                if ( nperyear.eq.12 ) then
                    if ( lprint ) then
                        print '(2a)','Month: ',corrmonths
                        if ( dump ) write(10,'(a,12i3)') '# Month: '
     +                       ,(1+mod(month+l-1,12),l=0,lsum-1)
                    endif
                else
                    if ( lprint ) then
                        print '(4a)','Starting month: ',months(month)
                    endif
                    call month2period(j1,nperyear,1)
                    call month2period(j2,nperyear,0)
                    if ( lprint ) then
                        print '(a,i3,a,i3,a)','(Periods: ',j1,'-',
     +                       j2,')'
                        if ( dump ) write(10,'(a,12i3)')
     +                       '# Starting month: ',month
                    endif
                endif
            else                ! lsel > 1
                j1 = month
                j2 = month+lsel-1
                if ( nperyear.eq.12 .or. lsum.eq.1 ) then
                    write(corrmonths,'(3a)')
     +                   months(j1),'-',
     +                   months(1+mod(j2+lsum-2,12))
                    if ( lprint ) then
                        print '(4a)','Months: ',corrmonths
                        if ( dump ) write(10,'(4a)') '# Months: '
     +                       ,corrmonths
                    end if
                else
                    write(corrmonths,'(3a)')
     +                   months(j1),'-',
     +                   months(1+mod(j2-1,12))
                    if ( lprint ) then
                        print '(4a)','Starting months: ',corrmonths
                        if ( dump ) write(10,'(4a)')
     +                       '# Starting months: ',corrmonths
                    end if
                end if
                if ( nperyear.ne.12 ) then
                    call month2period(j1,nperyear,1)
                    call month2period(j2,nperyear,0)
                    if ( j2.lt.j1 ) j2 = j2 + nperyear
                endif
            endif
        endif
        if ( lwrite ) then
            print *,'getj1j2: j1,j2 = ',j1,j2
        end if
        end
