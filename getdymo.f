        subroutine getdymo(dy,mo,firstmo,nperyear)
!
!       given the firstmo-period of the year out of nperyear periods per
!       year, computes the day and month
!
        implicit none
        integer dy,mo,firstmo,nperyear
        integer m,i,dpm(12),dpm365(12)
        logical lwrite
        data dpm    /31,29,31,30,31,30,31,31,30,31,30,31/
        data dpm365 /31,28,31,30,31,30,31,31,30,31,30,31/

        lwrite = .false.
        if ( nperyear.le.12 ) then
            dy = 1
            mo = firstmo
            return
        endif
        m = 1+mod(firstmo-1,nperyear)
        if ( m.le.0 ) m = m + nperyear
        if ( nperyear.eq.36 ) then
            dy = 1
            do i=1,(m-1)/3
                dy = dy + dpm(i)
            enddo
            dy = dy + 5 + 10*mod(m-1,3)
        elseif ( nperyear.le.366 ) then
            dy = nint(0.5 + (m-0.5)*nint(366./nperyear))
        else
            dy = nint(0.5 + (m-0.5)*366./nperyear)
        endif
        mo = 1
 400    continue
        if ( .false. .and. lwrite ) print *,'dy,mo = ',dy,mo
        if ( nperyear.eq.365 .or. nperyear.eq.73 ) then
            if ( dy.gt.dpm365(mo) ) then
                dy = dy - dpm365(mo)
                mo = mo + 1
                goto 400
            endif
        elseif ( nperyear.eq.360 ) then
            if ( dy.gt.30 ) then
                dy = dy - 30
                mo = mo + 1
                goto 400
            endif
        else
            if ( dy.gt.dpm(mo) ) then
                dy = dy - dpm(mo)
                mo = mo + 1
                goto 400
            endif
        endif
        if ( lwrite ) then
            print *,'getdymo: input: firstmo,nperyear = ',firstmo
     +           ,nperyear
            print *,'         output: dy,mo           = ',dy,mo
        end if
        if ( mo.le.0 .or. mo.gt.12 ) then
            write(0,*) 'getdymo: error: impossible month ',mo
            mo = 1
        endif
        end

        subroutine invgetdymo(dy,mo,firstmo,nperyear)
!
!       given dy and mo, compute in which period out of nperyear ones it
!       falls
!
        implicit none
        integer dy,mo,firstmo,nperyear
        integer m,i,dpm(12),dpm365(12)
        logical lwrite
        data dpm    /31,29,31,30,31,30,31,31,30,31,30,31/
        data dpm365 /31,28,31,30,31,30,31,31,30,31,30,31/
!
        if ( nperyear.le.12 ) then
            firstmo = mo
            dy = 1
        else if ( nperyear.eq.360 ) then
            firstmo = 30*(mo-1) + dy
        else if ( nperyear.eq.365 ) then
            firstmo = 0
            do m=1,mo-1
                firstmo = firstmo + dpm365(mo)
            end do
            firstmo = firstmo + dy
        else if ( nperyear.eq.366 ) then
            firstmo = 0
            do m=1,mo-1
                firstmo = firstmo + dpm(m)
            end do
            firstmo = firstmo + dy
        else
            write(0,*) 'invgetdymo: error: cannot yet hndle nperyear = '
     +           ,nperyear
            call abort
        end if
        end
