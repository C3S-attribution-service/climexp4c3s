        subroutine rkeepalive(i,n)
*
*       print out something every half minute to keep the httpd happy
*       secold is set at the last call to keepalive to figure out
*         when 30 seconds have passed
*       sec0 is set at the first call as a base for the wall time
*
        implicit none
        integer i,n
        integer iarray(3),sec,secold,sec0
        save secold,sec0
        real tarray(2)
        real etime
        data sec0 /99999/
*
        call itime(iarray)
        if ( sec0.eq.99999 ) then
            sec0 = iarray(3) + 60*iarray(2) + 60*60*iarray(1)
            secold = sec0
        endif
        sec = iarray(3) + 60*iarray(2) + 60*60*iarray(1)
        if ( sec.lt.secold ) then
           secold = secold - 24*60*60
           sec0 = sec0 - 24*60*60
        endif
        if ( sec-secold.gt.30 ) then
            secold = sec
            sec = sec - sec0
            iarray(1) = sec/(60*60)
            sec = sec - iarray(1)*60*60
            iarray(2) = sec/60
            sec = sec - iarray(2)*60
            iarray(3) = sec
            write(0,'(a,i8,a,i8,a,f8.1,a,i4,a,i2.2,a,i2.2,a)')
     +           'Still computing ',i,'/',n,' (CPU time ',etime(tarray)
     +           ,'s, wall time ',iarray(1),':',iarray(2),':',iarray(3),
     +           ')<p>'
        endif
        end
