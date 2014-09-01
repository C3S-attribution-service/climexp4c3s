        subroutine keepalive(i,n)
        implicit none
        integer i,n
        call keepalive1('Still computing, ',i,n)
        end
        subroutine keepalive1(string,i,n)
*
*       print out something every half minute to keep the httpd happy
*       secold is set at the last call to keepalive to figure out
*         when 30 seconds have passed
*       sec0 is set at the first call as a base for the wall time
*
        implicit none
        integer i,n
        character string*(*)
        integer iarray(8),sec,secold,sec0
        save secold,sec0
        real tarray(2)
        real etime
        data sec0 /99999/
*
        call date_and_time(values=iarray)
        if ( sec0.eq.99999 ) then
            sec0 = iarray(7) + 60*iarray(6) + 60*60*iarray(5)
            secold = sec0
        endif
        sec = iarray(7) + 60*iarray(6) + 60*60*iarray(5)
        if ( sec.lt.secold ) then
           secold = secold - 24*60*60
           sec0 = sec0 - 24*60*60
        endif
        !!!write(0,*) 'keepalive: delta(sec) = ',sec-secold
        if ( sec-secold.gt.30 ) then
            secold = sec
            sec = sec - sec0
            iarray(5) = sec/(60*60)
            sec = sec - iarray(5)*60*60
            iarray(6) = sec/60
            sec = sec - iarray(6)*60
            iarray(7) = sec
            write(0,'(2a,i8,a,i8,a,f8.1,a,i4,a,i2.2,a,i2.2,a)')
     +           trim(string),' ',i,'/',n,' (CPU time ',etime(tarray)
     +           ,'s, wall time ',iarray(5),':',iarray(6),':',iarray(7),
     +           ')<p>'
            flush(0)
        endif
        end
