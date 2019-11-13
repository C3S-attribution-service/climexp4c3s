subroutine keepalive(i,n)
    implicit none
    integer,intent(in) :: i,n
    call keepalive1('Still computing, ',i,n)
end subroutine keepalive

subroutine keepalive1(string,i,n)
    implicit none
    integer,intent(in) :: i,n
    character,intent(in) :: string*(*)
    call keepalive2(string,i,n,.false.)
end subroutine keepalive1

subroutine keepalive2(string,i,n,lforce)

!   print out something every half minute to keep the httpd happy
!   secold is set at the last call to keepalive to figure out when 30 seconds have passed
!   sec0 is set at the first call as a base for the wall time

    implicit none
    integer,intent(in) :: i,n
    character,intent(in) :: string*(*)
    logical,intent(in) :: lforce
    integer :: iarray(8),sec
    integer,save :: secold,sec0=99999
    character,save :: laststring*50
    real :: tarray(2)
    real :: etime
    if ( sec0 == 99999 ) then
        laststring = ' '
    end if

    call date_and_time(values=iarray)
    if ( sec0 == 99999 ) then
        sec0 = iarray(7) + 60*iarray(6) + 60*60*iarray(5)
        secold = sec0
    endif
    sec = iarray(7) + 60*iarray(6) + 60*60*iarray(5)
    if ( sec < secold ) then
        secold = secold - 24*60*60
        sec0 = sec0 - 24*60*60
    endif
!!!write(0,*) 'keepalive: delta(sec) = ',sec-secold
    if ( sec-secold > 30 .or. (lforce .and. laststring /= string) ) then
        secold = sec
        laststring = string
        sec = sec - sec0
        iarray(5) = sec/(60*60)
        sec = sec - iarray(5)*60*60
        iarray(6) = sec/60
        sec = sec - iarray(6)*60
        iarray(7) = sec
        if ( n < 0 ) then
            write(0,'(2a,i8,a,f8.1,a,i4,a,i2.2,a,i2.2,a)') &
                trim(string),' ',i,' (CPU time ',etime(tarray) &
                ,'s, wall time ',iarray(5),':',iarray(6),':',iarray(7), &
                ')<p>'
        else
            write(0,'(2a,i8,a,i8,a,f8.1,a,i4,a,i2.2,a,i2.2,a)') &
                trim(string),' ',i,'/',n,' (CPU time ',etime(tarray) &
                ,'s, wall time ',iarray(5),':',iarray(6),':',iarray(7), &
                ')<p>'
        end if
        if ( sec > 10*60*60 ) then ! no calculations should last more than 10 hours
            write(0,'(2a)') trim(string),': maximum run time is 10 hours, aborted'
            call exit(-1)
        end if
        flush(0)
    endif
end subroutine keepalive2
