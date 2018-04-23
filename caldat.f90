subroutine caldat(julian,month,day,year)
!
!1wrapper around conversion routine from https://stackoverflow.com/questions/39719011/convert-between-julian-and-gregorian-in-fortran
!
    implicit none
    integer :: julian,month,day,year
    integer :: hour,minute,second
    
    call gregorian(dble(julian),year,month,day,hour,minute,second)

end subroutine caldat

subroutine gregorian(jd,year,month,day,hour,minute,second)

    implicit none
    double precision, intent(in) :: jd
    integer, intent(out) :: year, month, day, hour, minute, second
    double precision :: jt
    integer :: i,j,k,l,n

    l = int(jd)+68569
    n = 4*l/146097
    l = l-(146097*n+3)/4
    i = 4000*(l+1)/1461001
    l = l-1461*i/4+31
    j = 80*l/2447
    k = l-2447*j/80
    l = j/11
    j = j+2-12*l
    i = 100*(n-49)+i+l

    year = i
    month = j
    day = k

    jt = dmod(jd,1.d0)*24.d0
    hour = int(jt)
    jt = dmod(jt,1.d0)*60.d0
    minute = int(jt)
    jt = dmod(jt,1.d0)*60.d0
    second = nint(jt)

    if (second == 60) then
        second = second-60
        minute = minute+1
    end if

end subroutine gregorian
