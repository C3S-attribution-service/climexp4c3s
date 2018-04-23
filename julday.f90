integer function julday(month,day,year)
!   wrapper for julian
    implicit none
    integer,intent(in) :: month,day,year
    integer,external :: julian
    julday = julian(year,month,day)
end function julday

integer function julian(year,month,day)
!   from https://stackoverflow.com/questions/39719011/convert-between-julian-and-gregorian-in-fortran%20!
    implicit none
    integer,intent(in) :: year,month,day
    julian = day-32075+1461*(year+4800+(month-14)/12)/4 + &
        367*(month-2-(month-14)/12*12)/12 - &
        3*((year+4900+(month-14)/12)/100)/4
end function julian
