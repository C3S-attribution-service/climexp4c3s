program testwetbulb
    implicit none
    real temp,dewtemp,pres,wetbulb
    real,external :: wetbulbdew
    logical lwrite
    lwrite = .true.
1   print *,'T,Tdew,P (Celsius,hPa)'
    read *,temp,dewtemp,pres
    !!!pres = 1000
    wetbulb = wetbulbdew(temp,pres,dewtemp,lwrite)
    print *,'Twetbulb = ',wetbulb
    goto 1
end program