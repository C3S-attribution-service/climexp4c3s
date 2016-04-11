real function wetbulbmax(tmax,temp,rh,pres,lwrite)
!
!   estimate daily max wet bulb temperature from Tmax, Tmean and RH
!
    implicit none
    real tmax,temp,rh,pres
    logical lwrite
    real es,e2
    real,external :: wetbulbvap

    es = 6.112 * exp(17.67 * temp / (temp + 243.5))
    e2 = rh*es/100
    if ( lwrite ) print *,'e2 = ',e2
    wetbulbmax = wetbulbvap(tmax,pres,e2,lwrite)

    if ( lwrite ) then
        print *,'T,Tx    = ',temp,tmax
        print *,'e,es,rh = ',e2,es,rh
        print *,'Tw      = ',wetbulbmax
    end if
end function

real function wetbulbdew(temp,pres,dewtemp,lwrite)
!
!   compute wet bulb temperature from temperature (Celsius), pressure (mb)
!   and dew point temperature
!
    implicit none
    real temp,pres,dewtemp
    logical lwrite
    real e2,es,rh
    real,external :: wetbulbvap
    
    if ( lwrite ) then
        print *,' wetbulbdew: temp,pres,dewtemp = ',temp,pres,dewtemp
    end if
    if ( temp > 1e33 .or. pres > 1e33 .or. dewtemp > 1e33 ) then
        wetbulbdew = 3e33
        return
    end if
    es = 6.112 * exp(17.67 * temp / (temp + 243.5))
    e2 = 6.112 * exp(17.67 * dewtemp / (dewtemp + 243.5))
    rh = e2/es*100
    if ( lwrite ) print *,'rh = ',rh
    wetbulbdew = wetbulbvap(temp,pres,e2,lwrite)

    if ( lwrite ) then
        print *,'T,Td    = ',temp,dewtemp
        print *,'e,es,rh = ',e2,es,rh
        print *,'Tw      = ',wetbulbdew
    end if
end function

real function wetbulbsim(temp,dewtemp,lwrite)
!
!   compute simplified wet bulb temperature from temperature (Celsius)
!   and dew point temperature
!
    implicit none
    real temp,dewtemp
    logical lwrite
    real e2,es,rh
    real,external :: wetbulbvap
    
    if ( lwrite ) then
        print *,' wetbulbdew: temp,dewtemp = ',temp,dewtemp
    end if
    if ( temp > 1e33 .or. dewtemp > 1e33 ) then
        wetbulbsim = 3e33
        return
    end if
    es = 6.112 * exp(17.67 * temp / (temp + 243.5))
    e2 = 6.112 * exp(17.67 * dewtemp / (dewtemp + 243.5))
    rh = e2/es*100
    if ( lwrite ) print *,'rh = ',rh
    wetbulbsim = 0.567*temp + 0.393*e2 + 3.94

    if ( lwrite ) then
        print *,'T,Td    = ',temp,dewtemp
        print *,'e,es,rh = ',e2,es,rh
        print *,'Tw      = ',wetbulbsim
    end if
end function

real function wetbulbvap(temp,pres,E2,lwrite)
!
!   compute wet bulb temperature from temperature (Celsius), pressure (mb)
!   and vapour pressure
!   based on javascript at http://www.srh.noaa.gov/epz/?n=wxcalc_rh
!   (not the most sophisticated solver I have ever seen.._
!
    implicit none
    real temp,pres,E2
    logical lwrite
    integer i,prevsign,cursign
    real Edifference,incr,Twguess,Ewguess,Eguess
!
!   init
!
    if ( lwrite ) then
        print *,' wetbulbvap: temp,pres,E2 = ',temp,pres,E2
    end if
    if ( temp > 1e33 .or. pres > 1e33 .or. E2 > 1e33 ) then
        wetbulbvap = 3e33
        return
    end if
    Edifference = 1
    prevsign = 1
    incr = 10
    Twguess = 0

    i = 0
    do while ( ( i < 10 .or. abs(Edifference) > 0.0001 ) .and. i < 1000 )
        i = i + 1
        Ewguess = 6.112 * exp((17.67 * Twguess) / (Twguess + 243.5))
        Eguess = Ewguess - pres * (temp - Twguess) * 0.00066 * (1 + (0.00115 * Twguess))
        if ( lwrite ) print *,'Ewguess,Eguess,Twguess,incr = ',Ewguess,Eguess,Twguess,incr
        Edifference = E2 - Eguess
        if ( Edifference == 0 ) then
            exit
        else if ( Edifference < 0 ) then
            cursign = -1
            if ( cursign /= prevsign ) then
                prevsign = cursign
                incr = incr/10
            end if
        else ! Edifference > 0
            cursign = +1
            if ( cursign /= prevsign ) then
                prevsign = cursign
                incr = incr/10
            end if
        end if
        Twguess = Twguess + incr * prevsign
    end do
    if ( i >=1000 ) then
        write(0,*) 'wetbulb: did not converge ',temp,pres,E2
    end if
    if ( lwrite ) print *,'wetbulb = ',Twguess
    wetbulbvap = Twguess
end function