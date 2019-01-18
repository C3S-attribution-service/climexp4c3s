program fillin

!   interpolate missing values linearly from previous and next value
!   used for the CO2 time series.

    implicit none
    include 'param.inc'
    integer :: mo,yr,mo1,yr1,mo2,yr2,k,m,nperyear,max
    real :: data(npermax,yrbeg:yrend),mean(npermax)
    character :: file*255,var*40,units*60,string*100
    logical :: lwrite
            
    lwrite = .false. 
    call get_command_argument(2,file)
    if ( file == ' ' ) go to 901
    call get_command_argument(1,string)
    read(string,*,err=901) max
    call readseries(file,data,npermax,yrbeg,yrend,nperyear, &
        var,units, .false. ,lwrite)
    call anomalclim(data,npermax,nperyear,yrbeg,yrend,yrbeg,yrend &
        ,mean)

    do yr=yrbeg,yrend
        do mo=1,nperyear
            mo1 = mo + 1
            call normon(mo1,yr,yr1,nperyear)
            if ( yr1 > yrend ) exit
            if ( data(mo,yr) < 1e33 .and. &
            data(mo1,yr1) > 1e33 ) then
                do m=1,max
                    mo2 = mo1 + m
                    call normon(mo2,yr1,yr2,nperyear)
                    if ( data(mo2,yr2) < 1e33 ) then
                        do k=1,m
                            mo1 = mo + k
                            call normon(mo1,yr,yr1,nperyear)
                            data(mo1,yr1) = ((m+1-k)*data(mo,yr) + k*data(mo2,yr2))/(m+1)
                            if ( lwrite ) write(0,*) 'filled in ',mo1,yr1,data(mo1,yr1)
                        end do
                        exit
                    end if
                end do
            end if
        end do
    end do
    do yr=yrbeg,yrend
        do mo=1,nperyear
            if ( data(mo,yr) < 1e33 .and. mean(mo) < 1e33 ) then
                data(mo,yr) = data(mo,yr) + mean(mo)
            else
                data(mo,yr) = 3e33
            end if
        end do
    end do
    call copyheader(file,6)
    print '(a)','# interpolated undefined anomalies linearly'
    call printdatfile(6,data,npermax,nperyear,yrbeg,yrend)
    goto 999
    901 continue
    write(0,*) 'usage: fillin max file'
    write(0,*) 'interpolates single missing values linearly in time'
    call exit(-1)
    999 continue
    END PROGRAM
