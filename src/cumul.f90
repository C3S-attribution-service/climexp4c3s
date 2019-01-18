program cumul
!
!   compute cumulative time series
!
    implicit none
    include 'param.inc'
    integer yr,mo,nperyear
    real data(npermax,yrbeg:yrend),cum(npermax,yrbeg:yrend),sum
    logical lwrite
    character file*1024,var*100,units*50
    
    lwrite = .false.
    call get_command_argument(1,file)
    call readseries(file,data,npermax,yrbeg,yrend,nperyear, &
        var,units,.false.,lwrite)
    sum = 0
    do yr=yrbeg,yrend
        do mo=1,npermax
            if ( data(mo,yr) < 1e33 ) then
                sum = sum + data(mo,yr)
                cum(mo,yr) = sum
            else
                cum(mo,yr) = 3e33
            end if
        end do
    end do
    print '(a)','# cumulative of'
    call copyheader(file,6)
    call printdatfile(6,cum,npermax,nperyear,yrbeg,yrend)
end program cumul