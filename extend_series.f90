program extendseries
!
!   extend the series by one month to enable real-time use fo teh global eman temperature,
!   assuming persistence
!
    implicit none
    integer yrbeg,yrend,npermax
    parameter(yrbeg=1500,yrend=2100,npermax=12)
    integer yr,mo,yr1,mo1,nperyear
    real data(npermax,yrbeg:yrend)
    character file*255,var*80,units*40
    logical lstandardunits,lwrite
    lwrite = .false.
    call get_command_argument(1,file)
    lstandardunits = .false.
    call readseries(file,data,npermax,yrbeg,yrend,nperyear,var,units,lstandardunits,lwrite)
    do yr=yrend,yrbeg,-1
        do mo=12,1,-1
            if ( data(mo,yr).lt.3e33 ) then
                mo1 = mo + 1
                call normon(mo1,yr,yr1,nperyear)
                data(mo1,yr1) = data(mo,yr)
                goto 800
            end if
        end do
    end do
    800 continue
    call copyheader(file,6)
    call printdatfile(6,data,npermax,nperyear,yrbeg,yrend)
end program