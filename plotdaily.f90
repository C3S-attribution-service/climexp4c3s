program plotdaily
!
!   make the graphics file for gnuplot to plot the last N days of a daily time series
!   plus climatology.
!
    implicit none
    include 'param.inc'
    include 'getopts.inc'
    integer nday,nperyear,mens,mens1,i,j,k,dy,mm,mo,yy,yr,yrlast,molast,dylast
    real,allocatable :: data(:,:),mean(:)
    character file*255,var*20,units*20,string*80,enddate*20
    integer iargc
    
    if ( iargc() < 3 ) then
        write(0,*) 'usage: plotdaily infile nday enddate [begin yr1 end yr2]'
        call exit(-1)
    end if
    
    call getarg(1,file)
    call getarg(2,string)
    call getarg(3,enddate)
    read(string,*) nday
    allocate(data(npermax,yrbeg:yrend))
    call readseries(file,data,npermax,yrbeg,yrend,nperyear, &
    &   var,units,lstandardunits,lwrite)
    mens1 = 0
    mens = 0
    call getopts(4,iargc(),nperyear,yrbeg,yrend,.true.,mens1,mens)

    allocate(mean(nperyear))
    call anomalclim(data,npermax,nperyear,yrbeg,yrend,yr1,yr2,mean)

    if ( enddate == 'last' ) then
        ! find last time with data
        do yr=yrend,yrbeg,-1
            do mo=nperyear,1,-1
                if ( data(mo,yr).lt.1e33 ) goto 800
            end do
        end do
        800 continue
        yrlast = yr
        molast = mo
    else
        read(enddate,'(i4,2i2.2)') yrlast,mo,dy
        call invgetdymo(dy,mo,molast,nperyear)
    end if

    print '(a,i5,2a)','# last ',nday,' of data in file ',trim(file)
    call copyheader(file,6)
    do k=nday-1,0,-1
        mm = molast-k
        call normon(mm,yrlast,yr,nperyear)
        call getdymo(dy,mo,mm,nperyear)
        if ( data(mm,yr).lt.1e33 .and. mean(mm).lt.1e33 ) then
            print '(i4,2i2.2,2f12.4)',yr,mo,dy,data(mm,yr)+mean(mm),mean(mm)
        else
            print '(a)'
        end if
    end do
end program