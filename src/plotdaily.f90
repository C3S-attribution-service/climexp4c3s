program plotdaily
!
!   make the graphics file for gnuplot to plot the last N days of a daily time series
!   plus climatology.
!   also works for monthly and annual data, so the name is a bit awkward.
!
    implicit none
    include 'param.inc'
    include 'getopts.inc'
    integer :: nday,nperyear,mens,mens1,i,j,k,dy,mm,mo,yy,yr,yrlast,molast,dylast,n
    real :: cumdata,cummean
    real,allocatable :: data(:,:),mean(:)
    logical :: cdf
    character :: file*255,var*20,units*20,string*80,enddate*20
    character :: lvar*120,svar*120,history*50000,metadata(2,100)*2000
    integer,external :: leap
    
    if ( command_argument_count() < 3 ) then
        write(0,*) 'usage: plotdaily infile nday enddate [cdf] [begin yr1 end yr2] [anom]'
        call exit(-1)
    end if
    
    call get_command_argument(1,file)
    call get_command_argument(2,string)
    read(string,*) nday
    call get_command_argument(3,enddate)
    allocate(data(npermax,yrbeg:yrend))
    lstandardunits = .true.
    lwrite = .false.
    call readseriesmeta(file,data,npermax,yrbeg,yrend,nperyear,var,units, &
        lvar,svar,history,metadata,lstandardunits,lwrite)
    mens1 = 0
    mens = 0
    n = 4
    cdf = .false.
    if ( command_argument_count() >= n ) then
        call get_command_argument(n,string)
        if ( string(1:3) == 'cdf' ) then
            n = n + 1
            cdf = .true.
        end if
    end if
    call getopts(n,command_argument_count(),nperyear,yrbeg,yrend,.true.,mens1,mens)

    allocate(mean(nperyear))
    if ( .not.anom ) then
        call anomalclim(data,npermax,nperyear,yrbeg,yrend,yr1,yr2,mean)
    else
        mean = 0
    end if

    if ( enddate == 'last' ) then
        ! find last time with data
        do yr=yrend,yrbeg,-1
            do mo=nperyear,1,-1
                if ( data(mo,yr).lt.1e33 ) goto 800
            end do
        end do
        write(0,*) 'plotdaily: error: cannot find any valid data'
        call exit(-1)
    800 continue
        yrlast = yr
        molast = mo
    else
        read(enddate,'(i4,2i2.2)') yrlast,mo,dy
        call invgetdymo(dy,mo,molast,nperyear)
    end if

    print '(a,i5,2a)','# last ',nday,' of data in file ',trim(file)
    call printvar(6,var,units,lvar)
    call printmetadata(6,file,' ',' ',history,metadata)
    if ( cdf ) then
        cumdata = 0
        cummean = 0
    end if
    do k=nday-1,-1,-1
        mm = molast-k
        call normon(mm,yrlast,yr,nperyear)
        call getdymo(dy,mo,mm,nperyear)
        if ( yr >= yrbeg .and. yr <= yrend ) then
            if ( data(mm,yr).lt.1e33 .and. mean(mm).lt.1e33 ) then
                if ( cdf ) then
                    cumdata = cumdata + data(mm,yr)+mean(mm)
                    cummean = cummean + mean(mm)
                else
                    cumdata = data(mm,yr)+mean(mm)
                    cummean = mean(mm)
                end if
                print '(i4,2i2.2,2g14.6)',yr,mo,dy,cumdata,cummean
            else if ( nperyear /= 366 .or. leap(yr) == 2 .or. mm /= 60 ) then
                print '(a)'
            end if
        end if
    end do
end program