program month2lead

!   convert files with indices sorted by analysis month
!   into files sorted by lead time

    implicit none
    integer :: yrbeg,yrend,npermax
    parameter(yrbeg=1940,yrend=2020,npermax=12)
    integer :: i,j,k,mo,lead,yr,nplus,nperyear,mo1,yr1,leadmax,ana
    real :: data(12,yrbeg:yrend,12),newdata(12,yrbeg:yrend,0:11)
    character file*255,months(12)*3,var*40,units*20,lines(10)*255
    logical :: lstandardunits,lwrite
    data months /'jan','feb','mar','apr','may','jun', &
    'jul','aug','sep','oct','nov','dec'/
    lwrite = .false. 
    lstandardunits = .false. 

    if ( command_argument_count() < 1 ) then
        print *,'usage: month2lead files_+++_rest.dat'
        print *,'       replaces +++ with jan feb ... dec'
        print *,'       and produces files with lead time 0 ...'
    endif

    call get_command_argument(1,file)
    nplus = index(file,'+++')
    if ( nplus == 0 ) then
        write(0,*) 'error: expecting ''+++'' in file name'
        call exit(-1)
    endif

!       read series

    do ana=1,12
        file(nplus:nplus+2) = months(ana)
        call readseries(file,data(1,yrbeg,ana),npermax,yrbeg,yrend &
        ,nperyear,var,units,lstandardunits,lwrite)
    enddo
    lines = ' '
    open(1,file=file)
    do i=1,10
        read(1,'(a)',end=80) lines(i)
    enddo
    80 continue

!       reorder data

    leadmax=0
    do lead=0,11
        do mo=1,12
            do yr=yrbeg,yrend
                ana = mo - lead
                if ( ana <= 0 ) ana = ana + 12
                newdata(mo,yr,lead) = data(mo,yr,ana)
                if ( newdata(mo,yr,lead) < 1e33 ) then
                    leadmax = max(leadmax,lead)
                endif
            enddo
        enddo
    enddo

!       print data

    do lead=0,leadmax
        write(file(nplus:nplus+2),'(a,i2.2)') '+',lead
        open(1,file=file)
        do i=1,10
            if ( lines(i)(1:1) == '#' ) write(1,'(a)')trim(lines(i))
        enddo
        write(1,'(a,i2.2)') '# lead time ',lead
        call printdatfile(1,newdata(1,yrbeg,lead),npermax,nperyear, &
        yrbeg,yrend)
        close(1)
    enddo

!       also make seasonal averages

    do lead=0,leadmax-2
        do yr=yrbeg,yrend
            do mo=1,12
                mo1 = mo + 1
                call normon(mo1,yr,yr1,12)
                if ( yr1 > yrend ) then
                    newdata(mo,yr,lead) = 3e33
                elseif ( newdata(mo,yr,lead) < 1e33 .and. &
                    newdata(mo1,yr1,lead+1) < 1e33 ) then
                    newdata(mo,yr,lead) = newdata(mo,yr,lead) + &
                    newdata(mo1,yr1,lead+1)
                else
                    newdata(mo,yr,lead) = 3e33
                endif
                mo1 = mo + 2
                call normon(mo1,yr,yr1,12)
                if ( yr1 > yrend ) then
                    newdata(mo,yr,lead) = 3e33
                elseif ( newdata(mo,yr,lead) < 1e33 .and. &
                    newdata(mo1,yr1,lead+1) < 1e33 ) then
                    newdata(mo,yr,lead) = (newdata(mo,yr,lead) + &
                    newdata(mo1,yr1,lead+2))/3
                else
                    newdata(mo,yr,lead) = 3e33
                endif
            enddo
        enddo
    enddo

    i = index(file,'.dat')
    file(i:) = '_3m.dat'
    do lead=0,leadmax-2
        write(file(nplus:nplus+2),'(a,i2.2)') '+',lead
        open(1,file=file)
        do i=1,10
            if ( lines(i)(1:1) == '#' ) write(1,'(a)')trim(lines(i))
        enddo
        write(1,'(a,i2.2)') '# lead time ',lead
        write(1,'(a)') '# 3-month avrages '
        call printdatfile(1,newdata(1,yrbeg,lead),12,12,yrbeg,yrend)
        close(1)
    enddo
            
end program month2lead
