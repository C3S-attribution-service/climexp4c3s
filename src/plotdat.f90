program plotdat

!   convert 14-column data into a format suitable for gnuplot,
!   taking anomalies if requested

    implicit none
    include 'param.inc'
    include 'getopts.inc'
    integer :: i,j,k,l,n,year,loop,type,nperyear,iens,iarg,mens,mens1,ayr1,ayr2
    real :: mean(npermax),s,nextx,lasty
    real,allocatable :: data(:,:,:)
    logical :: lanomal,lskip,lexist,lopened
    character :: line*1023,ensfile*255,var*40,units*40,name*1023,title*100
    character :: lvar*120,svar*120,history*50000,metadata(2,100)*1000
    lwrite = .false. 
    lskip = .false. 
    call getenv('PLOTDAT_LWRITE',line)
    if ( line(1:1) == 't' .or. line(1:1) == 'T' ) then
        lwrite = .true. 
    endif

!   init

    if ( command_argument_count() == 0 ) then
        print *,'usage: plotdat [anom [yr1 yr2]] [normsd [[mon M ] ave N]] datafile'
        stop
    endif

    iarg = 1
    ayr1 = yrbeg
    ayr2 = yrend
    lanomal = .false. 
    iarg = 1
    if ( command_argument_count() > 1 ) then
        call get_command_argument(1,line)
        if ( line(1:4) == 'anom' ) then
            lanomal = .true. 
            if ( command_argument_count() < 3 ) then
                iarg = 2
            else
                call get_command_argument(2,line)
                read(line,*,err=901) ayr1
                ayr1 = max(yrbeg,ayr1)
                call get_command_argument(3,line)
                read(line,*,err=901) ayr2
                ayr2 = min(yrend,ayr2)
                iarg = 4
            endif
        else
            lanomal = .false. 
            iarg = 1
        end if
    endif
      
    call get_command_argument(command_argument_count(),line)
    call keepalive1('Reading data',0,0)
    lstandardunits = .false. ! because we cannot call getopts before we know mens1,mens...
    allocate(data(npermax,yrbeg:yrend,0:nensmax))
    call readensseriesmeta(line,data,npermax,yrbeg,yrend,nensmax &
        ,nperyear,mens1,mens,var,units,lvar,svar,history,metadata,lstandardunits,lwrite)
    call keepalive1('Reading data',0,0)
    call getopts(iarg,command_argument_count()-1,nperyear,yrbeg,yrend,.true.,mens1,mens)
    if ( index(line,'++') > 0 .or. index(line,'%%') > 0 ) then
        inquire(file=trim(line),exist=lexist)
        if ( .not. lexist .or. index(line,'.nc') == 0 ) then
            name = line
            call filloutens(line,0)
            inquire(file=trim(line),exist=lexist)
            if ( .not. lexist ) then
                line = name
                call filloutens(line,1)
            end if
        end if
    end if
    title = ' '
    if ( svar /= ' ' ) then
        do n=1,100
            if ( metadata(1,n) == ' ' ) exit
        end do
        if ( n <= 100 ) then
            metadata(1,n) = 'standard_name'
            metadata(2,n) = svar
        end if
    end if
    call copyheadermeta(line,6,title,history,metadata)
    if ( lnormsd ) then
        if ( m1 == 0 ) then
            print '(2a)','# taking relative anomalies wrt standard seasons'
        else
            print '(a,i3,a,i3)','# taking relative anomalies wrt ', &
                lsum,'-month season starting in month ',m1
        end if
        var = trim(var)//'rel'
        units = '1'
        lvar = 'relative '//trim(lvar)
    end if
    if ( var /= ' ' ) then
        print '(6a)','# ',trim(var),' [',trim(units),'] ',trim(lvar)
    end if
    if ( lstandardunits ) then
        do iens=mens1,mens
            call makestandardseries(data(1,yrbeg,iens), &
                npermax,yrbeg,yrend,nperyear,var,units,lwrite)
        end do
    endif

    if ( lanomal ) then
        if ( lensanom ) then
            call ensanomalclim(data,npermax,nperyear,yrbeg,yrend, &
                nens1,nens2,ayr1,ayr2,mean)
            if ( lnormsd ) then
                call takerelanom(data,mean,npermax,yrbeg,yrend, &
                    nens1,nens2,nperyear,m1,lsum,lwrite)
            end if
        else
            do iens=nens1,nens2
                call ensanomalclim(data,npermax,nperyear,yrbeg,yrend &
                    ,iens,iens,ayr1,ayr2,mean)
                if ( lnormsd ) then
                    call takerelanom(data,mean,npermax,yrbeg,yrend, &
                        iens,iens,nperyear,m1,lsum,lwrite)
                end if
            end do
        end if
    end if                  ! lanomal

!   print

    do iens=nens1,nens2
        call keepalive1('Writing data',iens-nens1+1,nens2-nens1+1)
        if ( nens2 > nens1 ) then
            print '(a,i4)','# ensemble member ',iens
            print '(a)'
        endif
        do i=yrbeg,yrend
            do j=1,nperyear
                k = j
                n = nperyear
                if ( nperyear == 366 .and. (mod(i,4) /= 0 .or. ( &
                     mod(i,100) == 0 .and. mod(i,400) /= 0)) ) then
                    n = 365
                    if ( j > 60 ) k = j-1
                    if ( j == 60 .and. data(j,i,iens) < 1e33 ) then
                        print '(a,i4,f14.6)','# skipping valid data on 29feb',i,data(j,i,iens)
                    endif
                endif
                if ( abs(data(j,i,iens)) < 1e33 ) then
                    print '(f10.4,g14.6)',i+(k-1.)/n,data(j,i,iens)
                    nextx = i+real(k)/n
                    lasty = data(j,i,iens)
                    lskip = .true. 
                elseif ( lskip .and. .not. ( nperyear == 366 .and. j == 60) ) then
!                   to counter peculiarity in gnuplot steps plotting
                    print '(f10.4,g14.6,a)',nextx,lasty,'# repeat last y to get nice gnuplot plot'
                    print '(a)'
                    lskip = .false. 
                endif
            enddo
        enddo
    enddo

    goto 999
901 write(0,*) 'plotdat: expecting year, found ',line
    call exit(-1)
999 continue
end program

