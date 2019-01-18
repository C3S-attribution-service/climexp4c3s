program quantiles_series

!   compute the data files necessary for the box-and-whisker plot on
!   the right side of the time series plot in the Atlas

    implicit none
    integer,parameter :: yrbeg=1850,yrend=2100,nmodmax=50,nensmax=100,nscenmax=5
    integer :: imod,nmod,iens,nens(nmodmax),iiens(nensmax,nmodmax),i,j &
        ,n,m,yr,iret,ifirst,iscen,nperyear,iperiod,mon,ave &
        ,offset,p,yr1r,yr2r,yr1s,yr2s,nperiod,pid,jens
    real :: series(12,yrbeg:yrend,nensmax,nmodmax)
    real :: seriesclim(nensmax,nmodmax),seriesscen(nensmax,nmodmax,9)
    real :: seriesmean(9,nscenmax),seriesquant(-6:6,9,nscenmax), &
        seriesdecvar(9,nscenmax),mean(12)
    real :: s,s1,s2
    logical :: lwrite,lexist,standardunits,lfirst=.true.
    character :: scenarios(nscenmax)*7
    character :: models(nmodmax)*100,shortmodel*20
    character :: string*128,format*20,vars*80,units*30,file*255 &
        ,line*255,var*40,region_subregion*100,season*10 &
        ,oneall*3,inpattern*1023,quantfile*1023,spid*10 &
        ,rcplist*30,command*1023,filelist*10000
    character :: lvars*120,svars*120,metadata(2,100)*2000,history*50000
    integer :: getpid
    data scenarios /'rcp26','rcp45','rcp60','rcp85','sresa1b'/
    lwrite = .false. 
    call getenv('LWRITE',string)
    call tolower(string)
    if ( string == 'true' ) lwrite = .true. 

!   find list of model files

    call get_command_argument(1,var)
    call get_command_argument(2,region_subregion)
    call get_command_argument(3,oneall)
    if ( oneall /= 'one' .and. oneall /= 'all' ) then
        write(0,*) 'usage: quantiles_series var region_subregion '// &
            'one|all mon N ave M yr1r yr2r yr1s yr2s '// &
            'rcplist inpattern quantfile'
        write(0,*) 'yr1r=yr2r=0 signifies no anomalies'
        call exit(-1)
    end if
    call get_command_argument(4,string)
    if ( string(1:3) /= 'mon' ) then
        write(0,*) 'error: expecting mon'
        call exit(-1)
    end if
    call get_command_argument(5,string)
    read(string,*) mon
    if ( mon < 1 .or. mon > 12 ) then
        write(0,*) 'error: expecting 1 <= mon <= 12, not ',mon
        call exit(-1)
    end if
    call get_command_argument(6,string)
    if ( string(1:3) /= 'ave' ) then
        write(0,*) 'error: expecting ave'
        call exit(-1)
    end if
    call get_command_argument(7,string)
    read(string,*) ave
    if ( ave < 1 .or. ave > 12 ) then
        write(0,*) 'error: expecting 1 <= ave <= 12, not ',ave
        call exit(-1)
    end if
    call get_command_argument(8,string)
    read(string,*) yr1r
    call get_command_argument(9,string)
    read(string,*) yr2r
    call get_command_argument(10,string)
    read(string,*) yr1s
    call get_command_argument(11,string)
    read(string,*) yr2s
    if ( (yr1r < 1800 .and. yr1r /= 0) .or. yr1s < 1800 .or. &
         (yr2r > 2300 .and. yr2r /= 0) .or. yr2s > 2300 ) then
        write(0,*) 'quantiles_series: error: yrs out of range: ',yr1r,yr2r,yr1s,yr2s
        call exit(-1)
    end if
    if ( yr1r > yr2r ) then
        write(0,*) 'quantiles_series: error: yr1r > yr2r: ',yr1r,yr2r
        call exit(-1)
    end if
    if ( yr1s > yr2s ) then
        write(0,*) 'quantiles_series: error: yr1s > yr2s: ',yr1s,yr2s
        call exit(-1)
    end if
    call get_command_argument(12,rcplist)
    call get_command_argument(13,inpattern)
    call get_command_argument(14,quantfile)

    pid = getpid()
    write(spid,'(i10.10)') pid

    write(season,'(a,i2.2,a,i2.2)') 'mon',mon,'ave',ave
    nperiod = 1
    if ( mon+ave-1 > 12 ) then
        offset = -1
    else
        offset = 0
    end if

    do iscen=1,nscenmax
        if ( lwrite ) print *,scenarios(iscen),' ',rcplist,index(rcplist,scenarios(iscen))
        if ( index(rcplist,trim(scenarios(iscen))) == 0 ) cycle
        print *,trim(scenarios(iscen))
        if ( iscen <= 4 ) then
            if ( oneall == 'all' ) then
                command = &
                    'ls '//trim(inpattern)// &
                    ' | fgrep '//trim(scenarios(iscen))// &
                    ' > /tmp/list_of_all_files_'// &
                    trim(scenarios(iscen))//'_'//trim(spid)//'.txt'
                call mysystem(trim(command),iret)
            else
                command = &
                    'ls '//trim(inpattern)// &
                    ' | fgrep '//trim(scenarios(iscen))// &
                    ' | fgrep -v EC-EARTH | '// &
                    'fgrep -v HadGEM2-ES > /tmp/list_of_one_files_'// &
                    trim(scenarios(iscen))//'_'//trim(spid)//'.txt'
            !!!print *,trim(command)
                call mysystem(trim(command),iret)
            !!! SPECIAL CASES, UPDATE WHEN MORE ENSEMBLE MEMBERS BECOME AVAILABLE
                if ( trim(scenarios(iscen)) /= 'rcp60' ) then
                    command = &
                        'ls '//trim(inpattern)// &
                        ' | fgrep '//trim(scenarios(iscen))// &
                        ' | fgrep EC-EARTH '// &
                        '| fgrep r8i1p1 >> /tmp/list_of_one_files_'// &
                        trim(scenarios(iscen))//'_'//trim(spid)//'.txt'
                !!!print *,trim(command)
                end if
                call mysystem(trim(command),iret)
                command = &
                    'ls '//trim(inpattern)// &
                    ' | fgrep '//trim(scenarios(iscen))// &
                    ' | fgrep HadGEM2-ES '// &
                    '| fgrep r2i1p1 >> /tmp/list_of_one_files_'// &
                    trim(scenarios(iscen))//'_'//trim(spid)//'.txt'
            !!!print *,trim(command)
                call mysystem(trim(command),iret)
            end if
        else
            command = &
                'ls '//trim(inpattern)// &
                ' | fgrep '//trim(scenarios(iscen))// &
                ' > /tmp/list_of_'//oneall//'_files_'// &
                trim(scenarios(iscen))//'_'//trim(spid)//'.txt'
            call mysystem(trim(command),iret)
        end if
        if ( iret /= 0 ) then
            write(0,*) 'error in system call ',iret
        !!!call exit(-1)
        end if
        nmod = 0
        nens = 0
        iens = 0
        open(1,file='/tmp/list_of_'//oneall//'_files_' &
            //trim(scenarios(iscen))//'_'//trim(spid)//'.txt',status='old')
      1 continue
        read(1,'(a)',end=2,err=2) line
        if ( index(line,'modmean') /= 0 .or. &
             index(line,'onemean') /= 0 .or. &
             index(line,'_ave') /= 0 ) go to 1
        ifirst = index(inpattern,'*')
        if ( index(inpattern,'CORDEX') == 0 ) then
            if ( iscen <= 4 ) then
                i = ifirst + index(line(ifirst:),'_') - 2
            else
                i = ifirst + index(line(ifirst:),'_144') - 2
                j = i
                if ( line(i-2:i-2) == '_' .and. &
                ichar(line(i:i)) >= ichar('0') .and. &
                ichar(line(i:i)) <= ichar('9') .and. &
                ichar(line(i-1:i-1)) >= ichar('0') .and. &
                ichar(line(i-1:i-1)) <= ichar('9') ) then
                    i = i-3
                end if
            end if
            string = line(ifirst:i)
            iens = iens + 1
            !!! YET ANOTHER SPECIAL CASE; I ASSUME ONLY GISS HAS DIFFERENT p's
            if ( string(1:4) == 'GISS' .and. index(string,'CC') == 0 ) then
                i = index(line,'i1p') + 3
                string = trim(string)//'_p'//line(i:i)
            end if
            do imod=1,nmod
                if ( trim(models(imod)) == trim(string) ) then
                    nens(imod) = nens(imod) + 1
                    call getensnumber(line,ifirst, &
                    iiens(nens(imod),imod),lwrite)
                    goto 1
                end if
            end do
        else
            ! CORDEX
            i = ifirst + index(line(ifirst:),'_v') - 2
            string = line(ifirst:i)
        end if
        nmod = nmod + 1
        models(nmod) = string
        nens(nmod) = 1
        if ( index(inpattern,'CORDEX') == 0 ) then
            call getensnumber(line,ifirst,iiens(nens(nmod),nmod),lwrite)
            if ( nmod > 1 .and. lwrite ) print *,'         with ',nens(nmod-1),' realisations ',iens
        end if
        if ( lwrite ) print *,'Found model ',nmod,' ',trim(models(nmod))
        goto 1
      2 continue
        close(1,status='delete')
        if ( nmod > 1 ) then
            if ( lwrite) print *,'         with ',nens(nmod-1),' realisations ',iens
        end if
        if ( oneall == 'all' .and. index(inpattern,'CORDEX') == 0 ) then
            print *,'Found ',nmod,' models with in all ',iens,' ensemble members<br>'
        else
            print *,'Found ',nmod,' models, using one ensemble member per model<br>'
        end if

!       read data
    
        series = 3e33
        do imod=1,nmod
            do iens=1,nens(imod)
                if ( oneall == 'one') then
                    jens = 1
                    if ( models(imod) == 'EC-EARTH' ) then
                        if ( iiens(iens,imod) /= 8 ) cycle
                    else if ( models(imod) == 'HadGEM2-ES' ) then
                        if ( iiens(iens,imod) /= 2 ) cycle
                    else
                        if ( iiens(iens,imod) /= 1 ) cycle
                    end if
                else
                    jens = iens
                end if
                if ( iiens(iens,imod) < 10 ) then
                    format = '(5a,i1,a,i1,5a)'
                else
                    format= '(5a,i2,a,i1,5a)'
                end if
                if ( index(inpattern,'CORDEX') /= 0 ) then
                    write(file,'(5a)') inpattern(:ifirst-1),trim(models(imod)), &
                        '_v1_mon_195101-209912_latlon_',trim(region_subregion),'.dat'
                    i = index(file,'rcp')
                    file(i:i+4) = trim(scenarios(iscen)) !! quick dirty hack
                    inquire(file=file,exist=lexist)
                    if ( .not. lexist ) then
                        ! some MOHC-forced runs only go to november
                        i = index(file,'209912')
                        file(i:i+5) = '209911'
                    end if
                else if ( iscen <= 4 ) then
                    shortmodel = models(imod)
                    i = index(shortmodel,'_p')
                    if ( i /= 0 ) then
                        read(shortmodel(i+2:i+2),'(i1)') p
                        shortmodel = shortmodel(1:i-1)
                    else
                        p = 1
                    end if
                                            
                    write(file,format) inpattern(:ifirst-1),trim(shortmodel),'_' &
                        ,trim(scenarios(iscen)),'_r',iiens(iens &
                        ,imod),'i1p',p,'_',trim(region_subregion),'.dat'
                    i = index(file,'rcp')
                    file(i:i+4) = trim(scenarios(iscen)) !! quick dirty hack
                else
                    write(file,'(3a,i2.2,3a)') inpattern(:ifirst-1),trim(models(imod)) &
                        ,'_',iiens(iens,imod)-1,'_144_',trim(region_subregion),'.dat'
                    inquire(file=trim(file),exist=lexist)
                    if ( .not. lexist ) then
                        write(file,'(5a)') inpattern(:ifirst-1),trim(models(imod)),'_144_' &
                            ,trim(region_subregion),'.dat'
                    end if
                end if
                if ( lwrite ) print *,'reading ',trim(file),'<br>'
                standardunits = .true.
                if ( lfirst ) then
                    lfirst = .false.
                    call readseriesmeta(file,series(1,yrbeg,jens,imod),12, &
                        yrbeg,yrend,nperyear,vars,units,lvars,svars,history,metadata, &
                        standardunits,.false.)
                    filelist = file
                else
                    call readseries(file,series(1,yrbeg,jens,imod),12 &
                        ,yrbeg,yrend,nperyear,vars,units,standardunits,.false.)
                    filelist = trim(filelist)//' '//trim(file)
                end if
                if ( index(var,'rel ') > 0 .and. yr1r > 0 ) then
                    call ensanomalclim(series(1,yrbeg,jens,imod),12 &
                        ,nperyear,yrbeg,yrend,0,0,yr1r,yr2r,mean)
                    call takerelanom(series(1,yrbeg,jens,imod),mean, &
                        12,yrbeg,yrend,0,0,nperyear,mon,ave,lwrite)
                end if
            end do
        end do
        if ( oneall == 'one' ) then
            do imod=1,nmod
                nens(imod) = 1
            end do
        end if
        if ( lwrite ) then
            print *,'metadata'
            do i=1,100
                if ( metadata(1,i) == ' ' ) exit
                print *,trim(metadata(1,i)),' :: ',trim(metadata(2,i))
            end do
            print *,'history :: ',trim(history)
            print *,'files :: ',trim(filelist)
        end if

!       compute means
    
        do imod=1,nmod
            do iens=1,nens(imod)
                call sumit(series(1,yrbeg,iens,imod),12,nperyear,yrbeg,yrend,ave,'v')
            end do
        end do
        seriesclim = 0
        seriesscen = 0
        iperiod = 1
        do imod=1,nmod
            do iens=1,nens(imod)
                n = 0
                if ( yr1r > 0 ) then
                    do yr=yr1r+offset,yr2r+offset
                        if ( series(mon,yr,iens,imod) < 1e30 ) then
                            n = n + 1
                            seriesclim(iens,imod) = &
                            seriesclim(iens,imod) &
                            + series(mon,yr,iens,imod)
                        end if
                    end do
                    if ( n == 0 ) then
                        print *,'weird, n=0 for ',models(imod)
                        do yr=yr1r+offset,yr2r+offset
                            print *,'series(',mon,yr,iens,imod,') = ',series(mon,yr,iens,imod)
                        end do
                    else
                        seriesclim(iens,imod) = &
                        seriesclim(iens,imod)/n
                    end if
                else
                    seriesclim(iens,imod) = 0
                end if
                n = 0
                do yr=yr1s+offset,yr2s+offset
                    if ( series(mon,yr,iens,imod) < 1e30 ) then
                        n = n + 1
                        seriesscen(iens,imod,iperiod) = &
                        seriesscen(iens,imod,iperiod) + series(mon,yr,iens,imod)
                    end if
                end do
                if ( n == 0 ) then
                    if ( lwrite ) print *,'weird, n = 0'
                    seriesscen(iens,imod,iperiod) = 3e33
                else
                    seriesscen(iens,imod,iperiod) = &
                    seriesscen(iens,imod,iperiod)/n - seriesclim(iens,imod)
                end if
            end do
        end do
    
!       compute s.d. of 20-yr means from inter-ensemble spread
    
        do iperiod=nperiod,1,-1
            seriesdecvar(iperiod,iscen) = 0
            m = 0
            do imod=1,nmod
                if ( nens(imod) > 2 ) then
                    m = m + 1
                    s = 0
                    s1 = 0
                    n = 0
                    do iens=1,nens(imod)
                        n = n + 1
                        s1 = s1 + seriesscen(iens,imod,iperiod)
                    end do
                    s1 = s1/n
                    s2 = 0
                    n = 0
                    do iens=1,nens(imod)
                        n = n + 1
                        s2 = s2 + &
                        (seriesscen(iens,imod,iperiod)-s1)**2
                    end do
                    seriesdecvar(iperiod,iscen) = &
                    seriesdecvar(iperiod,iscen) + s2/(n-1)
                end if
            end do
            if ( m > 0 ) then
                seriesdecvar(iperiod,iscen) = &
                seriesdecvar(iperiod,iscen)/m
            else
                seriesdecvar(iperiod,iscen) = 3e33
            end if
            if ( lwrite ) print *,'seriesdecvar(',iperiod,iscen,') = ',seriesdecvar(iperiod,iscen)
        end do
    
!       compute mean & quantiles
    
        do iperiod=nperiod,1,-1
            call gethistmean(seriesmean(iperiod,iscen), &
                seriesscen(1,1,iperiod),nens,nmod,nensmax,lwrite)
            call gethistquant(seriesquant(-6,iperiod,iscen), &
                seriesscen(1,1,iperiod),nens,nmod,nensmax,lwrite)
        end do
    end do                  ! scenario

    open(1,file=trim(quantfile))
    call printmetadata(1,' ',' ','quantiles of '//trim(vars)//' ('//trim(lvars)//')',history,metadata)
    write(1,'(2a)') '# files :: ',trim(filelist)
    write(1,'(2a)') '# year mean min 2.5% 5% 10% 17% 25% ', &
        '50% 75% 83% 90% 95% 97.5% max sd'
    do iscen=1,nscenmax
        write(1,'(2a)') '# ',scenarios(iscen)
        do iperiod=nperiod,1,-1
            write(1,'(i4,50f12.4)') yr2s,seriesmean(iperiod,iscen), &
            (seriesquant(i,iperiod,iscen),i=-6,6),sqrt(seriesdecvar(1,iscen))
        end do
        write(1,'(a)')
        write(1,'(a)')
    end do
end program quantiles_series

subroutine getensnumber(string,ifirst,n,lwrite)
    implicit none
    integer,intent(in) :: ifirst
    integer,intent(out) :: n
    character,intent(in) :: string*(*)
    logical,intent(in) :: lwrite
    integer :: i,j
    i = index(string,'sresa1b')
    if ( i /= 0 ) then      ! SRES A1b runs
        i = ifirst + index(string(ifirst:),'_144') - 2
        j = i
        if ( string(i-2:i-2) == '_' .and. &
             ichar(string(i:i)) >= ichar('0') .and. &
             ichar(string(i:i)) <= ichar('9') .and. &
             ichar(string(i-1:i-1)) >= ichar('0') .and. &
             ichar(string(i-1:i-1)) <= ichar('9') ) then
            i = i-3
        end if
        if ( j /= i ) then
            read(string(i+2:j),*) n
        else
            n = 0
        end if
        n = n + 1
    else                    ! RCP runs
        i = index(string,'_r', .true. )
        if ( i == 0 ) then
            write(0,*) 'error: cannot find "_r" in ',trim(string)
            n = 0
            return
        end if
        i = i+2
        j = i + index(string(i:),'i') - 2
!!!     print *,'searching for ensemble number in ',trim(string),i,j,string(i:j)
        if ( j < i ) then
            write(0,*) 'getensnumber: error: cannot find "i" in ',trim(string(i:))
            call exit(-1)
        end if
        read(string(i:j),*) n
    end if
    if ( lwrite ) print *,'found ensemble number ',n,' in ',trim(string)
end subroutine getensnumber

subroutine gethistmean(mean,data,nens,nmod,nensmax,lwrite)

!   Compute the mean of a set of values.
!   All models are weighted equally

    implicit none
    integer,intent(in) :: nmod,nensmax
    integer,intent(in) :: nens(nmod)
    real,intent(in) :: data(nensmax,nmod)
    real,intent(out) :: mean
    logical,intent(in) :: lwrite
    integer :: iens,imod,nmodmax
    real :: s

    if ( lwrite ) then
        print *,'gethistmean: input'
        do imod=1,nmod
            print '(i3,100f7.3)',imod,(data(iens,imod),iens=1,nens(imod))
        end do
    end if
    s = 0
    do imod=1,nmod
        do iens=1,nens(imod)
            if ( data(iens,imod) > 100 ) then
                write(0,*) 'gethistmean: warning: data(',iens,imod,') = ',data(iens,imod)
            end if
            s = s + data(iens,imod)/nens(imod)
        end do
    end do
    mean = s/nmod
    if ( lwrite ) print '(a,f7.3)','mean = ',mean
end subroutine gethistmean

subroutine gethistquant(quant,data,nens,nmod,nensmax,lwrite)

!   Compute some qunatiles of a set of values:
!   min, 5%, 25%, 50%, 75%, 95%, max
!   All models are weighted equally

    implicit none
    integer,intent(in) :: nmod,nensmax
    integer,intent(in) :: nens(nmod)
    real,intent(in) :: data(nensmax,nmod)
    real,intent(out) :: quant(-6:6)
    logical,intent(in) :: lwrite
    integer :: iens,imod,i,n,nnens(1024),nmodmax
    real :: point(1024)

!   just convert the data to the same format as quantiles_field uses

    n = 0
    do imod=1,nmod
        do iens=1,nens(imod)
            n = n + 1
            if ( n > 1024 ) then
                write(0,*) 'gethistquant: error: array too small ',1024
                call exit(-1)
            end if
            point(n) = data(iens,imod)
            nnens(n) = nens(imod)
        end do
    end do

    call getweightedquant(point,nnens,n,nmod,quant,lwrite)

end subroutine gethistquant
