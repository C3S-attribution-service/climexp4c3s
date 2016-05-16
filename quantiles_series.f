        program quantiles_series
!
!       compute the data files necessary for the box-and-whisker plot on
!       the right side of the time series plot
!
        implicit none
        integer yrbeg,yrend,nmodmax,nensmax,nscenmax
        parameter (yrbeg=1850,yrend=2100,nmodmax=50,nensmax=100,
     +       nscenmax=5)
        integer imod,nmod,iens,nens(nmodmax),iiens(nensmax,nmodmax),i,j
     +       ,n,m,yr,iret,ifirst,iscen,nperyear,iperiod,mon,ave
     +       ,offset,p,yr1r,yr2r,yr1s,yr2s,nperiod,pid,jens
        real series(12,yrbeg:yrend,nensmax,nmodmax)
        real seriesclim(nensmax,nmodmax)
     +       ,seriesscen(nensmax,nmodmax,9)
        real seriesmean(9,nscenmax),seriesquant(-6:6,9,nscenmax),
     +       seriesdecvar(9,nscenmax),mean(12)
        real s,s1,s2
        logical lwrite,lexist,standardunits
        character scenarios(nscenmax)*7
        character models(nmodmax)*100,shortmodel*20
        character string*128,format*20,vars*80,units*30,file*255
     +       ,line*255,var*40,region_subregion*100,season*10
     +       ,oneall*3,inpattern*1023,quantfile*1023,spid*10
     +       ,rcplist*30,command*1023
        integer getpid
        data scenarios /'rcp26','rcp45','rcp60','rcp85','sresa1b'/
        lwrite = .false.
        call getenv('LWRITE',string)
        call tolower(string)
        if ( string.eq.'true' ) lwrite = .true.
!
!       find list of model files
!
        call getarg(1,var)
        call getarg(2,region_subregion)
        call getarg(3,oneall)
        if ( oneall.ne.'one' .and. oneall.ne.'all' ) then
            write(0,*) 'usage: quantiles_series var region_subregion '//
     +           'one|all mon N ave M yr1r yr2r yr1s yr2s '//
     +           'rcplist inpattern quantfile'
            write(0,*) 'yr1r=yr2r=0 signifies no anomalies'
            call abort
        end if
        call getarg(4,string)
        if ( string(1:3).ne.'mon' ) then
                write(0,*) 'error: expecting mon'
                call abort
        end if
        call getarg(5,string)
        read(string,*) mon
        if ( mon.lt.1 .or. mon.gt.12 ) then
                write(0,*) 'error: expecting 1 <= mon <= 12, not ',mon
                call abort
        end if
        call getarg(6,string)
        if ( string(1:3).ne.'ave' ) then
                write(0,*) 'error: expecting ave'
                call abort
        end if
        call getarg(7,string)
        read(string,*) ave
        if ( ave.lt.1 .or. ave.gt.12 ) then
                write(0,*) 'error: expecting 1 <= ave <= 12, not ',ave
                call abort
        end if
        call getarg(8,string)
        read(string,*) yr1r
        call getarg(9,string)
        read(string,*) yr2r
        call getarg(10,string)
        read(string,*) yr1s
        call getarg(11,string)
        read(string,*) yr2s
        if ( (yr1r.lt.1800 .and. yr1r.ne.0) .or. yr1s.lt.1800 .or.
     +       (yr2r.gt.2300 .and. yr2r.ne.0) .or. yr2s.gt.2300 )
     +       then
            write(0,*) 'quantiles_series: error: yrs out of range: ',
     +           yr1r,yr2r,yr1s,yr2s
            call abort
        end if
        if ( yr1r.gt.yr2r ) then
            write(0,*) 'quantiles_series: error: yr1r.gt.yr2r: ',
     +           yr1r,yr2r
            call abort
        end if
        if ( yr1s.gt.yr2s ) then
            write(0,*) 'quantiles_series: error: yr1s.gt.yr2s: ',
     +           yr1s,yr2s
            call abort
        end if
        call getarg(12,rcplist)
        call getarg(13,inpattern)
        call getarg(14,quantfile)

        pid = getpid()
        write(spid,'(i10.10)') pid

        write(season,'(a,i2.2,a,i2.2)') 'mon',mon,'ave',ave
                nperiod = 1
                if ( mon+ave-1.gt.12 ) then
                        offset = -1
                else
                offset = 0
        end if

        do iscen=1,nscenmax
            if ( lwrite ) print *,scenarios(iscen),' ',rcplist,
     +           index(rcplist,scenarios(iscen))
            if ( index(rcplist,trim(scenarios(iscen))).eq.0 ) cycle
            print *,trim(scenarios(iscen))
            if ( iscen.le.4 ) then
                if ( oneall.eq.'all' ) then
                    command = 
     +                   'ls '//trim(inpattern)//
     +                   ' | fgrep '//trim(scenarios(iscen))//
     +                   ' > /tmp/list_of_all_files_'//
     +                   trim(scenarios(iscen))//'_'//trim(spid)//
     +                   '.txt'
                    call mysystem(trim(command),iret)
                else
                    command = 
     +                   'ls '//trim(inpattern)//
     +                   ' | fgrep '//trim(scenarios(iscen))//
     +                   ' | fgrep -v EC-EARTH | '//
     +                   'fgrep -v HadGEM2-ES > /tmp/list_of_one_files_'
     +                   //trim(scenarios(iscen))//'_'//trim(spid)//
     +                   '.txt'
                    !!!print *,trim(command)
                    call mysystem(trim(command),iret)
!!! SPECIAL CASES, UPDATE WHEN MORE ENSEMBLE MEMBERS BECOME AVAILABLE
                    if ( trim(scenarios(iscen)).ne.'rcp60' ) then
                        command = 
     +                   'ls '//trim(inpattern)//
     +                   ' | fgrep '//trim(scenarios(iscen))//
     +                   ' | fgrep EC-EARTH '//
     +                   '| fgrep r8i1p1 >> /tmp/list_of_one_files_'
     +                   //trim(scenarios(iscen))//'_'//trim(spid)//
     +                   '.txt'
                    !!!print *,trim(command)
                    end if
                    call mysystem(trim(command),iret)
                    command = 
     +                   'ls '//trim(inpattern)//
     +                   ' | fgrep '//trim(scenarios(iscen))//
     +                   ' | fgrep HadGEM2-ES '//
     +                   '| fgrep r2i1p1 >> /tmp/list_of_one_files_'
     +                   //trim(scenarios(iscen))//'_'//trim(spid)//
     +                   '.txt'
                    !!!print *,trim(command)
                    call mysystem(trim(command),iret)
                end if
            else
                command = 
     +               'ls '//trim(inpattern)//
     +               ' | fgrep '//trim(scenarios(iscen))//
     +               ' > /tmp/list_of_'//oneall//'_files_'//
     +               trim(scenarios(iscen))//'_'//trim(spid)//
     +               '.txt'
                call mysystem(trim(command),iret)
            end if
            if ( iret.ne.0 ) then
                write(0,*) 'error in system call ',iret
                !!!call abort
            end if
            nmod = 0
            nens = 0
            iens = 0
            open(1,file='/tmp/list_of_'//oneall//'_files_'
     +           //trim(scenarios(iscen))//'_'//trim(spid)//
     +           '.txt',status='old')
 1          continue
            read(1,'(a)',end=2,err=2) line
            if ( index(line,'modmean').ne.0 .or.
     +           index(line,'onemean').ne.0 .or.
     +           index(line,'_ave').ne.0 ) goto 1
            ifirst = index(inpattern,'*')
            if ( index(inpattern,'CORDEX').eq.0 ) then
                if ( iscen.le.4 ) then
                    i = ifirst + index(line(ifirst:),'_') - 2
                else
                    i = ifirst + index(line(ifirst:),'_144') - 2
                    j = i
                    if ( line(i-2:i-2).eq.'_' .and.
     +                   ichar(line(i:i)).ge.ichar('0') .and.
     +                   ichar(line(i:i)).le.ichar('9') .and.
     +                   ichar(line(i-1:i-1)).ge.ichar('0') .and.
     +                   ichar(line(i-1:i-1)).le.ichar('9') ) then
                        i = i-3
                    end if
                end if
                string = line(ifirst:i)
                iens = iens + 1
!!! YET ANOTHER SPECIAL CASE; I ASSUME ONLY GISS HAS DIFFERENT p's
                if ( string(1:4).eq.'GISS' .and. index(string,'CC').eq.0
     +               ) then
                    i = index(line,'i1p') + 3
                    string = trim(string)//'_p'//line(i:i)
                end if
                do imod=1,nmod
                    if ( trim(models(imod)).eq.trim(string) ) then
                        nens(imod) = nens(imod) + 1
                        call getensnumber(line,ifirst,
     +                       iiens(nens(imod),imod),lwrite)
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
            if ( index(inpattern,'CORDEX').eq.0 ) then
                call getensnumber(line,ifirst,iiens(nens(nmod),nmod),
     +               lwrite)
                if ( nmod.gt.1 .and. lwrite ) print *,
     +               '         with ',nens(nmod-1),' realisations ',iens
            end if
            if ( lwrite ) print *,'Found model ',nmod,' ',
     +           trim(models(nmod))
            goto 1
 2          continue
            close(1,status='delete')
            if ( nmod.gt.1 ) then
                if ( lwrite) print *,'         with ',nens(nmod-1),
     +               ' realisations ',iens
            end if
            if ( oneall.eq.'all' .and. index(inpattern,'CORDEX').eq.0 )
     +           then
                print *,'Found ',nmod,' models with in all ',iens,
     +               ' ensemble members<br>'
            else
                print *,'Found ',nmod,' models, using one ensemble ',
     +               ' member per model<br>'
            end if
!
!           read data
!
            series = 3e33
            do imod=1,nmod
                do iens=1,nens(imod)
                    if ( oneall.eq.'one') then
                        jens = 1
                        if ( models(imod).eq.'EC-EARTH' ) then
                            if ( iiens(iens,imod).ne.8 ) cycle
                        else if ( models(imod).eq.'HadGEM2-ES' ) then
                            if ( iiens(iens,imod).ne.2 ) cycle
                        else
                            if ( iiens(iens,imod).ne.1 ) cycle
                        end if
                    else
                        jens = iens
                    end if
                    if ( iiens(iens,imod).lt.10 ) then
                        format = '(5a,i1,a,i1,5a)'
                    else
                        format= '(5a,i2,a,i1,5a)'
                    end if
                    if ( index(inpattern,'CORDEX').ne.0 ) then
                        write(file,'(5a)') inpattern(:ifirst-1)
     +                       ,trim(models(imod)),
     +                       '_v1_mon_195101-209912_latlon_',
     +                       trim(region_subregion),'.dat'
                        i = index(file,'rcp')
                        file(i:i+4) = trim(scenarios(iscen)) !! quick dirty hack
                        inquire(file=file,exist=lexist)
                        if ( .not.lexist ) then
                            ! some MOHC-forced runs only go to november
                            i = index(file,'209912')
                            file(i:i+5) = '209911'
                        end if
                    else if ( iscen.le.4 ) then
                        shortmodel = models(imod)
                        i = index(shortmodel,'_p')
                        if ( i.ne.0 ) then
                            read(shortmodel(i+2:i+2),'(i1)') p
                            shortmodel = shortmodel(1:i-1)
                        else
                            p = 1
                        end if
                        
                        write(file,format) inpattern(:ifirst-1)
     +                       ,trim(shortmodel),'_'
     +                       ,trim(scenarios(iscen)),'_r',iiens(iens
     +                       ,imod),'i1p',p,'_',trim(region_subregion)
     +                       ,'.dat'
                        i = index(file,'rcp')
                        file(i:i+4) = trim(scenarios(iscen)) !! quick dirty hack
                    else
                        write(file,'(3a,i2.2,3a)') inpattern(:ifirst-1)
     +                       ,trim(models(imod))
     +                       ,'_',iiens(iens,imod)-1,'_144_'
     +                       ,trim(region_subregion),'.dat'
                        inquire(file=trim(file),exist=lexist)
                        if ( .not.lexist ) then
                            write(file,'(5a)') inpattern(:ifirst-1)
     +                           ,trim(models(imod)),'_144_'
     +                           ,trim(region_subregion),'.dat'
                        end if
                    end if
                    if ( lwrite ) print *,'reading ',trim(file),'<br>'
                    standardunits = .true.
                    call readseries(file,series(1,yrbeg,jens,imod),12
     +                   ,yrbeg,yrend,nperyear,vars,units,standardunits
     +                   ,.false.)
                    if ( index(var,'rel ').gt.0 .and. yr1r.gt.0 ) then
                        call ensanomalclim(series(1,yrbeg,jens,imod),12
     +                       ,nperyear,yrbeg,yrend,0,0,yr1r,yr2r,mean)
                        call takerelanom(series(1,yrbeg,jens,imod),mean,
     +                       12,yrbeg,yrend,0,0,nperyear,mon,ave,lwrite)
                    end if
                end do
            end do
            if ( oneall.eq.'one' ) then
                do imod=1,nmod
                    nens(imod) = 1
                end do
            end if
!
!           compute means
!
            do imod=1,nmod
                do iens=1,nens(imod)
                    call sumit(series(1,yrbeg,iens,imod),12,nperyear
     +                   ,yrbeg,yrend,ave,'v')
                end do
            end do
            seriesclim = 0
            seriesscen = 0
            iperiod = 1
            do imod=1,nmod
                do iens=1,nens(imod)
                    n = 0
                    if ( yr1r.gt.0 ) then
                        do yr=yr1r+offset,yr2r+offset
                            if ( series(mon,yr,iens,imod).lt.1e30 ) then
                                n = n + 1
                                seriesclim(iens,imod) =
     +                               seriesclim(iens,imod)
     +                               + series(mon,yr,iens,imod)
                            end if
                        end do
                        if ( n.eq.0 ) then
                            print *,'weird, n=0 for ',models(imod)
                            do yr=yr1r+offset,yr2r+offset
                                print *,'series(',mon,yr,iens,imod,
     +                               ') = ',series(mon,yr,iens,imod)
                            end do
                        else
                            seriesclim(iens,imod) = 
     +                           seriesclim(iens,imod)/n
                        end if
                    else
                        seriesclim(iens,imod) = 0
                    end if
                    n = 0
                    do yr=yr1s+offset,yr2s+offset
                        if ( series(mon,yr,iens,imod).lt.1e30 )
     +                       then
                            n = n + 1
                            seriesscen(iens,imod,iperiod) =
     +                           seriesscen(iens,imod,iperiod)
     +                           + series(mon,yr,iens,imod)
                        end if
                    end do
                    if ( n.eq.0 ) then
                        if ( lwrite ) print *,'weird, n = 0'
                         seriesscen(iens,imod,iperiod) = 3e33
                    else
                         seriesscen(iens,imod,iperiod) = 
     +                        seriesscen(iens,imod,iperiod)/n
     +                        - seriesclim(iens,imod)
                     end if
                end do
            end do
!
!           compute s.d. of 20-yr means from inter-ensemble spread
!
            do iperiod=nperiod,1,-1
                seriesdecvar(iperiod,iscen) = 0
                m = 0
                do imod=1,nmod
                    if ( nens(imod).gt.2 ) then
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
                            s2 = s2 +
     +                           (seriesscen(iens,imod,iperiod)-s1)**2
                        end do
                        seriesdecvar(iperiod,iscen) =
     +                       seriesdecvar(iperiod,iscen) + s2/(n-1)
                    end if
                end do
                if ( m.gt.0 ) then
                        seriesdecvar(iperiod,iscen) =
     +                  seriesdecvar(iperiod,iscen)/m
                        else
                                seriesdecvar(iperiod,iscen) = 3e33
                        end if
                if ( lwrite ) print *,'seriesdecvar(',iperiod,iscen,
     +               ') = ',seriesdecvar(iperiod,iscen)
            end do
!
!           compute mean & quantiles
!
            do iperiod=nperiod,1,-1
                call gethistmean(seriesmean(iperiod,iscen),
     +               seriesscen(1,1,iperiod),nens,nmod,nensmax,lwrite)
                call gethistquant(seriesquant(-6,iperiod,iscen),
     +               seriesscen(1,1,iperiod),nens,nmod,nensmax,lwrite)
            end do
        end do                  ! scenario

        open(1,file=trim(quantfile))
        write(1,'(10a)') '# generated by quantiles_series ',trim(var)
     +       ,' ',trim(region_subregion),' ',trim(season),' '
     +       ,oneall
        write(1,'(2a)') '# year mean min 2.5% 5% 10% 17% 25% ',
     +       '50% 75% 83% 90% 95% 97.5% max sd'
        do iscen=1,nscenmax
            write(1,'(2a)') '# ',scenarios(iscen)
            do iperiod=nperiod,1,-1
                write(1,'(i4,50f12.4)') yr2s,
     +               seriesmean(iperiod,iscen),
     +               (seriesquant(i,iperiod,iscen),i=-6,6)
     +               ,sqrt(seriesdecvar(1,iscen))
            end do
            write(1,'(a)')
            write(1,'(a)')
        end do
        end

        subroutine getensnumber(string,ifirst,n,lwrite)
        implicit none
        integer ifirst,n
        character*(*) string
        logical lwrite
        integer i,j
        i = index(string,'sresa1b')
        if ( i.ne.0 ) then      ! SRES A1b runs
            i = ifirst + index(string(ifirst:),'_144') - 2
            j = i
            if ( string(i-2:i-2).eq.'_' .and.
     +           ichar(string(i:i)).ge.ichar('0') .and.
     +           ichar(string(i:i)).le.ichar('9') .and.
     +           ichar(string(i-1:i-1)).ge.ichar('0') .and.
     +           ichar(string(i-1:i-1)).le.ichar('9') ) then
                i = i-3
            end if
            if ( j.ne.i ) then
                read(string(i+2:j),*) n
            else
                n = 0
            end if
            n = n + 1
       else                    ! RCP runs
            i = index(string,'_r',.true.)
            if ( i.eq.0 ) then
                write(0,*) 'error: cannot find "_r" in ',trim(string)
                n = 0
                return
            end if
            i = i+2
            j = i + index(string(i:),'i') - 2
!!!            print *,'searching for ensemble number in ',trim(string),
!!!     +           i,j,string(i:j)
            if ( j.lt.i ) then
                write(0,*) 'getensnumber: error: cannot find "i" in ',
     +               trim(string(i:))
                call abort
            end if
            read(string(i:j),*) n
        end if
        if ( lwrite ) print *,'found ensemble number ',n,' in '
     +       ,trim(string)
        end

        subroutine gethistmean(mean,data,nens,nmod,nensmax,lwrite)
!
!       Compute the mean of a set of values.
!       All models are weighted equally
!
        implicit none
        integer nmod,nensmax,nmodmax
        integer nens(nmod)
        real mean,data(nensmax,nmod)
        logical lwrite
        integer iens,imod
        real s

        if ( lwrite ) then
            print *,'gethistmean: input'
            do imod=1,nmod
                print '(i3,100f7.3)',
     +               imod,(data(iens,imod),iens=1,nens(imod))
            end do
        end if
        s = 0
        do imod=1,nmod
            do iens=1,nens(imod)
                if ( data(iens,imod).gt.100 ) then
                    write(0,*) 'gethistmean: warning: data(',iens,imod,
     +                   ') = ',data(iens,imod)
                end if
                s = s + data(iens,imod)/nens(imod)
            end do
        end do
        mean = s/nmod
        if ( lwrite ) print '(a,f7.3)','mean = ',mean
        end

        subroutine gethistquant(quant,data,nens,nmod,nensmax,lwrite)
!
!       Compute some qunatiles of a set of values:
!       min, 5%, 25%, 50%, 75%, 95%, max
!       All models are weighted equally
!
        implicit none
        integer nmod,nensmax,nmodmax
        integer nens(nmod)
        real quant(-6:6),data(nensmax,nmod)
        logical lwrite
        integer iens,imod,i,n,nnens(1024)
        real point(1024)
!
!       just convert the data to the same format as quantiles_field uses
!
        n = 0
        do imod=1,nmod
            do iens=1,nens(imod)
                n = n + 1
                point(n) = data(iens,imod)
                nnens(n) = nens(imod)
            end do
        end do

        call getweightedquant(point,nnens,n,nmod,quant,lwrite)

        end
