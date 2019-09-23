program getmomentsfield
!
!       compute the first four moments of a field, or a percentile, or
!       the min or max, and write them out as GrADS files
!       9-nov-2006 added GPD fits
!
        use lsdata
        implicit none
        include 'params.h'
        include 'netcdf.inc'
        include 'getopts.inc'
        integer,parameter :: recfa4=4,recfac=4
        real,parameter :: absent=3e33
        integer,parameter :: nvarmax=14,ntmax=13,nyears=1000
        integer nx,ny,nz,nt,firstyr,lastyr,firstmo,nvars, &
     &       ivars(2,nvarmax),jvars(6,nvarmax),ncid,endian, &
     &       status,nperyear,mens,mens1,years(nyears),iyear
        integer jx,jy,jz,i,j,jj,j1,j2,k,m,n,month,yr,imoment,ldir, &
     &       iens,yrstart,yrstop,f,irec,year,ntvarid,itimeaxis(ntmax), &
     &       itype,iyeartype
        logical lexist
        real xx(nxmax),yy(nymax),zz(nzmax),undef,xxls(nxmax),yyls(nymax) &
     &       ,a,b,xi,t(10),t25(10),t975(10),tx,tx1,tx25,tx975,sx,sy,s,xtreme
        real xmom(5),var,perc,z
        real,allocatable :: field(:,:,:,:,:,:),res(:,:,:,:,:)
        real,allocatable :: fxy(:,:,:),ddata(:)
        character infile*255,datfile*255,outfile*255,line*255,lz(3)*20 &
     &       ,vars(nvarmax)*15,lvars(nvarmax)*255,svars(nvarmax)*80 &
     &       ,title*255,history*20000,units(nvarmax)*20, &
     &       cell_methods(nvarmax)*100,lsmasktype*4,ltime*120,metadata(2,100)*2000
        character yesno*1,dir*255,string*15,saveunits*20,format*10,assume*5
        integer rindex
!
!       process command line
!
        n = command_argument_count()
        if ( n.lt.3 ) then
            write(0,*) 'usage: getmomentsfield infile.[ctl|nc] [1234] ' &
     &            //'[month m[:n] [sum|ave|max|min|sel n] [log|sqrt] ' &
     &            //'[minfac r] [minnum n] [begin yr] [end yr] ' &
     &            //'[lt cut] [gt cut] [diff [nyr]] [detrend] ' &
     &            //'outfile.[ctl|nc]'
            call exit(-1)
        end if
        lwrite = .false.
        do i=1,n
            call get_command_argument(i,line)
            if ( line.eq.'debug' .or. line.eq.'lwrite' ) then
                lwrite = .true.
            end if
        end do
        call get_command_argument(1,infile)
        call getmetadata(infile,mens1,mens,ncid,datfile,nxmax,nx &
     &       ,xx,nymax,ny,yy,nzmax,nz,zz,lz,nt,nperyear,firstyr,firstmo &
     &       ,ltime,undef,endian,title,history,nvarmax,nvars,vars,jvars &
     &       ,lvars,svars,units,cell_methods,metadata,lwrite)
        call getlastyr(firstyr,firstmo,nt,nperyear,lastyr)
!       co-ordinate with adjustunits, adjustvar
!       this catches the time of year of min, max occurrence
        itype = 0
        if ( vars(1)(1:5).eq.'time_' .and. nperyear.eq.1 ) then
            if ( units(1).eq.'dy' ) then
                itype = 366
            else if ( units(1).eq.'mo' ) then
                itype = 12
            else if ( units(1).eq.'season' ) then
                itype = 4
            else if ( units(1).eq.'period' ) then
                write(0,*) 'getomomentsfield: error: cannot determine ', &
     &               'length of period'
                call exit(-1)
            end if
            call add_metadata('method','getmomentsfield: mean was taken over a circle',metadata)
            iyeartype = 1 ! the default year is defined Jan-Dec
            i = index(history,'daily2longer')
            if ( i > 0 ) then
                ! search for the second space after the command
                i = i + index(history(i:),' ')
                i = i + index(history(i:),' ')
                if ( history(i:i+1) == '-1' ) iyeartype = -1 ! the year was defined from July to June
            end if
        end if
!       process arguments
        call get_command_argument(2,line)
        if ( line(1:4).eq.'perc' .or. line(1:6).eq.'pot_rt' .or. &
     &       line(1:6).eq.'gev_rt' .or. line(1:7).eq.'norm_rt' .or. &
     &       line(1:6).eq.'norm_z' .or. line(1:4).eq.'rank' ) then
            i = 4
        else
            i = 3
        end if
        call getlsmask(i,lsmasktype,nxmax,xxls,nymax,yyls,lwrite)
        if ( lsmasktype.ne.'all' ) then
            call checkgridequal(nx,ny,xx,yy,nxls,nyls,xxls,yyls)
        end if
        call getopts(i,n-1,nperyear,yrbeg,yrend,.true.,mens1,mens)
        call get_command_argument(1,infile)
        if ( index(infile,'%').gt.0 .or. &
     &       index(infile,'++').gt.0 ) then
            write(0,*) 'Using ensemble members ',nens1,' to ',nens2 &
     &           ,'<br>'
        end if
        yr1 = max(yr1,firstyr)
        yr2 = min(yr2,lastyr)
        firstyr = yr1
        lastyr = yr2
!
!       allocate fields
!
        if ( lwrite ) print *,'allocating ',nx*ny*nz*nperyear* &
     &       (lastyr-firstyr+1)*(nens2-nens1+1),' reals'
        allocate(field(nx,ny,nz,nperyear,firstyr:lastyr,nens1:nens2))
        allocate(res(nx,ny,nz,0:12,nvarmax))
        allocate(fxy(nperyear,yr1:yr2,0:nens2))
        allocate(ddata(1+nperyear*(yr2-yr1+1)*(nens2-nens1+1)))
        call get_command_argument(2,line)
        year = 0
        if (  ichar(line(1:1)).ge.ichar('0') .and. &
     &        ichar(line(1:1)).le.ichar('9') ) then
            read(line,*,err=901) imoment
            if ( imoment.lt.-2 .or. imoment.gt.4 ) goto 901
        elseif ( line(1:4).eq.'adev' ) then
            imoment = 0
        elseif ( line(1:3).eq.'ave' .or. line(1:4).eq.'mean' ) then
            imoment = 1
        elseif ( line(1:3).eq.'sdm' .or. line(1:6).eq.'s.d./m' ) then
            imoment = -2
        elseif ( line(1:2).eq.'sd' .or. line(1:4).eq.'s.d.' .or. &
     &           line(1:4).eq.'norm' ) then
            imoment = 2
            if ( line(1:7).eq.'norm_rt' .or. line(1:6).eq.'norm_z' ) &
     &           then
                call get_command_argument(3,line)
                read(line,*,err=905) year
            end if
        elseif ( line(1:4).eq.'skew' ) then
            imoment = 3
        elseif ( line(1:4).eq.'kurt' .or. line(1:4).eq.'curt' ) then
            imoment = 4
        elseif ( line(1:4).eq.'perc' ) then
            imoment = -1
            call get_command_argument(3,line)
            read(line,*,err=902) perc
        elseif ( line(1:3).eq.'max' ) then
            imoment = +100
        elseif ( line(1:3).eq.'min' ) then
            imoment = -100
        elseif ( line(1:3).eq.'pot' ) then
            if ( line(1:6).eq.'pot_rt' ) then
                call get_command_argument(3,line)
                read(line,*,err=905) year
            end if
            assume = 'shift'
            imoment = 200
        elseif ( line(1:3).eq.'gev' ) then
            if ( line(1:6).eq.'gev_rt' ) then
                call get_command_argument(3,line)
                read(line,*,err=905) year
            end if
            imoment = 300
        elseif ( line(1:4).eq.'rank' ) then
            call get_command_argument(3,line)
            read(line,*,err=905) year
            imoment = 1000
        elseif ( line(1:5).eq.'timex' ) then
            imoment = 1001
        else
            goto 901
        end if
        if ( lag1.ne.0 .or. lag2.ne.0 ) print * &
     &        ,'getmomentsfield: lags do not make sense'
        if ( dump ) print *,'getmomentsfield: dump not supported'
        if ( plot ) print *,'getmomentsfield: plot not supported'
        if ( lks ) print *,'getmomentsfield: K-S not supported'
        if ( lconting ) print *,'getmomentsfield: contingency '// &
     &        'tables not supported'
        do i=1,indxuse
            if ( lincl(i) ) print *,'getmomentsfield: what do ', &
     &          'you mean with ',strindx(i),'?'
        end do
        if ( itype.ne.0 .and. imoment.ne.1 ) then
            write(0,*) 'getmomemntsfield: error: can only compute ', &
     &           'mean of cyclical data, not ',trim(line)
            call exit(-1)
        end if
!       range of years
        yr1 = max(yr1,firstyr)
        yr2 = min(yr2,firstyr + (firstmo+nt-2)/nperyear)
        yrstart = yr2
        yrstop  = yr1
!
!       read field, change absent values to our convention
!
        do iens=nens1,nens2
            call keepalive1('Reading ensemble member', &
     &           iens-nens1+1,nens2-nens1+1)
            if ( ncid.eq.-1 ) then
                call get_command_argument(1,infile)
                if ( index(infile,'%').gt.0 .or. &
     &               index(infile,'++').gt.0 ) then
                    call filloutens(infile,iens)
                    if ( lwrite ) print *,'calling parsectl on ', &
     &                   trim(infile)
                end if
                call parsectl(infile,datfile,nxmax,nx,xx,nymax,ny,yy &
     &               ,nzmax,nz,zz,nt,nperyear,f,firstmo,undef &
     &               ,endian,title,1,nvars,vars,ivars,lvars,units)
                if (lwrite) print '(2a)','# looking for ',trim(datfile)
                inquire(file=datfile,exist=lexist)
                if ( .not.lexist ) then
                    if (lwrite) print '(3a)','# looking for ' &
     &                   ,trim(datfile),'.gz'
                    inquire(file=trim(datfile)//'.gz',exist=lexist)
                    if ( .not.lexist ) then
                        nens2 = iens-1
                        if ( nens2.ge.nens1 ) then
                            write(0,*) 'Found ensemble 0 to ',nens2,'<br>'
                            goto 5
                        else
                            write(0,*) 'Cannot locate file ',trim(datfile)
                            call exit(-1)
                        end if
                    end if
                end if
                if ( lwrite ) then
                    print *,'opening file ',trim(datfile)
                end if
                call zreaddatfile(datfile,field(1,1,1,1,firstyr,iens), &
     &                nx,ny,nz,nx,ny,nz,nperyear,firstyr,lastyr, &
     &                f,firstmo,nt,undef,endian,lwrite,yr1,yr2,1,1 &
     &                )
            else
                call get_command_argument(1,infile)
                if ( index(infile,'%').gt.0 .or. &
     &               index(infile,'++').gt.0 ) then
                    call filloutens(infile,iens)
                    if ( lwrite ) print *,'calling parsenc on ',trim(infile)
                    status = nf_open(infile,nf_nowrite,ncid)
                end if
                call parsenc(infile,ncid,nxmax,nx,xx,nymax,ny,yy &
     &               ,nzmax,nz,zz,nt,nperyear,f,firstmo &
     &               ,undef,title,1,nvars,vars,jvars,lvars,units)
                call zreadncfile(ncid,field(1,1,1,1,firstyr,iens) &
     &               ,nx,ny,nz,nx,ny,nz,nperyear,firstyr,lastyr,f &
     &               ,firstmo,nt,undef,lwrite,yr1,yr2,jvars)
            end if
        end do
    5   continue
!
!       apply land/sea mask
!
        call applylsmask(field,lsmask,nx,ny,nz,nperyear,firstyr,lastyr, &
     &       nens1,nens2,lsmasktype,lwrite)
!
!       open output file
!
        call get_command_argument(command_argument_count(),outfile)
        inquire(file=outfile,exist=lexist)
        if ( lexist ) then
            print *,'output file ',outfile(1:index(outfile,' ')-1), &
     &            ' already exists, overwrite? [y/n]'
            read(*,'(a)') yesno
            if (  yesno.ne.'y' .and. yesno.ne.'Y' .and. &
     &            yesno.ne.'j' .and. yesno.ne.'J' ) then
                stop
            end if            
            open(2,file=outfile)
            close(2,status='delete')
        end if
        if ( index(outfile,'.ctl').ne.0 ) then
            i = index(outfile,'.ctl')
            if ( i.ne.0 ) then
                datfile = outfile(:i-1)//'.grd'
            else
                datfile = outfile
            end if
            open(unit=2,file=datfile,form='unformatted',access='direct' &
     &            ,recl=recfac*nx*ny*nz,err=920)
        end if
!
!       compute minfac if it has not been set explicitly
!
        if ( minfac.lt.0 .and. minnum.lt.0 ) then
!           heuristic, gives 0.25 for 150 yrs, 0.5 for 50 yrs, 0.75 for 20yrs
            minfac = max(0.1, &
     &            min(0.6, &
     &            1.5-log(1+real(min(nt,nperyear*(yr2-yr1+1))-1) &
     &            /nperyear)/4))
        end if
        write(0,'(a,i2,a)') 'Requiring at least ', &
     &            nint(100*minfac),'% valid points'
!
!       loop over gridpoints
!
        do jz=1,nz
            do jy=1,ny
                call keepalive1('Computing latitude ', &
     &               jy+(jz-1)*ny,ny*nz)
                do jx=1,nx
                    do month=0,min(12,nperyear)
                        res(jx,jy,jz,month,1:nvarmax) = 3e33
                    end do
                    xyear = 3e33
!
!                   create 1-D series from field
!
                    n = 0
                    do iens=nens1,nens2
                        do i=yr1,yr2
                            do j=1,nperyear
                                fxy(j,i,iens) = &
     &                                  field(jx,jy,jz,j,i,iens)
                                if ( fxy(j,i,iens).lt.0.9*absent ) &
     &                                  n = n+1
                            end do
                        end do
                    end do
                    if ( n.lt.3 ) then
                        if ( lwrite ) print '(a,3i5)', &
     &                          'no valid points at ',jx,jy,jz
                        goto 800
                    end if
                    do iens=nens1,nens2
!
!                       sum
!
                        if ( lsum.gt.1 ) then
                            call sumit(fxy(1,yr1,iens),nperyear, &
     &                          nperyear,yr1,yr2,lsum,oper)
                        end if
!
!                       log,sqrt
!
                        if ( logscale ) then
                            do i=yr1,yr2
                                do j=1,nperyear
                                    if ( fxy(j,i,iens).lt.1e33 .and. &
     &                                   fxy(j,i,iens).gt.0 ) &
     &                                  then
                                        fxy(j,i,iens) = &
     &                                          log10(fxy(j,i,iens))
                                    else
                                        fxy(j,i,iens) = 3e33
                                    end if
                                end do
                            end do
                        end if
                        if ( sqrtscale ) then
                            do i=yr1,yr2
                                do j=1,nperyear
                                    if ( fxy(j,i,iens).lt.1e33 .and. &
     &                                   fxy(j,i,iens).ge.0 ) &
     &                                  then
                                        fxy(j,i,iens) = &
     &                                          sqrt(fxy(j,i,iens))
                                    else
                                        fxy(j,i,iens) = 3e33
                                    end if
                                end do
                            end do
                        end if
                        if ( squarescale ) then
                            do i=yr1,yr2
                                do j=1,nperyear
                                    if ( fxy(j,i,iens).lt.1e16 ) then
                                        fxy(j,i,iens) = &
     &                                          fxy(j,i,iens)**2
                                    else
                                        fxy(j,i,iens) = 3e33
                                    end if
                                end do
                            end do
                        end if
                        if ( cubescale ) then
                            do i=yr1,yr2
                                do j=1,nperyear
                                    if ( fxy(j,i,iens).lt.1e11 ) then
                                        fxy(j,i,iens) = &
     &                                          fxy(j,i,iens)**3
                                    else
                                        fxy(j,i,iens) = 3e33
                                    end if
                                end do
                            end do
                        end if
                        if ( twothirdscale ) then
                            do i=yr1,yr2
                                do j=1,nperyear
                                    if ( fxy(j,i,iens).lt.1e33 .and. &
     &                                   fxy(j,i,iens).ge.0 ) &
     &                                  then
                                        fxy(j,i,iens) = &
     &                                          fxy(j,i,iens)**(2./3.)
                                    else
                                        fxy(j,i,iens) = 3e33
                                    end if
                                end do
                            end do
                        end if
!
!                       detrend
!
                        if ( ldetrend ) then
                            if ( lwrite ) print *,'Detrending field'
                            call detrend(fxy(1,yr1,iens),nperyear, &
     &                           nperyear,yr1,yr2,yr1,yr2,m1,m2,lsel &
     &                           )
                        end if
!
!                       differentiate
!
                        if ( ndiff.ne.0 ) then
                            if ( lwrite ) print *,'Taking differences'
                            call diffit(fxy(1,yr1,iens),nperyear, &
     &                          nperyear,yr1,yr2,ndiff)
                        end if
!
!                       anomalies
!
                        if ( anom .or. lsel.gt.1 .and. ndiff.le.0 ) &
     &                      then
                            call anomal(fxy(1,yr1,iens),nperyear, &
     &                          nperyear,yr1,yr2,yr1,yr2)
                        end if
                    end do
!
!                   anomalies wrt ensemble mean
!
                    if ( nens2.gt.nens1 .and. lensanom ) then
                        call anomalensemble(fxy,nperyear,nperyear, &
     &                       yr1,yr2,yr1,yr2,nens1,nens2)
                    end if
!
!                   normalize to s.d.
!
                    do iens=nens1,nens2
                        if ( lnormsd ) then
                            call normsd(fxy(1,yr1,iens),nperyear, &
     &                          nperyear,yr1,yr2,yr1,yr2)
                        end if
                    end do       ! iens
!
!                   get moments or other properties
!
                    do month=m1,m2
                        call getj1j2(j1,j2,month,nperyear,.false.)
!
!                       fill linear arrays without absent values
!                       and compute moment
!
                        n = 0
                        do iens=nens1,nens2
                            do yr=yr1-1,yr2
                                do jj=j1,j2
                                    j = jj
                                    call normon(j,yr,i,nperyear)
                                    if ( yr.eq.year ) then
                                        if ( nens1.ne.nens2 ) then
                                            write(0,*)'getmomentsfield' &
     &                                           //': error: cannot ' &
     &                                           //'handle ensembles'
                                            write(*,*)'getmomentsfield' &
     &                                           //': error: cannot ' &
     &                                           //'handle ensembles'
                                            call exit(-1)
                                        end if
                                        if ( j1.ne.j2 ) then
                                            write(0,*)'getmomentsfield' &
     &                                           //': error: can only ' &
     &                                           //'handle annual data'
                                            write(*,*)'getmomentsfield' &
     &                                           //': error: can only ' &
     &                                           //'handle annual data'
                                            call exit(-1)
                                        end if
                                        xyear = fxy(j,i,iens)
                                        fxy(j,i,iens) = absent
                                    end if
                                    if ( i.lt.yr1 .or.i.gt.yr2 ) &
     &                                  goto 710
                                    if (  fxy(j,i,iens).lt.absent/3.and. &
     &                                    fxy(j,i,iens).lt.maxindx .and. &
     &                                    fxy(j,i,iens).gt.minindx) then
                                        n = n+1
                                        ddata(n) = fxy(j,i,iens)
                                        yrstart = min(yrstart,i)
                                        yrstop  = max(yrstop,i)
                                    end if
  710                               continue
                                end do
                            end do
                        end do
                        if ( month.eq.0 .and. &
     &                        n.lt.minfac*min(nt,nperyear*(yr2-yr1+1)) &
     &                        .or. &
     &                        month.ne.0 .and. &
     &                        n.lt.minfac*min(nt/nperyear,yr2-yr1+1) &
     &                        .or. &
     &                        n.lt.minnum ) then
                            if ( lwrite ) print '(a,3i5,i3,a,2i6)' &
     &                          ,'not enough valid points at ',jx,jy,jz &
     &                          ,month,': ',n,nt
                            goto 790
                        end if
!
                        m = month-m1
                        if ( itype.gt.0 ) then
                            ! remember, only mean is supported
                            nvars = 1
                            sx = 0
                            sy = 0
                            do i=1,n
                                sx = sx + cos(8*atan(1.)*ddata(i)/itype)
                                sy = sy + sin(8*atan(1.)*ddata(i)/itype)
                            end do
                            s = itype*atan2(sy,sx)/(8*atan(1.))
                            if ( lwrite ) print *,sx/n,sy/n,s
                            if ( s.lt.0 ) s = s + itype
                            if ( iyeartype == -1 .and. s < itype/2 ) s = s + itype
                            res(jx,jy,jz,m,1) = s
                        else if ( imoment.eq.-1 ) then
!
!                           sort results and find perc-th percentile
!                           using the routine already written for the
!                           contingency table
!
                            call nrsort(n,ddata)
                            call getcut(res(jx,jy,jz,m,1),perc,n,ddata)
                            nvars = 1
                        elseif ( abs(imoment).le.100 ) then
!
!                           call Numerical Recipes routine
!
                            call moment(ddata,n,xmom(1),xmom(5),xmom(2), &
     &                          var,xmom(3),xmom(4))
                            do i=1,5
                                res(jx,jy,jz,m,i) = xmom(i)
                            end do
                            if ( xmom(1).ne.0 ) then
                                res(jx,jy,jz,m,6) = xmom(2)/xmom(1)
                            else
                                res(jx,jy,jz,m,6) = 3e33
                            end if
                            res(jx,jy,jz,m,7) = minval(ddata(1:n))
                            res(jx,jy,jz,m,8) = maxval(ddata(1:n))
                            nvars = 8
                            if ( year.ne.0 ) then
                                nvars = 10
                                if ( xyear.lt.1e33 ) then
                                    if ( .not.lchangesign ) then
                                        call fitgau(ddata,n, &
     &                                       xmom(1),xmom(2),a,b, &
     &                                       minindx,maxindx,3,j1,j2, &
     &                                       year,xyear,t,t25,t975, &
     &                                       tx,tx25,tx975,confidenceinterval, &
     &                                       .false.,.false.,lweb,.false.,lwrite)
                                    else
                                        do i=1,n
                                            ddata(i) = -ddata(i)
                                        end do
                                        if ( xyear < 1e33 ) then
                                            xyear = -xyear
                                        end if
                                        call fitgau(ddata,n, &
     &                                       -xmom(1),xmom(2),a,b, &
     &                                       minindx,maxindx,3,j1,j2, &
     &                                       year,xyear,t,t25,t975, &
     &                                       tx,tx25,tx975,confidenceinterval, &
     &                                       .false.,.false.,lweb,.true.,lwrite)
                                    end if
                                    res(jx,jy,jz,m,9) = tx
                                    res(jx,jy,jz,m,10) = &
     &                                   (xyear-xmom(1))/xmom(2)
                                else
                                    res(jx,jy,jz,m,9) = 3e33
                                    res(jx,jy,jz,m,10) = 3e33
                                end if
                            end if
                        elseif ( imoment.eq.200 ) then
!
!                           GPD fit requested
!
                            if ( lchangesign ) then
                                do i=1,n
                                    ddata(i) = -ddata(i)
                                end do
                                if ( xyear < 1e33 ) then
                                    xyear = -xyear
                                end if
                            end if
                            if ( pmindata.le.0 .or. pmindata.ge.100 ) &
     &                           then
                                write(0,*) 'getmomentsfield: error: '// &
     &                               'threshold invalid: ',pmindata
                                call exit(-1)
                            end if
                            call moment(ddata,n,xmom(1),xmom(5),xmom(2), &
     &                          var,xmom(3),xmom(4))
                            call fitgpd(ddata,n,xmom(1),xmom(2),b,xi, &
     &                           j1,j2,lweb,2,lchangesign,pmindata, &
     &                           mindata,year,xyear,t,t25,t975, &
     &                           tx,tx25,tx975,restrain,assume,confidenceinterval, &
     &                           .false.,.false.,lwrite)
                            res(jx,jy,jz,m,1) = b
                            res(jx,jy,jz,m,2) = xi
                            do i=1,10
                                res(jx,jy,jz,m,i+2) = t(i)
                            end do
                            nvars = 12
                            if ( year.ne.0 ) then
                                nvars = 13
                                if ( xyear.lt.3e33 ) then
                                    res(jx,jy,jz,m,13) = tx
                                else
                                    res(jx,jy,jz,m,13) = 3e33
                                end if
                            end if
                        elseif ( imoment.eq.300 ) then
!
!                           GEV fit requested
!
                            if ( lchangesign ) then
                                do i=1,n
                                    ddata(i) = -ddata(i)
                                end do
                                if ( xyear < 1e33 ) then
                                    xyear = -xyear
                                end if
                            end if
                            call moment(ddata,n,xmom(1),xmom(5),xmom(2), &
     &                          var,xmom(3),xmom(4))
                            call fitgev(ddata,n,xmom(1),xmom(2),a,b,xi, &
     &                           j1,j2,lweb,4,lchangesign,year,xyear,t &
     &                           ,t25,t975,tx,tx25,tx975,restrain &
     &                           ,confidenceinterval,.false.,.false.,lwrite)
                            res(jx,jy,jz,m,1) = a
                            res(jx,jy,jz,m,2) = b
                            res(jx,jy,jz,m,3) = xi
                            do i=1,10
                                res(jx,jy,jz,m,i+3) = t(i)
                            end do
                            nvars = 13
                            if ( year.ne.0 ) then
                                nvars = 14
                                if ( xyear.lt.3e33 ) then
                                    res(jx,jy,jz,m,14) = tx
                                else
                                    res(jx,jy,jz,m,14) = 3e33
                                end if
                            end if
                        elseif ( imoment.eq.1000 ) then
!
!                           rank requested
!
                            if ( xyear.gt.1e33 ) then
                                res(jx,jy,jz,m,1) = 3e33
                            else
                                if ( lchangesign ) then
                                    do i=1,n
                                        ddata(i) = -ddata(i)
                                    end do
                                    if ( xyear < 1e33 ) then
                                        xyear = -xyear
                                    end if
                                end if
                                call nrsort(n,ddata)
                                do i=1,n
                                    if ( ddata(n-i+1).lt.xyear ) exit
                                end do
                                res(jx,jy,jz,m,1) = real(i)
                            end if
                            nvars = 1
                        elseif ( imoment.eq.1001 ) then
                            call nrsort(n,ddata)
                            if ( lchangesign ) then
                                xtreme = ddata(1)
                            else
                                xtreme = ddata(n)
                            end if
                            years = -9999
                            iyear = 0
                            do iens=nens1,nens2
                                do yr=yrstart-1,yrstop
                                    do jj=j1,j2
                                        j = jj
                                        call normon(j,yr,i,nperyear)
                                        if ( i.ge.yrstart .and. i.le.yrstop ) then
                                            if ( fxy(j,i,iens).eq.xtreme ) then
                                                if ( iyear.lt.nyears ) iyear = iyear + 1
                                                years(iyear) = i
                                            end if
                                        end if
                                    end do
                                end do
                            end do
                            if ( iyear.eq.0 ) then
                                res(jx,jy,jz,m,1) = 3e33
                            else if ( iyear.eq.1 ) then
                                res(jx,jy,jz,m,1) = years(iyear)
                            else
                                ! take random choice of equal years
                                call random_number(s)
                                i = 1 + int(iyear*s)
                                res(jx,jy,jz,m,1) = years(i)
                            end if
                            nvars = 1
                        else
                            write(0,*) 'error: unknown imoment ',imoment
                            call exit(-1)
                        end if
                        if ( lwrite ) then
                            do i=1,nvars
                                print '(a,3i5,2i3,a,g23.6)', &
     &                               'res(',jx,jy,jz,m,i,') = ', &
     &                               res(jx,jy,jz,m,i)
                            end do
                        end if    
  790                   continue    ! valid point/month
                    end do           ! month
  800               continue        ! valid point
                end do               ! nx
            end do                   ! ny
        end do                       ! nz
!
!       convert to standard units
!
        if ( lstandardunits ) then
            saveunits = units(1)
            do month=m1,m2
                m = month-m1
                do i=1,nvars
                    if ( .not. ( abs(imoment).le.100 .and. &
     &                   (i.eq.3 .or. i.eq.4 .or. i.eq.6 .or. i.eq.9 &
     &                   .or. i.eq.10 ) ) &
     &                   .and. .not. ( imoment.eq.200 .and. i.eq.2 ) &
     &                  .and. .not. imoment.eq.1000 .and. .not. imoment.eq.1001 ) then
                        units(i) = saveunits
                        call makestandardfield(res(1,1,1,m,i),nx,ny,nz,1,0,0 &
     &                      ,nx,ny,nz,1,0,0,vars(1),units(i),lwrite)
                    else
                        units(i) = "1"
                    end if
                end do
            end do
        else
            if ( nvars.eq.8 .or. nvars.eq.9 ) then
                units(3:9) = "1"
                units(2) = units(1)
                units(5) = units(1)
                units(7) = units(1)
                units(8) = units(1)
            elseif ( nvars.eq.12 ) then
                units(2) = "1"
                units(3:12) = units(1)
            end if
        end if
!
!       write output field in GrADS or netcdf format
!
        call savestartstop(yrstart,yrstop)
        call getenv('DIR',dir)
        ldir = len_trim(dir)
        if ( ldir.eq.0 ) ldir=1
        if ( dir(ldir:ldir).ne.'/' ) then
            ldir = ldir + 1
            dir(ldir:ldir) = '/'
        end if
        svars(2:) = ' '
        cell_methods(2:) = cell_methods(1)
        if ( imoment.eq.-1 ) then
            vars(1) = 'perc'
            write(lvars(1),'(f5.2,a)') perc,'% percentile'
        elseif ( imoment.eq.1000 ) then
            vars(1) = 'rank'
            write(lvars(1),'(a,i4,a)') 'rank of year ', &
     &               year,' in the context of the other years'
            units(1) = '1'
        elseif ( imoment.eq.1001 ) then
            vars(1) = 'timex'
            lvars(1) = 'year of extreme'
            units(1) = '1'
        elseif ( abs(imoment).le.100 ) then
            vars(1) = 'mean'
            lvars(1) = 'mean'
            vars(2) = 'sd'
            lvars(2) = 'standard deviation'
            vars(3) = 'skew'
            lvars(3) = 'skewness'
            units(3) = '1'
            vars(4) = 'kurt'
            lvars(4) = 'kurtosis'
            units(4) = '1'
            vars(5) = 'adev'
            lvars(5) = 'absolute deviation'
            vars(6) = 'sdm'
            lvars(6) = 'standard deviation / mean'
            units(6) = '1'
            vars(7) = 'min'
            lvars(7) = 'minimum'
            vars(8) = 'max'
            lvars(8) = 'maximum'
            if ( year.ne.0 ) then
                call check4(year)
                write(vars(9),'(a,i4.4)') 'norm_rt_',year
                write(lvars(9),'(a,i4,a)') 'return time of year ', &
     &               year,' in the context of the other years'
                units(9) = 'yr'
                write(vars(10),'(a,i4.4)') 'z_',year
                write(lvars(10),'(a,i4,a)') 'z-value of year ', &
     &               year,' in the context of the other years'
                units(10) = '1'
            end if
        elseif ( imoment.eq.200 ) then
            svars(1:2) = ' '
            vars(1) = 'pot_scale'
            lvars(1) = 'scale parameter sigma of GPD fit'
            vars(2) = 'pot_shape'
            lvars(2) = 'shape parameter xi of GPD fit'
            units(2) = '1'
            do i=1,4
                write(format,'(a,i1,a)') '(a,i',i+1,',a)'
                write(vars(3*i),format) 't',10**i
                if ( j1.eq.j2 ) then
                    write(lvars(3*i),format) 'return value at ', &
     &                   10**i,' years'
                else
                    write(lvars(3*i),'(a)') &
     &                   'not quite sure what this number means'
                end if
                if ( i.eq.4 ) cycle
                write(vars(3*i+1),format) 't',2*10**i
                if ( j1.eq.j2 ) then
                    write(lvars(3*i+1),format) 'return value at ', &
     &                   2*10**i,' years'
                else
                    write(lvars(3*i+1),'(a)') &
     &                   'not quite sure what this number means'
                end if
                write(vars(3*i+2),format) 't',5*10**i
                if ( j1.eq.j2 ) then
                    write(lvars(3*i+2),format) 'return value at ', &
     &                   5*10**i,' years'
                else
                    write(lvars(3*i+2),'(a)') &
     &                   'not quite sure what this number means'
                end if
            end do
            if ( year.ne.0 ) then
                call check4(year)
                write(vars(13),'(a,i4.4)') 'pot_rt_',year
                write(lvars(13),'(a,i4,a)') 'return time of year ', &
     &               year,' in the context of the other years'
                units(13) = 'yr'
            end if
        elseif ( imoment.eq.300 ) then
            svars(1:3) = ' '
            vars(1) = 'gev_pos'
            lvars(1) = 'position parameter mu of GEV fit'
            vars(2) = 'gev_scale'
            lvars(2) = 'scale parameter sigma of GEV fit'
            vars(3) = 'gev_shape'
            lvars(3) = 'shape parameter xi of GEV fit'
            units(3) = '1'
            do i=1,4
                write(format,'(a,i1,a)') '(a,i',i+1,',a)'
                write(vars(3*i+1),format) 't',10**i
                if ( j1.eq.j2 ) then
                    write(lvars(3*i+1),format) 'return value at ', &
     &                   10**i,' years'
                else
                    write(lvars(3*i+1),'(a)') &
     &                   'not quite sure what this number means'
                end if
                if ( i.eq.4 ) cycle
                write(vars(3*i+2),format) 't',2*10**i
                if ( j1.eq.j2 ) then
                    write(lvars(3*i+2),format) 'return value at ', &
     &                   2*10**i,' years'
                else
                    write(lvars(3*i+2),'(a)') &
     &                   'not quite sure what this number means'
                end if
                write(vars(3*i+3),format) 't',5*10**i
                if ( j1.eq.j2 ) then
                    write(lvars(3*i+3),format) 'return value at ', &
     &                   5*10**i,' years'
                else
                    write(lvars(3*i+3),'(a)') &
     &                   'not quite sure what this number means'
                end if
            end do
            if ( year.ne.0 ) then
                call check4(year)
                write(vars(14),'(a,i4.4)') 'gev_rt_',year
                write(lvars(14),'(a,i4,a)') 'return time of year ', &
     &               year,' in the context of the other years'
                units(14) = 'yr'
            end if
        else
            write(0,*) 'getmomentsfield: error: imoment = ',imoment
            call exit(-1)
        end if
        do i=1,nvars
            ivars(1,i) = nz
            ivars(2,i) = 99
        end do
!       give correlations dates in 0-1
        if ( m1.eq.0 ) then
            i = 0
            j = nperyear
        else
            i = 1
            j = m1
        end if
        if ( index(outfile,'.ctl').ne.0 ) then
            if ( lwrite ) print '(a)','# writing ctl file'
            call args2title(title)
            call writectl(outfile,datfile,nx,xx,ny,yy,nz,zz &
     &           ,1+(m2-m1),nperyear,i,j,absent,title,nvars,vars,ivars &
     &           ,lvars,units)
            print '(a)','# writing output'
            irec = 0
            do month=m1,m2
                m = month-m1
                do i=1,nvars
                    irec = irec + 1
                    write(2,rec=irec) (((res(jx,jy,jz,m,i),jx=1,nx),jy=1 &
     &                   ,ny),jz=1,nz)
                end do
            end do
            close(2)
        else
            if ( lwrite ) print '(a)','# writing netcdf metadata'
            title = 'statistical properties of '//title
            call enswritenc(outfile,ncid,ntvarid,itimeaxis,ntmax,nx,xx &
     &           ,ny,yy,nz,zz,lz,1+(m2-m1),nperyear,i,j,ltime,absent &
     &           ,title,history,nvars,vars,ivars,lvars,svars,units &
     &           ,cell_methods,metadata,0,0)
            do month=m1,m2
                m = month-m1
                do i=1,nvars
                    if ( lwrite ) print '(a,2i3)','# writing netcdf data',m,i
                    call writencslice(ncid,ntvarid,itimeaxis,ntmax, &
     &                   ivars(1,i),res(1,1,1,m,i),nx,ny,nz,nx,ny,nz, &
     &                   m+1,1)
                end do
            end do
            i = nf_close(ncid)  ! do not forget to close file!!!!
        end if
!
!       error messages
!
        goto 999
  901   print *,'getmomentsfield: error reading moment [1-4] from ',trim(line)
        call exit(-1)
  902   print *,'getmomentsfield: error reading percentile from ',trim(line)
        call exit(-1)
  903   print *,'error reading date from file ',trim(line),' at record ',k
        call exit(-1)
  904   print *,'error cannot locate field file file ',trim(line)
        call exit(-1)
  905   print *,'getmomentsfield: error reading year from ',trim(line)
        call exit(-1)
  920   print *,'error cannot open new correlations file ',trim(datfile)
        call exit(-1)
  999   continue
 end program

 subroutine check4(year)
        implicit none
        integer year
        if ( year.gt.9999 .or. year.lt.-999 ) then
            write(0,*) 'getmomentsfield: error: year = ',year
            write(*,*) 'getmomentsfield: error: year = ',year
            call exit(-1)
        end if
end subroutine
