program eof

!       Get the EOFs of a time-varying field, the PCs are computed seperately.
!       Closely follows the Bausteine routine "eof" by Koos Verbeek, but uses
!       Lapack rather than Eispack for greater efficiency.  Does not use the
!       time dimension yet to compensate for this efficiency gain :-(

    use lsdata
    implicit none
    include 'params.h'
    include 'netcdf.inc'
    include 'getopts.inc'
    integer,parameter :: recfa4=4
    integer,parameter :: nvarmax=1,nyrmax=79,neofmax=500,ntmax=24*366
    integer,parameter :: nxymax=180*100,lwork=648108032,liwork=90003 ! ssyevd
    integer :: i,ii,i1,i2,j,jj,j1,j2,jx,jy,k,m,n,n1,n2,nx,ny,nz,nt &
        ,month,yr,mo,nperyear,neigen,ldir,nxf,nyf
    integer :: firstyr,firstmo,nvars,ivars(2,neofmax),jvars(6,neofmax) &
        ,ncid,endian,status,iarg,mens1,mens,fyr,lyr,nnmax,ntvarid,itimeaxis(ntmax)
    integer :: x1,x2,y1,y2,nn(nxmax,nymax),nxy,indx(2,nxymax),isupeof(2*neofmax)
    integer :: info,nrec,mm1,mm2,mm,ncovs,iens,retval &
        ,yrstart,yrstop,lun
    integer,allocatable :: iwork(:)
    real :: xx(nxmax),yy(nymax),zz(nzmax),undef,xxls(nxmax),yyls(nymax)
    real :: fac,sum,wgt,wx(nxmax),wy(nymax),dum,eigen(nxymax),eps,a,b,t,s
    real,allocatable :: cov(:,:),work(:)
    real,allocatable :: field(:,:,:,:,:)
    real,allocatable :: fxy(:,:,:),pc(:,:,:,:),eofxy(:,:,:),pattern(:,:,:,:)
    real*4 :: tarray(2)
    character title*255,vars(neofmax)*60,lvars(neofmax)*200, &
        units(neofmax)*20,lsmasktype*4,outfile*255,infile*255, &
        datfile*255,line*255,yesno*1,dir*255,string*20, &
        fieldtitle*255,lvar*200,unit*80,cell_methods(neofmax)*100, &
        history*50000,ltime*120,lz(3)*20,svars(neofmax)*80, &
        metadata(2,100)*2000,FORM_field*250
    logical :: lexist,xrev,yrev,xwrap,compute_work_sizes
    integer :: getpid,putenv,system
    real*4 :: etime

    lwrite = .false. 
    if ( command_argument_count() < 2 ) then
        print *,'usage: eof field.[ctl|nc] [n] '// &
            '[lsmask maskfile all|lnd|sea] '// &
            '[normalize maxspace|varspace|vartime] '// &
            '[month m[:n] [lag n[:m]] '// &
            '[sum|ave|max|min|sel n] '// &
            '[log|sqrt|rank] '// &
            '[minfac r] [minnum n] [begin yr] [end yr] '// &
            '[diff [nyr]] [detrend] outfield.[ctl|nc]'
        call exit(-1)
    end if
    call get_command_argument(1,infile)
    call getmetadata(infile,mens1,mens,ncid,datfile,nxmax,nx &
        ,xx,nymax,ny,yy,nzmax,nz,zz,lz,nt,nperyear,firstyr,firstmo &
        ,ltime,undef,endian,title,history,nvarmax,nvars,vars,jvars &
        ,lvars,svars,units,cell_methods,metadata,lwrite)
    if ( nz > 1 ) then
        write(0,*) 'eof: error: 3D EOFs not yet ready'
        call exit(-1)
    end if
    if ( index(infile,'%') > 0 .or. &
    index(infile,'++') > 0 ) then
        ensemble = .true. 
    else
        ensemble = .false. 
    end if
    call get_command_argument(2,line)
    if (  ichar(line(1:1)) >= ichar('0') .and. &
    ichar(line(1:1)) <= ichar('9') ) then
        iarg = 3
    else
        iarg = 2
    end if
    call getlsmask(iarg,lsmasktype,nxmax,xxls,nymax,yyls,lwrite)
    if ( lsmasktype /= 'all' ) then
        call checkgridequal(nx,ny,xx,yy,nxls,nyls,xxls,yyls)
    end if
    call getopts(iarg,command_argument_count()-1,nperyear,yrbeg,yrend,.true.,mens1,mens)
    !!!print *,'m1,m2 = ',m1,m2
    call get_command_argument(1,infile)
    lyr = firstyr + (nt+firstmo-1)/nperyear
    lyr = min(yr2,lyr)
    fyr = max(yr1,firstyr)
    if ( lwrite ) print *,'allocating field(',nx,ny,nperyear,fyr,lyr,nens2,')'
    allocate(field(nx,ny,nperyear,fyr:lyr,0:nens2))

    call get_command_argument(2,line)
    if (  ichar(line(1:1)) >= ichar('0') .and. &
    ichar(line(1:1)) <= ichar('9') ) then
        read(line,*,err=901) neigen
        if ( neigen > neofmax ) then
            write(0,*) 'eof: cannot show more than ',neofmax,' EOFs'
            call exit(-1)
        end if
        iarg = 3
    else
        neigen = 4
        iarg = 2
    end if
    print '(a,i3,a,f8.2)','Computing ',neigen, ' EOFs, time: ',etime(tarray)
    if ( lag1 /= 0 .or. lag2 /= 0 ) print *,'eof: lags do not make sense'
    if ( dump ) print *,'eof: dump not supported'
    if ( plot ) print *,'eof: plot not supported'
    if ( lks ) print *,'eof: K-S not supported'
    if ( lconting ) print *,'eof: contingency tables not supported'
    do i=1,indxuse
        if ( lincl(i) ) print *,'eof: what do you mean with ' &
        ,strindx(i),'?'
    end do
!   range of field
    nxf = nx
    nyf = ny
!   range of years
    yr1 = max(yr1,firstyr,firstyr - (min(lag1,lag2)+11)/nperyear)
    yr2 = min(yr2,firstyr + (firstmo+nt-1)/nperyear, &
        firstyr + (firstmo+nt-1)/nperyear - (max(lag1,lag2)-11)/nperyear)
    yrstart = yr2
    yrstop  = yr1
    allocate(pc(nperyear,yr1:yr2,0:nens2,neigen))
    pc = 3e33
    allocate(fxy(nperyear,yr1:yr2,0:nens2))
    allocate(eofxy(nxf,nyf,neigen))
    allocate(pattern(nxf,nyf,nperyear,0:1))

!   read kill info file (in the climexp) and add own PID

    call killfile(string,line,outfile,0)
    call readfield(ncid,infile,datfile,field,nx,ny,nz &
        ,nperyear,fyr,lyr,nens1,nens2,nx,ny,nz,nperyear,yr1,yr2 &
        ,firstyr,firstmo,nt,undef,endian,vars,units,lstandardunits &
        ,lwrite)

!   apply land/sea mask

    call applylsmask(field,lsmask,nx,ny,nz,nperyear,fyr,lyr, &
        nens1,nens2,lsmasktype,lwrite)

!   open output file

    call get_command_argument(command_argument_count(),outfile)
    inquire(file=outfile,exist=lexist)
    if ( lexist ) then
        print *,'output file ',outfile(1:index(outfile,' ')-1), &
            ' already exists, overwrite? [y/n]'
        read(*,'(a)') yesno
        if (  yesno /= 'y' .and. yesno /= 'Y' .and. &
              yesno /= 'j' .and. yesno /= 'J' ) then
            stop
        end if
        open(2,file=outfile)
        close(2,status='delete')
    end if

!   save variables

    fieldtitle = title
    lvar = lvars(1)
    unit = units(1)

!   compute minfac if it has not been set explicitly

    if ( minfac < 0 .and. minnum < 0 ) then
    !           heuristic, gives 0.25 for 150 yrs, 0.5 for 50 yrs, 0.75 for 20yrs
        minfac = max(0.1, &
                 min(0.6, &
                1.5-log(1+real(min(nt,nperyear*(yr2-yr1+1))-1)/nperyear)/4))
    end if
    write(0,'(a,i2,a)') 'Requiring at least ',nint(100*minfac),'% valid points<br>'

!	get boundaries in grid points

    call getxyprop(xx,nx,yy,ny,xrev,yrev,xwrap)
    call getlatlonwindow(lat1,lat2,lon1,lon2,xx,nx,xwrap,avex,yy,ny &
        ,avey,x1,x2,y1,y2,lwrite)

!       average, cut out window - everything to make the arrays smaller

    write(0,'(a,f8.2,a)') 'Averaging, shifting, cutting, time: ',etime(tarray),'<br>'
    call enscutoutwindow(x1,x2,y1,y2,xx,nx,xwrap,xrev,avex,yy,ny &
        ,avey,wx,wy,field,nxf,nyf,nens1,nens2,nperyear,fyr,lyr,yr1,yr2,lwrite)

!   loop over gridpoints carrying out all requested operations

    print '(a,f8.2)','Processing timeseries options, time: ',etime(tarray)
    nnmax = 0
    do jy=1,ny
        do jx=1,nx
            do iens=nens1,nens2
            
!               create 1-D series from field
            
                do i=yr1,yr2
                    do j=1,nperyear
                        fxy(j,i,iens) = field(jx,jy,j,i,iens)
                    end do
                end do
                if ( lwrite .and. jx == (nx+1)/2 .and. &
                     jy == (ny+1)/2 .and. iens == nens1 ) then
                    print *,'to begin with: '
                    call printdatfile(6,fxy(1,yr1,iens),nperyear,nperyear,yr1,yr2)
                end if
            
!               sum
            
                if ( lsum > 1 ) then
                    call sumit(fxy(1,yr1,iens),nperyear,nperyear,yr1,yr2,lsum,oper)
                end if
                if ( lwrite .and. jx == (nx+1)/2 .and. &
                     jy == (ny+1)/2 .and. iens == nens1 ) then
                    print *,'after sumit: '
                    call printdatfile(6,fxy(1,yr1,iens),nperyear,nperyear,yr1,yr2)
                end if
            
!               enforce datacuts
            
                do i=yr1,yr2
                    do j=1,nperyear
                        if ( fxy(j,i,iens) > maxindx .or. fxy(j,i,iens) < minindx ) then
                            fxy(j,i,iens) = 3e33
                        end if
                    end do
                end do

!               sqrt,log

                if ( logscale ) then
                    call takelog(fxy(1,yrbeg,iens),nperyear,nperyear,yr1,yr2)
                end if
                if ( sqrtscale ) then
                    call takesqrt(fxy(1,yrbeg,iens),nperyear,nperyear,yr1,yr2)
                end if
            
!               detrend
            
                if ( ldetrend ) then
                    if ( lwrite ) print *,'Detrending field'
                    call detrend(fxy(1,yr1,iens),nperyear,nperyear,yr1,yr2,yr1,yr2,m1,m2,lsel)
                end if
            
!               differentiate
            
                if ( ndiff /= 0 ) then
                    if ( lwrite ) print *,'Taking differences'
                    call ndiffit(fxy(1,yr1,iens),nperyear,nperyear,yr1,yr2,ndiff,minfacsum)
                end if
                if ( lwrite .and. jx == (nx+1)/2 .and. &
                     jy == (ny+1)/2 .and. iens == nens1 ) then
                    print *,'after diffit: '
                    call printdatfile(6,fxy(1,yr1,iens),nperyear,nperyear,yr1,yr2)
                end if
            end do
        
!           normalize to s.d.
        
            if ( lnormsd ) then ! implies anomal
                call ensnormsd(fxy,nperyear,nperyear,yr1,yr2,nens1,nens2,yr1,yr2)
            else if ( anom .or. lsel == 1 .and. m1 /= 0  ) then
                call ensanomal(fxy,nperyear,nperyear,yr1,yr2,nens1,nens2,yr1,yr2)
            else if ( m1 == m2 ) then
                call getj1j2(mm1,mm2,m1,nperyear, .false. )
                call enssubtractmean(fxy,nperyear,nperyear,yr1,yr2,nens1,nens2,yr1,yr2,mm1,mm2)
            else
                write(0,*) 'eof: can only compute EOFs including ', &
                    ' the seasonal cycle for a single starting month'
                call exit(-1)
            end if
        
!           anomalies wrt ensemble mean
        
            if ( .not. composite .and. nens2 > nens1 .and. lensanom ) then
                call anomalensemble(fxy,nperyear,nperyear,yr1,yr2,yr1,yr2,nens1,nens2)
            end if
        
!           check max number of entries in any single month
        
            do j=1,nperyear
                n = 0
                do iens=nens1,nens2
                    do i=yr1,yr2
                        if ( fxy(j,i,iens) <= 1e33 ) then
                            n = n + 1
                        end if
                    end do
                end do
                nnmax = max(nnmax,n)
            end do
        
!           copy back
        
            nn(jx,jy) = 0
            do iens=nens1,nens2
                do i=yr1,yr2
                    do j=1,nperyear
                        field(jx,jy,j,i,iens) = fxy(j,i,iens)
                        if ( abs(field(jx,jy,j,i,iens)) > 1e20 .and. &
                             field(jx,jy,j,i,iens) < 1e33 ) then
                            print *,'eof: internal error: ',jx,jy,j,i,field(jx,jy,j,i,iens)
                            call exit(-1)
                        end if
                        if ( field(jx,jy,j,i,iens) < 1e33 ) then
                            nn(jx,jy) = nn(jx,jy) + 1
                            yrstart = min(yrstart,i)
                            yrstop  = max(yrstop,i)
                        end if
                    end do
                end do
            end do           ! iens
        
!           flag as invalid if too few points are present
        
            if ( lwrite) print *,'nn(',jx,jy,') = ', &
                nn(jx,jy),xx(jx),wx(jx),yy(jy),wy(jy)
            if ( nn(jx,jy) < max(4,minnum) ) then
                if ( lwrite .and. nn(jx,jy) > 0 ) &
                    print *,'set to zero ',nn(jx,jy),minnum
                nn(jx,jy) = 0
            end if
        end do               ! jx
    end do                   ! jy

!   check whether the problem will fit, set index array n -> i,j

    nxy = 0
    do jy=1,ny
        do jx=1,nx
            if ( nn(jx,jy) > 0 ) then
                if ( nn(jx,jy) < minfac*nnmax ) then
                    if ( lwrite .and. nn(jx,jy) > 0 ) &
                        print *,'set to zero ',nn(jx,jy),minfac,nnmax
                    nn(jx,jy) = 0
                else
                    nxy = nxy + 1
                    if ( nxy <= nxymax ) then
                        indx(1,nxy) = jx
                        indx(2,nxy) = jy
                    end if
                end if
            end if
        end do
    end do
    if ( nxy > nxymax ) goto 901
    write(0,'(a,i5,a,i5,a,f8.2,a)') 'Found approx ',nxy,'/',nx*ny, &
        ' valid points, time: ',etime(tarray),'<br>'
    if ( nxy < neigen ) then
        neigen = nxy
        write(0,*) 'eof: warning: reduced number of eigenvalues to:',nxy
    end if

!	check whether the job is not too large

!	some data points, all month 0 (together)
!       petruk: 350MHz PII, pgf77 -O, ssyevr
!	               nt   nx*ny    nxy    cov  eigen
!	Jones SLP    1476     540    472    22s     1.3s
!	Jones&Parker 1692    2592   1835   365s   102s
!         begin 1925  852    2592   1916   204s   117s
!         ave 2x2    1692     648    474    25s     1.7s
!	  1925, 2x1   852    1296   1010    58s    17s
!       Parker SLP   1488    2664   2664 839/818s   too long
!               2x1  1488    1260   1260   188s    25.6s
!               2x2  1488     648    648    50s     3.4s
!       Broeikas, 194MHz R10000, f77 -O, ssyevr
!       Parker SLP   1488    2664   2664  4216s
!               2x1  1488    1260   1260  1040s    42.6s
!               2x2  1488     648    648   233s     4.2s
!       Zuidzee, 733 MHz PIV, ssyevd
!       Parker SLP   1488    2664   2664   401s >1500s
!               2x1  1488    1260   1260  89.5s    49.8s
!               2x2  1488     648    648  23.4s     6.9s
! pdate 2008
!       bima, 2.16 GHz Intel Core 2 Duo ssyevd, gfortran+vecLib
!       HadSLP2r     1895    2664   2664   183s    13.2s
!               2x1  1895    1260   1260  43.8s     1.4s
!               2x2  1895     648    648  11.1s     0.3s
!       bhlclim, 3.4 GHz Intel Xeon (2x) pgf90
!       HadSLP2r     1895    2664   2664   320s     127s
!               2x1  1895    1260   1260  97.5s    14.5s
!               2x2  1895     648    648  19.9s     2.3s
!       zuidzee, 2.7 GHz Intel Xeon pgf90
!       HadSLP2r     1895    2664   2664   968s     165s
!               2x1  1895    1260   1260   265s    39.7s
!               2x2  1895     648    648  68.2s     2.7s
!       superstorm 1.8 GHZ AMD Opteron (2x) g95
!       HadSLP2r     1895    2664   2664   383s    66.5s
!               2x1  1895    1260   1260  88.8s    17.9s
!               2x2  1895     648    648  22.5s     1.9s
!       equator 3.1 GHz Intel Xeon pgf90
!       HadSLP2r     1895    2664   2664   794s     142s
!               2x1  1895    1260   1260   178s    14.9s
!               2x2  1895     648    648  45.6s     2.2s

!	so the following formula seems to hold
!		cov = A*nt*nxy**2 + B*nxy**3
!	with for the 350MHz PII A =  80ns, B = 17ns
!	and for the broeikas    A = 440ns, B = 21ns ssyevr
!       and for the zuidzee     A =  40ns, B = 25ns ssyevd
!       for the new zuidzee     A =  70ns, B = 9ns
!       for the bhlclim         A =  24ns, B = 7ns
!       for the bima            A =  14ns, B = 0.7ns
    if ( m1 == 0 ) then
        n = nt
    else
        n = (yr2-yr1+1)*(m2+lsel-m1)
    end if
!**	a = 80e-9
!**	b = 17e-9
    a = 4*40e-9
    b = 25e-9
    t = a*n*(nens2-nens1+1)*real(nxy)**2 + b*real(nxy)**3
    t = t/3 ! for the new server
!!!        write(0,'(a,2i6)') '# n,nxy = ',n,nxy
    write(0,'(a,f8.0,a)') 'Estimated time to completion: ',t ,'s<p>'

!   open output file

    do i=1,neigen
        if ( i < 10 ) then
            write(vars(i),'(a,i1)') 'eof',i
        elseif ( i < 100 ) then
            write(vars(i),'(a,i2)') 'eof',i
        elseif ( i < 1000 ) then
            write(vars(i),'(a,i3)') 'eof',i
        else
            write(vars(i),'(a,i4)') 'eof',i
        end if
        ivars(1,i) = 1
        ivars(2,i) = 99
        write(lvars(i),'(a,i4,2a)') 'EOF #',i,' of'
        lvars(i) = trim(lvars(i))//' '//trim(lvar)
        svars = ' '
        if ( normalization > 0 ) then
            units(i) = '1'
        else
            units(i) = unit
        end if
        if ( i > 1 ) cell_methods(i) = cell_methods(1)
    end do
    call getenv(FORM_field,FORM_field)
!   adjust x axis to what it was before we started playing
    if ( lon1 /= 0 .or. lon2 /= 360 ) then
        if ( abs(xx(1)-360-lon1) < abs(xx(1)-lon1) ) then
            do i=1,nx
                xx(i) = xx(i) - 360
            end do
        elseif ( abs(xx(1)+360-lon1) < abs(xx(1)-lon1) ) then
            do i=1,nx
                xx(i) = xx(i) + 360
            end do
        end if
    end if
!   give EOF dates in 2000-2001
    if ( m1 == 0 ) then
        i = 2000
        j = nperyear
    else
        i = 2001
        j = m1
    end if
    k = index(outfile,'.ctl')
    if ( k /= 0 ) then
        datfile = outfile(:k-1)//'.grd'
        open(2,file=datfile,form='unformatted',access='direct' &
            ,recl=recfa4*nx*ny,err=920)
        nrec = 0
        ncid = 0
        call args2title(title)
        call writectl(outfile,datfile,nx,xx,ny,yy,1,zz,1+(m2-m1) &
            ,nperyear,i,j,3e33,title,neigen,vars,ivars,lvars,units)
    else
        undef = 3e33
        call enswritenc(outfile,ncid,ntvarid,itimeaxis,ntmax,nx,xx,ny &
            ,yy,nz,zz,lz,1+(m2-m1),nperyear,i,j,ltime,undef,title &
            ,history,neigen,vars,ivars,lvars,svars,units,cell_methods &
            ,metadata,0,0)
    end if

!   compute covariance matrix

    allocate(cov(nxymax,nxymax))
    ncovs = 0
    do month=m1,m2
!       memory corruption or bug in gfortran - but it happens...
        if ( month > m2 ) exit
        call getj1j2(mm1,mm2,month,nperyear, .false. )
        write(0,'(a,i3,a,i3,a,f8.2,a)') &
            'Computing covariance for periods ',mm1,'-', &
            1+mod(mm2-1,nperyear),', time: ',etime(tarray),'<p>'
        do n2=1,nxy
            call keepalive1('Computing covariance matrix entry',n2,nxy)
            i2 = indx(1,n2)
            j2 = indx(2,n2)
            if ( lwrite ) print '(a,i5,a,i5,2i4)', &
            'row ',n2,'/',nxy,i2,j2
            do n1=n2,nxy
                i1 = indx(1,n1)
                j1 = indx(2,n1)
                sum = 0
                n = 0
                do yr=yr1-1,yr2
                    do mm=mm1,mm2
                        m = mm
                        i = yr
                        if ( m > nperyear ) then
                            m = m-nperyear
                            i = i+1
                        end if
                        if ( i < yr1 .or. i > yr2 ) cycle
                        do iens=nens1,nens2
                            if ( field(i1,j1,m,i,iens) < 1e33 .and. &
                                 field(i2,j2,m,i,iens) < 1e33 ) then
                                sum = sum + field(i1,j1,m,i,iens)* &
                                field(i2,j2,m,i,iens)
                                n = n + 1
                            end if
                        end do
                    end do
                end do
                if ( n > 3 &
                     .and. n > minfac*nnmax*(mm2-mm1) &
                     .and. n > minnum ) then
                    if ( wx(i1)*wx(i2)*wy(j1)*wy(j2) < 0 ) then
                        write(0,*) 'eof: error: negative weights ' &
                            ,i1,wx(i1),i2,wx(i2),j1,wy(j1),j2,wy(j2)
                        call exit(-1)
                    end if
                    cov(n1,n2) = sqrt(wx(i1)*wx(i2)*wy(j1)*wy(j2))*sum/n
                    ncovs = ncovs + 1
                else
                    if ( lwrite ) print &
                        '(a,4i4,a,2i5,a,i5,f6.0,f5.2,i4,2i3,2i5)', &
                        'not enough points at ',i1,j1,i2,j2,'(' &
                        ,n1,n2,'): ',n,minfac*nnmax*(mm2-mm1) &
                        /nperyear,minfac,minnum,mm1,mm2,yr1,yr2
                    cov(n1,n2) = 0 ! best guess...
                end if
                if ( n1 /= n2 ) then
                    cov(n2,n1) = cov(n1,n2)
                end if
            end do
        end do
        if ( lwrite ) then
            print *,'Covariance matrix:'
            do j=1,nxy
                print '(20f12.2)',(cov(i,j),i=1,min(nxy,20))
            end do
        end if
        if ( ncovs <= 2 ) then
            print *,'no points with valid covariances'
            write(0,*)'no points with valid covariances'
            call exit(-1)
        end if
    
!       compute eigenvectors - LAPACK routine, see manpage.
    
        dum = 0
    eps = 1e-3
        write(0,'(a,i8,a,f8.2,a)')'Computing eigenvalues, nxy= ' &
            ,nxy,', time: ',etime(tarray),'<p>'
        i = lwork
        j = liwork

!       estimate time needed for routine and start keepalive job

        call getenv('DIR',dir)
        ldir = len_trim(dir)
        if ( ldir <= 1 ) then
            dir = '/usr/people/oldenbor/climexp/'
            ldir= len_trim(dir)
        end if
        if ( dir(ldir:ldir) /= '/' ) then
            ldir = ldir + 1
            dir(ldir:ldir) = '/'
        end if
        if ( lweb ) then
            s = b*real(nxy)**3/30 ! half minutes
            write(string,'(i19)') getpid()
            write(line,'(2a,i3,3a)') dir(1:len_trim(dir)) &
                ,'./stillcomputing.cgi ',nint(s),' ', &
                trim(string),'&'
            call mysystem(line,retval)
            if ( retval /= 0 ) then
                write(0,*) 'eof: error: ',trim(line),' failed: ',retval
            end if
        end if
    
!       and call lapack
    
        if ( lwrite) print '(a)','Calling LAPACK routine ssyevd'
        compute_work_sizes = .false. ! run once if nxymax changes
        if ( compute_work_sizes ) then
            i = -1
            j = -1
            nxy = nxymax
        end if
        allocate(iwork(liwork))
        allocate(work(lwork))
        call ssyevd('V','U',nxy,cov,nxymax,eigen,work,i,iwork,j,info)
        if ( i == -1 .or. j == -1 ) then
            print *,'for nxymax =  ',nxy
            print *,'set lwork to  ',work(1)
            print *,'set liwork to ',iwork(1)
            call exit(-1)
        end if
        if ( info /= 0 ) then
            if ( info < 0 ) then
                write(0,*) 'lapack ssyevr: argument ',-info,' illegal'
            else
                write(0,*) 'lapack ssyevr: internal error ',info
            end if
            call exit(-1)
        end if
    
!       treat eigenvalues, note that they are returned in *as*cending order
    
        sum = 0
        do i=1,nxy
            sum = sum + max(0.,eigen(i))
            if ( lwrite .and. i > nxy-30 ) print *,i,eigen(i),sum
        end do
        write(0,'(a,f8.2)')'Eigenvalues, time: ',etime(tarray)
        write(0,'(a)')'<table class="realtable" '// &
            'border=0 cellpadding=0 cellspacing=0>'// &
            '<tr><th>#</th><th>eigenvalue</th>'
        write(0,'(a)')'<th>explained variance</th>'// &
            '<th>cumulative</th></tr>'
        s = 0
        do i=nxy,nxy-neigen+1,-1
            s = s + eigen(i)
            write(0,'(a,i3,a,g16.5,a,f8.2,a,f8.2,a)') &
                '<tr><td>',nxy-i+1,'</td><td>',eigen(i) &
                ,'</td><td align=right>',eigen(i)/sum*100 &
                ,'%</td><td align=right>',s/sum*100 &
                ,'%</td></tr>'
        end do
        write(0,'(a)')'</table>'
    
        print '(a,f8.2)','Eigenvectors (EOFs), time: ',etime(tarray)
        do i=1,neigen
!           set to absent (we will miss holes)
            do j1=1,ny
                do i1=1,nx
                    eofxy(i1,j1,i) = 3e33
                end do
            end do
!           translate back to latlon grid
            do n1=1,nxy
                i1 = indx(1,n1)
                j1 = indx(2,n1)
                if ( cov(n1,nxy-i+1) /= 0 ) then
                    if ( lwrite ) print '(a,2i6,a,i4,a,i4,a,3g12.3)' &
                        ,'cov(',n1,nxy-i+1,'),wx(',i1,'),wy(',j1 &
                        ,') = ',cov(n1,nxy-i+1),wx(i1),wy(j1)
                    eofxy(i1,j1,i) = &
                    cov(n1,nxy-i+1)/sqrt(wx(i1)*wy(j1))
                    if ( lwrite ) print '(a,3i5,a,g12.6)','eofxy(' &
                        ,i1,j1,i,') = ',eofxy(i1,j1,i)
                end if
            end do
        end do
        if ( normalization > 0 ) then
            call normeof1(eofxy,nxf,nyf,nx,ny,neigen,normalization)
        end if
        do i=1,neigen
!           construct a pattern file suitable for use in project
            pattern = 3e33
            if ( month /= 0 ) then
                pattern(1:nx,1:ny,month,1) = eofxy(1:nx,1:ny,i)
            else
                pattern(1:nx,1:ny,nperyear,0) = eofxy(1:nx,1:ny,i)
            end if
            call project(pc(1,yr1,0,i),nperyear, &
                yr1,yr2,nens1,nens2,xx,nx,yy,ny, &
                field,pattern,nxf,nyf,nperyear,fyr,lyr, &
                month,minfac,lnomissing,(m1 /= 0),anom,lwrite)
        end do
        if ( normalization < 0 ) then
            call normeof2(pc,eofxy,nperyear,yr1,yr2,mm1,mm2, &
                nens1,nens2,neigen,nxf,nyf,nx,ny, &
                normalization,lnomissing,lwrite)
        end if
    
!       write out pattern
    
        if ( ncid == 0 ) then
            do i=1,neigen
                nrec = nrec + 1
                write(2,rec=nrec) ((eofxy(i1,j1,i),i1=1,nx),j1=1,ny)
            end do
        else
            nrec = nrec + 1
            do i=1,neigen
                call writencslice(ncid,ntvarid,itimeaxis,ntmax,ivars(1,i), &
                    eofxy(1,1,i),nxf,nyf,nz,nx,ny,nz,nrec,1)
            end do
        end if
    end do
    if ( ncid /= 0 ) status = nf_close(ncid)

!   write out time series

    if ( ncid == 0 ) then
        datfile(index(datfile,'.grd'):) = ' '
    else
        datfile = outfile(:index(outfile,'.nc')-1)
    end if
    k = len_trim(datfile)
    do iens=nens1,nens2
        do i=1,neigen
            if ( ensemble ) then
                write(datfile(k+1:),'(a,i2.2,a,i2.2,a)') '_',i,'_',iens,'.dat'
            else
                write(datfile(k+1:),'(a,i2.2,a)') '_',i,'.dat'
            end if
            call rsunit(lun)
            open(lun,file=datfile)
            if ( fieldtitle /= ' ' ) then
                fieldtitle = infile
            end if
            call printmetadata(lun,infile,FORM_field,fieldtitle,history,metadata)
            if ( normalization > 0 ) then
                if ( i < 10 ) then
                    write(lun,'(a,i1,3a,i1,2a)') '# PC',i,' [',trim(unit), &
                    '] PC ',i,' of ',trim(fieldtitle)
                else
                    write(lun,'(a,i2,3a,i2,2a)') '# PC',i,' [',trim(unit), &
                    '] PC ',i,' of ',trim(fieldtitle)
                end if
            else
                if ( i < 10 ) then
                    write(lun,'(a,i1,3a)') '# PC',i,' [1]'
                else
                    write(lun,'(a,i2,3a)') '# PC',i,' [1]'
                end if
            end if
!           spread value over whole season
            if ( m1 == m2 .and. m1 > 0 .and. lsum > 1 ) then
                do yr=yr1,yr2
                    do m=2,lsum
                        mo = m1 + m - 1
                        call normon(mo,yr,ii,nperyear)
                        if ( ii <= yr2 ) then
                            pc(mo,ii,iens,i) = pc(m1,yr,iens,i)
                        end if
                    end do
                end do
            end if
            call printdatfile(lun,pc(1,yr1,iens,i),nperyear,nperyear,yr1,yr2)
            close(lun)
        end do
    end do

!   finito

    call savestartstop(yrstart,yrstop)
    write(string,'(i19)') getpid()
    write(line,'(4a)') dir(1:len_trim(dir)),'./kill_stillcomputing.cgi ',trim(string)
    call mysystem(line,retval)
    if ( retval /= 0 ) then
        write(0,*) 'eof: error: ',trim(line),' failed: ',retval
    end if
    print '(a,f8.2)','finished, time: ',etime(tarray)

!	error messages

    return
901	write(0,*) 'eof: error: covariance array not big enough'
    write(0,*) '     available size: ',nxymax,', need ',nxy
    write(0,*) '     try selecting a smaller region or averaging'
    write(*,*) 'eof: error: covariance array not big enough'
    write(*,*) '     available size: ',nxymax,', need ',nxy
    write(*,*) '     try selecting a smaller region or averaging'
call exit(-1)
920 print *,'error cannot open new EOF file ',datfile(1:index(datfile,' ')-1),recfa4,nx,ny
    call exit(-1)
end program eof
