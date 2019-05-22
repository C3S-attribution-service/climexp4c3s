program correlatefield
!
!   program to correlate a field series to a point series 
!   to give fields of correlation coefficients, probabilities that
!   these are significant, and the fit coefficients a, b and their
!   errors.
!   aug-2002 added support for ensembles; ensemble field x scalar
!   time series, and ensemble field x ensemble time series the
!   third option is not (yet?) implemented
!
    use lsdata
    implicit none
    include 'params.h'
    include 'netcdf.inc'
    include 'getopts.inc'
    integer,parameter :: nvarmax=50,nyrmax=50,nlevmax=1,mensmax=1,mpermax=366
    integer,parameter :: ntmax=1000,ndatmax=2*mpermax*(yrend-yrbeg+1)*10
    integer,parameter :: nmc=1000
    integer,parameter :: recfa4=4
    real,parameter :: absent=3e33
    integer :: i,j,n,nx,ny,nz,nt,firstyr,firstmo,nvars,iarg,             &
         ivars(2,nvarmax),ncid,jvars(6,nvarmax),endian,status,        &
         nperyear,lastyr,mens,mens1
    integer :: jj,k,kk,nn,lag,jx,jy,jz,yr,month,j1,j2,m,mm,mo,ii &
         ,l,ldir,ntvarid,itimeaxis(ntmax)    &
         ,nrec,iens,jens,ndup(0:mpermax),validens(nensmax)     &
         ,nens2series,iens2,imens(0:1),nold,yrstart,yrstop   &
         ,fyr,mdatmax,irec,ntp,ndiffn,nmetadata
    integer,allocatable :: yrmo(:,:)
    real,allocatable :: field(:,:,:,:,:,:),r(:,:,:,:),prob(:,:,:,:), &
         a(:,:,:,:),b(:,:,:,:),da(:,:,:,:),db(:,:,:,:),             &
         a1(:,:,:,:),da1(:,:,:,:),cov(:,:,:,:),relregr(:,:,:,:),    &
         drelregr(:,:,:,:),                                         &
         rmin(:,:,:,:),rmax(:,:,:,:),zdif(:,:,:,:),rprob(:,:,:,:),  &
         xn(:,:,:,:)
    real,allocatable :: data(:,:,:),mcdata(:,:,:),fxy(:,:,:)
    real,allocatable :: aaa1(:,:,:,:,:),bbb1(:,:,:,:,:),            &
         aaa(:),bbb(:),field2(:,:,:,:,:,:)
    real :: xx(nxmax),yy(nymax),zz(nzmax),undef,xxls(nxmax),yyls(nymax)
    real,allocatable :: ddata(:),dindx(:),dddata(:)
    real :: adata,sxx,aindx,syy,sxy,df,d,zd,z,probd,sig(1),chi2, &
         q,sum,fac,aa,daa,bb,dbb,dresult(-2:2),         &
         results(nmc),rmins(nmc),rmaxs(nmc),zdifs(nmc),dum,zold,    &
         sxxold,alpha,xrand,s
    logical :: lexist,ensseries,llwrite
    logical,allocatable :: lnewyr(:,:),lfirst(:)
    character :: title*512,vars(nvarmax)*60,lvars(nvarmax)*128,        &
          units(nvarmax)*60,lsmasktype*4
    character :: invars(nvarmax)*60,inlvars(nvarmax)*128,intitle*255,inunits(nvarmax)*60
    character :: line*80,yesno*1,string*10,file*255,infile*255,        &
         datfile*255,outfile*255,dir*255,ensfile*255,var*60,unit*60, &
         lvar*120,svar*120,tmpunits*60,tmpvars*60,string1*42,string2*42
    character :: lz(3)*20,svars(100)*100,ltime*120,history*50000,serieshistory*50000, &
        cell_methods(100)*100,metadata(2,100)*2000,seriesmetadata(2,100)*2000, &
        seriestitle*100
    integer :: rindex
!
!   check arguments
!
    n = command_argument_count()
    if ( lwrite ) print *,'correlatefield: called with ',n,' arguments'
    if ( n < 3 ) then
        print *,'usage: correlatefield field.[ctl|nc] '//           &
             'series [month m[:n] [[ens]anom] [lag n[:m]] '//       &
             '[sum|ave|max|min|sel n] [sum2 n] '//                  &
             '[log|sqrt|rank] '//                                   &
             '[minfac r] [minnum n] [begin yr] [end yr] '//         &
             '[lt cut] [gt cut] [diff [nyr]] [detrend] '//          &
             '[runcorr|runregr nyr dummyfile '//                    &
             '[random series|field]] [noise white|red]'//           &
              'outfield'
        print *,'       ensemble input is denoted by %% in the name'
        call exit(-1)
    end if
    call get_command_argument(1,infile)
    call getmetadata(infile,mens1,mens,ncid,datfile,nxmax,nx &
        ,xx,nymax,ny,yy,nzmax,nz,zz,lz,nt,nperyear,firstyr,firstmo &
        ,ltime,undef,endian,title,history,nvarmax,nvars,vars,ivars &
        ,lvars,svars,units,cell_methods,metadata,lwrite)
    iarg = 2
    call getlastyr(firstyr,firstmo,nt,nperyear,lastyr)
    call add_varnames_metadata(vars(1),lvars(1),svars(1),metadata,'field_variable')
!
!   process arguments
!
    j=3
    call getlsmask(j,lsmasktype,nxmax,xxls,nymax,yyls,lwrite)
    if ( lsmasktype /= 'all' ) then
        call checkgridequal(nx,ny,xx,yy,nxls,nyls,xxls,yyls)
    end if
    if ( index(infile,'%') > 0 .or. index(infile,'++') > 0 ) then
        ensemble = .true.
    else
        ensemble = .false.
    end if
    call getopts(j,n-1,nperyear,yrbeg,yrend,.true.,mens1,mens)
    if ( ensemble ) write(0,*) 'Using ensemble members ',nens1      &
          ,' to ',nens2,'<br>'
    call get_command_argument(2,infile)
!
!   allocate arrays
!
    allocate(field(nx,ny,nz,nperyear,firstyr:lastyr,nens1:nens2))
    allocate(lsmask(nx,ny))
    allocate(r(nx,ny,nz,0:nperyear))
    allocate(prob(nx,ny,nz,0:nperyear))
    allocate(a(nx,ny,nz,0:nperyear))
    allocate(b(nx,ny,nz,0:nperyear))
    allocate(da(nx,ny,nz,0:nperyear))
    allocate(db(nx,ny,nz,0:nperyear))
    allocate(cov(nx,ny,nz,0:nperyear))
    allocate(a1(nx,ny,nz,0:nperyear))
    allocate(da1(nx,ny,nz,0:nperyear))
    allocate(relregr(nx,ny,nz,0:nperyear))
    allocate(drelregr(nx,ny,nz,0:nperyear))
    allocate(rmin(nx,ny,nz,0:nperyear))
    allocate(rmax(nx,ny,nz,0:nperyear))
    allocate(zdif(nx,ny,nz,0:nperyear))
    allocate(rprob(nx,ny,nz,0:nperyear))
    allocate(xn(nx,ny,nz,0:nperyear))
!
    intitle = title
    invars = vars
    intitle = title
    inunits = units
    n = command_argument_count()
    if ( lag1 < 0 ) print *,'(point leading field)'
    if ( lag2 > 0 ) print *,'(field leading point)'
    if ( dump ) write(0,*)'correlatefield: dump not supported'
    if ( plot ) write(0,*)'correlatefield: plot not supported'
    if ( lks ) write(0,*)'correlatefield: K-S not supported'
    if ( lconting ) write(0,*)'correlatefield: contingency tables not supported'
    do i=1,indxuse
        if ( lincl(i) ) write(0,*)'correlatefield: what do you mean with ',strindx(i),'?'
    end do
    if ( composite ) then
        if ( maxdata > 1e33 .and. mindata < -1e33 .and.           &
              (pmaxdata <= 0 .or. pmaxdata >= 100) .and.            &
              (pmindata <= 0 .or. pmindata >= 100) ) then
            write(0,*) 'correlatefield: cannot make composites '//  &
                 'without cut-offs ',maxdata,mindata,pmaxdata       &
                 ,pmindata
            write(*,*) 'correlatefield: cannot make composites '//  &
                  'without cut-offs',maxdata,mindata,pmaxdata       &
                 ,pmindata
            call exit(-1)
        end if
    end if
    call get_command_argument(n,outfile)
    inquire(file=outfile,exist=lexist)
    if ( lexist ) then
        print *,'output file ',trim(outfile),' already exists, overwrite? [y/n]'
        read(*,'(a)') yesno
        if (  yesno /= 'y' .and. yesno /= 'Y' .and.                 &
              yesno /= 'j' .and. yesno /= 'J' ) then
            call exit(-1)
        end if
        open(1,file=outfile)
        close(1,status='delete')
    end if
!   save time on the initialization - but not too much.
    yr1 = max(yr1,firstyr,firstyr - (min(lag1,lag2)+nperyear-1)     &
          /nperyear)
    yr2 = min(yr2,firstyr + (firstmo+nt-2)/nperyear,                &
          firstyr + (firstmo+nt-2)/nperyear - (max(lag1,lag2)       &
          -nperyear+1)/nperyear)
    if ( nyrwindow > 0 ) then
        write(0,'(2a,i6,a)') 'p-values are computed'               &
             ,' against a ',nmc,' sample Monte Carlo<br>'
    end if
!
!   read kill info file (in the climexp) and add own PID
!
    call killfile(dir,title,infile,0)
!
!   allocate time series arrays
!
    allocate(data(mpermax,yrbeg:yrend,0:nensmax))
    allocate(mcdata(mpermax,yrbeg:yrend,0:nensmax))
    allocate(fxy(mpermax,yrbeg:yrend,0:nensmax))
    allocate(lnewyr(mpermax,yrbeg:yrend))
    if ( ncrossvalidate > 0 ) then
        mdatmax = nperyear*(lastyr-firstyr+1)*(nens2-nens1+1)
        allocate(aaa(mdatmax))
        allocate(bbb(mdatmax))
        allocate(aaa1(nx,ny,nz,nperyear,yr1:yr2))
        allocate(bbb1(nx,ny,nz,nperyear,yr1:yr2))
        aaa1 = 3e33
        bbb1 = 3e33
    end if
    allocate(yrmo(3,ndatmax),lfirst(ndatmax))
    allocate(ddata(ndatmax),dindx(ndatmax),dddata(ndatmax))
!
!   init
!
    print *,'init'
    nold = -1
    zold = 3e33
    sxxold = 3e33
    if ( lwrite ) print *,'filling fields with absent'
    if ( m1 /= m2 .and. lag1 /= lag2 ) then
        print *,'Sorry, can only handle either lags varying or'// &
              ' months varying, not both'
        print *,'(months:',m1,m2,', lags:',lag1,lag2,')'
        call exit(-1)
    end if
    if ( lag2-lag1 > nperyear ) then
        print *,'Sorry, can only store ',nperyear+1,            &
             ' fields maximum'
        call exit(-1)
    end if
    do iens=nens1,nens2
        if ( lwrite ) print *,'calling makeabsent data',iens
        call makeabsent(data(1,yrbeg,iens),mpermax,yrbeg,yrend)
        if ( lwrite ) print *,'calling makeabsent fxy ',iens
        call makeabsent(fxy(1,yrbeg,iens),mpermax,yrbeg,yrend)
    end do
!
!   read field, change absent values to our convention
!   if subtracting fields, read the whole field; otherwise just the part needed
!
    if ( lsubtract ) then
        yr1a = firstyr
        yr2a = lastyr
    else
        yr1a = yr1
        yr2a = yr2
    end if
    do iens=nens1,nens2
        call keepalive1('Reading ensemble member',iens-nens1+1,nens2-nens1+1)
        if ( ncid == -1 ) then
            call get_command_argument(1,infile)
            if ( ensemble ) then
                call filloutens(infile,iens)
            end if
            if ( lwrite ) print *,'calling parsectl on ',trim(infile)
            call parsectl(infile,datfile,nxmax,nx,xx,nymax,ny,yy    &
                 ,nzmax,nz,zz,nt,nperyear,fyr,firstmo,undef         &
                 ,endian,title,1,nvars,vars,ivars,lvars,units)
            if (lwrite) print '(2a)','# looking for ',trim(datfile)
            inquire(file=datfile,exist=lexist)
            if ( .not.lexist ) then
                if (lwrite) print '(3a)','# looking for ',trim(datfile),'.gz'
                inquire(file=trim(datfile)//'.gz',exist=lexist)
                if ( .not.lexist ) then
                    nens2 = iens-1
                    if ( nens2 >= nens1 ) then
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
            call zreaddatfile(datfile,field(1,1,1,1,firstyr,iens),  &
                  nx,ny,nz,nx,ny,nz,nperyear,firstyr,lastyr,        &
                 fyr,firstmo,nt,undef,endian,lwrite,yr1a,yr2a,      &
                 1,1)
        else
            call get_command_argument(1,infile)
            if ( ensemble ) then
                call filloutens(infile,iens)
            end if
            if ( lwrite ) print *,'calling parsenc on ',trim(infile)
            status = nf_open(infile,nf_nowrite,ncid)
            call parsenc(infile,ncid,nxmax,nx,xx,nymax,ny,yy        &
                 ,nzmax,nz,zz,nt,nperyear,fyr,firstmo,              &
                 undef,title,1,nvars,vars,jvars,lvars,units)
            call zreadncfile(ncid,field(1,1,1,1,firstyr,iens),      &
                 nx,ny,nz,nx,ny,nz,nperyear,firstyr,lastyr,fyr,     &
                 firstmo,nt,undef,lwrite,yr1a,yr2a,jvars)
        end if
    end do
5   continue
    if ( lstandardunits ) then
!       convert to standard units
        do iens=nens1,nens2
            tmpunits = units(1) ! they are otherwise adjusted
            call makestandardfield(field(1,1,1,1,firstyr,iens),nx,ny &
                 ,nz,nperyear,firstyr,lastyr,nx,ny,nz,nperyear      &
                 ,max(firstyr,yr1a),min(lastyr,yr2a),vars(1)        &
                 ,tmpunits,lwrite)
        end do
        units(1) = tmpunits
        if ( lwrite ) then
            print *,'correlatefield: just after standard units'
            print *,'field(',(nx+1)/2,(ny+1)/2,(nz+1)/2,firstmo     &
                 ,max(firstyr,yr1),') = ',field((nx+1)/2,(ny+1)/2   &
                 ,(nz+1)/2,firstmo,max(firstyr,yr1),nens1)
        end if
    end if
!
!   apply land/sea mask
!
    call applylsmask(field,lsmask,nx,ny,nz,nperyear,firstyr,lastyr, &
         nens1,nens2,lsmasktype,lwrite)
!
!   read series
!
    call get_command_argument(iarg,file)
    if ( index(file,'%%') == 0 .and. index(file,'++') == 0 ) then
        ensseries = .FALSE.
        nens2series = nens1
        print *,'reading file ',file(1:index(file,' ')-1)
        call readseriesmeta(file,data(1,yrbeg,nens1),mpermax,yrbeg,yrend &
             ,n,var,unit,lvar,svar,serieshistory,seriesmetadata,lstandardunits,lwrite)
    else
        ensseries = .TRUE.
        if ( ensemble ) then
            nens2series = nens2
        else
            write(0,*) 'Cannot correlate a single field with an '// &
                  'ensemble series yet'
            call exit(-1)
        end if
        do iens=nens1,nens2series
            ensfile=file
            call filloutens(ensfile,iens)
            if ( lwrite ) write(0,*) 'looking for file ',trim(ensfile)
            inquire(file=ensfile,exist=lexist)
            if ( .not.lexist ) goto 10
            print *,'reading file ',ensfile(1:index(ensfile,' ')-1)
            if ( iens == nens1 ) then
                call readseriesmeta(ensfile,data(1,yrbeg,iens),mpermax,yrbeg &
                    ,yrend,n,var,unit,lvar,svar,serieshistory,seriesmetadata, &
                    lstandardunits,lwrite)
            else
                call readseries(ensfile,data(1,yrbeg,iens),mpermax,yrbeg &
                    ,yrend,n,var,unit,lstandardunits,lwrite)
            end if
        end do
        goto 11
10      continue
        if ( iens <= nens1 ) then
            write(0,*) 'error: cannot find time series'
            call exit(-1)
        end if
        write(0,*) 'Using series ensemble members ',nens1,' to '    &
              ,iens-1,'<br>'
        nens2series = iens - 1
        if ( nens2 < 0 ) then
            write(0,*)                                              &
                  'correlatefield: error: could not find ensemble ' &
                  ,trim(file),trim(ensfile)
            call exit(-1)
        end if
11      continue
    end if
    call add_varnames_metadata(var,lvar,svar,seriesmetadata,'variable')

    if ( n /= nperyear ) then
        write(0,*) 'correlatefield: error: cannot interpolate '//   &
              'in time (yet)',nperyear,n
        write(*,*) 'correlatefield: error: cannot interpolate '//   &
              'in time (yet)',nperyear,n
        call exit(-1)
    end if
    do iens=nens1,nens2series
!
!       take monthly anomalies
!
        if ( mdiff > 0 ) then
            if ( lwrite ) print *,'taking monthly anomalies'
            call mdiffit(data(1,yrbeg,iens),mpermax,nperyear,yrbeg  &
                  ,yrend,mdiff)
        end if
!
!       sum series
!
        if ( lsum > 1 ) then
            if ( lwrite ) print *,'taking sum'
            call sumit(data(1,yrbeg,iens),mpermax,nperyear,yrbeg    &
                  ,yrend,lsum,oper)
        end if
!
!       log, sqrt
!
        if ( logscale ) then
            call takelog(data(1,yrbeg,iens),mpermax,nperyear,yrbeg  &
                  ,yrend)
        end if
        if ( sqrtscale ) then
            call takesqrt(data(1,yrbeg,iens),mpermax,nperyear,yrbeg &
                  ,yrend)
        end if
!
!       detrend data
!
        if ( ldetrend ) then
            i = len_trim(file)
            if ( index(file,'time') == 0 ) then
                if ( lwrite ) print *,'detrending series'
                call detrend(data(1,yrbeg,iens),mpermax,nperyear,   &
                     yrbeg,yrend,yr1,yr2,m1,m2,lsel)
            end if
        end if
!
!       differentiate data
!
        if ( ndiff /= 0 ) then
            if ( lwrite ) print *,'Taking differences/averaging'
            call ndiffit(data(1,yrbeg,iens),mpermax,nperyear,yrbeg  &
                  ,yrend,ndiff,minfacsum)
        end if
        if ( ndiff2 /= 0 ) then
            if ( lwrite ) print *,'Taking differences/averaging2'
            call ndiffit(data(1,yrbeg,iens),mpermax,nperyear,yrbeg  &
                  ,yrend,ndiff2,minfacsum)
        end if
    end do
!
!   composites
!
    if ( pmindata > 0 .and. pmindata < 100 .or.                   &
         pmaxdata > 0 .and. pmaxdata < 100 ) then
        if ( fix2 ) then
            n = 0
        else
            if ( lag1 == lag2 ) then
                n = lag1
            else
                write(0,*) 'correlatefield: error: can only handle', &
                     ' a single lag in composite analysis'
                call exit(-1)
            end if
        end if
        if ( pmindata > 0 .and. pmindata < 100 ) then
            call getj1j2(j1,j2,m1,nperyear,.false.)
            call getenscutoff(mindata,pmindata,data,mpermax          &
                 ,nperyear,yrbeg,yrend,nensmax,nens1,nens2series,yr1 &
                 ,yr2,j1,j2,n)
            write(0,'(a,f6.2,a,g14.6,a)') 'Converted ',             &
                 pmindata,'% to ',mindata,'<br>'
        end if
        if ( pmaxdata > 0 .and. pmaxdata < 100 ) then
            call getj1j2(j1,j2,m1,nperyear,.false.)
            call getenscutoff(maxdata,pmaxdata,data,mpermax          &
                 ,nperyear,yrbeg,yrend,nensmax,nens1,nens2series,yr1 &
                 ,yr2,j1,j2,n)
            write(0,'(a,f6.2,a,g14.6,a)') 'Converted ',             &
                 pmaxdata,'% to ',maxdata,'<br>'
        end if
    end if
    if ( composite ) then
        call whencomposite(data,mpermax,yrbeg,yrend,nensmax         &
             ,nens2series,lnewyr,nperyear)
    end if
!
!   anomalies - necessary if we consider more than one month
!
    if ( .not.composite .and. anom ) then
        if ( lwrite ) print *,'Taking anomalies'
        do iens=nens1,nens2
            call anomal(data(1,yrbeg,iens),mpermax,nperyear,yrbeg   &
                  ,yrend,yr1,yr2)
        end do
    end if
!
!   ensemble anomalies
!
    if ( .not.composite .and. nens2series > nens1 .and. lensanom ) then
        if ( lwrite ) print *,'Taking anomalies wrt ensemble'
        call anomalensemble(data,mpermax,nperyear,yrbeg,yrend,      &
             yr1,yr2,nens1,nens2series)
    end if
!
!   adjust yr1,yr2 to not be longer than the time series
!
    call adjustyr(yr1,yr2,data(1,yrbeg,nens1),mpermax,nperyear,yrbeg,yrend)
    call adjustyr(yr2,yr1,data(1,yrbeg,nens1),mpermax,nperyear,yrbeg,yrend)
!
!   compute number of time steps with data (<=nt) for minfac cuts
!
    ntp = 0
    do yr=yr1,yr2
        do mo=1,nperyear
            do iens=nens1,nens2
                do jz=1,nz
                    do jy=1,ny
                        do jx=1,nx
                            if ( field(jx,jy,jz,mo,yr,iens) < 1e33 ) then
                                ntp = ntp + 1
                                goto 101
                            end if
                        end do
                    end do
                end do
            end do
101         continue
        end do
    end do
!   check for ntp a multiple of 365 while nperyear = 366
    s = ntp/real(nperyear)
    if ( abs(s/nint(s)-1) < 0.01 ) then
        ntp = nperyear*nint(s)
    end if
    if ( lwrite ) print *,'nt,ntp = ',nt,ntp
!
!   compute minfac if it has not been set explicitly
!
    if ( minfac < 0 .and. minnum < 0 ) then ! not used in the web scripts
!       heuristic, gives 0.25 for 150 yrs, 0.5 for 50 yrs, 0.75 for 20yrs
        minfac = max(0.1,                                           &
             min(0.6,                                               &
             (1.5-log(1+real(min(ntp,nperyear*(yr2-yr1+1))-1)  &
             /nperyear)/4)))
    end if
    if ( minfac > 0 ) then
        write(0,'(a,i2,a)') 'Requiring at least ',nint(100*minfac),'% valid points<br>'
    end if
!   compute percentage valid points in data
    if ( maxdata < 1e33 .or. mindata > -1e33 ) then
        n = 0
        nn = 0
        call getj1j2(j1,j2,m1,nperyear,.false.)
        do yr=yr1-1,yr2+1
            do jj=j1,j2
                if ( fix2 ) then
                    j = jj+lag1
                else
                    j = jj
                end if
                call normon(j,yr,i,nperyear)
                if ( i >= yr1 .and. i <= yr2 ) then
                    if ( data(j,i,nens1) < 1e33 ) then
                        nn = nn + 1
                        if ( ((data(j,i,nens1) <= maxdata) .and.    &
                                (data(j,i,nens1) >= mindata)) .eqv.    &
                                (maxdata >= mindata) ) then
                            n = n + 1
                        end if
                    end if
                end if
            end do
        end do
        !!!write(0,*) '@@@ n,nn = ',n,nn,'<br>'
        ! adjust for valid points in data, includinhg 
        minfac = minfac*real(n)/real(nn) 
        !!!write(0,*) '@@@ minfac = ',minfac
    end if
!
!   loop over grid points
!
    print *,'correlating'
    yrstart = yr2
    yrstop  = yr1
    do jz=1,nz
        do jy=1,ny
            do jx=1,nx
                if ( nyrwindow == 0 ) then
                    if ( mod(jx,10) == 1 ) then
                        call keepalive1('Computing correlations for latitude ', &
                         jy+(jz-1)*ny,ny*nz)
                    end if
                else
                    call keepalive1('Computing correlations for point ', &
                         jx+(jy-1)*nx+(jy-1)*(jz-1)*nx*ny,nx*ny*nz)
                end if
                do month=0,nperyear
                    r(jx,jy,jz,month) = absent
                    prob(jx,jy,jz,month) = absent
                    a(jx,jy,jz,month) = absent
                    b(jx,jy,jz,month) = absent
                    da(jx,jy,jz,month) = absent
                    db(jx,jy,jz,month) = absent
                    relregr(jx,jy,jz,month) = absent
                    drelregr(jx,jy,jz,month) = absent
                    xn(jx,jy,jz,month) = absent
                end do
!
!               create 1-D series from field
!
                n = 0
                do iens=nens1,nens2
                    do i=yr1,yr2
                        do j=1,nperyear
                            fxy(j,i,iens) = field(jx,jy,jz,j,i,iens)
                            if ( fxy(j,i,iens) < 0.9*absent ) then
                                n = n+1
                            end if
                        end do
                    end do
                    if ( n < 3 ) then
                        if ( lwrite ) print '(a,3i5,a,3f7.2,a,i5,a)'    &
                             ,'not enough valid points at ',jx,jy,jz    &
                             ,' (',xx(jx),yy(jy),zz(jz),'): ',n         &
                             ,' < 3'
                        goto 800
                    end if
!
!                   take monthly anomalies
!
                    if ( mdiff2 > 0 ) then
                        call mdiffit(fxy(1,yrbeg,iens),mpermax          &
                              ,nperyear,yrbeg,yrend,mdiff2)
                    end if
!       
!                       sum
!
                    if ( lsum2 > 1 ) then
                        call sumit(fxy(1,yrbeg,iens),mpermax            &
                              ,nperyear,yrbeg,yrend,lsum2,'v')
                    end if
!
!                   log, sqrt
!
                    if ( logfield ) then
                        call takelog(fxy(1,yrbeg,iens),mpermax          &
                             ,nperyear,yrbeg,yrend)
                    end if
                    if ( sqrtfield ) then
                        call takesqrt(fxy(1,yrbeg,iens),mpermax         &
                             ,nperyear,yrbeg,yrend)
                    end if
!
!                   detrend
!
                    if ( ldetrend ) then
                        if ( lwrite ) print *,'Detrending field'
                        if ( lag1 == 0 .and. lag2 == 0 .or. m1 == 0     &
                              .or.lsel == 12 ) then
                            call detrend(fxy(1,yrbeg,iens),mpermax      &
                                 ,nperyear,yrbeg,yrend,yr1,yr2,m1,m2    &
                                 ,lsel)
                        else
                            call detrend(fxy(1,yrbeg,iens),mpermax      &
                                 ,nperyear,yrbeg,yrend,yr1,yr2,1        &
                                 ,nperyear,lsel)
                        end if
                    end if
!       
!                   differentiate
!
                    if ( ndiff /= 0 ) then
                        if ( lwrite ) print *,'Taking differences/averaging'
                        call ndiffit(fxy(1,yrbeg,iens),mpermax          &
                             ,nperyear,yrbeg,yrend,ndiff,minfacsum)
                    end if
                    if ( ndiff2 /= 0 ) then
                        if ( lwrite ) print *,'Taking differences/averaging2'
                        call ndiffit(fxy(1,yrbeg,iens),mpermax          &
                             ,nperyear,yrbeg,yrend,ndiff2,minfacsum)
                    end if
!       
!                   anomalies
!
                    if ( anom .or. lsel > 1 .and. ndiff <= 0 .and.     &
                          ndiff2 <= 0 .or. composite ) then
                        if ( lwrite ) print *,'Taking anomalies'
                        call anomal(fxy(1,yrbeg,iens),mpermax           &
                              ,nperyear,yrbeg,yrend,yr1,yr2)
                    end if
                end do       ! loop over ensemble members
!
!               anomalies wrt ensemble mean
!
                if ( nens2 > nens1 .and. lensanom ) then
                    if ( lwrite ) print *,'Taking anomalies wrt ensemble mean'
                    call anomalensemble(fxy,mpermax,nperyear,           &
                         yrbeg,yrend,yr1,yr2,nens1,nens2)
                end if
!
!               copy ensemble members so that there is the same
!               number of valid ones at every time step
!       
                do j=0,nperyear
                    ndup(j) = 0
                end do
                if ( ensemble ) then
!**                        call makeensfull(ndup,nperyear,fxy,mpermax,yrbeg
!**     +                        ,yrend,nens1,nens2,validens,lwrite)
                    do j=1,nperyear
                        do i=yrbeg,yrend
                            mens = 0
                            do iens=nens1,nens2
                                if ( fxy(j,i,iens) < 1e33 ) then
                                    mens = mens + 1
                                    validens(mens) = iens
                                end if
                            end do
                            if (  mens > 1 .and. mens <= (nens2-nens1) ) then
!                                   some - but not all - are undefined
!                                   copy as many whole sets in as fit
                                if ( lwrite ) print                     &
                                      '(a,i4,i2,a,i4,a,2i4)','At ',i    &
                                      ,j,' only found ',mens            &
                                      ,' ensemble members out of '      &
                                      ,nens1,nens2
                                k = (nens2-nens1+1)/mens
                                iens = nens1
                                do l=1,k-1
                                    do jens=1,mens
300                                     continue
                                        if ( fxy(j,i,iens) < 1e33 ) then
                                            if ( iens >= nens2 ) then
                                                write(0,*) 'error: ',iens,nens2
                                                do k=nens1,nens2
                                                    write(0,*) fxy(j,i,k)
                                                end do
                                                call exit(-1)
                                            end if
                                            iens = iens + 1
                                            goto 300
                                        end if
                                        if ( .FALSE. .and. lwrite ) then
                                              print *,'copying ',validens(jens),' to ' &
                                              ,iens,fxy(j,i,validens(jens)),fxy(j,i,iens)
                                        end if
                                        ndup(j) = ndup(j) + 1
                                        fxy(j,i,iens) = fxy(j,i,validens(jens))
                                        if ( ensseries ) then
                                            data(j,i,iens) = data(j,i,validens(jens))
                                        end if
                                    end do
                                end do
!                               and fill out the rest with random members
                                do l=nens1+k*mens,nens2
!                                   check...
                                    do k=1,mens
                                        if ( validens(k) >= 0 ) goto 305
                                    end do
                                    write(0,*) 'error: no more members'
                                    call exit(-1)
305                                 continue
310                                 continue
                                    call random_number(xrand)
                                    k = 1 + int(mens*xrand)
                                    if ( k <= 0 .or. k > mens .or. validens(k) < 0 ) goto 310
320                                 continue
                                    if ( fxy(j,i,iens) < 1e33 ) then
                                        iens = iens + 1
                                        goto 320
                                    end if
                                    if ( .FALSE. .and. lwrite ) then
                                          print *,'copying member '     &
                                          ,validens(k),' to ',iens      &
                                          ,fxy(j,i,validens(k))         &
                                          ,fxy(j,i,iens)
                                    end if
                                    ndup(j) = ndup(j) + 1
                                    fxy(j,i,iens) = fxy(j,i,validens(k))
                                    if ( ensseries ) then
                                        data(j,i,iens) = data(j,i,validens(k))
                                    end if
                                    validens(k) = -1
                                end do
                            end if
                        end do
                    end do
                end if
!
!               correlate!
!
                do month=m1,m2
                    call getj1j2(j1,j2,month,nperyear,.false.)
                    if ( ensemble ) then
                        ndup(0) = 0
                        do j=j1,j2
                            ndup(0) = ndup(0) + ndup(j)
                        end do
                    end if
                    do lag=lag1,lag2
!       
!                       fill linear arrays without absent values and compute r
!
                        n = 0
                        do yr=yr1-1,yr2
                            do jj=j1,j2
                                if ( fix2 ) then
                                    j = jj+lag
                                else
                                    j = jj
                                end if
                                call normon(j,yr,i,nperyear)
                                if ( i < yr1 .or.i > yr2 ) goto 710
                                m = j-lag
                                call normon(m,i,ii,nperyear)
                                if ( ii < yrbeg .or.ii > yrend ) goto 710
                                if ( yr == yr1 .and. jx == nx/2 .and. jy == ny/2 .and. jz == 1 ) then
                                    print '(a,i4,i5,i4,i5,a)','Correlating months/days', &
                                        j,i,m,ii,' of point and field'
                                end if
                                do iens=nens1,nens2
                                    if ( ensseries ) then
                                        iens2 = iens
                                    else
                                        iens2 = nens1
                                    end if
                                    !!!print *,'@@@  fxy,maxindx,minindx = ',fxy(m,ii,iens),maxindx,minindx
                                    !!!print *,'@@@ data,maxdata,mindata = ',data(j,i,iens2),maxdata,mindata
                                    if ( fxy(m,ii,iens) < 1e33 .and. ((        &
                                          ( fxy(m,ii,iens) <= maxindx ) .and.    &
                                          ( fxy(m,ii,iens) >= minindx )) .eqv.   & 
                                          ( maxindx >= minindx ) )                & 
                                          .and.                                 &
                                          data(j,i,iens2) < 1e33 .and. ((      &
                                          ( data(j,i,iens2) <= maxdata ) .and.    &
                                          ( data(j,i,iens2) >= mindata )) .eqv.   &
                                          ( maxdata >= mindata ) ) ) then
                                        !!!print *,'@@@ OK'
                                        n = n+1
                                        if ( n > ndatmax ) goto 909
                                        ddata(n) = data(j,i,iens2)
                                        dindx(n) = fxy(m,ii,iens)
!!!                                            print *,'@data(',j,i,iens2,
!!!     +                                           ') = ',data(j,i,iens2)
!!!                                            print *,'@ fxy(',m,ii,iens,
!!!     +                                           ') = ',fxy(m,ii,iens)
                                        yrmo(1,n) = i
                                        yrmo(2,n) = j
                                        yrmo(3,n) = iens
                                        yrstart = min(yrstart,i,ii)
                                        yrstop  = max(yrstop,i,ii)
                                    end if ! valid point
                                end do ! iens
710                             continue
                            end do ! jj
                        end do ! yr
                        if ( m1 /= m2 ) then
                            m = month-m1
                        elseif ( .not.fix2 ) then
                            m = lag2-lag
                        else
                            m = lag-lag1
                        end if
                        if ( jx == nx/2 .and. jy == ny/2 .and. jz == 1 ) then
                            print '(a,i2)','and storing it at ',m
                        end if
                        xn(jx,jy,jz,m) = n
                        ndiffn = 1 + max(0,-ndiff)
                        if ( month == 0 .and. n < minfac*              &
                             min(ntp/ndiffn,nperyear*(yr2-yr1+1))       &
                             .or. month /= 0 .and. n < minfac*         &
                             min(ntp/nperyear/ndiffn,yr2-yr1+1)         &
                             .or. n < minnum ) then
                            if ( lwrite ) print '(a,3i5,a,3f7.2,a,2i3,a,2i6)' &
                                 ,'not enough valid points at ',jx      &
                                 ,jy,jz,' (',xx(jx),yy(jy),zz(jz)       &
                                 ,')',month,lag,': ',n,ntp
                            goto 790
                        end if
                        if ( ncrossvalidate > 0 ) then
                            if ( n > mdatmax ) then
                                write(0,*) 'error: mdatmax<n ',mdatmax,n
                                call exit(-1)
                            end if
                        end if
                        call fitcross(dindx,ddata,n,sig,0,              &
                             aa,a(jx,jy,jz,m),daa,da(jx,jy,jz,m)        &
                             ,chi2,q,ncrossvalidate,aaa,bbb,.true.)
                        if ( ncrossvalidate > 0 ) then
                            do i=1,n
                                aaa1(jx,jy,jz,yrmo(2,i),yrmo(1,i)) = aaa(i)
                                bbb1(jx,jy,jz,yrmo(2,i),yrmo(1,i)) = bbb(i)
                            end do
                        end if
                        call fitcross(ddata,dindx,n,sig,0,              &
                             bb,b(jx,jy,jz,m),dbb,db(jx,jy,jz,m)        &
                             ,chi2,q,ncrossvalidate,aaa,bbb,.false.)
                        noisetype = 1
                        if ( ensseries ) then
                            imens(0) = nens2
                        else
                            imens(0) = 0
                        end if
                        imens(1) = nens2
                        if ( .false. ) then
                            do k=1,n
                                print *,k,dindx(k),ddata(k)
                            end do
                        end if
                        if ( index(file,'/time') /= 0 .or. file(1:4) == 'time') then
                            alpha = 0
                        else if ( fitfunc > 1 ) then
                            write(0,*) 'cirrelatefield: warning: ',     &
                                'serial correlations are not taken',    &
                                ' into account for non-linear fit'
                            alpha = 0
                        else
                            call getred(alpha,j1,j2,lag,1,nperyear      &    
                                 ,imens,1,data,fxy,mpermax,yrbeg        &
                                 ,yrend,nensmax,bb,b(jx,jy,jz,m))
                        end if
                        if ( lwrite ) then
                            print *,'j1,j2 = ',j1,j2
                            call getdf(df,month,n,ndup(0),decor,        &
                                 max(lsum,lsum2),max(1,1-ndiff,         &
                                 1-ndiff2),nperyear)
                            print *,'df was',df,n,ndup(0)
                        end if
                        if ( alpha == 1 .or. n-ndup(0) <= 0 ) then
                            df = -2
                        elseif ( alpha > 2/sqrt(real(n-ndup(0)))       &
                                 .and.alpha > exp(-1.) ) then
                            df = n/(1 - 1/(2*log(alpha))) - 2
                        else
                            df = n - 2
                        end if
                        if ( lwrite ) then
                            print *,'df  is',df,alpha
                        end if
                        if ( df < 1 ) then
                            df = 3e33
                        end if
                        if ( df < 1 .or. df > 1e33 ) then
                            r(jx,jy,jz,m) = 3e33
                            prob(jx,jy,jz,m) = 3e33
                            a(jx,jy,jz,m) = 3e33
                            da(jx,jy,jz,m) = 3e33
                            b(jx,jy,jz,m) = 3e33
                            db(jx,jy,jz,m) = 3e33
                            a1(jx,jy,jz,m) = 3e33
                            da1(jx,jy,jz,m) = 3e33
                            if ( lwrite ) print '(a,f5.2,3i5)'          &
                                 ,'error: df < 1: ',df,jx,jy,jz
                            goto 790
                        end if
                        if ( lfitnoise ) then
                            write(0,*) 'correlatefield: error: option fitnoise no longer supported'
                            write(0,*) 'correlatefield: error: option fitnoise no longer supported'
                            call exit(-1)
                        elseif ( lrank ) then
                            if ( lwrite ) then
                                if ( month == 0 ) then
                                    sum = max(lsum,lsum2) + decor
                                else
                                    sum = 1 + (max(lsum,lsum2)-1)/nperyear + decor/nperyear
                                end if
                                print *,'sum was: ',sum
                            end if
                            sum = n/(df+2)
                            if ( lwrite ) then
                                print *,'sum  is: ',sum,n,df
                            end if
                            call spearx(ddata,dindx,n,ddata,dindx,d     &
                                  ,zd,probd,r(jx,jy,jz,m),              &
                                  prob(jx,jy,jz,m),sum,adata,sxx        &
                                  ,aindx,syy)
                        elseif ( composite ) then
                            call makecomposite(ddata,dindx,n,           &
                                 b(jx,jy,jz,m),db(jx,jy,jz,m),          &
                                 prob(jx,jy,jz,m),df,lwrite)
                        elseif ( fitfunc == 1 ) then
                             if ( df <= 0 ) then
                                if ( lwrite ) print '(a,f5.2,3i5)'      &
                                  ,'error: df <= 0: ',df,jx,jy,jz
                                goto 790
                            end if
!                               adjust error estimates from fit()
                            da(jx,jy,jz,m) = da(jx,jy,jz,m)*sqrt((n-2)/df)
                            db(jx,jy,jz,m) = db(jx,jy,jz,m)*sqrt((n-2)/df)
                            call pearsncross(ddata,dindx,n,             &
                                 r(jx,jy,jz,m),prob(jx,jy,jz,m),z,      &
                                 adata,sxx,aindx,syy,sxy,df             &
                                 ,ncrossvalidate)
                            if ( sxx == 0 .or. syy == 0 ) then
                                r(jx,jy,jz,m) = absent
                                prob(jx,jy,jz,m) = absent
                            else
                                cov(jx,jy,jz,m) = sxy/(n-1)
                            end if
                            if ( abs(aindx) > 1e-33 .and. aindx < 1e33 .and. &
                                     b(jx,jy,jz,m) < 1e33 ) then
                                relregr(jx,jy,jz,m) = b(jx,jy,jz,m)/aindx
                                drelregr(jx,jy,jz,m) = db(jx,jy,jz,m)/aindx
                            else
                                relregr(jx,jy,jz,m) = 3e33
                                drelregr(jx,jy,jz,m) = 3e33
                            end if
                            if ( nyrwindow > 0 ) then
                                if ( ensseries ) then
                                    imens(0) = nens2
                                else
                                    imens(0) = 0
                                end if
                                imens(1) = nens2
                                if ( nfittime /= 0 ) goto 901
                                call getruncorr(dindx,ddata,lfirst      &
                                     ,dddata,ndatmax,j1,j2,lag,1        &
                                     ,month,nperyear,imens,1,fxy        &
                                     ,data,mpermax,yrbeg,yrend          &
                                     ,nensmax,ndup(0)                   &
                                     ,' ',.false.,.false.               &
                                     ,rmin(jx,jy,jz,m)                  &
                                     ,rmax(jx,jy,jz,m)                  &
                                     ,zdif(jx,jy,jz,m))
                                if ( zdif(jx,jy,jz,m) > 1e30 ) then
                                    rprob(jx,jy,jz,m) = absent
                                    goto 750
                                end if
!
!                                   Monte Carlo to determine significance
!
                                call filllinarray(dindx,ddata           &
                                     ,lfirst,dddata,ndatmax,n,j1        &
                                     ,j2,lag,1,nperyear,imens,1         &
                                     ,fxy,data,mpermax,yrbeg            &
                                     ,yrend,nensmax,-999,-999,yrmo)
                                if ( lwrite ) print *,'correlate: a,b,sd,r,sde '    &
                                     ,a(jx,jy,jz,m),b(jx,jy,jz,m),sqrt(sxx/(n-1))   &
                                     ,r(jx,jy,jz,m),sqrt(sxx/(n-1)*(1-r(jx,jy,jz,m)**2))
!
!                                   copy the data to a scratch array
!                                   to indicate where the holes are
!
                                if ( ensseries ) then
                                    do iens=nens1,nens2
                                        do i=yr1,yr2
                                            do j=1,12
                                                mcdata(j,i,iens) = data(j,i,iens)
                                            end do
                                        end do
                                    end do
                                else
                                    do i=yr1,yr2
                                        do j=1,12
                                            mcdata(j,i,0) = data(j,i,0)
                                        end do
                                    end do
                                end if
                                do i=1,nmc
                                    call makemcseries(mcdata,fxy        &
                                         ,mpermax,yrbeg,yrend           &
                                         ,nensmax,1,nperyear,1          &
                                         ,lag,j1,j2,imens,adata         &
                                         ,sxx,aindx,syy,sxy             &
                                         ,alpha,n)
                                    if ( lwrite ) then
                                        lwrite = .false.
                                        call filllinarray(dindx,ddata,lfirst        &
                                             ,dddata,ndatmax,n,j1,j2,lag,1          &
                                             ,nperyear,imens,1,fxy,mcdata,mpermax   &
                                             ,yrbeg,yrend,nensmax,-999,-999,yrmo)
                                        call printcorr(dindx,ddata,lfirst           &
                                             ,dddata,yrmo,n,ndup(0)                 &
                                             ,j1,j2,month,nperyear,lag,' ',         &
                                             .false.,.false.,results(i),dresult,dum)
                                        lwrite = .true.
                                    end if
                                    call getruncorr(dindx,ddata,lfirst,dddata,ndatmax   &
                                         ,j1,j2,lag,1,month,nperyear,imens,1,fxy        &
                                         ,mcdata,mpermax,yrbeg,yrend,nensmax,ndup(0)    &
                                         ,' ',.false.,.false.,rmins(i),rmaxs(i),zdifs(i))
                                    if ( lwrite ) print *,'zdif = ',zdifs(i)
                                end do
                                if ( lwrite ) call getsign('result '            &
                                     ,r(jx,jy,jz,m),results,nmc,1               &
                                     ,rprob(jx,jy,jz,m),lwrite)
                                call getsign('zdif',zdif(jx,jy,jz,m)            &
                                     ,zdifs,nmc,1,rprob(jx,jy,jz,m),            &
                                     lwrite.or..true.)
750                             continue
                            end if
                        else ! fitfunc>1, nonlinear fit
                            if ( fitfunc < 1 .or. fitfunc > 3 ) then
                                write(0,*) 'correlatefield: error: '            &
                                     //'cannot handle polynomial '//            &
                                     'of degree ',fitfunc,' yet'
                                write(*,*) 'correlatefield: error: '            &
                                     //'cannot handle polynomial '//            &
                                     'of degree ',fitfunc,' yet'
                                call exit(-1)
                            end if
                            call getdf(df,month,n,ndup(0),decor,max(lsum,lsum2), &
                                 max(1,1-ndiff,1-ndiff2),nperyear)
                            if ( df <= 0 ) then
                                if ( lwrite ) print '(a,f5.2,3i5)'              &
                                  ,'error: df <= 0: ',df,jx,jy,jz
                                goto 790
                            end if
                            if ( fitfunc == 2 ) then
                                call fit2(ddata,dindx,n,a(jx,jy,jz,m),b(jx,jy,jz,m),    &
                                     a1(jx,jy,jz,m),da(jx,jy,jz,m),db(jx,jy,jz,m),      &
                                     da1(jx,jy,jz,m),r(jx,jy,jz,m),prob(jx,jy,jz,m),    &
                                     df,ncrossvalidate,lwrite)
                            else
                                write(0,*) 'correlatefield: error: cubic fitting not yet ready'
                                write(*,*) 'correlatefield: error: cubic fitting not yet ready'
                                call exit(-1)
                            end if
                        end if
                        if ( lwrite .and. fitfunc == 1 ) then
                            print '(a,3i5,2i3,a,i6,a,2f8.4,4f12.4)','point ',   &
                              jx,jy,jz,month,lag,' OK (',n,'): ',               &
                              r(jx,jy,jz,m),prob(jx,jy,jz,m),                   &
                              a(jx,jy,jz,m),da(jx,jy,jz,m),                     &
                              b(jx,jy,jz,m),db(jx,jy,jz,m)
                        end if
                        if ( lwrite .and. fitfunc > 1 ) then
                            print '(a,3i5,2i3,a,i6,a,2f8.4,6f12.4)','point ',   &
                              jx,jy,jz,month,lag,' OK (',n,'): ',               &
                              r(jx,jy,jz,m),prob(jx,jy,jz,m),                   &
                              a(jx,jy,jz,m),da(jx,jy,jz,m),                     &
                              b(jx,jy,jz,m),db(jx,jy,jz,m),                     &
                              a1(jx,jy,jz,m),da1(jx,jy,jz,m)
                        end if
790                     continue ! valid point/month
                    end do   ! lag
                end do       ! month
800             continue     ! valid point
            end do           ! nx
        end do               ! ny
    end do                   ! nz
    if ( index(outfile,'.ctl') /= 0 ) then
        i = index(outfile,'.ctl')
        datfile = outfile(:i-1)//'.grd'
        open(2,file=datfile,form='unformatted',access='direct'  &
              ,recl=recfa4*nx*ny*nz,err=920)
    elseif ( lsubtract ) then
        write(*,*) 'correlatefield: error: netcdf output not yet ready'
        write(0,*) 'correlatefield: error: netcdf output not yet ready'
        call exit(-1)
    end if
    if ( .not.lsubtract ) then
        call getenv('DIR',dir)
        ldir = len_trim(dir)
        if ( ldir == 0 ) ldir=1
        if ( dir(ldir:ldir) /= '/' ) then
            ldir = ldir + 1
            dir(ldir:ldir) = '/'
        end if
        if ( title == ' ' ) call get_command_argument(1,title)
        title = 'Linear correlations and regressions of '//trim(title)//' and '//trim(file)
        seriestitle = ' '
        call merge_metadata(metadata,nmetadata,seriesmetadata,seriestitle,serieshistory,'series_')
        tmpunits = units(1)
        tmpvars = vars(1)
        units = ' '
        svars = ' '
        if ( composite ) then
            nvars = 7
            vars(1) = 'unknown3'
            lvars(1) = 'unknown3'
            vars(2) = 'prob'
            lvars(2) = 'p-value'
            vars(3) = 'unknown'
            lvars(3) = 'nothing'
            vars(4) = 'composite'
            lvars(4) = 'composite'
            vars(5) = 'nothing'
            lvars(5) = 'error on unknown'
            vars(6) = 'errorcomp'
            lvars(6) = 'error on composite'
            vars(7) = 'unknown2'
            lvars(7) = 'unknown2'
        elseif ( lfitnoise ) then
            nvars = 9
            vars(1) = 'skeweps'
            vars(2) = 'dskeweps'
            vars(3) = 'a2'
            vars(4) = 'bn'
            vars(5) = 'da2'
            vars(6) = 'dbn'
            vars(7) = 'sdeps'
            vars(8) = 'a1'
            vars(9) = 'da1'
            lvars(1) = 'skew noise term'
            lvars(2) = 'error on skew noise term'
            lvars(3) = 'nonlinear term a2'
            lvars(4) = 'multiplicative noise term'
            lvars(5) = 'error on nonlinear term a2'
            lvars(6) = 'error on multiplicative noise term'
            lvars(7) = 'standard deviation of noise'
            lvars(8) = 'linear term a2'
            lvars(9) = 'error on linear term a2'
        else
            if ( fitfunc == 1 ) then
                nvars = 9
            else
                nvars = 8
            end if
            vars(1) = 'corr'
            if ( lrank ) then
                lvars(1) = 'rank correlation'
            else
                lvars(1) = 'correlation'
            end if
            units(1) = '1'
            vars(2) = 'prob'
            lvars(2) = 'p-value'
            units(2) = '1'
            vars(3) = 'regr1'
            lvars(3) = 'regression of series '//trim(var)//' on field '//trim(tmpvars)
            call divideunits(units(3),unit,tmpunits)
            vars(4) = 'regr'
            lvars(4) = 'regression of field '//trim(tmpvars)//' on series '//trim(var)
            call divideunits(units(4),tmpunits,unit)
            vars(5) = 'd'//vars(3)
            lvars(5) = 'error on '//lvars(3)
            units(5) = units(3)
            vars(6) = 'd'//vars(4)
            lvars(6) = 'error on '//lvars(4)
            units(6) = units(4)
            if ( fitfunc == 1 ) then
                vars(7) = 'n'
                vars(8) = 'relregr'
                vars(9) = 'drelregr'
                lvars(7) = 'number of valid points'
                lvars(8) = 'relative regression'
                lvars(9) = 'error on relative regression'
                units(7) = '1'
                units(8) = '1'
                units(9) = '1'
            elseif ( fitfunc == 2 ) then
                vars(7) = 'quad'
                lvars(7) = 'coefficient of quadratic term'
                vars(8) = 'errorquad'
                lvars(8) = 'error on coefficient of quadratic term'
            else
                write(0,*) 'error: fitfunc = ',fitfunc,' not yet supported'
                write(*,*) 'error: fitfunc = ',fitfunc,' not yet supported'
                call exit(-1)
            end if
        end if
        if ( nyrwindow > 0 .and. irunvar > 0 ) then
            if ( irunvar == 1 ) then
                vars(nvars+1) = 'rmin'
                vars(nvars+2) = 'rmax'
                vars(nvars+3) = 'zdif'
            elseif ( irunvar == 2 ) then
                vars(nvars+1) = 'bmin'
                vars(nvars+2) = 'bmax'
                vars(nvars+3) = 'bdif'
            else
                write(0,*) 'correlatefield: unknown irunvar ',irunvar
                write(*,*) 'correlatefield: unknown irunvar ',irunvar
                call exit(-1)
            end if
            vars(nvars+4) = 'rprob'
            lvars(nvars+1) = 'minimum running correlation'
            lvars(nvars+2) = 'maximum running correlation'
            lvars(nvars+3) = 'difference between min and max'
            lvars(nvars+4) = 'probability of difference'
            nvars = nvars + 4
        end if
        do i=1,nvars
            if ( nz > 1 ) then
                ivars(1,i) = nz
            else
                ivars(1,i) = 0
            end if
            ivars(2,i) = 99
        end do
!           give correlations dates in 0-1
        if ( m1 == 0 ) then
            i = 0
        else
            i = 1
        end if
        j = m1-lag2
        if ( j <= 0 ) then
            j = j + nperyear*(1-j/nperyear)
        elseif ( j > nperyear ) then
            j = j - nperyear*((j-1)/nperyear)
        end if
        if ( nperyear > 12 ) then
!               rough approximation
            j = 1 + 12*(j-1)/nperyear
        end if
        if ( index(outfile,'.ctl') /= 0 ) then
            call writectl(outfile,datfile,nx,xx,ny,yy,nz,zz         &
                 ,1+(m2-m1)+(lag2-lag1),min(12,nperyear),i,j,3e33   &
                 ,title,nvars,vars,ivars,lvars,units)
        else
            call enswritenc(outfile,ncid,ntvarid,itimeaxis,ntmax,nx,xx &
                ,ny,yy,nz,zz,lz,1+(m2-m1)+(lag2-lag1),min(12,nperyear) &
                ,i,j,ltime,3e33,title,history,nvars,vars,ivars,lvars,svars &
                ,units,cell_methods,metadata,0,0)
        end if
!
!       write output field in GrADS or netCDF format
!
        print *,'writing output'
        do lag=lag2,lag1,-1
            do mo=m1,m2
                if ( m1 /= m2 ) then
                    m = mo-m1
                else
                    m = lag2-lag
                end if
                if ( index(outfile,'.ctl') /= 0 ) then
                    if ( lwrite ) then
                        print *,'writing records ',nvars*m+1,'-'    &
                              ,nvars*(m+1),' of fields ',m,' of size ',nx*ny*nz*recfa4
                        do jz=1,nz
                            do jy=1,ny
                                do jx = 1,nx
                                    if ( abs(r(jx,jy,jz,m)) <= 1 ) then
                                        i = i + 1
                                        d = d + abs(r(jx,jy,jz,m))
                                    end if
                                end do
                            end do
                        end do
                        if ( i > 0 ) then
                            print *,'there are ',i,' valid values in record ',m &
                                 ,' with mean value ',d/i
                        else
                            print *,'there are ',i,' valid values in record ',m
                        end if
                    end if
                    write(2,rec=nvars*m+1) (((r(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                    write(2,rec=nvars*m+2) (((prob(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                    write(2,rec=nvars*m+3) (((a(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                    write(2,rec=nvars*m+4) (((b(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                    write(2,rec=nvars*m+5) (((da(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                    write(2,rec=nvars*m+6) (((db(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                    if ( fitfunc == 1 ) then
                        write(2,rec=nvars*m+7) (((xn(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                        if ( lfitnoise ) then
                            write(2,rec=nvars*m+8) (((a1(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                            write(2,rec=nvars*m+9) (((da1(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                        elseif ( .not.composite ) then
                            write(2,rec=nvars*m+8) (((relregr(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                            write(2,rec=nvars*m+9) (((drelregr(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                        end if
                    else
                        write(2,rec=nvars*m+7) (((a1(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                        write(2,rec=nvars*m+8) (((da1(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                    end if
                    if ( nyrwindow > 0 ) then
                        write(2,rec=nvars*m+9) (((rmin(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                        write(2,rec=nvars*m+10) (((rmax(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                        write(2,rec=nvars*m+11) (((zdif(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                        write(2,rec=nvars*m+12) (((rprob(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                    end if
                else
!                       netCDF file
                    call writencslice(ncid,0,0,0,ivars(1,1),r(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                    call writencslice(ncid,0,0,0,ivars(1,2),prob(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                    call writencslice(ncid,0,0,0,ivars(1,3),a(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                    call writencslice(ncid,0,0,0,ivars(1,4),b(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                    call writencslice(ncid,0,0,0,ivars(1,5),da(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                    call writencslice(ncid,0,0,0,ivars(1,6),db(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                    call writencslice(ncid,0,0,0,ivars(1,7),xn(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                    if ( fitfunc == 1 ) then
                        if ( lfitnoise ) then
                            call writencslice(ncid,0,0,0,ivars(1,8),a1(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                            call writencslice(ncid,0,0,0,ivars(1,9),da1(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                        elseif ( .not.composite ) then
                            call writencslice(ncid,0,0,0,ivars(1,8),relregr(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                            call writencslice(ncid,0,0,0,ivars(1,9),drelregr(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                        end if
                    else
                        call writencslice(ncid,0,0,0,ivars(1,7),a1(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                        call writencslice(ncid,0,0,0,ivars(1,8),da1(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                    end if
                    if ( nyrwindow > 0 ) then
                       call writencslice(ncid,0,0,0,ivars(1,nvars-3),rmin(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                       call writencslice(ncid,0,0,0,ivars(1,nvars-2),rmax(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                       call writencslice(ncid,0,0,0,ivars(1,nvars-1),zdif(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                       call writencslice(ncid,0,0,0,ivars(1,nvars),rprob(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                    end if
                end if
            end do
        end do
        if ( index(outfile,'.ctl') == 0 ) then
            i = nf_close(ncid)
        end if
        if ( ncrossvalidate > 0 ) then
            title = 'data for subtractions of cross-validated regression'
            vars(1) = 'a'
            lvars(1) = 'constant in cross-validated regression'
            units(1) = inunits(1)
            open(1,file=trim(bbfile)//'_a.ctl')
            close(1,status='delete')
            call writectl(trim(bbfile)//'_a.ctl',                   &
                 trim(bbfile)//'_a.grd',nx,xx,ny,yy,nz,zz,          &
                 nperyear*(yr2-yr1+1),nperyear,yr1,1,3e33,title,    &
                 1,vars,ivars,lvars,units)
            open(1,file=trim(bbfile)//'_a.grd',access='direct',     &
                 form='unformatted',recl=nx*ny*nz*recfa4)
            irec = 0
            do i=yr1,yr2
                do j=1,nperyear
                    irec = irec + 1
                    write(1,rec=irec) (((aaa1(jx,jy,jz,j,i),jx=1,nx),jy=1,ny),jz=1,nz)
                end do
            end do
            close(1)
            vars(1) = 'b'
            lvars(1) = 'regr in cross-validated regression'
            units(1) = units(4)
            open(1,file=trim(bbfile)//'_b.ctl')
            close(1,status='delete')
            call writectl(trim(bbfile)//'_b.ctl',                   &
                 trim(bbfile)//'_b.grd',nx,xx,ny,yy,nz,zz,          &
                 nperyear*(yr2-yr1+1),nperyear,yr1,1,3e33,title,    &
                 1,vars,ivars,lvars,units)
            open(1,file=trim(bbfile)//'_b.grd',access='direct',     &
                 form='unformatted',recl=nx*ny*nz*recfa4)
            irec = 0
            do i=yr1,yr2
                do j=1,nperyear
                    irec = irec + 1
                    write(1,rec=irec) (((bbb1(jx,jy,jz,j,i),jx=1,nx),jy=1,ny),jz=1,nz)
                end do
            end do
            close(1)
        end if
    elseif ( ensemble .and. nens1 /= nens2 ) then ! subtract
        print *,'correlatefield: subtraction for ensembles not yet ready'
        call exit(-1)
    else                    ! subtract
!
!       subtract best fit from field and give this as output
!
        if ( lwrite ) print *,'Subtracting best fit from field'
        nrec = 0
        allocate(field2(nx,ny,nz,nperyear,firstyr:lastyr,nens1:nens2))
        field2 = 3e33
        do yr=firstyr,lastyr
!               skip month 0, which is everything together, if there is more
            if ( m2 /= 0 ) then
                mm = max(1,m1)
            else
                mm = m1
            end if
            do month=mm,m2
                if ( month == 0 ) then
                    j1 = 1
                    j2 = nperyear
                else
                    j1 = month
                    j2 = month + lsel - 1
                end if
                do lag=lag1,lag2
                    do jj=j1,j2
                        if ( fix2 ) then
                            j = jj+lag
                        else
                            j = jj
                        end if
                        call normon(j,yr,i,nperyear)
                        if ( i < yrbeg .or.i > yrend ) goto 910
                        if ( m1 /= m2 ) then
                            m = month-m1
                        elseif ( .not.fix2 ) then
                            m = lag2-lag
                        else
                            m = lag-lag1
                        end if
                        do kk=j-lag,j-lag+lsum-1
                            k = kk
                            call normon(k,i,ii,nperyear)
                            if ( ii < firstyr .or.ii > lastyr ) goto 910
                            if ( yr == yr1 ) then
                                print '(a,i3,i5,i3,i5,a)' &
                                     ,'Subtracting months',j,i,k,ii,' of point and field'
                            end if
                            do jz=1,nz
                                do jy=1,ny
                                    do jx=1,nx
                                        if ( b(jx,jy,jz,m) < 1e33 .and.                &
                                             field(jx,jy,jz,k,ii,nens1) < 1e33 .and.   &
                                             data(j,i,nens1) < 1e33 ) then
                                            if ( .false. .and. &
                                                 jx == (nx+1)/2.and.jy == (ny+1)/2.and.nz == (nz+1)/2 ) then
                                                print *,'field was ',field(jx,jy,jz,k,ii,nens1),k,ii,nens1
                                                print *,' b,data = ',b(jx,jy,jz,m),data(j,i,nens1)
                                            end if
                                            field2(jx,jy,jz,k,ii,nens1) = field(jx,jy,jz,k,ii,nens1) &
                                                 - b(jx,jy,jz,m)*data(j,i,nens1)
                                            if ( .false. .and. &
                                                 jx == (nx+1)/2.and.jy == (ny+1)/2.and.nz == (nz+1)/2 ) then
                                                print *,'field2 is ',field2(jx,jy,jz,k,ii,nens1)
                                            end if
                                        end if
                                    end do ! jy
                                end do ! jx
                            end do ! jz
                        end do ! k
910                     continue
                    end do   ! jj
                end do       ! lag
            end do           ! month
            do mo=1,nperyear
                nrec = nrec + 1
                if ( lwrite ) print *,'Writing new field for ',yr,mo
                write(2,rec=nrec) (((field2(jx,jy,jz,mo,yr,nens1),      &
                     jx=1,nx),jy=1,ny),jz=1,nz)
            end do
        end do               ! yr
        if ( index(outfile,'.ctl') /= 0 ) then
            title = trim(intitle)//' with the effect of '//     &
                  file(1+rindex(file,'/'):index(file,' ')-1)//  &
                  ' linearly subtracted'
            nvars = 1
            vars(1) = invars(1)
            ivars(1,1) = nz
            ivars(2,1) = 99
            lvars(1) = trim(inlvars(1))//' with '//             &
                 file(1+rindex(file,'/'):index(file,' ')-1)//   &
                 ' linearly subtracted'
            units(1) = inunits(1)
            call writectl(outfile,datfile,nx,xx,ny,yy,nz,zz     &
                  ,nrec,nperyear,firstyr,1,3e33,title,          &
                  nvars,vars,ivars,lvars,units)
        end if
    end if                   ! subtract
    close(2)
    call savestartstop(yrstart,yrstop)
!
!   error messages
!
    goto 999
901 write(0,*) 'correlatefield: error: fittime not yet implemented'
    call exit(-1)
909 write(0,*) 'correlatefield: error: array too small ',n
    call exit(-1)
920 write(0,*) 'correlatefield: error cannot open new '//       &
         'correlations file ',datfile(1:index(datfile,' ')-1)
    call exit(-1)
999 continue
end program correlatefield
