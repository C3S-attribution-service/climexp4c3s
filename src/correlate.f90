    program correlate
!       correlate station parameters from www.ncdc.noaa.gov or other
!       sources with all kinds of indices.  Currently included are:
!           SOI SLP index, from 1866 to 1997, from Phil Jones via Marleen
!           4 NINO SST indices, from 1950 to present, from
!               http://nic.fb4.noaa.gov:80/data/cddb/
!           NAO SLP index, from 1865 to 1996
!           SIDC sunspot index, from Zuerich,  http://www.astroinfo.ch/sunspot
!           SIDC sunspot cycle length, from Zuerich
!           time (in years)
!           any other file in .dat or .txt format (max 10)

!       GJvO hack, nov-1997, revised dec-1997, jan-1998, dec 1998, oct 1999,
!       jan-2000, added ensembles sep-2002, added running correlatons
!       2004, finally implemented ensembles correctly 2005
!
    implicit none
    integer,parameter :: nmc=1000,nvarmax=1
    include 'param.inc'
    include 'getopts.inc'
    integer :: i,ii,j,jj,k,l,m,n,j1,j2,if,im,ip,jm,jp,lag,year,month &
        ,yr,n1,ks,ks1,ks2,nperyear,ldir,iens,jens,ilen,jindx &
        ,ndup(0:npermax),ndum(0:npermax),validens(nensmax), &
        imens1(0:indxmx),imens(0:indxmx),yr1s,yr2s,nunequal &
        ,iunequal,mo,nfac(indxmx),yrstart,yrstop,mdata,nx,ny,nz,nt &
        ,nvars,nu,yrbg,yred,nensmx,iindex,lastmeta,init
    integer,allocatable :: yrmo(:,:),iarray(:)
    real,allocatable :: data(:,:,:),indx(:,:,:,:), &
        mcdata(:,:,:),mcindx(:,:,:,:), &
        ddata(:),dddata(:),dindx(:),aa(:),bb(:), &
        aa1(:,:,:),bb1(:,:,:)
    real :: anino(8),val12(npermax),adata,aindx,sxx,syy,yrmin(25), &
        yrmax(25),xmin,xmax,rise,fall,slength(25),ayr,r,prob, &
        absent,z,probd,sxy,df, &
        ss(nensmax*(indxmx+1)),minval, &
        result,dresult(-2:2),a(2),da(2,2),results(nmc),sig(1), &
        chi2,q,rmin,rmins(nmc),rmax,rmaxs(nmc),zdif,zdifs(nmc), &
        addfac(npermax,indxmx),sign,signmin,signmax,signdif, &
        alpha
    logical :: lfirstzero,lexist,lboot,ensindex,lprint,xwrap, &
        laddfile(indxmx),lbb1allocated
    logical,allocatable :: lfirst(:)
    parameter (absent=3e33)
    character :: line*1024,file*1024,string*10,dir*256,ensfile*1024, &
        newunits*60,dum1*40,dum2*40,varorg*80,unitsorg*60
    character :: var*80,units*60,lvar*120,svar*120,history*50000,metadata(2,100)*2000
    character :: ivar(indxmx)*80,iunits(indxmx)*60,ilvar(indxmx)*120,isvar(indxmx)*120, &
        ihistory(indxmx)*50000,imetadata(2,100,indxmx)*1000

    character(4) :: runs(3,2)
    data runs /'rmin','rmax','zdif', 'bmin','bmax','bdif'/

!   check arguments

    lwrite = .false. 
    n = command_argument_count()
    if ( n < 1 ) then
        print *,'usage: correlate datafile index [lag n[:m]] '// &
            '[month m[:n] [sum|ave|max|min|sel m] [log|sqrt|rank]' &
            //' [begin yr] [end yr] [ks plot|cut-off] [detrend]', &
            '[diff [nyr]] [plot file] [lt maxindex] [gt minindex]' &
            //'[decor n] [runcorr|runregr nyr outfile [random '// &
            'series|index] [noise white|red]]'
            call exit(-1)
    endif
    if ( indxmx /= 19 ) then
        write(0,*) 'expecting indxmx=19, not ',indxmx
        call exit(-1)
    endif

!   do not make them much bigger then needed...

    call get_command_argument(1,file)
    call getfileunits(file,nx,ny,nz,nt,nperyear,nvarmax,nvars,var,units,newunits, &
        lvar,svar,xwrap,lwrite)
    if ( nperyear <= 0 ) then
        write(0,*) 'correlate: error: the file ',trim(file),' is not readable'
        call exit(-1)
    end if
!       do not let the arrays get too large
    if ( nperyear <= 12 ) then
        yrbg = yrbeg
        nensmx = nensmax
    else if ( nperyear <= 366 ) then
        yrbg = 1700
        nensmx = 10
    else
        yrbg = 1900
        nensmx = 1
    end if
    yred = yrend
    mdata = nperyear*(yred-yrbg+1)*(nensmx+1)
    if ( lwrite ) print *,'allocating ',mdata,' reals: ',nperyear,yrbg,yred,nensmx
    allocate(data(nperyear,yrbg:yred,0:nensmx))
    if ( lwrite ) print *,'allocating ',mdata*indxmx,' reals: ',nperyear,yrbg,yred,nensmx,indxmx
    allocate(indx(nperyear,yrbg:yred,0:nensmx,indxmx))
    allocate(ddata(mdata))
    allocate(dddata(mdata))
    allocate(dindx(mdata))
    allocate(lfirst(mdata))
    allocate(iarray(mdata))
    if ( lwrite ) print *,'allocating ',mdata*3,' reals'
    allocate(yrmo(3,mdata))

!   init

    imens1 = 0
    imens = 0
    yrstart = yred
    yrstop  = yrbg
    call getenv('DIR',dir)
    ldir = len_trim(dir)
    if ( ldir <= 1 ) then
        call getenv('HOME',dir)
        ldir = len_trim(dir)
        dir = dir(1:len_trim(dir))//'/climexp/'
        ldir = len_trim(dir)
    elseif ( dir(ldir:ldir) /= '/' ) then
        ldir = ldir + 1
        dir(ldir:ldir) = '/'
    endif
    iunequal = 0
    nunequal = 1

!   read data from station file downloaded from
!   http://www.ncdc.noaa.gov/ghcn/ghcnV1.CLIMVIS.html

    call get_command_argument(1,file)
    if ( lwrite ) print *,'reading data file ',trim(file)
    call readensseriesmeta(file,data,nperyear,yrbg,yred,nensmx, &
        n,imens1(0),imens(0),var,units,lvar,svar,history,metadata, &
        lstandardunits,lwrite)
    if ( n /= nperyear ) call exit(-1)
    if ( imens(0) > 0 ) ensemble = .true. 

!   process options

    n = command_argument_count()
    call getopts(2,n,nperyear,yrbg,yred, .true. ,imens1(0),imens(0))
    if ( lstandardunits ) then
        ! not right first time...
        if ( lwrite ) print *,'calling makestandardseries'
        varorg = var
        unitsorg = units
        do iens = nens1,nens2
            var = varorg
            units = unitsorg
            call makestandardseries(data(1,yrbg,iens),nperyear, &
                yrbg,yred,nperyear,var,units,lwrite)
        end do
    end if
    if ( minfac < 0 ) minfac = 0.75
    if ( imens(0) > 0 ) then
        print '(a,i4,a,i4,a)','# taking ensemble members ',nens1,' to ',nens2,' of series'
        ensemble = .true. 
    endif
    if ( lsubtract ) write(0,*) 'correlate: subtract not yet implemented'
    do i=1,indxmx
        if ( lincl(i) ) goto 100
    enddo
!   demand at least one series...
    write(0,*) 'correlate: please select an index to correlate with'
    call exit(-1)
100 continue
!**        if ( lsum.gt.1 .and. lsel.gt.1 ) goto 915
    if ( dump ) then
        write(10,'(6a)') '# ',trim(var),' [',trim(units),']'
        if ( logscale ) write(10,'(a)') '# logarithmic plot'
        if ( sqrtscale ) write(10,'(a)') '# sqrt plot'
    endif
    if ( ncrossvalidate > 0 ) then
        allocate(aa(mdata))
        allocate(bb(mdata))
    end if

!   get SOI,NINO,NAO (really old set-up)

    do k=1,indxmx
        call makeabsent(indx(1,yrbg,0,k),nperyear,yrbg,yred)
        imens1(k) = 0
        imens(k) = 0
    enddo

!   get the SOI data

    if ( lincl(1) ) then
        if ( lwrite ) print *,'reading SOI file soi.dat'
        call readseriesmeta(trim(dir)//'CRUData/soi.dat', &
            indx(1,yrbg,0,1),nperyear,yrbg,yred,n,ivar(1),iunits(1),ilvar(1),isvar(1), &
            ihistory(1),imetadata(1,1,1),.false.,lwrite)
        if ( n /= nperyear ) goto 916
    endif

!   get NINO indices from ERSST v4

    if ( lincl(2) ) then
        if ( lwrite ) print *,'reading NCDCData/ersst_nino12a.dat'
        call readseriesmeta(trim(dir)//'NCDCData/ersst_nino12.dat', &
            indx(1,yrbg,0,2),nperyear,yrbg,yred,n,ivar(1),iunits(2),ilvar(2),isvar(2), &
            ihistory(2),imetadata(1,1,2),.false.,lwrite)
        if ( n /= nperyear ) goto 916
    end if
    if ( lincl(3) ) then
        if ( lwrite ) print *,'reading NCDCData/ersst_nino3a.dat'
        call readseriesmeta(trim(dir)//'NCDCData/ersst_nino3a.dat', &
            indx(1,yrbg,0,3),nperyear,yrbg,yred,n,ivar(3),iunits(3),ilvar(3),isvar(3), &
            ihistory(3),imetadata(1,1,3),.false.,lwrite)
        if ( n /= nperyear ) goto 916
    end if
    if ( lincl(4) ) then
        if ( lwrite ) print *,'reading NCDCData/ersst_nino4a.dat'
        call readseriesmeta(trim(dir)//'NCDCData/ersst_nino4a.dat', &
            indx(1,yrbg,0,4),nperyear,yrbg,yred,n,ivar(4),iunits(4),ilvar(4),isvar(4), &
            ihistory(4),imetadata(1,1,4),.false.,lwrite)
        if ( n /= nperyear ) goto 916
    end if
    if ( lincl(5) ) then
        if ( lwrite ) print *,'reading NCDCData/ersst_nino3.4a.dat'
        call readseriesmeta(trim(dir)//'NCDCData/ersst_nino3.4a.dat', &
            indx(1,yrbg,0,5),nperyear,yrbg,yred,n,ivar(5),iunits(5),ilvar(5),isvar(5), &
            ihistory(5),imetadata(1,1,5),.false.,lwrite)
        if ( n /= nperyear ) goto 916
    end if

!   and the NAO data

    if ( lincl(6) ) then
        call readseriesmeta(trim(dir)//'CRUData/nao.dat', &
            indx(1,yrbg,0,6),nperyear,yrbg,yred,n,ivar(6),iunits(6),ilvar(6),isvar(6), &
            ihistory(6),imetadata(1,1,6),.false.,lwrite)
        if ( nperyear /= n ) goto 916
    endif

!   and the sunspot data, 1749-1991

    if ( lincl(7) ) then
        call readseriesmeta(trim(dir)//'SIDCData/sunspots.dat', &
            indx(1,yrbg,0,7),nperyear,yrbg,yred,n,ivar(7),iunits(7),ilvar(7),isvar(7), &
            ihistory(7),imetadata(1,1,7),.false.,lwrite)
        if ( nperyear /= n ) goto 916
    endif
    if ( lincl(8) ) then
        write(0,*) 'this option is no longer supported'
        write(*,*) 'this option is no longer supported'
        call exit(-1)
    endif
    if ( lincl(9) ) then
        ! time
        ivar(9) = 'time'
        iunits(9) = 'yr'
        isvar(9) = 'time'
        ilvar(9) = 'years from 2000'
        ihistory(9) = ' '
        imetadata(:,:,9) = ' '
        do i=yrbg,yred
            do j=1,nperyear
                indx(j,i,0,9) = i-2000 + (j-0.5)/nperyear
            enddo
        enddo
    endif

!   get my own file

    do i=10,indxuse
        if ( lwrite ) print *,'Opening index file ' ,trim(indexfiles(i))
        if ( indexfiles(i) == ' ' ) then
            write(0,*) 'correlate: error: empty file name ',i
            call exit(-1)
        end if
        call readensseriesmeta(indexfiles(i),indx(1,yrbg,0,i),nperyear &
            ,yrbg,yred,nensmx,n,imens1(i),imens(i),ivar(i),iunits(i),ilvar(i),isvar(i), &
            ihistory(i),imetadata(1,1,i),lstandardunits,lwrite)
        if ( lwrite ) then
            print *,'nperyear from this file = ',n
            print *,'metadata ',i
            do j=1,100
                if ( imetadata(1,j,i) == ' ' ) exit
                print *,j,trim(imetadata(1,j,i)),' :: ',trim(imetadata(2,j,i))
            end do
        end if
        if ( n /= nperyear ) goto 916
        if ( imens(i) > 0 ) then
            ensindex = .true. 
            if ( .not. ensemble ) then
                nens1 = max(nens1,imens1(i))
                if ( nens2 > 0 ) then
                    nens2 = min(nens2,imens(i))
                else
                    nens2 = imens(i)
                endif
                write(*,'(a,i2,a,i3,a,i3)') '# Found ensemble ',i &
                    ,' from ',imens1(i),' to ',imens(i)
            else
                if ( imens(i) /= imens(0) ) then
                    iunequal = iunequal + 1
                    if ( iunequal >= nunequal ) then
                        nunequal = 2*nunequal
                        write(0,'(a,2i3,a)') 'warning: unequal sizes ensembles: ', &
                            imens(0),imens(i),' using smallest'
                        write(*,'(a,2i3,a)') '# warning: unequal sizes ensembles: ', &
                            imens(0),imens(i),' using smallest'
                    end if
                    nens2 = min(imens(i),imens(0))
                    imens(0) = nens2
                    imens(i) = nens2
                endif
                imens(i) = min(nens2,imens(i))
                do k=1,i-1
                    if ( lincl(k) .and. &
                    imens(k) > 0 .and. &
                    imens(i) /= imens(k) ) then
                        iunequal = iunequal + 1
                        if ( iunequal >= nunequal ) then
                            nunequal = 2*nunequal
                            write(0,'(a,2i3,a)') 'warning: unequal ensembles: ' &
                                ,imens(k),imens(i),' using smallest'
                            write(*,'(a,2i3,a)') '# warning: unequal ensembles: ' &
                                ,imens(k),imens(i),' using smallest'
                        end if
                        imens(k) = min(imens(i),imens(k))
                        imens(i) = imens(k)
                        nens2 = imens(k)
                        if ( ensemble ) imens(0) = imens(k)
                    endif
                enddo
            endif
        endif
    enddo

!   merge metadata

    do lastmeta = 1,100
        if ( metadata(1,lastmeta) == ' ' ) exit
    end do
    lastmeta = lastmeta - 1
    iindex = 0
    do i=1,indxuse
        if ( lincl(i) .or. i >= 10 ) iindex = iindex + 1
    end do
    do i=1,indxuse
        if ( lincl(i) .or. i >= 10 ) then
            if ( iindex == 1) then
                string = 'index_'
            else
                write(string,'(a,i2.2,a)') 'index_',i,'_'
            end if
            call merge_metadata(metadata,lastmeta,imetadata(1,1,i),' ',ihistory(i),trim(string))
        end if
    end do
    if ( dump ) then
        call printmetadata(10,file,' ','correlation analysis of',history,metadata)
    end if

!   am I being called as addseries?

    call get_command_argument(0,line)
    if ( index(line,'addseries ') /= 0 ) then
        if ( lwrite ) print *,'adding series'
        if ( nens2 > 0 ) then
            write(0,*) 'error: cannot add ensembles yet'
            call exit(-1)
        endif
        if ( lag1 /= lag2 ) then
            write(0,*) 'error: cannot use range of lags',lag1,lag2
            call exit(-1)
        endif
        lag = lag1
    
!       construct name
    
        i=1
        do k=1,indxuse
            if ( lincl(k) ) then
                plotfile(i:) = strindx(k)
                j = index(strindx(k),'.')
                if ( j > 0 ) then
                    i = i+j-1
                else
                    i = i + len_trim(strindx(k))
                endif
                plotfile(i:) = '+'
                i = i+1
            endif
        enddo
        i = i-1
        plotfile(i:i) = ' '
        print '(a)',plotfile(1:i-1)
    
!       construct file name
    
        call get_command_argument(1,plotfile)
        do i=len(plotfile),1,-1
            if ( plotfile(i:i) == '/' ) goto 510
        enddo
    510 continue
        i = i+1
        plotfile(i:) = 'iadded'
        i = i+6
        ii = 0
    520 continue
        ii = ii + 1
        if ( ii < 10 ) then
            write(plotfile(i:),'(i1,a)') ii,'.dat'
        elseif ( ii < 100 ) then
            write(plotfile(i:),'(i2,a)') ii,'.dat'
        elseif ( ii < 1000 ) then
            write(plotfile(i:),'(i3,a)') ii,'.dat'
        else
            write(0,*) 'addseries: error: cannot open output file ' &
                ,plotfile(1:len_trim(plotfile))
            write(*,*) 'addseries: error: cannot open output file ' &
                ,plotfile(1:len_trim(plotfile))
            call exit(-1)
        endif
        open(99,file=plotfile,status='new',err=520)
        print '(a)',plotfile(1:len_trim(plotfile))
    
!       get the coefficients from the environment
!       I assume they are stored as FORM_a1,FORM_a2,...
    
        jindx = 0
        laddfile = .false. 
        nfac = 1
        do k=1,indxuse
            if ( lincl(k) ) then
                if ( k == 10 ) then
                    ! assume the first user-defined file is the main one...
                    call copyheader(indexfiles(k),99)
                end if
                jindx = jindx + 1
                if ( jindx < 10 ) then
                    write(string,'(a,i1)') 'FORM_a',jindx
                    i = 7
                else
                    write(string,'(a,i2)') 'FORM_a',jindx
                    i = 8
                endif
                call getenv(string(1:i),line)
                if ( line == ' ' ) then
                    write(0,*) 'addseries: error: environment variable ',string(1:i), &
                        ' corresponding to ',trim(strindx(k)),' not set'
                    call exit(-1)
                end if
                if ( lwrite ) print *,'found env ',string(1:i)  ,' with content ',trim(line)
                if ( line(1:4) /= 'file' ) then
                ! numerical values
                    nfac(k) = 1
                    do i=1,len(line)-1
                        if ( line(i:i) == '+' ) line(i:i) = ' '
                        if ( line(i:i+1) == '--' ) &
                        line(i:i+1) = '  '
                        if ( line(i:i) == ':' ) then
                            line(i:i) = ' '
                            nfac(k) = nfac(k) + 1
                        endif
                    enddo
                    if ( nfac(k) == nperyear+1 ) nfac(k) = nperyear
                    if ( nfac(k) /= 1 .and. nfac(k) /= nperyear ) then
                        write(0,*) 'addseries: error: can only use ' &
                            //'#addfac 1 or ',nperyear,', not ',nfac(k)
                        write(*,*) 'addseries: error: can only use ' &
                            //'#addfac 1 or ',nperyear,', not ',nfac(k)
                        call exit(-1)
                    endif
                    read(line,*,err=600,end=600) &
                    (addfac(j,k),j=1,nfac(k))
                    if ( k == 10 .and. oper == '+' .and. lsum > 1 ) then
                        do j=1,nfac(k)
                            addfac(j,k) = addfac(j,k)*lsum
                        enddo
                    endif
                    if ( lwrite ) then
                        print *,'adding ',strindx(k),' with weight ',(addfac(j,k),j=1,nfac(k))
                        print *,'using lag ',lag,fix2
                    endif
                    write(99,'(3a,366f10.6)') '# added ',strindx(k),' with weight', &
                        (addfac(j,k),j=1,nfac(k))
                else
                    !   the coefficients are in a file
                    if ( .not. allocated(bb1) ) then
                        allocate(aa1(nperyear,yrbg:yred,indxmx))
                        allocate(bb1(nperyear,yrbg:yred,indxmx))
                        aa1 = 3e33
                        bb1 = 3e33
                        ilen = index(bbfile,'.dat') - 1
                        if ( ilen <= 0 ) ilen = len_trim(bbfile)
                        call readseries(bbfile(:ilen)//'_a.dat', &
                            aa1(1,yrbg,k),nperyear,yrbg,yred,i, &
                            dum1,dum2,.false.,lwrite)
                        if ( i /= nperyear ) then
                            write(0,*) 'error: found nperyear = ',i, &
                                ' in ',trim(line),', expected ',nperyear
                            call exit(-1)
                        end if
                    end if
                    if ( jindx < 10 ) then
                        write(bbfile(ilen+1:),'(a,i1,a)') '_b',jindx,'.dat'
                    else
                        write(bbfile(ilen+1:),'(a,i2,a)') '_b',jindx,'.dat'
                    end if
                    call readseries(bbfile,bb1(1,yrbg,k),nperyear &
                        ,yrbg,yred,i,dum1,dum2,.false.,lwrite)
                    write(99,'(4a)') '# added ',strindx(k), &
                        ' with weights from file ',trim(bbfile)
                    laddfile(k) = .true. 
!                   if there is only one value defined, set all
!                   values for the year afterwards to that value
                    n = 0
                    m = -1
                    do j=1,nperyear
                        do i=yrbg,yred
                            if ( bb1(j,i,k) < 1e33 ) then
                                n = n + 1
                                m = j
                                exit
                            end if
                        end do
                    end do
                    if ( n == 1 .and. nperyear > 1 ) then
                        do i=yrbg,yred
                            do j=1,nperyear
                                if ( j > m ) then
                                    aa1(j,i,k) = aa1(m,i,k)
                                    bb1(j,i,k) = bb1(m,i,k)
                                elseif ( j < m .and. i > 1 ) then
                                    aa1(j,i,k) = aa1(m,i-1,k)
                                    bb1(j,i,k) = bb1(m,i-1,k)
                                end if
                            end do
                        end do
                        n = nperyear
                    end if
                    nfac(k) = n
                end if
                if ( k /= 10 .and. lag /= 0 ) then
                    write(99,'(a,i3,a,i3)') '# used lag ',lag,' for index ',k
                endif
            600 continue
            else
                addfac(1,k) = 0
            endif
        enddo
!
!       transform time series if requested
!
        if ( logscale ) then
            if ( lwrite ) print '(a,2i3)','# Taking log of series '
            call takelog(indx(1,yrbg,0,10),nperyear,nperyear,yrbg,yred)
        endif
        if ( logfield ) then
            do k=1,indxuse
                if ( lincl(k) .and. k /= 10 ) then
                    if ( lwrite ) print '(a,2i3)','# Taking log of index ',k
                    call takelog(indx(1,yrbg,iens,k),nperyear,nperyear,yrbg,yred)
                endif
            enddo
        endif
        if ( sqrtscale ) then
            if ( lwrite ) print '(a,2i3)','# Taking sqrt of series '
            call takesqrt(indx(1,yrbg,0,10),nperyear,nperyear,yrbg,yred)
        endif
        if ( sqrtfield ) then
            do k=1,indxuse
                if ( lincl(k) .and. k /= 10 ) then
                    if ( lwrite ) print '(a,2i3)','# Taking sqrt of index ',k
                    call takesqrt(indx(1,yrbg,iens,k),nperyear,nperyear,yrbg,yred)
                endif
            enddo
        endif
!
!       take anomalies of indices to preserve mean of series
!
        do k=1,indxuse
            if ( lincl(k) .and. k /= 10 .and. .not.laddfile(k) ) then
                call anomal(indx(1,yrbg,iens,k),nperyear,nperyear,yrbg,yred,yrbg,yred)
            end if
        end do        
!
!       perform the subtraction
!
        do yr=yrbg,yred
            n = 0
            do jj=1,nperyear
                if ( fix2 ) then
                    j = jj+lag
                else
                    j = jj
                endif
                call normon(j,yr,i,nperyear)
                if ( i < yrbg .or. i > yred ) cycle
                m = j-lag
                call normon(m,i,ii,nperyear)
                if ( ii < yrbg .or. ii > yred ) cycle
                data(j,i,0) = 0
                do k=1,indxuse
                    if ( lincl(k) .and. ( laddfile(k) .or. &
                         addfac(min(j,nfac(k)),k) /= 0 ) ) then
                        if ( k == 10 ) then ! no lag
                            if ( data(j,i,0) < 0.9*absent .and. &
                                 indx(j,i,0,k) < 0.9*absent) then
                                if ( laddfile(k) ) then
                                    if ( aa1(j,i,k) < 1e33 .and. &
                                    bb1(j,i,k) < 1e33 ) then
                                        data(j,i,0) = data(j,i,0) - aa1(j,i,k) - bb1(j,i,k)*indx(j,i,0,k)
                                    else
                                        data(j,i,0) = absent
                                    end if
                                else
                                    data(j,i,0) = data(j,i,0) + addfac(min(j,nfac(k)),k)*indx(j,i,0,k)
                                end if
                            else
                                data(j,i,0) = absent
                            endif
                        else
                            if ( data(j,i,0) < 0.9*absent .and. &
                                 indx(m,ii,0,k) < 0.9*absent) then
                                if ( laddfile(k) ) then
                                    if ( aa1(j,i,k) < 1e33 .and. &
                                         bb1(j,i,k) < 1e33 ) then
                                        data(j,i,0) = data(j,i,0) - aa1(j,i,k) - bb1(j,i,k)*indx(m,ii,0,k)
                                    else
                                        data(j,i,0) = absent
                                    end if
                                else
                                    data(j,i,0) = data(j,i,0) + addfac(min(j,nfac(k)),k)*indx(m,ii,0,k)
                                end if
                            else
                                data(j,i,0) = absent
                            endif
                        endif
                    endif
                enddo
            end do ! nperyear
        end do ! yr
!
!       transform back (only one ensemble member)
!
        if ( logscale ) then
            call takeexp(data,nperyear,nperyear,yrbg,yred)
        endif
        if ( sqrtscale ) then
            call takesquare(data,nperyear,nperyear,yrbg,yred)
        end if
!
!       output
!
        call printmetadata(99,file,' ','time series with regressions subtractief',history,metadata)
        call printdatfile(99,data,nperyear,nperyear,yrbg,yred)
        goto 999
    endif

!   monthly diffs

    if ( mdiff > 0 ) then
        if ( lwrite ) print '(a)','# Taking monthly anomalies of series'
        do iens=imens1(0),imens(0)
            call mdiffit(data(1,yrbg,iens),nperyear,nperyear,yrbg,yred,mdiff)
        enddo
    endif
    if ( mdiff2 > 0 ) then
        do k=1,indxuse
            if ( lincl(k) ) then
                if ( lwrite ) print '(a,i3)','# Taking monthly anomalies of index ',k
                do iens=imens1(k),imens(k)
                    call mdiffit(indx(1,yrbg,iens,k),nperyear,nperyear,yrbg,yred,mdiff2)
                enddo
            endif
        enddo
    endif

!   sum

    if ( lsum > 1 ) then
        if ( lwrite ) print '(a,i3)','# Summing series ',lsum
        do iens=imens1(0),imens(0)
            call sumit(data(1,yrbg,iens),nperyear,nperyear,yrbg,yred,lsum,oper)
        enddo
    endif
    if ( lsum2 > 1 ) then
        do k=1,indxuse
            if ( lincl(k) ) then
                if ( lwrite ) print '(a,2i3)','# Averaging index ',k,lsum2
                do iens=imens1(k),imens(k)
                    call sumit(indx(1,yrbg,iens,k),nperyear,nperyear,yrbg,yred,lsum2,'v')
                enddo
            endif
        enddo
    endif

!   logscale

    if ( logscale ) then
        do iens=imens1(0),imens(0)
            if ( lwrite ) print '(a,2i3)','# Taking log of series '
            call takelog(data(1,yrbg,iens),nperyear,nperyear,yrbg,yred)
        enddo
    endif
    if ( logfield ) then
        do k=1,indxuse
            if ( lincl(k) ) then
                if ( lwrite ) print '(a,2i3)','# Taking log of index ',k
                do iens=imens1(k),imens(k)
                    call takelog(indx(1,yrbg,iens,k),nperyear &
                    ,nperyear,yrbg,yred)
                enddo
            endif
        enddo
    endif

!   sqrtscale

    if ( sqrtscale ) then
        do iens=imens1(0),imens(0)
            if ( lwrite ) print '(a,2i3)','# Taking sqrt of series '
            call takesqrt(data(1,yrbg,iens),nperyear,nperyear,yrbg,yred)
        enddo
    endif
    if ( sqrtfield ) then
        do k=1,indxuse
            if ( lincl(k) ) then
                if ( lwrite ) print '(a,2i3)' &
                ,'# Taking sqrt of index ',k
                do iens=imens1(k),imens(k)
                    call takesqrt(indx(1,yrbg,iens,k),nperyear,nperyear,yrbg,yred)
                enddo
            endif
        enddo
    endif

!   find begin,end of valid data

    if ( lwrite ) print *,'adjusting yr1,yr2 from ',yr1,yr2
    do yr=yr1,yr2
        do iens=imens1(0),imens(0)
            do mo=1,nperyear
                if ( data(mo,yr,iens) < 1e33 ) then
                    goto 710
                endif
            enddo
        enddo
    enddo
710 continue
    if ( fix2 ) then
        yr1 = max(yr1,yr+min(0,(lag1-nperyear+1)/nperyear))
    else
        yr1 = max(yr1,yr)
    endif
    do yr=yr1,yr2
        do k=1,indxuse
            if ( lincl(k) ) then
                do iens=imens1(k),imens(k)
                    do mo=1,nperyear
                        if ( indx(mo,yr,iens,k) < 1e33 ) then
                            goto 720
                        endif
                    enddo
                enddo
            endif
        enddo
    enddo
720 continue
    if ( fix2 ) then
        yr1 = max(yr1,yr)
    else
        yr1 = max(yr1,yr-max(0,(lag2+nperyear-1)/nperyear))
    endif

    do yr=yr2,yr1,-1
        do iens=imens1(0),imens(0)
            do mo=1,nperyear
                if ( data(mo,yr,iens) < 1e33 ) then
                    goto 730
                endif
            enddo
        enddo
    enddo
    print '(a)','# found no valid data in data'
730 continue
    if ( fix2 ) then
        yr2 = min(yr2,yr+max(0,(lag2+nperyear-1)/nperyear))
    else
        yr2 = min(yr2,yr)
    endif
    do yr=yr2,yr1,-1
        do k=1,indxuse
            if ( lincl(k) ) then
                do iens=imens1(k),imens(k)
                    do mo=1,nperyear
                        if ( indx(mo,yr,iens,k) < 1e33 ) then
                            goto 740
                        endif
                    enddo
                enddo
            endif
        enddo
    enddo
    print '(a,i3)','# found no valid data in index ',k
740 continue
    if ( fix2 ) then
        yr2 = max(yr2,yr)
    else
        yr2 = max(yr2,yr-min(0,(lag1-nperyear+1)/nperyear))
    endif
    if ( lwrite ) print *,'                    to ',yr1,yr2

!   detrending

    if ( ldetrend ) then
        if ( lwrite ) print *,'Detrending series'
        do iens=imens1(0),imens(0)
            call detrend(data(1,yrbg,iens),nperyear,nperyear, yrbg,yred,yr1,yr2,m1,m2,lsel)
        enddo
        do k=1,indxuse
            if ( lincl(k) .and. k /= 9 ) then
                if ( lwrite ) print *,'Detrending index ',k
                do iens=imens1(k),imens(k)
                    if ( lag1 == 0 .and. lag2 == 0 .or. m1 == 0 .or. &
                    lsel == 12 ) then
                        call detrend(indx(1,yrbg,iens,k),nperyear, &
                            nperyear,yrbg,yred,yr1,yr2,m1,m2,lsel)
                    else
                        call detrend(indx(1,yrbg,iens,k),nperyear, &
                            nperyear,yrbg,yred,yr1,yr2,1,nperyear,lsel)
                    endif
                enddo
            endif
        enddo
    endif
    if ( ndiff /= 0 ) then
        if ( lwrite ) print *,'Taking differences - series'
        do iens=imens1(0),imens(0)
            call ndiffit(data(1,yrbg,iens),nperyear,nperyear,yrbg,yred,ndiff,minfacsum)
        enddo
        if ( ndiff < 0 .and. lnooverlap ) then
            call dooverlap(data,nperyear,yrbg,yred,imens1(0),imens(0),nperyear,ndiff)
        end if
        do k=1,indxuse
            if ( lincl(k) .and. k /= 9 ) then
                if ( lwrite ) print * &
                ,'Taking differences - index ',k
                do iens=imens1(k),imens(k)
                    call ndiffit(indx(1,yrbg,iens,k),nperyear,nperyear,yrbg,yred,ndiff,minfacsum)
                enddo
            endif
        enddo
    endif

!   bias corrections

    n = 0
    do k=1,indxuse
        if ( lincl(k) ) n = n + 1
    end do
    if ( n == 1 ) then
    ! no different bias corrections for different models yet
        if ( debias == 1 ) then
            print '(a)','# correcting for bias in mean'
            call debiasmean(indx(1,yrbg,0,indxuse),data,nperyear, &
                nperyear,yrbg,yred,yr1,yr2,imens1(0),imens(0),var,lwrite)
        elseif ( debias == 2 ) then
            print '(a)','# correcting for bias in mean and variance'
            call debiasvar(indx(1,yrbg,0,indxuse),data,nperyear, &
                nperyear,yrbg,yred,yr1,yr2,imens1(0),imens(0),var,lwrite)
        elseif ( debias == 3 ) then
            print '(a)','# correcting for bias in distribution'
            call debiasall(indx(1,yrbg,0,indxuse),data,nperyear, &
                nperyear,yrbg,yred,yr1,yr2,imens1(0),imens(0),var,lwrite)
        endif
        if ( debias > 0 ) then
            ncrossvalidate = 0 ! is already included
        end if
    end if

!   anomalies - necessary if we fit time derivative with relaxation

    if ( anom .or. nfittime > 0 .and. ndiff <= 0 ) then
        if ( lwrite ) print *,'Taking anomalies - data'
        do iens=imens1(0),imens(0)
            call anomal(data(1,yrbg,iens),nperyear,nperyear,yrbg,yred,yr1,yr2)
        enddo
        do k=1,indxuse
            if ( lincl(k) ) then
                if ( lwrite ) print *,'Taking anomalies - index ',k
                do iens=imens1(k),imens(k)
                    call anomal(indx(1,yrbg,iens,k),nperyear,nperyear,yrbg,yred,yr1,yr2)
                enddo
            endif
        enddo
    endif

!   anomalies wrt ensemble mean

    if ( lensanom .and. imens1(0) /= imens(0) ) then
        if ( lwrite ) print *,'Taking anomalies wrt ensemble mean - data'
        call anomalensemble(data,nperyear,nperyear,yrbg, &
            yred,yr1,yr2,max(nens1,imens1(0)),min(nens2,imens(0)))
    endif
    do k=1,indxuse
        if ( lincl(k) .and. lensanom .and. imens1(k) /= imens(k)) &
        then
            if ( lwrite ) print *,'Taking anomalies wrt ensemble mean - index ',k
            call anomalensemble(indx(1,yrbg,0,k),nperyear,nperyear, &
                yrbg,yred,yr1,yr2,max(nens1,imens1(k)),min(nens2,imens(k)))
        endif
    enddo

!   copy ensemble members so that there is the same
!   number of valid ones at every time step

    do j=0,nperyear
        ndup(j) = 0
    enddo
    if ( ensemble .and. lmakeensfull ) then
        do k=1,indxmx
            if ( ensemble .and. imens(k) > 0 ) then
                do yr=yrbg,yred
                    do mo=1,nperyear
                        do iens=nens1,nens2
                            if ( data(mo,yr,iens) < 1e33 .neqv. &
                            indx(mo,yr,iens,k) < 1e33 ) then
                                if ( lwrite) print *,'correlate: warning: inconsistent ', &
                                'ensembles ',mo,yr,iens,k,data(mo,yr,iens),indx(mo,yr,iens,k)
                                data(mo,yr,iens) = 3e33
                                indx(mo,yr,iens,k) = 3e33
                            endif
                        enddo
                    enddo
                enddo
            endif
        enddo
        if ( lwrite ) print *,'correlate: calling makeensfull for data'
        call makeensfull(ndup,nperyear,data,nperyear,yrbg,yred &
            ,nens1,nens2,validens,lwrite)
        call random_number(prob)
    endif
    if ( ensindex .and. lmakeensfull ) then
        do k=1,indxmx
            if ( imens(k) > 0 ) then
                call makeensfull(ndum,nperyear,indx(1,yrbg,0,k) &
                    ,nperyear,yrbg,yred,nens1,nens2,validens,lwrite)
                call random_number(probd)
                if ( prob /= probd ) then
                    write(0,*) 'error: ensembles not filled identically: ',prob,probd
                endif
            endif
        enddo
    endif

!   correlate!

    call printcorrheader
    do month=m1,m2
        call getj1j2(j1,j2,month,nperyear, .true. )
        lboot = lbootstrap
        if ( ensemble ) then
            call checkdup(lboot,ndup,nperyear,j1,j2)
        endif
    
!       fill linear arrays without absent values and compute r
        do lag=lag1,lag2
            do k=1,indxuse
                if ( .not. lincl(k) ) goto 800
                call perc2cut(lag,k,j1,j2,nperyear,imens1,imens &
                    ,indxmx,indx,data,nperyear,yrbg,yred,nensmx)
                if ( lwrite ) print *,'index ',strindx(k),k
                if ( dump ) write(10,'(2a)') '# ',strindx(k)
                call filllinarray(dindx,ddata,lfirst,dddata,mdata,n &
                    ,j1,j2,lag,k,nperyear,imens,indxmx,indx,data &
                    ,nperyear,yrbg,yred,nensmx,yrstart,yrstop,yrmo)
                if ( n < minnum ) then
                    if ( lwrite ) print *,'correlate: not enough points: ',n,minnum
                    goto 800
                endif
                call printcorr(dindx,ddata,lfirst,dddata,yrmo,n &
                    ,ndup(0),j1,j2,month,nperyear,lag,strindx(k) &
                    ,lboot, .true. ,result,dresult,prob)
            
!               running correlations
            
                if ( nyrwindow > 0 ) then
                    call getruncorr(dindx,ddata,lfirst,dddata,mdata &
                        ,j1,j2,lag,k,month,nperyear,imens,indxmx &
                        ,indx,data,nperyear,yrbg,yred,nensmx &
                        ,ndup(0),strindx(k),lboot,.true. &
                        ,rmin,rmax,zdif)
                
!                   Monte Carlo
                
                    call filllinarray(dindx,ddata,lfirst,dddata &
                        ,mdata,n,j1,j2,lag,k,nperyear,imens,indxmx &
                        ,indx,data,nperyear,yrbg,yred,nensmx &
                        ,-999,-999,yrmo)
                    if ( month == 0 ) then
                        df = (n-ndup(0))/(max(lsum,lsum2) + decor)/real(max(1,1-ndiff)) - 2
                    else
                        df = (n-ndup(0))/(1+(max(lsum,lsum2)-1) &
                            /nperyear+ decor/nperyear)/real(max(1,1-ndiff))- 2
                    endif
                    call pearsncross(ddata,dindx,n,r,prob,z,adata &
                        ,sxx,aindx,syy,sxy,df,ncrossvalidate)
                    call fitcross(dindx,ddata,n,sig,0,a(2),a(1) &
                        ,da(2,2),da(1,1),chi2,q,ncrossvalidate,aa &
                        ,bb,.false.)
                    if ( lag /= 0 ) then
                        write(0,*) 'correlate: error: cannot handle running lagged correlations yet'
                        write(*,*) 'correlate: error: cannot handle running lagged correlations yet'
                        call exit(-1)
                    endif
                    call getred(alpha,j1,j2,lag,k,nperyear,imens &
                        ,indxmx,indx,data,nperyear,yrbg,yred,nensmx,a(2),a(1))
                    if ( lwrite ) then
                        print *,'correlate: b,a,sd,r,sdeff,alpha = ' &
                            ,a,sqrt(sxx/(n-1)),r,sqrt(sxx/(n-1)*(1-r**2)),alpha
                    endif
                    allocate(mcdata(nperyear,yrbg:yred,0:nensmx))
                    allocate(mcindx(nperyear,yrbg:yred,0:nensmx,indxmx))
                    do i=1,nmc
                        call makemcseries(mcdata,mcindx,nperyear &
                            ,yrbg,yred,nensmx,indxmx,nperyear,k,lag &
                            ,j1,j2,imens,adata,sxx,aindx,syy,sxy &
                            ,alpha,n)
                        call filllinarray(dindx,ddata,lfirst,dddata &
                            ,mdata,n,j1,j2,lag,k,nperyear,imens &
                            ,indxmx,mcindx,mcdata,nperyear,yrbg &
                            ,yred,nensmx,-999,-999,yrmo)
                        if ( lwrite ) then
                            lwrite = .false. 
                            call printcorr(dindx,ddata,lfirst,dddata &
                                ,yrmo,n,ndup(0),j1,j2,month &
                                ,nperyear,lag,strindx(k),.false., &
                                .true.,results(i),dresult,prob)
                            lwrite = .true. 
                        endif
                        if ( .false. .and. i <= 10 ) then
                            lprint = .true. 
                            write(14,'(a)')
                            write(14,'(a)')
                        else
                            lprint = .false. 
                        endif
                        call getruncorr(dindx,ddata,lfirst,dddata &
                            ,mdata,j1,j2,lag,k,month,nperyear,imens &
                            ,indxmx,mcindx,mcdata,nperyear,yrbg &
                            ,yred,nensmx,ndup(0),strindx(k), &
                            .false.,lprint,rmins(i),rmaxs(i),zdifs(i))
                        if ( lwrite ) print *,'zdif = ',zdifs(i)
                        call keepalive(i,nmc)
                    enddo
                    if ( irunvar > 0 ) then
                        if ( lwrite ) call getsign('result ',result,results,nmc,1,sign,.true.)
                        if ( lweb ) then
                            print '(a,i6,a)','Significances are computed against a ',nmc,' sample Monte Carlo<br>'
                        else
                            print '(a,i6,a)','# significances are computed against a ',nmc,' sample Monte Carlo'
                        endif
                        call getsign(runs(1,irunvar),rmin,rmins,nmc,-1,signmin, .true. )
                        call getsign(runs(2,irunvar),rmax,rmaxs,nmc,1,signmax, .true. )
                        call getsign(runs(3,irunvar),zdif,zdifs,nmc,1,signdif, .true. )
                        write(14,'(3a,f6.2,a,e14.6,a)') '# ',runs(1,irunvar),' = ',rmin,' P ',1-signmin
                        write(14,'(3a,f6.2,a,e14.6,a)') '# ',runs(2,irunvar),' = ',rmax,' P ',1-signmax
                        write(14,'(3a,f6.2,a,e14.6,a)') '# ',runs(3,irunvar),' = ',zdif,' P ',1-signdif
                    end if
                endif       ! running correlations
            800 continue
            enddo           ! k (index)
        
!           new type dumpfile
        
            if ( dump ) then
                n = 0
                do k=1,indxuse
                    if ( lincl(k) ) n = n+1
                enddo
                if ( nfittime > 0 ) n = n+1
                write(line,'(a,i2,a)') '(',n+1,'g14.6,i5,i4,i5,i4)'
                if ( lrank ) then
!                   replace values by their quantile in the PDF
                    init = 0
                    do k=1,indxuse
                        if ( lincl(k) ) then
                            call filllinarray(dindx,ddata,lfirst,dddata,mdata,n &
                                ,j1,j2,lag,k,nperyear,imens,indxmx,indx,data &
                                ,nperyear,yrbg,yred,nensmx,yrstart,yrstop,yrmo)
                            if ( init == 0 ) then
                                ! data
                                call ffsort(ddata,iarray,n)
                                do i=1,n
                                    yr = yrmo(1,iarray(i))
                                    mo = yrmo(2,iarray(i))
                                    iens = yrmo(3,iarray(i))
                                    data(mo,yr,iens) = i
                                end do
                                init = 1
                            end if
                            ! index
                            call ffsort(dindx,iarray,n)
                            do i=1,n
                                yr = yrmo(1,iarray(i))
                                mo = yrmo(2,iarray(i))
                                m = mo - lag
                                call normon(m,yr,ii,nperyear)
                                if ( imens(k) > 0 ) then
                                    iens = yrmo(3,i)
                                    indx(m,ii,iens,k) = i
                                else
                                    indx(m,ii,0,k) = i
                                end if
                            end do
                        end if
                    end do
                end if
                do iens=nens1,nens2
                    if ( nens2 > nens1 ) then
                        write(10,'(a)')
                        write(10,'(a,i4)')'# ensemble member ',iens
                    endif
                    do yr=yr1-1,yr2
                        do jj=j1,j2
                            if ( fix2 ) then
                                j = jj+lag
                            else
                                j = jj
                            endif
                            call normon(j,yr,i,nperyear)
                            if ( i < yr1 .or. i > yr2 ) cycle
                            m = j-lag
                            call normon(m,i,ii,nperyear)
                            if ( ii < yrbg .or. ii > yred ) cycle
                            n = 0
                            do k=1,indxuse
                                if ( lincl(k) ) then
                                    n = n + 1
                                    if ( imens(k) > 0 ) then
                                        ss(n) = indx(m,ii,iens,k)
                                    else
                                        ss(n) = indx(m,ii,0,k)
                                    endif
                                endif
                            enddo
                            n = n + 1
                            ss(n) = data(j,i,iens)
                            if ( nfittime > 0 ) then
                                n = n + 1
                                ss(n) = indx(j,i,iens,indxmx)
                            endif
                            l = 0
                            do k=1,n-1
                                if ( ss(k) > 0.9*absent .or. &
                                    .not. lconting .and. .not. ( &
                                    (ss(k) <= maxindx) .eqv. &
                                    (ss(k) >= minindx) .eqv. &
                                    (maxindx >= minindx) ) ) then
                                    ss(k) = -999.9
                                    l = l + 1
                                endif
                            enddo
                            if ( l < n-1 .and. &
                                ss(n) < 0.9*absent .and. &
                                ( lconting .or. ( &
                                (ss(n) <= maxdata) .eqv. &
                                (ss(n) >= mindata) .eqv. &
                                (maxdata >= mindata) ) ) ) then
                                write(10,line) (ss(k),k=1,n),i,j,ii,m
                            endif
                        enddo ! month
                    enddo   ! year
                enddo       ! ensemble
            endif           ! dump
        enddo               ! lag
    enddo                   ! month
    call printcorrfooter
    call savestartstop(yrstart,yrstop)

!   error messages

    goto 999
902 print *,'error reading NINO data'
    call exit(-1)
904 print *,'error reading sunspot data'
    print *,line
    call exit(-1)
905 print *,'error reading lag value ',line
    call exit(-1)
906 print *,'error reading sum value ',line
    call exit(-1)
907 print *,'error reading begin value ',line
    call exit(-1)
908 print *,'error reading end value ',line
    call exit(-1)
909 print *,'error reading maximum index value ',line
    call exit(-1)
910 print *,'error reading minimum index value ',line
    call exit(-1)
911 print *,'error reading minimum ks cut-off or "plot" ',line
    call exit(-1)
912 print *,'error reading month ',line
    call exit(-1)
913 print *,'error reading ldecorrelation length ',line
    call exit(-1)
914 print *,'error reading number of months to select ',line
    call exit(-1)
915 print *,'error: both summing and selecting does not make sense'
    call exit(-1)
916 print *,'error: cannot interpolate in time (yet) ',n,nperyear
    call exit(-1)
999 continue
end program

subroutine checkdup(lboot,ndup,nperyear,j1,j2)
    implicit none
    logical :: lboot
    integer :: nperyear,ndup(0:nperyear),j1,j2
    integer :: j
    ndup(0) = 0
    do j=j1,j2
        ndup(0) = ndup(0) + ndup(j)
    enddo
    if ( ndup(0) > 0 ) then
        lboot = .false. 
        print '(a,i10,a)','# cannot compute bootstrap with ',ndup(0),' duplicates'
    endif
end subroutine checkdup
