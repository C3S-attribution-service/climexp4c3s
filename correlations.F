        program correlate
*  #[ comment:
*
*       correlate station parameters from www.ncdc.noaa.gov or other
*       sources with all kinds of indices.  Currently included are:
*           SOI SLP index, from 1866 to 1997, from Phil Jones via Marleen
*           4 NINO SST indices, from 1950 to present, from
*               http://nic.fb4.noaa.gov:80/data/cddb/
*           NAO SLP index, from 1865 to 1996
*           SIDC sunspot index, from Zuerich,  http://www.astroinfo.ch/sunspot
*           SIDC sunspot cycle length, from Zuerich
*           time (in years)
*           any other file in .dat or .txt format (max 10)
*
*       GJvO hack, nov-1997, revised dec-1997, jan-1998, dec 1998, oct 1999,
*       jan-2000, added ensembles sep-2002, added running correlatons
*       2004, finally implemented ensembles correctly 2005
*
*  #] comment:
*  #[ declarations:
        implicit none
        integer nmc,nvarmax
        parameter(nmc=1000,nvarmax=1)
#include "param.inc"
#include "getopts.inc"
        integer i,ii,j,jj,k,l,m,n,j1,j2,if,im,ip,jm,jp,lag,year,month
     +       ,yr,n1,ks,ks1,ks2,nperyear,ldir,iens,jens,ilen,jindx
     +       ,ndup(0:npermax),ndum(0:npermax),validens(nensmax),
     +       imens1(0:indxmx),imens(0:indxmx),yr1s,yr2s,nunequal
     +       ,iunequal,mo,nfac(indxmx),yrstart,yrstop,mdata,nx,ny,nz,nt
     +       ,nvars,nu,yrbg,yred,nensmx
        integer,allocatable :: yrmo(:,:)
        real,allocatable :: data(:,:,:),indx(:,:,:,:),
     +       mcdata(:,:,:),mcindx(:,:,:,:),
     +       ddata(:),dddata(:),dindx(:),aa(:),bb(:),
     +       aa1(:,:,:),bb1(:,:,:)
        real anino(8),val12(npermax),adata,aindx,sxx,syy,yrmin(25),
     +       yrmax(25),xmin,xmax,rise,fall,slength(25),ayr,r,prob,
     +       absent,z,probd,sxy,df,
     +       ss(nensmax*(indxmx+1)),filter(100),minval,
     +       result,dresult(-2:2),a(2),da(2,2),results(nmc),sig(1),
     +       chi2,q,rmin,rmins(nmc),rmax,rmaxs(nmc),zdif,zdifs(nmc),
     +       addfac(npermax,indxmx),sign,signmin,signmax,signdif,
     +       alpha
	logical lfirstzero,lexist,lboot,ensindex,lprint,xwrap,
     +       laddfile(indxmx),lbb1allocated
        logical,allocatable :: lfirst(:)
        parameter (absent=3e33)
        character line*256,string*10,dir*256,ensfile*256,var*40,units*60
     +       ,ivar*40,iunits*20,newunits*60,dum1*40,dum2*40,varorg*40
     +       ,unitsorg*60
        character*4 runs(3,2)
        integer iargc,llen
        data runs /'rmin','rmax','zdif', 'bmin','bmax','bdif'/
*  #] declarations:
*  #[ check arguments:
*
*       check arguments
*       
        lwrite = .FALSE.
        n = iargc()
        if ( n.lt.1 ) then
            print *,'usage: correlate datafile index [lag n[:m]] '//
     +            '[month m[:n] [sum|ave|max|min|sel m] [log|sqrt|rank]'
     +            //' [begin yr] [end yr] [ks plot|cut-off] [detrend]',
     +            '[diff [nyr]] [plot file] [lt maxindex] [gt minindex]'
     +            //'[decor n] [runcorr|runregr nyr outfile [random '//
     +           'series|index] [noise white|red]]'
            stop
        endif
        if ( indxmx.ne.19 ) then
            write(0,*) 'expecting indxmx=19, not ',indxmx
            call exit(-1)
        endif
        if ( lwrite ) print *,'putting arrays to absent'
        if ( nfittime.gt.0 ) then
*           the last value is a heuristic based on advise in NumRec
            call savgol(filter,100,nfittime,nfittime,1,min(nfittime,4))
        endif
*  #] check arguments:
*  #[ allocate big arrays:
*
*       do not make them much bigger then needed...
*
        call getarg(1,line)
        call getfileunits(line,nx,ny,nz,nt,nperyear,nvarmax,nvars,var
     +       ,units,newunits,xwrap,lwrite)
        if ( nperyear.le.0 ) then
            write(0,*) 'correlate: error: the file ',trim(line),
     +           ' is not readable'
            call exit(-1)
        end if
!       do not let the arrays get too large
        if ( nperyear.le.12 ) then
            yrbg = yrbeg
            nensmx = nensmax
       else if ( nperyear.le.366 ) then
            yrbg = 1700
            nensmx = 10
        else
            yrbg = 1900
            nensmx = 1
        end if
        yred = yrend
        mdata = nperyear*(yred-yrbg+1)*(nensmx+1)
        if ( lwrite ) print *,'allocating ',mdata,' reals: ',
     +       nperyear,yrbg,yred,nensmx
        allocate(data(nperyear,yrbg:yred,0:nensmx))
        if ( lwrite ) print *,'allocating ',mdata*indxmx,' reals: ',
     +       nperyear,yrbg,yred,nensmx,indxmx
        allocate(indx(nperyear,yrbg:yred,0:nensmx,indxmx))
        allocate(ddata(mdata))
        allocate(dddata(mdata))
        allocate(dindx(mdata))
        allocate(lfirst(mdata))
        if ( lwrite ) print *,'allocating ',mdata*2,' reals'
        allocate(yrmo(2,mdata))
*
*  #] allocate big arrays:
*  #[ init:
*
*       init
*       
        imens1 = 0
        imens = 0
        yrstart = yred
        yrstop  = yrbg
        call getenv('DIR',dir)
        ldir = llen(dir)
        if ( ldir.le.1 ) then
            call getenv('HOME',dir)
            ldir = llen(dir)
            dir = dir(1:llen(dir))//'/climexp/'
            ldir = llen(dir)
        elseif ( dir(ldir:ldir).ne.'/' ) then
            ldir = ldir + 1
            dir(ldir:ldir) = '/'
        endif
        iunequal = 0
        nunequal = 1
*  #] init:
*  #[ read data:
*
*       read data from station file downloaded from
*       http://www.ncdc.noaa.gov/ghcn/ghcnV1.CLIMVIS.html
*
        call getarg(1,line)
        if ( lwrite ) print *,'reading data file ',trim(line)
        call readensseries(line,data,nperyear,yrbg,yred,nensmx
     +       ,n,imens1(0),imens(0),var,units,lstandardunits
     +       ,lwrite)
        if ( n.ne.nperyear ) call exit(-1)
        if ( imens(0).gt.0 ) ensemble = .TRUE.
*
*  #] read data:
*  #[ process options:
*
*       process options
*
        n = iargc()
        call getopts(2,n,nperyear,yrbg,yred,.true.,imens1(0),imens(0))
        if ( lstandardunits ) then
            ! not right first time...
            if ( lwrite ) print *,'calling makestandardseries'
            varorg = var
            unitsorg = units
            do iens = nens1,nens2
                var = varorg
                units = unitsorg
                call makestandardseries(data(1,yrbg,iens),nperyear,
     +               yrbg,yred,nperyear,var,units,lwrite)
            end do
        end if
        if ( minfac.lt.0 ) minfac = 0.75
        if ( imens(0).gt.0 ) then
            print '(a,i4,a,i4,a)','# taking ensemble members ',nens1
     +            ,' to ',nens2,' of series'
            ensemble = .TRUE.
        endif
	if ( lsubtract ) write(0,*) 'correlate: subtract not yet '//
     +        'implemented'
        do i=1,indxmx
            if ( lincl(i) ) goto 100
        enddo
!       demand at least one series...
        write(0,*) 'correlate: please select an index to correlate with'
        call exit(-1)
  100   continue
***        if ( lsum.gt.1 .and. lsel.gt.1 ) goto 915
        if ( dump ) then
            write(10,'(2a)') '# ',line(1:index(line,' '))
            write(10,'(6a)') '# ',trim(var),' [',trim(units),']'
            if ( logscale ) write(10,'(a)') '# logarithmic plot'
            if ( sqrtscale ) write(10,'(a)') '# sqrt plot'
        endif
        if ( ncrossvalidate.gt.0 ) then
            allocate(aa(mdata))
            allocate(bb(mdata))
        end if
*  #] process options:
*  #[ get SOI,NINO,NAO:
*
        do k=1,indxmx
            call makeabsent(indx(1,yrbg,0,k),nperyear,yrbg,yred)
            imens1(k) = 0
            imens(k) = 0
        enddo
*
*       get the SOI data
*
        if ( lincl(1) ) then
            if ( lwrite ) print *,'reading SOI file soi.dat'
            call readseries(trim(dir)//'CRUData/soi.dat',
     +           indx(1,yrbg,0,1),nperyear,yrbg,yred,n,var,units
     +           ,.false.,lwrite)
            if ( n.ne.nperyear ) goto 916
        endif
*
*       get NINO indices from ERSST v4
*
        if ( lincl(2) ) then
            if ( lwrite ) print *,'reading NCDCData/ersst_nino12a.dat'
            call readseries(trim(dir)//'NCDCData/ersst_nino12.dat'
     +           ,indx(1,yrbg,0,2),nperyear,yrbg,yred,n,var,units
     +           ,.false.,lwrite)
            if ( n.ne.nperyear ) goto 916
        end if
        if ( lincl(3) ) then
            if ( lwrite ) print *,'reading NCDCData/ersst_nino3a.dat'
            call readseries(trim(dir)//'NCDCData/ersst_nino3a.dat'
     +           ,indx(1,yrbg,0,3),nperyear,yrbg,yred,n,var,units
     +           ,.false.,lwrite)
            if ( n.ne.nperyear ) goto 916
        end if
        if ( lincl(4) ) then
            if ( lwrite ) print *,'reading NCDCData/ersst_nino4a.dat'
            call readseries(trim(dir)//'NCDCData/ersst_nino4a.dat'
     +           ,indx(1,yrbg,0,4),nperyear,yrbg,yred,n,var,units
     +           ,.false.,lwrite)
            if ( n.ne.nperyear ) goto 916
        end if
        if ( lincl(5) ) then
            if ( lwrite ) print *,'reading NCDCData/ersst_nino3.4a.dat'
            call readseries(trim(dir)//'NCDCData/ersst_nino3.4a.dat'
     +           ,indx(1,yrbg,0,5),nperyear,yrbg,yred,n,var,units
     +           ,.false.,lwrite)
            if ( n.ne.nperyear ) goto 916
        end if
*
*       and the NAO data
*
        if ( lincl(6) ) then
            call readseries(trim(dir)//'CRUData/nao.dat',
     +           indx(1,yrbg,0,6),nperyear,yrbg,yred,n,var,units
     +           ,.false.,lwrite)
            if ( nperyear.ne.n ) goto 916
        endif
*  #] get SOI,NINO,NAO:
*  #[ get sunspot data:
*
*       and the sunspot data, 1749-1991
*
        if ( lincl(7) ) then
            call readseries(trim(dir)//'SIDCData/sunspots.dat',
     +           indx(1,yrbg,0,7),nperyear,yrbg,yred,n,var,units
     +           ,.false.,lwrite)
            if ( nperyear.ne.n ) goto 916
        endif
        if ( lincl(8) ) then
            n = 12
            if ( nperyear.ne.12 ) goto 916
            if ( lwrite ) print *,'Opening file sunhist_2.html'
            open(1,file=trim(dir)//'SIDCData/sunhist_2.html',status
     +            ='old')
*           data is in the format
*           no yrmin valmin yrmax valmax rise fall length
            do i=1,10000
                read(1,'(a)',err=904,end=500) line
                if ( line(1:2).eq.'  ' .and.
     +                (line(3:3).ge.'0' .and. line(3:3).le.'9' .or.
     +                 line(3:3).eq.' ') .and.
     +                (line(4:4).ge.'0' .and. line(4:4).le.'9') ) then
                    read(line,*,err=904) j
                    if ( j.lt.22 ) then
                        read(line,*,err=904) j,yrmin(j),xmin,yrmax(j)
     +                        ,xmax,rise,fall,slength(j)
                    else
                        read(line,*,err=904) j,yrmin(j),xmin,yrmax(j)
     +                        ,xmax
                    endif
                    if ( lwrite ) print *,'read ',j,yrmin(j),xmin
     +                    ,yrmax(j),xmax,rise,fall,slength(j)
                endif
            enddo
  500       continue
            j = 1
            do year=yrbg,yred
                do month=1,12
                    ayr = year + (month-0.5)/12
                    if ( j.le.22 ) then
                        if ( ayr.gt.yrmin(j) ) then
                            j = j + 1
                        endif
                        if ( j.gt.1 .and. j.le.22 ) then
                            indx(month,year,0,8) = slength(j-1)
                            if ( lwrite ) print *,year,month,slength(j-1
     +                            )
                        endif
                    endif
                enddo
            enddo
            close(1)
        endif
        if ( lincl(9) ) then
            do i=yrbg,yred
                do j=1,nperyear
                    indx(j,i,0,9) = i-2000 + (j-0.5)/nperyear
                enddo
            enddo
        endif
*  #] get sunspot data:
*  #[ get user file:
*
*       get my own file
*
        do i=10,indxuse
            if ( lwrite ) print *,'Opening index file '
     +           ,trim(indexfiles(i))
            if ( indexfiles(i).eq.' ' ) then
                write(0,*) 'correlate: error: empty file name ',i
                call exit(-1)
            end if
            call readensseries(indexfiles(i),indx(1,yrbg,0,i),nperyear
     +           ,yrbg,yred,nensmx,n,imens1(i),imens(i),ivar,iunits
     +           ,lstandardunits,lwrite)
            if ( lwrite ) print *,'nperyear from this file = ',n
            if ( n.ne.nperyear ) goto 916
            if ( imens(i).gt.0 ) then
                ensindex = .TRUE.
                if ( .not.ensemble ) then
                    nens1 = max(nens1,imens1(i))
                    if ( nens2.gt.0 ) then
                        nens2 = min(nens2,imens(i))
                    else
                        nens2 = imens(i)
                    endif
                    write(*,'(a,i2,a,i3,a,i3)') '# Found ensemble ',i
     +                   ,' from ',imens1(i),' to ',imens(i)
                else
                    if ( imens(i).ne.imens(0) ) then
                        iunequal = iunequal + 1
                        if ( iunequal.ge.nunequal ) then
                            nunequal = 2*nunequal
                            write(0,'(a,2i3,a)')
     +                           'warning: unequal sizes ensembles: ',
     +                           imens(0),imens(i),' using smallest'
                            write(*,'(a,2i3,a)')
     +                           '# warning: unequal sizes ensembles: ',
     +                           imens(0),imens(i),' using smallest'
                        end if
                        nens2 = min(imens(i),imens(0))
                        imens(0) = nens2
                        imens(i) = nens2
                    endif
                    imens(i) = min(nens2,imens(i))
                    do k=1,i-1
                        if ( lincl(k) .and. 
     +                       imens(k).gt.0 .and. 
     +                       imens(i).ne.imens(k) ) then
                            iunequal = iunequal + 1
                            if ( iunequal.ge.nunequal ) then
                                nunequal = 2*nunequal
                                write(0,'(a,2i3,a)')
     +                               'warning: unequal ensembles: '
     +                               ,imens(k),imens(i)
     +                               ,' using smallest'
                                write(*,'(a,2i3,a)')
     +                               '# warning: unequal ensembles: '
     +                               ,imens(k),imens(i)
     +                               ,' using smallest'
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
*  #] get user file:
*  #[ addseries:
*
*       am I being called as addseries?
*
        call getarg(0,line)
        if ( index(line,'addseries ').ne.0 ) then
            if ( lwrite ) print *,'adding series'
            if ( nens2.gt.0 ) then
                write(0,*) 'error: cannot add ensembles yet'
                call exit(-1)
            endif
            if ( lag1.ne.lag2 ) then
                write(0,*) 'error: cannot use range of lags',lag1,lag2
                call exit(-1)
            endif
            lag = lag1
*
*           construct name
*
            i=1
            do k=1,indxuse
                if ( lincl(k) ) then
                    plotfile(i:) = strindx(k)
                    j = index(strindx(k),'.')
                    if ( j.gt.0 ) then
                        i = i+j-1
                    else
                        i = i + llen(strindx(k))
                    endif
                    plotfile(i:) = '+'
                    i = i+1
                endif
            enddo
            i = i-1
            plotfile(i:i) = ' '
            print '(a)',plotfile(1:i-1)
*
*           construct file name
*
            call getarg(1,plotfile)
            do i=len(plotfile),1,-1
                if ( plotfile(i:i).eq.'/' ) goto 510
            enddo
  510       continue
            i = i+1
            plotfile(i:) = 'iadded'
            i = i+6
            ii = 0
  520       continue
            ii = ii + 1
            if ( ii.lt.10 ) then
                write(plotfile(i:),'(i1,a)') ii,'.dat'
            elseif ( ii.lt.100 ) then
                write(plotfile(i:),'(i2,a)') ii,'.dat'
            elseif ( ii.lt.1000 ) then
                write(plotfile(i:),'(i3,a)') ii,'.dat'
            else
                write(0,*) 'addseries: error: cannot open output file '
     +               ,plotfile(1:llen(plotfile))
                write(*,*) 'addseries: error: cannot open output file '
     +               ,plotfile(1:llen(plotfile))
                call exit(-1)
            endif
            open(99,file=plotfile,status='new',err=520)
            print '(a)',plotfile(1:llen(plotfile))
*
*           get the coefficients from the environment
*           I assume they are stored as FORM_a1,FORM_a2,...
*
            jindx = 0
            laddfile = .false.
            nfac = 1
            do k=1,indxuse
                if ( lincl(k) ) then
                    !!!print *,k,strindx(k),jindx
                    if ( k.eq.10 ) then
                        ! assume the first user-defined file is the main one...
                        call copyheader(indexfiles(k),99)
                    end if
                    jindx = jindx + 1
                    if ( jindx.lt.10 ) then
                        write(string,'(a,i1)') 'FORM_a',jindx
                        i = 7
                    else
                        write(string,'(a,i2)') 'FORM_a',jindx
                        i = 8
                    endif
                    call getenv(string(1:i),line)
                    if ( line.eq.' ' ) then
                        write(0,*) 'addseries: error: environment '//
     +                       'variable ',string(1:i),
     +                       ' corresponding to ',trim(strindx(k)),
     +                       ' not set'
                        call exit(-1)
                    end if
                    if ( lwrite ) print *,'found env ',string(1:i)
     +                    ,' with content ',line
                    if ( line(1:4).ne.'file' ) then
                        ! numerical values
                        nfac(k) = 1
                        do i=1,len(line)-1
                            if ( line(i:i).eq.'+' ) line(i:i) = ' '
                            if ( line(i:i+1).eq.'--' )
     +                           line(i:i+1) = '  '
                            if ( line(i:i).eq.':' ) then
                                line(i:i) = ' '
                                nfac(k) = nfac(k) + 1
                            endif
                        enddo
                        if ( nfac(k).eq.nperyear+1 ) nfac(k) = nperyear
                        if ( nfac(k).ne.1 .and. nfac(k).ne.nperyear )
     +                       then
                            write(0,*) 'addseries: error: can only use '
     +                           //'#addfac 1 or ',nperyear,', not '
     +                           ,nfac(k)
                            write(*,*) 'addseries: error: can only use '
     +                           //'#addfac 1 or ',nperyear,', not '
     +                           ,nfac(k)
                            call exit(-1)
                        endif
                        read(line,*,err=600,end=600)
     +                       (addfac(j,k),j=1,nfac(k))
                        if ( k.eq.10 .and. oper.eq.'+' .and. lsum.gt.1 )
     +                       then
                            do j=1,nfac(k)
                                addfac(j,k) = addfac(j,k)*lsum
                            enddo
                        endif
                        if ( lwrite ) then
                            print *,'adding ',strindx(k),' with weight '
     +                           ,(addfac(j,k),j=1,nfac(k))
                            print *,'using lag ',lag,fix2
                        endif
                        write(99,'(3a,366f10.6)') '# added '
     +                       ,strindx(k),' with weight',(addfac(j,k),j
     +                       =1,nfac(k))
                    else
!                       the coefficients are in a file
                        if ( .not.allocated(bb1) ) then
                            allocate(aa1(nperyear,yrbg:yred,indxmx))
                            allocate(bb1(nperyear,yrbg:yred,indxmx))
                            aa1 = 3e33
                            bb1 = 3e33
                            ilen = index(bbfile,'.dat') - 1
                            if ( ilen.le.0 ) ilen = len_trim(bbfile)
                            call readseries(bbfile(:ilen)//'_a.dat',
     +                           aa1(1,yrbg,k),nperyear,yrbg,yred,i,
     +                           dum1,dum2,.false.,lwrite)
                            if ( i.ne.nperyear ) then
                                write(0,*) 'error: found nperyear = ',i,
     +                               ' in ',trim(line),', expected '
     +                               ,nperyear
                                call exit(-1)
                            end if
                        end if
                        if ( jindx.lt.10 ) then
                            write(bbfile(ilen+1:),'(a,i1,a)') 
     +                           '_b',jindx,'.dat'
                        else
                            write(bbfile(ilen+1:),'(a,i2,a)') 
     +                           '_b',jindx,'.dat'
                        end if
                        call readseries(bbfile,bb1(1,yrbg,k),nperyear
     +                       ,yrbg,yred,i,dum1,dum2,.false.,lwrite)
                        write(99,'(4a)') '# added ',strindx(k),
     +                       ' with weights from file ',trim(bbfile)
                        laddfile(k) = .true.
!                       if there is only one value defined, set all
!                       values for the year afterwards to that value
                        n = 0
                        m = -1
                        do j=1,nperyear
                            do i=yrbg,yred
                                if ( bb1(j,i,k).lt.1e33 ) then
                                    n = n + 1
                                    m = j
                                    exit
                                end if
                            end do
                        end do
                        if ( n.eq.1 .and. nperyear.gt.1 ) then
                            do i=yrbg,yred
                                do j=1,nperyear
                                    if ( j.gt.m ) then
                                        aa1(j,i,k) = aa1(m,i,k)
                                        bb1(j,i,k) = bb1(m,i,k)
                                    elseif ( j.lt.m .and. i.gt.1 ) then
                                        aa1(j,i,k) = aa1(m,i-1,k)
                                        bb1(j,i,k) = bb1(m,i-1,k)
                                    end if
                                end do
                            end do
                            n = nperyear
                        end if
                        nfac(k) = n
                    end if
                    if ( k.ne.10 .and. lag.ne.0 ) then
                        write(99,'(a,i3,a,i3)') '# used lag ',lag
     +                       ,' for index ',k
                    endif
***                    write(99,'(a,i3)') '# ',k
  600               continue
                else
                    addfac(1,k) = 0
                endif
            enddo
            do yr=yrbg,yred
                n = 0
                do jj=1,nperyear
                    if ( fix2 ) then
                        j = jj+lag
                    else
                        j = jj
                    endif
                    call normon(j,yr,i,nperyear)
                    if ( i.lt.yrbg .or.i.gt.yred ) cycle
                    m = j-lag
                    call normon(m,i,ii,nperyear)
                    if ( ii.lt.yrbg .or.ii.gt.yred ) cycle
                    data(j,i,0) = 0
                    do k=1,indxuse
                        if ( lincl(k) .and. ( laddfile(k) .or. 
     +                       addfac(min(j,nfac(k)),k).ne.0 ) ) then
                            if ( k.eq.10 ) then ! no lag
                                if ( data(j,i,0).lt.0.9*absent .and.
     +                               indx(j,i,0,k).lt.0.9*absent) then
                                    if ( laddfile(k) ) then
                                        if ( aa1(j,i,k).lt.1e33 .and.
     +                                       bb1(j,i,k).lt.1e33 ) then
                                            print *,'data(',j,i,
     +                                           ') was ',data(j,i,0)
                                            data(j,i,0) = data(j,i,0) -
     +                                           aa1(j,i,k) - bb1(j,i,k)
     +                                           *indx(j,i,0,k)
                                            print *,'aa1,bb1(',j,i,k,
     +                                           ') = ',aa1(j,i,k),
     +                                           bb1(j,i,k)
                                            print *,'data(',j,i,
     +                                           ')  is ',data(j,i,0)
                                        else
                                            data(j,i,0) = absent
                                        end if
                                    else
                                        data(j,i,0) = data(j,i,0) +
     +                                       addfac(min(j,nfac(k)),k)
     +                                       *indx(j,i,0,k)
                                    end if
                                else
                                    data(j,i,0) = absent
                                endif
                            else
                                if ( data(j,i,0).lt.0.9*absent .and.
     +                               indx(m,ii,0,k).lt.0.9*absent) then
                                    if ( laddfile(k) ) then
                                        if ( aa1(j,i,k).lt.1e33 .and.
     +                                       bb1(j,i,k).lt.1e33 ) then
                                            data(j,i,0) = data(j,i,0) -
     +                                           aa1(j,i,k) - bb1(j,i,k)
     +                                           *indx(m,ii,0,k)
                                        else
                                            data(j,i,0) = absent
                                        end if
                                    else
                                        data(j,i,0) = data(j,i,0) +
     +                                       addfac(min(j,nfac(k)),k)
     +                                       *indx(m,ii,0,k)
                                    end if
                                else
                                    data(j,i,0) = absent
                                endif
                            endif
                        endif
                    enddo
                    if ( data(j,i,0).lt.0.9*absent ) then
                        n = n + 1
                    else
                        data(j,i,0) = -999.9
                    endif
                enddo
                if ( n.gt.0 ) then
                    write(99,'(i5,366f12.4)') i,(data(j,i,0),j=1
     +                    ,nperyear)
                endif
            enddo
            goto 999
        endif
*  #] addseries:
*  #[ manipulate time series:
*
*       monthly diffs
*
        if ( mdiff.gt.0 ) then
            if ( lwrite ) print '(a)'
     +           ,'# Taking monthly anomalies of series'
            do iens=imens1(0),imens(0)
                call mdiffit(data(1,yrbg,iens),nperyear,nperyear,yrbg
     +                ,yred,mdiff)
            enddo
        endif
        if ( mdiff2.gt.0 ) then
            do k=1,indxuse
                if ( lincl(k) ) then
                    if ( lwrite ) print '(a,i3)',
     +                    '# Taking monthly anomalies of index ',k
                    do iens=imens1(k),imens(k)
                        call mdiffit(indx(1,yrbg,iens,k),nperyear
     +                        ,nperyear,yrbg,yred,mdiff2)
                    enddo
                endif
            enddo
        endif
*
*       sum
*
        if ( lsum.gt.1 ) then
            if ( lwrite ) print '(a,i3)','# Summing series ',lsum
            do iens=imens1(0),imens(0)
                call sumit(data(1,yrbg,iens),nperyear,nperyear,yrbg
     +                ,yred,lsum,oper)
            enddo
        endif
        if ( lsum2.gt.1 ) then
            do k=1,indxuse
                if ( lincl(k) ) then
                    if ( lwrite ) print '(a,2i3)','# Averaging index ',k
     +                    ,lsum2
                    do iens=imens1(k),imens(k)
                        call sumit(indx(1,yrbg,iens,k),nperyear
     +                       ,nperyear,yrbg,yred,lsum2,'v')
                    enddo
                endif
            enddo
        endif
*
*       logscale
*
        if ( logscale ) then
            do iens=imens1(0),imens(0)
                if ( lwrite ) print '(a,2i3)','# Taking log of series '
                call takelog(data(1,yrbg,iens),nperyear,nperyear,yrbg
     +                ,yred)
            enddo
        endif
        if ( logfield ) then
            do k=1,indxuse
                if ( lincl(k) ) then
                    if ( lwrite ) print '(a,2i3)'
     +                   ,'# Taking log of index ',k
                    do iens=imens1(k),imens(k)
                        call takelog(indx(1,yrbg,iens,k),nperyear
     +                       ,nperyear,yrbg,yred)
                    enddo
                endif
            enddo
        endif
*
*       sqrtscale
*
        if ( sqrtscale ) then
            do iens=imens1(0),imens(0)
                if ( lwrite ) print '(a,2i3)','# Taking sqrt of series '
                call takesqrt(data(1,yrbg,iens),nperyear,nperyear,yrbg
     +                ,yred)
            enddo
        endif
        if ( sqrtfield ) then
            do k=1,indxuse
                if ( lincl(k) ) then
                    if ( lwrite ) print '(a,2i3)'
     +                   ,'# Taking sqrt of index ',k
                    do iens=imens1(k),imens(k)
                        call takesqrt(indx(1,yrbg,iens,k),nperyear
     +                       ,nperyear,yrbg,yred)
                    enddo
                endif
            enddo
        endif
*
*       find begin,end of valid data
*
        if ( lwrite ) print *,'adjusting yr1,yr2 from ',yr1,yr2
        do yr=yr1,yr2
            do iens=imens1(0),imens(0)
                do mo=1,nperyear
                    if ( data(mo,yr,iens).lt.1e33 ) then
                        goto 710
                    endif
                enddo
            enddo
        enddo
 710    continue
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
                            if ( indx(mo,yr,iens,k).lt.1e33 ) then
                                goto 720
                            endif
                        enddo
                    enddo
                endif
            enddo
        enddo
 720    continue
        if ( fix2 ) then
            yr1 = max(yr1,yr)
        else
            yr1 = max(yr1,yr-max(0,(lag2+nperyear-1)/nperyear))
        endif

        do yr=yr2,yr1,-1
            do iens=imens1(0),imens(0)
                do mo=1,nperyear
                    if ( data(mo,yr,iens).lt.1e33 ) then
                        goto 730
                    endif
                enddo      
            enddo
        enddo
        print '(a)','# found no valid data in data'
 730    continue
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
                            if ( indx(mo,yr,iens,k).lt.1e33 ) then
                                goto 740
                            endif
                        enddo
                    enddo
                endif
            enddo
        enddo
        print '(a,i3)','# found no valid data in index ',k
 740    continue
        if ( fix2 ) then
            yr2 = max(yr2,yr)
        else
            yr2 = max(yr2,yr-min(0,(lag1-nperyear+1)/nperyear))
        endif
        if ( lwrite ) print *,'                    to ',yr1,yr2
*
*       detrending
*
        if ( ldetrend ) then
            if ( lwrite ) print *,'Detrending series'
            do iens=imens1(0),imens(0)
                call detrend(data(1,yrbg,iens),nperyear,nperyear,
     +               yrbg,yred,yr1,yr2,m1,m2,lsel)
            enddo
            do k=1,indxuse
                if ( lincl(k) .and. k.ne.9 ) then
                    if ( lwrite ) print *,'Detrending index ',k
                    do iens=imens1(k),imens(k)
                        if ( lag1.eq.0 .and. lag2.eq.0 .or. m1.eq.0 .or.
     +                        lsel.eq.12 ) then
                            call detrend(indx(1,yrbg,iens,k),nperyear,
     +                           nperyear,yrbg,yred,yr1,yr2,m1,m2,lsel
     +                           )
                        else
                            call detrend(indx(1,yrbg,iens,k),nperyear,
     +                           nperyear,yrbg,yred,yr1,yr2,1,nperyear
     +                           ,lsel)
                        endif
                    enddo
                endif
            enddo
        endif
        if ( ndiff.ne.0 ) then
            if ( lwrite ) print *,'Taking differences - series'
            do iens=imens1(0),imens(0)
                call ndiffit(data(1,yrbg,iens),nperyear,nperyear,yrbg
     +                ,yred,ndiff,minfacsum)
            enddo
            if ( ndiff.lt.0 .and. lnooverlap ) then
                call dooverlap(data,nperyear,yrbg,yred,imens1(0)
     +               ,imens(0),nperyear,ndiff)
            end if
            do k=1,indxuse
                if ( lincl(k) .and. k.ne.9 ) then
                    if ( lwrite ) print *
     +                    ,'Taking differences - index ',k
                    do iens=imens1(k),imens(k)
                        call ndiffit(indx(1,yrbg,iens,k),nperyear
     +                       ,nperyear,yrbg,yred,ndiff,minfacsum)
                    enddo
                endif
            enddo
        endif
*
*       bias corrections
*
        n = 0
        do k=1,indxuse
            if ( lincl(k) ) n = n + 1
        end do
        if ( n.eq.1 ) then
            ! no different bias corrections for different models yet
            if ( debias.eq.1 ) then
                print '(a)','# correcting for bias in mean'
                call debiasmean(indx(1,yrbg,0,indxuse),data,nperyear,
     +               nperyear,yrbg,yred,yr1,yr2,imens1(0),imens(0),var,
     +               lwrite)
            elseif ( debias.eq.2 ) then
                print '(a)','# correcting for bias in mean and variance'
                call debiasvar(indx(1,yrbg,0,indxuse),data,nperyear,
     +               nperyear,yrbg,yred,yr1,yr2,imens1(0),imens(0),var,
     +               lwrite)
            elseif ( debias.eq.3 ) then
                print '(a)','# correcting for bias in distribution'
                call debiasall(indx(1,yrbg,0,indxuse),data,nperyear,
     +               nperyear,yrbg,yred,yr1,yr2,imens1(0),imens(0),var,
     +               lwrite)
            endif
            if ( debias.gt.0 ) then
                ncrossvalidate = 0 ! is already included
            end if
        end if
*
*       anomalies - necessary if we fit time derivative with relaxation
*
        if ( anom .or. nfittime.gt.0 .and. ndiff.le.0 ) then
            if ( lwrite ) print *,'Taking anomalies - data'
            do iens=imens1(0),imens(0)
                call anomal(data(1,yrbg,iens),nperyear,nperyear,yrbg
     +                ,yred,yr1,yr2)
            enddo
            do k=1,indxuse
                if ( lincl(k) ) then
                    if ( lwrite ) print *
     +                    ,'Taking anomalies - index ',k
                    do iens=imens1(k),imens(k)
                        call anomal(indx(1,yrbg,iens,k),nperyear
     +                        ,nperyear,yrbg,yred,yr1,yr2)
                    enddo
                endif
            enddo
        endif
!
!       anomalies wrt ensemble mean
!
        if ( lensanom .and. imens1(0).ne.imens(0) ) then
            if ( lwrite ) print *
     +           ,'Taking anomalies wrt ensemble mean - data'
            call anomalensemble(data,nperyear,nperyear,yrbg,
     +           yred,yr1,yr2,max(nens1,imens1(0)),min(nens2,imens(0)))
        endif
        do k=1,indxuse
            if ( lincl(k) .and. lensanom .and. imens1(k).ne.imens(k))
     +           then
                if ( lwrite ) print *
     +               ,'Taking anomalies wrt ensemble mean - index ',k
                call anomalensemble(indx(1,yrbg,0,k),nperyear,nperyear,
     +               yrbg,yred,yr1,yr2,max(nens1,imens1(k)),
     +               min(nens2,imens(k)))
            endif
        enddo
*
*       copy ensemble members so that there is the same
*       number of valid ones at every time step
*
        do j=0,nperyear
            ndup(j) = 0
        enddo
        if ( ensemble .and. lmakeensfull ) then
	        do k=1,indxmx
    	        if ( ensemble .and. imens(k).gt.0 ) then                
        	        do yr=yrbg,yred
            	        do mo=1,nperyear
                	        do iens=nens1,nens2                            
                    	        if ( data(mo,yr,iens).lt.1e33 .neqv.
     +                  	          indx(mo,yr,iens,k).lt.1e33 ) then
                            	    if ( lwrite) print *,'correlate: ',
     +                          	      'warning: inconsistent ',
     +									  'ensembles ',mo,yr
     +                              	  ,iens,k,data(mo,yr,iens)
     +               	                  ,indx(mo,yr,iens,k)
                        	        data(mo,yr,iens) = 3e33
                            	    indx(mo,yr,iens,k) = 3e33
	                            endif
    	                    enddo
        	            enddo
            	    enddo
    	        endif
        	enddo
            if ( lwrite ) print *,
     +          'correlate: calling makeensfull for data'
            call makeensfull(ndup,nperyear,data,nperyear,yrbg,yred
     +            ,nens1,nens2,validens,lwrite)
            call random_number(prob)
        endif
        if ( ensindex .and. lmakeensfull ) then
            do k=1,indxmx
                if ( imens(k).gt.0 ) then
                    call makeensfull(ndum,nperyear,indx(1,yrbg,0,k)
     +                   ,nperyear,yrbg,yred,nens1,nens2,validens
     +                   ,lwrite)
		    call random_number(probd)
		    if ( prob.ne.probd ) then
			write(0,*) 'error: ensembles not filled '//
     +				'identically: ',prob,probd
		    endif
                endif
            enddo
        endif
*  #] manipulate time series:
*  #[ correlate:
*
*       correlate!
*
        call printcorrheader
        do month=m1,m2
            call getj1j2(j1,j2,month,nperyear,.true.)
            lboot = lbootstrap
            if ( ensemble ) then
                call checkdup(lboot,ndup,nperyear,j1,j2)
            endif
*
*           fill linear arrays without absent values and compute r
            do lag=lag1,lag2
                do k=1,indxuse
                    if ( .not.lincl(k) ) goto 800
                    call perc2cut(lag,k,j1,j2,nperyear,imens1,imens
     +                   ,indxmx,indx,data,nperyear,yrbg,yred,nensmx)
                    if ( lwrite ) print *,'index ',strindx(k),k
                    if ( dump ) write(10,'(2a)') '# ',strindx(k)
                    call filllinarray(dindx,ddata,lfirst,dddata,mdata,n
     +                   ,j1,j2,lag,k,nperyear,imens,indxmx,indx,data
     +                   ,nperyear,yrbg,yred,nensmx,filter,yrstart
     +                   ,yrstop,yrmo)
                    if ( n.lt.minnum ) then
                        if ( lwrite ) print *
     +                        ,'correlate: not enough points: ',n,minnum
                        goto 800
                    endif
                    call printcorr(dindx,ddata,lfirst,dddata,yrmo,n
     +                   ,ndup(0),j1,j2,month,nperyear,lag,strindx(k)
     +                   ,lboot,.true.,result,dresult,prob)
*
*                   running correlations
*
                    if ( nyrwindow.gt.0 ) then
                        call getruncorr(dindx,ddata,lfirst,dddata,mdata
     +                       ,j1,j2,lag,k,month,nperyear,imens,indxmx
     +                       ,indx,data,nperyear,yrbg,yred,nensmx
     +                       ,ndup(0),filter,strindx(k),lboot,.true.
     +                       ,rmin,rmax,zdif)
*
*                       Monte Carlo
*
                        call filllinarray(dindx,ddata,lfirst,dddata
     +                       ,mdata,n,j1,j2,lag,k,nperyear,imens,indxmx
     +                       ,indx,data,nperyear,yrbg,yred,nensmx
     +                       ,filter,-999,-999,yrmo)
                        if ( month.eq.0 ) then
                            df = (n-ndup(0))/(max(lsum,lsum2) + decor)/
     +                            real(max(1,1-ndiff)) - 2
                        else
                            df = (n-ndup(0))/(1+(max(lsum,lsum2)-1)
     +                            /nperyear+ decor/nperyear)/real(max(1
     +                            ,1-ndiff))- 2
                        endif
                        call pearsncross(ddata,dindx,n,r,prob,z,adata
     +                       ,sxx,aindx,syy,sxy,df,ncrossvalidate)
                        call fitcross(dindx,ddata,n,sig,0,a(2),a(1)
     +                       ,da(2,2),da(1,1),chi2,q,ncrossvalidate,aa
     +                       ,bb,.false.)
                        if ( lag.ne.0 ) then
                            write(0,*) 'correlate: error: cannot hand'//
     +                           'le running lagged correlations yet'
                            write(*,*) 'correlate: error: cannot hand'//
     +                           'le running lagged correlations yet'
                            call exit(-1)
                        endif
                        call getred(alpha,j1,j2,lag,k,nperyear,imens
     +                       ,indxmx,indx,data,nperyear,yrbg,yred
     +                       ,nensmx,a(2),a(1))
                        if ( lwrite ) then
                            print *,'correlate: b,a,sd,r,sdeff,alpha = '
     +                           ,a,sqrt(sxx/(n-1)),r,sqrt(sxx/(n-1)*(1
     +                           -r**2)),alpha
                        endif
                        allocate(mcdata(nperyear,yrbg:yred,0:nensmx))
                        allocate(mcindx(nperyear,yrbg:yred,0:nensmx,
     +                       indxmx))
                        do i=1,nmc
                            call makemcseries(mcdata,mcindx,nperyear
     +                           ,yrbg,yred,nensmx,indxmx,nperyear,k,lag
     +                           ,j1,j2,imens,adata,sxx,aindx,syy,sxy
     +                           ,alpha,n)
                            call filllinarray(dindx,ddata,lfirst,dddata
     +                           ,mdata,n,j1,j2,lag,k,nperyear,imens
     +                           ,indxmx,mcindx,mcdata,nperyear,yrbg
     +                           ,yred,nensmx,filter,-999,-999,yrmo)
                            if ( lwrite ) then
                                lwrite = .false.
                                call printcorr(dindx,ddata,lfirst,dddata
     +                               ,yrmo,n,ndup(0),j1,j2,month
     +                               ,nperyear,lag,strindx(k),.false.,
     +                               .true.,results(i),dresult,prob)
                                lwrite = .true.
                            endif
                            if ( .false. .and. i.le.10 ) then
                               lprint = .true.
                               write(14,'(a)')
                               write(14,'(a)')
                            else
                               lprint = .false.
                            endif
                            call getruncorr(dindx,ddata,lfirst,dddata
     +                           ,mdata,j1,j2,lag,k,month,nperyear,imens
     +                           ,indxmx,mcindx,mcdata,nperyear,yrbg
     +                           ,yred,nensmx,ndup(0),filter,strindx(k),
     +                           .false.,lprint,rmins(i),rmaxs(i)
     +                           ,zdifs(i))
                            if ( lwrite ) print *,'zdif = ',zdifs(i)
                            call keepalive(i,nmc)
                        enddo
                        if ( lwrite ) call getsign('result ',result
     +                        ,results,nmc,1,sign,.true.)
                        if ( lweb ) then
                            print '(2a,i6,a)'
     +                           ,'Significances are computed'
     +                           ,' against a ',nmc
     +                           ,' sample Monte Carlo<br>'
                        else
                            print '(2a,i6,a)'
     +                           ,'# significances are computed'
     +                           ,' against a ',nmc
     +                           ,' sample Monte Carlo'
                        endif
                        call getsign(runs(1,irunvar),rmin,rmins,nmc,-1
     +                       ,signmin,.true.)
                        call getsign(runs(2,irunvar),rmax,rmaxs,nmc,1
     +                       ,signmax,.true.)
                        call getsign(runs(3,irunvar),zdif,zdifs,nmc,1
     +                       ,signdif,.true.)
                        write(14,'(3a,f6.2,a,e14.6,a)') '# ',
     +                       runs(1,irunvar),' = ',rmin,' P ',1-signmin
                        write(14,'(3a,f6.2,a,e14.6,a)') '# ',
     +                       runs(2,irunvar),' = ',rmax,' P ',1-signmax
                        write(14,'(3a,f6.2,a,e14.6,a)') '# ',
     +                       runs(3,irunvar),' = ',zdif,' P ',1-signdif
                    endif       ! running correlations
  800               continue
                enddo           ! k (index)
*
*               new type dumpfile
*
                if ( dump ) then
                    n = 0
                    do k=1,indxuse
                        if ( lincl(k) ) n = n+1
                    enddo
                    if ( nfittime.gt.0 ) n = n+1
                    write(line,'(a,i2,a)') '(',n+1,'g14.6,i5,i4,i5,i4)'
                    if ( lrank ) then
!                       replace values by their quantile in the PDF
                        print '(a)'
     +                       ,'# TODO: replace values by quantiles'
                    end if
                    do iens=nens1,nens2
                        if ( nens2.gt.nens1 ) then
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
                                if ( i.lt.yr1 .or.i.gt.yr2 ) goto
     +                                810
                                m = j-lag
                                call normon(m,i,ii,nperyear)
                                if ( ii.lt.yrbg .or.ii.gt.yred ) goto
     +                                810
                                n = 0
                                do k=1,indxuse
                                    if ( lincl(k) ) then
                                        n = n + 1
                                        if ( imens(k).gt.0 ) then
                                            ss(n) = indx(m,ii,iens,k)
                                        else
                                            ss(n) = indx(m,ii,0,k)
                                        endif
                                    endif
                                enddo
                                n = n + 1
                                ss(n) = data(j,i,iens)
                                if ( nfittime.gt.0 ) then
                                    n = n + 1
                                    ss(n) = indx(j,i,iens,indxmx)
                                endif
                                l = 0
                                do k=1,n-1
                                    if ( ss(k).gt.0.9*absent .or.
     +                                    .not.lconting .and. .not.(
     +                                    (ss(k).le.maxindx) .eqv.
     +                                    (ss(k).ge.minindx) .eqv.
     +                                    (maxindx.ge.minindx) ) ) then
                                        ss(k) = -999.9
                                        l = l + 1
                                    endif
                                enddo
                                if ( l.lt.n-1 .and.
     +                                ss(n).lt.0.9*absent .and.
     +                                ( lconting .or. (
     +                                (ss(n).le.maxdata) .eqv.
     +                                (ss(n).ge.mindata) .eqv.
     +                                (maxdata.ge.mindata) ) ) ) then
                                    write(10,line) (ss(k),k=1,n),i,j,ii
     +                                    ,m
                                endif
  810                           continue
                            enddo ! month
                        enddo   ! year
                    enddo       ! ensemble
                endif           ! dump
            enddo               ! lag
        enddo                   ! month
        call printcorrfooter
        call savestartstop(yrstart,yrstop)
*  #] correlate:
*  #[ error messages:
*
*       error messages
*
        goto 999
  902   print *,'error reading NINO data'
        call exit(-1)
  904   print *,'error reading sunspot data'
        print *,line
        call exit(-1)
  905   print *,'error reading lag value ',line
        call exit(-1)
  906   print *,'error reading sum value ',line
        call exit(-1)
  907   print *,'error reading begin value ',line
        call exit(-1)
  908   print *,'error reading end value ',line
        call exit(-1)
  909   print *,'error reading maximum index value ',line
        call exit(-1)
  910   print *,'error reading minimum index value ',line
        call exit(-1)
  911   print *,'error reading minimum ks cut-off or "plot" ',line
        call exit(-1)
  912   print *,'error reading month ',line
        call exit(-1)
  913   print *,'error reading ldecorrelation length ',line
        call exit(-1)
  914   print *,'error reading number of months to select ',line
        call exit(-1)
  915   print *,'error: both summing and selecting does not make sense'
        call exit(-1)
 916    print *,'error: cannot interpolate in time (yet) ',n,nperyear
        call exit(-1)
  999   continue
*  #] error messages:
        end
*  #[ checkdup:        
        subroutine checkdup(lboot,ndup,nperyear,j1,j2)
        implicit none
        logical lboot
        integer nperyear,ndup(0:nperyear),j1,j2
        integer j
        ndup(0) = 0
        do j=j1,j2
            ndup(0) = ndup(0) + ndup(j)
        enddo
        if ( ndup(0).gt.0 ) then
            lboot = .FALSE.
            print '(a,i10,a)','# cannot compute bootstrap with '
     +            ,ndup(0),' duplicates'
        endif
        end
*  #] checkdup:        
