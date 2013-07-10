*  #[ correlatefield:
        program correlatefield
*
*       program to correlate a field series to a point series 
*       to give fields of correlation coefficients, probabilities that
*       these are significant, and the fit coefficients a, b and their
*       errors.
*       aug-2002 added support for ensembles; ensemble field x scalar
*       time series, and ensemble field x ensemble time series the
*       third option is not (yet?) implemented
*
        use lsdata
        implicit none
#include "recfac.h"
#include "params.h"
#include "netcdf.inc"
#include "getopts.inc"
        integer nvarmax,nyrmax,nlevmax,mensmax
        parameter (nvarmax=1,nyrmax=50,nlevmax=1,mensmax=1)
        integer j,n,nx,ny,nz,nt,firstyr,firstmo,nvars,iarg,
     +       ivars(2,nvarmax),ncid,jvars(6,nvmax),endian,status,
     +       nperyear,lastyr,mens,mens1
        logical lexist
        real,allocatable :: field(:,:,:,:,:,:),
     +       r(:,:,:,:),prob(:,:,:,:),
     +       a(:,:,:,:),b(:,:,:,:),da(:,:,:,:),db(:,:,:,:),
     +       a1(:,:,:,:),da1(:,:,:,:),cov(:,:,:,:),relregr(:,:,:,:),
     +       rmin(:,:,:,:),rmax(:,:,:,:),zdif(:,:,:,:),rprob(:,:,:,:),
     +       xn(:,:,:,:)
        real xx(nxmax),yy(nymax),zz(nzmax),undef,xxls(nxmax),yyls(nymax)
        character title*512,vars(nvarmax)*60,lvars(nvarmax)*128,
     +        units(nvarmax)*60,lsmasktype*4
        character infile*255,datfile*255
        integer iargc,llen
*
*       check arguments
*
        n = iargc()
        if ( lwrite ) print *,'correlatefield: called with ',n
     +        ,' arguments'
        if ( n.lt.3 ) then
            print *,'usage: correlatefield field.[ctl|nc] '//
     +           'series [month m[:n] [[ens]anom] [lag n[:m]] '//
     +           '[sum|ave|max|min|sel n] [sum2 n] '//
     +           '[log|sqrt|rank] '//
     +           '[minfac r] [minnum n] [begin yr] [end yr] '//
     +           '[lt cut] [gt cut] [diff [nyr]] [detrend] '//
     +           '[runcorr|runregr nyr dummyfile '//
     +           '[random series|field]] [noise white|red]'//
     +            'outfield'
            print *,'       ensemble input is denoted by %% in the name'
            stop
        endif
        call getarg(1,infile)
        mens1 = 0
        if ( index(infile,'%').gt.0 .or. index(infile,'++').gt.0 ) then
            ensemble = .true.
            call filloutens(infile,0)
            inquire(file=infile,exist=lexist)
            if ( .not.lexist ) then
                call getarg(1,infile)
                mens1 = 1
                call filloutens(infile,1)
            endif
        else
            ensemble = .false.
            mens1 = 0
            mens = 0
        endif
        if ( lwrite ) print *,'correlatefield: nf_opening file '
     +        ,infile(1:llen(infile))
        status = nf_open(infile,nf_nowrite,ncid)
        if ( status.ne.nf_noerr ) then
            if ( lwrite ) print *,'correlatefield: calling parsectl'
            call parsectl(infile,datfile,nxmax,nx,xx,nymax,ny,yy,nzmax
     +            ,nz,zz,nt,nperyear,firstyr,firstmo,undef,endian,title
     +            ,1,nvars,vars,ivars,lvars,units)
            nz = max(1,ivars(1,1))
            ncid = -1
            if ( ensemble ) then
                do mens=1,nensmax
                    call getarg(1,infile)
                    call filloutens(infile,mens)
                    inquire(file=infile,exist=lexist)
                    if ( .not.lexist ) goto 100
                enddo
  100           continue
                mens = mens - 1
                call getarg(1,infile)
                write(0,*) 'located ',mens-mens1+1
     +               ,' ensemble members of ',trim(infile),'<br>'
            endif
        else
            if ( lwrite ) print *,'correlatefield: calling parsenc'
     +            ,ncid
            call parsenc(infile,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax
     +            ,nz,zz,nt,nperyear,firstyr,firstmo,undef,title,1
     +            ,nvars,vars,jvars,lvars,units)
            if ( ensemble ) then
                do mens=1,nensmax
                    call getarg(1,infile)
                    call filloutens(infile,mens)
                    status = nf_open(infile,nf_nowrite,ncid)
                    if ( status.ne.nf_noerr ) goto 200
                enddo
  200           continue
                mens = mens - 1
                write(0,*) 'located ',mens-mens1+1
     +               ,' ensemble members<br>'
            endif
        endif
        iarg = 2
        lastyr = firstyr + (firstmo+nt-2)/nperyear
*       process arguments
        j=3
        call getlsmask(j,lsmasktype,nxmax,xxls,nymax,yyls,lwrite)
        if ( lsmasktype.ne.'all' ) then
            call checkgridequal(nx,ny,xx,yy,nxls,nyls,xxls,yyls)
        endif
        call getopts(j,n-1,nperyear,yrbeg,yrend,.true.,mens1,mens)
        if ( ensemble ) write(0,*) 'Using ensemble members ',nens1
     +        ,' to ',nens2,'<br>'
        call getarg(2,infile)
        if ( .not.ensemble .and. 
     +       (index(infile,'%').gt.0 .or. index(infile,'++').gt.0) )
     +       then
***            write(0,*) 'infile = ',infile
***            write(*,*) 'infile = ',infile
            ensemble = .true.
        endif
*
*       allocate arrays
*
        allocate(field(nx,ny,nz,nperyear,firstyr:lastyr,nens1:nens2))
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
        allocate(rmin(nx,ny,nz,0:nperyear))
        allocate(rmax(nx,ny,nz,0:nperyear))
        allocate(zdif(nx,ny,nz,0:nperyear))
        allocate(rprob(nx,ny,nz,0:nperyear))
        allocate(xn(nx,ny,nz,0:nperyear))
*
*       superfluous now
*
        if ( lwrite ) print *,'correlatefield: calling cfield'
        call cfield(datfile,ncid,field,lsmask,nx,xx,ny,yy,nz,zz,nt
     +       ,nperyear,firstyr,lastyr,firstmo,undef,endian,jvars,iarg,r
     +       ,prob,a,b,da,db,cov,a1,da1,relregr,rmin,rmax,zdif,rprob,xn
     +       ,vars,lvars,units,title,lsmasktype)
        end
*  #] correlatefield:
*  #[ cfield:
        subroutine cfield(datfile,ncid,field,lsmask,nx,xx,ny,yy,nz,zz,nt
     +       ,nperyear,firstyr,lastyr,firstmo,undef,endian,jvars,iarg,r
     +       ,prob,a,b,da,db,cov,a1,da1,relregr,rmin,rmax,zdif,rprob,xn
     +       ,invars,inlvars,inunits,intitle,lsmasktype)
*
*       break to use the array field compactly and conserve RAM
*
        implicit none
#include "recfac.h"
#include "params.h"
#include "getopts.inc"
#include "netcdf.inc"
        integer nvarmax,mpermax
        parameter (nvarmax=13,mpermax=366)
        integer nmc
        parameter(nmc=1000)
*
        integer ncid,nx,ny,nz,nt,nperyear,firstyr,lastyr,firstmo,iarg,
     +        jvars(6,nvmax),endian
        real field(nx,ny,nz,nperyear,firstyr:lastyr,nens1:nens2),
     +       lsmask(nx,ny),xx(nx),yy(nx),zz(nz),undef,
     +       r(nx,ny,nz,0:nperyear),
     +       prob(nx,ny,nz,0:nperyear),
     +       a(nx,ny,nz,0:nperyear),
     +       b(nx,ny,nz,0:nperyear),
     +       da(nx,ny,nz,0:nperyear),
     +       db(nx,ny,nz,0:nperyear),
     +       cov(nx,ny,nz,0:nperyear),
     +       a1(nx,ny,nz,0:nperyear),
     +       da1(nx,ny,nz,0:nperyear),
     +       relregr(nx,ny,nz,0:nperyear),
     +       rmin(nx,ny,nz,0:nperyear),
     +       rmax(nx,ny,nz,0:nperyear),
     +       zdif(nx,ny,nz,0:nperyear),
     +       rprob(nx,ny,nz,0:nperyear),
     +       xn(nx,ny,nz,0:nperyear)
        character datfile*(*),invars(1)*(*),inlvars(1)*(*),intitle*(*),
     +       inunits(1)*(*),lsmasktype*4
*
        integer ntmax,ndatmax
        parameter(ntmax=1000,ndatmax=2*mpermax*(yrend-yrbeg+1)*10)
        integer i,j,jj,k,kk,n,nn,lag,jx,jy,jz,yr,month,j1,j2,m,mm,mo,ii
     +       ,l,nvars,ivars(2,nvarmax),ldir,ntvarid,itimeaxis(ntmax)
     +       ,nrec,iens,jens,ndup(0:mpermax),mens,validens(nensmax)
     +       ,nens2series,iens2,status,imens(0:1),nold,yrstart,yrstop
     +       ,fyr,yrmo(2,ndatmax),mdatmax,irec,ntp,ndiffn
        real,allocatable :: data(:,:,:),mcdata(:,:,:),fxy(:,:,:)
        real,allocatable :: aaa1(:,:,:,:,:),bbb1(:,:,:,:,:),
     +       aaa(:),bbb(:),field2(:,:,:,:,:,:)
        real ddata(ndatmax),dindx(ndatmax),dddata(ndatmax),
     +       adata,sxx,aindx,syy,sxy,df,d,zd,z,probd,absent,sig(1),chi2,
     +       q,sum,fac,filter(100),aa,daa,bb,dbb,dresult(-2:2),
     +       results(nmc),rmins(nmc),rmaxs(nmc),zdifs(nmc),dum,zold,
     +       sxxold,alpha,xrand,s
        parameter (absent=3e33)
        logical lexist,ensseries,lfirst(ndatmax),llwrite
        logical,allocatable :: lnewyr(:,:)
        character line*80,yesno*1,string*10,file*255,infile*255,
     +       outfile*255,title*250,vars(nvarmax)*60,lvars(nvarmax)*128,
     +       dir*255,ensfile*255,units(nvarmax)*60,var*60,unit*60,
     +       tmpunits*60,tmpvars*60,string1*42,string2*42
        integer iargc,llen,rindex
        external llen,rindex
*
        n = iargc()
        if ( lag1.lt.0 ) print *,'(point leading field)'
        if ( lag2.gt.0 ) print *,'(field leading point)'
	if ( dump ) write(0,*)'correlatefield: dump not supported'
	if ( plot ) write(0,*)'correlatefield: plot not supported'
	if ( lks ) write(0,*)'correlatefield: K-S not supported'
	if ( lconting ) write(0,*)'correlatefield: contingency '//
     +        'tables not supported'
	do i=1,indxuse
	    if ( lincl(i) ) write(0,*)'correlatefield: what do ',
     +	    	'you mean with ',strindx(i),'?'
	enddo
        if ( composite ) then
            if ( maxdata.gt.1e33 .and. mindata.lt.-1e33 .and.
     +            (pmaxdata.le.0 .or. pmaxdata.ge.100) .and.
     +            (pmindata.le.0 .or. pmindata.ge.100) ) then
                write(0,*) 'correlatefield: cannot make composites '//
     +               'without cut-offs ',maxdata,mindata,pmaxdata
     +               ,pmindata
                write(*,*) 'correlatefield: cannot make composites '//
     +                'without cut-offs',maxdata,mindata,pmaxdata
     +               ,pmindata
                call abort
            endif
        endif
        call getarg(n,outfile)
        inquire(file=outfile,exist=lexist)
        if ( lexist ) then
            print *,'output file ',outfile(1:index(outfile,' ')-1),
     +            ' already exists, overwrite? [y/n]'
            read(*,'(a)') yesno
            if (  yesno.ne.'y' .and. yesno.ne.'Y' .and. 
     +            yesno.ne.'j' .and. yesno.ne.'J' ) then
                stop
            endif
            open(1,file=outfile)
            close(1,status='delete')
        endif
*       save time on the initialization - but not too much.
        yr1 = max(yr1,firstyr,firstyr - (min(lag1,lag2)+nperyear-1)
     +        /nperyear)
        yr2 = min(yr2,firstyr + (firstmo+nt-2)/nperyear,
     +        firstyr + (firstmo+nt-2)/nperyear - (max(lag1,lag2)
     +        -nperyear+1)/nperyear)
        if ( nyrwindow.gt.0 ) then
            write(0,'(2a,i6,a)') ,'p-values are computed'
     +           ,' against a ',nmc,' sample Monte Carlo<br>'
        endif
*
*       read kill info file (in the climexp) and add own PID
*
        call killfile(dir,title,infile,0)
*
*       allocate time series arrays
*
        allocate(data(mpermax,yrbeg:yrend,0:nensmax))
        allocate(mcdata(mpermax,yrbeg:yrend,0:nensmax))
        allocate(fxy(mpermax,yrbeg:yrend,0:nensmax))
        allocate(lnewyr(mpermax,yrbeg:yrend))
        if ( ncrossvalidate.gt.0 ) then
            mdatmax = nperyear*(lastyr-firstyr+1)*(nens2-nens1+1)
            allocate(aaa(mdatmax))
            allocate(bbb(mdatmax))
            allocate(aaa1(nx,ny,nz,nperyear,yr1:yr2))
            allocate(bbb1(nx,ny,nz,nperyear,yr1:yr2))
            aaa1 = 3e33
            bbb1 = 3e33
        end if
*
*       init
*
	print *,'init'
        nold = -1
        zold = 3e33
        sxxold = 3e33
	if ( lwrite ) print *,'filling fields with absent'
        if ( m1.ne.m2 .and. lag1.ne.lag2 ) then
            print *,'Sorry, can only handle either lags varying or'//
     +            ' months varying, not both'
            print *,'(months:',m1,m2,', lags:',lag1,lag2,')'
            call abort
        endif
        if ( lag2-lag1.gt.nperyear ) then
            print *,'Sorry, can only store ',nperyear+1,
     +           ' fields maximum'
            call abort
        endif
        do iens=nens1,nens2
            if ( lwrite ) print *,'calling makeabsent data',iens
            call makeabsent(data(1,yrbeg,iens),mpermax,yrbeg,yrend)
            if ( lwrite ) print *,'calling makeabsent fxy ',iens
            call makeabsent(fxy(1,yrbeg,iens),mpermax,yrbeg,yrend)
        enddo
*
*       read field, change absent values to our convention
*       if subtracting fields, read the whole field; otherwise just the
*       part needed
*
        if ( lsubtract ) then
            yr1a = firstyr
            yr2a = lastyr
        else
            yr1a = yr1
            yr2a = yr2
        endif
        do iens=nens1,nens2
            call keepalive(iens-nens1+1,nens2-nens1+1)
            if ( ncid.eq.-1 ) then
                call getarg(1,infile)
                if ( ensemble ) then
                    call filloutens(infile,iens)
                endif
                if ( lwrite ) print *,'calling parsectl on ',
     +               trim(infile)
                call parsectl(infile,datfile,nxmax,nx,xx,nymax,ny,yy
     +               ,nzmax,nz,zz,nt,nperyear,fyr,firstmo,undef
     +               ,endian,title,1,nvars,vars,ivars,lvars,units)
                if (lwrite) print '(2a)','# looking for ',trim(datfile)
                inquire(file=datfile,exist=lexist)
                if ( .not.lexist ) then
                    if (lwrite) print '(3a)','# looking for '
     +                   ,trim(datfile),'.gz'
                    inquire(file=trim(datfile)//'.gz',exist=lexist)
                    if ( .not.lexist ) then
                        nens2 = iens-1
                        if ( nens2.ge.nens1 ) then
                            write(0,*) 'Found ensemble 0 to ',nens2
     +                           ,'<br>'
                            goto 5
                        else
                            write(0,*) 'Cannot locate file '
     +                            ,datfile(1:llen(dir))
                            call abort
                        endif
                    endif
                endif
                if ( lwrite ) then
                    print *,'opening file ',trim(datfile)
                endif
                call zreaddatfile(datfile,field(1,1,1,1,firstyr,iens),
     +                nx,ny,nz,nx,ny,nz,nperyear,firstyr,lastyr,
     +               fyr,firstmo,nt,undef,endian,lwrite,yr1a,yr2a,
     +               1,1)
            else
                call getarg(1,infile)
                if ( ensemble ) then
                    call filloutens(infile,iens)
                endif
                if ( lwrite ) print *,'calling parsenc on ',
     +               infile(1:llen(infile))
                status = nf_open(infile,nf_nowrite,ncid)
                call parsenc(infile,ncid,nxmax,nx,xx,nymax,ny,yy
     +               ,nzmax,nz,zz,nt,nperyear,fyr,firstmo,
     +               undef,title,1,nvars,vars,jvars,lvars,units)
                call zreadncfile(ncid,field(1,1,1,1,firstyr,iens),
     +               nx,ny,nz,nx,ny,nz,nperyear,firstyr,lastyr,fyr,
     +               firstmo,nt,undef,lwrite,yr1a,yr2a,jvars)
            endif
        enddo
 5      continue
        if ( lstandardunits ) then
*           convert to standard units
            do iens=nens1,nens2
                tmpunits = units(1) ! they are otherwise adjusted
                call makestandardfield(field(1,1,1,1,firstyr,iens),nx,ny
     +               ,nz,nperyear,firstyr,lastyr,nx,ny,nz,nperyear
     +               ,max(firstyr,yr1a),min(lastyr,yr2a),vars(1)
     +               ,tmpunits,lwrite)
            enddo
            units(1) = tmpunits
            if ( lwrite ) then
                print *,'correlatefield: just after standard units'
                print *,'field(',(nx+1)/2,(ny+1)/2,(nz+1)/2,firstmo
     +               ,max(firstyr,yr1),') = ',field((nx+1)/2,(ny+1)/2
     +               ,(nz+1)/2,firstmo,max(firstyr,yr1),nens1)
            endif
        endif
*
*       apply land/sea mask
*
        call applylsmask(field,lsmask,nx,ny,nz,nperyear,firstyr,lastyr,
     +       nens1,nens2,lsmasktype,lwrite)
*
*       read series
*
        call getarg(iarg,file)
        if ( index(file,'%%').eq.0 .and. index(file,'++').eq.0 ) then
            ensseries = .FALSE.
            nens2series = nens1
            print *,'reading file ',file(1:index(file,' ')-1)
            call readseries(file,data(1,yrbeg,nens1),mpermax,yrbeg,yrend
     +           ,n,var,unit,lstandardunits,lwrite)
        else
            ensseries = .TRUE.
            if ( ensemble ) then
                nens2series = nens2
            else
                write(0,*) 'Cannot correlate a single field with an '//
     +                'ensemble series yet'
                call abort
            endif
            do iens=nens1,nens2series
                ensfile=file
                call filloutens(ensfile,iens)
                if ( lwrite ) write(0,*) 'looking for file '
     +                ,ensfile(1:llen(ensfile))
                inquire(file=ensfile,exist=lexist)
                if ( .not.lexist ) goto 10
                print *,'reading file ',ensfile(1:index(ensfile,' ')-1)
                call readseries(ensfile,data(1,yrbeg,iens),mpermax,yrbeg
     +               ,yrend,n,var,unit,lstandardunits,lwrite)
            enddo
            goto 11
   10       continue
            if ( iens.le.nens1 ) then
                write(0,*) 'error: cannot find time series'
                call abort
            endif
            write(0,*) 'Using series ensemble members ',nens1,' to '
     +            ,iens-1,'<br>'
            nens2series = iens - 1
            if ( nens2.lt.0 ) then
                write(0,*)
     +                'correlatefield: error: could not find ensemble '
     +                ,file(1:llen(file)),ensfile(1:llen(ensfile))
                call abort
            endif
   11       continue
        endif

        if ( n.ne.nperyear ) then
            write(0,*) 'correlatefield: error: cannot interpolate '//
     +            'in time (yet)',nperyear,n
            write(*,*) 'correlatefield: error: cannot interpolate '//
     +            'in time (yet)',nperyear,n
            call abort
        endif
        do iens=nens1,nens2series
*
*           take monthly anomalies
*
            if ( mdiff.gt.0 ) then
                if ( lwrite ) print *,'taking monthly anomalies'
                call mdiffit(data(1,yrbeg,iens),mpermax,nperyear,yrbeg
     +                ,yrend,mdiff)
            endif
*
*           sum series
*
            if ( lsum.gt.1 ) then
                if ( lwrite ) print *,'taking sum'
                call sumit(data(1,yrbeg,iens),mpermax,nperyear,yrbeg
     +                ,yrend,lsum,oper)
            endif
*
*           log, sqrt
*
            if ( logscale ) then
                call takelog(data(1,yrbeg,iens),mpermax,nperyear,yrbeg
     +                ,yrend)
            endif
            if ( sqrtscale ) then
                call takesqrt(data(1,yrbeg,iens),mpermax,nperyear,yrbeg
     +                ,yrend)
            endif
*
*           detrend data
*
            if ( ldetrend ) then
                i = llen(file)
                if ( index(file,'time').eq.0 ) then
                    if ( lwrite ) print *,'detrending series'
                    call detrend(data(1,yrbeg,iens),mpermax,nperyear,
     +                   yrbeg,yrend,yr1,yr2,m1,m2,lsel)
                endif
            endif
*
*           differentiate data
*
            if ( ndiff.ne.0 ) then
                if ( lwrite ) print *,'Taking differences/averaging'
                call ndiffit(data(1,yrbeg,iens),mpermax,nperyear,yrbeg
     +                ,yrend,ndiff,minfacsum)
            endif
            if ( ndiff2.ne.0 ) then
                if ( lwrite ) print *,'Taking differences/averaging2'
                call ndiffit(data(1,yrbeg,iens),mpermax,nperyear,yrbeg
     +                ,yrend,ndiff2,minfacsum)
            endif
        enddo
*
*       composites
*
        if ( pmindata.gt.0 .and. pmindata.lt.100 .or.
     +       pmaxdata.gt.0 .and. pmaxdata.lt.100 ) then
            if ( fix2 ) then
                n = 0
            else
                if ( lag1.eq.lag2 ) then
                    n = lag1
                else
                    write(0,*) 'correlatefield: error: can only handle',
     +                   ' a single lag in composite analysis'
                    call abort
                end if
            end if
            if ( pmindata.gt.0 .and. pmindata.lt.100 ) then
                call getj1j2(j1,j2,m1,nperyear,.false.)
                call getenscutoff(mindata,pmindata,data,mpermax
     +               ,nperyear,yrbeg,yrend,nensmax,nens1,nens2series,yr1
     +               ,yr2,j1,j2,n)
                write(0,'(a,f6.2,a,g14.6,a)') 'Converted ',
     +               pmindata,'% to ',mindata,'<br>'
            endif
            if ( pmaxdata.gt.0 .and. pmaxdata.lt.100 ) then
                call getj1j2(j1,j2,m1,nperyear,.false.)
                call getenscutoff(maxdata,pmaxdata,data,mpermax
     +               ,nperyear,yrbeg,yrend,nensmax,nens1,nens2series,yr1
     +               ,yr2,j1,j2,n)
                write(0,'(a,f6.2,a,g14.6,a)') 'Converted ',
     +               pmaxdata,'% to ',maxdata,'<br>'
            endif
        end if
        if ( composite ) then
            call whencomposite(data,mpermax,yrbeg,yrend,nensmax
     +           ,nens2series,lnewyr,nperyear)
        endif
*
*       anomalies - necessary if we consider more than one month
*
        if ( .not.composite .and. anom ) then
            if ( lwrite ) print *,'Taking anomalies'
            do iens=nens1,nens2
                call anomal(data(1,yrbeg,iens),mpermax,nperyear,yrbeg
     +                ,yrend,yr1,yr2)
            enddo
        endif
*
*       ensemble anomalies
*
        if ( .not.composite .and. nens2series.gt.nens1 .and. lensanom )
     +        then
            if ( lwrite ) print *,'Taking anomalies wrt ensemble'
            call anomalensemble(data,mpermax,nperyear,yrbeg,yrend,
     +           yr1,yr2,nens1,nens2series)
        endif
*
*       adjust yr1,yr2 to not be longer than the time series
*
        call adjustyr(yr1,yr2,data(1,yrbeg,nens1),mpermax,nperyear,yrbeg
     +        ,yrend)
        call adjustyr(yr2,yr1,data(1,yrbeg,nens1),mpermax,nperyear,yrbeg
     +        ,yrend)
!
!       compute number of time steps with data (<=nt) for minfac cuts
!
        ntp = 0
        do yr=yr1,yr2
            do mo=1,nperyear
                do iens=nens1,nens2
                    do jz=1,nz
                        do jy=1,ny
                            do jx=1,nx
                                if ( field(jx,jy,jz,mo,yr,iens).lt.1e33
     +                               ) then
                                    ntp = ntp + 1
                                    goto 101
                                end if
                            end do
                        end do
                    end do
                end do
 101            continue
            end do
        end do
!       check for ntp a multiple of 365 while nperyaer = 366
        s = ntp/real(nperyear)
        if ( abs(s/nint(s)-1).lt.0.01 ) then
            ntp = nperyear*nint(s)
        end if
        if ( lwrite ) print *,'nt,ntp = ',nt,ntp
*
*       compute minfac if it has not been set explicitly
*
        if ( minfac.lt.0 .and. minnum.lt.0 ) then
*           compute percentage valid points in data
            n = 0
            nn = 0
            call getj1j2(j1,j2,m1,nperyear,.false.)
            do yr=yr1-1,yr2
                do jj=j1,j2
                    if ( fix2 ) then
                        j = jj+lag1
                    else
                        j = jj
                    endif
                    call normon(j,yr,i,nperyear)
                    if ( i.ge.yr1 .and. i.le.yr2 ) then
                        nn = nn + 1
                        if (  data(j,i,nens1).lt.1e33 .and. (
     +                        (data(j,i,nens1).le.maxdata) .eqv.
     +                        (data(j,i,nens1).ge.mindata) .eqv.
     +                        (maxdata.ge.mindata) ) ) then
                            n = n + 1
                        endif
                    endif
                enddo
            enddo
*           heuristic, gives 0.25 for 150 yrs, 0.5 for 50 yrs, 0.75 for 20yrs
            minfac = max(0.1,
     +           min(0.6,
     +           (1.5-log(1+real(min(ntp,nperyear*(yr2-yr1+1))-1)*n/nn
     +           /nperyear)/4)*n/nn))
        endif
        if ( minfac.gt.0 ) then
            write(0,'(a,i2,a)') 'Requiring at least ',
     +            nint(100*minfac),'% valid points<br>'
        endif
*
*       loop over grid points
*
	print *,'correlating'
        yrstart = yr2
        yrstop  = yr1
        do jz=1,nz
            do jy=1,ny
                do jx=1,nx
                    if ( nyrwindow.eq.0 ) then
                        if ( mod(jx,10).eq.1 ) 
     +                       call keepalive(jy+(jz-1)*ny,ny*nz)
                    else
                        call keepalive(jx+(jy-1)*nx+(jy-1)*(jz-1)*nx*ny
     +                       ,nx*ny*nz)
                    endif
                    do month=0,nperyear
                        r(jx,jy,jz,month) = absent
                        prob(jx,jy,jz,month) = absent
                        a(jx,jy,jz,month) = absent
                        b(jx,jy,jz,month) = absent
                        da(jx,jy,jz,month) = absent
                        db(jx,jy,jz,month) = absent
                        relregr(jx,jy,jz,month) = absent
                        xn(jx,jy,jz,month) = absent
                    enddo
*
*                   create 1-D series from field
*
                    n = 0
                    do iens=nens1,nens2
                        do i=yr1,yr2
                            do j=1,nperyear
                                fxy(j,i,iens) = field(jx,jy,jz,j,i,iens)
                                if ( fxy(j,i,iens).lt.0.9*absent ) 
     +                                n = n+1
                            enddo
                        enddo
                        if ( n.lt.3 ) then
                            if ( lwrite ) print '(a,3i5,a,3f7.2,a,i5,a)'
     +                           ,'not enough valid points at ',jx,jy,jz
     +                           ,' (',xx(jx),yy(jy),zz(jz),'): ',n
     +                           ,' < 3'
                            goto 800
                        endif
*
*                       take monthly anomalies
*
                        if ( mdiff2.gt.0 ) then
                            call mdiffit(fxy(1,yrbeg,iens),mpermax
     +                            ,nperyear,yrbeg,yrend,mdiff2)
                        endif
*       
*                       sum
*
                        if ( lsum2.gt.1 ) then
                            call sumit(fxy(1,yrbeg,iens),mpermax
     +                            ,nperyear,yrbeg,yrend,lsum2,'v')
                        endif
*
*                       log, sqrt
*
                        if ( logfield ) then
                            call takelog(fxy(1,yrbeg,iens),mpermax
     +                           ,nperyear,yrbeg,yrend)
                        endif
                        if ( sqrtfield ) then
                            call takesqrt(fxy(1,yrbeg,iens),mpermax
     +                           ,nperyear,yrbeg,yrend)
                        endif
*
*                       detrend
*
                        if ( ldetrend ) then
                            if ( lwrite ) print *,'Detrending field'
                            if ( lag1.eq.0 .and. lag2.eq.0 .or. m1.eq.0
     +                            .or.lsel.eq.12 ) then
                                call detrend(fxy(1,yrbeg,iens),mpermax
     +                               ,nperyear,yrbeg,yrend,yr1,yr2,m1,m2
     +                               ,lsel)
                            else
                                call detrend(fxy(1,yrbeg,iens),mpermax
     +                               ,nperyear,yrbeg,yrend,yr1,yr2,1
     +                               ,nperyear,lsel)
                            endif
                        endif
*       
*                       differentiate
*
                        if ( ndiff.ne.0 ) then
                            if ( lwrite ) print *
     +                            ,'Taking differences/averaging'
                            call ndiffit(fxy(1,yrbeg,iens),mpermax
     +                           ,nperyear,yrbeg,yrend,ndiff,minfacsum)
                        endif
                        if ( ndiff2.ne.0 ) then
                            if ( lwrite ) print *
     +                            ,'Taking differences/averaging2'
                            call ndiffit(fxy(1,yrbeg,iens),mpermax
     +                           ,nperyear,yrbeg,yrend,ndiff2,minfacsum)
                        endif
*       
*                       anomalies
*
                        if ( anom .or. lsel.gt.1 .and. ndiff.le.0 .and.
     +                        ndiff2.le.0 .or. composite ) then
                            if ( lwrite ) print *,'Taking anomalies'
                            call anomal(fxy(1,yrbeg,iens),mpermax
     +                            ,nperyear,yrbeg,yrend,yr1,yr2)
                        endif
                    enddo       ! loop over ensemble members
!
!                   anomalies wrt ensemble mean
!
                    if ( nens2.gt.nens1 .and. lensanom ) then
                        if ( lwrite ) print *
     +                       ,'Taking anomalies wrt ensemble mean'
                        call anomalensemble(fxy,mpermax,nperyear,
     +                       yrbeg,yrend,yr1,yr2,nens1,nens2)
                    endif
*
*                   copy ensemble members so that there is the same
*                   number of valid ones at every time step
*       
                    do j=0,nperyear
                        ndup(j) = 0
                    enddo
                    if ( ensemble ) then
***                        call makeensfull(ndup,nperyear,fxy,mpermax,yrbeg
***     +                        ,yrend,nens1,nens2,validens,lwrite)
                        do j=1,nperyear
                            do i=yrbeg,yrend
                                mens = 0
                                do iens=nens1,nens2
                                    if ( fxy(j,i,iens).lt.1e33 ) then
                                        mens = mens + 1
                                        validens(mens) = iens
                                    endif
                                enddo
                                if (  mens.gt.1 .and. 
     +                                mens.le.(nens2-nens1) ) then
*                                   some - but not all - are undefined
*                                   copy as many whole sets in as fit
                                    if ( lwrite ) print
     +                                    '(a,i4,i2,a,i4,a,2i4)','At ',i
     +                                    ,j,' only found ',mens
     +                                    ,' ensemble members out of '
     +                                    ,nens1,nens2
                                    k = (nens2-nens1+1)/mens
                                    iens = nens1
                                    do l=1,k-1
                                        do jens=1,mens
  300                                       continue
                                            if ( fxy(j,i,iens).lt.1e33
     +                                            ) then
                                                if ( iens.ge.nens2 )
     +                                                then
                                                    write(0,*) 'error: '
     +                                                    ,iens,nens2
                                                    do k=nens1,nens2
                                                        write(0,*) fxy(j
     +                                                        ,i,k)
                                                    enddo
                                                    call abort
                                                endif
                                                iens = iens + 1
                                                goto 300
                                            endif
                                            if ( .FALSE. .and. lwrite )
     +                                            print*,'copying '
     +                                            ,validens(jens),' to '
     +                                            ,iens,fxy(j,i
     +                                            ,validens(jens))
     +                                            ,fxy(j,i,iens)
                                            ndup(j) = ndup(j) + 1
                                            fxy(j,i,iens) = 
     +                                            fxy(j,i,validens(jens)
     +                                            )
                                            if ( ensseries ) then
                                                data(j,i,iens) =
     +                                                data(j,i
     +                                                ,validens(jens))
                                            endif
                                        enddo
                                    enddo
*                                   and fill out the rest with random members
                                    do l=nens1+k*mens,nens2
*                                       check...
                                        do k=1,mens
                                            if ( validens(k).ge.0 ) 
     +                                            goto 305
                                        enddo
                                        write(0,*)
     +                                        'error: no more members'
                                        call abort
  305                                   continue
  310                                   continue
                                        call random_number(xrand)
                                        k = 1 + int(mens*xrand)
                                        if ( k.le.0 .or. k.gt.mens .or.
     +                                        validens(k).lt.0 ) goto
     +                                        310
  320                                   continue
                                        if ( fxy(j,i,iens).lt.1e33 )
     +                                        then
                                            iens = iens + 1
                                            goto 320
                                        endif
                                        if ( .FALSE. .and. lwrite )
     +                                        print *,'copying member '
     +                                        ,validens(k),' to ',iens
     +                                        ,fxy(j,i,validens(k))
     +                                        ,fxy(j,i,iens)
                                        ndup(j) = ndup(j) + 1
                                        fxy(j,i,iens) = fxy(j,i
     +                                        ,validens(k))
                                        if ( ensseries ) then
                                            data(j,i,iens) = 
     +                                            data(j,i,validens(k))
                                        endif
                                        validens(k) = -1
                                    enddo
                                endif
                            enddo
                        enddo
                    endif
*
*                   correlate!
*
                    do month=m1,m2
                        call getj1j2(j1,j2,month,nperyear,.false.)
                        if ( ensemble ) then
                            ndup(0) = 0
                            do j=j1,j2
                                ndup(0) = ndup(0) + ndup(j)
                            enddo
                        endif
                        do lag=lag1,lag2
*       
*                           fill linear arrays without absent values and
*                           compute r
*
                            n = 0
                            do yr=yr1-1,yr2
                                do jj=j1,j2
                                    if ( fix2 ) then
                                        j = jj+lag
                                    else
                                        j = jj
                                    endif
                                    call normon(j,yr,i,nperyear)
                                    if ( i.lt.yr1 .or.i.gt.yr2 )
     +                                    goto 710
                                    m = j-lag
                                    call normon(m,i,ii,nperyear)
                                    if ( ii.lt.yrbeg .or.ii.gt.yrend )
     +                                    goto 710
                                    if ( yr.eq.yr1 .and. jx.eq.nx/2 .and
     +                                    .jy.eq.ny/2 .and. jz.eq.1 )
     +                                    then
                                        print '(a,i4,i5,i4,i5,a)'
     +                                       ,'Correlating months/days',
     +                                       j,i,m,ii
     +                                       ,' of point and field'
                                    endif
                                    do iens=nens1,nens2
                                        if ( ensseries ) then
                                            iens2 = iens
                                        else
                                            iens2 = nens1
                                        endif
                                        if ( fxy(m,ii,iens).lt.1e33
     +                                        .and. ( 
     +                                        (fxy(m,ii,iens).le.maxindx
     +                                        ) .eqv. 
     +                                        (fxy(m,ii,iens).ge.minindx
     +                                        ) .eqv. 
     +                                        (maxindx.ge.minindx) ) 
     +                                        .and. 
     +                                        data(j,i,iens2).lt.1e33
     +                                        .and.((
     +                                        data(j,i,iens2).le.maxdata
     +                                        ).eqv.(
     +                                        data(j,i,iens2).ge.mindata
     +                                        ).eqv.(
     +                                        maxdata.ge.mindata) )
     +                                        ) then
                                            n = n+1
                                            if ( n.gt.ndatmax ) goto 909
                                            ddata(n) = data(j,i,iens2)
                                            dindx(n) = fxy(m,ii,iens)
                                            yrmo(1,n) = i
                                            yrmo(2,n) = j
                                            yrstart = min(yrstart,i,ii)
                                            yrstop  = max(yrstop,i,ii)
                                        endif ! valid point
                                    enddo ! iens
  710                               continue
                                enddo ! jj
                            enddo ! yr
                            if ( m1.ne.m2 ) then
                                m = month-m1
                            elseif ( .not.fix2 ) then
                                m = lag2-lag
                            else
                                m = lag-lag1
                            endif
                            if ( jx.eq.nx/2 .and. jy.eq.ny/2 .and. 
     +                            jz.eq.1 ) then
                                print '(a,i2)','and storing it at ',m
                            endif
                            xn(jx,jy,jz,m) = n
                            ndiffn = 1 + max(0,-ndiff)
                            if ( month.eq.0 .and. n.lt.minfac*
     +                           min(ntp/ndiffn,nperyear*(yr2-yr1+1))
     +                           .or. month.ne.0 .and. n.lt.minfac*
     +                           min(ntp/nperyear/ndiffn,yr2-yr1+1)
     +                           .or. n.lt.minnum ) then
                                if ( lwrite ) print
     +                               '(a,3i5,a,3f7.2,a,2i3,a,2i6)'
     +                               ,'not enough valid points at ',jx
     +                               ,jy,jz,' (',xx(jx),yy(jy),zz(jz)
     +                               ,')',month,lag,': ',n,ntp
                                goto 790
                            endif
                            if ( ncrossvalidate.gt.0 ) then
                                if ( n.gt.mdatmax ) then
                                    write(0,*) 'error: mdatmax<n '
     +                                   ,mdatmax,n
                                    call abort
                                end if
                            end if
                            call fitcross(dindx,ddata,n,sig,0,
     +                           aa,a(jx,jy,jz,m),daa,da(jx,jy,jz,m)
     +                           ,chi2,q,ncrossvalidate,aaa,bbb,.true.)
                            if ( ncrossvalidate.gt.0 ) then
                                do i=1,n
                                    aaa1(jx,jy,jz,yrmo(2,i),yrmo(1,i)) =
     +                                   aaa(i)
                                    bbb1(jx,jy,jz,yrmo(2,i),yrmo(1,i)) =
     +                                   bbb(i)
                                end do
                            end if
                            call fitcross(ddata,dindx,n,sig,0,
     +                           bb,b(jx,jy,jz,m),dbb,db(jx,jy,jz,m)
     +                           ,chi2,q,ncrossvalidate,aaa,bbb,.false.)
                            noisetype = 1
                            if ( ensseries ) then
                                imens(0) = nens2
                            else
                                imens(0) = 0
                            endif
                            imens(1) = nens2
                            if ( .false. ) then
                                do k=1,n
                                    print *,k,dindx(k),ddata(k)
                                end do
                            end if
                            if ( index(file,'/time').ne.0 .or.
     +                           file(1:4).eq.'time') then
                                alpha = 0
                            else if ( fitfunc.gt.1 ) then
                            	write(0,*) 'cirrelatefield: warning: ',
     +                              'serial correlations are not taken'
     +				    ,' into account for non-linear fit'
				alpha = 0
                            else
                                call getred(alpha,j1,j2,lag,1,nperyear
     +                               ,imens,1,data,fxy,mpermax,yrbeg
     +                               ,yrend,nensmax,bb,b(jx,jy,jz,m))
                            end if
                            if ( lwrite ) then
                                print *,'j1,j2 = ',j1,j2
                                call getdf(df,month,n,ndup(0),decor,
     +                               max(lsum,lsum2),max(1,1-ndiff,
     +                               1-ndiff2),nperyear)
                                print *,'df was',df,n,ndup(0)
                            endif
                            if ( alpha.eq.1 .or. n-ndup(0).le.0 ) then
                                df = -2
                            elseif ( alpha.gt.2/sqrt(real(n-ndup(0))) 
     +                               .and.alpha.gt.exp(-1.) ) then
                                df = n/(1 - 1/(2*log(alpha))) - 2
                            else
                                df = n - 2
                            endif
                            if ( lwrite ) then
                                print *,'df  is',df,alpha
                            endif
                            if ( df.lt.1 ) then
                                df = 3e33
                            end if
                            if ( df.lt.1 ) then
                                r(jx,jy,jz,m) = 3e33
                                prob(jx,jy,jz,m) = 3e33
                                a(jx,jy,jz,m) = 3e33
                                da(jx,jy,jz,m) = 3e33
                                b(jx,jy,jz,m) = 3e33
                                db(jx,jy,jz,m) = 3e33
                                a1(jx,jy,jz,m) = 3e33
                                da1(jx,jy,jz,m) = 3e33
                                if ( lwrite ) print '(a,f5.2,3i5)'
     +                               ,'error: df < 1: ',df,jx,jy,jz
                                goto 790
                            endif
                            if ( lfitnoise ) then
                                call fitnoisemodel(ddata,dindx,n,
     +                               a1(jx,jy,jz,m),da1(jx,jy,jz,m),
     +                               a(jx,jy,jz,m),da(jx,jy,jz,m),
     +                               b(jx,jy,jz,m),db(jx,jy,jz,m),
     +                               r(jx,jy,jz,m),prob(jx,jy,jz,m),
     +                               cov(jx,jy,jz,m),lbootstrap,lwrite)
                            elseif ( lrank ) then
                                if ( lwrite ) then
                                    if ( month.eq.0 ) then
                                        sum = max(lsum,lsum2) + decor
                                    else
                                        sum = 1 + (max(lsum,lsum2)-1)
     +                                       /nperyear +decor/nperyear
                                    endif
                                    print *,'sum was: ',sum
                                endif
                                sum = n/(df+2)
                                if ( lwrite ) then
                                    print *,'sum  is: ',sum
                                endif
                                call spearx(ddata,dindx,n,ddata,dindx,d
     +                                ,zd,probd,r(jx,jy,jz,m),
     +                                prob(jx,jy,jz,m),sum,adata,sxx
     +                                ,aindx,syy)
                            elseif ( composite ) then
                                call makecomposite(ddata,dindx,n,
     +                               b(jx,jy,jz,m),db(jx,jy,jz,m),
     +                               prob(jx,jy,jz,m),df,lwrite)
                            elseif ( fitfunc.eq.1 ) then
                                 if ( df.le.0 ) then
                                    if ( lwrite ) print '(a,f5.2,3i5)'
     +                                ,'error: df <= 0: ',df,jx,jy,jz
                                    goto 790
                                endif
!                               adjust error estimates from fit()
                                da(jx,jy,jz,m) = da(jx,jy,jz,m)
     +                               *sqrt((n-2)/df)
                                db(jx,jy,jz,m) = db(jx,jy,jz,m)
     +                               *sqrt((n-2)/df)
                                call pearsncross(ddata,dindx,n,
     +                               r(jx,jy,jz,m),prob(jx,jy,jz,m),z,
     +                               adata,sxx,aindx,syy,sxy,df
     +                               ,ncrossvalidate)
                                if ( sxx.eq.0 .or. syy.eq.0 ) then
                                    r(jx,jy,jz,m) = absent
                                    prob(jx,jy,jz,m) = absent
                                else
                                    cov(jx,jy,jz,m) = sxy/(n-1)
                                endif
                                if ( abs(aindx).gt.1e-33 .and. 
     +                                   aindx.lt.1e33 .and.
     +                                   b(jx,jy,jz,m).lt.1e33 ) then
                                    relregr(jx,jy,jz,m) = 
     +                                   b(jx,jy,jz,m)/aindx
                                else
                                    relregr(jx,jy,jz,m) = 3e33
                                endif
                                if ( nyrwindow.gt.0 ) then
                                    if ( ensseries ) then
                                        imens(0) = nens2
                                    else
                                        imens(0) = 0
                                    endif
                                    imens(1) = nens2
                                    if ( nfittime.ne.0 ) goto 901
                                    call getruncorr(dindx,ddata,lfirst
     +                                   ,dddata,ndatmax,j1,j2,lag,1
     +                                   ,month,nperyear,imens,1,fxy
     +                                   ,data,mpermax,yrbeg,yrend
     +                                   ,nensmax,ndup(0),filter
     +                                   ,' ',.false.,.false.
     +                                   ,rmin(jx,jy,jz,m)
     +                                   ,rmax(jx,jy,jz,m)
     +                                   ,zdif(jx,jy,jz,m))
                                    if ( zdif(jx,jy,jz,m).gt.1e30 ) then
                                        rprob(jx,jy,jz,m) = absent
                                        goto 750
                                    endif
*
*                                   Monte Carlo to determine significance
*
                                    call filllinarray(dindx,ddata
     +                                   ,lfirst,dddata,ndatmax,n,j1
     +                                   ,j2,lag,1,nperyear,imens,1
     +                                   ,fxy,data,mpermax,yrbeg
     +                                   ,yrend,nensmax,filter,-999,-999
     +                                   ,yrmo)
                                    if ( lwrite ) print *
     +                                   ,'correlate: a,b,sd,r,sde '
     +                                   ,a(jx,jy,jz,m)
     +                                   ,b(jx,jy,jz,m)
     +                                   ,sqrt(sxx/(n-1))
     +                                   ,r(jx,jy,jz,m)
     +                                   ,sqrt(sxx/(n-1)*(1-r(jx,jy
     +                                   ,jz,m)**2))
*
*                                   copy the data to a scratch array
*                                   to indicate where the holes are
*
                                    if ( ensseries ) then
                                        do iens=nens1,nens2
                                            do i=yr1,yr2
                                                do j=1,12
                                                    mcdata(j,i,iens) =
     +                                                   data(j,i,iens)
                                                enddo
                                            enddo
                                        enddo
                                    else
                                        do i=yr1,yr2
                                            do j=1,12
                                                mcdata(j,i,0) = 
     +                                               data(j,i,0)
                                            enddo
                                        enddo
                                    endif
                                    do i=1,nmc
                                        call makemcseries(mcdata,fxy
     +                                       ,mpermax,yrbeg,yrend
     +                                       ,nensmax,1,nperyear,1
     +                                       ,lag,j1,j2,imens,adata
     +                                       ,sxx,aindx,syy,sxy
     +                                       ,alpha,n)
                                        if ( lwrite ) then
                                            lwrite = .false.
                                            call filllinarray(dindx
     +                                           ,ddata,lfirst
     +                                           ,dddata,ndatmax,n
     +                                           ,j1,j2,lag,1
     +                                           ,nperyear,imens,1
     +                                           ,fxy,mcdata,mpermax
     +                                           ,yrbeg,yrend
     +                                           ,nensmax,filter,-999,
     +                                           -999,yrmo)
                                            call printcorr(dindx
     +                                           ,ddata,lfirst
     +                                           ,dddata,yrmo,n,ndup(0)
     +                                           ,j1,j2,month
     +                                           ,nperyear,lag,' ',
     +                                           .false.,.false.
     +                                           ,results(i),dresult
     +                                           ,dum)
                                            lwrite = .true.
                                        endif
                                        call getruncorr(dindx,ddata
     +                                       ,lfirst,dddata,ndatmax
     +                                       ,j1,j2,lag,1,month
     +                                       ,nperyear,imens,1,fxy
     +                                       ,mcdata,mpermax,yrbeg
     +                                       ,yrend,nensmax,ndup(0)
     +                                       ,filter,' ',.false.,
     +                                       .false.,rmins(i)
     +                                       ,rmaxs(i),zdifs(i))
                                        if ( lwrite ) print *,'zdif = '
     +                                       ,zdifs(i)
                                    enddo
                                    if ( lwrite ) call getsign('result '
     +                                   ,r(jx,jy,jz,m),results,nmc,1
     +                                   ,rprob(jx,jy,jz,m),lwrite)
                                    call getsign('zdif',zdif(jx,jy,jz,m)
     +                                   ,zdifs,nmc,1,rprob(jx,jy,jz,m),
     +                                   lwrite.or..true.)
 750                                continue
                                endif
                            else ! fitfunc>1, nonlinear fit
                                if ( fitfunc.lt.1 .or. fitfunc.gt.3 )
     +                               then
                                    write(0,*) 'correlatefield: error: '
     +                                   //'cannot handle polynomial '//
     +                                   'of degree ',fitfunc,' yet'
                                    write(*,*) 'correlatefield: error: '
     +                                   //'cannot handle polynomial '//
     +                                   'of degree ',fitfunc,' yet'
                                    call abort
                                endif
                                call getdf(df,month,n,ndup(0),decor,
     +                               max(lsum,lsum2),
     +                               max(1,1-ndiff,1-ndiff2),nperyear)
                                if ( df.le.0 ) then
                                    if ( lwrite ) print '(a,f5.2,3i5)'
     +                                ,'error: df <= 0: ',df,jx,jy,jz
                                    goto 790
                                endif
                                if ( fitfunc.eq.2 ) then
                                    call fit2(ddata,dindx,n,
     +                                   a(jx,jy,jz,m),b(jx,jy,jz,m),
     +                                   a1(jx,jy,jz,m),
     +                                   da(jx,jy,jz,m),db(jx,jy,jz,m),
     +                                   da1(jx,jy,jz,m),
     +                                   r(jx,jy,jz,m),prob(jx,jy,jz,m),
     +                                   df,ncrossvalidate,lwrite)
                                else
                                    write(0,*) 'correlatefield: error:'/
     +                                   /' cubic fitting ot yet ready'
                                    write(*,*) 'correlatefield: error:'/
     +                                   /' cubic fitting ot yet ready'
                                    call abort
                                endif
                            endif
                            if ( lwrite .and. fitfunc.eq.1 )
     +                          print '(a,3i5,2i3,a,i6,a,2f8.4,4f12.4)'
     +                           ,'point ',
     +                            jx,jy,jz,month,lag,' OK (',n,'): ',
     +                            r(jx,jy,jz,m),prob(jx,jy,jz,m),
     +                            a(jx,jy,jz,m),da(jx,jy,jz,m),
     +                            b(jx,jy,jz,m),db(jx,jy,jz,m)
                            if ( lwrite .and. fitfunc.gt.1 )
     +                          print '(a,3i5,2i3,a,i6,a,2f8.4,6f12.4)'
     +                           ,'point ',
     +                            jx,jy,jz,month,lag,' OK (',n,'): ',
     +                            r(jx,jy,jz,m),prob(jx,jy,jz,m),
     +                            a(jx,jy,jz,m),da(jx,jy,jz,m),
     +                            b(jx,jy,jz,m),db(jx,jy,jz,m),
     +                            a1(jx,jy,jz,m),da1(jx,jy,jz,m)
 790                        continue ! valid point/month
                        enddo   ! lag
                    enddo       ! month
  800               continue    ! valid point
                enddo           ! nx
            enddo               ! ny
        enddo                   ! nz
        if ( index(outfile,'.ctl').ne.0 ) then
            i = index(outfile,'.ctl')
            datfile = outfile(:i-1)//'.grd'
            open(2,file=datfile,form='unformatted',access='direct'
     +            ,recl=recfac*nx*ny*nz,err=920)
        elseif ( lsubtract ) then
            write(*,*)
     +           'correlatefield: error: netcdf output not yet ready'
            write(0,*)
     +           'correlatefield: error: netcdf output not yet ready'
            call abort
        endif
        if ( .not.lsubtract ) then
            call getenv('DIR',dir)
            ldir = llen(dir)
            if ( ldir.eq.0 ) ldir=1
            if ( dir(ldir:ldir).ne.'/' ) then
                ldir = ldir + 1
                dir(ldir:ldir) = '/'
            endif
            call args2title(title)
            tmpunits = units(1)
            tmpvars = vars(1)
            units = ' '
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
                nvars = 8
                vars(1) = 'corr'
                if ( lrank ) then
                    lvars(1) = 'rank correlation'
                else
                    lvars(1) = 'correlation'
                endif
                units(1) = '1'
                vars(2) = 'prob'
                lvars(2) = 'p-value'
                units(2) = '1'
                vars(3) = 'regr1'
                lvars(3) = 'regression of series '//trim(var)//
     +               ' on field '//trim(tmpvars)
                call divideunits(units(3),unit,tmpunits)
                vars(4) = 'regr'
                lvars(4) = 'regression of field '//trim(tmpvars)//
     +               ' on series '//trim(var)
                call divideunits(units(4),tmpunits,unit)
                vars(5) = 'd'//vars(3)
                lvars(5) = 'error on '//lvars(3)
                units(5) = units(3)
                vars(6) = 'd'//vars(4)
                lvars(6) = 'error on '//lvars(4)
                units(6) = units(4)
                if ( fitfunc.eq.1 ) then
                    vars(7) = 'n'
                    vars(8) = 'relregr'
                    lvars(7) = 'number of valid points'
                    lvars(8) = 'relative regression'
                    units(7) = '1'
                    units(8) = '1'
                elseif ( fitfunc.eq.2 ) then
                    vars(7) = 'quad'
                    lvars(7) = 'coefficient of quadratic term'
                    vars(8) = 'errorquad'
                    lvars(8) = 'error on coefficient of quadratic term'
                else
                    write(0,*) 'error: fitfunc = ',fitfunc
     +                   ,' not yet supported'
                    write(*,*) 'error: fitfunc = ',fitfunc
     +                   ,' not yet supported'
                    call abort
                endif
            endif
            if ( nyrwindow.gt.0 ) then
                if ( irunvar.eq.1 ) then
                    vars(nvars+1) = 'rmin'
                    vars(nvars+2) = 'rmax'
                    vars(nvars+3) = 'zdif'
                elseif ( irunvar.eq.2 ) then
                    vars(nvars+1) = 'bmin'
                    vars(nvars+2) = 'bmax'
                    vars(nvars+3) = 'bdif'
                else
                    write(0,*) 'correlatefield: unknown irunvar '
     +                   ,irunvar
                    write(*,*) 'correlatefield: unknown irunvar '
     +                   ,irunvar
                    call abort
                endif
                vars(nvars+4) = 'rprob'
                lvars(nvars+1) = 'minimum running correlation'
                lvars(nvars+2) = 'maximum running correlation'
                lvars(nvars+3) = 'difference between min and max'
                lvars(nvars+4) = 'probability of difference'
                nvars = nvars + 4
            endif
            do i=1,nvars
                if ( nz.gt.1 ) then
                    ivars(1,i) = nz
                else
                    ivars(1,i) = 0
                endif
                ivars(2,i) = 99
            enddo
*           give correlations dates in 0-1
            if ( m1.eq.0 ) then
                i = 0
            else
                i = 1
            endif
            j = m1-lag2
            if ( j.le.0 ) then
                j = j + nperyear*(1-j/nperyear)
            elseif ( j.gt.nperyear ) then
                j = j - nperyear*((j-1)/nperyear)
            endif
            if ( nperyear.gt.12 ) then
!               rough approximation
                j = 1 + 12*(j-1)/nperyear
            end if
            if ( index(outfile,'.ctl').ne.0 ) then
                call writectl(outfile,datfile,nx,xx,ny,yy,nz,zz
     +               ,1+(m2-m1)+(lag2-lag1),min(12,nperyear),i,j,3e33
     +               ,title,nvars,vars,ivars,lvars,units)
            else
                call writenc(outfile,ncid,ntvarid,itimeaxis,ntmax,nx,xx
     +               ,ny,yy,nz,zz,1+(m2-m1)+(lag2-lag1),min(12,nperyear)
     +               ,i,j,3e33,title,nvars,vars,ivars,lvars,units,0,0)
            endif
*
*           write output field in GrADS or netCDF format
*
            print *,'writing output'
            do lag=lag2,lag1,-1
                do mo=m1,m2
                    if ( m1.ne.m2 ) then
                        m = mo-m1
                    else
                        m = lag2-lag
                    endif
                    if ( index(outfile,'.ctl').ne.0 ) then
                        if ( lwrite ) then
                            print *,'writing records ',nvars*m+1,'-'
     +                            ,nvars*(m+1),' of fields ',m
     +                            ,' of size ',nx*ny*nz*recfa4
                            do jz=1,nz
                                do jy=1,ny
                                    do jx = 1,nx
                                        if ( abs(r(jx,jy,jz,m)).le.1 )
     +                                        then
                                            i = i + 1
                                            d = d + abs(r(jx,jy,jz,m))
                                        endif
                                    enddo
                                enddo
                            enddo
                            if ( i.gt.0 ) then
                                print *,'there are ',i
     +                               ,' valid values in record ',m
     +                               ,' with mean value ',d/i
                            else
                                print *,'there are ',i
     +                               ,' valid values in record ',m
                            endif
                        endif
                        write(2,rec=nvars*m+1)
     +                    (((r(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                        write(2,rec=nvars*m+2)
     +                    (((prob(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                        write(2,rec=nvars*m+3)
     +                    (((a(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                        write(2,rec=nvars*m+4)
     +                    (((b(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                        write(2,rec=nvars*m+5)
     +                    (((da(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                        write(2,rec=nvars*m+6)
     +                    (((db(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                        if ( fitfunc.eq.1 ) then
                            write(2,rec=nvars*m+7)
     +                    (((xn(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                            if ( lfitnoise ) then
                                write(2,rec=nvars*m+8)
     +                    (((a1(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                                write(2,rec=nvars*m+9)
     +                    (((da1(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                            elseif ( .not.composite ) then
                                write(2,rec=nvars*m+8)
     +                 (((relregr(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                            endif
                        else
                            write(2,rec=nvars*m+7)
     +                    (((a1(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                            write(2,rec=nvars*m+8)
     +                    (((da1(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                        endif
                        if ( nyrwindow.gt.0 ) then
                            write(2,rec=nvars*m+9)
     +                    (((rmin(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                            write(2,rec=nvars*m+10)
     +                    (((rmax(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                            write(2,rec=nvars*m+11)
     +                    (((zdif(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                            write(2,rec=nvars*m+12)
     +                   (((rprob(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                        endif
                    else
*                       netCDF file
                        call writencslice(ncid,0,0,0,ivars(1,1)
     +                       ,r(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                        call writencslice(ncid,0,0,0,ivars(1,2)
     +                       ,prob(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                        call writencslice(ncid,0,0,0,ivars(1,3)
     +                       ,a(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                        call writencslice(ncid,0,0,0,ivars(1,4)
     +                       ,b(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                        call writencslice(ncid,0,0,0,ivars(1,5)
     +                       ,da(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                        call writencslice(ncid,0,0,0,ivars(1,6)
     +                       ,db(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                        call writencslice(ncid,0,0,0,ivars(1,7)
     +                       ,xn(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                        if ( fitfunc.eq.1 ) then
                            if ( lfitnoise ) then
                                call writencslice(ncid,0,0,0,ivars(1,8)
     +                           ,a1(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                                call writencslice(ncid,0,0,0,ivars(1,9)
     +                           ,da1(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                            elseif ( .not.composite ) then
                                call writencslice(ncid,0,0,0,ivars(1,8)
     +                           ,relregr(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1
     +                           ,1)
                            endif
                        else
                            call writencslice(ncid,0,0,0,ivars(1,7)
     +                           ,a1(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                            call writencslice(ncid,0,0,0,ivars(1,8)
     +                           ,da1(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                        endif
                        if ( nyrwindow.gt.0 ) then
                           call writencslice(ncid,0,0,0,ivars(1,nvars-3)
     +                           ,rmin(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                           call writencslice(ncid,0,0,0,ivars(1,nvars-2)
     +                          ,rmax(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                           call writencslice(ncid,0,0,0,ivars(1,nvars-1)
     +                          ,zdif(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                           call writencslice(ncid,0,0,0,ivars(1,nvars)
     +                          ,rprob(1,1,1,m),nx,ny,nz,nx,ny,nz,m+1,1)
                        endif
                    endif
                enddo
            enddo
            if ( index(outfile,'.ctl').eq.0 ) then
                i = nf_close(ncid)
            endif
            if ( ncrossvalidate.gt.0 ) then
                title = 'data for subtractions of '//
     +               'cross-validated regression'
                vars(1) = 'a'
                lvars(1) = 'constant in cross-validated regression'
                units(1) = inunits(1)
                open(1,file=trim(bbfile)//'_a.ctl')
                close(1,status='delete')
                call writectl(trim(bbfile)//'_a.ctl',
     +               trim(bbfile)//'_a.grd',nx,xx,ny,yy,nz,zz,
     +               nperyear*(yr2-yr1+1),nperyear,yr1,1,3e33,title,
     +               1,vars,ivars,lvars,units)
                open(1,file=trim(bbfile)//'_a.grd',access='direct',
     +               form='unformatted',recl=nx*ny*nz*recfa4)
                irec = 0
                do i=yr1,yr2
                    do j=1,nperyear
                        irec = irec + 1
                        write(1,rec=irec) (((aaa1(jx,jy,jz,j,i),
     +                       jx=1,nx),jy=1,ny),jz=1,nz)
                    end do
                end do
                close(1)
                vars(1) = 'b'
                lvars(1) = 'regr in cross-validated regression'
                units(1) = units(4)
                open(1,file=trim(bbfile)//'_b.ctl')
                close(1,status='delete')
                call writectl(trim(bbfile)//'_b.ctl',
     +               trim(bbfile)//'_b.grd',nx,xx,ny,yy,nz,zz,
     +               nperyear*(yr2-yr1+1),nperyear,yr1,1,3e33,title,
     +               1,vars,ivars,lvars,units)
                open(1,file=trim(bbfile)//'_b.grd',access='direct',
     +               form='unformatted',recl=nx*ny*nz*recfa4)
                irec = 0
                do i=yr1,yr2
                    do j=1,nperyear
                        irec = irec + 1
                        write(1,rec=irec) (((bbb1(jx,jy,jz,j,i),
     +                       jx=1,nx),jy=1,ny),jz=1,nz)
                    end do
                end do
                close(1)
            end if
        elseif ( ensemble .and. nens1.ne.nens2 ) then ! subtract
            print *,'correlatefield: subtraction for ensembles '//
     +            'not yet ready'
            call abort
        else                    ! subtract
*
*           subtract best fit from field and give this as output
*
            if ( lwrite ) print *,'Subtracting best fit from field'
            nrec = 0
            allocate(field2(nx,ny,nz,nperyear,firstyr:lastyr,nens1:nens2
     +           ))
            field2 = 3e33
            do yr=firstyr,lastyr
*               skip month 0, which is everything together, if there is more
                if ( m2.ne.0 ) then
                    mm = max(1,m1)
                else
                    mm = m1
                endif
                do month=mm,m2
                    if ( month.eq.0 ) then
                        j1 = 1
                        j2 = nperyear
                    else
                        j1 = month
                        j2 = month + lsel - 1
                    endif
                    do lag=lag1,lag2
                        do jj=j1,j2
                            if ( fix2 ) then
                                j = jj+lag
                            else
                                j = jj
                            endif
                            call normon(j,yr,i,nperyear)
                            if ( i.lt.yrbeg .or.i.gt.yrend )
     +                            goto 910
                            if ( m1.ne.m2 ) then
                                m = month-m1
                            elseif ( .not.fix2 ) then
                                m = lag2-lag
                            else
                                m = lag-lag1
                            endif
                            do kk=j-lag,j-lag+lsum-1
                                k = kk
                                call normon(k,i,ii,nperyear)
                                if ( ii.lt.firstyr .or.ii.gt.lastyr )
     +                               goto 910
                                if ( yr.eq.yr1 ) then
                                    print '(a,i3,i5,i3,i5,a)'
     +                                   ,'Subtracting months',j,i,k,ii
     +                                   ,' of point and field'
                                endif
                                do jz=1,nz
                                    do jy=1,ny
                                        do jx=1,nx
                                            if ( b(jx,jy,jz,m).lt.1e33
     +                                           .and.
     +                                           field(jx,jy,jz,k,ii,
     +                                           nens1).lt.1e33
     +                                           .and.
     +                                           data(j,i,nens1).lt.
     +                                           1e33 ) then
                                                if ( .false. .and.
     +                                               jx.eq.(nx+1)/2.and.
     +                                               jy.eq.(ny+1)/2.and.
     +                                               nz.eq.(nz+1)/2 )
     +                                               then
                                                    print *,'field was '
     +                                                   ,field(jx,jy,jz
     +                                                   ,k,ii,nens1),k
     +                                                   ,ii,nens1
                                                    print *,' b,data = '
     +                                                   ,b(jx,jy,jz,m)
     +                                                   ,data(j,i,nens1
     +                                                   )
                                                end if
                                                field2(jx,jy,jz,k,ii
     +                                               ,nens1) = field(jx
     +                                               ,jy,jz,k,ii,nens1)
     +                                               - b(jx,jy,jz,m)
     +                                               *data(j,i,nens1)
                                                if ( .false. .and.
     +                                               jx.eq.(nx+1)/2.and.
     +                                               jy.eq.(ny+1)/2.and.
     +                                               nz.eq.(nz+1)/2 )
     +                                               then
                                                    print *,'field2 is '
     +                                                   ,field2(jx,jy
     +                                                   ,jz,k,ii,nens1)
                                                end if
                                            endif
                                        enddo ! jy
                                    enddo ! jx
                                enddo ! jz
                            enddo ! k
  910                       continue
                        enddo   ! jj
                    enddo       ! lag
                enddo           ! month
                do mo=1,nperyear
                    nrec = nrec + 1
                    if ( lwrite ) print *,'Writing new field for ',yr,mo
                    write(2,rec=nrec)
     +                   (((field2(jx,jy,jz,mo,yr,nens1),
     +                   jx=1,nx),jy=1,ny),jz=1,nz)
                enddo
            enddo               ! yr
            if ( index(outfile,'.ctl').ne.0 ) then
                title = intitle(1:llen(intitle))//
     +                ' with the effect of '//
     +                file(1+rindex(file,'/'):index(file,' ')-1)//
     +                ' linearly subtracted'
                nvars = 1
                vars(1) = invars(1)
                ivars(1,1) = nz
                ivars(2,1) = 99
                lvars(1) = trim(inlvars(1))//' with '//
     +               file(1+rindex(file,'/'):index(file,' ')-1)//
     +               ' linearly subtracted'
                units(1) = inunits(1)
                call writectl(outfile,datfile,nx,xx,ny,yy,nz,zz
     +                ,nrec,nperyear,firstyr,1,3e33,title,
     +                nvars,vars,ivars,lvars,units)
            endif
        endif                   ! subtract
        close(2)
        call savestartstop(yrstart,yrstop)
*
*       error messages
*
        goto 999
 901    write(0,*) 'correlatefield: error: fittime not yet implemented'
        call abort
 909    write(0,*) 'correlatefield: error: array too small ',n
        call abort
 920    write(0,*) 'correlatefield: error cannot open new '//
     +       'correlations file ',datfile(1:index(datfile,' ')-1)
        call abort
 999    continue
        end
*  #] cfield:
