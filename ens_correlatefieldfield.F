        program correlatefieldfield
*
*       program to correlate a field series to another field series on
*       the same grid to give fields of correlation coefficients,
*       probabilities that these are significant, and the fit
*       coefficients a, b and their errors.
*       when the option 'subtract' is given, the output consists of 
*       the first field minus the best fit to the second one
*       If one of the two fields is 3D and one 2D, the 2D field is 
*       used on all levels.
*
        implicit none
#include "params.h"
#include "recfac.h"
#include "netcdf.inc"
        integer nyrmax,nlevmax
        parameter(nyrmax=70,nlevmax=1)
        integer j,n,ncid1,ncid2,nx1,ny1,nz1,nt1,nper1,firstyr1,firstmo1,
     +        nx2,ny2,nz2,nt2,nper2,firstyr2,firstmo2,nvars,ivars1(2,1)
     +        ,ivars2(2,1),jvars1(6,nvmax),jvars2(6,nvmax),endian1
     +        ,endian2,status,nxf,nyf,nzf,nperyear,firstyr,lastyr,
     +        mens1,mens2
        real  xx1(nxmax),yy1(nymax),zz1(nzmax),
     +        xx2(nxmax),yy2(nymax),zz2(nzmax),
     +        field1(nxmax,nymax,12,nyrmax),
     +        field2(nxmax,nymax,12,nyrmax),
     +        u1,u2,
     +        r(nxmax,nymax,nlevmax,0:npermax),
     +        prob(nxmax,nymax,nlevmax,0:npermax),
     +        a(nxmax,nymax,nlevmax,0:npermax),
     +        b(nxmax,nymax,nlevmax,0:npermax),
     +        da(nxmax,nymax,nlevmax,0:npermax),
     +        db(nxmax,nymax,nlevmax,0:npermax),
     +        mean2(nxmax,nymax,nlevmax,0:npermax)
        logical lwrite,lexist
*
*	for f2c
	common /ccorrelatefield/ field1,field2,r,mean2,a,b,da,db,prob
*
        character line*256,datfile1*256,datfile2*256,
     +        title1*256,title2*256,vars(1)*5,lvars(1)*40,
     +        units1(1)*10,units2(1)*10
        integer iargc,llen
        external iargc,getarg,llen
*
*       check arguments
*       
        lwrite = .FALSE.
        n = iargc()
        if ( n.lt.3 ) then
            print *,'usage: correlatefieldfield '//
     +            'field1.[ctl|nc] field2.[ctl|nc] '//
     +            '[lag n[:m]] [sum|ave|max|min|sel n] '//
     +            '[log|sqrt|rank] '//
     +            '[minfac r] [minnum n] [begin yr] [end yr] '//
     +            '[lt cut] [gt cut] [diff [nyr]] [detrend] '//
     +            '[subtract] outfield'
            stop
        endif
        call getarg(1,line)
        call getmens(line,mens1,nensmax)
        if ( lwrite ) print *,'correlatefieldfield: nf_opening file '
     +        ,line(1:llen(line))
        status = nf_open(line,nf_nowrite,ncid1)
        if ( status.ne.nf_noerr ) then
            call parsectl(line,datfile1,nxmax,nx1,xx1,nymax,ny1,yy1
     +            ,nzmax,nz1,zz1,nt1,nper1,firstyr1,firstmo1,u1,endian1
     +            ,title1,1,nvars,vars,ivars1,lvars,units1)
            nz1 = max(1,ivars1(1,1))
            ncid1 = -1
        else
            datfile1 = line
            call parsenc(line,ncid1,nxmax,nx1,xx1,nymax,ny1,yy1
     +            ,nzmax,nz1,zz1,nt1,nper1,firstyr1,firstmo1,u1,title1,1
     +            ,nvars,vars,jvars1,lvars,units1)
        endif
*
        call getarg(2,line)
        call getmens(line,mens2,nensmax)
        if ( lwrite ) print *,'correlatefieldfield: nf_opening file '
     +        ,line(1:llen(line))
        status = nf_open(line,nf_nowrite,ncid2)
        if ( status.ne.nf_noerr ) then
            call parsectl(line,datfile2,nxmax,nx2,xx2,nymax,ny2,yy2
     +            ,nzmax,nz2,zz2,nt2,nper2,firstyr2,firstmo2,u2,endian2
     +            ,title2,1,nvars,vars,ivars2,lvars,units2)
            nz2 = max(1,ivars2(1,1))
            ncid2 = -1
        else
            datfile2 = line
            call parsenc(line,ncid2,nxmax,nx2,xx2,nymax,ny2,yy2
     +            ,nzmax,nz2,zz2,nt2,nper2,firstyr2,firstmo2,u2,title2,1
     +            ,nvars,vars,jvars2,lvars,units2)
        endif
*
        nxf = max(nx1,nx2)
        nyf = max(ny1,ny2)
        nzf = max(nz1,nz2)
        firstyr = max(firstyr1,firstyr2)
        nperyear = max(nper1,nper2)
        lastyr = min(firstyr1 + (nt1-1)/nperyear,firstyr2 + (nt2-1)
     +        /nperyear)
        if ( nxf*nyf*nzf*nperyear*(lastyr-firstyr+1).gt.nxmax*nymax*12
     +        *nyrmax*(1+max(mens1,mens2)) ) then
            write(0,*) 'correlatefieldfield: error: fields too large: '
     +            ,nxf,nyf,nzf,firstyr,yrend,nxmax,nymax,nyrmax
     +            ,1+max(mens1,mens2)
            write(*,*) 'correlatefieldfield: error: fields too large: '
     +            ,nxf,nyf,nzf,firstyr,yrend,nxmax,nymax,nyrmax
     +            ,1+max(mens1,mens2)
            call abort
        endif
        if ( nxf*nyf*nzf .gt. nxmax*nymax*nlevmax ) then
            write(0,*) 'correlatefieldfield: error: fields too large: '
     +            ,nxf,nyf,nzf,nxmax,nymax,nlevmax
            write(*,*) 'correlatefieldfield: error: fields too large: '
     +            ,nxf,nyf,nzf,nxmax,nymax,nlevmax
            call abort
        endif
        if ( nper1.ne.nper2 ) then
            write(0,*) 'correlatefield: error: cannot handle different'
     +            //' time scales yet',nper1,nper2
            write(*,*) 'correlatefield: error: cannot handle different'
     +            //' time scales yet',nper1,nper2
            call abort            
        endif
        if ( mens1.gt.0 .and. mens2.gt.0 .and. mens1.ne.mens2 ) then
            write(0,*) 'correlatefield: warning: found different number'
     +            //'of ensemble members: ',mens1,mens2
            write(*,*) 'correlatefield: warning: found different number'
     +            //'of ensemble members: ',mens1,mens2
        endif
*
        call cfieldfield(nxf,nyf,nzf,nperyear,firstyr,
     +        datfile1,ncid1,field1,nx1,xx1,ny1,yy1,nz1,zz1,nt1,nper1
     +        ,firstyr1,firstmo1,u1,endian1,jvars1,mens1,
     +        datfile2,ncid2,field2,nx2,xx2,ny2,yy2,nz2,zz2,nt2,nper2
     +        ,firstyr2,firstmo2,u2,endian2,jvars2,mens2,
     +        r,prob,a,b,da,db,mean2)
*
        end

        subroutine cfieldfield(nxf,nyf,nzf,nperyear,firstyr,
     +        datfile1,ncid1,field1,nx1,xx1,ny1,yy1,nz1,zz1,nt1,nper1
     +        ,firstyr1,firstmo1,u1,endian1,jvars1,mens1,
     +        datfile2,ncid2,field2,nx2,xx2,ny2,yy2,nz2,zz2,nt2,nper2
     +        ,firstyr2,firstmo2,u2,endian2,jvars2,mens2,
     +        r,prob,a,b,da,db,mean2)
*
*       break to use the arrays field1,field2,... compactly and conserve RAM
*
        implicit none
#include "params.h"
#include "recfac.h"
#include "getopts.inc"
*
        integer nxf,nyf,nzf,nperyear,firstyr,
     +        ncid1,nx1,ny1,nz1,nt1,nper1,firstyr1,firstmo1,endian1
     +        ,jvars1(6,nvmax),mens1,
     +        ncid2,nx2,ny2,nz2,nt2,nper2,firstyr2,firstmo2,endian2
     +        ,jvars2(6,nvmax),mens2
        real field1(nxf,nyf,nzf,nperyear,firstyr:yrend,0:mens1),
     +        field2(nxf,nyf,nzf,nperyear,firstyr:yrend,0:mens2),
     +        xx1(nx1),yy1(ny1),zz1(nz1),xx2(nx2),yy2(ny2),zz2(nz2),
     +        u1,u2,
     +        r(nxf,nyf,nzf,0:npermax),
     +        prob(nxf,nyf,nzf,0:npermax),
     +        a(nxf,nyf,nzf,0:npermax),
     +        b(nxf,nyf,nzf,0:npermax),
     +        da(nxf,nyf,nzf,0:npermax),
     +        db(nxf,nyf,nzf,0:npermax),
     +        mean2(nxf,nyf,nzf,0:npermax)
        character*(*) datfile1,datfile2
*       
        integer ntmax
        parameter(ntmax=1000)
        integer i,j,jj,k,n,lag,jx,jy,jz,mo,j1,j2,m,ii,l,yr,if,jm,jp,im
     +        ,ip,lastyr,nt,nrec,nvars,ivars(2,6),ldir,nx,ny,nz,ncid
     +        ,ntvarid,itimeaxis(ntmax)
        real  fxy1(npermax,yrbeg:yrend,0:nensmax),
     +        fxy2(npermax,yrbeg:yrend,0:nensmax),
     +        ddata(ndata),dindx(ndata),dddata(ndata),ddindx(ndata)
     +        ,adata,sxx,aindx,syy,sxy,df,d,zd,z,probd,absent,chi2,q,sum
     +        ,fac,xx(nxmax),yy(nxmax),zz(nzmax),aa(2),x(ndata)
     +        ,sig(ndata),u(ndata,2),v(2,2),w(2),daa(2,2),filter(100)
        common /c_findx/ dddata,dindx
*
*	not only for f2c
        common /c_cfieldfield/ fxy1,fxy2
        parameter (absent=3e33)
        logical lexist
        character line*80,yesno*1,string*10,file*255,outfile*255,
     +        datfile*255,title*255,vars(6)*10,lvars(6)*40,dir*255,
     +        units(6)*10
        integer iargc,llen
        external iargc,getarg,llen,findx
*       
*       check if arrays big enough
        if ( nperyear.gt.npermax ) then
            write(0,*)'correlatefieldfield: recompile with npermax '
     +            ,nperyear,npermax
            write(*,*)'correlatefieldfield: recompile with npermax '
     +            ,nperyear,npermax
            call abort
        endif
*
*       save time on the initialization - but not too much.
        lastyr = min(firstyr1 + (nt1-1)/nperyear,firstyr2 + (nt2-1)
     +        /nperyear)
        nt = nperyear*(lastyr-firstyr+1)
        n = iargc()
        call getopts(3,n-1,nperyear,yrbeg,yrend,.true.,mens?)
        if ( lag1.lt.0 ) print *,'(field1 leading field2)'
        if ( lag2.gt.0 ) print *,'(field2 leading field1)'
	if ( dump ) write(0,*)'correlatefieldfield: dump not supported'
	if ( plot ) write(0,*)'correlatefieldfield: plot not supported'
	if ( lks ) write(0,*)'correlatefieldfield: K-S not supported'
	if ( lconting ) write(0,*)'correlatefieldfield: contingency '//
     +        'tables not supported'
        if ( composite ) write(0,*)'composites not yet supported'
	do i=1,indxuse
	    if ( lincl(i) ) write(0,*)'correlatefieldfield: what do ',
     +            'you mean with ',strindx(i),'?'
	enddo
        if ( nens2.gt.max(mens1,mens2) ) then
            write(0,*) 'There are only enesemble members up to '
     +            ,max(mens1,mens2),', using that instead of ',nens2
            nens2 = max(mens1,mens2)
        endif
        mens1 = min(mens1,nens2)
        mens2 = min(mens2,nens2)
        yr1 = max(yr1,firstyr,firstyr - (min(lag1,lag2)+nperyear-1)
     +        /nperyear)
        yr2 = min(yr2,lastyr,lastyr - (max(lag1,lag2)-nperyear+1)
     +        /nperyear)
        if ( lwrite ) then
            print *,'cfieldfield: correlating ',datfile1(1:llen(datfile1
     +            )),llen(datfile1),len(datfile1)
            print *,'                    with ',datfile2(1:llen(datfile2
     +            )),llen(datfile2),len(datfile2)
            print *,'years: ',yr1,yr2
        endif
*
*       init
*
        if ( m1.ne.m2 .and. lag1.ne.lag2 ) then
            print *,'Sorry, can only handle either lags varying or'//
     +            ' months varying, not both'
            print *,'(months:',m1,m2,', lags:',lag1,lag2,')'
            call abort
        endif
        if ( lag2-lag1.gt.npermax ) then
            print *,'Sorry, can only store ',npermax+1,' fields maximum'
            call abort
        endif
        if ( lag1.ne.lag2 .and. lsubtract ) then
            print *,'Sorry, can only subtract one lag'
            call abort
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
	print *,'init'
**	if ( lwrite ) print *,'filling fields with absent'
**        call makefieldabsent(field1,nx1,ny1,yrbeg,yrend,yr1,yr2)
**        call makefieldabsent(field2,nx2,ny2,yrbeg,yrend,yr1,yr2)
        do iens=0,nens2
            call makeabsent(fxy1,npermax,yrbeg,yrend,iens)
            call makeabsent(fxy2,npermax,yrbeg,yrend,iens)
        enddo
*
*       compute minfac if it has not been set explicitly
*
        if ( minfac.lt.0 .and. minnum.lt.0 ) then
*           heuristic, gives 0.25 for 150 yrs, 0.5 for 50 yrs, 0.75 for 20yrs
            minfac = max(0.1,
     +            min(0.6,
     +            1.5-log(1+real(min(nt,nperyear*(yr2-yr1+1))-1)
     +            /nperyear)/4))
        endif
        write(0,'(a,i2,a)') 'Requiring at least ',
     +            nint(100*minfac),'% valid points<p>'
*       
*       compute derivative coefficients if these are needed
*       
        if ( nfittime.gt.0 ) then
*           the last value is a heuristic based on advise in NumRec
            call savgol(filter,100,nfittime,nfittime,1,min(nfittime,4))
        endif
*
*       read fields
*
        call keepalive(1,2)
        if ( mens1.ge.0 ) then
            do iens=nens1,nens2
                datfile = datfile1
                call filloutend(datfile,iens)
                if ( ncid1.eq.-1 ) then
                    call zreaddatfile(datfile,field1(1,1,1,1,firstyr
     +                    ,iens),nxf,nyf,nzf,nx1,ny1,nz1,nperyear
     +                    ,firstyr,yrend,firstyr1,firstmo1,nt1,u1
     +                    ,endian1,lwrite,yr1,yr2,1,1)
                else
                    write(0,*) 'Cannot handle ensembles of netCDF '//
     +                    'files (yet)'
                    call abort
                endif
            enddo
        else
            if ( ncid1.eq.-1 ) then
                call zreaddatfile(datfile1,field1,nxf,nyf,nzf,nx1,ny1
     +                ,nz1,nperyear,firstyr,yrend,firstyr1,firstmo1,nt1
     +                ,u1,endian1,lwrite,yr1,yr2,1,1)
            else
                if ( nzf.ne.1 ) then
                    write(0,*) 'cannot read 3D netCDF files yet'
                    call abort
                endif
                call readncfile(ncid1,field1,nxf,nyf,nx1,ny1,nperyear
     +                ,firstyr,yrend,firstyr1,firstmo1,nt1,u1,lwrite,yr1
     +                ,yr2,jvars1)
            endif
        endif
        if ( lwrite ) then
            print *,'field1 @ 0,60N'
            call dump060(xx1,yy1,zz1,field1,nxf,nyf,nzf,nx1,ny1,nz1
     +            ,nperyear,firstyr,yrend)
        endif
        call keepalive(2,2)
        if ( mens2.ge.0 ) then
            do iens=nens1,nens2
                datfile = datfile2
                call filloutend(datfile,iens)
                if ( ncid2.eq.-1 ) then
                    call zreaddatfile(datfile2,field2(1,1,1,1,firstyr
     +                    ,iens),nxf,nyf,nzf,nx2,ny2,nz2,nperyear
     +                    ,firstyr,yrend,firstyr2,firstmo2,nt2,u2
     +                    ,endian2,lwrite,yr1,yr2,1,1)
                else
                    write(0,*) 'correlatefieldfield: cannot handle '//
     +                    'ensembles of netCDF files (yet)'
                    call abort
                endif
            enddo
        else
            if ( ncid2.eq.-1 ) then
                call zreaddatfile(datfile2,field2,nxf,nyf,nzf,nx2,ny2
     +                ,nz2,nperyear,firstyr,yrend,firstyr2,firstmo2,nt2
     +                ,u2,endian2,lwrite,yr1,yr2,1,1)
            else
                if ( nzf.ne.1 ) then
                    write(0,*) 'correlatefieldfield: cannot handle 3D '/
     +                    /'netCDF files yet'
                    call abort
                endif
                call readncfile(ncid2,field2,nxf,nyf,nx2,ny2,nperyear
     +                ,firstyr,yrend,firstyr2,firstmo2,nt2,u2,lwrite,yr1
     +                ,yr2,jvars2)
            endif
        endif
        if ( lwrite ) then
            print *,'field1 @ 0,60N'
            call dump060(xx2,yy2,zz2,field2,nxf,nyf,nzf,nx2,ny2,nz2
     +            ,nperyear,firstyr,yrend)
        endif
        call keepalive(0,0)
*
*       interpolate fields to common grid
*
        do iens=0,max(0,mens1)
            if ( lwrite ) print *,'calling zinterpu'
            call zinterpu(
     +            field1(1,1,1,1,yr1,iens),zz1,nx1,ny1,nz1,
     +            field2(1,1,1,1,yr1,iens),zz2,nx2,ny2,nz2,
     +            zz,nz,yr1,yr2,yr1,yr2,nxf,nyf,nzf,nperyear,lwrite)
            if ( lwrite ) print *,'calling xyinterpu'
            call xyinterpu(
     +            field1(1,1,1,1,yr1,iens),xx1,nx1,yy1,ny1,
     +            field2(1,1,1,1,yr1,iens),xx2,nx2,yy2,ny2,
     +            xx,nx,yy,ny,yr1,yr2,yr1,yr2,nxf,nyf,nzf,nz,
     +            nperyear,intertype,lwrite)
        enddo
*
*       loop over grid points
*
	print *,'correlating'
        do jz=1,nz
            do jy=1,ny
                call keepalive(jy+(jz-1)*ny,ny*nz)
                do jx=1,nx
                    do mo=0,nperyear
                        r(jx,jy,jz,mo) = absent
                        prob(jx,jy,jz,mo) = absent
                        a(jx,jy,jz,mo) = absent
                        b(jx,jy,jz,mo) = absent
                        da(jx,jy,jz,mo) = absent
                        db(jx,jy,jz,mo) = absent
                    enddo
*
*                   create 1-D series from fields
*
                    n = 0
                    do iens=0,max(0,mens1)
                        do i=yr1,yr2
                            do j=1,nperyear
                                fxy1(j,i,iens) = field1(jx,jy,jz,j,i
     +                                ,iens)
                                if ( fxy1(j,i,iens).lt.0.9*absent ) n =
     +                                n + 1
                            enddo
                        enddo
                    enddo
                    if ( n.lt.3 ) then
                        if ( lwrite ) print '(a,4i5)'
     +                        ,'not enough valid points in field1 at '
     +                        ,jx,jy,jz,n
                        goto 800
                    endif
                    n = 0
                    do iens=0,max(0,mens2)
                        do i=yr1,yr2
                            do j=1,nperyear
                                fxy2(j,i,iens) = field2(jx,jy,jz,j,i
     +                                ,iens)
                                if ( fxy2(j,i,iens).lt.0.9*absent ) n =
     +                                n + 1
                            enddo
                        enddo
                    enddo
                    if ( n.lt.3 ) then
                        if ( lwrite ) print '(a,4i5)'
     +                        ,'not enough valid points in field2 at '
     +                        ,jx,jy,jz,n
                        goto 800
                    endif
                    do iens=0,max(0,mens1)
*
*                   take monthly anomalies
*
                        if ( mdiff.gt.0 ) then
                            call mdiffit(fxy1(1,yrbeg,iens),npermax
     +                            ,nperyear,yrbeg,yrend,mdiff)
                        endif
*
*                   log, sqrt
*
                        if ( logscale ) then
                            call takelog(fxy1(1,yrbeg,iens),npermax
     +                            ,nperyear,yrbeg,yrend)
                        endif
                        if ( sqrtscale ) then
                            call takesqrt(fxy1(1,yrbeg,iens),npermax
     +                            ,nperyear,yrbeg,yrend)
                        endif
*
*                   sum
*
                        if ( lsum.gt.1 ) then
                            call sumit(fxy1(1,yrbeg,iens),npermax
     +                            ,nperyear,yrbeg,yrend,lsum,'v')
                        endif
*
*                   detrend
*
                        if ( ldetrend ) then
                            if ( lwrite ) print *,'Detrending field'
                            if ( lag1.eq.0 .and. lag2.eq.0 .or. m1.eq.0
     +                            .or.lsel.eq.nperyear ) then
                                call detrend(fxy1(1,yr1,iens),npermax
     +                                ,nperyear,yr1,yr2,m1,m2,lsel)
                            else
                                call detrend(fxy1(1,yr1),npermax
     +                                ,nperyear,yr1,yr2,m1,m2,lsel)
                            endif
                        endif
*
*                   differentiate
*
                        if ( ndiff.gt.0 ) then
                            if ( lwrite ) print *,'Taking differences'
                            call diffit(fxy1(1,yrbeg,iens),npermax
     +                            ,nperyear,yrbeg,yrend,ndiff)
                        endif
*
*                   anomalies - necessary if we consider more than one month
*
                        if ( anom .or. (lsel.gt.1 .or. nfittime.gt.0 )
     +                        .and.ndiff.eq.0 ) then
                            if ( lwrite ) print *,'Taking anomalies'
                            call anomal(fxy1(1,yrbeg,iens),npermax
     +                            ,nperyear,yrbeg,yrend,yr1,yr2)
                        endif
                    enddo
*
*                   copy ensemble members so that there is the same
*                   number of valid ones at every time step
*       
                    call makeensfull(@@)
                    do iens=0,max(0,mens2)
                        if ( mdiff2.gt.0 ) then
                            call mdiffit(fxy2(1,yrbeg,iens),npermax
     +                            ,nperyear,yrbeg,yrend,mdiff2)
                        endif
                        if ( lsum2.gt.1 ) then
                            call sumit(fxy2(1,yrbeg,iens),npermax
     +                            ,nperyear,yrbeg,yrend,lsum2,'v')
                        endif
                        if ( ldetrend ) then
                            if ( lag1.eq.0 .and. lag2.eq.0 .or. m1.eq.0
     +                            .or.lsel.eq.nperyear ) then
                                call detrend(fxy2(1,yr1,iens),npermax
     +                                ,nperyear,yr1,yr2,m1,m2,lsel)
                            else
                                call detrend(fxy2(1,yr1,iens),npermax
     +                                ,nperyear,yr1,yr2,m1,m2,lsel)
                            endif
                        endif
                        if ( ndiff.gt.0 ) then
                            call diffit(fxy2(1,yrbeg,iens),npermax
     +                            ,nperyear,yrbeg,yrend,ndiff)
                        endif
                        if ( anom .or. (lsel.gt.1 .or. nfittime.gt.0 )
     +                        .and.ndiff.eq.0 ) then
                            call anomal(fxy2(1,yrbeg,iens),npermax
     +                            ,nperyear,yrbeg,yrend,yr1,yr2)
                        endif
                    enddo
*       
*                   correlate   !
*
                    do mo=m1,m2
                        if ( mo.eq.0 ) then
                            j1 = 1
                            j2 = nperyear
                        else
                            j1 = mo
                            j2 = mo + lsel - 1
                        endif
                        do lag=lag1,lag2
*
*                           fill linear arrays without absent values 
*                           and compute r
*
                            n = 0
                            mean2(jx,jy,jz,mo) = 0
                            do yr=yr1-1,yr2
                                do jj=j1,j2
                                    j = jj
                                    call normon(j,yr,i,nperyear)
                                    if ( i.lt.yr1 .or.i.gt.yr2 )
     +                                    goto 710
                                    m = j-lag
                                    call normon(m,i,ii,nperyear)
                                    if ( ii.lt.yr1 .or.ii.gt.yr2 )
     +                                    goto 710
                                    if ( i.eq.(yr1+yr2)/2 .and. 
     +                                    jx.eq.nx/2.and.jy.eq.ny/2 
     +                                    .and. jz.eq.1 ) then
                                        print '(a,i3,i5,i3,i5,a)'
     +                                        ,'Correlating months',
     +                                        j,i,m,ii
     +                                        ,' of field1 and field2'
                                    endif
                                    if (  fxy1(j,i).lt.1e33 .and. (
     +                                    (fxy1(j,i).le.maxdata) .eqv.
     +                                    (fxy1(j,i).ge.mindata) .eqv.
     +                                    (maxdata.ge.mindata) ) .and. 
     +                                    fxy2(m,ii).lt.1e33 .and. (
     +                                    (fxy2(m,ii).le.maxindx) .eqv.
     +                                    (fxy2(m,ii).ge.minindx) .eqv.
     +                                    (maxindx.ge.minindx) ) ) then
                                        n = n+1
                                        ddata(n) = fxy1(j,i)
                                        dindx(n) = fxy2(m,ii)
                                        if ( nfittime.gt.0 ) then
                                            dddata(n) = 0
                                            do if=1,nfittime
                                                jm = j-if
                                                call normon(jm,i,im,
     +                                                nperyear)
                                                jp = j+if
                                                call normon(jp,i,ip,
     +                                                nperyear)
                                                if (  ip.gt.yrend .or. 
     +                                                im.lt.yrbeg ) then
                                                    n = n - 1
                                                    goto 710
                                                endif
                                                if ( fxy1(jp,ip).gt.1e33
     +                                               .or.
     +                                               fxy1(jm,im).gt.1e33
     +                                                ) then
                                                    n = n - 1
                                                    goto 710
                                                endif
                                                dddata(n) = dddata(n) +
     +                                                filter(1+if)*(
     +                                                fxy1(jm,im)
     +                                                - fxy1(jp,ip))
                                            enddo ! ip=1,nfittime
                                        endif ! nfittime
                                        mean2(jx,jy,jz,mo) = 
     +                                    mean2(jx,jy,jz,mo) + fxy2(j,i)
                                    endif
  710                               continue
                                enddo
                            enddo
                            if ( lwrite ) then
                                if ( mo.eq.0 ) then
                                    print *,'Comparing n=',n
     +                                    ,' with minfac*N = ',minfac
     +                                    ,min(nt,nperyear*(yr2-yr1+1))
                                else
                                    print *,'Comparing n=',n
     +                                    ,' with minfac*N = ',minfac
     +                                    ,min(nt/nperyear,yr2-yr1+1)
     +                                    *lsel
                                endif
                            endif
                            if (  mo.eq.0 .and. n.lt.minfac*
     +                            min(nt,nperyear*(yr2-yr1+1)) .or. 
     +                            mo.ne.0 .and.n.lt.minfac*
     +                            min(nt/nperyear,yr2-yr1+1)*lsel .or.
     +                            n.lt.minnum ) then
                                if ( lwrite ) print '(a,3i5,2i3,a,2i6)'
     +                                ,'not enough valid points at ',
     +                                jx,jy,jz,mo,lag,': ',n,nt
                                goto 790
                            endif
                            mean2(jx,jy,jz,mo) = mean2(jx,jy,jz,mo)/n
                            if ( m1.ne.m2 ) then
                                m = mo-m1
                            else
                                m = lag2-lag
                            endif
                            if ( jx.eq.nx/2 .and. jy.eq.ny/2 .and.
     +                            jz.eq.1 ) then
                                print '(a,i2)','and storing it at ',m
                            endif
                            if ( nfittime.gt.0 ) then
                                if ( n.le.4 ) then
                                    print*,'not enough data for timefit'
                                    goto 800
                                endif
*       
*                               fit to data = a*d(data)/dt + b*indx
*       
                                do j=1,n
                                    x(j) = j
                                enddo
                                do j=1,n
                                    sig(j) = 1
                                enddo
                                call svdfit(x,ddata,sig,n,aa,2,u,v,w,
     +                                ndata,2,chi2,findx)
                                call svdvar(v,2,2,w,daa,2)
                                a(jx,jy,jz,m) = aa(1)
                                b(jx,jy,jz,m) = aa(2)
                                da(jx,jy,jz,m) = sqrt(daa(1,1)/chi2)
                                db(jx,jy,jz,m) = sqrt(daa(2,2)/chi2)
                                do j=1,n
                                    x(j) = aa(1)*dddata(j) + 
     +                                    aa(2)*dindx(j)
                                enddo
                                if ( mo.eq.0 ) then
                                    df = n/(max(lsum,lsum2) + decor)/
     +                                    real(max(1,1-ndiff)) - 3
                                else
                                    df = n/(1 + (max(lsum,lsum2)-1)
     +                                    /nperyear+decor/nperyear)
     +                                    /real(max(1,1-ndiff)) - 3
                                endif
                                call pearsncross(ddata,x,n,r(jx,jy,jz,m)
     +                               ,prob(jx,jy,jz,m),z,adata,sxx,aindx
     +                               ,syy,sxy,df,ncrossvalidate)
*
                            elseif ( lrank ) then
                                if ( mo.eq.0 ) then
                                    sum = max(lsum,lsum2) + decor
                                else
                                    sum = 1 + (max(lsum,lsum2)-1)
     +                                    /nperyear +decor/nperyear
                                endif
                                call spearx(dindx,ddata,n,ddindx,dddata,
     +                                d,zd,probd,r(jx,jy,jz,m),
     +                                prob(jx,jy,jz,m),sum,adata,sxx
     +                                ,aindx,syy)
                            else
                                if ( mo.eq.0 ) then
                                    df = n/(max(lsum,lsum2) + decor) - 2
                                else
                                    df = n/(1 + (max(lsum,lsum2)-1)
     +                                    /nperyear+ decor/nperyear) - 2
                                endif
                                if ( df.le.0 ) then
                                    if ( lwrite ) print '(a,f5.2,3i5)'
     +                                    ,'error: df <= 0: ',df,jx,jy
     +                                    ,jz
                                    goto 790
                                endif
                                call pearsncross(dindx,ddata,n,
     +                               r(jx,jy,jz,m),prob(jx,jy,jz,m),z
     +                               ,aindx,syy,adata,sxx,sxy,df
     +                               ,ncrossvalidate)
                                if ( sxx.eq.0 .or. syy.eq.0 ) then
                                    r(jx,jy,jz,m) = absent
                                    prob(jx,jy,jz,m) = absent
                                else
                                    call fit(dindx,ddata,n,sig,0,
     +                                    a(jx,jy,jz,m),b(jx,jy,jz,m),
     +                                    da(jx,jy,jz,m),db(jx,jy,jz,m),
     +                                    chi2,q)
                                endif
                            endif
                            if ( abs(prob(jx,jy,jz,m)).lt.0.9*absent )
     +                            then
                                prob(jx,jy,jz,m) = 1-prob(jx,jy,jz,m)
                            else
                                prob(jx,jy,jz,m) = absent
                            endif
                            if (lwrite) print '(a,3i5,2i3,a,i6,a,6f9.4)'
     +                            ,'point ',jx,jy,jz,mo,lag,' OK (',n
     +                            ,'): ',r(jx,jy,jz,m),prob(jx,jy,jz,m)
     +                            ,a(jx,jy,jz,m),da(jx,jy,jz,m),b(jx,jy
     +                            ,jz,m),db(jx,jy,jz,m)
  790                       continue ! valid point/month
                        enddo   ! lag
                    enddo       ! month
  800               continue    ! valid point
                enddo           ! nx
            enddo               ! ny
        enddo                   ! nz
        if ( index(outfile,'.ctl').ne.0 ) then
            i = index(outfile,'.ctl')
            datfile = outfile(:i-1)//'.dat'
            open(unit=2,file=datfile,form='unformatted',access='direct'
     +            ,recl=recfac*nx*ny*nz,err=920)
        endif
        if ( .not.lsubtract ) then
            call getenv('DIR',dir)
            ldir = llen(dir)
            if ( ldir.eq.0 ) ldir=1
            if ( dir(ldir:ldir).ne.'/' ) then
                ldir = ldir + 1
                dir(ldir:ldir) = '/'
            endif
            title = ' '
            n = 1
            do i=0,iargc()-1
                call getarg(i,line)
                if ( i.le.2 ) then
                    do j=len(line),1,-1
                        if ( line(j:j).eq.'/' ) goto 810
                    enddo
  810               continue
                    j = j + 1
                    title(n:) = line(j:)
                else
                    title(n:) = line
                endif
                n = min(llen(title) + 2,len(title))
            enddo
            nvars = 6
            vars(1) = 'corr'
            vars(2) = 'sign'
            if ( nfittime.eq.0 ) then
                vars(3) = 'intercept'
                vars(4) = 'slope'
                vars(5) = 'errorintercept'
                vars(6) = 'errorslope'
            else
                vars(3) = 'relaxation'
                vars(4) = 'relation'
                vars(5) = 'drelaxation'
                vars(6) = 'drelation'
            endif
            do i=1,nvars
                ivars(1,i) = nz
                ivars(2,i) = 99
            enddo
            lvars(1) = 'correlation'
            lvars(2) = 'significance'
            if ( nfittime.eq.0) then
                lvars(3) = 'intercept'
                lvars(4) = 'slope'
                lvars(5) = 'error on intercept'
                lvars(6) = 'error on slope'
            else
                lvars(3) = 'relaxation'
                lvars(4) = 'error on relaxation'
                lvars(5) = 'relation'
                lvars(6) = 'error on relation'
            endif                        
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
            if ( index(outfile,'.ctl').ne.0 ) then
                call writectl(outfile,datfile,nx,xx,ny,yy,nz,zz
     +                ,1+(m2-m1)+(lag2-lag1),nperyear,i,j,3e33,title
     +                ,nvars,vars,ivars,lvars,units)
            else
                units = ' '
                call writenc(outfile,ncid,ntvarid,itimeaxis,ntmax,nx,xx
     +                ,ny,yy,nz,zz,1+(m2-m1)+(lag2-lag1),nperyear,i,j
     +                ,3e33,title,nvars,vars,ivars,lvarsi,units,0,0)
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
                            print *,'writing records ',6*m+1,'-',6*m+6
     +                            ,' of fields ',m,' of size ',nx*ny*nz
     +                            *recfa4
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
                            print *,'there are ',i
     +                            ,' valid values in record ',m
     +                            ,' with mean value ',d/i
                        endif
                        write(2,rec=6*m+1)
     +                    (((r(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                        write(2,rec=6*m+2)
     +                    (((prob(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                        write(2,rec=6*m+3)
     +                    (((a(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                        write(2,rec=6*m+4)
     +                    (((b(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                        write(2,rec=6*m+5)
     +                    (((da(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                        write(2,rec=6*m+6)
     +                    (((db(jx,jy,jz,m),jx=1,nx),jy=1,ny),jz=1,nz)
                    else
*                       netCDF file
                        call writencslice(ncid,0,0,0,ivars(1,1)
     +                       ,r(1,1,1,m),nxf,nyf,nzf,nx,ny,nz,m+1,1)
                        call writencslice(ncid,0,0,0,ivars(1,2)
     +                       ,prob(1,1,1,m),nxf,nyf,nzf,nx,ny,nz,m+1,1)
                        call writencslice(ncid,0,0,0,ivars(1,3)
     +                       ,a(1,1,1,m),nxf,nyf,nzf,nx,ny,nz,m+1,1)
                        call writencslice(ncid,0,0,0,ivars(1,4)
     +                       ,b(1,1,1,m),nxf,nyf,nzf,nx,ny,nz,m+1,1)
                        call writencslice(ncid,0,0,0,ivars(1,5)
     +                       ,da(1,1,1,m),nxf,nyf,nzf,nx,ny,nz,m+1,1)
                        call writencslice(ncid,0,0,0,ivars(1,6)
     +                       ,db(1,1,1,m),nxf,nyf,nzf,nx,ny,nz,m+1,1)
                    endif
                enddo
            enddo
            if ( index(outfile,'.ctl').eq.0 ) then
                call nf_close(ncid)
            endif
        else                    ! subtract
*
*           subtract best fit from field1 and give this as output
*
            nrec = 0
            fac = 1
            if ( lsum.gt.1 .and. oper.eq.'+' ) fac = fac/lsum
            do yr=firstyr1,firstyr1 + (nt1-1)/nperyear
                if ( yr.lt.firstyr2 .or. yr.gt.firstyr2 + (nt2-1)
     +                /nperyear ) then
                    if ( lwrite ) print *,'Invalid field at ',yr
                    do mo=1,nperyear
                        do jz=1,nz
                            do jy=1,ny
                                do jx=1,nx
                                    field1(jx,jy,jz,mo,yr) = 3e33
                                enddo
                            enddo
                        enddo
                    enddo
                else            ! overlap with field2
                    do mo=1,nperyear
                        if ( mo.ge.m1 .and. mo.le.m2 ) then
                            m = mo-m1
                        elseif (  mo.gt.m2 .and. mo.lt.m2+lsel .or.
     +                            mo+nperyear.gt.m2 .and. mo.lt.m2+lsel
     +                            ) then
                            m = m2-m1
                        else
                            m = -1
                        endif
                        if ( lwrite ) print *,'Using m=',m,' at mo=',mo
                        do jz=1,nz
                            do jy=1,ny
                                do jx=1,nx
                                    if ( m.ge.0 ) then
                                        if ( r(jx,jy,jz,m).lt.1e33 .and.
     +                                        field1(jx,jy,jz,mo,yr)
     +                                        .lt.1e33 .and.
     +                                        field2(jx,jy,jz,mo,yr)
     +                                        .lt.1e33 ) then
                                            field1(jx,jy,jz,mo,yr) = 
     +                                        field1(jx,jy,jz,mo,yr) 
     +                                        - b(jx,jy,jz,m)*fac*
     +                                        (field2(jx,jy,jz,mo,yr))
***     +                                        - mean2(jx,jy,jz,mo))
                                        else ! invalid point
                                            field1(jx,jy,jz,mo,yr) =
     +                                            3e33
                                        endif
                                    else ! invalid point
                                        field1(jx,jy,jz,mo,yr) = 3e33
                                    endif
                                enddo ! jy
                            enddo ! jx
                        enddo   ! jz
                    enddo       ! m
                endif           ! year OK
                do mo=1,nperyear
                    nrec = nrec + 1
                    if ( lwrite ) print *,'Writing new field1 for ',yr
     +                    ,mo
                    write(2,rec=nrec) 
     +                    (((field1(jx,jy,jz,mo,yr),
     +                    jx=1,nx),jy=1,ny),jz=1,nz)
                enddo
            enddo               ! yr
            if ( index(outfile,'.ctl').ne.0 ) then
                title = 'Data of '//
     +                datfile1(1:index(datfile1,' ')-1)//
     +                ' with the effect of '//
     +                datfile2(1:index(datfile2,' ')-1)//
     +                ' subtracted'
                nvars = 1
                vars(1) = 'var'
                ivars(1,1) = 0
                ivars(2,1) = 99
                lvars(1) = 'unknown, should be fixed'
                zz(1) = 0
                call writectl(outfile,datfile,nx,xx,ny,yy,1,zz
     +                ,nt1,nperyear,firstyr1,1,3e33,title,
     +                nvars,vars,ivars,lvars,units)
            endif
        endif                   ! subtract
        close(2)
*
*       error messages
*
        goto 999
  903   print *,'error reading date from file ',line(1:index(line,' ')-1
     +        ),' at record ',k
        call abort
  920   print *,'error cannot open new correlations file '
     +        ,datfile(1:index(datfile,' ')-1),' with record length '
     +        ,recfac*nx*ny*nz,recfac,nx,ny,nz
        call abort
  999   continue
        end
        
        subroutine findx(xi,f,n)
*       
*       used by the multiple-parameter fitting routine (lfittime)
*       
        implicit none
#include "params.h"
        integer n
        real xi,f(n)
        real dddata(ndata),dindx(ndata)
        common /c_findx/ dddata,dindx
        integer i,j
*
        if ( n.ne.2 ) goto 901
        i = nint(xi)
        if ( abs(xi-i).gt.0.01 ) goto 902
        f(1) = dddata(i)
        f(2) = dindx(i)
        return
  901   print *,'findx: should be called with n=2, not ',n
        call abort
  902   print *,'findx: wrong input! ',xi
        call abort
        end

        subroutine dump060(xx,yy,zz,field,nxf,nyf,nzf,nx,ny,nz
     +            ,nperyear,firstyr,lastyr)
*       
*       dumps the field at 0,60N
*       
        implicit none
#include "params.h"
        integer nxf,nyf,nzf,nx,ny,nz,nperyear,firstyr,lastyr
        real xx(nx),yy(ny),zz(nz),
     +        field(nxf,nyf,nzf,nperyear,firstyr:lastyr)
        integer x1,x2,y1,y2,i,j,yr,mo
        real lon1,lat1,lon2,lat2,lon1c,lat1c,lon2c,lat2c,
     +        data(npermax,yrbeg:yrend)
        
        lon1 = 0
        lat1 = 60
        lon2 = 0
        lat2 = 60
        call getlonwindow(lon1,lon2,x1,x2,xx,nx,lon1c,lon2c)
        call getlatwindow(lat1,lat2,y1,y2,yy,ny,lat1c,lat2c)
        if ( lon1c.gt.1e33 .or. lat1c.ge.1e33 ) then
            x1 = 1
            lon1c = xx(1)
            y1 = 1
            lat1c = yy(1)
        endif
        print *,'cutting out longitude ',x1,x2,lon1c,lon2c
        print *,'cutting out latitude  ',y1,y2,lat1c,lat2c
        call makeabsent (data,npermax,yrbeg,yrend)
        do yr=max(yrbeg,firstyr),yrend
            do mo=1,nperyear
                data(mo,yr) = field(x1,y1,1,mo,yr)
            enddo
        enddo
        call printdatfile(6,data,npermax,nperyear,yrbeg,yrend)
        end
