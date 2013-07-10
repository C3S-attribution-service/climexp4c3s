        program addfields
*
*       written to add convective and large-scale precipitation fields
*       copied from correlatefiledfield
*
        implicit none
#include "params.h"
#include "recfac.h"
#include "netcdf.inc"
        integer nyrmax
        parameter(nyrmax=75)
        integer j,n,ncid1,ncid2,nx1,ny1,nz1,nt1,nper1,firstyr1,firstmo1,
     +        nx2,ny2,nz2,nt2,nper2,firstyr2,firstmo2,nvars,ivars(2,1)
     +        ,jvars1(6,nvmax),jvars2(6,nvmax),endian1,endian2,status
     +        ,nxf,nyf,nperyear,firstyr
        real  xx1(nxmax),yy1(nymax),zz1(nzmax),
     +        xx2(nxmax),yy2(nymax),zz2(nzmax),
     +        field1(nxmax,nymax,12,nyrmax),
     +        field2(nxmax,nymax,12,nyrmax),
     +        u1,u2
        logical lwrite
*
*	for f2c
	common /ccorrelatefield/ field1,field2
*
        character line*256,datfile1*256,datfile2*256,
     +        title1*256,title2*256,vars(1)*5,lvars(1)*40,
     +        units1(1)*10,units2(1)*10
        integer iargc,llen
*
*       check arguments
*       
        lwrite = .FALSE.
        n = iargc()
        if ( n.lt.3 ) then
            print *,'usage: correlatefieldfield '//
     +            'field1.[ctl|nc] field2.[ctl|nc] '//
     +            '[lag n[:m]] [sum|ave|max|min|sel n] '//
     +            '[log|sqrt|rank] [minfac r] [begin yr] [end yr] '//
     +            '[lt cut] [gt cut] [diff [nyr]] [detrend] '//
     +            '[subtract] outfield'
            stop
        endif
        call getarg(1,line)
        if ( lwrite ) print *,'correlatefieldfield: nf_opening file '
     +        ,line(1:llen(line))
        status = nf_open(line,nf_nowrite,ncid1)
        if ( status.ne.nf_noerr ) then
            call parsectl(line,datfile1,nxmax,nx1,xx1,nymax,ny1,yy1
     +            ,nzmax,nz1,zz1,nt1,nper1,firstyr1,firstmo1,u1,endian1
     +            ,title1,1,nvars,vars,ivars,lvars,units1)
            ncid1 = -1
        else
            datfile1 = line
            call parsenc(line,ncid1,nxmax,nx1,xx1,nymax,ny1,yy1
     +            ,nzmax,nz1,zz1,nt1,nper1,firstyr1,firstmo1,u1,title1,1
     +            ,nvars,vars,jvars1,lvars,units1)
        endif
*
        call getarg(2,line)
        if ( lwrite ) print *,'correlatefieldfield: nf_opening file '
     +        ,line(1:llen(line))
        status = nf_open(line,nf_nowrite,ncid2)
        if ( status.ne.nf_noerr ) then
            call parsectl(line,datfile2,nxmax,nx2,xx2,nymax,ny2,yy2
     +            ,nzmax,nz2,zz2,nt2,nper2,firstyr2,firstmo2,u2,endian2
     +            ,title2,1,nvars,vars,ivars,lvars,units2)
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
        firstyr = max(firstyr1,firstyr2)
        nperyear = max(nper1,nper2)
        if ( nxf*nyf*nperyear*(yrend-firstyr+1).gt.nxmax*nymax*12*nyrmax
     +        ) then
            write(0,*) 'correlatefieldfield: error: fields too large: '
     +            ,nxf,nyf,firstyr,yrend,nxmax,nymax,nyrmax
            write(*,*) 'correlatefieldfield: error: fields too large: '
     +            ,nxf,nyf,firstyr,yrend,nxmax,nymax,nyrmax
            call abort
        endif
        if ( nper1.ne.nper2 ) then
            write(0,*) 'correlatefield: error: cannot handle different'
     +            //' time scales yet',nper1,nper2
            write(*,*) 'correlatefield: error: cannot handle different'
     +            //' time scales yet',nper1,nper2
            call abort            
        endif
*
        call cfieldfield(nxf,nyf,nperyear,firstyr,
     +        datfile1,ncid1,field1,nx1,xx1,ny1,yy1,nt1,nper1,firstyr1
     +        ,firstmo1,u1,endian1,jvars1,
     +        datfile2,ncid2,field2,nx2,xx2,ny2,yy2,nt2,nper2,firstyr2
     +        ,firstmo2,u2,endian2,jvars2)
*
        end

        subroutine cfieldfield(nxf,nyf,nperyear,firstyr,
     +        datfile1,ncid1,field1,nx1,xx1,ny1,yy1,nt1,nper1,firstyr1
     +        ,firstmo1,u1,endian1,jvars1,
     +        datfile2,ncid2,field2,nx2,xx2,ny2,yy2,nt2,nper2,firstyr2
     +        ,firstmo2,u2,endian2,jvars2)
*
*       break to use the arrays field1,field2 compactly and conserve RAM
*
        implicit none
#include "params.h"
#include "recfac.h"
#include "getopts.inc"
*
        integer nxf,nyf,nperyear,firstyr,ncid1,nx1,ny1,nt1,nper1
     +        ,firstyr1,firstmo1,endian1,jvars1(6,nvmax),ncid2,nx2,ny2
     +        ,nt2,nper2,firstyr2,firstmo2,endian2,jvars2(6,nvmax)
        real field1(nxf,nyf,nperyear,firstyr:yrend),field2(nxf,nyf
     +        ,nperyear,firstyr:yrend),xx1(nx1),yy1(ny1),xx2(nx2)
     +        ,yy2(ny2),u1,u2
        character*(*) datfile1,datfile2
*
        integer i,j,jj,k,n,lag,jx,jy,mo,j1,j2,m,ii,l,yr
     +        ,lastyr,nt,nrec,nvars,ivars(2,6),ldir,nx,ny
        real fxy1(npermax,yrbeg:yrend),fxy2(npermax,yrbeg:yrend),
     +        r(nxmax,nymax,0:npermax),mean2(nxmax,nymax,0:npermax),
     +        a(nxmax,nymax,0:npermax),b(nxmax,nymax,0:npermax),
     +        da(nxmax,nymax,0:npermax),db(nxmax,nymax,0:npermax),
     +        prob(nxmax,nymax,0:npermax),ddata(npermax*(yrend-yrbeg+1
     +        )),dindx(npermax*(yrend-yrbeg+1)),adata,sxx,aindx,syy,df
     +        ,d,zd,z,probd,absent,sig(1),chi2,q,sum,fac,xx(nxmax)
     +        ,yy(nxmax),zz(1)
*
*	not only for f2c
        common /c_cfieldfield/ fxy1,fxy2,r,mean2,a,b,da,db,prob,ddata
     +        ,dindx
        parameter (absent=3e33)
        logical lexist
        character line*80,yesno*1,string*10,file*255,outfile*255,
     +        datfile*255,title*255,vars(6)*10,lvars(6)*40,dir*255
        integer iargc,llen
        external iargc,getarg,llen
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
        call getopts(3,n-1,nperyear,yrbeg,yrend,.true.,0,0)
        if ( lag1.lt.0 ) print *,'(field2 leading field1)'
        if ( lag2.gt.0 ) print *,'(field1 leading field2)'
	if ( dump ) print *,'correlatefieldfield: dump not supported'
	if ( plot ) print *,'correlatefieldfield: plot not supported'
	if ( lks ) print *,'correlatefieldfield: K-S not supported'
	if ( lconting ) print *,'correlatefieldfield: contingency '//
     +        'tables not supported'
	do i=1,indxuse
	    if ( lincl(i) ) print *,'correlatefieldfield: what do ',
     +            'you mean with ',strindx(i),'?'
	enddo
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
        call makeabsent(fxy1,npermax,yrbeg,yrend)
        call makeabsent(fxy2,npermax,yrbeg,yrend)
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
*       read fields
*
        call keepalive(1,2)
        if ( ncid1.eq.-1 ) then
            call readdatfile(datfile1,field1,nxf,nyf,nx1,ny1,nperyear
     +            ,firstyr,yrend,firstyr1,firstmo1,nt1,u1,endian1,lwrite
     +            ,yr1,yr2,1,1)
        else
            call readncfile(ncid1,field1,nxf,nyf,nx1,ny1,nperyear
     +            ,firstyr,yrend,firstyr1,firstmo1,nt1,u1,lwrite,yr1,yr2
     +            ,jvars1)
        endif
        call keepalive(2,2)
        if ( ncid2.eq.-1 ) then
            call readdatfile(datfile2,field2,nxf,nyf,nx2,ny2,nperyear
     +            ,firstyr,yrend,firstyr2,firstmo2,nt2,u2,endian2,lwrite
     +            ,yr1,yr2,1,1)
        else
            call readncfile(ncid2,field2,nxf,nyf,nx2,ny2,nperyear
     +            ,firstyr,yrend,firstyr2,firstmo2,nt2,u2,lwrite,yr1,yr2
     +            ,jvars2)
        endif
        call keepalive(0,0)
*
*       interpolate fields to common grid
*
        call interpu(
     +        field1(1,1,1,yr1),xx1,yy1,nx1,ny1,
     +        field2(1,1,1,yr1),xx2,nx2,yy2,ny2,
     +        xx,nx,yy,ny,yr1,yr2,yr1,yr2,nxf,nyf,nperyear,intertype,
     +        lwrite)
*
*       add fields
*
        do i=yr1,yr2
            do j=1,12
                do jy=1,jy
                    do jx=1,nx
                        field1(jx,jy,j,i) = field1(jx,jy,j,i) +
     +                        field2(jx,jy,j,i)
                    enddo
                enddo
            enddo
        enddo
*       
*       output
*
        open(unit=2,file=datfile,form='unformatted',access='direct'
     +            ,recl=recfac*nx*ny,err=920)

        call writectl(outfile,datfile,nx,xx,ny,yy,1,zz
     +                ,1+(m2-m1)+(lag2-lag1),nperyear,i,j,3e33,title,
     +                nvars,vars,ivars,lvars,units)

        end
