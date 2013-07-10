        program svd
*
*       program make an SVD of two fields.
*       this version only accepts a single region and interpolates
*       the fields to a common grid, in orer to be able to output
*       a single outfile with the two fields.
*
        implicit none
#include "params.h"
#include "getopts.inc"
#include "recfac.h"
#include "netcdf.inc"
        integer nyrmax,nlevmax
        parameter(nyrmax=40,nlevmax=1)
        integer j,n,ncid1,ncid2,nx1,ny1,nz1,nt1,nper1,firstyr1,firstmo1,
     +        nx2,ny2,nz2,nt2,nper2,firstyr2,firstmo2,nvars,ivars1(2,1)
     +        ,ivars2(2,1),jvars1(6,nvmax),jvars2(6,nvmax),endian1
     +        ,endian2,status,nxf,nyf,nzf,nperyear,firstyr,lastyr,ncid
     +        ,mens
        real  xx1(nxmax),yy1(nymax),zz1(nzmax),
     +        xx2(nxmax),yy2(nymax),zz2(nzmax),
     +        field1(nxmax,nymax,12,nyrmax),
     +        field2(nxmax,nymax,12,nyrmax),
     +        u1,u2
*
*	for f2c
	common /ccorrelatefield/ field1,field2
*
        character line*256,datfile1*256,datfile2*256,
     +        title1*256,title2*256,vars(1)*5,lvars(1)*40,units(1)*10
        logical lexist
        integer iargc,llen
*
*       check arguments
*       
        lwrite = .FALSE.
        n = iargc()
        if ( n.lt.3 ) then
            print *,'usage: svd '//
     +            'field1.[ctl|nc] field2.[ctl|nc] '//
     +            '[lag n[:m]] [sum|ave|max|min|sel n] '//
     +            '[log|sqrt|rank] '//
     +            '[minfac r] [minnum n] [begin yr] [end yr] '//
     +            '[lt cut] [gt cut] [diff [nyr]] [detrend] '//
     +            '[subtract] outfield'
            stop
        endif
        call getarg(1,line)
        if ( index(line,'%').gt.0 .or. index(line,'++').gt.0 ) then
            ensemble = .true.
            call filloutens(line,0)
        else
            ensemble = .false.
            mens = 0
        endif
        if ( lwrite ) print *,'svd: nf_opening file '
     +        ,line(1:llen(line))
        status = nf_open(line,nf_nowrite,ncid1)
        if ( status.ne.nf_noerr ) then
            call parsectl(line,datfile1,nxmax,nx1,xx1,nymax,ny1,yy1
     +            ,nzmax,nz1,zz1,nt1,nper1,firstyr1,firstmo1,u1,endian1
     +            ,title1,1,nvars,vars,ivars1,lvars,units)
            nz1 = max(1,ivars1(1,1))
            ncid1 = -1
            if ( ensemble ) then
                do mens=1,nensmax
                    call getarg(1,line)
                    call filloutens(line,mens)
                    inquire(file=line,exist=lexist)
                    if ( .not.lexist ) goto 100
                enddo
  100           continue
                mens = mens - 1
                write(0,*) 'located ',mens+1,' ensemble members<br>'
            endif
        else
            datfile1 = line
            call parsenc(line,ncid1,nxmax,nx1,xx1,nymax,ny1,yy1
     +            ,nzmax,nz1,zz1,nt1,nper1,firstyr1,firstmo1,u1,title1,1
     +            ,nvars,vars,jvars1,lvars,units)
            if ( ensemble ) then
                do mens=1,nensmax
                    call getarg(1,line)
                    call filloutens(line,mens)
                    status = nf_open(line,nf_nowrite,ncid)
                    if ( status.ne.nf_noerr ) goto 200
                enddo
  200           continue
                mens = mens - 1
                write(0,*) 'located ',mens+1,' ensemble members<br>'
            endif
        endif
*
        call getarg(2,line)
        if ( ensemble ) then
           call filloutens(line,0)
        endif
        if ( lwrite ) print *,'svd: nf_opening file '
     +        ,line(1:llen(line))
        status = nf_open(line,nf_nowrite,ncid2)
        if ( status.ne.nf_noerr ) then
            call parsectl(line,datfile2,nxmax,nx2,xx2,nymax,ny2,yy2
     +            ,nzmax,nz2,zz2,nt2,nper2,firstyr2,firstmo2,u2,endian2
     +            ,title2,1,nvars,vars,ivars2,lvars,units)
            nz2 = max(1,ivars2(1,1))
            ncid2 = -1
        else
            datfile2 = line
            call parsenc(line,ncid2,nxmax,nx2,xx2,nymax,ny2,yy2
     +            ,nzmax,nz2,zz2,nt2,nper2,firstyr2,firstmo2,u2,title2,1
     +            ,nvars,vars,jvars2,lvars,units)
        endif
*
        nxf = max(nx1,nx2)
        nyf = max(ny1,ny2)
        nzf = max(nz1,nz2)
        firstyr = max(firstyr1,firstyr2)
        nperyear = max(nper1,nper2)
        lastyr = min(firstyr1 + (nt1-1)/nperyear,firstyr2 + (nt2-1)
     +        /nperyear)
        if ( nxf*nyf*nzf*nperyear*(lastyr-firstyr+1)*(mens+1).gt.
     +       nxmax*nymax*12*nyrmax) then
            write(0,*) 'svd: error: fields too large: ',
     +            nxf*nyf*nzf*nperyear*(lastyr-firstyr+1)*(mens+1),
     +            nxmax*nymax*12*nyrmax
***     +            ,nxf,nyf,nzf,firstyr,yrend,nxmax,nymax,nyrmax
            write(*,*) 'svd: error: fields too large: '
     +           ,nxf,nyf,nzf,firstyr,lastyr,mens,nxmax,nymax,nyrmax
            call abort
        endif
*
        call csvd(mens,nxf,nyf,nzf,nperyear,firstyr,lastyr,
     +        datfile1,ncid1,field1,nx1,xx1,ny1,yy1,nz1,zz1,nt1,nper1
     +        ,firstyr1,firstmo1,u1,endian1,jvars1,
     +        datfile2,ncid2,field2,nx2,xx2,ny2,yy2,nz2,zz2,nt2,nper2
     +        ,firstyr2,firstmo2,u2,endian2,jvars2)
*
        end

        subroutine csvd(mens,nxf,nyf,nzf,nperyear,firstyr,lastyr,
     +     datfile1,ncid1,field1,nx1,xx1,ny1,yy1,nz1,zz1,nt1,nper1
     +     ,firstyr1,firstmo1,u1,endian1,jvars1,
     +     datfile2,ncid2,field2,nx2,xx2,ny2,yy2,nz2,zz2,nt2,nper2
     +     ,firstyr2,firstmo2,u2,endian2,jvars2)
*
*       break to use the arrays field1,field2,... compactly and conserve RAM
*
        implicit none
#include "params.h"
#include "recfac.h"
#include "getopts.inc"
#include "netcdf.inc"
*
        integer mens,nxf,nyf,nzf,nperyear,firstyr,lastyr,
     +        ncid1,nx1,ny1,nz1,nt1,nper1,firstyr1,firstmo1,endian1,
     +        jvars1(6,nvmax),
     +        ncid2,nx2,ny2,nz2,nt2,nper2,firstyr2,firstmo2,endian2,
     +        jvars2(6,nvmax)
        real field1(nxf,nyf,nzf,nperyear,firstyr:lastyr,0:mens),
     +       field2(nxf,nyf,nzf,nperyear,firstyr:lastyr,0:mens),
     +       xx1(nx1),yy1(ny1),zz1(nz1),xx2(nx2),yy2(ny2),zz2(nz2),
     +       u1,u2
        character*(*) datfile1,datfile2
*       
        integer nxymax,lwork,liwork
        parameter(nxymax=120*51,lwork=107208721,liwork=36603)
        integer i,j,jj,k,m,n,jx,jy,jz,jx1,jy1,jz1,jx2,jy2,jz2,mo,i1,i2
     +       ,j1,j2,n1,n2,nn1,nn2,ii,l,yr,nt,nrec,nvars,ivars(2,7),ldir
     +       ,nx,ny,nz,ncid,ntvarid,iens,status,x1,x2,y1,y2
        integer ijk1(3,nxymax),ijk2(3,nxymax)
        real xx(nxmax),yy(nxmax),zz(nzmax),wx(nxmax),wy(nymax),
     +       tarray(2),s,f1,f2
        real cov(nxymax,nxymax)
        logical lexist,xrev,yrev,xwrap
        character line*80,yesno*1,string*10,file*255,outfile*255,
     +       datfile*255,title*255,vars(7)*10,lvars(7)*40,dir*255,
     +       units(7)*10
        integer iargc,llen
        real etime
        external llen,findx
*       
*       check if arrays big enough
        if ( nperyear.gt.npermax ) then
            write(0,*)'svd: recompile with npermax '
     +           ,nperyear,npermax
            write(*,*)'svd: recompile with npermax '
     +           ,nperyear,npermax
            call abort
        endif
*
*       save time on the initialization - but not too much.
        nt = nperyear*(lastyr-firstyr+1)
        n = iargc()
        call getopts(3,n-1,nperyear,yrbeg,yrend,.true.,0,0)
        if ( lag1.lt.0 ) print *,'(field1 leading field2)'
        if ( lag2.gt.0 ) print *,'(field2 leading field1)'
        if ( dump ) write(0,*)'svd: dump not supported'
        if ( plot ) write(0,*)'svd: plot not supported'
        if ( lks ) write(0,*)'svd: K-S not supported'
        if ( lconting ) write(0,*)'svd: contingency '//
     +       'tables not supported'
        if ( composite ) write(0,*)'composites not yet supported'
        do i=1,indxuse
            if ( lincl(i) ) write(0,*)'svd: what do ',
     +           'you mean with ',strindx(i),'?'
        enddo
        yr1 = max(yr1,firstyr,firstyr - (min(lag1,lag2)+nperyear-1)
     +       /nperyear)
        yr2 = min(yr2,lastyr,lastyr - (max(lag1,lag2)-nperyear+1)
     +       /nperyear)
        if ( lwrite ) then
            print *,'csvd: correlating ',datfile1(1:llen(datfile1
     +           )),llen(datfile1),len(datfile1)
            print *,'                    with ',datfile2(1:llen(datfile2
     +           )),llen(datfile2),len(datfile2)
            print *,'years: ',yr1,yr2
        endif
*
*       init
*
        if ( m1.ne.m2 .and. lag1.ne.lag2 ) then
            print *,'Sorry, can only handle either lags varying or'//
     +           ' months varying, not both'
            print *,'(months:',m1,m2,', lags:',lag1,lag2,')'
            call abort
        endif
        if ( lag2.ne.lag1 ) then
            print *,'Sorry, only makes sense for lag1=lag2'
            call abort
        endif
        call getarg(n,outfile)
        inquire(file=outfile,exist=lexist)
        if ( lexist ) then
            print *,'output file ',outfile(1:index(outfile,' ')-1),
     +           ' already exists, overwrite? [y/n]'
            read(*,'(a)') yesno
            if (  yesno.ne.'y' .and. yesno.ne.'Y' .and. 
     +           yesno.ne.'j' .and. yesno.ne.'J' ) then
                stop
            endif
            open(1,file=outfile)
            close(1,status='delete')
        endif
        print *,'init'
*       
*       compute minfac if it has not been set explicitly
*       
        if ( minfac.lt.0 .and. minnum.lt.0 ) then
*           heuristic, gives 0.25 for 150 yrs, 0.5 for 50 yrs, 0.75 for
*           20yrs
            minfac = max(0.1,
     +           min(0.6,
     +           1.5-log(1+real(min(nt,nperyear*(yr2-yr1+1))-1)
     +           /nperyear)/4))
        endif
        write(0,'(a,i2,a)') 'Requiring at least ',
     +       nint(100*minfac),'% valid points<p>'
*       
*       read fields
*       
        call keepalive(1,2)
        if ( ncid1.eq.-1 ) then
            if ( ensemble ) then
                write(0,*) 'Using ensemble members ',nens1,' to ',
     +               nens2,'<p>'
                call getarg(1,dir)
            endif
            do iens=nens1,nens2
                call keepalive(iens-nens1,2*(nens2-nens1+1))
                if ( iens.gt.nens1 ) then
                    outfile=dir
                    call filloutens(outfile,iens)
                    open(1,file=outfile,status='old')
                    read(1,'(a)') line
                    j = index(line,'^')
                    if ( j.ne.0 ) then
                        do i=len(outfile),1,-1
                            if ( outfile(i:i).eq.'/' ) goto 110
                        enddo
  110                   continue
                        datfile1 = outfile(1:i)//line(j+1:)
                    else
                        datfile1 = line(6:)
                    endif
                endif
                call zreaddatfile(datfile1,field1(1,1,1,1,firstyr,iens)
     +               ,nxf,nyf,nzf,nx1,ny1,nz1,nperyear,firstyr,lastyr
     +               ,firstyr1,firstmo1,nt1,u1,endian1,lwrite,yr1,yr2,1
     +               ,1)
            enddo
        else
            if ( nzf.ne.1 ) then
                write(0,*) 'cannot read 3D netCDF files yet'
                call abort
            endif
            if ( ensemble ) then
                write(0,*) 'Using ensemble members ',nens1,' to ',nens2
     +               ,'<p>'
            endif
            do iens=nens1,nens2
                call keepalive(iens-nens1,2*(nens2-nens1+1))
                if ( ensemble ) then
                    call getarg(1,file)
                    call filloutens(file,iens)
                    status = nf_open(file,nf_nowrite,ncid)
                    call parsenc(file,ncid1,nxmax,nx1,xx1,nymax,ny1
     +                   ,yy1,nzmax,nz1,zz1,nt1,nperyear,firstyr1
     +                   ,firstmo1,u1,title,1,nvars,vars,jvars1,lvars
     +                   ,units)
                endif
                call readncfile(ncid1,field1(1,1,1,1,firstyr,iens)
     +               ,nxf,nyf,nx1,ny1,nperyear,firstyr,lastyr,firstyr1
     +               ,firstmo1,nt1,u1,lwrite,yr1,yr2,jvars1)
            enddo
        endif
        if ( lwrite ) then
            print *,'field1 @ 0,60N'
            call dump060(xx1,yy1,zz1,field1,nxf,nyf,nzf,nx1,ny1,nz1
     +            ,nperyear,firstyr,lastyr)
        endif
        if ( ncid2.eq.-1 ) then
            if ( ensemble ) then
                call getarg(2,dir)
            endif
            do iens=nens1,nens2
                call keepalive(iens-nens1+nens2-nens1+1,
     +               2*(nens2-nens1+1))
                if ( iens.gt.nens1 ) then
                    file=dir
                    call filloutens(file,iens)
                    open(1,file=file,status='old')
                    read(1,'(a)') line
                    j = index(line,'^')
                    if ( j.ne.0 ) then
                        do i=len(file),1,-1
                            if ( file(i:i).eq.'/' ) goto 120
                        enddo
  120                   continue
                        datfile2 = file(1:i)//line(j+1:)
                    else
                        datfile2 = line(6:)
                    endif
                endif
                call zreaddatfile(datfile2,field2(1,1,1,1,firstyr,iens)
     +               ,nxf,nyf,nzf,nx2,ny2,nz2,nperyear,firstyr,lastyr
     +               ,firstyr2,firstmo2,nt2,u2,endian2,lwrite,yr1,yr2,1
     +               ,1)
            enddo
        else
            if ( nzf.ne.1 ) then
                write(0,*) 'cannot handle 3D netCDF files yet'
                call abort
            endif
            do iens=nens1,nens2
                call keepalive(iens-nens1,2*(nens2-nens1+1))
                if ( ensemble ) then
                    call getarg(2,file)
                    call filloutens(file,iens)
                    status = nf_open(file,nf_nowrite,ncid)
                    call parsenc(file,ncid2,nxmax,nx2,xx2,nymax,ny2
     +                   ,yy2,nzmax,nz2,zz2,nt1,nperyear,firstyr2
     +                   ,firstmo2,u2,title,1,nvars,vars,jvars2,lvars
     +                   ,units)
                endif
                call readncfile(ncid2,field2(1,1,1,1,firstyr,iens)
     +               ,nxf,nyf,nx2,ny2,nperyear,firstyr,lastyr,firstyr2
     +               ,firstmo2,nt2,u2,lwrite,yr1,yr2,jvars2)
            enddo
        endif
        if ( lwrite ) then
            print *,'field2 @ 0,60N'
            call dump060(xx2,yy2,zz2,field2,nxf,nyf,nzf,nx2,ny2,nz2
     +            ,nperyear,firstyr,lastyr)
        endif
        call keepalive(0,0)
*
*       interpolate fields to common grid
*
        if ( lwrite ) print *,'calling zinterpu'
        do iens=0,mens
           call zinterpu(
     +          field1(1,1,1,1,yr1,iens),zz1,nx1,ny1,nz1,
     +          field2(1,1,1,1,yr1,iens),zz2,nx2,ny2,nz2,
     +          zz,nz,yr1,yr2,yr1,yr2,nxf,nyf,nzf,nperyear,lwrite)
        enddo
        if ( lwrite ) print *,'calling xyinterpu'
        call ensxyinterpu(
     +       field1(1,1,1,1,yr1,0),xx1,nx1,yy1,ny1,
     +       field2(1,1,1,1,yr1,0),xx2,nx2,yy2,ny2,
     +       xx,nx,yy,ny,yr1,yr2,yr1,yr2,nxf,nyf,nzf,mens,nz,
     +       nperyear,intertype,lwrite)
*
*	get boundaries in grid points
*
        call getxyprop(xx,nx,yy,ny,xrev,yrev,xwrap)
        call getlatlonwindow(lat1,lat2,lon1,lon2,xx,nx,xwrap,avex,yy,ny
     +        ,avey,x1,x2,y1,y2,lwrite)
*
*       average, cut out window - everything to make the arrays smaller
*       
        write(0,'(a,f8.2,a)') 'Averaging, shifting, cutting, time: '
     +        ,etime(tarray),'<br>'
*       save a copy of xx,yy
        nx1 = nx
        do i=1,nx           
           xx1(i) = xx(i)
        enddo
        ny1 = ny
        do i=1,ny
           yy1(i) = yy(i)
        enddo
        call enscutoutwindow(x1,x2,y1,y2,xx1,nx1,xwrap,xrev,avex,yy1,ny1
     +       ,avey,wx,wy,field1,nxf,nyf,nens1,nens2,nperyear,
     +       firstyr,lastyr,yr1,yr2,lwrite)
        call enscutoutwindow(x1,x2,y1,y2,xx,nx,xwrap,xrev,avex,yy,ny
     +       ,avey,wx,wy,field2,nxf,nyf,nens1,nens2,nperyear,
     +       firstyr,lastyr,yr1,yr2,lwrite)
*
*       compute covariance matrix
*
       print '(a,f8.2)','Processing timeseries options, time: '
     +        ,etime(tarray)
        call processoptions(field1,nx,ny,nz,nperyear,firstyr,lastyr,mens
     +       ,1)
        call processoptions(field2,nx,ny,nz,nperyear,firstyr,lastyr,mens
     +       ,2)
*
*       covariance
*
        do n1=1,nxymax
            do n2=1,nxymax
                cov(n1,n2) = 3e33
            enddo
        enddo
        n1 = 0
        n2 = 0
        do jx1=1,nx
            do jy1=1,ny
                do jz1=1,jz
                    do jx2=1,nx
                        do jy2=1,ny
                            do jz2=1,jz
                                s = 0
                                n = 0
                                mo=m1
                                do iens=nens1,nens2
                                    do yr=yr1,yr2
                                        if ( mo.eq.0 ) then
                                            j1 = 1
                                            j2 = nperyear
                                        else
                                            j1 = mo
                                            j2 = mo + lsel - 1
                                        endif
                                        do jj=j1,j2
                                            j = jj
                                            call normon(j,yr,i,nperyear)
                                            if ( i.lt.firstyr .or.
     +                                           i.gt.lastyr ) 
     +                                           goto 710
                                            m = j-lag1
                                            call normon(m,i,ii,nperyear)
                                            if ( ii.lt.firstyr .or.
     +                                           ii.gt.lastyr ) 
     +                                           goto 710
                                            f1 = field1(jx1,jy1,jz1,j,i
     +                                           ,iens)
                                            f2 = field1(jx2,jy2,jz2,m,ii
     +                                           ,iens)
                                            if ( f1.lt.1e33 .and. (
     +                                           ((f1.le.maxdata).eqv.
     +                                           (f1.ge.mindata)) .eqv.
     +                                           (maxdata.ge.mindata) )
     +                                           .and.
     +                                           f2.lt.1e33 .and. (
     +                                           ((f2.le.maxindx).eqv.
     +                                           (f2.ge.minindx) .eqv.
     +                                           (maxindx.ge.minindx) )
     +                                           ) ) then
                                                n = n + 1
                                                s = s + f1*f2
                                            endif
 710                                        continue
                                        enddo ! jj
                                    enddo ! yr
                                enddo ! iens
                                if ( lwrite ) then
                                    if ( mo.eq.0 ) then
                                        print *,'Comparing n=',n
     +                                       ,' with minfac*N = ',minfac
     +                                       ,min(nt,nperyear*(yr2-yr1+1
     +                                       ))
                                    else
                                        print *,'Comparing n=',n
     +                                       ,' with minfac*N = ',minfac
     +                                       ,min(nt/nperyear,yr2-yr1+1)
     +                                       *lsel
                                    endif
                                endif
                                if (  mo.eq.0 .and. n.lt.minfac*
     +                               min(nt,nperyear*(yr2-yr1+1)) .or. 
     +                               mo.ne.0 .and.n.lt.minfac*
     +                               min(nt/nperyear,yr2-yr1+1)*lsel .or
     +                               .n.lt.minnum ) then
                                    if ( lwrite ) print
     +                                   '(a,3i5,2i3,a,2i6)'
     +                                   ,'not enough valid points at '
     +                                   ,jx1,jy1,jz1,jx2,jy2,jz2,': ',n
     +                                   ,nt
                                    goto 790
                                endif
                                if ( n1.eq.0 ) then
                                    n1 = 1
                                    ijk1(1,n1) = jx1
                                    ijk1(2,n1) = jy1
                                    ijk1(3,n1) = jz1
                                elseif ( ijk1(1,n1).ne.jx1 .or.
     +                               ijk1(2,n1).ne.jy1 .or.
     +                               ijk1(3,n1).ne.jz1 ) then
                                    n1 = n1 + 1
                                    if ( n1.gt.nxymax ) then
                                        write(0,*) 'svd: error: too '//
     +                                       'many points in field1'
                                        write(*,*) 'svd: error: too '//
     +                                       'many points in field1'
                                        call abort
                                    endif
                                    ijk1(1,n1) = jx1
                                    ijk1(2,n1) = jy1
                                    ijk1(3,n1) = jz1
                                endif
                                if ( n2.eq.0 ) then
                                    n2 = 1
                                    ijk2(1,n2) = jx2
                                    ijk2(2,n2) = jy2
                                    ijk2(3,n2) = jz2
                                else
                                    do i2=n2,1,-1
                                        if ( ijk2(1,i).eq.jx2 .and.
     +                                       ijk2(2,i).eq.jy2 .and.
     +                                       ijk2(3,i).eq.jz2 ) goto 780
                                    enddo
                                    n2 = n2 + 1
                                    if ( n2.gt.nxymax ) then
                                        write(0,*) 'svd: error: too '//
     +                                       'many points in field2'
                                        write(*,*) 'svd: error: too '//
     +                                       'many points in field2'
                                        call abort
                                    endif
                                    ijk2(1,n2) = jx2
                                    ijk2(2,n2) = jy2
                                    ijk2(3,n2) = jz2
                                    i2 = n2
 780                                continue
                                endif
                                cov(n1,i2) = s/(n-1)
                                if (lwrite) print
     +                               '(a,3i5,2i3,a,i6,a,6f9.4)','point '
     +                               ,jx1,jy1,jz1,jx2,jy2,jz2,' OK (',n,
     +                               '): cov(',n1,i2,') = ',cov(n1,i2)
 790                            continue ! valid point
                            enddo ! nz2
                        enddo   ! ny2
                    enddo       ! nx2
                enddo           ! nz1
            enddo               ! ny1
        enddo                   ! nx1
*
*       I now have a covariance matrix with invalid data
*       get rid of the rows and columns with too many undefines
*       and put the rest of the undefineds to zero (I do not know a
*       better value)
*
 801    continue
        do i1=1,n1
            do i2=1,n2
                if ( cov(i1,i2).gt.1e30 ) then
                    nn1 = 0
                    do j1=1,n1
                        if ( cov(j1,i2).gt.1e33 ) nn1 = nn1 + 1
                    enddo
                    nn2 = 0
                    do j2=1,n2
                        if ( cov(j1,i2).gt.1e33 ) nn2 = nn2 + 1
                    enddo
                    if ( nn1.gt.n1/2 .or. nn2.gt.n2/2 ) then
                        if ( nn1/real(n1).gt.nn2/real(n2) ) then
                            if ( lwrite ) print *,'deleting row ',i1
     +                           ,' with ',nn1,'/',n1,' undefs'
                            n1 = n1 - 1
                            do j1=i1,n1
                                do k=1,3
                                    ijk1(k,j1) = ijk1(k,j1+1)
                                enddo
                                do j2=1,n2
                                    cov(j1,j2) = cov(j1+1,j2)
                                enddo
                            enddo
                        else
                            if ( lwrite ) print *,'deleting col ',i2
     +                           ,' with ',nn2,'/',n2,' undefs'
                            n2 = n2 - 1
                            do j2=i2,n2
                                do k=1,3
                                    ijk2(k,j2) = ijk2(k,j2+2)
                                enddo
                                do j1=1,n1
                                    cov(j1,j2) = cov(j1,j2+1)
                                enddo
                            enddo
                        endif
                        goto 801 ! for safety
                    endif
                endif
            enddo
        enddo
        n = 0
        do i1=1,n1
            do i2=1,n2
                if ( cov(i1,i2).gt.1e30 ) then
                    n = n + 1
                    cov(i1,i2) = 0
                endif
            enddo
        enddo
        print '(a,i6,a)','svd: set ',n,' undefined covariances to zero'
        print *,'still to do: compute SVD...'
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
        call getlonwindow(lon1,lon2,x1,x2,xx,nx,lon1c,lon2c,.false.)
        call getlatwindow(lat1,lat2,y1,y2,yy,ny,lat1c,lat2c,.false.)
        if ( lon1c.gt.1e33 .or. lat1c.ge.1e33 ) then
            x1 = 1
            lon1c = xx(1)
            y1 = 1
            lat1c = yy(1)
        endif
        print *,'cutting out longitude ',x1,x2,lon1c,lon2c
        print *,'cutting out latitude  ',y1,y2,lat1c,lat2c
        call makeabsent (data,npermax,yrbeg,yrend)
        do yr=max(yrbeg,firstyr),min(lastyr,yrend)
            do mo=1,nperyear
                data(mo,yr) = field(x1,y1,1,mo,yr)
            enddo
        enddo
        call printdatfile(6,data,npermax,nperyear,yrbeg,yrend)
        end

      subroutine processoptions(field,nx,ny,nz,nperyear,firstyr,lastyr
     +     ,mens,i12)
      implicit none
#include "params.h"
#include "getopts.inc"
      integer nx,ny,nz,nperyear,firstyr,lastyr,mens,i12
      real field(nx,ny,nz,nperyear,firstyr:lastyr,0:mens)
      integer jx,jy,jz,iens,i,j
      real fxy(npermax,yrbeg:yrend)
      call makeabsent(fxy,npermax,yrbeg,yrend)
      do jz=1,nz
            do jy=1,ny
                call keepalive(jy+(jz-1)*ny,ny*nz)
                do jx=1,nx
                    do iens=nens1,nens2
*
*                       create 1-D series from field
*
                        call makeabsent(fxy,npermax,yrbeg,yrend)
                        do i=yr1,yr2
                            do j=1,nperyear
                                if ( field(jx,jy,jz,j,i,iens).lt.1e30 )
     +                               then
                                    fxy(j,i) = field(jx,jy,jz,j,i,iens)
                                endif
                            enddo
                        enddo
*                       
*                       take monthly anomalies
*                       
                        if ( i12.eq.1 .and. mdiff.gt.0 ) then
                            call mdiffit(fxy,npermax,nperyear,yrbeg
     +                           ,yrend,mdiff)
                        endif
                        if ( i12.eq.2 .and. mdiff2.gt.0 ) then
                            call mdiffit(fxy,npermax,nperyear,yrbeg
     +                           ,yrend,mdiff2)
                        endif
*                       
*                       sum
*
                        if ( i12.eq.1 .and. lsum.gt.1 ) then
                            call sumit(fxy,npermax,nperyear,yrbeg,yrend
     +                           ,lsum,'v')
                        endif
                        if ( i12.eq.2 .and. lsum2.gt.1 ) then
                            call sumit(fxy,npermax,nperyear,yrbeg,yrend
     +                           ,lsum2,'v')
                        endif
*
*                       log, sqrt
*
                        if ( i12.eq.1 .and. logscale .or. 
     +                       i12.eq.2 .and. logfield ) then
                            call takelog(fxy,npermax,nperyear,yrbeg
     +                           ,yrend)
                        endif
                        if ( i12.eq.1 .and. sqrtscale .or. 
     +                       i12.eq.2 .and. logfield ) then
                            call takesqrt(fxy,npermax,nperyear,yrbeg
     +                           ,yrend)
                        endif
*
*                       detrend
*
                        if ( ldetrend ) then
                            if ( lwrite ) print *,'Detrending field'
                            if ( lag1.eq.0 .and. lag2.eq.0 .or. m1.eq.0
     +                           .or.lsel.eq.nperyear ) then
                                call detrend(fxy,npermax,nperyear,
     +                               yrbeg,yrend,yr1,yr2,m1,m2,lsel)
                            else
                                call detrend(fxy,npermax,nperyear,
     +                               yrbeg,yrend,yr1,yr2,1,12,lsel)
                            endif
                        endif
*
*                       differentiate
*
                        if ( ndiff.gt.0 ) then
                            if ( lwrite ) print *,'Taking differences'
                            call diffit(fxy,npermax,nperyear,yrbeg
     +                           ,yrend,ndiff)
                        endif
*       
*                       normalize to s.d.
*
                        if ( lnormsd ) then
                            call normsd(fxy,npermax,nperyear,yrbeg,yrend
     +                           ,yr1,yr2)
                        else
                            call anomal(fxy,npermax,nperyear,yrbeg
     +                           ,yrend,yr1,yr2)
                        endif
*
*                       copy back to field
*
                        do i=yr1,yr2
                            do j=1,nperyear
                                field(jx,jy,jz,j,i,iens) = fxy(i,j)
                            enddo
                        enddo
                    enddo       ! iens
                enddo           ! jz
            enddo               ! jy
        enddo                   ! jx
        end
