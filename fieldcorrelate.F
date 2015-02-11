*  #[ fieldcorrelate:
        program fieldcorrelate
*
*       program to compute a time series of filed correlations
*       between two fields
*
        implicit none
#include "params.h"
#include "recfac.h"
#include "netcdf.inc"
        integer nyrmax
        parameter(nyrmax=50)
        integer j,n,ncid1,ncid2,nx1,ny1,nz1,nt1,nper1,firstyr1,firstmo1,
     +       nx2,ny2,nz2,nt2,nper2,firstyr2,firstmo2,nvars,ivars1(2,1),
     +       ivars2(2,1),jvars1(6,nvmax),jvars2(6,nvmax),endian1,
     +       endian2,status,nxf,nyf,nperyear,firstyr,lastyr
        real  xx1(nxmax),yy1(nymax),zz1(nzmax),
     +        xx2(nxmax),yy2(nymax),zz2(nzmax),
     +        u1,u2
        real,allocatable :: field1(:,:,:,:),field2(:,:,:,:)
        logical lwrite
*
        character line*256,datfile1*256,datfile2*256,
     +       title1*256,title2*256,vars1(1)*20,lvars1(1)*80,
     +       vars2(1)*20,lvars2(1)*80,units1(1)*20,units2(1)*20
        integer iargc,llen
*
*       check arguments
*
        allocate(field1(nxmax,nymax,12,nyrmax))
        allocate(field2(nxmax,nymax,12,nyrmax))
        lwrite = .FALSE.
        n = iargc()
        if ( n.lt.4 ) then
            print *,'usage: fieldcorrelate '//
     +            'field1.[ctl|nc] field2.[ctl|nc] '//
     +            'corr|sign|regr|intercept '//
     +            '[lag n '//
     +            '[log|sqrt|rank] [sum|ave|max|min|sel n] '//
     +            '[minfac r] [begin yr] [end yr] '
***     +            '[lt cut] [gt cut] [diff [npoints]] '
            stop
        endif
        call getarg(1,line)
        if ( lwrite ) print *,'fieldcorrelate: nf_opening file '
     +        ,line(1:llen(line))
        status = nf_open(line,nf_nowrite,ncid1)
        if ( status.ne.nf_noerr ) then
            call parsectl(line,datfile1,nxmax,nx1,xx1,nymax,ny1,yy1
     +            ,nzmax,nz1,zz1,nt1,nper1,firstyr1,firstmo1,u1,endian1
     +            ,title1,1,nvars,vars1,ivars1,lvars1,units1)
            ncid1 = -1
        else
            datfile1 = line
            call parsenc(line,ncid1,nxmax,nx1,xx1,nymax,ny1,yy1
     +            ,nzmax,nz1,zz1,nt1,nper1,firstyr1,firstmo1,u1,title1,1
     +            ,nvars,vars1,jvars1,lvars1,units1)
        endif
*
        call getarg(2,line)
        if ( lwrite ) print *,'fieldcorrelate: nf_opening file '
     +        ,line(1:llen(line))
        status = nf_open(line,nf_nowrite,ncid2)
        if ( status.ne.nf_noerr ) then
            call parsectl(line,datfile2,nxmax,nx2,xx2,nymax,ny2,yy2
     +            ,nzmax,nz2,zz2,nt2,nper2,firstyr2,firstmo2,u2,endian2
     +            ,title2,1,nvars,vars2,ivars2,lvars2,units2)
            ncid2 = -1
        else
            datfile2 = line
            call parsenc(line,ncid2,nxmax,nx2,xx2,nymax,ny2,yy2
     +            ,nzmax,nz2,zz2,nt2,nper2,firstyr2,firstmo2,u2,title2,1
     +            ,nvars,vars2,jvars2,lvars2,units2)
        endif
*
        nxf = max(nx1,nx2)
        nyf = max(ny1,ny2)
        firstyr = max(firstyr1,firstyr2)
        nperyear = max(nper1,nper2)
        lastyr = min(firstyr1 + (nt1-1)/nperyear,firstyr2 + (nt2-1)
     +        /nperyear)
        if ( nxf*nyf*nperyear*(lastyr-firstyr+1).gt.nxmax*nymax*12
     +        *nyrmax) then
            write(0,*) 'fieldcorrelate: error: fields too large: '
     +            ,nxf,nyf,firstyr,yrend,nxmax,nymax,nyrmax
            write(*,*) 'fieldcorrelate: error: fields too large: '
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
     +        ,firstmo1,u1,endian1,jvars1,vars1,units1,
     +        datfile2,ncid2,field2,nx2,xx2,ny2,yy2,nt2,nper2,firstyr2
     +        ,firstmo2,u2,endian2,jvars2,vars2,units2)
*
        end
*  #] fieldcorrelate:
*  #[ cfieldfield:
        subroutine cfieldfield(nxf,nyf,nperyear,firstyr,
     +        datfile1,ncid1,field1,nx1,xx1,ny1,yy1,nt1,nper1,firstyr1
     +        ,firstmo1,u1,endian1,jvars1,vars1,units1,
     +        datfile2,ncid2,field2,nx2,xx2,ny2,yy2,nt2,nper2,firstyr2
     +        ,firstmo2,u2,endian2,jvars2,vars2,units2)
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
     +        ,nperyear,firstyr:yrend),xx1(nxmax),yy1(nymax),xx2(nxmax)
     +        ,yy2(nymax),u1,u2
        character*(*) datfile1,datfile2,vars1(1)*(*),units1(1)*(*),
     +       vars2(1)*(*),units2(1)*(*)
*
        integer ntmax
        parameter(ntmax=1000)
        integer i,j,k,n,jx,jy,mo,m,yr,y,lastyr,nt,nrec,nvars,ivars(2,6)
     +        ,ldir,nx,ny,ncid,iout,x1,x2,y1,y2
        logical xrev1,xwrap1,yrev1,xrev2,xwrap2,yrev2,equal
        real r(npermax,yrbeg:yrend),prob(npermax,yrbeg:yrend),
     +       a(npermax,yrbeg:yrend),b(npermax,yrbeg:yrend),
     +       da(npermax,yrbeg:yrend),db(npermax,yrbeg:yrend),
     +       f1(ndata),f2(ndata),w1(ndata),w2(ndata),
     +       df1(ndata),df2(ndata),sig,
     +       adata,sxx,aindx,syy,sxy,df,d,zd,z,probd,absent,chi2,q,sum,
     +       fac,wx1(ndata),wy1(ndata),wx2(ndata),wy2(ndata),xx(nxmax),
     +       yy(nymax),wx(nxmax),wy(nxmax),rmse,bias,var1,var2
        character output*4
        common /c_findx/ df1,f2
*
        parameter (absent=3e33)
        logical lexist
        character line*80,yesno*1,string*10,file*255,outfile*255,
     +        datfile*255,title*255,vars(6)*10,lvars(6)*40,dir*255
        integer iargc,llen
*
*       check if arrays big enough
        if ( nxf*nyf.gt.nxmax*nymax ) then
            write(0,*)'fieldcorrelate: recompile with nxmax*nymax '
     +            ,nxf*nyf
            write(*,*)'fieldcorrelate: recompile with nxmax*nymax '
     +            ,nxf*nyf
            call abort
        endif
*
*       save time on the initialization - but not too much.
        lastyr = min(firstyr1 + (nt1-1)/nperyear,firstyr2 + (nt2-1)
     +        /nperyear)
        nt = nperyear*(lastyr-firstyr+1)
        n = iargc()
        call getarg(3,output)
        if ( output.ne.'corr' .and. output.ne.'sign' .and.
     +       output.ne.'regr' .and. output.ne.'inte' .and.
     +       output.ne.'rmse' .and. output.ne.'bias' .and.
     +       output.ne.'rvar' .and. output.ne.'var1' .and.
     +       output.ne.'nrms' .and. output.ne.'nbia' .and.
     +       output.ne.'nvar' .and. output.ne.'var2' ) goto 902
        call getopts(4,n-1,nperyear,yrbeg,yrend,.true.,0,0)
        if ( lag1.lt.0 ) print *,'(field1 leading field2)'
        if ( lag1.ne.lag2) then
            write(0,*) 'fieldcorrelate: cannot handle lag range'
            lag2 = lag1
        endif
        if ( dump ) write(0,*) 'fieldcorrelate: dump not supported'
        if ( plot ) write(0,*) 'fieldcorrelate: plot not supported'
        if ( lks ) write(0,*) 'fieldcorrelate: K-S not supported'
        if ( mdiff.ne.0 .or. mdiff2.ne.0 ) write(0,*) 
     +      'fieldcorrelate: monthly anomalies not supported'
        if ( lconting ) write(0,*) 'fieldcorrelate: contingency '//
     +        'tables not supported'
        do i=1,indxuse
            if ( lincl(i) ) write(0,*) 'fieldcorrelate: what do ',
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
        iF ( output.ne.'rmse' .and. output.ne.'bias' .and.
     +       output.ne.'var1' .and. output.ne.'var2' .and. 
     +       output.ne.'nrms' .and. output.ne.'nbia'.and. 
     +       output.ne.'nvar' .and. output.ne.'rvar' ) then
            lrmse = .false.
        else
            lrmse = .true.
        endif
*
*       init
*
*       compute minfac if it has not been set explicitly
*
        if ( minfac.lt.0 ) then
*           heuristic, gives 0.25 for 150 yrs, 0.5 for 50 yrs, 0.75 for 20yrs
            minfac = 0.40
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
        if ( lstandardunits ) then
            call makestandardfield(field1,nxf,nyf,1,nperyear,firstyr
     +           ,yrend,nx1,ny1,1,nperyear,yr1,yr2,vars1(1),units1(1)
     +           ,lwrite)
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
        if ( lstandardunits ) then
            call makestandardfield(field2,nxf,nyf,1,nperyear,firstyr
     +           ,yrend,nx2,ny2,1,nperyear,yr1,yr2,vars2(1),units2(1)
     +           ,lwrite)
        endif
*
*       cut out region of interest
*
        call getxyprop(xx1,nx1,yy1,ny1,xrev1,yrev1,xwrap1)
        call getlatlonwindow(lat1,lat2,lon1,lon2,xx1,nx1,xwrap1,1,yy1
     +        ,ny1,1,x1,x2,y1,y2,lwrite)
        call enscutoutwindow(x1,x2,y1,y2,xx1,nx1,xwrap1,xrev1,1,yy1,ny1
     +       ,1,wx1,wy1,field1,nxf,nyf,0,0,nperyear,firstyr,yrend,yr1
     +       ,yr2,lwrite)
        call getxyprop(xx2,nx2,yy2,ny2,xrev2,yrev2,xwrap2)
        call getlatlonwindow(lat1,lat2,lon1,lon2,xx2,nx2,xwrap2,1,yy2
     +        ,ny2,1,x1,x2,y1,y2,lwrite)
        call enscutoutwindow(x1,x2,y1,y2,xx2,nx2,xwrap2,xrev2,1,yy2,ny2
     +       ,1,wx2,wy2,field2,nxf,nyf,0,0,nperyear,firstyr,yrend,yr1
     +       ,yr2,lwrite)
*
*       interpolate fields to common grid
*
        call interpu(
     +        field1(1,1,1,yr1),xx1,yy1,nx1,ny1,
     +        field2(1,1,1,yr1),xx2,yy2,nx2,ny2,
     +        xx,nx,yy,ny,yr1,yr2,yr1,yr2,nxf,nyf,nperyear,intertype,
     +        lwrite)
        call copyweights(xx,xx1,xx2,wx,wx1,wx2,nx,nx1,nx2)
        call copyweights(yy,yy1,yy2,wy,wy1,wy2,ny,ny1,ny2)
*
*       differentiate
*
        if ( ndiff.gt.0 ) then
            if ( lwrite ) print *,'Taking differences'
            write(0,*) 'not yet implemented'
            call abort
        endif
*
*       take anomalies
*
        if ( anom ) then
            if ( .TRUE..or.lwrite ) write(0,*) 'Taking anomalies<p>'
            call fieldanomal(field1(1,1,1,yr1),nxf,nyf,nperyear,yr1,yr2
     +            ,nx,ny)
            call fieldanomal(field2(1,1,1,yr1),nxf,nyf,nperyear,yr1,yr2
     +            ,nx,ny)
        endif
*
*       loop over time
*
        write(0,*) 'Correlating<p>'
        call makeabsent(r,npermax,yrbeg,yrend)
        call makeabsent(prob,npermax,yrbeg,yrend)
        call makeabsent(a,npermax,yrbeg,yrend)
        call makeabsent(b,npermax,yrbeg,yrend)
        call makeabsent(da,npermax,yrbeg,yrend)
        call makeabsent(db,npermax,yrbeg,yrend)
        call printdatfile(6,r,npermax,nperyear,yrbeg,yrend)
        do yr=yr1,yr2
            call keepalive(yr-yr1+1,yr2-yr1+1)
            do mo=1,nperyear
*
*               create 1-D series from fields
*
                n = 0
                do jy=1,ny
                    do jx=1,nx
                        n = n + 1
                        f1(n) = 0
                        w1(n) = 0
                        do i=1,lsum
                            m=mo+i-1
                            call normon(m,yr,y,nperyear)
                            if ( y.gt.yr2 ) then
                                f1(n) = 3e33
                                w1(n) = 3e33
                            else
                                f1(n) = f1(n) + field1(jx,jy,m,y)
                                w1(n) = w1(n) + wx(jx)*wy(jy)
                            endif
                        enddo
                        f2(n) = 0
                        w2(n) = 0
                        do i=1,lsum2
                            m=mo+i-1-lag1
                            call normon(m,yr,y,nperyear)
                            if ( y.gt.yr2 .or.y.lt.yr1 ) then
                                f2(n) = 3e33
                                w2(n) = 3e33
                            else
                                f2(n) = f2(n) + field2(jx,jy,m,y)
                                w2(n) = w2(n) + wx(jx)*wy(jy)
                            endif
                        enddo
                        if ( f1(n).gt.1e30 .or. f2(n).gt.1e30 )
     +                        then
                            n = n - 1
                        endif
                    enddo
                enddo
                if ( n.lt.3 .or. n.lt.minfac*nx*ny ) then
                    if ( lwrite ) print '(a,i3,i5,i8,f10.0)',
     +                    'not enough valid points in field1 ',mo,yr,n,
     +                    minfac*nx*ny
                    goto 790
                endif
                if ( lrmse ) then
                    if ( lwrite ) print *,'computing bias,rmse,vars'
*
*                   compute RMS error, bias, ratio variances
*
                    call getrms(f1,w1,f2,w2,n,bias,adata,var1,aindx,var2
     +                   ,rmse)
                    r(mo,yr) = rmse
                    a(mo,yr) = bias
                    b(mo,yr) = var1
                    db(mo,yr)= var2
                    if ( var2.gt.0 ) then
                        prob(mo,yr) = rmse/sqrt(var2)
                        da(mo,yr) = bias/sqrt(var2)
                    endif
*
*                   correlate!
*
                elseif ( lrank ) then
                    if ( lwrite ) print *,'computing rank correlation'
                    sum = 1 + decor/nperyear
                    call spearx(f1,f2,n,df2,df1,d,zd,probd,r(mo,yr),
     +                    prob(mo,yr),sum,adata,sxx,aindx,syy)
                else
                    if ( lwrite ) print *,'computing correlation,fit'
                    df = n/(1 + decor/nperyear) - 2
                    if ( df.le.0 ) then
                        if ( lwrite ) print '(a,f5.2,2i5)'
     +                        ,'error: df <= 0: ',df,jx,jy
                        goto 790
                    endif
                    if ( anom ) then
                        call apearsnxx(f2,f1,n,r(mo,yr),prob(mo,yr),z,
     +                          aindx,syy,adata,sxx,sxy,df)
                    else
                        call pearsncross(f2,f1,n,r(mo,yr),prob(mo,yr),z,
     +                       aindx,syy,adata,sxx,sxy,df,ncrossvalidate)
                    endif
                    if ( sxx.eq.0 .or. syy.eq.0 ) then
                        r(mo,yr) = absent
                        prob(mo,yr) = absent
                    else
                        call fit(f2,f1,n,sig,0,
     +                        a(mo,yr),b(mo,yr),
     +                        da(mo,yr),db(mo,yr),chi2,q)
                    endif
                endif
                if ( output.eq.'sign' .and. 
     +               abs(prob(mo,yr)).lt.0.9*absent ) then
                    prob(mo,yr) = 1-prob(mo,yr)
                endif
                if ( lwrite ) then
                    if ( lrmse ) then
                        print '(a,i3,i5,i3,a,i6,a,6g12.4)'
     +                       ,'point ',mo,yr,lag1,' OK (',n,'): '
     +                       ,r(mo,yr),prob(mo,yr),a(mo,yr),da(mo,yr)
     +                       ,b(mo,yr),db(mo,yr)
                    else
                        print '(a,i3,i5,i3,a,i6,a,6f9.4)'
     +                       ,'point ',mo,yr,lag1,' OK (',n,'): '
     +                       ,r(mo,yr),prob(mo,yr),a(mo,yr),da(mo,yr)
     +                       ,b(mo,yr),db(mo,yr)
                    endif
                endif
  790           continue        ! valid point
            enddo               ! mo
        enddo                   ! yr
        call getarg(iargc(),outfile)
        open(1,file=outfile,err=920)
        write(1,'(4a)') '# fieldcorrelations of '
     +       ,datfile1(1:llen(datfile1)),' and '
     +       ,datfile2(1:llen(datfile2))
        write(1,'(a,f6.2,a,f6.2)') '# longitude ',lon1,':',lon2
        write(1,'(a,f6.2,a,f6.2)') '# latitude  ',lat1,':',lat2
        write(1,'(a)') '# KNMI Climate Explorer, http://climexp.knmi.nl'
        if ( output.eq.'corr' ) then
            write(1,'(a)') '# correlations'
            call printdatfile(1,r,npermax,nperyear,yrbeg,yrend)
        elseif ( output.eq.'sign' ) then
            write(1,'(a)') '# significances correlations'
            call printdatfile(1,prob,npermax,nperyear,yrbeg,yrend)
        elseif ( output.eq.'slop' .or. output.eq.'regr' ) then
            write(1,'(a)') '# regression'
            call printdatfile(1,b,npermax,nperyear,yrbeg,yrend)
        elseif ( output.eq.'inte' ) then
            write(1,'(a)') '# intercept'
            call printdatfile(1,a,npermax,nperyear,yrbeg,yrend)
        elseif ( output.eq.'rmse' ) then
            write(1,'(a)') '# rms error (w/o bias)'
            call printdatfile(1,r,npermax,nperyear,yrbeg,yrend)
        elseif ( output.eq.'nrms' ) then
            write(1,'(a)') '# rms error/standard deviation field 2'
            call printdatfile(1,prob,npermax,nperyear,yrbeg,yrend)
        elseif ( output.eq.'bias' ) then
            write(1,'(a)') '# bias (field1-field2)'
            call printdatfile(1,a,npermax,nperyear,yrbeg,yrend)
        elseif ( output.eq.'nbia' ) then
            write(1,'(a)') '# bias/standard deviation field 2'
            call printdatfile(1,da,npermax,nperyear,yrbeg,yrend)
        elseif ( output.eq.'var1' ) then
            write(1,'(a)') '# variances field1'
            call printdatfile(1,b,npermax,nperyear,yrbeg,yrend)
        elseif ( output.eq.'var2' ) then
            write(1,'(a)') '# variances field2'
            call printdatfile(1,db,npermax,nperyear,yrbeg,yrend)
        elseif ( output.eq.'rvar' ) then
            write(1,'(a)') '# ratio variances var1/var2'
            do yr=yr1,yr2
                do m=1,nperyear
                    r(mo,yr) = b(mo,yr)/db(mo,yr)
                enddo
            enddo
            call printdatfile(1,r,npermax,nperyear,yrbeg,yrend)
        elseif ( output.eq.'nvar' ) then
            write(1,'(a)') '# var1/var2-1 or 1-var2/var1'
            do yr=yr1,yr2
                do m=1,nperyear
                    if ( b(mo,yr).lt.db(mo,yr) ) then
                        r(mo,yr) = 1 - db(mo,yr)/b(mo,yr)
                    else
                        r(mo,yr) = b(mo,yr)/db(mo,yr) - 1
                    endif
                enddo
            enddo
            call printdatfile(1,r,npermax,nperyear,yrbeg,yrend)
        else
            goto 902
        endif
        close(1)
*
*       error messages
*
        goto 999
  902   write(*,*) 'fieldcorrelate: error: expecting '//
     +            'corr|prob|regr|intercept'
        call abort
  903   print *,'error reading date from file ',line(1:index(line,' ')-1
     +        ),' at record ',k
        call abort
  920   print *,'error cannot open new correlations file '
     +        ,outfile(1:index(outfile,' ')-1)
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
*  #] cfieldfield:

        subroutine copyweights(xx,xx1,xx2,wx,wx1,wx2,nx,nx1,nx2)
        integer nx,nx1,nx2
        real xx(nx),xx1(nx),xx2(nx),wx(nx),wx1(nx),wx2(nx)
        integer i
        logical equal
        
        equal = .true.
        do i=1,nx
            if ( xx(i).ne.xx1(i) .and. 
     +           xx(i).ne.xx1(i)+360 .and.
     +           xx(i).ne.xx1(i)-360 ) then
                equal = .false.
                exit
            end if
        end do
        if ( equal ) then
            do i=1,nx
                wx(i) = wx1(i)
            end do
        else
            equal = .true.
            do i=1,nx
                if ( xx(i).ne.xx1(nx+1-i) .and.
     +               xx(i).ne.xx1(nx+1-i)+360 .and.
     +               xx(i).ne.xx1(nx+1-i)-360 ) then
                    equal = .false.
                    exit
                end if
            end do
            if ( equal ) then
                do i=1,nx
                    wx(i) = wx1(nx+1-i)
                end do
            else
                equal = .true.
                do i=1,nx
                    if ( xx(i).ne.xx2(i) .and.
     +                   xx(i).ne.xx2(i)+360 .and.
     +                   xx(i).ne.xx2(i)-360 ) then
                        equal = .false.
                        exit
                    end if
                end do
                if ( equal ) then
                    do i=1,nx
                        wx(i) = wx2(i)
                    end do
                else
                    equal = .true.
                    do i=1,nx
                        if ( xx(i).ne.xx2(nx+1-i) .and.
     +                       xx(i).ne.xx2(nx+1-i)+360 .and.
     +                       xx(i).ne.xx2(nx+1-i)-360 ) then
                            equal = .false.
                            exit
                        end if
                    end do
                    if ( equal ) then
                        do i=1,nx
                            wx(i) = wx2(nx+1-i)
                        end do
                    else
                        write(0,*) 'copyweights: error: cannot find ',
     +                       'grid back'
                        write(0,*) 'xx  = ',(xx(i),i=1,nx)
                        write(0,*) 'xx1 = ',(xx1(i),i=1,nx1)
                        write(0,*) 'xx2 = ',(xx2(i),i=1,nx2)
                        call abort
                    end if
                end if
            end if
        end if
        end
