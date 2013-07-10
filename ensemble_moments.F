        program ensemble_moments
*
*       compute the requested moments of an ensemble and write them out
*       as GrADS files
*
        implicit none
#include "params.h"
#include "netcdf.inc"
#include "getopts.inc"
        integer nyrmax,nlevmax,mensmax
        parameter (nyrmax=121,nlevmax=1,mensmax=1)
        integer n,nx,ny,nz,nt,firstyr,lastyr,firstmo,nvars,
     +          ivars(2,nvmax),jvars(6,nvmax),ncid,endian,
     +          status,nperyear
        real xx(nxmax),yy(nymax),zz(nzmax),undef,
     +        field(nxmax,nymax,12,nyrmax,mensmax),
     +        res(nxmax,nymax,12,nyrmax)
        character infile*255,datfile*255,outfile*255
     +        ,vars(nvmax)*10,lvars(nvmax)*40,title*255,
     +        units(nvmax)*10
        integer iargc,llen
*
*       for f2c
        common /ccorrelatefield/ field,res
*
*       process command line
*
        n = iargc()
        if ( n.lt.3 ) then
            write(0,*) 'usage: ensemble_moments infile.[ctl|nc] [1234] '
     +            //'outfile.[ctl|nc]'
            stop
        endif
        call getarg(1,infile)
        if ( index(infile,'%').gt.0 .or. index(infile,'++').gt.0 ) then
            ensemble = .true.
            call filloutens(infile,0)
            inquire(file=infile,exist=lexist)
            if ( .not.lexist ) then
                mens1 = 1
                call filloutens(infile,1)
            else
                mens1 = 0
            endif
        else
            write(0,*) 'only for ensembles, not ',infile
            call abort
        endif
        if ( lwrite ) print *,'getmomentsfield: nf_opening file '
     +        ,infile(1:llen(infile))
        status = nf_open(infile,nf_nowrite,ncid)
        if ( status.ne.nf_noerr ) then
            call parsectl(infile,datfile,nxmax,nx,xx,nymax,ny,yy,nzmax
     +            ,nz,zz,nt,nperyear,firstyr,firstmo,undef,endian,title
     +            ,1,nvars,vars,ivars,lvars,units)
            ncid = -1
            if ( ensemble ) then
                do mens=1,nensmax
                    call getarg(1,line)
                    call filloutens(line,mens)
                    inquire(file=line,exist=lexist)
                    if ( .not.lexist ) goto 100
                enddo
  100           continue
                mens = mens - 1
                write(*,*) 'located ',mens-mens1+1
     +               ,' ensemble members<br>'
            endif
        else
            call parsenc(infile,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax
     +            ,nz,zz,nt,nperyear,firstyr,firstmo,undef,title,1,nvars
     +            ,vars,jvars,lvars,units)
            if ( ensemble ) then
                do mens=1,nensmax
                    call getarg(1,line)
                    call filloutens(line,mens)
                    status = nf_open(line,nf_nowrite,ncid)
                    if ( status.ne.nf_noerr ) goto 200
                enddo
  200           continue
                mens = mens - 1
                write(*,*) 'located ',mens-mens1+1
     +               ,' ensemble members<br>'
            endif
        endif
        lastyr = firstyr + (firstmo+nt-2)/nperyear
*       process arguments
        call getopts(3,n-1,nperyear,yrbeg,yrend,.true.,mens1,mens)
        call getarg(1,infile)
*
*       check dimensions
*
        if ( nx*ny*nz*nperyear*(lastyr-firstyr+1)*(nens2-nens1+1)
     +        .gt.nxmax*nymax*12*nyrmax*mensmax ) then
            write(0,*) 'correlatefield: error: field too large '
            write(0,*) '  nx       = ',nx
            write(0,*) '  ny       = ',ny
            write(0,*) '  nz       = ',nz
            write(0,*) '  nperyear = ',nperyear
            write(0,*) '  years    = ',firstyr,lastyr
            write(0,*) '  ensemble = ',nens1,nens2
            write(0,*) 'total request',
     +            nx*ny*nz*nperyear*(lastyr-firstyr+1)*(nens2-nens1+1)
            write(0,*) 'available    ',
     +            nxmax*nymax*12*nyrmax*mensmax
            write(*,*) 'correlatefield: error: field too large '
            write(*,*) '  nx       = ',nx
            write(*,*) '  ny       = ',ny
            write(*,*) '  nz       = ',nz
            write(*,*) '  nperyear = ',nperyear
            write(*,*) '  years    = ',firstyr,lastyr
            write(*,*) '  ensemble = ',nens1,nens2
            write(*,*) 'total request',
     +            nx*ny*nz*nperyear*(lastyr-firstyr+1)*(nens2-nens1+1)
            write(*,*) 'available    ',
     +            nxmax*nymax*12*nyrmax*mensmax
            call abort
        endif
        if ( nx*ny*nz .gt. nxmax*nymax*nlevmax ) then
            write(0,*) 'correlatefield: error: fields too large: '
     +            ,nx,ny,nz,nxmax,nymax,nlevmax
            write(*,*) 'correlatefield: error: fields too large: '
     +            ,nx,ny,nz,nxmax,nymax,nlevmax
            call abort
        endif
*
*       to save on RAM usage
*
        call gfield(datfile,ncid,field,res,nx,xx,ny,yy,nz,zz,nt
     +        ,nperyear,firstyr,lastyr,firstmo,undef,endian,jvars)
        end

        subroutine gfield(datfile,ncid,field,res,nx,xx,ny,yy,nz,zz,nt
     +        ,nperyear,firstyr,lastyr,firstmo,undef,endian,jvars)
*
*       break to use field() compactly
*
        implicit none
#include "recfac.h"
#include "params.h"
#include "getopts.inc"
        real absent
        parameter (absent=3e33)
*
        integer ncid,endian,nx,ny,nz,nt,nperyear,firstyr,lastyr,
     +          firstmo,jvars(6,nvmax)
        real field(nx,ny,nz,nperyear,firstyr:lastyr,nens1:nens2),
     +          res(nx,ny,nz,nperyear,firstyr:lastyr),
     +          undef,xx(nx),yy(ny),zz(nz)
        character datfile*(*)
*
        integer jx,jy,jz,i,j,jj,j1,j2,k,m,n,month,yr,imoment,ldir,
     +          nvars,ivars(2,nvmax),iens
        real fxy(npermax,yrbeg:yrend,0:nensmax),
     +        xmom(0:4),var,ddata(npermax*(yrend-yrbeg+1)),perc
        logical lexist
        character outfile*255,line*255,yesno*1,dir*255,string*10
     +        ,vars(nvmax)*10,lvars(nvmax)*40,title*255
        integer iargc,llen,rindex
*
        call getarg(2,line)
        if (  ichar(line(1:1)).ge.ichar('0') .and. 
     +        ichar(line(1:1)).le.ichar('9') ) then
            read(line,*,err=901) imoment
            if ( imoment.lt.-2 .or. imoment.gt.4 ) goto 901
        elseif ( line(1:4).eq.'adev' ) then
            imoment = 0
        elseif ( line(1:3).eq.'ave' .or. line(1:4).eq.'mean' ) then
            imoment = 1
        elseif ( line(1:3).eq.'sdm' .or. line(1:6).eq.'s.d./m' ) then
            imoment = -2
        elseif ( line(1:2).eq.'sd' .or. line(1:4).eq.'s.d.' ) then
            imoment = 2
        elseif ( line(1:4).eq.'skew' ) then
            imoment = 3
        elseif ( line(1:4).eq.'kurt' .or. line(1:4).eq.'curt' ) then
            imoment = 4
        elseif ( line(1:4).eq.'perc' ) then
            imoment = -1
            call getarg(3,line)
            read(line,*,err=902) perc
        else
            goto 901
        endif
        if ( lag1.ne.0 .or. lag2.ne.0 ) print *
     +        ,'getmomentsfield: lags do not make sense'
        if ( dump ) print *,'getmomentsfield: dump not supported'
        if ( plot ) print *,'getmomentsfield: plot not supported'
        if ( lks ) print *,'getmomentsfield: K-S not supported'
        if ( lconting ) print *,'getmomentsfield: contingency '//
     +        'tables not supported'
        do i=1,indxuse
            if ( lincl(i) ) print *,'getmomentsfield: what do ',
     +          'you mean with ',strindx(i),'?'
        enddo
*       range of years
        yr1 = max(yr1,firstyr)
        yr2 = min(yr2,firstyr + (firstmo+nt-2)/nperyear)
*
*       read field, change absent values to our convention
*
        if ( ensemble ) then
*           put the %% back in datfile...
            if ( nens2.lt.10 ) then
                i = 1
            elseif ( nens2.lt.100 ) then
                i = 2
            elseif ( nens2.lt.1000 ) then
                i = 3
            else
                write(0,*)
     +                'correlatefield: cannot handle ensembles up to '
     +                ,nens2,' yet'
                call abort
            endif
            string = '0000000'
            j = rindex(datfile,string(1:i))
            if ( j.eq.0 ) then
                write(0,*) 'correlatefield: error: cannot find '
     +                ,string(1:i),' in ',datfile(1:llen(datfile))
                call abort
            endif
            do k=j,j+i-1
                datfile(k:k) = '%'
            enddo
        endif
        do iens=nens1,nens2
            call keepalive(iens-nens1+1,nens2-nens1+1)
            if ( ncid.eq.-1 ) then
                dir=datfile
                if ( ensemble ) call filloutens(dir,iens)
                print *,'looking for '//dir(1:llen(dir))
                inquire(file=dir,exist=lexist)
                if ( .not.lexist ) then
                    print *,'looking for '//dir(1:llen(dir))//'.gz'
                    inquire(file=dir(1:llen(dir))//'.gz',exist=lexist)
                    if ( .not.lexist ) then
                        nens2 = iens-1
                        if ( nens2.ge.nens1 ) then
                            write(0,*) 'Found ensemble 0 to ',nens2
     +                            ,'<br>'
                            goto 5
                        else
                            write(0,*) 'Cannot locate file '
     +                            ,dir(1:llen(dir))
                            call abort
                        endif
                    endif
                endif
                if ( lwrite ) then
                    print *,'opening file ',dir(1:llen(dir))
                endif
                call zreaddatfile(dir,field(1,1,1,1,firstyr,iens),
     +                nx,ny,nz,nx,ny,nz,nperyear,firstyr,lastyr,
     +                firstyr,firstmo,nt,undef,endian,lwrite,yr1,yr2,1,1
     +                )
            else
                if ( nz.ne.1 ) then
                    write(0,*) 'cannot read 3D netCDF files yet'
                    call abort
                endif
                if ( ensemble ) then
                    write(0,*)
     +                    'cannot handle ensembles of netcdf files yet'
                    call abort
                endif
                call readncfile(ncid,field,nx,ny,nx,ny,nperyear,firstyr
     +                ,lastyr,firstyr,firstmo,nt,undef,lwrite,yr1,yr2
     +                ,jvars)
            endif
        enddo
    5   continue
*       open output file
*
        call getarg(iargc(),outfile)
        inquire(file=outfile,exist=lexist)
        if ( lexist ) then
            print *,'output file ',outfile(1:index(outfile,' ')-1),
     +            ' already exists, overwrite? [y/n]'
            read(*,'(a)') yesno
            if (  yesno.ne.'y' .and. yesno.ne.'Y' .and. 
     +            yesno.ne.'j' .and. yesno.ne.'J' ) then
                stop
            endif            
            open(2,file=outfile)
            close(2,status='delete')
        endif
        if ( index(outfile,'.nc').ne.0 ) then
            print *,'netCDF output not yet ready'
            stop
        else
            i = index(outfile,'.ctl')
            if ( i.ne.0 ) then
                datfile = outfile(:i-1)//'.dat'
            else
                datfile = outfile
            endif
            open(unit=2,file=datfile,form='unformatted',access='direct'
     +            ,recl=recfac*nx*ny*nz,err=920)
        endif
*
*       loop over gridpoints
*
        do jz=1,nz
            do jy=1,ny
                call keepalive(jy,ny)
                do jx=1,nx
                    do month=0,nperyear
                        res(jx,jy,jz,month) = absent
                    enddo
*
*                   create 1-D series from field
*
                    n = 0
                    do iens=nens1,nens2
                        do i=yr1,yr2
                            do j=1,nperyear
                                fxy(j,i,iens) = 
     +                                  field(jx,jy,jz,j,i,iens)
                                if ( fxy(j,i,iens).lt.0.9*absent )
     +                                  n = n+1
                            enddo
                        enddo
                    enddo
                    if ( n.lt.3 ) then
                        if ( lwrite ) print '(a,2i5)',
     +                          'no valid points at ',jx,jy
                        goto 800
                    endif
                    do iens=nens1,nens2
                        do i=yrbeg,yr1-1
                            do j=1,nperyear
                                fxy(j,i,iens) = 3e33
                            enddo
                        enddo
                        do i=yr2+1,yrend
                            do j=1,nperyear
                                fxy(j,i,iens) = 3e33
                            enddo
                        enddo
*
*                       sum
*
                        if ( lsum.gt.1 ) then
                            call sumit(fxy(1,yrbeg,iens),npermax,
     +                          nperyear,yrbeg,yrend,lsum,oper)
                        endif
*
*                       log,sqrt
*
                        if ( logscale ) then
                            do i=yr1,yr2
                                do j=1,nperyear
                                    if ( fxy(j,i,iens).lt.1e33 .and.
     +                                   fxy(j,i,iens).gt.0 )
     +                                  then
                                        fxy(j,i,iens) =
     +                                          log10(fxy(j,i,iens))
                                    else
                                        fxy(j,i,iens) = 3e33
                                    endif
                                enddo
                            enddo
                        endif
                        if ( sqrtscale ) then
                            do i=yr1,yr2
                                do j=1,nperyear
                                    if ( fxy(j,i,iens).lt.1e33 .and. 
     +                                   fxy(j,i,iens).ge.0 )
     +                                  then
                                        fxy(j,i,iens) =
     =                                          sqrt(fxy(j,i,iens))
                                    else
                                        fxy(j,i,iens) = 3e33
                                    endif
                                enddo
                            enddo
                        endif
*
*                       detrend
*
                        if ( ldetrend ) then
                            if ( lwrite ) print *,'Detrending field'
                            call detrend(fxy(1,yrbeg,iens),npermax,
     +                           nperyear,yrbeg,yrend,yr1,yr2,m1,m2,lsel
     +                           )
                        endif
*
*                       differentiate
*
                        if ( ndiff.ne.0 ) then
                            if ( lwrite ) print *,'Taking differences'
                            call diffit(fxy(1,yrbeg,iens),npermax,
     +                          nperyear,yrbeg,yrend,ndiff)
                        endif
*
*                       anomalies
*
                        if ( anom .or. lsel.gt.1 .and. ndiff.le.0 )
     +                      then
                            call anomal(fxy(1,yrbeg,iens),npermax,
     +                          nperyear,yrbeg,yrend,yr1,yr2)
                        endif
*
*                       normalize to s.d.
*
                        if ( lnormsd ) then
                            call normsd(fxy(1,yrbeg,iens),npermax,
     +                          nperyear,yrbeg,yrend,yr1,yr2)
                        endif
                    enddo       ! iens
*
*                   get moments
*
                    do month=m1,m2
                        if ( month.eq.0 ) then
                            j1 = 1
                            j2 = nperyear
                        else
                            j1 = month
                            j2 = month + lsel - 1
                        endif
*
*                       fill linear arrays without absent values
*                       and compute moment
*
                        n = 0
                        do iens=nens1,nens2
                            do yr=yr1-1,yr2
                                do jj=j1,j2
                                    j = jj
                                    call normon(j,yr,i,nperyear)
                                    if ( i.lt.yr1 .or.i.gt.yr2 )
     +                                  goto 710
                                    if (  fxy(j,i,iens).lt.absent/3.and. 
     +                                    fxy(j,i,iens).lt.maxindx .and. 
     +                                    fxy(j,i,iens).gt.minindx) then
                                        n = n+1
                                        ddata(n) = fxy(j,i,iens)
                                    endif
  710                               continue
                                enddo
                            enddo
                        enddo
                        if ( month.eq.0 .and. 
     +                        n.lt.minfac*min(nt,nperyear*(yr2-yr1+1))
     +                        .or. 
     +                        month.ne.0 .and.
     +                        n.lt.minfac*min(nt/nperyear,yr2-yr1+1)
     +                        .or.
     +                        n.lt.minnum ) then
                            if ( lwrite ) print '(a,2i5,2i3,a,2i6)'
     +                          ,'not enough valid points at ',jx,jy
     +                          ,month,': ',n,nt
                            goto 790
                        endif
*
                        m = month-m1
                        if ( imoment.ne.-1 ) then
*
*                           call Numerical Recipes routine
*
                            call moment(ddata,n,xmom(1),xmom(0),xmom(2),
     +                          var,xmom(3),xmom(4))
                            if ( imoment.ge.0 ) then
                                res(jx,jy,jz,m) = xmom(imoment)
                            elseif ( imoment.eq.-2 ) then
                                if ( xmom(1).ne.0 ) then
                                    res(jx,jy,jz,m) = xmom(2)/xmom(1)
                                else
                                    res(jx,jy,jz,m) = 3e33
                                endif
                            else
                                write(0,*)
     +                              'getmomentfield: internal error:',
     +                              ' imoment = ',imoment
                                call abort
                            endif
                        else
*
*                           sort results and find perc-th percentile
*                           using the routine already written for the
*                           contingency table
*
                            call nrsort(n,ddata)
                            call getcut(res(jx,jy,jz,m),perc,n,ddata)
                        endif
  790                   continue    ! valid point/month
                    enddo           ! month
  800               continue        ! valid point
                enddo               ! nx
            enddo                   ! ny
        enddo                       ! nz
*
*       write output field in GrADS format
*
        print *,'writing output'
        do month=m1,m2
            m = month-m1
            if ( lwrite ) print *,'writing records ',6*m+1,'-',6*m+6
     +            ,' of fields ',m,' of size ',nx*ny*recfa4
            write(2,rec=m+1) (((res(jx,jy,jz,m),jx=1,nx),jy=1,ny),
     +          jz=1,nz)
        enddo
        close(2)
        if ( index(outfile,'.ctl').ne.0 ) then
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
                if ( line(1:ldir).eq.dir(1:ldir) ) then
                    title(n:) = line(ldir+1:)
                else
                    title(n:) = line
                endif
                n = llen(title) + 2
            enddo
            nvars = 1
            if ( imoment.eq.-1 ) then
                vars(1) = 'perc'
                write(lvars(1),'(f5.2,a)') perc,'% percentile'
            elseif ( imoment.eq.0 ) then
                vars(1) = 'adev'
                lvars(1) = 'absolute deviation'
            elseif ( imoment.eq.1 ) then
                vars(1) = 'ave'
                lvars(1) = 'mean'
            elseif ( imoment.eq.2 ) then
                vars(1) = 'sd'
                lvars(1) = 'standard deviation'
            elseif ( imoment.eq.-2 ) then
                vars(1) = 'sdm'
                lvars(1) = 'standard deviation / mean'
            elseif ( imoment.eq.3 ) then
                vars(1) = 'skew'
                lvars(1) = 'skewness'
            elseif ( imoment.eq.4 ) then
                vars(1) = 'kurt'
                lvars(1) = 'kurtosis'
            else
                write(0,*) 'getmomentsfield: error: imoment = ',imoment
                call abort
            endif
            do i=1,nvars
                ivars(1,i) = 0
                ivars(2,i) = 99
            enddo
*           give correlations dates in 0-1
            if ( m1.eq.0 ) then
                i = 0
                j = 12
            else
                i = 1
                j = m1
            endif
            call writectl(outfile,datfile,nx,xx,ny,yy,nz,zz
     +            ,1+(m2-m1),nperyear,i,j,absent,title,nvars,vars,ivars
     +            ,lvars,units)
        endif
*
*       error messages
*
        goto 999
  901   print *,'getmomentsfield: error reading moment [1-4] from '
     +        ,line(1:llen(line))
        call abort
  902   print *,'getmomentsfield: error reading percentile from '
     +        ,line(1:llen(line))
        call abort
  903   print *,'error reading date from file ',line(1:index(line,' ')-1
     +        ),' at record ',k
        call abort
  904   print *,'error cannot locate field file file ',line(1:index(line
     +        ,' ')-1)
        call abort
  920   print *,'error cannot open new correlations file '
     +        ,datfile(1:index(datfile,' ')-1)
        call abort
  999   continue
        end
