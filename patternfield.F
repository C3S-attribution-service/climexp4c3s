        program patternfield
*
*       computes the projection of a pattern onto a field, which gives a
*       time series on stdout
*
        implicit none
#include "params_s.h"
#include "netcdf.inc"
*	maximum number of variables in the pattern file
        integer nvarmax,nyrmax
        parameter (nvarmax=80,nyrmax=85)
        integer i,j,n,ncid1,ncid2,nx1,ny1,nz1,nt1,nper1,firstyr1
     +        ,firstmo1,nx2,ny2,nz2,nt2,nper2,firstyr2,firstmo2,nvars
     +        ,ivars(2,nvarmax),jvars1(6,nvmax),jvars2(6,nvarmax),ivar
     +        ,month,endian1,endian2,status,nperyear,nxf,nyf
        real  xx1(nxmax),yy1(nymax),zz1(nzmax),
     +        xx2(nxmax),yy2(nymax),zz2(nzmax),
     +        field1(nxmax,nymax,12,nyrmax),
     +        field2(nxmax,nymax,npermax,2),
     +        u1,u2
        logical lwrite
*
*	for f2c
	common /cpatternfield/ field1,field2
*
        character line*80,datfile1*256,datfile2*256,variable*20,
     +        title1*256,title2*256
        character vars1(nvarmax)*20,lvars1(nvarmax)*80,
     +       units1(nvarmax)*20
        character vars2(nvarmax)*20,lvars2(nvarmax)*80,
     +       units2(nvarmax)*20
        integer iargc,llen
        lwrite = .FALSE.
*
*       check arguments
*
        n = iargc()
        if ( n.lt.4 ) then
            print *,'usage: patternfield field.[ctl|nc] '//
     +            'pattern.[ctl|nc] variable month [minfac n]'
            stop
        endif
        call killfile(title1,datfile1,datfile2,0)
        call getarg(iargc(),line)
        call tolower(line)
        if ( line(1:5).eq.'debug' .or. line(1:6).eq.'lwrite' ) then
            lwrite = .true.
            print *,'turned debug output on'
        end if
        call getarg(1,line)
        if ( lwrite ) print *,'patternfield: nf_opening file '
     +        ,line(1:llen(line))
        status = nf_open(line,nf_nowrite,ncid1)
        if ( status.ne.nf_noerr ) then
            call parsectl(line,datfile1,nxmax,nx1,xx1,nymax,ny1,yy1
     +            ,nzmax,nz1,zz1,nt1,nper1,firstyr1,firstmo1,u1,endian1
     +            ,title1,1,nvars,vars1,ivars,lvars1,units1)
            ncid1 = -1
        else
            datfile1 = line
            call parsenc(line,ncid1,nxmax,nx1,xx1,nymax,ny1,yy1
     +            ,nzmax,nz1,zz1,nt1,nper1,firstyr1,firstmo1,u1,title1,1
     +            ,nvars,vars1,jvars1,lvars1,units1)
        endif
*
        call getarg(2,line)
        if ( lwrite ) print *,'patternfield: nf_opening file '
     +        ,line(1:llen(line))
        status = nf_open(line,nf_nowrite,ncid2)
        if ( status.ne.nf_noerr ) then
            call parsectl(line,datfile2,nxmax,nx2,xx2,nymax,ny2,yy2
     +            ,nzmax,nz2,zz2,nt2,nper2,firstyr2,firstmo2,u2
     +           ,endian2,title2,nvarmax,nvars,vars2,ivars,lvars2,
     +            units2)
            ncid2 = -1
        else
            datfile2 = line
            call parsenc(line,ncid2,nxmax,nx2,xx2,nymax,ny2,yy2
     +            ,nzmax,nz2,zz2,nt2,nper2,firstyr2,firstmo2,u2,title2
     +            ,nvarmax,nvars,vars2,jvars2,lvars2,units2)
        endif
        if ( nper1.ne.nper2 ) then
            write(0,*) 'patternfield: error: time scales do not agree '
     +            ,nper1,nper2
            call abort
        endif
        nperyear = max(nper1,nper2)
*
        call getarg(3,variable)
        do ivar=1,nvars
            if ( vars2(ivar).eq.variable ) then
                goto 100
            endif
        enddo
        write(0,*) 'patternfield: cannot locate '
     +        ,variable(1:llen(variable)),' in pattern file '
     +        ,line(1:llen(line))
        write(0,*) 'I only have ',(vars2(ivar),ivar=1,nvars)
        call abort
  100   continue
        if ( ncid2.ge.0 ) then
*           make sure the variable is the first one in the jvar array
            if ( ivar.gt.1 ) then
                do i=1,5
                    jvars2(i,1) = jvars2(i,ivar)
                enddo
            endif
        endif
        if ( lwrite ) print *,'located ',variable(1:llen(variable)),ivar
        call getarg(4,line)
        read(line,*,err=903) month
        if ( month.lt.0 .or. month.gt.nperyear ) goto 903
        if ( lwrite ) print *,'picking pattern for month ',month
        do i=len(datfile1),1,-1
            if ( datfile1(i:i).eq.'/' ) goto 200
        enddo
  200   continue
        i = i + 1
        print '(3a,i2,4a)','# patternfield: projecting variable '
     +        ,variable(1:llen(variable)),', month ',month
     +        ,' of pattern ',datfile2(1:llen(datfile2)),' on field '
     +        ,datfile1(1:llen(datfile1))
        print '(2a)','# ',title1(1:llen(title1))
        print '(2a)','# ',title2(1:llen(title2))
        do j=1,2
            print '(a)','#'
        enddo
        nxf = max(nx1,nx2)
        nyf = max(ny1,ny2)
        if ( nxf*nyf*nt1.gt.nxmax*nymax*12*nyrmax ) then
            write(0,*) 'patternfield: error: field too large: ',nxf,nyf
     +            ,nt1,' vs ',nxmax,nymax,12,nyrmax
            write(*,*) 'patternfield: error: field too large: ',nxf,nyf
     +            ,nt1,' vs ',nxmax,nymax,12,nyrmax
            write(0,*) 'patternfield: error: field too large: ',nxf*nyf
     +            *nt1,' vs ',nxmax*nymax*12*nyrmax
            write(*,*) 'patternfield: error: field too large: ',nxf*nyf
     +            *nt1,' vs ',nxmax*nymax*12*nyrmax
            call abort
        endif
        if ( nxf*nyf*nper2.gt.nxmax*nymax*npermax ) then
            write(0,*) 'patternfield: error: pattern too large: ',nxf
     +            ,nyf,nper2,' vs ',nxmax,nymax,npermax
            write(*,*) 'patternfield: error: pattern too large: ',nxf
     +            ,nyf,nper2,' vs ',nxmax,nymax,npermax
            call abort
        endif
        call patfield(nperyear,nxf,nyf,nvarmax,ivar
     +        ,nvars,month,datfile1,ncid1,field1,nx1,xx1,ny1,yy1,nt1
     +        ,firstyr1,firstmo1,u1,endian1,jvars1,datfile2,ncid2,field2
     +        ,nx2,xx2,ny2,yy2,nt2,firstyr2,firstmo2,u2,endian2,jvars2
     +       ,vars1,units1,lwrite)
*
*       error messages
        goto 999
  902   write(0,*) 'error: firstyr1,lastyr1 = ',firstyr1,firstyr1 + (nt1
     +        -1)/nperyear,' should be between ',yrbeg,' and ',yrend
        write(0,*) '       recompile if this is too restrictive'
        call abort
  903   write(0,*)'error: month should be between 0 and ',nperyear
     +        ,', not ',line(1:llen(line))
        call abort
  999   continue
        end

        subroutine patfield(nperyear,nxf,nyf,nvarmax,ivar,nvars,month,
     +        datfile1,ncid1,field1,nx1,xx1,ny1,yy1,nt1,firstyr1
     +        ,firstmo1,u1,endian1,jvars1,
     +        datfile2,ncid2,field2,nx2,xx2,ny2,yy2,nt2,firstyr2
     +        ,firstmo2,u2,endian2,jvars2,vars,units,lwrite)
*
*       break to use the arrays field1,field2 compactly and conserve RAM
*
        implicit none
#include "params.h"
#include "recfac.h"
*
        integer nperyear,nxf,nyf,nvarmax,ivar,nvars,month,
     +      ncid1,nx1,ny1,nt1,firstyr1,firstmo1,endian1,jvars1(6,nvmax),
     +      ncid2,nx2,ny2,nt2,firstyr2,firstmo2,endian2,
     +       jvars2(6,nvarmax)
        real field1(nxf,nyf,nperyear,firstyr1:yrend),
     +        field2(nxf,nyf,nperyear,0:1)
     +        ,xx1(nx1),yy1(ny1),xx2(nx2),yy2(ny2),u1,u2
        character*(*) datfile1,datfile2
        character vars(nvarmax)*(*),units(nvarmax)*(*)
        logical lwrite
*       
        integer yr1,yr2,nx,ny,nz,i,j,n,yr,mo,iskip
        real xx(nxmax),yy(nymax),zz(nymax),minfac
        real,allocatable :: var(:,:)
        logical anom,lstandardunits
        character*80 line
        integer iargc,llen
        anom = .false.
        lstandardunits = .false.
        
        minfac = 0.5
        do i=5,iargc()-1
            if ( iskip.gt.0 ) then
                iskip = iskip - 1
                cycle
            endif
            call getarg(i,line)
            if ( line(1:6).eq.'minfac' ) then
                call getarg(6,line)
                iskip = 1
                read(line,*,err=904) minfac
                if ( minfac.gt.1 ) minfac=minfac/100
            elseif ( line(1:5).eq.'debug' .or. line(1:6).eq.'lwrite' ) 
     +               then
                lwrite = .true.
                print *,'turned debug output on'
            elseif ( line(1:4).eq.'stan' ) then
                lstandardunits = .true.
                print '(a)','# converting to standard units'
            else
                write(0,*) 'error: do not understand argument ',line
                call abort
            endif
        enddo
*
*       read fields
*
        call keepalive(0,3)
        yr1 = firstyr1
        yr2 = min(yrend,firstyr1 + (firstmo1+nt1-1)/nperyear)
        if ( ncid1.eq.-1 ) then
            call readdatfile(datfile1,field1,nxf,nyf,nx1,ny1,nperyear
     +            ,firstyr1,yrend,firstyr1,firstmo1,nt1,u1,endian1
     +            ,lwrite,yr1,yr2,1,1)
        else
            if ( lwrite ) print *,'calling readncfile '
     +          ,datfile1(1:llen(datfile1)),ncid1,nxf,nyf,nx1,ny1
     +          ,nperyear,firstyr1,yrend,firstyr1,firstmo1,nt1,u1,lwrite
     +          ,yr1,yr2,jvars1
            call readncfile(ncid1,field1,nxf,nyf,nx1,ny1,nperyear
     +            ,firstyr1,yrend,firstyr1,firstmo1,nt1,u1,lwrite,yr1
     +            ,yr2,jvars1)
        endif
        if ( lstandardunits ) then
*           convert to standard units
            call makestandardfield(field1,nxf,nyf,1
     +           ,nperyear,firstyr1,yrend,nx1,ny1,1,nperyear
     +           ,firstyr1,yrend,vars(1),units(1),lwrite)
            if ( lwrite ) then
                print *,'patternfield: just after standard units'
                print *,'field1(',(nx1+1)/2,(ny1+1)/2,firstmo1
     +               ,firstyr1,') = ',field1((nx1+1)/2,(ny1+1)/2
     +               ,firstmo1,firstyr1)
            endif
        endif
        call keepalive(1,3)
        if ( ncid2.eq.-1 ) then
            call readdatfile(datfile2,field2,nxf,nyf,nx2,ny2,nperyear,0
     +            ,1,firstyr2,firstmo2,nt2,u2,endian2,lwrite,0,1,ivar
     +            ,nvars)
        else
            if ( lwrite ) print *,'calling readncfile '
     +          ,datfile2(1:llen(datfile2)),ncid2,nxf,nyf,nx2,ny2
     +          ,nperyear,yrbeg,yrend,firstyr2,firstmo2,nt2,u2,lwrite
     +          ,yr1,yr2,jvars2
            call readncfile(ncid2,field2,nxf,nyf,nx2,ny2,nperyear,0,1
     +            ,firstyr2,firstmo2,nt2,u2,lwrite,0,1,jvars2)
        endif
        call keepalive(2,3)
*
*       interpolate fields to common grid
*       
        call interpu(
     +        field1,xx1,yy1,nx1,ny1,
     +        field2,xx2,yy2,nx2,ny2,
     +        xx,nx,yy,ny,firstyr1,yr2,0,1,nxf,nyf,nperyear,0,lwrite)
        call keepalive(3,3)
        allocate(var(nperyear,yr1:yr2))
        call project(var,nperyear,yr1,yr2,0,0,xx,nx,yy,ny,
     +       field1,field2,nxf,nyf,nperyear,firstyr1,yrend,
     +       month,minfac,.false.,.false.,anom,lwrite)
        call printdatfile(6,var,nperyear,nperyear,yr1,yr2)
        return
*
*       error messages
*        
  901   write(0,*) 'error: interpolation not yet ready'
        write(0,*) xx1
        write(0,*) xx2
        call abort
  902   write(0,*) 'error: interpolation not yet ready'
        write(0,*) yy1
        write(0,*) yy2
        call abort
  904   write(0,*)'error: cannot read minfac from ',line(1:llen(line))
        call abort
        end
