        program quantiles_field
!
!       Compute quantiles of a set of fields using all ensemble members
!       but weighing each model equally.
!
        implicit none
#include "params.h"
#include "netcdf.inc"
        integer nvarmax,ntmax
        parameter (nvarmax=20,ntmax=12)
        integer i,j,nens,iens,fyr,yr1,yr2,imod,nmod,nskip,nensmin
        integer ix,iy,iz,mo,lastmo,ntvarid,itimeaxis(ntmax),it,status
        integer mens1,mens,ncid,nx,ny,nz,nt,nperyear,ivar
     +       ,firstyr,firstmo,endian,nvars,jvars(6,nvarmax)
     +       ,ivars(2,nvarmax)
        integer,allocatable :: nnens(:),nnmod(:)
        real xx(nxmax),yy(nymax),zz(nzmax),undef
        real,allocatable :: field(:,:,:,:,:,:),quantiles(:,:,:,:,:)
     +       ,sdmod(:,:,:,:,:),sd(:,:,:,:),sdsd(:,:,:,:)
        character infile*255,outfile*255,datfile*255,lz(3)*20,ltime*100
     +       ,title*255,history*20000,vars(nvarmax)*40
     +       ,lvars(nvarmax)*140,svars(nvarmax)*100,units(nvarmax)*100
     +       ,cell_methods(nvarmax)*100
        character var*20,model*20,exp*20,rip*20,yrs*20,oldmodel*20
     +       ,oldexp*5
        character models(500)*20,inunts*40,p*2,oldp*2
        logical lwrite,lstandardunits,lcomputesd
        integer iargc

        lcomputesd = .false.
        call getarg(iargc(),outfile)
        if ( index(outfile,'CMIP5ext').ne.0 ) lcomputesd = .true.
        lwrite = .false. ! .true.
        lstandardunits = .false.
        nensmin = 3             ! minimum number of ensemble members to base a sd on

        if ( iargc().lt.3 ) then
            write(0,*) 'usage: quantiles_field infile1 infile2 ...
     +           outfile'
            write(0,*) 'computes a few standrd quantiles in a '//
     +           'multi-model ensemble. Each file should correspond '//
     +           'to one ensemble member.'
            call abort
        end if
!
!       read data
!
        nens = iargc() - 1
        nskip  = 0
        do iens=1,nens
            call getarg(iens,infile)
            call getmetadata(infile,mens1,mens,ncid,datfile,nxmax,nx
     +           ,xx,nymax,ny,yy,nzmax,nz,zz,lz,nt,nperyear,firstyr
     +           ,firstmo,ltime,undef,endian,title,history,nvarmax,nvars
     +           ,vars,jvars,lvars,svars,units,cell_methods,lwrite)
            if ( nt+firstmo-1.gt.nperyear ) then
                write(0,*) 'quantiles_field: error: can only handle one'
     +               ,' year of data in each field, not ',firstmo,+nt-1
                call abort
            end if
            if ( iens.eq.1 ) then
                fyr = firstyr
                yr1 = firstyr
                yr2 = firstyr
                allocate(field(nx,ny,nz,nperyear,fyr:fyr,nens))
                allocate(nnmod(nens))
                allocate(nnens(nens))
                allocate(quantiles(nx,ny,nz,nperyear,-6:6))
                nmod = 0
                nnmod = 0
                oldmodel = ' '
                oldp = ' '
                oldexp = ' '
            end if
            call parsename(infile,var,model,exp,rip,yrs,lwrite)
            if ( model.eq.'modmean' .or. model.eq.'onemean' ) then
                write(0,*) 'skipping ',trim(model)
                nskip = nskip + 1
                cycle
            end if
            i = index(rip,'p')
            p = rip(i+1:)
            if ( model.ne.oldmodel .or. p.ne.oldp .or. exp.ne.oldexp )
     +           then
                nmod = nmod + 1
                nnmod(nmod) = 0
                models(nmod) = model
                oldmodel = model
                oldp = p
                oldexp = exp
            end if
            nnmod(nmod) = nnmod(nmod) + 1
            call readfield(ncid,infile,datfile,field(1,1,1,1,fyr,iens)
     +           ,nx,ny,nz,nperyear,fyr,fyr,0,0,nx,ny,nz,nperyear,yr1
     +           ,yr2,firstyr,firstmo,nt,undef,endian,vars,units
     +           ,lstandardunits,lwrite)
        end do
        if ( lcomputesd ) then
            allocate(sdmod(nx,ny,nz,nperyear,nmod))
            allocate(sd(nx,ny,nz,nperyear))
            allocate(sdsd(nx,ny,nz,nperyear))
        end if
        nens = nens - nskip
        lastmo = firstmo + nt - 1
!
!       how many ensemble members per model?
!
        print *,'Found ',nmod,' model/experiment combinations with '//
     +       'in all ',nens,' ensemble members<br>'
        iens = 0
        do imod=1,nmod
            if ( lwrite) print *,imod,nnmod(imod),' ',trim(models(imod))
            do i=1,nnmod(imod)
                iens = iens + 1
                nnens(iens) = nnmod(imod)
            end do
        end do
        if ( iens.ne.nens ) then
            write(0,*) 'quantiles_field: error: something went wrong '
     +           ,iens,nens
            call abort
        end if
!
!       compute min,max,quantiles
!
        call getweightedquantile(field,nnens,nx,ny,nz,nperyear,nens
     +       ,firstmo,lastmo,nmod,quantiles,lwrite)

        if ( lcomputesd ) then
!
!           compute natural variability
!
            sdmod = 3e33
            iens = 0
            do imod=1,nmod
                iens = iens + nnmod(imod)
                if ( nnmod(imod).ge.nensmin ) then
                    if ( lwrite ) print *,'using model ',models(imod),
     +                   ' to compute natural variability'
                    call getnaturalvariability(field(1,1,1,1,fyr,iens)
     +                   ,nx,ny,nz,nperyear,nnmod(imod),firstmo,lastmo
     +                   ,sdmod(1,1,1,1,imod),lwrite)
                end if
            end do
            call averagefields(sdmod,nx,ny,nz,nperyear,nmod,
     +           firstmo,lastmo,sd,sdsd,lwrite)
        end if
!
!       write output to file
!
        nvars = 0
        call definevariable(nvars,vars,lvars,'min'
     +       ,'minimum of all ensemble members')
        call definevariable(nvars,vars,lvars,'p025','2.5% quantile')
        call definevariable(nvars,vars,lvars,'p05','5% quantile')
        call definevariable(nvars,vars,lvars,'p10','10% quantile')
        call definevariable(nvars,vars,lvars,'p17','17% quantile')
        call definevariable(nvars,vars,lvars,'p25','25% quantile')
        call definevariable(nvars,vars,lvars,'p50','50% quantile')
        call definevariable(nvars,vars,lvars,'p75','75% quantile')
        call definevariable(nvars,vars,lvars,'p83','83% quantile')
        call definevariable(nvars,vars,lvars,'p90','90% quantile')
        call definevariable(nvars,vars,lvars,'p95','95% quantile')
        call definevariable(nvars,vars,lvars,'p975','97.5% quantile')
        call definevariable(nvars,vars,lvars,'max'
     +       ,'maximum of all ensemble members')
        if ( lcomputesd ) then
            call definevariable(nvars,vars,lvars,'sd'
     +           ,'mean standard deviation from intra-model spread')
            call definevariable(nvars,vars,lvars,'sdsd'
     +           ,'intra-model standard deviation of standard '//
     +           'deviation from intra-model spread')
        end if
        do ivar=1,nvars
            units(ivar) = units(1)
            ivars(1,ivar) = nz
            ivars(2,ivar) = 99
        end do

        title = 'Quantiles of the multi-model ensemble'
        call getarg(iargc(),outfile)
        call writenc(outfile,ncid,ntvarid,itimeaxis,ntmax,nx,xx,ny,yy
     +       ,nz,zz,nt,nperyear,firstyr,firstmo,3e33,title,nvars
     +       ,vars,ivars,lvars,units,0,0)
        do it=1,nt
            do ivar=1,nvars
                if ( ivar.le.13 ) then
                    call writencslice(ncid,ntvarid,itimeaxis,ntmax
     +                   ,ivars(1,ivar),quantiles(1,1,1,firstmo+it-1
     +                   ,ivar-7),nx,ny,nz,nx,ny,nz,it,1)
                else if ( lcomputesd .and. ivar.eq.nvars-1 ) then
                    call writencslice(ncid,ntvarid,itimeaxis,ntmax
     +                   ,ivars(1,ivar),sd(1,1,1,firstmo+it-1),nx,ny,nz
     +                   ,nx,ny,nz,it,1)
                else if ( lcomputesd .and. ivar.eq.nvars ) then
                    call writencslice(ncid,ntvarid,itimeaxis,ntmax
     +                   ,ivars(1,ivar),sdsd(1,1,1,firstmo+it-1),
     +                   nx,ny,nz,nx,ny,nz,it,1)
                else
                    write(0,*) 'quantiles_field: error: unknown ivar ='
     +                   ,ivar
                end if
            end do
        end do
        status = nf_close(ncid)
!
        end program

        subroutine parsename(file,var,model,exp,rip,yrs,lwrite)
!
!       parse a standard file name
!
        implicit none
        integer i0,i
        character file*(*),var*(*),model*(*),exp*(*),rip*(*),yrs*(*)
        character measure*10,type*10
        logical lwrite
        i0 = 1 + index(file,'/',.true.)
        if ( lwrite ) print *,'file name = ',trim(file(i0:))
        i = i0 + index(file(i0:),'_') - 1
        i0 = i + 1
        call parsename_next(file,i0,i,measure,'measure',lwrite)
        call parsename_next(file,i0,i,var,'var',lwrite)
        call parsename_next(file,i0,i,type,'type',lwrite)
        if ( type.eq.'A1' ) then ! CMIP3
            i = index(file,'_144') - 1
            if ( ichar(file(i:i)).ge.ichar('0') .and.
     +           ichar(file(i:i)).le.ichar('9') .and.
     +           ichar(file(i-1:i-1)).ge.ichar('0') .and.
     +           ichar(file(i-1:i-1)).le.ichar('9') .and.
     +           file(i-2:i-2).eq.'_' ) then
                rip = file(i-1:i-1)
                i = i-3
            else
                rip = ' '
            end if
            model = file(i0:i)
            exp = 'sresa1b'
            yrs = ' '
        else ! CMIP5
            call parsename_next(file,i0,i,model,'model',lwrite)
            call parsename_next(file,i0,i,exp,'exp',lwrite)
            call parsename_next(file,i0,i,rip,'rip',lwrite)
            call parsename_next(file,i0,i,yrs,'yrs',lwrite)
        end if
!
        end subroutine

        subroutine parsename_next(file,i0,i,var,name,lwrite)
        implicit none
        integer i0,i
        character file*(*),var*(*),name*(*)
        logical lwrite

        i = i0 + index(file(i0:),'_') - 1
        if ( i.eq.i0 ) then
            write(0,*) 'parsename: error: cannot find ',name,' in ',
     +           trim(file),i0
            call abort
        end if
        var = file(i0:i-1)
        if ( lwrite ) print *,trim(name),' = ',trim(var)
        i0 = i + 1
!
        end subroutine

        subroutine getweightedquantile(field,nnens,nx,ny,nz,nperyear
     +       ,nens,firstmo,lastmo,nmod,quantiles,lwrite)
!
!       compute the weighted quantile of field using the number of
!       ensemble members in nn as weighting information
!
        implicit none
        integer nnens(nens),nx,ny,nz,nperyear,nens,firstmo,lastmo,nmod
        real field(nx,ny,nz,nperyear,1,nens)
        real quantiles(nx,ny,nz,nperyear,-6:6)
        logical lwrite
        integer ix,iy,iz,mo,iens,iquant
        real quant(-6:6)
        real,allocatable :: point(:)
!
        allocate(point(nens))
        if ( lwrite ) print *,'getweightedquantile: ',nx,ny,nz,nperyear
     +       ,nens
        do mo=firstmo,lastmo
            do iz=1,nz
                do iy=1,ny
                    do ix=1,nx
                        do iens=1,nens
                            point(iens) = field(ix,iy,iz,mo,1,iens)
                        end do
                        call getweightedquant(point,nnens
     +                       ,nens,nmod,quant,lwrite)
                        do iquant = -6,6
                            quantiles(ix,iy,iz,mo,iquant) =
     +                           quant(iquant)
                        end do
                    end do
                end do
            end do
        end do
        deallocate(point)
!               
        end subroutine

        subroutine getnaturalvariability(field,nx,ny,nz,nperyear,nens
     +       ,firstmo,lastmo,sdmod,lwrite)
!
!       compute an estimate of the natural variability from the spread
!       of the nens ensemble members of one model
!
        implicit none
        integer nx,ny,nz,nperyear,nens,firstmo,lastmo
        real field(nx,ny,nz,nperyear,1,nens),sdmod(nx,ny,nz,nperyear)
        logical lwrite
        integer ix,iy,iz,mo,iens,n
        real ave,adev,sdev,var,skew,curt
        real,allocatable :: xx(:)
        
        if ( nens.lt.3 ) then
            write(0,*) 'getnaturalvariability: error: need more than '
     +           ,nens,' ensemble members'
            sdmod = 3e33
            return
        end if
        allocate(xx(nens))
        do mo=firstmo,lastmo
            do iz=1,nz
                do iy=1,ny
                    do ix=1,nx
                        n = 0
                        do iens=1,nens
                            if ( field(ix,iy,iz,mo,1,iens).lt.1e33 )
     +                           then
                                n = n + 1
                                xx(n) = field(ix,iy,iz,mo,1,iens)
                            end if
                        end do
                        ! Numerical recipes routine
                        if ( n.ge.2 ) then
                            call moment(xx,n,ave,adev,sdev,var,skew
     +                           ,curt)
                            sdmod(ix,iy,iz,mo) = sdev
                        else
                            sdmod(ix,iy,iz,mo) = 3e33
                        end if
                    end do
                end do
            end do
        end do
        deallocate(xx)
!
        end subroutine

        subroutine averagefields(sdmod,nx,ny,nz,nperyear,nmod
     +       ,firstmo,lastmo,sd,sdsd,lwrite)
!
!       average the fields in sdmod (one per model) with equal weights
!       also computes the sd of the sd.
!
        implicit none
        integer nx,ny,nz,nperyear,nmod,firstmo,lastmo
        real sdmod(nx,ny,nz,nperyear,nmod),sd(nx,ny,nz,nperyear),
     +       sdsd(nx,ny,nz,nperyear)
        logical lwrite
        integer ix,iy,iz,mo,imod,n
        real s,s2
        
        do mo=firstmo,lastmo
            do iz=1,nz
                do iy=1,ny
                    do ix=1,nx
                        n = 0
                        s = 0
                        do imod=1,nmod
                            if ( sdmod(ix,iy,iz,mo,imod).lt.1e33 ) then
                                n = n + 1
                                s = s + sdmod(ix,iy,iz,mo,imod)
                            end if
                        end do
                        if ( n.lt.1 ) then
                            sd(ix,iy,iz,mo) = 3e33
                            sdsd(ix,iy,iz,mo) = 3e33
                        else
                            s = s/n
                            sd(ix,iy,iz,mo) = s
                            s2 = 0
                            do imod=1,nmod
                                if ( sdmod(ix,iy,iz,mo,imod).lt.1e33 )
     +                               then
                                    s2 = s2 +
     +                                   (sdmod(ix,iy,iz,mo,imod)-s)**2
                                end if
                            end do
                            if ( n.lt.2 ) then
                                sdsd(ix,iy,iz,mo) = 3e33
                            else
                                sdsd(ix,iy,iz,mo) = sqrt(s2/(n-1))
                            end if
                        end if
                    end do
                end do
            end do
        end do

        end subroutine

        subroutine definevariable(nvars,vars,lvars,name,lname)
        implicit none
        integer nvars
        character vars(nvars+1)*(*),lvars(nvars+1)*(*),
     +       name*(*),lname*(*)
        nvars = nvars + 1
        vars(nvars) = name
        lvars(nvars) = lname
        !!!print *,'vars(',nvars,') = ',trim(name)
        end subroutine

        
