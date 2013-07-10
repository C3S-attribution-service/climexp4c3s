        program convertmetadata
*
*       convert metadat from a netcdf object (most often over opendap)
*       into a form that my climexp scripts can use
*
        implicit none
#include "netcdf.inc"
#include "params.h"
        integer ntmax,nvarmax
        parameter (ntmax=60000,nvarmax=99)
        integer ncid,nx,ny,nz,nt,nperyear,firstyr
     +       ,firstmo,nens1,nens2,ntvars,ivars(6,nvarmax),iperyear
        real xx(nxmax),yy(nymax),zz(nzmax),undef
        character file*256,title*1024,vars(nvarmax)*40
     +       ,lvars(nvarmax)*80,units*40,oper*4,variable*20
     +       ,list*500,sourcelist(0:200,3)*60,creftimes(ntmax)*5
*       local variables
        integer status,ndims,nvars,ngatts,unlimdimid,varid,xtype
     +       ,ndimvar,dimids(nf_max_var_dims),natts,dimid,len,ix,iy,iz
     +       ,it,ie,i,j,k,n,l,iens,start(2),count(2),fmo,dy,mo,yr
     +       ,jul0,valid,invalid,nsource,isource,l1,l2,stringlength
        real*8 tt(ntmax),ttu(ntmax),tt1
        real x,xmin
        character name*(nf_max_name),dimname*(nf_max_name),months(12)*3,
     +       format*30,clwrite*1
        logical lwrite,foundtime,foundvalid
        integer llen,leap,iargc,julday
        data months
     +       /'jan','feb','mar','apr','may','jun'
     +       ,'jul','aug','sep','oct','nov','dec'/
*
        if ( iargc().lt.3 ) then
            print *,'usage: convertmetadata file list|convert'
            print *,'       variable [list]'
            print *,'       list: give list of values of variable, '//
     +           'possibly restricted to list'
            print *,'       convert: list all the indices for which '//
     +           'variable is equal to list'
            stop
        endif
        call getarg(1,file)
        do i=1,len(file)
            if ( file(i:i).eq.char(13) ) file(i:i) = ' '
        end do
        call getarg(2,oper)
        call getarg(3,variable)
        call getarg(4,list)
        foundtime = .false.
        lwrite = .false.
        call getenv('CONVERTMETADATA_LWRITE',clwrite)
        if ( index(clwrite,'T') + index(clwrite,'t') .gt.0 ) then
            lwrite = .true.
            print *,'convertmetadata: debug output requested'
        endif
        nx = 1
        xx(1) = 0
        ny = 1
        yy(1) = 0
        nz = 1
        zz(1) = 0
        ix = -1
        iy = -1
        iz = -1
        it = -1
        ie = -1
*
*       open file
*
        if ( lwrite ) print *,'convertmetadata: opening file '
     +       ,trim(file)
        status = nf_open(trim(file),nf_nowrite,ncid)
        if ( status.ne.nf_noerr ) call handle_err(status,file)
        if ( lwrite ) print *,'convertmetadata: opened with ncid = '
     +       ,ncid
        call gettitle(ncid,title,lwrite)
        call getnumbers(ncid,ndims,nvars,ngatts,unlimdimid,lwrite)
        call getdims(ncid,ndims,ix,nx,nxmax,iy,ny,nymax,iz,nz,nzmax,it
     +       ,nt,ntmax,ie,nens1,nens2,nensmax,lwrite)
*
*       loop over variables
*
        nsource = 0
        ntvars = 0
        do varid=1,nvars
*           get dimensions of variable
            status = nf_inq_var(ncid,varid,name,xtype,ndimvar,dimids
     +            ,natts)
            if ( status.ne.nf_noerr ) call handle_err(status
     +            ,'nf_inq_var')
            if ( lwrite ) then
                print *,'convertmetadata: variable: ',varid
                print *,'         name:     ',trim(name)
                print *,'         dims:     ',ndimvar,':',
     +                (dimids(i),i=1,ndimvar)
                print *,'         natts:    ',natts
            endif
            if ( index(name,'_bnd').ne.0 ) then
                if ( lwrite ) print *
     +               ,'convertmetadata: disregarding boundary '
     +               ,trim(name)
                cycle
            endif
*           what kind of variable do we have?
            if ( variable.eq.'variables' ) then
                if ( oper.eq.'list' ) then
!                   list all time-varying variables except uninteresting ones
                    if ( name.eq.'time' .or. name.eq.'reftime' .or. 
     +                   name.eq.'leadtime' ) then
                        if ( lwrite ) print *,'skipping ',trim(name)
                        cycle
                    endif
                    do i=1,ndimvar
                        if ( dimids(i).eq.it ) then
                            call gettextattopt(ncid,varid,'long_name'
     +                           ,dimname,lwrite)
                            if ( dimname.eq.' ' ) dimname = name
                            call gettextattopt(ncid,varid,'units',units
     +                           ,lwrite)
!                           there are GDS servers that do not have units...
                            if ( units.eq.' ' ) then
                                l1 = index(dimname,'[')
                                l2 = index(dimname,']')
                                if ( l1.ne.0 .and. l2.gt.l1 ) then
                                    if ( l2.gt.l1+1 ) then
                                        units = dimname(l1+1:l2-1)
                                    endif
                                    dimname(l1:l2) = ' '
                                endif
                            endif
!                           get rid of stars in dimname - these are
!                           expanded
                            l1 = index(dimname,'*',.true.)
                            if ( dimname(l1-1:l1).eq.'**' ) then
                                dimname = dimname(l1+2:)
                            else
                                dimname(l1:l1) = ' '
                            endif
                            print '(8a)',trim(name)
     +                           ,' (',trim(dimname)
     +                           ,') [',trim(units),']'
                        endif
                    enddo
                else
                    write(0,*) 'convertmetadata: cannot handle oper '
     +                   ,oper
                    write(*,*) 'convertmetadata: cannot handle oper '
     +                   ,oper
                    call abort
                endif
            elseif ( name.eq.variable .or. variable.eq.'ids' ) then
                if ( variable.eq.'level' .or. 
     +               variable.eq.'lev' .or. 
     +               variable.eq.'depth' .or. 
     +               variable.eq.'longitude' .or.
     +               variable.eq.'lon' .or.
     +               variable.eq.'latitude' .or.
     +               variable.eq.'lat' ) then
*                   first check whether the variable has any levels
*                   associated with it
                    do i=1,nvars
                        status = nf_inq_var(ncid,i,name,xtype,ndimvar
     +                       ,dimids,natts)
                        if ( status.ne.nf_noerr ) call handle_err(status
     +                       ,'nf_inq_var 2')
                        if ( lwrite ) then
                            print *,'comparing ',trim(name)
                            print *,'     with ',trim(list)
                        endif
                        if ( name.eq.list ) then
                            do j=1,ndimvar
                                if ( variable(1:3).eq.'lev' .and.
     +                               dimids(j).eq.iz .or.
     +                               variable.eq.'depth' .and.
     +                               dimids(j).eq.iz .or.
     +                               variable(1:3).eq.'lon' .and.
     +                               dimids(j).eq.ix .or.
     +                               variable(1:3).eq.'lat' .and.
     +                               dimids(j).eq.iy ) then
                                    goto 700
                                endif
                            enddo
!                           no levels, no output
                            goto 999
                        endif
                    enddo
                    write(0,*) 'convertmetadata: error: could not'//
     +                   ' find variable ',trim(list)
                    write(*,*) 'convertmetadata: error: could not'//
     +                   ' find variable ',trim(list)
                    call abort
 700                continue
                    status = nf_inq_var(ncid,varid,name,xtype,ndimvar
     +                   ,dimids,natts)
                    if ( lwrite ) print *,'found ',trim(variable)
                    if ( variable.eq.'level'     .and. dimids(1).ne.iz 
     +                   .or.
     +                   variable.eq.'depth'     .and. dimids(1).ne.iz 
     +                   .or.
     +                   variable.eq.'longitude' .and. dimids(1).ne.ix
     +                   .or.
     +                   variable.eq.'atitude'   .and. dimids(1).ne.iy
     +                   ) then
                        write(0,*) 'convertmetadata: error: expected '//
     +                       'first dimension of level to be '
     +                       ,trim(variable)
                        write(*,*) 'convertmetadata: error: expected '//
     +                       'first dimension of level to be '
     +                       ,trim(variable),iz,dimids
                        call abort
                    endif
!                   read data
                    status = nf_get_var_real(ncid,varid,xx)
                    if ( status.ne.nf_noerr ) call
     +                   handle_err(status,'nf_get_var_double levels')
                    if ( variable(1:3).eq.'lev' .or. 
     +                   variable.eq.'depth' ) then
                        n  = nz
                    elseif ( variable(1:3).eq.'lon' ) then
                        n = nx
                    elseif ( variable(1:3).eq.'lat' ) then
                        n = ny
                    endif
                    if ( lwrite ) then
                        do i=1,n
                            print *,i,xx(i)
                        enddo
                    endif
                    if ( oper.eq.'list' ) then
                        do i=1,n
                            print '(f15.2)',xx(i)
                        enddo
                        goto 999
                    else
                        write(0,*)
     +                       'convertmetadata: cannot handle oper ',oper
                        write(*,*)
     +                       'convertmetadata: cannot handle oper ',oper
                        call abort
                    endif
                elseif ( ( variable.eq.'ids' .or. variable.eq.name )
     +                   .and.
     +                   ( name.eq.'source' .or. 
     +                   name.eq.'experiment_id' .or.
     +                   name.eq.'institution' ) ) then
*                   translate ensemble numbers to model names and vv
                    nsource = nsource + 1
                    if ( name.eq.'experiment_id' ) then
                        isource = 1
                        stringlength=4 ! hard-coded
                    elseif ( name.eq.'source' ) then
                        isource = 2
                        stringlength=60 ! hard-coded
                    elseif ( name.eq.'institution' ) then
                        isource = 3
                        stringlength=15 ! hard-coded
                    endif
                    if ( lwrite ) then
                        print *,'found ',trim(variable)
                        print *,'ie,ndimvar,dimids = ',ie,ndimvar
     +                       ,dimids(1:ndimvar)
                        print *,'isource = ',isource
                    end if
                    if ( dimids(ndimvar).ne.ie ) then
                        write(0,*) 'convertmetadata: error: '//
     +                       'expected first dimension of source '//
     +                       'to be ensemble'
                        write(*,*) 'convertmetadata: error: '//
     +                       'expected first dimension of source '//
     +                       'to be ensemble',ie,dimids(1:2)
                        call abort
                    endif
                    if ( ndimvar.eq.2 ) then
!                       read data
                        count(1) = stringlength
                        count(2) = 1
                        do iens=nens1,nens2
                            start(1) = 1
                            start(2) = 1 + iens
                            if ( lwrite ) then
                                print *,'calling nf_get_vara_text with '
                                print *,'ncid,varid = ',ncid,varid
                                print *,'start = ',start(1:2)
                                print *,'count = ',count(1:2)
                            end if
                            status = nf_get_vara_text(ncid,varid
     +                           ,start,count,sourcelist(iens,isource))
                            if ( status.ne.nf_noerr ) call
     +                           handle_err(status,'nf_get_vara_text')
                            if ( lwrite ) print *,iens
     +                           ,trim(sourcelist(iens,isource))
                        enddo
                    else if ( ndimvar.eq.1 ) then ! new way to translate opendap into netcdf3
                        count(1) = 1 ! length of character array
                        do iens=nens1,nens2
                            start(1) = 1 + iens
                            status = nf_get_vara_text(ncid,varid
     +                           ,start,count,sourcelist(iens,isource))
                            if ( status.ne.nf_noerr ) call
     +                           handle_err(status,'nf_get_vara_text')
                            if ( lwrite ) print *,iens
     +                           ,trim(sourcelist(iens,isource))
                        enddo                        
                    else
                        write(0,*) 'expected one or two dimensions, '//
     +                       'not ',ndimvar
                        call abort
                    endif
!                   first get rid of the spaces
                    do iens=nens1,nens2
                        do i=1,len_trim(sourcelist(iens,isource))
                            if ( sourcelist(iens,isource)(i:i)
     +                           .eq.char(0) )
     +                           sourcelist(iens,isource)(i:i) = ' '
                        enddo
!                       special cases -- too long and boring
                        i = index(sourcelist(iens,isource)
     +                       ,', System')
                        if ( i.gt.0 ) sourcelist(iens,isource)(i:) = ' '
                        i = index(sourcelist(iens,isource)
     +                       ,', ENSEMBLES')
                        if ( i.gt.0 ) sourcelist(iens,isource)(i:) = ' '
                        do i=1,len_trim(sourcelist(iens,isource))
                            if ( sourcelist(iens,isource)(i:i).eq.' '
     +                           .or.
     +                           sourcelist(iens,isource)(i:i).eq.'+' )
     +                           sourcelist(iens,isource)(i:i) = '_'
                        enddo
                    enddo
                    if ( variable.eq.'ids' .and. nsource.lt.3 ) cycle
                    if ( oper.eq.'list' ) then
                        if ( nsource.eq.3 ) isource = 1
!                       find and print unique strings - a quadratic algorithm
                        do iens=nens1,nens2
                            do i=nens1,iens-1
                                k = 0
                                do j=isource,isource + nsource - 1
                                    if ( sourcelist(i,j)
     +                                   .eq.sourcelist(iens,j) )
     +                                   k = k + 1
                                enddo
                                if ( k.eq.nsource ) goto 701
                            enddo
                            print '(6a)',(trim(sourcelist(iens,i)),' ',
     +                           i=isource,isource+nsource-1)
 701                        continue
                        enddo
                        goto 999
                    elseif ( oper.eq.'conv' ) then
!                       find and print ensemble number corresponding to list
                        do iens=nens1,nens2
                            if ( lwrite ) then
                                print *,'comparing list  ',trim(list)
                                print *,'with sourcelist '
     +                               ,trim(sourcelist(iens,isource))
                            endif
                            if ( list.eq.'all' .or. 
     +                           list.eq.sourcelist(iens,isource) ) then
                                print '(i4)',iens
                            endif
                        enddo
                    else
                        write(0,*) 'convertmetadata: error: '//
     +                       'unknown operation ',oper
                        write(*,*) 'convertmetadata: error: '//
     +                       'unknown operation ',oper
                        call abort
                    endif
                elseif ( variable.eq.'reftime' ) then
*                   translate forecast reference time into dates
                    call getreftime(ncid,varid,tt,nt,firstmo,firstyr
     +                   ,nperyear,iperyear,lwrite)
                    if ( oper.eq.'list' ) then
!                       find and print unique strings
                        if ( nperyear.gt.12 ) then
                            write(0,*) 'convertmetadata: cannot handle'
     +                           //' nperyear = ',nperyear
                            write(*,*) 'convertmetadata: cannot handle'
     +                           //' nperyear = ',nperyear
                            call abort
                        endif
                        do i=1,12
                            do j=1,nperyear
                                if ( mod(firstmo+(j-1)*(12/nperyear),12)
     +                               .eq.mod(i,12) ) then
                                    print '(2a)','1',months(i)
                                    exit
                                endif
                            enddo
                        enddo
                        goto 999
                    elseif ( oper.eq.'conv' ) then
!                       find and print ensemble number corresponding to list
                        do fmo=1,12
                            if ( index(list,months(fmo)).ne.0 ) goto 723
                        enddo
                        write(0,*) 'convertmetadata: error: cannot '//
     +                       'locate month in ',trim(list)
                        write(*,*) 'convertmetadata: error: cannot '//
     +                       'locate month in ',trim(list)
                        call abort
 723                    continue
                        if ( lwrite ) print *,'searching for refmonths '
     +                       ,fmo
!                       set all time steps that do not match to undefined
                        if ( iperyear.eq.366 ) then
                            call getdymo(dy,mo,firstmo,nperyear)
                            jul0 = julday(mo,dy,firstyr)
                            tt1 = tt(1)
                            do it=1,nt
                                i = nint(tt(it) - tt1)
                                call caldat(jul0+i,mo,dy,yr)
                                if ( mo.ne.fmo ) then
                                    if ( lwrite ) print *,it,'mo = ',mo
     +                                   ,' set tt to undef',tt(it)
                                    tt(it) = 3e33
                                else
                                    if ( lwrite ) print *,it,'mo = ',mo
     +                                   ,' OK ',tt(it)
                                endif
                            enddo
                        elseif ( iperyear.le.12 .and. nperyear.le.12 )
     +                           then
                            do it=1,nt
                                i = nint((tt(it)-tt(1))
     +                               *nperyear/iperyear) + firstmo
                                if ( mod(i,nperyear).ne.
     +                               mod(fmo,nperyear) ) then
                                    tt(it) = 3e33
                                endif
                            enddo
                        elseif ( nperyear.gt.1 ) then
                            write(0,*) 'convertmetadata: error: '//
     +                           'cannot handle iperyear = '
     +                           ,iperyear,' yet'
                            write(*,*) 'convertmetadata: error: '//
     +                           'cannot handle iperyear = '
     +                           ,iperyear,' yet'
                            call abort
                        endif
!                       construct a set that describes the rest
!                       not using strides, as tehse do not work
                        foundvalid = .false.
                        invalid=0
                        do
!                           find valid entry
                            do valid=invalid+1,nt
                                if ( tt(valid).lt.1e33 ) exit
                            enddo
                            if ( valid.gt.nt ) then
                                if ( .not.foundvalid )
     +                               then
                                    write(0,*) 'metadataconvert: '//
     +                                   'error: no valid data'
                                    goto 999
                                endif
                                exit
                            endif
                            if ( lwrite ) print *,'valid data at ',valid
                            foundvalid = .true.
!                           find invalid
                            do invalid=valid+1,nt
                                if ( tt(invalid).gt.1e33 ) exit
                            enddo
                            if ( lwrite ) print *,'invalid data at '
     +                           ,invalid
!                           remember that the offset of ncks is 0-based
                            write(format,'(a,i1,a,i1,a)') '(i',1+
     +                           int(log10(max(1.,real(valid-1)))),
     +                           ',a,i',
     +                           1+int(log10(real(invalid-2))),
     +                           ')'
                            if ( lwrite ) print *,'format = '
     +                           ,trim(format)
                            print format,valid-1,',',invalid-2
                            if ( invalid.gt.nt ) exit
                        end do
                        goto 999
                    else
                        write(0,*) 'convertmetadata: error: '//
     +                       'unknown operation ',oper
                        write(*,*) 'convertmetadata: error: '//
     +                       'unknown operation ',oper
                        call abort
                    endif
                elseif ( variable.eq.'leadtime' ) then
*                   translate forecast reference time into dates
                    call getleadtime(ncid,varid,tt,nt,dimname,lwrite)
!                   convert to months
                    if ( dimname(1:4).eq.'days' .and. tt(2)-tt(1).ge.28
     +                   .and. tt(2)-tt(1).le.31 ) then
                        dimname = 'months'
                        n = 6
                        do it=1,nt
                            tt(it) = nint((tt(it)-15)/30)
                        enddo
                    endif
                    if ( oper.eq.'list' ) then
!                       find and print unique strings
!                       get unique time steps - quadratic algorithm, but that's OK
                        j = 1
                        ttu(1) = tt(1)
                        do it=2,nt
                            do i=1,j
                                if ( tt(it).eq.ttu(i) ) goto 711
                            enddo
                            j = j + 1
                            ttu(j) = tt(it)
 711                        continue
                        enddo
                        do i=1,j
                            print '(i3,a)',nint(ttu(i)),dimname(1:n)
                        enddo
                        goto 999
                    elseif ( oper.eq.'conv' ) then
!                       find and print ensemble number corresponding to list
                        write(0,*) 'convertmetadata: error: '//
     +                       'cannot handle convert yet'
                        call abort
                    else
                        write(0,*) 'convertmetadata: error: '//
     +                       'unknown operation ',oper
                        call abort
                    endif
                elseif ( variable.ne.'ids' ) then
                    write(*,*) 'convertmetadata: error: cannot handle '
     +                   ,trim(name),' yet'
                    write(0,*) 'convertmetadata: error: cannot handle '
     +                   ,trim(name),' yet'
                    call abort
                endif
            endif               ! variable-recognition case
            goto 899
 801        if ( lwrite ) 
     +           write(0,*) 'convertmetadata: warning: disregarding '//
     +           'variable with strange dimensions ',name(1:llen(name))
            ntvars = ntvars - 1
 899        continue
        enddo
        goto 999
* 
*       errors
*
 902    write(0,*) 'convertmetadata: error: found variable with '//
     +       'different dimensions ',name(1:llen(name))
        call abort
 999    continue
        end
