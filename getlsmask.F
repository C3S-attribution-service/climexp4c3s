        module lsdata
        integer :: nxls,nyls
        real,allocatable :: lsmask(:,:)
        end module lsdata

        subroutine getlsmask(startopts,lsmasktype,nxmax,xxls,nymax,yyls
     +       ,lwrite)
!
!       read lsmask if options (startopts,startopts+1,startopts+2)
!       are (lsm* file land|sea|notl|nots|all)
!
        use lsdata
        implicit none
#include "netcdf.inc"
#include "recfac.h"
        integer startopts,nxmax,nymax
        real xxls(nxmax),yyls(nymax)
        character lsmasktype*4
        logical lwrite
        integer ncid,nzls,nt,nperyear,firstyr,firstmo,nvars,jvars(6,1)
     +       ,ivars(2,1),endian,status,jx,jy
        real zzls(1),undef,scale,offset
        character string*10,file*255,datfile*255,title*255,vars(1)*20
     +       ,lvars(1)*80,units(1)*40
        integer,external :: get_endian

        call getarg(startopts,string)
        lsmasktype = 'all'
        if ( string(1:3).eq.'lsm' ) then
            call getarg(startopts+1,file)
            call getarg(startopts+2,lsmasktype)
            startopts = startopts + 3
            if ( lsmasktype.ne.'all' .and.
     +           lsmasktype.ne.'sea' .and.
     +           lsmasktype.ne.'nots' .and.
     +           lsmasktype.ne.'land' .and.
     +           lsmasktype.ne.'notl' .and.
     +           lsmasktype.ne.'5lan' .and.
     +           lsmasktype.ne.'5sea' ) then
                write(0,*) 'getlsmask: error: lsmasktype should be '//
     +               'all, land, sea, notl, nots, 5lan or 5sea '//
     +               'but I found ',lsmasktype
                write(*,*) 'getlsmask: error: lsmasktype should be '//
     +               'all, land, sea, notl, nots, 5lan or 5sea '//
     +               'but I found ',lsmasktype
            endif
            if ( lsmasktype.eq.'all' ) return
            print '(4a)','# using only ',lsmasktype,
     +           ' points from LS mask in ',trim(file)
            status = nf_open(file,nf_nowrite,ncid)
            if ( status.ne.nf_noerr ) then
                ncid = -1
                call parsectl(file,datfile,nxmax,nxls,xxls,nymax,nyls
     +               ,yyls,1,nzls,zzls,nt,nperyear,firstyr,firstmo,undef
     +               ,endian,title,1,nvars,vars,ivars,lvars,units)
            else
                if ( lwrite ) print *,'calling parsenc on ',trim(file)
                call parsenc(file,ncid,nxmax,nxls,xxls,nymax,nyls,yyls,1
     +               ,nzls,zzls,nt,nperyear,firstyr,firstmo,undef,title
     +               ,1,nvars,vars,jvars,lvars,units)
            endif
*           I do not (yet?) check that the grids are equal to the data
*           grids later on
            allocate(lsmask(nxls,nyls))
            if ( ncid.eq.-1 ) then
                open(1,file=datfile,form='unformatted',access='direct',
     +               recl=recfa4*nxls*nyls,status='old')
                read(1,rec=1) lsmask
                if ( endian*get_endian().eq.-1 ) then
                    call swapbyte4(lsmask,nxls*nyls)
                endif
            else
                call readonencfield(ncid,jvars,lsmask,nxls,nyls,lwrite)
            endif
            call checklsmask(lsmask,nxls,nyls,.false.)
        endif
        end
