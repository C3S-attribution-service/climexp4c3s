module lsdata
    integer :: nxls,nyls
    real,allocatable :: lsmask(:,:)
end module lsdata

subroutine getlsmask(startopts,lsmasktype,nxmax,xxls,nymax,yyls,lwrite)

!   read lsmask if options (startopts,startopts+1,startopts+2)
!   are (lsm* file land|sea|notl|nots|all)

    use lsdata
    implicit none
    include 'netcdf.inc'
    integer,parameter :: recfa4=4
    integer :: startopts,nxmax,nymax
    real :: xxls(nxmax),yyls(nymax)
    character :: lsmasktype*4
    logical :: lwrite
    integer :: ncid,nzls,nt,nperyear,firstyr,firstmo,nvars,jvars(6,1) &
        ,ivars(2,1),endian,status,jx,jy
    real :: zzls(1),undef,scale,offset
    character string*10,file*255,datfile*255,title*255,vars(1)*20 &
        ,lvars(1)*80,units(1)*40,longlsmasktype*100
    integer,external :: get_endian

    call get_command_argument(startopts,string)
    lsmasktype = 'all'
    if ( string(1:3) == 'lsm' ) then
        call get_command_argument(startopts+1,file)
        call get_command_argument(startopts+2,lsmasktype)
        startopts = startopts + 3
        if ( lsmasktype == 'all' ) return
        longlsmasktype = lsmasktype
        if ( lsmasktype == '5lan' ) then
            longlsmasktype = 'more than 50% land'
        else if ( lsmasktype == '5sea' ) then
            longlsmasktype = 'more than 50% sea'
        else if ( lsmasktype == 'land' ) then
            longlsmasktype = '100% land'
        else if ( lsmasktype == 'sea' ) then
            longlsmasktype = '100% sea'
        else if ( lsmasktype == 'notl' ) then
            longlsmasktype = 'not 100% land'
        else if ( lsmasktype == 'nots' ) then
            longlsmasktype = 'not 100% sea'
        else
            write(0,*) 'getlsmask: error: lsmasktype should be '// &
                'all, land, sea, notl, nots, 5lan or 5sea '// &
                'but I found ',lsmasktype
            write(*,*) 'getlsmask: error: lsmasktype should be '// &
                'all, land, sea, notl, nots, 5lan or 5sea '// &
                'but I found ',lsmasktype
        end if
        print '(4a)','# using only ',trim(longlsmasktype),' points from land/sea mask in ',trim(file)
        status = nf_open(file,nf_nowrite,ncid)
        if ( status /= nf_noerr ) then
            ncid = -1
            call parsectl(file,datfile,nxmax,nxls,xxls,nymax,nyls &
                ,yyls,1,nzls,zzls,nt,nperyear,firstyr,firstmo,undef &
                ,endian,title,1,nvars,vars,ivars,lvars,units)
        else
            if ( lwrite ) print *,'calling parsenc on ',trim(file)
            call parsenc(file,ncid,nxmax,nxls,xxls,nymax,nyls,yyls,1 &
                ,nzls,zzls,nt,nperyear,firstyr,firstmo,undef,title &
                ,1,nvars,vars,jvars,lvars,units)
        endif
!       I do not (yet?) check that the grids are equal to the data grids later on
        allocate(lsmask(nxls,nyls))
        if ( ncid == -1 ) then
            open(1,file=datfile,form='unformatted',access='direct', &
            recl=recfa4*nxls*nyls,status='old')
            read(1,rec=1) lsmask
            if ( endian*get_endian() == -1 ) then
                call swapbyte4(lsmask,nxls*nyls)
            endif
        else
            call readonencfield(ncid,jvars,lsmask,nxls,nyls,lwrite)
        endif
        call checklsmask(lsmask,nxls,nyls, .false. )
    endif
end subroutine getlsmask
