program getunits
!
!       return the units a netcdf, grads or text file
!       for use in the web interface
!
    implicit none
    integer nvarmax
    parameter(nvarmax=100)
    integer iargc
    integer nx,ny,nz,nt,nperyear,nvars,i,j,nens,off,nargs,iarg,nt1
    logical xwrap,lwrite,lexist
    character file*255,var(nvarmax)*40,units(nvarmax)*60, &
 &       newunits(nvarmax)*60,lvar(nvarmax)*120,svar(nvarmax)*120,string*80,ensfile*255
    if ( iargc().lt.1 ) then
        print *,'usage: getunits file [file2 ..] [debug]'
        stop
    endif
    call getarg(1,file)
    lwrite = .false.
    nargs = iargc()
    if ( nargs.gt.1 ) then
        call getarg(nargs,string)
        if ( string.eq.'debug' .or. string.eq.'lwrite' ) then
            lwrite = .true.
            print *,'getunits: turned on debug printing'
            nargs = nargs - 1
        endif
    endif
    i = index(file,'%%')
    j = index(file,'++')
    if ( i.ne.0 .or. j.ne.0 ) then
        if ( lwrite ) print *,'counting ensemble members'
        off = 0
        do nens=0,999
            ensfile = file
            call filloutens(ensfile,nens)
            inquire(file=trim(ensfile),exist=lexist)
            if ( .not.lexist ) then
                if ( nens.eq.0 ) then
                    off = 1
                else
                    exit
                end if
            end if
        end do
        nens = nens - off
        call printassignment('NENS',nens)
    end if
    nt = 0
    do iarg=1,nargs
        call getarg(iarg,file)
        call getfileunits(file,nx,ny,nz,nt1,nperyear,nvarmax,nvars,var &
        &  ,units,newunits,lvar,svar,xwrap,lwrite)
        nt = nt + nt1
    end do
    call printassignment('NX',nx)
    call printassignment('NY',ny)
    call printassignment('NZ',nz)
    call printassignment('NT',nt)
    call printassignment('NPERYEAR',nperyear)
    if ( xwrap ) print '(a)','XWRAP=true'
    print '(3a)','VAR="',trim(var(1)),'"'
    print '(3a)','UNITS="',trim(units(1)),'"'
    print '(3a)','NEWUNITS="',trim(newunits(1)),'"'
    if ( lvar(1) /= ' '  ) print '(3a)','LVAR="',trim(lvar(1)),'"'
    if ( svar(1) /= ' '  ) print '(3a)','SVAR="',trim(svar(1)),'"'
end program

subroutine printassignment(name,value)
    implicit none
    integer value
    character name*(*)

    if ( value.lt.10 ) then
        print '(2a,i1)',trim(name),'=',value
    elseif ( value.lt.100 ) then
        print '(2a,i2)',trim(name),'=',value
    elseif ( value.lt.1000 ) then
        print '(2a,i3)',trim(name),'=',value
    elseif ( value.lt.10000 ) then
        print '(2a,i4)',trim(name),'=',value
    elseif ( value.lt.100000 ) then
        print '(2a,i5)',trim(name),'=',value
    else
        print '(2a,i6)',trim(name),'=',value
    endif
end subroutine