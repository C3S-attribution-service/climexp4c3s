program getunits
!
!       return the units a netcdf, grads or text file
!       for use in the web interface
!
    implicit none
    include 'params.h'
    integer,parameter :: nvarmax=100,mtmax=24*366*200
    integer :: nx,ny,nz,nt,nperyear,nvars,i,j,nens,off,nargs,iarg,nt1,ncid,firstmo,firstyr, &
        ivars(6,nvarmax),nens1,nens2,ndpm
    logical :: xwrap,lwrite,lexist,tdefined(mtmax),xrev,yrev
    real :: xx(nxmax),yy(nymax),zz(nzmax),undef,mean,offset,slope
    character :: file*255,var(nvarmax)*40,units(nvarmax)*60,title*500, &
        newunits(nvarmax)*60,lvar(nvarmax)*120,svar(nvarmax)*120,string*80,ensfile*255, &
        lz(3)*20,ltime*120,cell_methods(100)*100,history*50000,metadata(2,100)*1000
 
    if ( command_argument_count() < 1 ) then
        print *,'usage: getunits file [file2 ..] [debug]'
        stop
    endif
    call get_command_argument(1,file)
    lwrite = .false.
    nargs = command_argument_count()
    if ( nargs > 1 ) then
        call get_command_argument(nargs,string)
        if ( string == 'debug' .or. string == 'lwrite' ) then
            lwrite = .true.
            print *,'getunits: turned on debug printing'
            nargs = nargs - 1
        endif
    endif
    i = index(file,'%%')
    j = index(file,'++')
    if ( i /= 0 .or. j /= 0 ) then
        if ( lwrite ) print *,'counting ensemble members'
        off = 0
        do nens=0,999
            ensfile = file
            call filloutens(ensfile,nens)
            inquire(file=trim(ensfile),exist=lexist)
            if ( .not.lexist ) then
                if ( nens == 0 ) then
                    off = 1
                else
                    exit
                end if
            end if
        end do
        nens = nens - off
        call printassignment('NENS',nens)
    end if
    i = index(file,'@@')
    if ( i > 0 ) then
        ! an internal ensemble dimension, 
        ! only defined for .nc files, the .dat file has been derived from that
        i = index(file,'.dat')
        if ( i /= 0 ) then
            file(i:) = '.nc'
        end if
        ncid = 0
        call ensparsenc(file,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax &
            ,nz,zz,lz,nt,nperyear,firstyr,firstmo,ltime,tdefined,mtmax &
            ,nens1,nens2,undef,title,history,nvmax,nvars,var,ivars &
            ,lvar,svar,units,cell_methods,metadata)
        call printassignment('NENS',1+nens2)
        if ( units(1) == 'K' ) then
            mean = 273.13
        else
            mean = 0
        endif
        call makestandardunits(mean,nperyear,var,units(1) &
            ,newunits(1),offset,slope,ndpm,lwrite)
        call getxyprop(xx,nx,yy,ny,xrev,yrev,xwrap)
    else
        nt = 0
        do iarg=1,nargs
            call get_command_argument(iarg,file)
            call getfileunits(file,nx,ny,nz,nt1,nperyear,nvarmax,nvars,var &
               ,units,newunits,lvar,svar,xwrap,lwrite)
            nt = nt + nt1
        end do
    end if
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

    if ( value < 10 ) then
        print '(2a,i1)',trim(name),'=',value
    elseif ( value < 100 ) then
        print '(2a,i2)',trim(name),'=',value
    elseif ( value < 1000 ) then
        print '(2a,i3)',trim(name),'=',value
    elseif ( value < 10000 ) then
        print '(2a,i4)',trim(name),'=',value
    elseif ( value < 100000 ) then
        print '(2a,i5)',trim(name),'=',value
    else
        print '(2a,i6)',trim(name),'=',value
    endif
end subroutine