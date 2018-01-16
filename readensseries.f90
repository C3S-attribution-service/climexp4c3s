subroutine readensseries(file,data,npermax,yrbeg,yrend,nensmax &
    ,nperyear,mens1,mens,var,units,lstandardunits,lwrite)
    implicit none
    integer :: npermax,yrbeg,yrend,nperyear,nensmax,mens1,mens
    real :: data(npermax,yrbeg:yrend,0:nensmax)
    logical :: lstandardunits,lwrite
    character :: file*(*),var*(*),units*(*)
    character ::  lvar*120,svar*120,history*50000,metadata(2,100)*1000
    
    call readensseriesmeta(file,data,npermax,yrbeg,yrend,nensmax &
        ,nperyear,mens1,mens,var,units,lvar,svar,history,metadata,lstandardunits,lwrite)

end subroutine readensseries

subroutine readensseriesmeta(file,data,npermax,yrbeg,yrend,nensmax &
    ,nperyear,mens1,mens,var,units,lvar,svar,history,metadata,lstandardunits,lwrite)

!   read a data file or an ensemble of data files in the array data
!   also returns the number of ensemble members in mens, 0=not an ensemble

    implicit none
    integer :: npermax,yrbeg,yrend,nperyear,nensmax,mens1,mens
    real :: data(npermax,yrbeg:yrend,0:nensmax)
    logical :: lstandardunits,lwrite
    character :: file*(*),var*(*),units*(*)
    character ::  lvar*(*),svar*(*),history*(*),metadata(2,100)*(*)
    integer :: iens
    logical :: lexist,lfirst
    character :: ensfile*255,line*10

    if ( lwrite ) then
        print *,'readensseries: file = ',trim(file)
        print *,'               npermax,yrbeg,yrend = ',npermax &
        ,yrbeg,yrend
    endif

!   ensembles are denoted by a '%' or '++' in the file name

    if ( index(file,'%') == 0 .AND. index(file,'++') == 0 ) then
        mens = 0
        mens1 = 0
        if ( lwrite ) print *,'not an ensemble, calling readseries'
        call readseriesmeta(file,data(1,yrbeg,0),npermax,yrbeg,yrend &
            ,nperyear,var,units,lvar,svar,history,metadata,lstandardunits,lwrite)
    else
        mens = -1
        mens1 = 0
        lfirst = .true.
        do iens=0,nensmax
            ensfile = file
            call filloutens(ensfile,iens)
            !!!if ( lwrite) print *,'looking for file ',trim(ensfile)
            inquire(file=ensfile,exist=lexist)
            if ( .NOT. lexist ) then
                if ( mens == -1 ) then
                    mens1 = iens + 1
                end if
            else
                mens = iens
            end if
        end do
        if ( mens < 0 ) then
            write(0,*) 'readensseries: error: could not find ensemble ' &
                ,trim(file),trim(ensfile)
            call abort
        end if
        do iens=mens1,mens
            ensfile = file
            call filloutens(ensfile,iens)
            inquire(file=ensfile,exist=lexist)
            if ( .not. lexist ) then
                data(:,:,iens) = 3e33
            else
                if ( lwrite) print *,'reading file ',trim(ensfile)
                call keepalive1('Reading ensemble member',iens,mens)
                if ( lfirst ) then
                    lfirst = .false.
                    call readseriesmeta(ensfile,data(1,yrbeg,iens),npermax &
                        ,yrbeg,yrend,nperyear,var,units,lstandardunits &
                        ,lwrite)
                else
                    call readseries(ensfile,data(1,yrbeg,iens),npermax &
                        ,yrbeg,yrend,nperyear,var,units,lvar,svar,history &
                        ,metadata,lstandardunits,lwrite)
                end if
            end if
        end do
    end if
end subroutine readensseriesmeta
