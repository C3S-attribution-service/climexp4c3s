program seriesensanomal

!   compute the anomaly series or ensemble

    implicit none
    include 'param.inc'
    integer :: nperyear,mens1,mens,iens,nens1
    real,allocatable :: data(:,:,:)
    logical :: lwrite,lstandardunits,lexist
    character :: file*255,ensfile*255,var*40,units*40,line*255
    character :: lvar*120,svar*120,history*50000,metadata(2,100)*1000
    lwrite = .false.
    lstandardunits = .false. 

    if ( command_argument_count() /= 3 ) then
        print *,'usage: seriesensanomal ensmember ensfile dummy'
        print *,'       ensmember <0 denotes the full ensemble'
        call exit(-1)
    endif

    call get_command_argument(1,ensfile)
    read(ensfile,*) nens1
    call get_command_argument(2,ensfile)

    allocate(data(npermax,yrbeg:yrend,0:nensmax))
    call readensseriesmeta(ensfile,data,npermax,yrbeg,yrend,nensmax &
        ,nperyear,mens1,mens,var,units,lvar,svar,history,metadata &
        ,lstandardunits,lwrite)
    call anomalensemble(data,npermax,nperyear,yrbeg,yrend, &
        yrbeg,yrend,mens1,mens)
    write(6,'(a)') '# anomalies with respect to the ensemble mean'
    call printvar(6,var,units,lvar)
    if ( nens1 >= 0 ) then
        ! one ensemble member
        write(6,'(a,i4)') '# ensemble member ',nens1
        call filloutens(ensfile,nens1)
        call printmetadata(6,ensfile,' ',' ',history,metadata)
        call printdatfile(6,data(1,yrbeg,nens1),npermax,nperyear,yrbeg,yrend)
    else
        ! whole ensemble
        call printmetadata(6,ensfile,' ',' ',history,metadata)
        do iens=mens1,mens
            write(6,'(a,i4)') '# ensemble member ',iens
            write(6,'(a)')
            call printdatfile(6,data(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend)
        enddo
    end if
end program seriesensanomal