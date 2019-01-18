program filteryearseries

!   Apply hi- or lo-pass filters to a time series in standard
!   format and return the resulting series.

    implicit none
    include 'param.inc'
    include 'getopts.inc'
    integer :: i,j,yr,mo,nyr,nperyear,n,mens,mens1
    real,allocatable :: data(:,:)
    character :: file*1023,line*128,hilo*2,filtertype*12,var*40,units*60
    character :: lvar*120,svar*120,history*50000,metadata(2,100)*2000

    lwrite = .false. 
    if ( command_argument_count() < 4 ) then
        write(0,*) 'usage: filteryearseries hi|lo filtertype nyr file [minfac N]'
        call exit(-1)
    endif
    call get_command_argument(1,hilo)
    if ( hilo /= 'hi' .and. hilo /= 'lo' ) then
        write(0,*) 'filterseries: error: say hi or lo, not ',hilo
        call exit(-1)
    endif
    call get_command_argument(2,filtertype)
    call get_command_argument(3,line)
    read(line,*,err=901) nyr
    call get_command_argument(4,file)
    allocate(data(npermax,yrbeg:yrend))
    call readseriesmeta(file,data,npermax,yrbeg,yrend,nperyear,var,units,lvar,svar, &
        history,metadata,.false.,lwrite)
    call printvar(6,var,units,lvar)
    call printmetadata(6,file,' ',' ',history,metadata)
    write(6,'(3a,i4,3a)') '# operation :: each calendar month ',hilo, &
        '-pass filtered with a ',nyr,'-yr ',trim(filtertype),' filter'
    call getopts(5,command_argument_count()-1,nperyear,yrbeg,yrend, .true. ,mens1,mens)
    if ( minfac <= 0 ) minfac = 0.75

    if ( filtertype == 'running-mean' .or. filtertype(1:3) == 'box' ) then
        if ( hilo == 'hi' ) then
            call hipass(data,npermax,nperyear,yrbeg,yrend,nyr,minfac)
        else
            call ndiffit(data,npermax,nperyear,yrbeg,yrend,-nyr+1,minfacsum)
            call shiftseries(data,npermax,nperyear,yrbeg,yrend,-nperyear*((nyr-1)/2))
        endif
    else if ( filtertype == 'loess1' .or. filtertype == 'loess2' ) then
        call myloess(data,npermax,nperyear,yrbeg,yrend,nyr/2 &
            ,minfac,filtertype,hilo,'year','gaussian',lwrite)
    else
        write(0,*) 'filterseries: error: filtertype ',filtertype,' not yet implemented'
        call exit(-1)
    endif

    call printdatfile(6,data,npermax,nperyear,yrbeg,yrend)

    goto 999
901 write(0,*) 'filteryearseries: expcting an integer, not ',trim(line)
    call exit(-1)
902 write(0,*) 'filteryearseries: expcting a number, not ',trim(line)
    call exit(-1)
999 continue
end program filteryearseries
