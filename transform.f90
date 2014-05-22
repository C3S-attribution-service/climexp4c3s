program transform
!
!   tranformation program belonging to the KNMI'14 scenarios
!   Based on the R versions by Alexander Bakker
!
    implicit none
    integer npermax,yrbeg,yrend
    parameter(npermax=366,yrbeg=1800,yrend=2200)
    integer horizon,nperyear
    real series(npermax,yrbeg:yrend),newseries(npermax,yrbeg,yrend),deltas(5,12)
    character variable*10,scenario*2,region*3,infile*255,var*40,units*80
    logical lwrite
    integer iargc
    lwrite = .false.

    if ( iargc().lt.5 ) then
        write(0,*) 'usage: transform infile variable horizon scenario region'
        write(0,*) '       gives transformed time series following the KNMI''14 scenario'
        stop
    end if

    call get_parameters(infile,variable,horizon,scenario,region)
    call read_deltas(variable,horizon,scenario,region,deltas)
    call readseries(infile,series,npermax,yrbeg,yrend,nperyear, &
&       var,units,.false.,lwrite)
    if ( nperyear /= 366 ) then
        write(0,*) 'transform: error: only for daily data'
        write(*,*) 'transform: error: only for daily data'
        call abort
    end if 
    if ( variable == 'tg' .or. variable == 'tn' .or. variable == 'tx' ) then
        call lintransform(series,npermax,yrbeg,yrend,deltas,newseries)
    else if ( variable == 'rr' ) then
        call exptranform(series,npermax,yrbeg,yrend,deltas,newseries)
    else
        write(0,*) 'transform: error: unknow value for varibale, expecting ', &
&           'tg,tn,tx or rr, found ',variable
        call abort
    end if
    print '(2a)','# Transformed following the KNMI''14 scenario ',trim(scenario)
    print '(4a)','# for the region ',trim(region),' and for variable ',trim(variable)
    call copyheader(infile,6)
    call printdatfile(6,newseries,npermax,nperyear,yrbeg,yrend)

end program transform

subroutine get_parameter(infile,variable,horizon,scenario,region)
!
!   read parameters from the command line
!
    implicit none
    integer horizon
    character infile*(*), variable*(*), scenario*(*), region*(*)
    character string*80
    
    call getarg(1,infile)
    call getarg(2,variable)
    call getarg(3,string)
    read(string,*) horizon
    call getarg(4,scenario)
    call getarg(5,region)
    
end subroutine get_parameters

subroutine read_deltas(variable,horizon,scenario,region,deltas)
!
!   read the deltas from file for the years 2030, 2050 and 2085
!   and interpolate to the requested horizon
!
    implicit none
    integer nregions
    parameter(nregions=7,nyears=3)
    integer horizon,yrs(nyears)
    real deltas(5,12)
    character variable*(*),scenario*(*),region*(*)
    integer i,j,iyr,yr,iregion,mo
    character file*255,line*255,scen*2
    real delta0(5,12,nyears),a
    yrs = (/2030,2050,2085/)

    delta0 = 3e33
    do iyr=1,nyears
        yr = yrs(iyr)
        if ( yr.eq.2030 ) then
            scen = '__'
        else
            scen = scenario
        end if
        write(file,'(5a,i4,a)') 'KNMI14/deltas-KNMI14__',trim(variable),'__', &
&           trim(scen),'__',yr,'.txt'
        open(1,file=trim(file),status='old',err=900)
        read(1,'(a)') line
        do
            read(1,'(a)',end=100,err=901) line
            if ( line(1:3).eq.region ) then
                read(line(4:),*) mo
                read(line(4:),*) i,(delta0(j,mo,iyr),j=1,5)
            end if
        end do
100     continue
    end do iyr
    
    if ( horizon.lt.yrs(1) ) then
        write(0,*) 'get_deltas: error: horizon should be equal to or larger than ',yrs(1)
        write(*,*) 'transform: error: horizon should be equal to or larger than ',yrs(1)
        call abort
    end if
    if ( horizon.gt.yrs(nyears) ) then
        write(0,*) 'get_deltas: error: horizon should be equal to or smaller than ',yrs(nyears)
        write(*,*) 'transform: error: horizon should be equal to or smaller than ',yrs(nyears)
        call abort
    end if
    do iyr=1,nyears-1
        if ( horizon >= yrs(iyr) .and. horizon < yrs(ir+1) ) then
            ! linear interpolation
            a = real(horizon-yrs(iyr)/real(yrs(iyr+1)-yrs(iyr))
            do mo=1,12
                do j=1,5
                    deltas(j,mo) = (1-a)*delta0(j,mo,iyr) + a*delta0(j,mo,iyr+1)
                end do j
            end do mo
        end if
    end do iyr
    return

900 write(0,*) 'read_deltas: error: cannot locate file ',trim(file)
    call abort
    
901 write(0,*) 'read_deltas: error reading frome file ',trim(file)
    write(0,*) 'read_deltas: around line ',trim(line)
    call abort
    
end read_deltas

subroutine 