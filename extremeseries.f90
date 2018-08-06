program extremeseries
!
!   compute extreme indices (ETCCDI, ECA&D) from daily data
!   input: daily time series
!   output: yearly/halfyearly/seasonal/monthly time series
!
    implicit none
    include 'param.inc'
    include 'getopts.inc'
    integer mpermax
    parameter(mpermax=366)
    integer nperyear,nperyearnew,yr,mo,dy,i,j,n,itype,nperyear2
    integer*2 nn(mpermax)
    real olddata(mpermax,yrbeg:yrend)
    real,allocatable :: newdata(:,:)
    character file*512,string*512,climdex*20,var*60,units*40,lvar*120,newunits*10

    lwrite = .false.
    lstandardunits = .true.
!
    if ( command_argument_count().lt.3 ) then
        print *,'usage: extremeseries infile nperyearnew extremeindex [gt percentile%]'
        call exit(-1)
    endif
!
!   read data
!
    call get_command_argument(1,file)
    call readseries(file,olddata,mpermax,yrbeg,yrend,nperyear,var,units,lstandardunits,lwrite)
!
!   read operation
!
    call get_command_argument(2,string)
    read(string,*,err=901) nperyearnew
    if ( abs(nperyearnew).gt.12 .or. nperyearnew.lt.1 ) then
        write(0,*) 'extremeseries: error: nperyearnew = ',nperyearnew,' not yet supported'
        write(*,*) 'extremeseries: error: nperyearnew = ',nperyearnew,' not yet supported'
        call abort
    endif
    call get_command_argument(3,climdex)
    call getopts(4,command_argument_count(),nperyear,yrbeg,yrend,.true.,0,0)
    allocate(newdata(nperyearnew,yrbeg:yrend))
!
!   perform operation
!
    call climexp2extreme(olddata,newdata,npermax,yrbeg,yrend,yr1,yr2  &
 &       ,nperyear,nperyearnew,minindx,var,units,newunits,climdex     &
 &       ,lvar,lwrite)
!
!   print out new series
!
    call copyheader(file,6)
    print '(6a)','# ',trim(climdex),' [',trim(newunits),'] ',trim(lvar)
    if ( climdex.eq.'Rnnmm' ) then
        print '(a,f6.1,a)','# using threshold ',minindx,' mm/dy'
    end if
    call printdatfile(6,newdata,abs(nperyearnew),abs(nperyearnew),yrbeg,yrend)
!
!   error messages
!
    goto 999
901 write(0,*) 'extremeseries: expecting nperyearnew, not ',string
    call abort
902 write(0,*) 'extremeseries: expecting value[%|p], not ',string
    call abort
999 continue
end program