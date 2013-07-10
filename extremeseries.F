        program extremeseries
*
*       compute extreme indices (ETCCDI, ECA&D) from daily data
*       input: daily time series
*       output: yearly/halfyearly/seasonal/monthly time series
*
        implicit none
#include "param.inc"
#include "getopts.inc"
        integer mpermax
        parameter(mpermax=366)
        integer nperyear,nperyearnew,yr,mo,dy,i,j,n,itype,nperyear2
        integer*2 nn(mpermax)
        real olddata(mpermax,yrbeg:yrend)
        real,allocatable :: newdata(:,:)
        character file*512,string*512,climdex*20,var*60,units*40,
     +       lvar*120,newunits*10
        integer iargc
        lwrite = .false.
        lstandardunits = .true.
*
        if ( iargc().lt.3 ) then
            print *,'usage: extremeseries infile nperyearnew'//
     +           ' extremeindex'
            stop
        endif
*
*       read data
*
        call getarg(1,file)
        call readseries(file,olddata,mpermax,yrbeg,yrend,nperyear,var
     +       ,units,lstandardunits,lwrite)
*
*       read operation
*
        call getarg(2,string)
        read(string,*,err=901) nperyearnew
        if ( abs(nperyearnew).gt.12 .or. nperyearnew.lt.1 ) then
            write(0,*) 'daily2longer: error: nperyearnew = ',nperyearnew
     +           ,' not yet supported'
            write(*,*) 'daily2longer: error: nperyearnew = ',nperyearnew
     +           ,' not yet supported'
            call abort
        endif
        call getarg(3,climdex)
        call getopts(4,iargc(),nperyear,yrbeg,yrend,.true.,0,0)
        allocate(newdata(nperyearnew,yrbeg:yrend))
*
*       perform operation
*
        call climexp2extreme(olddata,newdata,npermax,yrbeg,yrend,yr1,yr2
     +       ,nperyear,nperyearnew,var,units,newunits,climdex,lvar
     +       ,lwrite)
!
!       print out new series
!
        call copyheader(file,6)
        print '(6a)','# ',trim(climdex),' [',trim(newunits),'] ',
     +       trim(lvar)
        call printdatfile(6,newdata,abs(nperyearnew),abs(nperyearnew),
     +       yrbeg,yrend)
*
*       error messages
*
        goto 999
  901   write(0,*) 'daily2longer: expecting nperyearnew, not ',string
        call abort
  902   write(0,*) 'daily2longer: expecting value[%|p], not ',string
        call abort
  999   continue
        end
