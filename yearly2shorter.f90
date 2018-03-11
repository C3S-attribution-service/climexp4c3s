program yearly2shorter
!
!   convert a yearly time series to a monthly one
!
    implicit none
    include 'param.inc'
    integer :: yr,mo,dy,nperyear,npernew,m1,n,k,i,nfac
    real :: data(npermax,yrbeg:yrend),newdata(npermax,yrbeg:yrend)
    character :: file*1023,var*80,units*120,lvar*120,svar*120,history*10000,metadata(2,100)*2000
    logical :: lwrite,lstandardunits
    integer :: iargc

    lwrite = .false.
    lstandardunits = .false.
    if ( iargc().lt.2 ) then
        print *,'usage: yearly2shorter infile.dat nperyearnew [mon n ave|sum m]'
        stop
    end if
    call getarg(1,file)
    call readseriesmeta(file,data,npermax,yrbeg,yrend,nperyear,var,units,lvar,svar, &
        history,metadata,lstandardunits,lwrite)
    call printvar(6,var,units,lvar)
    call copyheadermeta(file,6,' ',history,metadata)
    call getarg(2,file)
    read(file,*,err=901) npernew
    if ( iargc().gt.2 ) then
        call getarg(3,file)
        if ( file(1:3).ne.'mon' ) goto 902
        call getarg(4,file)
        read(file,*,err=903) m1
        call getarg(6,file)
        read(file,*,err=905) n
        call getarg(5,file)
        if ( file(1:3).eq.'ave' ) then
            nfac = 1
        else if ( file(1:3).eq.'sum' ) then
            nfac = n
        else
            goto 904
        end if
    else
        if ( nperyear /= 4 ) then
            m1 = 1
        else
            m1 = 0
        end if
        n = npernew
        nfac = 1
    end if
    call annual2shorter(data,npermax,yrbeg,yrend,nperyear, &
 &       newdata,npermax,yrbeg,yrend,npernew,m1,n,nfac,lwrite)
    call printdatfile(6,newdata,npermax,npernew,yrbeg,yrend)
    goto 999
901 write(0,*) 'yearly2shorter: error reading npernew from ',trim(file)
    call exit(-1)
902 write(0,*) 'yearly2shorter: error: expecting ''month'', not ',trim(file)
    call exit(-1)
903 write(0,*) 'yearly2shorter: error reading first month from ',trim(file)
    call exit(-1)
904 write(0,*) 'yearly2shorter: error: expecting ''ave|sum'', not ',trim(file)
    call exit(-1)
905 write(0,*) 'yearly2shorter: error reading number of months from ',trim(file)
    call exit(-1)
999 continue
end program
