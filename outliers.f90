program outliers
!
!   detect outliers in a comparison of two time series
!
    implicit none
    include 'param.inc'
    integer dy,mo,yr,i,nperyear,nperyear2,n
    real xx(npermax,yrbeg:yrend),yy(npermax,yrbeg:yrend),sd
    logical lstandardunits,lwrite
    character var*40,units*40,xfile*255,yfile*255,string*20
    integer iargc
    
    if ( iargc() < 2 ) then
        write(0,*) 'usage: outliers xfile yfile [debig]'
        call exit(-1)
    end if
    call getarg(1,xfile)
    call getarg(2,yfile)
    lwrite = .false.
    if ( iargc() == 3 ) then
        call getarg(3,string)
        if ( string(1:5) == 'debug' .or. string (1:6) == 'lwrite' ) then
            lwrite = .true.
        end if
    end if
    lstandardunits = .true.
    call readseries(xfile,xx,npermax,yrbeg,yrend,nperyear,var,units,lstandardunits,lwrite)
    call readseries(yfile,yy,npermax,yrbeg,yrend,nperyear2,var,units,lstandardunits,lwrite)
    if ( nperyear /= nperyear2 ) then
        write(0,*) 'outliers: error: different time scales',nperyear,nperyear2
        call exit(-1)
    end if
    do yr=yrbeg,yrend
        do mo=1,npermax
            if ( xx(mo,yr) < 1e33 .and. yy(mo,yr) < 1e33 ) then
                xx(mo,yr) = xx(mo,yr) - yy(mo,yr) ! could include fit later on
            else
                xx(mo,yr) = 3e33
            end if
        end do
    end do
    call getseriesmoment(2,xx,npermax,yrbeg,yrend,nperyear,yrbeg,yrend,sd)
    do yr=yrbeg,yrend
        do mo=1,npermax
            if ( xx(mo,yr) < 1e33 .and. abs(xx(mo,yr)) > 5*sd ) then
                if ( nperyear.lt.360 ) then
                    print *,'outlier: ',yr,mo
                else
                    call getdymo(dy,i,mo,nperyear)
                    print *,'outlier: ',yr,i,dy
                end if
                print *,'outlier: ',xx(mo,yr),xx(mo,yr)+yy(mo,yr),yy(mo,yr)
            end if
        end do
    end do
 end program