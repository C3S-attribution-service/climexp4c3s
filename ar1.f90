program ar1
!
!   generate a an AR(1) time series
!   X(mo,yr) = a1*X(mo-1,yr) + a2*X(mo,yr-1) + noise
!
    implicit none
    include 'param.inc'
    integer nseedmax
    integer yr,mo,nperyear,i,n,iarray(8)
    integer,allocatable :: iseed(:)
    real a1,a2,data(npermax,yrbeg:yrend),out(npermax,yrbeg:yrend)
    real ave(npermax),adev(npermax),sdev(npermax),var(npermax),skew(npermax),curt(npermax)
    real aa1(npermax),aa2(npermax),s1,s2
    character string*1023,varname*80,units*40
    logical lstandardunits,lwrite
    integer iargc
    real gasdev
    lwrite = .false.

    if ( iargc() > 4 .or. iargc() < 1 ) then
        print *,'Usage: ar1 file|nperyear a1 a2 [dummy]'
        stop
    endif
!
!   load random number seed
!   Due to a compiler bug in pgf90 this has to happen first
!
    call random_seed(size=n)
    allocate(iseed(n))
    call date_and_time(values=iarray)
    do i=7,1,-1
        iarray(i) = iarray(i)*iarray(i+1)
    enddo
    do i=1,n
        iseed(i) = iarray(8-mod(i,8))
    enddo
    call random_seed(put=iseed)
!
!   process arguments
!
    if ( iargc() >= 3 ) then
        call getarg(1,string)
        read(string,*,err=901) nperyear
        if ( nperyear < 1 .or. nperyear > npermax ) goto 901
        call getarg(2,string)
        read(string,*,err=902) a1
        if ( abs(a1) > 1 ) goto 902
        call getarg(3,string)
        read(string,*,err=903) a2
        if ( abs(a2) > 1 ) goto 903
        if ( abs(a1+a2) > 1 ) goto 904
        ave = 0
        var = 1
        aa1 = a1
        aa2 = a2
    else
        a1 = -999
        a2 = -999
        call getarg(1,string)
        lstandardunits = .true.
        call readseries(string,data,npermax,yrbeg,yrend,nperyear,varname,units,lstandardunits,lwrite)
        call detrend(data,npermax,yrbeg,yrend,yrbeg,yrend,1,nperyear,1)
        call seriesmoment(data,npermax,yrbeg,yrend,nperyear,yrbeg,yrend,ave,adev,sdev,var,skew,curt)
        call seriesautocor(data,npermax,yrbeg,yrend,nperyear,yrbeg,yrend,ave,var,aa1,aa2)
    end if

!   make correct autocorrelation structure with mean zero and sd 1

    do yr=yrbeg,yrend
        do mo=1,nperyear
            out(mo,yr) = gasdev(iseed)
            if ( mo > 1 .and. aa1(mo) < 1e33 ) then
                out(mo,yr) = (1-aa1(mo))*out(mo,yr) + aa1(mo)*out(mo-1,yr)
            elseif ( yr > yrbeg .and. aa1(mo) < 1e33 ) then
                out(mo,yr) = (1-aa1(mo))*out(mo,yr) + aa1(mo)*out(nperyear,yr-1)
            endif
            if ( yr > yrbeg .and. aa2(mO) < 1e33 ) then
                out(mo,yr) = (1-aa2(mo))*out(mo,yr) + aa2(mo)*out(mo,yr-1)
            endif
        enddo
    enddo
    do yr=yrbeg,yrend
        do mo=1,nperyear
            if ( aa1(mO) < 1e33 ) then
                out(mo,yr) = out(mo,yr)*sqrt(1-aa1(mo)**2)/(1-aa1(mo))
            end if
            if ( aa2(mO) < 1e33 ) then
                out(mo,yr) = out(mo,yr)*sqrt(1-aa2(mo)**2)/(1-aa2(mo))
            end if
        enddo
    enddo

!   normalise and copy undefs if required

    if ( a1 == -999 ) then
        do yr=yrbeg,yrend
            do mo=1,nperyear
                if ( data(mo,yr) < 1e33 ) then
                    out(mo,yr) = ave(mo) + sqrt(var(mo))*out(mo,yr)
                else
                    out(mo,yr) = 3e33
                end if
            enddo
        enddo
    end if

!   compute some statistics

    n = 0
    s1 = 0
    s2 = 0
    do yr=yrbeg,yrend
        do mo=1,nperyear
            if ( out(mo,yr) < 1e33 ) then
                n = n + 1
                s1 = s1 + out(mo,yr)
                s2 = s2 + out(mo,yr)**2
            end if
        enddo
    enddo
    write(6,'(a,100i16)') '# AR(1) noise time series, iseed = ',iseed
    if ( a1 /= -999 ) then
        write(6,'(a,f10.4)') '# a1 = ',a1
        write(6,'(a,f10.4)') '# a2 = ',a2
    else
        write(6,'(2a)') '# based on detrended file ',trim(string)
        write(6,'(6a)') '# ',trim(varname),' [',trim(units),'] noise'
        write(6,'(a,2000i10)')   '# day = ',(mo,mo=1,nperyear)
        write(6,'(a,2000f10.4)') '# ave = ',(ave(mo),mo=1,nperyear)
        write(6,'(a,2000f10.4)') '# sd  = ',(sqrt(var(mo)),mo=1,nperyear)
        write(6,'(a,2000f10.4)') '# a1  = ',(aa1(mo),mo=1,nperyear)
        write(6,'(a,2000f10.4)') '# a2  = ',(aa2(mo),mo=1,nperyear)
    end if
    write(6,'(a,f10.4)') '# mean = ',s1/n
    write(6,'(a,f10.4)') '# s.d. = ',sqrt(abs(s2/n - (s1/n)**2))
    call printdatfile(6,out,npermax,nperyear,yrbeg,yrend)
    goto 999
    
901 continue
    write(0,*) 'ar1: error: nperyear should be between 1 and ',npermax,', not ',trim(string)
    call exit(-1)
902 continue
    write(0,*) 'ar1: error: a1 should be between -1 and 1, not ',trim(string)
    call exit(-1)
903 continue
    write(0,*) 'ar1: error: a2 should be between -1 and 1, not ',trim(string)
    call exit(-1)
904 continue
    write(0,*) 'ar2: error: a1+a2 should be between -1 and 1, not ',a1+a2
    call exit(-1)
999 continue
end program
