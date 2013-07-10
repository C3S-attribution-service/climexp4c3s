        program ar1
*
*       generate a an AR(1) time series
*       X(mo,yr) = a1*X(mo-1,yr) + a2*X(mo,yr-1) + noise
*
        implicit none
#include "param.inc"
        integer nseedmax
        parameter (nseedmax=34)
        integer yr,mo,nperyear,i,n,iseed(nseedmax),iarray(8)
        real a1,a2,data(npermax,yrbeg:yrend),s1,s2
        character string*32
        integer iargc
        real gasdev
*
        if ( iargc().lt.3 ) then
            print *,'Usage: ar1 nperyear a1 a2'
            stop
        endif
        call getarg(1,string)
        read(string,*,err=901) nperyear
        if ( nperyear.lt.1 .or. nperyear.gt.npermax ) goto 901
        call getarg(2,string)
        read(string,*,err=902) a1
        if ( abs(a1).gt.1 ) goto 902
        call getarg(3,string)
        read(string,*,err=903) a2
        if ( abs(a1+a2).gt.1 ) goto 903
*
*       load random number seed
*
        call random_seed(size=n)
        if ( n.gt.nseedmax ) then
            write(0,*) 'ar1: error: n > 10 ',n,nseedmax
            write(*,*) 'ar1: error: n > 10 ',n,nseedmax
            call abort
        endif
        call date_and_time(values=iarray)
        do i=7,1,-1
            iarray(i) = iarray(i)*iarray(i+1)
        enddo
        do i=1,n
            iseed(i) = iarray(8-mod(i,8))
        enddo
        call random_seed(put=iseed)
*
        n = 0
        s1 = 0
        s2 = 0
        do yr=yrbeg,yrend
            do mo=1,nperyear
                data(mo,yr) = gasdev(iseed)
                if ( mo.gt.1 ) then
                    data(mo,yr) = data(mo,yr) + a1*data(mo-1,yr)
                elseif ( yr.gt.yrbeg ) then
                    data(mo,yr) = data(mo,yr) + a1*data(nperyear,yr-1)
                endif
                if ( yr.gt.yrbeg ) then
                    data(mo,yr) = data(mo,yr) + a2*data(mo,yr-1)
                endif
                n = n + 1
                s1 = s1 + data(mo,yr)
                s2 = s2 + data(mo,yr)**2
            enddo
        enddo
        write(6,'(a,100i16)') '# AR(1) noise time series, iseed = '
     +       ,iseed
        write(6,'(a,f10.4)') '# a1 = ',a1
        write(6,'(a,f10.4)') '# a2 = ',a2
        write(6,'(a,f10.4)') '# mean = ',s1/n
        write(6,'(a,f10.4)') '# var. = ',s2/n - (s1/n)**2
        call printdatfile(6,data,npermax,nperyear,yrbeg,yrend)
        goto 999
  901   continue
        print *,'error: nperyear should be between 1 and ',npermax
     +        ,', not ',string
        call abort
  902   continue
        print *,'error: a1 should be between -1 and 1, not ',string
        call abort
  903   continue
        print *,'error: a2 should be between -1 and 1, not ',string
        call abort
  999   continue
        end
