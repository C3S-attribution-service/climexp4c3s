        subroutine whencomposite(data,mpermax,yrbeg,yrend,nensmax
     +       ,nens2series,lnewyr,nperyear)
!
!       determine which time slices to use for the composite
!       the boudnaries in data are given in maxdata/mindata
!       the timeslices are marked by replacing the value of data
!       by -1,0,+1
!
        implicit none
#include "getopts.inc"
        integer nens2series,mpermax,yrbeg,yrend,nensmax,nperyear
        real data(mpermax,yrbeg:yrend,0:nensmax)
        logical lnewyr(mpermax,yrbeg:yrend)
        integer iens,yr,mo,mm,mmold,moold,yrold,d,m
        real ddata,ddataold
        logical lfirst
        integer,save :: init
        data init /0/
!
        if ( init.eq.0 ) then
            init = 1
            lnewyr = .true.
            write(0,'(a)') 'For composite taking dates: '
        endif
!
        mmold = -1
        lfirst = .true.
        do yr=yrbeg,yrend
            do mo=1,nperyear
                do iens=nens1,nens2series
                    ddata = data(mo,yr,iens)
                    if ( data(mo,yr,iens).gt.1e33 ) then
*                       do nothing
                    elseif ( maxdata.lt.1e33 ) then
                        if ( mindata.gt.-1e33 ) then
                            if ( maxdata.gt.mindata ) then
*                               interval composite
                                if ( data(mo,yr,iens).gt.maxdata.or.
     +                               data(mo,yr,iens).lt.mindata )
     +                               then
                                    data(mo,yr,iens) = 0
                                else
                                    data(mo,yr,iens) = 1
                                endif
                            else
*                               signed composite
                                if ( data(mo,yr,iens).gt.maxdata )
     +                               then
                                    data(mo,yr,iens) = +1
                                elseif ( data(mo,yr,iens).lt.mindata
     +                                   ) then
                                    data(mo,yr,iens) = -1
                                else
                                    data(mo,yr,iens) = 0
                                endif
                            endif
                        else
*                           unsigned composite
                            if ( data(mo,yr,iens).lt.maxdata ) then
                                data(mo,yr,iens) = +1
                            else
                                data(mo,yr,iens) = 0
                            endif
                        endif
                    elseif ( mindata.gt.-1e33 ) then
                        if ( data(mo,yr,iens).gt.mindata ) then
                            data(mo,yr,iens) = +1
                        else
                            data(mo,yr,iens) = 0
                        endif
                    endif
                    if ( data(mo,yr,iens).lt.1e33 .and. 
     +                   data(mo,yr,iens).ne.0 .and. lnewyr(mo,yr) )
     +                   then
!
!                       print times used for composite
!
                        lnewyr(mo,yr) = .false.
                        mm = nperyear*(yr-yrbeg) + mo
                        if ( nperyear.eq.1 ) then
                            write(0,'(i4,f8.3,a)') yr,ddata,' '
                        elseif ( nperyear.le.12 ) then
                            if ( mmold.gt.0 .and. mm.ne.mmold+1 ) then
                                if ( .not.lfirst ) then
                                    yrold = yrbeg + (mmold-1)/nperyear
                                    moold = 1 + mod(mmold-1,nperyear)
                                    write(0,'(a)') '-'
                                    call printmoyrval(moold,yrold
     +                                   ,ddataold)
                                endif
                                mmold = -1
                            endif
                            if ( mmold.eq.-1 ) then
                                call printmoyrval(mo,yr,ddata)
                                lfirst = .true.
                            else
                                lfirst = .false.
                            endif
                            mmold = mm
                            ddataold = ddata
                        else
                            if ( mmold.gt.0 .and. mm.ne.mmold+1 ) then
                                if ( .not.lfirst ) then
                                    yrold = yrbeg + (mmold-1)/nperyear
                                    moold = 1 + mod(mmold-1,nperyear)
                                    write(0,'(a)') '-'
                                    call printdymoyrval(moold,yrold
     +                                   ,ddataold,nperyear)
                                endif
                                mmold = -1
                            endif
                            if ( mmold.eq.-1 ) then
                                call printdymoyrval(mo,yr,ddata,nperyear
     +                               )
                                lfirst = .true.
                            else
                                lfirst = .false.
                            endif
                            mmold = mm
                            ddataold = ddata
                        endif
                    endif
                enddo
            enddo
        enddo
!
!       last bit
!
        if ( nperyear.eq.1 ) then
!           no intervals
        else if (  nperyear.le.12 ) then
            if ( .not.lfirst ) then
                yrold = yrbeg + (mmold-1)/nperyear
                moold = 1 + mod(mmold-1,nperyear)
                write(0,'(a)') '-'
                call printmoyrval(moold,yrold,ddataold)
            endif
        else 
            if ( .not.lfirst ) then
                yrold = yrbeg + (mmold-1)/nperyear
                moold = 1 + mod(mmold-1,nperyear)
                write(0,'(a)') '-'
                call printdymoyrval(moold,yrold,ddataold,nperyear)
            endif
        endif

*       avoid cutting off ...
        mindata = -3e33
        maxdata = +3e33
        end

        subroutine printmoyrval(mo,yr,ddata)
        implicit none
        integer mo,yr
        real ddata
        write(0,'(i2,a,i4,a,f8.3,a)') mo,'-',yr,' (',ddata,') '
        end

        subroutine printdymoyrval(nr,yr,ddata,nperyear)
        implicit none
        integer nr,yr,nperyear
        real ddata
        integer d,m
        character months(12)*3
        data months /'jan','feb','mar','apr','may','jun',
     +               'jul','aug','sep','oct','nov','dec'/
        call getdymo(d,m,nr,nperyear)
        write(0,'(i2,a,i4,a,f8.3,a)') d,months(m),yr,' (',ddata,') '
        end

        subroutine makecomposite(ddata,dindx,n,b,db,prob,df,lwrite)
!
!       compute a composite of dindx based on the values of
!       ddata (-1,0,+1).  The number of degrees of freedom is input,
!       this may be different from n-1 because of serial correlations.
!       Output: composite b, error on composite db, p-value prob
!       I assume that dindx is already an anomaly field, so I compute
!       the mean of ddata*dindx
!       Destroys ddata
!
        implicit none
        integer n
        real ddata(n),dindx(n),b,db,prob,df
        logical lwrite
        integer i,n1,n0
        real s1,s2
        real,external :: erfcc

        s1 = 0
        do i=1,n
            s1 = s1 + dindx(i)
        end do
        s1 = s1/n
        if ( lwrite) print *,'makecomposite: mean was : ',s1,n
        dindx = dindx - s1

        n0 = 0
        n1 = 0
        do i=1,n
            if ( abs(ddata(i)).ne.1 .and. ddata(i).ne.0 ) then
                write(0,*)
     +               'makecomposite: error: expecting -1,0,+1, not '
     +               ,ddata(i)
            end if
            if ( ddata(i).ne.0 ) then
                n1 = n1 + 1
                ddata(n1) = ddata(i)*dindx(i)
            else
                n0 = n0 + 1
            endif
        end do
        call getmoment(1,ddata,n1,s1)
        ddata = ddata - s1
        call getmoment(2,ddata,n1,s2)
        b = s1
        if ( n0.eq.0 .or. n1.eq.0 .or. df.le.0 ) then
            db = 3e33
            prob = 3e33
        else
            db = sqrt(s2*(1./n1+1./n0)*(n-2)/df)
!           assume Gaussian for prob
            if ( db.gt.0 ) then
                prob = erfcc(abs(b/db)/sqrt(2.))
            else
                prob = 3e33
            end if
        end if
        if ( lwrite ) then
            print *,'makecomposite: composite: ',b,n1
            print *,'makecomposite: errorcomp: ',db
            print *,'makecomposite: p-value:   ',prob
        end if
        end
