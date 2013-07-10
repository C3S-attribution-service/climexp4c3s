        program average_ensemble
!
!       compute the ensemble average of a set of time series
!
        implicit none
#include "param.inc"
        integer i,mo,yr,iens,mens1,mens,nperyear,iarg
        integer nn(1:npermax,yrbeg:yrend),nens1,nens2
        real data(1:npermax,yrbeg:yrend,0:nensmax),
     +       mean(1:npermax,yrbeg:yrend)
        logical lstandardunits,lwrite
        character file*255,var*20,units*20,string*80,oper*3
        integer iargc
        
        if ( iargc().lt.2 ) then
            print *,'usage: average_ensemble ensfile [ens n1 n2] dummy'
            call abort
        endif

        lstandardunits = .false.
        lwrite = .false.
        call getarg(1,file)
        nens1 = 0
        nens2 = 999
        oper = 'mea'
        iarg = 2
        do while ( iarg.le.iargc() )
            call getarg(iarg,string)
            if ( string.eq.'ens' ) then
                call getarg(iarg+1,string)
                read(string,*) nens1
                call getarg(iarg+2,string)
                read(string,*) nens2
                iarg = iarg + 3
            else if ( string.eq.'mean' ) then
                oper = 'mea'
                iarg = iarg + 1
            else if ( string.eq.'min' ) then
                oper = 'min'
                iarg = iarg + 1
            else if ( string.eq.'max' ) then
                oper = 'max'
                iarg = iarg + 1
            else if ( string.eq.'debug' .or. string.eq.'lwrite' ) then
                lwrite = .true.
                print *,'turned on debug output'
                iarg = iarg + 1
            else
                iarg = iarg + 1
            end if
        end do
        call readensseries(file,data,npermax,yrbeg,yrend,nensmax
     +        ,nperyear,mens1,mens,var,units,lstandardunits,lwrite)
        nens1 = max(nens1,mens1)
        nens2 = min(nens2,mens)

        nn(1:nperyear,yrbeg:yrend) = 0
        if ( oper.eq.'mea' ) then
            mean(1:nperyear,yrbeg:yrend) = 0
        else if ( oper.eq.'min' ) then
            mean(1:nperyear,yrbeg:yrend) = 3e33
        else if ( oper.eq.'max' ) then
            mean(1:nperyear,yrbeg:yrend) = -3e33
        else
            write(0,*) 'average_ensemble: error: unknown operation '
     +           ,oper
            call abort
        end if
        do iens=nens1,nens2
            do yr=yrbeg,yrend
                do mo=1,nperyear
                    if ( data(mo,yr,iens).lt.1e33 ) then
                        nn(mo,yr) = nn(mo,yr) + 1
                        if ( oper.eq.'mea' ) then
                            mean(mo,yr) = mean(mo,yr) + data(mo,yr,iens)
                        else if ( oper.eq.'min' ) then
                            mean(mo,yr) =
     +                           min(mean(mo,yr),data(mo,yr,iens))
                        else if ( oper.eq.'max' ) then
                            mean(mo,yr) =
     +                           max(mean(mo,yr),data(mo,yr,iens))
                        else
                            write(0,*) 'average_ensemble: error: '//
     +                           'unknown operation 2 ',oper
                            call abort
                        end if
                    endif
                enddo
            enddo
        enddo
        do yr=yrbeg,yrend
            do mo=1,nperyear
                if ( nn(mo,yr).gt.1 ) then
                    if ( oper.eq.'mea' ) then
                        mean(mo,yr) = mean(mo,yr)/nn(mo,yr)
                    end if
                else
                    mean(mo,yr) = 3e33
                endif
            enddo
        enddo

        print '(4a)','# ensemble ',oper,' of file ',trim(file)
        print '(a,i3,a,i3)','# using ensemble members ',nens1,' to '
     +       ,nens2
        call filloutens(file,mens1)
        call copyheader(file,6)
        call printdatfile(6,mean,npermax,nperyear,yrbeg,yrend)

        end
