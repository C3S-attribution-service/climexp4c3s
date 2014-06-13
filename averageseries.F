        program addseries
*
*       add a few timeseries under the assumption that they represent the
*       same signal with independent noise and a station- and seasondependent
*       multplicative correction factor, which is assumed 1 for the first
*       series.
*
        implicit none
#include "param.inc"
        integer nmax,nyearmin,mpermax
        parameter(nmax=1000,nyearmin=10,mpermax=366)
        integer i,j,k,l,n,nseries,nperyear,nl,yr1,yr2
        real w
        real,allocatable :: data(:,:,:),means(:),weights(:,:)
        character string(nmax)*80,firstline*10000,var*40,units*20,
     +       weigh*10
        logical lwrite,lstandardunits,lweigh
        integer iargc,llen
        yr1=1850
        yr2=2100

        lwrite = .false.
        lstandardunits = .false.
        nseries = iargc() - 1 
        if ( nseries.lt.2 ) then
            write(0,*) 'usage: averageseries const|weight series1 '//
     +           '... seriesN'
            stop
        end if
        if ( nseries.gt.nmax ) then
            print *,'averageseries: error: recompile with nmax='
     +            ,nseries
            call abort
        endif
        allocate(data(mpermax,yr1:yr2,0:nseries))
        allocate(means(mpermax))
        allocate(weights(mpermax,nmax))
*
*       read data
*
        call getarg(0,firstline)
        nl = llen(firstline) + 1
        call getarg(1,weigh)
        if ( weigh(1:5).eq.'const' ) then
            lweigh = .false.
        else if ( weigh(1:5).eq.'weigh' ) then
            lweigh = .true.
        else
            write(0,*) 'error: expecting const|weigh, not ',weigh
            stop
        end if
        do k=1,nseries
            call getarg(k+1,string(k))
            if ( nl+1.lt.len(firstline) ) then
                firstline(nl+1:) = string(k)
                nl = llen(firstline) + 1
            endif
            write(0,*) 'reading file ',k,' ',trim(string(k))
            call readseries(string(k),data(1,yr1,k),mpermax,
     +           yr1,yr2,nperyear,var,units,lstandardunits,lwrite)
        enddo
*
*       compute means of first series (reference values)
*
        if ( lweigh ) then
            do j=1,nperyear
                n = 0
                means(j) = 0
                do i=yr1,yr2
                    if ( data(j,i,1).lt.1e33 ) then
                        n = n + 1
                        means(j) = means(j) + data(j,i,1)
                    endif
                enddo
                if ( n.ge.nyearmin ) then
                    means(j) = means(j)/n
                    weights(j,1) = 1
                else
                    print *,'error: not enough data for month ',j
     +                   ,' in reference series'
                    call abort
                endif
            enddo
            if ( nperyear.ge.360 ) then
                call averageseries_smooth(means,nperyear,5)
            end if
*           
*           compute weights
*           
            do k=2,nseries
                do j=1,nperyear
                    n = 0
                    weights(j,k) = 0
                    do i=yr1,yr2
                        if ( data(j,i,k).lt.1e33 ) then
                            n = n + 1
                            weights(j,k) = weights(j,k) + data(j,i,k)
                        endif
                    enddo
                    if ( n.ge.nyearmin ) then
                        weights(j,k) = means(j)/(weights(j,k)/n)
                    else
                        weights(j,k) = 3e33
                    endif
                enddo
                if ( nperyear.ge.360 ) then
                    ! smooth twice with a 5-day running mean filter
                    call averageseries_smooth(weights(1,k),nperyear,5)
                end if
            enddo
            open(1,file='weights.dat')
            write(1,'(a)') '# weights for the summation'
            do k=1,nseries
                write(1,'(a20,366f6.2)') string(k),
     +               (weights(j,k),j=1,nperyear)
            enddo
            close(1)
        else
            weights = 1
        end if
*
*       sum
*
        do i=yr1,yr2
            do j=1,nperyear
                data(j,i,0) = 0
                w = 0
                do k=1,nseries
                    if ( data(j,i,k).lt.1e33 .and. weights(j,k).lt.1e33
     +                    ) then
                        w = w + weights(j,k)
                        data(j,i,0) = data(j,i,0) +
     +                        weights(j,k)*data(j,i,k)
                    endif
                enddo
                if ( w.gt.0 ) then
                    data(j,i,0) = data(j,i,0)/w
                else
                    data(j,i,0) = 3e33
                endif
            enddo
        enddo
*
*       output
*
        print '(a)','# addseries: added the time series'
        print '(2a)','# ',trim(firstline)
        if ( lweigh ) print '(a)','# with weights in weights.dat'
        call copyheader(string(1),6)
        call printdatfile(6,data(1,yr1,0),mpermax,nperyear,
     +       yr1,yr2)
        end

        subroutine averageseries_smooth(weights,nperyear,nsmooth)
!
!       smooth weights by applying an  N-day runnng mean twice (i.e., a trinagle)
!
        implicit none
        integer nperyear,nsmooth
        real weights(nperyear),weights1(nperyear)
        real,allocatable :: array(:)

        allocate(array(nperyear))
        call averageseries_smooth1(weights,array,nperyear,nsmooth)
        call averageseries_smooth1(array,weights,nperyear,nsmooth)
        deallocate(array)
        end

        subroutine averageseries_smooth1(weights,array,nperyear,nsmooth)
!
!       put an 5-day smoothed version of weights into array
!
        implicit none
        integer nperyear,nsmooth
        real weights(nperyear),array(nperyear)
        integer j,k,l,n

        do j=1,nperyear
            n = 0
            array(j) = 0
            do k=-nsmooth/2,nsmooth/2
                l = j+k
                if ( l.lt.1 ) l = l + nperyear
                if ( l.gt.nperyear) l = l - nperyear
                if ( weights(l).lt.1e33 ) then
                    n = n + 1
                    array(j) = array(j) + weights(l)
                end if
            end do
            array(j) = array(j)/n
        end do
        end
