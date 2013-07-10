        subroutine diffit(data,npermax,nperyear,yrbeg,yrend,ndiff)
*
*       old entry point
*
        implicit none
        integer npermax,nperyear,yrbeg,yrend,ndiff
        real data(npermax,yrbeg:yrend)
        real minfacsum
        minfacsum = 0.75
        call ndiffit(data,npermax,nperyear,yrbeg,yrend,ndiff,minfacsum)
        end

        subroutine ndiffit(data,npermax,nperyear,yrbeg,yrend,ndiff
     +       ,minfacsum)
*
*       Take anomaly with respect to ndiff previous years, less when
*       undefined or near beginning of series (high-pass filter)
*       When ndiff < 0, take the sum of the n previous years
*       (poor man's low-pass filer).
*       7-dec-2001: only accept the point when at least 3/4 is valid data
*       29 november 2006: made this user-chosable.
*
        implicit none
        integer npermax,nperyear,yrbeg,yrend,ndiff
        real data(npermax,yrbeg:yrend),minfacsum
        integer i,j,k,n,m
        real absent,s
        parameter (absent=3e33)
        logical lwrite
        parameter(lwrite=.false.)
*
        if ( lwrite) then
            print *,'ndiffit: ndiff,minfacsum = ',ndiff,minfacsum
            print *,'         yrbeg,yrend  = ',yrbeg,yrend
        end if
        if ( minfacsum.lt.0 ) minfacsum = 0.75
        if ( ndiff.gt.0 ) then
            do i=yrend,yrbeg,-1
                do j=1,nperyear
                    if ( data(j,i).lt.0.9*absent ) then 
                        n = 0
                        s = 0
                        do k=max(i-abs(ndiff),yrbeg),i-1
                            if ( data(j,k).lt.0.9*absent ) then
                                s = s + data(j,k)
                                n = n + 1
                            endif
                        enddo
                        if ( n.ge.minfacsum*abs(ndiff) ) then
                            data(j,i) = data(j,i) - s/n
                        else
                            data(j,i) = absent
                        endif
                    else
                        data(j,i) = absent
                    endif
                enddo
            enddo
        elseif ( ndiff.lt.0 ) then
            do i=yrend,yrbeg,-1
                m = 0
                do j=1,nperyear
                    n = 0
                    s = 0
                    do k=max(i-abs(ndiff),yrbeg),i
                        if ( data(j,k).lt.0.9*absent ) then
                            s = s + data(j,k)
                            n = n + 1
                        endif
                    enddo
                    !!!print *,'ndiff: ',i,j,n
                    if ( n.ge.minfacsum*(1-ndiff) ) then
                        data(j,i) = s/n
                        m = m + 1
                    else
                        data(j,i) = absent
                    endif
                enddo
                if ( lwrite .and. m.gt.0 ) then
                    print *,i,(data(j,i),j=1,nperyear)
                end if
            enddo
        endif
*
        end

        subroutine dooverlap(fcst,npermax,yrbeg,yrend,nens1,nens2
     +       ,nperyear,ndiff)
        implicit none
        integer npermax,yrbeg,yrend,nens1,nens2,nperyear,ndiff
        real fcst(npermax,yrbeg:yrend,0:nens2)
        integer yr,mo,i,j,iens
        logical lvalid
!
!       make sure that the summed regions do not overlap
!
        do mo=1,nperyear
            yr = yrbeg
 101        continue
            do i=yr,yrend
!               search for a mo,yr with at least one valid point
                lvalid = .false.
                do iens=nens1,nens2
                    if ( fcst(mo,i,iens).lt.1e30 ) then
                        lvalid = .true.
                        exit
                    end if
                end do
                if ( lvalid .and. i.lt.yrend ) then
!                   put the next -ndiff years to undefined
                    do j=i+1,min(yrend,i-ndiff)
                        do iens=nens1,nens2
                            fcst(mo,j,iens) = 3e33
                        end do
                    end do
!                   and continue after that block
                    yr = j
                    goto 101
                end if
            end do
        end do
        end
