        subroutine shiftseries(data,npermax,nperyear,yrbeg,yrend,nshift)
*       
*       shift the time series over nshift units (normally months)
*
        implicit none
        integer npermax,nperyear,yrbeg,yrend,nshift
        real data(npermax,yrbeg:yrend)
        integer i,j,ii,k
*       
        if ( nshift.gt.0 ) then
            do i=yrend,yrbeg,-1
                do j=nperyear,1,-1
                    k = j + nshift
                    call normon(k,i,ii,nperyear)
                    if ( ii.le.yrend ) then
                        data(k,ii) = data(j,i)
                        data(j,i) = 3e33
                    endif
                enddo
            enddo
        elseif ( nshift.lt.0 ) then
            do i=yrbeg,yrend
                do j=1,nperyear
                    k = j + nshift
                    call normon(k,i,ii,nperyear)
                    if ( ii.ge.yrbeg ) then
                        data(k,ii) = data(j,i)
                        data(j,i) = 3e33
                    endif
                enddo
            enddo
        endif
        end

        subroutine shiftseriesyear(series,
     +      npermax,nperyear,yrbeg,yrend,dyr,lwrite)
!
!       shift the dates forward by dyr years
!
        implicit none
        integer npermax,nperyear,yrbeg,yrend,dyr
        real series(npermax,yrbeg:yrend)
        logical lwrite
        integer j,dy,yr,d,dtemp
        logical allundef
        integer leap
        external leap
!
        dtemp = 0
        if ( dyr.eq.0 ) then
            return
        else if ( dyr.gt.0 ) then
            if ( nperyear.ne.366 ) then
                do yr=yrend-dyr,yrbeg,-1
                    do dy=1,nperyear
                        series(dy,yr+dyr) = series(dy,yr)
                        series(dy,yr) = 3e33
                    end do
                end do
            else
                ! leap years...
                d = 0
                allundef = .true.
                do yr=yrend-dyr,yrbeg,-1
                    do dy=nperyear,1,-1
                        ! start with a straight copy
                        if ( allundef .and. series(dy,yr).lt.1e33 ) then
                            allundef = .false.
                            d = 0
                        end if
                        if ( dy.eq.60 ) then
                            if ( leap(yr).eq.2 .and. leap(yr+dyr).ne.2 ) 
     +                           then
                                d = d - 1
                                series(dy,yr+dyr) = 3e33
                                if ( lwrite .and. 
     +                               series(dy-3,yr).lt.1e33 ) then
                                    print *,'different leap years, '//
     +                                   'changing d to ',d
                                    print *,'series(',dy,yr+dyr,
     +                                   ') = 3e33'
                                end if
                                if ( d.eq.0 ) dtemp = 1
                                if ( d.lt.-1 .and.
     +                               series(dy-3,yr).lt.1e33 ) then
                                    write(0,*) 'shiftseriesyear: unte'//
     +                                   'sted case of two leap years',d
                                end if
                            end if
                            if ( leap(yr).ne.2 .and. leap(yr+dyr).eq.2 )
     +                           then
                                d = d + 1
                                if ( lwrite .and.
     +                               series(dy-3,yr).lt.1e33 ) then
                                    print *,'different leap years, '//
     +                                   'changing d to ',d
                                end if
                                if ( d.gt.1 .and.
     +                               series(dy-3,yr).lt.1e33 ) then
                                    write(0,*) 'shiftseriesyear: unte'//
     +                                   'sted case of two leap years',d
                                end if
                                cycle
                            end if
                        end if
                        if ( dy+d.lt.1 ) then
                            j = dy+d+nperyear
                            series(j,yr+dyr-1)=series(dy,yr)
                            if ( lwrite .and. series(dy,yr).lt.1e33 ) 
     +                           print *,'copying from ',dy,yr,
     +                           ' to ',j,yr+dyr-1
                        else if ( dy+d.gt.nperyear ) then
                            j = dy+d-nperyear
                            series(j,yr+dyr+1)=series(dy,yr)
                            if ( lwrite .and. series(dy,yr).lt.1e33 ) 
     +                           print *,'copying from ',dy,yr,
     +                           ' to ',j,yr+dyr+1
                        else
                            j = dy+d+dtemp
                            if ( j.lt.1 .or. j.gt.366 ) then
                                write(0,*) 'shiftseries: j out of '//
     +                               'bounds: ',j,dy,d,dtemp
                            end if
                            if ( j.eq.60 .and. dy.ne.60 
     +                           .and. leap(yr+dyr).eq.1 ) then
                                series(j,yr+dyr) = series(dy-1,yr)
                            else
                                series(j,yr+dyr) = series(dy,yr)
                            end if
                            if ( j.lt.3 .or. j.gt.363 .or. 
     +                           abs(j-60).lt.3 ) then
                                if ( lwrite .and. series(dy,yr).lt.1e33)
     +                               print *,'copying from ',
     +                               dy,yr,' to ',j,yr+dyr
                            end if
                        end if
                        dtemp = 0
                        series(dy,yr) = 3e33
                    end do
                end do
            end if
        else if ( dyr.lt.0 ) then
            if ( nperyear.ne.366 ) then
                do yr=yrbeg-dyr,yrend
                    do dy=1,nperyear
                        series(dy,yr+dyr) = series(dy,yr)
                        series(dy,yr) = 3e33
                    end do
                end do
            else
                ! leap years...
                d = 0
                allundef = .true.
                do yr=yrbeg-dyr,yrend
                    do dy=1,nperyear
                        ! start with a straight copy
                        if ( allundef .and. series(dy,yr).lt.1e33 ) then
                        allundef = .false.
                            d = 0
                        end if
                        if ( dy.eq.60 ) then
                            if ( leap(yr).eq.2 .and. leap(yr+dyr).ne.2 ) 
     +                           then
                                d = d - 1
                                series(dy,yr+dyr) = 3e33
                                if ( lwrite .and. 
     +                               series(dy+3,yr).lt.1e33 ) then
                                    print *,'different leap years, '//
     +                                   'changing d to ',d
                                    print *,'series(',dy,yr+dyr,
     +                                   ') = 3e33'
                                end if
                                if ( d.eq.0 ) dtemp = 1
                            end if
                            if ( leap(yr).ne.2 .and. leap(yr+dyr).eq.2 )
     +                           then
                                d = d + 1
                                if ( lwrite .and.
     +                               series(dy+3,yr).lt.1e33 ) then
                                    print *,'different leap years, '//
     +                                   'changing d to ',d
                                end if
                                cycle
                            end if
                        end if
                        if ( dy-d.lt.1 ) then
                            j = dy-d+nperyear
                            series(j,yr+dyr-1)=series(dy,yr)
                            if ( lwrite .and. series(dy,yr).lt.1e33 ) 
     +                           print *,'copying from ',dy,yr,
     +                           ' to ',j,yr+dyr-1
                        else if ( dy-d.gt.nperyear ) then
                            j = dy-d-nperyear
                            series(j,yr+dyr+1)=series(dy,yr)
                            if ( lwrite .and. series(dy,yr).lt.1e33 ) 
     +                           print *,'copying from ',dy,yr,
     +                           ' to ',j,yr+dyr+1
                        else
                            j = dy-d-dtemp
                            if ( j.eq.60 .and. dy.ne.60
     +                           .and. leap(yr+dyr).eq.1 ) then
                                series(j,yr+dyr) = series(dy-1,yr)
                            else
                                series(j,yr+dyr) = series(dy,yr)
                            end if
                            if ( j.lt.3 .or. j.gt.363 .or. 
     +                           abs(j-60).lt.3 ) then
                                if ( lwrite .and. series(dy,yr).lt.1e33)
     +                               print *,'copying from ',
     +                               dy,yr,' to ',j,yr+dyr
                            end if
                        end if
                        dtemp = 0
                        series(dy,yr) = 3e33
                    end do
                end do
            end if
        end if
        end subroutine
