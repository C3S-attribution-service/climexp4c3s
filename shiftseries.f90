subroutine shiftseries(data,npermax,nperyear,yrbeg,yrend,nshift)

!   shift the time series over nshift units (normally months)

    implicit none
    integer,intent(in) :: npermax,nperyear,yrbeg,yrend,nshift
    real,intent(inout) :: data(npermax,yrbeg:yrend)
    integer :: i,j,ii,k

    if ( nshift > 0 ) then
        do i=yrend,yrbeg,-1
            do j=nperyear,1,-1
                k = j + nshift
                call normon(k,i,ii,nperyear)
                if ( ii <= yrend ) then
                    data(k,ii) = data(j,i)
                    data(j,i) = 3e33
                endif
            enddo
        enddo
    elseif ( nshift < 0 ) then
        do i=yrbeg,yrend
            do j=1,nperyear
                k = j + nshift
                call normon(k,i,ii,nperyear)
                if ( ii >= yrbeg ) then
                    data(k,ii) = data(j,i)
                    data(j,i) = 3e33
                endif
            enddo
        enddo
    endif
end subroutine shiftseries

subroutine shiftseriesyear(series,npermax,nperyear,yrbeg,yrend,dyr,lwrite)

!   shift the dates forward by dyr years

    implicit none
    integer,intent(in) :: npermax,nperyear,yrbeg,yrend,dyr
    real,intent(inout) :: series(npermax,yrbeg:yrend)
    logical,intent(in) :: lwrite
    integer :: k,j,dy,yr,d,dtemp
    logical :: allundef
    integer :: leap
    external leap

!!!lwrite = .true.
    dtemp = 0
    if ( dyr == 0 ) then
        return
    else if ( dyr > 0 ) then
        if ( nperyear /= 366 ) then
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
                    ! start with a straight copy if the final year is less than 2100
                    if ( allundef .and. series(dy,yr) < 1e33 ) then
                        allundef = .false. 
                        if ( yr+dyr < 2100 ) then
                            d = 0
                        else
                            ! but with an offset of 1 if the year is above 2100,
                            ! which is not a leap year
                            d = 1
                        end if
                    end if
                    if ( dy == 60 ) then
                        if ( leap(yr) == 2 .and. leap(yr+dyr) /= 2 ) then
                            d = d - 1
                            series(dy,yr+dyr) = 3e33
                            if ( lwrite .and. &
                            series(dy-3,yr) < 1e33 ) then
                                print *,'different leap years, changing d to ',d
                                print *,'series(',dy,yr+dyr,') = 3e33'
                            end if
                            if ( d == 0 ) dtemp = 1
                            if ( d < -1 .and. series(dy-3,yr) < 1e33 ) then
                                write(0,*) 'shiftseriesyear: untested case of two leap years',d
                            end if
                        end if
                        if ( leap(yr) /= 2 .and. leap(yr+dyr) == 2 ) then
                            d = d + 1
                            if ( lwrite .and. &
                            series(dy-3,yr) < 1e33 ) then
                                print *,'different leap years, changing d to ',d
                            end if
                            if ( d > 1 .and. series(dy-3,yr) < 1e33 ) then
                                write(0,*) 'shiftseriesyear: untested case of two leap years',d
                            end if
                            cycle
                        end if
                    end if
                    k = dy
                    if ( dy+d < 1 ) then
                        j = dy+d+nperyear
                        series(j,yr+dyr-1)=series(dy,yr)
                        if ( lwrite .and. series(dy,yr) < 1e33 ) &
                            print *,'copying from ',dy,yr,' to ',j,yr+dyr-1
                    else if ( dy+d > nperyear ) then
                        j = dy+d-nperyear
                        series(j,yr+dyr+1)=series(dy,yr)
                        if ( lwrite .and. series(dy,yr) < 1e33 ) &
                            print *,'copying from ',dy,yr,' to ',j,yr+dyr+1
                    else
                        j = dy+d+dtemp
                        if ( j == 60 .and. k /= 60 &
                         .and. leap(yr+dyr) == 1 ) then
                            j = j+d
                        end if
                        if ( j /= 60 .and. k == 60 &
                         .and. leap(yr) == 1 ) then
                            k = k-d
                        end if
                        if ( series(k,yr) < 1e33 ) then
                            series(j,yr+dyr) = series(k,yr)
                            if ( lwrite .and. (j < 3 .or. j > 363 .or. abs(j-60) < 3 ) ) then
                                print *,'copying from ',k,yr,' to ',j,yr+dyr
                            end if
                        end if
                    end if
                    dtemp = 0
                    series(k,yr) = 3e33
                end do
            end do
        end if
    else if ( dyr < 0 ) then
        if ( nperyear /= 366 ) then
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
                    if ( allundef .and. series(dy,yr) < 1e33 ) then
                        allundef = .false. 
                        d = 0
                    end if
                    if ( dy == 60 ) then
                        if ( leap(yr) == 2 .and. leap(yr+dyr) /= 2 ) then
                            d = d - 1
                            series(dy,yr+dyr) = 3e33
                            if ( lwrite .and. series(dy+3,yr) < 1e33 ) then
                                print *,'different leap years, changing d to ',d
                                print *,'series(',dy,yr+dyr,') = 3e33'
                            end if
                            if ( d == 0 ) dtemp = 1
                        end if
                        if ( leap(yr) /= 2 .and. leap(yr+dyr) == 2 ) then
                            d = d + 1
                            if ( lwrite .and. series(dy+3,yr) < 1e33 ) then
                                print *,'different leap years, changing d to ',d
                            end if
                            cycle
                        end if
                    end if
                    if ( dy-d < 1 ) then
                        j = dy-d+nperyear
                        series(j,yr+dyr-1)=series(dy,yr)
                        if ( lwrite .and. series(dy,yr) < 1e33 ) &
                            print *,'copying from ',dy,yr,' to ',j,yr+dyr-1
                    else if ( dy-d > nperyear ) then
                        j = dy-d-nperyear
                        series(j,yr+dyr+1)=series(dy,yr)
                        if ( lwrite .and. series(dy,yr) < 1e33 ) &
                            print *,'copying from ',dy,yr,' to ',j,yr+dyr+1
                    else
                        j = dy-d-dtemp
                        if ( j == 60 .and. dy /= 60 .and. leap(yr+dyr) == 1 ) then
                            series(j,yr+dyr) = series(dy-1,yr)
                        else
                            series(j,yr+dyr) = series(dy,yr)
                        end if
                        if ( j < 3 .or. j > 363 .or. abs(j-60) < 3 ) then
                            if ( lwrite .and. series(dy,yr) < 1e33) &
                                print *,'copying from ',dy,yr,' to ',j,yr+dyr
                        end if
                    end if
                    dtemp = 0
                    series(dy,yr) = 3e33
                end do
            end do
        end if
    end if
end subroutine shiftseriesyear
