        subroutine takerelanom(data,mean,npermax,yrbeg,yrend,nens1,nens2
     +       ,nperyear,lsum,lwrite)
!       take relative anomalies.
!       the variable lsum indocates the length of the seaon over
!       which the climatology should be avraged: annual, bi
!       -annual,
!       seasonal, monthly.
        implicit none
        integer npermax,yrbeg,yrend,nens1,nens2,nperyear,lsum
        real data(npermax,yrbeg:yrend,0:nens2),mean(npermax)
        logical lwrite
        integer i,j,iens
        real s
!
        if ( lsum.eq.12 ) then
            s = sum(mean(1:nperyear))
            mean = s/nperyear
        else if ( lsum.eq.nperyear/2 ) then
            if ( nperyear.eq.12 ) then
                s = sum(mean(4:9))
                mean(4:9) = s/6
                s = sum(mean(1:3)) + sum(mean(10:12))
                mean(1:3) = s/6
                mean(10:12) = s/6
            else
                goto 903
            end if
        else if ( lsum.eq.3 ) then
            if ( nperyear.eq.12 ) then
                do i=2,4
                    s = sum(mean(3*(i-1):3*i-1))
                    mean(3*(i-1):3*i-1) = s/3
                end do
                s = mean(12) + mean(1) + mean(2)
                mean(12) = s/3
                mean(1) = s/3
                mean(2) = s/3
            else if ( nperyear.eq.4 ) then
                if ( lwrite ) print *,'do nothing'
            else
                goto 903
            end if
        else if ( lsum.eq.12 ) then
            if ( nperyear.eq.12 ) then
                if ( lwrite ) print *,'do nothing'
            else
                goto 903
            end if
        else
            goto 903
        end if
        do iens=nens1,nens2
            do i=yrbeg,yrend
                do j=1,nperyear
                    if ( data(j,i,iens).lt.1e30 .and. mean(j).ne.0 )
     +                   then
                        data(j,i,iens) = data(j,i,iens)/mean(j)
                    else
                        data(j,i,iens) = 3e33
                    end if
                end do
            end do
        end do
        goto 999
 903    write(0,*) 'plotdat: error: cannot take relative anoamlies with'
     +       ,'nperyear = ',nperyear,' and lsum = ',lsum,' yet'
        call abort
 999    continue
        end
